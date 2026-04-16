#%%
"""
UC1 - Semantic proximity predicts cross-layer quantitative coherence

Hypothesis: Gene pairs connected in the KG by specific semantic relationships
(e.g., "causally influences", "molecularly interacts with") show higher
quantitative coherence across omic layers than random gene pairs.

Architecture:
1. Query ArangoDB for gene-gene pairs connected by each predicate type
2. For each pair, compute Pearson correlation of mRNA, protein, and
   cross-layer mRNA-protein vectors across BRCA samples
3. Compare semantically connected pairs vs size-matched random pairs
4. Stratify by predicate type

Visualizations:
- Boxplot: correlation distributions by predicate type vs random
- Network: KG subgraph with edges colored by mRNA-protein correlation strength
- Scatter: semantic distance (hop count) vs cross-layer correlation
"""

#%% Imports
import os
import sys
import json
import random
import warnings
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx

# Add scripts dir to path for local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

from arangodb_utils import setup_arangodb_connection
from query_utils import (
    OMIC_CONFIG,
    get_gene_expression_matrix,
    get_protein_abundance,
    resolve_protein_position,
    list_samples_with_complete_omics,
    get_gene_vector,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
DB_NAME = "PKT_test10000"
COHORT = "TCGA-BRCA"

# Predicate groups to test — each requires a different extraction strategy
# because the KG encodes relationships through different class_code topologies:
#
# DIRECT gene-gene:
#   "genetically interacts with"  EntrezID <-> EntrezID
#
# VIA PROTEIN (PR->PR, resolve each PR to its gene via "has gene product"):
#   "molecularly interacts with"  PR <-> PR
#
# VIA PROTEIN-COMPOUND (PR<->CHEBI, resolve PR to gene):
#   "interacts with"              CHEBI <-> PR  (drug-target)
#
# CO-PARTICIPATION (genes sharing the same intermediate node):
#   "has participant"             R-HSA -> EntrezID (pathway co-membership)
#   "has phenotype"               MONDO -> HP (disease-phenotype, needs gene->disease)
#   "causes or contributes to"    EntrezID -> MONDO (gene pairs sharing disease)
#
# CAUSAL VIA SNP (dbSNP->EntrezID, but SNP maps back to same gene = self-loop):
#   "causally influences"         dbSNP -> EntrezID  (excluded: no gene PAIRS)
#
# Extraction strategy key:
#   "direct"       = both endpoints are EntrezID
#   "via_protein"  = PR<->PR, resolve to gene
#   "co_pathway"   = genes co-participating in same R-HSA pathway
#   "co_disease"   = genes contributing to same MONDO disease
PREDICATE_GROUPS = {
    "genetically interacts with": {
        "predicates": ["genetically interacts with"],
        "strategy": "direct",
    },
    "molecularly interacts with": {
        "predicates": ["molecularly interacts with"],
        "strategy": "via_protein",
    },
    # NOTE: "interacts with" excluded — edges are CHEBI<->PR (drug-target),
    # not gene-gene. CHEBI compounds don't resolve to EntrezID genes.
    "co-pathway (has participant)": {
        "predicates": ["has participant"],
        "strategy": "co_pathway",
    },
    "co-disease": {
        "predicates": ["causes or contributes to condition"],
        "strategy": "co_disease",
    },
}

# How many gene pairs to sample per predicate group
MAX_PAIRS_PER_PREDICATE = 500
# How many random pairs for baseline
N_RANDOM_PAIRS = 500
# Minimum samples required for correlation
MIN_SAMPLES_CORR = 20
# Random seed for reproducibility
SEED = 42

# Output directory
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "results", "uc1")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# STEP 1: Extract gene pairs from KG by predicate type
# =============================================================================

def _dedup_pairs(results: List[Dict]) -> List[Tuple[str, str, str, str]]:
    """Deduplicate gene pairs, removing self-loops."""
    pairs = []
    seen = set()
    for r in results:
        a, b = r["gene_a"], r["gene_b"]
        if a == b or not a or not b:
            continue
        key = tuple(sorted([a, b]))
        if key not in seen:
            seen.add(key)
            pairs.append((a, b, r.get("gene_a_label", ""), r.get("gene_b_label", "")))
    return pairs


def get_gene_pairs_direct(db, predicate_labels: List[str],
                           limit: int = 5000) -> List[Tuple[str, str, str, str]]:
    """
    Strategy 'direct': both endpoints are EntrezID gene nodes.
    Works for: genetically interacts with.
    """
    aql = """
    FOR e IN edges
        FILTER e.predicate_label IN @predicates
        LET src = DOCUMENT(e._from)
        LET tgt = DOCUMENT(e._to)
        FILTER src.class_code == "EntrezID" AND tgt.class_code == "EntrezID"
        LIMIT @limit
        RETURN {
            gene_a: src.entity_id, gene_b: tgt.entity_id,
            gene_a_label: src.label, gene_b_label: tgt.label
        }
    """
    results = list(db.aql.execute(aql, bind_vars={
        "predicates": predicate_labels, "limit": limit,
    }))
    return _dedup_pairs(results)


def _build_pr_to_gene_map(db) -> Dict[str, Dict]:
    """
    Build a lookup: PR node _key -> {entity_id, label} of its EntrezID gene.
    Single bulk AQL query — avoids per-edge traversal that causes timeouts.
    """
    aql = """
    FOR e IN edges
        FILTER e.predicate_label IN ["gene product of", "has gene product"]
        LET src = DOCUMENT(e._from)
        LET tgt = DOCUMENT(e._to)
        LET pr_node  = (src.class_code == "PR") ? src : tgt
        LET gene_node = (src.class_code == "EntrezID") ? src : tgt
        FILTER pr_node.class_code == "PR" AND gene_node.class_code == "EntrezID"
        RETURN DISTINCT {pr_key: pr_node._key, gene_id: gene_node.entity_id, gene_label: gene_node.label}
    """
    results = list(db.aql.execute(aql))
    return {r["pr_key"]: {"entity_id": r["gene_id"], "label": r["gene_label"]} for r in results}


# Module-level cache so we only build the map once per run
_PR_TO_GENE_CACHE: Optional[Dict] = None


def get_gene_pairs_via_protein(db, predicate_labels: List[str],
                                limit: int = 5000) -> List[Tuple[str, str, str, str]]:
    """
    Strategy 'via_protein': edges connect PR<->PR (or CHEBI<->PR).
    Two-step batch approach (avoids per-edge traversal timeout):
    1. Build global PR->gene lookup (one query, cached)
    2. Fetch PR-involving edges, resolve both endpoints in Python
    Works for: molecularly interacts with, interacts with.
    """
    global _PR_TO_GENE_CACHE
    if _PR_TO_GENE_CACHE is None:
        print("    Building PR -> gene lookup map...")
        _PR_TO_GENE_CACHE = _build_pr_to_gene_map(db)
        print(f"    PR->gene map: {len(_PR_TO_GENE_CACHE)} proteins mapped")
    pr_to_gene = _PR_TO_GENE_CACHE

    # Fetch edges where at least one endpoint is PR
    aql = """
    FOR e IN edges
        FILTER e.predicate_label IN @predicates
        LET src = DOCUMENT(e._from)
        LET tgt = DOCUMENT(e._to)
        FILTER src.class_code == "PR" OR tgt.class_code == "PR"
        LIMIT @limit
        RETURN {
            src_key: src._key, src_cc: src.class_code,
            src_eid: src.entity_id, src_label: src.label,
            tgt_key: tgt._key, tgt_cc: tgt.class_code,
            tgt_eid: tgt.entity_id, tgt_label: tgt.label
        }
    """
    edges = list(db.aql.execute(aql, bind_vars={
        "predicates": predicate_labels, "limit": limit * 3,
    }))
    print(f"    Fetched {len(edges)} edges with PR endpoints")

    # Resolve each endpoint to gene via lookup
    results = []
    for e in edges:
        if e["src_cc"] == "EntrezID":
            gene_a = {"entity_id": e["src_eid"], "label": e["src_label"]}
        elif e["src_cc"] == "PR":
            gene_a = pr_to_gene.get(e["src_key"])
        else:
            gene_a = None

        if e["tgt_cc"] == "EntrezID":
            gene_b = {"entity_id": e["tgt_eid"], "label": e["tgt_label"]}
        elif e["tgt_cc"] == "PR":
            gene_b = pr_to_gene.get(e["tgt_key"])
        else:
            gene_b = None

        if gene_a and gene_b and gene_a["entity_id"] != gene_b["entity_id"]:
            results.append({
                "gene_a": gene_a["entity_id"], "gene_b": gene_b["entity_id"],
                "gene_a_label": gene_a.get("label", ""), "gene_b_label": gene_b.get("label", ""),
            })

    return _dedup_pairs(results)


def get_gene_pairs_co_pathway(db, predicate_labels: List[str],
                               limit: int = 5000,
                               max_pathway_size: int = 50) -> List[Tuple[str, str, str, str]]:
    """
    Strategy 'co_pathway': genes that co-participate in the same R-HSA pathway.
    Edges are R-HSA -[has participant]-> EntrezID.
    We group genes by pathway, then form pairs within each pathway.
    max_pathway_size caps very large pathways to avoid combinatorial explosion.
    """
    aql = """
    FOR e IN edges
        FILTER e.predicate_label IN @predicates
        LET pathway = DOCUMENT(e._from)
        LET gene = DOCUMENT(e._to)
        FILTER pathway.class_code == "R-HSA" AND gene.class_code == "EntrezID"
        COLLECT pw = pathway._key INTO gene_docs = gene
        LET unique_genes = (
            FOR g IN gene_docs
                COLLECT gid = g.entity_id INTO labels = g.label
                RETURN {entity_id: gid, label: FIRST(labels)}
        )
        FILTER LENGTH(unique_genes) >= 2 AND LENGTH(unique_genes) <= @maxSize
        RETURN {pathway: pw, genes: unique_genes}
    """
    pathways = list(db.aql.execute(aql, bind_vars={
        "predicates": predicate_labels,
        "maxSize": max_pathway_size,
    }))

    all_pairs = []
    for pw in pathways:
        genes = pw["genes"]
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                all_pairs.append({
                    "gene_a": genes[i]["entity_id"],
                    "gene_b": genes[j]["entity_id"],
                    "gene_a_label": genes[i].get("label", ""),
                    "gene_b_label": genes[j].get("label", ""),
                })
                if len(all_pairs) >= limit:
                    return _dedup_pairs(all_pairs)
    return _dedup_pairs(all_pairs)


def get_gene_pairs_co_disease(db, predicate_labels: List[str],
                               limit: int = 5000,
                               max_disease_genes: int = 30) -> List[Tuple[str, str, str, str]]:
    """
    Strategy 'co_disease': genes that both cause/contribute to the same MONDO disease.
    Edges are EntrezID -[causes or contributes to condition]-> MONDO.
    Group genes by disease, form pairs.
    """
    aql = """
    FOR e IN edges
        FILTER e.predicate_label IN @predicates
        LET gene = DOCUMENT(e._from)
        LET disease = DOCUMENT(e._to)
        FILTER gene.class_code == "EntrezID" AND disease.class_code == "MONDO"
        COLLECT dis = disease._key INTO gene_docs = gene
        LET unique_genes = (
            FOR g IN gene_docs
                COLLECT gid = g.entity_id INTO labels = g.label
                RETURN {entity_id: gid, label: FIRST(labels)}
        )
        FILTER LENGTH(unique_genes) >= 2 AND LENGTH(unique_genes) <= @maxSize
        RETURN {disease: dis, genes: unique_genes}
    """
    diseases = list(db.aql.execute(aql, bind_vars={
        "predicates": predicate_labels,
        "maxSize": max_disease_genes,
    }))

    all_pairs = []
    for d in diseases:
        genes = d["genes"]
        for i in range(len(genes)):
            for j in range(i + 1, len(genes)):
                all_pairs.append({
                    "gene_a": genes[i]["entity_id"],
                    "gene_b": genes[j]["entity_id"],
                    "gene_a_label": genes[i].get("label", ""),
                    "gene_b_label": genes[j].get("label", ""),
                })
                if len(all_pairs) >= limit:
                    return _dedup_pairs(all_pairs)
    return _dedup_pairs(all_pairs)


def get_gene_pairs_by_predicate(db, predicate_config: Dict,
                                 limit: int = 5000) -> List[Tuple[str, str, str, str]]:
    """
    Dispatch to the correct extraction strategy based on predicate topology.
    """
    strategy = predicate_config["strategy"]
    predicates = predicate_config["predicates"]

    if strategy == "direct":
        return get_gene_pairs_direct(db, predicates, limit)
    elif strategy == "via_protein":
        return get_gene_pairs_via_protein(db, predicates, limit)
    elif strategy == "co_pathway":
        return get_gene_pairs_co_pathway(db, predicates, limit)
    elif strategy == "co_disease":
        return get_gene_pairs_co_disease(db, predicates, limit)
    else:
        print(f"  WARNING: unknown strategy '{strategy}', skipping")
        return []


def get_gene_pairs_by_hop_distance(db, gene_key: str, max_hops: int = 3,
                                    limit_per_hop: int = 200) -> Dict[int, List[str]]:
    """
    For a given gene node, find other EntrezID genes reachable at 1, 2, 3 hops.
    Returns dict {hop_distance: [entrez_ids]}.
    """
    results_by_hop = {}
    for h in range(1, max_hops + 1):
        aql = f"""
        FOR v IN {h}..{h} ANY CONCAT("nodes/", @geneKey) edges
            FILTER v.class_code == "EntrezID"
            FILTER v._key != @geneKey
            LIMIT @limit
            RETURN DISTINCT v.entity_id
        """
        res = list(db.aql.execute(aql, bind_vars={
            "geneKey": gene_key, "limit": limit_per_hop
        }))
        results_by_hop[h] = res
    return results_by_hop


# =============================================================================
# STEP 2: Build expression / protein matrices for gene pairs
# =============================================================================

def build_paired_vectors(db, gene_pairs: List[Tuple[str, str]],
                         cohort: str,
                         sample_ids: List[str] = None
                         ) -> pd.DataFrame:
    """
    For each gene pair, retrieve mRNA expression vectors and protein abundance
    vectors across samples. Returns a DataFrame with columns:
    gene_a, gene_b, corr_mrna, corr_protein, corr_cross (mRNA_a vs protein_b),
    n_samples_mrna, n_samples_protein.
    """
    # Collect all unique entrez IDs needed
    all_genes = sorted(set(g for pair in gene_pairs for g in pair[:2]))
    print(f"  Fetching expression matrix for {len(all_genes)} genes across {len(sample_ids) if sample_ids else 'all'} samples...")

    # Batch fetch expression matrix
    expr_data = get_gene_expression_matrix(
        db, cohort, gene_ids=all_genes,
        id_type="entrez_id", value_type="tpm",
        sample_ids=sample_ids
    )

    expr_df = None
    if expr_data["samples"] and expr_data["genes"]:
        expr_df = pd.DataFrame(
            expr_data["matrix"],
            index=expr_data["samples"],
            columns=expr_data["genes"]
        ).apply(pd.to_numeric, errors="coerce")
        print(f"  Expression matrix: {expr_df.shape[0]} samples x {expr_df.shape[1]} genes")

    # Fetch protein data gene by gene (protein index uses gene_symbol, not entrez_id
    # but some mappings use entrez_id - try both)
    print(f"  Fetching protein abundance for {len(all_genes)} genes...")
    prot_records = {}
    for eid in all_genes:
        data = get_protein_abundance(db, cohort, entrez_id=eid, sample_ids=sample_ids)
        if data:
            prot_records[eid] = {d["sample_id"]: d["abundance"] for d in data}

    prot_df = None
    if prot_records:
        prot_df = pd.DataFrame(prot_records).apply(pd.to_numeric, errors="coerce")
        print(f"  Protein matrix: {prot_df.shape[0]} samples x {prot_df.shape[1]} genes")

    # Compute correlations for each pair
    rows = []
    for gene_a, gene_b, *labels in gene_pairs:
        row = {
            "gene_a": gene_a,
            "gene_b": gene_b,
            "gene_a_label": labels[0] if labels else "",
            "gene_b_label": labels[1] if len(labels) > 1 else "",
            "corr_mrna": np.nan,
            "corr_protein": np.nan,
            "corr_cross_ab": np.nan,  # mRNA_a vs protein_b
            "corr_cross_ba": np.nan,  # mRNA_b vs protein_a
            "n_samples_mrna": 0,
            "n_samples_protein": 0,
        }

        # mRNA-mRNA correlation
        if expr_df is not None and gene_a in expr_df.columns and gene_b in expr_df.columns:
            va = expr_df[gene_a].dropna()
            vb = expr_df[gene_b].dropna()
            common = va.index.intersection(vb.index)
            if len(common) >= MIN_SAMPLES_CORR:
                r, _ = stats.pearsonr(va.loc[common], vb.loc[common])
                row["corr_mrna"] = r
                row["n_samples_mrna"] = len(common)

        # Protein-protein correlation
        if prot_df is not None and gene_a in prot_df.columns and gene_b in prot_df.columns:
            va = prot_df[gene_a].dropna()
            vb = prot_df[gene_b].dropna()
            common = va.index.intersection(vb.index)
            if len(common) >= MIN_SAMPLES_CORR:
                r, _ = stats.pearsonr(va.loc[common], vb.loc[common])
                row["corr_protein"] = r
                row["n_samples_protein"] = len(common)

        # Cross-layer: mRNA_a vs protein_b
        if expr_df is not None and prot_df is not None:
            if gene_a in expr_df.columns and gene_b in prot_df.columns:
                va = expr_df[gene_a].dropna()
                vb = prot_df[gene_b].dropna()
                common = va.index.intersection(vb.index)
                if len(common) >= MIN_SAMPLES_CORR:
                    r, _ = stats.pearsonr(va.loc[common], vb.loc[common])
                    row["corr_cross_ab"] = r

            if gene_b in expr_df.columns and gene_a in prot_df.columns:
                va = expr_df[gene_b].dropna()
                vb = prot_df[gene_a].dropna()
                common = va.index.intersection(vb.index)
                if len(common) >= MIN_SAMPLES_CORR:
                    r, _ = stats.pearsonr(va.loc[common], vb.loc[common])
                    row["corr_cross_ba"] = r

        rows.append(row)

    return pd.DataFrame(rows)


# =============================================================================
# STEP 3: Generate random gene pairs as baseline
# =============================================================================

def generate_random_gene_pairs(all_gene_ids: List[str], n_pairs: int,
                                seed: int = 42) -> List[Tuple[str, str]]:
    """Generate random pairs of gene IDs (no self-pairs, no duplicates)."""
    rng = random.Random(seed)
    pairs = set()
    attempts = 0
    max_attempts = n_pairs * 10
    while len(pairs) < n_pairs and attempts < max_attempts:
        a, b = rng.sample(all_gene_ids, 2)
        key = tuple(sorted([a, b]))
        pairs.add(key)
        attempts += 1
    return [(a, b, "", "") for a, b in pairs]


# =============================================================================
# STEP 4: Hop-distance analysis
# =============================================================================

def compute_hop_distance_correlations(db, cohort: str,
                                       seed_genes: List[str],
                                       sample_ids: List[str],
                                       max_hops: int = 3,
                                       pairs_per_hop: int = 100,
                                       seed: int = 42) -> pd.DataFrame:
    """
    For a set of seed genes, find genes at 1, 2, 3 hops in the KG.
    Compute mRNA correlations for pairs at each hop distance.
    """
    rng = random.Random(seed)
    hop_pairs = defaultdict(list)

    for gene_key in seed_genes:
        neighbors_by_hop = get_gene_pairs_by_hop_distance(
            db, gene_key, max_hops=max_hops, limit_per_hop=200
        )
        for hop, neighbor_ids in neighbors_by_hop.items():
            for nid in neighbor_ids:
                if nid != gene_key:
                    hop_pairs[hop].append((gene_key, nid, "", ""))

    # Sample pairs per hop
    rows = []
    for hop in range(1, max_hops + 1):
        available = hop_pairs.get(hop, [])
        if not available:
            continue
        sampled = rng.sample(available, min(pairs_per_hop, len(available)))
        print(f"  Hop {hop}: {len(sampled)} pairs (from {len(available)} available)")

        df = build_paired_vectors(db, sampled, cohort, sample_ids)
        df["hop_distance"] = hop
        rows.append(df)

    if rows:
        return pd.concat(rows, ignore_index=True)
    return pd.DataFrame()


# =============================================================================
# VISUALIZATIONS
# =============================================================================

def plot_correlation_boxplots(results_by_group: Dict[str, pd.DataFrame],
                               corr_column: str = "corr_mrna",
                               title: str = "mRNA-mRNA Correlation by Semantic Relation",
                               output_path: str = None):
    """
    Boxplot comparing correlation distributions across predicate groups vs random.
    """
    plot_data = []
    group_order = []

    for group_name, df in results_by_group.items():
        vals = df[corr_column].dropna()
        if len(vals) == 0:
            continue
        for v in vals:
            plot_data.append({"Relation": group_name, "Correlation": v})
        group_order.append(group_name)

    if not plot_data:
        print(f"  No data for {corr_column}, skipping plot.")
        return

    plot_df = pd.DataFrame(plot_data)

    # Order: random first, then by median correlation descending
    medians = plot_df.groupby("Relation")["Correlation"].median()
    if "random" in medians.index:
        others = medians.drop("random").sort_values(ascending=False).index.tolist()
        order = ["random"] + others
    else:
        order = medians.sort_values(ascending=False).index.tolist()

    fig, ax = plt.subplots(figsize=(14, 7))

    # Color palette: grey for random, gradient for others
    palette = {}
    cmap = plt.cm.RdYlGn
    for i, name in enumerate(order):
        if name == "random":
            palette[name] = "#CCCCCC"
        else:
            palette[name] = cmap(0.3 + 0.7 * i / max(len(order) - 2, 1))

    sns.boxplot(data=plot_df, x="Relation", y="Correlation", order=order,
                palette=palette, ax=ax, showfliers=False, width=0.6)
    sns.stripplot(data=plot_df, x="Relation", y="Correlation", order=order,
                  color="black", alpha=0.15, size=2, ax=ax, jitter=True)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.set_xlabel("KG Semantic Relation", fontsize=12)
    ax.set_ylabel(f"Pearson r ({corr_column.replace('corr_', '')})", fontsize=12)
    ax.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)

    # Add sample sizes
    for i, name in enumerate(order):
        subset = plot_df[plot_df["Relation"] == name]["Correlation"].dropna()
        n = len(subset)
        med = subset.median()
        ax.text(i, ax.get_ylim()[1] * 0.95, f"n={n}\nmed={med:.3f}",
                ha="center", va="top", fontsize=8, color="dimgray")

    plt.xticks(rotation=35, ha="right")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def plot_hop_distance_scatter(hop_df: pd.DataFrame,
                               corr_column: str = "corr_mrna",
                               title: str = "KG Hop Distance vs Cross-Layer Correlation",
                               output_path: str = None):
    """
    Scatter + boxplot: semantic distance (hops) vs correlation.
    """
    df = hop_df.dropna(subset=[corr_column])
    if df.empty:
        print(f"  No data for hop distance plot ({corr_column}), skipping.")
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6),
                                     gridspec_kw={"width_ratios": [2, 1]})

    # Scatter with jitter
    jitter = np.random.default_rng(42).uniform(-0.2, 0.2, size=len(df))
    colors = {1: "#2ca02c", 2: "#ff7f0e", 3: "#d62728"}
    for hop in sorted(df["hop_distance"].unique()):
        subset = df[df["hop_distance"] == hop]
        j = jitter[:len(subset)]
        ax1.scatter(subset["hop_distance"] + j[:len(subset)],
                    subset[corr_column],
                    alpha=0.3, s=20, color=colors.get(hop, "grey"),
                    label=f"Hop {hop} (n={len(subset)})")

    # Trend line: median per hop
    medians = df.groupby("hop_distance")[corr_column].median()
    ax1.plot(medians.index, medians.values, "k-o", linewidth=2, markersize=8,
             label="Median", zorder=5)

    ax1.set_xlabel("Semantic Distance (KG hops)", fontsize=12)
    ax1.set_ylabel(f"Pearson r ({corr_column.replace('corr_', '')})", fontsize=12)
    ax1.set_title(title, fontsize=13, fontweight="bold")
    ax1.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax1.legend(fontsize=9)
    ax1.set_xticks(sorted(df["hop_distance"].unique()))

    # Boxplot
    sns.boxplot(data=df, x="hop_distance", y=corr_column, ax=ax2,
                palette=colors, showfliers=False)
    ax2.set_xlabel("Hop Distance", fontsize=12)
    ax2.set_ylabel("")
    ax2.set_title("Distribution by Hop", fontsize=12)
    ax2.axhline(y=0, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def plot_network_subgraph(results_df: pd.DataFrame,
                           predicate_name: str,
                           corr_column: str = "corr_mrna",
                           top_n: int = 80,
                           title: str = None,
                           output_path: str = None):
    """
    Network visualization: KG subgraph with edges colored by correlation strength.
    Shows top_n pairs with valid correlations for a given predicate.
    """
    df = results_df.dropna(subset=[corr_column])
    if df.empty:
        print(f"  No data for network plot, skipping.")
        return

    # Take top-N by absolute correlation
    df = df.reindex(df[corr_column].abs().sort_values(ascending=False).index).head(top_n)

    G = nx.Graph()
    for _, row in df.iterrows():
        label_a = row.get("gene_a_label", row["gene_a"]) or row["gene_a"]
        label_b = row.get("gene_b_label", row["gene_b"]) or row["gene_b"]
        # Use short labels (gene symbol)
        short_a = label_a.split(" (")[0] if " (" in str(label_a) else str(label_a)
        short_b = label_b.split(" (")[0] if " (" in str(label_b) else str(label_b)

        G.add_node(short_a, entrez=row["gene_a"])
        G.add_node(short_b, entrez=row["gene_b"])
        G.add_edge(short_a, short_b, weight=row[corr_column])

    if len(G.nodes) == 0:
        return

    fig, ax = plt.subplots(figsize=(14, 12))

    pos = nx.spring_layout(G, seed=42, k=2.0 / np.sqrt(len(G.nodes)), iterations=60)

    # Edge colors based on correlation
    edge_weights = [G[u][v]["weight"] for u, v in G.edges()]
    edge_colors = plt.cm.RdBu([(w + 1) / 2 for w in edge_weights])  # map [-1,1] to [0,1]
    edge_widths = [max(0.5, abs(w) * 3) for w in edge_weights]

    # Node sizes by degree
    degrees = dict(G.degree())
    node_sizes = [100 + degrees[n] * 40 for n in G.nodes()]

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color=edge_colors,
                           width=edge_widths, alpha=0.7)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_size=node_sizes,
                           node_color="#4A90D9", alpha=0.85, edgecolors="white",
                           linewidths=0.5)

    # Labels only for high-degree nodes to avoid clutter
    degree_threshold = np.percentile(list(degrees.values()), 70)
    labels = {n: n for n, d in degrees.items() if d >= degree_threshold}
    nx.draw_networkx_labels(G, pos, labels, font_size=7, font_weight="bold", ax=ax)

    # Colorbar
    sm = plt.cm.ScalarMappable(cmap=plt.cm.RdBu, norm=plt.Normalize(-1, 1))
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.6, pad=0.02)
    cbar.set_label(f"Pearson r ({corr_column.replace('corr_', '')})", fontsize=10)

    title = title or f"KG Subgraph: {predicate_name}\n(edges colored by {corr_column})"
    ax.set_title(title, fontsize=13, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def plot_summary_heatmap(summary_stats: pd.DataFrame, output_path: str = None):
    """
    Heatmap: predicate type vs correlation metric (median values).
    """
    metrics = ["corr_mrna", "corr_protein", "corr_cross_ab"]
    metric_labels = ["mRNA-mRNA", "Protein-Protein", "mRNA-Protein"]

    data = []
    for _, row in summary_stats.iterrows():
        for m, ml in zip(metrics, metric_labels):
            col = f"median_{m}"
            if col in row:
                data.append({
                    "Relation": row["predicate_group"],
                    "Metric": ml,
                    "Median r": row[col],
                })
    if not data:
        return

    hm_df = pd.DataFrame(data).pivot(index="Relation", columns="Metric", values="Median r")

    # Order by mRNA-mRNA median
    if "mRNA-mRNA" in hm_df.columns:
        hm_df = hm_df.sort_values("mRNA-mRNA", ascending=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.heatmap(hm_df, annot=True, fmt=".3f", cmap="RdBu", center=0,
                linewidths=0.5, ax=ax, vmin=-0.3, vmax=0.3)
    ax.set_title("Median Correlation by Semantic Relation and Omic Layer",
                 fontsize=13, fontweight="bold")
    ax.set_ylabel("")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def print_statistical_tests(results_by_group: Dict[str, pd.DataFrame],
                             corr_column: str = "corr_mrna"):
    """
    Mann-Whitney U test: each predicate group vs random.
    """
    random_vals = results_by_group.get("random", pd.DataFrame())[corr_column].dropna()
    if len(random_vals) < 5:
        print("  Not enough random pairs for statistical testing.")
        return

    print(f"\n{'='*70}")
    print(f"  Mann-Whitney U test: each predicate vs random ({corr_column})")
    print(f"{'='*70}")
    print(f"  {'Predicate':<35} {'n':>5}  {'Median':>8}  {'U stat':>10}  {'p-value':>10}  {'Effect size r':>12}")
    print(f"  {'-'*35} {'-'*5}  {'-'*8}  {'-'*10}  {'-'*10}  {'-'*12}")

    random_med = random_vals.median()
    print(f"  {'random (baseline)':<35} {len(random_vals):>5}  {random_med:>8.4f}")

    for group_name, df in results_by_group.items():
        if group_name == "random":
            continue
        vals = df[corr_column].dropna()
        if len(vals) < 5:
            continue
        u_stat, p_val = stats.mannwhitneyu(vals, random_vals, alternative="greater")
        # Rank-biserial correlation as effect size (derived directly from U)
        n1, n2 = len(vals), len(random_vals)
        effect_r = 2 * u_stat / (n1 * n2) - 1  # ranges from -1 to +1

        med = vals.median()
        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else ""
        print(f"  {group_name:<35} {len(vals):>5}  {med:>8.4f}  {u_stat:>10.0f}  {p_val:>10.2e}  {effect_r:>10.4f}  {sig}")


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_uc1(db_name: str = DB_NAME, cohort: str = COHORT):
    """
    Execute the full UC1 analysis pipeline.
    """
    print("=" * 70)
    print("  UC1: Semantic Proximity Predicts Cross-Layer Quantitative Coherence")
    print("=" * 70)

    # --- Connect to ArangoDB ---
    db = setup_arangodb_connection(db_name)
    if db is None:
        raise RuntimeError("Failed to connect to ArangoDB")

    # --- Get list of samples with complete omics ---
    print("\n[1/6] Finding samples with complete omics...")
    # Note: specimen_type in SAMPLES is "Solid Tissue", not "Primary Tumor".
    # The tumor/normal distinction is in the sample_type field.
    # First get all samples with both expression + protein, then filter to tumors.
    all_sample_ids = list_samples_with_complete_omics(
        db, cohort=cohort,
        omic_types=["gene_expression", "protein"],
        specimen_type=None,  # no filter at this level
    )
    print(f"  Found {len(all_sample_ids)} samples with both expression + protein data")

    # Filter to Primary Tumor via SAMPLES.sample_type
    tumor_aql = """
    FOR s IN SAMPLES
        FILTER s.sample_type == "Primary Tumor"
        RETURN s._key
    """
    tumor_keys = set(db.aql.execute(tumor_aql))
    sample_ids = [s for s in all_sample_ids if s in tumor_keys]
    print(f"  Of which {len(sample_ids)} are Primary Tumor samples")

    if len(sample_ids) < MIN_SAMPLES_CORR:
        print(f"  ERROR: Need at least {MIN_SAMPLES_CORR} samples, got {len(sample_ids)}")
        return

    # --- Extract gene pairs by predicate ---
    print("\n[2/6] Extracting gene pairs from KG by predicate type...")
    results_by_group = {}
    all_connected_genes = set()

    for group_name, pred_config in PREDICATE_GROUPS.items():
        print(f"\n  --- {group_name} (strategy: {pred_config['strategy']}) ---")
        pairs = get_gene_pairs_by_predicate(db, pred_config, limit=MAX_PAIRS_PER_PREDICATE * 3)
        print(f"  Found {len(pairs)} unique gene-gene pairs")

        if not pairs:
            continue

        # Sample down if too many
        if len(pairs) > MAX_PAIRS_PER_PREDICATE:
            rng = random.Random(SEED)
            pairs = rng.sample(pairs, MAX_PAIRS_PER_PREDICATE)
            print(f"  Sampled down to {len(pairs)} pairs")

        for a, b, *_ in pairs:
            all_connected_genes.add(a)
            all_connected_genes.add(b)

        # Compute correlations
        df = build_paired_vectors(db, pairs, cohort, sample_ids)
        df["predicate_group"] = group_name
        results_by_group[group_name] = df
        valid = df["corr_mrna"].notna().sum()
        print(f"  Valid mRNA correlations: {valid}/{len(df)}")

    # --- Generate random baseline ---
    print("\n[3/6] Generating random gene pair baseline...")
    # Get all gene entrez IDs available in expression data
    gene_vector = get_gene_vector(db, cohort, id_type="entrez_id", omic_type="gene_expression")
    valid_genes = [g for g in gene_vector if g and g != ""]
    print(f"  Total genes in expression index: {len(valid_genes)}")

    random_pairs = generate_random_gene_pairs(valid_genes, N_RANDOM_PAIRS, seed=SEED)
    random_df = build_paired_vectors(db, random_pairs, cohort, sample_ids)
    random_df["predicate_group"] = "random"
    results_by_group["random"] = random_df
    valid = random_df["corr_mrna"].notna().sum()
    print(f"  Valid random mRNA correlations: {valid}/{len(random_df)}")

    # --- Combine all results ---
    all_results = pd.concat(results_by_group.values(), ignore_index=True)
    all_results.to_csv(os.path.join(OUTPUT_DIR, "uc1_correlation_results.csv"), index=False)
    print(f"\n  Saved all results to {OUTPUT_DIR}/uc1_correlation_results.csv")

    # --- Summary statistics ---
    print("\n[4/6] Computing summary statistics...")
    summary_rows = []
    for group_name, df in results_by_group.items():
        row = {
            "predicate_group": group_name,
            "n_pairs": len(df),
            "n_valid_mrna": df["corr_mrna"].notna().sum(),
            "n_valid_protein": df["corr_protein"].notna().sum(),
            "n_valid_cross": df["corr_cross_ab"].notna().sum(),
            "median_corr_mrna": df["corr_mrna"].median(),
            "mean_corr_mrna": df["corr_mrna"].mean(),
            "std_corr_mrna": df["corr_mrna"].std(),
            "median_corr_protein": df["corr_protein"].median(),
            "mean_corr_protein": df["corr_protein"].mean(),
            "median_corr_cross_ab": df["corr_cross_ab"].median(),
            "mean_corr_cross_ab": df["corr_cross_ab"].mean(),
        }
        summary_rows.append(row)

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "uc1_summary_stats.csv"), index=False)
    print(summary_df.to_string(index=False))

    # --- Statistical tests ---
    print_statistical_tests(results_by_group, "corr_mrna")
    print_statistical_tests(results_by_group, "corr_protein")
    print_statistical_tests(results_by_group, "corr_cross_ab")

    # --- Visualizations ---
    print("\n[5/6] Generating visualizations...")

    # 1. Boxplot: mRNA-mRNA correlation by predicate
    plot_correlation_boxplots(
        results_by_group, corr_column="corr_mrna",
        title="UC1: mRNA-mRNA Correlation by KG Semantic Relation (TCGA-BRCA)",
        output_path=os.path.join(OUTPUT_DIR, "uc1_boxplot_mrna.png")
    )

    # 2. Boxplot: Protein-protein correlation by predicate
    plot_correlation_boxplots(
        results_by_group, corr_column="corr_protein",
        title="UC1: Protein-Protein Correlation by KG Semantic Relation (TCGA-BRCA)",
        output_path=os.path.join(OUTPUT_DIR, "uc1_boxplot_protein.png")
    )

    # 3. Boxplot: Cross-layer mRNA-protein
    plot_correlation_boxplots(
        results_by_group, corr_column="corr_cross_ab",
        title="UC1: Cross-Layer mRNA-Protein Correlation by KG Semantic Relation (TCGA-BRCA)",
        output_path=os.path.join(OUTPUT_DIR, "uc1_boxplot_cross_layer.png")
    )

    # 4. Summary heatmap
    plot_summary_heatmap(
        summary_df,
        output_path=os.path.join(OUTPUT_DIR, "uc1_heatmap_summary.png")
    )

    # 5. Network subgraph for top predicate group
    best_group = summary_df.loc[
        summary_df["predicate_group"] != "random", "median_corr_mrna"
    ].idxmax()
    best_name = summary_df.loc[best_group, "predicate_group"]
    if best_name in results_by_group:
        plot_network_subgraph(
            results_by_group[best_name],
            predicate_name=best_name,
            corr_column="corr_mrna",
            top_n=80,
            output_path=os.path.join(OUTPUT_DIR, f"uc1_network_{best_name.replace(' ', '_')}.png")
        )

    # --- Hop-distance analysis ---
    print("\n[6/6] Computing hop-distance analysis...")
    # Pick seed genes: high-degree genes in our results
    seed_candidates = list(all_connected_genes)[:30]
    if seed_candidates:
        hop_df = compute_hop_distance_correlations(
            db, cohort, seed_genes=seed_candidates[:15],
            sample_ids=sample_ids, max_hops=3,
            pairs_per_hop=150, seed=SEED
        )
        if not hop_df.empty:
            hop_df.to_csv(os.path.join(OUTPUT_DIR, "uc1_hop_distance_results.csv"), index=False)
            plot_hop_distance_scatter(
                hop_df, corr_column="corr_mrna",
                title="UC1: KG Semantic Distance vs mRNA Correlation (TCGA-BRCA)",
                output_path=os.path.join(OUTPUT_DIR, "uc1_scatter_hop_distance.png")
            )

    print("\n" + "=" * 70)
    print("  UC1 analysis complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print("=" * 70)

    return results_by_group, summary_df


#%% Entry point
if __name__ == "__main__":
    if "ipykernel" not in sys.modules:
        run_uc1()
    else:
        # Interactive: run with defaults
        results_by_group, summary_df = run_uc1()
