#%%
"""
UC2 - CNV-expression-protein discordance with semantic explanation

Hypothesis: Genes showing discordance between DNA copy number amplification,
mRNA expression, and protein abundance are enriched for specific semantic
properties in the Knowledge Graph (GO terms, pathways, interaction types)
that explain dosage compensation mechanisms.

Architecture:
1. For each gene in the RPPA protein index (~487 proteins):
   - Compute Pearson correlation CNV -> mRNA across BRCA samples
   - Compute Pearson correlation mRNA -> protein across BRCA samples
2. Classify genes:
   - Concordant: both correlations high (CNV drives mRNA drives protein)
   - Transcriptional discordant: CNV-mRNA low (compensation at transcription)
   - Post-transcriptional discordant: CNV-mRNA high but mRNA-protein low
3. For discordant genes: query KG for GO terms, pathways, predicates
4. Enrichment analysis: are discordant genes enriched in specific GO processes?

Visualizations:
- Scatter plot 2D: CNV-mRNA correlation (x) vs mRNA-protein correlation (y)
- Network: discordant genes in KG, colored by GO process type
- Heatmap: genes x omic layers for discordant vs concordant subsets
"""

#%% Imports
import os
import sys
import json
import warnings
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import fisher_exact
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import networkx as nx

# Add scripts dir to path for local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

from arangodb_utils import setup_arangodb_connection
from query_utils import (
    OMIC_CONFIG,
    get_gene_expression_by_gene,
    get_cnv_by_gene,
    get_protein_abundance,
    get_gene_vector,
    list_samples_with_complete_omics,
    resolve_cnv_position,
    resolve_gene_position,
    resolve_protein_position,
    get_node_neighbors,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
DB_NAME = "PKT_test10000"
COHORT = "TCGA-BRCA"

# Correlation thresholds for classification
# Genes with |r| above this are considered "correlated" for that layer transition
CORR_THRESHOLD_HIGH = 0.3
# Genes with |r| below this are considered "uncorrelated" / compensated
CORR_THRESHOLD_LOW = 0.15

# Minimum samples required for correlation
MIN_SAMPLES_CORR = 20

# Random seed
SEED = 42

# Output directory
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "results", "uc2")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# STEP 1: Get protein index genes and build multi-omic matrices
# =============================================================================

def get_protein_index_genes(db, cohort: str) -> List[Dict]:
    """
    Retrieve the full protein (RPPA) index for a cohort.
    Returns list of dicts with entrez_id, gene_symbol, position, etc.
    """
    config = OMIC_CONFIG["protein"]
    index_key = f"{config['index_key_prefix']}_{cohort}"

    aql = f"""
    FOR doc IN {config['collection']}
        FILTER doc.data_type == @indexDataType
        FILTER doc._key == @indexKey
        RETURN doc.{config['mapping_field']}
    """
    result = list(db.aql.execute(aql, bind_vars={
        "indexKey": index_key,
        "indexDataType": config["index_data_type"],
    }))

    if not result or not result[0]:
        print("  WARNING: No protein index found")
        return []

    mappings = result[0]
    print(f"  Protein index: {len(mappings)} entries")
    return mappings


def build_cross_layer_matrices(db, cohort: str, protein_genes: List[Dict],
                                sample_ids: List[str]
                                ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Build aligned CNV, mRNA, and protein matrices for the protein-index genes.
    Returns three DataFrames (samples x genes) with entrez_id as column names.
    """
    # Collect unique entrez_ids from protein index
    gene_info = {}
    for m in protein_genes:
        eid = m.get("entrez_id")
        sym = m.get("gene_symbol", "")
        if eid:
            gene_info[eid] = sym

    entrez_ids = sorted(gene_info.keys())
    print(f"  Unique entrez IDs in protein index: {len(entrez_ids)}")

    # --- CNV matrix (gene by gene, no batch function available) ---
    print(f"  Fetching CNV data for {len(entrez_ids)} genes...")
    cnv_records = {}
    for eid in entrez_ids:
        data = get_cnv_by_gene(db, cohort, entrez_id=eid, sample_ids=sample_ids)
        if data:
            cnv_records[eid] = {d["sample_id"]: d["copy_number"] for d in data}

    cnv_df = pd.DataFrame(cnv_records).apply(pd.to_numeric, errors="coerce") if cnv_records else pd.DataFrame()
    print(f"  CNV matrix: {cnv_df.shape[0]} samples x {cnv_df.shape[1]} genes")

    # --- mRNA expression matrix (gene by gene for consistency) ---
    print(f"  Fetching mRNA expression data for {len(entrez_ids)} genes...")
    expr_records = {}
    for eid in entrez_ids:
        data = get_gene_expression_by_gene(db, cohort, entrez_id=eid,
                                            value_type="tpm", sample_ids=sample_ids)
        if data:
            expr_records[eid] = {d["sample_id"]: d["value"] for d in data}

    expr_df = pd.DataFrame(expr_records).apply(pd.to_numeric, errors="coerce") if expr_records else pd.DataFrame()
    print(f"  mRNA matrix: {expr_df.shape[0]} samples x {expr_df.shape[1]} genes")

    # --- Protein abundance matrix (gene by gene) ---
    print(f"  Fetching protein abundance data for {len(entrez_ids)} genes...")
    prot_records = {}
    for eid in entrez_ids:
        data = get_protein_abundance(db, cohort, entrez_id=eid, sample_ids=sample_ids)
        if data:
            prot_records[eid] = {d["sample_id"]: d["abundance"] for d in data}

    prot_df = pd.DataFrame(prot_records).apply(pd.to_numeric, errors="coerce") if prot_records else pd.DataFrame()
    print(f"  Protein matrix: {prot_df.shape[0]} samples x {prot_df.shape[1]} genes")

    return cnv_df, expr_df, prot_df


# =============================================================================
# STEP 2: Compute per-gene cross-layer correlations
# =============================================================================

def compute_cross_layer_correlations(cnv_df: pd.DataFrame,
                                      expr_df: pd.DataFrame,
                                      prot_df: pd.DataFrame,
                                      min_samples: int = MIN_SAMPLES_CORR
                                      ) -> pd.DataFrame:
    """
    For each gene present in all three matrices, compute:
    - corr_cnv_mrna: Pearson r(CNV, mRNA) across samples
    - corr_mrna_prot: Pearson r(mRNA, protein) across samples
    - corr_cnv_prot: Pearson r(CNV, protein) across samples (end-to-end)
    """
    # Find genes present in all three matrices
    common_genes = sorted(
        set(cnv_df.columns) & set(expr_df.columns) & set(prot_df.columns)
    )
    print(f"  Genes with all three layers: {len(common_genes)}")

    rows = []
    for gene in common_genes:
        cnv_vals = cnv_df[gene].dropna()
        expr_vals = expr_df[gene].dropna()
        prot_vals = prot_df[gene].dropna()

        row = {
            "entrez_id": gene,
            "corr_cnv_mrna": np.nan,
            "pval_cnv_mrna": np.nan,
            "corr_mrna_prot": np.nan,
            "pval_mrna_prot": np.nan,
            "corr_cnv_prot": np.nan,
            "pval_cnv_prot": np.nan,
            "n_samples_cnv_mrna": 0,
            "n_samples_mrna_prot": 0,
            "n_samples_cnv_prot": 0,
        }

        # CNV -> mRNA
        common_cm = cnv_vals.index.intersection(expr_vals.index)
        if len(common_cm) >= min_samples:
            r, p = stats.pearsonr(cnv_vals.loc[common_cm], expr_vals.loc[common_cm])
            row["corr_cnv_mrna"] = r
            row["pval_cnv_mrna"] = p
            row["n_samples_cnv_mrna"] = len(common_cm)

        # mRNA -> Protein
        common_mp = expr_vals.index.intersection(prot_vals.index)
        if len(common_mp) >= min_samples:
            r, p = stats.pearsonr(expr_vals.loc[common_mp], prot_vals.loc[common_mp])
            row["corr_mrna_prot"] = r
            row["pval_mrna_prot"] = p
            row["n_samples_mrna_prot"] = len(common_mp)

        # CNV -> Protein (end-to-end)
        common_cp = cnv_vals.index.intersection(prot_vals.index)
        if len(common_cp) >= min_samples:
            r, p = stats.pearsonr(cnv_vals.loc[common_cp], prot_vals.loc[common_cp])
            row["corr_cnv_prot"] = r
            row["pval_cnv_prot"] = p
            row["n_samples_cnv_prot"] = len(common_cp)

        rows.append(row)

    return pd.DataFrame(rows)


def classify_discordance(corr_df: pd.DataFrame,
                          high_thresh: float = CORR_THRESHOLD_HIGH,
                          low_thresh: float = CORR_THRESHOLD_LOW
                          ) -> pd.DataFrame:
    """
    Classify each gene into discordance categories:
    - concordant: CNV-mRNA high AND mRNA-protein high
    - transcriptional_discordant: CNV-mRNA low (compensation at transcription level)
    - post_transcriptional_discordant: CNV-mRNA high but mRNA-protein low
    - mixed_discordant: both transitions show low correlation
    """
    df = corr_df.copy()

    # Drop genes with missing correlations
    df = df.dropna(subset=["corr_cnv_mrna", "corr_mrna_prot"])

    conditions = []
    for _, row in df.iterrows():
        cm = row["corr_cnv_mrna"]
        mp = row["corr_mrna_prot"]

        if cm >= high_thresh and mp >= high_thresh:
            conditions.append("concordant")
        elif cm < low_thresh and mp >= high_thresh:
            conditions.append("transcriptional_discordant")
        elif cm >= high_thresh and mp < low_thresh:
            conditions.append("post_transcriptional_discordant")
        elif cm < low_thresh and mp < low_thresh:
            conditions.append("mixed_discordant")
        else:
            conditions.append("intermediate")

    df["category"] = conditions
    return df


# =============================================================================
# STEP 3: KG semantic enrichment for discordant genes
# =============================================================================

def get_gene_go_terms(db, entrez_id: str, limit: int = 50) -> List[Dict]:
    """
    Query KG for GO terms associated with a gene (via its EntrezID node).
    Traverses: EntrezID -> (any edge) -> GO node
    Also traverses through protein nodes: EntrezID -> PR -> GO
    """
    aql = """
    LET direct_go = (
        FOR v, e IN 1..2 ANY CONCAT("nodes/", @entrezId) edges
            FILTER v.class_code == "GO"
            LIMIT @limit
            RETURN DISTINCT {
                go_id: v.entity_id,
                go_label: v.label,
                go_key: v._key,
                predicate: e.predicate_label,
                bioentity_type: v.bioentity_type
            }
    )
    RETURN direct_go
    """
    result = list(db.aql.execute(aql, bind_vars={
        "entrezId": entrez_id, "limit": limit,
    }))
    return result[0] if result else []


def get_gene_pathways(db, entrez_id: str, limit: int = 50) -> List[Dict]:
    """
    Query KG for Reactome pathways (R-HSA) associated with a gene.
    """
    aql = """
    FOR v, e IN 1..2 ANY CONCAT("nodes/", @entrezId) edges
        FILTER v.class_code == "R-HSA"
        LIMIT @limit
        RETURN DISTINCT {
            pathway_id: v.entity_id,
            pathway_label: v.label,
            pathway_key: v._key,
            predicate: e.predicate_label
        }
    """
    result = list(db.aql.execute(aql, bind_vars={
        "entrezId": entrez_id, "limit": limit,
    }))
    return result


def get_gene_interaction_predicates(db, entrez_id: str, limit: int = 100) -> List[Dict]:
    """
    Get all predicates (edge types) connected to a gene node.
    """
    aql = """
    FOR v, e IN 1..1 ANY CONCAT("nodes/", @entrezId) edges
        LIMIT @limit
        RETURN {
            predicate: e.predicate_label,
            neighbor_class: v.class_code,
            neighbor_label: v.label
        }
    """
    result = list(db.aql.execute(aql, bind_vars={
        "entrezId": entrez_id, "limit": limit,
    }))
    return result


def enrich_genes_with_kg(db, gene_ids: List[str],
                          label: str = "group") -> Dict:
    """
    For a list of entrez_ids, collect GO terms, pathways, and interaction predicates.
    Returns aggregated counts (number of UNIQUE GENES having each term, not
    occurrence counts) and per-gene details.
    """
    # Counters track how many genes have each term (not how many times it appears)
    go_counter = Counter()
    pathway_counter = Counter()
    predicate_counter = Counter()
    per_gene = {}

    for eid in gene_ids:
        gene_data = {"go_terms": [], "pathways": [], "predicates": []}

        # GO terms — deduplicate per gene
        go_terms = get_gene_go_terms(db, eid)
        gene_go_labels = set()
        for g in go_terms:
            go_label = g.get("go_label", "unknown")
            gene_go_labels.add(go_label)
            gene_data["go_terms"].append(go_label)
        for go_label in gene_go_labels:
            go_counter[go_label] += 1

        # Pathways — deduplicate per gene
        pathways = get_gene_pathways(db, eid)
        gene_pw_labels = set()
        for p in pathways:
            pw_label = p.get("pathway_label", "unknown")
            gene_pw_labels.add(pw_label)
            gene_data["pathways"].append(pw_label)
        for pw_label in gene_pw_labels:
            pathway_counter[pw_label] += 1

        # Interaction predicates — deduplicate per gene
        interactions = get_gene_interaction_predicates(db, eid)
        gene_pred_labels = set()
        for i in interactions:
            pred = i.get("predicate", "unknown")
            gene_pred_labels.add(pred)
            gene_data["predicates"].append(pred)
        for pred in gene_pred_labels:
            predicate_counter[pred] += 1

        per_gene[eid] = gene_data

    return {
        "label": label,
        "n_genes": len(gene_ids),
        "go_term_counts": go_counter,
        "pathway_counts": pathway_counter,
        "predicate_counts": predicate_counter,
        "per_gene": per_gene,
    }


def compare_enrichment(discordant_enrichment: Dict,
                        concordant_enrichment: Dict,
                        term_type: str = "go_term_counts",
                        min_count: int = 3) -> pd.DataFrame:
    """
    Compare term frequencies between discordant and concordant gene sets.
    Compute fold enrichment and Fisher's exact test p-value.
    """
    disc_counts = discordant_enrichment[term_type]
    conc_counts = concordant_enrichment[term_type]
    n_disc = discordant_enrichment["n_genes"]
    n_conc = concordant_enrichment["n_genes"]

    all_terms = set(disc_counts.keys()) | set(conc_counts.keys())

    rows = []
    for term in all_terms:
        d = disc_counts.get(term, 0)
        c = conc_counts.get(term, 0)

        if d + c < min_count:
            continue

        # Fisher's exact test: is 'term' enriched in discordant vs concordant?
        # Contingency table (gene counts, not occurrence counts):
        #                  has_term   no_term
        # discordant         d        n_disc - d
        # concordant         c        n_conc - c
        # Clamp: d/c should not exceed n_disc/n_conc after dedup fix,
        # but guard against edge cases
        d_clamped = min(d, n_disc)
        c_clamped = min(c, n_conc)
        table = [[d_clamped, n_disc - d_clamped],
                  [c_clamped, n_conc - c_clamped]]
        try:
            odds, pval = fisher_exact(table, alternative="greater")
        except ValueError:
            odds, pval = np.nan, np.nan

        # Fold enrichment (frequency in discordant / frequency in concordant)
        freq_disc = d / n_disc if n_disc > 0 else 0
        freq_conc = c / n_conc if n_conc > 0 else 0
        fold = freq_disc / freq_conc if freq_conc > 0 else np.inf

        rows.append({
            "term": term,
            "count_discordant": d,
            "count_concordant": c,
            "freq_discordant": freq_disc,
            "freq_concordant": freq_conc,
            "fold_enrichment": fold,
            "fisher_pval": pval,
            "odds_ratio": odds,
        })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values("fisher_pval")
    return df


# =============================================================================
# STEP 4: VISUALIZATIONS
# =============================================================================

def plot_discordance_scatter(corr_df: pd.DataFrame, output_path: str = None):
    """
    Scatter plot: CNV-mRNA correlation (x) vs mRNA-protein correlation (y).
    Points colored by discordance category, with quadrant shading.
    """
    df = corr_df.dropna(subset=["corr_cnv_mrna", "corr_mrna_prot"]).copy()
    if df.empty:
        print("  No data for scatter plot, skipping.")
        return

    fig, ax = plt.subplots(figsize=(12, 10))

    # Quadrant shading
    ax.axhline(y=CORR_THRESHOLD_HIGH, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axvline(x=CORR_THRESHOLD_HIGH, color="grey", linestyle="--", linewidth=0.8, alpha=0.5)
    ax.axhline(y=CORR_THRESHOLD_LOW, color="grey", linestyle=":", linewidth=0.6, alpha=0.3)
    ax.axvline(x=CORR_THRESHOLD_LOW, color="grey", linestyle=":", linewidth=0.6, alpha=0.3)

    # Light background quadrants
    xlim = (-0.6, 1.0)
    ylim = (-0.6, 1.0)
    # Concordant quadrant (top-right): green
    ax.fill_between([CORR_THRESHOLD_HIGH, xlim[1]], CORR_THRESHOLD_HIGH, ylim[1],
                     alpha=0.06, color="green")
    # Post-transcriptional discordant (bottom-right): red
    ax.fill_between([CORR_THRESHOLD_HIGH, xlim[1]], ylim[0], CORR_THRESHOLD_LOW,
                     alpha=0.06, color="red")
    # Transcriptional discordant (top-left): orange
    ax.fill_between([xlim[0], CORR_THRESHOLD_LOW], CORR_THRESHOLD_HIGH, ylim[1],
                     alpha=0.06, color="orange")
    # Mixed discordant (bottom-left): purple
    ax.fill_between([xlim[0], CORR_THRESHOLD_LOW], ylim[0], CORR_THRESHOLD_LOW,
                     alpha=0.06, color="purple")

    # Color mapping
    cat_colors = {
        "concordant": "#2ca02c",
        "transcriptional_discordant": "#ff7f0e",
        "post_transcriptional_discordant": "#d62728",
        "mixed_discordant": "#9467bd",
        "intermediate": "#BBBBBB",
    }
    cat_labels = {
        "concordant": "Concordant",
        "transcriptional_discordant": "Transcriptional discordant",
        "post_transcriptional_discordant": "Post-transcriptional discordant",
        "mixed_discordant": "Mixed discordant",
        "intermediate": "Intermediate",
    }

    for cat, color in cat_colors.items():
        subset = df[df["category"] == cat]
        if subset.empty:
            continue
        ax.scatter(subset["corr_cnv_mrna"], subset["corr_mrna_prot"],
                   c=color, label=f"{cat_labels[cat]} (n={len(subset)})",
                   s=40, alpha=0.7, edgecolors="white", linewidths=0.3)

    # Label notable genes (most extreme discordant)
    for cat in ["post_transcriptional_discordant", "transcriptional_discordant", "mixed_discordant"]:
        subset = df[df["category"] == cat]
        if subset.empty:
            continue
        # Label top 5 most extreme by distance from diagonal
        subset = subset.copy()
        subset["dist"] = abs(subset["corr_cnv_mrna"] - subset["corr_mrna_prot"])
        top = subset.nlargest(5, "dist")
        for _, row in top.iterrows():
            label = row.get("gene_symbol", row["entrez_id"])
            ax.annotate(label, (row["corr_cnv_mrna"], row["corr_mrna_prot"]),
                        fontsize=7, alpha=0.8, ha="left",
                        xytext=(4, 4), textcoords="offset points")

    ax.set_xlabel("Pearson r (CNV → mRNA)", fontsize=13)
    ax.set_ylabel("Pearson r (mRNA → Protein)", fontsize=13)
    ax.set_title("UC2: Cross-Layer Discordance — CNV vs mRNA vs Protein (TCGA-BRCA)",
                 fontsize=14, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9, framealpha=0.9)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    # Diagonal reference
    ax.plot(xlim, ylim, "k--", alpha=0.2, linewidth=0.8)

    # Quadrant labels
    ax.text(0.85, 0.85, "CONCORDANT", transform=ax.transAxes,
            fontsize=9, color="green", alpha=0.5, ha="center", fontweight="bold")
    ax.text(0.15, 0.85, "TRANSCRIPTIONAL\nDISCORDANCE", transform=ax.transAxes,
            fontsize=8, color="orange", alpha=0.5, ha="center", fontweight="bold")
    ax.text(0.85, 0.12, "POST-TRANSCRIPTIONAL\nDISCORDANCE", transform=ax.transAxes,
            fontsize=8, color="red", alpha=0.5, ha="center", fontweight="bold")
    ax.text(0.15, 0.12, "MIXED\nDISCORDANCE", transform=ax.transAxes,
            fontsize=8, color="purple", alpha=0.5, ha="center", fontweight="bold")

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def plot_discordant_network(db, discordant_genes: pd.DataFrame,
                             enrichment_data: Dict,
                             top_n_genes: int = 40,
                             output_path: str = None):
    """
    Network: discordant genes in KG, colored by most frequent GO process type.
    Edges from KG (1-hop) between the discordant genes.
    """
    if discordant_genes.empty:
        print("  No discordant genes for network plot, skipping.")
        return

    genes = discordant_genes.head(top_n_genes)
    gene_ids = genes["entrez_id"].tolist()

    # Assign dominant GO term category per gene
    per_gene = enrichment_data.get("per_gene", {})

    # Simplify GO labels to broad categories
    def categorize_go(go_labels):
        categories = {
            "translation": ["translation", "ribosom", "peptide biosynthetic"],
            "protein modification": ["ubiquitin", "phosphoryl", "acetyl", "sumoyl",
                                      "protein modification", "post-translational"],
            "proteolysis": ["proteolysis", "proteasom", "peptidase", "degradation"],
            "transcription regulation": ["transcription", "RNA polymerase",
                                          "chromatin", "histone"],
            "signal transduction": ["signal", "kinase activity", "receptor",
                                     "GTPase", "MAPK"],
            "cell cycle": ["cell cycle", "mitotic", "cell division", "checkpoint"],
            "apoptosis": ["apoptot", "programmed cell death", "caspase"],
            "metabolism": ["metabol", "biosynthe", "catabol", "oxidation-reduction"],
        }
        for go_label in go_labels:
            go_lower = go_label.lower()
            for cat, keywords in categories.items():
                if any(kw in go_lower for kw in keywords):
                    return cat
        return "other"

    gene_categories = {}
    for eid in gene_ids:
        go_terms = per_gene.get(eid, {}).get("go_terms", [])
        gene_categories[eid] = categorize_go(go_terms)

    # Build network: query KG edges between these genes
    # Use entrez_ids as _keys (EntrezID nodes have _key = entrez_id)
    node_keys = gene_ids
    try:
        aql = """
        FOR nk IN @nodeKeys
            FOR v, e IN 1..1 ANY CONCAT("nodes/", nk) edges
                FILTER v._key IN @nodeKeys
                RETURN DISTINCT {
                    source: nk,
                    target: v._key,
                    predicate: e.predicate_label
                }
        """
        edges = list(db.aql.execute(aql, bind_vars={"nodeKeys": node_keys}))
    except Exception as exc:
        print(f"  Warning: KG edge query failed ({exc}), building node-only network.")
        edges = []

    G = nx.Graph()

    # Add nodes
    for _, row in genes.iterrows():
        eid = row["entrez_id"]
        label = row.get("gene_symbol", eid)
        G.add_node(eid, label=label, category=gene_categories.get(eid, "other"),
                   corr_cnv_mrna=row.get("corr_cnv_mrna", 0),
                   corr_mrna_prot=row.get("corr_mrna_prot", 0))

    # Add edges
    seen = set()
    for e in edges:
        key = tuple(sorted([e["source"], e["target"]]))
        if key not in seen and e["source"] != e["target"]:
            seen.add(key)
            G.add_edge(e["source"], e["target"], predicate=e["predicate"])

    if len(G.nodes) == 0:
        print("  Empty network, skipping.")
        return

    fig, ax = plt.subplots(figsize=(14, 12))

    # Category colors
    cat_palette = {
        "translation": "#e41a1c",
        "protein modification": "#377eb8",
        "proteolysis": "#4daf4a",
        "transcription regulation": "#984ea3",
        "signal transduction": "#ff7f00",
        "cell cycle": "#a65628",
        "apoptosis": "#f781bf",
        "metabolism": "#999999",
        "other": "#CCCCCC",
    }

    pos = nx.spring_layout(G, seed=42, k=2.5 / max(np.sqrt(len(G.nodes)), 1), iterations=80)

    # Node colors by GO category
    node_colors = [cat_palette.get(G.nodes[n].get("category", "other"), "#CCCCCC")
                   for n in G.nodes()]
    node_sizes = [200 + G.degree(n) * 60 for n in G.nodes()]

    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, edge_color="grey", width=0.8)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85,
                           edgecolors="white", linewidths=0.5)

    # Labels
    labels = {n: G.nodes[n].get("label", n) for n in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=7, font_weight="bold", ax=ax)

    # Legend
    present_cats = set(gene_categories.values())
    legend_patches = [mpatches.Patch(color=cat_palette.get(c, "#CCCCCC"), label=c)
                      for c in sorted(present_cats) if c in cat_palette]
    ax.legend(handles=legend_patches, loc="upper left", fontsize=8,
              title="Dominant GO Process", title_fontsize=9)

    ax.set_title(f"UC2: Discordant Genes in KG — Colored by GO Process\n"
                 f"({len(G.nodes)} genes, {len(G.edges)} KG edges)",
                 fontsize=13, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


def plot_layer_heatmap(corr_df: pd.DataFrame, gene_info: Dict[str, str],
                        top_n: int = 50, output_path: str = None):
    """
    Heatmap: genes x omic correlation metrics, grouped by discordance category.
    Shows the three correlation values (CNV-mRNA, mRNA-Prot, CNV-Prot) per gene.
    """
    df = corr_df.dropna(subset=["corr_cnv_mrna", "corr_mrna_prot"]).copy()
    if df.empty:
        print("  No data for heatmap, skipping.")
        return

    # Take top genes per category by extremeness
    selected = []
    for cat in ["concordant", "post_transcriptional_discordant",
                "transcriptional_discordant", "mixed_discordant"]:
        subset = df[df["category"] == cat].copy()
        if subset.empty:
            continue
        # Sort by how "characteristic" they are of the category
        if cat == "concordant":
            subset["score"] = subset["corr_cnv_mrna"] + subset["corr_mrna_prot"]
        elif cat == "post_transcriptional_discordant":
            subset["score"] = subset["corr_cnv_mrna"] - subset["corr_mrna_prot"]
        elif cat == "transcriptional_discordant":
            subset["score"] = subset["corr_mrna_prot"] - subset["corr_cnv_mrna"]
        else:
            subset["score"] = -(subset["corr_cnv_mrna"] + subset["corr_mrna_prot"])
        selected.append(subset.nlargest(min(top_n // 4, len(subset)), "score"))

    if not selected:
        return

    plot_df = pd.concat(selected)

    # Build heatmap matrix
    metrics = ["corr_cnv_mrna", "corr_mrna_prot", "corr_cnv_prot"]
    metric_labels = ["CNV → mRNA", "mRNA → Protein", "CNV → Protein"]

    # Use gene symbols as row labels where available
    plot_df = plot_df.reset_index(drop=True)
    labels = []
    seen = set()
    for _, row in plot_df.iterrows():
        sym = gene_info.get(row["entrez_id"], row["entrez_id"])
        # Deduplicate labels
        if sym in seen:
            sym = f"{sym} ({row['entrez_id']})"
        seen.add(sym)
        labels.append(sym)
    plot_df["gene_label"] = labels

    hm_data = plot_df[metrics].copy()
    hm_data.index = plot_df["gene_label"]
    hm_data.columns = metric_labels

    # Row colors by category
    cat_colors_map = {
        "concordant": "#2ca02c",
        "transcriptional_discordant": "#ff7f0e",
        "post_transcriptional_discordant": "#d62728",
        "mixed_discordant": "#9467bd",
        "intermediate": "#BBBBBB",
    }
    row_colors = pd.Series(
        [cat_colors_map.get(c, "#BBBBBB") for c in plot_df["category"]],
        index=hm_data.index,
        name="Category"
    )

    fig_height = max(8, len(hm_data) * 0.25)
    try:
        g = sns.clustermap(hm_data, row_cluster=False, col_cluster=False,
                            cmap="RdBu_r", center=0, vmin=-0.5, vmax=0.8,
                            figsize=(8, fig_height), linewidths=0.3,
                            row_colors=row_colors, annot=True, fmt=".2f",
                            annot_kws={"size": 7},
                            cbar_kws={"label": "Pearson r", "shrink": 0.5})

        g.ax_heatmap.set_title("UC2: Cross-Layer Correlations by Discordance Category",
                                fontsize=13, fontweight="bold", pad=20)
        g.ax_heatmap.set_ylabel("")

        # Add category legend
        legend_patches = [mpatches.Patch(color=c, label=cat.replace("_", " ").title())
                          for cat, c in cat_colors_map.items()
                          if cat in plot_df["category"].values]
        g.ax_heatmap.legend(handles=legend_patches, loc="upper right",
                             bbox_to_anchor=(1.3, 1.0), fontsize=8, title="Category")

        if output_path:
            g.savefig(output_path, dpi=200, bbox_inches="tight")
            print(f"  Saved: {output_path}")
        plt.show()
        return g
    except Exception as exc:
        print(f"  Warning: clustermap failed ({exc}), falling back to simple heatmap.")
        fig, ax = plt.subplots(figsize=(8, fig_height))
        sns.heatmap(hm_data, cmap="RdBu_r", center=0, vmin=-0.5, vmax=0.8,
                    annot=True, fmt=".2f", ax=ax, linewidths=0.3)
        ax.set_title("UC2: Cross-Layer Correlations by Discordance Category",
                      fontsize=13, fontweight="bold")
        plt.tight_layout()
        if output_path:
            fig.savefig(output_path, dpi=200, bbox_inches="tight")
            print(f"  Saved: {output_path}")
        plt.show()
        return fig


def plot_publication_figure(corr_df: pd.DataFrame,
                             go_enrichment_df: pd.DataFrame,
                             pred_enrichment_df: pd.DataFrame,
                             output_path: str = None):
    """
    Composite publication figure for UC2 (Bioinformatics journal style).

    Layout (double-column, ~180 mm wide):
        ┌──────────────────────────────┬─────────────────────────────┐
        │                              │  B  GO Term Enrichment      │
        │  A  Cross-layer scatter      │                             │
        │                              ├─────────────────────────────┤
        │                              │  C  KG Predicate Enrichment │
        └──────────────────────────────┴─────────────────────────────┘

    Panel A (left, ~60%): scatter CNV-mRNA vs mRNA-Protein (main result).
    Panel B (top-right):  top GO term enrichment (discordant vs concordant).
    Panel C (bottom-right): top KG predicate enrichment.
    """
    df = corr_df.dropna(subset=["corr_cnv_mrna", "corr_mrna_prot"]).copy()
    if df.empty:
        print("  No data for publication figure, skipping.")
        return

    # Bioinformatics-style typography: sans-serif, small sizes, thin lines
    rc_backup = {k: plt.rcParams[k] for k in [
        "font.family", "font.size", "axes.labelsize", "axes.titlesize",
        "xtick.labelsize", "ytick.labelsize", "legend.fontsize",
        "axes.linewidth", "xtick.major.width", "ytick.major.width",
    ]}
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 6.5,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
    })

    # Figure size: ~180 mm wide (double column) x ~135 mm tall
    fig = plt.figure(figsize=(7.2, 5.4))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.45, 1.0],
        height_ratios=[1.0, 1.0],
        wspace=0.55, hspace=0.75,
        left=0.07, right=0.985, top=0.90, bottom=0.09,
    )
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])

    # =========================================================================
    # Panel A — Scatter cross-layer
    # =========================================================================
    xlim = (-0.6, 1.0)
    ylim = (-0.6, 1.0)

    cat_colors = {
        "concordant": "#2ca02c",
        "transcriptional_discordant": "#ff7f0e",
        "post_transcriptional_discordant": "#d62728",
        "mixed_discordant": "#9467bd",
        "intermediate": "#BBBBBB",
    }
    cat_labels = {
        "concordant": "Concordant",
        "transcriptional_discordant": "Transcr. discordant",
        "post_transcriptional_discordant": "Post-transcr. discordant",
        "mixed_discordant": "Mixed discordant",
        "intermediate": "Intermediate",
    }

    ax_a.fill_between([CORR_THRESHOLD_HIGH, xlim[1]], CORR_THRESHOLD_HIGH, ylim[1],
                      alpha=0.05, color="green", linewidth=0)
    ax_a.fill_between([CORR_THRESHOLD_HIGH, xlim[1]], ylim[0], CORR_THRESHOLD_LOW,
                      alpha=0.05, color="red", linewidth=0)
    ax_a.fill_between([xlim[0], CORR_THRESHOLD_LOW], CORR_THRESHOLD_HIGH, ylim[1],
                      alpha=0.05, color="orange", linewidth=0)
    ax_a.fill_between([xlim[0], CORR_THRESHOLD_LOW], ylim[0], CORR_THRESHOLD_LOW,
                      alpha=0.05, color="purple", linewidth=0)

    ax_a.axhline(y=CORR_THRESHOLD_HIGH, color="grey", linestyle="--", linewidth=0.5, alpha=0.5)
    ax_a.axvline(x=CORR_THRESHOLD_HIGH, color="grey", linestyle="--", linewidth=0.5, alpha=0.5)
    ax_a.axhline(y=CORR_THRESHOLD_LOW, color="grey", linestyle=":", linewidth=0.4, alpha=0.4)
    ax_a.axvline(x=CORR_THRESHOLD_LOW, color="grey", linestyle=":", linewidth=0.4, alpha=0.4)
    ax_a.plot(xlim, ylim, "k--", alpha=0.2, linewidth=0.5)

    for cat, color in cat_colors.items():
        subset = df[df["category"] == cat]
        if subset.empty:
            continue
        ax_a.scatter(subset["corr_cnv_mrna"], subset["corr_mrna_prot"],
                     c=color, label=f"{cat_labels[cat]} (n={len(subset)})",
                     s=14, alpha=0.75, edgecolors="white", linewidths=0.25)

    # Label extreme discordant genes — alternate offsets to avoid overlap
    labelled = set()
    annotations = []  # (x, y, label)
    for cat in ["post_transcriptional_discordant", "transcriptional_discordant", "mixed_discordant"]:
        subset = df[df["category"] == cat].copy()
        if subset.empty:
            continue
        subset["dist"] = abs(subset["corr_cnv_mrna"] - subset["corr_mrna_prot"])
        top = subset.nlargest(2, "dist")
        for _, row in top.iterrows():
            label = row.get("gene_symbol", row["entrez_id"])
            if not isinstance(label, str) or not label or label in labelled:
                continue
            labelled.add(label)
            annotations.append((row["corr_cnv_mrna"], row["corr_mrna_prot"], label))

    # Place labels with offsets that alternate to reduce collisions
    offsets = [(5, 4), (5, -8), (-5, 5), (-5, -8), (6, 9), (6, -10)]
    for i, (x, y, label) in enumerate(annotations):
        dx, dy = offsets[i % len(offsets)]
        ha = "left" if dx >= 0 else "right"
        ax_a.annotate(label, (x, y),
                      fontsize=6, alpha=0.95, ha=ha, fontstyle="italic",
                      xytext=(dx, dy), textcoords="offset points",
                      arrowprops=dict(arrowstyle="-", lw=0.3, color="grey", alpha=0.5))

    ax_a.set_xlabel("Pearson r (CNV → mRNA)")
    ax_a.set_ylabel("Pearson r (mRNA → Protein)")
    ax_a.set_xlim(xlim)
    ax_a.set_ylim(ylim)
    ax_a.set_xticks(np.arange(-0.6, 1.01, 0.2))
    ax_a.set_yticks(np.arange(-0.6, 1.01, 0.2))
    ax_a.tick_params(direction="out", length=3)
    for spine in ["top", "right"]:
        ax_a.spines[spine].set_visible(False)

    # Quadrant annotations — corner-anchored to avoid covering data
    ax_a.text(0.985, 0.985, "Concordant", transform=ax_a.transAxes,
              fontsize=6.5, color="#2ca02c", alpha=0.85, ha="right",
              va="top", fontweight="bold")
    ax_a.text(0.015, 0.985, "Transcr.\ndiscordance", transform=ax_a.transAxes,
              fontsize=6.5, color="#ff7f0e", alpha=0.85, ha="left",
              va="top", fontweight="bold")
    ax_a.text(0.985, 0.015, "Post-transcr.\ndiscordance", transform=ax_a.transAxes,
              fontsize=6.5, color="#d62728", alpha=0.85, ha="right",
              va="bottom", fontweight="bold")
    ax_a.text(0.015, 0.015, "Mixed\ndiscordance", transform=ax_a.transAxes,
              fontsize=6.5, color="#9467bd", alpha=0.85, ha="left",
              va="bottom", fontweight="bold")

    # Legend below the panel, horizontal, no overlap with data
    ax_a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.10),
                fontsize=6, framealpha=0.0, ncol=3,
                handletextpad=0.3, borderpad=0.2,
                columnspacing=0.9, labelspacing=0.3)

    # =========================================================================
    # Helper for B/C — horizontal bar with frequency + significance markers
    # =========================================================================
    def _enrichment_panel(ax, enrichment_df, top_n, ylabel_prefix,
                          max_term_chars=42, show_legend=False):
        if enrichment_df is None or enrichment_df.empty:
            ax.text(0.5, 0.5, "No enrichment data",
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=7, color="grey")
            ax.set_xticks([]); ax.set_yticks([])
            return None

        df_e = enrichment_df.copy()
        # Pick significant first, otherwise top by p-value
        sig = df_e[df_e["fisher_pval"] < 0.05]
        if len(sig) >= 3:
            sel = sig.nsmallest(top_n, "fisher_pval")
        else:
            sel = df_e.nsmallest(top_n, "fisher_pval")

        # Order top -> bottom by significance
        sel = sel.sort_values("fisher_pval", ascending=True).head(top_n)

        terms = sel["term"].astype(str).tolist()
        terms_short = [t if len(t) <= max_term_chars else t[:max_term_chars - 1] + "…"
                       for t in terms]
        freq_disc = sel["freq_discordant"].values
        freq_conc = sel["freq_concordant"].values
        pvals = sel["fisher_pval"].values

        y = np.arange(len(sel))
        h = 0.38
        # Two grouped horizontal bars per term
        b_disc = ax.barh(y - h/2, freq_disc, height=h, color="#d62728",
                         alpha=0.85, edgecolor="white", linewidth=0.3,
                         label="Discordant")
        b_conc = ax.barh(y + h/2, freq_conc, height=h, color="#2ca02c",
                         alpha=0.85, edgecolor="white", linewidth=0.3,
                         label="Concordant")

        # Significance markers at end of discordant bar
        xmax = max(max(freq_disc) if len(freq_disc) else 0,
                   max(freq_conc) if len(freq_conc) else 0)
        for yi, p, fd in zip(y, pvals, freq_disc):
            if p < 0.001:
                star = "***"
            elif p < 0.01:
                star = "**"
            elif p < 0.05:
                star = "*"
            else:
                star = ""
            if star:
                ax.text(fd + xmax * 0.015, yi - h/2, star,
                        va="center", ha="left", fontsize=7,
                        fontweight="bold", color="#333333")

        ax.set_yticks(y)
        ax.set_yticklabels(terms_short, fontsize=6.5)
        ax.invert_yaxis()
        ax.set_xlabel("Gene frequency", fontsize=7.5)
        ax.set_xlim(0, xmax * 1.18 if xmax > 0 else 1)
        ax.tick_params(direction="out", length=2.5)
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        if show_legend:
            ax.legend(loc="lower right", fontsize=6, framealpha=0.92,
                      handletextpad=0.4, borderpad=0.3, labelspacing=0.25,
                      frameon=True, edgecolor="lightgrey")
        return (b_disc, b_conc)

    # =========================================================================
    # Panel B — GO term enrichment
    # =========================================================================
    bars_b = _enrichment_panel(ax_b, go_enrichment_df, top_n=6,
                                ylabel_prefix="GO", max_term_chars=42,
                                show_legend=False)
    ax_b.set_title("Top GO biological processes", fontsize=8, pad=6)

    # =========================================================================
    # Panel C — KG predicate enrichment (carries shared legend for B+C)
    # =========================================================================
    bars_c = _enrichment_panel(ax_c, pred_enrichment_df, top_n=6,
                                ylabel_prefix="Predicate", max_term_chars=42,
                                show_legend=True)
    ax_c.set_title("Top KG predicates", fontsize=8, pad=6)

    # =========================================================================
    # Panel labels (A, B, C) — Bioinformatics style: bold, top-left corner
    # =========================================================================
    ax_a.text(-0.05, 1.02, "A", transform=ax_a.transAxes,
              fontsize=12, fontweight="bold", va="bottom", ha="left")
    ax_b.text(-0.32, 1.08, "B", transform=ax_b.transAxes,
              fontsize=12, fontweight="bold", va="bottom", ha="left")
    ax_c.text(-0.32, 1.08, "C", transform=ax_c.transAxes,
              fontsize=12, fontweight="bold", va="bottom", ha="left")

    # Global suptitle — concise, above panel labels
    fig.suptitle(
        "Cross-layer discordance and KG-driven semantic enrichment in TCGA-BRCA",
        fontsize=9, fontweight="bold", y=0.998,
    )

    if output_path:
        # Save both PNG (300 dpi) and PDF (vector) for publication
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        pdf_path = os.path.splitext(output_path)[0] + ".pdf"
        try:
            fig.savefig(pdf_path, bbox_inches="tight")
            print(f"  Saved: {output_path}")
            print(f"  Saved: {pdf_path}")
        except Exception as exc:
            print(f"  Saved: {output_path}  (PDF export skipped: {exc})")

    plt.show()
    plt.rcParams.update(rc_backup)
    return fig


def plot_enrichment_bars(enrichment_df: pd.DataFrame,
                          term_type_label: str = "GO Term",
                          top_n: int = 20,
                          output_path: str = None):
    """
    Bar chart of top enriched terms in discordant vs concordant genes.
    """
    if enrichment_df.empty:
        print(f"  No enrichment data for {term_type_label}, skipping.")
        return

    # Filter significant and take top by fold enrichment
    df = enrichment_df[enrichment_df["fisher_pval"] < 0.1].copy()
    if df.empty:
        df = enrichment_df.head(top_n)
    else:
        df = df.nsmallest(top_n, "fisher_pval")

    # Truncate long term names
    df["term_short"] = df["term"].apply(
        lambda x: (x[:50] + "...") if len(str(x)) > 53 else str(x)
    )

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, max(6, len(df) * 0.35)),
                                     gridspec_kw={"width_ratios": [2, 1]})

    # Left: frequency comparison
    y_pos = range(len(df))
    ax1.barh(y_pos, df["freq_discordant"], height=0.4, align="center",
             color="#d62728", alpha=0.8, label="Discordant")
    ax1.barh([y + 0.4 for y in y_pos], df["freq_concordant"], height=0.4,
             align="center", color="#2ca02c", alpha=0.8, label="Concordant")
    ax1.set_yticks([y + 0.2 for y in y_pos])
    ax1.set_yticklabels(df["term_short"], fontsize=8)
    ax1.set_xlabel("Frequency (genes with term / total genes)", fontsize=10)
    ax1.set_title(f"Top {term_type_label} Enrichment:\nDiscordant vs Concordant", fontsize=12)
    ax1.legend(fontsize=9)
    ax1.invert_yaxis()

    # Right: -log10(p-value)
    df["neg_log_p"] = -np.log10(df["fisher_pval"].clip(lower=1e-20))
    ax2.barh(y_pos, df["neg_log_p"], color="#4A90D9", alpha=0.8)
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels([])
    ax2.set_xlabel("-log10(p-value)", fontsize=10)
    ax2.set_title("Statistical Significance", fontsize=12)
    ax2.axvline(x=-np.log10(0.05), color="red", linestyle="--", linewidth=0.8,
                label="p=0.05")
    ax2.legend(fontsize=8)
    ax2.invert_yaxis()

    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_uc2(db_name: str = DB_NAME, cohort: str = COHORT):
    """
    Execute the full UC2 analysis pipeline.
    """
    print("=" * 70)
    print("  UC2: CNV-Expression-Protein Discordance with Semantic Explanation")
    print("=" * 70)

    # --- Connect to ArangoDB ---
    db = setup_arangodb_connection(db_name)
    if db is None:
        raise RuntimeError("Failed to connect to ArangoDB")

    # --- Get samples with complete omics (CNV + expression + protein) ---
    print("\n[1/7] Finding samples with complete tri-layer omics...")
    all_sample_ids = list_samples_with_complete_omics(
        db, cohort=cohort,
        omic_types=["gene_expression", "cnv", "protein"],
        specimen_type=None,
    )
    print(f"  Found {len(all_sample_ids)} samples with CNV + expression + protein")

    # Filter to Primary Tumor
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

    # --- Get protein index genes ---
    print("\n[2/7] Retrieving protein (RPPA) index genes...")
    protein_genes = get_protein_index_genes(db, cohort)
    if not protein_genes:
        print("  ERROR: No protein index genes found")
        return

    # Build gene_info lookup: entrez_id -> gene_symbol
    gene_info = {}
    for m in protein_genes:
        eid = m.get("entrez_id")
        sym = m.get("gene_symbol", "")
        if eid:
            gene_info[eid] = sym if sym else eid

    # --- Build multi-omic matrices ---
    print("\n[3/7] Building cross-layer matrices (CNV, mRNA, Protein)...")
    cnv_df, expr_df, prot_df = build_cross_layer_matrices(
        db, cohort, protein_genes, sample_ids
    )

    # --- Compute per-gene cross-layer correlations ---
    print("\n[4/7] Computing per-gene cross-layer correlations...")
    corr_df = compute_cross_layer_correlations(cnv_df, expr_df, prot_df)
    print(f"  Computed correlations for {len(corr_df)} genes")

    # Add gene symbols
    corr_df["gene_symbol"] = corr_df["entrez_id"].map(gene_info)

    # Classify discordance
    corr_df = classify_discordance(corr_df)

    # Summary
    cat_counts = corr_df["category"].value_counts()
    print("\n  Discordance classification:")
    for cat, count in cat_counts.items():
        print(f"    {cat}: {count}")

    # Save correlation results
    corr_df.to_csv(os.path.join(OUTPUT_DIR, "uc2_correlation_results.csv"), index=False)
    print(f"\n  Saved: {OUTPUT_DIR}/uc2_correlation_results.csv")

    # --- Summary statistics ---
    print("\n[5/7] Computing summary statistics...")
    summary_rows = []
    for cat in corr_df["category"].unique():
        subset = corr_df[corr_df["category"] == cat]
        summary_rows.append({
            "category": cat,
            "n_genes": len(subset),
            "median_cnv_mrna": subset["corr_cnv_mrna"].median(),
            "mean_cnv_mrna": subset["corr_cnv_mrna"].mean(),
            "median_mrna_prot": subset["corr_mrna_prot"].median(),
            "mean_mrna_prot": subset["corr_mrna_prot"].mean(),
            "median_cnv_prot": subset["corr_cnv_prot"].median(),
            "mean_cnv_prot": subset["corr_cnv_prot"].mean(),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(os.path.join(OUTPUT_DIR, "uc2_summary_stats.csv"), index=False)
    print(summary_df.to_string(index=False))

    # Mann-Whitney: discordant categories vs concordant
    concordant_vals_cm = corr_df[corr_df["category"] == "concordant"]["corr_cnv_mrna"].dropna()
    concordant_vals_mp = corr_df[corr_df["category"] == "concordant"]["corr_mrna_prot"].dropna()
    print(f"\n  Statistical tests (discordant vs concordant):")
    print(f"  {'Category':<40} {'Metric':<15} {'U stat':>10} {'p-value':>12}")
    print(f"  {'-'*40} {'-'*15} {'-'*10} {'-'*12}")
    for cat in ["transcriptional_discordant", "post_transcriptional_discordant", "mixed_discordant"]:
        subset = corr_df[corr_df["category"] == cat]
        if len(subset) < 5:
            continue
        vals_cm = subset["corr_cnv_mrna"].dropna()
        vals_mp = subset["corr_mrna_prot"].dropna()
        if len(vals_cm) >= 5 and len(concordant_vals_cm) >= 5:
            u, p = stats.mannwhitneyu(vals_cm, concordant_vals_cm, alternative="less")
            print(f"  {cat:<40} {'CNV-mRNA':<15} {u:>10.0f} {p:>12.2e}")
        if len(vals_mp) >= 5 and len(concordant_vals_mp) >= 5:
            u, p = stats.mannwhitneyu(vals_mp, concordant_vals_mp, alternative="less")
            print(f"  {'':<40} {'mRNA-Prot':<15} {u:>10.0f} {p:>12.2e}")

    # --- KG semantic enrichment ---
    print("\n[6/7] Querying KG for semantic enrichment of discordant genes...")

    # Combine all discordant categories
    discordant_ids = corr_df[
        corr_df["category"].isin([
            "transcriptional_discordant",
            "post_transcriptional_discordant",
            "mixed_discordant",
        ])
    ]["entrez_id"].tolist()

    concordant_ids = corr_df[
        corr_df["category"] == "concordant"
    ]["entrez_id"].tolist()

    print(f"  Discordant genes: {len(discordant_ids)}")
    print(f"  Concordant genes: {len(concordant_ids)}")

    disc_enrichment = enrich_genes_with_kg(db, discordant_ids, label="discordant")
    conc_enrichment = enrich_genes_with_kg(db, concordant_ids, label="concordant")

    # Compare GO term enrichment
    go_enrichment = compare_enrichment(disc_enrichment, conc_enrichment,
                                        term_type="go_term_counts", min_count=2)
    go_enrichment.to_csv(os.path.join(OUTPUT_DIR, "uc2_go_enrichment.csv"), index=False)

    # Compare pathway enrichment
    pw_enrichment = compare_enrichment(disc_enrichment, conc_enrichment,
                                        term_type="pathway_counts", min_count=2)
    pw_enrichment.to_csv(os.path.join(OUTPUT_DIR, "uc2_pathway_enrichment.csv"), index=False)

    # Compare predicate enrichment
    pred_enrichment = compare_enrichment(disc_enrichment, conc_enrichment,
                                          term_type="predicate_counts", min_count=2)
    pred_enrichment.to_csv(os.path.join(OUTPUT_DIR, "uc2_predicate_enrichment.csv"), index=False)

    print(f"\n  Top enriched GO terms in discordant genes:")
    if not go_enrichment.empty:
        for _, row in go_enrichment.head(10).iterrows():
            sig = "*" if row["fisher_pval"] < 0.05 else ""
            print(f"    {row['term'][:60]:<60} fold={row['fold_enrichment']:.2f}  "
                  f"p={row['fisher_pval']:.3e} {sig}")

    print(f"\n  Top enriched predicates in discordant genes:")
    if not pred_enrichment.empty:
        for _, row in pred_enrichment.head(10).iterrows():
            sig = "*" if row["fisher_pval"] < 0.05 else ""
            print(f"    {row['term'][:60]:<60} fold={row['fold_enrichment']:.2f}  "
                  f"p={row['fisher_pval']:.3e} {sig}")

    # --- Visualizations ---
    print("\n[7/7] Generating visualizations...")

    # Each plot is wrapped in try/except so one failure doesn't block the rest

    # 1. Scatter plot: CNV-mRNA vs mRNA-Protein
    try:
        print(f"\n  Building scatter plot of CNV-mRNA vs mRNA-Protein correlations...")
        plot_discordance_scatter(
            corr_df,
            output_path=os.path.join(OUTPUT_DIR, "uc2_scatter_discordance.png")
        )
    except Exception as exc:
        print(f"  ERROR in scatter plot: {exc}")

    # 2. Network: discordant genes colored by GO process
    try:
        print(f"\n  Building network plot for discordant genes in KG...")
        all_discordant = corr_df[
            corr_df["category"].isin([
                "transcriptional_discordant",
                "post_transcriptional_discordant",
                "mixed_discordant",
            ])
        ].copy()
        if not all_discordant.empty:
            print(f"    Plotting top {min(40, len(all_discordant))} discordant genes in KG network...")
            plot_discordant_network(
                db, all_discordant, disc_enrichment,
                top_n_genes=40,
                output_path=os.path.join(OUTPUT_DIR, "uc2_network_discordant.png")
            )
    except Exception as exc:
        print(f"  ERROR in network plot: {exc}")

    # 3. Heatmap: genes x layers by category
    try:
        print(f"\n  Building heatmap of layer correlations...")
        plot_layer_heatmap(
            corr_df, gene_info, top_n=60,
            output_path=os.path.join(OUTPUT_DIR, "uc2_heatmap_layers.png")
        )
    except Exception as exc:
        print(f"  ERROR in heatmap: {exc}")

    # 4. Enrichment bar charts
    for enrichment_data, term_label, top_n_terms, fname in [
        (go_enrichment, "GO Term", 20, "uc2_enrichment_go.png"),
        (pred_enrichment, "KG Predicate", 15, "uc2_enrichment_predicates.png"),
        (pw_enrichment, "Pathway", 15, "uc2_enrichment_pathways.png"),
    ]:
        try:
            print(f"\n  Building enrichment bar chart for {term_label}...")
            plot_enrichment_bars(
                enrichment_data, term_type_label=term_label,
                top_n=top_n_terms,
                output_path=os.path.join(OUTPUT_DIR, fname)
            )
        except Exception as exc:
            print(f"  ERROR in {term_label} enrichment plot: {exc}")

    # 5. Publication multi-panel figure (Bioinformatics-style)
    try:
        print(f"\n  Building publication multi-panel figure (A|B,C)...")
        plot_publication_figure(
            corr_df, go_enrichment, pred_enrichment,
            output_path=os.path.join(OUTPUT_DIR, "uc2_figure_publication.png"),
        )
    except Exception as exc:
        print(f"  ERROR in publication figure: {exc}")

    print("\n" + "=" * 70)
    print("  UC2 analysis complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print("=" * 70)

    return corr_df, summary_df, disc_enrichment, conc_enrichment


def run_uc2_plots_only():
    """
    Skip the analysis and regenerate all plots from saved CSV files in OUTPUT_DIR.
    Requires: uc2_correlation_results.csv, uc2_go_enrichment.csv,
              uc2_pathway_enrichment.csv, uc2_predicate_enrichment.csv

    NOTE: The network plot (uc2_network_discordant.png) requires a live ArangoDB
    connection to query KG edges and cannot be regenerated in skip-analysis mode.
    All other plots are regenerated from CSVs only.
    """
    print("=" * 70)
    print("  UC2: Regenerating plots from saved results (skip-analysis mode)")
    print("=" * 70)

    corr_path = os.path.join(OUTPUT_DIR, "uc2_correlation_results.csv")
    go_path = os.path.join(OUTPUT_DIR, "uc2_go_enrichment.csv")
    pw_path = os.path.join(OUTPUT_DIR, "uc2_pathway_enrichment.csv")
    pred_path = os.path.join(OUTPUT_DIR, "uc2_predicate_enrichment.csv")

    missing = [p for p in [corr_path] if not os.path.exists(p)]
    if missing:
        print(f"  ERROR: Missing required files:\n  " + "\n  ".join(missing))
        print("  Run without --skip-analysis first to generate the data.")
        return

    print(f"  Loading: {corr_path}")
    corr_df = pd.read_csv(corr_path)

    # gene_info: entrez_id -> gene_symbol (from the saved CSV)
    gene_info = {}
    if "gene_symbol" in corr_df.columns:
        gene_info = dict(zip(corr_df["entrez_id"].astype(str),
                             corr_df["gene_symbol"].fillna("").astype(str)))

    print(f"  Genes loaded: {len(corr_df)}")
    if "category" in corr_df.columns:
        print(f"  Categories: {corr_df['category'].value_counts().to_dict()}")

    print("\n  Generating visualizations...")

    # 1. Scatter: CNV-mRNA vs mRNA-Protein
    try:
        plot_discordance_scatter(
            corr_df,
            output_path=os.path.join(OUTPUT_DIR, "uc2_scatter_discordance.png")
        )
    except Exception as exc:
        print(f"  ERROR in scatter plot: {exc}")

    # 2. Network: SKIPPED — requires live ArangoDB connection
    print("\n  NOTE: Network plot skipped (requires live ArangoDB connection).")
    print("        Existing uc2_network_discordant.png is preserved.")

    # 3. Heatmap: genes x layers by category
    try:
        plot_layer_heatmap(
            corr_df, gene_info, top_n=60,
            output_path=os.path.join(OUTPUT_DIR, "uc2_heatmap_layers.png")
        )
    except Exception as exc:
        print(f"  ERROR in heatmap: {exc}")

    # 4. Enrichment bar charts (from saved CSVs)
    enrichment_files = [
        (go_path,   "GO Term",       20, "uc2_enrichment_go.png"),
        (pred_path, "KG Predicate",  15, "uc2_enrichment_predicates.png"),
        (pw_path,   "Pathway",       15, "uc2_enrichment_pathways.png"),
    ]
    loaded_enrichment = {}
    for fpath, term_label, top_n_terms, fname in enrichment_files:
        if not os.path.exists(fpath):
            print(f"  WARNING: {fpath} not found, skipping {term_label} enrichment plot.")
            continue
        try:
            df = pd.read_csv(fpath)
            loaded_enrichment[term_label] = df
            plot_enrichment_bars(
                df, term_type_label=term_label,
                top_n=top_n_terms,
                output_path=os.path.join(OUTPUT_DIR, fname)
            )
        except Exception as exc:
            print(f"  ERROR in {term_label} enrichment plot: {exc}")

    # 5. Publication multi-panel figure (Bioinformatics-style)
    try:
        print(f"\n  Building publication multi-panel figure (A|B,C)...")
        plot_publication_figure(
            corr_df,
            loaded_enrichment.get("GO Term", pd.DataFrame()),
            loaded_enrichment.get("KG Predicate", pd.DataFrame()),
            output_path=os.path.join(OUTPUT_DIR, "uc2_figure_publication.png"),
        )
    except Exception as exc:
        print(f"  ERROR in publication figure: {exc}")

    print("\n" + "=" * 70)
    print("  UC2 plots regenerated!")
    print(f"  Output: {OUTPUT_DIR}")
    print("=" * 70)


#%% Entry point
if __name__ == "__main__":
    skip_analysis = "--skip-analysis" in sys.argv

    if "ipykernel" not in sys.modules:
        if skip_analysis:
            run_uc2_plots_only()
        else:
            run_uc2()
    else:
        # Interactive: set skip_analysis = True to regenerate plots only
        if skip_analysis:
            run_uc2_plots_only()
        else:
            corr_df, summary_df, disc_enrichment, conc_enrichment = run_uc2()
