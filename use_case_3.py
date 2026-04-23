#%%
"""
UC3 - Phenotype-anchored trans-omic subgraph extraction

Hypothesis: Starting from a clinically relevant phenotype node in the KG
(HP ontology), a semantically grounded gene set can be derived via graph
traversal that shows coherent and interpretable multi-omic patterns across
TCGA-BRCA samples — patterns that a purely statistical gene set would lack.

Architecture:
1. Seed a phenotype node (HP ontology) in the KG
2. Traverse 2-hop: phenotype -> genes (via "has phenotype"/"phenotype of"/
   "causes or contributes to condition") -> proteins (via "has gene product")
3. Extract multi-omic profiles (TPM, CNV, methylation beta, RPPA) for the
   semantically grounded gene set across all BRCA tumor samples
4. Visualize the trans-omic map of the phenotype

Visualizations:
- Semantic subgraph: phenotype -> genes -> proteins -> GO terms, colored by
  mean expression in BRCA
- Multi-omic heatmap: genes (rows) x samples (columns), annotated by layer
- Radar/spider chart: mean multi-omic profile per layer for the gene set
- Sankey diagram: phenotype -> genes -> pathways -> omic layers
"""

#%% Imports
import os
import sys
import json
import random
import warnings
from collections import defaultdict, Counter
from typing import Dict, List, Optional, Tuple, Set

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from matplotlib.colors import Normalize, TwoSlopeNorm
import seaborn as sns
import networkx as nx

# Add scripts dir to path for local imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "scripts"))

from arangodb_utils import setup_arangodb_connection
from query_utils import (
    OMIC_CONFIG,
    get_gene_expression_matrix,
    get_gene_expression_by_gene,
    get_cnv_by_gene,
    get_methylation_by_gene,
    get_protein_abundance,
    resolve_protein_position,
    list_samples_with_complete_omics,
    get_gene_vector,
    get_node_neighbors,
    search_nodes,
    get_node_by_key,
)

warnings.filterwarnings("ignore", category=RuntimeWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================
DB_NAME = "PKT_test10000"
COHORT = "TCGA-BRCA"

# Seed phenotype: Breast carcinoma (HP:0003002)
# Alternative seeds can be tried: HP:0002671 (Basal cell carcinoma),
# HP:0010788 (Testicular neoplasm), or any HP term relevant to breast.
SEED_HP_ENTITY_ID = "HP_0003002"  # _key format in nodes collection

# Traversal configuration — gene set composition
# The gene set mixes hop-1 (direct phenotype→gene links) and hop-2
# (phenotype→disease→gene links) to capture both GWAS-style associations
# AND causal driver genes reachable only through disease intermediaries.
# Breast-specific MONDO diseases are prioritized at hop 2.
MAX_GENES_HOP1 = 100     # direct phenotype→gene associations
MAX_GENES_HOP2 = 100     # genes via disease intermediaries (drivers live here)
MAX_GENES_IN_SET = 200    # total cap
MIN_GENES_IN_SET = 10     # minimum to proceed with analysis

# Random baseline for statistical comparison
N_RANDOM_GENES = 200      # size-matched random gene set (expression / CNV / methylation)
N_RANDOM_PROTEIN = 100    # protein random set — total RPPA coverage ~300 genes, so 200 is unreasonable
SEED = 42

# Sample filtering
MIN_SAMPLES = 20

# Omic layers to extract
OMIC_LAYERS = ["gene_expression", "cnv", "methylation", "protein"]

# CNV data note: values in the CNV collection are absolute copy numbers
# (diploid = 2). For visualization and statistics, we transform to
# log2(CN/2) so that 0 = diploid, positive = gain, negative = loss.

# Output directory
OUTPUT_DIR = os.path.join(os.path.dirname(__file__), "results", "uc3")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# =============================================================================
# STEP 1: KG Traversal — Phenotype-anchored gene set extraction
# =============================================================================

def find_phenotype_seed(db, hp_key: str) -> Optional[Dict]:
    """
    Find the HP phenotype seed node in the KG.
    Tries by _key first, then by entity_id, then by label search.
    """
    # Try direct key lookup
    node = get_node_by_key(db, hp_key)
    if node and node.get("class_code") == "HP":
        return node

    # Try by entity_id
    aql = """
    FOR n IN nodes
        FILTER n.class_code == "HP"
        FILTER n.entity_id == @entityId OR n._key == @entityId
        LIMIT 1
        RETURN n
    """
    results = list(db.aql.execute(aql, bind_vars={"entityId": hp_key}))
    if results:
        return results[0]

    # Try label search for "breast carcinoma"
    aql = """
    FOR n IN nodes
        FILTER n.class_code == "HP"
        FILTER CONTAINS(LOWER(n.label), "breast carcinoma")
           OR CONTAINS(LOWER(n.label), "breast neoplasm")
        LIMIT 5
        RETURN n
    """
    results = list(db.aql.execute(aql, bind_vars={}))
    if results:
        print(f"  Found {len(results)} HP nodes matching breast carcinoma:")
        for r in results:
            print(f"    {r['_key']}: {r.get('label', 'N/A')}")
        return results[0]

    return None


def traverse_phenotype_to_genes(db, seed_key: str,
                                 max_hop1: int = MAX_GENES_HOP1,
                                 max_hop2: int = MAX_GENES_HOP2) -> Dict:
    """
    Traverse KG from phenotype seed to build a semantically grounded gene set.

    Two-phase strategy to capture both association-level and causal-driver genes:

    Phase 1 (hop 1): Genes directly linked to the phenotype via
        "causes or contributes to condition" — typically GWAS-level associations.
        Capped at max_hop1 to leave room for drivers.

    Phase 2 (hop 2): Phenotype -> MONDO diseases -> genes via
        "causes or contributes to condition" / "disease has basis in dysfunction of".
        Breast-specific diseases are prioritized. This is where canonical driver
        genes (BRCA1, TP53, PTEN, PIK3CA, etc.) are found — they connect to
        HP:Breast carcinoma indirectly through disease entities like
        "hereditary breast carcinoma" or "breast-ovarian cancer familial".

    Returns dict with:
    - genes: {entrez_id: {label, kg_key, hop_distance, path_predicates, ...}}
    - diseases: {mondo_key: {label, kg_key, hop_distance, is_breast_specific}}
    - proteins: {pr_key: {label, ...}}
    - pathways: {rhsa_key: {label, ...}}
    - subgraph_edges: list of (src_key, tgt_key, predicate_label) tuples
    """
    genes = {}
    diseases = {}
    proteins = {}
    pathways = {}
    subgraph_edges = []

    # ── Phase 1: direct phenotype → gene links (hop 1) ──────────────────
    print("  Phase 1: direct phenotype -> gene links...")
    aql_hop1 = """
    FOR v, e IN 1..1 ANY CONCAT("nodes/", @seedKey) edges
        FILTER e.predicate_label IN ["causes or contributes to condition",
                                     "causally influences", "causally influenced by"]
        RETURN {vertex: v, edge: e}
    """
    hop1_results = list(db.aql.execute(aql_hop1, bind_vars={"seedKey": seed_key}))

    for r in hop1_results:
        v, e = r["vertex"], r["edge"]
        cc = v.get("class_code", "")
        src_key = e.get("_from", "").replace("nodes/", "")
        tgt_key = e.get("_to", "").replace("nodes/", "")
        pred = e.get("predicate_label", "")
        subgraph_edges.append((src_key, tgt_key, pred))

        if cc == "EntrezID":
            eid = v.get("entity_id", v.get("_key"))
            if eid and eid not in genes:
                genes[eid] = {
                    "label": v.get("label", ""),
                    "kg_key": v.get("_key"),
                    "hop_distance": 1,
                    "path_predicates": [pred],
                    "class_code": cc,
                    "source": "direct",
                }
        elif cc == "MONDO":
            diseases[v.get("_key")] = {
                "label": v.get("label", ""),
                "kg_key": v.get("_key"),
                "hop_distance": 1,
                "is_breast_specific": False,
            }

    hop1_total = len([g for g in genes.values() if g["source"] == "direct"])
    print(f"    Found {hop1_total} genes at hop 1")

    # Cap hop-1 genes
    if hop1_total > max_hop1:
        rng = random.Random(SEED)
        hop1_ids = [eid for eid, g in genes.items() if g["source"] == "direct"]
        keep = set(rng.sample(hop1_ids, max_hop1))
        genes = {eid: g for eid, g in genes.items() if eid in keep}
        print(f"    Sampled down to {max_hop1} hop-1 genes")

    # ── Phase 2: phenotype → MONDO diseases → genes (hop 2) ─────────────
    # First collect MONDO diseases linked to the phenotype
    print("  Phase 2: phenotype -> diseases -> genes (driver genes)...")
    aql_mondo = """
    FOR v, e IN 1..1 ANY CONCAT("nodes/", @seedKey) edges
        FILTER v.class_code == "MONDO"
        FILTER e.predicate_label IN ["has phenotype", "phenotype of"]
        RETURN DISTINCT {key: v._key, label: v.label, id: v._id}
    """
    mondo_results = list(db.aql.execute(aql_mondo, bind_vars={"seedKey": seed_key}))
    print(f"    Found {len(mondo_results)} linked MONDO diseases")

    # Classify diseases: breast-specific vs. other
    breast_mondo = []
    other_mondo = []
    for m in mondo_results:
        label_lower = (m.get("label") or "").lower()
        is_breast = ("breast" in label_lower or "mammary" in label_lower
                     or "ovarian cancer" in label_lower)
        diseases[m["key"]] = {
            "label": m.get("label", ""),
            "kg_key": m["key"],
            "hop_distance": 1,
            "is_breast_specific": is_breast,
        }
        subgraph_edges.append((seed_key, m["key"], "has phenotype"))
        if is_breast:
            breast_mondo.append(m)
        else:
            other_mondo.append(m)

    print(f"    Breast-specific diseases: {len(breast_mondo)}, other: {len(other_mondo)}")

    # Extract genes from diseases — breast-specific first, then others
    hop2_genes = {}
    for disease_batch, batch_label in [(breast_mondo, "breast-specific"),
                                        (other_mondo, "other")]:
        if not disease_batch:
            continue
        mondo_keys = [m["key"] for m in disease_batch]
        aql_d2g = """
        FOR mk IN @mondoKeys
            FOR v, e IN 1..1 ANY CONCAT("nodes/", mk) edges
                FILTER v.class_code == "EntrezID"
                FILTER e.predicate_label IN ["causes or contributes to condition",
                                              "disease has basis in dysfunction of"]
                RETURN DISTINCT {
                    gene_eid: v.entity_id, gene_key: v._key,
                    gene_label: v.label, disease_key: mk,
                    pred: e.predicate_label
                }
        """
        d2g_results = list(db.aql.execute(aql_d2g, bind_vars={"mondoKeys": mondo_keys}))

        for r in d2g_results:
            eid = r["gene_eid"]
            if eid and eid not in genes and eid not in hop2_genes:
                disease_info = diseases.get(r["disease_key"], {})
                hop2_genes[eid] = {
                    "label": r["gene_label"],
                    "kg_key": r["gene_key"],
                    "hop_distance": 2,
                    "path_predicates": ["has phenotype", r["pred"]],
                    "class_code": "EntrezID",
                    "source": f"hop2_{batch_label}",
                    "via_disease": disease_info.get("label", ""),
                    "via_disease_key": r["disease_key"],
                }
                subgraph_edges.append((r["disease_key"], r["gene_key"], r["pred"]))

        print(f"    {batch_label}: {len([g for g in hop2_genes.values() if g['source'] == f'hop2_{batch_label}'])} new genes")

    # Add hop-2 genes (cap at max_hop2, breast-specific prioritized)
    breast_hop2 = {k: v for k, v in hop2_genes.items()
                   if v["source"] == "hop2_breast-specific"}
    other_hop2 = {k: v for k, v in hop2_genes.items()
                  if v["source"] == "hop2_other"}

    # Add all breast-specific hop-2 genes first
    added_hop2 = 0
    for eid, info in breast_hop2.items():
        if added_hop2 >= max_hop2:
            break
        genes[eid] = info
        added_hop2 += 1

    # Fill remaining with other hop-2 genes
    for eid, info in other_hop2.items():
        if added_hop2 >= max_hop2:
            break
        genes[eid] = info
        added_hop2 += 1

    print(f"    Added {added_hop2} hop-2 genes (total gene set: {len(genes)})")

    # ── Collect proteins and pathways from 1-hop around gene nodes ───────
    gene_keys = [g["kg_key"] for g in genes.values() if g.get("kg_key")]
    if gene_keys:
        aql_ctx = """
        FOR gk IN @geneKeys
            FOR v, e IN 1..1 ANY CONCAT("nodes/", gk) edges
                FILTER e.predicate_label IN ["has gene product", "gene product of",
                                              "has participant", "participates in"]
                FILTER v.class_code IN ["PR", "R-HSA"]
                LIMIT 2000
                RETURN DISTINCT {vertex: v, gene_key: gk, pred: e.predicate_label}
        """
        ctx_results = list(db.aql.execute(aql_ctx, bind_vars={"geneKeys": gene_keys[:200]}))
        for r in ctx_results:
            v = r["vertex"]
            cc = v.get("class_code", "")
            if cc == "PR":
                proteins[v["_key"]] = {
                    "label": v.get("label", ""),
                    "kg_key": v["_key"],
                    "entity_id": v.get("entity_id", ""),
                    "hop_distance": 2,
                }
            elif cc == "R-HSA":
                pathways[v["_key"]] = {
                    "label": v.get("label", ""),
                    "kg_key": v["_key"],
                    "hop_distance": 2,
                }
            subgraph_edges.append((r["gene_key"], v["_key"], r["pred"]))

    # Summary
    n_hop1 = len([g for g in genes.values() if g.get("hop_distance") == 1])
    n_hop2 = len([g for g in genes.values() if g.get("hop_distance") == 2])
    n_breast_d = len([d for d in diseases.values() if d.get("is_breast_specific")])
    print(f"\n  Final gene set: {len(genes)} genes "
          f"({n_hop1} hop-1 + {n_hop2} hop-2)")
    print(f"  Diseases: {len(diseases)} total ({n_breast_d} breast-specific)")
    print(f"  Proteins: {len(proteins)}, Pathways: {len(pathways)}")

    return {
        "genes": genes,
        "proteins": proteins,
        "diseases": diseases,
        "pathways": pathways,
        "subgraph_edges": subgraph_edges,
    }


# =============================================================================
# STEP 2: Multi-omic data extraction for the gene set
# =============================================================================

def extract_multiomic_profiles(db, cohort: str,
                                gene_set: Dict,
                                sample_ids: List[str]) -> Dict[str, pd.DataFrame]:
    """
    Extract multi-omic data matrices for the gene set across samples.

    Returns dict of DataFrames:
    - "expression": samples x genes (TPM)
    - "cnv": samples x genes (copy number)
    - "methylation": samples x genes (mean beta across probes)
    - "protein": samples x genes (RPPA abundance)
    """
    entrez_ids = sorted(gene_set.keys())
    gene_labels = {eid: gene_set[eid].get("label", eid).split(" (")[0]
                   for eid in entrez_ids}

    print(f"  Extracting multi-omic data for {len(entrez_ids)} genes "
          f"across {len(sample_ids)} samples...")

    # --- Gene Expression (TPM) ---
    print("    Gene expression (TPM)...")
    expr_data = get_gene_expression_matrix(
        db, cohort, gene_ids=entrez_ids,
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
        print(f"      {expr_df.shape[0]} samples x {expr_df.shape[1]} genes")

    # --- CNV (absolute copy number -> log2 ratio) ---
    # Raw values are absolute copy numbers (2 = diploid normal).
    # We transform to log2(CN / 2) so that 0 = diploid, >0 = gain, <0 = loss.
    # This makes CNV comparable across genes and biologically interpretable.
    print("    CNV (copy number -> log2 ratio)...")
    cnv_records = {}
    for eid in entrez_ids:
        data = get_cnv_by_gene(db, cohort, entrez_id=eid, sample_ids=sample_ids)
        if data:
            cnv_records[eid] = {d["sample_id"]: d["copy_number"] for d in data}
    cnv_df = None
    if cnv_records:
        cnv_df = pd.DataFrame(cnv_records).apply(pd.to_numeric, errors="coerce")
        # Transform: log2(CN / 2), clamp CN >= 0.01 to avoid log(0)
        cnv_df = np.log2(cnv_df.clip(lower=0.01) / 2.0)
        print(f"      {cnv_df.shape[0]} samples x {cnv_df.shape[1]} genes "
              f"(log2 ratio, range [{cnv_df.min().min():.2f}, {cnv_df.max().max():.2f}])")

    # --- Methylation (mean beta per gene) ---
    print("    Methylation (beta)...")
    meth_records = {}
    # Need HGNC symbols for methylation lookup
    for eid in entrez_ids:
        symbol = gene_labels.get(eid, "")
        if not symbol or symbol == eid:
            continue
        data = get_methylation_by_gene(db, cohort, gene_symbol=symbol,
                                        sample_ids=sample_ids)
        if data and data.get("probes") and data.get("samples"):
            for s in data["samples"]:
                sid = s["sample_id"]
                vals = [v for v in s["values"] if v is not None]
                if vals:
                    meth_records.setdefault(eid, {})[sid] = np.mean(vals)
    meth_df = None
    if meth_records:
        meth_df = pd.DataFrame(meth_records).apply(pd.to_numeric, errors="coerce")
        print(f"      {meth_df.shape[0]} samples x {meth_df.shape[1]} genes")

    # --- Protein (RPPA) ---
    print("    Protein (RPPA abundance)...")
    prot_records = {}
    for eid in entrez_ids:
        data = get_protein_abundance(db, cohort, entrez_id=eid,
                                      sample_ids=sample_ids)
        if data:
            prot_records[eid] = {d["sample_id"]: d["abundance"] for d in data}
    prot_df = None
    if prot_records:
        prot_df = pd.DataFrame(prot_records).apply(pd.to_numeric, errors="coerce")
        print(f"      {prot_df.shape[0]} samples x {prot_df.shape[1]} genes")

    return {
        "expression": expr_df,
        "cnv": cnv_df,
        "methylation": meth_df,
        "protein": prot_df,
    }


# =============================================================================
# STEP 3: Visualization — Semantic Subgraph
# =============================================================================

def plot_semantic_subgraph(traversal_result: Dict,
                           seed_node: Dict,
                           mean_expression: Optional[pd.Series] = None,
                           title: str = "UC3: Phenotype-Anchored Semantic Subgraph",
                           output_path: str = None,
                           max_nodes: int = 120):
    """
    Network visualization of the semantic subgraph:
    phenotype -> genes -> proteins -> GO terms
    Nodes colored by mean expression in BRCA (if available).
    """
    G = nx.DiGraph()

    # Add seed phenotype node
    seed_key = seed_node["_key"]
    seed_label = seed_node.get("label", seed_key)
    G.add_node(seed_key, label=seed_label, node_type="phenotype",
               mean_expr=0)

    genes = traversal_result["genes"]
    proteins = traversal_result["proteins"]
    diseases = traversal_result["diseases"]
    pathways = traversal_result["pathways"]

    # Add gene nodes
    for eid, info in list(genes.items())[:max_nodes]:
        short_label = info["label"].split(" (")[0] if info["label"] else eid
        mean_val = 0
        if mean_expression is not None and eid in mean_expression.index:
            mean_val = mean_expression[eid]
        G.add_node(info["kg_key"], label=short_label, node_type="gene",
                   mean_expr=mean_val, entrez_id=eid)

    # Add disease nodes (subset)
    for dk, dinfo in list(diseases.items())[:15]:
        short_label = dinfo["label"][:30] if dinfo["label"] else dk
        G.add_node(dk, label=short_label, node_type="disease", mean_expr=0)

    # Add pathway nodes (subset)
    for pk, pinfo in list(pathways.items())[:15]:
        short_label = pinfo["label"][:30] if pinfo["label"] else pk
        G.add_node(pk, label=short_label, node_type="pathway", mean_expr=0)

    # Add protein nodes (subset)
    for pr_key, pinfo in list(proteins.items())[:20]:
        short_label = pinfo["label"][:20] if pinfo["label"] else pr_key
        G.add_node(pr_key, label=short_label, node_type="protein", mean_expr=0)

    # Add edges (only those connecting nodes in G)
    node_keys_in_G = set(G.nodes())
    for src, tgt, pred in traversal_result["subgraph_edges"]:
        if src in node_keys_in_G and tgt in node_keys_in_G:
            G.add_edge(src, tgt, predicate=pred)

    if len(G.edges()) == 0:
        # Fallback: add edges from seed to diseases/genes found at hop 1
        for dk in diseases:
            if dk in node_keys_in_G:
                G.add_edge(seed_key, dk, predicate="has phenotype")
        for eid, info in genes.items():
            gk = info["kg_key"]
            if gk in node_keys_in_G and info["hop_distance"] <= 2:
                G.add_edge(seed_key, gk, predicate="linked")

    if len(G.nodes()) == 0:
        print("  No nodes to plot for semantic subgraph, skipping.")
        return

    fig, ax = plt.subplots(figsize=(18, 14))

    # Layout
    pos = nx.spring_layout(G, seed=42, k=2.5 / np.sqrt(max(len(G.nodes()), 1)),
                           iterations=80)

    # Node colors by type
    type_colors = {
        "phenotype": "#E74C3C",
        "gene": "#3498DB",
        "protein": "#2ECC71",
        "disease": "#F39C12",
        "pathway": "#9B59B6",
    }
    node_types = nx.get_node_attributes(G, "node_type")
    node_colors = [type_colors.get(node_types.get(n, "gene"), "#CCCCCC")
                   for n in G.nodes()]

    # Node sizes: phenotype large, genes medium, others small
    node_sizes = []
    for n in G.nodes():
        nt = node_types.get(n, "gene")
        if nt == "phenotype":
            node_sizes.append(800)
        elif nt == "gene":
            node_sizes.append(200)
        else:
            node_sizes.append(120)

    # If we have expression data, modulate gene node intensity
    if mean_expression is not None:
        mean_exprs = nx.get_node_attributes(G, "mean_expr")
        valid_exprs = [v for v in mean_exprs.values() if v != 0]
        if valid_exprs:
            vmin, vmax = np.percentile(valid_exprs, [5, 95])
            norm = Normalize(vmin=max(vmin, 0.01), vmax=max(vmax, 1))
            for i, n in enumerate(G.nodes()):
                if node_types.get(n) == "gene" and mean_exprs.get(n, 0) > 0:
                    intensity = norm(mean_exprs[n])
                    node_colors[i] = plt.cm.Blues(0.3 + 0.7 * intensity)

    # Draw
    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#CCCCCC",
                           alpha=0.4, arrows=True, arrowsize=8,
                           connectionstyle="arc3,rad=0.1")
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=node_sizes, alpha=0.85,
                           edgecolors="white", linewidths=0.5)

    # Labels for high-degree or important nodes
    degrees = dict(G.degree())
    degree_thresh = np.percentile(list(degrees.values()), 60) if degrees else 0
    labels_dict = {}
    node_labels = nx.get_node_attributes(G, "label")
    for n in G.nodes():
        nt = node_types.get(n, "")
        if nt == "phenotype" or degrees.get(n, 0) >= degree_thresh:
            labels_dict[n] = node_labels.get(n, n)[:20]
    nx.draw_networkx_labels(G, pos, labels_dict, font_size=7,
                            font_weight="bold", ax=ax)

    # Legend
    legend_handles = [mpatches.Patch(color=c, label=t.capitalize())
                      for t, c in type_colors.items()]
    ax.legend(handles=legend_handles, loc="upper left", fontsize=9,
              framealpha=0.8)

    ax.set_title(title, fontsize=14, fontweight="bold")
    ax.axis("off")
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


# =============================================================================
# STEP 4: Visualization — Multi-omic Heatmap
# =============================================================================

def plot_multiomic_heatmap(omic_dfs: Dict[str, pd.DataFrame],
                            gene_set: Dict,
                            title: str = "UC3: Multi-Omic Heatmap (Phenotype Gene Set)",
                            output_path: str = None,
                            max_genes: int = 60,
                            max_samples: int = 80):
    """
    Clustered heatmap: genes (rows) x samples (columns), with side bars
    for each omic layer (TPM, CNV, methylation beta, RPPA).
    """
    # Collect available layers
    available = {}
    for layer_name, df in omic_dfs.items():
        if df is not None and not df.empty:
            available[layer_name] = df

    if not available:
        print("  No omic data available for heatmap, skipping.")
        return

    # Find common genes and samples
    all_genes_per_layer = [set(df.columns) for df in available.values()]
    common_genes = sorted(set.intersection(*all_genes_per_layer)) if all_genes_per_layer else []

    if not common_genes:
        # Use union and fill NaN
        common_genes = sorted(set.union(*all_genes_per_layer))

    # Limit genes
    if len(common_genes) > max_genes:
        # Prefer genes with data across more layers
        gene_coverage = {}
        for g in common_genes:
            gene_coverage[g] = sum(1 for df in available.values()
                                   if g in df.columns and df[g].notna().any())
        common_genes = sorted(gene_coverage.keys(),
                              key=lambda g: -gene_coverage[g])[:max_genes]

    # Gene labels
    gene_labels = {}
    for g in common_genes:
        info = gene_set.get(g, {})
        lbl = info.get("label", g)
        gene_labels[g] = lbl.split(" (")[0] if " (" in str(lbl) else str(lbl)

    n_layers = len(available)
    fig, axes = plt.subplots(1, n_layers, figsize=(5 * n_layers, max(8, len(common_genes) * 0.18)),
                              sharey=True)
    if n_layers == 1:
        axes = [axes]

    layer_cmaps = {
        "expression": "Reds",
        "cnv": "RdBu_r",
        "methylation": "YlGn",
        "protein": "PuBu",
    }
    layer_titles = {
        "expression": "Gene Expression (TPM)",
        "cnv": "Copy Number Variation",
        "methylation": "Methylation (Beta)",
        "protein": "Protein Abundance (RPPA)",
    }

    for ax_idx, (layer_name, df) in enumerate(available.items()):
        ax = axes[ax_idx]

        # Filter to common genes present in this layer
        genes_in_layer = [g for g in common_genes if g in df.columns]
        if not genes_in_layer:
            ax.set_title(layer_titles.get(layer_name, layer_name))
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes)
            continue

        sub_df = df[genes_in_layer].dropna(how="all")
        if sub_df.empty:
            continue

        # Limit samples
        if sub_df.shape[0] > max_samples:
            sub_df = sub_df.sample(n=max_samples, random_state=42)

        # Z-score normalization for visualization
        sub_z = sub_df.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x)
        sub_z = sub_z.clip(-3, 3)

        # Rename columns to gene symbols
        rename_map = {g: gene_labels.get(g, g) for g in genes_in_layer}
        sub_z = sub_z.rename(columns=rename_map)

        cmap = layer_cmaps.get(layer_name, "viridis")
        if layer_name == "cnv":
            center = 0
        else:
            center = None

        sns.heatmap(sub_z.T, ax=ax, cmap=cmap, center=center,
                    xticklabels=False,
                    yticklabels=(ax_idx == 0),
                    cbar_kws={"shrink": 0.5, "label": "z-score"})

        ax.set_title(layer_titles.get(layer_name, layer_name),
                     fontsize=11, fontweight="bold")
        ax.set_xlabel(f"Samples (n={sub_df.shape[0]})", fontsize=9)
        if ax_idx == 0:
            ax.set_ylabel("Genes", fontsize=10)

    fig.suptitle(title, fontsize=14, fontweight="bold", y=1.02)
    plt.tight_layout()

    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


# =============================================================================
# STEP 5: Visualization — Radar / Spider Chart  [DEPRECATED — kept for reference]
# The radar chart is commented out because its per-axis normalization is
# relative (each axis divided by max(phenotype, random)), making the visual
# shape sensitive to the random baseline, which itself varies across runs
# due to non-deterministic ArangoDB query ordering. Replaced by
# plot_omic_barplot_comparison() and plot_bootstrap_profile() below.
# =============================================================================

# def plot_radar_chart(omic_dfs: Dict[str, pd.DataFrame],
#                       gene_set: Dict,
#                       random_omic_dfs: Optional[Dict[str, pd.DataFrame]] = None,
#                       random_gene_set: Optional[Dict] = None,
#                       title: str = "UC3: Multi-Omic Radar Profile",
#                       output_path: str = None):
#     """
#     Radar chart comparing the phenotype gene set vs a real random gene set
#     across multiple omic dimensions.
#
#     Axes:
#     - Per-layer mean signal (log2(TPM+1), |CNV log2 ratio|, beta, |RPPA|)
#     - Per-layer cross-sample variance (biological variability)
#     - Omic coverage fraction
#     """
#     metric_labels = []
#     phenotype_values = []
#     baseline_values = []
#
#     # Signal intensity metrics
#     layer_metrics = {
#         "expression": ("Mean\nlog2(TPM+1)", lambda df: np.log2(df.mean(axis=0).mean() + 1)),
#         "cnv": ("Mean\n|CNV log2R|", lambda df: df.abs().mean(axis=0).mean()),
#         "methylation": ("Mean\nBeta", lambda df: df.mean(axis=0).mean()),
#         "protein": ("Mean\n|RPPA|", lambda df: df.abs().mean(axis=0).mean()),
#     }
#
#     for layer_name, (label, fn) in layer_metrics.items():
#         pheno_df = omic_dfs.get(layer_name)
#         rand_df = random_omic_dfs.get(layer_name) if random_omic_dfs else None
#         if pheno_df is not None and not pheno_df.empty:
#             metric_labels.append(label)
#             phenotype_values.append(fn(pheno_df))
#             baseline_values.append(fn(rand_df) if rand_df is not None and not rand_df.empty else 0)
#
#     # Cross-sample variance per layer (captures biological heterogeneity)
#     for layer_name in ["expression", "cnv"]:
#         pheno_df = omic_dfs.get(layer_name)
#         rand_df = random_omic_dfs.get(layer_name) if random_omic_dfs else None
#         if pheno_df is not None and not pheno_df.empty:
#             metric_labels.append(f"{layer_name.capitalize()}\nVariance")
#             phenotype_values.append(pheno_df.var(axis=0).mean())
#             baseline_values.append(
#                 rand_df.var(axis=0).mean() if rand_df is not None and not rand_df.empty else 0)
#
#     # Coverage metrics
#     n_pheno = max(len(gene_set), 1)
#     n_random = max(len(random_gene_set), 1) if random_gene_set else 1
#     for layer_name in ["expression", "protein"]:
#         pheno_df = omic_dfs.get(layer_name)
#         rand_df = random_omic_dfs.get(layer_name) if random_omic_dfs else None
#         if pheno_df is not None and not pheno_df.empty:
#             metric_labels.append(f"{layer_name.capitalize()}\nCoverage")
#             phenotype_values.append(pheno_df.shape[1] / n_pheno)
#             baseline_values.append(
#                 rand_df.shape[1] / n_random if rand_df is not None and not rand_df.empty else 0)
#
#     if len(metric_labels) < 3:
#         print("  Not enough metrics for radar chart, skipping.")
#         return
#
#     # Normalize to [0, 1] range
#     max_vals = [max(abs(p), abs(b), 0.001) for p, b in
#                 zip(phenotype_values, baseline_values)]
#     pheno_norm = [p / m for p, m in zip(phenotype_values, max_vals)]
#     base_norm = [b / m for b, m in zip(baseline_values, max_vals)]
#
#     N = len(metric_labels)
#     angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
#     angles += angles[:1]
#     pheno_norm += pheno_norm[:1]
#     base_norm += base_norm[:1]
#
#     fig, ax = plt.subplots(figsize=(10, 10), subplot_kw=dict(polar=True))
#
#     ax.fill(angles, pheno_norm, alpha=0.25, color="#E74C3C",
#             label="Phenotype Gene Set (semantic)")
#     ax.plot(angles, pheno_norm, "o-", linewidth=2, color="#E74C3C",
#             markersize=6)
#
#     ax.fill(angles, base_norm, alpha=0.15, color="#3498DB",
#             label="Random Gene Set (size-matched)")
#     ax.plot(angles, base_norm, "o--", linewidth=1.5, color="#3498DB",
#             markersize=5)
#
#     ax.set_xticks(angles[:-1])
#     ax.set_xticklabels(metric_labels, fontsize=9)
#     ax.set_ylim(0, 1.3)
#     ax.set_title(title, fontsize=14, fontweight="bold", pad=30)
#     ax.legend(loc="upper right", bbox_to_anchor=(1.35, 1.1), fontsize=10)
#
#     plt.tight_layout()
#     if output_path:
#         fig.savefig(output_path, dpi=200, bbox_inches="tight")
#         print(f"  Saved: {output_path}")
#     plt.show()
#     return fig


# =============================================================================
# STEP 5a: Visualization — Omic Barplot Comparison (replaces radar)
# =============================================================================

def plot_omic_barplot_comparison(omic_dfs: Dict[str, pd.DataFrame],
                                  random_omic_dfs: Dict[str, pd.DataFrame],
                                  gene_set: Dict,
                                  random_gene_set: Dict,
                                  title: str = "UC3: Multi-Omic Signal — Phenotype vs Random",
                                  output_path: str = None) -> Optional[object]:
    """
    Grouped barplot comparing phenotype gene set vs random baseline on four
    metrics (one per omic layer), with error bars (±1 SEM of the group mean).

    Each bar shows the mean per-gene signal across the gene set (not
    normalized) so axes have interpretable biological units. Error bars
    are ±1 SEM = SD/√n of the per-gene values → they represent the
    precision of the group mean, not the cross-gene spread (which would
    be ±1 SD and would dwarf any real phenotype-vs-random difference
    because cross-gene variability is intrinsically much larger than
    the mean shift between two gene sets). A Mann-Whitney U two-sided
    p-value on the per-gene distributions is annotated above each pair.

    Metric per layer (one value per gene = cross-sample mean):
      expression  — mean log2(TPM+1)
      cnv         — mean |log2 CNV ratio|
      methylation — mean beta
      protein     — mean |RPPA|
    """
    layer_cfg = {
        "expression":  ("Expression\nlog2(TPM+1)",  "#C0392B",
                        lambda df: np.log2(df.mean(axis=0) + 1)),
        "cnv":         ("CNV\n|log2 ratio|",         "#2980B9",
                        lambda df: df.abs().mean(axis=0)),
        "methylation": ("Methylation\nMean β",        "#27AE60",
                        lambda df: df.mean(axis=0)),
        "protein":     ("Protein\n|RPPA|",            "#8E44AD",
                        lambda df: df.abs().mean(axis=0)),
    }

    rows = []
    for ln, (xlabel, color, fn) in layer_cfg.items():
        p_df = omic_dfs.get(ln)
        r_df = random_omic_dfs.get(ln)
        if p_df is None or p_df.empty:
            continue
        p_vals = fn(p_df).dropna().values
        r_vals = fn(r_df).dropna().values if (r_df is not None and not r_df.empty) else np.array([])
        if len(p_vals) == 0:
            continue
        u_stat, p_val = (stats.mannwhitneyu(p_vals, r_vals, alternative="two-sided")
                         if len(r_vals) >= 5 else (np.nan, np.nan))

        p_mean = p_vals.mean()
        r_mean = r_vals.mean() if len(r_vals) else np.nan
        # SEM = SD / sqrt(n) — precision of the group mean, not spread
        p_sem = p_vals.std(ddof=1) / np.sqrt(len(p_vals)) if len(p_vals) > 1 else 0
        r_sem = (r_vals.std(ddof=1) / np.sqrt(len(r_vals))
                 if len(r_vals) > 1 else 0)

        rows.append({
            "layer": ln, "label": xlabel, "color": color,
            "pheno_mean": p_mean, "pheno_sem": p_sem,
            "rand_mean":  r_mean, "rand_sem":  r_sem,
            "pheno_vals": p_vals, "rand_vals": r_vals,
            "p_value": p_val, "n_pheno": len(p_vals), "n_rand": len(r_vals),
        })

    if not rows:
        print("  No data for barplot comparison, skipping.")
        return None

    n = len(rows)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 4.8), sharey=False)
    if n == 1:
        axes = [axes]

    x = np.array([0, 1])
    bar_w = 0.55

    for ax, row in zip(axes, rows):
        color = row["color"]
        pheno_col = color
        rand_col  = "#95A5A6"

        means = [row["pheno_mean"], row["rand_mean"]]
        sems  = [row["pheno_sem"],  row["rand_sem"]]

        # Replace NaN mean (no random data) with 0 so the bar doesn't vanish
        plot_means = [m if not np.isnan(m) else 0 for m in means]
        plot_sems  = [s if not np.isnan(s) else 0 for s in sems]

        ax.bar(x, plot_means, width=bar_w, yerr=plot_sems, capsize=4,
               color=[pheno_col, rand_col],
               error_kw={"linewidth": 1.0, "ecolor": "#555"},
               edgecolor="white", linewidth=0.5, alpha=0.88)

        # Overlay individual per-gene points (strip/jitter) for transparency
        # on the true cross-gene distribution
        for i, vals in enumerate([row["pheno_vals"], row["rand_vals"]]):
            if len(vals) == 0:
                continue
            jitter = np.random.default_rng(SEED + i).normal(0, 0.055, len(vals))
            ax.scatter(np.full(len(vals), x[i]) + jitter, vals,
                       s=5, alpha=0.25, color="#333", zorder=3, linewidths=0)

        ax.set_xticks(x)
        ax.set_xticklabels(["Phenotype\n(semantic)", "Random\n(size-matched)"],
                           fontsize=8)
        ax.set_title(row["label"], fontsize=9, fontweight="bold", pad=4)
        ax.set_ylabel("Mean signal ± SEM", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=7)

        # Significance annotation
        p = row["p_value"]
        if not np.isnan(p):
            sig = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
            # Place annotation above the tallest point, with a bracket
            all_top = np.concatenate([row["pheno_vals"], row["rand_vals"]]) \
                        if len(row["rand_vals"]) else row["pheno_vals"]
            y_top = float(np.nanmax(all_top)) if len(all_top) else max(plot_means)
            y_top = y_top * 1.06
            ax.plot([0, 1], [y_top, y_top], color="#555", linewidth=0.8)
            ax.annotate(f"p={p:.2e}  {sig}",
                        xy=(0.5, y_top * 1.02),
                        ha="center", va="bottom", fontsize=7.5,
                        color="#333333")
        ax.text(0.5, -0.22,
                f"n={row['n_pheno']} / {row['n_rand']} genes",
                transform=ax.transAxes, ha="center", fontsize=6.5, color="grey")

    fig.suptitle(title, fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


# =============================================================================
# STEP 5b: Visualization — Bootstrap Profile (replaces radar)
# =============================================================================

# Number of bootstrap resamples for the null distribution
N_BOOTSTRAP = 50


def plot_bootstrap_profile(omic_dfs: Dict[str, pd.DataFrame],
                            valid_genes: List[str],
                            gene_set: Dict,
                            random_omic_dfs: Optional[Dict[str, pd.DataFrame]] = None,
                            title: str = "UC3: Multi-Omic Signal vs Bootstrap Null",
                            output_path: str = None,
                            n_bootstrap: int = N_BOOTSTRAP,
                            boot_seed: int = SEED) -> Optional[object]:
    """
    For each omic layer, compare the *per-gene* signal distribution of the
    phenotype gene set against a bootstrap null distribution built by pooling
    the per-gene values from N_BOOTSTRAP random gene sets of the same size
    (sampled from valid_genes, excluding phenotype genes).

    Rationale — why per-gene values and not mean-of-means:
    The mean-of-a-200-gene-set has very small variance (law of large numbers):
    50 resamples give a ~constant null with no visible spread, and the
    phenotype value looks like a single point next to another single point.
    Instead, we show the *per-gene* distributions: red violin = 200 per-gene
    phenotype values; grey violin = ~200×50 = ~10k per-gene values pooled
    from all 50 random resamples. The null then reflects the true
    cross-gene variability, and the Mann-Whitney test compares the two
    distributions directly.

    The plot is deterministic (boot_seed fixes the RNG).

    Metric per gene (one scalar per gene, cross-sample aggregate):
      expression  — log2(mean TPM + 1)
      cnv         — mean |log2 CNV ratio|
      methylation — mean beta
      protein     — mean |RPPA|
    """
    layer_cfg = {
        "expression":  ("Expression\nlog2(TPM+1)",  "#C0392B",
                        lambda df, cols: np.log2(df[cols].mean(axis=0) + 1).dropna().values
                        if len(cols) else np.array([])),
        "cnv":         ("CNV\n|log2 ratio|",         "#2980B9",
                        lambda df, cols: df[cols].abs().mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
        "methylation": ("Methylation\nMean β",        "#27AE60",
                        lambda df, cols: df[cols].mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
        "protein":     ("Protein\n|RPPA|",            "#8E44AD",
                        lambda df, cols: df[cols].abs().mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
    }

    pheno_ids = set(gene_set.keys())
    pool = [g for g in valid_genes if g and g not in pheno_ids]
    n_sample = len(gene_set)

    active_layers = [(ln, cfg) for ln, cfg in layer_cfg.items()
                     if omic_dfs.get(ln) is not None and not omic_dfs[ln].empty]
    if not active_layers:
        print("  No data for bootstrap profile, skipping.")
        return None

    # Phenotype per-gene values
    pheno_vals = {}
    for ln, (_, _, fn) in active_layers:
        df = omic_dfs[ln]
        cols = [g for g in gene_set if g in df.columns]
        pheno_vals[ln] = fn(df, cols)

    # Null — preferred: from random_omic_dfs columns (correct size-matched
    # random set with data for each layer). Legacy fallback: resample from
    # valid_genes (only works if omic_dfs contains data beyond phenotype set,
    # which is NOT the case here — leaves null empty).
    rng = np.random.default_rng(boot_seed)
    null_pools = {ln: np.array([]) for ln, _ in active_layers}
    null_set_means = {ln: [] for ln, _ in active_layers}

    if random_omic_dfs is not None:
        for ln, (_, _, fn) in active_layers:
            rand_df = random_omic_dfs.get(ln)
            if rand_df is None or rand_df.empty:
                continue
            rand_cols = list(rand_df.columns)
            vals = fn(rand_df, rand_cols)
            null_pools[ln] = vals
            # For z-score: bootstrap-resample set means from the null pool
            if len(vals) > 1:
                n_target = min(n_sample, len(vals))
                means = [float(rng.choice(vals, size=n_target, replace=True).mean())
                         for _ in range(n_bootstrap)]
                null_set_means[ln] = means
    else:
        # Legacy (broken) path — kept for backward compatibility
        _per_resample = {ln: [] for ln, _ in active_layers}
        for _ in range(n_bootstrap):
            samp = rng.choice(pool, size=min(n_sample, len(pool)),
                              replace=False).tolist()
            for ln, (_, _, fn) in active_layers:
                df = omic_dfs[ln]
                cols = [g for g in samp if g in df.columns]
                if cols:
                    vals = fn(df, cols)
                    _per_resample[ln].append(vals)
                    if len(vals):
                        null_set_means[ln].append(float(vals.mean()))
        null_pools = {ln: (np.concatenate(v) if v else np.array([]))
                      for ln, v in _per_resample.items()}

    n = len(active_layers)
    fig, axes = plt.subplots(1, n, figsize=(3.2 * n, 5.0), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, (ln, (label, color, _)) in zip(axes, active_layers):
        pheno = pheno_vals[ln]
        null  = null_pools[ln]
        set_means = np.array(null_set_means[ln])

        # Violin — null (grey, left) and phenotype (colored, right)
        positions, datasets, colors, edges = [], [], [], []
        if len(null) > 1:
            positions.append(0); datasets.append(null)
            colors.append("#BDC3C7"); edges.append("#7F8C8D")
        if len(pheno) > 1:
            positions.append(1); datasets.append(pheno)
            colors.append(color);   edges.append(color)

        if datasets:
            parts = ax.violinplot(datasets, positions=positions, widths=0.75,
                                  showmedians=True, showextrema=False)
            for pc, c, e in zip(parts["bodies"], colors, edges):
                pc.set_facecolor(c); pc.set_alpha(0.55)
                pc.set_edgecolor(e); pc.set_linewidth(0.8)
            parts["cmedians"].set_color("#333"); parts["cmedians"].set_linewidth(1.0)

        # Overlay boxplot stats (Q1/Q3)
        for pos, data in zip(positions, datasets):
            q1, q3 = np.percentile(data, [25, 75])
            ax.plot([pos - 0.12, pos + 0.12], [q1, q1], color="#333", lw=0.8)
            ax.plot([pos - 0.12, pos + 0.12], [q3, q3], color="#333", lw=0.8)

        # Mann-Whitney on per-gene distributions (phenotype vs pooled null)
        p_val, sig, mean_diff = np.nan, "", np.nan
        if len(pheno) >= 5 and len(null) >= 5:
            u_stat, p_val = stats.mannwhitneyu(pheno, null, alternative="two-sided")
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else \
                  "*" if p_val < 0.05 else "ns"
            # Z against resample-set-means distribution (null of group means)
            if len(set_means) > 1 and set_means.std() > 0:
                z = (float(pheno.mean()) - set_means.mean()) / set_means.std()
            else:
                z = np.nan
        else:
            z = np.nan

        ax.set_xticks([0, 1])
        ax.set_xticklabels([f"Null\n(pooled,\n{n_bootstrap} sets)",
                            "Phenotype\n(semantic)"],
                           fontsize=7.2)
        ax.set_title(label, fontsize=9, fontweight="bold", pad=6)
        ax.set_ylabel("Per-gene signal", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(axis="y", labelsize=7)

        # Significance bracket
        if not np.isnan(p_val):
            all_vals = np.concatenate([null, pheno]) if len(null) else pheno
            y_top = float(np.nanpercentile(all_vals, 98)) * 1.05
            ax.plot([0, 1], [y_top, y_top], color="#555", lw=0.8)
            z_str = f", z={z:.2f}" if not np.isnan(z) else ""
            ax.annotate(f"p={p_val:.2e}  {sig}{z_str}",
                        xy=(0.5, y_top * 1.02),
                        ha="center", va="bottom", fontsize=7, color="#333")

        ax.text(0.5, -0.20,
                f"n_pheno={len(pheno)}, n_null={len(null)} "
                f"({n_bootstrap} resamples)",
                transform=ax.transAxes, ha="center", fontsize=6.5, color="grey")

    fig.suptitle(title, fontsize=11, fontweight="bold", y=1.02)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=200, bbox_inches="tight")
        print(f"  Saved: {output_path}")
    plt.show()
    return fig


# =============================================================================
# STEP 6: Visualization — Sankey Diagram
# =============================================================================

def plot_sankey_diagram(traversal_result: Dict,
                        omic_dfs: Dict[str, pd.DataFrame],
                        seed_node: Dict,
                        gene_set: Dict,
                        title: str = "UC3: Phenotype -> Genes -> Pathways -> Omic Layers",
                        output_path: str = None):
    """
    Sankey diagram using Plotly (standard approach for alluvial/Sankey plots).

    Flow: Phenotype -> Diseases/Direct -> Top Genes -> Omic Layers
    Link width proportional to mean gene expression.
    Exported to PNG via kaleido at 200 DPI.
    """
    try:
        import plotly.graph_objects as go
        import plotly.io as pio
    except ImportError:
        print("  plotly not available — install with: pip install plotly kaleido")
        return

    seed_label = seed_node.get("label", "Phenotype")[:35]
    diseases = traversal_result.get("diseases", {})

    # ── Top genes by mean expression ──────────────────────────────────────
    expr_df = omic_dfs.get("expression")
    if expr_df is not None and not expr_df.empty:
        mean_expr_series = expr_df.mean(axis=0).sort_values(ascending=False)
        top_genes = mean_expr_series.head(20).index.tolist()
    else:
        top_genes = list(gene_set.keys())[:20]
        mean_expr_series = pd.Series(dtype=float)

    # ── Disease list: breast-specific first, then fill to 8 ──────────────
    breast_diseases = [(k, v) for k, v in diseases.items()
                       if v.get("is_breast_specific")]
    other_diseases  = [(k, v) for k, v in diseases.items()
                       if not v.get("is_breast_specific")]
    disease_list = breast_diseases[:5] + other_diseases[:max(0, 7 - len(breast_diseases[:5]))]
    # "Direct association" node for hop-1 genes
    disease_list.insert(0, ("direct_assoc", {
        "label": "Direct association", "is_breast_specific": False,
    }))
    disease_list = disease_list[:9]

    # ── Active omic layers ────────────────────────────────────────────────
    layer_meta = {
        "expression": {"label": "Gene Expression (TPM)",    "color": "rgba(231,76,60,0.85)"},
        "cnv":        {"label": "Copy Number Variation",    "color": "rgba(52,152,219,0.85)"},
        "methylation":{"label": "DNA Methylation (Beta)",   "color": "rgba(46,204,113,0.85)"},
        "protein":    {"label": "Protein Abundance (RPPA)", "color": "rgba(155,89,182,0.85)"},
    }
    active_layers = [l for l in layer_meta
                     if omic_dfs.get(l) is not None and not omic_dfs[l].empty]

    # ── Build node index ──────────────────────────────────────────────────
    # Node order: [phenotype] + [diseases] + [genes] + [omic_layers]
    nodes = []   # list of dicts: label, color, x, y

    # Node 0: phenotype
    nodes.append({"label": seed_label, "color": "rgba(192,57,43,0.9)",
                  "x": 0.01, "y": 0.5})

    # Disease nodes
    d_idx = {}
    for i, (dk, dinfo) in enumerate(disease_list):
        dlabel = (dinfo.get("label") or dk)[:30]
        is_breast = dinfo.get("is_breast_specific", False)
        color = "rgba(211,84,0,0.85)" if is_breast else "rgba(243,156,18,0.75)"
        y_pos = 0.05 + i * (0.9 / max(len(disease_list) - 1, 1))
        nodes.append({"label": dlabel, "color": color,
                      "x": 0.28, "y": round(y_pos, 3)})
        d_idx[dk] = len(nodes) - 1

    # Gene nodes
    g_idx = {}
    for i, eid in enumerate(top_genes):
        info = gene_set.get(eid, {})
        glabel = info.get("label", eid)
        glabel = glabel.split(" (")[0] if " (" in str(glabel) else str(glabel)
        expr_val = mean_expr_series.get(eid, 0) if len(mean_expr_series) else 0
        intensity = min(0.95, max(0.35, np.log2(float(expr_val) + 1) / 14))
        r = int(30 + (1 - intensity) * 80)
        g = int(100 + intensity * 100)
        b = int(180 + intensity * 55)
        color = f"rgba({r},{g},{b},0.82)"
        y_pos = 0.02 + i * (0.96 / max(len(top_genes) - 1, 1))
        nodes.append({"label": glabel[:18], "color": color,
                      "x": 0.55, "y": round(y_pos, 3)})
        g_idx[eid] = len(nodes) - 1

    # Omic layer nodes
    l_idx = {}
    for i, ln in enumerate(active_layers):
        lm = layer_meta[ln]
        df = omic_dfs[ln]
        covered = len([g for g in top_genes if g in df.columns])
        label = f"{lm['label']}<br>({covered}/{len(top_genes)} genes)"
        y_pos = 0.1 + i * (0.8 / max(len(active_layers) - 1, 1))
        nodes.append({"label": label, "color": lm["color"],
                      "x": 0.92, "y": round(y_pos, 3)})
        l_idx[ln] = len(nodes) - 1

    # ── Build links ───────────────────────────────────────────────────────
    sources, targets, values, link_colors = [], [], [], []

    # Phenotype -> diseases (equal weight)
    for dk, _ in disease_list:
        if dk in d_idx:
            sources.append(0)
            targets.append(d_idx[dk])
            values.append(10)
            link_colors.append("rgba(192,57,43,0.25)")

    # Diseases -> genes
    for eid in top_genes:
        if eid not in g_idx:
            continue
        info = gene_set.get(eid, {})
        via_key = info.get("via_disease_key")
        src_node = d_idx.get(via_key) if via_key and via_key in d_idx \
                   else d_idx.get("direct_assoc")
        if src_node is None:
            src_node = 1  # fallback to first disease node
        expr_val = float(mean_expr_series.get(eid, 1)) if len(mean_expr_series) else 1
        weight = max(1, min(20, expr_val * 0.8))
        is_hop2 = (info.get("via_disease_key") is not None)
        lc = "rgba(211,84,0,0.30)" if is_hop2 else "rgba(52,152,219,0.30)"
        sources.append(src_node)
        targets.append(g_idx[eid])
        values.append(weight)
        link_colors.append(lc)

    # Genes -> omic layers
    omic_link_colors = {
        "expression": "rgba(231,76,60,0.22)",
        "cnv":        "rgba(52,152,219,0.22)",
        "methylation":"rgba(46,204,113,0.22)",
        "protein":    "rgba(155,89,182,0.22)",
    }
    for ln in active_layers:
        df = omic_dfs[ln]
        if ln not in l_idx:
            continue
        for eid in top_genes:
            if eid not in g_idx:
                continue
            if eid in df.columns:
                mean_val = float(df[eid].mean()) if not df[eid].isna().all() else 0
                weight = max(1, min(15, abs(mean_val) * 0.7 + 1))
                sources.append(g_idx[eid])
                targets.append(l_idx[ln])
                values.append(weight)
                link_colors.append(omic_link_colors.get(ln, "rgba(150,150,150,0.2)"))

    # ── Assemble plotly Sankey ─────────────────────────────────────────────
    node_labels  = [n["label"]  for n in nodes]
    node_colors  = [n["color"]  for n in nodes]
    node_x       = [n["x"]      for n in nodes]
    node_y       = [n["y"]      for n in nodes]

    fig = go.Figure(go.Sankey(
        arrangement="snap",
        node=dict(
            pad=18,
            thickness=22,
            line=dict(color="white", width=0.8),
            label=node_labels,
            color=node_colors,
            x=node_x,
            y=node_y,
            hovertemplate="%{label}<extra></extra>",
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values,
            color=link_colors,
            hovertemplate="  %{source.label} → %{target.label}<extra></extra>",
        ),
    ))

    fig.update_layout(
        title=dict(text=f"<b>{title}</b>", font=dict(size=16, color="#2c3e50"),
                   x=0.5, xanchor="center", y=0.97),
        font=dict(family="Arial, sans-serif", size=11, color="#2c3e50"),
        paper_bgcolor="white",
        plot_bgcolor="white",
        width=1400,
        height=820,
        margin=dict(l=30, r=30, t=70, b=30),
    )

    # ── Export ────────────────────────────────────────────────────────────
    if output_path:
        try:
            pio.write_image(fig, output_path, format="png", scale=2,
                            width=1400, height=820)
            print(f"  Saved: {output_path}")
        except Exception as e:
            print(f"  WARNING: kaleido PNG export failed ({e}). "
                  "Falling back to HTML save.")
            html_path = output_path.replace(".png", ".html")
            fig.write_html(html_path)
            print(f"  Saved HTML fallback: {html_path}")
    try:
        fig.show()
    except Exception:
        pass
    return fig


# =============================================================================
# STEP 6b: Publication-ready composite figures (Bioinformatics style)
# =============================================================================

# Bioinformatics journal figure conventions:
#   single-column =  86 mm  = 3.39 in
#   1.5-column    = 120 mm  = 4.72 in
#   double-column = 178 mm  = 7.00 in
#   Fonts: Arial/Helvetica, panel labels 10-11 pt bold, axes 7-8 pt
#   DPI ≥ 300 for raster, prefer vector (PDF) when possible.

_BIOINF_RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
    "legend.title_fontsize": 7.5,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "savefig.dpi": 300,
    "figure.dpi": 300,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def _panel_label(ax, letter, x=-0.09, y=1.04, fig=None, anchor=None):
    """
    Add a bold panel letter (A, B, ...) on the top-left of an axes.
    If `fig` is given along with `anchor=(fx, fy)` in figure coordinates,
    the letter is placed in figure space (robust for polar/complex axes).
    """
    if fig is not None and anchor is not None:
        fig.text(anchor[0], anchor[1], letter,
                 fontsize=11, fontweight="bold", va="top", ha="left")
    else:
        ax.text(x, y, letter, transform=ax.transAxes,
                fontsize=11, fontweight="bold", va="bottom", ha="left")


def _prep_heatmap_data(omic_dfs, gene_set, max_genes=40, max_samples=60):
    """Shared heatmap data preparation for composite figures."""
    available = {n: df for n, df in omic_dfs.items()
                 if df is not None and not df.empty}
    if not available:
        return None, None

    all_cols = [set(df.columns) for df in available.values()]
    common = sorted(set.intersection(*all_cols)) if all_cols else []
    if not common:
        common = sorted(set.union(*all_cols))

    # Rank by expression (if available) then coverage
    expr = available.get("expression")
    if expr is not None:
        means = expr.mean(axis=0)
        common = [g for g in common if g in means.index]
        common = sorted(common, key=lambda g: -means.get(g, 0))
    if len(common) > max_genes:
        common = common[:max_genes]

    labels = {}
    for g in common:
        info = gene_set.get(g, {})
        lbl = info.get("label", g) if isinstance(info, dict) else g
        labels[g] = str(lbl).split(" (")[0]

    return available, (common, labels)


def _draw_heatmap_panel(ax, df, genes, gene_labels, cmap, center, title,
                        show_ylabels=True, max_samples=60, cbar_label="z-score"):
    """Draw a single omic heatmap panel (genes × samples, z-scored)."""
    cols = [g for g in genes if g in df.columns]
    if not cols:
        ax.set_title(title, fontsize=8.5, fontweight="bold")
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        return

    sub = df[cols].dropna(how="all")
    if sub.empty:
        ax.set_title(title, fontsize=8.5, fontweight="bold")
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="grey")
        ax.set_xticks([]); ax.set_yticks([])
        return

    if sub.shape[0] > max_samples:
        sub = sub.sample(n=max_samples, random_state=42)

    z = sub.apply(lambda x: (x - x.mean()) / x.std() if x.std() > 0 else x)
    z = z.clip(-3, 3)
    z = z.rename(columns={g: gene_labels.get(g, g) for g in cols})

    sns.heatmap(z.T, ax=ax, cmap=cmap, center=center,
                xticklabels=False,
                yticklabels=show_ylabels,
                cbar_kws={"shrink": 0.55, "label": cbar_label,
                          "pad": 0.02, "aspect": 20})
    cbar = ax.collections[0].colorbar
    if cbar is not None:
        cbar.ax.tick_params(labelsize=6)
        cbar.set_label(cbar_label, fontsize=7)
    ax.set_title(title, fontsize=8.5, fontweight="bold", pad=4)
    ax.set_xlabel(f"Samples (n={sub.shape[0]})", fontsize=7)
    if show_ylabels:
        ax.set_ylabel("Genes", fontsize=8)
        ax.tick_params(axis="y", labelsize=5.5)
    else:
        ax.set_ylabel("")


def _build_radar_metrics(omic_dfs, gene_set, random_omic_dfs, random_gene_set):
    """Compute normalized radar axes shared by both composite figures."""
    metric_labels, pheno_vals, base_vals = [], [], []

    layer_metrics = {
        "expression": ("Mean\nlog2(TPM+1)",
                       lambda df: np.log2(df.mean(axis=0).mean() + 1)),
        "cnv":        ("Mean\n|CNV log2R|",
                       lambda df: df.abs().mean(axis=0).mean()),
        "methylation":("Mean\nβ",
                       lambda df: df.mean(axis=0).mean()),
        "protein":    ("Mean\n|RPPA|",
                       lambda df: df.abs().mean(axis=0).mean()),
    }
    for ln, (lab, fn) in layer_metrics.items():
        p = omic_dfs.get(ln); r = random_omic_dfs.get(ln) if random_omic_dfs else None
        if p is not None and not p.empty:
            metric_labels.append(lab)
            pheno_vals.append(fn(p))
            base_vals.append(fn(r) if r is not None and not r.empty else 0)

    for ln in ["expression", "cnv"]:
        p = omic_dfs.get(ln); r = random_omic_dfs.get(ln) if random_omic_dfs else None
        if p is not None and not p.empty:
            metric_labels.append(f"{ln.capitalize()[:4]}.\nVariance")
            pheno_vals.append(p.var(axis=0).mean())
            base_vals.append(r.var(axis=0).mean() if r is not None and not r.empty else 0)

    n_p = max(len(gene_set), 1)
    n_r = max(len(random_gene_set), 1) if random_gene_set else 1
    for ln in ["expression", "protein"]:
        p = omic_dfs.get(ln); r = random_omic_dfs.get(ln) if random_omic_dfs else None
        if p is not None and not p.empty:
            metric_labels.append(f"{ln.capitalize()[:4]}.\nCoverage")
            pheno_vals.append(p.shape[1] / n_p)
            base_vals.append(r.shape[1] / n_r if r is not None and not r.empty else 0)

    if len(metric_labels) < 3:
        return None

    maxes = [max(abs(a), abs(b), 0.001) for a, b in zip(pheno_vals, base_vals)]
    pn = [a / m for a, m in zip(pheno_vals, maxes)]
    bn = [b / m for b, m in zip(base_vals, maxes)]
    return metric_labels, pn, bn


def _draw_radar_panel(ax, radar_data):
    """Draw the radar chart onto a polar axes."""
    metric_labels, pn, bn = radar_data
    N = len(metric_labels)
    angles = np.linspace(0, 2 * np.pi, N, endpoint=False).tolist()
    angles_cl = angles + angles[:1]
    pn_cl = pn + pn[:1]; bn_cl = bn + bn[:1]

    ax.fill(angles_cl, pn_cl, alpha=0.28, color="#C0392B")
    ax.plot(angles_cl, pn_cl, "o-", linewidth=1.4, color="#C0392B",
            markersize=3.5, label="Phenotype (semantic)")
    ax.fill(angles_cl, bn_cl, alpha=0.18, color="#2980B9")
    ax.plot(angles_cl, bn_cl, "s--", linewidth=1.2, color="#2980B9",
            markersize=3.5, label="Random baseline")

    ax.set_xticks(angles)
    ax.set_xticklabels(metric_labels, fontsize=6.5)
    ax.set_ylim(0, 1.15)
    ax.set_yticks([0.25, 0.5, 0.75, 1.0])
    ax.set_yticklabels(["0.25", "0.50", "0.75", "1.00"], fontsize=5.5,
                       color="grey")
    ax.tick_params(axis="x", pad=1)
    ax.grid(linewidth=0.4, alpha=0.6)
    ax.spines["polar"].set_linewidth(0.5)
    ax.legend(loc="upper right", bbox_to_anchor=(1.32, 1.08),
              fontsize=6.5, frameon=False, handlelength=1.5)


def _build_bootstrap_panel_data(omic_dfs, gene_set, valid_genes,
                                 random_omic_dfs=None,
                                 n_bootstrap=N_BOOTSTRAP, boot_seed=SEED):
    """
    Build per-gene null pools for each active omic layer.

    Preferred path (random_omic_dfs provided): null pool = per-gene values
    from random_omic_dfs columns directly. This is the correct approach
    when the phenotype omic_dfs were extracted for phenotype genes only
    (their columns don't include random genes, so resampling from
    valid_genes and filtering by omic_dfs.columns yields empty cols).

    Legacy fallback path (random_omic_dfs=None): resample from valid_genes
    and filter by omic_dfs.columns — only works if omic_dfs contains data
    for a large gene pool (not just phenotype).

    Returns a dict {layer: {pheno_vals, null_pool, label, color}} where
    pheno_vals is an array of per-gene values for the phenotype set.
    """
    layer_cfg = {
        "expression":  ("Expression log2(TPM+1)", "#C0392B",
                        lambda df, cols: np.log2(df[cols].mean(axis=0) + 1).dropna().values
                        if len(cols) else np.array([])),
        "cnv":         ("CNV |log2R|",             "#2980B9",
                        lambda df, cols: df[cols].abs().mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
        "methylation": ("Methylation β",            "#27AE60",
                        lambda df, cols: df[cols].mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
        "protein":     ("Protein |RPPA|",           "#8E44AD",
                        lambda df, cols: df[cols].abs().mean(axis=0).dropna().values
                        if len(cols) else np.array([])),
    }

    pheno_ids = set(gene_set.keys())
    pool = [g for g in valid_genes if g and g not in pheno_ids]
    n_sample = len(gene_set)
    rng = np.random.default_rng(boot_seed)

    result = {}
    for ln, (label, color, fn) in layer_cfg.items():
        df = omic_dfs.get(ln)
        if df is None or df.empty:
            continue
        pheno_cols = [g for g in pheno_ids if g in df.columns]
        if not pheno_cols:
            continue
        pheno_vals = fn(df, pheno_cols)

        # --- Preferred path: null pool from random_omic_dfs ---
        if random_omic_dfs is not None:
            rand_df = random_omic_dfs.get(ln)
            if rand_df is not None and not rand_df.empty:
                rand_cols = list(rand_df.columns)
                null_pool = fn(rand_df, rand_cols)
            else:
                null_pool = np.array([])
        else:
            # --- Legacy fallback: resample from valid_genes (only works if
            # omic_dfs contains data for genes beyond phenotype set) ---
            null_per_resample = []
            for _ in range(n_bootstrap):
                samp = rng.choice(pool, size=min(n_sample, len(pool)),
                                  replace=False).tolist()
                cols = [g for g in samp if g in df.columns]
                if cols:
                    vals = fn(df, cols)
                    if len(vals):
                        null_per_resample.append(vals)
            null_pool = (np.concatenate(null_per_resample)
                         if null_per_resample else np.array([]))

        result[ln] = {"pheno_vals": pheno_vals, "null_pool": null_pool,
                      "label": label, "color": color}
    return result


def _draw_bootstrap_panel(ax, boot_data: dict):
    """
    Draw a side-by-side violin panel per omic layer: null (grey) vs
    phenotype (colored), each built from per-gene values.

    A Mann-Whitney U p-value is annotated above each layer group.
    """
    if not boot_data:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="grey", fontsize=7)
        ax.axis("off")
        return

    layers = list(boot_data.keys())
    n = len(layers)
    group_w = 1.0
    offset  = 0.24

    all_top = []
    for i, ln in enumerate(layers):
        d = boot_data[ln]
        null  = d["null_pool"]
        pheno = d["pheno_vals"]
        color = d["color"]

        xc = i * group_w
        positions, datasets, facecolors, edges = [], [], [], []
        if len(null) > 1:
            positions.append(xc - offset); datasets.append(null)
            facecolors.append("#D5D8DC"); edges.append("#95A5A6")
        if len(pheno) > 1:
            positions.append(xc + offset); datasets.append(pheno)
            facecolors.append(color); edges.append(color)

        if datasets:
            parts = ax.violinplot(datasets, positions=positions, widths=0.42,
                                  showmedians=True, showextrema=False)
            for pc, fc, ec in zip(parts["bodies"], facecolors, edges):
                pc.set_facecolor(fc); pc.set_alpha(0.55)
                pc.set_edgecolor(ec); pc.set_linewidth(0.6)
            parts["cmedians"].set_color("#333"); parts["cmedians"].set_linewidth(0.8)

        # Mann-Whitney per-gene, annotate
        if len(null) >= 5 and len(pheno) >= 5:
            u_stat, p_val = stats.mannwhitneyu(pheno, null, alternative="two-sided")
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else \
                  "*" if p_val < 0.05 else "ns"
            y_top = float(np.nanpercentile(np.concatenate([null, pheno]), 98))
            all_top.append(y_top)
            ax.plot([xc - offset, xc + offset], [y_top * 1.02, y_top * 1.02],
                    color="#555", lw=0.5)
            ax.text(xc, y_top * 1.04, sig,
                    ha="center", va="bottom", fontsize=6, color="#333",
                    fontweight="bold")

    ax.set_xticks([i * group_w for i in range(n)])
    ax.set_xticklabels([boot_data[ln]["label"] for ln in layers],
                       fontsize=6, rotation=0)
    ax.set_ylabel("Per-gene signal", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=6)

    # Legend
    handles = [
        mpatches.Patch(facecolor="#D5D8DC", edgecolor="#95A5A6",
                       alpha=0.7, label="Null (bootstrap)"),
        mpatches.Patch(facecolor="#C0392B", edgecolor="#C0392B",
                       alpha=0.55, label="Phenotype"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=5.5,
              frameon=False, handlelength=1.2)


def _draw_boxplot_panel(ax, boot_data: dict):
    """
    Boxplot variant of _draw_bootstrap_panel. Side-by-side per layer:
    null (grey) vs phenotype (colored). Mann-Whitney significance annotated.
    """
    if not boot_data:
        ax.text(0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, color="grey", fontsize=7)
        ax.axis("off")
        return

    layers = list(boot_data.keys())
    n = len(layers)
    group_w = 1.0
    offset  = 0.24

    for i, ln in enumerate(layers):
        d = boot_data[ln]
        null  = d["null_pool"]
        pheno = d["pheno_vals"]
        color = d["color"]
        xc = i * group_w

        if len(null) > 1:
            bp1 = ax.boxplot([null], positions=[xc - offset], widths=0.38,
                             patch_artist=True, showfliers=False,
                             medianprops=dict(color="#333", linewidth=0.8),
                             whiskerprops=dict(color="#95A5A6", linewidth=0.6),
                             capprops=dict(color="#95A5A6", linewidth=0.6),
                             boxprops=dict(linewidth=0.6))
            for patch in bp1["boxes"]:
                patch.set_facecolor("#D5D8DC"); patch.set_alpha(0.65)
                patch.set_edgecolor("#95A5A6")
        if len(pheno) > 1:
            bp2 = ax.boxplot([pheno], positions=[xc + offset], widths=0.38,
                             patch_artist=True, showfliers=False,
                             medianprops=dict(color="#333", linewidth=0.8),
                             whiskerprops=dict(color=color, linewidth=0.6),
                             capprops=dict(color=color, linewidth=0.6),
                             boxprops=dict(linewidth=0.6))
            for patch in bp2["boxes"]:
                patch.set_facecolor(color); patch.set_alpha(0.55)
                patch.set_edgecolor(color)

        if len(null) >= 5 and len(pheno) >= 5:
            u_stat, p_val = stats.mannwhitneyu(pheno, null, alternative="two-sided")
            sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else \
                  "*" if p_val < 0.05 else "ns"
            y_top = float(np.nanpercentile(np.concatenate([null, pheno]), 98))
            ax.plot([xc - offset, xc + offset], [y_top * 1.02, y_top * 1.02],
                    color="#555", lw=0.5)
            ax.text(xc, y_top * 1.04, sig,
                    ha="center", va="bottom", fontsize=6, color="#333",
                    fontweight="bold")

    ax.set_xticks([i * group_w for i in range(n)])
    ax.set_xticklabels([boot_data[ln]["label"] for ln in layers],
                       fontsize=6, rotation=0)
    ax.set_xlim(-0.6, (n - 1) * group_w + 0.6)
    ax.set_ylabel("Per-gene signal", fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    ax.tick_params(axis="y", labelsize=6)

    handles = [
        mpatches.Patch(facecolor="#D5D8DC", edgecolor="#95A5A6",
                       alpha=0.7, label="Random baseline"),
        mpatches.Patch(facecolor="#C0392B", edgecolor="#C0392B",
                       alpha=0.55, label="Phenotype"),
    ]
    ax.legend(handles=handles, loc="upper right", fontsize=5.5,
              frameon=False, handlelength=1.2)


def _draw_subgraph_panel(ax, traversal, seed_node, mean_expression,
                          max_gene_labels=18, max_nodes=90):
    """Draw a compact semantic subgraph suited for a publication panel."""
    G = nx.DiGraph()
    seed_key = seed_node["_key"]
    seed_label = seed_node.get("label", seed_key)
    G.add_node(seed_key, label=seed_label, node_type="phenotype", mean_expr=0)

    genes = traversal["genes"]; diseases = traversal["diseases"]
    proteins = traversal["proteins"]; pathways = traversal["pathways"]

    # Rank genes by expression for label prioritization
    if mean_expression is not None:
        sorted_g = sorted(genes.items(),
                          key=lambda kv: -float(mean_expression.get(kv[0], 0)))
    else:
        sorted_g = list(genes.items())

    top_label_keys = set()
    for eid, info in sorted_g[:max_nodes]:
        kg_key = info.get("kg_key") if isinstance(info, dict) else eid
        short = (info.get("label", eid) if isinstance(info, dict) else eid)
        short = str(short).split(" (")[0]
        mv = float(mean_expression.get(eid, 0)) if mean_expression is not None else 0
        G.add_node(kg_key, label=short, node_type="gene",
                   mean_expr=mv, entrez_id=eid)
        if len(top_label_keys) < max_gene_labels:
            top_label_keys.add(kg_key)

    # Prioritize breast-specific diseases
    d_items = sorted(diseases.items(),
                     key=lambda kv: (not kv[1].get("is_breast_specific"),
                                     kv[1].get("label", "")))
    for dk, dinfo in d_items[:8]:
        lbl = (dinfo.get("label") or dk)[:18]
        G.add_node(dk, label=lbl, node_type="disease", mean_expr=0)

    for pk, pinfo in list(pathways.items())[:6]:
        lbl = (pinfo.get("label") or pk)[:14]
        G.add_node(pk, label=lbl, node_type="pathway", mean_expr=0)
    for prk, prinfo in list(proteins.items())[:8]:
        lbl = (prinfo.get("label") or prk)[:10]
        G.add_node(prk, label=lbl, node_type="protein", mean_expr=0)

    in_G = set(G.nodes())
    for src, tgt, pred in traversal["subgraph_edges"]:
        if src in in_G and tgt in in_G:
            G.add_edge(src, tgt, predicate=pred)
    if len(G.edges()) == 0:
        for dk in diseases:
            if dk in in_G: G.add_edge(seed_key, dk, predicate="has phenotype")
        for eid, info in genes.items():
            gk = info.get("kg_key") if isinstance(info, dict) else eid
            if gk in in_G: G.add_edge(seed_key, gk, predicate="linked")

    if len(G.nodes()) == 0:
        ax.text(0.5, 0.5, "No subgraph", ha="center", va="center",
                transform=ax.transAxes, color="grey")
        ax.axis("off")
        return

    pos = nx.spring_layout(G, seed=42,
                           k=2.2 / np.sqrt(max(len(G.nodes()), 1)),
                           iterations=80)

    type_colors = {"phenotype": "#C0392B", "gene": "#2980B9",
                   "protein": "#27AE60", "disease": "#E67E22",
                   "pathway": "#8E44AD"}
    node_types = nx.get_node_attributes(G, "node_type")
    node_colors = [type_colors.get(node_types.get(n, "gene"), "#BBBBBB")
                   for n in G.nodes()]
    sizes = []
    for n in G.nodes():
        nt = node_types.get(n, "gene")
        sizes.append(260 if nt == "phenotype" else
                     90 if nt == "gene" else 55)

    if mean_expression is not None:
        mes = nx.get_node_attributes(G, "mean_expr")
        vals = [v for v in mes.values() if v > 0]
        if vals:
            vmin, vmax = np.percentile(vals, [5, 95])
            norm = Normalize(vmin=max(vmin, 0.01), vmax=max(vmax, 1))
            for i, n in enumerate(G.nodes()):
                if node_types.get(n) == "gene" and mes.get(n, 0) > 0:
                    node_colors[i] = plt.cm.Blues(0.35 + 0.6 * norm(mes[n]))

    nx.draw_networkx_edges(G, pos, ax=ax, edge_color="#B8B8B8",
                           alpha=0.45, width=0.35, arrows=False)
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=sizes, alpha=0.9,
                           edgecolors="white", linewidths=0.4)

    labels = {}
    node_labels = nx.get_node_attributes(G, "label")
    for n in G.nodes():
        nt = node_types.get(n, "")
        if nt in ("phenotype", "disease") or n in top_label_keys:
            labels[n] = str(node_labels.get(n, n))[:16]
    nx.draw_networkx_labels(G, pos, labels, font_size=5, ax=ax)

    legend_handles = [mpatches.Patch(color=c, label=t.capitalize())
                      for t, c in type_colors.items()]
    ax.legend(handles=legend_handles, loc="lower left",
              fontsize=6.5, frameon=False, handlelength=1,
              ncol=1, borderpad=0.2)
    ax.axis("off")


def _draw_sankey_panel(ax, traversal, omic_dfs, gene_set, seed_node,
                        top_n_genes=14):
    """
    Lightweight matplotlib Sankey-style flow diagram (phenotype → disease
    class → gene → omic). Uses rectangles + quadratic Bezier bands so the
    figure stays vector-friendly and integrates cleanly with the composite.
    """
    ax.axis("off")
    ax.set_xlim(-0.02, 1.05); ax.set_ylim(-0.02, 1.02)

    diseases = traversal.get("diseases", {})
    expr = omic_dfs.get("expression")
    if expr is not None and not expr.empty:
        mean_expr = expr.mean(axis=0).sort_values(ascending=False)
        top_genes = [g for g in mean_expr.index if g in gene_set][:top_n_genes]
    else:
        mean_expr = pd.Series(dtype=float)
        top_genes = list(gene_set.keys())[:top_n_genes]

    active_layers = [l for l in ["expression", "cnv", "methylation", "protein"]
                     if omic_dfs.get(l) is not None and not omic_dfs[l].empty]
    layer_disp = {"expression": "Expr", "cnv": "CNV",
                  "methylation": "Methyl", "protein": "Protein"}
    layer_color = {"expression": "#C0392B", "cnv": "#2980B9",
                   "methylation": "#27AE60", "protein": "#8E44AD"}

    # Column 0: Phenotype (single)
    # Column 1: Disease classes (breast-specific vs other vs direct)
    # Column 2: Top genes
    # Column 3: Omic layers
    col_x = [0.08, 0.35, 0.65, 0.93]

    # Build disease class nodes (3 fixed buckets)
    n_breast = sum(1 for d in diseases.values() if d.get("is_breast_specific"))
    n_other  = sum(1 for d in diseases.values() if not d.get("is_breast_specific"))
    bucket_items = [
        ("direct",  f"Direct\nassoc.",       "#F39C12"),
        ("breast",  f"Breast-specific\ndiseases (n={n_breast})", "#D35400"),
        ("other",   f"Other\ndiseases (n={n_other})",   "#F5B041"),
    ]

    # Node layout: evenly spaced vertically within column
    def column_positions(n, top_margin=0.05, bottom_margin=0.05):
        if n == 1: return [0.5]
        span = 1 - top_margin - bottom_margin
        return [1 - top_margin - i * span / (n - 1) for i in range(n)]

    bucket_y = column_positions(3)
    bucket_pos = {bid: bucket_y[i] for i, (bid, *_rest) in enumerate(bucket_items)}

    gene_y = column_positions(len(top_genes)) if top_genes else []
    layer_y = column_positions(len(active_layers)) if active_layers else []

    # Rectangle sizes
    W = 0.085
    H_seed = 0.14; H_bucket = 0.09; H_gene = 0.048; H_layer = 0.085

    def rect(ax, x, y, w, h, color, label=None, label_color="black",
             fontsize=6.5, fontweight="normal"):
        ax.add_patch(mpatches.FancyBboxPatch(
            (x - w / 2, y - h / 2), w, h,
            boxstyle="round,pad=0.002,rounding_size=0.008",
            linewidth=0.4, edgecolor="white", facecolor=color, alpha=0.92))
        if label is not None:
            ax.text(x, y, label, ha="center", va="center",
                    fontsize=fontsize, color=label_color,
                    fontweight=fontweight, zorder=5)

    def flow(ax, x0, y0, x1, y1, color, alpha=0.25, lw_scale=1.0):
        """Bezier band between two nodes."""
        import matplotlib.path as mpath
        import matplotlib.patches as mpatch
        Path = mpath.Path
        verts = [(x0, y0),
                 ((x0 + x1) / 2, y0),
                 ((x0 + x1) / 2, y1),
                 (x1, y1)]
        codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4]
        p = mpatch.PathPatch(Path(verts, codes), fc="none",
                             ec=color, alpha=alpha,
                             linewidth=1.0 * lw_scale)
        ax.add_patch(p)

    # --- Draw nodes ---
    # Phenotype
    seed_lbl = (seed_node.get("label", "Phenotype"))[:14]
    rect(ax, col_x[0], 0.5, W * 1.05, H_seed, "#C0392B",
         label=seed_lbl, label_color="white",
         fontsize=6.5, fontweight="bold")

    # Disease buckets
    bucket_color = {bid: col for bid, _, col in bucket_items}
    for (bid, blab, col), y in zip(bucket_items, bucket_y):
        rect(ax, col_x[1], y, W, H_bucket, col, label=blab,
             label_color="white", fontsize=6, fontweight="bold")

    # Genes — color by expression
    gene_color_map = {}
    if len(mean_expr):
        valid = mean_expr[mean_expr > 0]
        if len(valid):
            vmin, vmax = np.percentile(valid, [5, 95])
            norm = Normalize(vmin=max(vmin, 0.01), vmax=max(vmax, 1))
            for g in top_genes:
                ev = mean_expr.get(g, 0)
                if ev > 0:
                    gene_color_map[g] = plt.cm.Blues(0.4 + 0.55 * norm(ev))
                else:
                    gene_color_map[g] = "#BDC3C7"
    else:
        for g in top_genes: gene_color_map[g] = "#BDC3C7"

    gene_pos = {}
    for g, y in zip(top_genes, gene_y):
        info = gene_set.get(g, {})
        lbl = info.get("label", g) if isinstance(info, dict) else g
        lbl = str(lbl).split(" (")[0][:8]
        rect(ax, col_x[2], y, W * 0.9, H_gene, gene_color_map.get(g, "#BDC3C7"),
             label=lbl, label_color="black", fontsize=5.2,
             fontweight="normal")
        gene_pos[g] = y

    # Omic layers
    layer_pos = {}
    for ln, y in zip(active_layers, layer_y):
        df = omic_dfs[ln]
        covered = sum(1 for g in top_genes if g in df.columns)
        lbl = f"{layer_disp[ln]}\n({covered}/{len(top_genes)})"
        rect(ax, col_x[3], y, W * 1.1, H_layer, layer_color[ln],
             label=lbl, label_color="white", fontsize=5.8,
             fontweight="bold")
        layer_pos[ln] = y

    # --- Draw flows ---
    # Phenotype → disease buckets
    for bid, y in bucket_pos.items():
        flow(ax, col_x[0] + W / 2, 0.5, col_x[1] - W / 2, y,
             bucket_color[bid], alpha=0.35, lw_scale=3.0)

    # Disease bucket → gene
    for g, y in gene_pos.items():
        info = gene_set.get(g, {})
        via = info.get("via_disease_key") if isinstance(info, dict) else None
        if via is None:
            bid = "direct"
        else:
            d = diseases.get(via, {})
            bid = "breast" if d.get("is_breast_specific") else "other"
        flow(ax, col_x[1] + W / 2, bucket_pos[bid],
             col_x[2] - W * 0.45, y,
             bucket_color[bid], alpha=0.30, lw_scale=1.3)

    # Gene → omic layer (coverage-weighted)
    for ln, y_l in layer_pos.items():
        df = omic_dfs[ln]
        for g, y_g in gene_pos.items():
            if g in df.columns and not df[g].isna().all():
                flow(ax, col_x[2] + W * 0.45, y_g,
                     col_x[3] - W * 0.55, y_l,
                     layer_color[ln], alpha=0.22, lw_scale=1.0)

    # Column headers (kept inside the axes, not overlapping the pane title)
    headers = ["Phenotype", "Disease ctx.", "Top genes", "Omic layers"]
    for x, h in zip(col_x, headers):
        ax.text(x, 0.97, h, ha="center", va="top",
                fontsize=6.5, fontweight="bold", color="#2C3E50")


def plot_composite_figure_2panel(omic_dfs, gene_set,
                                  valid_genes=None,
                                  random_omic_dfs=None,
                                  seed_node=None,
                                  title=None,
                                  output_path=None):
    """
    Bioinformatics-style 2-panel figure:
        A — Multi-omic heatmap (4 subpanels stacked horizontally)
        B — Bootstrap null comparison (phenotype dot vs violin null distribution)

    Panel B replaces the former radar chart. The bootstrap null is built by
    resampling N_BOOTSTRAP random gene sets from valid_genes (the expression
    index), computing the same per-layer metric each time, and displaying the
    resulting distribution as a violin+IQR box. The phenotype value is the
    colored dot. Z-scores are annotated above each layer column.

    Designed for double-column layout (178 mm wide).
    """
    prep = _prep_heatmap_data(omic_dfs, gene_set, max_genes=35, max_samples=60)
    if prep is None or prep[0] is None:
        print("  [composite-2p] No omic data, skipping.")
        return None
    available, (common_genes, gene_labels) = prep

    vgenes = valid_genes or []
    boot_data = _build_bootstrap_panel_data(omic_dfs, gene_set, vgenes,
                                             random_omic_dfs=random_omic_dfs)

    with plt.rc_context(_BIOINF_RC):
        fig = plt.figure(figsize=(7.2, 4.9))

        outer = gridspec.GridSpec(1, 2, width_ratios=[1.7, 1.0],
                                  wspace=0.32, left=0.07, right=0.98,
                                  top=0.83, bottom=0.12)

        # ── Panel A: heatmap strip ──
        n_layers = len(available)
        inner_a = gridspec.GridSpecFromSubplotSpec(
            1, n_layers, subplot_spec=outer[0],
            wspace=0.18)

        layer_cmaps = {"expression": "Reds", "cnv": "RdBu_r",
                       "methylation": "YlGn", "protein": "PuBu"}
        layer_titles = {"expression": "Expression\n(TPM)",
                        "cnv": "CNV\n(log2 ratio)",
                        "methylation": "Methylation\n(β)",
                        "protein": "Protein\n(RPPA)"}

        for i, (ln, df) in enumerate(available.items()):
            ax = fig.add_subplot(inner_a[0, i])
            _draw_heatmap_panel(
                ax, df, common_genes, gene_labels,
                cmap=layer_cmaps.get(ln, "viridis"),
                center=0 if ln == "cnv" else None,
                title=layer_titles.get(ln, ln),
                show_ylabels=(i == 0),
                max_samples=60)

        # ── Panel B: bootstrap null comparison ──
        ax_b = fig.add_subplot(outer[1])
        _draw_bootstrap_panel(ax_b, boot_data)
        ax_b.set_title("Multi-omic signal\nvs bootstrap null",
                       fontsize=8.5, fontweight="bold", pad=6)

        _panel_label(None, "A", fig=fig, anchor=(0.015, 0.93))
        _panel_label(None, "B", fig=fig, anchor=(0.595, 0.93))

        if title is None:
            seed_lab = seed_node.get("label", "Phenotype") if seed_node else "Phenotype"
            title = f"Phenotype-anchored multi-omic signature — {seed_lab}"
        fig.suptitle(title, fontsize=10.5, fontweight="bold", y=0.97)

        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            pdf_path = output_path.rsplit(".", 1)[0] + ".pdf"
            try:
                fig.savefig(pdf_path, bbox_inches="tight")
            except Exception:
                pass
            print(f"  Saved: {output_path}")
        plt.show()
        return fig


def plot_composite_figure_4panel(traversal, omic_dfs, gene_set, seed_node,
                                  valid_genes=None,
                                  random_omic_dfs=None,
                                  mean_expression=None,
                                  title=None,
                                  output_path=None):
    """
    Bioinformatics-style 4-panel figure:
        A — Semantic subgraph (KG traversal topology)
        B — Multi-omic heatmap (4 layers, z-scored)
        C — Bootstrap null comparison (replaces radar; robust across runs)
        D — Sankey-like flow (phenotype → disease → gene → omic)

    Panel C uses _build_bootstrap_panel_data / _draw_bootstrap_panel:
    N_BOOTSTRAP random gene sets are drawn from valid_genes (columns already
    present in the omic DataFrames — no DB queries), producing a violin null
    distribution per layer. The phenotype dot is overlaid with a z-score.

    Designed for double-column layout (178 mm wide).
    """
    prep = _prep_heatmap_data(omic_dfs, gene_set, max_genes=28, max_samples=60)
    if prep is None or prep[0] is None:
        print("  [composite-4p] No omic data, skipping.")
        return None
    available, (common_genes, gene_labels) = prep

    vgenes = valid_genes or []
    boot_data = _build_bootstrap_panel_data(omic_dfs, gene_set, vgenes,
                                             random_omic_dfs=random_omic_dfs)

    with plt.rc_context(_BIOINF_RC):
        fig = plt.figure(figsize=(7.2, 8.8))

        outer = gridspec.GridSpec(
            2, 2,
            width_ratios=[1.15, 1.0],
            height_ratios=[1.0, 1.0],
            hspace=0.40, wspace=0.24,
            left=0.06, right=0.985,
            top=0.92, bottom=0.04)

        # ── Panel A: semantic subgraph (top-left) ──
        ax_a = fig.add_subplot(outer[0, 0])
        _draw_subgraph_panel(ax_a, traversal, seed_node, mean_expression,
                              max_gene_labels=12, max_nodes=70)
        ax_a.set_title("Semantic subgraph (KG traversal)",
                       fontsize=9, fontweight="bold", pad=4)

        # ── Panel B: heatmap strip (top-right) ──
        n_layers = len(available)
        inner_b = gridspec.GridSpecFromSubplotSpec(
            1, n_layers, subplot_spec=outer[0, 1],
            wspace=0.22)

        layer_cmaps = {"expression": "Reds", "cnv": "RdBu_r",
                       "methylation": "YlGn", "protein": "PuBu"}
        layer_titles = {"expression": "Expr.", "cnv": "CNV",
                        "methylation": "Methyl.", "protein": "Prot."}

        for i, (ln, df) in enumerate(available.items()):
            ax = fig.add_subplot(inner_b[0, i])
            _draw_heatmap_panel(
                ax, df, common_genes, gene_labels,
                cmap=layer_cmaps.get(ln, "viridis"),
                center=0 if ln == "cnv" else None,
                title=layer_titles.get(ln, ln),
                show_ylabels=(i == 0),
                max_samples=60)
            ax.set_xlabel("")
        # Shared xlabel for the whole strip
        bpos = outer[0, 1].get_position(fig)
        fig.text(bpos.x0 + bpos.width / 2,
                 bpos.y0 - 0.005,
                 "Samples (n=60)", ha="center", fontsize=7.2)

        # ── Panel C: bootstrap null comparison (bottom-left) ──
        ax_c = fig.add_subplot(outer[1, 0])
        _draw_bootstrap_panel(ax_c, boot_data)
        ax_c.set_title("Multi-omic signal vs bootstrap null",
                       fontsize=9, fontweight="bold", pad=4)

        # ── Panel D: Sankey flow (bottom-right) ──
        ax_d = fig.add_subplot(outer[1, 1])
        _draw_sankey_panel(ax_d, traversal, omic_dfs, gene_set, seed_node,
                            top_n_genes=12)
        ax_d.set_title("Trans-omic flow diagram",
                       fontsize=9, fontweight="bold", pad=4)

        # Panel letters — placed in figure coords above each quadrant so they
        # never collide with titles. Derive positions from gridspec bboxes.
        def _fig_anchor(sp, dx=-0.01, dy=0.012):
            p = sp.get_position(fig)
            return (p.x0 + dx, p.y1 + dy)

        _panel_label(None, "A", fig=fig, anchor=_fig_anchor(outer[0, 0]))
        _panel_label(None, "B", fig=fig, anchor=_fig_anchor(outer[0, 1]))
        _panel_label(None, "C", fig=fig, anchor=_fig_anchor(outer[1, 0]))
        _panel_label(None, "D", fig=fig, anchor=_fig_anchor(outer[1, 1]))

        if title is None:
            seed_lab = seed_node.get("label", "Phenotype")
            title = (f"Phenotype-anchored trans-omic signature "
                     f"(HP: {seed_lab})")
        fig.suptitle(title, fontsize=10.5, fontweight="bold", y=0.975)

        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            pdf_path = output_path.rsplit(".", 1)[0] + ".pdf"
            try:
                fig.savefig(pdf_path, bbox_inches="tight")
            except Exception:
                pass
            print(f"  Saved: {output_path}")
        plt.show()
        return fig


def _plot_composite_3panel_stacked_impl(traversal, omic_dfs, gene_set, seed_node,
                                         valid_genes, random_omic_dfs,
                                         mean_expression, title, output_path,
                                         null_panel_drawer, null_panel_label):
    """
    Shared impl for 3-panel stacked publication figure.
        A (top)    — Semantic subgraph (wide horizontal)
        B (middle) — Multi-omic heatmap strip (4 subpanels)
        C (bottom) — Null-vs-phenotype comparison (violin OR boxplot)
    """
    prep = _prep_heatmap_data(omic_dfs, gene_set, max_genes=30, max_samples=60)
    if prep is None or prep[0] is None:
        print("  [composite-3p] No omic data, skipping.")
        return None
    available, (common_genes, gene_labels) = prep

    vgenes = valid_genes or []
    boot_data = _build_bootstrap_panel_data(omic_dfs, gene_set, vgenes,
                                             random_omic_dfs=random_omic_dfs)

    with plt.rc_context(_BIOINF_RC):
        fig = plt.figure(figsize=(7.0, 10.5))

        outer = gridspec.GridSpec(
            3, 1,
            height_ratios=[1.1, 0.85, 0.85],
            hspace=0.42,
            left=0.08, right=0.97,
            top=0.95, bottom=0.05)

        # ── Panel A: semantic subgraph (wide top) ──
        ax_a = fig.add_subplot(outer[0, 0])
        _draw_subgraph_panel(ax_a, traversal, seed_node, mean_expression,
                              max_gene_labels=16, max_nodes=90)
        ax_a.set_title("Semantic subgraph (KG traversal)",
                       fontsize=9, fontweight="bold", pad=4)

        # ── Panel B: multi-omic heatmap strip (4 subpanels) ──
        n_layers = len(available)
        inner_b = gridspec.GridSpecFromSubplotSpec(
            1, n_layers, subplot_spec=outer[1, 0],
            wspace=0.22)

        layer_cmaps = {"expression": "Reds", "cnv": "RdBu_r",
                       "methylation": "YlGn", "protein": "PuBu"}
        layer_titles = {"expression": "Expression\n(TPM)",
                        "cnv": "CNV\n(log2 ratio)",
                        "methylation": "Methylation\n(β)",
                        "protein": "Protein\n(RPPA)"}

        for i, (ln, df) in enumerate(available.items()):
            ax = fig.add_subplot(inner_b[0, i])
            _draw_heatmap_panel(
                ax, df, common_genes, gene_labels,
                cmap=layer_cmaps.get(ln, "viridis"),
                center=0 if ln == "cnv" else None,
                title=layer_titles.get(ln, ln),
                show_ylabels=(i == 0),
                max_samples=60)

        # ── Panel C: null comparison (violin or boxplot) ──
        ax_c = fig.add_subplot(outer[2, 0])
        null_panel_drawer(ax_c, boot_data)
        ax_c.set_title(null_panel_label,
                       fontsize=9, fontweight="bold", pad=4)

        def _fig_anchor(sp, dx=-0.02, dy=0.012):
            p = sp.get_position(fig)
            return (p.x0 + dx, p.y1 + dy)

        _panel_label(None, "A", fig=fig, anchor=_fig_anchor(outer[0, 0]))
        _panel_label(None, "B", fig=fig, anchor=_fig_anchor(outer[1, 0]))
        _panel_label(None, "C", fig=fig, anchor=_fig_anchor(outer[2, 0]))

        if title is None:
            seed_lab = seed_node.get("label", "Phenotype") if seed_node else "Phenotype"
            title = (f"Phenotype-anchored trans-omic signature — {seed_lab}")
        fig.suptitle(title, fontsize=10.5, fontweight="bold", y=0.985)

        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches="tight")
            pdf_path = output_path.rsplit(".", 1)[0] + ".pdf"
            try:
                fig.savefig(pdf_path, bbox_inches="tight")
            except Exception:
                pass
            print(f"  Saved: {output_path}")
        plt.show()
        return fig


def plot_composite_figure_3panel_stacked(traversal, omic_dfs, gene_set, seed_node,
                                          valid_genes=None,
                                          random_omic_dfs=None,
                                          mean_expression=None,
                                          title=None,
                                          output_path=None):
    """3-panel stacked publication figure (violin null, panel C)."""
    return _plot_composite_3panel_stacked_impl(
        traversal, omic_dfs, gene_set, seed_node,
        valid_genes, random_omic_dfs, mean_expression, title, output_path,
        null_panel_drawer=_draw_bootstrap_panel,
        null_panel_label="Multi-omic signal vs bootstrap null")


def plot_composite_figure_3panel_stacked_box(traversal, omic_dfs, gene_set, seed_node,
                                              valid_genes=None,
                                              random_omic_dfs=None,
                                              mean_expression=None,
                                              title=None,
                                              output_path=None):
    """3-panel stacked publication figure (boxplot null, panel C)."""
    return _plot_composite_3panel_stacked_impl(
        traversal, omic_dfs, gene_set, seed_node,
        valid_genes, random_omic_dfs, mean_expression, title, output_path,
        null_panel_drawer=_draw_boxplot_panel,
        null_panel_label="Multi-omic signal — phenotype vs random baseline")


# =============================================================================
# STEP 7: Summary statistics
# =============================================================================

def compute_summary_stats(omic_dfs: Dict[str, pd.DataFrame],
                           gene_set: Dict) -> pd.DataFrame:
    """
    Compute summary statistics per omic layer for the phenotype gene set.
    """
    rows = []
    for layer_name, df in omic_dfs.items():
        if df is None or df.empty:
            continue

        gene_coverage = len([g for g in gene_set if g in df.columns])
        values = df.values.flatten()
        values = values[~np.isnan(values)]

        rows.append({
            "omic_layer": layer_name,
            "n_genes_covered": gene_coverage,
            "n_samples": df.shape[0],
            "mean": np.mean(values) if len(values) > 0 else np.nan,
            "median": np.median(values) if len(values) > 0 else np.nan,
            "std": np.std(values) if len(values) > 0 else np.nan,
            "min": np.min(values) if len(values) > 0 else np.nan,
            "max": np.max(values) if len(values) > 0 else np.nan,
            "pct_nonzero": (np.sum(values != 0) / len(values) * 100)
            if len(values) > 0 else 0,
        })

    return pd.DataFrame(rows)


# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

def run_uc3(db_name: str = DB_NAME, cohort: str = COHORT,
            seed_hp_key: str = SEED_HP_ENTITY_ID):
    """
    Execute the full UC3 analysis pipeline.
    """
    print("=" * 70)
    print("  UC3: Phenotype-Anchored Trans-Omic Subgraph Extraction")
    print("=" * 70)

    # --- Connect to ArangoDB ---
    db = setup_arangodb_connection(db_name)
    if db is None:
        raise RuntimeError("Failed to connect to ArangoDB")

    # --- Step 1: Find phenotype seed ---
    print("\n[1/8] Finding phenotype seed node...")
    seed_node = find_phenotype_seed(db, seed_hp_key)
    if seed_node is None:
        print("  ERROR: Could not find phenotype seed node.")
        print("  Trying broader search for HP nodes related to breast...")
        aql = """
        FOR n IN nodes
            FILTER n.class_code == "HP"
            FILTER CONTAINS(LOWER(n.label), "breast")
                OR CONTAINS(LOWER(n.label), "neoplasm")
            LIMIT 10
            RETURN {key: n._key, label: n.label, entity_id: n.entity_id}
        """
        candidates = list(db.aql.execute(aql))
        if candidates:
            print("  Available HP nodes:")
            for c in candidates:
                print(f"    {c['key']}: {c['label']} ({c.get('entity_id', '')})")
            # Use the first candidate
            seed_node = get_node_by_key(db, candidates[0]["key"])
        else:
            raise RuntimeError("No suitable HP phenotype node found in KG")

    print(f"  Seed: {seed_node['_key']} - {seed_node.get('label', 'N/A')}")
    print(f"         class_code: {seed_node.get('class_code')}, "
          f"entity_id: {seed_node.get('entity_id')}")

    # --- Step 2: KG Traversal (two-phase: hop-1 + hop-2 via diseases) ---
    print(f"\n[2/8] Traversing KG from phenotype (hop-1 + hop-2)...")
    traversal = traverse_phenotype_to_genes(
        db, seed_node["_key"],
        max_hop1=MAX_GENES_HOP1,
        max_hop2=MAX_GENES_HOP2
    )
    gene_set = traversal["genes"]

    if len(gene_set) < MIN_GENES_IN_SET:
        print(f"  WARNING: Only {len(gene_set)} genes found. "
              "Results may be limited.")

    # Save gene set
    gene_set_path = os.path.join(OUTPUT_DIR, "uc3_gene_set.json")
    with open(gene_set_path, "w", encoding="utf-8") as f:
        json.dump(gene_set, f, indent=2, default=str)
    print(f"  Saved gene set to {gene_set_path}")

    # --- Step 3: Get tumor samples ---
    print("\n[3/8] Finding BRCA tumor samples with multi-omic data...")
    all_sample_ids = list_samples_with_complete_omics(
        db, cohort=cohort,
        omic_types=["gene_expression", "protein"],
        specimen_type=None,
    )
    print(f"  Found {len(all_sample_ids)} samples with expression + protein")

    # Filter to Primary Tumor
    tumor_aql = """
    FOR s IN SAMPLES
        FILTER s.sample_type == "Primary Tumor"
        RETURN s._key
    """
    tumor_keys = set(db.aql.execute(tumor_aql))
    sample_ids = [s for s in all_sample_ids if s in tumor_keys]
    print(f"  Of which {len(sample_ids)} are Primary Tumor samples")

    if len(sample_ids) < MIN_SAMPLES:
        print(f"  WARNING: Only {len(sample_ids)} samples (minimum {MIN_SAMPLES})")

    # --- Step 4: Extract multi-omic profiles ---
    print(f"\n[4/8] Extracting multi-omic profiles for {len(gene_set)} genes...")
    omic_dfs = extract_multiomic_profiles(db, cohort, gene_set, sample_ids)

    # --- Step 5: Per-layer random baselines ---
    # Each omic layer has different gene coverage:
    #   expression / CNV : genome-wide (~20k+ genes)
    #   methylation      : needs HGNC symbol lookup (not all genes have array probes)
    #   protein (RPPA)   : ~300 genes total — sampling 200 from expression pool yields ~3 hits
    #
    # Fix 1 (methylation): build entrez→HGNC mapping so random genes carry proper labels;
    #   extract_multiomic_profiles skips methylation when label == entrez_id.
    # Fix 2 (protein): sample random genes exclusively from the protein-covered pool.
    print(f"\n[5/8] Building per-layer random gene set baselines...")
    gene_vector = get_gene_vector(db, cohort, id_type="entrez_id",
                                   omic_type="gene_expression")
    valid_genes = [g for g in gene_vector if g and g != ""]
    print(f"  Total genes in expression index: {len(valid_genes)}")

    # Entrez → HGNC symbol mapping (both vectors are position-sorted, so zip is safe)
    hgnc_vector = get_gene_vector(db, cohort, id_type="hgnc_symbol",
                                   omic_type="gene_expression")
    entrez_to_hgnc = {e: h for e, h in zip(valid_genes, hgnc_vector)
                      if e and h and e != h}

    # --- Expression / CNV / Methylation random set ---
    # Sampled from the expression pool; HGNC labels enable methylation lookup.
    rng = random.Random(SEED)
    expr_pool = [g for g in valid_genes if g not in gene_set]
    random_gene_ids = rng.sample(expr_pool, min(N_RANDOM_GENES, len(expr_pool)))
    random_gene_set = {eid: {"label": entrez_to_hgnc.get(eid, eid)}
                       for eid in random_gene_ids}

    print(f"  Extracting expr/CNV/methylation for {len(random_gene_ids)} random genes...")
    random_omic_dfs = extract_multiomic_profiles(db, cohort, random_gene_set, sample_ids)

    # --- Protein random set ---
    # Sample from the protein index (RPPA coverage ~300 genes).
    # N_RANDOM_PROTEIN is capped well below total coverage so sampling is valid.
    print(f"  Building protein random baseline from protein index (n≤{N_RANDOM_PROTEIN})...")
    protein_gene_vector = get_gene_vector(db, cohort, id_type="entrez_id",
                                           omic_type="protein")
    valid_protein_genes = [g for g in protein_gene_vector if g and g != ""]
    print(f"  Total genes in protein index: {len(valid_protein_genes)}")
    protein_pool = [g for g in valid_protein_genes if g not in gene_set]
    n_random_prot = min(N_RANDOM_PROTEIN, len(protein_pool))
    random_protein_ids = rng.sample(protein_pool, n_random_prot)
    random_protein_set = {eid: {"label": entrez_to_hgnc.get(eid, eid)}
                          for eid in random_protein_ids}

    print(f"  Extracting protein data for {n_random_prot} random protein-covered genes...")
    random_prot_omic = extract_multiomic_profiles(db, cohort, random_protein_set, sample_ids)
    random_omic_dfs["protein"] = random_prot_omic.get("protein")

    # --- Methylation random set ---
    # Methylation array has partial coverage (~37% of expression genes).
    # Sampling from expression pool → only ~74/200 random genes have methylation data.
    # Fix: sample directly from the methylation-covered gene pool (HGNC symbols with
    # at least one probe in the array) to get a random baseline size-matched to pheno.
    print(f"  Building methylation random baseline from methylation index (n≤{N_RANDOM_GENES})...")
    meth_config = OMIC_CONFIG["methylation"]
    meth_index_key = f"{meth_config['index_key_prefix']}_{cohort}"
    meth_aql = f"""
    FOR doc IN {meth_config['collection']}
        FILTER doc.data_type == @indexDataType
        FILTER doc._key == @indexKey
        RETURN UNIQUE(FLATTEN(FOR pm IN doc.probe_mappings RETURN pm.gene_symbols))
    """
    meth_result = list(db.aql.execute(meth_aql, bind_vars={
        "indexKey": meth_index_key,
        "indexDataType": meth_config["index_data_type"]}))
    meth_symbols = meth_result[0] if meth_result and meth_result[0] else []
    print(f"  Total gene symbols with methylation coverage: {len(meth_symbols)}")

    # Reverse HGNC → Entrez mapping (restrict to genes also in expression index)
    hgnc_to_entrez = {h: e for e, h in entrez_to_hgnc.items()}
    meth_pool = [hgnc_to_entrez[s] for s in meth_symbols
                 if s in hgnc_to_entrez and hgnc_to_entrez[s] not in gene_set]
    n_random_meth = min(N_RANDOM_GENES, len(meth_pool))
    print(f"  Methylation-covered pool after excluding phenotype: {len(meth_pool)}")
    random_meth_ids = rng.sample(meth_pool, n_random_meth)

    # Inline methylation extraction (avoids re-querying expr/CNV/protein)
    print(f"  Extracting methylation data for {n_random_meth} methylation-covered random genes...")
    meth_records_rand = {}
    for eid in random_meth_ids:
        symbol = entrez_to_hgnc.get(eid)
        if not symbol:
            continue
        data = get_methylation_by_gene(db, cohort, gene_symbol=symbol,
                                        sample_ids=sample_ids)
        if data and data.get("probes") and data.get("samples"):
            for s in data["samples"]:
                sid = s["sample_id"]
                vals = [v for v in s["values"] if v is not None]
                if vals:
                    meth_records_rand.setdefault(eid, {})[sid] = np.mean(vals)
    if meth_records_rand:
        rand_meth_df = pd.DataFrame(meth_records_rand).apply(pd.to_numeric, errors="coerce")
        random_omic_dfs["methylation"] = rand_meth_df
        print(f"      {rand_meth_df.shape[0]} samples x {rand_meth_df.shape[1]} genes")

    # --- Step 6: Summary statistics + statistical tests ---
    print("\n[6/8] Computing summary statistics and statistical tests...")
    summary_df = compute_summary_stats(omic_dfs, gene_set)
    summary_path = os.path.join(OUTPUT_DIR, "uc3_summary_stats.csv")
    summary_df.to_csv(summary_path, index=False)
    print(summary_df.to_string(index=False))

    # Statistical comparison: phenotype gene set vs random
    print("\n  --- Phenotype Gene Set vs Random: Statistical Comparison ---")
    comparison_rows = []
    for layer_name in ["expression", "cnv", "methylation", "protein"]:
        pheno_df = omic_dfs.get(layer_name)
        rand_df = random_omic_dfs.get(layer_name)
        if pheno_df is None or rand_df is None:
            continue
        if pheno_df.empty or rand_df.empty:
            continue

        # Compare mean absolute signal per gene across samples
        pheno_means = pheno_df.mean(axis=0).dropna()
        rand_means = rand_df.mean(axis=0).dropna()

        if len(pheno_means) < 5 or len(rand_means) < 5:
            continue

        u_stat, p_val = stats.mannwhitneyu(
            pheno_means.abs(), rand_means.abs(), alternative="two-sided"
        )
        effect_r = 2 * u_stat / (len(pheno_means) * len(rand_means)) - 1

        # Cross-gene variance comparison
        pheno_var = pheno_df.var(axis=0).mean()
        rand_var = rand_df.var(axis=0).mean()

        row = {
            "layer": layer_name,
            "n_pheno_genes": len(pheno_means),
            "n_random_genes": len(rand_means),
            "pheno_mean_abs_signal": pheno_means.abs().mean(),
            "random_mean_abs_signal": rand_means.abs().mean(),
            "pheno_cross_gene_var": pheno_var,
            "random_cross_gene_var": rand_var,
            "mannwhitney_U": u_stat,
            "p_value": p_val,
            "effect_size_r": effect_r,
        }
        comparison_rows.append(row)

        sig = "***" if p_val < 0.001 else "**" if p_val < 0.01 else "*" if p_val < 0.05 else "ns"
        print(f"  {layer_name:15s}  pheno={pheno_means.abs().mean():.4f}  "
              f"random={rand_means.abs().mean():.4f}  "
              f"U={u_stat:.0f}  p={p_val:.2e}  r={effect_r:.3f}  {sig}")

    comparison_df = pd.DataFrame(comparison_rows)
    if not comparison_df.empty:
        comp_path = os.path.join(OUTPUT_DIR, "uc3_phenotype_vs_random.csv")
        comparison_df.to_csv(comp_path, index=False)
        print(f"  Saved comparison: {comp_path}")

    # --- Step 7: Visualizations ---
    print("\n[7/8] Generating visualizations...")

    # Mean expression for coloring
    expr_df = omic_dfs.get("expression")
    mean_expression = None
    if expr_df is not None and not expr_df.empty:
        mean_expression = expr_df.mean(axis=0)

    # 7a. Semantic subgraph
    print("  [7a] Semantic subgraph...")
    plot_semantic_subgraph(
        traversal, seed_node,
        mean_expression=mean_expression,
        title=f"UC3: Semantic Subgraph from {seed_node.get('label', 'Phenotype')[:40]}",
        output_path=os.path.join(OUTPUT_DIR, "uc3_semantic_subgraph.png"),
    )

    # 7b. Multi-omic heatmap
    print("  [7b] Multi-omic heatmap...")
    plot_multiomic_heatmap(
        omic_dfs, gene_set,
        title=f"UC3: Multi-Omic Heatmap — {seed_node.get('label', 'Phenotype')[:40]} Gene Set",
        output_path=os.path.join(OUTPUT_DIR, "uc3_multiomic_heatmap.png"),
    )

    # 7c. Omic barplot comparison (phenotype vs single random set, per-layer)
    print("  [7c] Omic barplot comparison...")
    plot_omic_barplot_comparison(
        omic_dfs, random_omic_dfs, gene_set, random_gene_set,
        title=f"UC3: Multi-Omic Signal — {seed_node.get('label', 'Phenotype')[:40]}",
        output_path=os.path.join(OUTPUT_DIR, "uc3_barplot_comparison.png"),
    )

    # 7c2. Bootstrap profile (phenotype vs null from per-layer random_omic_dfs)
    print("  [7c2] Bootstrap null profile...")
    plot_bootstrap_profile(
        omic_dfs, valid_genes, gene_set,
        random_omic_dfs=random_omic_dfs,
        title=f"UC3: Multi-Omic Signal vs Bootstrap Null — "
              f"{seed_node.get('label', 'Phenotype')[:40]}",
        output_path=os.path.join(OUTPUT_DIR, "uc3_bootstrap_profile.png"),
    )

    # 7d. Sankey diagram
    print("  [7d] Sankey diagram...")
    plot_sankey_diagram(
        traversal, omic_dfs, seed_node, gene_set,
        title=f"UC3: {seed_node.get('label', 'Phenotype')[:30]} -> Genes -> Omic Layers",
        output_path=os.path.join(OUTPUT_DIR, "uc3_sankey_diagram.png"),
    )

    # 7e. Publication-ready composite figures (Bioinformatics style)
    # Both composite figures now use the bootstrap null panel (C/B) instead of radar.
    print("  [7e] Composite figure (2-panel, Bioinformatics style)...")
    plot_composite_figure_2panel(
        omic_dfs, gene_set,
        valid_genes=valid_genes,
        random_omic_dfs=random_omic_dfs,
        seed_node=seed_node,
        output_path=os.path.join(OUTPUT_DIR, "uc3_figure_composite_2panel.png"),
    )
    print("  [7f] Composite figure (4-panel, Bioinformatics style)...")
    plot_composite_figure_4panel(
        traversal, omic_dfs, gene_set, seed_node,
        valid_genes=valid_genes,
        random_omic_dfs=random_omic_dfs,
        mean_expression=mean_expression,
        output_path=os.path.join(OUTPUT_DIR, "uc3_figure_composite_4panel.png"),
    )

    # 7g/7h. Final publication figure — 3-panel stacked (violin + boxplot variants)
    print("  [7g] Composite figure (3-panel stacked, violin null)...")
    plot_composite_figure_3panel_stacked(
        traversal, omic_dfs, gene_set, seed_node,
        valid_genes=valid_genes,
        random_omic_dfs=random_omic_dfs,
        mean_expression=mean_expression,
        output_path=os.path.join(OUTPUT_DIR, "uc3_figure_3panel_stacked.png"),
    )
    print("  [7h] Composite figure (3-panel stacked, boxplot null)...")
    plot_composite_figure_3panel_stacked_box(
        traversal, omic_dfs, gene_set, seed_node,
        valid_genes=valid_genes,
        random_omic_dfs=random_omic_dfs,
        mean_expression=mean_expression,
        output_path=os.path.join(OUTPUT_DIR, "uc3_figure_3panel_stacked_box.png"),
    )

    # --- Step 8: Save full results ---
    print("\n[8/8] Saving results...")

    # Save omic matrices
    for layer_name, df in omic_dfs.items():
        if df is not None and not df.empty:
            path = os.path.join(OUTPUT_DIR, f"uc3_{layer_name}_matrix.csv")
            df.to_csv(path)
            print(f"  Saved {layer_name} matrix: {path}")

    # Save traversal metadata
    n_hop1 = len([g for g in gene_set.values() if g.get("hop_distance") == 1])
    n_hop2 = len([g for g in gene_set.values() if g.get("hop_distance") == 2])
    n_breast_d = len([d for d in traversal["diseases"].values()
                      if d.get("is_breast_specific")])
    traversal_meta = {
        "seed_node": {
            "key": seed_node["_key"],
            "label": seed_node.get("label"),
            "class_code": seed_node.get("class_code"),
            "entity_id": seed_node.get("entity_id"),
        },
        "gene_set_size": len(gene_set),
        "n_hop1_genes": n_hop1,
        "n_hop2_genes": n_hop2,
        "n_proteins": len(traversal["proteins"]),
        "n_diseases": len(traversal["diseases"]),
        "n_breast_specific_diseases": n_breast_d,
        "n_pathways": len(traversal["pathways"]),
        "n_samples": len(sample_ids),
        "omic_coverage": {
            layer: {"n_genes": df.shape[1] if df is not None else 0,
                    "n_samples": df.shape[0] if df is not None else 0}
            for layer, df in omic_dfs.items()
        },
    }
    meta_path = os.path.join(OUTPUT_DIR, "uc3_traversal_metadata.json")
    with open(meta_path, "w", encoding="utf-8") as f:
        json.dump(traversal_meta, f, indent=2, default=str)
    print(f"  Saved metadata: {meta_path}")

    # Save full traversal structure (diseases, proteins, pathways, subgraph_edges)
    # These are needed to reproduce the semantic subgraph and Sankey plots exactly.
    traversal_full = {
        "diseases": traversal["diseases"],
        "proteins": traversal["proteins"],
        "pathways": traversal["pathways"],
        "subgraph_edges": traversal["subgraph_edges"],
    }
    traversal_full_path = os.path.join(OUTPUT_DIR, "uc3_traversal_full.json")
    with open(traversal_full_path, "w", encoding="utf-8") as f:
        json.dump(traversal_full, f, indent=2, default=str)
    print(f"  Saved full traversal: {traversal_full_path}")

    print("\n" + "=" * 70)
    print("  UC3 analysis complete!")
    print(f"  Results saved to: {OUTPUT_DIR}")
    print(f"  Gene set: {len(gene_set)} genes ({n_hop1} hop-1 + {n_hop2} hop-2) "
          f"from '{seed_node.get('label', 'N/A')}'")
    print("=" * 70)

    return traversal, omic_dfs, summary_df


def run_uc3_plots_only():
    """
    Skip the analysis and regenerate all plots from saved files in OUTPUT_DIR.

    Requires:
    - uc3_gene_set.json
    - uc3_traversal_metadata.json
    - uc3_traversal_full.json  (diseases, proteins, pathways, subgraph_edges)
    - uc3_expression_matrix.csv, uc3_cnv_matrix.csv,
      uc3_methylation_matrix.csv, uc3_protein_matrix.csv
    - uc3_phenotype_vs_random.csv (for radar random baseline stats)
    """
    print("=" * 70)
    print("  UC3: Regenerating plots from saved results (skip-analysis mode)")
    print("=" * 70)

    gene_set_path = os.path.join(OUTPUT_DIR, "uc3_gene_set.json")
    meta_path = os.path.join(OUTPUT_DIR, "uc3_traversal_metadata.json")
    traversal_full_path = os.path.join(OUTPUT_DIR, "uc3_traversal_full.json")

    missing = [p for p in [gene_set_path] if not os.path.exists(p)]
    if missing:
        print(f"  ERROR: Missing required files:\n  " + "\n  ".join(missing))
        print("  Run without --skip-analysis first to generate the data.")
        return

    # --- Load gene set ---
    print(f"  Loading: {gene_set_path}")
    with open(gene_set_path, encoding="utf-8") as f:
        gene_set = json.load(f)
    print(f"  Gene set: {len(gene_set)} genes")

    # --- Load traversal metadata (seed node info) ---
    traversal_meta = {}
    if os.path.exists(meta_path):
        print(f"  Loading: {meta_path}")
        with open(meta_path, encoding="utf-8") as f:
            traversal_meta = json.load(f)

    seed_info = traversal_meta.get("seed_node", {})
    seed_node = {
        "_key": seed_info.get("key", "HP_seed"),
        "label": seed_info.get("label", "Breast Carcinoma"),
        "class_code": seed_info.get("class_code", "HP"),
        "entity_id": seed_info.get("entity_id", ""),
    }
    print(f"  Seed: {seed_node['_key']} - {seed_node['label']}")

    # --- Load full traversal (diseases, proteins, pathways, subgraph_edges) ---
    if os.path.exists(traversal_full_path):
        print(f"  Loading: {traversal_full_path}")
        with open(traversal_full_path, encoding="utf-8") as f:
            traversal_full = json.load(f)
        # subgraph_edges is stored as list of lists, convert back to list of tuples
        subgraph_edges = [tuple(e) for e in traversal_full.get("subgraph_edges", [])]
        traversal = {
            "genes": gene_set,
            "proteins": traversal_full.get("proteins", {}),
            "diseases": traversal_full.get("diseases", {}),
            "pathways": traversal_full.get("pathways", {}),
            "subgraph_edges": subgraph_edges,
        }
        print(f"  Diseases: {len(traversal['diseases'])}, "
              f"Proteins: {len(traversal['proteins'])}, "
              f"Pathways: {len(traversal['pathways'])}, "
              f"Edges: {len(subgraph_edges)}")
    else:
        # Fallback: reconstruct minimally from gene_set (no disease/pathway nodes)
        print(f"  WARNING: {traversal_full_path} not found — "
              "subgraph/Sankey will show genes only (no disease/pathway nodes).")
        seed_key = seed_node["_key"]
        subgraph_edges = []
        for eid, info in gene_set.items():
            gk = info.get("kg_key", eid)
            hop = info.get("hop_distance", 1)
            if hop == 1:
                subgraph_edges.append((seed_key, gk,
                                       info.get("path_predicates", ["linked"])[-1]))
            elif hop == 2:
                via_key = info.get("via_disease_key", seed_key)
                subgraph_edges.append((via_key, gk,
                                       info.get("path_predicates", ["linked"])[-1]))
        traversal = {
            "genes": gene_set,
            "proteins": {},
            "diseases": {},
            "pathways": {},
            "subgraph_edges": subgraph_edges,
        }

    # --- Load omic matrices ---
    omic_dfs = {}
    for layer_name in ["expression", "cnv", "methylation", "protein"]:
        fpath = os.path.join(OUTPUT_DIR, f"uc3_{layer_name}_matrix.csv")
        if os.path.exists(fpath):
            print(f"  Loading: {fpath}")
            df = pd.read_csv(fpath, index_col=0)
            df.columns = df.columns.astype(str)
            omic_dfs[layer_name] = df
        else:
            print(f"  WARNING: {fpath} not found, {layer_name} layer skipped.")
            omic_dfs[layer_name] = None

    print(f"\n  Generating visualizations...")

    # 7a. Semantic subgraph
    print("  [7a] Semantic subgraph...")
    expr_df = omic_dfs.get("expression")
    mean_expression = expr_df.mean(axis=0) if expr_df is not None and not expr_df.empty else None
    try:
        plot_semantic_subgraph(
            traversal, seed_node,
            mean_expression=mean_expression,
            title=f"UC3: Semantic Subgraph from {seed_node.get('label', 'Phenotype')[:40]}",
            output_path=os.path.join(OUTPUT_DIR, "uc3_semantic_subgraph.png"),
        )
    except Exception as exc:
        print(f"  ERROR in semantic subgraph: {exc}")

    # 7b. Multi-omic heatmap
    print("  [7b] Multi-omic heatmap...")
    try:
        plot_multiomic_heatmap(
            omic_dfs, gene_set,
            title=f"UC3: Multi-Omic Heatmap — {seed_node.get('label', 'Phenotype')[:40]} Gene Set",
            output_path=os.path.join(OUTPUT_DIR, "uc3_multiomic_heatmap.png"),
        )
    except Exception as exc:
        print(f"  ERROR in multiomic heatmap: {exc}")

    # 7c. Omic barplot comparison (needs a random set — reconstruct from CSV)
    print("  [7c] Omic barplot comparison...")
    comp_path = os.path.join(OUTPUT_DIR, "uc3_phenotype_vs_random.csv")
    random_omic_dfs = {}; random_gene_set_plot = {}
    if os.path.exists(comp_path):
        comp_df = pd.read_csv(comp_path)
        for _, row in comp_df.iterrows():
            layer = row["layer"]
            n_genes = int(row.get("n_random_genes", 1))
            mean_abs = float(row.get("random_mean_abs_signal", 0))
            rand_cols = [f"rand_{i}" for i in range(n_genes)]
            random_omic_dfs[layer] = pd.DataFrame(
                np.full((1, n_genes), mean_abs), columns=rand_cols)
            for c in rand_cols:
                random_gene_set_plot[c] = {"label": c}
    try:
        plot_omic_barplot_comparison(
            omic_dfs, random_omic_dfs, gene_set, random_gene_set_plot,
            title=f"UC3: Multi-Omic Signal — {seed_node.get('label', 'Phenotype')[:40]}",
            output_path=os.path.join(OUTPUT_DIR, "uc3_barplot_comparison.png"),
        )
    except Exception as exc:
        print(f"  ERROR in barplot comparison: {exc}")

    # 7c2. Bootstrap profile (uses columns in existing omic DFs — no DB queries)
    print("  [7c2] Bootstrap null profile...")
    valid_genes_local = (list(omic_dfs["expression"].columns)
                         if omic_dfs.get("expression") is not None else [])
    try:
        plot_bootstrap_profile(
            omic_dfs, valid_genes_local, gene_set,
            random_omic_dfs=random_omic_dfs,
            title=f"UC3: Multi-Omic Signal vs Bootstrap Null — "
                  f"{seed_node.get('label', 'Phenotype')[:40]}",
            output_path=os.path.join(OUTPUT_DIR, "uc3_bootstrap_profile.png"),
        )
    except Exception as exc:
        print(f"  ERROR in bootstrap profile: {exc}")

    # 7d. Sankey diagram
    print("  [7d] Sankey diagram...")
    try:
        plot_sankey_diagram(
            traversal, omic_dfs, seed_node, gene_set,
            title=f"UC3: {seed_node.get('label', 'Phenotype')[:30]} -> Genes -> Omic Layers",
            output_path=os.path.join(OUTPUT_DIR, "uc3_sankey_diagram.png"),
        )
    except Exception as exc:
        print(f"  ERROR in Sankey diagram: {exc}")

    # 7e/7f. Composite figures — bootstrap panel instead of radar
    print("  [7e] Composite figure (2-panel, Bioinformatics style)...")
    try:
        plot_composite_figure_2panel(
            omic_dfs, gene_set,
            valid_genes=valid_genes_local,
            random_omic_dfs=random_omic_dfs,
            seed_node=seed_node,
            output_path=os.path.join(OUTPUT_DIR, "uc3_figure_composite_2panel.png"),
        )
    except Exception as exc:
        print(f"  ERROR in 2-panel composite: {exc}")

    print("  [7f] Composite figure (4-panel, Bioinformatics style)...")
    try:
        plot_composite_figure_4panel(
            traversal, omic_dfs, gene_set, seed_node,
            valid_genes=valid_genes_local,
            random_omic_dfs=random_omic_dfs,
            mean_expression=mean_expression,
            output_path=os.path.join(OUTPUT_DIR, "uc3_figure_composite_4panel.png"),
        )
    except Exception as exc:
        print(f"  ERROR in 4-panel composite: {exc}")

    print("  [7g] Composite figure (3-panel stacked, violin null)...")
    try:
        plot_composite_figure_3panel_stacked(
            traversal, omic_dfs, gene_set, seed_node,
            valid_genes=valid_genes_local,
            random_omic_dfs=random_omic_dfs,
            mean_expression=mean_expression,
            output_path=os.path.join(OUTPUT_DIR, "uc3_figure_3panel_stacked.png"),
        )
    except Exception as exc:
        print(f"  ERROR in 3-panel stacked (violin): {exc}")

    print("  [7h] Composite figure (3-panel stacked, boxplot null)...")
    try:
        plot_composite_figure_3panel_stacked_box(
            traversal, omic_dfs, gene_set, seed_node,
            valid_genes=valid_genes_local,
            random_omic_dfs=random_omic_dfs,
            mean_expression=mean_expression,
            output_path=os.path.join(OUTPUT_DIR, "uc3_figure_3panel_stacked_box.png"),
        )
    except Exception as exc:
        print(f"  ERROR in 3-panel stacked (box): {exc}")

    print("\n" + "=" * 70)
    print("  UC3 plots regenerated!")
    print(f"  Output: {OUTPUT_DIR}")
    print("=" * 70)


#%% Entry point
if __name__ == "__main__":
    skip_analysis = "--skip-analysis" in sys.argv

    if "ipykernel" not in sys.modules:
        if skip_analysis:
            run_uc3_plots_only()
        else:
            run_uc3()
    else:
        # Interactive: set skip_analysis = True to regenerate plots only
        if skip_analysis:
            run_uc3_plots_only()
        else:
            traversal, omic_dfs, summary_df = run_uc3()
