#%%
"""
Command-line helper to build a transomic property-graph JSON for one sample.

Usage example:
    python build_transomic_network.py \
        --db-name PKT_test10000 \
        --cohort TCGA-BRCA \
        --sample-id TCGA-XX-YYYY \
        --specimen-type "Primary Tumor" \
        --value-abs-threshold 1.0 \
        --top-n 500 \
        --predicate-filter "has gene product" "interacts with"

If --sample-id is omitted, the first sample having all requested omic types
(and optional specimen-type) is selected automatically.
"""

"""
Uso
Notebook (rete densa per topologia):
graph_obj, saved_path = build_transomic_network_for_sample(
    db_name=DB_NAME, cohort=COHORT, sample_id=SAMPLE_ID,
    value_threshold_per_omic=DEFAULT_VALUE_THRESHOLDS,
    predicate_filter=BIOLOGICAL_PREDICATES,
    edge_limit=200000,
)

CLI:
python build_transomic_network.py --cohort TCGA-BRCA --sample-id ... \
    --use-default-thresholds --use-biological-predicates --edge-limit 200000

Override per-omica:
--value-threshold-per-omic gene_expression=2.0 protein=0.8


Per la versione snella ML/visualizzazione basta ridurre edge_limit (es. 10000) o alzare le soglie.
 Il value_abs_threshold globale resta come fallback per omiche non specificate.
"""


"""
DEV

% rete completa, zoom su caso, insight 
% visualizzazione di path multiomici non facilmente  riconoscibili nei db scollegato
% visualizzazione del pezzo, poi arrichimato semantico
% statistiche 
[complessità della rete]
- generale
- specifica
- arricchita semanticamte
- dati statistici grazie all'integrazione semantica

- dicussione : cosa si può fare con il sistema, scalabile, aperto etc

% template della rivista --- limiti di sottomissime



"""


#%% Imports
import argparse
import json
import os
import sys
from typing import Dict, List, Optional, Tuple

from arangodb_utils import setup_arangodb_connection
from query_utils import (
    OMIC_CONFIG,
    build_transomic_property_graph,
    list_samples_with_complete_omics,
)

# =============================================================================
# BIOLOGICAL PRESETS
# =============================================================================
# Per-omic |value| thresholds. Scales differ across omics so a single global
# threshold is not appropriate. Values here are defensible starting points for
# bulk TCGA-like data; tune per study.
#
# - gene_expression (TPM): drop near-zero expression to remove noise/off-genes.
# - cnv (relative copy number log2 ratio): keep >|0.2| to flag gain/loss.
# - methylation (beta in [0,1]): keep probes far from 0.5 (|beta-0.5| proxy not
#   used; here we simply keep |beta|>=0.2, meaning not fully unmethylated. For
#   differential methylation use the per-study preferred scale).
# - mirna (expression, log or RPM depending on dataset): drop near-zero.
# - protein (RPPA normalized log2 ratio ~[-3,3]): keep |value|>=0.5.
DEFAULT_VALUE_THRESHOLDS: Dict[str, float] = {
    "gene_expression": 1.0,
    "cnv": 0.2,
    "methylation": 0.2,
    "mirna": 1.0,
    "protein": 0.5,
}

# Biologically-meaningful KG predicates for topological analysis.
# Focused on molecular interactions, regulation, gene-product, transcription,
# and disease/phenotype links. Drops pure ontological/structural edges
# (type, part_of, has_quality, etc.) that inflate the graph without carrying
# sample-specific biology.
BIOLOGICAL_PREDICATES: List[str] = [
    # Interactions
    "molecularly interacts with",
    "interacts with",
    "genetically interacts with",
    # Regulation
    "regulates (processual)",
    "positively regulates",
    "negatively regulates",
    "causally influences",
    "causally influenced by",
    # Transcription / translation
    "transcribed to",
    "transcribed from",
    "ribosomal translation of",
    "ribosomally translates to",
    # Gene product
    "has gene product",
    "gene product of",
    "has_gene_template",
    # Participation / function
    "participates in",
    "has participant",
    "has function",
    "function of",
    # Disease / phenotype
    "has phenotype",
    "phenotype of",
    "causes or contributes to condition",
    "is substance that treats",
    "is treated by substance",
    "disease has basis in dysfunction of",
]

# ---------------------------------------------------------------------------
# MOLECULAR_INTERACTION_PREDICATES
# Lightweight set for network topology / graph-ML tasks.
# Focus: direct molecular wiring between genes, transcripts, proteins.
# Rationale: these are the highest-confidence, most reciprocal edges in the
# PKT (interaction pairs are symmetric) and are the backbone of PPI / GRN
# analyses. Excludes disease/phenotype context to keep the graph tight.
# Expected edge volume: medium (~interactions 1.5M + gene-product 40K +
# transcription 450K → heavily reduced after node filtering).
# Good for: centrality, community detection, GNN on molecular circuits.
# ---------------------------------------------------------------------------
MOLECULAR_INTERACTION_PREDICATES: List[str] = [
    # Direct molecular interactions — highest volume, most informative
    "molecularly interacts with",       # 1384824 edges
    "interacts with",                   # 176638
    "genetically interacts with",       # 3386
    # Transcription / translation — connects genes to transcripts to proteins
    "transcribed to",                   # 182692
    "transcribed from",                 # 182692
    "ribosomal translation of",         # 44205
    "ribosomally translates to",        # 44205
    # Gene → protein mapping
    "has gene product",                 # 19521
    "gene product of",                  # 19521
    "has_gene_template",                # 19845
    # Causal regulation (directional signal flow)
    "causally influences",              # 145129
    "causally influenced by",           # 145129
    "positively regulates",             # 3115
    "negatively regulates",             # 3126
]

# ---------------------------------------------------------------------------
# DISEASE_MECHANISM_PREDICATES
# Lightweight set for clinical / translational use cases.
# Focus: gene–disease, gene–phenotype, drug–target links only.
# Rationale: minimal set that connects omic alterations directly to clinical
# entities. Strips all molecular wiring (interactions, transcription) to
# highlight only the disease-relevant semantic layer.
# Expected edge volume: low (~disease causation 83K + phenotype 856K +
# treatment 558K → compact after node filtering).
# Good for: patient stratification, drug repurposing, pathway-to-disease
# mapping, explainability in clinical ML models.
# ---------------------------------------------------------------------------
DISEASE_MECHANISM_PREDICATES: List[str] = [
    # Gene/protein → disease causation
    "causes or contributes to condition",   # 83530
    "disease has basis in dysfunction of",  # 4594
    # Gene/protein → phenotype
    "has phenotype",                        # 428374
    "phenotype of",                         # 428374
    # Drug → disease treatment
    "is substance that treats",             # 279052
    "is treated by substance",              # 279052
    # Disease mechanism detail
    "disease has location",                 # 2757
    "disease has feature",                  # 1105
]


"""
## I tre PREDICATE_FILTERS a confronto

| **Categoria**                  | **BIOLOGICAL_PREDICATES**                  | **MOLECULAR_INTERACTION_PREDICATES** | **DISEASE_MECHANISM_PREDICATES** |
|-------------------------------|--------------------------------------------|-------------------------------------|----------------------------------|
| **Predicati**                 | ~25                                        | 15                                  | 8                                |
| **Focus**                     | Biologia molecolare + clinica              | Solo cablaggio molecolare           | Solo gene→malattia→farmaco       |
| **Edge stimati (PKT)**        | medio-alti                                 | medi                                | bassi                            |
| **Usa per**                   | analisi topologica completa                | PPI/GRN, centrality, GNN            | stratificazione pazienti, drug repurposing |
| **CLI flag**                  | `--use-biological-predicates`              | `--use-molecular-predicates`        | `--use-disease-predicates`       |

### Descrizioni dettagliate

**MOLECULAR_INTERACTION_PREDICATES**  
Taglia tutto il layer clinico e tiene solo il cablaggio molecolare diretto: interazioni, trascrizione, traduzione, regolazione causale.  
*Ideale per*: analisi di rete molecolare (betweenness centrality, community detection, GNN su pathway).

**DISEASE_MECHANISM_PREDICATES**  
Taglia tutto il layer molecolare e tiene solo i link gene→fenotipo→malattia→farmaco.  
*Il grafo più snello*: connette le alterazioni omiche del campione direttamente alle entità cliniche rilevanti, senza rumore di interazioni molecolari.
"""

# from arangodb_utils import *
# # --- FUNZIONI DI BASE DI ARANGODB ---
# db_connection = setup_arangodb_connection("PKT_test10000")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build a transomic property graph for a sample")
    parser.add_argument("--db-name", default="PKT_test10000", help="ArangoDB database name")
    parser.add_argument("--cohort", required=True, help="Cohort identifier (e.g., TCGA-BRCA)")
    parser.add_argument("--sample-id", help="Sample _key; if omitted auto-pick one with complete omics")
    parser.add_argument(
        "--specimen-type",
        help="Optional specimen_type filter when auto-selecting a sample (e.g., Primary Tumor, Solid Tissue Normal)",
    )
    parser.add_argument(
        "--omic-types",
        nargs="*",
        default=list(OMIC_CONFIG.keys()),
        help="Subset of omic types to include (default: all)",
    )
    parser.add_argument(
        "--value-abs-threshold",
        type=float,
        default=None,
        help="Global |value| threshold (fallback for omics not in --value-threshold-per-omic)",
    )
    parser.add_argument(
        "--value-threshold-per-omic",
        nargs="*",
        default=None,
        metavar="OMIC=THRESHOLD",
        help="Per-omic |value| thresholds, e.g. gene_expression=1.0 cnv=0.2",
    )
    parser.add_argument(
        "--use-default-thresholds",
        action="store_true",
        help="Apply DEFAULT_VALUE_THRESHOLDS preset for all omics",
    )
    parser.add_argument(
        "--top-n",
        type=int,
        default=None,
        help="Keep only the first N features per omic after filtering",
    )
    parser.add_argument(
        "--predicate-filter",
        nargs="*",
        default=None,
        help="Explicit list of predicate_label values to retain in KG edges (overrides presets)",
    )
    parser.add_argument(
        "--use-biological-predicates",
        action="store_true",
        default=True,
        help="Apply BIOLOGICAL_PREDICATES preset (default, ~25 predicates)",
    )
    parser.add_argument(
        "--use-molecular-predicates",
        action="store_true",
        help="Apply MOLECULAR_INTERACTION_PREDICATES preset instead of default",
    )
    parser.add_argument(
        "--use-disease-predicates",
        action="store_true",
        help="Apply DISEASE_MECHANISM_PREDICATES preset instead of default",
    )
    parser.add_argument(
        "--no-predicate-filter",
        action="store_true",
        help="Disable predicate filtering — include ALL KG edges (warning: very large graphs)",
    )
    parser.add_argument(
        "--edge-limit",
        type=int,
        default=None,
        help="Optional hard cap on KG edges (AQL LIMIT, no ranking). Default: no cap.",
    )
    parser.add_argument(
        "--star-expansion",
        action="store_true",
        help="Return edges touching at least one mapped node (default: induced subgraph — both endpoints mapped)",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output path for JSON; defaults to ../transomic-networks/transomic_graph_<cohort>_<sample>.json",
    )
    return parser.parse_args()



def pick_sample(db, cohort: str, omic_types: List[str], specimen_type: Optional[str]) -> Optional[str]:
    candidates = list_samples_with_complete_omics(db, cohort=cohort, omic_types=omic_types, specimen_type=specimen_type)
    return candidates[0] if candidates else None


#%% Notebook-friendly helper
def build_transomic_network_for_sample(
    db_name: str,
    cohort: str,
    sample_id: str,
    omic_types: Optional[List[str]] = None,
    value_abs_threshold: Optional[float] = None,
    value_threshold_per_omic: Optional[Dict[str, float]] = None,
    top_n: Optional[int] = None,
    predicate_filter: Optional[List[str]] = None,
    edge_limit: Optional[int] = None,
    induced_subgraph: bool = False,
    output_path: Optional[str] = None,
    tag_flag: Optional[str] = None,
) -> Tuple[dict, str]:
    """
    Build and persist a transomic property-graph for a single sample.

    Key filters:
    - value_threshold_per_omic: dict {omic: threshold} — preferred way to drop
      low-signal features. Scales differ across omics so prefer this over
      value_abs_threshold.
    - predicate_filter: list of predicate_label values to retain in KG edges.
      Default: BIOLOGICAL_PREDICATES. Pass None to include all KG edges.
    - edge_limit: hard cap (AQL LIMIT, no ranking). For dense analysis raise
      this well above the default 5000.
    """
    db = setup_arangodb_connection(db_name)
    if db is None:
        raise RuntimeError("Failed to connect to ArangoDB")

    omic_types = omic_types or list(OMIC_CONFIG.keys())

    graph = build_transomic_property_graph(
        db_connection=db,
        cohort=cohort,
        sample_id=sample_id,
        omic_types=omic_types,
        value_abs_threshold=value_abs_threshold,
        value_threshold_per_omic=value_threshold_per_omic,
        top_n=top_n,
        predicate_filter=predicate_filter,
        edge_limit=edge_limit,
        induced_subgraph=induced_subgraph,
        include_edges=True,
    )

    out_path = output_path
    if out_path is None:
        safe_sample = sample_id.replace("/", "-")
        flag = "_thresholded" if value_threshold_per_omic or value_abs_threshold else "_full"
        flag += f"_{tag_flag}" if tag_flag else ""

        out_path = os.path.join("..", "transomic-networks", f"transomic_graph_{cohort}_{safe_sample}{flag}.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2)

    return graph, out_path


def _parse_per_omic_thresholds(items: Optional[List[str]]) -> Optional[Dict[str, float]]:
    if not items:
        return None
    out: Dict[str, float] = {}
    for item in items:
        if "=" not in item:
            raise ValueError(f"Invalid --value-threshold-per-omic entry: {item!r} (expected OMIC=FLOAT)")
        k, v = item.split("=", 1)
        out[k.strip()] = float(v)
    return out


#%% CLI glue remains for standalone execution
def main():
    args = parse_args()

    db = setup_arangodb_connection(args.db_name)
    if db is None:
        raise SystemExit("Failed to connect to ArangoDB")

    sample_id = args.sample_id
    if not sample_id:
        sample_id = pick_sample(db, cohort=args.cohort, omic_types=args.omic_types, specimen_type=args.specimen_type)
        if not sample_id:
            raise SystemExit("No sample found with complete omics for requested settings")
        print(f"Auto-selected sample: {sample_id}")

    per_omic = _parse_per_omic_thresholds(args.value_threshold_per_omic)
    if args.use_default_thresholds:
        merged = dict(DEFAULT_VALUE_THRESHOLDS)
        if per_omic:
            merged.update(per_omic)
        per_omic = merged

    predicate_filter = args.predicate_filter
    if args.no_predicate_filter:
        predicate_filter = None
    elif args.use_disease_predicates:
        predicate_filter = list(DISEASE_MECHANISM_PREDICATES) + (predicate_filter or [])
    elif args.use_molecular_predicates:
        predicate_filter = list(MOLECULAR_INTERACTION_PREDICATES) + (predicate_filter or [])
    else:
        # Default: BIOLOGICAL_PREDICATES (also when --use-biological-predicates is set)
        predicate_filter = list(BIOLOGICAL_PREDICATES) + (predicate_filter or [])

    graph = build_transomic_property_graph(
        db_connection=db,
        cohort=args.cohort,
        sample_id=sample_id,
        omic_types=args.omic_types,
        value_abs_threshold=args.value_abs_threshold,
        value_threshold_per_omic=per_omic,
        top_n=args.top_n,
        predicate_filter=predicate_filter,
        edge_limit=args.edge_limit,
        induced_subgraph=not args.star_expansion,
        include_edges=True,
    )

    out_path = args.output
    if out_path is None:
        safe_sample = sample_id.replace("/", "-")
        out_path = os.path.join("..", "transomic-networks", f"transomic_graph_{args.cohort}_{safe_sample}.json")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(graph, f, indent=2)

    print(f"Saved transomic graph -> {out_path}")


#%% Quick-start cell for notebooks (edit variables then run)

# {"_key": "TCGA-E9-A1N8-01A", 
# {"_key": "TCGA-E9-A1N8-11A", 
{
  "01": "Primary Solid Tumor",
  "02": "Recurrent Solid Tumor",
  "03": "Primary Blood Derived Cancer – Peripheral Blood",
  "10": "Blood Derived Normal",
  "11": "Solid Tissue Normal",
  "06": "Metastatic",
  "07": "Additional New Primary"
}
#%% Entry point for CLI execution

if __name__ == "__main__":
    # Avoid argparse errors when executed inside notebooks (%run) by skipping CLI
    # if ipykernel is present. Users can call build_transomic_network_for_sample instead.
    if "ipykernel" not in sys.modules:
        print("Running as CLI script...")
        main()

    else:
        print("Running inside interactive environment; skipping CLI. Call build_transomic_network_for_sample() directly with desired parameters.")
        # Uncomment and execute in a notebook/VS Code interactive window
        DB_NAME = "PKT_test10000"
        COHORT = "TCGA-BRCA"
        Case = "TCGA-BH-A1F2" #TCGA-GM-A2DD
        Sample = "01"
        SAMPLE_ID =f"{Case}-{Sample}A"  # required: single sample only

        for predicate_set, tag in [
            (BIOLOGICAL_PREDICATES, "biological"),
            (MOLECULAR_INTERACTION_PREDICATES, "molecular"),
            (DISEASE_MECHANISM_PREDICATES, "disease"),
        ]:
            graph_obj, saved_path = build_transomic_network_for_sample(
                db_name=DB_NAME,
                cohort=COHORT,
                sample_id=SAMPLE_ID,
                omic_types=list(OMIC_CONFIG.keys()),
                value_threshold_per_omic=DEFAULT_VALUE_THRESHOLDS,
                predicate_filter=predicate_set,
                top_n=None,
                tag_flag=tag,
            )
            print(f"Graph saved to {saved_path}")
            print(f"Nodes: {len(graph_obj['nodes'])}, Edges: {len(graph_obj['edges'])}")

#%%
"""
Thresholded vs non trheholded
Graph saved to ..
Nodes: 78012, 
-----------------------------------
Graph saved to ..
Nodes: 123651, 
"""

"""
biological + clinical predicates
Nodes: 78012, Edges: 986792

molecular interaction predicates
Nodes: 78012, Edges: 737289

disease mechanism predicates
Nodes: 78012, Edges: 41457
"""

# %%
