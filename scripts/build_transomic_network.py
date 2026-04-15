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
        help="List of predicate_label values to retain in KG edges",
    )
    parser.add_argument(
        "--use-biological-predicates",
        action="store_true",
        help="Apply BIOLOGICAL_PREDICATES preset as predicate_filter",
    )
    parser.add_argument(
        "--edge-limit",
        type=int,
        default=None,
        help="Optional hard cap on KG edges (AQL LIMIT, no ranking). Default: no cap.",
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
    output_path: Optional[str] = None,
) -> Tuple[dict, str]:
    """
    Build and persist a transomic property-graph for a single sample.

    Key filters:
    - value_threshold_per_omic: dict {omic: threshold} — preferred way to drop
      low-signal features. Scales differ across omics so prefer this over
      value_abs_threshold.
    - predicate_filter: list of predicate_label values to retain in KG edges.
      Use BIOLOGICAL_PREDICATES for a topology-focused, biologically-relevant
      subgraph.
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
        include_edges=True,
    )

    out_path = output_path
    if out_path is None:
        safe_sample = sample_id.replace("/", "-")
        out_path = os.path.join("..", "transomic-networks", f"transomic_graph_{cohort}_{safe_sample}.json")
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
    if args.use_biological_predicates:
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
        main()

    else:

        # Uncomment and execute in a notebook/VS Code interactive window
        DB_NAME = "PKT_test10000"
        COHORT = "TCGA-BRCA"
        Case = "TCGA-BH-A1F2" #TCGA-GM-A2DD
        Sample = "01"
        SAMPLE_ID =f"{Case}-{Sample}A"  # required: single sample only

        graph_obj, saved_path = build_transomic_network_for_sample(
            db_name=DB_NAME,
            cohort=COHORT,
            sample_id=SAMPLE_ID,
            omic_types=list(OMIC_CONFIG.keys()),
            # value_threshold_per_omic=DEFAULT_VALUE_THRESHOLDS,
            predicate_filter=BIOLOGICAL_PREDICATES,
            top_n=None,
        )
        print(f"Graph saved to {saved_path}")
        print(f"Nodes: {len(graph_obj['nodes'])}, Edges: {len(graph_obj['edges'])}")


# %%
