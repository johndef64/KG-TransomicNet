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

#%% Imports
import argparse
import json
import os
import sys
from typing import List, Optional, Tuple

from arangodb_utils import setup_arangodb_connection
from query_utils import (
    OMIC_CONFIG,
    build_transomic_property_graph,
    list_samples_with_complete_omics,
)

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
        help="Keep features with |value| >= threshold",
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
        help="Optional list of predicate_label values to retain in KG edges",
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
    top_n: Optional[int] = None,
    predicate_filter: Optional[List[str]] = None,
    output_path: Optional[str] = None,
) -> Tuple[dict, str]:
    """
    Build and optionally persist a transomic property-graph for a single sample.

    Designed for notebook (#%%) execution: returns the graph object and the path used.
    The caller must provide a specific sample_id; no cohort-wide batching is done here.
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
        top_n=top_n,
        predicate_filter=predicate_filter,
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

    graph = build_transomic_property_graph(
        db_connection=db,
        cohort=args.cohort,
        sample_id=sample_id,
        omic_types=args.omic_types,
        value_abs_threshold=args.value_abs_threshold,
        top_n=args.top_n,
        predicate_filter=args.predicate_filter,
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

#%%
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
    value_abs_threshold=None,
    top_n=None,
    predicate_filter=None,
)
print(f"Graph saved to {saved_path}")

#%% Entry point for CLI execution
if __name__ == "__main__":
    # Avoid argparse errors when executed inside notebooks (%run) by skipping CLI
    # if ipykernel is present. Users can call build_transomic_network_for_sample instead.
    if "ipykernel" not in sys.modules:
        main()


# %%
