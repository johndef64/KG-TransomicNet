"""
to_omic_matrices.py
===================
Export the omic layers of one cohort as feature x sample matrices, in the layout
of the files distributed by UCSC Xena: first column = feature identifier (with
the original header), one column per sample, tab-separated, gzip-compressed.
Downstream preprocessing can then read the export exactly as it would read the
original download.

For each layer of the layout file the adapter runs extraction/layer_index.aql
(row identifiers, in index order) and extraction/layer_vectors.aql (one vector
per sample); an optional label table is written from extraction/sample_labels.aql.
Every sample of the cohort is exported: intersections across layers, filtering
and normalisation are left to the downstream pipeline.

Layouts: params/omic_layers_pancan_atlas.json (PanCanAtlas view),
params/omic_layers_gdc.json (a GDC cohort; override it with --cohort).

Usage
-----
    python queries/adapters/to_omic_matrices.py --layout queries/params/omic_layers_pancan_atlas.json \
        --out-dir raw/
    python queries/adapters/to_omic_matrices.py --layout queries/params/omic_layers_gdc.json \
        --cohort TCGA-LUAD --layers star_tpm.tsv.gz --out-dir luad/
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import QUERIES_DIR, add_connection_args, connect_from_args, read_json, run  # noqa: E402

Q = QUERIES_DIR / "extraction"


def export_layer(db, layer, cohort, out_dir):
    index = layer["index"].format(cohort=cohort)
    feats = [r["id"] for r in run(db, Q / "layer_index.aql", {
        "index": index, "mapping": layer["mapping"], "id_field": layer["id_field"]})]
    if not feats:
        print(f"  [{layer['file']}] no index {index}: skipped")
        return None
    cols, empty = {}, []
    for r in run(db, Q / "layer_vectors.aql", {
            "@collection": layer["collection"], "cohort": cohort, "data_type": layer["data_type"],
            "value_field": layer["value_field"], "positions": None}, batch_size=100):
        vals = r["values"] or []
        if len(vals) != len(feats):
            raise SystemExit(f"{layer['file']}: sample {r['sample_id']} has {len(vals)} values, "
                             f"index has {len(feats)} features")
        if all(v is None for v in vals):
            empty.append(r["sample_id"])      # a vector with no measured value is not a sample
            continue
        cols[r["sample_id"]] = vals
    if empty:
        print(f"  [{layer['file']}] skipped {len(empty)} vector(s) without any value: {empty[:5]}")
    df = pd.DataFrame(cols, index=pd.Index(feats, name=layer["header"]))
    df = df[sorted(df.columns)]
    df.to_csv(out_dir / layer["file"], sep="\t", na_rep="NA", compression="gzip", lineterminator="\n")
    print(f"  [{layer['file']}] {df.shape[0]:,} features x {df.shape[1]:,} samples")
    return df.shape


def export_labels(db, spec, cohort, out_dir):
    rows = [(r["sample_id"], r["label"]) for r in run(db, Q / "sample_labels.aql", {
        "cohort": cohort, "label_field": spec["label_field"]})]
    pd.DataFrame(rows, columns=spec["columns"]).to_csv(
        out_dir / spec["file"], sep="\t", index=False, compression="gzip", lineterminator="\n")
    print(f"  [{spec['file']}] {len(rows):,} labelled samples")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--layout", required=True, help="layer layout (JSON)")
    p.add_argument("--cohort", help="override the cohort of the layout")
    p.add_argument("--layers", help="comma-separated subset of layer files to export")
    p.add_argument("--no-labels", action="store_true")
    p.add_argument("--out-dir", default="omic_matrices")
    add_connection_args(p)
    args = p.parse_args()

    spec = read_json(args.layout)
    cohort = args.cohort or spec["cohort"]
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    wanted = set(args.layers.split(",")) if args.layers else None
    db = connect_from_args(args)

    print(f"[to_omic_matrices] cohort {cohort}")
    for layer in spec["layers"]:
        if wanted is None or layer["file"] in wanted:
            export_layer(db, layer, cohort, out)
    if spec.get("labels") and not args.no_labels:
        export_labels(db, spec["labels"], cohort, out)
    print(f"[to_omic_matrices] -> {out}")


if __name__ == "__main__":
    main()
