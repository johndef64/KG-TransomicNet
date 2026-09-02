"""
Storage footprint of the deployed KG-TransomicNet instance.

Read-only. Reports, per collection, document count and on-disk size of
documents and indexes, grouped into the semantic backbone and the five
quantitative layers.

    python scripts/benchmark/bench_storage.py
    python scripts/benchmark/bench_storage.py --db MY_DATABASE
"""
from __future__ import annotations

import argparse

from bench_common import (add_db_argument, collection_figures, connect,
                          list_collections, write_csv, write_env)

BACKBONE = {"nodes", "edges"}
LAYERS = {
    "GENE_EXPRESSION": "RNA-Seq",
    "CNV": "CNV",
    "MIRNA": "miRNA",
    "METHYLATION": "Methylation",
    "PROTEIN": "RPPA",
}
SUPPORT = {"GENES", "SAMPLES", "PROJECTS", "CASES"}


def layer_breakdown(db, collection: str) -> dict:
    """Split a quantitative collection into index vs vector documents."""
    q = f"""
    FOR d IN {collection}
      COLLECT t = d.data_type WITH COUNT INTO n
      RETURN {{t: t, n: n}}
    """
    out = {"n_index": 0, "n_vector": 0}
    for row in db.aql.execute(q):
        if row["t"] and "index" in row["t"]:
            out["n_index"] += row["n"]
        elif row["t"] and "vector" in row["t"]:
            out["n_vector"] += row["n"]
    return out


def main() -> None:
    args = add_db_argument(argparse.ArgumentParser()).parse_args()
    db = connect(args.db)
    write_env()

    rows = []
    for name in list_collections(args.db):
        fig = collection_figures(args.db, name)
        if name in BACKBONE:
            group = "semantic_backbone"
        elif name in LAYERS:
            group = "quantitative_layer"
        else:
            group = "support"

        row = {
            "collection": name,
            "group": group,
            "layer": LAYERS.get(name, ""),
            "n_documents": fig["count"],
            "documents_mb": round(fig["documents_bytes"] / 1e6, 2),
            "indexes_mb": round(fig["indexes_bytes"] / 1e6, 2),
            "total_mb": round((fig["documents_bytes"] + fig["indexes_bytes"]) / 1e6, 2),
            "n_index_docs": "",
            "n_vector_docs": "",
        }
        if name in LAYERS:
            br = layer_breakdown(db, name)
            row["n_index_docs"] = br["n_index"]
            row["n_vector_docs"] = br["n_vector"]
        rows.append(row)

    for r in rows:
        print(f"  {r['collection']:20s} {r['group']:18s} "
              f"{r['n_documents']:>12,d} docs  {r['total_mb']:>9.1f} MB")

    tot = sum(r["total_mb"] for r in rows)
    bb = sum(r["total_mb"] for r in rows if r["group"] == "semantic_backbone")
    qt = sum(r["total_mb"] for r in rows if r["group"] == "quantitative_layer")
    print(f"\n  backbone {bb:,.1f} MB | quantitative {qt:,.1f} MB | total {tot:,.1f} MB")

    write_csv(rows, "storage_footprint.csv")


if __name__ == "__main__":
    main()
