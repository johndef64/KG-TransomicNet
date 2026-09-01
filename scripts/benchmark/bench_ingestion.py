"""
Ingestion benchmark: transaction-size sweep and scaling in the number of samples.

Sweep A  backbone-shaped documents (small, many)  -- real `nodes` documents
Sweep B  vector documents (large, few)            -- real per-sample expression vectors
Scaling  ingestion time and storage as the number of per-sample vector
         documents grows, which is the O(N_samples + N_cohorts) claim

Sweep A calibrates the 1,000-document transaction unit reported in
Section 3.5, which until now was stated without supporting measurements.

Documents are read from the production instance and re-inserted into a scratch
database, so the measured throughput is on real document shapes and sizes.
The production database is only read from.

    python scripts/benchmark/bench_ingestion.py
"""
from __future__ import annotations

import argparse
import time

from bench_common import (SCRATCH_DB, collection_figures, connect, drop_scratch,
                          make_scratch, write_csv, write_env)

COHORT = "TCGA-BRCA"
BATCH_SIZES = [100, 500, 1000, 2500, 5000, 10000]
VECTOR_BATCH_SIZES = [1, 5, 10, 25, 50]


def fetch_nodes(db, n: int) -> list[dict]:
    q = "FOR d IN nodes LIMIT @n RETURN UNSET(d, '_id', '_rev')"
    return list(db.aql.execute(q, bind_vars={"n": n}))


def fetch_vectors(db, n: int) -> list[dict]:
    q = """
    FOR d IN GENE_EXPRESSION
      FILTER d.data_type == "gene_expression_vector" AND d.cohort == @c
      LIMIT @n
      RETURN {_key: d._key, sample_id: d.sample_id, cohort: d.cohort,
              data_type: d.data_type, n_genes: d.n_genes,
              expression_index_ref: d.expression_index_ref,
              values_tpm: d.values_tpm}
    """
    return list(db.aql.execute(q, bind_vars={"c": COHORT, "n": n}))


def ingest(sdb, name: str, docs: list[dict], batch: int) -> float:
    if sdb.has_collection(name):
        sdb.delete_collection(name)
    col = sdb.create_collection(name)
    t0 = time.perf_counter()
    for i in range(0, len(docs), batch):
        col.insert_many(docs[i:i + batch])
    return time.perf_counter() - t0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-nodes", type=int, default=100_000)
    ap.add_argument("--n-vectors", type=int, default=300)
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()

    write_env()
    db = connect()

    print(f"[fetching {args.n_nodes:,} backbone documents]")
    nodes = fetch_nodes(db, args.n_nodes)
    print(f"[fetching {args.n_vectors} per-sample vector documents]")
    vectors = fetch_vectors(db, args.n_vectors)
    vec_floats = len(vectors[0]["values_tpm"]) if vectors else 0
    print(f"  vector documents carry {vec_floats:,} values each")

    sdb = make_scratch()
    rows = []

    print("\n[Sweep A: backbone-shaped documents]")
    for b in BATCH_SIZES:
        dt = ingest(sdb, "sweepA", nodes, b)
        fig = collection_figures(SCRATCH_DB, "sweepA")
        rows.append({"sweep": "A_backbone_nodes", "batch_size": b,
                     "n_documents": len(nodes), "elapsed_s": round(dt, 2),
                     "docs_per_s": round(len(nodes) / dt, 1),
                     "storage_mb": round(fig["documents_bytes"] / 1e6, 2)})
        print(f"  batch {b:>6,d}  {dt:7.2f} s  {len(nodes)/dt:>10,.0f} docs/s")

    print("\n[Sweep B: per-sample vector documents]")
    for b in VECTOR_BATCH_SIZES:
        dt = ingest(sdb, "sweepB", vectors, b)
        fig = collection_figures(SCRATCH_DB, "sweepB")
        mb = fig["documents_bytes"] / 1e6
        rows.append({"sweep": "B_vector_documents", "batch_size": b,
                     "n_documents": len(vectors), "elapsed_s": round(dt, 2),
                     "docs_per_s": round(len(vectors) / dt, 1),
                     "storage_mb": round(mb, 2)})
        print(f"  batch {b:>6,d}  {dt:7.2f} s  {len(vectors)/dt:>10,.1f} docs/s"
              f"  {mb/dt:>7.1f} MB/s")

    write_csv(rows, "ingestion_batch_sweep.csv")

    print("\n[Scaling: vector documents vs ingestion time and storage]")
    scal = []
    steps = [n for n in (50, 100, 150, 200, 250, 300) if n <= len(vectors)]
    for n in steps:
        dt = ingest(sdb, "scaling", vectors[:n], 25)
        fig = collection_figures(SCRATCH_DB, "scaling")
        scal.append({"n_sample_documents": n, "n_documents_total": fig["count"],
                     "elapsed_s": round(dt, 2),
                     "storage_mb": round(fig["documents_bytes"] / 1e6, 2),
                     "ms_per_document": round(dt * 1000 / n, 2),
                     "mb_per_document": round(fig["documents_bytes"] / 1e6 / n, 3)})
        print(f"  {n:>4d} samples  {dt:6.2f} s  "
              f"{fig['documents_bytes']/1e6:7.1f} MB  "
              f"{dt*1000/n:6.1f} ms/doc")
    write_csv(scal, "ingestion_scaling.csv")

    if not args.keep:
        drop_scratch()
        print("  scratch database removed")


if __name__ == "__main__":
    main()
