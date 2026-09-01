"""
Storage-layout comparison: the proposed vector+index model against the
alternative representations available for the same information content.

  S1  measurement-as-node      one document per (feature, sample) measurement
  S2  measurement-as-property  one document per feature, array over samples
  S3  KG + external matrix     values in an external columnar file (Parquet),
                               joined application-side  <- the de-facto baseline
  S6  vector + index           proposed: one index document per cohort,
                               one dense vector document per sample

All four hold exactly the same values: G features x S samples of TCGA-BRCA
gene expression. Reported per layout: ingestion time, on-disk size, and the
latency of three access patterns.

  A  sample slice    all G features for one sample
  B  feature slice   one feature across all S samples
  C  gene-set slice  k features across all S samples  (the UC2/UC3 pattern)

Everything runs in a scratch database that is dropped at the end; the
production database is only read from.

    python scripts/benchmark/bench_schema.py [--genes 2000] [--keep]
"""
from __future__ import annotations

import argparse
import shutil
import statistics
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from bench_common import (SCRATCH_DB, collection_figures, connect, drop_scratch,
                          make_scratch, write_csv, write_env)

COHORT = "TCGA-BRCA"
TMP_DIR = Path(__file__).resolve().parent / "_tmp_parquet"


def fetch_matrix(db, n_genes: int):
    """Pull a G x S slice of the BRCA expression layer from the live instance."""
    idx = db.aql.execute("""
        LET idx = DOCUMENT(CONCAT("GENE_EXPRESSION/expr_index_", @c))
        FOR m IN idx.gene_mappings
          FILTER m.entrez_id != null
          LIMIT @g
          RETURN {e: m.entrez_id, p: m.position, s: m.hgnc_symbol}
    """, bind_vars={"c": COHORT, "g": n_genes})
    mappings = list(idx)
    positions = [m["p"] for m in mappings]

    rows = db.aql.execute("""
        FOR s IN GENE_EXPRESSION
          FILTER s.data_type == "gene_expression_vector" AND s.cohort == @c
          RETURN {sample: s.sample_id, v: (FOR p IN @pos RETURN s.values_tpm[p])}
    """, bind_vars={"c": COHORT, "pos": positions})
    samples, values = [], []
    for r in rows:
        samples.append(r["sample"])
        values.append([0.0 if x is None else float(x) for x in r["v"]])
    return mappings, samples, values


def timed(fn, repeats: int = 5, warmup: int = 1):
    for _ in range(warmup):
        fn()
    lat = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        lat.append((time.perf_counter() - t0) * 1000.0)
    lat.sort()
    return round(statistics.median(lat), 2)


# --- S6: vector + index ----------------------------------------------------
def build_s6(sdb, mappings, samples, values):
    col = sdb.create_collection("S6_expression")
    col.add_index({"type": "persistent", "fields": ["data_type"]})
    t0 = time.perf_counter()
    col.insert({"_key": "idx", "data_type": "index",
                "gene_mappings": [{"position": i, "entrez_id": m["e"]}
                                  for i, m in enumerate(mappings)]})
    docs = [{"_key": s, "data_type": "vector", "index_ref": "idx", "values": v}
            for s, v in zip(samples, values)]
    col.insert_many(docs)
    return round(time.perf_counter() - t0, 2)


def query_s6(sdb, n_genes, sample_key, gene_pos, kset):
    A = lambda: list(sdb.aql.execute(
        "RETURN DOCUMENT('S6_expression/@@k').values", bind_vars={}) ) # placeholder
    def a():
        list(sdb.aql.execute("FOR d IN S6_expression FILTER d._key == @k RETURN d.values",
                             bind_vars={"k": sample_key}))
    def b():
        list(sdb.aql.execute("""
            FOR d IN S6_expression FILTER d.data_type == "vector"
              RETURN d.values[@p]""", bind_vars={"p": gene_pos}))
    def c():
        list(sdb.aql.execute("""
            FOR d IN S6_expression FILTER d.data_type == "vector"
              RETURN (FOR p IN @ps RETURN d.values[p])""", bind_vars={"ps": kset}))
    return timed(a), timed(b), timed(c)


# --- S1: measurement-as-node ----------------------------------------------
def build_s1(sdb, mappings, samples, values):
    col = sdb.create_collection("S1_measurements")
    col.add_index({"type": "persistent", "fields": ["sample"]})
    col.add_index({"type": "persistent", "fields": ["entrez_id"]})
    t0 = time.perf_counter()
    batch = []
    for si, s in enumerate(samples):
        for gi, m in enumerate(mappings):
            batch.append({"sample": s, "entrez_id": m["e"], "value": values[si][gi]})
            if len(batch) >= 50000:
                col.insert_many(batch)
                batch = []
    if batch:
        col.insert_many(batch)
    return round(time.perf_counter() - t0, 2)


def query_s1(sdb, sample_key, entrez, kset_ids):
    def a():
        list(sdb.aql.execute("FOR d IN S1_measurements FILTER d.sample == @s RETURN d.value",
                             bind_vars={"s": sample_key}))
    def b():
        list(sdb.aql.execute("FOR d IN S1_measurements FILTER d.entrez_id == @e RETURN d.value",
                             bind_vars={"e": entrez}))
    def c():
        list(sdb.aql.execute("""
            FOR d IN S1_measurements FILTER d.entrez_id IN @es
              RETURN {s: d.sample, v: d.value}""", bind_vars={"es": kset_ids}))
    return timed(a, repeats=3), timed(b, repeats=3), timed(c, repeats=3)


# --- S2: measurement-as-property ------------------------------------------
def build_s2(sdb, mappings, samples, values):
    col = sdb.create_collection("S2_features")
    t0 = time.perf_counter()
    docs = [{"_key": f"g{i}", "entrez_id": m["e"], "samples": samples,
             "values": [values[si][i] for si in range(len(samples))]}
            for i, m in enumerate(mappings)]
    col.insert_many(docs)
    return round(time.perf_counter() - t0, 2)


def query_s2(sdb, sample_pos, gene_key, kset_keys):
    def a():
        list(sdb.aql.execute("FOR d IN S2_features RETURN d.values[@p]",
                             bind_vars={"p": sample_pos}))
    def b():
        list(sdb.aql.execute("FOR d IN S2_features FILTER d._key == @k RETURN d.values",
                             bind_vars={"k": gene_key}))
    def c():
        list(sdb.aql.execute("FOR d IN S2_features FILTER d._key IN @ks RETURN d.values",
                             bind_vars={"ks": kset_keys}))
    return timed(a), timed(b), timed(c)


# --- S3: KG + external matrix (Parquet, application-side join) -------------
def build_s3(mappings, samples, values):
    TMP_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    cols = {"sample": pa.array(samples)}
    for gi, m in enumerate(mappings):
        cols[m["e"]] = pa.array([values[si][gi] for si in range(len(samples))])
    pq.write_table(pa.table(cols), TMP_DIR / "expr.parquet", compression="snappy")
    return round(time.perf_counter() - t0, 2)


def query_s3(samples, entrez, kset_ids, sample_key):
    path = TMP_DIR / "expr.parquet"
    def a():
        t = pq.read_table(path, filters=[("sample", "=", sample_key)])
        t.to_pydict()
    def b():
        pq.read_table(path, columns=["sample", entrez]).to_pydict()
    def c():
        pq.read_table(path, columns=["sample"] + kset_ids).to_pydict()
    return timed(a), timed(b), timed(c)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes", type=int, default=2000)
    ap.add_argument("--kset", type=int, default=200)
    ap.add_argument("--keep", action="store_true", help="do not drop the scratch DB")
    args = ap.parse_args()

    write_env()
    db = connect()
    print(f"[fetching {args.genes} genes x all BRCA samples from the live instance]")
    mappings, samples, values = fetch_matrix(db, args.genes)
    G, S = len(mappings), len(samples)
    print(f"  matrix: {G} features x {S} samples = {G*S:,} measurements")

    sdb = make_scratch()
    kset_idx = list(range(0, G, max(1, G // args.kset)))[:args.kset]
    kset_ids = [mappings[i]["e"] for i in kset_idx]
    kset_keys = [f"g{i}" for i in kset_idx]

    rows = []

    print("[S6 vector + index]")
    t = build_s6(sdb, mappings, samples, values)
    a, b, c = query_s6(sdb, G, samples[0], kset_idx[0], kset_idx)
    f = collection_figures(SCRATCH_DB, "S6_expression")
    rows.append({"layout": "S6_vector_index", "n_documents": f["count"],
                 "ingest_s": t, "storage_mb": round((f["documents_bytes"]+f["indexes_bytes"])/1e6, 2),
                 "A_sample_slice_ms": a, "B_feature_slice_ms": b, "C_geneset_slice_ms": c})
    print(f"  docs={f['count']:,} ingest={t}s A={a}ms B={b}ms C={c}ms")

    print("[S2 measurement-as-property]")
    t = build_s2(sdb, mappings, samples, values)
    a, b, c = query_s2(sdb, 0, kset_keys[0], kset_keys)
    f = collection_figures(SCRATCH_DB, "S2_features")
    rows.append({"layout": "S2_measurement_as_property", "n_documents": f["count"],
                 "ingest_s": t, "storage_mb": round((f["documents_bytes"]+f["indexes_bytes"])/1e6, 2),
                 "A_sample_slice_ms": a, "B_feature_slice_ms": b, "C_geneset_slice_ms": c})
    print(f"  docs={f['count']:,} ingest={t}s A={a}ms B={b}ms C={c}ms")

    print("[S3 KG + external matrix (Parquet)]")
    t = build_s3(mappings, samples, values)
    a, b, c = query_s3(samples, kset_ids[0], kset_ids, samples[0])
    size = (TMP_DIR / "expr.parquet").stat().st_size
    rows.append({"layout": "S3_kg_plus_external_matrix", "n_documents": 1,
                 "ingest_s": t, "storage_mb": round(size/1e6, 2),
                 "A_sample_slice_ms": a, "B_feature_slice_ms": b, "C_geneset_slice_ms": c})
    print(f"  parquet ingest={t}s A={a}ms B={b}ms C={c}ms")

    print("[S1 measurement-as-node]  (slowest, be patient)")
    t = build_s1(sdb, mappings, samples, values)
    a, b, c = query_s1(sdb, samples[0], kset_ids[0], kset_ids)
    f = collection_figures(SCRATCH_DB, "S1_measurements")
    rows.append({"layout": "S1_measurement_as_node", "n_documents": f["count"],
                 "ingest_s": t, "storage_mb": round((f["documents_bytes"]+f["indexes_bytes"])/1e6, 2),
                 "A_sample_slice_ms": a, "B_feature_slice_ms": b, "C_geneset_slice_ms": c})
    print(f"  docs={f['count']:,} ingest={t}s A={a}ms B={b}ms C={c}ms")

    for r in rows:
        r["n_features"] = G
        r["n_samples"] = S
        r["n_measurements"] = G * S
    write_csv(rows, "schema_comparison.csv")

    if not args.keep:
        drop_scratch()
        shutil.rmtree(TMP_DIR, ignore_errors=True)
        print("  scratch database and temporary files removed")


if __name__ == "__main__":
    main()
