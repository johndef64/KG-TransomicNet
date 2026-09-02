# Database-engineering benchmarks

Measurements backing the performance-evaluation section of the manuscript.
All scripts write CSV to `results/benchmark/` together with an
`environment.json` recording the machine and the ArangoDB build, without which
the timings are not interpretable.

## Requirements

```
pip install python-arango pyarrow psutil
```

A running ArangoDB instance holding the deployed database. Connection
parameters are at the top of `bench_common.py` (defaults: `localhost:8529`,
database `PKT_main`). Every script accepts `--db` to target a differently
named instance without editing the sources:

```bash
python scripts/benchmark/bench_storage.py --db MY_DATABASE
```

## Scripts

| Script | Produces | Touches production? |
|---|---|---|
| `bench_storage.py` | `storage_footprint.csv` | read-only |
| `bench_query.py` | `query_latency.csv`, `index_build_time.csv` | read-only, unless `--create-indexes` |
| `bench_ingestion.py` | `ingestion_batch_sweep.csv`, `ingestion_scaling.csv` | read-only (writes to a scratch database) |
| `bench_schema.py` | `schema_comparison.csv` | read-only (writes to a scratch database) |

Write experiments run in a scratch database (`BENCH_scratch`) that is created
and dropped by the script. Pass `--keep` to inspect it afterwards.

## Running

```bash
python scripts/benchmark/bench_storage.py
python scripts/benchmark/bench_query.py                      # as-deployed indexes
python scripts/benchmark/bench_query.py --create-indexes     # + Section 3.5 ablation
python scripts/benchmark/bench_ingestion.py
python scripts/benchmark/bench_schema.py --genes 2000 --kset 200
```

`bench_query.py --create-indexes` creates persistent secondary indexes on
`nodes.bioentity_type`, `nodes.class_code`, `edges.predicate_label` and
`edges.predicate_class_code`, then re-runs the suite to produce a
with/without-index comparison. Index creation is additive;
`--drop-created-indexes` removes exactly the indexes the script created and
nothing else.

## What is measured

**Storage** — per-collection document count and on-disk size of documents and
indexes, split into semantic backbone, quantitative layers and support
collections. Sizes are the RocksDB estimates ArangoDB exposes through
`/_api/collection/{name}/figures`; they lag compaction, so small collections
measured immediately after a bulk insert can be reported a few per cent off.

**Query latency** — three retrieval modalities, mirroring the three use cases:

- `Q1` semantic step: predicate-stratified gene-pair retrieval,
- `Q2` quantitative step: multi-layer per-gene alignment through `*_index`,
- `Q3` end-to-end: phenotype-anchored traversal followed by vector slicing
  (the pattern printed as Listing 1 of the manuscript).

Each is executed once before any warm-up (`cold_first_ms`, a proxy for a cold
block cache) and then ten times after two warm-up runs, reported as
mean/p50/p95/min/max.

**Ingestion** — transaction-size sweep on two document shapes: backbone-shaped
documents (small and numerous) and per-sample vector documents (large and few).
Sweep A calibrates the 1,000-document transaction unit reported in Section 3.5.
The scaling run measures ingestion time and storage against the number of
per-sample vector documents.

**Schema comparison** — the same `G x S` block of TCGA-BRCA expression values
materialised under four layouts, each measured for ingestion time, on-disk size
and three access patterns:

| Layout | Description |
|---|---|
| `S1` | measurement-as-node: one document per (feature, sample) |
| `S2` | measurement-as-property: one document per feature, array over samples |
| `S3` | KG + external matrix: values in Parquet, joined application-side |
| `S6` | proposed: one index document per cohort, one dense vector per sample |

Access patterns: `A` all features for one sample, `B` one feature across all
samples, `C` k features across all samples (the UC2/UC3 pattern).

## Schema note

Backbone nodes carry `entity_id` together with `class_code`
(`class_code == "EntrezID"` identifies Entrez genes). There is no `entrez_id`
attribute on `nodes`; the join key against the omic `*_index` documents is
`entity_id`. Queries in these scripts follow the deployed schema.
