"""
Query-latency benchmark for the three retrieval modalities of the framework.

Q1  semantic step      predicate-stratified gene-pair retrieval          (UC1)
Q2  quantitative step  multi-layer per-gene alignment via *_index        (UC2)
Q3  end-to-end         phenotype-anchored traversal + vector slicing     (UC3)

Q3 is the two-step pattern printed as Listing 1 of the manuscript.

By default the script is READ-ONLY. With --create-indexes it additionally
creates the secondary indexes described in Section 3.5 of the manuscript and
re-runs the suite, producing a with/without-index ablation. Index creation is
additive and reversible (--drop-created-indexes removes exactly those).

    python scripts/benchmark/bench_query.py
    python scripts/benchmark/bench_query.py --create-indexes
    python scripts/benchmark/bench_query.py --drop-created-indexes

NOTE on the schema: backbone nodes carry `entity_id` + `class_code`
(class_code == 'EntrezID' identifies Entrez genes). There is no `entrez_id`
attribute on `nodes`; the join key against the omic indexes is `entity_id`.
"""
from __future__ import annotations

import argparse
import time

from bench_common import (add_db_argument, connect, time_query, write_csv,
                          write_env)

COHORT = "TCGA-BRCA"
REPEATS = 10

# Indexes described in Section 3.5 of the manuscript. Created only on request.
DECLARED_INDEXES = [
    ("nodes", ["bioentity_type"]),
    ("nodes", ["class_code"]),
    ("edges", ["predicate_label"]),
    ("edges", ["predicate_class_code"]),
]

# --- Q1: semantic step -----------------------------------------------------
# In the PKT instance-based OWLNETS build only `genetically interacts with`
# connects two Entrez genes directly; `molecularly interacts with` relates
# proteins and `has participant` relates pathways to their members. UC1 gene
# pairs are therefore obtained in two ways, and both are benchmarked:
#
#   Q1a  direct     seed gene panel -> predicate-filtered gene neighbours
#   Q1b  2-hop      seed gene panel -> shared intermediary -> gene
#
# Both are anchored on a bounded seed panel, which is the access pattern UC1
# actually uses (the hop-distance analysis of Section 4.3.1).
Q1A = """
FOR g IN nodes
  FILTER g.class_code == "EntrezID"
  LIMIT @nseed
  FOR v, e IN 1..1 ANY g edges
    FILTER e.predicate_label == @pred AND v.class_code == "EntrezID"
    RETURN {a: g.entity_id, b: v.entity_id}
"""

Q1B = """
FOR g IN nodes
  FILTER g.class_code == "EntrezID"
  LIMIT @nseed
  FOR m IN 1..1 ANY g edges
    FILTER m.bioentity_type == @mid
    FOR g2 IN 1..1 ANY m edges
      FILTER g2.class_code == "EntrezID" AND g2._key != g._key
      LIMIT @lim
      RETURN {a: g.entity_id, b: g2.entity_id}
"""

# --- Q2: quantitative step, three layers -----------------------------------
Q2 = """
LET pidx = DOCUMENT(CONCAT("PROTEIN/protein_index_", @cohort))
LET targets = (FOR m IN pidx.protein_mappings
                 FILTER m.entrez_id != null RETURN DISTINCT m.entrez_id)

LET eidx = DOCUMENT(CONCAT("GENE_EXPRESSION/expr_index_", @cohort))
LET epos = (FOR m IN eidx.gene_mappings
              FILTER m.entrez_id IN targets
              RETURN {e: m.entrez_id, p: m.position})

LET cidx = DOCUMENT(CONCAT("CNV/cnv_index_", @cohort))
LET cpos = (FOR m IN cidx.gene_mappings
              FILTER m.entrez_id IN targets
              RETURN {e: m.entrez_id, p: m.position})

LET expr = (FOR s IN GENE_EXPRESSION
              FILTER s.data_type == "gene_expression_vector" AND s.cohort == @cohort
              LIMIT @nsamples
              RETURN {sample: s.sample_id,
                      v: (FOR x IN epos RETURN s.values_tpm[x.p])})

LET cnv = (FOR s IN CNV
             FILTER s.data_type == "cnv_vector" AND s.cohort == @cohort
             LIMIT @nsamples
             RETURN {sample: s.sample_id,
                     v: (FOR x IN cpos RETURN s.values_copy_number[x.p])})

RETURN {n_targets: LENGTH(targets), n_expr: LENGTH(expr), n_cnv: LENGTH(cnv)}
"""

# --- Q3: end-to-end, Listing 1 of the manuscript ---------------------------
Q3 = """
LET targets = (
  FOR v IN 1..@depth ANY @seed edges
    OPTIONS {bfs: true, uniqueVertices: "global"}
    FILTER v.bioentity_type == "gene" AND v.class_code == "EntrezID"
    RETURN DISTINCT v.entity_id)

LET idx = DOCUMENT(CONCAT("GENE_EXPRESSION/expr_index_", @cohort))
LET pos = (FOR m IN idx.gene_mappings
             FILTER m.entrez_id IN targets
             RETURN {e: m.entrez_id, p: m.position})

LET sliced = (
  FOR s IN GENE_EXPRESSION
    FILTER s.data_type == "gene_expression_vector" AND s.cohort == @cohort
    LIMIT @nsamples
    RETURN {sample: s.sample_id,
            tpm: (FOR x IN pos RETURN [x.e, s.values_tpm[x.p]])})

RETURN {n_targets: LENGTH(targets), n_pos: LENGTH(pos), n_samples: LENGTH(sliced)}
"""

CASES = [
    ("Q1a_semantic_direct", "genetically interacts with / 200 seeds", Q1A,
     {"pred": "genetically interacts with", "nseed": 200}),
    ("Q1b_semantic_2hop", "co-pathway / 50 seeds", Q1B,
     {"mid": "pathway", "nseed": 50, "lim": 5000}),
    ("Q1b_semantic_2hop", "co-protein / 50 seeds", Q1B,
     {"mid": "protein", "nseed": 50, "lim": 5000}),
    ("Q2_quantitative", "3-layer alignment / 878 samples", Q2,
     {"cohort": COHORT, "nsamples": 878}),
    ("Q3_end_to_end", "HP:0003002 depth=1", Q3,
     {"seed": "nodes/HP_0003002", "depth": 1, "cohort": COHORT, "nsamples": 878}),
    ("Q3_end_to_end", "HP:0003002 depth=2", Q3,
     {"seed": "nodes/HP_0003002", "depth": 2, "cohort": COHORT, "nsamples": 878}),
]


def existing_index(db, collection, fields):
    for idx in db.collection(collection).indexes():
        if idx.get("fields") == fields and idx["type"] in ("persistent", "hash"):
            return idx
    return None


def create_declared_indexes(db):
    created = []
    for coll, fields in DECLARED_INDEXES:
        if existing_index(db, coll, fields):
            print(f"  index on {coll}{fields} already present, skipping")
            continue
        print(f"  creating persistent index on {coll}{fields} ...", flush=True)
        t0 = time.perf_counter()
        db.collection(coll).add_index(
            {"type": "persistent", "fields": fields, "sparse": False,
             "name": f"bench_{coll}_{'_'.join(fields)}"})
        dt = time.perf_counter() - t0
        print(f"    done in {dt:.1f}s")
        created.append((coll, fields, round(dt, 2)))
    return created


def drop_created_indexes(db):
    for coll, fields in DECLARED_INDEXES:
        name = f"bench_{coll}_{'_'.join(fields)}"
        for idx in db.collection(coll).indexes():
            if idx.get("name") == name:
                db.collection(coll).delete_index(idx["id"])
                print(f"  dropped {name}")


def run_suite(db, index_state: str) -> list[dict]:
    rows = []
    for qid, variant, aql, bind in CASES:
        # first execution before any warm-up: proxy for a cold block cache
        t0 = time.perf_counter()
        list(db.aql.execute(aql, bind_vars=bind))
        cold_ms = round((time.perf_counter() - t0) * 1000.0, 2)

        stats = time_query(db, aql, bind, repeats=REPEATS, warmup=2)
        row = {"query": qid, "variant": variant, "index_state": index_state,
               "cold_first_ms": cold_ms, **stats}
        rows.append(row)
        print(f"  {qid:16s} {variant:32s} [{index_state:12s}] "
              f"cold {cold_ms:>9.1f} ms | warm p50 {stats['p50_ms']:>9.1f} ms")
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--create-indexes", action="store_true",
                    help="create the Section 3.5 indexes and re-run the suite")
    ap.add_argument("--drop-created-indexes", action="store_true",
                    help="remove only the indexes this script created, then exit")
    add_db_argument(ap)
    args = ap.parse_args()

    db = connect(args.db)

    if args.drop_created_indexes:
        drop_created_indexes(db)
        return

    write_env()
    print("\n[baseline: indexes as deployed]")
    rows = run_suite(db, "as_deployed")

    if args.create_indexes:
        print("\n[creating Section 3.5 secondary indexes]")
        created = create_declared_indexes(db)
        if created:
            write_csv([{"collection": c, "fields": "|".join(f), "build_s": s}
                       for c, f, s in created], "index_build_time.csv")
        print("\n[with declared indexes]")
        rows += run_suite(db, "with_declared_idx")

    write_csv(rows, "query_latency.csv")


if __name__ == "__main__":
    main()
