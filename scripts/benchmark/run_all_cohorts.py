"""
Per-cohort ingestion of the full TCGA/TARGET set, with measurement.

For every study not yet loaded, runs the three existing pipeline stages and
records what each one cost:

    download_omics.py            --studies STUDY
    build_omics_collections.py   --studies STUDY
    load_omics_collections_to_arangodb.py --studies STUDY   <- timed

After a successful load the study's raw downloads and intermediate JSON
collections are deleted, so peak disk stays at roughly one cohort's working set
plus the growing database instead of the ~29 GB a naive download-all/build-all
/load-all pass would need.

Results are appended to results/benchmark/cohort_ingestion.csv after every
cohort, so an interrupted run loses nothing and re-running resumes where it
stopped. TCGA-BRCA is skipped: it is already loaded and its local files are
left untouched.

    python scripts/benchmark/run_all_cohorts.py            # all remaining
    python scripts/benchmark/run_all_cohorts.py --dry-run  # show the plan
    python scripts/benchmark/run_all_cohorts.py --studies TCGA-CHOL TCGA-UVM
    python scripts/benchmark/run_all_cohorts.py --no-cleanup
"""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from bench_common import (DB_NAME, RESULTS_DIR, add_db_argument,
                          collection_figures, connect)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = PROJECT_ROOT / "scripts"
OMICS_ROOT = PROJECT_ROOT / "data" / "omics"
COLLECTIONS_ROOT = PROJECT_ROOT / "data" / "arangodb_collections"
CSV_PATH = RESULTS_DIR / "cohort_ingestion.csv"
LOG_PATH = RESULTS_DIR / "cohort_ingestion.log"

ALREADY_LOADED = {"TCGA-BRCA"}

LAYER_COLLECTIONS = {
    "GENE_EXPRESSION": "gene_expression_vector",
    "CNV": "cnv_vector",
    "MIRNA": "mirna_vector",
    "PROTEIN": "protein_vector",
    "METHYLATION": "methylation_vector",
}

FIELDS = ["study", "status", "download_s", "build_s", "load_s",
          "n_expression", "n_cnv", "n_mirna", "n_protein", "n_methylation",
          "n_samples_total", "db_delta_mb", "db_total_mb", "timestamp", "note"]


def log(msg: str) -> None:
    line = f"[{datetime.now():%H:%M:%S}] {msg}"
    print(line, flush=True)
    LOG_PATH.parent.mkdir(parents=True, exist_ok=True)
    with LOG_PATH.open("a", encoding="utf-8") as fh:
        fh.write(line + "\n")


def study_lists() -> list[str]:
    """Read the canonical study lists from the pipeline itself."""
    if str(SCRIPTS) not in sys.path:
        sys.path.insert(0, str(SCRIPTS))          # the loader imports siblings
    from download_omics import TARGET_STUDIES, TCGA_STUDIES
    return list(TCGA_STUDIES) + list(TARGET_STUDIES)


def run_stage(script: str, study: str, extra: list[str] | None = None,
              timeout: int = 7200) -> tuple[bool, float, str]:
    """Run one pipeline stage for one study. Returns (ok, elapsed_s, tail)."""
    cmd = [sys.executable, str(SCRIPTS / script), "--studies", study]
    if extra:
        cmd += extra
    t0 = time.perf_counter()
    try:
        p = subprocess.run(cmd, cwd=PROJECT_ROOT, capture_output=True,
                           text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return False, time.perf_counter() - t0, "timeout"
    dt = time.perf_counter() - t0
    if p.returncode != 0:
        tail = (p.stderr or p.stdout or "").strip().splitlines()
        return False, dt, " | ".join(tail[-3:])[:300]
    return True, dt, ""


def db_total_mb(db_name: str = DB_NAME) -> float:
    total = 0.0
    for c in list(LAYER_COLLECTIONS) + ["SAMPLES", "CASES", "PROJECTS", "GENES"]:
        try:
            f = collection_figures(db_name, c)
            total += (f["documents_bytes"] + f["indexes_bytes"]) / 1e6
        except Exception:
            pass
    return round(total, 2)


def layer_counts(db, study: str) -> dict:
    out = {}
    for coll, dtype in LAYER_COLLECTIONS.items():
        q = f"""RETURN COUNT(FOR d IN {coll}
                 FILTER d.cohort == @s AND d.data_type == @t RETURN 1)"""
        try:
            out[coll] = list(db.aql.execute(q, bind_vars={"s": study, "t": dtype}))[0]
        except Exception:
            out[coll] = 0
    return out


def done_studies() -> set[str]:
    if not CSV_PATH.exists():
        return set()
    with CSV_PATH.open(encoding="utf-8") as fh:
        return {r["study"] for r in csv.DictReader(fh) if r.get("status") == "ok"}


def append_row(row: dict) -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    new = not CSV_PATH.exists()
    with CSV_PATH.open("a", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        if new:
            w.writeheader()
        w.writerow(row)


def cleanup(study: str) -> None:
    for path in (OMICS_ROOT / study, COLLECTIONS_ROOT / study):
        if path.exists():
            shutil.rmtree(path, ignore_errors=True)
            log(f"    removed {path.relative_to(PROJECT_ROOT)}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--studies", nargs="+", help="restrict to these studies")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-cleanup", action="store_true",
                    help="keep raw downloads and intermediate JSON")
    add_db_argument(ap)
    args = ap.parse_args()

    all_studies = args.studies or study_lists()
    done = done_studies()
    todo = [s for s in all_studies if s not in ALREADY_LOADED and s not in done]

    log(f"{len(all_studies)} studies known | {len(done)} already recorded | "
        f"{len(todo)} to process")
    if args.dry_run:
        for s in todo:
            print("  ", s)
        return

    db = connect(args.db)
    before_total = db_total_mb(args.db)
    log(f"database quantitative footprint at start: {before_total:,.1f} MB")

    for i, study in enumerate(todo, 1):
        log(f"[{i}/{len(todo)}] {study}")
        row = {k: "" for k in FIELDS}
        row.update(study=study, timestamp=datetime.now().isoformat(timespec="seconds"))
        mb_before = db_total_mb()

        ok, dt, err = run_stage("download_omics.py", study)
        row["download_s"] = round(dt, 1)
        if not ok:
            log(f"    download FAILED: {err}")
            row.update(status="download_failed", note=err)
            append_row(row)
            continue
        log(f"    downloaded in {dt:.0f}s")

        ok, dt, err = run_stage("build_omics_collections.py", study)
        row["build_s"] = round(dt, 1)
        # build_omics_collections.py logs "[OK]" and exits 0 even when a study
        # raised, so the exit code alone is not evidence: check for output.
        built = sorted((COLLECTIONS_ROOT / study).glob(f"*_samples_{study}.json")) \
            if (COLLECTIONS_ROOT / study).exists() else []
        if not ok or not built:
            reason = err or "no *_samples_*.json produced"
            log(f"    build FAILED: {reason}")
            row.update(status="build_failed", note=reason)
            append_row(row)
            if not args.no_cleanup:
                cleanup(study)
            continue
        log(f"    built in {dt:.0f}s ({len(built)} layer files)")

        # Pass the database name explicitly rather than relying on the loader's
        # own default, so both stay in step if either is pointed elsewhere.
        ok, dt, err = run_stage("load_omics_collections_to_arangodb.py", study,
                                extra=["--db", args.db])
        row["load_s"] = round(dt, 1)
        if not ok:
            log(f"    load FAILED: {err}")
            row.update(status="load_failed", note=err)
            append_row(row)
            if not args.no_cleanup:
                cleanup(study)
            continue

        counts = layer_counts(db, study)
        if sum(counts.values()) == 0:
            reason = "load reported success but no vector documents present"
            log(f"    load FAILED: {reason}")
            row.update(status="load_empty", note=reason)
            append_row(row)
            if not args.no_cleanup:
                cleanup(study)
            continue

        mb_after = db_total_mb()
        row.update(status="ok",
                   n_expression=counts["GENE_EXPRESSION"], n_cnv=counts["CNV"],
                   n_mirna=counts["MIRNA"], n_protein=counts["PROTEIN"],
                   n_methylation=counts["METHYLATION"],
                   n_samples_total=sum(counts.values()),
                   db_delta_mb=round(mb_after - mb_before, 2),
                   db_total_mb=mb_after)
        log(f"    loaded in {dt:.0f}s | {sum(counts.values())} vectors | "
            f"+{mb_after - mb_before:,.1f} MB (total {mb_after:,.1f} MB)")
        append_row(row)

        if not args.no_cleanup:
            cleanup(study)

    log(f"finished. footprint {before_total:,.1f} -> {db_total_mb():,.1f} MB")
    log(f"per-cohort measurements in {CSV_PATH}")


if __name__ == "__main__":
    main()
