"""
load_omics_collections_to_arangodb.py
=====================================
Load the per-study omic JSON collections produced by
scripts/build_omics_collections.py into ArangoDB.

Each requested study is read from data/arangodb_collections/<STUDY>/ and pushed
into the following unified collections (one shared schema across all studies):

  Semantic / metadata layer (loaded once, idempotent via _key):
    GENES, PROJECTS, SAMPLES, CASES

  Unified quantitative collections (index + sample documents combined,
  distinguished by the `data_type` field):
    GENE_EXPRESSION, CNV, MIRNA, PROTEIN, METHYLATION

Loading is incremental: existing documents are kept unless --replace-existing
is passed (which drops the affected collections before reload).

Usage examples (mirror the CLI of scripts/download_omics.py)
------------------------------------------------------------
    # Default: TCGA-BRCA, all layers, into db PKT_main
    python scripts/load_omics_collections_to_arangodb.py

    # Specific studies
    python scripts/load_omics_collections_to_arangodb.py --studies TCGA-BRCA TCGA-LUAD

    # Full cohort
    python scripts/load_omics_collections_to_arangodb.py --cohort tcga
    python scripts/load_omics_collections_to_arangodb.py --cohort target
    python scripts/load_omics_collections_to_arangodb.py --cohort all

    # Restrict to specific layers
    python scripts/load_omics_collections_to_arangodb.py --layers gene_expression mirna

    # Drop and reload everything (DESTRUCTIVE)
    python scripts/load_omics_collections_to_arangodb.py --replace-existing

    # Custom ArangoDB endpoint / credentials
    python scripts/load_omics_collections_to_arangodb.py \\
        --host http://my-server:8529 --user myuser --password mypass --db MyDB

Default connection settings are defined at the top of this file
(ARANGODB_HOSTS / ARANGODB_USER / ARANGODB_PASSWORD).
"""

import argparse
import json
import logging
import math
import os
import sys
from pathlib import Path
from typing import List

from arango import ArangoClient

# Reuse the study catalogue from download_omics
from download_omics import TCGA_STUDIES, TARGET_STUDIES

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parent
DEFAULT_COLLECTIONS_ROOT = PROJECT_ROOT / "data" / "arangodb_collections"
DEFAULT_DB_NAME = "PKT_main"

# Default ArangoDB connection (kept here for backwards compatibility; can be
# overridden via the CLI or by configuring scripts/arangodb_utils.py).
ARANGODB_HOSTS = "http://localhost:8529"
ARANGODB_USER = "root"
ARANGODB_PASSWORD = "avocadodb"

BATCH_SIZE = 1000

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)

# Map between layer name (CLI) and (unified collection, index file, samples file template).
LAYER_SPEC = {
    "gene_expression": ("GENE_EXPRESSION", "gene_expression_index.json", "gene_expression_samples_{study}.json"),
    "cnv":             ("CNV",             "cnv_index.json",             "cnv_samples_{study}.json"),
    "mirna":           ("MIRNA",           "mirna_index.json",           "mirna_samples_{study}.json"),
    "protein":         ("PROTEIN",         "protein_index.json",         "protein_samples_{study}.json"),
    "methylation":     ("METHYLATION",     "methylation_index.json",     "methylation_samples_{study}.json"),
}
ALL_LAYERS = list(LAYER_SPEC.keys())

# Metadata files always loaded when present
METADATA_FILES = {
    "projects.json": "PROJECTS",
    "samples.json":  "SAMPLES",
    "cases.json":    "CASES",
}


# ---------------------------------------------------------------------------
# Connection / collection management
# ---------------------------------------------------------------------------

def setup_arangodb_connection(db_name: str,
                              host: str = ARANGODB_HOSTS,
                              user: str = ARANGODB_USER,
                              password: str = ARANGODB_PASSWORD):
    """Connect to ArangoDB and ensure the target database exists."""
    try:
        client = ArangoClient(hosts=host)
        sys_db = client.db("_system", username=user, password=password)
        if not sys_db.has_database(db_name):
            sys_db.create_database(db_name)
            logger.info(f"Created database: {db_name}")
        else:
            logger.info(f"Database {db_name} already exists")
        db = client.db(db_name, username=user, password=password)
        logger.info(f"Connected to {host} / database '{db_name}'")
        return db
    except Exception as exc:
        logger.error(f"Failed to connect to ArangoDB at {host}: {exc}")
        return None


def create_tcga_collections(db):
    """Create the unified omic collections + semantic/metadata collections and indexes."""
    omic_collections = [spec[0] for spec in LAYER_SPEC.values()]
    other_collections = ["GENES", "PROJECTS", "SAMPLES", "CASES"]

    for cname in omic_collections + other_collections:
        if not db.has_collection(cname):
            db.create_collection(cname)
            logger.info(f"Created collection: {cname}")
        else:
            logger.info(f"Collection {cname} already exists ({db.collection(cname).count()} docs)")

    try:
        db.collection("GENES").add_index(
            {"type": "hash", "fields": ["symbol"], "unique": False, "sparse": True}
        )
        db.collection("SAMPLES").add_index(
            {"type": "hash", "fields": ["submitter_id"], "unique": False, "sparse": True}
        )
        for cname in omic_collections:
            col = db.collection(cname)
            col.add_index({"type": "hash", "fields": ["data_type"], "unique": False, "sparse": False})
            col.add_index({"type": "hash", "fields": ["sample_id"], "unique": False, "sparse": True})
            col.add_index({"type": "hash", "fields": ["cohort"],    "unique": False, "sparse": True})
        logger.info("Created indexes on omic / metadata collections")
    except Exception as exc:
        logger.warning(f"Could not create some indexes: {exc}")


# ---------------------------------------------------------------------------
# JSON loading / sanitisation
# ---------------------------------------------------------------------------

def _sanitize_value(val):
    """Recursively replace NaN/Inf with None so ArangoDB accepts the documents."""
    if isinstance(val, float):
        if math.isnan(val) or math.isinf(val):
            return None
        return val
    if isinstance(val, list):
        return [_sanitize_value(v) for v in val]
    if isinstance(val, dict):
        return {k: _sanitize_value(v) for k, v in val.items()}
    return val


def _load_json_lines(path: str):
    """Load a JSON file in either array or JSON Lines format and sanitise NaN/Inf."""
    with open(path, "r", encoding="utf-8") as fh:
        first = fh.read(1)
        fh.seek(0)
        if not first:
            return []
        if first == "[":
            docs = json.load(fh)
        else:
            docs = [json.loads(line) for line in fh if line.strip()]
    return [_sanitize_value(d) for d in docs]


def _insert_incremental(collection, docs):
    """Insert documents in batches, counting inserts / skipped duplicates / failures."""
    total = len(docs)
    if total == 0:
        logger.info("  no documents to insert")
        return

    inserted = skipped = failed = 0
    sample_errors = []

    for i in range(0, total, BATCH_SIZE):
        batch = docs[i:i + BATCH_SIZE]
        try:
            result = collection.insert_many(batch, overwrite=False, silent=False)
            for r in result:
                if isinstance(r, dict):
                    if r.get("_key") and not r.get("error"):
                        inserted += 1
                    elif r.get("errorNum", 0) == 1210:  # unique constraint violated
                        skipped += 1
                    else:
                        failed += 1
                        if len(sample_errors) < 3:
                            sample_errors.append(r)
                elif hasattr(r, "error_code"):
                    if r.error_code == 1210:
                        skipped += 1
                    else:
                        failed += 1
                        if len(sample_errors) < 3:
                            sample_errors.append(str(r))
                else:
                    inserted += 1
            if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= total:
                logger.info(f"  processed {min(i + BATCH_SIZE, total)}/{total}")
        except Exception as exc:
            failed += len(batch)
            logger.warning(f"  batch {i // BATCH_SIZE + 1} failed: {exc}")

    logger.info(f"  -> inserted {inserted}, skipped(duplicate) {skipped}, failed {failed}")
    if sample_errors:
        logger.warning(f"  sample errors: {sample_errors[:3]}")


# ---------------------------------------------------------------------------
# Per-study loader
# ---------------------------------------------------------------------------

def load_study(db, study: str, layers: List[str], collections_root: Path,
               do_semantic: bool, do_metadata: bool, replace_existing: bool):
    """Load all artifacts for one study from collections_root/<study>/."""
    study_dir = collections_root / study
    if not study_dir.exists():
        logger.warning(f"[{study}] collection directory not found: {study_dir} -- skipping")
        return

    logger.info("=" * 60)
    logger.info(f"Loading {study}  from  {study_dir}")
    logger.info("=" * 60)

    if replace_existing:
        logger.warning("--replace-existing is set: dropping affected collections first")
        affected = ["GENES"] + [LAYER_SPEC[l][0] for l in layers] + ["PROJECTS", "SAMPLES", "CASES"]
        for cname in affected:
            if db.has_collection(cname):
                db.delete_collection(cname)
                logger.info(f"  dropped {cname}")
        create_tcga_collections(db)

    # 1. GENES (semantic layer)
    if do_semantic:
        genes_path = study_dir / "genes.json"
        if genes_path.exists():
            logger.info(f"[{study}] Loading GENES from {genes_path.name}")
            _insert_incremental(db.collection("GENES"), _load_json_lines(str(genes_path)))
        else:
            logger.info(f"[{study}] genes.json not found, skipping")

    # 2. Unified omic collections (index + samples)
    for layer in layers:
        cname, index_fname, samples_template = LAYER_SPEC[layer]
        col = db.collection(cname)

        index_path = study_dir / index_fname
        if index_path.exists():
            logger.info(f"[{study}] Loading {cname} index from {index_fname}")
            _insert_incremental(col, _load_json_lines(str(index_path)))
        else:
            logger.info(f"[{study}] {index_fname} not found, skipping index")

        samples_path = study_dir / samples_template.format(study=study)
        if samples_path.exists():
            logger.info(f"[{study}] Loading {cname} samples from {samples_path.name}")
            _insert_incremental(col, _load_json_lines(str(samples_path)))
        else:
            logger.info(f"[{study}] {samples_path.name} not found, skipping samples")

    # 3. Metadata
    if do_metadata:
        for fname, cname in METADATA_FILES.items():
            path = study_dir / fname
            if not path.exists():
                logger.info(f"[{study}] {fname} not found, skipping")
                continue
            logger.info(f"[{study}] Loading {cname} from {fname}")
            _insert_incremental(db.collection(cname), _load_json_lines(str(path)))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Load per-study omic JSON collections into ArangoDB.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Usage")[1] if "Usage" in __doc__ else "",
    )
    group = p.add_mutually_exclusive_group()
    group.add_argument("--studies", nargs="+", metavar="STUDY",
                       help="One or more study IDs to load (e.g. TCGA-BRCA TARGET-AML).")
    group.add_argument("--cohort", choices=["tcga", "target", "all"],
                       help="Load an entire cohort. 'all' = TCGA + TARGET.")
    p.add_argument("--layers", nargs="+", choices=ALL_LAYERS, default=ALL_LAYERS,
                   help=f"Quantitative layers to load (default: all of {ALL_LAYERS}).")
    p.add_argument("--db", default=DEFAULT_DB_NAME,
                   help=f"ArangoDB database name (default: {DEFAULT_DB_NAME}).")
    p.add_argument("--host", default=ARANGODB_HOSTS,
                   help=f"ArangoDB host URL (default: {ARANGODB_HOSTS}).")
    p.add_argument("--user", default=ARANGODB_USER,
                   help=f"ArangoDB user (default: {ARANGODB_USER}).")
    p.add_argument("--password", default=ARANGODB_PASSWORD,
                   help="ArangoDB password (default: same as the project default).")
    p.add_argument("--collections-root", type=Path, default=DEFAULT_COLLECTIONS_ROOT,
                   help="Root directory holding the per-study JSON collections.")
    p.add_argument("--no-semantic", action="store_true",
                   help="Skip loading the semantic layer (GENES).")
    p.add_argument("--no-metadata", action="store_true",
                   help="Skip loading the metadata layer (PROJECTS, SAMPLES, CASES).")
    p.add_argument("--replace-existing", action="store_true",
                   help="Drop affected collections before loading (DESTRUCTIVE).")
    return p.parse_args()


def resolve_studies(args) -> List[str]:
    if args.studies:
        return args.studies
    if args.cohort == "tcga":   return TCGA_STUDIES
    if args.cohort == "target": return TARGET_STUDIES
    if args.cohort == "all":    return TCGA_STUDIES + TARGET_STUDIES
    return ["TCGA-BRCA"]


def main():
    args = parse_args()
    studies = resolve_studies(args)
    if not args.studies and not args.cohort:
        logger.info("No --studies or --cohort provided; defaulting to TCGA-BRCA.")

    logger.info("=" * 60)
    logger.info("load_omics_collections_to_arangodb")
    logger.info(f"  host             : {args.host}")
    logger.info(f"  user             : {args.user}")
    logger.info(f"  database         : {args.db}")
    logger.info(f"  collections root : {args.collections_root}")
    logger.info(f"  studies          : {', '.join(studies)}")
    logger.info(f"  layers           : {', '.join(args.layers)}")
    logger.info(f"  replace_existing : {args.replace_existing}")
    logger.info("=" * 60)

    db = setup_arangodb_connection(args.db, host=args.host, user=args.user, password=args.password)
    if db is None:
        sys.exit("[ERROR] Could not connect to ArangoDB.")

    create_tcga_collections(db)

    for study in studies:
        try:
            load_study(
                db, study, args.layers, args.collections_root,
                do_semantic=not args.no_semantic,
                do_metadata=not args.no_metadata,
                replace_existing=args.replace_existing,
            )
            # --replace-existing only triggers a drop on the first study; subsequent
            # studies share the same target collections.
            args.replace_existing = False
        except Exception as exc:
            logger.exception(f"[{study}] load failed: {exc}")

    # Final summary
    logger.info("\n--- Collections in database ---")
    for col_info in db.collections():
        if col_info["name"].startswith("_"):
            continue
        cnt = db.collection(col_info["name"]).count()
        logger.info(f"  {col_info['name']:<20} {cnt:>12} documents")

    logger.info("[OK] load_omics_collections_to_arangodb finished.")


if __name__ == "__main__":
    main()
