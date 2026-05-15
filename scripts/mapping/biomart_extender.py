"""
biomart_extender.py
===================
Reproduce `data/mappings/biomart_gene_mappings_filled_extended.tsv` starting
from the raw BioMart harvest and a few auxiliary tables.

This script consolidates the historical four-step pipeline
(`biomart_fixer.py` -> `biomart_fixer_v2.py` -> `biomart_fixer_v2_part2.py`
plus the manual ENST merge) into a single, runnable CLI.

Pipeline
--------
Input:
    data/mappings/biomart_gene_mappings.tsv                (raw BioMart harvest)
    data/mappings/gencode.v36.annotation.gtf.gene.probemap (auto-downloaded by build_omics_collections)
    data/mappings/biomart_ensg_enst_mapping.tsv            (from biomart_maps_ensg_enst_mapping.py;
                                                            shipped as .zip)

Steps:
    1. Fill missing `hgnc_symbol` via GENCODE probemap.
    2. Fill missing `entrez_id` via self-lookup (hgnc_symbol -> entrez_id) using
       rows where both are present.
    3. (Optional, slow) Fill remaining `entrez_id` by querying NCBI Entrez for
       each orphan HGNC symbol; checkpoints to disk so the run can resume.
    4. Merge with `biomart_ensg_enst_mapping.tsv` to add `gene_id` and
       `transcript_id` columns, producing the final 15-column table.

Output:
    data/mappings/biomart_gene_mappings_filled_extended.tsv

Usage
-----
    # Quick run (steps 1, 2, 4 — no NCBI scraping; ~30s)
    python scripts/mapping/biomart_extender.py --skip-ncbi

    # Full run with NCBI Entrez fallback (recommended; takes ~30-60 min)
    python scripts/mapping/biomart_extender.py --ncbi-email you@example.com

    # Use a checkpoint so an interrupted NCBI run can resume
    python scripts/mapping/biomart_extender.py --ncbi-email you@example.com \\
        --ncbi-checkpoint data/mappings/temp/ncbi_entrez_lookup.tsv

Notes
-----
* `--ncbi-email` is required by NCBI's API terms of service for step 3.
* The version of `_extended.tsv` shipped with the repository was produced
  with this script after a multi-hour NCBI run; if you don't need full
  coverage, `--skip-ncbi` is fine for most use cases.
"""

import argparse
import re
import sys
import time
from pathlib import Path

import pandas as pd
from tqdm import tqdm

REPO_ROOT = Path(__file__).resolve().parents[2]
MAPS_DIR = REPO_ROOT / "data" / "mappings"

DEFAULT_BIOMART_RAW   = MAPS_DIR / "biomart_gene_mappings.tsv"
DEFAULT_PROBEMAP      = MAPS_DIR / "gencode.v36.annotation.gtf.gene.probemap"
DEFAULT_ENSG_ENST     = MAPS_DIR / "biomart_ensg_enst_mapping.tsv"
DEFAULT_OUTPUT        = MAPS_DIR / "biomart_gene_mappings_filled_extended.tsv"
DEFAULT_CHECKPOINT    = MAPS_DIR / "temp" / "ncbi_entrez_lookup_checkpoint.tsv"
DEFAULT_NOTFOUND      = MAPS_DIR / "temp" / "ncbi_notfound_symbols.tsv"

# Pattern matching HGNC-like symbols that contain a version suffix
# (e.g. "AC123456.1") -- these never resolve in NCBI Gene so we skip them.
_VERSIONED_SYMBOL = re.compile(r"^[A-Z][A-Z0-9]*\..*")


# ---------------------------------------------------------------------------
# Step 1: fill hgnc_symbol from GENCODE probemap
# ---------------------------------------------------------------------------

def fill_hgnc_from_gencode(biomart: pd.DataFrame, probemap_path: Path) -> pd.DataFrame:
    print(f"\n[1/4] Filling missing hgnc_symbol from {probemap_path.name} ...")
    probemap = pd.read_csv(probemap_path, sep="\t")
    probemap["id_no_version"] = probemap["id"].str.split(".").str[0]

    before = biomart["hgnc_symbol"].isna().sum()
    merged = biomart.merge(
        probemap[["id_no_version", "gene"]],
        left_on="gene_stable_id", right_on="id_no_version",
        how="left", suffixes=("", "_gencode"),
    )
    merged["hgnc_symbol"] = merged["hgnc_symbol"].fillna(merged["gene"])
    merged = merged.drop(columns=["id_no_version", "gene"])
    after = merged["hgnc_symbol"].isna().sum()
    print(f"      hgnc_symbol missing: {before:,} -> {after:,}  (filled {before - after:,})")
    return merged


# ---------------------------------------------------------------------------
# Step 2: self-lookup hgnc -> entrez (uses rows that have both)
# ---------------------------------------------------------------------------

def fill_entrez_self_lookup(biomart: pd.DataFrame) -> pd.DataFrame:
    print("\n[2/4] Filling entrez_id via self-lookup (hgnc_symbol -> entrez_id) ...")
    table = biomart[biomart["hgnc_symbol"].notna() & biomart["entrez_id"].notna()]
    lookup = dict(zip(table["hgnc_symbol"], table["entrez_id"]))
    print(f"      built lookup with {len(lookup):,} hgnc->entrez pairs")

    before = biomart["entrez_id"].isna().sum()
    mask = biomart["entrez_id"].isna() & biomart["hgnc_symbol"].notna()
    biomart.loc[mask, "entrez_id"] = biomart.loc[mask, "hgnc_symbol"].map(lookup)
    after = biomart["entrez_id"].isna().sum()
    print(f"      entrez_id missing: {before:,} -> {after:,}  (filled {before - after:,})")
    return biomart


# ---------------------------------------------------------------------------
# Step 3 (optional): NCBI Entrez API fallback
# ---------------------------------------------------------------------------

def _query_ncbi_entrez(symbol: str, max_retries: int = 3):
    """Query NCBI Gene for `symbol` (Homo sapiens). Returns the Entrez ID or None."""
    from Bio import Entrez  # local import: only required for --no-skip-ncbi
    for attempt in range(max_retries):
        try:
            handle = Entrez.esearch(
                db="gene",
                term=f"{symbol}[Gene Name] AND Homo sapiens[Organism]",
                retmax=1,
            )
            record = Entrez.read(handle)
            handle.close()
            return record["IdList"][0] if record["IdList"] else None
        except Exception as exc:
            if attempt < max_retries - 1:
                time.sleep(1)
                continue
            print(f"      [warn] giving up on {symbol}: {exc}")
            return None
    return None


def fill_entrez_via_ncbi(biomart: pd.DataFrame,
                        ncbi_email: str,
                        checkpoint_path: Path,
                        notfound_path: Path,
                        request_delay: float = 0.4,
                        checkpoint_every: int = 50) -> pd.DataFrame:
    """Fill remaining entrez_id by hitting NCBI Gene; checkpoints to disk."""
    try:
        from Bio import Entrez
    except ImportError:
        sys.exit("[ERROR] Biopython is required for the NCBI step.\n"
                 "        Install it with `pip install biopython` or pass --skip-ncbi.")
    Entrez.email = ncbi_email

    print(f"\n[3/4] Filling remaining entrez_id via NCBI Entrez API (email: {ncbi_email})")
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)

    # Resume from checkpoint
    lookup, notfound = {}, set()
    if checkpoint_path.exists():
        df_chk = pd.read_csv(checkpoint_path, sep="\t", dtype=str)
        lookup = dict(zip(df_chk["hgnc_symbol"], df_chk["entrez_id"]))
        print(f"      resumed checkpoint: {len(lookup):,} symbols already resolved")
    if notfound_path.exists():
        df_nf = pd.read_csv(notfound_path, sep="\t", dtype=str)
        notfound = set(df_nf["hgnc_symbol"].dropna().tolist())
        print(f"      resumed not-found list: {len(notfound):,} symbols")

    # Build worklist: missing entrez_id, has hgnc_symbol, not seen before
    todo_mask = biomart["entrez_id"].isna() & biomart["hgnc_symbol"].notna()
    symbols_todo = (
        biomart.loc[todo_mask, "hgnc_symbol"]
        .dropna().drop_duplicates().tolist()
    )
    symbols_todo = [s for s in symbols_todo if s not in lookup and s not in notfound]
    print(f"      symbols to query: {len(symbols_todo):,}")

    for idx, symbol in enumerate(tqdm(symbols_todo, desc="NCBI queries", ncols=80), 1):
        # Skip versioned/synthetic symbols that NCBI doesn't know about
        if _VERSIONED_SYMBOL.match(symbol):
            notfound.add(symbol)
            continue
        entrez_id = _query_ncbi_entrez(symbol)
        if entrez_id:
            lookup[symbol] = entrez_id
        else:
            notfound.add(symbol)
        # checkpoint
        if idx % checkpoint_every == 0:
            pd.DataFrame(list(lookup.items()),
                         columns=["hgnc_symbol", "entrez_id"]
                         ).to_csv(checkpoint_path, sep="\t", index=False)
            pd.DataFrame(sorted(notfound), columns=["hgnc_symbol"]
                         ).to_csv(notfound_path, sep="\t", index=False)
        # NCBI rate limit: <= 3 req/s
        time.sleep(request_delay)

    # Final checkpoint
    pd.DataFrame(list(lookup.items()), columns=["hgnc_symbol", "entrez_id"]
                 ).to_csv(checkpoint_path, sep="\t", index=False)
    pd.DataFrame(sorted(notfound), columns=["hgnc_symbol"]
                 ).to_csv(notfound_path, sep="\t", index=False)
    print(f"      NCBI lookup: {len(lookup):,} resolved, {len(notfound):,} not found")

    # Apply lookup to biomart
    before = biomart["entrez_id"].isna().sum()
    mask = biomart["entrez_id"].isna() & biomart["hgnc_symbol"].isin(lookup)
    biomart.loc[mask, "entrez_id"] = biomart.loc[mask, "hgnc_symbol"].map(lookup)
    after = biomart["entrez_id"].isna().sum()
    print(f"      entrez_id missing: {before:,} -> {after:,}  (filled {before - after:,})")
    return biomart


# ---------------------------------------------------------------------------
# Step 4: merge with ENSG -> ENST table
# ---------------------------------------------------------------------------

def attach_transcripts(biomart: pd.DataFrame, ensg_enst_path: Path) -> pd.DataFrame:
    print(f"\n[4/4] Attaching transcripts from {ensg_enst_path.name} ...")
    enst = pd.read_csv(ensg_enst_path, sep="\t", dtype=str)
    # Expected columns: gene_id, transcript_id (and maybe more)
    if "gene_id" not in enst.columns or "transcript_id" not in enst.columns:
        sys.exit(f"[ERROR] {ensg_enst_path} must contain 'gene_id' and 'transcript_id'.")

    # Aggregate transcript_id as a list per gene_id (mirrors the schema expected
    # by build_omics_collections.create_lookup_dicts -> ensg_to_ensts).
    enst_grouped = (
        enst.groupby("gene_id")["transcript_id"]
        .apply(lambda s: str(s.dropna().unique().tolist()))
        .reset_index()
    )

    merged = biomart.merge(
        enst_grouped, left_on="gene_stable_id", right_on="gene_id",
        how="left",
    )
    n_with_tx = merged["transcript_id"].notna().sum()
    print(f"      genes with at least one transcript: {n_with_tx:,} / {len(merged):,}")
    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _resolve_input(path: Path, name: str) -> Path:
    """If `path` is missing but `path.zip` exists, unzip in place; then return path."""
    if path.exists():
        return path
    zip_path = path.with_suffix(path.suffix + ".zip")
    if zip_path.exists():
        import zipfile
        print(f"      extracting {zip_path.name} ...")
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(path.parent)
        if path.exists():
            return path
    sys.exit(f"[ERROR] Required input not found: {path}\n"
             f"        (nor its .zip counterpart). "
             f"Hint: run scripts/mapping/biomart_maps{'' if name != 'ENSG->ENST' else '_ensg_enst_mapping'}.py first.")


def parse_args():
    p = argparse.ArgumentParser(
        description="Rebuild biomart_gene_mappings_filled_extended.tsv end-to-end.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--biomart-raw", type=Path, default=DEFAULT_BIOMART_RAW,
                   help="Raw BioMart harvest (default: data/mappings/biomart_gene_mappings.tsv).")
    p.add_argument("--probemap", type=Path, default=DEFAULT_PROBEMAP,
                   help="GENCODE gene probemap (auto-downloaded by build_omics_collections).")
    p.add_argument("--ensg-enst", type=Path, default=DEFAULT_ENSG_ENST,
                   help="ENSG->ENST mapping (auto-unzipped from .zip if needed).")
    p.add_argument("--output", type=Path, default=DEFAULT_OUTPUT,
                   help="Output TSV path (default: data/mappings/biomart_gene_mappings_filled_extended.tsv).")
    p.add_argument("--skip-ncbi", action="store_true",
                   help="Skip the NCBI Entrez fallback step (fast, lower coverage).")
    p.add_argument("--ncbi-email", default=None,
                   help="Email address sent to NCBI (required unless --skip-ncbi).")
    p.add_argument("--ncbi-checkpoint", type=Path, default=DEFAULT_CHECKPOINT,
                   help="Checkpoint TSV for resuming an interrupted NCBI run.")
    p.add_argument("--ncbi-notfound", type=Path, default=DEFAULT_NOTFOUND,
                   help="TSV listing symbols NCBI couldn't resolve (skipped on re-runs).")
    p.add_argument("--ncbi-delay", type=float, default=0.4,
                   help="Seconds between NCBI requests (NCBI limit: 3/s; default: 0.4).")
    return p.parse_args()


def main():
    args = parse_args()
    if not args.skip_ncbi and not args.ncbi_email:
        sys.exit("[ERROR] --ncbi-email is required for the NCBI step. "
                 "Pass --skip-ncbi if you want to skip it.")

    print("=" * 60)
    print("biomart_extender")
    print(f"  input (raw)   : {args.biomart_raw}")
    print(f"  probemap      : {args.probemap}")
    print(f"  ENSG->ENST    : {args.ensg_enst}")
    print(f"  output        : {args.output}")
    print(f"  ncbi step     : {'skipped' if args.skip_ncbi else 'enabled'}")
    print("=" * 60)

    biomart_raw  = _resolve_input(args.biomart_raw, "BioMart raw")
    probemap     = _resolve_input(args.probemap, "GENCODE probemap")
    ensg_enst    = _resolve_input(args.ensg_enst, "ENSG->ENST")

    print(f"\nReading {biomart_raw.name} ...")
    biomart = pd.read_csv(biomart_raw, sep="\t")
    print(f"  loaded {len(biomart):,} rows, columns: {list(biomart.columns)}")

    biomart = fill_hgnc_from_gencode(biomart, probemap)
    biomart = fill_entrez_self_lookup(biomart)
    if not args.skip_ncbi:
        biomart = fill_entrez_via_ncbi(
            biomart,
            ncbi_email=args.ncbi_email,
            checkpoint_path=args.ncbi_checkpoint,
            notfound_path=args.ncbi_notfound,
            request_delay=args.ncbi_delay,
        )
    biomart = attach_transcripts(biomart, ensg_enst)

    # Reorder columns to match the canonical shipped schema
    canonical_order = [
        "gene_stable_id", "gene_stable_id_version", "gene_type", "gene_description",
        "chromosome", "gene_start_bp", "gene_end_bp", "strand",
        "hgnc_symbol", "mirbase_id", "uniprot_swissprot_id", "uniprot_isoform_id",
        "entrez_id", "gene_id", "transcript_id",
    ]
    present = [c for c in canonical_order if c in biomart.columns]
    extras  = [c for c in biomart.columns if c not in canonical_order]
    biomart = biomart[present + extras]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    biomart.to_csv(args.output, sep="\t", index=False)
    print("\n" + "=" * 60)
    print(f"[OK] wrote {args.output}")
    print(f"     rows: {len(biomart):,}  columns: {len(biomart.columns)}")
    print("=" * 60)


if __name__ == "__main__":
    main()
