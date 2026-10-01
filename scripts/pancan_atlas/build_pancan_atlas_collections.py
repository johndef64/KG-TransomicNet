"""
build_pancan_atlas_collections.py
=================================
Build the ArangoDB-ready JSON collections of the PanCanAtlas view (view P)
from the matrices fetched by scripts/pancan_atlas/download_pancan_atlas.py.

View P is a second, pan-cancer view of the quantitative layer. It sits next
to the per-study GDC view and never modifies it:

  * it is written to its own folder, data/arangodb_collections/PANCAN-ATLAS/;
  * every document carries cohort = "PANCAN-ATLAS", so the cohort-keyed
    queries of query_utils.py address it as one more cohort
    (index keys expr_index_PANCAN-ATLAS, cnv_index_PANCAN-ATLAS, ...);
  * vector keys are 15-character PanCanAtlas barcodes (TCGA-XX-XXXX-01),
    while GDC vectors use 16-character barcodes (TCGA-XX-XXXX-01A), so the
    two views can share the unified collections without key collisions.

Output files (names follow LAYER_SPEC in load_omics_collections_to_arangodb.py,
so the existing loader reads them with --studies PANCAN-ATLAS):

  gene_expression_index.json   gene_expression_samples_PANCAN-ATLAS.json
  cnv_index.json               cnv_samples_PANCAN-ATLAS.json
  mirna_index.json             mirna_samples_PANCAN-ATLAS.json
  samples.json                 (one document per sample, iCluster label)
  projects.json                (one document, the PANCAN-ATLAS cohort)

Gene identifiers. Expression and CNV rows are HGNC symbols (plus a few bare
Entrez IDs). They are resolved to Entrez with the HGNC conversion table, using
approved symbols only: the same rule used by the downstream classifiers, so the
gene vocabulary of view P is the one their feature matrices are built on.

Value fields.
  expression  values_ebpp          (log2(norm_value+1), EB++ adjusted)
  cnv         values_copy_number   (GISTIC2 gene-level, continuous; see
                                    analysis_method / value_scale)
  mirna       values_expression    (EB-adjusted, mature miRNAs)
Expression is read with value_type="ebpp"; CNV and miRNA use the field names
already queried by query_utils.py.

Usage
-----
    # Full build (about 5-6 GB of JSON; run it on the build machine)
    python scripts/pancan_atlas/build_pancan_atlas_collections.py

    # Quick test on the first 50 samples of each matrix, then check the output
    python scripts/pancan_atlas/build_pancan_atlas_collections.py --limit-samples 50 --verify

    # Only some layers
    python scripts/pancan_atlas/build_pancan_atlas_collections.py --layers gene_expression mirna
"""

import argparse
import json
import logging
import math
import random
import re
import zipfile
from pathlib import Path
from typing import Dict, Iterable, Optional

import sys

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
SCRIPTS = SCRIPT_DIR.parent
PROJECT_ROOT = SCRIPTS.parent
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))          # reuse the GDC builder of scripts/

from build_omics_collections import CollectionBuilder  # noqa: E402
RAW_DIR = PROJECT_ROOT / "data" / "omics" / "PANCAN-ATLAS"
MAPS_ROOT = PROJECT_ROOT / "data" / "mappings"
DEFAULT_OUTPUT_ROOT = PROJECT_ROOT / "data" / "arangodb_collections"
HGNC_TABLE = MAPS_ROOT / "hgcn_mart_conversion_table.zip"
MIRNA_MAP = MAPS_ROOT / "mirna_hgcn_map.tsv"

COHORT = "PANCAN-ATLAS"
SOURCE_RELEASE = "PanCanAtlas 2018"

logger = logging.getLogger(__name__)

# layer -> raw file, decimals kept from the source, collection file stem
LAYERS = {
    "gene_expression": {"raw": "expression.tsv.gz", "decimals": 2, "stem": "gene_expression"},
    "cnv":             {"raw": "cnv.tsv.gz",        "decimals": 3, "stem": "cnv"},
    "mirna":           {"raw": "mirna.tsv.gz",      "decimals": 2, "stem": "mirna"},
}


# ---------------------------------------------------------------------------
# Inputs
# ---------------------------------------------------------------------------

def read_matrix(path: Path, limit_samples: Optional[int] = None) -> pd.DataFrame:
    """Feature x sample matrix, float32, features as index (source order kept)."""
    header = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    cols = header if limit_samples is None else header[: limit_samples + 1]
    dtypes = {c: np.float32 for c in cols[1:]}
    dtypes[cols[0]] = str
    df = pd.read_csv(path, sep="\t", usecols=cols, dtype=dtypes, index_col=0,
                     na_values=["NA", "NaN", ""], low_memory=False)
    logger.info(f"  {path.name}: {df.shape[0]:,} features x {df.shape[1]:,} samples")
    return df


def load_hgnc_table(path: Path = HGNC_TABLE) -> Dict[str, Dict[str, str]]:
    """Approved symbol -> {entrez_id, ensembl_id, hgnc_id} (unique symbols)."""
    conv = pd.read_csv(path, sep="\t", compression="zip", dtype=str)
    conv["NCBI gene ID"] = pd.to_numeric(conv["NCBI gene ID"], errors="coerce").astype("Int64")
    conv = conv.dropna(subset=["NCBI gene ID", "Approved symbol"]).drop_duplicates(subset=["Approved symbol"])
    return {
        row["Approved symbol"]: {
            "entrez_id":  str(int(row["NCBI gene ID"])),
            "ensembl_id": row["Ensembl gene ID"] if pd.notna(row["Ensembl gene ID"]) else None,
            "hgnc_id":    row["HGNC ID"],
        }
        for _, row in conv.iterrows()
    }


_ARM_RE = re.compile(r"-(3p|5p)$")


def mature_to_precursor(mature_ids: Iterable[str], path: Path = MIRNA_MAP) -> Dict[str, str]:
    """Mature miRNA -> HGNC symbol of its precursor, when the match is unique.

    The arm suffix is dropped and miR is lowered to mir (hsa-miR-21-5p ->
    hsa-mir-21); the result must be a precursor id of the miRNA map. Mature
    forms shared by several precursors (hsa-let-7a-5p) stay unresolved.
    """
    mp = pd.read_csv(path, sep="\t", dtype=str)
    pre2hgnc = dict(zip(mp["miRNA_ID"], mp["hgcn_id"]))
    out = {}
    for m in mature_ids:
        pre = _ARM_RE.sub("", m).replace("-miR-", "-mir-")
        if pre in pre2hgnc and pd.notna(pre2hgnc[pre]):
            out[m] = pre2hgnc[pre]
    return out


# ---------------------------------------------------------------------------
# Builders
# ---------------------------------------------------------------------------

class PanCanAtlasViewBuilder(CollectionBuilder):
    """Index and vector documents of view P, written as streamed JSON lines."""

    def __init__(self, output_dir: Path, hgnc: Dict[str, Dict[str, str]]):
        super().__init__(COHORT, mapping_dfs={}, lookups={}, output_dir=output_dir)
        self.hgnc = hgnc

    def save_stream(self, documents: Iterable[Dict], collection_name: str) -> int:
        """Like save_collection, but consumes a generator (vectors do not fit in RAM as dicts)."""
        out = self.output_dir / f"{collection_name}.json"
        n = 0
        with open(out, "w", encoding="utf-8") as fh:
            for doc in documents:
                json.dump(doc, fh, ensure_ascii=False, default=self._json_default)
                fh.write("\n")
                n += 1
        logger.info(f"  -> wrote {n:,} documents to {out.name}")
        return n

    # -- index documents ---------------------------------------------------

    def _gene_mappings(self, row_ids) -> list:
        maps = []
        for i, rid in enumerate(row_ids):
            hit = self.hgnc.get(rid)
            maps.append({
                "position":      i,
                "source_id":     rid,
                "hgnc_symbol":   rid if hit else None,
                "entrez_id":     hit["entrez_id"] if hit else None,
                "gene_id_base":  hit["ensembl_id"] if hit else None,
                "gene_ref":      f"genes/{hit['ensembl_id']}" if hit and hit["ensembl_id"] else None,
            })
        return maps

    def expression_index(self, df: pd.DataFrame) -> Dict:
        maps = self._gene_mappings(df.index.tolist())
        return {
            "_key": f"expr_index_{COHORT}", "cohort": COHORT,
            "data_type": "gene_expression_index", "n_genes": len(maps),
            "platform": "IlluminaHiSeq RNASeqV2", "normalization": "EB++ batch-adjusted",
            "unit": "log2(norm_value+1)", "source_release": SOURCE_RELEASE,
            "mapping_rule": "HGNC approved symbol -> NCBI gene ID",
            "gene_mappings": maps,
            "description": "Gene position mapping for PanCanAtlas expression vectors.",
        }

    def cnv_index(self, df: pd.DataFrame) -> Dict:
        maps = self._gene_mappings(df.index.tolist())
        return {
            "_key": f"cnv_index_{COHORT}", "cohort": COHORT,
            "data_type": "cnv_index", "n_genes": len(maps),
            "analysis_method": "GISTIC2", "value_scale": "gene-level, continuous (all_data_by_genes)",
            "source_release": SOURCE_RELEASE,
            "mapping_rule": "HGNC approved symbol -> NCBI gene ID",
            "gene_mappings": maps,
            "description": "Gene position mapping for PanCanAtlas CNV vectors.",
        }

    def mirna_index(self, df: pd.DataFrame) -> Dict:
        ids = df.index.tolist()
        pre = mature_to_precursor(ids)
        maps = []
        for i, mid in enumerate(ids):
            sym = pre.get(mid)
            hit = self.hgnc.get(sym) if sym else None
            maps.append({
                "position":    i,
                "mirna_id":    mid,
                "mirbase_id":  mid,
                "form":        "mature",
                "hgnc_symbol": sym,
                "entrez_id":   hit["entrez_id"] if hit else None,
                "description": f"miRNA {mid}",
            })
        return {
            "_key": f"mirna_index_{COHORT}", "cohort": COHORT,
            "data_type": "mirna_index", "n_mirnas": len(maps),
            "platform": "Illumina miRNA-seq", "normalization": "EB batch-adjusted",
            "source_release": SOURCE_RELEASE,
            "mapping_rule": "mature -> unique precursor (mirna_hgcn_map) -> NCBI gene ID",
            "mirna_mappings": maps,
            "description": "miRNA position mapping for PanCanAtlas miRNA vectors.",
        }

    # -- vector documents --------------------------------------------------

    @staticmethod
    def _column(df: pd.DataFrame, sid: str, decimals: int) -> list:
        col = np.round(df[sid].to_numpy(dtype=np.float64), decimals)
        values = col.tolist()
        if np.isnan(col).any():
            values = [None if (v is None or math.isnan(v)) else v for v in values]
        return values

    def vectors(self, layer: str, df: pd.DataFrame) -> Iterable[Dict]:
        dec = LAYERS[layer]["decimals"]
        for sid in df.columns:
            base = {"_key": sid, "sample_id": sid, "cohort": COHORT, "source_release": SOURCE_RELEASE}
            if layer == "gene_expression":
                yield {**base, "data_type": "gene_expression_vector",
                       "expression_index_ref": f"expr_index_{COHORT}",
                       "platform": "IlluminaHiSeq RNASeqV2", "n_genes": len(df),
                       "values_ebpp": self._column(df, sid, dec),
                       "normalization": "EB++ batch-adjusted, log2(norm_value+1)"}
            elif layer == "cnv":
                yield {**base, "data_type": "cnv_vector",
                       "cnv_index_ref": f"cnv_index_{COHORT}",
                       "analysis_method": "GISTIC2", "n_genes": len(df),
                       "values_copy_number": self._column(df, sid, dec)}
            elif layer == "mirna":
                yield {**base, "data_type": "mirna_vector",
                       "mirna_index_ref": f"mirna_index_{COHORT}",
                       "platform": "Illumina miRNA-seq", "n_mirnas": len(df),
                       "values_expression": self._column(df, sid, dec)}

    # -- metadata ----------------------------------------------------------

    def samples(self, sample_sets: Dict[str, set], icluster: pd.Series) -> list:
        all_ids = sorted(set().union(*sample_sets.values()) | set(icluster.index))
        docs = []
        for sid in all_ids:
            label = icluster.get(sid)
            docs.append({
                "_key": sid, "sample_id": sid, "submitter_id": sid[:12],
                "case_ref": f"cases/{sid[:12]}", "project_ref": f"projects/{COHORT}",
                "sample_type_code": sid[13:15], "cohort": COHORT,
                "layers": [l for l, s in sample_sets.items() if sid in s],
                "icluster_k28": label if isinstance(label, str) else None,
                "source_release": SOURCE_RELEASE, "entity_type": "sample",
            })
        return docs

    def project(self, n_samples: int) -> list:
        return [{
            "_key": COHORT, "name": "TCGA PanCanAtlas (harmonised pan-cancer view)",
            "program": "TCGA", "source_release": SOURCE_RELEASE,
            "entity_type": "project", "n_samples": n_samples,
        }]


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def verify_layer(out_dir: Path, layer: str, df: pd.DataFrame, n_check: int = 5) -> bool:
    """Re-read some written vectors and compare them with the source matrix."""
    field = {"gene_expression": "values_ebpp", "cnv": "values_copy_number",
             "mirna": "values_expression"}[layer]
    path = out_dir / f"{LAYERS[layer]['stem']}_samples_{COHORT}.json"
    wanted = set(random.Random(0).sample(list(df.columns), min(n_check, df.shape[1])))
    tol = 0.5 * 10 ** (-LAYERS[layer]["decimals"]) + 1e-6
    ok = True
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            doc = json.loads(line)
            if doc["_key"] not in wanted:
                continue
            got = np.array([np.nan if v is None else v for v in doc[field]], dtype=np.float64)
            ref = df[doc["_key"]].to_numpy(dtype=np.float64)
            same_nan = np.array_equal(np.isnan(got), np.isnan(ref))
            diff = np.nanmax(np.abs(got - ref)) if got.size else 0.0
            if not (same_nan and diff <= tol and len(got) == len(ref)):
                ok = False
                logger.error(f"  [verify] {layer} {doc['_key']}: max diff {diff:.4g}, NaN pattern equal {same_nan}")
    logger.info(f"  [verify] {layer}: {'OK' if ok else 'FAILED'} on {len(wanted)} samples")
    return ok


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--raw-dir", type=Path, default=RAW_DIR)
    p.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT_ROOT)
    p.add_argument("--layers", nargs="+", choices=list(LAYERS), default=list(LAYERS))
    p.add_argument("--limit-samples", type=int, default=None,
                   help="Use only the first N samples of each matrix (test build)")
    p.add_argument("--verify", action="store_true",
                   help="Re-read written vectors and compare them with the source")
    return p.parse_args()


def main():
    logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
    args = parse_args()
    out_dir = args.output_root / COHORT
    if not HGNC_TABLE.exists():
        raise SystemExit(f"HGNC conversion table not found: {HGNC_TABLE}")
    builder = PanCanAtlasViewBuilder(out_dir, load_hgnc_table())
    logger.info(f"Building view P ({SOURCE_RELEASE}) -> {out_dir}")

    sample_sets, all_ok = {}, True
    for layer in args.layers:
        raw = args.raw_dir / LAYERS[layer]["raw"]
        if not raw.exists():
            logger.warning(f"[{layer}] raw file missing: {raw} (run download_pancan_atlas.py)")
            continue
        logger.info(f"[{layer}] reading {raw.name}")
        df = read_matrix(raw, args.limit_samples)
        index_doc = {"gene_expression": builder.expression_index,
                     "cnv": builder.cnv_index, "mirna": builder.mirna_index}[layer](df)
        builder.save_collection([index_doc], f"{LAYERS[layer]['stem']}_index")
        builder.save_stream(builder.vectors(layer, df), f"{LAYERS[layer]['stem']}_samples_{COHORT}")
        sample_sets[layer] = set(df.columns)
        if args.verify:
            all_ok &= verify_layer(out_dir, layer, df)
        del df

    clinical = args.raw_dir / "clinical.tsv.gz"
    icluster = pd.Series(dtype=str)
    if clinical.exists():
        cl = pd.read_csv(clinical, sep="\t", dtype=str)
        icluster = cl.set_index(cl.columns[0])[cl.columns[1]]
        if args.limit_samples is not None:
            keep = set().union(*sample_sets.values()) if sample_sets else set()
            icluster = icluster[icluster.index.isin(keep)]
    samples = builder.samples(sample_sets, icluster) if sample_sets else []
    builder.save_collection(samples, "samples")
    builder.save_collection(builder.project(len(samples)), "projects")

    if args.verify and not all_ok:
        raise SystemExit("Verification failed")


if __name__ == "__main__":
    main()
