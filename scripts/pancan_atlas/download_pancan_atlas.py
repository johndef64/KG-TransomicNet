"""
download_pancan_atlas.py
========================
Download the TCGA PanCanAtlas omic matrices that form the *PanCanAtlas view*
(view P) of the quantitative layer, next to the per-study GDC view built by
scripts/download_omics.py.

View P is an addition: it does not touch any GDC file or collection. Its raw
inputs are the four Xena matrices of the PanCanAtlas release (2018):

  expression  EB++ batch-adjusted RNA-seq (IlluminaHiSeq RNASeqV2), gene-level
  mirna       EB-adjusted miRNA-seq, mature miRNAs
  cnv         GISTIC2 gene-level copy number (all_data_by_genes)
  clinical    iCluster k=28 tumour-sample assignments (subtype labels)

Files are stored in data/omics/PANCAN-ATLAS/ as <key>.tsv.gz, together with
a manifest (manifest.json) that records URL, size and SHA-256 of each file.

Usage
-----
    # Download all four files
    python scripts/pancan_atlas/download_pancan_atlas.py

    # Only some of them
    python scripts/pancan_atlas/download_pancan_atlas.py --keys expression,clinical

    # Reuse a copy already on disk (same file names), then write the manifest
    python scripts/pancan_atlas/download_pancan_atlas.py --from-dir D:/path/to/tcga_pancan
"""

import argparse
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

import requests

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[1]
DEFAULT_OUT = PROJECT_ROOT / "data" / "omics" / "PANCAN-ATLAS"

PANCAN_ATLAS_HUB = "https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/"
TCGA_XENA_HUB = "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/"

SOURCE_RELEASE = "PanCanAtlas 2018"

# key -> (hub, Xena file name). The key is also the local file stem.
PANCAN_ATLAS_FILES = {
    "expression": (PANCAN_ATLAS_HUB, "EB%2B%2BAdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz"),
    "mirna":      (PANCAN_ATLAS_HUB, "pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz"),
    "cnv":        (TCGA_XENA_HUB,    "TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz"),
    "clinical":   (PANCAN_ATLAS_HUB, "TCGA_PanCan33_iCluster_k28_tumor.gz"),
}


def _sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def _download(url: str, dest: Path) -> None:
    tmp = dest.with_suffix(dest.suffix + ".part")
    with requests.get(url, stream=True, timeout=60) as r:
        r.raise_for_status()
        with open(tmp, "wb") as fh:
            for block in r.iter_content(chunk_size=1 << 20):
                fh.write(block)
    tmp.replace(dest)


def fetch(out_dir: Path, keys=None, from_dir: Path = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest = {"source_release": SOURCE_RELEASE, "files": {}}
    for key, (hub, fname) in PANCAN_ATLAS_FILES.items():
        if keys and key not in keys:
            continue
        url = hub + fname
        dest = out_dir / f"{key}.tsv.gz"
        if dest.exists():
            print(f"[{key}] already present: {dest}")
        elif from_dir is not None and (from_dir / dest.name).exists():
            print(f"[{key}] copying from {from_dir / dest.name}")
            shutil.copy2(from_dir / dest.name, dest)
        else:
            print(f"[{key}] downloading {url}")
            _download(url, dest)
        manifest["files"][key] = {
            "url": url,
            "local_file": dest.name,
            "size_bytes": dest.stat().st_size,
            "sha256": _sha256(dest),
        }
        print(f"[{key}] {dest.stat().st_size:,} bytes")
    manifest["written_at"] = datetime.now(timezone.utc).isoformat(timespec="seconds")

    manifest_path = out_dir / "manifest.json"
    if manifest_path.exists():
        old = json.loads(manifest_path.read_text(encoding="utf-8"))
        old.get("files", {}).update(manifest["files"])
        manifest["files"] = old["files"]
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"Manifest -> {manifest_path}")
    return manifest


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Output directory")
    p.add_argument("--keys", default=None,
                   help=f"Comma-separated subset of {','.join(PANCAN_ATLAS_FILES)}")
    p.add_argument("--from-dir", type=Path, default=None,
                   help="Directory holding <key>.tsv.gz copies to reuse instead of downloading")
    return p.parse_args()


def main():
    args = parse_args()
    keys = args.keys.split(",") if args.keys else None
    fetch(args.out, keys, args.from_dir)


if __name__ == "__main__":
    main()
