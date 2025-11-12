"""
Essential script to create small sample CSV files from MetaKG entities and triples.

- Defaults:
  - Entities:  `metakg_entities.csv`
  - Triples:   `metakg_triples-001.csv`
  - Output:    same folder, files named `sample_entities_<n>.csv` and `sample_triples_<n>.csv`
  - Sample N:  1000 rows (change with `-n`)
  - Seed:      42 (change with `--seed`)

Implementation uses standard library only. By default it uses fast head sampling
(first N rows) to avoid scanning huge files. Optionally, use uniform random
sampling with reservoir sampling via `--mode random`.
"""

from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path
import random
from typing import Iterable, List, Optional, Tuple


def detect_dialect(file_path: Path, default=csv.excel, enable_sniff: bool = False) -> csv.Dialect:
	"""Optionally sniff CSV dialect; by default return excel (comma-delimited).

	Sniffing can be slow on very large files, so it's disabled unless requested.
	"""
	if not enable_sniff:
		return default
	try:
		with open(file_path, "r", encoding="utf-8", errors="replace") as f:
			sample = f.read(64_000)
		sniffer = csv.Sniffer()
		return sniffer.sniff(sample)
	except Exception:
		return default


def reservoir_sample_rows(
	file_path: Path,
	n: int,
	seed: int = 42,
	enable_sniff: bool = False,
) -> Tuple[List[str], List[List[str]]]:
	"""Reservoir sample n data rows (excluding header) from a CSV file.

	Returns: (header, rows)
	"""
	if n <= 0:
		raise ValueError("n must be > 0")

	dialect = detect_dialect(file_path, enable_sniff=enable_sniff)
	rng = random.Random(seed)

	header: Optional[List[str]] = None
	reservoir: List[List[str]] = []
	count = 0

	with open(file_path, "r", encoding="utf-8", errors="replace", newline="") as f:
		reader = csv.reader(f, dialect)
		try:
			header = next(reader)
		except StopIteration:
			# empty file
			return [], []
		for row in reader:
			if count < n:
				reservoir.append(row)
			else:
				j = rng.randint(0, count)
				if j < n:
					reservoir[j] = row
			count += 1

	return header or [], reservoir


def head_sample_rows(
	file_path: Path,
	n: int,
	enable_sniff: bool = False,
) -> Tuple[List[str], List[List[str]]]:
	"""Read header and the first n rows quickly (no full scan)."""
	if n <= 0:
		raise ValueError("n must be > 0")
	dialect = detect_dialect(file_path, enable_sniff=enable_sniff)
	rows: List[List[str]] = []
	with open(file_path, "r", encoding="utf-8", errors="replace", newline="") as f:
		reader = csv.reader(f, dialect)
		try:
			header = next(reader)
		except StopIteration:
			return [], []
		for i, row in enumerate(reader):
			if i >= n:
				break
			rows.append(row)
	return header, rows


def write_csv(file_path: Path, header: List[str], rows: List[List[str]], dialect: csv.Dialect | None = None) -> None:
	file_path.parent.mkdir(parents=True, exist_ok=True)
	if dialect is None:
		dialect = csv.excel
	with open(file_path, "w", encoding="utf-8", newline="") as f:
		writer = csv.writer(f, dialect)
		if header:
			writer.writerow(header)
		writer.writerows(rows)


def make_sample(
	src: Path,
	out: Path,
	n: int,
	seed: int,
	mode: str = "head",
	sniff: bool = False,
) -> None:
	if not src.exists():
		print(f"[WARN] Input not found, skipping: {src}")
		return
	print(f"[INFO] Sampling {n} rows from {src.name} (mode={mode}, seed={seed})…")
	if mode == "random":
		header, rows = reservoir_sample_rows(src, n=n, seed=seed, enable_sniff=sniff)
	else:
		header, rows = head_sample_rows(src, n=n, enable_sniff=sniff)
	write_csv(out, header, rows)
	print(f"[OK] Wrote {len(rows)} rows to {out}")


def parse_args() -> argparse.Namespace:
	here = Path(__file__).resolve().parent
	parser = argparse.ArgumentParser(
		description="Create small sample CSVs from MetaKG entities and triples (essential)."
	)
	parser.add_argument(
		"-n", "--num-rows", type=int, default=1000, help="Number of rows to sample from each input"
	)
	parser.add_argument("--seed", type=int, default=42, help="Random seed for sampling")
	parser.add_argument(
		"--entities",
		type=Path,
		default=here / "metakg_entities.csv",
		help="Path to entities CSV (default: metakg_entities.csv in this folder)",
	)
	parser.add_argument(
		"--triples",
		type=Path,
		default=here / "metakg_triples-001.csv",
		help="Path to triples CSV (default: metakg_triples-001.csv in this folder)",
	)
	parser.add_argument(
		"--outdir",
		type=Path,
		default=here,
		help="Output directory (default: same as this script)",
	)
	parser.add_argument(
		"--mode",
		choices=["head", "random"],
		default="head",
		help="Sampling mode: 'head' (fast, first N rows) or 'random' (uniform reservoir sampling, scans full file)",
	)
	parser.add_argument(
		"--sniff",
		action="store_true",
		help="Enable CSV dialect sniffing (can be slow on large files).",
	)
	parser.add_argument(
		"--entities-only",
		action="store_true",
		help="Only create entities sample",
	)
	parser.add_argument(
		"--triples-only",
		action="store_true",
		help="Only create triples sample",
	)
	return parser.parse_args()


def main() -> None:
	args = parse_args()
	n = max(1, int(args.num_rows))
	seed = int(args.seed)

	outdir: Path = args.outdir
	outdir.mkdir(parents=True, exist_ok=True)

	do_entities = True
	do_triples = True
	if args.entities_only and not args.triples_only:
		do_entities, do_triples = True, False
	elif args.triples_only and not args.entities_only:
		do_entities, do_triples = False, True

	if do_entities:
		entities_in: Path = args.entities
		entities_out = outdir / f"sample_entities_{n}.csv"
		make_sample(entities_in, entities_out, n=n, seed=seed, mode=args.mode, sniff=args.sniff)

	if do_triples:
		triples_in: Path = args.triples
		triples_out = outdir / f"sample_triples_{n}.csv"
		make_sample(triples_in, triples_out, n=n, seed=seed, mode=args.mode, sniff=args.sniff)

	print("[DONE]")


if __name__ == "__main__":
	main()
