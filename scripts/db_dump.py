#!/usr/bin/env python3
"""
Create, verify, restore and fetch a dump of the KG-TransomicNet ArangoDB
instance, including the quantitative layers.

Why this exists rather than a bare arangodump call: a dump of this database can
fail *silently*. A run that exits 0, whose gzip streams are intact and whose
last line is valid JSON, can still be missing documents. Restoring it yields a
database that looks healthy and quietly under-reports a layer. Every command
here therefore counts documents and refuses to declare success on a mismatch.

    # create a dump and verify it against the live database
    python scripts/db_dump.py dump --db PKT_main --out /path/to/dump

    # check a dump you already have (no server needed)
    python scripts/db_dump.py verify --input /path/to/dump

    # restore into a database, then verify what landed
    python scripts/db_dump.py restore --input /path/to/dump --db PKT_main --create

    # fetch the published dump from Hugging Face
    python scripts/db_dump.py download --out /path/to/dump

Notes
-----
* Do not write a dump into a cloud-synced folder (Google Drive, Dropbox,
  OneDrive). The sync client rewrites files under the writer's feet and is a
  known cause of silent truncation. --out refuses such paths unless you pass
  --allow-synced-dir.
* The full instance is ~11.4 GB restored; the compressed dump is ~3.5 GB.
"""
from __future__ import annotations

import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

MANIFEST = "kg_transomicnet_manifest.json"

# Collection groups, so the backbone and the heavy layers can be handled apart.
BACKBONE = ["nodes", "edges"]
QUANTITATIVE = ["GENE_EXPRESSION", "CNV", "MIRNA", "PROTEIN", "METHYLATION"]
SUPPORT = ["GENES", "SAMPLES", "CASES", "PROJECTS"]
ALL_COLLECTIONS = BACKBONE + QUANTITATIVE + SUPPORT

HF_REPO = "johndef64/KG-TransomicNet"

SYNCED_MARKERS = ("google drive", "googledrive", "dropbox", "onedrive",
                  "icloud", "\\my drive", "/my drive")


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------

def log(msg: str) -> None:
    print(msg, flush=True)


def die(msg: str) -> "NoReturn":  # noqa: F821
    print(f"[ERROR] {msg}", file=sys.stderr, flush=True)
    sys.exit(1)


def find_tool(name: str, explicit: str | None) -> str:
    """Locate arangodump/arangorestore: --arango-bin, then PATH, then defaults."""
    if explicit:
        cand = Path(explicit) / name
        for p in (cand, cand.with_suffix(".exe")):
            if p.is_file():
                return str(p)
        die(f"{name} not found under {explicit}")

    found = shutil.which(name)
    if found:
        return found

    patterns = [
        "C:/Program Files/ArangoDB3*/usr/bin",
        "D:/Program Files/ArangoDB3*/usr/bin",
        "/usr/bin", "/usr/local/bin", "/opt/arangodb/bin",
    ]
    import glob
    for pat in patterns:
        for d in sorted(glob.glob(pat), reverse=True):
            for p in (Path(d) / name, Path(d) / f"{name}.exe"):
                if p.is_file():
                    return str(p)
    die(f"{name} not found. Put it on PATH or pass --arango-bin <dir>.")


def is_synced_dir(path: Path) -> bool:
    p = str(path.resolve()).lower()
    return any(m in p for m in SYNCED_MARKERS)


def connect(db_name: str, host: str, user: str, password: str):
    from arango import ArangoClient
    from arango.http import DefaultHTTPClient
    client = ArangoClient(hosts=host,
                          http_client=DefaultHTTPClient(request_timeout=3600))
    return client.db(db_name, username=user, password=password)


def live_counts(db, collections: list[str]) -> dict[str, int]:
    out = {}
    for c in collections:
        if db.has_collection(c):
            out[c] = db.collection(c).count()
    return out


def count_gz_lines(path: Path) -> int:
    """Documents in a dump data file: arangodump writes one JSON object per line."""
    n = 0
    with gzip.open(path, "rb") as fh:
        while True:
            chunk = fh.read(32 << 20)
            if not chunk:
                break
            n += chunk.count(b"\n")
    return n


def count_plain_lines(path: Path) -> int:
    n = 0
    with open(path, "rb") as fh:
        while True:
            chunk = fh.read(32 << 20)
            if not chunk:
                break
            n += chunk.count(b"\n")
    return n


def dump_files(dump_dir: Path) -> dict[str, list[Path]]:
    """Map collection name -> its data files (arangodump may write several)."""
    out: dict[str, list[Path]] = {}
    for p in sorted(dump_dir.iterdir()):
        name = p.name
        if ".data.json" not in name:
            continue
        coll = name.split(".data.json")[0]
        # strip the md5 suffix arangodump appends: NAME_<32 hex>
        if "_" in coll:
            head, _, tail = coll.rpartition("_")
            if len(tail) == 32 and all(ch in "0123456789abcdef" for ch in tail):
                coll = head
        out.setdefault(coll, []).append(p)
    return out


def count_dump(dump_dir: Path) -> dict[str, int]:
    counts = {}
    for coll, files in dump_files(dump_dir).items():
        total = 0
        for f in files:
            total += (count_gz_lines(f) if f.suffix == ".gz"
                      else count_plain_lines(f))
        counts[coll] = total
    return counts


def compare(expected: dict[str, int], actual: dict[str, int]) -> list[str]:
    """Return a list of human-readable problems; empty means all good."""
    problems = []
    for coll in sorted(set(expected) | set(actual)):
        exp, act = expected.get(coll), actual.get(coll)
        if exp is None:
            problems.append(f"{coll}: present in dump ({act:,}) but not expected")
        elif act is None:
            problems.append(f"{coll}: expected {exp:,} documents, absent from dump")
        elif exp != act:
            problems.append(f"{coll}: expected {exp:,}, found {act:,} "
                            f"({act - exp:+,})")
    return problems


def report_table(expected: dict[str, int], actual: dict[str, int]) -> None:
    log(f"  {'collection':<20} {'expected':>14} {'in dump':>14}   ")
    for coll in sorted(set(expected) | set(actual)):
        exp = expected.get(coll)
        act = actual.get(coll)
        mark = "ok" if exp == act else "MISMATCH"
        log(f"  {coll:<20} {(f'{exp:,}' if exp is not None else '-'):>14} "
            f"{(f'{act:,}' if act is not None else '-'):>14}   {mark}")


def run(cmd: list[str]) -> int:
    log("  $ " + " ".join(f'"{c}"' if " " in c else c for c in cmd[:6]) + " ...")
    t0 = time.perf_counter()
    p = subprocess.run(cmd)
    log(f"  finished in {time.perf_counter() - t0:,.0f}s (exit {p.returncode})")
    return p.returncode


# --------------------------------------------------------------------------
# commands
# --------------------------------------------------------------------------

def cmd_dump(args) -> int:
    out = Path(args.out)
    if is_synced_dir(out) and not args.allow_synced_dir:
        die(f"{out} looks like a cloud-synced folder. The sync client can "
            f"rewrite files while arangodump writes them, which silently "
            f"truncates the dump. Choose a local path, or pass "
            f"--allow-synced-dir if you are sure.")

    collections = args.collections or ALL_COLLECTIONS
    db = connect(args.db, args.host, args.user, args.password)
    expected = live_counts(db, collections)
    missing = [c for c in collections if c not in expected]
    if missing:
        die(f"not in database {args.db}: {', '.join(missing)}")

    log(f"dumping {len(expected)} collections from {args.db}")
    for c, n in sorted(expected.items()):
        log(f"  {c:<20} {n:>14,} documents")

    if out.exists() and args.overwrite:
        shutil.rmtree(out)
    out.mkdir(parents=True, exist_ok=True)

    cmd = [find_tool("arangodump", args.arango_bin),
           # ArangoDB derives its config filename from argv[0]; on Windows
           # shutil.which returns "...arangodump.EXE", which sends it looking
           # for "arangodump.EXE.conf" and aborting. We need no config file.
           "--configuration", "none",
           "--server.endpoint", args.host.replace("http://", "tcp://"),
           "--server.username", args.user,
           "--server.password", args.password,
           "--server.database", args.db,
           "--output-directory", str(out),
           "--compress-output", "true",
           "--overwrite", "true",
           "--threads", str(args.threads)]
    for c in collections:
        cmd += ["--collection", c]

    rc = run(cmd)
    if rc != 0:
        die(f"arangodump exited {rc}")

    log("verifying (counting documents in the dump — this reads every file)")
    actual = count_dump(out)
    report_table(expected, actual)
    problems = compare(expected, actual)

    manifest = {
        "source_database": args.db,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "arangodump": find_tool("arangodump", args.arango_bin),
        "collections": {c: {"documents": expected[c]} for c in sorted(expected)},
        "verified": not problems,
    }
    (out / MANIFEST).write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    if problems:
        for p in problems:
            log(f"  ! {p}")
        die("dump is INCOMPLETE — do not publish it. Retry with fewer threads "
            "(--threads 1) and a local, non-synced output directory.")

    total = sum(f.stat().st_size for f in out.iterdir() if f.is_file())
    log(f"dump verified complete: {sum(expected.values()):,} documents, "
        f"{total / 1e9:.2f} GB on disk")
    log(f"manifest written to {out / MANIFEST}")
    return 0


def cmd_verify(args) -> int:
    inp = find_dump_dir(Path(args.input))
    man_path = inp / MANIFEST
    if man_path.is_file():
        expected = {c: v["documents"] for c, v in
                    json.loads(man_path.read_text(encoding="utf-8"))["collections"].items()}
        log(f"expected counts from {MANIFEST}")
    elif args.db:
        db = connect(args.db, args.host, args.user, args.password)
        # Without a manifest there is no record of what the dump was meant to
        # hold, so --collections says which ones to require.
        wanted = args.collections or ALL_COLLECTIONS
        expected = live_counts(db, wanted)
        log(f"expected counts from live database {args.db} "
            f"({len(expected)} collections)")
    else:
        die(f"no {MANIFEST} in {inp}; pass --db to compare against a live database")

    actual = count_dump(inp)
    report_table(expected, actual)
    problems = compare(expected, actual)
    if problems:
        for p in problems:
            log(f"  ! {p}")
        die("dump is INCOMPLETE")
    log(f"dump verified complete: {sum(actual.values()):,} documents")
    return 0


def find_dump_dir(root: Path) -> Path:
    """Locate the dump inside `root`, which may be a parent (e.g. a download).

    A dump directory is the one holding the *.structure.json files, so a repo
    that keeps the dump in a subfolder restores just as well as a flat one.
    """
    if not root.is_dir():
        die(f"{root} is not a directory")
    if any(root.glob("*.structure.json")):
        return root
    candidates = sorted({p.parent for p in root.rglob("*.structure.json")})
    if not candidates:
        die(f"no dump found under {root} (no *.structure.json anywhere)")
    if len(candidates) > 1:
        names = ", ".join(str(c.relative_to(root)) for c in candidates)
        die(f"several dumps found under {root}: {names}. Point --input at one.")
    log(f"  dump found at {candidates[0]}")
    return candidates[0]


def cmd_restore(args) -> int:
    inp = find_dump_dir(Path(args.input))

    if args.create:
        from arango import ArangoClient
        sys_db = ArangoClient(hosts=args.host).db("_system", username=args.user,
                                                  password=args.password)
        if not sys_db.has_database(args.db):
            sys_db.create_database(args.db)
            log(f"created database {args.db}")

    cmd = [find_tool("arangorestore", args.arango_bin),
           "--configuration", "none",   # see the note in cmd_dump
           "--server.endpoint", args.host.replace("http://", "tcp://"),
           "--server.username", args.user,
           "--server.password", args.password,
           "--server.database", args.db,
           "--input-directory", str(inp),
           "--threads", str(args.threads)]
    rc = run(cmd)
    if rc != 0:
        die(f"arangorestore exited {rc}")

    man_path = inp / MANIFEST
    if not man_path.is_file():
        log(f"[warning] no {MANIFEST} in the dump; cannot verify what landed")
        return 0

    expected = {c: v["documents"] for c, v in
                json.loads(man_path.read_text(encoding="utf-8"))["collections"].items()}
    db = connect(args.db, args.host, args.user, args.password)
    actual = live_counts(db, list(expected))
    report_table(expected, actual)
    problems = compare(expected, actual)
    if problems:
        for p in problems:
            log(f"  ! {p}")
        die(f"restore into {args.db} is INCOMPLETE")
    log(f"restore verified: {sum(actual.values()):,} documents in {args.db}")
    return 0


def cmd_download(args) -> int:
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    try:
        from huggingface_hub import snapshot_download
    except ImportError:
        die("huggingface_hub is not installed: pip install huggingface_hub")

    log(f"downloading {args.repo} -> {out}")
    snapshot_download(repo_id=args.repo, repo_type="dataset",
                      local_dir=str(out), allow_patterns=args.include or None)

    log("download complete; verifying")
    dump_dir = find_dump_dir(out)
    args.input = str(dump_dir)
    args.db = None
    args.collections = None
    rc = cmd_verify(args)
    if rc == 0:
        log("")
        log("restore it with:")
        log(f"  python scripts/db_dump.py restore --input {dump_dir} "
            f"--db PKT_main --create")
    return rc


# --------------------------------------------------------------------------

def main() -> int:
    try:
        import arangodb_utils as au
        d_host, d_user, d_pass = (au.arangodb_hosts, au.arangodb_user,
                                  au.arangodb_password)
        d_db = au.db_name
    except Exception:
        d_host, d_user, d_pass, d_db = ("http://localhost:8529", "root", "", "PKT_main")

    p = argparse.ArgumentParser(
        description="Create, verify, restore or download a KG-TransomicNet database dump.",
        formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)

    def common(sp, need_db=True):
        sp.add_argument("--host", default=d_host)
        sp.add_argument("--user", default=d_user)
        sp.add_argument("--password", default=os.environ.get("ARANGO_PASSWORD", d_pass))
        sp.add_argument("--arango-bin", default=None,
                        help="directory holding arangodump/arangorestore")
        sp.add_argument("--threads", type=int, default=2,
                        help="parallel collections (default: 2; lower is safer "
                             "for very large documents)")
        if need_db:
            sp.add_argument("--db", default=d_db,
                            help=f"ArangoDB database name (default: {d_db}).")

    sp = sub.add_parser("dump", help="create a dump and verify it")
    common(sp)
    sp.add_argument("--out", required=True)
    sp.add_argument("--collections", nargs="+",
                    help=f"default: all ({len(ALL_COLLECTIONS)} collections)")
    sp.add_argument("--overwrite", action="store_true",
                    help="delete the output directory first")
    sp.add_argument("--allow-synced-dir", action="store_true")
    sp.set_defaults(func=cmd_dump)

    sp = sub.add_parser("verify", help="count documents in an existing dump")
    common(sp, need_db=False)
    sp.add_argument("--input", required=True)
    sp.add_argument("--db", default=None,
                    help="compare against this live database instead of the manifest")
    sp.add_argument("--collections", nargs="+",
                    help="with --db: the collections the dump should hold "
                         "(default: all)")
    sp.set_defaults(func=cmd_verify)

    sp = sub.add_parser("restore", help="restore a dump and verify what landed")
    common(sp)
    sp.add_argument("--input", required=True)
    sp.add_argument("--create", action="store_true",
                    help="create the target database if it does not exist")
    sp.set_defaults(func=cmd_restore)

    sp = sub.add_parser("download", help="fetch the published dump from Hugging Face")
    common(sp, need_db=False)
    sp.add_argument("--out", required=True)
    sp.add_argument("--repo", default=HF_REPO)
    sp.add_argument("--include", nargs="+",
                    help="restrict to these path patterns")
    sp.set_defaults(func=cmd_download)

    args = p.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
