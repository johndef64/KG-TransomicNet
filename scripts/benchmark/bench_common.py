"""
Shared helpers for the KG-TransomicNet database-engineering benchmarks.

All benchmark scripts write CSV to results/benchmark/ and never modify the
production database except where explicitly stated (see bench_query.py,
--create-indexes).

Usage:
    from bench_common import connect, time_query, write_csv, ...
"""
from __future__ import annotations

import csv
import json
import platform
import statistics
import time
from pathlib import Path

import requests
from arango import ArangoClient
from arango.http import DefaultHTTPClient

# --- connection ------------------------------------------------------------
HOST = "http://localhost:8529"
USER = "root"
PASSWORD = "avocadodb"
DB_NAME = "PKT_TransomicNet_v3"
SCRATCH_DB = "BENCH_scratch"

HTTP_TIMEOUT = 3600  # bulk pulls and index builds exceed the 60 s default

RESULTS_DIR = Path(__file__).resolve().parents[2] / "results" / "benchmark"


def client() -> ArangoClient:
    return ArangoClient(hosts=HOST,
                        http_client=DefaultHTTPClient(request_timeout=HTTP_TIMEOUT))


def connect(db_name: str = DB_NAME):
    return client().db(db_name, username=USER, password=PASSWORD)


def sys_db():
    return client().db("_system", username=USER, password=PASSWORD)


def make_scratch(name: str = SCRATCH_DB, drop_first: bool = True):
    """Create (optionally recreating) a scratch database for write experiments."""
    s = sys_db()
    if drop_first and s.has_database(name):
        s.delete_database(name)
    if not s.has_database(name):
        s.create_database(name)
    return connect(name)


def drop_scratch(name: str = SCRATCH_DB):
    s = sys_db()
    if s.has_database(name):
        s.delete_database(name)


# --- storage ---------------------------------------------------------------
def collection_figures(db_name: str, collection: str) -> dict:
    """Return count, documentsSize, indexesSize (bytes) via the REST API.

    python-arango 8.x dropped Collection.figures(), so we call the endpoint
    directly. Sizes are RocksDB estimates, which is what ArangoDB exposes.
    """
    url = f"{HOST}/_db/{db_name}/_api/collection/{collection}/figures"
    r = requests.get(url, auth=(USER, PASSWORD), timeout=120)
    r.raise_for_status()
    j = r.json()
    return {
        "count": j["count"],
        "documents_bytes": j["figures"]["documentsSize"],
        "indexes_bytes": j["figures"]["indexes"]["size"],
    }


def list_collections(db_name: str) -> list[str]:
    url = f"{HOST}/_db/{db_name}/_api/collection?excludeSystem=true"
    r = requests.get(url, auth=(USER, PASSWORD), timeout=120)
    r.raise_for_status()
    return sorted(c["name"] for c in r.json()["result"])


# --- timing ----------------------------------------------------------------
def time_query(db, aql: str, bind_vars: dict | None = None,
               repeats: int = 10, warmup: int = 2) -> dict:
    """Execute an AQL query repeatedly and return latency statistics in ms.

    `warmup` executions are discarded (they populate the RocksDB block cache),
    so the reported figures are warm-cache latencies. Cold-cache figures are
    obtained separately in bench_query.py by restarting the cache.
    """
    for _ in range(warmup):
        list(db.aql.execute(aql, bind_vars=bind_vars))

    lat = []
    n_rows = 0
    for _ in range(repeats):
        t0 = time.perf_counter()
        rows = list(db.aql.execute(aql, bind_vars=bind_vars))
        lat.append((time.perf_counter() - t0) * 1000.0)
        n_rows = len(rows)

    lat.sort()
    return {
        "n_rows": n_rows,
        "repeats": repeats,
        "mean_ms": round(statistics.fmean(lat), 2),
        "p50_ms": round(statistics.median(lat), 2),
        "p95_ms": round(lat[min(len(lat) - 1, int(0.95 * len(lat)))], 2),
        "min_ms": round(lat[0], 2),
        "max_ms": round(lat[-1], 2),
    }


# --- output ----------------------------------------------------------------
def write_csv(rows: list[dict], filename: str) -> Path:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    path = RESULTS_DIR / filename
    if not rows:
        raise ValueError(f"no rows to write for {filename}")
    fields = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"  -> wrote {path} ({len(rows)} rows)")
    return path


def write_env(filename: str = "environment.json") -> Path:
    """Record the hardware/software environment: required to interpret timings."""
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    try:
        import psutil
        mem_gb = round(psutil.virtual_memory().total / 1e9, 1)
        cpus = psutil.cpu_count(logical=True)
        cpus_phys = psutil.cpu_count(logical=False)
    except Exception:
        mem_gb, cpus, cpus_phys = None, None, None

    ver = requests.get(f"{HOST}/_api/version", auth=(USER, PASSWORD), timeout=30).json()
    env = {
        "arangodb_version": ver.get("version"),
        "arangodb_license": ver.get("license"),
        "storage_engine": "rocksdb",
        "python": platform.python_version(),
        "platform": platform.platform(),
        "processor": platform.processor(),
        "cpu_logical": cpus,
        "cpu_physical": cpus_phys,
        "ram_gb": mem_gb,
        "timestamp": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    path = RESULTS_DIR / filename
    path.write_text(json.dumps(env, indent=2), encoding="utf-8")
    print(f"  -> wrote {path}")
    return path
