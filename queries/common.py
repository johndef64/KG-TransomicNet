"""
common.py
=========
Shared helpers for the query collection: database connection and execution of
an .aql file with bind variables.

Connection defaults (host, user, password, database name) are read from
scripts/arangodb_utils.py, so the query collection and the rest of the
repository never drift apart. Every helper accepts explicit overrides.
"""

import json
import sys
from pathlib import Path
from typing import Dict, Iterator, Optional

from arango import ArangoClient

QUERIES_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = QUERIES_DIR.parent
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

if str(SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPTS_DIR))
import arangodb_utils as _cfg  # noqa: E402  (defaults only; nothing is executed on import)

DEFAULT_DB = _cfg.db_name
DEFAULT_HOST = _cfg.arangodb_hosts
DEFAULT_USER = _cfg.arangodb_user
DEFAULT_PASSWORD = _cfg.arangodb_password

# Extraction queries scan the whole edge collection (about 11 M documents):
# allow them to run for a while before the client gives up.
DEFAULT_TIMEOUT_S = 3600


def connect(db_name: str = DEFAULT_DB, host: str = DEFAULT_HOST, user: str = DEFAULT_USER,
            password: str = DEFAULT_PASSWORD, timeout: int = DEFAULT_TIMEOUT_S):
    """Connect to an EXISTING database (never creates one, as in arangodb_utils)."""
    client = ArangoClient(hosts=host, request_timeout=timeout)
    sys_db = client.db("_system", username=user, password=password)
    if not sys_db.has_database(db_name):
        raise SystemExit(f"[ERROR] Database '{db_name}' does not exist on {host}; "
                         f"pass --db with the name of an existing database.")
    return client.db(db_name, username=user, password=password)


def load_aql(path) -> str:
    """Read an .aql file; '//' comment lines are kept (AQL accepts them)."""
    p = Path(path)
    if not p.is_absolute() and not p.exists():
        p = QUERIES_DIR / p
    return p.read_text(encoding="utf-8")


def run(db, aql_path, bind_vars: Optional[Dict] = None, batch_size: int = 1_000,
        timeout: int = DEFAULT_TIMEOUT_S) -> Iterator:
    """Execute an .aql file and stream its results.

    Keys starting with "_" (e.g. "_comment" in a parameter file) are documentation
    and are not passed to the server, which rejects undeclared bind variables.
    """
    bind = {k: v for k, v in (bind_vars or {}).items() if not k.startswith("_")}
    # Plain (non-streaming) cursor, as in scripts/query_utils.py; max_runtime also
    # stops the query on the server if the client goes away.
    cursor = db.aql.execute(load_aql(aql_path), bind_vars=bind,
                            batch_size=batch_size, ttl=timeout, max_runtime=timeout)
    yield from cursor


def add_connection_args(parser):
    """Standard --db/--host/--user/--password options for every CLI here."""
    parser.add_argument("--db", default=DEFAULT_DB, help=f"database name (default: {DEFAULT_DB})")
    parser.add_argument("--host", default=DEFAULT_HOST)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--password", default=DEFAULT_PASSWORD)
    return parser


def connect_from_args(args):
    return connect(args.db, args.host, args.user, args.password)


def read_json(path):
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)
