"""
run_query.py
============
Run any .aql file of this collection with bind variables and save the result.

Bind variables come from a JSON file (--bind), from key=value pairs (--set,
values parsed as JSON when possible), or both (--set wins).

Examples
--------
    # Two-step retrieval on TCGA-BRCA, printed to screen
    python queries/run_query.py queries/examples/01_two_step_retrieval.aql \
        --set cohort=TCGA-BRCA --set seed=nodes/HP_0003002 --limit 5

    # Typed subgraph with a routing table, saved as JSON lines
    python queries/run_query.py queries/extraction/typed_relation_subgraph.aql \
        --bind queries/params/typed_relations.json --out triples.jsonl
"""

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import add_connection_args, connect_from_args, read_json, run  # noqa: E402


def parse_set(pairs):
    out = {}
    for pair in pairs or []:
        key, _, raw = pair.partition("=")
        try:
            out[key] = json.loads(raw)
        except json.JSONDecodeError:
            out[key] = raw
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("aql", help="path to the .aql file")
    p.add_argument("--bind", help="JSON file with bind variables")
    p.add_argument("--set", action="append", metavar="KEY=VALUE", help="single bind variable")
    p.add_argument("--out", help="output file (.jsonl); default: print to screen")
    p.add_argument("--limit", type=int, help="stop after N results")
    add_connection_args(p)
    args = p.parse_args()

    bind = read_json(args.bind) if args.bind else {}
    bind.update(parse_set(args.set))
    db = connect_from_args(args)

    sink = open(args.out, "w", encoding="utf-8") if args.out else sys.stdout
    n = 0
    try:
        for row in run(db, args.aql, bind):
            sink.write(json.dumps(row, ensure_ascii=False) + "\n")
            n += 1
            if args.limit and n >= args.limit:
                break
    finally:
        if args.out:
            sink.close()
    print(f"[run_query] {n:,} results", file=sys.stderr)


if __name__ == "__main__":
    main()
