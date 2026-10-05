"""
to_typed_triples.py
===================
Export a type-constrained relational subgraph of the backbone as a typed-triple
TSV, the input format of knowledge-graph embedding and relational GNN frameworks:

    head <TAB> interaction <TAB> tail <TAB> source <TAB> type
    entity = "<Type>::<entity_id>"      (node type = prefix before "::")

The subgraph comes from extraction/typed_relation_subgraph.aql with a routing
table (default: params/typed_relations.json). Relations from other resources
(e.g. a curated drug--target layer) can be appended with --extra: they are
external tables keyed on the same identifiers, and their rows keep their own
`source` label.

Rows are written sorted (relation, then head, then tail), so two exports of the
same database are byte-identical: a positional train/test split stays stable.

Usage
-----
    python queries/adapters/to_typed_triples.py --out triples.tsv
    python queries/adapters/to_typed_triples.py --relations PPI,GENE_PRODUCT --zip
    python queries/adapters/to_typed_triples.py --extra extra_layer.tsv --out graph.tsv
"""

import argparse
import csv
import io
import sys
import zipfile
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import QUERIES_DIR, add_connection_args, connect_from_args, read_json, run  # noqa: E402

QUERY = QUERIES_DIR / "extraction" / "typed_relation_subgraph.aql"
DEFAULT_ROUTING = QUERIES_DIR / "params" / "typed_relations.json"
SOURCE = "PKT"

# bioentity_type -> entity prefix
TYPE_PREFIX = {
    "chemical": "Compound", "protein": "Protein", "gene": "Gene", "go": "GO",
    "pathway": "Pathway", "disease": "Disease", "phenotype": "Phenotype",
    "variant": "Variant",
}


def extract(db, routing_file):
    params = read_json(routing_file)
    edges = defaultdict(set)
    for r in run(db, QUERY, {"routing": params["routing"], "predicates": params["predicates"]}):
        h = f"{TYPE_PREFIX.get(r['head_type'], r['head_type'])}::{r['head_id']}"
        t = f"{TYPE_PREFIX.get(r['tail_type'], r['tail_type'])}::{r['tail_id']}"
        edges[r["relation"]].add((h, t))
    return edges


def load_extra(paths, edges, rel_source):
    """Append external relations: TSV with head, interaction, tail[, source]."""
    for path in paths or []:
        with open(path, encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rel = row["interaction"]
                edges[rel].add((row["head"], row["tail"]))
                rel_source[rel] = row.get("source") or Path(path).stem


def write(out, edges, relations, rel_source, as_zip):
    buf = io.StringIO()
    w = csv.writer(buf, delimiter="\t", lineterminator="\n")
    w.writerow(["head", "interaction", "tail", "source", "type"])
    n = 0
    for rel in sorted(relations):
        for h, t in sorted(edges.get(rel, ())):
            w.writerow([h, rel, t, rel_source.get(rel, SOURCE),
                        f"{h.split('::', 1)[0]}-{t.split('::', 1)[0]}"])
            n += 1
    out = Path(out)
    if as_zip:
        with zipfile.ZipFile(out.with_suffix(out.suffix + ".zip"), "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr(out.name, buf.getvalue())
    else:
        out.write_bytes(buf.getvalue().encode("utf-8"))   # "\n" line ends on every OS
    return n


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--routing", default=str(DEFAULT_ROUTING), help="routing table (JSON)")
    p.add_argument("--relations", help="comma-separated subset to write (default: all)")
    p.add_argument("--extra", action="append", help="external relation TSV to append")
    p.add_argument("--out", default="typed_triples.tsv")
    p.add_argument("--zip", action="store_true", help="write <out>.zip instead of <out>")
    add_connection_args(p)
    args = p.parse_args()

    db = connect_from_args(args)
    edges = extract(db, args.routing)
    rel_source = {}
    load_extra(args.extra, edges, rel_source)
    relations = args.relations.split(",") if args.relations else list(edges)
    n = write(args.out, edges, relations, rel_source, args.zip)

    for rel in sorted(relations):
        print(f"  {rel:<24} {len(edges.get(rel, ())):>9,}  ({rel_source.get(rel, SOURCE)})")
    print(f"[to_typed_triples] {n:,} triples -> {args.out}{'.zip' if args.zip else ''}")


if __name__ == "__main__":
    main()
