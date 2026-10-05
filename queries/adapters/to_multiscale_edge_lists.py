"""
to_multiscale_edge_lists.py
===========================
Export the higher biological scales of a gene panel as node and edge lists, the
input of heterogeneous graph builders (one CSV per node type, one per relation):

    node_gene.csv        symbol, idx
    node_pathway.csv     reactome_id, idx
    node_GO_term.csv     go_id, idx
    node_disease.csv     mondo_id, idx
    edge_gene_member_of_pathway.csv       gene_idx, gene, pathway_idx, reactome_id
    edge_gene_annotated_with_GO_term.csv  gene_idx, gene, GO_idx, go_id
    edge_GO_term_is_a_GO_term.csv         GO_src_idx, go_src, GO_dst_idx, go_dst
    edge_gene_associated_with_disease.csv gene_idx, gene, disease_idx, mondo_id

Queries: extraction/gene_scale_edges.aql (gene -> pathway, gene -> disease),
extraction/gene_ontology_annotation.aql (gene -> GO through the protein bridge,
one gene per protein, minimum support) and extraction/ontology_closure.aql
(is_a / part_of among the retained terms).

The panel is a CSV with a gene column (HGNC symbols by default). Symbols are
resolved to Entrez with data/mappings/hgcn_mart_conversion_table.zip (approved
symbols); node_gene.csv keeps the panel order. Indices of the other node types
follow the sorted identifiers. Interaction layers between genes (protein
interactions, regulatory edges) are not part of the backbone and are not
produced here.

Usage
-----
    python queries/adapters/to_multiscale_edge_lists.py --genes panel.csv --out-dir hetero/
    python queries/adapters/to_multiscale_edge_lists.py --genes panel.csv --go-min-support 3 \
        --no-disease --out-dir hetero_no_disease/
"""

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from common import PROJECT_ROOT, QUERIES_DIR, add_connection_args, connect_from_args, run  # noqa: E402

Q = QUERIES_DIR / "extraction"
HGNC_TABLE = PROJECT_ROOT / "data" / "mappings" / "hgcn_mart_conversion_table.zip"

SCALES = {
    "pathway": {"name": "pathway", "target_cc": "R-HSA", "predicates": ["participates in"]},
    "disease": {"name": "disease", "target_cc": "MONDO",
                "predicates": ["causes or contributes to condition"]},
}
GO_PREDICATES = ["participates in", "has function", "located_in"]
GO_HIERARCHY = ["type", "part_of"]


def symbol_to_entrez(path=HGNC_TABLE):
    conv = pd.read_csv(path, sep="\t", compression="zip")
    conv["NCBI gene ID"] = pd.to_numeric(conv["NCBI gene ID"], errors="coerce").astype("Int64")
    conv = conv.dropna(subset=["NCBI gene ID", "Approved symbol"]).drop_duplicates(subset=["Approved symbol"])
    return {s: str(int(e)) for s, e in zip(conv["Approved symbol"], conv["NCBI gene ID"])}


def load_panel(path, column, id_type):
    df = pd.read_csv(path)
    col = column or next(c for c in ("symbol", "gene", "entrez_id") if c in df.columns)
    ids = df[col].astype(str).tolist()
    if id_type == "entrez":
        return ids, {e: e for e in ids}
    s2e = symbol_to_entrez()
    kept = [s for s in ids if s in s2e]
    if len(kept) < len(ids):
        print(f"[panel] {len(ids) - len(kept)} symbols without Entrez dropped (kept {len(kept)})")
    return kept, {s: s2e[s] for s in kept}


def vocab(ids):
    v = sorted(set(ids))
    return v, {x: i for i, x in enumerate(v)}


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="CSV with the gene panel")
    p.add_argument("--column", help="gene column (default: symbol | gene | entrez_id)")
    p.add_argument("--id-type", choices=["symbol", "entrez"], default="symbol")
    p.add_argument("--go-min-support", type=int, default=3)
    p.add_argument("--all-genes-per-protein", action="store_true",
                   help="assign a protein to every linked gene (default: one gene per protein)")
    p.add_argument("--no-pathway", action="store_true")
    p.add_argument("--no-go", action="store_true")
    p.add_argument("--no-disease", action="store_true")
    p.add_argument("--out-dir", default="multiscale_edges")
    add_connection_args(p)
    args = p.parse_args()

    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    genes, sym2e = load_panel(args.genes, args.column, args.id_type)
    entrez2idx = {}
    for i, g in enumerate(genes):
        entrez2idx[sym2e[g]] = i              # same as a dict comprehension: last symbol wins
    idx2sym = dict(enumerate(genes))
    panel = list(entrez2idx)
    db = connect_from_args(args)

    pd.DataFrame({"symbol": genes, "idx": range(len(genes))}).to_csv(out / "node_gene.csv", index=False, lineterminator="\n")
    summary = [("gene", "node", len(genes))]

    scales = [SCALES[s] for s in ("pathway", "disease") if not getattr(args, f"no_{s}")]
    rows = list(run(db, Q / "gene_scale_edges.aql", {"panel": panel, "scales": scales})) if scales else []
    for name, node_file, id_col, edge_file in (
            ("pathway", "node_pathway.csv", "reactome_id", "edge_gene_member_of_pathway.csv"),
            ("disease", "node_disease.csv", "mondo_id", "edge_gene_associated_with_disease.csv")):
        if getattr(args, f"no_{name}"):
            continue
        pairs = [(entrez2idx[r["gene"]], r["target"]) for r in rows if r["scale"] == name]
        v, x2i = vocab(t for _, t in pairs)
        pd.DataFrame({id_col: v, "idx": range(len(v))}).to_csv(out / node_file, index=False, lineterminator="\n")
        pd.DataFrame([(gi, idx2sym[gi], x2i[t], t) for gi, t in pairs],
                     columns=["gene_idx", "gene", f"{name}_idx", id_col]).to_csv(out / edge_file, index=False, lineterminator="\n")
        summary += [(name, "node", len(v)), (name, "edge", len(pairs))]

    if not args.no_go:
        ann = list(run(db, Q / "gene_ontology_annotation.aql", {
            "panel": panel, "annotation_predicates": GO_PREDICATES,
            "min_support": args.go_min_support,
            "one_gene_per_protein": not args.all_genes_per_protein}))
        pairs = [(entrez2idx[r["gene"]], r["term"]) for r in ann]
        v, x2i = vocab(t for _, t in pairs)
        pd.DataFrame({"go_id": v, "idx": range(len(v))}).to_csv(out / "node_GO_term.csv", index=False, lineterminator="\n")
        pd.DataFrame([(gi, idx2sym[gi], x2i[t], t) for gi, t in pairs],
                     columns=["gene_idx", "gene", "GO_idx", "go_id"]).to_csv(
            out / "edge_gene_annotated_with_GO_term.csv", index=False, lineterminator="\n")
        hier = list(run(db, Q / "ontology_closure.aql", {
            "terms": v, "class_code": "GO", "hierarchy_predicates": GO_HIERARCHY}))
        hier = sorted((x2i[h["src"]], h["src"], x2i[h["dst"]], h["dst"]) for h in hier)
        pd.DataFrame(hier, columns=["GO_src_idx", "go_src", "GO_dst_idx", "go_dst"]).to_csv(
            out / "edge_GO_term_is_a_GO_term.csv", index=False, lineterminator="\n")
        summary += [("GO_term", "node", len(v)), ("annotated_with", "edge", len(pairs)),
                    ("is_a", "edge", len(hier))]

    pd.DataFrame(summary, columns=["name", "kind", "count"]).to_csv(out / "summary.csv", index=False, lineterminator="\n")
    for name, kind, n in summary:
        print(f"  {kind:<5} {name:<16} {n:>7,}")
    print(f"[to_multiscale_edge_lists] -> {out}")


if __name__ == "__main__":
    main()
