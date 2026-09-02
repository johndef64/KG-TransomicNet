#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
spotcheck_mapping.py — D4 Mapping-accuracy spot-check (TCGA-BRCA cohort).

Audits a random sample of feature -> backbone-node mappings, drawn from the
per-layer `*_index_TCGA-BRCA` documents in ArangoDB, against AUTHORITATIVE
EXTERNAL identifier sources (BioMart, TCPA RPPA500, miRBase/HGNC, Illumina HM27
manifest).  This is deliberately an *external* comparison: re-running the build
pipeline and comparing it against itself would be tautological (it would measure
determinism, not accuracy).  See docs/TODO-D4-BRCA.md.

Definition of an ERROR
----------------------
A feature whose index mapping resolves it to a node / cross-reference that does
NOT match the authoritative identifier of that feature in the external source:
  * RNA-Seq / CNV : the index `entrez_id` recorded for an ENSG disagrees with the
                    Entrez ID BioMart assigns to that same ENSG.
  * RPPA          : the index (entrez_id, gene_symbol) for an antibody disagrees
                    with the TCPA RPPA500 antibody->gene metadata.
  * miRNA         : the index `hgnc_symbol` / `mirbase_id` for a miRNA disagrees
                    with the miRBase->HGNC reference for that miRNA accession.
  * Methylation   : the gene(s) the index associates with a CpG probe disagree
                    with the gene(s) the HM27 manifest assigns to that probe.

A mapping is also flagged when the backbone node it should resolve to does not
exist (for gene/protein layers the EntrezID node has `_key == entrez_id`).

Output
------
  * validation/spotcheck_results.csv  — per-feature audit log (every sampled row)
  * validation/spotcheck_summary.csv  — per-layer + Overall (N, errors, rate)
  * a LaTeX-ready table printed to stdout for Supplementary Table S3.

Reproducibility: fixed SEED (default 42); the sampled positions are deterministic
given the index contents.

Usage
-----
  conda run -n gnn python validation/spotcheck_mapping.py
  conda run -n gnn python validation/spotcheck_mapping.py --n 100 --seed 42 \
      --db PKT_main --cohort TCGA-BRCA
"""

import argparse
import ast
import csv
import os
import random
import sys

import pandas as pd

# --- repo paths -------------------------------------------------------------
HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(HERE)
SCRIPTS = os.path.join(REPO, "scripts")
MAPPINGS = os.path.join(REPO, "data", "mappings")
sys.path.insert(0, SCRIPTS)

from arango import ArangoClient  # noqa: E402


# ===========================================================================
# Identifier normalisation helpers
# ===========================================================================
def norm_entrez(value):
    """Normalise an Entrez ID to a bare integer string.

    The index stores it as a string ('7105'); the BioMart table stores it as a
    float-looking string ('7105.0') and the TCPA table as a plain int. Return
    None for empty / NaN.
    """
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() in {"nan", "none", "<na>"}:
        return None
    try:
        return str(int(float(s)))
    except (ValueError, TypeError):
        return s


def norm_symbol(value):
    """Normalise an HGNC symbol: upper-case, strip. The KG node labels drop the
    '-' from miRNA symbols (MIRLET7A-1 -> MIRLET7A1); we compare on the raw
    symbol and, for miRNA, also on the dash-stripped form."""
    if value is None:
        return None
    s = str(value).strip().upper()
    return s or None


# ===========================================================================
# ArangoDB
# ===========================================================================
def connect(db_name, host, user, password):
    client = ArangoClient(hosts=host)
    return client.db(db_name, username=user, password=password)


def load_index_mappings(db, collection, index_key, mapping_field):
    """Return the list of mapping dicts from a *_index document."""
    aql = f"""
    FOR doc IN {collection}
        FILTER doc._key == @key
        RETURN doc.{mapping_field}
    """
    res = list(db.aql.execute(aql, bind_vars={"key": index_key}))
    if not res or not res[0]:
        raise RuntimeError(
            f"Index '{index_key}' not found / empty in collection '{collection}'."
        )
    return res[0]


def entrez_nodes_exist(db, entrez_ids):
    """Return the subset of entrez_ids that exist as EntrezID backbone nodes
    (gene nodes have _key == entrez_id)."""
    keys = [e for e in entrez_ids if e]
    if not keys:
        return set()
    aql = """
    FOR n IN nodes
        FILTER n._key IN @keys
        FILTER n.class_code == "EntrezID"
        RETURN n._key
    """
    return set(db.aql.execute(aql, bind_vars={"keys": list(set(keys))}))


# ===========================================================================
# External ground-truth tables
# ===========================================================================
def load_biomart():
    """ENSG (base) -> set of authoritative Entrez IDs (BioMart extended).

    A given ENSG can legitimately carry several Entrez IDs (the table is exploded
    over Entrez/transcript). The index records ONE entrez_id per ENSG; we count a
    match if the index value is among the BioMart Entrez IDs for that ENSG.
    """
    path = os.path.join(MAPPINGS, "biomart_gene_mappings_filled_extended.tsv")
    df = pd.read_csv(path, sep="\t", usecols=["gene_stable_id", "entrez_id",
                                              "hgnc_symbol"], dtype=str)
    ensg2entrez = {}
    ensg2symbol = {}
    for ensg, sub in df.groupby("gene_stable_id"):
        eset = {norm_entrez(v) for v in sub["entrez_id"]}
        eset.discard(None)
        ensg2entrez[ensg] = eset
        syms = {norm_symbol(v) for v in sub["hgnc_symbol"]}
        syms.discard(None)
        ensg2symbol[ensg] = syms
    return ensg2entrez, ensg2symbol


def load_tcpa():
    """RPPA antibody (peptide id) -> authoritative (entrez, symbol)."""
    path = os.path.join(MAPPINGS, "TCPA_RPPA500_metadata_mapping.tsv")
    df = pd.read_csv(path, sep="\t", dtype=str)
    by_id = {}
    for _, row in df.iterrows():
        pid = str(row["RPPA_Protein_ID"]).strip().upper()
        by_id[pid] = {
            "entrez": norm_entrez(row.get("Entrez_Gene_ID")),
            "symbol": norm_symbol(row.get("Gene_Symbol")),
        }
    return by_id


def load_mirna_ref():
    """miRNA accession (hsa-...) -> authoritative HGNC symbol (mirna_hgcn_map)."""
    path = os.path.join(MAPPINGS, "mirna_hgcn_map.tsv")
    df = pd.read_csv(path, sep="\t", dtype=str)
    by_acc = {}
    for _, row in df.iterrows():
        acc = str(row["miRNA_ID"]).strip().lower()
        by_acc[acc] = norm_symbol(row["hgcn_id"])
    return by_acc


def load_hm27_manifest():
    """CpG probe id -> authoritative set of gene symbol(s) from the Illumina
    HM27 manifest (GDC Xena mirror).

    The manifest stores multiple genes for a probe as a SINGLE comma-separated
    cell (e.g. "AC092720.2,JPH3"); we split it into a set so it can be compared
    against the index, which stores the same genes as a list."""
    path = os.path.join(MAPPINGS, "HM27.hg38.manifest.gencode.v36.probeMap")
    df = pd.read_csv(path, sep="\t", dtype=str)
    # columns: #id, gene, chrom, chromStart, chromEnd, strand
    id_col = df.columns[0]
    by_probe = {}
    for _, row in df.iterrows():
        probe = str(row[id_col]).strip()
        cell = row["gene"]
        by_probe.setdefault(probe, set())
        if cell is not None and str(cell).strip().lower() not in {"", "nan", "none"}:
            for part in str(cell).split(","):
                gene = norm_symbol(part)
                if gene:
                    by_probe[probe].add(gene)
    return by_probe


# ===========================================================================
# Per-layer audits
# ===========================================================================
def sample_positions(mappings, n, rng):
    """Deterministically sample up to n mappings."""
    if len(mappings) <= n:
        return list(mappings)
    idx = rng.sample(range(len(mappings)), n)
    return [mappings[i] for i in idx]


def audit_gene_layer(db, mappings, n, rng, ensg2entrez, layer_name):
    """RNA-Seq / CNV: check index entrez_id for each ENSG against BioMart, and
    check the EntrezID backbone node exists."""
    sample = sample_positions(mappings, n, rng)
    entrez_ids = [norm_entrez(m.get("entrez_id")) for m in sample]
    existing_nodes = entrez_nodes_exist(db, entrez_ids)

    rows = []
    for m in sample:
        ensg = m.get("gene_id_base")
        idx_entrez = norm_entrez(m.get("entrez_id"))
        truth = ensg2entrez.get(ensg, set())
        node_ok = idx_entrez in existing_nodes if idx_entrez else False

        if not truth:
            # ENSG absent from BioMart ground truth -> cannot adjudicate; treat
            # as "unverifiable", excluded from N and error count (logged).
            verdict, note = "unverifiable", "ENSG not in BioMart reference"
        elif idx_entrez is None:
            verdict, note = "error", "no entrez_id in index (BioMart has one)"
        elif idx_entrez in truth:
            verdict = "ok"
            note = "" if node_ok else "xref ok but EntrezID node missing"
        else:
            verdict = "error"
            note = f"index entrez {idx_entrez} not in BioMart {sorted(truth)}"

        rows.append({
            "layer": layer_name, "feature": ensg, "index_value": idx_entrez,
            "authoritative": ";".join(sorted(truth)) if truth else "",
            "node_exists": node_ok, "verdict": verdict, "note": note,
        })
    return rows


def audit_protein_layer(db, mappings, n, rng, tcpa):
    sample = sample_positions(mappings, n, rng)
    entrez_ids = [norm_entrez(m.get("entrez_id")) for m in sample]
    existing_nodes = entrez_nodes_exist(db, entrez_ids)

    rows = []
    for m in sample:
        pid = str(m.get("peptide_target", "")).strip().upper()
        idx_entrez = norm_entrez(m.get("entrez_id"))
        idx_symbol = norm_symbol(m.get("gene_symbol"))
        truth = tcpa.get(pid)
        node_ok = idx_entrez in existing_nodes if idx_entrez else False

        if truth is None:
            verdict, note = "unverifiable", "antibody not in TCPA reference"
        else:
            t_entrez, t_symbol = truth["entrez"], truth["symbol"]
            entrez_match = (idx_entrez == t_entrez) if t_entrez else None
            symbol_match = (idx_symbol == t_symbol) if t_symbol else None
            # match if either authoritative key agrees (RPPA antibodies map to a
            # gene product; entrez and symbol are redundant keys for the gene).
            checks = [c for c in (entrez_match, symbol_match) if c is not None]
            if checks and any(checks):
                verdict = "ok"
                note = "" if node_ok else "xref ok but EntrezID node missing"
            else:
                verdict = "error"
                note = (f"index ({idx_entrez},{idx_symbol}) vs "
                        f"TCPA ({t_entrez},{t_symbol})")
        rows.append({
            "layer": "RPPA", "feature": pid, "index_value": f"{idx_entrez}/{idx_symbol}",
            "authoritative": (f"{truth['entrez']}/{truth['symbol']}"
                              if truth else ""),
            "node_exists": node_ok, "verdict": verdict, "note": note,
        })
    return rows


def audit_mirna_layer(mappings, n, rng, mirna_ref):
    sample = sample_positions(mappings, n, rng)
    rows = []
    for m in sample:
        acc = str(m.get("mirna_id", "")).strip().lower()
        idx_symbol = norm_symbol(m.get("hgnc_symbol"))
        truth = mirna_ref.get(acc)

        if truth is None:
            verdict, note = "unverifiable", "miRNA accession not in HGNC reference"
        else:
            # KG drops '-' from symbols; compare on both raw and dash-stripped.
            idx_variants = {idx_symbol, (idx_symbol or "").replace("-", "")}
            truth_variants = {truth, truth.replace("-", "")}
            if idx_variants & truth_variants:
                verdict, note = "ok", ""
            else:
                verdict = "error"
                note = f"index {idx_symbol} vs miRBase/HGNC {truth}"
        rows.append({
            "layer": "miRNA", "feature": acc, "index_value": idx_symbol,
            "authoritative": truth or "", "node_exists": "",
            "verdict": verdict, "note": note,
        })
    return rows


def audit_methylation_layer(mappings, n, rng, hm27):
    sample = sample_positions(mappings, n, rng)
    rows = []
    for m in sample:
        probe = str(m.get("probe_id", "")).strip()
        idx_genes = {norm_symbol(g) for g in (m.get("gene_symbols") or [])}
        idx_genes.discard(None)
        truth = hm27.get(probe)

        if truth is None:
            verdict, note = "unverifiable", "probe not in HM27 manifest"
        elif not truth:
            # manifest lists probe but with no gene (intergenic) — match if index
            # also has no gene.
            verdict = "ok" if not idx_genes else "error"
            note = "" if not idx_genes else "manifest intergenic, index has gene"
        elif idx_genes & truth:
            verdict, note = "ok", ""
        elif not idx_genes:
            verdict, note = "error", f"index has no gene, manifest {sorted(truth)}"
        else:
            verdict = "error"
            note = f"index {sorted(idx_genes)} vs manifest {sorted(truth)}"
        rows.append({
            "layer": "Methylation", "feature": probe,
            "index_value": ";".join(sorted(idx_genes)),
            "authoritative": ";".join(sorted(truth)) if truth else "",
            "node_exists": "", "verdict": verdict, "note": note,
        })
    return rows


# ===========================================================================
# Driver
# ===========================================================================
LAYER_LABELS = {
    "RNA-Seq": "RNA-Seq (Ensembl/Entrez)",
    "CNV": "CNV (Entrez)",
    "miRNA": "miRNA (miRBase)",
    "Methylation": "Methylation (CpG/gene)",
    "RPPA": "RPPA (Entrez/HGNC)",
}
LAYER_ORDER = ["RNA-Seq", "CNV", "miRNA", "Methylation", "RPPA"]


def summarise(all_rows):
    """Build the per-layer + Overall summary.

    N_audited counts only VERIFIABLE mappings (verdict ok/error) — i.e. features
    that actually carry a mapping to an authoritative external identifier. Rows
    marked 'unverifiable' are features with no authoritative cross-reference to
    audit (e.g. non-coding ENSGs with no Entrez ID, or antibodies absent from
    the TCPA reference); they are a COMPLETENESS statistic, not accuracy, and are
    reported separately as n_sampled / coverage."""
    summary = []
    tot_n = tot_err = tot_unv = 0
    for layer in LAYER_ORDER:
        rows = [r for r in all_rows if r["layer"] == layer]
        verifiable = [r for r in rows if r["verdict"] in {"ok", "error"}]
        errors = [r for r in verifiable if r["verdict"] == "error"]
        unver = [r for r in rows if r["verdict"] == "unverifiable"]
        n_sampled = len(rows)
        n = len(verifiable)
        e = len(errors)
        rate = (100.0 * e / n) if n else 0.0
        coverage = (100.0 * n / n_sampled) if n_sampled else 0.0
        examples = "; ".join(
            f"{r['feature']}: {r['note']}" for r in errors[:2]
        ) or "no errors"
        summary.append({
            "layer": layer, "label": LAYER_LABELS[layer],
            "n_sampled": n_sampled, "n_audited": n, "n_errors": e,
            "error_rate_pct": round(rate, 2),
            "n_unverifiable": len(unver),
            "coverage_pct": round(coverage, 1),
            "examples": examples,
        })
        tot_n += n
        tot_err += e
        tot_unv += len(unver)
    tot_sampled = tot_n + tot_unv
    overall_rate = (100.0 * tot_err / tot_n) if tot_n else 0.0
    overall_cov = (100.0 * tot_n / tot_sampled) if tot_sampled else 0.0
    summary.append({
        "layer": "Overall", "label": "Overall",
        "n_sampled": tot_sampled, "n_audited": tot_n, "n_errors": tot_err,
        "error_rate_pct": round(overall_rate, 2),
        "n_unverifiable": tot_unv,
        "coverage_pct": round(overall_cov, 1),
        "examples": "---",
    })
    return summary


def print_latex(summary):
    print("\n% ---- Supplementary Table S3 (filled) ----")
    print(r"\begin{tabular}{@{}lrrrl@{}}")
    print(r"\toprule")
    print(r"\textbf{Layer} & \textbf{\shortstack{$N$\\audited}} & "
          r"\textbf{\shortstack{$N$\\errors}} & "
          r"\textbf{\shortstack{Error\\rate}} & "
          r"\textbf{Notes / example errors} \\ \midrule")
    for s in summary:
        if s["layer"] == "Overall":
            print(r"\midrule")
            print(rf"\textbf{{Overall}} & {s['n_audited']} & {s['n_errors']} & "
                  rf"{s['error_rate_pct']:.1f}\% & --- \\ \bottomrule")
            continue
        if s["n_errors"] == 0:
            note = "no errors"
        else:
            note = s["examples"]
        note = note.replace("_", r"\_").replace("%", r"\%").replace("&", r"\&")
        if len(note) > 60:
            note = note[:57] + "..."
        print(rf"{s['label']} & {s['n_audited']} & {s['n_errors']} & "
              rf"{s['error_rate_pct']:.1f}\% & {note} \\")
    print(r"\end{tabular}")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--n", type=int, default=100, help="audited mappings per layer")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--db", default="PKT_main")
    ap.add_argument("--cohort", default="TCGA-BRCA")
    ap.add_argument("--host", default="http://localhost:8529")
    ap.add_argument("--user", default="root")
    ap.add_argument("--password", default="avocadodb")
    ap.add_argument("--out-dir", default=HERE)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    print(f"[D4 spot-check] db={args.db} cohort={args.cohort} "
          f"N={args.n}/layer seed={args.seed}")

    db = connect(args.db, args.host, args.user, args.password)

    print("[1/3] Loading external ground-truth tables...")
    ensg2entrez, _ = load_biomart()
    tcpa = load_tcpa()
    mirna_ref = load_mirna_ref()
    hm27 = load_hm27_manifest()
    print(f"      BioMart ENSGs: {len(ensg2entrez)} | TCPA antibodies: {len(tcpa)} | "
          f"miRNA refs: {len(mirna_ref)} | HM27 probes: {len(hm27)}")

    co = args.cohort
    layer_specs = [
        ("RNA-Seq", "GENE_EXPRESSION", f"expr_index_{co}", "gene_mappings"),
        ("CNV", "CNV", f"cnv_index_{co}", "gene_mappings"),
        ("miRNA", "MIRNA", f"mirna_index_{co}", "mirna_mappings"),
        ("Methylation", "METHYLATION", f"methylation_index_{co}", "probe_mappings"),
        ("RPPA", "PROTEIN", f"protein_index_{co}", "protein_mappings"),
    ]

    print("[2/3] Auditing layers...")
    all_rows = []
    for layer, coll, key, field in layer_specs:
        maps = load_index_mappings(db, coll, key, field)
        if layer in ("RNA-Seq", "CNV"):
            rows = audit_gene_layer(db, maps, args.n, rng, ensg2entrez, layer)
        elif layer == "RPPA":
            rows = audit_protein_layer(db, maps, args.n, rng, tcpa)
        elif layer == "miRNA":
            rows = audit_mirna_layer(maps, args.n, rng, mirna_ref)
        else:
            rows = audit_methylation_layer(maps, args.n, rng, hm27)
        all_rows.extend(rows)
        nv = sum(1 for r in rows if r["verdict"] in {"ok", "error"})
        ne = sum(1 for r in rows if r["verdict"] == "error")
        print(f"      {layer:12s}: {len(maps):6d} mappings, "
              f"{nv} verifiable, {ne} errors")

    print("[3/3] Writing results...")
    os.makedirs(args.out_dir, exist_ok=True)
    detail_path = os.path.join(args.out_dir, "spotcheck_results.csv")
    with open(detail_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["layer", "feature", "index_value",
                                          "authoritative", "node_exists",
                                          "verdict", "note"])
        w.writeheader()
        w.writerows(all_rows)

    summary = summarise(all_rows)
    summary_path = os.path.join(args.out_dir, "spotcheck_summary.csv")
    with open(summary_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=["layer", "label", "n_sampled",
                                          "n_audited", "n_errors",
                                          "error_rate_pct", "n_unverifiable",
                                          "coverage_pct", "examples"])
        w.writeheader()
        w.writerows(summary)

    print("\n=== SUMMARY (seed={}, N={}/layer sampled) ===".format(args.seed, args.n))
    print(f"{'Layer':<26}{'Smpl':>5}{'Aud':>5}{'Err':>5}{'Rate%':>8}{'Cov%':>7}")
    for s in summary:
        print(f"{s['label']:<26}{s['n_sampled']:>5}{s['n_audited']:>5}"
              f"{s['n_errors']:>5}{s['error_rate_pct']:>8.2f}"
              f"{s['coverage_pct']:>7.1f}")

    print_latex(summary)
    print(f"\nWrote: {detail_path}")
    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
