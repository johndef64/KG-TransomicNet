# `data/mappings/` — Identifier-mapping tables

This folder hosts the controlled mapping tables used by the build pipeline to
join the **semantic layer** (PheKnowLator backbone) with the **quantitative
layer** (TCGA / TARGET omic matrices) by gene, transcript, miRNA, CpG probe,
RPPA antibody, dbSNP rsID and MONDO disease identifiers.

All consumers go through `scripts/build_omics_collections.py:load_mappings()`,
which automatically:

* unzips a `*.tsv.zip` archive on first use if the `.tsv` is not present;
* downloads probemap files from the [GDC Xena hub](https://gdc-hub.s3.us-east-1.amazonaws.com/download/) if they are missing.

So the typical user only needs `git pull` + run the build script.

---

## File inventory

Legend:
* **shipped (.tsv)** — small enough to commit as plain text
* **shipped (.zip)** — committed compressed; auto-extracted on first use
* **auto-download** — fetched from the GDC hub on first use (not committed)

| File | Size (raw) | Ship mode | Produced by | Used by |
|---|---:|---|---|---|
| `biomart_gene_mappings_filled_extended.tsv` | 36 MB | shipped (.zip → 4.7 MB) | `scripts/mapping/biomart_extender.py` | `build_omics_collections` semantic+index layers |
| `biomart_ensg_enst_mapping.tsv` | 110 MB | shipped (.zip → 12.6 MB) | `scripts/mapping/biomart_maps_ensg_enst_mapping.py` | `build_omics_collections` CNV index |
| `rsid_checkpoint.tsv` | 5.1 MB | shipped (.zip → 2.1 MB) | `scripts/mapping/dbsnp_maps.py` | (annotation of dbSNP nodes; future use) |
| `mirna_hgcn_map.tsv` | 40 KB | shipped (.zip → 10 KB) | `scripts/omics_utils.py::make_mirna_hgcn_map()` | `build_omics_collections` miRNA index |
| `TCPA_RPPA500_metadata_mapping.tsv` | 40 KB | shipped (.tsv) | manual download from [TCPA portal](https://tcpaportal.org/tcpa/) + manual cleanup | `build_omics_collections` protein index |
| `mondo_to_primarydisease.tsv` | <1 KB | shipped (.tsv) | `scripts/mapping/mondo_mapping_manual.py` (semi-manual) | `build_omics_collections` PROJECT collection |
| `HM27.hg38.manifest.gencode.v36.probeMap` | 1.3 MB | auto-download | UCSC Xena (GDC hub mirror) | `build_omics_collections` methylation index |
| `HM450.hg38.manifest.gencode.v36.probeMap` | 22 MB | auto-download | UCSC Xena (GDC hub mirror) | `build_omics_collections` methylation index (HM450 cohorts) |
| `gencode.v36.annotation.gtf.gene.probemap` | ~18 MB | auto-download | UCSC Xena (GDC hub mirror) | reserved for future use |
| `monarch/` | various | shipped (.tsv.gz) | downloaded from [Monarch Initiative](https://data.monarchinitiative.org/) | future enrichment of gene–disease and gene–phenotype edges |
| `pkt/` | small | shipped | `scripts/kg_utils/pkt_nodelabel-processing.py` + manual curation | offline reference, PKT node-type catalogue |

The `backups/` and `dev/` subfolders contain previous versions / scratch
artifacts; they are **not** read by any pipeline script and are kept only for
provenance.

---

## Provenance details

### `biomart_gene_mappings_filled_extended.tsv` — the core BioMart table

This is the canonical gene-centric join table used by every quantitative
layer. It has 15 columns: `gene_stable_id`, `gene_stable_id_version`,
`gene_type`, `gene_description`, `chromosome`, `gene_start_bp`, `gene_end_bp`,
`strand`, `hgnc_symbol`, `mirbase_id`, `uniprot_swissprot_id`,
`uniprot_isoform_id`, `entrez_id`, `gene_id`, `transcript_id` (Python literal
list).

It is rebuilt end-to-end by
[`scripts/mapping/biomart_extender.py`](../../scripts/mapping/biomart_extender.py),
which consolidates four historical steps into one CLI:

1. **Initial harvest** (input) — `data/mappings/biomart_gene_mappings.tsv`,
   produced by [`biomart_maps.py`](../../scripts/mapping/biomart_maps.py)
   querying BioMart chromosome by chromosome (`hsapiens_gene_ensembl`,
   Ensembl v109 / GRCh38 / GENCODE v36).

2. **HGNC fill** — missing `hgnc_symbol` values are filled from the GENCODE
   gene probemap (auto-downloaded by the build pipeline).

3. **Entrez fill (two passes)**:
   * a fast self-lookup using rows where both `hgnc_symbol` and `entrez_id`
     are present;
   * an optional NCBI Entrez fallback for the remaining orphans, with on-disk
     checkpoint so the run can be resumed (slow: ~30-60 min, hits NCBI's
     3 req/s limit).

4. **Transcript merge** — joins with `biomart_ensg_enst_mapping.tsv`
   (from [`biomart_maps_ensg_enst_mapping.py`](../../scripts/mapping/biomart_maps_ensg_enst_mapping.py))
   to add `gene_id` and `transcript_id` columns, producing the final 15-column
   table.

**Reproducing it locally:**

```bash
# Quick run (skip NCBI; ~30 s, slightly lower entrez_id coverage)
python scripts/mapping/biomart_extender.py --skip-ncbi

# Full run with NCBI Entrez fallback (recommended; ~30-60 min)
python scripts/mapping/biomart_extender.py --ncbi-email you@example.com
```

The shipped `.zip` is the output of a full run we performed for the paper.
Most users will never need to rerun this script.

### `biomart_ensg_enst_mapping.tsv` — gene → transcript mapping

Produced by [`scripts/mapping/biomart_maps_ensg_enst_mapping.py`](../../scripts/mapping/biomart_maps_ensg_enst_mapping.py)
(pybiomart, queries `hsapiens_gene_ensembl` per chromosome). Two columns:
`gene_id`, `transcript_id`. Used by `build_omics_collections` to attach the
list of ENSTs to each ENSG in the CNV index. Long but fully scripted; could
be regenerated in ~5 minutes if BioMart is reachable.

### `rsid_checkpoint.tsv` — dbSNP variant annotations

Produced by [`scripts/mapping/dbsnp_maps.py`](../../scripts/mapping/dbsnp_maps.py).
For every rsID node in the PKT NodeLabels (~50k rsIDs), it queries
`clinicaltables.nlm.nih.gov/api/snps/v3/search` to obtain chromosome,
position, alleles and gene. Uses `ThreadPoolExecutor` and writes a
**checkpoint after every batch** so the run can resume after interruption.

**Reproducibility note.** The NLM API is rate-limited, the full run takes
many hours, and individual queries occasionally fail; the shipped checkpoint
captures the successful subset that was used for the paper.

### `mirna_hgcn_map.tsv` — miRBase ↔ HGNC

Produced by [`scripts/omics_utils.py::make_mirna_hgcn_map()`](../../scripts/omics_utils.py).
A trivial regex on the `miRNA_ID` column of any TCGA miRNA matrix:
`hsa-mir-XXX` → `MIRXXX`, `hsa-let-XXX` → `MIRLETXXX`, uppercased.
Could be re-derived in 1 second, but is shipped for convenience.

### `mondo_to_primarydisease.tsv` — TCGA disease → MONDO

Produced by [`scripts/mapping/mondo_mapping_manual.py`](../../scripts/mapping/mondo_mapping_manual.py).
33 rows, one per TCGA study disease type. The MONDO IDs were assigned
manually after a one-to-one inspection of MONDO labels against the TCGA
`disease_type.project` values. Plain text, ships uncompressed.

### `TCPA_RPPA500_metadata_mapping.tsv` — RPPA antibody → gene

Downloaded manually from the [TCPA portal antibody metadata](https://tcpaportal.org/tcpa/)
and lightly post-processed to expose four useful columns
(`RPPA_Protein_ID`, `Entrez_Gene_ID`, `Gene_Symbol`, `Protein_Type`).
No script reproduces it — it is treated as a third-party reference table.

### `HM27` / `HM450` / `gencode` probemaps

Public GDC hub assets, mirrored on the same S3 endpoint TCGA omic matrices
are served from. They are **never committed** (gitignored as `*.probemap`)
and are downloaded transparently on first call to `load_mappings()`. URLs are
defined in `scripts/download_omics.py:PROBEMAP_URLS`.

### `monarch/` — Monarch Initiative associations

Static dumps from `data.monarchinitiative.org`, downloaded with
[`data/mappings/monarch/read_monarch.py`](monarch/read_monarch.py).
Currently not consumed by the build pipeline — kept for future enrichment
of gene↔disease, gene↔phenotype and gene↔pathway edges.

### `pkt/` — PKT node-type reference

Offline reference materials (PDFs, node-type CSVs). The complete CSV was
built by `scripts/kg_utils/pkt_nodelabel-processing.py` (now folded into
`scripts/download_pkt.py`). Used only for documentation / debugging.

---

## How `load_mappings()` resolves each file

For every required mapping `fname` under `data/mappings/`:

1. If `data/mappings/fname` exists → use it.
2. Else if `data/mappings/fname.zip` exists → unzip in place, use the result.
3. Else if `fname` is a known GDC probemap → download from
   `https://gdc-hub.s3.us-east-1.amazonaws.com/download/<fname>`.
4. Else → log a WARNING and continue (the missing mapping is treated as
   optional unless `biomart_genes` is missing, in which case the build
   aborts with a clear error message pointing back to this README).

---

## Decision matrix — committed `.zip` vs scripted rebuild

For the paper data-availability statement we keep the **`.zip`-with-shipping**
strategy for the four heavy/curated files. Rationale:

* `biomart_gene_mappings_filled_extended.tsv` can now be fully regenerated by
  `biomart_extender.py`, but the run takes ~30-60 minutes (NCBI rate limit)
  and the .zip preserves the exact snapshot used in the paper.
* `rsid_checkpoint.tsv` requires multi-hour API scraping that's fragile.
* `biomart_ensg_enst_mapping.tsv` could be regenerated, but reproducibility
  is better served by pinning the exact snapshot we used than by relying on
  BioMart serving identical results months later.
* `mirna_hgcn_map.tsv` is trivial to regenerate but is shipped for symmetry.

The `scripts/mapping/` and `data/mappings/dev/` folders document the
provenance for transparency and reviewer scrutiny; users do not need to run
them.

---

## What ships under `scripts/mapping/`

Only the scripts that are **runnable in isolation** and that document the
provenance of a shipped mapping are committed to the repository:

| Script | Status | Reproduces |
|---|---|---|
| [`biomart_maps.py`](../../scripts/mapping/biomart_maps.py) | exploratory but functional | `biomart_gene_mappings.tsv` (per-chromosome BioMart harvest, slow) |
| [`biomart_maps_ensg_enst_mapping.py`](../../scripts/mapping/biomart_maps_ensg_enst_mapping.py) | clean, CLI `main()` | `biomart_ensg_enst_mapping.tsv` |
| [`biomart_extender.py`](../../scripts/mapping/biomart_extender.py) | clean, CLI with `--skip-ncbi` | `biomart_gene_mappings_filled_extended.tsv` (full end-to-end pipeline) |
| [`dbsnp_maps.py`](../../scripts/mapping/dbsnp_maps.py) | clean, CLI `main()` with checkpoint resume | `rsid_checkpoint.tsv` |
| [`mondo_mapping_manual.py`](../../scripts/mapping/mondo_mapping_manual.py) | notebook-style, kept for documentation only | (manual assignments behind `mondo_to_primarydisease.tsv`) |

The following scripts are part of the historical provenance but are **not
committed** (listed in `.gitignore`); they have been **superseded by
`biomart_extender.py`** and are kept locally only for archival purposes:

| Script (local only) | Superseded by |
|---|---|
| `biomart_fixer.py` | `biomart_extender.py` step 1-2 |
| `biomart_fixer_v2.py` | `biomart_extender.py` step 1-3 |
| `biomart_fixer_v2_part2.py` | `biomart_extender.py` step 4 |
| `create_mappings_collection.py` | (scratch fragment, no replacement needed) |
