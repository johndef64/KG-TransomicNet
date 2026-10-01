# PanCanAtlas view (view P)

The quantitative layer of KG-TransomicNet holds omics in the *vector + index*
model, one cohort at a time. The per-study GDC data (STAR-TPM, ASCAT3, miRNA-seq,
RPPA, HM27) form the **GDC view**. The **PanCanAtlas view** adds the harmonised
pan-cancer release of TCGA (PanCanAtlas 2018) as one more cohort, `PANCAN-ATLAS`,
in the same collections and with the same document schema.

Why a second view: the PanCanAtlas matrices are batch-corrected across cancer
types (EB++), and the integrated iCluster subtypes of Hoadley et al. (2018) were
derived on exactly these data. Pan-cancer tasks that use those labels therefore
need these values, not the per-study GDC ones.

The view is purely additive. It writes its own folder and its own cohort key and
never modifies GDC files or GDC documents.

## Sources

| Layer | File (Xena) | Hub |
|---|---|---|
| RNA-seq, EB++ adjusted, log2(norm_value+1) | `EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena.gz` | tcga-pancan-atlas-hub |
| miRNA, EB adjusted, mature forms | `pancanMiRs_EBadjOnProtocolPlatformWithoutRepsWithUnCorrectMiRs_08_04_16.xena.gz` | tcga-pancan-atlas-hub |
| CNV, GISTIC2 gene-level, continuous | `TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz` | tcga-xena-hub |
| iCluster k=28 assignments (27 classes) | `TCGA_PanCan33_iCluster_k28_tumor.gz` | tcga-pancan-atlas-hub |

## Commands

```bash
# 1. Download (about 810 MB) into data/omics/PANCAN-ATLAS/, with a SHA-256 manifest
python scripts/pancan_atlas/download_pancan_atlas.py
python scripts/pancan_atlas/download_pancan_atlas.py --from-dir /path/with/copies   # reuse local copies

# 2. Build the JSON collections into data/arangodb_collections/PANCAN-ATLAS/
python scripts/pancan_atlas/build_pancan_atlas_collections.py --limit-samples 50 --verify   # quick test
python scripts/pancan_atlas/build_pancan_atlas_collections.py                               # full, about 3.4 GB

# 3. (Optional) Load with the standard loader, into a database of your choice
python scripts/load_omics_collections_to_arangodb.py --studies PANCAN-ATLAS \
       --layers gene_expression cnv mirna --no-semantic --db <DB>
```

Never pass `--replace-existing` when loading view P: that flag drops whole
collections, GDC documents included.

## Documents

| Layer | Index `_key` | Vector value field | Notes |
|---|---|---|---|
| gene_expression | `expr_index_PANCAN-ATLAS` | `values_ebpp` | read with `value_type="ebpp"` |
| cnv | `cnv_index_PANCAN-ATLAS` | `values_copy_number` | `analysis_method: "GISTIC2"`, `value_scale` gives the scale |
| mirna | `mirna_index_PANCAN-ATLAS` | `values_expression` | mature miRNA ids |

- Vector `_key` is the 15-character PanCanAtlas barcode (`TCGA-XX-XXXX-01`);
  GDC vectors use 16 characters, so the two views never collide.
- `samples.json` holds one document per sample with its layers and its
  `icluster_k28` label; `projects.json` holds the `PANCAN-ATLAS` cohort.
- Every document carries `source_release: "PanCanAtlas 2018"`.
- The existing query helpers in `query_utils.py` work unchanged with
  `cohort="PANCAN-ATLAS"`.

Gene identifiers are resolved with `data/mappings/hgcn_mart_conversion_table.zip`
(HGNC approved symbol → NCBI gene ID). Obsolete symbols and the few bare
Entrez rows of the expression matrix stay unresolved. A mature miRNA is mapped
only when it belongs to a single precursor.

## Statistics of the view

Computed from the full matrices without loading them. Coverage is the share of
features that resolve to a gene node of the backbone (780,753 nodes). "In GDC
view" counts samples whose 15-character barcode prefix is also loaded in the
GDC view (16,938 distinct samples).

| Layer | Samples | Patients | Features | Resolved to Entrez | In backbone | Coverage | In GDC view |
|---|---:|---:|---:|---:|---:|---:|---:|
| RNA-seq (EB++) | 11,069 | 10,274 | 20,531 | 17,133 | 16,988 | 0.827 | 11,063 |
| CNV (GISTIC2) | 10,845 | 10,845 | 24,776 | 21,470 | 19,520 | 0.788 | 10,785 |
| miRNA (EB) | 10,824 | 10,113 | 743 | 641 | 637 | 0.857 | 10,824 |

Samples with all three layers and an iCluster label: **9,355** (27 classes,
31 TCGA projects), all of them also present in the GDC view.
