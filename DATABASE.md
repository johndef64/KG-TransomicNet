# KG-TransomicNet database structure

Schema of the materialised ArangoDB instance: collections, document fields and the
query pattern that joins the quantitative layers to the knowledge graph. Field lists
and examples are taken from the published instance.

- **Database name:** `PKT_main` by default; every script accepts `--db` (see the main
  [readme](readme.md#database-name))
- **Graph model:** property graph (document collections plus one edge collection)

## Collections at a glance

| Collection | Layer | Documents | Content |
|---|---|---:|---|
| `nodes` | semantic | 780,753 | ontologically grounded biological entities (PheKnowLator v3.0.2) |
| `edges` | semantic | 11,082,103 | typed relations between `nodes` |
| `GENES` | shared metadata | 61,855 | Ensembl gene records shared by all cohorts |
| `PROJECTS` | metadata | 42 | TCGA and TARGET projects |
| `CASES` | metadata | 16,086 | patients: demographics, diagnoses, survival |
| `SAMPLES` | metadata | 19,287 | biospecimens (barcode, type, case, project) |
| `GENE_EXPRESSION` | quantitative | 15,475 | expression vectors + 42 cohort indexes |
| `CNV` | quantitative | 10,665 | copy-number vectors + 33 cohort indexes |
| `MIRNA` | quantitative | 13,441 | miRNA vectors + 38 cohort indexes |
| `PROTEIN` | quantitative | 7,936 | RPPA vectors + 32 cohort indexes |
| `METHYLATION` | quantitative | 3,150 | methylation vectors + 13 cohort indexes |

Total: 12,010,793 documents.

## Semantic layer

### `nodes`

| Field | Type | Description |
|---|---|---|
| `_key` | string | node identifier, equal to `entity_id` (e.g. `80724`, `PR_Q6JQN1`, `MONDO_0021282`) |
| `uri` | string | canonical URI of the entity |
| `namespace` | string | URI host (e.g. `www.ncbi.nlm.nih.gov`, `purl.obolibrary.org`) |
| `entity_id` | string | identifier in the source namespace |
| `class_code` | string | ontology or database prefix (e.g. `EntrezID`, `PR`, `GO`, `MONDO`, `HP`, `R-HSA`, `dbSNP`, `ENST`, `CHEBI`) |
| `bioentity_type` | string | coarse entity type (see below) |
| `label` | string | human-readable label (e.g. `ACAD10 (human)`) |
| `description` | string | textual definition |
| `synonym` | string | alternative names |
| `source` | string | source database or ontology (e.g. `NCBI Entrez Gene`) |
| `source_type` | string | `Database` or `Ontology` |
| `integer_id` | integer | internal numeric identifier |

`bioentity_type` takes 14 values: `gene`, `protein`, `rna`, `variant`, `chemical`,
`disease`, `phenotype`, `go`, `pathway`, `anatomy`, `cell`, `catalyst`, `organism`,
`unknown`.

There are no dedicated mapping fields on the nodes: an omic feature is joined to a node
through `class_code` and `entity_id`.

| Target entity | `class_code` | `entity_id` format | Example |
|---|---|---|---|
| gene | `EntrezID` | NCBI gene ID | `7157` |
| protein | `PR` | `PR_` + UniProt accession | `PR_P04637` |
| transcript | `ENST` | Ensembl transcript ID | `ENST00000363558` |
| variant | `dbSNP` | rsID | `rs201867256` |
| disease | `MONDO` | `MONDO_` + number | `MONDO_0021282` |
| phenotype | `HP` | `HP_` + number | `HP_0025607` |
| GO term | `GO` | `GO_` + number | `GO_0005737` |
| pathway | `R-HSA` | Reactome stable ID | `R-HSA-162582` |

### `edges`

| Field | Type | Description |
|---|---|---|
| `_key` | string | edge identifier (e.g. `edge_0`) |
| `_from`, `_to` | string | `nodes/<_key>` of source and target |
| `source_uri`, `target_uri` | string | URIs of source and target |
| `predicate_uri` | string | relation URI (e.g. `http://purl.obolibrary.org/obo/RO_0002205`) |
| `predicate_label` | string | relation label (e.g. `has gene product`, `participates in`, `causes or contributes to condition`) |
| `predicate_class_code` | string | ontology of the predicate (e.g. `RO`) |
| `predicate_source` | string | source ontology name (e.g. `Relation Ontology`) |
| `predicate_bioentity_type` | string | semantic type attached to the predicate |

## Metadata collections

### `PROJECTS`

`_key` (project id, e.g. `TCGA-BRCA`), `name`, `program` (`TCGA` or `TARGET`),
`primary_site`, `disease_types` (array), `n_cases`, `n_samples`, `entity_type`
(`project`).

### `CASES`

`_key` (case submitter id, e.g. `TCGA-3C-AAAU`), `project_ref` (`projects/<id>`),
`disease_type`, `primary_site`, `entity_type` (`case`), and three nested objects:

- `demographic`: `gender`, `race`, `ethnicity`, `vital_status`, `age_at_index`,
  `days_to_birth`, `year_of_birth`, `year_of_death`, `days_to_death`
- `diagnoses`: diagnosis records when available (may be null)
- `survival`: `overall_survival_time`, `overall_survival_status`

### `SAMPLES`

`_key` (sample barcode, e.g. `TCGA-BH-A0W3-01A`), `submitter_id` (case id),
`case_ref` (`cases/<id>`), `project_ref` (`projects/<id>`), `sample_type`
(e.g. `Primary Tumor`), `sample_type_id`, `tissue_type`, `tumor_descriptor`,
`specimen_type`, `composition`, `preservation_method`, `days_to_collection`,
`entity_type` (`sample`).

### `GENES`

`_key` (Ensembl gene ID without version, e.g. `ENSG00000141510`), `gene_stable_id`,
`gene_stable_id_version`, `hgnc_symbol`, `entrez_id`, `uniprot_id`, `mirbase_id`,
`gene_type`, `gene_description`, `chromosome`, `gene_start_bp`, `gene_end_bp`,
`strand`, `transcript_ids` (array), `source` (`Ensembl_BioMart`), `source_version`
(`GRCh38.v36`), `entity_type` and `bioentity_type` (`gene`).

Join to the knowledge graph: `GENES.entrez_id` = `nodes._key` of a node with
`class_code == "EntrezID"`.

## Quantitative layer: vector + index

Each omic collection holds two document types, distinguished by `data_type`:

1. **vector** documents (`*_vector`), one per sample, with the measurements as dense
   arrays (`values_*`);
2. **index** documents (`*_index`), one per cohort, mapping each array position to a
   biological feature.

A vector document points to its index through `*_index_ref`; index keys follow the
pattern `<prefix>_<cohort>`.

| Collection | Index `_key` | Mapping field | Vector value fields | Platform / method |
|---|---|---|---|---|
| `GENE_EXPRESSION` | `expr_index_<cohort>` | `gene_mappings` | `values_tpm`, `values_fpkm`, `values_counts` | STAR, 60,660 genes |
| `CNV` | `cnv_index_<cohort>` | `gene_mappings` | `values_copy_number` | ASCAT3 gene-level, 60,623 genes |
| `MIRNA` | `mirna_index_<cohort>` | `mirna_mappings` | `values_expression` | miRNA-seq, 1,881 miRNAs |
| `PROTEIN` | `protein_index_<cohort>` | `protein_mappings` | `values_abundance` | RPPA, 487 targets |
| `METHYLATION` | `methylation_index_<cohort>` | `probe_mappings` | `values_beta` | Illumina HM27, 27,578 probes |

### Vector documents

Common fields: `_key` and `sample_id` (sample barcode), `cohort` (project id),
`data_type`, `*_index_ref`, the feature count (`n_genes`, `n_mirnas`, `n_proteins` or
`n_probes`) and the value arrays. Layer-specific fields: `platform` (expression,
miRNA, protein, methylation), `normalization` (expression), `analysis_method` (CNV),
`value_range` and `description` (methylation). Missing measurements are `null`.

### Index documents

Common fields: `_key`, `cohort`, `data_type`, feature count, `description`, plus
`platform`, `genome_version` or `analysis_method` where relevant. Entries of the
mapping arrays:

| Mapping | Entry fields |
|---|---|
| `gene_mappings` (expression) | `position`, `gene_id_ensembl`, `gene_id_base`, `gene_ref`, `entrez_id`, `hgnc_symbol` |
| `gene_mappings` (CNV) | `position`, `gene_id_ensembl`, `gene_id_base`, `gene_ref`, `entrez_id`, `enst_ids` |
| `mirna_mappings` | `position`, `mirna_id`, `mirbase_id`, `hgnc_symbol`, `description` |
| `protein_mappings` | `position`, `peptide_target`, `entrez_id`, `gene_symbol`, `protein_type` |
| `probe_mappings` | `position`, `probe_id`, `chromosome`, `genomic_start`, `genomic_end`, `strand`, `gene_symbols`, `gene_ids`, `gene_refs` |

Example of a `gene_mappings` entry:

```json
{"position": 0, "gene_id_ensembl": "ENSG00000000003.15", "gene_id_base": "ENSG00000000003",
 "gene_ref": "genes/ENSG00000000003", "entrez_id": "7105", "hgnc_symbol": "TSPAN6"}
```

### Optional PanCanAtlas view

The same collections can also host the harmonised TCGA PanCanAtlas 2018 release as the
cohort `PANCAN-ATLAS` (index keys `expr_index_PANCAN-ATLAS`, `cnv_index_PANCAN-ATLAS`,
`mirna_index_PANCAN-ATLAS`; expression values in `values_ebpp`). It is not part of the
published instance; see [scripts/pancan_atlas/README.md](scripts/pancan_atlas/README.md).

## Joining the layers

| Omic layer | Path to the knowledge graph |
|---|---|
| expression, CNV | `gene_mappings[].entrez_id` → `nodes` with `class_code == "EntrezID"` |
| protein | `protein_mappings[].entrez_id` → gene node; the protein node is reached through the `has gene product` edge |
| miRNA | `mirna_mappings[].hgnc_symbol` → `GENES` → `entrez_id` → gene node |
| methylation | `probe_mappings[].gene_refs` → `GENES` → `entrez_id` → gene node |

**Example: TPM of one gene across the samples of a cohort**

```aql
LET idx = DOCUMENT("GENE_EXPRESSION/expr_index_TCGA-BRCA")
LET pos = FIRST(FOR m IN idx.gene_mappings FILTER m.hgnc_symbol == "TP53" RETURN m.position)
FOR s IN GENE_EXPRESSION
  FILTER s.data_type == "gene_expression_vector" AND s.cohort == "TCGA-BRCA"
  RETURN { sample_id: s.sample_id, tpm: s.values_tpm[pos] }
```

**Example: pathways of the same gene in the knowledge graph**

```aql
FOR g IN nodes
  FILTER g.class_code == "EntrezID" AND g.entity_id == "7157"
  FOR v, e IN 1..1 OUTBOUND g edges
    FILTER v.class_code == "R-HSA" AND e.predicate_label == "participates in"
    RETURN { pathway: v._key, label: v.label }
```

The semantic collections (`nodes`, `edges`) are a stable, reusable layer; the omic
collections carry the per-sample quantitative evidence in compact vectors; the index
documents are the bridge between array positions and the entities of the graph.

More queries, and adapters that export task-specific datasets, are in
[queries/README.md](queries/README.md).
