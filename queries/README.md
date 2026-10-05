# Query collection

AQL queries over the integrated database, and small adapters that turn their
results into the input files of downstream analyses. Everything here only
**reads** the database.

```
queries/
├── common.py, run_query.py     connection helpers; run any .aql with bind variables
├── examples/                   the data model at work (four short queries)
├── extraction/                 building blocks for task-specific datasets
├── params/                     routing tables and layer layouts used by the adapters
└── adapters/                   query results -> files for downstream frameworks
```

Connection defaults (host, credentials, database `PKT_main`) come from
`scripts/arangodb_utils.py`; every script accepts `--db`, `--host`, `--user`,
`--password`.

## Examples

| Query | What it shows | Main bind variables |
|---|---|---|
| `01_two_step_retrieval.aql` | semantic seed → cohort index → per-sample vectors | `seed`, `cohort`, `value_field` |
| `02_predicate_stratified_pairs.aql` | gene pairs that share a node through one predicate (e.g. a pathway) | `predicates`, `shared_cc`, `shared_side`, `min_size`, `max_size` |
| `03_cross_layer_profile.aql` | one gene on CNV, mRNA and protein, sample by sample | `entrez`, `cohort` |
| `04_phenotype_anchored_set.aql` | genes reached from a phenotype, with their layer coverage | `phenotype`, `cohort` |

```bash
python queries/run_query.py queries/examples/01_two_step_retrieval.aql \
    --set seed=nodes/HP_0003002 --set cohort=TCGA-BRCA --set value_field=values_tpm --limit 3
python queries/run_query.py queries/examples/02_predicate_stratified_pairs.aql \
    --set 'predicates=["has participant"]' --set shared_cc=R-HSA --set shared_side=source \
    --set min_size=2 --set max_size=50 --out pairs.jsonl
```

Traversals are predicate-constrained on purpose: an unconstrained multi-hop
expansion reaches the class hubs of the ontology and grows without bound.

## Extraction queries

| Query | Returns |
|---|---|
| `typed_relation_subgraph.aql` | edges kept and renamed by (source type, predicate, target type); symmetric relations once per unordered pair |
| `gene_scale_edges.aql` | gene → higher-scale edges (pathway, disease, ...) for a gene panel, admitted by (target class, predicate) |
| `gene_ontology_annotation.aql` | gene → GO terms through the protein bridge (both bridge directions; one gene per protein; minimum support) |
| `ontology_closure.aql` | is_a / part_of edges among a given set of terms |
| `layer_index.aql` | ordered feature identifiers of one omic layer of one cohort |
| `layer_vectors.aql` | per-sample vectors of one layer (whole or sliced) |
| `sample_labels.aql` | a sample-level label stored on `SAMPLES` |
| `labelled_cohort.aql` | samples with expression, CNV and miRNA (and a label) |
| `aligned_feature_axis.aql` | features shared by two layers, with their positions |

Each file starts with a comment block listing its bind variables.

## Adapters

| Adapter | Output | Queries used |
|---|---|---|
| `to_typed_triples.py` | `head / interaction / tail / source / type` TSV, entities as `Type::id` | `typed_relation_subgraph` |
| `to_multiscale_edge_lists.py` | `node_*.csv` and `edge_*.csv` per node type and relation | `gene_scale_edges`, `gene_ontology_annotation`, `ontology_closure` |
| `to_omic_matrices.py` | feature × sample matrices in the Xena file layout (`.tsv.gz`), plus a label table | `layer_index`, `layer_vectors`, `sample_labels` |

```bash
python queries/adapters/to_typed_triples.py --out graph.tsv
python queries/adapters/to_multiscale_edge_lists.py --genes panel.csv --go-min-support 3 --out-dir hetero/
python queries/adapters/to_omic_matrices.py --layout queries/params/omic_layers_gdc.json \
    --cohort TCGA-BRCA --out-dir brca/
python queries/adapters/to_omic_matrices.py --layout queries/params/omic_layers_pancan_atlas.json \
    --out-dir pancan/          # needs the PanCanAtlas view (scripts/pancan_atlas/README.md)
```

What the adapters do **not** do: interaction layers that are not part of the
backbone (protein interaction networks, regulatory networks, curated
drug--target tables) are external inputs; `to_typed_triples.py --extra` appends
such a table, keyed on the same identifiers. Feature selection, normalisation and
cohort filtering belong to the downstream pipeline.

Exports are deterministic: rows are written in a fixed order, so two runs on the
same database produce identical files.
