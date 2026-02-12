## ArangoDB PKT_test10000 Database Structure

### Database Configuration

- **Database name**: PKT_test10000
- **Graph model**: Property Graph

***

### Knowledge Graph Collections

#### Collection: **nodes**

**Purpose**: Semantic layer containing ontologically grounded biological entities

**Properties**:

- `_key` (string): Unique node identifier (e.g., "PR_Q9H609", "79177")
- `uri` (string): Canonical URI of the biological entity
- `namespace` (string): Source namespace (e.g., "purl.obolibrary.org", "www.ncbi.nlm.nih.gov")
- `entity_id` (string): Original entity identifier from source
- `class_code` (string): Ontology or database class code (e.g., "PR", "EntrezID", "MONDO", "HP")
- `label` (string): Human-readable label
- `bioentity_type` (string): Entity type classifier — **gene**, **protein**, **variant**, **disease**, **phenotype**, **go**, etc.
- `description` (string): Textual description of the entity
- `synonym` (string): Alternative names (pipe-separated)
- `source` (string): Source database or ontology name
- `source_type` (string): Type of source ("Ontology", "Database")
- `integer_id` (integer): Internal numeric identifier

**Mapping keys for omics integration**:

- `entrez_id` (string): Links to gene entities (bioentity_type = "gene")
- `uniprot_id` (string): Links to protein entities (bioentity_type = "protein")
- `rsid` (string): Links to variant entities (bioentity_type = "variant")
- `mondo_id` (string): Links to disease entities (bioentity_type = "disease")
- `hp_id` (string): Links to phenotype entities (bioentity_type = "phenotype")

***

#### Collection: **edges**

**Purpose**: Semantic relationships between biological entities

**Properties**:

- `_from` (string): Source node reference (collection/key format)
- `_to` (string): Target node reference (collection/key format)
- `edge_id` (string): Unique edge identifier
- `source_uri` (string): URI of source entity
- `target_uri` (string): URI of target entity
- `predicate_uri` (string): Relationship type URI (from ontologies)
- `predicate_label` (string): Human-readable relationship label (e.g., "gene product of", "associated with")
- `predicate_class_code` (string): Ontology code for predicate (e.g., "RO")
- `predicate_bioentity_type` (string): Semantic type of the relationship
- `predicate_source` (string): Source ontology of the predicate

***

### Omics Data Collections

All omics collections follow a **unified vector + index architecture** where each collection contains two document types:

1. **Vector documents** (`*_vector`): Per-s
ample arrays of quantitative values
2. **Index documents** (`*_index`): Position-to-feature mappings for interpreting vectors

| Collection      | Scope / Content                                          | Key Fields                                                                                                              | Typical Use                                                                      |
| --------------- | -------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| PROJECTS        | TCGA project metadata (e.g., TCGA-BRCA, TCGA-LUAD)       | `_key`= project_id                                                                                                      | Navigate cohorts, descriptions, project ↔ study mappings                         |
| SAMPLES         | TCGA sample metadata (barcode, case, project, type)      | `_key`= sample_id                                                                                                       | Filter sample subsets by type, cohort, case; link to quantitative layers         |
| GENES           | Gene metadata shared across all studies                  | `_key`= Ensembl base (e.g., ENSG00000141510)                                                                            | Semantic join from quantitative layers to gene info (symbol, biotype, etc.)      |
| GENE_EXPRESSION | Sample expression vectors + index                        | Vector: `_key`=sample_id, `values_tpm[]`, `expression_index_ref`<br>Index: `_key`=expr_index_TCGA-\<STUDY\>, `gene_mappings[]` | Build gene × sample matrices; retrieve gene order for interpreting vectors       |
| CNV             | Sample CNV vectors + index                               | Vector: `_key`=sample_id, `values_copy_number[]`, `cnv_index_ref`<br>Index: `_key`=cnv_index_TCGA-\<STUDY\>, `gene_mappings[]` | CNV analysis by cohort/case; get gene→position mapping in vector                 |
| METHYLATION     | Sample methylation vectors + index                       | Vector: `_key`=sample_id, `values_beta[]`, `methylation_index_ref`<br>Index: `_key`=methylation_index_TCGA-\<STUDY\>, `probe_mappings[]` | Methylation analysis by cohort/case; probe→position mapping                      |
| MIRNA           | Sample miRNA expression vectors + index                  | Vector: `_key`=sample_id, `values_expression[]`, `mirna_index_ref`<br>Index: `_key`=mirna_index_TCGA-\<STUDY\>, `mirna_mappings[]` | miRNA analysis and multi-omic integration; miRNA→position mapping                |
| PROTEIN         | Sample proteomic vectors + index                         | Vector: `_key`=sample_id, `values_abundance[]`, `protein_index_ref`<br>Index: `_key`=protein_index_TCGA-\<STUDY\>, `protein_mappings[]` | Proteomic layer for multi-omic models; peptide→position mapping                  |

***

### Supporting Collections

#### Collection: **PROJECTS**

**Purpose**: TCGA project-level metadata

**Properties**:

- `_key` (string): TCGA project identifier (e.g., "TCGA-BRCA")
- `name` (string): Project full name (e.g., "Breast Invasive Carcinoma")
- `program` (string): Program name (e.g., "TCGA")
- `primary_site` (string): Primary anatomical site (e.g., "Breast")
- `disease_types` (array): List of disease type classifications
- `entity_type` (string): Entity type ("project")
- `n_cases` (integer): Number of cases in project
- `n_samples` (integer): Number of samples in project

***

#### Collection: **SAMPLES**

**Purpose**: Biological specimen and tissue sample metadata

**Properties**:

- `_key` (string): Sample barcode identifier (e.g., "TCGA-BH-A0W3-01A")
- `submitter_id` (string): Case submitter ID (e.g., "TCGA-BH-A0W3")
- `case_ref` (string): Reference to case (e.g., "cases/TCGA-BH-A0W3")
- `project_ref` (string): Reference to project (e.g., "projects/TCGA-BRCA")
- `sample_type` (string): Type of sample (e.g., "Primary Tumor")
- `sample_type_id` (integer): Numeric sample type code
- `tissue_type` (string): Tissue classification (e.g., "Tumor")
- `tumor_descriptor` (string): Tumor descriptor (e.g., "Primary")
- `specimen_type` (string): Specimen processing type (e.g., "Solid Tissue")
- `composition` (string): Sample composition
- `preservation_method` (string): Preservation protocol (e.g., "OCT")
- `days_to_collection` (integer): Days from diagnosis to sample collection
- `entity_type` (string): Entity type ("sample")

***

#### Collection: **GENES**

**Purpose**: Central gene metadata repository shared across all studies

**Properties**:

- `_key` (string): Ensembl gene ID base (e.g., "ENSG00000290825")
- `entity_type` (string): Entity type ("gene")
- `bioentity_type` (string): Bio-entity type ("gene")
- `gene_stable_id` (string): Ensembl stable ID without version
- `gene_stable_id_version` (string): Ensembl stable ID with version (e.g., "ENSG00000290825.2")
- `hgnc_symbol` (string): HGNC gene symbol (e.g., "TP53")
- `entrez_id` (string): Entrez Gene ID
- `uniprot_id` (string): UniProt protein ID
- `mirbase_id` (string): miRBase ID (for miRNA genes)
- `gene_type` (string): Gene biotype (e.g., "protein_coding", "lncRNA")
- `gene_description` (string): Gene description from source
- `chromosome` (string): Chromosome location
- `gene_start_bp` (string): Gene start position (bp)
- `gene_end_bp` (string): Gene end position (bp)
- `strand` (string): DNA strand ("1" or "-1")
- `transcript_ids` (array): List of associated transcript IDs
- `source` (string): Data source (e.g., "Ensembl_BioMart")
- `source_version` (string): Source version (e.g., "GRCh38.v36")

**KG mapping**: `entrez_id` → nodes(`_key`) where bioentity_type = "gene"

***

### Omics Vector Collections (Detailed)

#### Collection: **GENE_EXPRESSION**

**Purpose**: Quantitative gene expression measurements from RNA-seq

**Document type: `gene_expression_vector`**

- `_key` (string): Sample identifier (e.g., "TCGA-D8-A146-01A")
- `sample_id` (string): Sample identifier
- `cohort` (string): TCGA cohort (e.g., "TCGA-BRCA")
- `data_type` (string): "gene_expression_vector"
- `expression_index_ref` (string): Reference to index document (e.g., "expr_index_TCGA-BRCA")
- `platform` (string): Analysis platform (e.g., "STAR")
- `n_genes` (integer): Number of genes in vector
- `values_counts` (array): Raw read counts (optional)
- `values_fpkm` (array): FPKM normalized values (optional)
- `values_tpm` (array): TPM normalized values

**Document type: `gene_expression_index`**

- `_key` (string): Index identifier (e.g., "expr_index_TCGA-BRCA")
- `cohort` (string): TCGA cohort
- `data_type` (string): "gene_expression_index"
- `n_genes` (integer): Number of genes
- `platform` (string): Analysis platform
- `genome_version` (string): Genome assembly version (e.g., "GRCh38.v36")
- `gene_mappings` (array): Position-to-gene mappings
  - `position` (integer): Index in value arrays
  - `gene_id_ensembl` (string): Ensembl ID with version
  - `gene_id_base` (string): Ensembl ID base
  - `gene_ref` (string): Reference to GENES collection
  - `entrez_id` (string): Entrez Gene ID
  - `hgnc_symbol` (string): HGNC symbol

***

#### Collection: **CNV**

**Purpose**: Copy Number Variation data per sample

**Document type: `cnv_vector`**

- `_key` (string): Sample identifier
- `sample_id` (string): Sample identifier
- `cohort` (string): TCGA cohort
- `data_type` (string): "cnv_vector"
- `cnv_index_ref` (string): Reference to index document
- `analysis_method` (string): CNV calling method (e.g., "ASCAT3")
- `n_genes` (integer): Number of genes
- `values_copy_number` (array): Copy number values per gene (nullable)

**Document type: `cnv_index`**

- `_key` (string): Index identifier (e.g., "cnv_index_TCGA-BRCA")
- `cohort` (string): TCGA cohort
- `data_type` (string): "cnv_index"
- `n_genes` (integer): Number of genes
- `analysis_method` (string): CNV calling method
- `gene_mappings` (array): Position-to-gene mappings
  - `position` (integer): Index in value arrays
  - `gene_id_ensembl` (string): Ensembl ID with version
  - `gene_id_base` (string): Ensembl ID base
  - `gene_ref` (string): Reference to GENES collection
  - `entrez_id` (string): Entrez Gene ID
  - `enst_ids` (array): Associated transcript IDs

***

#### Collection: **METHYLATION**

**Purpose**: DNA methylation beta values at CpG sites

**Document type: `methylation_vector`**

- `_key` (string): Sample identifier
- `sample_id` (string): Sample identifier
- `cohort` (string): TCGA cohort
- `data_type` (string): "methylation_vector"
- `methylation_index_ref` (string): Reference to index document
- `platform` (string): Methylation array platform (e.g., "Illumina HumanMethylation27")
- `n_probes` (integer): Number of CpG probes
- `values_beta` (array): Beta values [0.0, 1.0] representing methylation levels
- `value_range` (string): Value range description
- `description` (string): Data description

**Document type: `methylation_index`**

- `_key` (string): Index identifier (e.g., "methylation_index_TCGA-BRCA")
- `cohort` (string): TCGA cohort
- `data_type` (string): "methylation_index"
- `n_probes` (integer): Number of probes
- `platform` (string): Array platform
- `genome_version` (string): Genome assembly version
- `probe_mappings` (array): Position-to-probe mappings
  - `position` (integer): Index in value arrays
  - `probe_id` (string): CpG probe identifier (e.g., "cg00000292")
  - `chromosome` (string): Chromosome location
  - `genomic_start` (integer): Genomic start position
  - `genomic_end` (integer): Genomic end position
  - `strand` (string): DNA strand
  - `gene_symbols` (array): Associated gene symbols
  - `gene_ids` (array): Associated Ensembl gene IDs
  - `gene_refs` (array): References to GENES collection

***

#### Collection: **MIRNA**

**Purpose**: miRNA expression measurements

**Document type: `mirna_vector`**

- `_key` (string): Sample identifier
- `sample_id` (string): Sample identifier
- `cohort` (string): TCGA cohort
- `data_type` (string): "mirna_vector"
- `mirna_index_ref` (string): Reference to index document
- `platform` (string): Sequencing platform (e.g., "Illumina")
- `n_mirnas` (integer): Number of miRNAs
- `values_expression` (array): Expression values (log2 RPM)

**Document type: `mirna_index`**

- `_key` (string): Index identifier (e.g., "mirna_index_TCGA-BRCA")
- `cohort` (string): TCGA cohort
- `data_type` (string): "mirna_index"
- `n_mirnas` (integer): Number of miRNAs
- `platform` (string): Sequencing platform
- `mirna_mappings` (array): Position-to-miRNA mappings
  - `position` (integer): Index in value arrays
  - `mirna_id` (string): miRNA identifier (e.g., "hsa-let-7a-1")
  - `mirbase_id` (string): miRBase accession
  - `hgnc_symbol` (string): HGNC symbol (e.g., "MIRLET7A-1")
  - `description` (string): miRNA description

***

#### Collection: **PROTEIN**

**Purpose**: Reverse Phase Protein Array (RPPA) abundance measurements

**Document type: `protein_vector`**

- `_key` (string): Sample identifier
- `sample_id` (string): Sample identifier
- `cohort` (string): TCGA cohort
- `data_type` (string): "protein_vector"
- `protein_index_ref` (string): Reference to index document
- `platform` (string): Proteomic platform ("RPPA")
- `n_proteins` (integer): Number of proteins/antibodies
- `values_abundance` (array): Normalized protein abundance values

**Document type: `protein_index`**

- `_key` (string): Index identifier (e.g., "protein_index_TCGA-BRCA")
- `cohort` (string): TCGA cohort
- `data_type` (string): "protein_index"
- `n_proteins` (integer): Number of proteins
- `platform` (string): Proteomic platform
- `protein_mappings` (array): Position-to-protein mappings
  - `position` (integer): Index in value arrays
  - `peptide_target` (string): Antibody target name (e.g., "1433BETA")
  - `entrez_id` (string): Entrez Gene ID
  - `gene_symbol` (string): HGNC gene symbol
  - `protein_type` (string): Measurement type (e.g., "Total", "Phosphorylated")

***

### Integration Architecture

The **unified vector + index architecture** enables efficient storage and retrieval of multi-omics data:

1. **Vector documents** store dense arrays of values per sample, minimizing document count
2. **Index documents** provide the semantic mapping from array positions to biological features
3. **References** (`*_ref` fields) link samples to their corresponding indices

**Query Example: Retrieve TPM of a specific gene for all samples**

```aql
// Retrieve TPM of TP53 for all BRCA samples
LET gene_symbol = "TP53"
LET index_doc = DOCUMENT("GENE_EXPRESSION/expr_index_TCGA-BRCA")
LET gene_pos = FIRST(
  FOR m IN index_doc.gene_mappings
    FILTER m.hgnc_symbol == gene_symbol
    RETURN m.position
)
FOR sample IN GENE_EXPRESSION
  FILTER sample.data_type == "gene_expression_vector"
  FILTER sample.cohort == "TCGA-BRCA"
  RETURN {
    sample_id: sample._key,
    tpm: sample.values_tpm[gene_pos]
  }
```

**KG Integration Mapping Strategy**:

1. **Gene expression/CNV** → Match via `entrez_id` in gene_mappings → nodes where bioentity_type = "gene"
2. **Protein expression** → Match via `entrez_id`/`gene_symbol` in protein_mappings → nodes
3. **miRNA expression** → Match via `mirbase_id` in mirna_mappings → nodes
4. **Methylation** → Match via `gene_refs` in probe_mappings → GENES → nodes

This decoupled architecture follows the **Information Representation and Reuse (IRR) paradigm**, where:

- The `nodes` and `edges` collections form a **stable, reusable semantic layer**
- Omics collections contain **dynamic quantitative evidence** in efficient vector format
- Index documents provide **semantic bridges** between quantitative arrays and biological annotations
- Ontological grounding ensures **semantic interoperability** across heterogeneous data sources

***
