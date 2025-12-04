## ArangoDB PKT_test10000 Database Structure

### Database Configuration

- **Database name**: KG_TransomicNet_test10000
- **Graph model**: Property Graph

***

### Core Collections

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

#### Collection: **GENE_EXPRESSION**

**Purpose**: Quantitative gene expression measurements from RNA-seq experiments

**Data sources**: `star_counts`, `star_fpkm`, `star_tpm`

**Properties**:

- `_key` (string): Unique identifier (format: `expr_{sample_id}_{gene_id}`)
- `sample_id` (string): Biological sample identifier (e.g., "TCGA-D8-A1XU-01A")
- `gene_id` (string): Ensembl gene identifier with version (e.g., "ENSG00000223972.5")
- `entrez_id` (string): **Mapping key to nodes** — links to node with bioentity_type = "gene"
- `star_counts` (float): Raw read counts from STAR aligner
- `star_fpkm` (float): Fragments Per Kilobase Million normalized expression
- `star_tpm` (float): Transcripts Per Million normalized expression
- `platform` (string): Sequencing/analysis platform (e.g., "STAR")

**KG mapping**: `entrez_id` → nodes(`_key`) where bioentity_type = "gene"

***

#### Collection: **PROTEIN_EXPRESSION**

**Purpose**: Quantitative protein abundance measurements from proteomic assays

**Data source**: `protein` (RPPA)

**Properties**:

- `_key` (string): Unique identifier (format: `prot_{sample_id}_{rppa_id}`)
- `sample_id` (string): Biological sample identifier
- `rppa_id` (string): Reverse Phase Protein Array antibody identifier
- `uniprot_id` (string): **Mapping key to nodes** — links to node with bioentity_type = "protein"
- `entrez_id` (string): Secondary mapping to corresponding gene node
- `abundance` (float): Normalized protein abundance value
- `platform` (string): Proteomic platform (e.g., "RPPA")

**KG mapping**: `uniprot_id` → nodes(`_key`) where bioentity_type = "protein"

***

#### Collection: **VARIANTS**

**Purpose**: Somatic mutation and genetic variant data

**Data source**: `somaticmutation_wxs` (Whole Exome Sequencing)

**Properties**:

- `_key` (string): Unique identifier (format: `mut_{sample_id}_{gene_symbol}_{position}`)
- `sample_id` (string): Biological sample identifier
- `gene_id` (string): Ensembl gene identifier
- `entrez_id` (string): Mapping to gene node
- `rsid` (string): **Mapping key to nodes** — links to node with bioentity_type = "variant"
- `gene_symbol` (string): HGNC gene symbol (e.g., "TP53")
- `chromosome` (string): Chromosome location (e.g., "chr17")
- `start` (integer): Genomic start position
- `ref` (string): Reference allele
- `alt` (string): Alternate allele
- `amino_acid_change` (string): Protein-level mutation annotation (e.g., "R175H")
- `effect` (string): Variant effect classification (e.g., "missense_variant")
- `dna_vaf` (float): Variant Allele Frequency (DNA level)
- `callers` (array): List of variant calling algorithms used

**KG mapping**: `rsid` → nodes(`_key`) where bioentity_type = "variant"

***

#### Collection: **DISEASE_ASSOCIATIONS**

**Purpose**: Clinical disease annotations and diagnoses linked to samples

**Data source**: `clinical` (TCGA diagnoses)

**Properties**:

- `_key` (string): Unique identifier (format: `disease_{sample_id}_{mondo_id}`)
- `sample_id` (string): Biological sample identifier
- `case_id` (string): Patient/case identifier
- `mondo_id` (string): **Mapping key to nodes** — links to node with bioentity_type = "disease"
- `primary_diagnosis` (string): Clinical diagnosis text
- `disease_type` (string): Disease category
- `tissue_or_organ_of_origin` (string): Anatomical site
- `morphology` (string): ICD-O morphology code
- `icd_10_code` (string): ICD-10 classification code
- `tumor_grade` (string): Histological grade
- `ajcc_pathologic_stage` (string): AJCC staging (e.g., "Stage IIA")
- `ajcc_pathologic_t` (string): Tumor size category
- `ajcc_pathologic_n` (string): Lymph node involvement
- `ajcc_pathologic_m` (string): Metastasis status

**KG mapping**: `mondo_id` → nodes(`_key`) where bioentity_type = "disease"

***

#### Collection: **PHENOTYPE_ASSOCIATIONS**

**Purpose**: Phenotypic observations and clinical features linked to samples

**Data source**: `clinical` (derived from demographics, exposures, treatments)

**Properties**:

- `_key` (string): Unique identifier (format: `phenotype_{sample_id}_{hp_id}`)
- `sample_id` (string): Biological sample identifier
- `case_id` (string): Patient/case identifier
- `hp_id` (string): **Mapping key to nodes** — links to node with bioentity_type = "phenotype"
- `phenotype_label` (string): Human Phenotype Ontology term label
- `vital_status` (string): Patient vital status (e.g., "Alive", "Deceased")
- `age_at_diagnosis_years` (float): Patient age at diagnosis
- `gender` (string): Patient gender
- `race` (string): Patient race
- `ethnicity` (string): Patient ethnicity
- `treatment_type` (string): Type of treatment received
- `treatment_or_therapy` (string): Whether treatment was administered
- `overall_survival_time` (float): Survival time in years
- `overall_survival_status` (integer): Survival status (0 = alive, 1 = deceased)

**KG mapping**: `hp_id` → nodes(`_key`) where bioentity_type = "phenotype"

***

### Supporting Collections

#### Collection: **SAMPLES**

**Purpose**: Biological specimens and tissue samples

**Properties**:

- `_key` (string): Sample barcode identifier (e.g., "TCGA-D8-A1XU-01A")
- `case_id` (string): Reference to patient/case
- `sample_id` (string): UUID identifier
- `sample_type` (string): Type of sample (e.g., "Primary Tumor")
- `sample_type_id` (integer): Numeric sample type code
- `tissue_type` (string): Tissue classification
- `specimen_type` (string): Specimen processing type
- `preservation_method` (string): Sample preservation protocol


#### Collection: **CASES**

**Purpose**: Patient/case-level clinical and demographic data

**Properties**:

- `_key` (string): Case submitter ID (e.g., "TCGA-BH-A0W3")
- `id` (string): UUID identifier
- `tcga_project` (string): TCGA project code
- `primary_site` (string): Primary anatomical site
- `demographics` (object): Nested demographic information
- `diagnoses` (array): Array of diagnosis records
- `treatments` (array): Array of treatment records
- `survival` (object): Survival outcome data

***

### Integration Architecture

The **node_id mapping strategy** enables seamless integration between omics data collections and the semantic knowledge graph:

1. **Gene expression/methylation/miRNA** → Match on `entrez_id`
2. **Protein expression** → Match on `uniprot_id`
3. **Variants/mutations** → Match on `rsid`
4. **Disease associations** → Match on `mondo_id`
5. **Phenotype associations** → Match on `hp_id`

This decoupled architecture follows the **Information Representation and Reuse (IRR) paradigm**, where:

- The `nodes` and `edges` collections form a **stable, reusable semantic layer**
- Omics collections contain **dynamic quantitative evidence** linked via controlled identifiers
- Ontological grounding ensures **semantic interoperability** across heterogeneous data sources[^1]

***
