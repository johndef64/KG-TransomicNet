# TCGA Omics Data Loader for ArangoDB

## Overview
This script (`step3.1_load_to_arangodb.py`) loads TCGA omics datasets into ArangoDB collections, handling both study-independent and study-specific data appropriately.

## Data Structure

### Study-Independent Collections (Incremental)
These collections accumulate data across all studies and skip duplicates:

1. **genes** - Gene metadata mapping
   - Contains: gene IDs (Ensembl, Entrez), HGNC symbols
   - Skip rule: If `_key` exists, skip insertion
   - Source: `genes.json`

2. **Index Collections** - Position mappings for omics vectors
   - `expression_index` - Gene expression position mapping
   - `cnv_index` - Copy number variation position mapping
   - `mirna_index` - miRNA position mapping
   - `protein_index` - Protein expression position mapping
   - Skip rule: If `_key` exists, skip insertion
   - Source: `*_index.json`

### Study-Specific Collections
These collections contain data specific to each TCGA study:

3. **Metadata Collections**
   - `projects` - TCGA project/study metadata
   - `cases` - Patient case information
   - `samples` - Sample metadata (tumor/normal tissue info)

4. **Omics Data Collections**
   - `gene_expression_samples` - RNA-seq expression vectors
   - `cnv_samples` - Copy number variation vectors
   - `mirna_samples` - miRNA expression vectors
   - `protein_samples` - Protein expression vectors
   - Format: Each document contains sample_id and data vector

## Key Features

### 1. Incremental Loading
- Study-independent data (genes, indexes) is loaded only once
- Duplicate entries are automatically skipped based on `_key`
- Multiple studies can be loaded sequentially without data duplication

### 2. Batch Processing
- Data is imported in batches of 1000 documents
- Progress is reported every 10 batches
- Failed batches are logged but don't stop the process

### 3. Flexible File Format Support
- Handles both JSON arrays and JSONL (JSON Lines) format
- Automatically detects format based on file content
- JSONL is used for large sample data files

### 4. Error Handling
- Missing files are reported but don't stop the import
- Failed batches are counted and reported
- Comprehensive error messages for debugging

## Usage

### Basic Usage
```python
# Set the study to load
STUDY = 'BRCA'  # Breast Invasive Carcinoma

# Set database name
db_name = 'PKT_test10000'

# Run the script
python step3.1_load_to_arangodb.py
```

### Loading Multiple Studies
To load multiple TCGA studies sequentially:

```python
studies = ['BRCA', 'LUAD', 'COAD', 'PRAD']

for study in studies:
    STUDY = study
    data_dir = f"../data/arangodb_collections/TCGA-{STUDY}/"
    # Run import process
    import_data_to_arangodb(db)
```

## Configuration

### Database Settings
```python
db_name = 'PKT_test10000'           # Database name
arangodb_hosts = 'http://localhost:8529'
arangodb_user = 'root'
arangodb_password = 'avocadodb'
```

### Study Settings
```python
STUDY = 'BRCA'  # Current study code
data_dir = f"../data/arangodb_collections/TCGA-{STUDY}/"
```

### Batch Size
```python
BATCH_SIZE = 1000  # Documents per batch
```

## Indexes Created

### Genes Collection
- `gene_id_base` (hash index)
- `hgnc_symbol` (hash index)
- `entrez_id` (hash index)

### Index Collections
- `position` (hash index)
- `gene_id_base` (hash index)

### Samples Collection
- `submitter_id` (hash index)
- `case_ref` (hash index)
- `project_ref` (hash index)

### Omics Data Collections
- `sample_id` (hash index)

## Data Files Expected

For each TCGA study (e.g., TCGA-BRCA), the following files should be present:

```
data/arangodb_collections/TCGA-BRCA/
├── genes.json                              # Study-independent
├── expression_index.json                    # Study-independent
├── cnv_index.json                          # Study-independent
├── mirna_index.json                        # Study-independent
├── protein_index.json                      # Study-independent
├── projects.json                           # Study-specific
├── cases.json                              # Study-specific
├── samples.json                            # Study-specific
├── gene_expression_samples_TCGA-BRCA.json  # Study-specific
├── cnv_samples_TCGA-BRCA.json              # Study-specific
├── mirna_samples_TCGA-BRCA.json            # Study-specific
└── protein_samples_TCGA-BRCA.json          # Study-specific
```

## Output

The script provides:
1. Progress updates during loading
2. Count of documents imported/skipped/failed per collection
3. Final statistics showing total documents in each collection

Example output:
```
============================================================
IMPORTING TCGA-BRCA DATA
============================================================

📂 Loading genes from genes.json...
  Loaded 60660 documents
  Processed 60660/60660 documents...
✓ Imported 60660 documents to genes (skipped: 0, failed: 0)

📂 Loading expression_index from expression_index.json...
  Loaded 60660 documents
  Processed 60660/60660 documents...
✓ Imported 60660 documents to expression_index (skipped: 0, failed: 0)

...

============================================================
✅ Import completed!
   Total imported: 185432
   Total skipped: 0
============================================================

============================================================
COLLECTION STATISTICS
============================================================
  genes: 60,660 documents
  expression_index: 60,660 documents
  cnv_index: 60,660 documents
  mirna_index: 60,660 documents
  protein_index: 60,660 documents
  projects: 1 documents
  cases: 1,098 documents
  samples: 1,255 documents
  gene_expression_samples: 1,098 documents
  cnv_samples: 1,098 documents
  mirna_samples: 1,098 documents
  protein_samples: 1,098 documents
============================================================
```

## Query Examples

### Get all samples for a project
```python
samples = list(db.aql.execute(
    'FOR s IN samples FILTER s.project_ref == @project RETURN s',
    bind_vars={'project': 'projects/TCGA-BRCA'}
))
```

### Get gene expression for a specific sample
```python
expr_data = list(db.aql.execute(
    'FOR e IN gene_expression_samples FILTER e.sample_id == @sample_id RETURN e',
    bind_vars={'sample_id': 'TCGA-BH-A0W3-01A'}
))
```

### Get genes by symbol
```python
genes = list(db.aql.execute(
    'FOR g IN genes FILTER g.hgnc_symbol == @symbol RETURN g',
    bind_vars={'symbol': 'TP53'}
))
```

### Join sample with expression data
```python
results = list(db.aql.execute('''
    FOR s IN samples
        FILTER s.project_ref == @project
        FOR e IN gene_expression_samples
            FILTER e.sample_id == s._key
            RETURN MERGE(s, {expression: e})
''', bind_vars={'project': 'projects/TCGA-BRCA'}))
```

## Troubleshooting

### Connection Issues
- Ensure ArangoDB is running on localhost:8529
- Check username/password credentials
- Verify database permissions

### Missing Files
- Check that all required JSON files exist in data directory
- Verify file names match the expected pattern
- Check file permissions

### Memory Issues
- Reduce BATCH_SIZE if memory errors occur
- Process files in smaller chunks
- Consider increasing ArangoDB memory limits

### Duplicate Data
- Study-independent collections automatically skip duplicates
- Use `skip_duplicates=True` for custom collections
- Check `_key` values for uniqueness

## Notes

1. **Large Files**: Gene metadata (genes.json) can be >50MB. The script handles streaming where possible.

2. **JSONL Format**: Sample data files use JSONL (one JSON per line) for efficiency with large datasets.

3. **Incremental Design**: The script is designed to be run multiple times for different studies without data duplication.

4. **Index Performance**: Indexes are created after data import for better performance.

5. **Error Recovery**: Failed batches don't stop the import process, allowing partial data recovery.
