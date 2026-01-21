#%%
"""
"""

from sqlalchemy import null
from arangodb_utils import *
import pandas as pd
# --- FUNZIONI DI BASE DI ARANGODB ---
db_connection = setup_arangodb_connection("PKT_test10000")

#%%  Load PKT NODELABELS as ref ===============================
pkt_nodelabels_path = "../data/pkt/builds/v3.0.2/PKT_NodeLabels_with_metadata_v3.0.2.csv"
pkt_nodelabels = pd.read_csv(pkt_nodelabels_path, dtype=str)
pkt_nodelabels['bioentity_type'].value_counts()
#%%
# serch in "db_connection" for collection "nodes"in nodes containing "tumor or cancer ion labels"
nodes = get_nodes_by_pattern(db_connection, "nodes", "label", "%tumor%")
nodes += get_nodes_by_pattern(db_connection, "nodes", "label", "%cancer%")
print(f"Found {len(nodes)} nodes with 'tumor' or 'cancer' in labels")
# example of nodes properties
for node in nodes[:5]:
    for key, value in node.items():
        print(f"{key}: {value}")
    print("-----")
    
#%%
nodes = get_nodes_by_pattern(db_connection, "nodes", "class_code", "%HP%")
print(f"Found {len(nodes)} nodes with class_code containing 'HP'")
# example of nodes properties
for node in nodes[:5]:
    for key, value in node.items():
        print(f"{key}: {value}")
    print("-----")

#%%
# serch "vital" in HP labels
nodes_query = get_nodes_by_pattern(db_connection, "nodes", "label", "% female%")
for n in nodes_query:
    print(f"{n['class_code']}: {n['label']}")
    
#%%

"""
in teoria sono mappabili su PKT
bioentity_type --> node property:
- gene --> entrez_id
- protein --> uniprot_id
- variant --> rsid
- disease --> mondo_id
- phenotype --> hp_id

KG key maps
gene expression data:
entrez_id  --> node with bioentity_type = gene

protein expression data:
uniprot_id --> node with bioentity_type = protein

variant data:
rsid --> node with bioentity_type = variant

disease associations:
mondo_id --> node with bioentity_type = disease

phenotype associations:
hp_id --> node with bioentity_type = phenotype


ArangoDB structure:
- database: PKT_test10000

Knowledge Graph in PKT_test10000:
- collection: nodes
  properties: _key, uri, namespace, entity_id, class_code, label, bioentity_type, description, synonym, source, source_type, integer_id, entrez_id, uniprot_id, rsid, mondo_id, hp_id

- collection: edges
  properties: _from, _to, edge_id, source_uri, target_uri, predicate_uri, predicate_label, predicate_class_code, predicate_bioentity_type, predicate_source

Omics Data Collections in PKT_test10000: 

| Collection              | Scope / contenuto principale                                | Chiave e campi chiave                                         | Uso tipico in query                                                                   |
| ----------------------- | ----------------------------------------------------------- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| PROJECTS                | Metadata TCGA di progetto (per es. TCGA-BRCA, TCGA-LUAD, …) | _key= project_id TCGA                                         | Navigare coorti, descrizioni e mapping project ↔ studi.                               |
| SAMPLES                 | Metadata campioni TCGA (barcode, case, project, type, ecc.) | _key= sample_id TCGA                                          | Filtrare subset di campioni per tipo, coorte, case e collegare ai layer quantitativi. |
| GENES                   | Metadata di geni, condivisi tra tutti gli studi             | _key= Ensembl base; es.ENSG00000141510                        | Join semantico dai layer quantitativi a info gene (symbol, biotype, ecc.).            |
| GENE_EXPRESSION_SAMPLES | Vettori di espressione per campione (TPM / counts / FPKM)   | _key=sample_id; campi:values_tpm,expr_index_ref,cohort        | Costruire matrici gene × sample via join conEXPRESSION_INDEXeGENES.                   |
| GENE_EXPRESSION_INDEX   | Indice globale dei geni per vettori di espressione          | _key=expr_index_TCGA-<STUDY>                                  | Recuperare l’ordine dei geni per interpretare i vettori diGENE_EXPRESSION_SAMPLES.    |
| CNV_SAMPLES             | Vettori CNV per campione                                    | _key=sample_id; campi:values_cnv,cnv_index_ref,cohort         | Analisi CNV per coorte / case, export verso ML.                                       |
| CNV_INDEX               | Indice feature per CNV (gene / segmenti)                    | _key=cnv_index_TCGA-<STUDY>                                   | Ottenere posizione feature CNV nei vettori diCNV_SAMPLES.                             |
| METHYLATION_SAMPLES     | Vettori di metilazione per campione                         | _key=sample_id; campi:values_methylation,methylation_index_ref,cohort | Analisi metilazione per coorte / case, export verso ML.                               |
| METHYLATION_INDEX       | Indice feature per metilazione (CpG sites
| MIRNA_SAMPLES           | Vettori di espressione miRNA per campione                   | _key=sample_id; campi:values_mirna,mirna_index_ref,cohort     | Analisi miRNA e integrazione multi-omica.                                             |
| MIRNA_INDEX             | Indice feature per miRNA                                    | _key=mirna_index_TCGA-<STUDY>                                 | Interpretare i vettori miRNA inMIRNA_SAMPLES.                                         |
| PROTEIN_SAMPLES         | Vettori proteomici per campione                             | _key=sample_id; campi:values_protein,protein_index_ref,cohort | Layer proteomico per modelli multi-omici.                                             |
| PROTEIN_INDEX           | Indice feature per proteomica (peptidi/proteine)            | _key=protein_index_TCGA-<STUDY>                               | Interpretare i vettori proteici inPROTEIN_SAMPLES.                                    |

in IDEXES the value "position" refers to the index in the omic data vector per sample ("_SAMPLES"), example "values_tpm" that corresponds to the gene.)

"""
def get_data_structure_info(collection_name):
    """
    get dict atrigutes kyes and values  from on entry in collection_name
    """
    aql_query = f"""
    FOR doc IN {collection_name}
      LIMIT 1
      RETURN doc
    """
    result = db_connection.aql.execute(aql_query)
    for doc in result:
        return doc

""" PROJECTS document example:
{'_key': 'TCGA-BRCA',
 '_id': 'PROJECTS/TCGA-BRCA',
 '_rev': '_kujDl82---',
 'name': 'Breast Invasive Carcinoma',
 'program': 'TCGA',
 'primary_site': 'Breast',
 'disease_types': ['Squamous Cell Neoplasms',
  'Adnexal and Skin Appendage Neoplasms',
  'Epithelial Neoplasms, NOS',
  'Complex Epithelial Neoplasms',
  'Fibroepithelial Neoplasms',
  'Cystic, Mucinous and Serous Neoplasms',
  'Basal Cell Neoplasms',
  'Adenomas and Adenocarcinomas',
  'Ductal and Lobular Neoplasms'],
 'entity_type': 'project',
 'n_cases': 1098,
 'n_samples': 1255}
"""
# get_data_structure_info("PROJECTS")

""" SAMPLES document example:
{'_key': 'TCGA-BH-A0W3-01A',
 '_id': 'SAMPLES/TCGA-BH-A0W3-01A',
 '_rev': '_kujMAnC---',
 'submitter_id': 'TCGA-BH-A0W3',
 'case_ref': 'cases/TCGA-BH-A0W3',
 'project_ref': 'projects/TCGA-BRCA',
 'sample_type': 'Primary Tumor',
 'sample_type_id': 1,
 'tissue_type': 'Tumor',
 'tumor_descriptor': 'Primary',
 'specimen_type': 'Solid Tissue',
 'composition': 'Not Reported',
 'preservation_method': 'OCT',
 'days_to_collection': 85,
 'entity_type': 'sample'}
"""
# get_data_structure_info("SAMPLES")

""" GENES  document example:
{'_key': 'ENSG00000290825',
 '_id': 'GENES/ENSG00000290825',
 '_rev': '_kujDS_q---',
 'entity_type': 'gene',
 'bioentity_type': 'gene',
 'gene_stable_id': 'ENSG00000290825',
 'gene_stable_id_version': 'ENSG00000290825.2',
 'hgnc_symbol': None,
 'entrez_id': '727856',
 'uniprot_id': None,
 'mirbase_id': None,
 'gene_type': 'lncRNA',
 'gene_description': 'DEAD/H-box helicase 11 like 16 (pseudogene) [Source:NCBI gene (formerly Entrezgene);Acc:727856]',
 'chromosome': '1',
 'gene_start_bp': '11121',
 'gene_end_bp': '24894',
 'strand': '1',
 'transcript_ids': ['ENST00000832824',
  'ENST00000832825',
  'ENST00000832826',
  'ENST00000832827',
  'ENST00000832828',
  'ENST00000832829',
  'ENST00000832830',
  'ENST00000832837',
...
  'ENST00000832848',
  'ENST00000832849',
  'ENST00000832823'],
 'source': 'Ensembl_BioMart',
 'source_version': 'GRCh38.v36'}
"""
# get_data_structure_info("GENES")

""" GENE_EXPRESSION_SAMPLES document example:
{'_key': 'TCGA-D8-A146-01A',
 '_id': 'GENE_EXPRESSION_SAMPLES/TCGA-D8-A146-01A',
 '_rev': '_kuh2f5O---',
 'sample_id': 'TCGA-D8-A146-01A',
 'cohort': 'TCGA-BRCA',
 'data_type': 'gene_expression_vector',
 'expression_index_ref': 'expression_index/expr_index_TCGA-BRCA',
 'platform': 'STAR',
 'n_genes': 60660,
 'values_counts': None,
 'values_fpkm': None,
 'values_tpm': [5.662037403353558,
  3.37609586203262,
  6.8601396907537975,
  4.400551591403899,
  2.84516885873,
  ...,
  ...,
  ...
  ]
}
"""
# get_data_structure_info("GENE_EXPRESSION_SAMPLES")

"""  GENE_EXPRESSION_INDEX document example:

  {'position': 167,
   'gene_id_ensembl': 'ENSG00000007080.11',
   'gene_id_base': 'ENSG00000007080',
   'gene_ref': 'genes/ENSG00000007080',
   'entrez_id': '115098',
   'hgnc_symbol': 'CCDC124'},
  {'position': 168,
   'gene_id_ensembl': 'ENSG00000007129.18',
   'gene_id_base': 'ENSG00000007129',
   'gene_ref': 'genes/ENSG00000007129',
   'entrez_id': '90273',
   'hgnc_symbol': 'CEACAM21'},
  {'position': 169,
   'gene_id_ensembl': 'ENSG00000007168.14',
   'gene_id_base': 'ENSG00000007168',
   'gene_ref': 'genes/ENSG00000007168',
   'entrez_id': '5048',
   'hgnc_symbol': 'PAFAH1B1'},
  {...}
"""
# get_data_structure_info("GENE_EXPRESSION_INDEX")

""" CNV_SAMPLES document example:
{'_key': 'TCGA-D8-A1XU-01A',
 '_id': 'CNV_SAMPLES/TCGA-D8-A1XU-01A',
 '_rev': '_k7ty2ii---',
 'sample_id': 'TCGA-D8-A1XU-01A',
 'cohort': 'TCGA-BRCA',
 'data_type': 'cnv_vector',
 'cnv_index_ref': 'cnv_index/cnv_index_TCGA-BRCA',
 'analysis_method': 'ASCAT3',
 'n_genes': 60623,
 'values_copy_number': [None,
  None,
  None,
  None,
  2,
  2,
  2,
  2,
...
  2,
  2,
  2,
  2,
  ...]}
"""
# get_data_structure_info("CNV_SAMPLES")

""" CNV_INDEX document example:
{'_key': 'cnv_index_TCGA-BRCA',
 '_id': 'CNV_INDEX/cnv_index_TCGA-BRCA',
 '_rev': '_k7uFnie---',
 'cohort': 'TCGA-BRCA',
 'data_type': 'cnv_index',
 'n_genes': 60623,
 'analysis_method': 'ASCAT3',
 'gene_mappings': [{'position': 0,
   'gene_id_ensembl': 'ENSG00000223972.5',
   'gene_id_base': 'ENSG00000223972',
   'gene_ref': 'genes/ENSG00000223972',
   'entrez_id': '100287102',
   'enst_ids': ['ENST00000450305']},
  {'position': 1,
   'gene_id_ensembl': 'ENSG00000227232.5',
   'gene_id_base': 'ENSG00000227232',
   'gene_ref': 'genes/ENSG00000227232',
   'entrez_id': '653635',
   'enst_ids': ['ENST00000488147']},
  {'position': 2,
   'gene_id_ensembl': 'ENSG00000278267.1',
   'gene_id_base': 'ENSG00000278267',
   'gene_ref': 'genes/ENSG00000278267',
   'entrez_id': '102466751',
   'enst_ids': ['ENST00000619216']},
...
   'gene_ref': 'genes/ENSG00000230523',
   'entrez_id': '101928460',
   'enst_ids': ['ENST00000417651']},
  ...],
 'description': 'Gene position mapping for CNV vectors in sample documents'}
"""
# get_data_structure_info("CNV_INDEX")

"""METHYLATION_SAMPLES document example:
{'_key': 'TCGA-A7-A0CD-01A',
 '_id': 'METHYLATION_SAMPLES/TCGA-A7-A0CD-01A',
 '_rev': '_kujMAE----',
 'sample_id': 'TCGA-A7-A0CD-01A',
 'cohort': 'TCGA-BRCA',
 'data_type': 'methylation_vector',
 'methylation_index_ref': 'methylation_index/methylation_index_TCGA-BRCA',
 'platform': 'Illumina HumanMethylation27',
 'n_probes': 27578,
 'values_beta': [0.891275578180925,
  0.0337566988200663,
  0.789521222963392,
  0.619850753240135,
  0.0619325703831766,
  0.0218987623260666,
  0.985018975396962,
  0.0132840494835158,
  0.0123023732109782,
  None,
  0.90177723238615,
  0.0161257521439586,
  0.190583229352712,
  0.0201717551329457,
  0.10779725647655,
  0.0204887059820779,
...
  0.320087517917168,
  0.0741297816754439,
  ...],
 'value_range': '[0.0, 1.0]',
 'description': 'Beta values representing methylation levels at CpG sites'}
"""
# get_data_structure_info("METHYLATION_SAMPLES")

""" METHYLATION_INDEX document example:
{'_key': 'methylation_index_TCGA-BRCA',
 '_id': 'METHYLATION_INDEX/methylation_index_TCGA-BRCA',
 '_rev': '_k7rvuzu---',
 'cohort': 'TCGA-BRCA',
 'data_type': 'methylation_index',
 'n_probes': 27578,
 'platform': 'Illumina HumanMethylation27',
 'genome_version': 'GRCh38',
 'probe_mappings': [{'position': 0,
   'probe_id': 'cg00000292',
   'chromosome': 'chr16',
   'genomic_start': 28878778,
   'genomic_end': 28878780,
   'strand': '-',
   'gene_symbols': ['AC009093.11', 'ATP2A1'],
   'gene_ids': ['ENSG00000196296'],
   'gene_refs': ['genes/ENSG00000196296']},
  {'position': 1,
   'probe_id': 'cg00002426',
   'chromosome': 'chr3',
   'genomic_start': 57757815,
   'genomic_end': 57757817,
   'strand': '-',
   'gene_symbols': ['SLMAP'],
   'gene_ids': ['ENSG00000163681'],
...
   'gene_symbols': ['SYTL4'],
   'gene_ids': ['ENSG00000102362'],
   'gene_refs': ['genes/ENSG00000102362']},
  ...],
 'description': 'CpG probe position mapping for methylation beta value vectors in sample documents'}
"""
# get_data_structure_info("METHYLATION_INDEX")

""" MIRNA_SAMPLES document example:
{'_key': 'TCGA-D8-A146-01A',
 '_id': 'MIRNA_SAMPLES/TCGA-D8-A146-01A',
 '_rev': '_kuh3bmC---',
 'sample_id': 'TCGA-D8-A146-01A',
 'cohort': 'TCGA-BRCA',
 'data_type': 'mirna_vector',
 'mirna_index_ref': 'mirna_index/mirna_index_TCGA-BRCA',
 'platform': 'Illumina',
 'n_mirnas': 1881,
 'values_expression': [12.815167145771248,
  12.816050718114727,
  12.855015959831322,
  14.730158401240104,
  11.169628343075129,
  8.859832599955457,
  9.505781518958475,
...
  0.395469768620962,
  0.395469768620962,
  0,
  0,
  ...]}
"""
# get_data_structure_info("MIRNA_SAMPLES")

""" MIRNA_INDEX document example:
{'_key': 'mirna_index_TCGA-BRCA',
 '_id': 'MIRNA_INDEX/mirna_index_TCGA-BRCA',
 '_rev': '_k7uFnqq---',
 'cohort': 'TCGA-BRCA',
 'data_type': 'mirna_index',
 'n_mirnas': 1881,
 'platform': 'Illumina',
 'mirna_mappings': [{'position': 0,
   'mirna_id': 'hsa-let-7a-1',
   'mirbase_id': 'hsa-let-7a-1',
   'hgnc_symbol': 'MIRLET7A-1',
   'description': 'miRNA hsa-let-7a-1'},
  {'position': 1,
   'mirna_id': 'hsa-let-7a-2',
   'mirbase_id': 'hsa-let-7a-2',
   'hgnc_symbol': 'MIRLET7A-2',
   'description': 'miRNA hsa-let-7a-2'},
  {'position': 2,
   'mirna_id': 'hsa-let-7a-3',
   'mirbase_id': 'hsa-let-7a-3',
   'hgnc_symbol': 'MIRLET7A-3',
   'description': 'miRNA hsa-let-7a-3'},
  {'position': 3,
   'mirna_id': 'hsa-let-7b',
   'mirbase_id': 'hsa-let-7b',
...
   'mirbase_id': 'hsa-mir-4717',
   'hgnc_symbol': 'MIR4717',
   'description': 'miRNA hsa-mir-4717'},
  ...],
 'description': 'miRNA position mapping for expression vectors in sample documents'}
"""
# get_data_structure_info("MIRNA_INDEX")

""" PROTEIN_SAMPLES document example:
{'_key': 'TCGA-WT-AB41-01A',
 '_id': 'PROTEIN_SAMPLES/TCGA-WT-AB41-01A',
 '_rev': '_k7tzEcK---',
 'sample_id': 'TCGA-WT-AB41-01A',
 'cohort': 'TCGA-BRCA',
 'data_type': 'protein_vector',
 'protein_index_ref': 'protein_index/protein_index_TCGA-BRCA',
 'platform': 'RPPA',
 'n_proteins': 487,
 'values_abundance': [0.06515,
  0.051282,
  0.23858,
  0.62398,
  -0.39313,
  0.7282,
  -0.15681,
...
  -0.32683,
  0.1320667,
  -0.1346183,
  1.442465,
  0]}
"""
# get_data_structure_info("PROTEIN_SAMPLES")

""" PROTEIN_INDEX document example:
{'_key': 'protein_index_TCGA-BRCA',
 '_id': 'PROTEIN_INDEX/protein_index_TCGA-BRCA',
 '_rev': '_k7uFnru---',
 'cohort': 'TCGA-BRCA',
 'data_type': 'protein_index',
 'n_proteins': 487,
 'platform': 'RPPA',
 'protein_mappings': [{'position': 0,
   'peptide_target': '1433BETA',
   'entrez_id': '7529',
   'gene_symbol': 'YWHAB',
   'protein_type': 'Total'},
  {'position': 1,
   'peptide_target': '1433EPSILON',
   'entrez_id': '7531',
   'gene_symbol': 'YWHAE',
   'protein_type': 'Total'},
  {'position': 2,
   'peptide_target': '1433ZETA',
   'entrez_id': '7534',
   'gene_symbol': 'YWHAZ',
   'protein_type': 'Total'},
  {'position': 3,
   'peptide_target': '4EBP1',
   'entrez_id': '1977',
...
   'peptide_target': 'ZEB1',
   'entrez_id': '9708',
   'gene_symbol': 'ZEB1',
   'protein_type': 'Total'}],
 'description': 'Protein/peptide position mapping for RPPA vectors in sample documents'}
"""
# get_data_structure_info("PROTEIN_INDEX")

#%%
"""
### Strategie Query AQL Raccomandate

Per estrarre sottografi multi-omici integrati:

Data structure in PKT_test10000:
# 4. Collezione GENE_EXPRESSION_SAMPLES (normalizzato)
{"_key": "expr_TCGA-D8-A1XU-01A_ENSG00000223972.5",
  "sample_id": "TCGA-D8-A1XU-01A",
  "gene_id": "ENSG00000223972.5",
  "star_counts": 11.737669863177452,
  "star_fpkm": 5.123456,
  "star_tpm": 6.543210,
  "platform": "STAR",
  "analysis_date": "2015-03-15"
}


```aql
// Query trans-omica per pathway-specific expression
FOR gene IN genes
  FILTER gene.entrez_id IN @pathway_genes
  FOR expr IN gene_expression_TCGA_BRCA
    FILTER expr.gene_ref == CONCAT("genes/", gene._key)
    FILTER expr.sample_id IN @tumor_samples
    RETURN {
      gene: gene.hgnc_symbol,
      sample: expr.sample_id,
      expression: expr.tpm
    }
```
"""

def query_gene_expression_by_pathway(db_connection, 
                                     pathway_genes, 
                                     tumor_samples):
    aql_query = """
    FOR gene IN genes
      FILTER gene.entrez_id IN @pathway_genes
      FOR expr IN gene_expression_TCGA_BRCA
        FILTER expr.gene_ref == CONCAT("genes/", gene._key)
        FILTER expr.sample_id IN @tumor_samples
        RETURN {
          gene: gene.hgnc_symbol,
          sample: expr.sample_id,
          expression: expr.tpm
        }
    """
    bind_vars = {
        "pathway_genes": pathway_genes,
        "tumor_samples": tumor_samples
    }
    result = db_connection.aql.execute(aql_query, bind_vars=bind_vars)
    return list(result)

# example usage
pathway_genes = ["7157", "7422", "1956"]  # Example Entrez IDs for TP53, EGFR, AKT1
tumor_samples = ["TCGA-A1-A0SB-01A", "TCGA-BH-A0B3-01A"]  # Example sample IDs
expression_data = query_gene_expression_by_pathway(db_connection, pathway_genes, tumor_samples)
for record in expression_data:
    print(record)
    

"""

## Architettura v4: Sample-Centric Vector Strategy

La nuova versione risolve completamente il problema di scalabilità attraverso una ristrutturazione radicale del layer quantitativo, mantenendo inalterato il layer semantico.[3][4]

### Struttura Multi-Layer

**1. Semantic Layer (invariato)**
- `genes`: ~20K nodi con annotazioni ontologiche complete
- `projects`: metadati studio TCGA

**2. Quantitative Index Layer (nuovo)**
- `expression_index`: un documento per coorte contenente l'array ordinato `gene_mappings` che mappa posizione→gene_id
- `cnv_index`, `mirna_index`, `protein_index`: analoghi per altri layer omici

**3. Quantitative Vector Layer (completamente ridisegnato)**
- `gene_expression_samples_{COHORT}`: **~1K documenti** invece di ~20M
- Ogni documento = un sample con vettori paralleli `values_tpm`, `values_fpkm`, `values_counts`
- Riferimento esplicito a `expression_index_ref` per interpretazione

### Riduzione Scalabilità

La strategia v4 ottiene una **riduzione di 20.000× nel numero di documenti** rispetto alla v3:
- v3: 20.000 geni × 1.000 samples = 20.000.000 documenti
- v4: 1.000 samples = 1.000 documenti + 1 indice

Storage stimato: ~100MB vs ~10GB per coorte.[5][3]

### Query Pattern Ottimizzate

**Accesso sample-specific** (caso d'uso primario per ML/GNN):
```aql
FOR s IN gene_expression_samples_TCGA_BRCA
  FILTER s.sample_id == "TCGA-A1-A0SB-01A"
  RETURN s.values_tpm  // O(1) lookup
```

**Estrazione gene-specific cross-cohort**:
```aql
LET idx = DOCUMENT("expression_index/expr_index_TCGA-BRCA")
LET pos = FIRST(FOR g IN idx.gene_mappings 
                FILTER g.hgnc_symbol == "TP53" 
                RETURN g.position)
FOR s IN gene_expression_samples_TCGA_BRCA
  RETURN {sample: s.sample_id, tp53: s.values_tpm[pos]}
```

**Costruzione matrici features per ML**:
```aql
FOR s IN gene_expression_samples_TCGA_BRCA
  RETURN {sample_id: s.sample_id, features: s.values_tpm}
```

Questa architettura è compatibile con le funzionalità di vector indexing native di ArangoDB 3.12+ per query di similarità, utili per identificare campioni con profili trascrizionali simili.[4][3][5]

### Integrazione con il Framework IRR

La strategia v4 mantiene la separazione concettuale del paradigma IRR:[2]
- **Representation**: il layer semantico (`genes`, ontologie) rimane application-agnostic
- **Reuse**: i vettori quantitativi sono riutilizzabili per molteplici analisi (GNN, network analysis, survival modeling)
- **Quantitative Evidence**: separata topologicamente dalla conoscenza stabile, permettendo aggiornamenti incrementali per nuove coorti senza duplicare annotazioni ontologiche

"""