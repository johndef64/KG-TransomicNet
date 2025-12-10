

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

"""
in teoria sono mappabili su PKT
bioentity_type --> node property:
- gene --> entrez_id
- protein --> uniprot_id
- variant --> rsid
- disease --> mondo_id
- phenotype --> hp_id
"""

# serch in db_connection for collection "nodes"in nodes containing "tumor or cancer ion labels"
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

Knowledge Graph:
- collection: nodes
  properties: _key, uri, namespace, entity_id, class_code, label, bioentity_type, description, synonym, source, source_type, integer_id, entrez_id, uniprot_id, rsid, mondo_id, hp_id

- collection: edges
  properties: _from, _to, edge_id, source_uri, target_uri, predicate_uri, predicate_label, predicate_class_code, predicate_bioentity_type, predicate_source

Omics Data Collections:

  | Collection              | Scope / contenuto principale                                | Chiave e campi chiave                                         | Uso tipico in query                                                                   |
| ----------------------- | ----------------------------------------------------------- | ------------------------------------------------------------- | ------------------------------------------------------------------------------------- |
| GENES                   | Metadata di geni, condivisi tra tutti gli studi             | _key= Ensembl base; es.ENSG00000141510                        | Join semantico dai layer quantitativi a info gene (symbol, biotype, ecc.).            |
| EXPRESSION_INDEX        | Indice globale dei geni per vettori di espressione          | _key=expr_index_TCGA-<STUDY>                                  | Recuperare l’ordine dei geni per interpretare i vettori diGENE_EXPRESSION_SAMPLES.    |
| CNV_INDEX               | Indice feature per CNV (gene / segmenti)                    | _key=cnv_index_TCGA-<STUDY>                                   | Ottenere posizione feature CNV nei vettori diCNV_SAMPLES.                             |
| MIRNA_INDEX             | Indice feature per miRNA                                    | _key=mirna_index_TCGA-<STUDY>                                 | Interpretare i vettori miRNA inMIRNA_SAMPLES.                                         |
| PROTEIN_INDEX           | Indice feature per proteomica (peptidi/proteine)            | _key=protein_index_TCGA-<STUDY>                               | Interpretare i vettori proteici inPROTEIN_SAMPLES.                                    |
| PROJECTS                | Metadata TCGA di progetto (per es. TCGA-BRCA, TCGA-LUAD, …) | _key= project_id TCGA                                         | Navigare coorti, descrizioni e mapping project ↔ studi.                               |
| SAMPLES                 | Metadata campioni TCGA (barcode, case, project, type, ecc.) | _key= sample_id TCGA                                          | Filtrare subset di campioni per tipo, coorte, case e collegare ai layer quantitativi. |
| GENE_EXPRESSION_SAMPLES | Vettori di espressione per campione (TPM / counts / FPKM)   | _key=sample_id; campi:values_tpm,expr_index_ref,cohort        | Costruire matrici gene × sample via join conEXPRESSION_INDEXeGENES.                   |
| CNV_SAMPLES             | Vettori CNV per campione                                    | _key=sample_id; campi:values_cnv,cnv_index_ref,cohort         | Analisi CNV per coorte / case, export verso ML.                                       |
| MIRNA_SAMPLES           | Vettori di espressione miRNA per campione                   | _key=sample_id; campi:values_mirna,mirna_index_ref,cohort     | Analisi miRNA e integrazione multi-omica.                                             |
| PROTEIN_SAMPLES         | Vettori proteomici per campione                             | _key=sample_id; campi:values_protein,protein_index_ref,cohort | Layer proteomico per modelli multi-omici.                                             |

"""

"""
### Strategie Query AQL Raccomandate

Per estrarre sottografi multi-omici integrati:

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