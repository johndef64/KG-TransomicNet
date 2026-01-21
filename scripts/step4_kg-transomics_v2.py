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

def test_query(query_aql):
    result = db_connection.aql.execute(query_aql)
    for doc in result:
        return doc
    
# for data structure see "db_structure_omics_data.md"


#%%
"""
## Pattern di mapping KG↔omiche
Nel KG (`nodes`) il mapping suggerito è: gene→`entrezid`, protein→`uniprotid`, variant→`rsid`, disease→`mondoid`, phenotype→`hpid`, e questo è il ponte principale per collegare evidenza quantitativa ai nodi semantici. 
Nei layer omici TCGA il mapping si appoggia agli “index documents” per coorte: `GENEEXPRESSIONINDEX.genemappings` porta a `GENES` via `generef`, CNV analogo (`CNVINDEX.genemappings`), miRNA (`MIRNAINDEX.mirnamappings`), proteomica (`PROTEININDEX.proteinmappings`), metilazione (`METHYLATIONINDEX.probemappings` con `generefs`). 

Pattern AQL consigliati (riusabili in tutte le query successive):
- **Risoluzione di un gene KG → posizione indice (coorte)** (evita join “per riga” sui vettori):  
```aql
LET idx = DOCUMENT(@@exprIndexCol, @exprIndexKey)   // es. @@exprIndexCol="GENEEXPRESSIONINDEX"
LET pos = FIRST(
  FOR m IN idx.genemappings
    FILTER m.entrezid == @entrezId OR m.hgncsymbol == @hgnc
    RETURN m.position
)
RETURN pos
```
"""

#- **Risoluzione di un gene KG → posizione indice (coorte)** (evita join “per riga” sui vettori): 
AQL = """
LET idx = DOCUMENT(@@exprIndexCol, @@exprIndexKey)   // es. @@exprIndexCol="GENEEXPRESSIONINDEX"
LET pos = FIRST(
  FOR m IN idx.genemappings
    FILTER m.entrezid == @entrezId OR m.hgncsymbol == @hgnc
    RETURN m.position
)
RETURN pos
"""

AQL_BASE = """
FOR gene IN GENES
  RETURN {
    gene_id: gene.entrez_id,
    hgnc_symbol: gene.hgnc_symbol
  }
"""

def test_query(query_aql):
    result = db_connection.aql.execute(query_aql)
    for doc in result:
        return doc

test_query(AQL_BASE)

#data la struttyra di questi dati, devi sviuluppare un siostema AQL python che permetta il recupero dei valori quantitativi dagli array dei samples per la variabile inserita esempio entrex_id

resolve_gene_to_position(db_connection, "TCGA-BRCA", entrez_id="727856")  # TP53
# test_query(AQL_BASE)
#%%


# ...existing code...

from omics_query_utils import (
    resolve_gene_position,
    get_gene_expression_by_gene,
    get_gene_expression_matrix,
    get_cnv_by_gene,
    get_methylation_by_gene,
    get_mirna_expression,
    get_protein_abundance,
    get_multiomics_by_gene,
    get_sample_full_profile
)

#%% Test: Resolve gene position
pos = resolve_gene_position(db_connection, "TCGA-BRCA", entrez_id="7157")
print(f"TP53 position in expression vector: {pos}")

#%% Test: Get TP53 expression across all BRCA samples
tp53_expr = get_gene_expression_by_gene(
    db_connection, 
    cohort="TCGA-BRCA",
    hgnc_symbol="TP53",
    value_type="tpm"
)
print(f"Found expression for {len(tp53_expr)} samples")
for r in tp53_expr[:5]:
    print(f"  {r['sample_id']}: {r['value']:.2f} TPM")

#%% Test: Get expression matrix for multiple genes
pathway_genes = ["7157", "7422", "1956"]  # TP53, VEGFA, EGFR
matrix_data = get_gene_expression_matrix(
    db_connection,
    cohort="TCGA-BRCA",
    gene_ids=pathway_genes,
    id_type="entrez_id",
    value_type="tpm"
)
print(f"Matrix: {len(matrix_data['samples'])} samples x {len(matrix_data['genes'])} genes")

# show the matriox as dataframe
df_matrix = pd.DataFrame(matrix_data,
                        index=matrix_data['samples'],
                        columns=matrix_data['genes'])
df_matrix


#%% Test: Get multi-omics data for TP53
tp53_multiomics = get_multiomics_by_gene(
    db_connection,
    cohort="TCGA-BRCA",
    hgnc_symbol="TP53",
    entrez_id="7157"
)
print(f"TP53 multi-omics summary:")
print(f"  Expression samples: {len(tp53_multiomics['expression'])}")
print(f"  CNV samples: {len(tp53_multiomics['cnv'])}")
print(f"  Methylation probes: {len(tp53_multiomics['methylation'].get('probes', []))}")
print(f"  Protein samples: {len(tp53_multiomics['protein'])}")

#%% Test: Get full profile for a single sample
sample_profile = get_sample_full_profile(
    db_connection, 
    sample_id="TCGA-D8-A146-01A", 
    cohort="TCGA-BRCA"
)
print(f"Sample profile keys: {sample_profile.keys()}")














#%%


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

#%%

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