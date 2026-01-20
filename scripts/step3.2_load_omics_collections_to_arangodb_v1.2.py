#%%

from arango import ArangoClient
import json
import traceback
import os
import math

STUDY = 'BRCA' # Change this to load different TCGA studies

# --- CONFIGURAZIONE ---
db_name = 'PKT_test10000'
# db_name = 'PKT_transomics_v1'
arangodb_hosts = 'http://localhost:8529'
arangodb_user = 'root'
arangodb_password = 'avocadodb'

data_dir = f"../data/arangodb_collections/TCGA-{STUDY}/"

print(os.getcwd())

# show all file in data dir
print("Data files in directory:")
file_list = []
for f in os.listdir(data_dir):
    print(" -", f)
    file_list.append(f)

# STUDY = 'XXXX' (whatever you set it to)
sorted_files = [
    # genes
    'genes.json',
    # expression
    'gene_expression_index.json',
    f'gene_expression_samples_TCGA-{STUDY}.json',
    # copy-number
    'cnv_index.json',
    f'cnv_samples_TCGA-{STUDY}.json',
    # miRNA
    'mirna_index.json',
    f'mirna_samples_TCGA-{STUDY}.json',
    # protein
    'protein_index.json',
    f'protein_samples_TCGA-{STUDY}.json',
    # methylation (NEW!)
    'methylation_index.json',
    f'methylation_samples_TCGA-{STUDY}.json',
    # projects
    'projects.json',
    # samples
    'samples.json',
]

BATCH_SIZE = 1000 # Dimensione del batch per le operazioni di inserimento

# --- FUNZIONI DI BASE DI ARANGODB ---

def setup_arangodb_connection(db_name: str):
    """Setup ArangoDB connection and create database if needed."""
    try:
        client = ArangoClient(hosts=arangodb_hosts)
        # Connetti prima al database di sistema
        sys_db = client.db('_system', username=arangodb_user, password=arangodb_password)
        
        # Crea database se non esiste
        if not sys_db.has_database(db_name):
            sys_db.create_database(db_name)
            print(f"✓ Created database: {db_name}")
        else:
            print(f"Database {db_name} already exists")
        
        # Connetti al database di destinazione
        db = client.db(db_name, username=arangodb_user, password=arangodb_password)
        print(f"✓ Connected to database: {db_name}")
        return db
    
    except Exception as e:
        print(f"✗ Error setting up ArangoDB connection. Is ArangoDB running on {arangodb_hosts}? Error: {e}")
        return None

def create_tcga_collections(db):
    """
    Crea le collection per i dataset TCGA.
    Nessuna collection viene droppata: carichi incrementale.
    """
    if db is None:
        return
    
    try:
        collections_vertex = [
            "GENES",
            "GENE_EXPRESSION_INDEX",
            "CNV_INDEX",
            "MIRNA_INDEX",
            "PROTEIN_INDEX",
            "METHYLATION_INDEX", # NEW!
            "PROJECTS",
            "SAMPLES",
            "GENE_EXPRESSION_SAMPLES",
            "CNV_SAMPLES",
            "MIRNA_SAMPLES",
            "PROTEIN_SAMPLES",
            "METHYLATION_SAMPLES", # NEW!
        ]
        
        for cname in collections_vertex:
            if not db.has_collection(cname):
                db.create_collection(cname)
                print(f"✓ Created collection: {cname}")
            else:
                print(f"Collection {cname} already exists ({db.collection(cname).count()} documents)")
        
        # Indici essenziali (oltre all'indice implicito su _key)
        try:
            genes_col = db.collection("GENES")
            genes_col.add_index({
                "type": "hash",
                "fields": ["symbol"],
                "unique": False,
                "sparse": True
            })
            
            samples_col = db.collection("SAMPLES")
            samples_col.add_index({
                "type": "hash",
                "fields": ["sample_id"],
                "unique": True,
                "sparse": False
            })
            
            ge_col = db.collection("GENE_EXPRESSION_SAMPLES")
            ge_col.add_index({
                "type": "hash",
                "fields": ["sample_id"],
                "unique": True,
                "sparse": False
            })
            ge_col.add_index({
                "type": "hash",
                "fields": ["cohort"],
                "unique": False,
                "sparse": True
            })
            
            # Aggiungi indici simili per CNV/MIRNA/PROTEIN/METHYLATION
            for col_name in ["CNV_SAMPLES", "MIRNA_SAMPLES", "PROTEIN_SAMPLES", "METHYLATION_SAMPLES"]:
                col = db.collection(col_name)
                col.add_index({
                    "type": "hash",
                    "fields": ["sample_id"],
                    "unique": True,
                    "sparse": False
                })
                col.add_index({
                    "type": "hash",
                    "fields": ["cohort"],
                    "unique": False,
                    "sparse": True
                })
            
            print("✓ Created additional indexes")
        except Exception as e:
            print(f"Warning: Could not create some indexes: {e}")
    
    except Exception as e:
        print(f"✗ Error creating TCGA collections: {e}")
        raise

# --- UTILS PER CARICARE E INSERIRE JSON ---

def sanitize_value(val):
    """
    Recursively sanitize values in nested structures.
    Converts NaN, Infinity, -Infinity to None (null in JSON).
    """
    if isinstance(val, float):
        if math.isnan(val) or math.isinf(val):
            return None
        return val
    elif isinstance(val, list):
        return [sanitize_value(v) for v in val]
    elif isinstance(val, dict):
        return {k: sanitize_value(v) for k, v in val.items()}
    else:
        return val

def _load_json_lines(path: str):
    """
    Carica un file JSON che può essere:
    - array JSON unico
    - JSON Lines (un documento JSON per riga)
    
    Sanitizza automaticamente valori NaN/Inf.
    """
    with open(path, "r", encoding="utf-8") as f:
        first = f.read(1)
        f.seek(0)
        
        if not first:
            return []
        
        if first == "[":
            # Array JSON
            docs = json.load(f)
        else:
            # JSON Lines
            docs = []
            for line in f:
                line = line.strip()
                if not line:
                    continue
                docs.append(json.loads(line))
    
    # Sanitize all documents
    return [sanitize_value(doc) for doc in docs]

def _insert_incremental(collection, docs, unique_key: str = "_key"):
    """
    Inserimento incrementale.
    Assumiamo che ogni documento abbia _key unico e usiamo insert_many(overwrite=False).
    Questo evita di sovrascrivere documenti già presenti.
    """
    total = len(docs)
    if total == 0:
        print("  No documents to insert.")
        return
    
    inserted = 0
    failed = 0
    
    for i in range(0, total, BATCH_SIZE):
        batch = docs[i:i + BATCH_SIZE]
        try:
            collection.insert_many(batch, overwrite=False, silent=True)
            inserted += len(batch)
            if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= total:
                print(f"  Inserted {min(i + BATCH_SIZE, total)}/{total}")
        except Exception as e:
            failed += len(batch)
            print(f"  ⚠ Batch {i // BATCH_SIZE + 1} failed: {e}")
    
    print(f"  → Inserted {inserted}, failed {failed}")

# --- IMPORTAZIONE DATASETS TCGA ---

def import_tcga_datasets(db, 
                         PARTIAL_LOAD_LIST=[],
                         REPLACE_EXISTING=False):
    """
    Carica:
    - GENES (global, dedup su _key: insert senza overwrite)
    - *_INDEX (globali, uno per study/coorte, incrementali)
    - *_SAMPLES (quantitativo per-sample, study-specific ma sempre per _key)
    - PROJECTS, SAMPLES (meta, incrementali)
    """
    if db is None:
        print("✗ No database connection available")
        return
    
    if REPLACE_EXISTING:
        print("\n⚠ REPLACE_EXISTING is True: Dropping existing collections before load.")
        for col_name in sorted_files:
            collection_name = col_name.split('.')[0].upper()
            if db.has_collection(collection_name):
                db.delete_collection(collection_name)
                print(f"  Dropped collection: {collection_name}")
        create_tcga_collections(db)
    
    # 1. GENES (global, no overwrite; entries con stesso _key non vengono ricaricate)
    genes_path = os.path.join(data_dir, "genes.json")
    if os.path.exists(genes_path):
        print(f"\n📂 Loading GENES from {genes_path}...")
        genes_docs = _load_json_lines(genes_path)
        genes_col = db.collection("GENES")
        _insert_incremental(genes_col, genes_docs, unique_key="_key")
    else:
        print("GENES file not found, skipping.")
    
    # 2. Index collections (global, uno per study/coorte)
    collection_names_map = {
        "index_map" :{
            "gene_expression_index.json": "GENE_EXPRESSION_INDEX",
            "cnv_index.json": "CNV_INDEX",
            "mirna_index.json": "MIRNA_INDEX",
            "protein_index.json": "PROTEIN_INDEX",
            "methylation_index.json": "METHYLATION_INDEX" # NEW!
        },
        "sample_map" :{
            f"gene_expression_samples_TCGA-{STUDY}.json": "GENE_EXPRESSION_SAMPLES",
            f"cnv_samples_TCGA-{STUDY}.json": "CNV_SAMPLES",
            f"mirna_samples_TCGA-{STUDY}.json": "MIRNA_SAMPLES",
            f"protein_samples_TCGA-{STUDY}.json": "PROTEIN_SAMPLES",
            f"methylation_samples_TCGA-{STUDY}.json": "METHYLATION_SAMPLES" # NEW!
    }
    }

    index_map = collection_names_map["index_map"]
    sample_map = collection_names_map["sample_map"]

    if len(PARTIAL_LOAD_LIST) > 0:
        # keep only specified collections that contatins elements in PARTIAL_LOAD_LIST
        index_map = {k:v for k,v in index_map.items() if any(item in k for item in PARTIAL_LOAD_LIST)}
        sample_map = {k:v for k,v in sample_map.items() if any(item in k for item in PARTIAL_LOAD_LIST)}

    
    # Load Omics collections
    for fname, colname in index_map.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            print(f"{fname} not found, skipping.")
            continue
        print(f"\n📂 Loading {colname} from {path}...")
        docs = _load_json_lines(path)
        col = db.collection(colname)
        _insert_incremental(col, docs, unique_key="_key")
    
    for fname, colname in sample_map.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            print(f"{fname} not found, skipping.")
            continue
        print(f"\n📂 Loading {colname} (study {STUDY}) from {path}...")
        docs = _load_json_lines(path)
        col = db.collection(colname)
        _insert_incremental(col, docs, unique_key="_key")
    
    # 4. PROJECTS, SAMPLES (meta, global/incrementale)
    meta_map = {
        "projects.json": "PROJECTS",
        "samples.json": "SAMPLES",
    }
    
    for fname, colname in meta_map.items():
        path = os.path.join(data_dir, fname)
        if not os.path.exists(path):
            print(f"{fname} not found, skipping.")
            continue
        print(f"\n📂 Loading {colname} from {path}...")
        docs = _load_json_lines(path)
        col = db.collection(colname)
        _insert_incremental(col, docs, unique_key="_key")

#%%
db = setup_arangodb_connection(db_name)
create_tcga_collections(db)
#%%
if __name__ == "__main__":
    print("=" * 50)
    print("ARANGODB IMPORT SCRIPT - TCGA DATASETS (v4.2)")
    print("=" * 50)
    
    # 1. Connessione e Configurazione
    db = setup_arangodb_connection(db_name)
    
    if db:
        # 2. Creazione collection TCGA (no drop, incrementale)
        create_tcga_collections(db)
        
        # 3. Importazione Dati TCGA
        import_tcga_datasets(db, PARTIAL_LOAD_LIST=[
            "gene_expression",
            "cnv",
            "mirna", 
            "protein",
            "methylation"],
            REPLACE_EXISTING=True)  
        
        print("\n✓ TCGA datasets import completed.")
    else:
        print("✗ Failed to connect to ArangoDB. Please check if ArangoDB is running.")

# %%
# db = setup_arangodb_connection(db_name)

# check collections and number of documents
if db:
    # show db create_tcga_collections
    print("\n--- Collections in database ---")
    for col_name in db.collections():
        if col_name['name'].startswith('_'):
            continue
        collection = db.collection(col_name['name'])
        count = db.collection(col_name['name']).count()
        print(f"Collection: {col_name['name']}, Documents: {count}")


#%%