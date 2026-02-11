#%%

"""
TCGA Omics Data Loader for ArangoDB - Version 3

This version consolidates index and sample data into unified collections:
- GENE_EXPRESSION (contains both index and per-sample data)
- CNV (contains both index and per-sample data)
- MIRNA (contains both index and per-sample data)
- PROTEIN (contains both index and per-sample data)
- METHYLATION (contains both index and per-sample data)

The "data_type" field in each document indicates whether it's an index or per-sample data.
"""

from arango import ArangoClient
import json
import traceback
import os
import math

STUDY = 'BRCA'  # Change this to load different TCGA studies

# --- CONFIGURATION ---
db_name = 'PKT_test10000'
# db_name = 'PKT_transomics_v1'
arangodb_hosts = 'http://localhost:8529'
arangodb_user = 'root'
arangodb_password = 'avocadodb'

data_dir = f"../data/arangodb_collections/TCGA-{STUDY}/"

print(os.getcwd())

# Show all files in data dir
print("Data files in directory:")
file_list = []
for f in os.listdir(data_dir):
    print(" -", f)
    file_list.append(f)

# Mapping of JSON files to unified collections
# Each omic type has both index and samples loaded into the same collection
OMIC_TYPES = {
    "gene_expression": {
        "collection": "GENE_EXPRESSION",
        "index_file": "gene_expression_index.json",
        "samples_file": f"gene_expression_samples_TCGA-{STUDY}.json"
    },
    "cnv": {
        "collection": "CNV",
        "index_file": "cnv_index.json",
        "samples_file": f"cnv_samples_TCGA-{STUDY}.json"
    },
    "mirna": {
        "collection": "MIRNA",
        "index_file": "mirna_index.json",
        "samples_file": f"mirna_samples_TCGA-{STUDY}.json"
    },
    "protein": {
        "collection": "PROTEIN",
        "index_file": "protein_index.json",
        "samples_file": f"protein_samples_TCGA-{STUDY}.json"
    },
    "methylation": {
        "collection": "METHYLATION",
        "index_file": "methylation_index.json",
        "samples_file": f"methylation_samples_TCGA-{STUDY}.json"
    }
}

BATCH_SIZE = 1000  # Batch size for insert operations

# --- ARANGODB BASE FUNCTIONS ---

def setup_arangodb_connection(db_name: str):
    """Setup ArangoDB connection and create database if needed."""
    try:
        client = ArangoClient(hosts=arangodb_hosts)
        # Connect to system database first
        sys_db = client.db('_system', username=arangodb_user, password=arangodb_password)
        
        # Create database if not exists
        if not sys_db.has_database(db_name):
            sys_db.create_database(db_name)
            print(f"✓ Created database: {db_name}")
        else:
            print(f"Database {db_name} already exists")
        
        # Connect to target database
        db = client.db(db_name, username=arangodb_user, password=arangodb_password)
        print(f"✓ Connected to database: {db_name}")
        return db
    
    except Exception as e:
        print(f"✗ Error setting up ArangoDB connection. Is ArangoDB running on {arangodb_hosts}? Error: {e}")
        return None


def create_tcga_collections(db):
    """
    Create collections for TCGA datasets.
    Uses unified collections for each omic type (index + samples combined).
    No collection is dropped: incremental loading.
    """
    if db is None:
        return
    
    try:
        # Unified omic collections (each contains both index and sample data)
        omic_collections = [
            "GENE_EXPRESSION",
            "CNV",
            "MIRNA",
            "PROTEIN",
            "METHYLATION"
        ]
        
        # Other vertex collections
        other_collections = [
            "GENES",
            "PROJECTS",
            "SAMPLES",
            "CASES"
        ]
        
        all_collections = omic_collections + other_collections
        
        for cname in all_collections:
            if not db.has_collection(cname):
                db.create_collection(cname)
                print(f"✓ Created collection: {cname}")
            else:
                print(f"Collection {cname} already exists ({db.collection(cname).count()} documents)")
        
        # Create essential indexes
        try:
            # GENES collection index
            genes_col = db.collection("GENES")
            genes_col.add_index({
                "type": "hash",
                "fields": ["symbol"],
                "unique": False,
                "sparse": True
            })
            
            # SAMPLES collection index - use submitter_id instead of sample_id
            # Note: sample_id field is not present in samples.json, documents use _key
            samples_col = db.collection("SAMPLES")
            samples_col.add_index({
                "type": "hash",
                "fields": ["submitter_id"],
                "unique": False,
                "sparse": True
            })
            
            # Create indexes for unified omic collections
            for col_name in omic_collections:
                col = db.collection(col_name)
                
                # Index on data_type to distinguish between index and sample documents
                col.add_index({
                    "type": "hash",
                    "fields": ["data_type"],
                    "unique": False,
                    "sparse": False
                })
                
                # Index on sample_id for sample documents
                col.add_index({
                    "type": "hash",
                    "fields": ["sample_id"],
                    "unique": False,  # Not unique because index documents don't have sample_id
                    "sparse": True
                })
                
                # Index on cohort
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


# --- UTILITIES FOR LOADING AND INSERTING JSON ---

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
    Load a JSON file that can be:
    - Single JSON array
    - JSON Lines (one JSON document per line)
    
    Automatically sanitizes NaN/Inf values.
    """
    with open(path, "r", encoding="utf-8") as f:
        first = f.read(1)
        f.seek(0)
        
        if not first:
            return []
        
        if first == "[":
            # JSON array
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
    Incremental insertion.
    Assumes each document has a unique _key and uses insert_many(overwrite=False).
    This avoids overwriting existing documents.
    """
    from arango.exceptions import DocumentInsertError
    
    total = len(docs)
    if total == 0:
        print("  No documents to insert.")
        return
    
    inserted = 0
    skipped = 0  # Documents that already exist (duplicates)
    failed = 0
    errors_sample = []
    
    for i in range(0, total, BATCH_SIZE):
        batch = docs[i:i + BATCH_SIZE]
        try:
            # Use silent=False to get actual results
            result = collection.insert_many(batch, overwrite=False, silent=False)
            
            # Count actual insertions vs errors
            for r in result:
                # Check if result item is a dict (success) or an error object
                if isinstance(r, dict):
                    if r.get('_key') and not r.get('error'):
                        inserted += 1
                    else:
                        # This is an error in dict form
                        error_num = r.get('errorNum', 0)
                        if error_num == 1210:  # Unique constraint violated (duplicate)
                            skipped += 1
                        else:
                            failed += 1
                            if len(errors_sample) < 3:
                                errors_sample.append(r)
                elif hasattr(r, 'error_code'):
                    # It's a DocumentInsertError object
                    if r.error_code == 1210:  # Unique constraint violated (duplicate)
                        skipped += 1
                    else:
                        failed += 1
                        if len(errors_sample) < 3:
                            errors_sample.append(str(r))
                else:
                    # Unknown format, count as inserted
                    inserted += 1
            
            if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= total:
                print(f"  Processed {min(i + BATCH_SIZE, total)}/{total}")
                
        except Exception as e:
            failed += len(batch)
            print(f"  ⚠ Batch {i // BATCH_SIZE + 1} failed: {e}")
    
    print(f"  → Inserted {inserted}, skipped (already exist) {skipped}, failed {failed}")
    if errors_sample:
        print(f"  Sample errors: {errors_sample[:3]}")


# --- TCGA DATASETS IMPORT ---

def import_tcga_datasets(db, 
                         PARTIAL_LOAD_LIST=[],
                         REPLACE_EXISTING=False):
    """
    Load TCGA data into unified collections:
    - GENES (global, dedup on _key: insert without overwrite)
    - GENE_EXPRESSION, CNV, MIRNA, PROTEIN, METHYLATION 
      (unified collections containing both index and per-sample data)
    - PROJECTS, SAMPLES (metadata, incremental)
    
    The data_type field in each document indicates whether it's an index or sample data.
    """
    if db is None:
        print("✗ No database connection available")
        return
    
    if REPLACE_EXISTING:
        print("\n⚠ REPLACE_EXISTING is True: Dropping existing collections before load.")
        collections_to_drop = [
            "GENES", 
            "GENE_EXPRESSION", "CNV", "MIRNA", "PROTEIN", "METHYLATION",
            "PROJECTS", "SAMPLES"
        ]
        for col_name in collections_to_drop:
            if db.has_collection(col_name):
                db.delete_collection(col_name)
                print(f"  Dropped collection: {col_name}")
        create_tcga_collections(db)
    
    # 1. GENES (global, no overwrite; entries with same _key are not reloaded)
    genes_path = os.path.join(data_dir, "genes.json")
    if os.path.exists(genes_path):
        print(f"\n📂 Loading GENES from {genes_path}...")
        genes_docs = _load_json_lines(genes_path)
        genes_col = db.collection("GENES")
        _insert_incremental(genes_col, genes_docs, unique_key="_key")
    else:
        print("GENES file not found, skipping.")
    
    # 2. Load unified omic collections (index + samples combined)
    # Filter based on PARTIAL_LOAD_LIST if provided
    omic_types_to_load = OMIC_TYPES
    if len(PARTIAL_LOAD_LIST) > 0:
        omic_types_to_load = {
            k: v for k, v in OMIC_TYPES.items() 
            if any(item in k for item in PARTIAL_LOAD_LIST)
        }
    
    for omic_name, omic_config in omic_types_to_load.items():
        collection_name = omic_config["collection"]
        index_file = omic_config["index_file"]
        samples_file = omic_config["samples_file"]
        
        print(f"\n{'='*50}")
        print(f"📊 Loading {collection_name} (unified collection)")
        print(f"{'='*50}")
        
        col = db.collection(collection_name)
        
        # Load index data
        index_path = os.path.join(data_dir, index_file)
        if os.path.exists(index_path):
            print(f"\n  📂 Loading INDEX data from {index_file}...")
            index_docs = _load_json_lines(index_path)
            print(f"     Found {len(index_docs)} index documents")
            _insert_incremental(col, index_docs, unique_key="_key")
        else:
            print(f"  ⚠ {index_file} not found, skipping index data.")
        
        # Load samples data
        samples_path = os.path.join(data_dir, samples_file)
        if os.path.exists(samples_path):
            print(f"\n  📂 Loading SAMPLES data from {samples_file}...")
            samples_docs = _load_json_lines(samples_path)
            print(f"     Found {len(samples_docs)} sample documents")
            _insert_incremental(col, samples_docs, unique_key="_key")
        else:
            print(f"  ⚠ {samples_file} not found, skipping samples data.")
    
    # 3. PROJECTS, SAMPLES, CASES (metadata, global/incremental)
    meta_map = {
        "projects.json": "PROJECTS",
        "samples.json": "SAMPLES",
        "cases.json": "CASES" 
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
    print("=" * 60)
    print("ARANGODB IMPORT SCRIPT - TCGA DATASETS (v3 - Unified Collections)")
    print("=" * 60)
    
    # 1. Connection and Configuration
    db = setup_arangodb_connection(db_name)
    
    if db:
        # 2. Create TCGA collections (no drop, incremental)
        create_tcga_collections(db)
        
        # 3. Import TCGA Data
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
# Check collections and number of documents
if db:
    print("\n--- Collections in database ---")
    for col_name in db.collections():
        if col_name['name'].startswith('_'):
            continue
        collection = db.collection(col_name['name'])
        count = db.collection(col_name['name']).count()
        print(f"Collection: {col_name['name']}, Documents: {count}")

# %%


# collection fixer

# # Elimina la collection per rimuovere l'indice vecchio
# db.delete_collection("SAMPLES")

# # Ricrea la collection (questo creerà il nuovo indice corretto)
# create_tcga_collections(db)

# # Ricarica i dati
# path = "../data/arangodb_collections/TCGA-BRCA/samples.json"
# docs = _load_json_lines(path)
# print(f"Documenti da caricare: {len(docs)}")
# col = db.collection("SAMPLES")
# _insert_incremental(col, docs)

# # Verifica
# print(f"Documenti caricati: {col.count()}")

#%%
# Load Cases
# create_tcga_collections(db)
# path = "../data/arangodb_collections/TCGA-BRCA/cases.json"
# docs = _load_json_lines(path)
# print(f"Documenti da caricare: {len(docs)}")
# col = db.collection("CASES")
# _insert_incremental(col, docs)