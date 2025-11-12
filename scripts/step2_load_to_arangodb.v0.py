#%%
from arango import ArangoClient
import json
import traceback

# --- CONFIGURAZIONE ---
db_name = 'PKT_test10000'
arangodb_hosts = 'http://localhost:8529'
arangodb_user = 'root'
arangodb_password = 'avocadodb'
data_dir = "../temp_dir/"
NODES_FILE = f"{data_dir}nodes.json"
EDGES_FILE = f"{data_dir}edges.json"
BATCH_SIZE = 1000 # Dimensione del batch per le operazioni di inserimento

# --- FUNZIONI DI BASE DI ARANGODB ---

def setup_arangodb_connection(db_name):
    """Setup ArangoDB connection and create database if needed"""
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
        return db
        
    except Exception as e:
        print(f"✗ Error setting up ArangoDB connection. Is ArangoDB running on {arangodb_hosts}? Error: {e}")
        return None

def create_arangodb_collections(db_connection):
    """Create nodes (vertex) and edges (edge) collections"""
    if db_connection is None: return

    try:
        # Rimuovi le collezioni esistenti
        for collection_name in ['nodes', 'edges']:
            if db_connection.has_collection(collection_name):
                db_connection.delete_collection(collection_name)
                print(f"Deleted existing collection: {collection_name}")
        
        # Crea la collezione di vertici per i nodi
        db_connection.create_collection('nodes')
        print("✓ Created nodes (vertex) collection")
        
        # Crea la collezione di archi per le relazioni
        db_connection.create_collection('edges', edge=True)
        print("✓ Created edges (edge) collection")
        
        # Crea indici
        try:
            db_connection.collection('nodes').add_index({'fields': ['label'], 'type': 'hash'})
            db_connection.collection('nodes').add_index({'fields': ['entity_id'], 'type': 'hash'})
            db_connection.collection('edges').add_index({'fields': ['relationship'], 'type': 'hash'})
            print("✓ Created indexes for better query performance")
        except Exception as e:
            print(f"Warning: Could not create indexes: {e}")
            
    except Exception as e:
        print(f"✗ Error creating ArangoDB collections: {e}")
        sys.exit(1)


# --- FUNZIONE DI IMPORTAZIONE DATI ---

def import_data_to_arangodb(db_connection):
    """Load JSON files and import data into ArangoDB"""
    
    if db_connection is None: return

    try:
        print(f"\nLoading data from {NODES_FILE} and {EDGES_FILE}...")
        
        # 1. Carica Nodi
        with open(NODES_FILE, 'r', encoding='utf-8') as f:
            nodes_data = json.load(f)
        
        # Mappa per trovare la chiave (_key) ArangoDB dal node_id (URI)
        # Questo è fondamentale per costruire gli archi
        node_key_map = {node['node_id']: node['_key'] for node in nodes_data}
        
        nodes_collection = db_connection.collection('nodes')
        print(f"Inserting {len(nodes_data)} nodes into ArangoDB...")
        
        # Inserisci i nodi
        try:
            nodes_collection.insert_many(nodes_data, overwrite=True, silent=True)
            print(f"✓ Successfully inserted {len(nodes_data)} nodes.")
        except Exception as e:
            print(f"✗ Error inserting all nodes in batch. Trying individual insertion: {e}")
            for doc in nodes_data:
                try:
                    nodes_collection.insert(doc, overwrite=True)
                except Exception as single_error:
                    print(f"Failed to insert node {doc.get('entity_id', 'unknown')}: {single_error}")
            
        # 2. Carica Archi
        with open(EDGES_FILE, 'r', encoding='utf-8') as f:
            edges_data = json.load(f)
            
        edges_collection = db_connection.collection('edges')
        print(f"Preparing and inserting {len(edges_data)} edges into ArangoDB...")
        
        edge_docs = []
        for i, edge in enumerate(edges_data):
            # ArangoDB edge format richiede _from e _to
            source_key = node_key_map.get(edge['source_uri'])
            target_key = node_key_map.get(edge['target_uri'])
            
            if source_key and target_key:
                edge_doc = edge.copy()
                edge_doc['_from'] = f"nodes/{source_key}"
                edge_doc['_to'] = f"nodes/{target_key}"
                edge_doc['_key'] = f"edge_{i}" # Un _key univoco non basato sull'URI
                edge_docs.append(edge_doc)
            else:
                # Questo non dovrebbe succedere se lo script 1 è stato eseguito correttamente
                print(f"Warning: Missing key for edge {i} (Source URI: {edge['source_uri']}, Target URI: {edge['target_uri']}). Skipping.")

        # Inserisci gli archi
        try:
            edges_collection.insert_many(edge_docs, overwrite=True, silent=True)
            print(f"✓ Successfully inserted {len(edge_docs)} edges.")
        except Exception as e:
            print(f"✗ Error inserting all edges in batch. Trying individual insertion: {e}")
            for doc in edge_docs:
                try:
                    edges_collection.insert(doc, overwrite=True)
                except Exception as single_error:
                    print(f"Failed to insert edge {doc.get('_key', 'unknown')}: {single_error}")

        print(f"\nArangoDB import completed! Inserted {nodes_collection.count()} nodes and {edges_collection.count()} edges.")

    except Exception as e:
        print(f"✗ Fatal error during data import: {e}")
        traceback.print_exc()

def check_collections_data(db_connection):
    """Check and print the number of documents in each collection"""
    if db_connection is None: return

    try:
        for collection_name in ['nodes', 'edges']:
            if db_connection.has_collection(collection_name):
                collection = db_connection.collection(collection_name)
                count = collection.count()
                print(f"Collection '{collection_name}' has {count} documents.")
            else:
                print(f"Collection '{collection_name}' does not exist.")
    except Exception as e:
        print(f"✗ Error checking collections data: {e}")

#%%
db = setup_arangodb_connection(db_name)
#%%

# --- ESECUZIONE PRINCIPALE ---

if __name__ == "__main__":
    print("="*50)
    print("ARANGODB IMPORT SCRIPT (STEP 2)")
    print("="*50)
    
    # 1. Connessione e Configurazione
    db = setup_arangodb_connection(db_name)
    if db:
        create_arangodb_collections(db)
        
        # 2. Importazione Dati
        import_data_to_arangodb(db)
# %%
check_collections_data(db)