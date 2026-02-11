#%%
from arango import ArangoClient
import json
import traceback
import sys
import re

# --- CONFIGURAZIONE ---
db_name = 'PKT_test10000'
# db_name = 'PKT_transomics_v1'
arangodb_hosts = 'http://localhost:8529'
arangodb_user = 'root'
arangodb_password = 'avocadodb'
data_dir = "../temp_dir/"
NODES_FILE = f"{data_dir}nodes.json"
EDGES_FILE = f"{data_dir}edges.json"
BATCH_SIZE = 1000 # Dimensione del batch per le operazioni di inserimento

"""
The structure of the JSON data to be loaded are this

edge.json
[
  {
    "edge_id": "edge_1518276",
    "source_uri": "http://purl.obolibrary.org/obo/PR_Q9H609",
    "target_uri": "http://www.ncbi.nlm.nih.gov/gene/79177",
    "predicate_uri": "http://purl.obolibrary.org/obo/RO_0002204",
    "predicate_label": "gene product of",
    "predicate_class_code": "RO",
    "predicate_bioentity_type": "phenotype",
    "predicate_source": "Relation Ontology"
  },
  {
    "edge_id": "edge_7453030",
    "source_uri": "http://purl.obolibrary.org/obo/CLO_0023083",
    "target_uri": "http://purl.obolibrary.org/obo/CLO_0000652",
    "predicate_uri": "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
    "predicate_label": "type",
    "predicate_class_code": "RDF",
    "predicate_bioentity_type": "go",
    "predicate_source": "Resource Description Framework"
  },
  {

node.json
[
  {
    "_key": "PR_Q9H609",
    "uri": "http://purl.obolibrary.org/obo/PR_Q9H609",
    "namespace": "purl.obolibrary.org",
    "entity_id": "PR_Q9H609",
    "class_code": "PR",
    "label": "zinc finger protein 576 (human)",
    "bioentity_type": "protein",
    "description": "A zinc finger protein 576 that is encoded in the genome of human.",
    "synonym": "hZNF576|ZNF576",
    "source": "Protein Ontology",
    "source_type": "Ontology",
    "integer_id": 39530
  },
  {
    "_key": "79177",
    "uri": "http://www.ncbi.nlm.nih.gov/gene/79177",
    "namespace": "www.ncbi.nlm.nih.gov",
    "entity_id": "79177",
    "class_code": "EntrezID",
    "label": "ZNF576 (human)",
    "bioentity_type": "gene",
    "description": "A protein coding gene ZNF576 in human.",
    "synonym": "",
    "source": "NCBI Entrez Gene",
    "source_type": "Database",
    "integer_id": 216413
  },
"""


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
        print(f"✓ Connected to database: {db_name}")

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

# you must implement this function that is missing correctly, the main gial i sto load nodes edges and their propertis correcly i aragnodb 

def import_data_to_arangodb(db_connection):
    """Load JSON files and import data into ArangoDB"""
    if db_connection is None:
        print("✗ No database connection available")
        return
    
    try:
        
        # Load nodes
        print(f"\n📂 Loading nodes from {NODES_FILE}...")
        with open(NODES_FILE, 'r', encoding='utf-8') as f:
            nodes_data = json.load(f)

        # Load edges
        print(f"📂 Loading edges from {EDGES_FILE}...")
        with open(EDGES_FILE, 'r', encoding='utf-8') as f:
            edges_data = json.load(f)
        
        print(f"✓ Loaded {len(nodes_data)} nodes and {len(edges_data)} edges from JSON files")
        
        # Get collections
        nodes_collection = db_connection.collection('nodes')
        edges_collection = db_connection.collection('edges')

          
        # Import nodes in batches
        print(f"\n🔄 Importing nodes in batches of {BATCH_SIZE}...")
        nodes_imported = 0
        nodes_failed = 0
        
        for i in range(0, len(nodes_data), BATCH_SIZE):
            batch = nodes_data[i:i + BATCH_SIZE]
            try:
                result = nodes_collection.insert_many(batch, overwrite=False, silent=True)
                nodes_imported += len(batch)
                if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= len(nodes_data):
                    print(f"  Processed {min(i + BATCH_SIZE, len(nodes_data))}/{len(nodes_data)} nodes...")
            except Exception as e:
                nodes_failed += len(batch)
                print(f"  ⚠ Warning: Batch {i//BATCH_SIZE + 1} failed: {e}")
        
        print(f"✓ Imported {nodes_imported} nodes (failed: {nodes_failed})")
            
        # Import edges in batches
        print(f"\n🔄 Importing edges in batches of {BATCH_SIZE}...")
        edges_imported = 0
        edges_failed = 0
        
        for i in range(0, len(edges_data), BATCH_SIZE):
            batch = edges_data[i:i + BATCH_SIZE]
            
            # Transform edges to ArangoDB format while preserving ALL properties
            transformed_batch = []
            for edge in batch:
                # Extract entity_id from URIs (last component after /)
                source_key = edge['source_uri'].split('/')[-1].split('?t=')[-1]  # Handle query params
                target_key = edge['target_uri'].split('/')[-1].split('?t=')[-1]

                # Fix ENST identifiers if they are in the format Summary?t=ENST00000XXXXX
    
                
                # Create edge document with _from/_to and ALL original properties
                edge_doc = {
                    '_key': edge['edge_id'],
                    '_from': f"nodes/{source_key}",
                    '_to': f"nodes/{target_key}",
                    # Preserve ALL edge properties from the original JSON
                    **{k: v for k, v in edge.items() if k != 'edge_id'}
                }
                transformed_batch.append(edge_doc)
            
            try:
                result = edges_collection.insert_many(transformed_batch, overwrite=False, silent=True)
                edges_imported += len(transformed_batch)
                if (i // BATCH_SIZE + 1) % 10 == 0 or i + BATCH_SIZE >= len(edges_data):
                    print(f"  Processed {min(i + BATCH_SIZE, len(edges_data))}/{len(edges_data)} edges...")
            except Exception as e:
                edges_failed += len(transformed_batch)
                print(f"  ⚠ Warning: Batch {i//BATCH_SIZE + 1} failed: {e}")
        
        print(f"✓ Imported {edges_imported} edges (failed: {edges_failed})")
        print(f"\n✅ Import completed successfully!")
        
    except FileNotFoundError as e:
        print(f"✗ Error: Could not find data files: {e}")
        print(f"  Please ensure {NODES_FILE} and {EDGES_FILE} exist")
    except json.JSONDecodeError as e:
        print(f"✗ Error: Invalid JSON format: {e}")
    except Exception as e:
        print(f"✗ Unexpected error during import: {e}")
        traceback.print_exc()


def _extract_enst_id(uri):
    """Return ENST identifier found after the Summary query parameter in a URI."""
    if not uri:
        return None

    marker = "Summary?t="
    candidate = uri.split(marker, 1)[-1] if marker in uri else uri
    candidate = candidate.split('&', 1)[0]
    candidate = candidate.split('#', 1)[0]
    candidate = candidate.split('/', 1)[0]

    match = re.search(r"ENST\d+", candidate)
    return match.group(0) if match else None

def fix_summary_edges(db_connection):
    """Replace placeholder nodes/Summary references in edges with the actual ENST ids."""
    if db_connection is None:
        print("✗ No database connection available for ENST edge fix")
        return

    print("\n🔧 Fixing edges referencing nodes/Summary ...")
    query = """
    FOR edge IN edges
        FILTER edge._from == "nodes/Summary" OR edge._to == "nodes/Summary"
        RETURN {
            _key: edge._key,
            _from: edge._from,
            _to: edge._to,
            source_uri: edge.source_uri,
            target_uri: edge.target_uri
        }
    """

    edges_to_fix = list(db_connection.aql.execute(query))
    print(f"Found {len(edges_to_fix)} edges with placeholder references to nodes/Summary")
    if not edges_to_fix:
        print("✓ No edges with placeholder references detected")
        return

    edges_collection = db_connection.collection('edges')
    updated = 0
    skipped = 0
    from tqdm import tqdm
    for edge in tqdm(edges_to_fix, desc="Fixing edges"):
        update_doc = {'_key': edge['_key']}
        missing_identifier = False

        if edge['_from'] == 'nodes/Summary':
            enst_id = _extract_enst_id(edge.get('source_uri')) or _extract_enst_id(edge.get('target_uri'))
            if enst_id:
                update_doc['_from'] = f"nodes/{enst_id}"
            else:
                missing_identifier = True

        if edge['_to'] == 'nodes/Summary':
            enst_id = _extract_enst_id(edge.get('target_uri')) or _extract_enst_id(edge.get('source_uri'))
            if enst_id:
                update_doc['_to'] = f"nodes/{enst_id}"
            else:
                missing_identifier = True

        if len(update_doc) == 1 or missing_identifier:
            skipped += 1
            if missing_identifier:
                print(f"  ⚠ Skipped edge {edge['_key']} - unable to determine ENST identifier")
            continue

        try:
            edges_collection.update(update_doc)
            updated += 1
        except Exception as e:
            skipped += 1
            print(f"  ⚠ Failed to update edge {edge['_key']}: {e}")

    print(f"✓ Fixed {updated}/{len(edges_to_fix)} edges (skipped: {skipped})")

from arangodb_utils import *

#%%
# db = setup_arangodb_connection(db_name)

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

        # 2b. Correggi gli edge ENST con riferimenti a nodes/Summary
        # fix_summary_edges(db)
        
        # 3. Verifica finale
        print("\nFinal verification:")
        get_collections_data(db)

        # 4. Visualizza grafo esempio di query
        query = 'PR_Q9H609'
        subgraph = get_node_centric_graph(db, node_key=query)
        # subgraph = get_node_centric_graph(db, filters={'label': 'TP53'}, edge_limit=100)
        plot_subgraph(subgraph, layout='kamada_kawai', highlight_center=True, save_path=f'{query}_subgraph.png')
    else:
        print("✗ Failed to connect to ArangoDB. Please check if ArangoDB is running.")

#%%
if __name__ == "__main__":
    import arangodb_utils
    from importlib import reload
    reload(arangodb_utils)
    from arangodb_utils import get_node_centric_graph
    node = get_nodes_by_pattern(db, pattern="ENST00000362302", collection_name='nodes', property_name='entity_id')
    print(node)
#%%

# %%

#%%
# Fix ENST edges in the database (example for ENST00000362302)
"""
in nodes:
_id:nodes/ENST00000541341
_rev:_kneFK2K--F
_key:ENST00000541341
{
  "uri": "https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000541341",
  "namespace": "uswest.ensembl.org",
  "entity_id": "ENST00000541341",
  "class_code": "ENST",
  "label": "CHFR-219",
  "bioentity_type": "rna",
  "description": "Transcript CHFR-219 is classified as type 'retained_intron'.",
  "synonym": "",
  "source": "Ensembl Transcript",
  "source_type": "Database",
  "integer_id": 13
}

in edges:
_id:edges/edge_4
_rev:_kneF6ZS--B
_key:edge_4
_from:nodes/Summary
_to:nodes/SO_0000110
{
  "source_uri": "https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000541341",
  "target_uri": "http://purl.obolibrary.org/obo/SO_0000110",
  "predicate_uri": "http://www.w3.org/1999/02/22-rdf-syntax-ns#type",
  "predicate_label": "type",
  "predicate_class_code": "RDF",
  "predicate_bioentity_type": "go",
  "predicate_source": "Resource Description Framework"
}


questo _from:nodes/Summary va corretto nella collezione "edges" 
in ogni caso in cui _from o _to sia "nodes/Summary", questo va corretto, Summary va sostituito col valore reale di "ENST00000XXXXX"
il valore reale fi enst si trova nel relativo source o target_uri dopo "Summary?t="

scrivi un funzoine, essenziale funzionale che opera questa correzione nella collezione "edges"
"""

# %%
# 4. Visualizza il grafo
G = visualize_random_graph(
    db, 
    sample_size=200,  # Campiona 200 nodi per grafi grandi
    layout='spring'   # Usa layout spring (altre opzioni: 'circular', 'kamada_kawai')
)
# %%

# 3. Crea il grafo
graph = create_arango_graph(db, graph_name='PKT_graph')

# open_arango_web_viewer(graph_name='PKT_graph', db_name='PKT_test10000')

G = visualize_arango_graph(db, graph_name='PKT_graph', 
                           depth=10, 
                           limit=500, output_dir='.')


# %%
