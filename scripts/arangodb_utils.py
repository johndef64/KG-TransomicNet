#%%
import os
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
OMICS_PATH = "../data/omics/"


# --- FUNZIONI DI BASE DI ARANGODB ---

def setup_arangodb_connection(db_name = db_name):
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

# manage arando db databases
def list_databases():
    """List all databases in ArangoDB"""
    try:
        client = ArangoClient(hosts=arangodb_hosts)
        sys_db = client.db('_system', username=arangodb_user, password=arangodb_password)
        databases = sys_db.databases()
        print("Available databases:")
        for db_name in databases:
            print(f" - {db_name}")
    except Exception as e:
        print(f"✗ Error listing databases: {e}") 

def delete_database(db_name):
    """Delete a database in ArangoDB"""
    try:
        client = ArangoClient(hosts=arangodb_hosts)
        sys_db = client.db('_system', username=arangodb_user, password=arangodb_password)
        if sys_db.has_database(db_name):
            sys_db.delete_database(db_name)
            print(f"✓ Deleted database: {db_name}")
        else:
            print(f"⚠ Database '{db_name}' does not exist.")
    except Exception as e:
        print(f"✗ Error deleting database: {e}")

def dump_database(db_name, output_dir="./dumps/", include_system=False):
    """
    Dump an ArangoDB database to JSON files
    
    Args:
        db_name: Name of the database to dump
        output_dir: Directory where to save the dump files
        include_system: Whether to include system collections
    
    Returns:
        Path to the dump directory
    """
    try:
        import os
        import json
        from datetime import datetime
        
        # Setup connection
        db = setup_arangodb_connection(db_name)
        if db is None:
            print(f"✗ Failed to connect to database '{db_name}'")
            return None
        
        # Create dump directory with timestamp
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        dump_path = os.path.join(output_dir, f"{db_name}_{timestamp}")
        os.makedirs(dump_path, exist_ok=True)
        
        print(f"\n📦 Dumping database '{db_name}' to {dump_path}")
        
        # Get all collections
        collections = db.collections()
        
        # Metadata file
        metadata = {
            "database_name": db_name,
            "dump_date": timestamp,
            "collections": []
        }
        
        total_docs = 0
        
        for coll_info in collections:
            coll_name = coll_info['name']
            
            # Skip system collections if not requested
            if coll_name.startswith('_') and not include_system:
                continue
            
            print(f"\n  Processing collection: {coll_name}")
            
            try:
                collection = db.collection(coll_name)
                count = collection.count()
                
                print(f"    Documents: {count}")
                
                # Export collection data
                coll_file = os.path.join(dump_path, f"{coll_name}.json")
                
                # Get all documents
                aql = f"FOR doc IN {coll_name} RETURN doc"
                cursor = db.aql.execute(aql, batch_size=1000)
                
                documents = []
                batch_count = 0
                for doc in cursor:
                    documents.append(doc)
                    batch_count += 1
                    if batch_count % 1000 == 0:
                        print(f"    Exported {batch_count} documents...")
                
                # Save to file
                with open(coll_file, 'w', encoding='utf-8') as f:
                    json.dump(documents, f, indent=2, ensure_ascii=False)
                
                print(f"    ✓ Saved {len(documents)} documents to {coll_name}.json")
                
                # Add to metadata
                metadata["collections"].append({
                    "name": coll_name,
                    "type": coll_info['type'],
                    "count": count
                })
                
                total_docs += count
                
            except Exception as e:
                print(f"    ✗ Error exporting collection {coll_name}: {e}")
        
        # Check for graphs and export their definitions
        try:
            graphs = db.graphs()
            if graphs:
                print(f"\n  Processing {len(graphs)} graph(s)")
                graph_definitions = []
                
                for graph in graphs:
                    graph_info = {
                        "name": graph['name'],
                        "edge_definitions": graph['edge_definitions'],
                        "orphan_collections": graph.get('orphan_collections', [])
                    }
                    graph_definitions.append(graph_info)
                    print(f"    ✓ Exported graph: {graph['name']}")
                
                # Save graph definitions
                graphs_file = os.path.join(dump_path, "_graphs.json")
                with open(graphs_file, 'w', encoding='utf-8') as f:
                    json.dump(graph_definitions, f, indent=2)
                
                metadata["graphs"] = graph_definitions
        except Exception as e:
            print(f"  ✗ Error exporting graphs: {e}")
        
        # Save metadata
        metadata_file = os.path.join(dump_path, "_metadata.json")
        with open(metadata_file, 'w', encoding='utf-8') as f:
            json.dump(metadata, f, indent=2)
        
        print(f"\n✅ Database dump completed!")
        print(f"   Total documents: {total_docs}")
        print(f"   Total collections: {len(metadata['collections'])}")
        print(f"   Dump location: {dump_path}")
        
        return dump_path
        
    except Exception as e:
        print(f"✗ Error dumping database: {e}")
        traceback.print_exc()
        return None

def restore_database_from_dump(dump_path, target_db_name=None, overwrite=False):
    """
    Restore an ArangoDB database from dump files
    
    Args:
        dump_path: Path to the dump directory
        target_db_name: Name for the restored database (None = use original name)
        overwrite: If True, delete existing database before restoring
    
    Returns:
        Database connection object
    """
    try:
        import os
        import json
        
        # Check if dump path exists
        if not os.path.exists(dump_path):
            print(f"✗ Dump path not found: {dump_path}")
            return None
        
        # Load metadata
        metadata_file = os.path.join(dump_path, "_metadata.json")
        if not os.path.exists(metadata_file):
            print(f"✗ Metadata file not found: {metadata_file}")
            return None
        
        with open(metadata_file, 'r', encoding='utf-8') as f:
            metadata = json.load(f)
        
        # Determine database name
        original_db_name = metadata['database_name']
        db_name = target_db_name if target_db_name else original_db_name
        
        print(f"\n📥 Restoring database '{db_name}' from {dump_path}")
        print(f"   Original database: {original_db_name}")
        print(f"   Dump date: {metadata['dump_date']}")
        
        # Check if database exists
        client = ArangoClient(hosts=arangodb_hosts)
        sys_db = client.db('_system', username=arangodb_user, password=arangodb_password)
        
        if sys_db.has_database(db_name):
            if overwrite:
                print(f"   ⚠ Deleting existing database '{db_name}'...")
                sys_db.delete_database(db_name)
            else:
                print(f"   ✗ Database '{db_name}' already exists. Use overwrite=True to replace.")
                return None
        
        # Create database
        sys_db.create_database(db_name)
        print(f"   ✓ Created database: {db_name}")
        
        # Connect to new database
        db = client.db(db_name, username=arangodb_user, password=arangodb_password)
        
        # Restore collections
        total_restored = 0
        
        for coll_info in metadata['collections']:
            coll_name = coll_info['name']
            coll_type = coll_info['type']
            expected_count = coll_info['count']
            
            print(f"\n  Restoring collection: {coll_name}")
            print(f"    Type: {'edge' if coll_type == 3 else 'document'}")
            print(f"    Expected documents: {expected_count}")
            
            try:
                # Create collection
                if coll_type == 3:  # Edge collection
                    collection = db.create_collection(coll_name, edge=True)
                else:  # Document collection
                    collection = db.create_collection(coll_name)
                
                # Load data file
                coll_file = os.path.join(dump_path, f"{coll_name}.json")
                if not os.path.exists(coll_file):
                    print(f"    ⚠ Data file not found: {coll_file}")
                    continue
                
                with open(coll_file, 'r', encoding='utf-8') as f:
                    documents = json.load(f)
                
                # Insert documents in batches
                batch_size = 1000
                inserted = 0
                
                for i in range(0, len(documents), batch_size):
                    batch = documents[i:i+batch_size]
                    collection.insert_many(batch, overwrite=True)
                    inserted += len(batch)
                    if inserted % 5000 == 0:
                        print(f"    Inserted {inserted}/{len(documents)} documents...")
                
                actual_count = collection.count()
                print(f"    ✓ Restored {actual_count} documents")
                total_restored += actual_count
                
                if actual_count != expected_count:
                    print(f"    ⚠ Warning: Expected {expected_count} but restored {actual_count}")
                
            except Exception as e:
                print(f"    ✗ Error restoring collection {coll_name}: {e}")
        
        # Restore graphs if present
        if 'graphs' in metadata and metadata['graphs']:
            print(f"\n  Restoring {len(metadata['graphs'])} graph(s)")
            
            for graph_info in metadata['graphs']:
                try:
                    graph_name = graph_info['name']
                    
                    # Check if graph already exists
                    if db.has_graph(graph_name):
                        print(f"    ⚠ Graph '{graph_name}' already exists, skipping...")
                        continue
                    
                    # Create graph
                    db.create_graph(
                        name=graph_name,
                        edge_definitions=graph_info['edge_definitions'],
                        orphan_collections=graph_info.get('orphan_collections', [])
                    )
                    print(f"    ✓ Restored graph: {graph_name}")
                    
                except Exception as e:
                    print(f"    ✗ Error restoring graph {graph_name}: {e}")
        
        print(f"\n✅ Database restore completed!")
        print(f"   Total documents restored: {total_restored}")
        print(f"   Total collections: {len(metadata['collections'])}")
        
        return db
        
    except Exception as e:
        print(f"✗ Error restoring database: {e}")
        traceback.print_exc()
        return None

def list_dumps(dumps_dir="./dumps/"):
    """
    List all available database dumps
    
    Args:
        dumps_dir: Directory containing dumps
    
    Returns:
        List of dump information dictionaries
    """
    try:
        import os
        import json
        
        if not os.path.exists(dumps_dir):
            print(f"✗ Dumps directory not found: {dumps_dir}")
            return []
        
        dumps = []
        
        for item in os.listdir(dumps_dir):
            item_path = os.path.join(dumps_dir, item)
            if os.path.isdir(item_path):
                metadata_file = os.path.join(item_path, "_metadata.json")
                if os.path.exists(metadata_file):
                    with open(metadata_file, 'r', encoding='utf-8') as f:
                        metadata = json.load(f)
                    
                    dump_info = {
                        "path": item_path,
                        "database_name": metadata['database_name'],
                        "dump_date": metadata['dump_date'],
                        "collections_count": len(metadata['collections']),
                        "graphs_count": len(metadata.get('graphs', []))
                    }
                    dumps.append(dump_info)
        
        if dumps:
            print(f"📦 Found {len(dumps)} dump(s) in {dumps_dir}:")
            for dump in sorted(dumps, key=lambda x: x['dump_date'], reverse=True):
                print(f"\n  Database: {dump['database_name']}")
                print(f"    Date: {dump['dump_date']}")
                print(f"    Collections: {dump['collections_count']}")
                print(f"    Graphs: {dump['graphs_count']}")
                print(f"    Path: {dump['path']}")
        else:
            print(f"No dumps found in {dumps_dir}")
        
        return dumps
        
    except Exception as e:
        print(f"✗ Error listing dumps: {e}")
        traceback.print_exc()
        return []




# --- FUNZIONI AGGIUNTIVE ---
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


# funzioni per caricare multi omics datasets
# setup_arangodb_connection()
#%%

# --- FUNZIONI DI RETRIEVAL NODI ---

def get_nodes_by_ids(db, collection_name, ids_list):
    """
    Recupera nodi da una collezione in base a una lista di IDs
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        ids_list: Lista di IDs dei nodi da recuperare
    
    Returns:
        Lista di nodi che corrispondono agli IDs
    """
    try:
        aql = f"""
        FOR doc IN {collection_name}
            FILTER doc._key IN @ids
            RETURN doc
        """
        cursor = db.aql.execute(aql, bind_vars={'ids': ids_list})
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes with specified IDs")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def get_nodes_by_property(db, collection_name, property_name, property_value):
    """
    Recupera nodi da una collezione in base a una proprietà specifica
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà da filtrare
        property_value: Valore della proprietà
    
    Returns:
        Lista di nodi che corrispondono al criterio
    """
    try:
        aql = f"""
        FOR doc IN {collection_name}
            FILTER doc.{property_name} == @value
            RETURN doc
        """
        cursor = db.aql.execute(aql, bind_vars={'value': property_value})
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes with {property_name}={property_value}")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def get_nodes_by_multiple_properties(db, collection_name, filters):
    """
    Recupera nodi in base a multipli criteri
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        filters: Dizionario con coppie property_name: property_value
    
    Returns:
        Lista di nodi che corrispondono a tutti i criteri
    """
    try:
        filter_conditions = []
        bind_vars = {}
        
        for i, (prop, value) in enumerate(filters.items()):
            filter_conditions.append(f"doc.{prop} == @value{i}")
            bind_vars[f'value{i}'] = value
        
        filter_str = " AND ".join(filter_conditions)
        
        aql = f"""
        FOR doc IN {collection_name}
            FILTER {filter_str}
            RETURN doc
        """
        
        cursor = db.aql.execute(aql, bind_vars=bind_vars)
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes matching all criteria")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def get_nodes_by_property_range(db, collection_name, property_name, min_value=None, max_value=None):
    """
    Recupera nodi con proprietà in un range specificato
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà numerica
        min_value: Valore minimo (opzionale)
        max_value: Valore massimo (opzionale)
    
    Returns:
        Lista di nodi nel range specificato
    """
    try:
        conditions = []
        bind_vars = {}
        
        if min_value is not None:
            conditions.append(f"doc.{property_name} >= @min_val")
            bind_vars['min_val'] = min_value
        
        if max_value is not None:
            conditions.append(f"doc.{property_name} <= @max_val")
            bind_vars['max_val'] = max_value
        
        filter_str = " AND ".join(conditions) if conditions else "true"
        
        aql = f"""
        FOR doc IN {collection_name}
            FILTER {filter_str}
            RETURN doc
        """
        
        cursor = db.aql.execute(aql, bind_vars=bind_vars)
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes in range")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def get_nodes_by_pattern(db, collection_name, property_name, pattern, case_sensitive=False):
    """
    Recupera nodi con proprietà che corrispondono a un pattern (LIKE)
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà stringa
        pattern: Pattern da cercare (es. "%gene%")
        case_sensitive: Se True, ricerca case-sensitive
    
    Returns:
        Lista di nodi che corrispondono al pattern

    Patters wilde cards in AQL:
        % - rappresenta zero o più caratteri
        _ - rappresenta un singolo carattere 

    Patters examples:
        "%gene%" - contiene 'gene'
        "gene%"  - inizia con 'gene'
        "%gene"  - termina con 'gene'
        "_gen%"  - inizia con qualsiasi carattere, seguito da 'gen'
        "%gen_"  - termina con qualsiasi carattere, preceduto da 'gen'
        "_gen_"  - contiene 'gen' con qualsiasi carattere prima e dopo
        "%gen%123" - contiene 'gen' seguito da '123'
        "_gen_123" - contiene 'gen' con qualsiasi carattere prima e dopo, seguito da '123'
    """
    try:
        like_func = "LIKE" if case_sensitive else "LIKE"
        
        aql = f"""
        FOR doc IN {collection_name}
            FILTER {like_func}(doc.{property_name}, @pattern, true)
            RETURN doc
        """
        
        cursor = db.aql.execute(aql, bind_vars={'pattern': pattern})
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes matching pattern '{pattern}'")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []



def get_nodes_by_property_list(db, collection_name, property_name, values_list):
    """
    Recupera nodi dove la proprietà è presente in una lista di valori
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà
        values_list: Lista di valori da cercare
    
    Returns:
        Lista di nodi con proprietà nella lista
    """
    try:
        aql = f"""
        FOR doc IN {collection_name}
            FILTER doc.{property_name} IN @values
            RETURN doc
        """
        
        cursor = db.aql.execute(aql, bind_vars={'values': values_list})
        results = [doc for doc in cursor]
        print(f"✓ Found {len(results)} nodes with {property_name} in list")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def get_all_nodes_with_properties(db, collection_name, properties_list=None):
    """
    Recupera tutti i nodi, opzionalmente proiettando solo alcune proprietà
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        properties_list: Lista di proprietà da restituire (None = tutte)
    
    Returns:
        Lista di nodi
    """
    try:
        if properties_list:
            return_obj = "{" + ", ".join([f"{p}: doc.{p}" for p in properties_list]) + "}"
            aql = f"""
            FOR doc IN {collection_name}
                RETURN {return_obj}
            """
        else:
            aql = f"""
            FOR doc IN {collection_name}
                RETURN doc
            """
        
        cursor = db.aql.execute(aql)
        results = [doc for doc in cursor]
        print(f"✓ Retrieved {len(results)} nodes from {collection_name}")
        return results
    except Exception as e:
        print(f"✗ Error retrieving nodes: {e}")
        traceback.print_exc()
        return []

def count_nodes_by_property(db, collection_name, property_name, property_value):
    """
    Conta i nodi che corrispondono a una proprietà specifica
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà
        property_value: Valore della proprietà
    
    Returns:
        Numero di nodi che corrispondono
    """
    try:
        aql = f"""
        FOR doc IN {collection_name}
            FILTER doc.{property_name} == @value
            COLLECT WITH COUNT INTO length
            RETURN length
        """
        
        cursor = db.aql.execute(aql, bind_vars={'value': property_value})
        count = next(cursor, 0)
        print(f"✓ Count: {count} nodes with {property_name}={property_value}")
        return count
    except Exception as e:
        print(f"✗ Error counting nodes: {e}")
        traceback.print_exc()
        return 0

def get_distinct_property_values(db, collection_name, property_name):
    """
    Ottiene tutti i valori distinti di una proprietà
    
    Args:
        db: Database ArangoDB
        collection_name: Nome della collezione
        property_name: Nome della proprietà
    
    Returns:
        Lista di valori distinti
    """
    try:
        aql = f"""
        FOR doc IN {collection_name}
            RETURN DISTINCT doc.{property_name}
        """
        
        cursor = db.aql.execute(aql)
        values = [val for val in cursor if val is not None]
        print(f"✓ Found {len(values)} distinct values for {property_name}")
        return values
    except Exception as e:
        print(f"✗ Error getting distinct values: {e}")
        traceback.print_exc()
        return []



# --- FUNZIONI GRAFI ---
def create_arango_graph(db_connection, graph_name='PKT_graph'):
    """Create an ArangoDB named graph from nodes and edges collections"""
    if db_connection is None:
        print("✗ No database connection available")
        return None
    
    try:
        # Check if graph already exists
        if db_connection.has_graph(graph_name):
            print(f"⚠ Graph '{graph_name}' already exists. Deleting it...")
            db_connection.delete_graph(graph_name, drop_collections=False)
        
        # Define edge definitions for the graph
        edge_definitions = [
            {
                'edge_collection': 'edges',
                'from_vertex_collections': ['nodes'],
                'to_vertex_collections': ['nodes']
            }
        ]
        
        # Create the graph
        print(f"🔨 Creating graph '{graph_name}'...")
        graph = db_connection.create_graph(
            name=graph_name,
            edge_definitions=edge_definitions
        )
        
        print(f"✅ Graph '{graph_name}' created successfully!")
        print(f"   - Vertex collections: nodes")
        print(f"   - Edge collections: edges")
        
        # Get some statistics
        nodes_count = db_connection.collection('nodes').count()
        edges_count = db_connection.collection('edges').count()
        print(f"\n📊 Graph statistics:")
        print(f"   - Total nodes: {nodes_count}")
        print(f"   - Total edges: {edges_count}")
        
        return graph
        
    except Exception as e:
        print(f"✗ Error creating graph: {e}")
        traceback.print_exc()
        return None

def visualize_random_graph(db_connection, 
                          sample_size=100, layout='spring'):
    """Visualize ArangoDB graph using NetworkX and matplotlib
    
    Args:
        db_connection: ArangoDB database connection
        graph_name: Name of the graph to visualize
        sample_size: Number of nodes to sample (for large graphs)
        layout: Layout algorithm ('spring', 'circular', 'kamada_kawai', 'random')
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        from collections import defaultdict
        
        if db_connection is None:
            print("✗ No database connection available")
            return
        
        print(f"\n📊 Visualizing graph ...")
        
        # Get collections
        nodes_collection = db_connection.collection('nodes')
        edges_collection = db_connection.collection('edges')
        
        total_nodes = nodes_collection.count()
        total_edges = edges_collection.count()
        
        print(f"   Total nodes: {total_nodes}, Total edges: {total_edges}")
        
        # Sample nodes if graph is too large
        if total_nodes > sample_size:
            print(f"   Sampling {sample_size} nodes for visualization...")
            aql_nodes = f"""
                FOR node IN nodes
                LIMIT {sample_size}
                RETURN node
            """
        else:
            aql_nodes = """
                FOR node IN nodes
                RETURN node
            """
        
        # Fetch nodes
        cursor_nodes = db_connection.aql.execute(aql_nodes)
        nodes_list = list(cursor_nodes)
        sampled_keys = {node['_key'] for node in nodes_list}
        
        # Fetch edges connected to sampled nodes
        aql_edges = f"""
            FOR edge IN edges
            FILTER SPLIT(edge._from, '/')[1] IN {list(sampled_keys)}
            AND SPLIT(edge._to, '/')[1] IN {list(sampled_keys)}
            RETURN edge
        """
        cursor_edges = db_connection.aql.execute(aql_edges)
        edges_list = list(cursor_edges)
        
        print(f"   Loaded {len(nodes_list)} nodes and {len(edges_list)} edges")
        
        # Create NetworkX graph
        G = nx.DiGraph()
        
        # Add nodes with attributes
        node_types = defaultdict(list)
        for node in nodes_list:
            node_id = node['_key']
            G.add_node(node_id, **{k: v for k, v in node.items() if k not in ['_key', '_id', '_rev']})
            node_label = node.get('label', 'unknown')
            node_types[node_label].append(node_id)
        
        # Add edges with attributes
        edge_types = defaultdict(int)
        for edge in edges_list:
            source = edge['_from'].split('/')[1]
            target = edge['_to'].split('/')[1]
            relationship = edge.get('relationship', 'unknown')
            G.add_edge(source, target, **{k: v for k, v in edge.items() 
                                         if k not in ['_key', '_id', '_rev', '_from', '_to']})
            edge_types[relationship] += 1
        
        # Print graph statistics
        print(f"\n📈 Graph statistics:")
        print(f"   - Nodes by type:")
        for node_type, nodes in sorted(node_types.items(), key=lambda x: len(x[1]), reverse=True):
            print(f"     • {node_type}: {len(nodes)}")
        print(f"   - Edges by relationship:")
        for edge_type, count in sorted(edge_types.items(), key=lambda x: x[1], reverse=True):
            print(f"     • {edge_type}: {count}")
        
        # Create visualization
        plt.figure(figsize=(16, 12))
        
        # Choose layout
        if layout == 'spring':
            pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        else:
            pos = nx.random_layout(G)
        
        # Define colors for different node types
        color_map = {
            'protein': '#FF6B6B',
            'gene': '#4ECDC4',
            'chemical': '#95E1D3',
            'disease': '#F38181',
            'pathway': '#FEE191',
            'transcript': '#AA96DA',
            'metabolite': '#FCBAD3',
            'unknown': '#CCCCCC'
        }
        
        # Color nodes by type
        node_colors = []
        for node in G.nodes():
            node_data = nodes_list[[n['_key'] for n in nodes_list].index(node)]
            node_label = node_data.get('label', 'unknown')
            node_colors.append(color_map.get(node_label, '#CCCCCC'))
        
        # Draw the graph
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                              node_size=300, alpha=0.8)
        nx.draw_networkx_edges(G, pos, edge_color='gray', 
                              arrows=True, arrowsize=10, 
                              alpha=0.3, width=0.5)
        
        # Add labels for small graphs
        if len(G.nodes()) < 50:
            labels = {node: node[:15] for node in G.nodes()}  # Truncate long labels
            nx.draw_networkx_labels(G, pos, labels, font_size=8)
        
        plt.title(f"Knowledge Graph Visualization: \n"
                 f"({len(G.nodes())} nodes, {len(G.edges())} edges)", 
                 fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()
        
        # Save figure
        output_file = f"random_graph_visualization.png"
        # plt.savefig(output_file, dpi=300, bbox_inches='tight')
        # print(f"\n✅ Visualization saved to: {output_file}")
        
        plt.show()
        
        return G
        
    except ImportError as e:
        print(f"✗ Missing required libraries. Install with:")
        print(f"   pip install networkx matplotlib")
    except Exception as e:
        print(f"✗ Error visualizing graph: {e}")
        traceback.print_exc()
        return None

def open_arango_web_viewer(graph_name='PKT_graph', host='localhost', port=8529, db_name='PKT_test10000'):
    """Open ArangoDB web interface graph viewer in browser"""
    import webbrowser
    
    # Construct URL for graph viewer
    url = f"http://{host}:{port}/_db/{db_name}/_admin/aardvark/index.html#graph/{graph_name}"
    
    print(f"🌐 Opening ArangoDB graph viewer in browser...")
    print(f"   URL: {url}")
    print(f"\n📌 Instructions:")
    print(f"   1. Login with your credentials (root/avocadodb)")
    print(f"   2. The graph viewer will show your '{graph_name}' graph")
    print(f"   3. You can:")
    print(f"      - Click nodes to explore")
    print(f"      - Drag nodes to rearrange")
    print(f"      - Change depth and limit in settings")
    print(f"      - Take screenshots")
    print(f"      - Search specific nodes")
    
    webbrowser.open(url)

def visualize_arango_graph(db_connection, graph_name='PKT_graph', 
                          start_node=None, depth=2, limit=250, 
                          layout='spring', output_dir='./'):
    """Visualize a saved ArangoDB graph by querying it
    
    Args:
        db_connection: ArangoDB database connection
        graph_name: Name of the saved graph in ArangoDB
        start_node: Starting node key for traversal (None = random)
        depth: Traversal depth
        limit: Maximum nodes to retrieve
        layout: Layout algorithm
        output_dir: Output directory for saving visualization
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        from collections import defaultdict
        import os
        import re
        
        if db_connection is None:
            print("✗ No database connection available")
            return
        
        # Check if graph exists
        if not db_connection.has_graph(graph_name):
            print(f"✗ Graph '{graph_name}' does not exist in database")
            print(f"   Available graphs: {db_connection.graphs()}")
            return
        
        print(f"\n📊 Retrieving graph '{graph_name}' from ArangoDB...")
        
        # Get the saved graph
        graph = db_connection.graph(graph_name)
        
        # Get graph info
        edge_definitions = graph.edge_definitions()
        vertex_collections = graph.vertex_collections()
        
        print(f"   Vertex collections: {vertex_collections}")
        print(f"   Edge definitions: {len(edge_definitions)}")
        
        # If no start node, pick a random one
        if start_node is None:
            aql_random = f"""
                FOR v IN {vertex_collections[0]}
                LIMIT 1
                RETURN v._key
            """
            cursor = db_connection.aql.execute(aql_random)
            start_node = list(cursor)[0]
        
        print(f"   Start node: {start_node}")
        print(f"   Depth: {depth}, Limit: {limit}")
        
        # Traverse the saved graph
        aql_traverse = f"""
            FOR v, e, p IN 1..{depth} ANY '{vertex_collections[0]}/{start_node}' GRAPH '{graph_name}'
            OPTIONS {{uniqueVertices: 'global', bfs: true}}
            LIMIT {limit}
            RETURN {{vertex: v, edge: e, path: p}}
        """
        
        cursor = db_connection.aql.execute(aql_traverse)
        traversal_results = list(cursor)
        
        print(f"   Retrieved {len(traversal_results)} traversal results")
        
        # Extract unique nodes and edges
        nodes_dict = {}
        edges_list = []
        
        # Add start node
        start_doc = graph.vertex_collection(vertex_collections[0]).get(start_node)
        nodes_dict[start_node] = start_doc
        
        for result in traversal_results:
            vertex = result['vertex']
            edge = result['edge']
            
            if vertex:
                nodes_dict[vertex['_key']] = vertex
            
            if edge:
                edges_list.append(edge)
        
        print(f"   Unique nodes: {len(nodes_dict)}, Unique edges: {len(edges_list)}")
        
        # Create NetworkX graph
        G = nx.DiGraph()
        
        # Add nodes with attributes
        node_types = defaultdict(list)
        for node_key, node_data in nodes_dict.items():
            G.add_node(node_key, **{k: v for k, v in node_data.items() 
                                    if k not in ['_key', '_id', '_rev']})
            node_label = node_data.get('label', 'unknown')
            node_types[node_label].append(node_key)
        
        # Add edges with attributes
        edge_types = defaultdict(int)
        for edge in edges_list:
            source = edge['_from'].split('/')[1]
            target = edge['_to'].split('/')[1]
            if source in nodes_dict and target in nodes_dict:
                relationship = edge.get('relationship', 'unknown')
                G.add_edge(source, target, **{k: v for k, v in edge.items() 
                                             if k not in ['_key', '_id', '_rev', '_from', '_to']})
                edge_types[relationship] += 1
        
        # Print statistics
        print(f"\n📈 Graph statistics:")
        print(f"   - Nodes by type (top 10):")
        for node_type, nodes in sorted(node_types.items(), key=lambda x: len(x[1]), reverse=True)[:10]:
            print(f"     • {node_type}: {len(nodes)}")
        print(f"   - Edges by relationship:")
        for edge_type, count in sorted(edge_types.items(), key=lambda x: x[1], reverse=True)[:10]:
            print(f"     • {edge_type}: {count}")
        
        # Visualize
        plt.figure(figsize=(16, 12))
        
        # Layout
        if layout == 'spring':
            pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        else:
            pos = nx.random_layout(G)
        
        # Colors
        color_map = {
            'protein': '#FF6B6B',
            'gene': '#4ECDC4',
            'chemical': '#95E1D3',
            'disease': '#F38181',
            'pathway': '#FEE191',
            'transcript': '#AA96DA',
            'metabolite': '#FCBAD3',
            'unknown': '#CCCCCC'
        }
        
        node_colors = []
        for node in G.nodes():
            node_data = nodes_dict[node]
            node_label = node_data.get('label', 'unknown')
            node_colors.append(color_map.get(node_label, '#CCCCCC'))
        
        # Draw
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                              node_size=300, alpha=0.8)
        nx.draw_networkx_edges(G, pos, edge_color='gray', 
                              arrows=True, arrowsize=10, 
                              alpha=0.3, width=0.5)
        
        if len(G.nodes()) < 50:
            labels = {node: node[:15] for node in G.nodes()}
            nx.draw_networkx_labels(G, pos, labels, font_size=8)
        
        plt.title(f"Saved Graph: {graph_name}\n({len(G.nodes())} nodes, {len(G.edges())} edges)", 
                 fontsize=16, fontweight='bold')
        plt.axis('off')
        plt.tight_layout()
        
        # Save
        # safe_name = re.sub(r'[<>:"/\\|?*]', '_', graph_name)
        # output_file = os.path.join(output_dir, f"{safe_name}_from_db.png")
        # os.makedirs(output_dir, exist_ok=True)
        
        # plt.savefig(output_file, dpi=300, bbox_inches='tight')
        # print(f"\n✅ Visualization saved to: {output_file}")
        
        plt.show()
        
        return G
        
    except Exception as e:
        print(f"✗ Error: {e}")
        traceback.print_exc()
        return None



#%%
list_databases()
#%%
# rename_database('PKT_test10000', 'PKT_transomics_v1')
# #%%
# name = "PKT_transomics_v1"
# delete_database(name)
# #%%
# dump_database("PKT_test10000", include_system=False)
# #%%
# restore_database_from_dump("./dumps/PKT_test10000_20240605_123456", target_db_name="PKT_transomics_v1", overwrite=True)
#%%

# --- MAIN DI TEST ---

if __name__ == "__main__":
    print("="*80)
    print("TESTING ARANGODB RETRIEVAL FUNCTIONS")
    print("="*80)
    
    # Connessione al database
    print("\n1. Testing database connection...")
    db = setup_arangodb_connection()
    
    if db is None:
        print("Failed to connect to database. Exiting.")
        exit(1)
    
    # Lista delle collezioni disponibili
    print("\n2. Available collections:")
    collections = db.collections()
    for coll in collections:
        if not coll['name'].startswith('_'):  # Escludi collezioni di sistema
            print(f"   - {coll['name']} (type: {coll['type']})")
    
    # Seleziona una collezione di nodi per i test
    # Modifica questo nome in base alle tue collezioni
    test_collection = "nodes"  # Cambia con il nome della tua collezione
    
    if not db.has_collection(test_collection):
        print(f"\n⚠ Collection '{test_collection}' not found. Please update test_collection variable.")
        print("Available collections:", [c['name'] for c in collections if not c['name'].startswith('_')])
        exit(1)
    
    print(f"\n3. Using collection: {test_collection}")
    
    # Test 1: Recupera tutti i nodi (limitato a 5 per il test)
    print("\n" + "="*80)
    print("TEST 1: Get all nodes (limited to 5)")
    print("="*80)
    try:
        aql = f"FOR doc IN {test_collection} LIMIT 5 RETURN doc"
        cursor = db.aql.execute(aql)
        sample_nodes = [doc for doc in cursor]
        print(f"Sample of {len(sample_nodes)} nodes:")
        for i, node in enumerate(sample_nodes, 1):
            print(f"\nNode {i}:")
            for key, value in node.items():
                if key != '_rev':  # Skip revision field
                    print(f"  {key}: {value}")
    except Exception as e:
        print(f"Error: {e}")
        sample_nodes = []
    
    # Test 2: Ottieni valori distinti di una proprietà
    print("\n" + "="*80)
    print("TEST 2: Get distinct property values")
    print("="*80)
    if sample_nodes:
        # Trova la prima proprietà disponibile (escluse quelle di sistema)
        test_property = None
        for key in sample_nodes[0].keys():
            if not key.startswith('_') and key not in ['id', 'key']:
                test_property = key
                break
        
        if test_property:
            distinct_values = get_distinct_property_values(db, test_collection, test_property)
            print(f"First 10 distinct values of '{test_property}':")
            for val in distinct_values[:10]:
                print(f"  - {val}")
    
    # Test 3: Cerca nodi per proprietà specifica
    print("\n" + "="*80)
    print("TEST 3: Get nodes by property")
    print("="*80)
    if sample_nodes and test_property:
        test_value = sample_nodes[0].get(test_property)
        if test_value:
            results = get_nodes_by_property(db, test_collection, test_property, test_value)
            print(f"Found {len(results)} nodes with {test_property}='{test_value}'")
    
    # Test 4: Cerca nodi per multipli criteri
    print("\n" + "="*80)
    print("TEST 4: Get nodes by multiple properties")
    print("="*80)
    if sample_nodes and len(sample_nodes[0].keys()) >= 2:
        filters = {}
        for key in list(sample_nodes[0].keys())[:2]:
            if not key.startswith('_'):
                filters[key] = sample_nodes[0][key]
        
        if filters:
            print(f"Searching with filters: {filters}")
            results = get_nodes_by_multiple_properties(db, test_collection, filters)
            if results:
                print(f"First result: {results[0]}")
    
    # Test 5: Conta nodi per proprietà
    print("\n" + "="*80)
    print("TEST 5: Count nodes by property")
    print("="*80)
    if sample_nodes and test_property:
        test_value = sample_nodes[0].get(test_property)
        if test_value:
            count = count_nodes_by_property(db, test_collection, test_property, test_value)
    
    # Test 6: Cerca nodi per lista di valori
    print("\n" + "="*80)
    print("TEST 6: Get nodes by property list")
    print("="*80)
    if sample_nodes and test_property:
        values_list = [node.get(test_property) for node in sample_nodes[:3] if node.get(test_property)]
        if values_list:
            print(f"Searching for nodes with {test_property} in: {values_list}")
            results = get_nodes_by_property_list(db, test_collection, test_property, values_list)
    
    # Test 7: Cerca nodi per pattern (solo se proprietà è stringa)
    print("\n" + "="*80)
    print("TEST 7: Get nodes by pattern")
    print("="*80)
    if sample_nodes:
        # Trova una proprietà stringa
        string_prop = None
        string_value = None
        for key, value in sample_nodes[0].items():
            if isinstance(value, str) and not key.startswith('_') and len(value) > 2:
                string_prop = key
                string_value = value
                break
        
        if string_prop and string_value:
            # Cerca pattern con i primi 3 caratteri
            pattern = f"%{string_value[:3]}%"
            print(f"Searching for pattern '{pattern}' in property '{string_prop}'")
            results = get_nodes_by_pattern(db, test_collection, string_prop, pattern)
    
    # Test 8: Cerca nodi per range (solo se proprietà è numerica)
    print("\n" + "="*80)
    print("TEST 8: Get nodes by property range")
    print("="*80)
    if sample_nodes:
        # Trova una proprietà numerica
        numeric_prop = None
        numeric_value = None
        for key, value in sample_nodes[0].items():
            if isinstance(value, (int, float)) and not key.startswith('_'):
                numeric_prop = key
                numeric_value = value
                break
        
        if numeric_prop and numeric_value is not None:
            min_val = numeric_value - 10
            max_val = numeric_value + 10
            print(f"Searching for {numeric_prop} between {min_val} and {max_val}")
            results = get_nodes_by_property_range(db, test_collection, numeric_prop, min_val, max_val)
    
    # Test 9: Recupera nodi con proiezione
    print("\n" + "="*80)
    print("TEST 9: Get all nodes with property projection")
    print("="*80)
    if sample_nodes:
        # Prendi le prime 3 proprietà non di sistema
        props = [k for k in sample_nodes[0].keys() if not k.startswith('_')][:3]
        if props:
            print(f"Retrieving only properties: {props}")
            results = get_all_nodes_with_properties(db, test_collection, props)
            if results:
                print(f"First result: {results[0]}")
    
    print("\n" + "="*80)
    print("ALL TESTS COMPLETED!")
    print("="*80)

### ADDITIONAL TESTS (PKT-Transomics) ######

# retrieval function list
retrieval_functions_bundle = [
    get_nodes_by_property,
    get_nodes_by_multiple_properties,
    get_nodes_by_property_range,
    get_nodes_by_pattern,
    get_nodes_by_property_list,
    get_all_nodes_with_properties,
    count_nodes_by_property,
    get_distinct_property_values
]

#%%
if __name__ == "__main__":
    # Additional test: Get nodes by property 'class_code'='PR'
    import time
    test_collection = "nodes"
    test_property = "class_code"
    test_value = "PR"
    print("\n" + "="*80)
    print("ADDITIONAL TEST: Get nodes by property 'class_code'='PR'")
    print("="*80)
    if db.has_collection(test_collection):
        # time the function
        start_time = time.time()
        results = get_nodes_by_property(db, test_collection, test_property, test_value)
        end_time = time.time()
        print(f"Found {len(results)} nodes with {test_property}='{test_value}'")
        print(f"First 3 results:")
        for node in results[:3]:
            print("",node["_key"], node.get("label", "N/A")) 
        print(f"Time taken: {round(end_time - start_time, 2)} seconds")
#%%
if __name__ == "__main__":
    # test pattern matching
    pattern = "protein"
    string_prop = "label"
    test_collection = "nodes"
    print(f"Searching for pattern '{pattern}' in property '{string_prop}'")
    results = get_nodes_by_pattern(db, test_collection, string_prop, pattern)
    print(f"Found {len(results)} nodes with {string_prop}='{pattern}'")
    print(f"First 3 results:")
    for node in results[:3]:
        print("",node["_key"], node.get("label", "N/A")) 

#%%
# serche for _key contains "H7C1B8"
patterns_to_test = """Q96HL8
"""
patterns_to_test = patterns_to_test.strip().split("\n")

if __name__ == "__main__":
    for pattern in patterns_to_test:
        pattern = f"%{pattern}%"
        string_prop = "_key"
        test_collection = "nodes"
        print(f"Searching for pattern '{pattern}' in property '{string_prop}'")
        results = get_nodes_by_pattern(db, test_collection, string_prop, pattern, case_sensitive=True)
        print(f"Found {len(results)} nodes with {string_prop}='{pattern}'")
        print(f"Results:")
        for node in results:
            print("",node["_key"], node.get("label", "N/A"))

#%%
######## EXAMPLE DATA STRUCTURE ########


"""
The structure of the COLLECTIONs data example
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

