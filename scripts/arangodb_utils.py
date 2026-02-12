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

if __name__ == "__main__":
    db = setup_arangodb_connection()

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
def get_collections_data(db_connection, pattern=None, show_summary = True):
    """Check and print the number of documents in each collection
    
    Args:
        db_connection: ArangoDB database connection
        pattern: Optional regex pattern to filter collection names (e.g. "^node" or "edge.*")
    
    Returns:
        List of dicts with collection info: [{"name": str, "count": int, "type": str}, ...]
    """
    import re
    if db_connection is None: return []

    try:
        collections = db_connection.collections()
        
        print(f"\n📊 Collections in database:")
        total_docs = 0
        result_list = []
        
        for coll_info in collections:
            coll_name = coll_info['name']
            
            # Skip system collections
            if coll_name.startswith('_'):
                continue
            
            # Apply regex filter if provided
            if pattern:
                if not re.search(pattern, coll_name, re.IGNORECASE):
                    continue
            
            collection = db_connection.collection(coll_name)
            count = collection.count()
            coll_type = "edge" if coll_info['type'] == 3 else "document"
            
            print(f"  • {coll_name}: {count} documents ({coll_type})")
            total_docs += count
            result_list.append({"name": coll_name, "count": count, "type": coll_type})
        
        if show_summary:
            print(f"\n  Total: {len(result_list)} collections, {total_docs} documents")
            if pattern:
                print(f"  (filtered by pattern: '{pattern}')")

        collection_names = [coll['name'] for coll in result_list]
        
        return collection_names
            
    except Exception as e:
        print(f"✗ Error checking collections data: {e}")
        return []

def delete_collection(db_connection, collection_name, confirm=True):
    """Delete a collection completely from the database
    
    Args:
        db_connection: ArangoDB database connection
        collection_name: Name of the collection to delete
        confirm: If True, ask for confirmation before deleting
    
    Returns:
        True if deleted successfully, False otherwise
    """
    if db_connection is None:
        print("✗ No database connection available")
        return False
    
    try:
        if not db_connection.has_collection(collection_name):
            print(f"⚠ Collection '{collection_name}' does not exist.")
            return False
        
        # Get collection info before deletion
        collection = db_connection.collection(collection_name)
        doc_count = collection.count()
        
        if confirm:
            print(f"\n⚠ WARNING: You are about to delete collection '{collection_name}'")
            print(f"  This collection contains {doc_count} documents.")
            response = input("  Type 'yes' to confirm deletion: ")
            if response.lower() != 'yes':
                print("  Deletion cancelled.")
                return False
        
        # Delete the collection
        db_connection.delete_collection(collection_name)
        print(f"✓ Collection '{collection_name}' deleted successfully ({doc_count} documents removed)")
        return True
        
    except Exception as e:
        print(f"✗ Error deleting collection '{collection_name}': {e}")
        traceback.print_exc()
        return False

# funzioni per caricare multi omics datasets
if False:
    db_connection = setup_arangodb_connection()
    collections = get_collections_data(db_connection)
    collections_to_delete = [coll for coll in collections if coll.endswith(('_INDEX', '_SAMPLES'))]
    collections_to_delete
    # for coll in collections_to_delete:
        # # !!! ATTENZIONE: questa operazione è irreversibile, assicurati di avere un backup prima di eseguirla !!!
        # delete_collection(db_connection, coll, confirm=False)
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
def create_arango_graph(db_connection, graph_name='PKT_graph',
                        edge_collection='edges', 
                        from_vertex_collections=['nodes'],
                        to_vertex_collections=['nodes']):
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
                'edge_collection': edge_collection,
                'from_vertex_collections': from_vertex_collections,
                'to_vertex_collections': to_vertex_collections
            }
        ]
        
        # Create the graph
        print(f"🔨 Creating graph '{graph_name}'...")
        graph = db_connection.create_graph(
            name=graph_name,
            edge_definitions=edge_definitions
        )
        
        print(f"✅ Graph '{graph_name}' created successfully!")
        print(f"   - Vertex collections: {from_vertex_collections}")
        print(f"   - Edge collections: {edge_collection}")
        
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

def _resolve_center_node(db_connection, node_key=None, filters=None, node_collection='nodes'):
    """Helper to locate the center node either by explicit key or by property filters."""
    if node_key:
        nodes_coll = db_connection.collection(node_collection)
        return nodes_coll.get(node_key)

    if filters:
        filter_conditions = []
        bind_vars = {'@node_coll': node_collection}
        for idx, (prop, value) in enumerate(filters.items()):
            filter_conditions.append(f"doc.{prop} == @value{idx}")
            bind_vars[f"value{idx}"] = value
        condition_str = " AND ".join(filter_conditions) if filter_conditions else "true"
        aql = f"""
            FOR doc IN @@node_coll
                FILTER {condition_str}
                LIMIT 1
                RETURN doc
        """
        cursor = db_connection.aql.execute(aql, bind_vars=bind_vars)
        return next(cursor, None)

    return None


def get_node_centric_graph(db_connection, node_key=None, filters=None, edge_limit=500,
                           node_collection='nodes', edge_collection='edges'):
    """Return a node-centric subgraph consisting of the center node, its 1-hop neighbors, and connecting edges."""
    if db_connection is None:
        print("✗ No database connection available")
        return None

    if not node_key and not filters:
        print("✗ Provide either node_key or filters to identify the center node")
        return None

    center_node = _resolve_center_node(db_connection, node_key=node_key, filters=filters, node_collection=node_collection)
    if not center_node:
        print("✗ Center node not found with provided parameters")
        return None

    node_key = center_node['_key']
    start_vertex = f"{node_collection}/{node_key}"

    edge_query = """
        FOR edge IN @@edge_coll
            FILTER edge._from == @start_vertex OR edge._to == @start_vertex
            LIMIT @edge_limit
            RETURN edge
    """
    edges_cursor = db_connection.aql.execute(edge_query, bind_vars={
        '@edge_coll': edge_collection,
        'start_vertex': start_vertex, 
        'edge_limit': edge_limit
    })
    edges_list = list(edges_cursor)

    neighbor_keys = set()
    for edge in edges_list:
        neighbor_keys.add(edge['_from'].split('/')[1])
        neighbor_keys.add(edge['_to'].split('/')[1])
    neighbor_keys.discard(node_key)

    neighbors = []
    if neighbor_keys:
        neighbors_query = """
            FOR doc IN @@node_coll
                FILTER doc._key IN @keys
                RETURN doc
        """
        neighbors_cursor = db_connection.aql.execute(neighbors_query, bind_vars={
            '@node_coll': node_collection,
            'keys': list(neighbor_keys)
        })
        neighbors = list(neighbors_cursor)

    result = {
        'center_node': center_node,
        'neighbors': neighbors,
        'edges': edges_list
    }

    print(f"✓ Found center node {node_key} with {len(neighbors)} neighbors and {len(edges_list)} edges")
    return result


def traverse_node_centric_graph(db_connection, node_key=None, filters=None, depth=2, direction='ANY', limit=1000, graph_name='PKT_graph'):
    """Perform a multi-hop traversal from a center node using the saved graph and return the induced subgraph."""
    if db_connection is None:
        print("✗ No database connection available")
        return None

    if not node_key and not filters:
        print("✗ Provide either node_key or filters to identify the center node")
        return None

    direction_token = direction.upper()
    if direction_token not in {'ANY', 'INBOUND', 'OUTBOUND'}:
        print("✗ direction must be one of: 'ANY', 'INBOUND', 'OUTBOUND'")
        return None

    if depth < 1:
        print("✗ depth must be >= 1")
        return None

    if not graph_name.replace('_', '').isalnum():
        print("✗ graph_name contains invalid characters")
        return None

    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist. Create it with create_arango_graph() first.")
        return None

    center_node = _resolve_center_node(db_connection, node_key=node_key, filters=filters)
    if not center_node:
        print("✗ Center node not found with provided parameters")
        return None

    node_key = center_node['_key']
    start_vertex = f"nodes/{node_key}"

    traversal_query = f"""
        FOR v, e IN 1..{depth} {direction_token} @start_vertex GRAPH '{graph_name}'
            OPTIONS {{uniqueVertices: 'global', bfs: true}}
            LIMIT @limit
            RETURN {{vertex: v, edge: e}}
    """

    cursor = db_connection.aql.execute(traversal_query, bind_vars={'start_vertex': start_vertex, 'limit': limit})

    nodes_dict = {node_key: center_node}
    edges_list = []

    for item in cursor:
        vertex = item.get('vertex')
        edge = item.get('edge')
        if vertex and vertex['_key'] not in nodes_dict:
            nodes_dict[vertex['_key']] = vertex
        if edge:
            edges_list.append(edge)

    neighbors = [node for key, node in nodes_dict.items() if key != node_key]

    result = {
        'center_node': center_node,
        'nodes': list(nodes_dict.values()),
        'neighbors': neighbors,
        'edges': edges_list,
        'depth': depth
    }

    print(f"✓ Traversed depth {depth} from {node_key}: {len(nodes_dict)} nodes, {len(edges_list)} edges (limit {limit})")
    return result


# --- GRAPH QUERY FUNCTIONS ---

def find_shortest_path(db_connection, start_key, end_key, direction='ANY',
                       weight_attribute=None, graph_name='PKT_graph'):
    """
    Trova il percorso più breve tra due nodi.
    
    Args:
        db_connection: Connessione ArangoDB
        start_key: _key del nodo di partenza
        end_key: _key del nodo di arrivo
        direction: 'ANY', 'OUTBOUND', 'INBOUND'
        weight_attribute: Nome attributo peso per weighted shortest path (None = unweighted)
        graph_name: Nome del grafo
    
    Returns:
        Dict con 'vertices', 'edges', 'distance' o None se non trovato
    """
    if db_connection is None:
        print("✗ No database connection available")
        return None
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return None
    
    direction = direction.upper()
    if direction not in {'ANY', 'OUTBOUND', 'INBOUND'}:
        print("✗ direction must be 'ANY', 'OUTBOUND', or 'INBOUND'")
        return None
    
    start_vertex = f"nodes/{start_key}"
    end_vertex = f"nodes/{end_key}"
    
    try:
        if weight_attribute:
            query = f"""
                FOR v, e IN {direction} SHORTEST_PATH @start TO @end GRAPH '{graph_name}'
                    OPTIONS {{weightAttribute: @weight}}
                    RETURN {{vertex: v, edge: e}}
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'end': end_vertex, 'weight': weight_attribute
            })
        else:
            query = f"""
                FOR v, e IN {direction} SHORTEST_PATH @start TO @end GRAPH '{graph_name}'
                    RETURN {{vertex: v, edge: e}}
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'end': end_vertex
            })
        
        results = list(cursor)
        
        if not results:
            print(f"✗ No path found between {start_key} and {end_key}")
            return None
        
        vertices = [r['vertex'] for r in results if r['vertex']]
        edges = [r['edge'] for r in results if r['edge']]
        
        result = {
            'vertices': vertices,
            'edges': edges,
            'distance': len(edges),
            'start': start_key,
            'end': end_key
        }
        
        print(f"✓ Found path: {start_key} → {end_key} (distance: {len(edges)} hops)")
        return result
        
    except Exception as e:
        print(f"✗ Error finding shortest path: {e}")
        traceback.print_exc()
        return None


def find_k_shortest_paths(db_connection, start_key, end_key, k=5, direction='OUTBOUND',
                          weight_attribute=None, graph_name='PKT_graph'):
    """
    Trova i k percorsi più brevi tra due nodi.
    
    Args:
        db_connection: Connessione ArangoDB
        start_key: _key del nodo di partenza
        end_key: _key del nodo di arrivo
        k: Numero di percorsi da restituire
        direction: 'OUTBOUND' o 'INBOUND' (K_SHORTEST_PATHS non supporta ANY)
        weight_attribute: Nome attributo peso (None = unweighted)
        graph_name: Nome del grafo
    
    Returns:
        Lista di path, ognuno con 'vertices', 'edges', 'weight'
    """
    if db_connection is None:
        print("✗ No database connection available")
        return []
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return []
    
    direction = direction.upper()
    if direction not in {'OUTBOUND', 'INBOUND'}:
        print("✗ direction must be 'OUTBOUND' or 'INBOUND' for K_SHORTEST_PATHS")
        return []
    
    start_vertex = f"nodes/{start_key}"
    end_vertex = f"nodes/{end_key}"
    
    try:
        if weight_attribute:
            query = f"""
                FOR path IN {direction} K_SHORTEST_PATHS @start TO @end GRAPH '{graph_name}'
                    OPTIONS {{weightAttribute: @weight}}
                    LIMIT @k
                    RETURN path
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'end': end_vertex, 'weight': weight_attribute, 'k': k
            })
        else:
            query = f"""
                FOR path IN {direction} K_SHORTEST_PATHS @start TO @end GRAPH '{graph_name}'
                    LIMIT @k
                    RETURN path
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'end': end_vertex, 'k': k
            })
        
        paths = list(cursor)
        
        results = []
        for i, path in enumerate(paths):
            results.append({
                'rank': i + 1,
                'vertices': path.get('vertices', []),
                'edges': path.get('edges', []),
                'weight': path.get('weight', len(path.get('edges', [])))
            })
        
        print(f"✓ Found {len(results)} paths between {start_key} and {end_key}")
        return results
        
    except Exception as e:
        print(f"✗ Error finding k shortest paths: {e}")
        traceback.print_exc()
        return []


def find_all_shortest_paths(db_connection, start_key, end_key, direction='OUTBOUND',
                            graph_name='PKT_graph'):
    """
    Trova tutti i percorsi minimi (stessa lunghezza) tra due nodi.
    
    Args:
        db_connection: Connessione ArangoDB
        start_key: _key del nodo di partenza
        end_key: _key del nodo di arrivo
        direction: 'OUTBOUND' o 'INBOUND'
        graph_name: Nome del grafo
    
    Returns:
        Lista di path con 'vertices' e 'edges'
    """
    if db_connection is None:
        print("✗ No database connection available")
        return []
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return []
    
    direction = direction.upper()
    if direction not in {'OUTBOUND', 'INBOUND'}:
        print("✗ direction must be 'OUTBOUND' or 'INBOUND' for ALL_SHORTEST_PATHS")
        return []
    
    start_vertex = f"nodes/{start_key}"
    end_vertex = f"nodes/{end_key}"
    
    try:
        query = f"""
            FOR path IN {direction} ALL_SHORTEST_PATHS @start TO @end GRAPH '{graph_name}'
                RETURN path
        """
        cursor = db_connection.aql.execute(query, bind_vars={
            'start': start_vertex, 'end': end_vertex
        })
        
        paths = list(cursor)
        
        results = []
        for path in paths:
            results.append({
                'vertices': path.get('vertices', []),
                'edges': path.get('edges', [])
            })
        
        if results:
            print(f"✓ Found {len(results)} shortest paths (all same length: {len(results[0]['edges'])} hops)")
        else:
            print(f"✗ No paths found between {start_key} and {end_key}")
        
        return results
        
    except Exception as e:
        print(f"✗ Error finding all shortest paths: {e}")
        traceback.print_exc()
        return []


def find_pattern(db_connection, pattern_types, edge_filters=None, direction='OUTBOUND',
                 limit=100, graph_name='PKT_graph'):
    """
    Trova motivi (pattern) nel grafo, es. gene → protein → disease.
    
    Args:
        db_connection: Connessione ArangoDB
        pattern_types: Lista di bioentity_type da matchare in sequenza
                       Es: ['gene', 'protein', 'disease']
        edge_filters: Dict opzionale con filtri per gli edge
                      Es: {'predicate_label': 'regulates'}
        direction: Direzione traversal ('OUTBOUND', 'INBOUND', 'ANY')
        limit: Numero massimo di pattern da restituire
        graph_name: Nome del grafo
    
    Returns:
        Lista di pattern trovati, ognuno con i nodi matchati
    """
    if db_connection is None:
        print("✗ No database connection available")
        return []
    
    if not pattern_types or len(pattern_types) < 2:
        print("✗ pattern_types must have at least 2 elements")
        return []
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return []
    
    direction = direction.upper()
    if direction not in {'OUTBOUND', 'INBOUND', 'ANY'}:
        print("✗ direction must be 'OUTBOUND', 'INBOUND', or 'ANY'")
        return []
    
    try:
        # Costruisci la query dinamicamente
        query_parts = []
        return_fields = []
        
        # Primo nodo
        query_parts.append(f"FOR v0 IN nodes")
        query_parts.append(f"    FILTER v0.bioentity_type == @type0")
        return_fields.append("v0")
        
        # Nodi successivi con traversal
        for i in range(1, len(pattern_types)):
            edge_var = f"e{i-1}"
            node_var = f"v{i}"
            
            query_parts.append(f"    FOR {node_var}, {edge_var} IN 1..1 {direction} v{i-1} GRAPH '{graph_name}'")
            query_parts.append(f"        FILTER {node_var}.bioentity_type == @type{i}")
            
            # Aggiungi filtri sugli edge se specificati
            if edge_filters:
                for key, value in edge_filters.items():
                    query_parts.append(f"        FILTER {edge_var}.{key} == @edge_{key}")
            
            return_fields.append(node_var)
            return_fields.append(edge_var)
        
        # Limit e return
        query_parts.append(f"    LIMIT @limit")
        return_obj = "{" + ", ".join([f"{f}: {f}" for f in return_fields]) + "}"
        query_parts.append(f"    RETURN {return_obj}")
        
        query = "\n".join(query_parts)
        
        # Bind vars
        bind_vars = {'limit': limit}
        for i, ptype in enumerate(pattern_types):
            bind_vars[f'type{i}'] = ptype
        if edge_filters:
            for key, value in edge_filters.items():
                bind_vars[f'edge_{key}'] = value
        
        cursor = db_connection.aql.execute(query, bind_vars=bind_vars)
        results = list(cursor)
        
        # Formatta i risultati
        patterns = []
        for r in results:
            pattern = {'nodes': [], 'edges': []}
            for i in range(len(pattern_types)):
                node = r.get(f'v{i}')
                if node:
                    pattern['nodes'].append(node)
                if i > 0:
                    edge = r.get(f'e{i-1}')
                    if edge:
                        pattern['edges'].append(edge)
            patterns.append(pattern)
        
        pattern_str = " → ".join(pattern_types)
        print(f"✓ Found {len(patterns)} patterns matching: {pattern_str}")
        return patterns
        
    except Exception as e:
        print(f"✗ Error finding patterns: {e}")
        traceback.print_exc()
        return []


def find_common_neighbors(db_connection, node_key_a, node_key_b, direction='ANY',
                          filter_type=None, graph_name='PKT_graph'):
    """
    Trova i vicini comuni tra due nodi.
    
    Args:
        db_connection: Connessione ArangoDB
        node_key_a: _key del primo nodo
        node_key_b: _key del secondo nodo
        direction: Direzione ('ANY', 'OUTBOUND', 'INBOUND')
        filter_type: Filtra per bioentity_type (es. 'protein')
        graph_name: Nome del grafo
    
    Returns:
        Dict con 'common_neighbors', 'only_a', 'only_b', 'jaccard_similarity'
    """
    if db_connection is None:
        print("✗ No database connection available")
        return None
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return None
    
    direction = direction.upper()
    if direction not in {'ANY', 'OUTBOUND', 'INBOUND'}:
        print("✗ direction must be 'ANY', 'OUTBOUND', or 'INBOUND'")
        return None
    
    vertex_a = f"nodes/{node_key_a}"
    vertex_b = f"nodes/{node_key_b}"
    
    try:
        if filter_type:
            query = f"""
                LET neighbors_a = (
                    FOR v IN 1..1 {direction} @vertex_a GRAPH '{graph_name}'
                        FILTER v.bioentity_type == @filter_type
                        RETURN v
                )
                LET neighbors_b = (
                    FOR v IN 1..1 {direction} @vertex_b GRAPH '{graph_name}'
                        FILTER v.bioentity_type == @filter_type
                        RETURN v
                )
                LET keys_a = neighbors_a[*]._key
                LET keys_b = neighbors_b[*]._key
                LET common_keys = INTERSECTION(keys_a, keys_b)
                LET only_a_keys = MINUS(keys_a, keys_b)
                LET only_b_keys = MINUS(keys_b, keys_a)
                RETURN {{
                    common: (FOR n IN neighbors_a FILTER n._key IN common_keys RETURN n),
                    only_a: (FOR n IN neighbors_a FILTER n._key IN only_a_keys RETURN n),
                    only_b: (FOR n IN neighbors_b FILTER n._key IN only_b_keys RETURN n),
                    count_a: LENGTH(keys_a),
                    count_b: LENGTH(keys_b),
                    count_common: LENGTH(common_keys)
                }}
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'vertex_a': vertex_a, 'vertex_b': vertex_b, 'filter_type': filter_type
            })
        else:
            query = f"""
                LET neighbors_a = (
                    FOR v IN 1..1 {direction} @vertex_a GRAPH '{graph_name}'
                        RETURN v
                )
                LET neighbors_b = (
                    FOR v IN 1..1 {direction} @vertex_b GRAPH '{graph_name}'
                        RETURN v
                )
                LET keys_a = neighbors_a[*]._key
                LET keys_b = neighbors_b[*]._key
                LET common_keys = INTERSECTION(keys_a, keys_b)
                LET only_a_keys = MINUS(keys_a, keys_b)
                LET only_b_keys = MINUS(keys_b, keys_a)
                RETURN {{
                    common: (FOR n IN neighbors_a FILTER n._key IN common_keys RETURN n),
                    only_a: (FOR n IN neighbors_a FILTER n._key IN only_a_keys RETURN n),
                    only_b: (FOR n IN neighbors_b FILTER n._key IN only_b_keys RETURN n),
                    count_a: LENGTH(keys_a),
                    count_b: LENGTH(keys_b),
                    count_common: LENGTH(common_keys)
                }}
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'vertex_a': vertex_a, 'vertex_b': vertex_b
            })
        
        result = next(cursor, None)
        
        if not result:
            print(f"✗ Could not compute common neighbors")
            return None
        
        # Calcola Jaccard similarity
        union_size = result['count_a'] + result['count_b'] - result['count_common']
        jaccard = result['count_common'] / union_size if union_size > 0 else 0.0
        
        output = {
            'node_a': node_key_a,
            'node_b': node_key_b,
            'common_neighbors': result['common'],
            'only_a_neighbors': result['only_a'],
            'only_b_neighbors': result['only_b'],
            'count_common': result['count_common'],
            'count_a': result['count_a'],
            'count_b': result['count_b'],
            'jaccard_similarity': round(jaccard, 4)
        }
        
        print(f"✓ {node_key_a} and {node_key_b}: {result['count_common']} common neighbors (Jaccard: {jaccard:.3f})")
        return output
        
    except Exception as e:
        print(f"✗ Error finding common neighbors: {e}")
        traceback.print_exc()
        return None


def get_node_neighbors(db_connection, node_key, direction='ANY', depth=1,
                       filter_type=None, limit=500, graph_name='PKT_graph'):
    """
    Ottiene i vicini di un nodo a una certa profondità.
    
    Args:
        db_connection: Connessione ArangoDB
        node_key: _key del nodo centrale
        direction: 'ANY', 'OUTBOUND', 'INBOUND'
        depth: Profondità del vicinato (1 = vicini diretti)
        filter_type: Filtra per bioentity_type
        limit: Numero massimo di vicini
        graph_name: Nome del grafo
    
    Returns:
        Lista di nodi vicini
    """
    if db_connection is None:
        print("✗ No database connection available")
        return []
    
    if not db_connection.has_graph(graph_name):
        print(f"✗ Graph '{graph_name}' does not exist")
        return []
    
    direction = direction.upper()
    start_vertex = f"nodes/{node_key}"
    
    try:
        if filter_type:
            query = f"""
                FOR v IN 1..@depth {direction} @start GRAPH '{graph_name}'
                    FILTER v.bioentity_type == @filter_type
                    LIMIT @limit
                    RETURN DISTINCT v
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'depth': depth, 'filter_type': filter_type, 'limit': limit
            })
        else:
            query = f"""
                FOR v IN 1..@depth {direction} @start GRAPH '{graph_name}'
                    LIMIT @limit
                    RETURN DISTINCT v
            """
            cursor = db_connection.aql.execute(query, bind_vars={
                'start': start_vertex, 'depth': depth, 'limit': limit
            })
        
        neighbors = list(cursor)
        
        type_info = f" (type: {filter_type})" if filter_type else ""
        print(f"✓ Found {len(neighbors)} neighbors of {node_key} at depth {depth}{type_info}")
        return neighbors
        
    except Exception as e:
        print(f"✗ Error getting neighbors: {e}")
        traceback.print_exc()
        return []


def plot_subgraph(subgraph, layout='spring', node_size=400, figsize=(14, 10), 
                  show_labels=True, highlight_center=True, save_path=None):
    """
    Visualizza un sottografo restituito da get_node_centric_graph o traverse_node_centric_graph.
    
    Args:
        subgraph: Dizionario con 'center_node', 'neighbors'/'nodes', 'edges'
        layout: Algoritmo di layout ('spring', 'circular', 'kamada_kawai', 'shell')
        node_size: Dimensione dei nodi
        figsize: Dimensione della figura (width, height)
        show_labels: Se True, mostra le etichette sui nodi
        highlight_center: Se True, evidenzia il nodo centrale
        save_path: Percorso per salvare l'immagine (None = non salva)
    
    Returns:
        NetworkX DiGraph object
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        from collections import defaultdict
        
        if subgraph is None:
            print("✗ Subgraph is None")
            return None
        
        center_node = subgraph.get('center_node')
        edges_list = subgraph.get('edges', [])
        
        # Supporta sia 'nodes' (da traverse) che 'neighbors' (da get_node_centric)
        if 'nodes' in subgraph:
            all_nodes = subgraph['nodes']
        else:
            all_nodes = [center_node] + subgraph.get('neighbors', [])
        
        if not all_nodes:
            print("✗ No nodes in subgraph")
            return None
        
        # Crea grafo NetworkX
        G = nx.DiGraph()
        
        # Mappa nodi per lookup veloce
        nodes_map = {node['_key']: node for node in all_nodes}
        
        # Aggiungi nodi con attributi
        bioentity_types = defaultdict(list)
        for node in all_nodes:
            node_key = node['_key']
            attrs = {k: v for k, v in node.items() if not k.startswith('_')}
            G.add_node(node_key, **attrs)
            bioentity_types[node.get('bioentity_type', 'unknown')].append(node_key)
        
        # Aggiungi archi
        edge_labels_map = defaultdict(int)
        for edge in edges_list:
            source = edge['_from'].split('/')[1]
            target = edge['_to'].split('/')[1]
            if source in nodes_map and target in nodes_map:
                attrs = {k: v for k, v in edge.items() if not k.startswith('_')}
                G.add_edge(source, target, **attrs)
                edge_labels_map[edge.get('predicate_label', 'related')] += 1
        
        # Layout
        if layout == 'spring':
            pos = nx.spring_layout(G, k=1.5, iterations=50, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            pos = nx.kamada_kawai_layout(G)
        elif layout == 'shell':
            # Metti il centro al centro, i vicini attorno
            center_key = center_node['_key'] if center_node else None
            shells = [[center_key]] if center_key else []
            other_nodes = [n for n in G.nodes() if n != center_key]
            if other_nodes:
                shells.append(other_nodes)
            pos = nx.shell_layout(G, nlist=shells if len(shells) > 1 else None)
        else:
            pos = nx.random_layout(G)
        
        # Colori per bioentity type
        color_map = {
            'protein': '#FF6B6B',
            'gene': '#4ECDC4',
            'rna': '#AA96DA',
            'chemical': '#95E1D3',
            'disease': '#F38181',
            'pathway': '#FEE191',
            'phenotype': '#FFB6B9',
            'metabolite': '#FCBAD3',
            'go': '#A8D8EA',
            'unknown': '#CCCCCC'
        }
        
        node_colors = []
        for node_key in G.nodes():
            node_data = nodes_map.get(node_key, {})
            btype = node_data.get('bioentity_type', 'unknown')
            node_colors.append(color_map.get(btype, '#CCCCCC'))
        
        # Plot
        plt.figure(figsize=figsize)
        
        # Disegna tutti i nodi
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                               node_size=node_size, alpha=0.85)
        
        # Evidenzia nodo centrale
        if highlight_center and center_node:
            center_key = center_node['_key']
            if center_key in G.nodes():
                nx.draw_networkx_nodes(G, pos, nodelist=[center_key],
                                       node_color='gold', node_size=node_size * 1.8,
                                       edgecolors='black', linewidths=2)
        
        # Disegna archi
        nx.draw_networkx_edges(G, pos, edge_color='gray', 
                               arrows=True, arrowsize=12, 
                               alpha=0.5, width=0.8,
                               connectionstyle="arc3,rad=0.1")
        
        # Etichette
        if show_labels and len(G.nodes()) <= 50:
            labels = {}
            for node_key in G.nodes():
                node_data = nodes_map.get(node_key, {})
                label = node_data.get('label', node_key)
                labels[node_key] = label[:20] + '...' if len(str(label)) > 20 else label
            nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
        
        # Legenda
        legend_elements = []
        from matplotlib.patches import Patch
        for btype, nodes in sorted(bioentity_types.items(), key=lambda x: -len(x[1])):
            color = color_map.get(btype, '#CCCCCC')
            legend_elements.append(Patch(facecolor=color, label=f'{btype} ({len(nodes)})'))
        
        if legend_elements:
            plt.legend(handles=legend_elements, loc='upper left', fontsize=8)
        
        # Titolo
        center_label = center_node.get('label', center_node.get('_key', '?')) if center_node else '?'
        depth_info = f", depth={subgraph.get('depth', 1)}" if 'depth' in subgraph else ''
        plt.title(f"Subgraph: {center_label}\n({len(G.nodes())} nodes, {len(G.edges())} edges{depth_info})",
                  fontsize=14, fontweight='bold')
        
        plt.axis('off')
        plt.tight_layout()
        
        # Salva se richiesto
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved plot to: {save_path}")
        
        plt.show()
        
        # Stampa statistiche
        print(f"\n📈 Subgraph statistics:")
        print(f"   Nodes by bioentity_type:")
        for btype, nodes in sorted(bioentity_types.items(), key=lambda x: -len(x[1])):
            print(f"     • {btype}: {len(nodes)}")
        print(f"   Edges by predicate_label (top 5):")
        for label, count in sorted(edge_labels_map.items(), key=lambda x: -x[1])[:5]:
            print(f"     • {label}: {count}")
        
        return G
        
    except ImportError:
        print("✗ Missing libraries. Install with: pip install networkx matplotlib")
        return None
    except Exception as e:
        print(f"✗ Error plotting subgraph: {e}")
        traceback.print_exc()
        return None


def plot_graph_result(result, layout='kamada_kawai', node_size=400, figsize=(14, 10),
                      show_labels=True, highlight_endpoints=True, save_path=None,
                      title=None, path_index=0):
    """
    Funzione universale per plottare output di tutte le funzioni graph query.
    
    Rileva automaticamente il tipo di risultato:
    - Shortest path: dict con 'vertices' e 'edges'
    - K shortest paths / All shortest paths: lista di path
    - Pattern matching: lista di pattern con 'nodes' e 'edges'
    - Common neighbors: dict con 'common_neighbors', 'only_a_neighbors', etc.
    - Node neighbors: lista di nodi
    - Subgraph (traversal): dict con 'center_node', 'neighbors'/'nodes', 'edges'
    
    Args:
        result: Output di una delle funzioni graph query
        layout: 'spring', 'circular', 'kamada_kawai', 'shell', 'path'
        node_size: Dimensione nodi
        figsize: Dimensione figura
        show_labels: Mostra etichette
        highlight_endpoints: Evidenzia nodi start/end per path
        save_path: Percorso per salvare immagine
        title: Titolo custom (None = auto)
        path_index: Per k_shortest_paths, quale path plottare (0 = primo)
    
    Returns:
        NetworkX DiGraph o None
    
    Esempi:
        # Shortest path
        path = find_shortest_path(db, 'A', 'B')
        plot_graph_result(path)
        
        # K shortest paths (plotta il terzo path)
        paths = find_k_shortest_paths(db, 'A', 'B', k=5)
        plot_graph_result(paths, path_index=2)
        
        # Pattern matching
        patterns = find_pattern(db, ['gene', 'protein', 'disease'])
        plot_graph_result(patterns)  # plotta tutti i pattern sovrapposti
        
        # Common neighbors
        common = find_common_neighbors(db, 'A', 'B')
        plot_graph_result(common)
        
        # Node neighbors
        neighbors = get_node_neighbors(db, 'A', depth=2)
        plot_graph_result(neighbors)
        
        # Traversal subgraph
        subgraph = traverse_node_centric_graph(db, node_key='A', depth=2)
        plot_graph_result(subgraph)
    """
    try:
        import networkx as nx
        import matplotlib.pyplot as plt
        from collections import defaultdict
        from matplotlib.patches import Patch
        
        if result is None:
            print("✗ Result is None")
            return None
        
        # --- DETECT RESULT TYPE ---
        result_type = None
        nodes_list = []
        edges_list = []
        special_nodes = {}  # per evidenziare nodi speciali
        auto_title = "Graph Result"
        
        # Caso 1: Shortest path singolo
        if isinstance(result, dict) and 'vertices' in result and 'start' in result:
            result_type = 'shortest_path'
            nodes_list = result.get('vertices', [])
            edges_list = result.get('edges', [])
            special_nodes['start'] = result.get('start')
            special_nodes['end'] = result.get('end')
            auto_title = f"Shortest Path: {result.get('start')} → {result.get('end')} ({result.get('distance', '?')} hops)"
        
        # Caso 2: K shortest paths o All shortest paths (lista di path)
        elif isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) and 'vertices' in result[0]:
            result_type = 'k_paths'
            if path_index >= len(result):
                print(f"✗ path_index {path_index} out of range (max {len(result)-1})")
                path_index = 0
            selected_path = result[path_index]
            nodes_list = selected_path.get('vertices', [])
            edges_list = selected_path.get('edges', [])
            if nodes_list:
                special_nodes['start'] = nodes_list[0].get('_key') if isinstance(nodes_list[0], dict) else None
                special_nodes['end'] = nodes_list[-1].get('_key') if isinstance(nodes_list[-1], dict) else None
            rank = selected_path.get('rank', path_index + 1)
            weight = selected_path.get('weight', len(edges_list))
            auto_title = f"Path {rank}/{len(result)} (weight: {weight})"
        
        # Caso 3: Pattern matching (lista di pattern con 'nodes')
        elif isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) and 'nodes' in result[0]:
            result_type = 'patterns'
            # Combina tutti i pattern in un grafo
            all_nodes = {}
            all_edges = []
            for pattern in result:
                for node in pattern.get('nodes', []):
                    all_nodes[node['_key']] = node
                all_edges.extend(pattern.get('edges', []))
            nodes_list = list(all_nodes.values())
            edges_list = all_edges
            auto_title = f"Pattern Matching ({len(result)} patterns found)"
        
        # Caso 4: Common neighbors
        elif isinstance(result, dict) and 'common_neighbors' in result:
            result_type = 'common_neighbors'
            # Crea nodi fittizi per A e B
            fake_node_a = {'_key': result['node_a'], 'bioentity_type': 'query', 'label': result['node_a']}
            fake_node_b = {'_key': result['node_b'], 'bioentity_type': 'query', 'label': result['node_b']}
            nodes_list = [fake_node_a, fake_node_b] + result.get('common_neighbors', [])
            # Aggiungi anche solo_a e solo_b se presenti
            nodes_list.extend(result.get('only_a_neighbors', []))
            nodes_list.extend(result.get('only_b_neighbors', []))
            # Non abbiamo edges reali, ma possiamo creare connessioni concettuali
            edges_list = []
            special_nodes['start'] = result['node_a']
            special_nodes['end'] = result['node_b']
            jaccard = result.get('jaccard_similarity', 0)
            auto_title = f"Common Neighbors: {result['node_a']} ∩ {result['node_b']}\n(Jaccard: {jaccard:.3f}, Common: {result.get('count_common', 0)})"
        
        # Caso 5: Lista di nodi (get_node_neighbors)
        elif isinstance(result, list) and len(result) > 0 and isinstance(result[0], dict) and '_key' in result[0] and 'vertices' not in result[0]:
            result_type = 'neighbors_list'
            nodes_list = result
            edges_list = []
            auto_title = f"Node Neighbors ({len(result)} nodes)"
        
        # Caso 6: Subgraph da traversal (ha center_node)
        elif isinstance(result, dict) and 'center_node' in result:
            # Usa la funzione esistente plot_subgraph
            return plot_subgraph(result, layout=layout, node_size=node_size, figsize=figsize,
                                 show_labels=show_labels, highlight_center=highlight_endpoints,
                                 save_path=save_path)
        
        else:
            print("✗ Cannot detect result type. Supported: shortest_path, k_paths, patterns, common_neighbors, neighbors_list, subgraph")
            return None
        
        if not nodes_list:
            print("✗ No nodes to plot")
            return None
        
        # --- BUILD NETWORKX GRAPH ---
        G = nx.DiGraph()
        
        # Mappa nodi
        nodes_map = {}
        for node in nodes_list:
            if isinstance(node, dict) and '_key' in node:
                nodes_map[node['_key']] = node
        
        # Aggiungi nodi
        bioentity_types = defaultdict(list)
        for node_key, node in nodes_map.items():
            attrs = {k: v for k, v in node.items() if not k.startswith('_')}
            G.add_node(node_key, **attrs)
            bioentity_types[node.get('bioentity_type', 'unknown')].append(node_key)
        
        # Aggiungi archi
        edge_labels_map = defaultdict(int)
        for edge in edges_list:
            if isinstance(edge, dict) and '_from' in edge and '_to' in edge:
                source = edge['_from'].split('/')[1]
                target = edge['_to'].split('/')[1]
                if source in nodes_map and target in nodes_map:
                    attrs = {k: v for k, v in edge.items() if not k.startswith('_')}
                    G.add_edge(source, target, **attrs)
                    edge_labels_map[edge.get('predicate_label', 'related')] += 1
        
        # --- LAYOUT ---
        if layout == 'path' and result_type in ('shortest_path', 'k_paths'):
            # Layout lineare per path
            pos = {}
            for i, node_key in enumerate(G.nodes()):
                pos[node_key] = (i, 0)
        elif layout == 'spring':
            pos = nx.spring_layout(G, k=1.5, iterations=50, seed=42)
        elif layout == 'circular':
            pos = nx.circular_layout(G)
        elif layout == 'kamada_kawai':
            try:
                pos = nx.kamada_kawai_layout(G)
            except:
                pos = nx.spring_layout(G, seed=42)
        elif layout == 'shell':
            pos = nx.shell_layout(G)
        else:
            pos = nx.spring_layout(G, seed=42)
        
        # --- COLORI ---
        color_map = {
            'protein': '#FF6B6B',
            'gene': '#4ECDC4',
            'rna': '#AA96DA',
            'chemical': '#95E1D3',
            'disease': '#F38181',
            'pathway': '#FEE191',
            'phenotype': '#FFB6B9',
            'metabolite': '#FCBAD3',
            'go': '#A8D8EA',
            'query': '#FFD700',
            'unknown': '#CCCCCC'
        }
        
        node_colors = []
        for node_key in G.nodes():
            node_data = nodes_map.get(node_key, {})
            btype = node_data.get('bioentity_type', 'unknown')
            node_colors.append(color_map.get(btype, '#CCCCCC'))
        
        # --- PLOT ---
        plt.figure(figsize=figsize)
        
        # Nodi base
        nx.draw_networkx_nodes(G, pos, node_color=node_colors, 
                               node_size=node_size, alpha=0.85)
        
        # Evidenzia endpoint
        if highlight_endpoints:
            highlight_list = []
            for key in ['start', 'end']:
                node_key = special_nodes.get(key)
                if node_key and node_key in G.nodes():
                    highlight_list.append(node_key)
            if highlight_list:
                colors = ['#00FF00', '#FF0000'][:len(highlight_list)]  # verde=start, rosso=end
                for i, node_key in enumerate(highlight_list):
                    nx.draw_networkx_nodes(G, pos, nodelist=[node_key],
                                           node_color=colors[i], node_size=node_size * 1.5,
                                           edgecolors='black', linewidths=2)
        
        # Archi
        if G.edges():
            nx.draw_networkx_edges(G, pos, edge_color='gray', 
                                   arrows=True, arrowsize=12, 
                                   alpha=0.6, width=1.2,
                                   connectionstyle="arc3,rad=0.1")
        
        # Etichette
        if show_labels and len(G.nodes()) <= 60:
            labels = {}
            for node_key in G.nodes():
                node_data = nodes_map.get(node_key, {})
                label = node_data.get('label', node_key)
                labels[node_key] = str(label)[:18] + '...' if len(str(label)) > 18 else str(label)
            nx.draw_networkx_labels(G, pos, labels, font_size=8, font_weight='bold')
        
        # Legenda
        legend_elements = []
        for btype, nodes in sorted(bioentity_types.items(), key=lambda x: -len(x[1])):
            color = color_map.get(btype, '#CCCCCC')
            legend_elements.append(Patch(facecolor=color, label=f'{btype} ({len(nodes)})'))
        
        if highlight_endpoints and special_nodes.get('start'):
            legend_elements.append(Patch(facecolor='#00FF00', label='Start'))
        if highlight_endpoints and special_nodes.get('end'):
            legend_elements.append(Patch(facecolor='#FF0000', label='End'))
        
        if legend_elements:
            plt.legend(handles=legend_elements, loc='upper left', fontsize=8)
        
        # Titolo
        final_title = title if title else auto_title
        plt.title(f"{final_title}\n({len(G.nodes())} nodes, {len(G.edges())} edges)",
                  fontsize=13, fontweight='bold')
        
        plt.axis('off')
        plt.tight_layout()
        
        if save_path:
            plt.savefig(save_path, dpi=300, bbox_inches='tight')
            print(f"✓ Saved plot to: {save_path}")
        
        plt.show()
        
        # Stats
        print(f"\n📈 Result statistics:")
        print(f"   Type: {result_type}")
        print(f"   Nodes by bioentity_type:")
        for btype, nodes in sorted(bioentity_types.items(), key=lambda x: -len(x[1])):
            print(f"     • {btype}: {len(nodes)}")
        if edge_labels_map:
            print(f"   Edges by predicate_label:")
            for label, count in sorted(edge_labels_map.items(), key=lambda x: -x[1])[:5]:
                print(f"     • {label}: {count}")
        
        return G
        
    except ImportError:
        print("✗ Missing libraries. Install with: pip install networkx matplotlib")
        return None
    except Exception as e:
        print(f"✗ Error plotting result: {e}")
        traceback.print_exc()
        return None

# --- ESEMPI DI UTILIZZO DELLE NUOVE FUNZIONI ---
# 
# Tutte le funzioni graph query possono essere visualizzate con plot_graph_result()
# che rileva automaticamente il tipo di output e lo visualizza appropriatamente.
#
# # Shortest Path
# path = find_shortest_path(db, 'PR_Q9H609', 'DOID_162', direction='ANY')
# plot_graph_result(path)  # auto-detect shortest_path
# plot_graph_result(path, layout='path', save_path='shortest_path.png')  # layout lineare
# 
# # K Shortest Paths
# paths = find_k_shortest_paths(db, 'PR_Q9H609', 'DOID_162', k=10, direction='OUTBOUND')
# plot_graph_result(paths)  # plotta il primo path
# plot_graph_result(paths, path_index=2)  # plotta il terzo path
# 
# # All Shortest Paths
# all_paths = find_all_shortest_paths(db, 'PR_Q9H609', 'DOID_162')
# plot_graph_result(all_paths, path_index=0)
# 
# # Pattern Matching: gene → protein → disease
# patterns = find_pattern(db, ['gene', 'protein', 'disease'], limit=50)
# plot_graph_result(patterns)  # combina tutti i pattern in un grafo
# 
# # Common Neighbors
# common = find_common_neighbors(db, 'PR_Q9H609', 'PR_P12345')
# plot_graph_result(common)  # mostra nodi A, B e vicini comuni
# 
# # Node Neighbors (lista di nodi)
# neighbors = get_node_neighbors(db, 'PR_Q9H609', depth=2, filter_type='gene')
# plot_graph_result(neighbors)
#
# # Subgraph/Traversal (usa center_node, quindi rileva automaticamente)
# subgraph = traverse_node_centric_graph(db, node_key='PR_Q9H609', depth=2)
# plot_graph_result(subgraph)  # equivalente a plot_subgraph()

#%%

#%% Test Graph Retrieval and Visualization

if __name__ == "__main__":
    # --- ESEMPI DI UTILIZZO ---
    # Esempio 1: get_node_centric_graph (1-hop)
    query = 'PR_Q9H609'
    subgraph = get_node_centric_graph(db, node_key=query)
    # subgraph = get_node_centric_graph(db, filters={'label': 'TP53'}, edge_limit=100)
    plot_subgraph(subgraph, layout='kamada_kawai', highlight_center=True, save_path=f'{query}_subgraph.png')

#%%
if __name__ == "__main__":
    # Esempio 2: traverse_node_centric_graph (multi-hop)
    db = setup_arangodb_connection()
    query = 'PR_Q9H609'
    create_arango_graph(db, graph_name='PKT_graph')  # crea il grafo se non esiste
    subgraph = traverse_node_centric_graph(db, node_key=query, depth=4, limit=50)
    # subgraph = traverse_node_centric_graph(db, filters={'class_code': 'EntrezID'}, depth=3, direction='OUTBOUND')
    plot_subgraph(subgraph, layout='spring', highlight_center=True, save_path=f'{query}_traversal.png')

#%%
if __name__ == "__main__":
    # Esempio 3: Visualizzare il sottografo
    db = setup_arangodb_connection()
    query = 'PR_Q9H609'
    G = plot_subgraph(subgraph, layout='spring', highlight_center=True)
    
    # oppure con traversal multi-hop
    subgraph = traverse_node_centric_graph(db, node_key='PR_Q9H609', depth=2, limit=100)
    G = plot_subgraph(subgraph, layout='shell', save_path='subgraph.png')


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

if __name__ == "__main__":
    # serche for _key contains "H7C1B8"
    patterns_to_test = "Q96HL8"
    patterns_to_test = patterns_to_test.strip().split("\n")

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

"""
Oltre al **traversal** già implementato, ArangoDB supporta diverse query su grafi:

### 1. **Shortest Path** (percorso più breve)
```aql
FOR v, e IN OUTBOUND SHORTEST_PATH 'nodes/A' TO 'nodes/B' GRAPH 'PKT_graph'
    RETURN {vertex: v, edge: e}
```

### 2. **K Shortest Paths** (k percorsi più brevi)
```aql
FOR path IN OUTBOUND K_SHORTEST_PATHS 'nodes/A' TO 'nodes/B' GRAPH 'PKT_graph'
    LIMIT 5
    RETURN path
```

### 3. **K Paths** (tutti i percorsi fino a k, senza ottimizzazione peso)
```aql
FOR path IN 1..5 OUTBOUND K_PATHS 'nodes/A' TO 'nodes/B' GRAPH 'PKT_graph'
    RETURN path
```

### 4. **All Shortest Paths** (tutti i percorsi minimi)
```aql
FOR path IN OUTBOUND ALL_SHORTEST_PATHS 'nodes/A' TO 'nodes/B' GRAPH 'PKT_graph'
    RETURN path
```

### 5. **Pattern Matching** (subgraph isomorphism)
```aql
FOR v1 IN nodes
    FILTER v1.bioentity_type == 'gene'
    FOR v2, e1 IN 1..1 OUTBOUND v1 GRAPH 'PKT_graph'
        FILTER v2.bioentity_type == 'protein'
        FOR v3, e2 IN 1..1 OUTBOUND v2 GRAPH 'PKT_graph'
            FILTER v3.bioentity_type == 'disease'
            RETURN {gene: v1, protein: v2, disease: v3}
```

### 6. **Neighbor Query** (vicini diretti)
```aql
FOR v IN 1..1 ANY 'nodes/TP53' GRAPH 'PKT_graph'
    RETURN v
```

### 7. **Common Neighbors** (vicini comuni tra due nodi)
```aql
LET neighbors_a = (FOR v IN 1..1 ANY 'nodes/A' GRAPH 'PKT_graph' RETURN v._key)
LET neighbors_b = (FOR v IN 1..1 ANY 'nodes/B' GRAPH 'PKT_graph' RETURN v._key)
RETURN INTERSECTION(neighbors_a, neighbors_b)
```

### 8. **Prune/Filter durante traversal**
```aql
FOR v, e, p IN 1..3 OUTBOUND 'nodes/start' GRAPH 'PKT_graph'
    PRUNE v.bioentity_type == 'disease'  -- ferma il ramo qui
    FILTER e.predicate_label == 'regulates'
    RETURN p
```

### 9. **Weighted Shortest Path**
```aql
FOR v, e IN OUTBOUND SHORTEST_PATH 'nodes/A' TO 'nodes/B' GRAPH 'PKT_graph'
    OPTIONS {weightAttribute: 'weight'}
    RETURN {v, e}
```

### 10. **Centrality / PageRank** (tramite Pregel o SmartGraphs)
```javascript
// Via API o arangosh
var pregel = require("@arangodb/pregel");
pregel.start("pagerank", "PKT_graph", {maxGSS: 100, resultField: "rank"});
```

---

**Vuoi che implementi una o più di queste come funzioni Python nel tuo arangodb_utils.py?** Le più utili per knowledge graph biomedici sono tipicamente:
- **Shortest Path** (trovare connessioni tra entità)
- **Pattern Matching** (trovare motivi gene→protein→disease)
- **Common Neighbors** (entità che collegano due nodi)
"""