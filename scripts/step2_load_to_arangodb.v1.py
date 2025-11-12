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

"""
The structure of the JSON data to be loaded are this

edge.json
[
  {
    "edge_id": "edge_0",
    "source_uri": "http://purl.obolibrary.org/obo/PR_P04234",
    "target_uri": "http://purl.obolibrary.org/obo/PR_P17693-1",
    "relationship": "molecularly_interacts_with",
    "predicate_uri": "http://purl.obolibrary.org/obo/RO_0002436"
  },
  {
    "edge_id": "edge_1",
    "source_uri": "http://purl.obolibrary.org/obo/CHEBI_30185",
    "target_uri": "https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000578112",
    "relationship": "interacts_with",
    "predicate_uri": "http://purl.obolibrary.org/obo/RO_0002434"
  },

node.json
[
  {
    "_key": "PR_P04234",
    "node_id": "http://purl.obolibrary.org/obo/PR_P04234",
    "label": "protein",
    "namespace": "purl.obolibrary.org",
    "entity_id": "PR_P04234",
    "full_uri": "http://purl.obolibrary.org/obo/PR_P04234"
  },
  {
    "_key": "PR_P17693-1",
    "node_id": "http://purl.obolibrary.org/obo/PR_P17693-1",
    "label": "protein",
    "namespace": "purl.obolibrary.org",
    "entity_id": "PR_P17693-1",
    "full_uri": "http://purl.obolibrary.org/obo/PR_P17693-1"
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

#you must implement this function that is missing correctly, the main gial i sto load nodes edges and their propertis correcly i aragnodb 

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
                source_key = edge['source_uri'].split('/')[-1].split('?')[0]  # Handle query params
                target_key = edge['target_uri'].split('/')[-1].split('?')[0]
                
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

# --- FUNZIONI GRAFICO ---
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
        
        # 3. Verifica finale
        print("\nFinal verification:")
        check_collections_data(db)

        # 4. Visualizza il grafo
        G = visualize_random_graph(db)
    else:
        print("✗ Failed to connect to ArangoDB. Please check if ArangoDB is running.")
        
# %%


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
