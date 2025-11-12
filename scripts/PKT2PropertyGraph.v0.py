#%%
# TEST ARANGO DB CONNECTION
from arango import ArangoClient
db_name = 'PheKnowLator_REDUX_NO-GO_test10000'

def test_arango_connection():
    """Test connection to ArangoDB"""
    try:
        client = ArangoClient(hosts='http://localhost:8529')
        sys_db = client.db('_system', username='root', password='avocadodb')
        
        # Check if the database exists
        if sys_db.has_database(db_name):
            print(f"✓ Successfully connected to ArangoDB database: {db_name}")
        else:
            print(f"✗ Database {db_name} does not exist")
        
    except Exception as e:
        print(f"✗ Error connecting to ArangoDB: {e}")
test_arango_connection()

# create a new database if it does not exist
def create_arango_database(db_name):
    """Create a new ArangoDB database if it does not exist"""
    try:
        client = ArangoClient(hosts='http://localhost:8529')
        sys_db = client.db('_system', username='root', password='avocadodb')
        
        if not sys_db.has_database(db_name):
            sys_db.create_database(db_name)
            print(f"✓ Created database: {db_name}")
        else:
            print(f"Database {db_name} already exists")
        
    except Exception as e:
        print(f"✗ Error creating database: {e}")
    
create_arango_database(db_name)
test_arango_connection()

#%%
import pandas as pd
file = "PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt.tar.gz"
pkt_file = r'..\data\pkt\builds\\'+file
pkt_df = pd.read_csv(pkt_file)
pkt_df.head().to_csv()
# pkt_df.sample(10000).to_csv(r'..\temp_dir\sample_10000.csv', index=False)
# pkt_df = pd.read_csv(r'..\temp_dir\sample_10000.csv')
#%%
pkt_df.head()#.to_csv()
#%%
# pkt_df data structure (Tabular URI):
#  ,subject,predicate,object
# 0,http://purl.obolibrary.org/obo/PR_Q7Z4G1,http://purl.obolibrary.org/obo/RO_0002436,http://purl.obolibrary.org/obo/PR_Q96L50
# 1,http://purl.obolibrary.org/obo/CHEBI_7553,http://purl.obolibrary.org/obo/RO_0002434,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000586592
# 2,http://purl.obolibrary.org/obo/CHEBI_46024,http://purl.obolibrary.org/obo/RO_0002434,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000468105
# 3,http://purl.obolibrary.org/obo/PR_Q8N2Z9,http://purl.obolibrary.org/obo/RO_0002436,http://purl.obolibrary.org/obo/PR_Q9NVI1-3
# 4,http://www.ncbi.nlm.nih.gov/gene/145773,http://purl.obolibrary.org/obo/RO_0002511,https://uswest.ensembl.org/Homo_sapiens/Transcript/Summary?t=ENST00000557914


# Effective conversion of RDF dataset pkt_df into a Property graph

import json
import re
from collections import defaultdict
from urllib.parse import urlparse

def extract_entity_info(uri):
    """Extract namespace, entity type, and ID from URI"""
    parsed = urlparse(uri)
    
    # Extract namespace
    if parsed.netloc:
        namespace = parsed.netloc
    else:
        namespace = "unknown"
    
    # Extract entity type and ID
    path = parsed.path
    entity_id = uri.split('/')[-1]  # Get the last part
    
    # Load BFO mappings if not already cached
    if not hasattr(extract_entity_info, '_bfo_mappings'):
        extract_entity_info._bfo_mappings = load_bfo_mappings_from_owl()
    
    bfo_mappings = extract_entity_info._bfo_mappings
    
    # Determine entity type based on namespace and patterns
    entity_type = "unknown"
    if "obo/PR_" in uri:
        entity_type = "protein"
    elif "obo/CHEBI_" in uri:
        entity_type = "chemical"
    elif "obo/RO_" in uri:
        entity_type = "relation"
    elif "obo/BFO_" in uri:
        # Use BFO mappings for BFO entities
        bfo_code = entity_id
        entity_type = bfo_mappings.get(bfo_code, f"bfo_{bfo_code}")
    elif "gene/" in uri:
        entity_type = "gene"
    elif "ensembl.org" in uri:
        entity_type = "transcript"
    elif "obo/GO_" in uri:
        entity_type = "gene_ontology"
    elif "obo/HP_" in uri:
        entity_type = "phenotype"
    elif "obo/CL_" in uri:
        entity_type = "cell_type"
    elif "obo/UBERON_" in uri:
        entity_type = "anatomy"
    
    return {
        'namespace': namespace,
        'entity_type': entity_type,
        'entity_id': entity_id,
        'full_uri': uri
    }

def load_ro_mappings_from_owl(owl_file_path="../data/owl/ro.owl"):
    """Load RO relation mappings from ro.owl file"""
    import xml.etree.ElementTree as ET
    
    try:
        # Parse the OWL file
        tree = ET.parse(owl_file_path)
        root = tree.getroot()
        
        # Define namespaces
        namespaces = {
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
            'owl': 'http://www.w3.org/2002/07/owl#',
            'obo': 'http://purl.obolibrary.org/obo/'
        }
        
        relation_mappings = {}
        
        # Find all ObjectProperty declarations
        for obj_prop in root.findall('.//owl:ObjectProperty', namespaces):
            about = obj_prop.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
            
            if about and 'RO_' in about:
                # Extract RO code
                ro_code = about.split('/')[-1]
                
                # Find rdfs:label for this property
                label_elem = obj_prop.find('.//rdfs:label', namespaces)
                if label_elem is not None:
                    label = label_elem.text
                    if label:
                        # Clean up the label - convert to snake_case
                        clean_label = label.lower().replace(' ', '_').replace('-', '_')
                        relation_mappings[ro_code] = clean_label
        
        print(f"Loaded {len(relation_mappings)} RO mappings from {owl_file_path}")
        return relation_mappings
        
    except Exception as e:
        print(f"Error loading RO mappings from OWL file: {e}")
        # Fallback to hardcoded mappings
        return {
            'RO_0002436': 'molecularly_interacts_with',
            'RO_0002434': 'interacts_with',
            'RO_0002511': 'inheres_in',
            'RO_0002212': 'negatively_regulates',
            'RO_0002213': 'positively_regulates',
            'RO_0002211': 'regulates',
            'RO_0002331': 'involved_in',
            'RO_0002202': 'develops_from',
            'RO_0002220': 'adjacent_to'
        }

def load_bfo_mappings_from_owl(owl_file_path="../data/owl/bfo_classes_only.owl"):
    """Load BFO class mappings from bfo_classes_only.owl file"""
    import xml.etree.ElementTree as ET
    
    try:
        # Parse the OWL file
        tree = ET.parse(owl_file_path)
        root = tree.getroot()
        
        # Define namespaces
        namespaces = {
            'rdf': 'http://www.w3.org/1999/02/22-rdf-syntax-ns#',
            'rdfs': 'http://www.w3.org/2000/01/rdf-schema#',
            'owl': 'http://www.w3.org/2002/07/owl#',
            'obo': 'http://purl.obolibrary.org/obo/'
        }
        
        class_mappings = {}
        
        # Find all Class declarations for BFO
        for owl_class in root.findall('.//owl:Class', namespaces):
            about = owl_class.get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}about')
            
            if about and 'BFO_' in about:
                # Extract BFO code
                bfo_code = about.split('/')[-1]
                
                # Find BFO_0000179 annotation (elucidation/label) for this class
                label_elem = owl_class.find('.//obo:BFO_0000179', namespaces)
                if label_elem is not None:
                    label = label_elem.text
                    if label:
                        # Clean up the label - convert to snake_case
                        clean_label = label.lower().replace(' ', '_').replace('-', '_')
                        class_mappings[bfo_code] = clean_label
                else:
                    # Fallback to rdfs:label if BFO_0000179 not found
                    label_elem = owl_class.find('.//rdfs:label', namespaces)
                    if label_elem is not None:
                        label = label_elem.text
                        if label:
                            clean_label = label.lower().replace(' ', '_').replace('-', '_')
                            class_mappings[bfo_code] = clean_label
        
        print(f"Loaded {len(class_mappings)} BFO class mappings from {owl_file_path}")
        return class_mappings
        
    except Exception as e:
        print(f"Error loading BFO mappings from OWL file: {e}")
        # Fallback to hardcoded mappings for common BFO classes
        return {
            'BFO_0000001': 'entity',
            'BFO_0000002': 'continuant',
            'BFO_0000003': 'occurrent',
            'BFO_0000004': 'independent_continuant',
            'BFO_0000006': 'spatial_region',
            'BFO_0000008': 'temporal_region',
            'BFO_0000015': 'process',
            'BFO_0000016': 'disposition',
            'BFO_0000017': 'realizable_entity',
            'BFO_0000019': 'quality',
            'BFO_0000020': 'specifically_dependent_continuant',
            'BFO_0000023': 'role',
            'BFO_0000024': 'fiat_object_part',
            'BFO_0000027': 'object_aggregate',
            'BFO_0000029': 'site',
            'BFO_0000030': 'object',
            'BFO_0000031': 'generically_dependent_continuant',
            'BFO_0000034': 'function',
            'BFO_0000035': 'process_boundary',
            'BFO_0000038': 'one_dimensional_temporal_region',
            'BFO_0000040': 'material_entity',
            'BFO_0000141': 'immaterial_entity'
        }

def extract_predicate_info(predicate_uri):
    """Extract relation type from predicate URI using complete RO ontology"""
    # Load RO mappings from OWL file (cached)
    if not hasattr(extract_predicate_info, '_ro_mappings'):
        extract_predicate_info._ro_mappings = load_ro_mappings_from_owl()
    
    relation_mappings = extract_predicate_info._ro_mappings
    
    # Extract RO code from predicate
    ro_match = re.search(r'RO_(\d+)', predicate_uri)
    if ro_match:
        ro_code = f"RO_{ro_match.group(1)}"
        return relation_mappings.get(ro_code, ro_code)
    
    return predicate_uri.split('/')[-1]

# Convert RDF triples to Property Graph format
def convert_rdf_to_property_graph(df):
    """Convert RDF dataframe to property graph with nodes and edges"""
    
    # Extract unique nodes (subjects and objects)
    nodes = {}
    edges = []
    
    print("Converting RDF triples to property graph...")
    
    for idx, row in df.iterrows():
        subject = row['subject']
        predicate = row['predicate']
        obj = row['object']
        
        # Process subject node
        if subject not in nodes:
            subject_info = extract_entity_info(subject)
            nodes[subject] = {
                'node_id': subject,
                'label': subject_info['entity_type'],
                'namespace': subject_info['namespace'],
                'entity_id': subject_info['entity_id'],
                'full_uri': subject_info['full_uri']
            }
        
        # Process object node
        if obj not in nodes:
            object_info = extract_entity_info(obj)
            nodes[obj] = {
                'node_id': obj,
                'label': object_info['entity_type'],
                'namespace': object_info['namespace'],
                'entity_id': object_info['entity_id'],
                'full_uri': object_info['full_uri']
            }
        
        # Create edge
        edge = {
            'edge_id': f"edge_{idx}",
            'source': subject,
            'target': obj,
            'relationship': extract_predicate_info(predicate),
            'predicate_uri': predicate
        }
        edges.append(edge)
        
        if idx % 10000 == 0:
            print(f"Processed {idx} triples...")
    
    return list(nodes.values()), edges

# Perform the conversion
nodes, edges = convert_rdf_to_property_graph(pkt_df)

print(f"Conversion completed!")
print(f"Total nodes: {len(nodes)}")
print(f"Total edges: {len(edges)}")

# Create summary statistics
node_types = defaultdict(int)
edge_types = defaultdict(int)

for node in nodes:
    node_types[node['label']] += 1

for edge in edges:
    edge_types[edge['relationship']] += 1

print("\nNode type distribution:")
for node_type, count in sorted(node_types.items()):
    print(f"  {node_type}: {count}")

print("\nEdge type distribution:")
for edge_type, count in sorted(edge_types.items()):
    print(f"  {edge_type}: {count}")

#%%
# from urllib.parse import urlparse

# # Extract namespaces from subjects
# def extract_namespace(uri):
#     parsed = urlparse(uri)
#     if parsed.netloc:
#         return parsed.netloc
#     elif parsed.path:
#         # For URIs like http://purl.obolibrary.org/obo/PR_Q7Z4G1
#         parts = parsed.path.split('/')
#         if len(parts) >= 2:
#             return '/'.join(parts[:-1])
#     return uri

# # Extract all unique namespaces from subjects
# subject_namespaces = set(pkt_df['subject'].apply(extract_namespace))
# # predicate_namespace
# # object_namespace
# #%%
# print("Subject namespaces:")
# for ns in sorted(subject_namespaces):
#     print(ns)



#%%
# ArangoDB connection

from arango import ArangoClient

def setup_arangodb_connection(db_name):
    """Setup ArangoDB connection and create database if needed"""
    try:
        client = ArangoClient(hosts='http://localhost:8529')
        
        # Connect to system database first
        sys_db = client.db('_system', username='root', password='avocadodb')
        
        # Create database if it doesn't exist
        if not sys_db.has_database(db_name):
            sys_db.create_database(db_name)
            print(f"Created database: {db_name}")
        else:
            print(f"Database {db_name} already exists")
        
        # Connect to the target database
        db = client.db(db_name, username='root', password='avocadodb')
        return db
        
    except Exception as e:
        print(f"Error setting up ArangoDB connection: {e}")
        return None

# Setup the database connection
db = setup_arangodb_connection(db_name)

# Export property graph to different formats

def export_to_json(nodes, edges, output_dir="../temp_dir/"):
    """Export property graph to JSON format"""
    
    # Create comprehensive JSON structure
    property_graph = {
        "metadata": {
            "total_nodes": len(nodes),
            "total_edges": len(edges),
            "node_types": dict(node_types),
            "edge_types": dict(edge_types),
            "export_format": "property_graph_json"
        },
        "nodes": nodes,
        "edges": edges
    }
    
    # Export complete graph
    with open(f"{output_dir}property_graph.json", 'w', encoding='utf-8') as f:
        json.dump(property_graph, f, indent=2, ensure_ascii=False)
    
    # Export nodes separately
    with open(f"{output_dir}nodes.json", 'w', encoding='utf-8') as f:
        json.dump(nodes, f, indent=2, ensure_ascii=False)
    
    # Export edges separately
    with open(f"{output_dir}edges.json", 'w', encoding='utf-8') as f:
        json.dump(edges, f, indent=2, ensure_ascii=False)
    
    print(f"JSON exports completed in {output_dir}")

def export_to_tsv(nodes, edges, output_dir="../temp_dir/"):
    """Export property graph to TSV format"""
    
    # Export nodes to TSV
    nodes_df = pd.DataFrame(nodes)
    nodes_df.to_csv(f"{output_dir}nodes.tsv", sep='\t', index=False)
    
    # Export edges to TSV
    edges_df = pd.DataFrame(edges)
    edges_df.to_csv(f"{output_dir}edges.tsv", sep='\t', index=False)
    
    # Create summary statistics file
    stats_data = {
        'metric': ['total_nodes', 'total_edges'] + [f'nodes_{k}' for k in node_types.keys()] + [f'edges_{k}' for k in edge_types.keys()],
        'value': [len(nodes), len(edges)] + list(node_types.values()) + list(edge_types.values())
    }
    stats_df = pd.DataFrame(stats_data)
    stats_df.to_csv(f"{output_dir}graph_statistics.tsv", sep='\t', index=False)
    
    print(f"TSV exports completed in {output_dir}")

def export_to_cypher(nodes, edges, output_dir="../temp_dir/"):
    """Export property graph to Cypher statements for Neo4j"""
    
    with open(f"{output_dir}property_graph.cypher", 'w', encoding='utf-8') as f:
        # Write node creation statements
        f.write("// Create nodes\n")
        for node in nodes:
            label = node['label'].replace(' ', '_').replace('-', '_')
            entity_id = node['entity_id'].replace("'", "\\'")
            namespace = node['namespace'].replace("'", "\\'")
            full_uri = node['full_uri'].replace("'", "\\'")
            node_id = node['node_id'].replace("'", "\\'")
            
            f.write(f"CREATE (n:{label} {{")
            f.write(f"node_id: '{node_id}', ")
            f.write(f"entity_id: '{entity_id}', ")
            f.write(f"namespace: '{namespace}', ")
            f.write(f"full_uri: '{full_uri}'")
            f.write("});\n")
        
        f.write("\n// Create relationships\n")
        for edge in edges:
            relationship = edge['relationship'].replace(' ', '_').replace('-', '_').upper()
            source_id = edge['source'].replace("'", "\\'")
            target_id = edge['target'].replace("'", "\\'")
            predicate_uri = edge['predicate_uri'].replace("'", "\\'")
            
            f.write(f"MATCH (a {{node_id: '{source_id}'}}) ")
            f.write(f"MATCH (b {{node_id: '{target_id}'}}) ")
            f.write(f"CREATE (a)-[:{relationship} {{")
            f.write(f"edge_id: '{edge['edge_id']}', ")
            f.write(f"predicate_uri: '{predicate_uri}'")
            f.write("}]->(b);\n")
    
    print(f"Cypher export completed in {output_dir}")

# ArangoDB import functionality
#pip install python-arango --upgrade
def export_to_arangodb(nodes, edges, db_connection=None):
    """Export property graph to ArangoDB"""
    
    if db_connection is None:
        print("ArangoDB connection not provided. Skipping ArangoDB export.")
        return
    
    try:
        # Clear existing collections if they exist
        for collection_name in ['nodes', 'edges']:
            if db_connection.has_collection(collection_name):
                db_connection.delete_collection(collection_name)
                print(f"Deleted existing collection: {collection_name}")
        
        # Create vertex collection for nodes
        nodes_collection = db_connection.create_collection('nodes')
        print("Created nodes collection")
        
        # Create edge collection for relationships
        edges_collection = db_connection.create_collection('edges', edge=True)
        print("Created edges collection")
        
        # Insert nodes in batches
        print("Inserting nodes into ArangoDB...")
        batch_size = 1000
        node_docs = []
        
        for i, node in enumerate(nodes):
            # ArangoDB requires _key field
            node_doc = node.copy()
            # Create a safe _key from the entity_id
            safe_key = node['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_')
            node_doc['_key'] = safe_key
            node_docs.append(node_doc)
            
            # Insert in batches
            if len(node_docs) >= batch_size:
                try:
                    nodes_collection.insert_many(node_docs, overwrite=True)
                    print(f"Inserted {i+1} nodes...")
                    node_docs = []
                except Exception as e:
                    print(f"Error inserting nodes batch: {e}")
                    # Try inserting one by one for this batch
                    for doc in node_docs:
                        try:
                            nodes_collection.insert(doc, overwrite=True)
                        except Exception as single_error:
                            print(f"Failed to insert node {doc.get('entity_id', 'unknown')}: {single_error}")
                    node_docs = []
        
        # Insert remaining nodes
        if node_docs:
            try:
                nodes_collection.insert_many(node_docs, overwrite=True)
                print(f"Inserted final batch of {len(node_docs)} nodes")
            except Exception as e:
                print(f"Error inserting final nodes batch: {e}")
                for doc in node_docs:
                    try:
                        nodes_collection.insert(doc, overwrite=True)
                    except Exception as single_error:
                        print(f"Failed to insert node {doc.get('entity_id', 'unknown')}: {single_error}")
        
        # Insert edges in batches
        print("Inserting edges into ArangoDB...")
        edge_docs = []
        
        for i, edge in enumerate(edges):
            edge_doc = edge.copy()
            
            # Create safe keys for source and target
            source_key = nodes[next(j for j, n in enumerate(nodes) if n['node_id'] == edge['source'])]['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_')
            target_key = nodes[next(j for j, n in enumerate(nodes) if n['node_id'] == edge['target'])]['entity_id'].replace(':', '_').replace('/', '_').replace('?', '_').replace('=', '_')
            
            # ArangoDB edge format requires _from and _to
            edge_doc['_from'] = f"nodes/{source_key}"
            edge_doc['_to'] = f"nodes/{target_key}"
            edge_doc['_key'] = f"edge_{i}"
            edge_docs.append(edge_doc)
            
            # Insert in batches
            if len(edge_docs) >= batch_size:
                try:
                    edges_collection.insert_many(edge_docs, overwrite=True)
                    print(f"Inserted {i+1} edges...")
                    edge_docs = []
                except Exception as e:
                    print(f"Error inserting edges batch: {e}")
                    # Try inserting one by one for this batch
                    for doc in edge_docs:
                        try:
                            edges_collection.insert(doc, overwrite=True)
                        except Exception as single_error:
                            print(f"Failed to insert edge {doc.get('_key', 'unknown')}: {single_error}")
                    edge_docs = []
        
        # Insert remaining edges
        if edge_docs:
            try:
                edges_collection.insert_many(edge_docs, overwrite=True)
                print(f"Inserted final batch of {len(edge_docs)} edges")
            except Exception as e:
                print(f"Error inserting final edges batch: {e}")
                for doc in edge_docs:
                    try:
                        edges_collection.insert(doc, overwrite=True)
                    except Exception as single_error:
                        print(f"Failed to insert edge {doc.get('_key', 'unknown')}: {single_error}")
        
        print(f"ArangoDB export completed! Inserted {len(nodes)} nodes and {len(edges)} edges.")
        
        # Create some useful indexes
        try:
            nodes_collection.add_index({'fields': ['label'], 'type': 'hash'})
            nodes_collection.add_index({'fields': ['namespace'], 'type': 'hash'})
            edges_collection.add_index({'fields': ['relationship'], 'type': 'hash'})
            print("Created indexes for better query performance")
        except Exception as e:
            print(f"Warning: Could not create indexes: {e}")
        
    except Exception as e:
        print(f"Error exporting to ArangoDB: {e}")
        import traceback
        traceback.print_exc()

# Export in all formats
print("\nExporting property graph...")
# export_to_json(nodes, edges)
export_to_tsv(nodes, edges)
# export_to_cypher(nodes, edges)

# Export to ArangoDB
# print("\nExporting to ArangoDB...")
# try:
#     export_to_arangodb(nodes, edges, db)
#     print("✓ ArangoDB export completed successfully!")
# except Exception as e:
#     print(f"✗ ArangoDB export failed: {e}")
#     print("Make sure ArangoDB is running and the database exists.")

# Display sample of converted data
print("\nSample nodes:")
for i, node in enumerate(nodes[:5]):
    print(f"  {i+1}. {node}")

print("\nSample edges:")
for i, edge in enumerate(edges[:5]):
    print(f"  {i+1}. {edge}")

# Create schema documentation

def create_schema_documentation(nodes, edges, output_dir="../temp_dir/"):
    """Create comprehensive schema documentation for the property graph"""
    
    schema_doc = {
        "property_graph_schema": {
            "description": "Property graph schema for PheKnowLator RDF dataset conversion with RO and BFO ontology mappings",
            "version": "1.0",
            "ontologies_used": {
                "RO": "Relations Ontology - for relationship mapping",
                "BFO": "Basic Formal Ontology - for entity type classification"
            },
            "node_schema": {
                "description": "Schema for graph nodes (entities)",
                "properties": {
                    "node_id": {
                        "type": "string",
                        "description": "Unique identifier for the node (original URI)",
                        "required": True
                    },
                    "label": {
                        "type": "string", 
                        "description": "Node type/category (includes BFO-mapped entity types)",
                        "required": True,
                        "possible_values": list(node_types.keys())
                    },
                    "namespace": {
                        "type": "string",
                        "description": "Namespace/domain of the entity",
                        "required": True
                    },
                    "entity_id": {
                        "type": "string",
                        "description": "Entity identifier extracted from URI",
                        "required": True
                    },
                    "full_uri": {
                        "type": "string",
                        "description": "Complete original URI",
                        "required": True
                    }
                }
            },
            "edge_schema": {
                "description": "Schema for graph edges (relationships)",
                "properties": {
                    "edge_id": {
                        "type": "string",
                        "description": "Unique identifier for the edge",
                        "required": True
                    },
                    "source": {
                        "type": "string",
                        "description": "Source node URI",
                        "required": True
                    },
                    "target": {
                        "type": "string", 
                        "description": "Target node URI",
                        "required": True
                    },
                    "relationship": {
                        "type": "string",
                        "description": "Type of relationship (RO-mapped)",
                        "required": True,
                        "possible_values": list(edge_types.keys())
                    },
                    "predicate_uri": {
                        "type": "string",
                        "description": "Original predicate URI from RDF",
                        "required": True
                    }
                }
            },
            "statistics": {
                "total_nodes": len(nodes),
                "total_edges": len(edges),
                "node_type_counts": dict(node_types),
                "edge_type_counts": dict(edge_types)
            }
        }
    }
    
    # Save schema as JSON
    with open(f"{output_dir}property_graph_schema.json", 'w', encoding='utf-8') as f:
        json.dump(schema_doc, f, indent=2, ensure_ascii=False)
    
    # Create human-readable README
    readme_content = f"""# Property Graph Schema Documentation

## Overview
This property graph was converted from the PheKnowLator RDF dataset, containing biological entities and their relationships.
The conversion uses mappings from:
- **RO (Relations Ontology)**: For relationship type mapping
- **BFO (Basic Formal Ontology)**: For entity type classification

## Statistics
- **Total Nodes**: {len(nodes):,}
- **Total Edges**: {len(edges):,}

## Node Types
{chr(10).join([f"- **{node_type}**: {count:,} nodes" for node_type, count in sorted(node_types.items())])}

## Edge Types
{chr(10).join([f"- **{edge_type}**: {count:,} relationships" for edge_type, count in sorted(edge_types.items())])}

## File Formats

### JSON Format
- `property_graph.json`: Complete graph structure
- `nodes.json`: Node data only
- `edges.json`: Edge data only
- `property_graph_schema.json`: Schema documentation

### TSV Format
- `nodes.tsv`: Node data in tabular format
- `edges.tsv`: Edge data in tabular format
- `graph_statistics.tsv`: Summary statistics

### Cypher Format
- `property_graph.cypher`: Neo4j import statements

### Ontology Mappings
- `ro_mappings.json` / `ro_mappings.tsv`: RO relationship mappings
- `bfo_mappings.json` / `bfo_mappings.tsv`: BFO entity type mappings

## Usage Examples

### Loading in Python
```python
import json
import pandas as pd

# Load JSON format
with open('property_graph.json', 'r') as f:
    graph = json.load(f)

# Load TSV format
nodes_df = pd.read_csv('nodes.tsv', sep='\\t')
edges_df = pd.read_csv('edges.tsv', sep='\\t')

# Load ontology mappings
with open('ro_mappings.json', 'r') as f:
    ro_mappings = json.load(f)

with open('bfo_mappings.json', 'r') as f:
    bfo_mappings = json.load(f)
```

### Neo4j Import
```cypher
// Run the cypher file
:source property_graph.cypher
```

### ArangoDB Import
Use the provided ArangoDB export function in the conversion script.

## Ontology Information

### BFO Entity Types
BFO (Basic Formal Ontology) provides a standardized way to classify entities:
- **entity**: Top-level entity
- **continuant**: Entities that persist through time
- **occurrent**: Entities that happen or occur
- **material_entity**: Physical entities
- **process**: Temporal activities
- And many more specific classifications

### RO Relationship Types  
RO (Relations Ontology) provides standardized relationship types:
- **interacts_with**: General interaction
- **molecularly_interacts_with**: Molecular-level interaction
- **regulates**: Regulatory relationship
- **part_of**: Structural relationship
- And many more specific relationships
"""
    
    with open(f"{output_dir}README.md", 'w', encoding='utf-8') as f:
        f.write(readme_content)
    
    print(f"Schema documentation created in {output_dir}")

# Generate schema documentation
create_schema_documentation(nodes, edges)

print("\n" + "="*50)
print("PROPERTY GRAPH CONVERSION COMPLETED")
print("="*50)
print(f"✓ Converted {len(pkt_df)} RDF triples")
print(f"✓ Generated {len(nodes)} unique nodes")
print(f"✓ Generated {len(edges)} relationships")
print(f"✓ Applied RO ontology mapping for relationships")
print(f"✓ Applied BFO ontology mapping for entity types")
print(f"✓ Exported to JSON, TSV, and Cypher formats")
print(f"✓ Created ontology mapping files (RO + BFO)")
print(f"✓ Created schema documentation")
print(f"✓ Files saved in ../temp_dir/")
print("="*50)

#%%
# Test the RO and BFO mappings loading
print("Testing RO mappings extraction...")
ro_mappings = load_ro_mappings_from_owl()

# Display some sample mappings
print(f"\nSample RO mappings loaded ({len(ro_mappings)} total):")
for i, (ro_code, label) in enumerate(sorted(ro_mappings.items())[:10]):
    print(f"  {ro_code}: {label}")

if len(ro_mappings) > 10:
    print(f"  ... and {len(ro_mappings) - 10} more")

print("\nTesting BFO mappings extraction...")
bfo_mappings = load_bfo_mappings_from_owl()

# Display some sample BFO mappings
print(f"\nSample BFO mappings loaded ({len(bfo_mappings)} total):")
for i, (bfo_code, label) in enumerate(sorted(bfo_mappings.items())[:10]):
    print(f"  {bfo_code}: {label}")

if len(bfo_mappings) > 10:
    print(f"  ... and {len(bfo_mappings) - 10} more")

# Test extraction with sample predicates from your data
sample_predicates = [
    'http://purl.obolibrary.org/obo/RO_0002436',
    'http://purl.obolibrary.org/obo/RO_0002434', 
    'http://purl.obolibrary.org/obo/RO_0002511'
]

print("\nTesting predicate extraction:")
for pred in sample_predicates:
    result = extract_predicate_info(pred)
    print(f"  {pred} -> {result}")

# Test BFO entity extraction
sample_bfo_entities = [
    'http://purl.obolibrary.org/obo/BFO_0000001',
    'http://purl.obolibrary.org/obo/BFO_0000002',
    'http://purl.obolibrary.org/obo/BFO_0000040'
]

print("\nTesting BFO entity extraction:")
for entity in sample_bfo_entities:
    result = extract_entity_info(entity)
    print(f"  {entity} -> {result['entity_type']} ({result['entity_id']})")

#%%
# Save RO mappings to file for future reference
def save_ro_mappings(mappings, output_dir="../temp_dir/"):
    """Save RO mappings to JSON and TSV files"""
    
    # Save as JSON
    with open(f"{output_dir}ro_mappings.json", 'w', encoding='utf-8') as f:
        json.dump(mappings, f, indent=2, ensure_ascii=False)
    
    # Save as TSV for easy viewing
    mappings_df = pd.DataFrame([
        {'ro_code': k, 'label': v} 
        for k, v in sorted(mappings.items())
    ])
    mappings_df.to_csv(f"{output_dir}ro_mappings.tsv", sep='\t', index=False)
    
    print(f"RO mappings saved to {output_dir}ro_mappings.json and {output_dir}ro_mappings.tsv")

def save_bfo_mappings(mappings, output_dir="../temp_dir/"):
    """Save BFO mappings to JSON and TSV files"""
    
    # Save as JSON
    with open(f"{output_dir}bfo_mappings.json", 'w', encoding='utf-8') as f:
        json.dump(mappings, f, indent=2, ensure_ascii=False)
    
    # Save as TSV for easy viewing
    mappings_df = pd.DataFrame([
        {'bfo_code': k, 'label': v} 
        for k, v in sorted(mappings.items())
    ])
    mappings_df.to_csv(f"{output_dir}bfo_mappings.tsv", sep='\t', index=False)
    
    print(f"BFO mappings saved to {output_dir}bfo_mappings.json and {output_dir}bfo_mappings.tsv")

# Save the loaded RO and BFO mappings
save_ro_mappings(ro_mappings)
save_bfo_mappings(bfo_mappings)

#%%
# ArangoDB verification queries
def verify_arangodb_data(db_connection):
    """Verify that data was loaded correctly into ArangoDB"""
    
    if db_connection is None:
        print("No ArangoDB connection available for verification")
        return
    
    try:
        print("\n" + "="*50)
        print("ARANGODB DATA VERIFICATION")
        print("="*50)
        
        # Count nodes by type
        aql_query = """
        FOR doc IN nodes
        COLLECT label = doc.label WITH COUNT INTO count
        RETURN {label: label, count: count}
        """
        result = db_connection.aql.execute(aql_query)
        
        print("\nNode counts by type:")
        for item in result:
            print(f"  {item['label']}: {item['count']}")
        
        # Count edges by relationship type
        aql_query = """
        FOR doc IN edges
        COLLECT relationship = doc.relationship WITH COUNT INTO count
        RETURN {relationship: relationship, count: count}
        """
        result = db_connection.aql.execute(aql_query)
        
        print("\nEdge counts by relationship type:")
        for item in result:
            print(f"  {item['relationship']}: {item['count']}")
        
        # Sample query: find proteins that interact with chemicals
        aql_query = """
        FOR edge IN edges
        FILTER edge.relationship == "interacts_with"
        FOR source IN nodes
        FILTER source._key == SPLIT(edge._from, "/")[1]
        FILTER source.label == "protein"
        FOR target IN nodes
        FILTER target._key == SPLIT(edge._to, "/")[1]
        FILTER target.label == "chemical"
        LIMIT 5
        RETURN {
            protein: source.entity_id,
            chemical: target.entity_id,
            relationship: edge.relationship
        }
        """
        result = db_connection.aql.execute(aql_query)
        
        print("\nSample protein-chemical interactions:")
        for item in result:
            print(f"  {item['protein']} {item['relationship']} {item['chemical']}")
        
        print("\n✓ ArangoDB data verification completed")
        
    except Exception as e:
        print(f"Error verifying ArangoDB data: {e}")

# Run verification after export
if db is not None:
    verify_arangodb_data(db)

#%%
