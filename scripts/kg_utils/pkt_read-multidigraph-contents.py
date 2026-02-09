# this script read and show sample of contents of this multidigraph
#%%
import networkx as nx
import pickle
import tarfile
import os
from  pkt_utils import *
set_working_directory()

#%%
# get the file name fron the working this, thew only file that contqaint "NetworkxMultiDiGraph"
file_name = [f for f in os.listdir() if "NetworkxMultiDiGraph" in f][0]


print(f"Loading MultiDiGraph from: {file_name}")
print("=" * 80)

# Extract and load the graph
with tarfile.open(file_name, 'r:gz') as tar:
    # Get the first member (the .gpickle file)
    members = tar.getmembers()
    print(f"Archive contains: {[m.name for m in members]}")
    
    # Extract and load the pickle file
    extracted_file = tar.extractfile(members[0])
    graph = pickle.load(extracted_file)

print(f"\nGraph type: {type(graph)}")
print(f"Number of nodes: {graph.number_of_nodes()}")
print(f"Number of edges: {graph.number_of_edges()}")
print("=" * 80)

# Sample nodes
print("\nSample of 10 nodes:")
print("-" * 80)
for i, node in enumerate(list(graph.nodes())[:10]):
    node_data = graph.nodes[node]
    print(f"{i+1}. {node}")
    if node_data:
        print(f"   Attributes: {node_data}")
    print()

# Sample edges
print("\nSample of 10 edges:")
print("-" * 80)
for i, (source, target, key, data) in enumerate(list(graph.edges(keys=True, data=True))[:10]):
    print(f"{i+1}. {source} -> {target}")
    print(f"   Key: {key}")
    print(f"   Attributes: {data}")
    print()

# Graph statistics
print("\nGraph Statistics:")
print("-" * 80)
print(f"Is directed: {graph.is_directed()}")
print(f"Is multigraph: {graph.is_multigraph()}")

# Degree statistics (for first 10 nodes)
print("\nDegree statistics for first 10 nodes:")
print("-" * 80)
for i, node in enumerate(list(graph.nodes())[:10]):
    in_degree = graph.in_degree(node)
    out_degree = graph.out_degree(node)
    print(f"{i+1}. {node}")
    print(f"   In-degree: {in_degree}, Out-degree: {out_degree}")

# Edge types/relations (if available in edge attributes)
print("\nEdge relation types (sample):")
print("-" * 80)
edge_types = set()
for _, _, _, data in graph.edges(keys=True, data=True):
    if 'label' in data:
        edge_types.add(data['label'])
    elif 'relation' in data:
        edge_types.add(data['relation'])
    elif 'type' in data:
        edge_types.add(data['type'])
    if len(edge_types) >= 20:  # Limit to first 20 unique types
        break

for i, edge_type in enumerate(sorted(edge_types)[:20], 1):
    print(f"{i}. {edge_type}")

#%%
# Find the top 1000 hub nodes with most connections and save to CSV
import pandas as pd

print("\nFinding top 1000 hub nodes...")
print("=" * 80)

# Calculate total degree (in + out) for each node
node_degrees = []
for node in graph.nodes():
    in_deg = graph.in_degree(node)
    out_deg = graph.out_degree(node)
    total_deg = in_deg + out_deg
    
    # Get node attributes
    node_attrs = dict(graph.nodes[node])
    
    node_info = {
        'node_id': node,
        'in_degree': in_deg,
        'out_degree': out_deg,
        'total_degree': total_deg,
        **node_attrs  # Unpack all node attributes
    }
    node_degrees.append(node_info)

# Sort by total degree descending and get top 1000
top_hubs = sorted(node_degrees, key=lambda x: x['total_degree'], reverse=True)[:1000]

# Create DataFrame
df_hubs = pd.DataFrame(top_hubs)

# Save to CSV
output_file = "top_1000_hub_nodes.csv"
df_hubs.to_csv(output_file, index=False)

print(f"Top 1000 hub nodes saved to: {output_file}")
print(f"Columns in CSV: {list(df_hubs.columns)}")

# Display summary of top 10 hubs
print("\nTop 10 hub nodes summary:")
print("-" * 80)
for i, hub in enumerate(top_hubs[:10], 1):
    print(f"{i}. {hub['node_id']}")
    print(f"   Total degree: {hub['total_degree']} (In: {hub['in_degree']}, Out: {hub['out_degree']})")
    # Show a few key attributes if they exist
    for attr in ['name', 'label', 'type', 'category']:
        if attr in hub and hub[attr]:
            print(f"   {attr}: {hub[attr]}")
    print()

print(f"\nAll {len(top_hubs)} hub nodes saved to {output_file}")
print("=" * 80)

print("\n" + "=" * 80)
print("Done!")

#%%

# show the saved top 1000 hub nodes
import pandas as pd
df_hubs = pd.read_csv("top_1000_hub_nodes.csv")
df_hubs
#%%
# load nodeLabel and show line of "http://purl.obolibrary.org/obo/SO_0000673"
import pandas as pd

NodeLabels = pd.read_csv("PKT_NodeLabels_with_metadata_v3.0.2.csv")
#%%

entity_uri = "http://purl.obolibrary.org/obo/SO_0000673"
NodeLabels[NodeLabels['entity_uri'] == f"<{entity_uri}>"] 
#%%
# merge NodeLabels for only nodes in top 1000 hubs
df_hubs['<node_id>'] = df_hubs['node_id'].apply(lambda x: f"<{x}>")
merged_df = pd.merge(df_hubs, NodeLabels, left_on='<node_id>', right_on='entity_uri', how='left')
merged_df.drop(columns=['<node_id>'], inplace=True)
merged_df
#%%
merged_df.columns

small_merged_df = merged_df[['node_id', 'in_degree', 'out_degree','label', 'class_code', "bioentity_type",'description/definition']]
small_merged_df[small_merged_df['bioentity_type']=="rna"]
#%%
# save merged df to csv
merged_df#.to_csv("top_1000_hub_nodes_with_labels.csv", index=False)



#%%  ====================================================================
#load PKT csv knwoledge graph
import pandas as pd
import tarfile
file_name="PKT.nt.tar.gz"

def read_kg_file(file_path):
    print(f"Reading file ID  {file_path}")    
    with tarfile.open(file_path, 'r:gz') as tar:
        members = tar.getmembers()
        extracted_file = tar.extractfile(members[0])
        df = pd.read_csv(extracted_file, sep=' ', header=None, names=['subject', 'predicate', 'object', '.'])
        return df[['subject', 'predicate', 'object']]


PKT = read_kg_file(file_name)
PKT.head()
#%%
# cerca in PKT dqtaframe le righe che contengono "http://purl.obolibrary.org/obo/SO_0000673"
entiry = "SO_0001217"#"SO_0000831" # SO_0000673
# http://purl.obolibrary.org/obo/SO_0001217
# http://purl.obolibrary.org/obo/NCBITaxon_9606
# http://purl.obolibrary.org/obo/CHEBI_23367
entity_uri = "<http://purl.obolibrary.org/obo/" + entiry + ">"
PKT_entity = PKT[(PKT['subject'] == entity_uri) | (PKT['object'] == entity_uri)]
NodeLabels[NodeLabels['entity_uri'] == entity_uri]
#%%
PKT_entity#.head()
#%%
#====================================================================
# get all "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>" relations
# this type of relatioa are used to define the class of an entity
rdf_type = "<http://www.w3.org/1999/02/22-rdf-syntax-ns#type>"

PKT_rdf = PKT[PKT['predicate'] == rdf_type]
# print a sentece with NodeLabels.label for subject and object
for i, row in PKT_rdf.head(10).iterrows():
    subject_label = NodeLabels[NodeLabels['entity_uri'] == row['subject']]['label'].values
    object_label = NodeLabels[NodeLabels['entity_uri'] == row['object']]['label'].values
    subject_label_str = subject_label[0] if len(subject_label) > 0 else "Unknown"
    object_label_str = object_label[0] if len(object_label) > 0 else "Unknown"
    print(f" {subject_label_str} --rdf:type--> {object_label_str}\n")
#%%
print ( f" numeber of unique objects with rdf:type: {PKT_rdf['object'].nunique() }")
NodeLabelRDFob = NodeLabels[NodeLabels['entity_uri'].isin(PKT_rdf['object'])]
NodeLabelRDFob
#%%

print ( f" numeber of unique objects with rdf:type: {PKT_rdf['subject'].nunique() }")
NodeLabelRDFsub = NodeLabels[NodeLabels['entity_uri'].isin(PKT_rdf['subject'])]
NodeLabelRDFsub
#%%
#====================================================================
NodeLabels[NodeLabels['entity_type'] == "NODES"]["entity_uri"].nunique()

# get a subset of PKT with only nodes that dont have rdf:type relation
PKT_no_rdf_type = PKT[~PKT['subject'].isin(PKT_rdf['subject']) & ~PKT['object'].isin(PKT_rdf['subject'])]
# show them in NodeLabels
NodeLabels_no_rdf_type = NodeLabels[NodeLabels['entity_uri'].isin(PKT_no_rdf_type['subject']) | NodeLabels['entity_uri'].isin(PKT_no_rdf_type['object'])]
NodeLabels_no_rdf_type
#%%
# show bioentity_types count of NodeLabels_no_rdf_type
NodeLabels_no_rdf_type['bioentity_type'].value_counts()
#%%
NodeLabels['bioentity_type'].value_counts()

#%%

# # mostra 5 edge che connettono "http://purl.obolibrary.org/obo/SO_0000673"
# entity_uri = "http://purl.obolibrary.org/obo/SO_0000673"
# edges_with_entity = []
# for u, v, k, data in graph.edges(keys=True, data=True):
#     if u == entity_uri or v == entity_uri:
#         edges_with_entity.append((u, v, k, data))
#     if len(edges_with_entity) >= 5:
#         break
# #%%
# for i, (u, v, k, data) in enumerate(edges_with_entity, 1):
#     print(f"{i}. {u} -> {v}")
#     print(f"   Key: {k}")
#     print(f"   Attributes: {data}")
#     print()
