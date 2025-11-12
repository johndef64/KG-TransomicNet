#%%
# import needed libraries
import networkx as nx
import os

from rdflib import Graph, Namespace, URIRef, BNode, Literal
from tqdm import tqdm

from pkt_kg.utils import *  # provides access to helper functions

import pickle
def save_graph_as_gpickle(graph, filename):
    # This function saves a MultiDiGraph as a gpickle file
    with open(filename, 'wb') as file:
        pickle.dump(graph, file)
def load_graph_from_gpickle(filename):
    # This function loads a MultiDiGraph from a gpickle file
    with open(filename, 'rb') as file:
        graph = pickle.load(file)
    return graph
def save_graph_as_graphml(graph, filename):
    # This function saves the graph in GraphML format
    nx.write_graphml(graph, filename)

# import built-in namespaces
from rdflib.namespace import OWL, RDF, RDFS

# create namespaces
obo = Namespace('http://purl.obolibrary.org/obo/')
entrez = Namespace('http://www.ncbi.nlm.nih.gov/gene/')

type(obo)
obo.Ladies
# %%


# Set global variables
write_location = '../temp_dir/'
if not os.path.exists(write_location):
    os.mkdir(write_location)

# download data to the temp_dir directory
data_urls = [
    # 'http://purl.obolibrary.org/obo/vo.owl',
    # 'https://zenodo.org/records/10055990/files/PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt.tar.gz',#?download=1'
    # 'http://purl.obolibrary.org/obo/ro.owl',
    # 'http://purl.obolibrary.org/obo/chebi.owl',
]

for url in data_urls:
    file_name = url.split('/')[-1]
    if not os.path.exists(write_location + file_name):
        data_downloader(url, write_location)

go_graph = Graph()
go_graph.parse(write_location + 'go.owl', format='xml')

#%%
# get GO:GO relations
# go_relations = set()
# for s, p, o in tqdm(go_graph.triples((None, None, None))):
#     if "GO_" in o and "GO_" in s:
#          go_relations.add((s, p , o))    
#     go_relations.add((s, p ,o))

# %%
# Approach 1: SPARQL query -- note the use of OWL and RDF built-in namespaces

def get_class_query(graph):
    class_query =  graph.query(
        """SELECT DISTINCT ?class
        WHERE {
          ?class rdf:type owl:Class .}
        """, initNs={'rdf': RDF,
                     'owl': OWL}) 
    return class_query

class_query = (get_class_query(go_graph))
print('There are {} owl:Class objects'.format(len(class_query)))
#%%
# Approach 2: Iterate over RDFLib graph object
my_graph = go_graph
classes = list(my_graph.subjects(RDF.type, OWL.Class))
objects = list(my_graph.subjects(RDF.type, OWL.ObjectProperty))
triples = list(my_graph.triples((None, None, None)))

print('There are {} triples, {} owl:Class objects, and {} owl:ObjectProperties objects'.format(len(triples), len(classes), len(objects)))

#%%
# You can use RDFLib to query and get labels from URIs

import rdflib

def get_label_from_uri(uri, graph):
    # SPARQL query for getting the rdfs:label from the URI
    query = f"""
    SELECT ?label
    WHERE {{
      <{uri}> rdfs:label ?label .
    }}
    """
    # Execute the query and return the label if exists
    results = graph.query(query)
    for row in results:
        return str(row.label) # return label as string
    return None

def get_subjects(graph):
    # Use a set to store subjects as they must be unique
    subjects = set()
    for subj in graph.subjects():
        subjects.add(subj)
    return list(subjects)

subjects = get_subjects(go_graph)

uri = rdflib.URIRef('http://purl.obolibrary.org/obo/RO_HOM0000012') # replace with your URI
# urin = rdflib.BNode('N0f9fb0ecfe714c9294dc368d610c077f')
# get_label_from_uri(urin, ro_graph)

labels = []
for uri in tqdm(subjects):
    label = get_label_from_uri(uri, go_graph)
    labels.append(label)
    #print(label)
#%%
labels
#%%
obo = Namespace('http://purl.obolibrary.org/obo/')
obo.GO_0002186
# %%
# get labels for the URI VO_0002186 -- see two different ways to reference a URI
# using obo namespace
# go_0002186_label = go_graph.label(obo.GO_0002186)
# print(str(go_0002186_label))

# get all triples that VO_0002186 participates in
GO_0002186_triples = list(go_graph.triples((obo.GO_0002186, None, None))) +\
                     list(go_graph.triples((None, None, obo.GO_0002186)))

for s, p, o in tqdm(GO_0002186_triples[:100]):
#   print(s, p, o)
    print(str(s).split('/')[-1], 
          str(p).split('/')[-1], 
          str(o).split('/')[-1] + '\n')  # converting entities to str
# %%

#### Obtaining Detailed Network Statistics  

# convert RDFLib graph to Networkx MultiDiGraph
def rdfgraph2multidigraph(rdf_graph):
    nx_graph = nx.MultiDiGraph()
    
    for s, p, o in tqdm(rdf_graph):
        nx_graph.add_edge(s, o, **{'key': p})
    return nx_graph

nx_graph = rdfgraph2multidigraph(go_graph)

# get the number of nodes, edges, and self-loops
def get_nun_node_edges_loops(nx_graph):
    nodes = nx.number_of_nodes(nx_graph)
    edges = nx.number_of_edges(nx_graph)
    self_loops = nx.number_of_selfloops(nx_graph)
    # get degree information
    avg_degree = float(edges) / nodes
    # get network density
    density = nx.density(nx_graph)
    
    print('There are {} nodes, {} edges, and {} self-loop(s)'.format(nodes, edges, self_loops))
    print('The Average Degree is {}'.format(avg_degree))
    print('The density of the graph is: {}'.format(density))
    print('\n')
    return nodes, edges, self_loops, avg_degree

nodes, edges, self_loops, avg_degree = get_nun_node_edges_loops(nx_graph)
# get_nun_node_edges_loops(nx_ppi_graph)

#%%
# get 5 nodes with the highest degress
def get_highests_degree(nx_graph, top = 5, bottom=0):
    top = top + bottom 
    my_list = [(str(x[0]), x[1]) for x in nx_graph.degree]
    n_deg = sorted(my_list, key=lambda x: x[1], reverse=1)[bottom:top-1]
    
    for x in n_deg:
        print('{} (degree={})'.format(x[0], x[1]))
    return n_deg

get_highests_degree(nx_graph, 5)
# get_highests_degree(nx_ppi_graph,10, 10000)

#%%

# get connected components -- have to convert MultiDiGraph to undirected graph
def get_connected_components(nx_graph):
    nx_graph_und = nx_graph.to_undirected()
    
    # get connected components
    components = sorted(list(nx.connected_components(nx_graph_und)), key=len, reverse=True)
    cc_content = {x: str(len(components[x])) + ' nodes: ' + ' | '.join(components[x]) if len(components[x]) < 500
                  else len(components[x]) for x in range(len(components))}
    
    for k, v in tqdm(cc_content.items()):
        if isinstance(v, int):
            print('COMPONENT: {} Consists of {} nodes'.format(str(k), str(v)))
        else:
            print('\nCOMPONENT: {} Consists of the following nodes:'.format(str(k)))
            for node in v.split(': ')[-1].split(' | '):
                print(node)
    return cc_content

cc_content = get_connected_components(nx_graph)
print('There are {} connected components'.format(len(cc_content)))

#%%
# get shortest path from GO_0002186
# go_0002186_path = nx.single_source_shortest_path(nx_graph, obo.GO_0002186)
go_0002186_path = nx.single_source_shortest_path(nx_graph, obo.GO_0002186)

for k, v in go_0002186_path.items():
    if k != obo.GO_0002186:
        print('\n{} - {} Path:'.format(str(obo.GO_0002186).split('/')[-1], str(k).split('/')[-1]))
        for i in v:
            print(i)
# %%

# save multidigraph version of graph
#nx.write_gpickle(nx_graph, write_location + 'vo_NetworkxMultiDiGraph.gpickle')
save_graph_as_gpickle(nx_graph, write_location + 'vo_NetworkxMultiDiGraph.gpickle')

# Example usage:
# Define file path
file_path = write_location + 'go_NetworkxMultiDiGraph.gpickle'

# Load graph from a gpickle file
nx_graph = load_graph_from_gpickle(file_path)

# read in multidigraph version of graph
#nx_graph = nx.read_gpickle(write_location + 'vo_NetworkxMultiDiGraph.gpickle')

def save_rdf_graph(graph, file_path, file_name='vo_graph_data', format='nt'):
    graph.serialize(file_path + f"{file_name}.{format}", format=format)

# save vo_graph as `ntriple` format
save_rdf_graph(go_graph, write_location, format='nt')
# save vo_graph as OWL (i.e. RDF/XML)
save_rdf_graph(go_graph,write_location, format='xml')


###########################
# integrate ontos

# download data to the temp_dir directory
data_urls = [
    #'http://purl.obolibrary.org/obo/ro.owl',
    'http://purl.obolibrary.org/obo/go.owl',
    #'http://purl.obolibrary.org/obo/chebi.owl',
    #'http://purl.obolibrary.org/obo/pr.owl',
]

for url in data_urls:
    file_name = url.split('/')[-1]
    if not os.path.exists(write_location + file_name):
        data_downloader(url, write_location)

#%%
import os
from rdflib import Graph, URIRef
write_location = '../temp_dir/'
# Create a graph and parse the RDF data/ontology files.
g = Graph()
# Here, you can parse from local files or directly from an online source if available.
# Replace these with your file paths or URLs.
ontos =[
    "ro.owl",
    "go.owl",
    #"chebi.owl", 
    #"pr.owl", 
]
def add_ontos_to_graph(g, ontos, root_path = write_location):
    for onto in ontos:
        print(f"Loading {onto}")
        onto_path = root_path+onto
        if not os.path.exists(onto_path):
            onto_path = f"http://purl.obolibrary.org/obo/{onto}"
            
        g.parse(onto_path, format="xml") # Example ontology