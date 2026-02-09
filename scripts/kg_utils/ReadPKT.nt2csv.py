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
write_location = '../data/pkt/'
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

#%%

pkt_graph = Graph()
pkt_build = "PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt"
# pkt_build = "PheKnowLator_v2.0.0_REDUX_NO-GO.nt"

file_path = write_location + "builds/" + pkt_build
# file_path = write_location + "pkt_DrugGene_subgraph.nt"
pkt_graph.parse(file_path, format='nt')
# %%


# Nome del file RDF
write_location = '../data/ptk/'
pkt_build = 'PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt'
# pkt_build = 'PheKnowLator_v2.0.0_REDUX_NO-GO.nt'
file_path = write_location + "builds/" + pkt_build

# Numero di righe da visualizzare
# num_lines = 10

# Leggere e stampare le prime 100 righe del file
# def show_file_head(file_path, num_lines=10):
#     with open(file_path, 'r') as file:
#         for i, line in enumerate(file):
#             if i < num_lines:
#                 print(line.strip())
#             else:
#                 break
# show_file_head(file_path, num_lines=100)

# %%
# Redux Grapoh, purge ontos

###python
# Define the input and output file paths
file_path = '../temp_dir/PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt'
output_file = '../temp_dir/PheKnowLator_v2.0.0_REDUX_NO-GO.nt'

# List of strings to remove lines containing them
expell = [
    "http://purl.obolibrary.org/obo/MONDO",
    "http://purl.obolibrary.org/obo/HP",
    "http://purl.obolibrary.org/obo/VO",
    "http://purl.obolibrary.org/obo/UBERON",
    
    
    "https://www.ncbi.nlm.nih.gov/snp/rs",
    "http://purl.obolibrary.org/obo/SO",
          
    "http://purl.obolibrary.org/obo/CLO",
    "http://purl.obolibrary.org/obo/CL",
    
    "http://purl.obolibrary.org/obo/GO"]

# Read lines from the input file
with open(file_path, 'r') as file:
    lines = file.readlines()

# Filter out lines containing any of the strings in expell
filtered_lines = [line for line in lines if not any(s in line for s in expell)]

# Write the filtered lines to the output file
with open(output_file, 'w') as file:
    file.writelines(filtered_lines)

#%%

# Efficient Parsing loading of the graph

from rdflib import Graph
import os
from tqdm import tqdm
write_location = '../data/ptk/'

def read_in_chunks(file_path, chunk_size=10000):
    with open(file_path, "r") as file:
        chunk = []
        for line in file:
            chunk.append(line)
            if len(chunk) >= chunk_size:
                yield "".join(chunk)
                chunk = []
        # Yield any remaining lines
        if chunk:
            yield "".join(chunk)


def parse_in_chunks(file_path, chunk_size=10000):
    """Parse file in parallel using multiprocessing."""
    chunks = read_in_chunks(file_path, chunk_size=chunk_size)
    num_chunks = 0
    # Printing the chunks before processing
    for chunk in chunks:
        num_chunks += 1
    print(f"Number of chunks to be processed: {num_chunks}")

    g = Graph()
    with open(file_path, "r") as file:
        lines = []
        #for i, line in enumerate(file):
        #for i, line in enumerate(tqdm(file, desc="Processing chunks")):
        for i, line in enumerate(tqdm(file, total=num_chunks, desc="Processing chunks")):
            lines.append(line)
            if i % chunk_size == 0:
                #print(f"Parsing chunk {i // chunk_size}")
                g.parse(data="\n".join(lines), format="nt")
                lines = []
        # Parse remaining lines
        g.parse(data="\n".join(lines), format="nt")
    return g

# Usa parse_in_chunks per caricare i tuoi dati
write_location = '../data/ptk/'
pkt_build = 'PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt'
# pkt_build = "PheKnowLator_v2.0.0_REDUX_REDUX.nt"
# pkt_build = "PheKnowLator_v2.0.0_REDUX_NO-GO.nt"
file_path = write_location + "builds/" + pkt_build
pkt_graph = parse_in_chunks(file_path)


"""
CPU times: user 8min 22s, sys: 5.35 s, total: 8min 27s
Wall time: 8min 27s
-----
[REDUX on my OLD VERSION pc --> RAM satura]
CPU times: total: 10min 36s
Wall time: 12min 12s
-----
[REDUX on my pc]
CPU times: total: 6min 52s
Wall time: 8min 4s
"""
#%%
pkt_graph = None
#%%
# load file.nt with pandas

import pandas as pd
def load_nt_to_dataframe(file_path):
    # Leggi solo le prime tre colonne, ignorando il punto finale
    df = pd.read_csv(
        file_path,
        sep=' ',
        header=None,
        usecols=[0, 1, 2],
        names=['subject', 'predicate', 'object'],
        engine='python',
        quoting=3
    )
    # Rimuovi il punto finale dalla colonna object, se presente
    df['object'] = df['object'].str.rstrip('.')
    return df

write_location = '../data/ptk/'
pkt_build = 'PheKnowLator_v2.0.0_full_instance_relationsOnly_noOWL_OWLNETS.nt'

# Example usage
file_path = write_location + "builds/" + pkt_build
pkt_df = load_nt_to_dataframe(file_path)


# Display the first few rows of the DataFrame
#%%

# define pkt_df
def remove_angle_brackets(df):
    return df.replace({'<|>': ''}, regex=True)
pkt_df = remove_angle_brackets(pkt_df)
pkt_df.head()

#%%
# save dataframa as csv zipped
def save_dataframe_as_csv(df, file_path):
    # Save the DataFrame to a CSV file with gzip compression
    df.to_csv(file_path, index=False, compression='gzip')
# Example usage
output_file_path = write_location + "builds/" + pkt_build.replace(".nt",'.csv.gz')
save_dataframe_as_csv(pkt_df, output_file_path)

#%%
# Load Graph from CSV file
def load_graph_from_csv(file_path):
    # Load the DataFrame from the CSV file
    df = pd.read_csv(file_path, compression='gzip')
    
    # Create a new RDF graph
    g = Graph()
    
    # Iterate through the DataFrame and add triples to the graph
    for _, row in df.iterrows():
        g.add((URIRef(row['subject']), URIRef(row['predicate']), URIRef(row['object'])))
    
    return g
# Example usage
file_path = output_file_path  # Path to the CSV file
pkt_graph = load_graph_from_csv(file_path)
# Display the number of triples in the graph
print(f"Number of triples in the graph: {len(pkt_graph)}")


#%%
# get basic descriptive statistics
nodes = set(list(pkt_graph.subjects()) + list(pkt_graph.objects()))
rels = set(list(pkt_graph.predicates()))

# print stats
print('Graph Stats: {} triples, {} nodes, {} predicates'.format(len(pkt_graph), len(nodes), len(rels)))

#%%
# DEFINE NAMESPACE: Proteins
prot = Namespace('http://purl.obolibrary.org/obo/PR_')  
# Iterare attraverso tutti gli oggetti e filtrare quelli che appartengono al namespace specificato
pkt_proteins = [obj for obj in pkt_graph.objects() if obj.startswith(prot)]
pkt_proteins

#%%

# Using mygene package to convert UniProt ID to gene symbol

from mygene import MyGeneInfo

# Initialize MyGeneInfo
mg = MyGeneInfo()

def uniprot_to_gene_symbol(uniprot_id):
    # Query for converting UniProt ID to gene symbol
    result = mg.query(uniprot_id, fields='symbol', species='human')
    if result['hits']:
        return result['hits'][0]['symbol']
    return None

# Example usage
uniprot_id = "P12345" # replace with your uniprot_id
gene_symbol = uniprot_to_gene_symbol(uniprot_id)
print(f"Gene symbol for {uniprot_id} is {gene_symbol}")
# Using the UniProt API and BioPython to map Uniprot IDs to gene symbols

#%%
#pip install bio_utils 
from bio_utils import *

# Example usage
uniprot_id = "P12345"
for obj in pkt_proteins[:10]:
    uniprot_id = obj.split('/')[-1].replace('PR_','')
    gene_symbol = uniprot_to_gene_symbol(uniprot_id)
    print(f"Gene symbol for UniProt ID {uniprot_id} is {gene_symbol}")

#%%
from rdflib.namespace import OWL, RDF, RDFS
# get all owl:ObjectProperty objects
pkt_graph_props = list(pkt_graph.subjects(RDF.type, OWL.ObjectProperty))
#pkt_graph_props = list(pkt_graph.subjects(predicate_uri))
for p in tqdm(pkt_graph_props[:100]):
    print(p)

#### Create Subgraph ###

def filters_edge_data(graph, pred_namespace, subj_namespace, obj_namespace):
    """Method takes an input RDFLib graph and filters it using an the input pred_namespace, subj_namespace and obj_namespace variables.
    
    Args:
        graph: An RDFLib Graph object.
        pred_namespace: A URIRef object containing a Relation Ontology relation. 
        subj_namespace: A URIRef object containing information to filter subjects by.
        obj_namespace : A URIRef object containing information to filter objects by.
    
    Returns:
        filtered_edges: 
    """  
     
    filtered_triples = []

    for s, p, o in tqdm(graph):
        #if p == pred_namespace:
        if str(p).startswith(str(pred_namespace)):
            if str(s).startswith(str(subj_namespace)) and str(o).startswith(str(obj_namespace)):
                filtered_triples += [(s, p, o)]
    
    return filtered_triples

#%%

# get gene-drug edges - then keep subjects with 
drug_gene_triples = filters_edge_data(graph=pkt_graph,
                                      pred_namespace=obo.RO_0002434, #'interacts with'
                                      subj_namespace=obo.CHEBI,
                                      obj_namespace=entrez)

print('There are {} drug-gene edges'.format(len(drug_gene_triples)))

#%%
# get drug-disease edges
drug_disease_triples = filters_edge_data(graph=pkt_graph,
                                         pred_namespace=obo.RO_0002606,
                                         subj_namespace=obo.CHEBI,
                                         obj_namespace=obo.MONDO)

print('There are {} drug-disease edges'.format(len(drug_disease_triples)))


#%%
# get gene-disease edges
gene_disease_triples = filters_edge_data(graph=pkt_graph,
                                         pred_namespace=obo.RO_0003302,
                                         subj_namespace=entrez,
                                         obj_namespace=obo.MONDO)

print('There are {} gene-disease edges'.format(len(gene_disease_triples)))

#%%
# combine triples into single graph
filtered_edges = drug_gene_triples + drug_disease_triples + gene_disease_triples
gene_drug_disease_graph = adds_edges_to_graph(Graph(), filtered_edges)

print('The drug-gene-disease Subgraph contains {} triples'.format(len(gene_drug_disease_graph)))
#%%


#### GET TRIPLES: build subgraphs
from myRDFLIB_utils import *
# GET from OBJ FOR "GO" 
subjects, predicates, objects = get_unique_triple_value_cont(pkt_graph, 'GO')
# GET from OBJ FOR "PR" 
subjects, predicates, objects = get_unique_triple_value_cont(pkt_graph, 'PR')



from itertools import islice
import pandas as pd
# Create a dataframe from the first 10 elements from islice
data = list(islice(zip(subjects, predicates, objects), 10))
df = pd.DataFrame(data, columns=['Subject', 'Predicate', 'Object'])
df
# GET NAMESPACES
obo    = Namespace('http://purl.obolibrary.org/obo/')
entrez = Namespace('http://www.ncbi.nlm.nih.gov/gene/')
trans  = Namespace('https://uswest.ensembl.org/Homo_sapiens/Transcript/') #Transcript
embl   = Namespace('http://www.ebi.ac.uk/')  # diseases
ract   = Namespace( 'https://reactome.org/content/detail/')  # Pathways
prot   = Namespace('http://purl.obolibrary.org/obo/PR_')  # Proteins
# GET ALL
subjects, predicates, objects = get_unique_triple_value(pkt_graph)
# for pred in predicates:
#     print(pred)
# Create a dataframe of all the grpah
data = list(zip(subjects, predicates, objects))
pkg_df = pd.DataFrame(data, columns=['Subject', 'Predicate', 'Object'])
#%%
len(predicates),  len(objects), len(subjects)
#entrez
def get_namespaces(subjects):
    subjects_namespaces = set()
    for i in subjects:
        subj = i.replace(i.split('/')[-1], '')
        subjects_namespaces.add(subj)
    return subjects_namespaces
    
sub_namespaces = get_namespaces(subjects)
obj_namespaces = get_namespaces(objects)
obj_namespaces

#%%
def get_predicates_for_subj_obj(graph, namespace):
    # This function prints predicates associated with the given namespace in an RDFLib graph
    predicates_for_sub = set()
    predicates_for_obj = set()
    for subj, pred, obj in tqdm(graph):
        if namespace in subj:  # check if subject includes the namespace
            predicates_for_sub.add(pred)
        if namespace in obj:     
            predicates_for_obj.add(pred)
    return predicates_for_sub, predicates_for_obj
# Show unique predicates associated with the 'entrez' namespace

namespace_pred_sub, namespace_pred_obj = get_predicates_for_subj_obj(pkt_graph, prot)
for pred in namespace_pred_obj:
    print(pred)
#%%

def get_subjects_for_namespace(graph, namespace):
    # This function returns subjects associated with the given namespace in an RDFLib graph
    subjects = set()
    for subj, pred, obj in tqdm(graph):
        if namespace in subj:  # check if subject includes the namespace
            subjects.add(subj)

    return subjects

namespace_entries = get_subjects_for_namespace(pkt_graph, ract)
namespace_entries

# ..............
#%%
import rdflib
def get_objects_for_subject_predicate(graph, subject_ns, predicate):
    # This function retrieves objects' namespaces for a given subject namespace and predicate
    object_namespaces = set()
    objects = set()
    for subj, pred, obj in tqdm(graph):
        if subject_ns in subj and pred == predicate:
            # Get the namespace of the object
            obj_namespace = rdflib.namespace.split_uri(obj)[0]
            object_namespaces.add(obj_namespace)
            objects.add(obj)
    return object_namespaces, objects

object_namespaces, objects= get_objects_for_subject_predicate(pkt_graph, prot, obo.RO_0002436)
objects
#%%
# get gene-drug edges - then keep subjects with 
relation = obo.RO_0002434
#relation = obo.RO_0002436 #molecularly interacts with
relation = obo.RO_0002511 
relation = obo.RO_0000056
relation = obo.RO_0002436
entity_entity_triples = filters_edge_data(graph=pkt_graph,
                                         relation=relation,
                                         #subj=entrez,
                                         #obj=ract
                                         subj=prot,
                                         obj=prot
                                         )

print('There are {} entity-entity edges'.format(len(entity_entity_triples)))
#%%
entity_entity_graph = adds_edges_to_graph(Graph(), entity_entity_triples)

# convert RDFLib graph to Networkx MultiDiGraph
nx_graph_entity_entity = nx.MultiDiGraph()

for s, p, o in tqdm(entity_entity_graph):
    nx_graph_entity_entity.add_edge(s, o, **{'key': p})

label = 'ProtProt'
nt_file      = write_location + f'pkt_{label}_subgraph.nt'
gpickle_file = write_location + f'pkt_{label}_NetworkxMultiDiGraph.gpickle'
graphml_file = write_location + f'pkt_{label}_NetworkxMultiDiGraph.graphml'
# save subgraph
entity_entity_graph.serialize(nt_file, format='nt')
#nx.write_gpickle(nx_graph_drug_gene, write_location + 'pkt_DrugGene_NetworkxMultiDiGraph.gpickle')
#%%
entity_entity_graph = Graph()
entity_entity_graph.parse(nt_file, format='nt')

entity_entity_graph = adds_edges_to_graph(Graph(), drug_gene_triples)

# convert RDFLib graph to Networkx MultiDiGraph
nx_graph_entity_entity = nx.MultiDiGraph()

for s, p, o in tqdm(drug_gene_graph):
    nx_graph_entity_entity.add_edge(s, o, **{'key': p})

# save multidigraph version of graph
#nx.write_gpickle(nx_graph, write_location + 'vo_NetworkxMultiDiGraph.gpickle')
save_graph_as_gpickle(nx_graph_entity_entity, gpickle_file)

# save graphml file
nx_graph_entity_entity = load_graph_from_gpickle(gpickle_file)
save_graph_as_graphml(nx_graph_entity_entity, graphml_file)
#%%
#### Convert Subgraph to NetworkX MultiDiGraph

# convert RDFLib graph to Networkx MultiDiGraph
nx_graph_dgd = nx.MultiDiGraph()

for s, p, o in tqdm(gene_drug_disease_graph):
    nx_graph_dgd.add_edge(s, o, **{'key': p})

#%%

##### Get Graph Descriptives


# get the number of nodes, edges, and self-loops
nodes = nx.number_of_nodes(nx_graph_dgd)
edges = nx.number_of_edges(nx_graph_dgd)
avg_degree = float(edges) / nodes

print('There are {} nodes, {} edges, and has an average degree of {}'.format(nodes, edges, avg_degree))


# get 5 nodes with the highest degress
n_deg = sorted([(str(x[0]), x[1]) for x in  nx_graph_dgd.degree], key=lambda x: x[1], reverse=1)[:6]

for x in n_deg:
    print('{} (degree={})'.format(x[0], x[1]))
# get connected components -- have to convert MultiDiGraph to undirected graph
nx_graph_und = nx_graph_dgd.to_undirected()

# get connected components
components = sorted(list(nx.connected_components(nx_graph_und)), key=len, reverse=True)
cc_content = {x: str(len(components[x])) + ' nodes: ' + ' | '.join(components[x]) if len(components[x]) < 500
              else len(components[x]) for x in range(len(components))}

for k, v in cc_content.items():
    if isinstance(v, int):
        print('COMPONENT: {} Consists of {} nodes'.format(str(k), str(v)))
    else:
        print('\nCOMPONENT: {} Consists of the following nodes:'.format(str(k)))
        for node in v.split(': ')[-1].split(' | '):
            print(node)

#%%
#### Explore Graph


# perform bidirectional search for path between -- epilepsy and valporic acid
nx.bidirectional_shortest_path(nx_graph_dgd, obo.CHEBI_39867, obo.MONDO_0005027)
# look at shortest path between 
shortest_paths = nx.shortest_path(nx_graph_dgd, source=obo.CHEBI_39867)

for k, v in shortest_paths.items():
    if k != obo.CHEBI_39867:
        print('\n{} - {} Path:'.format(str(obo.CHEBI_39867).split('/')[-1], str(k).split('/')[-1]))
        for i in v:
            print(i)

#%%
#### Save Output
# save subgraph
gene_drug_disease_graph.serialize(write_location + 'pkt_DrugGeneDisease_subgraph.nt', format='nt')
nx.write_gpickle(nx_graph_dgd, write_location + 'pkt_DrugGeneDisease_NetworkxMultiDiGraph.gpickle')