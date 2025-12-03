
#%% read .ttl file
filename = 'biolink-model.owl.ttl'

# load libraries
from rdflib import Graph

# load data
g = Graph()
g.parse(filename, format='turtle')

# print number of triples
print(f"Number of triples: {len(g)}")

rows = []
# optionally, iterate through triples
for subject, predicate, obj in g:
    print(f"{subject} {predicate} {obj}")
    subject = subject.split('/')[-1]
    predicate = predicate.split('/')[-1]
    object_ = obj.split('/')[-1]
    row = {'Subject': subject, 'Predicate': predicate, 'Object': object_}
    rows.append(row)

df = pd.DataFrame(rows)
df

#%%
# filtraggio dei black nodes
from rdflib import BNode, Graph
import pandas as pd 
filename = 'biolink-model.owl.ttl'
# load data
g = Graph()
g.parse(filename, format='turtle')

rows = []
# optionally, iterate through triples
for subject, predicate, obj in g:
    # Salta i blank nodes
    if isinstance(subject, BNode) or isinstance(obj, BNode):
        continue
    
    print(f"{subject} {predicate} {obj}")
    subject = str(subject).split('/')[-1]
    predicate = str(predicate).split('/')[-1]
    object_ = str(obj).split('/')[-1]
    row = {'Subject': subject, 'Predicate': predicate, 'Object': object_}
    rows.append(row)

df = pd.DataFrame(rows)
#%%
df["Subject"].value_counts()
df["Object"].value_counts()
df["Predicate"].value_counts()

# correct this to get unique predicates
df["pred_type"] = df["Predicate"].apply(lambda x: x.split('#')[0])
df["pred_type"].value_counts()

#%%
df_rdf = df[df["pred_type"] == "rdf-schema"]
df_rdf.Predicate.value_counts().to_clipboard()

df[df["Predicate"].str.contains("subPropertyOf")]
# df[df["Subject"].str.contains("agent")]
df[df["Object"].str.contains("related_to")]
df[df["Subject"].str.contains("contributes_to")]

#%%

small_graph= df[df["Subject"].str.contains("contributes_to")]

import pygraphviz as pgv

def show_graph(df_graph):
    # Filtra le righe con valori vuoti o None
    df_graph = df_graph[
        (df_graph['Subject'].notna()) & 
        (df_graph['Object'].notna()) & 
        (df_graph['Predicate'].notna()) &
        (df_graph['Subject'] != '') & 
        (df_graph['Object'] != '') & 
        (df_graph['Predicate'] != '')
    ].copy()
    
    if df_graph.empty:
        print("No valid data to visualize")
        return None

    G = pgv.AGraph(directed=True, strict=True)

    # 4. Set Default Node & Edge Attributes
    # --------------------------------------------------------

    G.graph_attr.update({
        'rankdir': 'LR',
        'splines': 'curved',
        'ranksep': '8.0',       # AUMENTA MOLTO questo valore per espandere orizzontalmente
        'nodesep': '0.1',       # Riduci questo per compattare verticalmente
        'overlap': 'false',
        'fontname': 'Arial',
        'bgcolor': 'white',
        'dpi': '150',           # Risoluzione
        'margin': '0.5',        # Margini
        'ratio': 'fill'         # Riempi lo spazio disponibile
    })

    G.node_attr.update({
        'shape': 'point',
        'width': '0.1',
        'height': '0.1',
        'style': 'filled',
        'fillcolor': '#555555',
        'color': '#555555',
        'fontname': 'Arial',
        'fontsize': '20',        # Font leggermente più piccolo
    })

    G.edge_attr.update({
        'color': '#888888',
        'arrowhead': 'none',
        'penwidth': '0.8'
    })


    for index, row in df_graph.iterrows():
        # Ruse the triple graph dataframe
        G.add_edge(str(row['Subject']), str(row['Object']), label=str(row['Predicate']))
    # Layout and Render
    G.layout(prog='dot')
    import os
    output_path = os.path.join(os.getcwd(), 'small_graph.png')
    G.draw(output_path)
    print(f"Graph saved to: {output_path}")
    # show grqaph
    import matplotlib.pyplot as plt
    img = plt.imread(output_path)
    plt.figure(figsize=(10, 10))
    plt.imshow(img)
    plt.axis('off')
    plt.show()
    return G

G_small = show_graph(small_graph)

#%%
df_core = df[df["pred_type"] == "core"]
df_core
#%%
# find ubjects or object containing related_to
related_to_entities = set()
for subj, pred, obj in g.triples((None, None, None)):
    if str(subj).endswith("related_to"):
        related_to_entities.add(subj)
        related_to_entities.add(obj)
print("Entities related to 'related_to':")
for entity in related_to_entities:
    print(entity)
#%%

# the sub graph from "related_to" 
import pandas as pd

rows = []
related_to_triples = g.triples((None, None, None))
for subj, pred, obj in related_to_triples:
    if str(subj).endswith("related_to"):
        print(f"{subj.split('/')[-1]} {pred.split('/')[-1]} {obj.split('/')[-1]}")
        # get sub dataframe
        subject = subj.split('/')[-1]
        predicate = pred.split('/')[-1]
        object_ = obj.split('/')[-1]
        row = {'Subject': subject, 'Predicate': predicate, 'Object': object_}
        rows.append(row)

df_related_to = pd.DataFrame(rows)
df_related_to
#%%