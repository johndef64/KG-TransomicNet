#%%
import pandas as pd
import pygraphviz as pgv
import os

# 1. Load and Preprocess Data
# ---------------------------------------------------------
filename = 'BIOLINK.csv'

# Handle potential missing file for the sake of the script running
if not os.path.exists(filename):
    print(f"Error: {filename} not found. Please ensure the file is in the same directory.")
    exit()

df = pd.read_csv(filename)



# Select relevant columns
real_cols = ['Class ID', 'Preferred Label', 'Parents']
df_real = df[real_cols].copy()


# Helper function to extract a short name from a URI or return the label
def get_node_label(row):
    if pd.notna(row['Preferred Label']):
        return row['Preferred Label']
    # Fallback to extracting the last part of the URI if label is missing
    return str(row['Class ID']).split('/')[-1]

def get_id(uri_string):
    # Extract unique ID from URI (e.g., "PathologicalAnatomicalExposure")
    if pd.isna(uri_string):
        return "Unknown"
    return str(uri_string).split('/')[-1]

# 2. Initialize Graphviz AGraph
# ---------------------------------------------------------
# strict=True removes duplicate edges
G = pgv.AGraph(directed=True, strict=True)

# 3. Set Graph Attributes for "Dendrogram" Style
# ---------------------------------------------------------
G.graph_attr.update({
    'rankdir': 'LR',        # Left-to-Right layout (Horizontal)
    'splines': 'curved',    # Curved lines like the screenshot
    'ranksep': '20',       # Increase separation between ranks (horizontal space for labels)
    'nodesep': '0.1',       # Vertical separation between nodes
    'overlap': 'false',     # Prevent overlap
    'fontname': 'Arial',    # Clean font
    'bgcolor': 'white'
})

# 4. Set Default Node & Edge Attributes
# ---------------------------------------------------------
G.node_attr.update({
    'shape': 'point',       # Small dot for the node
    'width':  '0.2',         # Width of the dot
    'height': '0.2',        # Height of the dot
    'style': 'filled',
    'fillcolor': '#555555', # Dark gray/black dots
    'color': '#555555',
    'fontname': 'Arial',
    'fontsize': '50',
})

G.edge_attr.update({
    'color': '#888888',     # Grey edges
    'arrowhead': 'none',    # No arrows (typical for dendrograms)
    'penwidth': '1'       # Thin lines
})

# 5. Build the Graph
# ---------------------------------------------------------
for index, row in df_real.iterrows():
    # Current Node
    node_id = get_id(row['Class ID'])
    node_label = get_node_label(row)
    
    # Add the node with an external label (xlabel) to mimic the text floating next to the dot
    # Note: We use 'xlabel' so the dot remains a dot, and text is placed beside it.
    # For specific alignment, sometimes 'label' with shape='plaintext' is used, 
    # but 'shape=point' + 'xlabel' matches the screenshot's "dot+text" style best.
    G.add_node(node_id, xlabel=f"  {node_label}", labelloc="c") 

    # Parents
    if pd.notna(row['Parents']):
        parents = str(row['Parents']).split("|")
        for parent_uri in parents:
            parent_id = get_id(parent_uri)
            # Add edge from Parent -> Child
            G.add_edge(parent_id, node_id)

# 6. Output
# ---------------------------------------------------------
compact_timestamp = pd.Timestamp.now().strftime("%Y%m%d%H%M%S")
output_filename = f'{compact_timestamp}_biolink_hierarchy.png'
print(f"Generating {output_filename}...")
# Layout using 'dot' engine (hierarchical)
G.layout(prog='dot')
G.draw(output_filename)
print("Done!    ")

# %%
# Funzione ottimizzata per ottenere la gerarchia dei discendenti
def get_descendants_hierarchy(df, root_filter, max_depth=3):
    """
    Ottiene tutti i discendenti di un nodo fino a una profondità massima.
    
    Parameters:
    -----------
    df : pd.DataFrame
        DataFrame con colonne 'Class ID' e 'Parents'
    root_filter : str
        Stringa da cercare nei genitori (es. 'NamedThing')
    max_depth : int
        Numero di livelli di discendenti da includere
    
    Returns:
    --------
    pd.DataFrame
        DataFrame con tutti i nodi filtrati
    """
    # Trova nodi radice
    current_level = df[df['Parents'].str.contains(root_filter, na=False)].copy()
    all_descendants = [current_level]
    
    # Itera per ogni livello di profondità
    for _ in range(max_depth):
        current_ids = current_level['Class ID'].tolist()
        if not current_ids:
            break
        
        # Trova figli del livello corrente
        next_level = df[df['Parents'].isin(current_ids)].copy()
        if next_level.empty:
            break
            
        all_descendants.append(next_level)
        current_level = next_level
    
    # Combina tutti i livelli
    result = pd.concat(all_descendants, ignore_index=True).drop_duplicates()
    return result


# show all entity at the same level of NamedThing
# get the parent of NamedThing
parent_of_namedthing = df_real[df_real['Class ID'].str.contains('NamedThing', na=False)]['Parents'].values
print("Parent(s) of NamedThing:", parent_of_namedthing)


# Crea il subset con la gerarchia di NamedThing (3 livelli di profondità)
SUBSETS = ['NamedThing', 'Entity', 'Association']
SUBSET_NAME = 'NamedThing'  # Cambia qui per scegliere il subset
df_subset = get_descendants_hierarchy(df_real, SUBSET_NAME, max_depth=10)


from mygraphviz_v1 import build_pgv_graph, render_pgv_graph
from datetime import datetime

preset = "expanded"
preset = "boxes"
A_xlabel = build_pgv_graph(df_subset, 
                           preset=preset, 
                           show_labels=True, use_xlabel=True)
compact_timestamp = datetime.now().strftime("%Y%m%d%H%M%S")

render_pgv_graph(A_xlabel, f"{compact_timestamp}_pgv_ontology_{preset}.png")




#%%
# 4. Set Default Node & Edge Attributes
# --------------------------------------------------------
G = pgv.AGraph(directed=True, strict=True)
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


for index, row in df_subset.iterrows():
    # Root nodes labels must be added

    # Current Node
    node_id = get_id(row['Class ID'])
    node_label = get_node_label(row)
    
    # Add the node with an external label (xlabel) to mimic the text floating next to the dot
    # Note: We use 'xlabel' so the dot remains a dot, and text is placed beside it.
    # For specific alignment, sometimes 'label' with shape='plaintext' is used, 
    # but 'shape=point' + 'xlabel' matches the screenshot's "dot+text" style best.
    G.add_node(node_id, xlabel=f"  {node_label}", labelloc="c") 

    # Parents
    if pd.notna(row['Parents']):
        parents = str(row['Parents']).split("|")
        for parent_uri in parents:
            parent_id = get_id(parent_uri)
            # Add edge from Parent -> Child
            G.add_edge(parent_id, node_id)

output_filename = f'biolink_{SUBSET_NAME.lower()}_hierarchy.png'
print(f"Generating {output_filename}...")

# Layout using 'dot' engine (hierarchical)
G.layout(prog='dot')
G.draw(output_filename, format='png')



#%% ==========================================================================
# Interactive Visualization with PyVis

import pandas as pd
import networkx as nx
from pyvis.network import Network
import re
import os

# 1. Caricamento e Pulizia Dati
# ---------------------------------------------------------
filename = 'biolink.csv'
if not os.path.exists(filename):
    print(f"Errore: {filename} non trovato.")
    exit()

df = pd.read_csv(filename)

# Funzione per pulire gli ID (gestisce sia URL puri che formato markdown [url](url))
def clean_id(text):
    if pd.isna(text): return "Unknown"
    text = str(text)
    # Se è formato markdown [url](url), estrai l'url tra parentesi tonde
    match = re.search(r'\((.*?)\)', text)
    if match:
        uri = match.group(1)
    else:
        # Rimuovi eventuali parentesi quadre residue
        uri = text.replace('[', '').replace(']', '')
    
    # Ritorna l'ultimo pezzo dell'URI (es. "NamedThing")
    return uri.split('/')[-1]

# Funzione per pulire i genitori
def get_parents(parents_str):
    if pd.isna(parents_str): return []
    # Gestisce separatore pipe, pulisce ogni genitore
    return [clean_id(p) for p in parents_str.split('|')]

# 2. Costruzione del Grafo NetworkX
# ---------------------------------------------------------
G = nx.DiGraph()

for index, row in df.iterrows():
    node_id = clean_id(row['Class ID'])
    label = row['Preferred Label'] if pd.notna(row['Preferred Label']) else node_id
    definition = str(row['Definitions']) if pd.notna(row['Definitions']) else "No definition available."
    
    # Aggiungi nodo con metadati HTML per il tooltip (title)
    # 'label' è ciò che si vede nel nodo, 'title' è il tooltip
    tooltip_html = f"<b>{label}</b><br><i>{node_id}</i><br><br>{definition}"
    
    G.add_node(node_id, label=label, title=tooltip_html, group='biolink')
    
    # Aggiungi archi
    parents = get_parents(row['Parents'])
    for parent in parents:
        if parent: # evita stringhe vuote
            G.add_edge(parent, node_id)

# 3. Configurazione PyVis Interattiva
# ---------------------------------------------------------
# height='90vh' usa quasi tutto lo schermo verticale
net = Network(height='85vh', width='100%', bgcolor='#ffffff', font_color='black', directed=True)

# Carica il grafo NetworkX
net.from_nx(G)

# CONFIGURAZIONE AVANZATA DELLE OPZIONI
# Impostiamo un layout gerarchico (hierarchical) Left-to-Right per imitare il dendrogramma
options = {
    "layout": {
        "hierarchical": {
            "enabled": True,
            "direction": "LR",        # Left to Right
            "sortMethod": "directed", # Dispone in base alla direzione delle frecce
            "levelSeparation": 250,   # Spazio orizzontale tra i livelli
            "nodeSpacing": 150        # Spazio verticale tra i nodi
        }
    },
    "nodes": {
        "shape": "dot",
        "size": 10,
        "font": {"size": 14, "face": "Arial"},
        "color": {"border": "#2B7CE9", "background": "#97C2FC"}
    },
    "edges": {
        "color": {"color": "#848484", "highlight": "#ff0000"},
        "smooth": {"type": "cubicBezier", "forceDirection": "horizontal", "roundness": 0.4},
        "arrows": {"to": {"enabled": True, "scaleFactor": 0.5}}
    },
    "physics": {
        "hierarchicalRepulsion": {
            "centralGravity": 0.0,
            "springLength": 100,
            "springConstant": 0.01,
            "nodeDistance": 150,
            "damping": 0.09
        },
        "solver": "hierarchicalRepulsion"
    },
    "interaction": {
        "hover": True,
        "navigationButtons": True,
        "keyboard": True
    }
}

# Converti il dizionario opzioni in JSON stringa valido per JS
import json
net.set_options(json.dumps(options))

# Aggiungi menu di filtro (opzionale, utile per cercare nodi specifici)
# net.show_buttons(filter_=['physics', 'layout']) # Scommenta per mostrare i controlli UI completi

# 4. Salvataggio e Apertura
# ---------------------------------------------------------
output_file = "biolink_interactive.html"
net.save_graph(output_file)

print(f"Generato {output_file}. Aprilo nel tuo browser.")

# Opzionale: apre automaticamente il file se sei in locale
try:
    import webbrowser
    webbrowser.open(output_file)
except:
    pass


