#%%

# graph_plots.py - REFACTORED VERSION

import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

# Opzionali (usa try/except così lo script non esplode se mancano)
try:
    import pygraphviz as pgv
except ImportError:
    pgv = None

try:
    import plotly.graph_objects as go
except ImportError:
    go = None

try:
    from pyvis.network import Network
except ImportError:
    Network = None

# =============================================================================
# CARICAMENTO DATI E COSTRUZIONE GRAFO
# =============================================================================

def load_biolink_graph(
    filename="BIOLINK.csv",
    class_col="Class ID",
    label_col="Preferred Label",
    parents_col="Parents",
):
    """
    Carica BIOLINK.csv e restituisce:
    - df_real: dataframe con colonne (Class ID, Preferred Label, Parents)
    - G: DiGraph NetworkX con archi Parent -> Child
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} non trovato nella cartella corrente.")
    
    df = pd.read_csv(filename)
    real_cols = [class_col, label_col, parents_col]
    df_real = df[real_cols].copy()
    
    G = nx.DiGraph()
    for _, row in df_real.iterrows():
        class_id = row[class_col]
        parents = (
            str(row[parents_col]).split("|")
            if pd.notna(row[parents_col])
            else []
        )
        for parent in parents:
            if parent:
                G.add_edge(parent, class_id)
    
    return df_real, G


def compute_levels(G):
    """
    Calcola livelli gerarchici (come in read.py).
    """
    levels = {node: 0 for node in G.nodes()}
    roots = [n for n in G.nodes() if G.in_degree(n) == 0]
    
    for root in roots:
        try:
            distances = nx.single_source_shortest_path_length(G, root)
            for node, dist in distances.items():
                levels[node] = max(levels[node], dist)
        except Exception:
            pass
    
    return levels


# =============================================================================
# HELPER FUNCTION PER NETWORKX PLOTTING (ELIMINA RIDONDANZE)
# =============================================================================

def _draw_nx_graph(
    G,
    pos,
    title,
    out_file=None,
    figsize=(16, 12),
    with_labels=False,
    node_size=30,
    node_color="lightblue",
    edge_color="gray",
    alpha=0.6,
    arrowsize=8,
    font_size=6,
    label_dict=None,
):
    """
    Helper function per disegnare grafi NetworkX con matplotlib.
    Centralizza la logica comune ed elimina duplicazioni.
    
    Parameters:
    -----------
    G : nx.Graph
        Il grafo da visualizzare
    pos : dict
        Posizioni dei nodi {node: (x, y)}
    title : str
        Titolo del grafico
    with_labels : bool
        Se True, mostra le labels dei nodi
    label_dict : dict, optional
        Dizionario custom {node_id: label_text} per le labels
    """
    plt.figure(figsize=figsize)
    
    # Prepara labels se richieste
    if with_labels and label_dict is None:
        # Default: usa l'ultima parte dell'URI troncata a 20 caratteri
        label_dict = {n: str(n).split("/")[-1][:20] for n in G.nodes()}
    
    nx.draw(
        G,
        pos,
        with_labels=with_labels,
        labels=label_dict if with_labels else None,
        node_size=node_size,
        node_color=node_color,
        edge_color=edge_color,
        alpha=alpha,
        arrows=True,
        arrowsize=arrowsize,
        font_size=font_size,
    )
    
    plt.title(title)
    plt.tight_layout()
    
    if out_file:
        plt.savefig(out_file, dpi=80, bbox_inches="tight")
    
    plt.show()


# =============================================================================
# NETWORKX + MATPLOTLIB PLOTTING
# =============================================================================

def plot_nx_spring(
    G, 
    out_file=None, 
    figsize=(12, 12), 
    with_labels=False,
    label_dict=None
):
    """
    Layout spring (force-directed).
    """
    pos = nx.spring_layout(G, k=0.15, iterations=20)
    _draw_nx_graph(
        G,
        pos,
        title="Ontology - Spring Layout",
        out_file=out_file,
        figsize=figsize,
        with_labels=with_labels,
        label_dict=label_dict,
        node_size=10,
        arrowsize=5,
    )


def plot_nx_multipartite(
    G, 
    levels=None, 
    out_file=None, 
    figsize=(20, 12),
    with_labels=False,
    label_dict=None
):
    """
    Layout multipartite per livelli.
    """
    if levels is None:
        levels = compute_levels(G)
    
    # Assegna i livelli come attributi dei nodi
    for node, level in levels.items():
        G.nodes[node]['subset'] = level
    
    pos = nx.multipartite_layout(G, subset_key='subset')
    _draw_nx_graph(
        G,
        pos,
        title="Ontology Hierarchy - Multipartite Layout",
        out_file=out_file,
        figsize=figsize,
        with_labels=with_labels,
        label_dict=label_dict,
        node_size=20,
        node_color="lightblue",
        arrowsize=10,
    )


def plot_nx_circular(
    G, 
    out_file=None, 
    figsize=(16, 16),
    with_labels=False,
    label_dict=None
):
    """
    Layout circolare.
    """
    pos = nx.circular_layout(G)
    _draw_nx_graph(
        G,
        pos,
        title="Ontology - Circular Layout",
        out_file=out_file,
        figsize=figsize,
        with_labels=with_labels,
        label_dict=label_dict,
        node_size=30,
        node_color="lightgreen",
        arrowsize=8,
    )


def plot_nx_shell(
    G, 
    levels=None, 
    out_file=None, 
    figsize=(20, 20),
    with_labels=False,
    label_dict=None
):
    """
    Layout a shell (concentrico per livelli).
    """
    if levels is None:
        levels = compute_levels(G)
    
    shells = []
    for level in sorted(set(levels.values())):
        shell = [n for n, l in levels.items() if l == level]
        if shell:
            shells.append(shell)
    
    pos = nx.shell_layout(G, nlist=shells)
    _draw_nx_graph(
        G,
        pos,
        title="Ontology - Shell Layout (Concentric by Levels)",
        out_file=out_file,
        figsize=figsize,
        with_labels=with_labels,
        label_dict=label_dict,
        node_size=30,
        node_color="lightcoral",
        arrowsize=8,
    )


def plot_nx_kamada_kawai(
    G, 
    out_file=None, 
    figsize=(20, 15), 
    with_labels=False,
    label_dict=None
):
    """
    Layout Kamada-Kawai.
    """
    pos = nx.kamada_kawai_layout(G)
    _draw_nx_graph(
        G,
        pos,
        title="Ontology - Kamada-Kawai Layout",
        out_file=out_file,
        figsize=figsize,
        with_labels=with_labels,
        label_dict=label_dict,
        node_size=25,
        node_color="lightseagreen",
        font_size=6,
    )


def plot_nx_graphviz_layout(
    G, 
    out_file=None, 
    figsize=(20, 15), 
    prog="dot",
    with_labels=False,
    label_dict=None
):
    """
    Layout gerarchico usando graphviz_layout (nx_agraph).
    Richiede pygraphviz installato.
    """
    try:
        pos = nx.nx_agraph.graphviz_layout(G, prog=prog)
        _draw_nx_graph(
            G,
            pos,
            title=f"Ontology Hierarchy - Graphviz {prog} Layout",
            out_file=out_file,
            figsize=figsize,
            with_labels=with_labels,
            label_dict=label_dict,
            node_size=20,
            node_color="lightblue",
            arrowsize=10,
        )
    except Exception as e:
        print(f"Graphviz non disponibile o errore in graphviz_layout: {e}")


# =============================================================================
# PYGRAPHVIZ PLOTTING - CONFIGURAZIONE UNIFICATA
# =============================================================================

# Configurazioni predefinite unificate (elimina ridondanze)
PGV_PRESETS = {
    "basic": {
        'graph': {
            'rankdir': 'LR',
            'splines': 'curved',
            'ranksep': '2.5',
            'nodesep': '0.4',
            'overlap': 'false',
            'fontname': 'Arial',
            'bgcolor': 'white',
            'dpi': '80',
            'margin': '0.5',
            'ratio': 'fill'
        },
        'node': {
            'shape': 'point',
            'width': '0.1',
            'height': '0.1',
            'style': 'filled',
            'fillcolor': '#555555',
            'color': '#555555',
            'fontname': 'Arial',
            'fontsize': '10',
        },
        'edge': {
            'color': '#888888',
            'arrowhead': 'none',
            'penwidth': '0.8'
        }
    },
    "expanded": {
        'graph': {
            'rankdir': 'LR',
            'splines': 'curved',
            'ranksep': '8.0',  # Più espanso orizzontalmente
            'nodesep': '0.1',  # Più compatto verticalmente
            'overlap': 'false',
            'fontname': 'Arial',
            'bgcolor': 'white',
            'dpi': '80',
            'margin': '0.5',
            'ratio': 'fill'
        },
        'node': {
            'shape': 'point',
            'width': '0.1',
            'height': '0.1',
            'style': 'filled',
            'fillcolor': '#555555',
            'color': '#555555',
            'fontname': 'Arial',
            'fontsize': '20',
        },
        'edge': {
            'color': '#888888',
            'arrowhead': 'none',
            'penwidth': '0.8'
        }
    },
    "boxes": {
        'graph': {
            'rankdir': 'LR',
            'splines': 'curved',
            'ranksep': '8.0',  # Più espanso orizzontalmente
            'nodesep': '0.1',  # Più compatto verticalmente
            'overlap': 'false',
            'fontname': 'Arial',
            'bgcolor': 'white',
            'dpi': '80',
            'margin': '0.5',
            'ratio': 'fill'
        },
        'node': {
            'shape': 'box',  # Diverso per boxes
            'style': 'filled',
            'fillcolor': '#E8F4F8',
            'color': '#2B7CE9',
            'fontname': 'Arial',
            'fontsize': '10',
        },
        'edge': {
            'color': '#888888',
            'arrowhead': 'vee',
            'penwidth': '0.8'
        }
    }
}


def apply_pgv_style(A, preset="basic", **custom_attrs):
    """
    Applica uno stile predefinito a un AGraph PyGraphviz.
    Versione semplificata che elimina ridondanze.
    
    Parameters:
    -----------
    A : pygraphviz.AGraph
        Il grafo da stilizzare
    preset : str
        Nome del preset ('basic', 'expanded', 'boxes')
    **custom_attrs : dict
        Attributi custom che sovrascrivono il preset
        Es: graph={'rankdir': 'TB'}, node={'shape': 'circle'}
    """
    if pgv is None:
        raise ImportError("pygraphviz non è installato.")
    
    if preset not in PGV_PRESETS:
        raise ValueError(f"Preset '{preset}' non riconosciuto. Usa: {list(PGV_PRESETS.keys())}")
    
    style = PGV_PRESETS[preset]
    
    # Applica attributi con override custom
    A.graph_attr.update(custom_attrs.get('graph', style['graph']))
    A.node_attr.update(custom_attrs.get('node', style['node']))
    A.edge_attr.update(custom_attrs.get('edge', style['edge']))
    
    return A


def build_pgv_graph(
    df_real,
    class_col="Class ID",
    label_col="Preferred Label",
    parents_col="Parents",
    preset="basic",
    show_labels=False,
    use_xlabel=False,
):
    """
    Costruisce un AGraph PyGraphviz unificato.
    Elimina le 3 funzioni build_pgv_* separate.
    
    Parameters:
    -----------
    show_labels : bool
        Se True, mostra le labels (come 'label' o 'xlabel' a seconda di use_xlabel)
    use_xlabel : bool
        Se True usa xlabel (label esterna), altrimenti label (dentro il nodo)
    """
    if pgv is None:
        raise ImportError("pygraphviz non è installato.")
    
    def get_id(uri_string):
        if pd.isna(uri_string):
            return "Unknown"
        return str(uri_string).split("/")[-1]
    
    def get_label(row):
        if pd.notna(row[label_col]):
            label = str(row[label_col])[:30]  # Tronca a 30 caratteri
        else:
            label = get_id(row[class_col])
        return label
    
    A = pgv.AGraph(directed=True, strict=True)
    A = apply_pgv_style(A, preset=preset)
    
    for _, row in df_real.iterrows():
        node_id = get_id(row[class_col])
        
        if show_labels:
            label_text = get_label(row)
            if use_xlabel:
                A.add_node(node_id, xlabel=f" {label_text}")
            else:
                A.add_node(node_id, label=label_text)
        else:
            A.add_node(node_id)
        
        # Aggiungi archi
        if pd.notna(row[parents_col]):
            parents = str(row[parents_col]).split("|")
            for parent_uri in parents:
                if parent_uri:
                    parent_id = get_id(parent_uri)
                    A.add_edge(parent_id, node_id)
    
    return A


def render_pgv_graph(A, out_file="graph.png", prog="dot"):
    """
    Esegue layout e draw su un AGraph.
    """
    if pgv is None:
        raise ImportError("pygraphviz non è installato.")
    
    A.layout(prog=prog)
    A.draw(out_file)
    print(f"Salvato: {out_file}")


# =============================================================================
# PLOTLY SUNBURST / TREEMAP
# =============================================================================

def plot_plotly_sunburst(
    df_real,
    class_col="Class ID",
    parents_col="Parents",
    out_html="ontology_sunburst.html"
):
    """
    Sunburst interattivo (Plotly).
    """
    if go is None:
        raise ImportError("plotly non è installato.")
    
    labels = []
    parents = []
    values = []
    
    for _, row in df_real.iterrows():
        class_id = str(row[class_col]).split("/")[-1]
        parent_list = (
            str(row[parents_col]).split("|")
            if pd.notna(row[parents_col])
            else [""]
        )
        for parent in parent_list:
            parent_short = parent.split("/")[-1] if parent else ""
            labels.append(class_id)
            parents.append(parent_short)
            values.append(1)
    
    df_tree = pd.DataFrame(
        {"labels": labels, "parents": parents, "values": values}
    ).drop_duplicates()
    
    fig = go.Figure(
        go.Sunburst(
            labels=df_tree["labels"],
            parents=df_tree["parents"],
            values=df_tree["values"],
            branchvalues="total",
        )
    )
    
    fig.update_layout(title="Ontology Sunburst", width=1000, height=1000)
    fig.write_html(out_html)
    print(f"Sunburst salvato in {out_html}")


# =============================================================================
# PYVIS INTERACTIVE GRAPH
# =============================================================================

def build_pyvis_graph(
    filename="BIOLINK.csv",
    class_col="Class ID",
    label_col="Preferred Label",
    def_col="Definitions",
    parents_col="Parents",
    out_html="ontology_interactive.html",
):
    """
    Grafo interattivo con PyVis (layout gerarchico LR).
    """
    if Network is None:
        raise ImportError("pyvis non è installato.")
    
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} non trovato.")
    
    import re
    
    df = pd.read_csv(filename)
    
    def clean_id(text):
        if pd.isna(text):
            return "Unknown"
        text = str(text)
        # Gestisce eventuale formato markdown [url](url)
        match = re.search(r"\((.*?)\)", text)
        if match:
            uri = match.group(1)
        else:
            uri = text.replace("[", "").replace("]", "")
        return uri.split("/")[-1]
    
    def get_parents(parents_str):
        if pd.isna(parents_str):
            return []
        return [clean_id(p) for p in str(parents_str).split("|")]
    
    G = nx.DiGraph()
    
    for _, row in df.iterrows():
        node_id = clean_id(row[class_col])
        label = (
            row[label_col]
            if pd.notna(row[label_col])
            else node_id
        )
        definition = (
            str(row[def_col])
            if def_col in df.columns and pd.notna(row[def_col])
            else "No definition available."
        )
        
        tooltip_html = f"<b>{label}</b><br>ID: {node_id}<br>{definition}"
        
        G.add_node(
            node_id,
            label=label,
            title=tooltip_html,
            group="biolink",
        )
        
        parents = get_parents(row[parents_col]) if parents_col in df.columns else []
        for parent in parents:
            if parent:
                G.add_edge(parent, node_id)
    
    net = Network(
        height="85vh", width="100%", bgcolor="#ffffff", font_color="black", directed=True
    )
    
    net.from_nx(G)
    
    import json
    options = {
        "layout": {
            "hierarchical": {
                "enabled": True,
                "direction": "LR",
                "sortMethod": "directed",
                "levelSeparation": 250,
                "nodeSpacing": 80,
            }
        },
        "nodes": {
            "shape": "dot",
            "size": 10,
            "font": {"size": 14, "face": "Arial"},
            "color": {
                "border": "#2B7CE9",
                "background": "#97C2FC",
            },
        },
        "edges": {
            "color": {"color": "#848484", "highlight": "#ff0000"},
            "smooth": {
                "type": "cubicBezier",
                "forceDirection": "horizontal",
                "roundness": 0.4,
            },
            "arrows": {"to": {"enabled": True, "scaleFactor": 0.5}},
        },
        "physics": {
            "hierarchicalRepulsion": {
                "centralGravity": 0.0,
                "springLength": 100,
                "springConstant": 0.01,
                "nodeDistance": 80,
                "damping": 0.09,
            },
            "solver": "hierarchicalRepulsion",
        },
        "interaction": {
            "hover": True,
            "navigationButtons": True,
            "keyboard": True,
        },
    }
    
    net.set_options(json.dumps(options))
    net.save_graph(out_html)
    print(f"Grafo interattivo salvato in {out_html}")


# =============================================================================
# ESEMPIO DI USO (main)
# =============================================================================

if __name__ == "__main__":
    # Carica dati e grafo base
    df_real, G = load_biolink_graph()
    
    # NetworkX - ORA CON LABELS ABILITATE!
    if False:
        plot_nx_spring(G, out_file="nx_ontology_spring.png", with_labels=True)
        
        levels = compute_levels(G)
        plot_nx_multipartite(G, levels, out_file="nx_ontology_multipartite.png", with_labels=True)
        plot_nx_circular(G, out_file="nx_ontology_circular.png", with_labels=True)
        plot_nx_shell(G, levels, out_file="nx_ontology_shell.png", with_labels=True)
        plot_nx_kamada_kawai(G, out_file="nx_ontology_kamada_kawai.png", with_labels=True)
        plot_nx_graphviz_layout(G, out_file="nx_ontology_graphviz_dot.png", prog="dot", with_labels=True)
    
    # PyGraphviz - ORA UNIFICATO CON LABELS
    if pgv is not None:
        # Versione senza labels (punti minimal)
        A_minimal = build_pgv_graph(df_real, preset="basic", show_labels=False)
        render_pgv_graph(A_minimal, "pgv_ontology_minimal.png")
        
        # Versione con xlabel (label esterne)
        A_xlabel = build_pgv_graph(df_real, preset="expanded", show_labels=True, use_xlabel=True)
        render_pgv_graph(A_xlabel, "pgv_ontology_xlabel.png")
        
        # Versione con boxes e label interne
        A_boxes = build_pgv_graph(df_real, preset="boxes", show_labels=True, use_xlabel=False)
        render_pgv_graph(A_boxes, "pgv_ontology_boxes.png")
    
    # Plotly Sunburst
    if go is not None:
        plot_plotly_sunburst(df_real)
    
    # PyVis interactive
    if Network is not None:
        build_pyvis_graph()

# %%
