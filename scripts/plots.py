#%%
import plotly.graph_objects as go
import networkx as nx
import math, json
import numpy as np
from collections import defaultdict


def _repel_nodes(pos, node_sizes_px, plot_w=1450, plot_h=980,
                 iterations=200, k_repel=0.004):
    """
    Post-process layout positions to push overlapping nodes apart.
    node_sizes_px: dict {node: size in px}
    Positions are in layout units; we convert node radii to the same space.
    """
    nodes = list(pos.keys())
    coords = np.array([pos[n] for n in nodes], dtype=float)

    xs = coords[:, 0]; ys = coords[:, 1]
    x_range = xs.max() - xs.min() or 1
    y_range = ys.max() - ys.min() or 1

    # radius of each node in layout units
    radii = np.array([
        ((node_sizes_px[n] / 2 / (plot_w * 0.85) * x_range) +
         (node_sizes_px[n] / 2 / (plot_h * 0.85) * y_range)) / 2
        for n in nodes
    ])

    for _ in range(iterations):
        for i in range(len(nodes)):
            for j in range(i + 1, len(nodes)):
                dx = coords[i, 0] - coords[j, 0]
                dy = coords[i, 1] - coords[j, 1]
                dist = math.sqrt(dx * dx + dy * dy) or 1e-9
                min_dist = (radii[i] + radii[j]) * 1.25  # 25% padding
                if dist < min_dist:
                    overlap = (min_dist - dist) / 2
                    nx_ = dx / dist; ny_ = dy / dist
                    coords[i, 0] += nx_ * overlap * k_repel * 100
                    coords[i, 1] += ny_ * overlap * k_repel * 100
                    coords[j, 0] -= nx_ * overlap * k_repel * 100
                    coords[j, 1] -= ny_ * overlap * k_repel * 100

    return {n: (coords[i, 0], coords[i, 1]) for i, n in enumerate(nodes)}


def plot_kg_schema(schema,
                   output_path="kg_schema_plot.png",
                   show=True,
                   filter_schema=True,
                   min_degree=10):
    """
    Plotta lo schema generico del KG come grafo diretto.

    Args:
        schema: lista di dict con chiavi:
                source_class_code, predicate_class_code,
                target_class_code, count
                (output di extract_generic_schema_from_files)
        output_path: path dove salvare il PNG
        show: se True apre il plot nel browser (interattivo)
        min_degree: soglia minima — edge con count aggregato < min_degree
                    e nodi con degree totale < min_degree vengono esclusi
    """

    # ── Costruisce grafo ───────────────────────────────────────────────────
    G = nx.DiGraph()
    all_nodes = set()

    if filter_schema:
        schema = [r for r in schema if r['count'] > 100]
        TO_REMOVE = [
            "NCBITaxon",
            "anatomy",
            "UBERON",
            "CLO",
            "ECTO",
            "ENVO",
            # "CL",
            "MOD",
            "GNO",
            "MPATH",
            "CARO",
            ]
        for code in TO_REMOVE:
            schema = [r for r in schema if r['source_class_code'] != code and r['target_class_code'] != code]

    # ── Aggrega edge counts per coppia (s, t) ─────────────────────────────
    edge_counts = defaultdict(int)
    for r in schema:
        s, t = r['source_class_code'], r['target_class_code']
        if s != t:
            edge_counts[(s, t)] += r['count']

    # drop edges below min_degree
    edge_counts = {k: v for k, v in edge_counts.items() if v >= min_degree}

    # drop nodes whose total degree < min_degree
    node_degree_check = defaultdict(int)
    for (s, t), c in edge_counts.items():
        node_degree_check[s] += c
        node_degree_check[t] += c
    active_nodes = {n for n, d in node_degree_check.items() if d >= min_degree}
    edge_counts = {(s, t): c for (s, t), c in edge_counts.items()
                   if s in active_nodes and t in active_nodes}

    # ricostruisce schema filtrato per self-loop e metadati predicato
    schema_filtered = [r for r in schema
                       if (r['source_class_code'] == r['target_class_code'] and
                           r['source_class_code'] in active_nodes)
                       or (r['source_class_code'], r['target_class_code']) in edge_counts]

    for r in schema_filtered:
        all_nodes.add(r['source_class_code'])
        all_nodes.add(r['target_class_code'])

    for (s, t), c in edge_counts.items():
        all_nodes.add(s); all_nodes.add(t)
        if G.has_edge(s, t):
            G[s][t]['count'] += c
        else:
            G.add_edge(s, t, count=c)

    schema = schema_filtered

    node_weighted_degree = defaultdict(int)
    for r in schema:
        node_weighted_degree[r['source_class_code']] += r['count']
        node_weighted_degree[r['target_class_code']] += r['count']
    max_deg = max(node_weighted_degree.values())

    pos = nx.kamada_kawai_layout(G, scale=1.6)
    _sizes_px = {n: 30 + 58 * (node_weighted_degree[n] / max_deg) ** 0.42 for n in all_nodes}
    pos = _repel_nodes(pos, _sizes_px)

    # ── Colori automatici per nodi e predicati ─────────────────────────────
    palette = [
        "#26C6DA", "#FFA726", "#AB47BC", "#EF5350",
        "#66BB6A", "#42A5F5", "#FFCA28", "#FF7043",
        "#F06292", "#80CBC4", "#A5D6A7", "#CE93D8",
    ]
    unique_nodes = sorted(all_nodes)
    unique_preds = sorted(set(r['predicate_class_code'] for r in schema))
    node_color_map = {n: palette[i % len(palette)] for i, n in enumerate(unique_nodes)}
    pred_color_map = {p: palette[i % len(palette)] for i, p in enumerate(unique_preds)}

    # ── Edge traces per predicato ──────────────────────────────────────────
    pred_edge_map = defaultdict(list)
    for r in schema:
        pred_edge_map[r['predicate_class_code']].append(r)

    all_traces = []
    annotations = []

    for idx, (pred, pedges) in enumerate(pred_edge_map.items()):
        color = pred_color_map[pred]
        ex, ey = [], []
        for r in pedges:
            s, t = r['source_class_code'], r['target_class_code']
            if s == t:
                continue
            x0, y0 = pos[s]; x1, y1 = pos[t]
            ox = (y1 - y0) * 0.15; oy = -(x1 - x0) * 0.15
            mx, my = (x0+x1)/2 + ox, (y0+y1)/2 + oy
            ex += [x0, mx, x1, None]
            ey += [y0, my, y1, None]
        all_traces.append(go.Scatter(
            x=ex, y=ey, mode='lines',
            line=dict(width=2.5, color=color), opacity=0.55,
            name=pred, legendgroup=f"pred_{pred}",
            legendgrouptitle=dict(text="<b>Predicates</b>") if idx == 0 else None,
            hoverinfo='none'
        ))

    # Frecce direzionali — si fermano al bordo del nodo target
    # Calcola il raggio di ogni nodo in unità di layout (serve la scala del plot)
    # node_radius_lu = raggio nodo in unità di layout proporzionale alla dimensione px
    plot_w, plot_h = 1450, 980
    xs_all = [pos[n][0] for n in pos]; ys_all = [pos[n][1] for n in pos]
    x_range = max(xs_all) - min(xs_all) or 1
    y_range = max(ys_all) - min(ys_all) or 1

    def node_radius_lu(node):
        """Raggio approssimato in unità di layout per il nodo dato."""
        size_px = 30 + 58 * (node_weighted_degree[node] / max_deg) ** 0.42
        # converti px → unità layout (metà dimensione marker / metà asse in px)
        rx = size_px / 2 / (plot_w * 0.85) * x_range
        ry = size_px / 2 / (plot_h * 0.85) * y_range
        return (rx + ry) / 2  # raggio medio

    for r in schema:
        s, t = r['source_class_code'], r['target_class_code']
        if s == t:
            continue
        color = pred_color_map[r['predicate_class_code']]
        x0, y0 = pos[s]; x1, y1 = pos[t]
        dx, dy = x1-x0, y1-y0
        L = math.sqrt(dx**2 + dy**2) or 1
        # arretra il punto di partenza dal bordo del nodo source
        src_r = node_radius_lu(s)
        ax = x0 + dx/L * src_r
        ay = y0 + dy/L * src_r
        # arretra il punto di arrivo al bordo del nodo target
        tgt_r = node_radius_lu(t)
        tx = x1 - dx/L * tgt_r
        ty = y1 - dy/L * tgt_r
        annotations.append(dict(
            ax=ax, ay=ay,
            x=tx, y=ty,
            xref='x', yref='y', axref='x', ayref='y',
            showarrow=True, arrowhead=2,
            arrowsize=1.5, arrowwidth=2.2,
            arrowcolor=color, opacity=0.9
        ))

    # Label conteggio edge
    elx, ely, elt, elh = [], [], [], []
    for r in schema:
        s, t = r['source_class_code'], r['target_class_code']
        if s == t:
            continue
        x0, y0 = pos[s]; x1, y1 = pos[t]
        ox = (y1 - y0) * 0.15; oy = -(x1 - x0) * 0.15
        elx.append((x0+x1)/2 + ox); ely.append((y0+y1)/2 + oy)
        c = r['count']
        elt.append(f"{c/1000:.1f}k" if c >= 1000 else str(c))
        elh.append(f"<b>{s}</b> –[{r['predicate_class_code']}]→ <b>{t}</b><br>{c:,} edges")
    all_traces.append(go.Scatter(
        x=elx, y=ely, mode='text', text=elt,
        textfont=dict(size=10, color='rgba(30,30,30,0.85)', family='Courier New'),
        hovertext=elh, hoverinfo='text', showlegend=False
    ))

    # Self-loop badges
    for r in schema:
        if r['source_class_code'] != r['target_class_code']:
            continue
        x, y = pos[r['source_class_code']]
        annotations.append(dict(
            x=x+0.08, y=y+0.2, xref='x', yref='y',
            text=f"↺ {r['predicate_class_code']}  {r['count']/1000:.1f}k",
            showarrow=False,
            font=dict(size=11, color=pred_color_map[r['predicate_class_code']]),
            bgcolor='rgba(240,240,240,0.92)',
            bordercolor='rgba(0,0,0,0.2)',
            borderwidth=1, borderpad=4,
        ))

    # Legenda nodi
    for i, node in enumerate(unique_nodes):
        all_traces.append(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=13, color=node_color_map[node],
                        line=dict(width=2, color='rgba(255,255,255,0.3)')),
            name=node, legendgroup=f"node_{node}",
            legendgrouptitle=dict(text="<b>Node class_code</b>") if i == 0 else None,
            showlegend=True
        ))

    # Nodi principali
    nx_l, ny_l, nt_l, nh_l, ns_l, nc_l = [], [], [], [], [], []
    for node in all_nodes:
        x, y = pos[node]
        nx_l.append(x); ny_l.append(y); nt_l.append(node)
        deg = node_weighted_degree[node]
        out_c = sum(r['count'] for r in schema if r['source_class_code'] == node)
        in_c  = sum(r['count'] for r in schema if r['target_class_code'] == node)
        nh_l.append(f"<b>{node}</b><br>Weighted degree: {deg:,}<br>Out: {out_c:,}  In: {in_c:,}")
        ns_l.append(30 + 58 * (deg / max_deg) ** 0.42)
        nc_l.append(node_color_map[node])
    all_traces.append(go.Scatter(
        x=nx_l, y=ny_l, mode='markers+text',
        marker=dict(size=ns_l, color=nc_l,
                    line=dict(width=2.0, color='rgba(0,0,0,0.35)'), opacity=0.95),
        text=nt_l, textposition='middle center',
        textfont=dict(size=14, color='#111111', family='Arial Black'),
        hovertext=nh_l, hoverinfo='text', showlegend=False
    ))

    fig = go.Figure(
        data=all_traces,
        layout=go.Layout(
            title=dict(
                text=(
                    "Knowledge Graph Schema  —  class_code topology<br>"
                    "<span style='font-size:14px;font-weight:normal;color:#555555'>"
                    "Node size ∝ weighted degree  ·  Edge labels = count  ·  Arrows = direction"
                    "</span>"
                ),
                x=0.5, xanchor='center', font=dict(size=22, color='#111111')
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white', paper_bgcolor='white',
            font=dict(color='#111111', family='Arial'),
            hovermode='closest',
            legend=dict(
                bgcolor='rgba(245,245,245,0.95)',
                bordercolor='rgba(0,0,0,0.15)', borderwidth=1,
                font=dict(size=12.5, color='#111111'), tracegroupgap=10,
                x=1.01, y=1.0, xanchor='left', yanchor='top',
            ),
            annotations=annotations,
            width=1450, height=980,
            margin=dict(l=40, r=210, t=120, b=40),
        )
    )

    fig.write_image(output_path, scale=2)
    print(f"✓ Schema plot salvato in {output_path}")
    if show:
        fig.show()
    return fig


def plot_kg_schema_one_type(schema,
                            output_path="kg_schema_plot.png",
                            show=True,
                            min_count=0,
                            min_degree=10,
                            nodes_to_remove=None):
    """
    Plot KG schema with a single predicate type (e.g. only "RO" or only "RDF").
    Edges are drawn as straight arrows with width ∝ log(count).
    Nodes are colored by class and sized by weighted degree.

    Args:
        schema:           list of dicts with keys:
                          source_class_code, predicate_class_code,
                          target_class_code, count
        output_path:      where to save the PNG
        show:             if True opens interactive plot in browser
        min_count:        drop edges with count <= this threshold
        min_degree:       drop edges whose aggregated count < this value;
                          also drop nodes whose total degree (sum of all
                          connected edge counts) < this value
        nodes_to_remove:  list of class_codes to exclude (default list applied if None)
    """
    if nodes_to_remove is None:
        nodes_to_remove = [
            "NCBITaxon", "anatomy", "UBERON", "CLO", "ECTO",
            "ENVO", "MOD", "GNO", "MPATH", "CARO",
        ]

    schema = [r for r in schema if r['count'] > min_count]
    for code in nodes_to_remove:
        schema = [r for r in schema
                  if r['source_class_code'] != code and r['target_class_code'] != code]

    pred_type = schema[0]['predicate_class_code'] if schema else "?"

    # ── Build graph (aggregate counts for same s→t pair) ──────────────────
    G = nx.DiGraph()
    edge_counts = defaultdict(int)
    for r in schema:
        s, t = r['source_class_code'], r['target_class_code']
        if s != t:
            edge_counts[(s, t)] += r['count']

    # drop edges below min_degree
    edge_counts = {k: v for k, v in edge_counts.items() if v >= min_degree}

    # drop nodes whose total degree < min_degree
    node_degree = defaultdict(int)
    for (s, t), c in edge_counts.items():
        node_degree[s] += c
        node_degree[t] += c
    active_nodes = {n for n, d in node_degree.items() if d >= min_degree}
    edge_counts = {(s, t): c for (s, t), c in edge_counts.items()
                   if s in active_nodes and t in active_nodes}

    for (s, t), c in edge_counts.items():
        G.add_edge(s, t, count=c)

    all_nodes = set(G.nodes())
    if not all_nodes:
        print("No nodes after filtering.")
        return None

    # ── Weighted degree ────────────────────────────────────────────────────
    node_weighted_degree = defaultdict(int)
    for (s, t), c in edge_counts.items():
        node_weighted_degree[s] += c
        node_weighted_degree[t] += c
    max_deg = max(node_weighted_degree.values()) or 1

    pos = nx.kamada_kawai_layout(G, scale=2.0)
    _sizes_px = {n: 28 + 60 * (node_weighted_degree[n] / max_deg) ** 0.42 for n in all_nodes}
    pos = _repel_nodes(pos, _sizes_px)

    # ── Node colors (one distinct color per node class) ────────────────────
    palette = [
        "#E53935", "#8E24AA", "#1E88E5", "#00897B", "#F4511E",
        "#3949AB", "#00ACC1", "#43A047", "#FB8C00", "#6D4C41",
        "#D81B60", "#039BE5", "#7CB342", "#FFB300", "#546E7A",
        "#EF9A9A", "#CE93D8", "#90CAF9", "#80CBC4", "#FFCC80",
    ]
    unique_nodes = sorted(all_nodes)
    node_color_map = {n: palette[i % len(palette)] for i, n in enumerate(unique_nodes)}

    # ── Layout geometry helpers ────────────────────────────────────────────
    plot_w, plot_h = 1450, 980
    xs_all = [pos[n][0] for n in pos]
    ys_all = [pos[n][1] for n in pos]
    x_range = max(xs_all) - min(xs_all) or 1
    y_range = max(ys_all) - min(ys_all) or 1

    def node_size_px(node):
        return 28 + 60 * (node_weighted_degree[node] / max_deg) ** 0.42

    def node_radius_lu(node):
        sz = node_size_px(node)
        rx = sz / 2 / (plot_w * 0.85) * x_range
        ry = sz / 2 / (plot_h * 0.85) * y_range
        return (rx + ry) / 2

    max_count = max(edge_counts.values()) or 1

    # ── Annotations: one arrow per edge ───────────────────────────────────
    annotations = []
    elx, ely, elt, elh = [], [], [], []

    for (s, t), c in edge_counts.items():
        x0, y0 = pos[s]; x1, y1 = pos[t]
        dx, dy = x1 - x0, y1 - y0
        L = math.sqrt(dx**2 + dy**2) or 1

        ax = x0 + dx / L * node_radius_lu(s)
        ay = y0 + dy / L * node_radius_lu(s)
        tx = x1 - dx / L * node_radius_lu(t)
        ty = y1 - dy / L * node_radius_lu(t)

        # arrow width scaled by log count
        w = 1.2 + 3.0 * math.log1p(c) / math.log1p(max_count)
        color = node_color_map[s]

        annotations.append(dict(
            ax=ax, ay=ay, x=tx, y=ty,
            xref='x', yref='y', axref='x', ayref='y',
            showarrow=True, arrowhead=3,
            arrowsize=1.2, arrowwidth=w,
            arrowcolor=color, opacity=0.75,
        ))

        # edge count label at midpoint
        elx.append((x0 + x1) / 2)
        ely.append((y0 + y1) / 2)
        elt.append(f"{c/1000:.1f}k" if c >= 1000 else str(c))
        elh.append(f"<b>{s}</b> →[{pred_type}]→ <b>{t}</b><br>{c:,} edges")

    # ── Self-loop badges ───────────────────────────────────────────────────
    for r in schema:
        if r['source_class_code'] != r['target_class_code']:
            continue
        node = r['source_class_code']
        if node not in pos:  # node filtered out (no non-self edges)
            continue
        x, y = pos[node]
        c = r['count']
        annotations.append(dict(
            x=x + 0.1, y=y + 0.22, xref='x', yref='y',
            text=f"↺ {c/1000:.1f}k" if c >= 1000 else f"↺ {c}",
            showarrow=False,
            font=dict(size=11, color=node_color_map[r['source_class_code']]),
            bgcolor='rgba(240,240,240,0.92)',
            bordercolor='rgba(0,0,0,0.18)',
            borderwidth=1, borderpad=4,
        ))

    # ── Traces ────────────────────────────────────────────────────────────
    all_traces = []

    # edge count labels
    all_traces.append(go.Scatter(
        x=elx, y=ely, mode='text', text=elt,
        textfont=dict(size=10, color='rgba(40,40,40,0.8)', family='Courier New'),
        hovertext=elh, hoverinfo='text', showlegend=False,
    ))

    # legend entries (one per node)
    for i, node in enumerate(unique_nodes):
        all_traces.append(go.Scatter(
            x=[None], y=[None], mode='markers',
            marker=dict(size=12, color=node_color_map[node],
                        line=dict(width=1.5, color='rgba(0,0,0,0.2)')),
            name=node,
            legendgrouptitle=dict(text="<b>Node class</b>") if i == 0 else None,
            showlegend=True,
        ))

    # nodes
    nx_l, ny_l, nt_l, nh_l, ns_l, nc_l = [], [], [], [], [], []
    for node in all_nodes:
        x, y = pos[node]
        nx_l.append(x); ny_l.append(y); nt_l.append(node)
        deg = node_weighted_degree[node]
        out_c = sum(c for (s, t), c in edge_counts.items() if s == node)
        in_c  = sum(c for (s, t), c in edge_counts.items() if t == node)
        nh_l.append(f"<b>{node}</b><br>Weighted degree: {deg:,}<br>Out: {out_c:,}  In: {in_c:,}")
        ns_l.append(node_size_px(node))
        nc_l.append(node_color_map[node])

    all_traces.append(go.Scatter(
        x=nx_l, y=ny_l, mode='markers+text',
        marker=dict(size=ns_l, color=nc_l,
                    line=dict(width=2.0, color='rgba(0,0,0,0.3)'), opacity=0.95),
        text=nt_l, textposition='middle center',
        textfont=dict(size=13, color='#111111', family='Arial Black'),
        hovertext=nh_l, hoverinfo='text', showlegend=False,
    ))

    fig = go.Figure(
        data=all_traces,
        layout=go.Layout(
            title=dict(
                text=(
                    f"KG Schema  —  predicate: <b>{pred_type}</b><br>"
                    "<span style='font-size:13px;font-weight:normal;color:#666'>"
                    "Node size ∝ weighted degree  ·  Arrow width ∝ log(count)  ·  Color = source node"
                    "</span>"
                ),
                x=0.5, xanchor='center', font=dict(size=21, color='#111111'),
            ),
            xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
            plot_bgcolor='white', paper_bgcolor='white',
            font=dict(color='#111111', family='Arial'),
            hovermode='closest',
            legend=dict(
                bgcolor='rgba(245,245,245,0.95)',
                bordercolor='rgba(0,0,0,0.15)', borderwidth=1,
                font=dict(size=12, color='#111111'),
                x=1.01, y=1.0, xanchor='left', yanchor='top',
            ),
            annotations=annotations,
            width=1450, height=980,
            margin=dict(l=40, r=200, t=110, b=40),
        )
    )

    fig.write_image(output_path, scale=2)
    print(f"✓ Schema plot salvato in {output_path}")
    if show:
        fig.show()
    # print min degree warning
    if min_degree > 0:
        print(f"⚠️  Nodes with total degree < {min_degree} were removed.")
    return fig
