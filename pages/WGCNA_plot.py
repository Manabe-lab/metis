import streamlit as st
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import seaborn as sns
from pathlib import Path
import tempfile
import io
from adjustText import adjust_text
import math

# ForceAtlas2 layout implementation - use ForceAtlas2 package only
FA2_AVAILABLE = False

# Try ForceAtlas2Py package (works with Python 3.11+)
try:
    import ForceAtlas2 as FA2Py
    FA2_AVAILABLE = True
except ImportError:
    FA2Py = None

# Try python-igraph as fallback
IGRAPH_AVAILABLE = False
if not FA2_AVAILABLE:
    try:
        import igraph as ig
        IGRAPH_AVAILABLE = True
    except ImportError:
        pass

# Matplotlib用のエッジカラー（タプル形式）
EDGE_COLORS = {
    'regular': (0.4, 0.4, 0.4, 0.3),      # 濃いグレー
    'hub': (1.0, 0.08, 0.58, 0.4),        # 濃いdeeppink
    'highlight': (1.0, 0.55, 0.0, 0.4)     # 濃いorange
}

# SVG用のエッジカラー（文字列形式）
SVG_EDGE_COLORS = {
    'regular': "rgba(100,100,100,0.3)",    # 濃いグレー
    'hub': "rgba(255,20,147,0.4)",         # 濃いdeeppink
    'highlight': "rgba(255,140,0,0.4)"     # 濃いorange
}

# WGCNAで使用される可能性のある全ての色のマッピング
# R WGCNAのlabels2colors関数から抽出
WGCNA_COLORS = {
    # Standard colors
    'white': '#FFFFFF',
    'black': '#000000',
    'grey': '#808080',
    'red': '#FF0000',
    'green': '#00FF00',
    'blue': '#0000FF',
    'yellow': '#FFFF00',
    'magenta': '#FF00FF',
    'cyan': '#00FFFF',
    'brown': '#A52A2A',
    'pink': '#FFC0CB',
    'turquoise': '#40E0D0',
    'midnightblue': '#191970',
    'purple': '#800080',
    'greenyellow': '#ADFF2F',
    'tan': '#D2B48C',
    'salmon': '#FA8072',
    'lightcyan': '#E0FFFF',
    'royalblue': '#4169E1',
    'darkred': '#8B0000',
    'darkgreen': '#006400',
    'darkturquoise': '#00CED1',
    'darkgrey': '#404040',
    'orange': '#FFA500',
    'darkorange': '#FF8C00',
    'lightgreen': '#90EE90',
    'darkmagenta': '#8B008B',
    'darkolivegreen': '#556B2F',
    'lightsteelblue': '#B0C4DE',
    'darkviolet': '#9400D3',
    'skyblue': '#87CEEB',
    'saddlebrown': '#8B4513',
    'steelblue': '#4682B4',
    'paleturquoise': '#AFEEEE',
    'violet': '#EE82EE',
    'darkcyan': '#008B8B',
    'lightyellow': '#FFFFE0',
    'thistle': '#D8BFD8',
    'plum': '#DDA0DD',
    'palevioletred': '#DB7093',
    'orchid': '#DA70D6',
    'sienna': '#A0522D',
    'coral': '#FF7F50',
    'chocolate': '#D2691E',
    'maroon': '#800000',
    'ivory': '#FFFFF0',
    'floralwhite': '#FFFAF0',
    'peru': '#CD853F',
    'gold': '#FFD700',
    'cornsilk': '#FFF8DC',
    'firebrick': '#B22222',
    'forestgreen': '#228B22',
    'indianred': '#CD5C5C',
    'khaki': '#F0E68C',
    'lawngreen': '#7CFC00',
    'lightblue': '#ADD8E6',
    'lightpink': '#FFB6C1',
    'mediumpurple': '#9370DB',
    'navajowhite': '#FFDEAD',
    'navy': '#000080',
    'olive': '#808000',
    'orangered': '#FF4500',
    'rosybrown': '#BC8F8F',
    'seagreen': '#2E8B57',
    'tomato': '#FF6347',
    'wheat': '#F5DEB3',
    # Grey scale colors
    'darkgrey': '#404040',
    'dimgrey': '#696969',
    'lightgrey': '#A9A9A9',
    'silver': '#C0C0C0',
    # Additional colors for large networks
    'antiquewhite': '#FAEBD7',
    'aqua': '#00FFFF',
    'aquamarine': '#7FFFD4',
    'azure': '#F0FFFF',
    'beige': '#F5F5DC',
    'bisque': '#FFE4C4',
    'blanchedalmond': '#FFEBCD',
    'burlywood': '#DEB887',
    'cadetblue': '#5F9EA0',
    'chartreuse': '#7FFF00',
    'cornflowerblue': '#6495ED',
    'crimson': '#DC143C',
    'cyan2': '#00FFFF',
    'cyan3': '#00CED1',
    'darkblue': '#00008B',
    'darkgoldenrod': '#B8860B',
    'darkkhaki': '#BDB76B',
    'darkseagreen': '#8FBC8F',
    'darkslateblue': '#483D8B',
    'darkslategray': '#2F4F4F',
    'deeppink': '#FF1493',
    'deepskyblue': '#00BFFF',
    'dodgerblue': '#1E90FF',
    'hotpink': '#FF69B4',
    'lavender': '#E6E6FA',
    'lavenderblush': '#FFF0F5',
    'lemonchiffon': '#FFFACD',
    'lime': '#00FF00',
    'limegreen': '#32CD32',
    'mediumblue': '#0000CD',
    'mediumorchid': '#BA55D3',
    'mediumseagreen': '#3CB371',
    'mediumslateblue': '#7B68EE',
    'mediumspringgreen': '#00FA9A',
    'mediumturquoise': '#48D1CC',
    'mediumvioletred': '#C71585',
    'mintcream': '#F5FFFA',
    'mistyrose': '#FFE4E1',
    'moccasin': '#FFE4B5',
    'oldlace': '#FDF5E6',
    'olivedrab': '#6B8E23',
    'palegoldenrod': '#EEE8AA',
    'palegreen': '#98FB98',
    'peachpuff': '#FFDAB9',
    'powderblue': '#B0E0E6',
    'rebeccapurple': '#663399',
    'sandybrown': '#F4A460',
    'seashell': '#FFF5EE',
    'slateblue': '#6A5ACD',
    'slategray': '#708090',
    'springgreen': '#00FF7F',
    'teal': '#008080',
    'yellowgreen': '#9ACD32'
}


@st.cache_data
def generate_network_plots(tom_array, gene_names, module_colors, threshold, 
                         selected_modules, n_hubs, highlight_indices, 
                         minimum_component_size, plot_type):
    """Generate and cache plots based on plot type"""
    if plot_type == "Interactive":
        # インタラクティブプロット生成
        html_content = plot_network_svg(
            tom_array, gene_names, module_colors, threshold,
            selected_modules, n_hubs, highlight_indices, minimum_component_size
        )
        # 静的プロット生成（PDFとPNG用）
        fig_static = plot_network(
            tom_array, gene_names, module_colors, threshold,
            selected_modules, n_hubs, highlight_indices, minimum_component_size
        )
        return html_content, fig_static
    
    elif plot_type == "Force-Directed":
        # Force-directedプロット生成
        fig = plot_network(
            tom_array, gene_names, module_colors, threshold,
            selected_modules, n_hubs, highlight_indices, minimum_component_size
        )
        return fig
    
    elif plot_type == "ForceAtlas2":
        # ForceAtlas2プロット生成
        fig = plot_network_forceatlas2(
            tom_array, gene_names, module_colors, threshold,
            selected_modules, n_hubs, highlight_indices, minimum_component_size
        )
        return fig
    
    elif plot_type == "Circular Layout":
        # Circular Layoutプロット生成
        fig = plot_circular_layout(
            tom_array, gene_names, module_colors, 
            selected_modules, threshold
        )
        return fig
    
    elif plot_type == "Hive Plot":
        # Hive Plotプロット生成
        fig = plot_hive_layout(
            tom_array, gene_names, module_colors, 
            selected_modules, threshold,
            max_nodes=500
        )
        return fig

def plot_network_svg(tom_array, gene_names, module_colors=None, threshold=0.1, 
                     selected_modules=None, n_hubs=10, highlight_nodes=None,
                     minimum_component_size=20):
    """Create interactive network visualization with hover that matches force-directed appearance"""
    if tom_array is None:
        return None
    
    # Create graph
    G = nx.Graph()
    
    # Get hub genes before filtering nodes
    hub_indices = get_hub_genes(tom_array, module_colors, selected_modules, n_hubs)
    
    # Filter nodes by module if specified
    if selected_modules is not None and module_colors is not None:
        nodes_to_keep = [i for i, color in enumerate(module_colors) 
                        if color in selected_modules]
    else:
        nodes_to_keep = range(len(tom_array))
    
    # Add edges that meet the threshold with weights
    for i in nodes_to_keep:
        for j in nodes_to_keep:
            if i < j and tom_array[i,j] > threshold:
                G.add_edge(i, j, weight=1/tom_array[i,j], tom_value=tom_array[i,j])
    
    # Process components
    components = list(nx.connected_components(G))
    if components:
        significant_components = [comp for comp in components 
                               if len(comp) >= minimum_component_size]
        if significant_components:
            largest_component = max(significant_components, key=len)
            G = G.subgraph(largest_component).copy()
            st.write(f"Original network size: {len(nodes_to_keep)} nodes")
            st.write(f"Largest component size: {len(largest_component)} nodes")
        else:
            st.warning(f"No components with size >= {minimum_component_size} nodes found.")
            return None
    
    # Calculate layout using TOM-based weights
    pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), 
                          iterations=50, 
                          weight='weight',
                          seed=42)
    
    # Scale positions to SVG size
    scale = 700  # Slightly smaller than 800 to leave margin
    margin = 50
    
    # Find position ranges
    pos_array = np.array(list(pos.values()))
    x_min, y_min = pos_array.min(axis=0)
    x_max, y_max = pos_array.max(axis=0)
    
    def scale_pos(pos):
        x = (pos[0] - x_min) / (x_max - x_min) * scale + margin
        y = (pos[1] - y_min) / (y_max - y_min) * scale + margin
        return x, y

    # Create HTML content with white background
    html_content = f'''
    <div id="network-container" style="position: relative; width: 800px; height: 800px;">
        <div id="tooltip" style="display: none; position: absolute; background: white; 
             padding: 5px; border: 1px solid black; border-radius: 5px; pointer-events: none;"></div>
        <svg width="800" height="800" xmlns="http://www.w3.org/2000/svg">
        <rect width="800" height="800" fill="white"/>
    '''

    # Prepare nodes and edges like force-directed
    nodes = list(G.nodes())
    if module_colors is not None:
        colors = [get_color(module_colors[n]) for n in nodes]
    else:
        colors = ['skyblue'] * len(nodes)

    # Node classification
    regular_nodes = []
    important_nodes = []
    regular_colors = []
    important_colors = []
    regular_sizes = []
    important_sizes = []
    important_edges = []

    # Size settings (match force-directed)
    REGULAR_SIZE = 7     # 通常のノード
    HIGHLIGHT_SIZE = 12   # ハイライトされたノード
    HUB_SIZE = 16        # ハブノード

    # Categorize nodes
    for node in nodes:
        if node in hub_indices:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HUB_SIZE)
            important_edges.append(EDGE_COLORS['hub'])
        elif highlight_nodes and node in highlight_nodes:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HIGHLIGHT_SIZE)
            important_edges.append(EDGE_COLORS['highlight'])
        else:
            regular_nodes.append(node)
            regular_colors.append(colors[nodes.index(node)])
            regular_sizes.append(REGULAR_SIZE)

    # Draw edges with TOM-based width
    for (u, v, data) in G.edges(data=True):
        start = scale_pos(pos[u])
        end = scale_pos(pos[v])
        tom_value = data['tom_value']
        edge_width = tom_value * 5  # TOM値に基づいてエッジの太さを設定

        if u in hub_indices or v in hub_indices:
            edge_color = SVG_EDGE_COLORS['hub']  # ここを修正
        elif highlight_nodes and (u in highlight_nodes or v in highlight_nodes):
            edge_color = SVG_EDGE_COLORS['highlight']  # ここを修正
        else:
            edge_color = SVG_EDGE_COLORS['regular']  # ここを修正
        
        html_content += f'<line class="edge" x1="{start[0]}" y1="{start[1]}" x2="{end[0]}" y2="{end[1]}" ' \
                       f'stroke="{edge_color}" stroke-width="{edge_width}"/>'
    
    # Draw regular nodes first
    for node, color, size in zip(regular_nodes, regular_colors, regular_sizes):
        x, y = scale_pos(pos[node])
        html_content += f'''<circle class="node" cx="{x}" cy="{y}" r="{size/2}"
                           fill="{color}" stroke="none"
                           data-gene="{gene_names[node]}" data-type="Regular"
                           data-connections="{len(list(G.neighbors(node)))}"
                           data-tom-sum="{sum(tom_array[node, n] for n in G.neighbors(node)):.3f}"/>'''

    # Draw important nodes with border
    for node, color, size, edge_color in zip(important_nodes, important_colors, important_sizes, important_edges):
        x, y = scale_pos(pos[node])
        node_type = "Hub" if node in hub_indices else "Highlighted"
        html_content += f'''<circle class="node" cx="{x}" cy="{y}" r="{size/2}"
                           fill="{color}" stroke="rgba{edge_color}" stroke-width="2"
                           data-gene="{gene_names[node]}" data-type="{node_type}"
                           data-connections="{len(list(G.neighbors(node)))}"
                           data-tom-sum="{sum(tom_array[node, n] for n in G.neighbors(node)):.3f}"/>'''


    # Add labels for important nodes with white background
    for node in important_nodes:
        x, y = scale_pos(pos[node])
        # テキストの長さに基づいて背景の幅を調整
        text_length = len(gene_names[node])
        rect_width = text_length * 6  # 文字数に基づいて幅を調整
        rect_height = 14  # 高さを小さく
        html_content += f'''
            <g class="label-group">
                <rect x="{x-rect_width/2}" y="{y-rect_height/2}" 
                      width="{rect_width}" height="{rect_height}" 
                      fill="white" fill-opacity="0.7" rx="3"/>
                <text x="{x}" y="{y+4}" text-anchor="middle" 
                      font-family="Arial" font-size="10px">
                    {gene_names[node]}
                </text>
            </g>
        '''



    # Add title
    if hub_indices:
        connected_hubs = [i for i in hub_indices if i in nodes]
        if connected_hubs:
            hub_info = "\\nHub genes:\\n" + ", ".join([gene_names[i] for i in connected_hubs])
        else:
            hub_info = "\\nNo connected hub genes"
    else:
        hub_info = ""

    html_content += f'''
        <text x="400" y="30" text-anchor="middle" font-family="Arial" font-size="12px">
            Network Visualization (threshold > {threshold:.3f})
        </text>
        <text x="400" y="50" text-anchor="middle" font-family="Arial" font-size="12px">
            Showing modules: {", ".join(selected_modules) if selected_modules else "all"}
        </text>
        <text x="400" y="70" text-anchor="middle" font-family="Arial" font-size="12px">
            {len(G.nodes())} nodes, {len(G.edges())} edges{hub_info}
        </text>
    '''

    # Close SVG
    html_content += '</svg>'
    
    # Add JavaScript for hover functionality
    html_content += '''
        <script>
            document.addEventListener("DOMContentLoaded", function() {
                const tooltip = document.getElementById("tooltip");
                const nodes = document.getElementsByClassName("node");
                const container = document.getElementById("network-container");
                
                Array.from(nodes).forEach(node => {
                    node.addEventListener("mouseover", function(e) {
                        const gene = this.getAttribute("data-gene");
                        const type = this.getAttribute("data-type");
                        const connections = this.getAttribute("data-connections");
                        const tomSum = this.getAttribute("data-tom-sum");
                        tooltip.innerHTML = `Gene: ${gene}<br>` +
                                          `Type: ${type}<br>` +
                                          `Connections: ${connections}<br>` +
                                          `Total TOM: ${tomSum}`;
                        tooltip.style.display = "block";
                        tooltip.style.left = (e.clientX - container.getBoundingClientRect().left + 10) + "px";
                        tooltip.style.top = (e.clientY - container.getBoundingClientRect().top + 10) + "px";
                    });
                    
                    node.addEventListener("mousemove", function(e) {
                        tooltip.style.left = (e.clientX - container.getBoundingClientRect().left + 10) + "px";
                        tooltip.style.top = (e.clientY - container.getBoundingClientRect().top + 10) + "px";
                    });
                    
                    node.addEventListener("mouseout", function() {
                        tooltip.style.display = "none";
                    });
                });
            });
        </script>
    </div>
    '''
    
    return html_content

def plot_circular_layout(tom_array, gene_names, module_colors=None, selected_modules=None, threshold=0.1):
    """
    Create a Circular Layout visualization of the network
    
    Parameters:
    - tom_array: Topological Overlap Matrix
    - gene_names: List of gene names
    - module_colors: List of module colors for each gene
    - selected_modules: List of modules to visualize
    - threshold: Minimum TOM value to create an edge
    """
    # Create graph
    G = nx.Graph()
    
    # Add nodes and edges based on TOM threshold
    for i in range(len(tom_array)):
        for j in range(i+1, len(tom_array)):
            if tom_array[i,j] > threshold:
                G.add_edge(i, j, weight=tom_array[i,j])
    
    # Filter nodes by module if specified
    if selected_modules is not None and module_colors is not None:
        nodes_to_keep = [i for i, color in enumerate(module_colors) 
                         if color in selected_modules]
        G = G.subgraph(nodes_to_keep).copy()
    
    # Compute node centrality
    centrality = nx.eigenvector_centrality(G)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Compute circular layout
    pos = nx.circular_layout(G)
    
    # Draw edges first (with weight-based thickness)
    for (u, v, d) in G.edges(data=True):
        plt.plot([pos[u][0], pos[v][0]], 
                 [pos[u][1], pos[v][1]], 
                 color='lightgrey', 
                 alpha=0.3,
                 linewidth=d['weight']*3)
    
    # Draw nodes
    node_sizes = [centrality[node]*5000 for node in G.nodes()]
    
    # Color nodes based on module
    node_colors = []
    for node in G.nodes():
        if module_colors is not None and len(module_colors) > 0 and node < len(module_colors):
            node_colors.append(module_colors[node])
        else:
            node_colors.append('grey')
    
    nx.draw_networkx_nodes(G, pos, 
                            node_color=node_colors, 
                            node_size=node_sizes, 
                            alpha=0.7)
    
    # Add labels for top centrality nodes
    top_nodes = sorted(centrality, key=centrality.get, reverse=True)[:10]
    labels = {node: gene_names[node] for node in top_nodes}
    nx.draw_networkx_labels(G, pos, labels, font_size=8)
    
    plt.title('Circular Layout Network Visualization')
    plt.axis('off')
    plt.tight_layout()
    return fig


def get_color(module_color):
    """モジュール色名をmatplotlibの色に変換"""
    return WGCNA_COLORS.get(module_color, '#808080')  # 未知の色の場合はグレーを返す

@st.cache_data()
def calculate_tom_distribution(tom_array, module_colors=None, selected_modules=None):
    """Calculate TOM distribution data"""
    if selected_modules and module_colors is not None:
        # Get indices for selected modules
        selected_indices = [i for i, color in enumerate(module_colors) 
                          if color in selected_modules]
        # Get submatrix for selected modules
        sub_tom = tom_array[np.ix_(selected_indices, selected_indices)]
        # Get upper triangle values
        upper_tri = sub_tom[np.triu_indices_from(sub_tom, k=1)]
    else:
        upper_tri = tom_array[np.triu_indices_from(tom_array, k=1)]
    
    # Calculate statistics
    stats = {
        "values": upper_tri[upper_tri > 0],
        "percentiles": {
            "90th": np.percentile(upper_tri, 90),
            "95th": np.percentile(upper_tri, 95),
            "99th": np.percentile(upper_tri, 99)
        },
        "mean": np.mean(upper_tri)
    }
    return stats

def plot_tom_distribution_from_stats(stats, selected_modules=None):
    """Plot TOM distribution from pre-calculated statistics"""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot distribution
    sns.histplot(stats["values"], bins=50, ax=ax)
    ax.set_xlabel('TOM value')
    ax.set_ylabel('Count')
    ax.set_title('Distribution of TOM values' + 
                (f' for modules: {", ".join(selected_modules)}' if selected_modules else ''))
    
    # Add vertical lines for percentiles
    colors = ['r', 'g', 'b']
    percentiles = ['90th', '95th', '99th']
    for p, c in zip(percentiles, colors):
        val = stats["percentiles"][p]
        ax.axvline(x=val, color=c, linestyle='--', 
                  label=f'{p} percentile: {val:.3f}')
    
    ax.legend()
    return fig

@st.cache_data()
def get_hub_genes(tom_array, module_colors=None, selected_modules=None, n_hubs=10):
    """Get hub genes based on connectivity"""
    if selected_modules and module_colors is not None:
        # Get indices for selected modules
        node_indices = [i for i, color in enumerate(module_colors) 
                       if color in selected_modules]
        sub_tom = tom_array[np.ix_(node_indices, node_indices)]
        connectivity = np.sum(sub_tom, axis=1)
        if n_hubs == 0:
            return []
        else:
            # Get top n hub genes
            hub_indices = np.argsort(connectivity)[-n_hubs:]
            # Map back to original indices
            return [node_indices[i] for i in hub_indices]
    else:
        connectivity = np.sum(tom_array, axis=1)
        return np.argsort(connectivity)[-n_hubs:]

def plot_network(tom_array, gene_names, module_colors=None, threshold=0.1, 
               selected_modules=None, n_hubs=10, highlight_nodes=None,
               minimum_component_size=20):
    """Plot network visualization with hub genes and highlighted nodes"""
    if tom_array is None:
        return None
    
    # Create graph
    G = nx.Graph()
    
    # Get hub genes before filtering nodes
    hub_indices = get_hub_genes(tom_array, module_colors, selected_modules, n_hubs)
    
    # Filter nodes by module if specified
    if selected_modules is not None and module_colors is not None:
        nodes_to_keep = [i for i, color in enumerate(module_colors) 
                        if color in selected_modules]
    else:
        nodes_to_keep = range(len(tom_array))
    
    # Add edges that meet the threshold with weights
    for i in nodes_to_keep:
        for j in nodes_to_keep:
            if i < j and tom_array[i,j] > threshold:
                # TOM値の逆数を重みとして使用（高いTOM = 短い距離）
                G.add_edge(i, j, weight=1/tom_array[i,j], tom_value=tom_array[i,j])
    
    # Process components
    components = list(nx.connected_components(G))
    if components:
        significant_components = [comp for comp in components 
                               if len(comp) >= minimum_component_size]
        
        if significant_components:
            largest_component = max(significant_components, key=len)
            G = G.subgraph(largest_component).copy()
            
            st.write(f"Original network size: {len(nodes_to_keep)} nodes")
            st.write(f"Largest component size: {len(largest_component)} nodes")
        else:
            st.warning(f"No components with size >= {minimum_component_size} nodes found.")
            return None
    
    # Create plot with white background
    fig, ax = plt.subplots(figsize=(15, 15), facecolor='white')
    ax.set_facecolor('white')
    
    # Calculate layout using TOM-based weights
    pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), 
                          iterations=50, 
                          weight='weight',
                          seed=42)
    
    # Draw nodes
    nodes = list(G.nodes())
    if module_colors is not None:
        colors = [get_color(module_colors[n]) for n in nodes]
    else:
        colors = ['skyblue'] * len(nodes)

    # Node classification
    regular_nodes = []
    important_nodes = []
    regular_colors = []
    important_colors = []
    regular_sizes = []
    important_sizes = []
    important_edges = []

    # Size settings
    REGULAR_SIZE = 30     # 通常のノード
    HIGHLIGHT_SIZE = 45  # ハイライトされたノード
    HUB_SIZE = 60       # ハブノード



    # ノードの分類
    for node in nodes:
        if node in hub_indices:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HUB_SIZE)
            important_edges.append(EDGE_COLORS['hub'])
        elif highlight_nodes and node in highlight_nodes:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HIGHLIGHT_SIZE)
            important_edges.append(EDGE_COLORS['highlight'])
        else:
            regular_nodes.append(node)
            regular_colors.append(colors[nodes.index(node)])
            regular_sizes.append(REGULAR_SIZE)

    # エッジを分類して描画
    regular_edge_list = []
    hub_edge_list = []
    highlight_edge_list = []
    
    # エッジ幅も保存
    regular_widths = []
    hub_widths = []
    highlight_widths = []
    
    for (u, v, data) in G.edges(data=True):
        tom_value = data['tom_value']
        edge_width = tom_value * 5  # TOM値に基づいてエッジの太さを設定
        
        if u in hub_indices or v in hub_indices:
            hub_edge_list.append((u, v))
            hub_widths.append(edge_width)
        elif (highlight_nodes and 
             (u in highlight_nodes or v in highlight_nodes)):
            highlight_edge_list.append((u, v))
            highlight_widths.append(edge_width)
        else:
            regular_edge_list.append((u, v))
            regular_widths.append(edge_width)
   
    # Draw edges
    if regular_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=regular_edge_list,
                            width=regular_widths,
                            edge_color=[EDGE_COLORS['regular'][:3]] * len(regular_edge_list),  # RGBのみ
                            alpha=EDGE_COLORS['regular'][3])  # alphaを別途設定

    if hub_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=hub_edge_list,
                            width=hub_widths,
                            edge_color=[EDGE_COLORS['hub'][:3]] * len(hub_edge_list),
                            alpha=EDGE_COLORS['hub'][3])

    if highlight_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=highlight_edge_list,
                            width=highlight_widths,
                            edge_color=[EDGE_COLORS['highlight'][:3]] * len(highlight_edge_list),
                            alpha=EDGE_COLORS['highlight'][3])



    # Draw regular nodes first
    nx.draw_networkx_nodes(G, pos,
                        nodelist=regular_nodes,
                        node_color=regular_colors,
                        node_size=regular_sizes,
                        edgecolors='none',  # これは'none'のままで良い
                        linewidths=0)

    # Draw important nodes (hubs and highlighted) with border
    if important_nodes:
        # ノードごとにedgecolorsを設定
        edge_colors = []
        for node in important_nodes:
            if node in hub_indices:
                edge_colors.append(EDGE_COLORS['hub'])
            else:  # highlight node
                edge_colors.append(EDGE_COLORS['highlight'])
                
        nx.draw_networkx_nodes(G, pos,
                           nodelist=important_nodes,
                           node_color=important_colors,
                           node_size=important_sizes,
                           edgecolors=edge_colors,  # リストで色を指定
                           linewidths=2)
    # Add labels for important nodes
    texts = []
    for node in important_nodes:
        x, y = pos[node]
        texts.append(plt.text(x, y, gene_names[node],
                    bbox=dict(facecolor='white',
                           edgecolor='none',
                           alpha=0.7,
                           pad=0.5),
                    fontsize=10))

    # Adjust label positions
    adjust_text(texts,
               arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5),
               expand_points=(1.5, 1.5))
    
    # Title and hub genes info
    if hub_indices:
        connected_hubs = [i for i in hub_indices if i in nodes]
        if connected_hubs:
            hub_info = "\nHub genes:\n" + ", ".join([gene_names[i] for i in connected_hubs])
        else:
            hub_info = "\nNo connected hub genes"
    else:
        hub_info = ""
    
    plt.title(f'Network Visualization (threshold > {threshold:.3f})\n'
            f'Showing modules: {", ".join(selected_modules) if selected_modules else "all"}\n'
            f'{len(G.nodes())} nodes, {len(G.edges())} edges'
            f'{hub_info}')
    
    plt.axis('off')
    return fig


def plot_network_forceatlas2(tom_array, gene_names, module_colors=None, threshold=0.1, 
                             selected_modules=None, n_hubs=10, highlight_nodes=None,
                             minimum_component_size=20):
    """Plot network visualization using ForceAtlas2 layout"""
    if tom_array is None:
        return None
    
    # Create graph
    G = nx.Graph()
    
    # Get hub genes before filtering nodes
    hub_indices = get_hub_genes(tom_array, module_colors, selected_modules, n_hubs)
    
    # Filter nodes by module if specified
    if selected_modules is not None and module_colors is not None:
        nodes_to_keep = [i for i, color in enumerate(module_colors) 
                        if color in selected_modules]
    else:
        nodes_to_keep = range(len(tom_array))
    
    # Add edges that meet the threshold with weights
    for i in nodes_to_keep:
        for j in nodes_to_keep:
            if i < j and tom_array[i,j] > threshold:
                # ForceAtlas2では重みが大きいほど強い引力を表すため、TOM値をそのまま使用
                G.add_edge(i, j, weight=tom_array[i,j], tom_value=tom_array[i,j])
    
    # Process components
    components = list(nx.connected_components(G))
    if components:
        significant_components = [comp for comp in components 
                               if len(comp) >= minimum_component_size]
        
        if significant_components:
            largest_component = max(significant_components, key=len)
            G = G.subgraph(largest_component).copy()
            
            st.write(f"Original network size: {len(nodes_to_keep)} nodes")
            st.write(f"Largest component size: {len(largest_component)} nodes")
        else:
            st.warning(f"No components with size >= {minimum_component_size} nodes found.")
            return None
    
    # Create plot with white background
    fig, ax = plt.subplots(figsize=(15, 15), facecolor='white')
    ax.set_facecolor('white')
    
    # Calculate ForceAtlas2 layout
    if FA2_AVAILABLE:
        try:
            # Create ForceAtlas2 layout using ForceAtlas2Py package
            pos = FA2Py.forceatlas2_networkx_layout(
                G, 
                pos=None,
                iterations=100,
                # ForceAtlas2 parameters
                outboundAttractionDistribution=True,
                linLogMode=False,
                adjustSizes=False,
                edgeWeightInfluence=1.0,
                jitterTolerance=1.0,
                barnesHutOptimize=True,
                barnesHutTheta=1.2,
                multiThreaded=False,
                scalingRatio=2.0,
                strongGravityMode=False,
                gravity=1.0,
                verbose=False
            )
            st.success("✅ Using ForceAtlas2 algorithm")
            
        except Exception as e:
            st.warning(f"⚠️ ForceAtlas2 failed ({str(e)}), using spring layout approximation")
            pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), 
                                  iterations=100, seed=42)
                                  
    elif IGRAPH_AVAILABLE:
        try:
            nodes = list(G.nodes())
            edges = [(nodes.index(u), nodes.index(v)) for u, v in G.edges()]
            weights = [data['weight'] for u, v, data in G.edges(data=True)]
            
            # Create igraph graph
            ig_graph = ig.Graph(n=len(nodes), edges=edges)
            ig_graph.es['weight'] = weights
            
            # Use igraph's forceatlas2 layout
            layout = ig_graph.layout_forceatlas2(
                iterations=100,
                linlog=False,
                pos=None,
                nohubs=False,
                weight_attr='weight'
            )
            
            # Convert to networkx format
            pos = {nodes[i]: (layout[i][0], layout[i][1]) for i in range(len(nodes))}
            st.success("✅ Using ForceAtlas2 algorithm (igraph)")
            
        except Exception as e:
            st.warning(f"⚠️ igraph ForceAtlas2 failed ({str(e)}), using spring layout approximation")
            pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), 
                                  iterations=100, seed=42)
    else:
        # Use spring layout as ForceAtlas2 approximation
        st.info("📦 ForceAtlas2 package not available, using spring layout approximation")
        st.info("💡 Install with: pip install ForceAtlas2")
        pos = nx.spring_layout(G, k=1/np.sqrt(len(G.nodes())), 
                              iterations=100, seed=42)
    
    # Draw nodes (same as regular plot_network function)
    nodes = list(G.nodes())
    if module_colors is not None:
        colors = [get_color(module_colors[n]) for n in nodes]
    else:
        colors = ['skyblue'] * len(nodes)

    # Node classification
    regular_nodes = []
    important_nodes = []
    regular_colors = []
    important_colors = []
    regular_sizes = []
    important_sizes = []
    important_edges = []

    # Size settings
    REGULAR_SIZE = 30
    HIGHLIGHT_SIZE = 45
    HUB_SIZE = 60

    # Categorize nodes
    for node in nodes:
        if node in hub_indices:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HUB_SIZE)
            important_edges.append(EDGE_COLORS['hub'])
        elif highlight_nodes and node in highlight_nodes:
            important_nodes.append(node)
            important_colors.append(colors[nodes.index(node)])
            important_sizes.append(HIGHLIGHT_SIZE)
            important_edges.append(EDGE_COLORS['highlight'])
        else:
            regular_nodes.append(node)
            regular_colors.append(colors[nodes.index(node)])
            regular_sizes.append(REGULAR_SIZE)

    # Draw edges with classification (same as regular plot)
    regular_edge_list = []
    hub_edge_list = []
    highlight_edge_list = []
    
    regular_widths = []
    hub_widths = []
    highlight_widths = []
    
    for (u, v, data) in G.edges(data=True):
        tom_value = data['tom_value']
        edge_width = tom_value * 5
        
        if u in hub_indices or v in hub_indices:
            hub_edge_list.append((u, v))
            hub_widths.append(edge_width)
        elif (highlight_nodes and 
             (u in highlight_nodes or v in highlight_nodes)):
            highlight_edge_list.append((u, v))
            highlight_widths.append(edge_width)
        else:
            regular_edge_list.append((u, v))
            regular_widths.append(edge_width)

    # Draw edges
    if regular_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=regular_edge_list,
                            width=regular_widths,
                            edge_color=[EDGE_COLORS['regular'][:3]] * len(regular_edge_list),
                            alpha=EDGE_COLORS['regular'][3])

    if hub_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=hub_edge_list,
                            width=hub_widths,
                            edge_color=[EDGE_COLORS['hub'][:3]] * len(hub_edge_list),
                            alpha=EDGE_COLORS['hub'][3])

    if highlight_edge_list:
        nx.draw_networkx_edges(G, pos,
                            edgelist=highlight_edge_list,
                            width=highlight_widths,
                            edge_color=[EDGE_COLORS['highlight'][:3]] * len(highlight_edge_list),
                            alpha=EDGE_COLORS['highlight'][3])

    # Draw regular nodes first
    nx.draw_networkx_nodes(G, pos,
                        nodelist=regular_nodes,
                        node_color=regular_colors,
                        node_size=regular_sizes,
                        edgecolors='none',
                        linewidths=0)

    # Draw important nodes with border
    if important_nodes:
        edge_colors = []
        for node in important_nodes:
            if node in hub_indices:
                edge_colors.append(EDGE_COLORS['hub'])
            else:
                edge_colors.append(EDGE_COLORS['highlight'])
                
        nx.draw_networkx_nodes(G, pos,
                           nodelist=important_nodes,
                           node_color=important_colors,
                           node_size=important_sizes,
                           edgecolors=edge_colors,
                           linewidths=2)

    # Add labels for important nodes
    texts = []
    for node in important_nodes:
        x, y = pos[node]
        texts.append(plt.text(x, y, gene_names[node],
                    bbox=dict(facecolor='white',
                           edgecolor='none',
                           alpha=0.7,
                           pad=0.5),
                    fontsize=10))

    # Adjust label positions
    if texts:  # Only adjust if there are labels
        adjust_text(texts,
                   arrowprops=dict(arrowstyle='-', color='gray', alpha=0.5),
                   expand_points=(1.5, 1.5))
    
    # Title and info
    layout_info = "ForceAtlas2" if FA2_AVAILABLE else "Spring Layout (ForceAtlas2 approximation)"
    
    if hub_indices:
        connected_hubs = [i for i in hub_indices if i in nodes]
        if connected_hubs:
            hub_info = "\nHub genes:\n" + ", ".join([gene_names[i] for i in connected_hubs])
        else:
            hub_info = "\nNo connected hub genes"
    else:
        hub_info = ""
    
    plt.title(f'Network Visualization - {layout_info} (threshold > {threshold:.3f})\n'
            f'Showing modules: {", ".join(selected_modules) if selected_modules else "all"}\n'
            f'{len(G.nodes())} nodes, {len(G.edges())} edges'
            f'{hub_info}')
    
    plt.axis('off')
    return fig


def main():
    st.title('WGCNA Network Visualization')
    st.markdown("#### Upload WGCNA object ---.p file")
    uploaded_file = st.file_uploader("Upload WGCNA object SCENIC---.p file", type=['p'], label_visibility = 'collapsed')
    
    if uploaded_file is not None:
        try:
            with tempfile.TemporaryDirectory() as tmp_dir:
                temp_path = Path(tmp_dir) / "temp.p"
                with open(temp_path, 'wb') as f:
                    f.write(uploaded_file.getvalue())
                
                with open(temp_path, 'rb') as f:
                    pywgcna_obj = pickle.load(f)
            
            if hasattr(pywgcna_obj, 'TOM'):
                st.success("File loaded successfully!")
                
                # Get TOM and gene names
                if isinstance(pywgcna_obj.TOM, pd.DataFrame):
                    tom_array = pywgcna_obj.TOM.values
                    gene_names = pywgcna_obj.TOM.index.tolist()
                else:
                    tom_array = pywgcna_obj.TOM
                    # 遺伝子名はdatExpr.var.indexから取得
                    gene_names = pywgcna_obj.datExpr.var.index.tolist()
                
                # Get module colors
                module_colors = None
                if hasattr(pywgcna_obj, 'datExpr'):
                    if hasattr(pywgcna_obj.datExpr, 'var'):
                        if 'moduleColors' in pywgcna_obj.datExpr.var:
                            module_colors = pywgcna_obj.datExpr.var['moduleColors'].values
                            # Print debug information
                            st.write("Module colors found:")
                            st.write(pd.Series(module_colors).value_counts())


                # Module selection
                if module_colors is not None:
                    unique_modules = sorted(set(module_colors))
                    st.subheader("Module Selection")
                    show_all = st.checkbox("Show all modules", value=True)
                    if not show_all:
                        selected_modules = st.multiselect(
                            "Select modules",
                            unique_modules,
                            default=[unique_modules[0]] if unique_modules else None
                        )
                    else:
                        selected_modules = unique_modules


                # TOM Distribution と Statistics を選択されたモジュールに基づいて計算
                # 計算部分（キャッシュされる）
                tom_stats = calculate_tom_distribution(tom_array, module_colors, 
                                                    selected_modules if not show_all else None)
                
                # TOM Distribution
                st.header("TOM Value Distribution")
                fig_dist = plot_tom_distribution_from_stats(tom_stats, 
                                                          selected_modules if not show_all else None)
                st.pyplot(fig_dist)
                
                # TOM Statistics
                st.subheader("TOM Matrix Statistics")
                stats_display = {
                    "Total nodes": tom_array.shape[0],
                    "Mean TOM value": f"{tom_stats['mean']:.4f}",
                    "90th percentile": f"{tom_stats['percentiles']['90th']:.4f}",
                    "95th percentile": f"{tom_stats['percentiles']['95th']:.4f}",
                    "99th percentile": f"{tom_stats['percentiles']['99th']:.4f}"
                }
                st.write(stats_display)
                

                                  # Hub genes and node selection
                
                # 選択されたモジュールの遺伝子のみを表示
                if selected_modules and not show_all:
                    # 選択されたモジュールの遺伝子インデックスを取得
                    selected_indices = [i for i, color in enumerate(module_colors) 
                                      if color in selected_modules]
                    # 選択可能な遺伝子名リストを作成
                    available_genes = [gene_names[i] for i in selected_indices]
                else:
                    # すべてのモジュールが選択されている場合は全遺伝子を表示
                    available_genes = gene_names


                # Settings form
                with st.form("network_settings"):
                    
                    if selected_modules and not show_all:
                        # 選択されたモジュールのTOMサブマトリックスを取得
                        selected_indices = [i for i, color in enumerate(module_colors) 
                                          if color in selected_modules]
                        sub_tom = tom_array[np.ix_(selected_indices, selected_indices)]
                        
                        # 選択されたモジュールのデータに基づいてthresholdを設定
                        threshold = st.number_input(
                            "TOM similarity threshold",
                            min_value=0.0000,
                            max_value=1.0000,
                            value=float(tom_stats['percentiles']['90th']),  # 選択モジュールの90パーセンタイル
                            format="%.4f",
                            step=0.0100
                        )
                    else:
                        # 全体のデータに基づいてthresholdを設定
                        threshold = st.number_input(
                            "TOM similarity threshold",
                            min_value=0.0000,
                            max_value=1.0000,
                            value=float(tom_stats['percentiles']['90th']),
                            format="%.4f",
                            step=0.0100
                        )




                    st.subheader("Node Highlighting")
                    n_hubs = st.number_input("Number of hub genes to show", 
                                       min_value=0, max_value=50, value=10)

                    # Gene selection for highlighting
                    if module_colors is not None:
                        # 選択されたモジュールの遺伝子のみを表示
                        if selected_modules and not show_all:
                            selected_indices = [i for i, color in enumerate(module_colors) 
                                              if color in selected_modules]
                            available_genes = [gene_names[i] for i in selected_indices]
                        else:
                            available_genes = gene_names
                        
                        selected_genes = st.multiselect(
                            "Select genes to highlight",
                            available_genes,
                            key="highlight_genes"
                        )
                        highlight_indices = [gene_names.index(gene) for gene in selected_genes]
                    else:
                        highlight_indices = None

                                                                # Visualization type selection
                    vis_type = st.selectbox(
                        "Select Network Visualization Type",
                        ["Force-Directed", "ForceAtlas2", "Interactive", "Circular Layout", "Hive Plot"],
                        key="vis_type_select",index=2
                    )
                    if not FA2_AVAILABLE and vis_type == "ForceAtlas2":
                        st.warning("📝 **ForceAtlas2 Note**: Using spring_layout approximation.")
                        st.info("💡 For true ForceAtlas2, install: `pip install ForceAtlas2` or `pip install python-igraph`")
                    st.write("Hive plot is extremely slow....")

                    submitted = st.form_submit_button("Update Network")


                if submitted:
                    st.header("Network Visualization")
                    
                    # プロット生成（キャッシュされる）
                    if vis_type == "Interactive":
                        html_content, fig_static = generate_network_plots(
                            tom_array, gene_names, module_colors, threshold,
                            selected_modules, n_hubs, highlight_indices,
                            minimum_component_size=20, plot_type=vis_type
                        )
                        
                        if html_content and fig_static:
                            st.components.v1.html(html_content, height=800, width=800)
                            
                            download_container = st.container()
                            with download_container:
                                st.download_button(
                                    "Download Interactive HTML",
                                    html_content,
                                    "interactive_network.html",
                                    "text/html"
                                )
                                
                                col1, col2 = st.columns(2)
                                # Save as PDF
                                pdf = io.BytesIO()
                                fig_static.savefig(pdf, format='pdf', bbox_inches='tight')
                                col1.download_button(
                                    "Download PDF",
                                    pdf.getvalue(),
                                    "network.pdf",
                                    "application/pdf"
                                )
                                
                                # Save as PNG
                                png = io.BytesIO()
                                fig_static.savefig(png, format='png', dpi=300, bbox_inches='tight')
                                col2.download_button(
                                    "Download PNG",
                                    png.getvalue(),
                                    "network.png",
                                    "image/png"
                                )
                    
                    else:  # その他のプロットタイプ
                        fig = generate_network_plots(
                            tom_array, gene_names, module_colors, threshold,
                            selected_modules, n_hubs, highlight_indices,
                            minimum_component_size=20, plot_type=vis_type
                        )
                        
                        if fig:
                            st.pyplot(fig)
                            col1, col2 = st.columns(2)
                            
                            # Save as PDF
                            pdf = io.BytesIO()
                            fig.savefig(pdf, format='pdf', bbox_inches='tight')
                            col1.download_button(
                                "Download PDF",
                                pdf.getvalue(),
                                f"{vis_type.lower().replace(' ', '_')}_network.pdf",
                                "application/pdf"
                            )
                            
                            # Save as PNG
                            png = io.BytesIO()
                            fig.savefig(png, format='png', dpi=300, bbox_inches='tight')
                            col2.download_button(
                                "Download PNG",
                                png.getvalue(),
                                f"{vis_type.lower().replace(' ', '_')}_network.png",
                                "image/png"
                            )


      
        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
            import traceback
            st.write("Error traceback:", traceback.format_exc())


if __name__ == "__main__":
    main()
