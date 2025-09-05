import streamlit as st
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from umap import UMAP
from sklearn.preprocessing import RobustScaler
import networkx as nx
import pickle
import io
import tempfile
import os
from PyWGCNA import PyWGCNA
from adjustText import adjust_text

st.sidebar.title("Options")

@st.cache_data
def read_pywgcna(uploaded_file):
    with tempfile.NamedTemporaryFile(delete=False, suffix='.p') as temp_file:
        temp_file_path = temp_file.name
        temp_file.write(uploaded_file.getbuffer())
    
    try:
        pyWGCNA_df = PyWGCNA.readWGCNA(temp_file_path)
    finally:
        os.unlink(temp_file_path)  # Delete the temporary file
    
    return pyWGCNA_df

@st.cache_data
def process_pywgcna_data(_pyWGCNA_df):
    module_colors = _pyWGCNA_df.datExpr.var.moduleColors
    module_counts = module_colors.value_counts().sort_index()
    
    module_info = pd.DataFrame({
        'Module': module_counts.index,
        'Color': module_counts.index,
        'Gene Count': module_counts.values
    })
    
    modules = ['All'] + sorted(module_colors.unique())
    
    return module_colors, module_info, modules

@st.cache_data
def get_hub_gene_edges(tom, hub_genes, percentile=99):
    edges = []
    for hub_gene in hub_genes:
        for connected_gene in hub_genes:
            if hub_gene != connected_gene:
                weight = tom.loc[hub_gene, connected_gene]
                edges.append((hub_gene, connected_gene, weight))
    
    # Sort edges by weight and select top percentile
    edges.sort(key=lambda x: x[2], reverse=True)
    num_edges = int(len(edges) * (100 - percentile) / 100)
    return edges[:num_edges]

@st.cache_data
def calculate_kme(tom_values, module_genes_indices):
    return np.mean(tom_values[:, module_genes_indices], axis=1)


@st.cache_data
def run_module_umap(_pyWGCNA_df, module_colors, n_neighbors=30, min_dist=0.5, random_state=42):
    tom = _pyWGCNA_df.TOM
    if isinstance(tom, np.ndarray):
        tom = pd.DataFrame(tom, 
                           index=_pyWGCNA_df.datExpr.var.index, 
                           columns=_pyWGCNA_df.datExpr.var.index)
    # Calculate kME
    if not hasattr(_pyWGCNA_df, 'kME'):
        kME = pd.DataFrame(index=tom.index, columns=module_colors.unique())
        for color in module_colors.unique():
            module_genes = module_colors[module_colors == color].index
            module_genes_indices = [tom.index.get_loc(gene) for gene in module_genes]
            kME[color] = np.mean(tom.values[:, module_genes_indices], axis=1)
    else:
        kME = _pyWGCNA_df.kME
    
    # Run UMAP
    umap_model = UMAP(n_neighbors=n_neighbors, min_dist=min_dist, random_state=random_state)
    umap_embedding = umap_model.fit_transform(RobustScaler().fit_transform(tom))
    
    # Create DataFrame with UMAP results
    umap_df = pd.DataFrame(umap_embedding, columns=['UMAP1', 'UMAP2'], index=tom.index)
    umap_df['color'] = module_colors
    umap_df['kME'] = kME.max(axis=1)  # Use max kME across all modules
    
    return umap_df, tom


def plot_module_umap_network(umap_df, tom, point_size=20, size_scale=1, min_point_size=5,
                             edge_alpha=0.25, label_hubs=2, keep_grey_edges=False, edge_percentile=99,
                             min_edge_width=0.1, max_edge_width=2, show_edges=True, fig_size=(20, 16),
                             gene_label_font_size=8, selected_module=None):
    fig = plt.figure(figsize=fig_size)
    
    gs = fig.add_gridspec(1, 2, width_ratios=[3, 1])
    ax = fig.add_subplot(gs[0, 0])
    legend_ax = fig.add_subplot(gs[0, 1])

    if selected_module and selected_module != 'All':
        umap_df = umap_df[umap_df['color'] == selected_module]

    sizes = np.maximum(umap_df['kME'] * size_scale * point_size, min_point_size)

    hub_genes = []
    for module in umap_df['color'].unique():
        if module == 'grey' and not keep_grey_edges:
            continue
        module_df = umap_df[umap_df['color'] == module].sort_values('kME', ascending=False)
        hub_genes.extend(module_df.index[:label_hubs].tolist())

    if show_edges:
        edges = get_hub_gene_edges(tom, hub_genes, percentile=edge_percentile)
        G = nx.Graph()
        G.add_weighted_edges_from(edges)
        
        if G.number_of_edges() > 0:
            weights = [G[u][v]['weight'] for u, v in G.edges()]
            min_weight, max_weight = min(weights), max(weights)
            edge_widths = [min_edge_width + (max_edge_width - min_edge_width) * 
                           (w - min_weight) / (max_weight - min_weight) for w in weights]
            pos = {gene: (umap_df.loc[gene, 'UMAP1'], umap_df.loc[gene, 'UMAP2']) for gene in G.nodes() if gene in umap_df.index}
            nx.draw_networkx_edges(G, pos, ax=ax, alpha=edge_alpha, width=edge_widths)
        else:
            st.warning("No edges to display for the selected module(s).")

    # Draw non-hub genes
    non_hub_genes = umap_df.index[~umap_df.index.isin(hub_genes)]
    scatter = ax.scatter(
        umap_df.loc[non_hub_genes, 'UMAP1'], 
        umap_df.loc[non_hub_genes, 'UMAP2'], 
        c=umap_df.loc[non_hub_genes, 'color'], 
        s=sizes[~umap_df.index.isin(hub_genes)],
        alpha=0.7
    )

    # Draw hub genes with hollow circles
    texts = []
    for gene in hub_genes:
        if gene in umap_df.index:
            row = umap_df.loc[gene]
            ax.scatter(row['UMAP1'], row['UMAP2'], 
                       s=sizes[umap_df.index.get_loc(gene)] * 1.5,  # Make hub genes slightly larger
                       facecolors='none', 
                       edgecolors=row['color'], 
                       linewidth=2, 
                       alpha=1)
            texts.append(ax.text(row['UMAP1'], row['UMAP2'], gene, 
                                 fontsize=gene_label_font_size, 
                                 fontstyle='italic', 
                                 fontweight='bold'))

    # Adjust text positions to avoid overlaps
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='->', color='black', lw=0.5))

    title = 'WGCNA Module UMAP Visualization'
    if selected_module and selected_module != 'All':
        title += f' - {selected_module} Module'
    elif show_edges:
        title += f'\nwith Top {100-edge_percentile}% Hub Gene Connections'
    ax.set_title(title, fontsize=16, pad=20)
    ax.set_xlabel('UMAP1', fontsize=12)
    ax.set_ylabel('UMAP2', fontsize=12)
    ax.tick_params(axis='both', which='major', labelsize=10)

    # Remove frame
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)

    # Set aspect ratio to 'equal'
    ax.set_aspect('equal', adjustable='datalim')

    unique_colors = sorted(umap_df['color'].unique())
    legend_elements = [plt.scatter([], [], c=color, label=color, s=30) for color in unique_colors]
    # Add a legend element for hub genes
    legend_elements.append(plt.scatter([], [], facecolors='none', edgecolors='black', s=45, label='Hub Gene'))
    legend = legend_ax.legend(handles=legend_elements, title='Modules', loc='center left', 
                     fontsize=8, title_fontsize=10, bbox_to_anchor=(0, 0.5))
    legend_ax.axis('off')

    # Remove legend frame
    legend.get_frame().set_linewidth(0.0)

    plt.tight_layout()
    return fig

def main():
    st.title('WGCNA Module UMAP Visualization')

    # Initialize session state
    if 'pyWGCNA_df' not in st.session_state:
        st.session_state.pyWGCNA_df = None
        st.session_state.module_colors = None
        st.session_state.module_info = None
        st.session_state.modules = None
    if 'uploaded_file' not in st.session_state:
        st.session_state.uploaded_file = None

    uploaded_file = st.file_uploader("Upload pyWGCNA object (.p file)", type="p")
    
    if uploaded_file is not None:
        # Load pyWGCNA_df only if it hasn't been loaded yet or if a new file is uploaded
        if st.session_state.pyWGCNA_df is None or uploaded_file != st.session_state.uploaded_file:
            try:
                st.session_state.pyWGCNA_df = read_pywgcna(uploaded_file)
                st.session_state.module_colors, st.session_state.module_info, st.session_state.modules = process_pywgcna_data(st.session_state.pyWGCNA_df)
                st.session_state.uploaded_file = uploaded_file
            except Exception as e:
                st.error(f"Error loading the file: {str(e)}")
                return

        st.sidebar.header('Visualization Parameters')
        show_edges = st.sidebar.checkbox('Show Edges', value=True)
        
        point_size = st.sidebar.number_input('Base Point Size', value=20, min_value=1)
        size_scale = st.sidebar.number_input('Size Scale', value=2.0, min_value=0.1)
        min_point_size = st.sidebar.number_input('Minimum Point Size', value=5, min_value=1)
        gene_label_font_size = st.sidebar.number_input('Gene Label Font Size', value=8, min_value=1)
        label_hubs = st.sidebar.number_input('Number of Hub Genes to Label per Module', value=2, min_value=1)
        
        if show_edges:
            edge_alpha = st.sidebar.number_input('Edge Alpha', value=0.25, min_value=0.0, max_value=1.0)
            keep_grey_edges = st.sidebar.checkbox('Keep Grey Edges', value=False)
            edge_percentile = st.sidebar.number_input('Edge Percentile', value=99.0, min_value=0.0, max_value=100.0)
            min_edge_width = st.sidebar.number_input('Minimum Edge Width', value=0.1, min_value=0.01)
            max_edge_width = st.sidebar.number_input('Maximum Edge Width', value=2.0, min_value=0.1)

        fig_width = st.sidebar.number_input('Figure Width', value=8, min_value=1)
        fig_height = st.sidebar.number_input('Figure Height', value=6, min_value=1)

        st.sidebar.subheader('Download Options')
        file_format = st.sidebar.radio('Select file format:', ('PNG', 'PDF'))

        if st.session_state.module_info is not None:
            st.subheader("Module Information")
            st.dataframe(st.session_state.module_info)

            selected_module = st.selectbox('Select Module to Display', st.session_state.modules)

            if st.button('Run'):
                umap_results, tom = run_module_umap(
                                                        st.session_state.pyWGCNA_df, 
                                                        st.session_state.module_colors, 
                                                        n_neighbors=30, 
                                                        min_dist=0.5
                                                    )

                fig = plot_module_umap_network(
                    umap_results, 
                    tom,
                    point_size=point_size, 
                    size_scale=size_scale, 
                    min_point_size=min_point_size,
                    edge_alpha=edge_alpha if show_edges else 0,
                    label_hubs=label_hubs,
                    keep_grey_edges=keep_grey_edges if show_edges else False,
                    edge_percentile=edge_percentile if show_edges else 100,
                    min_edge_width=min_edge_width if show_edges else 0,
                    max_edge_width=max_edge_width if show_edges else 0,
                    show_edges=show_edges,
                    fig_size=(fig_width, fig_height),
                    gene_label_font_size=gene_label_font_size,
                    selected_module=selected_module if selected_module != 'All' else None
                )

                st.pyplot(fig)

                buf = io.BytesIO()
                if file_format == 'PNG':
                    plt.savefig(buf, format='png', dpi=300, bbox_inches='tight')
                    file_extension = 'png'
                    mime = "image/png"
                else:
                    plt.savefig(buf, format='pdf', bbox_inches='tight')
                    file_extension = 'pdf'
                    mime = "application/pdf"
                buf.seek(0)
                
                st.download_button(
                    label="Download Figure",
                    data=buf,
                    file_name=f"wgcna_umap_visualization.{file_extension}",
                    mime=mime
                )
    else:
        st.info("Please upload a pyWGCNA object file to begin.")

if __name__ == '__main__':
    main()
