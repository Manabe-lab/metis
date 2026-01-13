"""
Dynamo Visualization
Visualize Dynamo analysis results with advanced vector field plots
"""

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
from helper_func import clear_old_directories, clear_old_files
from streamlit_sortables import sort_items

@st.cache_data
def create_cell_color_mapping(cell_list, palette_name):
    """
    Create a consistent mapping between cell names/clusters and colors

    Parameters
    ----------
    cell_list : list
        List of cell names/cluster names
    palette_name : str
        Name of the discrete color palette to use

    Returns
    -------
    dict
        Dictionary with cell names/cluster names as keys and colors as values
    """
    n_cells = len(cell_list)
    base_palette = sns.color_palette(palette_name)
    base_n = len(base_palette)

    if n_cells <= base_n:
        colors = base_palette[:n_cells]
    else:
        colors = sns.color_palette(palette_name, n_colors=n_cells)

    return {cell: color for cell, color in zip(cell_list, colors)}

import hashlib

def _get_file_hash(file_obj):
    """Get MD5 hash of file content for cache key"""
    content = file_obj.getvalue()
    return hashlib.md5(content).hexdigest()

@st.cache_resource
def load_h5ad_file(_file_obj, _file_hash):
    """Load h5ad file with caching using cache_resource (no deep copy)"""
    import io
    temp_path = os.path.join("temp", f"h5ad_vis_temp_{round(time.time())}.h5ad")
    with open(temp_path, "wb") as f:
        f.write(_file_obj.getvalue())
    adata = sc.read_h5ad(temp_path)
    os.remove(temp_path)
    return adata

@st.cache_data
def extract_metadata_from_h5ad(file_bytes):
    """
    Extract metadata from h5ad file with caching

    Returns a dictionary containing:
    - vecfld_keys: list of Vector Field keys
    - vecfld_bases: list of basis names from Vector Fields
    - available_embeddings: list of available embeddings
    - categorical_cols: list of categorical column names
    - continuous_cols: list of continuous column names
    - geometry_cols: list of geometry column names
    - geometry_by_basis: dict mapping basis to geometry columns
    - dynamics_genes: list of genes for dynamics
    """
    import tempfile

    # Load adata temporarily
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    adata = sc.read_h5ad(tmp_path)
    os.unlink(tmp_path)

    # Extract Vector Field information
    vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
    vecfld_bases = [k.replace('VecFld_', '') for k in vecfld_keys]

    # Get available embeddings
    available_embeddings = [key.replace('X_', '') for key in adata.obsm.keys() if key.startswith('X_')]

    # Get available color options
    categorical_cols = [col for col in adata.obs.columns
                       if adata.obs[col].dtype.name == 'category' or
                       adata.obs[col].nunique() < 50]
    continuous_cols = [col for col in adata.obs.columns
                      if adata.obs[col].dtype.name in ['float64', 'float32', 'int64', 'int32']]

    # Find geometry columns
    geometry_cols = [col for col in adata.obs.columns
                    if any(key in col for key in ['speed', 'divergence', 'curl', 'acceleration', 'curvature'])]

    # Group geometry columns by basis
    geometry_by_basis = {}
    feature_types = ['speed', 'divergence', 'curl', 'acceleration', 'curvature']

    for col in geometry_cols:
        matched_basis = None

        # Try to match against each Vector Field basis
        for basis_name in vecfld_bases:
            if f'_{basis_name}' in col:
                matched_basis = basis_name
                break
            elif col in feature_types and len(vecfld_bases) == 1:
                matched_basis = vecfld_bases[0]
                break

        # If still no match, try to find the feature base and extract suffix
        if not matched_basis:
            for feature_type in feature_types:
                if col.startswith(feature_type):
                    if col == feature_type:
                        if vecfld_bases:
                            matched_basis = vecfld_bases[0]
                    elif col.startswith(f'{feature_type}_'):
                        suffix = col[len(feature_type)+1:]
                        if suffix in vecfld_bases:
                            matched_basis = suffix
                    break

        # Add to appropriate group
        if matched_basis:
            if matched_basis not in geometry_by_basis:
                geometry_by_basis[matched_basis] = []
            geometry_by_basis[matched_basis].append(col)

    # Get genes for phase portraits
    if 'use_for_dynamics' in adata.var.columns:
        dynamics_genes = list(adata.var_names[adata.var.use_for_dynamics])
    else:
        dynamics_genes = list(adata.var_names[:100])

    return {
        'vecfld_keys': vecfld_keys,
        'vecfld_bases': vecfld_bases,
        'available_embeddings': available_embeddings,
        'categorical_cols': categorical_cols,
        'continuous_cols': continuous_cols,
        'geometry_cols': geometry_cols,
        'geometry_by_basis': geometry_by_basis,
        'dynamics_genes': dynamics_genes
    }

# Import dynamo
try:
    import dynamo as dyn
    DYNAMO_AVAILABLE = True
except ImportError:
    DYNAMO_AVAILABLE = False

st.set_page_config(page_title="Dynamo Visualization", page_icon="🌊", layout="wide")

st.title("🌊 Dynamo Visualization")

if not DYNAMO_AVAILABLE:
    st.error("""
    ❌ **Dynamo is not installed**

    Please install Dynamo using:
    ```bash
    pip install dynamo-release
    ```

    See: https://github.com/aristoteleo/dynamo-release
    """)
    st.stop()

st.markdown("""
Visualize Dynamo analysis results using various methods.

### Visualization Types
- **Streamline plot**: Represent vector field flow with smooth curves
- **Cell-wise vectors**: Display vectors as arrows for each cell
- **Grid vectors**: Vector field on a grid
- **Phase portraits**: Spliced/unspliced dynamics for each gene
- **Topography**: Potential landscape (energy landscape) - UMAP required
- **Geometric features panel**: Visualize speed, divergence, acceleration, curvature
  - **Multiple basis support**: If Geometry is computed for multiple bases, panels are automatically generated for all bases

### Support for Multiple Basis Analysis
When Vector Field/Geometry is computed for multiple bases (e.g., `rna.mnn.umap`, `rna.pca.umap`) in Dynamo Analysis:
- Vector Field visualization: Select each basis to visualize individually
- Geometric features panel: **Automatically generates panels for all bases** and displays them side by side (for comparison)

### References
- [Qiu et al. (2022) "Mapping transcriptomic vector fields of single cells" Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4)
- [Dynamo Documentation](https://dynamo-release.readthedocs.io/)
""")

# Initialize session state
if "dynamo_vis_temp_dir" not in st.session_state:
    dynamo_vis_temp_dir = os.path.join("temp", f"dynamo_vis_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(dynamo_vis_temp_dir, exist_ok=True)
    st.session_state.dynamo_vis_temp_dir = dynamo_vis_temp_dir
else:
    dynamo_vis_temp_dir = st.session_state.dynamo_vis_temp_dir

# ========================================
# Upload file
# ========================================
st.header("Step 1: Upload Dynamo result")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file (Dynamo result)",
    type=['h5ad'],
    key="dynamo_vis_h5ad_upload",
    help="h5ad file generated by Dynamo analysis app"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load data with caching (use hash for cache key to detect file changes)
    with st.spinner("Loading data..."):
        file_hash = _get_file_hash(uploaded_h5ad)
        adata = load_h5ad_file(uploaded_h5ad, file_hash)
        st.info(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Extract metadata with caching (speeds up repeated interactions)
    file_bytes = uploaded_h5ad.getvalue()
    metadata = extract_metadata_from_h5ad(file_bytes)

    # Unpack metadata
    vecfld_keys = metadata['vecfld_keys']
    vecfld_bases = metadata['vecfld_bases']
    available_embeddings = metadata['available_embeddings']
    categorical_cols = metadata['categorical_cols']
    continuous_cols = metadata['continuous_cols']
    geometry_cols = metadata['geometry_cols']
    geometry_by_basis = metadata['geometry_by_basis']
    dynamics_genes = metadata['dynamics_genes']

    # Check for required data
    has_velocity = 'velocity_S' in adata.layers or 'velocity' in adata.layers
    has_vector_field = len(vecfld_keys) > 0

    if not has_velocity:
        st.warning("⚠️ Velocity layer (velocity_S) not found. Some visualizations may not work.")

    if not has_vector_field:
        st.warning("⚠️ Vector field not found. Vector field visualizations will not be available.")
    else:
        with st.expander("📊 Vector Field Information", expanded=False):
            st.write(f"✓ Found vector fields: {', '.join(vecfld_keys)}")
            st.write(f"✓ Vector field available for bases: {', '.join(vecfld_bases)}")

    if len(available_embeddings) == 0:
        st.error("❌ No embeddings found in data. Cannot proceed with visualization.")
        st.stop()
    else:
        with st.expander("🗺️ Embedding Information", expanded=False):
            st.write(f"✓ Available embeddings: {', '.join(available_embeddings)}")

    # Display geometry information
    if geometry_cols:
        with st.expander("⚡ Geometry Features Information", expanded=False):
            st.write(f"✓ Found {len(geometry_cols)} geometry columns: {', '.join(geometry_cols[:10])}{'...' if len(geometry_cols) > 10 else ''}")
            if geometry_by_basis:
                st.write(f"✓ Detected geometry features for {len(geometry_by_basis)} bases: {', '.join(geometry_by_basis.keys())}")
            else:
                st.write("⚠️ Could not group geometry columns by basis")

    # ========================================
    # Visualization options
    # ========================================
    st.header("Step 2: Select visualization type")

    viz_type = st.selectbox(
        "Visualization type",
        [
            "Streamline plot",
            "Cell-wise vectors",
            "Grid vectors",
            "Phase portraits",
            "Topography (potential landscape)",
            "Geometric features panel",
            "Basic UMAP"
        ],
        help="Select visualization type"
    )

    # ========================================
    # Common parameters
    # ========================================
    st.subheader("Visualization parameters")

    with st.expander("Parameter Guide", expanded=False):
        st.markdown("""
        ### Color options
        - **Categorical**: Clusters, cell types, etc.
        - **Continuous**: Gene expression, pseudotime, etc.
        - **Geometric features**: speed, divergence, acceleration, curvature, etc.

        ### Coloring
        - **Discrete colormap**: For categorical variables (clusters, cell types, etc.)
        - **Continuous colormap**: For continuous variables (gene expression, pseudotime, etc.)

        ### Cluster Settings (Sidebar)
        - **Change cluster order**: Change cluster display order by drag and drop
        - Order changes are also reflected in color assignments

        ### Basis
        - Embedding used for visualization
        - **Visualizations requiring Vector Field** (Streamline, Cell-wise vectors, Grid vectors, Topography):
          - Only bases with computed Vector Field can be selected
          - If Vector Field is computed for multiple bases in Dynamo Analysis, all will be shown as options
        - **Geometric features panel**:
          - If Geometry is computed for multiple bases, **panels are automatically generated for all bases**
          - Independent panels are displayed for each basis (for comparison)
        - **Other visualizations** (Basic UMAP):
          - All embeddings present in the data can be selected
        - Common embeddings: umap, mnn_umap, harmony_umap, pca, tsne, draw_graph_fa, etc.
        - UMAP is recommended (if available)
        - **UMAP is required for Topography analysis**: Topography is difficult to analyze in high-dimensional space,
          so it only works accurately with UMAP-based vector fields

        ### Plot parameters
        - **figsize**: Figure size (width, height in inches)
        - **pointsize**: Size of data points
        - **alpha**: Transparency (0-1)
        - **quiver_length**: Vector length (cell-wise/grid vectors)
        - **quiver_size**: Vector thickness (cell-wise/grid vectors)
        """)

    # Sidebar for colormap selection
    with st.sidebar:
        st.markdown("### Visualization Options")

        colormap_discrete = st.selectbox(
            "Colormap (discrete):",
            ["tab10", "Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
             "tab20b", "tab20c", "Pastel1", "Pastel2", "Accent"],
            index=0,
            help="Color palette for categorical variables"
        )

        colormap_continuous = st.selectbox(
            "Colormap (continuous):",
            ["viridis", "plasma", "inferno", "magma", "cividis",
             "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "Blues", "Greens", "Greys"],
            index=0,
            help="Color palette for continuous variables"
        )

    col1, col2, col3 = st.columns(3)

    # Initialize default values
    color_col = None
    palette = None
    show_cluster_labels = False

    with col1:
        if viz_type in ["Streamline plot", "Cell-wise vectors", "Grid vectors", "Topography (potential landscape)", "Basic UMAP"]:
            # For Topography, add Potential/Pseudotime option
            if viz_type == "Topography (potential landscape)":
                color_type = st.radio("Color by", ["Categorical", "Continuous", "Geometric features", "Potential/Pseudotime (FP)"], index=0)
                # Store current color type in session state
                st.session_state['dynamo_color_type'] = color_type
            else:
                color_type = st.radio("Color by", ["Categorical", "Continuous", "Geometric features"], index=0)
                st.session_state['dynamo_color_type'] = color_type

            if color_type == "Categorical" and categorical_cols:
                default_idx = 0
                for i, col in enumerate(categorical_cols):
                    if any(key in col.lower() for key in ['cluster', 'cell_type', 'celltype', 'leiden', 'louvain']):
                        default_idx = i
                        break

                color_col = st.selectbox("Color column", categorical_cols, index=default_idx)

                # Cluster order management for categorical variables
                if color_col:
                    cluster_list = adata.obs[color_col].cat.categories.tolist() if adata.obs[color_col].dtype.name == 'category' else sorted(adata.obs[color_col].unique().tolist())

                    # Initialize sorted order if not exists
                    if 'dynamo_vis_sorted_order' not in st.session_state:
                        st.session_state.dynamo_vis_sorted_order = cluster_list.copy()
                        st.session_state.dynamo_vis_prev_color_col = color_col
                    elif st.session_state.get('dynamo_vis_prev_color_col') != color_col:
                        st.session_state.dynamo_vis_sorted_order = cluster_list.copy()
                        st.session_state.dynamo_vis_prev_color_col = color_col
                    else:
                        sorted_order = st.session_state.get('dynamo_vis_sorted_order')
                        # Ensure sorted_order includes all current clusters
                        missing_clusters = [c for c in cluster_list if c not in sorted_order]
                        if missing_clusters:
                            sorted_order = sorted_order + missing_clusters
                            st.session_state.dynamo_vis_sorted_order = sorted_order

                    # Color map management
                    if 'dynamo_vis_cluster_color_map' not in st.session_state or \
                       st.session_state.get('dynamo_vis_current_cmap', '') != colormap_discrete or \
                       st.session_state.get('dynamo_vis_prev_color_col') != color_col:
                        st.session_state.dynamo_vis_cluster_color_map = create_cell_color_mapping(
                            st.session_state.dynamo_vis_sorted_order, colormap_discrete
                        )
                        st.session_state.dynamo_vis_current_cmap = colormap_discrete
                    else:
                        # Ensure existing color map has all clusters
                        existing_map = st.session_state.dynamo_vis_cluster_color_map
                        missing_in_map = [c for c in cluster_list if c not in existing_map]
                        if missing_in_map:
                            st.session_state.dynamo_vis_cluster_color_map = create_cell_color_mapping(
                                st.session_state.dynamo_vis_sorted_order, colormap_discrete
                            )

                    # Cluster order sorting UI in sidebar
                    with st.sidebar:
                        st.markdown("---")
                        st.markdown("### Cluster Settings")

                        show_cluster_labels = st.checkbox(
                            "Show cluster labels on plot",
                            value=True,
                            help="Display labels on each cluster (when off, only legend is shown)"
                        )

                        sort_clusters = st.checkbox("Change cluster order?")
                        if sort_clusters:
                            with st.form("dynamo_cluster_sorter"):
                                sorted_order_new = sort_items(st.session_state.dynamo_vis_sorted_order.copy())
                                submitted_sort = st.form_submit_button("Done sorting")
                            if submitted_sort:
                                st.session_state.dynamo_vis_sorted_order = sorted_order_new
                                current_cmap = st.session_state.get('dynamo_vis_current_cmap', colormap_discrete)
                                st.session_state.dynamo_vis_cluster_color_map = create_cell_color_mapping(
                                    sorted_order_new, current_cmap
                                )
                                st.success("✓ Cluster order updated!")

                    palette = st.session_state.get('dynamo_vis_cluster_color_map')
                else:
                    palette = None

            elif color_type == "Continuous" and continuous_cols:
                color_col = st.selectbox("Color column", continuous_cols)
                palette = None  # Use continuous colormap instead
                show_cluster_labels = False  # Not applicable for continuous

            elif color_type == "Geometric features" and geometry_cols:
                color_col = st.selectbox("Geometric feature", geometry_cols)
                palette = None  # Use continuous colormap for geometry features
                show_cluster_labels = False  # Not applicable for geometric features

            elif color_type == "Potential/Pseudotime (FP)":
                # Check if potential_fp and pseudotime_fp exist
                has_potential_fp = 'potential_fp' in adata.obs.columns
                has_pseudotime_fp = 'pseudotime_fp' in adata.obs.columns

                topo_feature = st.radio(
                    "Select feature:",
                    ["potential_fp", "pseudotime_fp"],
                    index=0,
                    help="Compute and visualize potential and pseudotime using Fokker-Planck method"
                )

                # Store the selected feature for later computation
                st.session_state['selected_topo_feature'] = topo_feature
                color_col = topo_feature  # Set color_col for immediate use
                palette = None
                show_cluster_labels = False

                if not has_potential_fp or not has_pseudotime_fp:
                    st.info(f"Will compute {topo_feature} when 'Generate Visualization' is clicked")
                else:
                    st.success(f"{topo_feature} has already been computed")

            else:
                color_col = None
                palette = None
                show_cluster_labels = False
                st.warning("No suitable columns found for selected color type")

    with col2:
        # For vector field visualizations, show all embeddings if Vector Field exists
        if viz_type in ["Streamline plot", "Cell-wise vectors", "Grid vectors", "Topography (potential landscape)"]:
            if has_vector_field:
                # For Topography, ONLY show 2D embeddings with Vector Field
                if viz_type == "Topography (potential landscape)":
                    # Topography requires BOTH:
                    # 1. A 2D embedding
                    # 2. A Vector Field computed in that same 2D space
                    # Cannot use high-dimensional embeddings because Vector Field is also high-dimensional
                    available_bases_for_viz = []

                    for emb in available_embeddings:
                        vecfld_key = f'VecFld_{emb}'
                        if vecfld_key in adata.uns:
                            emb_key = f'X_{emb}'
                            if emb_key in adata.obsm:
                                n_dims = adata.obsm[emb_key].shape[1]
                                if n_dims == 2:  # ONLY 2D embeddings
                                    available_bases_for_viz.append(emb)

                    if not available_bases_for_viz:
                        st.error("No **2D embedding** with computed Vector Field found.")
                        st.warning(f"""
                        **Requirements for Topography analysis:**
                        - 2D embedding (UMAP, tSNE, etc.)
                        - Vector Field computed in that 2D space

                        **Current status:**
                        - Available Vector Fields: {', '.join(vecfld_bases) if vecfld_bases else 'None'}

                        **Solutions:**
                        1. Compute Vector Field using **2D embedding like UMAP** in Dynamo Analysis
                        2. Or use other visualization types (Streamline plot, etc.)

                        **Note:** Topography does not work with high-dimensional embeddings (PCA, etc.).
                        Vector Field is also computed in high-dimensional space, so it cannot be reduced to 2D.
                        """)
                        basis = None
                    else:
                        # Prefer 'umap' for 2D visualization
                        default_basis_idx = 0
                        umap_like = [b for b in available_bases_for_viz if 'umap' in b.lower()]
                        if umap_like:
                            default_basis_idx = available_bases_for_viz.index(umap_like[0])

                        basis = st.selectbox(
                            "Basis (2D embedding with Vector Field)",
                            available_bases_for_viz,
                            index=default_basis_idx,
                            help="For Topography analysis: 2D embedding + Vector Field required"
                        )

                        st.caption(f"✓ 2D embedding with Vector Field: {basis}")
                else:
                    # For other visualizations, show all embeddings
                    available_bases_for_viz = available_embeddings

                    # Prefer 'umap' for 2D visualization, else first available
                    default_basis_idx = 0
                    umap_like = [b for b in available_bases_for_viz if 'umap' in b.lower()]
                    if umap_like:
                        default_basis_idx = available_bases_for_viz.index(umap_like[0])

                    basis = st.selectbox(
                        "Basis (embedding for visualization)",
                        available_bases_for_viz,
                        index=default_basis_idx,
                        help="Embedding used for visualization. If Vector Field exists, any embedding can be used for visualization.\n\n"
                             f"Vector Field computed: {', '.join(vecfld_bases)}"
                    )

                    st.caption(f"ℹ️ Vector Field: {', '.join(vecfld_bases)}")
            else:
                st.error("Vector Field not found. Please compute it in Dynamo Analysis.")
                basis = None
        elif viz_type == "Geometric features panel":
            # Geometric features panel automatically generates panels for ALL computed bases
            # No basis selection needed
            basis = None  # Not used for this visualization type
            if geometry_by_basis:
                st.info(f"**Automatic panel generation mode**\n\nAutomatically generating panels for all computed bases ({len(geometry_by_basis)}):\n- {', '.join(geometry_by_basis.keys())}")
            else:
                st.warning("No bases with computed Geometry found")
        else:
            # For other visualizations (Basic UMAP), show all embeddings
            default_basis_idx = 0
            if 'umap' in available_embeddings:
                default_basis_idx = available_embeddings.index('umap')
            elif 'pca' in available_embeddings:
                default_basis_idx = available_embeddings.index('pca')

            basis = st.selectbox(
                "Basis",
                available_embeddings,
                index=default_basis_idx,
                help="Embedding used for visualization (only those present in data are shown)"
            )

        fig_width = st.slider("Figure width", 2, 20, 6)
        fig_height = st.slider("Figure height", 2, 20, 4)

    with col3:
        pointsize = st.slider("Point size", 0.01, 2.0, 0.5, 0.05)
        alpha = st.slider("Transparency (alpha)", 0.0, 1.0, 0.7, 0.05)

        if viz_type in ["Cell-wise vectors", "Grid vectors"]:
            quiver_length = st.slider("Quiver length", 1, 20, 6)
            quiver_size = st.slider("Quiver size", 1, 20, 6)

    # Type-specific parameters
    if viz_type == "Phase portraits":
        st.subheader("Phase portrait parameters")

        n_genes_select = st.slider("Number of genes to plot", 1, 20, 4)

        # Gene selection method
        gene_select_method = st.radio(
            "Gene selection",
            ["Use dynamics genes (auto)", "Manual selection"],
            index=0
        )

        if gene_select_method == "Manual selection":
            all_genes = list(adata.var_names)
            selected_genes = st.multiselect(
                "Select genes",
                all_genes,
                default=all_genes[:n_genes_select] if len(all_genes) >= n_genes_select else all_genes,
                help="Select genes for phase portraits"
            )
        else:
            selected_genes = dynamics_genes[:n_genes_select]

        st.info(f"Will plot phase portraits for: {', '.join(selected_genes[:10])}{'...' if len(selected_genes) > 10 else ''}")

    # ========================================
    # Generate visualization
    # ========================================
    if st.button("🎨 Generate Visualization", type="primary"):
        st.header("Step 3: Visualization")

        try:
            # Clear all previous matplotlib figures to avoid overlap issues
            plt.close('all')

            # Set Dynamo style
            dyn.configuration.set_figure_params('dynamo', background='white')

            if viz_type == "Streamline plot":
                st.subheader("Streamline Plot")

                with st.spinner("Generating streamline plot..."):
                    # Explicitly close all previous figures to avoid overlap
                    plt.close('all')
                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    if color_col:
                        # Determine legend display based on user preference
                        if palette and show_cluster_labels:
                            legend_setting = 'on data'
                        elif palette and not show_cluster_labels:
                            legend_setting = 'right margin'
                        else:
                            legend_setting = False

                        plot_kwargs = {
                            'adata': adata,
                            'color': [color_col],
                            'basis': basis,
                            'show_legend': legend_setting,
                            'show_arrowed_spines': False,
                            'pointsize': pointsize,
                            'alpha': alpha,
                            'ax': ax,
                            'save_show_or_return': 'return'
                        }
                        # Add palette for categorical or cmap for continuous
                        if palette:
                            plot_kwargs['palette'] = palette
                        else:
                            plot_kwargs['cmap'] = colormap_continuous

                        dyn.pl.streamline_plot(**plot_kwargs)
                    else:
                        dyn.pl.streamline_plot(
                            adata,
                            basis=basis,
                            show_arrowed_spines=False,
                            pointsize=pointsize,
                            alpha=alpha,
                            ax=ax,
                            save_show_or_return='return'
                        )

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_vis_temp_dir, "streamline_plot.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="streamline_plot.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        try:
                            fig_path_pdf = os.path.join(dynamo_vis_temp_dir, "streamline_plot.pdf")
                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                            with open(fig_path_pdf, "rb") as f:
                                st.download_button(
                                    "⬇️ Download PDF",
                                    f,
                                    file_name="streamline_plot.pdf",
                                    mime="application/pdf"
                                )
                        except ValueError as e:
                            st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                        except Exception as e:
                            st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

            elif viz_type == "Cell-wise vectors":
                st.subheader("Cell-wise Vectors")

                with st.spinner("Generating cell-wise vectors..."):
                    # Explicitly close all previous figures to avoid overlap
                    plt.close('all')
                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    if color_col:
                        # Determine legend display based on user preference
                        if palette and show_cluster_labels:
                            legend_setting = 'on data'
                        elif palette and not show_cluster_labels:
                            legend_setting = 'right margin'
                        else:
                            legend_setting = False

                        plot_kwargs = {
                            'adata': adata,
                            'color': [color_col],
                            'basis': basis,
                            'show_legend': legend_setting,
                            'quiver_length': quiver_length,
                            'quiver_size': quiver_size,
                            'pointsize': pointsize,
                            'alpha': alpha,
                            'show_arrowed_spines': False,
                            'ax': ax,
                            'save_show_or_return': 'return'
                        }
                        if palette:
                            plot_kwargs['palette'] = palette
                        else:
                            plot_kwargs['cmap'] = colormap_continuous

                        dyn.pl.cell_wise_vectors(**plot_kwargs)
                    else:
                        dyn.pl.cell_wise_vectors(
                            adata,
                            basis=basis,
                            quiver_length=quiver_length,
                            quiver_size=quiver_size,
                            pointsize=pointsize,
                            alpha=alpha,
                            show_arrowed_spines=False,
                            ax=ax,
                            save_show_or_return='return'
                        )

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_vis_temp_dir, "cell_wise_vectors.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="cell_wise_vectors.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        try:
                            fig_path_pdf = os.path.join(dynamo_vis_temp_dir, "cell_wise_vectors.pdf")
                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                            with open(fig_path_pdf, "rb") as f:
                                st.download_button(
                                    "⬇️ Download PDF",
                                    f,
                                    file_name="cell_wise_vectors.pdf",
                                    mime="application/pdf"
                                )
                        except ValueError as e:
                            st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                        except Exception as e:
                            st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

            elif viz_type == "Grid vectors":
                st.subheader("Grid Vectors")

                with st.spinner("Generating grid vectors..."):
                    # Explicitly close all previous figures to avoid overlap
                    plt.close('all')
                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    if color_col:
                        # Determine legend display based on user preference
                        if palette and show_cluster_labels:
                            legend_setting = 'on data'
                        elif palette and not show_cluster_labels:
                            legend_setting = 'right margin'
                        else:
                            legend_setting = False

                        plot_kwargs = {
                            'adata': adata,
                            'color': color_col,
                            'basis': basis,
                            'show_legend': legend_setting,
                            'quiver_length': quiver_length,
                            'quiver_size': quiver_size,
                            'ax': ax,
                            'save_show_or_return': 'return'
                        }
                        # Add palette for categorical or cmap for continuous
                        if palette:
                            plot_kwargs['palette'] = palette
                        else:
                            plot_kwargs['cmap'] = colormap_continuous

                        dyn.pl.grid_vectors(**plot_kwargs)
                    else:
                        dyn.pl.grid_vectors(
                            adata,
                            basis=basis,
                            quiver_length=quiver_length,
                            quiver_size=quiver_size,
                            ax=ax,
                            save_show_or_return='return'
                        )

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_vis_temp_dir, "grid_vectors.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="grid_vectors.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        try:
                            fig_path_pdf = os.path.join(dynamo_vis_temp_dir, "grid_vectors.pdf")
                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                            with open(fig_path_pdf, "rb") as f:
                                st.download_button(
                                    "⬇️ Download PDF",
                                    f,
                                    file_name="grid_vectors.pdf",
                                    mime="application/pdf"
                                )
                        except ValueError as e:
                            st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                        except Exception as e:
                            st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

            elif viz_type == "Phase portraits":
                st.subheader("Phase Portraits")

                with st.spinner(f"Generating phase portraits for {len(selected_genes)} genes..."):
                    # Explicitly close all previous figures to avoid overlap
                    plt.close('all')

                    # Calculate grid layout for display
                    n_genes = len(selected_genes)
                    n_cols = min(4, n_genes)
                    n_rows = (n_genes + n_cols - 1) // n_cols

                    # Dynamo's phase_portraits creates its own figure, so we let it handle the layout
                    if color_col:
                        fig = dyn.pl.phase_portraits(
                            adata,
                            genes=selected_genes,
                            color=color_col,
                            figsize=(fig_width * n_cols, fig_height * n_rows),
                            ncols=n_cols,
                            save_show_or_return='return'
                        )
                    else:
                        fig = dyn.pl.phase_portraits(
                            adata,
                            genes=selected_genes,
                            figsize=(fig_width * n_cols, fig_height * n_rows),
                            ncols=n_cols,
                            save_show_or_return='return'
                        )

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_vis_temp_dir, "phase_portraits.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="phase_portraits.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        try:
                            fig_path_pdf = os.path.join(dynamo_vis_temp_dir, "phase_portraits.pdf")
                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                            with open(fig_path_pdf, "rb") as f:
                                st.download_button(
                                    "⬇️ Download PDF",
                                    f,
                                    file_name="phase_portraits.pdf",
                                    mime="application/pdf"
                                )
                        except ValueError as e:
                            st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                        except Exception as e:
                            st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

            elif viz_type == "Topography (potential landscape)":
                st.subheader("Topography / Potential Landscape")

                # basis is guaranteed to have Vector Field (filtered in selection)
                if basis is None:
                    st.warning("No 2D embedding with computed Vector Field found. Please compute it in Dynamo Analysis.")
                else:
                    st.info(f"✓ Using: **{basis}** (Vector Field: VecFld_{basis})")

                    # Show embedding info
                    embedding_key = f'X_{basis}'
                    if embedding_key in adata.obsm:
                        emb_shape = adata.obsm[embedding_key].shape
                        st.caption(f"📍 Embedding shape: {emb_shape}")

                    # Check if potential_fp/pseudotime_fp was selected and needs computation
                    # ONLY process if Color type is "Potential/Pseudotime (FP)"
                    selected_topo_feature = st.session_state.get('selected_topo_feature', None)
                    current_color_type = st.session_state.get('dynamo_color_type', None)

                    if current_color_type == "Potential/Pseudotime (FP)" and selected_topo_feature and selected_topo_feature != "None":
                        has_potential_fp = 'potential_fp' in adata.obs.columns
                        has_pseudotime_fp = 'pseudotime_fp' in adata.obs.columns
                        has_fp_transition_rate = 'fp_transition_rate' in adata.obsp

                        # Check if both are already computed
                        if has_potential_fp and has_pseudotime_fp:
                            st.success(f"{selected_topo_feature} has already been computed. Using existing values.")
                            color_col = selected_topo_feature
                        # Compute if not already computed
                        elif not has_potential_fp or not has_pseudotime_fp:
                            with st.spinner(f"Computing {selected_topo_feature} using Fokker-Planck method..."):
                                try:
                                    # Compute velocities using Dynamo's standard function
                                    st.caption("Computing cell velocities for transition matrix...")

                                    # Use Dynamo's cell_velocities to create transition matrix
                                    dyn.tl.cell_velocities(adata, basis=basis)

                                    st.caption("Cell velocities computed")

                                    # Now compute potential and divergence using ddhodge
                                    st.caption("Computing potential using Hodge decomposition...")

                                    from dynamo.ext import ddhodge

                                    ddhodge(
                                        adata,
                                        basis=basis,
                                        adjmethod='graphize_vecfld',
                                        enforce=True
                                    )

                                    st.caption("Hodge decomposition completed")

                                    # Check what columns were added by ddhodge
                                    potential_cols = [col for col in adata.obs.columns if 'potential' in col.lower()]
                                    divergence_cols = [col for col in adata.obs.columns if 'divergence' in col.lower() or 'div' in col.lower()]

                                    st.caption(f"Columns with 'potential': {potential_cols}")
                                    st.caption(f"Columns with 'divergence': {divergence_cols}")

                                    # Try different possible column names
                                    potential_col = None
                                    if 'potential' in adata.obs.columns:
                                        potential_col = 'potential'
                                    elif f'potential_{basis}' in adata.obs.columns:
                                        potential_col = f'potential_{basis}'
                                    elif potential_cols:
                                        potential_col = potential_cols[0]

                                    if potential_col:
                                        adata.obs["potential_fp"] = adata.obs[potential_col]
                                        adata.obs["pseudotime_fp"] = -adata.obs[potential_col]
                                        st.caption(f"Using '{potential_col}' for potential_fp and pseudotime_fp")
                                    else:
                                        st.error(f"""
                                        ❌ ddhodge did not generate 'potential' column

                                        Available obs columns: {list(adata.obs.columns)[-20:]}
                                        """)
                                        st.stop()

                                    st.success(f"✅ {selected_topo_feature} computed successfully!")

                                except ImportError as e:
                                    st.error(f"""
                                    ❌ **Import Error**

                                    Required Dynamo modules not found: {str(e)}

                                    Please ensure you have the latest version of dynamo-release installed.
                                    """)
                                    st.stop()
                                except Exception as e:
                                    st.error(f"""
                                    ❌ **Computation Error**

                                    {str(e)}

                                    Please check that your data contains valid Vector Field.
                                    """)
                                    st.exception(e)
                                    st.stop()

                        # Override color_col with the selected feature ONLY if FP was selected
                        if current_color_type == "Potential/Pseudotime (FP)":
                            color_col = selected_topo_feature
                            st.info(f"✓ Using {selected_topo_feature} for coloring")

                if basis is not None:
                    with st.spinner("Generating topography plot..."):
                        # basis is guaranteed to be 2D (filtered in selection above)
                        # Explicitly close all previous figures to avoid overlap
                        plt.close('all')
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        # Try with frontier=True first, fallback to frontier=False if it fails
                        try_frontier = True
                        topography_success = False
                        retry_count = 0
                        max_retries = 2  # Try with frontier=True, then frontier=False

                        while not topography_success and retry_count < max_retries:
                            try:
                                if color_col:
                                    # Determine legend display based on user preference
                                    if palette and show_cluster_labels:
                                        legend_setting = 'on data'
                                    elif palette and not show_cluster_labels:
                                        legend_setting = 'right margin'
                                    else:
                                        legend_setting = False

                                    plot_kwargs = {
                                        'adata': adata,
                                        'basis': basis,  # Visualization embedding and Vector Field basis must match
                                        'fps_basis': basis,  # Same as visualization basis
                                        'background': 'white',
                                        'color': [color_col],  # Wrap single color in list
                                        'streamline_color': 'black',
                                        'show_legend': legend_setting,
                                        'frontier': try_frontier,
                                        'pointsize': pointsize,
                                        'alpha': alpha,
                                        'ax': ax,
                                        'save_show_or_return': 'return'
                                    }
                                    # Add palette for categorical or cmap for continuous
                                    if palette:
                                        plot_kwargs['palette'] = palette
                                    else:
                                        plot_kwargs['cmap'] = colormap_continuous

                                    dyn.pl.topography(**plot_kwargs)
                                else:
                                    dyn.pl.topography(
                                        adata,
                                        basis=basis,  # Visualization embedding and Vector Field basis must match
                                        fps_basis=basis,  # Same as visualization basis
                                        background='white',
                                        streamline_color='black',
                                        frontier=try_frontier,
                                        pointsize=pointsize,
                                        alpha=alpha,
                                        ax=ax,
                                        save_show_or_return='return'
                                    )
                                topography_success = True
                                if not try_frontier:
                                    st.info("Fixed points analysis was skipped")
                            except (TypeError, AttributeError) as e:
                                retry_count += 1
                                if try_frontier and "NoneType" in str(e) and retry_count < max_retries:
                                    # Retry without frontier (fixed points)
                                    st.warning("Error occurred in fixed points analysis. Retrying with frontier=False...")
                                    try_frontier = False
                                    plt.close('all')
                                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
                                else:
                                    # Failed even with frontier=False
                                    st.error(f"""
                                    Topography analysis failed (basis: {basis})

                                    **Error:** {str(e)}

                                    **Possible causes:**
                                    - There may be issues with Vector Field computation
                                    - Data dimensions may be inappropriate

                                    **Solutions:**
                                    1. Try other visualization types (Streamline plot, etc.)
                                    2. Recompute Vector Field in Dynamo Analysis
                                    """)
                                    topography_success = False  # Skip topography generation
                                    break

                        # Only generate outputs if topography was successful
                        if topography_success:
                            # Generate PNG and PDF bytes (save to session_state to prevent re-rendering on download)
                            from matplotlib.backends.backend_agg import FigureCanvasAgg
                            from matplotlib.backends.backend_pdf import FigureCanvasPdf
                            import io

                            # PNG
                            canvas_png = FigureCanvasAgg(fig)
                            buf_png = io.BytesIO()
                            canvas_png.print_figure(buf_png, dpi=300, bbox_inches='tight', pad_inches=0.1, format='png')
                            buf_png.seek(0)
                            st.session_state['topography_png'] = buf_png.getvalue()

                            # PDF - conditional settings based on color type
                            try:
                                buf_pdf = io.BytesIO()

                                # For FP features, use default bbox to preserve layer alignment
                                # For other color types, use tight bbox for cleaner output
                                if current_color_type == "Potential/Pseudotime (FP)":
                                    # Don't use bbox_inches='tight' as it breaks the layer alignment for FP
                                    fig.savefig(buf_pdf, format='pdf', dpi=300)
                                else:
                                    # Use tight bbox for other color types
                                    fig.savefig(buf_pdf, format='pdf', bbox_inches='tight', pad_inches=0.1, dpi=300)

                                buf_pdf.seek(0)
                                st.session_state['topography_pdf'] = buf_pdf.getvalue()
                                st.session_state['topography_pdf_error'] = None
                            except Exception as e:
                                st.session_state['topography_pdf'] = None
                                st.session_state['topography_pdf_error'] = str(e)

                            st.pyplot(fig)

                            # Save as PNG and PDF (use cached bytes from session_state)
                            col_dl1, col_dl2 = st.columns(2)

                            with col_dl1:
                                if 'topography_png' in st.session_state:
                                    st.download_button(
                                        "⬇️ Download PNG",
                                        st.session_state['topography_png'],
                                        file_name="topography.png",
                                        mime="image/png"
                                    )
                                else:
                                    st.error("PNG save failed")
                                    st.info("The figure is displayed in the browser. Please use a screenshot.")

                            with col_dl2:
                                if st.session_state.get('topography_pdf') is not None:
                                    st.download_button(
                                        "⬇️ Download PDF",
                                        st.session_state['topography_pdf'],
                                        file_name="topography.pdf",
                                        mime="application/pdf"
                                    )
                                elif st.session_state.get('topography_pdf_error') is not None:
                                    st.warning(f"PDF save failed: {st.session_state['topography_pdf_error']}\n\nPlease use the PNG version.")
                                else:
                                    st.warning("PDF save failed\n\nPlease use the PNG version.")

            elif viz_type == "Geometric features panel":
                st.subheader("Geometric Features Panel")

                if not geometry_cols:
                    st.error("❌ No geometric features found in data. Please run Dynamo Analysis with 'Compute geometric features' enabled.")
                elif not geometry_by_basis:
                    st.error("❌ Could not detect basis information for geometric features.")
                else:
                    # Display info about available bases
                    st.info(f"✓ Found geometric features for {len(geometry_by_basis)} basis: {', '.join(geometry_by_basis.keys())}")

                    # Generate panel for each basis
                    for basis_idx, (geom_basis, geom_cols_for_basis) in enumerate(geometry_by_basis.items()):
                        st.markdown(f"### Basis: **{geom_basis}**")

                        with st.spinner(f"Generating geometric features panel for {geom_basis}..."):
                            # Explicitly close all previous figures to avoid overlap
                            plt.close('all')

                            # Create 2x2 panel
                            fig, axes = plt.subplots(2, 2, figsize=(fig_width * 2, fig_height * 2), constrained_layout=True)

                            # Speed
                            speed_cols = [col for col in geom_cols_for_basis if 'speed' in col]
                            if speed_cols:
                                speed_col = speed_cols[0]
                                dyn.pl.cell_wise_vectors(
                                    adata,
                                    color=speed_col,
                                    basis=geom_basis,
                                    pointsize=pointsize,
                                    alpha=alpha,
                                    quiver_length=6,
                                    quiver_size=6,
                                    ax=axes[0, 0],
                                    save_show_or_return='return'
                                )
                                axes[0, 0].set_title(f'Speed ({geom_basis})', fontsize=14)

                            # Divergence
                            div_cols = [col for col in geom_cols_for_basis if 'divergence' in col]
                            if div_cols:
                                div_col = div_cols[0]
                                dyn.pl.grid_vectors(
                                    adata,
                                    color=div_col,
                                    basis=geom_basis,
                                    quiver_length=12,
                                    quiver_size=12,
                                    ax=axes[0, 1],
                                    save_show_or_return='return'
                                )
                                axes[0, 1].set_title(f'Divergence ({geom_basis})', fontsize=14)

                            # Acceleration
                            accel_cols = [col for col in geom_cols_for_basis if 'acceleration' in col]
                            if accel_cols:
                                accel_col = accel_cols[0]
                                dyn.pl.streamline_plot(
                                    adata,
                                    color=accel_col,
                                    basis=geom_basis,
                                    ax=axes[1, 0],
                                    save_show_or_return='return'
                                )
                                axes[1, 0].set_title(f'Acceleration ({geom_basis})', fontsize=14)

                            # Curvature
                            curv_cols = [col for col in geom_cols_for_basis if 'curvature' in col]
                            if curv_cols:
                                curv_col = curv_cols[0]
                                dyn.pl.streamline_plot(
                                    adata,
                                    color=curv_col,
                                    basis=geom_basis,
                                    ax=axes[1, 1],
                                    save_show_or_return='return'
                                )
                                axes[1, 1].set_title(f'Curvature ({geom_basis})', fontsize=14)

                            st.pyplot(fig)

                            # Save as PNG and PDF (with basis name in filename)
                            col_dl1, col_dl2 = st.columns(2)

                            with col_dl1:
                                fig_path_png = os.path.join(dynamo_vis_temp_dir, f"geometric_features_panel_{geom_basis}.png")
                                fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                                with open(fig_path_png, "rb") as f:
                                    st.download_button(
                                        f"⬇️ Download PNG ({geom_basis})",
                                        f,
                                        file_name=f"geometric_features_panel_{geom_basis}.png",
                                        mime="image/png",
                                        key=f"png_geom_{basis_idx}"
                                    )

                            with col_dl2:
                                try:
                                    fig_path_pdf = os.path.join(dynamo_vis_temp_dir, f"geometric_features_panel_{geom_basis}.pdf")
                                    fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                                    with open(fig_path_pdf, "rb") as f:
                                        st.download_button(
                                            f"⬇️ Download PDF ({geom_basis})",
                                            f,
                                            file_name=f"geometric_features_panel_{geom_basis}.pdf",
                                            mime="application/pdf",
                                            key=f"pdf_geom_{basis_idx}"
                                        )
                                except ValueError as e:
                                    st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                                except Exception as e:
                                    st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

                            # Add separator between bases
                            if basis_idx < len(geometry_by_basis) - 1:
                                st.markdown("---")

            elif viz_type == "Basic UMAP":
                st.subheader("Basic UMAP")

                with st.spinner("Generating UMAP plot..."):
                    # Explicitly close all previous figures to avoid overlap
                    plt.close('all')
                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    if color_col:
                        # Determine legend display based on user preference
                        if palette and show_cluster_labels:
                            legend_setting = 'on data'
                        elif palette and not show_cluster_labels:
                            legend_setting = 'right margin'
                        else:
                            legend_setting = False

                        dyn.pl.umap(
                            adata,
                            color=color_col,
                            pointsize=pointsize,
                            alpha=alpha,
                            show_legend=legend_setting,
                            ax=ax,
                            save_show_or_return='return'
                        )
                    else:
                        dyn.pl.umap(
                            adata,
                            pointsize=pointsize,
                            alpha=alpha,
                            ax=ax,
                            save_show_or_return='return'
                        )

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_vis_temp_dir, "umap.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="umap.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        try:
                            fig_path_pdf = os.path.join(dynamo_vis_temp_dir, "umap.pdf")
                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                            with open(fig_path_pdf, "rb") as f:
                                st.download_button(
                                    "⬇️ Download PDF",
                                    f,
                                    file_name="umap.pdf",
                                    mime="application/pdf"
                                )
                        except ValueError as e:
                            st.warning(f"PDF save failed (contains non-finite numbers): {str(e)}\n\nPlease use the PNG version.")
                        except Exception as e:
                            st.warning(f"PDF save failed: {str(e)}\n\nPlease use the PNG version.")

            st.success("Visualization generated successfully!")

        except Exception as e:
            st.error(f"❌ Error during visualization: {str(e)}")
            st.exception(e)

else:
    st.info("Upload an h5ad file generated by Dynamo analysis app to get started")
