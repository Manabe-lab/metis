"""
Dynamo LAP (Least Action Path) Analysis
Compute optimal transition paths between cell states using vector field
"""

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
from pathlib import Path
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

# Import dynamo
try:
    import dynamo as dyn
    DYNAMO_AVAILABLE = True
except ImportError:
    DYNAMO_AVAILABLE = False

st.set_page_config(page_title="Dynamo LAP", page_icon="🛤️", layout="wide")

st.title("🛤️ Least Action Path Analysis")

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
Least Action Path (LAP) analysis computes the optimal transition path between cell states
in the vector field. This is based on the principle of least action from physics, where
cells follow paths that minimize the "action" (a measure of the energetic cost of transition).

**Prerequisites**:
- ✅ Vector field computed
- ✅ Transition matrix computed
- ✅ Cell type annotations

**Reference**: [Dynamo Documentation](https://dynamo-release.readthedocs.io/)
""")

@st.cache_data
def load_h5ad_file(file_bytes, filename):
    """Load h5ad file with caching"""
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    adata = sc.read_h5ad(tmp_path)
    os.unlink(tmp_path)

    return adata

@st.cache_data
def load_tfs(species="auto"):
    """
    Load transcription factors list from local file or dynamo sample data
    Returns list of TF gene symbols (case insensitive for mouse compatibility)

    Parameters
    ----------
    species : str
        "human", "mouse", or "auto" (default). If auto, tries local file first.

    Returns
    -------
    list
        List of TF gene symbols (lowercase for case-insensitive matching)
    """
    # Try local file first (streamlit directory)
    local_tf_paths = [
        str(Path(__file__).resolve().parent.parent / "human_tfs.txt"),
        str(Path(__file__).resolve().parent.parent / "mouse_tfs.txt"),
        "./human_tfs.txt",
        "./mouse_tfs.txt",
    ]

    for tf_path in local_tf_paths:
        if os.path.exists(tf_path):
            try:
                tfs_df = pd.read_csv(tf_path, sep="\t")
                # Get gene symbols from 'Symbol' column or first column
                if "Symbol" in tfs_df.columns:
                    tfs_names = list(tfs_df["Symbol"])
                else:
                    tfs_names = list(tfs_df.iloc[:, 0])

                # Convert to lowercase for case-insensitive matching (important for mouse!)
                tfs_names_lower = [tf.lower() for tf in tfs_names if isinstance(tf, str)]
                st.info(f"✓ Loaded {len(tfs_names_lower)} TFs from {tf_path}")
                return tfs_names_lower
            except Exception as e:
                st.warning(f"Could not load TFs from {tf_path}: {str(e)}")
                continue

    # Fallback: try dynamo's built-in function
    if species == "human" or species == "auto":
        try:
            tfs_df = dyn.sample_data.human_tfs()
            tfs_names = list(tfs_df["Symbol"]) if "Symbol" in tfs_df.columns else list(tfs_df.iloc[:, 0])
            tfs_names_lower = [tf.lower() for tf in tfs_names if isinstance(tf, str)]
            st.info(f"✓ Loaded {len(tfs_names_lower)} TFs from dynamo (human)")
            return tfs_names_lower
        except Exception as e:
            st.warning(f"Could not load TFs from dynamo: {str(e)}")

    st.error("❌ Could not load TFs data from any source")
    return None

def save_lap_to_adata(adata, transition_graph, selected_cell_types, representative_cells,
                      lap_basis, vis_basis, group_key):
    """
    Save all LAP analysis results to AnnData object

    Parameters
    ----------
    adata : AnnData
        AnnData object to save results to
    transition_graph : dict
        Dictionary containing all transition results
    selected_cell_types : list
        List of cell types analyzed
    representative_cells : dict
        Dictionary mapping cell type to representative cell name
    lap_basis : str
        Basis used for LAP calculation
    vis_basis : str
        Basis used for visualization
    group_key : str
        Key in adata.obs for cell type annotations
    """

    # 1. Save all transition results
    # Helper function to convert data to h5ad-compatible types
    def convert_for_h5ad(obj):
        """Recursively convert objects to h5ad-compatible types"""
        if obj is None:
            return None
        elif isinstance(obj, pd.Index):
            # Convert pandas Index to list
            return obj.tolist()
        elif isinstance(obj, (np.ndarray,)):
            return obj.tolist()  # numpy arrays to lists
        elif isinstance(obj, dict):
            return {str(k): convert_for_h5ad(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            return [convert_for_h5ad(item) for item in obj]
        elif isinstance(obj, (np.integer, np.floating)):
            return float(obj)
        elif isinstance(obj, (int, float, str, bool)):
            return obj
        else:
            # For complex objects, try to convert to string or skip
            try:
                return str(obj)
            except:
                return None

    lap_transitions = {}
    for key, data in transition_graph.items():
        # Convert ranking DataFrame to dict for h5ad compatibility
        ranking_dict = None
        if data.get('ranking') is not None:
            ranking_dict = {
                'genes': [str(g) for g in data['ranking']['all'].tolist()],
                'scores': [float(s) for s in data['ranking']['all_values'].tolist()]
            }

        # Convert LAP_dict and LAP_vis_dict to h5ad-compatible formats
        lap_dict_converted = convert_for_h5ad(data.get('LAP_dict'))
        lap_vis_dict_converted = convert_for_h5ad(data.get('LAP_vis_dict'))

        lap_transitions[str(key)] = {
            'action': float(data['action']) if data['action'] is not None else None,
            'time': float(data['time']) if data['time'] is not None else None,
            'ranking': ranking_dict,
            'LAP_dict': lap_dict_converted,
            'LAP_vis_dict': lap_vis_dict_converted
        }

    adata.uns['LAP_transitions'] = lap_transitions

    # 2. Save Action and Time matrices
    action_matrix = pd.DataFrame(
        index=selected_cell_types,
        columns=selected_cell_types,
        dtype=float
    )
    time_matrix = pd.DataFrame(
        index=selected_cell_types,
        columns=selected_cell_types,
        dtype=float
    )

    for key, data in transition_graph.items():
        start, end = key.split('->')
        if data['action'] is not None:
            action_matrix.loc[start, end] = data['action']
        if data['time'] is not None:
            time_matrix.loc[start, end] = data['time']

    # Convert to dict for h5ad compatibility
    # Ensure all values are plain Python types, not pandas/numpy types
    adata.uns['LAP_action_matrix'] = {
        'data': [[float(v) if not pd.isna(v) else None for v in row]
                 for row in action_matrix.values],
        'index': [str(idx) for idx in action_matrix.index],
        'columns': [str(col) for col in action_matrix.columns]
    }
    adata.uns['LAP_time_matrix'] = {
        'data': [[float(v) if not pd.isna(v) else None for v in row]
                 for row in time_matrix.values],
        'index': [str(idx) for idx in time_matrix.index],
        'columns': [str(col) for col in time_matrix.columns]
    }

    # 3. Save representative cells info
    # Mark representative cells in adata.obs
    adata.obs['LAP_representative'] = 'No'
    for cell_type, cell_name in representative_cells.items():
        if cell_name in adata.obs_names:
            adata.obs.loc[cell_name, 'LAP_representative'] = cell_type

    # Also save as dict in uns (convert to plain strings)
    adata.uns['LAP_representative_cells'] = {str(k): str(v) for k, v in representative_cells.items()}

    # 4. Save metadata (ensure all types are h5ad-compatible)
    adata.uns['LAP_metadata'] = {
        'lap_basis': str(lap_basis),
        'vis_basis': str(vis_basis),
        'group_key': str(group_key),
        'cell_types': [str(ct) for ct in selected_cell_types],
        'n_transitions': int(len(transition_graph))
    }

    return adata

def check_prerequisites(adata, basis):
    """Check if preprocessing is complete"""
    issues = []

    # Check vector field
    vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
    if not vecfld_keys:
        issues.append("❌ Vector field not found (.uns['VecFld_*'])")

    # Check transition matrix
    transition_keys = [k for k in adata.obsp.keys() if 'transition_matrix' in k]
    if not transition_keys:
        issues.append("❌ No transition matrix found in .obsp")

    # Check embedding
    if f'X_{basis}' not in adata.obsm:
        issues.append(f"❌ Embedding 'X_{basis}' not found")

    return issues

# Initialize session state
if "dynamo_lap_temp_dir" not in st.session_state:
    dynamo_lap_temp_dir = os.path.join("temp", f"dynamo_lap_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(dynamo_lap_temp_dir, exist_ok=True)
    st.session_state.dynamo_lap_temp_dir = dynamo_lap_temp_dir
else:
    dynamo_lap_temp_dir = st.session_state.dynamo_lap_temp_dir

if "lap_computed" not in st.session_state:
    st.session_state.lap_computed = False

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload Dynamo result")

st.markdown("""
**Required data**:
- ✅ Vector field computed (`adata.uns['VecFld_*']`)
- ✅ Transition matrix computed (`adata.obsp`)
- ✅ Cell type annotations (`adata.obs`)
- ✅ Embedding (UMAP, etc.)

→ Use the h5ad file generated from **Dynamo Analysis** app
""")

uploaded_h5ad = st.file_uploader(
    "Upload Dynamo result (h5ad)",
    type=['h5ad'],
    key="dynamo_lap_h5ad_upload",
    help="h5ad file generated by Dynamo Analysis app"
)

if uploaded_h5ad is not None:
    # Load h5ad file
    adata = load_h5ad_file(uploaded_h5ad.getvalue(), uploaded_h5ad.name)
    st.success(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Check required data
    st.subheader("Data validation")

    col1, col2 = st.columns(2)

    with col1:
        # Check for Vector Field keys
        vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
        if len(vecfld_keys) > 0:
            st.success("✓ Vector field available")
            vecfld_bases = [k.replace('VecFld_', '') for k in vecfld_keys]
            st.info(f"Vector Fields: {', '.join(vecfld_bases)}")
        else:
            st.error("❌ Vector field not found")
            st.warning("Please run Dynamo Analysis first")

    with col2:
        # Check for embeddings
        embedding_keys = [k for k in adata.obsm.keys() if k.startswith('X_')]
        if len(embedding_keys) > 0:
            st.success(f"✓ Embeddings available")
            available_bases = [k.replace('X_', '') for k in embedding_keys]
            st.info(f"Bases: {', '.join(available_bases)}")
        else:
            st.error("❌ No embeddings found")

    # Check if analysis can proceed
    can_proceed = (len(vecfld_keys) > 0 and len(embedding_keys) > 0)

    if not can_proceed:
        st.error("""
        ❌ **Cannot proceed with LAP analysis**

        This file is missing required data from Dynamo analysis.

        Required data:
        - ✅ Vector field
        - ✅ Embedding (UMAP, etc.)
        - ✅ Transition matrix

        Please run vector field computation in the Dynamo Analysis app first.
        """)
        st.stop()

    # ========================================
    # Sidebar: Colormap selection
    # ========================================
    with st.sidebar:
        st.markdown("### Visualization Options")

        colormap_discrete = st.selectbox(
            "Colormap (Discrete):",
            ["tab10", "Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
             "tab20b", "tab20c", "Pastel1", "Pastel2", "Accent"],
            index=0,
            help="Color palette for categorical variables"
        )

        colormap_continuous = st.selectbox(
            "Colormap (Continuous):",
            ["viridis", "plasma", "inferno", "magma", "cividis",
             "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "Blues", "Greens", "Greys"],
            index=0,
            help="Color palette for continuous variables"
        )

    # ========================================
    # Step 2: Configuration
    # ========================================
    st.header("Step 2: LAP Configuration")

    col1, col2 = st.columns(2)

    with col1:
        # Extract basis names from Vector Fields
        all_vf_bases = [k.replace('VecFld_', '') for k in vecfld_keys]

        # Filter out UMAP/tSNE/etc for LAP calculation (high-dimensional spaces only)
        high_dim_keywords = ['umap', 'tsne', 'phate', 'trimap']
        vf_bases_for_lap = [b for b in all_vf_bases
                            if not any(kw in b.lower() for kw in high_dim_keywords)]

        if not vf_bases_for_lap:
            st.error("❌ No high-dimensional Vector Field found")
            st.warning("""
            **LAP calculation requires high-dimensional space:**
            - ✅ PCA, MNN, etc.
            - ❌ UMAP, tSNE, etc. (low-dimensional embeddings not supported)

            Please compute Vector Field with a high-dimensional basis in Dynamo Analysis.
            """)
            st.stop()

        st.info(f"""
        **Available high-dimensional bases for LAP:**
        {', '.join(vf_bases_for_lap)}
        """)

        # Basis selection (high-dimensional space only)
        lap_basis = st.selectbox(
            "Vector Field basis for LAP calculation:",
            vf_bases_for_lap,
            help="LAP calculation is performed in high-dimensional space (PCA, MNN, etc.)"
        )

        # Check prerequisites
        issues = check_prerequisites(adata, lap_basis)

        if issues:
            st.error("**⚠️ Prerequisites not met:**")
            for issue in issues:
                st.markdown(issue)
            st.info("""
            **Solution**: Please run the following in the Dynamo Analysis app:
            1. Vector Field computation
            2. Cell velocities computation
            3. Transition matrix computation
            """)
            st.stop()
        else:
            st.success("✅ Prerequisites met!")

    with col2:
        # Visualization basis selection (UMAP, etc.)
        st.info("""
        **Visualization embedding:**
        Used for results visualization
        """)

        vis_basis = st.selectbox(
            "Embedding basis for visualization:",
            available_bases,
            index=available_bases.index('umap') if 'umap' in available_bases else 0,
            help="Embedding used for results visualization (UMAP recommended)"
        )

        # Transition matrix key selection
        transition_keys = [k for k in adata.obsp.keys() if 'transition_matrix' in k]

        if not transition_keys:
            st.error("❌ No transition matrix found")
            st.stop()

        adj_key = st.selectbox(
            "Transition matrix key:",
            transition_keys,
            help="Key for the transition matrix"
        )

    # ========================================
    # Cell Type Selection
    # ========================================
    st.subheader("Cell Type Selection")

    # Group key selection
    group_options = [c for c in adata.obs.columns
                if adata.obs[c].dtype.name in ['category', 'object']]

    if not group_options:
        st.error("❌ No categorical columns found in .obs")
        st.stop()

    group_key = st.selectbox(
        "Group key (cell type annotation):",
        group_options,
        index=group_options.index('cell_type') if 'cell_type' in group_options else 0
    )

    cell_types = sorted(adata.obs[group_key].unique().tolist())

    # Streamline plot preview (using visualization basis)
    with st.expander("📊 Preview Vector Field", expanded=False):
        st.caption(f"Visualization embedding: {vis_basis}")
        if st.button("Generate Streamline Plot"):
            with st.spinner("Generating plot..."):
                plt.close('all')
                fig, ax = plt.subplots(figsize=(10, 8))

                try:
                    # Get cluster list and create color mapping
                    cluster_list = adata.obs[group_key].cat.categories.tolist() if adata.obs[group_key].dtype.name == 'category' else sorted(adata.obs[group_key].unique().tolist())
                    palette = create_cell_color_mapping(cluster_list, colormap_discrete)

                    dyn.pl.streamline_plot(
                        adata,
                        color=[group_key],
                        basis=vis_basis,
                        palette=palette,
                        ax=ax,
                        save_show_or_return="return"
                    )

                    st.pyplot(fig)
                    plt.close()
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.exception(e)

    # Cell type selection
    selected_cell_types = st.multiselect(
        "Choose at least 2 cell types:",
        cell_types,
        default=cell_types[:min(3, len(cell_types))],
        help="LAP will be computed for all cell type pairs"
    )

    if len(selected_cell_types) < 2:
        st.warning("⚠️ Please select at least 2 cell types")
        st.stop()

    # Representative cell selection method
    st.subheader("Representative Cell Selection")

    selection_method = st.radio(
        "Method:",
        ["Centroid-nearest", "First cell", "Random sample"],
        help="""
        - Centroid-nearest: Cell closest to the centroid (recommended)
        - First cell: First cell in the group
        - Random sample: Random sampling
        """
    )

    # LAP settings
    st.subheader("LAP Parameters")

    col_p1, col_p2 = st.columns(2)

    with col_p1:
        em_steps = st.number_input("EM steps:", 1, 10, 2,
                                  help="Number of EM optimization steps")

    with col_p2:
        min_lap_t = st.checkbox("Minimize initial T", value=False,
                               help="Start from minimum transition time")

    # ========================================
    # Step 3: Run LAP
    # ========================================
    st.header("Step 3: Run LAP Analysis")

    st.markdown(f"""
    **Configuration Summary:**
    - LAP calculation basis: `{lap_basis}` (high-dimensional space)
    - Visualization basis: `{vis_basis}` (for visualization)
    - Adjacency key: `{adj_key}`
    - Cell types: {len(selected_cell_types)}
    - Total transitions: {len(selected_cell_types) * (len(selected_cell_types) - 1)}
    - Selection method: {selection_method}
    - EM steps: {em_steps}
    """)

    if st.button("🚀 Compute LAP", type="primary"):
        with st.spinner("Computing LAP... This may take several minutes..."):
            try:
                # Representative cell selection
                representative_cells = {}

                progress_bar = st.progress(0)
                status_text = st.empty()

                # Step 1: Select representative cell for each cell type
                status_text.text("Selecting representative cells...")

                for ct in selected_cell_types:
                    ct_mask = adata.obs[group_key] == ct
                    ct_cells = adata.obs_names[ct_mask]

                    if selection_method == "Centroid-nearest":
                        # Select cell closest to centroid (computed in high-dimensional space)
                        ct_coords = adata.obsm[f'X_{lap_basis}'][ct_mask]
                        centroid = ct_coords.mean(axis=0)
                        distances = np.linalg.norm(ct_coords - centroid, axis=1)
                        nearest_idx = np.argmin(distances)
                        representative_cells[ct] = ct_cells[nearest_idx]

                    elif selection_method == "First cell":
                        representative_cells[ct] = ct_cells[0]

                    else:  # Random sample
                        representative_cells[ct] = np.random.choice(ct_cells)

                st.info(f"✅ Selected representative cells: {representative_cells}")

                # CRITICAL: Set use_for_pca=True for all genes before LAP calculation
                # This must be done BEFORE the loop and restored AFTER the loop
                original_use_for_pca = None
                if 'use_for_pca' in adata.var:
                    original_use_for_pca = adata.var['use_for_pca'].copy()
                    n_hvg = original_use_for_pca.sum()
                    if n_hvg < len(adata.var):
                        st.info(f"""
                        ⚙️ **Adjusting use_for_pca for LAP calculation:**
                        - Original: {n_hvg} genes (HVG only)
                        - Updated: {len(adata.var)} genes (all genes)

                        This is necessary because LAP projects back to gene space.
                        Will be restored after all LAP calculations complete.
                        """)
                        # IMPORTANT: Use [:] to modify array in-place, not replace with scalar
                        adata.var['use_for_pca'][:] = True

                # Step 2: LAP computation
                transition_graph_dict = {}
                n_transitions = len(selected_cell_types) * (len(selected_cell_types) - 1)
                current = 0

                try:
                    for i, start_ct in enumerate(selected_cell_types):
                        for j, end_ct in enumerate(selected_cell_types):
                            if start_ct != end_ct:
                                status_text.text(f"Computing {start_ct} → {end_ct}")

                                try:
                                    # Step 1: LAP computation on visualization basis (UMAP, etc.) - optional
                                    # This is only for saving the visualization path
                                    vis_lap_success = False
                                    try:
                                        # Check if Vector Field exists for vis_basis
                                        if f'VecFld_{vis_basis}' in adata.uns:
                                            st.caption(f"  Computing LAP on {vis_basis} for visualization...")

                                            # Find appropriate adjacency key
                                            vis_adj_candidates = [
                                                f'X_{vis_basis}_distances',
                                                f'{vis_basis}_distances',
                                                f'cosine_transition_matrix',
                                                adj_key
                                            ]
                                            vis_adj_key = None
                                            for key in vis_adj_candidates:
                                                if key in adata.obsp:
                                                    vis_adj_key = key
                                                    break

                                            if vis_adj_key:
                                                dyn.pd.least_action(
                                                    adata,
                                                    init_cells=[representative_cells[start_ct]],
                                                    target_cells=[representative_cells[end_ct]],
                                                    basis=vis_basis,
                                                    adj_key=vis_adj_key,
                                                    vf_key='VecFld',
                                                    min_lap_t=min_lap_t if current == 0 else False,
                                                    EM_steps=em_steps
                                                )
                                                vis_lap_success = True
                                    except Exception as e_vis:
                                        st.caption(f"  ⚠️ Skipping {vis_basis} LAP (visualization only): {str(e_vis)[:50]}")

                                    # Step 2: LAP computation in high-dimensional space (PCA, etc.) - required
                                    # This is the main LAP calculation
                                    st.caption(f"  Computing LAP on {lap_basis} for analysis...")
                                    lap_result = dyn.pd.least_action(
                                        adata,
                                        init_cells=[representative_cells[start_ct]],
                                        target_cells=[representative_cells[end_ct]],
                                        basis=lap_basis,
                                        adj_key=adj_key,
                                        vf_key='VecFld',
                                        min_lap_t=min_lap_t if (current == 0 and not vis_lap_success) else False,
                                        EM_steps=em_steps
                                    )

                                    # CRITICAL: Extract LAP data immediately after calculation
                                    # Deep copy to preserve data for later use
                                    import copy
                                    lap_dict = copy.deepcopy(adata.uns.get(f'LAP_{lap_basis}', {}))
                                    lap_vis_dict = copy.deepcopy(adata.uns.get(f'LAP_{vis_basis}', {}))

                                    # Validate LAP_dict structure
                                    if lap_dict:
                                        st.caption(f"  ✓ LAP_{lap_basis} keys: {list(lap_dict.keys())}")
                                        if 'prediction' in lap_dict:
                                            n_pred = len(lap_dict['prediction'])
                                            st.caption(f"  ✓ Saved {n_pred} prediction(s)")
                                            # Show shape of first prediction
                                            if n_pred > 0:
                                                pred_shape = lap_dict['prediction'][0].shape if hasattr(lap_dict['prediction'][0], 'shape') else len(lap_dict['prediction'][0])
                                                st.caption(f"  ✓ First prediction shape: {pred_shape}")
                                    else:
                                        st.caption(f"  ⚠️ No LAP_{lap_basis} data found in adata.uns")

                                    # Gene trajectory and ranking (computed from high-dimensional LAP)
                                    from dynamo.prediction import GeneTrajectory

                                    # Gene trajectory calculation
                                    # Fix: Use full gene mean when PCs shape doesn't match pca_mean
                                    if adata.uns['PCs'].shape[0] != adata.uns['pca_mean'].shape[0]:
                                        # Calculate mean for all genes
                                        if hasattr(adata.X, 'toarray'):
                                            full_mean = np.array(adata.X.mean(axis=0)).flatten()
                                        else:
                                            full_mean = adata.X.mean(axis=0)
                                        # Use full gene mean for GeneTrajectory
                                        gtraj = GeneTrajectory(adata, mean=full_mean)
                                    else:
                                        gtraj = GeneTrajectory(adata)
                                    gtraj.from_pca(lap_result.X, t=lap_result.t)
                                    gtraj.calc_msd()
                                    ranking = dyn.vf.rank_genes(adata, 'traj_msd', output_values=True)

                                    key = f"{start_ct}->{end_ct}"
                                    transition_graph_dict[key] = {
                                        'lap': lap_result,
                                        'action': lap_result.action_t()[-1] if hasattr(lap_result, 'action_t') else None,
                                        'time': lap_result.t[-1] if hasattr(lap_result, 't') else None,
                                        'LAP_dict': lap_dict,
                                        'LAP_vis_dict': lap_vis_dict,
                                        'ranking': ranking,
                                        'gtraj': gtraj,
                                    }

                                except Exception as e:
                                    st.warning(f"⚠️ Failed for {start_ct}→{end_ct}: {str(e)}")

                                current += 1
                                progress_bar.progress(current / n_transitions)

                    # After all transitions complete, save to session state
                    st.session_state.transition_graph = transition_graph_dict
                    st.session_state.lap_computed = True
                    st.session_state.selected_cell_types = selected_cell_types
                    st.session_state.representative_cells = representative_cells
                    st.session_state.lap_basis = lap_basis
                    st.session_state.vis_basis = vis_basis
                    st.session_state.lap_group_key = group_key

                    # Save LAP results to AnnData
                    status_text.text("Saving LAP results to AnnData...")
                    try:
                        adata = save_lap_to_adata(
                            adata,
                            transition_graph_dict,
                            selected_cell_types,
                            representative_cells,
                            lap_basis,
                            vis_basis,
                            group_key
                        )
                        st.session_state.lap_adata = adata
                        st.info("✅ LAP results saved to AnnData object")
                    except Exception as e:
                        st.warning(f"⚠️ Failed to save LAP results to AnnData: {str(e)}")
                        st.session_state.lap_adata = adata

                    progress_bar.empty()
                    status_text.empty()

                    st.success(f"✅ LAP computed for {len(transition_graph_dict)} transitions!")

                except Exception as e:
                    st.error(f"❌ Error: {str(e)}")
                    st.exception(e)
                finally:
                    # CRITICAL: Restore original use_for_pca after all LAP calculations
                    if original_use_for_pca is not None:
                        adata.var['use_for_pca'] = original_use_for_pca
                        st.info("✅ Restored original use_for_pca values")

            except Exception as e_outer:
                st.error(f"❌ Unexpected error during LAP computation: {str(e_outer)}")
                st.exception(e_outer)

    # ========================================
    # Step 4: Results
    # ========================================
    if st.session_state.lap_computed:
        st.markdown("---")
        st.header("Step 4: Results")

        adata_display = st.session_state.lap_adata
        lap_basis_display = st.session_state.lap_basis
        vis_basis_display = st.session_state.vis_basis
        group_key_display = st.session_state.lap_group_key

        st.info(f"""
        **Computed LAP paths:**
        - Calculation basis: `{lap_basis_display}` (high-dimensional space)
        - Visualization basis: `{vis_basis_display}` (for visualization)
        """)

        # Action/Time matrices
        st.subheader("Transition Matrices")

        cell_types_computed = st.session_state.selected_cell_types
        action_matrix = pd.DataFrame(
            index=cell_types_computed,
            columns=cell_types_computed,
            dtype=float
        )
        time_matrix = pd.DataFrame(
            index=cell_types_computed,
            columns=cell_types_computed,
            dtype=float
        )

        for key, data in st.session_state.transition_graph.items():
            start, end = key.split('->')
            if data['action'] is not None:
                action_matrix.loc[start, end] = data['action']
            if data['time'] is not None:
                time_matrix.loc[start, end] = data['time']

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("**Action Matrix**")
            plt.close('all')
            fig_action, ax_action = plt.subplots(figsize=(8, 7))
            sns.heatmap(action_matrix.astype(float), annot=True, fmt='.2f', ax=ax_action, cmap='viridis')
            ax_action.set_title("Action (energetic cost)")
            st.pyplot(fig_action)

            # Save PNG and PDF
            col_dl1, col_dl2 = st.columns(2)

            with col_dl1:
                fig_path_png = os.path.join(dynamo_lap_temp_dir, "action_matrix.png")
                fig_action.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                with open(fig_path_png, "rb") as f:
                    st.download_button(
                        "⬇️ Download PNG",
                        f,
                        file_name="action_matrix.png",
                        mime="image/png",
                        key="png_action"
                    )

            with col_dl2:
                try:
                    fig_path_pdf = os.path.join(dynamo_lap_temp_dir, "action_matrix.pdf")
                    fig_action.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                    with open(fig_path_pdf, "rb") as f:
                        st.download_button(
                            "⬇️ Download PDF",
                            f,
                            file_name="action_matrix.pdf",
                            mime="application/pdf",
                            key="pdf_action"
                        )
                except Exception as e:
                    st.warning(f"⚠️ Failed to save PDF: {str(e)}")

            plt.close(fig_action)

        with col2:
            st.markdown("**Time Matrix**")
            plt.close('all')
            fig_time, ax_time = plt.subplots(figsize=(8, 7))
            sns.heatmap(time_matrix.astype(float), annot=True, fmt='.2f', ax=ax_time, cmap='plasma')
            ax_time.set_title("Transition Time")
            st.pyplot(fig_time)

            # Save PNG and PDF
            col_dl1, col_dl2 = st.columns(2)

            with col_dl1:
                fig_path_png = os.path.join(dynamo_lap_temp_dir, "time_matrix.png")
                fig_time.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                with open(fig_path_png, "rb") as f:
                    st.download_button(
                        "⬇️ Download PNG",
                        f,
                        file_name="time_matrix.png",
                        mime="image/png",
                        key="png_time"
                    )

            with col_dl2:
                try:
                    fig_path_pdf = os.path.join(dynamo_lap_temp_dir, "time_matrix.pdf")
                    fig_time.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                    with open(fig_path_pdf, "rb") as f:
                        st.download_button(
                            "⬇️ Download PDF",
                            f,
                            file_name="time_matrix.pdf",
                            mime="application/pdf",
                            key="pdf_time"
                        )
                except Exception as e:
                    st.warning(f"⚠️ Failed to save PDF: {str(e)}")

            plt.close(fig_time)

        # Display matrices as tables
        st.subheader("Data Tables")

        col_t1, col_t2 = st.columns(2)

        with col_t1:
            st.markdown("**Action Matrix (CSV)**")
            st.dataframe(action_matrix)

        with col_t2:
            st.markdown("**Time Matrix (CSV)**")
            st.dataframe(time_matrix)

        # ========================================
        # Analysis Results Tabs
        # ========================================
        st.markdown("---")

        # Check if transition_graph has data
        transition_keys_list = list(st.session_state.transition_graph.keys())

        if len(transition_keys_list) == 0:
            st.warning("⚠️ No LAP results available. Please run LAP computation first.")
        else:
            st.subheader("Visualization")

            # Create tabs for different analysis types
            tab1, tab2, tab3 = st.tabs(["🎨 LAP Path Visualization", "📊 Gene Ranking", "🔥 Kinetic Heatmap"])

            # ========================================
            # Tab 1: LAP Path Visualization
            # ========================================
            with tab1:
                st.info(f"Displaying LAP paths on {vis_basis_display}")

                # Transition selection
                n_paths = st.slider(
                    "Number of transitions to visualize:",
                    min_value=1,
                    max_value=min(len(transition_keys_list), 10),
                    value=min(3, len(transition_keys_list)),
                    help="Display multiple LAP paths simultaneously"
                )

                selected_paths = []
                for i in range(n_paths):
                    path = st.selectbox(
                        f"Transition {i+1}:",
                        transition_keys_list,
                        index=min(i, len(transition_keys_list)-1),
                        key=f"path_select_{i}"
                    )
                    selected_paths.append(path)

                if st.button("🎨 Plot LAP Paths"):
                    with st.spinner("Generating LAP path visualization..."):
                        try:
                            from dynamo.plot.utils import map2color

                            plt.close('all')
                            fig, ax = plt.subplots(figsize=(10, 8))

                            # Base streamline plot
                            cluster_list = adata_display.obs[group_key_display].cat.categories.tolist() if adata_display.obs[group_key_display].dtype.name == 'category' else sorted(adata_display.obs[group_key_display].unique().tolist())
                            palette = create_cell_color_mapping(cluster_list, colormap_discrete)

                            dyn.pl.streamline_plot(
                                adata_display,
                                color=[group_key_display],
                                basis=vis_basis_display,
                                palette=palette,
                                ax=ax,
                                save_show_or_return="return"
                            )

                            # Overlay LAP paths
                            for path in selected_paths:
                                lap_vis_dict = st.session_state.transition_graph[path].get('LAP_vis_dict')

                                if lap_vis_dict is not None and 'prediction' in lap_vis_dict:
                                    predictions = lap_vis_dict['prediction']
                                    actions = lap_vis_dict['action']

                                    for prediction, action in zip(predictions, actions):
                                        # Plot path
                                        ax.plot(*prediction[:, [0, 1]].T, c='black', linewidth=2, alpha=0.7)
                                        # Plot points colored by action
                                        colors = map2color(action)
                                        ax.scatter(*prediction[:, [0, 1]].T, c=colors, s=30, alpha=0.8, edgecolors='black', linewidths=0.5)
                                else:
                                    st.warning(f"⚠️ LAP path data not available for {path}")

                            ax.set_title(f"LAP Paths: {', '.join(selected_paths)}", fontsize=14)

                            st.pyplot(fig)

                            # Save
                            col_dl1, col_dl2 = st.columns(2)

                            with col_dl1:
                                fig_path_png = os.path.join(dynamo_lap_temp_dir, "lap_paths.png")
                                fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                                with open(fig_path_png, "rb") as f:
                                    st.download_button(
                                        "⬇️ Download PNG",
                                        f,
                                        file_name="lap_paths.png",
                                        mime="image/png",
                                        key="png_lap_paths"
                                    )

                            with col_dl2:
                                try:
                                    fig_path_pdf = os.path.join(dynamo_lap_temp_dir, "lap_paths.pdf")
                                    fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                                    with open(fig_path_pdf, "rb") as f:
                                        st.download_button(
                                            "⬇️ Download PDF",
                                            f,
                                            file_name="lap_paths.pdf",
                                            mime="application/pdf",
                                            key="pdf_lap_paths"
                                        )
                                except Exception as e:
                                    st.warning(f"⚠️ Failed to save PDF: {str(e)}")

                            plt.close(fig)

                        except Exception as e:
                            st.error(f"Error: {str(e)}")
                            st.exception(e)

            # ========================================
            # Tab 2: Gene Ranking
            # ========================================
            with tab2:
                st.subheader("Gene Ranking (Mean Squared Displacement)")

                transition_for_genes = st.selectbox(
                    "Select transition:",
                    transition_keys_list,
                    index=0 if len(transition_keys_list) > 0 else None,
                    key="gene_rank_transition"
                )

                top_n_genes = st.slider("Top N genes:", 5, 50, 10)

                if st.button("📊 Plot Gene Rankings"):
                    if transition_for_genes is None:
                        st.error("⚠️ Please select a transition")
                    else:
                        with st.spinner("Generating gene ranking plot..."):
                            try:
                                ranking = st.session_state.transition_graph[transition_for_genes].get('ranking')

                                if ranking is not None:
                                    plt.close('all')
                                    fig, ax = plt.subplots(figsize=(10, 6))

                                    top_genes = ranking.head(top_n_genes)

                                    sns.barplot(
                                        y='all',
                                        x='all_values',
                                        data=top_genes,
                                        ax=ax,
                                        color='steelblue'
                                    )

                                    ax.set_title(f"Top {top_n_genes} genes for {transition_for_genes}", fontsize=14)
                                    ax.set_xlabel("Ranking Score (MSD)", fontsize=12)
                                    ax.set_ylabel("Genes", fontsize=12)

                                    st.pyplot(fig)

                                    # Save
                                    col_dl1, col_dl2 = st.columns(2)

                                    with col_dl1:
                                        fig_path_png = os.path.join(dynamo_lap_temp_dir, "gene_ranking.png")
                                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                                        with open(fig_path_png, "rb") as f:
                                            st.download_button(
                                                "⬇️ Download PNG",
                                                f,
                                                file_name="gene_ranking.png",
                                                mime="image/png",
                                                key="png_gene_rank"
                                            )

                                    with col_dl2:
                                        try:
                                            fig_path_pdf = os.path.join(dynamo_lap_temp_dir, "gene_ranking.pdf")
                                            fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                                            with open(fig_path_pdf, "rb") as f:
                                                st.download_button(
                                                    "⬇️ Download PDF",
                                                    f,
                                                    file_name="gene_ranking.pdf",
                                                    mime="application/pdf",
                                                    key="pdf_gene_rank"
                                                )
                                        except Exception as e:
                                            st.warning(f"⚠️ Failed to save PDF: {str(e)}")

                                    plt.close(fig)

                                    # Show table
                                    st.markdown("**Top genes table:**")
                                    st.dataframe(top_genes)

                                else:
                                    st.error("❌ Gene ranking data not available for this transition")

                            except Exception as e:
                                st.error(f"Error: {str(e)}")
                                st.exception(e)

            # ========================================
            # Tab 3: Kinetic Heatmap
            # ========================================
            with tab3:
                st.subheader("Kinetic Heatmap")

                st.markdown("""
                Heatmap of gene expression dynamics along the LAP path
                """)

                # Heatmap options in a form
                with st.form("heatmap_options_form"):
                    transition_for_heatmap = st.selectbox(
                        "Select transition:",
                        transition_keys_list,
                        key="heatmap_transition"
                    )

                    col_h1, col_h2, col_h3 = st.columns(3)

                    with col_h1:
                        n_genes_heatmap = st.slider("Number of genes:", 10, 200, 50)

                    with col_h2:
                        genes_to_use = st.radio(
                            "Gene selection:",
                            ["Top ranked genes", "Use for transition (if available)", "Custom gene list"],
                            help="Method for selecting genes"
                        )

                    with col_h3:
                        filter_tfs = st.checkbox(
                            "Filter for TFs only",
                            value=False,
                            help="Same as Shiny app: Filter only TFs from use_for_transition genes"
                        )

                    # Custom gene list input (always shown, used only when "Custom gene list" is selected)
                    custom_gene_input = st.text_area(
                        "Custom gene list (if using 'Custom gene list' option above):",
                        height=100,
                        placeholder="e.g., Gene1 Gene2, Gene3\tGene4\nGene5",
                        help="Gene names separated by space/comma/tab/newline. Comparison is case insensitive."
                    )

                    # Additional heatmap options
                    st.markdown("#### Heatmap Options")
                    col_opt1, col_opt2, col_opt3, col_opt4 = st.columns(4)

                    with col_opt1:
                        gene_order_method = st.selectbox(
                            "Gene ordering:",
                            ["Peak time", "Hierarchical clustering", "As selected"],
                            help="Peak time: Order by gene expression peak time\nHierarchical: Hierarchical clustering\nAs selected: Original ranking order"
                        )

                    with col_opt2:
                        heatmap_cmap = st.selectbox(
                            "Colormap:",
                            ["bwr", "RdBu_r", "coolwarm", "viridis", "plasma"],
                            help="Heatmap color scheme"
                        )

                    with col_opt3:
                        normalization_method = st.radio(
                            "Normalization:",
                            ["Z-score", "Standard scale (Min-Max)"],
                            index=0,
                            help="Z-score: Normalize to mean 0, standard deviation 1\nStandard scale: Dynamo default, scale to 0-1"
                        )

                    with col_opt4:
                        show_gene_labels = st.checkbox(
                            "Show gene names",
                            value=True,
                            help="Display gene names on Y-axis"
                        )

                    # Figure size options
                    st.markdown("#### Figure Size")
                    col_size1, col_size2 = st.columns(2)

                    with col_size1:
                        fig_width = st.slider(
                            "Figure width:",
                            min_value=6,
                            max_value=24,
                            value=12,
                            step=1,
                            help="Heatmap width (inches)"
                        )

                    with col_size2:
                        fig_height = st.slider(
                            "Figure height:",
                            min_value=4,
                            max_value=30,
                            value=10,
                            step=1,
                            help="Heatmap height (inches). Increase for large number of genes"
                        )

                    # Submit button
                    submit_heatmap = st.form_submit_button("🔥 Generate Kinetic Heatmap", type="primary")

                if submit_heatmap:
                    with st.spinner("Generating kinetic heatmap..."):
                        try:
                            transition_data = st.session_state.transition_graph[transition_for_heatmap]
                            ranking = transition_data.get('ranking')

                            if ranking is not None:
                                            # Copy adata and add LAP data (following Shiny app approach)
                                            # Original Shiny code:
                                            #   _adata = adata.copy()
                                            #   _adata.uns["LAP_umap"] = transition_graph()[path]["LAP_umap"]
                                            #   _adata.uns["LAP_pca"] = transition_graph()[path]["LAP_pca"]
                                            _adata = adata_display.copy()

                                            # Check if LAP_dict is available
                                            lap_dict = transition_data.get('LAP_dict')
                                            lap_vis_dict = transition_data.get('LAP_vis_dict')

                                            if not lap_dict:
                                                st.error("❌ LAP data not available for this transition. Please recompute LAP.")
                                                st.stop()

                                            # Validate LAP_dict structure
                                            if 'prediction' not in lap_dict or 'action' not in lap_dict:
                                                st.error("❌ Invalid LAP data structure. Missing 'prediction' or 'action' keys.")
                                                st.info(f"Available keys: {list(lap_dict.keys())}")
                                                st.stop()

                                            # Add LAP dictionaries to adata.uns (matching Shiny app structure)
                                            # Shiny saves both LAP_umap and LAP_pca in adata.uns
                                            # IMPORTANT: Dynamo stores trajectory data as lists when multiple trajectories exist
                                            # For single trajectory, unwrap only specific keys: prediction, action, t, mftp, exprs
                                            # Also convert pandas Index to list for h5ad compatibility

                                            trajectory_keys = ['prediction', 'action', 't', 'mftp', 'exprs']

                                            lap_dict_unwrapped = {}
                                            for key, value in lap_dict.items():
                                                if key in trajectory_keys and isinstance(value, list) and len(value) == 1:
                                                    # Single trajectory stored as list - unwrap it
                                                    lap_dict_unwrapped[key] = value[0]
                                                elif key == 'genes' and hasattr(value, 'tolist'):
                                                    # Convert pandas Index to list for h5ad compatibility
                                                    lap_dict_unwrapped[key] = value.tolist()
                                                else:
                                                    lap_dict_unwrapped[key] = value

                                            _adata.uns[f'LAP_{lap_basis_display}'] = lap_dict_unwrapped

                                            if lap_vis_dict:
                                                lap_vis_dict_unwrapped = {}
                                                for key, value in lap_vis_dict.items():
                                                    if key in trajectory_keys and isinstance(value, list) and len(value) == 1:
                                                        lap_vis_dict_unwrapped[key] = value[0]
                                                    elif key == 'genes' and hasattr(value, 'tolist'):
                                                        lap_vis_dict_unwrapped[key] = value.tolist()
                                                    else:
                                                        lap_vis_dict_unwrapped[key] = value
                                                _adata.uns[f'LAP_{vis_basis_display}'] = lap_vis_dict_unwrapped

                                            # Select genes following Shiny app implementation
                                            if genes_to_use == "Top ranked genes":
                                                genes_for_heatmap = ranking['all'].head(n_genes_heatmap).tolist()
                                            elif genes_to_use == "Custom gene list":
                                                # Parse custom gene list input
                                                if custom_gene_input:
                                                    import re
                                                    # Split by space, comma, tab, or newline
                                                    input_genes = re.split(r'[\s,\t\n]+', custom_gene_input.strip())
                                                    # Remove empty strings
                                                    input_genes = [g.strip() for g in input_genes if g.strip()]

                                                    if len(input_genes) == 0:
                                                        st.error("❌ No genes found in custom input. Please enter gene names.")
                                                        st.stop()

                                                    # Case insensitive matching
                                                    # Create mapping: lowercase -> original gene name in dataset
                                                    var_names_lower_map = {g.lower(): g for g in _adata.var_names}

                                                    # Match input genes (case insensitive)
                                                    genes_for_heatmap = []
                                                    not_found = []
                                                    for input_gene in input_genes:
                                                        gene_lower = input_gene.lower()
                                                        if gene_lower in var_names_lower_map:
                                                            # Use the original gene name from dataset
                                                            genes_for_heatmap.append(var_names_lower_map[gene_lower])
                                                        else:
                                                            not_found.append(input_gene)

                                                    if len(genes_for_heatmap) == 0:
                                                        st.error(f"❌ None of the input genes were found in the dataset.")
                                                        if len(not_found) > 0:
                                                            st.info(f"Not found: {', '.join(not_found[:10])}" +
                                                                   (f" ... and {len(not_found)-10} more" if len(not_found) > 10 else ""))
                                                        st.stop()

                                                    # Show warning if some genes were not found
                                                    if len(not_found) > 0:
                                                        st.warning(f"⚠️ {len(not_found)} gene(s) not found in dataset: {', '.join(not_found[:5])}" +
                                                                 (f" ... and {len(not_found)-5} more" if len(not_found) > 5 else ""))

                                                    st.success(f"✓ Matched {len(genes_for_heatmap)} genes from custom list (case insensitive)")
                                                else:
                                                    st.error("❌ Please enter gene names in the text area.")
                                                    st.stop()
                                            else:  # "Use for transition (if available)"
                                                # Use genes marked for transition (Shiny app approach)
                                                if 'use_for_transition' in _adata.var.columns:
                                                    use_for_trans_genes = _adata.var_names[_adata.var.use_for_transition]
                                                    genes_for_heatmap = use_for_trans_genes[:n_genes_heatmap].tolist()
                                                else:
                                                    # Fallback to ranked genes
                                                    genes_for_heatmap = ranking['all'].head(n_genes_heatmap).tolist()
                                                    st.warning("use_for_transition not found, using ranked genes instead")

                                            # Apply TF filter if requested (after gene selection)
                                            if filter_tfs:
                                                # Load TFs list (case insensitive)
                                                tfs_names_lower = load_tfs()
                                                if tfs_names_lower is not None:
                                                    # Filter for TFs only with case-insensitive matching
                                                    genes_before_filter = len(genes_for_heatmap)
                                                    genes_for_heatmap = [gene for gene in genes_for_heatmap
                                                                         if gene.lower() in tfs_names_lower]
                                                    if len(genes_for_heatmap) == 0:
                                                        st.error("❌ No TFs found in selected genes. Please disable TF filter or select different genes.")
                                                        st.stop()
                                                else:
                                                    st.warning("⚠️ Could not load TFs, TF filter not applied")

                                            plt.close('all')

                                            # Try using dynamo's kinetic_heatmap first (following Shiny app implementation)
                                            fig = None
                                            try:
                                                st.caption("Attempting to generate heatmap using dynamo.pl.kinetic_heatmap (Shiny method)...")
                                                # Generate kinetic heatmap following Shiny app parameters exactly
                                                sns.set(font_scale=0.8)
                                                sns_heatmap = dyn.pl.kinetic_heatmap(
                                                    _adata,
                                                    basis=lap_basis_display,
                                                    mode='lap',
                                                    genes=genes_for_heatmap,
                                                    project_back_to_high_dim=True,
                                                    traj_ind=0,  # Use first (and only) trajectory for single transition
                                                    save_show_or_return='return',
                                                    color_map=heatmap_cmap,
                                                    transpose=False,
                                                    xticklabels=False,  # Match Shiny app
                                                    yticklabels=True,   # Match Shiny app
                                                )

                                                # Shiny app code: plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels())
                                                plt.setp(sns_heatmap.ax_heatmap.yaxis.get_majorticklabels())
                                                plt.tight_layout()
                                                fig = plt.gcf()  # Get current figure for saving
                                                st.pyplot(fig)
                                                st.success("✅ Kinetic heatmap generated using Dynamo's built-in method!")
                                            except IndexError as idx_err:
                                                st.warning(f"⚠️ dynamo.pl.kinetic_heatmap failed: {str(idx_err)}")
                                                st.info("""
                                                **Note**: The Shiny app does not have a fallback for this error.
                                                Attempting custom alternative heatmap generation using GeneTrajectory...
                                                """)

                                                # Alternative implementation matching Dynamo's kinetic_heatmap processing
                                                # This uses GeneTrajectory directly with Dynamo-compatible processing
                                                try:
                                                    gtraj = transition_data.get('gtraj')
                                                    if gtraj is None:
                                                        st.error("❌ GeneTrajectory object not available")
                                                        raise ValueError("GeneTrajectory not found")

                                                    # Get time values for trajectory
                                                    lap_dict_for_time = transition_data.get('LAP_dict')
                                                    if 't' in lap_dict_for_time:
                                                        t_values = lap_dict_for_time['t']
                                                        if isinstance(t_values, list):
                                                            t_values = t_values[0] if len(t_values) > 0 else np.arange(gtraj.X.shape[0])
                                                    else:
                                                        t_values = np.arange(gtraj.X.shape[0])

                                                    # Get gene expression along trajectory
                                                    # gtraj.X contains the gene expression trajectory
                                                    gene_indices = [_adata.var_names.get_loc(g) for g in genes_for_heatmap if g in _adata.var_names]
                                                    valid_genes = [g for g in genes_for_heatmap if g in _adata.var_names]

                                                    if len(gene_indices) == 0:
                                                        st.error("❌ No valid genes found in dataset")
                                                        raise ValueError("No valid genes")

                                                    st.info(f"Using {len(valid_genes)} valid genes (out of {len(genes_for_heatmap)} requested)")

                                                    # Extract expression data (trajectory_points x genes)
                                                    expr_traj = gtraj.X[:, gene_indices]

                                                    # STEP 1: Apply LOESS smoothing (matching Dynamo's lowess_smoother)
                                                    # Import Dynamo's smoother
                                                    from dynamo.plot.time_series import lowess_smoother

                                                    st.caption("Applying LOESS smoothing (n_convolve=30)...")
                                                    # Transpose for smoother: genes x time_points
                                                    expr_traj_smooth = lowess_smoother(
                                                        t_values,
                                                        expr_traj.T,
                                                        spaced_num=None,  # No resampling
                                                        n_convolve=30     # Dynamo default
                                                    )
                                                    # Result: genes x time_points

                                                    # STEP 2: Normalization
                                                    # Matching Dynamo's standard_scale (lines 343-346)
                                                    if normalization_method == "Z-score":
                                                        # Z-score normalization
                                                        from scipy import stats
                                                        expr_traj_norm = stats.zscore(expr_traj_smooth, axis=1, ddof=1, nan_policy='omit')
                                                        cbar_label = 'Z-score'
                                                        vmin, vmax = -2, 2
                                                    else:  # "Standard scale (Min-Max)"
                                                        # Standard scale (Min-Max) - Dynamo default
                                                        expr_min = np.min(expr_traj_smooth, axis=1, keepdims=True)
                                                        expr_ptp = np.ptp(expr_traj_smooth, axis=1, keepdims=True)  # peak-to-peak
                                                        expr_traj_norm = (expr_traj_smooth - expr_min) / (expr_ptp + 1e-10)
                                                        cbar_label = 'Expression (scaled)'
                                                        vmin, vmax = 0, 1

                                                    # STEP 3: Gene ordering
                                                    if gene_order_method == "Peak time":
                                                        # Sort by peak time (matching Dynamo's "maximum" method)
                                                        peak_indices = np.argmax(expr_traj_norm, axis=1)
                                                        sort_order = np.argsort(peak_indices)
                                                        st.caption("Ordering genes by peak expression time...")
                                                    elif gene_order_method == "Hierarchical clustering":
                                                        # Hierarchical clustering
                                                        from scipy.cluster.hierarchy import linkage, dendrogram
                                                        linkage_matrix = linkage(expr_traj_norm, method='ward')
                                                        dendro = dendrogram(linkage_matrix, no_plot=True)
                                                        sort_order = dendro['leaves']
                                                        st.caption("Ordering genes by hierarchical clustering...")
                                                    else:  # "As selected"
                                                        # Keep original ranking order
                                                        sort_order = np.arange(len(valid_genes))
                                                        st.caption("Using original gene ranking order...")

                                                    # Apply sorting
                                                    expr_traj_sorted = expr_traj_norm[sort_order, :]
                                                    genes_sorted = [valid_genes[i] for i in sort_order]

                                                    # STEP 4: Create heatmap using seaborn (matching Dynamo style)
                                                    # Create DataFrame for seaborn
                                                    n_points = expr_traj_sorted.shape[1]
                                                    time_labels = [f'{t_values[i]:.2f}' for i in range(n_points)]
                                                    df_heatmap = pd.DataFrame(
                                                        expr_traj_sorted,
                                                        index=genes_sorted,
                                                        columns=time_labels
                                                    )

                                                    # Choose between clustermap (with dendrogram) or regular heatmap
                                                    if gene_order_method == "Hierarchical clustering":
                                                        # Use clustermap with dendrogram (matching Dynamo's approach)
                                                        st.caption("Creating clustermap with dendrogram...")

                                                        # Use user-specified figure size
                                                        figsize = (fig_width, fig_height)

                                                        g = sns.clustermap(
                                                            df_heatmap,
                                                            cmap=heatmap_cmap,
                                                            vmin=vmin,
                                                            vmax=vmax,
                                                            row_cluster=True,   # Cluster genes
                                                            col_cluster=False,  # Don't cluster time (preserve order)
                                                            figsize=figsize,
                                                            yticklabels=show_gene_labels,
                                                            xticklabels=False,  # Too many time points
                                                            cbar_kws={'label': cbar_label},
                                                            dendrogram_ratio=(0.15, 0.0),  # Only row dendrogram
                                                            cbar_pos=(0.02, 0.8, 0.03, 0.15),  # Colorbar position
                                                            method='ward'  # Linkage method
                                                        )

                                                        # Set labels on heatmap axes
                                                        g.ax_heatmap.set_xlabel('LAP trajectory time', fontsize=12)
                                                        g.ax_heatmap.set_ylabel('Genes', fontsize=12)

                                                        # Add title
                                                        g.fig.suptitle(
                                                            f'Kinetic Heatmap: {transition_for_heatmap}',
                                                            fontsize=14,
                                                            fontweight='bold',
                                                            y=0.98
                                                        )

                                                        # Add X-axis time labels at key positions
                                                        tick_positions = [0, n_points//4, n_points//2, 3*n_points//4, n_points-1]
                                                        g.ax_heatmap.set_xticks(tick_positions)
                                                        g.ax_heatmap.set_xticklabels([
                                                            f'{t_values[0]:.2f}',
                                                            f'{t_values[n_points//4]:.2f}',
                                                            f'{t_values[n_points//2]:.2f}',
                                                            f'{t_values[3*n_points//4]:.2f}',
                                                            f'{t_values[-1]:.2f}'
                                                        ], rotation=0)

                                                        fig = g.fig

                                                    else:
                                                        # Use regular heatmap (no dendrogram)
                                                        st.caption("Creating heatmap...")

                                                        # Use user-specified figure size
                                                        figsize = (fig_width, fig_height)
                                                        fig, ax = plt.subplots(figsize=figsize)

                                                        # Create heatmap with seaborn
                                                        sns.heatmap(
                                                            df_heatmap,
                                                            cmap=heatmap_cmap,
                                                            vmin=vmin,
                                                            vmax=vmax,
                                                            yticklabels=show_gene_labels,
                                                            xticklabels=False,  # Too many time points
                                                            cbar_kws={'label': cbar_label},
                                                            ax=ax
                                                        )

                                                        # Set labels
                                                        ax.set_xlabel('LAP trajectory time', fontsize=12)
                                                        ax.set_ylabel('Genes', fontsize=12)
                                                        ax.set_title(
                                                            f'Kinetic Heatmap: {transition_for_heatmap}',
                                                            fontsize=14,
                                                            fontweight='bold'
                                                        )

                                                        # Add X-axis time labels at key positions
                                                        tick_positions = [0, n_points//4, n_points//2, 3*n_points//4, n_points-1]
                                                        ax.set_xticks(tick_positions)
                                                        ax.set_xticklabels([
                                                            f'{t_values[0]:.2f}',
                                                            f'{t_values[n_points//4]:.2f}',
                                                            f'{t_values[n_points//2]:.2f}',
                                                            f'{t_values[3*n_points//4]:.2f}',
                                                            f'{t_values[-1]:.2f}'
                                                        ], rotation=0)

                                                        plt.tight_layout()

                                                    st.pyplot(fig)
                                                    st.success("✅ Kinetic heatmap generated with Dynamo-compatible processing (using seaborn)!")

                                                except Exception as alt_err:
                                                    st.error(f"❌ Alternative heatmap generation also failed: {str(alt_err)}")
                                                    st.exception(alt_err)

                                                    # Show detailed LAP dict structure for debugging
                                                    with st.expander("🔍 Detailed LAP_dict structure"):
                                                        for key, value in lap_dict.items():
                                                            if isinstance(value, (list, np.ndarray)):
                                                                st.write(f"- {key}: {type(value).__name__} of length {len(value)}")
                                                                if len(value) > 0:
                                                                    st.write(f"  First element: {type(value[0])}, shape: {getattr(value[0], 'shape', 'N/A')}")
                                                            else:
                                                                st.write(f"- {key}: {type(value).__name__} = {value}")
                                                    raise
                                            except Exception as e:
                                                st.error(f"❌ Unexpected error: {str(e)}")
                                                st.exception(e)
                                                raise

                                            # Save
                                            col_dl1, col_dl2 = st.columns(2)

                                            with col_dl1:
                                                fig_path_png = os.path.join(dynamo_lap_temp_dir, "kinetic_heatmap.png")
                                                fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                                                with open(fig_path_png, "rb") as f:
                                                    st.download_button(
                                                        "⬇️ Download PNG",
                                                        f,
                                                        file_name="kinetic_heatmap.png",
                                                        mime="image/png",
                                                        key="png_heatmap"
                                                    )

                                            with col_dl2:
                                                try:
                                                    fig_path_pdf = os.path.join(dynamo_lap_temp_dir, "kinetic_heatmap.pdf")
                                                    fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                                                    with open(fig_path_pdf, "rb") as f:
                                                        st.download_button(
                                                            "⬇️ Download PDF",
                                                            f,
                                                            file_name="kinetic_heatmap.pdf",
                                                            mime="application/pdf",
                                                            key="pdf_heatmap"
                                                        )
                                                except Exception as e:
                                                    st.warning(f"⚠️ Failed to save PDF: {str(e)}")

                                            plt.close(fig)

                            else:
                                st.error("❌ Ranking data not available for this transition")

                        except Exception as e:
                            st.error(f"Error generating kinetic heatmap: {str(e)}")
                            st.exception(e)

        # ========================================
        # Export Results
        # ========================================
        st.markdown("---")
        st.subheader("📥 Export Results")

        st.markdown("""
        Export LAP analysis results as TSV files (tab-delimited).

        **Files saved (bundled in ZIP):**
        - ✅ `action_matrix.tsv` - Action matrix
        - ✅ `time_matrix.tsv` - Time matrix
        - ✅ `gene_rankings/` - Gene rankings for each transition
        - ✅ `representative_cells.tsv` - Representative cell information
        - ✅ `metadata.tsv` - LAP analysis metadata
        """)

        # Show what's in the AnnData object
        with st.expander("🔍 View saved LAP data in AnnData"):
            if 'LAP_metadata' in adata_display.uns:
                st.markdown("**LAP Metadata:**")
                st.json(adata_display.uns['LAP_metadata'])

            if 'LAP_representative_cells' in adata_display.uns:
                st.markdown("**Representative Cells:**")
                st.json(adata_display.uns['LAP_representative_cells'])

            if 'LAP_transitions' in adata_display.uns:
                st.markdown(f"**Number of transitions saved:** {len(adata_display.uns['LAP_transitions'])}")
                st.markdown("**Transition keys:**")
                st.write(list(adata_display.uns['LAP_transitions'].keys()))

            if 'LAP_representative' in adata_display.obs.columns:
                rep_cells = adata_display.obs[adata_display.obs['LAP_representative'] != 'No']
                st.markdown(f"**Representative cells marked in .obs:** {len(rep_cells)}")

        # Export button
        col_exp1, col_exp2 = st.columns([2, 1])

        with col_exp1:
            export_filename = st.text_input(
                "Output filename:",
                value="dynamo_lap_results.zip",
                help="Save TSV files bundled in a ZIP archive"
            )

        with col_exp2:
            st.write("")  # spacing
            st.write("")  # spacing
            if st.button("💾 Export TSV (ZIP)", type="primary", use_container_width=True):
                with st.spinner("Exporting LAP results to TSV files..."):
                    try:
                        import tempfile
                        import zipfile
                        from io import BytesIO

                        # Create a temporary directory
                        with tempfile.TemporaryDirectory() as tmpdir:
                            # 1. Export Action Matrix
                            st.caption("Exporting action_matrix.tsv...")
                            action_path = os.path.join(tmpdir, "action_matrix.tsv")
                            action_matrix.to_csv(action_path, sep='\t')

                            # 2. Export Time Matrix
                            st.caption("Exporting time_matrix.tsv...")
                            time_path = os.path.join(tmpdir, "time_matrix.tsv")
                            time_matrix.to_csv(time_path, sep='\t')

                            # 3. Export Gene Rankings (create subdirectory)
                            rankings_dir = os.path.join(tmpdir, "gene_rankings")
                            os.makedirs(rankings_dir, exist_ok=True)

                            st.caption("Exporting gene rankings...")
                            for transition_key, transition_data in st.session_state.transition_graph.items():
                                ranking = transition_data.get('ranking')
                                if ranking is not None:
                                    # Sanitize filename (replace -> with _to_)
                                    safe_name = transition_key.replace('->', '_to_').replace(' ', '_')
                                    ranking_path = os.path.join(rankings_dir, f"{safe_name}.tsv")
                                    ranking.to_csv(ranking_path, sep='\t')

                            # 4. Export Representative Cells
                            st.caption("Exporting representative_cells.tsv...")
                            rep_cells_df = pd.DataFrame([
                                {'cell_type': ct, 'representative_cell': cell}
                                for ct, cell in st.session_state.representative_cells.items()
                            ])
                            rep_cells_path = os.path.join(tmpdir, "representative_cells.tsv")
                            rep_cells_df.to_csv(rep_cells_path, sep='\t', index=False)

                            # 5. Export Metadata
                            st.caption("Exporting metadata.tsv...")
                            metadata_df = pd.DataFrame([
                                {'parameter': 'lap_basis', 'value': lap_basis_display},
                                {'parameter': 'vis_basis', 'value': vis_basis_display},
                                {'parameter': 'group_key', 'value': group_key_display},
                                {'parameter': 'n_cell_types', 'value': len(cell_types_computed)},
                                {'parameter': 'n_transitions', 'value': len(st.session_state.transition_graph)},
                                {'parameter': 'cell_types', 'value': ', '.join(cell_types_computed)}
                            ])
                            metadata_path = os.path.join(tmpdir, "metadata.tsv")
                            metadata_df.to_csv(metadata_path, sep='\t', index=False)

                            # Create ZIP file
                            st.caption("Creating ZIP archive...")
                            zip_buffer = BytesIO()
                            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zipf:
                                # Add all files to ZIP
                                zipf.write(action_path, "action_matrix.tsv")
                                zipf.write(time_path, "time_matrix.tsv")
                                zipf.write(rep_cells_path, "representative_cells.tsv")
                                zipf.write(metadata_path, "metadata.tsv")

                                # Add gene ranking files
                                for root, dirs, files in os.walk(rankings_dir):
                                    for file in files:
                                        file_path = os.path.join(root, file)
                                        arcname = os.path.join("gene_rankings", file)
                                        zipf.write(file_path, arcname)

                            zip_buffer.seek(0)
                            zip_bytes = zip_buffer.read()

                        # Provide download button
                        st.success("✅ TSV files created successfully!")
                        st.download_button(
                            label="⬇️ Download ZIP file",
                            data=zip_bytes,
                            file_name=export_filename,
                            mime="application/zip",
                            key="download_zip"
                        )

                        st.info(f"""
                        **Export summary:**
                        - File size: {len(zip_bytes) / (1024**2):.2f} MB
                        - Files included:
                          - action_matrix.tsv
                          - time_matrix.tsv
                          - representative_cells.tsv
                          - metadata.tsv
                          - gene_rankings/ ({len(st.session_state.transition_graph)} files)
                        - Total: {5 + len(st.session_state.transition_graph)} files
                        """)

                    except Exception as e:
                        st.error(f"❌ Failed to export TSV files: {str(e)}")
                        st.exception(e)

else:
    st.info("👆 Please upload an h5ad file with completed Dynamo analysis to begin.")
