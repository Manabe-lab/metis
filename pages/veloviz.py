"""
VeloViz: RNA Velocity-Informed Embedding
Create 2D embeddings informed by RNA velocity using VeloViz (R package)
"""

import streamlit as st
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import io
import time
import tempfile
from scipy.sparse import issparse
from helper_func import clear_old_directories, clear_old_files

# ========================================
# Cached functions
# ========================================

@st.cache_data
def check_r_veloviz_installed():
    """Check if R and veloviz are installed"""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr

        # Check if veloviz is installed
        ro.r('library(veloviz)')
        return True, "VeloViz is available"
    except Exception as e:
        return False, str(e)


def run_veloviz(curr_matrix, proj_matrix, gene_names, cell_names,
                n_pcs=20, k=30, similarity_threshold=0.2, distance_weight=1.0,
                normalize_depth=True, use_ods_genes=True, seed=42):
    """
    Run VeloViz using rpy2 with numpy arrays

    Parameters
    ----------
    curr_matrix : np.ndarray
        Current expression matrix (cells x genes), dense
    proj_matrix : np.ndarray
        Projected expression matrix (cells x genes), dense
    gene_names : list
        Gene names for rownames
    cell_names : list
        Cell names for colnames
    n_pcs : int
        Number of PCs to use
    k : int
        Number of neighbors
    similarity_threshold : float
        Similarity threshold for graph construction
    distance_weight : float
        Weight for distance component
    normalize_depth : bool
        Whether to normalize by depth
    use_ods_genes : bool
        Whether to use overdispersed genes
    seed : int
        Random seed

    Returns
    -------
    np.ndarray
        VeloViz embedding coordinates (cells x 2)
    """
    import rpy2.robjects as ro
    from rpy2.robjects import numpy2ri
    from rpy2.robjects.vectors import StrVector

    # Activate numpy conversion
    numpy2ri.activate()

    try:
        # Transfer numpy arrays to R (cells x genes)
        ro.globalenv['curr_py'] = curr_matrix
        ro.globalenv['proj_py'] = proj_matrix
        ro.globalenv['gene_names'] = StrVector(gene_names)
        ro.globalenv['cell_names'] = StrVector(cell_names)
        ro.globalenv['n_pcs'] = n_pcs
        ro.globalenv['k_neighbors'] = k
        ro.globalenv['sim_thresh'] = similarity_threshold
        ro.globalenv['dist_weight'] = distance_weight
        ro.globalenv['norm_depth'] = normalize_depth
        ro.globalenv['use_ods'] = use_ods_genes
        ro.globalenv['seed_val'] = seed

        # Run VeloViz in R
        ro.r('''
        library(veloviz)
        set.seed(seed_val)

        # Transpose matrices (R expects genes x cells)
        curr_t <- t(curr_py)
        proj_t <- t(proj_py)

        # Set row and column names
        rownames(curr_t) <- gene_names
        colnames(curr_t) <- cell_names
        rownames(proj_t) <- gene_names
        colnames(proj_t) <- cell_names

        # Run buildVeloviz
        veloviz_result <- buildVeloviz(
            curr = curr_t,
            proj = proj_t,
            normalize.depth = norm_depth,
            use.ods.genes = use_ods,
            alpha = 0.05,
            pca = TRUE,
            center = TRUE,
            scale = TRUE,
            nPCs = n_pcs,
            k = k_neighbors,
            similarity.threshold = sim_thresh,
            distance.weight = dist_weight,
            distance.threshold = 1,
            weighted = TRUE,
            verbose = TRUE
        )

        # Extract coordinates
        veloviz_coords <- veloviz_result$fdg_coords
        ''')

        # Get coordinates back to Python
        coords = np.array(ro.globalenv['veloviz_coords'])

        return coords

    finally:
        # Deactivate to avoid conflicts
        numpy2ri.deactivate()

        # Clean up R environment
        try:
            ro.r('rm(list = c("curr_py", "proj_py", "curr_t", "proj_t", "gene_names", "cell_names", "veloviz_result", "veloviz_coords"))')
        except:
            pass


@st.cache_data
def create_color_mapping(categories, palette_name='tab20'):
    """Create consistent color mapping for categories"""
    unique_cats = sorted(categories.unique())
    n_cats = len(unique_cats)

    if n_cats <= 20:
        colors = sns.color_palette(palette_name, n_colors=max(n_cats, 1))
    else:
        colors = sns.color_palette('husl', n_colors=n_cats)

    return {cat: colors[i] for i, cat in enumerate(unique_cats)}


# ========================================
# Main App
# ========================================

st.set_page_config(page_title="VeloViz", page_icon="🌀", layout="wide")

st.title("🌀 VeloViz: RNA Velocity-Informed Embedding")
st.markdown("""
Create 2D embedding considering RNA velocity. VeloViz uses both current expression state and predicted future expression state to
construct embedding space for cell trajectory visualization.

### Features
- Unlike conventional UMAP/tSNE, reflects velocity information in embedding
- Cell differentiation direction is more clearly visualized
- Can be calculated directly from scVelo analysis results

### References
- [Atta et al. (2021) "VeloViz: RNA-velocity informed embeddings for visualizing cellular trajectories" Bioinformatics](https://doi.org/10.1093/bioinformatics/btab653)
- [VeloViz GitHub](https://github.com/JEFworks-Lab/veloviz)
""")

# ========================================
# Check R/VeloViz availability
# ========================================
veloviz_available, veloviz_message = check_r_veloviz_installed()

if not veloviz_available:
    st.error(f"""
    ❌ **VeloViz (R package) is not available**

    Error: {veloviz_message}

    ### Installation
    Run the following in R:
    ```r
    if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
    BiocManager::install("veloviz")
    ```

    Also, rpy2 is required in Python:
    ```bash
    pip install rpy2
    ```
    """)
    st.stop()

st.success("✓ VeloViz (R) is available")

# ========================================
# Initialize session state
# ========================================
if "veloviz_temp_dir" not in st.session_state:
    veloviz_temp_dir = os.path.join("temp", f"veloviz_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(veloviz_temp_dir, exist_ok=True)
    st.session_state.veloviz_temp_dir = veloviz_temp_dir

if "veloviz_adata" not in st.session_state:
    st.session_state.veloviz_adata = None

if "veloviz_coords" not in st.session_state:
    st.session_state.veloviz_coords = None

if "veloviz_uploaded_file" not in st.session_state:
    st.session_state.veloviz_uploaded_file = None

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload h5ad file")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file with velocity data",
    type=['h5ad'],
    key="veloviz_h5ad_upload",
    help="scVelo-analyzed h5ad file (velocity layer required)"
)

if uploaded_h5ad is not None:
    # Check if new file uploaded
    if st.session_state.veloviz_uploaded_file != uploaded_h5ad.name:
        with st.spinner("Loading data..."):
            # Save and load file
            temp_h5ad_path = os.path.join(st.session_state.veloviz_temp_dir, "input.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            try:
                adata = sc.read_h5ad(temp_h5ad_path)
                st.session_state.veloviz_adata = adata
                st.session_state.veloviz_uploaded_file = uploaded_h5ad.name
                st.session_state.veloviz_coords = None  # Reset coords on new upload
                st.success(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")
            except Exception as e:
                st.error(f"❌ Failed to load file: {e}")
                st.stop()
    else:
        st.success(f"✓ File already loaded: {uploaded_h5ad.name}")

    adata = st.session_state.veloviz_adata

    # ========================================
    # Step 2: Data validation
    # ========================================
    st.header("Step 2: Data check")

    # Check required layers
    col1, col2, col3 = st.columns(3)

    with col1:
        has_velocity = 'velocity' in adata.layers
        if has_velocity:
            st.success("✓ velocity layer")
        else:
            st.error("✗ velocity layer missing")

    with col2:
        has_spliced = 'spliced' in adata.layers or 'Ms' in adata.layers
        spliced_key = 'Ms' if 'Ms' in adata.layers else 'spliced'
        if has_spliced:
            st.success(f"✓ {spliced_key} layer")
        else:
            st.error("✗ spliced/Ms layer missing")

    with col3:
        # Check for existing embeddings
        embeddings = [k for k in adata.obsm.keys() if k.startswith('X_')]
        if embeddings:
            st.success(f"✓ Embeddings: {', '.join(embeddings)}")
        else:
            st.warning("⚠ No existing embeddings")

    if not has_velocity or not has_spliced:
        st.error("""
        ❌ **Required data missing**

        VeloViz requires velocity and spliced (or Ms) layers.
        Please upload h5ad file analyzed by scVelo analysis app.
        """)
        st.stop()

    # Show data info
    st.info(f"""
    **Data summary:**
    - Cells: {adata.n_obs:,}
    - Genes: {adata.n_vars:,}
    - Layers: {', '.join(adata.layers.keys())}
    """)

    # ========================================
    # Step 3: Select columns
    # ========================================
    st.header("Step 3: Select columns for visualization")

    col1, col2 = st.columns(2)

    with col1:
        # Cluster column for coloring
        categorical_cols = [col for col in adata.obs.columns
                          if adata.obs[col].dtype.name == 'category' or adata.obs[col].nunique() < 50]

        default_cluster = 'leiden' if 'leiden' in categorical_cols else (
            'louvain' if 'louvain' in categorical_cols else (
                categorical_cols[0] if categorical_cols else None
            )
        )

        cluster_col = st.selectbox(
            "Cluster column for coloring:",
            options=categorical_cols,
            index=categorical_cols.index(default_cluster) if default_cluster in categorical_cols else 0,
            key="veloviz_cluster_col"
        )

    with col2:
        # Color palette
        palette_options = ['tab20', 'tab10', 'Set1', 'Set2', 'Set3', 'Paired', 'husl']
        color_palette = st.selectbox(
            "Color palette:",
            options=palette_options,
            index=0,
            key="veloviz_palette"
        )

    # ========================================
    # Step 4: VeloViz parameters
    # ========================================
    st.header("Step 4: VeloViz parameters")

    with st.expander("📊 Parameter settings", expanded=True):
        col1, col2 = st.columns(2)

        with col1:
            n_pcs = st.slider(
                "Number of PCs (nPCs):",
                min_value=5, max_value=100, value=20, step=5,
                help="Number of principal components. Higher is more detailed but increases computation cost"
            )

            k_neighbors = st.slider(
                "Number of neighbors (k):",
                min_value=5, max_value=100, value=30, step=5,
                help="Number of neighbors per cell. Larger gives smoother embedding"
            )

            similarity_threshold = st.slider(
                "Similarity threshold:",
                min_value=0.0, max_value=1.0, value=0.2, step=0.05,
                help="Similarity threshold with velocity direction. Higher is stricter"
            )

        with col2:
            distance_weight = st.slider(
                "Distance weight:",
                min_value=0.0, max_value=2.0, value=1.0, step=0.1,
                help="Weight of distance component"
            )

            normalize_depth = st.checkbox(
                "Normalize depth",
                value=True,
                help="Normalize by depth"
            )

            use_ods_genes = st.checkbox(
                "Use overdispersed genes",
                value=True,
                help="Use highly variable genes"
            )

            seed = st.number_input(
                "Random seed:",
                min_value=0, max_value=99999, value=42,
                help="Random seed for reproducibility"
            )

    # ========================================
    # Step 5: Run VeloViz
    # ========================================
    st.header("Step 5: Run VeloViz")

    # Estimate computation (rough)
    n_cells = adata.n_obs
    estimated_time = "< 1 min" if n_cells < 5000 else ("1-5 min" if n_cells < 20000 else "5+ min")
    st.info(f"⏱️ Estimated time: {estimated_time} ({n_cells:,} cells)")

    if st.button("🚀 Run VeloViz", type="primary", key="run_veloviz"):
        with st.spinner("Running VeloViz... This may take a few minutes."):
            try:
                # Prepare data
                # Get current expression (spliced or Ms)
                if 'Ms' in adata.layers:
                    curr = adata.layers['Ms']
                else:
                    curr = adata.layers['spliced']

                # Get velocity
                velocity = adata.layers['velocity']

                # Convert to dense if sparse
                if issparse(curr):
                    curr = curr.toarray()
                if issparse(velocity):
                    velocity = velocity.toarray()

                # Handle NaN in velocity
                velocity = np.nan_to_num(velocity, nan=0.0)

                # Calculate projected expression
                proj = curr + velocity

                # Ensure non-negative
                proj = np.clip(proj, 0, None)

                # Get gene and cell names
                gene_names = adata.var_names.tolist()
                cell_names = adata.obs_names.tolist()

                st.info(f"Data prepared: curr shape = {curr.shape}, proj shape = {proj.shape}, genes = {len(gene_names)}, cells = {len(cell_names)}")

                # Run VeloViz
                coords = run_veloviz(
                    curr_matrix=curr.astype(np.float64),
                    proj_matrix=proj.astype(np.float64),
                    gene_names=gene_names,
                    cell_names=cell_names,
                    n_pcs=n_pcs,
                    k=k_neighbors,
                    similarity_threshold=similarity_threshold,
                    distance_weight=distance_weight,
                    normalize_depth=normalize_depth,
                    use_ods_genes=use_ods_genes,
                    seed=int(seed)
                )

                # Store results
                st.session_state.veloviz_coords = coords

                # Add to adata
                adata.obsm['X_veloviz'] = coords
                st.session_state.veloviz_adata = adata

                st.success(f"✓ VeloViz completed! Embedding shape: {coords.shape}")

            except Exception as e:
                st.error(f"""
                ❌ **VeloViz failed**

                Error: {str(e)}

                **Troubleshooting:**
                1. Check if data has many NaN values
                2. Check if cell count is too low (recommended: >500 cells)
                3. Check if velocity calculation was performed correctly
                """)
                import traceback
                st.code(traceback.format_exc())

    # ========================================
    # Step 6: Visualization
    # ========================================
    if st.session_state.veloviz_coords is not None:
        st.header("Step 6: Visualization")

        coords = st.session_state.veloviz_coords
        adata = st.session_state.veloviz_adata

        # Create color mapping
        color_map = create_color_mapping(adata.obs[cluster_col], color_palette)
        colors = [color_map[cat] for cat in adata.obs[cluster_col]]

        # Plot options
        col1, col2 = st.columns(2)
        with col1:
            point_size = st.slider("Point size:", 1, 50, 10, key="veloviz_point_size")
        with col2:
            alpha = st.slider("Transparency:", 0.1, 1.0, 0.7, key="veloviz_alpha")

        # Create figure (VeloViz only)
        fig, ax = plt.subplots(figsize=(10, 8))

        # VeloViz plot
        scatter = ax.scatter(
            coords[:, 0], coords[:, 1],
            c=colors, s=point_size, alpha=alpha
        )
        ax.set_xlabel("VeloViz 1")
        ax.set_ylabel("VeloViz 2")
        ax.set_title(f"VeloViz Embedding (colored by {cluster_col})")
        ax.set_aspect('equal', 'box')

        # Add legend
        unique_cats = sorted(adata.obs[cluster_col].unique())
        legend_handles = [plt.scatter([], [], c=[color_map[cat]], label=cat, s=50) for cat in unique_cats]
        fig.legend(handles=legend_handles, loc='center right', bbox_to_anchor=(1.15, 0.5), title=cluster_col)

        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

        # ========================================
        # Optional: Velocity stream on VeloViz
        # ========================================
        st.subheader("Velocity stream on VeloViz embedding")

        if st.checkbox("Show velocity stream", key="veloviz_show_stream"):
            with st.spinner("Computing velocity stream..."):
                try:
                    # Temporarily set the embedding
                    fig_stream, ax_stream = plt.subplots(figsize=(10, 8))

                    # Use scvelo to plot velocity stream on veloviz embedding
                    scv.pl.velocity_embedding_stream(
                        adata,
                        basis='veloviz',
                        color=cluster_col,
                        ax=ax_stream,
                        legend_loc='right margin',
                        show=False
                    )

                    st.pyplot(fig_stream)
                    plt.close()

                except Exception as e:
                    st.warning(f"Could not compute velocity stream: {e}")

        # ========================================
        # Step 7: Download
        # ========================================
        st.header("Step 7: Download results")

        col1, col2 = st.columns(2)

        with col1:
            # Download h5ad with embedding
            output_path = os.path.join(st.session_state.veloviz_temp_dir, "output_veloviz.h5ad")
            adata.write_h5ad(output_path)

            with open(output_path, "rb") as f:
                st.download_button(
                    label="📥 Download h5ad with VeloViz embedding",
                    data=f.read(),
                    file_name=f"{uploaded_h5ad.name.replace('.h5ad', '')}_VeloViz.h5ad",
                    mime="application/octet-stream",
                    key="download_veloviz_h5ad"
                )

        with col2:
            # Download coordinates as CSV
            coords_df = pd.DataFrame(
                coords,
                columns=['VeloViz_1', 'VeloViz_2'],
                index=adata.obs_names
            )
            coords_df[cluster_col] = adata.obs[cluster_col].values

            csv_buffer = io.StringIO()
            coords_df.to_csv(csv_buffer)

            st.download_button(
                label="📥 Download coordinates (CSV)",
                data=csv_buffer.getvalue(),
                file_name="coordinates_VeloViz.csv",
                mime="text/csv",
                key="download_veloviz_csv"
            )

        # Download plot
        st.subheader("Download plot")

        col1, col2 = st.columns(2)
        with col1:
            plot_format = st.selectbox("Format:", ["PDF", "PNG", "SVG"], key="veloviz_plot_format")
        with col2:
            plot_dpi = st.slider("DPI (for PNG):", 100, 600, 300, key="veloviz_plot_dpi")

        # Regenerate plot for download (VeloViz only)
        fig_download, ax_download = plt.subplots(figsize=(10, 8))

        ax_download.scatter(coords[:, 0], coords[:, 1], c=colors, s=point_size, alpha=alpha)
        ax_download.set_xlabel("VeloViz 1")
        ax_download.set_ylabel("VeloViz 2")
        ax_download.set_title(f"VeloViz Embedding (colored by {cluster_col})")
        ax_download.set_aspect('equal', 'box')

        legend_handles = [plt.scatter([], [], c=[color_map[cat]], label=cat, s=50) for cat in unique_cats]
        fig_download.legend(handles=legend_handles, loc='center right', bbox_to_anchor=(1.15, 0.5), title=cluster_col)
        plt.tight_layout()

        # Save to buffer
        buf = io.BytesIO()
        fig_download.savefig(buf, format=plot_format.lower(), dpi=plot_dpi, bbox_inches='tight')
        buf.seek(0)
        plt.close(fig_download)

        st.download_button(
            label=f"📥 Download plot ({plot_format})",
            data=buf.getvalue(),
            file_name=f"plot_VeloViz.{plot_format.lower()}",
            mime=f"image/{plot_format.lower()}" if plot_format != "PDF" else "application/pdf",
            key="download_veloviz_plot"
        )

else:
    st.info("👆 Upload an h5ad file with velocity data to get started.")
