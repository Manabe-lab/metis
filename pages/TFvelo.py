"""
TFvelo Web Interface
Copyright (C) 2024

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

This is a Streamlit web application for TFvelo analysis.
Original TFvelo: https://github.com/xiaoyeye/TFvelo
Citation: Li, Jiachen, et al. "TFvelo: gene regulation inspired RNA velocity estimation." bioRxiv (2023)
"""

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import tempfile
import shutil
from datetime import datetime
import traceback
import subprocess
import sys
import io
# from contextlib import redirect_stdout, redirect_stderr  # Not needed

# TFvelo specific imports
try:
    import scvelo as scv
    import anndata as ad
    from scipy import sparse
except ImportError as e:
    st.error(f"Required package not installed: {e}")
    st.stop()

# Import TFvelo with relative path
import sys
import os
current_dir = os.path.dirname(os.path.abspath(__file__))
tfvelo_path = os.path.join(os.path.dirname(current_dir), 'TFvelo')
sys.path.append(tfvelo_path)
try:
    import TFvelo as TFv
    TFVELO_AVAILABLE = True
    st.success("✅ TFvelo package loaded successfully")
except ImportError as e:
    TFVELO_AVAILABLE = False
    st.error(f"❌ TFvelo package not available: {e}")

# Import additional required packages for demo functions
import scipy.stats

# TFvelo demo missing functions (from TFvelo_demo.ipynb)
def get_sort_positions(arr):
    """Helper function for get_sort_t - normalizes array positions"""
    positions = np.argsort(np.argsort(arr))
    positions_normed = positions/(len(arr)-1)
    return positions_normed

def get_sort_t(adata):
    """
    Advanced temporal analysis function from TFvelo demo
    Creates latent_t layer and non_blank_gene annotations
    """
    t = adata.layers['fit_t_raw'].copy()
    normed_t = adata.layers['fit_t_raw'].copy()
    n_bins = 20
    n_cells, n_genes = adata.shape
    sort_t = np.zeros([n_cells, n_genes])
    non_blank_gene = np.zeros(n_genes, dtype=int)
    hist_all, bins_all = np.zeros([n_genes, n_bins]), np.zeros([n_genes, n_bins+1])

    for i in range(n_genes):
        gene_name = adata.var_names[i]
        tmp = t[:,i].copy()
        if np.isnan(tmp).sum():
            non_blank_gene[i] = 1
            continue
        hist, bins = np.histogram(tmp, bins=n_bins)
        hist_all[i], bins_all[i] = hist, bins
        if not (0 in list(hist)):
            if (tmp.min() < 0.1) and (tmp.max() > 0.8):
                blank_start_bin_id = np.argmin(hist)
                blank_end_bin_id = blank_start_bin_id
                non_blank_gene[i] = 1
                blank_start_bin = bins[blank_start_bin_id]
                blank_end_bin = bins[blank_end_bin_id]
                tmp = (tmp < blank_start_bin)*1 + tmp
            else:
                blank_end_bin = tmp.min()
        else:
            blank_start_bin_id = list(hist).index(0)
            for j in range(blank_start_bin_id+1, len(hist)):
                if hist[j] > 0:
                    blank_end_bin_id = j
                    break
            blank_start_bin = bins[blank_start_bin_id]
            blank_end_bin = bins[blank_end_bin_id]
            tmp = (tmp < blank_start_bin)*1 + tmp

        t[:,i] = tmp
        tmp = tmp - blank_end_bin
        tmp = tmp/tmp.max()
        normed_t[:,i] = tmp
        sort_t[:,i] = get_sort_positions(tmp)

    adata.layers['latent_t'] = sort_t.copy()
    adata.var['non_blank_gene'] = non_blank_gene.copy()
    return adata

def get_metric_pseudotime(adata, t_key='latent_t'):
    """
    Calculate Spearman correlation between latent time and pseudotime
    """
    n_cells, n_genes = adata.shape
    adata.var['spearmanr_pseudotime'] = 0.0
    for i in range(n_genes):
        correlation, _ = scipy.stats.spearmanr(adata.layers[t_key][:,i], adata.obs['velocity_pseudotime'])
        adata.var['spearmanr_pseudotime'][i] = correlation
    return adata

def apply_demo_gene_filtering(adata, loss_percent_thres=50, spearmanr_thres=0.8):
    """
    Apply the complete TFvelo demo gene filtering pipeline
    """
    adata_filtered = adata.copy()

    # Filter 1: Loss percentile threshold (50th percentile from demo)
    if 'min_loss' in adata_filtered.var:
        thres_loss = np.percentile(adata_filtered.var['min_loss'], loss_percent_thres)
        n_genes_before = adata_filtered.n_vars
        adata_filtered = adata_filtered[:, adata_filtered.var['min_loss'] < thres_loss]
        print(f"Loss filter ({loss_percent_thres}th percentile): {n_genes_before} → {adata_filtered.n_vars} genes")

    # Filter 2: Cell count threshold (10% of total cells from demo)
    if 'n_cells' in adata_filtered.var:
        thres_n_cells = adata_filtered.X.shape[0] * 0.1
        n_genes_before = adata_filtered.n_vars
        adata_filtered = adata_filtered[:, adata_filtered.var['n_cells'] > thres_n_cells]
        print(f"Cell count filter (10%): {n_genes_before} → {adata_filtered.n_vars} genes")

    # Filter 3: Non-blank gene filter
    if 'non_blank_gene' in adata_filtered.var:
        n_genes_before = adata_filtered.n_vars
        adata_filtered = adata_filtered[:, adata_filtered.var['non_blank_gene'] == 0]
        print(f"Non-blank gene filter: {n_genes_before} → {adata_filtered.n_vars} genes")

    # Filter 4: Spearman correlation filter (0.8 threshold from demo)
    if 'spearmanr_pseudotime' in adata_filtered.var:
        n_genes_before = adata_filtered.n_vars
        adata_filtered = adata_filtered[:, adata_filtered.var['spearmanr_pseudotime'] > spearmanr_thres]
        print(f"Spearman correlation filter (>{spearmanr_thres}): {n_genes_before} → {adata_filtered.n_vars} genes")

    return adata_filtered

# Page configuration
st.set_page_config(
    page_title="TFvelo Analysis",
    page_icon="🧬",
    layout="wide"
)

# Title and description
st.title("🧬 TFvelo Analysis")
st.markdown("""
**Transcription Factor-based Gene Regulatory Network Analysis**

TFvelo analyzes single-cell dynamics using transcription factor regulatory networks without requiring RNA velocity (spliced/unspliced) data.
Based on the method from [Qiu et al., Nature Methods 2022](https://github.com/xiaoyeye/TFvelo)
""")

# Initialize session state
if 'adata' not in st.session_state:
    st.session_state.adata = None
if 'tfvelo_results' not in st.session_state:
    st.session_state.tfvelo_results = {}
if 'temp_dir' not in st.session_state:
    st.session_state.temp_dir = None

# Data Input Section
st.header("📁 Data Input")

# Data source selection (moved from visualization section)
data_source = st.radio(
    "Choose data source:",
    ["Upload new data", "Load saved TFvelo analysis for visualization"],
    help="Upload a new h5ad file for analysis or load previously saved TFvelo results for visualization"
)

if data_source == "Upload new data":
    # File upload for new analysis
    uploaded_file = st.file_uploader(
        "Upload h5ad file",
        type=['h5ad'],
        help="Upload a preprocessed AnnData object with RNA velocity information"
    )
elif data_source == "Load saved TFvelo analysis for visualization":
    # Load saved analysis option - upload only
    st.subheader("📤 Upload TFvelo Analysis File")
    uploaded_file = st.file_uploader(
        "Choose a TFvelo analysis h5ad file",
        type=['h5ad'],
        help="Upload a TFvelo analysis result file (.h5ad format)",
        key="tfvelo_upload"
    )

    if uploaded_file is not None:
        try:
            # Save uploaded file temporarily
            temp_dir = "/tmp/tfvelo_uploads"
            os.makedirs(temp_dir, exist_ok=True)
            temp_path = os.path.join(temp_dir, uploaded_file.name)

            with open(temp_path, "wb") as f:
                f.write(uploaded_file.getbuffer())

            # Load the uploaded file
            import anndata as ad
            adata_tfvelo = ad.read_h5ad(temp_path)

            # Store in session state for visualization
            st.session_state.adata_tfvelo = adata_tfvelo
            st.session_state.tfvelo_results = {'completed': True}
            st.session_state.data_source = "saved_analysis"

            # Also store as regular adata for tab display
            st.session_state.adata = adata_tfvelo

            st.success(f"✅ Loaded TFvelo analysis: {uploaded_file.name}")
            st.info(f"Data shape: {adata_tfvelo.shape}")

            # Show basic info about uploaded data
            with st.expander("📊 Analysis Data Summary", expanded=True):
                col_a, col_b = st.columns(2)
                with col_a:
                    st.metric("Cells", adata_tfvelo.n_obs)
                    st.metric("Genes", adata_tfvelo.n_vars)
                with col_b:
                    if 'velocity' in adata_tfvelo.layers:
                        st.metric("Velocity Layer", "✅ Available")
                    if 'n_TFs' in adata_tfvelo.var.columns:
                        st.metric("Total TFs", len(adata_tfvelo.uns.get('all_TFs', [])))

            # Redirect to Visualization tab instead of showing duplicate visualization
            st.success("✅ TFvelo analysis loaded successfully!")
            st.info("💡 **Scroll down** to see the analysis tabs ('📊 Data Overview', '📈 Visualization', etc.).")

            # Skip duplicate visualization - use Visualization tab instead
            if False:  # Disabled duplicate visualization
                col1, col2 = st.columns(2)

                with col1:
                    selected_embedding = st.selectbox(
                        "Select embedding for projection",
                        available_embeddings,
                        format_func=lambda x: x.replace('X_', '').upper(),
                        help="Choose which dimensionality reduction to use for visualization"
                    )

                with col2:
                    # Select coloring variable
                    color_options = []

                    # Add TFvelo-specific metrics
                    if 'velocity_pseudotime' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_pseudotime')
                    if 'velocity_length' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_length')
                    if 'velocity_confidence' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_confidence')

                    # Add cluster information
                    if 'clusters' in adata_tfvelo.obs.columns:
                        color_options.append('clusters')

                    # Add other metadata
                    for col in adata_tfvelo.obs.columns:
                        if col not in color_options:
                            color_options.append(col)

                    color_by = st.selectbox(
                        "Color by",
                        color_options,
                        index=0 if 'velocity_pseudotime' in color_options else 0,
                        help="Choose variable for coloring cells"
                    )

                # Visualization options
                st.subheader("Visualization Options")

                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    arrow_length = st.slider("Arrow length", 1, 10, 4, key="upload_arrow_length")
                with col2:
                    arrow_size = st.slider("Arrow size", 0.5, 3.0, 1.0, key="upload_arrow_size")
                with col3:
                    density = st.slider("Stream density", 0.5, 4.0, 2.0, key="upload_density")
                with col4:
                    point_size = st.slider("Cell size", 10, 1000, 100, key="upload_point_size")

                # Advanced options
                with st.expander("Advanced Visualization Options"):
                    col1, col2, col3, col4 = st.columns(4)

                    with col1:
                        alpha_value = st.slider("Cell transparency", 0.1, 1.0, 0.7, 0.1, key="upload_alpha")

                    with col2:
                        plot_height = st.slider("Plot height", 4, 12, 6, key="upload_height")

                    with col3:
                        # Determine if selected color is categorical or continuous
                        is_categorical = False
                        if color_by in adata_tfvelo.obs.columns:
                            dtype = adata_tfvelo.obs[color_by].dtype
                            is_categorical = dtype == 'object' or dtype.name == 'category'

                        # Colormap selection based on data type
                        if is_categorical:
                            colormap = st.selectbox("Colormap",
                                                   ['tab20', 'tab10', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Paired', 'Accent'],
                                                   index=0,
                                                   key="upload_colormap")
                        else:
                            colormap = st.selectbox("Colormap",
                                                   ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'coolwarm', 'RdBu_r', 'seismic'],
                                                   index=0,
                                                   key="upload_colormap2")

                    with col4:
                        dpi = st.slider("DPI (resolution)", 50, 300, 100, key="upload_dpi")

                # Check if velocity information exists
                has_velocity = 'velocity' in adata_tfvelo.layers or 'velocities' in adata_tfvelo.layers

                if not has_velocity:
                    st.warning("⚠️ No velocity information found in the uploaded data.")
                    st.info("The uploaded file should contain velocity calculations. Available layers: " +
                           str(list(adata_tfvelo.layers.keys()) if adata_tfvelo.layers else "None"))

                    # Still show basic embedding plot without velocity
                    st.subheader("Data Projection (without velocity)")

                    import matplotlib.pyplot as plt
                    import scanpy as sc

                    fig, ax = plt.subplots(1, 1, figsize=(plot_height, plot_height))

                    # Use scanpy for basic plotting
                    sc.pl.embedding(
                        adata_tfvelo,
                        basis=selected_embedding.replace('X_', ''),
                        color=color_by,
                        ax=ax,
                        size=point_size,
                        alpha=alpha_value,
                        palette=colormap if is_categorical else None,
                        cmap=colormap if not is_categorical else None,
                        show=False,
                        legend_loc='right margin' if is_categorical else 'none',
                        title=f"{selected_embedding.replace('X_', '').upper()} - {color_by}"
                    )
                    ax.grid(False)

                    plt.tight_layout()
                    st.pyplot(fig, dpi=dpi)

                else:
                    # Create visualizations with velocity
                    st.subheader("Velocity Projections")

                    # Prepare figure
                    import matplotlib.pyplot as plt
                    import scvelo as scv

                    fig, axes = plt.subplots(1, 2, figsize=(plot_height*2, plot_height))

                    # Left plot: Velocity stream with manual colormap application
                    try:
                        import numpy as np

                        # Get embedding coordinates
                        embedding_key = selected_embedding
                        X_emb = adata_tfvelo.obsm[embedding_key]

                        # Get color data
                        if color_by in adata_tfvelo.obs.columns:
                            color_data = adata_tfvelo.obs[color_by].values

                            # Handle categorical vs continuous data
                            if is_categorical:
                                # Use categorical colors
                                unique_vals = np.unique(color_data)
                                colors = plt.cm.get_cmap(colormap, len(unique_vals))
                                color_map = {val: colors(i) for i, val in enumerate(unique_vals)}
                                point_colors = [color_map[val] for val in color_data]
                            else:
                                # Use continuous colormap
                                norm = plt.Normalize(vmin=np.min(color_data), vmax=np.max(color_data))
                                cmap = plt.cm.get_cmap(colormap)
                                point_colors = cmap(norm(color_data))
                        else:
                            point_colors = 'gray'

                        # Plot scatter points with proper colors
                        scatter = axes[0].scatter(
                            X_emb[:, 0], X_emb[:, 1],
                            c=point_colors,
                            s=point_size,
                            alpha=alpha_value,
                            edgecolors='none'
                        )

                        # Then overlay velocity stream (arrows only)
                        try:
                            # Clear the axis to avoid overlap issues
                            axes[0].clear()

                            # Re-plot scatter with colors
                            axes[0].scatter(
                                X_emb[:, 0], X_emb[:, 1],
                                c=point_colors,
                                s=point_size,
                                alpha=alpha_value,
                                edgecolors='none'
                            )

                            # Try to overlay velocity stream
                            scv.pl.velocity_embedding_stream(
                                adata_tfvelo,
                                basis=selected_embedding.replace('X_', ''),
                                ax=axes[0],
                                color=color_by,  # Try with color specified
                                arrow_length=arrow_length,
                                arrow_size=arrow_size,
                                density=density,
                                show=False,
                                linewidth=1.5,
                                alpha=0.6
                            )
                        except:
                            # If overlay fails, just keep the colored scatter plot
                            pass

                        axes[0].set_title(f"Velocity Stream ({selected_embedding.replace('X_', '').upper()})")
                        axes[0].grid(False)
                    except Exception as e:
                        st.warning(f"Could not create velocity stream plot: {str(e)}")
                        # Fallback to basic plot
                        sc.pl.embedding(
                            adata_tfvelo,
                            basis=selected_embedding.replace('X_', ''),
                            color=color_by,
                            ax=axes[0],
                            size=point_size,
                            alpha=alpha_value,
                            palette=colormap if is_categorical else None,
                            cmap=colormap if not is_categorical else None,
                            show=False,
                            legend_loc='right margin' if is_categorical else 'none',
                            title=f"{selected_embedding.replace('X_', '').upper()} - {color_by}"
                        )
                        axes[0].grid(False)

                    # Right plot: Velocity arrows with colormap
                    try:
                        scv.pl.velocity_embedding(
                            adata_tfvelo,
                            basis=selected_embedding.replace('X_', ''),
                            ax=axes[1],
                            color=color_by,
                            size=point_size,
                            alpha=alpha_value,
                            arrow_length=arrow_length*2,
                            arrow_size=arrow_size*1.5,
                            palette=colormap if is_categorical else None,
                            cmap=colormap if not is_categorical else None,
                            show=False,
                            legend_loc='right margin' if is_categorical else 'right',
                            title=f"Velocity Pseudotime ({selected_embedding.replace('X_', '').upper()})"
                        )
                        axes[1].grid(False)
                    except Exception as e:
                        # Fallback to basic plot
                        sc.pl.embedding(
                            adata_tfvelo,
                            basis=selected_embedding.replace('X_', ''),
                            color=color_by,
                            ax=axes[1],
                            size=point_size,
                            alpha=alpha_value,
                            palette=colormap if is_categorical else None,
                            cmap=colormap if not is_categorical else None,
                            show=False,
                            legend_loc='none',
                            title=f"{selected_embedding.replace('X_', '').upper()} (no velocity)"
                        )
                        axes[1].grid(False)

                    plt.tight_layout()
                    st.pyplot(fig, dpi=dpi)

                # Download options for plots
                col1, col2 = st.columns(2)
                with col1:
                    if st.button("💾 Download Plot as PNG", key="upload_download_png"):
                        fig.savefig('/tmp/tfvelo_plot.png', dpi=300, bbox_inches='tight')
                        with open('/tmp/tfvelo_plot.png', 'rb') as f:
                            st.download_button(
                                label="📥 Download PNG",
                                data=f.read(),
                                file_name="tfvelo_visualization.png",
                                mime='image/png',
                                key="upload_dl_png"
                            )

                with col2:
                    if st.button("💾 Download Plot as PDF", key="upload_download_pdf"):
                        try:
                            fig.savefig('/tmp/tfvelo_plot.pdf', format='pdf', bbox_inches='tight')
                            with open('/tmp/tfvelo_plot.pdf', 'rb') as f:
                                st.download_button(
                                    label="📥 Download PDF",
                                    data=f.read(),
                                    file_name="tfvelo_visualization.pdf",
                                    mime='application/pdf',
                                    key="upload_dl_pdf"
                                )
                        except Exception as e:
                            st.error(f"PDF export failed: {str(e)}")
            # else block removed - was incorrectly always running due to "if False:"

        except Exception as e:
            st.error(f"Error loading TFvelo analysis file: {str(e)}")

# Continue with new data upload processing only if "Upload new data" is selected
uploaded_file = None if data_source == "Load saved TFvelo analysis for visualization" else uploaded_file

if uploaded_file is not None:
    try:
        # Create temporary directory
        if st.session_state.temp_dir is None:
            st.session_state.temp_dir = tempfile.mkdtemp(prefix="tfvelo_")

        # Save uploaded file
        temp_path = os.path.join(st.session_state.temp_dir, uploaded_file.name)
        with open(temp_path, 'wb') as f:
            f.write(uploaded_file.getbuffer())

        # Load AnnData object
        with st.spinner("Loading h5ad file..."):
            adata = sc.read_h5ad(temp_path)
            st.session_state.adata = adata

        col1, col2, col3 = st.columns(3)
        with col1:
            st.success(f"✅ Loaded data: {adata.n_obs} cells × {adata.n_vars} genes")

        with col2:
            # Display available embeddings
            if adata.obsm:
                st.info("**Available embeddings:**")
                for key in adata.obsm.keys():
                    st.text(f"  • {key}")

        with col3:
            # Display available layers with more details in expandable section
            with st.expander("📋 Available layers", expanded=False):
                if adata.layers:
                    for key in adata.layers.keys():
                        layer_shape = adata.layers[key].shape
                        layer_type = type(adata.layers[key]).__name__
                        st.text(f"• {key} ({layer_type}, {layer_shape})")
                else:
                    st.text("• None")

            # Also show if raw data exists
            if hasattr(adata, 'raw') and adata.raw is not None:
                raw_shape = adata.raw.X.shape
                raw_type = type(adata.raw.X).__name__
                st.text(f"  • raw.X ({raw_type}, {raw_shape})")

        # Data layer selection
        st.markdown("### 📊 Data Layer Selection")

        # Add information about common layer types
        with st.expander("ℹ️ Common layer types", expanded=False):
            st.markdown("""
            **Common h5ad layers:**
            - `X` - Main expression matrix (usually normalized)
            - `raw.X` - Raw count matrix (if preserved)
            - `counts` - Raw counts
            - `logcounts` - Log-transformed counts
            - `scaled` - Scaled/standardized data
            - `pearson_residuals` - Pearson residuals (SCTransform)
            - `normalized` - Normalized expression values
            - `batch_corrected` - Batch-corrected expression
            """)


        # Create options for data source
        data_options = ["adata.X (current)"]
        if hasattr(adata, 'raw') and adata.raw is not None:
            data_options.append("adata.raw.X (raw counts)")
        if adata.layers:
            for layer_name in adata.layers.keys():
                layer_shape = adata.layers[layer_name].shape
                data_options.append(f"adata.layers['{layer_name}'] ({layer_shape})")

        selected_data_source = st.selectbox(
            "Select data source:",
            data_options,
            index=0,
            help="Choose which data to use for TFvelo analysis"
        )

        # Check if spliced/unspliced layers are selected
        is_rna_velocity_mode = False
        if "spliced" in selected_data_source.lower() or "unspliced" in selected_data_source.lower():
            is_rna_velocity_mode = True
            st.info("🧬 **RNA Velocity mode detected**: Using spliced/unspliced data (paper demo compatible)")
            data_type = "Raw counts (apply full normalization)"  # Auto-set for RNA velocity
        else:
            # Data type selection (only show when not in RNA velocity mode)
            data_type = st.radio(
                "Data type:",
                ["Raw counts (apply full normalization)", "Normalized data (skip log transform)"],
                index=0,
                help="Select whether your data is raw counts or already normalized"
            )

        # Data source and type validation warnings (immediate feedback)
        try:
            # Get preview data for validation
            if selected_data_source == "adata.X (current)":
                validation_data = adata.X
            elif selected_data_source == "adata.raw.X (raw counts)":
                validation_data = adata.raw.X
            else:
                # Extract layer name from the string
                if "adata.layers['" in selected_data_source:
                    if "] (" in selected_data_source:
                        layer_name = selected_data_source.split("'")[1].split("'")[0]
                    else:
                        layer_name = selected_data_source.split("'")[1]
                    validation_data = adata.layers[layer_name]
                else:
                    validation_data = adata.X

            # Sample data for quick validation
            if hasattr(validation_data, 'toarray'):
                sample_flat = validation_data[:100, :100].toarray().flatten() if validation_data.shape[0] > 100 else validation_data.toarray().flatten()
            else:
                sample_flat = validation_data[:100, :100].flatten() if validation_data.shape[0] > 100 else validation_data.flatten()

            non_zero_sample = sample_flat[sample_flat > 0]
            if len(non_zero_sample) > 0:
                has_decimals = np.any(non_zero_sample % 1 != 0)
                max_val = np.max(non_zero_sample)

                # Adjust for radio button text difference
                normalized_text = "Normalized data (skip log transform)"

                # Show validation warnings
                if has_decimals and data_type == "Raw counts (apply full normalization)":
                    st.warning("⚠️ Data contains decimal values but is marked as raw counts. Consider selecting 'Normalized data'.")

                if not has_decimals and max_val > 50 and data_type == normalized_text:
                    st.warning("⚠️ Data appears to be raw counts (integers, high values) but is marked as normalized. Consider selecting 'Raw counts'.")

                if has_decimals and max_val < 20 and data_type == "Raw counts (apply full normalization)":
                    st.warning("⚠️ Data appears to be log-transformed (decimals, low values) but is marked as raw counts. Consider selecting 'Normalized data'.")

                if selected_data_source == "adata.raw.X (raw counts)" and data_type == normalized_text:
                    st.warning("⚠️ Selected raw.X as data source but marked as normalized data. This combination may be inconsistent.")

                if "layers" in selected_data_source and "raw" in selected_data_source.lower() and data_type == normalized_text:
                    st.warning("⚠️ Selected a 'raw' layer but marked as normalized data. Check if this combination is intended.")

        except Exception:
            pass  # Skip validation if data cannot be accessed

        # Display data preview (always shown)
        st.subheader("Data Preview")
        if selected_data_source == "adata.X (current)":
            preview_data = adata.X
        elif selected_data_source == "adata.raw.X (raw counts)":
            preview_data = adata.raw.X
        else:
            # Extract layer name from the string (handle both old and new format)
            if "adata.layers['" in selected_data_source:
                if "] (" in selected_data_source:
                    # New format with shape info: "adata.layers['layer_name'] (shape)"
                    layer_name = selected_data_source.split("'")[1].split("'")[0]
                else:
                    # Old format: "adata.layers['layer_name']"
                    layer_name = selected_data_source.split("'")[1]
                preview_data = adata.layers[layer_name]
            else:
                st.error(f"Cannot parse data source: {selected_data_source}")
                preview_data = adata.X

        # Create preview DataFrame
        if hasattr(preview_data, 'toarray'):
            preview_array = preview_data[:5, :8].toarray()
        else:
            preview_array = preview_data[:5, :8]

        preview_df = pd.DataFrame(
            preview_array,
            index=adata.obs_names[:5],
            columns=adata.var_names[:8]
        )
        st.dataframe(preview_df)

        # Display data characteristics (informational only)
        if hasattr(preview_data, 'toarray'):
            data_flat = preview_data.data if hasattr(preview_data, 'data') else preview_data.toarray().flatten()
        else:
            data_flat = preview_data.flatten()

        non_zero_data = data_flat[data_flat > 0]
        if len(non_zero_data) > 0:
            has_decimals = np.any(non_zero_data % 1 != 0)
            max_val = np.max(non_zero_data)
            st.write(f"🔍 **Data characteristics:** Max value: {max_val:.2f}, Contains decimals: {has_decimals}")

    except Exception as e:
        st.error(f"Error loading file: {str(e)}")
        st.session_state.adata = None

st.divider()

# Main content area
if st.session_state.adata is not None:
    adata = st.session_state.adata

    # Create tabs for different analyses
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "📊 Data Overview",
        "⚙️ TFvelo Setup & Analysis",
        "🔬 Gene Filtering",
        "📈 Visualization",
        "💾 Download Results"
    ])

    # Tab 1: Data Overview
    with tab1:
        st.header("Data Overview")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Number of cells", adata.n_obs)
        with col2:
            st.metric("Number of genes", adata.n_vars)
        with col3:
            st.metric("Number of observations",
                     len(adata.obs.columns) if adata.obs is not None else 0)

        # Display metadata
        if adata.obs is not None and len(adata.obs.columns) > 0:
            st.subheader("Cell Metadata")
            st.dataframe(adata.obs.head(10))

        # Display gene information
        if adata.var is not None and len(adata.var.columns) > 0:
            st.subheader("Gene Information")
            st.dataframe(adata.var.head(10))

    # Tab 2: TFvelo Setup
    with tab2:
        st.header("TFvelo Configuration")

        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Basic Parameters")

            # Variable genes selection
            n_top_genes = st.number_input(
                "Number of variable genes",
                min_value=500,
                max_value=5000,
                value=2000,
                step=500,
                help="Number of most variable genes to use for analysis. More genes = more TFs available but slower computation (2000 genes = ~20min, 4000 genes = ~40min)"
            )

            # Number of CPUs
            n_cpus = st.number_input(
                "Number of CPUs",
                min_value=1,
                max_value=os.cpu_count(),
                value=min(4, os.cpu_count()),
                help="Number of CPUs to use for parallel processing"
            )

            # Neighbor selection
            n_neighbors = st.number_input(
                "Number of neighbors",
                min_value=5,
                max_value=100,
                value=30,
                help="Number of neighbors for KNN graph"
            )

            # TF database selection (multiple)
            tf_databases = st.multiselect(
                "TF Databases",
                ["ENCODE", "ChEA"],
                default=["ENCODE", "ChEA"],
                help="Select transcription factor databases to combine"
            )

            # General TF filtering parameters
            st.markdown("**General TF Filtering**")
            min_expression = st.number_input(
                "Min expression threshold",
                min_value=0.0,
                max_value=1.0,
                value=0.1,
                step=0.01,
                help="Minimum expression for TF-target correlation (default: 0.1)"
            )
            min_cells = st.number_input(
                "Min cells (general)",
                min_value=1,
                max_value=50,
                value=2,
                help="Minimum cells where both TF and target are expressed (default: 2)"
            )

            # Important TF options
            st.markdown("**Important TFs (Optional)**")
            force_include_tfs_input = st.text_area(
                "Force include TFs",
                value="",
                help="Enter important TF names (one per line or comma-separated). These TFs use separate filtering thresholds.",
                placeholder="e.g., FOXP3\nBCL6\nGATA3"
            )

            # Parse force_include_tfs
            if force_include_tfs_input.strip():
                # Handle both comma and newline separated
                force_include_tfs = [tf.strip() for tf in force_include_tfs_input.replace(',', '\n').split('\n') if tf.strip()]
            else:
                force_include_tfs = None

            # Show filtering options for important TFs only if some are specified
            if force_include_tfs:
                st.info(f"Important TFs: {', '.join(force_include_tfs)}")
                min_expression_force = st.number_input(
                    "Min expression (important TFs)",
                    min_value=0.0,
                    max_value=1.0,
                    value=0.0,
                    step=0.01,
                    help="Expression threshold for important TFs (default: 0.0 = any expression)"
                )
                min_cells_force = st.number_input(
                    "Min cells (important TFs)",
                    min_value=1,
                    max_value=50,
                    value=1,
                    help="Minimum cells for important TF-target pairs (default: 1)"
                )
            else:
                min_expression_force = 0.0
                min_cells_force = 1

        with col2:
            st.subheader("TFvelo Parameters")

            # TFvelo demo parameters
            # Memory-aware core usage based on dataset size and available memory
            try:
                available_mem_gb = int(os.popen("free -g | awk '/^Mem:/{print $7}'").read().strip())
            except:
                available_mem_gb = 64  # Default assumption

            # Estimate memory per worker based on adata size
            n_cells = adata.n_obs
            n_genes = adata.n_vars
            # Each worker copies adata: ~5 layers × n_cells × n_genes × 8 bytes
            # Plus varm matrices: ~5 × n_genes × n_genes × 8 bytes (before cleanup)
            # Plus overhead
            adata_mem_gb = (n_cells * n_genes * 5 * 8 + n_genes * n_genes * 5 * 8) / 1e9
            worker_mem_gb = adata_mem_gb + 1  # +1GB overhead per worker

            # Reserve memory for main process (result matrices + session state)
            main_process_gb = (n_cells * n_genes * 8 * 4 / 1e9) + 10  # 8 result matrices (float32) + 10GB buffer

            # Calculate optimal cores
            usable_mem_gb = available_mem_gb - main_process_gb
            mem_based_cores = max(1, int(usable_mem_gb / worker_mem_gb))
            optimal_cores = min(32, mem_based_cores, os.cpu_count() // 2)

            n_jobs = st.number_input(
                "Number of jobs",
                min_value=1,
                max_value=os.cpu_count(),
                value=optimal_cores,
                help=f"Dataset: {n_cells:,} cells × {n_genes:,} genes | "
                     f"Worker mem: ~{worker_mem_gb:.1f}GB | "
                     f"Available: {available_mem_gb}GB → Recommended: {optimal_cores} cores"
            )

            # Always use all preprocessed genes
            var_names = "all"

            init_weight_method = st.selectbox(
                "Weight initialization",
                ["correlation"],
                index=0,
                help="Method to initialize weights"
            )

            WX_method = st.selectbox(
                "Optimization method",
                ["lsq_linear", "LS", "LASSO", "Ridge", "constant", "LS_constrant"],
                index=0,
                help="Method for weight optimization"
            )

            WX_thres = st.number_input(
                "Weight threshold",
                min_value=1,
                max_value=100,
                value=20,
                help="Threshold for weights (default: 20)"
            )

            max_n_TF = st.number_input(
                "Max number of TFs",
                min_value=1,
                max_value=200,
                value=99,
                help="Maximum number of TFs per gene (default: 99)"
            )

            max_iter = st.number_input(
                "Max iterations",
                min_value=1,
                max_value=100,
                value=20,
                help="Maximum iterations in EM algorithm (default: 20)"
            )


            n_time_points = st.number_input(
                "Time points",
                min_value=100,
                max_value=2000,
                value=1000,
                help="Number of time points (default: 1000)"
            )


        st.subheader("Preprocessing Options")

        # Check data status for TFvelo
        col1, col2 = st.columns(2)
        with col1:
            if adata.X is not None:
                st.success("✅ Expression matrix found")
            else:
                st.warning("⚠️ Expression matrix not found")

        with col2:
            # Display gene count
            if adata.var is not None:
                st.success(f"✅ {adata.n_vars} genes available")


        # TFvelo preprocessing options
        st.subheader("TFvelo Preprocessing")

        if not TFVELO_AVAILABLE:
            st.error("TFvelo package not available. Please fix the installation.")
            st.stop()

        st.info(f"Selected databases: {', '.join(tf_databases)}")

        # Combined TFvelo workflow button
        if st.button("🚀 Run Complete TFvelo Analysis", type="primary"):
            if not tf_databases:
                st.error("Please select at least one TF database")
                st.stop()

            # Phase 1: Preprocessing
            with st.spinner("Phase 1/2: Running TFvelo preprocessing..."):
                try:
                    # Copy data for processing
                    adata_processed = adata.copy()

                    # Step 0: Apply selected data layer and type settings
                    st.info(f"Using data source: {selected_data_source}")
                    st.info(f"Data type: {data_type}")

                    # Set the X matrix and use_raw flag based on selected data source
                    # Note: For pancreas dataset, official demo uses use_raw=0 (default) even with raw counts
                    # This means: start from raw counts, normalize, then use normalized data for analysis
                    if selected_data_source == "adata.X (current)":
                        use_raw_flag = 0  # Using processed data (official demo default)
                        pass  # Already using adata.X
                    elif selected_data_source == "adata.raw.X (raw counts)":
                        use_raw_flag = 0  # Official demo default - use normalized data for analysis
                        adata_processed.X = adata_processed.raw.X.copy()
                    else:
                        use_raw_flag = 0  # Using layer data (official demo default)
                        # Extract layer name from the string (handle both old and new format)
                        if "adata.layers['" in selected_data_source:
                            if "] (" in selected_data_source:
                                # New format with shape info: "adata.layers['layer_name'] (shape)"
                                layer_name = selected_data_source.split("'")[1].split("'")[0]
                            else:
                                # Old format: "adata.layers['layer_name']"
                                layer_name = selected_data_source.split("'")[1]
                            adata_processed.X = adata_processed.layers[layer_name].copy()
                        else:
                            st.error(f"Cannot parse data source: {selected_data_source}")
                            # Fallback to current X
                            pass

                    # Step 1: Initial data setup and gene filtering (following demo order)
                    adata_processed.var_names_make_unique()
                    adata_processed.obs_names_make_unique()

                    # Store original gene names
                    adata_processed.uns['genes_all'] = np.array(adata_processed.var_names)

                    # Step 2: Data normalization based on type
                    if data_type == "Raw counts (apply full normalization)":
                        st.info("Applying TFvelo normalization (raw counts → normalized)...")

                        # Set up layers for TFvelo based on data layer selection
                        if is_rna_velocity_mode:
                            # RNA Velocity mode: Use spliced + unspliced → total (paper demo method)
                            st.info("📊 Using spliced + unspliced data (paper demo method)")
                            available_layers = list(adata_processed.layers.keys()) if adata_processed.layers else []
                            if "spliced" in available_layers and "unspliced" in available_layers:
                                spliced = adata_processed.layers["spliced"]
                                unspliced = adata_processed.layers["unspliced"]
                                if hasattr(spliced, 'todense'):
                                    spliced = spliced.todense()
                                if hasattr(unspliced, 'todense'):
                                    unspliced = unspliced.todense()
                                adata_processed.layers['total'] = spliced + unspliced
                            else:
                                st.error("RNA Velocity mode selected but spliced/unspliced layers not found!")
                                st.stop()
                        else:
                            # Standard mode: use selected data source
                            st.info("📊 Using selected data source")
                            adata_processed.layers['total'] = adata_processed.X.copy()
                        adata_processed.layers['total_raw'] = adata_processed.layers['total'].copy()

                        # Apply gene and cell filtering first (following demo order)
                        n_cells, n_genes = adata_processed.X.shape
                        st.info(f"Applying gene filtering (min_cells: {int(n_cells/50)})...")
                        sc.pp.filter_genes(adata_processed, min_cells=int(n_cells/50))

                        st.info(f"Applying cell filtering (min_genes: {int(n_genes/50)})...")
                        sc.pp.filter_cells(adata_processed, min_genes=int(n_genes/50))

                        # Apply TFvelo's normalization pipeline
                        st.info(f"Applying TFvelo normalization pipeline with {n_top_genes} variable genes...")
                        TFv.pp.filter_and_normalize(adata_processed, min_shared_counts=20, n_top_genes=n_top_genes, log=True)
                        adata_processed.X = adata_processed.layers['total'].copy()

                    else:
                        st.info("Using normalized data (skipping log transform)...")

                        # Set up layers for TFvelo based on data layer selection
                        if is_rna_velocity_mode:
                            # RNA Velocity mode: Use spliced + unspliced → total (paper demo method)
                            st.info("📊 Using spliced + unspliced data (paper demo method)")
                            available_layers = list(adata_processed.layers.keys()) if adata_processed.layers else []
                            if "spliced" in available_layers and "unspliced" in available_layers:
                                spliced = adata_processed.layers["spliced"]
                                unspliced = adata_processed.layers["unspliced"]
                                if hasattr(spliced, 'todense'):
                                    spliced = spliced.todense()
                                if hasattr(unspliced, 'todense'):
                                    unspliced = unspliced.todense()
                                adata_processed.layers['total'] = spliced + unspliced
                            else:
                                st.error("RNA Velocity mode selected but spliced/unspliced layers not found!")
                                st.stop()
                        else:
                            # Standard mode: use selected data source
                            st.info("📊 Using selected data source")
                            adata_processed.layers['total'] = adata_processed.X.copy()
                        adata_processed.layers['total_raw'] = adata_processed.layers['total'].copy()

                        # Apply gene and cell filtering first (following demo order)
                        n_cells, n_genes = adata_processed.X.shape
                        st.info(f"Applying gene filtering (min_cells: {int(n_cells/50)})...")
                        sc.pp.filter_genes(adata_processed, min_cells=int(n_cells/50))

                        st.info(f"Applying cell filtering (min_genes: {int(n_genes/50)})...")
                        sc.pp.filter_cells(adata_processed, min_genes=int(n_genes/50))

                        # Apply normalization without log transform (high variable gene selection only)
                        st.info(f"Applying gene selection with {n_top_genes} variable genes (no log transform)...")
                        TFv.pp.filter_and_normalize(adata_processed, min_shared_counts=20, n_top_genes=n_top_genes, log=False)
                        adata_processed.X = adata_processed.layers['total'].copy()

                    # Step 3: Convert gene names to uppercase (following demo order after normalization)
                    st.info("Converting gene names to uppercase...")
                    gene_names = []
                    for tmp in adata_processed.var_names:
                        gene_names.append(tmp.upper())
                    adata_processed.var_names = gene_names
                    adata_processed.var_names_make_unique()
                    adata_processed.obs_names_make_unique()

                    # Step 4: moments will automatically create M_total layer

                    # Step 4: Run TFvelo moments (equivalent to demo line 78, creates M_total automatically)
                    st.info("Computing moments...")
                    TFv.pp.moments(adata_processed, n_pcs=30, n_neighbors=n_neighbors)

                    # Step 5: Get TFs (equivalent to demo line 80)
                    st.info(f"Loading TF databases: {tf_databases}")
                    if force_include_tfs:
                        st.info(f"Important TFs: {force_include_tfs} (min_expr={min_expression_force}, min_cells={min_cells_force})")

                    try:
                        TFv.pp.get_TFs(
                            adata_processed,
                            databases=tf_databases,
                            min_expression=min_expression,
                            min_cells=min_cells,
                            force_include_tfs=force_include_tfs,
                            min_expression_force=min_expression_force,
                            min_cells_force=min_cells_force
                        )
                        st.success("✅ TF databases loaded successfully")
                    finally:
                        pass

                    # Step 6: Store processed data
                    adata_processed.uns['genes_pp'] = np.array(adata_processed.var_names)
                    st.session_state.adata_processed = adata_processed

                    st.success("✅ TFvelo preprocessing completed")
                    st.info(f"Data shape: {adata_processed.n_obs} cells × {adata_processed.n_vars} genes")

                    # Display TF information
                    if 'n_TFs' in adata_processed.var.columns:
                        total_tf_relations = adata_processed.var['n_TFs'].sum()
                        st.info(f"Total TF-target relationships: {total_tf_relations}")

                except Exception as e:
                    st.error(f"Error during preprocessing: {str(e)}")
                    st.text(traceback.format_exc())
                    st.stop()

            # Phase 2: Run dynamics recovery immediately after preprocessing
            with st.spinner("Phase 2/2: Running TFvelo dynamics recovery (this may take several minutes)..."):
                try:
                    # TFvelo main analysis following demo exactly
                    try:
                        st.info("Starting dynamics recovery...")

                        # Automatically determine use_raw based on data type selection
                        if data_type == "Raw counts (apply full normalization)":
                            use_raw_param = 0  # Use normalized data after preprocessing
                            st.info("🔢 Using normalized data (use_raw=0) - raw counts will be normalized during preprocessing")
                        else:
                            use_raw_param = 0  # Use normalized data
                            st.info("🔢 Using normalized data (use_raw=0) - data is already normalized")

                        # Simple progress display
                        progress_container = st.container()
                        progress_bar = progress_container.progress(0)
                        status_text = progress_container.empty()
                        gene_info = progress_container.empty()

                        # Initialize display
                        n_genes = adata_processed.n_vars
                        status_text.text("🔄 Starting TFvelo dynamics recovery...")
                        gene_info.text(f"Will process {n_genes} genes")
                        progress_bar.progress(0.1)

                        # Show detailed processing information
                        status_text.text("Running TFvelo dynamics recovery (this may take several minutes)...")
                        gene_info.text(f"Processing {n_genes} genes using {n_jobs} cores...")
                        progress_bar.progress(0.3)

                        # Provide timing estimate based on genes and cells
                        n_cells = adata_processed.n_obs
                        # Empirical estimates:
                        # - Demo pancreas data: ~2 hours
                        # - Reference (2000 genes × 300 cells): ~20 minutes
                        # Time scales mainly with genes and complexity of data
                        # Note: n_top_genes from UI determines gene count after preprocessing

                        if n_genes > 10000 and n_cells > 1000:
                            # Large dataset like pancreas demo
                            estimated_hours = max(1, (n_genes / 10000) * 2)
                            if estimated_hours >= 1:
                                st.info(f"⏱️ Estimated time: ~{estimated_hours:.1f} hours for {n_genes} genes × {n_cells} cells")
                            else:
                                estimated_minutes = int(estimated_hours * 60)
                                st.info(f"⏱️ Estimated time: ~{estimated_minutes} minutes for {n_genes} genes × {n_cells} cells")
                        else:
                            # Smaller dataset - base on 2000 genes × 300 cells = 20 minutes
                            base_time = (n_genes / 2000) * 20  # Scale from reference: 2000 genes = 20 min
                            cell_factor = max(1, n_cells / 300)  # Scale from reference: 300 cells
                            estimated_minutes = max(1, int(base_time * cell_factor))
                            st.info(f"⏱️ Estimated time: ~{estimated_minutes} minutes for {n_genes} genes × {n_cells} cells")

                        # Show processing details
                        with st.expander("Processing Details", expanded=False):
                            st.write(f"- **Genes to process**: {n_genes}")
                            st.write(f"- **Cells to process**: {n_cells}")
                            st.write(f"- **Total complexity**: {n_genes * n_cells:,} gene-cell pairs")
                            st.write(f"- **CPU cores**: {n_jobs}")
                            st.write(f"- **Max iterations per gene**: {max_iter}")
                            st.write(f"- **Using data**: {'Raw counts' if use_raw_flag else 'Normalized'}")
                            st.write("- **Note**: Elapsed time tracking available during processing")

                        # Real progress tracking by monitoring TFvelo's temporary files/state
                        import time
                        import threading
                        import os

                        # Simple status display without false progress indication
                        def monitor_processing():
                            """Monitor processing without false progress simulation"""
                            start_time = time.time()

                            while progress_state['running']:
                                elapsed = time.time() - start_time
                                elapsed_min = int(elapsed // 60)
                                elapsed_sec = int(elapsed % 60)

                                status_text.text(f"🔄 TFvelo processing... ({elapsed_min:02d}:{elapsed_sec:02d})")
                                gene_info.text(f"Processing {n_genes} genes with {n_jobs} cores")

                                time.sleep(1)

                        # Initialize progress state
                        progress_state = {'running': True}

                        # Start monitoring thread
                        monitor_thread = threading.Thread(target=monitor_processing, daemon=True)
                        monitor_thread.start()

                        # Suppress TFvelo debug output
                        import logging
                        import sys
                        from io import StringIO

                        # Capture and suppress stdout/stderr during TFvelo execution
                        old_stdout = sys.stdout
                        old_stderr = sys.stderr
                        sys.stdout = StringIO()
                        sys.stderr = StringIO()

                        # Set logging to WARNING to suppress DEBUG messages
                        logging.getLogger().setLevel(logging.WARNING)

                        try:
                            # Run TFvelo
                            flag = TFv.tl.recover_dynamics(
                            adata_processed,
                            n_jobs=n_jobs,
                            max_iter=max_iter,
                            var_names=var_names,
                            WX_method=WX_method,
                            WX_thres=WX_thres,
                            max_n_TF=max_n_TF,
                            n_top_genes=n_top_genes,
                            fit_scaling=True,
                            use_raw=use_raw_flag,
                            init_weight_method=init_weight_method,
                            n_time_points=n_time_points
                        )
                        finally:
                            # Restore stdout/stderr
                            sys.stdout = old_stdout
                            sys.stderr = old_stderr

                        # Stop progress simulation
                        progress_state['running'] = False

                        # Complete progress
                        progress_bar.progress(1.0)
                        status_text.text("✅ TFvelo dynamics recovery completed!")
                        gene_info.text(f"Successfully processed all {n_genes} genes")

                        if flag == False:
                            st.error("TFvelo dynamics recovery failed")
                            st.stop()

                        # Handle highly variable genes metadata
                        if 'highly_variable_genes' in adata_processed.var.keys():
                            # Convert boolean to string for saving
                            if adata_processed.var['highly_variable_genes'][0] in [True, False]:
                                adata_processed.var['highly_variable_genes'] = adata_processed.var['highly_variable_genes'].map({True: 'True', False: 'False'})

                        st.success("✅ TFvelo dynamics recovery completed successfully")

                        # TFvelo post-processing (following demo exactly)
                        st.info("🔄 Running TFvelo post-processing...")

                        # 1. Create velocity layer (from TFvelo_analysis_demo.py line 136)
                        if 'velo_hat' in adata_processed.layers and 'fit_scaling_y' in adata_processed.var:
                            n_cells = adata_processed.shape[0]
                            expanded_scaling_y = np.expand_dims(np.array(adata_processed.var['fit_scaling_y']), 0).repeat(n_cells, axis=0)
                            adata_processed.layers['velocity'] = adata_processed.layers['velo_hat'] / expanded_scaling_y
                            st.success("✅ Velocity layer created")

                        # 2. Calculate loss values
                        if 'loss' in adata_processed.varm:
                            losses = adata_processed.varm['loss'].copy()
                            losses[np.isnan(losses)] = 1e6
                            adata_processed.var['min_loss'] = losses.min(1)
                            st.success("✅ Loss metrics calculated")

                        # 3. Generate PCA and UMAP if not present
                        if 'X_pca' not in adata_processed.obsm.keys():
                            sc.tl.pca(adata_processed, n_comps=50, svd_solver='arpack')
                            st.success("✅ PCA computed")

                        if 'X_umap' not in adata_processed.obsm.keys():
                            sc.tl.umap(adata_processed)
                            st.success("✅ UMAP computed")

                        # 4. Calculate velocity graph and pseudotime
                        if 'velocity' in adata_processed.layers:
                            TFv.tl.velocity_graph(adata_processed, basis=None, vkey='velocity', xkey='M_total')
                            TFv.tl.velocity_pseudotime(adata_processed, vkey='velocity', modality='M_total')
                            st.success("✅ Velocity graph and pseudotime computed")

                        # 5. Apply TFvelo demo functions (following exact demo methodology)
                        if 'fit_t_raw' in adata_processed.layers:
                            st.info("🔄 Applying TFvelo demo functions...")

                            # Apply get_sort_t function from demo
                            adata_processed = get_sort_t(adata_processed)
                            st.success("✅ Latent time analysis completed (get_sort_t)")

                            # Calculate Spearman correlation with pseudotime
                            if 'velocity_pseudotime' in adata_processed.obs:
                                adata_processed = get_metric_pseudotime(adata_processed)
                                st.success("✅ Pseudotime correlation metrics calculated")

                        # Memory optimization: clear large varm matrices no longer needed
                        # These are only used during get_TFs() and recover_dynamics()
                        large_varm_keys = ['TFs', 'TFs_id', 'TFs_times', 'TFs_correlation', 'knockTF_Log2FC']
                        cleared_size = 0
                        for key in large_varm_keys:
                            if key in adata_processed.varm:
                                cleared_size += adata_processed.varm[key].nbytes if hasattr(adata_processed.varm[key], 'nbytes') else 0
                                del adata_processed.varm[key]
                        if cleared_size > 0:
                            st.info(f"🧹 Cleared {cleared_size / 1e9:.1f}GB of intermediate TF matrices from memory")
                        import gc
                        gc.collect()

                        # Store unfiltered data in session_state for later filtering adjustment
                        st.session_state.adata_tfvelo_unfiltered = adata_processed.copy()
                        st.success("✅ TFvelo processing completed (before gene filtering)")
                        st.info("📊 Unfiltered data saved to session_state. You can adjust filtering parameters later.")

                        # Store results
                        st.session_state.adata_tfvelo = adata_processed
                        st.session_state.tfvelo_results = {
                            'completed': True,
                            'tf_databases': tf_databases,
                            'parameters': {
                                'n_jobs': n_jobs,
                                'max_iter': max_iter,
                                'var_names': var_names,
                                'WX_method': WX_method,
                                'WX_thres': WX_thres,
                                'max_n_TF': max_n_TF,
                                'n_top_genes': n_top_genes,
                                'init_weight_method': init_weight_method,
                                'n_time_points': n_time_points,
                                'use_raw': 0
                            },
                            'tf_filtering': {
                                'min_expression': min_expression,
                                'min_cells': min_cells,
                                'force_include_tfs': force_include_tfs,
                                'min_expression_force': min_expression_force,
                                'min_cells_force': min_cells_force
                            }
                        }

                        # Display results summary
                        st.subheader("Analysis Results")
                        if 'velocity' in adata_processed.layers:
                            st.success("✅ Velocity computed")
                        if 'latent_time' in adata_processed.obs:
                            st.success("✅ Latent time computed")

                        # Show some statistics
                        if 'n_TFs' in adata_processed.var.columns:
                            total_tf_relations = adata_processed.var['n_TFs'].sum()
                            st.info(f"Total TF-target relationships used: {total_tf_relations}")

                        st.info("💡 Go to the 'Visualization' tab to explore results")


                    finally:
                        pass

                except Exception as e:
                    st.error(f"Error during TFvelo analysis: {str(e)}")
                    st.text(traceback.format_exc())


    # Tab 3: Gene Filtering
    with tab3:
        st.header("🔬 Gene Filtering")

        # Check if unfiltered data is available
        if 'adata_tfvelo_unfiltered' not in st.session_state:
            st.warning("⚠️ No unfiltered data available. Please run the analysis in the 'TFvelo Setup & Analysis' tab first.")
            st.info("If you loaded a saved analysis file, Gene Filtering is not available. Please use the Visualization tab.")
            # Don't use st.stop() - it blocks other tabs from rendering
        else:
            adata_unfiltered = st.session_state.adata_tfvelo_unfiltered

            st.info(f"📊 **Unfiltered data**: {adata_unfiltered.n_vars} genes, {adata_unfiltered.n_obs} cells")

            st.markdown("""
            **📝 About default values**: All filtering parameters below use the **official TFvelo default values**.
            - `loss_percent_thres` = 50
            - `cell_percent` = 10%
            - `non_blank_gene` filter = ON
            - `spearmanr_thres` = 0.8

            If the number of genes remaining is too low, consider relaxing `spearmanr_thres` in particular.
            """)

            # Filtering parameters
            st.subheader("Filtering Parameters")

            col1, col2 = st.columns(2)

            with col1:
                # Filter 1: Loss percentile
                if 'min_loss' in adata_unfiltered.var:
                    loss_percentile = st.slider(
                        "Loss percentile threshold (loss_percent_thres) (%)",
                        0, 100, 50,
                        help="""TFvelo default=50.

**Meaning**: Filters genes based on TFvelo model goodness-of-fit (loss function value).
Genes with smaller loss have better model fit for gene dynamics.
Retains genes with loss values below this percentile.

Example: At 50%, genes with loss in the bottom 50% (= top 50% model fit) are retained."""
                    )
                    thres_loss = np.percentile(adata_unfiltered.var['min_loss'], loss_percentile)
                    genes_pass_loss = (adata_unfiltered.var['min_loss'] < thres_loss).sum()
                    st.metric("Loss filter", f"{genes_pass_loss} / {adata_unfiltered.n_vars} genes")

                # Filter 2: Cell count threshold
                if 'n_cells' in adata_unfiltered.var:
                    cell_percent = st.slider(
                        "Minimum cell percentage (%) (cell_percent)",
                        0, 100, 10,
                        help="""TFvelo default=10%.

**Meaning**: Filters genes based on the proportion of cells expressing the gene.
Genes expressed in too few cells are excluded due to low statistical reliability.

Example: At 10%, only genes expressed in at least 10% of all cells are retained."""
                    )
                    thres_n_cells = adata_unfiltered.n_obs * (cell_percent / 100.0)
                    genes_pass_cells = (adata_unfiltered.var['n_cells'] > thres_n_cells).sum()
                    st.metric("Cell count filter", f"{genes_pass_cells} / {adata_unfiltered.n_vars} genes")

            with col2:
                # Filter 3: Non-blank genes
                apply_blank_filter = st.checkbox(
                    "Non-blank gene filter (non_blank_gene)",
                    value=True,
                    help="""TFvelo default=ON.

**Meaning**: A filter that excludes genes with large gaps (blanks) in expression.
Retains only genes with continuously changing expression along pseudotime.
Genes with unnatural gaps in expression patterns have lower velocity estimation reliability.

non_blank_gene == 0: High-quality genes with continuous expression
non_blank_gene != 0: Genes with expression gaps (excluded)"""
                )
                if apply_blank_filter and 'non_blank_gene' in adata_unfiltered.var:
                    genes_pass_blank = (adata_unfiltered.var['non_blank_gene'] == 0).sum()
                    st.metric("Non-blank filter", f"{genes_pass_blank} / {adata_unfiltered.n_vars} genes")

                # Filter 4: Pseudotime correlation (Spearman)
                if 'spearmanr_pseudotime' in adata_unfiltered.var:
                    corr_threshold = st.slider(
                        "Spearman correlation threshold (spearmanr_thres)",
                        0.0, 1.0, 0.8,
                        step=0.05,
                        help="""TFvelo default=0.8.

**Meaning**: Filters based on Spearman correlation between gene-specific time and pseudotime.
Genes with higher correlation change expression more regularly along the cell differentiation trajectory.
Only genes with correlation above this threshold are retained.

Example: At 0.8, only genes strongly correlated with pseudotime are retained (strict condition).
If the number of genes is too low, consider relaxing to around 0.5-0.6."""
                    )
                    genes_pass_corr = (adata_unfiltered.var['spearmanr_pseudotime'] > corr_threshold).sum()
                    st.metric("Correlation filter", f"{genes_pass_corr} / {adata_unfiltered.n_vars} genes")

            # Apply filters button
            st.markdown("---")
            if st.button("🔄 Apply Filters", type="primary"):
                with st.spinner("Applying filters..."):
                    adata_filtered = adata_unfiltered.copy()

                    # Apply each filter
                    filter_log = []

                    # Filter 1: Loss
                    if 'min_loss' in adata_filtered.var:
                        n_before = adata_filtered.n_vars
                        adata_filtered = adata_filtered[:, adata_filtered.var['min_loss'] < thres_loss]
                        filter_log.append(f"Loss filter ({loss_percentile}%): {n_before} → {adata_filtered.n_vars} genes")

                    # Filter 2: Cell count
                    if 'n_cells' in adata_filtered.var:
                        n_before = adata_filtered.n_vars
                        adata_filtered = adata_filtered[:, adata_filtered.var['n_cells'] > thres_n_cells]
                        filter_log.append(f"Cell count filter (>{cell_percent}%): {n_before} → {adata_filtered.n_vars} genes")

                    # Filter 3: Non-blank
                    if apply_blank_filter and 'non_blank_gene' in adata_filtered.var:
                        n_before = adata_filtered.n_vars
                        adata_filtered = adata_filtered[:, adata_filtered.var['non_blank_gene'] == 0]
                        filter_log.append(f"Non-blank filter: {n_before} → {adata_filtered.n_vars} genes")

                    # Filter 4: Pseudotime correlation
                    if 'spearmanr_pseudotime' in adata_filtered.var:
                        n_before = adata_filtered.n_vars
                        adata_filtered = adata_filtered[:, adata_filtered.var['spearmanr_pseudotime'] > corr_threshold]
                        filter_log.append(f"Correlation filter (>{corr_threshold}): {n_before} → {adata_filtered.n_vars} genes")

                    # Display filtering results
                    st.subheader("Filtering Results")
                    for log in filter_log:
                        st.info(log)

                    # Check if genes remain
                    if adata_filtered.n_vars == 0:
                        st.error("❌ **No genes remaining after filtering!**")
                        st.error("Please relax the filtering parameters.")
                    else:
                        st.success(f"✅ **Final result**: {adata_filtered.n_vars} genes, {adata_filtered.n_obs} cells")

                        # Re-compute velocity graph
                        with st.spinner("Re-computing velocity graph..."):
                            try:
                                import TFvelo as TFv
                                TFv.tl.velocity_graph(adata_filtered, basis=None, vkey='velocity', xkey='M_total')
                                st.success("✅ Velocity graph computed successfully")

                                # Store filtered data
                                st.session_state.adata_tfvelo = adata_filtered
                                st.session_state.tfvelo_results['completed'] = True

                                st.balloons()
                                st.success("🎉 Filtering complete! Check the results in the Visualization tab.")

                            except Exception as e:
                                st.error(f"❌ Velocity graph computation failed: {str(e)}")
                                st.text(traceback.format_exc())

    # Tab 4: Visualization (renamed from tab3)
    with tab4:
        st.header("TFvelo Visualization")

        # Check if TFvelo analysis is available (from current session or loaded data)
        if 'adata_tfvelo' not in st.session_state or 'completed' not in st.session_state.tfvelo_results:
            st.warning("🔍 No TFvelo analysis available for visualization.")
            st.info("Please either:")
            st.markdown("1. **Load saved analysis**: Go to Data Input → 'Load saved TFvelo analysis for visualization'")
            st.markdown("2. **Run new analysis**: Upload data → Complete analysis in 'TFvelo Setup & Analysis' tab")
            st.markdown("3. **Apply gene filtering**: Go to 'Gene Filtering' tab")
            st.stop()

        adata_tfvelo = st.session_state.adata_tfvelo


        if 'completed' in st.session_state.tfvelo_results:
            # Select embedding for visualization
            available_embeddings = [key for key in adata_tfvelo.obsm.keys() if key.startswith('X_')]

            if available_embeddings:
                col1, col2 = st.columns(2)

                with col1:
                    selected_embedding = st.selectbox(
                        "Select embedding for projection",
                        available_embeddings,
                        format_func=lambda x: x.replace('X_', '').upper(),
                        help="Choose which dimensionality reduction to use for visualization"
                    )

                with col2:
                    # Select coloring variable (TFvelo specific)
                    color_options = []

                    # Add TFvelo-specific metrics
                    if 'velocity_pseudotime' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_pseudotime')
                    if 'velocity_length' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_length')
                    if 'velocity_confidence' in adata_tfvelo.obs.columns:
                        color_options.append('velocity_confidence')

                    # Add cluster information
                    if 'clusters' in adata_tfvelo.obs.columns:
                        color_options.append('clusters')

                    # Add other metadata
                    for col in adata_tfvelo.obs.columns:
                        if col not in color_options:
                            color_options.append(col)

                    color_by = st.selectbox(
                        "Color by",
                        color_options,
                        index=0 if 'velocity_pseudotime' in color_options else 0,
                        help="Choose variable for coloring cells"
                    )

                # Visualization options
                st.subheader("Visualization Options")

                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    arrow_length = st.slider("Arrow length", 1, 10, 4)
                with col2:
                    arrow_size = st.slider("Arrow size", 0.5, 3.0, 1.0)
                with col3:
                    density = st.slider("Stream density", 0.5, 4.0, 2.0)
                with col4:
                    point_size = st.slider("Cell size", 10, 1000, 100, help="Size of individual cells in scatter plots (points²)")

                # Advanced visualization options
                with st.expander("Advanced Visualization Options"):
                    col1, col2, col3, col4 = st.columns(4)

                    with col1:
                        alpha_value = st.slider("Cell transparency", 0.1, 1.0, 0.7, 0.1, help="Transparency of cells (lower = more transparent)")

                    with col2:
                        plot_height = st.slider("Plot height", 4, 12, 6, help="Height of the plot in inches")

                    with col3:
                        # Determine if selected color is categorical or continuous
                        if color_by in adata_tfvelo.obs.columns:
                            is_categorical = adata_tfvelo.obs[color_by].dtype == 'category' or adata_tfvelo.obs[color_by].dtype == 'object'

                            if is_categorical:
                                # Discrete/categorical colormap options
                                colormap_options = ['tab10', 'tab20', 'Set1', 'Set2', 'Set3', 'Pastel1', 'Pastel2', 'Dark2', 'Accent']
                                colormap_help = "Colormap for categorical data"
                            else:
                                # Continuous colormap options with reverse variants
                                base_colormaps = ['viridis', 'plasma', 'inferno', 'magma', 'cividis', 'Blues', 'Reds', 'Greens',
                                                'Oranges', 'Purples', 'BuGn', 'BuPu', 'GnBu', 'OrRd', 'PuBu', 'PuRd',
                                                'YlGn', 'YlOrRd', 'RdYlBu', 'RdYlGn', 'Spectral', 'coolwarm', 'bwr',
                                                'seismic', 'RdBu', 'PiYG', 'PRGn', 'BrBG', 'RdGy', 'PuOr', 'gnuplot']
                                colormap_options = []
                                for cmap in base_colormaps:
                                    colormap_options.append(cmap)
                                    colormap_options.append(cmap + '_r')  # Add reversed version
                                colormap_help = "Colormap for continuous data (_r suffix for reversed)"
                        else:
                            base_colormaps = ['viridis', 'plasma', 'inferno', 'magma']
                            colormap_options = []
                            for cmap in base_colormaps:
                                colormap_options.append(cmap)
                                colormap_options.append(cmap + '_r')
                            colormap_help = "Default colormap (_r suffix for reversed)"

                        selected_colormap = st.selectbox(
                            "Colormap",
                            colormap_options,
                            index=0,
                            help=colormap_help
                        )

                    with col4:
                        # Removed show_colorbar option - always show colorbar
                        st.write("Colorbar: Always shown")

                # Additional Analyses Section (removed velocity confidence/length as not in TFvelo demo)
                # These options are commented out as TFvelo demo doesn't compute velocity confidence/length
                # If needed, can be re-enabled with proper velocity layer computation

                # Generate plots
                if st.button("Generate Plots"):
                    # Check available data
                    with st.expander("📊 Available Data Layers", expanded=False):
                        col1, col2 = st.columns(2)

                        with col1:
                            st.write("**Layers:**")
                            for layer in adata_tfvelo.layers.keys():
                                st.write(f"- {layer}")

                        with col2:
                            st.write("**Observations:**")
                            velocity_cols = [col for col in adata_tfvelo.obs.columns if 'velocity' in col.lower()]
                            for col in velocity_cols[:10]:  # Show first 10
                                st.write(f"- {col}")

                    # Try TFvelo-specific plotting
                    try:
                        import matplotlib.pyplot as plt
                        fig, axes = plt.subplots(1, 2, figsize=(16, plot_height))

                        # Set tight layout parameters to prevent clipping
                        plt.subplots_adjust(left=0.05, right=0.88, bottom=0.1, top=0.9, wspace=0.25)

                        basis = selected_embedding.replace('X_', '')

                        # Plot 1: TFvelo velocity stream with colormap fix only
                        if 'velocity' in adata_tfvelo.layers:
                            try:
                                # Apply colormap for both categorical and continuous data
                                cmap_param = selected_colormap if color_by in adata_tfvelo.obs.columns else None

                                scv.pl.velocity_embedding_stream(
                                    adata_tfvelo,
                                    basis=basis,
                                    ax=axes[0],
                                    density=density,
                                    arrow_size=arrow_size,
                                    show=False,
                                    color=color_by if color_by in adata_tfvelo.obs.columns else None,
                                    cmap=cmap_param,
                                )

                                # Fix colors for categorical data to match right plot
                                for collection in axes[0].collections:
                                    if hasattr(collection, 'get_offsets') and len(collection.get_offsets()) > 0:
                                        # Fix colors for categorical data using selected colormap
                                        if color_by in adata_tfvelo.obs.columns:
                                            color_data = adata_tfvelo.obs[color_by]
                                            if hasattr(color_data, 'dtype') and (color_data.dtype == 'category' or color_data.dtype == 'object'):
                                                # Categorical data - apply colormap manually (use same order as right plot)
                                                unique_vals = color_data.unique()
                                                try:
                                                    colors = getattr(plt.cm, selected_colormap)(np.linspace(0, 1, len(unique_vals)))
                                                except AttributeError:
                                                    colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))

                                                # Create color mapping for each cell
                                                color_dict = dict(zip(unique_vals, colors))
                                                point_colors = [color_dict[val] for val in color_data]

                                                # Apply colors to collection
                                                if hasattr(collection, 'set_facecolors'):
                                                    collection.set_facecolors(point_colors)
                                                elif hasattr(collection, 'set_color'):
                                                    collection.set_color(point_colors)

                                axes[0].set_title(f"Velocity Stream ({basis.upper()})")

                            except Exception as e:
                                # Final fallback: Simple embedding with consistent manual approach
                                embedding = adata_tfvelo.obsm[selected_embedding]

                                # Use same coloring approach as main implementation
                                if color_by in adata_tfvelo.obs.columns:
                                    color_data = adata_tfvelo.obs[color_by]
                                    if color_data.dtype == 'category' or color_data.dtype == 'object':
                                        # Categorical coloring - same method as main
                                        unique_vals = color_data.unique()
                                        try:
                                            colors = getattr(plt.cm, selected_colormap)(np.linspace(0, 1, len(unique_vals)))
                                        except AttributeError:
                                            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))

                                        for i, val in enumerate(unique_vals):
                                            mask = color_data == val
                                            axes[0].scatter(embedding[mask, 0], embedding[mask, 1],
                                                          c=[colors[i]], s=point_size, alpha=alpha_value, label=val)
                                    else:
                                        # Continuous coloring
                                        axes[0].scatter(embedding[:, 0], embedding[:, 1],
                                                      c=color_data, cmap=selected_colormap, s=point_size, alpha=alpha_value)
                                else:
                                    # Default coloring
                                    axes[0].scatter(embedding[:, 0], embedding[:, 1],
                                                  c='lightgray', s=point_size, alpha=alpha_value)

                                axes[0].set_title(f"Embedding ({basis.upper()})")
                        else:
                            # No velocity layer - show embedding only (use same manual approach)
                            embedding = adata_tfvelo.obsm[selected_embedding]

                            # Use same coloring approach as main implementation
                            if color_by in adata_tfvelo.obs.columns:
                                color_data = adata_tfvelo.obs[color_by]
                                if color_data.dtype == 'category' or color_data.dtype == 'object':
                                    # Categorical coloring - same method as main
                                    unique_vals = color_data.unique()
                                    try:
                                        colors = getattr(plt.cm, selected_colormap)(np.linspace(0, 1, len(unique_vals)))
                                    except AttributeError:
                                        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))

                                    for i, val in enumerate(unique_vals):
                                        mask = color_data == val
                                        axes[0].scatter(embedding[mask, 0], embedding[mask, 1],
                                                      c=[colors[i]], s=point_size, alpha=alpha_value, label=val)
                                else:
                                    # Continuous coloring
                                    axes[0].scatter(embedding[:, 0], embedding[:, 1],
                                                  c=color_data, cmap=selected_colormap, s=point_size, alpha=alpha_value)
                            else:
                                # Default coloring
                                axes[0].scatter(embedding[:, 0], embedding[:, 1],
                                              c='lightgray', s=point_size, alpha=alpha_value)

                            axes[0].set_title(f"Embedding ({basis.upper()})")

                        # Plot 2: Pseudotime or selected coloring
                        try:
                            if color_by == 'velocity_pseudotime' and 'velocity_pseudotime' in adata_tfvelo.obs.columns:
                                # Use TFvelo's scatter plot for pseudotime
                                try:
                                    # Try to disable automatic colorbar in TFvelo
                                    try:
                                        TFv.pl.scatter(
                                            adata_tfvelo,
                                            basis=basis,
                                            color='velocity_pseudotime',
                                            cmap=selected_colormap,
                                            ax=axes[1],
                                            colorbar=False,  # Try to disable automatic colorbar
                                            show=False
                                        )
                                    except TypeError:
                                        # If colorbar parameter not supported, use without it
                                        TFv.pl.scatter(
                                            adata_tfvelo,
                                            basis=basis,
                                            color='velocity_pseudotime',
                                            cmap=selected_colormap,
                                            ax=axes[1],
                                            show=False
                                        )

                                    # Manually adjust scatter points for size and alpha
                                    for collection in axes[1].collections:
                                        if hasattr(collection, 'set_sizes'):
                                            collection.set_sizes([point_size] * len(collection.get_offsets()))
                                        if hasattr(collection, 'set_alpha'):
                                            collection.set_alpha(alpha_value)

                                    # Always add colorbar
                                    # Find the scatter collection for colorbar
                                    scatter_collection = None
                                    for collection in axes[1].collections:
                                        if hasattr(collection, 'get_array') and collection.get_array() is not None:
                                            scatter_collection = collection
                                            break
                                    if scatter_collection is not None:
                                        try:
                                            cbar = plt.colorbar(scatter_collection, ax=axes[1])
                                            cbar.set_label('velocity_pseudotime')
                                        except Exception:
                                            pass
                                except Exception:
                                    # Fallback to manual plotting
                                    embedding = adata_tfvelo.obsm[selected_embedding]
                                    scatter = axes[1].scatter(embedding[:, 0], embedding[:, 1],
                                                            c=adata_tfvelo.obs['velocity_pseudotime'],
                                                            cmap=selected_colormap, s=point_size, alpha=alpha_value)
                                    try:
                                        cbar = plt.colorbar(scatter, ax=axes[1])
                                        cbar.set_label('velocity_pseudotime')
                                    except Exception:
                                        pass
                                axes[1].set_title(f"Velocity Pseudotime ({basis.upper()})")
                                # Remove plot frame/spines and ticks
                                axes[1].spines['top'].set_visible(False)
                                axes[1].spines['right'].set_visible(False)
                                axes[1].spines['bottom'].set_visible(False)
                                axes[1].spines['left'].set_visible(False)
                                axes[1].set_xticks([])
                                axes[1].set_yticks([])
                            else:
                                # Generic scatter plot
                                embedding = adata_tfvelo.obsm[selected_embedding]
                                if color_by in adata_tfvelo.obs.columns:
                                    if adata_tfvelo.obs[color_by].dtype == 'category' or adata_tfvelo.obs[color_by].dtype == 'object':
                                        # Categorical coloring
                                        unique_vals = adata_tfvelo.obs[color_by].unique()
                                        try:
                                            colors = getattr(plt.cm, selected_colormap)(np.linspace(0, 1, len(unique_vals)))
                                        except AttributeError:
                                            # Fallback to default colormap if selected one doesn't exist
                                            colors = plt.cm.tab10(np.linspace(0, 1, len(unique_vals)))
                                        for i, val in enumerate(unique_vals):
                                            mask = adata_tfvelo.obs[color_by] == val
                                            axes[1].scatter(embedding[mask, 0], embedding[mask, 1],
                                                          c=[colors[i]], label=val, s=point_size, alpha=alpha_value)
                                        # Remove legend box frame
                                        legend = axes[1].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
                                        legend.set_frame_on(False)
                                    else:
                                        # Continuous coloring
                                        scatter = axes[1].scatter(embedding[:, 0], embedding[:, 1],
                                                                c=adata_tfvelo.obs[color_by],
                                                                cmap=selected_colormap, s=point_size, alpha=alpha_value)
                                        try:
                                            cbar = plt.colorbar(scatter, ax=axes[1])
                                            cbar.set_label(color_by)
                                        except Exception as cbar_error:
                                            # Skip colorbar if it fails
                                            pass
                                    axes[1].set_title(f"{color_by} ({basis.upper()})")
                                    # Remove plot frame/spines and ticks
                                    axes[1].spines['top'].set_visible(False)
                                    axes[1].spines['right'].set_visible(False)
                                    axes[1].spines['bottom'].set_visible(False)
                                    axes[1].spines['left'].set_visible(False)
                                    axes[1].set_xticks([])
                                    axes[1].set_yticks([])
                                else:
                                    # Default embedding
                                    axes[1].scatter(embedding[:, 0], embedding[:, 1],
                                                  c='orange', s=point_size, alpha=alpha_value)
                                    axes[1].set_title(f"Embedding ({basis.upper()})")
                                    # Remove plot frame/spines and ticks
                                    axes[1].spines['top'].set_visible(False)
                                    axes[1].spines['right'].set_visible(False)
                                    axes[1].spines['bottom'].set_visible(False)
                                    axes[1].spines['left'].set_visible(False)
                                    axes[1].set_xticks([])
                                    axes[1].set_yticks([])
                        except Exception as e:
                            axes[1].text(0.5, 0.5, f"Plot error: {str(e)[:30]}",
                                       transform=axes[1].transAxes, ha='center', va='center')

                        # Adjust aspect ratio and margins
                        for ax in axes:
                            # Use auto aspect ratio instead of equal to reduce vertical stretching
                            ax.set_aspect('auto')
                            # Add some margin to prevent clipping
                            xlim = ax.get_xlim()
                            ylim = ax.get_ylim()
                            x_margin = (xlim[1] - xlim[0]) * 0.03
                            y_margin = (ylim[1] - ylim[0]) * 0.03
                            ax.set_xlim(xlim[0] - x_margin, xlim[1] + x_margin)
                            ax.set_ylim(ylim[0] - y_margin, ylim[1] + y_margin)

                        plt.tight_layout()
                        st.pyplot(fig, bbox_inches='tight', dpi=150)

                        # Add download buttons for the plot
                        import io
                        col1, col2 = st.columns(2)

                        with col1:
                            # PNG download
                            buffer_png = io.BytesIO()
                            fig.savefig(buffer_png, format='png', bbox_inches='tight', dpi=300, facecolor='white')
                            buffer_png.seek(0)
                            st.download_button(
                                label="📥 Download Plot (PNG)",
                                data=buffer_png.getvalue(),
                                file_name=f"tfvelo_plot_{datetime.now().strftime('%Y%m%d_%H%M%S')}.png",
                                mime="image/png"
                            )


                    except Exception as e:
                        st.error(f"Plotting error: {str(e)}")
                        st.text(f"Full error: {traceback.format_exc()}")

                    # Additional velocity analyses removed - not computed in TFvelo demo
                    # TFvelo focuses on transcription factor dynamics rather than velocity confidence/length metrics

                # TF Activity Visualization
                if 'filtered_tf_targets' in st.session_state.tfvelo_results:
                    st.subheader("TF Activity Analysis")

                    tf_targets_used = st.session_state.tfvelo_results['filtered_tf_targets']

                    col1, col2 = st.columns(2)

                    with col1:
                        if st.checkbox("Show TF Activity Heatmap"):
                            tf_activity_cols = [col for col in adata.obs.columns if col.startswith('TF_activity_')]

                            if len(tf_activity_cols) > 0:
                                fig, ax = plt.subplots(figsize=(8, 6))

                                # Create heatmap of TF activities
                                tf_activity_data = adata.obs[tf_activity_cols]
                                tf_names = [col.replace('TF_activity_', '') for col in tf_activity_cols]

                                # Sort cells by first TF activity for visualization
                                sort_idx = np.argsort(tf_activity_data.iloc[:, 0])
                                sorted_data = tf_activity_data.iloc[sort_idx]

                                im = ax.imshow(sorted_data.T, aspect='auto', cmap='viridis')
                                ax.set_yticks(range(len(tf_names)))
                                ax.set_yticklabels(tf_names)
                                ax.set_xlabel("Cells (sorted by first TF activity)")
                                ax.set_ylabel("Transcription Factors")
                                ax.set_title("TF Activity Scores")

                                plt.colorbar(im, ax=ax, label="Activity Score")
                                plt.tight_layout()
                                st.pyplot(fig)

                    with col2:
                        if st.checkbox("Show TF Target Networks"):
                            selected_tf_for_network = st.selectbox(
                                "Select TF for network view:",
                                list(tf_targets_used.keys()),
                                help="Choose a TF to visualize its target network"
                            )

                            if selected_tf_for_network:
                                targets = tf_targets_used[selected_tf_for_network]
                                st.info(f"**{selected_tf_for_network}** has {len(targets)} target genes in the dataset:")
                                st.text(", ".join(targets[:20]))  # Show first 20 targets
                                if len(targets) > 20:
                                    st.text(f"... and {len(targets) - 20} more")

            else:
                st.warning("No embeddings found in the data. Please ensure your data has dimensionality reductions (UMAP, t-SNE, PCA, etc.)")
        else:
            st.info("Please run the analysis first in the 'Run Analysis' tab")

    # Tab 5: Download Results
    with tab5:
        st.header("Download Analysis Results")

        # Check if TFvelo analysis has been completed
        if 'adata_tfvelo' not in st.session_state or 'completed' not in st.session_state.tfvelo_results:
            st.warning("🔍 No TFvelo analysis available for download.")
            st.info("Please either:")
            st.markdown("1. **Load saved analysis**: Go to Data Input → 'Load saved TFvelo analysis for visualization'")
            st.markdown("2. **Run new analysis**: Upload data → Complete analysis in 'TFvelo Setup & Analysis' tab")
            st.stop()

        # Option to select filtered or unfiltered data
        data_options = []
        if 'adata_tfvelo' in st.session_state:
            data_options.append("Filtered (adata_tfvelo)")
        if 'adata_tfvelo_unfiltered' in st.session_state:
            data_options.append("Unfiltered (adata_tfvelo_unfiltered)")

        if len(data_options) > 1:
            selected_data = st.radio(
                "Select data to download",
                data_options,
                index=0,  # Default to filtered
                horizontal=True,
                help="Filtered: after gene filtering | Unfiltered: before gene filtering"
            )
            if "Unfiltered" in selected_data:
                adata_tfvelo = st.session_state.adata_tfvelo_unfiltered
                st.info("📋 Using **unfiltered** data (before gene filtering)")
            else:
                adata_tfvelo = st.session_state.adata_tfvelo
        else:
            adata_tfvelo = st.session_state.adata_tfvelo

        if 'completed' in st.session_state.tfvelo_results:
            # Analysis summary
            st.subheader("📊 Analysis Summary")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total Cells", adata_tfvelo.n_obs)
            with col2:
                st.metric("Total Genes", adata_tfvelo.n_vars)
            with col3:
                if 'n_TFs' in adata_tfvelo.var.columns:
                    total_tf_relations = adata_tfvelo.var['n_TFs'].sum()
                    st.metric("TF-Target Relations", int(total_tf_relations))
                else:
                    st.metric("TF-Target Relations", "N/A")

            # Available data layers
            st.subheader("📋 Available Data Layers")
            layers_info = []
            layers_info.append(f"**Expression matrix (X)**: {adata_tfvelo.X.shape}")
            if adata_tfvelo.layers:
                for layer_name, layer_data in adata_tfvelo.layers.items():
                    layers_info.append(f"**{layer_name}**: {layer_data.shape}")
            if adata_tfvelo.obsm:
                for obsm_name in adata_tfvelo.obsm.keys():
                    layers_info.append(f"**{obsm_name}** (embedding): {adata_tfvelo.obsm[obsm_name].shape}")

            for info in layers_info:
                st.markdown(f"• {info}")

            # Download options
            st.subheader("💾 Download Options")
            col1, col2 = st.columns(2)

            with col1:
                if st.button("📦 Download Complete Analysis (.h5ad)", help="Download the complete analysis results as h5ad file", key="tab4_download_h5ad"):
                    try:
                        import tempfile
                        from datetime import datetime

                        # Generate filename with timestamp
                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                        filename = f"tfvelo_analysis_{timestamp}.h5ad"

                        # Create temporary file for download
                        temp_dir = tempfile.mkdtemp()
                        temp_path = os.path.join(temp_dir, filename)

                        # Convert boolean columns to string for saving (comprehensive method)
                        adata_to_save = adata_tfvelo.copy()

                        # Convert var columns
                        for key in list(adata_to_save.var.keys()):
                            try:
                                col_data = adata_to_save.var[key]
                                if hasattr(col_data, 'dtype'):
                                    # Convert boolean types
                                    if col_data.dtype == bool or str(col_data.dtype) == 'bool':
                                        adata_to_save.var[key] = col_data.astype(str)
                                    # Convert object types that contain booleans
                                    elif col_data.dtype == 'object' and len(col_data) > 0:
                                        # Check if any values are boolean
                                        sample_vals = col_data.dropna().head(10).tolist()
                                        if any(isinstance(v, (bool, np.bool_)) for v in sample_vals):
                                            adata_to_save.var[key] = col_data.astype(str)
                                    # Convert numpy boolean types
                                    elif 'bool' in str(col_data.dtype):
                                        adata_to_save.var[key] = col_data.astype(str)
                            except Exception as convert_error:
                                st.warning(f"Could not convert var column {key}: {str(convert_error)}")
                                # Try a simple string conversion as fallback
                                try:
                                    adata_to_save.var[key] = adata_to_save.var[key].astype(str)
                                except:
                                    pass

                        # Convert obs columns too
                        for key in list(adata_to_save.obs.keys()):
                            try:
                                col_data = adata_to_save.obs[key]
                                if hasattr(col_data, 'dtype'):
                                    if col_data.dtype == bool or str(col_data.dtype) == 'bool' or 'bool' in str(col_data.dtype):
                                        adata_to_save.obs[key] = col_data.astype(str)
                                    elif col_data.dtype == 'object' and len(col_data) > 0:
                                        sample_vals = col_data.dropna().head(10).tolist()
                                        if any(isinstance(v, (bool, np.bool_)) for v in sample_vals):
                                            adata_to_save.obs[key] = col_data.astype(str)
                            except Exception as convert_error:
                                st.warning(f"Could not convert obs column {key}: {str(convert_error)}")
                                try:
                                    adata_to_save.obs[key] = adata_to_save.obs[key].astype(str)
                                except:
                                    pass

                        # Save the data to temp file
                        adata_to_save.write(temp_path)

                        # Read file for download
                        with open(temp_path, 'rb') as f:
                            st.download_button(
                                label=f"📥 Download {filename}",
                                data=f.read(),
                                file_name=filename,
                                mime='application/octet-stream',
                                key="download_h5ad"
                            )

                        # Clean up temp file
                        os.remove(temp_path)
                        os.rmdir(temp_dir)
                        st.success("✅ Download prepared!")

                    except Exception as save_error:
                        st.error(f"Failed to prepare download: {str(save_error)}")

            with col2:
                if st.button("📋 Download Summary Report (.txt)", help="Download analysis summary as text file", key="tab4_download_report"):
                    try:
                        from datetime import datetime
                        import tempfile

                        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

                        # Create report content
                        report_content = []
                        report_content.append("=== TFVELO ANALYSIS REPORT ===\n")
                        report_content.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                        # Basic statistics
                        report_content.append("=== BASIC STATISTICS ===\n")
                        report_content.append(f"Total cells: {adata_tfvelo.n_obs}\n")
                        report_content.append(f"Total genes: {adata_tfvelo.n_vars}\n")

                        # TF information
                        if 'n_TFs' in adata_tfvelo.var.columns:
                            total_tfs = int(adata_tfvelo.var['n_TFs'].sum())
                            avg_tfs = adata_tfvelo.var['n_TFs'].mean()
                            report_content.append(f"Total TF-target relationships: {total_tfs}\n")
                            report_content.append(f"Average TFs per gene: {avg_tfs:.2f}\n")

                        # Available layers
                        report_content.append("\n=== AVAILABLE DATA LAYERS ===\n")
                        report_content.append(f"Expression matrix (X): {adata_tfvelo.X.shape}\n")
                        if adata_tfvelo.layers:
                            for layer_name, layer_data in adata_tfvelo.layers.items():
                                report_content.append(f"{layer_name}: {layer_data.shape}\n")

                        # Embeddings
                        if adata_tfvelo.obsm:
                            report_content.append("\n=== AVAILABLE EMBEDDINGS ===\n")
                            for obsm_name in adata_tfvelo.obsm.keys():
                                report_content.append(f"{obsm_name}: {adata_tfvelo.obsm[obsm_name].shape}\n")

                        # Gene information
                        if adata_tfvelo.var.columns.tolist():
                            report_content.append("\n=== GENE METADATA COLUMNS ===\n")
                            for col in adata_tfvelo.var.columns:
                                report_content.append(f"- {col}\n")

                        # Cell information
                        if adata_tfvelo.obs.columns.tolist():
                            report_content.append("\n=== CELL METADATA COLUMNS ===\n")
                            for col in adata_tfvelo.obs.columns:
                                report_content.append(f"- {col}\n")

                        # Create report string
                        report_text = "".join(report_content)

                        # Provide download
                        filename = f"tfvelo_report_{timestamp}.txt"
                        st.download_button(
                            label=f"📄 Download {filename}",
                            data=report_text,
                            file_name=filename,
                            mime='text/plain',
                            key="download_report"
                        )
                        st.success("✅ Report prepared!")

                    except Exception as report_error:
                        st.error(f"Failed to generate report: {str(report_error)}")

            # Additional download options
            st.subheader("🔍 Export Specific Data")

            st.info("💡 **Gene Expression Matrix**: Contains the preprocessed expression data (adata.X) used in TFvelo analysis - includes normalized, filtered genes and cells.")
            st.info("🎯 **TF Weights Matrix**: The core TFvelo result showing TF regulatory weights for each gene (positive=activation, negative=repression).")

            if st.button("📊 Download Gene Expression Matrix (.tsv)", help="Download the processed expression matrix used in TFvelo analysis", key="tab4_download_expression"):
                try:
                    import pandas as pd
                    from datetime import datetime

                    # Convert to DataFrame
                    if hasattr(adata_tfvelo.X, 'todense'):
                        expression_df = pd.DataFrame(
                            adata_tfvelo.X.todense(),
                            index=adata_tfvelo.obs_names,
                            columns=adata_tfvelo.var_names
                        )
                    else:
                        expression_df = pd.DataFrame(
                            adata_tfvelo.X,
                            index=adata_tfvelo.obs_names,
                            columns=adata_tfvelo.var_names
                        )

                    # Convert to TSV
                    tsv_string = expression_df.to_csv(sep='\t')

                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"expression_matrix_{timestamp}.tsv"

                    st.download_button(
                        label=f"📥 Download {filename}",
                        data=tsv_string,
                        file_name=filename,
                        mime='text/tab-separated-values',
                        key="download_expression"
                    )
                    st.success(f"✅ TFvelo processed expression matrix prepared! ({expression_df.shape[0]} cells × {expression_df.shape[1]} genes)")

                except Exception as tsv_error:
                    st.error(f"Failed to export expression matrix: {str(tsv_error)}")

            # TF Weights Matrix download
            if 'fit_weights' in adata_tfvelo.varm and st.button("🎯 Download TF Weights Matrix (.tsv)", help="Download TF regulatory weights matrix (core TFvelo result)", key="tab4_download_tf_weights"):
                try:
                    import pandas as pd
                    from datetime import datetime

                    # Get actual TF names from TFvelo data structures
                    tf_names = []
                    n_tfs = adata_tfvelo.varm['fit_weights'].shape[1]

                    # Method 0: Try to get TF names using correct weight matrix ordering
                    # Based on TFvelo source: TFs_x = adata.layers[layer2use][:, adata.varm['TFs_id'][gene_id][:self.n_TFs]]
                    if 'all_TFs' in adata_tfvelo.uns and 'TFs_id' in adata_tfvelo.varm:
                        try:
                            all_tfs_list = adata_tfvelo.uns['all_TFs']  # All TFs in processing order
                            tf_ids_matrix = adata_tfvelo.varm['TFs_id']  # TF indices for each gene
                            tf_names_matrix = adata_tfvelo.varm['TFs'] if 'TFs' in adata_tfvelo.varm else None


                            # Find a representative gene to get TF ordering (use first gene with maximum TFs)
                            max_tfs_per_gene = [np.sum(row[:n_tfs] >= 0) for row in tf_ids_matrix]
                            representative_gene_idx = np.argmax(max_tfs_per_gene)

                            # Get TF indices for this gene in the order used for weight calculation
                            tf_indices = tf_ids_matrix[representative_gene_idx][:n_tfs]

                            # Map indices back to TF names using all_TFs list
                            tf_names = []
                            for i, tf_idx in enumerate(tf_indices):
                                if tf_idx >= 0 and tf_idx < len(all_tfs_list):
                                    tf_name = all_tfs_list[int(tf_idx)]
                                    tf_names.append(tf_name)

                            if len(tf_names) > 0:
                                st.success(f"🎯 Found {len(tf_names)} TF names using TFs_id matrix")

                        except Exception as tfs_id_error:
                            st.warning(f"Could not extract TF names from TFs_id matrix: {str(tfs_id_error)}")

                            # Fallback to simple uns['all_TFs'] order
                            try:
                                all_tfs_ordered = adata_tfvelo.uns['all_TFs']
                                tf_names = all_tfs_ordered[:n_tfs]
                                st.info(f"🎯 Using fallback from uns['all_TFs']")
                            except Exception as uns_error:
                                pass

                    # Method 1: Try to get TF names from adata.varm['TFs'] - FALLBACK
                    if not tf_names and 'TFs' in adata_tfvelo.varm:
                        try:
                            # TFvelo stores weights in the same order as TFs are processed
                            # We need to maintain the correspondence between weight columns and TF order

                            # Get TF order by examining the TFs matrix structure
                            tf_matrix = adata_tfvelo.varm['TFs']

                            # Method 1a: Try to get consistent TF ordering from first gene with most TFs
                            max_tf_gene_idx = np.argmax([np.sum(row != '') if hasattr(row, '__iter__') else 1 for row in tf_matrix])
                            first_gene_tfs = tf_matrix[max_tf_gene_idx]

                            tf_names = []
                            if hasattr(first_gene_tfs, '__iter__') and not isinstance(first_gene_tfs, str):
                                for tf in first_gene_tfs[:n_tfs]:
                                    if tf and str(tf) != 'nan' and str(tf) != '':
                                        tf_names.append(str(tf))

                            # Ensure we have the right number of TF names
                            if len(tf_names) >= n_tfs:
                                tf_names = tf_names[:n_tfs]

                        except Exception as tf_extract_error:
                            st.warning(f"Could not extract TF names from varm['TFs']: {str(tf_extract_error)}")

                    # Method 2: Try session state filtered targets
                    if not tf_names and 'filtered_tf_targets' in st.session_state.tfvelo_results:
                        try:
                            tf_targets = st.session_state.tfvelo_results['filtered_tf_targets']
                            unique_tfs = sorted(set([target['TF'] for target in tf_targets]))
                            tf_names = unique_tfs[:n_tfs]
                            st.info(f"🎯 Found {len(tf_names)} TF names from session results")
                        except Exception as session_error:
                            st.warning(f"Could not extract TF names from session: {str(session_error)}")

                    # Method 3: Check if TF names are stored in var columns
                    if not tf_names:
                        tf_cols = [col for col in adata_tfvelo.var.columns if 'TF' in col and 'name' in col.lower()]
                        if tf_cols:
                            try:
                                tf_names = adata_tfvelo.var[tf_cols[0]].dropna().unique()[:n_tfs].tolist()
                                st.info(f"🎯 Found {len(tf_names)} TF names from var columns")
                            except Exception:
                                pass

                    # Method 4: Try to extract from database files used during preprocessing
                    if not tf_names:
                        try:
                            # Check if TF databases info is stored somewhere
                            if hasattr(adata_tfvelo, 'uns') and 'tf_databases' in adata_tfvelo.uns:
                                tf_db_info = adata_tfvelo.uns['tf_databases']
                                st.info(f"Found TF database info: {tf_db_info}")

                            # Check session state for TF database information
                            if 'tfvelo_results' in st.session_state and 'tf_databases' in st.session_state.tfvelo_results:
                                used_databases = st.session_state.tfvelo_results['tf_databases']
                                st.info(f"Used TF databases: {used_databases}")
                        except Exception:
                            pass

                    # Create proper weight matrix with correct TF indexing for each gene
                    st.info("🔧 Creating gene-specific TF weight matrix...")

                    # Get necessary data
                    weights_matrix = adata_tfvelo.varm['fit_weights']
                    tf_names_matrix = adata_tfvelo.varm['TFs']
                    gene_names = adata_tfvelo.var_names

                    # Create a list to store all TF-weight entries
                    tf_weight_entries = []

                    # Process each gene individually
                    for gene_idx, gene_name in enumerate(gene_names):
                        gene_tf_names = tf_names_matrix[gene_idx]
                        gene_weights = weights_matrix[gene_idx]

                        # Extract valid TF entries (non-blank, non-zero weight)
                        for col_idx in range(len(gene_tf_names)):
                            tf_name = gene_tf_names[col_idx]
                            weight = gene_weights[col_idx]

                            # Only include non-blank TFs with non-zero weights
                            if (tf_name != 'blank' and
                                str(tf_name) != 'nan' and
                                str(tf_name) != '' and
                                abs(weight) > 1e-10):  # Consider very small weights as effectively zero

                                tf_weight_entries.append({
                                    'Gene': gene_name,
                                    'TF': str(tf_name),
                                    'Weight': weight,
                                    'Gene_Index': gene_idx,
                                    'Column_Index': col_idx
                                })

                    # Create DataFrame from entries
                    tf_weight_df = pd.DataFrame(tf_weight_entries)

                    if len(tf_weight_df) > 0:
                        # Create pivot table: genes as rows, TFs as columns
                        weights_df = tf_weight_df.pivot_table(
                            index='Gene',
                            columns='TF',
                            values='Weight',
                            fill_value=0.0
                        )

                        # Sort columns (TFs) alphabetically for consistency
                        weights_df = weights_df.reindex(sorted(weights_df.columns), axis=1)

                        # Ensure gene order matches original
                        weights_df = weights_df.reindex(gene_names, fill_value=0.0)

                        st.success(f"✅ Created proper weight matrix: {weights_df.shape[0]} genes × {weights_df.shape[1]} TFs")

                        # Show summary statistics
                        with st.expander("📊 Weight Matrix Summary"):
                            st.write(f"**Dimensions:** {weights_df.shape[0]} genes × {weights_df.shape[1]} TFs")
                            st.write(f"**Total non-zero entries:** {len(tf_weight_df)}")
                            st.write(f"**Unique TFs:** {len(weights_df.columns)}")
                            st.write(f"**TFs per gene (avg):** {len(tf_weight_df) / len(gene_names):.1f}")

                            # Top/Bottom TFs by total weight sum
                            tf_weight_sums = weights_df.sum(axis=0).sort_values(ascending=False)
                            st.write("**Top 10 TFs by total weight sum:**")
                            for tf, total_weight in tf_weight_sums.head(10).items():
                                st.write(f"  - {tf}: {total_weight:.3f}")

                            st.write("**Bottom 5 TFs by total weight sum:**")
                            for tf, total_weight in tf_weight_sums.tail(5).items():
                                st.write(f"  - {tf}: {total_weight:.3f}")

                            # Top/Bottom TF-gene combinations by weight value
                            tf_gene_weights = []
                            for gene in weights_df.index:
                                for tf in weights_df.columns:
                                    weight = weights_df.loc[gene, tf]
                                    if abs(weight) > 1e-10:  # Only non-zero weights
                                        tf_gene_weights.append({
                                            'TF': tf,
                                            'Gene': gene,
                                            'Weight': weight
                                        })

                            tf_gene_df = pd.DataFrame(tf_gene_weights)
                            if len(tf_gene_df) > 0:
                                tf_gene_sorted = tf_gene_df.sort_values('Weight', ascending=False)

                                st.write("**Top 10 TF-Gene combinations (highest activation):**")
                                for _, row in tf_gene_sorted.head(10).iterrows():
                                    st.write(f"  - {row['TF']} → {row['Gene']}: {row['Weight']:.6f}")

                                st.write("**Bottom 10 TF-Gene combinations (strongest repression):**")
                                for _, row in tf_gene_sorted.tail(10).iterrows():
                                    st.write(f"  - {row['TF']} → {row['Gene']}: {row['Weight']:.6f}")

                            # Show sample of the matrix (top weight TFs)
                            top_weight_tfs = tf_weight_sums.head(10).index
                            st.write("**Sample of weight matrix (first 5 genes, top 10 TFs by weight sum):**")
                            sample_df = weights_df.loc[:, top_weight_tfs].iloc[:5]
                            st.dataframe(sample_df)

                    else:
                        st.error("❌ No valid TF-weight entries found!")
                        # Fallback to original method
                        tf_names = [f"TF_{i+1}" for i in range(n_tfs)]
                        weights_df = pd.DataFrame(
                            adata_tfvelo.varm['fit_weights'],
                            index=adata_tfvelo.var_names,
                            columns=tf_names
                        )

                    # Convert to TSV
                    tsv_string = weights_df.to_csv(sep='\t')

                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"tf_weights_matrix_{timestamp}.tsv"

                    st.download_button(
                        label=f"📥 Download {filename}",
                        data=tsv_string,
                        file_name=filename,
                        mime='text/tab-separated-values',
                        key="download_tf_weights"
                    )
                    # Show TF names information
                    if any('TF_' in name for name in tf_names):
                        st.warning(f"⚠️ TF weights matrix prepared with generic names! ({weights_df.shape[0]} genes × {weights_df.shape[1]} TFs)")

                        # Show debug information to help find TF names
                        with st.expander("🔍 Debug: Available data structures", expanded=False):
                            st.write("**Available varm keys:**", list(adata_tfvelo.varm.keys()))
                            st.write("**Available var columns:**", list(adata_tfvelo.var.columns))
                            st.write("**Available uns keys:**", list(adata_tfvelo.uns.keys()) if hasattr(adata_tfvelo, 'uns') else "No uns")

                            if 'TFs' in adata_tfvelo.varm:
                                tf_matrix = adata_tfvelo.varm['TFs']
                                st.write(f"**TFs varm shape:**", tf_matrix.shape)
                                st.write(f"**fit_weights shape:**", adata_tfvelo.varm['fit_weights'].shape)

                                # Show first few genes TF assignments
                                st.write("**TF assignments (first 3 genes):**")
                                for i in range(min(3, len(tf_matrix))):
                                    gene_name = adata_tfvelo.var_names[i]
                                    gene_tfs = tf_matrix[i]
                                    if hasattr(gene_tfs, '__iter__'):
                                        tfs_list = [str(tf) for tf in gene_tfs if tf and str(tf) != 'nan' and str(tf) != ''][:5]
                                        st.write(f"  - {gene_name}: {tfs_list}")
                                    else:
                                        st.write(f"  - {gene_name}: {gene_tfs}")

                                # Show max TF gene
                                if len(tf_matrix) > 0:
                                    max_tf_counts = [np.sum(row != '') if hasattr(row, '__iter__') else 1 for row in tf_matrix]
                                    max_idx = np.argmax(max_tf_counts)
                                    st.write(f"**Gene with most TFs:** {adata_tfvelo.var_names[max_idx]} ({max_tf_counts[max_idx]} TFs)")

                            # Show weight matrix info
                            if 'fit_weights' in adata_tfvelo.varm:
                                weights = adata_tfvelo.varm['fit_weights']
                                st.write(f"**Weight statistics:**")
                                st.write(f"  - Shape: {weights.shape}")
                                st.write(f"  - Range: [{weights.min():.3f}, {weights.max():.3f}]")
                                st.write(f"  - Non-zero weights: {np.sum(weights != 0)} / {weights.size}")

                        with st.expander("⚠️ TF Name-Weight Correspondence Issue", expanded=False):
                            st.warning("Using generic TF names (TF_1, TF_2, ...) because actual TF names could not be reliably mapped to weight matrix columns.")
                            st.write("**Possible causes:**")
                            st.write("- TF names stored in different format than expected")
                            st.write("- Weight matrix column order doesn't match TF discovery order")
                            st.write("- Missing TF metadata in analysis results")
                            st.write("**Impact:** Weight values are correct, but column headers are generic")
                    else:
                        st.success(f"✅ TF weights matrix prepared with actual TF names! ({weights_df.shape[0]} genes × {weights_df.shape[1]} TFs)")

                except Exception as weights_error:
                    st.error(f"Failed to export TF weights matrix: {str(weights_error)}")

            if 'velocity' in adata_tfvelo.layers and st.button("🏃 Download Velocity Data (.tsv)", help="Download velocity estimates as TSV", key="tab4_download_velocity"):
                try:
                    import pandas as pd
                    from datetime import datetime

                    # Convert velocity layer to DataFrame
                    if hasattr(adata_tfvelo.layers['velocity'], 'todense'):
                        velocity_df = pd.DataFrame(
                            adata_tfvelo.layers['velocity'].todense(),
                            index=adata_tfvelo.obs_names,
                            columns=adata_tfvelo.var_names
                        )
                    else:
                        velocity_df = pd.DataFrame(
                            adata_tfvelo.layers['velocity'],
                            index=adata_tfvelo.obs_names,
                            columns=adata_tfvelo.var_names
                        )

                    # Convert to TSV
                    tsv_string = velocity_df.to_csv(sep='\t')

                    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
                    filename = f"velocity_data_{timestamp}.tsv"

                    st.download_button(
                        label=f"📥 Download {filename}",
                        data=tsv_string,
                        file_name=filename,
                        mime='text/tab-separated-values',
                        key="download_velocity"
                    )
                    st.success(f"✅ Velocity data prepared! ({velocity_df.shape[0]} cells × {velocity_df.shape[1]} genes)")

                except Exception as vel_error:
                    st.error(f"Failed to export velocity data: {str(vel_error)}")
        else:
            st.info("Please complete the TFvelo analysis first to enable downloads.")

else:
    # No data loaded
    st.info("👈 Please upload an h5ad file in the sidebar to begin analysis")

    # Show example usage
    with st.expander("📖 Usage Guide"):
        st.markdown("""
        ### How to use TFvelo Analysis

        1. **Prepare your data**: Ensure your h5ad file contains:
           - Spliced and unspliced count matrices
           - Dimensionality reductions (UMAP, t-SNE, PCA)
           - Cell type annotations (optional)

        2. **Upload file**: Use the sidebar to upload your h5ad file

        3. **Configure settings**: Adjust parameters in the TFvelo Setup tab

        4. **Run analysis**: Execute the velocity analysis

        5. **Visualize results**: Explore velocity streams and gene-specific patterns

        6. **Export**: Save your results for further analysis

        ### About TFvelo

        TFvelo improves upon standard RNA velocity by incorporating gene regulatory
        network information, providing more accurate cell fate predictions.

        Reference: [Xing et al., TFvelo](https://github.com/xiaoyeye/TFvelo)
        """)

# Cleanup function
def cleanup_temp_dir():
    if st.session_state.temp_dir and os.path.exists(st.session_state.temp_dir):
        shutil.rmtree(st.session_state.temp_dir)
        st.session_state.temp_dir = None

# Register cleanup
import atexit
atexit.register(cleanup_temp_dir)