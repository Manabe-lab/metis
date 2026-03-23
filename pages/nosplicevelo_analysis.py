"""
noSpliceVelo Analysis (Subprocess version)
RNA velocity inference without unspliced/spliced separation

This app calls noSpliceVelo via subprocess using a separate conda environment.
Reference: https://github.com/Tarun-Mahajan/noSpliceVelo
"""

import streamlit as st
import scanpy as sc
import numpy as np
import os
import subprocess
import tempfile
import time
from pathlib import Path

st.set_page_config(
    page_title="noSpliceVelo Analysis",
    page_icon="🔮",
    layout="wide"
)

st.title("🔮 noSpliceVelo Analysis")

st.markdown("""
**RNA velocity inference without unspliced/spliced transcript separation**

noSpliceVelo estimates RNA velocity without requiring spliced/unspliced separation.
It infers the directionality of gene expression (upregulation/downregulation) from
the variance-mean relationship.

**Reference**: [Mahajan et al. "noSpliceVelo: RNA velocity without splicing information"](https://github.com/Tarun-Mahajan/noSpliceVelo)

### Features
- **No loom file required**: Works with standard scRNA-seq count data only
- **Two-stage model**: SCVIModified (VAE) -> noSpliceVelo (velocity)
- **Variance-mean principle**: "variance always leads, mean follows"
""")

st.warning("""
**Warning: Computation time can be long**

noSpliceVelo trains two-stage deep learning models,
which can take tens of minutes to several hours even with GPU.
""")

# Sidebar - Environment configuration
with st.sidebar:
    st.header("Environment Settings")

    # noSpliceVelo environment Python path
    default_nosplicevelo_python = str(Path.home() / "anaconda3" / "envs" / "nosplicevelo_env" / "bin" / "python")

    nosplicevelo_python = st.text_input(
        "noSpliceVelo Python path",
        value=default_nosplicevelo_python,
        help="Path to Python interpreter in noSpliceVelo conda environment"
    )

    # Check if the path exists
    if os.path.exists(nosplicevelo_python):
        st.success("Python path exists")
    else:
        st.warning("Python path not found")
        st.markdown("""
        **To create noSpliceVelo environment:**
        ```bash
        # Clone repository
        git clone https://github.com/Tarun-Mahajan/noSpliceVelo.git
        cd noSpliceVelo

        # Create environment (GPU recommended)
        conda env create -f environment_gpu.yml
        # Or for CPU only:
        # conda env create -f environment_cpu.yml

        conda activate nosplicevelo_env
        ```
        """)

    st.markdown("---")
    st.header("About noSpliceVelo")
    st.markdown("""
    **Principle:**
    Variance always leads and mean follows
    in gene expression dynamics.

    **Advantages:**
    - No loom file needed
    - Works with standard scRNA-seq
    - Handles complex dynamics
    - Deep learning based
    """)

# Get script path
script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
nosplicevelo_script = os.path.join(script_dir, "scripts", "run_nosplicevelo.py")

if not os.path.exists(nosplicevelo_script):
    st.error(f"noSpliceVelo script not found: {nosplicevelo_script}")
    st.info("Please ensure the script exists at the specified path.")
    st.stop()

# File upload
st.header("Step 1: Upload Data")

st.info("""
**Input data requirements:**
- h5ad file (standard scRNA-seq count data)
- spliced/unspliced layers are **not required**
- Cluster/cell type annotations are required
""")

uploaded_file = st.file_uploader(
    "Upload h5ad file",
    type=["h5ad"],
    help="AnnData file with count matrix and cell annotations"
)

if uploaded_file is not None:
    # Create temp directory
    temp_dir = tempfile.mkdtemp(prefix="nosplicevelo_")
    input_path = os.path.join(temp_dir, "input.h5ad")
    output_path = os.path.join(temp_dir, "output.h5ad")

    # Generate output filename
    original_name = uploaded_file.name
    if original_name.endswith('.h5ad'):
        output_filename = original_name[:-5] + ".noSpliceVelo.h5ad"
    else:
        output_filename = original_name + ".noSpliceVelo.h5ad"

    # Save uploaded file
    with open(input_path, "wb") as f:
        f.write(uploaded_file.getvalue())

    # Load data for preview
    with st.spinner("Loading data..."):
        adata = sc.read_h5ad(input_path)

    st.success(f"Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Show data info
    col1, col2, col3 = st.columns(3)
    with col1:
        st.subheader("Layers")
        if adata.layers:
            for layer in list(adata.layers.keys())[:10]:
                st.write(f"- {layer}")
        else:
            st.write("- X (main matrix)")

    with col2:
        st.subheader("Embeddings")
        embeddings = [k for k in adata.obsm.keys() if k.startswith("X_")]
        if embeddings:
            for emb in embeddings[:5]:
                st.write(f"- {emb}")
        else:
            st.warning("No embeddings (will be computed)")

    with col3:
        st.subheader("Metadata")
        categorical_cols = [col for col in adata.obs.columns
                          if adata.obs[col].dtype.name == 'category' or
                          adata.obs[col].dtype == 'object']
        if categorical_cols:
            for col in categorical_cols[:10]:
                st.write(f"- {col}")
        else:
            st.warning("No categorical columns")

    # Note about spliced/unspliced
    has_spliced = 'spliced' in adata.layers
    has_unspliced = 'unspliced' in adata.layers

    if has_spliced and has_unspliced:
        st.info("spliced/unspliced layers detected, but noSpliceVelo does not use them.")
    else:
        st.success("spliced/unspliced layers are not required (a feature of noSpliceVelo)")

    # Settings
    st.header("Step 2: Settings")

    # Get cluster/celltype columns
    cluster_cols = [col for col in adata.obs.columns
                   if adata.obs[col].dtype.name == 'category' or
                   adata.obs[col].dtype == 'object']

    if not cluster_cols:
        st.error("No categorical columns found for cell type labels.")
        st.error("noSpliceVelo requires cell type/cluster annotations.")
        st.stop()

    # Find default cluster column
    priority_names = ['celltype', 'cell_type', 'leiden', 'louvain', 'clusters', 'seurat_clusters']
    default_col = cluster_cols[0]
    for name in priority_names:
        if name in cluster_cols:
            default_col = name
            break

    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Basic Settings")

        # Label column selection
        label_col = st.selectbox(
            "Cell type/cluster column (required)",
            cluster_cols,
            index=cluster_cols.index(default_col) if default_col in cluster_cols else 0,
            help="Column containing cell type or cluster annotations"
        )

        unique_labels = adata.obs[label_col].unique()
        st.info(f"Found {len(unique_labels)} unique labels")

        # Number of highly variable genes
        n_top_genes = st.number_input(
            "Number of HVGs",
            value=2000,
            min_value=500,
            max_value=5000,
            step=500,
            help="Number of highly variable genes to use"
        )

        # GPU selection
        st.subheader("Computation Device")
        from gpu_utils import render_simple_gpu_selector
        gpu_value = render_simple_gpu_selector(key_prefix="nosplicevelo")

    with col2:
        st.subheader("Model Parameters")

        # Model architecture
        n_hidden = st.number_input(
            "Hidden layer dimension",
            value=128,
            min_value=32,
            max_value=512,
            step=32,
            help="Dimension of hidden layers in neural networks"
        )

        n_latent = st.number_input(
            "Latent dimension",
            value=10,
            min_value=5,
            max_value=50,
            step=5,
            help="Dimension of latent space"
        )

        tmax = st.number_input(
            "Maximum time (tmax)",
            value=24.0,
            min_value=1.0,
            max_value=100.0,
            step=1.0,
            help="Maximum time for velocity model"
        )

        # Training parameters
        st.subheader("Training")

        with st.expander("Training Parameters", expanded=False):
            batch_size = st.number_input(
                "Batch size",
                value=512,
                min_value=64,
                max_value=2048,
                step=64,
                help="Batch size for training"
            )

            scvi_epochs = st.number_input(
                "Stage 1 max epochs (SCVIModified)",
                value=10000,
                min_value=1000,
                max_value=50000,
                step=1000,
                help="Maximum epochs for first stage (VAE)"
            )

            velo_epochs = st.number_input(
                "Stage 2 max epochs (noSpliceVelo)",
                value=100000,
                min_value=10000,
                max_value=500000,
                step=10000,
                help="Maximum epochs for second stage (velocity)"
            )

            early_patience = st.number_input(
                "Stage 1 early stopping patience",
                value=100,
                min_value=10,
                max_value=500,
                step=10
            )

            velo_patience = st.number_input(
                "Stage 2 early stopping patience",
                value=1500,
                min_value=100,
                max_value=5000,
                step=100
            )

    # Run analysis
    st.header("Step 3: Run Analysis")

    # Estimate computation time
    n_cells = adata.n_obs
    device_str = "CPU" if gpu_value == -1 else f"GPU:{gpu_value}"

    st.info(f"""
    **Configuration Summary:**
    - Cells: {n_cells:,}
    - Genes: {adata.n_vars:,} (will filter to {n_top_genes} HVGs)
    - Label column: {label_col}
    - Device: {device_str}
    - Training: Stage1 {scvi_epochs:,} epochs + Stage2 {velo_epochs:,} epochs

    **Estimated time:** {device_str} mode, depending on early stopping
    """)

    if st.button("Run noSpliceVelo Analysis", type="primary"):
        if not os.path.exists(nosplicevelo_python):
            st.error(f"noSpliceVelo Python not found: {nosplicevelo_python}")
            st.error("Please configure the correct path in the sidebar.")
            st.stop()

        progress_bar = st.progress(0)
        status_text = st.empty()

        # Build command
        cmd = [
            nosplicevelo_python,
            nosplicevelo_script,
            input_path,
            output_path,
            "--label", label_col,
            "--n-top-genes", str(n_top_genes),
            "--gpu", str(gpu_value),
            "--n-hidden", str(n_hidden),
            "--n-latent", str(n_latent),
            "--tmax", str(tmax),
            "--batch-size", str(batch_size),
            "--scvi-epochs", str(scvi_epochs),
            "--velo-epochs", str(velo_epochs),
            "--early-stopping-patience", str(early_patience),
            "--velo-patience", str(velo_patience),
        ]

        status_text.text("Starting noSpliceVelo analysis...")
        progress_bar.progress(5)

        with st.expander("Command", expanded=False):
            st.code(' '.join(cmd))

        # Run subprocess
        try:
            log_lines = []
            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1
            )

            # Stream output
            with st.expander("noSpliceVelo Log", expanded=True):
                log_placeholder = st.empty()

                for line in iter(process.stdout.readline, ''):
                    log_lines.append(line.strip())
                    log_placeholder.code('\n'.join(log_lines[-50:]))

                    # Update progress based on log
                    line_lower = line.lower()
                    if "stage 0" in line_lower or "preprocessing" in line_lower:
                        progress_bar.progress(10)
                        status_text.text("Stage 0: Preprocessing...")
                    elif "stage 1" in line_lower:
                        progress_bar.progress(20)
                        status_text.text("Stage 1: Training SCVIModified (VAE)...")
                    elif "scvimodified training complete" in line_lower:
                        progress_bar.progress(40)
                        status_text.text("Stage 1 complete, extracting parameters...")
                    elif "stage 2" in line_lower:
                        progress_bar.progress(50)
                        status_text.text("Stage 2: Training noSpliceVelo...")
                    elif "nosplicevelo training complete" in line_lower:
                        progress_bar.progress(80)
                        status_text.text("Stage 2 complete, post-processing...")
                    elif "stage 3" in line_lower:
                        progress_bar.progress(85)
                        status_text.text("Stage 3: Computing velocity graph...")
                    elif "saving output" in line_lower:
                        progress_bar.progress(95)
                        status_text.text("Saving results...")

                process.wait()

            if process.returncode != 0:
                st.error(f"noSpliceVelo failed with return code {process.returncode}")
                st.error("Check the log above for details.")

                # Common error hints
                if any("cuda" in line.lower() or "gpu" in line.lower() for line in log_lines):
                    st.warning("""
                    **GPU Error Detected**

                    Try:
                    1. Use CPU mode (GPU: -1)
                    2. Check CUDA installation
                    3. Reduce batch size
                    """)
                st.stop()

            progress_bar.progress(100)
            status_text.text("Complete!")

            # Load result
            if os.path.exists(output_path):
                st.success("noSpliceVelo analysis complete!")

                result_adata = sc.read_h5ad(output_path)
                st.session_state.nosplicevelo_result = result_adata

                # Show result info
                st.subheader("Result Summary")
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.metric("Cells", result_adata.n_obs)
                    st.metric("Genes", result_adata.n_vars)
                with col2:
                    st.write("**Velocity layers:**")
                    for layer in result_adata.layers.keys():
                        if 'velocity' in layer.lower():
                            st.write(f"- {layer}")
                with col3:
                    st.write("**Computed:**")
                    if 'velocity_pseudotime' in result_adata.obs:
                        st.write("- Velocity pseudotime")
                    if 'velocity_confidence' in result_adata.obs:
                        st.write("- Velocity confidence")
                    if 'X_umap' in result_adata.obsm:
                        st.write("- UMAP embedding")

                # Download
                st.header("Step 4: Download Results")

                with open(output_path, "rb") as f:
                    st.download_button(
                        label=f"Download {output_filename}",
                        data=f.read(),
                        file_name=output_filename,
                        mime="application/octet-stream",
                        type="primary"
                    )

                st.info("""
                **Result contains:**
                - `velocity`: noSpliceVelo velocity estimates
                - `velocity_nosplicevelo`: Copy of velocity
                - `mu_nosplicevelo`: Mean expression estimates
                - `velocity_graph`: Velocity transition matrix
                - `velocity_pseudotime`: Pseudotime estimates

                **Next steps:**
                - Load in **scVelo visualization** to visualize velocity streams
                - Compare with scVelo/UniTVelo results if available
                """)

            else:
                st.error("Output file not found. Check the log for errors.")

        except Exception as e:
            st.error(f"Error running noSpliceVelo: {str(e)}")
            import traceback
            st.code(traceback.format_exc())

else:
    st.info("Please upload an h5ad file to begin")

    st.markdown("""
    ### Advantages of noSpliceVelo

    | Feature | scVelo | noSpliceVelo |
    |------|--------|--------------|
    | Input data | spliced + unspliced | Count data only |
    | Loom file | Required | **Not required** |
    | Principle | Transcriptional dynamics model | Variance-mean relationship |

    ### Workflow

    1. **Prepare h5ad file** from standard scRNA-seq analysis
    2. Add cluster/cell type annotations
    3. **Run velocity analysis** with this app
    """)

# Footer
st.markdown("---")
st.caption("""
**Note**: noSpliceVelo runs in a separate conda environment via subprocess.
Configure the Python path in the sidebar if needed.
Computation may take several hours depending on dataset size and hardware.
""")
