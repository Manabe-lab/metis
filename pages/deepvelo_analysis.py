"""
DeepVelo Analysis (Subprocess version)
Deep learning-based RNA velocity analysis using DeepVelo

This app calls DeepVelo via subprocess using a separate conda environment.
"""

import streamlit as st
import scanpy as sc
import scvelo as scv
import numpy as np
import os
import subprocess
import tempfile
import time
from pathlib import Path

st.set_page_config(
    page_title="DeepVelo Analysis",
    page_icon="🧠",
    layout="wide"
)

st.title("🧠 DeepVelo Analysis")
st.markdown("""
Deep learning-based RNA velocity estimation using DeepVelo.

**Reference**: [Gao et al. (2024) Nature Communications](https://www.nature.com/articles/s41467-024-51278-6)

This app runs DeepVelo in a separate conda environment via subprocess.
""")

# Sidebar - Environment configuration
with st.sidebar:
    st.header("⚙️ Environment Settings")

    # DeepVelo environment Python path
    default_deepvelo_python = str(Path.home() / "anaconda3" / "envs" / "deepvelo_env" / "bin" / "python")

    deepvelo_python = st.text_input(
        "DeepVelo Python path",
        value=default_deepvelo_python,
        help="Path to Python interpreter in DeepVelo conda environment"
    )

    # Check if the path exists
    if os.path.exists(deepvelo_python):
        st.success("✅ Python path exists")
    else:
        st.warning("⚠️ Python path not found")
        st.markdown("""
        **To create DeepVelo environment:**
        ```bash
        conda create -n deepvelo_env python=3.10
        conda activate deepvelo_env
        pip install deepvelo scvelo scanpy
        pip install torch dgl
        ```
        """)

    st.markdown("---")
    st.header("About DeepVelo")
    st.markdown("""
    DeepVelo uses graph neural networks to estimate RNA velocity.

    **Advantages:**
    - More accurate velocity estimation
    - Better handles complex dynamics
    - Robust to noise
    """)

# Get script path
script_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
deepvelo_script = os.path.join(script_dir, "scripts", "run_deepvelo.py")

if not os.path.exists(deepvelo_script):
    st.error(f"❌ DeepVelo script not found: {deepvelo_script}")
    st.stop()

# File upload
st.header("Step 1: Upload Data")

input_type = st.radio(
    "Input data type",
    ["scVelo preprocessed h5ad", "Raw h5ad (with spliced/unspliced)"],
    help="Choose whether your data needs preprocessing or is already processed by scVelo"
)

uploaded_file = st.file_uploader(
    "Upload h5ad file",
    type=["h5ad"],
    help="AnnData file with spliced and unspliced counts"
)

if uploaded_file is not None:
    # Create temp directory
    temp_dir = tempfile.mkdtemp(prefix="deepvelo_")
    input_path = os.path.join(temp_dir, "input.h5ad")
    output_path = os.path.join(temp_dir, "output.h5ad")

    # Generate output filename: original_name.DeepVelo.h5ad
    original_name = uploaded_file.name
    if original_name.endswith('.h5ad'):
        output_filename = original_name[:-5] + ".DeepVelo.h5ad"
    else:
        output_filename = original_name + ".DeepVelo.h5ad"

    # Save uploaded file
    with open(input_path, "wb") as f:
        f.write(uploaded_file.getvalue())

    # Load data for preview
    with st.spinner("Loading data..."):
        adata = sc.read_h5ad(input_path)

    st.success(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    # Show data info
    col1, col2 = st.columns(2)
    with col1:
        st.subheader("Layers")
        if adata.layers:
            for layer in adata.layers.keys():
                st.write(f"- {layer}")
        else:
            st.warning("No layers found")

    with col2:
        st.subheader("Embeddings")
        embeddings = [k for k in adata.obsm.keys() if k.startswith("X_")]
        if embeddings:
            for emb in embeddings:
                st.write(f"- {emb}")
        else:
            st.warning("No embeddings found")

    # Check required data
    has_spliced = 'spliced' in adata.layers
    has_unspliced = 'unspliced' in adata.layers

    if not has_spliced or not has_unspliced:
        st.error("❌ Missing spliced/unspliced layers. Please provide data with RNA velocity information.")
        st.stop()

    # Preprocessing options
    st.header("Step 2: Settings")

    need_preprocess = input_type == "Raw h5ad (with spliced/unspliced)"

    col1, col2 = st.columns(2)
    with col1:
        if need_preprocess:
            st.info("Data will be preprocessed before DeepVelo analysis")
            min_shared_counts = st.number_input("Min shared counts", value=20, min_value=1)
            n_top_genes = st.number_input("Number of top genes", value=2000, min_value=100)
        else:
            st.info("Using existing preprocessing")
            min_shared_counts = 20
            n_top_genes = 2000

    with col2:
        n_pcs = st.number_input("Number of PCs", value=30, min_value=5)
        n_neighbors = st.number_input("Number of neighbors", value=30, min_value=5)

    # Run analysis
    st.header("Step 3: Run Analysis")

    if st.button("🚀 Run DeepVelo Analysis", type="primary"):
        if not os.path.exists(deepvelo_python):
            st.error(f"❌ DeepVelo Python not found: {deepvelo_python}")
            st.error("Please configure the correct path in the sidebar.")
            st.stop()

        progress_bar = st.progress(0)
        status_text = st.empty()
        log_container = st.empty()

        # Build command
        cmd = [
            deepvelo_python,
            deepvelo_script,
            input_path,
            output_path,
            "--n-pcs", str(n_pcs),
            "--n-neighbors", str(n_neighbors),
        ]

        if need_preprocess:
            cmd.extend([
                "--preprocess",
                "--min-shared-counts", str(min_shared_counts),
                "--n-top-genes", str(n_top_genes),
            ])

        status_text.text("Starting DeepVelo analysis...")
        progress_bar.progress(10)

        st.info(f"Running command:\n```\n{' '.join(cmd)}\n```")

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
            with st.expander("📋 DeepVelo Log", expanded=True):
                log_placeholder = st.empty()

                for line in iter(process.stdout.readline, ''):
                    log_lines.append(line.strip())
                    log_placeholder.code('\n'.join(log_lines[-50:]))  # Show last 50 lines

                    # Update progress based on log
                    if "preprocessing" in line.lower():
                        progress_bar.progress(30)
                        status_text.text("Preprocessing...")
                    elif "deepvelo velocity" in line.lower():
                        progress_bar.progress(50)
                        status_text.text("Running DeepVelo (this may take several minutes)...")
                    elif "velocity graph" in line.lower():
                        progress_bar.progress(80)
                        status_text.text("Computing velocity graph...")
                    elif "saving output" in line.lower():
                        progress_bar.progress(90)
                        status_text.text("Saving results...")

                process.wait()

            if process.returncode != 0:
                st.error(f"❌ DeepVelo failed with return code {process.returncode}")
                st.error("Check the log above for details.")
                st.stop()

            progress_bar.progress(100)
            status_text.text("Complete!")

            # Load result
            if os.path.exists(output_path):
                st.success("🎉 DeepVelo analysis complete!")

                result_adata = sc.read_h5ad(output_path)
                st.session_state.deepvelo_result = result_adata

                # Show result info
                st.subheader("Result Summary")
                col1, col2 = st.columns(2)
                with col1:
                    st.write(f"**Cells:** {result_adata.n_obs}")
                    st.write(f"**Genes:** {result_adata.n_vars}")
                with col2:
                    st.write("**Velocity layers:**")
                    for layer in result_adata.layers.keys():
                        if 'velocity' in layer.lower():
                            st.write(f"- {layer}")

                # Download
                st.header("Step 4: Download Results")

                with open(output_path, "rb") as f:
                    st.download_button(
                        label=f"📥 Download {output_filename}",
                        data=f.read(),
                        file_name=output_filename,
                        mime="application/octet-stream"
                    )

                st.info("""
                **Result contains:**
                - `velocity`: DeepVelo velocity estimates
                - `velocity_deepvelo`: Copy of DeepVelo velocity
                - `velocity_graph`: Velocity transition matrix
                - `velocity_confidence`: Confidence scores

                Load this file in **scVelo visualization** to compare with scVelo results.
                """)

            else:
                st.error("❌ Output file not found. Check the log for errors.")

        except Exception as e:
            st.error(f"❌ Error running DeepVelo: {str(e)}")
            import traceback
            st.code(traceback.format_exc())

else:
    st.info("👆 Please upload an h5ad file to begin")

# Footer
st.markdown("---")
st.caption("""
**Note**: DeepVelo runs in a separate conda environment via subprocess.
Configure the Python path in the sidebar if needed.
""")
