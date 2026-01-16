"""
COMPASS - Metabolic Flux Analysis
Characterizing metabolic pathway activity states in single cells
"""

import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import os
import tempfile
import time
import subprocess
import shutil
from helper_func import clear_old_directories, clear_old_files, check_species_index

st.set_page_config(page_title="COMPASS Analysis", page_icon="🧬", layout="wide")

st.title("🧬 COMPASS - Metabolic Flux Analysis")
st.markdown("""
[COMPASS](https://github.com/wagnerlab-berkeley/Compass) is an algorithm to characterize
the metabolic state of cells based on single-cell RNA-seq and flux balance analysis (FBA).

### Features
- Estimates metabolic reaction activity from gene expression
- Supports both human and mouse data
- Uses genome-scale metabolic models (RECON2, Human1, Mouse1)

### Requirements
- **Gurobi License**: COMPASS requires a Gurobi WLS license (free for academic use)
- Install: `pip install git+https://github.com/wagnerlab-berkeley/Compass.git`

### References
- Wagner et al. (2021) "Metabolic modeling of single Th17 cells reveals regulators of autoimmunity" Cell
""")

# Check if COMPASS is installed
def check_compass_installed():
    """Check if COMPASS is available in PATH"""
    result = subprocess.run(['which', 'compass'], capture_output=True, text=True)
    return result.returncode == 0

# Initialize session state
if "compass_temp_dir" not in st.session_state:
    compass_temp_dir = os.path.join("temp", f"compass_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(compass_temp_dir, exist_ok=True)
    st.session_state.compass_temp_dir = compass_temp_dir
else:
    compass_temp_dir = st.session_state.compass_temp_dir

if "compass_complete" not in st.session_state:
    st.session_state.compass_complete = False

# Check COMPASS installation
compass_available = check_compass_installed()

if not compass_available:
    st.error("""
    ⚠️ **COMPASS is not installed or not in PATH**

    Please install COMPASS:
    ```bash
    pip install git+https://github.com/wagnerlab-berkeley/Compass.git
    ```

    And set up Gurobi license:
    ```bash
    compass --set-license <PATH_TO_LICENSE>
    ```

    Academic licenses are available free at: https://www.gurobi.com/academia/academic-program-and-licenses/
    """)

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload h5ad file")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file (single-cell RNA-seq data)",
    type=['h5ad'],
    key="compass_h5ad_upload",
    help="AnnData h5ad file with normalized gene expression (CPM/TPM recommended)"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load and preview data
    if ("compass_adata_info" not in st.session_state or
        st.session_state.get("compass_uploaded_file") != uploaded_h5ad.name):

        with st.spinner("Reading file..."):
            # Save uploaded file
            temp_h5ad_path = os.path.join(compass_temp_dir, "input.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            # Read metadata
            adata = sc.read_h5ad(temp_h5ad_path)

            # Detect species from gene names
            gene_list = adata.var_names.tolist()
            detected_species_idx = check_species_index(gene_list)

            st.session_state.compass_adata_info = {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'obs_columns': list(adata.obs.columns),
                'gene_sample': gene_list[:10],
                'detected_species_idx': detected_species_idx
            }
            st.session_state.compass_uploaded_file = uploaded_h5ad.name
            st.session_state.compass_h5ad_path = temp_h5ad_path

    info = st.session_state.compass_adata_info

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Cells", info['n_cells'])
    with col2:
        st.metric("Genes", info['n_genes'])
    with col3:
        detected = "Human" if info['detected_species_idx'] == 1 else "Mouse"
        st.metric("Detected Species", detected)

    with st.expander("Sample genes"):
        st.write(info['gene_sample'])

    # ========================================
    # Step 2: Configure analysis parameters
    # ========================================
    st.header("Step 2: Configure analysis")

    with st.expander("📚 Parameter Guide", expanded=False):
        st.markdown("""
        ### Species
        - **Human (homo_sapiens)**: Uses RECON2 or Human1 metabolic model
        - **Mouse (mus_musculus)**: Uses Mouse1 metabolic model (MGI gene symbols)

        ### Metabolic Model
        - **RECON2**: Classic human metabolic reconstruction (default)
        - **Human1**: Newer human model from Metabolic Atlas
        - **Mouse1**: Mouse model from Metabolic Atlas

        ### Number of Processes
        - Parallel processing for faster computation
        - Recommended: number of CPU cores available

        ### Output
        - **reactions.tsv**: Reaction penalty scores (lower = more likely active)
        - **Compass score**: Alignment of cell's metabolic program with reaction flux
        """)

    with st.form("compass_params_form"):
        col1, col2 = st.columns(2)

        with col1:
            species = st.selectbox(
                "Species",
                options=["homo_sapiens", "mus_musculus"],
                index=1 - info['detected_species_idx'],  # Auto-select based on detection
                help="Species for metabolic model selection"
            )

            # Show gene symbol format hint
            if species == "homo_sapiens":
                st.caption("Gene symbols: HGNC format (e.g., GAPDH, ACTB)")
            else:
                st.caption("Gene symbols: MGI format (e.g., Gapdh, Actb)")

        with col2:
            num_processes = st.number_input(
                "Number of processes",
                min_value=1,
                max_value=32,
                value=4,
                help="Number of parallel processes (recommended: number of CPU cores)"
            )

        col1, col2 = st.columns(2)

        with col1:
            metabolic_model = st.selectbox(
                "Metabolic model",
                options=["RECON2 (default)", "Human1", "Mouse1"],
                index=0 if species == "homo_sapiens" else 2,
                help="Genome-scale metabolic model to use"
            )

        with col2:
            output_prefix = st.text_input(
                "Output prefix",
                value="compass_result",
                help="Prefix for output files"
            )

        # Advanced options
        st.subheader("Advanced options")

        col1, col2 = st.columns(2)

        with col1:
            calc_metabolites = st.checkbox(
                "Calculate metabolite scores",
                value=False,
                help="Also calculate metabolite uptake/secretion scores (slower)"
            )

        with col2:
            test_mode = st.checkbox(
                "Test mode (subset of reactions)",
                value=False,
                help="Run on subset of reactions for quick testing"
            )

        st.markdown("---")

        submit_analysis = st.form_submit_button(
            "🧬 Run COMPASS Analysis",
            type="primary",
            disabled=not compass_available
        )

    # ========================================
    # Step 3: Run analysis
    # ========================================
    if submit_analysis:
        st.header("Step 3: Running analysis")

        with st.spinner("Running COMPASS..."):
            progress_bar = st.progress(0)
            status_text = st.empty()
            log_container = st.empty()

            try:
                status_text.text("Preparing data...")
                progress_bar.progress(10)

                # Create output directory
                output_dir = os.path.join(compass_temp_dir, "output")
                os.makedirs(output_dir, exist_ok=True)

                # Build COMPASS command
                cmd = [
                    "compass",
                    "--data", st.session_state.compass_h5ad_path,
                    "--num-processes", str(num_processes),
                    "--species", species,
                    "--output-dir", output_dir
                ]

                # Add model option if not default
                if "Human1" in metabolic_model:
                    cmd.extend(["--model", "Human1"])
                elif "Mouse1" in metabolic_model:
                    cmd.extend(["--model", "Mouse1"])

                # Add optional flags
                if calc_metabolites:
                    cmd.append("--calc-metabolites")

                if test_mode:
                    cmd.append("--test-mode")

                status_text.text(f"Running: {' '.join(cmd)}")
                progress_bar.progress(20)

                # Run COMPASS
                st.info(f"Command: `{' '.join(cmd)}`")

                process = subprocess.Popen(
                    cmd,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                    text=True,
                    bufsize=1
                )

                log_lines = []
                for line in process.stdout:
                    log_lines.append(line.strip())
                    # Keep last 20 lines
                    if len(log_lines) > 20:
                        log_lines = log_lines[-20:]
                    log_container.code('\n'.join(log_lines), language='text')

                process.wait()

                if process.returncode != 0:
                    st.error(f"COMPASS failed with return code {process.returncode}")
                    st.stop()

                progress_bar.progress(80)
                status_text.text("Loading results...")

                # Check output files
                reactions_file = os.path.join(output_dir, "reactions.tsv")

                if not os.path.exists(reactions_file):
                    st.error(f"Output file not found: {reactions_file}")
                    st.write("Available files:", os.listdir(output_dir))
                    st.stop()

                # Load results
                reactions_df = pd.read_csv(reactions_file, sep='\t', index_col=0)

                st.session_state.compass_reactions = reactions_df
                st.session_state.compass_output_dir = output_dir
                st.session_state.compass_complete = True

                # Load metabolites if calculated
                metabolites_file = os.path.join(output_dir, "metabolites.tsv")
                if os.path.exists(metabolites_file):
                    metabolites_df = pd.read_csv(metabolites_file, sep='\t', index_col=0)
                    st.session_state.compass_metabolites = metabolites_df

                progress_bar.progress(100)
                status_text.text("Analysis complete!")

                st.success(f"""
                ✅ **COMPASS analysis completed!**

                - Reactions analyzed: {len(reactions_df)}
                - Cells processed: {reactions_df.shape[1]}
                """)

            except FileNotFoundError as e:
                st.error(f"COMPASS not found. Please install it first: {e}")
            except Exception as e:
                st.error(f"Error during analysis: {str(e)}")
                st.exception(e)

    # ========================================
    # Step 4: Results and Download
    # ========================================
    if st.session_state.compass_complete:
        st.header("Step 4: Results")

        reactions_df = st.session_state.compass_reactions

        # Summary statistics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Total reactions", len(reactions_df))
        with col2:
            st.metric("Cells", reactions_df.shape[1])
        with col3:
            # Count reactions with variation
            reaction_var = reactions_df.var(axis=1)
            active_reactions = (reaction_var > 0).sum()
            st.metric("Variable reactions", active_reactions)

        # Show top variable reactions
        st.subheader("Top variable reactions")

        reaction_stats = pd.DataFrame({
            'mean': reactions_df.mean(axis=1),
            'std': reactions_df.std(axis=1),
            'min': reactions_df.min(axis=1),
            'max': reactions_df.max(axis=1)
        })
        reaction_stats = reaction_stats.sort_values('std', ascending=False)

        st.dataframe(reaction_stats.head(50))

        # Reaction score distribution
        st.subheader("Score distribution")

        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 4))

        # Mean scores histogram
        axes[0].hist(reactions_df.mean(axis=1), bins=50, edgecolor='black', alpha=0.7)
        axes[0].set_xlabel('Mean COMPASS score')
        axes[0].set_ylabel('Number of reactions')
        axes[0].set_title('Distribution of mean reaction scores')

        # Variance histogram
        axes[1].hist(reactions_df.var(axis=1), bins=50, edgecolor='black', alpha=0.7)
        axes[1].set_xlabel('Variance')
        axes[1].set_ylabel('Number of reactions')
        axes[1].set_title('Distribution of reaction score variance')

        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

        # Preview data
        with st.expander("Preview reaction scores"):
            st.dataframe(reactions_df.head(100))

        # Metabolites if available
        if 'compass_metabolites' in st.session_state:
            st.subheader("Metabolite scores")
            metabolites_df = st.session_state.compass_metabolites
            st.write(f"Metabolites analyzed: {len(metabolites_df)}")

            with st.expander("Preview metabolite scores"):
                st.dataframe(metabolites_df.head(50))

        # Download section
        st.subheader("Download results")

        col1, col2, col3 = st.columns(3)

        with col1:
            # Download reactions.tsv
            reactions_csv = reactions_df.to_csv(sep='\t')
            st.download_button(
                label="⬇️ Download reactions.tsv",
                data=reactions_csv,
                file_name=f"{output_prefix}_reactions.tsv",
                mime="text/tab-separated-values",
                type="primary"
            )

        with col2:
            # Download reaction statistics
            stats_csv = reaction_stats.to_csv(sep='\t')
            st.download_button(
                label="⬇️ Download statistics.tsv",
                data=stats_csv,
                file_name=f"{output_prefix}_statistics.tsv",
                mime="text/tab-separated-values"
            )

        with col3:
            if 'compass_metabolites' in st.session_state:
                metabolites_csv = st.session_state.compass_metabolites.to_csv(sep='\t')
                st.download_button(
                    label="⬇️ Download metabolites.tsv",
                    data=metabolites_csv,
                    file_name=f"{output_prefix}_metabolites.tsv",
                    mime="text/tab-separated-values"
                )

        # Usage instructions
        st.markdown("""
        ---
        ### Interpreting Results

        **COMPASS scores** represent the "penalty" for carrying flux through each reaction:
        - **Lower scores** = reaction is more likely to be active
        - **Higher scores** = reaction is less likely to be active

        **Recommended analysis:**
        1. Identify differentially active reactions between cell types/conditions
        2. Map reactions to metabolic pathways
        3. Correlate reaction scores with gene expression

        **Resources:**
        - [Virtual Metabolic Human](https://www.vmh.life/) - RECON2 reaction information
        - [Metabolic Atlas](https://metabolicatlas.org/) - Human1/Mouse1 models
        """)

else:
    st.info("👆 Upload an h5ad file to begin COMPASS analysis")

    # Show installation instructions
    st.markdown("""
    ### Getting Started

    1. **Install COMPASS:**
    ```bash
    pip install numpy
    pip install git+https://github.com/wagnerlab-berkeley/Compass.git
    ```

    2. **Get Gurobi license** (free for academic use):
       - Visit: https://www.gurobi.com/academia/academic-program-and-licenses/
       - Download license file
       - Register: `compass --set-license <PATH_TO_LICENSE>`

    3. **Prepare data:**
       - h5ad file with normalized gene expression (CPM/TPM)
       - Gene symbols: HGNC (human) or MGI (mouse)

    4. **Run analysis** and download reaction scores
    """)

# Footer
st.markdown("---")
st.markdown("""
**Reference:** Wagner, A., et al. (2021). Metabolic modeling of single Th17 cells reveals
regulators of autoimmunity. Cell, 184(16), 4168-4185.

**Documentation:** [COMPASS Docs](https://compass-wagnerlab.readthedocs.io/)
""")
