"""
scFEA - Single-cell Flux Estimation Analysis
Graph neural network model to estimate cell-wise metabolic flux
"""

import streamlit as st
import scanpy as sc
import pandas as pd
import numpy as np
import os
import sys
import tempfile
import time
import subprocess
import shutil
from pathlib import Path
from helper_func import clear_old_directories, clear_old_files, check_species_index

st.set_page_config(page_title="scFEA Analysis", page_icon="🔬", layout="wide")

st.title("🔬 scFEA - Single-cell Flux Estimation Analysis")
st.markdown("""
[scFEA](https://github.com/changwn/scFEA) is a graph neural network model to estimate
cell-wise metabolic flux using single-cell RNA-seq data.

### Features
- **Fast**: Estimation via Graph Neural Network (significantly faster than COMPASS)
- **No license required**: Open source (no Gurobi license needed)
- Outputs both Metabolic flux and Metabolite balance
- Supports Human / Mouse

### vs COMPASS
| Item | scFEA | COMPASS |
|------|-------|---------|
| Method | Graph Neural Network | Linear Programming |
| Speed | **Minutes to tens of minutes** | Hours to days |
| License | Not required | Gurobi required |
| Output | Flux + Balance | Reaction scores |
| GPU | Recommended (10k+ cells) | Not required |

### References
- Alghamdi et al. (2021) "A graph neural network model to estimate cell-wise metabolic flux using single-cell RNA-seq data" Genome Research
""")

with st.expander("📖 Difference between Flux and Balance", expanded=False):
    st.markdown("""
### Flux (Metabolic Flux)
- **Definition**: The **reaction rate** of each metabolic module (reaction pathway)
- **Unit**: Relative values (comparable between cells)
- **Meaning**: How actively a metabolic pathway is operating
- **Example**: High glycolysis flux = active conversion of glucose to pyruvate

```
Glucose → G6P → G3P → Pyruvate → Lactate
         ↑ Each arrow is a module, its reaction rate is flux
```

### Balance (Metabolite Balance/Stress)
- **Definition**: The **difference between inflow and outflow** (imbalance) of each intermediate metabolite
- **Calculation**: Balance = Sum(inflow flux) - Sum(outflow flux)
- **Meaning**:
  - **Positive value** = Metabolite accumulation (production > consumption)
  - **Negative value** = Metabolite depletion (production < consumption)
  - **Near zero** = Equilibrium state

### When to use which
| Purpose | Data to use |
|------|-----------|
| Compare metabolic pathway activity | **Flux** (168 modules) |
| Detect metabolite accumulation/depletion | **Balance** (70 metabolites) |
| Assess metabolic stress | **Balance** |
| Detect Warburg effect | **Flux** (lactate production flux increased) |
""")

with st.expander("📋 Supermodule (Metabolic Pathway Category) List", expanded=False):
    st.markdown("""
Module numbers are ordered by metabolic pathway category (Supermodule).

| Module Range | ID | Metabolic Pathway | Number of Modules |
|-----------|:--:|---------|:-----------:|
| M_1 - M_14 | 1 | **Glycolysis + TCA cycle** | 14 |
| M_15 - M_32 | 2 | Serine metabolism | 18 |
| M_33 | 3 | Pentose phosphate pathway | 1 |
| M_34 - M_35 | 4 | Fatty acids metabolism | 2 |
| M_36 - M_40 | 5 | Aspartate metabolism | 5 |
| M_41 - M_45 | 6 | Beta-alanine metabolism | 5 |
| M_46 - M_47 | 7 | Propionyl-CoA metabolism | 2 |
| M_48 - M_52 | 8 | Glutamate metabolism | 5 |
| M_53 - M_60 | 9 | Leucine + Valine + Isoleucine (branched-chain amino acids) | 8 |
| M_61 - M_68 | 10 | Urea cycle | 8 |
| M_69 - M_70 | 11 | Spermine metabolism | 2 |
| M_71 - M_105 | 12 | **Transporters** | 35 |
| M_106 - M_110 | 13 | Hyaluronic acid synthesis | 5 |
| M_111 | 14 | Glycogen synthesis | 1 |
| M_112 | 15 | Glycosaminoglycan synthesis | 1 |
| M_113 - M_124 | 16 | N-linked glycan synthesis | 12 |
| M_125 - M_128 | 17 | O-linked glycan synthesis | 4 |
| M_129 - M_131 | 18 | Sialic acid synthesis | 3 |
| M_132 | 19 | Glycan synthesis | 1 |
| M_133 - M_149 | 20 | Purine synthesis | 17 |
| M_150 - M_166 | 21 | Pyrimidine synthesis | 17 |
| M_167 - M_169 | 22 | Steroid hormone synthesis | 3 |
""")

# scFEA directory
SCFEA_DIR = Path(__file__).resolve().parent.parent / "scFEA"
SCFEA_SRC = SCFEA_DIR / "src"
SCFEA_DATA = SCFEA_DIR / "data"

# Model configurations
MODEL_CONFIGS = {
    "Human Complete (168 modules)": {
        "moduleGene": "module_gene_m168.csv",
        "stoichiometry": "cmMat_c70_m168.csv",
        "cName": "cName_c70_m168.csv",
        "species": "human",
        "description": "All metabolic modules (168) - for comprehensive analysis"
    },
    "Mouse Complete (168 modules)": {
        "moduleGene": "module_gene_complete_mouse_m168.csv",
        "stoichiometry": "cmMat_complete_mouse_c70_m168.csv",
        "cName": "cName_complete_mouse_c70_m168.csv",
        "species": "mouse",
        "description": "Mouse all metabolic modules (168)"
    },
    "Glutaminolysis 1 (23 modules)": {
        "moduleGene": "module_gene_glutaminolysis1_m23.csv",
        "stoichiometry": "cmMat_glutaminolysis1_c17_m23.csv",
        "cName": "cName_glutaminolysis1_c17_m23.csv",
        "species": "human",
        "description": "Focused on glutamine metabolism (fast)"
    },
    "Glutaminolysis 2 (20 modules)": {
        "moduleGene": "module_gene_glutaminolysis2_m20.csv",
        "stoichiometry": "cmMat_glutaminolysis2_c16_m20.csv",
        "cName": "cName_glutaminolysis2_c16_m20.csv",
        "species": "human",
        "description": "Glutamine metabolism variant"
    },
    "Iron Ion (15 modules)": {
        "moduleGene": "module_gene_iron_m15.csv",
        "stoichiometry": "cmMat_iron_c8_m15.csv",
        "cName": "cName_iron_c8_m15.csv",
        "species": "human",
        "description": "Focused on iron ion metabolism (fast)"
    }
}


def check_scfea_installed():
    """Check if scFEA and dependencies are available"""
    checks = {
        "scFEA directory": SCFEA_DIR.exists(),
        "scFEA.py": (SCFEA_SRC / "scFEA.py").exists(),
        "Model data": (SCFEA_DATA / "module_gene_m168.csv").exists()
    }
    return checks


def check_pytorch_installed():
    """Check if PyTorch is installed"""
    try:
        import torch
        return True, torch.__version__, torch.cuda.is_available()
    except ImportError:
        return False, None, False


def check_magic_installed():
    """Check if MAGIC is installed"""
    try:
        import magic
        return True
    except ImportError:
        return False


# Initialize session state
if "scfea_temp_dir" not in st.session_state:
    scfea_temp_dir = os.path.join("temp", f"scfea_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(scfea_temp_dir, exist_ok=True)
    st.session_state.scfea_temp_dir = scfea_temp_dir
else:
    scfea_temp_dir = st.session_state.scfea_temp_dir

if "scfea_complete" not in st.session_state:
    st.session_state.scfea_complete = False

# Check installation
scfea_checks = check_scfea_installed()
pytorch_installed, pytorch_version, cuda_available = check_pytorch_installed()
magic_installed = check_magic_installed()

# Show installation status
with st.expander("📋 Installation Status", expanded=not all(scfea_checks.values()) or not pytorch_installed):
    col1, col2 = st.columns(2)

    with col1:
        st.markdown("**scFEA Files:**")
        for check_name, status in scfea_checks.items():
            icon = "✅" if status else "❌"
            st.write(f"{icon} {check_name}")

    with col2:
        st.markdown("**Dependencies:**")
        if pytorch_installed:
            cuda_status = "with CUDA" if cuda_available else "CPU only"
            st.write(f"✅ PyTorch {pytorch_version} ({cuda_status})")
        else:
            st.write("❌ PyTorch not installed")

        if magic_installed:
            st.write("✅ MAGIC (imputation)")
        else:
            st.write("⚠️ MAGIC not installed (imputation disabled)")

    if not pytorch_installed:
        st.error("""
        **PyTorch is required for scFEA.**

        Install with:
        ```bash
        # CPU only
        pip install torch

        # With CUDA (recommended for large datasets)
        pip install torch --index-url https://download.pytorch.org/whl/cu118
        ```
        """)

    if not magic_installed:
        st.warning("""
        **MAGIC is optional but recommended for 10x data.**

        Install with:
        ```bash
        pip install magic-impute
        ```
        """)

if not all(scfea_checks.values()):
    st.error(f"scFEA is not properly installed at {SCFEA_DIR}")
    st.stop()

if not pytorch_installed:
    st.error("PyTorch is required. Please install it first.")
    st.stop()

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload h5ad file")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file (single-cell RNA-seq data)",
    type=['h5ad'],
    key="scfea_h5ad_upload",
    help="AnnData h5ad file with gene expression data (raw or normalized counts)"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load and preview data
    if ("scfea_adata_info" not in st.session_state or
        st.session_state.get("scfea_uploaded_file") != uploaded_h5ad.name):

        with st.spinner("Reading file..."):
            # Save uploaded file
            temp_h5ad_path = os.path.join(scfea_temp_dir, "input.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            # Read metadata
            adata = sc.read_h5ad(temp_h5ad_path)

            # Detect species from gene names
            gene_list = adata.var_names.tolist()
            detected_species_idx = check_species_index(gene_list)

            st.session_state.scfea_adata_info = {
                'n_cells': adata.n_obs,
                'n_genes': adata.n_vars,
                'obs_columns': list(adata.obs.columns),
                'gene_sample': gene_list[:10],
                'detected_species_idx': detected_species_idx
            }
            st.session_state.scfea_uploaded_file = uploaded_h5ad.name
            st.session_state.scfea_h5ad_path = temp_h5ad_path

    info = st.session_state.scfea_adata_info

    col1, col2, col3, col4 = st.columns(4)
    with col1:
        st.metric("Cells", info['n_cells'])
    with col2:
        st.metric("Genes", info['n_genes'])
    with col3:
        detected = "Human" if info['detected_species_idx'] == 1 else "Mouse"
        st.metric("Detected Species", detected)
    with col4:
        gpu_status = "GPU" if cuda_available else "CPU"
        st.metric("Compute", gpu_status)

    with st.expander("Sample genes"):
        st.write(info['gene_sample'])

    # ========================================
    # Raw Counts Layer Selection & CPM Normalization
    # ========================================
    st.subheader("📊 Raw Counts Layer Selection")
    st.markdown("""
    **scFEA recommends CPM/TPM-normalized data.**

    1. **Select the layer** containing raw counts
    2. The selected layer will be automatically **CPM-normalized** (no log transformation)
    3. If max > 50, scFEA internally applies log2 transformation
    """)

    # Load adata to check layers
    adata = sc.read_h5ad(st.session_state.scfea_h5ad_path)

    # Show current .X statistics
    if hasattr(adata.X, 'toarray'):
        x_data = adata.X.toarray()
    else:
        x_data = adata.X

    x_max = float(np.max(x_data))
    x_min = float(np.min(x_data))
    x_mean = float(np.mean(x_data))
    has_decimals = not np.allclose(x_data, np.round(x_data))
    x_sparsity = float((x_data == 0).sum() / x_data.size * 100)

    col1, col2, col3, col4, col5 = st.columns(5)
    with col1:
        st.metric(".X max", f"{x_max:.2f}")
    with col2:
        st.metric(".X min", f"{x_min:.2f}")
    with col3:
        st.metric(".X mean", f"{x_mean:.2f}")
    with col4:
        decimal_status = "Yes" if has_decimals else "No"
        st.metric("Has decimals", decimal_status)
    with col5:
        st.metric("Sparsity", f"{x_sparsity:.1f}%")

    # Build layer options
    layer_options = [
        "adata.X \u2192 CPM normalize (if raw counts)",
        "adata.X (use as-is, if already CPM/TPM normalized)"
    ]
    if hasattr(adata, 'raw') and adata.raw is not None:
        layer_options.insert(0, "adata.raw.X \u2192 CPM normalize")
    if adata.layers:
        for layer_name in adata.layers.keys():
            layer_options.insert(0, f"adata.layers['{layer_name}'] \u2192 CPM normalize")

    # Show available layers info
    with st.expander("\U0001f4cb Layer information", expanded=False):
        st.write("**adata.X**")
        st.write(f"  - Shape: {adata.X.shape}")
        st.write(f"  - Max: {x_max:.2f}, Sparsity: {x_sparsity:.1f}%, Has decimals: {has_decimals}")

        if hasattr(adata, 'raw') and adata.raw is not None:
            if hasattr(adata.raw.X, 'toarray'):
                raw_data = adata.raw.X.toarray()
            else:
                raw_data = adata.raw.X
            raw_max = float(np.max(raw_data))
            raw_sparsity = float((raw_data == 0).sum() / raw_data.size * 100)
            raw_decimals = not np.allclose(raw_data, np.round(raw_data))
            st.write("**adata.raw.X**")
            st.write(f"  - Shape: {adata.raw.X.shape}")
            st.write(f"  - Max: {raw_max:.2f}, Sparsity: {raw_sparsity:.1f}%, Has decimals: {raw_decimals}")

        if adata.layers:
            for layer_name in adata.layers.keys():
                layer_data = adata.layers[layer_name]
                if hasattr(layer_data, 'toarray'):
                    layer_arr = layer_data.toarray()
                else:
                    layer_arr = layer_data
                layer_max = float(np.max(layer_arr))
                layer_sparsity = float((layer_arr == 0).sum() / layer_arr.size * 100)
                layer_decimals = not np.allclose(layer_arr, np.round(layer_arr))
                st.write(f"**adata.layers['{layer_name}']**")
                st.write(f"  - Max: {layer_max:.2f}, Sparsity: {layer_sparsity:.1f}%, Has decimals: {layer_decimals}")

        st.markdown("""
        **Guidelines:**
        - Max > 100 and integers = Raw counts (CPM normalization needed)
        - Max < 50 and has decimals = Already normalized
        - Max approx. 10,000 - 1,000,000 = CPM/TPM normalized
        """)

    # Layer selection
    selected_layer = st.selectbox(
        "Select data source:",
        options=layer_options,
        index=0,
        key="scfea_layer_selection",
        help="Select raw counts to apply CPM normalization. If already CPM/TPM normalized, select 'use as-is'"
    )

    # Show warning/info based on selection
    if "use as-is" in selected_layer:
        if x_max > 100 and not has_decimals:
            st.warning("""
            \u26a0\ufe0f **Warning:** adata.X may contain raw counts (max > 100).

            If raw counts, please select "adata.X \u2192 CPM normalize" instead.
            """)
        else:
            st.info("\u2713 Using adata.X as-is (CPM normalization will not be applied)")
    elif "CPM normalize" in selected_layer:
        _arrow = "\u2192"
        _layer_name = selected_layer.split(f' {_arrow}')[0]
        st.success(f"""
        \u2705 **{_layer_name}** {_arrow} **CPM normalization**

        CPM normalization (target_sum=1,000,000, no log transformation) will be applied automatically before scFEA execution.
        """)

    # Calculate sparsity for selected layer (for MAGIC recommendation)
    if "use as-is" in selected_layer:
        sel_sparsity = x_sparsity
    elif "adata.raw.X" in selected_layer:
        sel_data = adata.raw.X.toarray() if hasattr(adata.raw.X, 'toarray') else adata.raw.X
        sel_sparsity = float((sel_data == 0).sum() / sel_data.size * 100)
    elif "adata.layers['" in selected_layer:
        layer_name = selected_layer.split("'")[1]
        sel_data = adata.layers[layer_name]
        if hasattr(sel_data, 'toarray'):
            sel_data = sel_data.toarray()
        sel_sparsity = float((sel_data == 0).sum() / sel_data.size * 100)
    else:
        sel_sparsity = x_sparsity

    # MAGIC recommendation based on sparsity
    if sel_sparsity > 70:
        st.info(f"\U0001f4a1 Sparsity {sel_sparsity:.1f}% (>70%): MAGIC imputation is recommended (e.g., 10x data)")
        magic_recommended = True
    elif sel_sparsity > 50:
        st.info(f"\u2139\ufe0f Sparsity {sel_sparsity:.1f}% (50-70%): MAGIC imputation is optional")
        magic_recommended = False
    else:
        st.info(f"\u2139\ufe0f Sparsity {sel_sparsity:.1f}% (<50%): MAGIC imputation is not needed (e.g., full-length data)")
        magic_recommended = False

    st.session_state.scfea_selected_layer = selected_layer
    st.session_state.scfea_sparsity = sel_sparsity
    st.session_state.scfea_magic_recommended = magic_recommended

    # Clean up temporary adata
    del adata

    # ========================================
    # Step 2: Configure analysis parameters
    # ========================================
    st.header("Step 2: Configure analysis")

    with st.expander("📚 Parameter Guide", expanded=False):
        st.markdown("""
        ### Metabolic Model
        | Model | Modules | Purpose |
        |-------|---------|------|
        | **Human Complete** | 168 | Comprehensive metabolic analysis |
        | **Mouse Complete** | 168 | For mouse data |
        | **Glutaminolysis** | 20-23 | Glutamine metabolism focused (fast) |
        | **Iron Ion** | 15 | Iron metabolism focused (fast) |

        ### MAGIC Imputation
        - Recommended for data with high dropout rates such as 10x
        - Completes missing values in sparse data
        - Increases computation time, so may be unnecessary for small datasets

        ### Training Epochs
        - Default: 100
        - Increase if not converging (200-500)
        - Decrease for faster computation (50)

        ### Output
        - **Flux matrix**: Flux values for each metabolic module per cell
        - **Balance matrix**: Metabolite stress values per cell
        """)

    with st.form("scfea_params_form"):
        col1, col2 = st.columns(2)

        with col1:
            # Filter models by detected species
            detected_species = "human" if info['detected_species_idx'] == 1 else "mouse"
            available_models = [
                name for name, config in MODEL_CONFIGS.items()
                if config["species"] == detected_species or config["species"] == "human"
            ]

            # Set default to the matching Complete model
            if detected_species == "mouse":
                default_model = "Mouse Complete (168 modules)"
            else:
                default_model = "Human Complete (168 modules)"
            default_idx = available_models.index(default_model) if default_model in available_models else 0

            selected_model = st.selectbox(
                "Metabolic Model",
                options=available_models,
                index=default_idx,
                help="Default model is set based on species auto-detected from gene names"
            )

            # Show model description
            model_config = MODEL_CONFIGS[selected_model]
            st.caption(model_config["description"])

        with col2:
            train_epochs = st.number_input(
                "Training Epochs",
                min_value=10,
                max_value=1000,
                value=100,
                step=10,
                help="Number of training iterations for the neural network"
            )

        col1, col2 = st.columns(2)

        with col1:
            # Get MAGIC recommendation from session state
            magic_recommended = st.session_state.get("scfea_magic_recommended", False)
            sparsity = st.session_state.get("scfea_sparsity", 0)

            # Build help text with sparsity info
            if magic_installed:
                if sparsity > 70:
                    help_text = f"Sparsity {sparsity:.1f}% -> **Recommended** (e.g., 10x data)"
                elif sparsity > 50:
                    help_text = f"Sparsity {sparsity:.1f}% -> Optional"
                else:
                    help_text = f"Sparsity {sparsity:.1f}% -> Not needed (e.g., full-length data)"
            else:
                help_text = "MAGIC is not installed"

            use_imputation = st.checkbox(
                "Enable MAGIC Imputation",
                value=magic_recommended and magic_installed,
                disabled=not magic_installed,
                help=help_text
            )

        with col2:
            output_prefix = st.text_input(
                "Output prefix",
                value="scfea_result",
                help="Prefix for output files"
            )

        # Estimate time
        st.markdown("---")
        n_cells = info['n_cells']
        n_modules = int(selected_model.split("(")[1].split(" ")[0]) if "(" in selected_model else 168

        if cuda_available:
            estimated_time = "a few minutes" if n_cells < 5000 else "tens of minutes"
        else:
            estimated_time = "a few to tens of minutes" if n_cells < 5000 else "tens of minutes to 1 hour"

        st.info(f"""
        **Estimated computation time**: {estimated_time}
        - Cells: {n_cells:,}
        - Modules: {n_modules}
        - Device: {"GPU" if cuda_available else "CPU"}
        - Epochs: {train_epochs}
        """)

        submitted = st.form_submit_button("🚀 Run scFEA", type="primary")

    # ========================================
    # Step 3: Run analysis
    # ========================================
    if submitted:
        st.header("Step 3: Running scFEA")

        # Prepare input data
        with st.spinner("Preparing input data..."):
            # Load adata
            adata = sc.read_h5ad(st.session_state.scfea_h5ad_path)

            # Get selected layer
            selected_layer = st.session_state.get("scfea_selected_layer", "adata.X \u2192 CPM normalize (if raw counts)")

            st.info(f"\U0001f4ca **Data source:** {selected_layer}")

            # Check if CPM normalization is needed
            apply_cpm = "CPM normalize" in selected_layer

            # Apply selected layer to adata.X
            if "adata.raw.X" in selected_layer:
                if hasattr(adata, 'raw') and adata.raw is not None:
                    adata = adata.raw.to_adata()
                    st.info("\u2713 Copied adata.raw.X to .X")
                else:
                    st.error("adata.raw.X is not available")
                    st.stop()
            elif "adata.layers['" in selected_layer:
                layer_name = selected_layer.split("'")[1]
                if layer_name in adata.layers:
                    adata.X = adata.layers[layer_name].copy()
                    st.info(f"\u2713 Copied adata.layers['{layer_name}'] to .X")
                else:
                    st.error(f"Layer '{layer_name}' not found")
                    st.stop()
            elif "use as-is" in selected_layer:
                st.info("\u2713 Using adata.X as-is (no CPM normalization)")
            # else: "adata.X -> CPM normalize" - use adata.X, apply CPM below

            # Apply CPM normalization if selected
            if apply_cpm:
                import matplotlib.pyplot as plt

                # Store total counts before normalization
                if hasattr(adata.X, 'toarray'):
                    total_counts_before = np.array(adata.X.sum(axis=1)).flatten()
                else:
                    total_counts_before = np.array(adata.X.sum(axis=1)).flatten()

                # Apply CPM normalization (no log transform)
                sc.pp.normalize_total(adata, target_sum=1e6)

                st.success("""
                \u2705 **CPM normalization applied**

                - target_sum = 1,000,000 (CPM)
                - No log transformation (scFEA internally applies log2 transformation when max > 50)
                """)

                # Show distribution plot
                with st.expander("\U0001f4ca Data distribution (before/after normalization)", expanded=False):
                    fig, axes = plt.subplots(1, 2, figsize=(10, 4))

                    # Plot 1: Total counts before normalization
                    axes[0].hist(total_counts_before, bins=30, alpha=0.7, color='steelblue')
                    axes[0].set_xlabel('Total counts')
                    axes[0].set_ylabel('Number of cells')
                    axes[0].set_title('Before CPM normalization')
                    axes[0].ticklabel_format(style='scientific', axis='x', scilimits=(0, 0))

                    # Plot 2: Expression distribution after normalization
                    if hasattr(adata.X, 'toarray'):
                        X_sample = adata.X[:, :1000].toarray().flatten()
                    else:
                        X_sample = adata.X[:, :1000].flatten()
                    X_nonzero = X_sample[X_sample > 0]
                    axes[1].hist(np.log10(X_nonzero + 1), bins=50, alpha=0.7, color='forestgreen')
                    axes[1].set_xlabel('log10(CPM + 1)')
                    axes[1].set_ylabel('Frequency')
                    axes[1].set_title('After CPM normalization')

                    plt.tight_layout()
                    st.pyplot(fig)
                    plt.close()

            data_matrix = adata.X

            # Get expression matrix (genes x cells)
            # scFEA expects genes as rows, cells as columns
            if hasattr(data_matrix, 'toarray'):
                expr_matrix = pd.DataFrame(
                    data_matrix.toarray().T,
                    index=adata.var_names,
                    columns=adata.obs_names
                )
            else:
                expr_matrix = pd.DataFrame(
                    data_matrix.T,
                    index=adata.var_names,
                    columns=adata.obs_names
                )

            # Save as CSV
            input_csv_path = os.path.join(scfea_temp_dir, "input_expression.csv")
            expr_matrix.to_csv(input_csv_path)
            st.success(f"\u2713 Expression matrix saved: {expr_matrix.shape[0]} genes x {expr_matrix.shape[1]} cells")

        # Create output directory
        output_dir = os.path.join(scfea_temp_dir, "output")
        os.makedirs(output_dir, exist_ok=True)

        # Prepare output file names
        flux_file = os.path.join(output_dir, f"{output_prefix}_flux.csv")
        balance_file = os.path.join(output_dir, f"{output_prefix}_balance.csv")

        # Build command
        model_config = MODEL_CONFIGS[selected_model]

        # Run scFEA using Python directly (not subprocess)
        progress_bar = st.progress(0)
        status_text = st.empty()

        try:
            # Add scFEA src to path
            if str(SCFEA_SRC) not in sys.path:
                sys.path.insert(0, str(SCFEA_SRC))

            status_text.text("Loading scFEA modules...")
            progress_bar.progress(10)

            # Import required modules
            import torch
            from torch.autograd import Variable
            import warnings

            # Import scFEA components
            from ClassFlux import FLUX
            from util import pearsonr
            from DatasetFlux import MyDataset

            # scFEA hyperparameters
            LEARN_RATE = 0.008
            LAMB_BA = 1
            LAMB_NG = 1
            LAMB_CELL = 1
            LAMB_MOD = 1e-2

            status_text.text("Loading data...")
            progress_bar.progress(20)

            # Device selection
            device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
            st.info(f"Using device: {device}")

            # Read expression data
            geneExpr = pd.read_csv(input_csv_path, index_col=0)
            geneExpr = geneExpr.T  # cells x genes
            geneExpr = geneExpr * 1.0

            # MAGIC imputation
            if use_imputation and magic_installed:
                status_text.text("Applying MAGIC imputation...")
                progress_bar.progress(25)
                import magic as magic_pkg
                magic_operator = magic_pkg.MAGIC()
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    geneExpr = magic_operator.fit_transform(geneExpr)

            # Log transform if needed
            if geneExpr.max().max() > 50:
                geneExpr = (geneExpr + 1).apply(np.log2)

            geneExprSum = geneExpr.sum(axis=1)
            stand = geneExprSum.mean()
            geneExprScale = geneExprSum / stand
            geneExprScale = torch.FloatTensor(geneExprScale.values).to(device)

            BATCH_SIZE = geneExpr.shape[0]

            status_text.text("Loading metabolic model...")
            progress_bar.progress(30)

            # Load module genes
            moduleGene = pd.read_csv(
                SCFEA_DATA / model_config["moduleGene"],
                sep=',',
                index_col=0
            )
            moduleLen = np.array([moduleGene.iloc[i, :].notna().sum() for i in range(moduleGene.shape[0])])

            # Find overlapping genes
            module_gene_all = set()
            for i in range(moduleGene.shape[0]):
                for j in range(moduleGene.shape[1]):
                    if pd.notna(moduleGene.iloc[i, j]):
                        module_gene_all.add(moduleGene.iloc[i, j])

            data_gene_all = set(geneExpr.columns)
            gene_overlap = sorted(list(data_gene_all.intersection(module_gene_all)))

            n_overlap = len(gene_overlap)
            n_total_module_genes = len(module_gene_all)
            overlap_pct = n_overlap / n_total_module_genes * 100

            st.info(f"Gene overlap: {n_overlap}/{n_total_module_genes} ({overlap_pct:.1f}%) metabolic genes found in data")

            if overlap_pct < 30:
                st.warning("⚠️ Low gene overlap. Results may be unreliable. Check if gene symbols match the expected format.")

            # Load stoichiometry matrix
            cmMat = pd.read_csv(SCFEA_DATA / model_config["stoichiometry"], sep=',', header=None)
            cmMat = torch.FloatTensor(cmMat.values).to(device)

            # Load compound names
            cName = pd.read_csv(SCFEA_DATA / model_config["cName"], sep=',', header=0)
            cName = cName.columns

            status_text.text("Processing data...")
            progress_bar.progress(40)

            # Process expression data
            geneExpr = geneExpr[gene_overlap]
            gene_names = geneExpr.columns
            cell_names = geneExpr.index.astype(str)
            n_modules = moduleGene.shape[0]
            n_genes = len(gene_names)
            n_cells = len(cell_names)
            n_comps = cmMat.shape[0]

            # Build module-gene expression matrix
            geneExprDf_list = []
            emptyNode = []

            for i in range(n_modules):
                genes = moduleGene.iloc[i, :].values.astype(str)
                genes = [g for g in genes if g != 'nan']
                if not genes:
                    emptyNode.append(i)
                    continue
                temp = geneExpr.copy()
                temp.loc[:, [g for g in gene_names if g not in genes]] = 0
                temp = temp.T
                temp['Module_Gene'] = ['%02d_%s' % (i, g) for g in gene_names]
                geneExprDf_list.append(temp)

            geneExprDf = pd.concat(geneExprDf_list, ignore_index=True)
            geneExprDf.index = geneExprDf['Module_Gene']
            geneExprDf.drop('Module_Gene', axis='columns', inplace=True)
            X = torch.FloatTensor(geneExprDf.values.T).to(device)

            # Prepare module scale
            df = geneExprDf.copy()
            df.index = [int(i.split('_')[0]) for i in df.index]
            module_scale = df.groupby(df.index).sum().T
            module_scale = torch.FloatTensor(module_scale.values / moduleLen).to(device)

            status_text.text("Initializing neural network...")
            progress_bar.progress(50)

            # Initialize network
            torch.manual_seed(16)
            net = FLUX(X, n_modules, f_in=n_genes, f_out=1).to(device)
            optimizer = torch.optim.Adam(net.parameters(), lr=LEARN_RATE)

            # Dataloader
            dataloader_params = {
                'batch_size': BATCH_SIZE,
                'shuffle': False,
                'num_workers': 0,
                'pin_memory': False
            }

            dataSet = MyDataset(X, geneExprScale, module_scale)
            train_loader = torch.utils.data.DataLoader(dataset=dataSet, **dataloader_params)

            # Training
            status_text.text("Training neural network...")

            # Loss function (inline to avoid import issues)
            def myLoss(m, c, lamb1=0.2, lamb2=0.2, lamb3=0.2, lamb4=0.2, geneScale=None, moduleScale=None):
                total1 = torch.sum(torch.pow(c, 2), dim=1)
                error = torch.abs(m) - m
                total2 = torch.sum(error, dim=1)
                diff = torch.pow(torch.sum(m, dim=1) - geneScale, 2)
                if (diff > 0).sum().item() == m.shape[0]:
                    total3 = torch.pow(diff, 0.5)
                else:
                    total3 = diff

                if lamb4 > 0:
                    corr = torch.ones(m.shape[0], device=m.device)
                    for i in range(m.shape[0]):
                        corr[i] = pearsonr(m[i, :], moduleScale[i, :])
                    corr = torch.abs(corr)
                    total4 = torch.ones(m.shape[0], device=m.device) - corr
                else:
                    total4 = torch.zeros(m.shape[0], device=m.device)

                loss = torch.sum(lamb1 * total1 + lamb2 * total2 + lamb3 * total3 + lamb4 * total4)
                return loss, torch.sum(lamb1 * total1), torch.sum(lamb2 * total2), torch.sum(lamb3 * total3), torch.sum(lamb4 * total4)

            loss_history = []
            net.train()
            start_time = time.time()

            for epoch in range(train_epochs):
                loss = 0
                for i, (X_batch, X_scale, m_scale) in enumerate(train_loader):
                    X_batch = Variable(X_batch.float().to(device))
                    X_scale_batch = Variable(X_scale.float().to(device))
                    m_scale_batch = Variable(m_scale.float().to(device))

                    out_m_batch, out_c_batch = net(X_batch, n_modules, n_genes, n_comps, cmMat)
                    loss_batch, _, _, _, _ = myLoss(
                        out_m_batch, out_c_batch,
                        lamb1=LAMB_BA, lamb2=LAMB_NG, lamb3=LAMB_CELL, lamb4=LAMB_MOD,
                        geneScale=X_scale_batch, moduleScale=m_scale_batch
                    )

                    optimizer.zero_grad()
                    loss_batch.backward()
                    optimizer.step()

                    loss += loss_batch.cpu().data.numpy()

                loss_history.append(loss)

                # Update progress
                progress = 50 + int((epoch + 1) / train_epochs * 40)
                progress_bar.progress(progress)
                status_text.text(f"Training: Epoch {epoch + 1}/{train_epochs}, Loss: {loss:.4f}")

            training_time = time.time() - start_time
            st.success(f"✓ Training completed in {training_time:.1f} seconds")

            status_text.text("Generating predictions...")
            progress_bar.progress(95)

            # Prediction
            dataloader_params['batch_size'] = 1
            dataSet = MyDataset(X, geneExprScale, module_scale)
            test_loader = torch.utils.data.DataLoader(dataset=dataSet, **dataloader_params)

            fluxStatuTest = np.zeros((n_cells, n_modules), dtype='f')
            balanceStatus = np.zeros((n_cells, n_comps), dtype='f')

            net.eval()
            with torch.no_grad():
                for i, (X_batch, X_scale, _) in enumerate(test_loader):
                    X_batch = Variable(X_batch.float().to(device))
                    out_m_batch, out_c_batch = net(X_batch, n_modules, n_genes, n_comps, cmMat)
                    fluxStatuTest[i, :] = out_m_batch.cpu().numpy()
                    balanceStatus[i, :] = out_c_batch.cpu().numpy()

            # Save results
            status_text.text("Saving results...")

            # Load module information for descriptive names
            module_info_path = SCFEA_DATA / "Human_M168_information.symbols.csv"
            if module_info_path.exists():
                module_info = pd.read_csv(module_info_path)
                # Create descriptive names: M_1_Glucose_to_G6P
                module_name_map = {}
                for _, row in module_info.iterrows():
                    module_id = row['Unnamed: 0']
                    in_name = str(row['Compound_IN_name']).replace(' ', '_').replace('+', '_')
                    out_name = str(row['Compound_OUT_name']).replace(' ', '_').replace('+', '_')
                    module_name_map[module_id] = f"{module_id}_{in_name}_to_{out_name}"
            else:
                module_name_map = None

            # Flux matrix
            flux_df = pd.DataFrame(fluxStatuTest)
            flux_df.columns = moduleGene.index
            # Apply descriptive names if available
            if module_name_map:
                flux_df.columns = [module_name_map.get(col, col) for col in flux_df.columns]
            flux_df.index = geneExpr.index.tolist()
            flux_df.to_csv(flux_file)

            # Balance matrix
            balance_df = pd.DataFrame(balanceStatus)
            balance_df.index = flux_df.index
            balance_df.columns = cName
            balance_df.to_csv(balance_file)

            progress_bar.progress(100)
            status_text.text("Complete!")

            st.session_state.scfea_complete = True
            st.session_state.scfea_flux_file = flux_file
            st.session_state.scfea_balance_file = balance_file
            st.session_state.scfea_loss_history = loss_history
            st.session_state.scfea_flux_df = flux_df
            st.session_state.scfea_balance_df = balance_df

            st.success("✅ scFEA analysis completed!")

        except Exception as e:
            st.error(f"Error running scFEA: {str(e)}")
            import traceback
            st.code(traceback.format_exc())
            st.stop()

    # ========================================
    # Step 4: Results
    # ========================================
    if st.session_state.get("scfea_complete", False):
        st.header("Step 4: Results")

        flux_df = st.session_state.scfea_flux_df
        balance_df = st.session_state.scfea_balance_df
        loss_history = st.session_state.scfea_loss_history

        # Summary metrics
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Cells", flux_df.shape[0])
        with col2:
            st.metric("Flux Modules", flux_df.shape[1])
        with col3:
            st.metric("Metabolites", balance_df.shape[1])

        # Loss plot
        st.subheader("Training Loss")
        import plotly.express as px
        fig = px.line(
            x=list(range(1, len(loss_history) + 1)),
            y=loss_history,
            labels={'x': 'Epoch', 'y': 'Loss'},
            title='Training Loss Curve'
        )
        st.plotly_chart(fig, use_container_width=True)

        # Flux matrix preview
        st.subheader("Metabolic Flux Matrix")
        st.dataframe(flux_df.head(20), use_container_width=True)

        # Flux heatmap
        st.subheader("Flux Heatmap (Top 50 cells)")
        import plotly.express as px
        fig = px.imshow(
            flux_df.head(50).T,
            labels=dict(x="Cells", y="Modules", color="Flux"),
            title="Metabolic Flux",
            aspect="auto"
        )
        st.plotly_chart(fig, use_container_width=True)

        # Balance matrix preview
        st.subheader("Metabolite Balance Matrix")
        st.dataframe(balance_df.head(20), use_container_width=True)

        # Download buttons
        st.subheader("Download Results")

        # Create cell_metadata from original h5ad
        adata_orig = sc.read_h5ad(st.session_state.scfea_h5ad_path)
        cell_metadata = adata_orig.obs.copy()
        cell_metadata.index.name = 'cell_id'
        # Align to flux_df cells
        common_cells = [c for c in flux_df.index if c in cell_metadata.index]
        cell_metadata = cell_metadata.loc[common_cells]
        cell_metadata_csv = cell_metadata.to_csv()
        del adata_orig

        import zipfile
        from io import BytesIO

        cell_metadata_tsv = cell_metadata.to_csv(sep='\t')

        zip_buffer = BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
            zf.writestr(f"{output_prefix}_flux.tsv", flux_df.to_csv(sep='\t'))
            zf.writestr(f"{output_prefix}_balance.tsv", balance_df.to_csv(sep='\t'))
            zf.writestr(f"{output_prefix}_cell_metadata.tsv", cell_metadata_tsv)
        zip_buffer.seek(0)

        st.download_button(
            label="📥 Download Results (ZIP)",
            data=zip_buffer.getvalue(),
            file_name=f"{output_prefix}_results.zip",
            mime="application/zip"
        )

        st.markdown(f"""
        **ZIP contains:**
        - `{output_prefix}_flux.tsv` - Metabolic flux matrix (cells x modules)
        - `{output_prefix}_balance.tsv` - Metabolite balance matrix (cells x metabolites)
        - `{output_prefix}_cell_metadata.tsv` - Cell metadata (for scFEA postprocess)
        """)

        # Export as h5ad files
        st.subheader("Export as AnnData (h5ad)")

        st.markdown("""
        scFEA output (Flux/Balance) contains relative values comparable between cells.
        Library size normalization (TPM/CPM conversion, etc.) is not required.

        Below, two h5ad files are exported with Flux/Balance values stored in `adata.X`.
        Following the official tutorial (scflux.org), analysis proceeds in order: scaling -> PCA -> clustering.
        """)

        if st.button("Generate h5ad files"):
            with st.spinner("Creating h5ad files..."):
                # Get original cell metadata
                adata_orig = sc.read_h5ad(st.session_state.scfea_h5ad_path)

                # Load module info for supermodule annotations
                import anndata as ad
                module_info_path = SCFEA_DATA / "Human_M168_information.symbols.csv"
                if module_info_path.exists():
                    module_info = pd.read_csv(module_info_path)
                    module_info.set_index('Unnamed: 0', inplace=True)
                else:
                    module_info = None

                # Supermodule name mapping
                supermodule_names = {
                    1: "Glycolysis_TCA", 2: "Serine", 3: "Pentose_phosphate",
                    4: "Fatty_acids", 5: "Aspartate", 6: "Beta_alanine",
                    7: "Propionyl_CoA", 8: "Glutamate", 9: "BCAA",
                    10: "Urea_cycle", 11: "Spermine", 12: "Transporters",
                    13: "Hyaluronic_acid", 14: "Glycogen", 15: "Glycosaminoglycan",
                    16: "N_glycan", 17: "O_glycan", 18: "Sialic_acid",
                    19: "Glycan", 20: "Purine", 21: "Pyrimidine", 22: "Steroid_hormone"
                }

                # Create Flux AnnData
                adata_flux = ad.AnnData(
                    X=flux_df.values,
                    obs=adata_orig.obs.copy(),
                    var=pd.DataFrame(index=flux_df.columns)
                )
                adata_flux.var_names = flux_df.columns
                adata_flux.obs_names = flux_df.index

                # Add supermodule info to var
                if module_info is not None:
                    # Extract original module ID (M_1, M_2, etc.) from descriptive name
                    original_ids = [col.split('_')[0] + '_' + col.split('_')[1] for col in flux_df.columns]
                    supermodule_ids = [module_info.loc[mid, 'Supermodule_id'] if mid in module_info.index else None for mid in original_ids]
                    supermodule_labels = [supermodule_names.get(sid, '') if sid else '' for sid in supermodule_ids]
                    adata_flux.var['supermodule_id'] = supermodule_ids
                    adata_flux.var['supermodule_name'] = supermodule_labels

                # Copy original embeddings (PCA, UMAP, etc.)
                if adata_orig.obsm:
                    for key in adata_orig.obsm.keys():
                        adata_flux.obsm[key] = adata_orig.obsm[key].copy()
                # Mark all modules as highly variable for PCA
                adata_flux.var['highly_variable'] = True
                adata_flux.uns['data_type'] = 'scFEA_flux'
                adata_flux.uns['description'] = 'Metabolic flux values from scFEA. All modules marked as highly_variable.'

                # Create Balance AnnData
                adata_balance = ad.AnnData(
                    X=balance_df.values,
                    obs=adata_orig.obs.copy(),
                    var=pd.DataFrame(index=balance_df.columns)
                )
                adata_balance.var_names = balance_df.columns
                adata_balance.obs_names = balance_df.index
                # Copy original embeddings (PCA, UMAP, etc.)
                if adata_orig.obsm:
                    for key in adata_orig.obsm.keys():
                        adata_balance.obsm[key] = adata_orig.obsm[key].copy()
                # Mark all metabolites as highly variable for PCA
                adata_balance.var['highly_variable'] = True
                adata_balance.uns['data_type'] = 'scFEA_balance'
                adata_balance.uns['description'] = 'Metabolite balance values from scFEA. All metabolites marked as highly_variable.'

                # Save files
                flux_h5ad_path = os.path.join(scfea_temp_dir, f"{output_prefix}_flux.h5ad")
                balance_h5ad_path = os.path.join(scfea_temp_dir, f"{output_prefix}_balance.h5ad")

                adata_flux.write_h5ad(flux_h5ad_path)
                adata_balance.write_h5ad(balance_h5ad_path)

                st.success("✓ h5ad files created")

                st.info("""
                **Analysis method in SCALA (common for Flux/Balance)**:
                1. Load the h5ad file
                2. Select "Scale all genes"
                3. Check "Use all scaled data (not just variable features)"
                4. PCA -> UMAP -> Clustering

                * Flux: 168 modules, Balance: 70 metabolites (for the Complete model)
                * FindVariableFeatures is not needed (all features are used)
                """)

                st.info("""
                **Analysis method with Scanpy/Python**:
                ```python
                import scanpy as sc
                adata = sc.read_h5ad('xxx_flux.h5ad')  # or xxx_balance.h5ad

                # Scale all features -> PCA -> clustering
                sc.pp.scale(adata)
                sc.pp.pca(adata, n_comps=min(10, adata.n_vars-1))
                sc.pp.neighbors(adata)
                sc.tl.umap(adata)
                sc.tl.leiden(adata)
                ```
                """)

                # Create zip file with h5ad files and cell metadata
                import zipfile
                zip_path = os.path.join(scfea_temp_dir, f"{output_prefix}_scfea.zip")
                with zipfile.ZipFile(zip_path, 'w', zipfile.ZIP_DEFLATED) as zf:
                    zf.write(flux_h5ad_path, f"{output_prefix}_flux.h5ad")
                    zf.write(balance_h5ad_path, f"{output_prefix}_balance.h5ad")
                    zf.writestr(f"{output_prefix}_cell_metadata.tsv", cell_metadata_tsv)

                # Download button
                with open(zip_path, 'rb') as f:
                    st.download_button(
                        label="📥 Download h5ad files (ZIP)",
                        data=f.read(),
                        file_name=f"{output_prefix}_scfea.zip",
                        mime="application/zip",
                        help="Contains flux.h5ad and balance.h5ad"
                    )

else:
    st.info("👆 Upload an h5ad file to start analysis")
