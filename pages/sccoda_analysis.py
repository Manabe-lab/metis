"""
scCODA - Compositional Data Analysis for Single-Cell
Perform differential composition analysis using scCODA (pertpy)
"""

import streamlit as st
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import time
import zipfile
import io
from matplotlib.backends.backend_pdf import PdfPages
from helper_func import clear_old_directories, clear_old_files

# Import pertpy for scCODA
try:
    import pertpy as pt
    PERTPY_AVAILABLE = True
except ImportError:
    PERTPY_AVAILABLE = False

st.set_page_config(page_title="scCODA Analysis", page_icon="📊", layout="wide")

st.title("📊 scCODA - Compositional Data Analysis")

if not PERTPY_AVAILABLE:
    st.error("""
    ❌ **pertpy is not installed**

    Please install pertpy using:
    ```bash
    pip install pertpy
    ```

    See: https://pertpy.readthedocs.io/
    """)
    st.stop()

st.markdown("""
Statistically detect **cell type composition changes** in single-cell data using scCODA.

### scCODA Features
- **Compositional Data support**: Properly handles compositional constraints of ratio data
- **Hierarchical Bayesian model**: Estimates cell type proportion changes between groups
- **Sparse estimation**: Detects only cell types that actually change
- **Reference cell type**: Can select cell type assumed to be unchanged
- **Bayesian statistics**: Significance determined by Inclusion probability, not p-value or FDR

### Workflow
1. **h5ad file upload**: Data containing cell type, sample, and condition information
2. **Parameter setting**: Select cell type column, sample column, condition column, reference cell type
3. **Run scCODA**: Differential composition analysis by Bayesian estimation
4. **Display results**: Visualization and download of significantly changed cell types

### References
- [Büttner et al. (2021) "scCODA is a Bayesian model for compositional single-cell data analysis" Nature Communications](https://www.nature.com/articles/s41467-021-27150-6)
- [pertpy Documentation - scCODA tutorial](https://pertpy.readthedocs.io/en/latest/tutorials/notebooks/sccoda_extended.html)
- [scCODA GitHub](https://github.com/theislab/scCODA)
""")

# Initialize session state
if "sccoda_temp_dir" not in st.session_state:
    sccoda_temp_dir = os.path.join("temp", f"sccoda_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(sccoda_temp_dir, exist_ok=True)
    st.session_state.sccoda_temp_dir = sccoda_temp_dir
else:
    sccoda_temp_dir = st.session_state.sccoda_temp_dir

if "sccoda_adata" not in st.session_state:
    st.session_state.sccoda_adata = None

if "sccoda_results" not in st.session_state:
    st.session_state.sccoda_results = None

if "sccoda_plots" not in st.session_state:
    st.session_state.sccoda_plots = {
        "step2": [],  # Step 2 plots
        "step3": []   # Step 3 plots
    }

# ========================================
# Step 1: Upload h5ad file
# ========================================
st.header("Step 1: Upload h5ad file")

st.markdown("""
### Required Information
h5ad file (`adata`) requires the following information:

- **Cell type information**: Included as column in `adata.obs` (e.g., `'cell_type'`, `'clusters'`, `'leiden'`)
- **Sample information**: Included as column in `adata.obs` (e.g., `'sample'`, `'batch'`, `'replicate'`)
- **Condition information**: Included as column in `adata.obs` (e.g., `'condition'`, `'group'`, `'treatment'`)

💡 **Tip**: Data analyzed with Seurat should be converted to `anndata` format in Python before uploading.
""")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file",
    type=["h5ad"],
    help="h5ad file containing cell type, sample, and condition information"
)

@st.cache_data(show_spinner=False)
def load_h5ad_cached(file_content, file_name):
    """Load h5ad file with caching based on file content hash"""
    import hashlib
    import tempfile

    # Create temp file
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp.write(file_content)
        tmp_path = tmp.name

    # Load anndata
    adata = sc.read_h5ad(tmp_path)

    # Clean up
    os.unlink(tmp_path)

    return adata

if uploaded_h5ad is not None:
    # Get file content for caching
    file_content = uploaded_h5ad.read()
    uploaded_h5ad.seek(0)  # Reset file pointer

    # Check if we need to reload (new file)
    import hashlib
    file_hash = hashlib.md5(file_content).hexdigest()

    if "sccoda_file_hash" not in st.session_state or st.session_state.sccoda_file_hash != file_hash:
        with st.spinner("Loading h5ad file..."):
            try:
                adata = load_h5ad_cached(file_content, uploaded_h5ad.name)
                st.session_state.sccoda_adata = adata
                st.session_state.sccoda_file_hash = file_hash
                # Reset parameters when new file is loaded
                if "sccoda_params" in st.session_state:
                    del st.session_state.sccoda_params
                # Reset results
                st.session_state.sccoda_results = None

                st.success(f"✅ h5ad file loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
            except Exception as e:
                st.error(f"❌ h5ad file loading error: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
                st.stop()
    else:
        adata = st.session_state.sccoda_adata
        st.success(f"✅ h5ad file (cached): {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Display data summary
    adata = st.session_state.sccoda_adata

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Total cells", f"{adata.shape[0]:,}")
    with col2:
        st.metric("Genes", f"{adata.shape[1]:,}")
    with col3:
        if hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
            st.metric("Metadata columns", len(adata.obs.columns))

    # Show available metadata columns
    with st.expander("📋 Available metadata columns"):
        st.write("**Columns in adata.obs:**")
        st.dataframe(pd.DataFrame({
            "Column name": adata.obs.columns.tolist(),
            "Data type": [str(dtype) for dtype in adata.obs.dtypes],
            "Unique values": [adata.obs[col].nunique() for col in adata.obs.columns]
        }))

    # Preview data
    with st.expander("🔍 Data preview (first 10 rows)"):
        st.dataframe(adata.obs.head(10))

# ========================================
# Step 2: Configure scCODA parameters
# ========================================
if st.session_state.sccoda_adata is not None:
    st.header("Step 2: Configure scCODA parameters")

    adata = st.session_state.sccoda_adata
    obs_columns = adata.obs.columns.tolist()

    # Filtering: Only columns with <50 unique values are candidates
    filtered_columns = [col for col in obs_columns if adata.obs[col].nunique() < 50]

    if len(filtered_columns) == 0:
        st.error("❌ No columns with <50 unique values found. Please check metadata columns.")
        st.stop()

    # Default column selection function
    def find_default_celltype_col(columns, obs_df):
        """Find default cell type column"""
        # 1. Columns containing ident (excluding orig.ident)
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 2. Columns containing type
        type_cols = [col for col in columns if 'type' in col.lower()]
        if type_cols:
            return type_cols[0]
        # 3. Default is first column
        return columns[0]

    def find_default_sample_col(columns, obs_df):
        """Find default sample column"""
        # 1. Prioritize orig.ident
        if 'orig.ident' in columns:
            return 'orig.ident'
        # 2. Columns containing ident (other than orig.ident)
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 3. Columns containing sample
        sample_cols = [col for col in columns if 'sample' in col.lower()]
        if sample_cols:
            return sample_cols[0]
        # 4. Default is first column
        return columns[0]

    def find_default_condition_col(columns, obs_df):
        """Find default condition column"""
        # Find columns containing condition, stim, KO
        keywords = ['condition', 'stim', 'ko']
        for keyword in keywords:
            matching_cols = [col for col in columns if keyword in col.lower()]
            if matching_cols:
                return matching_cols[0]
        # Default is second column (first is likely used for sample)
        return columns[min(1, len(columns)-1)]

    # Determine default columns
    default_celltype = find_default_celltype_col(filtered_columns, adata.obs)
    default_sample = find_default_sample_col(filtered_columns, adata.obs)
    default_condition = find_default_condition_col(filtered_columns, adata.obs)

    # Get index
    celltype_idx = filtered_columns.index(default_celltype) if default_celltype in filtered_columns else 0
    sample_idx = filtered_columns.index(default_sample) if default_sample in filtered_columns else 0
    condition_idx = filtered_columns.index(default_condition) if default_condition in filtered_columns else min(1, len(filtered_columns)-1)

    # Initialize session state for form values
    if "sccoda_params" not in st.session_state:
        st.session_state.sccoda_params = {
            "celltype_col": default_celltype,
            "sample_col": default_sample,
            "condition_col": default_condition,
            "reference_celltype": "Auto select",
            "inclusion_threshold": 0.95
        }

    # Get current values from session state for form defaults
    current_celltype = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    current_sample = st.session_state.sccoda_params.get("sample_col", default_sample)
    current_condition = st.session_state.sccoda_params.get("condition_col", default_condition)
    current_reference = st.session_state.sccoda_params.get("reference_celltype", "Auto select")
    current_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Calculate indices for selectboxes based on current values
    form_celltype_idx = filtered_columns.index(current_celltype) if current_celltype in filtered_columns else celltype_idx
    form_sample_idx = filtered_columns.index(current_sample) if current_sample in filtered_columns else sample_idx
    form_condition_idx = filtered_columns.index(current_condition) if current_condition in filtered_columns else condition_idx

    # Form for parameter selection
    st.info("💡 After selecting parameters, click "✅ Confirm Parameters" button")
    with st.form("sccoda_params_form"):
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("Data Column Selection")

            # Cell type column
            celltype_col = st.selectbox(
                "Cell type column",
                filtered_columns,
                index=form_celltype_idx,
                help="Column containing cell type or cluster information"
            )

            # Sample column
            sample_col = st.selectbox(
                "Sample column",
                filtered_columns,
                index=form_sample_idx,
                help="Column containing biological replicate (sample) information"
            )

            # Condition column
            condition_col = st.selectbox(
                "Condition column (groups to compare)",
                filtered_columns,
                index=form_condition_idx,
                help="Conditions to compare (e.g., Control vs Treatment)"
            )

        with col2:
            st.subheader("Model Parameters")

            # Reference cell type - use celltypes from current celltype column
            current_celltype_for_ref = filtered_columns[form_celltype_idx]
            all_celltypes = adata.obs[current_celltype_for_ref].unique().tolist()
            ref_options = ["Auto select"] + all_celltypes
            ref_idx = ref_options.index(current_reference) if current_reference in ref_options else 0
            reference_celltype = st.selectbox(
                "Reference cell type (assumed unchanged)",
                ref_options,
                index=ref_idx,
                help="Cell type assumed unchanged. With Auto select, the most stable cell type is chosen."
            )

            # Inclusion probability threshold (used instead of FDR in scCODA)
            inclusion_threshold = st.slider(
                "Inclusion probability threshold",
                min_value=0.5,
                max_value=1.0,
                value=current_threshold,
                step=0.05,
                help="Threshold for significant effect. 0.95 means effect is non-zero in 95% of MCMC samples"
            )

        # Submit button
        params_submitted = st.form_submit_button("✅ Confirm Parameters", type="primary")

        if params_submitted:
            st.session_state.sccoda_params = {
                "celltype_col": celltype_col,
                "sample_col": sample_col,
                "condition_col": condition_col,
                "reference_celltype": reference_celltype,
                "inclusion_threshold": inclusion_threshold
            }
            st.success("✅ Parameters confirmed")

    # Use confirmed parameters or defaults
    celltype_col = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    sample_col = st.session_state.sccoda_params.get("sample_col", default_sample)
    condition_col = st.session_state.sccoda_params.get("condition_col", default_condition)
    reference_celltype = st.session_state.sccoda_params.get("reference_celltype", "Auto select")
    inclusion_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Formula input (outside form for flexibility)
    st.markdown("---")
    formula = st.text_input(
        "Model formula (optional)",
        value=condition_col,
        help="Default is condition column only. Specify 'condition + batch' to include additional covariates"
    )

    # Display confirmed parameters
    st.markdown("### 📌 Confirmed Parameters")
    param_col1, param_col2, param_col3 = st.columns(3)
    with param_col1:
        st.info(f"**Cell type column:** `{celltype_col}`")
        st.info(f"**Sample column:** `{sample_col}`")
    with param_col2:
        st.info(f"**Condition column:** `{condition_col}`")
        st.info(f"**Reference cell type:** `{reference_celltype}`")
    with param_col3:
        st.info(f"**Model formula:** `{formula}`")
        st.info(f"**Inclusion threshold:** `{inclusion_threshold}`")

    # Show data summary
    with st.expander("📊 Summary of selected data (confirmed)"):
        st.markdown(f"""
        **Cell type column:** `{celltype_col}`
        - Unique cell types: {adata.obs[celltype_col].nunique()}
        - Cell types: {', '.join(map(str, adata.obs[celltype_col].unique()[:10].tolist()))}{'...' if adata.obs[celltype_col].nunique() > 10 else ''}

        **Sample column:** `{sample_col}`
        - Number of samples: {adata.obs[sample_col].nunique()}
        - Samples: {', '.join(map(str, adata.obs[sample_col].unique().tolist()))}

        **Condition column:** `{condition_col}`
        - Number of conditions: {adata.obs[condition_col].nunique()}
        - Conditions: {', '.join(map(str, adata.obs[condition_col].unique().tolist()))}
        """)

        # Cross-tabulation
        st.markdown("**Sample × Condition cross tabulation:**")
        crosstab = pd.crosstab(adata.obs[sample_col], adata.obs[condition_col])
        st.dataframe(crosstab)

    # Visualize cell type distributions
    st.markdown("---")
    if st.button("📊 Set conditions and check data distribution", type="secondary"):
        with st.spinner("Visualizing data distribution..."):
            try:
                # Prepare data for visualization
                sccoda_vis = pt.tl.Sccoda()
                adata_vis = sccoda_vis.load(
                    adata,
                    type="cell_level",
                    cell_type_identifier=celltype_col,
                    sample_identifier=sample_col,
                    covariate_obs=[condition_col]
                )

                st.subheader("📈 Cell type distribution")

                # Reset Step 2 plots
                st.session_state.sccoda_plots["step2"] = []

                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("#### Cell type abundance (by condition)")
                    try:
                        import matplotlib.pyplot as plt
                        fig = sccoda_vis.plot_boxplots(
                            adata_vis,
                            modality_key="coda",
                            feature_name=condition_col,
                            add_dots=True
                        )
                        if fig is None:
                            fig = plt.gcf()
                        st.pyplot(fig)
                        # Save figure for PDF
                        st.session_state.sccoda_plots["step2"].append(("Cell type abundance (by condition)", fig))
                    except Exception as e:
                        st.warning(f"Failed to create boxplot: {str(e)}")

                with col2:
                    st.markdown("#### Cell type composition (by sample)")
                    try:
                        import matplotlib.pyplot as plt
                        fig = sccoda_vis.plot_stacked_barplot(
                            adata_vis,
                            modality_key="coda",
                            feature_name=sample_col
                        )
                        if fig is None:
                            fig = plt.gcf()
                        st.pyplot(fig)
                        # Save figure for PDF
                        st.session_state.sccoda_plots["step2"].append(("Cell type composition (by sample)", fig))
                    except Exception as e:
                        st.warning(f"Failed to create stacked bar plot: {str(e)}")

                # Additional plot: composition by condition (cell-level calculation)
                st.markdown("#### Cell type composition (by condition, cell-level)")
                st.caption("* Count all cells to calculate ratio (same method as original barplot)")
                try:
                    import matplotlib.pyplot as plt

                    # Calculate cell-level proportions (same as original barplot)
                    cell_counts = pd.crosstab(
                        adata.obs[condition_col],
                        adata.obs[celltype_col],
                        normalize='index'
                    ) * 100

                    # Create stacked bar plot
                    fig, ax = plt.subplots(figsize=(8, 6))
                    cell_counts.plot(kind='bar', stacked=True, ax=ax, width=0.8)
                    ax.set_ylabel('Proportion (%)')
                    ax.set_xlabel(condition_col)
                    ax.set_title(f'Cell Type Composition by {condition_col}')
                    ax.legend(title=celltype_col, bbox_to_anchor=(1.05, 1), loc='upper left')
                    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                    plt.tight_layout()

                    st.pyplot(fig)
                    # Save figure for PDF
                    st.session_state.sccoda_plots["step2"].append(("Cell type composition (by condition, cell-level)", fig))
                except Exception as e:
                    st.warning(f"Failed to create cell-level composition plot: {str(e)}")

                # scCODA sample-level average plot (for comparison)
                with st.expander("📊 Show sample-level mean (scCODA method)"):
                    st.caption("* Calculate ratio for each sample, then average within condition (method used in statistical analysis)")
                    try:
                        import matplotlib.pyplot as plt
                        fig = sccoda_vis.plot_stacked_barplot(
                            adata_vis,
                            modality_key="coda",
                            feature_name=condition_col
                        )
                        if fig is None:
                            fig = plt.gcf()
                        st.pyplot(fig)
                        # Save figure for PDF
                        st.session_state.sccoda_plots["step2"].append(("Cell type composition (by condition, sample mean)", fig))
                    except Exception as e:
                        st.warning(f"Failed to create sample-level plot: {str(e)}")

                # Download plots as PDF ZIP
                if len(st.session_state.sccoda_plots["step2"]) > 0:
                    st.markdown("---")
                    st.markdown("### 💾 Download plots")

                    # Create PDF in memory
                    pdf_buffer = io.BytesIO()
                    with PdfPages(pdf_buffer) as pdf:
                        for title, fig in st.session_state.sccoda_plots["step2"]:
                            pdf.savefig(fig, bbox_inches='tight')

                    pdf_buffer.seek(0)

                    st.download_button(
                        label="📥 Download data distribution plots as PDF",
                        data=pdf_buffer,
                        file_name="sccoda_step2_plots.pdf",
                        mime="application/pdf"
                    )

            except Exception as e:
                st.error(f"❌ Visualization error: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # ========================================
    # Step 3: Run scCODA
    # ========================================
    st.header("Step 3: Run scCODA")

    # Initialize result session states
    if "sccoda_effect_df" not in st.session_state:
        st.session_state.sccoda_effect_df = None
    if "sccoda_intercept_df" not in st.session_state:
        st.session_state.sccoda_intercept_df = None
    if "sccoda_cell_counts" not in st.session_state:
        st.session_state.sccoda_cell_counts = None
    if "sccoda_result_path" not in st.session_state:
        st.session_state.sccoda_result_path = None
    if "sccoda_analysis_params" not in st.session_state:
        st.session_state.sccoda_analysis_params = None

    if st.button("🚀 Run scCODA analysis", type="primary"):
        with st.spinner("Running scCODA..."):
            try:
                # Prepare data for scCODA
                st.info("📊 Preparing data...")

                # Create scCODA object
                sccoda = pt.tl.Sccoda()

                # Load data into scCODA format
                adata_sccoda = sccoda.load(
                    adata,
                    type="cell_level",
                    cell_type_identifier=celltype_col,
                    sample_identifier=sample_col,
                    covariate_obs=[condition_col]
                )

                st.success(f"✅ Created scCODA data: {adata_sccoda.shape}")

                # Fix: Copy condition column to coda modality if missing
                if hasattr(adata_sccoda, 'mod') and 'coda' in adata_sccoda.mod:
                    coda_mod = adata_sccoda.mod['coda']
                    if condition_col not in coda_mod.obs.columns:
                        rna_mod = adata_sccoda.mod['rna']
                        if condition_col in rna_mod.obs.columns and sample_col in rna_mod.obs.columns:
                            sample_condition = rna_mod.obs.groupby(sample_col)[condition_col].first()
                            coda_mod.obs[condition_col] = coda_mod.obs.index.map(sample_condition)
                            st.info(f"📝 Copied '{condition_col}' to coda modality")

                # Set reference cell type
                if reference_celltype != "Auto select":
                    st.info(f"📌 Reference cell types: {reference_celltype}")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type=reference_celltype
                    )
                else:
                    st.info("🔄 Auto-selecting reference cell type...")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type="automatic"
                    )

                # Run scCODA
                st.info("🔬 Running Bayesian estimation (may take several minutes)...")
                sccoda.run_nuts(
                    adata_sccoda,
                    num_warmup=1000,
                    num_samples=1000
                )

                st.success("✅ scCODA analysis completed!")

                # Store results in session state
                st.session_state.sccoda_results = adata_sccoda

                # Get and store effect_df
                try:
                    effect_df = sccoda.get_effect_df(adata_sccoda)
                    st.session_state.sccoda_effect_df = effect_df
                except Exception as e:
                    st.warning(f"Failed to get effect_df: {str(e)}")
                    st.session_state.sccoda_effect_df = None

                # Get and store intercept_df
                try:
                    intercept_df = sccoda.get_intercept_df(adata_sccoda)
                    st.session_state.sccoda_intercept_df = intercept_df
                except:
                    st.session_state.sccoda_intercept_df = None

                # Calculate and store cell-level proportions
                try:
                    cell_counts = pd.crosstab(
                        adata.obs[condition_col],
                        adata.obs[celltype_col],
                        normalize='index'
                    ) * 100
                    st.session_state.sccoda_cell_counts = cell_counts
                except:
                    st.session_state.sccoda_cell_counts = None

                # Save h5ad and store path
                result_path = os.path.join(sccoda_temp_dir, "sccoda_results.h5ad")
                adata_sccoda.write(result_path)
                st.session_state.sccoda_result_path = result_path

                # Store analysis parameters for display
                st.session_state.sccoda_analysis_params = {
                    "celltype_col": celltype_col,
                    "sample_col": sample_col,
                    "condition_col": condition_col,
                    "reference_celltype": reference_celltype,
                    "formula": formula,
                    "inclusion_threshold": inclusion_threshold
                }

                # Reset plots
                st.session_state.sccoda_plots["step3"] = []

                st.rerun()

            except Exception as e:
                st.error(f"❌ scCODA analysis error: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # ========================================
    # Display Results (outside button block - persists after rerun)
    # ========================================
    if st.session_state.sccoda_effect_df is not None:
        effect_df = st.session_state.sccoda_effect_df
        params = st.session_state.sccoda_analysis_params or {}

        st.subheader("📈 Analysis Results")

        # Effect summary
        st.markdown("### Effects by cell type")

        st.info("""
        **Interpreting scCODA results:**
        - **Final Parameter**: Posterior mean estimate of effect (0 = not significant, non-zero = significant)
        - **Inclusion probability**: Probability that effect is non-zero (higher = more certain change)
        - **log2-fold change**: Log2-transformed fold change
        - **HDI**: Highest Density Interval (Bayesian confidence interval)
        """)

        st.dataframe(effect_df, use_container_width=True)

        # Identify significant cell types
        if 'Final Parameter' in effect_df.columns:
            significant = effect_df[effect_df['Final Parameter'] != 0]

            if 'Inclusion probability' in effect_df.columns:
                thresh = params.get("inclusion_threshold", 0.95)
                significant_high_prob = significant[significant['Inclusion probability'] >= thresh]

                if len(significant_high_prob) > 0:
                    st.success(f"✨ **{len(significant_high_prob)}Significant changes detected in  cell types** (Inclusion prob ≥ {thresh})")
                    st.markdown("#### Cell types with significant changes")
                    st.dataframe(significant_high_prob, use_container_width=True)
                else:
                    st.warning(f"Inclusion probability ≥ {thresh}  - no cell types found")
                    if len(significant) > 0:
                        st.info(f"However, {len(significant)} cell types with changes detected (below threshold)")
            else:
                if len(significant) > 0:
                    st.success(f"✨ **{len(significant)}Significant changes detected in  cell types**")
                    st.dataframe(significant, use_container_width=True)
                else:
                    st.info("No cell types with significant changes detected")

        # Intercept parameters
        if st.session_state.sccoda_intercept_df is not None:
            with st.expander("📊 Intercept parameters"):
                st.dataframe(st.session_state.sccoda_intercept_df, use_container_width=True)

        # Visualization
        st.subheader("📊 Visualization")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("#### Cell type proportion changes")
            if 'Final Parameter' in effect_df.columns and (effect_df['Final Parameter'] != 0).any():
                st.info("Effect plot for cell types with significant changes")
                # Note: Cannot regenerate pertpy plot without sccoda object
                # Would need to store the figure in session_state during analysis
            else:
                st.info("No effect plot displayed because no significant changes")

        with col2:
            st.markdown("#### Cell type composition (by condition, cell-level)")
            if st.session_state.sccoda_cell_counts is not None:
                try:
                    import matplotlib.pyplot as plt
                    cell_counts = st.session_state.sccoda_cell_counts
                    cond_col = params.get("condition_col", "condition")
                    ct_col = params.get("celltype_col", "celltype")

                    fig, ax = plt.subplots(figsize=(8, 6))
                    cell_counts.plot(kind='bar', stacked=True, ax=ax, width=0.8)
                    ax.set_ylabel('Proportion (%)')
                    ax.set_xlabel(cond_col)
                    ax.set_title(f'Cell Type Composition by {cond_col}')
                    ax.legend(title=ct_col, bbox_to_anchor=(1.05, 1), loc='upper left')
                    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
                    plt.tight_layout()

                    st.pyplot(fig)
                    plt.close(fig)
                except Exception as e:
                    st.warning(f"Failed to create composition plot: {str(e)}")

        # Download results
        st.subheader("💾 Download results")

        col1, col2, col3 = st.columns(3)
        with col1:
            csv = effect_df.to_csv(index=True)
            st.download_button(
                label="📥 Download as CSV",
                data=csv,
                file_name="sccoda_results.csv",
                mime="text/csv",
                key="download_csv"
            )
        with col2:
            tsv = effect_df.to_csv(index=True, sep='\t')
            st.download_button(
                label="📥 Download as TSV",
                data=tsv,
                file_name="sccoda_results.tsv",
                mime="text/tab-separated-values",
                key="download_tsv"
            )
        with col3:
            if st.session_state.sccoda_result_path and os.path.exists(st.session_state.sccoda_result_path):
                with open(st.session_state.sccoda_result_path, "rb") as f:
                    st.download_button(
                        label="📥 Download as h5ad",
                        data=f,
                        file_name="sccoda_results.h5ad",
                        mime="application/octet-stream",
                        key="download_h5ad"
                    )

        # Clear results button
        if st.button("🗑️ Clear results", type="secondary"):
            st.session_state.sccoda_effect_df = None
            st.session_state.sccoda_intercept_df = None
            st.session_state.sccoda_cell_counts = None
            st.session_state.sccoda_result_path = None
            st.session_state.sccoda_analysis_params = None
            st.session_state.sccoda_results = None
            st.rerun()

# ========================================
# Additional information
# ========================================
with st.expander("ℹ️ About scCODA"):
    st.markdown("""
    ## What is scCODA?

    scCODA (single-cell Compositional Data Analysis) is a hierarchical Bayesian model for statistically detecting **cell type composition changes** in
    single-cell data.

    ### Key Features

    1. **Compositional data support**
       - Cell type proportions have compositional constraint (sum=1)
       - Conventional statistical methods are inappropriate
       - scCODA is specialized for compositional data

    2. **Sparse estimation**
       - Detects only cell types that actually change
       - Reduces false positives

    3. **Hierarchical model**
       - Accounts for variation between samples
       - Distinguishes biological and technical variation

    4. **Reference cell type**
       - Select cell type assumed unchanged
       - Evaluate changes in other cell types relative to this

    5. **Bayesian inference**
       - Uses **Inclusion probability** instead of p-value or FDR
       - Inclusion probability is proportion of MCMC samples where effect is non-zero
       - Usually 0.95+ is significant (95% probability of change)

    ### When to use?

    - ✅ **Want to compare cell type proportions between 2+ conditions**
    - ✅ **Have multiple biological replicates**
    - ✅ **Need statistical analysis considering compositional data characteristics**

    ### Comparison with other methods

    | Method | Data | Sample variation | Compositional | Sparse | Statistical method |
    |------|--------|--------------|---------|-------------|---------|
    | **scCODA** | Count | ✅ | ✅ | ✅ | Bayesian (Inclusion prob) |
    | **Propeller** | Count | ✅ | ❌ | ❌ | Frequentist (p-value) |
    | **Beta-binomial** | Count | ✅ | ❌ | ❌ | Frequentist (p-value/LRT) |
    | **χ² test** | Count | ❌ | ❌ | ❌ | Frequentist (p-value) |

    ### References

    - Büttner, M., Ostner, J., Müller, C.L. et al. (2021)
      "scCODA is a Bayesian model for compositional single-cell data analysis"
      *Nature Communications* **12**, 6876
      [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)
    """)
