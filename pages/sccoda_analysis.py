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
scCODAtheuseite、shinguruseruDatatookeru**serutaipustructbecomeratioofchangeize**theStatisticalaltocheckoutshimasu。

### scCODAofspecfeature
- **matchbecomeData (Compositional Data) topairrespond**: ratiorateDataofmatchbecomecontrolaboutthefitcuttoprocproc
- **tierlayeralbeizumoderu**: GroupbetweenofserutaipuratioratechangeizetheEstimation
- **supa-suEstimation**: realocctochangeizeshiteexistserutaipuofmithecheckout
- **refrefserutaipu**: changeizeshiteinotandtempsetdoserutaiputheSelectpossible
- **beizuStatistical**: pvalyaFDRwithisnaku、Inclusion probabilitywithSignificantnaturethejudgeset

### wa-kufuro-
1. **h5adFileUpload**: serutaipu,Sample,ConditioninfoinfotheincludemuData
2. **ParameterSettings**: serutaipucol、Samplecol、Conditioncol、refrefserutaiputheSelect
3. **scCODARun**: beizuEstimationtoyorudiffnextalstructbecomeratioAnalysis
4. **ResultDisplay**: SignificanttochangeizeshitaserutaipuofVisualizationandDownload

### refthinktextdedic
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
        "step2": [],  # Step 2ofpuroto
        "step3": []   # Step 3ofpuroto
    }

# ========================================
# Step 1: Upload h5ad file
# ========================================
st.header("Step 1: Upload h5ad file")

st.markdown("""
### requirednainfoinfo
h5adFile (`adata`) toisorbelowofinfoinfoisrequiredwithsu：

- **serutaipuinfoinfo**: `adata.obs` tocolandshiteincludemareru（Example: `'cell_type'`, `'clusters'`, `'leiden'`）
- **Sampleinfoinfo**: `adata.obs` tocolandshiteincludemareru（Example: `'sample'`, `'batch'`, `'replicate'`）
- **Conditioninfoinfo**: `adata.obs` tocolandshiteincludemareru（Example: `'condition'`, `'group'`, `'treatment'`）

💡 **hinto**: SeuratwithAnalysisshitaDatais、Pythonwith`anndata`shapeformattochangechangeshitefromUploadshitekudasai。
""")

uploaded_h5ad = st.file_uploader(
    "h5adFiletheUpload",
    type=["h5ad"],
    help="serutaipu、Sample、Conditioninfoinfotheincludemuh5adFile"
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
        with st.spinner("h5adFiletheLoadingmid..."):
            try:
                adata = load_h5ad_cached(file_content, uploaded_h5ad.name)
                st.session_state.sccoda_adata = adata
                st.session_state.sccoda_file_hash = file_hash
                # Reset parameters when new file is loaded
                if "sccoda_params" in st.session_state:
                    del st.session_state.sccoda_params
                # Reset results
                st.session_state.sccoda_results = None

                st.success(f"✅ h5adFiletheLoadingmashita: {adata.shape[0]} cells × {adata.shape[1]} genes")
            except Exception as e:
                st.error(f"❌ h5adFileofLoadingError: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
                st.stop()
    else:
        adata = st.session_state.sccoda_adata
        st.success(f"✅ h5adFile（kiyashiyudonemi）: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Display data summary
    adata = st.session_state.sccoda_adata

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("totalCellnum", f"{adata.shape[0]:,}")
    with col2:
        st.metric("Genenum", f"{adata.shape[1]:,}")
    with col3:
        if hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
            st.metric("metaDatacolnum", len(adata.obs.columns))

    # Show available metadata columns
    with st.expander("📋 useusepossiblenametaDatacol"):
        st.write("**adata.obs toincludemarerucol:**")
        st.dataframe(pd.DataFrame({
            "colname": adata.obs.columns.tolist(),
            "Datatype": [str(dtype) for dtype in adata.obs.dtypes],
            "yuni-kuvalnum": [adata.obs[col].nunique() for col in adata.obs.columns]
        }))

    # Preview data
    with st.expander("🔍 Datapurebiyu- (firsthead10row)"):
        st.dataframe(adata.obs.head(10))

# ========================================
# Step 2: Configure scCODA parameters
# ========================================
if st.session_state.sccoda_adata is not None:
    st.header("Step 2: Configure scCODA parameters")

    adata = st.session_state.sccoda_adata
    obs_columns = adata.obs.columns.tolist()

    # Filtering: needelemnum50not yetfullofcolofmithe候suppanddo
    filtered_columns = [col for col in obs_columns if adata.obs[col].nunique() < 50]

    if len(filtered_columns) == 0:
        st.error("❌ needelemnumis50not yetfullofcolisviewtsukarimasen。metaDataofcoltheConfirmshitekudasai。")
        st.stop()

    # defuorutocolofSelectrelnum
    def find_default_celltype_col(columns, obs_df):
        """serutaipucolofdefaulttheviewtsukeru"""
        # 1. identtheincludemucol（orig.identtheremoveku）
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 2. typetheincludemucol
        type_cols = [col for col in columns if 'type' in col.lower()]
        if type_cols:
            return type_cols[0]
        # 3. defuorutoismostfirstofcol
        return columns[0]

    def find_default_sample_col(columns, obs_df):
        """Samplecolofdefaulttheviewtsukeru"""
        # 1. orig.identthepriorfirst
        if 'orig.ident' in columns:
            return 'orig.ident'
        # 2. identtheincludemucol（orig.identorout）
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 3. sampletheincludemucol
        sample_cols = [col for col in columns if 'sample' in col.lower()]
        if sample_cols:
            return sample_cols[0]
        # 4. defuorutoismostfirstofcol
        return columns[0]

    def find_default_condition_col(columns, obs_df):
        """Conditioncolofdefaulttheviewtsukeru"""
        # condition, stim, KOtheincludemucolthesearchsu
        keywords = ['condition', 'stim', 'ko']
        for keyword in keywords:
            matching_cols = [col for col in columns if keyword in col.lower()]
            if matching_cols:
                return matching_cols[0]
        # defuorutois2numidxofcol（mostfirstisSampleandshiteusewarerupossiblenatureishighi）
        return columns[min(1, len(columns)-1)]

    # defuorutocolthedetset
    default_celltype = find_default_celltype_col(filtered_columns, adata.obs)
    default_sample = find_default_sample_col(filtered_columns, adata.obs)
    default_condition = find_default_condition_col(filtered_columns, adata.obs)

    # indekusuthegetget
    celltype_idx = filtered_columns.index(default_celltype) if default_celltype in filtered_columns else 0
    sample_idx = filtered_columns.index(default_sample) if default_sample in filtered_columns else 0
    condition_idx = filtered_columns.index(default_condition) if default_condition in filtered_columns else min(1, len(filtered_columns)-1)

    # Initialize session state for form values
    if "sccoda_params" not in st.session_state:
        st.session_state.sccoda_params = {
            "celltype_col": default_celltype,
            "sample_col": default_sample,
            "condition_col": default_condition,
            "reference_celltype": "selfmoveSelect",
            "inclusion_threshold": 0.95
        }

    # Get current values from session state for form defaults
    current_celltype = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    current_sample = st.session_state.sccoda_params.get("sample_col", default_sample)
    current_condition = st.session_state.sccoda_params.get("condition_col", default_condition)
    current_reference = st.session_state.sccoda_params.get("reference_celltype", "selfmoveSelect")
    current_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Calculate indices for selectboxes based on current values
    form_celltype_idx = filtered_columns.index(current_celltype) if current_celltype in filtered_columns else celltype_idx
    form_sample_idx = filtered_columns.index(current_sample) if current_sample in filtered_columns else sample_idx
    form_condition_idx = filtered_columns.index(current_condition) if current_condition in filtered_columns else condition_idx

    # Form for parameter selection
    st.info("💡 ParametertheSelectafter、「✅ Parameterthecertainset」botanthekurikushitekudasai")
    with st.form("sccoda_params_form"):
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("DatacolofSelect")

            # Cell type column
            celltype_col = st.selectbox(
                "serutaipucol",
                filtered_columns,
                index=form_celltype_idx,
                help="serutaipualsoisClusterinfoinfotheincludemucol"
            )

            # Sample column
            sample_col = st.selectbox(
                "Samplecol",
                filtered_columns,
                index=form_sample_idx,
                help="genthinglearnalrepurike-to（Sample）infoinfotheincludemucol"
            )

            # Condition column
            condition_col = st.selectbox(
                "Conditioncol（ComparisondoGroup）",
                filtered_columns,
                index=form_condition_idx,
                help="ComparisonshitaiCondition（Example: Control vs Treatment）"
            )

        with col2:
            st.subheader("moderuParameter")

            # Reference cell type - use celltypes from current celltype column
            current_celltype_for_ref = filtered_columns[form_celltype_idx]
            all_celltypes = adata.obs[current_celltype_for_ref].unique().tolist()
            ref_options = ["selfmoveSelect"] + all_celltypes
            ref_idx = ref_options.index(current_reference) if current_reference in ref_options else 0
            reference_celltype = st.selectbox(
                "refrefserutaipu（changeizeshinotandtempset）",
                ref_options,
                index=ref_idx,
                help="changeizeshiteinotandtempsetdoserutaipu。selfmoveSelectofplacematch、mostmosafesetshitaserutaipuisselbaremasu。"
            )

            # Inclusion probability threshold (scCODAwithisFDRofsubwaritouseuse)
            inclusion_threshold = st.slider(
                "Inclusion probabilityThreshold",
                min_value=0.5,
                max_value=1.0,
                value=current_threshold,
                step=0.05,
                help="effectresultisSignificantandjudgesetsareruforofThreshold。0.95isMCMCSampleof95%withkoofeffectresultis0withnotkoandthemeanmeanshimasu"
            )

        # Submit button
        params_submitted = st.form_submit_button("✅ Parameterthecertainset", type="primary")

        if params_submitted:
            st.session_state.sccoda_params = {
                "celltype_col": celltype_col,
                "sample_col": sample_col,
                "condition_col": condition_col,
                "reference_celltype": reference_celltype,
                "inclusion_threshold": inclusion_threshold
            }
            st.success("✅ Parameterthecertainsetshimashita")

    # Use confirmed parameters or defaults
    celltype_col = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    sample_col = st.session_state.sccoda_params.get("sample_col", default_sample)
    condition_col = st.session_state.sccoda_params.get("condition_col", default_condition)
    reference_celltype = st.session_state.sccoda_params.get("reference_celltype", "selfmoveSelect")
    inclusion_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Formula input (outside form for flexibility)
    st.markdown("---")
    formula = st.text_input(
        "moderuformat（Option）",
        value=condition_col,
        help="defuorutoisConditioncolofmi。addaddofcochangeamounttheincludemeruplacematchis 'condition + batch' ofliketopointset"
    )

    # Display confirmed parameters
    st.markdown("### 📌 certainsetdonemiParameter")
    param_col1, param_col2, param_col3 = st.columns(3)
    with param_col1:
        st.info(f"**serutaipucol:** `{celltype_col}`")
        st.info(f"**Samplecol:** `{sample_col}`")
    with param_col2:
        st.info(f"**Conditioncol:** `{condition_col}`")
        st.info(f"**refrefserutaipu:** `{reference_celltype}`")
    with param_col3:
        st.info(f"**moderuformat:** `{formula}`")
        st.info(f"**InclusionThreshold:** `{inclusion_threshold}`")

    # Show data summary
    with st.expander("📊 SelectshitaDataofneedabout（certainsetdonemi）"):
        st.markdown(f"""
        **serutaipucol:** `{celltype_col}`
        - yuni-kunaserutaipunum: {adata.obs[celltype_col].nunique()}
        - serutaipu: {', '.join(map(str, adata.obs[celltype_col].unique()[:10].tolist()))}{'...' if adata.obs[celltype_col].nunique() > 10 else ''}

        **Samplecol:** `{sample_col}`
        - Samplenum: {adata.obs[sample_col].nunique()}
        - Sample: {', '.join(map(str, adata.obs[sample_col].unique().tolist()))}

        **Conditioncol:** `{condition_col}`
        - Conditionnum: {adata.obs[condition_col].nunique()}
        - Condition: {', '.join(map(str, adata.obs[condition_col].unique().tolist()))}
        """)

        # Cross-tabulation
        st.markdown("**Sample×Conditionofkurosugathercalc:**")
        crosstab = pd.crosstab(adata.obs[sample_col], adata.obs[condition_col])
        st.dataframe(crosstab)

    # Visualize cell type distributions
    st.markdown("---")
    if st.button("📊 ConditiontheSettingsshiteDatadivdisttheConfirm", type="secondary"):
        with st.spinner("DatadivdisttheVisualizationmid..."):
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

                st.subheader("📈 Cell typeofdivdist")

                # Reset Step 2 plots
                st.session_state.sccoda_plots["step2"] = []

                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("#### serutaipuexistatamount（Conditionsep）")
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
                        st.session_state.sccoda_plots["step2"].append(("serutaipuexistatamount（Conditionsep）", fig))
                    except Exception as e:
                        st.warning(f"bokusupurotoofmakebecometofailfail: {str(e)}")

                with col2:
                    st.markdown("#### serutaipustructbecomeratio（Samplesep）")
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
                        st.session_state.sccoda_plots["step2"].append(("serutaipustructbecomeratio（Samplesep）", fig))
                    except Exception as e:
                        st.warning(f"sutakuba-purotoofmakebecometofailfail: {str(e)}")

                # Additional plot: composition by condition (cell-level calculation)
                st.markdown("#### serutaipustructbecomeratio（Conditionsep,Cellreberugathercalc）")
                st.caption("※ allCellthekauntoshiteratioratetheCalculation（sourceofBarplotandsamewaymethod）")
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
                    st.session_state.sccoda_plots["step2"].append(("serutaipustructbecomeratio（Conditionsep,Cellreberu）", fig))
                except Exception as e:
                    st.warning(f"Cellreberustructbecomeratiopurotoofmakebecometofailfail: {str(e)}")

                # scCODA sample-level average plot (for comparison)
                with st.expander("📊 Samplereberuflatavg（scCODAwayformat）theDisplay"):
                    st.caption("※ eachSampleofratioratetheCalculationafter、Conditioninwithflatavg（StatisticalAnalysiswithuseusedowaymethod）")
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
                        st.session_state.sccoda_plots["step2"].append(("serutaipustructbecomeratio（Conditionsep,Sampleflatavg）", fig))
                    except Exception as e:
                        st.warning(f"Samplereberupurotoofmakebecometofailfail: {str(e)}")

                # Download plots as PDF ZIP
                if len(st.session_state.sccoda_plots["step2"]) > 0:
                    st.markdown("---")
                    st.markdown("### 💾 purotoofDownload")

                    # Create PDF in memory
                    pdf_buffer = io.BytesIO()
                    with PdfPages(pdf_buffer) as pdf:
                        for title, fig in st.session_state.sccoda_plots["step2"]:
                            pdf.savefig(fig, bbox_inches='tight')

                    pdf_buffer.seek(0)

                    st.download_button(
                        label="📥 DatadivdistpurotothePDFwithDownload",
                        data=pdf_buffer,
                        file_name="sccoda_step2_plots.pdf",
                        mime="application/pdf"
                    )

            except Exception as e:
                st.error(f"❌ VisualizationError: {str(e)}")
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

    if st.button("🚀 scCODAAnalysistheRun", type="primary"):
        with st.spinner("scCODAtheRunmid..."):
            try:
                # Prepare data for scCODA
                st.info("📊 Datathelevelprepmid...")

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

                st.success(f"✅ scCODADatathemakebecome: {adata_sccoda.shape}")

                # Fix: Copy condition column to coda modality if missing
                if hasattr(adata_sccoda, 'mod') and 'coda' in adata_sccoda.mod:
                    coda_mod = adata_sccoda.mod['coda']
                    if condition_col not in coda_mod.obs.columns:
                        rna_mod = adata_sccoda.mod['rna']
                        if condition_col in rna_mod.obs.columns and sample_col in rna_mod.obs.columns:
                            sample_condition = rna_mod.obs.groupby(sample_col)[condition_col].first()
                            coda_mod.obs[condition_col] = coda_mod.obs.index.map(sample_condition)
                            st.info(f"📝 '{condition_col}' the coda modality tokopi-shimashita")

                # Set reference cell type
                if reference_celltype != "selfmoveSelect":
                    st.info(f"📌 refrefserutaipu: {reference_celltype}")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type=reference_celltype
                    )
                else:
                    st.info("🔄 refrefserutaiputheselfmoveSelectmid...")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type="automatic"
                    )

                # Run scCODA
                st.info("🔬 beizuEstimationtheRunmid（numdivkakaruplacematchisarimasu）...")
                sccoda.run_nuts(
                    adata_sccoda,
                    num_warmup=1000,
                    num_samples=1000
                )

                st.success("✅ scCODAAnalysisisCompleteshimashita！")

                # Store results in session state
                st.session_state.sccoda_results = adata_sccoda

                # Get and store effect_df
                try:
                    effect_df = sccoda.get_effect_df(adata_sccoda)
                    st.session_state.sccoda_effect_df = effect_df
                except Exception as e:
                    st.warning(f"effect_dfofgetgettofailfail: {str(e)}")
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
                st.error(f"❌ scCODAAnalysisError: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # ========================================
    # Display Results (outside button block - persists after rerun)
    # ========================================
    if st.session_state.sccoda_effect_df is not None:
        effect_df = st.session_state.sccoda_effect_df
        params = st.session_state.sccoda_analysis_params or {}

        st.subheader("📈 AnalysisResult")

        # Effect summary
        st.markdown("### serutaipugoandofeffectresult")

        st.info("""
        **scCODAResultofsolveinterp:**
        - **Final Parameter**: effectresultofthingafterflatavgEstimationval（0naraSignificantwithnot、non0naraSignificant）
        - **Inclusion probability**: effectresultis0withnotcertainrate（highihodocertainrealtochangeizeshiteexist）
        - **log2-fold change**: log2changechangesaretatimesratechangeize
        - **HDI**: highdensedegreetrustuseareabetween（beizuveroftrustrelyareabetween）
        """)

        st.dataframe(effect_df, use_container_width=True)

        # Identify significant cell types
        if 'Final Parameter' in effect_df.columns:
            significant = effect_df[effect_df['Final Parameter'] != 0]

            if 'Inclusion probability' in effect_df.columns:
                thresh = params.get("inclusion_threshold", 0.95)
                significant_high_prob = significant[significant['Inclusion probability'] >= thresh]

                if len(significant_high_prob) > 0:
                    st.success(f"✨ **{len(significant_high_prob)}pieceofserutaipuwithSignificantnachangeizeischeckoutsaremashita** (Inclusion prob ≥ {thresh})")
                    st.markdown("#### Significantnachangeizetheshowsuserutaipu")
                    st.dataframe(significant_high_prob, use_container_width=True)
                else:
                    st.warning(f"Inclusion probability ≥ {thresh} thefulltasuserutaipuisarimasenwithshita")
                    if len(significant) > 0:
                        st.info(f"however、{len(significant)}pieceofserutaipuwithchangeizeischeckoutsareteimasu（Thresholdorbelow）")
            else:
                if len(significant) > 0:
                    st.success(f"✨ **{len(significant)}pieceofserutaipuwithSignificantnachangeizeischeckoutsaremashita**")
                    st.dataframe(significant, use_container_width=True)
                else:
                    st.info("Significantnachangeizeischeckoutsaretaserutaipuisarimasenwithshita")

        # Intercept parameters
        if st.session_state.sccoda_intercept_df is not None:
            with st.expander("📊 cut片Parameter"):
                st.dataframe(st.session_state.sccoda_intercept_df, use_container_width=True)

        # Visualization
        st.subheader("📊 Visualization")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("#### serutaipuratiorateofchangeize")
            if 'Final Parameter' in effect_df.columns and (effect_df['Final Parameter'] != 0).any():
                st.info("Significantnachangeizeisexistserutaipuofeffectresultpuroto")
                # Note: Cannot regenerate pertpy plot without sccoda object
                # Would need to store the figure in session_state during analysis
            else:
                st.info("Significantnachangeizeisnotfor、effectresultpurotoisDisplaysaremasen")

        with col2:
            st.markdown("#### serutaipustructbecomeratio（Conditionsep,Cellreberu）")
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
                    st.warning(f"structbecomeratiopurotoofmakebecometofailfail: {str(e)}")

        # Download results
        st.subheader("💾 ResultofDownload")

        col1, col2, col3 = st.columns(3)
        with col1:
            csv = effect_df.to_csv(index=True)
            st.download_button(
                label="📥 CSVwithDownload",
                data=csv,
                file_name="sccoda_results.csv",
                mime="text/csv",
                key="download_csv"
            )
        with col2:
            tsv = effect_df.to_csv(index=True, sep='\t')
            st.download_button(
                label="📥 TSVwithDownload",
                data=tsv,
                file_name="sccoda_results.tsv",
                mime="text/tab-separated-values",
                key="download_tsv"
            )
        with col3:
            if st.session_state.sccoda_result_path and os.path.exists(st.session_state.sccoda_result_path):
                with open(st.session_state.sccoda_result_path, "rb") as f:
                    st.download_button(
                        label="📥 h5adwithDownload",
                        data=f,
                        file_name="sccoda_results.h5ad",
                        mime="application/octet-stream",
                        key="download_h5ad"
                    )

        # Clear results button
        if st.button("🗑️ Resultthekuria", type="secondary"):
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
with st.expander("ℹ️ scCODAtotsuite"):
    st.markdown("""
    ## scCODAandis？

    scCODA (single-cell Compositional Data Analysis) is、shinguruseruDatatookeru**serutaipustructbecomeratioofchangeize**the
    Statisticalaltocheckoutdoforoftierlayerbeizumoderuwithsu。

    ### mainnaspecfeature

    1. **matchbecomeDataheofpairrespond**
       - serutaipuratiorateismatchbecomecontrolabout（totalsum=1）theholdtsu
       - passnormalofStatisticalhandmethodwithisnotfitcut
       - scCODAismatchbecomeDatatospecizeshitahandmethod

    2. **supa-suEstimation**
       - realocctochangeizeshiteexistserutaipuofmithecheckout
       - falsesunnaturethereducerasu

    3. **tierlayermoderu**
       - Samplebetweenofchangemovethethinkconsider
       - genthinglearnalchangemoveand技術alchangemovetheareasep

    4. **refrefserutaipu**
       - changeizeshiteinotandtempsetdoserutaiputheSelect
       - otherofserutaipuofchangeizethekorethebaseleveltoevalval

    5. **beizuStatisticaltoyorujudgeset**
       - pvalyaFDRwithisnaku**Inclusion probability**theuseuse
       - Inclusion probabilityisMCMCSampleofuchieffectresultis0withnotratiomatch
       - passnormal0.95oruptheSignificantanddo（95%ofcertainratewithchangeizeshiteexist）

    ### itsuuseubekika？

    - ✅ **2grouporupofConditionbetweenwithserutaipuratioratetheComparisonshitai**
    - ✅ **multinumofgenthinglearnalrepurike-toisexist**
    - ✅ **matchbecomeDataofspecnaturethethinkconsidershitaStatisticalAnalysisisrequired**

    ### otherofhandmethodandofComparison

    | handmethod | Data | Samplebetweenchangemove | matchbecomecontrolabout | supa-suEstimation | Statisticalhandmethod |
    |------|--------|--------------|---------|-------------|---------|
    | **scCODA** | Count | ✅ | ✅ | ✅ | beizu (Inclusion prob) |
    | **Propeller** | Count | ✅ | ❌ | ❌ | freqdegreelogic (pval) |
    | **Beta-binomial** | Count | ✅ | ❌ | ❌ | freqdegreelogic (pval/LRT) |
    | **χ² test** | Count | ❌ | ❌ | ❌ | freqdegreelogic (pval) |

    ### refthinktextdedic

    - Büttner, M., Ostner, J., Müller, C.L. et al. (2021)
      "scCODA is a Bayesian model for compositional single-cell data analysis"
      *Nature Communications* **12**, 6876
      [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)
    """)
