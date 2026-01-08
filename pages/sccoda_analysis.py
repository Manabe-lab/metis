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
scCODAを用いて、シングルセルデータにおける**セルタイプ構成比の変化**を統計的に検出します。

### scCODAの特徴
- **合成データ (Compositional Data) に対応**: 比率データの合成制約を適切に処理
- **階層的ベイズモデル**: グループ間のセルタイプ比率変化を推定
- **スパース推定**: 実際に変化しているセルタイプのみを検出
- **参照セルタイプ**: 変化していないと仮定するセルタイプを選択可能
- **ベイズ統計**: p値やFDRではなく、Inclusion probabilityで有意性を判定

### ワークフロー
1. **h5adファイルアップロード**: セルタイプ・サンプル・条件情報を含むデータ
2. **パラメータ設定**: セルタイプ列、サンプル列、条件列、参照セルタイプを選択
3. **scCODA実行**: ベイズ推定による差次的構成比解析
4. **結果表示**: 有意に変化したセルタイプの可視化とダウンロード

### 参考文献
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
        "step2": [],  # Step 2のプロット
        "step3": []   # Step 3のプロット
    }

# ========================================
# Step 1: Upload h5ad file
# ========================================
st.header("Step 1: Upload h5ad file")

st.markdown("""
### 必要な情報
h5adファイル (`adata`) には以下の情報が必要です：

- **セルタイプ情報**: `adata.obs` に列として含まれる（例: `'cell_type'`, `'clusters'`, `'leiden'`）
- **サンプル情報**: `adata.obs` に列として含まれる（例: `'sample'`, `'batch'`, `'replicate'`）
- **条件情報**: `adata.obs` に列として含まれる（例: `'condition'`, `'group'`, `'treatment'`）

💡 **ヒント**: Seuratで解析したデータは、Pythonで`anndata`形式に変換してからアップロードしてください。
""")

uploaded_h5ad = st.file_uploader(
    "h5adファイルをアップロード",
    type=["h5ad"],
    help="セルタイプ、サンプル、条件情報を含むh5adファイル"
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
        with st.spinner("h5adファイルを読み込み中..."):
            try:
                adata = load_h5ad_cached(file_content, uploaded_h5ad.name)
                st.session_state.sccoda_adata = adata
                st.session_state.sccoda_file_hash = file_hash
                # Reset parameters when new file is loaded
                if "sccoda_params" in st.session_state:
                    del st.session_state.sccoda_params
                # Reset results
                st.session_state.sccoda_results = None

                st.success(f"✅ h5adファイルを読み込みました: {adata.shape[0]} cells × {adata.shape[1]} genes")
            except Exception as e:
                st.error(f"❌ h5adファイルの読み込みエラー: {str(e)}")
                import traceback
                st.code(traceback.format_exc())
                st.stop()
    else:
        adata = st.session_state.sccoda_adata
        st.success(f"✅ h5adファイル（キャッシュ済み）: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # Display data summary
    adata = st.session_state.sccoda_adata

    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("総細胞数", f"{adata.shape[0]:,}")
    with col2:
        st.metric("遺伝子数", f"{adata.shape[1]:,}")
    with col3:
        if hasattr(adata, 'obs') and len(adata.obs.columns) > 0:
            st.metric("メタデータ列数", len(adata.obs.columns))

    # Show available metadata columns
    with st.expander("📋 利用可能なメタデータ列"):
        st.write("**adata.obs に含まれる列:**")
        st.dataframe(pd.DataFrame({
            "列名": adata.obs.columns.tolist(),
            "データ型": [str(dtype) for dtype in adata.obs.dtypes],
            "ユニーク値数": [adata.obs[col].nunique() for col in adata.obs.columns]
        }))

    # Preview data
    with st.expander("🔍 データプレビュー (先頭10行)"):
        st.dataframe(adata.obs.head(10))

# ========================================
# Step 2: Configure scCODA parameters
# ========================================
if st.session_state.sccoda_adata is not None:
    st.header("Step 2: Configure scCODA parameters")

    adata = st.session_state.sccoda_adata
    obs_columns = adata.obs.columns.tolist()

    # フィルタリング: 要素数50未満の列のみを候補とする
    filtered_columns = [col for col in obs_columns if adata.obs[col].nunique() < 50]

    if len(filtered_columns) == 0:
        st.error("❌ 要素数が50未満の列が見つかりません。メタデータの列を確認してください。")
        st.stop()

    # デフォルト列の選択関数
    def find_default_celltype_col(columns, obs_df):
        """セルタイプ列のdefaultを見つける"""
        # 1. identを含む列（orig.identを除く）
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 2. typeを含む列
        type_cols = [col for col in columns if 'type' in col.lower()]
        if type_cols:
            return type_cols[0]
        # 3. デフォルトは最初の列
        return columns[0]

    def find_default_sample_col(columns, obs_df):
        """サンプル列のdefaultを見つける"""
        # 1. orig.identを優先
        if 'orig.ident' in columns:
            return 'orig.ident'
        # 2. identを含む列（orig.ident以外）
        ident_cols = [col for col in columns if 'ident' in col.lower() and col != 'orig.ident']
        if ident_cols:
            return ident_cols[0]
        # 3. sampleを含む列
        sample_cols = [col for col in columns if 'sample' in col.lower()]
        if sample_cols:
            return sample_cols[0]
        # 4. デフォルトは最初の列
        return columns[0]

    def find_default_condition_col(columns, obs_df):
        """条件列のdefaultを見つける"""
        # condition, stim, KOを含む列を探す
        keywords = ['condition', 'stim', 'ko']
        for keyword in keywords:
            matching_cols = [col for col in columns if keyword in col.lower()]
            if matching_cols:
                return matching_cols[0]
        # デフォルトは2番目の列（最初はサンプルとして使われる可能性が高い）
        return columns[min(1, len(columns)-1)]

    # デフォルト列を決定
    default_celltype = find_default_celltype_col(filtered_columns, adata.obs)
    default_sample = find_default_sample_col(filtered_columns, adata.obs)
    default_condition = find_default_condition_col(filtered_columns, adata.obs)

    # インデックスを取得
    celltype_idx = filtered_columns.index(default_celltype) if default_celltype in filtered_columns else 0
    sample_idx = filtered_columns.index(default_sample) if default_sample in filtered_columns else 0
    condition_idx = filtered_columns.index(default_condition) if default_condition in filtered_columns else min(1, len(filtered_columns)-1)

    # Initialize session state for form values
    if "sccoda_params" not in st.session_state:
        st.session_state.sccoda_params = {
            "celltype_col": default_celltype,
            "sample_col": default_sample,
            "condition_col": default_condition,
            "reference_celltype": "自動選択",
            "inclusion_threshold": 0.95
        }

    # Get current values from session state for form defaults
    current_celltype = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    current_sample = st.session_state.sccoda_params.get("sample_col", default_sample)
    current_condition = st.session_state.sccoda_params.get("condition_col", default_condition)
    current_reference = st.session_state.sccoda_params.get("reference_celltype", "自動選択")
    current_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Calculate indices for selectboxes based on current values
    form_celltype_idx = filtered_columns.index(current_celltype) if current_celltype in filtered_columns else celltype_idx
    form_sample_idx = filtered_columns.index(current_sample) if current_sample in filtered_columns else sample_idx
    form_condition_idx = filtered_columns.index(current_condition) if current_condition in filtered_columns else condition_idx

    # Form for parameter selection
    st.info("💡 パラメータを選択後、「✅ パラメータを確定」ボタンをクリックしてください")
    with st.form("sccoda_params_form"):
        col1, col2 = st.columns(2)

        with col1:
            st.subheader("データ列の選択")

            # Cell type column
            celltype_col = st.selectbox(
                "セルタイプ列",
                filtered_columns,
                index=form_celltype_idx,
                help="セルタイプまたはクラスター情報を含む列"
            )

            # Sample column
            sample_col = st.selectbox(
                "サンプル列",
                filtered_columns,
                index=form_sample_idx,
                help="生物学的レプリケート（サンプル）情報を含む列"
            )

            # Condition column
            condition_col = st.selectbox(
                "条件列（比較するグループ）",
                filtered_columns,
                index=form_condition_idx,
                help="比較したい条件（例: Control vs Treatment）"
            )

        with col2:
            st.subheader("モデルパラメータ")

            # Reference cell type - use celltypes from current celltype column
            current_celltype_for_ref = filtered_columns[form_celltype_idx]
            all_celltypes = adata.obs[current_celltype_for_ref].unique().tolist()
            ref_options = ["自動選択"] + all_celltypes
            ref_idx = ref_options.index(current_reference) if current_reference in ref_options else 0
            reference_celltype = st.selectbox(
                "参照セルタイプ（変化しないと仮定）",
                ref_options,
                index=ref_idx,
                help="変化していないと仮定するセルタイプ。自動選択の場合、最も安定したセルタイプが選ばれます。"
            )

            # Inclusion probability threshold (scCODAではFDRの代わりに使用)
            inclusion_threshold = st.slider(
                "Inclusion probability閾値",
                min_value=0.5,
                max_value=1.0,
                value=current_threshold,
                step=0.05,
                help="効果が有意と判定されるための閾値。0.95はMCMCサンプルの95%でこの効果が0でないことを意味します"
            )

        # Submit button
        params_submitted = st.form_submit_button("✅ パラメータを確定", type="primary")

        if params_submitted:
            st.session_state.sccoda_params = {
                "celltype_col": celltype_col,
                "sample_col": sample_col,
                "condition_col": condition_col,
                "reference_celltype": reference_celltype,
                "inclusion_threshold": inclusion_threshold
            }
            st.success("✅ パラメータを確定しました")

    # Use confirmed parameters or defaults
    celltype_col = st.session_state.sccoda_params.get("celltype_col", default_celltype)
    sample_col = st.session_state.sccoda_params.get("sample_col", default_sample)
    condition_col = st.session_state.sccoda_params.get("condition_col", default_condition)
    reference_celltype = st.session_state.sccoda_params.get("reference_celltype", "自動選択")
    inclusion_threshold = st.session_state.sccoda_params.get("inclusion_threshold", 0.95)

    # Formula input (outside form for flexibility)
    st.markdown("---")
    formula = st.text_input(
        "モデル式（オプション）",
        value=condition_col,
        help="デフォルトは条件列のみ。追加の共変量を含める場合は 'condition + batch' のように指定"
    )

    # Display confirmed parameters
    st.markdown("### 📌 確定済みパラメータ")
    param_col1, param_col2, param_col3 = st.columns(3)
    with param_col1:
        st.info(f"**セルタイプ列:** `{celltype_col}`")
        st.info(f"**サンプル列:** `{sample_col}`")
    with param_col2:
        st.info(f"**条件列:** `{condition_col}`")
        st.info(f"**参照セルタイプ:** `{reference_celltype}`")
    with param_col3:
        st.info(f"**モデル式:** `{formula}`")
        st.info(f"**Inclusion閾値:** `{inclusion_threshold}`")

    # Show data summary
    with st.expander("📊 選択したデータの要約（確定済み）"):
        st.markdown(f"""
        **セルタイプ列:** `{celltype_col}`
        - ユニークなセルタイプ数: {adata.obs[celltype_col].nunique()}
        - セルタイプ: {', '.join(map(str, adata.obs[celltype_col].unique()[:10].tolist()))}{'...' if adata.obs[celltype_col].nunique() > 10 else ''}

        **サンプル列:** `{sample_col}`
        - サンプル数: {adata.obs[sample_col].nunique()}
        - サンプル: {', '.join(map(str, adata.obs[sample_col].unique().tolist()))}

        **条件列:** `{condition_col}`
        - 条件数: {adata.obs[condition_col].nunique()}
        - 条件: {', '.join(map(str, adata.obs[condition_col].unique().tolist()))}
        """)

        # Cross-tabulation
        st.markdown("**サンプル×条件のクロス集計:**")
        crosstab = pd.crosstab(adata.obs[sample_col], adata.obs[condition_col])
        st.dataframe(crosstab)

    # Visualize cell type distributions
    st.markdown("---")
    if st.button("📊 条件を設定してデータ分布を確認", type="secondary"):
        with st.spinner("データ分布を可視化中..."):
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

                st.subheader("📈 細胞タイプの分布")

                # Reset Step 2 plots
                st.session_state.sccoda_plots["step2"] = []

                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("#### セルタイプ存在量（条件別）")
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
                        st.session_state.sccoda_plots["step2"].append(("セルタイプ存在量（条件別）", fig))
                    except Exception as e:
                        st.warning(f"ボックスプロットの作成に失敗: {str(e)}")

                with col2:
                    st.markdown("#### セルタイプ構成比（サンプル別）")
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
                        st.session_state.sccoda_plots["step2"].append(("セルタイプ構成比（サンプル別）", fig))
                    except Exception as e:
                        st.warning(f"スタックバープロットの作成に失敗: {str(e)}")

                # Additional plot: composition by condition (cell-level calculation)
                st.markdown("#### セルタイプ構成比（条件別・細胞レベル集計）")
                st.caption("※ 全細胞をカウントして比率を計算（元のBarplotと同じ方法）")
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
                    st.session_state.sccoda_plots["step2"].append(("セルタイプ構成比（条件別・細胞レベル）", fig))
                except Exception as e:
                    st.warning(f"細胞レベル構成比プロットの作成に失敗: {str(e)}")

                # scCODA sample-level average plot (for comparison)
                with st.expander("📊 サンプルレベル平均（scCODA方式）を表示"):
                    st.caption("※ 各サンプルの比率を計算後、条件内で平均（統計解析で使用する方法）")
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
                        st.session_state.sccoda_plots["step2"].append(("セルタイプ構成比（条件別・サンプル平均）", fig))
                    except Exception as e:
                        st.warning(f"サンプルレベルプロットの作成に失敗: {str(e)}")

                # Download plots as PDF ZIP
                if len(st.session_state.sccoda_plots["step2"]) > 0:
                    st.markdown("---")
                    st.markdown("### 💾 プロットのダウンロード")

                    # Create PDF in memory
                    pdf_buffer = io.BytesIO()
                    with PdfPages(pdf_buffer) as pdf:
                        for title, fig in st.session_state.sccoda_plots["step2"]:
                            pdf.savefig(fig, bbox_inches='tight')

                    pdf_buffer.seek(0)

                    st.download_button(
                        label="📥 データ分布プロットをPDFでダウンロード",
                        data=pdf_buffer,
                        file_name="sccoda_step2_plots.pdf",
                        mime="application/pdf"
                    )

            except Exception as e:
                st.error(f"❌ 可視化エラー: {str(e)}")
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

    if st.button("🚀 scCODA解析を実行", type="primary"):
        with st.spinner("scCODAを実行中..."):
            try:
                # Prepare data for scCODA
                st.info("📊 データを準備中...")

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

                st.success(f"✅ scCODAデータを作成: {adata_sccoda.shape}")

                # Fix: Copy condition column to coda modality if missing
                if hasattr(adata_sccoda, 'mod') and 'coda' in adata_sccoda.mod:
                    coda_mod = adata_sccoda.mod['coda']
                    if condition_col not in coda_mod.obs.columns:
                        rna_mod = adata_sccoda.mod['rna']
                        if condition_col in rna_mod.obs.columns and sample_col in rna_mod.obs.columns:
                            sample_condition = rna_mod.obs.groupby(sample_col)[condition_col].first()
                            coda_mod.obs[condition_col] = coda_mod.obs.index.map(sample_condition)
                            st.info(f"📝 '{condition_col}' を coda modality にコピーしました")

                # Set reference cell type
                if reference_celltype != "自動選択":
                    st.info(f"📌 参照セルタイプ: {reference_celltype}")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type=reference_celltype
                    )
                else:
                    st.info("🔄 参照セルタイプを自動選択中...")
                    sccoda.prepare(
                        adata_sccoda,
                        formula=formula,
                        reference_cell_type="automatic"
                    )

                # Run scCODA
                st.info("🔬 ベイズ推定を実行中（数分かかる場合があります）...")
                sccoda.run_nuts(
                    adata_sccoda,
                    num_warmup=1000,
                    num_samples=1000
                )

                st.success("✅ scCODA解析が完了しました！")

                # Store results in session state
                st.session_state.sccoda_results = adata_sccoda

                # Get and store effect_df
                try:
                    effect_df = sccoda.get_effect_df(adata_sccoda)
                    st.session_state.sccoda_effect_df = effect_df
                except Exception as e:
                    st.warning(f"effect_dfの取得に失敗: {str(e)}")
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
                st.error(f"❌ scCODA解析エラー: {str(e)}")
                import traceback
                st.code(traceback.format_exc())

    # ========================================
    # Display Results (outside button block - persists after rerun)
    # ========================================
    if st.session_state.sccoda_effect_df is not None:
        effect_df = st.session_state.sccoda_effect_df
        params = st.session_state.sccoda_analysis_params or {}

        st.subheader("📈 解析結果")

        # Effect summary
        st.markdown("### セルタイプごとの効果")

        st.info("""
        **scCODA結果の解釈:**
        - **Final Parameter**: 効果の事後平均推定値（0なら有意でない、非0なら有意）
        - **Inclusion probability**: 効果が0でない確率（高いほど確実に変化している）
        - **log2-fold change**: log2変換された倍率変化
        - **HDI**: 高密度信用区間（ベイズ版の信頼区間）
        """)

        st.dataframe(effect_df, use_container_width=True)

        # Identify significant cell types
        if 'Final Parameter' in effect_df.columns:
            significant = effect_df[effect_df['Final Parameter'] != 0]

            if 'Inclusion probability' in effect_df.columns:
                thresh = params.get("inclusion_threshold", 0.95)
                significant_high_prob = significant[significant['Inclusion probability'] >= thresh]

                if len(significant_high_prob) > 0:
                    st.success(f"✨ **{len(significant_high_prob)}個のセルタイプで有意な変化が検出されました** (Inclusion prob ≥ {thresh})")
                    st.markdown("#### 有意な変化を示すセルタイプ")
                    st.dataframe(significant_high_prob, use_container_width=True)
                else:
                    st.warning(f"Inclusion probability ≥ {thresh} を満たすセルタイプはありませんでした")
                    if len(significant) > 0:
                        st.info(f"ただし、{len(significant)}個のセルタイプで変化が検出されています（閾値以下）")
            else:
                if len(significant) > 0:
                    st.success(f"✨ **{len(significant)}個のセルタイプで有意な変化が検出されました**")
                    st.dataframe(significant, use_container_width=True)
                else:
                    st.info("有意な変化が検出されたセルタイプはありませんでした")

        # Intercept parameters
        if st.session_state.sccoda_intercept_df is not None:
            with st.expander("📊 切片パラメータ"):
                st.dataframe(st.session_state.sccoda_intercept_df, use_container_width=True)

        # Visualization
        st.subheader("📊 可視化")

        col1, col2 = st.columns(2)

        with col1:
            st.markdown("#### セルタイプ比率の変化")
            if 'Final Parameter' in effect_df.columns and (effect_df['Final Parameter'] != 0).any():
                st.info("有意な変化があるセルタイプの効果プロット")
                # Note: Cannot regenerate pertpy plot without sccoda object
                # Would need to store the figure in session_state during analysis
            else:
                st.info("有意な変化がないため、効果プロットは表示されません")

        with col2:
            st.markdown("#### セルタイプ構成比（条件別・細胞レベル）")
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
                    st.warning(f"構成比プロットの作成に失敗: {str(e)}")

        # Download results
        st.subheader("💾 結果のダウンロード")

        col1, col2, col3 = st.columns(3)
        with col1:
            csv = effect_df.to_csv(index=True)
            st.download_button(
                label="📥 CSVでダウンロード",
                data=csv,
                file_name="sccoda_results.csv",
                mime="text/csv",
                key="download_csv"
            )
        with col2:
            tsv = effect_df.to_csv(index=True, sep='\t')
            st.download_button(
                label="📥 TSVでダウンロード",
                data=tsv,
                file_name="sccoda_results.tsv",
                mime="text/tab-separated-values",
                key="download_tsv"
            )
        with col3:
            if st.session_state.sccoda_result_path and os.path.exists(st.session_state.sccoda_result_path):
                with open(st.session_state.sccoda_result_path, "rb") as f:
                    st.download_button(
                        label="📥 h5adでダウンロード",
                        data=f,
                        file_name="sccoda_results.h5ad",
                        mime="application/octet-stream",
                        key="download_h5ad"
                    )

        # Clear results button
        if st.button("🗑️ 結果をクリア", type="secondary"):
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
with st.expander("ℹ️ scCODAについて"):
    st.markdown("""
    ## scCODAとは？

    scCODA (single-cell Compositional Data Analysis) は、シングルセルデータにおける**セルタイプ構成比の変化**を
    統計的に検出するための階層ベイズモデルです。

    ### 主な特徴

    1. **合成データへの対応**
       - セルタイプ比率は合成制約（総和=1）を持つ
       - 通常の統計手法では不適切
       - scCODAは合成データに特化した手法

    2. **スパース推定**
       - 実際に変化しているセルタイプのみを検出
       - 偽陽性を減らす

    3. **階層モデル**
       - サンプル間の変動を考慮
       - 生物学的変動と技術的変動を区別

    4. **参照セルタイプ**
       - 変化していないと仮定するセルタイプを選択
       - 他のセルタイプの変化をこれを基準に評価

    5. **ベイズ統計による判定**
       - p値やFDRではなく**Inclusion probability**を使用
       - Inclusion probabilityはMCMCサンプルのうち効果が0でない割合
       - 通常0.95以上を有意とする（95%の確率で変化している）

    ### いつ使うべきか？

    - ✅ **2群以上の条件間でセルタイプ比率を比較したい**
    - ✅ **複数の生物学的レプリケートがある**
    - ✅ **合成データの特性を考慮した統計解析が必要**

    ### 他の手法との比較

    | 手法 | データ | サンプル間変動 | 合成制約 | スパース推定 | 統計手法 |
    |------|--------|--------------|---------|-------------|---------|
    | **scCODA** | Count | ✅ | ✅ | ✅ | ベイズ (Inclusion prob) |
    | **Propeller** | Count | ✅ | ❌ | ❌ | 頻度論 (p値) |
    | **Beta-binomial** | Count | ✅ | ❌ | ❌ | 頻度論 (p値/LRT) |
    | **χ² test** | Count | ❌ | ❌ | ❌ | 頻度論 (p値) |

    ### 参考文献

    - Büttner, M., Ostner, J., Müller, C.L. et al. (2021)
      "scCODA is a Bayesian model for compositional single-cell data analysis"
      *Nature Communications* **12**, 6876
      [doi:10.1038/s41467-021-27150-6](https://doi.org/10.1038/s41467-021-27150-6)
    """)
