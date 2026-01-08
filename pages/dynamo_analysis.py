"""
Dynamo RNA Velocity Analysis
Perform advanced RNA velocity and vector field analysis using Dynamo
"""

import streamlit as st
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd
import os
import io
import tempfile
import time
from helper_func import clear_old_directories, clear_old_files


def parse_gene_list(genes_str):
    """
    Parse gene list from string with multiple delimiters (space, comma, tab, newline).
    Case-insensitive. Returns list of unique gene names (uppercase).
    """
    if not genes_str or len(genes_str.strip()) == 0:
        return []

    # Remove quotes
    genes_str = genes_str.replace("'", "")
    genes_str = genes_str.replace('"', "")

    # Split by space first
    gene_list = genes_str.split(' ')
    gene_list = list(filter(lambda a: a != '', gene_list))  # Remove empty strings

    # Split by comma
    if ',' in genes_str:
        gene_list = sum([x.split(',') for x in gene_list], [])

    # Split by tab
    if '\t' in genes_str:
        gene_list = sum([x.split('\t') for x in gene_list], [])

    # Split by newline
    if '\n' in genes_str:
        gene_list = sum([x.split('\n') for x in gene_list], [])

    # Remove empty strings and strip whitespace
    gene_list = [g.strip() for g in gene_list if g.strip()]

    # Convert to uppercase for case-insensitive matching
    gene_list = [g.upper() for g in gene_list]

    # Remove duplicates while preserving order
    gene_list = sorted(set(gene_list), key=gene_list.index)

    return gene_list

# Import dynamo
try:
    import dynamo as dyn
    DYNAMO_AVAILABLE = True
except ImportError:
    DYNAMO_AVAILABLE = False

st.set_page_config(page_title="Dynamo Analysis", page_icon="🌊", layout="wide")

st.title("🌊 Dynamo RNA Velocity Analysis")

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
Dynamoを用いて高度なRNA velocity解析とベクトル場解析を実行します。

### ワークフロー
1. **ファイル読み込み**: h5ad (optional) + loom (spliced/unspliced matrices)
2. **前処理** (オプション): scVelo解析済みh5adの場合はスキップ推奨
3. **Dynamics推定**: RNA dynamicsのモデリング（モーメント計算を含む）
4. **Vector Field**: 連続的なベクトル場の再構築
5. **幾何学的解析**: Speed, curl, divergence, accelerationの計算
6. **結果保存**: 解析結果をh5adファイルでダウンロード

### Dynamoの特徴
- **Vector field reconstruction**: 連続的な速度場による軌跡予測
- **Differential geometry**: Jacobian, acceleration, curvatureなどの力学的解析
- **Least action paths**: 細胞状態間の最適遷移経路の推定
- **In silico perturbation**: 遺伝子摂動の効果予測
- **Regulatory network**: 転写制御ネットワークの推論

### 参考
- [Qiu et al. (2022) "Mapping transcriptomic vector fields of single cells" Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4)
- [Dynamo Documentation](https://dynamo-release.readthedocs.io/)
- [Dynamo GitHub](https://github.com/aristoteleo/dynamo-release)
""")

# Initialize session state
if "dynamo_temp_dir" not in st.session_state:
    dynamo_temp_dir = os.path.join("temp", f"dynamo_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(dynamo_temp_dir, exist_ok=True)
    st.session_state.dynamo_temp_dir = dynamo_temp_dir
else:
    dynamo_temp_dir = st.session_state.dynamo_temp_dir

if "dynamo_complete" not in st.session_state:
    st.session_state.dynamo_complete = False

# ========================================
# Step 1: Upload files
# ========================================
st.header("Step 1: Upload files")

st.markdown("""
### ファイルの準備

#### 🌟 推奨ワークフロー（デフォルト）

**h5adファイル + loomファイル** の両方をアップロード:

1. **h5adファイル**: scVelo解析済みファイル
   - ✅ UMAP embedding (`adata.obsm`に`umap`を含むキー、例: `X_umap`, `umap`, `X_umap_2d`)
   - ✅ クラスタリング結果 (`adata.obs['leiden']`, `adata.obs['louvain']`, `adata.obs['clusters']`等)
   - ✅ QC/正規化済みデータ
   - → これらの既存情報を活用し、前処理をスキップ

2. **loomファイル**: Data filtering appで**フィルタリング済み**のloomファイル
   - spliced/unspliced行列を含む
   - h5adと同じ細胞IDを持つ

#### ⚠️ loomファイルのみから解析する場合

- **Data filtering appでフィルタリングしていない生loomファイル**を使用
- Dynamoの前処理（QC、正規化、特徴選択、PCA）を実行
- UMAP/クラスター情報は新規に計算
""")

col1, col2 = st.columns(2)

with col1:
    uploaded_h5ad = st.file_uploader(
        "Upload h5ad file (推奨)",
        type=['h5ad'],
        key="dynamo_h5ad_upload",
        help="scVelo解析済みh5adファイル（UMAP/クラスター/QC済み）"
    )

with col2:
    uploaded_loom = st.file_uploader(
        "Upload loom file (必須)",
        type=['loom'],
        key="dynamo_loom_upload",
        help="Data filtering appでフィルタリング済みのloomファイル（または生loom）"
    )

# Cache for file loading - use file content hash for proper cache invalidation
import hashlib

def _get_file_hash(file_obj):
    """Get MD5 hash of file content for cache key"""
    content = file_obj.getvalue()
    return hashlib.md5(content).hexdigest()

@st.cache_data
def load_h5ad_file(_file_obj, _file_hash):
    """Load h5ad file with caching based on file content hash"""
    import io
    temp_path = os.path.join("temp", f"h5ad_temp_{round(time.time())}.h5ad")
    with open(temp_path, "wb") as f:
        f.write(_file_obj.getvalue())
    adata = sc.read_h5ad(temp_path)
    os.remove(temp_path)
    return adata

@st.cache_data
def load_loom_file(_file_obj, _file_hash):
    """Load loom file with caching based on file content hash"""
    import io
    temp_path = os.path.join("temp", f"loom_temp_{round(time.time())}.loom")
    with open(temp_path, "wb") as f:
        f.write(_file_obj.getvalue())
    adata = ad.read_loom(temp_path)
    os.remove(temp_path)
    return adata

# Can proceed with just loom file
if uploaded_loom is not None:
    if uploaded_h5ad is not None:
        st.success("✓ 両方のファイルがアップロードされました - h5adのUMAP/クラスター情報を使用します（推奨）")

        # Load h5ad with cache
        adata_check = load_h5ad_file(uploaded_h5ad, _get_file_hash(uploaded_h5ad))

        # Extract available embeddings (X_umap, X_pca, etc.)
        available_embeddings = [key.replace('X_', '') for key in adata_check.obsm.keys() if key.startswith('X_')]
        if len(available_embeddings) == 0:
            available_embeddings = ['umap', 'pca']  # fallback
            st.warning("⚠️ h5adにembeddingが見つかりません。デフォルトのumap/pcaを使用します")
        else:
            st.info(f"✓ 利用可能なembeddings: {', '.join(available_embeddings)}")
    else:
        st.warning("⚠️ loomファイルのみ - Dynamoで新規にQC/正規化/UMAP計算を実行します（フィルタリングしていない生loomを使用してください）")
        available_embeddings = ['umap', 'pca']  # Default for loom-only
        adata_check = None

    # ========================================
    # Step 2: Basis Selection (outside form)
    # ========================================
    st.header("Step 2: Basis Selection")

    if uploaded_h5ad is not None:
        st.markdown("""
        **Vector Field計算用のbasis選択:**
        - 高次元遺伝子発現空間（PCA, MNN, Harmony等、50次元程度）を選択
        - Perturbation解析、Jacobian計算に使用

        ⚠️ **Perturbation解析の制約:**
        - **使用可能**: PCA、scVI、Harmony (PCAベース)
          - これらはローディング行列を持ち、遺伝子発現空間への逆変換が可能
        - **使用不可**: FastMNN (mnn)
          - FastMNNは新しい補正空間を作成し、ローディング行列を持たない
          - 遺伝子レベルの解釈ができないため、Perturbation解析では使えません
        - 💡 Perturbation解析を行う場合は、PCAまたはscVIベースのVector Fieldを必ず計算してください
        """)

        # Check for gene expression space embeddings
        pca_candidates = [b for b in available_embeddings if 'pca' in b.lower()]

        # Set default: prefer MNN/Harmony corrected spaces
        if pca_candidates:
            mnn_harmony_pca = [b for b in pca_candidates if 'mnn' in b.lower() or 'harmony' in b.lower()]
            default_basis = mnn_harmony_pca[0] if mnn_harmony_pca else pca_candidates[0]
        else:
            default_basis = available_embeddings[0]

        selected_basis = st.selectbox(
            "Basis for Vector Field computation:",
            available_embeddings,
            index=available_embeddings.index(default_basis),
            help="遺伝子発現空間 (pca, mnn, harmony等) の選択を推奨"
        )

        st.caption(f"💡 選択: **{selected_basis}** - Vector Fieldはこのbasisと、後で選択する2D embedding、新規X_dynamo.umapで計算されます")

        # Select corresponding UMAP/2D embedding
        st.markdown("---")
        st.markdown("**2D embedding選択（可視化用）:**")

        def levenshtein_distance(s1, s2):
            """Calculate Levenshtein distance between two strings"""
            if len(s1) < len(s2):
                return levenshtein_distance(s2, s1)

            if len(s2) == 0:
                return len(s1)

            previous_row = range(len(s2) + 1)
            for i, c1 in enumerate(s1):
                current_row = [i + 1]
                for j, c2 in enumerate(s2):
                    insertions = previous_row[j + 1] + 1
                    deletions = current_row[j] + 1
                    substitutions = previous_row[j] + (c1 != c2)
                    current_row.append(min(insertions, deletions, substitutions))
                previous_row = current_row

            return previous_row[-1]

        # Find best matching 2D embedding by Levenshtein distance
        umap_candidates = [b for b in available_embeddings if 'umap' in b.lower() or 'tsne' in b.lower()]

        # Calculate distances and find closest (excluding selected_basis itself)
        if umap_candidates:
            valid_candidates = [u for u in umap_candidates if u != selected_basis]

            if valid_candidates:
                # Calculate Levenshtein distance for each candidate
                distances = [(u, levenshtein_distance(selected_basis.lower(), u.lower()))
                             for u in valid_candidates]
                # Sort by distance (ascending) - closest match first
                distances.sort(key=lambda x: x[1])
                default_umap = distances[0][0]
            elif umap_candidates:
                # If all UMAP candidates equal selected_basis, use first
                default_umap = umap_candidates[0]
            else:
                default_umap = available_embeddings[0] if available_embeddings else 'umap'
        else:
            default_umap = available_embeddings[0] if available_embeddings else 'umap'

        selected_umap = st.selectbox(
            "2D embedding for visualization:",
            available_embeddings,
            index=available_embeddings.index(default_umap) if default_umap in available_embeddings else 0,
            help="可視化用の2D embedding。Topography解析等で使用。元の名前で保持され、X_umapにもコピーされます。"
        )

        st.info(f"""
        ✓ **Vector Field計算対象（3つ）:**
        1. **{selected_umap}** - 選択した2D embedding（元の名前で保持、X_umapにもコピー）
        2. **{selected_basis}** - 高次元遺伝子発現空間（Perturbation解析用）
        3. **X_dynamo.umap** - 新規計算（KNN整合、後で作成）
        """)
    else:
        selected_basis = None
        selected_umap = None
        st.info("loomファイルのみ: Vector Fieldは新規計算されるUMAPで実行されます")

    # ========================================
    # Step 3: Preprocessing options (outside form)
    # ========================================
    st.header("Step 3: Preprocessing options")

    col1, col2 = st.columns(2)
    with col1:
        skip_preprocess = st.checkbox(
            "Skip preprocessing",
            value=True,
            help="既にQC/正規化/PCA済みのh5adの場合はチェック（推奨）"
        )

    if skip_preprocess:
        with col2:
            recompute_hvg = st.checkbox(
                "Recompute HVGs",
                value=False,
                help="既存の正規化データを使いつつ、HVGだけ再計算する"
            )
    else:
        recompute_hvg = False

    # ========================================
    # Step 4: Configure analysis parameters (form)
    # ========================================
    st.header("Step 4: Configure analysis parameters")

    with st.expander("📚 Parameter Guide", expanded=False):
        st.markdown("""
        ### Preprocessing
        - **Skip preprocessing**: 既にQC/正規化/PCA済みのh5adファイルの場合はチェック
          - scVelo解析済みh5ad → スキップ推奨（既存の正規化/PCAを保持しつつDynamo用の最小限の構造を設定）
          - 生loomファイルのみ → 前処理を実行（QC、正規化、特徴選択、PCAを全て実行）
          - **注**: Dynamoは内部データ構造の設定が必要なため、完全にスキップはできません

        - **recipe** (前処理実行時): 前処理レシピ
          - **monocle** (推奨): Monocle3スタイル（QC、正規化、特徴選択、PCA）
          - **seurat**: Seuratスタイル

        - **Cell cycle score**: 細胞周期ステージの推定（前処理の有無に関わらず実行可能）
          - 用途: `dyn.pl.phase_diagram()`での細胞周期による色分け
          - 結果: `adata.obs`に細胞周期関連列が追加
          - **注**: h5adファイルから遺伝子名を取得（マウス/ヒト両対応）。h5adなしの場合は失敗する可能性あり（オプションなので解析は継続）

        ### Dynamics estimation (dyn.tl.dynamics)
        - **model**: Dynamicsモデル
          - **stochastic** (推奨): 確率的モデル、最も一般的
          - **deterministic**: 決定論的モデル
          - **degradation**: Degradation-onlyモデル
        - **est_method**: 推定法（stochasticの場合）
          - **gmm**: Gaussian mixture model（デフォルト、推奨）
          - **negbin**: Negative binomialモデル

        ### Vector Field (dyn.vf.VectorField)
        - **bases**: 計算基底（複数選択可）
          - **umap**: UMAPベースのベクトル場
            - 📊 **Topography可視化用** - curl計算、地形図可視化に最適
            - 2D embeddingが必須
          - **pca**: PCAベースのベクトル場
            - 🧬 **Perturbation解析用** - 遺伝子発現空間での定量解析
            - In silico perturbation、Jacobian解析、ranking/differential解析に必須
            - ⚠️ 後でPerturbation解析を行う場合は必ず選択してください
          - **mnn_umap, harmony_umap等**: バッチ補正後のUMAP（Topography用）
          - 💡 **推奨**: UMAPとPCA両方を選択（全ての解析に対応）
          - 各basisごとに独立してVector Fieldが保存されます（`VecFld_{basis}`）
          - 複数選択した場合、全てのbasisでVector Fieldが計算されます
        - **M**: サンプリングポイント数（デフォルト: 1000）
        - **pot_curl_div**: Potential, curl, divergenceの計算

        ### Differential Geometry
        - **speed**: 速度の大きさ
        - **curl**: 回転場（渦度）
        - **divergence**: 発散場（増殖/分化）
        - **acceleration**: 加速度
        - **curvature**: 曲率

        Qiu et al. (2022) Cell を参照
        """)

    with st.form("dynamo_params_form"):
        # HVG parameters (dynamic display based on preprocessing options)
        if not skip_preprocess or recompute_hvg:
            st.subheader("HVG parameters")

            n_top_genes = st.slider(
                "Number of highly variable genes",
                min_value=500,
                max_value=10000,
                value=2000,
                step=500,
                help="高変動遺伝子の数（デフォルト: 2000）"
            )

            # Force include specific genes
            st.markdown("##### Force include genes (comma, space, tab, CR separated):")
            force_genes = st.text_area(
                "Force include genes",
                value="",
                height=100,
                label_visibility='collapsed',
                help="必ずHVGに含めたい遺伝子名（複数の区切り文字対応: , space tab CR）"
            )
        elif skip_preprocess and not recompute_hvg:
            # Only option to add specific genes
            st.subheader("Additional genes")
            st.markdown("##### Add specific genes to existing HVGs (comma, space, tab, CR separated):")
            force_genes_only = st.text_area(
                "Add specific genes",
                value="",
                height=100,
                label_visibility='collapsed',
                help="既存のHVGに追加したい遺伝子名（複数の区切り文字対応: , space tab CR）"
            )
            n_top_genes = 2000  # Not used
            force_genes = ""
        else:
            n_top_genes = 2000
            force_genes = ""
            force_genes_only = ""

        # PCA components parameter
        st.subheader("PCA parameters")

        # Get the dimension of the selected embedding (if available)
        max_pca_components = 100  # Default max
        existing_embedding_dim = None
        if adata_check is not None and 'selected_basis' in dir() and selected_basis:
            embedding_key = f"X_{selected_basis}"
            if embedding_key in adata_check.obsm:
                existing_embedding_dim = adata_check.obsm[embedding_key].shape[1]
                max_pca_components = existing_embedding_dim

        n_pca_components = st.slider(
            "Number of PCA components",
            min_value=10,
            max_value=max_pca_components,
            value=min(30, max_pca_components),  # dynamo default is 30
            step=5,
            help="PCA成分数。dynamoデフォルトは30。Perturbation解析を行う場合、embeddingの次元数と一致させる必要があります。"
        )

        if existing_embedding_dim is not None:
            if n_pca_components == existing_embedding_dim:
                st.success(f"✓ PCA成分数 ({n_pca_components}) が既存のembedding (X_{selected_basis}) と一致しています")
            else:
                st.warning(f"⚠️ 既存のembedding (X_{selected_basis}) は {existing_embedding_dim} 次元です。Perturbation解析には {existing_embedding_dim} に設定してください。")
        else:
            st.info(f"💡 PCA成分数: {n_pca_components} (dynamoデフォルト: 30)")

        # Set preprocessing recipe based on skip_preprocess
        if not skip_preprocess:
            col1, col2 = st.columns(2)
            with col1:
                preprocess_recipe = st.selectbox(
                    "Preprocessing recipe",
                    ["monocle", "seurat"],
                    index=0,
                    help="前処理レシピ。monocleが推奨"
                )
            with col2:
                if preprocess_recipe == "monocle":
                    st.info("✨ **Monocle recipe** (Recommended)")
        else:
            preprocess_recipe = None
            force_genes_only = force_genes_only if not recompute_hvg else ""

        # 細胞周期スコアは前処理の有無に関わらず計算可能
        col1, col2 = st.columns(2)
        with col1 if skip_preprocess else col2:
            cell_cycle_score = st.checkbox(
                "Compute cell cycle score",
                value=True,
                help="細胞周期スコアを計算（phase diagram可視化用）"
            )

        st.subheader("Dynamics estimation")

        col1, col2 = st.columns(2)
        with col1:
            dynamics_model = st.selectbox(
                "Dynamics model",
                ["stochastic", "deterministic", "degradation"],
                index=0,
                help="Dynamicsモデル。stochasticが推奨"
            )

            if dynamics_model == "stochastic":
                st.info("✨ **Stochastic model** (Default / Recommended)")

        with col2:
            est_method = "gmm"  # Default
            if dynamics_model == "stochastic":
                est_method = st.selectbox(
                    "Estimation method",
                    ["gmm", "negbin"],
                    index=0,
                    help="推定方法。gmmが推奨"
                )
            else:
                st.info(f"Method: {dynamics_model} (fixed)")

            # Get number of available CPUs
            import multiprocessing
            n_cpus = multiprocessing.cpu_count()

            # Set default: 75% of available cores, max 24
            default_jobs = min(24, int(n_cpus * 0.75))

            n_jobs = st.number_input(
                "n_jobs (parallel cores)",
                min_value=1,
                max_value=n_cpus,
                value=default_jobs,
                help=f"並列計算のコア数。利用可能: {n_cpus}コア"
            )

        st.subheader("Vector field reconstruction")

        col1, col2, col3 = st.columns(3)
        with col1:
            compute_vector_field = st.checkbox(
                "Compute vector field",
                value=True,
                help="ベクトル場の再構築（推奨）"
            )

            if compute_vector_field:
                st.info("✨ **Enabled** (Recommended)")

        with col2:
            vf_bases = []
            if compute_vector_field:
                # Determine bases for Vector Field computation
                if uploaded_h5ad is not None and selected_basis and selected_umap:
                    # h5ad uploaded: use selected basis, selected umap, and dynamo.umap (will be computed)
                    vf_bases = [selected_umap, selected_basis, 'dynamo.umap']
                    st.info(f"""
                    ✓ **計算対象 (3 bases):**
                    1. {selected_umap}
                    2. {selected_basis}
                    3. dynamo.umap (後で計算)
                    """)
                else:
                    # loom-only: use dynamo.umap only
                    vf_bases = ['dynamo.umap']
                    st.info("✓ **計算対象:** dynamo.umap (新規計算)")

        with col3:
            vf_M = 1000
            if compute_vector_field:
                vf_M = st.number_input(
                    "Sampling points (M)",
                    min_value=100,
                    max_value=5000,
                    value=1000,
                    step=100,
                    help="サンプリングポイント数"
                )

        st.subheader("Differential geometry analysis")

        compute_geometry = st.checkbox(
            "Compute geometric features",
            value=True,
            help="Speed, curl, divergence, accelerationなどを計算"
        )

        if compute_geometry:
            if compute_vector_field and vf_bases:
                st.info(f"✨ **Enabled** - Vector Field basesで選択した全てのbasis ({len(vf_bases)}個) で計算します")
                st.caption(f"Geometry計算対象: {', '.join(vf_bases)}")
            else:
                st.warning("⚠️ Vector Fieldが有効でない場合、Geometry analysisは実行されません")

        st.markdown("---")

        check_config = st.form_submit_button("✓ Check Configuration", type="primary")

    # ========================================
    # Configuration Check
    # ========================================
    if check_config:
        st.markdown("---")
        st.subheader("🔍 Configuration Check")

        validation_passed = True

        # Check 1: Basis validation (h5ad only)
        if uploaded_h5ad is not None:
            if selected_basis:
                basis_key = f'X_{selected_basis}'
                if adata_check and basis_key not in adata_check.obsm:
                    st.error(f"""
                    ❌ **Basis validation failed**

                    Selected basis '{selected_basis}' not found in adata.obsm

                    **Available embeddings:** {', '.join(available_embeddings)}
                    """)
                    validation_passed = False
                else:
                    if adata_check:
                        basis_shape = adata_check.obsm[basis_key].shape
                        st.success(f"✅ Basis '{selected_basis}' validated (shape: {basis_shape})")

            if selected_umap:
                umap_key = f'X_{selected_umap}'
                if adata_check and umap_key not in adata_check.obsm:
                    st.error(f"""
                    ❌ **UMAP validation failed**

                    Selected UMAP '{selected_umap}' not found in adata.obsm

                    **Available embeddings:** {', '.join(available_embeddings)}
                    """)
                    validation_passed = False
                else:
                    if adata_check:
                        umap_shape = adata_check.obsm[umap_key].shape
                        st.success(f"✅ UMAP '{selected_umap}' validated (shape: {umap_shape})")

        # Check 2: Vector Field bases
        if compute_vector_field:
            if vf_bases:
                st.success(f"✅ Vector Field bases: {', '.join(vf_bases)}")
            else:
                st.error("❌ Vector Field is enabled but no basis selected")
                validation_passed = False

        # Check 3: File compatibility
        st.success(f"✅ Files uploaded: h5ad={uploaded_h5ad is not None}, loom={uploaded_loom is not None}")

        # Summary
        if validation_passed:
            st.markdown("---")
            st.success("✅ **All checks passed!** Configuration is valid.")

            # Store configuration in session state
            st.session_state.config_validated = True
            # Parse gene lists
            force_genes_list = parse_gene_list(force_genes) if not skip_preprocess or recompute_hvg else parse_gene_list(force_genes) if skip_preprocess and recompute_hvg else []
            force_genes_only_list = parse_gene_list(force_genes_only) if skip_preprocess and not recompute_hvg else []

            st.session_state.config = {
                'selected_basis': selected_basis,
                'selected_umap': selected_umap,
                'skip_preprocess': skip_preprocess,
                'preprocess_recipe': preprocess_recipe,
                'recompute_hvg': recompute_hvg if skip_preprocess else False,
                'n_top_genes': n_top_genes,
                'force_genes': force_genes_list,
                'force_genes_only': force_genes_only_list,
                'cell_cycle_score': cell_cycle_score,
                'dynamics_model': dynamics_model,
                'est_method': est_method,
                'n_jobs': n_jobs,
                'compute_vector_field': compute_vector_field,
                'vf_bases': vf_bases,
                'vf_M': vf_M,
                'compute_geometry': compute_geometry
            }

            # Show parsed genes if any
            if force_genes_list:
                st.info(f"✓ Force include genes: {len(force_genes_list)} genes - {force_genes_list[:5]}{'...' if len(force_genes_list) > 5 else ''}")
            if force_genes_only_list:
                st.info(f"✓ Add to existing HVGs: {len(force_genes_only_list)} genes - {force_genes_only_list[:5]}{'...' if len(force_genes_only_list) > 5 else ''}")

            st.info("👇 Click **Run Analysis** below to start the computation")
        else:
            st.error("❌ **Configuration check failed.** Please fix the issues above and check again.")
            st.session_state.config_validated = False

    # ========================================
    # Step 3: Run analysis
    # ========================================
    if st.session_state.get('config_validated', False):
        st.markdown("---")

        run_analysis = st.button("🌊 Run Analysis", type="primary", use_container_width=True)

        if run_analysis:
            st.header("Step 5: Running analysis")

            # Retrieve configuration
            config = st.session_state.config
            selected_basis = config['selected_basis']
            selected_umap = config['selected_umap']
            skip_preprocess = config['skip_preprocess']
            preprocess_recipe = config['preprocess_recipe']
            recompute_hvg = config.get('recompute_hvg', False)
            n_top_genes = config.get('n_top_genes', 2000)  # Default to 2000 for backward compatibility
            cell_cycle_score = config['cell_cycle_score']
            dynamics_model = config['dynamics_model']
            est_method = config['est_method']
            n_jobs = config['n_jobs']
            compute_vector_field = config['compute_vector_field']
            vf_bases = config['vf_bases']
            vf_M = config['vf_M']
            compute_geometry = config['compute_geometry']

            with st.spinner("Loading and processing files..."):
                progress_bar = st.progress(0)
                status_text = st.empty()

                try:
                    # Load h5ad if provided (use as base)
                    if uploaded_h5ad is not None:
                        status_text.text("Loading h5ad file...")
                        progress_bar.progress(10)

                        adata = load_h5ad_file(uploaded_h5ad, _get_file_hash(uploaded_h5ad))
                        st.info(f"✓ Loaded h5ad: {adata.n_obs} cells, {adata.n_vars} genes")

                        # Load loom
                        status_text.text("Loading loom file...")
                        progress_bar.progress(15)

                        ldata = load_loom_file(uploaded_loom, _get_file_hash(uploaded_loom))
                        st.info(f"✓ Loaded loom: {ldata.n_obs} cells, {ldata.n_vars} genes")

                        # Check for spliced/unspliced layers in loom
                        if 'spliced' not in ldata.layers or 'unspliced' not in ldata.layers:
                            st.error("""
                            ❌ **Missing spliced/unspliced layers in loom**

                            Dynamoの解析にはspliced/unspliced行列が必要です。
                            velocytoまたはData filtering appで生成されたloomファイルをアップロードしてください。
                            """)
                            st.stop()

                        # Check for ambiguous layer (optional but recommended for Dynamo)
                        available_layers = list(ldata.layers.keys())
                        st.info(f"✓ Available layers in loom: {', '.join(available_layers)}")

                        if 'ambiguous' in ldata.layers:
                            st.success("✓ Ambiguous reads layer found - Dynamoはより正確な推定が可能です")
                        else:
                            st.warning("⚠️ Ambiguous layer not found - 通常はvelocytoで生成されます（オプション）")

                        # Check cell overlap
                        status_text.text("Checking cell correspondence...")
                        progress_bar.progress(18)

                        h5ad_cells = set(adata.obs_names)
                        loom_cells = set(ldata.obs_names)
                        common_cells = h5ad_cells & loom_cells

                        if len(common_cells) == 0:
                            st.error(f"""
                            ❌ **No common cells found**

                            h5ad cells: {len(h5ad_cells)}
                            loom cells: {len(loom_cells)}

                            Cell IDが一致しません。同じデータセットから生成されたファイルをアップロードしてください。
                            """)
                            st.stop()

                        st.info(f"✓ Found {len(common_cells)} common cells")

                        # Filter h5ad to common cells
                        status_text.text("Filtering to common cells...")
                        progress_bar.progress(20)

                        adata = adata[list(common_cells), :].copy()
                        st.success(f"✓ Filtered h5ad to {adata.n_obs} cells")

                        # Merge using scvelo's merge function (adds spliced/unspliced from loom to h5ad)
                        status_text.text("Merging h5ad with loom layers...")
                        progress_bar.progress(22)

                        # Backup varm and uns before merge (scvelo merge may not preserve them)
                        # Also backup gene names for proper reindexing after merge
                        genes_before_merge = list(adata.var_names)
                        varm_backup = {k: v.copy() for k, v in adata.varm.items()} if len(adata.varm) > 0 else {}
                        uns_backup = {}
                        for k in list(adata.uns.keys()):
                            if any(kw in k.lower() for kw in ['pca', 'mnn', 'harmony', 'mean', 'gene', 'stdev']):
                                uns_backup[k] = adata.uns[k].copy() if hasattr(adata.uns[k], 'copy') else adata.uns[k]

                        if varm_backup:
                            st.info(f"📦 Backed up varm: {list(varm_backup.keys())} ({len(genes_before_merge)} genes)")

                        import scvelo as scv
                        adata = scv.utils.merge(adata, ldata)
                        st.success(f"✓ Merged data: {adata.n_obs} cells, {adata.n_vars} genes")

                        # Restore varm (PCA loadings) after merge
                        # Use gene names to properly reindex if genes were filtered
                        if varm_backup:
                            import pandas as pd
                            genes_after_merge = list(adata.var_names)
                            n_genes_after = len(genes_after_merge)
                            restored_varm = []

                            for key, mat in varm_backup.items():
                                n_genes_backup = mat.shape[0]
                                if n_genes_backup == n_genes_after and genes_before_merge == genes_after_merge:
                                    # Same genes in same order, can use directly
                                    adata.varm[key] = mat
                                    restored_varm.append(key)
                                else:
                                    # Genes were filtered or reordered - need to reindex
                                    st.info(f"⚠️ {key}: reindexing {n_genes_backup} → {n_genes_after} genes")

                                    # Create DataFrame with gene names as index for proper reindexing
                                    df_backup = pd.DataFrame(mat, index=genes_before_merge)

                                    # Filter to genes that exist after merge
                                    common_genes = [g for g in genes_after_merge if g in genes_before_merge]

                                    if len(common_genes) > 0:
                                        # Reindex to match current gene order
                                        mat_reindexed = df_backup.loc[common_genes].values

                                        # If some genes are missing, we need to fill with zeros
                                        if len(common_genes) < n_genes_after:
                                            # Create full matrix with zeros
                                            full_mat = np.zeros((n_genes_after, mat.shape[1]))
                                            # Fill in values for genes that have loadings
                                            gene_to_idx = {g: i for i, g in enumerate(genes_after_merge)}
                                            for i, g in enumerate(common_genes):
                                                full_mat[gene_to_idx[g], :] = mat_reindexed[i, :]
                                            adata.varm[key] = full_mat
                                        else:
                                            adata.varm[key] = mat_reindexed
                                        restored_varm.append(key)
                                    else:
                                        st.warning(f"⚠️ Could not restore {key}: no common genes")

                            if restored_varm:
                                st.success(f"✓ Restored varm: {restored_varm}")

                        # Restore uns (PCA means, etc.) after merge
                        if uns_backup:
                            for key, val in uns_backup.items():
                                adata.uns[key] = val
                            st.success(f"✓ Restored uns: {list(uns_backup.keys())}")

                        # Clean up old scVelo metadata to avoid conflicts with Dynamo
                        old_metadata_keys = []
                        if 'recover_dynamics' in adata.uns:
                            old_metadata_keys.append('recover_dynamics')
                        if 'velocity_params' in adata.uns:
                            old_metadata_keys.append('velocity_params')

                        for key in old_metadata_keys:
                            del adata.uns[key]

                        if old_metadata_keys:
                            st.info(f"✓ Removed old scVelo metadata: {old_metadata_keys}")

                    else:
                        # Load only loom (no h5ad provided)
                        status_text.text("Loading loom file...")
                        progress_bar.progress(10)

                        adata = load_loom_file(uploaded_loom, _get_file_hash(uploaded_loom))
                        st.info(f"✓ Loaded loom: {adata.n_obs} cells, {adata.n_vars} genes")

                        # Check for spliced/unspliced layers
                        if 'spliced' not in adata.layers or 'unspliced' not in adata.layers:
                            st.error("""
                            ❌ **Missing spliced/unspliced layers**

                            Dynamoの解析にはspliced/unspliced行列が必要です。
                            velocytoまたはData filtering appで生成されたloomファイルをアップロードしてください。
                            """)
                            st.stop()

                        # Check for ambiguous layer (optional but recommended for Dynamo)
                        available_layers = list(adata.layers.keys())
                        st.info(f"✓ Available layers in loom: {', '.join(available_layers)}")

                        if 'ambiguous' in adata.layers:
                            st.success("✓ Ambiguous reads layer found - Dynamoはより正確な推定が可能です")
                        else:
                            st.warning("⚠️ Ambiguous layer not found - 通常はvelocytoで生成されます（オプション）")

                    st.success(f"✓ Found spliced/unspliced layers")

                    # Preprocessing
                    if not skip_preprocess:
                        status_text.text("Preprocessing with Dynamo...")
                        progress_bar.progress(25)

                        preprocessor = dyn.pp.Preprocessor(
                            cell_cycle_score_enable=False,  # Will compute separately if needed
                            pca_kwargs={'n_pca_components': n_pca_components}
                        )
                        preprocessor.preprocess_adata(
                            adata,
                            recipe=preprocess_recipe,
                            n_top_genes=n_top_genes
                        )
                        st.info(f"✓ PCA computed with {n_pca_components} components")
                        st.success(f"✓ Preprocessing complete (recipe={preprocess_recipe})")
                    else:
                        # Skip preprocessing - but run full preprocessing then restore Seurat data
                        status_text.text("Backing up Seurat-processed data...")
                        progress_bar.progress(20)

                        # Backup Seurat-processed data (X matrix and all embeddings)
                        X_backup = adata.X.copy()
                        obsm_backup = {}
                        for key in adata.obsm.keys():
                            obsm_backup[key] = adata.obsm[key].copy()

                        st.info(f"✓ Backed up X matrix and {len(obsm_backup)} embeddings")

                        # Run full preprocessing (for spliced/unspliced layers)
                        status_text.text("Preprocessing with Dynamo (spliced/unspliced layers)...")
                        progress_bar.progress(22)

                        # Save cell names before preprocessing
                        cells_before = set(adata.obs_names)
                        n_cells_before = len(cells_before)

                        preprocessor = dyn.pp.Preprocessor(
                            cell_cycle_score_enable=False,
                            pca_kwargs={'n_pca_components': n_pca_components}
                        )
                        preprocessor.preprocess_adata(
                            adata,
                            recipe=preprocess_recipe if preprocess_recipe else "monocle"
                        )
                        st.info(f"✓ PCA computed with {n_pca_components} components (for spliced/unspliced layers)")

                        # Check which cells were filtered
                        cells_after = set(adata.obs_names)
                        n_cells_after = len(cells_after)
                        filtered_cells = cells_before - cells_after

                        st.info(f"✓ Preprocessing complete (spliced/unspliced fully processed)")
                        st.info(f"  Cells before: {n_cells_before}, after: {n_cells_after}")

                        # Reselect HVGs if n_top_genes differs from default (2000)
                        if n_top_genes != 2000:
                            status_text.text(f"Reselecting HVGs ({n_top_genes} genes)...")
                            import scanpy as sc
                            sc.pp.highly_variable_genes(
                                adata,
                                n_top_genes=n_top_genes,
                                flavor='seurat_v3'
                            )
                            # Filter to HVGs
                            adata = adata[:, adata.var['highly_variable']].copy()
                            st.info(f"  ✓ Reselected HVGs: {adata.n_vars} genes")

                        if len(filtered_cells) > 0:
                            st.warning(f"  {len(filtered_cells)} cells were filtered during preprocessing")
                            st.info(f"  Filtered cells: {list(filtered_cells)[:10]}...")  # Show first 10

                        # Restore Seurat-processed X matrix and all embeddings
                        status_text.text("Restoring Seurat-processed data...")
                        progress_bar.progress(25)

                        # Check if cell count changed
                        if X_backup.shape[0] != adata.shape[0]:
                            st.warning(f"Cell count mismatch: backup has {X_backup.shape[0]}, adata has {adata.shape[0]}")
                            st.info("Filtering backup data to match surviving cells...")

                            # Get the cells that survived preprocessing
                            surviving_cells = adata.obs_names

                            # Filter X_backup to only include surviving cells
                            # Assuming X_backup is a numpy array/sparse matrix
                            # We need to find the indices of surviving cells in the original data
                            original_cells = list(cells_before)
                            surviving_indices = [i for i, cell in enumerate(original_cells) if cell in surviving_cells]

                            st.info(f"  Surviving cells: {len(surviving_cells)}/{len(original_cells)}")

                            # Filter X_backup
                            X_backup_filtered = X_backup[surviving_indices, :]
                            adata.X = X_backup_filtered

                            # Filter obsm_backup
                            for key, value in obsm_backup.items():
                                obsm_backup[key] = value[surviving_indices, :]
                                adata.obsm[key] = obsm_backup[key]

                            st.success("✓ Velocity layers (spliced/unspliced) fully preprocessed")
                            st.info("✓ Seurat X matrix and embeddings filtered to match surviving cells")
                        else:
                            # No cell filtering occurred
                            adata.X = X_backup
                            for key, value in obsm_backup.items():
                                adata.obsm[key] = value

                            st.success("✓ Velocity layers (spliced/unspliced) fully preprocessed")
                            st.info("✓ Seurat X matrix and all embeddings preserved")

                        # Recompute HVGs if requested (using existing normalized data)
                        if recompute_hvg:
                            status_text.text(f"Recomputing HVGs ({n_top_genes} genes)...")
                            progress_bar.progress(27)

                            import scanpy as sc
                            # Use scanpy to compute HVGs on the normalized data (adata.X)
                            sc.pp.highly_variable_genes(
                                adata,
                                n_top_genes=n_top_genes,
                                flavor='seurat_v3',
                                layer=None  # Use adata.X (normalized data)
                            )

                            # Filter to HVGs
                            adata = adata[:, adata.var['highly_variable']].copy()

                            st.success(f"✓ Recomputed HVGs: {adata.n_vars} genes selected from {n_top_genes} requested")
                            st.info("  Using existing normalized data (Seurat X matrix preserved)")

                    # Compute cell cycle score if requested (independent of preprocessing)
                    if cell_cycle_score:
                        status_text.text("Computing cell cycle score...")
                        try:
                            from dynamo.preprocessing.cell_cycle import cell_cycle_scores
                            cell_cycle_scores(adata)
                            st.success("✓ Cell cycle score computed")
                        except Exception as e:
                            st.warning(f"⚠️ Cell cycle score computation skipped: {str(e)}\n(This is optional - analysis will continue)")

                    # Dynamics estimation
                    status_text.text(f"Estimating dynamics ({dynamics_model} model)...")
                    progress_bar.progress(40)

                    st.warning(f"⏳ Dynamics estimation may take several minutes... (using {n_jobs} cores)")

                    if dynamics_model == "stochastic":
                        dyn.tl.dynamics(
                            adata,
                            model='stochastic',
                            est_method=est_method,
                            cores=n_jobs
                        )
                    else:
                        dyn.tl.dynamics(
                            adata,
                            model=dynamics_model,
                            cores=n_jobs
                        )

                    st.success(f"✓ Dynamics estimated ({dynamics_model} model)")

                    # Copy PCA loadings from varm to uns for Dynamo perturbation compatibility
                    status_text.text("Preparing PCA loadings for Dynamo perturbation...")
                    progress_bar.progress(58)

                    if uploaded_h5ad is not None and selected_basis:
                        # Check if PCA loadings exist in varm
                        varm_key = f'PCs_{selected_basis}'

                        if varm_key in adata.varm.keys():
                            # Copy to uns (Dynamo perturbation expects loadings in uns, not varm)
                            adata.uns['PCs'] = adata.varm[varm_key].copy()
                            adata.uns[varm_key] = adata.varm[varm_key].copy()

                            st.success(f"✓ Copied PCA loadings: varm['{varm_key}'] → uns['PCs'] (Dynamo perturbation compatible)")
                            st.info(f"  Shape: {adata.varm[varm_key].shape} (genes × PCs)")

                            # Set use_for_pca flag (required by Dynamo perturbation)
                            # Mark genes that have non-zero loadings as used for PCA
                            loadings = adata.varm[varm_key]
                            pca_genes_mask = (np.abs(loadings).sum(axis=1) > 0)
                            adata.var['use_for_pca'] = pca_genes_mask

                            n_pca_genes = pca_genes_mask.sum()
                            st.info(f"✓ Set use_for_pca flag: {n_pca_genes} / {adata.n_vars} genes used in PCA")

                            # Also copy pca_mean if it exists in uns (from seurat2ann)
                            pca_mean_key = f'{selected_basis}_mean'
                            if pca_mean_key in adata.uns.keys():
                                adata.uns['pca_mean'] = adata.uns[pca_mean_key].copy()
                                st.info(f"✓ Copied PCA mean: uns['{pca_mean_key}'] → uns['pca_mean']")
                        else:
                            st.warning(f"⚠️ PCA loadings not found in varm['{varm_key}']")
                            st.info("  Perturbation analysis may not work without PCA loadings")
                            st.info("  Ensure h5ad was exported with 'for_cellxgene = FALSE' from SCALA")

                    # Dimensionality reduction - preserve selected_umap and prepare for X_dynamo.umap
                    status_text.text("Managing UMAP embeddings...")
                    progress_bar.progress(60)

                    if uploaded_h5ad is not None and selected_umap:
                        # h5ad uploaded with selected_umap - preserve it and copy to X_umap
                        selected_umap_key = f'X_{selected_umap}'

                        if selected_umap_key in adata.obsm:
                            # selected_umap is already preserved with its original name
                            st.info(f"✓ Selected UMAP '{selected_umap}' preserved: {selected_umap_key}")

                            # Copy to X_umap for Dynamo compatibility (required for Perturbation analysis)
                            adata.obsm['X_umap'] = adata.obsm[selected_umap_key].copy()
                            st.info(f"  Copied {selected_umap_key} → X_umap (for Dynamo/Perturbation compatibility)")

                            # New UMAP (X_dynamo.umap) will be computed after neighbors
                            st.info("  New UMAP (X_dynamo.umap) will be computed after neighbors (KNN consistent)")
                        else:
                            st.warning(f"⚠️ Selected UMAP {selected_umap_key} not found in adata.obsm")
                            st.info("  Will create X_umap and X_dynamo.umap from scratch")
                    else:
                        # loom-only: compute UMAP if needed (after neighbors)
                        st.info("✓ UMAP (X_dynamo.umap) will be computed after neighbors are established")

                    # Neighbors computation - check existing and compute with correct use_rep
                    status_text.text("Checking/computing neighbors...")
                    progress_bar.progress(65)

                    # Check if neighbors exist
                    has_connectivities = 'connectivities' in adata.obsp
                    has_distances = 'distances' in adata.obsp

                    if has_connectivities and has_distances:
                        # Neighbors exist - preserve them
                        st.success("✓ Using existing neighbors (connectivities & distances preserved)")
                        if 'neighbors' in adata.uns:
                            existing_use_rep = adata.uns['neighbors']['params'].get('use_rep', 'unknown')
                            existing_n_neighbors = adata.uns['neighbors']['params'].get('n_neighbors', 'unknown')
                            st.info(f"  Existing neighbors: use_rep='{existing_use_rep}', n_neighbors={existing_n_neighbors}")
                    else:
                        # Compute neighbors with correct use_rep
                        st.warning("⚠️ Neighbors not found - computing with selected basis")

                        # Determine use_rep from selected_basis
                        if selected_basis:
                            use_rep_key = f'X_{selected_basis}'
                        else:
                            # Fallback: use PCA if available
                            if 'X_pca' in adata.obsm:
                                use_rep_key = 'X_pca'
                            else:
                                use_rep_key = None

                        if use_rep_key and use_rep_key in adata.obsm:
                            st.info(f"  Computing neighbors: use_rep='{use_rep_key}', n_neighbors=30, metric='euclidean'")

                            try:
                                sc.pp.neighbors(
                                    adata,
                                    use_rep=use_rep_key,
                                    n_neighbors=30,
                                    metric='euclidean',
                                    method='umap',
                                    random_state=0
                                )
                                st.success(f"✓ Neighbors computed using {use_rep_key}")

                                # Update flags after computation
                                has_connectivities = True
                                has_distances = True
                            except Exception as e:
                                st.error(f"❌ Failed to compute neighbors: {str(e)}")
                                has_connectivities = False
                                has_distances = False
                        else:
                            st.error(f"❌ Cannot compute neighbors - {use_rep_key} not found in adata.obsm")
                            has_connectivities = False
                            has_distances = False

                    # Compute new UMAP (X_dynamo.umap) consistent with KNN
                    status_text.text("Computing new UMAP (X_dynamo.umap)...")
                    progress_bar.progress(68)

                    if has_connectivities and has_distances:
                        try:
                            # Step 1: Backup existing X_umap if it exists (to preserve original)
                            existing_x_umap = None
                            if 'X_umap' in adata.obsm:
                                existing_x_umap = adata.obsm['X_umap'].copy()

                            # Step 2: Compute new UMAP using the neighbors graph
                            # This will create/overwrite X_umap
                            sc.tl.umap(adata, random_state=0)

                            # Step 3: Save newly computed UMAP as X_dynamo.umap
                            if 'X_umap' in adata.obsm:
                                adata.obsm['X_dynamo.umap'] = adata.obsm['X_umap'].copy()
                                st.success("✓ New UMAP computed → X_dynamo.umap (consistent with KNN)")

                                # Step 4: Restore original X_umap if it existed
                                if existing_x_umap is not None:
                                    adata.obsm['X_umap'] = existing_x_umap
                                    st.info("  Original X_umap restored")
                                else:
                                    st.info("  X_umap created (new UMAP also saved as X_dynamo.umap)")
                            else:
                                st.warning("⚠️ UMAP computation did not create X_umap")
                        except Exception as e:
                            st.warning(f"⚠️ Could not compute new UMAP: {str(e)}")
                    else:
                        st.warning("⚠️ Cannot compute UMAP - neighbors not available")

                    # Cell velocities - only compute general velocities if NOT using per-basis Vector Field
                    if not (compute_vector_field and vf_bases):
                        status_text.text("Computing cell velocities...")
                        progress_bar.progress(70)

                        dyn.tl.cell_velocities(
                            adata,
                            method='pearson',
                            other_kernels_dict={'transform': 'sqrt'}
                        )
                        st.success("✓ Cell velocities computed")

                        # Cell-wise confidence
                        status_text.text("Computing cell-wise confidence...")
                        progress_bar.progress(75)

                        dyn.tl.cell_wise_confidence(adata)
                        st.success("✓ Cell-wise confidence computed")
                    else:
                        st.info("✓ Skipping general cell velocities - will compute per-basis in Vector Field loop")
                        progress_bar.progress(75)

                    # Vector field reconstruction
                    if compute_vector_field and vf_bases:
                        st.warning(f"⏳ Vector field reconstruction for {len(vf_bases)} bases may take several minutes...")

                        # Neighbors should already be computed above - just verify
                        has_connectivities = 'connectivities' in adata.obsp
                        has_distances = 'distances' in adata.obsp
                        neighbors_available = has_connectivities and has_distances

                        if not neighbors_available:
                            st.error("❌ Neighbors not available - cannot proceed with Vector Field reconstruction")
                            st.warning("This should not happen - neighbors should have been computed earlier")
                        else:
                            st.info("✓ Neighbors available for Vector Field reconstruction")

                        for i, vf_basis in enumerate(vf_bases):
                            # Check if embedding exists
                            status_text.text(f"Checking embedding for {vf_basis} [{i+1}/{len(vf_bases)}]...")
                            progress_bar.progress(80 + (i * 5 // len(vf_bases)))  # Progress 80-85

                            try:
                                embedding_key = f'X_{vf_basis}'
                                if embedding_key not in adata.obsm:
                                    st.warning(f"⚠️ Embedding {embedding_key} not found, skipping {vf_basis}")
                                    continue

                                # Neighbors already computed - proceed with cell velocities
                                status_text.text(f"Computing cell velocities for {vf_basis} [{i+1}/{len(vf_bases)}]...")

                                dyn.tl.cell_velocities(
                                    adata,
                                    basis=vf_basis,
                                    method='pearson',
                                    other_kernels_dict={'transform': 'sqrt'}
                                )
                                st.info(f"✓ Cell velocities computed for {vf_basis}")
                            except Exception as e:
                                st.warning(f"⚠️ Could not compute cell velocities for {vf_basis}: {str(e)}")
                                continue  # Skip this basis if velocity computation fails

                            # Now compute vector field
                            status_text.text(f"Reconstructing vector field ({vf_basis}) [{i+1}/{len(vf_bases)}]...")

                            try:
                                dyn.vf.VectorField(
                                    adata,
                                    basis=vf_basis,
                                    M=vf_M,
                                    pot_curl_div=True
                                )
                                st.success(f"✓ Vector field reconstructed on {vf_basis}")

                                # Debug: Check what Vector Field keys were actually created
                                expected_key = f'VecFld_{vf_basis}'
                                if expected_key in adata.uns:
                                    st.info(f"🔍 Debug: Vector Field saved as **{expected_key}**")
                                else:
                                    # List all VecFld keys to see what was actually created
                                    vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
                                    st.warning(f"⚠️ Debug: Expected key '{expected_key}' not found. Available: {vecfld_keys}")

                                # Save basis-specific PCA loadings for perturbation analysis
                                # This ensures PCs dimensions match the embedding used for VectorField
                                # Skip low-dimensional spaces (UMAP, tSNE, etc.) - they're not useful for perturbation
                                low_dim_keywords = ['umap', 'tsne', 'phate', 'trimap', 'draw_graph', 'diffmap']
                                is_low_dim_basis = any(kw in vf_basis.lower() for kw in low_dim_keywords)

                                embedding_key = f'X_{vf_basis}'
                                if embedding_key in adata.obsm and not is_low_dim_basis:
                                    n_dims = adata.obsm[embedding_key].shape[1]
                                    pcs_key = f'PCs_{vf_basis}'
                                    mean_key = f'{vf_basis}_mean'

                                    # Check if basis-specific PCs already exist with correct dimensions
                                    # First check varm (Seurat conversion stores loadings there)
                                    # Then check uns (dynamo stores loadings there)
                                    pcs_exists = False
                                    if pcs_key in adata.varm and adata.varm[pcs_key].shape[1] == n_dims:
                                        st.success(f"✓ Found existing PCA loadings in varm: {pcs_key} ({adata.varm[pcs_key].shape})")
                                        pcs_exists = True
                                    elif pcs_key in adata.uns and adata.uns[pcs_key].shape[1] == n_dims:
                                        st.success(f"✓ Found existing PCA loadings in uns: {pcs_key} ({adata.uns[pcs_key].shape})")
                                        pcs_exists = True

                                    if not pcs_exists:
                                        st.info(f"📊 Computing PCA loadings for {vf_basis} ({n_dims} components)...")
                                        try:
                                            from sklearn.decomposition import PCA as sklearn_PCA

                                            # Get HVG genes
                                            hvg_mask = adata.var.use_for_pca.values
                                            n_hvg = hvg_mask.sum()

                                            # Get expression data for HVG genes
                                            if hasattr(adata.X, 'toarray'):
                                                X_hvg = adata.X[:, hvg_mask].toarray()
                                            else:
                                                X_hvg = adata.X[:, hvg_mask]

                                            # Fit PCA with matching dimensions
                                            n_components = min(n_dims, n_hvg - 1)
                                            pca_model = sklearn_PCA(n_components=n_components, random_state=0)
                                            pca_model.fit(X_hvg)

                                            # Save loadings and mean
                                            adata.uns[pcs_key] = pca_model.components_.T  # (n_genes, n_components)
                                            adata.uns[mean_key] = pca_model.mean_

                                            st.success(f"✓ Saved {pcs_key}: {adata.uns[pcs_key].shape}, {mean_key}: {adata.uns[mean_key].shape}")
                                        except Exception as e:
                                            st.warning(f"⚠️ Could not compute PCA loadings for {vf_basis}: {str(e)}")
                            except Exception as e:
                                st.error(f"❌ Vector field reconstruction failed for {vf_basis}: {str(e)}")

                        progress_bar.progress(85)

                        # Compute cell-wise confidence after Vector Field reconstruction
                        # (only if it wasn't computed earlier)
                        if 'velocity_confidence' not in adata.obs:
                            status_text.text("Computing cell-wise confidence...")
                            try:
                                dyn.tl.cell_wise_confidence(adata)
                                st.success("✓ Cell-wise confidence computed")
                            except Exception as e:
                                st.warning(f"⚠️ Could not compute cell-wise confidence: {str(e)}")

                    # Differential geometry
                    if compute_geometry and compute_vector_field and vf_bases:
                        st.warning(f"⏳ Computing geometric features for {len(vf_bases)} bases may take several minutes...")

                        for j, geom_basis in enumerate(vf_bases):
                            status_text.text(f"Computing geometric features for {geom_basis} [{j+1}/{len(vf_bases)}]...")
                            progress_bar.progress(85 + (j * 10 // len(vf_bases)))  # Progress 85-95

                            # Vector Field should already exist (computed in previous step)
                            # But double-check
                            if f'VecFld_{geom_basis}' not in adata.uns:
                                st.warning(f"⚠️ Vector Field not found for {geom_basis}, skipping geometry")
                                continue

                            # Speed
                            try:
                                dyn.vf.speed(adata, basis=geom_basis)
                                st.info(f"✓ Computed speed on {geom_basis}")
                            except Exception as e:
                                st.warning(f"Could not compute speed for {geom_basis}: {str(e)}")

                            # Curl (only for 2D embeddings)
                            try:
                                # Check if geom_basis is 2D
                                embedding_key = f'X_{geom_basis}'
                                if embedding_key in adata.obsm and adata.obsm[embedding_key].shape[1] == 2:
                                    dyn.vf.curl(adata, basis=geom_basis)
                                    st.info(f"✓ Computed curl on {geom_basis}")
                                else:
                                    st.info(f"⚠️ Curl skipped for {geom_basis} (requires 2D embedding)")
                            except Exception as e:
                                st.warning(f"Could not compute curl for {geom_basis}: {str(e)}")

                            # Divergence
                            try:
                                dyn.vf.divergence(adata, basis=geom_basis)
                                st.info(f"✓ Computed divergence on {geom_basis}")
                            except Exception as e:
                                st.warning(f"Could not compute divergence for {geom_basis}: {str(e)}")

                            # Acceleration
                            try:
                                dyn.vf.acceleration(adata, basis=geom_basis)
                                st.info(f"✓ Computed acceleration on {geom_basis}")
                            except Exception as e:
                                st.warning(f"Could not compute acceleration for {geom_basis}: {str(e)}")

                            # Curvature
                            try:
                                dyn.vf.curvature(adata, basis=geom_basis)
                                st.info(f"✓ Computed curvature on {geom_basis}")
                            except Exception as e:
                                st.warning(f"Could not compute curvature for {geom_basis}: {str(e)}")

                        st.success(f"✓ Geometric features computed for {len(vf_bases)} bases")

                    # Save result
                    status_text.text("Saving results...")
                    progress_bar.progress(95)

                    output_path = os.path.join(dynamo_temp_dir, "dynamo_result.h5ad")
                    if os.path.exists(output_path):
                        os.remove(output_path)

                    adata.write_h5ad(output_path, compression="gzip")

                    # Store in session state
                    st.session_state.dynamo_result_path = output_path
                    st.session_state.dynamo_adata = adata
                    st.session_state.dynamo_model = dynamics_model
                    st.session_state.dynamo_complete = True

                    progress_bar.progress(100)
                    status_text.text("Analysis complete!")

                    st.success("""
                    ✅ **Analysis completed successfully!**

                    Computed results:
                    - RNA dynamics
                    - Cell velocities
                    - Vector field (if enabled)
                    - Geometric features (if enabled)
                    """)

                    # Show summary statistics
                    st.subheader("Summary Statistics")

                    summary_col1, summary_col2, summary_col3 = st.columns(3)

                    with summary_col1:
                        st.metric("Total cells", adata.n_obs)
                        st.metric("Total genes", adata.n_vars)

                    with summary_col2:
                        if 'velocity_confidence' in adata.obs:
                            mean_confidence = adata.obs['velocity_confidence'].mean()
                            st.metric("Mean velocity confidence", f"{mean_confidence:.4f}")

                        if 'speed' in adata.obs:
                            mean_speed = adata.obs['speed'].mean()
                            st.metric("Mean speed", f"{mean_speed:.4f}")

                    with summary_col3:
                        if compute_vector_field and vf_bases:
                            st.metric("Vector field bases", f"{len(vf_bases)} bases")
                            st.write(f"**Bases**: {', '.join(vf_bases)}")
                            st.metric("Sampling points", vf_M)

                    # Show available data
                    with st.expander("📊 Available data in result", expanded=False):
                        st.write("**Layers:**")
                        st.write(list(adata.layers.keys()))

                        st.write("**Obs (cell metadata):**")
                        velocity_cols = [col for col in adata.obs.columns if 'velocity' in col or 'speed' in col or 'curl' in col or 'div' in col or 'accel' in col or 'curv' in col]
                        st.write(velocity_cols)

                        st.write("**Obsm (embeddings):**")
                        st.write(list(adata.obsm.keys()))

                        st.write("**Uns (Vector Field storage):**")
                        # Show all Vector Field related keys
                        vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
                        if vecfld_keys:
                            st.write(f"✓ Vector Field keys found: **{', '.join(vecfld_keys)}**")
                            # Also show which embeddings these correspond to
                            vecfld_bases = [k.replace('VecFld_', '') for k in vecfld_keys]
                            st.write(f"   Bases: {', '.join(vecfld_bases)}")

                            # Check if embeddings exist for each Vector Field
                            for vf_base in vecfld_bases:
                                emb_key = f'X_{vf_base}'
                                if emb_key in adata.obsm:
                                    emb_shape = adata.obsm[emb_key].shape
                                    st.write(f"   ✓ {emb_key}: shape={emb_shape}")
                                else:
                                    st.write(f"   ❌ {emb_key}: NOT FOUND (mismatch!)")
                        else:
                            st.write("❌ No Vector Field data found")

                except Exception as e:
                    st.error(f"❌ Error during analysis: {str(e)}")
                    st.exception(e)
                    st.session_state.dynamo_complete = False

                finally:
                    progress_bar.empty()
                    status_text.empty()

    # ========================================
    # Step 4: Download results
    # ========================================
    if st.session_state.dynamo_complete:
        st.header("Step 6: Download results")

        with open(st.session_state.dynamo_result_path, "rb") as f:
            result_bytes = f.read()

        # Create filename
        if uploaded_h5ad is not None:
            h5ad_basename = os.path.splitext(uploaded_h5ad.name)[0]
        else:
            loom_basename = os.path.splitext(uploaded_loom.name)[0]
            h5ad_basename = loom_basename

        dynamics_model = st.session_state.get('dynamo_model', 'stochastic')
        output_filename = f"{h5ad_basename}.Dynamo.{dynamics_model}.h5ad"

        st.download_button(
            label="⬇️ Download Dynamo result (h5ad)",
            data=result_bytes,
            file_name=output_filename,
            mime="application/octet-stream",
            type="primary"
        )

        st.info("""
        ### 次のステップ

        ダウンロードしたh5adファイルには以下が含まれています：
        - RNA dynamics parameters
        - Velocity vectors (`adata.layers['velocity']`)
        - Vector field (`adata.uns['VecFld']`)（計算した場合）
        - Geometric features: speed, curl, divergence, acceleration, curvature

        **📌 UMAP embeddingsについて（重要）:**
        - **元のUMAP** (例: `X_rna.mnn.umap`): Seurat/元の解析で作成されたもの
        - **dynamo.umap** (`X_dynamo.umap`): 今回新しく計算されたKNN graphから作成

        ⚠️ **注意**: 元のUMAPでVector Fieldを可視化する場合、UMAP座標と速度ベクトルの整合性が若干ズレる可能性があります。
        これは、元のUMAPが元のneighbor graphに基づいており、今回計算された速度ベクトルは新しいneighbor graphに基づいているためです。

        💡 **推奨**: KNN graphと速度が完全に整合する`dynamo.umap`を使用すると、より正確な可視化が可能です。

        #### 🌟 推奨: Streamlitアプリで可視化

        **Dynamo Visualizationアプリ**を使ってインタラクティブに可視化：
        1. ダウンロードしたh5adファイルをアップロード
        2. Streamline plot、Cell-wise vectors、Phase portraits等を簡単に作成
        3. パラメータを調整してリアルタイムで結果を確認

        その後、**Dynamo Perturbationアプリ**で高度な解析：
        - In silico perturbation（遺伝子摂動シミュレーション）
        - Least action paths（最適遷移経路）
        - Jacobian analysis（制御ネットワーク推論）

        #### Pythonでの可視化例:
        ```python
        import dynamo as dyn
        import scanpy as sc

        # Load result
        adata = sc.read_h5ad('result.Dynamo.h5ad')

        # Visualizations
        # Streamline plot (vector field)
        dyn.pl.streamline_plot(adata, color=['clusters'], basis='umap')

        # Cell-wise vectors
        dyn.pl.cell_wise_vectors(adata, color=['clusters'], basis='umap')

        # Phase portraits
        dyn.pl.phase_portraits(adata, genes=['gene1', 'gene2'])

        # Geometric features
        dyn.pl.cell_wise_vectors(adata, color='speed', basis='umap')
        dyn.pl.cell_wise_vectors(adata, color='divergence', basis='umap')
        dyn.pl.cell_wise_vectors(adata, color='acceleration', basis='umap')

        # Kinetic parameters
        dyn.pl.kinetic_heatmap(adata, genes=['gene1', 'gene2'])
        ```

        #### 高度な解析 (Pythonで実行):
        ```python
        # Least action paths (細胞状態間の最適経路)
        dyn.pd.least_action(adata, [start_cells], [end_cells])

        # In silico perturbation (遺伝子摂動予測)
        dyn.pd.perturbation(adata, genes=['gene1'], expression=[0])

        # Regulatory network inference
        dyn.vf.jacobian(adata, regulators=['TF1', 'TF2'], effectors=['gene1', 'gene2'])
        ```

        #### その他のDynamo機能:
        - **Vector field analysis**: 軌跡予測、固定点解析
        - **Differential geometry**: Jacobian matrix、Hessian解析
        - **Cell fate prediction**: Attractor state identification
        - **Perturbation simulation**: 遺伝子ノックアウト/過剰発現シミュレーション

        詳細は [Dynamo Documentation](https://dynamo-release.readthedocs.io/) を参照してください。
        """)

else:
    st.info("👆 loomファイル（spliced/unspliced matrices必須）をアップロードして開始してください。h5adファイルは任意（メタデータ用）です。")
