"""
Dynamo Perturbation Analysis
Perform in silico perturbations and least action path analysis using Dynamo
"""

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import io
import tempfile
import time
from helper_func import clear_old_directories, clear_old_files

# Import dynamo
try:
    import dynamo as dyn
    DYNAMO_AVAILABLE = True
except ImportError:
    DYNAMO_AVAILABLE = False

st.set_page_config(page_title="Dynamo Prediction and Network", page_icon="🧬", layout="wide")

st.title("🧬 Dynamo Prediction and Network")

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
Dynamoを用いて遺伝子摂動シミュレーションと細胞運命予測を実行します。

### ワークフロー
1. **ファイル読み込み**: Dynamo解析済みh5adファイル（vector field計算済み）
2. **摂動シミュレーション**: 遺伝子のノックアウト/過剰発現の効果を予測
3. **最適経路解析**: 細胞状態間の最適遷移経路（Least Action Path）を計算
4. **制御ネットワーク**: Jacobian解析による遺伝子制御ネットワークの推論

### Dynamoの摂動解析機能
- **In silico perturbation**: ベクトル場を用いた遺伝子摂動の効果予測
- **Least action paths**: 細胞状態間の最も確率の高い遷移経路
- **Jacobian analysis**: 遺伝子間の制御関係の定量化
- **Fate prediction**: 細胞運命の予測と可視化

### 参考
- [Qiu et al. (2022) "Mapping transcriptomic vector fields of single cells" Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4)
- [Dynamo Documentation - Perturbation](https://dynamo-release.readthedocs.io/en/latest/API.html#module-dynamo.prediction)
""")

@st.cache_data
def load_h5ad_file(file_bytes, filename):
    """
    Load h5ad file with caching
    """
    import tempfile
    with tempfile.NamedTemporaryFile(delete=False, suffix=".h5ad") as tmp:
        tmp.write(file_bytes)
        tmp_path = tmp.name

    adata = sc.read_h5ad(tmp_path)
    os.unlink(tmp_path)  # Clean up temp file

    return adata

# Initialize session state
if "dynamo_pert_temp_dir" not in st.session_state:
    dynamo_pert_temp_dir = os.path.join("temp", f"dynamo_pert_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(dynamo_pert_temp_dir, exist_ok=True)
    st.session_state.dynamo_pert_temp_dir = dynamo_pert_temp_dir
else:
    dynamo_pert_temp_dir = st.session_state.dynamo_pert_temp_dir

if "dynamo_pert_complete" not in st.session_state:
    st.session_state.dynamo_pert_complete = False

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload Dynamo result")

st.markdown("""
### 必要なファイル

**Dynamo解析済みh5adファイル**をアップロード:
- ✅ Vector field計算済み (`adata.uns['VecFld_*']` - basis ごとに保存)
  - 例: `VecFld_rna.pca` (PCA basis)
  - 例: `VecFld_rna.mnn.umap` (UMAP basis)
- ✅ RNA dynamics推定済み
- ✅ Cell velocities計算済み
- ✅ Embedding (UMAP, PCA等)

→ Dynamo Analysisアプリで生成されたh5adファイルを使用してください

**⚠️ 重要**: Perturbation解析には**PCA basis**のVector Fieldが推奨されます（遺伝子発現空間での定量解析のため）
""")

uploaded_h5ad = st.file_uploader(
    "Upload Dynamo result (h5ad)",
    type=['h5ad'],
    key="dynamo_pert_h5ad_upload",
    help="Dynamo Analysisアプリで生成されたh5adファイル"
)

if uploaded_h5ad is not None:
    # Load file with caching
    with st.spinner("Loading file..."):
        adata = load_h5ad_file(uploaded_h5ad.getvalue(), uploaded_h5ad.name)
        st.success(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

        # Check required data
        st.subheader("Data validation")

        col1, col2, col3 = st.columns(3)

        with col1:
            # Check for Vector Field keys (can be VecFld or VecFld_<basis>)
            vecfld_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
            if len(vecfld_keys) > 0:
                st.success("✓ Vector field available")
                # Extract basis names
                vf_bases = [k.replace('VecFld_', '') if k.startswith('VecFld_') else 'default' for k in vecfld_keys]
                st.info(f"Bases: {', '.join(vf_bases)}")
            else:
                st.error("❌ Vector field not found")
                st.warning("Please run Dynamo Analysis first")

        with col2:
            # Check for velocity layers (can be velocity or velocity_<basis>)
            velocity_keys = [k for k in adata.layers.keys() if 'velocity' in k.lower()]
            if len(velocity_keys) > 0:
                st.success("✓ Velocity available")
                st.info(f"Keys: {', '.join(velocity_keys)}")
            else:
                st.error("❌ Velocity not found")

        with col3:
            umap_keys = [k for k in adata.obsm.keys() if 'umap' in k.lower()]
            if len(umap_keys) > 0:
                st.success(f"✓ UMAP available")
                st.info(f"Keys: {', '.join(umap_keys)}")
            else:
                st.error("❌ UMAP not found")

        # Check if analysis can proceed
        can_proceed = (len(vecfld_keys) > 0 and
                      len(velocity_keys) > 0 and
                      len(umap_keys) > 0)

        if not can_proceed:
            st.error("""
            ❌ **Cannot proceed with perturbation analysis**

            このファイルにはDynamo解析の必須データが不足しています。
            Dynamo Analysisアプリで先にvector field計算を実行してください。
            """)
            st.stop()

        # Get categorical columns for cell selection
        categorical_cols = [col for col in adata.obs.columns
                           if adata.obs[col].dtype.name == 'category' or
                              (adata.obs[col].dtype == 'object' and adata.obs[col].nunique() < 50)]

        # Get highly variable genes
        if 'highly_variable' in adata.var.columns:
            hvg_genes = adata.var_names[adata.var['highly_variable']].tolist()
        else:
            # Use all genes
            hvg_genes = adata.var_names.tolist()

        # ========================================
        # Step 2: Select analysis type
        # ========================================
        st.header("Step 2: Select analysis type")

        with st.expander("📚 Analysis Guide", expanded=False):
            st.markdown("""
            ### In silico Perturbation
            - **遺伝子摂動**: 特定遺伝子の発現を変化させた際の細胞状態への影響を予測
            - **Perturbation types**:
              - **Knockout**: 遺伝子発現を0に設定
              - **Overexpression**: 遺伝子発現を増加
              - **Custom**: 任意の発現値を設定
            - **結果**: ベクトル場の変化、細胞運命の変化

            ### Least Action Path
            - **最適経路**: 開始状態から終了状態への最も確率の高い遷移経路
            - **用途**: 細胞分化・転写再プログラミングの経路予測
            - **計算**: ベクトル場に基づくaction integral最小化

            ### Jacobian Analysis (Regulatory Network)
            - **Jacobian行列**: 遺伝子間の瞬間的な制御関係
            - **Regulators**: 制御因子（転写因子など）
            - **Effectors**: 被制御因子（標的遺伝子）
            - **結果**: 遺伝子制御ネットワークの定量的マップ

            詳細は [Dynamo Documentation](https://dynamo-release.readthedocs.io/) を参照
            """)

        analysis_type = st.selectbox(
            "Select analysis type",
            ["In silico Perturbation", "Least Action Path", "Jacobian Analysis (Regulatory Network)"],
            help="実行する解析タイプを選択"
        )

        # ========================================
        # In silico Perturbation
        # ========================================
        if analysis_type == "In silico Perturbation":
            st.subheader("In silico Perturbation Settings")

            with st.form("perturbation_form"):
                # Basis selection for Vector Field
                if len(vf_bases) > 0:
                    selected_basis = st.selectbox(
                        "Vector Field basis for perturbation analysis",
                        options=vf_bases,
                        help="Perturbation/Jacobian計算に使用するVector Field basis\n\n"
                             "⚠️ 推奨: PCA basis (pca, rna.pca, mnn.pca等)\n"
                             "• 遺伝子発現空間での定量的解析に最適\n"
                             "• Jacobian行列、ranking/differential解析が可能\n\n"
                             "UMAP basis: 可視化には適していますが、Perturbation解析には非推奨"
                    )

                    # Show warning if UMAP basis is selected
                    if 'umap' in selected_basis.lower() and 'pca' not in selected_basis.lower():
                        st.warning("⚠️ UMAP basisが選択されています。Perturbation解析には**PCA basis**（pca, rna.pca等）を推奨します。")
                    elif 'pca' in selected_basis.lower():
                        st.success("✓ PCA basisが選択されています - Perturbation解析に最適です")
                else:
                    st.error("No Vector Field bases available")
                    st.stop()

                col1, col2 = st.columns(2)

                with col1:
                    # Select genes to perturb
                    st.markdown("#### Genes to perturb")

                    # Gene search
                    gene_search = st.text_input(
                        "Search genes (comma-separated)",
                        help="遺伝子名をカンマ区切りで入力（例: Gata1, Gata2, Tal1）\n全遺伝子を検索可能（大文字小文字を区別しません）"
                    )

                    if gene_search:
                        searched_genes = [g.strip() for g in gene_search.split(',')]
                        # Create case-insensitive mapping
                        gene_name_map = {gene.lower(): gene for gene in adata.var_names}
                        # Filter genes that exist in data (case-insensitive)
                        perturb_genes = [gene_name_map[g.lower()] for g in searched_genes if g.lower() in gene_name_map]

                        if len(perturb_genes) > 0:
                            st.success(f"✓ Found {len(perturb_genes)} genes: {', '.join(perturb_genes)}")
                        else:
                            st.warning("⚠️ No matching genes found")
                            perturb_genes = []
                    else:
                        perturb_genes = []

                    # Or select from list
                    if len(perturb_genes) == 0:
                        perturb_genes = st.multiselect(
                            "Or select from gene list",
                            options=hvg_genes[:500],  # Show top 500 genes
                            default=[],
                            help="リストから選択（Highly Variable Genesの上位500遺伝子のみ）"
                        )

                with col2:
                    # Perturbation type
                    st.markdown("#### Perturbation type")

                    pert_type = st.selectbox(
                        "Perturbation method",
                        ["Knockout", "Overexpression", "Custom expression"],
                        help="遺伝子摂動の方法"
                    )

                    if pert_type == "Knockout":
                        expression_value = 0.0
                        st.info("Expression will be set to 0 (knockout)")
                    elif pert_type == "Overexpression":
                        expression_fold = st.slider(
                            "Overexpression fold change",
                            min_value=2.0,
                            max_value=10.0,
                            value=5.0,
                            step=0.5,
                            help="発現倍率（現在の発現値に対する倍率）"
                        )
                        st.info(f"Expression will be multiplied by {expression_fold}x")
                    else:  # Custom
                        expression_value = st.number_input(
                            "Custom expression value",
                            min_value=0.0,
                            max_value=100.0,
                            value=10.0,
                            step=1.0,
                            help="設定する発現値"
                        )

                # Prediction settings
                st.markdown("#### Prediction settings")

                col1, col2, col3 = st.columns(3)

                with col1:
                    t_end = st.number_input(
                        "Prediction time (t_end)",
                        min_value=1,
                        max_value=100,
                        value=20,
                        help="予測時間の長さ"
                    )

                with col2:
                    # Select cell population
                    if len(categorical_cols) > 0:
                        cell_group_col = st.selectbox(
                            "Cell group column",
                            options=['All cells'] + categorical_cols,
                            help="解析対象の細胞グループを選択するための列"
                        )

                        if cell_group_col != 'All cells':
                            cell_groups = adata.obs[cell_group_col].cat.categories.tolist()
                            selected_groups = st.multiselect(
                                f"Select {cell_group_col}",
                                options=cell_groups,
                                default=cell_groups[:1] if len(cell_groups) > 0 else [],
                                help="解析対象の細胞グループ"
                            )
                        else:
                            cell_group_col = None
                            selected_groups = []
                    else:
                        cell_group_col = None
                        selected_groups = []
                        st.info("No cell grouping available - will use all cells")

                with col3:
                    n_cells = st.number_input(
                        "Number of cells to predict",
                        min_value=10,
                        max_value=1000,
                        value=100,
                        help="予測する細胞数（ランダムサンプリング）"
                    )

                run_perturbation = st.form_submit_button("🧬 Run Perturbation", type="primary")

            if run_perturbation:
                if len(perturb_genes) == 0:
                    st.error("❌ Please select genes to perturb")
                else:
                    st.header("Step 3: Running perturbation analysis")

                    with st.spinner("Running in silico perturbation..."):
                        try:
                            # Prepare expression dict
                            if pert_type == "Knockout":
                                expression_dict = {gene: 0.0 for gene in perturb_genes}
                            elif pert_type == "Overexpression":
                                # Get current mean expression
                                mean_expr = adata[:, perturb_genes].X.mean(axis=0)
                                if hasattr(mean_expr, 'A1'):
                                    mean_expr = mean_expr.A1
                                expression_dict = {gene: mean_expr[i] * expression_fold
                                                 for i, gene in enumerate(perturb_genes)}
                            else:  # Custom
                                expression_dict = {gene: expression_value for gene in perturb_genes}

                            st.info(f"Perturbation targets: {expression_dict}")

                            # Select cells
                            if cell_group_col and len(selected_groups) > 0:
                                cell_mask = adata.obs[cell_group_col].isin(selected_groups)
                                cell_indices = np.where(cell_mask)[0]
                            else:
                                cell_indices = np.arange(adata.n_obs)

                            # Random sample
                            if len(cell_indices) > n_cells:
                                cell_indices = np.random.choice(cell_indices, n_cells, replace=False)

                            st.info(f"Predicting {len(cell_indices)} cells")

                            # Run perturbation prediction
                            st.warning("⏳ This may take several minutes...")

                            # Step 1: Compute perturbation effect
                            st.subheader("1. Computing perturbation effect")

                            progress_bar = st.progress(0)
                            status_text = st.empty()

                            try:
                                status_text.text("Computing perturbation effect...")
                                progress_bar.progress(10)

                                # Prepare expression list
                                if pert_type == "Knockout":
                                    expression_values = [0.0] * len(perturb_genes)
                                elif pert_type == "Overexpression":
                                    mean_expr = adata[:, perturb_genes].X.mean(axis=0)
                                    if hasattr(mean_expr, 'A1'):
                                        mean_expr = mean_expr.A1
                                    expression_values = [mean_expr[i] * expression_fold for i in range(len(perturb_genes))]
                                else:  # Custom
                                    expression_values = [expression_value] * len(perturb_genes)

                                st.info(f"Perturbation: {dict(zip(perturb_genes, expression_values))}")

                                # Debug: Show all Vector Field keys
                                all_vf_keys = [k for k in adata.uns.keys() if k.startswith('VecFld')]
                                st.info(f"🔍 Debug: Available Vector Field keys in adata.uns: **{all_vf_keys}**")
                                st.info(f"🔍 Debug: Selected basis: **{selected_basis}**")

                                # Verify Vector Field exists
                                expected_vf_key = f'VecFld_{selected_basis}'
                                if expected_vf_key not in adata.uns:
                                    st.error(f"""
                                    ❌ **Vector Field not found for basis: {selected_basis}**

                                    Expected key: `{expected_vf_key}`
                                    Available Vector Field keys: {all_vf_keys}

                                    Please ensure you selected a basis that was computed in Dynamo Analysis.
                                    """)
                                    st.stop()

                                st.success(f"✓ Found Vector Field: {expected_vf_key}")

                                # CRITICAL WORKAROUND: Copy Vector Field AND PCA data to work around Dynamo bug
                                # Dynamo's perturbation() has a bug: it internally calls jacobian(basis='pca')
                                # regardless of the basis parameter passed to perturbation()
                                # We need to copy:
                                # 1. VecFld_<basis> → VecFld_pca
                                # 2. X_<basis> → X_pca
                                # 3. PCs_<basis> → PCs (if exists)
                                # 4. pca_mean_<basis> → pca_mean (if exists)

                                st.warning(f"⚠️ Copying Vector Field and PCA data to work around Dynamo bug")
                                st.info(f"Selected basis: {selected_basis} → Copying to: pca")

                                # 1. Copy Vector Field to VecFld_pca
                                adata.uns['VecFld_pca'] = adata.uns[expected_vf_key].copy()
                                st.success(f"✓ Copied {expected_vf_key} → VecFld_pca")

                                # 2. Copy embedding X_<basis> → X_pca
                                embedding_key = f'X_{selected_basis}'
                                if embedding_key in adata.obsm:
                                    adata.obsm['X_pca'] = adata.obsm[embedding_key].copy()
                                    st.success(f"✓ Copied {embedding_key} → X_pca")
                                else:
                                    st.warning(f"⚠️ Embedding {embedding_key} not found")

                                # 3. Try to extract PCs from Vector Field object
                                pcs_extracted = False
                                vf_dict = adata.uns[expected_vf_key]

                                # Check if Vector Field contains PCs information
                                if 'PCs' in vf_dict:
                                    adata.uns['PCs'] = vf_dict['PCs'].copy()
                                    st.success(f"✓ Extracted PCs from {expected_vf_key}")
                                    pcs_extracted = True
                                elif 'pca_dict' in vf_dict and 'PCs' in vf_dict['pca_dict']:
                                    adata.uns['PCs'] = vf_dict['pca_dict']['PCs'].copy()
                                    st.success(f"✓ Extracted PCs from {expected_vf_key}['pca_dict']")
                                    pcs_extracted = True
                                else:
                                    # Try standard keys
                                    pcs_key = f'PCs_{selected_basis}'
                                    if pcs_key in adata.uns:
                                        adata.uns['PCs'] = adata.uns[pcs_key].copy()
                                        st.success(f"✓ Copied {pcs_key} → PCs")
                                        pcs_extracted = True
                                    elif pcs_key in adata.varm:
                                        adata.uns['PCs'] = adata.varm[pcs_key].copy()
                                        st.success(f"✓ Copied varm[{pcs_key}] → uns['PCs']")
                                        pcs_extracted = True

                                if not pcs_extracted:
                                    # Last resort: check if existing PCs can be truncated
                                    if 'PCs' in adata.uns or 'PCs' in adata.varm:
                                        existing_pcs = adata.uns.get('PCs', adata.varm.get('PCs', None))
                                        if existing_pcs is not None:
                                            n_dims = adata.obsm[embedding_key].shape[1]
                                            if existing_pcs.shape[1] >= n_dims:
                                                adata.uns['PCs'] = existing_pcs[:, :n_dims].copy()
                                                st.warning(f"⚠️ Truncated existing PCs to {n_dims} dimensions to match {selected_basis}")
                                                pcs_extracted = True

                                if not pcs_extracted:
                                    st.error(f"❌ Cannot find PCs for {selected_basis} - perturbation may fail")

                                # 4. Try to extract pca_mean from Vector Field object
                                mean_extracted = False
                                if 'pca_mean' in vf_dict:
                                    adata.uns['pca_mean'] = vf_dict['pca_mean'].copy()
                                    st.success(f"✓ Extracted pca_mean from {expected_vf_key}")
                                    mean_extracted = True
                                elif 'pca_dict' in vf_dict and 'pca_mean' in vf_dict['pca_dict']:
                                    adata.uns['pca_mean'] = vf_dict['pca_dict']['pca_mean'].copy()
                                    st.success(f"✓ Extracted pca_mean from {expected_vf_key}['pca_dict']")
                                    mean_extracted = True
                                else:
                                    # Try standard key
                                    mean_key = f'pca_mean_{selected_basis}'
                                    if mean_key in adata.uns:
                                        adata.uns['pca_mean'] = adata.uns[mean_key].copy()
                                        st.success(f"✓ Copied {mean_key} → pca_mean")
                                        mean_extracted = True

                                if not mean_extracted:
                                    st.warning(f"⚠️ Mean for {selected_basis} not found - will use zeros")

                                # Also copy to default VecFld for safety
                                adata.uns['VecFld'] = adata.uns[expected_vf_key].copy()
                                st.success(f"✓ Copied {expected_vf_key} → VecFld (default)")

                                # Get basis for embedding visualization
                                umap_key = [k for k in adata.obsm.keys() if 'umap' in k.lower()][0]
                                basis_name = umap_key.replace('X_', '')

                                # Use default VecFld (no basis parameter)
                                dynamo_basis = None
                                st.info("Using default Vector Field (VecFld) - copied from selected basis")

                                progress_bar.progress(20)

                                # Skip pre-computing Jacobian since we're using default VecFld
                                # The perturbation() function will compute it internally
                                progress_bar.progress(25)

                                # Compute perturbation effect (instantaneous)
                                # API: dyn.pd.perturbation(adata, genes, expression, basis, emb_basis, ...)
                                # Since we copied Vector Field to 'VecFld', we don't specify basis parameter
                                status_text.text("Computing perturbation effect vector...")

                                # Use default VecFld (no basis parameter)
                                dyn.pd.perturbation(
                                    adata,
                                    genes=perturb_genes,
                                    expression=expression_values,
                                    emb_basis=basis_name  # Embedding basis for projection
                                )

                                progress_bar.progress(40)
                                st.success("✓ Perturbation effect computed!")

                                # Step 2: Predict cell fate trajectories
                                # API: dyn.pd.fate(adata, init_cells, basis, t_end, direction, interpolation_num, ...)
                                # Since we're using default VecFld, use basis_name for visualization
                                st.subheader("2. Predicting cell fate trajectories")
                                status_text.text(f"Predicting trajectories for {len(cell_indices)} cells...")

                                # Run fate prediction on selected cells
                                # Use basis_name for embedding visualization
                                dyn.pd.fate(
                                    adata,
                                    init_cells=cell_indices,
                                    basis=basis_name,  # Embedding basis for visualization
                                    t_end=t_end,
                                    direction='forward',
                                    interpolation_num=100
                                )

                                progress_bar.progress(60)
                                status_text.text("✓ Trajectory prediction complete")

                                st.success("✓ Cell fate trajectories computed!")

                                # Visualize trajectories
                                st.subheader("Predicted Trajectories")

                                fig, ax = plt.subplots(figsize=(10, 8))

                                # Plot all cells
                                coords = adata.obsm[umap_key]
                                ax.scatter(coords[:, 0], coords[:, 1],
                                         c='lightgray', s=10, alpha=0.3, label='All cells')

                                # Highlight initial cells
                                selected_coords = adata[cell_indices].obsm[umap_key]
                                ax.scatter(selected_coords[:, 0], selected_coords[:, 1],
                                         c='blue', s=50, alpha=0.6, label='Initial cells')

                                # Plot fate trajectories
                                # Check multiple possible keys where fate might store results
                                fate_key = f'fate_{basis_name}'
                                trajectories_found = False

                                if fate_key in adata.uns:
                                    fate_data = adata.uns[fate_key]
                                    if isinstance(fate_data, (list, np.ndarray)):
                                        # Plot trajectories for first few cells
                                        n_show = min(10, len(cell_indices), len(fate_data))
                                        for i in range(n_show):
                                            if i < len(fate_data):
                                                traj = fate_data[i]
                                                if traj is not None and len(traj) > 0:
                                                    ax.plot(traj[:, 0], traj[:, 1],
                                                          c='red', linewidth=2, alpha=0.5)
                                                    trajectories_found = True
                                        if trajectories_found:
                                            st.info(f"Showing trajectories for {n_show} cells")

                                # Also check for prediction results in obsm
                                elif f'X_{basis_name}_fate' in adata.obsm:
                                    st.info(f"Fate prediction results stored in adata.obsm['X_{basis_name}_fate']")
                                    trajectories_found = True

                                if not trajectories_found:
                                    st.warning(f"Trajectory data not found in expected locations. Check adata.uns['{fate_key}'] or adata.obsm for results.")
                                    st.info("Fate prediction may have succeeded but stored results in a different format.")

                                ax.set_title(f"Predicted Cell Fate after {pert_type}\n({', '.join(perturb_genes)})")
                                ax.set_xlabel(f"{basis_name.upper()} 1")
                                ax.set_ylabel(f"{basis_name.upper()} 2")
                                ax.legend()

                                st.pyplot(fig)
                                plt.close()

                                progress_bar.progress(70)

                                # Step 3: Optionally compute Jacobian
                                st.subheader("3. Computing regulatory network (Jacobian)")
                                status_text.text("Computing Jacobian matrix...")

                                try:
                                    # Use default VecFld (already copied from selected basis)
                                    st.info("Using default Vector Field (VecFld)")
                                    dyn.vf.jacobian(
                                        adata,
                                        regulators=perturb_genes,
                                        effectors=hvg_genes[:100],  # Top 100 genes as effectors
                                        cell_idx=cell_indices[:100]  # Compute for subset of cells
                                    )

                                    progress_bar.progress(90)
                                    st.success("✓ Jacobian computation complete")

                                    # Show Jacobian heatmap
                                    if 'jacobian' in adata.uns:
                                        st.subheader("Regulatory Network (Jacobian)")

                                        fig, ax = plt.subplots(figsize=(12, 8))

                                        # Get Jacobian matrix
                                        jac = adata.uns['jacobian']['jacobian_gene']

                                        # Plot heatmap
                                        max_show = min(20, jac.shape[0], jac.shape[1])
                                        sns.heatmap(jac[:max_show, :max_show],
                                                  cmap='RdBu_r',
                                                  center=0,
                                                  xticklabels=adata.uns['jacobian']['effectors'][:max_show],
                                                  yticklabels=adata.uns['jacobian']['regulators'][:max_show],
                                                  ax=ax)
                                        ax.set_title("Jacobian Matrix (Regulators → Effectors)")
                                        ax.set_xlabel("Effector genes")
                                        ax.set_ylabel("Regulator genes")

                                        st.pyplot(fig)
                                        plt.close()

                                except Exception as e:
                                    st.warning(f"Jacobian computation failed: {str(e)}")
                                    st.info("Trajectory prediction succeeded but Jacobian failed")

                                progress_bar.progress(100)
                                status_text.text("✓ Analysis complete")

                                # Store results
                                st.session_state.dynamo_pert_adata = adata
                                st.session_state.dynamo_pert_genes = perturb_genes
                                st.session_state.dynamo_pert_type = pert_type
                                st.session_state.dynamo_analysis_type = "perturbation"
                                st.session_state.dynamo_pert_complete = True

                                st.success("""
                                ✅ **Perturbation analysis completed!**

                                Results:
                                - Perturbation trajectories computed
                                - Regulatory network (Jacobian) analyzed
                                - Results saved to adata
                                """)

                            except Exception as e:
                                st.error(f"Perturbation computation failed: {str(e)}")
                                st.info("This may happen if the vector field is incomplete or incompatible")
                                st.exception(e)

                        except Exception as e:
                            st.error(f"❌ Perturbation analysis failed: {str(e)}")
                            st.exception(e)

        # ========================================
        # Least Action Path
        # ========================================
        elif analysis_type == "Least Action Path":
            st.subheader("Least Action Path Settings")

            st.markdown("""
            最適遷移経路を計算するには、開始細胞と終了細胞を指定する必要があります。
            """)

            # Cell group column selection OUTSIDE the form
            if len(categorical_cols) > 0:
                lap_group_col = st.selectbox(
                    "Cluster identity column",
                    options=categorical_cols,
                    help="クラスター/細胞タイプ等のカテゴリカル列"
                )

                lap_groups = adata.obs[lap_group_col].cat.categories.tolist()
            else:
                st.warning("No categorical columns available for cell selection")
                lap_group_col = None

            # Only show form if we have a valid column
            if lap_group_col:
                with st.form("lap_form"):
                    col1, col2 = st.columns(2)

                    with col1:
                        st.markdown("#### Start cluster")
                        start_group = st.selectbox(
                            "Select start cluster",
                            options=lap_groups,
                            key="lap_start_group",
                            help="開始クラスター"
                        )

                    with col2:
                        st.markdown("#### End cluster")
                        end_group = st.selectbox(
                            "Select end cluster",
                            options=lap_groups,
                            key="lap_end_group",
                            index=min(1, len(lap_groups)-1),
                            help="終了クラスター"
                        )

                    # Additional settings
                    st.markdown("#### Path calculation settings")
                    col1, col2, col3 = st.columns(3)

                    with col1:
                        n_paths = st.number_input(
                            "Number of paths",
                            min_value=1,
                            max_value=20,
                            value=5,
                            help="計算する経路の数"
                        )

                    with col2:
                        lap_basis = st.selectbox(
                            "Basis for calculation",
                            options=['pca', 'umap'],
                            index=0,
                            help="経路計算の基底"
                        )

                    with col3:
                        n_steps = st.number_input(
                            "Number of steps",
                            min_value=20,
                            max_value=200,
                            value=100,
                            help="経路の分割数"
                        )

                    run_lap = st.form_submit_button("🛤️ Calculate Least Action Path", type="primary")
            else:
                run_lap = False

            if run_lap:
                if lap_group_col is None:
                    st.error("❌ Cell selection requires categorical columns")
                else:
                    st.header("Step 3: Calculating least action path")

                    with st.spinner("Computing optimal trajectory..."):
                        try:
                            # Get start and end cell indices
                            start_cells = adata.obs_names[adata.obs[lap_group_col] == start_group]
                            end_cells = adata.obs_names[adata.obs[lap_group_col] == end_group]

                            st.info(f"Start cells: {len(start_cells)}, End cells: {len(end_cells)}")

                            # Select representative cells
                            n_start = min(n_paths, len(start_cells))
                            n_end = min(n_paths, len(end_cells))

                            start_cells_selected = np.random.choice(start_cells, n_start, replace=False)
                            end_cells_selected = np.random.choice(end_cells, n_end, replace=False)

                            st.warning("⏳ This may take several minutes...")

                            # Progress tracking
                            progress_bar = st.progress(0)
                            status_text = st.empty()

                            status_text.text("Preparing least action path calculation...")
                            progress_bar.progress(10)

                            # Get start and end cell indices (integers)
                            start_indices = [np.where(adata.obs_names == cell)[0][0] for cell in start_cells_selected]
                            end_indices = [np.where(adata.obs_names == cell)[0][0] for cell in end_cells_selected]

                            st.info(f"""
                            **Least Action Path Configuration:**
                            - Start: {start_group} ({n_start} cells)
                            - End: {end_group} ({n_end} cells)
                            - Basis: {lap_basis}
                            - Steps: {n_steps}
                            - Paths: {n_paths}
                            """)

                            progress_bar.progress(20)

                            # Compute least action paths
                            status_text.text("Computing least action paths...")

                            try:
                                # Run least action path calculation
                                dyn.pd.least_action(
                                    adata,
                                    init_cells=start_indices,
                                    target_cells=end_indices,
                                    basis=lap_basis,
                                    num_t=n_steps,
                                    adj_key='distances'  # Use nearest neighbor distances
                                )

                                progress_bar.progress(80)
                                status_text.text("✓ Path computation complete")

                                st.success("✓ Least action paths computed!")

                                # Visualize paths
                                st.subheader("Optimal Transition Paths")

                                # Get basis key
                                if lap_basis == 'umap':
                                    basis_key = [k for k in adata.obsm.keys() if 'umap' in k.lower()][0]
                                elif lap_basis == 'pca':
                                    basis_key = 'X_pca'
                                else:
                                    basis_key = f'X_{lap_basis}'

                                fig, ax = plt.subplots(figsize=(12, 10))

                                # Plot all cells
                                coords = adata.obsm[basis_key]
                                ax.scatter(coords[:, 0], coords[:, 1],
                                         c='lightgray', s=10, alpha=0.3, label='All cells')

                                # Highlight start and end clusters
                                start_mask = adata.obs[lap_group_col] == start_group
                                end_mask = adata.obs[lap_group_col] == end_group

                                start_coords = adata[start_mask].obsm[basis_key]
                                end_coords = adata[end_mask].obsm[basis_key]

                                ax.scatter(start_coords[:, 0], start_coords[:, 1],
                                         c='blue', s=50, alpha=0.6, label=f'Start: {start_group}')
                                ax.scatter(end_coords[:, 0], end_coords[:, 1],
                                         c='red', s=50, alpha=0.6, label=f'End: {end_group}')

                                # Plot paths
                                if 'LAP' in adata.uns:
                                    lap_data = adata.uns['LAP']

                                    # Plot paths
                                    for i in range(min(n_paths, len(start_indices))):
                                        try:
                                            path_key = f'LAP_{i}'
                                            if path_key in lap_data:
                                                path = lap_data[path_key]
                                                ax.plot(path[:, 0], path[:, 1],
                                                      c='orange', linewidth=2, alpha=0.6)
                                        except:
                                            continue

                                    st.info(f"Computed {min(n_paths, len(start_indices))} optimal paths")
                                else:
                                    st.warning("Path data not found in expected format")

                                ax.set_title(f"Least Action Paths: {start_group} → {end_group}")
                                ax.set_xlabel(f"{lap_basis.upper()} 1")
                                ax.set_ylabel(f"{lap_basis.upper()} 2")
                                ax.legend()

                                st.pyplot(fig)
                                plt.close()

                                progress_bar.progress(100)
                                status_text.text("✓ Analysis complete")

                                # Store results
                                st.session_state.dynamo_pert_adata = adata
                                st.session_state.dynamo_lap_start = start_group
                                st.session_state.dynamo_lap_end = end_group
                                st.session_state.dynamo_analysis_type = "lap"
                                st.session_state.dynamo_pert_complete = True

                                st.success("""
                                ✅ **Least action path analysis completed!**

                                Results:
                                - Optimal transition paths computed
                                - Paths represent the most probable cell state transitions
                                - Results saved to adata.uns['LAP']
                                """)

                            except Exception as e:
                                st.error(f"Least action path computation failed: {str(e)}")
                                st.info("This may happen if:")
                                st.markdown("""
                                - The vector field is incomplete
                                - The basis doesn't have sufficient dimensions
                                - Start and end states are too similar
                                """)
                                st.exception(e)

                        except Exception as e:
                            st.error(f"❌ Least action path calculation failed: {str(e)}")
                            st.exception(e)

        # ========================================
        # Jacobian Analysis
        # ========================================
        elif analysis_type == "Jacobian Analysis (Regulatory Network)":
            st.subheader("Jacobian Analysis Settings")

            st.markdown("""
            Jacobian行列を計算して遺伝子制御ネットワークを推論します。
            """)

            with st.form("jacobian_form"):
                # Basis selection for Vector Field
                if len(vf_bases) > 0:
                    jac_selected_basis = st.selectbox(
                        "Vector Field basis for Jacobian analysis",
                        options=vf_bases,
                        help="Jacobian計算に使用するVector Field basis\n\n"
                             "⚠️ 推奨: PCA basis (pca, rna.pca, mnn.pca等)\n"
                             "• 遺伝子発現空間での定量的解析に最適\n"
                             "• 制御ネットワークの正確な推論が可能\n\n"
                             "UMAP basis: 可視化には適していますが、Jacobian解析には非推奨"
                    )

                    # Show warning if UMAP basis is selected
                    if 'umap' in jac_selected_basis.lower() and 'pca' not in jac_selected_basis.lower():
                        st.warning("⚠️ UMAP basisが選択されています。Jacobian解析には**PCA basis**（pca, rna.pca等）を推奨します。")
                    elif 'pca' in jac_selected_basis.lower():
                        st.success("✓ PCA basisが選択されています - Jacobian解析に最適です")
                else:
                    st.error("No Vector Field bases available")
                    st.stop()

                col1, col2 = st.columns(2)

                with col1:
                    st.markdown("#### Regulator genes")

                    reg_search = st.text_input(
                        "Search regulator genes (comma-separated)",
                        help="制御因子（転写因子など）をカンマ区切りで入力\n全遺伝子を検索可能（大文字小文字を区別しません）"
                    )

                    if reg_search:
                        searched_regs = [g.strip() for g in reg_search.split(',')]
                        # Create case-insensitive mapping
                        gene_name_map = {gene.lower(): gene for gene in adata.var_names}
                        # Filter genes that exist in data (case-insensitive)
                        regulators = [gene_name_map[g.lower()] for g in searched_regs if g.lower() in gene_name_map]

                        if len(regulators) > 0:
                            st.success(f"✓ Found {len(regulators)} regulators: {', '.join(regulators)}")
                        else:
                            st.warning("⚠️ No matching genes found")
                            regulators = []
                    else:
                        regulators = []

                    if len(regulators) == 0:
                        regulators = st.multiselect(
                            "Or select from gene list",
                            options=hvg_genes[:200],
                            default=[],
                            key="jac_regulator_multiselect",
                            help="リストから選択（Highly Variable Genesの上位200遺伝子のみ）"
                        )

                with col2:
                    st.markdown("#### Effector genes")

                    eff_search = st.text_input(
                        "Search effector genes (comma-separated)",
                        help="被制御因子（標的遺伝子）をカンマ区切りで入力\n全遺伝子を検索可能（大文字小文字を区別しません）"
                    )

                    if eff_search:
                        searched_effs = [g.strip() for g in eff_search.split(',')]
                        # Create case-insensitive mapping
                        gene_name_map = {gene.lower(): gene for gene in adata.var_names}
                        # Filter genes that exist in data (case-insensitive)
                        effectors = [gene_name_map[g.lower()] for g in searched_effs if g.lower() in gene_name_map]

                        if len(effectors) > 0:
                            st.success(f"✓ Found {len(effectors)} effectors: {', '.join(effectors)}")
                        else:
                            st.warning("⚠️ No matching genes found")
                            effectors = []
                    else:
                        effectors = []

                    if len(effectors) == 0:
                        use_all_genes = st.checkbox(
                            "Use all genes as effectors",
                            value=False,
                            help="全遺伝子を被制御因子として使用"
                        )

                        if use_all_genes:
                            effectors = adata.var_names.tolist()
                            st.info(f"Using all {len(effectors)} genes as effectors")
                        else:
                            effectors = st.multiselect(
                                "Or select from gene list",
                                options=hvg_genes[:200],
                                default=[],
                                key="jac_effector_multiselect",
                                help="リストから選択（Highly Variable Genesの上位200遺伝子のみ）"
                            )

                # Cell selection
                col1, col2 = st.columns(2)

                with col1:
                    if len(categorical_cols) > 0:
                        jac_cell_col = st.selectbox(
                            "Cell group for analysis",
                            options=['All cells'] + categorical_cols,
                            help="解析対象の細胞グループ"
                        )

                        if jac_cell_col != 'All cells':
                            jac_groups = adata.obs[jac_cell_col].cat.categories.tolist()
                            jac_selected_groups = st.multiselect(
                                f"Select {jac_cell_col}",
                                options=jac_groups,
                                default=jac_groups,
                                help="解析対象グループ"
                            )
                        else:
                            jac_cell_col = None
                            jac_selected_groups = []
                    else:
                        jac_cell_col = None
                        jac_selected_groups = []

                with col2:
                    max_cells = st.number_input(
                        "Maximum cells for computation",
                        min_value=100,
                        max_value=10000,
                        value=1000,
                        help="Jacobian計算に使用する最大細胞数"
                    )

                run_jacobian = st.form_submit_button("🔗 Compute Jacobian", type="primary")

            if run_jacobian:
                if len(regulators) == 0 or len(effectors) == 0:
                    st.error("❌ Please select both regulators and effectors")
                else:
                    st.header("Step 3: Computing Jacobian matrix")

                    with st.spinner("Computing regulatory network..."):
                        try:
                            # Select cells
                            if jac_cell_col and len(jac_selected_groups) > 0:
                                cell_mask = adata.obs[jac_cell_col].isin(jac_selected_groups)
                                cell_indices = np.where(cell_mask)[0]
                            else:
                                cell_indices = np.arange(adata.n_obs)

                            # Sample if too many
                            if len(cell_indices) > max_cells:
                                cell_indices = np.random.choice(cell_indices, max_cells, replace=False)

                            st.info(f"Computing Jacobian for {len(regulators)} regulators × {len(effectors)} effectors on {len(cell_indices)} cells")

                            st.warning("⏳ This may take several minutes...")

                            # Determine correct basis to use
                            if jac_selected_basis == 'default':
                                jac_dynamo_basis = None
                                st.info("Using default Vector Field (VecFld)")
                            else:
                                jac_dynamo_basis = jac_selected_basis
                                st.info(f"Using Vector Field: VecFld_{jac_dynamo_basis}")

                            # Compute Jacobian
                            if jac_dynamo_basis is not None:
                                dyn.vf.jacobian(
                                    adata,
                                    regulators=regulators,
                                    effectors=effectors,
                                    cell_idx=cell_indices,
                                    basis=jac_dynamo_basis  # CRITICAL: Specify Vector Field basis
                                )
                            else:
                                # Use default VecFld (no basis parameter)
                                dyn.vf.jacobian(
                                    adata,
                                    regulators=regulators,
                                    effectors=effectors,
                                    cell_idx=cell_indices
                                )

                            st.success("✓ Jacobian computation complete")

                            # Store results
                            st.session_state.dynamo_pert_adata = adata
                            st.session_state.dynamo_pert_regulators = regulators
                            st.session_state.dynamo_pert_effectors = effectors
                            st.session_state.dynamo_analysis_type = "jacobian"
                            st.session_state.dynamo_pert_complete = True

                            st.success("""
                            ✅ **Jacobian analysis completed!**

                            Results stored in adata.uns['jacobian']
                            """)

                            # Visualize Jacobian
                            if 'jacobian' in adata.uns:
                                st.subheader("Regulatory Network Heatmap")

                                jac_data = adata.uns['jacobian']
                                jac_matrix = jac_data['jacobian_gene']

                                # Limit visualization size
                                max_show = 30
                                n_reg_show = min(len(regulators), max_show)
                                n_eff_show = min(len(effectors), max_show)

                                fig, ax = plt.subplots(figsize=(max(10, n_eff_show * 0.5),
                                                               max(8, n_reg_show * 0.4)))

                                sns.heatmap(
                                    jac_matrix[:n_reg_show, :n_eff_show],
                                    cmap='RdBu_r',
                                    center=0,
                                    xticklabels=jac_data['effectors'][:n_eff_show],
                                    yticklabels=jac_data['regulators'][:n_reg_show],
                                    cbar_kws={'label': 'Jacobian value'},
                                    ax=ax
                                )

                                ax.set_title("Gene Regulatory Network (Jacobian Matrix)")
                                ax.set_xlabel("Effector genes")
                                ax.set_ylabel("Regulator genes")
                                plt.tight_layout()

                                st.pyplot(fig)
                                plt.close()

                                # Show top interactions
                                st.subheader("Top Regulatory Interactions")

                                # Flatten Jacobian and get top values
                                jac_flat = jac_matrix.flatten()
                                jac_abs = np.abs(jac_flat)
                                top_indices = np.argsort(jac_abs)[::-1][:20]

                                interactions = []
                                for idx in top_indices:
                                    reg_idx = idx // len(jac_data['effectors'])
                                    eff_idx = idx % len(jac_data['effectors'])

                                    interactions.append({
                                        'Regulator': jac_data['regulators'][reg_idx],
                                        'Effector': jac_data['effectors'][eff_idx],
                                        'Jacobian': jac_flat[idx],
                                        'Abs_Jacobian': jac_abs[idx]
                                    })

                                interactions_df = pd.DataFrame(interactions)
                                st.dataframe(interactions_df, use_container_width=True)

                                # Download interactions
                                csv = interactions_df.to_csv(index=False)
                                st.download_button(
                                    "⬇️ Download top interactions (CSV)",
                                    csv,
                                    "jacobian_top_interactions.csv",
                                    "text/csv"
                                )

                        except Exception as e:
                            st.error(f"❌ Jacobian computation failed: {str(e)}")
                            st.exception(e)

        # ========================================
        # Step 4: Download results
        # ========================================
        if st.session_state.dynamo_pert_complete:
            st.header("Step 4: Download results")

            # Save updated adata
            output_path = os.path.join(dynamo_pert_temp_dir, "perturbation_result.h5ad")

            # Remove existing file if it exists to avoid file locking issues
            if os.path.exists(output_path):
                try:
                    os.remove(output_path)
                except OSError as e:
                    st.warning(f"Could not remove existing file: {e}")

            st.session_state.dynamo_pert_adata.write_h5ad(output_path, compression="gzip")

            with open(output_path, "rb") as f:
                result_bytes = f.read()

            # Generate filename based on analysis type
            analysis_type = st.session_state.get('dynamo_analysis_type', 'perturbation')

            if analysis_type == "perturbation":
                # In silico Perturbation
                target_genes = st.session_state.get('dynamo_pert_genes', [])
                if target_genes and len(target_genes) > 0:
                    if len(target_genes) > 5:
                        gene_str = '_'.join(target_genes[:5]) + '_etc'
                    else:
                        gene_str = '_'.join(target_genes)
                    output_filename = f"perturbation.{gene_str}.h5ad"
                else:
                    output_filename = "perturbation.h5ad"

            elif analysis_type == "jacobian":
                # Jacobian Analysis
                regulators = st.session_state.get('dynamo_pert_regulators', [])
                if regulators and len(regulators) > 0:
                    if len(regulators) > 5:
                        gene_str = '_'.join(regulators[:5]) + '_etc'
                    else:
                        gene_str = '_'.join(regulators)
                    output_filename = f"jacobian.{gene_str}.h5ad"
                else:
                    output_filename = "jacobian.h5ad"

            elif analysis_type == "lap":
                # Least Action Path
                start = st.session_state.get('dynamo_lap_start', 'start')
                end = st.session_state.get('dynamo_lap_end', 'end')
                output_filename = f"LAP.{start}_to_{end}.h5ad"

            else:
                # Fallback
                output_filename = "dynamo_result.h5ad"

            st.download_button(
                label="⬇️ Download perturbation result (h5ad)",
                data=result_bytes,
                file_name=output_filename,
                mime="application/octet-stream",
                type="primary"
            )

            st.info("""
            ### 次のステップ

            ダウンロードしたh5adファイルには以下が含まれています：
            - Jacobian matrix (`adata.uns['jacobian']`)
            - Regulatory network information

            #### Pythonでの追加解析:
            ```python
            import dynamo as dyn
            import scanpy as sc

            # Load result
            adata = sc.read_h5ad('perturbation_result.h5ad')

            # Visualize regulatory network
            dyn.pl.jacobian_heatmap(adata)

            # Network analysis
            dyn.pl.jacobian_kinetics(adata, genes=['gene1', 'gene2'])

            # Compute in silico perturbation
            dyn.pd.perturbation(
                adata,
                genes=['gene1'],
                expression=[0],  # knockout
                t_end=20
            )

            # Visualize perturbation effects
            dyn.pl.perturbation(adata, genes=['gene1'])
            ```

            詳細は [Dynamo Documentation](https://dynamo-release.readthedocs.io/) を参照してください。
            """)

else:
    st.info("👆 Dynamo解析済みh5adファイルをアップロードして開始してください。")
