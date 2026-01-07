"""
scVelo Visualization
Visualize RNA velocity results from scVelo analysis
"""

import streamlit as st
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import io
import time
from helper_func import clear_old_directories, clear_old_files
from streamlit_sortables import sort_items

@st.cache_data
def create_cell_color_mapping(cell_list, palette_name):
    """
    細胞名/クラスターと色の一貫したマッピングを作成する関数

    Parameters
    ----------
    cell_list : list
        細胞名/クラスター名のリスト
    palette_name : str
        使用する離散カラーパレット名

    Returns
    -------
    dict
        細胞名/クラスター名をキー、色をバリューとする辞書
    """
    n_cells = len(cell_list)
    base_palette = sns.color_palette(palette_name)
    base_n = len(base_palette)

    if n_cells <= base_n:
        colors = base_palette[:n_cells]
    else:
        colors = sns.color_palette(palette_name, n_colors=n_cells)

    return {cell: color for cell, color in zip(cell_list, colors)}

st.set_page_config(page_title="scVelo Visualization", page_icon="🌊", layout="wide")

st.title("🌊 scVelo Visualization")
st.markdown("""
scVelo解析結果を可視化します。

### 可視化の種類
1. **Velocity embedding stream**: 連続的なvelocityフロー
2. **Velocity embedding grid**: グリッド状のvelocity矢印
3. **Phase portraits**: Spliced/unspliced遺伝子動態
4. **Velocity confidence**: Velocityの信頼度
5. **Velocity pseudotime**: Velocityベースの疑似時間
6. **Gene-specific velocity**: 特定遺伝子のvelocity
7. **PAGA**: Partition-based graph abstraction（クラスター間接続性）

**注:** Pseudotimeにおける遺伝子発現変化は **Pseudotime gene expression** appで可視化できます。

### 参考
- [scVelo Basics](https://scvelo.readthedocs.io/en/stable/VelocityBasics.html)
- [PAGA](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.paga.html)
""")

# Initialize session state
if "velocity_vis_temp_dir" not in st.session_state:
    velocity_vis_temp_dir = os.path.join("temp", f"velocity_vis_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(velocity_vis_temp_dir, exist_ok=True)
    st.session_state.velocity_vis_temp_dir = velocity_vis_temp_dir
else:
    velocity_vis_temp_dir = st.session_state.velocity_vis_temp_dir

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload scVelo result")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file (scVelo result)",
    type=['h5ad'],
    key="velocity_vis_h5ad_upload",
    help="scVelo analysisアプリで生成されたh5adファイル"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load data
    if ("velocity_vis_adata" not in st.session_state or
        st.session_state.get("velocity_vis_uploaded_file") != uploaded_h5ad.name):

        with st.spinner("Loading data..."):
            # Save uploaded file temporarily
            temp_h5ad_path = os.path.join(velocity_vis_temp_dir, "input.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            # Read data with error handling for malformed h5ad
            try:
                adata = sc.read_h5ad(temp_h5ad_path)
            except ValueError as e:
                if "Data must be 1-dimensional" in str(e):
                    st.warning("⚠️ Detected malformed metadata in h5ad file. Attempting to fix...")

                    # Try to read with backed mode and fix the issue
                    try:
                        import h5py
                        import anndata as ad

                        # Read with backed mode
                        adata_backed = ad.read_h5ad(temp_h5ad_path, backed='r')

                        # Create a new AnnData object with cleaned metadata
                        adata = ad.AnnData(
                            X=adata_backed.X[:],
                            obs=adata_backed.obs.copy(),
                            var=adata_backed.var.copy(),
                            uns=adata_backed.uns.copy() if len(adata_backed.uns) > 0 else {},
                            obsm=adata_backed.obsm.copy() if len(adata_backed.obsm) > 0 else {},
                            varm=adata_backed.varm.copy() if len(adata_backed.varm) > 0 else {},
                            layers={k: adata_backed.layers[k][:] for k in adata_backed.layers.keys()} if len(adata_backed.layers) > 0 else {}
                        )

                        # Clean obs columns - fix any 2D arrays
                        for col in adata.obs.columns:
                            if hasattr(adata.obs[col], 'shape') and len(adata.obs[col].shape) > 1:
                                if adata.obs[col].shape[1] == 1:
                                    adata.obs[col] = adata.obs[col].values.flatten()
                                else:
                                    st.warning(f"Removing malformed column '{col}' from metadata")
                                    adata.obs = adata.obs.drop(columns=[col])

                        # Clean var columns
                        for col in adata.var.columns:
                            if hasattr(adata.var[col], 'shape') and len(adata.var[col].shape) > 1:
                                if adata.var[col].shape[1] == 1:
                                    adata.var[col] = adata.var[col].values.flatten()
                                else:
                                    st.warning(f"Removing malformed column '{col}' from gene metadata")
                                    adata.var = adata.var.drop(columns=[col])

                        st.success("✓ Successfully cleaned and loaded h5ad file")

                    except Exception as fix_error:
                        st.error(f"""
                        ❌ **Failed to load h5ad file**

                        Error: {str(e)}

                        This h5ad file has malformed metadata (2D arrays in obs/var columns).

                        **Possible solutions:**
                        1. Re-run the velocity analysis to generate a new h5ad file
                        2. Open the file in Python and manually clean the metadata:

                        ```python
                        import scanpy as sc
                        import anndata as ad

                        # Load the file
                        adata = sc.read_h5ad('your_file.h5ad')

                        # Check for problematic columns
                        for col in adata.obs.columns:
                            if hasattr(adata.obs[col], 'shape'):
                                print(f"{{col}}: {{adata.obs[col].shape}}")

                        # Fix 2D columns by flattening
                        for col in adata.obs.columns:
                            if hasattr(adata.obs[col], 'shape') and len(adata.obs[col].shape) > 1:
                                adata.obs[col] = adata.obs[col].values.flatten()

                        # Save cleaned file
                        adata.write_h5ad('cleaned_file.h5ad')
                        ```

                        Attempted fix error: {str(fix_error)}
                        """)
                        st.stop()
                else:
                    raise

            st.session_state.velocity_vis_adata = adata
            st.session_state.velocity_vis_uploaded_file = uploaded_h5ad.name

            st.info(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    adata = st.session_state.velocity_vis_adata

    # Check required data (also accept _scvelo/_deepvelo suffixed versions)
    has_velocity = ('velocity' in adata.layers or
                    'velocity_scvelo' in adata.layers or
                    'velocity_deepvelo' in adata.layers)
    has_velocity_graph = ('velocity_graph' in adata.uns or
                          'velocity_graph_scvelo' in adata.uns or
                          'velocity_graph_deepvelo' in adata.uns)

    required_data = []
    if not has_velocity:
        required_data.append("velocity layer")
    if not has_velocity_graph:
        required_data.append("velocity_graph")

    if required_data:
        st.error(f"""
        ❌ **Missing required data: {', '.join(required_data)}**

        このファイルはscVelo解析結果ではないようです。
        scVelo analysisアプリで解析したh5adファイルをアップロードしてください。
        """)
        st.stop()

    # ========================================
    # Velocity Source Selection (scVelo vs DeepVelo)
    # ========================================
    def swap_velocity_data(adata, suffix):
        """
        Swap all velocity-related data to use the specified source.
        suffix: 'scvelo' or 'deepvelo'
        """
        # layers
        layer_keys = ['velocity', 'velocity_u']
        for key in layer_keys:
            src_key = f'{key}_{suffix}'
            if src_key in adata.layers:
                adata.layers[key] = adata.layers[src_key].copy()

        # uns (sparse matrices)
        uns_keys = ['velocity_graph', 'velocity_graph_neg', 'velocity_params']
        for key in uns_keys:
            src_key = f'{key}_{suffix}'
            if src_key in adata.uns:
                adata.uns[key] = adata.uns[src_key].copy()

        # obs (cell-level metrics)
        obs_keys = ['velocity_pseudotime', 'velocity_confidence', 'velocity_length',
                    'velocity_self_transition', 'velocity_confidence_transition']
        for key in obs_keys:
            src_key = f'{key}_{suffix}'
            if src_key in adata.obs.columns:
                adata.obs[key] = adata.obs[src_key].copy()

    # Check for available velocity sources
    velocity_sources = {}  # {display_name: suffix}
    if 'velocity_scvelo' in adata.layers:
        velocity_sources["scVelo (Dynamical)"] = 'scvelo'
    if 'velocity_deepvelo' in adata.layers:
        velocity_sources["DeepVelo (GNN)"] = 'deepvelo'

    # Show velocity source selector
    if len(velocity_sources) >= 1:
        st.subheader("🔄 Velocity source")

        # Show velocity statistics comparison
        with st.expander("📈 Velocity comparison statistics", expanded=False):
            stat_data = []
            for name, suffix in velocity_sources.items():
                v = adata.layers[f'velocity_{suffix}']
                non_nan_genes = np.sum(~np.all(np.isnan(v), axis=0))
                stat_data.append({
                    "Source": name,
                    "Mean": f"{np.nanmean(v):.6f}",
                    "Std": f"{np.nanstd(v):.6f}",
                    "Non-NaN genes": non_nan_genes
                })

            # Correlation between sources if both exist
            if len(velocity_sources) == 2:
                v1 = adata.layers['velocity_scvelo'].flatten()
                v2 = adata.layers['velocity_deepvelo'].flatten()
                mask = ~np.isnan(v1) & ~np.isnan(v2)
                if np.sum(mask) > 0:
                    corr = np.corrcoef(v1[mask], v2[mask])[0, 1]
                    st.metric("Correlation (scVelo vs DeepVelo)", f"{corr:.4f}")

            st.dataframe(pd.DataFrame(stat_data), hide_index=True)

        if len(velocity_sources) > 1:
            # Multiple sources - show selector
            col_vel1, col_vel2 = st.columns([2, 3])
            source_names = list(velocity_sources.keys())
            with col_vel1:
                selected_velocity = st.radio(
                    "Select velocity to visualize",
                    source_names,
                    index=0,  # Default to first (scVelo)
                    help="異なるvelocity計算結果を切り替えて比較できます"
                )
            with col_vel2:
                if "scVelo" in selected_velocity:
                    st.info("📊 **scVelo**: ODE-based dynamical model")
                elif "DeepVelo" in selected_velocity:
                    st.info("🧠 **DeepVelo**: Graph Neural Network model")

            # Swap ALL velocity-related data
            selected_suffix = velocity_sources[selected_velocity]
            swap_velocity_data(adata, selected_suffix)

            # Check if velocity_graph was swapped successfully
            if f'velocity_graph_{selected_suffix}' in adata.uns:
                st.success(f"✓ Using {selected_velocity} (all data swapped)")
            else:
                # velocity_graph for this source doesn't exist - need to recompute
                st.warning(f"⚠️ velocity_graph_{selected_suffix} not found, recomputing...")
                # Remove old velocity_graph to force recomputation
                if 'velocity_graph' in adata.uns:
                    del adata.uns['velocity_graph']
                if 'velocity_graph_neg' in adata.uns:
                    del adata.uns['velocity_graph_neg']
                with st.spinner("Computing velocity graph..."):
                    scv.tl.velocity_graph(adata)
                st.success(f"✓ Using {selected_velocity} (velocity_graph recomputed)")
        else:
            # Single source - use that source
            source_name = list(velocity_sources.keys())[0]
            selected_suffix = velocity_sources[source_name]
            swap_velocity_data(adata, selected_suffix)

            if "scVelo" in source_name:
                st.info("📊 Using **scVelo** velocity")
            elif "DeepVelo" in source_name:
                st.info("🧠 Using **DeepVelo** velocity")
    else:
        # No specific source found, use original velocity layer
        if 'velocity' in adata.layers:
            st.info("📊 Using velocity layer (default)")
        # Ensure velocity_graph exists
        if 'velocity_graph' not in adata.uns:
            with st.spinner("Computing velocity graph..."):
                scv.tl.velocity_graph(adata)

    # Get available embeddings (UMAP, tSNE, PCA, etc.)
    available_bases = []
    for key in adata.obsm.keys():
        if key.startswith('X_'):
            base_name = key.replace('X_', '')
            available_bases.append(base_name)

    if not available_bases:
        st.error("❌ No embeddings found! Please ensure your h5ad file contains UMAP or other embeddings.")
        st.stop()

    # Get available categorical columns for coloring
    categorical_columns = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == 'category' or adata.obs[col].nunique() < 50:
            categorical_columns.append(col)

    # Get available continuous columns
    continuous_columns = []
    for col in adata.obs.columns:
        if pd.api.types.is_numeric_dtype(adata.obs[col]) and adata.obs[col].nunique() > 10:
            continuous_columns.append(col)

    # Get gene list
    gene_list = adata.var_names.tolist()

    # ========================================
    # Step 2: Configure visualization
    # ========================================
    st.header("Step 2: Configure visualization")

    with st.expander("📚 Visualization Guide", expanded=False):
        st.markdown("""
        ### Embedding Basis（埋め込み空間）
        - **複数選択可能**: UMAP、tSNE、PCA等を同時に選択できます
        - 各embeddingに対して個別のプロットが生成されます
        - デフォルトでUMAPが選択されます（利用可能な場合）

        ### Coloring（色付け）
        - **Cluster identity**: メタデータから色付けする列を選択
        - **Discrete colormap**: カテゴリカル変数用（クラスター、細胞タイプ等）
          - 可視化タイプに応じて自動的に表示/非表示が切り替わります
        - **Continuous colormap**: 連続変数用（遺伝子発現、pseudotime等）
          - 可視化タイプに応じて自動的に表示/非表示が切り替わります
        - **Show cluster names**: 凡例の表示/非表示を切り替え

        ### Cluster Settings（サイドバー）
        - **Change cluster order**: クラスターの表示順序をドラッグ&ドロップで変更
        - 順序変更は色の割り当てにも反映されます

        ### Velocity Plots
        - **Stream plot**: 連続的なvelocityフロー
          - 全体的な流れと軌跡を把握
          - Advanced options: density, smoothing, min_mass
        - **Grid plot**: グリッド状のvelocity矢印
          - 個々の細胞のvelocity方向を詳細に把握
          - Advanced options: scale, min_mass, alpha (transparency)

        ### Advanced Options（詳細設定）
        - **Stream density**: ストリームラインの密度（0.5-3.0）
        - **Smoothing**: ストリームラインの滑らかさ（0.0-1.0）
        - **Min mass**: 矢印を表示する最小質量（閾値）
        - **Scale**: 矢印のスケール（0.1-2.0）
        - **Alpha**: 矢印の透明度（0.1-1.0）

        ### PAGA (Partition-based graph abstraction)
        - クラスター間の接続性を可視化
        - Velocity-directed PAGAはクラスター間の方向性のある遷移を表示
        - **PAGA graph**: クラスター間の接続グラフ
        - **PAGA compare**: PAGAグラフとUMAP/tSNEの比較表示
        - PAGAが計算されていない場合は表示されません

        ### Download（ダウンロード）
        - **PNG format**: 高解像度ビットマップ（300 DPI）
        - **PDF format**: ベクター形式（論文投稿・印刷に最適）
        - 生成された全てのプロットが両形式でダウンロード可能
        """)

    # Check if PAGA is available
    has_paga = 'paga' in adata.uns

    # Visualization type
    viz_types = [
        "Velocity embedding stream",
        "Velocity embedding grid",
        "Phase portraits",
        "Velocity confidence",
        "Velocity pseudotime",
        "Gene velocity",
    ]

    if has_paga:
        viz_types.append("PAGA")

    viz_type = st.selectbox(
        "Visualization type",
        viz_types,
        help="可視化のタイプを選択"
    )

    if not has_paga and viz_type == "PAGA":
        st.warning("⚠️ PAGA data not found. Please run PAGA computation in the analysis step.")

    # Colormap selection in sidebar - conditional based on viz_type
    with st.sidebar:
        st.markdown("### Visualization Options")

        # Determine if this viz type uses discrete or continuous colormap
        uses_discrete = viz_type in ["Velocity embedding stream", "Velocity embedding grid", "Gene velocity", "PAGA"]
        uses_continuous = viz_type in ["Phase portraits", "Velocity confidence", "Velocity pseudotime"]

        if uses_discrete:
            colormap_discrete = st.selectbox(
                "Colormap (離散カラーマップ):",
                ["tab10", "Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
                 "tab20b", "tab20c", "Pastel1", "Pastel2", "Accent"],
                index=0,
                help="カテゴリカル変数用のカラーパレット"
            )
        else:
            # Set default for non-discrete viz types
            colormap_discrete = "tab10"

        if uses_continuous:
            colormap_continuous = st.selectbox(
                "Colormap (連続カラーマップ):",
                ["viridis", "plasma", "inferno", "magma", "cividis",
                 "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "Blues", "Greens", "Greys"],
                index=0,
                help="連続変数用のカラーマップ"
            )
        else:
            # Set default for non-continuous viz types
            colormap_continuous = "viridis"

        st.markdown("---")

        # Legend/cluster name display toggle
        show_legend = st.checkbox(
            "Show cluster names (legend)",
            value=True,
            help="クラスター名/カテゴリ名の凡例を表示するかどうか"
        )

        st.markdown("---")
        st.markdown("### Advanced Display Options")

        # Dot size and alpha controls
        dot_size = st.slider(
            "Dot size (細胞のサイズ)",
            min_value=1,
            max_value=200,
            value=10,
            step=5,
            help="散布図上の細胞ドットのサイズ"
        )

        dot_alpha = st.slider(
            "Dot alpha (細胞の透明度)",
            min_value=0.1,
            max_value=1.0,
            value=0.8,
            step=0.1,
            help="細胞ドットの透明度（0.1=透明、1.0=不透明）"
        )

    st.markdown("---")

    # Common parameters
    col1, col2 = st.columns(2)

    with col1:
        # Embedding basis selection - allows single or multiple selections
        # Find umap embedding if exists, otherwise use first
        default_basis = None
        for base in available_bases:
            if 'umap' in base.lower():
                default_basis = base
                break
        if default_basis is None and available_bases:
            default_basis = available_bases[0]

        basis_selections = st.multiselect(
            "Embedding basis (可視化する埋め込み空間)",
            available_bases,
            default=[default_basis] if default_basis else [],
            help="可視化に使用する埋め込み空間を選択。複数選択すると各embeddingに対してプロットを生成します"
        )

        if not basis_selections:
            st.warning("⚠️ Please select at least one embedding")
            basis_selections = [default_basis] if default_basis else []

    with col2:
        if viz_type in ["Velocity embedding stream", "Velocity embedding grid", "Gene velocity"]:
            # Cluster identity selection (not for PAGA)
            cluster_identity = st.selectbox(
                "Cluster identity (color by)",
                ["None"] + categorical_columns,
                index=1 if len(categorical_columns) > 0 else 0,
                help="メタデータから色付けする列を選択"
            )

    # Split by option - placed outside col2 but below cluster selection
    if viz_type in ["Velocity embedding stream", "Velocity embedding grid", "Gene velocity"]:
        split_by_enabled = st.checkbox(
            "Split plots by metadata",
            value=False,
            help="メタデータの各カテゴリごとに個別のプロットを作成"
        )

        if split_by_enabled:
            split_by_column = st.selectbox(
                "Select column to split by",
                categorical_columns,
                help="この列の各カテゴリごとにプロットを分割"
            )
        else:
            split_by_column = None
    else:
        split_by_enabled = False
        split_by_column = None

    # Continue with cluster order and color settings
    if viz_type in ["Velocity embedding stream", "Velocity embedding grid", "Gene velocity"]:
        # If cluster identity is selected, show cluster order option
        if cluster_identity != "None":
            with st.sidebar:
                st.markdown("---")
                st.markdown("### Cluster Settings")

                # Get unique clusters from current data
                cluster_list = sorted(adata.obs[cluster_identity].unique().astype(str).tolist())

                # Initialize sorted_order in session state
                if "velocity_sorted_order" not in st.session_state or \
                   st.session_state.get('velocity_prev_cluster_identity') != cluster_identity:
                    st.session_state.velocity_sorted_order = cluster_list
                    st.session_state.velocity_prev_cluster_identity = cluster_identity
                    sorted_order = cluster_list
                else:
                    sorted_order = st.session_state.get('velocity_sorted_order')
                    # Ensure sorted_order includes all current clusters
                    missing_clusters = [c for c in cluster_list if c not in sorted_order]
                    if missing_clusters:
                        sorted_order = sorted_order + missing_clusters
                        st.session_state.velocity_sorted_order = sorted_order

                # Color map management - always ensure all clusters are included
                if 'velocity_cluster_color_map' not in st.session_state or \
                   st.session_state.get('velocity_current_cmap', '') != colormap_discrete or \
                   st.session_state.get('velocity_prev_cluster_identity') != cluster_identity:
                    # Create color map with all current clusters
                    st.session_state.velocity_cluster_color_map = create_cell_color_mapping(sorted_order, colormap_discrete)
                    st.session_state.velocity_current_cmap = colormap_discrete
                else:
                    # Ensure existing color map has all clusters
                    existing_map = st.session_state.velocity_cluster_color_map
                    missing_in_map = [c for c in cluster_list if c not in existing_map]
                    if missing_in_map:
                        # Add missing clusters to color map
                        all_clusters = sorted_order
                        st.session_state.velocity_cluster_color_map = create_cell_color_mapping(all_clusters, colormap_discrete)

                # Cluster order sorting
                sort_clusters = st.checkbox("Change cluster order?")
                if sort_clusters:
                    with st.form("cluster_sorter"):
                        sorted_order_new = sort_items(sorted_order.copy())
                        submitted_sort = st.form_submit_button("Done sorting")
                    if submitted_sort:
                        st.session_state.velocity_sorted_order = sorted_order_new
                        current_cmap = st.session_state.get('velocity_current_cmap', colormap_discrete)
                        st.session_state.velocity_cluster_color_map = create_cell_color_mapping(sorted_order_new, current_cmap)

            color_by = cluster_identity
            palette = st.session_state.get('velocity_cluster_color_map')
        else:
            color_by = None
            palette = None
            sorted_order = None
    elif viz_type == "PAGA":
        # PAGA-specific cluster customization
        # For PAGA, we use the cluster groups from PAGA computation (no user selection)
        if has_paga:
            # Get the cluster key used by PAGA
            paga_groups_key = adata.uns['paga']['groups']

            # Display info about which cluster column PAGA is using
            st.info(f"💡 PAGA uses '{paga_groups_key}' clusters (auto-detected from PAGA computation)")

            # Get unique clusters from PAGA groups
            if paga_groups_key in adata.obs.columns:
                cluster_list = list(adata.obs[paga_groups_key].cat.categories)

                with st.sidebar:
                    st.markdown("---")
                    st.markdown("### PAGA Cluster Settings")

                    # Initialize sorted_order in session state for PAGA
                    if "velocity_paga_sorted_order" not in st.session_state or \
                       st.session_state.get('velocity_paga_prev_cluster') != paga_groups_key:
                        st.session_state.velocity_paga_sorted_order = cluster_list
                        st.session_state.velocity_paga_prev_cluster = paga_groups_key
                        sorted_order = cluster_list
                    else:
                        sorted_order = st.session_state.get('velocity_paga_sorted_order')
                        # Ensure sorted_order includes all current clusters
                        missing_clusters = [c for c in cluster_list if c not in sorted_order]
                        if missing_clusters:
                            sorted_order = sorted_order + missing_clusters
                            st.session_state.velocity_paga_sorted_order = sorted_order

                    # Color map management for PAGA
                    if 'velocity_paga_cluster_color_map' not in st.session_state or \
                       st.session_state.get('velocity_paga_current_cmap', '') != colormap_discrete or \
                       st.session_state.get('velocity_paga_prev_cluster') != paga_groups_key:
                        # Create color map with all current clusters
                        st.session_state.velocity_paga_cluster_color_map = create_cell_color_mapping(sorted_order, colormap_discrete)
                        st.session_state.velocity_paga_current_cmap = colormap_discrete
                    else:
                        # Ensure existing color map has all clusters
                        existing_map = st.session_state.velocity_paga_cluster_color_map
                        missing_in_map = [c for c in cluster_list if c not in existing_map]
                        if missing_in_map:
                            # Add missing clusters to color map
                            all_clusters = sorted_order
                            st.session_state.velocity_paga_cluster_color_map = create_cell_color_mapping(all_clusters, colormap_discrete)

                    # Cluster order sorting
                    sort_clusters = st.checkbox("Change cluster order?")
                    if sort_clusters:
                        with st.form("paga_cluster_sorter"):
                            sorted_order_new = sort_items(sorted_order.copy())
                            submitted_sort = st.form_submit_button("Done sorting")
                        if submitted_sort:
                            st.session_state.velocity_paga_sorted_order = sorted_order_new
                            current_cmap = st.session_state.get('velocity_paga_current_cmap', colormap_discrete)
                            st.session_state.velocity_paga_cluster_color_map = create_cell_color_mapping(sorted_order_new, current_cmap)

                color_by = paga_groups_key
                palette = st.session_state.get('velocity_paga_cluster_color_map')
            else:
                st.warning(f"⚠️ PAGA cluster column '{paga_groups_key}' not found in data")
                color_by = None
                palette = None
                sorted_order = None
        else:
            color_by = None
            palette = None
            sorted_order = None

    # Type-specific parameters
    if viz_type == "Velocity embedding stream":
        st.subheader("Stream plot parameters")
        col1, col2 = st.columns(2)
        with col1:
            arrow_size = st.slider("Arrow size", 0.5, 5.0, 2.0, 0.1)
        with col2:
            linewidth = st.slider("Line width", 0.5, 3.0, 1.5, 0.1)

        # Advanced options
        with st.expander("Advanced stream options"):
            col1, col2, col3 = st.columns(3)
            with col1:
                stream_density = st.slider("Stream density", 0.5, 3.0, 1.0, 0.1,
                                          help="ストリームラインの密度")
            with col2:
                stream_smooth = st.slider("Smoothing", 0.0, 1.0, 0.5, 0.1,
                                         help="ストリームラインの滑らかさ")
            with col3:
                stream_min_mass = st.slider("Min mass", 0.0, 10.0, 1.0, 0.5,
                                           help="矢印を表示する最小質量（閾値）")

    elif viz_type == "Velocity embedding grid":
        st.subheader("Grid plot parameters")
        col1, col2, col3 = st.columns(3)
        with col1:
            arrow_size = st.slider("Arrow size", 0.5, 5.0, 2.0, 0.1)
        with col2:
            density = st.slider("Grid density", 0.2, 2.0, 0.8, 0.1)
        with col3:
            arrow_length = st.slider("Arrow length", 1, 10, 3, 1)

        # Advanced options
        with st.expander("Advanced grid options"):
            col1, col2, col3 = st.columns(3)
            with col1:
                grid_scale = st.slider("Scale", 0.1, 2.0, 1.0, 0.1,
                                      help="矢印のスケール")
            with col2:
                grid_min_mass = st.slider("Min mass", 0.0, 10.0, 1.0, 0.5,
                                         help="矢印を表示する最小質量（閾値）")
            with col3:
                grid_alpha = st.slider("Alpha (transparency)", 0.1, 1.0, 1.0, 0.1,
                                      help="矢印の透明度")

    elif viz_type == "Phase portraits":
        st.subheader("Phase portrait parameters")

        # Gene selection
        n_genes = st.slider("Number of genes to show", 1, 20, 6, 1)

        gene_selection_method = st.radio(
            "Gene selection method",
            ["Top velocity genes", "Manual selection"],
            help="Top velocity genes: 最も高いvelocityを持つ遺伝子を自動選択"
        )

        if gene_selection_method == "Manual selection":
            selected_genes = st.multiselect(
                "Select genes",
                gene_list,
                max_selections=20,
                help="表示する遺伝子を選択（最大20個）"
            )
        else:
            selected_genes = None

    elif viz_type == "Gene velocity":
        st.subheader("Gene velocity parameters")
        selected_genes = st.multiselect(
            "Select genes",
            gene_list,
            max_selections=10,
            help="表示する遺伝子を選択（最大10個）"
        )

    elif viz_type == "PAGA":
        st.subheader("PAGA visualization parameters")

        paga_plot_type = st.radio(
            "Plot type",
            ["PAGA graph", "PAGA compare"],
            help="PAGA graph: クラスター間接続グラフのみ\nPAGA compare: PAGAグラフとUMAP/tSNEの比較"
        )

        # Common parameters for both plot types
        col1, col2, col3 = st.columns(3)
        with col1:
            paga_threshold = st.slider(
                "Edge threshold",
                0.0, 1.0, 0.01, 0.01,
                help="表示するエッジの閾値（低い値ほど多くの接続を表示）"
            )
        with col2:
            paga_node_size_scale = st.slider(
                "Node size scale",
                0.5, 5.0, 1.0, 0.5,
                help="ノードサイズのスケール"
            )
        with col3:
            paga_arrowsize = st.slider(
                "Arrow size",
                0, 30, 15, 1,
                help="方向性を示す矢印のサイズ（0で矢印なし）"
            )

    # Figure size
    st.subheader("Figure settings")
    col1, col2 = st.columns(2)
    with col1:
        fig_width = st.slider("Figure width", 4.0, 20.0, 10.0, 0.5)
    with col2:
        fig_height = st.slider("Figure height", 4.0, 20.0, 8.0, 0.5)

    # ========================================
    # Step 3: Generate visualization
    # ========================================
    if st.button("🎨 Generate Visualization", type="primary"):
        st.header("Step 3: Visualization")

        with st.spinner("Generating plots..."):
            try:
                # Set scvelo settings
                scv.settings.figdir = velocity_vis_temp_dir
                scv.settings.set_figure_params('scvelo', dpi=100, dpi_save=300,
                                              frameon=False, transparent=True)

                fig_files_png = []
                fig_files_pdf = []

                if viz_type == "Velocity embedding stream":
                    st.subheader("Velocity Embedding Stream")

                    # Determine split categories
                    if split_by_enabled and split_by_column:
                        split_categories = sorted(adata.obs[split_by_column].unique().astype(str).tolist())
                    else:
                        split_categories = [None]

                    # Generate plot for each selected embedding
                    for b in basis_selections:
                        for split_cat in split_categories:
                            # Subset data if splitting
                            if split_cat is not None:
                                adata_subset = adata[adata.obs[split_by_column] == split_cat].copy()
                                plot_title = f"Velocity Stream ({b}) - {split_by_column}: {split_cat}"
                                file_suffix = f"{b}_{split_by_column}_{split_cat}"
                            else:
                                adata_subset = adata
                                plot_title = f"Velocity Stream ({b})"
                                file_suffix = b

                            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                            scv.pl.velocity_embedding_stream(
                                adata_subset,
                                basis=b,
                                color=color_by if color_by not in ["None", None] else None,
                                palette=palette if palette else None,
                                size=dot_size,
                                alpha=dot_alpha,
                                arrow_size=arrow_size,
                                linewidth=linewidth,
                                density=stream_density,
                                smooth=stream_smooth,
                                min_mass=stream_min_mass,
                                legend_loc='right margin' if show_legend else 'none',
                                ax=ax,
                                show=False,
                                title=plot_title
                            )

                            # Save as PNG
                            fig_path_png = os.path.join(velocity_vis_temp_dir, f"velocity_stream_{file_suffix}.png")
                            fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                            fig_files_png.append(fig_path_png)

                            # Save as PDF (with error handling for non-finite values)
                            fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"velocity_stream_{file_suffix}.pdf")
                            try:
                                fig.savefig(fig_path_pdf, bbox_inches='tight')
                                fig_files_pdf.append(fig_path_pdf)
                            except ValueError as pdf_error:
                                if "finite" in str(pdf_error):
                                    st.warning(f"⚠️ PDF保存をスキップ ({file_suffix}): データに無限大/NaNが含まれています")
                                else:
                                    raise

                            st.pyplot(fig)
                            plt.close(fig)

                elif viz_type == "Velocity embedding grid":
                    st.subheader("Velocity Embedding Grid")

                    # Determine split categories
                    if split_by_enabled and split_by_column:
                        split_categories = sorted(adata.obs[split_by_column].unique().astype(str).tolist())
                    else:
                        split_categories = [None]

                    # Generate plot for each selected embedding
                    for b in basis_selections:
                        for split_cat in split_categories:
                            # Subset data if splitting
                            if split_cat is not None:
                                adata_subset = adata[adata.obs[split_by_column] == split_cat].copy()
                                plot_title = f"Velocity Grid ({b}) - {split_by_column}: {split_cat}"
                                file_suffix = f"{b}_{split_by_column}_{split_cat}"
                            else:
                                adata_subset = adata
                                plot_title = f"Velocity Grid ({b})"
                                file_suffix = b

                            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                            scv.pl.velocity_embedding_grid(
                                adata_subset,
                                basis=b,
                                color=color_by if color_by not in ["None", None] else None,
                                palette=palette if palette else None,
                                size=dot_size,
                                arrow_size=arrow_size,
                                density=density,
                                arrow_length=arrow_length,
                                scale=grid_scale,
                                min_mass=grid_min_mass,
                                alpha=grid_alpha,
                                legend_loc='right margin' if show_legend else 'none',
                                ax=ax,
                                show=False,
                                title=plot_title
                            )

                            # Save as PNG
                            fig_path_png = os.path.join(velocity_vis_temp_dir, f"velocity_grid_{file_suffix}.png")
                            fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                            fig_files_png.append(fig_path_png)

                            # Save as PDF (with error handling for non-finite values)
                            fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"velocity_grid_{file_suffix}.pdf")
                            try:
                                fig.savefig(fig_path_pdf, bbox_inches='tight')
                                fig_files_pdf.append(fig_path_pdf)
                            except ValueError as pdf_error:
                                if "finite" in str(pdf_error):
                                    st.warning(f"⚠️ PDF保存をスキップ ({file_suffix}): データに無限大/NaNが含まれています")
                                else:
                                    raise

                            st.pyplot(fig)
                            plt.close(fig)

                elif viz_type == "Phase portraits":
                    st.subheader("Phase Portraits")

                    if gene_selection_method == "Manual selection" and selected_genes:
                        genes_to_plot = selected_genes
                    else:
                        # Get top velocity genes
                        try:
                            # Try rank_velocity_genes (may fail with older scvelo + newer pandas)
                            scv.tl.rank_velocity_genes(adata, groupby=None, n_genes=n_genes)
                            genes_to_plot = adata.var.index[adata.var['velocity_genes'].argsort()[-n_genes:]].tolist()
                        except (AttributeError, ValueError) as e:
                            st.error(f"""
                            ❌ **rank_velocity_genes failed**

                            **Error:** {str(e)}

                            **原因:** scvelo v0.3.4 と pandas 2.0+ の互換性問題
                            - `AttributeError: property 'categories' of 'Categorical' object has no setter`

                            **解決方法:**
                            1. `rank_velocity_genes.py` を修正
                               - Line 173-177: `cat.reorder_categories()` → `cat.set_categories()`
                               - Line 191-197: `cat.categories = ...` → `cat.rename_categories()`
                            2. または "Manual selection" で遺伝子を手動選択してください
                            """)
                            genes_to_plot = []

                    if not genes_to_plot:
                        st.warning("⚠️ No genes selected!")
                    else:
                        st.write(f"**Showing {len(genes_to_plot)} genes:**")
                        st.write(", ".join(genes_to_plot))

                        # Tutorial method: scv.pl.velocity() for phase portraits with dynamical model lines
                        # This displays spliced vs unspliced with transcription/degradation lines
                        scv.pl.velocity(
                            adata,
                            var_names=genes_to_plot,
                            ncols=3,
                            figsize=(fig_width, fig_height),
                            show=False
                        )

                        # Get the current figure that scvelo created
                        fig = plt.gcf()

                        # Save as PNG
                        fig_path_png = os.path.join(velocity_vis_temp_dir, "phase_portraits.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        # Save as PDF (with error handling for non-finite values)
                        fig_path_pdf = os.path.join(velocity_vis_temp_dir, "phase_portraits.pdf")
                        try:
                            fig.savefig(fig_path_pdf, bbox_inches='tight')
                            fig_files_pdf.append(fig_path_pdf)
                        except ValueError as pdf_error:
                            if "finite" in str(pdf_error):
                                st.warning("⚠️ PDF保存をスキップ: データに無限大/NaNが含まれています")
                            else:
                                raise

                        st.pyplot(fig)
                        plt.close(fig)

                elif viz_type == "Velocity confidence":
                    st.subheader("Velocity Confidence")

                    if 'velocity_confidence' not in adata.obs:
                        st.warning("Computing velocity confidence...")
                        scv.tl.velocity_confidence(adata)

                    # Generate plot for each selected embedding
                    for b in basis_selections:
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        scv.pl.scatter(
                            adata,
                            basis=b,
                            color='velocity_confidence',
                            cmap='coolwarm',
                            size=dot_size,
                            alpha=dot_alpha,
                            ax=ax,
                            show=False,
                            title=f"Velocity Confidence ({b})"
                        )

                        # Save as PNG
                        fig_path_png = os.path.join(velocity_vis_temp_dir, f"velocity_confidence_{b}.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        # Save as PDF (with error handling for non-finite values)
                        fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"velocity_confidence_{b}.pdf")
                        try:
                            fig.savefig(fig_path_pdf, bbox_inches='tight')
                            fig_files_pdf.append(fig_path_pdf)
                        except ValueError as pdf_error:
                            if "finite" in str(pdf_error):
                                st.warning(f"⚠️ PDF保存をスキップ ({b}): データに無限大/NaNが含まれています")
                            else:
                                raise

                        st.pyplot(fig)
                        plt.close(fig)

                    # Show statistics
                    st.write("**Confidence statistics:**")
                    st.write(adata.obs['velocity_confidence'].describe())

                elif viz_type == "Velocity pseudotime":
                    st.subheader("Velocity Pseudotime")

                    if 'velocity_pseudotime' not in adata.obs:
                        st.warning("Computing velocity pseudotime...")
                        try:
                            scv.tl.velocity_pseudotime(adata)
                        except Exception as e:
                            st.error(f"❌ Could not compute velocity pseudotime: {str(e)}")
                            st.stop()

                    # Generate plot for each selected embedding
                    for b in basis_selections:
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        scv.pl.scatter(
                            adata,
                            basis=b,
                            color='velocity_pseudotime',
                            cmap='viridis',
                            size=dot_size,
                            alpha=dot_alpha,
                            ax=ax,
                            show=False,
                            title=f"Velocity Pseudotime ({b})"
                        )

                        # Save as PNG
                        fig_path_png = os.path.join(velocity_vis_temp_dir, f"velocity_pseudotime_{b}.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        # Save as PDF (with error handling for non-finite values)
                        fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"velocity_pseudotime_{b}.pdf")
                        try:
                            fig.savefig(fig_path_pdf, bbox_inches='tight')
                            fig_files_pdf.append(fig_path_pdf)
                        except ValueError as pdf_error:
                            if "finite" in str(pdf_error):
                                st.warning(f"⚠️ PDF保存をスキップ ({b}): データに無限大/NaNが含まれています")
                            else:
                                raise

                        st.pyplot(fig)
                        plt.close(fig)

                    # Show statistics
                    st.write("**Pseudotime statistics:**")
                    st.write(adata.obs['velocity_pseudotime'].describe())

                elif viz_type == "Gene velocity":
                    st.subheader("Gene-specific Velocity")

                    if not selected_genes:
                        st.warning("⚠️ Please select at least one gene!")
                    else:
                        # Determine split categories
                        if split_by_enabled and split_by_column:
                            split_categories = sorted(adata.obs[split_by_column].unique().astype(str).tolist())
                        else:
                            split_categories = [None]

                        # Generate plots for each selected embedding and gene
                        for b in basis_selections:
                            st.markdown(f"### Embedding: {b}")

                            for split_cat in split_categories:
                                # Subset data if splitting
                                if split_cat is not None:
                                    adata_subset = adata[adata.obs[split_by_column] == split_cat].copy()
                                    st.markdown(f"#### {split_by_column}: {split_cat}")
                                else:
                                    adata_subset = adata

                                for gene in selected_genes:
                                    if gene not in adata_subset.var_names:
                                        st.warning(f"⚠️ Gene '{gene}' not found in dataset")
                                        continue

                                    if split_cat is not None:
                                        st.markdown(f"**{gene}**")
                                    else:
                                        st.markdown(f"**{gene}**")

                                    # Use tutorial method: scv.pl.velocity automatically creates phase portrait + velocity embedding
                                    # Let scvelo create the figure automatically
                                    scv.pl.velocity(
                                        adata_subset,
                                        var_names=[gene],
                                        basis=b,
                                        color=color_by if color_by not in ["None", None] else None,
                                        palette=palette if palette else None,
                                        size=dot_size,
                                        alpha=dot_alpha,
                                        figsize=(fig_width, fig_height//2),
                                        legend_loc='best' if show_legend else 'none',
                                        show=False
                                    )

                                    # Get the current figure that scvelo created
                                    fig = plt.gcf()

                                    # Adjust layout
                                    plt.tight_layout()

                                    # Save as PNG
                                    if split_cat is not None:
                                        file_suffix = f"{gene}_{b}_{split_by_column}_{split_cat}"
                                    else:
                                        file_suffix = f"{gene}_{b}"

                                    fig_path_png = os.path.join(velocity_vis_temp_dir, f"gene_velocity_{file_suffix}.png")
                                    fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                                    fig_files_png.append(fig_path_png)

                                    # Save as PDF (with error handling for non-finite values)
                                    fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"gene_velocity_{file_suffix}.pdf")
                                    try:
                                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                                        fig_files_pdf.append(fig_path_pdf)
                                    except ValueError as pdf_error:
                                        if "finite" in str(pdf_error):
                                            st.warning(f"⚠️ PDF保存をスキップ ({file_suffix}): データに無限大/NaNが含まれています")
                                        else:
                                            raise

                                st.pyplot(fig)
                                plt.close(fig)

                elif viz_type == "PAGA":
                    st.subheader("PAGA Visualization")

                    if 'paga' not in adata.uns:
                        st.error("❌ PAGA data not found. Please run PAGA computation in the velocity analysis step.")
                    else:
                        # color_by and palette are already set correctly in configuration section
                        # They are guaranteed to use PAGA's groups if PAGA is selected

                        if paga_plot_type == "PAGA graph":
                            st.markdown("### PAGA Graph")

                            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                            # Plot PAGA graph with directed arrows
                            scv.pl.paga(
                                adata,
                                color=color_by if color_by else None,
                                palette=palette if palette else None,
                                threshold=paga_threshold,
                                node_size_scale=paga_node_size_scale,
                                arrowsize=paga_arrowsize,
                                legend_loc='on data' if show_legend else 'none',
                                ax=ax,
                                show=False,
                                title="PAGA Graph"
                            )

                            # Save as PNG
                            fig_path_png = os.path.join(velocity_vis_temp_dir, "paga_graph.png")
                            fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                            fig_files_png.append(fig_path_png)

                            # Save as PDF (with error handling for non-finite values)
                            fig_path_pdf = os.path.join(velocity_vis_temp_dir, "paga_graph.pdf")
                            try:
                                fig.savefig(fig_path_pdf, bbox_inches='tight')
                                fig_files_pdf.append(fig_path_pdf)
                            except ValueError as pdf_error:
                                if "finite" in str(pdf_error):
                                    st.warning("⚠️ PDF保存をスキップ: データに無限大/NaNが含まれています")
                                else:
                                    raise

                            st.pyplot(fig)
                            plt.close(fig)

                        else:  # PAGA compare
                            st.markdown("### PAGA Compare")

                            # Generate plot for each selected embedding
                            for b in basis_selections:
                                # Let scVelo create the figure - don't pre-create axes
                                # This avoids internal axes handling issues in scvelo
                                plt.figure(figsize=(fig_width*2, fig_height))

                                scv.pl.paga(
                                    adata,
                                    basis=b,
                                    color=color_by if color_by else None,
                                    palette=palette if palette else None,
                                    threshold=paga_threshold,
                                    node_size_scale=paga_node_size_scale,
                                    arrowsize=paga_arrowsize,
                                    size=200,
                                    legend_loc='on data' if show_legend else 'none',
                                    show=False
                                )

                                # Get the current figure
                                fig = plt.gcf()

                                # Adjust layout
                                plt.tight_layout()

                                # Save as PNG
                                fig_path_png = os.path.join(velocity_vis_temp_dir, f"paga_compare_{b}.png")
                                fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                                fig_files_png.append(fig_path_png)

                                # Save as PDF (with error handling for non-finite values)
                                fig_path_pdf = os.path.join(velocity_vis_temp_dir, f"paga_compare_{b}.pdf")
                                try:
                                    fig.savefig(fig_path_pdf, bbox_inches='tight')
                                    fig_files_pdf.append(fig_path_pdf)
                                except ValueError as pdf_error:
                                    if "finite" in str(pdf_error):
                                        st.warning(f"⚠️ PDF保存をスキップ ({b}): データに無限大/NaNが含まれています")
                                    else:
                                        raise

                                st.pyplot(fig)
                                plt.close(fig)

                            # Show connectivity information
                            st.write("**PAGA Connectivity:**")
                            if 'connectivities' in adata.uns['paga']:
                                conn_df = pd.DataFrame(
                                    adata.uns['paga']['connectivities'].toarray(),
                                    index=adata.obs[adata.uns['paga']['groups']].cat.categories,
                                    columns=adata.obs[adata.uns['paga']['groups']].cat.categories
                                )
                                # Round to 3 decimal places (avoids jinja2 dependency)
                                st.dataframe(conn_df.round(3))

                # Download section
                if fig_files_png:
                    st.success(f"✅ Generated {len(fig_files_png)} plot(s)")

                    st.subheader("Download plots")

                    # Create two columns for PNG and PDF downloads
                    col_png, col_pdf = st.columns(2)

                    with col_png:
                        st.markdown("**PNG Format**")
                        for fig_file in fig_files_png:
                            file_name = os.path.basename(fig_file)
                            with open(fig_file, "rb") as f:
                                st.download_button(
                                    label=f"⬇️ {file_name}",
                                    data=f.read(),
                                    file_name=file_name,
                                    mime="image/png",
                                    key=f"png_{file_name}"
                                )

                    with col_pdf:
                        st.markdown("**PDF Format**")
                        for fig_file in fig_files_pdf:
                            file_name = os.path.basename(fig_file)
                            with open(fig_file, "rb") as f:
                                st.download_button(
                                    label=f"⬇️ {file_name}",
                                    data=f.read(),
                                    file_name=file_name,
                                    mime="application/pdf",
                                    key=f"pdf_{file_name}"
                                )

            except Exception as e:
                st.error(f"❌ Error during visualization: {str(e)}")
                st.exception(e)

else:
    st.info("👆 scVelo解析結果のh5adファイルをアップロードして開始してください")
