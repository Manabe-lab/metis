"""
Pseudotime Gene Expression Visualization
Visualize gene expression trends along pseudotime
"""

import streamlit as st
import scanpy as sc
import scvelo as scv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
import re
from scipy.interpolate import UnivariateSpline
from scipy.ndimage import zoom
from helper_func import clear_old_directories, clear_old_files
from streamlit_sortables import sort_items
from mpl_toolkits.axes_grid1 import make_axes_locatable

@st.cache_data
def create_cell_color_mapping(cell_list, palette_name):
    """
    クラスター/細胞タイプと色の一貫したマッピングを作成する関数

    Parameters
    ----------
    cell_list : list
        クラスター名/細胞タイプ名のリスト
    palette_name : str
        使用する離散カラーパレット名

    Returns
    -------
    dict
        クラスター/細胞タイプ名をキー、RGBタプルを値とする辞書
    """
    n_cells = len(cell_list)

    # 指定されたパレットから色を取得
    base_palette = sns.color_palette(palette_name)

    # パレットの色数が足りない場合は拡張
    if n_cells <= len(base_palette):
        colors = base_palette[:n_cells]
    else:
        colors = sns.color_palette(palette_name, n_colors=n_cells)

    # 辞書を作成
    color_dict = {cell: color for cell, color in zip(cell_list, colors)}

    return color_dict

def add_cluster_annotation_bar(ax, pseudotime_sorted, cluster_sorted, bins, cluster_color_dict, reverse_pseudotime=False, height_inches=0.15):
    """
    プロット下部（x軸ラベルの下）にクラスターアノテーションバーを追加する関数

    axes_dividerを使用してx軸の下にアノテーションバーを配置する

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        プロットのaxesオブジェクト
    pseudotime_sorted : np.ndarray
        ソート済みpseudotime値
    cluster_sorted : np.ndarray
        ソート済みクラスター割り当て
    bins : np.ndarray
        ビンの境界
    cluster_color_dict : dict
        クラスター名をキー、色を値とする辞書
    reverse_pseudotime : bool
        pseudotimeの方向を反転するか
    height_inches : float
        アノテーションバーの高さ（インチ単位）
    """
    # Calculate dominant cluster for each bin
    bin_indices = np.digitize(pseudotime_sorted, bins)
    bin_cluster_labels = []

    for i in range(1, len(bins)):
        mask = bin_indices == i
        if np.sum(mask) > 0:
            bin_cluster_values = cluster_sorted[mask]
            unique, counts = np.unique(bin_cluster_values, return_counts=True)
            dominant_cluster = str(unique[np.argmax(counts)])
            bin_cluster_labels.append(dominant_cluster)
        else:
            bin_cluster_labels.append(None)

    # Reverse if needed
    if reverse_pseudotime:
        bin_cluster_labels = bin_cluster_labels[::-1]

    # Create divider for existing axes
    divider = make_axes_locatable(ax)

    # Add new axes below the x-axis with specified height
    ax_anno = divider.append_axes("bottom", size=height_inches, pad=0.3, sharex=ax)

    # Plot colored rectangles
    for i in range(len(bins)-1):
        cluster = bin_cluster_labels[i]
        if cluster is not None:
            color = cluster_color_dict.get(cluster, (0.5, 0.5, 0.5))
            ax_anno.barh(0, bins[i+1] - bins[i], left=bins[i], height=1,
                        color=color, edgecolor='white', linewidth=0.5)

    # Set limits and remove axes decorations
    ax_anno.set_xlim(pseudotime_sorted.min(), pseudotime_sorted.max())
    ax_anno.set_ylim(-0.5, 0.5)
    ax_anno.set_xticks([])
    ax_anno.set_yticks([])
    ax_anno.spines['top'].set_visible(False)
    ax_anno.spines['bottom'].set_visible(False)
    ax_anno.spines['left'].set_visible(False)
    ax_anno.spines['right'].set_visible(False)

    if reverse_pseudotime:
        ax_anno.invert_xaxis()

    return ax_anno

st.set_page_config(page_title="Pseudotime Gene Expression", page_icon="📈", layout="wide")

st.title("📈 Pseudotime Gene Expression Visualization")
st.markdown("""
Pseudotimeに沿った遺伝子発現トレンドを可視化します。

### 可視化の種類
1. **Gene expression - Line plot**: Pseudotimeに沿った遺伝子発現の変化を線グラフで表示
   - **Individual genes**: 各遺伝子ごとに個別のプロット
   - **Overlay genes**: 全遺伝子を1つのグラフに重ね合わせ
2. **Gene expression - Heatmap**: 複数遺伝子の発現パターンをヒートマップで表示
3. **Cluster density**: Pseudotimeに沿ったクラスター構成の変化を可視化
   - **Line plot**: 各クラスターの密度変化
   - **Stacked area**: クラスター構成の推移（100%積み上げ）
   - **Heatmap**: クラスター密度の2次元表示

### 使用できるPseudotime
- **Velocity pseudotime** (scVelo) - RNA velocity解析から計算される疑似時間
- **Latent time** (scVelo dynamical model) - 動的モデルに基づく潜在時間
- **DPT pseudotime** (Diffusion pseudotime) - 拡散疑似時間
- **Fate probabilities** (CellRank) - 各終末状態への運命確率（`fate_prob_{state}`列）
  - 横軸 = 0（確率低）〜 1（確率高）
  - 特定の細胞運命への分化過程の遺伝子発現変化を可視化
- その他のカスタムpseudotime列

### 参考
- [CellRank Gene Trends](https://cellrank.readthedocs.io/en/stable/api.html#plotting)
""")

# Initialize session state
if "pseudotime_vis_temp_dir" not in st.session_state:
    pseudotime_vis_temp_dir = os.path.join("temp", f"pseudotime_vis_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(pseudotime_vis_temp_dir, exist_ok=True)
    st.session_state.pseudotime_vis_temp_dir = pseudotime_vis_temp_dir
else:
    pseudotime_vis_temp_dir = st.session_state.pseudotime_vis_temp_dir

# ========================================
# Step 1: Upload file
# ========================================
st.header("Step 1: Upload h5ad file")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file",
    type=['h5ad'],
    key="pseudotime_vis_h5ad_upload",
    help="scVelo/CellRank解析済みのh5adファイル"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load data
    if ("pseudotime_vis_adata" not in st.session_state or
        st.session_state.get("pseudotime_vis_uploaded_file") != uploaded_h5ad.name):

        with st.spinner("Loading data..."):
            # Save uploaded file temporarily
            temp_h5ad_path = os.path.join(pseudotime_vis_temp_dir, "input.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            # Read data
            try:
                adata = sc.read_h5ad(temp_h5ad_path)
            except ValueError as e:
                if "Data must be 1-dimensional" in str(e):
                    st.warning("⚠️ Detected malformed metadata. Attempting to fix...")
                    import anndata as ad
                    adata_backed = ad.read_h5ad(temp_h5ad_path, backed='r')
                    adata = ad.AnnData(
                        X=adata_backed.X[:],
                        obs=adata_backed.obs.copy(),
                        var=adata_backed.var.copy(),
                        uns=adata_backed.uns.copy() if len(adata_backed.uns) > 0 else {},
                        obsm=adata_backed.obsm.copy() if len(adata_backed.obsm) > 0 else {},
                        layers={k: adata_backed.layers[k][:] for k in adata_backed.layers.keys()} if len(adata_backed.layers) > 0 else {}
                    )
                    for col in adata.obs.columns:
                        if hasattr(adata.obs[col], 'shape') and len(adata.obs[col].shape) > 1:
                            if adata.obs[col].shape[1] == 1:
                                adata.obs[col] = adata.obs[col].values.flatten()
                            else:
                                adata.obs = adata.obs.drop(columns=[col])
                    st.success("✓ Successfully cleaned and loaded h5ad file")
                else:
                    raise

            st.session_state.pseudotime_vis_adata = adata
            st.session_state.pseudotime_vis_uploaded_file = uploaded_h5ad.name

            st.info(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    adata = st.session_state.pseudotime_vis_adata

    # Get available pseudotime columns
    pseudotime_cols = []
    priority_names = ['velocity_pseudotime', 'latent_time', 'dpt_pseudotime', 'pseudotime']

    # Check for priority names first
    for name in priority_names:
        if name in adata.obs.columns:
            pseudotime_cols.append(name)

    # Add other numeric columns that might be pseudotime
    for col in adata.obs.columns:
        if pd.api.types.is_numeric_dtype(adata.obs[col]) and col not in pseudotime_cols:
            col_lower = col.lower()
            # Include time-based, fate probability, and similar columns
            if ('time' in col_lower or 'pseudotime' in col_lower or 'latent' in col_lower or
                ('fate' in col_lower and 'prob' in col_lower)):
                pseudotime_cols.append(col)

    if not pseudotime_cols:
        st.error("""
        ❌ **No pseudotime information found!**

        This file does not contain pseudotime data. Please run one of the following:
        - scVelo velocity analysis (generates 'velocity_pseudotime')
        - scVelo dynamical model (generates 'latent_time')
        - Scanpy DPT (generates 'dpt_pseudotime')
        """)
        st.stop()

    # Get gene list
    gene_list = adata.var_names.tolist()

    # Get categorical columns for cluster identity
    categorical_columns = []
    for col in adata.obs.columns:
        if adata.obs[col].dtype.name == 'category' or (adata.obs[col].dtype == 'object' and adata.obs[col].nunique() < 50):
            categorical_columns.append(col)

    # ========================================
    # Step 2: Configure visualization
    # ========================================
    st.header("Step 2: Configure visualization")

    # Main analysis type selection
    analysis_type = st.radio(
        "Analysis type",
        ["Gene expression", "Cluster density"],
        help="遺伝子発現 or クラスター密度を選択"
    )

    st.markdown("---")

    # Visualization type selection (depends on analysis type)
    if analysis_type == "Gene expression":
        viz_type = st.radio(
            "Visualization type",
            ["Line plot", "Heatmap"],
            help="線グラフ or ヒートマップを選択"
        )
    else:  # Cluster density
        viz_type = st.radio(
            "Visualization type",
            ["Line plot", "Stacked area", "Heatmap"],
            help="線グラフ、積み上げエリア、またはヒートマップを選択"
        )

    st.markdown("---")

    # Common parameters
    st.subheader("Pseudotime selection")

    selected_pseudotime = st.selectbox(
        "Select pseudotime",
        pseudotime_cols,
        help="使用するpseudotimeを選択"
    )

    # Show pseudotime statistics
    with st.expander("📊 Pseudotime statistics"):
        pt_stats = adata.obs[selected_pseudotime].describe()
        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Min", f"{pt_stats['min']:.3f}")
        with col2:
            st.metric("Mean", f"{pt_stats['mean']:.3f}")
        with col3:
            st.metric("Max", f"{pt_stats['max']:.3f}")
        with col4:
            st.metric("Valid cells", f"{adata.obs[selected_pseudotime].notna().sum()}")

    # Sidebar options for pseudotime
    st.sidebar.markdown("### Pseudotime Options")
    reverse_pseudotime = st.sidebar.checkbox(
        "Reverse pseudotime direction",
        value=False,
        help="Pseudotimeの方向を逆転（終点→始点）"
    )

    log1p_pseudotime = st.sidebar.checkbox(
        "Log1p transform pseudotime (x-axis)",
        value=False,
        help="Pseudotimeをlog1p変換してx軸に表示（log(1+x)）。Pseudotimeの分布が歪んでいる場合に有効"
    )

    st.markdown("---")

    # Selection based on analysis type
    if analysis_type == "Gene expression":
        st.subheader("Gene selection")

        gene_input_method = st.radio(
            "Gene input method",
            ["Multi-select", "Text input"],
            help="遺伝子の選択方法"
        )

        if gene_input_method == "Multi-select":
            selected_genes = st.multiselect(
                "Select genes",
                gene_list,
                max_selections=50,
                help="表示する遺伝子を選択（最大50個）"
            )
        else:
            gene_text_input = st.text_area(
                "Enter gene names",
                height=150,
                help="遺伝子名をスペース、コンマ、タブ、または改行で区切って入力"
            )

            if gene_text_input:
                # Parse gene names from text
                # Split by comma, space, tab, or newline
                genes_raw = re.split(r'[,\s\t\n]+', gene_text_input.strip())
                # Remove empty strings
                genes_raw = [g.strip() for g in genes_raw if g.strip()]

                # Create case-insensitive mapping (lowercase -> original gene name)
                gene_map_lower = {g.lower(): g for g in gene_list}

                # Check which genes are in the dataset (case-insensitive)
                selected_genes = []
                missing_genes = []
                for g in genes_raw:
                    g_lower = g.lower()
                    if g_lower in gene_map_lower:
                        selected_genes.append(gene_map_lower[g_lower])
                    else:
                        missing_genes.append(g)

                st.info(f"✓ Found {len(selected_genes)} genes in dataset")

                if missing_genes:
                    st.warning(f"⚠️ {len(missing_genes)} genes not found: {', '.join(missing_genes[:10])}" +
                              (f" ... and {len(missing_genes)-10} more" if len(missing_genes) > 10 else ""))
            else:
                selected_genes = []

        if not selected_genes:
            st.warning("⚠️ Please select at least one gene")

    else:  # Cluster density
        st.subheader("Cluster selection")

        if not categorical_columns:
            st.error("❌ No categorical columns found in metadata!")
            st.stop()

        cluster_key = st.selectbox(
            "Select cluster column",
            categorical_columns,
            help="クラスター情報を含むmetadata列を選択"
        )

        # Get unique clusters
        clusters = adata.obs[cluster_key].unique()
        if hasattr(clusters, 'categories'):
            clusters = clusters.categories.tolist()
        else:
            clusters = sorted([str(c) for c in clusters if pd.notna(c)])

        st.info(f"✓ Found {len(clusters)} clusters: {', '.join(map(str, clusters[:10]))}" +
               (f" ... and {len(clusters)-10} more" if len(clusters) > 10 else ""))

        # For cluster density, we'll use all clusters by default
        selected_genes = None  # Not used for cluster density

    st.markdown("---")

    # Type-specific parameters
    if analysis_type == "Gene expression" and viz_type == "Line plot":
        st.subheader("Line plot parameters")

        line_plot_mode = st.radio(
            "Plot mode",
            ["Individual genes", "Overlay genes"],
            help="Individual: 各遺伝子ごとに個別のプロット\nOverlay: 全遺伝子を1つのグラフに重ね合わせ"
        )

        col1, col2, col3 = st.columns(3)
        with col1:
            n_bins = st.slider(
                "Number of bins",
                5, 100, 10, 5,
                help="Pseudotimeをビン分割する数（少ない方が滑らかになります）"
            )
        with col2:
            show_scatter = st.checkbox(
                "Show scatter points",
                value=False,
                help="個々の細胞を散布図として表示"
            )
        with col3:
            smoothing_mode = st.selectbox(
                "Smoothing mode",
                ["Strong (smooth)", "Interpolation (detailed)"],
                index=0,
                help="スムージング強度\nStrong: 平均値と信頼区間を強力にスムージング（推奨）\nInterpolation: ビン平均を完全に補間（詳細だがノイジー）"
            )

        col4, col5, col6 = st.columns(3)
        with col4:
            show_ci = st.checkbox(
                "Show confidence interval",
                value=True,
                help="95%信頼区間を表示"
            )
        with col5:
            use_zscore_line = st.checkbox(
                "Z-score normalization",
                value=False,
                help="各遺伝子をZ-score正規化（平均0、標準偏差1）"
            )
        with col6:
            line_color = st.selectbox(
                "Color scheme",
                ["Red (#E64B35)", "Blue (#4DBBD5)", "Green (#00A087)", "Purple (#3C5488)", "Orange (#F39B7F)"],
                index=0,
                help="線の色（Individual genes用）"
            )

        # Cluster split option
        st.markdown("---")
        split_by_cluster = st.checkbox(
            "Split by cluster",
            value=False,
            help="クラスター別に遺伝子発現トレンドを分けて表示"
        )

        if split_by_cluster:
            if not categorical_columns:
                st.error("❌ No categorical columns found in metadata!")
                split_by_cluster = False
            else:
                col_cluster, col_palette = st.columns(2)
                with col_cluster:
                    split_cluster_key = st.selectbox(
                        "Select cluster column",
                        categorical_columns,
                        help="クラスター情報を含むmetadata列を選択"
                    )
                with col_palette:
                    split_cluster_palette = st.selectbox(
                        "Cluster color palette",
                        ["husl", "tab10", "Set2", "Set3", "Paired", "Dark2"],
                        index=0,
                        help="クラスターの色分けパレット"
                    )

        # Cluster annotation bar option
        st.markdown("---")
        st.markdown("### Cluster Annotation Bar")
        show_cluster_annotation_line = st.checkbox(
            "Show cluster annotation bar",
            value=False,
            help="グラフ下部に各pseudotime範囲の主要クラスターをカラーバーで表示",
            key="line_plot_cluster_annotation"
        )

        if show_cluster_annotation_line:
            if not categorical_columns:
                st.error("❌ No categorical columns found in metadata!")
                show_cluster_annotation_line = False
            else:
                col_ann_cluster, col_ann_palette = st.columns(2)
                with col_ann_cluster:
                    annotation_cluster_key = st.selectbox(
                        "Cluster identity for annotation",
                        categorical_columns,
                        help="アノテーション用のクラスター列を選択",
                        key="line_annotation_cluster_key"
                    )
                with col_ann_palette:
                    annotation_cluster_palette = st.selectbox(
                        "Annotation color palette",
                        ["tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Paired"],
                        index=0,
                        help="クラスターアノテーションの色分けパレット",
                        key="line_annotation_palette"
                    )

                # Get unique clusters for annotation
                annotation_clusters = adata.obs[annotation_cluster_key].unique()
                if hasattr(annotation_clusters, 'categories'):
                    annotation_cluster_list = annotation_clusters.categories.tolist()
                else:
                    annotation_cluster_list = sorted([str(c) for c in annotation_clusters if pd.notna(c)])

                # Create color mapping for annotation
                annotation_cluster_colors = sns.color_palette(annotation_cluster_palette, len(annotation_cluster_list))
                annotation_cluster_color_dict = {str(cluster): color for cluster, color in zip(annotation_cluster_list, annotation_cluster_colors)}

        if line_plot_mode == "Individual genes":
            ncols = st.slider(
                "Number of columns",
                1, 5, 3, 1,
                help="サブプロットの列数"
            )
        else:
            col_legend, col_palette = st.columns(2)
            with col_legend:
                show_legend = st.checkbox(
                    "Show legend",
                    value=True,
                    help="遺伝子名の凡例を表示"
                )
            with col_palette:
                overlay_palette = st.selectbox(
                    "Color palette",
                    ["husl", "tab10", "Set2", "Set3", "Paired", "Dark2"],
                    index=0,
                    help="遺伝子の色分けパレット"
                )

    elif analysis_type == "Gene expression" and viz_type == "Heatmap":
        st.subheader("Heatmap parameters")

        col1, col2, col3 = st.columns(3)
        with col1:
            n_bins = st.slider(
                "Number of bins",
                10, 200, 100, 10,
                help="Pseudotimeをビン分割する数"
            )
        with col2:
            use_clustermap = st.checkbox(
                "Hierarchical clustering",
                value=True,
                help="Dendrogramを表示する階層的クラスタリングを使用"
            )
        with col3:
            use_interpolation = st.checkbox(
                "Smooth interpolation",
                value=False,
                help="ビン間を視覚的に補間（bilinear）"
            )

        if use_clustermap:
            col4, col5 = st.columns(2)
            with col4:
                cluster_rows = st.checkbox(
                    "Cluster genes (rows)",
                    value=True,
                    help="遺伝子（行）を階層的クラスタリング"
                )
            with col5:
                linkage_method = st.selectbox(
                    "Linkage method",
                    ["ward", "average", "complete", "single"],
                    index=0,
                    help="階層的クラスタリングの連結法"
                )

            # Pseudotime columns are never clustered (time order must be preserved)
            cluster_cols = False
        else:
            cluster_genes = st.checkbox(
                "Cluster genes",
                value=True,
                help="遺伝子を発現パターンで階層的クラスタリング（dendrogramなし）"
            )

        col6, col7, col8 = st.columns(3)
        with col6:
            use_zscore = st.checkbox(
                "Z-score normalization",
                value=True,
                help="各遺伝子をZ-score正規化（平均0、標準偏差1）"
            )
        with col7:
            heatmap_cmap = st.selectbox(
                "Colormap",
                ["viridis", "plasma", "inferno", "magma", "RdYlBu_r", "coolwarm", "bwr"],
                index=0,
                help="ヒートマップのカラーマップ"
            )
        with col8:
            show_gene_names = st.checkbox(
                "Show gene names",
                value=True,
                help="遺伝子名をY軸に表示"
            )

        # Cluster annotation
        st.markdown("### Cluster Annotation (上部カラーバー)")
        show_cluster_annotation = st.checkbox(
            "Show cluster annotation",
            value=False,
            help="Heatmap上部に各binの主要クラスターを表示"
        )

        if show_cluster_annotation:
            col7, col8 = st.columns(2)
            with col7:
                cluster_annotation_key = st.selectbox(
                    "Cluster identity",
                    categorical_columns,
                    help="アノテーション用のクラスター列を選択"
                )
            with col8:
                cluster_palette = st.selectbox(
                    "Cluster colormap",
                    ["tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Paired"],
                    index=0,
                    help="クラスターの色分けパレット"
                )

            # Get unique clusters for annotation
            clusters_for_annotation = adata.obs[cluster_annotation_key].unique()
            if hasattr(clusters_for_annotation, 'categories'):
                cluster_list = clusters_for_annotation.categories.tolist()
            else:
                cluster_list = sorted([str(c) for c in clusters_for_annotation if pd.notna(c)])

            # Initialize cluster order and color mapping in session state
            if "pseudotime_heatmap_cluster_order" not in st.session_state or \
               st.session_state.get("pseudotime_heatmap_cluster_key") != cluster_annotation_key:
                st.session_state.pseudotime_heatmap_cluster_order = cluster_list
                st.session_state.pseudotime_heatmap_cluster_key = cluster_annotation_key
                sorted_order = cluster_list
            else:
                sorted_order = st.session_state.get('pseudotime_heatmap_cluster_order')
                # Ensure sorted_order includes all current clusters
                missing_clusters = [c for c in cluster_list if c not in sorted_order]
                if missing_clusters:
                    sorted_order = sorted_order + missing_clusters
                    st.session_state.pseudotime_heatmap_cluster_order = sorted_order

            # Initialize or update cluster color mapping
            if 'pseudotime_heatmap_cluster_color_map' not in st.session_state or \
               st.session_state.get("pseudotime_heatmap_cluster_palette") != cluster_palette or \
               st.session_state.get("pseudotime_heatmap_cluster_key") != cluster_annotation_key:
                st.session_state.pseudotime_heatmap_cluster_color_map = create_cell_color_mapping(sorted_order, cluster_palette)
                st.session_state.pseudotime_heatmap_cluster_palette = cluster_palette
            else:
                # Check if there are new clusters
                existing_map = st.session_state.pseudotime_heatmap_cluster_color_map
                if set(sorted_order) != set(existing_map.keys()):
                    # Add missing clusters to color map
                    st.session_state.pseudotime_heatmap_cluster_color_map = create_cell_color_mapping(sorted_order, cluster_palette)

            # Cluster order sorting
            sort_clusters = st.checkbox("Change cluster order?")

            if sort_clusters:
                sorted_order_new = sort_items(sorted_order.copy())

                # Update session state when order changes
                st.session_state.pseudotime_heatmap_cluster_order = sorted_order_new
                current_cmap = st.session_state.get("pseudotime_heatmap_cluster_palette", cluster_palette)
                st.session_state.pseudotime_heatmap_cluster_color_map = create_cell_color_mapping(sorted_order_new, current_cmap)

                st.info(f"✓ Cluster order updated: {', '.join(sorted_order_new)}")

            cluster_color_map = st.session_state.get('pseudotime_heatmap_cluster_color_map')
        else:
            cluster_annotation_key = None
            cluster_color_map = None

    elif analysis_type == "Cluster density":
        st.subheader("Cluster density parameters")

        col1, col2 = st.columns(2)
        with col1:
            n_bins = st.slider(
                "Number of bins",
                5, 100, 10, 5,
                help="Pseudotimeをビン分割する数（少ない方が滑らかになります）"
            )
        with col2:
            density_type = st.selectbox(
                "Density calculation",
                ["Proportion (%)", "Cell count"],
                help="Proportion: 各ビンにおける各クラスターの割合\nCell count: 各ビンにおける細胞数"
            )

        # Normalization option for cluster size
        normalize_cluster_size = st.checkbox(
            "Normalize by cluster size",
            value=False,
            help="各クラスターの細胞総数を同じにして比較。クラスターサイズの違いによる影響を除外し、相対的な密度変化を可視化します。"
        )

        if viz_type == "Heatmap":
            heatmap_cmap = st.selectbox(
                "Colormap",
                ["viridis", "plasma", "inferno", "magma", "YlOrRd", "Blues", "Greens"],
                index=0,
                help="ヒートマップのカラーマップ"
            )

        # Cluster annotation bar option for Line plot and Stacked area
        if viz_type in ["Line plot", "Stacked area"]:
            st.markdown("---")
            st.markdown("### Cluster Annotation Bar")
            show_cluster_annotation_density = st.checkbox(
                "Show cluster annotation bar",
                value=False,
                help="グラフ下部に各pseudotime範囲の主要クラスターをカラーバーで表示",
                key="density_plot_cluster_annotation"
            )

            if show_cluster_annotation_density:
                col_ann_cluster, col_ann_palette = st.columns(2)
                with col_ann_cluster:
                    annotation_cluster_key_density = st.selectbox(
                        "Cluster identity for annotation",
                        categorical_columns,
                        help="アノテーション用のクラスター列を選択（Cluster densityで使用中のものと異なる列も選択可能）",
                        key="density_annotation_cluster_key"
                    )
                with col_ann_palette:
                    annotation_cluster_palette_density = st.selectbox(
                        "Annotation color palette",
                        ["tab20", "tab20b", "tab20c", "Set1", "Set2", "Set3", "Paired"],
                        index=0,
                        help="クラスターアノテーションの色分けパレット",
                        key="density_annotation_palette"
                    )

                # Get unique clusters for annotation
                annotation_clusters_density = adata.obs[annotation_cluster_key_density].unique()
                if hasattr(annotation_clusters_density, 'categories'):
                    annotation_cluster_list_density = annotation_clusters_density.categories.tolist()
                else:
                    annotation_cluster_list_density = sorted([str(c) for c in annotation_clusters_density if pd.notna(c)])

                # Create color mapping for annotation
                annotation_cluster_colors_density = sns.color_palette(annotation_cluster_palette_density, len(annotation_cluster_list_density))
                annotation_cluster_color_dict_density = {str(cluster): color for cluster, color in zip(annotation_cluster_list_density, annotation_cluster_colors_density)}

    # Figure size
    st.subheader("Figure settings")
    col1, col2 = st.columns(2)
    with col1:
        fig_width = st.slider("Figure width", 4, 20, 12, 1)
    with col2:
        fig_height = st.slider("Figure height", 2, 30, 8, 1)

    # ========================================
    # Step 3: Generate visualization
    # ========================================
    if st.button("🎨 Generate Visualization", type="primary"):
        if analysis_type == "Gene expression" and not selected_genes:
            st.error("❌ Please select at least one gene")
            st.stop()

        st.header("Step 3: Visualization")

        with st.spinner("Generating plots..."):
            try:
                fig_files_png = []
                fig_files_pdf = []

                if analysis_type == "Gene expression":
                    # Get pseudotime values
                    pseudotime = adata.obs[selected_pseudotime].values

                    # Filter out NaN pseudotime
                    valid_mask = ~np.isnan(pseudotime)
                    pseudotime = pseudotime[valid_mask]

                    # Get expression data
                    if hasattr(adata.X, 'toarray'):
                        expr_data = adata.X.toarray()[valid_mask, :]
                    else:
                        expr_data = adata.X[valid_mask, :]

                    # Get gene indices
                    gene_indices = [adata.var_names.get_loc(g) for g in selected_genes]
                    gene_expr = expr_data[:, gene_indices]

                    # Sort by pseudotime
                    sort_idx = np.argsort(pseudotime)
                    pseudotime_sorted = pseudotime[sort_idx]
                    gene_expr_sorted = gene_expr[sort_idx, :]

                    # Apply log1p transformation if requested
                    if log1p_pseudotime:
                        pseudotime_sorted = np.log1p(pseudotime_sorted)
                        pseudotime_label = f'{selected_pseudotime} (log1p)'
                    else:
                        pseudotime_label = selected_pseudotime

                if analysis_type == "Gene expression" and viz_type == "Line plot":
                    # Map color selection to hex code
                    color_map = {
                        "Red (#E64B35)": "#E64B35",
                        "Blue (#4DBBD5)": "#4DBBD5",
                        "Green (#00A087)": "#00A087",
                        "Purple (#3C5488)": "#3C5488",
                        "Orange (#F39B7F)": "#F39B7F"
                    }
                    selected_color = color_map[line_color]

                    # Get cluster information if split by cluster
                    if split_by_cluster:
                        cluster_assignments = adata.obs[split_cluster_key].values[valid_mask]
                        cluster_assignments_sorted = cluster_assignments[sort_idx]

                        # Get unique clusters
                        unique_clusters = adata.obs[split_cluster_key].unique()
                        if hasattr(unique_clusters, 'categories'):
                            unique_clusters = unique_clusters.categories.tolist()
                        else:
                            unique_clusters = sorted([str(c) for c in unique_clusters if pd.notna(c)])

                        # Get cluster colors
                        cluster_colors = sns.color_palette(split_cluster_palette, len(unique_clusters))
                        cluster_color_dict = {str(cluster): color for cluster, color in zip(unique_clusters, cluster_colors)}

                    if line_plot_mode == "Individual genes":
                        st.subheader("Gene Expression Trends (Individual)")

                        # Calculate grid layout
                        n_genes = len(selected_genes)
                        nrows = int(np.ceil(n_genes / ncols))

                        fig, axes = plt.subplots(nrows, ncols, figsize=(fig_width, fig_height))
                        if n_genes == 1:
                            axes = np.array([axes])
                        axes = axes.flatten()

                        for idx, gene in enumerate(selected_genes):
                            ax = axes[idx]
                            gene_idx = gene_indices[idx]
                            expr = gene_expr_sorted[:, idx]

                            # Apply Z-score normalization if requested
                            if use_zscore_line:
                                expr = (expr - np.mean(expr)) / (np.std(expr) + 1e-10)

                            if split_by_cluster:
                                # Plot cluster-specific lines
                                for cluster in unique_clusters:
                                    cluster_mask = cluster_assignments_sorted == cluster
                                    if np.sum(cluster_mask) == 0:
                                        continue

                                    cluster_pseudotime = pseudotime_sorted[cluster_mask]
                                    cluster_expr = expr[cluster_mask]

                                    # Plot scatter if requested
                                    if show_scatter:
                                        ax.scatter(cluster_pseudotime, cluster_expr,
                                                 alpha=0.1, s=1, color=cluster_color_dict[str(cluster)])

                                    # Bin and calculate mean and standard error
                                    bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                                    bin_indices_cluster = np.digitize(cluster_pseudotime, bins)

                                    bin_means = []
                                    bin_centers = []
                                    bin_sems = []
                                    for i in range(1, len(bins)):
                                        mask = bin_indices_cluster == i
                                        if np.sum(mask) > 0:
                                            bin_centers.append((bins[i-1] + bins[i]) / 2)
                                            bin_means.append(np.mean(cluster_expr[mask]))
                                            if np.sum(mask) > 1:
                                                bin_sems.append(np.std(cluster_expr[mask]) / np.sqrt(np.sum(mask)))
                                            else:
                                                bin_sems.append(0)

                                    if len(bin_centers) > 3:
                                        # Smooth the binned data with spline
                                        if smoothing_mode == "Strong (smooth)":
                                            mean_smoothing = len(bin_centers) * 0.5
                                            sem_smoothing = len(bin_centers) * 2.0
                                        else:
                                            mean_smoothing = 0
                                            sem_smoothing = 0

                                        mean_spline = UnivariateSpline(bin_centers, bin_means, s=mean_smoothing, k=3)
                                        x_smooth = np.linspace(min(bin_centers), max(bin_centers), 300)
                                        y_smooth = mean_spline(x_smooth)

                                        sem_spline = UnivariateSpline(bin_centers, bin_sems, s=sem_smoothing, k=3)
                                        sem_smooth = sem_spline(x_smooth)

                                        if show_ci:
                                            ax.fill_between(x_smooth, y_smooth - 1.96*sem_smooth,
                                                           y_smooth + 1.96*sem_smooth,
                                                           alpha=0.2, color=cluster_color_dict[str(cluster)])

                                        ax.plot(x_smooth, y_smooth, linewidth=2,
                                               color=cluster_color_dict[str(cluster)],
                                               label=str(cluster), alpha=0.8)
                                    elif len(bin_centers) > 0:
                                        ax.plot(bin_centers, bin_means, linewidth=2,
                                               color=cluster_color_dict[str(cluster)],
                                               label=str(cluster), alpha=0.8)

                                # Add legend for clusters
                                ax.legend(title=split_cluster_key, fontsize=6, frameon=False)

                            else:
                                # Original single-line plot (no cluster split)
                                # Plot scatter if requested
                                if show_scatter:
                                    ax.scatter(pseudotime_sorted, expr, alpha=0.1, s=1, color='gray')

                                # Bin and calculate mean and standard error
                                bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                                bin_indices = np.digitize(pseudotime_sorted, bins)

                                bin_means = []
                                bin_centers = []
                                bin_sems = []  # Standard error of mean
                                for i in range(1, len(bins)):
                                    mask = bin_indices == i
                                    if np.sum(mask) > 0:
                                        bin_centers.append((bins[i-1] + bins[i]) / 2)
                                        bin_means.append(np.mean(expr[mask]))
                                        if np.sum(mask) > 1:
                                            bin_sems.append(np.std(expr[mask]) / np.sqrt(np.sum(mask)))
                                        else:
                                            bin_sems.append(0)

                                if len(bin_centers) > 3:
                                    # Smooth the binned data with spline
                                    # Choose smoothing parameters based on mode
                                    if smoothing_mode == "Strong (smooth)":
                                        # Strong smoothing for publication-quality plots
                                        mean_smoothing = len(bin_centers) * 0.5
                                        sem_smoothing = len(bin_centers) * 2.0
                                    else:  # Interpolation (detailed)
                                        # Perfect interpolation through all bin points
                                        mean_smoothing = 0
                                        sem_smoothing = 0

                                    mean_spline = UnivariateSpline(bin_centers, bin_means, s=mean_smoothing, k=3)
                                    x_smooth = np.linspace(min(bin_centers), max(bin_centers), 300)
                                    y_smooth = mean_spline(x_smooth)

                                    # Smooth the standard errors
                                    sem_spline = UnivariateSpline(bin_centers, bin_sems, s=sem_smoothing, k=3)
                                    sem_smooth = sem_spline(x_smooth)

                                    # Plot confidence interval (95% CI) if requested
                                    if show_ci:
                                        ax.fill_between(x_smooth, y_smooth - 1.96*sem_smooth,
                                                       y_smooth + 1.96*sem_smooth,
                                                       alpha=0.2, color=selected_color)

                                    # Plot smoothed line
                                    ax.plot(x_smooth, y_smooth, linewidth=2, color=selected_color)
                                elif len(bin_centers) > 0:
                                    # Fall back to simple line if too few bins
                                    ax.plot(bin_centers, bin_means, linewidth=2, color=selected_color)

                            ax.set_xlabel(f'{pseudotime_label}')
                            ax.set_ylabel('Z-score' if use_zscore_line else 'Expression')
                            ax.set_title(gene, fontweight='bold')
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            ax.grid(True, alpha=0.3)

                            # Invert x-axis if pseudotime direction is reversed
                            if reverse_pseudotime:
                                ax.invert_xaxis()

                        # Remove empty subplots
                        for idx in range(n_genes, len(axes)):
                            fig.delaxes(axes[idx])

                        plt.tight_layout()

                        # Add cluster annotation bar if requested (after tight_layout)
                        if show_cluster_annotation_line:
                            annotation_cluster_sorted = adata.obs[annotation_cluster_key].values[valid_mask][sort_idx]
                            bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                            # Add annotation bar to all subplots
                            for idx in range(n_genes):
                                add_cluster_annotation_bar(axes[idx], pseudotime_sorted, annotation_cluster_sorted,
                                                          bins, annotation_cluster_color_dict, reverse_pseudotime, height_inches=0.2)

                        # Save
                        fig_path_png = os.path.join(pseudotime_vis_temp_dir, "gene_trends_individual.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "gene_trends_individual.pdf")
                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                        fig_files_pdf.append(fig_path_pdf)

                        st.pyplot(fig)
                        plt.close(fig)

                    else:  # Overlay genes
                        st.subheader("Gene Expression Trends (Overlay)")

                        if split_by_cluster:
                            # Create subplots for each cluster
                            n_clusters = len(unique_clusters)
                            ncols_cluster = min(3, n_clusters)
                            nrows_cluster = int(np.ceil(n_clusters / ncols_cluster))

                            fig, axes = plt.subplots(nrows_cluster, ncols_cluster,
                                                    figsize=(fig_width, fig_height))
                            if n_clusters == 1:
                                axes = np.array([axes])
                            if n_clusters > 1:
                                axes = axes.flatten()

                            # Use selected color palette for genes
                            gene_colors = sns.color_palette(overlay_palette, len(selected_genes))

                            for cluster_idx, cluster in enumerate(unique_clusters):
                                ax = axes[cluster_idx] if n_clusters > 1 else axes[0]
                                cluster_mask = cluster_assignments_sorted == cluster

                                if np.sum(cluster_mask) == 0:
                                    ax.text(0.5, 0.5, f'No cells in {cluster}',
                                           ha='center', va='center', transform=ax.transAxes)
                                    ax.set_xticks([])
                                    ax.set_yticks([])
                                    continue

                                cluster_pseudotime = pseudotime_sorted[cluster_mask]

                                for gene_idx, gene in enumerate(selected_genes):
                                    expr = gene_expr_sorted[:, gene_idx]
                                    cluster_expr = expr[cluster_mask]

                                    # Apply Z-score normalization if requested
                                    if use_zscore_line:
                                        cluster_expr = (cluster_expr - np.mean(cluster_expr)) / (np.std(cluster_expr) + 1e-10)

                                    # Bin and calculate mean and standard error
                                    bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                                    bin_indices_cluster = np.digitize(cluster_pseudotime, bins)

                                    bin_means = []
                                    bin_centers = []
                                    bin_sems = []
                                    for i in range(1, len(bins)):
                                        mask = bin_indices_cluster == i
                                        if np.sum(mask) > 0:
                                            bin_centers.append((bins[i-1] + bins[i]) / 2)
                                            bin_means.append(np.mean(cluster_expr[mask]))
                                            if np.sum(mask) > 1:
                                                bin_sems.append(np.std(cluster_expr[mask]) / np.sqrt(np.sum(mask)))
                                            else:
                                                bin_sems.append(0)

                                    if len(bin_centers) > 3:
                                        # Smooth the binned data with spline
                                        if smoothing_mode == "Strong (smooth)":
                                            mean_smoothing = len(bin_centers) * 0.5
                                            sem_smoothing = len(bin_centers) * 2.0
                                        else:
                                            mean_smoothing = 0
                                            sem_smoothing = 0

                                        mean_spline = UnivariateSpline(bin_centers, bin_means, s=mean_smoothing, k=3)
                                        x_smooth = np.linspace(min(bin_centers), max(bin_centers), 300)
                                        y_smooth = mean_spline(x_smooth)

                                        sem_spline = UnivariateSpline(bin_centers, bin_sems, s=sem_smoothing, k=3)
                                        sem_smooth = sem_spline(x_smooth)

                                        if show_ci:
                                            ax.fill_between(x_smooth, y_smooth - 1.96*sem_smooth,
                                                           y_smooth + 1.96*sem_smooth,
                                                           alpha=0.2, color=gene_colors[gene_idx])

                                        ax.plot(x_smooth, y_smooth, linewidth=2,
                                               label=gene, color=gene_colors[gene_idx], alpha=0.8)
                                    elif len(bin_centers) > 0:
                                        ax.plot(bin_centers, bin_means, linewidth=2,
                                               label=gene, color=gene_colors[gene_idx], alpha=0.8)

                                ax.set_xlabel(f'{pseudotime_label}', fontsize=10)
                                ax.set_ylabel('Z-score' if use_zscore_line else 'Expression', fontsize=10)
                                ax.set_title(f'{cluster}', fontweight='bold', fontsize=12)
                                ax.spines['top'].set_visible(False)
                                ax.spines['right'].set_visible(False)
                                ax.grid(True, alpha=0.3)

                                if reverse_pseudotime:
                                    ax.invert_xaxis()

                                if show_legend:
                                    ax.legend(fontsize=6, frameon=False)

                            # Remove empty subplots
                            for idx in range(n_clusters, len(axes)):
                                fig.delaxes(axes[idx])

                        else:
                            # Original single overlay plot (no cluster split)
                            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                            # Use selected color palette
                            colors = sns.color_palette(overlay_palette, len(selected_genes))

                            for idx, gene in enumerate(selected_genes):
                                expr = gene_expr_sorted[:, idx]

                                # Apply Z-score normalization if requested
                                if use_zscore_line:
                                    expr = (expr - np.mean(expr)) / (np.std(expr) + 1e-10)

                                # Bin and calculate mean and standard error
                                bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                                bin_indices = np.digitize(pseudotime_sorted, bins)

                                bin_means = []
                                bin_centers = []
                                bin_sems = []  # Standard error of mean
                                for i in range(1, len(bins)):
                                    mask = bin_indices == i
                                    if np.sum(mask) > 0:
                                        bin_centers.append((bins[i-1] + bins[i]) / 2)
                                        bin_means.append(np.mean(expr[mask]))
                                        if np.sum(mask) > 1:
                                            bin_sems.append(np.std(expr[mask]) / np.sqrt(np.sum(mask)))
                                        else:
                                            bin_sems.append(0)

                                if len(bin_centers) > 3:
                                    # Smooth the binned data with spline
                                    # Choose smoothing parameters based on mode
                                    if smoothing_mode == "Strong (smooth)":
                                        # Strong smoothing for publication-quality plots
                                        mean_smoothing = len(bin_centers) * 0.5
                                        sem_smoothing = len(bin_centers) * 2.0
                                    else:  # Interpolation (detailed)
                                        # Perfect interpolation through all bin points
                                        mean_smoothing = 0
                                        sem_smoothing = 0

                                    mean_spline = UnivariateSpline(bin_centers, bin_means, s=mean_smoothing, k=3)
                                    x_smooth = np.linspace(min(bin_centers), max(bin_centers), 300)
                                    y_smooth = mean_spline(x_smooth)

                                    # Smooth the standard errors
                                    sem_spline = UnivariateSpline(bin_centers, bin_sems, s=sem_smoothing, k=3)
                                    sem_smooth = sem_spline(x_smooth)

                                    # Plot confidence interval (95% CI) if requested
                                    if show_ci:
                                        ax.fill_between(x_smooth, y_smooth - 1.96*sem_smooth,
                                                       y_smooth + 1.96*sem_smooth,
                                                       alpha=0.2, color=colors[idx])

                                    # Plot smoothed line
                                    ax.plot(x_smooth, y_smooth, linewidth=2,
                                           label=gene, color=colors[idx], alpha=0.8)
                                elif len(bin_centers) > 0:
                                    # Fall back to simple line if too few bins
                                    ax.plot(bin_centers, bin_means, linewidth=2,
                                           label=gene, color=colors[idx], alpha=0.8)

                            ax.set_xlabel(f'{pseudotime_label}', fontsize=12)
                            ax.set_ylabel('Z-score' if use_zscore_line else 'Expression', fontsize=12)
                            ax.set_title(f'Gene Expression Trends ({len(selected_genes)} genes)',
                                        fontweight='bold', fontsize=14)
                            ax.spines['top'].set_visible(False)
                            ax.spines['right'].set_visible(False)
                            ax.grid(True, alpha=0.3)

                            # Invert x-axis if pseudotime direction is reversed
                            if reverse_pseudotime:
                                ax.invert_xaxis()

                            if show_legend:
                                ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left',
                                         frameon=False, fontsize=8)

                        plt.tight_layout()

                        # Add cluster annotation bar if requested (after tight_layout)
                        if show_cluster_annotation_line and not split_by_cluster:
                            annotation_cluster_sorted = adata.obs[annotation_cluster_key].values[valid_mask][sort_idx]
                            bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                            add_cluster_annotation_bar(ax, pseudotime_sorted, annotation_cluster_sorted,
                                                      bins, annotation_cluster_color_dict, reverse_pseudotime, height_inches=0.2)

                        # Save
                        fig_path_png = os.path.join(pseudotime_vis_temp_dir, "gene_trends_overlay.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "gene_trends_overlay.pdf")
                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                        fig_files_pdf.append(fig_path_pdf)

                        st.pyplot(fig)
                        plt.close(fig)

                elif analysis_type == "Gene expression" and viz_type == "Heatmap":
                    st.subheader("Gene Expression Heatmap")

                    # Bin expression data
                    bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                    bin_indices = np.digitize(pseudotime_sorted, bins)

                    heatmap_data = np.zeros((len(selected_genes), n_bins-1))

                    for i in range(1, len(bins)):
                        mask = bin_indices == i
                        if np.sum(mask) > 0:
                            heatmap_data[:, i-1] = np.mean(gene_expr_sorted[mask, :], axis=0)

                    # Normalize each gene (row) if requested
                    if use_zscore:
                        heatmap_data_norm = (heatmap_data - heatmap_data.mean(axis=1, keepdims=True)) / (heatmap_data.std(axis=1, keepdims=True) + 1e-10)
                        cbar_label = 'Z-score'
                        center_value = 0
                    else:
                        heatmap_data_norm = heatmap_data
                        cbar_label = 'Expression'
                        center_value = None

                    # Apply smooth interpolation if requested
                    if use_interpolation:
                        # Upsample using bilinear interpolation (zoom by factor of 3)
                        zoom_factor = 3
                        heatmap_data_interp = zoom(heatmap_data_norm, (1, zoom_factor), order=1)
                        # Generate interpolated column labels
                        interp_bins = np.linspace(bins[0], bins[-1], heatmap_data_interp.shape[1] + 1)
                        interp_columns = [f"{interp_bins[i]:.2f}" for i in range(len(interp_bins)-1)]
                    else:
                        heatmap_data_interp = heatmap_data_norm
                        interp_columns = [f"{bins[i]:.2f}" for i in range(len(bins)-1)]

                    # Create DataFrame
                    heatmap_df = pd.DataFrame(
                        heatmap_data_interp,
                        index=selected_genes,
                        columns=interp_columns
                    )

                    # Calculate bin annotations (dominant cluster for each bin)
                    bin_cluster_annotation = None
                    if show_cluster_annotation and cluster_annotation_key is not None:
                        bin_cluster_labels = []

                        for i in range(1, len(bins)):
                            mask = bin_indices == i
                            if np.sum(mask) > 0:
                                # Get cells in this bin
                                bin_cluster_values = adata.obs[cluster_annotation_key].values[valid_mask][sort_idx][mask]
                                # Find most common cluster (mode)
                                unique, counts = np.unique(bin_cluster_values, return_counts=True)
                                dominant_cluster = str(unique[np.argmax(counts)])
                                bin_cluster_labels.append(dominant_cluster)
                            else:
                                bin_cluster_labels.append("N/A")

                        # Reverse cluster labels if pseudotime direction is reversed
                        if reverse_pseudotime:
                            bin_cluster_labels = bin_cluster_labels[::-1]

                        # If interpolation is enabled, expand cluster labels to match interpolated columns
                        if use_interpolation:
                            # Repeat each cluster label for the interpolated bins
                            zoom_factor = 3
                            bin_cluster_labels_interp = []
                            for label in bin_cluster_labels:
                                bin_cluster_labels_interp.extend([label] * zoom_factor)
                            # Trim to match exact column count
                            bin_cluster_labels_interp = bin_cluster_labels_interp[:len(heatmap_df.columns)]
                        else:
                            bin_cluster_labels_interp = bin_cluster_labels

                        # Create color list for col_colors (use DataFrame to add label)
                        bin_colors = [cluster_color_map.get(str(c), (0.5, 0.5, 0.5)) for c in bin_cluster_labels_interp]
                        bin_cluster_annotation = pd.DataFrame({
                            "Cluster": bin_colors
                        }, index=heatmap_df.columns)

                    # Plot with hierarchical clustering or standard heatmap
                    if use_clustermap:
                        # Use clustermap with dendrograms
                        g = sns.clustermap(
                            heatmap_df,
                            cmap=heatmap_cmap,
                            center=center_value,
                            cbar_kws={'label': cbar_label},
                            xticklabels=False if n_bins > 50 else True,
                            yticklabels=show_gene_names,
                            row_cluster=cluster_rows,
                            col_cluster=cluster_cols,
                            method=linkage_method,
                            figsize=(fig_width, fig_height),
                            dendrogram_ratio=(0.15, 0.15),
                            cbar_pos=(0.02, 0.8, 0.03, 0.15),
                            col_colors=bin_cluster_annotation if bin_cluster_annotation is not None else None
                        )

                        g.ax_heatmap.set_xlabel(f'{pseudotime_label}', fontsize=12)
                        g.ax_heatmap.set_ylabel('Genes', fontsize=12)
                        g.fig.suptitle(f'Gene Expression Heatmap ({len(selected_genes)} genes)',
                                      fontweight='bold', fontsize=14, y=0.98)

                        # Add x-axis labels at start, middle, end (for clustermap)
                        if n_bins <= 50:
                            g.ax_heatmap.set_xticks([0, len(bins)//2, len(bins)-1])
                            g.ax_heatmap.set_xticklabels([f"{bins[0]:.2f}", f"{bins[len(bins)//2]:.2f}", f"{bins[-1]:.2f}"])

                        # Invert x-axis if pseudotime direction is reversed
                        if reverse_pseudotime:
                            g.ax_heatmap.invert_xaxis()

                        fig = g.fig
                    else:
                        # Standard heatmap with optional manual clustering
                        if cluster_genes and len(selected_genes) > 1:
                            from scipy.cluster.hierarchy import linkage, dendrogram
                            linkage_matrix = linkage(heatmap_data_norm, method='ward')
                            dendro = dendrogram(linkage_matrix, no_plot=True)
                            gene_order = dendro['leaves']
                            heatmap_df = heatmap_df.iloc[gene_order, :]

                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        sns.heatmap(
                            heatmap_df,
                            cmap=heatmap_cmap,
                            center=center_value,
                            cbar_kws={'label': cbar_label},
                            xticklabels=False if n_bins > 50 else True,
                            yticklabels=show_gene_names,
                            ax=ax
                        )

                        ax.set_xlabel(f'{pseudotime_label}', fontsize=12)
                        ax.set_ylabel('Genes', fontsize=12)
                        ax.set_title(f'Gene Expression Heatmap ({len(selected_genes)} genes)',
                                    fontweight='bold', fontsize=14)

                        # Add x-axis labels at start, middle, end
                        ax.set_xticks([0, len(bins)//2, len(bins)-1])
                        ax.set_xticklabels([f"{bins[0]:.2f}", f"{bins[len(bins)//2]:.2f}", f"{bins[-1]:.2f}"])

                        # Invert x-axis if pseudotime direction is reversed
                        if reverse_pseudotime:
                            ax.invert_xaxis()

                        plt.tight_layout()

                    # Save
                    fig_path_png = os.path.join(pseudotime_vis_temp_dir, "gene_expression_heatmap.png")
                    fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                    fig_files_png.append(fig_path_png)

                    fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "gene_expression_heatmap.pdf")
                    fig.savefig(fig_path_pdf, bbox_inches='tight')
                    fig_files_pdf.append(fig_path_pdf)

                    st.pyplot(fig)
                    plt.close(fig)

                # ========================================
                # Cluster Density Visualization
                # ========================================
                elif analysis_type == "Cluster density":
                    # Get pseudotime values
                    pseudotime = adata.obs[selected_pseudotime].values

                    # Filter out NaN pseudotime
                    valid_mask = ~np.isnan(pseudotime)
                    pseudotime_valid = pseudotime[valid_mask]

                    # Get cluster assignments for valid cells
                    cluster_assignments = adata.obs[cluster_key].values[valid_mask]

                    # Sort by pseudotime
                    sort_idx = np.argsort(pseudotime_valid)
                    pseudotime_sorted = pseudotime_valid[sort_idx]
                    cluster_sorted = cluster_assignments[sort_idx]

                    # Apply log1p transformation if requested
                    if log1p_pseudotime:
                        pseudotime_sorted = np.log1p(pseudotime_sorted)
                        pseudotime_label = f'{selected_pseudotime} (log1p)'
                    else:
                        pseudotime_label = selected_pseudotime

                    # Create bins
                    bins = np.linspace(pseudotime_sorted.min(), pseudotime_sorted.max(), n_bins)
                    bin_indices = np.digitize(pseudotime_sorted, bins)

                    # Calculate cluster density for each bin
                    cluster_density = np.zeros((len(clusters), n_bins-1))

                    for i in range(1, len(bins)):
                        bin_mask = bin_indices == i
                        if np.sum(bin_mask) > 0:
                            bin_clusters = cluster_sorted[bin_mask]
                            # Count cells in each cluster
                            for j, cluster in enumerate(clusters):
                                count = np.sum(bin_clusters == cluster)
                                cluster_density[j, i-1] = count

                    # Normalize by cluster size if requested
                    if normalize_cluster_size:
                        # Calculate total cells per cluster across all bins
                        cluster_totals = cluster_density.sum(axis=1, keepdims=True)
                        cluster_totals[cluster_totals == 0] = 1  # Avoid division by zero
                        # Normalize each cluster to have the same total (sum to 1)
                        cluster_density = cluster_density / cluster_totals
                        # Scale back to make values more readable (multiply by mean cluster size)
                        mean_cluster_size = np.mean(cluster_totals)
                        cluster_density = cluster_density * mean_cluster_size

                    # Convert to proportion if requested
                    if density_type == "Proportion (%)":
                        # Calculate proportion for each bin
                        bin_totals = cluster_density.sum(axis=0)
                        bin_totals[bin_totals == 0] = 1  # Avoid division by zero
                        cluster_density = (cluster_density / bin_totals) * 100

                    # Generate visualization based on type
                    if viz_type == "Line plot":
                        st.subheader("Cluster Density - Line Plot")

                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        # Use color palette
                        colors = sns.color_palette("husl", len(clusters))

                        # Bin centers for x-axis
                        bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]

                        for j, cluster in enumerate(clusters):
                            # Apply spline smoothing on binned data
                            if len(bin_centers) > 3:
                                # Smooth with spline interpolation (with smoothing factor)
                                smoothing_factor = len(bin_centers) * 0.5
                                spline = UnivariateSpline(bin_centers, cluster_density[j, :], s=smoothing_factor, k=3)
                                x_smooth = np.linspace(min(bin_centers), max(bin_centers), 300)
                                y_smooth = spline(x_smooth)

                                # Plot smoothed line
                                ax.plot(x_smooth, y_smooth, linewidth=2,
                                       label=str(cluster), color=colors[j], alpha=0.8)
                            elif len(bin_centers) > 0:
                                # Fall back to straight line if too few points
                                ax.plot(bin_centers, cluster_density[j, :],
                                       linewidth=2, label=str(cluster),
                                       color=colors[j], alpha=0.8)

                        ax.set_xlabel(f'{pseudotime_label}', fontsize=12)
                        if density_type == "Proportion (%)":
                            ylabel = 'Proportion (%)' + (' (normalized)' if normalize_cluster_size else '')
                            ax.set_ylabel(ylabel, fontsize=12)
                            ax.set_ylim([0, 100])
                        else:
                            ylabel = 'Relative cell count' if normalize_cluster_size else 'Cell count'
                            ax.set_ylabel(ylabel, fontsize=12)

                        ax.set_title(f'Cluster Density along Pseudotime',
                                    fontweight='bold', fontsize=14)
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.grid(True, alpha=0.3)
                        ax.legend(title=cluster_key, bbox_to_anchor=(1.05, 1),
                                 loc='upper left', frameon=False)

                        # Invert x-axis if pseudotime direction is reversed
                        if reverse_pseudotime:
                            ax.invert_xaxis()

                        plt.tight_layout()

                        # Add cluster annotation bar if requested (after tight_layout)
                        if show_cluster_annotation_density:
                            annotation_cluster_sorted_density = adata.obs[annotation_cluster_key_density].values[valid_mask][sort_idx]
                            add_cluster_annotation_bar(ax, pseudotime_sorted, annotation_cluster_sorted_density,
                                                      bins, annotation_cluster_color_dict_density, reverse_pseudotime, height_inches=0.2)

                        # Save
                        fig_path_png = os.path.join(pseudotime_vis_temp_dir, "cluster_density_lineplot.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "cluster_density_lineplot.pdf")
                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                        fig_files_pdf.append(fig_path_pdf)

                        st.pyplot(fig)
                        plt.close(fig)

                    elif viz_type == "Stacked area":
                        st.subheader("Cluster Density - Stacked Area")

                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        # Use color palette
                        colors = sns.color_palette("husl", len(clusters))

                        # Bin centers for x-axis
                        bin_centers = [(bins[i] + bins[i+1]) / 2 for i in range(len(bins)-1)]

                        # For stacked area, always use proportions (convert if needed)
                        if density_type == "Cell count":
                            # Convert to proportions for stacking
                            cluster_density_prop = cluster_density.copy()
                            bin_totals = cluster_density_prop.sum(axis=0)
                            bin_totals[bin_totals == 0] = 1
                            cluster_density_prop = (cluster_density_prop / bin_totals) * 100
                        else:
                            cluster_density_prop = cluster_density

                        # Create stacked area plot
                        ax.stackplot(bin_centers, cluster_density_prop,
                                    labels=[str(c) for c in clusters],
                                    colors=colors, alpha=0.8)

                        ax.set_xlabel(f'{pseudotime_label}', fontsize=12)
                        ax.set_ylabel('Proportion (%)', fontsize=12)
                        ax.set_ylim([0, 100])
                        ax.set_title(f'Cluster Composition along Pseudotime',
                                    fontweight='bold', fontsize=14)
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                        ax.grid(True, alpha=0.3, axis='y')
                        ax.legend(title=cluster_key, bbox_to_anchor=(1.05, 1),
                                 loc='upper left', frameon=False)

                        # Invert x-axis if pseudotime direction is reversed
                        if reverse_pseudotime:
                            ax.invert_xaxis()

                        plt.tight_layout()

                        # Add cluster annotation bar if requested (after tight_layout)
                        if show_cluster_annotation_density:
                            annotation_cluster_sorted_density = adata.obs[annotation_cluster_key_density].values[valid_mask][sort_idx]
                            add_cluster_annotation_bar(ax, pseudotime_sorted, annotation_cluster_sorted_density,
                                                      bins, annotation_cluster_color_dict_density, reverse_pseudotime, height_inches=0.2)

                        # Save
                        fig_path_png = os.path.join(pseudotime_vis_temp_dir, "cluster_density_stacked.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "cluster_density_stacked.pdf")
                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                        fig_files_pdf.append(fig_path_pdf)

                        st.pyplot(fig)
                        plt.close(fig)

                    elif viz_type == "Heatmap":
                        st.subheader("Cluster Density - Heatmap")

                        # Create DataFrame for heatmap
                        bin_labels = [f"{bins[i]:.2f}" for i in range(len(bins)-1)]
                        heatmap_df = pd.DataFrame(
                            cluster_density,
                            index=[str(c) for c in clusters],
                            columns=bin_labels
                        )

                        # Plot
                        fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                        # Determine colorbar label based on density type and normalization
                        if density_type == "Proportion (%)":
                            cbar_label = 'Proportion (%)' + (' (normalized)' if normalize_cluster_size else '')
                        else:
                            cbar_label = 'Relative cell count' if normalize_cluster_size else 'Cell count'

                        sns.heatmap(
                            heatmap_df,
                            cmap=heatmap_cmap,
                            cbar_kws={'label': cbar_label},
                            xticklabels=False if n_bins > 50 else True,
                            yticklabels=True,
                            ax=ax
                        )

                        ax.set_xlabel(f'{pseudotime_label}', fontsize=12)
                        ax.set_ylabel(cluster_key, fontsize=12)
                        ax.set_title(f'Cluster Density Heatmap',
                                    fontweight='bold', fontsize=14)

                        # Add x-axis labels at start, middle, end
                        ax.set_xticks([0, len(bins)//2, len(bins)-1])
                        ax.set_xticklabels([f"{bins[0]:.2f}",
                                          f"{bins[len(bins)//2]:.2f}",
                                          f"{bins[-1]:.2f}"])

                        # Invert x-axis if pseudotime direction is reversed
                        if reverse_pseudotime:
                            ax.invert_xaxis()

                        plt.tight_layout()

                        # Save
                        fig_path_png = os.path.join(pseudotime_vis_temp_dir, "cluster_density_heatmap.png")
                        fig.savefig(fig_path_png, bbox_inches='tight', dpi=300)
                        fig_files_png.append(fig_path_png)

                        fig_path_pdf = os.path.join(pseudotime_vis_temp_dir, "cluster_density_heatmap.pdf")
                        fig.savefig(fig_path_pdf, bbox_inches='tight')
                        fig_files_pdf.append(fig_path_pdf)

                        st.pyplot(fig)
                        plt.close(fig)

                # Download section
                if fig_files_png:
                    st.success(f"✅ Generated visualization")

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
    st.info("👆 h5adファイルをアップロードして開始してください")
