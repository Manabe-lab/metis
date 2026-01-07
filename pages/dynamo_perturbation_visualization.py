"""
Dynamo Perturbation Visualization
Visualize perturbation analysis results including Jacobian matrices and regulatory networks
"""

import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import time
from helper_func import clear_old_directories, clear_old_files

# Import network analysis libraries
try:
    import networkx as nx
    NETWORKX_AVAILABLE = True
except ImportError:
    NETWORKX_AVAILABLE = False

# Import dynamo
try:
    import dynamo as dyn
    DYNAMO_AVAILABLE = True
except ImportError:
    DYNAMO_AVAILABLE = False

st.set_page_config(page_title="Dynamo Perturbation Visualization", page_icon="🔬", layout="wide")

st.title("🔬 Dynamo Perturbation Visualization")

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
**Dynamo Prediction and Network**アプリで計算したperturbation解析結果を可視化します。

### 可視化タイプ
- **Jacobian heatmap**: 遺伝子制御ネットワークのヒートマップ
- **Regulatory network graph**: 遺伝子制御ネットワークのネットワーク図
- **Top regulatory interactions**: 上位の制御関係のランキング
- **Jacobian statistics**: Jacobian行列の統計情報

### 必要なデータ
- **Dynamo Prediction and Network**アプリで生成されたh5adファイル
- Jacobian行列が計算済み（`adata.uns['jacobian']`）

### 参考
- [Qiu et al. (2022) "Mapping transcriptomic vector fields of single cells" Cell](https://www.cell.com/cell/fulltext/S0092-8674(21)01577-4)
- [Dynamo Documentation - Perturbation](https://dynamo-release.readthedocs.io/en/latest/API.html#module-dynamo.prediction)
""")

# Initialize session state
if "dynamo_pert_vis_temp_dir" not in st.session_state:
    dynamo_pert_vis_temp_dir = os.path.join("temp", f"dynamo_pert_vis_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(dynamo_pert_vis_temp_dir, exist_ok=True)
    st.session_state.dynamo_pert_vis_temp_dir = dynamo_pert_vis_temp_dir
else:
    dynamo_pert_vis_temp_dir = st.session_state.dynamo_pert_vis_temp_dir

# ========================================
# Upload file
# ========================================
st.header("Step 1: Upload perturbation result")

uploaded_h5ad = st.file_uploader(
    "Upload h5ad file (Perturbation result)",
    type=['h5ad'],
    key="dynamo_pert_vis_h5ad_upload",
    help="Dynamo Prediction and Networkアプリ（dynamo_perturbation.py）で生成されたh5adファイル"
)

if uploaded_h5ad is not None:
    st.success("✓ File uploaded")

    # Load data
    if ("dynamo_pert_vis_adata" not in st.session_state or
        st.session_state.get("dynamo_pert_vis_uploaded_file") != uploaded_h5ad.name):

        with st.spinner("Loading data..."):
            # Save uploaded file temporarily
            temp_h5ad_path = os.path.join(dynamo_pert_vis_temp_dir, "perturbation_result.h5ad")
            with open(temp_h5ad_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            # Read h5ad
            adata = sc.read_h5ad(temp_h5ad_path)

            st.session_state.dynamo_pert_vis_adata = adata
            st.session_state.dynamo_pert_vis_uploaded_file = uploaded_h5ad.name

            st.info(f"✓ Loaded: {adata.n_obs} cells, {adata.n_vars} genes")

    adata = st.session_state.dynamo_pert_vis_adata

    # Check for required data
    has_jacobian = 'jacobian' in adata.uns

    if not has_jacobian:
        st.error("""
        ❌ **Jacobian matrix not found**

        このファイルにはJacobian行列が含まれていません。

        **Dynamo Prediction and Network**アプリ（pages/dynamo_perturbation.py）で先にJacobian解析を実行してください：
        1. Dynamo解析済みh5adファイルをアップロード
        2. "Jacobian Analysis (Regulatory Network)" を選択
        3. Regulator genesとEffector genesを設定
        4. "Compute Jacobian" をクリック
        5. 結果のh5adファイルをダウンロード
        6. このアプリでダウンロードしたファイルをアップロード
        """)
        st.stop()

    st.success("✓ Jacobian matrix found")

    # Extract Jacobian information
    jac_data = adata.uns['jacobian']
    jac_matrix = jac_data['jacobian_gene']
    regulators = jac_data['regulators']
    effectors = jac_data['effectors']

    st.info(f"✓ Jacobian matrix: {len(regulators)} regulators × {len(effectors)} effectors")

    # ========================================
    # Visualization options
    # ========================================
    st.header("Step 2: Select visualization type")

    viz_type = st.selectbox(
        "Visualization type",
        [
            "Jacobian heatmap",
            "Regulatory network graph",
            "Top regulatory interactions",
            "Jacobian statistics"
        ],
        help="可視化のタイプを選択"
    )

    # ========================================
    # Visualization parameters
    # ========================================
    st.subheader("Visualization parameters")

    with st.expander("📚 Parameter Guide", expanded=False):
        st.markdown("""
        ### Jacobian heatmap
        - 遺伝子制御ネットワークをヒートマップで表示
        - 赤: 正の制御（活性化）
        - 青: 負の制御（抑制）
        - 表示する遺伝子数を調整可能

        ### Regulatory network graph
        - 遺伝子制御ネットワークをグラフで表示
        - ノード: 遺伝子
        - エッジ: 制御関係（太さは制御強度）
        - 上位の制御関係のみ表示（閾値で調整）

        ### Top regulatory interactions
        - 制御強度の上位ランキング
        - 正の制御と負の制御を分けて表示
        - バープロットで可視化

        ### Jacobian statistics
        - Jacobian行列の統計情報
        - 分布、相関、主成分分析など
        """)

    # Sidebar for visualization options
    with st.sidebar:
        st.markdown("### Visualization Options")

        colormap_continuous = st.selectbox(
            "Colormap (diverging):",
            ["RdBu_r", "RdYlBu_r", "seismic", "coolwarm", "bwr", "PiYG", "PRGn"],
            index=0,
            help="Jacobianヒートマップ用のdivergingカラーマップ"
        )

    col1, col2, col3 = st.columns(3)

    with col1:
        if viz_type in ["Jacobian heatmap", "Regulatory network graph", "Top regulatory interactions"]:
            n_top_genes = st.slider(
                "Number of top genes to show",
                min_value=5,
                max_value=min(50, len(regulators), len(effectors)),
                value=min(20, len(regulators), len(effectors)),
                help="表示する遺伝子数（上位N個）"
            )

        if viz_type == "Regulatory network graph":
            interaction_threshold = st.slider(
                "Interaction threshold (absolute Jacobian value)",
                min_value=0.0,
                max_value=float(np.abs(jac_matrix).max()),
                value=float(np.percentile(np.abs(jac_matrix), 95)),
                help="表示する制御関係の閾値（絶対値）"
            )

    with col2:
        fig_width = st.slider("Figure width", 2, 20, 10)
        fig_height = st.slider("Figure height", 2, 20, 8)

    with col3:
        if viz_type == "Regulatory network graph":
            node_size = st.slider("Node size", 100, 2000, 500)
            edge_width_scale = st.slider("Edge width scale", 0.1, 5.0, 1.0, 0.1)

    # ========================================
    # Generate visualization
    # ========================================
    if st.button("🎨 Generate Visualization", type="primary"):
        st.header("Step 3: Visualization")

        try:
            # Clear all previous matplotlib figures
            plt.close('all')

            if viz_type == "Jacobian heatmap":
                st.subheader("Jacobian Matrix Heatmap")

                with st.spinner("Generating heatmap..."):
                    # Limit to top N genes
                    n_reg_show = min(n_top_genes, len(regulators))
                    n_eff_show = min(n_top_genes, len(effectors))

                    # Create figure
                    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                    # Plot heatmap
                    sns.heatmap(
                        jac_matrix[:n_reg_show, :n_eff_show],
                        cmap=colormap_continuous,
                        center=0,
                        xticklabels=effectors[:n_eff_show],
                        yticklabels=regulators[:n_reg_show],
                        cbar_kws={'label': 'Jacobian value'},
                        ax=ax
                    )

                    ax.set_title("Gene Regulatory Network (Jacobian Matrix)", fontsize=16)
                    ax.set_xlabel("Effector genes", fontsize=12)
                    ax.set_ylabel("Regulator genes", fontsize=12)
                    plt.xticks(rotation=45, ha='right')
                    plt.yticks(rotation=0)
                    plt.tight_layout()

                    st.pyplot(fig)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_pert_vis_temp_dir, "jacobian_heatmap.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="jacobian_heatmap.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        fig_path_pdf = os.path.join(dynamo_pert_vis_temp_dir, "jacobian_heatmap.pdf")
                        fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                        with open(fig_path_pdf, "rb") as f:
                            st.download_button(
                                "⬇️ Download PDF",
                                f,
                                file_name="jacobian_heatmap.pdf",
                                mime="application/pdf"
                            )

            elif viz_type == "Regulatory network graph":
                st.subheader("Regulatory Network Graph")

                if not NETWORKX_AVAILABLE:
                    st.error("❌ NetworkX is not installed. Please install it using: `pip install networkx`")
                else:
                    with st.spinner("Generating network graph..."):
                        # Create network graph
                        G = nx.DiGraph()

                        # Add nodes
                        all_genes = list(set(regulators) | set(effectors))
                        G.add_nodes_from(all_genes)

                        # Add edges based on threshold
                        edges_to_add = []
                        for i, reg in enumerate(regulators):
                            for j, eff in enumerate(effectors):
                                jac_value = jac_matrix[i, j]
                                if abs(jac_value) >= interaction_threshold:
                                    edges_to_add.append((reg, eff, jac_value))

                        st.info(f"Adding {len(edges_to_add)} edges (threshold: {interaction_threshold:.3f})")

                        for reg, eff, weight in edges_to_add:
                            G.add_edge(reg, eff, weight=weight)

                        if len(G.edges()) == 0:
                            st.warning("⚠️ No edges to display. Try lowering the interaction threshold.")
                        else:
                            # Create figure
                            fig, ax = plt.subplots(figsize=(fig_width, fig_height))

                            # Layout
                            pos = nx.spring_layout(G, k=1, iterations=50, seed=42)

                            # Draw nodes
                            node_colors = ['lightblue' if node in regulators else 'lightgreen' for node in G.nodes()]
                            nx.draw_networkx_nodes(
                                G, pos,
                                node_color=node_colors,
                                node_size=node_size,
                                alpha=0.9,
                                ax=ax
                            )

                            # Draw edges
                            edge_weights = [abs(G[u][v]['weight']) for u, v in G.edges()]
                            edge_widths = [w * edge_width_scale for w in edge_weights]
                            edge_colors = ['red' if G[u][v]['weight'] > 0 else 'blue' for u, v in G.edges()]

                            nx.draw_networkx_edges(
                                G, pos,
                                width=edge_widths,
                                edge_color=edge_colors,
                                alpha=0.5,
                                arrows=True,
                                arrowsize=20,
                                ax=ax
                            )

                            # Draw labels
                            nx.draw_networkx_labels(
                                G, pos,
                                font_size=8,
                                font_color='black',
                                ax=ax
                            )

                            ax.set_title("Regulatory Network Graph", fontsize=16)
                            ax.axis('off')

                            # Add legend
                            from matplotlib.patches import Patch
                            legend_elements = [
                                Patch(facecolor='lightblue', label='Regulators'),
                                Patch(facecolor='lightgreen', label='Effectors'),
                                Patch(facecolor='red', alpha=0.5, label='Activation (+)'),
                                Patch(facecolor='blue', alpha=0.5, label='Repression (-)')
                            ]
                            ax.legend(handles=legend_elements, loc='upper right')

                            plt.tight_layout()
                            st.pyplot(fig)

                            # Network statistics
                            st.markdown("### Network Statistics")
                            col_stat1, col_stat2, col_stat3 = st.columns(3)
                            with col_stat1:
                                st.metric("Nodes", len(G.nodes()))
                            with col_stat2:
                                st.metric("Edges", len(G.edges()))
                            with col_stat3:
                                avg_degree = sum(dict(G.degree()).values()) / len(G.nodes()) if len(G.nodes()) > 0 else 0
                                st.metric("Avg Degree", f"{avg_degree:.2f}")

                            # Save as PNG and PDF
                            col_dl1, col_dl2 = st.columns(2)

                            with col_dl1:
                                fig_path_png = os.path.join(dynamo_pert_vis_temp_dir, "regulatory_network.png")
                                fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                                with open(fig_path_png, "rb") as f:
                                    st.download_button(
                                        "⬇️ Download PNG",
                                        f,
                                        file_name="regulatory_network.png",
                                        mime="image/png"
                                    )

                            with col_dl2:
                                fig_path_pdf = os.path.join(dynamo_pert_vis_temp_dir, "regulatory_network.pdf")
                                fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                                with open(fig_path_pdf, "rb") as f:
                                    st.download_button(
                                        "⬇️ Download PDF",
                                        f,
                                        file_name="regulatory_network.pdf",
                                        mime="application/pdf"
                                    )

            elif viz_type == "Top regulatory interactions":
                st.subheader("Top Regulatory Interactions")

                with st.spinner("Analyzing top interactions..."):
                    # Flatten Jacobian and get top values
                    jac_flat = jac_matrix.flatten()
                    jac_abs = np.abs(jac_flat)

                    # Get top positive and negative interactions
                    n_top = min(n_top_genes, len(jac_flat))

                    # Top positive
                    top_pos_indices = np.argsort(jac_flat)[::-1][:n_top]
                    top_pos_interactions = []
                    for idx in top_pos_indices:
                        reg_idx = idx // len(effectors)
                        eff_idx = idx % len(effectors)
                        if jac_flat[idx] > 0:
                            top_pos_interactions.append({
                                'Regulator': regulators[reg_idx],
                                'Effector': effectors[eff_idx],
                                'Jacobian': jac_flat[idx]
                            })

                    # Top negative
                    top_neg_indices = np.argsort(jac_flat)[:n_top]
                    top_neg_interactions = []
                    for idx in top_neg_indices:
                        reg_idx = idx // len(effectors)
                        eff_idx = idx % len(effectors)
                        if jac_flat[idx] < 0:
                            top_neg_interactions.append({
                                'Regulator': regulators[reg_idx],
                                'Effector': effectors[eff_idx],
                                'Jacobian': jac_flat[idx]
                            })

                    # Create figure with two subplots
                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(fig_width, fig_height))

                    # Plot top positive
                    if top_pos_interactions:
                        df_pos = pd.DataFrame(top_pos_interactions)
                        df_pos['Interaction'] = df_pos['Regulator'] + ' → ' + df_pos['Effector']
                        df_pos = df_pos.sort_values('Jacobian')

                        ax1.barh(df_pos['Interaction'], df_pos['Jacobian'], color='red', alpha=0.7)
                        ax1.set_xlabel('Jacobian value', fontsize=12)
                        ax1.set_title('Top Positive Regulations (Activation)', fontsize=14)
                        ax1.grid(axis='x', alpha=0.3)

                    # Plot top negative
                    if top_neg_interactions:
                        df_neg = pd.DataFrame(top_neg_interactions)
                        df_neg['Interaction'] = df_neg['Regulator'] + ' → ' + df_neg['Effector']
                        df_neg = df_neg.sort_values('Jacobian', ascending=False)

                        ax2.barh(df_neg['Interaction'], df_neg['Jacobian'], color='blue', alpha=0.7)
                        ax2.set_xlabel('Jacobian value', fontsize=12)
                        ax2.set_title('Top Negative Regulations (Repression)', fontsize=14)
                        ax2.grid(axis='x', alpha=0.3)

                    plt.tight_layout()
                    st.pyplot(fig)

                    # Show data tables
                    col_table1, col_table2 = st.columns(2)
                    with col_table1:
                        st.markdown("#### Top Activations")
                        if top_pos_interactions:
                            st.dataframe(pd.DataFrame(top_pos_interactions), use_container_width=True)
                        else:
                            st.info("No positive interactions found")

                    with col_table2:
                        st.markdown("#### Top Repressions")
                        if top_neg_interactions:
                            st.dataframe(pd.DataFrame(top_neg_interactions), use_container_width=True)
                        else:
                            st.info("No negative interactions found")

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_pert_vis_temp_dir, "top_interactions.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download PNG",
                                f,
                                file_name="top_interactions.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        fig_path_pdf = os.path.join(dynamo_pert_vis_temp_dir, "top_interactions.pdf")
                        fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                        with open(fig_path_pdf, "rb") as f:
                            st.download_button(
                                "⬇️ Download PDF",
                                f,
                                file_name="top_interactions.pdf",
                                mime="application/pdf"
                            )

            elif viz_type == "Jacobian statistics":
                st.subheader("Jacobian Matrix Statistics")

                with st.spinner("Computing statistics..."):
                    # Basic statistics
                    st.markdown("### Basic Statistics")
                    col_stat1, col_stat2, col_stat3, col_stat4 = st.columns(4)

                    with col_stat1:
                        st.metric("Mean", f"{jac_matrix.mean():.4f}")
                    with col_stat2:
                        st.metric("Std Dev", f"{jac_matrix.std():.4f}")
                    with col_stat3:
                        st.metric("Min", f"{jac_matrix.min():.4f}")
                    with col_stat4:
                        st.metric("Max", f"{jac_matrix.max():.4f}")

                    # Distribution plot
                    fig, axes = plt.subplots(2, 2, figsize=(fig_width, fig_height))

                    # Histogram
                    axes[0, 0].hist(jac_matrix.flatten(), bins=50, color='steelblue', alpha=0.7, edgecolor='black')
                    axes[0, 0].set_xlabel('Jacobian value')
                    axes[0, 0].set_ylabel('Frequency')
                    axes[0, 0].set_title('Distribution of Jacobian Values')
                    axes[0, 0].axvline(0, color='red', linestyle='--', linewidth=2)
                    axes[0, 0].grid(alpha=0.3)

                    # Absolute value distribution
                    axes[0, 1].hist(np.abs(jac_matrix.flatten()), bins=50, color='coral', alpha=0.7, edgecolor='black')
                    axes[0, 1].set_xlabel('Absolute Jacobian value')
                    axes[0, 1].set_ylabel('Frequency')
                    axes[0, 1].set_title('Distribution of Absolute Values')
                    axes[0, 1].grid(alpha=0.3)

                    # Regulator strength (row-wise sum)
                    regulator_strength = np.abs(jac_matrix).sum(axis=1)
                    top_reg_idx = np.argsort(regulator_strength)[::-1][:20]
                    axes[1, 0].barh(
                        [regulators[i] for i in top_reg_idx],
                        regulator_strength[top_reg_idx],
                        color='lightblue',
                        edgecolor='black'
                    )
                    axes[1, 0].set_xlabel('Total regulatory strength')
                    axes[1, 0].set_title('Top 20 Regulators by Strength')
                    axes[1, 0].grid(axis='x', alpha=0.3)

                    # Effector sensitivity (column-wise sum)
                    effector_sensitivity = np.abs(jac_matrix).sum(axis=0)
                    top_eff_idx = np.argsort(effector_sensitivity)[::-1][:20]
                    axes[1, 1].barh(
                        [effectors[i] for i in top_eff_idx],
                        effector_sensitivity[top_eff_idx],
                        color='lightgreen',
                        edgecolor='black'
                    )
                    axes[1, 1].set_xlabel('Total sensitivity')
                    axes[1, 1].set_title('Top 20 Effectors by Sensitivity')
                    axes[1, 1].grid(axis='x', alpha=0.3)

                    plt.tight_layout()
                    st.pyplot(fig)

                    # Correlation analysis
                    st.markdown("### Regulator-Regulator Correlation")
                    st.info("Showing correlation between regulatory patterns of different regulators")

                    # Compute correlation between regulators
                    reg_corr = np.corrcoef(jac_matrix)

                    fig_corr, ax_corr = plt.subplots(figsize=(fig_width, fig_width))
                    n_show = min(30, len(regulators))
                    sns.heatmap(
                        reg_corr[:n_show, :n_show],
                        cmap='coolwarm',
                        center=0,
                        xticklabels=regulators[:n_show],
                        yticklabels=regulators[:n_show],
                        ax=ax_corr
                    )
                    ax_corr.set_title(f'Regulator-Regulator Correlation (top {n_show})')
                    plt.xticks(rotation=45, ha='right')
                    plt.yticks(rotation=0)
                    plt.tight_layout()
                    st.pyplot(fig_corr)

                    # Save as PNG and PDF
                    col_dl1, col_dl2 = st.columns(2)

                    with col_dl1:
                        fig_path_png = os.path.join(dynamo_pert_vis_temp_dir, "jacobian_statistics.png")
                        fig.savefig(fig_path_png, dpi=300, bbox_inches='tight')
                        with open(fig_path_png, "rb") as f:
                            st.download_button(
                                "⬇️ Download Statistics PNG",
                                f,
                                file_name="jacobian_statistics.png",
                                mime="image/png"
                            )

                    with col_dl2:
                        fig_path_pdf = os.path.join(dynamo_pert_vis_temp_dir, "jacobian_statistics.pdf")
                        fig.savefig(fig_path_pdf, format='pdf', bbox_inches='tight')
                        with open(fig_path_pdf, "rb") as f:
                            st.download_button(
                                "⬇️ Download Statistics PDF",
                                f,
                                file_name="jacobian_statistics.pdf",
                                mime="application/pdf"
                            )

            st.success("✅ Visualization generated successfully!")

        except Exception as e:
            st.error(f"❌ Error during visualization: {str(e)}")
            st.exception(e)

else:
    st.info("👆 **Dynamo Prediction and Network**アプリ（pages/dynamo_perturbation.py）で生成されたh5adファイルをアップロードして開始してください")
