"""
scVelo RNA Velocity Analysis
Perform RNA velocity analysis using scVelo or DeepVelo
"""

import streamlit as st
import scanpy as sc
import anndata as ad
import scvelo as scv
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

st.set_page_config(page_title="scVelo Analysis", page_icon="🚀", layout="wide")

st.title("🚀 RNA Velocity Analysis")
st.markdown("""
scVeloまたはDeepVeloを用いてRNA velocity解析を実行します。

### ワークフロー
1. **ファイル読み込み**: h5ad + loom ファイル
2. **マージ**: scvelo.utils.merge でデータ統合
3. **前処理**: フィルタリング、正規化、moments計算
4. **Velocity計算**: scVelo (Stochastic/Dynamical/Deterministic) または DeepVelo
5. **結果保存**: 解析結果をh5adファイルでダウンロード

### 参考
- [Bergen et al. (2020) "Generalizing RNA velocity to transient cell states through dynamical modeling" Nature Biotechnology](https://www.nature.com/articles/s41587-020-0591-3)
- [Gao et al. (2024) "DeepVelo: Deep learning enables accurate estimation of single-cell RNA velocity" Nature Communications](https://www.nature.com/articles/s41467-024-51278-6)
""")

# Initialize session state
if "velocity_analysis_temp_dir" not in st.session_state:
    velocity_analysis_temp_dir = os.path.join("temp", f"scvelo_{round(time.time())}")
    os.makedirs("temp", exist_ok=True)
    clear_old_directories("temp")
    clear_old_files("temp")
    os.makedirs(velocity_analysis_temp_dir, exist_ok=True)
    st.session_state.velocity_analysis_temp_dir = velocity_analysis_temp_dir
else:
    velocity_analysis_temp_dir = st.session_state.velocity_analysis_temp_dir

if "analysis_complete" not in st.session_state:
    st.session_state.analysis_complete = False

# ========================================
# Step 1: Upload files
# ========================================
st.header("Step 1: Upload files")

col1, col2 = st.columns(2)

with col1:
    uploaded_h5ad = st.file_uploader(
        "Upload h5ad file",
        type=['h5ad'],
        key="scvelo_h5ad_upload",
        help="Seurat/Scanpy解析済みのh5adファイル"
    )

with col2:
    uploaded_loom = st.file_uploader(
        "Upload loom file",
        type=['loom'],
        key="scvelo_loom_upload",
        help="velocytoまたはData filtering appで生成されたloomファイル"
    )

if uploaded_h5ad is not None and uploaded_loom is not None:
    st.success("✓ Both files uploaded")

    # Pre-load h5ad metadata to get categorical columns for PAGA
    if "velocity_h5ad_categorical_cols" not in st.session_state or \
       st.session_state.get("velocity_h5ad_name") != uploaded_h5ad.name:

        with st.spinner("Reading h5ad metadata..."):
            temp_preview_path = os.path.join(velocity_analysis_temp_dir, "preview.h5ad")
            with open(temp_preview_path, "wb") as f:
                f.write(uploaded_h5ad.read())

            adata_preview = sc.read_h5ad(temp_preview_path)

            # Get categorical columns
            categorical_cols = [col for col in adata_preview.obs.columns
                               if adata_preview.obs[col].dtype.name == 'category']

            st.session_state.velocity_h5ad_categorical_cols = categorical_cols
            st.session_state.velocity_h5ad_name = uploaded_h5ad.name

            del adata_preview  # Free memory

    categorical_cols = st.session_state.velocity_h5ad_categorical_cols

    # ========================================
    # Step 2: Configure parameters
    # ========================================
    st.header("Step 2: Configure analysis parameters")

    with st.expander("📚 Parameter Guide", expanded=False):
        st.markdown("""
        ### filter_and_normalize
        - **min_shared_counts**: 最小共有カウント数（デフォルト: 20）
        - **n_top_genes**: 使用する高変動遺伝子数（デフォルト: 2000）

        ### moments
        - **n_pcs**: 使用する主成分数（デフォルト: 30）
        - **n_neighbors**: 近傍細胞数（デフォルト: 30）

        ### velocity
        - **mode**: 計算モード
          - **dynamical** (推奨): 最も正確、転写動態の完全なモデリング。原著論文で推奨
          - **stochastic**: 最速、簡易的な解析
          - **deterministic**: 中間的

        - **n_jobs**: 並列計算に使用するCPUコア数（dynamicalモードのみ有効）
          - 32コアシステムでは24コア程度が推奨（メモリ使用量とのバランス）

        ### PAGA (Partition-based graph abstraction)
        - クラスター間の接続性を計算し、分化軌跡を抽象化
        - velocity graphと組み合わせることで、directed PAGAが計算されます
        - クラスター情報（adata.obsのカテゴリカル変数）が必要

        Bergen et al. (2020) Nature Biotechnology では dynamical モードが推奨されています。
        """)

    with st.form("analysis_params_form"):
        st.subheader("Preprocessing parameters")

        col1, col2 = st.columns(2)
        with col1:
            min_shared_counts = st.number_input(
                "min_shared_counts",
                min_value=1,
                max_value=100,
                value=20,
                help="遺伝子と細胞の両方でフィルタリングに使用する最小カウント数"
            )

            n_top_genes = st.number_input(
                "n_top_genes",
                min_value=500,
                max_value=5000,
                value=2000,
                step=100,
                help="高変動遺伝子の数"
            )

            # Retain specific genes (like scVelo's retain_genes parameter)
            st.markdown("##### Retain genes (comma, space, tab, CR separated):")
            retain_genes_str = st.text_area(
                "Retain genes",
                value="",
                height=100,
                label_visibility='collapsed',
                help="HVG選択に関係なく必ず保持する遺伝子名（複数の区切り文字対応: , space tab CR）"
            )

        with col2:
            n_pcs = st.number_input(
                "n_pcs",
                min_value=10,
                max_value=100,
                value=30,
                help="主成分分析で使用する成分数"
            )

            n_neighbors = st.number_input(
                "n_neighbors",
                min_value=5,
                max_value=100,
                value=30,
                help="近傍グラフ構築に使用する近傍細胞数"
            )

        st.subheader("Velocity calculation")

        col1, col2 = st.columns([2, 1])

        with col1:
            # Build velocity mode options (scVelo modes only)
            velocity_options = ["dynamical", "stochastic", "deterministic"]

            velocity_mode = st.selectbox(
                "scVelo mode",
                velocity_options,
                index=0,
                help="scVelo計算モード。dynamicalが推奨（Bergen et al. 2020）"
            )

            if velocity_mode == "dynamical":
                st.info("✨ **Dynamical mode** (Default / Recommended)")

            st.caption("💡 DeepVelo is available in a separate app (DeepVelo analysis)")

        with col2:
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
                help=f"並列計算のコア数（dynamicalモードのみ有効）。利用可能: {n_cpus}コア"
            )

        st.subheader("PAGA (optional)")

        compute_paga = st.checkbox(
            "Compute PAGA",
            value=True,
            help="Partition-based graph abstractionを計算（クラスター情報が必要）"
        )

        if compute_paga:
            st.info("✨ **Compute PAGA** (Default / Recommended)")

        paga_cluster_key = None
        if compute_paga:
            if not categorical_cols:
                st.warning("⚠️ No categorical columns found in h5ad metadata. PAGA requires cluster information.")
                st.info("Your h5ad file should contain cluster information (e.g., 'louvain', 'leiden', 'clusters')")
            else:
                # Find default cluster column (prioritize common names)
                default_cluster = None
                priority_names = ['leiden', 'louvain', 'clusters', 'celltype', 'cell_type', 'seurat_clusters']

                for name in priority_names:
                    if name in categorical_cols:
                        default_cluster = name
                        break

                if default_cluster is None:
                    default_cluster = categorical_cols[0]

                default_index = categorical_cols.index(default_cluster)

                paga_cluster_key = st.selectbox(
                    "Cluster column for PAGA",
                    categorical_cols,
                    index=default_index,
                    help="PAGAの計算に使用するクラスター情報を選択"
                )

                st.info(f"💡 Using '{paga_cluster_key}' for PAGA computation")

        st.markdown("---")

        submit_analysis = st.form_submit_button("🚀 Run Analysis", type="primary")

    # ========================================
    # Step 3: Run analysis
    # ========================================
    if submit_analysis:
        st.header("Step 3: Running analysis")

        with st.spinner("Loading and processing files..."):
            progress_bar = st.progress(0)
            status_text = st.empty()

            try:
                # Load h5ad (use preview.h5ad if it exists, otherwise save new one)
                status_text.text("Loading h5ad file...")
                progress_bar.progress(10)

                temp_h5ad_path = os.path.join(velocity_analysis_temp_dir, "preview.h5ad")
                if not os.path.exists(temp_h5ad_path):
                    with open(temp_h5ad_path, "wb") as f:
                        f.write(uploaded_h5ad.read())

                adata = sc.read_h5ad(temp_h5ad_path)
                st.info(f"✓ Loaded h5ad: {adata.n_obs} cells, {adata.n_vars} genes")

                # Load loom
                status_text.text("Loading loom file...")
                progress_bar.progress(20)

                temp_loom_path = os.path.join(velocity_analysis_temp_dir, "input.loom")
                with open(temp_loom_path, "wb") as f:
                    f.write(uploaded_loom.read())

                ldata = ad.read_loom(temp_loom_path)
                st.info(f"✓ Loaded loom: {ldata.n_obs} cells, {ldata.n_vars} genes")

                # Check cell overlap
                status_text.text("Checking cell correspondence...")
                progress_bar.progress(25)

                h5ad_cells = set(adata.obs_names)
                loom_cells = set(ldata.obs_names)
                common_cells = h5ad_cells & loom_cells

                st.write(f"""
                **Cell correspondence:**
                - h5ad cells: {len(h5ad_cells)}
                - loom cells: {len(loom_cells)}
                - Common cells: {len(common_cells)}
                """)

                if len(common_cells) == 0:
                    st.error("""
                    ❌ **No matching cells found between h5ad and loom files!**

                    This can happen if:
                    - Cell IDs don't match between files
                    - Wrong loom file was uploaded
                    - Data filtering step was not performed correctly

                    Please check your files and ensure the loom file was generated
                    from the same dataset as the h5ad file.
                    """)
                    st.stop()

                match_percentage = (len(common_cells) / len(h5ad_cells)) * 100
                if match_percentage < 50:
                    st.warning(f"""
                    ⚠️ **Low cell overlap: {match_percentage:.1f}%**

                    Only {len(common_cells)} out of {len(h5ad_cells)} cells in h5ad
                    are present in the loom file. This may indicate:
                    - Incorrect data filtering
                    - Mismatched files
                    - Expected behavior if you intentionally filtered cells
                    """)

                # Filter h5ad to only include cells present in loom
                status_text.text("Filtering h5ad to match loom cells...")
                progress_bar.progress(28)

                adata = adata[list(common_cells), :].copy()
                st.success(f"✓ Filtered h5ad to {adata.n_obs} cells matching loom")

                # Merge
                status_text.text("Merging datasets...")
                progress_bar.progress(30)

                adata = scv.utils.merge(adata, ldata)
                st.success(f"✓ Merged data: {adata.n_obs} cells, {adata.n_vars} genes")

                # Show velocity layer info
                if 'spliced' in adata.layers and 'unspliced' in adata.layers:
                    n_spliced = np.sum(adata.layers['spliced'] > 0)
                    n_unspliced = np.sum(adata.layers['unspliced'] > 0)
                    st.write(f"Spliced counts: {n_spliced}, Unspliced counts: {n_unspliced}")

                # Filter and normalize
                status_text.text("Filtering and normalizing...")
                progress_bar.progress(40)

                # Parse retain genes
                retain_genes_list = parse_gene_list(retain_genes_str)

                # Match gene names case-insensitively
                if retain_genes_list:
                    # Get actual gene names from adata (case-sensitive)
                    all_genes = adata.var_names.tolist()
                    all_genes_upper = [g.upper() for g in all_genes]

                    # Find matching genes
                    retain_genes_matched = []
                    for gene in retain_genes_list:
                        if gene in all_genes_upper:
                            idx = all_genes_upper.index(gene)
                            retain_genes_matched.append(all_genes[idx])

                    st.info(f"ℹ️ Retain genes: {len(retain_genes_matched)} / {len(retain_genes_list)} found")
                    if len(retain_genes_matched) < len(retain_genes_list):
                        missing = set(retain_genes_list) - set([g.upper() for g in retain_genes_matched])
                        st.warning(f"⚠️ Not found: {list(missing)[:10]}")
                else:
                    retain_genes_matched = None

                scv.pp.filter_and_normalize(
                    adata,
                    min_shared_counts=min_shared_counts,
                    n_top_genes=n_top_genes,
                    retain_genes=retain_genes_matched
                )

                if retain_genes_matched:
                    st.success(f"✓ Filtered and normalized (top {n_top_genes} genes + {len(retain_genes_matched)} retained)")
                else:
                    st.success(f"✓ Filtered and normalized (top {n_top_genes} genes)")

                # Moments
                status_text.text("Computing moments...")
                progress_bar.progress(55)

                scv.pp.moments(
                    adata,
                    n_pcs=n_pcs,
                    n_neighbors=n_neighbors
                )
                st.success(f"✓ Computed moments (n_pcs={n_pcs}, n_neighbors={n_neighbors})")

                # Velocity
                status_text.text(f"Computing velocity ({velocity_mode} mode)...")
                progress_bar.progress(70)

                # Run scVelo first
                if velocity_mode == "dynamical":
                    st.warning(f"⏳ Dynamical mode may take several minutes... (using {n_jobs} cores)")
                    try:
                        scv.tl.recover_dynamics(adata, n_jobs=n_jobs)
                    except Exception as parallel_error:
                        if "pickle" in str(parallel_error).lower() or "serialize" in str(parallel_error).lower() or "BrokenProcessPool" in str(type(parallel_error).__name__):
                            st.warning(f"⚠️ Parallel processing failed. Retrying with single core (n_jobs=1)...")
                            scv.tl.recover_dynamics(adata, n_jobs=1)
                        else:
                            raise
                    scv.tl.velocity(adata, mode="dynamical")
                else:
                    scv.tl.velocity(adata, mode=velocity_mode)

                # Save scVelo velocity result
                adata.layers['velocity_scvelo'] = adata.layers['velocity'].copy()
                st.success(f"✓ scVelo velocity computed ({velocity_mode} mode) → saved as 'velocity_scvelo'")

                # Velocity graph
                status_text.text("Computing velocity graph...")
                progress_bar.progress(85)

                scv.tl.velocity_graph(adata)
                st.success("✓ Velocity graph computed")

                # Compute additional metrics BEFORE PAGA (needed for directed PAGA)
                status_text.text("Computing velocity metrics...")
                progress_bar.progress(85)

                # Velocity confidence (also computes velocity_length automatically)
                scv.tl.velocity_confidence(adata)
                st.info("✓ Computed velocity confidence and velocity length")

                # Velocity pseudotime (needed for directed PAGA arrows)
                paga_use_time_prior = False
                try:
                    scv.tl.velocity_pseudotime(adata)
                    st.info("✓ Computed velocity pseudotime")
                    paga_use_time_prior = True
                except Exception as e:
                    st.warning(f"Could not compute velocity pseudotime: {str(e)}")
                    st.info("Will compute PAGA without time prior")

                # PAGA (if requested) - compute AFTER velocity_pseudotime
                if compute_paga:
                    status_text.text("Computing PAGA with directed transitions...")
                    progress_bar.progress(92)

                    if paga_cluster_key is None or paga_cluster_key not in adata.obs.columns:
                        st.warning("⚠️ Selected cluster column not found in data. Skipping PAGA computation.")
                    else:
                        # Compute PAGA with user-selected cluster column (tutorial method)
                        try:
                            # IMPORTANT: Fix for scvelo bug - ensure neighbors structure is correct
                            # This is required by scVelo tutorial
                            if 'neighbors' in adata.uns:
                                if 'distances' in adata.obsp:
                                    adata.uns['neighbors']['distances'] = adata.obsp['distances']
                                if 'connectivities' in adata.obsp:
                                    adata.uns['neighbors']['connectivities'] = adata.obsp['connectivities']
                                st.info("✓ Fixed neighbors structure for PAGA computation")

                            # First compute standard PAGA connectivity with scanpy
                            sc.tl.paga(adata, groups=paga_cluster_key)
                            st.info("✓ Computed PAGA connectivity")

                            # Compute velocity-directed PAGA with scvelo (tutorial method)
                            scv.tl.paga(adata, groups=paga_cluster_key)

                            if 'transitions_confidence' in adata.uns.get('paga', {}):
                                trans_conf = adata.uns['paga']['transitions_confidence']
                                st.success(f"✓ PAGA computed with scvelo (directed transitions)")
                                st.info(f"   Transitions matrix shape: {trans_conf.shape}, max: {trans_conf.max():.4f}")
                            else:
                                st.warning("⚠️ scvelo PAGA computed but transitions_confidence not generated")

                        except Exception as paga_error:
                            error_msg = str(paga_error)
                            st.error(f"❌ Could not compute PAGA: {error_msg}")

                            # Check if this is the known scvelo bug
                            if "mismatching number of index arrays" in error_msg:
                                st.warning("""
                                **⚠️ 既知のscveloバグが検出されました**

                                このエラーは scvelo v0.3.4 以前に存在する既知のバグです。
                                GitHub Issue: https://github.com/theislab/scvelo/issues/1241
                                Fix PR: https://github.com/theislab/scvelo/pull/1308 (まだマージされていません)
                                """)

                                st.info("""
                                **🔧 修正方法:**

                                以下のファイルを編集してバグを修正できます：

                                **ファイル:** `scvelo/tools/paga.py` (約50-52行目)

                                **修正前:**
                                ```python
                                if len(edges) > 0:
                                    return csr_matrix((weights, zip(*edges)), shape=shape)
                                ```

                                **修正後:**
                                ```python
                                if len(edges) > 0:
                                    rows, cols = zip(*edges)
                                    return csr_matrix((weights, (rows, cols)), shape=shape)
                                ```

                                修正後、このアプリを再起動してください。
                                """)

                                # Show file location
                                try:
                                    import scvelo
                                    import os
                                    paga_file = os.path.join(os.path.dirname(scvelo.__file__), 'tools', 'paga.py')
                                    st.code(f"File location:\n{paga_file}", language="text")
                                except:
                                    pass

                            # Show full traceback in expander
                            with st.expander("🔍 Detailed error traceback"):
                                import traceback
                                st.text(traceback.format_exc())

                # Save result
                status_text.text("Saving results...")
                progress_bar.progress(98)

                output_path = os.path.join(velocity_analysis_temp_dir, "velocity_result.h5ad")
                if os.path.exists(output_path):
                    os.remove(output_path)

                adata.write_h5ad(output_path, compression="gzip")

                # Store in session state
                st.session_state.velocity_result_path = output_path
                st.session_state.velocity_adata = adata
                st.session_state.velocity_mode = velocity_mode  # Store mode for filename
                st.session_state.analysis_complete = True

                progress_bar.progress(100)
                status_text.text("Analysis complete!")

                # Build success message
                success_msg = """
                ✅ **Analysis completed successfully!**

                Computed layers and embeddings:
                - Velocity vectors
                - Velocity graph
                - Velocity confidence
                - Velocity pseudotime (if available)
                """

                if compute_paga and 'paga' in adata.uns:
                    success_msg += "- PAGA (partition-based graph abstraction)\n"

                st.success(success_msg)

                # Show summary statistics
                st.subheader("Summary Statistics")

                summary_col1, summary_col2, summary_col3 = st.columns(3)

                with summary_col1:
                    st.metric("Total cells", adata.n_obs)
                    st.metric("Total genes", adata.n_vars)

                with summary_col2:
                    if 'velocity_length' in adata.obs:
                        mean_velocity = adata.obs['velocity_length'].mean()
                        st.metric("Mean velocity length", f"{mean_velocity:.4f}")

                    if 'velocity_confidence' in adata.obs:
                        mean_confidence = adata.obs['velocity_confidence'].mean()
                        st.metric("Mean velocity confidence", f"{mean_confidence:.4f}")

                with summary_col3:
                    if 'velocity_pseudotime' in adata.obs:
                        st.metric("Velocity pseudotime range",
                                f"{adata.obs['velocity_pseudotime'].min():.3f} - {adata.obs['velocity_pseudotime'].max():.3f}")

                # Show available data
                with st.expander("📊 Available data in result", expanded=False):
                    st.write("**Layers:**")
                    st.write(list(adata.layers.keys()))

                    st.write("**Obs (cell metadata):**")
                    st.write(adata.obs.columns.tolist())

                    st.write("**Obsm (embeddings):**")
                    st.write(list(adata.obsm.keys()))

            except Exception as e:
                st.error(f"❌ Error during analysis: {str(e)}")
                st.exception(e)
                st.session_state.analysis_complete = False

            finally:
                progress_bar.empty()
                status_text.empty()

    # ========================================
    # Step 4: Download results
    # ========================================
    if st.session_state.analysis_complete:
        st.header("Step 4: Download results")

        with open(st.session_state.velocity_result_path, "rb") as f:
            result_bytes = f.read()

        # Create filename based on input h5ad and velocity mode
        h5ad_basename = os.path.splitext(uploaded_h5ad.name)[0]
        velocity_mode = st.session_state.get('velocity_mode', 'dynamical')
        output_filename = f"{h5ad_basename}.scVelo.{velocity_mode}.h5ad"

        st.download_button(
            label="⬇️ Download velocity result (h5ad)",
            data=result_bytes,
            file_name=output_filename,
            mime="application/octet-stream",
            type="primary"
        )

        st.info("""
        ### 次のステップ

        ダウンロードしたh5adファイルには以下が含まれています：
        - Velocity vectors (`adata.layers['velocity']`)
        - Velocity graph (`adata.uns['velocity_graph']`)
        - Velocity confidence (`adata.obs['velocity_confidence']`)
        - Velocity pseudotime (`adata.obs['velocity_pseudotime']`)

        #### 推奨: scVelo Visualizationアプリで可視化

        1. **scVelo Visualization** ページに移動
        2. このh5adファイルをアップロード
        3. 以下の可視化が可能：
           - Velocity stream plots
           - Velocity grid plots
           - Phase portraits
           - Velocity embedding plots

        #### Pythonでのvisualization例:
        ```python
        import scvelo as scv
        import scanpy as sc

        # Load result
        adata = sc.read_h5ad('result_velocity.h5ad')

        # Visualizations
        scv.pl.velocity_embedding_stream(adata, basis='umap')
        scv.pl.velocity_embedding_grid(adata, basis='umap')
        scv.pl.velocity_graph(adata)
        scv.pl.velocity(adata, var_names=['gene1', 'gene2'])
        ```

        #### その他のオプション:
        - **CellRank analysis**: 細胞運命と系統解析
        - **Pseudotime gene expression**: 遺伝子発現トレンドの可視化
        """)

else:
    st.info("👆 Please upload both h5ad and loom files to start")
