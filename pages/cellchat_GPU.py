import scanpy as sc
import numpy as np
import streamlit as st
import pandas as pd
import os
import pickle
import re
import time
import matplotlib
matplotlib.use("Agg") # PDFの透明化対策 interactiveではうまく行かない?
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
#from scipy import sparse
import scipy.sparse
#from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
import logging
import traceback
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr
from joblib import Parallel, delayed
from tqdm import tqdm
#print("rpy2 version:", rpy2.__version__)
import scipy
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import plot
from numba import njit, prange
# import scipy.stats as stats
from typing import Union, Optional, List
from streamlit_sortables import sort_items
from matplotlib.patches import FancyArrowPatch
from matplotlib.path import Path
import math
import gc


from pages.cellchat_vis import (
    netVisual_circle, netVisual_circle_individual, netVisual_chord_by_pathway,
    netVisual_aggregate, netVisual_chord, netVisual_heatmap,
    netAnalysis_signalingRole_network, plotGeneExpression)



# ロギング設定を最初に行う
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("CellChat")


def identify_overexpressed_genes_optimized(
    adata,
    group_by: str = None,
    idents_use: list = None,
    invert: bool = False,
    features_name: str = "features",
    only_pos: bool = True,
    features: list = None,
    return_object: bool = True,
    thresh_pct: float = 0,
    thresh_fc: float = 0,
    thresh_p: float = 0.05,
    do_de: bool = True,
    do_fast: bool = True,
    min_cells: int = 10,
    use_gpu: bool = False,
    batch_size: int = 5000,
    n_jobs: int = -1,
):
    """
    高速化した過剰発現遺伝子同定関数
    
    Parameters
    ----------
    [既存のパラメータ]
    use_gpu : bool, optional
        GPUを使用するかどうか（一部の計算のみ対応）
    batch_size : int, optional
        メモリ効率のためのバッチサイズ
    n_jobs : int, optional
        並列処理のジョブ数。-1は全コア使用
        
    Returns
    -------
    dict
        "features": 過剰発現遺伝子のリスト
        "markers_all": 詳細なマーカー情報のデータフレーム
    """
    import pandas as pd
    import numpy as np
    import scipy.sparse
    import time
    from sklearn.utils import parallel_backend
    from joblib import Parallel, delayed
    import os
    
    # 処理開始時間記録
    start_time = time.time()
    
    if not hasattr(pd.DataFrame, "iteritems"):
        pd.DataFrame.iteritems = pd.DataFrame.items
    
    # GPUセットアップ（使用時のみ）
    has_gpu = False
    if use_gpu:
        try:
            import cupy as cp
            from cupyx.scipy import sparse as cusparse
            has_gpu = True
            print("GPU acceleration enabled for compatible operations")
        except ImportError:
            print("GPU acceleration requested but cupy not available, using CPU only")
    
    # n_jobs設定の調整
    if n_jobs == -1:
        n_jobs = os.cpu_count() or 4
    
    # 自動変換を有効化
    numpy2ri.activate()
    pandas2ri.activate()
    presto = rpackages.importr("presto")
    
    # データ行列の準備 - バッチ処理による最適化
    if scipy.sparse.issparse(adata.X):
        if has_gpu:
            try:
                # スパース行列をCuPyでGPUに転送し処理
                X_gpu = cusparse.csr_matrix(adata.X)
                X = cp.asnumpy(X_gpu.toarray())
                del X_gpu  # GPUメモリ解放
                cp.get_default_memory_pool().free_all_blocks()
            except Exception as e:
                print(f"GPU sparse matrix conversion failed: {str(e)}, using CPU")
                # 大規模行列はバッチ処理
                if adata.X.shape[0] > batch_size:
                    X = np.vstack([adata.X[i:min(i+batch_size, adata.X.shape[0])].toarray() 
                                  for i in range(0, adata.X.shape[0], batch_size)])
                else:
                    X = adata.X.toarray()
        else:
            # メモリ効率のためにバッチ処理
            if adata.X.shape[0] > batch_size:
                X = np.vstack([adata.X[i:min(i+batch_size, adata.X.shape[0])].toarray() 
                              for i in range(0, adata.X.shape[0], batch_size)])
            else:
                X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    
    # セルラベルの取得 - 最適化なし
    labels = adata.obs[group_by] if group_by else adata.obs.iloc[:, 0]
    labels = pd.Categorical(labels)
    
    # 遺伝子リストの取得
    features_use = list(adata.var_names) if features is None else list(set(features) & set(adata.var_names))
    
    # 高速インデックス変換のための辞書作成
    gene_to_idx = {gene: idx for idx, gene in enumerate(adata.var_names)}
    feature_indices = [gene_to_idx[gene] for gene in features_use]
    
    # 過剰発現解析用の行列を作成 - 転置処理の最適化
    X_gene = X.T  # 転置は全体で一度だけ実行
    data_use = X_gene[feature_indices, :]
    
    # 各クラスタのラベル
    level_use = list(labels.categories)
    
    if do_de:
        if do_fast:
            # 前処理の経過時間記録
            prep_time = time.time() - start_time
            print(f"Preprocessing completed in {prep_time:.2f}s")
            
            # データフレーム変換の最適化
            # Rとの変換を効率化するためPandasを直接使用
            data_use_df = pd.DataFrame(data_use, index=features_use, columns=adata.obs_names)
            
            # R変換の最適化
            with localconverter(ro.default_converter + pandas2ri.converter):
                # labelsとgroups_useのR変換
                r_labels = ro.FactorVector(labels.astype(str).tolist())
                r_groups = ro.StrVector(level_use)
                
                # presto呼び出し - これはR依存なので高速化余地は少ない
                print(f"Running wilcoxauc on {len(features_use)} genes and {len(level_use)} groups")
                r_start = time.time()
                res = presto.wilcoxauc(data_use_df, r_labels, groups_use=r_groups)
                r_time = time.time() - r_start
                print(f"wilcoxauc completed in {r_time:.2f}s")
                
                # 結果変換の最適化
                genes_de = pd.DataFrame(ro.conversion.rpy2py(res))
            
            # 後処理の最適化 - インプレース操作と効率的なフィルタリング
            # 列名変更
            genes_de.rename(columns={
                "group": "clusters",
                "feature": "features",
                "pval": "pvalues",
                "logFC": "logFC",
                "pct_in": "pct.1",
                "pct_out": "pct.2",
                "padj": "pvalues.adj"
            }, inplace=True)
            
            # 効率的な演算
            genes_de["logFC_abs"] = genes_de["logFC"].abs()
            genes_de["pct.max"] = genes_de[["pct.1", "pct.2"]].max(axis=1)
            
            # マスクベースのフィルタリング - 単一の演算で効率的に
            mask = ((genes_de["pvalues"] < thresh_p) & 
                   (genes_de["logFC_abs"] >= thresh_fc) & 
                   (genes_de["pct.max"] > thresh_pct))
            markers_all = genes_de[mask].sort_values("pvalues")
            
        else:
            # 標準Python実装の最適化 - 並列処理を適用
            def process_cluster(cluster, features_use, data_use, labels):
                """単一クラスターの処理を並列化するための関数"""
                cluster_markers = []
                cluster_mask = (labels == cluster)
                cluster_data = data_use[:, cluster_mask]
                other_data = data_use[:, ~cluster_mask]
                
                # バッチ処理で効率化
                batch_size = 500  # 遺伝子をバッチで処理
                for b_start in range(0, len(features_use), batch_size):
                    b_end = min(b_start + batch_size, len(features_use))
                    batch_genes = features_use[b_start:b_end]
                    batch_indices = list(range(b_start, b_end))
                    
                    # バッチ内の統計量を効率的に計算
                    pct1_batch = np.mean(cluster_data[batch_indices, :] > 0, axis=1)
                    pct2_batch = np.mean(other_data[batch_indices, :] > 0, axis=1)
                    
                    # Mann-Whitney U テストはまとめて計算できないので、遺伝子ごとに
                    for i, gene_idx in enumerate(batch_indices):
                        gene = batch_genes[i]
                        pct1 = pct1_batch[i]
                        pct2 = pct2_batch[i]
                        
                        # 統計テスト
                        stat, p_val = scipy.stats.mannwhitneyu(
                            cluster_data[gene_idx, :], 
                            other_data[gene_idx, :], 
                            alternative='two-sided'
                        )
                        
                        # 発現差の計算
                        log_fc = np.log(
                            np.mean(cluster_data[gene_idx, :] + 1e-10) / 
                            np.mean(other_data[gene_idx, :] + 1e-10)
                        )
                        
                        cluster_markers.append({
                            "clusters": cluster,
                            "features": gene,
                            "pvalues": p_val,
                            "logFC": log_fc,
                            "pct.1": pct1,
                            "pct.2": pct2
                        })
                
                return cluster_markers
            
            # 並列処理で各クラスターを処理
            print(f"Running parallel Mann-Whitney U tests for {len(level_use)} clusters")
            mwu_start = time.time()
            with parallel_backend("loky", n_jobs=n_jobs):
                all_markers = Parallel()(
                    delayed(process_cluster)(
                        cluster, features_use, data_use, labels
                    ) for cluster in level_use
                )
            mwu_time = time.time() - mwu_start
            print(f"Mann-Whitney U tests completed in {mwu_time:.2f}s")
            
            # 結果をフラット化して1つのデータフレームに
            markers = [marker for cluster_markers in all_markers for marker in cluster_markers]
            markers_all = pd.DataFrame(markers)
            
            # 後処理
            markers_all["logFC_abs"] = markers_all["logFC"].abs()
            markers_all["pct.max"] = markers_all[["pct.1", "pct.2"]].max(axis=1)
            mask = ((markers_all["pvalues"] < thresh_p) & 
                   (markers_all["logFC_abs"] >= thresh_fc) & 
                   (markers_all["pct.max"] > thresh_pct))
            markers_all = markers_all[mask].sort_values("pvalues")
        
        features_sig = markers_all["features"].unique()
    else:
        # do.DE=False モードの最適化 - 効率的な集計
        if has_gpu and data_use.size > 1_000_000:  # 大規模データのみGPUを使用
            try:
                # CuPyを使った高速集計
                data_use_gpu = cp.array(data_use)
                nCells = cp.sum(data_use_gpu > 0, axis=1).get()
                del data_use_gpu
                cp.get_default_memory_pool().free_all_blocks()
            except Exception as e:
                print(f"GPU counting failed: {str(e)}, using CPU")
                nCells = (data_use > 0).sum(axis=1)
        else:
            # 標準的な集計 - NumPyは十分に最適化済み
            nCells = (data_use > 0).sum(axis=1)
        
        # 効率的なデータフレーム作成
        markers_all = pd.DataFrame({
            "features": features_use,
            "nCells": nCells
        })
        markers_all = markers_all[markers_all["nCells"] >= min_cells]
        features_sig = markers_all["features"].unique()
    
    # 総実行時間の計算
    total_time = time.time() - start_time
    print(f"identify_overexpressed_genes completed in {total_time:.2f}s")
    print(f"{len(features_sig)} genes passed filtering")
    
    return {"features": features_sig, "markers_all": markers_all}

def cellchat_analysis_optimized(
    adata,
    groupby,
    db,
    use_layer=None,
    min_cells=10,
    expr_prop=0.1,
    pseudocount=1.0,
    trim_threshold=0.05,
    k=0.5,
    n=1,
    type_mean="triMean",
    raw_use=True,
    population_size=False,
    nboot=100,
    seed=12345,
    n_jobs=8,
    key_added="cellchat_res",
    trim=0.1,
    apply_pval_filter=True,
    features=None,
    optimize_memory=True,   # 記憶域最適化の有効化
    use_gpu=True,          # GPUを使用するかどうか
    gpu_precision='float64', # GPU計算の精度 ('float64' 推奨)
    validate_gpu_results=True, # CPU結果とGPU結果を比較検証するか
    validation_tolerance=1e-5, # 検証時の許容誤差
    gpu_memory_limit=0.8,    # GPU使用メモリの上限（0-1の割合）
    hybrid_mode=True         # ハイブリッドモード（一部の演算のみGPU使用）
):
    """
    GPUサポートを追加した、最適化されたCellChatアルゴリズム
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    groupby : str
        Column name in .obs containing cell type/cluster information
    db : dict
        CellChatDB dictionary (obtained from get_cellchatdb_from_r())
    
    [従来のパラメータは同様]
    
    optimize_memory : bool, optional
        Whether to use memory optimization techniques. Default is True
    use_gpu : bool, optional
        GPUを使用するかどうか。Default is False
    gpu_precision : str, optional
        GPU計算の精度 ('float64' または 'float32')。Default is 'float64'
    validate_gpu_results : bool, optional
        GPU計算結果をCPU結果と比較して検証するか。Default is True
    validation_tolerance : float, optional
        GPU/CPU結果比較時の許容誤差。Default is 1e-5
    gpu_memory_limit : float, optional
        GPUメモリ使用量の上限（0-1の割合）。Default is 0.8
    hybrid_mode : bool, optional
        ハイブリッドモード（一部の演算のみGPU使用）。Default is True
        
    Returns
    -------
    dict
        Cell-cell communication analysis results, identical to the original cellchat_analysis
    """

    has_gpu = False
    cp = None
    
    if use_gpu:
        try:
            import cupy as cp
            has_gpu = True
            # GPUメモリ使用量を制限
            if gpu_memory_limit > 0 and gpu_memory_limit <= 1.0:
                mempool = cp.get_default_memory_pool()
                mempool.set_limit(cp.cuda.Device().mem_info[1] * gpu_memory_limit)
            # 精度の設定
            gpu_dtype = cp.float64 if gpu_precision == 'float64' else cp.float32
            logger_info = f"GPU acceleration enabled using {gpu_precision} precision"
            st.write("GPU available")
        except ImportError:
            logger_info = "GPU acceleration requested but cupy not available, falling back to CPU"
            use_gpu = False
            st.write("GPU is not ready")
    else:
        logger_info = "Using CPU for calculations (GPU not enabled)"
    
    # Initialize logger
    logger = logging.getLogger("CellChat-GPU-Optimized")
    logger.info(logger_info)
    
    try:
        start_time = time.time()
        
        # Progress display function
        if 'st' in globals():
            try:
                progress_bar = st.progress(0)
                status_area = st.empty()
                
                def show_progress(progress, message):
                    progress_bar.progress(progress)
                    status_area.text(message)
            except:
                def show_progress(progress, message):
                    logger.info(f"{int(progress * 100)}% - {message}")
        else:
            def show_progress(progress, message):
                if int(progress * 100) % 10 == 0:  # Only log every 10%
                    logger.info(f"{int(progress * 100)}% - {message}")
        
        show_progress(0.0, "Starting CellChat analysis...")
        
        # GPU結果検証用のヘルパー関数
        def validate_results(cpu_result, gpu_result, name=""):
            """GPU計算結果とCPU計算結果を比較検証する"""
            if not has_gpu or not validate_gpu_results:
                return True
            
            # CPUデータに変換
            if hasattr(gpu_result, 'get'):
                gpu_result_cpu = gpu_result.get()
            else:
                gpu_result_cpu = gpu_result
                
            # データ形状が同じか確認
            if isinstance(cpu_result, np.ndarray) and isinstance(gpu_result_cpu, np.ndarray):
                if cpu_result.shape != gpu_result_cpu.shape:
                    logger.warning(f"Shape mismatch in {name}: CPU {cpu_result.shape} vs GPU {gpu_result_cpu.shape}")
                    return False
                
                # 最大差異を計算
                abs_diff = np.abs(cpu_result - gpu_result_cpu)
                max_diff = np.max(abs_diff) if abs_diff.size > 0 else 0
                mean_diff = np.mean(abs_diff) if abs_diff.size > 0 else 0
                
                # 許容範囲内か確認
                if max_diff > validation_tolerance:
                    logger.warning(f"Results diverge in {name}: Max diff {max_diff}, Mean diff {mean_diff}")
                    if max_diff > validation_tolerance * 10:
                        logger.warning("Large difference detected, falling back to CPU result")
                        return False
                    logger.info(f"Difference within acceptable range, continuing with GPU result")
                
                logger.info(f"Validation successful for {name}: Max diff {max_diff}, Mean diff {mean_diff}")
                return True
            
            return True  # 比較できない場合はTrueを返す
        
        # Set random seed - exact same as original
        np.random.seed(seed)
        if has_gpu:
            cp.random.seed(seed)
                
        # Import GPU-optimized helper functions if available
        if has_gpu and use_gpu:
            try:
                from pages.gpu_cellchat_helpers import (
                    compute_hill_outer_gpu,
                    precompute_gene_expressions_gpu,
                    geometricMean_gpu
                )
            except ImportError:
                logger.warning("GPU helper functions not available, falling back to CPU implementations")
        show_progress(0.05, "Filgering genes...")
        # Preprocess data - exactly as original
        if features is not None:
            logger.info(f"Using provided gene list: {len(features)} genes")
            features = [f for f in features if f in adata.var_names]
            adata_filtered = adata[:, features].copy()
        else:
            adata_filtered = preprocess_data(adata, groupby, min_cells=min_cells, thresh_pct=expr_prop)
        
        show_progress(0.1, "Data preprocessing complete")
        
        # Get database components - exact same references
        resource = db.get('interaction', pd.DataFrame())
        complex_input = db.get('complex', pd.DataFrame())
        cofactor_input = db.get('cofactor', pd.DataFrame())
        gene_info = db.get('geneInfo', pd.DataFrame())
        
        # Basic validation - same as original
        if resource.empty:
            raise ValueError("Empty interaction information in CellChatDB.")
        for required_col in ['ligand', 'receptor']:
            if required_col not in resource.columns:
                raise ValueError(f"Required column '{required_col}' not found in resource dataframe")
        
        # Get expression data - exact same approach
        if use_layer is not None and use_layer in adata_filtered.layers:
            logger.info(f"Using layer '{use_layer}'")
            X = adata_filtered.layers[use_layer]
        else:
            logger.info("Using default X matrix")
            X = adata_filtered.X
        
        # Convert sparse matrix to dense - efficient approach with GPU options
        if scipy.sparse.issparse(X):
            data_size = X.shape[0] * X.shape[1]
            
            # GPUを使用する場合、大きいスパース行列の場合はGPUで直接変換
            if has_gpu and use_gpu and not hybrid_mode and data_size < cp.cuda.Device().mem_info[0] * 0.3:
                try:
                    import cupyx.scipy.sparse as cusparse
                    logger.info("Converting sparse matrix on GPU")
                    X_gpu = cusparse.csr_matrix(X)
                    X = cp.asnumpy(X_gpu.toarray())
                    del X_gpu
                    cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    st.warning(f"GPU sparse conversion failed: {str(e)}, falling back to CPU")
                    if optimize_memory and X.shape[0] > 5000:
                        chunk_size = 5000
                        result = []
                        for i in range(0, X.shape[0], chunk_size):
                            end = min(i + chunk_size, X.shape[0])
                            result.append(X[i:end].toarray())
                        X = np.vstack(result)
                    else:
                        X = X.toarray()
            else:
                # CPUで変換
                if optimize_memory and X.shape[0] > 5000:
                    chunk_size = 5000
                    result = []
                    for i in range(0, X.shape[0], chunk_size):
                        end = min(i + chunk_size, X.shape[0])
                        result.append(X[i:end].toarray())
                    X = np.vstack(result)
                else:
                    X = X.toarray()
        
        # Get cell type labels - exact same process
        cell_labels = adata_filtered.obs[groupby].copy()
        cell_types = np.array(sorted(cell_labels.unique()))
        
        logger.info(f"Number of cell types: {len(cell_types)}")
        if len(cell_types) < 2:
            raise ValueError(f"Only {len(cell_types)} cell type found. At least 2 required.")
        
        # 正規化 - CPU版とGPU版の両方を実装
        # CPU版のデータを常に保持
        data_use_cpu = X / np.max(X).astype(np.float64)
        
        if has_gpu and use_gpu and not hybrid_mode:
            # GPU版
            X_gpu = cp.array(X, dtype=gpu_dtype)
            max_val_gpu = cp.max(X_gpu).astype(gpu_dtype)
            data_use_gpu = X_gpu / max_val_gpu
            
            # 検証
            if validate_gpu_results:
                valid = validate_results(data_use_cpu, data_use_gpu, "data normalization")
                if not valid:
                    # GPUの結果を使用しない場合
                    data_use = data_use_cpu
                    del X_gpu, data_use_gpu
                    cp.get_default_memory_pool().free_all_blocks()
                else:
                    data_use = data_use_gpu
            else:
                data_use = data_use_gpu
        else:
            # CPU版のみ
            data_use = data_use_cpu.copy()
        
        nC = data_use.shape[0] if isinstance(data_use, np.ndarray) else data_use.shape[0]
        
        # Free original data for large datasets but keep data_use and data_use_cpu
        if optimize_memory:
            del X
            gc.collect()
            if has_gpu and use_gpu and not hybrid_mode:
                # Keep only necessary GPU arrays
                if 'X_gpu' in locals():
                    del X_gpu
                cp.get_default_memory_pool().free_all_blocks()
        
        # Define mean function - GPU対応バージョン
        def FunMean_gpu(x):
            # GPU用の実装
            if has_gpu and use_gpu and isinstance(x, cp.ndarray):
                x_no_nan = x[~cp.isnan(x)]
                if len(x_no_nan) == 0:
                    return cp.nan
                q = cp.quantile(x_no_nan, cp.array([0.25, 0.5, 0.5, 0.75]))
                return cp.mean(q)
            else:
                # CPUにフォールバック
                x_no_nan = x[~np.isnan(x)]
                if len(x_no_nan) == 0:
                    return np.nan
                q = np.quantile(x_no_nan, [0.25, 0.5, 0.5, 0.75], method='linear')
                return np.mean(q)
        
        # 元のFunMean実装を保持
        if type_mean == "triMean":
            def FunMean(x):
                x_no_nan = x[~np.isnan(x)]
                if len(x_no_nan) == 0:
                    return np.nan
                q = np.quantile(x_no_nan, [0.25, 0.5, 0.5, 0.75], method='linear')
                return np.mean(q)
        elif type_mean == "truncatedMean":
            from scipy.stats import trim_mean
            def FunMean(x):
                return trim_mean(x, proportiontocut=trim)
        elif type_mean == "median":
            def FunMean(x):
                return np.median(x)
        else:
            def FunMean(x):
                return np.mean(x)
        
        show_progress(0.15, "Calculating cell type averages...")
        
        # Optimization: Precompute cell indexes to avoid repeated searches
        cell_counts = {}
        cell_type_indices = {}
        for cell_type in cell_types:
            indices = np.where(cell_labels == cell_type)[0]
            cell_counts[cell_type] = len(indices)
            cell_type_indices[cell_type] = indices
        
        # Identify valid cell types
        valid_cell_types = [ct for ct in cell_types if cell_counts.get(ct, 0) >= min_cells]
        if len(valid_cell_types) < 2:
            raise ValueError(f"Less than 2 cell types have at least {min_cells} cells")
            
        # GPU対応の平均発現計算
        data_use_avg_dict = {}
        
        # GPUを使用する場合
        if has_gpu and use_gpu and not hybrid_mode and isinstance(data_use, cp.ndarray):
            # 各細胞タイプごとに平均発現を計算（GPU版）
            for cell_type in valid_cell_types:
                indices = cell_type_indices[cell_type]
                indices_gpu = cp.array(indices)
                
                # 大きなデータセットの場合はバッチ処理
                if optimize_memory and data_use.shape[1] > 10000:
                    gene_chunk_size = 1000
                    n_genes = data_use.shape[1]
                    avg_expr = cp.zeros(n_genes, dtype=gpu_dtype)
                    
                    for i in range(0, n_genes, gene_chunk_size):
                        end = min(i + gene_chunk_size, n_genes)
                        data_subset = data_use[indices_gpu, i:end]
                        # GPU用の平均計算関数を適用
                        temp_result = cp.apply_along_axis(FunMean_gpu, 0, data_subset)
                        avg_expr[i:end] = temp_result
                else:
                    data_subset = data_use[indices_gpu]
                    avg_expr = cp.apply_along_axis(FunMean_gpu, 0, data_subset)
                
                # 結果をCPUに戻す
                if validate_gpu_results:
                    # CPU版でも計算して検証
                    cpu_indices = np.array(indices)
                    cpu_data_subset = data_use_cpu[cpu_indices]
                    cpu_avg_expr = np.apply_along_axis(FunMean, 0, cpu_data_subset)
                    
                    # 検証
                    valid = validate_results(cpu_avg_expr, avg_expr, f"cell type {cell_type} mean expression")
                    if not valid:
                        st.warning(f"Using CPU result for cell type {cell_type}")
                        data_use_avg_dict[cell_type] = cpu_avg_expr
                    else:
                        data_use_avg_dict[cell_type] = cp.asnumpy(avg_expr)
                        st.write("Using GPU results...")
                else:
                    data_use_avg_dict[cell_type] = cp.asnumpy(avg_expr)
        else:
            # CPU版の計算
            for cell_type in valid_cell_types:
                indices = cell_type_indices[cell_type]
                
                # 最適化: 大きなデータセットは分割処理
                if optimize_memory and data_use_cpu.shape[1] > 10000:
                    gene_chunk_size = 1000
                    n_genes = data_use_cpu.shape[1]
                    avg_expr = np.zeros(n_genes)
                    
                    for i in range(0, n_genes, gene_chunk_size):
                        end = min(i + gene_chunk_size, n_genes)
                        data_subset = data_use_cpu[indices, i:end]
                        avg_expr[i:end] = np.apply_along_axis(FunMean, 0, data_subset)
                else:
                    data_subset = data_use_cpu[indices]
                    avg_expr = np.apply_along_axis(FunMean, 0, data_subset)
                
                data_use_avg_dict[cell_type] = avg_expr
        
        # Convert to DataFrame
        data_use_avg_df = pd.DataFrame(data_use_avg_dict, index=adata_filtered.var_names)
        
        show_progress(0.2, "Computing ligand-receptor expression...")
        
        # Compute ligand-receptor expression levels - original function for consistency
        gene_to_index = {gene: i for i, gene in enumerate(adata_filtered.var_names)}
        
        # Use original functions for exact results
        dataLavg = computeExpr_LR(resource['ligand'].values, data_use_avg_df, complex_input)
        dataRavg = computeExpr_LR(resource['receptor'].values, data_use_avg_df, complex_input)
        
        # Account for co-receptors - using exact same functions
        dataRavg_co_A_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "A")
        dataRavg_co_I_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "I")
        dataRavg = dataRavg * dataRavg_co_A_receptor / dataRavg_co_I_receptor
        
        # Apply population size effect if enabled - exact same calculation
        if population_size:
            # Calculate proportion of cells for each cell type
            cell_proportions = np.array([np.sum(cell_labels == ct) for ct in valid_cell_types]) / nC
            # Use same proportion for all LR pairs
            dataLavg2 = np.tile(cell_proportions, (len(resource), 1))
            dataRavg2 = dataLavg2
        else:
            dataLavg2 = np.ones((len(resource), len(valid_cell_types)))
            dataRavg2 = np.ones((len(resource), len(valid_cell_types)))
        
        # Identify agonist and antagonist indices - exactly as original
        index_agonist = np.where(resource['agonist'].notna() & (resource['agonist'] != ""))[0] if 'agonist' in resource.columns else []
        index_antagonist = np.where(resource['antagonist'].notna() & (resource['antagonist'] != ""))[0] if 'antagonist' in resource.columns else []
        
        show_progress(0.3, "Preparing permutation data...")
        
        # Prepare permutation data with GPU support
        permutation = np.zeros((nC, nboot), dtype=np.int32)
        
        if has_gpu and use_gpu and not hybrid_mode:
            # GPUでの生成
            for i in range(nboot):
                # GPUで乱数生成し、CPUに転送
                perm_gpu = cp.random.permutation(nC)
                permutation[:, i] = cp.asnumpy(perm_gpu)
        else:
            # CPU版
            for i in range(nboot):
                permutation[:, i] = np.random.permutation(nC)
        
        # Precompute gene expressions with GPU support
        if nboot > 0:
            show_progress(0.35, "Precomputing gene expression for permutations...")
            # Always pass data_use_cpu to precompute_gene_expressions as a fallback
            data_for_precompute = data_use_cpu
            
            # Choose between GPU and CPU implementation
            if has_gpu and use_gpu and not hybrid_mode:
                try:
                    # GPU実装を使用
                    all_gene_expr = precompute_gene_expressions_gpu(
                        data_for_precompute, 
                        cell_labels, 
                        permutation, 
                        valid_cell_types, 
                        FunMean_gpu, 
                        nboot, 
                        n_jobs,
                        gpu_dtype=gpu_dtype
                    )
                    
                    # 検証（サンプルデータで）
                    if validate_gpu_results:
                        # CPUでも一部計算して比較
                        sample_size = min(10, nboot)
                        cpu_result = precompute_gene_expressions(
                            data_for_precompute,
                            cell_labels,
                            permutation[:, :sample_size],
                            valid_cell_types,
                            FunMean,
                            sample_size,
                            n_jobs
                        )
                        
                        # 最初の数個のブートストラップ結果だけ比較
                        valid = validate_results(
                            cpu_result[:, :, :3], 
                            all_gene_expr[:, :, :3], 
                            "gene expression precomputation"
                        )
                        
                        if not valid:
                            st.warning("GPU precomputation results diverged, falling back to CPU")
                            st.warning("Bootstrap GPU analysis is not valid. Using CPU")
                            # 全体をCPUで再計算
                            all_gene_expr = precompute_gene_expressions(
                                data_for_precompute,
                                cell_labels,
                                permutation,
                                valid_cell_types,
                                FunMean,
                                nboot,
                                n_jobs
                            )
                        else:
                            st.write("GPU analysis is validated.")
                except Exception as e:
                    st.warning(f"GPU gene expression precomputation failed: {str(e)}")
                    # CPUにフォールバック
                    st.warning("GPU analysis is not valid. Using CPU")
                    all_gene_expr = precompute_gene_expressions(
                        data_for_precompute,
                        cell_labels,
                        permutation,
                        valid_cell_types,
                        FunMean,
                        nboot,
                        n_jobs
                    )
            else:
                # CPU版を使用
                st.warning("GPU analysis is not valid. Using CPU")
                all_gene_expr = precompute_gene_expressions(
                    data_for_precompute,
                    cell_labels,
                    permutation,
                    valid_cell_types,
                    FunMean,
                    nboot,
                    n_jobs
                )
        
        # 元のdata_useとdata_use_cpuが不要なら解放（all_gene_exprを計算した後）
        if optimize_memory:
            if use_gpu:
                if isinstance(data_use, cp.ndarray):
                    del data_use
                    cp.get_default_memory_pool().free_all_blocks()
            else:
                del data_use
            # data_use_cpuはまだ必要な場合があるので保持
            gc.collect()
        
        show_progress(0.4, "Setting up interaction matrices...")
        
        # Prepare complex mapping in advance - optimization that doesn't change results
        complex_mapping = {}
        if not complex_input.empty:
            for complex_name in complex_input.index:
                # Get complex subunits
                subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[complex_name, subunits_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in gene_to_index]
                
                if subunits:
                    # Store gene indices for subunits
                    complex_mapping[complex_name] = [gene_to_index[s] for s in subunits]
        
        # Initialize probability and p-value matrices
        numCluster = len(valid_cell_types)
        nLR = len(resource)
        Prob = np.zeros((numCluster, numCluster, nLR))
        Pval = np.zeros((numCluster, numCluster, nLR))
        
        # Get ligand and receptor gene indices - optimization that doesn't change results
        ligand_indices = []
        receptor_indices = []
        
        for i in range(nLR):
            ligand = resource['ligand'].iloc[i]
            receptor = resource['receptor'].iloc[i]
            
            # Determine if single gene or complex
            if isinstance(ligand, str) and ligand in gene_to_index:
                # Single gene
                ligand_indices.append((i, [gene_to_index[ligand]], False))
            elif isinstance(ligand, str) and ligand in complex_mapping:
                # Complex
                ligand_indices.append((i, complex_mapping[ligand], True))
            else:
                # Unknown gene/complex
                ligand_indices.append((i, [], None))
            
            if isinstance(receptor, str) and receptor in gene_to_index:
                receptor_indices.append((i, [gene_to_index[receptor]], False))
            elif isinstance(receptor, str) and receptor in complex_mapping:
                receptor_indices.append((i, complex_mapping[receptor], True))
            else:
                receptor_indices.append((i, [], None))
        
        show_progress(0.45, "Computing interaction probabilities...")
        
        # Process LR pairs with GPU optimization when possible
        # Set up parallel processing for batch computing
        from joblib import Parallel, delayed
        
        # Function for computing permutation batch - GPU-hybrid approach
        def compute_permutation_batch(batch_indices, ligand_info, receptor_info, numCluster, n, k, 
                                      population_size, cell_labels, permutation, valid_cell_types, nC,
                                      use_gpu_compute=False):
            batch_results = np.zeros((numCluster * numCluster, len(batch_indices)))
            
            for idx, j in enumerate(batch_indices):
                # Get ligand expression
                lr_i, ligand_gene_indices, is_ligand_complex = ligand_info
                
                if not is_ligand_complex:
                    # Single gene
                    if ligand_gene_indices:
                        ligand_idx = ligand_gene_indices[0]
                        dataLavgB = all_gene_expr[ligand_idx, :, j].reshape(1, -1)
                    else:
                        dataLavgB = np.zeros((1, numCluster))
                else:
                    # Complex - calculate geometric mean
                    expr_values = np.array([all_gene_expr[idx, :, j] for idx in ligand_gene_indices])
                    if len(expr_values) > 0:
                        # Log transform, average, then transform back
                        log_values = np.log(expr_values + 1e-10)
                        dataLavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                    else:
                        dataLavgB = np.zeros((1, numCluster))
                
                # Get receptor expression
                lr_i, receptor_gene_indices, is_receptor_complex = receptor_info
                
                if not is_receptor_complex:
                    if receptor_gene_indices:
                        receptor_idx = receptor_gene_indices[0]
                        dataRavgB = all_gene_expr[receptor_idx, :, j].reshape(1, -1)
                    else:
                        dataRavgB = np.zeros((1, numCluster))
                else:
                    expr_values = np.array([all_gene_expr[idx, :, j] for idx in receptor_gene_indices])
                    if len(expr_values) > 0:
                        log_values = np.log(expr_values + 1e-10)
                        dataRavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                    else:
                        dataRavgB = np.zeros((1, numCluster))
                
                # Calculate interaction probability - using GPU if available and requested
                if use_gpu_compute and has_gpu and use_gpu:
                    try:
                        dataLavgB_gpu = cp.array(dataLavgB[0, :], dtype=gpu_dtype)
                        dataRavgB_gpu = cp.array(dataRavgB[0, :], dtype=gpu_dtype)
                        
                        # GPU calculation
                        dataLRB_gpu = cp.outer(dataLavgB_gpu, dataRavgB_gpu)
                        P1_boot_gpu = dataLRB_gpu**n / (k**n + dataLRB_gpu**n)
                        
                        # Convert back to CPU
                        P1_boot = cp.asnumpy(P1_boot_gpu)
                        
                        # Cleanup
                        del dataLavgB_gpu, dataRavgB_gpu, dataLRB_gpu, P1_boot_gpu
                        # Free memory
                        if idx % 10 == 0:
                            cp.get_default_memory_pool().free_all_blocks()
                    except Exception as e:
                        # Fallback to CPU on error
                        st.warning(f"GPU calculation failed in permutation {j}: {str(e)}. Using CPU.")
                        dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                        P1_boot = dataLRB**n / (k**n + dataLRB**n)
                else:
                    # CPU calculation
                    dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                    P1_boot = dataLRB**n / (k**n + dataLRB**n)
                
                # Use simplified approach for permutation (same as original)
                P2_boot = np.ones((numCluster, numCluster))
                P3_boot = np.ones((numCluster, numCluster))
                
                # Population size effect
                P4_boot = np.ones((numCluster, numCluster))
                if population_size:
                    group_boot = cell_labels.values[permutation[:, j]]
                    cell_proportions_boot = np.array([np.sum(group_boot == ct) for ct in valid_cell_types]) / nC
                    P4_boot = np.outer(cell_proportions_boot, cell_proportions_boot)
                
                # Final probability
                Pboot_result = P1_boot * P2_boot * P3_boot * P4_boot
                batch_results[:, idx] = Pboot_result.flatten()
            
            return batch_results
        
        # Process each LR pair with GPU optimization for hill function calculation
        for i in range(nLR):
            # Show periodic progress without spamming
            if i % max(1, nLR // 10) == 0:
                show_progress(0.45 + 0.3 * i / nLR, f"Processing LR pairs {i+1}/{nLR}...")
                
            # Hill function to compute probability - GPU or CPU implementation
            if has_gpu and use_gpu and not hybrid_mode:
                try:
                    # GPU implementation
                    dataLavg_gpu = cp.array(dataLavg[i, :], dtype=gpu_dtype)
                    dataRavg_gpu = cp.array(dataRavg[i, :], dtype=gpu_dtype)
                    
                    # Call GPU-optimized function
                    P1_gpu = compute_hill_outer_gpu(dataLavg_gpu, dataRavg_gpu, k, n)
                    
                    # Validate with CPU
                    if validate_gpu_results:
                        P1_cpu = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
                        valid = validate_results(P1_cpu, P1_gpu, f"hill function for LR pair {i}")
                        
                        if valid:
                            P1 = cp.asnumpy(P1_gpu)
                        else:
                            P1 = P1_cpu
                    else:
                        P1 = cp.asnumpy(P1_gpu)
                    
                    # Cleanup
                    del dataLavg_gpu, dataRavg_gpu, P1_gpu
                    if i % 50 == 0:  # 定期的にメモリ解放
                        cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    logger.warning(f"GPU Hill function failed for LR pair {i}: {str(e)}. Using CPU.")
                    P1 = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
            else:
                # CPU implementation
                P1 = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
            
            # Agonist effect - use original function
            P2 = np.ones((numCluster, numCluster))
            if i in index_agonist:
                data_agonist = computeExpr_agonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P2 = np.outer(data_agonist, data_agonist)
            
            # Antagonist effect - use original function
            P3 = np.ones((numCluster, numCluster))
            if i in index_antagonist:
                data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P3 = np.outer(data_antagonist, data_antagonist)
            
            # Population size effect - exact same calculation
            P4 = np.ones((numCluster, numCluster))
            if population_size:
                P4 = np.outer(dataLavg2[i, :], dataRavg2[i, :])
            
            # Final probability - exact same calculation
            Pnull = P1 * P2 * P3 * P4
            Prob[:, :, i] = Pnull
            
            # Skip p-value calculation if no interaction
            if np.sum(Pnull) == 0 or nboot == 0:
                Pval[:, :, i] = 1
                continue
            
            Pnull_vec = Pnull.flatten()
            
            # Get ligand and receptor info
            ligand_info = ligand_indices[i]
            receptor_info = receptor_indices[i]
            
            # Skip if expression can't be retrieved
            if ligand_info[2] is None or receptor_info[2] is None:
                Pval[:, :, i] = 1
                continue
            
            # Initialize bootstrap probabilities
            Pboot = np.zeros((numCluster * numCluster, nboot))
            
            # Optimization: Use batch processing
            batch_size = min(20, nboot)  # Adjust batch size for efficiency 
            
            # Process in batches with parallel computing
            n_jobs_effective = min(n_jobs, os.cpu_count() or 1, nboot)
            
            # 各バッチでのGPU使用を制御
            # ヒートメモリ問題を避けるため、一部のバッチのみGPU使用
            use_gpu_for_batches = has_gpu and use_gpu and not hybrid_mode
            
            for b_start in range(0, nboot, batch_size):
                b_end = min(b_start + batch_size, nboot)
                batch_indices = list(range(b_start, b_end))
                
                if n_jobs_effective > 1 and len(batch_indices) > 1:
                    # 並列処理 - 各ワーカーにGPU使用フラグを渡す
                    batch_results_list = Parallel(n_jobs=n_jobs_effective, backend="loky")(
                        delayed(compute_permutation_batch)(
                            [j], 
                            ligand_info, 
                            receptor_info,
                            numCluster,
                            n,
                            k,
                            population_size,
                            cell_labels,
                            permutation,
                            valid_cell_types,
                            nC,
                            use_gpu_compute=use_gpu_for_batches and (j % 5 == 0)  # 5バッチごとにGPU使用
                        ) for j in batch_indices
                    )
                    # Combine results
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results_list[j_idx][:, 0]
                else:
                    # Single-threaded processing
                    batch_results = compute_permutation_batch(
                        batch_indices, 
                        ligand_info, 
                        receptor_info,
                        numCluster,
                        n,
                        k,
                        population_size,
                        cell_labels,
                        permutation,
                        valid_cell_types,
                        nC,
                        use_gpu_compute=use_gpu_for_batches
                    )
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results[:, j_idx]
                
                # GPUメモリ解放
                if has_gpu and use_gpu and not hybrid_mode and b_start % 100 == 0:
                    cp.get_default_memory_pool().free_all_blocks()
            
            # Calculate p-values - GPUでの効率的な計算
            if has_gpu and use_gpu and not hybrid_mode:
                try:
                    # GPUでp値計算
                    Pnull_vec_gpu = cp.array(Pnull_vec, dtype=gpu_dtype)
                    Pboot_gpu = cp.array(Pboot, dtype=gpu_dtype)
                    
                    # 閾値を超える値のカウント
                    nReject_gpu = cp.sum(Pboot_gpu > cp.expand_dims(Pnull_vec_gpu, 1), axis=1)
                    p_gpu = nReject_gpu / nboot
                    
                    # CPUに戻す
                    p = cp.asnumpy(p_gpu)
                    
                    # 検証
                    if validate_gpu_results:
                        nReject_cpu = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
                        p_cpu = nReject_cpu / nboot
                        
                        valid = validate_results(p_cpu, p, f"p-value calculation for LR pair {i}")
                        if not valid:
                            st.warning(f"Using CPU p-values for LR pair {i}")
                            p = p_cpu
                    
                    # メモリ解放
                    del Pnull_vec_gpu, Pboot_gpu, nReject_gpu, p_gpu
                    if i % 50 == 0:
                        cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    st.warning(f"GPU p-value calculation failed for LR pair {i}: {str(e)}. Using CPU.")
                    nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
                    p = nReject / nboot
            else:
                # CPU版
                nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
                p = nReject / nboot
            
            Pval[:, :, i] = p.reshape(numCluster, numCluster)
        
        show_progress(0.75, "Finalizing probability calculations...")
        
        # Set p-values to 1 where probability is 0
        Pval[Prob == 0] = 1
        
        # Apply p-value filtering if enabled - exact same approach
        if apply_pval_filter:
            Prob[Pval > trim_threshold] = 0
        
        # Set dimension names
        dimnames = [list(valid_cell_types), list(valid_cell_types), list(resource.index)]
        
        # Compute pathway-level communication probabilities
        show_progress(0.8, "Computing pathway-level communication...")
        netP = computeCommunProbPathway({"prob": Prob, "pval": Pval}, resource, 
                                       thresh=trim_threshold, apply_pval_filter=apply_pval_filter)
        
        # Compute centrality
        show_progress(0.85, "Computing network centrality...")
        
        # Compute centrality metrics for pathways
        if netP["pathways"] is not None and len(netP["pathways"]) > 0:
            netP["centr"] = {}
            for p_idx in range(len(netP["pathways"])):
                pathway_prob = netP["prob"][:, :, p_idx]
                pathway_prob_3d = np.expand_dims(pathway_prob, axis=2)
                netP["centr"][p_idx] = netAnalysis_computeCentrality(pathway_prob_3d)[0]
        
        # Compute centrality for aggregate network
        prob_sum = np.sum(Prob, axis=2)
        prob_sum_3d = np.expand_dims(prob_sum, axis=2)
        net_centr = netAnalysis_computeCentrality(prob_sum_3d)[0]
        
        # Compute aggregate network
        show_progress(0.9, "Computing aggregate network...")
        net_summary = aggregateCell_Cell_Communication({"prob": Prob, "pval": Pval}, 
                                                      valid_cell_types, pval_threshold=0.05, 
                                                      apply_pval_filter=apply_pval_filter)
        
        # Prepare results DataFrame
        show_progress(0.95, "Preparing results...")
        results_data = {
            'source': [],
            'target': [],
            'interaction_name': [],
            'ligand': [],
            'receptor': [],
            'prob': [],
            'pval': []
        }
        
        # Include only the significant interactions in results dataframe - exact same as original
        for i in range(len(valid_cell_types)):
            for j in range(len(valid_cell_types)):
                for k in range(nLR):
                    if Prob[i, j, k] > 0:
                        results_data['source'].append(valid_cell_types[i])
                        results_data['target'].append(valid_cell_types[j])
                        results_data['interaction_name'].append(resource.index[k])
                        results_data['ligand'].append(resource['ligand'].iloc[k])
                        results_data['receptor'].append(resource['receptor'].iloc[k])
                        results_data['prob'].append(Prob[i, j, k])
                        results_data['pval'].append(Pval[i, j, k])
        
        results_df = pd.DataFrame(results_data)
        
        # Clean up progress indicators if using streamlit
        if 'progress_bar' in locals() and progress_bar is not None:
            try:
                progress_bar.empty()
            except:
                pass
        if 'status_area' in locals() and status_area is not None:
            try:
                status_area.empty()
            except:
                pass
        
        # 最終的なメモリ解放
        if has_gpu and use_gpu:
            cp.get_default_memory_pool().free_all_blocks()
        
        # Calculate total execution time
        execution_time = time.time() - start_time
        logger.info(f"CellChat analysis completed in {execution_time:.2f} seconds")
        
        # Return the same structure as original cellchat_analysis
        return {
            'adata': adata_filtered,
            'results': results_df,
            'net': {
                "prob": Prob,
                "pval": Pval,
                "dimnames": dimnames,
                "centr": net_centr
            },
            'netP': netP,
            'network': net_summary,
            'groupby': groupby 
        }
    
    except Exception as e:
        logger.error(f"Error in GPU-optimized CellChat analysis: {str(e)}")
        logger.error(traceback.format_exc())
        
        # Clean up progress indicators if using streamlit
        if 'progress_bar' in locals() and progress_bar is not None:
            try:
                progress_bar.empty()
            except:
                pass
        if 'status_area' in locals() and status_area is not None:
            try:
                status_area.empty()
            except:
                pass
            
        # GPUメモリ解放
        if has_gpu and use_gpu:
            try:
                cp.get_default_memory_pool().free_all_blocks()
            except:
                pass
            
        return {'error': str(e), 'traceback': traceback.format_exc()}

def sanitize_filename(filename, max_length=20):
    """ファイル名を安全にして最大長さを制限する"""
    # 特殊文字を削除または置換
    filename = re.sub(r'[\\/*?:"<>|]', "_", filename)
    # 長さを制限
    if len(filename) > max_length:
        base, ext = os.path.splitext(filename)
        filename = base[:max_length] + ext
    return filename

def save_cellchat_result(result, uploaded_filename, selected_types, output_dir="temp"):
    """CellChatの結果をpickleファイルとして保存"""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # ファイル名を安全にして長さを制限
    safe_filename = sanitize_filename(uploaded_filename, 20)
    
    # シグナルタイプを短い形式に
    signal_types = "_".join([sig_type[:3] for sig_type in selected_types])
    
    # 保存ファイル名を作成
    save_filename = f"cellchat_{safe_filename}_{signal_types}.pkl"
    save_path = os.path.join(output_dir, save_filename)
    
    # 結果を保存
    with open(save_path, 'wb') as f:
        pickle.dump(result, f)
    
    return save_path

def load_cellchat_result(filepath):
    """保存したCellChat結果を読み込む"""
    with open(filepath, 'rb') as f:
        result = pickle.load(f)
    return result

def identify_overexpressed_genes(
    adata,
    group_by: str = None,
    idents_use: list = None,
    invert: bool = False,
    features_name: str = "features",
    only_pos: bool = True,
    features: list = None,
    return_object: bool = True,
    thresh_pct: float = 0,
    thresh_fc: float = 0,
    thresh_p: float = 0.05,
    do_de: bool = True,
    do_fast: bool = True,
    min_cells: int = 10,
):
    if not hasattr(pd.DataFrame, "iteritems"): #iteritemsがない場合の対応
        pd.DataFrame.iteritems = pd.DataFrame.items

    # 自動変換を有効化
    numpy2ri.activate()
    pandas2ri.activate()
    presto = rpackages.importr("presto")
    # データ行列の準備（ここでは genes x cells を前提とする）
    if scipy.sparse.issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = np.asarray(adata.X)
    # ※ AnnData では X は通常 cells×genes なので、必要に応じて転置してください
    # 例: X = X.T  # これで X の行が遺伝子、列が細胞となる

    # セルラベルの取得
    labels = adata.obs[group_by] if group_by else adata.obs.iloc[:, 0]
    labels = pd.Categorical(labels)
    
    # 遺伝子リストの取得
    features_use = list(adata.var_names) if features is None else list(set(features) & set(adata.var_names))
    
    # 過剰発現解析用の行列を作成（ここでは X の行は遺伝子、列は細胞と仮定）
   # data_use = X[np.array([adata.var_names.get_loc(gene) for gene in features_use]), :]
    X_gene = X.T 
    data_use = X_gene[np.array([adata.var_names.get_loc(gene) for gene in features_use]), :]

    
    # 各クラスタのラベル
    level_use = list(labels.categories)
    
    if do_de:
        if do_fast:
            # R の presto::wilcoxauc を呼び出す
            # data_use: 行が遺伝子、列が細胞のデータ行列
            # labels: 細胞ごとのグループ（長さは細胞数）
            with localconverter(ro.default_converter + pandas2ri.converter):
                # data_use_df として pandas.DataFrame に変換
                data_use_df = pd.DataFrame(data_use, index=features_use, columns=adata.obs_names)
                # R では行名は features、列名はサンプル名
            # labels を R の factor に変換
            r_labels = ro.FactorVector(labels.astype(str).tolist())
            # groups_use パラメータには、対象とするグループのベクトルを渡す（ここでは level_use をそのまま渡す）
            r_groups = ro.StrVector(level_use)
            
            # presto の wilcoxauc 呼び出し：
            # R の呼び出しは: wilcoxauc(x, group, groups_use = level_use)
            # x はデータ行列（features x cells）
            res = presto.wilcoxauc(data_use_df, r_labels, groups_use=r_groups)
            
            # res は R の DataFrame になるので、pandas に変換
            with localconverter(ro.default_converter + pandas2ri.converter):
                genes_de = pd.DataFrame(ro.conversion.rpy2py(res))
            
            # R版と同様に、列名の変更（例： "pval" -> "pvalues", "logFC" -> "logFC" など）
            # 以下は例として調整。実際は presto の返す列名に合わせて調整してください。
            genes_de.rename(columns={
                "group": "clusters",
                "feature": "features",
                "pval": "pvalues",
                "logFC": "logFC",
                "pct_in": "pct.1",
                "pct_out": "pct.2",
                "padj": "pvalues.adj"
            }, inplace=True)
            genes_de["logFC_abs"] = genes_de["logFC"].abs()
            # R版では pct.max は pct.1 と pct.2 の最大値（R ではパーセンテージとして扱う場合もあるため注意）
            genes_de["pct.max"] = genes_de[["pct.1", "pct.2"]].max(axis=1)
            
            # フィルタリング（ここでは R版同様 thresh_p, thresh_fc, thresh_pct を利用）
            # ※ R版では thresh_pct は 0.1 を 10% として比較している可能性があるため、ここは 0～1 のスケールと仮定
            markers_all = genes_de[(genes_de["pvalues"] < thresh_p) &
                                    (genes_de["logFC_abs"] >= thresh_fc) &
                                    (genes_de["pct.max"] > thresh_pct)]
            markers_all.sort_values("pvalues", inplace=True)
            
        else:
            # 従来の Python 実装（scipy の Mann–Whitney U 検定を使用）
            markers = []
            for cluster in level_use:
                cluster_mask = (labels == cluster)
                cluster_data = data_use[:, cluster_mask]
                other_data = data_use[:, ~cluster_mask]
                for i, gene in enumerate(features_use):
                    pct1 = np.mean(cluster_data[i, :] > 0)
                    pct2 = np.mean(other_data[i, :] > 0)
                    stat, p_val = scipy.stats.mannwhitneyu(cluster_data[i, :], other_data[i, :], alternative='two-sided')
                    log_fc = np.log(np.mean(cluster_data[i, :] + 1e-10) / np.mean(other_data[i, :] + 1e-10))
                    markers.append({
                        "clusters": cluster,
                        "features": gene,
                        "pvalues": p_val,
                        "logFC": log_fc,
                        "pct.1": pct1,
                        "pct.2": pct2
                    })
            markers_all = pd.DataFrame(markers)
            markers_all["logFC_abs"] = markers_all["logFC"].abs()
            markers_all["pct.max"] = markers_all[["pct.1", "pct.2"]].max(axis=1)
            markers_all = markers_all[(markers_all["pvalues"] < thresh_p) &
                                      (markers_all["logFC_abs"] >= thresh_fc) &
                                      (markers_all["pct.max"] > thresh_pct)]
            markers_all.sort_values("pvalues", inplace=True)
        
        features_sig = markers_all["features"].unique()
    else:
        # do.DE=False の場合：各遺伝子が min_cells 以上発現しているかで選別
        nCells = (data_use > 0).sum(axis=1)
        markers_all = pd.DataFrame({
            "features": features_use,
            "nCells": nCells
        })
        markers_all = markers_all[markers_all["nCells"] >= min_cells]
        features_sig = markers_all["features"].unique()
    
    # 結果を辞書として返す
  #  st.write(f"{len(features_sig)} genes passed filtering")
    return {"features": features_sig, "markers_all": markers_all}


@njit
def compute_hill_outer(dataL, dataR, k, n):
    numCluster = dataL.shape[0]
    result = np.empty((numCluster, numCluster), dtype=np.float64)
    kn = k ** n
    for i in prange(numCluster):
        for j in range(numCluster):
            x = dataL[i] * dataR[j]
            xn = x ** n
            result[i, j] = xn / (kn + xn)
    return result

# Python実装のCellChatDB処理をデバッグするコード
# get_cellchatdb_from_r関数の後に以下の関数を追加

def debug_cellchatdb(db):
    """CellChatDBの内容を詳細に出力して検証"""
    print("==== CellChatDB検証 ====")
    
    if 'interaction' in db:
        interaction = db['interaction']
        print(f"リガンド-レセプターペア数: {len(interaction)}")
        if len(interaction) > 0:
            print("最初の5ペア:")
            cols_to_show = ['ligand', 'receptor']
            if 'interaction_name' in interaction.columns:
                cols_to_show.insert(0, 'interaction_name')
            if 'pathway_name' in interaction.columns:
                cols_to_show.append('pathway_name')
            print(interaction[cols_to_show].head(5))
            
            # リガンド・レセプターの種類を確認
            ligands = interaction['ligand'].unique()
            receptors = interaction['receptor'].unique()
            
            print(f"\nユニークなリガンド数: {len(ligands)}")
            print(f"リガンドの例 (最初の10): {', '.join(map(str, ligands[:10]))}")
            
            print(f"\nユニークなレセプター数: {len(receptors)}")
            print(f"レセプターの例 (最初の10): {', '.join(map(str, receptors[:10]))}")
            
            # リガンド・レセプターの型を確認
            print(f"\nリガンドのデータ型: {type(interaction['ligand'].iloc[0])}")
            print(f"レセプターのデータ型: {type(interaction['receptor'].iloc[0])}")
            
            # 複合体のリガンド・レセプターの確認
            complex_ligands = [l for l in ligands if isinstance(l, str) and '_' in l]
            complex_receptors = [r for r in receptors if isinstance(r, str) and '_' in r]
            
            print(f"\n複合体リガンド数: {len(complex_ligands)}")
            if complex_ligands:
                print(f"複合体リガンドの例: {', '.join(complex_ligands[:5])}")
                
            print(f"複合体レセプター数: {len(complex_receptors)}")
            if complex_receptors:
                print(f"複合体レセプターの例: {', '.join(complex_receptors[:5])}")

            # データベース取得後、構造を詳細に確認
            print("DB構造:")
            for key in db:
                print(f"{key}: {type(db[key])}, 行数: {len(db[key])}")

            print("リガンドの最初の5つの例:")
            print(resource['ligand'].head(5).tolist())
            print("リガンドのデータ型:", resource['ligand'].dtype)
    else:
        print("CellChatDBにinteractionが含まれていません")
    
    # complex, cofactorの確認
    if 'complex' in db and not db['complex'].empty:
        print(f"\n複合体情報数: {len(db['complex'])}")
        print("複合体の例 (最初の3):")
        print(db['complex'].head(3))
    else:
        print("\n複合体情報が含まれていないか空です")
        
    if 'cofactor' in db and not db['cofactor'].empty:
        print(f"\n補助因子情報数: {len(db['cofactor'])}")
        print("補助因子の例 (最初の3):")
        print(db['cofactor'].head(3))
    else:
        print("\n補助因子情報が含まれていないか空です")
        
    print("==== CellChatDB検証完了 ====\n")
    
    return db  # 元のDBをそのまま返して元の処理を継続可能に

# cellchat_analysis関数内でのデバッグコードを追加（主要な計算ステップに挿入）

def debug_expression_data(adata, groupby):
    """発現データの検証"""
    print("==== 発現データ検証 ====")
    
    # 細胞タイプの確認
    cell_types = sorted(adata.obs[groupby].unique())
    print(f"細胞タイプ数: {len(cell_types)}")
    print(f"細胞タイプ: {', '.join(map(str, cell_types))}")
    
    # 遺伝子数の確認
    gene_count = len(adata.var_names)
    print(f"遺伝子数: {gene_count}")
    
    # 発現データの基本統計
    X = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    print("\n発現データの基本統計:")
    print(f"最小値: {np.min(X)}")
    print(f"最大値: {np.max(X)}")
    print(f"平均値: {np.mean(X)}")
    print(f"0より大きい値の割合: {np.mean(X > 0)}")
    temp_df = pd.DataFrame(
    adata.X[:5,:8].toarray() if scipy.sparse.issparse(adata.X) else adata.X[:5,:8],
    index=adata.obs_names[:5],
    columns=adata.var_names[:8]
    )
    st.write('adata.X')
    st.dataframe(temp_df) 
    
    # 一部の遺伝子の平均発現の確認
    print("\n例として最初の5遺伝子の平均発現:")
    selected_genes = adata.var_names[:5]
    
    for cell_type in cell_types[:3]:  # 最初の3つの細胞タイプのみ
        print(f"細胞タイプ: {cell_type}")
        cell_indices = np.where(adata.obs[groupby] == cell_type)[0]
        for i, gene in enumerate(selected_genes):
            gene_idx = adata.var_names.get_loc(gene)
            gene_expr = X[cell_indices, gene_idx]
            print(f"{gene} 平均発現: {np.mean(gene_expr)}")
            print(f"{gene} 発現率: {np.mean(gene_expr > 0)}")
        print()
    
    print("==== 発現データ検証完了 ====\n")

def debug_lr_calculations(resource, dataLavg, dataRavg, Prob, Pval, valid_cell_types, k=0.5, n=1):
    """リガンド-レセプター計算の詳細な検証"""
    print("==== LR計算検証 ====")
    
    # リガンド・レセプター発現値の統計
    print("リガンド発現値の統計:")
    print(f"  最小値: {np.min(dataLavg)}")
    print(f"  最大値: {np.max(dataLavg)}")
    print(f"  平均値: {np.mean(dataLavg)}")
    print(f"  0より大きい値の割合: {np.mean(dataLavg > 0)}")
    
    print("\nレセプター発現値の統計:")
    print(f"  最小値: {np.min(dataRavg)}")
    print(f"  最大値: {np.max(dataRavg)}")
    print(f"  平均値: {np.mean(dataRavg)}")
    print(f"  0より大きい値の割合: {np.mean(dataRavg > 0)}")
    
    # 具体的なLRペアの例
    print("\n具体的なLRペアの計算例:")
    for lr_idx in range(min(5, len(resource))):
        ligand = resource['ligand'].iloc[lr_idx]
        receptor = resource['receptor'].iloc[lr_idx]
        print(f"LRペア {lr_idx+1}: {ligand}-{receptor}")
        
        # 各細胞タイプでの発現
        for i, source in enumerate(valid_cell_types[:3]):  # 最初の3つの送信側
            for j, target in enumerate(valid_cell_types[:3]):  # 最初の3つの受信側
                # リガンド・レセプター発現値
                l_expr = dataLavg[lr_idx, i]
                r_expr = dataRavg[lr_idx, j]
                
                # リガンド・レセプター相互作用
                dataLR = l_expr * r_expr
                
                # ヒル関数による確率計算
                expected_prob = dataLR**n / (k**n + dataLR**n)
                actual_prob = Prob[i, j, lr_idx]
                p_value = Pval[i, j, lr_idx]
                
                if l_expr > 0 and r_expr > 0:
                    print(f"  {source} -> {target}:")
                    print(f"    リガンド発現: {l_expr:.6f}")
                    print(f"    レセプター発現: {r_expr:.6f}")
                    print(f"    L*R = {dataLR:.6f}")
                    print(f"    期待される確率: {expected_prob:.6f}")
                    print(f"    実際の確率: {actual_prob:.6f}")
                    print(f"    P値: {p_value:.6f}")
                    print()
    
    # 相互作用数のカウント
    total_interactions = np.sum(Prob > 0)
    print(f"総相互作用数 (確率 > 0): {total_interactions}")
    
    sig_interactions_05 = np.sum((Prob > 0) & (Pval <= 0.05))
    print(f"有意な相互作用数 (P ≤ 0.05): {sig_interactions_05}")
    
    sig_interactions_10 = np.sum((Prob > 0) & (Pval <= 0.1))
    print(f"有意な相互作用数 (P ≤ 0.1): {sig_interactions_10}")
    
    # 細胞タイプペアごとのインタラクション数
    interaction_counts = np.zeros((len(valid_cell_types), len(valid_cell_types)), dtype=int)
    for i in range(len(valid_cell_types)):
        for j in range(len(valid_cell_types)):
            interaction_counts[i, j] = np.sum((Prob[i, j, :] > 0) & (Pval[i, j, :] <= 0.05))
    
    print("\n細胞タイプペアごとのインタラクション数 (P ≤ 0.05):")
    interaction_counts_df = pd.DataFrame(interaction_counts, 
                                         index=valid_cell_types, 
                                         columns=valid_cell_types)
    print(interaction_counts_df)
    
    # TSVファイルとして保存
    interaction_counts_df.to_csv("python_interaction_number.tsv", sep="\t")
    
    # ヒル関数パラメータの確認
    print("\nヒル関数パラメータ確認")
    print(f"現在のk値: {k}")
    print(f"現在のn値: {n}")
    
    print("==== LR計算検証完了 ====\n")
    
# 計算過程の検証用関数
def verify_hill_function(k=0.5, n=1):
    """ヒル関数の計算を検証"""
    print("==== ヒル関数検証 ====")
    
    # サンプル値での計算
    test_values = [0, 0.001, 0.01, 0.1, 0.25, 0.5, 0.75, 1.0, 2.0, 5.0, 10.0]
    
    print(f"ヒル関数: y = x^n / (k^n + x^n) [k={k}, n={n}]")
    print("\n入力値 -> 出力値")
    for x in test_values:
        y = x**n / (k**n + x**n)
        print(f"{x:.4f} -> {y:.6f}")
    
    print("==== ヒル関数検証完了 ====\n")

@st.cache_data
def get_cellchatdb_from_r(species="human"):
    """
    RのCellChatパッケージからCellChatDBを取得する（rpy2 3.5.1に対応）
    
    Parameters
    ----------
    species : str
        'human' または 'mouse' を指定してデータベースを選択
        
    Returns
    -------
    dict
        'interaction', 'complex', 'cofactor', 'geneInfo' を含むCellChatDBの辞書（各値は Pandas DataFrame）
    """
    try:
        # 自動変換の有効化
        pandas2ri.activate()
        # baseパッケージは importr を使って読み込む
        base = importr("base")
        # CellChatパッケージのロード
        ro.r("library(CellChat)")
        
        # 対象のデータセットをロード
        if species.lower() == "human":
            r_command = """
            data("CellChatDB.human")
            CellChatDB.human <- CellChatDB.human
            """
            db_name = "CellChatDB.human"
        elif species.lower() == "mouse":
            r_command = """
            data("CellChatDB.mouse")
            CellChatDB.mouse <- CellChatDB.mouse
            """
            db_name = "CellChatDB.mouse"
        else:
            raise ValueError("species は 'human' または 'mouse' を指定してください")
        
        # Rコマンドの実行
        ro.r(r_command)
        
        db = {}
        components = ['interaction', 'complex', 'cofactor', 'geneInfo']
        for comp in components:
            # 各コンポーネントを取得（自動変換により Pandas DataFrame になっているはず）
            r_data = ro.r(f"{db_name}${comp}")
            # もし r_data が Pandas DataFrame ならそのまま使う
            if isinstance(r_data, pd.DataFrame):
                py_data = r_data
            else:
                # 変換を行う（通常はここには来ないはずです）
                with ro.default_converter + pandas2ri.converter as cv:
                    py_data = ro.conversion.rpy2py(r_data)
            db[comp] = py_data


        # get_cellchatdb_from_r関数内での修正
        # データフレーム取得後、型を明示的に変換
        if 'ligand' in db['interaction'].columns:
            db['interaction']['ligand'] = db['interaction']['ligand'].astype(str)
        if 'receptor' in db['interaction'].columns:
            db['interaction']['receptor'] = db['interaction']['receptor'].astype(str)



        return db
    except Exception as e:
        import traceback
        print("CellChatDBの取得中にエラーが発生しました: " + str(e))
        traceback.print_exc()
        return None


def debug_cellchatdb(db):
    """CellChatDBの内容を詳細に出力して検証"""
    print("==== CellChatDB検証 ====")
    
    if 'interaction' in db:
        interaction = db['interaction']
        print(f"リガンド-レセプターペア数: {len(interaction)}")
        if len(interaction) > 0:
            print("最初の5ペア:")
            cols_to_show = ['ligand', 'receptor']
            if 'interaction_name' in interaction.columns:
                cols_to_show.insert(0, 'interaction_name')
            if 'pathway_name' in interaction.columns:
                cols_to_show.append('pathway_name')
            print(interaction[cols_to_show].head(5))
            
            # リガンド・レセプターの種類を確認
            ligands = interaction['ligand'].unique()
            receptors = interaction['receptor'].unique()
            
            print(f"\nユニークなリガンド数: {len(ligands)}")
            print(f"リガンドの例 (最初の10): {', '.join(map(str, ligands[:10]))}")
            
            print(f"\nユニークなレセプター数: {len(receptors)}")
            print(f"レセプターの例 (最初の10): {', '.join(map(str, receptors[:10]))}")
            
            # リガンド・レセプターの型を確認
            print(f"\nリガンドのデータ型: {type(interaction['ligand'].iloc[0])}")
            print(f"レセプターのデータ型: {type(interaction['receptor'].iloc[0])}")
            
            # 複合体のリガンド・レセプターの確認
            complex_ligands = [l for l in ligands if isinstance(l, str) and '_' in l]
            complex_receptors = [r for r in receptors if isinstance(r, str) and '_' in r]
            
            print(f"\n複合体リガンド数: {len(complex_ligands)}")
            if complex_ligands:
                print(f"複合体リガンドの例: {', '.join(complex_ligands[:5])}")
                
            print(f"複合体レセプター数: {len(complex_receptors)}")
            if complex_receptors:
                print(f"複合体レセプターの例: {', '.join(complex_receptors[:5])}")
    else:
        print("CellChatDBにinteractionが含まれていません")
    
    # complex, cofactorの確認
    if 'complex' in db and not db['complex'].empty:
        print(f"\n複合体情報数: {len(db['complex'])}")
        print("複合体の例 (最初の3):")
        print(db['complex'].head(3))
    else:
        print("\n複合体情報が含まれていないか空です")
        
    if 'cofactor' in db and not db['cofactor'].empty:
        print(f"\n補助因子情報数: {len(db['cofactor'])}")
        print("補助因子の例 (最初の3):")
        print(db['cofactor'].head(3))
    else:
        print("\n補助因子情報が含まれていないか空です")
        
    print("==== CellChatDB検証完了 ====\n")
    
    return db  # 元のDBをそのまま返して元の処理を継続可能に

# ディレクトリ処理用ヘルパー関数
@st.cache_data
def clear_old_directories(path):
    import shutil
    import os
    from datetime import datetime, timedelta
    
    now = datetime.now()
    for dir_name in os.listdir(path):
        dir_path = os.path.join(path, dir_name)
        if os.path.isdir(dir_path):
            try:
                # ディレクトリ名からタイムスタンプを取得
                timestamp = float(dir_name)
                dir_time = datetime.fromtimestamp(timestamp)
                # 24時間以上前のディレクトリを削除
                if now - dir_time > timedelta(hours=24):
                    shutil.rmtree(dir_path)
            except:
                pass

@st.cache_data
def clear_old_files(path):
    import os
    from datetime import datetime, timedelta
    
    now = datetime.now()
    for file_name in os.listdir(path):
        file_path = os.path.join(path, file_name)
        if os.path.isfile(file_path):
            file_time = datetime.fromtimestamp(os.path.getmtime(file_path))
            # 24時間以上前のファイルを削除
            if now - file_time > timedelta(hours=24):
                os.remove(file_path)

@st.cache_data
def check_species_index(gene_list):
    """遺伝子シンボルの大文字出現頻度からヒトかマウスかを推測する"""
    if not gene_list:
        return 0  # 空のリストの場合はデフォルトでヒト
    
    # サンプリング（遺伝子リストが大きい場合）
    sample_genes = gene_list[:500] if len(gene_list) > 500 else gene_list
    
    # 大文字の出現頻度を計算
    uppercase_ratios = []
    for gene in sample_genes:
        if not gene or not isinstance(gene, str):
            continue
        uppercase_count = sum(1 for char in gene if char.isupper())
        ratio = uppercase_count / len(gene) if len(gene) > 0 else 0
        uppercase_ratios.append(ratio)
    
    # 平均値を計算
    avg_uppercase_ratio = sum(uppercase_ratios) / len(uppercase_ratios) if uppercase_ratios else 0
    
    # 結果をログに出力
    print(f"遺伝子シンボルの大文字出現頻度: {avg_uppercase_ratio:.2f}")
    
    # ヒト：大文字が多い（BRCA1）、マウス：小文字が多い（Brca1）
    return 1 if avg_uppercase_ratio > 0.5 else 0  # 0.5未満ならマウス、それ以上ならヒト

def find_first_index_or_zero(lst, elements):
    for element in elements:
        try:
            return lst.index(element)
        except ValueError:
            continue
    return 0

@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    return adata

# 最適化された行列処理関数
def optimize_matrix_operations(X, adata):
    """
    行列操作を最適化する関数
    """
    logger.info(f"X形状: {X.shape}, adata.obs_names長さ: {len(adata.obs_names)}, adata.var_names長さ: {len(adata.var_names)}")
    
    # スパース行列の変換（一度に変換してメモリ効率を向上）
    if scipy.sparse.issparse(X):
        logger.info("スパース行列を密行列に変換しています")
        
        # メモリ効率のために部分的に変換
        chunk_size = 5000  # 適切なチャンクサイズを選択
        
        if X.shape[0] > chunk_size:
            logger.info(f"大きな行列を{chunk_size}行ずつ変換します")
            result = []
            for i in range(0, X.shape[0], chunk_size):
                end = min(i + chunk_size, X.shape[0])
                result.append(X[i:end].toarray())
            X = np.vstack(result)
        else:
            X = X.toarray()
    
    return X


# 集計関数の最適化
def calculate_mean_expression_optimized(data_use, cell_labels, cell_types, min_cells, FunMean):
    """
    細胞タイプごとの平均発現量を効率的に計算
    """
    data_use_avg_dict = {}
    cell_counts = {}
    
    # 前処理：各細胞タイプのインデックスをまとめて取得
    cell_type_indices = {}
    for cell_type in cell_types:
        indices = np.where(cell_labels == cell_type)[0]
        cell_counts[cell_type] = len(indices)
        cell_type_indices[cell_type] = indices
    
    # 有効な細胞タイプだけを処理
    valid_cell_types = [ct for ct in cell_types if cell_counts.get(ct, 0) >= min_cells]
    
    # 各有効な細胞タイプについて平均発現量を計算
    for cell_type in valid_cell_types:
        indices = cell_type_indices[cell_type]
        # メモリ効率のために一度に処理する遺伝子数を制限
        gene_chunk_size = 1000
        n_genes = data_use.shape[1]
        avg_expr = np.zeros(n_genes)
        
        for i in range(0, n_genes, gene_chunk_size):
            end = min(i + gene_chunk_size, n_genes)
            data_subset = data_use[indices, i:end]
            avg_expr[i:end] = np.apply_along_axis(FunMean, 0, data_subset)
        
        data_use_avg_dict[cell_type] = avg_expr
    
    return data_use_avg_dict, cell_counts, valid_cell_types

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
import logging

logger = logging.getLogger("CellChat")

def apply_mean_function(data_subset, fun_type='triMean', trim=0.1):
    """平均計算関数（Numbaなし）"""
    if fun_type == 'triMean':
        # 四分位数の計算と平均（Numba未使用）
        q1 = np.quantile(data_subset, 0.25, axis=0)
        q2 = np.quantile(data_subset, 0.5, axis=0)  # median
        q3 = np.quantile(data_subset, 0.75, axis=0)
        return (q1 + 2*q2 + q3) / 4
    elif fun_type == 'truncatedMean':
        # トリム平均
        return np.apply_along_axis(lambda x: np.mean(x, axis=0), 0, data_subset)
    else:
        # デフォルトは単純平均
        return np.mean(data_subset, axis=0)

def process_single_permutation(data_use, cluster_indices, valid_cell_types, fun_type='triMean', trim=0.1):
    """単一permutationの全クラスタと全遺伝子の平均発現を計算（Numbaなし）"""
    n_genes = data_use.shape[1]
    numCluster = len(valid_cell_types)
    result = np.zeros((n_genes, numCluster), dtype=np.float32)
    
    # 各クラスターについて処理
    for ct_idx in range(numCluster):
        cells = cluster_indices[ct_idx]
        if len(cells) > 0:
            data_subset = data_use[cells]
            result[:, ct_idx] = apply_mean_function(data_subset, fun_type, trim)
    
    return result

def process_permutation_batch(batch_indices, data_use, cell_labels, permutation, valid_cell_types, 
                             fun_type='triMean', trim=0.1):
    """permutationのバッチ処理関数"""
    n_genes = data_use.shape[1]
    numCluster = len(valid_cell_types)
    results = np.zeros((n_genes, numCluster, len(batch_indices)), dtype=np.float32)
    
    for idx, j in enumerate(batch_indices):
        # j番目のPermutation後の細胞タイプ
        group_boot = cell_labels.values[permutation[:, j]]
        
        # 各クラスターに属する細胞のインデックスを取得
        cluster_indices = [np.where(group_boot == ct)[0] for ct in valid_cell_types]
        
        # すべての遺伝子と細胞タイプの平均発現を計算
        results[:, :, idx] = process_single_permutation(
            data_use, cluster_indices, valid_cell_types, fun_type, trim
        )
    
    return results

def precompute_gene_expressions(data_use, cell_labels, permutation, valid_cell_types, FunMean, nboot, n_jobs=32):
    """
    事前に全遺伝子×全クラスター×全permutationの平均発現をまとめて計算
    （Numbaなしのバージョン）
    
    Parameters
    ----------
    data_use : np.ndarray
        正規化された発現行列 (cells × genes)
    cell_labels : pd.Series
        細胞タイプラベル
    permutation : np.ndarray
        置換テスト用のインデックス配列 (cells × nboot)
    valid_cell_types : list
        有効な細胞タイプのリスト
    FunMean : function
        平均計算関数
    nboot : int
        置換テスト回数
    n_jobs : int
        並列処理で使用するコア数
        
    Returns
    -------
    np.ndarray
        事前計算された平均発現配列 (genes × clusters × nboot)
    list
        単純化された置換クラス平均アクセス用のリスト
    """
    logger.info("全遺伝子×全クラスター×全置換の平均発現を事前計算中...")
    
    # FunMeanから計算タイプを決定（この部分は簡略化）
    fun_type = 'triMean'  # デフォルト
    trim = 0.1  # デフォルト値
    
    numCluster = len(valid_cell_types)
    n_genes = data_use.shape[1]
    
    # メモリ使用量を制御するためのバッチ処理
    batch_size = min(50, nboot)  # 一度に処理するpermutationの数を調整
    
    # 全遺伝子×全クラスター×全permutationの平均発現を格納する配列
    all_gene_expr = np.zeros((n_genes, numCluster, nboot), dtype=np.float32)
    
    # バッチに分割して並列処理
    all_batches = [list(range(i, min(i+batch_size, nboot))) for i in range(0, nboot, batch_size)]
    
    # コア数を制限
    n_jobs_effective = min(n_jobs, len(all_batches), 32)  # 最大8コアに制限
    
    if n_jobs_effective > 1:
        # 並列処理で計算
        logger.info(f"{n_jobs_effective}コアで並列処理を実行中...")
        results = Parallel(n_jobs=n_jobs_effective)(
            delayed(process_permutation_batch)(
                batch_indices, data_use, cell_labels, permutation, valid_cell_types, fun_type, trim
            ) for batch_indices in all_batches
        )
        
        # 結果の統合
        for batch_idx, batch_indices in enumerate(all_batches):
            for idx, j in enumerate(batch_indices):
                if j < nboot:  # 範囲チェック
                    all_gene_expr[:, :, j] = results[batch_idx][:, :, idx]
    else:
        # シングルコアで処理
        logger.info("シングルコアで処理を実行中...")
        for j in range(nboot):
            if j % 10 == 0:
                logger.info(f"Permutation {j+1}/{nboot} を処理中...")
            
            # j番目のPermutation後の細胞タイプ
            group_boot = cell_labels.values[permutation[:, j]]
            
            # 各クラスターに属する細胞のインデックスを取得
            cluster_indices = [np.where(group_boot == ct)[0] for ct in valid_cell_types]
            
            # すべての遺伝子と細胞タイプの平均発現を計算
            all_gene_expr[:, :, j] = process_single_permutation(
                data_use, cluster_indices, valid_cell_types, fun_type, trim
            )
    
    logger.info("All mean gene expression is calclulated.")
    
    return all_gene_expr


def precompute_complex_mapping(complex_input, gene_to_index):
    """
    複合体とその構成遺伝子のインデックスマッピングを事前計算
    """
    complex_mapping = {}
    
    if complex_input.empty:
        return complex_mapping
    
    for complex_name in complex_input.index:
        # 複合体のサブユニットを取得
        subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
        subunits = complex_input.loc[complex_name, subunits_cols].dropna().astype(str)
        subunits = [s for s in subunits if s != "" and s in gene_to_index]
        
        if subunits:
            # サブユニットの遺伝子インデックスを保存
            complex_mapping[complex_name] = [gene_to_index[s] for s in subunits]
    
    return complex_mapping


def geometricMean(expr_values):
    """
    Rの実装に合わせた幾何平均計算
    """
    # 値がゼロの場合の対処（極小値に置換）
    epsilon = 1e-10
    expr_values = np.maximum(expr_values, epsilon)
    
    if expr_values.ndim == 1:
        # 1次元配列の場合
        return np.exp(np.mean(np.log(expr_values)))
    else:
        # 2次元配列の場合（列ごとに計算）
        return np.exp(np.mean(np.log(expr_values), axis=0))

@st.cache_data
def computeExpr_coreceptor(cofactor_input, data_use, pairLRsig, type_coreceptor):
    """
    リガンド-レセプターの相互作用における共受容体効果をモデル化
    """
    if cofactor_input.empty or pairLRsig.empty:
        return np.ones((len(pairLRsig), data_use.shape[1]))
    
    coreceptor_col = 'co_A_receptor' if type_coreceptor == "A" else 'co_I_receptor'
    
    if coreceptor_col not in pairLRsig.columns:
        return np.ones((len(pairLRsig), data_use.shape[1]))
    
    coreceptor_all = pairLRsig[coreceptor_col].values
    index_coreceptor = np.where((pd.notna(coreceptor_all)) & (coreceptor_all != ""))[0]
    numCluster = data_use.shape[1]
    data_coreceptor = np.ones((len(coreceptor_all), numCluster))
    
    if len(index_coreceptor) > 0:
        for idx in index_coreceptor:
            coreceptor = coreceptor_all[idx]
            if coreceptor in cofactor_input.index:
                # 補助因子を取得
                cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
                cofactors = cofactor_input.loc[coreceptor, cofactor_cols].dropna().astype(str)
                # cofactors = [c for c in cofactors if c != "" and c in data_use.index]
                cofactors = [str(c) for c in cofactors if c != "" and str(c) in data_use.index]
                
                if len(cofactors) == 1:
                    data_coreceptor[idx] = 1 + data_use.loc[cofactors[0]].values
                elif len(cofactors) > 1:
                    # 1 + 各共受容体の発現の積
                    prod = np.ones(numCluster)
                    for c in cofactors:
                        prod *= (1 + data_use.loc[c].values)
                    data_coreceptor[idx] = prod
    
    return data_coreceptor

@st.cache_data
def computeExpr_agonist(data_use, pairLRsig, cofactor_input, index_agonist, Kh, n):
    """
    アゴニストがリガンド-レセプター相互作用に与える効果をモデル化
    """
    if cofactor_input.empty or pairLRsig.empty or 'agonist' not in pairLRsig.columns:
        return np.ones(data_use.shape[1])
    
    agonist = pairLRsig['agonist'].iloc[index_agonist]
    if pd.isna(agonist) or agonist == "" or agonist not in cofactor_input.index:
        return np.ones(data_use.shape[1])
    
    # アゴニスト遺伝子を取得
    cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
    agonist_genes = cofactor_input.loc[agonist, cofactor_cols].dropna().astype(str)
    agonist_genes = [g for g in agonist_genes if g != "" and g in data_use.index]
    
    if len(agonist_genes) == 1:
        data_avg = data_use.loc[agonist_genes[0]].values
        data_agonist = 1 + data_avg**n / (Kh**n + data_avg**n)
    elif len(agonist_genes) > 1:
        # 各アゴニスト遺伝子の効果の積
        data_agonist = np.ones(data_use.shape[1])
        for g in agonist_genes:
            data_avg = data_use.loc[g].values
            data_agonist *= (1 + data_avg**n / (Kh**n + data_avg**n))
    else:
        data_agonist = np.ones(data_use.shape[1])
    
    return data_agonist

@st.cache_data
def computeExpr_antagonist(data_use, pairLRsig, cofactor_input, index_antagonist, Kh, n):
    """
    アンタゴニストがリガンド-レセプター相互作用に与える効果をモデル化
    """
    if cofactor_input.empty or pairLRsig.empty or 'antagonist' not in pairLRsig.columns:
        return np.ones(data_use.shape[1])
    
    antagonist = pairLRsig['antagonist'].iloc[index_antagonist]
    if pd.isna(antagonist) or antagonist == "" or antagonist not in cofactor_input.index:
        return np.ones(data_use.shape[1])
    
    # アンタゴニスト遺伝子を取得
    cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
    antagonist_genes = cofactor_input.loc[antagonist, cofactor_cols].dropna().astype(str)
    antagonist_genes = [g for g in antagonist_genes if g != "" and g in data_use.index]
    
    if len(antagonist_genes) == 1:
        data_avg = data_use.loc[antagonist_genes[0]].values
        data_antagonist = Kh**n / (Kh**n + data_avg**n)
    elif len(antagonist_genes) > 1:
        # 各アンタゴニスト遺伝子の効果の積
        data_antagonist = np.ones(data_use.shape[1])
        for g in antagonist_genes:
            data_avg = data_use.loc[g].values
            data_antagonist *= Kh**n / (Kh**n + data_avg**n)
    else:
        data_antagonist = np.ones(data_use.shape[1])
    
    return data_antagonist

@st.cache_data
def computeCommunProbPathway(net, pairLR_use, thresh=0.05, apply_pval_filter=True):
    """
    シグナル経路レベルでの通信確率を計算
    
    Parameters
    ----------
    net : dict
        通信確率と有意性を含む辞書
    pairLR_use : pd.DataFrame
        リガンド-レセプターペア情報
    thresh : float, optional
        有意性判定のp値閾値
    apply_pval_filter : bool, optional
        p値フィルタリングを適用するかどうか。デフォルトはTrue
        
    Returns
    -------
    dict
        パスウェイレベルでの通信確率情報
    """
    prob = net["prob"].copy()
    pval = net["pval"]
    
    # apply_pval_filterパラメータに基づいてp値でのフィルタリングを適用
    if apply_pval_filter:
        prob[pval > thresh] = 0
    
    # パスウェイ情報がない場合は処理をスキップ
    if 'pathway_name' not in pairLR_use.columns:
        return {
            "pathways": [],
            "prob": np.zeros((prob.shape[0], prob.shape[1], 0))
        }
    
    # シグナル経路ごとに集計
    pathways = pairLR_use['pathway_name'].dropna().unique()
    prob_pathways = np.zeros((prob.shape[0], prob.shape[1], len(pathways)))
    
    for i, pathway in enumerate(pathways):
        idx = np.where(pairLR_use['pathway_name'] == pathway)[0]
        prob_pathways[:, :, i] = np.sum(prob[:, :, idx], axis=2)
    
    # 総通信強度に基づいて並べ替え
    pathway_sums = np.sum(prob_pathways, axis=(0, 1))
    significant_idx = np.where(pathway_sums > 0)[0]
    
    if len(significant_idx) == 0:
        return {
            "pathways": [],
            "prob": np.zeros((prob.shape[0], prob.shape[1], 0))
        }
    
    sort_idx = significant_idx[np.argsort(-pathway_sums[significant_idx])]
    
    pathways_sig = pathways[sort_idx]
    prob_pathways_sig = prob_pathways[:, :, sort_idx]
    
    return {
        "pathways": pathways_sig,
        "prob": prob_pathways_sig
    }


@st.cache_data
def computeExpr_complex(complex_input, data_use, complex_genes):
    """複合体の発現を計算する関数（Rの実装に対応）"""
    result = np.zeros((len(complex_genes), data_use.shape[1]))
    
    for i, complex_gene in enumerate(complex_genes):
        if complex_gene in complex_input.index:
            # 複合体のサブユニットを取得
            subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
            subunits = complex_input.loc[complex_gene, subunits_cols].dropna().astype(str)
            subunits = [s for s in subunits if s != "" and s in data_use.index]
            
            if len(subunits) > 0:
                # 発現値の取得
                expr_values = data_use.loc[subunits].values
                # 幾何平均の計算
                result[i] = geometricMean(expr_values)
    
    return result

@st.cache_data
def computeExpr_LR(geneLR, data_use, complex_input):
    """リガンドまたはレセプターの発現を計算する関数（Rの実装に対応）"""
    geneLR = [str(gene) for gene in geneLR]  # 明示的に文字列に変換
    print(f"First 5 geneLR values: {geneLR[:5]}")
    print(f"Data_use.index type: {type(data_use.index)}")
    print(f"First 5 gene names in data_use: {list(data_use.index)[:5]}")
    print(f"Gene name types in data_use: {[type(g) for g in data_use.index[:5]]}")

    nLR = len(geneLR)
    numCluster = data_use.shape[1]
    
    # 単一遺伝子の処理（Rと同様の方法）
    index_singleL = [i for i, gene in enumerate(geneLR) if gene in data_use.index]
    dataLavg = np.zeros((nLR, numCluster))
    
    if index_singleL:
        # 単一遺伝子の発現をまとめて取得
        gene_indices = [geneLR[i] for i in index_singleL if geneLR[i] in data_use.index]
        if gene_indices:
            dataL1avg = data_use.loc[gene_indices].values
            for idx, gene_idx in enumerate(index_singleL):
                if idx < len(dataL1avg):
                    dataLavg[gene_idx] = dataL1avg[idx]
    
    # 複合体の処理（Rと同様の方法）
    index_complexL = [i for i in range(nLR) if i not in index_singleL]
    if index_complexL and not complex_input.empty:
        complex_genes = [geneLR[i] for i in index_complexL]
        data_complex = computeExpr_complex(complex_input, data_use, complex_genes)
        for idx, complex_idx in enumerate(index_complexL):
            dataLavg[complex_idx] = data_complex[idx]
    
    return dataLavg





@st.cache_data
def netAnalysis_computeCentrality(prob):
    """
    ネットワーク中心性指標を計算 (改良版)
    
    Parameters
    ----------
    prob : numpy.ndarray
        通信確率行列 (形状: [cell_types, cell_types, interactions])
    debug_mode : bool, optional
        デバッグモード
        
    Returns
    -------
    centrality : dict
        各インタラクション（ネットワーク）ごとの中心性指標のディクショナリ
    """
    import numpy as np
    import traceback
    
    # 入力が2次元の場合は3次元に拡張（軸エラー回避）
    if prob.ndim == 2:
        prob = np.expand_dims(prob, axis=2)
    
    # rpy2の設定と必要なRパッケージのインポート
    r_available = False
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import numpy2ri
        
        # 自動変換の有効化
        numpy2ri.activate()
        
        # Rライブラリのインポート
        base = importr('base')
        sna = importr('sna')
        
        # igraphパッケージの読み込み
        try:
            ro.r('library(igraph, quietly=TRUE)')
            r_available = True
            print("Rライブラリの読み込み成功: igraph, sna")
        except Exception as e:
            print(f"igraph読み込みエラー: {str(e)}")
            r_available = False
            
    except Exception as e:
        print(f"Rライブラリの読み込みエラー: {str(e)}")
        print(traceback.format_exc())
        print("NetworkXでフォールバックします。")
    
    # NetworkXを使用した計算
  #  centrality_nx = netAnalysis_computeCentrality_nx(prob)
    

    try:
        # Rの関数を定義
        ro.r('''
        computeCentralityLocal <- function(net) {
          centr <- list()
          G <- igraph::graph_from_adjacency_matrix(net, mode = "directed", weighted = TRUE)
          centr$outdeg_unweighted <- rowSums(net > 0)
          centr$indeg_unweighted <- colSums(net > 0)
          centr$outdeg <- igraph::strength(G, mode="out")
          centr$indeg <- igraph::strength(G, mode="in")
          centr$hub <- igraph::hub_score(G)$vector
          centr$authority <- igraph::authority_score(G)$vector
          centr$eigen <- igraph::eigen_centrality(G)$vector
          centr$page_rank <- igraph::page_rank(G)$vector
          igraph::E(G)$weight <- 1/igraph::E(G)$weight
          centr$betweenness <- igraph::betweenness(G)
          centr$flowbet <- tryCatch({
            sna::flowbet(net)
          }, error = function(e) {
            rep(0, nrow(net))
          })
          centr$info <- tryCatch({
            sna::infocent(net, diag = TRUE, rescale = TRUE, cmode = "lower")
          }, error = function(e) {
            rep(0, nrow(net))
          })
          return(centr)
        }
        ''')

        centrality_r = {}
        for i in range(prob.shape[2]):
            try:
                net_mat = prob[:, :, i]
                
                # NumPy行列をRの行列に変換
                r_mat = ro.r.matrix(net_mat, nrow=net_mat.shape[0], ncol=net_mat.shape[1])
                
                # R関数を実行
                r_centr = ro.r['computeCentralityLocal'](r_mat)
                
                # Rの結果をPythonの辞書に変換
                centrality_r[i] = {}
                for metric in ["outdeg", "indeg", "outdeg_unweighted", "indeg_unweighted", 
                              "hub", "authority", "eigen", "page_rank", "betweenness", "flowbet", "info"]:
                    try:
                        if metric in r_centr.names:
                            centrality_r[i][metric] = np.array(r_centr.rx2(metric))
                    except Exception as e:
                        print(f"指標 {metric} の取得中にエラー: {str(e)}")
                        centrality_r[i][metric] = np.zeros(net_mat.shape[0])  # デフォルト値を設定
            
            except Exception as e:
                print(f"R計算でエラー (interaction {i}): {str(e)}")
                centrality_r[i] = centrality_nx[i]  # NetworkXの結果を使用
        
        # デバッグモードの場合、結果を比較
  #      if debug_mode:
  #          compare_centrality_results(centrality_r, centrality_nx)
        
        return centrality_r
    except Exception as e:
        print(f"Rでの中心性計算に失敗しました: {str(e)}")
        print(traceback.format_exc())
        print("NetworkXの結果を使用します。")

@st.cache_data
def aggregateCell_Cell_Communication(net, valid_cell_types, pval_threshold=0.05, apply_pval_filter=True):
    """
    細胞間通信ネットワークを集計
    
    Parameters
    ----------
    net : dict
        通信確率と有意性を含む辞書
    valid_cell_types : list
        有効な細胞タイプのリスト
    pval_threshold : float, optional
        有意性判定のp値閾値。デフォルトは0.05
    apply_pval_filter : bool, optional
        p値フィルタリングを適用するかどうか。デフォルトはTrue
        
    Returns
    -------
    dict
        集計されたネットワーク情報
    """
    prob = net["prob"]
    pval = net["pval"]
    
    # apply_pval_filterパラメータに基づいてp値でのフィルタリングを適用
    sig_prob = prob.copy()
    if apply_pval_filter:
        sig_prob[pval > pval_threshold] = 0
    
    # 各細胞タイプペア間の総インタラクション強度
    strength_matrix = np.sum(sig_prob, axis=2)
    
    # インタラクション数
    count_matrix = np.sum(sig_prob > 0, axis=2)
    
    # 各リガンド-レセプターペアの寄与度
    lr_contribution = np.sum(np.sum(sig_prob, axis=0), axis=0)
    
    # 各細胞タイプの送受信総量
    outgoing = np.sum(strength_matrix, axis=1)
    incoming = np.sum(strength_matrix, axis=0)
    
    # 行列をDataFrameに変換
    strength_df = pd.DataFrame(strength_matrix, index=valid_cell_types, columns=valid_cell_types)
    count_df = pd.DataFrame(count_matrix, index=valid_cell_types, columns=valid_cell_types)
    
    # ネットワーク指標
    try:
        network_centrality = calculate_network_centrality({"strength_matrix": strength_df})
    except Exception as e:
        print(f"ネットワーク中心性計算エラー: {str(e)}")
        import traceback
        print(traceback.format_exc())
        network_centrality = pd.DataFrame()

    
    # 行列の統計情報
    strength_stats = {
        'min': float(np.min(strength_matrix)),
        'max': float(np.max(strength_matrix)),
        'mean': float(np.mean(strength_matrix)),
        'std': float(np.std(strength_matrix)),
        'non_zero': int(np.sum(strength_matrix > 0))
    }
    
    return {
        'strength_matrix': strength_df,
        'count_matrix': count_df,
        'lr_contribution': lr_contribution,
        'outgoing': pd.Series(outgoing, index=valid_cell_types),
        'incoming': pd.Series(incoming, index=valid_cell_types),
        'network_centrality': network_centrality,
        'strength_stats': strength_stats
    }





def plot_interaction_heatmap(matrix, order=None, title="Cell-Cell Interaction Strength", figsize=(10, 8)):
    """
    細胞間相互作用のヒートマップ描画関数
    
    Parameters
    ----------
    matrix : pd.DataFrame or numpy.ndarray
        相互作用強度/カウント行列
    title : str
        プロットのタイトル
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:

        if order is not None and isinstance(matrix, pd.DataFrame):
            matrix = matrix.reindex(order, axis=0).reindex(order, axis=1)

        fig, ax = plt.subplots(figsize=figsize)
        
        # Handle different matrix types
        is_pandas = hasattr(matrix, 'empty')
        
        # Check if matrix is empty
        if (is_pandas and matrix.empty) or (not is_pandas and matrix.size == 0):
            ax.text(0.5, 0.5, "データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title(title)
            return fig
            
        # Get values appropriately based on matrix type
        values = matrix.values.flatten() if is_pandas else matrix.flatten()
        all_zeros = np.allclose(values, 0, rtol=1e-5, atol=1e-8)
        all_same = np.allclose(values, values[0], rtol=1e-5, atol=1e-8) if len(values) > 0 else True
        
        if len(matrix) == 0:
            ax.text(0.5, 0.5, "データがありません", ha='center', va='center')
            ax.axis('off')
        elif all_zeros:
            # すべてゼロの場合
            # For numpy arrays, convert to DataFrame for heatmap
            if not is_pandas:
                matrix = pd.DataFrame(matrix)
            
            sns.heatmap(
                matrix, 
                cmap="YlOrRd", 
                square=True, 
                linewidths=0.5,
                annot=True,
                fmt=".6f", 
                ax=ax
            )
            plt.title(title)
            ax.text(0.5, -0.1, "すべての値がゼロです。パラメータを調整してください。", 
                    ha='center', va='center', transform=ax.transAxes, 
                    fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
        elif all_same and not all_zeros:
            # すべて同じ非ゼロ値
            if not is_pandas:
                matrix = pd.DataFrame(matrix)
                
            sns.heatmap(
                matrix, 
                cmap="YlOrRd", 
                square=True, 
                linewidths=0.5,
                annot=True,
                fmt=".6f", 
                ax=ax
            )
            plt.title(title)
            ax.text(0.5, -0.1, "すべての値が同じです。パラメータを調整してください。", 
                    ha='center', va='center', transform=ax.transAxes, 
                    fontsize=12, bbox=dict(facecolor='white', alpha=0.8))
        else:
            # 通常のヒートマップ表示
            # 値の範囲を確認
            vmin = np.min(values)
            vmax = np.max(values)
            value_range = vmax - vmin
            
            # 値の範囲が非常に小さい場合
            if value_range < 1e-6:
                logger.warning(f"値の範囲が非常に小さいです: {value_range:.2e}")
                
            # 適切なカラーマップを選択
            if value_range < 0.01:
                logger.info("値の範囲が小さいため、カラーマップを調整します")
                cmap = "YlOrRd"
                vmin = None
                vmax = None
            else:
                cmap = "YlOrRd"
                vmin = None
                vmax = None
            
            # Convert numpy array to DataFrame for seaborn
            if not is_pandas:
                matrix = pd.DataFrame(matrix)
            
            sns.heatmap(
                matrix, 
                cmap=cmap, 
                square=True, 
                linewidths=0.5,
                annot=True,  # 値を表示
                fmt=".6f",   # 6桁で表示
                ax=ax,
                vmin=vmin,
                vmax=vmax
            )
            
            plt.title(title)
        
        plt.tight_layout()
        return fig
    except Exception as e:
        logger.error(f"ヒートマップ描画エラー: {str(e)}")
        logger.error(traceback.format_exc())
        # エラーメッセージを含む空のプロットを返す
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"ヒートマップ作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title(title)
        return fig

def debug_interaction_matrix(matrix, title="Interaction Matrix Debug"):
    """
    相互作用行列のデバッグ情報を出力
    
    Parameters
    ----------
    matrix : pd.DataFrame or numpy.ndarray
        相互作用行列
    title : str
        タイトル
    """
    print(f"===== {title} =====")
    
    # Check if matrix is empty - handle both DataFrame and ndarray
    if hasattr(matrix, 'empty'):
        # For pandas DataFrame
        if matrix.empty:
            print("空の行列です")
            return
        values = matrix.values
        print(f"Shape: {matrix.shape}")
        print(f"Min value: {np.min(values)}")
        print(f"Max value: {np.max(values)}")
        print(f"Mean value: {np.mean(values)}")
        print(f"Std value: {np.std(values)}")
        print(f"Non-zero entries: {np.count_nonzero(values)} / {matrix.size}")
        print(f"Sample values:\n{matrix.iloc[:min(3, matrix.shape[0]), :min(3, matrix.shape[1])]}")
    else:
        # For numpy ndarray
        if matrix.size == 0:
            print("空の行列です")
            return
        print(f"Shape: {matrix.shape}")
        print(f"Min value: {np.min(matrix)}")
        print(f"Max value: {np.max(matrix)}")
        print(f"Mean value: {np.mean(matrix)}")
        print(f"Std value: {np.std(matrix)}")
        print(f"Non-zero entries: {np.count_nonzero(matrix)} / {matrix.size}")
        print(f"Sample values:\n{matrix[:min(3, matrix.shape[0]), :min(3, matrix.shape[1])]}")
    
    print("================\n")


def plot_circle_communication(network_summary, title="Cell-Cell Communication Network", figsize=(10, 10)):
    """
    サーキュラープロットで細胞間通信を描画
    
    Parameters
    ----------
    network_summary : dict
        ネットワーク集計結果
    title : str
        プロットのタイトル
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:
        fig, ax = plt.subplots(figsize=figsize)
        
        # 強度行列を取得
        if 'strength_matrix' not in network_summary or isinstance(network_summary['strength_matrix'], pd.DataFrame) and network_summary['strength_matrix'].empty:
            ax.text(0.5, 0.5, "データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title(title)
            return fig
            
        matrix = network_summary['strength_matrix']
        
        # 行列を正規化
        if isinstance(matrix, pd.DataFrame):
            if matrix.sum().sum() > 0:
                norm_matrix = matrix / matrix.max().max()
            else:
                norm_matrix = matrix
        else:
            if np.sum(matrix) > 0:
                norm_matrix = matrix / np.max(matrix)
            else:
                norm_matrix = matrix
        
        # カラーマップを設定
        colors = plt.cm.Set3(np.linspace(0, 1, len(matrix)))
        
        # サーキュラープロットの描画（簡易版）
        # 実際のCellChatではcirclizeなどのパッケージを使ってより洗練されたプロットを作成
        plt.title(title)
        plt.axis('equal')
        plt.axis('off')
        
        # 簡易的なサーキュラープロットのダミー実装
        theta = np.linspace(0, 2*np.pi, len(matrix)+1)[:-1]
        
        # ノードを配置
        x = np.cos(theta)
        y = np.sin(theta)
        
        # ノードの描画
        if isinstance(matrix, pd.DataFrame):
            labels = matrix.index
        else:
            labels = [f"Cluster{i}" for i in range(len(matrix))]
            
        for i, (xi, yi, label) in enumerate(zip(x, y, labels)):
            ax.scatter(xi, yi, s=300, color=colors[i], edgecolor='black', zorder=10)
            ax.text(xi*1.15, yi*1.15, label, ha='center', va='center', fontsize=12, fontweight='bold')
        
        # エッジの描画
        for i, sender in enumerate(range(len(matrix))):
            for j, receiver in enumerate(range(len(matrix))):
                # DataFrameとNumpyの両方に対応
                if isinstance(matrix, pd.DataFrame):
                    sender_label = matrix.index[i]
                    receiver_label = matrix.columns[j]
                    value = matrix.loc[sender_label, receiver_label]
                    if i != j and value > 0:
                        strength = norm_matrix.loc[sender_label, receiver_label]
                else:
                    if i != j and matrix[i, j] > 0:
                        strength = norm_matrix.iloc[i, j]
                    else:
                        continue
                
                # ベジェ曲線で弧を描画
                xi, yi = x[i], y[i]
                xj, yj = x[j], y[j]
                # 中間点を少しずらす
                xm = (xi + xj) / 2
                ym = (yi + yj) / 2
                # 中心から離す
                dx = xm
                dy = ym
                d = np.sqrt(dx**2 + dy**2)
                xm += dx / d * 0.3
                ym += dy / d * 0.3
                
                # 太さは強度に比例
                width = 1 + 4 * strength
                ax.plot([xi, xm, xj], [yi, ym, yj], 'gray', linewidth=width, alpha=0.6)
                
                # 矢印を描画
                ax.arrow(xm, ym, (xj-xm)*0.2, (yj-ym)*0.2, width=0.01*strength, 
                        head_width=0.05*width, head_length=0.1*width, 
                        fc='black', ec='black', zorder=5)
        
        return fig
    except Exception as e:
        logger.error(f"サーキュラープロット作成エラー: {str(e)}")
        logger.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"サーキュラープロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title(title)
        return fig




def plot_dot_lr_network(results_df, source_cells, target_cells, top_n=20, pval_threshold=0.05, figsize=(12, 10)):
    """
    ドットプロットで細胞間のLR相互作用を描画
    
    Parameters
    ----------
    results_df : pd.DataFrame
        CellChatの計算結果
    source_cells : list
        送信側の細胞タイプリスト
    target_cells : list
        受信側の細胞タイプリスト
    top_n : int
        表示する上位ペアの数
    pval_threshold : float
        有意とみなすP値の閾値
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:
        if results_df.empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効な相互作用データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
            
        # 選択した細胞タイプをフィルタリング
        filtered_df = results_df[
            (results_df['source'].isin(source_cells)) & 
            (results_df['target'].isin(target_cells)) &
            (results_df['pval'] <= pval_threshold)
        ]
        
        if len(filtered_df) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"指定した条件で有効な相互作用がありません。\n選択した細胞タイプ間では有意な相互作用が検出されませんでした。\nP値閾値を上げるか、別の細胞タイプを選択してください。", 
                   ha='center', va='center', wrap=True)
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
        
        # LRペアごとに集計
        if 'interaction_prob_normalized' not in filtered_df.columns:
            filtered_df['interaction_prob_normalized'] = filtered_df['prob']
            
        lr_summary = filtered_df.groupby(['ligand', 'receptor'])['interaction_prob_normalized'].sum().reset_index()
        lr_summary = lr_summary.sort_values('interaction_prob_normalized', ascending=False).head(top_n)
        
        # リガンド-レセプターペアを結合
        lr_summary['interaction'] = lr_summary['ligand'] + '-' + lr_summary['receptor']
        
        # 各ペアの細胞タイプペアごとの強度を計算
        dot_data = []
        
        for lr_pair in lr_summary[['ligand', 'receptor']].itertuples(index=False):
            ligand, receptor = lr_pair
            
            for source in source_cells:
                for target in target_cells:
                    # このLRペアとこの細胞タイプペアに対応する行を抽出
                    subset = filtered_df[
                        (filtered_df['ligand'] == ligand) & 
                        (filtered_df['receptor'] == receptor) &
                        (filtered_df['source'] == source) & 
                        (filtered_df['target'] == target)
                    ]
                    
                    if len(subset) > 0:
                        row = subset.iloc[0]
                        dot_data.append({
                            'interaction': f"{ligand}-{receptor}",
                            'source': source,
                            'target': target,
                            'strength': row['interaction_prob_normalized'],
                            'pvalue': row['pval']
                        })
        
        dot_df = pd.DataFrame(dot_data)
        
        if len(dot_df) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "ドットプロット用のデータがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
        
        # ソート
        unique_interactions = lr_summary['interaction'].tolist()
        
        # プロット作成
        fig, ax = plt.subplots(figsize=figsize)
        
        # ドットプロットのデータ準備
        plot_data = {}
        
        for source in source_cells:
            for target in target_cells:
                cell_pair = f"{source}->{target}"
                plot_data[cell_pair] = {}
                
                source_target_df = dot_df[(dot_df['source'] == source) & (dot_df['target'] == target)]
                
                for interaction in unique_interactions:
                    interaction_df = source_target_df[source_target_df['interaction'] == interaction]
                    
                    if len(interaction_df) > 0:
                        plot_data[cell_pair][interaction] = {
                            'strength': interaction_df['strength'].values[0],
                            'pvalue': interaction_df['pvalue'].values[0]
                        }
                    else:
                        plot_data[cell_pair][interaction] = {
                            'strength': 0,
                            'pvalue': 1
                        }
        
        # 各セルペアに対してプロット
        cell_pairs = list(plot_data.keys())
        n_pairs = len(cell_pairs)
        
        if n_pairs == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "ドットプロット用のデータがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Interactions')
            return fig
            
        # 縦にサブプロットを作成
        fig, axes = plt.subplots(n_pairs, 1, figsize=(figsize[0], figsize[1] * n_pairs / 3), sharex=True)
        
        if n_pairs == 1:
            axes = [axes]
        
        for i, (cell_pair, ax) in enumerate(zip(cell_pairs, axes)):
            strengths = []
            sizes = []
            interactions = []
            
            for interaction in unique_interactions:
                data = plot_data[cell_pair][interaction]
                strengths.append(data['strength'])
                # サイズは-log10(pvalue)に比例
                sizes.append(max(20, -np.log10(data['pvalue'] + 1e-10) * 20))
                interactions.append(interaction)
            
            # ゼロでない強度のみプロット
            non_zero = [i for i, s in enumerate(strengths) if s > 0]
            
            if non_zero:
                # カラーマップの設定
                max_strength = max(np.array(strengths)[non_zero]) if non_zero else 1
                colors = plt.cm.YlOrRd(np.array(strengths)[non_zero] / max_strength)
                
                # 散布図プロット
                scatter = ax.scatter(
                    [interactions[i] for i in non_zero],
                    [0] * len(non_zero),
                    s=[sizes[i] for i in non_zero],
                    c=colors,
                    alpha=0.7,
                    edgecolors='gray'
                )
                
                # Y軸ラベルをセルペア名に
                ax.set_ylabel(cell_pair, fontsize=12)
                
                # Y軸目盛りを非表示
                ax.set_yticks([])
            else:
                ax.text(0.5, 0, f"有意な相互作用なし", ha='center', va='center')
                ax.set_ylabel(cell_pair, fontsize=12)
                ax.set_yticks([])
            
            # 最後のサブプロットのみX軸ラベルを回転
            if i == len(axes) - 1:
                plt.xticks(rotation=90, ha='right')
                
        plt.suptitle('Ligand-Receptor Interactions Between Cell Types', fontsize=14)
        plt.tight_layout()
        plt.subplots_adjust(top=0.95)
        
        return fig
    except Exception as e:
        logger.error(f"ドットプロット作成エラー: {str(e)}")
        logger.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"ドットプロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Ligand-Receptor Interactions')
        return fig

def plot_aggregated_network(network_summary, threshold=0.01, figsize=(12, 10)):
    """
    集計されたネットワークをグラフ形式で描画
    
    Parameters
    ----------
    network_summary : dict
        ネットワーク集計結果
    threshold : float
        表示する最小相互作用強度
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:
        # 強度行列を取得
        if 'strength_matrix' not in network_summary or isinstance(network_summary['strength_matrix'], pd.DataFrame) and network_summary['strength_matrix'].empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効な相互作用データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Cell-Cell Communication Network')
            return fig
            
        matrix = network_summary['strength_matrix']
        
        # しきい値以下をフィルタリング
        if isinstance(matrix, pd.DataFrame):
            filtered_matrix = matrix.copy()
            filtered_matrix[filtered_matrix < threshold] = 0
        else:
            filtered_matrix = matrix.copy()
            filtered_matrix[filtered_matrix < threshold] = 0
        
        # グラフを作成
        if isinstance(filtered_matrix, pd.DataFrame):
            G = nx.from_pandas_adjacency(filtered_matrix, create_using=nx.DiGraph())
        else:
            # 安定のため、行列を numpy 配列に変換してからグラフ生成
            adj_mat = np.array(filtered_matrix)
            G = nx.from_numpy_array(adj_mat, create_using=nx.DiGraph())
        # 孤立ノードを削除
        G.remove_nodes_from(list(nx.isolates(G)))
        
        if len(G) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"閾値 {threshold} 以上の相互作用がありません。\n閾値を下げてみてください。", 
                   ha='center', va='center')
            ax.axis('off')
            plt.title('Cell-Cell Communication Network')
            return fig
        
        # エッジの重みをセット
        if isinstance(filtered_matrix, pd.DataFrame):
            for u, v, d in G.edges(data=True):
                d['weight'] = filtered_matrix.loc[u, v]
        else:
            for u, v, d in G.edges(data=True):
                d['weight'] = filtered_matrix[u, v]
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # ノードの色を設定
        node_color = plt.cm.Set3(np.linspace(0, 1, len(G)))[:len(G)]
        
        # ノードの位置を決定
        pos = nx.spring_layout(G, seed=42)
        
        # エッジの太さを重みに比例させる
        if isinstance(filtered_matrix, pd.DataFrame):
            max_weight = filtered_matrix.max().max()
        else:
            max_weight = np.max(filtered_matrix)
            
        if max_weight > 0:
            edge_weights = [d['weight'] * 5 / max_weight for u, v, d in G.edges(data=True)]
        else:
            edge_weights = [1.0 for u, v, d in G.edges(data=True)]
        
        # ノードを描画
        nx.draw_networkx_nodes(G, pos, node_size=700, node_color=node_color, alpha=0.8, ax=ax)
        
        # エッジを描画
        nx.draw_networkx_edges(
            G, pos, width=edge_weights, alpha=0.6, 
            edge_color='gray', arrows=True, 
            arrowstyle='->', arrowsize=20, ax=ax
        )
        
        # ラベルを描画
        nx.draw_networkx_labels(G, pos, font_size=12, font_family='sans-serif', ax=ax)

        # タイトルと調整
        plt.title('Cell-Cell Communication Network', fontsize=14)
        plt.axis('off')
        plt.tight_layout()
        
        return fig
    except Exception as e:
        logger.error(f"ネットワーク図作成エラー: {str(e)}")
        logger.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"ネットワーク図作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Cell-Cell Communication Network')
        return fig


def plot_pathway_enrichment(results_df, top_n=20, figsize=(12, 8)):
    """
    パスウェイエンリッチメント解析のシミュレーション
    
    実際のCellChatではシグナリングパスウェイの情報が含まれていますが、
    この実装では簡易的にリガンドのプレフィックスに基づいて模擬的なパスウェイ分類を行います。
    
    Parameters
    ----------
    results_df : pd.DataFrame
        CellChatの計算結果
    top_n : int
        表示する上位パスウェイの数
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:
        if results_df.empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効な相互作用データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Pathway Enrichment Analysis')
            return fig
            
        # リガンドプレフィックスからパスウェイを疑似的に決定
        def get_pathway(ligand):
            if not isinstance(ligand, str):
                return 'Other'
                
            if ligand.startswith(('IL', 'Il')):
                return 'Interleukin'
            elif ligand.startswith(('TNF', 'Tnf')):
                return 'TNF'
            elif ligand.startswith(('CCL', 'Ccl', 'CXCL', 'Cxcl')):
                return 'Chemokine'
            elif ligand.startswith(('WNT', 'Wnt')):
                return 'WNT'
            elif ligand.startswith(('TGF', 'Tgf')):
                return 'TGF-beta'
            elif ligand.startswith(('FGF', 'Fgf')):
                return 'FGF'
            elif ligand.startswith(('EGF', 'Egf')):
                return 'EGF'
            elif ligand.startswith(('PDGF', 'Pdgf')):
                return 'PDGF'
            elif ligand.startswith(('IGF', 'Igf')):
                return 'IGF'
            else:
                return 'Other'
        
        # 有意な相互作用をフィルタリング
        if 'pval' in results_df.columns:
            sig_df = results_df[results_df['pval'] <= 0.05].copy()
        else:
            sig_df = results_df.copy()
        
        if len(sig_df) == 0:
            # 有意な相互作用がなければ全データを使用
            sig_df = results_df.copy()
            
        # パスウェイ情報を追加
        sig_df['pathway'] = sig_df['ligand'].apply(get_pathway)
        
        # パスウェイごとに集計
        if 'interaction_prob_normalized' in sig_df.columns:
            pathway_strength = sig_df.groupby('pathway')['interaction_prob_normalized'].sum().reset_index()
        else:
            pathway_strength = sig_df.groupby('pathway')['prob'].sum().reset_index()
            pathway_strength.rename(columns={'prob': 'interaction_prob_normalized'}, inplace=True)
            
        pathway_strength = pathway_strength.sort_values('interaction_prob_normalized', ascending=False).head(top_n)
        
        if len(pathway_strength) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なパスウェイデータがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Pathway Enrichment Analysis')
            return fig
            
        fig, ax = plt.subplots(figsize=figsize)
        
        # 棒グラフの描画
        bars = ax.barh(
            pathway_strength['pathway'],
            pathway_strength['interaction_prob_normalized'],
            color=plt.cm.Spectral(np.linspace(0, 1, len(pathway_strength)))
        )
        
        ax.set_xlabel('Aggregated Interaction Strength')
        ax.set_title('Pathway Enrichment Analysis')
        
        plt.tight_layout()
        
        return fig
    except Exception as e:
        logger.error(f"パスウェイ解析プロット作成エラー: {str(e)}")
        logger.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"パスウェイ解析プロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Pathway Enrichment Analysis')
        return fig


def preprocess_data(
    adata, 
    group_by, 
    db=None,
    min_cells=10, 
    thresh_pct=0.1,
    thresh_p = 0.05,
    use_gpu=True
):
    """
    CellChatのデータ前処理を行う関数
    
    Parameters:
    ----------
    adata : AnnData
        AnnDataオブジェクト
    group_by : str
        細胞グループを定義するobsのカラム名
    db : dict, optional
        CellChatデータベース
    min_cells : int
        各細胞グループに必要な最小細胞数
    thresh_pct : float
        発現と見なすための発現率の閾値（0-1の間）
    
    Returns:
    --------
    adata_filtered : AnnData
        フィルタリング済みデータ
    """
    st.write(f"original data: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # 1. 細胞グループの細胞数チェック
    cell_counts = adata.obs[group_by].value_counts()
    valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
    if len(valid_groups) < len(cell_counts):
        st.warning(f"{len(cell_counts) - len(valid_groups)}個の細胞グループが{min_cells}細胞未満のため除外されました")
    
    # 有効なグループの細胞だけを保持
    adata_filtered = adata[adata.obs[group_by].isin(valid_groups)].copy()
    
    # 2. 高発現遺伝子の同定


    overexpressed_genes_result = identify_overexpressed_genes_optimized(
        adata_filtered,
        group_by=group_by,
        do_de=True,
        do_fast=True,
        thresh_p=thresh_p,
        thresh_pct=thresh_pct,
        min_cells=min_cells,
        only_pos=True,
        use_gpu=use_gpu,  # 新しいパラメータ
        batch_size=5000,  # 新しいパラメータ
        n_jobs=n_cpus     # UIから渡されるCPU数
    )
    
    # 高発現遺伝子を抽出
    features_sig = overexpressed_genes_result['features']
    st.write(f"{len(features_sig)} genes passed filtering")
    
    # 3. 発現遺伝子のみを保持
    adata_filtered = adata_filtered[:, features_sig].copy()

    st.write(f"filtered data: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} genes")
    
    # 4. シグナリング関連遺伝子のみを抽出（データベースが提供されている場合）
    if db is not None:
        # リガンド、レセプター、および関連遺伝子を抽出
        signaling_genes = []
        
        # 単一遺伝子リガンド・レセプター
        signaling_genes.extend(db['interaction']['ligand'].astype(str).tolist())
        signaling_genes.extend(db['interaction']['receptor'].astype(str).tolist())
        
        # 複合体のサブユニットを抽出
        if not db['complex'].empty:
            for _, row in db['complex'].iterrows():
                for col in [c for c in db['complex'].columns if 'subunit' in c]:
                    if pd.notna(row[col]) and row[col] != "":
                        signaling_genes.append(str(row[col]))
        
        # 共受容体などを抽出
        if not db['cofactor'].empty:
            for _, row in db['cofactor'].iterrows():
                for col in [c for c in db['cofactor'].columns if 'cofactor' in c]:
                    if pd.notna(row[col]) and row[col] != "":
                        signaling_genes.append(str(row[col]))
        
        # 重複を削除し、実際にデータに存在する遺伝子のみを保持
        signaling_genes = list(set(signaling_genes))
        valid_signaling_genes = [g for g in signaling_genes if g in features_sig]
        
        print(f"{len(valid_signaling_genes)}個のシグナリング関連遺伝子がデータに存在")
        
        # シグナリング遺伝子のみを保持
        adata_filtered = adata_filtered[:, valid_signaling_genes].copy()
    
    print(f"前処理後: {adata_filtered.shape[0]}遺伝子 x {adata_filtered.shape[1]}細胞")
    
    return adata_filtered


# プロット用のユーティリティ関数
def scPalette(n):
    """
    Generate colors from a customed color palette
    
    Parameters:
    ----------
    n : int
        Number of colors
        
    Returns:
    -------
    colors : list
        A color palette for plotting
    """
    colorSpace = ['#E41A1C','#377EB8','#4DAF4A','#984EA3','#F29403','#F781BF','#BC9DCC',
                 '#A65628','#54B0E4','#222F75','#1B9E77','#B2DF8A','#E3BE00','#FB9A99',
                 '#E7298A','#910241','#00CDD1','#A6CEE3','#CE1261','#5E4FA2','#8CA77B',
                 '#00441B','#DEDC00','#DCF0B9','#8DD3C7','#999999']
    if n <= len(colorSpace):
        colors = colorSpace[:n]
    else:
        # colorRampPaletteに相当する関数
        from matplotlib.colors import LinearSegmentedColormap
        cmap = LinearSegmentedColormap.from_list('custom_cmap', colorSpace)
        colors = [cmap(i / n) for i in range(n)]
    return colors

def ggPalette(n):
    """
    Generate ggplot2 colors
    
    Parameters:
    ----------
    n : int
        Number of colors to generate
        
    Returns:
    -------
    colors : list
        A color palette for plotting
    """
    from matplotlib.colors import hsv_to_rgb
    hues = np.linspace(15, 375, n + 1) / 360
    saturations = np.ones(n + 1) * 0.65
    values = np.ones(n + 1) * 0.65
    hsv_colors = np.column_stack((hues, saturations, values))
    rgb_colors = [hsv_to_rgb(hsv) for hsv in hsv_colors]
    return rgb_colors[:n]

# コードダイアグラム可視化のための関数
def plot_chord_cell(net, color_use=None, group=None, cell_order=None, 
                   sources_use=None, targets_use=None, lab_cex=0.8,
                   small_gap=1, big_gap=10, remove_isolate=False, 
                   transparency=0.4, title_name=None):
    """
    Chord diagram for visualizing cell-cell communication
    
    Parameters
    ----------
    net : numpy.ndarray or pandas.DataFrame
        A weighted matrix or a dataframe with three columns defining the cell-cell communication network
    color_use : dict or list, optional
        Colors for the cell groups
    group : dict, optional
        A named group labels for making multiple-group Chord diagrams
    cell_order : list, optional
        A vector defining the cell type orders
    sources_use : list, optional
        A vector giving the index or the name of source cell groups
    targets_use : list, optional
        A vector giving the index or the name of target cell groups
    lab_cex : float, optional
        Font size for the text
    small_gap : int, optional
        Small gap between sectors
    big_gap : int, optional
        Gap between the different sets of sectors
    remove_isolate : bool, optional
        Whether remove the isolate nodes in the communication network
    transparency : float, optional
        Transparency of link colors
    title_name : str, optional
        Title of the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Please install plotly: pip install plotly")
        return None
    
    # 入力がDataFrameの場合、行列に変換
    if isinstance(net, pd.DataFrame):
        if all(c in net.columns for c in ["source", "target", "prob"]):
            cell_levels = sorted(list(set(net['source'].unique()) | set(net['target'].unique())))
            net_matrix = np.zeros((len(cell_levels), len(cell_levels)))
            for _, row in net.iterrows():
                src_idx = cell_levels.index(row['source'])
                tgt_idx = cell_levels.index(row['target'])
                net_matrix[src_idx, tgt_idx] = row['prob']
            net = net_matrix
            
    # ソース・ターゲットのフィルタリング
    if sources_use is not None or targets_use is not None:
        cells_level = np.arange(net.shape[0])
        df_net = pd.DataFrame()
        for i in range(net.shape[0]):
            for j in range(net.shape[1]):
                if net[i, j] > 0:
                    df_net = df_net.append({
                        'source': i,
                        'target': j,
                        'value': net[i, j]
                    }, ignore_index=True)
        
        if sources_use is not None:
            if isinstance(sources_use[0], int):
                df_net = df_net[df_net['source'].isin(sources_use)]
            else:
                indices = [cell_levels.index(s) for s in sources_use if s in cell_levels]
                df_net = df_net[df_net['source'].isin(indices)]
        
        if targets_use is not None:
            if isinstance(targets_use[0], int):
                df_net = df_net[df_net['target'].isin(targets_use)]
            else:
                indices = [cell_levels.index(t) for t in targets_use if t in cell_levels]
                df_net = df_net[df_net['target'].isin(indices)]
        
        # 行列に変換し直す
        if not df_net.empty:
            net = np.zeros((len(cell_levels), len(cell_levels)))
            for _, row in df_net.iterrows():
                net[int(row['source']), int(row['target'])] = row['value']
    
    # 孤立ノードの削除（オプション）
    if remove_isolate:
        idx1 = np.where(np.sum(net, axis=0) == 0)[0]
        idx2 = np.where(np.sum(net, axis=1) == 0)[0]
        idx = np.intersect1d(idx1, idx2)
        if len(idx) > 0:
            net = np.delete(np.delete(net, idx, axis=0), idx, axis=1)
    
    # カラー設定
    if color_use is None:
        color_use = scPalette(net.shape[0])
    
    # ラベル設定
    if cell_order is None:
        labels = [f"Cell{i+1}" for i in range(net.shape[0])]
    else:
        labels = cell_order
    
    # Plotlyを使用したコードダイアグラムの作成
    source_list = []
    target_list = []
    value_list = []
    color_list = []
    
    for i in range(net.shape[0]):
        for j in range(net.shape[1]):
            if net.iloc[i, j] > 0:
                source_list.append(i)
                target_list.append(j)
                value_list.append(net.iloc[i, j])
                # 色はソースの色を使用
                if isinstance(color_use, list):
                    color_list.append(color_use[i % len(color_use)])
                else:
                    color_list.append(color_use.get(i, '#377EB8'))
    
    # カラーを16進形式に変換（必要な場合）
    for i in range(len(color_list)):
        if isinstance(color_list[i], tuple) and len(color_list[i]) in [3, 4]:
            if max(color_list[i]) <= 1:
                color_list[i] = f'rgba({int(color_list[i][0]*255)},{int(color_list[i][1]*255)},{int(color_list[i][2]*255)},{transparency if len(color_list[i]) == 3 else color_list[i][3]})'
            else:
                color_list[i] = f'rgba({int(color_list[i][0])},{int(color_list[i][1])},{int(color_list[i][2])},{transparency if len(color_list[i]) == 3 else color_list[i][3]/255})'
    
    fig = go.Figure(data=[go.Sankey(
        node=dict(
            pad=20,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=source_list,
            target=target_list,
            value=value_list,
            color=color_list
        )
    )])
    
    fig.update_layout(
        title=title_name,
        font_size=lab_cex * 10
    )
    
    return fig

def plot_chord_gene(object, signaling=None, pairLR_use=None, slot_name="net", 
                   color_use=None, thresh=0.05, title_name=None):
    """
    Chord diagram for visualizing cell-cell communication for a set of ligands/receptors or signaling pathways
    
    Parameters
    ----------
    object : CellChat object
        CellChat object
    signaling : list, optional
        A character vector giving the name of signaling networks
    pairLR_use : DataFrame, optional
        A data frame defining the L-R pairs of interest
    slot_name : str, optional
        The slot name of object
    color_use : dict or list, optional
        Colors for the cell groups
    thresh : float, optional
        Threshold of the p-value for determining significant interaction
    title_name : str, optional
        Title of the plot
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("Please install plotly: pip install plotly")
        return None
    
    # プレースホルダー関数 - 実際には詳細な実装が必要
    # 以下はサンプル実装の概要
    
    # プロ部（prob）と閾値情報（pval）の取得
    prob = getattr(object, slot_name).get("prob", np.array([]))
    pval = getattr(object, slot_name).get("pval", np.array([]))
    
    if prob.size == 0 or pval.size == 0:
        print("No probability or p-value data available")
        return None
    
    # 閾値でフィルタリング
    prob[pval > thresh] = 0
    
    # シグナリングまたはL-Rペアのフィルタリング
    if signaling is not None:
        # シグナリング関連のL-Rペアを選択
        pass
    elif pairLR_use is not None:
        # 特定のL-Rペアを選択
        pass
    else:
        print("Either signaling or pairLR_use must be provided")
        return None
    
    # コードダイアグラムの作成
    # 実際の実装では、特定のシグナリングパスウェイまたはL-Rペアに基づいてネットワークデータを作成
    
    # サンプルとして、ランダムなネットワークデータを生成
    n_nodes = 10
    net = np.random.rand(n_nodes, n_nodes)
    net[net < 0.7] = 0
    
    # plot_chord_cell関数を使用して可視化
    return plot_chord_cell(net, color_use=color_use, title_name=title_name)

# リバープロットのための関数
def plot_river(object, slot_name="netP", pattern="outgoing", cutoff=0.5,
              sources_use=None, targets_use=None, signaling=None,
              color_use=None, color_use_pattern=None, color_use_signaling="grey50",
              do_order=False, main_title=None, font_size=2.5, font_size_title=12):
    """
    River plot showing the associations of latent patterns with cell groups and signaling pathways
    
    Parameters
    ----------
    object : CellChat object
        CellChat object
    slot_name : str, optional
        The slot name of object
    pattern : str, optional
        "outgoing" or "incoming"
    cutoff : float, optional
        The threshold for filtering out weak links
    sources_use : list, optional
        A vector giving the index or the name of source cell groups
    targets_use : list, optional
        A vector giving the index or the name of target cell groups
    signaling : list, optional
        A vector giving the name of signaling pathways
    color_use : dict or list, optional
        Colors for the cell groups
    color_use_pattern : dict or list, optional
        Colors for the patterns
    color_use_signaling : str or list, optional
        Colors for the signaling
    do_order : bool, optional
        Whether reorder the cell groups or signaling
    main_title : str, optional
        The title of plot
    font_size : float, optional
        Font size of the text
    font_size_title : float, optional
        Font size of the title
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    try:
        import plotly.graph_objects as go
        from plotly.subplots import make_subplots
    except ImportError:
        print("Please install plotly: pip install plotly")
        return None
    
    # パターン情報の取得
    res_pattern = getattr(object, slot_name, {}).get("pattern", {}).get(pattern, {})
    data1 = res_pattern.get("pattern", {}).get("cell", None)
    data2 = res_pattern.get("pattern", {}).get("signaling", None)
    
    if data1 is None:
        print(f"No pattern data found for {pattern}")
        return None
    
    # カラー設定
    if color_use_pattern is None:
        n_patterns = len(np.unique(data1["Pattern"]))
        colors = ggPalette(n_patterns * 2)
        if pattern == "outgoing":
            color_use_pattern = [colors[i] for i in range(0, n_patterns * 2, 2)]
        else:
            color_use_pattern = [colors[i] for i in range(1, n_patterns * 2, 2)]
    
    if main_title is None:
        if pattern == "outgoing":
            main_title = "Outgoing communication patterns of secreting cells"
        else:
            main_title = "Incoming communication patterns of target cells"
    
    # データ前処理
    data1_filtered = data1.copy()
    data1_filtered.loc[data1_filtered["Contribution"] < cutoff, "Contribution"] = 0
    
    if data2 is not None:
        data2_filtered = data2.copy()
        data2_filtered.loc[data2_filtered["Contribution"] < cutoff, "Contribution"] = 0
        
        # パターンとシグナリングパスウェイの関連を可視化
        # Part 1: CellGroup -> Pattern
        # Part 2: Pattern -> Signaling
        
        # Plotlyを使用したSankeyダイアグラム
        # ノードを準備
        cell_nodes = np.unique(data1_filtered["CellGroup"])
        pattern_nodes = np.unique(data1_filtered["Pattern"])
        signaling_nodes = np.unique(data2_filtered["Signaling"])
        
        all_nodes = np.concatenate([cell_nodes, pattern_nodes, signaling_nodes])
        
        # リンクを準備（CellGroup -> Pattern）
        source_indices = []
        target_indices = []
        values = []
        colors = []
        
        for _, row in data1_filtered[data1_filtered["Contribution"] > 0].iterrows():
            source_idx = np.where(all_nodes == row["CellGroup"])[0][0]
            target_idx = np.where(all_nodes == row["Pattern"])[0][0]
            source_indices.append(source_idx)
            target_indices.append(target_idx)
            values.append(row["Contribution"])
            
            # 色の設定
            if color_use is not None and row["CellGroup"] in color_use:
                c = color_use[row["CellGroup"]]
                if isinstance(c, tuple):
                    colors.append(f'rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.7)')
                else:
                    colors.append(c)
            else:
                colors.append('rgba(100,100,100,0.7)')  # デフォルト色
        
        # リンクを準備（Pattern -> Signaling）
        for _, row in data2_filtered[data2_filtered["Contribution"] > 0].iterrows():
            source_idx = np.where(all_nodes == row["Pattern"])[0][0]
            target_idx = np.where(all_nodes == row["Signaling"])[0][0]
            source_indices.append(source_idx)
            target_indices.append(target_idx)
            values.append(row["Contribution"])
            
            # 色の設定
            if color_use_pattern is not None:
                pattern_idx = list(pattern_nodes).index(row["Pattern"])
                c = color_use_pattern[pattern_idx % len(color_use_pattern)]
                if isinstance(c, tuple):
                    colors.append(f'rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.7)')
                else:
                    colors.append(c)
            else:
                colors.append('rgba(100,100,100,0.7)')  # デフォルト色
        
        # Sankeyダイアグラムを作成
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_nodes
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=colors
            )
        )])
        
        fig.update_layout(
            title=main_title,
            font=dict(size=font_size_title)
        )
        
        return fig
    else:
        # シグナリングデータがない場合、CellGroup -> Pattern のみを可視化
        cell_nodes = np.unique(data1_filtered["CellGroup"])
        pattern_nodes = np.unique(data1_filtered["Pattern"])
        
        all_nodes = np.concatenate([cell_nodes, pattern_nodes])
        
        # リンクを準備
        source_indices = []
        target_indices = []
        values = []
        colors = []
        
        for _, row in data1_filtered[data1_filtered["Contribution"] > 0].iterrows():
            source_idx = np.where(all_nodes == row["CellGroup"])[0][0]
            target_idx = np.where(all_nodes == row["Pattern"])[0][0]
            source_indices.append(source_idx)
            target_indices.append(target_idx)
            values.append(row["Contribution"])
            
            # 色の設定
            if color_use is not None and row["CellGroup"] in color_use:
                c = color_use[row["CellGroup"]]
                if isinstance(c, tuple):
                    colors.append(f'rgba({int(c[0]*255)},{int(c[1]*255)},{int(c[2]*255)},0.7)')
                else:
                    colors.append(c)
            else:
                colors.append('rgba(100,100,100,0.7)')  # デフォルト色
        
        # Sankeyダイアグラムを作成
        fig = go.Figure(data=[go.Sankey(
            node=dict(
                pad=15,
                thickness=20,
                line=dict(color="black", width=0.5),
                label=all_nodes
            ),
            link=dict(
                source=source_indices,
                target=target_indices,
                value=values,
                color=colors
            )
        )])
        
        fig.update_layout(
            title=main_title,
            font=dict(size=font_size_title)
        )
        
        return fig

# マニフォールド可視化のための関数
def plot_embedding(object, slot_name="netP", type="functional", color_use=None,
                  pathway_labeled=None, top_label=1, pathway_remove=None, pathway_remove_show=True,
                  dot_size=(2, 6), dot_alpha=0.5, label_size=8, title_name=None,
                  font_size=10, font_size_title=12, do_label=True, show_legend=True, show_axes=True):
    """
    2D visualization of the learned manifold of signaling networks
    
    Parameters
    ----------
    object : CellChat object
        CellChat object
    slot_name : str, optional
        The slot name of object
    type : str, optional
        "functional" or "structural"
    color_use : dict or list, optional
        Colors for the cell groups
    pathway_labeled : list, optional
        A vector giving the signaling names to show
    top_label : float, optional
        The fraction of signaling pathways to label
    pathway_remove : list, optional
        A vector defining the signaling to remove
    pathway_remove_show : bool, optional
        Whether show the removed signaling names
    dot_size : tuple, optional
        A range defining the size of the symbol
    dot_alpha : float, optional
        Transparency
    label_size : float, optional
        Font size of the text
    title_name : str, optional
        The title of plot
    font_size : float, optional
        Font size of the text
    font_size_title : float, optional
        Font size of the title
    do_label : bool, optional
        Whether label the each point
    show_legend : bool, optional
        Whether show the legend
    show_axes : bool, optional
        Whether show the axes
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    # 比較設定
    comparison = "single"
    comparison_name = comparison
    
    # データ取得
    try:
        Y = getattr(object, slot_name).similarity[type]["dr"][comparison_name]
        Groups = getattr(object, slot_name).similarity[type]["group"][comparison_name]
        prob = getattr(object, slot_name).prob
    except (AttributeError, KeyError):
        print(f"Required data not found in the object. Check that object.{slot_name}.similarity.{type} exists.")
        return None
    
    # パスウェイ除外処理
    if pathway_remove is None:
        try:
            similarity = getattr(object, slot_name).similarity[type]["matrix"][comparison_name]
            # 孤立したノードを特定（列の和が1）
            pathway_remove = [row for row, col_sum in zip(similarity.index, np.sum(similarity, axis=1).tolist()) if col_sum == 1]
        except (AttributeError, KeyError):
            pathway_remove = []
    
    # パスウェイの除外
    if len(pathway_remove) > 0:
        try:
            pathway_remove_idx = [i for i, path in enumerate(prob.shape[2]) if path in pathway_remove]
            if pathway_remove_idx:
                prob = np.delete(prob, pathway_remove_idx, axis=2)
        except (IndexError, AttributeError):
            print("Error when trying to remove pathways")
    
    # 確率の合計計算
    prob_sum = np.sum(prob, axis=(0, 1))
    
    # データフレーム作成
    pathways = np.array([str(i) for i in range(prob.shape[2])])  # 実際のパスウェイ名がない場合のデモ
    
    df = pd.DataFrame({
        'x': Y[:, 0],
        'y': Y[:, 1],
        'Commun.Prob.': prob_sum/np.max(prob_sum) if np.max(prob_sum) > 0 else prob_sum,
        'labels': pathways,
        'Groups': np.array(Groups, dtype=str)
    })
    
    # カラー設定
    if color_use is None:
        unique_groups = np.unique(Groups)
        color_use = dict(zip(unique_groups, ggPalette(len(unique_groups))))
    
    # プロット作成
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 各グループごとに散布図をプロット
    for g in np.unique(df['Groups']):
        group_df = df.loc[df['Groups'] == g].copy()
        scatter = ax.scatter(
            group_df['x'], 
            group_df['y'],
            s=group_df['Commun.Prob.'] * (dot_size[1] - dot_size[0]) + dot_size[0],
            c=[color_use[g] if isinstance(color_use, dict) else color_use[list(np.unique(df['Groups'])).index(g)]],
            alpha=dot_alpha,
            label=g
        )
    
    # ラベル追加
    if do_label:
        if pathway_labeled is None:
            if top_label < 1:
                # トップNのパスウェイを選択
                top_n = int(top_label * len(df))
                data_label = df.nlargest(top_n, 'Commun.Prob.')
            else:
                data_label = df
        else:
            data_label = df[df['labels'].isin(pathway_labeled)]
        
        for _, row in data_label.iterrows():
            ax.annotate(
                row['labels'],
                (row['x'], row['y']),
                fontsize=label_size,
                color=color_use[row['Groups']] if isinstance(color_use, dict) else color_use[list(np.unique(df['Groups'])).index(row['Groups'])],
                xytext=(5, 5),
                textcoords='offset points',
                ha='left',
                va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.3),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            )
    
    # 除外パスウェイのメッセージ
    if len(pathway_remove) > 0 and pathway_remove_show:
        pathway_text = ', '.join(pathway_remove)
        if len(pathway_text) > 100:
            pathway_text = pathway_text[:97] + "..."
        ax.annotate(
            f"Isolate pathways: {pathway_text}",
            xy=(0.01, 0.99),
            xycoords='axes fraction',
            ha='left',
            va='top',
            fontsize=label_size,
            style='italic'
        )
    
    # 凡例の設定
    if show_legend:
        ax.legend()
    
    # 軸の設定
    if not show_axes:
        ax.axis('off')
    else:
        ax.set_xlabel('Dim 1', fontsize=font_size)
        ax.set_ylabel('Dim 2', fontsize=font_size)
    
    # タイトルの設定
    if title_name:
        ax.set_title(title_name, fontsize=font_size_title)
    
    plt.tight_layout()
    return fig



def plot_embedding_zoom_in(object, slot_name="netP", type="functional", color_use=None,
                          pathway_remove=None, ncol=1, dot_size=(2, 6), 
                          label_size=8, dot_alpha=0.5, do_label=True, 
                          show_legend=False, show_axes=True):
    """
    Zoom into the 2D visualization of the learned manifold learning of the signaling networks
    
    Parameters
    ----------
    object : CellChat object
        CellChat object
    slot_name : str, optional
        The slot name of object
    type : str, optional
        "functional" or "structural"
    color_use : dict or list, optional
        Colors for the cell groups
    pathway_remove : list, optional
        A vector defining the signaling to remove
    ncol : int, optional
        Number of columns if plotting multiple groups
    dot_size : tuple, optional
        A range defining the size of the symbol
    label_size : float, optional
        Font size of the text
    dot_alpha : float, optional
        Transparency
    do_label : bool, optional
        Whether label the each point
    show_legend : bool, optional
        Whether show the legend
    show_axes : bool, optional
        Whether show the axes
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    # 比較設定
    comparison = "single"
    comparison_name = comparison
    
    # データ取得
    try:
        Y = getattr(object, slot_name).similarity[type]["dr"][comparison_name]
        clusters = getattr(object, slot_name).similarity[type]["group"][comparison_name]
        prob = getattr(object, slot_name).prob
    except (AttributeError, KeyError):
        print(f"Required data not found in the object. Check that object.{slot_name}.similarity.{type} exists.")
        return None
    
    # パスウェイ除外処理
    if pathway_remove is None:
        try:
            similarity = getattr(object, slot_name).similarity[type]["matrix"][comparison_name]
            # 孤立したノードを特定（列の和が1）
            pathway_remove = [row for row, col_sum in zip(similarity.index, np.sum(similarity, axis=1).tolist()) if col_sum == 1]
        except (AttributeError, KeyError):
            pathway_remove = []
    
    # パスウェイの除外
    if len(pathway_remove) > 0:
        try:
            pathway_remove_idx = [i for i, path in enumerate(prob.shape[2]) if path in pathway_remove]
            if pathway_remove_idx:
                prob = np.delete(prob, pathway_remove_idx, axis=2)
        except (IndexError, AttributeError):
            print("Error when trying to remove pathways")
    
    # 確率の合計計算
    prob_sum = np.sum(prob, axis=(0, 1))
    
    # データフレーム作成
    pathways = np.array([str(i) for i in range(prob.shape[2])])  # 実際のパスウェイ名がない場合のデモ
    
    df = pd.DataFrame({
        'x': Y[:, 0],
        'y': Y[:, 1],
        'Commun.Prob.': prob_sum/np.max(prob_sum) if np.max(prob_sum) > 0 else prob_sum,
        'labels': pathways,
        'clusters': np.array(clusters, dtype=str)
    })
    
    # カラー設定
    if color_use is None:
        unique_clusters = np.unique(clusters)
        color_use = dict(zip(unique_clusters, ggPalette(len(unique_clusters))))
    
    # 各クラスターでズームインしたプロットを作成
    cluster_ids = np.unique(df['clusters'])
    n_clusters = len(cluster_ids)
    
    # サブプロットのレイアウト計算
    nrow = int(np.ceil(n_clusters / ncol))
    
    fig, axes = plt.subplots(nrow, ncol, figsize=(10 * ncol, 8 * nrow))
    
    # 1つのクラスターの場合にaxesを配列に変換
    if n_clusters == 1:
        axes = np.array([axes])
    
    # 軸の調整
    if nrow > 1 and ncol > 1:
        axes = axes.flatten()
    
    for i, cluster_id in enumerate(cluster_ids):
        if i < len(axes):
            ax = axes[i]
            
            # クラスター固有のデータを抽出
            cluster_df = df[df['clusters'] == cluster_id]
            
            # 散布図をプロット
            scatter = ax.scatter(
                cluster_df['x'], 
                cluster_df['y'],
                s=cluster_df['Commun.Prob.'] * (dot_size[1] - dot_size[0]) + dot_size[0],
                c=[color_use[cluster_id] if isinstance(color_use, dict) else color_use[list(np.unique(df['clusters'])).index(cluster_id)]],
                alpha=dot_alpha
            )
            
            # タイトル設定
            ax.set_title(f"Group {cluster_id}")
            
            # ラベル追加
            if do_label:
                for _, row in cluster_df.iterrows():
                    ax.annotate(
                        row['labels'],
                        (row['x'], row['y']),
                        fontsize=label_size,
                        color=color_use[cluster_id] if isinstance(color_use, dict) else color_use[list(np.unique(df['clusters'])).index(cluster_id)],
                        xytext=(5, 5),
                        textcoords='offset points',
                        ha='left',
                        va='bottom'
                    )
            
            # 凡例の設定
            if show_legend:
                ax.legend()
            
            # 軸の設定
            if not show_axes:
                ax.axis('off')
    
    # 使用されていないサブプロットを削除
    for i in range(n_clusters, len(axes)):
        fig.delaxes(axes[i])
    
    plt.tight_layout()
    return fig

def plot_embedding_pairwise(object, slot_name="netP", type="functional", comparison=None,
                           color_use=None, point_shape=None, pathway_labeled=None, 
                           top_label=1, pathway_remove=None, pathway_remove_show=True,
                           dot_size=(2, 6), label_size=8, dot_alpha=0.5, title_name=None,
                           font_size=10, font_size_title=12, do_label=True, 
                           show_legend=True, show_axes=True):
    """
    2D visualization of the joint manifold learning of signaling networks from two datasets
    
    Parameters
    ----------
    object : CellChat object
        CellChat object
    slot_name : str, optional
        The slot name of object
    type : str, optional
        "functional" or "structural"
    comparison : list, optional
        A numerical vector giving the datasets for comparison
    color_use : dict or list, optional
        Colors for the cell groups
    point_shape : list, optional
        A numeric vector giving the point shapes
    pathway_labeled : list, optional
        A vector giving the signaling names to show
    top_label : float, optional
        The fraction of signaling pathways to label
    pathway_remove : list, optional
        A vector defining the signaling to remove
    pathway_remove_show : bool, optional
        Whether show the removed signaling names
    dot_size : tuple, optional
        A range defining the size of the symbol
    label_size : float, optional
        Font size of the text
    dot_alpha : float, optional
        Transparency
    title_name : str, optional
        The title of plot
    font_size : float, optional
        Font size of the text
    font_size_title : float, optional
        Font size of the title
    do_label : bool, optional
        Whether label the each point
    show_legend : bool, optional
        Whether show the legend
    show_axes : bool, optional
        Whether show the axes
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    # 比較設定
    if comparison is None:
        comparison = list(range(1, len(getattr(object, "meta", {}).get("datasets", [])) + 1))
    
    comparison_name = "-".join(map(str, comparison))
    
    print(f"2D visualization of signaling networks from datasets {', '.join(map(str, comparison))}")
    
    # データ取得
    try:
        Y = getattr(object, slot_name).similarity[type]["dr"][comparison_name]
        clusters = getattr(object, slot_name).similarity[type]["group"][comparison_name]
        
        # 異なるデータセットの名前
        object_names = [name for name in getattr(object, slot_name, {}).__dict__.keys() if name != "similarity"][comparison]
        
        # 確率データの取得
        prob = []
        for i in range(len(comparison)):
            object_net = getattr(getattr(object, slot_name), comparison[i])
            prob.append(object_net.prob)
    except (AttributeError, KeyError, IndexError):
        print(f"Required data not found in the object.")
        return None
    
    # ポイントの形状設定
    if point_shape is None:
        point_shape = ['o', 's', '^', 'D', 'v', '<', '>']
    
    # パスウェイ除外処理
    if pathway_remove is None:
        try:
            similarity = getattr(object, slot_name).similarity[type]["matrix"][comparison_name]
            # 孤立したノードを特定（列の和が1）
            pathway_remove = [row for row, col_sum in zip(similarity.index, np.sum(similarity, axis=1).tolist()) if col_sum == 1]
        except (AttributeError, KeyError):
            pathway_remove = []
    
    # パスウェイの除外とデータ集計
    prob_sum_each = []
    signaling_all = []
    
    for i in range(len(prob)):
        prob_i = prob[i]
        if len(pathway_remove) > 0:
            try:
                pathway_remove_idx = [j for j, path in enumerate(f"{prob_i.shape[2]}--{object_names[i]}") if path in pathway_remove]
                if pathway_remove_idx:
                    prob_i = np.delete(prob_i, pathway_remove_idx, axis=2)
            except (IndexError, AttributeError):
                print("Error when trying to remove pathways")
        
        # 確率の合計を計算
        prob_sum_i = np.sum(prob_i, axis=(0, 1))
        prob_sum_each.append(prob_sum_i)
        
        # シグナリング名とデータセット名を組み合わせ
        signaling_all.extend([f"{j}--{object_names[i]}" for j in range(len(prob_sum_i))])
    
    # すべての確率値を結合
    prob_sum = np.concatenate(prob_sum_each)
    
    # グループとラベルの取得
    group = [name.split("--")[1] for name in signaling_all]
    labels = [name.split("--")[0] for name in signaling_all]
    
    # データフレーム作成
    df = pd.DataFrame({
        'x': Y[:, 0],
        'y': Y[:, 1],
        'Commun.Prob.': prob_sum/np.max(prob_sum) if np.max(prob_sum) > 0 else prob_sum,
        'labels': labels,
        'clusters': np.array(clusters, dtype=str),
        'group': group
    })
    
    # カラー設定
    if color_use is None:
        unique_clusters = np.unique(clusters)
        color_use = dict(zip(unique_clusters, ggPalette(len(unique_clusters))))
    
    # プロット作成
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 各クラスターとグループごとに散布図をプロット
    for g in np.unique(df['group']):
        g_idx = df['group'] == g
        
        for c in np.unique(df['clusters']):
            g_c_idx = g_idx & (df['clusters'] == c)
            g_c_df = df[g_c_idx]
            
            if not g_c_df.empty:
                # 点の形状を選択
                marker = point_shape[list(np.unique(df['group'])).index(g) % len(point_shape)]
                
                # 点の色を選択
                if isinstance(color_use, dict):
                    color = color_use[c]
                else:
                    color = color_use[list(np.unique(df['clusters'])).index(c)]
                
                scatter = ax.scatter(
                    g_c_df['x'], 
                    g_c_df['y'],
                    s=g_c_df['Commun.Prob.'] * (dot_size[1] - dot_size[0]) + dot_size[0],
                    c=[color],
                    marker=marker,
                    alpha=dot_alpha,
                    label=f"{c} - {g}"
                )
    
    # ラベル追加
    if do_label:
        if pathway_labeled is None:
            if top_label < 1:
                # トップNのパスウェイを選択
                top_n = int(top_label * len(df))
                data_label = df.nlargest(top_n, 'Commun.Prob.')
            else:
                data_label = df
        else:
            data_label = df[df['labels'].isin(pathway_labeled)]
        
        for _, row in data_label.iterrows():
            # 点の色を選択
            if isinstance(color_use, dict):
                color = color_use[row['clusters']]
            else:
                color = color_use[list(np.unique(df['clusters'])).index(row['clusters'])]
            
            ax.annotate(
                row['labels'],
                (row['x'], row['y']),
                fontsize=label_size,
                color=color,
                xytext=(5, 5),
                textcoords='offset points',
                ha='left',
                va='bottom',
                bbox=dict(boxstyle='round,pad=0.3', fc='white', alpha=0.3),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0')
            )
    
    # 除外パスウェイのメッセージ
    if len(pathway_remove) > 0 and pathway_remove_show:
        pathway_text = ', '.join(pathway_remove)
        if len(pathway_text) > 100:
            pathway_text = pathway_text[:97] + "..."
        ax.annotate(
            f"Isolate pathways: {pathway_text}",
            xy=(0.01, 0.99),
            xycoords='axes fraction',
            ha='left',
            va='top',
            fontsize=label_size,
            style='italic'
        )
    
    # 凡例の設定
    if show_legend:
        ax.legend()
    
    # 軸の設定
    if not show_axes:
        ax.axis('off')
    else:
        ax.set_xlabel('Dim 1', fontsize=font_size)
        ax.set_ylabel('Dim 2', fontsize=font_size)
    
    # タイトルの設定
    if title_name:
        ax.set_title(title_name, fontsize=font_size_title)
    
    plt.tight_layout()
    return fig

@st.cache_data
def calculate_network_centrality(net):
    """
    ネットワーク中心性を計算する - netAnalysis_signalingRole_networkの計算法を使用
    
    Parameters
    ----------
    net : dict
        ネットワークデータ（strength_matrixを含む）
        
    Returns
    -------
    pd.DataFrame
        中心性計算結果
    """
    import networkx as nx
    import pandas as pd
    import numpy as np
    
    # strength_matrixを取得
    if 'strength_matrix' in net:
        if isinstance(net['strength_matrix'], pd.DataFrame):
            adj = net['strength_matrix'].values
            cell_types = net['strength_matrix'].index.tolist()
        else:
            adj = np.array(net['strength_matrix'])
            # cell_typesを取得
            if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
                cell_types = net['net']['dimnames'][0]
            else:
                cell_types = [f"Cell{i+1}" for i in range(adj.shape[0])]
    elif 'net' in net and 'prob' in net['net']:
        # net['net']['prob']を使用
        adj = np.array(net['net']['prob'])
        if adj.ndim == 3:  # 3次元配列なら2次元に集約
            adj = np.sum(adj, axis=2)
        
        if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            cell_types = net['net']['dimnames'][0]
        else:
            cell_types = [f"Cell{i+1}" for i in range(adj.shape[0])]
    else:
        print("有効な相互作用データが見つかりません")
        return pd.DataFrame()
    
    # グラフを作成
    G = nx.from_numpy_array(adj, create_using=nx.DiGraph)
    mapping = {i: cell_types[i] for i in range(len(cell_types))}
    G = nx.relabel_nodes(G, mapping)
    
    # 中心性指標の計算
    # 1. Sender & Receiver: 元のグラフで計算
    sender = dict(G.out_degree(weight='weight'))
    receiver = dict(G.in_degree(weight='weight'))
    
    # 2. Betweenness: 重み逆数変換したグラフで計算
    G_betweenness = G.copy()
    for u, v, d in G_betweenness.edges(data=True):
        if 'weight' in d and d['weight'] > 0:
            d['weight'] = 1.0 / d['weight']
    
    try:
        betweenness = nx.betweenness_centrality(G_betweenness, weight='weight')
    except:
        print("媒介中心性の計算に失敗しました")
        betweenness = {node: 0.0 for node in G.nodes()}
    
    # 3. Mediator & Influencer: 特殊関数で計算
    try:
        from pages.cellchat_vis import compute_flow_betweenness, compute_information_centrality
        mediator = compute_flow_betweenness(G)
        influencer = compute_information_centrality(G)
    except Exception as e:
        print(f"特殊中心性計算エラー: {str(e)}")
        # フォールバック: 標準的な中心性指標で代用
        mediator = {node: 0.0 for node in G.nodes()}
        influencer = {node: 0.0 for node in G.nodes()}
    
    # データフレームを作成
    cell_type_list = list(G.nodes())
    centrality_df = pd.DataFrame({
        'cell_type': cell_type_list,
        'degree': [sender[node] + receiver[node] for node in cell_type_list],
        'in_degree': [receiver[node] for node in cell_type_list],
        'out_degree': [sender[node] for node in cell_type_list],
        'betweenness': [betweenness[node] for node in cell_type_list],
        'flowbet': [mediator[node] for node in cell_type_list],  # Mediator
        'info': [influencer[node] for node in cell_type_list]    # Influencer
    })
    
    return centrality_df


# メイン処理
if __name__ == "__main__":
    st.set_page_config(page_title="CellChat-Python", page_icon="💬")
    st.sidebar.title("Options")
    if "cellchat_res" not in st.session_state:
        st.session_state.cellchat_res = None

    if "cellchat_temp_dir" not in st.session_state:
        st.session_state.cellchat_temp_dir = None
        cellchat_temp_dir = "temp/" + str(round(time.time()))
        if not os.path.exists('temp'):
            os.mkdir('temp')
        else:
            clear_old_directories("temp")
            clear_old_files("temp")
        os.mkdir(cellchat_temp_dir)
        st.session_state.cellchat_temp_dir = cellchat_temp_dir
    else:
        cellchat_temp_dir = st.session_state.cellchat_temp_dir

    st.markdown("### CellChat Python implemetation")

    uploaded_file = st.file_uploader("H5AD or cellchat result pkl file", type=['h5ad', 'pkl'], help="H5ADファイルで新規解析、PKLファイルで保存済み結果を読み込みます")

    if uploaded_file is not None:
        file_type = uploaded_file.name.split('.')[-1].lower()
        
        # PKLファイルの場合は保存済みCellChat結果を読み込む
        if file_type == 'pkl':
            try:
                # 保存済みファイルを読み込む
                bytes_data = uploaded_file.read()
                result = pickle.loads(bytes_data)
                st.session_state.cellchat_res = result
                cell_list = []
                if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 0:
                    cell_list = sorted(result['net']['dimnames'][0])
                print("cell_list")
                print(cell_list)
                if 'groupby' in result:
                    # グローバル変数のgroupbyを更新
                    groupby = result['groupby']
                st.success("CellChat結果を正常に読み込みました！可視化を開始できます")
                
                # アップロードされたファイル名を保存（保存ボタン用）
                if 'uploaded_filename' not in st.session_state:
                    st.session_state.uploaded_filename = uploaded_file.name
                
                # ここからは可視化タブに直接進む
                if st.session_state.cellchat_res is not None:
                    # cell_listの抽出 - 複数の方法を試みる
                    pass
                    
            except Exception as e:
                st.error(f"PKLファイルの読み込みに失敗しました: {str(e)}")
                st.error(traceback.format_exc())
        
        # H5ADファイルの場合は通常の解析フロー
        elif file_type == 'h5ad':
            current_file_id = f"{uploaded_file.name}_{uploaded_file.size}"

            if 'last_file_id' not in st.session_state:
                st.session_state.last_file_id = None

            if current_file_id != st.session_state.last_file_id:
                st.cache_data.clear()
                st.session_state.last_file_id = current_file_id
                st.success("新しいファイルが検出されました。キャッシュをクリアしました。")

            # アップロードされたファイル名を保存（保存ボタン用）
            st.session_state.uploaded_filename = uploaded_file.name

            adata = read_h5ad(uploaded_file)

            st.write("Uploaded data:")
            st.write(adata)
            temp_df = pd.DataFrame(
            adata.X[:5,:8].toarray() if scipy.sparse.issparse(adata.X) else adata.X[:5,:8],
            index=adata.obs_names[:5],
            columns=adata.var_names[:8]
            )
            st.write("adata.X (default data layer)")
            st.dataframe(temp_df) 


            meta = adata.obs.columns.to_list()

            for i in ['nFeature_RNA', 'nCount_RNA', 'percent.mt', 'Cell_id']:
                try:
                    meta.remove(i)
                except:
                    pass

            # 利用可能なレイヤーを表示
            available_layers = list(adata.layers.keys()) if hasattr(adata, 'layers') else []
            st.write(f"Additional data layers: {', '.join(available_layers) if available_layers else 'NA'}")

            groupby = st.selectbox("Cell type:", meta, 
                                 index=find_first_index_or_zero(meta, ["cell.ident", "seurat_clusters", "cell_type", "louvain"]))
            
            cell_list = sorted(adata.obs[groupby].cat.categories.to_list() if hasattr(adata.obs[groupby], 'cat') else sorted(adata.obs[groupby].unique()))
            st.write(", ".join(cell_list))

            with st.form("Basic settings:"):
                species = st.radio("Species:", ('mouse', 'human'), 
                                  index=check_species_index(adata.var.index.to_list()[:50]))
                
                data_layer = st.selectbox("Using layer:", 
                                         ['X (default)'] + available_layers,
                                         index=0)
                
                selected_types = st.multiselect(
                    "Signaling types to analyze:",
                    ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling"],
                    default=["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"]
                )     
                expr_prop = st.number_input('Min fraction of expressing cells:', 
                                           min_value=0.0, max_value=0.9, step=0.01, value=0.1)

                thresh_p = st.number_input('Threshold p value for overexpressed genes:', 
                                           min_value=0.0, max_value=0.9, step=0.01, value=0.05,
                                           help = "Wilcoxonによる変化遺伝子選択のthreshold")
                
           #     remove_complex = st.checkbox("リガンド/レセプター複合体を除外する (例: L17A_IL17F, IL17RA_IL17RC)?")
                
          #      n_perms = st.slider("置換検定回数:", 
           #                        min_value=20, max_value=500, step=20, value=100)
                n_perms = 100
                
                n_cpus = st.slider("Num CPUs:", 
                                  min_value=1, max_value=os.cpu_count(), step=1, value=os.cpu_count()-2)
                
                submitted_basic = st.form_submit_button("Apply settings")

            if submitted_basic:
                st.session_state.cellchat_res = None

            st.markdown("### Click [Apply settings] before run!")

            if st.button('Run calc') or (st.session_state.cellchat_res is not None):
                if st.session_state.cellchat_res is None:
                    # CellChatDBを取得
                    try:
                        start_cpu = time.time()
                        with st.spinner('Getting CellChatDB from R...'):
                            cellchatdb = get_cellchatdb_from_r(species=species)
                            # デバッグ出力を追加
                        cellchatdb = debug_cellchatdb(cellchatdb)
                        st.success(f"{species} CellChatDB is succcesfully obtained")
                    except Exception as e:
                        st.error(f"CellChatDBの取得中にエラーが発生しました: {str(e)}")
                        st.stop()
                    interaction = cellchatdb['interaction']
                    interaction_filtered = interaction[interaction['annotation'].isin(selected_types)]
                    cellchatdb['interaction'] = interaction_filtered
                    
                    st.write(f"LR pair number: {len(cellchatdb['interaction'])}")
                    st.write(cellchatdb['interaction'].head(5))
                    
                    # データレイヤーの設定
                    use_layer = None if data_layer == 'X (default)' else data_layer

                    with st.spinner('CellChat calculation...'):
                        result = cellchat_analysis_optimized(
                            adata,
                            groupby=groupby,
                            db=cellchatdb,  # CellChatDBを渡す
                            use_layer=use_layer,
                            min_cells=5,
                            expr_prop=expr_prop,
                         #   log_scale=log_transformed,
                            pseudocount=1.0,
                            k=0.5,
                            trim_threshold=0.05,  
                            nboot=n_perms,
                            seed=12345,
                            n_jobs=n_cpus,
                         #   remove_complex=remove_complex,
                            trim=0.1,
                            apply_pval_filter=False,
                            use_gpu=False,                # GPUを有効化
                            gpu_precision='float64',     # 倍精度計算を使用
                            validate_gpu_results=False,   # 結果検証を有効化
                            hybrid_mode=False           # ハイブリッドモードを使用
                        )
                        
                        st.session_state.cellchat_res = result
                        time_cpu = time.time() - start_cpu
                        st.write(f"Computing time: {round(time_cpu)}")
                    
                    if 'error' in result:
                        st.error(f"解析中にエラーが発生しました: {result['error']}")
                        st.code(result['traceback'])
                    else:
                        st.success('CellChat解析が完了しました！')
                
            # 以下は結果の表示部分
            
        result = st.session_state.cellchat_res

        if result is not None:
            
            if 'error' in result:
                st.error(f"解析中にエラーが発生しました: {result['error']}")
                if 'traceback' in result:
                    st.code(result['traceback'])
            else:
                # 結果表示
                st.subheader("1. CellChat results")
                
                st.write("相互作用テーブル (上位10件):")
                if 'results' in result and len(result['results']) > 0:
                    st.write(result['results'].sort_values('prob', ascending=False).head(10))
                    
                    # 全結果のダウンロード
                    csv = result['results'].to_csv(index=False)
                    st.download_button(
                        label="結果をCSVとしてダウンロード",
                        data=csv,
                        file_name='cellchat_results.csv',
                        mime='text/csv',
                    )
                else:
                    st.warning("有効な相互作用結果がありません。パラメータを調整して再試行してください。")
                    
                # デバッグ情報を表示
                with st.expander("デバッグ情報"):
                    if 'network' in result and 'strength_matrix' in result['network']:
                        debug_interaction_matrix(result['network']['strength_matrix'], "相互作用強度行列")
                    if 'network' in result and 'count_matrix' in result['network']:
                        debug_interaction_matrix(result['network']['count_matrix'], "相互作用数行列")
                
                # 可視化
                st.subheader("2. Visualization of communication network")
                plt.style.use('default')
                with st.sidebar:
                    colormap_options = [
                        "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "RdPu", "Purples", 
                        "PuRd", "Blues", "Greens", "YlGn", "YlGnBu", "GnBu", "Greys", "binary",
                        "viridis", "plasma", "inferno", "magma", "cividis", 
                        "viridis_r", "plasma_r", "inferno_r", "magma_r", "cividis_r"
                    ]
                    color_heatmap = st.selectbox(
                        "Heatmap color:",
                        colormap_options,
                        index=2
                    )
                    cmap_name = st.sidebar.selectbox(
                        "Select colormap for other plots:",
                        ["tab10","Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
                        "tab20b", "tab20c","Pastel1",
                        "Pastel2",  "Accent", "viridis", "plasma", "inferno", "magma"],
                        index=0
                    )
                    sort_cell = st.checkbox("Change cell order?")
                    if sort_cell:
                        with st.form("Sorter"):
                            sorted_order = sort_items(sorted(adata.obs[groupby].unique()))
                            submitted_sort = st.form_submit_button("Done sorting")
                    else:
                        sorted_order = None
                
                tabs = st.tabs([
                    "Heatmap", 
                    "Chord", 
                    "LR contribution", 
                    "Dot", 
                    "Centrality", 
                    "Role Scatter",
                    "Signal contribution",
                    "Circle",   # New tab for circle plots
                    "Role",   # New tab for signaling role analysis 
                    "Expression"     # New tab for gene expression analysis
                ])

                # Keep your existing tab implementations
                with tabs[0]:
                    st.markdown("#### Cell interaction heatmap")
                    
                    heatmap_type = st.radio(
                        "Heatmap type:",
                        ["Interaction strength", "Interaction number"],
                        horizontal=True
                    )
                    heatmap_annot = st.checkbox("Show values?", value=True)
                    
                    
                    if 'network' in result:
                        if heatmap_type == "Interaction strength":
                            matrix = result['network']['strength_matrix']
                            title = "Cell-Cell Interaction Strength"
                        else:
                            matrix = result['network']['count_matrix']
                            title = "Number of Significant Interactions"

                        # Use the enhanced heatmap visualization
                        measure = "weight" if heatmap_type == "Interaction strength" else "count"
                        fig_heatmap = netVisual_heatmap(
                            result, 
                            measure=measure, 
                            color_heatmap=color_heatmap, 
                            font_size=14, 
                            font_size_title=20,
                            sorted_order = sorted_order,
                            annot=heatmap_annot
                        )   
                        st.pyplot(fig_heatmap)
                        
                        # PDF saving and download unchanged
                        pdf_path = f"{cellchat_temp_dir}/heatmap.pdf"
                        fig_heatmap.savefig(pdf_path, bbox_inches='tight')
                        
                        with open(pdf_path, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                            st.download_button(
                                label="Download heatmap pdf",
                                data=pdf_bytes,
                                file_name=f'cellchat_{heatmap_type}_heatmap.pdf',
                                mime='application/octet-stream'
                            )
                    else:
                        st.warning("ネットワーク情報が利用できません")

                with tabs[1]:  # コードダイアグラム tab
                    st.markdown("#### Chord diagram")
                    
                    if 'netP' in result and 'pathways' in result['netP']:
                        # 選択肢: "Aggregate" または 個別パスウェイの複数選択
                        option_type = st.radio(
                            "Display type:",
                            ["Aggregate", "Specific pathways"],
                            horizontal=True
                        )
                        
                        if option_type == "Aggregate":
                            selected_pathway = "Aggregate"
                            # Interaction measure を選択するラジオボタン
                            diagram_measure = st.radio(
                                "Select value:",
                                ["Interaction weight", "Interaction number"],
                                horizontal=True
                            )
                            measure = "weight" if diagram_measure == "Interaction weight" else "count"
                        else:
                            # 複数のパスウェイを選択可能にする
                            pathways = list(result['netP']['pathways'])
                            selected_pathways = st.multiselect(
                                "Select pathways (multiple selection allowed):",
                                sorted(pathways),
                                default=[pathways[0]] if pathways else []
                            )
                            selected_pathway = selected_pathways  # リストとして渡す
                            measure = "weight"  # 特定パスウェイでは常にweight
                        
                        if st.button("Generate chord diagram"):
                            if option_type == "Specific pathways" and not selected_pathways:
                                st.warning("Please select at least one pathway.")
                            else:
                                with st.spinner("Making chord diagram..."):
                                    from pages.cellchat_vis import create_chord_diagram_pycirclize
                                    try:
                                        fig_chord = create_chord_diagram_pycirclize(
                                            net=result,
                                            pathway_name=selected_pathway,
                                            figsize=(10, 10),
                                            cmap_name=cmap_name,
                                            sorted_order=sorted_order,
                                            measure=measure
                                        )
                                        st.pyplot(fig_chord)
                                        
                                        # PDF として保存・ダウンロード
                                        import io
                                        pdf_buffer = io.BytesIO()
                                        fig_chord.savefig(pdf_buffer, format="pdf", bbox_inches="tight")
                                        pdf_data = pdf_buffer.getvalue()
                                        
                                        if option_type == "Aggregate":
                                            filename = 'cellchat_chord_aggregate.pdf'
                                        else:
                                            # 複数パスウェイの場合はパスウェイ数を表示
                                            if len(selected_pathway) <= 2:
                                                pathway_str = '_'.join(selected_pathway)
                                            else:
                                                pathway_str = f'{len(selected_pathway)}_pathways'
                                            filename = f'cellchat_chord_{pathway_str}.pdf'
                                        
                                        st.download_button(
                                            label="Download chord diagram",
                                            data=pdf_data,
                                            file_name=filename,
                                            mime='application/octet-stream'
                                        )
                                    except Exception as e:
                                        st.error(f"コードダイアグラム生成エラー: {str(e)}")
                                        import traceback
                                        st.error(traceback.format_exc())
                    else:
                        st.warning("結果にパスウェイ情報がありません。")
                
                
                with tabs[2]:
                    st.markdown("#### リガンド-レセプターペアの寄与度")
                    
                    st.info("""
                    このグラフは各リガンド-レセプターペアが全体のシグナリングネットワークにどの程度寄与しているかを示しています。
                    寄与度が高いペアほど、細胞間コミュニケーションにおいて重要な役割を果たしています。
                    可能な場合は、各ペアが属するシグナリングパスウェイも表示されます。
                    """)
                    
                    if 'network' in result and 'lr_contribution' in result['network']:
                        lr_contribution = result['network']['lr_contribution']
                        
                        if isinstance(lr_contribution, np.ndarray) and len(lr_contribution) > 0:
                            # DataFrameの情報からデータの長さを取得
                            lr_contrib_len = np.sum(lr_contribution > 0)
                            
                            # 表示するトップLRペア数の選択（最大値を調整）
                            top_n = st.slider("表示するトップLRペア数:", 
                                             min_value=1, 
                                             max_value=min(30, int(lr_contrib_len)) if lr_contrib_len > 0 else 10, 
                                             value=min(10, int(lr_contrib_len)) if lr_contrib_len > 0 else 10)
                            
                            # カラースキーム選択
                            color_scheme = st.selectbox(
                                "カラースキーム:", 
                                ["viridis", "plasma", "inferno", "magma", "cividis", "Blues", "Greens", "Oranges", "Reds"],
                                index=0
                            )
                            
                            # プロットの高さを動的に調整（ペア数に比例）
                            plot_height = 6 + 0.3 * top_n  # ペア数に応じて高さを調整
                            from pages.cellchat_vis import  plot_lr_contribution

                            # ネットワークサマリーにLR名情報を追加
                            network_data = result['network'].copy()
                            
                            # interaction_namesの取得
                            interaction_names = result['net']['dimnames'][2]
                            
                            # LRペア情報の辞書を作成
                            lr_pairs = {}
                            for i, interaction_name in enumerate(interaction_names):
                                # resultsデータフレームからLR情報を検索
                                matching_rows = result['results'][result['results']['interaction_name'] == interaction_name]
                                if not matching_rows.empty:
                                    row = matching_rows.iloc[0]
                                    lr_pairs[interaction_name] = f"{row['ligand']}-{row['receptor']}"
                            
                            # network_dataにLR情報を追加
                            network_data['lr_names'] = lr_pairs
                            network_data['interaction_names'] = interaction_names
                            
                            # LR寄与度のグラフを生成
                            fig_lr = plot_lr_contribution(
                                network_data, 
                                top_n=top_n, 
                                figsize=(12, plot_height)
                            )

                            
                            # 図のスタイルを更新
                            for ax in fig_lr.get_axes():
                                # カラーマップを変更
                                if hasattr(ax, 'patches') and len(ax.patches) > 0:
                                    cmap = plt.cm.get_cmap(color_scheme)
                                    colors = cmap(np.linspace(0.2, 0.9, len(ax.patches)))
                                    for i, patch in enumerate(ax.patches):
                                        patch.set_color(colors[i])
                                
                                # 背景色の設定
                             #   ax.set_facecolor('#f8f9fa')  # 薄いグレー
                                
                                # グリッド線の設定
                             #   ax.grid(axis='x', linestyle='--', alpha=0.3)
                                
                                # 枠線の設定
                                for spine in ax.spines.values():
                                    spine.set_color('#cccccc')
                            
                            st.pyplot(fig_lr)
                            
                            # PDF保存とダウンロード
                            pdf_path = f"{cellchat_temp_dir}/lr_contribution.pdf"
                            fig_lr.savefig(pdf_path, bbox_inches='tight')
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label="LR寄与度PDFをダウンロード",
                                    data=pdf_bytes,
                                    file_name='cellchat_lr_contribution.pdf',
                                    mime='application/octet-stream'
                                )
                            
                            # 寄与度上位のLRペア情報をテーブルでも表示（オプション）
                            if st.checkbox("寄与度上位のLRペア詳細を表示"):
                                try:
                                    # 寄与度データの取得・整形
                                    lr_data = []
                                    for i, val in enumerate(lr_contribution):
                                        if val > 0:
                                            # LRペア情報の取得試行
                                            lr_info = {"index": i, "contribution": val}
                                            
                                            # 結果データフレームからの情報取得
                                            if 'results' in result and not result['results'].empty:
                                                # インタラクション名の取得試行
                                                if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 2:
                                                    interaction_names = result['net']['dimnames'][2]
                                                    if i < len(interaction_names):
                                                        matching_rows = result['results'][
                                                            result['results']['interaction_name'] == interaction_names[i]
                                                        ]
                                                        if not matching_rows.empty:
                                                            row = matching_rows.iloc[0]
                                                            lr_info.update({
                                                                "interaction_name": row.get('interaction_name', f"LR_{i}"),
                                                                "ligand": row.get('ligand', ''),
                                                                "receptor": row.get('receptor', ''),
                                                                "pathway": row.get('pathway_name', '')
                                                            })
                                            
                                            # 情報が取得できなかった場合のフォールバック
                                            if 'interaction_name' not in lr_info:
                                                lr_info.update({
                                                    "interaction_name": f"LR_{i}",
                                                    "ligand": "不明",
                                                    "receptor": "不明",
                                                    "pathway": "不明"
                                                })
                                            
                                            lr_data.append(lr_info)
                                    
                                    # 寄与度で降順ソート
                                    lr_df = pd.DataFrame(lr_data).sort_values('contribution', ascending=False).head(top_n)
                                    
                                    # 列名を日本語に変更
                                    lr_df.columns = [
                                        'インデックス' if col == 'index' else
                                        '寄与度' if col == 'contribution' else
                                        'インタラクション名' if col == 'interaction_name' else
                                        'リガンド' if col == 'ligand' else
                                        'レセプター' if col == 'receptor' else
                                        'パスウェイ' if col == 'pathway' else col
                                        for col in lr_df.columns
                                    ]
                                    
                                    # 表示列の順序設定
                                    display_cols = ['インタラクション名', 'リガンド', 'レセプター', 'パスウェイ', '寄与度']
                                    display_cols = [col for col in display_cols if col in lr_df.columns]
                                    
                                    # テーブル表示
                                    st.dataframe(lr_df[display_cols])
                                
                                except Exception as e:
                                    st.error(f"LRペア詳細情報の取得エラー: {str(e)}")
                        else:
                            st.warning("LR寄与度データがありません")
                    else:
                        st.warning("LR寄与度情報が利用できません")
                
                with tabs[3]:
                    st.markdown("#### LR相互作用ドットプロット")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        source_cells = st.multiselect(
                            'リガンド発現細胞の選択:',
                            cell_list,
                            default=[cell_list[0]] if len(cell_list) > 0 else []
                        )
                    
                    with col2:
                        target_cells = st.multiselect(
                            'レセプター発現細胞の選択:',
                            cell_list,
                            default=[cell_list[0]] if len(cell_list) > 0 else []
                        )
                    
                    if 'results' in result:
                        pval_col = 'pval' if 'pval' in result['results'].columns else 'pvalue'
                        if pval_col in result['results'].columns:
                            sig_lr_count = len(result['results'][result['results'][pval_col] <= 0.05])
                        else:
                            sig_lr_count = len(result['results'])
                            
                        top_n_dot = st.slider("表示するトップLRペア数 (ドットプロット):", 
                                             min_value=1, 
                                             max_value=min(30, sig_lr_count) if sig_lr_count > 0 else 15, 
                                             value=min(15, sig_lr_count) if sig_lr_count > 0 else 15)
                        pval_threshold = st.slider("P値閾値:", min_value=0.001, max_value=0.1, value=0.05, step=0.001)
                        
                        if source_cells and target_cells:
                            fig_dot = plot_dot_lr_network(
                                result['results'],
                                source_cells,
                                target_cells,
                                top_n=top_n_dot,
                                pval_threshold=pval_threshold
                            )
                            
                            st.pyplot(fig_dot)
                            
                            # PDF保存とダウンロード
                            pdf_path = f"{cellchat_temp_dir}/dot_plot.pdf"
                            fig_dot.savefig(pdf_path, bbox_inches='tight')
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label="ドットプロットPDFをダウンロード",
                                    data=pdf_bytes,
                                    file_name='cellchat_dot_plot.pdf',
                                    mime='application/octet-stream'
                                )
                        else:
                            st.warning("送信側と受信側の細胞タイプを選択してください。")
                    else:
                        st.warning("結果情報が利用できません")
                
                with tabs[4]:
                    st.markdown("#### ネットワーク中心性")
                    
                    if 'network' in result and 'network_centrality' in result['network']:
                        from pages.cellchat_vis import plot_network_centrality
    #                    fig_centrality = plot_network_centrality(result['network'])

                      
                      #  if 'network_centrality' in result['network']:
                      #      fig_centrality = plot_network_centrality({'network_centrality': result['network']['network_centrality']})
                      #  else:
                      #      fig_centrality = plot_network_centrality({'network_centrality': pd.DataFrame()})
                      #  st.pyplot(fig_centrality)
                        
                        # PDF保存とダウンロード
                      #  pdf_path = f"{cellchat_temp_dir}/centrality.pdf"
                      #  fig_centrality.savefig(pdf_path, bbox_inches='tight')
                        
                      #  with open(pdf_path, "rb") as pdf_file:
                      #      pdf_bytes = pdf_file.read()
                      #      st.download_button(
                      #          label="中心性測定PDFをダウンロード",
                      #          data=pdf_bytes,
                      #          file_name='cellchat_centrality.pdf',
                      #          mime='application/octet-stream'
                      #      )
                        
                        st.subheader("In/Out activity")
                        
                        if 'outgoing' in result['network'] and 'incoming' in result['network']:
                            outgoing = result['network']['outgoing']
                            incoming = result['network']['incoming']
                            
                            # 送受信データをマージ
                            if isinstance(outgoing, pd.Series) and not outgoing.empty:
                                io_df = pd.DataFrame({
                                    'cell_type': outgoing.index,
                                    'outgoing': outgoing.values,
                                    'incoming': incoming.values
                                })
                                
                                st.write(io_df)
                                
                                fig, ax = plt.subplots(figsize=(10, 6))
                                
                                x = np.arange(len(io_df))
                                width = 0.35
                                
                                ax.bar(x - width/2, io_df['outgoing'], width, label='Outgoing')
                                ax.bar(x + width/2, io_df['incoming'], width, label='Incoming')
                                
                                ax.set_xticks(x)
                                ax.set_xticklabels(io_df['cell_type'], rotation=45, ha='right')
                                ax.legend()
                                
                                ax.set_title('Outgoing vs Incoming Communication Activity')
                                ax.set_ylabel('Interaction Strength')
                                
                                plt.tight_layout()
                                
                                st.pyplot(fig)


                                # PDF保存とダウンロード
                                pdf_path = f"{cellchat_temp_dir}/centrality.pdf"
                                fig.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="Download in/out PDF",
                                        data=pdf_bytes,
                                        file_name='cellchat_in_out_activity.pdf',
                                        mime='application/octet-stream'
                                    )

                            else:
                                st.warning("送受信活性データが空です")
                        else:
                            st.warning("送受信活性情報が利用できません")
                    else:
                        st.warning("ネットワーク中心性情報が利用できません")
                
                with tabs[5]:
                    st.markdown("#### シグナル役割分析")
                    
                    st.info("""
                    シグナル役割分析は、通信ネットワークにおける各細胞タイプの
                    送信者または受信者としての貢献度と全体的な影響力を示します。
                    """)
                    
                    # UIを2カラムに分割して配置
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        # パスウェイ選択（Aggregate または特定パスウェイ）
                        if 'netP' in result and 'pathways' in result['netP']:
                            pathways = result['netP']['pathways']
                            selected_path = st.selectbox(
                                "役割分析のパスウェイを選択:",
                                options=["Aggregate"] + sorted(list(pathways)),
                                index=0
                            )
                        else:
                            st.warning("結果にパスウェイ情報がありません。")
                            selected_path = "Aggregate"
                        
                        # X軸とY軸の中心性指標を選択
                        measure_options = ["outdeg", "indeg", "outdeg_unweighted", "indeg_unweighted"]
                        if 'netP' in result and 'centr' in result['netP'] and len(result['netP']['centr']) > 0:
                            first_key = 0 if isinstance(result['netP']['centr'], list) else list(result['netP']['centr'].keys())[0]
                            if isinstance(result['netP']['centr'][first_key], dict):
                                measure_options = list(result['netP']['centr'][first_key].keys())
                        
                        x_measure = st.selectbox("X axis:", options=measure_options, 
                                                index=measure_options.index("outdeg") if "outdeg" in measure_options else 0)
                        y_measure = st.selectbox("Y axis:", options=measure_options, 
                                                index=measure_options.index("indeg") if "indeg" in measure_options else 0)
                    
                    with col2:
                        # 表示カスタマイズオプション
                        do_label = st.checkbox("細胞タイプラベルを表示", value=True)
                        dot_alpha = st.slider("点の透明度:", min_value=0.1, max_value=1.0, value=0.7, step=0.1)
                        
                   #     # 細胞グループの入力（オプション）
                   #     cell_group_input = st.text_input(
                   #         "細胞グループ (オプション, 例: 'CellA:Group1,CellB:Group2'):",
                   #         value=""
                   #     )
                        
                        # 入力された細胞グループを解析
                        cell_groups = None
                   #     if cell_group_input.strip():
                   #         try:
                   #             cell_groups = dict(item.strip().split(":") for item in cell_group_input.split(","))
                   #             st.success(f"{len(cell_groups)}個の細胞グループが定義されました")
                   #         except:
                   #             st.warning("細胞グループの形式が正しくありません。'CellType:GroupName'形式でカンマ区切りで入力してください。")
                    
                    # 詳細設定用のエクスパンダー
                    with st.expander("Options"):
                        label_size = st.slider("ラベルサイズ:", min_value=4, max_value=16, value=8)
                        dot_size_min = st.slider("最小点サイズ:", min_value=1, max_value=10, value=2)
                        dot_size_max = st.slider("最大点サイズ:", min_value=dot_size_min+1, max_value=20, value=6)
                        
                        # Weight Min/Max設定
                        weight_minmax_input = st.text_input(
                            "Weight MinMax (最小値,最大値):",
                            value=""
                        )
                        weight_MinMax = None
                        if weight_minmax_input.strip():
                            try:
                                min_val, max_val = map(float, weight_minmax_input.split(","))
                                weight_MinMax = (min_val, max_val)
                            except:
                                st.warning("Weight MinMaxの形式が正しくありません。'最小値,最大値'の形式で入力してください。")
                        
                        # タイトルと軸ラベル
                        title = st.text_input("タイトル:", value="Signaling Role Analysis")
                        xlabel = st.text_input("X軸ラベル:", value="Outgoing interaction strength" if x_measure == "outdeg" else x_measure)
                        ylabel = st.text_input("Y軸ラベル:", value="Incoming interaction strength" if y_measure == "indeg" else y_measure)
                    
                    # 分析実行ボタン
                    if st.button("Generate signaling role scatter plot"):
                        with st.spinner("Calculating..."):
                            try:
                                from pages.cellchat_vis import netAnalysis_signalingRole_scatter
                                
                                # パスウェイパラメータ設定
                                pathway_param = None if selected_path == "Aggregate" else [selected_path]
                                
                                # 散布図を生成
                                fig_role = netAnalysis_signalingRole_scatter(
                                    net=result,
                                    signaling=pathway_param,
                                    color_use=cmap_name,
                                    slot_name="netP",
                                    group=cell_groups,
                                    weight_MinMax=weight_MinMax,
                                    dot_size=(dot_size_min, dot_size_max),
                                    label_size=label_size,
                                    dot_alpha=dot_alpha,
                                    x_measure=x_measure,
                                    y_measure=y_measure,
                                    xlabel=xlabel,
                                    ylabel=ylabel,
                                    title=title,
                                    do_label=do_label
                                )
                                
                                # 図を表示
                                st.pyplot(fig_role)
                                
                                # PDFの保存とダウンロード提供
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_{selected_path}.pdf"
                                fig_role.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="シグナル役割散布図をPDFでダウンロード",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_{selected_path}.pdf',
                                        mime='application/octet-stream'
                                    )
                            except Exception as e:
                                st.error(f"シグナル役割散布図エラー: {str(e)}")
                                st.error(traceback.format_exc())  # トレースバックを表示してデバッグを容易に
                            

                
                with tabs[6]:
                    st.markdown("#### シグナル貢献度分析")
                    
                    st.info("""
                    このヒートマップでは、特定の細胞グループの送信または受信シグナルに最も寄与するシグナル経路を示します。
                    カラーバーは、細胞グループ間のシグナリング経路の相対的な強度を表します（注: 値は行方向でスケーリングされています）。
                    上部のカラーバープロットは、ヒートマップに表示されるすべてのシグナリング経路を合計した細胞グループの総シグナリング強度を示します。
                    右側のグレーバープロットは、ヒートマップに表示されるすべての細胞グループを合計したシグナリング経路の総シグナリング強度を示します。
                    """)
                    from pages.cellchat_vis import netAnalysis_signalingRole_heatmap
                    
                    # パターン選択（送信/受信）
                    pattern = st.radio(
                        "シグナルパターン:",
                        ["outgoing", "incoming"],
                        horizontal=True,
                        help="outgoing: 細胞からの送信シグナル, incoming: 細胞への受信シグナル"
                    )
                    
                    # シグナル経路の選択
                    if 'netP' in result and 'pathways' in result['netP']:
                        pathway_options = result['netP']['pathways']
                        
                        selection_mode = st.radio(
                            "シグナル経路選択:",
                            ["すべて", "選択"],
                            horizontal=True
                        )
                        
                        if selection_mode == "選択":
                            selected_pathways = st.multiselect(
                                "分析する経路:",
                                options=sorted(pathway_options),
                                default=sorted(pathway_options)[:min(5, len(pathway_options))],
                                help="分析する特定のシグナル経路を選択してください。"
                            )
                            signaling_param = selected_pathways if selected_pathways else None
                        else:
                            signaling_param = None  # すべての経路を分析
                            st.info(f"すべての経路を分析します ({len(pathway_options)}個の経路)")
                    else:
                        st.warning("結果にパスウェイ情報がありません。")
                        signaling_param = None
                    
                    # 可視化設定
                    with st.expander("可視化設定"):
                        col1, col2 = st.columns(2)
                        
                        with col1:

                            
                            display_numbers = st.checkbox(
                                "Show values in heatmap",
                                value=False,
                                help="ヒートマップのセルに数値を表示します"
                            )
                        
                        with col2:
                            cluster_rows = st.checkbox(
                                "行をクラスタリング",
                                value=True,
                                help="階層的クラスタリングを使用して行を整理します"
                            )
                            
                            cluster_cols = st.checkbox(
                                "列をクラスタリング",
                                value=False,
                                help="階層的クラスタリングを使用して列を整理します"
                            )
                        
                        thresh = st.slider(
                            "有意性閾値:",
                            min_value=0.01,
                            max_value=0.5,
                            value=0.05,
                            step=0.01,
                            help="有意とみなす相互作用のp値閾値"
                        )
                        
                        title = st.text_input(
                            "カスタムタイトル:",
                            value="",
                            help="図のカスタムタイトル（空白の場合は自動生成）"
                        )
                        
                        width = st.slider("図の幅:", min_value=6, max_value=15, value=10, step=1)
                        height = st.slider("図の高さ:", min_value=4, max_value=12, value=8, step=1)
                    
                    # ヒートマップ生成ボタン
                    if st.button("シグナル貢献度ヒートマップを生成"):
                        with st.spinner("ヒートマップを生成中..."):
                            try:
                                # ヒートマップ生成
                                # Generate the heatmap
                                fig_heatmap = netAnalysis_signalingRole_heatmap(
                                    net=result,
                                    signaling=signaling_param,
                                    pattern=pattern,
                                    title=title if title else None,
                                    color_heatmap=color_heatmap,
                                    thresh=thresh,
                                    width=width,
                                    height=height,
                                    font_size=10,
                                    cluster_rows=cluster_rows,
                                    cluster_cols=cluster_cols,
                                    display_numbers=display_numbers,
                                    sorted_order=sorted_order
                                )
                                
                                # 図を表示
                                st.pyplot(fig_heatmap)
                                
                                # PDFとして保存・ダウンロード提供
                                pathway_str = "all_pathways" if signaling_param is None else "_".join(signaling_param[:3]) + (f"_plus{len(signaling_param)-3}" if len(signaling_param) > 3 else "")
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_heatmap_{pattern}_{pathway_str}.pdf"
                                fig_heatmap.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="シグナル貢献度ヒートマップをPDFでダウンロード",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_heatmap_{pattern}_{pathway_str}.pdf',
                                        mime='application/octet-stream'
                                    )
                            
                            except Exception as e:
                                st.error(f"ヒートマップの生成中にエラーが発生しました: {str(e)}")
                                st.error(traceback.format_exc())

                with tabs[7]:
                    st.markdown("#### Network circle plot")
                    
                    # シグナル選択（Aggregate と各パスウェイ）
                    if 'netP' in result and 'pathways' in result['netP']:
                        # 選択肢: "Aggregate" または 個別パスウェイの複数選択
                        option_type = st.radio(
                            "Pathway type:",
                            ["Aggregate", "Specific pathways"],
                            horizontal=True
                        )
                        
                        if option_type == "Aggregate":
                            selected_pathway = "Aggregate"
                            # Interaction measure を選択するラジオボタン
                            aggregate_measure = st.radio(
                                "Select measure:",
                                ["Interaction Weight", "Interaction Number"],
                                horizontal=True
                            )
                            measure_key = "weight" if aggregate_measure == "Interaction Weight" else "count"
                        else:
                            # 複数のパスウェイを選択可能にする
                            pathways = list(result['netP']['pathways'])
                            selected_pathways = st.multiselect(
                                "Select pathways (selected pathways are summed):",
                                pathways,
                                default=[pathways[0]] if pathways else []
                            )
                            selected_pathway = selected_pathways  # リストとして渡す
                            measure_key = "weight"  # 特定パスウェイでは常にweight
                    else:
                        st.warning("パスウェイ情報が利用できません")
                        selected_pathway = None
                        measure_key = "weight"

                    circle_type = st.radio(
                            "Circle plot type:",
                            ["Total network", "Individual cell-type"],
                            horizontal=True
                        )
                    
                    # 必要に応じて vertex_weight も取得（例：result['adata'].obs から各細胞の数）
                    vertex_weight = None
                    if hasattr(result['adata'].obs, 'value_counts'):
                        try:
                            cell_counts = result['adata'].obs[result['adata'].obs.columns[0]].value_counts()
                            vertex_weight = [cell_counts.get(ct, 1) for ct in result['net']['dimnames'][0]]
                        except Exception as e:
                            st.warning(f"vertex_weight の取得に失敗しました: {str(e)}")

                    if st.button("Generate circle plot"):
                        if option_type == "Specific pathways" and (not selected_pathways or len(selected_pathways) == 0):
                            st.warning("Please select at least one pathway.")
                        else:
                            try:
                                if option_type == "Aggregate":
                                    # 集計ネットワークの場合
                                    title = f"Cell-Cell {aggregate_measure}"
                                    with st.spinner("集計サークルプロットを生成中..."):
                                        if circle_type == "Total network":
                                            fig_circle = netVisual_circle(
                                                net=result,
                                                title_name=title,
                                                measure=measure_key,
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                arrow=True
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                measure=measure_key,
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                arrow=True,
                                                vertex_size_max=6
                                            )
                                else:
                                    # 個別シグナルまたは複数シグナルの場合
                                    if len(selected_pathways) == 1:
                                        title = f"Cell-Cell Interaction: {selected_pathways[0]}"
                                    else:
                                        if len(selected_pathways) <= 3:
                                            title = f"Combined Pathways: {', '.join(selected_pathways)}"
                                        else:
                                            title = f"Combined {len(selected_pathways)} Pathways"
                                    
                                    with st.spinner(f"サークルプロットを生成中... ({title})"):
                                        if circle_type == "Total network":
                                            fig_circle = netVisual_circle(
                                                net=result,
                                                title_name=title,
                                                pathway_name=selected_pathway,
                                                measure="weight",  # 個別の場合は通常 weight を使用
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                arrow=True
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                pathway_name=selected_pathway,
                                                measure="weight",  # 個別の場合は通常 weight を使用
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                arrow=True,
                                                vertex_size_max=6
                                            )
                                st.pyplot(fig_circle)
                                
                                # PDF 保存・ダウンロード
                                import io
                                pdf_buffer = io.BytesIO()
                                fig_circle.savefig(pdf_buffer, format="pdf", bbox_inches='tight')
                                pdf_data = pdf_buffer.getvalue()
                                
                                # ファイル名の設定
                                if option_type == "Aggregate":
                                    filename = 'cellchat_circle_aggregate.pdf'
                                else:
                                    # 複数パスウェイの場合はパスウェイ数を表示
                                    if len(selected_pathway) <= 2:
                                        pathway_str = '_'.join(selected_pathway)
                                    else:
                                        pathway_str = f'{len(selected_pathway)}_pathways'
                                    filename = f'cellchat_circle_{pathway_str}_{circle_type.replace(" ", "_")}.pdf'
                                
                                st.download_button(
                                    label="Download circle plot",
                                    data=pdf_data,
                                    file_name=filename,
                                    mime='application/octet-stream'
                                )
                            except Exception as e:
                                st.error(f"サークルプロット生成エラー: {str(e)}")
                                import traceback
                                st.error(traceback.format_exc())


                with tabs[8]:
                    st.markdown("#### シグナル役割分析")
                    
                    st.info("""
                    シグナル役割分析は、通信ネットワークにおける各細胞タイプの
                    送信者(Sender)、受信者(Receiver)、仲介者(Mediator)、影響者(Influencer)としての役割を示します。
                    ヒートマップの色は各細胞タイプの相対的な重要度を表しています（各行で最大値を1に正規化）。
                    """)
                    # 2カラムに分ける
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        if 'netP' in result and 'pathways' in result['netP']:
                            pathways = result['netP']['pathways']
                            selected_path = st.selectbox(
                                "Select pathways for heatmap:",
                                options=["Aggregate"] + sorted(list(pathways)),
                                index=0
                            )
                                                    # デバッグ用に詳細を出力
                            print("Pathways:", pathways)
                            print("Pathways type:", type(pathways))
                            print("Pathways length:", len(pathways))
                        show_value_role = st.checkbox("Show values in plot?", value=False)
                    with col2:
                        font_size = st.slider("Font size:", min_value=8, max_value=14, value=10)
                    
                        # ヒートマップのサイズ設定
                        width = st.slider("Fig width:", min_value=6, max_value=15, value=10)
                        height = st.slider("Fig height:", min_value=3, max_value=8, value=4)

                        pathways = result['netP']['pathways']
    

                    if st.button("Generate plot"):
                        with st.spinner("Calculating..."):
                            try:
                                from pages.cellchat_vis import netAnalysis_signalingRole_network
                                
                                # 役割分析の実行
                                fig_role = netAnalysis_signalingRole_network(
                                    result,  # オブジェクト全体を渡す
                                    signaling=selected_path,
                                    font_size=font_size,
                                    width=width,
                                    height=height,
                                    color_heatmap=color_heatmap,
                                    sorted_order=sorted_order,  # 細胞順序の指定
                                    show_value=show_value_role
                                )
                                
                                # 図を表示
                                st.pyplot(fig_role)
                                
                                # PDFとして保存
                                pdf_path = f"{cellchat_temp_dir}/signaling_role_{selected_path}.pdf"
                                fig_role.savefig(pdf_path, bbox_inches='tight')
                                
                                # PDFダウンロードボタンを提供
                                with open(pdf_path, "rb") as pdf_file:
                                    pdf_bytes = pdf_file.read()
                                    st.download_button(
                                        label="Download PDF",
                                        data=pdf_bytes,
                                        file_name=f'cellchat_signaling_role_{selected_path}.pdf',
                                        mime='application/octet-stream'
                                    )
                                    
                            except Exception as e:
                                st.error(f"シグナル役割分析エラー: {str(e)}")
                                st.error(traceback.format_exc())

                # Add the Gene Expression Analysis tab
                with tabs[9]:
                    st.markdown("#### 遺伝子発現分析")
                                    
                    st.info("""
                    特定のシグナリングパスウェイに関わる遺伝子の発現を細胞タイプ間で可視化します。
                    """)

                    if 'cellchatdb' not in locals() or cellchatdb is None:
                        gene_list = result['adata'].var_names.tolist()
                        species_index = check_species_index(gene_list)
                        species = 'human' if species_index == 1 else 'mouse'
                        try:
                            with st.spinner(f'{species} CellChatDBを取得中...'):
                                cellchatdb = get_cellchatdb_from_r(species=species)
                                st.success(f"{species} CellChatDBを正常に取得しました")
                        except Exception as e:
                            st.error(f"CellChatDBの取得中にエラーが発生しました: {str(e)}")
                            cellchatdb = None

                    # adataがない場合のアップロード処理を追加
#                    if 'adata' not in locals() or adata is None:
                        # AnnDataファイルのアップロード
                    #    uploaded_adata = st.file_uploader("AnnDataファイルをアップロード", type=['.h5ad'])
                        
                   #     if uploaded_adata is not None:
                   #         try:
                   #             # ファイルをAnnDataとして読み込み
                   #             adata = sc.read_h5ad(uploaded_adata)
                   #             st.write('adata')
                   #             st.write(adata)
                   #             
                                # check_species_index関数を使用して種を推定
                   #             gene_list = adata.var_names.tolist()
                   #             species_index = check_species_index(gene_list)
                   #             species = 'human' if species_index == 1 else 'mouse'
                                
                   #             st.write(f"推定された種: {species}")
                                
                                # CellChatDBを取得
                           #     try:
                           #         with st.spinner(f'{species} CellChatDBを取得中...'):
                           #             cellchatdb = get_cellchatdb_from_r(species=species)
                           #             st.success(f"{species} CellChatDBを正常に取得しました")
                           #     except Exception as e:
                           #         st.error(f"CellChatDBの取得中にエラーが発生しました: {str(e)}")
                           #         cellchatdb = None
                            
                     #       except Exception as e:
                     #           st.error(f"AnnDataの読み込みエラー: {str(e)}")
                     #           adata = None
                     #           cellchatdb = None


                    # パスウェイ選択（Aggregate または特定パスウェイ）
                    if ('netP' in result) and ('pathways' in result['netP']):
                        pathways = result['netP']['pathways']
                        selected_path_expr = st.selectbox(
                            "Pathways to show gene expression:",
                            options=sorted(list(pathways)),
                            index=0
                        )
                      
                        if st.button("遺伝子発現を分析"):
                            from pages.cellchat_vis import plotGeneExpression
                            with st.spinner("遺伝子発現を分析中..."):
                                # adataが存在しない場合は警告を表示して処理を中断
                                if 'adata' not in result or result['adata'] is None:
                                    st.warning("AnnDataが読み込まれていません。")
                                else:
                                    try:
                                        fig_expr = plotGeneExpression(
                                            result,
                                            signaling=selected_path_expr,
                                          #  adata=adata,
                                            cellchatdb=cellchatdb,  # CellChatDBを渡す
                                            group_by=groupby,
                                            cmap_name=cmap_name
                                        )
                                        
                                        st.pyplot(fig_expr)
                                        
                                        # Save to PDF
                                        pdf_path = f"{cellchat_temp_dir}/gene_expression_{selected_path_expr}.pdf"
                                        fig_expr.savefig(pdf_path, bbox_inches='tight')
                                        
                                        with open(pdf_path, "rb") as pdf_file:
                                            pdf_bytes = pdf_file.read()
                                            st.download_button(
                                                label="遺伝子発現分析PDFをダウンロード",
                                                data=pdf_bytes,
                                                file_name=f'cellchat_gene_expression_{selected_path_expr}.pdf',
                                                mime='application/octet-stream'
                                            )
                                    except Exception as e:
                                        st.error(f"遺伝子発現分析エラー: {str(e)}")
                    else:
                        st.warning("結果にパスウェイ情報がありません。")
                
                # 可視化セクションの最後、解析を再実行ボタンの直前に以下を追加:
                st.markdown("---")

                # 解析結果の保存ボタン
                if st.session_state.cellchat_res is not None:
                    if st.button('Prepare CellChat analysis result to download'):
                        try:
                            # ファイル名の作成
                            if uploaded_file:
                                base_filename = sanitize_filename(uploaded_file.name, 20)
                            else:
                                base_filename = "cellchat_result"
                            
                            # シグナルタイプ文字列の作成
                            signal_str = "_".join([s[:3] for s in selected_types])
                            
                            # 保存用データを準備
                            import pickle
                            output = pickle.dumps(st.session_state.cellchat_res)
                            
                            # ダウンロードボタンを表示
                            st.download_button(
                                label="Download the result file",
                                data=output,
                                file_name=f"cellchat_{base_filename}_{signal_str}.pkl",
                                mime="application/octet-stream"
                            )
                            
                            st.success("解析結果を保存しました。ダウンロードするにはボタンをクリックしてください。")
                        except Exception as e:
                            st.error(f"保存中にエラーが発生しました: {str(e)}")


                if st.button('Rerun analysis'):
                    st.session_state.cellchat_res = None
                    st.rerun()

    else:
        st.info("h5adファイルをアップロードして開始してください。")
    