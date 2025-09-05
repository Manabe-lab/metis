import streamlit as st

# Set page config - must be the first Streamlit command
st.set_page_config(page_title="CellChat Permutation Comparison", page_icon="💬", layout="wide")

import scanpy as sc
import numpy as np
import pandas as pd
import os
import pickle
import re
import time
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
import zipfile
import io
from joblib import Parallel, delayed
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm
import copy
import gc
from streamlit_sortables import sort_items
import logging
import traceback
import scipy
from pages.cellchat import (
    get_cellchatdb_from_r,  debug_cellchatdb, check_species_index, identify_overexpressed_genes, extractGene,
    preprocess_data,  computeExpr_LR,computeExpr_coreceptor, compute_hill_outer, computeExpr_agonist, 
    computeExpr_antagonist, precompute_gene_expressions, computeCommunProbPathway,
    aggregateCell_Cell_Communication
)

from pages.cellchat_vis import (
    netVisual_circle, netVisual_circle_individual, netVisual_chord_by_pathway,
    plotGeneExpression, plot_lr_contribution,
    netVisual_heatmap, create_chord_diagram_pycirclize, netAnalysis_signalingRole_scatter,
    netAnalysis_signalingChanges_scatter
)


logger = logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context")
logger.disabled = True
#import warnings
#warnings.filterwarnings("ignore", message=".*missing ScriptRunContext.*")

# GPU初期化用関数：gpu_precision が True なら float64、False なら float32 とする
def init_gpu(use_gpu, gpu_precision):
    has_gpu = False
    cp = None
    gpu_dtype = None
    if use_gpu:
        try:
            import cupy as cp
            has_gpu = True
            # gpu_precision が True なら float64、False なら float32
            gpu_dtype = cp.float64 if gpu_precision else cp.float32
            st.write("GPU available for permutation testing")
        except (ImportError, ModuleNotFoundError):
            st.write("GPU requested but cupy not available, falling back to CPU")
            use_gpu = False
            cp = None
    return use_gpu, has_gpu, cp, gpu_dtype


# メインのパーミュテーションテスト関数
def run_permutation_cellchat(adata1_filtered, adata2_filtered, groupby,
                        gene_use,
                        complex_input,
                        cofactor_input,
                        resource, 
                        cell_types, n_perm=100, n_jobs=4, 
                        use_gpu=False, gpu_precision=True, hybrid_mode=True, do_every_deg=True,
                        union_mode="individual",feature_list_g1=None, feature_list_g2=None, r_patcher=False, **kwargs):
    # GPU初期化：gpu_precision は bool 型
    use_gpu, has_gpu, cp, gpu_dtype = init_gpu(use_gpu, gpu_precision)


    # ２つのAnnDataを結合し、グループラベルを作成
    adata_combined = adata1_filtered.concatenate(adata2)
    group_labels = np.zeros(adata_combined.n_obs)
    group_labels[:adata1_filtered.n_obs] = 0  # Group 1
    group_labels[adata1_filtered.n_obs:] = 1  # Group 2



    # 実際の解析（実データ）の実行
    st.write("Calculating actual difference between groups...")
    with st.spinner("Running CellChat analysis for Group 1"):
        if feature_list_g1 is not None: #もしdo_every_degの場合
            kwargs["features"] = feature_list_g1
        result1 = cellchat_minimal_analysis_optimized(
            adata1_filtered, groupby=groupby,
            gene_use=gene_use,
            complex_input=complex_input,
            cofactor_input=cofactor_input,
            resource=resource, 
            use_gpu=use_gpu, gpu_precision=gpu_precision, 
            validate_gpu_results=False, **kwargs
        )
    with st.spinner("Running CellChat analysis for Group 2"):
        if feature_list_g2 is not None:
            kwargs["features"] = feature_list_g2
        result2 = cellchat_minimal_analysis_optimized(
            adata2_filtered, groupby=groupby,
            gene_use=gene_use,
            complex_input=complex_input,
            cofactor_input=cofactor_input,
            resource=resource, 
            use_gpu=use_gpu, gpu_precision=gpu_precision, 
            validate_gpu_results=False, **kwargs
        )
    actual_diff = calculate_actual_diff(result1, result2, cell_types)

    # パーミュテーション生成：オリジナルの実装
    if use_gpu and has_gpu:
        try:
            permutation_gpu = cp.stack([cp.random.permutation(adata_combined.n_obs) for _ in range(n_perm)], axis=1)
            permutation = cp.asnumpy(permutation_gpu)
            cp.get_default_memory_pool().free_all_blocks()
        except Exception as e:
            st.warning(f"GPU permutation generation failed: {str(e)}, falling back to CPU")
            permutation = np.array([np.random.permutation(adata_combined.n_obs) for _ in range(n_perm)]).T
    else:
        permutation = np.array([np.random.permutation(adata_combined.n_obs) for _ in range(n_perm)]).T

    # バッチ処理の設定
    batch_size = determine_optimal_batch_size(n_perm, n_jobs)
    num_batches = (n_perm + batch_size - 1) // batch_size
    st.write(f"Running {n_perm} permutations in {num_batches} batches using {n_jobs} cores...")
    st.write(f"Batch size: {batch_size} permutations per batch")

    # 進捗表示用プレースホルダー
    batch_progress_placeholder = st.empty()
    overall_progress_bar = st.progress(0)
    batch_progress_bar = st.progress(0)

    permutation_results = []
    valid_results_count = 0
    start_time = time.time()

    # 各バッチごとに並列処理でパーミュテーション解析を実行
    for batch_idx in range(num_batches):
        start_idx = batch_idx * batch_size
        end_idx = min((batch_idx + 1) * batch_size, n_perm)
        batch_count = end_idx - start_idx

        batch_progress_placeholder.text(f"Running batch {batch_idx+1}/{num_batches}: permutations {start_idx} to {end_idx-1}")
        batch_start_time = time.time()

        batch_results = Parallel(n_jobs=n_jobs, verbose=10)(
            delayed(run_permutation_batch_gpu)(
                batch_idx, start_idx + i, min(start_idx + i + 1, end_idx),
                adata_combined, group_labels, groupby,
                gene_use,
                complex_input,
                cofactor_input,
                resource, 
                cell_types,
                adata1_filtered.n_obs, permutation,
                use_gpu=use_gpu and (i % 4 == 0 if hybrid_mode else True),
                gpu_precision=gpu_precision,
                do_every_deg=do_every_deg,
                union_mode=union_mode,
                r_patcher=r_patcher,
                kwargs=kwargs
            ) for i in range(batch_count)
        )
        # バッチ結果をフラット化
        flattened_results = [result for sublist in batch_results for result in sublist if result is not None]
        valid_results_count += len(flattened_results)
        permutation_results.extend(flattened_results)

        overall_progress_bar.progress((batch_idx + 1) / num_batches)
        batch_progress_bar.progress(1.0)
        batch_elapsed_time = time.time() - batch_start_time
        batch_progress_placeholder.text(
            f"Completed batch {batch_idx+1}/{num_batches}: {len(flattened_results)}/{batch_count} valid permutations in {batch_elapsed_time:.1f} sec."
        )
        if has_gpu and use_gpu:
            # Clear memory pool more frequently
            cp.get_default_memory_pool().free_all_blocks()
            
        time.sleep(0.1)
        batch_progress_bar.progress(0.0)

    elapsed_time = time.time() - start_time
    overall_progress_bar.progress(1.0)
    st.success(f"Completed {valid_results_count} valid permutations out of {n_perm} total in {elapsed_time:.1f} seconds.")

    
    p_values = calculate_p_values_vectorized(actual_diff, permutation_results)
    adjusted_p_values = adjust_p_values(p_values)

    return {
        "actual_diff": actual_diff,
        "permutation_results": permutation_results,
        "p_values": p_values,
        "adjusted_p_values": adjusted_p_values,
        "result1": result1,
        "result2": result2
    }

# バッチごとにパーミュテーション解析を実行する関数
def run_permutation_batch_gpu(
    batch_idx, start_idx, end_idx,
    adata_combined,
    group_labels,
    groupby,
    gene_use,
    complex_input,
    cofactor_input,
    resource,  
    cell_types, n1, permutation,
    use_gpu=False,
    gpu_precision=True,
    do_every_deg=True,
    union_mode="individual",
    deg_params=None,
    r_patcher=False,
    kwargs=None,

):
    print(f"Starting batch {batch_idx}: permutations {start_idx} to {end_idx-1}")
    batch_results = []

    use_gpu, has_gpu, cp, gpu_dtype = init_gpu(use_gpu, gpu_precision)

    if kwargs is None:
        kwargs = {}
    if deg_params is None:
        deg_params = {}

    for i in range(start_idx, end_idx):
        try:
            pid = os.getpid()
            print(f"[Batch {batch_idx}, Permutation {i}] Starting in process {pid}")

            # (1) パーミュテーションで分割
            permuted_indices = permutation[:, i]
            permuted_labels = group_labels[permuted_indices]
            mask1 = (permuted_labels == 0)
            mask2 = (permuted_labels == 1)
            adata_perm1 = adata_combined[mask1].copy()
            adata_perm2 = adata_combined[mask2].copy()

            # (2) do_every_deg=Trueなら、DEG計算 → 遺伝子セット決定
            #     その後 cellchat_minimal_analysis_optimized に features=... を指定して呼ぶ
            #     do_every_deg=False なら、外部kwargsの features をそのままor None として解析
            if do_every_deg: #多分、不要　繰り返し
                # もし外部で 'features' が渡されていても無視する
                if "features" in kwargs:
                    del kwargs["features"]

                # DEG 計算
                deg_res1 = identify_overexpressed_genes(adata_perm1, group_by=groupby, **deg_params)
                deg_res2 = identify_overexpressed_genes(adata_perm2, group_by=groupby, **deg_params)

                feats1 = set(deg_res1["features"])
                feats2 = set(deg_res2["features"])

                if union_mode == "union":
                    final_feats_1 = final_feats_2 = feats1.union(feats2)
                elif union_mode == "intersection":
                    final_feats_1 = final_feats_2 = feats1.intersection(feats2)
                elif union_mode == "individual":
                    # それぞれ別々に解析する
                    final_feats_1 = feats1
                    final_feats_2 = feats2
                else:
                    raise ValueError(f"Invalid union_mode: {union_mode}")
                
                # ここで GPU ロジックを判定しつつ、1回だけ呼ぶ
                use_gpu_now = (use_gpu and has_gpu)
                cellchat_res1 = cellchat_minimal_analysis_optimized(
                    adata_perm1,
                    groupby=groupby,
                    gene_use=gene_use,
                    complex_input=complex_input,
                    cofactor_input=cofactor_input,
                    resource=resource,
                    use_gpu=use_gpu_now,
                    gpu_precision=gpu_precision,
                    features=list(final_feats_1),
                    r_patcher=r_patcher,
                    **kwargs
                )
                cellchat_res2 = cellchat_minimal_analysis_optimized(
                    adata_perm2,
                    groupby=groupby,
                    gene_use=gene_use,
                    complex_input=complex_input,
                    cofactor_input=cofactor_input,
                    resource=resource,
                    use_gpu=use_gpu_now,
                    gpu_precision=gpu_precision,
                    features=list(final_feats_2),
                    r_patcher=r_patcher,
                    **kwargs
                )
            else:
                # do_every_deg=False の場合、features は外部kwargs 任せ（None なら全遺伝子）
                use_gpu_now = (use_gpu and has_gpu)
                cellchat_res1 = cellchat_minimal_analysis_optimized(
                    adata_perm1,
                    groupby=groupby,
                    gene_use=gene_use,
                    complex_input=complex_input,
                    cofactor_input=cofactor_input,
                    resource=resource,
                  #  db=db,
                    use_gpu=use_gpu_now,
                    gpu_precision=gpu_precision,
                    r_patcher=r_patcher,
                    **kwargs
                )
                cellchat_res2 = cellchat_minimal_analysis_optimized(
                    adata_perm2,
                    groupby=groupby,
                    gene_use=gene_use,
                    complex_input=complex_input,
                    cofactor_input=cofactor_input,
                    resource=resource,
                    use_gpu=use_gpu_now,
                    gpu_precision=gpu_precision,
                    r_patcher=r_patcher,
                    **kwargs
                )

            # (3) 差分計算
            print(f"[Batch {batch_idx}, Permutation {i}] Calculating differences in process {pid}")
            perm_diff = calculate_actual_diff(cellchat_res1, cellchat_res2, cell_types)
            print(f"[Batch {batch_idx}, Permutation {i}] Completed in process {pid}")

            # (4) 結果を追加
            batch_results.append(perm_diff)

            # GPUメモリ解放
            if use_gpu and has_gpu:
                cp.get_default_memory_pool().free_all_blocks()

        except Exception as e:
            print(f"[Batch {batch_idx}, Permutation {i}] Error: {str(e)}")
            import traceback
            print(traceback.format_exc())

    if use_gpu and has_gpu:
        cp.get_default_memory_pool().free_all_blocks()

    return batch_results


def cellchat_minimal_analysis_optimized(
    adata,
    groupby,
    gene_use, #これはLR genesの情報
    complex_input,
    cofactor_input,
    resource,
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
    seed=1,
    n_jobs=8,
    key_added="cellchat_res",
    trim=0.1,
    apply_pval_filter=True,
    features=None,
    use_gpu=False,
    gpu_precision='float64',
    validate_gpu_results=False,
    gpu_memory_limit=0.8,
    optimize_memory=True,
    r_patcher=False
):
    """
    GPU-optimized minimal version of CellChat algorithm for permutation testing.
    
    Parameters
    ----------
    adata : AnnData
        AnnData object
    groupby : str
        Column name in .obs containing cell type/cluster information
    db : dict
        CellChatDB dictionary
    use_layer : str, optional
        Data layer to use
    [他の既存パラメータ]
    use_gpu : bool, optional
        Whether to use GPU acceleration. Default is True
    gpu_precision : str, optional
        Precision for GPU calculations ('float64' or 'float32'). Default is 'float64'
    validate_gpu_results : bool, optional
        Whether to validate GPU results against CPU. Default is False
    gpu_memory_limit : float, optional
        GPU memory usage limit (0-1). Default is 0.8
    optimize_memory : bool, optional
        Whether to optimize memory usage. Default is True
        
    Returns
    -------
    dict
        Minimal cell-cell communication analysis results
    """
    
    # Initialize logger
    logger = logging.getLogger("CellChat-Minimal-Optimized")
    logger.setLevel(logging.INFO)
    
    # GPU初期化
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
            gpu_dtype = cp.float64 if gpu_precision else cp.float32
            logger.info(f"GPU acceleration enabled using {gpu_precision} precision")
        except ImportError:
            st.warning("GPU acceleration requested but cupy not available, falling back to CPU")
            use_gpu = False
    
    try:
        # 乱数のシードを設定
        np.random.seed(seed)
        if has_gpu and use_gpu:
            cp.random.seed(seed)

  #      gene_use, resource, complex_input, cofactor_input, gene_info = extractGene(db)
        
        # GPU用ヘルパー関数のインポート（利用可能な場合）
        if has_gpu and use_gpu:
            try:
                from pages.gpu_cellchat_helpers import (
                    compute_hill_outer_gpu,
                    precompute_gene_expressions_gpu,
                    geometricMean_gpu
                )
            except ImportError:
                st.warning("GPU helper functions not available, falling back to CPU implementations")
        
        # データの前処理
        if features is not None:
            # 提供された遺伝子リストを使用
            logger.info(f"Using provided gene list: {len(features)} genes")
            # 存在する遺伝子のみを保持
            
            adata_filtered, resource_filtered = preprocess_data(adata, groupby, complex_input, gene_use=gene_use, min_cells=min_cells,
                thresh_pct=expr_prop, resource=resource, features =features)
            
        else:
            adata_filtered = preprocess_data(adata, groupby, gene_use=gene_use, min_cells=min_cells,
                                                    tthresh_pct=expr_prop, resource=resource, features =None)
        
        
        # 発現データの取得
        if use_layer is not None and use_layer in adata_filtered.layers:
            logger.info(f"Using layer '{use_layer}'")
            X = adata_filtered.layers[use_layer]
        else:
            logger.info("Using default X matrix")
            X = adata_filtered.X
        
        # スパース行列を密行列に変換（効率的なアプローチとGPUオプション）
        if scipy.sparse.issparse(X):
            data_size = X.shape[0] * X.shape[1]
            
            # GPUを使用する場合、大きいスパース行列は直接GPUで変換を試みる
            if has_gpu and use_gpu and data_size < cp.cuda.Device().mem_info[0] * 0.3:
                try:
                    import cupyx.scipy.sparse as cusparse
                    logger.info("Converting sparse matrix on GPU")
                    X_gpu = cusparse.csr_matrix(X)
                    X = cp.asnumpy(X_gpu.toarray())
                    del X_gpu
                    cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    logger.warning(f"GPU sparse conversion failed: {str(e)}, falling back to CPU")
                    # CPUでの効率的な変換
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
                # CPUでの効率的な変換
                if optimize_memory and X.shape[0] > 5000:
                    chunk_size = 5000
                    result = []
                    for i in range(0, X.shape[0], chunk_size):
                        end = min(i + chunk_size, X.shape[0])
                        result.append(X[i:end].toarray())
                    X = np.vstack(result)
                else:
                    X = X.toarray()
        
        # 細胞タイプラベルの取得
        cell_labels = adata_filtered.obs[groupby].copy()
        cell_types = np.array(sorted(cell_labels.unique()))
        
        logger.info(f"Number of cell types: {len(cell_types)}")
        if len(cell_types) < 2:
            raise ValueError(f"Only {len(cell_types)} cell type found. At least 2 required.")
        
        # CPUとGPUの両方でデータ正規化
      #  data_use_cpu = X / np.max(X).astype(np.float64)
        data_use_cpu = X.astype(np.float64)
        
        if has_gpu and use_gpu:
            # GPU版
            X_gpu = cp.array(X, dtype=gpu_dtype)
            max_val_gpu = cp.max(X_gpu).astype(gpu_dtype)
            data_use_gpu = X_gpu / max_val_gpu
            
            # 検証
            if validate_gpu_results:
                # GPU結果とCPU結果の比較
                data_use_gpu_cpu = cp.asnumpy(data_use_gpu)
                max_diff = np.max(np.abs(data_use_cpu - data_use_gpu_cpu))
                if max_diff > 1e-5:
                    logger.warning(f"GPU and CPU normalization results differ: {max_diff}. Using CPU.")
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
        
        # メモリ最適化のため、元データを解放
        if optimize_memory:
            del X
            gc.collect()
            if has_gpu and use_gpu:
                if 'X_gpu' in locals():
                    del X_gpu
                cp.get_default_memory_pool().free_all_blocks()
        
        # 平均関数の定義
        if type_mean == "triMean":
            # GPU対応のtrimeean
            def FunMean_gpu(x):
                if has_gpu and use_gpu and isinstance(x, cp.ndarray):
                    x_no_nan = x[~cp.isnan(x)]
                    if len(x_no_nan) == 0:
                        return cp.nan
                    q = cp.quantile(x_no_nan, cp.array([0.25, 0.5, 0.5, 0.75]))
                    return cp.mean(q)
                else:
                    x_no_nan = x[~np.isnan(x)]
                    if len(x_no_nan) == 0:
                        return np.nan
                    q = np.quantile(x_no_nan, [0.25, 0.5, 0.5, 0.75], method='linear')
                    return np.mean(q)
            
            # オリジナルのCPU版
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
            # GPU版も同様に定義可能
        elif type_mean == "median":
            def FunMean(x):
                return np.median(x)
        else:
            def FunMean(x):
                return np.mean(x)
        
        # 細胞タイプごとの平均発現量の計算（最適化）
        # 事前に細胞インデックスを計算して繰り返し検索を避ける
        cell_counts = {}
        cell_type_indices = {}
        for cell_type in cell_types:
            indices = np.where(cell_labels == cell_type)[0]
            cell_counts[cell_type] = len(indices)
            cell_type_indices[cell_type] = indices
        
        # 有効な細胞タイプを特定
    #    cell_types = [ct for ct in cell_types if cell_counts.get(ct, 0) >= min_cells]
    #    if len(cell_types) < 2:
    #        raise ValueError(f"Less than 2 cell types have at least {min_cells} cells.")
        
        # GPU対応の平均発現計算
        data_use_avg_dict = {}
        
        # GPUを使用する場合の効率的な平均発現計算
        if has_gpu and use_gpu and cp is not None:
            for cell_type in cell_types:
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
                        temp_result = cp.apply_along_axis(FunMean_gpu, 0, data_subset)
                        avg_expr[i:end] = temp_result
                else:
                    data_subset = data_use[indices_gpu]
                    avg_expr = cp.apply_along_axis(FunMean_gpu, 0, data_subset)
                
                # GPU結果の検証
                if validate_gpu_results:
                    cpu_indices = np.array(indices)
                    cpu_data_subset = data_use_cpu[cpu_indices]
                    cpu_avg_expr = np.apply_along_axis(FunMean, 0, cpu_data_subset)
                    
                    # CPU結果とGPU結果の比較
                    avg_expr_cpu = cp.asnumpy(avg_expr)
                    max_diff = np.max(np.abs(cpu_avg_expr - avg_expr_cpu))
                    if max_diff > 1e-5:
                        logger.warning(f"GPU and CPU mean calculation differ: {max_diff}. Using CPU.")
                        data_use_avg_dict[cell_type] = cpu_avg_expr
                    else:
                        data_use_avg_dict[cell_type] = cp.asnumpy(avg_expr)
                else:
                    data_use_avg_dict[cell_type] = cp.asnumpy(avg_expr)
        else:
            # CPU版の計算（最適化）
            for cell_type in cell_types:
                indices = cell_type_indices[cell_type]
                
                # 大きなデータセットは分割処理
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
        
        # DataFrameに変換
        data_use_avg_df = pd.DataFrame(data_use_avg_dict, index=adata_filtered.var_names)
        
        # 遺伝子名からインデックスへのマッピング
        gene_to_index = {gene: i for i, gene in enumerate(adata_filtered.var_names)}
        
        # リガンド・レセプター発現レベルの計算 - オリジナル関数を使用
        dataLavg = computeExpr_LR(resource['ligand'].values, data_use_avg_df, complex_input)
        dataRavg = computeExpr_LR(resource['receptor'].values, data_use_avg_df, complex_input)
        
        # 共受容体の計算 - オリジナル関数を使用
        dataRavg_co_A_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "A")
        dataRavg_co_I_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "I")
        dataRavg = dataRavg * dataRavg_co_A_receptor / dataRavg_co_I_receptor
        
        # 集団サイズ効果が有効な場合 - オリジナルと同じ計算
        if population_size:
            cell_proportions = np.array([np.sum(cell_labels == ct) for ct in cell_types]) / nC
            dataLavg2 = np.tile(cell_proportions, (len(resource), 1))
            dataRavg2 = dataLavg2
        else:
            dataLavg2 = np.ones((len(resource), len(cell_types)))
            dataRavg2 = np.ones((len(resource), len(cell_types)))
        
        # アゴニスト/アンタゴニストのインデックス特定 - オリジナルと同じ
        index_agonist = np.where(resource['agonist'].notna() & (resource['agonist'] != ""))[0] if 'agonist' in resource.columns else []
        index_antagonist = np.where(resource['antagonist'].notna() & (resource['antagonist'] != ""))[0] if 'antagonist' in resource.columns else []
        
        # 置換データの準備 - GPU対応
        permutation = np.zeros((nC, nboot), dtype=np.int32)
        
        # オリジナルの実装: シンプルなNumPy permutation
        if has_gpu and use_gpu:
            # GPUでの生成
            for i in range(nboot):
                perm_gpu = cp.random.permutation(nC)
                permutation[:, i] = cp.asnumpy(perm_gpu)
        else:
            # CPU版
            for i in range(nboot):
                permutation[:, i] = np.random.permutation(nC)
        
        # 遺伝子発現の事前計算 - GPU対応
        if nboot > 0:
            # GPUの利用可能性に応じて適切な関数を選択
            if has_gpu and use_gpu:
                try:
                    all_gene_expr = precompute_gene_expressions_gpu(
                        data_use_cpu,  # 常にCPUデータを渡す
                        cell_labels, 
                        permutation, 
                        cell_types, 
                        FunMean_gpu if has_gpu and use_gpu else FunMean, 
                        nboot, 
                        n_jobs,
                        gpu_dtype=gpu_dtype
                    )
                    
                    if validate_gpu_results:
                        # サンプルの一部でCPU計算と比較検証
                        sample_size = min(5, nboot)
                        cpu_result = precompute_gene_expressions(
                            data_use_cpu,
                            cell_labels,
                            permutation[:, :sample_size],
                            cell_types,
                            FunMean,
                            sample_size,
                            n_jobs
                        )
                        
                        # 最初のいくつかの結果を比較
                        max_diff = np.max(np.abs(cpu_result[:, :, :3] - all_gene_expr[:, :, :3]))
                        if max_diff > 1e-5:
                            logger.warning(f"GPU and CPU precomputation differ: {max_diff}. Using CPU.")
                            all_gene_expr = precompute_gene_expressions(
                                data_use_cpu,
                                cell_labels,
                                permutation,
                                cell_types,
                                FunMean,
                                nboot,
                                n_jobs
                            )
                except Exception as e:
                    logger.warning(f"GPU gene expression precomputation failed: {str(e)}. Using CPU.")
                    all_gene_expr = precompute_gene_expressions(
                        data_use_cpu,
                        cell_labels,
                        permutation,
                        cell_types,
                        FunMean,
                        nboot,
                        n_jobs
                    )
            else:
                # CPU版
                all_gene_expr = precompute_gene_expressions(
                    data_use_cpu,
                    cell_labels,
                    permutation,
                    cell_types,
                    FunMean,
                    nboot,
                    n_jobs
                )
        
            # メモリ最適化のための解放
            if optimize_memory:
                # エラー箇所の修正: cp が None でないことを確認
                if has_gpu and use_gpu and cp is not None and isinstance(data_use, cp.ndarray):
                    del data_use
                    cp.get_default_memory_pool().free_all_blocks()
                else:
                    # CP が None または data_use が cp.ndarray でない場合
                    if 'data_use' in locals():
                        del data_use
                # data_use_cpuは他の計算でまだ必要な場合は保持
                gc.collect()
        
        # 複合体マッピングの事前計算 - 結果に影響しない最適化
        complex_mapping = {}
        if not complex_input.empty:
            for complex_name in complex_input.index:
                subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[complex_name, subunits_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in gene_to_index]
                
                if subunits:
                    complex_mapping[complex_name] = [gene_to_index[s] for s in subunits]
        
        # 確率と有意性の行列を初期化
        numCluster = len(cell_types)
        nLR = len(resource)
        Prob = np.zeros((numCluster, numCluster, nLR))
        Pval = np.zeros((numCluster, numCluster, nLR))
        
        # リガンド・レセプター遺伝子のインデックス事前計算 - 結果に影響しない最適化
        ligand_indices = []
        receptor_indices = []
        
        for i in range(nLR):
            ligand = resource['ligand'].iloc[i]
            receptor = resource['receptor'].iloc[i]
            
            # 単一遺伝子か複合体かを判定
            if isinstance(ligand, str) and ligand in gene_to_index:
                ligand_indices.append((i, [gene_to_index[ligand]], False))
            elif isinstance(ligand, str) and ligand in complex_mapping:
                ligand_indices.append((i, complex_mapping[ligand], True))
            else:
                ligand_indices.append((i, [], None))
            
            if isinstance(receptor, str) and receptor in gene_to_index:
                receptor_indices.append((i, [gene_to_index[receptor]], False))
            elif isinstance(receptor, str) and receptor in complex_mapping:
                receptor_indices.append((i, complex_mapping[receptor], True))
            else:
                receptor_indices.append((i, [], None))
        
        # 各LRペアを処理 - GPU最適化
        for i in range(nLR):
            # Hill関数による相互作用確率の計算 - GPU or CPU
            if has_gpu and use_gpu:
                try:
                    # GPU実装
                    dataLavg_gpu = cp.array(dataLavg[i, :], dtype=gpu_dtype)
                    dataRavg_gpu = cp.array(dataRavg[i, :], dtype=gpu_dtype)
                    
                    # GPU用Hill関数の呼び出し
                    P1_gpu = compute_hill_outer_gpu(dataLavg_gpu, dataRavg_gpu, k, n, dtype=gpu_dtype)
                    
                    # CPU結果と比較検証（オプション）
                    if validate_gpu_results:
                        P1_cpu = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
                        P1_gpu_cpu = cp.asnumpy(P1_gpu)
                        max_diff = np.max(np.abs(P1_cpu - P1_gpu_cpu))
                        
                        if max_diff > 1e-5:
                            logger.warning(f"GPU and CPU Hill function differ: {max_diff}. Using CPU.")
                            P1 = P1_cpu
                        else:
                            P1 = cp.asnumpy(P1_gpu)
                    else:
                        P1 = cp.asnumpy(P1_gpu)
                    
                    # メモリ解放
                    del dataLavg_gpu, dataRavg_gpu, P1_gpu
                    if i % 50 == 0:  # 定期的にメモリ解放
                        cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    logger.warning(f"GPU Hill function failed: {str(e)}. Using CPU.")
                    P1 = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
            else:
                # CPU実装
                P1 = compute_hill_outer(dataLavg[i, :], dataRavg[i, :], k, n)
            
            # アゴニスト効果 - オリジナル関数を使用
            P2 = np.ones((numCluster, numCluster))
            if i in index_agonist:
                data_agonist = computeExpr_agonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P2 = np.outer(data_agonist, data_agonist)
            
            # アンタゴニスト効果 - オリジナル関数を使用
            P3 = np.ones((numCluster, numCluster))
            if i in index_antagonist:
                data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P3 = np.outer(data_antagonist, data_antagonist)
            
            # 集団サイズ効果 - オリジナルと同じ計算
            P4 = np.ones((numCluster, numCluster))
            if population_size:
                P4 = np.outer(dataLavg2[i, :], dataRavg2[i, :])
            
            # 最終確率 - オリジナルと同じ計算
            Pnull = P1 * P2 * P3 * P4
            Prob[:, :, i] = Pnull
            
            # 相互作用がない場合はp値計算をスキップ
            if np.sum(Pnull) == 0 or nboot == 0:
                Pval[:, :, i] = 1
                continue
            
            Pnull_vec = Pnull.flatten()
            
            # リガンド・レセプター情報の取得
            ligand_info = ligand_indices[i]
            receptor_info = receptor_indices[i]
            
            # 発現取得できない場合はスキップ
            if ligand_info[2] is None or receptor_info[2] is None:
                Pval[:, :, i] = 1
                continue
            
            # ブートストラップ確率の初期化
            Pboot = np.zeros((numCluster * numCluster, nboot))
            
            # メモリ効率化のためバッチサイズを調整
            batch_size = min(20, nboot)  
            
            # 並列処理で効率的に計算
            n_jobs_effective = min(n_jobs, os.cpu_count() or 1, nboot)
            
            # GPU使用を制御
            use_gpu_for_batches = has_gpu and use_gpu
            
            # バッチごとの処理
            from joblib import Parallel, delayed
            
            # 置換計算バッチ関数の定義
            def compute_permutation_batch(batch_indices, ligand_info, receptor_info,
                                         numCluster, n, k, population_size, cell_labels,
                                         permutation, cell_types, nC, use_gpu_compute=False):
                batch_results = np.zeros((numCluster * numCluster, len(batch_indices)))
                
                for idx, j in enumerate(batch_indices):
                    # リガンド発現の取得
                    lr_i, ligand_gene_indices, is_ligand_complex = ligand_info
                    
                    if not is_ligand_complex:
                        # 単一遺伝子
                        if ligand_gene_indices:
                            ligand_idx = ligand_gene_indices[0]
                            dataLavgB = all_gene_expr[ligand_idx, :, j].reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    else:
                        # 複合体 - 幾何平均計算
                        expr_values = np.array([all_gene_expr[idx, :, j] for idx in ligand_gene_indices])
                        if len(expr_values) > 0:
                            # 対数変換、平均算出、逆変換
                            log_values = np.log(expr_values + 1e-10)
                            dataLavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    
                    # レセプター発現の取得
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
                    
                    # 相互作用確率の計算 - GPU使用オプション
                    if use_gpu_compute and has_gpu:
                        try:
                            dataLavgB_gpu = cp.array(dataLavgB[0, :], dtype=gpu_dtype)
                            dataRavgB_gpu = cp.array(dataRavgB[0, :], dtype=gpu_dtype)
                            
                            # GPU計算
                            dataLRB_gpu = cp.outer(dataLavgB_gpu, dataRavgB_gpu)
                            P1_boot_gpu = dataLRB_gpu**n / (k**n + dataLRB_gpu**n)
                            
                            # CPUに戻す
                            P1_boot = cp.asnumpy(P1_boot_gpu)
                            
                            # クリーンアップ
                            del dataLavgB_gpu, dataRavgB_gpu, dataLRB_gpu, P1_boot_gpu
                            # メモリ解放
                            if idx % 10 == 0:
                                cp.get_default_memory_pool().free_all_blocks()
                        except Exception as e:
                            # エラー時はCPUにフォールバック
                            dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                            P1_boot = dataLRB**n / (k**n + dataLRB**n)
                    else:
                        # CPU計算
                        dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                        P1_boot = dataLRB**n / (k**n + dataLRB**n)
                    
                    # 置換用の簡略化アプローチ（オリジナルと同様）
                    P2_boot = np.ones((numCluster, numCluster))
                    P3_boot = np.ones((numCluster, numCluster))
                    
                    # 集団サイズ効果
                    P4_boot = np.ones((numCluster, numCluster))
                    if population_size:
                        group_boot = cell_labels.values[permutation[:, j]]
                        cell_proportions_boot = np.array([np.sum(group_boot == ct) for ct in cell_types]) / nC
                        P4_boot = np.outer(cell_proportions_boot, cell_proportions_boot)
                    
                    # 最終確率
                    Pboot_result = P1_boot * P2_boot * P3_boot * P4_boot
                    batch_results[:, idx] = Pboot_result.flatten()
                
                return batch_results
            
            # バッチ処理実行
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
                            cell_types,
                            nC,
                            use_gpu_compute=use_gpu_for_batches and (j % 5 == 0)  # 5バッチごとにGPU使用
                        ) for j in batch_indices
                    )
                    # 結果の結合
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results_list[j_idx][:, 0]
                else:
                    # 単一スレッド処理
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
                        cell_types,
                        nC,
                        use_gpu_compute=use_gpu_for_batches
                    )
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results[:, j_idx]
                
                # GPUメモリ解放
                if has_gpu and use_gpu and b_start % 100 == 0:
                    cp.get_default_memory_pool().free_all_blocks()
            
            # p値の計算 - GPU対応
            if has_gpu and use_gpu:
                try:
                    # GPU計算
                    Pnull_vec_gpu = cp.array(Pnull_vec, dtype=gpu_dtype)
                    Pboot_gpu = cp.array(Pboot, dtype=gpu_dtype)
                    
                    # 閾値を超える値のカウント
                    nReject_gpu = cp.sum(Pboot_gpu > cp.expand_dims(Pnull_vec_gpu, 1), axis=1)
                    p_gpu = nReject_gpu / nboot
                    
                    # CPU転送
                    p = cp.asnumpy(p_gpu)
                    
                    # 検証
                    if validate_gpu_results:
                        nReject_cpu = np.sum(Pboot > (np.expand_dims(Pnull_vec, 1) + 1.49e-8), axis=1)
                        p_cpu = nReject_cpu / nboot
                        
                        max_diff = np.max(np.abs(p_cpu - p))
                        if max_diff > 1e-5:
                            logger.warning(f"GPU and CPU p-values differ: {max_diff}. Using CPU.")
                            p = p_cpu
                    
                    # メモリ解放
                    del Pnull_vec_gpu, Pboot_gpu, nReject_gpu, p_gpu
                    if i % 50 == 0:
                        cp.get_default_memory_pool().free_all_blocks()
                except Exception as e:
                    logger.warning(f"GPU p-value calculation failed: {str(e)}. Using CPU.")
                    nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
                    p = nReject / nboot
            else:
                # CPU計算
                nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
                p = nReject / nboot
            
            Pval[:, :, i] = p.reshape(numCluster, numCluster)
        
        # 確率が0の場所はp値を1に設定
        Pval[Prob == 0] = 1
        
        # p値フィルタリングが有効なら適用

        if apply_pval_filter:
            if r_patcher:
                Prob[Pval >= (trim_threshold - 1.49e-8)] = 0
            else:
                Prob[Pval >= trim_threshold] = 0


        # Filter out communication for cell types with few cells
        cell_counts = {ct: np.sum(cell_labels == ct) for ct in cell_types}
        cell_excludes = [ct for ct, count in cell_counts.items() if count <= min_cells]

        if cell_excludes:
            print(f"The cell-cell communication related with the following cell groups are excluded due to the few number of cells: {cell_excludes}")
            
            # Get the indices of excluded cell types
            exclude_indices = [i for i, ct in enumerate(cell_types) if ct in cell_excludes]
            
            # Set entire rows and columns to zero for excluded cell types
            for idx in exclude_indices:
                Prob[idx, :, :] = 0  # Sender is excluded
                Prob[:, idx, :] = 0  # Receiver is excluded
            
            # Update p-values
            Pval[Prob == 0] = 1
        
        # 次元名の設定
        dimnames = [list(cell_types), list(cell_types), list(resource.index)]
        
        # パスウェイレベルでの通信確率を計算
        netP = computeCommunProbPathway({"prob": Prob, "pval": Pval}, resource, 
                                       thresh=trim_threshold, apply_pval_filter=apply_pval_filter)
        
        # 集計ネットワークの計算
        net_summary = aggregateCell_Cell_Communication({"prob": Prob, "pval": Pval}, 
                                                     cell_types, pval_threshold=0.05, 
                                                     apply_pval_filter=apply_pval_filter)
        
        logger.info("CellChat minimal analysis completed")
        
        # 結果データフレームの準備
        results_data = {
            'source': [],
            'target': [],
            'interaction_name': [],
            'ligand': [],
            'receptor': [],
            'prob': [],
            'pval': []
        }
        
        for i in range(len(cell_types)):
            for j in range(len(cell_types)):
                for k in range(nLR):
                    if Prob[i, j, k] > 0:
                        results_data['source'].append(cell_types[i])
                        results_data['target'].append(cell_types[j])
                        results_data['interaction_name'].append(resource.index[k])
                        results_data['ligand'].append(resource['ligand'].iloc[k])
                        results_data['receptor'].append(resource['receptor'].iloc[k])
                        results_data['prob'].append(Prob[i, j, k])
                        results_data['pval'].append(Pval[i, j, k])
        
        results_df = pd.DataFrame(results_data)
        
        # 最終的なGPUメモリ解放
        if has_gpu and use_gpu:
            cp.get_default_memory_pool().free_all_blocks()
        
        # 可視化に必要な最小構造を返す
        return {
            'adata': adata_filtered,
            'results': results_df,
            'net': {
                "prob": Prob,
                "pval": Pval,
                "dimnames": dimnames,
                "centr": None  # 中心性は計算しない（可視化には不要）
            },
            'netP': netP,
            'network': net_summary,
            'groupby': groupby 
        }
    
    except Exception as e:
        print(f"Error in minimal CellChat analysis: {str(e)}")
        print("Full traceback:")
        print(traceback.format_exc())
        logger.error(f"Error in minimal CellChat analysis: {str(e)}")
        logger.error(traceback.format_exc())
        
        # エラー時もGPUメモリ解放
        if has_gpu and use_gpu:
            cp.get_default_memory_pool().free_all_blocks()
            
        return {'error': str(e), 'traceback': traceback.format_exc()}



# Create temp directory for saving files
if "cellchat_temp_dir" not in st.session_state:
    st.session_state.cellchat_temp_dir = True
    cellchat_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    os.mkdir(cellchat_temp_dir)
    st.session_state.cellchat_temp_dir = cellchat_temp_dir
else:
    cellchat_temp_dir = st.session_state.cellchat_temp_dir

# Function for splitting AnnData by a grouping variable
def split_anndata(adata, split_key, group1, group2):
    """Split AnnData object into two based on a grouping variable"""
    # Create masks for each group
    mask1 = adata.obs[split_key] == group1
    mask2 = adata.obs[split_key] == group2
    
    # Create subsets
    adata1 = adata[mask1].copy()
    adata2 = adata[mask2].copy()
    
    return adata1, adata2

# Function for finding the first matching index
def find_first_index_or_default(lst, elements, default=0):
    """Find the first matching element in a list or return default index"""
    for element in elements:
        try:
            return lst.index(element)
        except ValueError:
            continue
    return default

def determine_optimal_batch_size(n_perms, n_jobs):
    """
    Determine optimal batch size based on number of permutations and available jobs
    
    Parameters
    ----------
    n_perms : int
        Total number of permutations
    n_jobs : int
        Number of parallel jobs
    
    Returns
    -------
    int
        Optimal batch size
    """
    # Aim for at least 4 batches to show progress, but no more than 20
#    min_batch_count = 4
#    max_batch_count = 20
    
    # Aim for batch size to be a multiple of n_jobs for optimal parallelization
#    ideal_batch_size = max(1, n_jobs * 2)
    
    # Ensure we have at least min_batch_count batches
#    batch_size = min(ideal_batch_size, n_perms // min_batch_count)
    
    # Ensure we don't have more than max_batch_count batches
#    batch_size = max(batch_size, n_perms // max_batch_count)
    
    # Ensure batch size is at least 1
#    batch_size = max(1, batch_size)

    # まずは n_perms を 5〜10バッチに分ける程度を基本方針に
    # ここでは 8 で割ってバッチ数を出す
    batch_size = n_perms // 8
    
    # さらに n_jobs を参考に、極端に小さなバッチや大きなバッチを回避
    # 例: バッチサイズが n_jobs * 2 より小さい場合は n_jobs * 2 に合わせる
    # （ジョブを分配しやすいバランス）
    batch_size = max(batch_size, n_jobs * 2)
    
    return max(1, batch_size)    
 #   return batch_size


def calculate_actual_diff(result1, result2, cell_types):
    """
    Calculate actual differences between two CellChat results
    
    Parameters
    ----------
    result1 : dict
        First CellChat result
    result2 : dict
        Second CellChat result
    cell_types : list
        List of cell types
        
    Returns
    -------
    dict
        Dictionary containing probability and weight differences
    """
    # Prepare dicts for probability and weight differences
    prob_diff = []
    weight_diff = []
    
    # Extract values from results
    if ('netP' in result1 and 'pathways' in result1['netP'] and
        'netP' in result2 and 'pathways' in result2['netP']):
        
        # Get all pathways from both results
        pathways1 = result1['netP']['pathways']
        pathways2 = result2['netP']['pathways']
        all_pathways = sorted(set(pathways1).union(set(pathways2)))
        
        # Get actual cell types from results to handle potential mismatches
        if 'net' in result1 and 'dimnames' in result1['net'] and len(result1['net']['dimnames']) > 0:
            cells1 = result1['net']['dimnames'][0]
        else:
            cells1 = []
            
        if 'net' in result2 and 'dimnames' in result2['net'] and len(result2['net']['dimnames']) > 0:
            cells2 = result2['net']['dimnames'][0]
        else:
            cells2 = []
        
        # Use intersection of cell types to avoid index errors
        common_cells = [ct for ct in cell_types if ct in cells1 and ct in cells2]
        
        # Print debug info
        print(f"Cell types from input: {cell_types}")
        print(f"Cell types from result1: {cells1}")
        print(f"Cell types from result2: {cells2}")
        print(f"Common cell types: {common_cells}")
        
        if not common_cells:
            st.warning("No common cell types found between the two results!")
            # Use all provided cell types as fallback
            common_cells = cell_types
        
        # For each pathway, get probability differences
        for pathway in all_pathways:
            pathway_idx1 = list(pathways1).index(pathway) if pathway in pathways1 else None
            pathway_idx2 = list(pathways2).index(pathway) if pathway in pathways2 else None
            
            # Get the indices within each result
            for source in common_cells:
                for target in common_cells:
                    # Get source and target indices in each result
                    if source in cells1 and target in cells1:
                        i1 = cells1.index(source)
                        j1 = cells1.index(target)
                    else:
                        i1 = j1 = None
                        
                    if source in cells2 and target in cells2:
                        i2 = cells2.index(source)
                        j2 = cells2.index(target)
                    else:
                        i2 = j2 = None
                    
                    # Get probability values
                    prob1 = (result1['netP']['prob'][i1, j1, pathway_idx1] 
                             if i1 is not None and j1 is not None and pathway_idx1 is not None else 0)
                    prob2 = (result2['netP']['prob'][i2, j2, pathway_idx2] 
                             if i2 is not None and j2 is not None and pathway_idx2 is not None else 0)
                    
                    prob_diff.append({
                        'source': source,
                        'target': target,
                        'pathway': pathway,
                        'actual_prob_diff': prob2 - prob1  # Group 2 - Group 1
                    })
    
    # Calculate weight differences
    if ('network' in result1 and 'strength_matrix' in result1['network'] and 
        'network' in result2 and 'strength_matrix' in result2['network']):
        
        matrix1 = result1['network']['strength_matrix']
        matrix2 = result2['network']['strength_matrix']
        
        # Use intersection of cell types
        common_cells = [ct for ct in cell_types if ct in matrix1.index and ct in matrix1.columns
                      and ct in matrix2.index and ct in matrix2.columns]
                      
        if not common_cells:
            st.warning("No common cell types found in strength matrices!")
            # Fallback to provided cell types with checks
            for source in cell_types:
                for target in cell_types:
                    if (source in matrix1.index and target in matrix1.columns and
                        source in matrix2.index and target in matrix2.columns):
                        
                        weight1 = matrix1.loc[source, target]
                        weight2 = matrix2.loc[source, target]
                        
                        weight_diff.append({
                            'source': source,
                            'target': target,
                            'actual_weight_diff': weight2 - weight1  # Group 2 - Group 1
                        })
        else:
            # Use common cells
            for source in common_cells:
                for target in common_cells:
                    weight1 = matrix1.loc[source, target]
                    weight2 = matrix2.loc[source, target]
                    
                    weight_diff.append({
                        'source': source,
                        'target': target,
                        'actual_weight_diff': weight2 - weight1  # Group 2 - Group 1
                    })
    
    return {
        'prob_diff': pd.DataFrame(prob_diff),
        'weight_diff': pd.DataFrame(weight_diff)
    }




def calculate_p_values_vectorized(actual_diff, perm_results):
    """
    行列演算を使用して p 値を計算する最適化版 (元の関数と同じ結果を保証)
    
    Parameters
    ----------
    actual_diff : dict
        実際の差分を含む辞書
    perm_results : list
        順列結果のリスト
        
    Returns
    -------
    dict
        確率と重みの p 値を含む辞書
    """
    # 準備: 順列数を取得
    n_perms = len(perm_results)
    if n_perms == 0:
        return {'p_values_prob': actual_diff['prob_diff'].copy(), 
                'p_values_weight': actual_diff['weight_diff'].copy()}
    
    # 確率差のp値計算
    prob_p_values = actual_diff['prob_diff'].copy()
    weight_p_values = actual_diff['weight_diff'].copy()
    
    # カウント用の辞書を初期化 (元の関数と同様に)
    prob_counts = {}
    for _, row in prob_p_values.iterrows():
        key = (row['source'], row['target'], row['pathway'])
        prob_counts[key] = {
            'greater_equal': 0,
            'less_equal': 0
        }
    
    weight_counts = {}
    for _, row in weight_p_values.iterrows():
        key = (row['source'], row['target'])
        weight_counts[key] = {
            'greater_equal': 0,
            'less_equal': 0
        }
    
    # 順列結果のカウント (行列演算で高速化)
    # 1. 各キーごとに順列値を収集
    prob_keys = list(prob_counts.keys())
    weight_keys = list(weight_counts.keys())
    
    # 順列結果から値を抽出し、行列形式に変換
    prob_perm_values = {key: np.zeros(n_perms) for key in prob_keys}
    weight_perm_values = {key: np.zeros(n_perms) for key in weight_keys}
    
    # 各順列結果から値を抽出
    for p_idx, perm_result in enumerate(perm_results):
        # 確率差
        for _, row in perm_result['prob_diff'].iterrows():
            key = (row['source'], row['target'], row['pathway'])
            if key in prob_perm_values:
                prob_perm_values[key][p_idx] = row['actual_prob_diff']
        
        # 重み差
        for _, row in perm_result['weight_diff'].iterrows():
            key = (row['source'], row['target'])
            if key in weight_perm_values:
                weight_perm_values[key][p_idx] = row['actual_weight_diff']
 
    # 2. 各キーについて条件を満たす順列値をカウント (行列演算)
    for key in prob_keys:
        perm_values = prob_perm_values[key]
        # 実際の差分値を取得
        actual_diff_val = next((row['actual_prob_diff'] for _, row in prob_p_values.iterrows() 
                            if (row['source'], row['target'], row['pathway']) == key), 0)
        
        # 実際の差分値との比較に変更
        if actual_diff_val >= 0:
            prob_counts[key]['greater_equal'] = np.sum(perm_values >= actual_diff_val)
        else:
            prob_counts[key]['less_equal'] = np.sum(perm_values <= actual_diff_val)
    
    # 重み差についても同様に修正
    for key in weight_keys:
        perm_values = weight_perm_values[key]
        actual_diff_val = next((row['actual_weight_diff'] for _, row in weight_p_values.iterrows() 
                            if (row['source'], row['target']) == key), 0)
        
        if actual_diff_val >= 0:
            weight_counts[key]['greater_equal'] = np.sum(perm_values >= actual_diff_val)
        else:
            weight_counts[key]['less_equal'] = np.sum(perm_values <= actual_diff_val)
    
    
    # 3. p値の計算 (元の関数と完全に同じロジック)
    # 確率のp値
    for i, row in prob_p_values.iterrows():
        key = (row['source'], row['target'], row['pathway'])
        actual_diff_val = row['actual_prob_diff']
        
        if actual_diff_val >= 0:
            p_val = prob_counts[key]['greater_equal'] / n_perms
        else:
            p_val = prob_counts[key]['less_equal'] / n_perms
        
        prob_p_values.at[i, 'p_value_prob'] = p_val
    
    # 重みのp値
    for i, row in weight_p_values.iterrows():
        key = (row['source'], row['target'])
        actual_diff_val = row['actual_weight_diff']
        
        if actual_diff_val >= 0:
            p_val = weight_counts[key]['greater_equal'] / n_perms
        else:
            p_val = weight_counts[key]['less_equal'] / n_perms
        
        weight_p_values.at[i, 'p_value_weight'] = p_val
    
    return {
        'p_values_prob': prob_p_values,
        'p_values_weight': weight_p_values
    }

def calculate_p_values(actual_diff, perm_results):
    """
    Calculate p-values based on permutation results
    
    Parameters
    ----------
    actual_diff : dict
        Dictionary containing actual differences
    perm_results : list
        List of permutation result dictionaries
        
    Returns
    -------
    dict
        Dictionary containing p-values for probability and weight differences
    """
    # Initialize DataFrames to store p-values
    prob_p_values = actual_diff['prob_diff'].copy()
    weight_p_values = actual_diff['weight_diff'].copy()
    
    # Add columns for p-values
    prob_p_values['p_value_prob'] = 1.0
    weight_p_values['p_value_weight'] = 1.0
    
    # Initialize counting dictionaries for calculating p-values
    prob_counts = {}
    for _, row in prob_p_values.iterrows():
        key = (row['source'], row['target'], row['pathway'])
        prob_counts[key] = {
            'greater_equal': 0,
            'less_equal': 0
        }
    
    weight_counts = {}
    for _, row in weight_p_values.iterrows():
        key = (row['source'], row['target'])
        weight_counts[key] = {
            'greater_equal': 0,
            'less_equal': 0
        }
    
    # Count permutation results
    for perm_result in perm_results:
        # Count for probability
        for _, row in perm_result['prob_diff'].iterrows():
            key = (row['source'], row['target'], row['pathway'])
            if key in prob_counts:
                perm_diff = row['actual_prob_diff']
                prob_counts[key]['greater_equal'] += (perm_diff >= 0)
                prob_counts[key]['less_equal'] += (perm_diff <= 0)
        
        # Count for weight
        for _, row in perm_result['weight_diff'].iterrows():
            key = (row['source'], row['target'])
            if key in weight_counts:
                perm_diff = row['actual_weight_diff']
                weight_counts[key]['greater_equal'] += (perm_diff >= 0)
                weight_counts[key]['less_equal'] += (perm_diff <= 0)
    
    # Calculate p-values
    n_perms = len(perm_results)
    
    # For probability
    for i, row in prob_p_values.iterrows():
        key = (row['source'], row['target'], row['pathway'])
        actual_diff_val = row['actual_prob_diff']
        
        if actual_diff_val >= 0:
            p_val = prob_counts[key]['greater_equal'] / n_perms
        else:
            p_val = prob_counts[key]['less_equal'] / n_perms
        
        prob_p_values.at[i, 'p_value_prob'] = p_val
    
    # For weight
    for i, row in weight_p_values.iterrows():
        key = (row['source'], row['target'])
        actual_diff_val = row['actual_weight_diff']
        
        if actual_diff_val >= 0:
            p_val = weight_counts[key]['greater_equal'] / n_perms
        else:
            p_val = weight_counts[key]['less_equal'] / n_perms
        
        weight_p_values.at[i, 'p_value_weight'] = p_val
    
    return {
        'p_values_prob': prob_p_values,
        'p_values_weight': weight_p_values
    }

def adjust_p_values(p_values):
    """
    Adjust p-values for multiple testing
    
    Parameters
    ----------
    p_values : dict
        Dictionary containing p-values
        
    Returns
    -------
    dict
        Dictionary containing adjusted p-values
    """

    # Add debug info
    print("p_values keys:", p_values.keys())
    print("p_values_prob shape:", p_values['p_values_prob'].shape)
    print("p_values_prob head:", p_values['p_values_prob'].head())
    # エラー処理を追加
    if 'p_value_prob' in p_values['p_values_prob'].columns:
        print("p_value_prob unique values:", p_values['p_values_prob']['p_value_prob'].unique())
    else:
        print("Warning: p_value_prob column not found. Available columns:", p_values['p_values_prob'].columns.tolist())
        # p_value_prob列が存在しない場合、デフォルト値で作成
        p_values['p_values_prob']['p_value_prob'] = 1.0

    # Adjust probability p-values
    prob_adjusted = p_values['p_values_prob'].copy()
    
    # Check if p_values are valid and not all the same value
    if prob_adjusted.empty or 'p_value_prob' not in prob_adjusted.columns or prob_adjusted['p_value_prob'].nunique() <= 1:
        # Handle the case where there's no valid p-values to adjust
        if 'p_value_prob' not in prob_adjusted.columns:
            prob_adjusted['p_value_prob'] = 1.0
        prob_adjusted['p_adj_bonferroni_prob'] = prob_adjusted['p_value_prob']
        prob_adjusted['p_adj_BH_prob'] = prob_adjusted['p_value_prob']
    else:
        prob_adjusted['p_adj_bonferroni_prob'] = multipletests(
            prob_adjusted['p_value_prob'], method='bonferroni'
        )[1]
        prob_adjusted['p_adj_BH_prob'] = multipletests(
            prob_adjusted['p_value_prob'], method='fdr_bh'
        )[1]
    
    # Do the same for weight p-values
    weight_adjusted = p_values['p_values_weight'].copy()
    
    # エラー処理を追加: p_value_weight列が存在しない場合
    if 'p_value_weight' not in weight_adjusted.columns:
        print("Warning: p_value_weight column not found. Creating with default values.")
        weight_adjusted['p_value_weight'] = 1.0
    
    if weight_adjusted.empty or weight_adjusted['p_value_weight'].nunique() <= 1:
        weight_adjusted['p_adj_bonferroni_weight'] = weight_adjusted['p_value_weight']
        weight_adjusted['p_adj_BH_weight'] = weight_adjusted['p_value_weight']
    else:
        weight_adjusted['p_adj_bonferroni_weight'] = multipletests(
            weight_adjusted['p_value_weight'], method='bonferroni'
        )[1]
        weight_adjusted['p_adj_BH_weight'] = multipletests(
            weight_adjusted['p_value_weight'], method='fdr_bh'
        )[1]
    
    return {
        'p_values_prob_adjusted': prob_adjusted,
        'p_values_weight_adjusted': weight_adjusted
    }

def plot_p_value_distribution(p_values_adjusted):
    """
    Plot the distribution of p-values
    
    Parameters
    ----------
    p_values_adjusted : dict
        Dictionary containing adjusted p-values
        
    Returns
    -------
    tuple
        Tuple of matplotlib figures
    """
    # Probability p-values
    fig_prob, ax_prob = plt.subplots(figsize=(10, 6))
    sns.histplot(p_values_adjusted['p_values_prob_adjusted']['p_value_prob'], 
                bins=50, color='blue', alpha=0.7, ax=ax_prob)
    ax_prob.set_title('Distribution of P-values (Probability)')
    ax_prob.set_xlabel('P-value')
    ax_prob.set_ylabel('Count')
    
    # Weight p-values
    fig_weight, ax_weight = plt.subplots(figsize=(10, 6))
    sns.histplot(p_values_adjusted['p_values_weight_adjusted']['p_value_weight'], 
                bins=50, color='red', alpha=0.7, ax=ax_weight)
    ax_weight.set_title('Distribution of P-values (Weight)')
    ax_weight.set_xlabel('P-value')
    ax_weight.set_ylabel('Count')
    
    return fig_prob, fig_weight

def filter_significant_results(p_values_adjusted, threshold=0.05):
    """
    Filter significant results based on adjusted p-values
    
    Parameters
    ----------
    p_values_adjusted : dict
        Dictionary containing adjusted p-values
    threshold : float
        P-value threshold
        
    Returns
    -------
    dict
        Dictionary containing significant results
    """
    # Filter probability results
    significant_prob = p_values_adjusted['p_values_prob_adjusted'][
        p_values_adjusted['p_values_prob_adjusted']['p_adj_BH_prob'] < threshold
    ].sort_values('p_adj_BH_prob')
    
    # Filter weight results
    significant_weight = p_values_adjusted['p_values_weight_adjusted'][
        p_values_adjusted['p_values_weight_adjusted']['p_adj_BH_weight'] < threshold
    ].sort_values('p_adj_BH_weight')
    
    return {
        'prob': significant_prob,
        'weight': significant_weight
    }

def create_refined_heatmap(data, diff_var, p_value_var, title=None, adj_p_value_var=None, use_adjusted_p=False, 
                          n_perm=1000, is_pathway=False, cell_types=None, cell_order=None, reverse_diff=False):
    """
    Create a refined heatmap showing differences and p-values with improved visualization
    
    Parameters
    ----------
    data : pd.DataFrame
        DataFrame containing difference and p-value data
    diff_var : str
        Column name for difference values
    p_value_var : str
        Column name for nominal p-values
    adj_p_value_var : str, optional
        Column name for adjusted p-values
    use_adjusted_p : bool, optional
        Whether to use adjusted p-values for dot size
    title : str
        Plot title
    n_perm : int
        Number of permutations performed (default=1000)
    is_pathway : bool
        Whether the data is pathway-specific
    cell_types : list
        List of cell types to include
    cell_order : list
        List specifying the order of cell types
    reverse_diff : bool
        Whether to reverse the direction of the difference values
        
    Returns
    -------
    matplotlib.figure.Figure
        The heatmap figure
    """
    # Create a copy of the data
    plot_data = data.copy()
    
    # Apply reverse_diff if needed with improved error checking
    if reverse_diff:
        # Ensure the column is numeric before reversing
        if pd.api.types.is_numeric_dtype(plot_data[diff_var]):
            plot_data[diff_var] = plot_data[diff_var].multiply(-1)
            print(f"Values in {diff_var} have been reversed successfully")
        else:
            print(f"Warning: Column {diff_var} is not numeric. Cannot reverse values.")
    
    # Debug info
    print(f"Data shape: {plot_data.shape}")
    print(f"P-value column: {p_value_var}")
    if adj_p_value_var:
        print(f"Adjusted P-value column: {adj_p_value_var}")
    print(f"Using adjusted P-values: {use_adjusted_p}")
    print(f"Difference direction reversed: {reverse_diff}")
    
    # Determine which p-value column to use
    p_value_to_use = adj_p_value_var if use_adjusted_p and adj_p_value_var else p_value_var
    
    # Determine the minimum achievable p-value based on permutation count
    if n_perm <= 10:
        min_p_value = 0.1
        max_log_p = 1  # -log10(0.1) = 1
        p_thresholds = [0.1]  # Only show 0.1 in legend
    elif n_perm <= 100:
        min_p_value = 0.01
        max_log_p = 2  # -log10(0.01) = 2
        p_thresholds = [0.1, 0.01]  # Show 0.1 and 0.01 in legend
    else:
        min_p_value = 0.001
        max_log_p = 3  # -log10(0.001) = 3
        p_thresholds = [0.1, 0.01, 0.001]  # Show all thresholds in legend
    
    print(f"Minimum p-value for {n_perm} permutations: {min_p_value}")
    print(f"P-value thresholds to display in legend: {p_thresholds}")
    
    # Use cell_order if provided, otherwise use cell_types
    if cell_order is not None:
        ordered_cells = cell_order
    else:
        ordered_cells = cell_types
    
    # Ensure we have all source and target combinations
    if cell_types is not None:
        # Create a grid of all source and target combinations
        all_combinations = []
        for source in ordered_cells:
            for target in ordered_cells:
                if source in cell_types and target in cell_types:
                    all_combinations.append({
                        'source': source,
                        'target': target
                    })
        
        all_combinations_df = pd.DataFrame(all_combinations)
        
        # Merge with the actual data
        plot_data = pd.merge(
            all_combinations_df, plot_data, on=['source', 'target'], how='left'
        )
        
        # Fill NA values
        plot_data[diff_var] = plot_data[diff_var].fillna(0)
        plot_data[p_value_to_use] = plot_data[p_value_to_use].fillna(1)
    
    # Create a pivot table for the heatmap
    pivot_data = plot_data.pivot_table(
        index='source', columns='target', values=diff_var, fill_value=0
    )
    
    # Create a pivot table for p-values
    pivot_pvals = plot_data.pivot_table(
        index='source', columns='target', values=p_value_to_use, fill_value=1
    )
    
    # Reorder rows and columns if custom order is provided
    if cell_order is not None:
        # Filter to include only cells that exist in the pivot table
        valid_ordered_cells = [cell for cell in ordered_cells if cell in pivot_data.index]
        
        # Reorder the rows and columns
        pivot_data = pivot_data.reindex(index=valid_ordered_cells, columns=valid_ordered_cells)
        pivot_pvals = pivot_pvals.reindex(index=valid_ordered_cells, columns=valid_ordered_cells)
    
    # Set up the figure with additional bottom space for legend
    fig, ax = plt.subplots(figsize=(10, 9))  # Increased height to accommodate legend
    
    # Calculate max absolute value for symmetric color scale
    vmax = np.max(np.abs(pivot_data.values))
    vmin = -vmax
    
    # Define custom colormap: teal for negative, white for zero, red for positive
    cmap = sns.diverging_palette(200, 10, s=99, l=55, sep=10, as_cmap=True)
    
    # Create the heatmap with light grey background for the grid
    ax.set_facecolor('#f5f5f5')  # Light grey background
    
    # Create the heatmap with a smaller colorbar
    cbar_kws = {
        'label': 'Difference',
        'shrink': 0.4,     # Make colorbar 40% of original height
        'aspect': 15,      # Lower aspect ratio makes colorbar wider
        'pad': 0.02        # Less space between heatmap and colorbar
    }
    
    # Generate the heatmap
    heatmap = sns.heatmap(
        pivot_data, 
        cmap=cmap, 
        center=0,
        vmin=vmin, 
        vmax=vmax,
        annot=False, 
        fmt='.2f',
        linewidths=0.5,
        linecolor='white',
        ax=ax,
        cbar_kws=cbar_kws
    )
    
    # Base size for dots (to indicate significance)
    base_size = 20  # Large base size for visibility
    
    # Add dots with continuous scaling based on p-values
    for i in range(pivot_data.shape[0]):
        for j in range(pivot_data.shape[1]):
            p_val = pivot_pvals.iloc[i, j]
            
            # Skip if invalid p-value
            if p_val >= 1 or pd.isna(p_val):
                continue
            
            # Calculate -log10(p) value capped at max_log_p
            log_p_val = min(-np.log10(p_val), max_log_p)
            
            # Continuous scaling from 0 to max_log_p, mapped to size 0.25 to 1.0
            # This ensures the smallest p-value gets the largest dot
            if log_p_val >= 0:
                # Normalize to range [0, 1] based on max_log_p
                normalized = log_p_val / max_log_p
                
                # Scale to range [0.25, 1.0] for dot sizes
                scaled_size = 0.25 + 0.75 * normalized
                
                # Calculate final dot size
                final_size = base_size * scaled_size
                
                # Add a black dot to indicate significance
                ax.plot(j + 0.5, i + 0.5, 'o', markersize=final_size,
                        color='black', alpha=1.0)
    
    # Rotate x-axis labels for better readability
    plt.xticks(rotation=90, ha='right')
    
    # Add axis labels
    ax.set_xlabel('Target Cell Type')
    ax.set_ylabel('Source Cell Type')
    
    # Add p-value type to title if using adjusted p-values
    if use_adjusted_p:
        title = f"{title} (BH-adjusted p-values)"
    else:
        title = f"{title} (nominal p-values)"
    
    # Set title
    ax.set_title(title, fontsize=14)
    
    # Create a legend for p-value thresholds
    from matplotlib.lines import Line2D
    legend_elements = []
    
    # Create legend elements for the p-value thresholds
    for p in p_thresholds:
        # Calculate -log10(p)
        log_p = -np.log10(p)
        
        # Normalize to range [0, 1] based on max_log_p
        normalized = log_p / max_log_p
        
        # Scale to range [0.25, 1.0] for dot sizes
        scaled_size = 0.25 + 0.75 * normalized
        
        # Calculate legend size
        legend_size = base_size * scaled_size
        
        legend_elements.append(
            Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                   markersize=legend_size, label=f'{p}')
        )
    
    # Add the legend below the heatmap in a horizontal layout
    p_legend = ax.legend(
        handles=legend_elements, 
        title="p-value", 
        loc='upper center',        # Position at upper center
        bbox_to_anchor=(0.5, -0.15),  # Position below the heatmap
        fontsize=9,                # Increased font size
        title_fontsize=10,         # Increased title font size
        frameon=True,
        framealpha=0.8,     
        ncol=len(p_thresholds),    # Horizontal layout with one column per p-value
        borderpad=0.3,      
        handletextpad=0.5,  
        columnspacing=1.5,         # More space between columns
        labelspacing=0.3    
    )
    
    # Make room for the colorbar and legends at bottom
    fig.tight_layout(rect=[0, 0.05, 0.88, 0.95])
    
    return fig


def suggest_optimal_n_jobs():
    """
    Suggest optimal number of jobs based on available CPU cores
    
    Returns
    -------
    int
        Suggested optimal number of jobs
    """
    total_cores = os.cpu_count()
    
    # For less than 8 cores, use all but 1
    if total_cores <= 8:
        return max(1, total_cores - 1)
    
    # For 8-32 cores, use half of available cores
    if total_cores <= 32:
        return max(8, total_cores // 2)
    
    # For more than 32 cores, use 16 cores
    return 16

def display_visualizations(perm_results, cell_types, group1, group2, significance_threshold=0.05):
    """
    Display visualizations based on permutation results
    
    Parameters
    ----------
    perm_results : dict
        Permutation analysis results
    cell_types : list
        List of cell types
    group1 : str
        Name of first group
    group2 : str
        Name of second group
    significance_threshold : float
        P-value threshold for significance
    
    Returns
    -------
    None
    """
    # Get number of permutations from the results (if available, otherwise use default)
    n_perms = 1000  # Default value
    if 'permutation_results' in perm_results:
        n_perms = len(perm_results['permutation_results'])
        st.write(f"Analysis based on {n_perms} permutations")
    
    # Create tabs for visualization
    results_tabs = st.tabs([
        "P-value Distribution", 
        "Significant Results",
        "Weight Heatmap",
        "Pathway Heatmaps"
    ])
    
    # Tab 1: P-value Distribution
    with results_tabs[0]:
        st.header("P-value Distribution")
        
        col1, col2 = st.columns(2)
        
        # Plot the p-value distribution
        fig_prob, fig_weight = plot_p_value_distribution(perm_results['adjusted_p_values'])
        
        with col1:
            st.pyplot(fig_prob)
            
            # Save and download
            prob_pdf_path = f"{cellchat_temp_dir}/pval_distribution_prob.pdf"
            fig_prob.savefig(prob_pdf_path, bbox_inches='tight')
            
            with open(prob_pdf_path, "rb") as pdf_file:
                st.download_button(
                    label="Download Probability P-value Distribution",
                    data=pdf_file.read(),
                    file_name="pval_distribution_prob.pdf",
                    mime="application/pdf"
                )
        
        with col2:
            st.pyplot(fig_weight)
            
            # Save and download
            weight_pdf_path = f"{cellchat_temp_dir}/pval_distribution_weight.pdf"
            fig_weight.savefig(weight_pdf_path, bbox_inches='tight')
            
            with open(weight_pdf_path, "rb") as pdf_file:
                st.download_button(
                    label="Download Weight P-value Distribution",
                    data=pdf_file.read(),
                    file_name="pval_distribution_weight.pdf",
                    mime="application/pdf"
                )
    
    # Tab 2: Significant Results
    with results_tabs[1]:
        st.header("Significant Results")
        
        # Significance threshold
        local_threshold = st.slider(
            "Significance threshold (adjusted p-value):",
            min_value=0.001,
            max_value=0.1,
            value=significance_threshold,
            step=0.001,
            key="significance_threshold_slider"
        )
        
        # Filter significant results
        significant_results = filter_significant_results(
            perm_results['adjusted_p_values'],
            threshold=local_threshold
        )
        
        # Display results in tabs
        sig_tabs = st.tabs(["Probability", "Weight"])
        
        with sig_tabs[0]:
            st.subheader(f"Significant Probability Differences (BH-adjusted p < {local_threshold})")
            
            if len(significant_results['prob']) > 0:
                st.write(f"Found {len(significant_results['prob'])} significant probability differences")
                
                # Display the top significant results
                st.dataframe(
                    significant_results['prob'][['source', 'target', 'pathway', 
                                              'actual_prob_diff', 'p_value_prob', 
                                              'p_adj_BH_prob']]
                )
                
                # Download as CSV
                csv = significant_results['prob'].to_csv(index=False)
                st.download_button(
                    label="Download significant probability results as CSV",
                    data=csv,
                    file_name="significant_probability_differences.csv",
                    mime="text/csv"
                )
            else:
                st.write("No significant probability differences found at this threshold.")
        
        with sig_tabs[1]:
            st.subheader(f"Significant Weight Differences (BH-adjusted p < {local_threshold})")
            
            if len(significant_results['weight']) > 0:
                st.write(f"Found {len(significant_results['weight'])} significant weight differences")
                
                # Display the top significant results
                st.dataframe(
                    significant_results['weight'][['source', 'target', 
                                               'actual_weight_diff', 'p_value_weight', 
                                               'p_adj_BH_weight']]
                )
                
                # Download as CSV
                csv = significant_results['weight'].to_csv(index=False)
                st.download_button(
                    label="Download significant weight results as CSV",
                    data=csv,
                    file_name="significant_weight_differences.csv",
                    mime="text/csv"
                )
            else:
                st.write("No significant weight differences found at this threshold.")
    
    # Tab 3: Weight Heatmap

    # Weight Heatmap タブ内でも同様の修正を行う
    with results_tabs[2]:
        st.header("Weight Difference Heatmap")
        
        # BH調整済みP値を使用するかどうかのチェックボックス
        use_adjusted_p_weight = st.checkbox("BH調整済みP値を使用", value=False, key="use_adjusted_p_weight")
        
        # 使用するP値列を設定
        p_val_col_weight = 'p_value_weight'
        adj_p_val_col_weight = 'p_adj_BH_weight'

        reverse_direction = st.checkbox("Switch control group?", value=False, 
                               help="Toggle to view Group2 vs Group1 or Group1 vs Group2")

     
        # Create the heatmap with custom cell order if available
        # Update title based on reverse_direction
        if reverse_direction:
            heatmap_title = f"Weight Differences: {group1} vs {group2}"
        else:
            heatmap_title = f"Weight Differences: {group2} vs {group1}"
            
        fig_weight_heatmap = create_refined_heatmap(
            perm_results['adjusted_p_values']['p_values_weight_adjusted'],
            'actual_weight_diff',
            p_val_col_weight,
            title=heatmap_title,
            adj_p_value_var=adj_p_val_col_weight,
            use_adjusted_p=use_adjusted_p_weight,
            cell_types=cell_types,
            cell_order=st.session_state.cell_type_order,
            n_perm=n_perms,
            reverse_diff=reverse_direction
        )
        
        st.pyplot(fig_weight_heatmap)
        
        # Save and download
        weight_heatmap_path = f"{cellchat_temp_dir}/weight_difference_heatmap{'_adjP' if use_adjusted_p_weight else ''}.pdf"
        fig_weight_heatmap.savefig(weight_heatmap_path, bbox_inches='tight')
        
        with open(weight_heatmap_path, "rb") as pdf_file:
            st.download_button(
                label="Download Weight Difference Heatmap",
                data=pdf_file.read(),
                file_name=f"weight_difference_heatmap{'_adjP' if use_adjusted_p_weight else ''}.pdf",
                mime="application/pdf"
            )

    
    with results_tabs[3]:
        st.header("Pathway-Specific Heatmaps")
        
        # Get all available pathways
        if 'result1' in perm_results and 'netP' in perm_results['result1'] and 'pathways' in perm_results['result1']['netP']:
            pathways1 = set(perm_results['result1']['netP']['pathways'])
        else:
            pathways1 = set()
            
        if 'result2' in perm_results and 'netP' in perm_results['result2'] and 'pathways' in perm_results['result2']['netP']:
            pathways2 = set(perm_results['result2']['netP']['pathways'])
        else:
            pathways2 = set()
            
        all_pathways = sorted(pathways1.union(pathways2))
        
        if not all_pathways:
            st.warning("No pathways available for heatmap visualization.")
        else:
            # Initialize session state for pathway selection if not exists
            if 'selected_pathway' not in st.session_state:
                st.session_state.selected_pathway = all_pathways[0]
                
            # P値表示切り替えのためのチェックボックスを追加
            use_adjusted_p = st.checkbox("BH調整済みP値を使用", value=False, key="use_adjusted_p_pathway")
                
            # Pathway selection with callback function to update session state
            def on_pathway_change():
                # This function is called when the selectbox value changes
                # The new value is already in st.session_state.pathway_selector
                st.session_state.selected_pathway = st.session_state.pathway_selector
                
            selected_pathway = st.selectbox(
                "Select pathway for heatmap:",
                all_pathways,
                key="pathway_selector",
                index=all_pathways.index(st.session_state.selected_pathway) if st.session_state.selected_pathway in all_pathways else 0,
                on_change=on_pathway_change
            )
            
            # Use the pathway from session state
            current_pathway = st.session_state.selected_pathway
            
            # Filter data for the selected pathway
            pathway_data = perm_results['adjusted_p_values']['p_values_prob_adjusted'][
                perm_results['adjusted_p_values']['p_values_prob_adjusted']['pathway'] == current_pathway
            ]
            
            # 使用するP値列を設定
            p_val_col = 'p_value_prob'
            adj_p_val_col = 'p_adj_BH_prob'
            reverse_direction = st.checkbox("Change control group", value=False, 
                               help="Toggle to view Group2 vs Group1 or Group1 vs Group2")

            # Create the heatmap with custom cell order if available
            # Update title based on reverse_direction
            if reverse_direction:
                pathway_title = f"Probability Differences - {current_pathway}: {group1} vs {group2}"
            else:
                pathway_title = f"Probability Differences - {current_pathway}: {group2} vs {group1}"
                
            fig_pathway = create_refined_heatmap(
                pathway_data,
                'actual_prob_diff',
                p_val_col,
                title=pathway_title,
                adj_p_value_var=adj_p_val_col,
                use_adjusted_p=use_adjusted_p,  
                is_pathway=True,
                cell_types=cell_types,
                cell_order=st.session_state.cell_type_order,
                n_perm=n_perms,
                reverse_diff=reverse_direction
            )
            
            st.pyplot(fig_pathway)
            
            # Save and download
            pathway_pdf_path = f"{cellchat_temp_dir}/pathway_{current_pathway}_heatmap{'_adjP' if use_adjusted_p else ''}.pdf"
            fig_pathway.savefig(pathway_pdf_path, bbox_inches='tight')
            
            with open(pathway_pdf_path, "rb") as pdf_file:
                st.download_button(
                    label=f"Download {current_pathway} Pathway Heatmap",
                    data=pdf_file.read(),
                    file_name=f"pathway_{current_pathway}_heatmap{'_adjP' if use_adjusted_p else ''}.pdf",
                    mime="application/pdf"
                )


if __name__ == "__main__":
    # Main title
    st.title("CellChat Permutation Comparison Analysis")

    # Initialize session state variables if they don't exist
    if 'analysis_parameters' not in st.session_state:
        st.session_state.analysis_parameters = {}
    if 'analysis_completed' not in st.session_state:
        st.session_state.analysis_completed = False
    if 'current_mode' not in st.session_state:
        st.session_state.current_mode = "setup"  # "setup", "running", or "visualization"
    if 'cell_type_order' not in st.session_state:
        st.session_state.cell_type_order = None

    # Create sidebar
    with st.sidebar:
        st.title("Settings")
        
        # Cell type order settings (only shown when cell types are available)
        if 'cell_types' in st.session_state:
            st.header("Cell Type Display Settings")
            sort_cell = st.checkbox("Change cell type order?", key="sort_cell_checkbox")
            
            if sort_cell:
                with st.form("Sorter"):
                    sorted_cell_types = sort_items(sorted(st.session_state.cell_types))
                    submitted_sort = st.form_submit_button("Done sorting")
                st.session_state.cell_type_order = sorted_cell_types
                    
            else:
                # Reset to default order
                st.session_state.cell_type_order = None
    
    # Main panel
    st.header("Data Settings")

    uploaded_file = st.file_uploader("Upload file", type=['h5ad', 'pkl'], 
                                 help="AnnData file for CellChat analysis or saved permutation results (.pkl)")
    if uploaded_file is not None:
        # Check file type and handle accordingly
        if uploaded_file.name.endswith('.pkl'):
            # Load pickle file with saved permutation results
            try:
                with st.spinner('Loading permutation results...'):
                    # Load the pickle data
                    pickled_data = pickle.load(uploaded_file)
                    
                    # Validate that the pickle contains all necessary information
                    required_keys = ['group1', 'group2', 'actual_diff', 'adjusted_p_values']
                    missing_keys = [key for key in required_keys if key not in pickled_data]
                    
                    if missing_keys:
                        st.error(f"Invalid permutation results file. Missing required data: {', '.join(missing_keys)}")
                        st.stop()
                    
                    # Extract data from the pickle file
                    group1 = pickled_data['group1']
                    group2 = pickled_data['group2']
                    perm_results = {
                        'actual_diff': pickled_data['actual_diff'],
                        'adjusted_p_values': pickled_data['adjusted_p_values']
                    }
                    
                    # If the pickle includes cell types, use them
                    if 'cell_types' in pickled_data:
                        cell_types = pickled_data['cell_types']
                    # Otherwise try to extract from actual_diff if possible
                    elif 'prob_diff' in pickled_data['actual_diff'] and len(pickled_data['actual_diff']['prob_diff']) > 0:
                        # Extract unique cell types from the source and target columns
                        cell_types_from_source = pickled_data['actual_diff']['prob_diff']['source'].unique()
                        cell_types_from_target = pickled_data['actual_diff']['prob_diff']['target'].unique()
                        cell_types = sorted(set(cell_types_from_source) | set(cell_types_from_target))
                    else:
                        st.error("Could not determine cell types from the pickle file.")
                        st.stop()
                    
                    # Additional data if available
                    if 'result1' in pickled_data:
                        perm_results['result1'] = pickled_data['result1']
                    if 'result2' in pickled_data:
                        perm_results['result2'] = pickled_data['result2']
                    if 'permutation_results' in pickled_data:
                        perm_results['permutation_results'] = pickled_data['permutation_results']
                    
                    # Get number of permutations if available
                    n_perms = pickled_data.get('n_perms', 
                              pickled_data.get('parameters', {}).get('n_perms', 100))
                    
                    # Set session state variables for visualization
                    st.session_state.perm_results = perm_results
                    st.session_state.cell_types = cell_types
                    st.session_state.group1 = group1
                    st.session_state.group2 = group2
                    
                    # Set analysis as completed and mode to visualization
                    st.session_state.analysis_completed = True
                    st.session_state.current_mode = "visualization"
                    st.session_state.show_visualization = True

                    has_pathway_info = True

                    
                    # Display info about loaded data
                    st.success(f"Successfully loaded permutation results: {group1} vs {group2}")
                    st.write(f"Cell types: {', '.join(cell_types[:5])}{'...' if len(cell_types) > 5 else ''}")
                    st.write(f"Number of permutations: {n_perms}")
                    
                    # Initialize significance threshold if not exists
                    if 'significance_threshold' not in st.session_state:
                        st.session_state.significance_threshold = 0.05
                    
                    # Display visualizations
                    display_visualizations(
                        perm_results=perm_results,
                        cell_types=cell_types,
                        group1=group1,
                        group2=group2,
                        significance_threshold=st.session_state.significance_threshold
                    )
                    
            except Exception as e:
                st.error(f"Error loading pickle file: {str(e)}")
                st.exception(e)
        
        else:
            # Handle h5ad file (existing code)
            try:
                if 'adata' not in st.session_state or st.session_state.get('uploaded_file_name') != uploaded_file.name:
                    with st.spinner('Reading AnnData file...'):
                        adata = sc.read_h5ad(uploaded_file)
                    st.session_state.adata = adata
                    st.session_state.uploaded_file_name = uploaded_file.name
                    # Reset analysis state when a new file is uploaded
                    st.session_state.analysis_completed = False
                    st.session_state.current_mode = "setup"
                    st.session_state.analysis_parameters = {}
                    if 'perm_results' in st.session_state:
                        del st.session_state.perm_results
                else:
                    adata = st.session_state.adata
                
                st.success(f"Successfully loaded file: {uploaded_file.name}")


            
                # Display summary info
                st.write(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
                
                # Show available categorical variables for cell type
                categorical_cols = [col for col in adata.obs.columns if isinstance(adata.obs[col].dtype, pd.CategoricalDtype)]
                other_cols = [col for col in adata.obs.columns if col not in categorical_cols]
                
                # Select cell type annotation column
                default_celltype_idx = find_first_index_or_default(
                    categorical_cols, 
                    ['cell.ident', 'cell_type', 'celltype', 'seurat_clusters'],
                    0
                )
                celltype_key = st.selectbox(
                    "Select cell type annotation column:",
                    categorical_cols,
                    index=min(default_celltype_idx, len(categorical_cols)-1) if categorical_cols else 0,
                    help="Column to use for cell type annotations"
                )
                cell_list = sorted(adata.obs[celltype_key].cat.categories.to_list() if hasattr(adata.obs[celltype_key], 'cat') else sorted(adata.obs[celltype_key].unique()))
                st.write(", ".join(cell_list))

                # Select group comparison column
                default_group_idx = find_first_index_or_default(
                    categorical_cols, 
                    ['orig.ident', 'stimulus', 'condition', 'treatment', 'sample'],
                    0
                )
                splitby_key = st.selectbox(
                    "Select grouping variable for comparison:",
                    categorical_cols,
                    index=min(default_group_idx, len(categorical_cols)-1) if categorical_cols else 0,
                    help="Column to split the dataset for comparison"
                )

                # Get available groups from the selected grouping variable
                if splitby_key:
                    available_groups = adata.obs[splitby_key].cat.categories.tolist() if hasattr(adata.obs[splitby_key], 'cat') else sorted(adata.obs[splitby_key].unique())
                    
                    # If there are exactly 2 groups, select them by default
                    if len(available_groups) == 2:
                        group1 = available_groups[0]
                        group2 = available_groups[1]
                    else:
                        # Otherwise let the user select
                        group1 = st.selectbox("Select first group", available_groups, index=0)
                        remaining_groups = [g for g in available_groups if g != group1]
                        group2 = st.selectbox("Select second group", remaining_groups, index=0)
                    
                    # Display cell counts for selected groups
                    g1_count = (adata.obs[splitby_key] == group1).sum()
                    g2_count = (adata.obs[splitby_key] == group2).sum()
                    st.write(f"Group 1 ({group1}): {g1_count} cells")
                    st.write(f"Group 2 ({group2}): {g2_count} cells")
                
                # Available data layers
                available_layers = list(adata.layers.keys()) if hasattr(adata, 'layers') else []
                st.write(f"Available data layers: {', '.join(available_layers) if available_layers else 'No additional layers'}")
                
                data_layer = st.selectbox(
                    "Select data layer:",
                    ['X (default)'] + available_layers,
                    index=0,
                    help="Data layer to use for analysis"
                )
                
                # Infer species from gene names
                species = st.radio("Species:", ('mouse', 'human'), 
                                  index=check_species_index(adata.var.index.to_list()[:50]))
                
                # Analysis settings 
                st.header("Analysis Settings")
                
                with st.expander("CellChat Analysis Settings", expanded=True):
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        selected_types = st.multiselect(
                            "Signaling types to analyze:",
                            ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling"],
                            default=["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"],
                            help="Types of cellular communication to include in the analysis"
                        )
                        
                        expr_prop = st.number_input(
                            'Min fraction of expressing cells:',
                            min_value=0.0, max_value=0.9, step=0.01, value=0.1,
                            help="Minimum fraction of cells expressing a gene to be considered"
                        )
                        
                        thresh_p = st.number_input(
                            'P-value threshold for differential expression:',
                            min_value=0.0, max_value=0.9, step=0.01, value=0.05,
                            help="Threshold for Wilcoxon test for overexpressed genes"
                        )
                        
                        min_cells = st.number_input(
                            'Minimum cells per group:',
                            min_value=5, max_value=50, step=5, value=10,
                            help="Minimum number of cells required in each group"
                        )

                        population_size = st.checkbox("population.size", value=False, help="細胞数が多いクラスター同士の通信ほど強く（確率が高めに）評価される。")

                        st.markdown("---")
                        fun_type = st.radio("Method to calculate average expression per cell group", ['triMean','truncatedMean'], index = 0, help="推定されるリガンド-受容体ペアの数は、細胞グループごとの平均遺伝子発現量を計算する方法に依存する。trimeanは他の方法よりも少ない相互作用を生成。CellChatはより強い相互作用を予測する性能が優れていることがわかっており、実験的検証のために相互作用を絞り込む上で非常に有用。trimeanは25%切断平均に近似しており、これは1つのグループで発現している細胞の割合が25%未満の場合、平均遺伝子発現量がゼロになることを意味する。truncatedMeanでは5%および10%切断平均等を設定できる")

                        st.markdown("*To reduce stringency for pathway selection, use truncatedMean.*")
                        trim = st.number_input("Trimming for truncated mean:", min_value=0.00, max_value=0.25, value=0.10)

                        st.markdown("---")                    
                        apply_pval_filter = st.checkbox("Filter interaction pathways by P value", value=True)

                        trim_threshold = st.number_input('Filtering P threshold:', 
                                                   min_value=0.01, max_value=0.20, step=0.01, value=0.05)

                       # r_patcher = st.checkbox("Fine-tune to match R calc results", value=False,
                        #    help="Rとの計算結果の違いをある程度吸収する。weightについてはチェックしなくても近似する。")
                        r_patcher = False


                  
                    with col2:
                        # Permutation settings
                        n_perms = st.radio(
                            "Number of permutations:",
                            [10, 50, 100, 500, 1000],
                            index=2,  # Default to 100
                            help="Number of permutations for statistical testing"
                        )
                        
                        do_every_deg = st.checkbox("Calc DEG every premutation step?", value = True,
                            help="permutationごとにDEGを求めて解析する。チェックを外した場合は、最初のDEG解析で求めた遺伝子に関してすべての計算を行う。celchatのpermutationを忠実に行う観点からは毎回DEGを求める。")

                        if do_every_deg:
                            union_mode  = st.radio(
                                "Use individual overexpressed genes, or union or intersection from both groups",
                                ['individual','union','intersection'],
                                index=0,
                                help="Identifies overexpressed genes individually or from both groups and uses their union/intersect for analysis.CellChatのtutorialでは個々にDEG等でフィルタリングしたデータを解析し、それらの遺伝子についてLRの解析を行います。そのため、一方で検討対象とならない遺伝子が生じる。また、個々の解析で有意とならないLRペアはinteraction scoreが0になる。"
                            )
                        else:
                            union_mode  = st.radio(
                                "Use individual overexpressed genes, or union or intersection from both groups",
                                ['union','intersection'],
                                index=0,
                                help="Identifies overexpressed genes individually or from both groups and uses their union/intersect for analysis.CellChatのtutorialでは個々にDEG等でフィルタリングしたデータを解析し、それらの遺伝子についてLRの解析を行います。そのため、一方で検討対象とならない遺伝子が生じます。また、個々の解析で有意とならないLRペアはinteraction scoreが0になります。"
                            )

                        optimal_n_jobs = suggest_optimal_n_jobs()
                     #   st.info(f"For best performance, we recommend using {optimal_n_jobs} CPU cores.")
                        cpu_n = os.cpu_count()
                        
                        n_cpus = st.slider(
                            "Number of CPUs:",
                            min_value=1, max_value=cpu_n, step=1, value=int(cpu_n/8-2),
                            help="Number of CPU cores to use for parallel processing"
                        )

                        use_gpu = False
                        gpu_hybrid = False
                        gpu_precision = False
                     #   gpu_hybrid = st.checkbox("Use GPU?", value=False)
                     #   if gpu_hybrid: #多分不要だが、記載しておく
                     #       use_gpu = True
                     #   else:
                     #       use_gpu = False

                      #  gpu_precision = st.checkbox("float64 calculation?", value=False)

                    
                    # Store current settings
                    current_settings = {
                        'celltype_key': celltype_key,
                        'splitby_key': splitby_key,
                        'group1': group1,
                        'group2': group2,
                        'data_layer': None if data_layer == 'X (default)' else data_layer,
                        'species': species,
                        'selected_types': selected_types,
                        'expr_prop': expr_prop,
                        'thresh_p': thresh_p,
                        'min_cells': min_cells,
                        'n_perms': n_perms,
                        'n_cpus': n_cpus,
                        'cell_list': cell_list
                    }
                    
                    # Add "Run Analysis" button only if needed
                    if not st.session_state.analysis_completed or current_settings != st.session_state.analysis_parameters:
                        # Save current settings for comparison
                        st.session_state.analysis_parameters = current_settings.copy()
                        
                        if st.button("Run Analysis", key="run_analysis"):
                            st.session_state.current_mode = "running"
                            st.rerun()
                
                # Conditional display based on current mode
                if st.session_state.current_mode == "running":
                    # Run the analysis with current parameters
                    try:
                        parameters = st.session_state.analysis_parameters
                        
                        # Extract parameters
                        celltype_key = parameters['celltype_key']
                        splitby_key = parameters['splitby_key']
                        group1 = parameters['group1']
                        group2 = parameters['group2']
                        data_layer = parameters['data_layer']
                        species = parameters['species']
                        selected_types = parameters['selected_types']
                        expr_prop = parameters['expr_prop']
                        thresh_p = parameters['thresh_p']
                        min_cells = parameters['min_cells']
                        n_perms = parameters['n_perms']
                        n_cpus = parameters['n_cpus']
                        cell_list = parameters['cell_list']
                        
                        # Split AnnData into two groups
                        with st.spinner(f"Splitting data into {group1} and {group2} groups..."):
                            adata1, adata2 = split_anndata(adata, splitby_key, group1, group2)
                            
                            # Check if either group has too few cells
                            if adata1.n_obs < min_cells:
                                st.error(f"Group '{group1}' has fewer than {min_cells} cells.")
                                st.stop()
                            if adata2.n_obs < min_cells:
                                st.error(f"Group '{group2}' has fewer than {min_cells} cells.")
                                st.stop()
                            
                            st.success(f"Successfully split data: {group1} ({adata1.n_obs} cells) and {group2} ({adata2.n_obs} cells)")

                        # Get CellChatDB
                        with st.spinner(f"Getting {species} CellChatDB..."):
                            cellchatdb = get_cellchatdb_from_r(species=species)
                            cellchatdb = debug_cellchatdb(cellchatdb)
                            st.success(f"{species} CellChatDB successfully obtained")
                            
                            # Filter CellChatDB to selected signaling types
                            interaction = cellchatdb['interaction']
                            interaction_filtered = interaction[interaction['annotation'].isin(selected_types)]
                            cellchatdb['interaction'] = interaction_filtered
                            
                            st.write(f"Using {len(cellchatdb['interaction'])} ligand-receptor pairs")

                            #ここでdbを展開しておく
                            gene_use, resource, complex_input, cofactor_input, gene_info = extractGene(cellchatdb)

                            # データベースの検証
                            if resource.empty:
                                raise ValueError("DBの相互作用情報(interaction)が空です。有効なCellChatDBか確認してください。")
                                
                            # 必要な列の存在を確認
                            for required_col in ['ligand', 'receptor']:
                                if required_col not in resource.columns:
                                    raise ValueError(f"リソースデータフレームには '{required_col}' 列が必要です")

                        # 1. 細胞グループの細胞数チェック
                        cell_counts = adata1.obs[celltype_key].value_counts()
                        valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
                        if len(valid_groups) < len(cell_counts):
                            st.warning(f"{len(cell_counts) - len(valid_groups)}個の細胞グループが{min_cells}細胞未満のため除外されました")                   
                        # 有効なグループの細胞だけを保持
                        adata1_filtered = adata1[adata1.obs[celltype_key].isin(valid_groups)].copy()
                        valid_signaling_genes = [g for g in gene_use if g in adata1_filtered.var_names]
                        # シグナリング遺伝子のみを保持（R版のsubsetData相当の処理）
                        adata1_filtered = adata1_filtered[:, valid_signaling_genes].copy()
                        st.write(f"LR gene data of data1: {adata1_filtered.shape[0]} cells x {adata1_filtered.shape[1]} genes")

                        cell_counts = adata2.obs[celltype_key].value_counts()
                        valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
                        if len(valid_groups) < len(cell_counts):
                            st.warning(f"{len(cell_counts) - len(valid_groups)}個の細胞グループが{min_cells}細胞未満のため除外されました")                   
                        # 有効なグループの細胞だけを保持
                        adata2_filtered = adata2[adata2.obs[celltype_key].isin(valid_groups)].copy()
                        valid_signaling_genes = [g for g in gene_use if g in adata2_filtered.var_names]
                        # シグナリング遺伝子のみを保持（R版のsubsetData相当の処理）
                        adata2_filtered = adata2_filtered[:, valid_signaling_genes].copy()
                        st.write(f"LR gene data of data2: {adata2_filtered.shape[0]} cells x {adata2_filtered.shape[1]} genes")


                        # (1) adata1 で DE解析 ここからすべてfilterずみでーたにする
                        with st.spinner("Performing DE analysis for Group1..."):
                            de_result_g1 = identify_overexpressed_genes(
                                adata1_filtered,
                                group_by=celltype_key,  # 例: adata1の中のcelltype(あるいは同じsplitby_key)
                                do_de=True,
                                do_fast=True,
                                thresh_p=thresh_p,
                                thresh_pct=expr_prop,
                                min_cells=min_cells
                            )
                            feature_list_g1 = set(de_result_g1["features"])
                            st.write(f"Group1 DE: found {len(feature_list_g1)} genes")

                        # (2) adata2 で DE解析
                        with st.spinner("Performing DE analysis for Group2..."):
                            de_result_g2 = identify_overexpressed_genes(
                                adata2_filtered,
                                group_by=celltype_key,  # 例: adata2の中でも同じキーを利用
                                do_de=True,
                                do_fast=True,
                                thresh_p=thresh_p,
                                thresh_pct=expr_prop,
                                min_cells=min_cells
                            )
                            feature_list_g2 = set(de_result_g2["features"])
                            st.write(f"Group2 DE: found {len(feature_list_g2)} genes")

                        # (3) unionを取る

                        if union_mode=="union":
                            feature_union = list(feature_list_g1.union(feature_list_g2))
                        elif union_mode=="intersection":
                            feature_union = list(feature_list_g1 & feature_list_g2)
                        else:
                            feature_union = []
                        feature_list = list(feature_union)
                        st.write(f"Union of Group1/Group2 DE genes => {len(feature_list)} genes total")

                        


                        
                        # Run permutation analysis
                        with st.spinner(f"Running permutation analysis with {n_perms} permutations..."):
                            # Set up CellChat parameters
                            cellchat_params = {
                                'use_layer': data_layer,
                                'min_cells': min_cells,
                                'expr_prop': expr_prop,
                                'pseudocount': 1.0,
                                'trim_threshold': trim_threshold,
                                'k': 0.5,
                                'n': 1,
                                'type_mean': fun_type,
                                'raw_use': True,
                                'population_size': population_size,
                                'nboot': 100,  # Internal boostrapping in CellChat
                                'seed': 12345,
                                'key_added': "cellchat_res",
                                'trim': trim,
                                'apply_pval_filter': apply_pval_filter,
                                'features': feature_list,  # unionを設定
                            }
                            
                            # Get the list of cell types
                            cell_types = cell_list
                            
                            # Run permutation test
                            perm_results = run_permutation_cellchat(
                                adata1_filtered=adata1_filtered,
                                adata2_filtered=adata2_filtered,
                                groupby=celltype_key,
                                gene_use=gene_use,
                                complex_input=complex_input,
                                cofactor_input=cofactor_input,
                                resource=resource, 
                                cell_types=cell_types,
                                n_perm=n_perms,
                                n_jobs=n_cpus,
                                use_gpu=True if gpu_hybrid else False,
                                gpu_precision=gpu_precision,
                                hybrid_mode=False,
                                do_every_deg=do_every_deg,
                                union_mode=union_mode,
                                feature_list_g1=list(feature_list_g1) if union_mode =="individual" else None,
                                feature_list_g2=list(feature_list_g1) if union_mode =="individual" else None,
                                r_patcher=r_patcher,
                                #hybrid_mode=True if gpu_hybrid else False,
                                **cellchat_params
                            )
                            
                            # Store results in session state
                            st.session_state.perm_results = perm_results
                            st.session_state.cell_types = cell_types
                            st.session_state.group1 = group1
                            st.session_state.group2 = group2
                            
                            # Mark analysis as completed
                            st.session_state.analysis_completed = True
                            st.session_state.current_mode = "visualization"
                            
                            st.success(f"Permutation analysis completed successfully with {n_perms} permutations!")
                            
                            # Force a rerun to switch to visualization mode
                            st.rerun()
                    
                    except Exception as e:
                        st.error(f"Error during permutation analysis: {str(e)}")
                        st.exception(e)
                        st.session_state.current_mode = "setup"
                
                elif st.session_state.current_mode == "visualization" or st.session_state.analysis_completed:

                    
                    # Add a separate button for visualizing results
                    st.markdown("---")
                    if st.button("Show Visualizations", key="show_viz") or 'show_visualization' in st.session_state:

                        # Add visualization controls
                        st.header("Visualization Settings")
                        
                        # Initialize significance threshold in session state if not exists
                        if 'significance_threshold' not in st.session_state:
                            st.session_state.significance_threshold = 0.05
                        
                        # Function to update threshold in session state
                        def update_threshold():
                            st.session_state.significance_threshold = st.session_state.threshold_slider
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            significance_threshold = st.slider(
                                "Default significance threshold:",
                                min_value=0.001,
                                max_value=0.1,
                                value=st.session_state.significance_threshold,
                                step=0.001,
                                key="threshold_slider",
                                on_change=update_threshold
                            )


                        if 'perm_results' in st.session_state:
                            # Display results
                            perm_results = st.session_state.perm_results
                            cell_types = st.session_state.cell_types
                            group1 = st.session_state.group1
                            group2 = st.session_state.group2
                            
                            # Set flag to indicate visualization is shown
                            st.session_state.show_visualization = True
                            
                            # Display visualizations
                            display_visualizations(
                                perm_results=perm_results,
                                cell_types=cell_types,
                                group1=group1,
                                group2=group2,
                                significance_threshold=st.session_state.significance_threshold
                            )
                            
                            # Offer to save the entire permutation results
                            st.markdown("---")
                            st.header("Save Permutation Results")
                            if st.button("Prepare permutation results for download"):

                                try:
                                    # Prepare data for download
                                    result_data = {
                                        'group1': group1,
                                        'group2': group2,
                                        'n_perms': st.session_state.analysis_parameters['n_perms'],
                                        'actual_diff': perm_results['actual_diff'],
                                        'adjusted_p_values': perm_results['adjusted_p_values'],
                                        'parameters': st.session_state.analysis_parameters,  # cellchat_paramsの代わりにanalysis_parametersを使用
                                        'timestamp': time.time(),
                                        'file_name': uploaded_file.name
                                    }
                                        # pathway 情報を確保するために result1 と result2 の netP 部分を保存
                                    if 'result1' in perm_results and 'netP' in perm_results['result1']:
                                        result_data['result1'] = {'netP': perm_results['result1']['netP']}
                                    
                                    if 'result2' in perm_results and 'netP' in perm_results['result2']:
                                        result_data['result2'] = {'netP': perm_results['result2']['netP']}
                                    # Serialize to file
                                    with st.spinner("Serializing permutation results..."):
                                        pickled_data = pickle.dumps(result_data)
                                        
                                        # Offer download
                                        st.download_button(
                                            label="Download permutation analysis results",
                                            data=pickled_data,
                                            file_name=f"cellchat_permutation_{group1}_vs_{group2}_{n_perms}perms.pkl",
                                            mime="application/octet-stream"
                                        )
                                        
                                        st.success("Permutation results prepared for download!")
                                except Exception as e:
                                    st.error(f"Error preparing results for download: {str(e)}")
                                    st.exception(e)                       
                            
                        else:
                            st.error("No results to visualize. Please run the analysis first.")
                    
                    # Add a button to reset and run a new analysis
                    st.markdown("---")
                    if st.button('Reset Analysis'):
                        # Reset all relevant session state variables
                        st.session_state.analysis_completed = False
                        st.session_state.current_mode = "setup"
                        if 'perm_results' in st.session_state:
                            del st.session_state.perm_results
                        if 'show_visualization' in st.session_state:
                            del st.session_state.show_visualization
                        if 'selected_pathway' in st.session_state:
                            del st.session_state.selected_pathway
                        if 'significance_threshold' in st.session_state:
                            del st.session_state.significance_threshold
                        st.rerun()
            
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")
                st.exception(e)