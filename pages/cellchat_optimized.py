import numpy as np
import pandas as pd
import os
import pickle
import re
import time
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
import scipy.sparse as sparse
from scipy.stats import zscore
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
import psutil
import platform
from typing import Union, Optional, List, Dict, Any, Tuple

# ロギング設定
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("CellChat")

def configure_high_performance_system():
    """
    高性能計算環境（大容量メモリと多数のCPUコア）向けの最適なパラメータを取得する
    
    Returns
    -------
    dict
        最適化設定のディクショナリ
    """
    # システムリソースの検出
    mem = psutil.virtual_memory()
    total_mem_gb = mem.total / (1024 ** 3)
    available_mem_gb = mem.available / (1024 ** 3)
    cpu_count = os.cpu_count() or 4
    
    # システム情報のログ出力
    system_info = (
        f"システム情報: {platform.system()} {platform.release()}, "
        f"物理メモリ: {total_mem_gb:.1f}GB (利用可能: {available_mem_gb:.1f}GB), "
        f"CPUコア数: {cpu_count}"
    )
    print(system_info)
    
    # 高性能システム用の最適化設定
    if total_mem_gb > 200 and cpu_count >= 32:  # 256GB、64コアクラスのシステム
        print(f"高性能システムを検出: {total_mem_gb:.1f}GB RAM, {cpu_count}コア")
        
        # メモリが豊富な場合の並列度とバッチサイズを調整
        if total_mem_gb >= 256:
            # 超大容量メモリシステム (256GB+)
            recommended_workers = min(48, cpu_count)
            batch_size = 500
            permutation_chunk = 1000
            lr_chunk_size = 200
            cell_chunk_size = 50000  # 大きな細胞数データ用
            precompute_all = True    # すべての事前計算を有効化
        elif total_mem_gb >= 128:
            # 大容量メモリシステム (128-256GB)
            recommended_workers = min(32, cpu_count)
            batch_size = 300
            permutation_chunk = 500
            lr_chunk_size = 150
            cell_chunk_size = 20000
            precompute_all = True
        else:
            # 中程度の大容量メモリシステム (64-128GB)
            recommended_workers = min(24, cpu_count)
            batch_size = 200
            permutation_chunk = 200
            lr_chunk_size = 100 
            cell_chunk_size = 10000
            precompute_all = False
            
        # 有効なメモリがある場合のみ最大数のワーカーを使用
        memory_per_worker = 4.0  # ワーカーあたりの必要メモリ推定値(GB)
        max_workers_by_memory = int(available_mem_gb / memory_per_worker)
        recommended_workers = min(recommended_workers, max_workers_by_memory)
        
        return {
            "system_type": "high_performance",
            "available_memory_gb": available_mem_gb,
            "total_memory_gb": total_mem_gb,
            "cpu_count": cpu_count,
            "recommended_workers": recommended_workers,
            "recommended_batch_size": batch_size,
            "max_permutation_chunk": permutation_chunk,
            "lr_chunk_size": lr_chunk_size,
            "cell_chunk_size": cell_chunk_size,
            "precompute_all": precompute_all,
            "use_memory_mapping": True,  # 大きなデータセット用メモリマッピング
            "optimize_sparse": True     # スパース行列最適化
        }
    elif total_mem_gb > 64 and cpu_count >= 16:
        # 中程度の高性能システム
        print(f"中程度の高性能システムを検出: {total_mem_gb:.1f}GB RAM, {cpu_count}コア")
        recommended_workers = min(16, cpu_count)
        memory_per_worker = 3.0  # ワーカーあたりの必要メモリ(GB)
        max_workers_by_memory = int(available_mem_gb / memory_per_worker)
        recommended_workers = min(recommended_workers, max_workers_by_memory)
        
        return {
            "system_type": "medium_performance",
            "available_memory_gb": available_mem_gb,
            "total_memory_gb": total_mem_gb,
            "cpu_count": cpu_count,
            "recommended_workers": recommended_workers,
            "recommended_batch_size": 100,
            "max_permutation_chunk": 100,
            "lr_chunk_size": 50,
            "cell_chunk_size": 5000,
            "precompute_all": False,
            "use_memory_mapping": True,
            "optimize_sparse": True
        }
    else:
        # 標準システム
        print(f"標準システムを検出: {total_mem_gb:.1f}GB RAM, {cpu_count}コア")
        memory_per_worker = 2.0  # ワーカーあたりの必要メモリ(GB)
        recommended_workers = max(1, min(8, int(available_mem_gb / memory_per_worker)))
        
        return {
            "system_type": "standard",
            "available_memory_gb": available_mem_gb,
            "total_memory_gb": total_mem_gb,
            "cpu_count": cpu_count,
            "recommended_workers": recommended_workers,
            "recommended_batch_size": int(min(50, max(10, available_mem_gb * 10))),
            "max_permutation_chunk": 50,
            "lr_chunk_size": 20,
            "cell_chunk_size": 1000,
            "precompute_all": False,
            "use_memory_mapping": False,
            "optimize_sparse": True
        }

def sanitize_filename(filename, max_length=20):
    """ファイル名を安全にして最大長さを制限する"""
    # 特殊文字を削除または置換
    filename = re.sub(r'[\\/*?:"<>|]', "_", filename)
    # 長さを制限
    if len(filename) > max_length:
        base, ext = os.path.splitext(filename)
        filename = base[:max_length] + ext
    return filename

def compute_hill_outer_vectorized(dataL, dataR, k, n):
    """
    ベクトル化されたHill関数の外積計算
    
    Parameters
    ----------
    dataL : numpy.ndarray
        リガンドの発現レベル
    dataR : numpy.ndarray
        レセプターの発現レベル
    k : float
        Hill係数パラメータk
    n : float
        Hill係数パラメータn
        
    Returns
    -------
    numpy.ndarray
        相互作用確率の行列
    """
    # ブロードキャスティング用に形状を変更
    L_reshaped = dataL.reshape(-1, 1)  # 列ベクトル
    R_reshaped = dataR.reshape(1, -1)  # 行ベクトル
    
    # ベクトル化された外積
    x = L_reshaped * R_reshaped
    
    # ベクトル化されたHill関数
    kn = k ** n
    xn = x ** n
    result = xn / (kn + xn)
    
    return result

def calculate_expression_efficiently(data_use, cell_indices, fun_type='triMean', trim=0.1):
    """
    効率的に細胞のサブセットの発現値を計算する
    
    Parameters
    ----------
    data_use : numpy.ndarray
        発現データ (cells × genes)
    cell_indices : numpy.ndarray
        含める細胞のインデックス
    fun_type : str
        平均計算の種類: 'triMean', 'truncatedMean', または 'mean'
    trim : float
        truncatedMeanのトリムパラメータ
        
    Returns
    -------
    numpy.ndarray
        遺伝子あたりの平均発現値
    """
    if len(cell_indices) == 0:
        return np.zeros(data_use.shape[1])
        
    data_subset = data_use[cell_indices]
    
    if fun_type == 'triMean':
        # 効率的な分位数計算 - 転置して性能向上
        data_t = data_subset.T
        q1 = np.nanquantile(data_t, 0.25, axis=1)
        q2 = np.nanquantile(data_t, 0.50, axis=1)
        q3 = np.nanquantile(data_t, 0.75, axis=1)
        return (q1 + 2*q2 + q3) / 4
    elif fun_type == 'truncatedMean':
        # nanmeanのパフォーマンス向上のために転置
        data_t = data_subset.T
        n_elements = data_t.shape[1]
        k = int(trim * n_elements)
        if k > 0:
            # 各行（遺伝子）をソート
            sorted_data = np.sort(data_t, axis=1)
            # 両端をトリム
            trimmed_data = sorted_data[:, k:-k] if 2*k < n_elements else sorted_data
            return np.nanmean(trimmed_data, axis=1)
        else:
            return np.nanmean(data_t, axis=1)
    else:
        # デフォルトで平均に
        return np.nanmean(data_subset, axis=0)

def process_permutation_batch_optimized(batch_indices, data_use, cell_labels, permutation, valid_cell_types, 
                                       fun_type='triMean', trim=0.1):
    """
    パーミュテーションのバッチを効率的に処理する
    
    Parameters
    ----------
    batch_indices : list
        処理するパーミュテーションのインデックス
    data_use : numpy.ndarray
        正規化された発現行列 (cells × genes)
    cell_labels : pandas.Series
        細胞タイプラベル
    permutation : numpy.ndarray
        パーミュテーションインデックス (cells × nboot)
    valid_cell_types : list
        有効な細胞タイプのリスト
    fun_type : str
        平均計算の種類
    trim : float
        トリムパラメータ
        
    Returns
    -------
    numpy.ndarray
        バッチ結果 (genes × clusters × batch_size)
    """
    n_genes = data_use.shape[1]
    numCluster = len(valid_cell_types)
    batch_size = len(batch_indices)
    results = np.zeros((n_genes, numCluster, batch_size), dtype=np.float32)
    
    # 各細胞タイプのインデックスを事前計算してルックアップパフォーマンスを向上
    cell_type_indices = {ct: np.where(cell_labels == ct)[0] for ct in valid_cell_types}
    
    for idx, j in enumerate(batch_indices):
        # パーミュテーションされた細胞ラベルを取得
        perm_indices = permutation[:, j]
        
        # 各細胞タイプを処理
        for ct_idx, ct in enumerate(valid_cell_types):
            # この細胞タイプのオリジナルインデックスを取得
            orig_indices = cell_type_indices[ct]
            
            # パーミュテーションされたインデックスを取得
            perm_ct_indices = np.sort(perm_indices[orig_indices])
            
            # この細胞タイプの平均発現を計算
            results[:, ct_idx, idx] = calculate_expression_efficiently(
                data_use, perm_ct_indices, fun_type, trim
            )
    
    return results

def precompute_gene_expressions_optimized(data_use, cell_labels, permutation, valid_cell_types, fun_type='triMean', nboot=100, n_jobs=8, trim=0.1, max_chunk_size=None):
    """
    遺伝子発現の最適化された事前計算バージョン
    
    Parameters
    ----------
    data_use : numpy.ndarray
        正規化された発現行列 (cells × genes)
    cell_labels : pandas.Series
        細胞タイプラベル
    permutation : numpy.ndarray
        パーミュテーションインデックス (cells × nboot)
    valid_cell_types : list
        有効な細胞タイプのリスト
    fun_type : str
        平均計算の種類
    nboot : int
        パーミュテーション数
    n_jobs : int
        並列ジョブ数
    trim : float
        トリムパラメータ
    max_chunk_size : int, optional
        最大チャンクサイズ。システム設定から決定される場合はNone
        
    Returns
    -------
    numpy.ndarray
        事前計算された遺伝子発現 (genes × clusters × nboot)
    """
    logger.info("すべてのパーミュテーション用に平均遺伝子発現を事前計算...")
    
    numCluster = len(valid_cell_types)
    n_genes = data_use.shape[1]
    
    # システムによって適切なバッチサイズを決定
    total_elements = n_genes * numCluster * nboot
    element_size = data_use.dtype.itemsize  # 各要素のサイズ（バイト）
    total_memory = total_elements * element_size
    
    if max_chunk_size is None:
        # メモリに基づいてバッチサイズを調整
        system_config = configure_high_performance_system()
        if system_config["system_type"] == "high_performance":
            batch_size = system_config["max_permutation_chunk"]
        else:
            # 推定メモリに基づいて調整
            if total_memory > 2e9:  # 2GB以上
                batch_size = min(20, max(5, int(2e9 / (n_genes * numCluster * element_size))))
            else:
                batch_size = min(50, nboot)
    else:
        batch_size = max_chunk_size
        
    logger.info(f"パーミュテーションにバッチサイズを使用: {batch_size}")
    
    # 結果配列を初期化
    all_gene_expr = np.zeros((n_genes, numCluster, nboot), dtype=np.float32)
    
    # バッチに分割
    all_batches = [list(range(i, min(i+batch_size, nboot))) for i in range(0, nboot, batch_size)]
    
    # CPU数とバッチ数に基づいてジョブ数を制限
    n_jobs_effective = min(n_jobs, len(all_batches), os.cpu_count() or 1)
    
    if n_jobs_effective > 1:
        # バッチを並列処理
        logger.info(f"パーミュテーションを{n_jobs_effective}ジョブで並列処理...")
        
        # 改良された進捗レポート
        progress_bar = tqdm(total=len(all_batches), desc="パーミュテーションバッチを処理中")
        
        results = Parallel(n_jobs=n_jobs_effective, backend="loky")(
            delayed(process_permutation_batch_optimized)(
                batch_indices, data_use, cell_labels, permutation, valid_cell_types, fun_type, trim
            ) for batch_indices in all_batches
        )
        
        # 結果を結合
        for batch_idx, batch_indices in enumerate(all_batches):
            for idx, j in enumerate(batch_indices):
                if j < nboot:  # 範囲チェック
                    all_gene_expr[:, :, j] = results[batch_idx][:, :, idx]
            progress_bar.update(1)
            
        progress_bar.close()
    else:
        # シングルコア処理
        logger.info("パーミュテーションを順次処理...")
        progress_bar = tqdm(total=nboot, desc="パーミュテーションを処理中")
        
        for j in range(nboot):
            # パーミュテーションされた細胞ラベルを取得
            group_boot = cell_labels.values[permutation[:, j]]
            
            # 各細胞タイプを処理
            for ct_idx, ct in enumerate(valid_cell_types):
                # この細胞タイプのインデックスを取得
                indices = np.where(group_boot == ct)[0]
                
                # 平均発現を計算
                all_gene_expr[:, ct_idx, j] = calculate_expression_efficiently(
                    data_use, indices, fun_type, trim
                )
                
            progress_bar.update(1)
            
        progress_bar.close()
    
    logger.info("平均遺伝子発現計算完了.")
    return all_gene_expr

def optimize_matrix_operations(X, adata, use_memory_mapping=False, cell_chunk_size=10000):
    """
    行列演算を最適化してパフォーマンスを向上させる
    
    Parameters
    ----------
    X : scipy.sparse.spmatrix or numpy.ndarray
        最適化する発現行列
    adata : AnnData
        データを含むAnnDataオブジェクト
    use_memory_mapping : bool
        メモリマッピングを使用するかどうか
    cell_chunk_size : int
        大きな行列を処理する際のチャンクサイズ
        
    Returns
    -------
    numpy.ndarray
        最適化された密行列
    """
    logger.info(f"行列形状: {X.shape}")
    
    # スパース行列の最適化
    if sparse.issparse(X):
        logger.info("スパース行列を密行列に変換...")
        
        # メモリマッピングオプション
        if use_memory_mapping and X.shape[0] * X.shape[1] > 1e8:
            import tempfile
            import os
            
            # 一時ファイルを作成
            with tempfile.NamedTemporaryFile(suffix='.dat', delete=False) as f:
                temp_file = f.name
            
            try:
                # メモリマップ配列を作成
                result = np.memmap(temp_file, dtype=np.float32, mode='w+', 
                                  shape=(X.shape[0], X.shape[1]))
                
                # チャンクで処理
                for i in range(0, X.shape[0], cell_chunk_size):
                    end = min(i + cell_chunk_size, X.shape[0])
                    result[i:end] = X[i:end].toarray()
                
                # メモリマップをディスクと同期
                result.flush()
                
                # 読み取り専用モードでマップを再オープン
                result = np.memmap(temp_file, dtype=np.float32, mode='r', 
                                  shape=(X.shape[0], X.shape[1]))
                
                return result
            except Exception as e:
                logger.error(f"メモリマッピングエラー: {str(e)}")
                logger.info("標準の処理方法にフォールバック")
                os.unlink(temp_file)  # 一時ファイルを削除
        
        # 大きな行列の場合、メモリ使用量を減らすためにチャンクで処理
        if X.shape[0] > cell_chunk_size or X.shape[1] > cell_chunk_size:
            # メモリ制限に基づいて最適なチャンクサイズを決定
            mem_limit = 1e9  # チャンクあたり1GB制限
            element_size = 4  # float32の場合4バイト
            
            # 行または列のどちらが大きいかによってチャンキング方式を決定
            if X.shape[0] > X.shape[1]:
                # 行でチャンク
                chunk_size = min(cell_chunk_size, int(mem_limit / (X.shape[1] * element_size)))
                logger.info(f"大きな行列を{chunk_size}行ずつチャンクで処理")
                
                result = np.zeros((X.shape[0], X.shape[1]), dtype=np.float32)
                for i in range(0, X.shape[0], chunk_size):
                    end = min(i + chunk_size, X.shape[0])
                    result[i:end] = X[i:end].toarray()
                    
                return result
            else:
                # 列でチャンク
                chunk_size = min(cell_chunk_size, int(mem_limit / (X.shape[0] * element_size)))
                logger.info(f"大きな行列を{chunk_size}列ずつチャンクで処理")
                
                result = np.zeros((X.shape[0], X.shape[1]), dtype=np.float32)
                for i in range(0, X.shape[1], chunk_size):
                    end = min(i + chunk_size, X.shape[1])
                    result[:, i:end] = X[:, i:end].toarray()
                    
                return result
        else:
            # 小さな行列の場合、直接変換
            return X.toarray().astype(np.float32)
    else:
        # すでに密行列の場合
        if X.dtype != np.float32:
            # 精度と節約のためにfloat32に変換
            return X.astype(np.float32)
        return X

def precompute_complex_mappings(complex_input, cofactor_input, resource, gene_to_index):
    """
    すべての複合体とコファクターのマッピングを事前計算して繰り返し検索を避ける
    
    Parameters
    ----------
    complex_input : pd.DataFrame
        CellChatDBからの複合体情報
    cofactor_input : pd.DataFrame
        CellChatDBからのコファクター情報
    resource : pd.DataFrame
        CellChatDBからの相互作用情報
    gene_to_index : dict
        遺伝子名からインデックスへのマッピング
        
    Returns
    -------
    dict
        複合体、コファクター、アゴニスト、アンタゴニストのマッピング情報
    """
    # 複合体マッピング
    complex_mapping = {}
    if not complex_input.empty:
        for complex_name in complex_input.index:
            subunits_cols = [col for col in complex_input.columns if 'subunit' in col]
            subunits = complex_input.loc[complex_name, subunits_cols].dropna().astype(str)
            subunits = [s for s in subunits if s != "" and s in gene_to_index]
            
            if subunits:
                complex_mapping[complex_name] = [gene_to_index[s] for s in subunits]
    
    # コファクターマッピング
    cofactor_mapping = {}
    if not cofactor_input.empty:
        for cofactor_name in cofactor_input.index:
            cofactor_cols = [col for col in cofactor_input.columns if 'cofactor' in col]
            cofactors = cofactor_input.loc[cofactor_name, cofactor_cols].dropna().astype(str)
            cofactors = [c for c in cofactors if c != "" and c in gene_to_index]
            
            if cofactors:
                cofactor_mapping[cofactor_name] = [gene_to_index[c] for c in cofactors]
    
    # アゴニストとアンタゴニストのマッピング
    agonist_mapping = {}
    antagonist_mapping = {}
    
    if 'agonist' in resource.columns:
        for i, agonist in enumerate(resource['agonist']):
            if pd.notna(agonist) and agonist != "" and agonist in cofactor_mapping:
                agonist_mapping[i] = cofactor_mapping[agonist]
    
    if 'antagonist' in resource.columns:
        for i, antagonist in enumerate(resource['antagonist']):
            if pd.notna(antagonist) and antagonist != "" and antagonist in cofactor_mapping:
                antagonist_mapping[i] = cofactor_mapping[antagonist]
    
    return {
        'complex': complex_mapping,
        'cofactor': cofactor_mapping,
        'agonist': agonist_mapping,
        'antagonist': antagonist_mapping
    }

def compute_pairwise_interaction_probabilities(dataLavg, dataRavg, k, n, index_agonist, index_antagonist, 
                                      data_use_avg_df, resource, cofactor_input, dataLavg2, dataRavg2,
                                      population_size, numCluster, nLR, lr_chunk_size=100):
    """
    すべてのペアワイズ相互作用確率をベクトル化された方法で計算
    
    Parameters
    ----------
    dataLavg : numpy.ndarray
        リガンド発現行列 (LR_pairs × clusters)
    dataRavg : numpy.ndarray
        レセプター発現行列 (LR_pairs × clusters)
    k : float
        Hill係数パラメータk
    n : float
        Hill係数パラメータn
    index_agonist : list
        アゴニスト相互作用のインデックス
    index_antagonist : list
        アンタゴニスト相互作用のインデックス
    data_use_avg_df : pd.DataFrame
        細胞タイプごとの平均発現データ
    resource : pd.DataFrame
        CellChatDBからの相互作用情報
    cofactor_input : pd.DataFrame
        CellChatDBからのコファクター情報
    dataLavg2 : numpy.ndarray
        人口調整済みリガンド発現 (LR_pairs × clusters)
    dataRavg2 : numpy.ndarray
        人口調整済みレセプター発現 (LR_pairs × clusters)
    population_size : bool
        人口サイズ効果を考慮するかどうか
    numCluster : int
        細胞クラスター数
    nLR : int
        LRペア数
    lr_chunk_size : int
        LRペア処理のバッチサイズ
        
    Returns
    -------
    numpy.ndarray
        相互作用確率行列 (clusters × clusters × LR_pairs)
    """
    # 結果配列を初期化
    Prob = np.zeros((numCluster, numCluster, nLR))
    
    # LRペアをバッチで処理
    batch_size = min(lr_chunk_size, nLR)
    
    for batch_start in range(0, nLR, batch_size):
        batch_end = min(batch_start + batch_size, nLR)
        batch_size_actual = batch_end - batch_start
        
        # バッチ確率を初期化
        P1_batch = np.zeros((batch_size_actual, numCluster, numCluster))
        P2_batch = np.ones((batch_size_actual, numCluster, numCluster))
        P3_batch = np.ones((batch_size_actual, numCluster, numCluster))
        P4_batch = np.ones((batch_size_actual, numCluster, numCluster))
        
        # バッチ内のすべてのLRペアに対してP1（Hill関数）を計算
        for batch_idx, i in enumerate(range(batch_start, batch_end)):
            # ベクトル化されたHill関数を使用
            P1_batch[batch_idx] = compute_hill_outer_vectorized(dataLavg[i], dataRavg[i], k, n)
        
        
        # 関連するLRペアのP3（アンタゴニスト効果）を計算
        for batch_idx, i in enumerate(range(batch_start, batch_end)):
            if i in index_antagonist:
                data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P3_batch[batch_idx] = np.outer(data_antagonist, data_antagonist)
        
        # 必要な場合はP4（人口サイズ効果）を計算
        if population_size:
            for batch_idx, i in enumerate(range(batch_start, batch_end)):
                P4_batch[batch_idx] = np.outer(dataLavg2[i], dataRavg2[i])
        
        # バッチの最終確率を計算
        for batch_idx, i in enumerate(range(batch_start, batch_end)):
            Prob[:, :, i] = P1_batch[batch_idx] * P2_batch[batch_idx] * P3_batch[batch_idx] * P4_batch[batch_idx]

def calculate_lr_interaction(
    i, dataLavg, dataRavg, dataLavg2, dataRavg2,
    index_agonist, index_antagonist, data_use_avg_df,
    resource, cofactor_input, k, n, population_size,
    permutation, cell_labels, valid_cell_types, all_gene_expr,
    ligand_indices, receptor_indices, nboot, numCluster
):
    """
    単一のLRペアの相互作用確率とp値を計算
    
    Parameters
    ----------
    i : int
        LRペアのインデックス
    その他のパラメータは計算に必要なすべてのデータ
    
    Returns
    -------
    tuple
        このLRペアの(Prob, Pval)
    """
    # 相互作用確率を計算
    dataLR = np.outer(dataLavg[i, :], dataRavg[i, :])
    
    # Hill関数
    P1 = compute_hill_outer_vectorized(dataLavg[i, :], dataRavg[i, :], k, n)
    
    # アゴニスト効果
    P2 = np.ones((numCluster, numCluster))
    if i in index_agonist:
        data_agonist = computeExpr_agonist(data_use_avg_df, resource, cofactor_input, i, k, n)
        P2 = np.outer(data_agonist, data_agonist)
    
    # アンタゴニスト効果
    P3 = np.ones((numCluster, numCluster))
    if i in index_antagonist:
        data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
        P3 = np.outer(data_antagonist, data_antagonist)
    
    # 人口サイズ効果
    P4 = np.ones((numCluster, numCluster))
    if population_size:
        P4 = np.outer(dataLavg2[i, :], dataRavg2[i, :])
    
    # 最終確率
    Pnull = P1 * P2 * P3 * P4
    
    # 相互作用がない場合はp値計算をスキップ
    if np.sum(Pnull) == 0:
        return Pnull, np.ones((numCluster, numCluster))
    
    Pnull_vec = Pnull.flatten()
    
    # リガンドとレセプターの情報を取得
    ligand_info = ligand_indices[i]
    receptor_info = receptor_indices[i]
    
    # 発現が取得できない場合はスキップ
    if ligand_info[2] is None or receptor_info[2] is None:
        return Pnull, np.ones((numCluster, numCluster))
    
    # ブートストラップ確率を計算
    Pboot = np.zeros((numCluster * numCluster, nboot))
    
    for j in range(nboot):
        # リガンド発現を取得
        lr_i, ligand_gene_indices, is_ligand_complex = ligand_info
        
        if not is_ligand_complex:
            # 単一遺伝子
            if ligand_gene_indices:
                ligand_idx = ligand_gene_indices[0]
                dataLavgB = all_gene_expr[ligand_idx, :, j].reshape(1, -1)
            else:
                dataLavgB = np.zeros((1, numCluster))
        else:
            # 複合体 - 幾何平均を計算
            expr_values = np.array([all_gene_expr[idx, :, j] for idx in ligand_gene_indices])
            if len(expr_values) > 0:
                # log変換して平均を取り、元に戻す
                log_values = np.log(expr_values + 1e-10)
                dataLavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
            else:
                dataLavgB = np.zeros((1, numCluster))
        
        # レセプター発現を取得
        lr_i, receptor_gene_indices, is_receptor_complex = receptor_info
        
        if not is_receptor_complex:
            if receptor_gene_indices:
                receptor_idx = receptor_gene_indices[0]
                dataRavgB = all_gene_expr[receptor_idx, :, j].reshape(1, -1)
            else:
                dataRavgB = np.zeros((1, numCluster))
        else:
            # 複合体 - 幾何平均を計算
            expr_values = np.array([all_gene_expr[idx, :, j] for idx in receptor_gene_indices])
            if len(expr_values) > 0:
                log_values = np.log(expr_values + 1e-10)
                dataRavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
            else:
                dataRavgB = np.zeros((1, numCluster))
        
        # 相互作用確率を計算
        dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
        P1_boot = dataLRB**n / (k**n + dataLRB**n)
        
        # 簡略化されたアゴニスト/アンタゴニスト効果
        P2_boot = np.ones((numCluster, numCluster))
        P3_boot = np.ones((numCluster, numCluster))
        
        # 人口サイズ効果
        P4_boot = np.ones((numCluster, numCluster))
        if population_size:
            group_boot = cell_labels.values[permutation[:, j]]
            cell_proportions_boot = np.array([np.sum(group_boot == ct) for ct in valid_cell_types]) / len(cell_labels)
            P4_boot = np.outer(cell_proportions_boot, cell_proportions_boot)
        
        # 最終確率
        Pboot_result = P1_boot * P2_boot * P3_boot * P4_boot
        Pboot[:, j] = Pboot_result.flatten()
    
    # p値を計算
    nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
    p = nReject / nboot
    
    return Pnull, p.reshape(numCluster, numCluster)uster)
                else:
                    # 複合体 - 幾何平均を計算
                    expr_values = np.array([all_gene_expr[idx, :, j] for idx in lig_indices])
                    if len(expr_values) > 0:
                        log_values = np.log(expr_values + 1e-10)
                        ligand_expr = np.exp(np.mean(log_values, axis=0))
                    else:
                        ligand_expr = np.zeros(numCluster)
                
                # レセプター発現を取得
                lr_i, rec_indices, is_receptor_complex = receptor_info
                
                if not is_receptor_complex:
                    if rec_indices:
                        receptor_idx = rec_indices[0]
                        receptor_expr = all_gene_expr[receptor_idx, :, j]
                    else:
                        receptor_expr = np.zeros(numCluster)
                else:
                    expr_values = np.array([all_gene_expr[idx, :, j] for idx in rec_indices])
                    if len(expr_values) > 0:
                        log_values = np.log(expr_values + 1e-10)
                        receptor_expr = np.exp(np.mean(log_values, axis=0))
                    else:
                        receptor_expr = np.zeros(numCluster)
                
                # 注目している細胞ペアに対してベクトル化計算
                for pair_idx in range(len(row_idx)):
                    r, c = row_idx[pair_idx], col_idx[pair_idx]
                    
                    # P1（Hill関数）を計算
                    dataLR = ligand_expr[r] * receptor_expr[c]
                    P1_boot = dataLR**n / (k**n + dataLR**n)
                    
                    # 人口サイズ効果（有効な場合）
                    P4_boot = 1.0
                    if population_size:
                        group_boot = cell_labels.values[permutation[:, j]]
                        cell_proportions_boot = np.array([np.sum(group_boot == ct) for ct in valid_cell_types]) / nC
                        P4_boot = cell_proportions_boot[r] * cell_proportions_boot[c]
                    
                    # このペアの結果を保存
                    Pboot[pair_idx, j] = P1_boot * P4_boot
            
            # p値を計算
            p_values = np.zeros(len(row_idx))
            for pair_idx in range(len(row_idx)):
                nReject = np.sum(Pboot[pair_idx, :] > Pnull_values[pair_idx])
                p_values[pair_idx] = nReject / nboot
            
            batch_results.append((i, row_idx, col_idx, p_values))
        
        return batch_results
    
    # バッチを処理
    if n_jobs_effective > 1:
        logger.info(f"p値計算を{n_jobs_effective}個のジョブで並列処理...")
        all_results = Parallel(n_jobs=n_jobs_effective)(
            delayed(process_batch)(batch) for batch in all_batches
        )
        
        # 結果を平坦化
        flat_results = [item for sublist in all_results for item in sublist]
    else:
        logger.info("p値計算を順次処理...")
        flat_results = process_batch(nonzero_indices)
    
    # Pval配列を更新
    for i, row_idx, col_idx, p_values in flat_results:
        for pair_idx in range(len(row_idx)):
            r, c = row_idx[pair_idx], col_idx[pair_idx]
            Pval[r, c, i] = p_values[pair_idx]
    
    return Pval