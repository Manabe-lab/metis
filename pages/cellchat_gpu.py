import scanpy as sc
import numpy as np
import streamlit as st
import pandas as pd
import pickle
import re
import os
import time
import matplotlib
matplotlib.use("cairo") # PDFの透明化対策 interactiveではうまく行かない?
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from matplotlib.colors import LinearSegmentedColormap
#from scipy import sparse
import scipy.sparse
#from scipy.stats import zscore
from statsmodels.stats.multitest import multipletests
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
import logging
import matplotlib.cm as cm
import shutil
import io
import zipfile
from datetime import datetime, timedelta
import gc

# GPU imports with fallback
try:
    import cupy as cp
    import cupyx.scipy.sparse as cusparse
    HAS_GPU = True
except ImportError:
    HAS_GPU = False
    cp = None
    cusparse = None

logger = logging.getLogger("streamlit.runtime.scriptrunner_utils.script_run_context")
logger.disabled = True

from pages.cellchat_vis import (
    netVisual_circle, netVisual_circle_individual, netVisual_chord_by_pathway,
    netVisual_aggregate, netVisual_chord, netVisual_heatmap,
    netAnalysis_signalingRole_network, plotGeneExpression)

def init_gpu(use_gpu: bool = True, precision: str = 'float32', gpu_memory_limit: float = 0.8):
    """
    GPU初期化とメモリ情報取得
    
    Parameters
    ----------
    use_gpu : bool
        GPU使用フラグ
    precision : str
        計算精度 ('float32' or 'float64')
    gpu_memory_limit : float
        GPU メモリ使用率の上限（0.0-1.0）
        
    Returns
    -------
    tuple
        (use_gpu, has_gpu, cp_module, gpu_dtype, memory_info)
    """
    if not use_gpu or not HAS_GPU:
        cpu_dtype = getattr(np, precision)
        return False, False, None, cpu_dtype, {}
    
    try:
        # GPU メモリ情報取得
        device = cp.cuda.Device()
        total_memory = device.mem_info[1]
        available_memory = device.mem_info[0]
        
        memory_info = {
            'total_gb': total_memory / 1e9,
            'available_gb': available_memory / 1e9,
            'usage_limit_gb': available_memory * gpu_memory_limit / 1e9
        }
        
        gpu_dtype = getattr(cp, precision)
        
        # テスト用の小さな配列で動作確認
        test_array = cp.array([1, 2, 3], dtype=gpu_dtype)
        del test_array
        
        print(f"GPU initialized: {memory_info['available_gb']:.1f}GB available")
        return True, True, cp, gpu_dtype, memory_info
        
    except Exception as e:
        print(f"GPU initialization failed: {e}, falling back to CPU")
        cpu_dtype = getattr(np, precision)
        return False, False, None, cpu_dtype, {}

def estimate_memory_usage(n_cells: int, n_genes: int, n_permutations: int, 
                         precision: str = 'float32') -> dict:
    """GPU/CPU メモリ使用量の推定"""
    element_size = 4 if precision == 'float32' else 8
    
    # 基本データ行列
    data_memory = n_cells * n_genes * element_size
    
    # Permutation行列
    perm_memory = n_cells * n_permutations * 4  # int32
    
    # 中間結果（cell type別平均など）
    n_types = min(50, max(5, n_cells // 100))  # 推定cell type数
    intermediate_memory = n_types * n_genes * element_size * 2
    
    # 結果行列
    result_memory = n_types * n_types * n_permutations * element_size
    
    total_memory = data_memory + perm_memory + intermediate_memory + result_memory
    
    return {
        'data': data_memory,
        'permutation': perm_memory,
        'intermediate': intermediate_memory,
        'result': result_memory,
        'total': total_memory,
        'total_gb': total_memory / 1e9
    }

def determine_optimal_batch_size(n_cells: int, n_genes: int, n_permutations: int,
                                use_gpu: bool, available_memory_gb: float = 8.0) -> int:
    """最適なバッチサイズの決定"""
    if not use_gpu:
        return min(50, n_permutations)  # CPU では小さめのバッチ
    
    # GPU メモリ制約を考慮
    memory_estimate = estimate_memory_usage(n_cells, n_genes, 100)
    available_memory = available_memory_gb * 1e9 * 0.8  # 80%を使用
    
    # 利用可能メモリに基づくバッチサイズ
    max_batch_size = int(available_memory / (memory_estimate['total'] / 100))
    optimal_batch_size = max(10, min(max_batch_size, n_permutations, 200))
    
    return optimal_batch_size

def get_r_permutation(nC, nboot, seed=1):
    """
    R言語のsample.int関数を使用して置換データを生成し、
    エラー処理を含めた堅牢な実装
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        
        # Rの基本パッケージをインポート
        base = importr('base')
        
        # Rで乱数シードを設定
        ro.r('set.seed({})'.format(seed))
        
        # R側でpermutationデータを生成（トラブルシューティング用のコメント付き）
        r_code = """
        nC <- {}
        nboot <- {}
        cat("Generating permutation in R with nC =", nC, "and nboot =", nboot, "\\n")
        permutation <- replicate(nboot, sample.int(nC, size = nC))
        cat("R permutation dimensions:", dim(permutation), "\\n")
        permutation
        """.format(nC, nboot)
        
        # Rコードを実行して結果を取得
        r_permutation = ro.r(r_code)
        
        # R行列をNumPy配列に変換
        permutation = np.array(r_permutation).astype(int)
        
        # Rは1ベースなので、Pythonの0ベースに変換
        permutation = permutation - 1
        
        # 配列の形状を確認し、必要に応じて転置
        if permutation.shape[0] == nboot and permutation.shape[1] == nC:
            permutation = permutation.T  # CellChatの期待する形式に合わせる
        
        print(f"Successfully generated R permutation with shape: {permutation.shape}")
        print(f"Permutation range: min={np.min(permutation)}, max={np.max(permutation)}")
        print(f"Expected range: 0 to {nC-1}")
        return permutation
        
    except Exception as e:
        print(f"Error generating R permutation: {str(e)}")
        print("Falling back to numpy permutation")
        
        # フォールバック: NumPyの乱数生成を使用
        np.random.seed(seed)
        permutation = np.zeros((nC, nboot), dtype=int)
        for i in range(nboot):
            permutation[:, i] = np.random.permutation(nC)
        return permutation

@st.cache_data
def create_cell_color_mapping(cell_list, palette_name):
    """
    細胞名と色の一貫したマッピングを作成する関数

    Parameters
    ----------
    cell_list : list
        細胞名のリスト
    palette_name : str
        使用する離散カラーパレット名

    Returns
    -------
    dict
        細胞名をキー、色をバリューとする辞書

    Note
    ----
    指定された離散パレットの場合、まずデフォルトの色数（base_palette の長さ）を取得し、
    細胞数がその数以下ならそのまま利用、超える場合は sns.color_palette() の n_colors に
    細胞数を指定して補間色を生成します。
    """
    n_cells = len(cell_list)
    # デフォルトの離散パレットを取得（例: "Set1"なら9色）
    base_palette = sns.color_palette(palette_name)
    base_n = len(base_palette)
    
    if n_cells <= base_n:
        # 基本パレットの色数以内なら、先頭から必要な数を利用
        colors = base_palette[:n_cells]
    else:
        # 細胞数が基本パレットの色数を超える場合、必要な数の色を生成（線形補間）
        colors = sns.color_palette(palette_name, n_colors=n_cells)
        
    # 細胞名と色を対応付けた辞書を返す
    return {cell: color for cell, color in zip(cell_list, colors)}

@st.cache_data
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
            
            # データの次元チェック
            print(f"data_use_df shape: {data_use_df.shape}, ラベル数: {len(labels)}")
            
            # 次元が一致しない場合、ファイル経由でフォールバック
            if data_use_df.shape[1] != len(labels):
                print("警告: データの次元不一致。ファイル経由でRに処理を依頼します")
                genes_de = run_presto_through_files(data_use_df, labels, level_use)
            else:
                try:
                    # 通常のrpy2経由の処理を試みる
                    # labels を R の factor に変換
                    r_labels = ro.FactorVector(labels.astype(str).tolist())
                    # groups_use パラメータには、対象とするグループのベクトルを渡す
                    r_groups = ro.StrVector(level_use)
                    
                    # presto の wilcoxauc 呼び出し
                    res = presto.wilcoxauc(data_use_df, r_labels, groups_use=r_groups)
                    
                    # res は R の DataFrame になるので、pandas に変換
                    with localconverter(ro.default_converter + pandas2ri.converter):
                        genes_de = pd.DataFrame(ro.conversion.rpy2py(res))
                except Exception as e:
                    print(f"rpy2経由の実行エラー: {str(e)}")
                    print("ファイル経由でRに処理を依頼します")
                    genes_de = run_presto_through_files(data_use_df, labels, level_use)
            
            # R版と同様に、列名の変更
            genes_de.rename(columns={
                "group": "clusters",
                "feature": "features",
                "pval": "pvalues",
                "logFC": "logFC",
                "pct_in": "pct.1",
                "pct_out": "pct.2",
                "padj": "pvalues.adj"
            }, inplace=True)
            genes_de.loc[:, "logFC_abs"]  = genes_de["logFC"].abs().copy()
            # pct.max は pct.1 と pct.2 の最大値
            genes_de.loc[:, "pct.max"] = genes_de[["pct.1", "pct.2"]].max(axis=1)
            
            # フィルタリング
            markers_all = genes_de[(genes_de["pvalues"] < thresh_p) &
                                    (genes_de["logFC_abs"] >= thresh_fc) &
                                    (genes_de["pct.max"] > thresh_pct * 100)].copy()
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
                                      (markers_all["pct.max"] > thresh_pct)].copy()
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
    print(f"{len(features_sig)} genes passed filtering")
    return {"features": features_sig, "markers_all": markers_all}


def run_presto_through_files(data_df, labels, groups_use):
    """
    ファイル経由でR prestoを実行する
    
    Parameters:
    -----------
    data_df : pandas.DataFrame
        遺伝子×細胞の行列
    labels : pandas.Series
        細胞のグループラベル
    groups_use : list
        使用するグループリスト
        
    Returns:
    --------
    pandas.DataFrame
        wilcoxauc結果
    """
    import os
    import tempfile
    import subprocess
    import uuid
    
    # 一時ディレクトリを作成
    temp_dir = tempfile.mkdtemp(prefix="presto_")
    
    try:
        # ユニークなセッションID
        session_id = str(uuid.uuid4())
        
        # ファイルパス
        data_path = os.path.join(temp_dir, f"data_{session_id}.csv")
        labels_path = os.path.join(temp_dir, f"labels_{session_id}.csv")
        groups_path = os.path.join(temp_dir, f"groups_{session_id}.csv")
        output_path = os.path.join(temp_dir, f"result_{session_id}.csv")
        
        # データをCSVとして保存
        data_df.to_csv(data_path)
        pd.DataFrame({"label": labels.astype(str)}).to_csv(labels_path, index=False)
        pd.DataFrame({"group": groups_use}).to_csv(groups_path, index=False)
        
        # Rスクリプトを作成
        r_script_path = os.path.join(temp_dir, f"script_{session_id}.R")
        r_code = f"""
        # 必要なライブラリをロード
        library(presto)
        
        # データの読み込み
        data <- read.csv("{data_path}", row.names=1)
        labels <- read.csv("{labels_path}")$label
        groups <- read.csv("{groups_path}")$group
        
        # データを確認
        print(dim(data))
        print(length(labels))
        print(groups)
        
        # データの前処理
        data_matrix <- as.matrix(data)
        
        # wilcoxauc の実行
        result <- presto::wilcoxauc(data_matrix, labels, groups_use=groups)
        
        # 結果を保存
        write.csv(result, "{output_path}", row.names=FALSE)
        
        print("処理完了")
        """
        
        with open(r_script_path, "w") as f:
            f.write(r_code)
        
        # Rスクリプトを実行
        try:
            subprocess.run(["Rscript", r_script_path], check=True)
        except subprocess.CalledProcessError as e:
            print(f"Rスクリプト実行エラー: {e}")
            raise
        
        # 結果を読み込む
        if os.path.exists(output_path):
            result_df = pd.read_csv(output_path)
            return result_df
        else:
            raise FileNotFoundError(f"結果ファイルが見つかりません: {output_path}")
            
    finally:
        # 一時ファイルのクリーンアップ（オプション）
        # shutil.rmtree(temp_dir)
        print(f"一時ファイルは {temp_dir} に保存されています")
    
    return None


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

    
def hill_function(x, k, n):
    """
    ベクトル化されたヒル関数の計算
    数式: y = x^n / (k^n + x^n)
    x: 入力配列
    k: ヒル定数
    n: ヒル係数
    """
    x_n = np.power(x, n)
    return x_n / (np.power(k, n) + x_n)

def compute_hill_outer_vectorized(dataL, dataR, k, n):
    """
    dataLとdataRの外積に対してヒル関数を適用する（ベクトル化版）
    dataL: 1次元配列（リガンド発現値）
    dataR: 1次元配列（レセプター発現値）
    k, n: ヒル関数のパラメータ
    戻り値: 計算された2次元配列
    """
    outer = np.outer(dataL, dataR)
    return hill_function(outer, k, n)

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
            db['interaction'].loc[:, 'ligand'] = db['interaction']['ligand'].astype(str)
        if 'receptor' in db['interaction'].columns:
            db['interaction'].loc[:, 'receptor'] = db['interaction']['receptor'].astype(str)



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


def calculate_mean_expression_optimized(data_use, cell_labels, cell_types, min_cells, FunMean):
    """
    Calculate mean expression for ALL cell types (matching R behavior)
    """
    data_use_avg_dict = {}
    cell_counts = {}
    
    # Get cell indices for all cell types
    cell_type_indices = {}
    for cell_type in cell_types:
        indices = np.where(cell_labels == cell_type)[0]
        cell_counts[cell_type] = len(indices)
        cell_type_indices[cell_type] = indices
    
    # Calculate for ALL cell types (don't filter yet)
    for cell_type in cell_types:
        indices = cell_type_indices[cell_type]
        # Perform calculation even for cell types with few cells
        gene_chunk_size = 1000
        n_genes = data_use.shape[1]
        avg_expr = np.zeros(n_genes)
        
        for i in range(0, n_genes, gene_chunk_size):
            end = min(i + gene_chunk_size, n_genes)
            data_subset = data_use[indices, i:end]
            avg_expr[i:end] = np.apply_along_axis(FunMean, 0, data_subset)
        
        data_use_avg_dict[cell_type] = avg_expr
    
    
    return data_use_avg_dict, cell_counts



# logger = logging.getLogger("CellChat")



def apply_mean_function(data_subset, fun_type='triMean', trim=0.1):
    """平均計算関数（Numbaなし）"""
    if fun_type == 'triMean':
        # 四分位数の計算と平均（Numba未使用）
        q1 = np.quantile(data_subset, 0.25, axis=0, method='linear')
        q2 = np.quantile(data_subset, 0.5, axis=0, method='linear')  # median
        q3 = np.quantile(data_subset, 0.75, axis=0, method='linear')
        return (q1 + 2*q2 + q3) / 4
    elif fun_type == 'truncatedMean':
        # トリム平均
        return np.apply_along_axis(lambda x: np.mean(x, axis=0), 0, data_subset)
    else:
        # デフォルトは単純平均
        return np.mean(data_subset, axis=0)

def process_single_permutation(data_use, cluster_indices, cell_types, fun_type='triMean', trim=0.1):
    """単一permutationの全クラスタと全遺伝子の平均発現を計算（Numbaなし）"""
    n_genes = data_use.shape[1]
    numCluster = len(cell_types)
    result = np.zeros((n_genes, numCluster), dtype=np.float32)
    
    # 各クラスターについて処理
    for ct_idx in range(numCluster):
        cells = cluster_indices[ct_idx]
        if len(cells) > 0:
            data_subset = data_use[cells]
            result[:, ct_idx] = apply_mean_function(data_subset, fun_type, trim)
    
    return result

def process_permutation_batch(batch_indices, data_use, cell_labels, permutation, cell_types, 
                             fun_type='triMean', trim=0.1):
    """permutationのバッチ処理関数"""
    n_genes = data_use.shape[1]
    numCluster = len(cell_types)
    results = np.zeros((n_genes, numCluster, len(batch_indices)), dtype=np.float32)
    
    for idx, j in enumerate(batch_indices):
        # j番目のPermutation後の細胞タイプ
        group_boot = cell_labels.values[permutation[:, j]]
        
        # 各クラスターに属する細胞のインデックスを取得
        cluster_indices = [np.where(group_boot == ct)[0] for ct in cell_types]
        
        # すべての遺伝子と細胞タイプの平均発現を計算
        results[:, :, idx] = process_single_permutation(
            data_use, cluster_indices, cell_types, fun_type, trim
        )
    
    return results

def precompute_gene_expressions(data_use, cell_labels, permutation, cell_types, FunMean, nboot, n_jobs=32):
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
    cell_types : list
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
    
    numCluster = len(cell_types)
    n_genes = data_use.shape[1]
    
    # メモリ使用量を制御するためのバッチ処理
    batch_size = min(10, nboot)  # 一度に処理するpermutationの数を調整
    
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
                batch_indices, data_use, cell_labels, permutation, cell_types, fun_type, trim
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
            cluster_indices = [np.where(group_boot == ct)[0] for ct in cell_types]
            
            # すべての遺伝子と細胞タイプの平均発現を計算
            all_gene_expr[:, :, j] = process_single_permutation(
                data_use, cluster_indices, cell_types, fun_type, trim
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

import pandas as pd

def check_gene_symbol(gene_set, gene_info):
    """
    gene_info のインデックスに含まれていない遺伝子があれば、警告を出力する（修正は行わない）
    """
    missing = [gene for gene in gene_set if gene not in gene_info.index]
    if missing:
        print("Warning: The following genes are not in geneInfo:", missing)

def extract_gene_subset(genes, complex_input, gene_info):
    """
    複合体データ（complex_input）に存在する場合、各遺伝子の行からサブユニット情報を抽出する。
    存在しなければ元の遺伝子名をそのまま返す。
    """
    extracted = []
    for gene in genes:
        if gene in complex_input.index:
            row = complex_input.loc[gene]
            # NaN や空文字列を除去してサブユニットを抽出
            subunits = [str(x).strip() for x in row if pd.notna(x) and str(x).strip() != ""]
            if subunits:
                extracted.extend(subunits)
            else:
                extracted.append(gene)
        else:
            extracted.append(gene)
    return list(set(extracted))

def extractGene(db):
    """
    CellChatDB から関与している遺伝子シンボルを抽出する関数
    （R版 extractGene の挙動に合わせる）

    Parameters
    ----------
    db : dict
        以下のキーを含む辞書:
          - "interaction": pandas DataFrame（必須）→ 列 "ligand", "receptor", "agonist", "antagonist", "co_A_receptor", "co_I_receptor を含む
          - "complex": pandas DataFrame（行名に複合体名、各列にサブユニット情報）
          - "cofactor": pandas DataFrame（行名がコファクター遺伝子、"cofactor" で始まる列がサブユニット情報）
          - "geneInfo": pandas DataFrame（公式な遺伝子シンボルをインデックスとして保持）
          
    Returns
    -------
    list
        重複を除いた、関与しているユニークな遺伝子シンボルのリスト
    """
    # 取得方法はユーザー側で以下のように行っています
    # resource = db.get('interaction', pd.DataFrame())
    # complex_input = db.get('complex', pd.DataFrame())
    # cofactor_input = db.get('cofactor', pd.DataFrame())
    # gene_info = db.get('geneInfo', pd.DataFrame())
    resource = db.get('interaction', pd.DataFrame())
    complex_input = db.get('complex', pd.DataFrame())
    cofactor_input = db.get('cofactor', pd.DataFrame())
    gene_info = db.get('geneInfo', pd.DataFrame())

    # 複合体およびコファクター内の全要素について公式遺伝子シンボルかをチェック（警告を出す）
    if not complex_input.empty:
        complex_genes = complex_input.values.flatten()
        check_gene_symbol(list(complex_genes), gene_info)
    if not cofactor_input.empty:
        cofactor_genes = cofactor_input.values.flatten()
        check_gene_symbol(list(cofactor_genes), gene_info)

    # resource から ligand と receptor のユニークな遺伝子を取得
    geneL = resource["ligand"].dropna().astype(str).unique().tolist() if "ligand" in resource.columns else []
    geneR = resource["receptor"].dropna().astype(str).unique().tolist() if "receptor" in resource.columns else []
    geneLR = geneL + geneR

    # resource 内で複合体データに含まれていない遺伝子についてチェック（警告のみ）
    if not complex_input.empty:
        genes_not_in_complex = [gene for gene in geneLR if gene not in complex_input.index]
    else:
        genes_not_in_complex = geneLR
    check_gene_symbol(genes_not_in_complex, gene_info)

    # 複合体情報に基づいて、ligand, receptor の各リストを更新
    if not complex_input.empty:
        geneL = extract_gene_subset(geneL, complex_input, gene_info)
        geneR = extract_gene_subset(geneR, complex_input, gene_info)
    geneLR = geneL + geneR

    # resource から agonist, antagonist, co_A_receptor, co_I_receptor 列を取得し、コファクターのリストを作成
    cofactor_list = []
    for col in ["agonist", "antagonist", "co_A_receptor", "co_I_receptor"]:
        if col in resource.columns:
            cofactor_list.extend(resource[col].dropna().astype(str).tolist())
    cofactor_list = list(set([g for g in cofactor_list if g != ""]))

    # cofactor_input から、行名が cofactor_list に一致する行の "cofactor" で始まる列を抽出
    gene_cofactor = []
    if not cofactor_input.empty:
        cofactor_cols = [col for col in cofactor_input.columns if col.startswith("cofactor")]
        matched_rows = cofactor_input.loc[cofactor_input.index.intersection(cofactor_list)]
        if not matched_rows.empty and cofactor_cols:
            cofactorsubunits = matched_rows[cofactor_cols]
            cofactorsubunits_v = cofactorsubunits.values.flatten()
            gene_cofactor = list(set([str(x).strip() for x in cofactorsubunits_v if pd.notna(x) and str(x).strip() != ""]))

    # 最終的に、interaction 由来の遺伝子と cofactor 由来の遺伝子を統合してユニークなリストを返す
    gene_use = list(set(geneLR + gene_cofactor))
    return gene_use, resource, complex_input, cofactor_input, gene_info

def identify_overexpressed_interactions(features_sig, gene_use, resource, complex_input=None):
    """
    R版のidentifyOverExpressedInteractionsに完全に一致させる実装
    """
    # 型の正規化
    features_sig = [str(f) for f in features_sig]
    gene_use = [str(g) for g in gene_use]
    
    # デバッグ出力
    print(f"DEG genes count: {len(features_sig)}")
    
    # 1. 複合体のフィルタリング - R版と完全に同じロジック
    expressed_complex_names = []
    complexes_with_subunits = {}
    
    if complex_input is not None and not complex_input.empty:
        for complex_name in complex_input.index:
            # R版と同じサブユニット取得ロジック
            subunit_cols = [col for col in complex_input.columns if "subunit" in col]
            subunits = [str(s) for s in complex_input.loc[complex_name, subunit_cols].values if isinstance(s, str) and s != ""]
            
            # R版のロジック：サブユニットのいずれかがDEGで、全てが解析対象遺伝子に含まれる
            has_deg_subunit = any(s in features_sig for s in subunits)
            all_in_gene_use = all(s in gene_use for s in subunits)
            
            if has_deg_subunit and all_in_gene_use:
                expressed_complex_names.append(complex_name)
                complexes_with_subunits[complex_name] = subunits
    
    print(f"Complex subunits with DEGs count: {len(expressed_complex_names)}")
    
    # 2. LRペアのフィルタリング
    valid_elements = set(features_sig).union(set(expressed_complex_names))
    valid_pairs = []
    
    for idx, row in resource.iterrows():
        ligand = str(row['ligand'])
        receptor = str(row['receptor'])
        
        # R版の all(unlist(pairLR[x,]) %in% c(features.sig, rownames(complexSubunits.sig)))
        if ligand in valid_elements and receptor in valid_elements:
            valid_pairs.append(idx)
    
    # フィルタリングされたLRペア
    resource_filtered = resource.loc[valid_pairs].copy() if valid_pairs else resource.iloc[0:0].copy()
    print(f"Filtered LR pairs: {len(resource_filtered)}")
    
    # 結果をCSVにエクスポート（R結果と比較用）
    pd.Series(features_sig).to_csv("py_features_sig.csv")
    pd.Series(expressed_complex_names).to_csv("py_complex_names.csv")
    resource_filtered.to_csv("py_filtered_pairs.csv")
    
    # 3. 関連遺伝子の収集
    lr_related_genes = set(features_sig)
    
    for _, row in resource_filtered.iterrows():
        ligand = str(row['ligand'])
        receptor = str(row['receptor'])
        
        if ligand in gene_use:
            lr_related_genes.add(ligand)
        if receptor in gene_use:
            lr_related_genes.add(receptor)
        
        if ligand in complexes_with_subunits:
            lr_related_genes.update(complexes_with_subunits[ligand])
        if receptor in complexes_with_subunits:
            lr_related_genes.update(complexes_with_subunits[receptor])
    
    return resource_filtered, list(lr_related_genes)



def cellchat_analysis(
    adata,
    groupby,
    #db,
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
   # debug_mode=False,
    trim=0.1,
    apply_pval_filter=True,  # デフォルトをTrueに変更
    features=None,  # 追加：特定の遺伝子リストを使用するオプション
    r_patcher=False,
    use_gpu=True,  # GPU使用フラグ
    gpu_precision='float32',  # GPU精度 ('float32' or 'float64')
    gpu_memory_limit=0.8,  # GPUメモリ使用率制限
    auto_batch_size=True,  # 自動バッチサイズ調整
    do_de=True
):
    """
    最適化されたCellChatアルゴリズム
    
    Parameters
    ----------
    adata : AnnData
        AnnDataオブジェクト
    groupby : str
        細胞タイプ/クラスターを含む.obsのカラム名
    db : dict
        CellChatDB辞書（get_cellchatdb_from_r()で取得したもの）
        'interaction', 'complex', 'cofactor', 'geneInfo'を含む
    
    [その他のパラメータは省略]
    apply_pval_filter : bool, optional
        有意水準に満たない相互作用を削除するかどうか。Default: True（Rの実装と一致）

    features : list, optional
        特定の遺伝子リストを使用する場合に指定。Noneの場合はoverexpressed_genesから特定。
    
    Returns
    -------
    dict
        細胞間通信解析結果
    """
    try:
        logger.info("CellChatアルゴリズムを実行中...")

        progress_bar = st.progress(0)
        status_area = st.empty()
        status_area.text("Initializing GPU and preparing data...")

        # GPU初期化
        use_gpu_actual, has_gpu, cp_module, dtype, memory_info = init_gpu(
            use_gpu=use_gpu, 
            precision=gpu_precision, 
            gpu_memory_limit=gpu_memory_limit
        )
        
        if use_gpu_actual:
            st.success(f"GPU initialized: {memory_info['available_gb']:.1f}GB available")
        else:
            st.info("Using CPU implementation")

        # シードを設定
        np.random.seed(seed)
        if use_gpu_actual and cp_module:
            cp_module.random.seed(seed)


        # 前処理ステップを追加
        if features is not None:
            st.write("Using union/intersection genes...")
            # 提供された遺伝子リストを使用
            logger.info(f"Using provided gene list: {len(features)} genes")
            # 存在する遺伝子のみを保持
          #  features_sig = [f for f in features if f in adata.var_names]
          #  st.write(f"DEG genes found in data: {len(features_sig)}")
            
            # 新しく定義した関数を使用してLRペアと関連遺伝子をフィルタリング

            adata_filtered, resource_filtered = preprocess_data(adata, groupby, complex_input, gene_use=gene_use, min_cells=min_cells,
                thresh_pct=expr_prop, resource=resource, features =features)
            
            if not resource_filtered.empty:
                st.write(resource_filtered.head(3))
        else:
            st.write("Feastures is None")
            # 従来通りoverexpressed_genesを特定
            adata_filtered, resource_filtered = preprocess_data(adata, groupby, complex_input, gene_use=gene_use, min_cells=min_cells,
                thresh_pct=expr_prop, resource=resource, features =None)
        
        
        #以後、reource_filterdをresourceとする
        resource=resource_filtered.copy()
        progress_bar.progress(0.05)
        
       
        # データ形式の検証
        if not hasattr(adata_filtered, 'obs') or not hasattr(adata_filtered, 'var'):
            raise ValueError("入力がAnnDataオブジェクトではありません")
            
        if groupby not in adata_filtered.obs.columns:
            raise ValueError(f"指定されたグループ列 '{groupby}' がadata.obsに存在しません")
        
        # 発現データの取得
        if use_layer is not None and use_layer in adata_filtered.layers:
            logger.info(f"レイヤー '{use_layer}' を使用します")
            X = adata_filtered.layers[use_layer]
        else:
            logger.info("デフォルトXマトリックスを使用します")
            X = adata_filtered.X
        
        # スパース行列の変換（最適化）
        X = optimize_matrix_operations(X, adata_filtered)

        progress_bar.progress(0.1)
        
        # 発現行列をDataFrameに変換
        expr_df = pd.DataFrame(X, index=adata_filtered.obs_names, columns=adata_filtered.var_names)
        
        # 細胞タイプラベルの取得と検証
        cell_labels = adata_filtered.obs[groupby].copy()
        cell_types = np.array(sorted(cell_labels.unique()))
        
        logger.info(f"細胞タイプ数: {len(cell_types)}")
        if len(cell_types) < 2:
            raise ValueError(f"細胞タイプが1つしかありません。少なくとも2つ必要です。")
        
        # メモリ使用量推定とGPU対応の判定
        n_cells, n_genes = X.shape
        memory_estimate = estimate_memory_usage(n_cells, n_genes, nboot, gpu_precision)
        st.info(f"Estimated memory usage: {memory_estimate['total_gb']:.2f}GB")
        
        # 大規模データの場合の判定
        large_data_threshold = 5_000_000  # 500万要素
        is_large_data = n_cells * n_genes > large_data_threshold
        
        if is_large_data and use_gpu_actual:
            st.info(f"Large dataset detected ({n_cells:,} × {n_genes:,}), using GPU acceleration")
        
        # データ型の決定とGPU転送
        if use_gpu_actual and is_large_data:
            # GPU用のデータ準備
            X_cpu = X.astype(dtype if hasattr(dtype, 'dtype') else np.float32)
            
            # GPUメモリ制限チェック
            if memory_estimate['total_gb'] > memory_info.get('usage_limit_gb', 8):
                st.warning(f"Memory usage ({memory_estimate['total_gb']:.1f}GB) exceeds limit, falling back to CPU")
                use_gpu_actual = False
                data_use = X_cpu / np.max(X_cpu)
            else:
                # GPU転送と正規化
                X_gpu = cp_module.array(X_cpu)
                max_val_gpu = cp_module.max(X_gpu)
                data_use_gpu = X_gpu / max_val_gpu
                data_use = cp_module.asnumpy(data_use_gpu)  # 必要に応じてCPUに戻す
                
                # GPUメモリ解放
                del X_gpu, data_use_gpu
                if hasattr(cp_module, 'get_default_memory_pool'):
                    cp_module.get_default_memory_pool().free_all_blocks()
        else:
            # CPU処理
            data_use = X.astype(dtype if hasattr(dtype, 'dtype') else np.float64)
            data_use = data_use / np.max(data_use)
        
        # デバッグ用：データの統計情報を出力
        print(f"Data statistics after max normalization - min: {np.min(data_use)}, max: {np.max(data_use)}, mean: {np.mean(data_use)}")
        print(f"GPU used for preprocessing: {use_gpu_actual and is_large_data}")

        nC = data_use.shape[0]
        progress_bar.progress(0.2)
        status_area.text("Filtering data...")
        # FunMeanの選択 - Rの実装と一致させる
        if type_mean == "triMean":
            def FunMean(x):
                # NANを明示的に除外（na.rm=TRUEに対応）
                x_no_nan = x[~np.isnan(x)]
                if len(x_no_nan) == 0:
                    return np.nan
                # Rと完全に一致させるため、Rの quantile 関数を直接使用
                try:
                    import rpy2.robjects as ro
                    from rpy2.robjects import numpy2ri
                    numpy2ri.activate()
                    r_quantile = ro.r['quantile']
                    r_mean = ro.r['mean']
                    
                    # RでtriMeanを計算
                    r_vector = ro.FloatVector(x_no_nan)
                    r_result = r_mean(r_quantile(r_vector, ro.FloatVector([0.25, 0.5, 0.5, 0.75])))
                    return float(r_result[0])
                except:
                    # フォールバック: Pythonの実装
                    q1 = np.percentile(x_no_nan, 25, interpolation='linear')
                    q2 = np.percentile(x_no_nan, 50, interpolation='linear')
                    q3 = np.percentile(x_no_nan, 75, interpolation='linear')
                    return (q1 + 2*q2 + q3) / 4
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
        
        # 各遺伝子の各細胞タイプにおける平均発現量を計算（最適化）
        data_use_avg_dict, cell_counts = calculate_mean_expression_optimized(
            data_use, cell_labels, cell_types, min_cells, FunMean
        )
        
        logger.info(f"細胞タイプごとの細胞数: {cell_counts}")
        
        # 有効な細胞タイプを確認
      #  if len(cell_types) < 2:
      #      raise ValueError(f"最小細胞数 {min_cells} を満たす細胞タイプが2つ未満です。min_cellsを減らしてください。")
        progress_bar.progress(0.3)
        status_area.text("Calculating LR expression levels...")
        # 平均発現量をDataFrameに変換
        data_use_avg_df = pd.DataFrame(data_use_avg_dict, index=adata_filtered.var_names)
        
        # 遺伝子をインデックスにマッピング
        gene_to_index = {gene: i for i, gene in enumerate(adata_filtered.var_names)}
        
        # リガンドとレセプターの平均発現量を計算
        logger.info("リガンドとレセプターの発現量を計算中...")
        dataLavg = computeExpr_LR(resource['ligand'].values, data_use_avg_df, complex_input)
        dataRavg = computeExpr_LR(resource['receptor'].values, data_use_avg_df, complex_input)
        
        # デバッグ用：発現値をエクスポート
        pd.DataFrame(dataLavg, columns=cell_types, index=resource.index).to_csv("py_dataLavg.csv")
        pd.DataFrame(dataRavg, columns=cell_types, index=resource.index).to_csv("py_dataRavg.csv")
        
        # 共活性化および共阻害受容体の効果を考慮
        dataRavg_co_A_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "A")
        dataRavg_co_I_receptor = computeExpr_coreceptor(cofactor_input, data_use_avg_df, resource, "I")
        dataRavg = dataRavg * dataRavg_co_A_receptor / dataRavg_co_I_receptor
        
        # 細胞数の効果を考慮
        if population_size:
            # 各細胞タイプの細胞数の割合を計算
            cell_proportions = np.array([np.sum(cell_labels == ct) for ct in cell_types]) / nC
            # 各リガンド・レセプターペアに対して同じ割合を使用
            dataLavg2 = np.tile(cell_proportions, (len(resource), 1))
            dataRavg2 = dataLavg2
        else:
            dataLavg2 = np.ones((len(resource), len(cell_types)))
            dataRavg2 = np.ones((len(resource), len(cell_types)))
        
        # アゴニストとアンタゴニストのインデックスを特定
        index_agonist = np.where(resource['agonist'].notna() & (resource['agonist'] != ""))[0] if 'agonist' in resource.columns else []
        index_antagonist = np.where(resource['antagonist'].notna() & (resource['antagonist'] != ""))[0] if 'antagonist' in resource.columns else []
        
        progress_bar.progress(0.4)
        status_area.text("Preparing permutation data...")
        # 置換検定のためのデータを準備
        permutation = np.zeros((nC, nboot), dtype=int)

        # バッチサイズの決定
        if auto_batch_size:
            batch_size = determine_optimal_batch_size(
                nC, n_genes, nboot, use_gpu_actual, 
                memory_info.get('available_gb', 8.0)
            )
        else:
            batch_size = min(100, nboot)
        
        st.info(f"Using batch size: {batch_size} for {nboot} permutations")

        # Permutation生成 - GPU最適化
        if use_gpu_actual and is_large_data and nboot > 50:
            try:
                st.info("Generating permutations on GPU...")
                # 大きなpermutationはバッチごとに生成
                permutation = np.zeros((nC, nboot), dtype=np.int32)
                
                for batch_start in range(0, nboot, batch_size):
                    batch_end = min(batch_start + batch_size, nboot)
                    batch_count = batch_end - batch_start
                    
                    # GPU上でpermutation生成
                    perm_gpu = cp_module.random.randint(0, nC, (nC, batch_count))
                    permutation[:, batch_start:batch_end] = cp_module.asnumpy(perm_gpu)
                    
                    del perm_gpu
                    if batch_start % 500 == 0:  # 定期的にメモリ解放
                        cp_module.get_default_memory_pool().free_all_blocks()
                
                print("Using GPU permutation generation")
                
            except Exception as e:
                print(f"GPU permutation failed, using R/Python fallback: {e}")
                # Rとの完全一致のためにRの乱数生成を使用
                try:
                    permutation = get_r_permutation(nC, nboot, seed=seed)
                    print("Using R permutation for exact compatibility")
                except Exception as e2:
                    print(f"R permutation also failed, using Python fallback: {e2}")
                    for i in range(nboot):
                        permutation[:, i] = np.random.permutation(nC)
        else:
            # 小規模データまたはCPUモード
            try:
                permutation = get_r_permutation(nC, nboot, seed=seed)
                print("Using R permutation for exact compatibility")
            except Exception as e:
                print(f"R permutation failed, using Python fallback: {e}")
                for i in range(nboot):
                    permutation[:, i] = np.random.permutation(nC)
        progress_bar.progress(0.5)
        
        # 全遺伝子の平均発現を事前計算（最適化）
        #all_gene_expr, data_use_avg_boot = precompute_gene_expressions(
        all_gene_expr = precompute_gene_expressions(
            data_use, cell_labels, permutation, cell_types, FunMean, nboot
        )
        
        # コミュニケーション確率と有意性を計算
        numCluster = len(cell_types)
        nLR = len(resource)
        Prob = np.zeros((numCluster, numCluster, nLR))
        Pval = np.zeros((numCluster, numCluster, nLR))
        
        logger.info(f"リガンド-レセプターペアの解析を開始: {nLR}ペア")
        progress_bar.progress(0.6)
        status_area.text("Starting LR pair analysis...")
        
        # 複合体とその構成遺伝子のマッピングを事前に計算
        complex_mapping = precompute_complex_mapping(complex_input, gene_to_index)
        
        # リガンド・レセプターの遺伝子インデックスを取得
        ligand_indices = []
        receptor_indices = []
        
        nLR = len(resource)
        for i in range(nLR):
            ligand = resource['ligand'].iloc[i]
            receptor = resource['receptor'].iloc[i]
            
            # 単一遺伝子か複合体かを判定
            if isinstance(ligand, str) and ligand in gene_to_index:
                # 単一遺伝子
                ligand_indices.append((i, [gene_to_index[ligand]], False))
            elif isinstance(ligand, str) and ligand in complex_mapping:
                # 複合体
                ligand_indices.append((i, complex_mapping[ligand], True))
            else:
                # 未知の遺伝子/複合体
                ligand_indices.append((i, [], None))
            
            if isinstance(receptor, str) and receptor in gene_to_index:
                receptor_indices.append((i, [gene_to_index[receptor]], False))
            elif isinstance(receptor, str) and receptor in complex_mapping:
                receptor_indices.append((i, complex_mapping[receptor], True))
            else:
                receptor_indices.append((i, [], None))
        
        # LRペアごとのループ
        st.write("Long loop processes...")
        progress_bar_lr = st.progress(0)
        for i in range(nLR):
            dataLR = np.outer(dataLavg[i, :], dataRavg[i, :])
            
            # ヒル関数を使用（Rのコードと一致）
            P1 = compute_hill_outer_vectorized(dataLavg[i, :], dataRavg[i, :], k, n)
            
            # アゴニスト効果を計算
            P2 = np.ones((numCluster, numCluster))
            if i in index_agonist:
                data_agonist = computeExpr_agonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P2 = np.outer(data_agonist, data_agonist)
            
            # アンタゴニスト効果を計算
            P3 = np.ones((numCluster, numCluster))
            if i in index_antagonist:
                data_antagonist = computeExpr_antagonist(data_use_avg_df, resource, cofactor_input, i, k, n)
                P3 = np.outer(data_antagonist, data_antagonist)
            
            # 細胞数効果を計算
            P4 = np.ones((numCluster, numCluster))
            if population_size:
                P4 = np.outer(dataLavg2[i, :], dataRavg2[i, :])
            
            # 最終的な確率
            Pnull = P1 * P2 * P3 * P4
            Prob[:, :, i] = Pnull
            
            # 置換検定によるP値計算
            if np.sum(Pnull) == 0:
                # 相互作用がない場合
                Pval[:, :, i] = 1
                continue
            
            Pnull_vec = Pnull.flatten()
            
            # このLRペアのリガンドとレセプターのインデックス情報を取得
            ligand_info = ligand_indices[i]
            receptor_info = receptor_indices[i]
            
            # 発現が取得できない場合はスキップ
            if ligand_info[2] is None or receptor_info[2] is None:
                Pval[:, :, i] = 1
                continue
            
            # 各置換でのブートストラップ確率を計算
            Pboot = np.zeros((numCluster * numCluster, nboot))
            
            # バッチ処理を導入してメモリ使用量を管理
            batch_size = min(20, nboot)
            
            # 並列処理の設定
            n_jobs_to_use = min(n_jobs, os.cpu_count() or 1)
            
            # 並列処理用関数

            def compute_permutation_batch_vectorized(batch_indices):
                batch_results = np.zeros((numCluster * numCluster, len(batch_indices)))
                for idx, j in enumerate(batch_indices):
                    lr_i, ligand_gene_indices, is_ligand_complex = ligand_info
                    # リガンドの発現値取得
                    if not is_ligand_complex:
                        if ligand_gene_indices:
                            ligand_idx = ligand_gene_indices[0]
                            dataLavgB = all_gene_expr[ligand_idx, :, j].reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    else:
                        expr_values = np.array([all_gene_expr[l_idx, :, j] for l_idx in ligand_gene_indices])
                        if expr_values.size > 0:
                            log_values = np.log(expr_values + 1e-10)
                            dataLavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                        else:
                            dataLavgB = np.zeros((1, numCluster))
                    # レセプターの発現値取得
                    # レセプターの発現取得
                    lr_i, receptor_gene_indices, is_receptor_complex = receptor_info
                    if not is_receptor_complex:
                        if receptor_gene_indices:
                            receptor_idx = receptor_gene_indices[0]
                            dataRavgB = all_gene_expr[receptor_idx, :, j].reshape(1, -1)
                        else:
                            dataRavgB = np.zeros((1, numCluster))
                    else:
                        expr_values = np.array([all_gene_expr[r_idx, :, j] for r_idx in receptor_gene_indices])
                        if expr_values.size > 0:
                            log_values = np.log(expr_values + 1e-10)
                            dataRavgB = np.exp(np.mean(log_values, axis=0)).reshape(1, -1)
                        else:
                            dataRavgB = np.zeros((1, numCluster))
                    # 外積を計算し、ベクトル化ヒル関数を適用
                    dataLRB = np.outer(dataLavgB[0, :], dataRavgB[0, :])
                    P1_boot = hill_function(dataLRB, k, n)
                    batch_results[:, idx] = P1_boot.flatten()
                return batch_results
            # バッチ処理で並列計算
            for b_start in range(0, nboot, batch_size):
                b_end = min(b_start + batch_size, nboot)
                batch_indices = list(range(b_start, b_end))
                
                if n_jobs_to_use > 1 and len(batch_indices) > 1:
                    # 並列処理
                    batch_results_list = Parallel(n_jobs=n_jobs_to_use, backend="loky")(
                        delayed(compute_permutation_batch_vectorized)([j]) for j in batch_indices
                    )
                    # 結果の結合
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results_list[j_idx][:, 0]
                else:
                    # 単一スレッド処理
                    batch_results = compute_permutation_batch_vectorized(batch_indices)
                    for j_idx, j in enumerate(batch_indices):
                        Pboot[:, j] = batch_results[:, j_idx]
            
            # p値の計算
            #nReject = np.sum(Pboot > np.expand_dims(Pnull_vec, 1), axis=1)
            # Rの実装に合わせて差分で比較
            if r_patcher:
                # 浮動小数点精度の問題を回避
                nReject = np.sum((Pboot - np.expand_dims(Pnull_vec, 1)) > 1e-10, axis=1)
            else:
                # Rと同じ: rowSums(Pboot - Pnull > 0)
                nReject = np.sum((Pboot - np.expand_dims(Pnull_vec, 1)) > 0, axis=1)
            p = nReject / nboot
            Pval[:, :, i] = p.reshape(numCluster, numCluster)
            progress_bar_lr.progress((i + 1) / nLR)
        

        progress_bar_lr.empty()
        progress_bar.progress(0.8)

        # プロブが0の場合のp値を1に設定 (これは常に行う)
        Pval[Prob == 0] = 1

        # P値フィルタリング前の状態を検証
        print(f"CellChat_analysis Before pval filtering - Prob > 0 count: {np.sum(Prob > 0)}")
        
        # apply_pval_filterパラメータに基づいてp値でのフィルタリングを適用
        if apply_pval_filter:
            if r_patcher:
                #Prob[(Pval >= trim_threshold) > 1.49e-8] = 0
                Prob[Pval >= (trim_threshold - 1.49e-8)] = 0
            else:
                Prob[Pval >= trim_threshold] = 0
        # P値フィルタリング後の状態を検証
        print(f"After pval filtering - Prob > 0 count: {np.sum(Prob > 0)}")


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

        # 結果に名前を設定
        dimnames = [list(cell_types), list(cell_types), list(resource.index)]
        
        # デバッグ用：相互作用数の詳細を出力
        print(f"\n=== Interaction Count Summary ===")
        print(f"Total LR pairs analyzed: {nLR}")
        print(f"LR pairs with Prob > 0 (before p-value filter): {np.sum(np.any(Prob > 0, axis=(0,1)))}")
        print(f"Total interactions with Prob > 0: {np.sum(Prob > 0)}")
        print(f"Total interactions with p-value < {trim_threshold}: {np.sum(Pval < trim_threshold)}")
        print(f"Total significant interactions (Prob > 0 AND p < {trim_threshold}): {np.sum((Prob > 0) & (Pval < trim_threshold))}")
        
        # 細胞タイプペアごとの相互作用数
        interaction_counts = np.zeros((len(cell_types), len(cell_types)))
        for i in range(len(cell_types)):
            for j in range(len(cell_types)):
                if apply_pval_filter:
                    interaction_counts[i, j] = np.sum((Prob[i, j, :] > 0) & (Pval[i, j, :] < trim_threshold))
                else:
                    interaction_counts[i, j] = np.sum(Prob[i, j, :] > 0)
        
        print(f"\nInteraction counts per cell type pair:")
        interaction_df = pd.DataFrame(interaction_counts, index=cell_types, columns=cell_types)
        print(interaction_df)
        print(f"================================\n")
        
        # シグナル経路レベルでの通信確率を計算
        netP = computeCommunProbPathway({"prob": Prob, "pval": Pval}, resource,
            thresh=trim_threshold, apply_pval_filter=apply_pval_filter, r_patcher=r_patcher)
        
        # 中心性指標を計算
        if netP["pathways"] is not None and len(netP["pathways"]) > 0:
            netP["centr"] = {}
            for p_idx in range(len(netP["pathways"])):
                # 各パスウェイの確率行列を取得
                pathway_prob = netP["prob"][:, :, p_idx]
                
                # 単一のパスウェイに対して中心性を計算
                # 3次元配列を作成して計算関数に渡す
                pathway_prob_3d = np.expand_dims(pathway_prob, axis=2)
                netP["centr"][p_idx] = netAnalysis_computeCentrality(pathway_prob_3d)[0]
        
        # 集計ネットワークの中心性も計算
        prob_sum = np.sum(Prob, axis=2)
        prob_sum_3d = np.expand_dims(prob_sum, axis=2)
        net_centr = netAnalysis_computeCentrality(prob_sum_3d)[0]
        
        # 集計ネットワークを計算
        net_summary = aggregateCell_Cell_Communication({"prob": Prob, "pval": Pval},
            cell_types, pval_threshold=trim_threshold, apply_pval_filter=apply_pval_filter,
            r_patcher=r_patcher)
        progress_bar.empty()
        logger.info("CellChat解析が完了しました")
        
        # 結果をデータフレームに変換
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
        
        # GPU メモリ解放
        if use_gpu_actual and cp_module:
            try:
                cp_module.get_default_memory_pool().free_all_blocks()
                print("GPU memory freed")
            except Exception as e:
                print(f"GPU memory cleanup warning: {e}")
        
        # メモリ使用量の報告
        if use_gpu_actual:
            st.success(f"GPU-accelerated analysis completed successfully!")
            st.info(f"Final memory usage: {memory_estimate['total_gb']:.2f}GB")
        
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
            'groupby': groupby,
            'gpu_info': {
                'gpu_used': use_gpu_actual,
                'memory_estimate': memory_estimate,
                'batch_size': batch_size if 'batch_size' in locals() else None
            }
        }
    
    except Exception as e:
        logger.error(f"エラーが発生しました: {str(e)}")
        logger.error(traceback.format_exc())
        
        # エラー時もGPUメモリ解放
        if 'use_gpu_actual' in locals() and use_gpu_actual and 'cp_module' in locals() and cp_module:
            try:
                cp_module.get_default_memory_pool().free_all_blocks()
                print("GPU memory freed after error")
            except:
                pass
        
        return {'error': str(e), 'traceback': traceback.format_exc()}


def geometricMean(expr_values):
    """
    Rの実装に合わせた幾何平均計算
    """
    if expr_values.ndim == 1:
        # 1次元配列の場合
        # Rと同じ動作: log(0) = -Inf, mean with -Inf = -Inf, exp(-Inf) = 0
        with np.errstate(divide='ignore'):
            log_values = np.log(expr_values)
        # -Infが含まれる場合、meanも-Infになる（Rの動作を再現）
        if np.any(np.isneginf(log_values)):
            return 0.0
        else:
            return np.exp(np.mean(log_values))
    else:
        # 2次元配列の場合（列ごとに計算）
        result = np.zeros(expr_values.shape[1])
        for i in range(expr_values.shape[1]):
            with np.errstate(divide='ignore'):
                log_values = np.log(expr_values[:, i])
            if np.any(np.isneginf(log_values)):
                result[i] = 0.0
            else:
                result[i] = np.exp(np.mean(log_values))
        return result

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
def computeCommunProbPathway(net, pairLR_use, thresh=0.05, apply_pval_filter=True, r_patcher=False):
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
    
    # P値フィルタリング前の状態を検証
    print(f"computeCommunProbPathway - Before filtering - prob > 0 count: {np.sum(prob > 0)}")

    # apply_pval_filterパラメータに基づいてp値でのフィルタリングを適用
    if apply_pval_filter:
        if r_patcher:
           # prob[(pval - thresh) >= 1.49e-8] = 0
           prob[pval >= (thresh-1.49e-8)] = 0
        else:
            prob[pval >= thresh] = 0
    # P値フィルタリング後の状態を検証
    print(f"computeCommunProbPathway - After filtering - prob > 0 count: {np.sum(prob > 0)}")
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
def aggregateCell_Cell_Communication(net, cell_types, pval_threshold=0.05, apply_pval_filter=True, r_patcher=False):
    """
    細胞間通信ネットワークを集計
    
    Parameters
    ----------
    net : dict
        通信確率と有意性を含む辞書
    cell_types : list
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
    
    # Python実装でフィルタリング前後の値を比較
    print(f"aggregateCell_Cell_Communication Before filtering - sig_prob > 0 count: {np.sum(prob > 0)}")
    print('pval_threshold')
    print(pval_threshold)
    # apply_pval_filterパラメータに基づいてp値でのフィルタリングを適用
    sig_prob = prob.copy()
    if apply_pval_filter:
        if r_patcher:
            sig_prob[(pval - pval_threshold) >= 1.49e-8] = 0
        else:
            sig_prob[pval >= pval_threshold] = 0
    print(f"aggregateCell_Cell_Communication After filtering - sig_prob > 0 count: {np.sum(sig_prob > 0)}")
    # 異なる閾値での試し
    test_thresholds = [0, 0.0001, 0.001, 0.005, 0.006, 0.007, 0.008, 0.009, 0.0095, 0.0099, 0.01, 0.02, 0.05, 0.1]
    for thresh in test_thresholds:
        count = np.sum(sig_prob > thresh)
        print(f"Count with threshold {thresh}: {count}")


    # strength_matrix と count_matrix の計算箇所に以下を追加
    print(f"sig_prob shape: {sig_prob.shape}")
    print(f"sig_prob min: {np.min(sig_prob)}, max: {np.max(sig_prob)}, sum: {np.sum(sig_prob)}")
    print(f"sig_prob > 0 count: {np.sum(sig_prob > 0)}")
    print(f"sig_prob > 1e-10 count: {np.sum(sig_prob > 1e-10)}")
    print(f"sig_prob > 0.01 count: {np.sum(sig_prob > 0.01)}")
    
    # 各細胞タイプペア間の総インタラクション強度
    strength_matrix = np.sum(sig_prob, axis=2)
    
    # インタラクション数
    print("countmatrix >0")
    count_matrix = np.sum(sig_prob > 0, axis=2) #これだとRに比べて数が多くなる
    print(count_matrix)
    print("countmatrix >1.49e-8")
    count_matrix = np.sum(sig_prob > 1.49e-8, axis=2)
    print(count_matrix)
    print("countmatrix >0.001") #r_patcherがあるときはこれでcleanup
    count_matrix = np.sum(sig_prob > 0.001, axis=2)
    print(count_matrix)
    if not r_patcher:
        count_matrix = np.sum(sig_prob > 0, axis=2)

    # 各リガンド-レセプターペアの寄与度
    lr_contribution = np.sum(np.sum(sig_prob, axis=0), axis=0)
    
    # 各細胞タイプの送受信総量
    outgoing = np.sum(strength_matrix, axis=1)
    incoming = np.sum(strength_matrix, axis=0)
    
    # 行列をDataFrameに変換
    strength_df = pd.DataFrame(strength_matrix, index=cell_types, columns=cell_types)
    count_df = pd.DataFrame(count_matrix, index=cell_types, columns=cell_types)
    
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
        'outgoing': pd.Series(outgoing, index=cell_types),
        'incoming': pd.Series(incoming, index=cell_types),
        'network_centrality': network_centrality,
        'strength_stats': strength_stats
    }


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
        st.error(f"サーキュラープロット作成エラー: {str(e)}")
        st.error(traceback.format_exc())
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
        st.error(f"ドットプロット作成エラー: {str(e)}")
        st.error(traceback.format_exc())
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
        st.error(f"ネットワーク図作成エラー: {str(e)}")
        st.error(traceback.format_exc())
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
        st.error(f"パスウェイ解析プロット作成エラー: {str(e)}")
        st.error(traceback.format_exc())
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"パスウェイ解析プロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Pathway Enrichment Analysis')
        return fig

def preprocess_data(
    adata, 
    group_by,
    complex_input,
    gene_use=None,
    min_cells=10, 
    thresh_pct=0,
    thresh_p=0.05,
    resource=None,
    features=None
):
    """
    CellChatのデータ前処理を行う関数（R版の処理順序に合わせた修正版）
    """
    st.write(f"original data: {adata.shape[0]} cells x {adata.shape[1]} genes")
    
    # 1. 細胞グループの細胞数チェック これはRでは最初には行わない
#    cell_counts = adata.obs[group_by].value_counts()
#    valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
#    if len(valid_groups) < len(cell_counts):
#        st.warning(f"{len(cell_counts) - len(valid_groups)}個の細胞グループが{min_cells}細胞未満のため除外されました")
    
    # 有効なグループの細胞だけを保持
#    adata_filtered = adata[adata.obs[group_by].isin(valid_groups)].copy()
    
    valid_signaling_genes = [g for g in gene_use if g in adata.var_names]
    st.write(f"LR genes found in data: {len(valid_signaling_genes)}")
    # シグナリング遺伝子のみを保持（R版のsubsetData相当の処理）
    adata_filtered = adata[:, valid_signaling_genes].copy()
    st.write(f"LR gene data: {adata_filtered.shape[0]} cells x {adata_filtered.shape[1]} genes")
    
    if features is None:
        # 3. 高発現遺伝子の同定（R版のidentifyOverExpressedGenes相当）
        overexpressed_genes_result = identify_overexpressed_genes(
            adata_filtered,  # すでにシグナリング遺伝子にフィルタリング済み
            group_by=group_by,
            thresh_pct=thresh_pct,
            thresh_p=thresh_p,
            only_pos=True,
            do_de=True,
            do_fast=True,
            min_cells=min_cells
        )
        
        # 高発現遺伝子を抽出
        features_sig = overexpressed_genes_result['features']
    else:
        st.write("Using features...")
        features_sig = features.copy()
    
    # 4. LRペアフィルタリング（R版のidentifyOverExpressedInteractions相当）
    # タプルを受け取る
    result = identify_overexpressed_interactions(
        features_sig=features_sig, 
        gene_use=valid_signaling_genes,
        resource=resource,
        complex_input=complex_input
    )
    
    # タプルを展開
    resource_filtered, lr_genes_from_function = result
    
    st.write(f"Filtered LR pairs: {len(resource_filtered)}")
    
    # リソースが空でない場合のみ表示
    if not resource_filtered.empty:
        st.write(resource_filtered.head(3))
        
    # 5. 発現差のある遺伝子＋フィルタリングされたLRペアに関連する遺伝子のみを保持
    lr_related_genes = set(features_sig)  # DEG遺伝子を基本として保持
    
    for _, row in resource_filtered.iterrows():  # ここを修正
        ligand, receptor = row['ligand'], row['receptor']
        # 単純な遺伝子の場合はそのまま追加
        if ligand in valid_signaling_genes:
            lr_related_genes.add(ligand)
        if receptor in valid_signaling_genes:
            lr_related_genes.add(receptor)
        
        # 複合体の場合はサブユニットを追加
        if complex_input is not None:
            if ligand in complex_input.index:
                subunit_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[ligand, subunit_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in valid_signaling_genes]
                lr_related_genes.update(subunits)
            
            if receptor in complex_input.index:
                subunit_cols = [col for col in complex_input.columns if 'subunit' in col]
                subunits = complex_input.loc[receptor, subunit_cols].dropna().astype(str)
                subunits = [s for s in subunits if s != "" and s in valid_signaling_genes]
                lr_related_genes.update(subunits)
    
    # 関連する遺伝子のみを保持
    lr_related_genes = list(lr_related_genes)
    adata_filtered = adata_filtered[:, lr_related_genes].copy()
    
    st.write(f"前処理後: {adata_filtered.shape[0]}細胞 x {adata_filtered.shape[1]}遺伝子")
    
    return adata_filtered, resource_filtered


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
        st.session_state.cellchat_temp_dir = True
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
    st.markdown("###### weight(strength)はRとほぼ一致。interaction countは相関するが甘め。")

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
                st.markdown("##### Cells:")
                st.write(", ".join(cell_list))
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
                col1, col2 = st.columns(2)
                
                with col1:
                    species = st.radio("Species:", ('mouse', 'human'), 
                                      index=check_species_index(adata.var.index.to_list()[:50]))
                    
                    data_layer = st.selectbox("Using layer:", 
                                             ['X (default)'] + available_layers,
                                             index=0)
                    
                    selected_types = st.multiselect(
                        "Signaling types to analyze (can choose multiple):",
                        ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling"],
                        default=["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"]
                    ) 
                    n_cpus = st.slider("Num CPUs:", 
                                  min_value=1, max_value=os.cpu_count(), step=1, value=os.cpu_count()-2)
                    r_patcher = st.checkbox("Fine-tune to match R calc results", value=False,
                        help="Rとの計算結果の違いをある程度吸収する。weightについてはチェックしなくても近似する。countはチェックしないと多くなるが、傾向は一致する。countのthresholdの基準はやや恣意的。")

                with col2:
                    st.markdown("#### Do not change the following parameters unless necessary.")
                    min_cells = st.number_input('Min cells in a cell-type:', 
                                               min_value=0, max_value=100, step=1, value=10,
                                               help="SCALAのcellchatではフィルタリングはoff。")

                    expr_prop = st.number_input('Min fraction of expressing cells:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.0)

                    thresh_p = st.number_input('Threshold p value for overexpressed genes:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.05,
                                               help = "Wilcoxonによる変化遺伝子選択のthreshold")

                    apply_pval_filter = st.checkbox("Filter pathways by P value", value=True)

                    trim_threshold = st.number_input('Filtering P threshold:', 
                                               min_value=0.01, max_value=0.20, step=0.01, value=0.05)
                
                    n_perms = st.slider("Permutation number:", 
                                min_value=20, max_value=500, step=20, value=100)

                    population_size = st.checkbox("population.size", value=False, help="細胞数が多いクラスター同士の通信ほど強く（確率が高めに）評価する場合。")
                
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
                    st.write(cellchatdb['interaction'].head(3))
                    
                    # データレイヤーの設定
                    use_layer = None if data_layer == 'X (default)' else data_layer

                    #ここでdbを展開しておく
                    gene_use, resource, complex_input, cofactor_input, gene_info = extractGene(cellchatdb)
                    # ここで、resource はフィルタリング済みの interaction（つまり、selected_types に該当するペア）のみを含むはずです。
                    st.write(f"Resource LR pair numbers:{len(resource)}")
                    print(f'gene_use {len(gene_use)}')
                    print(f'complex_input {len(complex_input)}')
                    print(f'cofactor_input {len(cofactor_input)}')
                    print(f'gene_info {len(gene_info)}')


                    # データベースの検証
                    if resource.empty:
                        raise ValueError("DBの相互作用情報(interaction)が空です。有効なCellChatDBか確認してください。")
                        
                    # 必要な列の存在を確認
                    for required_col in ['ligand', 'receptor']:
                        if required_col not in resource.columns:
                            raise ValueError(f"リソースデータフレームには '{required_col}' 列が必要です")

                    st.write(f"LR genes in db: {len(gene_use)}")

                    with st.spinner('CellChat calculation...'):
                        result = cellchat_analysis(
                            adata,
                            groupby=groupby,
                          #  db=cellchatdb,
                            gene_use=gene_use,
                            complex_input=complex_input,
                            cofactor_input=cofactor_input,
                            resource=resource,
                            use_layer=use_layer,
                            min_cells=min_cells,
                            expr_prop=expr_prop,
                         #   log_scale=log_transformed,
                            pseudocount=1.0,
                            k=0.5,
                            trim_threshold=trim_threshold,  
                            nboot=n_perms,
                            seed=1,
                            n_jobs=n_cpus,
                            trim=0.1,
                            apply_pval_filter=apply_pval_filter,
                            r_patcher=r_patcher,
                            population_size=population_size
                        )
                        
                        st.session_state.cellchat_res = result
                    
                    if 'error' in result:
                        st.error(f"解析中にエラーが発生しました: {result['error']}")
                        st.code(result['traceback'])
                    else:
                        st.success('CellChat解析が完了しました！')
                        time_cpu = time.time() - start_cpu
                        st.write(f"Computing time: {round(time_cpu)}")
                        print(f"Computing time: {round(time_cpu)}")
                
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
                    

                    # まず、結果のCSV（相互作用テーブル）を文字列に変換
                    csv_str = result['results'].to_csv(index=False, sep = '\t')

                    # 次に、strength_matrixとcount_matrixをTSVに変換
                    # DataFrameの場合
                    if isinstance(result['network']['strength_matrix'], pd.DataFrame):
                        strength_tsv = result['network']['strength_matrix'].to_csv(sep='\t', index=True)
                    else:
                        # numpy配列の場合、DataFrame化してからTSV化（ここでは仮にdimnamesがあると仮定）
                        strength_df = pd.DataFrame(result['network']['strength_matrix'], 
                                                   index=result['network']['dimnames'][0], 
                                                   columns=result['network']['dimnames'][1])
                        strength_tsv = strength_df.to_csv(sep='\t', index=True)

                    if isinstance(result['network']['count_matrix'], pd.DataFrame):
                        count_tsv = result['network']['count_matrix'].to_csv(sep='\t', index=True)
                    else:
                        count_df = pd.DataFrame(result['network']['count_matrix'], 
                                                index=result['network']['dimnames'][0], 
                                                columns=result['network']['dimnames'][1])
                        count_tsv = count_df.to_csv(sep='\t', index=True)

                    # メモリ上でZIPファイルを作成
                    zip_buffer = io.BytesIO()
                    with zipfile.ZipFile(zip_buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as zf:
                        zf.writestr("cellchat_results.tsv", csv_str)
                        zf.writestr("strength_matrix.tsv", strength_tsv)
                        zf.writestr("count_matrix.tsv", count_tsv)

                    # バッファのポインタを先頭に戻す
                    zip_buffer.seek(0)

                    # st.download_buttonでZIPファイルを提供
                    st.download_button(
                        label="結果をZIPとしてダウンロード",
                        data=zip_buffer,
                        file_name="cellchat_results.zip",
                        mime="application/zip"
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
                st.markdown("###### Visualization options are displayed at the bottom of the left side panel")
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
                    # カラーマップ選択部分の後に追加
                    if 'cellchat_res' in st.session_state and st.session_state.cellchat_res is not None:
                        result = st.session_state.cellchat_res
                        # カラーマップが変更された場合または初めての場合
                        if 'cell_color_map' not in st.session_state or st.session_state.get('current_cmap', '') != cmap_name:
                            st.session_state.cell_color_map = create_cell_color_mapping(cell_list, cmap_name)
                            st.session_state.current_cmap = cmap_name

                    sort_cell = st.checkbox("Change cell order?")
                    if sort_cell:
                        with st.form("Sorter"):
                            #sorted_order = sort_items(sorted(adata.obs[groupby].unique()))
                            sorted_order = sort_items(cell_list.copy())
                            submitted_sort = st.form_submit_button("Done sorting")
                    else:
                        sorted_order = None
                
                tabs = st.tabs([
                    "Heatmap", 
                    "Chord", 
                    "LR contrib", 
                    "Dot", 
                    "Centrality", 
                    "Role sctr",
                    "Sig contrib",
                    "Circle",   # New tab for circle plots
                    "Sig role",   # New tab for signaling role analysis 
                    "Expression"     # New tab for gene expression analysis
                ])

                # Keep your existing tab implementations
                with tabs[0]:
                    st.markdown("#### Cell interaction heatmap")
                    
                    col1, col2 = st.columns(2)
                    
                    with col1:
                        
                        heatmap_type = st.radio(
                            "Heatmap type:",
                            ["Interaction strength", "Interaction number"],
                            horizontal=True
                        )
                        heatmap_annot = st.checkbox("Show values?", value=True)
                        show_color_bar = st.checkbox("Show cell-type color bars?", value=False)

                    with col2:
                        heatmap_title = st.slider("Title font size:", min_value=0, max_value=30, value=20)
                        heatmap_font = st.slider("Font size:", min_value=0, max_value=30, value=14)
                        heatmap_x = st.slider("Fig width:", min_value=1.0, max_value=20.0, value=10.0)
                        heatmap_y = st.slider("Fig height:", min_value=1.0, max_value=20.0, value=8.0)


                    if st.button("Generate heatmap"):
                        
                        
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
                                font_size=heatmap_font, 
                                font_size_title=heatmap_title,
                                sorted_order=sorted_order,
                                annot=heatmap_annot,
                                color_use=st.session_state.get('cell_color_map', None),  # 色マッピングを追加
                                show_color_bar=show_color_bar,
                                figx=heatmap_x,
                                figy=heatmap_y

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

                        col1, col2 = st.columns(2)

                        with col1:
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
                        
                        with col2:
                            # 表示カスタマイズオプション
                            sector_space = st.number_input("Space between sectors (degree):", min_value=0, max_value=20, value=5)
                            alpha_edge = st.slider("Transparency of edge:", min_value=0.1, max_value=1.0, value=0.4, step=0.1)
                            edge_border_width = st.slider("Edge border width", min_value=0.0, max_value=1.0, value=0.4, step=0.1)
                           # show_edge_border = st.checkbox("Show edge border", value=True)
                            chord_x = st.slider("Fig X size:", min_value=4.0, max_value=20.0, value=10.0, step=0.5)
                            chord_y = st.slider("Fig Y size:", min_value=4.0, max_value=20.0, value=10.0, step=0.5)
                        
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
                                            color_use=st.session_state.get('cell_color_map', None),
                                            sorted_order=sorted_order,
                                            measure=measure,
                                            figsize=(chord_x, chord_y),
                                            space=sector_space,               # セクター間の空白 5度
                                            alpha_edge=alpha_edge,        # リンクの透明度 0.6
                                        #    show_edge_border=show_edge_border, # リンクの黒い輪郭を非表示
                                            edge_border_width=edge_border_width
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

                            if st.button("Generate LR contribution graph"):
                            
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

                        if st.button("Generate dotplot"):
                            
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

                        st.subheader("In/Out activity")
                        
                        if 'outgoing' in result['network'] and 'incoming' in result['network']:
                            outgoing = result['network']['outgoing']
                            incoming = result['network']['incoming']
                            
                            # 送受信データをマージ
                            if isinstance(outgoing, pd.Series) and not outgoing.empty:

                                if st.button("Generate activity plot"):
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
                        use_count = st.checkbox("Dot size reflects link number", value=True)
                        cell_groups = None
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
                      #  title=" " #なにか入れるとplot重なって見えなくなるが、Noneだとg画面でうまく表示されない
                        xlabel = st.text_input("X軸ラベル:", value="Outgoing interaction strength" if x_measure == "outdeg" else x_measure)
                        ylabel = st.text_input("Y軸ラベル:", value="Incoming interaction strength" if y_measure == "indeg" else y_measure)
                        role_x = st.slider("Fig width:", min_value=1, max_value=20, value=6)
                        role_y = st.slider("Fig height:", min_value=1, max_value=20, value=4)
                    
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
                                    color_use=st.session_state.get('cell_color_map', None),
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
                                    do_label=do_label,
                                    width=role_x,
                                    height=role_y,
                                    use_count=use_count
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
                                    sorted_order=sorted_order,
                                    cmap_name=cmap_name
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


                        col1, col2 = st.columns(2)

                        with col1:
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
                                    sorted(pathways),
                                    default=[pathways[0]] if pathways else []
                                )
                                selected_pathway = selected_pathways  # リストとして渡す
                                measure_key = "weight"  # 特定パスウェイでは常にweight
                            circle_type = st.radio(
                                    "Circle plot type:",
                                    ["Total network", "Individual cell-type"],
                                    horizontal=True
                                )

                        with col2:
                            # 表示カスタマイズオプション
                            edge_width_max = st.slider("Max edge width:", min_value=1, max_value=20, value=8, step=1)
                            circle_alpha_edge = st.slider("Edge transparency:", min_value=0.1, max_value=1.0, value=0.6, step=0.1)
                            vertex_size_max= st.slider("Max node size:", min_value=1, max_value=15, value=6, step=1)
                            show_vertex = st.checkbox("Node size reflects the number of cells", value=False)

                    else:
                        st.warning("パスウェイ情報が利用できません")
                        selected_pathway = None
                        measure_key = "weight"


                    
                    # 必要に応じて vertex_weight も取得（例：result['adata'].obs から各細胞の数）
                    vertex_weight = None
                    if hasattr(result['adata'].obs, 'value_counts') and show_vertex:
                        try:
                            cell_counts = result['adata'].obs[groupby].value_counts()
                          #  st.write(cell_counts)
                            vertex_weight = [cell_counts.get(ct, 1) for ct in result['net']['dimnames'][0]]
                          #  st.write(vertex_weight)
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
                                                pathway_name=selected_pathway,
                                                measure=measure_key,
                                                sorted_order=sorted_order,
                                                edge_width_max=edge_width_max,
                                                alpha_edge=circle_alpha_edge,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                color_use=st.session_state.get('cell_color_map', None)  # 色マッピングを追加
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                measure=measure_key,
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                sorted_order=sorted_order,
                                                edge_width_max=edge_width_max,
                                                alpha_edge=circle_alpha_edge,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None) 
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
                                                edge_width_max=edge_width_max,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                alpha_edge=circle_alpha_edge,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None) 
                                            )
                                        else:
                                            fig_circle = netVisual_circle_individual(
                                                net=result,
                                                pathway_name=selected_pathway,
                                                measure="weight",  # 個別の場合は通常 weight を使用
                                                title_name=title,
                                                cmap_name=cmap_name,
                                                edge_width_max=edge_width_max,
                                                sorted_order=sorted_order,
                                                vertex_weight=vertex_weight,
                                                vertex_size_max=vertex_size_max,
                                                alpha_edge=circle_alpha_edge,
                                                arrow=True,
                                                color_use=st.session_state.get('cell_color_map', None) 
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
                        hide_color_bar = not st.checkbox("Show cell color bar?", value=False)
                    with col2:
                        font_size = st.slider("Font size:", min_value=8, max_value=14, value=10)
                    
                        # ヒートマップのサイズ設定
                        width = st.slider("Fig width:", min_value=6, max_value=15, value=10)
                        height = st.slider("Fig height:", min_value=1.0, max_value=8.0, value=3.0, step=0.5)

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
                                    show_value=show_value_role,
                                    cmap_name=cmap_name,
                                    hide_color_bar=hide_color_bar,
                                    color_use=st.session_state.get('cell_color_map', None)
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
        st.info("h5adファイルかpklファイルをアップロードして開始してください。")
    