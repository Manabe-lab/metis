import liana as li
# needed for visualization and toy data
import scanpy as sc

from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean, scseqcomm
from liana.mt import rank_aggregate
from functools import lru_cache
import itertools
from plotnine import scale_y_discrete

import numpy as np
import streamlit as st
import pandas as pd
import numpy as np
import os
from helper_func import clear_old_directories
from helper_func import clear_old_files
from helper_func import check_species_index
import time
import matplotlib.pyplot as plt

import datetime

from scipy import stats
from statsmodels.stats.multitest import multipletests
from functools import partial
import concurrent.futures

# we import plotnine
import plotnine as p9
from plotnine import scale_color_cmap

st.set_page_config(page_title="Liana_comparison", page_icon="💬")

from plotnine import scale_y_discrete
import matplotlib.pyplot as plt

import gc  # Add this at the top with other imports

#from pycirclize import Circos


# Define the function dictionary
function_dict = {
    'singlecellsignalr': li.method.singlecellsignalr,
    'connectome': li.method.connectome,
    'cellphonedb': li.method.cellphonedb,
    'natmi': li.method.natmi,
    'logfc': li.method.logfc,
    'cellchat': li.method.cellchat,
    'geometric_mean': li.method.geometric_mean,
    'rank_aggregate': li.mt.rank_aggregate,
    'scseqcomm': li.method.scseqcomm
}

def compare_conditions(adata, condition_col, condition1, condition2, groupby, resource, expr_prop, func_name, res_name):
    """
    Compare LR interactions between two conditions
    
    Parameters:
    -----------
    adata : AnnData object
        Contains the single-cell data
    condition_col : str
        Column name in adata.obs containing condition information
    condition1, condition2 : str
        Names of conditions to compare
    groupby : str
        Column name for cell type information
    resource : pandas DataFrame
        LR database to use
    expr_prop : float
        Expression proportion threshold
    func_name : str
        Name of LIANA method to use (should match keys in function_dict)
    res_name : str
        Key for storing results
        
    Returns:
    --------
    dict containing results for both conditions and differential analysis
    """
    # Split data by condition
    adata1 = adata[adata.obs[condition_col] == condition1].copy()
    adata2 = adata[adata.obs[condition_col] == condition2].copy()
    
    # Convert method name to match function_dict keys if necessary
    if func_name == 'Rank_Aggregate':
        func_name = 'rank_aggregate'
    elif func_name == 'CellPhoneDB':
        func_name = 'cellphonedb'
    elif func_name == 'CellChat':
        func_name = 'cellchat'
    elif func_name == 'Connectome':
        func_name = 'connectome'
    elif func_name == 'SingleCellSignalR':
        func_name = 'singlecellsignalr'
    elif func_name == 'Geometric Mean':
        func_name = 'geometric_mean'
    
    # Run LIANA for each condition
    rankaggregate_options = {
        'groupby': groupby, 
        'resource': resource, 
        "expr_prop": expr_prop,
        "use_raw": False, 
        "key_added": res_name, 
        "n_jobs": 8
    }
    
    # Calculate for condition 1
    if func_name not in function_dict:
        raise ValueError(f"Unknown method: {func_name}. Available methods: {list(function_dict.keys())}")
    
    function_dict[func_name](adata1, **rankaggregate_options)
    res1 = adata1.uns[res_name].copy()
    res1['condition'] = condition1
    
    # Calculate for condition 2
    function_dict[func_name](adata2, **rankaggregate_options)
    res2 = adata2.uns[res_name].copy()
    res2['condition'] = condition2
    
    # Compare results
    merged_res = pd.concat([res1, res2])
    
    # Calculate differential scores
    diff_res = calculate_differential_scores(res1, res2)
    
    return {
        'condition1': res1,
        'condition2': res2,
        'merged': merged_res,
        'differential': diff_res
    }



def plot_volcano(stats_df, padj_threshold=0.05, log2fc_threshold=1):
    """
    Create volcano plot of differential interactions
    """
    # Check if we have any results
    if len(stats_df) == 0:
        st.warning("No statistical results available for plotting")
        return None
        
    # Create significance categories
    stats_df['significant'] = 'NS'
    stats_df.loc[(stats_df['padj'] < padj_threshold) & 
                (abs(stats_df['mean_diff']) > log2fc_threshold), 'significant'] = 'Significant'
    
    # Create plot
    plot = (p9.ggplot(stats_df, p9.aes(x='mean_diff', y='-np.log10(padj)', 
                                      color='significant'))
           + p9.geom_point(alpha=0.6)
           + p9.geom_hline(yintercept=-np.log10(padj_threshold), linetype='dashed')
           + p9.geom_vline(xintercept=[-log2fc_threshold, log2fc_threshold], linetype='dashed')
           + p9.scale_color_manual(values={'NS': 'grey', 'Significant': 'red'})
           + p9.labs(x='Log2 Fold Change', y='-log10(adjusted p-value)')
           + p9.theme_bw())
    
    return plot


def calculate_differential_scores(res1, res2):
    """
    Calculate differential scores between two conditions
    """
    # Create unique identifier for each interaction
    def create_interaction_id(df):
        return df.apply(lambda x: f"{x['source']}|{x['target']}|{x['ligand_complex']}|{x['receptor_complex']}", axis=1)
    
    res1 = res1.copy()
    res2 = res2.copy()
    
    res1['interaction_id'] = create_interaction_id(res1)
    res2['interaction_id'] = create_interaction_id(res2)
    
    # Merge results
    diff_res = res1.merge(res2, on='interaction_id', suffixes=('_1', '_2'))
    
    # Calculate log2 fold change of magnitude scores
    magnitude_cols = ['lr_means', 'expr_prod', 'lrscore', 'magnitude_rank', 'lr_gmeans', 'inter_score', 'lr_probs']
    for col in magnitude_cols:
        if col in diff_res.columns:
            col1 = f"{col}_1"
            col2 = f"{col}_2"
            # Add small constant to avoid division by zero
            diff_res[f"{col}_log2fc"] = np.log2((diff_res[col1] + 1e-10) / (diff_res[col2] + 1e-10))
    
    return diff_res

# UI modifications for streamlit
def add_comparison_ui():
    """Add UI elements for condition comparison"""
    st.markdown("### Condition Selection")
    condition_col = st.selectbox("Column containing condition information:", 
                               meta, 
                               index=find_first_index_or_zero(meta, ["condition", "treatment"]))
    
    conditions = sorted(adata.obs[condition_col].unique())
    condition1 = st.selectbox("Condition 1:", conditions, index=0)
    condition2 = st.selectbox("Condition 2:", conditions, index=min(1, len(conditions)-1))
    
    return condition_col, condition1, condition2

# Modified visualization functions for comparison
def plot_differential_dotplot(diff_res, magnitude_col, specificity_col, source_labels, target_labels, 
                            dot_filter=0.05, dot_color='viridis', figure_size=(8, 7)):
    """Plot differential dotplot comparing conditions"""
    return li.pl.dotplot(
        liana_res=diff_res,
        colour=f"{magnitude_col}_log2fc",
        size=specificity_col,
        source_labels=source_labels,
        target_labels=target_labels,
        figure_size=figure_size,
        filter_fun=lambda x: x[specificity_col] <= dot_filter
    )


def single_permutation(data_tuple):
    """
    単一のpermutationによる差分を計算する関数
    
    Parameters:
    -----------
    data_tuple : tuple
        (expr_matrix, condition_labels, seed)
    
    Returns:
    --------
    float : permutationによる差分
    """
    # データの展開
    expr_matrix, condition_labels, seed = data_tuple
    
    # シードの設定
    np.random.seed(seed)
    
    # ラベルのシャッフル
    shuffled_labels = np.random.permutation(condition_labels)
    
    # シャッフルしたラベルでの差分計算
    mean1 = np.mean(expr_matrix[shuffled_labels == 0])
    mean2 = np.mean(expr_matrix[shuffled_labels == 1])
    
    return mean2 - mean1  # condition2 - condition1 の差分を返す

def permutation_test_with_liana(adata, condition_col, condition1, condition2, groupby,
                              resource, expr_prop, func_name, res_name, n_permutations=1000):
    """
    Perform permutation test using full LIANA recalculation
    """
    # Calculate observed result
    observed_result = {}
    for condition in [condition1, condition2]:
        subset = adata[adata.obs[condition_col] == condition].copy()
        rankaggregate_options = {
            'groupby': groupby,
            'resource': resource,
            "expr_prop": expr_prop,
            "use_raw": False,
            "key_added": res_name,
            "n_jobs": 8
        }
        function_dict[func_name](subset, **rankaggregate_options)
        observed_result[condition] = subset.uns[res_name]
    
    # Prepare permutation data
    permutation_data = [(adata, condition_col, groupby, resource, expr_prop, 
                        func_name, res_name, i) for i in range(n_permutations)]
    
    # Run permutations in parallel
    permuted_results = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        futures = [executor.submit(single_permutation, data) for data in permutation_data]
        
        # Create progress bar
        progress_bar = st.progress(0.0)
        st.write(f"Running {n_permutations} permutations...")
        
        try:
            for i, future in enumerate(concurrent.futures.as_completed(futures)):
                try:
                    result = future.result()
                    if result is not None:
                        permuted_results.append(result)
                except Exception as e:
                    st.warning(f"Error in permutation: {e}")
                
                # Update progress
                progress = (i + 1) / n_permutations
                progress_bar.progress(progress)
        finally:
            progress_bar.empty()
    
    return observed_result, permuted_results

def calculate_permutation_pvalues(observed_result, permuted_results, magnitude_col):
    """
    Calculate p-values from permutation results
    """
    # Extract observed values for each interaction
    observed_values = {}
    for interaction in observed_result['condition1'].groupby(
        ['source', 'target', 'ligand_complex', 'receptor_complex']):
        key = tuple(interaction[0])
        if key in observed_result['condition2'].groups:
            obs_val1 = observed_result['condition1'].get_group(key)[magnitude_col].mean()
            obs_val2 = observed_result['condition2'].get_group(key)[magnitude_col].mean()
            observed_values[key] = obs_val2 - obs_val1  # KO - WT
    
    # Calculate p-values
    pvalues = {}
    for key, obs_diff in observed_values.items():
        # Extract permuted values for this interaction
        perm_diffs = []
        for perm_result in permuted_results:
            if key in perm_result.groups:
                group = perm_result.get_group(key)
                n_half = len(group) // 2
                perm_diff = (group[magnitude_col].iloc[n_half:].mean() - 
                           group[magnitude_col].iloc[:n_half].mean())
                perm_diffs.append(perm_diff)
        
        if perm_diffs:
            # Calculate two-sided p-value
            pvalue = np.mean(np.abs(perm_diffs) >= np.abs(obs_diff))
            pvalues[key] = {
                'mean_diff': obs_diff,
                'pvalue': pvalue
            }
    
    return pd.DataFrame.from_dict(pvalues, orient='index')

def process_single_permutation(perm_data):
    """
    単一のpermutationを処理する関数
    
    Parameters:
    -----------
    perm_data: tuple
        (adata, condition_col, groupby, resource, expr_prop, func_name, res_name, seed)
    """
    (adata_orig, condition_col, groupby, resource, expr_prop, 
     func_name_converted, res_name, magnitude_col, seed) = perm_data
    
    # データのコピーと乱数シードの設定
    np.random.seed(seed)
    adata_perm = adata_orig.copy()
    
    try:
        # ラベルのシャッフル
        adata_perm.obs[condition_col] = np.random.permutation(adata_perm.obs[condition_col].values)
        
        # データの分割
        perm_adata1 = adata_perm[adata_perm.obs[condition_col] == adata_orig.obs[condition_col].unique()[0]].copy()
        perm_adata2 = adata_perm[adata_perm.obs[condition_col] == adata_orig.obs[condition_col].unique()[1]].copy()
        
        # 両条件でLIANA実行
        results = []
        for perm_adata in [perm_adata1, perm_adata2]:
            rankaggregate_options = {
                'groupby': groupby,
                'resource': resource,
                "expr_prop": expr_prop,
                "use_raw": False,
                "key_added": f"{res_name}_perm_{seed}",
                "n_jobs": 1
            }
            function_dict[func_name_converted](perm_adata, **rankaggregate_options)
            results.append(perm_adata.uns[f"{res_name}_perm_{seed}"][magnitude_col].mean())
        
        # 差分を計算
        perm_diff = results[1] - results[0]
        return perm_diff
        
    except Exception as e:
        print(f"Error in permutation {seed}: {str(e)}")
        return None

def process_cell_pair_parallel(cell_pair_data):
    """各細胞ペアのpermutation testを並列実行"""
    name, group1, group2, magnitude_col, n_permutations, adata, condition_col, groupby, resource, expr_prop, func_name, res_name = cell_pair_data
    
    # Calculate observed difference
    mean_diff = np.mean(group2[magnitude_col]) - np.mean(group1[magnitude_col])
    
    # Convert method name
    func_name_converted = func_name.lower()
    if func_name == 'Rank_Aggregate':
        func_name_converted = 'rank_aggregate'
    elif func_name == 'CellPhoneDB':
        func_name_converted = 'cellphonedb'
    elif func_name == 'CellChat':
        func_name_converted = 'cellchat'
    elif func_name == 'Connectome':
        func_name_converted = 'connectome'
    elif func_name == 'SingleCellSignalR':
        func_name_converted = 'singlecellsignalr'
    elif func_name == 'Geometric Mean':
        func_name_converted = 'geometric_mean'
    
    # Permutationsの並列実行
    n_cpus = os.cpu_count()
    max_workers = max(1, n_cpus - 1)  # 1コアは他の処理用に残す
    
    perm_data_list = [
        (adata, condition_col, groupby, resource, expr_prop, 
         func_name_converted, res_name, magnitude_col, i) 
        for i in range(n_permutations)
    ]
    
    perm_diffs = []
    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(process_single_permutation, data) for data in perm_data_list]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
                if result is not None:
                    perm_diffs.append(result)
            except Exception as e:
                print(f"Error in permutation: {str(e)}")
    
    # Calculate p-value
    if perm_diffs:
        pvalue = np.mean(np.abs(perm_diffs) >= np.abs(mean_diff))
    else:
        pvalue = np.nan
        
    return {
        'source': name[0],
        'target': name[1],
        'ligand': name[2],
        'receptor': name[3],
        'mean_diff': mean_diff,
        'pvalue': pvalue
    }

def get_significant_pairs(adata, groupby, resource, expr_prop, func_name, res_name, threshold=0.05):
    """
    単一条件でのLR解析を行い、adjusted P値で有意なペアを抽出
    """
    # LIANA実行
    rankaggregate_options = {
        'groupby': groupby, 
        'resource': resource, 
        "expr_prop": expr_prop,
        "use_raw": False, 
        "key_added": res_name, 
        "n_jobs": 8
    }
    
    # Convert method name
    func_name_converted = func_name.lower()
    if func_name == 'Rank_Aggregate':
        func_name_converted = 'rank_aggregate'
    elif func_name == 'CellPhoneDB':
        func_name_converted = 'cellphonedb'
    elif func_name == 'CellChat':
        func_name_converted = 'cellchat'
    elif func_name == 'Connectome':
        func_name_converted = 'connectome'
    elif func_name == 'SingleCellSignalR':
        func_name_converted = 'singlecellsignalr'
    elif func_name == 'Geometric Mean':
        func_name_converted = 'geometric_mean'
    
    function_dict[func_name_converted](adata, **rankaggregate_options)
    
    # 有意性の判定に使用するカラムを決定
    specificity_col = pval_col[method_list.index(func_name)]
    if specificity_col == "None":
        specificity_col = lr_means_list[method_list.index(func_name)]
    
    # 結果を取得
    result = adata.uns[res_name].copy()
    
    # P値の処理とFDR補正
    if specificity_col in result.columns:
        # P値を数値として処理
        p_values = pd.to_numeric(result[specificity_col], errors='coerce')
        
        # NANとInfinity値の処理
        p_values = p_values.fillna(1.0)
        p_values = np.minimum(p_values, 1.0)  # 1以上の値を1に制限
        p_values = np.maximum(p_values, 0.0)  # 負の値を0に制限
        
        # FDR補正
        mt_results = multipletests(p_values, method='fdr_bh', alpha=threshold)
        padj = mt_results[1]  # インデックス1にadjusted p-valuesが格納されている
        result['padj'] = padj
        
        # デバッグ情報を表示
 #       st.write(f"Original P-values range: [{p_values.min():.2e}, {p_values.max():.2e}]")
 #       st.write(f"Adjusted P-values range: [{padj.min():.2e}, {padj.max():.2e}]")
        
        # 補正前後で有意なペアの数を比較
        n_sig_before = sum(p_values < threshold)
        n_sig_after = sum(padj < threshold)
 #       st.write(f"Significant pairs before FDR: {n_sig_before}")
 #       st.write(f"Significant pairs after FDR: {n_sig_after}")
        
        # Adjusted P値でフィルタリング
        significant = result[result['padj'] < threshold]
        st.write(f"Found {len(significant)} pairs significant at adjusted P < {threshold}")
        
        # ペアのセットを作成
        pairs = set(zip(
            significant['source'],
            significant['target'],
            significant['ligand_complex'],
            significant['receptor_complex']
        ))
    else:
        st.warning(f"Specificity column '{specificity_col}' not found in results")
        pairs = set()
    
    return pairs

def compare_conditions_with_stats(adata, condition_col, condition1, condition2, groupby, 
                                resource, expr_prop, func_name, res_name, n_permutations=1000,
                                initial_padj_threshold=0.05):
    """Compare LR interactions between conditions using parallel permutation test"""
    # 各条件のデータを分割
    adata1 = adata[adata.obs[condition_col] == condition1].copy()
    adata2 = adata[adata.obs[condition_col] == condition2].copy()
    
    # 各条件で有意なペアを取得
    with st.spinner("Analyzing condition 1..."):
        pairs1 = get_significant_pairs(
            adata1, groupby, resource, expr_prop, func_name, res_name, 
            threshold=initial_padj_threshold
        )
    
    with st.spinner("Analyzing condition 2..."):
        pairs2 = get_significant_pairs(
            adata2, groupby, resource, expr_prop, func_name, res_name,
            threshold=initial_padj_threshold
        )
    
    # 有意なペアの和集合を取得
    significant_pairs = pairs1.union(pairs2)
    st.write(f"Total unique significant pairs: {len(significant_pairs)}")
    
    if len(significant_pairs) == 0:
        st.warning("No significant pairs found in either condition. Try adjusting the adjusted P-value threshold.")
        return None

    # メモリ解放
    del adata1, adata2
    gc.collect()
    
    # 通常の比較結果を取得
    comparison_results = compare_conditions(adata, condition_col, condition1, condition2, 
                                         groupby, resource, expr_prop, func_name, res_name)
    
    # magnitude columnの決定
    magnitude_col = lr_means_list[method_list.index(func_name)]
    if magnitude_col == "None":
        magnitude_col = pval_col[method_list.index(func_name)]
    
    # 有意なペアのみをグループ化
    res1_groups = comparison_results['condition1'].groupby(
        ['source', 'target', 'ligand_complex', 'receptor_complex'])
    res2_groups = comparison_results['condition2'].groupby(
        ['source', 'target', 'ligand_complex', 'receptor_complex'])
    
    statistical_results = []
    total_pairs = len(significant_pairs)
    
    # プログレスバーの作成
    progress_bar = st.progress(0)
    progress_text = st.empty()
    
    try:
        # 各ペアについて処理
        for i, pair in enumerate(significant_pairs):
            # 進捗状況の表示
            progress_text.write(f"Processing pair {i+1}/{total_pairs}: {pair[0]} -> {pair[1]}")
            
            if (pair in res1_groups.groups and pair in res2_groups.groups):
                group1 = res1_groups.get_group(pair)
                group2 = res2_groups.get_group(pair)
                
                # Permutation testのためのデータ準備
                expr_matrix = np.concatenate([
                    group1[magnitude_col].values,
                    group2[magnitude_col].values
                ])
                condition_labels = np.array([0] * len(group1) + [1] * len(group2))
                
                # 観測された差分の計算
                mean_diff = np.mean(group2[magnitude_col]) - np.mean(group1[magnitude_col])
                
                # Permutation testの実行
                try:
                    # CPUコア数の取得
                    n_cpus = os.cpu_count()
                    max_workers = max(1, n_cpus - 1)  # 1コアは他の処理用に残す
                    
                    perm_data_list = [
                        (expr_matrix, condition_labels, i) 
                        for i in range(n_permutations)
                    ]
                    
                    perm_diffs = []
                    with concurrent.futures.ProcessPoolExecutor(max_workers=max_workers) as executor:
                        # 並列処理でpermutationを実行
                        futures = [executor.submit(single_permutation, data) for data in perm_data_list]
                        for future in concurrent.futures.as_completed(futures):
                            result = future.result()
                            if result is not None:
                                perm_diffs.append(result)
                    
                    if perm_diffs:
                        pvalue = np.mean(np.abs(perm_diffs) >= np.abs(mean_diff))
                    else:
                        pvalue = np.nan
                        
                except Exception as e:
                    st.warning(f"Error in permutation test for pair {pair}: {str(e)}")
                    pvalue = np.nan
                
                # 結果の保存
                result = {
                    'source': pair[0],
                    'target': pair[1],
                    'ligand': pair[2],
                    'receptor': pair[3],
                    'mean_diff': mean_diff,
                    'pvalue': pvalue
                }
                statistical_results.append(result)
            
            # 進捗バーの更新
            progress_bar.progress((i + 1) / total_pairs)
            
            # メモリの解放
            if (i + 1) % 10 == 0:
                gc.collect()
    
    finally:
        # クリーンアップ
        progress_bar.empty()
        progress_text.empty()
    
    # 結果のDataFrame作成
    stats_df = pd.DataFrame(statistical_results)
    
    # FDR補正
    if len(stats_df) > 0 and not stats_df['pvalue'].isna().all():
        # Fill NA values with 1 for FDR correction
        stats_df['pvalue'] = pd.to_numeric(stats_df['pvalue'], errors='coerce').fillna(1)
        stats_df['padj'] = multipletests(stats_df['pvalue'], method='fdr_bh')[1]
    
    # 最終結果の保存
    comparison_results['statistics'] = stats_df
    return comparison_results


 

# Streamlitインターフェースの拡張
def add_statistical_ui():
    """Add UI elements for statistical analysis"""
    st.sidebar.markdown("### Statistical Analysis Options")
    initial_padj_threshold = st.sidebar.number_input("Adj P threshold for initial analysis", 
                                           min_value=0.001, 
                                           max_value=1.000, 
                                           value=0.050, 
                                           step=0.001,
                                           format='%.3f')
    n_permutations = st.sidebar.number_input("Number of permutations", 
                                           min_value=10, 
                                           max_value=10000, 
                                           value=1000, 
                                           step=100)
    padj_threshold = st.sidebar.number_input("Adjusted p-value threshold", 
                                           min_value=0.0, 
                                           max_value=1.0, 
                                           value=0.05, 
                                           step=0.01)
    log2fc_threshold = st.sidebar.number_input("Log2 fold change threshold", 
                                             min_value=0.0, 
                                             value=0.0, 
                                             step=0.1)
    
    return n_permutations, padj_threshold, log2fc_threshold, initial_padj_threshold




@st.cache_resource
def load_geneinfo():
	geneinfo_human = pd.read_csv("db/nichenetr.db/geneinfo_human.tsv", sep='\t')
	geneinfo_2022 = pd.read_csv("db/nichenetr.db/geneinfo_2022.tsv", sep='\t')
	return geneinfo_human, geneinfo_2022

@st.cache_data
def read_h5ad(file):
	adata = sc.read_h5ad(file)
	return adata

@st.cache_data
def read_map_df():
	map_df = li.rs.get_hcop_orthologs(url='https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz',
			  columns=['human_symbol', 'mouse_symbol'],
			   # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
			   min_evidence=3
			   )
	# rename the columns to source and target, respectively for the original organism and the target organism
	map_df = map_df.rename(columns={'human_symbol':'source', 'mouse_symbol':'target'})
	return map_df


def find_first_index_or_zero(lst, elements):
	for element in elements:
		try:
			return lst.index(element)
		except ValueError:
			continue
	return 0

def has_nan_like(lst):
	nan_values = {'nan', 'NaN', 'NAN', 'n/a', 'N/A', '', ' '}
	return any(str(item).strip().lower() in nan_values for item in lst)

@st.cache_data
def convert_human_to_mouse_symbols(symbols, version=1): # nichenetrのRスクリプトをClaude3.5で変換
	if not isinstance(symbols, (list, pd.Series)):
		raise ValueError("symbols should be a list or pandas Series of human gene symbols")
	if version == 1:
		geneinfo = geneinfo_human
	elif version == 2:
		geneinfo = geneinfo_2022
	else:
		raise ValueError("version must be 1 or 2")
	unambiguous_mouse_genes = (
		geneinfo.dropna()
		.groupby('symbol_mouse').size()
		.reset_index(name='count')
		.query('count < 2')['symbol_mouse']
		.tolist()
	)
	ambiguous_mouse_genes = (
		geneinfo.dropna()
		.groupby('symbol_mouse').size()
		.reset_index(name='count')
		.query('count >= 2')['symbol_mouse']
		.tolist()
	)
	geneinfo_ambiguous_solved = geneinfo[
		(geneinfo['symbol_mouse'].isin(ambiguous_mouse_genes)) &
		(geneinfo['symbol'] == geneinfo['symbol_mouse'].str.upper())
	]
	geneinfo = pd.concat([
		geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
		geneinfo_ambiguous_solved
	]).dropna()
	humansymbol2mousesymbol = dict(zip(geneinfo['symbol'], geneinfo['symbol_mouse']))
	converted_symbols = [humansymbol2mousesymbol.get(symbol, np.nan) for symbol in symbols]
	return converted_symbols


@lru_cache(maxsize=None)
def cached_convert_human_to_mouse_symbols(symbol):
	result = convert_human_to_mouse_symbols([symbol])
	return result[0] if result else np.nan


@st.cache_data
def mouse_conversion(resource):
 #   st.write("Converting to mouse genes... This may take a moment...")

	def convert_symbol(symbol):
		if pd.isna(symbol):
			return np.nan
		if '_' in symbol:
			parts = symbol.split('_')
			converted = [cached_convert_human_to_mouse_symbols(part) for part in parts]
			return '_'.join(converted) if all(not pd.isna(part) for part in converted) else np.nan
		return cached_convert_human_to_mouse_symbols(symbol)

	# Vectorize the convert_symbol function
	vectorized_convert = np.vectorize(convert_symbol, otypes=[object])

	# Convert all symbols at once
	mouse_ligands = vectorized_convert(resource['ligand'])
	mouse_receptors = vectorized_convert(resource['receptor'])

	# Create a DataFrame with both human and mouse symbols
	result = pd.DataFrame({
		'ligand': resource['ligand'],
		'receptor': resource['receptor'],
		'mouse_ligand': mouse_ligands,
		'mouse_receptor': mouse_receptors
	})

	# Apply has_nan_like to each row
	has_nan = result[['mouse_ligand', 'mouse_receptor']].isna().any(axis=1)

	# Split into converted and non-converted
	nonconversion = result[has_nan]
	mouse = result[~has_nan][['mouse_ligand', 'mouse_receptor']].rename(
		columns={'mouse_ligand': 'ligand', 'mouse_receptor': 'receptor'}
	)

	# Drop duplicates from mouse DataFrame
	mouse = mouse.drop_duplicates()

	return mouse, nonconversion


@st.cache_data
def gnerate_interaction_matrix(df, LR_count_threshold):
	def calculate_interaction_metrics(df, source_cell, target_cell, LR_count_threshold):
		# Filter the data for the specific source and target cell types
			filtered_df = df[(df['source'] == source_cell) & (df['target'] == target_cell)]
			# Calculate the total interaction strength
			interaction_strength = filtered_df[magnitude_col].sum()
			# Count the number of unique ligand-receptor pairs
			lr_count = len(filtered_df[filtered_df[specificity_col]<LR_count_threshold])
			return interaction_strength, lr_count

	# Get unique cell types
	cell_types = set(df['source'].unique()) | set(df['target'].unique())

	# Calculate interaction metrics for all combinations
	results = []
	for source, target in itertools.product(cell_types, repeat=2):
		strength, count = calculate_interaction_metrics(df, source, target, LR_count_threshold)
		results.append({
			'source': source,
			'target': target,
			'interaction_strength': strength,
			'lr_pair_count': count
			})

	# Convert results to a DataFrame
	results_df = pd.DataFrame(results)

	# Sort results by interaction strength in descending order
	results_df = results_df.sort_values('interaction_strength', ascending=False)
	return results_df

@st.cache_data
def calc_lr(_adata,rankaggregate_options, func_name):
	function_dict[func_name](adata,**rankaggregate_options)
	return adata


if "liana_res" not in st.session_state:
	st.session_state.liana_res = []

#if "liana_basic_change" not in st.session_state:
#	st.session_state.liana_basic_change = False

if "liana_temp_dir" not in st.session_state:
	st.session_state.liana_temp_dir = True
	liana_temp_dir = "temp/" + str(round(time.time()))
	if not os.path.exists('temp'):
		os.mkdir('temp')
	else:
		clear_old_directories("temp")
		clear_old_files("temp")
	os.mkdir(liana_temp_dir)
	st.session_state.liana_temp_dir = liana_temp_dir
else:
	liana_temp_dir = st.session_state.liana_temp_dir

st.markdown("### L-R interaction calculation using LIANA+ in one condition")

function_dict = {
'singlecellsignalr': singlecellsignalr,
'connectome': connectome,
'cellphonedb': cellphonedb,
'natmi': natmi,
'logfc': logfc,
'cellchat': cellchat,
'geometric_mean': geometric_mean,
'rank_aggregate': rank_aggregate,
'scseqcomm': scseqcomm
}

method_list = ['CellPhoneDB','Connectome','log2FC','NATMI','SingleCellSignalR','Rank_Aggregate','GeometricMean','scSeqComm','CellChat']
func_list = ['cellphonedb','connectome','logfc','natmi','singlecellsignalr','rank_aggregate','geometric_mean','scseqcomm','cellchat']
lr_means_list = ['lr_means','expr_prod','None','expr_prod','lrscore','magnitude_rank','lr_gmeans','inter_score','lr_probs']
pval_col = ['cellphone_pvals','scaled_weight','lr_logfc','spec_weight','None','specificity_rank','gmean_pvals','None','cellchat_pvals']


#============ 新しいファイルをアップロードしたときは、cacheをclearする

def get_file_identifier(file):
	if file is not None:
		return f"{file.name}_{file.size}"
	return None


uploaded_file = st.file_uploader("Upload a h5ad file", type=['h5ad'])

if uploaded_file  is not None:
	current_file_id = get_file_identifier(uploaded_file)

	if 'last_file_id' not in st.session_state:
		st.session_state.last_file_id = None

	if current_file_id != st.session_state.last_file_id:
		st.cache_data.clear()
		st.cache_resource.clear()
		st.session_state.last_file_id = current_file_id
		st.success("新しいファイルが検出されました。キャッシュをクリアしました。")

#---------------

	adata = read_h5ad(uploaded_file)
	st.write("Uploaded data")
	st.write(adata.X[:3,:3])

	meta = adata.obs.columns.to_list()
	for i in ['nFeature_RNA','nCount_RNA','percent.mt', 'Cell_id']:
		try:
			meta.remove(i)
		except:
			pass
	submitted_basic = False

	groupby = st.selectbox("cell identity:", meta, index = find_first_index_or_zero(meta, ["cell.ident",
	"seurat_clusters"]))

	sampleby = st.selectbox("sample/condition:", meta, index = find_first_index_or_zero(meta, ["orig.ident",
	"sample","KO"]))

	cell_list = sorted(adata.obs[groupby].cat.categories.to_list())

	sample_list = sorted(adata.obs[sampleby].cat.categories.to_list())

	if len(sample_list) != 2:
		st.warning("The number of conditions is not 2.")
		st.stop()

	st.write(f"Cell types: {', '.join(cell_list)}")
	st.write(f"Samples/Conditions: {', '.join(sample_list)}")
	with st.form("Basic settings:"):
		species = st.radio("Species:", ('mouse','human'), index = check_species_index(adata.var.index.to_list()[:50]))
		method = st.selectbox(
		"Method:",
		('Rank_Aggregate','CellPhoneDB','CellChat','Connectome','log2FC','NATMI','SingleCellSignalR','Geometric Mean','scSeqComm'), index = 0)

		func_name = func_list[method_list.index(method)]

		if species == "mouse":
			db_default = 1
		else:
			db_default = 0
		db = st.selectbox("database:", ('consensus', 'mouseconsensus', 'cellphonedb','cellchatdb', 'cellchat_secreted_signaling', 'baccin2019', 'cellcall', 'cellinker', 'celltalkdb', 'connectomedb2020',  'embrace', 'guide2pharma', 'hpmr', 'icellnet', 'italk', 'kirouac2010', 'lrdb','ramilowski2015'), index = db_default)
		st.write("LR DBs are in human symbols except mouseconsensus. For mouse, they will be converted.")

		expr_prop = st.number_input('Threshold for min fraction of cells expressing the gene', min_value =0.0, max_value=0.9,
		step =0.01, value=0.1)
		remove_complex = st.checkbox("Remove ligand/receptor complexes (e.g., L17A_IL17F, IL17RA_IL17RC)?")
		st.write("Note that some LR pairs have only comlexes.")
		submitted_basic = st.form_submit_button("Change settings")
	if submitted_basic:
		st.session_state.liana_res = []

	if db == 'cellchat_secreted_signaling': #cellchat secreted signalingを追加
		resource = pd.read_csv("db/Cellchat_secretory.tsv", sep = '\t')
		if species == 'mouse':
			geneinfo_human, geneinfo_2022 = load_geneinfo()
			mouse, unconverted = mouse_conversion(resource)
			if len(unconverted) >0:
				st.write("Uncorverted LR")
				st.write(unconverted)
			st.write("Original LR numbers: " + str(len(resource)))
			st.write("Mouse LR numbers: " + str(len(mouse)))
			resource = mouse
	elif db != 'mouseconsensus' and species == 'mouse':
		geneinfo_human, geneinfo_2022 = load_geneinfo()
#		map_df = read_map_df()
		resource = li.rs.select_resource(db)
		mouse, unconverted = mouse_conversion(resource)
		if len(unconverted) >0:
			st.write("Uncorverted LR")
			st.write(unconverted)
		st.write("Original LR numbers: " + str(len(resource)))
		st.write("Mouse LR numbers: " + str(len(mouse)))
		resource = mouse
#		st.write(mouse.head(100))
	else:
		resource = li.rs.select_resource(db)
		st.write("LR numbers: " + str(len(resource)))

	if remove_complex:
		resource = resource[~resource['ligand'].str.contains('_') & ~resource['receptor'].str.contains('_')]
		st.write("Filtered LR numbers: " + str(len(resource)))

	if db == "consensus":
		res_name = "liana" + "_res"
	else:
		res_name = db + "_res"

	st.write(" ")

	n_permutations, padj_threshold, log2fc_threshold, initial_padj_threshold = add_statistical_ui()

	# Streamlitでの使用
# Streamlitでの使用
	if st.button('Run analysis'):
	    with st.spinner('Running LIANA analysis and permutation tests...'):
	        results = compare_conditions_with_stats(
	            adata, 
	            condition_col=sampleby,
	            condition1=sample_list[0],
	            condition2=sample_list[1],
	            groupby=groupby,
	            resource=resource,
	            expr_prop=expr_prop,
	            func_name=method,
	            res_name='liana_res',
	            n_permutations=n_permutations,
	            initial_padj_threshold = initial_padj_threshold
	        )

	    if results is None:
	        st.error("Analysis could not be completed due to insufficient significant pairs.")
	        st.stop()

	    if not isinstance(results, dict) or 'statistics' not in results:
	        st.error("Invalid results format. Please check the analysis parameters.")
	        st.stop()

	    if len(results['statistics']) == 0:
	        st.warning("No statistical results were generated. Try adjusting the significance thresholds.")
	        st.stop()


	    # Ensure p-values are numeric and perform FDR correction
	    stats_df = results['statistics'].copy()
	    if 'pvalue' in stats_df.columns:
	        # Fill NA values with 1 for FDR correction
	        stats_df['pvalue'] = pd.to_numeric(stats_df['pvalue'], errors='coerce').fillna(1)
	        # Perform FDR correction
	        stats_df['padj'] = multipletests(stats_df['pvalue'], method='fdr_bh')[1]
	    else:
	        st.error("No p-values found in results. Check the statistical analysis.")
	        st.stop()

	    st.write("### Analysis Results")
	    
	    # タブで結果を整理
	    tab1, tab2, tab3 = st.tabs(["Summary", "Significant Interactions", "Visualization"])
	    
	    with tab1:
	        st.write("#### Summary Statistics")
	        total_interactions = len(stats_df)
	        significant = len(stats_df[
	            (stats_df['padj'] < padj_threshold) & 
	            (abs(stats_df['mean_diff']) > log2fc_threshold)
	        ])
	        
	        col1, col2, col3 = st.columns(3)
	        with col1:
	            st.metric("Total Interactions", total_interactions)
	        with col2:
	            st.metric("Significant Interactions", significant)
	        if total_interactions > 0:  # Avoid division by zero
	            with col3:
	                st.metric("Significance Rate", f"{(significant/total_interactions*100):.1f}%")
	    
	    with tab2:
	        st.write("#### Significant Interactions")
	        # 有意な結果のフィルタリングと表示
	        significant_interactions = stats_df[
	            (stats_df['padj'] < padj_threshold) & 
	            (abs(stats_df['mean_diff']) > log2fc_threshold)
	        ].copy()
	        
	        if len(significant_interactions) > 0:
	            # mean_diffの符号に基づいて条件を追加
	            significant_interactions['regulation'] = significant_interactions['mean_diff'].apply(
	                lambda x: f"Up in {sample_list[1]}" if x > 0 else f"Up in {sample_list[0]}"
	            )
	            
	            # DataFrameの表示オプション
	            st.write(f"Showing interactions with:")
	            st.write(f"- Adjusted p-value < {padj_threshold}")
	            st.write(f"- |Log2 fold change| > {log2fc_threshold}")
	            
	            # ソートオプション
	            sort_columns = [col for col in ['padj', 'mean_diff', 'source', 'target', 'ligand', 'receptor'] 
	                          if col in significant_interactions.columns]
	            sort_by = st.selectbox("Sort by:", sort_columns)
	            significant_interactions = significant_interactions.sort_values(sort_by)
	            
	            # 表示する列の選択
	            display_columns = [col for col in 
	                ['source', 'target', 'ligand', 'receptor', 'mean_diff', 'pvalue', 'padj', 'regulation'] 
	                if col in significant_interactions.columns]
	            
	            # 結果の表示
	            st.dataframe(
	                significant_interactions[display_columns].style.format({
	                    'mean_diff': '{:.3f}',
	                    'pvalue': '{:.2e}',
	                    'padj': '{:.2e}'
	                })
	            )
	            
	            # 結果のダウンロードボタン
	            csv = significant_interactions.to_csv(index=False)
	            st.download_button(
	                label="Download significant interactions as CSV",
	                data=csv,
	                file_name="significant_interactions.csv",
	                mime="text/csv"
	            )
	        else:
	            st.warning("No significant interactions found with current thresholds.")
	    
	    with tab3:
	        st.write("#### Visualizations")
	        
	        # Volcano plot
	        st.write("##### Volcano Plot")
	        if len(stats_df) > 0:
	            volcano_plot = plot_volcano(
	                stats_df,
	                padj_threshold=padj_threshold,
	                log2fc_threshold=log2fc_threshold
	            )
	            if volcano_plot is not None:
	                st.pyplot(volcano_plot.draw())
	                
	                # Volcanoプロットの保存ボタン
	                pdf_path = f"{liana_temp_dir}/volcano_plot.pdf"
	                volcano_plot.save(pdf_path, bbox_inches='tight')
	                with open(pdf_path, "rb") as pdf_file:
	                    st.download_button(
	                        label="Download volcano plot (PDF)",
	                        data=pdf_file,
	                        file_name="volcano_plot.pdf",
	                        mime="application/pdf"
	                    )
	        else:
	            st.warning("Not enough data points to create visualization.")