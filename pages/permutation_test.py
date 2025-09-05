import streamlit as st
import pandas as pd
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector
import pyper
import os
import time
import sys
import re
from collections import Counter
from helper_func import clear_old_directories, clear_old_files, remove_after_space, remove_sample_num
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
import plotly.express as px
import plotly.graph_objects as go

st.set_page_config(page_title="Permutation_Test", page_icon="🔀", layout="wide")

st.title("🔀 Permutation Test Analysis")
st.markdown("### Freedman-Lane & Simple Permutation Tests")

# 読み込み関数（Calc_DESeq2.pyから流用）
@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

def clean_column_names_for_r(df):
    """Convert DataFrame column names to R-safe format"""
    original_columns = df.columns.tolist()
    cleaned_columns = []
    mapping = {}
    
    for col in original_columns:
        cleaned = re.sub(r'[^a-zA-Z0-9_]', '.', str(col))
        if cleaned and cleaned[0].isdigit():
            cleaned = 'X' + cleaned
        cleaned = re.sub(r'\.+', '.', cleaned)
        cleaned = cleaned.rstrip('.')
        
        cleaned_columns.append(cleaned)
        if col != cleaned:
            mapping[col] = cleaned
    
    df.columns = cleaned_columns
    return df, mapping, original_columns

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []    
    min_length = min(len(s) for s in strings)
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break            
    if suffix_length == 0:
        return strings        
    return [s[:-suffix_length] for s in strings]

def rename_duplicates(df):
    """Rename duplicate indices by adding _2, _3, etc."""
    lis = df.index.values
    counts = Counter()
    new_indices = []
    
    for x in lis:
        counts[x] += 1
        if counts[x] == 1:
            new_indices.append(x)
        else:
            new_indices.append(f"{x}_{counts[x]}")
    
    if len(lis) != len(set(lis)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("Original duplicated entries: " + str(len(lis) - len(set(lis))))
        df.index = new_indices
        st.write("Renamed indices with suffixes _2, _3, etc.")
        
    return df

# Permutation test functions
def simple_permutation_test(y, x, n_permutations=1000, stat_type='t_stat'):
    """
    Simple permutation test
    """
    observed_stat = calculate_statistic(y, x, stat_type)
    permuted_stats = []
    
    for _ in range(n_permutations):
        # ラベルをシャッフル
        x_permuted = np.random.permutation(x)
        perm_stat = calculate_statistic(y, x_permuted, stat_type)
        permuted_stats.append(perm_stat)
    
    permuted_stats = np.array(permuted_stats)
    
    # Two-tailed p-value with continuity correction
    if stat_type in ['t_stat', 'correlation']:
        p_value = (np.sum(np.abs(permuted_stats) >= np.abs(observed_stat)) + 1) / (len(permuted_stats) + 1)
    else:
        p_value = (np.sum(permuted_stats >= observed_stat) + 1) / (len(permuted_stats) + 1)
    
    return observed_stat, permuted_stats, p_value

def freedman_lane_test(y, x_full, x_reduced, n_permutations=1000, stat_type='t_stat'):
    """
    Freedman-Lane permutation test
    """
    # Reduced modelで回帰
    reg_reduced = LinearRegression()
    reg_reduced.fit(x_reduced, y)
    y_pred_reduced = reg_reduced.predict(x_reduced)
    
    # 残差計算
    residuals = y - y_pred_reduced
    
    # Full modelで観測統計量計算
    observed_stat = calculate_full_model_statistic(y, x_full, stat_type)
    
    permuted_stats = []
    
    for _ in range(n_permutations):
        # 残差をシャッフル
        residuals_permuted = np.random.permutation(residuals)
        
        # 新しい応答変数
        y_star = y_pred_reduced + residuals_permuted
        
        # Full modelで統計量計算
        perm_stat = calculate_full_model_statistic(y_star, x_full, stat_type)
        permuted_stats.append(perm_stat)
    
    permuted_stats = np.array(permuted_stats)
    
    # Two-tailed p-value with continuity correction
    if stat_type in ['t_stat', 'correlation']:
        p_value = (np.sum(np.abs(permuted_stats) >= np.abs(observed_stat)) + 1) / (len(permuted_stats) + 1)
    else:
        p_value = (np.sum(permuted_stats >= observed_stat) + 1) / (len(permuted_stats) + 1)
    
    return observed_stat, permuted_stats, p_value

def permuco_freedman_lane_test(y, group, batch, n_permutations=1000, seed=42):
    """
    Run permuco::lmperm() Freedman-Lane test
    """
    try:
        # Import R packages
        ro.r('library(permuco)')
        
        # Convert to R objects  
        ro.globalenv['y_r'] = ro.FloatVector(y)
        ro.globalenv['group_r'] = ro.IntVector(group.astype(int))
        ro.globalenv['batch_r'] = ro.IntVector(batch.astype(int))
        
        # Run permuco test
        ro.r(f'''
        set.seed({seed})
        
        df <- data.frame(
            y = y_r,
            group = as.factor(group_r),
            batch = as.factor(batch_r)
        )
        
        # Run Freedman-Lane test
        perm_result <- lmperm(y ~ batch + group, data = df, 
                             np = {n_permutations}, 
                             method = 'freedman_lane',
                             rnd_rotation = y ~ batch)
        
        # Extract group effect results
        table_data <- perm_result$table
        group_row <- which(rownames(table_data) == 'group1')
        
        if(length(group_row) > 0) {{
            permuco_t_stat <- table_data[group_row, 't value']
            permuco_p_value <- table_data[group_row, 'resampled Pr(>|t|)']
            permuco_estimate <- table_data[group_row, 'Estimate']
            permuco_std_error <- table_data[group_row, 'Std. Error']
        }} else {{
            # Fallback to last row
            last_row <- nrow(table_data)
            permuco_t_stat <- table_data[last_row, 't value']
            permuco_p_value <- table_data[last_row, 'resampled Pr(>|t|)']
            permuco_estimate <- table_data[last_row, 'Estimate']
            permuco_std_error <- table_data[last_row, 'Std. Error']
        }}
        
        # Store additional information
        permuco_table <- table_data
        extraction_success <- TRUE
        ''')
        
        # Check if extraction was successful
        success = ro.r('extraction_success')[0]
        
        if success:
            t_stat = float(ro.r('permuco_t_stat')[0])
            p_value = float(ro.r('permuco_p_value')[0])
            estimate = float(ro.r('permuco_estimate')[0])
            std_error = float(ro.r('permuco_std_error')[0])
            
            return {
                'success': True,
                't_statistic': t_stat,
                'p_value': p_value,
                'estimate': estimate,
                'std_error': std_error,
                'method': 'permuco'
            }
        else:
            return {'success': False, 'error': 'Failed to extract results'}
            
    except Exception as e:
        return {'success': False, 'error': str(e)}

def calculate_statistic(y, x, stat_type):
    """Calculate various statistics"""
    if stat_type == 't_stat':
        # Two-sample t-statistic
        unique_groups = np.unique(x)
        if len(unique_groups) == 2:
            group1 = y[x == unique_groups[0]]
            group2 = y[x == unique_groups[1]]
            return stats.ttest_ind(group1, group2)[0]
        else:
            return 0
    elif stat_type == 'mean_diff':
        # Mean difference
        unique_groups = np.unique(x)
        if len(unique_groups) == 2:
            group1 = y[x == unique_groups[0]]
            group2 = y[x == unique_groups[1]]
            return np.mean(group1) - np.mean(group2)
        else:
            return 0
    elif stat_type == 'correlation':
        # Correlation (for continuous variables)
        return stats.pearsonr(y, x)[0] if not np.isnan(x).any() else 0
    else:
        return 0

def calculate_full_model_statistic(y, x_full, stat_type):
    """Calculate statistics for full model with proper standard errors"""
    if stat_type == 't_stat' and x_full.shape[1] >= 2:
        # 正しいt統計量の計算
        reg = LinearRegression()
        reg.fit(x_full, y)
        
        # 予測値と残差
        y_pred = reg.predict(x_full)
        residuals = y - y_pred
        
        # 残差平方和と自由度
        n = len(y)
        p = x_full.shape[1]
        df_resid = n - p
        
        if df_resid <= 0:
            return 0
        
        # 残差標準誤差
        mse = np.sum(residuals**2) / df_resid
        
        # デザイン行列の逆行列の計算（数値安定性のためpinvを使用）
        try:
            xtx_inv = np.linalg.pinv(x_full.T @ x_full)
            # 最後の係数（興味のある変数）の標準誤差
            se_coef = np.sqrt(mse * xtx_inv[-1, -1])
            
            if se_coef == 0:
                return 0
            
            # t統計量 = 係数 / 標準誤差
            t_stat = reg.coef_[-1] / se_coef
            return t_stat
            
        except (np.linalg.LinAlgError, ValueError):
            # 数値的に不安定な場合は0を返す
            return 0
            
    elif stat_type == 'correlation' and x_full.shape[1] >= 2:
        try:
            corr, _ = stats.pearsonr(y, x_full[:, -1])
            return corr if not np.isnan(corr) else 0
        except:
            return 0
    else:
        return calculate_statistic(y, x_full[:, -1] if x_full.shape[1] > 0 else x_full, stat_type)

# UI部分
st.sidebar.title("Analysis Options")

# Analysis method selection
test_method = st.sidebar.radio(
    "Select test method:",
    ["Simple Permutation", "Freedman-Lane", "Freedman-Lane (permuco)"],
    index=2
)

# Permutation options
st.sidebar.markdown("### Permutation Options:")
n_permutations = st.sidebar.number_input(
    "Number of permutations", 
    min_value=100, max_value=10000, value=1000, step=100
)

# 再現性のためのシード設定
random_seed = st.sidebar.number_input(
    "Random seed (for reproducibility)",
    min_value=0, max_value=99999, value=42, step=1
)

stat_type = st.sidebar.radio(
    "Test statistic:",
    ["t_stat", "mean_diff", "correlation"],
    index=0
)

with st.sidebar.expander("ℹ️ About test methods"):
    st.markdown("""
    **Simple Permutation Test:**
    - ラベル（群分け）をランダムに置換
    - 群間の差が偶然によるものかを検定
    - 共変量を考慮しない
    
    **Freedman-Lane Method:**
    - 共変量（バッチ効果等）を保持したまま検定
    - Reduced model（共変量のみ）の残差を置換
    - より厳密な統計的推論が可能
    - サンプル数が少ない場合に特に有効
    
    **Freedman-Lane (permuco):**
    - 確立されたRパッケージpermucoを使用
    - 我々の実装との比較・検証用
    - 論文での標準的な参照実装
    - より多くの統計情報を提供
    
    **Test Statistics:**
    - t_stat: 正確なt統計量（係数/標準誤差）
    - mean_diff: 平均値差
    - correlation: 相関係数（連続変数用）
    
    **統計的保証:**
    - Type I エラー率の正確な制御
    - 有限サンプルでの正確性
    - 分布仮定に依存しない頑健性
    """)


# データアップロード部分（Calc_DESeq2.pyから流用）
df = None
use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Available data")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        file_name_head = st.session_state.uploaded_file_name

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "", ('auto', 'tsv','csv','excel'), index = 0, label_visibility = 'collapsed')
    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])

    if uploaded_file is not None:
        if file_type == 'auto':
            try:
                df = read_csv(uploaded_file, sep = None)
                df = df.set_index(df.columns[0])
            except:
                df = read_excel(uploaded_file)
                df = df.set_index(df.columns[0])
        elif file_type != 'excel':
            if file_type == 'csv':
                df = read_csv(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep = '\t')
            df = df.set_index(df.columns[0])
        else:
            df = read_excel(uploaded_file)
            df = df.set_index(df.columns[0])
        
        file_name_head = os.path.splitext(uploaded_file.name)[0]

if df is not None:
    st.write('Data shape:', df.shape)
    st.write(df.head())
    
    # データの前処理
    df = df.astype(float)
    df = rename_duplicates(df)
    
    # 基本統計
    st.write("### Basic Statistics")
    col1, col2 = st.columns(2)
    with col1:
        st.write(f"Genes: {df.shape[0]}")
        st.write(f"Samples: {df.shape[1]}")
    with col2:
        st.write(f"Missing values: {df.isnull().sum().sum()}")
        st.write(f"Zero values: {(df == 0).sum().sum()}")
    
    # グループとバッチの設定（Calc_DESeq2.pyから流用）
    condition = [str(i) for i in df.columns.tolist()]
    group_condition = remove_common_suffix(condition)
    group_condition = [remove_sample_num(x) for x in group_condition]
    group_condition = [str(x).rstrip('.') for x in group_condition]

    df_e = pd.DataFrame(group_condition, index=condition, columns=["Group"])
    df_b = pd.DataFrame({'Batch': [str(x) for x in group_condition]}, index=condition)

    st.markdown("### Experimental Design")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("#### Group Assignment:")
        edited_df_e = st.data_editor(df_e)
        condition_groups = edited_df_e.iloc[:,0].tolist()
        
    with col2:
        if test_method == "Simple Permutation":
            st.markdown("#### Simple Permutation Test")
            st.info("No covariates needed for simple permutation test")
            batch_groups = None
        else:
            st.markdown("#### Batch/Covariate Assignment:")
            edited_df_b = st.data_editor(df_b)
            batch_groups = [str(x) for x in edited_df_b['Batch'].tolist()]

    # 解析実行
    if st.button("🚀 Run Permutation Test"):
        # 再現性のためのシード設定
        np.random.seed(random_seed)
        
        # 群の数をチェック
        unique_groups = list(set(condition_groups))
        if len(unique_groups) < 2:
            st.error("At least 2 groups are required for analysis")
            sys.exit()
        
        st.write(f"Groups found: {unique_groups}")
        st.write(f"Group assignment: {dict(zip(condition, condition_groups))}")
        st.write(f"Using random seed: {random_seed} for reproducibility")
        
        # 各遺伝子について解析実行
        results_list = []
        
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # サンプル数制限（計算時間を考慮）
        n_genes = min(df.shape[0], 100)  # 最大100遺伝子まで
        if df.shape[0] > 100:
            st.warning(f"Analysis limited to first {n_genes} genes for computational efficiency")
        
        for i, gene in enumerate(df.index[:n_genes]):
            progress_bar.progress((i + 1) / n_genes)
            status_text.text(f'Processing gene {i+1}/{n_genes}: {gene}')
            
            y = df.loc[gene].values
            
            if test_method == "Simple Permutation":
                # 群ラベルを数値に変換
                x = np.array([unique_groups.index(g) for g in condition_groups])
                
                try:
                    obs_stat, perm_stats, p_val = simple_permutation_test(
                        y, x, n_permutations, stat_type
                    )
                    
                    results_list.append({
                        'Gene': gene,
                        'Observed_Statistic': obs_stat,
                        'P_value': p_val,
                        'Method': 'Simple_Permutation'
                    })
                except Exception as e:
                    st.warning(f"Error processing gene {gene}: {str(e)}")
                    continue
                    
            elif test_method == "Freedman-Lane":  # Our Freedman-Lane implementation
                # 共変量マトリックス作成
                unique_batches = list(set(batch_groups))
                
                # Full model: intercept + batch + group
                x_full = np.ones((len(condition_groups), 1))  # intercept
                
                # Batch variables (one-hot encoding)
                for batch in unique_batches[:-1]:  # 最後のバッチは参照群
                    batch_col = np.array([1 if b == batch else 0 for b in batch_groups]).reshape(-1, 1)
                    x_full = np.hstack([x_full, batch_col])
                
                # Group variable
                group_col = np.array([unique_groups.index(g) for g in condition_groups]).reshape(-1, 1)
                x_full = np.hstack([x_full, group_col])
                
                # Reduced model: intercept + batch only
                x_reduced = x_full[:, :-1]  # 群変数を除く
                
                try:
                    obs_stat, perm_stats, p_val = freedman_lane_test(
                        y, x_full, x_reduced, n_permutations, stat_type
                    )
                    
                    results_list.append({
                        'Gene': gene,
                        'Observed_Statistic': obs_stat,
                        'P_value': p_val,
                        'Method': 'Freedman_Lane'
                    })
                except Exception as e:
                    st.warning(f"Error processing gene {gene}: {str(e)}")
                    continue
                    
            else:  # Freedman-Lane (permuco)
                # permucoを使用
                group_numeric = np.array([unique_groups.index(g) for g in condition_groups])
                unique_batches = list(set(batch_groups))
                batch_numeric = np.array([unique_batches.index(b) for b in batch_groups])
                
                try:
                    result = permuco_freedman_lane_test(
                        y, group_numeric, batch_numeric, n_permutations, random_seed
                    )
                    
                    if result['success']:
                        results_list.append({
                            'Gene': gene,
                            'Observed_Statistic': result['t_statistic'],
                            'P_value': result['p_value'],
                            'Method': 'Freedman_Lane_permuco'
                        })
                    else:
                        st.warning(f"permuco failed for gene {gene}: {result.get('error', 'Unknown error')}")
                        continue
                        
                except Exception as e:
                    st.warning(f"Error processing gene {gene} with permuco: {str(e)}")
                    continue
        
        progress_bar.empty()
        status_text.empty()
        
        if results_list:
            # 結果をDataFrameに変換
            results_df = pd.DataFrame(results_list)
            
            # FDR補正
            from statsmodels.stats.multitest import fdrcorrection
            
            # デバッグ情報
            st.write(f"P-value range: {results_df['P_value'].min():.6f} - {results_df['P_value'].max():.6f}")
            st.write(f"Unique P-values: {results_df['P_value'].nunique()}")
            
            # FDR補正を実行
            rejected, fdr_corrected = fdrcorrection(results_df['P_value'], method='indep')
            results_df['FDR'] = fdr_corrected
            
            # FDR補正後の確認
            st.write(f"FDR range: {results_df['FDR'].min():.6f} - {results_df['FDR'].max():.6f}")
            st.write(f"Unique FDR values: {results_df['FDR'].nunique()}")
            
            # 結果表示
            st.markdown("### Results")
            st.write(f"Analyzed {len(results_df)} genes")
            
            # 統計サマリー
            col1, col2, col3 = st.columns(3)
            with col1:
                sig_genes = (results_df['P_value'] < 0.05).sum()
                st.metric("Significant genes (p < 0.05)", sig_genes)
            with col2:
                fdr_sig_genes = (results_df['FDR'] < 0.05).sum()
                st.metric("FDR significant genes", fdr_sig_genes)
            with col3:
                median_p = results_df['P_value'].median()
                st.metric("Median p-value", f"{median_p:.4f}")
            
            # 結果テーブル
            st.dataframe(
                results_df.sort_values('P_value').head(20),
                use_container_width=True
            )
            
            # 可視化
            st.markdown("### Visualization")
            
            # P-value histogram
            col1, col2 = st.columns(2)
            
            with col1:
                fig, ax = plt.subplots(figsize=(8, 6))
                ax.hist(results_df['P_value'], bins=20, edgecolor='black', alpha=0.7)
                ax.set_xlabel('P-value')
                ax.set_ylabel('Frequency')
                ax.set_title('P-value Distribution')
                ax.axvline(x=0.05, color='red', linestyle='--', label='p = 0.05')
                ax.legend()
                st.pyplot(fig)
            
            with col2:
                # QQ plot
                fig, ax = plt.subplots(figsize=(8, 6))
                stats.probplot(results_df['P_value'], dist="uniform", plot=ax)
                ax.set_title('Q-Q Plot (Uniform Distribution)')
                ax.set_xlabel('Theoretical Quantiles')
                ax.set_ylabel('Sample Quantiles')
                st.pyplot(fig)
            
            # Volcano plot equivalent (if we have effect sizes)
            if stat_type in ['mean_diff', 't_stat']:
                fig, ax = plt.subplots(figsize=(10, 6))
                ax.scatter(results_df['Observed_Statistic'], -np.log10(results_df['P_value']), 
                          alpha=0.6)
                ax.set_xlabel(f'Observed {stat_type}')
                ax.set_ylabel('-log10(P-value)')
                ax.set_title('Effect Size vs Significance')
                ax.axhline(y=-np.log10(0.05), color='red', linestyle='--', label='p = 0.05')
                ax.legend()
                st.pyplot(fig)
            
            # ダウンロード
            csv = results_df.to_csv(index=False)
            st.download_button(
                label="📥 Download Results",
                data=csv,
                file_name=f"{file_name_head}_{test_method}_results.csv",
                mime="text/csv"
            )
            
        else:
            st.error("No results generated. Please check your data and parameters.")