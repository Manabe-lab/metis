
import scanpy as sc
import memento
import pandas as pd
import numpy as np
import os
from helper_func import clear_old_directories
from helper_func import clear_old_files
import time
import streamlit as st
from scipy.sparse import csr_matrix, isspmatrix_csr, issparse
import shutil
from statsmodels.stats.multitest import multipletests
st.set_page_config(page_title="Memento", page_icon="💬")

import scipy as sp
from scipy import interpolate

def estimate(pv, m=None, verbose=False, lowmem=False, pi0=None):
    """
    Estimates q-values from p-values

    Args
    =====

    m: number of tests. If not specified m = pv.size
    verbose: print verbose messages? (default False)
    lowmem: use memory-efficient in-place algorithm
    pi0: if None, it's estimated as suggested in Storey and Tibshirani, 2003.
         For most GWAS this is not necessary, since pi0 is extremely likely to be
         1

    """
    assert(pv.min() >= 0 and pv.max() <= 1), "p-values should be between 0 and 1"

    original_shape = pv.shape
    pv = pv.ravel()  # flattens the array in place, more efficient than flatten()

    if m is None:
        m = float(len(pv))
    else:
        # the user has supplied an m
        m *= 1.0

    # if the number of hypotheses is small, just set pi0 to 1
    if len(pv) < 100 and pi0 is None:
        pi0 = 1.0
    elif pi0 is not None:
        pi0 = pi0
    else:
        # evaluate pi0 for different lambdas
        pi0 = []
        lam = np.arange(0, 0.90, 0.01)
        counts = np.array([(pv > i).sum() for i in np.arange(0, 0.9, 0.01)])
        for l in range(len(lam)):
            pi0.append(counts[l]/(m*(1-lam[l])))

        pi0 = np.array(pi0)

        # fit natural cubic spline
        tck = interpolate.splrep(lam, pi0, k=3)
        pi0 = interpolate.splev(lam[-1], tck)
        if verbose:
            print("qvalues pi0=%.3f, estimated proportion of null features " % pi0)

        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f) while estimating qvalues, setting it to 1" % pi0)
            pi0 = 1.0

    assert(pi0 >= 0 and pi0 <= 1), "pi0 is not between 0 and 1: %f" % pi0

    if lowmem:
        # low memory version, only uses 1 pv and 1 qv matrices
        qv = sp.zeros((len(pv),))
        last_pv = pv.argmax()
        qv[last_pv] = (pi0*pv[last_pv]*m)/float(m)
        pv[last_pv] = -sp.inf
        prev_qv = last_pv
        for i in range(int(len(pv))-2, -1, -1):
            cur_max = pv.argmax()
            qv_i = (pi0*m*pv[cur_max]/float(i+1))
            pv[cur_max] = -sp.inf
            qv_i1 = prev_qv
            qv[cur_max] = min(qv_i, qv_i1)
            prev_qv = qv[cur_max]

    else:
        p_ordered = np.argsort(pv)
        pv = pv[p_ordered]
        qv = pi0 * m/len(pv) * pv
        qv[-1] = min(qv[-1], 1.0)

        for i in range(len(pv)-2, -1, -1):
            qv[i] = min(pi0*m*pv[i]/(i+1.0), qv[i+1])

        # reorder qvalues
        qv_temp = qv.copy()
        qv = np.zeros_like(qv)
        qv[p_ordered] = qv_temp

    # reshape qvalues
    qv = qv.reshape(original_shape)

    return qv

def convert_covariates_to_numeric(adata, covariate_cols):
    """
    指定されたカラムを数値に変換します。
    
    Args:
        adata: AnnData object
        covariate_cols: 数値に変換したいカラムのリスト
    
    Returns:
        変換後のAnnData object
    """
    adata = adata.copy()
    
    for col in covariate_cols:
        # カラムの型を確認
        if adata.obs[col].dtype == 'object' or adata.obs[col].dtype == 'category':
            unique_values = adata.obs[col].unique()
            if len(unique_values) <= 100:
                # カテゴリ数が10以下の場合はラベルエンコーディング
                categories = pd.Categorical(adata.obs[col]).categories
                adata.obs[col] = pd.Categorical(adata.obs[col]).codes

            else:
                try:
                    # 10分位数に分割（重複を許可）
                    adata.obs[col] = pd.qcut(
                        pd.factorize(adata.obs[col])[0], 
                        q=10, 
                        labels=False,
                        duplicates='drop'  # 重複する値は削除
                    )
                    st.write(f"{col} was divided into quantiles")
                except:
                    # qcutが失敗した場合は単純なラベルエンコーディング
                    categories = pd.Categorical(adata.obs[col]).categories
                    adata.obs[col] = pd.Categorical(adata.obs[col]).codes
                    st.write(f"Fallback to label encoding for {col}")
                    for i, cat in enumerate(categories):
                        st.write(f"{cat} -> {i}")
        else:
            st.write(f"{col} is already numeric")
            
    return adata

def run_memento_by_group(adata, groupby, label_columns, num_cpus=12, num_boot=5000, min_perc_group=0.9):
    # グループごとの結果を格納する辞書
    results_dict = {}
    
    # グループごとにループ
    for group in adata.obs[groupby].unique():
        # グループでデータをサブセット
        group_adata = adata[adata.obs[groupby] == group].copy()
        
        try:
            # mementoの解析を実行

            if len(cov_column) > 0:

            #    adata.obs['capture_rate'] = capture_rate
                memento.setup_memento(group_adata, q_column='capture_rate')
                memento.create_groups(group_adata, label_columns=label_columns) 
                # 3. モーメント計算
                memento.compute_1d_moments(group_adata, min_perc_group=min_perc_group)           
                # 4. メタデータ準備
                sample_meta = memento.get_groups(group_adata)
                # The covariate DataFrame - pick the covariate columns
                cov_df = sample_meta[cov_column]
                st.write()

                # The treatment DataFrame - pick the treatment column
                treat_df = sample_meta[['condition_encoded']]

                memento.ht_1d_moments(
                    group_adata,
                    treatment=treat_df,
                    covariate=cov_df,
                    num_boot=num_boot,
                    num_cpus=num_cpus
                )
                result = memento.get_1d_ht_result(group_adata)   
            
            else:
                #wrapper functionをそのまま書く
             #   adata.obs['capture_rate'] = capture_rate
                memento.setup_memento(group_adata, q_column='capture_rate')
                memento.create_groups(group_adata, label_columns=['condition_encoded'])
                memento.compute_1d_moments(group_adata, min_perc_group=min_perc_group)
                sample_meta = memento.get_groups(group_adata)[['condition_encoded']]
                memento.ht_1d_moments(
                    group_adata, 
                    treatment=sample_meta,
                    num_boot=num_boot,
                    num_cpus=num_cpus)
                result = memento.get_1d_ht_result(group_adata)

            # 結果を辞書に保存
            results_dict[group] = result
            st.write(f"Completed analysis for {group}")
            
        except Exception as e:
            st.write(f"Error in group {group}: {str(e)}")
            continue
    
    return results_dict



def save_results_to_tsv(results_dict, memento_temp_dir, comparison_str, control_group, test_group, adata, groupby):
    if os.path.exists(memento_temp_dir + "/memento"):
        shutil.rmtree(memento_temp_dir + "/memento")
    os.mkdir(memento_temp_dir + "/memento")
    for celltype, result in results_dict.items():
        filename = memento_temp_dir + "/memento/" + celltype.replace('/', '_').replace(' ', '_') + comparison_str + '.tsv'
        result.set_index('gene', inplace=True)

        try:
            # その細胞タイプのデータを抽出
            if groupby == "asone":
                adata_subset = adata.copy()
            else:
                adata_subset = adata[adata.obs[groupby] == celltype].copy()
            
            # データの状態を確認
          #  st.write(f"\nDiagnostics for {celltype}:")
         #   value_counts = adata_subset.obs['condition_encoded'].value_counts()
          #  st.write("condition_encoded values:")
         #   st.write(value_counts)
            
            if len(adata_subset) == 0:
                st.write(f"No cells found for {celltype}")
                continue
            
            # フィルタリングと正規化
            cell_counts = adata_subset.X.sum(axis=1).A1
            adata_subset = adata_subset[cell_counts > 0].copy()
            
            if len(adata_subset) == 0:
                st.write(f"No valid cells after filtering for {celltype}")
                continue
                
            # 正規化
            adata_norm = adata_subset.copy()
            sc.pp.normalize_total(adata_norm)
            sc.pp.log1p(adata_norm)
            
            # 条件ごとの平均を計算
            cond_0 = adata_norm.obs['condition_encoded'] == 0
            cond_1 = adata_norm.obs['condition_encoded'] == 1
            
            if not any(cond_0) or not any(cond_1):
                st.write(f"Missing one or both conditions for {celltype}")
                st.write("Number of condition 0:", sum(cond_0))
                st.write("Number of condition 1:", sum(cond_1))
                continue
                
            mean_0 = adata_norm[cond_0].X.mean(axis=0).A1
            mean_1 = adata_norm[cond_1].X.mean(axis=0).A1
            
            result_genes = result.index.to_list()
          
            df_mean = pd.DataFrame(index = adata_norm.var_names)
            df_mean[f'mean_expr_{control_group}'] = mean_0
            df_mean[f'mean_expr_{test_group}'] = mean_1
            df_mean['log2FC'] = df_mean[f'mean_expr_{test_group}'] - df_mean[f'mean_expr_{control_group}']
            result['log2FC'] = df_mean.loc[result_genes,'log2FC']
            result[f'mean_expr_{control_group}'] = df_mean.loc[result_genes,f'mean_expr_{control_group}']
            result[f'mean_expr_{test_group}'] = df_mean.loc[result_genes,f'mean_expr_{test_group}']
            result = result.sort_values(by='de_pval', ascending = True)
            st.write(f"Saved results for {celltype}")
            st.write(result.head())
            # 保存
            result.to_csv(filename, sep='\t')
            
        except Exception as e:
            st.write(f"Error saving {celltype}: Detailed error - {type(e).__name__}: {str(e)}")


def prepare_covariates(adata, covariate_cols):
    """
    adata.obsから指定された列を取得し、必要に応じてダミー変数化します。
    
    Args:
        adata: AnnData object
        covariate_cols: covariateとして使用する列名のリスト
    
    Returns:
        pd.DataFrame: 処理済みのcovariate dataframe
    """
    cov_df = pd.DataFrame(index=adata.obs.index)
    
    for col in covariate_cols:
        # 列のデータ型を確認
        dtype = adata.obs[col].dtype
        
        # カテゴリカル変数（objectまたはcategory）の場合
        if dtype == 'object' or dtype == 'category':
            # ユニークな値の数を確認
            n_unique = len(adata.obs[col].unique())
            
            if n_unique == 2:
                # 2値カテゴリの場合は0/1に変換
                cov_df[col] = (adata.obs[col] == adata.obs[col].unique()[1]).astype(int)
            else:
                # 3値以上のカテゴリはダミー変数化
                dummy = pd.get_dummies(adata.obs[col], prefix=col)
                # 多重共線性を避けるため、最後の列を除外
                cov_df = pd.concat([cov_df, dummy.iloc[:, :-1]], axis=1)
        
        # 数値型の場合はそのまま使用
        else:
            cov_df[col] = adata.obs[col]
    
    # interceptを追加
    cov_df['intercept'] = 1
    
    return cov_df


@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    st.write("uploaded_file data")
    temp_df = pd.DataFrame(
    adata.raw.X[:5,:8].toarray() if issparse(adata.raw.X) else adata.raw.X[:5,:8],
    index=adata.obs_names[:5],
    columns=adata.var_names[:8]
    )
    st.dataframe(temp_df) 
    return adata

def find_first_index_or_zero(lst, elements):
    for element in elements:
        try:
            return lst.index(element)
        except ValueError:
            continue
    return 0



#============ 新しいファイルをアップロードしたときは、cacheをclearする

def get_file_identifier(file):
    if file is not None:
        return f"{file.name}_{file.size}"
    return None


# temp内に保存する
# --- Initialising SessionState ---
if "memento_temp_dir" not in st.session_state:
    memento_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
    if not os.path.exists(memento_temp_dir):
        os.mkdir(memento_temp_dir)
#    if not os.path.exists(memento_temp_dir + "/memento"):
#        os.mkdir(memento_temp_dir + "/memento")
    st.session_state.memento_temp_dir = memento_temp_dir
else:
    memento_temp_dir = st.session_state.memento_temp_dir


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


    meta = adata.obs.columns.to_list()
    for i in ['nFeature_RNA','nCount_RNA','percent.mt', 'Cell_id']:
        try:
            meta.remove(i)
        except:
            pass
    submitted_basic = False

    st.markdown("##### Select condition:")
    sampleby = st.selectbox("sample/condition:", meta, index = find_first_index_or_zero(meta, ["orig.ident",
    "sample","KO"]), label_visibility = 'collapsed')
    sample_list = sorted(adata.obs[sampleby].cat.categories.to_list())

    if len(sample_list) > 2:
        # 複数の条件から2つを選択するマルチセレクト
        selected_conditions = st.multiselect(
            "Select 2 conditions to compare:",
            sample_list,
            max_selections=2  # 最大2つまで選択可能
        )
        
        # 2つ選択されるまで待機
        if len(selected_conditions) != 2:
            st.warning("Please select exactly 2 conditions to compare.")
            st.stop()
            
        # 選択された2条件でデータをサブセット化
        mask = adata.obs[sampleby].isin(selected_conditions)
        adata = adata[mask].copy()
        
        # condition_encodedカラムの作成（1つ目の条件を0、2つ目を1とする）
        adata.obs['condition_encoded'] = (adata.obs[sampleby] == selected_conditions[1]).astype(int)
        
        # 選択された条件を表示
        st.write(f"Comparing: {selected_conditions[0]} and {selected_conditions[1]}")

        sample_list = selected_conditions
        
    elif len(sample_list) < 2:
        st.warning("At least 2 conditions are required.")
        st.stop()

    st.markdown("##### Select control group:")
    control_group = st.radio("Control group",[sample_list[0], sample_list[1]], index =0, label_visibility = 'collapsed')
    test_group = [i for i in sample_list if i != control_group][0]

    comparison_str = "_" + test_group + "_vs_" + control_group + "_"

    st.markdown("---")

    asone = st.checkbox("Treat all cells as one type?")

    if not asone:
        st.markdown("##### Select cell type:")
        groupby = st.selectbox("cell identity:", [x for x in meta if x != sampleby], index = find_first_index_or_zero(meta, ["cell.ident", "seurat_clusters"]), label_visibility = 'collapsed')
        cell_list = sorted(adata.obs[groupby].cat.categories.to_list())
    else:
        groupby = "asone"

    st.write(f"Conditions: {', '.join(sample_list)}")
    if not asone:
        st.write(f"Cell types: {', '.join(cell_list)}")

    st.markdown("---")

    cov_column = []
    label_columns = ['condition_encoded']
    add_covariates = st.checkbox("Add covariates (e.g., replicates, batch)")
    if add_covariates:
        cov_options = [x for x in adata.obs.columns.to_list() if x not in [sampleby, groupby]]
        cov_column = st.multiselect("Choose covariates", cov_options)
        #cov_df = prepare_covariates(adata, cov_column)        
        label_columns.extend(cov_column)

        # カテゴリカルをcovariatesように変換
        adata = convert_covariates_to_numeric(adata, covariate_cols=cov_column)
    
    use_saturation_col = st.checkbox("Use a column containing sequence saturation for each sample?")
    if use_saturation_col:
        saturation_col = st.selectbox("Select thec colum for seq saturation", meta)
        adata.obs['capture_rate'] = adata.obs[[saturation_col]] * 0.25

    else:
        st.markdown("Please see web_summay.html for sequence saturation")
        saturation = st.number_input("Sequence saturation", min_value=0.10, max_value=1.00, value = 1.00, step = 0.05)
        capture_rate=0.25 * saturation
        adata.obs['capture_rate'] = capture_rate

    min_perc_group = st.number_input("min_percent_group", min_value=0.01, max_value=1.00, value = 0.70, step = 0.05)
    st.write("The minimum fraction of cells in each group where a gene must be detected to be included in the analysis")
    st.write("The default value is 0.7.(May need to decrease to analyze low expression genes.)")
    st.write("sequence saturationを低下させるほど、またmin_percent_groupを低下させるほど、単純な比較よりFCが乖離する遺伝子が増える傾向")
    st.write("min_percent_group = 0.5までは大きな乖離はなさそうだが、データセットによって確認を")

    if st.button('Run analysis'):
        adata.obs['condition_encoded'] = (adata.obs[sampleby] == test_group).astype(int) #test_groupが1になる
        adata.X = adata.raw.X
        # adata.Xの形式を確認し、CSRでない場合のみ変換
        if not isspmatrix_csr(adata.X):
            adata.X = csr_matrix(adata.X)


        num_cpus=12
        num_boot=5000


        # 実行例（例えばcell.ident2でグループ分け）
        if not asone:
            results_by_celltype = run_memento_by_group(adata, groupby=groupby,  min_perc_group=min_perc_group, label_columns=label_columns)
        else:
            results_by_celltype = {}
                # mementoの解析を実行       
            if len(cov_column) > 0:

#                adata.obs['capture_rate'] = capture_rate
                memento.setup_memento(adata, q_column='capture_rate')
                memento.create_groups(adata, label_columns=label_columns) 
                # 3. モーメント計算
                memento.compute_1d_moments(adata, min_perc_group=min_perc_group)           
                # 4. メタデータ準備
                sample_meta = memento.get_groups(adata)
                # The covariate DataFrame - pick the covariate columns
                cov_df = sample_meta[cov_column]

                # The treatment DataFrame - pick the treatment column
                treat_df = sample_meta[['condition_encoded']]


                memento.ht_1d_moments(
                    adata,
                    treatment=treat_df,
                    covariate=cov_df,
                    num_boot=num_boot,
                    num_cpus=num_cpus
                )
                result = memento.get_1d_ht_result(adata)   
            
            else:
                #wrapper functionをそのまま書く
                adata = adata.copy().copy()
                adata.obs['capture_rate'] = capture_rate
                memento.setup_memento(adata, q_column='capture_rate')
                memento.create_groups(adata, label_columns=['condition_encoded'])
                memento.compute_1d_moments(adata, min_perc_group=min_perc_group)
                sample_meta = memento.get_groups(adata)[['condition_encoded']]
                memento.ht_1d_moments(
                    adata, 
                    treatment=sample_meta,
                    num_boot=num_boot,
                    num_cpus=num_cpus)
                result = memento.get_1d_ht_result(adata)


            results_by_celltype["all_cells"] = result
  #      st.write(results_by_celltype)

        for celltype, result in results_by_celltype.items():
            # DM genes の補正
            result['de_padj'] = multipletests(result.de_pval, method='fdr_bh')[1]
            result['de_STqval'] = estimate(result.de_pval)
            # DV genes の補正
            result['dv_padj'] = multipletests(result.dv_pval, method='fdr_bh')[1]

                        # DE p-valueでソート
            result.sort_values('de_pval', ascending=True, inplace=True)
        st.markdown("##### Results summary")
        # 結果の確認例
        for celltype, result in results_by_celltype.items():
            st.write(f"\nResults for {celltype}:")
            st.write(f"Number of significant DM genes (FDR q<0.01): {sum(result.de_padj < 0.01)}")
            st.write(f"Number of significant DV genes (FDR q<0.1): {sum(result.dv_padj < 0.1)}")

        st.markdown("""
### 
##### Differential Mean (DM) genes：
de_coef (mean effect size): 2群間での遺伝子発現の平均値の差の大きさを示します

単純な2群間比較ではln(FC) covariatesがある場合は単純なln(FC)ではない

正の値：処理群で発現が上昇

負の値：処理群で発現が減少

de_pval: 平均値の差に関する統計的有意性を示すp値（補正前）

de_padj: FDR

##### Differential Variability (DV) genes：
ve_coef (variability effect size): 2群間での遺伝子発現の変動（ばらつき）の差の大きさ

正の値：処理群で発現のばらつきが増加

負の値：処理群で発現のばらつきが減少

dv_pval: 発現変動の差に関する統計的有意性を示すp値（補正前）)

dv_padj: FDR

原著ではDMG (differentially mean expression genes)としてFDR<0.01を用いている。
DVG (Differentially variable genes) FDR<0.1
""")


        output_dir = os.path.join(memento_temp_dir, "memento")
        os.makedirs(output_dir, exist_ok=True)

        st.write(groupby)
        
        # 結果をTSVとして保存
        save_results_to_tsv(
            results_dict=results_by_celltype,
            memento_temp_dir=memento_temp_dir,
            comparison_str=comparison_str,
            control_group=control_group,
            test_group=test_group,
            adata=adata,
            groupby=groupby
        )
        # ZIPファイルのパス
        zip_path = os.path.join(memento_temp_dir, "memento.zip")
        

        # ZIPファイルの作成
        shutil.make_archive(
            os.path.join(memento_temp_dir, "memento"), 
            'zip',
            root_dir=output_dir
        )
        
        # ZIPファイルの読み込み
        with open(zip_path, "rb") as fp:
            zip_data = fp.read()
        
        # ダウンロードボタンの表示
        st.download_button(
            label="Download Results",
            data=zip_data,
            file_name="memento_results.zip",
            mime="application/zip"
        )
                
