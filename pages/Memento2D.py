
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
import itertools
import re

from statsmodels.stats.multitest import multipletests
st.set_page_config(page_title="Memento", page_icon="💬")

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
                # マッピングを表示
                st.write(f"Mapping for {col}:")
                for i, cat in enumerate(categories):
                    st.write(f"{cat} -> {i}")
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

def run_memento_by_group(adata, groupby,  num_cpus=12, num_boot=5000, min_perc_group=0.7):
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
                memento.setup_memento(adata, q_column='capture_rate')
                memento.create_groups(adata, label_columns=label_columns) 
                # 3. モーメント計算
                memento.compute_1d_moments(adata, min_perc_group=min_perc_group)           
                # 4. メタデータ準備
                sample_meta = memento.get_groups(adata)
                # The covariate DataFrame - pick the covariate columns
                cov_df = sample_meta[cov_column]
                st.write()

                # The treatment DataFrame - pick the treatment column
                treat_df = sample_meta[['condition_encoded']]

                memento.ht_1d_moments(
                    adata,
                    treatment=sample_meta,
                    covariate=cov_df,
                    num_boot=num_boot,
                    num_cpus=num_cpus
                )
                result = memento.get_1d_ht_result(adata)   
            
            else:
                #wrapper functionをそのまま書く
                adata = adata.copy().copy()
             #   adata.obs['capture_rate'] = capture_rate
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
        filename = memento_temp_dir + "/memento/" + celltype.replace('/', '_').replace(' ', '_') + comparison_str + '_2D.tsv'
        result.to_csv(filename, sep='\t', index=True)
        st.write(f"Saved results for {celltype}")
        st.write(result.head())

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
def show_gene_counts_by_min_perc(_adata, min_perc_values=[0.5, 0.7, 0.9]):
   """
   異なるmin_perc_groupの値での残る遺伝子数を表示
   """
   # 元のデータの遺伝子数
   total_genes = adata.shape[1]
   st.write(f"Total genes before filtering: {total_genes}")
   
   for min_perc in min_perc_values:
       # テスト用のadataをコピー
       test_adata = adata.copy()
       
       # mementoのセットアップ
       memento.setup_memento(test_adata, q_column='capture_rate')
       memento.create_groups(test_adata, label_columns=['condition_encoded'])
       
       # モーメント計算
       memento.compute_1d_moments(test_adata, min_perc_group=min_perc)
       
       # 残った遺伝子数を取得
       remaining_genes = len(test_adata.uns['memento']['test_genes'])
       percent_remaining = (remaining_genes / total_genes) * 100
       
       st.write(f"min_perc_group = {min_perc}:")
       st.write(f"  Remaining genes: {remaining_genes} ({percent_remaining:.1f}%)")


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
 #   if not os.path.exists(memento_temp_dir + "/memento"):
  #      os.mkdir(memento_temp_dir + "/memento")
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
    if len(sample_list) < 2:
        st.warning("At least 2 conditions are required.")
        st.stop()

    groupby = None
    bycell = st.checkbox("Split data by cell-type?")
    if bycell:
        st.markdown("##### Select cell type:")
        groupby = st.selectbox("cell identity:", [x for x in meta if x != sampleby], index = find_first_index_or_zero(meta, ["cell.ident", "seurat_clusters"]), label_visibility = 'collapsed')
        cell_list = sorted(adata.obs[groupby].cat.categories.to_list())

    cov_column = []
    cov_df = None
    label_columns = ['condition_encoded']
    add_covariates = st.checkbox("Add covariates (e.g., replicates, batch)")
    if add_covariates:
        cov_options = [x for s in adata.obs.columns.to_list() if x not in [sampleby, groupby]]
        cov_column = st.multiselect("Choose covariates", cov_options)
        #cov_df = prepare_covariates(adata, cov_column)
        
        label_columns.extend(cov_column)

        # カテゴリカルをcovariatesように変換
        adata = convert_covariates_to_numeric(adata, covariate_cols=cov_column)
 
    adata_genes = adata.var.index.tolist()

    st.markdown("##### Input genes of interest (comma, space, CR separated):")
    genes = st.text_input("genes",label_visibility = 'collapsed')
    gene_list = []
    if len(genes) > 0:
        raw_genes = re.split(r'[,\s]+', genes)
        # Remove spaces from each gene name and filter out empty strings
        gene_list = [re.sub(r'\s', '', gene) for gene in raw_genes if gene.strip()]

        genes = list(filter(lambda x: x != "", genes)) #空白を除く
        gene_list =sorted(set(gene_list), key=gene_list.index)
        gene_list = [x for x in gene_list if x in adata_genes]
        st.write(gene_list[:3])

    Second_genes = st.checkbox("Specify the genes for the 2nd group?")
    st.write("If diselected, all genes are analyzed.")
    if Second_genes:
        st.markdown("##### Input genes of interest (comma, space, CR separated):")
        genes2 = st.text_input("genes2",label_visibility = 'collapsed')
        gene_list2 = []
        if len(genes2) > 0:
            raw_genes2 = re.split(r'[,\s]+', genes2)
            # Remove spaces from each gene name and filter out empty strings
            gene_list2 = [re.sub(r'\s', '', gene) for gene in raw_genes2 if gene.strip()]

            genes2 = list(filter(lambda x: x != "", genes2)) #空白を除く
            gene_list2 =sorted(set(gene_list2), key=gene_list2.index)
            gene_list2 = [x for x in gene_list2 if x in adata_genes]
            st.write(gene_list2[:3])
    if len(gene_list) > 0:
        if Second_genes:
            gene_pairs = list(itertools.product(gene_list, gene_list2))
        else:
            gene_pairs = list(itertools.product(gene_list, adata_genes))
    else:
        st.stop()

    st.write(gene_pairs[:5])

    use_saturation_col = st.checkbox("Use a column containing sequence saturation for each sample?")
    if use_saturation_col:
        saturation_col = st.selcetbox("Select thec colum for seq saturation", meta)
        adata.obs['capture_rate'] = adata.obs[[saturation_col]] * 0.25

    else:
        st.markdown("Please see web_summay.html for sequence saturation")
        saturation = st.number_input("Sequence saturation", min_value=0.10, max_value=1.00, value = 1.00, step = 0.05)
        capture_rate=0.25 * saturation
        adata.obs['capture_rate'] = capture_rate

    min_perc_group = st.number_input("min_percent_group", min_value=0.01, max_value=1.00, value = 0.70, step = 0.05)
    st.write("The minimum fraction of cells in each group where a gene must be detected to be included in the analysis")
    st.write("The default value is 0.7. (May need to decrease to analyze low expression genes.)")

    if st.button('Run analysis'):
        comparison_str = "_" + test_group + "_vs_" + control_group + "_"
        adata.obs['condition_encoded'] = (adata.obs[sampleby] == test_group).astype(int)
        adata.X = adata.raw.X
        # adata.Xの形式を確認し、CSRでない場合のみ変換
        if not isspmatrix_csr(adata.X):
            adata.X = csr_matrix(adata.X)


        num_cpus=12
        num_boot=5000


        # 実行例（例えばcell.ident2でグループ分け）
        if not bycell:
            memento.setup_memento(adata, q_column='capture_rate')
            memento.create_groups(adata, label_columns=label_columns)
            memento.compute_1d_moments(adata, min_perc_group=min_perc_group)
            filtered_gene_pairs = [(a,b) for a,b in gene_pairs if a in adata.var.index and b in adata.var.index]
            memento.compute_2d_moments(adata, filtered_gene_pairs)
            sample_meta = memento.get_groups(adata)
            if add_covariates:
                cov_df = sample_meta[cov_column]
            else:
                cov_df = None
            treat_df = sample_meta[['condition_encoded']]
            memento.ht_2d_moments(
                adata, 
                treatment=treat_df,
                covariate=cov_df,
                num_boot=5000, 
                verbose=1,
                num_cpus=12)

            results_dict = {}

            results_dict["all_cells"] = memento.get_2d_ht_result(adata)
            results_dict["all_cells"] = results_dict["all_cells"].sort_values('corr_pval')
            st.write(results_dict["all_cells"].head())

        else:
            results_dict = {}
        
            # グループごとにループ
            for group in cell_list:
                # グループでデータをサブセット
                group_adata = adata[adata.obs[groupby] == group].copy()
                
                try:
                    st.write(f"Analyzing {group}")


                    memento.setup_memento(group_adata, q_column='capture_rate')
                    memento.create_groups(group_adata, label_columns=label_columns)
                    memento.compute_1d_moments(group_adata, min_perc_group=min_perc_group)
                    filtered_gene_pairs = [(a,b) for a,b in gene_pairs if a in group_adata.var.index and b in group_adata.var.index]
                    memento.compute_2d_moments(group_adata, filtered_gene_pairs)
                    sample_meta = memento.get_groups(group_adata)
                    treat_df = sample_meta[['condition_encoded']]
                    if add_covariates:
                        cov_df = sample_meta[cov_column]
                    else:
                        cov_df = None
                    memento.ht_2d_moments(
                        group_adata, 
                        treatment=treat_df,
                        covariate=cov_df,
                        num_boot=5000, 
                        verbose=1,
                        num_cpus=12)

                    results_dict = {}

                    results_dict[group]= memento.get_2d_ht_result(adata)

                    # 結果を辞書に保存
                    results_dict[group] = result
                    st.write(f"Completed analysis for {group}")
                    st.write(results_dict[group].head())
                    
                except Exception as e:
                    st.write(f"Error in group {group}: {str(e)}")
                    continue
            
        for celltype, result in results_dict.items():
            # DM genes の補正
            result['corr_padj'] = multipletests(result.corr_pval, method='fdr_bh')[1]
            result.sort_values('corr_pval', ascending=True, inplace=True)

        st.markdown("""
### 

""")


        output_dir = os.path.join(memento_temp_dir, "memento")
        os.makedirs(output_dir, exist_ok=True)

        st.write(groupby)
        
        # 結果をTSVとして保存
        save_results_to_tsv(
            results_dict=results_dict,
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
                


