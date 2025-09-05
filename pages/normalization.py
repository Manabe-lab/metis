import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
import time
from helper_func import clear_old_directories
from helper_func import clear_old_files
from io import StringIO


st.set_page_config(page_title="Count_normalization", page_icon="√")

def normalize_totalreads(df):
    return 10**6 * df / df.sum()

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

def remove_sample_num(i):
    i = re.sub(r'[_-][^-_]*\d+$', '', i)
    i = re.sub(r'\d+$', '', i)
    return i

def remove_after_space(i):
    m = re.match(r'([^\ ]+)(\ )+.+',i)
    if m is not None:
        return m.group(1)
    else:
        return i

def plot_mean_sd_comparison(df_before, df_after):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 各遺伝子（行）の平均と標準偏差を計算
    mean_before = df_before.mean(axis=1)
    sd_before = df_before.std(axis=1)
    mean_after = df_after.mean(axis=1)
    sd_after = df_after.std(axis=1)
    
    # 正規化前のプロット
    ax1.scatter(mean_before, sd_before, alpha=0.3, color='blue')
    ax1.set_title('Before Normalization', fontsize = 30)
    ax1.set_xlabel('Mean (log scale)', fontsize = 24)
    ax1.set_ylabel('Standard Deviation (log scale)', fontsize = 24)
    ax1.set_xscale('symlog')
    ax1.set_yscale('symlog')
    ax1.grid(True, linestyle='--', alpha=0.7)
    
    # 正規化後のプロット
    ax2.scatter(mean_after, sd_after, alpha=0.3, color='red')
    ax2.set_title('After Normalization', fontsize = 30)
    ax2.set_xlabel('Mean (log scale)', fontsize = 24)
    ax2.set_ylabel('Standard Deviation (log scale)', fontsize = 24)
    ax2.set_xscale('symlog')
    ax2.set_yscale('symlog')
    ax2.grid(True, linestyle='--', alpha=0.7)
    
    # 全体のタイトル
    fig.suptitle('Mean-SD Comparison Before and After Normalization', fontsize=36)
    
    # レイアウトの調整
    plt.tight_layout()
    
    return fig


# Rを呼び出すときはsession_stateを設定する
if 'use_R' not in st.session_state:
    st.session_state.use_R = False


if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #   st.write('file name kept')
    filename_add = ""



input_file_type = st.radio(
    "Data format:",
    ('auto','tsv','csv', 'excel', 'Homer'), key='auto')

input_file_format = 'row = gene'

if input_file_type != 'Homer':
    input_file_format = st.radio(
        "Data structure:",
        ('row = gene', 'column = gene'))

    if input_file_format == 'row = gene':
        st.markdown("""
    row = gene
    |  | Sample1 | Sample2 |
    | --- | --- | --- |
    | Gene1 |  |  |
    | Gene2 | ||

    """)
    else:
        st.markdown("""
    column = gene
    |  | Gene1 | Gene2 |
    | --- | --- | --- |
    | Sample1 |  |  |
    | Sample2 | ||

    """)

uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
if uploaded_file is not None:
    try:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, index_col = 0)
        elif input_file_type == "auto":
            df = read_csv(uploaded_file, index_col = 0, sep = None)
        elif input_file_type == "excel":
            df = read_excel(uploaded_file)
        else:
            df = read_csv(uploaded_file, sep = '\t', index_col = 0)
    except:
        df = read_excel(uploaded_file)
      #  st.write(df.head())


    content = df.columns.tolist()
    Gene_column = content[0]

    if "Annotation/Divergence" in content:
         # colnamesの変換
        search_word = r'([^\ \(]*)\ \(.*'

        for i in range(1, len(content)):
            match = re.search(search_word, content[i])
            if match:
                content[i] = match.group(1).replace(' ', '_')
        df.columns = content # 一旦名前を変更
        df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # excel 対応
        pattern = "([^|]*)"
        repatter = re.compile(pattern)
        f_annotation = lambda x: repatter.match(x).group(1)
        df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
      #  st.write(df.head())
        # annotation/divergence以前を除く
        df = df.loc[:,'Annotation/Divergence':]
      #  st.write(df.head())
        content = df.columns.tolist()
        content[0] = 'Gene'
        df.columns = content
        st.write("Converted Annotation/Divergence to gene symbols.")
        df.set_index("Gene", inplace = True)

    if input_file_format == 'column = gene':
        df = df.T

    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio( "NaN:",
        ('remove Nan containing genes', 'conver to 0' ), key='remove Nan containing genes')
        if convert_nan == "conver to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')


    st.write('Original:')
    st.write(df.iloc[:4,:7])
    df_sum = pd.DataFrame(df.sum())
    df_sum.columns = ['Total counts']
    fig, ax = plt.subplots()
    ax = sns.barplot(data=df_sum.T)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    ax.tick_params(axis='y', labelsize=8)
    ax.set_ylabel("Total counts", fontsize = 14)
    st.pyplot(fig)
    st.write(df_sum)

    if any(df.sum() == 0): # count 0の列を除く
        st.markdown('#### There are the samples that have zero counts.')
        st.write(", ".join(df.columns[df.sum()  == 0].to_list()))
        st.write('They are removed. Now data are:')
        df = df.drop(df.columns[df.sum()  == 0].to_list(), axis = 1)
        st.write(df.iloc[:4,:7])


        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']
        fig, ax = plt.subplots()
        ax = sns.barplot(data=df_sum.T)
#        ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
        ax.tick_params(axis='y', labelsize=8)

        ax.set_ylabel("Total counts", fontsize = 14)
        st.pyplot(fig )

    fig, ax = plt.subplots()
    ax = sns.boxplot(data = np.log1p(df))
    ax.set_ylabel("ln(x+1)", fontsize = 14)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
    ax.tick_params(axis='y', labelsize=8)
    st.pyplot(fig )


    show_cor = st.checkbox('Show correlation coeficient matrix?')
    if show_cor:
        correlation_coefficients = df.corr()
        fig_c, ax_c = plt.subplots() #この形式でないとエラーになる
        ax_c = sns.heatmap(correlation_coefficients, vmax=1, vmin=-1, cmap='seismic', square=True,
            annot=False, xticklabels=1, yticklabels=1)
        st.pyplot(fig_c)

    with st.sidebar:
        st.markdown("### Filter out weakly-expressed samples")
        st.markdown("#### Total count min:")
        min_sample_threshold = st.number_input("count total minimum", value = 0.0, label_visibility = 'collapsed')
       # min_sample_threshold = float(min_sample_threshold)
        st.markdown("---")
        st.markdown("##### Filter the genes > counts in all samples:")
        min_threshold = st.number_input("count minimum", value = 0.0, label_visibility = 'collapsed')
        min_threshold = float(min_threshold)
        st.markdown("##### Filter the genes > counts in at least n samples:")
        max_threshold = st.number_input("count max", value = 0.0, label_visibility = 'collapsed')
        max_n = st.number_input("count max", value = 1, label_visibility = 'collapsed')
        max_threshold = float(max_threshold)


    if min_sample_threshold > 0:
        org_len = df.shape[1]
        org_column = df.columns.tolist()
        df = df.loc[:, df.sum() > min_sample_threshold]
        st.write("Filtered after total min threshold:")
        st.write(df.head(3))
        st.write(f"From {org_len} samples to {df.shape[1]} samples.")
        if df.shape[1] < org_len:
            st.markdown(f"##### {', '.join(list(set(org_column) - set(df.columns.tolist())))} removed.")
            st.write(f'Remaining: {", ".join(df.columns.tolist())}')

    if min_threshold > 0:
        org_len = len(df)
        df = df[df.apply(min, axis=1) > min_threshold]
        st.write("Filtered after min threshold:")
        st.write(f"From {org_len} to {len(df)} genes.")

    if max_threshold > 0:
        org_len = len(df)
        df = df[(df > max_threshold).sum(axis=1)  >= max_n]
#           df = df[df.apply(max, axis=1) > max_threshold]
        st.write("Filtered after max threshold:")
        st.write(f"From {org_len} to {len(df)} genes.")

    with st.form("Set_method and transformation"):
        st.markdown("##### 1. Normalization method:")
        convert_to = st.radio(
        "",
        ('TMM','UQ','CTF', 'CUF', 'CPM', 'rlog', 'vst', 'RPKM to TPM', 'None'), key='TMM', label_visibility='collapsed')
        with st.expander("ℹ️ About normalization methods"):
            st.write("using rnanorm package: https://rnanorm.readthedocs.io/en/latest/index.html")
            st.write("For CTF and CUF: https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02568-9")

        st.markdown("##### 2. Log transformation:")
        log_transform = st.radio(
        "",
        ('None', 'log2(x+1)', 'ln(x+1)', 'log10(x+1)', 'asinh (hyperbolic arcsine transformation)'), key='None', label_visibility='collapsed')
        st.write("###### Log transformation does not apply to rlog or vst.")
        st.write("  ")
        st.markdown("##### 3. Standardization:")
        zscore = st.checkbox(" Z-score transformation", value = False)
        st.write("  ")

        submitted = st.form_submit_button("Set options and go")


    if submitted or st.session_state.use_R:

        if convert_to in ['rlog', 'vst']:
            log_transform = "None" # log_transofrmをキャンセル

        if convert_to == "CTF":
            from rnanorm import CTF
            df_conv = CTF().set_output(transform="pandas").fit_transform(df.T)
        if convert_to == "TMM":
            from rnanorm import TMM
            df_conv = TMM().set_output(transform="pandas").fit_transform(df.T)

        if convert_to == "CUF":
            from rnanorm import CUF
            df_conv = CUF().set_output(transform="pandas").fit_transform(df.T)

        if convert_to == "UQ":
            from rnanorm import UQ
            df_conv = UQ().set_output(transform="pandas").fit_transform(df.T)

        if convert_to == "CPM":
            from rnanorm import CPM
            df_conv = CPM().set_output(transform="pandas").fit_transform(df.T)

        if convert_to == "RPKM to TPM":
            df_conv = normalize_totalreads(df).T

        if convert_to == "None":
            df_conv = df.copy()

        if (input_file_format == 'row = gene') and (convert_to not in ['rlog', 'vst', 'None']):
            df_conv = df_conv.T

        rlog_finish = False
        if convert_to in ['rlog', 'vst']:
            st.session_state.use_R = True
            import rpy2.robjects as ro
            from rpy2.robjects.packages import importr
            from rpy2.robjects import pandas2ri
            from rpy2.robjects.vectors import StrVector
            import pyper
            r = pyper.R(use_pandas=True)
            f = ro.r("source('/home/cellxgene/streamlit/pages/deseq2_func.R')") # full pathが必要

            condition = [str(i) for i in df.columns.tolist()] #error防止
            group_condition = [remove_after_space(x) for x in condition] #スペース以降を除く
            group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字を除く
            df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])
            df_b = pd.DataFrame({'Batch': group_condition}, index=condition)
            df = df.astype(float)
            df = df.round(0)
            batch = None
            group = st.checkbox('Set group?', value=True)
            use_batch = st.checkbox('Set batch and use limma::removeBatchEffect?', value=False)
            if group or use_batch:
                with st.form("Set_groups and batch"):
                    if group:
                        edited_df_e = st.data_editor(df_e)
                        condition = edited_df_e.iloc[:,0].tolist()
                        st.write('Group: ' + '  '.join(condition))
                    else:
                        condition = df.columns.tolist()
                        batch = None

                    if use_batch:
                        edited_df_b = st.data_editor(df_b)
                        batch = edited_df_b['Batch'].tolist()
                        st.write('Group: ' + '  '.join(condition))
                        st.write('Batch: ' + '  '.join(batch))
                    else:
                        batch = None

                    group_submitted = st.form_submit_button("Set group/batch")

            if st.button('Run calc'):
                temp_dir = "temp/" + str(round(time.time()))

                if not os.path.exists('temp'):
                    os.mkdir('temp')
                else:
                    clear_old_directories("temp")
                    clear_old_files("temp")
                os.mkdir(temp_dir)

                r.assign('df',df)
                pyper_df_path = "saveRDS(df, '" + temp_dir + "/pyper_df.RDS')"
                r(pyper_df_path)
                read_pyper_df = "cts <- readRDS('" + temp_dir + "/pyper_df.RDS')"
                ro.r(read_pyper_df)


                #まずベクターに変換
                r_condition =  ro.StrVector(condition)
                ro.r.assign('condition', r_condition)
                ro.r.assign('temp_dir', temp_dir)


                ro.r("make_coldata2()")

                if convert_to == 'rlog':
                    if group:
                        ro.r('calc_rlog()')
                    else:
                        ro.r('calc_rlog_no_group()')
                else:
                    if group:
                        ro.r('calc_vst()')
                    else:
                        ro.r('calc_vst_no_group()')

                # df_conv = ro.conversion.rpy2py(rld) うまくいかない
                rld_path = temp_dir + '/rld.tsv'
                df_conv = pd.read_csv(rld_path, sep = '\t', header = 0)
                content = df_conv.columns.tolist()
                content[0] = 'Gene'
                df_conv.columns = content
                df_conv.set_index("Gene", inplace = True)
          #      os.unlink(temp_dir + "/rld.tsv")


                # Perform batch effect removal if batch information is provided
                if batch and any(batch):
                    from sklearn.decomposition import PCA
          #          st.write(batch)
#                    st.write(f'''batch <- c({', '.join(batch)})''')
                    batch_code = '"' + '", "'.join(batch) + '"'
                    st.write("Batch to R:")
                    st.write(batch_code)

                    ro.r(f'''batch <- c({batch_code})''')#うまく渡せないのでベタ書き

                    ro.r(f'''
                    library(limma)
                    print('batch')
                    print(batch)
                    rld_data <- as.matrix(read.table("{rld_path}", header=TRUE, row.names=1, sep="\t"))
                    batch_vector <- as.factor(batch)
                    rld_batch_removed <- removeBatchEffect(rld_data, batch=batch_vector)
                    write.table(rld_batch_removed, file=file.path(temp_dir, "rld_batch_removed.tsv"), sep="\t", quote=FALSE)
                    ''')
                    df_conv_batch_removed = pd.read_csv(temp_dir + '/rld_batch_removed.tsv', sep='\t', header=0, index_col=0)

                    # Create scatter plot to compare before and after batch effect removal
                    fig, ax = plt.subplots(figsize=(10, 10))
                    sns.scatterplot(x=df_conv.values.flatten(), y=df_conv_batch_removed.values.flatten(), alpha=0.1)
                    plt.xlabel('Before batch effect removal', fontsize=16)
                    plt.ylabel('After batch effect removal', fontsize=16)
                    plt.title('Comparison before and after batch effect removal')
                    st.pyplot(fig)


                    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
                    pca = PCA(n_components=2)
                    pca_result_before = pca.fit_transform(df_conv.T)
                    pca_result_after = pca.fit_transform(df_conv_batch_removed.T)
                    
                    sns.scatterplot(x=pca_result_before[:, 0], y=pca_result_before[:, 1], hue=batch, ax=ax1, s=400, alpha=0.5)
                    ax1.set_title('PCA Before Batch Effect Removal', fontsize=30)
                    ax1.set_xlabel('PC1', fontsize=24)
                    ax1.set_ylabel('PC2', fontsize=24)
                    
                    sns.scatterplot(x=pca_result_after[:, 0], y=pca_result_after[:, 1], hue=batch, ax=ax2, s=400, alpha=0.5)
                    ax2.set_title('PCA After Batch Effect Removal', fontsize=30)
                    ax2.set_xlabel('PC1', fontsize=24)
                    ax2.set_ylabel('PC2', fontsize=24)
                    
                    plt.tight_layout()
                    st.pyplot(fig)

                    df_conv = df_conv_batch_removed

                os.unlink(temp_dir + "/rld.tsv")

                if os.path.exists(temp_dir + "/rld_batch_removed.tsv"):
                    os.unlink(temp_dir + "/rld_batch_removed.tsv")

                rlog_finish = True
                log_transform = "None"

        else:
            st.session_state.use_R = False # rlog, vstでない

        if (convert_to not in ['rlog', 'vst']) or (rlog_finish):

            st.write('Converted:')
            st.write(df_conv.iloc[:4,:7])


            log_transform_word = ''
            if log_transform  != "None":
                if log_transform == 'asinh (hyperbolic arcsine transformation)':
                    df_conv = np.arcsinh(df_conv)
                    log_transform_word = ".asinh"
                if log_transform == 'log2(x+1)':
                    df_conv = np.log2(df_conv+1)
                    log_transform_word = ".log2"
                if log_transform == 'log10(x+1)':
                    df_conv = np.log10(df_conv+1)
                    log_transform_word = ".log10"
                if log_transform == 'ln(x+1)':
                    df_conv = np.log1p(df_conv)
                    log_transform_word = ".ln"

                st.write('Transformed:')
                st.write(df_conv.iloc[:4,:7])


            fig, ax = plt.subplots()
            if log_transform != "None":
                ax = sns.boxplot(data = df_conv)
                ax.set_ylabel("transformed value", fontsize = 14)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
                ax.tick_params(axis='y', labelsize=8)

            else:
                ax = sns.boxplot(data = np.log1p(df_conv))
                ax.set_ylabel("ln(transformed value+1)", fontsize = 14)
                ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
                ax.tick_params(axis='y', labelsize=8)

            ax.set_xticklabels(ax.get_xticklabels(), rotation=90, fontsize=8)
            ax.tick_params(axis='y', labelsize=8)
            st.pyplot(fig )
            
            # RLE plot for normalized data
            st.markdown("### RLE (Relative Log Expression) Plot - Normalized Data")
            
            with st.expander("ℹ️ RLE Plotについて"):
                st.info("""
                **RLE (Relative Log Expression) Plot**とは：
                
                各サンプルの発現量を、全サンプルの中央値と比較したプロットです。
                
                - **計算方法**: 各遺伝子について、全サンプルの中央値を計算し、各サンプルの値から中央値を引いた値をプロット
                - **理想的な状態**: すべてのボックスが0付近に中心を持ち、幅が狭い
                - **問題がある場合**: 特定のサンプルが0から大きく離れている、または幅が広い
                
                このプロットは、正規化の効果やバッチ効果の確認に有用です。
                """)
            # Use the normalized data (df_conv) for RLE calculation
            if log_transform != "None":
                # If already log-transformed, use the data as is
                log_norm_data = df_conv
            else:
                # If not log-transformed, apply log transformation for RLE calculation
                log_norm_data = np.log1p(df_conv)
            
            # Calculate median expression for each gene across all samples
            median_expression = log_norm_data.median(axis=1)
            
            # Calculate RLE values
            rle_data = log_norm_data.subtract(median_expression, axis=0)
            
            # Create RLE plot
            fig_rle, ax_rle = plt.subplots(figsize=(10, 6))
            ax_rle = sns.boxplot(data=rle_data)
            ax_rle.set_ylabel("RLE", fontsize=14)
            ax_rle.set_xlabel("Samples", fontsize=14)
            ax_rle.set_title("Relative Log Expression (RLE) Plot - After Normalization", fontsize=16)
            ax_rle.axhline(y=0, color='red', linestyle='--', alpha=0.5)
            ax_rle.set_xticklabels(ax_rle.get_xticklabels(), rotation=90, fontsize=8)
            ax_rle.tick_params(axis='y', labelsize=8)
            plt.tight_layout()
            st.pyplot(fig_rle)

            if show_cor:
                correlation_coefficients = df.corr()
                fig_c, ax_c = plt.subplots() #この形式でないとエラーになる
                ax_c = sns.heatmap(correlation_coefficients, vmax=1, vmin=-1, cmap='seismic', square=True,
                    annot=False, xticklabels=1, yticklabels=1)
                st.pyplot(fig_c)


            if convert_to == "RPKM to TPM":
                file_name = os.path.splitext(uploaded_file.name)[0] + log_transform_word +'.TPM.txt'
                st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0] + log_transform_word +'.TPM'
            elif convert_to == "None":
                file_name = os.path.splitext(uploaded_file.name)[0] + '.filtered' + log_transform_word + '.txt'
                st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0] + '.filtered' + log_transform_word
            else:
                file_name = os.path.splitext(uploaded_file.name)[0] + '.' + convert_to + log_transform_word + '.txt'
                st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0] + log_transform_word

            st.session_state.df = df_conv

                # Mean-SD Comparison Plotを追加
            fig_mean_sd = plot_mean_sd_comparison(df, df_conv)
            st.pyplot(fig_mean_sd)


            if zscore:
                center0_z= True
                df_z = df_conv.copy()
                m = df_z.mean(1)
                s = df_z.std(1)
                df_z = df_z.sub(m, axis=0).div(s, axis = 0)
                df_z = np.round(df_z, decimals=10)
                df_z = df_z.loc[~(df_z==0).all(axis=1)] #すべて0のrowを除く
                df_z = df_z.dropna(how='any', axis=0) #エラー対応
                df_conv = df_z
                st.markdown("#### Normalized data are standardized.")
                
                # Z-score標準化後のデータを表示
                st.write('Z-score standardized:')
                st.write(df_conv.iloc[:4,:7])
                
                # Z-score標準化後のボックスプロット
                st.markdown("##### Distribution after Z-score standardization")
                fig_zscore, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
                
                # ボックスプロット
                ax1 = sns.boxplot(data=df_conv, ax=ax1)
                ax1.set_ylabel("Z-score", fontsize=14)
                ax1.set_xlabel("Samples", fontsize=14)
                ax1.set_title("Z-score Distribution by Sample", fontsize=16)
                ax1.axhline(y=0, color='red', linestyle='--', alpha=0.5)
                ax1.axhline(y=1, color='gray', linestyle=':', alpha=0.5)
                ax1.axhline(y=-1, color='gray', linestyle=':', alpha=0.5)
                ax1.set_xticklabels(ax1.get_xticklabels(), rotation=90, fontsize=8)
                ax1.tick_params(axis='y', labelsize=8)
                
                # ヒストグラム（全データの分布）
                ax2.hist(df_conv.values.flatten(), bins=50, edgecolor='black', alpha=0.7)
                ax2.set_xlabel("Z-score", fontsize=14)
                ax2.set_ylabel("Frequency", fontsize=14)
                ax2.set_title("Overall Z-score Distribution", fontsize=16)
                ax2.axvline(x=0, color='red', linestyle='--', alpha=0.5, label='Mean=0')
                ax2.axvline(x=1, color='gray', linestyle=':', alpha=0.5, label='±1 SD')
                ax2.axvline(x=-1, color='gray', linestyle=':', alpha=0.5)
                ax2.legend()
                
                plt.tight_layout()
                st.pyplot(fig_zscore)
                
                # 統計情報を表示
                with st.expander("📊 Z-score statistics"):
                    col1, col2 = st.columns(2)
                    with col1:
                        st.write("**Per-sample statistics:**")
                        sample_stats = pd.DataFrame({
                            'Mean': df_conv.mean(axis=0),
                            'Std': df_conv.std(axis=0),
                            'Min': df_conv.min(axis=0),
                            'Max': df_conv.max(axis=0)
                        })
                        st.dataframe(sample_stats.round(3))
                    with col2:
                        st.write("**Overall statistics:**")
                        st.write(f"Global mean: {df_conv.values.mean():.6f}")
                        st.write(f"Global std: {df_conv.values.std():.6f}")
                        st.write(f"Min value: {df_conv.values.min():.3f}")
                        st.write(f"Max value: {df_conv.values.max():.3f}")

            csv = convert_df(df_conv)
            st.download_button(
               "Press to Download",
               csv,
               file_name,
               "text/csv",
               key='download-csv'
            )

