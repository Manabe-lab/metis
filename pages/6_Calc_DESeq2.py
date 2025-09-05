#!!!!!!!!!!!!!! pip install rpy2==3.5.1  新しいバージョンはエラーが出る

# 基本的にglobal変数で計算する。
# pythonからassgnされるのはglobal変数


import streamlit as st
import rpy2
import csv
import re
import os
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector
import pyper
import shutil
from PIL import Image
from helper_func import clear_old_directories, clear_old_files, remove_after_space, remove_sample_num
import time
import sys
from rpy2.robjects.conversion import localconverter

from collections import Counter

def clean_column_names_for_r(df):
    """
    Convert DataFrame column names to R-safe format
    """
    original_columns = df.columns.tolist()
    cleaned_columns = []
    mapping = {}
    
    for col in original_columns:
        # R安全な文字（英数字とアンダースコア）以外をピリオドに置換
        cleaned = re.sub(r'[^a-zA-Z0-9_]', '.', str(col))
        # 先頭が数字の場合はXを追加
        if cleaned and cleaned[0].isdigit():
            cleaned = 'X' + cleaned
        # 連続するピリオドを1つに
        cleaned = re.sub(r'\.+', '.', cleaned)
        # 末尾のピリオドを削除
        cleaned = cleaned.rstrip('.')
        
        cleaned_columns.append(cleaned)
        if col != cleaned:
            mapping[col] = cleaned
    
    # データフレームのカラム名を変更
    df.columns = cleaned_columns
    
    return df, mapping, original_columns

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []    
    # 最も短い文字列の長さを取得
    min_length = min(len(s) for s in strings)
    # 共通の末尾部分の長さを見つける
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break            
    # 共通の末尾部分が見つからない場合は元のリストを返す
    if suffix_length == 0:
        return strings        
    # 共通の末尾部分を削除して新しいリストを作成
    return [s[:-suffix_length] for s in strings]


def rename_duplicates(df):
    """
    Rename duplicate indices by adding _2, _3, etc. to subsequent occurrences   
    Args:
        df: pandas DataFrame
    Returns:
        DataFrame with renamed indices
    """
    # Get current index values
    lis = df.index.values
    
    # Count occurrences of each value
    counts = Counter()
    new_indices = []
    
    for x in lis:
        counts[x] += 1
        if counts[x] == 1:
            new_indices.append(x)
        else:
            new_indices.append(f"{x}_{counts[x]}")
    
    # Check if there were any duplicates
    if len(lis) != len(set(lis)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene names of subsequent occurrences have _2, _3, etc. at the end.")
        
        # Display which names were changed
        for name, count in counts.items():
            if count > 1:
                st.write(f"'{name}' appears {count} times → {name}, " + 
                        ", ".join([f"{name}_{i}" for i in range(2, count + 1)]))
    
    # Set new index
    df.index = new_indices
    return df

#March-1 Sept-1対応
def excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.write("There are Excel-autoconverted gene names")
                k = 1
            autoconvert_flag = True
            st.write("Converting " + i)
            if x.group(2) == "Mar":
                index_name[j] = "March" + x.group(1)
            elif x.group(2) == "Sep":
                index_name[j] = "Sept" + x.group(1)
        j += 1
    dfx.index = index_name
    return(dfx)


def check_excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep|Oct|Dec|Feb|Nov)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.markdown("#### There are Excel-autoconverted gene names")
                st.write("Gene names are not converted.")
                k = 1
            autoconvert_flag = True
            st.write(i)
#    return(dfx)

r = pyper.R(use_pandas=True)
f = ro.r("source('pages/deseq2_func.R')") # full pathが必要

st.set_page_config(page_title="Calculate DESeq2.", page_icon="📃")


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c



@st.cache_data
def calc_barplot(data, ylabel):
    fig, ax = plt.subplots()
    ax = sns.barplot(data=data)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    ax.set_ylabel(ylabel, fontsize = 14)
    return fig

st.markdown("### DESeq2, limma eBayes, beta regression, GLM")
st.sidebar.title("Options")
st.markdown("#### Options are displayed at the bottom of the left side panel")
with st.sidebar:
    st.markdown("### Analysis Method:")
    test_method = st.radio("Select analysis method:", 
                         ["DESeq2", "limma eBayes", "Beta Regression", 
                          "Generalized Linear Model (GLM)"], 
                         index=0)
    
    st.markdown("###### limma eBayes with logit transformation, beta regression and GLM with beta regression are for proportion data.")

# temp内に保存する
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    #古いdirecotryとファイルを削除する
    temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)

else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)


st.markdown("### DESeq2にはraw count dataを使う")

use_sf = False # size factorの使用

use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Available data")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        # Homer対応
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list(): # Row_nameを含むとき
            df = df.set_index('Row_name')
            df.index.name = "Gene"

uploaded_size_factors = None
if "use_custom_size_factors" not in st.session_state:
    st.session_state.use_custom_size_factors = False

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('auto', 'Homer','tsv','csv','excel'), index = 0, label_visibility = 'collapsed')
    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])
    use_sf = st.checkbox('Upload Size Factors (Optional)')
    if use_sf:
        uploaded_size_factors = st.file_uploader("Choose a size factors file (TSV format)", type=['tsv'])

    if uploaded_file is not None:

        if file_type == 'auto':
            try:
                df = read_csv(uploaded_file, sep = None)
                st.write("Uploaded file:")
                st.write(df.head())

                content = df.columns.tolist()
#                Gene_column = content[0]

                if "Annotation/Divergence" in content:
                     # colnamesの変換
                    search_word = '([^\ \(]*)\ \(.*'

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
                    st.write("Converted Annotation/Divergence to gene symbols.")
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content

         #       df.set_index("Gene", inplace = True)

            except:# excel
                df = read_excel(uploaded_file)
                content = df.columns.tolist()
                if "Annotation/Divergence" in content:
                     # colnamesの変換
                    search_word = '([^\ \(]*)\ \(.*'

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
                    # annotation/divergence以前を除く
                    df = df.loc[:,'Annotation/Divergence':]
                    content = df.columns.tolist()
                    content[0] = 'Gene'
                    df.columns = content
                    st.write("Converted Annotation/Divergence to gene symbols.")
                else:
                    colnames = df.columns.tolist()
                    colnames[0] = 'Gene'
                    df.columns = colnames

        elif file_type != 'excel':
            if file_type == 'csv':
                df = read_csv(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep = '\t')
            st.write("Original:")
            st.write(df.head())
            if file_type == 'Homer':
                df = df.iloc[:,7:]
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                # colnamesの変換
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(colnames)):
                    match = re.search(search_word, colnames[i])
                    if match:
                        colnames[i] = match.group(1).replace(' ', '_')
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                try:
                    df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)
                    df.columns = colnames
                except:
                    st.markdown("### File format error. Non-Homer file?")

            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames
        else: # excel
            df = read_excel(uploaded_file)
            content = df.columns.tolist()
            if "Annotation/Divergence" in content:
                 # colnamesの変換
                search_word = '([^\ \(]*)\ \(.*'

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
                # annotation/divergence以前を除く
                df = df.loc[:,'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                st.write("Converted Annotation/Divergence to gene symbols.")
            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames

        df = df.set_index('Gene')
        file_name_head = os.path.splitext(uploaded_file.name)[0]

    else:
        sys.exit(1)

if df is not None:
    st.write('Original gene number:  ' + str(len(df)))
    st.write(df.head())

    # floatに変換 誤射悟入
    df = df.astype(float)

    if test_method == "DESeq2": # DESeq2のときだけ整数化
        if not float.is_integer(df.iloc[:,0].sum()*1000):
            st.markdown("# It is likely that your data are normalized. Please upload unnormalized raw count data.")

        df = df.round(0)

    df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く

########## excel対応?
    st.write("All zero count genes are removed.")
    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio( "NaN:",
        ('remove Nan containing genes', 'conver to 0' ), key='remove Nan containing genes')
        if convert_nan == "conver to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')
############
    check_excel_autoconversion(df)


    if len(df.index.values) != len(set(df.index.values)):
#        st.markdown("#### There are duplicated rows. Converting the names...")
#        st.write("The gene name of the second occurrence has _2 at the end.")
#        lis = df.index.values
#        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]
        df = rename_duplicates(df)
    # ここに新しい処理を追加
    df, column_name_mapping, original_column_names = clean_column_names_for_r(df)

    # 変換があった場合は警告を表示
    if column_name_mapping:
        st.warning("⚠️ Special characters in sample names have been converted for R compatibility:")
        for orig, clean in column_name_mapping.items():
            st.write(f"  • '{orig}' → '{clean}'")

    st.write(df.head())
    total_count = pd.DataFrame(df.sum()[1:])
    total_count.columns= ['Total counts']
    large_var = False
    if max(total_count['Total counts']) > min(total_count['Total counts']) * 2:
        large_var = True
        st.markdown("### Large difference (>2x) in counts")
        st.write(f"Minimum total counts: {min(total_count['Total counts'])}")
        st.markdown("##### Low count samples can be filtered on the side panel.")
        import matplotlib.pyplot as plt
        import seaborn as sns
        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']


        f1 = calc_barplot(df_sum.T, ylabel = "Total counts")
        st.pyplot(f1)

        f2 = calc_barplot(np.log1p(df), ylabel = "ln(x+1)")
        st.pyplot(f2)

    with st.sidebar:
        # DESeq2の既存のオプション
        if test_method == 'DESeq2':
            st.markdown("##### FC adjustment method:")
            type = st.radio("", ('ashr','apeglm', 'normal'), label_visibility="collapsed")
            st.write('DESeq2のdefaultはapeglm。apeglmではreference groupを指定する必要がある。')

            st.markdown("##### FDR cutoff for independent filtering:")
            results_alpha = st.number_input("alpha", value = 0.05, max_value=0.20, min_value=0.00, label_visibility = 'collapsed')
            st.write("alpha in results func")  

            st.markdown("##### Batch correction:")
            sva = st.checkbox('SVA batch removal?')
            sva_calc = True
            if sva:
                sva_calc = st.checkbox('Calculate only 2 surrogate variables? Deselect if want to calculate up to the recommended number.', value = True)
                st.markdown("---")

            ruv = st.checkbox('RUV batch removal?')

            if ruv:
                RUV_alpha = st.number_input('P values threshold for control genes in RUV', min_value=0.0, max_value = 0.5, step = 0.05, value=0.2)
            else:
                RUV_alpha = 0.2

        # limma eBayesのオプション
        elif test_method == 'limma eBayes':
            limma_data = st.radio("Data type:",
                ["RNA-seq count", "Non-count data", "0-1 data (proportion, AUC etc) to logit transformation"],
                index=1)

            if limma_data == "RNA-seq count":
                apply_logit = False
                limma_count = True
                default_trend = True  # RNA-seq countの場合はtrendをデフォルトでTRUE
            elif limma_data == "Non-count data":
                apply_logit = False
                limma_count = False
                default_trend = False  # Non-count dataの場合はtrendをデフォルトでFALSE
            else:
                apply_logit = True
                limma_count = False
                default_trend = False  # 0-1 dataの場合はtrendをデフォルトでFALSE
            
            # trend and robust options
            st.markdown("##### Advanced eBayes options:")
            limma_trend = st.checkbox("Use trend", value=default_trend)
            limma_robust = st.checkbox("Use robust", value=True)
            
            with st.expander("ℹ️ About trend and robust options"):
                st.markdown("""
                **trend option:**
                - 分散と平均発現量の関係をモデル化します
                - 低発現遺伝子と高発現遺伝子で分散が異なる場合に有効
                - RNA-seqデータで、voomやlog変換後のデータに推奨
                - すべての遺伝子で一定の分散を仮定せず、より柔軟なモデルを適用
                
                **robust option:**
                - 外れ値に対してロバストな推定を行います
                - ベイズ推定で外れ値の影響を軽減
                - データにアーティファクトや外れ値が含まれる可能性がある場合に推奨
                - 計算時間は少し長くなりますが、より安定した結果を得られます
                
                **推奨される使用例:**
                - RNA-seqカウントデータ: trend=TRUE, robust=TRUE
                - 非カウントデータ（マイクロアレイ等）: trend=FALSE, robust=TRUE
                - 0-1データ（比率データ）: trend=FALSE, robust=TRUE
                - きれいな正規化済みデータ: 両方FALSE（デフォルト）
                
                **データタイプ別の推奨設定：**
                - **"RNA-seq count"**: trend=TRUE, robust=TRUE (自動設定)
                  - カウントデータは平均-分散関係があり、外れ値も存在する可能性が高い
                - **"Non-count data"**: trend=FALSE, robust=TRUE (自動設定)
                  - 正規化済みデータは平均-分散関係が安定化済み
                  - 外れ値に対するロバスト性は依然として有用
                - **"0-1 data"**: trend=FALSE, robust=TRUE (自動設定)
                  - 比率データは変換済みで平均-分散関係は安定
                  - ロバスト推定で外れ値の影響を軽減
                
                **注意：** 
                - trendとrobustの初期値は、選択したデータタイプに基づいて自動設定されます
                - robustオプションは常にTRUEがデフォルトです（より安定した結果）
                - 必要に応じて手動で設定を変更できます
                """)

        # Beta Regressionのオプション
        elif test_method == 'Beta Regression':
            st.markdown("### Beta Regression Options:")
            epsilon = st.number_input("Epsilon for boundary adjustment (0-1 data)", 
                                    min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
            st.markdown("#### Batch correction:")
            use_batch = st.checkbox('Include batch effect?', value=False)
            n_cores = st.slider("Parallel cores", min_value=1, 
                               max_value=os.cpu_count()-1, 
                               value=max(1, os.cpu_count()//2-4))

        # GAMのオプション
        elif test_method == 'Generalized Linear Model (GLM)':
            st.markdown("### GLM Options:")
            dist_family = st.radio("Probability distribution", 
                                  ["Beta (0-1)", "Gaussian", "Poisson", "Negative Binomial"],
                                  index=0)


            with st.expander("Explain models"):
                st.markdown("""
GLMをmgcvパッケージのgam()関数で実装

1. Beta分布 (Beta Distribution)
適用データ：0と1の間の値（正規化発現量、比率データ）  
リンク関数：通常ロジットリンク（logit）を使用  
推定：最尤法でパラメータ推定（Beta分布の形状パラメータα、β）  
検定：条件の係数についてのWald検定（係数/標準誤差）を実施  
多重検定補正：Benjamini-Hochberg法でFDRを制御  
Beta回帰とGAM-Beta実装は本質的に同じアルゴリズム  
batchを含まない場合は本質的には「0-1データに適応したANOVA」

2. ガウス分布 (Gaussian Distribution)  
適用データ：連続値で正規分布に従うデータ（正規化済みのログカウントなど）  
リンク関数：恒等リンク（identity）  
推定：最小二乗法（または最尤法）でパラメータ推定  
検定：t検定に基づいてグループ間の差を評価  
多重検定補正：Benjamini-Hochberg法でFDR制御  
標準的な線形モデル（ANOVA的アプローチ）  
batchを含まない場合は標準的なANOVAと一致（F検定によるグループ間差異の評価）

3. ポアソン分布 (Poisson Distribution)  
適用データ：単純なカウントデータ（過分散がない場合）  
リンク関数：対数リンク（log）  
推定：最尤法でパラメータ推定（λ：平均=分散）  
検定：尤度比検定または Wald検定  
多重検定補正：Benjamini-Hochberg法でFDR制御  

4. 負の二項分布 (Negative Binomial Distribution)  
適用データ：過分散のあるカウントデータ（RNA-seqデータなど）  
リンク関数：対数リンク（log）  
推定：最尤法でパラメータ推定（μ：平均、θ：過分散パラメータ）  
過分散推定：遺伝子ごとの分散をモデル化（DESeq2と類似）  
検定：尤度比検定または Wald検定  
多重検定補正：Benjamini-Hochberg法でFDR制御  
ポアソン分布の拡張で、variance > mean を許容
                """)

            if dist_family == "Beta (0-1)":
                dist_short = "beta"
                epsilon = st.number_input("Epsilon for boundary adjustment", 
                                       min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
            elif dist_family == "Gaussian":
                dist_short = "gaussian"
            elif dist_family == "Poisson":
                dist_short = "poisson"
            elif dist_family == "Negative Binomial":
                dist_short = "nb"
                nb_theta = st.number_input("過分散パラメータ (theta)", 
                                         min_value=0.1, max_value=100.0, value=10.0)

            st.markdown("#### Batch correction:")
            use_batch = st.checkbox('Include batch effect?', value=False)
            n_cores = st.slider("Parallel cores", min_value=1, 
                               max_value=os.cpu_count()-1, 
                               value=max(1, os.cpu_count()//2-4))

        # 共通のオプション
        independentFiltering = True  # デフォルト値を設定
        if test_method == 'DESeq2' or (test_method == 'limma eBayes' and limma_count):
            st.markdown("#### Filter out weakly-expressed genes before multiple test correction:",help = "independentFiltering default:TRUE 平均正規化カウントに基づいて遺伝子をフィルタリングし、多重検定補正の負担を減らすことで統計的検出力を向上させる")
            independentFiltering = st.checkbox('Yes', value=True)

        st.markdown("#### Additional filtering:")
        st.markdown("##### Filter the genes > counts in all samples:")
        min_threshold = st.number_input("count minimum", value=0, label_visibility='collapsed')
        min_threshold = int(min_threshold)
        
        st.markdown("##### Filter the genes > counts in at least one sample:")
        max_threshold = st.number_input("count max", value=0, label_visibility='collapsed')
        max_threshold = int(max_threshold)

        sample_threshold = 0
        remove_zero_samples = False

        if large_var:
            st.markdown("##### Filter the samples <= counts:")
            sample_threshold = st.number_input("Minimum total cout", value = 0, label_visibility = 'collapsed')

        # DESeq2以外の場合、ゼロカウントサンプル除去のオプションを表示
        if test_method != 'DESeq2':
            st.markdown("##### Remove samples with zero counts:")
            remove_zero_samples = st.checkbox("Remove samples with total count = 0", value=False)
        else:
            # DESeq2の場合は常に除去
            remove_zero_samples = True

        st.markdown("##### Output style:")
        deseq2_flag = st.radio( "Homer: positive in A vs B means up in B; DESeq2: positive in A vs B means up in A", ('Homer','DESeq2'), index = 1)
        if deseq2_flag =='DESeq2':
            deseq2 = True
        else:
            deseq2 = False

    # DESeq2の場合、または明示的に指定された場合にゼロカウントサンプルを除去
    if test_method == 'DESeq2' or remove_zero_samples:
        if any(df.sum() <= sample_threshold): # count 0の列を除く
            st.markdown('#### There are the samples that have counts <= ' + str(sample_threshold))
            st.write(", ".join(df.columns[df.sum() <= sample_threshold].to_list()))
            st.markdown('##### They are removed. Now data are:')
            df = df.drop(df.columns[df.sum() <= sample_threshold].to_list(), axis = 1)
            st.write(df.head())

    if min_threshold > 0:
        df = df[df.apply(min, axis=1) > min_threshold]
    if max_threshold > 0:
        df = df[df.apply(max, axis=1) > max_threshold]

    st.markdown(f'#### Filtered gene number: {str(len(df))}')

    condition = [str(i) for i in df.columns.tolist()] #error防止
    group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
#    group_condition = [remove_after_space(x) for x in condition] #スペース以降を除く
    group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字と_を除く
    group_condition = [x.rstrip('.') for x in group_condition] # .を除く


    df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])
    df_b = pd.DataFrame(condition, index =condition, columns = ["Batch"])

    batch_in = st.checkbox('Setting batch?')
    with st.form("input_groups and batch"):
        st.write('Set groups:')
    #    edited_df_e = st.experimental_data_editor(df_e)
        edited_df_e = st.data_editor(df_e)

        condition = edited_df_e.iloc[:,0].tolist()
        
        if batch_in:
            st.write('Set batch:')
    #        edited_df_b = st.experimental_data_editor(df_b)
            edited_df_b = st.data_editor(df_b)

        if batch_in:
            batch = edited_df_b.iloc[:,0].tolist()
            st.write('Batch: ' + '  '.join(batch))
        else:
            batch = ["No batch"] #batchがないとき

        submitted = st.form_submit_button("Submit")
    st.write('Group: ' + ' '.join(condition))

    if (len(condition) != len(df.columns)):
            st.write("The number of group name does not match the data.")

#    df_condition = pd.DataFrame(condition)
#    df_batch = pd.DataFrame(batch)

    ref_in = st.checkbox('Setting referece group?')
    if ref_in or (type == 'apeglm'):
        st.markdown("##### Select a group used as the referece:")
        ref_group = st.selectbox("", set(condition), label_visibility="collapsed")
    else:
        ref_group = condition[0]

    rld_calc = st.checkbox('Calculate rlog?', value = False )

    for i in df.columns.values:
        a = df.select_dtypes(exclude='number')
        if not a.empty:
            st.write("There is a non-nemeric value in ")
            st.write(a)

#    df = excel_autoconversion(df) # 1-Marなどの誤変換への対応


    st.markdown("---")
    # サイズファクターの設定と確認
    if use_sf and uploaded_size_factors is not None:
        size_factors_df = pd.read_csv(uploaded_size_factors, sep='\t', index_col=0)
        st.write("Uploaded size factors file:")
        st.write(size_factors_df.head())

        selected_column = st.selectbox("Select the column containing size factors:", size_factors_df.columns)

        if st.button("Confirm size factors"):
            if size_factors_df.index.tolist() == df.columns.tolist():
                size_factors = size_factors_df[selected_column]
                with localconverter(ro.default_converter + pandas2ri.converter):
                    ro.r.assign('custom_size_factors', ro.FloatVector(size_factors.values))
                ro.r('use_custom_size_factors <- TRUE')
                st.success(f"Size factors from column '{selected_column}' will be used in the analysis.")
                st.session_state.use_custom_size_factors = True

                # R変数の内容を確認
                display_r_variable("use_custom_size_factors")
                display_r_variable("custom_size_factors")
            else:
                st.error("The row names in the size factors file do not match the column names in the count data. Please check your file.")
                ro.r('use_custom_size_factors <- FALSE')
                ro.r('custom_size_factors <- NULL')
                st.session_state.use_custom_size_factors = False
        st.markdown("---")
    else:
        ro.r('use_custom_size_factors <- FALSE')
        ro.r('custom_size_factors <- NULL')
        st.session_state.use_custom_size_factors = False



    # 「Run DESeq2」ボタンを「Run Analysis」に変更
    if st.button('Run Analysis'):
        # 共通のセットアップ
        try:
            r.assign('df', df)
            pyper_df_path = "saveRDS(df, '" + temp_dir + "/pyper_df.RDS')"
            r(pyper_df_path)
            read_pyper_df = "cts <- readRDS('" + temp_dir + "/pyper_df.RDS')"
            ro.r(read_pyper_df)
            
            r_condition = ro.StrVector(condition)
            ro.r.assign('condition', r_condition)
            r_batch = ro.StrVector(batch)
            ro.r.assign('batch', r_batch)
            ro.r.assign('ref_in', ref_in)
            ro.r.assign('ref_group', ref_group)
            ro.r.assign('independentFiltering', independentFiltering)
            ro.r.assign('deseq2', deseq2)
            ro.r.assign('res_dir', res_dir)
            ro.r.assign('temp_dir', temp_dir)
            
            ro.r("make_coldata()")
            
        except Exception as e:
            st.markdown(f"## Error: {str(e)}")
            st.markdown("## rpy2 error. Reinstall: pip install rpy2==3.5.1")
            sys.exit(1)

        # メソッド固有の実行ロジック
        if test_method == 'DESeq2':
            # 既存のDESeq2コード
            ro.r.assign('sva', sva)
            ro.r.assign('ruv', ruv)
            ro.r.assign('RUV.alpha', RUV_alpha)
            ro.r.assign('type', type)
            ro.r.assign('sva_calc', sva_calc)
            ro.r.assign('rld_calc', rld_calc)
            ro.r.assign('results_alpha', results_alpha)
            
            with st.spinner('Calculating DESeq2...'):
                if batch == ['No batch']:
                    ro.r('calc_dds_nobatch()')
                else:
                    ro.r('calc_dds_batch()')
                py_res = ro.r('calc_deseq()')
            
            image = Image.open(res_dir + '/DispersionEstimates.png')
            st.image(image, caption='Despersion Estimates')
            
            with open(res_dir + '/DESeq2_output.txt') as f:
                out = f.readlines()
                for i in out:
                    if ("NULL" not in i) and ("results" not in i):
                        st.write(i)
            
            res_df = pd.read_csv(res_dir + '/DESeq2_res.tsv', sep='\t', index_col=0)
            res_df = res_df.fillna(1)
            st.dataframe(res_df)
            
            if sva:
                st.markdown("#### =======SVA=======")
                try:
                    with st.spinner('Preparing SVAseq...'):
                        sva_n = ro.r("sv_n <- calc_sva_n()")
                    st.write("Recommended number of SVA covariates: " + str(int(sva_n[0])))
                except:
                    pass

                                # 両方の方法で計算
                ro.r("sv_both <- calc_sva_n_both()")
                
                # 結果を取得
                sv_be = int(ro.r("svn_be")[0])
                sv_leek_result = ro.r("svn_leek")
                
                # Leek法の結果確認（NAの可能性があるため）
                if sv_leek_result[0] != ro.NA_Logical:
                    sv_leek = int(sv_leek_result[0])
                else:
                    sv_leek = None
            
                # 結果を表示
                col1, col2 = st.columns(2)
                
                with col1:
                    st.metric(
                        "BE method (default)", 
                        sv_be,
                        help="Buja-Eyuboglu method - tends to estimate more factors"
                    )
                    
                with col2:
                    if sv_leek is not None:
                        st.metric(
                            "Leek method", 
                            sv_leek,
                            help="More conservative estimation"
                        )
                    else:
                        st.metric("Leek method", "Failed")
                        

                # 警告やアドバイスを表示
                if sv_leek is not None:
                    if sv_be > 10:
                        st.warning(f"""
                        ⚠️ **High number of surrogate variables detected**
                        
                        - BE method suggests {sv_be} variables
                        - Leek method suggests {sv_leek} variables
                        
                        This may indicate:
                        - Strong batch effects in your data
                        - Potential overfitting risk
                        - Need for alternative approaches
                        """)
                        
                    if sv_be > sv_leek * 2:
                        st.info(f"""
                        💡 **Large discrepancy between methods**
                        
                        The BE estimate ({sv_be}) is more than twice the Leek estimate ({sv_leek}).
                        Consider:
                        1. Using the more conservative Leek estimate
                        2. Starting with fewer SVs and checking DEG retention
                        3. Investigating sample outliers
                        """)
                
                with st.spinner('Calculating SVAseq...'):
                    ro.r("calc_svseq()")
                
                with open(res_dir + '/SVA_output.txt') as f:
                    out = f.readlines()
                    for i in out:
                        if ("NULL" not in i)  and ("results" not in i):
                            st.write(i)
                
            if ruv:
                st.markdown("#### =======RUV=======")
                with st.spinner('Calculating RUVseq...'):
                    ro.r("calc_ruvseq()")
                
                with open(res_dir + '/RUVseq.txt') as f:
                    out = f.readlines()
                    for i in out:
                        if ("NULL" not in i)  and ("results" not in i):
                            st.write(i)
            
            # セッションにデータを残す
            if sva:
                st.session_state.deseq2 = read_csv(res_dir + "/SVA_res.tsv", sep='\t', index_col=0)
            elif ruv:
                st.session_state.deseq2 = read_csv(res_dir + "/RUV_res.tsv", sep='\t', index_col=0)
            else:
                st.session_state.deseq2 = read_csv(res_dir + "/DESeq2_res.tsv", sep='\t', index_col=0)
            
            file_name = file_name_head + "_DESeq2"
        
        elif test_method == 'limma eBayes':
            # limma eBayesの実行ロジック
            ro.r.assign('limma_count', limma_count)
            ro.r.assign('apply_logit', apply_logit)
            ro.r.assign('limma_trend', limma_trend)
            ro.r.assign('limma_robust', limma_robust)
            
            with st.spinner('Calculating limma eBayes...'):
                ro.r('calc_limma()')
            
            if limma_count:
                try:
                    image = Image.open(res_dir + '/voom_plot.png')
                    st.image(image, caption='Voom Mean-Variance Trend')
                except:
                    st.write("Voom plot not available")
            res_df = pd.read_csv(res_dir + '/limma_res.tsv', sep='\t', index_col=0)
            res_df = res_df.fillna(1)
            # adj.pvalueで終わる列を検索して各比較ごとに有意遺伝子数を表示
            adj_pvalue_columns = [col for col in res_df.columns if '.adj.pvalue' in col]
            st.markdown("### Significant genes (FDR < 0.05):")
            total_significant = 0
            for col in adj_pvalue_columns:
                comparison = col.split('.adj.pvalue')[0]  # 比較名を取得
                sig_count = (res_df[col] < 0.05).sum()
                total_significant += sig_count
                st.write(f"- {comparison}: {sig_count}")
            st.write(f"Total significant genes: {total_significant}")

            st.dataframe(res_df)
            
            st.session_state.deseq2 = res_df  # 結果をセッションに保存
            file_name = file_name_head + "_limma"
        
        elif test_method == 'Beta Regression':
            ro.r.assign('epsilon', epsilon)
            ro.r.assign('n_cores', n_cores)
            ro.r.assign('use_batch', use_batch)

            with st.spinner('Calculating Beta Regression...'):
                ro.r('calc_betareg()')

            # 結果の読み込みと表示
            try:
                res_df = pd.read_csv(res_dir + '/betareg_res.tsv', sep='\t', index_col=0)
                res_df = res_df.fillna(1)
                
                # adj.pvalueで終わる列を検索して各比較ごとに有意遺伝子数を表示
                adj_pvalue_columns = [col for col in res_df.columns if '.adj.pvalue' in col]
                st.markdown("### Significant genes (FDR < 0.05):")
                total_significant = 0
                for col in adj_pvalue_columns:
                    comparison = col.split('.adj.pvalue')[0]  # 比較名を取得
                    sig_count = (res_df[col] < 0.05).sum()
                    total_significant += sig_count
                    st.write(f"- {comparison}: {sig_count}")
                st.write(f"Total significant genes: {total_significant}")
                
                st.dataframe(res_df)
                
                st.session_state.deseq2 = res_df
            except Exception as e:
                st.error(f"Error processing Beta Regression results: {str(e)}")
                st.write("Check the R output for details.")

            file_name = file_name_head + "_betareg"
        
        elif test_method == 'Generalized Linear Model (GLM)':
            # GAMの実行ロジック
            ro.r.assign('dist_short', dist_short)
            ro.r.assign('use_batch', use_batch)
            ro.r.assign('n_cores', n_cores)
            
            if dist_short == "beta":
                ro.r.assign('epsilon', epsilon)
            elif dist_short == "nb":
                ro.r.assign('nb_theta', nb_theta)
            
            with st.spinner('Calculating GLM...'):
                ro.r('calc_gam()')
            
            # GAMの場合
            res_df = pd.read_csv(res_dir + f'/glm_{dist_short}_res.tsv', sep='\t', index_col=0)
            res_df = res_df.fillna(1)

            adj_pvalue_columns = [col for col in res_df.columns if '.adj.pvalue' in col]
            st.markdown("### Significant genes (FDR < 0.05):")
            total_significant = 0
            for col in adj_pvalue_columns:
                comparison = col.split('.adj.pvalue')[0]
                sig_count = (res_df[col] < 0.05).sum()
                total_significant += sig_count
                st.write(f"- {comparison}: {sig_count}")
            st.write(f"Total significant genes: {total_significant}")

            st.dataframe(res_df)
            
            st.session_state.deseq2 = res_df
            file_name = file_name_head + "_glm"
        
        # 結果のZIPアーカイブ作成とダウンロードボタン
        shutil.make_archive(file_name, format='zip', root_dir=res_dir)
        
        with open(file_name + ".zip", "rb") as fp:
            btn = st.download_button(
                label="Download Results",
                data=fp,
                file_name=file_name + ".zip",
                mime="zip"
            )
        
        try:
            os.remove(file_name + ".zip")
        except:
            pass



#　データを送る前にすべてゼロのデータは除くべき


# refが指定されているときはファイル名を調整する?
