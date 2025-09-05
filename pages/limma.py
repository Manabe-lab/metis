#!!!!!!!!!!!!!! pip install rpy2==3.5.1  新しいバージョンはエラーが出る

# 基本的にglobal変数で計算する。
# pythonからassgnされるのはglobal変数


import streamlit as st
import csv
import re
import os
import numpy as np
import pandas as pd
import shutil
from PIL import Image
from helper_func import clear_old_directories, clear_old_files, make_r_compatible_column_names, remove_after_space, remove_sample_num
import time
import sys
from rpy2 import robjects
import io
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from itertools import combinations
from rpy2.robjects import StrVector
from rpy2.robjects.lib import grdevices
import matplotlib.pyplot as plt
import pdf2image

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

def format_comparison(comparison_str):
    # group_factor を含む部分を除去
    comparison_str = re.sub(r'([-]?\d*\*)?group_factor\w+\s*', '', comparison_str)
    
    # 係数の文字列を解析
    coefficients = comparison_str.split()
    pos_group = None
    neg_group = None
    
    for coef in coefficients:
        if coef.startswith('1*'):
            pos_group = coef[2:]  # '1*' を除去
        elif coef.startswith('-1*'):
            neg_group = coef[3:]  # '-1*' を除去
    
    if pos_group and neg_group:
        return f"{pos_group} vs. {neg_group}"
    elif pos_group:
        return f"{pos_group} vs. Control"
    elif neg_group:
        return f"Control vs. {neg_group}"
    else:
        return "Comparison information not available"  # パースできない場合のメッセージ


def create_integer_contrasts(groups):
    # グループの数を取得
    n = len(groups)
    
    # すべての可能な2グループの組み合わせを生成
    group_pairs = list(combinations(groups, 2))
    
    # コントラストを格納するリスト
    contrasts = []
    
    for pair in group_pairs:
        contrast = np.zeros(n, dtype=int)  # dtype=intを追加して整数型に
        first_index = groups.index(pair[0])
        second_index = groups.index(pair[1])
        
        # 先のグループを-1、後のグループを+1に設定
        contrast[first_index] = -1
        contrast[second_index] = 1
        
        contrasts.append(contrast)
    
    # コントラスト行列を作成
    contrast_matrix = np.array(contrasts, dtype=int)  # dtype=intを追加
    
    # DataFrameに変換
    contrast_df = pd.DataFrame(contrast_matrix, 
                               columns=groups,
                               index=[f"{pair[1]}_vs_{pair[0]}" for pair in group_pairs])
    
    return contrast_df

def capture_r_output_as_dataframe(r_code):
    # R のデータフレームを Python の pandas DataFrame に変換する設定
    pandas2ri.activate()
    
    # R コードを実行し、結果を取得
    result = ro.r(r_code)
    
    # 結果が data.frame の場合、pandas DataFrame に変換
    if isinstance(result, ro.vectors.DataFrame):
        df = pandas2ri.rpy2py(result)
        return df
    else:
        # データフレームでない場合は文字列として返す
        return str(result)

def capture_r_output(r_code):
    # 標準出力をキャプチャするための設定
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()

    try:
        # Rコードを実行
        ro.r(r_code)
        # キャプチャした出力を取得
        output = sys.stdout.getvalue()
    finally:
        # 標準出力を元に戻す
        sys.stdout = old_stdout

    return output


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

st.set_page_config(page_title="Calculate limma-voom", page_icon="📃")


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


st.sidebar.title("Options")
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


st.markdown("### raw count dataを使う")

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

#uploaded_size_factors = None
#if "use_custom_size_factors" not in st.session_state:
#    st.session_state.use_custom_size_factors = False

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('auto', 'Homer','tsv','csv','excel'), index = 0, label_visibility = 'collapsed')
    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])
    use_sf = st.checkbox('Upload Size Factors (Optional)')
#    if use_sf:
#        uploaded_size_factors = st.file_uploader("Choose a size factors file (TSV format)", type=['tsv'])

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

        #特殊文字を除く
        df = make_r_compatible_column_names(df)

    else:
        sys.exit(1)

if df is not None:
    st.write('Original gene number:  ' + str(len(df)))
    st.write(df.head())

    # floatに変換 誤射悟入
    df = df.astype(float)

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
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene name of the second occurrence has _2 at the end.")
        lis = df.index.values
        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]

    #特殊文字を除く
    df = make_r_compatible_column_names(df)

    st.write(df.head())
    total_count = pd.DataFrame(df.sum()[1:])
    total_count.columns= ['Total counts']
    large_var = False
    if max(total_count['Total counts']) > min(total_count['Total counts']) * 2:
        large_var = True
        st.markdown("### Large difference (>2x) in counts")
        import matplotlib.pyplot as plt
        import seaborn as sns
        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']


        f1 = calc_barplot(df_sum.T, ylabel = "Total counts")
        st.pyplot(f1)

        f2 = calc_barplot(np.log1p(df), ylabel = "ln(x+1)")
        st.pyplot(f2)


    with st.sidebar:
        st.markdown("##### Filter out weakly-expressed genes:")
        independentFiltering = st.checkbox('Yes', value= True)
        st.markdown("""低発現遺伝子の除外はFDRの計算を改善する。
        filtered outされた遺伝子が多い場合、GSEA用に全遺伝子をランキングするためにはチェックを外す。""")
        min_threshold = 0
        max_threshold = 0
        if independentFiltering:

            st.markdown("##### Filter the genes > counts in all samples:")
            min_threshold = st.number_input("count minimum", value = 0, label_visibility = 'collapsed')
            min_threshold = int(min_threshold)
            st.markdown("##### Filter the genes > counts in at least one sample:")
            max_threshold = st.number_input("count max", value = 0, label_visibility = 'collapsed')
            max_threshold = int(max_threshold)

            sample_threshold = 0

        if large_var:
            st.markdown("##### Filter the samples <= counts:")
            sample_threshold = st.number_input("Minimum total cout", value = 0, label_visibility = 'collapsed')

        which_voom = st.radio("Voom method: ", ['voom', 'voomWithQualityWeights', 'voomLmFit'] , index=2)
        st.write("voomWithQualityWeights for a dataset containing low quality data")
        st.write("voomLmFit for sparse data (e.g., scRNA)")


    if min_threshold > 0:
        df = df[df.apply(min, axis=1) > min_threshold]
    if max_threshold > 0:
        df = df[df.apply(max, axis=1) > max_threshold]

    st.write('Filtered gene number:  ' + str(len(df)))

    if any(df.sum() <= sample_threshold): # count 0の列を除く
        st.markdown('#### There are the samples that have counts <= ' + str(sample_threshold))
        st.write(", ".join(df.columns[df.sum() <= sample_threshold].to_list()))
        st.write('They are removed. Now data are:')
        df = df.drop(df.columns[df.sum() <= sample_threshold].to_list(), axis = 1)
        st.write(df.head())

    #特殊文字を除く
    df.columns = df.columns.str.replace('[^A-Za-z0-9]+', '_')
    df.columns = df.columns.str.replace('-', '_')

    condition = [str(i) for i in df.columns.tolist()] #error防止
    group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
#    group_condition = [remove_after_space(x) for x in condition] #スペース以降を除く
    group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字と_を除く


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

#    ref_in = st.checkbox('Setting referece group?')

    for i in df.columns.values:
        a = df.select_dtypes(exclude='number')
        if not a.empty:
            st.write("There is a non-nemeric value in ")
            st.write(a)

#    df = excel_autoconversion(df) # 1-Marなどの誤変換への対応


    st.markdown("---")



    if st.button('Run limma-voom'):
        with st.spinner('Calculating limma-voom...'):
            df_path = temp_dir + "/df.tsv"
            
            r_condition = StrVector(condition)
            ro.r.assign('condition', r_condition)

            df.to_csv(df_path, sep = '\t')
            ro.r(f'''
            library(limma)
            library(edgeR)
            
            # Read data
            rawdata <- read.csv('{df_path}', sep = '\t', row.names = 1)
            
            # Create DGEList object
            y <- DGEList(counts=rawdata)
            
            # Normalize library sizes
            y <- calcNormFactors(y)
            
            # Create design matrix
            group_factor <- factor(condition)
            design <- model.matrix(~0 + group_factor)
            colnames(design) <- levels(group_factor)
            ''')


            # 選択された方法に基づいて異なるRコードを実行
            if which_voom != 'voomLmFit':
                ro.r(f'''
                # voom変換と品質重みの適用
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                v <- {which_voom}(y, design, plot=TRUE)
                dev.off()
                
                # 線形モデルのフィッティング
                fit <- lmFit(v, design)
                ''')
            else:  # edgeR::voomLmFit
                ro.r(f'''
                # voom変換と線形モデルフィッティングを一度に行う
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                fit <- voomLmFit(y, design, plot=TRUE)
                dev.off()
                ''')

            # Display the voom plot
            try:
                with open(ro.r('png_file')[0], "rb") as file:
                    st.image(file.read(), caption='Voom Plot', use_column_width=True)
            except Exception as e:
                st.error(f"Failed to display Voom plot: {str(e)}")

            contrast_df = create_integer_contrasts(list(dict.fromkeys(condition)))

            res = dict()
            for i in range(len(contrast_df)):
                l = map(str, contrast_df.iloc[i, :].to_list())
                c = ', '.join(l)
                ro.r(f'''
                contrast <- makeContrasts(c({c}), levels=design)
                fit2 <- contrasts.fit(fit, contrast)
                fit2 <- eBayes(fit2)
                top_table <- topTable(fit2, sort.by="P", n=Inf)
                ''')

                # Convert R data frame to pandas DataFrame
                with localconverter(ro.default_converter + pandas2ri.converter):
                    top_table_df = ro.conversion.rpy2py(ro.r('top_table'))
                
                top_table_df = top_table_df.reset_index()
                top_table_df = top_table_df.rename(columns={'index': 'gene'})

                comparison_str = f"{contrast_df.index[i]}"
                st.write(comparison_str)
                st.write(top_table_df.head())

                # Generate and display MA plot
                ro.r(f'''
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                limma::plotMA(fit2, main=paste("MA Plot -", "{comparison_str}"))
                dev.off()
                ''')
                
                try:
                    with open(ro.r('png_file')[0], "rb") as file:
                        st.image(file.read(), caption=f'MA Plot - {comparison_str}', use_column_width=True)
                except Exception as e:
                    st.error(f"Failed to display MA plot for {comparison_str}: {str(e)}")

                res[comparison_str] = top_table_df

            # Merge all results
            new_dfs = []
            for key, df in res.items():
                # Rename columns to include the comparison name, but keep 'gene' as is
                df = df.rename(columns={col: f"{key}.{col}" for col in df.columns if col != 'gene'})
                new_dfs.append(df)

            # Use 'gene' column as the key for merging
            merged_df = new_dfs[0]
            for df in new_dfs[1:]:
                merged_df = pd.merge(merged_df, df, on='gene', how='outer')

            # Set 'gene' as the index after merging
            merged_df = merged_df.set_index('gene')

            st.write(merged_df)

            # Convert DataFrame to TSV
            tsv = merged_df.to_csv(index=True, sep='\t')

            # Create download button
            st.download_button(
                label="Download data as TSV",
                data=tsv,
                file_name=file_name_head + ".limma-voom.tsv",
                mime="text/tab-separated-values"
            )