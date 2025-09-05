import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import io
import sys
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()

st.set_page_config(page_title="PCA.", page_icon="√")
st.markdown("### PCA")

# プロットスタイルを設定する関数
def set_plot_style():
    plt.style.use('default')  # デフォルトスタイルをリセット
    plt.rcParams['figure.facecolor'] = 'white'  # 図の背景色を白に
    plt.rcParams['axes.facecolor'] = 'white'    # プロット領域の背景色を白に
#    sns.set_style("white")                      # seabornのスタイルを白背景に
    
# グラフ生成前に呼び出す
set_plot_style()

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

def normalize_totalreads(df):
    return 10**6 * df / df.sum()

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_xl = pd.read_csv(file, index_col = index_col, header = 0, sep = sep)
    return df_xl

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def pca_fit(df, random_state = 0):
    pca = PCA(random_state=0)
    emb = pca.fit_transform(df)
    return emb, pca



@st.cache_data
def mds_fit(df, n_components=2, random_state=0):
    embedding = MDS(n_components=n_components, random_state=random_state)
    emb = embedding.fit_transform(df)
    return emb


@st.cache_data
def calc_rlog_no_group():
    ro.r('calc_rlog_no_group()')

@st.cache_data
def calc_rlog():
    ro.r('calc_rlog()')

def remove_sample_num(i):
    i = re.sub('[_-]\d+$', '', i)
    i = re.sub('\d+$', '', i)
    return i

@st.cache_data
def calc_pca(df):
    from sklearn.decomposition import PCA
    z  = scaler.fit_transform(df.T)
    df_z = pd.DataFrame(z.T, columns=df.columns.to_list(), index=df.index.to_list())
    pca = PCA(random_state=0)
    x_embedded = pca.fit_transform(df_z.T)
    pca_col_name = ["PC" + str(i+1) for i in range(x_embedded.shape[1])]
    df_pca = pd.DataFrame(x_embedded[:,:], index= df.columns.to_list(), columns=pca_col_name)
    return df_pca, pca, pca_col_name


@st.cache_data
def calc_tsne(df, perplexity, n_components):
    from sklearn.manifold import TSNE
    if perplexity >= color_num:   # perplexityはサンプル数より多くないといけない
        perplexity = color_num - 1
        st.write('perplexity is set to ' + str(perplexity))
    embedding = TSNE(n_components=n_components, random_state=0, perplexity = perplexity)
    x_embedded = embedding.fit_transform(df.T)
    return x_embedded

@st.cache_data
def calc_umap(df, n_neighbors, min_dist, n_components):
    import umap.umap_ as UMAP
    if n_neighbors >= color_num:   # perplexityはサンプル数より多くないといけない
        n_neighbors = color_num - 1
        st.write('n_neighbors is set to ' + str(n_neighbors))
    embedding = UMAP.UMAP(n_components=n_components, random_state=0, n_neighbors=n_neighbors, min_dist=min_dist, spread=1.0)
    x_embedded = embedding.fit_transform(df.T)
    return x_embedded


@st.cache_data
def calc_mds(df, n_components=2):
    from sklearn.manifold import MDS
    mds = MDS(n_components=n_components,  random_state=0)
    x_embedded = mds.fit_transform(df.T)
    return x_embedded

st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #   st.write('file name kept')
    filename_add = ""

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

if "input_file_format" not in st.session_state:
    st.session_state.input_file_format = 'row = gene'
else:
    input_file_format = st.session_state.input_file_format
    

if use_upload == 'Yes':
    input_file_type = st.radio(
        "Data format:",
        ('tsv','csv', 'excel', 'Homer'), key='tsv')
    if input_file_type != 'Homer':
        input_file_format = st.radio(
            "Data structure:",
            ('row = gene', 'column = gene'))

        if input_file_format == 'row = gene':
            st.markdown("""
        row = genes (or other measurement items)
        |  | Sample1 | Sample2 |
        | --- | --- | --- |
        | Gene1 |  |  |
        | Gene2 | ||

        """)
        else:
            st.markdown("""
        column = gene (or other measurement items)
        |  | Gene1 | Gene2 |
        | --- | --- | --- |
        | Sample1 |  |  |
        | Sample2 | ||

        """)
        st.session_state.input_file_format = input_file_format

    st.write('')


    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        nonexplicit_excel = False
        try:
            if input_file_type == "csv":
                df = read_csv(uploaded_file, index_col = 0)
            elif input_file_type == "excel" :
                df = read_excel(uploaded_file, index_col = 0) #cacheされるfunctionを使う
            elif input_file_type ==  'Homer-excel':
                df = read_excel(uploaded_file)
            elif input_file_type == "Homer":
                df = read_csv(uploaded_file, sep = '\t')
            else:
                df = read_csv(uploaded_file, sep = '\t', index_col = 0)
        except: # excel fileを明示していないとき
            df = read_excel(uploaded_file)
            nonexplicit_excel = True

    ##### Homerの変換部分
        if input_file_type == "Homer":
            df = df.iloc[:,7:]
            colnames = df.columns.tolist()
            colnames[0] = 'Gene'

            # colnamesの変換
            search_word = '([^\ \(]*)\ \(.*'

            for i in range(1, len(colnames)):
                match = re.search(search_word, colnames[i])
                if match:
                    colnames[i] = match.group(1).replace(' ', '_')

            df.columns = colnames
            df['Gene'] = df['Gene'].astype(str)
            pattern = "([^|]*)"
            repatter = re.compile(pattern)
            f_annotation = lambda x: repatter.match(x).group(1)
            df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)

            df.set_index("Gene", inplace = True)

        else: # excelでは明示的にしない。
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
                df.set_index("Gene", inplace = True)
            elif nonexplicit_excel: # excel fileを明示していないとき
                df.set_index(content[0], inplace = True)


    ###########

        if input_file_format == 'column = gene':
            df = df.T

        file_name_head = os.path.splitext(uploaded_file.name)[0]


        st.session_state.df = df
        st.session_state.uploaded_file_name = file_name_head

    else:
        sys.exit(1)

if df is not None:

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
    df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く
############
    st.write(str(len(df)) + ' genes')
    st.write(df.iloc[:3,:30])
    convert_to = "raw_data"
    normalize_df = st.checkbox('Normalize count data?')
    if normalize_df:
        convert_to = st.radio("Convert to:",
            ('TMM','UQ','CTF', 'CUF', 'CPM', 'rlog', 'RPKM to TPM'), key='TMM')
        log_transform = None

        if convert_to != "raw_data":

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

            if (input_file_format == 'row = gene') and (convert_to != 'rlog'):
                df_conv = df_conv.T

 #           rlog_finish = False

            if convert_to == 'rlog':

                import rpy2.robjects as ro
                from rpy2.robjects.packages import importr
                from rpy2.robjects import pandas2ri
                from rpy2.robjects.vectors import StrVector
                import pyper
                r = pyper.R(use_pandas=True)
                f = ro.r("source('/home/cellxgene/streamlit/pages/deseq2_func.R')") # full pathが必要
                condition = [str(i) for i in df.columns.tolist()]
                group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
                group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字を除く
                df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])
                df = df.astype(float)
                df = df.round(0)
              #  df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く
                group = st.checkbox('Set groups for rlog?')
                if group:
                    edited_df_e = st.data_editor(df_e)
                    condition = edited_df_e.iloc[:,0].tolist()
                    st.write('Group: ' + '  '.join(condition))
                else:
                    condition = df.columns.tolist()
 #               if st.button('Run rlog calc'):
                r.assign('df',df)
                r("saveRDS(df, 'pyper_df.RDS')")

                ro.r("cts <- readRDS('pyper_df.RDS')")
                #まずベクターに変換
                r_condition =  ro.StrVector(condition)
                ro.r.assign('condition', r_condition)
#                batch = ["No batch"]
#                r_batch =  ro.StrVector(batch)
                ro.r.assign('condition', r_condition)
                ro.r.assign('temp_dir', 'temp')

                ro.r("make_coldata2()")
                if group:
                    calc_rlog()
                else:
                    calc_rlog_no_group()
                # df_conv = ro.conversion.rpy2py(rld) うまくいかない
                df_conv = pd.read_csv('temp/rld.tsv', sep = '\t', header = 0)
                content = df_conv.columns.tolist()
                content[0] = 'Gene'
                df_conv.columns = content
                df_conv.set_index("Gene", inplace = True)
              #  unlink("temp/rld.tsv")

                rlog_finish = True
                log_transform = None

            df = df_conv


    more_filt = st.checkbox('More filtering options')
    if more_filt:
        f_inf = -float('inf')
        p_inf = float('inf')

        min_val =  float(st.text_input("All values of each gene are larger than",  value=f_inf))
        max_val =  float(st.text_input("Max value is larger than",  value=f_inf))
        delta_val =  float(st.text_input("Delta (max - min value) >",  value=0))
        fold_val =  float(st.text_input("Fold (max / min) >",  value=1))
        min_variance =  float(st.text_input("Minmum variance across samples > (e.g., 0.3)",  value=0))
        high_min_val =  float(st.text_input("All values of each gene are smaller than or equal",  value=p_inf))
        high_max_val =  float(st.text_input("Min value is smaller than or equal",  value=p_inf))
        top_n =  float(st.text_input("Top n in mean",  value=p_inf))

        if min_val != f_inf:
            df = df[df.apply(min, axis=1) > min_val]


        if max_val != f_inf:
            df =  df[df.apply(max, axis=1) > max_val] #ここがminだと、0が一つでもあれば削除される。

        if delta_val > 1:
            df = df[df.apply(max, axis=1) > df.apply(min, axis=1) + delta_val]


        if fold_val > 1:
            df = df[df.apply(max, axis=1) > df.apply(min, axis=1) * fold_val]

        if min_variance > 0:
            df = df.loc[(df.index[(df.var(axis=1) > min_variance)]),:]

        if high_min_val != p_inf:
            df = df[df.apply(max, axis=1) <= high_min_val]

        if high_max_val != p_inf:
            df =  df[df.apply(min, axis=1) <= high_max_val] #ここがminだと、0が一つでもあれば削除される。

        if top_n != p_inf:
            top_ix = df.mean(axis = 1).sort_values(ascending=False).head(10).index
            new_index = [x for x in df.index.to_list() if x in top_ix]
            df = df.loc[new_index,:]

        df = df.astype('float')
    log_df = False
    if convert_to != 'rlog':
        log_df = st.checkbox('Log transforamtion?')
        if log_df:
            log_transform = st.radio(
            "Log transformation:",
            ('None', 'asinh','log2(x+1)', 'ln(x+1)', 'log10(x+1)'), key='None')

    condition = [str(i) for i in df.columns.tolist()] #error防止
    label2show = df.columns.to_list()
    category = st.checkbox('Set color groups and/or change labels??', value=True)
    st.write("###### Uncheck this to set individual colors for each point.")
    if category:
        condition = [str(i) for i in condition] # https://stackoverflow.com/questions/69578431/how-to-fix-streamlitapiexception-expected-bytes-got-a-int-object-conver
        # エラーが出るためすべてリストの内容をstrに変換してからdfを作る
#        df_f = pd.DataFrame(df.columns.tolist(), index = df.columns.tolist() , columns = ["Group"]).copy()
        group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
        group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字を除く
        df_f = pd.DataFrame(group_condition, index = condition , columns = ["Group"])
        df_f["Label"] = condition
        with st.form("input_groups and Label"):
            edited_df_f = st.data_editor(df_f)
            submitted = st.form_submit_button("Submit")

        condition = edited_df_f.iloc[:,0].tolist()
        label2show = edited_df_f.iloc[:,1].tolist()
        st.write('Group: ' + '  '.join(condition))
        st.write('Label: ' + '  '.join(label2show))



    log_transform_word = ''
    if log_df:
        if log_transform  != "None":
            if log_transform == 'asinh':
                df = np.arcsinh(df)
                log_transform_word = ".asinh"
            if log_transform == 'log2(x+1)':
                df = np.log2(df+1)
                log_transform_word = ".log2"
            if log_transform == 'log10(x+1)':
                df = np.log10(df+1)
                log_transform_word = ".log10"
            if log_transform == 'ln(x+1)':
                df = np.log1p(df)
                log_transform_word = ".ln"



    subset_gene = st.checkbox('Use subset of genes (rows)?')
    if subset_gene:
        st.markdown("##### Genes (comma, space, CR separated):")
        genes = st.text_input("genes",label_visibility = 'collapsed')
        gene_list = []
        if len(genes) > 0:
            if ',' not in genes:
                gene_list = genes.split(' ')
            else:
                genes =  ''.join(genes.split()) #空白を除く
                gene_list = genes.split(',')
            genes = list(filter(lambda x: x != "", genes)) #空白を除く

            # intersectを取る
            gene_list = set(df.index.tolist()) & set(gene_list)
            df = df.loc[list(gene_list),:]
    st.markdown('---')

    with st.sidebar:
        reduc  = st.radio(
        "Reduction method:",
        ('PCA', 'tSNE','UMAP', "MDS"), key='PCA')
        reduc_dim = st.radio(
        "Dimension:",
        ('2D', '3D'), key='2D')
        if reduc_dim ==  '3D':
            n_components = 3
            dot_size = 3
        else:
            n_components = 2
            dot_size = 8
        if reduc == 'PCA':
            x_show = int(st.text_input("X dim:", value = 1))
            y_show = int(st.text_input("Y dim:", value = 2))
            if x_show > 15 or y_show > 15:
                st.write("Max dimension is 15")
            else:
                pca_x = 'PC' + str(x_show)
                pca_y = 'PC' + str(y_show)
            if reduc_dim == '3D':
                z_show = int(st.text_input("Z dim:", value = 3))
                if z_show > 15 :
                    st.write("Max dimension is 15")
                else:
                    pca_z = 'PC' + str(z_show)

        if reduc == 'tSNE' or reduc == 'UMAP':
            calc_z = st.checkbox('Convert to Z score?')
            if calc_z:
                z  = scaler.fit_transform(df.T)
                df = pd.DataFrame(z.T, columns=df.columns.to_list(), index=df.index.to_list())
        if reduc == 'tSNE':
            perplexity = float(st.text_input("Perplexity:", value = 30))
            st.markdown('See https://distill.pub/2016/misread-tsne/ and https://qiita.com/maskot1977/items/2213e33c31cfc5403bf6')
            st.write('perplexityは局所と全体的な特性のどちらを保存するかのバランスを決める。どれだけ近傍の点を考慮するかの指標。大きなデータセットには大きな値を用いる。典型的には5-50だが調整が必要。')
        if reduc == "UMAP":
            n_neighbors = st.number_input("n_neighbors:", min_value = 2, value = 15) # integer
            min_dist = st.number_input("min_dist:", min_value= 0.0, max_value=1.0, value = 0.1, step = 0.05)
            st.write('n_neighborsは各データポイントを埋め込む際に考慮される近隣の点の数。数値が大きいと全体的構造が強調され、小さいと局所構造が保存される。典型的には2-100。')
            st.write('次元圧縮後の点間の最短距離を示す。小さいと点が密集し、大きいと点が広がりトポロジカルな構造を保存する。')


#        cmap  = st.radio(
#        "Color map:",
#        ('None', 'Set1','Set2','Set3','tab10','tab20', 'random'), key='None')

        st.markdown('##### Plot attributes')
        width = float(st.text_input("Plot x size:", value = 600))
        height = float(st.text_input("Plot y size:", value = 500))
        marker_size = float(st.text_input("Dot size:", value = dot_size))
        show_text = st.checkbox('Show labels?')
        if show_text:
            textfont_size = float(st.text_input("Text font size:", value = 14))
        change_color = st.checkbox('Change color sequence?')
        color_choice = 'Plotly'
        if change_color:
            show_color = st.checkbox('Show color sequence?')
            if show_color:
                color_fig = plt.figure()
                color_fig = px.colors.qualitative.swatches()
                st.plotly_chart(color_fig)
            color_choice = st.selectbox('Color sequence:',
                ('Plotly','D3','G10','T10','Alphabet','Dark24','Light24','Set1','Pastel1','Dark2','Set2',
                'Pastel2','Set3','Antique','Bold','Pastel','Prism','Safe','Vivid'), key = 'Plotly')
        c_choice = getattr(px.colors.qualitative, color_choice)
        saveas  = st.radio("Save plot as:",('pdf', 'png','html'), key='PDF')


    st.write('Data to use:')
    st.write(df.iloc[:3,:30])

    st.write(str(len(df)) + ' genes')

    color_num = len(df.columns.to_list())
    # Create an in-memory buffer
    buffer = io.BytesIO()
    if st.button("Draw", type="primary"):
        if reduc == 'PCA':
            df_pca, pca, pca_col_name = calc_pca(df)

            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_pca, x= pca_x, y= pca_y, color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_pca, x= pca_x, y= pca_y, color=condition, text = label2show , color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
                st.write("Explained variance ratio:")
                st.write(pca_x + ": " + str(pca.explained_variance_ratio_[x_show - 1]))
                st.write(pca_y + ": " + str(pca.explained_variance_ratio_[y_show - 1]))
            else:
                if not show_text:
                    fig = px.scatter_3d(df_pca, x= pca_x, y= pca_y, z= pca_z,  color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_pca, x= pca_x, y= pca_y, z= pca_z, color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
                st.write(pca_x + ": " + str(pca.explained_variance_ratio_[x_show - 1]))
                st.write(pca_y + ": " + str(pca.explained_variance_ratio_[y_show - 1]))
                st.write(pca_z + ": " + str(pca.explained_variance_ratio_[z_show - 1]))
            st.plotly_chart(fig)

        elif reduc == 'tSNE':
            x_embedded = calc_tsne(df, perplexity,n_components)
            if reduc_dim == "2D":
                df_tsne = pd.DataFrame(x_embedded[:,:2], index= df.columns.to_list(), columns=["tSNE1","tSNE2"])
            else:
                df_tsne = pd.DataFrame(x_embedded[:,:3], index= df.columns.to_list(), columns=["tSNE1","tSNE2",'tSNE3'])

  #          if category:
  #              df_tsne.index = condition

            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_tsne, x= 'tSNE1', y= 'tSNE2', color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_tsne, x= 'tSNE1', y= 'tSNE2', color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            else:
                if not show_text:
                    fig = px.scatter_3d(df_tsne, x= 'tSNE1', y= 'tSNE2',z='tSNE3',  color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_tsne, x= 'tSNE1', y= 'tSNE2',z='tSNE3',  color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            st.plotly_chart(fig)

        elif reduc == 'UMAP':
            x_embedded = calc_umap(df, n_neighbors, min_dist, n_components)
#            import umap.umap_ as UMAP
#            if n_neighbors >= color_num:   # perplexityはサンプル数より多くないといけない
#                n_neighbors = color_num - 1
#                st.write('n_neighbors is set to ' + str(n_neighbors))
#            x_embedded = umap_fit(df.T, n_components=n_components, random_state=0, n_neighbors=n_neighbors, min_dist = min_dist )

            if reduc_dim == "2D":
                df_umap = pd.DataFrame(x_embedded[:,:2], index= df.columns.to_list(), columns=["UMAP1","UMAP2"])
            else:
                df_umap = pd.DataFrame(x_embedded[:,:3], index= df.columns.to_list(), columns=["UMAP1","UMAP2",'UMAP3'])

#            if category:
#                df_umap.index = condition
            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_umap, x= 'UMAP1', y= 'UMAP2', color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_umap, x= 'UMAP1', y= 'UMAP2', color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            else:
                if not show_text:
                    fig = px.scatter_3d(df_umap, x= 'UMAP1', y= 'UMAP2',z='UMAP3', color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_umap, x= 'UMAP1', y= 'UMAP2',z='UMAP3', color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)

            st.plotly_chart(fig)


        if reduc == 'MDS':
            x_embedded = calc_mds(df, n_components)
            if reduc_dim == "2D":
                df_umap = pd.DataFrame(x_embedded[:,:2], index= df.columns.to_list(), columns=["MDS1","MDS2"])
            else:
                df_umap = pd.DataFrame(x_embedded[:,:3], index= df.columns.to_list(), columns=["MDS1","MDS2",'MDS3'])

            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_umap, x= 'MDS1', y= 'MDS2', color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_umap, x= 'MDS1', y= 'MDS2', color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            else:
                if not show_text:
                    fig = px.scatter_3d(df_umap, x= 'MDS1', y= 'MDS2',z='MDS3', color=condition, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_umap, x= 'MDS1', y= 'MDS2',z='MDS3', color=condition, text = label2show, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)

            st.plotly_chart(fig)

        if saveas == 'html':
            #fig.write_html(file=buffer )
            fig.write_html('dummy.html' )
        else:
            fig.write_image(file=buffer, format=saveas, engine="kaleido" )
        file_name = file_name_head + '.' + reduc + '.' + saveas

        if saveas != 'html':
            st.download_button(label="Download plot",
                            data=buffer,
                            file_name=file_name,
                            mime='application/octet-stream')
        else:
            with open("dummy.html", "rb") as file:
                st.download_button(label="Download plot",
                            data=file,
                            file_name=file_name,
                            mime='text/html')

        if reduc == 'PCA':
            loadings = pca.components_.T * np.sqrt(pca.explained_variance_)
            loading_matrix = pd.DataFrame(loadings, columns=pca_col_name, index=df.index)
            loading_tsv = convert_df(loading_matrix)
            loading_file_name = file_name_head + '.LoadingsMatrix.tsv'
            st.download_button(label="Download loadings matrix",
                                data=loading_tsv,
                                file_name=loading_file_name,
                                mime='test/csv')


