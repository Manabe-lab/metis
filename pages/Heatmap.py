import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import shutil
import os
import sys
import matplotlib.colors as mcolors
from helper_func import mk_temp_dir
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage as scipy_linkage
from sklearn.decomposition import NMF

# フォント設定を追加
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42  # TrueTypeフォントを使用
plt.rcParams['ps.fonttype'] = 42   # TrueTypeフォントを使用


st.set_page_config(page_title="Heatmap", page_icon="🌡")

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t', header = None).encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header= 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
    return df_c

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #   st.write('file name kept')
    filename_add = ""


#https://discuss.streamlit.io/t/dynamically-created-multiple-checkbox/18273/2
def checkbox_container(data):
#    st.header('Select columns or rows')
#    new_data = st.text_input('Enter country Code to add')
    cols = st.columns(5)
#    if cols[0].button('Add Coutry'):
#        dummy_data.append(new_data)
    if cols[0].button('Select All'):
        for i in data:
            st.session_state['dynamic_checkbox_' + i] = True
        st.rerun()
    if cols[1].button('UnSelect All'):
        for i in data:
            st.session_state['dynamic_checkbox_' + i] = False
        st.rerun()
    for i in data:
        st.checkbox(i, key='dynamic_checkbox_' + i)

def get_selected_checkboxes():
    return [i.replace('dynamic_checkbox_','') for i in st.session_state.keys() if i.startswith('dynamic_checkbox_') and st.session_state[i]]


def create_custom_cmap(colors, name='custom'):
    if len(colors) == 2:
        # For two colors, create a linear gradient
        return mcolors.LinearSegmentedColormap.from_list(name, colors, N=256)
    elif len(colors) == 3:
        # For three colors, create two segments
        cmap1 = mcolors.LinearSegmentedColormap.from_list("cmap1", colors[:2], N=128)
        cmap2 = mcolors.LinearSegmentedColormap.from_list("cmap2", colors[1:], N=128)
        # Combine the two segments
        newcolors = np.vstack((cmap1(np.linspace(0, 1, 128)),
                               cmap2(np.linspace(0, 1, 128))))
        return mcolors.ListedColormap(newcolors, name=name)
    else:
        raise ValueError("Only 2 or 3 colors are supported for custom colormap creation.")

# Matplotlibのカラーマップをカテゴリごとに分類
colormap_categories = {
    'Sequential': [
        'viridis', 'plasma', 'inferno', 'magma', 'cividis',
        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
    ],
    'Sequential (2)': [
        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
        'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
        'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
    ],
    'Diverging': [
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
        'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
    ],
    'Cyclic': ['twilight', 'twilight_shifted', 'hsv'],
    'Qualitative': [
        'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
        'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c'
    ],
    'Miscellaneous': [
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
        'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral',
        'gist_ncar'
    ]
}

def show_colormap(cmap_name):
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    fig, ax = plt.subplots(figsize=(10, 1))
    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap_name))
    ax.set_axis_off()
    st.pyplot(fig)


@st.cache_data
def perform_clustering_computation(data, method, metric):
    st.write(method)
    st.write(metric)
    # データの正規化や前処理をここで行う（必要な場合）
    
    # 行のクラスタリング
    row_distances = pdist(data, metric)
    row_linkage = scipy_linkage(row_distances, method)
    
    # 列のクラスタリング
    col_distances = pdist(data.T, metric)
    col_linkage = scipy_linkage(col_distances, method)
    
    return row_linkage, col_linkage, row_distances, col_distances

def plot_clustermap(df, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c, xticklabels, yticklabels, annot, fmt, linewidths, linecolor, py_x_size, py_y_size, x_font_size, y_font_size):
    g = sns.clustermap(
        df,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        center=v_center,
        cmap=cmap,
        vmin=v_min,
        vmax=v_max,
        row_cluster=y_c,
        col_cluster=x_c,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        annot=annot,
        fmt=fmt,
        linewidths=linewidths,
        linecolor=linecolor,
        figsize=(py_x_size, py_y_size)
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=x_font_size)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
    return g

def clean_column_name(col):
    if isinstance(col, (int, float)):
        # 整数の場合はそのまま文字列に変換
        if col.is_integer():
            return str(int(col))
        # 小数の場合は、末尾の0を削除
        return f'{col:g}'
    elif isinstance(col, str):
        # 文字列の場合、数値に変換できるなら変換して処理
        try:
            num = float(col)
            if num.is_integer():
                return str(int(num))
            return f'{num:g}'
        except ValueError:
            # 数値に変換できない場合はそのまま返す
            return col
    else:
        return str(col)

st.markdown("### Heatmap")
st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")
# temp内に保存する
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir, res_dir = mk_temp_dir("Heatmap")
    st.session_state.temp_dir = temp_dir
else:
    temp_dir = st.session_state.temp_dir
    temp_dir, res_dir = mk_temp_dir("Heatmap", temp_dir)


use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Data in memory:")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
 #       st.write(st.session_state.uploaded_file_name)
        file_name_head = st.session_state.uploaded_file_name
        # Homer対応
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list(): # Row_nameを含むとき
            df = df.set_index('Row_name')
            df.index.name = "Gene"


if use_upload == 'Yes':
    input_file_type = st.radio(
        "Data format:",
        ('tsv','csv', 'excel'))
    #Gene_column = 0
    genome = st.checkbox("Genome occupancy data (eg., Homer's annotatePeaks.pl output)?")
    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])

    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, header = None, index_col = None)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', header=None, index_col = None)
        else:
            df = read_excel(uploaded_file, index_col = None, header = None)
        st.write("Uploaded data:")
        st.write(df.head(3))
        st.write('Data Dimension: '+str(df.shape))

        st.markdown("###### Data format should be genes as rows and sample as columns.")
        st.markdown("""
    |  | Sample1 | Sample2 |
    | --- | --- | --- |
    | Gene1 |  |  |
    | Gene2 | ||

    """)
        st.write("    ")

        transpose_df = st.checkbox('Transpose the data?')

        if transpose_df:
            df = df.T
        if not genome:
            df.columns = df.iloc[0,:].tolist()   # transposeすると狂うので、transposeした後にcolumnsを決める
        else: # homer annotate に対して登場順で頭に数字をつける
            org_col = df.iloc[0,:].tolist()
            from collections import defaultdict
            # カラム名の出現回数を追跡するための辞書
            column_counts = defaultdict(int)
            # 新しいカラム名を生成
            new_columns = [org_col[0]]  # 最初の列名はそのまま
            for col in org_col[1:]:
                col = clean_column_name(col)
                #col = str(col).rstrip('.0')
                column_counts[col] += 1
                if column_counts[col] > 1:
                    new_col_name = f"{column_counts[col]}_{col}"
                else:
                    new_col_name = col
                new_columns.append(new_col_name)

            # 新しいカラム名を設定
            df.columns = new_columns

        df = df.drop(0, axis = 0) # 1行目を列名にして除く
        content = df.columns.tolist()
        # Rの(0,0)が空白のデータへの対応==============================
        if isinstance(content[0], float) and np.isnan(content[0]):
            st.write("0,0isnan")
            content[0] = "NaN"
            df.columns = content
        #========================================================

        Gene_column = content[0]
        if "Annotation/Divergence" in content:
              # colnamesの変換
            search_word = '([^\ \(]+).*'

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
    #        df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
            # annotation/divergence以前を除く
            df = df.loc[:,'Annotation/Divergence':]
            content = df.columns.tolist()
            content[0] = 'Gene'
            df.columns = content
            Gene_column = "Gene"
            st.write("Converted Annotation/Divergence to gene symbols.")

        elif "Gene" in content:
            Gene_column =  "Gene"
        else:
            Gene_column =  st.selectbox('Select gene name column', content)

        df = df.set_index(Gene_column)
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        st.session_state.uploaded_file_name = file_name_head
        st.session_state.df = df

    else:
        sys.exit(1)
   ##### file upload

if df is not None:
    nonzero = st.checkbox('Remove all zero genes?', value=True)
    if nonzero:
        df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く
    df = df.dropna(how='any', axis=0)
    df = df.astype(float) #文字列が残る

    st.write(df.head(3))

    st.markdown('---')
    st.markdown("##### Filter and transform data?")
    calc_z = False
    center0_z = False  # Z-scoreのときはTrueにする
    howlog = 'No'
    Manip = st.checkbox('minip', label_visibility = 'collapsed')
    if Manip:
        f_inf = -float('inf')
        p_inf = float('inf')
        min_val = f_inf
        max_val = f_inf
        high_min_val = p_inf
        high_max_val = p_inf
        delta_val = 1
        fold_val = 1
        min_variance = 0
        top_n = p_inf
        more_filt = st.checkbox('Additional filtering options (e.g., minimal values, FC...)')
        if more_filt:
            min_val =  float(st.text_input("All values of each gene are larger than",  value=f_inf))
            max_val =  float(st.text_input("Max value is larger than",  value=f_inf))
            delta_val =  float(st.text_input("Delta (max - min value) >",  value=0))
            fold_val =  float(st.text_input("Fold (max / min) >",  value=1))
            min_variance =  float(st.text_input("Minmum variance across samples > (e.g., 0.3)",  value=0))
            high_min_val =  float(st.text_input("All values of each gene are smaller than or equal",  value=p_inf))
            high_max_val =  float(st.text_input("Min value is smaller than or equal",  value=p_inf))
            top_n =  float(st.text_input("Top n in mean",  value=p_inf))

        st.markdown("######   ")
        calc_div = st.checkbox('Divided by (e.g., 1000, 1000,000)?', value = False)
        if calc_div:
            div_unit =  int(st.text_input("Unit: ",  value=1))
            df = df/div_unit
        calc_log = st.checkbox('Log transformation?')
        if calc_log:
            howlog = st.radio('Method', ['log2+1', 'log2', 'loge+1', 'loge','log10+1','log10', 'asinh'])
        else:
            howlog ='No'

        st.markdown("######   ")

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

        if howlog == 'log2+1':
            df = np.log2(df+1)
        elif howlog == 'log2':
            df = np.log2(df)
        elif howlog == 'loge+1':
            df = np.log1p(df)
        elif howlog == 'loge':
            df = np.log(df)
        elif howlog == 'log10+1':
            df = np.log10(df+1)
        elif howlog == 'log10':
            df = np.log10(df)
        elif howlog == 'asinh':
            df = np.arcsinh(df)



    st.markdown("##### Data standardization (Z-score transformation) ?")
    calc_z = st.checkbox('Z-score?', label_visibility =  'collapsed')
    if calc_z:
        center0_z= True
        df_z = df.copy()
        m = df_z.mean(1)
        s = df_z.std(1)
        df_z = df_z.sub(m, axis=0).div(s, axis = 0)
        df_z = np.round(df_z, decimals=10)
        df_z = df_z.loc[~(df_z==0).all(axis=1)] #すべて0のrowを除く
        df_z = df_z.dropna(how='any', axis=0) #エラー対応
        df = df_z


    st.markdown('---')


    st.markdown("##### Use subset of genes (or rows)?")
    subset_gene = st.checkbox('Use subset of genes (rows)?',label_visibility =  'collapsed')
    if subset_gene:
        st.markdown("##### Genes (comma, semicolon, space, CR separated):")
        genes = st.text_input("genes", label_visibility='collapsed')
        keep_all = st.checkbox('Do not remove duplicated genes?')
        gene_list = []
        if len(genes) > 0:
            genes = genes.replace("'", "")
            genes = genes.replace('"', "")
            gene_list = genes.split(' ')  # まず空白で分離
            gene_list = list(filter(lambda a: a != '', gene_list))  # 空白のみを除く
            
            # 各区切り文字で分割して順序を保持
            if ',' in genes:
                gene_list = sum([x.split(',') for x in gene_list], [])
            if ';' in genes:
                gene_list = sum([x.split(';') for x in gene_list], [])
            if '\t' in genes:
                gene_list = sum([x.split('\t') for x in gene_list], [])
            if '\n' in genes:
                gene_list = sum([x.split('\n') for x in gene_list], [])
                
            # 重複を削除しながら順序を保持
            seen = set()
            ordered_unique_genes = []
            for gene in gene_list:
                gene_lower = gene.lower()
                if gene_lower not in seen:
                    seen.add(gene_lower)
                    ordered_unique_genes.append(gene)

            if keep_all:
                ordered_unique_genes = gene_list
            
            # データフレームのインデックスとマッチング（順序を保持）
            df_index_set = set(df.index.tolist())
            gene_subset = []
            for gene in ordered_unique_genes:
                if any(x.lower() == gene.lower() for x in df_index_set):
                    matching_gene = next(x for x in df_index_set if x.lower() == gene.lower())
                    gene_subset.append(matching_gene)
                    
            df = df.loc[gene_subset, :]

    st.markdown('---')

    df = df.astype(float)
    st.markdown("##### Cleaned data:")
 #   df.iloc[:3,:]
    st.write(df.head())
    st.write('Data Dimension: '+str(df.shape))
    st.markdown('---')

    with st.sidebar:
        only_minmax = st.checkbox('**Two colors?**')
        if only_minmax:#ここで簡単に2色に変更できるようにする
            v_center = None
            c_cmap = 'viridis'
            cmap = 'viridis'
        else:
            c_cmap = 'bwr'
            cmap = 'bwr'

        v_min = None
        v_max = None
        if calc_z:
            v_center = 0
        else:
            v_center = None


        annot  = False
        fmt = None
        linewidths= 0
        change_c = st.checkbox('**Change color map?**')
        if change_c:
            create_c = st.checkbox('Create custom color map?')
            if not create_c:
                category = st.selectbox("Select colormap category", options=list(colormap_categories.keys()))
                c_cmap = st.selectbox("Select colormap", options=colormap_categories[category])
                inv_c = st.checkbox('Inverse color map?')
                if inv_c:
                    c_cmap = c_cmap + "_r"
                show_colormap(c_cmap)

#            from matplotlib.colors import ListedColormap
            if create_c:
            #    two_c = st.checkbox('Make two-color map?')
                if only_minmax:
                    min_c = st.color_picker('Min color:', '#ffffff')
                    max_c = st.color_picker('Max color:', '#EE4B2B')
                    colors = [min_c, max_c]
                else:
                    min_c = st.color_picker('Min color:', '#0096FF')
                    center_c = st.color_picker('Center color:', '#ffffff')
                    max_c = st.color_picker('Max color:', '#EE4B2B')
                    colors = [min_c, center_c, max_c]
                c_cmap = create_custom_cmap(colors, name="custom")
                show_colormap(c_cmap)
            st.write(f"Selected colormap: {c_cmap}")
            cmap = c_cmap

            st.markdown('---')
        Plot_attr = st.checkbox('**Change plot range?**')
        if Plot_attr:
#            st.markdown("#### Plot range")
            center0 = st.checkbox('Set center as 0?', value = center0_z)
            if center0:
                v_center = 0
            else:
                v_center = None
            minmax = st.checkbox('Change min/center/max?')
            if minmax:
                st.write('min: ' + str(df.to_numpy().min()) + '    mean: ' + str(df.to_numpy().mean()) + "    max: " + str(df.to_numpy().max()) )
                if only_minmax:
                    v_min =  float(st.text_input("Min value",  value=df.to_numpy().min()))
                    v_max =  float(st.text_input("Max value",  value=df.to_numpy().max()))
                    v_center = None
                else:
                    v_min =  float(st.text_input("Min value",  value=df.to_numpy().min()))
                    v_center = float(st.text_input("Center",  value=df.to_numpy().mean()))
                    v_max =  float(st.text_input("Max value",  value=df.to_numpy().max()))

            st.markdown('---')

        annot_on = st.checkbox('Show value on each cell?')
        if annot_on:
            annot = True
            annot_digit = st.text_input("Number of decimal palces",  value=0)
            annot_digit = int(annot_digit)
    #        if annot_digit == 0:
    #            fmt = "d" #"d"はうまくいかない
    #        else:
            fmt = "." + str(annot_digit) + 'f'
    #        st.write(fmt)


        grid_on = st.checkbox('Show grid lines?')
        if grid_on:
            linewidths = st.number_input("Grid line width",  value=0.1, step=0.001, format="%.3f")
            line_col = st.selectbox("Grid line color", ('white', 'black','red','blue','gray','yellow'), index=0)
        else:
            line_col = 'white'
        
        # Boundary line options for clustering
        show_boundary = st.checkbox('Show cluster boundary lines?', value=False)
        if show_boundary:
            boundary_line = st.number_input('Boundary line width', value=0.5, step=0.1, format="%.1f")
            boundary_color = st.selectbox("Boundary line color", ('black', 'white', 'red','blue','gray','yellow'), index=0)
        else:
            boundary_line = 0.0
            boundary_color = 'black'

        st.markdown('---')


    ###### 現状plot sizeを変えられない
        py_x_size = 8
        py_y_size = 8

        st.markdown('#### Plot size')
        py_x_size = float(st.text_input("Plot x size:", value = 8))
        py_y_size = float(st.text_input("Plot y size:", value = 8))
        st.markdown('#### Font size')
    #    sns_font_scale = float(st.text_input("Font scale:", value = 1))
        x_font_size = float(st.text_input("Sample name (column) font size:", value = 12))
        y_font_size = float(st.text_input("Gene name (row) font size:", value = 12))

        xticklabels= "auto"
        yticklabels= "auto"
        x_all = st.checkbox("Show all sample (column) names?", value = False)
        y_all = st.checkbox("Show all gene (row) names? !!!Do not check this for a large number of genes!!!", value = False)
        if x_all:
            xticklabels=1

        if y_all:
            yticklabels=1

    st.markdown('##### Clustering:')
    if calc_z:
        # NMFは負の値を扱えないため、Z-score変換時は除外
        clustering_type = st.radio("Clustering:", ('Nonclustering','Hierarchical','k-means','x-means', 'g-means'), label_visibility='collapsed')
    else:
        clustering_type = st.radio("Clustering:", ('Nonclustering','Hierarchical','k-means','x-means', 'g-means', 'NMF'), label_visibility='collapsed')
    if clustering_type == 'k-means':
        from kneed import KneeLocator
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        from sklearn.preprocessing import StandardScaler
        st.markdown('##### k-means options:')
        elbow = st.checkbox("Draw elbow plot and determine K automaticllay?", value = False)
        if "k" not in st.session_state:
            st.session_state.k = 3
        else:
            k_number = st.session_state.k
        if elbow:
            if st.button('Generate Elbow Plot'):
                with st.spinner('Generating elbow plot...'):
                    try:
                        sse = []
                        K = range(1, 11)
                        for k in K:
                            kmeans = KMeans(n_clusters=k, init='k-means++', n_init='auto', random_state=42)
                            kmeans.fit(df)
                            sse.append(kmeans.inertia_)
                        fig, ax = plt.subplots()
                        ax.plot(K, sse, 'bx-')
                        ax.set_xlabel('k')
                        ax.set_ylabel('Sum of squared distances')
                        ax.set_title('Elbow Method For Optimal k')
                        st.pyplot(fig)

                        kl = KneeLocator(K, sse, curve="convex", direction="decreasing")
                        if kl.elbow:
                            st.write(f"Optimal K suggested by elbow method: {kl.elbow}")
                            k_number = kl.elbow
                            st.session_state.k = k_number
                        else:
                            st.write("Could not determine optimal K automatically. Please select K manually.")
                    except Exception as e:
                        st.error(f"An error occurred while generating the elbow plot: {str(e)}")

        else:
            k_number = int(st.number_input("K clusters:", value = 3))
            st.session_state.k = k_number
        st.write("K number: " + str(k_number))
    
    elif clustering_type == 'NMF':
        st.markdown('##### NMF options:')
        st.markdown('NMF works best with non-negative data. Negative values will be set to 0.')
        
        with st.expander("ℹ️ About NMF (Non-negative Matrix Factorization)"):
            st.markdown("""
            **NMF** decomposes your data matrix into two non-negative matrices:
            - **W matrix**: Gene loadings on components (genes × components)
            - **H matrix**: Component loadings on samples (components × samples)
            
            **When to use NMF:**
            - Best for count data (RNA-seq, ChIP-seq peaks)
            - Finding parts-based representations
            - Identifying additive biological processes
            - Topic modeling in gene expression
            
            **Parameter Guide:**
            - **Number of components**: Start with 2-5 for exploratory analysis. More components = finer granularity
            - **Initialization methods**:
              - `nndsvd`: Best for sparse data (recommended for most genomics data)
              - `nndsvda`: NNDSVD with zeros filled with small random values
              - `nndsvdar`: NNDSVD with zeros filled with the data average
              - `random`: Random initialization (less reproducible)
            - **Max iterations**: Usually 200 is sufficient; increase if not converged
            """)
        
        if "n_components" not in st.session_state:
            st.session_state.n_components = 3
        n_components = int(st.number_input("Number of components:", value = st.session_state.n_components, min_value=2, max_value=min(df.shape), 
                                         help="Number of latent components to extract. Each component represents a pattern in your data."))
        st.session_state.n_components = n_components
        
        nmf_init = st.selectbox("Initialization method:", ('nndsvd', 'nndsvda', 'nndsvdar', 'random'), index=0,
                               help="Algorithm for matrix initialization. 'nndsvd' is recommended for sparse genomics data.")
        
        nmf_max_iter = int(st.number_input("Max iterations:", value = 200, min_value=100, max_value=1000,
                                          help="Maximum number of iterations. Increase if the algorithm doesn't converge."))
        
        st.write(f"Number of components: {n_components}")

#    save_cluster = st.checkbox("Save cluster info?", value = False)
    st.markdown('---')
    y_c = False
    x_c = False
    if clustering_type == 'Hierarchical':
        import fastcluster
        y_c = st.checkbox("Cluster rows (Y axis)?", value = True)
        x_c = st.checkbox("Cluster colums (X axis)?", value = False)
        st.markdown('---')
        method_type = st.radio("Clustering method:", ('average','weighted', 'ward', 'median','single','centroid'))
        metric_type = st.radio("Clustering metric:", ('euclidean', 'seuclidean', 'sqeuclidean', 'minkowski', 'correlation', 'mahalanobis', 'cityblock', 'jaccard', 'jensenshannon'))
        st.markdown('---')


    save_type = st.radio("Save heatmap as: (Preparing PDF may take a time.)", ('png','pdf'))

    show_cor = st.checkbox('Show correlation coeficient matrix?')
    st.write("Clustering: " + clustering_type)
#    if "clustering" not in st.session_state:
#        st.session_state.clustering = False
    make_plot = st.button('Make plot')
    if make_plot:
        make_plot = False

        #両方Falseの場合はclustering_typeも変更する　#これがX-meansへいかさなくしている’
        if show_cor:
            correlation_coefficients = df.corr()
            fig_c, ax_c = plt.subplots() #この形式でないとエラーになる
            ax_c = sns.heatmap(correlation_coefficients, vmax=1, vmin=-1, cmap='seismic', square=True,
                annot=False, xticklabels=1, yticklabels=1)
            st.pyplot(fig_c)
            fig_c.savefig(res_dir + "/corrlation." + save_type, format=save_type)

        df_file_name = file_name_head + '.Data4Heatmap.tsv'
        df.to_csv(res_dir + "/" + df_file_name,sep= '\t')

        if clustering_type =="Nonclustering":
            #fig = plt.figure(figsize=(py_x_size, py_y_size))
#            sns.set(rc={'figure.figsize':(py_x_size, py_y_size)})
#            sns.set(font_scale=sns_font_scale)
#            fig = plt.figure()
#            g = sns.heatmap(df, cmap = cmap, center = v_center, vmin= v_min, vmax = v_max)
            g = sns.clustermap(df, center = v_center, cmap = cmap,
                    vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                    xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths, linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'Hierarchical':
            with st.spinner('Performing hierarchical clustering...'):
                # クラスタリング計算の実行（キャッシュ可能）
                row_linkage, col_linkage, row_distances, col_distances  = perform_clustering_computation(df.values, method_type, metric_type)
                # method typeがaverageなど

                # プロットの生成（キャッシュしない）
                g = plot_clustermap(
                    df, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c,
                    xticklabels, yticklabels, annot, fmt, linewidths, line_col,
                    py_x_size, py_y_size, x_font_size, y_font_size
                )
            st.pyplot(g)

        elif clustering_type == 'k-means':
            #if st.button('Calculate k-means'):
            with st.spinner('Performing k-means clustering...'):
                kmeans = KMeans(n_clusters=int(k_number), init ='k-means++', n_init='auto', random_state=42)
                clusters = kmeans.fit_predict(df)
                df2 = df.copy()
                df2["cluster"] = clusters
                df3 = pd.DataFrame(df2['cluster'], index = df2.index)
                st.write(df3.head(3))

                cluster_file_name = file_name_head + '.k-' + str(k_number) + '.tsv'

                df3.sort_values('cluster').to_csv(res_dir + "/" + cluster_file_name,sep= '\t')
                df2 = df2.sort_values('cluster')
                df2_sorted = df2.copy(deep=True)
                df2 = df2.drop('cluster', axis =1)
                g = sns.clustermap(df2, center = v_center, cmap = cmap,
                vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                         figsize = (py_x_size,py_y_size))
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
                if show_boundary:
                    ax = g.ax_heatmap
                    cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                    for boundary in cluster_boundaries[:-1]:
                        ax.axhline(y=boundary, color=boundary_color, linewidth=boundary_line)
                    # df_boundaryを作成
                    df_boundary = pd.DataFrame({
                        'position': cluster_boundaries
                    }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                    # sizeを計算（最初の行はpositionと同じ、それ以降は差分）
                    df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                    st.write(df_boundary)
                    df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_kmeans' + str(k_number) + '_ClusterSize.tsv',sep= '\t')



            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'x-means':
            st.write("Calculating x-means...")
            from pyclustering.cluster import xmeans
            initial_centers = xmeans.kmeans_plusplus_initializer(df, 2).initialize() # k=2以上で探索
            xm = xmeans.xmeans(df, initial_centers=initial_centers, )
            xm.process()
            clusters = xm.predict(df)
            df2 = df.copy()
            df2["cluster"] = clusters
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head + '.xmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            df2 = df2.sort_values('cluster')
            df2_sorted = df2.copy(deep=True)
            df2 = df2.drop('cluster', axis =1)
            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
            if boundary_line > 0:
                ax = g.ax_heatmap
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                for boundary in cluster_boundaries[:-1]:
                    ax.axhline(y=boundary, color=line_col, linewidth=boundary_line)
                # df_boundaryを作成
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                # sizeを計算（最初の行はpositionと同じ、それ以降は差分）
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_xmeans_ClusterSize.tsv',sep= '\t')
            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'g-means':
            from pyclustering.cluster import gmeans
            ar = df.to_numpy()
            with st.spinner('This takes a long time...'):
                initial_centers = gmeans.kmeans_plusplus_initializer(ar, 2).initialize() # k=2以上で探索
                gm = gmeans.gmeans(ar, initial_centers=initial_centers, )
                gm.process()
                clusters = gm.predict(ar)
            df2 = df.copy()
            df2["cluster"] = clusters
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head +  '.gmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            df2 = df2.sort_values('cluster')
            df2_sorted = df2.copy(deep=True)
            df2 = df2.drop('cluster', axis =1)
            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
            if boundary_line > 0:
                ax = g.ax_heatmap
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                for boundary in cluster_boundaries[:-1]:
                    ax.axhline(y=boundary, color=line_col, linewidth=boundary_line)
                # df_boundaryを作成
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                # sizeを計算（最初の行はpositionと同じ、それ以降は差分）
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_gmeans_ClusterSize.tsv',sep= '\t')

            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'NMF':
            with st.spinner('Performing NMF decomposition...'):
                # Ensure non-negative data for NMF
                df_nmf = df.copy()
                if df_nmf.min().min() < 0:
                    st.warning("Negative values detected. Setting them to 0 for NMF.")
                    df_nmf = df_nmf.clip(lower=0)
                
                # Perform NMF
                nmf_model = NMF(n_components=n_components, init=nmf_init, max_iter=nmf_max_iter, random_state=42)
                W = nmf_model.fit_transform(df_nmf)  # Document-topic matrix (genes x components)
                H = nmf_model.components_  # Topic-term matrix (components x samples)
                
                # Assign each gene to the component with the highest weight
                cluster_assignments = np.argmax(W, axis=1)
                
                # Add cluster information to dataframe
                df2 = df.copy()
                df2["cluster"] = cluster_assignments
                df3 = pd.DataFrame(df2['cluster'], index = df2.index)
                
                st.write("NMF Component Assignments:")
                st.dataframe(df3.head(10))
                
                # Save cluster assignments
                cluster_file_name = file_name_head + f'.NMF_{n_components}components.tsv'
                df3.sort_values('cluster').to_csv(res_dir + "/" + cluster_file_name, sep='\t')
                
                # Save W matrix (gene loadings on components)
                W_df = pd.DataFrame(W, index=df.index, columns=[f'Component_{i}' for i in range(n_components)])
                W_df.to_csv(res_dir + "/" + file_name_head + f'.NMF_W_matrix_{n_components}comp.tsv', sep='\t')
                
                # Save H matrix (component loadings on samples)
                H_df = pd.DataFrame(H, index=[f'Component_{i}' for i in range(n_components)], columns=df.columns)
                H_df.to_csv(res_dir + "/" + file_name_head + f'.NMF_H_matrix_{n_components}comp.tsv', sep='\t')
                
                # Sort genes by cluster for visualization
                df2 = df2.sort_values('cluster')
                df2_sorted = df2.copy(deep=True)
                df2 = df2.drop('cluster', axis=1)
                
                # Create heatmap
                g = sns.clustermap(df2, center = v_center, cmap = cmap,
                    vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                    xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, 
                    linewidths= linewidths, linecolor=line_col,
                    figsize = (py_x_size, py_y_size))
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
                
                # Add cluster boundaries
                if show_boundary:
                    ax = g.ax_heatmap
                    cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                    for boundary in cluster_boundaries[:-1]:
                        ax.axhline(y=boundary, color=boundary_color, linewidth=boundary_line)
                    
                    # Create boundary dataframe
                    df_boundary = pd.DataFrame({
                        'position': cluster_boundaries
                    }, index=pd.Index(range(len(cluster_boundaries)), name='component'))
                    df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                    st.write("Component sizes:")
                    st.write(df_boundary)
                    df_boundary.to_csv(res_dir + "/" + file_name_head + f'_NMF_{n_components}_ComponentSize.tsv', sep='\t')
                
                st.pyplot(g)
                
                # Display reconstruction error
                reconstruction_error = nmf_model.reconstruction_err_
                st.write(f"NMF Reconstruction Error: {reconstruction_error:.4f}")
                
                # Display W matrix as heatmap
                with st.expander("View Gene-Component Matrix (W)"):
                    st.markdown("**W matrix**: Shows how much each gene contributes to each component")
                    # W matrixのサイズに応じて図のサイズを調整
                    w_fig_height = min(max(len(df.index) * 0.15, 5), 20)  # 高さを遺伝子数に応じて調整（最小5、最大20）
                    fig_w, ax_w = plt.subplots(figsize=(n_components + 2, w_fig_height))
                    
                    # 上位遺伝子のみ表示するオプション
                    show_top_genes = st.checkbox("Show only top genes per component", value=True, key="w_matrix_top")
                    if show_top_genes:
                        top_n_genes = st.slider("Number of top genes per component:", min_value=5, max_value=50, value=20, key="w_matrix_topn")
                        # 各コンポーネントでトップN遺伝子を取得
                        top_gene_indices = set()
                        for comp in range(n_components):
                            top_indices = np.argsort(W[:, comp])[-top_n_genes:]
                            top_gene_indices.update(top_indices)
                        top_gene_indices = sorted(list(top_gene_indices))
                        W_display = W_df.iloc[top_gene_indices, :]
                        st.write(f"Showing top {top_n_genes} genes per component ({len(top_gene_indices)} unique genes total)")
                    else:
                        W_display = W_df
                        st.write(f"Showing all {len(W_df)} genes")
                    
                    sns.heatmap(W_display, cmap='YlOrRd', annot=False, xticklabels=1, 
                               yticklabels=1 if len(W_display) <= 50 else 'auto', ax=ax_w, cbar_kws={'label': 'Weight'})
                    ax_w.set_title(f"Gene loadings on NMF components (W matrix)")
                    ax_w.set_xlabel("Components")
                    ax_w.set_ylabel("Genes")
                    st.pyplot(fig_w)
                    
                    # 各コンポーネントのトップ遺伝子をテキストで表示
                    st.markdown("**Top genes per component:**")
                    for comp in range(n_components):
                        top_gene_idx = np.argsort(W[:, comp])[-10:][::-1]  # トップ10遺伝子
                        top_genes = df.index[top_gene_idx].tolist()
                        top_weights = W[top_gene_idx, comp]
                        st.write(f"Component_{comp}: {', '.join([f'{gene} ({weight:.3f})' for gene, weight in zip(top_genes[:5], top_weights[:5])])}")
                
                # Display H matrix as heatmap
                with st.expander("View Component-Sample Matrix (H)"):
                    st.markdown("**H matrix**: Shows how much each component is expressed in each sample")
                    fig_h, ax_h = plt.subplots(figsize=(py_x_size, n_components))
                    sns.heatmap(H_df, cmap='YlGnBu', annot=False, xticklabels=1, yticklabels=1, ax=ax_h, cbar_kws={'label': 'Weight'})
                    ax_h.set_title("NMF component expression across samples (H matrix)")
                    ax_h.set_xlabel("Samples")
                    ax_h.set_ylabel("Components")
                    st.pyplot(fig_h)
            
            st.markdown('---')


        if howlog == "No":
            logmethod = ""
        else:
            logmethod = "_" + howlog
        if calc_z:
            logmethod = logmethod + '.Z'
        if save_type == 'pdf':
            if clustering_type == 'k-means':
                file_name = file_name_head + logmethod +  '.k-' + str(k_number) + '.heatmap.pdf'
            elif clustering_type == 'NMF':
                file_name = file_name_head + logmethod + f'.NMF-{n_components}.heatmap.pdf'
            else:
                file_name = file_name_head + logmethod + '.heatmap.pdf'
        else:
            if clustering_type == 'k-means':
                file_name = file_name_head + logmethod +  '.k-' + str(k_number) + '.heatmap.png'
            elif clustering_type == 'NMF':
                file_name = file_name_head + logmethod + f'.NMF-{n_components}.heatmap.png'
            else:
                file_name = file_name_head + logmethod + '.heatmap.png'
        try:
            with st.spinner('Generating PDF figure file may take a time...'):
                g.savefig(res_dir  + "/" + file_name, format=save_type)
        except:
            pass
        else:
            shutil.make_archive(temp_dir + "/Heatmap", format='zip',root_dir= res_dir)
            with open(temp_dir + "/Heatmap.zip", "rb") as fp:
                btn = st.download_button(
                    label="Download Results?",
                data=fp,
                file_name= file_name_head  + logmethod + '_' + clustering_type + ".Heatmap.zip",
                mime = "zip"
                )
#            shutil.rmtree(temp_dir)
#            os.mkdir(temp_dir)