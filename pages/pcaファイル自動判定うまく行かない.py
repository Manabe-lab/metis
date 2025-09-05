import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import plotly.express as px
import io
from sklearn.preprocessing import StandardScaler
scaler = StandardScaler()
from detect_delimiter import detect

st.set_page_config(page_title="PCA.", page_icon="√")

def normalize_totalreads(df):
    return 10**6 * df / df.sum()

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',',  header = 0):
    df_xl = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
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
def tsne_fit(df,n_components= 2, random_state=0, perplexity = 30):
    embedding = TSNE(n_components=n_components, random_state=random_state, perplexity = perplexity)
    emb = embedding.fit_transform(df)
    return emb

@st.cache_data
def umap_fit(df,n_components=2, random_state=0, n_neighbors=15):
    embedding = UMAP.UMAP(n_components=n_components, random_state=random_state, n_neighbors=n_neighbors)
    emb = embedding.fit_transform(df)
    return emb


def homer_clean(df):
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
    else:
        df.set_index(content[0], inplace = True)
    return(df)

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #   st.write('file name kept')
    filename_add = ""

#input_file_type = st.radio(
#    "Data format:",
#    ('tsv','csv', 'excel', 'Homer'), key='tsv')

input_file_format = 'row = gene'


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

st.write('')

def get_delimiter(file_path, bytes = 4096):
    sniffer = csv.Sniffer()
    data = open(file_path, "r").read(bytes)
    delimiter = sniffer.sniff(data).delimiter
    return delimiter

uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
if uploaded_file is not None:
    try:
        df = read_csv(uploaded_file)
        st.write('Read CSV file.')
        if (df.shape[1]) == 1 or ('\t' in df.iloc[0,0]):
            st.write('Re-reading as TSV.')
            try:
                df = read_csv(uploaded_file,sep = '\t')
            except:
                st.write("Possible format error.")
    except:
        try:
            df = read_csv(uploaded_file,sep = '\t')
            st.write('Read TSV file.')
        except:
            df = read_excel(uploaded_file)
        else:
            st.write("Format error. Please reload the page and upload the file again.")
#    df = read_csv(uploaded_file, sep = '\t')
    st.write(df.head(3))


    df = homer_clean(df)


###########

    if input_file_format == 'column = gene':
        df = df.T

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

    normalize_df = st.checkbox('Normalize count data?')
    if normalize_df:
        convert_to = st.radio(
        "Convert to:",
    ('TMM','UQ','CTF', 'CUF', 'CPM', 'RPKM to TPM'), key='TMM')


    if normalize_df:
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

        if input_file_format == 'row = gene':
            df_conv = df_conv.T

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

    log_df = st.checkbox('Log transforamtion?')
    if log_df:
        log_transform = st.radio(
        "Log transformation:",
        ('None', 'asinh','log2(x+1)', 'ln(x+1)', 'log10(x+1)'), key='None')


    category = st.checkbox('Create categories for coloring?')
    if category:
        df_e = pd.DataFrame(df.columns.tolist(), index =df.columns.tolist() , columns = ["Group"]).copy()
        edited_df_e = st.data_editor(df_e)
        condition = edited_df_e.iloc[:,0].tolist()
        st.write('Group: ' + '  '.join(condition))



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
        ('PCA', 'tSNE','UMAP'), key='PCA')
        reduc_dim = st.radio(
        "Dimension:",
        ('2D', '3D'), key='2D')
        if reduc_dim ==  '3D':
            n_components = 3
        else:
            n_components = 2
        if reduc == 'PCA':
            x_show = int(st.text_input("X dim:", value = 1))
            y_show = int(st.text_input("Y dim:", value = 2))
            if x_show > 15 or y_show > 15:
                st.write("Max dimension is 15")
            else:
                pca_x = 'PCA' + str(x_show)
                pca_y = 'PCA' + str(y_show)
            if reduc_dim == '3D':
                z_show = int(st.text_input("Z dim:", value = 3))
                if z_show > 15 :
                    st.write("Max dimension is 15")
                else:
                    pca_z = 'PCA' + str(z_show)

        if reduc == 'tSNE' or reduc == 'UMAP':
            calc_z = st.checkbox('Convert to Z score?')
            if calc_z:
                z  = scaler.fit_transform(df.T)
                df = pd.DataFrame(z.T, columns=df.columns.to_list(), index=df.index.to_list())
        if reduc == 'tSNE':
            perplexity = float(st.text_input("Perplexity:", value = 30))
        if reduc == "UMAP":
            n_neighbors = int(st.text_input("n_neighbors:", value = 15)) # integer



#        cmap  = st.radio(
#        "Color map:",
#        ('None', 'Set1','Set2','Set3','tab10','tab20', 'random'), key='None')

        st.markdown('##### Plot attributes')
        width = float(st.text_input("Plot x size:", value = 600))
        height = float(st.text_input("Plot y size:", value = 500))
        marker_size = float(st.text_input("Dot size:", value = 6))
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
            from sklearn.decomposition import PCA
            z  = scaler.fit_transform(df.T)
            df_z = pd.DataFrame(z.T, columns=df.columns.to_list(), index=df.index.to_list())
            x_embedded, pca = pca_fit(df_z.T)
            pca_col_name = ["PCA" + str(i+1) for i in range(x_embedded.shape[1])]
            if category:
                df_pca = pd.DataFrame(x_embedded[:,:], index= condition, columns=pca_col_name)
            else:
                df_pca = pd.DataFrame(x_embedded[:,:], index= df.columns.to_list(), columns=pca_col_name)

            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_pca, x= pca_x, y= pca_y, color=df_pca.index, color_discrete_sequence=c_choice)

                else:
                    fig = px.scatter(df_pca, x= pca_x, y= pca_y, color=df_pca.index, text = df_pca.index , color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
                st.write("Explained variance ratio:")
                st.write(pca_x + ": " + str(pca.explained_variance_ratio_[x_show - 1]))
                st.write(pca_y + ": " + str(pca.explained_variance_ratio_[y_show - 1]))
            else:
                if not show_text:
                    fig = px.scatter_3d(df_pca, x= pca_x, y= pca_y, z= pca_z,  color=df_pca.index, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_pca, x= pca_x, y= pca_y, z= pca_z, color=df_pca.index, text = df_pca.index, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
                st.write(pca_x + ": " + str(pca.explained_variance_ratio_[x_show - 1]))
                st.write(pca_y + ": " + str(pca.explained_variance_ratio_[y_show - 1]))
                st.write(pca_z + ": " + str(pca.explained_variance_ratio_[z_show - 1]))
            st.plotly_chart(fig)

        elif reduc == 'tSNE':
            from sklearn.manifold import TSNE
            if perplexity >= color_num:   # perplexityはサンプル数より少ない
                perplexity = color_num - 1
                st.write('perplexity is set to ' + str(perplexity))
#            embedding = TSNE(n_components=n_components, random_state=0, perplexity = perplexity)
            x_embedded = tsne_fit(df.T,n_components=n_components, perplexity = perplexity )
            if reduc_dim == "2D":
                df_tsne = pd.DataFrame(x_embedded[:,:2], index= df.columns.to_list(), columns=["tSNE1","tSNE2"])
            else:
                df_tsne = pd.DataFrame(x_embedded[:,:3], index= df.columns.to_list(), columns=["tSNE1","tSNE2",'tSNE3'])

            if category:
                df_tsne.index = condition

            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_tsne, x= 'tSNE1', y= 'tSNE2', color=df_tsne.index, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_tsne, x= 'tSNE1', y= 'tSNE2', color=df_tsne.index, text = df_tsne.index, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            else:
                if not show_text:
                    fig = px.scatter_3d(df_tsne, x= 'tSNE1', y= 'tSNE2',z='tSNE3',  color=df_tsne.index, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_tsne, x= 'tSNE1', y= 'tSNE2',z='tSNE3',  color=df_tsne.index, text = df_tsne.index, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            st.plotly_chart(fig)

        elif reduc == 'UMAP':
            import umap.umap_ as UMAP
            if n_neighbors >= color_num:   # perplexityはサンプル数より多くないといけない
                n_neighbors = color_num - 1
                st.write('n_neighbors is set to ' + str(n_neighbors))
            x_embedded = umap_fit(df.T, n_components=n_components, random_state=0, n_neighbors=n_neighbors)
            if reduc_dim == "2D":
                df_umap = pd.DataFrame(x_embedded[:,:2], index= df.columns.to_list(), columns=["UMAP1","UMAP2"])
            else:
                df_umap = pd.DataFrame(x_embedded[:,:3], index= df.columns.to_list(), columns=["UMAP1","UMAP2",'UMAP3'])

            if category:
                df_umap.index = condition
            if reduc_dim == '2D':
                if not show_text:
                    fig = px.scatter(df_umap, x= 'UMAP1', y= 'UMAP2', color=df_umap.index, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter(df_umap, x= 'UMAP1', y= 'UMAP2', color=df_umap.index, text = df_umap.index, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)
            else:
                if not show_text:
                    fig = px.scatter_3d(df_umap, x= 'UMAP1', y= 'UMAP2',z='UMAP3', color=df_umap.index, color_discrete_sequence=c_choice)
                else:
                    fig = px.scatter_3d(df_umap, x= 'UMAP1', y= 'UMAP2',z='UMAP3', color=df_umap.index, text = df_umap.index, color_discrete_sequence=c_choice)
                    fig.update_traces(textposition='top center', textfont_size=textfont_size)
                fig.update_traces(marker_size = marker_size)
                fig.update_layout(width = width, height=height)

            st.plotly_chart(fig)

        if saveas == 'html':
            #fig.write_html(file=buffer )
            fig.write_html('dummy.html' )
        else:
            fig.write_image(file=buffer, format=saveas, engine="kaleido" )
        file_name = os.path.splitext(uploaded_file.name)[0] + '.' + reduc + '.' + saveas

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
            loading_file_name = os.path.splitext(uploaded_file.name)[0] + '.LoadingsMatrix.tsv'
            st.download_button(label="Download loadings matrix",
                                data=loading_tsv,
                                file_name=loading_file_name,
                                mime='test/csv')

#    if (cmap == 'Set1' and color_num >9 ) or (cmap == 'Set2' and color_num > 8 ) or  (cmap == 'Set3' and color_num > 8 ) or (cmap == 'tab10' and  color_num > 10 ) or (cmap == 'tab20' and  color_num > 20 ):
#        st.markdown("### Too few color variatios. Choose another colormap")
#    else:
#        if cmap == 'random':
#            colors = np.random.randn(len(color_num ))
#            cmap = 'viridis'
#        else:
#            colors = range(color_num)


           # fig, ax = plt.subplots(figsize=(width, height)) #例にあるように、このようにaxisを指定して、axisに描画しないとエラーになる
           # if cmap == "None":
           #     ax.scatter(x_embedded[:, 0], x_embedded[:, 1])
           # else:
           #     ax.scatter(x_embedded[:, 0], x_embedded[:, 1], c=colors, cmap = cmap )
           #     ax.legend(loc='upper right')
           # st.pyplot(fig )
