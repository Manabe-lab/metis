import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np

st.set_page_config(page_title="Log_and_Z_transformation.", page_icon="√")


@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep)
    return df_c

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')


def homer_clean(df):
    content = df.columns.tolist()
    Gene_column = content[0]
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
    return(df)

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #   st.write('file name kept')
    filename_add = ""

f_inf = -float('inf')
min_val = f_inf
max_val = f_inf
delta_val = 1
fold_val = 1
min_variance = 0

input_file_type = st.radio(
    "Data format:",
    ('tsv','csv', 'excel'))

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
    if input_file_type == "csv":
        df = read_csv(uploaded_file, index_col = 0)
    elif input_file_type == 'tsv':
        df = read_csv(uploaded_file, sep = '\t', index_col = 0)
    else:
        df = read_excel(uploaded_file, index_col = 0)

    # excel並びにHomer対応
    df = homer_clean(df)


    st.write('Original:')
    st.write(df.head())

    st.markdown("######   ")
    nonzero = st.checkbox('Remove all zero genes?')

    howlog = 'No'
    st.markdown("######   ")
    calc_log = st.checkbox('Log transformation?')
    if calc_log:
        howlog = st.radio('Method', ['No', 'asinh', 'log2+1', 'log2', 'loge+1', 'loge','log10+1','log10'])

    st.markdown("######   ")
    calc_z = st.checkbox('Z-score?')

    st.markdown("######   ")
    more_filt = st.checkbox('More filtering options')
    if more_filt:
        min_val =  float(st.text_input("All values of each gene are larger than",  value=f_inf))
        max_val =  float(st.text_input("Max value is larger than",  value=f_inf))
        delta_val =  float(st.text_input("Delta (max - min value) >",  value=0))
        fold_val =  float(st.text_input("Fold (max / min) >",  value=1))
        min_variance =  float(st.text_input("Minmum variance across samples > ", value=0))

    if input_file_format == 'column = gene':
        df = df.T
    if nonzero:
        df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く

    if min_val != f_inf:
        df = df[df.apply(min, axis=1) > min_val]
        if "min" not in filename_add:
            filename_add = '.min' + str(min_val) + filename_add

    if max_val != f_inf:
        df =  df[df.apply(max, axis=1) > max_val] #ここがminだと、0が一つでもあれば削除される。
        if "max" not in filename_add:
            filename_add = '.max' + str(max_val) + filename_add

    if delta_val > 1:
        df = df[df.apply(max, axis=1) > df.apply(min, axis=1) + delta_val]
        if "delta" not in filename_add:
            filename_add = '.delta' + str(delta_val) + filename_add

    if fold_val > 1:
        df = df[df.apply(max, axis=1) > df.apply(min, axis=1) * fold_val]
        if "fold" not in filename_add:
            filename_add = '.fold' + str(fold_val) + filename_add

    if min_variance > 0:
        df = df.loc[(df.index[(df.var(axis=1) > min_variance)]),:]

    if howlog == 'asinh':
        df = np.arcsinh(df)

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
        df = np.log2(df)

    if howlog != 'No':
        filename_add = '.log' + filename_add

    if calc_z:
        df_z = df.copy()
        m = df_z.mean(1)
        s = df_z.std(1)
        df_z = df_z.sub(m, axis=0).div(s, axis = 0)
        df_z = np.round(df_z, decimals=10)
        filename_add = filename_add + ".Z"
        df = df_z

    if input_file_format == 'column = gene':
        df = df.T


    file_name = os.path.splitext(uploaded_file.name)[0] + filename_add + '.txt'

    st.markdown('#### Modified:')
    st.write("Number of data: " + str(len(df)))
    st.write(df.head())

    st.session_state.df = df
    st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0]

    csv = convert_df(df)
    st.download_button(
      "Press to Download",
       csv,
       file_name,
       "text/csv",
       key='download-csv'
       )
