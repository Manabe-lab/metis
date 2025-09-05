import streamlit as st
import pandas as pd
import csv
import re
import os

st.set_page_config(page_title="Convert_Homer_and_DESeq2_output_to_data_only", page_icon="📃")

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

input_file_type = st.radio(
    "File format (tsv: tab-separeted values, csv: comma-separated values)",
    ('Homer', 'tsv','csv','Homer-excel'))

rlog = False
if input_file_type ==  'tsv' or input_file_type ==   'csv':
    rlog = st.checkbox('Extract rlog data from DESeq2 output')

refseq = False
if input_file_type == "Homer" or input_file_type == "Homer-excel":
    refseq = st.checkbox('Output RefSeq IDs instead of gene symbols for Homer data')

uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv','csv','xlsx','xls'])
if uploaded_file is not None:
    if input_file_type == 'csv':
        df = pd.read_csv(uploaded_file)
    elif input_file_type ==  'Homer-excel':
        df = pd.read_excel(uploaded_file,  header = 0)
    else:
        df = pd.read_csv(uploaded_file, sep = '\t')

    st.write(df.head())

    if input_file_type == "Homer" or input_file_type ==  "Homer-excel":
        if refseq:
            df = df.drop(df.columns[range(1,8)], axis=1)
            colnames = df.columns.tolist()
            colnames[0] = 'RefSeq'
        else:
            df = df.iloc[:,7:]
            colnames = df.columns.tolist()
            colnames[0] = 'Gene'
    else:
        colnames = df.columns.tolist()
        colnames[0] = 'Gene'

    # colnamesの変換
    search_word = '([^\ \(]*)\ \(.*'

    if rlog:
        for i in range(1, len(colnames)):
            if ('Fold Change' in colnames[i]) or ('log2FC' in colnames[i]) :
                colnames = colnames[:i]
                df = df.iloc[:,:i]
                break
            else:
                match = re.search(search_word, colnames[i])
                if match:
                    colnames[i] = match.group(1).replace(' ', '_')
    else:
        for i in range(1, len(colnames)):
            match = re.search(search_word, colnames[i])
            if match:
                colnames[i] = match.group(1).replace(' ', '_')

    df.columns = colnames
    if refseq:
        df['RefSeq'] = df['RefSeq'].astype(str) # excel対応
    else:
        df['Gene'] = df['Gene'].astype(str) # excel対応


    pattern = "([^|]*)"
    repatter = re.compile(pattern)
    f_annotation = lambda x: repatter.match(x).group(1)
    df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)


########## excel対応?
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

    df = df.set_index(df.columns[0])


    st.markdown("### Data file")
    st.write(df.head())

#    st.write(os.path.splitext(uploaded_file.name)[0])

    if rlog:
        file_name = os.path.splitext(uploaded_file.name)[0] + '.rlog.txt'
        st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0] + ".rlog"
    else:
        file_name = os.path.splitext(uploaded_file.name)[0] + '.data.txt'
        st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0]

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