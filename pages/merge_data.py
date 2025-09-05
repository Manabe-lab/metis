import streamlit as st
import pandas as pd
import os
import itertools
import collections

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep)
    return df_c

st.set_page_config(page_title="Data_file_merger")

# Title
st.title('Data File Merger')
st.write("Data file format must be:")
st.markdown("""
    row = gene
    |  | Sample1 | Sample2 |
    | --- | --- | --- |
    | Gene1 |  |  |
    | Gene2 | ||

    """)

st.write(" ")

input_file_type = st.radio(
    "Data format:",
    ('tsv','csv', 'excel'), key='tsv')

# Create a list to hold DataFrames
dfs = []

# Allow multiple files upload
uploaded_files = st.file_uploader("Upload Data Files", type=['txt','tsv','csv','xls','xlsx'], accept_multiple_files=True)

if uploaded_files:
    file_number = 0
    for uploaded_file in uploaded_files:
        if uploaded_file:
            if input_file_type == "csv":
                df = read_csv(uploaded_file, index_col = 0)
            elif input_file_type == "excel":
                df = read_excel(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep = '\t', index_col = 0)
            dfs.append(df)
            file_number += 1
    st.write(f"Uploaded {file_number} files.")
    st.write("First data")
    st.write(dfs[0].head())


null_handler = st.radio("How to handle null data", ["Fill 0", "Delete the genes containing null"], index = 0)

if st.button("Merge files"):

    # Concatenate all dataframes in the list
    if dfs:
        col_names = list(itertools.chain.from_iterable([x.columns.values for x in dfs]))
        dup_list = [k for k, v in collections.Counter(col_names).items() if v > 1]
        if len(dup_list) > 0:
            st.write("Error!!! There are duplicated columns.")
            st.write(dup_list)
        else:
            file_number = 1
            for x in dfs:
                if x.index.duplicated().sum() > 0:
                    st.write("There are duplicated rows in file: " + str(file_number))
                    st.write(x.index[x.index.duplicated()])
                file_number += 1

            merged_df = pd.concat(dfs,join='outer', axis =1)
            #merged_df = pd.concat(dfs, axis =1)

            # Show the merged dataframe
            if merged_df.isnull().values.sum() > 0:
                st.write("There are " + str(merged_df.isnull().values.sum()) + "missing data.")
                if null_handler == "Fill 0":
                    merged_df = merged_df.fillna(0)
                    st.write("Missing data are filled with 0.")
                else:
                    merged_df = merged_df.dropna()
                    st.write("Genes containing missing data are removed.")

            st.write('Merged data')
            if merged_df.shape[1] > 10:
                st.write(merged_df.iloc[:5,:10])
            else:
                st.write(merged_df.head())

            st.session_state.df = merged_df
            file_name = ""
            for uploaded_file in uploaded_files:
                if file_name != "":
                    file_name = file_name + "_" + os.path.splitext(uploaded_file.name)[0]
                else:
                    file_name = file_name +  os.path.splitext(uploaded_file.name)[0]

            st.session_state.uploaded_file_name = file_name
            st.write(file_name)

            csv = convert_df(merged_df)
            st.download_button(
               "Press to Download",
               csv,
               file_name + '.merged.tsv',
               "text/csv",
               key='download-csv'
            )
    else:
        st.write("No files uploaded")
