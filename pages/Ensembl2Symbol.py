import streamlit as st
import pandas as pd
import pickle
import sys
import os
import re
import csv
import collections


st.set_page_config(page_title="Convert Ensembl ID to Symbol.", page_icon="❱")
st.markdown("### Can also summarize duplicated genes.")

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

species = st.radio(
    "Species:",
    ('Mouse','Human','No conversion, summarize duplicates'))

canonical = st.checkbox('Use only canonical gene/transcript?', value = True)


data_format = st.radio(
    "Data format:",
    ('TSV tab-separated','CSV comma-separated'))


uploaded_file = st.file_uploader("Choose a file", type=['txt','csv','tsv'])
if uploaded_file is not None:
    if data_format == "CSV comma-separated":
        df = pd.read_csv(uploaded_file, index_col = 0)
    else:
        df = pd.read_csv(uploaded_file, sep = '\t', index_col = 0)


    st.markdown("### Original")
    st.write('data lengths: ' + str(len(df)))
    st.write(df.head())

    col_names = df.columns.tolist()
#    col_names[0] = "Gene"
#    df.columns = col_names
    df.index.name = "Gene"
    st.write(df.head())

    OutDataName = os.path.splitext(uploaded_file.name)[0] + '.Symbol.txt'

    if species != 'No conversion, summarize duplicates':
        if canonical:
            if species == 'Human':
                with open("/home/cellxgene/streamlit/data/Hg38_ENSEMBLCanonicalSmbol.pkl", "rb") as tf:
                    dic = pickle.load(tf)
            else:
                with open("/home/cellxgene/streamlit/data/ENSEMBLCanonicalSmbol.pkl", "rb") as tf:
                    dic = pickle.load(tf)
        else:
            if species == 'Human':
                with open("/home/cellxgene/streamlit/data/Human_ENSMUSGT2Symbol.pkl", "rb") as tf:
                    dic = pickle.load(tf)
            else:
                with open("/home/cellxgene/streamlit/data/ENSMUSGT2Symbol.pkl", "rb") as tf:
                    dic = pickle.load(tf)

        gene_list = df.index.tolist()

        ens_id = []
        converted_id = []
        p = re.compile(r'([^\.]*)') # transcript等の.1, .2等を除いてサーチする
        n= 0
        removed = []
        exist_id = []
        removed_genes = []
        for i in gene_list: # Homer Refseq ID
            gene = i
            i_word = p.match(i).group(1)
            try:
                new_id = dic[i_word]
            except: #結果がないときはgene nameでサーチ
                ens_id.append(i_word) # どちらも存在しないとき
                removed.append(i_word)
                converted_id.append(i_word)
                removed_genes.append([i_word,gene])
            else:
                ens_id.append(new_id)
                exist_id.append(i_word)
                converted_id.append(new_id)
        if len(df) != len(ens_id):
            print("The number of converted ids is different from the original.")
            sys.exit()

        st.markdown("### Unconverted genes:")
        st.write("        " + str(len(removed_genes)) + " genes.")
        df_removed = pd.DataFrame(removed_genes, columns = ['Ensembl ID','Symbol'])
        st.write(df_removed)

        df.index = converted_id

        df.index.name = "Gene"

        st.markdown("### Updated data")
        st.write(df.head())

    st.markdown("### Duplicated")
    dup_d = df.loc[df.index.duplicated(keep=False),:].sort_index()
    st.write("Dupllicated genes: " + str(len(set(dup_d.index))))
    st.write(dup_d)
    grouping = df.groupby(level = 0)

    remove_un = st.radio(
        "Remove the unconverted genes?",
        ('Yes','Keep all'))

    if remove_un == 'Yes':
        if species == 'Human':
            st.write(df.head())
            ensmus = ["ENSG" not in str(i) for i in df.index.tolist()]
            df = df[ensmus]

        else:
            ensmus = ["ENSMUS" not in str(i) for i in df.index.tolist()]
            df = df[ensmus]
            st.write(len(ensmus))

    st.write('data length: ' + str(len(df)))

    agg_method = st.radio(
        "Aggregaton method:",
        ('Max','Mean'))

    if st.button('You can aggregate duplicates'):
        dup_gene = set(df.index[df.index.duplicated(keep=False)])
        df_nodup = df[~df.index.duplicated(keep='first')] # indexの重複をもとに削除
        if agg_method == "Mean":
            df_mean = grouping.mean(numeric_only=True)
        else:
            df_mean = grouping.max(numeric_only=True)

        for i in dup_gene:
            # 列名を明示的に指定して値を代入
            df_nodup.loc[i, df_mean.columns] = df_mean.loc[i, df_mean.columns].values

        st.write(df_nodup.loc[list(dup_gene),:].sort_index())
        df = df_nodup

    st.markdown("""
---
"""        )

    st.write('data length: ' + str(len(df)))

    csv = convert_df(df)
    st.download_button(
       "Press to download updated data",
       csv,
       OutDataName,
       "text/csv",
       key='download-csv'
    )
