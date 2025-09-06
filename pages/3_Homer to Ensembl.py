import streamlit as st
import pandas as pd
import pickle
import sys
import os
import re
import csv
import collections


st.set_page_config(page_title="Convert symbols to Ensembl IDs.", page_icon="❱")

@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

species = st.radio(
    "Species:",
    ('Mouse','Human'))

data_type = st.radio(
    "Data type:",
    ('Homer','Generic data file'))

data_format = st.radio(
    "Data format:",
    ('TSV tab-separated','CSV comma-separated'))


uploaded_file = st.file_uploader("Choose a file", type=['txt','csv','tsv'])
if uploaded_file is not None:
    if data_format == "CSV comma-separated":
        df = pd.read_csv(uploaded_file, index_col = 0)
    else:
        df = pd.read_csv(uploaded_file, sep = '\t', index_col = 0)

    # Homerの場合はrefseq idをもとに変換する
    st.markdown("### Original")
    st.write('data lengths: ' + str(len(df)))
    st.write(df.head(3))

    if data_type == "Homer":
        df = df.iloc[:,6:]
        col_names = df.columns.tolist()
        col_names[0] = "Gene"
    else:
        df.index.name = "Gene"
    st.write(df.head(3))

    if species == 'Human':
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "human_synonym2symbol.pkl"), "rb") as tf:
            synonym = pickle.load(tf)
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "Human_refseq2ensembl_2023-4.pkl"), "rb") as tf:
            refseq = pickle.load(tf)
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "Human_Symbol2Ensembl_2023-4.pkl"), "rb") as tf:
            dic = pickle.load(tf)
    else:
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "mouse_synonym2symbol.pkl"), "rb") as tf:
            synonym = pickle.load(tf)
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "mouse_refseq2ensembl_2023-4.pkl"), "rb") as tf:
            refseq = pickle.load(tf)
        with open(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "Symbol2Ensembl_2023-4.pkl"), "rb") as tf:
            dic = pickle.load(tf)

    OutDataName = os.path.splitext(uploaded_file.name)[0] + '.Ensembl.txt'

    gene_list = df.index.tolist()

    ens_id = []
    converted_id = []
    p = re.compile(r'([^|]*)')
    n= 0
    removed = []
    exist_id = []
    removed_genes = []
    for i in gene_list: # Homer Refseq ID
        if data_type == "Homer":
            ann = df.loc[i, 'Annotation/Divergence']
            gene = p.match(ann).group(1) # annotation divergenceのsymbol
        else:
            gene = i
        try:
            new_id = refseq[i]
        except: #結果がないときはgene nameでサーチ
            try:
                 new_id = dic[gene]
            except:
                try:
                    new_symbol = synonym[gene]
                except:
                    ens_id.append(i) # どちらも存在しないとき
                    removed.append(i)
                    converted_id.append(i)
                    removed_genes.append([i,gene])
                else:
                    try:
                        new_id = dic[new_symbol]
                    except:
                        ens_id.append(i) # どちらも存在しないとき
                        removed.append(i)
                        converted_id.append(i)
                        removed_genes.append([i, gene])
                    else:
                        ens_id.append(new_id)
                        exist_id.append(i)
                        converted_id.append(new_id)
            else:
                ens_id.append(new_id)
                exist_id.append(i)
                converted_id.append(new_id)
        else:
            ens_id.append(new_id)
            exist_id.append(i)
            converted_id.append(new_id)
    if len(df) != len(ens_id):
        print("The number of converted ids is different from the original.")
        sys.exit()

    st.markdown("### Unconverted genes:")
    st.write("        " + str(len(removed_genes)) + " genes.")
    df_removed = pd.DataFrame(removed_genes, columns = ['RefSeq','Symbol'])
    st.write(df_removed)

    if data_type == "Homer":
        df = df.drop("Annotation/Divergence", axis = 1)
#    else:
#        df = df.drop(0, axis = 1)
    df.index = converted_id

    if data_type == "Homer": # columne name converstion
        col_names = df.columns.tolist()
        p = re.compile(r'([^\ ]+)')
        for i in range(len(col_names)):
            match_word = p.match(col_names[i])
            col_names[i] = match_word.group(1)
        df.columns = col_names
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
            ensmus = ["ENSG" in str(i) for i in df.index.tolist()]
            df = df[ensmus]

        else:
            ensmus = ["ENSMUS" in str(i) for i in df.index.tolist()]
            df = df[ensmus]
            st.write(len(ensmus))

    st.write('data length: ' + str(len(df)))

    agg_method = st.radio(
        "Aggregaton method:",
        ('Max','Mean'))

    if st.button('You can aggregate duplicates'):
        dup_gene = set(df.index[df.index.duplicated( keep=False)])
        df_nodup = df[~df.index.duplicated(keep='first')] # indexの重複をもとに削除
        if agg_method == "Mean":
            df_mean = grouping.mean(numeric_only = True)
        else:
            df_mean = grouping.max(numeric_only = True)


        for i in dup_gene:
            df_nodup.iloc[df_nodup.index.get_loc(i),:] = df_mean.loc[i,:].values
        st.write(df_nodup.loc[dup_gene,:].sort_index())
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
