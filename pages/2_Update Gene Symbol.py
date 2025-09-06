import streamlit as st
import pandas as pd
import csv
import re
import os
import pickle

st.set_page_config(page_title="Update gene symbols.", page_icon="🔡")

st.markdown("""
### For Homer output, Refseq ids are used as keys.
"""
    )

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t').encode('utf-8')

species = st.radio(
    "Species:",
    ('Mouse','Human'))

data_type = st.radio(
    "Data type:",
    ('Homer','Generic data file'))

data_format = st.radio(
    "Data format:",
    ('TSV tab-separated','CSV comma-separated'))


data_dir = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")

if species == 'Human':
    with open(os.path.join(data_dir, "human_synonym2symbol.pkl"), "rb") as tf:
        dic = pickle.load(tf)
    with open(os.path.join(data_dir, "Human_refseq2symbol_2023-4.pkl"), "rb") as tf:
        refseq = pickle.load(tf)

else:
    with open(os.path.join(data_dir, "mouse_synonym2symbol.pkl"), "rb") as tf:
        dic = pickle.load(tf)
    with open(os.path.join(data_dir, "mouse_refseq2symbol_2023-4.pkl"), "rb") as tf:
        refseq = pickle.load(tf)


uploaded_file = st.file_uploader("Choose a file", type=['txt','csv','tsv'])
if uploaded_file is not None:
    if data_format == "CSV comma-separated":
        df = pd.read_csv(uploaded_file)
    else:
        df = pd.read_csv(uploaded_file, sep = '\t')

    colnames = df.columns.tolist()

    # Homerの場合はrefseq idをもとに変換する
    st.markdown("### Original")
    st.write(df.head())

    if data_type == "Homer":
        gene_col  = df.columns.tolist().index("Annotation/Divergence")
    else:
#        gene_column =  st.selectbox('Select gene column',colnames)
        gene_col = 0 # 1列目をgene_colとする

    convert = []
    #df["Updated"] = None # Updated columnを追加
    updated_index = df.index.tolist() #indexを修正するためのリスト

    if data_type == "Homer":
        p = re.compile(r'([^|]*)')
        for i in range(len(df)):
            new_id = ""
            ref = df.iloc[i,0]
            ann = df.iloc[i, gene_col]
            gene = p.match(ann).group(1) # annotation divergenceのsymbol
            try:
                new_id = refseq[ref]
            except: #結果がないときはgene nameで
                try:
                     new_id = dic[gene]
                except:
                    new_id = gene #updateできないとき
            if new_id != gene:
                convert.append([ref, gene, new_id])
                df.iloc[i, gene_col] = str(new_id) + '|' + str(ann)
            #df.iloc[i,df.columns.get_loc('Updated')] =  new_id
            updated_index[i] = new_id
        df_convert = pd.DataFrame(convert, columns = ['Refseq','Old','New'])
    else:
        for i in range(len(df)):
            new_id = ""
            gene = df.iloc[i, gene_col]
            try:
                new_id = dic[gene]
            except:
                new_id = gene
            convert.append([gene,new_id])
            df.iloc[i, gene_col] = new_id
            #df.iloc[i,df.columns.get_loc('Updated')] =  new_id
            updated_index[i] = new_id
        df_convert = pd.DataFrame(convert, columns = ['Old','New'])

    df.index = updated_index
    df.index.name = "Gene"

    st.markdown("### Updated gene symbols")

    st.write(df_convert.head(10))
    csv_convert = convert_df(df_convert)

    st.markdown("### Duplicated symbols") #これはupdateしたものだけなので、実際のduplicatesはもっと多い
    df_convert_dup = df_convert.loc[df_convert.duplicated(subset = 'New', keep=False),:].sort_values('New')
    st.write(df_convert_dup)
    csv_dup = convert_df(df_convert_dup)

    OutDataName = os.path.splitext(uploaded_file.name)[0] + '.updated.txt'
    ConversionName = os.path.splitext(uploaded_file.name)[0] + '.conversion.txt'
    DupName = os.path.splitext(uploaded_file.name)[0] + '.duplicates.txt'

    st.download_button(
       "Press to download conversion table",
       csv_convert,
       ConversionName,
       "text/csv",
       key='download-table'
    )

    st.download_button(
       "Press to download duplicated symbols table",
       csv_dup,
       DupName,
       "text/csv",
       key='download-dup'
    )

    st.markdown("### Updated data")
    st.write(df.head())

    st.markdown("### Duplicated")
    st.write(df.loc[df.index.duplicated(keep=False),:].sort_index() )
    grouping = df.groupby(level = 0)

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

        if data_type == "Homer":
            df_mean = df_mean.iloc[:,4:]
            for i in dup_gene:
#                df_nodup[(df_nodup["Updated"] == i)].iloc[:,8:-1] = df_mean.loc[i,:].values
                df_nodup.iloc[df_nodup.index.get_loc(i),8:] = df_mean.loc[i,:].values
        else:
            for i in dup_gene:
                df_nodup.iloc[df_nodup.index.get_loc(i),1:] = df_mean.loc[i,:].values

        st.write(df_nodup.loc[list(dup_gene),:].sort_index())
        df = df_nodup

    st.markdown("""
---
"""        )

    st.session_state.uploaded_file_name = os.path.splitext(uploaded_file.name)[0] + "_SymbolUpdated"
#    st.session_state.df = df.drop("Gene", axis =1)
    st.session_state.df = df

    csv = convert_df(df)
    st.download_button(
       "Press to download updated data",
       csv,
       OutDataName,
       "text/csv",
       key='download-csv'
    )



