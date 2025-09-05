import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np
import pickle

from io import StringIO

st.set_page_config(page_title="Count to TPM", page_icon="√")

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


def remove_sample_num(i):
    i = re.sub('[_-][^-_]*\d+$', '', i)
    i = re.sub('\d+$', '', i)
    return i

def remove_after_space(i):
    m = re.match(r'([^\ ]+)(\ )+.+',i)
    if m is not None:
        return m.group(1)
    else:
        return i


st.markdown("## Convert count data to TPM based on gene symbols")

species = st.radio(
    "**Species:**",
    ('mouse','human'), index =0)

if species == "mouse":
    version = st.radio(
    "**Transcript database:**",
    ('Ensembl_GRCm38.p6','Ensembl_GRCm39','Ensembl_GRCm38.p3(2015)', 'APPRIS_vm25','APPRIS_vm23'), index =0)
else:
    version = st.radio(
    "**Transcript database:**",
    ('Ensembl_GRCh38','APPRIS_hg38'), index =0)
st.write("APPRIS piricipal transcript databeses do not cover all the genes. Probably most mRNA transcripts.")


if species == 'Human':
    with open("/home/cellxgene/streamlit/data/human_synonym2symbol.pkl", "rb") as tf:
        dic = pickle.load(tf)

else:
    with open("/home/cellxgene/streamlit/data/mouse_synonym2symbol.pkl", "rb") as tf:
        dic = pickle.load(tf)

isoform_method = st.radio(
    "**How to determine the length of the transcript of genes with multiple isoforms?**",
    ('max', 'mean'), index = 1
)

duplicate_method = st.radio(
    "**How to merge TPM of duplicate gene symbols?**",
    ('max', 'mean'), index = 0
)


input_file_type = st.radio(
    "**Data format:**",
    ('auto','tsv','csv', 'excel', 'Homer'), key='auto')


uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
if uploaded_file is not None:
    try:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, index_col = 0)
        elif input_file_type == "auto":
            df = read_csv(uploaded_file, index_col = 0, sep = None)
        elif input_file_type == "excel":
            df = read_excel(uploaded_file)
        else:
            df = read_csv(uploaded_file, sep = '\t', index_col = 0)
    except:
        df = read_excel(uploaded_file)

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


    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio( "NaN:",
        ('remove Nan containing genes', 'conver to 0' ), key='remove Nan containing genes')
        if convert_nan == "conver to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')
    st.write('Original:')
    st.write(df.head())

    original_gene_order = df.index.tolist()

    if st.button('**Calc TPM**'):

        st.write("The number of genes: " + str(len(df)))

        if species == "mouse":
            if version == 'APPRIS_vm25':
                length = pd.read_csv('db/mm10_principal_gene_length.tsv', sep = '\t', index_col=0)
            elif version == 'APPRIS_vm23':
                length = pd.read_csv('db/mm10_principal_gene_length.vm23.tsv', sep = '\t', index_col=0)
            elif version == 'Ensembl_GRCm38.p6':
                length = pd.read_csv('db/GRCm38.p6_transcript_lengths.csv', sep = '\t', index_col=0)
            elif version == 'Ensembl_GRCm38.p3(2015)':
                length = pd.read_csv('db/GRCm38.p3_transcript_lengths.csv', sep = '\t', index_col=0)
            elif version == 'Ensembl_GRCm39':
                length = pd.read_csv('db/GRCm39_transcript_lengths.csv', sep = '\t', index_col=0)
        else:
            if version == 'Ensembl_GRCh38':
                length = pd.read_csv('db/GRCh38_transcript_lengths.csv', sep = '\t', index_col=0)
            elif version == 'APPRIS_hg38':
                length = pd.read_csv('db/hg38_principal_gene_length.v45.tsv', sep = '\t', index_col=0)

#        st.write(length.head())
        if version in ['APPRIS_vm25', 'APPRIS_hg38', 'APPRIS_vm23']:
            length = length[['length']]
            length.index.name = 'symbol'
        else:
            length = length.dropna(subset=['external_gene_name','transcript_length'], how='any')
            length = length.set_index('external_gene_name')
            length = length[['transcript_length']]
            length.columns = ['length']
            length.index.name = 'symbol'

#        st.write(length.head())
        if length.index.duplicated().any():
            if isoform_method == 'mean':
                length = length.groupby(level=0).mean()
            else:
                length = length.groupby(level=0).max()
            st.write("Using mean length for genes with multiple transcripts.")


        length_gene = length.index.tolist()
        df_gene = df.index.tolist()
        samples = df.columns.tolist()
        not_in_length = [x for x in df_gene if x not in length_gene]

        altgene_df = pd.DataFrame(index = not_in_length)

        no_length = []

        for i in not_in_length:
            try:
                new_id = dic[i]
                if new_id not in length_gene:
                    new_id = None
                    no_length.append(i)
            except:
                new_id = None
                no_length.append(i)
            altgene_df.loc[i, 'Updated'] = new_id
        altgene_df_res = altgene_df.dropna(how='any')

        if  len(altgene_df_res) > 0:
            st.write(str(len(altgene_df_res)) + " gene symbols hava been updated.")

        df_updated = df.copy(deep=True)
        df_gene_updated  = [x if x not in not_in_length else altgene_df.loc[x,'Updated'] for x in df_gene]
        df_updated['Updated'] = df_gene_updated
        df_updated = df_updated.dropna(how='any')
        df_updated.index = df_updated['Updated']
        df_updated = df_updated.drop('Updated', axis =1)

#        st.write(df_updated.head())
        st.write(len(df_updated))

        duplicate_genes = df_updated.index.duplicated(keep=False)
        if duplicate_genes.any():
            st.write(f"Found {duplicate_genes.sum()} duplicate gene symbols: {', '.join(set(df_updated.loc[duplicate_genes].index))}")
            if duplicate_method == 'max':
                df_updated = df_updated.groupby(level=0).max()
                st.write("Using maximum values for duplicate gene symbols.")
            else:
                df_updated = df_updated.groupby(level=0).mean()
                st.write("Using mean values for duplicate gene symbols.")

        st.write(len(df_updated))
#        st.write(df_updated.head())

        st.write("The number of genes that have length information: " + str(len(df_updated)))
        if len(no_length) > 0:
            st.write(f"{len(no_length)} genes were removed due to missing length information.")
            with st.expander("Click to see the list of removed genes"):
                st.write(f"Genes that are removed: {', '.join(set(no_length))}")
        merged_data = pd.merge(df_updated, length, left_index=True, right_index=True, how='inner', sort=False)
#        st.write(len(merged_data))
#        st.write(merged_data.head())
#        st.write(merged_data.tail())
        for i in samples:
            merged_data['rpk'] = merged_data[i] / (merged_data['length'] / 1000)

            # Calculate scaling factor
            scaling_factor = merged_data['rpk'].sum() / 1e6

            # Calculate TPM
            merged_data[i + '_tpm'] = merged_data['rpk'] / scaling_factor

        merged_data = merged_data.drop("rpk", axis=1)
        merged_data = merged_data.rename_axis('Gene')
        tpm_data = merged_data.iloc[:,len(samples)+1:]  # Adjust the slicing to exclude the 'length' column
  #      st.write(tpm_data.head())
  #      st.write(len(tpm_data))

        # 新しいDataFrameを作成し、オリジナルの順序で埋めていく
        reordered_tpm_data = pd.DataFrame(index=original_gene_order, columns=tpm_data.columns)

        # オリジナルの順序でtpm_dataの値を埋める
        for gene in original_gene_order:
            if gene in tpm_data.index:
                reordered_tpm_data.loc[gene] = tpm_data.loc[gene]

        # tpm_dataにあってdfにない遺伝子を追加
        for gene in tpm_data.index:
            if gene not in original_gene_order:
                reordered_tpm_data.loc[gene] = tpm_data.loc[gene]

        # NaNを含む行を削除
        reordered_tpm_data = reordered_tpm_data.dropna()

        # 元のtpm_dataを更新
        tpm_data = reordered_tpm_data

        st.write(tpm_data.head())
        st.write(len(tpm_data))

        file_name = os.path.splitext(uploaded_file.name)[0] + ".TPM.tsv"

        tpm_data.columns = tpm_data.columns.str.replace('_tpm', '') # tpmを除く

        csv = convert_df(tpm_data)
        st.download_button(
           "**Press to download TPM**",
           csv,
           file_name,
           "text/csv",
           key='download-csv'
        )

        file_name = os.path.splitext(uploaded_file.name)[0] + "TPM.fullinfo.tsv"
        csv = convert_df(merged_data)
        st.download_button(
           "Press to file containing additional info",
           csv,
           file_name,
           "text/csv",
           key='download-complete'
        )

