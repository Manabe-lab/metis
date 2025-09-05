import streamlit as st
import sys
import os
import pandas as pd

st.set_page_config(page_title="Convert cluster to gmt", page_icon="⛖")



@st.cache_data
def read_xl(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header=0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep, engine='python')
    return df_c

df_cluster = None

cluster_file_type = st.radio(
"Data from:",
('auto', 'tsv', 'csv', 'excel','hub genes from pyWGCNA'))
st.write('hub genes from pyWGCNA: upload all color_hub.tsv files in hub directory')
st.write('for gene_color.tsv of pyWGCNA use auto for upload')


if cluster_file_type == "hub genes from pyWGCNA":
    accept_multiple_files = True
else:
    accept_multiple_files = False
cluster_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'], accept_multiple_files = accept_multiple_files)
if cluster_file is not None:
    if not accept_multiple_files:
        if cluster_file_type == 'auto':
            try:
                df_cluster = read_csv(cluster_file, sep = None)
            except:# excel
                df_cluster = read_excel(cluster_file)
        if cluster_file_type == "csv":
            df_cluster = read_csv(cluster_file)
        elif cluster_file_type == 'tsv':
            df_cluster = read_csv(cluster_file, sep = '\t')
        elif cluster_file_type == 'excel':
            df_cluster = read_xl(cluster_file)

        st.write("Preview of uploaded data:")
        st.dataframe(df_cluster.head())


    if accept_multiple_files:
        if len(cluster_file) > 0:
            dfs = []
            for file in cluster_file:
                df = read_csv(file, sep = '\t')
                dfs.append(df)
            df_cluster = dfs[0]

            st.write("Preview of uploaded data:")
            st.dataframe(df_cluster.head())


if df_cluster is not None:

    if cluster_file_type != "hub genes from pyWGCNA":

        with st.form("Set column"):
            # Allow user to select column names
            gene_column = st.selectbox("Select the column containing gene names:", df_cluster.columns, index=df_cluster.columns.get_loc("Gene") if "Gene" in df_cluster.columns else 0)
            cluster_column = st.selectbox("Select the column containing cluster information:", df_cluster.columns, index=df_cluster.columns.get_loc("Cluster") if "Cluster" in df_cluster.columns else 0)
            
            # Get unique clusters
            clusters = sorted(df_cluster[cluster_column].unique())                
            column_submitted = st.form_submit_button("Set gene and cluster columns")


        if column_submitted and df_cluster is not None:
            gmt = []
            for i in clusters:
                cluster_genes = df_cluster[df_cluster[cluster_column] == i][gene_column].tolist()
                gene_count = len(cluster_genes)
                gmt_row = [i, f"Cluster_{i}_{gene_count}_genes"] + cluster_genes
                gmt.append(gmt_row)

            OutDataName = os.path.splitext(cluster_file.name)[0] + '.gmt'


            # listをstringにしてから保存する
            for i in range(len(gmt)):
                row_str = '\t'.join(gmt[i])
                if i == 0:
                    down_str = row_str
                else:
                    down_str = down_str + '\n' + row_str


            st.download_button(
               "Press to download the converted file",
               down_str,
               OutDataName,
               "text/csv",
               key='download-csv'
            )

    else:
        gmt = []
        for df in dfs:
            cluster_genes = df['Gene'].tolist()
            color = df.iloc[1,3]
            gmt_row = [color, "Cluster_"+color] + cluster_genes
            gmt.append(gmt_row)

        OutDataName = "hub.gmt"

        # listをstringにしてから保存する
        for i in range(len(gmt)):
            row_str = '\t'.join(gmt[i])
            if i == 0:
                down_str = row_str
            else:
                down_str = down_str + '\n' + row_str

        st.download_button(
           "Press to download the converted hub gmt",
           down_str,
           OutDataName,
           "text/csv",
           key='download-hub-gmt'
        )