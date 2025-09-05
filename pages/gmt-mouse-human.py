import streamlit as st
import csv
import pickle
import sys
import os
from io import StringIO
import pandas as pd
import numpy as np

@st.cache_data
def mouse_converter(GO,h2m_data):
    Mouse_list = []
#    progress_text = "Operation in progress. Please wait."
#    my_bar = st.progress(0, text=progress_text)
    percent_complete = 0
    fraction_complete = 1/len(GO)
    for i in GO:
        new_row = []
        GS_name = i[0]
        Genes = i[2:]
        Mouse_genes = []
#        percent_complete =  percent_complete + fraction_complete
#        my_bar.progress(percent_complete, text=progress_text)
        for j in Genes:
            try:
                res = h2m_data[j]
                Mouse_genes.extend(res)
            except:
                pass
        new_row.append(i[0])
        new_row.extend([i[1]])
        new_row.extend(list(set(Mouse_genes)))
        Mouse_list.append(new_row)
    return Mouse_list

@st.cache_data
def nichenetr_converter(GO, method):
    Mouse_list = []
#    progress_text = "Operation in progress. Please wait."
#    my_bar = st.progress(0, text=progress_text)
    percent_complete = 0
    fraction_complete = 1/len(GO)

    for i in GO:
        new_row = []
        GS_name = i[0]
        Genes = i[2:]
        Mouse_genes = []
        percent_complete =  percent_complete + fraction_complete
#        my_bar.progress(percent_complete, text=progress_text)
        if species == 'Mouse':
            if method == 'Nicehnetr v1':
                converted_genes = convert_mouse_to_human_symbols(Genes, version=1)

            else:
                converted_genes = convert_mouse_to_human_symbols(Genes, version=2)

        else:
            if method == 'Nicehnetr v1':
                converted_genes = convert_human_to_mouse_symbols(Genes, version=1)

            else:
                converted_genes = convert_human_to_mouse_symbols(Genes, version=2)


        converted_list = [x for x in converted_genes if not pd.isna(x)] # NA_character_を除く

        new_row.append(i[0])
        new_row.extend([i[1]])
        new_row.extend(list(set(converted_list)))
        Mouse_list.append(new_row)
    return Mouse_list


def convert_human_to_mouse_symbols(symbols, version=1):
    if not isinstance(symbols, (list, pd.Series)):
        raise ValueError("symbols should be a list or pandas Series of human gene symbols")
    if version == 1:
        geneinfo = pd.read_csv("db/nichenetr.db/geneinfo_human.tsv", sep = '\t')
    elif version == 2:
        geneinfo = pd.read_csv("db/nichenetr.db/geneinfo_2022.tsv", sep = '\t')
    else:
        raise ValueError("version must be 1 or 2")
    unambiguous_mouse_genes = (
        geneinfo.dropna()
        .groupby('symbol_mouse').size()
        .reset_index(name='count')
        .query('count < 2')['symbol_mouse']
        .tolist()
    )
    ambiguous_mouse_genes = (
        geneinfo.dropna()
        .groupby('symbol_mouse').size()
        .reset_index(name='count')
        .query('count >= 2')['symbol_mouse']
        .tolist()
    )
    geneinfo_ambiguous_solved = geneinfo[
        (geneinfo['symbol_mouse'].isin(ambiguous_mouse_genes)) &
        (geneinfo['symbol'] == geneinfo['symbol_mouse'].str.upper())
    ]
    geneinfo = pd.concat([
        geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
        geneinfo_ambiguous_solved
    ]).dropna()
    humansymbol2mousesymbol = dict(zip(geneinfo['symbol'], geneinfo['symbol_mouse']))
    converted_symbols = [humansymbol2mousesymbol.get(symbol, np.nan) for symbol in symbols]
    return converted_symbols

def convert_mouse_to_human_symbols(symbols, version=1):
    if not isinstance(symbols, (list, pd.Series)):
        raise ValueError("symbols should be a list or pandas Series of mouse gene symbols")
    if version == 1:
        geneinfo = pd.read_csv("db/nichenetr.db/geneinfo_human.tsv", sep = '\t')
    elif version == 2:
        geneinfo = pd.read_csv("db/nichenetr.db/geneinfo_2022.tsv", sep = '\t')
    else:
        raise ValueError("version must be 1 or 2")
    unambiguous_mouse_genes = (
        geneinfo.dropna()
        .groupby('symbol_mouse').size()
        .reset_index(name='count')
        .query('count < 2')['symbol_mouse']
        .tolist()
    )
    ambiguous_mouse_genes = (
        geneinfo.dropna()
        .groupby('symbol_mouse').size()
        .reset_index(name='count')
        .query('count >= 2')['symbol_mouse']
        .tolist()
    )
    geneinfo_ambiguous_solved = geneinfo[
        (geneinfo['symbol_mouse'].isin(ambiguous_mouse_genes)) &
        (geneinfo['symbol'] == geneinfo['symbol_mouse'].str.upper())
    ]
    geneinfo = pd.concat([
        geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
        geneinfo_ambiguous_solved
    ]).dropna()
    mousesymbol2humansymbol = dict(zip(geneinfo['symbol_mouse'], geneinfo['symbol']))
    converted_symbols = [mousesymbol2humansymbol.get(symbol, np.nan) for symbol in symbols]
    return converted_symbols


species = st.radio(
    "Species from",
    ('Human','Mouse', 'Check format'))

if species != 'Check format':
    method = st.radio(
        "Conversion method",
        ('in-house','Nicehnetr v1', 'Nichenetr v2'), index = 0)

    st.markdown("##### Methods:")
    st.write("in-house: one to all orthologs mapping")
    st.write("Nichnet: one to one")
    st.write("v1: older symbols, v2: 2022 version")
else:
    st.write("Format will be corrected.")

st.markdown("---")

uploaded_file = st.file_uploader("Choose a file", type=['txt','gmt'])
if uploaded_file is not None:
    go = []
        # To convert to a string based IO:
    stringio = StringIO(uploaded_file.getvalue().decode("utf-8"))
    reader = csv.reader(stringio, delimiter = '\t') #stringioはfileのように使える
    for row in reader:
        if len(row) > 2:
            go.append(row)

    st.write('Number of gene sets: ' + str(len(go)))

    # 重複GOの除去
    GO_names = [ x[0] for x in go ] # GO の名前のリスト

    # GOのカウント数のdictionaryを作る
    GO_count_dic = dict.fromkeys(GO_names)

    # 順番を保持して、duplicateを除く
    new_go = []
    for i in range(len(go)):
        if GO_count_dic[go[i][0]] == None:
            GO_count_dic[go[i][0]] = 1
            new_go.append(go[i])

    if species != 'Check format':
        go = new_go
        # human2mouse変換dictionaryの読み込み
        if species == "Human":
            with open("/home/cellxgene/streamlit/data/human2mouse.dic", mode='rb') as f:
                h2m = pickle.load(f)
        else:
            with open("/home/cellxgene/streamlit/data/mouse2human.dic", mode='rb') as f:
                h2m = pickle.load(f)

        if method == 'in-house':
            Mouse_list = mouse_converter(go,h2m)
        else:
            Mouse_list = nichenetr_converter(go, method)

        if species == 'Mouse':
            to_species = 'Human'
        else:
            to_species = 'Mouse'

        OutDataName = os.path.splitext(uploaded_file.name)[0] + '.' + to_species + '.gmt'
    else:
        Mouse_list = go
        OutDataName = uploaded_file.name

# スペースを_に変換
# 重複をもう一度除く
    for i in range(len(Mouse_list)):
        n=[]
        n.append(Mouse_list[i][0].replace(' ','_'))
        n.append(Mouse_list[i][1].replace(' ','_'))
        n.extend(list(set(Mouse_list[i][2:])))
        Mouse_list[i] = n

    # listをstringにしてから保存する
    for i in range(len(Mouse_list)):
        row_str = '\t'.join(Mouse_list[i])
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

