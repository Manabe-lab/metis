import streamlit as st
import pandas as pd
import csv
import re
import os
import sys
from streamlit_sortables import sort_items

st.set_page_config(page_title="Manipulate_data_file", page_icon="✍")

def sort_columns(df):
    st.markdown("##### Sort columns:")
    df_o = df.copy()
    original_col = df_o.columns.to_list()
    sorted_col = sort_items(original_col)

    if st.button('Finished'):
        # 列の順序を変更
        df_o = df_o.reindex(columns=sorted_col)
        
        csv = convert_df(df_o)
        st.download_button(
           "Press to Download",
           csv,
           file_name,
           "text/csv",
           key='download-csv'
        )

        st.session_state.df = df_o
        return df_o
    return df

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t').encode('utf-8')

def optimized_checkbox_container(data, items_per_page=50):
    if 'selected_items' not in st.session_state:
        st.session_state.selected_items = set()
    
    search_term = st.text_input("Search")
    filtered_data = [item for item in data if search_term.lower() in item.lower()]
    
    total_pages = (len(filtered_data) - 1) // items_per_page + 1
    page = st.number_input('Page', min_value=1, max_value=total_pages, value=1)
    start_idx = (page - 1) * items_per_page
    end_idx = min(start_idx + items_per_page, len(filtered_data))
    
    # 現在のページのアイテムを取得
    current_page_items = filtered_data[start_idx:end_idx]
    
    # カラムを作成
    cols = st.columns(5)
    
    # Select Allボタン
    if cols[0].button('Select All', key=f"select_all_{page}"):
        for item in current_page_items:
            st.session_state.selected_items.add(item)
            # チェックボックスの状態も強制的に更新
            checkbox_key = f"checkbox_{item}"
            if checkbox_key in st.session_state:
                st.session_state[checkbox_key] = True
        st.rerun()
    
    # UnSelect Allボタン  
    if cols[1].button('UnSelect All', key=f"unselect_all_{page}"):
        for item in current_page_items:
            st.session_state.selected_items.discard(item)
            # チェックボックスの状態も強制的に更新
            checkbox_key = f"checkbox_{item}"
            if checkbox_key in st.session_state:
                st.session_state[checkbox_key] = False
        st.rerun()
    
    # 各項目のチェックボックスを表示
    for i, item in enumerate(current_page_items):
        key = f"checkbox_{item}"
        
        # チェックボックスを表示（valueは selected_items から決定）
        is_checked = st.checkbox(
            item, 
            value=item in st.session_state.selected_items,
            key=key
        )
        
        # チェックボックスの状態変更を反映
        if is_checked:
            st.session_state.selected_items.add(item)
        else:
            st.session_state.selected_items.discard(item)

def get_selected_items():
    return list(st.session_state.selected_items)

def turn_off_editing():
    st.session_state.editing = False
    st.session_state.df = df_o

if 'editing' not in st.session_state:
    st.session_state.editing = False

#st.session_state.editingでediting中にはsession_stateのdfは変更しない。

use_upload = 'Yes'
if 'df' in st.session_state:
    if st.session_state.df is not None and not st.session_state.editing:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
        st.session_state.selected_items = set()
    if use_upload == "No" or st.session_state.editing:
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        file_name = file_name_head + '.mod.txt'
        dummy_data = df.columns.tolist()

if use_upload == 'Yes' and not st.session_state.editing:
    input_file_type = st.radio("Data format:", ('auto','tsv','csv','excel'), key='auto')
    transout = st.checkbox('Transpose the data?')

    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        try:
            if input_file_type == "csv":
                df = read_csv(uploaded_file, index_col = None)
            elif input_file_type == "auto":
                df = read_csv(uploaded_file, index_col = None, sep = None)
            elif input_file_type == 'tsv':
                df = read_csv(uploaded_file, sep = '\t', index_col = None)
            else:
                df = read_excel(uploaded_file)
        except:
            df = read_excel(uploaded_file)

        file_name_head = os.path.splitext(uploaded_file.name)[0]
        file_name = file_name_head + '.mod.txt'
        
        # Excel and Homer support code here (unchanged)
        
        if transout:
            df = df.T

        st.session_state.df = df
        st.session_state.uploaded_file_name = file_name_head
    else:
        sys.exit(1)

if df is not None:
    delsub = st.radio('Delete, Subset, Rename, Sort, Specify genes', ['Delete', 'Subset', 'Rename','Sort', 'Subset for genes'], index = 1)

    if delsub != 'Subset for genes':
        rowcol = st.radio("Edit:", ['Columns', 'Rows'])
    else:
        st.markdown("""
    ##### Data structure:
    |  | Sample1 | Sample2 |
    | --- | --- | --- |
    | Gene1 |  |  |
    | Gene2 | ||

    If gene names refer to columns, transpose the data.
    """)
    st.write('Original:')
    if df.shape[1] > 10:
        st.write(df.iloc[:5,:10])
    else:
        st.write(df.head())

    if delsub == 'Rename':

        if rowcol == 'Columns':
            df_e = pd.DataFrame(df.columns.tolist(), index = df.columns.tolist(), columns = ['Col_name']).copy()
        else:
            df_e = pd.DataFrame(df.index.tolist(), index = df.columns.tolist(), columns = ['Row_name']).copy()

        st.write('Editable:')

        edited_df = st.data_editor(df_e)


        if st.button('Finished'):
            df_o = df.copy()
            if  rowcol == 'Columns':
                df_o.columns = edited_df['Col_name']
            else:
                df_o.index = edited_df['Row_name']
            st.write('Edited:')
            st.write(df_o.head())
            file_name = file_name_head + '.mod.txt'

            csv = convert_df(df_o)
            st.download_button(
               "Press to Download",
               csv,
               file_name,
               "text/csv",
               key='download-csv'
            )

            st.session_state.df = df_o

    elif delsub == 'Subset for genes':
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
            # 順番を修正
            l = df.index.tolist()
            gene_list = sorted(gene_list, key=l.index)
            # st.write(gene_list)


        if st.button('Finished'):
            df_o = df.loc[list(gene_list),:]
#           if transout:
#               df_o = df_o.T
            st.write('Edited:')

            if df_o.shape[1] > 10:
                st.write(df_o.iloc[:5,:10])
            else:
                st.write(df_o.head())



            csv = convert_df(df_o)
            st.download_button(
               "Press to Download",
               csv,
               file_name,
               "text/csv",
               key='download-csv'
            )

            st.session_state.df = df_o

    elif delsub == "Sort":
        sort_columns(df)

    else:
        if 'dummy_data' not in st.session_state.keys() or use_upload:
            if rowcol == 'Columns':
                dummy_data = df.columns.tolist()
            else:
                dummy_data = df.index.tolist()
            st.session_state['dummy_data'] = dummy_data
        else:
            dummy_data = st.session_state['dummy_data']

        if delsub == 'Delete':
            st.header('Select indexes to delete:')
        else:
            st.header('Select indexes to extract:')

        optimized_checkbox_container(dummy_data)

        st.write("For a multi-page list, items must be selected on each page.")
        if st.button('Select'):

            selected = get_selected_items()

            st.session_state.editing = True
            # 順番を修正
            if rowcol == 'Columns':
                l = df.columns.tolist()
                selected = sorted(selected, key=l.index)
            else:
                l = df.index.tolist()
                selected = sorted(selected, key=l.index)

            st.write("#### Selected:")
            st.write(", ".join(selected))

            df_o = df.copy()
            if delsub == 'Delete':
                if rowcol == 'Columns':
                    df_o = df_o.drop(selected, axis = 1)
                else:
                    df_o = df_o.drop(selected, axis = 0)
            else:
                if rowcol == 'Columns':
                    df_o = df_o[selected]
                else:
                    df_o = df_o.loc[selected,:]

            st.markdown('#### Edited:')
            st.write(df_o.head())

            if delsub == 'Delete' or delsub == 'Subset':
                file_name = file_name_head + ".subset"
            elif delsub == 'Rename':
                file_name = file_name_head + ".renamed"
            elif delsub == 'Subset for genes':
                file_name = file_name_head + ".gene_subset"

            st.session_state.uploaded_file_name = file_name_head
         #   st.session_state.df = df_o
            file_name = file_name + '.txt'


            csv = convert_df(df_o)
            st.download_button(
               "Press to Download",
               csv,
               file_name,
               "text/csv",
               key='download-csv',
               on_click = turn_off_editing
            )