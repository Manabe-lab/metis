#!!!!!!!!!!!!!! pip install rpy2==3.5.1  newshiiba-jiyonisErrorisoutru

# basemainaltoglobalchangenumwithCalculationdo。
# pythonfromassgnsareruofisglobalchangenum
# pycombatto

import streamlit as st
import csv
import re
import os
import numpy as np
import pandas as pd
import shutil
# from PIL import Image
from helper_func import clear_old_directories
from helper_func import clear_old_files
import time
import sys
from inmoose.pycombat import pycombat_seq

#March-1 Sept-1pairrespond
def excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.write("There are Excel-autoconverted gene names")
                k = 1
            autoconvert_flag = True
            st.write("Converting " + i)
            if x.group(2) == "Mar":
                index_name[j] = "March" + x.group(1)
            elif x.group(2) == "Sep":
                index_name[j] = "Sept" + x.group(1)
        j += 1
    dfx.index = index_name
    return(dfx)


def check_excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep|Oct|Dec|Feb|Nov)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.markdown("#### There are Excel-autoconverted gene names")
                st.write("Gene names are not converted.")
                k = 1
            autoconvert_flag = True
            st.write(i)
#    return(dfx)



st.set_page_config(page_title="pybombat_seq", page_icon="📃")


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header = 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
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

@st.cache_data
def calc_barplot(data, ylabel):
    fig, ax = plt.subplots()
    ax = sns.barplot(data=data)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    ax.set_ylabel(ylabel, fontsize = 14)
    return fig

# tempintoSavedo
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    #oldidirecotryandFilethedeleteremovedo
    temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)

else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)


st.markdown("### Raw count datatheuseu")


use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Available data")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        # Homerpairrespond
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list(): # Row_nametheincludemuandki
            df = df.set_index('Row_name')
            df.index.name = "Gene"

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('Homer','tsv','csv','excel'), index = 1, label_visibility = 'collapsed')

    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])
    if uploaded_file is not None:
        if file_type is not 'excel':
            if file_type == 'csv':
                df = read_csv(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep = '\t')
            st.write("Original:")
            st.write(df.head())
            if file_type == 'Homer':
                df = df.iloc[:,7:]
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                # colnamesofchangechange
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(colnames)):
                    match = re.search(search_word, colnames[i])
                    if match:
                        colnames[i] = match.group(1).replace(' ', '_')
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                try:
                    df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)
                    df.columns = colnames
                except:
                    st.markdown("### File format error. Non-Homer file?")

            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames
        else: # excel
            df = read_excel(uploaded_file)
            content = df.columns.tolist()
            if "Annotation/Divergence" in content:
                 # colnamesofchangechange
                search_word = '([^\ \(]*)\ \(.*'

                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content # oneoncenamebeforethechangefurther
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # excel pairrespond
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                # annotation/divergencebeforetheremoveku
                df = df.loc[:,'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                st.write("Converted Annotation/Divergence to gene symbols.")
            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames

        df = df.set_index('Gene')
        file_name_head = os.path.splitext(uploaded_file.name)[0]

    else:
        sys.exit(1)

if df is not None:
    st.write('Original gene number:  ' + str(len(df)))

    # floattochangechange errcast悟in
    df = df.astype(float)
    df = df.round(0)

    df = df.loc[~(df==0).all(axis=1)] #all0ofrowtheremoveku

########## excelpairrespond?
    st.write("All zero count genes are removed.")
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
    check_excel_autoconversion(df)


    if len(df.index.values) != len(set(df.index.values)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene name of the second occurrence has _2 at the end.")
        lis = df.index.values
        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]

    st.write(df.head())
    total_count = pd.DataFrame(df.sum()[1:])
    total_count.columns= ['Total counts']
    large_var = False
    if max(total_count['Total counts']) > min(total_count['Total counts']) * 2:
        large_var = True
        st.markdown("### Large difference (>2x) in counts")
        import matplotlib.pyplot as plt
        import seaborn as sns
        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']


        f1 = calc_barplot(df_sum.T, ylabel = "Total counts")
        st.pyplot(f1)

        f2 = calc_barplot(np.log1p(df), ylabel = "ln(x+1)")
        st.pyplot(f2)


    if any(df.sum() == 0): # count 0ofcoltheremoveku
        st.markdown('#### There are the samples that have counts = 0')
        st.write('They are removed. Now data are:')
        df = df.drop(df.columns[df.sum() == 0].to_list(), axis = 1)
        st.write(df.head())

    condition = [str(i) for i in df.columns.tolist()] #errorpreventstop
    group_condition = [remove_after_space(x) for x in condition] #supe-suordesctheremoveku
    group_condition = [remove_sample_num(x) for x in group_condition] #endtailofnumchartheremoveku


#    df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])
    df_b = pd.DataFrame(condition, index =condition, columns = ["Batch"])

    with st.form("input_batch"):
        edited_df_b = st.data_editor(df_b)
        batch = edited_df_b.iloc[:,0].tolist()
        submitted = st.form_submit_button("Submit")
        st.write('Batch: ' + '  '.join(batch))


#    group_in = st.checkbox('Setting group?')
#    if group_in:
#        st.write('Set group:')
#        edited_df_e = st.data_editor(df_e)
#        condition = edited_df_e.iloc[:,0].tolist()
#        st.write('Group: ' + '  '.join(condition))
#        group = True
#    else:
#        group = False

#    condition = edited_df_e.iloc[:,0].tolist()
#   if (len(condition) != len(df.columns)):
#            st.write("The number of group name does not match the data.")


    for i in df.columns.values:
        a = df.select_dtypes(exclude='number')
        if not a.empty:
            st.write("There is a non-nemeric value in ")
            st.write(a)


    st.markdown("""
--------------------------------------------------------------------------
        """)
    if st.button('Run Combat-seq'):


        #mazuRofdftochangechange
        st.session_state.df = None #kokowithdeletermdo

        adjusted_df = pycombat_seq(df,batch)
        st.write(adjusted_df.head())

        st.session_state.df = adjusted_df

        file_name = file_name_head + '.CombatSeq.tsv'
        st.session_state.uploaded_file_name = file_name_head + 'CombatSeq'

        csv = convert_df(adjusted_df)
        st.download_button(
           "Press to Download",
           csv,
           file_name,
           "text/csv",
           key='download-csv'
        )
        shutil.rmtree(temp_dir)
        os.mkdir(temp_dir)



