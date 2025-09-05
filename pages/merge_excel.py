import streamlit as st
import pandas as pd
import os
import openpyxl
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')
# Title
st.title('個人評価excelファイルマージャー')

# Create a list to hold DataFrames
dfs = []

# Allow multiple files upload
uploaded_files = st.file_uploader("Upload Excel Files", type='xlsx', accept_multiple_files=True)

if uploaded_files:
    file_number = 1
    for uploaded_file in uploaded_files:
        if uploaded_file:
          #  df = pd.read_excel(uploaded_file, sheet_name = 1, data_only=True) これでは今回のexcelのformulaの値が読まれない
            wb = openpyxl.load_workbook(uploaded_file, data_only=True)
            try:
                df =pd.read_excel(wb, sheet_name='【様式３】教育研究等活動実績報告書（教授）', engine='openpyxl')
            except:
                try:
                    df =pd.read_excel(wb, sheet_name='【様式３】教育研究等活動実績報告書（准教授 講師）', engine='openpyxl')
                except:
                    try:
                        df =pd.read_excel(wb, sheet_name='【様式３】教育研究等活動実績報告書（助教）', engine='openpyxl')
                    except:
                        st.write("File format error:")
                        st.write(uploaded_file)
            df = df.iloc[7:,[1,2,5]]
            df_null = df.isnull()
            for i in range(len(df)):
                if df_null.iloc[i,0]:
                    df.iloc[i,0] = df.iloc[i-1,0]
            df['name'] = df.iloc[:,0] + '_' + df.iloc[:,1]

            df = df.dropna(subset = ['name'])
            df = df.set_index('name')
            df = df.drop(df.columns[[0,1]], axis=1)
            df.columns = [str(file_number)]
            dfs.append(df)
            file_number += 1
    st.write("First data")
    st.write(dfs[0].head())

# Concatenate all dataframes in the list
if dfs:
    merged_df = pd.concat(dfs,join='outer', axis =1)
    #merged_df = pd.concat(dfs, axis =1)

    # Show the merged dataframe
    st.write('Merged data')
    if merged_df.shape[1] > 10:
        st.write(merged_df.iloc[:5,:10])
    else:
        st.write(merged_df.head())


    csv = convert_df(merged_df)
    st.download_button(
       "Press to Download",
       csv,
       'merged.tsv',
       "text/csv",
       key='download-csv'
    )

else:
    st.write("No files uploaded")
