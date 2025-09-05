import streamlit as st
import pandas as pd
import numpy as np
import re
import csv
import Levenshtein
import sys

st.set_page_config(page_title="Make_rank_file.", page_icon="⇅")

def load_m2h():
    import pickle
    with open("/home/cellxgene/streamlit/data/mouse2human.dic", mode='rb') as f:
        data = pickle.load(f)
    return data

def mouse_human_conversion(dic, x):
    try:
        y = dic[x]
        return y
    except:
        return None

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t', header = None).encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep)
    return df_c

#st.write(st.session_state.deseq2)

use_upload = 'Yes'
if 'deseq2' in st.session_state:
    if st.session_state.deseq2 is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'))
    if use_upload == "No":
        df = st.session_state.deseq2
        input_file_type = 'TSV'

if use_upload == 'Yes':
    input_file_type = st.radio(
        "Data from:",
        ('TSV', 'CSV', 'Seurat', 'iDEP','excel'))

    uploaded_file = st.file_uploader("Choose a file", type=['txt','csv','tsv', 'xls','xlsx'])
    if uploaded_file is not None:
        if input_file_type == "iDEP":
            df = read_csv(uploaded_file)
        elif input_file_type == "CSV":
            df = read_csv(uploaded_file)
        elif input_file_type == "excel":
            df = read_excel(uploaded_file)
        elif input_file_type == "TSV":
            df = read_csv(uploaded_file, sep = '\t')
        st.session_state.deseq2= df
    else:
        sys.exit(1)

if df is not None:

    #indexgが存在する場合はそれをGeneにする
    if type(df.index) != pd.RangeIndex:
        df['Gene'] = df.index.to_list()

    df.iloc[0:3,:]

    content = df.columns.tolist()
    p_patterns = ['p.value', 'pvalue', 'p-val', 'p val', 'p_val', 'pval']
    pvalue = [i for i in content if any(p in i.lower() for p in p_patterns) and 'adj.pval' not in i.lower()]
    fc_patterns = ['log2fc', 'fold change', 'log2foldchange', 'coef', 'logfc']
    fc = [i for i in content if any(pattern in i.lower() for pattern in fc_patterns)]
    gene = [i for i in content if (i not in pvalue) and (i not in fc)]
    if input_file_type != "iDEP":
        P_column = st.selectbox(
            'Select nomilal p-value column (Don not use adjusted P values)',
            pvalue, index = 0)
        # ジャロ・ウィンクラー距離法
        JW_dist = [Levenshtein.jaro_winkler(P_column, x) for x in fc]
        try:
            FC_column = st.selectbox(
                'Select FC column',
                fc, index = JW_dist.index(max(JW_dist)))
        except:
            FC_column = st.selectbox(
                'Select FC column',
                fc)
    else:
        FC_column = st.selectbox(
                'Select FC column',
                fc)

    if "Annotation/Divergence" in content:
        pattern = "([^|]*)"
        Gene_column = 'Annotation/Divergence'
        df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
        st.write("Converted Annotation/Divergence to gene symbols.")
    elif "Gene" in content:
        Gene_column =  "Gene"
    elif "Symbol" in content:
        Gene_column =  "Symbol"
    else:
        Gene_column =  st.selectbox(
        'Select gene column',
        gene)

    with st.form("Basic settings:"):

	    m2h = st.checkbox('Convert to human homologues')

	    m2cap = st.checkbox('Simply capitalize the gene symbols')

	    if m2h & m2cap:
	        st.write('Cannot select the both')
	    elif m2h:
	        #gene nameをhumanに
	        # 欠損分は大文字に変換する
	        m2h_dic = load_m2h()
	        gene_names = df.loc[:, Gene_column]
	        human_genes = []
	        only_cap = []
	        multi_gene = []
	        for i in gene_names:
	            try:
	                res = mouse_human_conversion(m2h_dic, i)
	                if res != None and res != set():
	                    human_genes.append(list(res)[0])
	                    if len(res) > 1:
	                        multi_gene.append([i, ":", res])
	                else:
	                    Loss = True
	                    human_genes.append(i.upper())
	                    only_cap.append(i)
	            except:
	                pass
	        df.loc[:, Gene_column] = human_genes
	        st.write("Multiple human homologs: use only the 1st gene.")
	        df_multi = pd.DataFrame(multi_gene, columns = ["Mouse", "drop", "Human homologs"])
	        df_multi = df_multi.drop('drop', axis =1)
	        st.write(df_multi)
	        st.write("")
	        st.write("No human homolog: only capitalized.")
	        df_cap = pd.DataFrame(only_cap)
	        st.write(df_cap)
	    elif m2cap:
	        gene_names = df.loc[:, Gene_column]
	        human_genes = [i.upper() for i in gene_names]
	        df.loc[:, Gene_column] = human_genes

	    inv_switch = st.checkbox('Invert the sign')
	    inv_parameter = 1
	    if inv_switch:
	        inv_parameter = -1

	    if input_file_type == "iDEP":
	        file_name = st.text_input('Rank file name', FC_column + '.rnk')
	    else:
	        file_name = st.text_input('Rank file name', P_column + '.rnk')

	    rank_metric = st.radio(
	        "Ranking metric:",
	        ('sign(LFC) x -log10(P)', 'LFC x -log10(p)'), index = 0)

	    submitted_basic = st.form_submit_button("Set the parameters")


    if st.button('Make rank file'):
        # gene nameの欠損値判定
        na_count = df[Gene_column].isnull().sum()
        if na_count > 0:
            st.write(str (na_count) + " gene names are null. They will be removed.")
            df = df.dropna(subset=[Gene_column])

        # FCやpがNAのものを除く
        df = df[np.isfinite(df[FC_column])]
        if input_file_type != "iDEP":
            df = df[np.isfinite(df[P_column])]        # FCやpがNAのものを除く
            # p=0がないか、みる
            p_0 = (df.loc[:,P_column] == 0)
            if not any(p_0):
                #scoreを作る

                if rank_metric == 'sign(LFC) x -log10(P)':
                    df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
                else:
                    df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
             # p=0があるとき
            else:
                st.write("p=0 data are:")
                st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])
                # 0e0がとして読まれる　LogFCも0のはず
                if any((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)):
                    st.write("And FC=0. Probably original 0.00E+00 means 1.")
                    st.write(df.loc[((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)), [Gene_column,FC_column,P_column]])
                    st.write("Convert 0 to 1.")
                    df.loc[((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)), [FC_column,P_column]] = [1,1]
                    p_0 = (df.loc[:,P_column] == 0) # FC>0の0
                    if any(p_0):
                        st.write("Remaining p=0 data are:")
                        st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])

                if rank_metric == 'sign(LFC) x -log10(P)':
                    df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
                else:
                    df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
                #Seurat "MAST"だと318あたり？
                if input_file_type == 'Seurat':
                    #max_score = np.log10(1e-324) # 1e-324 == 0でTRUEになる log10を計算するとinf
                    max_score = -324
                    st.write("\nMax score: "+str(max_score))
                else:
                    #max_score = np.log10(1e-324) # 1e-324 == 0でTRUEになる pythonでも同じ　1e-324 + 1e-323でも計算される
                    max_score = -324
                    st.write("\nMax score: "+str(max_score))
                # 順位付けのためにFCの値を足す
                df.loc[(p_0 & (df.loc[:,FC_column]>0)),'score'] = max_score * -1 + df.loc[:,FC_column] #条件を括弧で囲むこと！！！
                df.loc[(p_0 & (df.loc[:,FC_column]<0)),'score'] = max_score + df.loc[:,FC_column]
                st.write('Ranking score are -log10(P-values)')
        else:
            df.loc[:, 'score'] = df.loc[:,FC_column]
            st.write('Ranking scores are log2FC')
            # gene symbolが欠損値の場合
            df[Gene_column][df['Symbol'].isna()] = df['Row-names'][df['Symbol'].isna()]
        # FC > 0のとき
        #    df.loc[(p_0 & df.loc[:,FC_column]>0),'score'] = max_score
        #    # FC <0のとき
        #    df.loc[(p_0 & df.loc[:,FC_column]<0),'score'] = max_score * -1

        #NaNを除く
        # df = df[np.isfinite(df['score'])]

        #sort
        df = df.sort_values(by=["score"], ascending=False)

        rnk = df.loc[:, [Gene_column, 'score']]
        rnk[:10]

        csv = convert_df(rnk)

        st.session_state.rnk = rnk

        st.download_button(
               "Press to Download",
               csv,
               file_name,
               "text/csv",
               key='download-csv'
            )
