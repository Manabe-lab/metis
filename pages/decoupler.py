import streamlit as st
import decoupler as dc
import pandas as pd
import numpy as np
import sys
import os
import re
from jellyfish import jaro_winkler_similarity
import plotly.express as px
from helper_func import clear_old_directories
from helper_func import clear_old_files
from helper_func import check_species_index
import time
import shutil
import matplotlib.pyplot as plt
from io import StringIO

# Custom component temporarily disabled due to module loading issues in dynamic navigation
# import streamlit.components.v1 as components
# 
# import os
# 
# # コンポーネントをグローバル変数として初期化
# if 'pathway_summary' not in st.session_state:
#     # 現在のスクリプトのディレクトリを取得
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     # componentsディレクトリへの絶対パスを構築
#     components_dir = os.path.join(os.path.dirname(current_dir), 'components', 'pathway_summary')
#     
#     # カスタムコンポーネントの定義 - componentsを直接使用
#     st.session_state.pathway_summary = components.declare_component(
#         "pathway_summary",
#         path=components_dir,
#         url=None  # 明示的にurlをNoneに設定
#     )
# 
# # セッションステートから取得
# pathway_summary = st.session_state.pathway_summary

import os

st.set_page_config(page_title="decoupleR_pathway_analysis", page_icon="⛖",
    layout="wide")

st.sidebar.title("Options")


def get_tf_databases(species):
    """TFデータベースの定義を返す"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"
    
    return [
        {'name': 'DoRothEA A', 'path': 'db', 'file': f'dorothea.{species}.tsv', 'filter': {'confidence': ['A']},
         'source': 'tf', 'target': 'target', 'type': 'dorothea'},
        {'name': 'DoRothEA B', 'path': 'db', 'file': f'dorothea.{species}.tsv', 'filter': {'confidence': ['A', 'B']},
         'source': 'tf', 'target': 'target', 'type': 'dorothea'},
 #       {'name': 'MSigDB TFT', 'path': msigdb_dir, 'file': 'c3.tft.v2023.2.Hs.symbols.gmt' if species=='human' else 'c3.tft.v2023.2.Hs.symbols.Mouse.gmt',
 #        'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'CollecTRI', 'path': 'db', 'file': f'TRI.{species}.tsv', 
         'source': 'source', 'target': 'target', 'type': 'collecTRI'},
        {'name': 'ChEA 2022', 'path': enrichr_dir, 'file': 'ChEA_2022.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'ENCODE and ChEA', 'path': enrichr_dir, 'file': 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'ARCHS4 TF Coexp', 'path': enrichr_dir, 'file': 'ARCHS4_TFs_Coexp.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'TRRUST 2019', 'path': enrichr_dir, 'file': 'TRRUST_Transcription_Factors_2019.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'TRANSFAC and JASPAR PWMs', 'path': enrichr_dir, 'file': 'TRANSFAC_and_JASPAR_PWMs.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
    ]

def get_pathway_databases(species):
    """パスウェイデータベースの定義を返す"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"
    
    return [
            {'name': 'MSigDB Hallmark', 'path': msigdb_dir, 'file': 'h.all.v2023.2.Hs.symbols.gmt' if species=='human' else 'mh.all.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB CP', 'path': msigdb_dir, 'file': 'c2.cp.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB GOBP', 'path': msigdb_dir, 'file': 'c5.go.bp.v2023.2.Hs.symbols.gmt' if species=='human' else 'm5.go.bp.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB Reactome', 'path': msigdb_dir, 'file': 'c2.cp.reactome.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.reactome.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr Reactome 2022', 'path': enrichr_dir, 'file': 'Reactome_2022.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr GOBP 2023', 'path': enrichr_dir, 'file': 'GO_Biological_Process_2023.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB KEGG legacy', 'path': msigdb_dir, 'file': 'c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.kegg_legacy.v2023.2.Hs.symbols2Mouse.txt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr KEGG 2021', 'path': enrichr_dir, 'file': 'KEGG_2021_Human.gmt' if species=='human' else 'KEGG_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
        ]

def get_celltype_databases(species):
    """パスウェイデータベースの定義を返す"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"
    
    return [
            {'name': 'MSigDB cell type signature', 'path': msigdb_dir, 'file': 'c8.all.v2023.2.Hs.symbols.gmt' if species=='human' else 'm8.all.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'CellMarker 2024', 'path': enrichr_dir, 'file': 'CellMarker_2024.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'CellMarker augmented 2021', 'path': enrichr_dir, 'file': 'CellMarker_Augmented_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'Tabula Sapiens', 'path': enrichr_dir, 'file': 'Tabula_Sapiens.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'Tabula Muris', 'path': enrichr_dir, 'file': 'Tabula_Muris.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'PanglaoDB augmented', 'path': enrichr_dir, 'file': 'PanglaoDB_Augmented_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
        ]

def run_summary_analysis(databases, gene_list, n_background):
    """
    共通の解析実行関数
    """
    all_results = {}
    progress_bar = st.progress(0)
    
    # 2列レイアウトを作成
    col1, col2 = st.columns(2)
    
    # 各データベースで解析を実行
    for idx, db in enumerate(databases):
        try:
            filepath = os.path.join(db['path'], db['file'])
            source = db['source']
            target = db['target']
            
            if db.get('type') == 'dorothea':
                net = pd.read_csv(filepath, sep='\t')
                net = net[net['confidence'].isin(db['filter']['confidence'])]

            elif db.get('type') == 'CollecTRI':
                net = pd.read_csv(filepath, sep='\t')

            else:
                net = dc.read_gmt(filepath)
            
            # キャッシュされた関数を使用して解析実行
            ora_res = run_single_database_analysis(
                filepath, 
                gene_list, 
                source, 
                target, 
                n_background,
                db_type=db.get('type'),
                filter_dict=db.get('filter')
            )
            
            if ora_res is not None:
                # 結果を保存
                all_results[db['name']] = ora_res
                progress_bar.progress((idx + 1) / len(databases))
                
                # プロットを表示（2列に分ける）
                with col1 if idx % 2 == 0 else col2:
                    if len(ora_res) > 0:
                        display_results(db, ora_res)
            else:
                st.warning(f"Database file not found: {filepath}")
                
        except Exception as e:
            st.error(f"Error processing {db['name']}: {str(e)}")

    progress_bar.empty()
    return all_results

def display_results(db, ora_res):
    """
    共通の結果表示関数
    """
    with st.expander(f"📊 {db['name']}", expanded=True):
        # プロット作成
        fig, ax = plt.subplots(figsize=(6, 4))
        
        # 上位5つを取得し、-log10(FDR)でソート
        plot_data = ora_res.head(5).copy()
        plot_data['neg_log10_fdr'] = -np.log10(plot_data['FDR p-value'])
        plot_data = plot_data.sort_values('neg_log10_fdr', ascending=True)
        
        # バープロット作成と色付け
        y_pos = range(len(plot_data))
        scores = plot_data['neg_log10_fdr']
        cmap = plt.cm.Reds
        colors = cmap(scores / scores.max())
        bars = ax.barh(y_pos, scores)
        
        for bar, color in zip(bars, colors):
            bar.set_color(color)
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['Term'].str[:50])
        ax.set_xlabel('-log10(FDR)')
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()
        
        # 詳細表示ボタン
        if st.button(f"Analyze {db['name']} in detail"):
            display_detailed_analysis(ora_res)

def display_detailed_analysis(ora_res):
    """
    共通の詳細解析表示関数
    """
    try:
        fig = create_enrichment_plot(ora_res, '')
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error creating enrichment plot: {str(e)}")

    try:
        fig_dotplot = create_dotplot(ora_res)
        st.pyplot(fig_dotplot)
    except Exception as e:
        st.error(f"Error creating dot plot: {str(e)}")

    # 結果テーブルを表示
    st.markdown("### Results Table")
    st.dataframe(
        ora_res[['Term', 'FDR p-value', 'Features']].style.format({
            'FDR p-value': '{:.2e}'
        })
    )

@st.cache_data
def run_single_database_analysis(filepath, gene_list, source, target, n_background, db_type=None, filter_dict=None):
    """単一のデータベースに対する解析をキャッシュする"""
    if not os.path.exists(filepath):
        return None
    
    try:
        if db_type == 'dorothea':
            net = pd.read_csv(filepath, sep='\t')
            if filter_dict:
                for key, values in filter_dict.items():
                    net = net[net[key].isin(values)]

        elif db_type == "collecTRI":
            net = pd.read_csv(filepath, sep='\t')
        else:
            net = dc.read_gmt(filepath)
            
        ora_res = dc.get_ora_df(gene_list, net, source=source, target=target, 
                               n_background=n_background, verbose=False)
        
        return ora_res.sort_values('FDR p-value', ascending=True)
    except Exception as e:
        st.error(f"Error in run_single_database_analysis: {str(e)}")
        return None

@st.cache_data
def create_enrichment_plot(ora_res_detail, title):
    """エンリッチメントプロットを生成してキャッシュする"""
    enr_pvals = ora_res_detail[['FDR p-value']]
    enr_pvals.index = ora_res_detail['Term']
    enr_pvals.values[enr_pvals.values == 0] = np.min(enr_pvals.values[enr_pvals.values != 0])
    enr_pvals = -np.log10(enr_pvals)
    
    g = dc.plot_barplot(enr_pvals.T, 'FDR p-value', vertical=True, top=15,
                       figsize=(8, 6), return_fig=True)
    g.gca().invert_yaxis()
    return g

@st.cache_data
def create_dotplot(ora_res_detail):
    """ドットプロットを生成してキャッシュする"""
    ora_res_detail = ora_res_detail.copy()
    ora_res_detail['count'] = ora_res_detail['Features'].str.split(';').str.len()
    
    max_count = ora_res_detail['count'].max()
    scale = 4 / max_count if max_count <= 10 else (2 / max_count if max_count <= 50 else 1 / max_count)
    
    dotplot = dc.plot_dotplot(
        ora_res_detail.sort_values('Combined score', ascending=False).head(15),
        x='Combined score',
        y='Term',
        s='count',
        c='FDR p-value',
        scale=scale,
        figsize=(8, 6),
        return_fig=True
    )
    return dotplot

@st.cache_data
def run_pathway_analysis(gene_list, net, source, target, n_background):
    """パスウェイ解析を実行してキャッシュする"""
    ora_res_detail = dc.get_ora_df(gene_list, net, source=source, target=target, 
                                  n_background=n_background, verbose=False)
    return ora_res_detail.sort_values('FDR p-value', ascending=True)


# プロットスタイルを設定する関数
def set_plot_style():
    plt.style.use('default')  # デフォルトスタイルをリセット
    plt.rcParams['figure.facecolor'] = 'white'  # 図の背景色を白に
    plt.rcParams['axes.facecolor'] = 'white'    # プロット領域の背景色を白に
#    sns.set_style("white")                      # seabornのスタイルを白背景に
    
# グラフ生成前に呼び出す
set_plot_style()


def delete_file(file_list):
    for i in file_list:
        if os.path.exists(i):
            os.remove(i)

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


#@st.cache_data
#def read_csv(file, index_col=None, sep=',', header = 0):
#    df_xl = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
#    return df_xl

@st.cache_data
def run_method(method, mat, net, source, target, weight, verbose=True, min_n = 0):
    norm = None
    corr = None
    if method == 'ulm':
        score, pvalue = dc.run_ulm(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'consensus':
        score, pvalue = dc.run_consensus(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'mlm':
        score, pvalue = dc.run_mlm(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'wsum_norm':
        score, norm, corr, pvalue = dc.run_wsum(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
        score = norm # wsum z-scoreを使う
    elif method == 'viper':
        score, pvalue = dc.run_viper(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    if method == 'wsum_norm':
        return score, pvalue, norm, corr
    else:
        return score, pvalue

@st.cache_data
def run_GSEA_df(mat, stat='stat', net = 'net', source='source', target='target', min_n=5, seed=42):
    GSEA_res = dc.get_gsea_df(mat, stat=stat, net = net, source=source, target=target, min_n=min_n, seed=seed)
    return GSEA_res

@st.cache_data
def calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch):
    orig_len = len(df)
    df = df[np.isfinite(df[P_column]) & pd.notnull(df[P_column])]     # FCやpがNAのものを除く
    df = df[np.isfinite(df[FC_column]) & pd.notnull(df[FC_column])]    # FCやpがNAのものを除く
    if len(df) < orig_len:
        st.warning("The P or FC columns contain inf or NA")
    inv_parameter = 1
    if inv_switch:
        inv_parameter = -1
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
        # 0e0がとして読まれる LogFCも0のはず
        # FC=0かつp=0の遺伝子を特定
        problematic_mask = (df[FC_column] == 0) & (df[P_column] == 0)
        if any(problematic_mask):
            st.warning(f"Found {sum(problematic_mask)} genes with FC=0 and p=0. These will be excluded from analysis.")
            excluded_genes = df.loc[problematic_mask, Gene_column].tolist()
            st.write("Excluded genes:", ", ".join(excluded_genes[:10]), "..." if len(excluded_genes) > 10 else "")
            df = df[~problematic_mask]
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
            #max_score = np.log10(1e-324) # 1e-324 == 0でTRUEになる pythonでも同じ 1e-324 + 1e-323でも計算される
            max_score = -324
            st.write("\nMax score: "+str(max_score))
        # 順位付けのためにFCの値を足す
        df.loc[(p_0 & (df.loc[:,FC_column]>0)),'score'] = max_score * -1 + df.loc[:,FC_column]  * inv_parameter#条件を括弧で囲むこと！！！
        df.loc[(p_0 & (df.loc[:,FC_column]<0)),'score'] = max_score + df.loc[:,FC_column] * inv_parameter
        st.write('Ranking score are -log10(P-values)')
    return df['score'].to_frame() #DFに変換してから返す



def set_back_func():
    input_file_type = st.radio(
        "Data format of the file containing all gene names (e.g., DESeq2 output):",
        ('tsv','csv', 'excel'))
    uploaded_file = st.file_uploader("Upload a file containing gene names (e.g., gene list, DESeq2, Homer)", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, header = None, index_col = None)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', header=None, index_col = None)
        else:
            df = read_excel(uploaded_file, index_col = None, header = None)

        # もし1列のデータで最初にGeneがないとき
        if df.shape[1] == 1:
            bk_genes = df.iloc[:,1].values
            if bk_genes[0] == "Gene" or bk_genes[0] == "GENE":
                bk_genes = bk_genes[1:]

        else:
            df.columns = df.iloc[0,:].tolist()  # transposeすると狂うので、transposeした後にcolumnsを決める
            df = df.drop(0, axis = 0) # 1行目を列名にして除く

            st.write(df.head())
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
        #       df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
                # annotation/divergence以前を除く
                df = df.loc[:,'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                Gene_column = "Gene"
                st.write("Converted Annotation/Divergence to gene symbols.")

            elif "Gene" in content:
                Gene_column =  "Gene"
            else:
                Gene_column =  st.selectbox(
                'Select gene column',
                content)

            bk_genes = df[Gene_column].values
        n_background = len(set(df[Gene_column].to_list()))
        st.write('Background gene number: ' + str(n_background))
        return n_background

#####
# --- Initialising SessionState ---
if "decouplercalc" not in st.session_state:
      st.session_state.decouplercalc = False
if "ORA" not in st.session_state:
      st.session_state.ORA = False

# temp内に保存する
# --- Initialising SessionState ---
if "dc_temp_dir" not in st.session_state:
    st.session_state.dc_temp_dir = True
    dc_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
    os.mkdir(dc_temp_dir)
    st.session_state.dc_temp_dir = dc_temp_dir
else:
    dc_temp_dir = st.session_state.dc_temp_dir

# progeny 読み込みの制御

if "net" not in st.session_state:
    st.session_state.net = None

if "num_progeny" not in st.session_state:
    st.session_state.num_progeny = "500"

if 'up_list' not in st.session_state:
    st.session_state.up_list = None
if 'down_list' not in st.session_state:
    st.session_state.down_list = None
if 'P_column' not in st.session_state:
    st.session_state.P_column = None

#############################


st.markdown('''## Signal pathway and TF activity inference using decoupleR
##### For rank mode, upload DEG result file or rank file
###### e.g.) FD-unshrunk DESeq2 file
###### decopuler orinigally uses stat values from DESeq2.
---
''')

use_upload = "No"

Analysis_mode = st.radio(
    "##### Analysis mode:",
    ('Rank','Over-representation'), key='Rank')
st.markdown("---")


if Analysis_mode == "Rank":
    input_file_type = st.radio(
    "Data format of DESeq2 result or rank file:",
    ('tsv','csv', 'excel', 'rank'), index = 0)
    st.write("DESeq2 res: tsv, csv, excel; rank file: rank")

    rnk_input = False

    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx','rnk','rank'])
    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, index_col = 0)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', index_col = 0)
        elif input_file_type == 'rank':
            rnk_input = True
            df = read_csv(uploaded_file, sep = '\t', index_col = 0, header = None)
        else:
            df = read_xl(uploaded_file, index_col = 0)

        down_file_name = os.path.basename(uploaded_file.name)

        df.iloc[0:3,:]

else: #ORA
    ORA_mode = st.radio(
    "##### How to provide gene list:",
    ('DEG result file','Cluster info file', "Input genes"), key='DEG result file')
    n_background = None
    set_back = False

    if ORA_mode == 'Cluster info file':
        cluster_file_type = st.radio(
            "Data from:",
            ('auto', 'tsv', 'csv', 'excel'))
        cluster_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
        if cluster_file is not None:
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


            if "cluster_submitted" not in st.session_state:
            	st.session_state.cluster_submitted = False


            with st.form("Set column"):
                # Allow user to select column names
                gene_column = st.selectbox("Select the column containing gene names:", df_cluster.columns, index=df_cluster.columns.get_loc("Gene") if "Gene" in df_cluster.columns else 0)
                cluster_column = st.selectbox("Select the column containing cluster information:", df_cluster.columns, index=df_cluster.columns.get_loc("Cluster") if "Cluster" in df_cluster.columns else 0)
                
                # Get unique clusters
                clusters = sorted(df_cluster[cluster_column].unique())           
                column_submitted = st.form_submit_button("Set gene and cluster columns")

            with st.form("Choose cluster"):
                selected_cluster = st.multiselect("Cluster:", clusters)
                cluster_submitted = st.form_submit_button("Select clusters")
                st.session_state.cluster_submitted = True

                
            if st.session_state.cluster_submitted:
#                cluster_genes = df_cluster[df_cluster[cluster_column] == selected_cluster][gene_column].tolist()
                cluster_genes = df_cluster[df_cluster[cluster_column].isin(selected_cluster)][gene_column].tolist()
                gene_list = list(set(cluster_genes))
                st.write(f"{len(gene_list)} genes in cluster {selected_cluster}:")
                st.write('Gene list: ' + ' '.join(gene_list[:10]))
                st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")

                set_back = st.checkbox('Set background genes?', value=False)
                if set_back:
                    n_background = set_back_func()

    elif ORA_mode == 'Input genes':
        st.markdown("##### Or input genes (comma, space, CR separated):")
        genes = st.text_input("genes",label_visibility = 'collapsed')
        gene_list = []
        if len(genes) > 0:
#            if ',' not in genes:
#                gene_list = genes.split(' ')
#            else:
#                genes =  ''.join(genes.split()) #空白を除く
#                gene_list = genes.split(',')
            raw_genes = re.split(r'[,\s]+', genes)
            # Remove spaces from each gene name and filter out empty strings
            gene_list = [re.sub(r'\s', '', gene) for gene in raw_genes if gene.strip()]

            genes = list(filter(lambda x: x != "", genes)) #空白を除く
            gene_list =sorted(set(gene_list), key=gene_list.index)
            st.write(gene_list[:3])
            st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")

            set_back = st.checkbox('Set background genes?', value=False)
            if set_back:
                n_background = set_back_func()


    else: # use DEG file

        if "up" not in st.session_state:
            st.session_state.up = None
        if "down" not in st.session_state:
            st.session_state.down = None

        use_upload = 'Yes' # deseq2がないときはYes
        df_res = None
        if 'deseq2' in st.session_state:
            st.write("There is a deseq2 result in the cache. If you use it, do not upload a new file.")
            if st.session_state.deseq2 is not None:
                use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
            if use_upload == "No" and "df_res" not in st.session_state: #df_resを作る
                if 'df_res' not in st.session_state:
                    df_res = st.session_state.deseq2
                    df_res['Gene'] = df_res.index
                    if  "deseq2_uploaded_file_name" in st.session_state:
                        file_name_head = st.session_state.deseq2_uploaded_file_name
                    else:
                        file_name_head = "res"
                    input_file_type = 'tsv'
                    if "Row_name" in df_res.columns.to_list(): # Row_nameを含むとき
                        df_res = df_res.set_index('Row_name')
                        df_res.index.name = "Gene"
                    st.session_state.df_res = df_res # df_resに記録
            elif "df_res" in st.session_state:
                df_res = st.session_state.df_res
            else:
                st.write("something is wrong...")



        if use_upload == 'Yes': # DEGの結果を使うとき df_res
    #        st.session_state.deseq2 = None
            input_file_type = st.radio(
                "Data from:",
                ('auto', 'tsv', 'csv', 'excel'))
            uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
            if uploaded_file is not None:
                if input_file_type == 'auto':
                    try:
                        df_res = read_csv(uploaded_file, sep = None)
                    except:# excel
                        df_res = read_excel(uploaded_file)


                if input_file_type == "csv":
                    df_res = read_csv(uploaded_file)
                elif input_file_type == 'tsv':
                    df_res = read_csv(uploaded_file, sep = '\t')
                elif input_file_type == 'excel':
                    df_res = read_xl(uploaded_file)
            #    st.session_state.deseq2 = df_res
                st.session_state.df_res = df_res
                file_name_head = os.path.splitext(uploaded_file.name)[0]
                st.session_state.deseq2_uploaded_file_name = file_name_head
            #    if 'seurat_res' not in st.session_state: #sueratのしょりをしたらTrue
            #        st.session_state.seurat_res = False

            else:
                sys.exit(1)
           ##### file upload

        if df_res is not None:
            if 'gene_list' not in st.session_state:
                st.session_state.gene_list = None

            st.write(df_res.head(3))


#            if use_upload == 'Yes': #Seuratをクリックしたとき
    #        if 'seurat_res' not in st.session_state: #sueratのしょりをしたらTrue
    #            st.session_state.seurat_res = False
    #        elif not st.session_state.seurat_res: #まだseuratの処理をしていないとき
            seurat = st.checkbox('Seurat resuls?', value=False)
            if seurat:
#                df_res.columns = ['col', 'p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene']
#                df_res = df_res.drop("col", axis = 1)

            #    st.write(df_res.head())

                # ユニークなクラスターを特定
                clusters = df_res['cluster'].unique()
                #clusters = list(set(df_res['cluster']))
            #    st.write(clusters)
            #    st.write(f"Unique clusters: {clusters}")

                # 各クラスターのDataFrameを保持するリストを作成
                cluster_dfs = []

                for cluster in clusters:
                    # このクラスターのデータをフィルタリング
                    cluster_data = df_res[df_res['cluster'] == cluster].copy()

            #        st.write(f"Processing cluster {cluster}, shape: {cluster_data.shape}")
                    # 'gene' 列をインデックスに設定
                    cluster_data = cluster_data.set_index('gene')
                    # 'cluster' カラムを削除
                    cluster_data = cluster_data.drop('cluster', axis=1)
                    # カラム名を変更してクラスター番号を先頭に
                    new_columns = [f'{cluster}_{col}' for col in cluster_data.columns]
                    cluster_data.columns = new_columns
                    # クラスターDataFrameのリストに追加
                    cluster_dfs.append(cluster_data)

                # 全てのクラスターDataFrameをマージ
                result_df = pd.concat(cluster_dfs, axis=1)

#                # カラムを並べ替えて 'gene' を最初の列にする
#                cols = ['gene'] + [col for col in result_df.columns if col != 'gene']
#                result_df = result_df[cols]
                df_res = result_df
                df_res["Gene"] = df_res.index.to_list()
                st.write(df_res.head(3))
                #st.session_state.seurat_res = True こうすると最初から戻ってきたときに動かない

            content = df_res.columns.tolist()
            p_patterns = ['p.val', 'pvalue', 'p-val', 'p val', 'p_val', 'pval']
            pvalue = [i for i in content if any(p in i.lower() for p in p_patterns)] #and 'adj.pval' not in i.lower()]
            fc_patterns = ['log2fc', 'fold change', 'log2foldchange', 'coef', 'logfc']
            fc = [i for i in content if any(pattern in i.lower() for pattern in fc_patterns)]
            gene = [i for i in content if (i not in pvalue) and (i not in fc)]
            P_column = st.selectbox(
                'Select adjusted P-value column',
                pvalue)
            # ジャロ・ウィンクラー距離法
            #JW_dist = [Levenshtein.jaro_winkler(P_column, x) for x in fc]

            JW_dist = [jaro_winkler_similarity(P_column, x) for x in fc]
            try:
                FC_column = st.selectbox(
                    'Select FC column',
                    fc, index = JW_dist.index(max(JW_dist)))
            except:
                FC_column = st.selectbox(
                    'Select FC column', fc)

            file_name_add = P_column
            file_name_add = file_name_add.replace('.adj.pvalue','')
            file_name_add = file_name_add.replace(' adj. p-value','')
            file_name = file_name_add

            set_gene = st.checkbox("Specify gene column?", value=False)

            if set_gene:
                Gene_column =  st.selectbox(
                'Select gene column',
                gene)
            else:

                if "Annotation/Divergence" in content:
                    pattern = "([^|]*)"
                    Gene_column = 'Annotation/Divergence'
                    df_res.loc[:,'Annotation/Divergence'] = df_res.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
                    st.write("Converted Annotation/Divergence to gene symbols.")
                elif "Gene" in content:
                    Gene_column =  "Gene"
                else:
                    Gene_column =  st.selectbox(
                    'Select gene column',
                    gene)

            if use_upload == 'Yes':
                # indexをGeneにする
                df_res.index = df_res[Gene_column].tolist()

            if 'df_res' in st.session_state:
                st.session_state.df_res = df_res

            with st.form("Basic settings:"):
                p_or_top = st.radio("P threshold or top/bottom:", ('P_threshold','top', 'both'))
                up_or_down = st.radio("Up or down genes:", ('Up','Down'))
                sort_val = st.radio("Top based on:", ('FC','P value'), index = 0)
                p_thre = st.number_input('Threshold for adjusted P', min_value =0.000, max_value=1.000,
                step =0.002, value=0.050)
                FC_thre = st.number_input('Threshold for log FC', min_value =0.0, step =0.1, value=0.0)

                top_n = st.number_input('Number of top genes', min_value =1, step = 1, value=50)

                st.write("Top genes include only those that meet the adjusted P-value threshold.")
                set_back = st.checkbox('Set background genes?', value=False)
                st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")
                ssubmitted_basic = st.form_submit_button("Change the parameters")
                up=None
                down=None
                df_thre  = None

                if p_or_top == 'top': #topのみのときはp_threは1に強制的に変更する
                    p_thre = 1


                # "Basic settings:"フォームの後、df_threを作成する前に追加
                # FC=0かつp=0の遺伝子を確認・除外
                problematic_mask = (df_res[FC_column] == 0) & (df_res[P_column] == 0)
                if any(problematic_mask):
                    st.warning(f"Found {sum(problematic_mask)} genes with FC=0 and p=0 in the DEG results.")
                    excluded_genes = df_res.loc[problematic_mask, Gene_column].tolist()
                    st.write("These genes will be excluded from ORA analysis:", ', '.join(excluded_genes[:10]), 
                             "..." if len(excluded_genes) > 10 else "")
                    
                    # df_resから除外
                    df_res = df_res[~problematic_mask]

                df_thre = df_res.copy(deep=True)
                df_thre = df_thre[df_thre[P_column] < p_thre]
                df_thre = df_thre[abs(df_thre[FC_column]) > FC_thre]
                if sort_val == "P value":
                    df_thre =df_thre.sort_values(P_column, ascending=True)
                    up = df_thre[df_thre[FC_column]>0].index.to_list()
                    down = df_thre[df_thre[FC_column]<0].index.to_list()
                else:
                    df_thre = df_thre.sort_values(FC_column, ascending=False)
                    up = df_thre[df_thre[FC_column]>0].index.to_list()
                    df_thre = df_thre.sort_values(FC_column, ascending=True)
                    down = df_thre[df_thre[FC_column]<0].index.to_list()
                if p_or_top in ["top",  "both"]:
                    up = up[:top_n]
                    down = down[:top_n]
                up = sorted(set(up), key=up.index) #重複を除く
                down = sorted(set(down), key=down.index)

                st.session_state.up = up
                st.session_state.down = down
            try:
                if len(df_thre) > 0:
                    if up_or_down == "Up":
                        st.write(','.join(st.session_state.up))
                        st.write("Number of genes: " + str(len(up)))
                        gene_list = st.session_state.up
                    else:
                        st.write(','.join(st.session_state.down))
                        st.write("Number of genes: " + str(len(down)))
                        gene_list = st.session_state.down
                    st.session_state.gene_list = gene_list
            except:
                pass



        if set_back:
            #bk_genes = df_res[Gene_column].values
            n_background = len(set(df_res[Gene_column].to_list()))
            st.write('Background gene number: ' + str(n_background))

        gene_list = st.session_state.gene_list

        if p_or_top == 'top':
            gene_list_file = P_column + "." + up_or_down + "-" + str(top_n) + '.txt'
        else:
            gene_list_file = P_column + "." + up_or_down + ".p" + str(p_thre) + '.txt'

        tsv_string = '\n'.join(gene_list)


        st.download_button(
            label="Download the gene list",
            data=tsv_string,
            file_name=gene_list_file,
            mime="text/tsv"
        )



if 'df' in locals()  or 'gene_list' in locals() or 'df_res' in locals(): # dfかgenesを入力したとき


    if Analysis_mode == "Rank":
        generated_rnk = False
        use_stat = False
        rank_calc = False
        rank_metric = None
        if not rnk_input:
            rank_metric = st.radio(
                "Ranking metric:",
                ('sign(LFC) x -log10(P)', 'LFC x -log10(p)', "DESeq2 stat"), index = 0)
                # calculate stat value
            # indexをGeneカラムにコピー
            df['Gene'] = df.index
            # indexの名前を除く
            df.index.name = None
            content = df.columns.tolist()
        if rank_metric ==  "DESeq2 stat" and not rnk_input:
            statvalue = [i for i in content if ('stat' in i)]
            stat_column = st.selectbox('Select stat column', statvalue)

        elif not rnk_input:
            generated_rnk = True # decoupler内でrnk fileを作ったとき
            st.write("Select pvalue and logFC")
            p_patterns = ['p.val', 'pvalue', 'p-val', 'p val', 'p_val', 'pval']
            pvalue = [i for i in content if any(p in i.lower() for p in p_patterns) and 'adj.pval' not in i.lower()]
 #           pvalue = [i for i in content if (('pval' in i) or ('p-val' in i) or ('p val' in i) or ('p_val' in i) or ('pval' in i) or ('p.val' in i)) and ('adj.pval' not in i)]
            fc = [i for i in content if ('log2FC' in i) or ('Fold Change' in i) or ('log2FoldChange' in i) or ('coef' in i) or ('logFC' in i)]
            gene = [i for i in content if (i not in pvalue) and (i not in fc)]
            P_column = st.selectbox('Select P-value column', pvalue)
            stat_column = re.match(r'([^\.]+)', P_column).group(1) #名前を変更する
            # ジャロ・ウィンクラー距離法
            JW_dist = [jaro_winkler_similarity(P_column, x) for x in fc]
            try:
                FC_column = st.selectbox(
                    'Select FC column',
                    fc, index = JW_dist.index(max(JW_dist)))
            except:
                FC_column = st.selectbox(
                    'Select FC column', fc)

            if "Gene" in content:
                Gene_column =  "Gene"
            elif "Symbol" in content:
                Gene_column =  "Symbol"
            else:
                Gene_column =  st.selectbox(
                'Select gene column',
                gene)

            inv_switch = st.checkbox('Invert the sign')

            df_sub = df[[P_column, FC_column]]
            df_score = calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch)
            stat_column = 'score'
            mat = df_score.transpose()
            rank_calc = True
            st.write(mat.iloc[:,:10])
            df = df_score

            # P_columnが変更されたらファイルを消去する
            if st.session_state.P_column != P_column:
                shutil.rmtree(dc_temp_dir)
                os.mkdir(dc_temp_dir)
                st.session_state.P_column = P_column

        else: # rank file
            stat_column = 'Rank'
            df.columns = ['Rank']

        if not rank_calc:

            if list(df.index.duplicated()).count(True) > 0:
                st.markdown("#### There are duplicated genes.")
                st.write('Dupliated genes:' +  ', '.join(list(df[df.index.duplicated()].index)))
                st.write("The first instances will be kept.")
                st.markdown("---")
                df = df[~df.index.duplicated(keep='first')]


            try:
                mat = df[[stat_column]].T
                st.write(mat.iloc[:,:10])
            except:
                st.markdown("#### Error. Prerank file?")
                sys.exit(1)


#--------------------------------------
    st.markdown("---")
    if Analysis_mode == 'Rank':
        species = st.radio("Species:", ('mouse','human'), index = check_species_index(df.index.to_list()[:50]))
        path = st.radio(
        "Pathway:",
        ('PROGENy', 'CollecTRI', 'DoRothEA', 'mSigDB', 'Enrichr', 'Homemade', 'Your own GMT file'))
        st.write("""
        PROGENy: signaling pathway\n
        CollecTRI: TF targets\n
        DoRothEA: TF targets (subset of TRI)
        """)
    else:
        species = st.radio("Species:", ('mouse','human'), index = check_species_index(gene_list[:50]))

        path = st.radio(
            "Pathway:",
            ( 'Pathway_summary','TF_summary', 'Celltype_summary', 'CollecTRI', 'DoRothEA', 'mSigDB', 'Enrichr', 'Homemade', 'Your own GMT file'),
            index = 0)
        st.write("""
        PROGENy: signaling pathway\n
        CollecTRI: TF targets\n
        DoRothEA: TF targets (subset of TRI)
        """)

    if path == 'PROGENy':
        num_progeny = st.radio("Number of top genes:", ('500','2000','5000', 'all'), index = 0,
            help="各パスウェイの活性計算に使用するフットプリント遺伝子数。500（高信頼性）〜all（全遺伝子）。デフォルトの500推奨。")
        if st.button('Load PROGENy db') or (num_progeny != st.session_state.num_progeny):
            net = pd.read_csv('./db/progeny.' + species + "." + num_progeny + '.tsv', sep = '\t')
            source='source'
            target='target'
            weight = 'weight'
            st.session_state.net = net
            st.session_state.num_progeny = num_progeny

        elif st.session_state.net is None:
            st.stop()
        else:
            net = st.session_state.net
            source='source'
            target='target'
            weight = 'weight'

    elif path == 'CollecTRI':
        net = pd.read_csv('./db/TRI.' + species + '.tsv', sep = '\t')
        source='source'
        target='target'
        weight = 'weight'
    elif path == 'DoRothEA':
        conf = st.selectbox("Minimal confidence level:", ("A","B","C","D"), index = 1)
        st.write("B is recommended.")
        conf_levels = ['A']
        if conf == 'B':
            conf_levels = ['A','B']
        elif conf == 'C':
            conf_levels = ['A', 'B','C']
        elif conf == 'D':
            conf_levels = ['A', 'B','C', 'D']
        net = pd.read_csv('./db/dorothea.' + species + '.tsv', sep = '\t')
        source='tf'
        target='target'
        weight = 'mor'
        net = net[net['confidence'].isin( conf_levels)]


    elif path in ['TF_summary', 'Pathway_summary', 'Celltype_summary']:
        source = 'source'
        target = 'target'
        weight = None
        
        if path == 'TF_summary':
            st.write("### Transcription Factor Activity Summary")
            databases = get_tf_databases(species)  # TFデータベース定義
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()
            
        elif path == 'Pathway_summary':
            st.write("### Pathway Enrichment Summary")
            databases = get_pathway_databases(species)  # パスウェイデータベース定義
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()

        elif path == 'Celltype_summary':
            st.write("### Cell Type Summary")
            databases = get_celltype_databases(species) 
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()



    elif path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
        source='source'
        target='target'
        weight = None
        if path == 'mSigDB':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/mSigDB_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/mSigDB"
        elif path == 'Enrichr':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/enrichr_gmt_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/enrichr_gmt"
        elif path == 'Homemade':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt"

        files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
        files_file.sort()
        key_index = None
        net_list = []
        if path == 'mSigDB':
            key_index = [item for item in files_file if "h.all" in item]
        else:
            key_index = files_file[0]
        GO_name = st.multiselect('Select gene sets',files_file, default = key_index)
        if len(GO_name) > 0:
            for i in GO_name:
                net_list.append(dc.read_gmt(dir_path + "/" + i))
            net = pd.concat(net_list, ignore_index = True)

    else:
        uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt','gmt'])
        if uploaded_gmt is not None:
            source='source'
            target='target'
            weight = None
            GO_name = uploaded_gmt.name
            path = uploaded_gmt.name.replace('.gmt','')
            stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
            s = stringio.read()
            #スペースがあるとエラーになる
            s = s.replace(' ', '_')

            with open('temp.gmt', mode='w') as f:
                f.write(s)
            net = dc.read_gmt('temp.gmt')
    #        os.remove("temp.gmt")

        else:
            st.stop()

    if "net" in locals():
        st.write('Pathway data:')
        st.write(net.head(2))
    #-------------------------------------
        if Analysis_mode == 'Rank':
            st.markdown("##### For CollecTRI TF enrichment: ULM, PROGENy pathway enrichment: MLM")
            st.write("For other weightless datasets (e.g., gmt gene sets), ulm treats each gene set independently while mlm accounts for all of them at the same time during fitting. If your gene sets contain many common targets, it is better to use ulm, as mlm may not be able to run at all.")
            method_index = 2
            if path ==  'CollecTRI' or path == 'DoRothEA':
                method_index = 0
            method = st.radio("Analytical model:", ("ulm", "consensus", "mlm",  "wsum_norm", 'viper', 'GSEA'), index = method_index)
            st.write('https://decoupler-py.readthedocs.io/en/latest/api.html#methods')
            st.write("Consensus score is calculated from ulm, mlm and wsum_norm.")


            min_n = 5
            min_n = int(st.number_input('Minimum number of targets per TF/set', min_value =0, step =1, value=5))
            st.write("Gene sets that have fewer genes in the uploaded data than this number will be ignored.")



            if method != 'GSEA':
                if 'dc_res' not in st.session_state:
                    st.session_state.dc_res = None
                    st.session_state.dc_score = None

                if st.button('Run analysis'): # or not st.session_state.decouplercalc:
                    if method == 'wsum_norm':
                        score, pvalue, norm, corr  = run_method(method= 'wsum_norm', mat=mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
                        score = norm # wsum z-scoreを使う
                    else:
                        score, pvalue = run_method(method= method, mat=mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)

                    res = score.T
                    res['p-value'] = pvalue.iloc[0,:]
                    res['adj.p_value'] = dc.p_adjust_fdr(pvalue.iloc[0,:])
                    res.columns = ['Score','p-value','adj.p_value']

                    res = res.sort_values('Score', ascending=False)
                    st.dataframe(res.head())
                    st.session_state.dc_res = res
                    st.session_state.dc_score = score

                st.markdown("##### You must run the analysis again if the parameters are changed.")
                st.markdown("---")

                if st.session_state.dc_res is not None:
                    res = st.session_state.dc_res
                    score = st.session_state.dc_score

                    with st.sidebar:
                        st.markdown("#### Barplot parameters")
                        if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                            if path == 'Homemade':
                                GO_file_name = '.'.join(GO_name)
                            else:
                                GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                            bar_name_org =  GO_file_name + '_barplot'
                        elif path == 'DoRothEA':
                            bar_name_org =  path + '-' + conf + '_barplot'
                        else:
                            bar_name_org =  path + '_barplot'
                        bar_name_head = st.text_input("Barplot file name: ", value = bar_name_org)
                        bar_name = bar_name_head + ".pdf"
                        bar_top = int(st.number_input('How many top terms to show', min_value =1, step =1, value=15))
                        bar_vertical = st.checkbox("Vertical plot?", value=True)
                        bar_vc = st.checkbox("Change center value?", value=False)
                        if bar_vc:
                            bar_v_center = st.number_input('Center value', value=0.0)
                        else:
                            bar_v_center = None
                        bar_x_size = st.number_input('X size', min_value =1, value=8)
                        bar_y_size = st.number_input('Y size', min_value =1, value=6)

                    bar_plot_name = dc_temp_dir + "/" + bar_name
                    fig_bar = dc.plot_barplot(score, stat_column, top=bar_top, vertical=bar_vertical, figsize = (bar_x_size, bar_y_size),
                        vcenter = bar_v_center,  return_fig=True) # save = dc_temp_dir + "/" + bar_name,
                    fig_bar.gca().invert_yaxis()
                    fig_bar.savefig(dc_temp_dir + "/" + bar_name,bbox_inches='tight')
                    st.pyplot(fig_bar)
                #        dc.plot_barplot(score, stat_column, top=bar_top, vertical=bar_vertical, figsize = (bar_x_size, bar_y_size),
                #            vcenter = bar_v_center, save = dc_temp_dir + "/" + bar_name)


                    # DECOUPLERでプロットを作成
                    fig = dc.plot_barplot(score, stat_column, top=bar_top, vertical=bar_vertical, 
                                         figsize=(bar_x_size, bar_y_size), vcenter=bar_v_center, return_fig=True)
                    
                    # 既存のカラーバーを削除
                    # 全てのaxesを取得
                    axes = fig.axes
                    # 最後のaxes（カラーバー）を削除
                    if len(axes) > 1:  # メインプロット以外にaxesがある場合
                        fig.delaxes(axes[-1])
                    
                    # メインのaxesを取得
                    ax = axes[0]
                    
                    # プロットデータの準備
                    plot_data = res.head(bar_top)
                    log_padj = -np.log10(plot_data['adj.p_value'])
                    
                    # カラーマップの作成
                    cmap = plt.cm.Reds
                    
                    # バーの色を変更
                    for i, bar in enumerate(ax.containers[0]):
                        bar.set_color(cmap(log_padj.iloc[i] / log_padj.max()))
                    
                    # 新しいカラーバーの追加
                    sm = plt.cm.ScalarMappable(cmap=cmap)
                    sm.set_array(log_padj)
                    cbar = fig.colorbar(sm, ax=ax)
                    cbar.set_label('-log10(adj P-value)')
                    
                    ax.invert_yaxis()
                    
                    # プロットの保存と表示
                    plt.tight_layout()
                    fig.savefig(dc_temp_dir + "/FDR_" + bar_name, bbox_inches='tight')
                    st.pyplot(fig)

                    if path == 'PROGENy' or path == 'CollecTRI':
                        tf_list = score.columns.to_list()
                        tf = st.selectbox('Select pathway to visualize', tf_list)
     #                   st.set_option('deprecation.showPyplotGlobalUse', False)

                        with st.sidebar:
                            st.markdown("#### Pathway plot parameters")
                            tf_name_head = st.text_input("TF file name: ", value = tf + '_targets')
                            tf_name = tf_name_head + ".pdf"
                            tf_top = int(st.number_input('How many top TFs to show', min_value =1, step =1, value=15))
                            tf_x_size = st.number_input('TF X size', min_value =1, value=8)
                            tf_y_size = st.number_input('TF Y size', min_value =1, value=6)


                        fig_targets = dc.plot_targets(df, stat=stat_column, source_name=tf, net=net, top=tf_top, figsize = (tf_x_size, tf_y_size),
                        save = dc_temp_dir + "/" + tf_name, return_fig=True)
                        st.pyplot(fig_targets)
                        if path == 'CollecTRI':
                            st.write("Weight shows positive/negative targets of the TF")

                    if generated_rnk:
                        add_head = P_column + '.' # rankファイルの名前をつける
                    else:
                        add_head = ""
                    if path == 'mSigDB' or path == 'Enrichr' or path == "Homemade":
                        if path == 'Homemade':
                            GO_file_name = ''.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method + "." + GO_file_name + '.tsv'
                        zip_name = down_file_name+ "." + add_head +  method + "." + GO_file_name
                    else:
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method+ "." + path + '.tsv'
                        zip_name = down_file_name + "." +add_head +  method + "." + path

                    if st.button('Prepare enrichment files to download'):
                        score.to_csv(out_file_name, sep = '\t')
                        shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)

                        with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                            btn = st.download_button(
                                label="Download Results of enrichment",
                            data=fp,
                            file_name=zip_name + ".zip",
                            mime = "zip"
                            )

            else: # GSEA
                if 'dc_gsea' not in st.session_state:
                    st.session_state.dc_gsea = None

                if st.button('Run GSEA analysis'):
                    GSEA_res = run_GSEA_df(mat.T, stat_column, net = net, source=source, target=target, min_n=min_n, seed=42)
                    st.session_state.dc_gsea = GSEA_res

                st.markdown("##### You must run the analysis again if the parameters are changed.")
                st.markdown("---")

                if st.session_state.dc_gsea is not None:
                    GSEA_res = st.session_state.dc_gsea

                    GO_terms = GSEA_res["Term"].to_list()
                    GSEA_res = GSEA_res.sort_values("FDR p-value", ascending = True)
                    nes_show = st.selectbox('Select term to visualize', GO_terms)

                    with st.sidebar:
                        st.markdown("#### GSEA enrichment plot parameters")
                        gsea_title_size = st.number_input('GSEA: title font size', min_value =1, value=14)
                        gsea_legend_size = st.number_input('GSEA: legend font size', min_value =1, value=12)
                        gsea_legend_x = st.number_input('GSEA: legend X pos', min_value =0.00, max_value = 1.00, value=0.2)
                        gsea_legend_y = st.number_input('GSEA: legend y pos', min_value =0.00, max_value = 1.00, value=0.5)
                        st.write('The left bottom is (0, 0). (0-1)')
                        gsea_x_size = st.number_input('GSEA: X size', min_value =1, value=5)
                        gsea_y_size = st.number_input('GSEA: Y size', min_value =1, value=5)

                    fig, d =  dc.plot_running_score(mat.T, stat_column, net = net, source=source, target=target,
                        set_name=nes_show, cmap='RdBu_r', figsize=(gsea_x_size, gsea_y_size), dpi=100, return_fig=True, save=None)
                    # fig, axのtuple objectがもどる dは遺伝子名

                    nes = GSEA_res[GSEA_res['Term']==nes_show]['NES'].iloc[-1]
                    fdr = GSEA_res[GSEA_res['Term']==nes_show]['FDR p-value'].iloc[-1]

                    s = "NES: " + str(nes) + "\nFDR:" + str(fdr)
                    #plt.text(len(mat.T)/10, 18, s, fontsize=14)
                    plt.figtext(gsea_legend_x, gsea_legend_y, s, fontsize = gsea_legend_size) # plotの中のcoordinateで書き込まれる。左下から。0-1
                    # axisを取り出して、タイトルを修正する
                    gsea_title = nes_show.replace("_", " ")
                    fig.axes[0].set_title(gsea_title, wrap=True, fontsize= gsea_title_size)
                    fig.axes[0].set_ylabel("Enrichment Score")

                    st.pyplot(fig)
                    fig.savefig(dc_temp_dir + "/"+ nes_show + '.pdf', format='pdf')

                    st.dataframe(GSEA_res)
                    if generated_rnk:
                        add_head = P_column + '.' # rankファイルの名前をつける
                    else:
                        add_head = ""
                    if path == 'mSigDB' or path == 'Enrichr' or path == "Homemade":
                        if path == 'Homemade':
                            GO_file_name = ''.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method + "." + GO_file_name + '.tsv'
                        zip_name = down_file_name+ "." + add_head +  method + "." + GO_file_name
                    else:
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method+ "." + path + '.tsv'
                        zip_name = down_file_name + "." +add_head +  method + "." + path

                    if st.button('Prepare GSEA files to download'):
                        GSEA_res.to_csv(out_file_name, sep = '\t')
                        shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)
                        with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                            btn = st.download_button(
                                label="Download GSEA Results",
                            data=fp,
                            file_name=zip_name + ".zip",
                            mime = "zip"
                            )


        else: # ORA

            if "ORA" not in st.session_state:
                st.session_state.ORA = False
            if "ORA_res" not in st.session_state:
                st.session_state.ORA_res = None

            if st.button('Run ORA analysis') or not st.session_state.ORA: # ボタンが押されない場合は再計算しない。
                try:
                    ORA_res = dc.get_ora_df(gene_list, net, source=source, target=target, n_background=n_background, verbose=False)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.markdown("#### Error. Possible wrong choice of gene sets or species.")
                    sys.exit(1)
                ORA_res = ORA_res.sort_values('FDR p-value', ascending = True)
                st.session_state.ORA = True
                st.session_state.ORA_res = ORA_res

            st.markdown("##### You must run the analysis again if the parameters are changed.")
            st.markdown("---")

            if st.session_state.ORA_res is not None and st.session_state.ORA:
                ORA_res = st.session_state.ORA_res

                # Set 0s to min p-value
                enr_pvals = ORA_res[['FDR p-value']]
                enr_pvals.index = ORA_res['Term']
                enr_pvals.values[enr_pvals.values == 0] = np.min(enr_pvals.values[enr_pvals.values != 0])
                enr_pvals = -np.log10(enr_pvals)

                with st.sidebar:
                    st.markdown("#### ORA barplot parameters")
                    if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                        if path == 'Homemade':
                            GO_file_name = '-'.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        bar_name_org =  GO_file_name + '_barplot'
                    elif path == 'DoRothEA':
                        bar_name_org =  path + '-' + conf + '_barplot'
                    else:
                        bar_name_org =  path + '_barplot'

                    if ORA_mode == 'Cluster info file' and st.session_state.cluster_submitted:
                        bar_name_org = bar_name_org + "." + "_".join(map(str, selected_cluster))
                    bar_name_head = st.text_input("ORA: Barplot file name: ", value = bar_name_org)
                    bar_name = bar_name_head + ".pdf"
                    bar_top = int(st.number_input('ORA: How many top terms to show', min_value =1, step =1, value=15))
                    bar_vertical = st.checkbox("ORA: Vertical plot?", value=True)
                    bar_vc = st.checkbox("ORA: Change center value?", value=False)
                    if bar_vc:
                        bar_v_center = st.number_input('ORA: Center value', value=0.0)
                    else:
                        bar_v_center = None
                    bar_x_size = st.number_input('ORA: X size', min_value =1, value=8)
                    bar_y_size = st.number_input('ORA: Y size', min_value =1, value=6)

                try:
                    fig_bar2 = dc.plot_barplot(enr_pvals.T, 'FDR p-value', vertical=bar_vertical, top=bar_top,
                        figsize = (bar_x_size, bar_y_size), vcenter = bar_v_center, return_fig=True)
                    fig_bar2.gca().invert_yaxis()
                    st.pyplot(fig_bar2)
                    st.markdown("#### Activity = -log10(adjP)")
                    st.markdown("###### -log10(0.05) = 1.301") 
                except Exception as e:
                    st.error(f"Error: {str(e)}")  
                    st.write("Probably little difference in FDR.")
                    st.markdown("#### The follwoing graph likely has no use!")
                    st.markdown("##### -log10(0.05) = 1.301")
                    st.write(enr_pvals)
                    vmn = enr_pvals['FDR p-value'].min()
                    vmx = enr_pvals['FDR p-value'].max()
                    vc = enr_pvals['FDR p-value'].mean()
                    fig_bar3 = dc.plot_barplot(enr_pvals.T, 'FDR p-value', vertical=bar_vertical, top=bar_top,
                        figsize = (bar_x_size, bar_y_size), vmin = vmn, vmax=vmx, vcenter =vc, return_fig=True)
                    fig_bar3.gca().invert_yaxis()
                    st.pyplot(fig_bar3)
                    st.markdown("##### Activity = -log10(adjP)")




                # Log-transform
                enr_pvals = -np.log10(enr_pvals)
                # countを追加
                ORA_res['count'] = ORA_res['Features'].str.split(';').str.len()
                st.dataframe(ORA_res)

                # 最大カウント値に基づいてscaleを計算
                max_count = ORA_res['count'].max()

                # 経験則として、小さいカウントの場合は大きめのscale、大きいカウントの場合は小さめのscaleを使用
                if max_count <= 10:
                    scale = 4 / max_count  # 小さいカウントの場合、より大きなscale
                elif max_count <= 50:
                    scale = 2 / max_count
                else:
                    scale = 1 / max_count   # 大きいカウントの場合、より小さなscale


                # dot plot
                try:
                    fig_dotplot2 = dc.plot_dotplot(
                        ORA_res.sort_values('Combined score', ascending=False).head(bar_top),
                        x='Combined score',
                        y='Term',
                        s='count',
                        c='FDR p-value',
                        scale = scale,
                        figsize = (bar_x_size, bar_y_size),  return_fig=True
                    )
                    st.pyplot(fig_dotplot2)
                    st.write("Combined score = -log10(P)")
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.write("Cannnot generate the dot plot")


                if use_upload == "Yes": # file_nameに追加
                    file_name_head = os.path.splitext(uploaded_file.name)[0]
                    file_name_add = file_name_head[:12] + "__"
                else:
                    file_name_add = ""
                if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                    if path == 'Homemade':
                        GO_file_name = '-'.join(GO_name)
                    else:
                        GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                    out_file_name = dc_temp_dir + "/ORA_" + file_name_add + GO_file_name + '.tsv'
                    zip_name = "ORA_" +file_name_add +  GO_file_name
                else:
                    out_file_name = dc_temp_dir + "/ORA_" + file_name_add + path + '.tsv'
                    zip_name = "ORA_" + file_name_add + path

                if ORA_mode == 'Cluster info file' and st.session_state.cluster_submitted:
                    out_file_name =  os.path.splitext(out_file_name)[0] + "." + "_".join(map(str, selected_cluster)) + '.tsv'

                st.write(out_file_name)
                st.write(zip_name)


                if st.button('Prepare ORA files to download'):
                    g.savefig( dc_temp_dir + "/" + bar_name, bbox_inches='tight')
                    dotplot.savefig( dc_temp_dir + "/dot_" + bar_name, bbox_inches='tight')
                    ORA_res.to_csv(out_file_name, sep = '\t')
                    shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)
                    with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                        btn = st.download_button(
                            label="Download Results of ORA",
                        data=fp,
                        file_name=zip_name + ".zip",
                        mime = "zip",
                        on_click = delete_file([out_file_name, dc_temp_dir + "/" + bar_name])#ダウンロードするとファイルを消す
                        )

            else:
                st.write("Click 'Run ORA analysis' to proceed.")