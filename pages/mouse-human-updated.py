import streamlit as st
import csv
import re
import pickle
import sys
import os
import pandas as pd
import itertools
import shutil
import numpy as np

# データディレクトリのパス
DATA_DIR = "/home/cellxgene/streamlit/data"

def load_m2h():
    import pickle
    with open(f"{DATA_DIR}/mouse2human.dic", mode='rb') as f:
        data = pickle.load(f)
    return data

def load_h2m():
    import pickle
    with open(f"{DATA_DIR}/human2mouse.dic", mode='rb') as f:
        data = pickle.load(f)
    return data

def mouse_human_conversion(dic, x):
    try:
        y = dic[x]
        return y
    except:
        return None

# 新しい変換テーブルをロード
@st.cache_data
def load_conversion_tables():
    """変換テーブルをロード"""
    tables = {}
    
    # HomoloGene + Ensembl Consensus (original)
    try:
        tables['consensus_original'] = pd.read_csv(f"{DATA_DIR}/consensus_orthologs_one_to_one.csv")
    except:
        st.warning("consensus_orthologs_one_to_one.csv not found")
        tables['consensus_original'] = None
    
    # HomoloGene + Ensembl Consensus (corrected)
    try:
        tables['consensus_corrected'] = pd.read_csv(f"{DATA_DIR}/consensus_orthologs_one_to_one_corrected.csv")
    except:
        st.warning("consensus_orthologs_one_to_one_corrected.csv not found")
        tables['consensus_corrected'] = None
    
    # NichenetR v1 (original)
    try:
        tables['nichenetr_v1_original'] = pd.read_csv("db/nichenetr.db/nichenetr_v1_original.csv")
    except:
        st.warning("db/nichenetr.db/nichenetr_v1_original.csv not found")
        tables['nichenetr_v1_original'] = None
    
    # NichenetR v1 (corrected)
    try:
        tables['nichenetr_v1_corrected'] = pd.read_csv("db/nichenetr.db/nichenetr_v1_corrected.csv")
    except:
        st.warning("db/nichenetr.db/nichenetr_v1_corrected.csv not found")
        tables['nichenetr_v1_corrected'] = None
    
    # NichenetR v2 (original)
    try:
        tables['nichenetr_v2_original'] = pd.read_csv("db/nichenetr.db/nichenetr_v2_original.csv")
    except:
        st.warning("db/nichenetr.db/nichenetr_v2_original.csv not found")
        tables['nichenetr_v2_original'] = None
    
    # NichenetR v2 (corrected)
    try:
        tables['nichenetr_v2_corrected'] = pd.read_csv("db/nichenetr.db/nichenetr_v2_corrected.csv")
    except:
        st.warning("db/nichenetr.db/nichenetr_v2_corrected.csv not found")
        tables['nichenetr_v2_corrected'] = None
    
    return tables

@st.cache_data
def convert_human_to_mouse_symbols(symbols, table_key='nichenetr_v1_original'):
    """ヒト遺伝子をマウス遺伝子に変換"""
    tables = load_conversion_tables()
    
    if table_key not in tables or tables[table_key] is None:
        st.error(f"Table {table_key} not available")
        return [np.nan] * len(symbols)
    
    geneinfo = tables[table_key]
    
    if table_key.startswith('consensus'):
        # Consensusテーブルの場合
        mapping_dict = dict(zip(geneinfo['human_symbol'], geneinfo['mouse_symbol']))
    else:
        # NichenetRテーブルの場合（従来のロジック）
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
        geneinfo_processed = pd.concat([
            geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
            geneinfo_ambiguous_solved
        ]).dropna()
        mapping_dict = dict(zip(geneinfo_processed['symbol'], geneinfo_processed['symbol_mouse']))
    
    converted_symbols = [mapping_dict.get(symbol, np.nan) for symbol in symbols]
    return converted_symbols

@st.cache_data
def convert_mouse_to_human_symbols(symbols, table_key='nichenetr_v1_original'):
    """マウス遺伝子をヒト遺伝子に変換"""
    tables = load_conversion_tables()
    
    if table_key not in tables or tables[table_key] is None:
        st.error(f"Table {table_key} not available")
        return [np.nan] * len(symbols)
    
    geneinfo = tables[table_key]
    
    if table_key.startswith('consensus'):
        # Consensusテーブルの場合
        mapping_dict = dict(zip(geneinfo['mouse_symbol'], geneinfo['human_symbol']))
    else:
        # NichenetRテーブルの場合（従来のロジック）
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
        geneinfo_processed = pd.concat([
            geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
            geneinfo_ambiguous_solved
        ]).dropna()
        mapping_dict = dict(zip(geneinfo_processed['symbol_mouse'], geneinfo_processed['symbol']))
    
    converted_symbols = [mapping_dict.get(symbol, np.nan) for symbol in symbols]
    return converted_symbols

# 初期化
if os.path.exists('res'):
    shutil.rmtree("res")
os.mkdir("res")

# UI
st.title("Mouse-Human Gene Conversion Tool")

species = st.radio(
    "Species from",
    ('Mouse','Human'), index=1)

method = st.radio(
    "Conversion method",
    (
        'in-house',
        'NichenetR v1 (original)', 
        'NichenetR v1 (corrected)',
        'NichenetR v2 (original)',
        'NichenetR v2 (corrected)',
        'Consensus (original)',
        'Consensus (corrected)'
    ), 
    index=0
)

# ヘルプセクション
with st.expander("📖 変換方法の詳細"):
    st.markdown("""
    ## 変換テーブルについて
    
    ### **in-house**
    - one to all orthologs mapping
    - 従来のピクルファイルベース
    
    ### **NichenetR**
    NCBIのEntrez Geneデータベースを基にした系統発生学的オルソログマッピング
    - **v1**: geneinfo_human（18,985エントリー）
    - **v2**: geneinfo_2022（122,797エントリー、より包括的）
    - **original**: 元のNichenetRデータ
    - **corrected**: 機能的に重要な遺伝子を修正
    
    ### **Consensus**
    HomoloGeneとEnsembl Comparaの両データベースで合意されたマッピング
    - **original**: 系統発生学的マッピングのみ
    - **corrected**: 機能的に重要な遺伝子を修正
    
    ## 修正された遺伝子（corrected版）
    
    研究で一般的に使われる機能的対応関係に基づいて以下の遺伝子を修正：
    
    | マウス | 修正前 | 修正後 | 根拠 |
    |--------|-------|--------|------|
    | **Ccl2** | CCL13 → **CCL2** | MCP-1として同じ機能 |
    | **Ccl3** | CCL3L3 → **CCL3** | MIP-1αとして同じ機能 |
    | **Cxcl1** | CXCL2 → **CXCL1** | KC/GRO-αとして同じ機能 |
    
    ### 修正の理由
    - **系統発生学的オルソログ**: 進化的に最も近い遺伝子
    - **機能的オルソログ**: 同じ生物学的機能を持つ遺伝子
    
    多くの免疫学・炎症研究では機能的対応関係（Ccl2=CCL2）が使われているため、
    研究目的に応じて適切なバージョンを選択してください。
    
    ### 推奨用途
    - **免疫学・炎症研究**: corrected版
    - **進化・系統解析**: original版
    - **一般的な発現解析**: corrected版（研究との整合性）
    """)

# ファイルアップロード
uploaded_file = st.file_uploader("Choose a CSV file", type="csv")

if uploaded_file is not None:
    # ファイル読み込み
    df = pd.read_csv(uploaded_file)
    st.write("Uploaded data preview:")
    st.write(df.head())
    
    # 遺伝子名カラムの選択
    gene_column = st.selectbox("Select gene name column:", df.columns)
    
    if st.button("Convert"):
        # 変換テーブルの選択
        table_key_map = {
            'in-house': None,
            'NichenetR v1 (original)': 'nichenetr_v1_original',
            'NichenetR v1 (corrected)': 'nichenetr_v1_corrected',
            'NichenetR v2 (original)': 'nichenetr_v2_original',
            'NichenetR v2 (corrected)': 'nichenetr_v2_corrected',
            'Consensus (original)': 'consensus_original',
            'Consensus (corrected)': 'consensus_corrected'
        }
        
        if method == 'in-house':
            # 従来のin-house方法
            if species == 'Mouse':
                dic = load_m2h()
                df['converted'] = df[gene_column].apply(lambda x: mouse_human_conversion(dic, x))
            else:
                dic = load_h2m()
                df['converted'] = df[gene_column].apply(lambda x: mouse_human_conversion(dic, x))
        else:
            # 新しい変換テーブル
            table_key = table_key_map[method]
            gene_list = df[gene_column].tolist()
            
            if species == 'Mouse':
                converted = convert_mouse_to_human_symbols(gene_list, table_key)
            else:
                converted = convert_human_to_mouse_symbols(gene_list, table_key)
            
            df['converted'] = converted
        
        # 結果表示
        st.write("Conversion results:")
        st.write(df)
        
        # 統計情報
        total_genes = len(df)
        converted_genes = df['converted'].notna().sum()
        conversion_rate = (converted_genes / total_genes) * 100
        
        st.write(f"**Conversion Statistics:**")
        st.write(f"- Total genes: {total_genes}")
        st.write(f"- Successfully converted: {converted_genes}")
        st.write(f"- Conversion rate: {conversion_rate:.1f}%")
        
        # ダウンロード
        csv = df.to_csv(index=False)
        st.download_button(
            label="Download converted data as CSV",
            data=csv,
            file_name=f'converted_{method.replace(" ", "_")}_{species.lower()}.csv',
            mime='text/csv'
        )

# サンプル変換テスト
with st.expander("🧪 サンプル変換テスト"):
    st.markdown("### 修正された遺伝子の変換例")
    
    sample_genes = {
        'Mouse': ['Ccl2', 'Ccl3', 'Cxcl1', 'Actb', 'Il6'],
        'Human': ['CCL2', 'CCL3', 'CXCL1', 'ACTB', 'IL6']
    }
    
    test_species = st.selectbox("Test species:", ['Mouse', 'Human'])
    test_method = st.selectbox("Test method:", [
        'NichenetR v1 (original)', 'NichenetR v1 (corrected)',
        'Consensus (original)', 'Consensus (corrected)'
    ])
    
    if st.button("Test conversion"):
        table_key = table_key_map[test_method]
        test_genes = sample_genes[test_species]
        
        if test_species == 'Mouse':
            results = convert_mouse_to_human_symbols(test_genes, table_key)
        else:
            results = convert_human_to_mouse_symbols(test_genes, table_key)
        
        test_df = pd.DataFrame({
            f'{test_species} Gene': test_genes,
            'Converted': results
        })
        
        st.write(f"**{test_method}** conversion results:")
        st.write(test_df)