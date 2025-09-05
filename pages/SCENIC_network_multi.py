import streamlit as st
import pandas as pd
import networkx as nx
from pyvis.network import Network
import tempfile
import os
from collections import defaultdict
import numpy as np
from typing import Dict, Set, Tuple, List, Optional
import streamlit.components.v1 as components
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import stats
import base64
import requests
import pandas as pd
import networkx as nx
from typing import Dict, List, Optional, Tuple
import requests
import time
import requests
import pandas as pd
import networkx as nx
from typing import Dict, List, Set, Tuple
import time
from collections import defaultdict



def get_string_interactions(proteins: List[str], species: str = "9606", score_threshold: float = 0.4) -> Dict[str, Dict[str, float]]:
    """STRING-dbからPPIデータを取得する"""
    base_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    
    # 種の情報を表示
    species_name = "Human" if species == "9606" else "Mouse"
    st.write(f"Fetching {species_name} protein interaction data from STRING-db")
    
    # シンボルの正規化
    normalized_proteins = [p.strip() for p in proteins]
    
    params = {
        "identifiers": "%0d".join(normalized_proteins),
        "species": species,
        "required_score": int(score_threshold * 1000),
        "network_type": "physical",
        "caller_identity": "your_app_name"
    }
    
    try:
        api_url = f"{base_url}/{output_format}/{method}"
        
        response = requests.post(api_url, data=params)
        response.raise_for_status()
        
        # レスポンスの解析
        interactions = defaultdict(dict)
        processed_lines = 0
        found_interactions = 0
        
        # 全ての相互作用を記録
        all_interactions = []
        
        for line in response.text.strip().split("\n"):
            if line:
                parts = line.split()
                if len(parts) >= 13:
                    try:
                        protein1, protein2 = parts[2], parts[3]
                        combined_score = float(parts[12])
                        all_interactions.append((protein1, protein2, combined_score))
                        processed_lines += 1
                    except (ValueError, IndexError) as e:
                        continue
        
        # デバッグ情報
        st.write("Debug: All found interactions:")
        for p1, p2, score in all_interactions[:10]:  # 最初の10個を表示
            st.write(f"{p1} - {p2}: {score}")
        
        # 相互作用の追加（大文字小文字を無視）
        for protein1, protein2, score in all_interactions:
            p1_upper = protein1.upper()
            p2_upper = protein2.upper()
            
            # 入力タンパク質リストとの照合（大文字小文字を無視）
            if p1_upper in [p.upper() for p in normalized_proteins] and \
               p2_upper in [p.upper() for p in normalized_proteins]:
                # オリジナルの大文字小文字を保持
                orig_p1 = next(p for p in normalized_proteins if p.upper() == p1_upper)
                orig_p2 = next(p for p in normalized_proteins if p.upper() == p2_upper)
                
                interactions[orig_p1][orig_p2] = score
                interactions[orig_p2][orig_p1] = score  # 双方向の相互作用を追加
                found_interactions += 1
        
        st.write(f"Processed {processed_lines} interactions from STRING-db")
        st.write(f"Found {found_interactions} relevant interactions between input proteins")
        
        if not interactions:
            st.warning(f"No {species_name} protein interactions found in STRING-db for the input proteins.")
        else:
            st.success(f"Successfully retrieved {species_name} protein interactions for {len(interactions)} proteins.")
            # 見つかった相互作用の例を表示
            st.write("Example interactions:")
            example_interactions = list(interactions.items())[:3]
            for protein1, partners in example_interactions:
                for protein2, score in list(partners.items())[:2]:
                    st.write(f"{protein1} - {protein2}: {score:.3f}")
        
        return dict(interactions)
        
    except requests.exceptions.RequestException as e:
        st.error(f"Error fetching {species_name} STRING-db data: {str(e)}")
        return {}
    except Exception as e:
        st.error(f"Unexpected error processing {species_name} STRING-db data: {str(e)}")
        return {}

def convert_to_string_symbols(proteins: List[str], species: str = "9606") -> Dict[str, str]:
    """
    遺伝子シンボルをSTRING-db形式に変換
    """
    # シンボルの正規化
    symbol_map = {}
    for protein in proteins:
        # 基本的なクリーニング（大文字変換は必要に応じて）
        cleaned = protein.strip()
        symbol_map[protein] = cleaned
    
    return symbol_map

def add_string_db_options(st, species: str) -> Dict:
    """STRING-db関連のオプションをUIに追加"""
    st.subheader("STRING-db Integration")
    
    with st.expander("STRING-db Settings", expanded=True):
        # STRING-dbフィルタリングの有効化
        enable_string = st.checkbox(
            "Enable STRING-db Filtering",
            help="Filter network edges based on STRING protein-protein interactions"
        )
        
        if enable_string:
            # スコア閾値の設定
            score_threshold = st.slider(
                "Minimum STRING interaction score",
                min_value=0.0,
                max_value=1.0,
                value=0.4,
                step=0.1,
                help="Filter interactions based on STRING combined score"
            )
            
            # 生物種の表示
            species_id = "9606" if species == "human" else "10090"
            st.info(f"Using STRING-db data for {species.capitalize()} (Taxonomy ID: {species_id})")
            
            return {
                "enable": enable_string,
                "score_threshold": score_threshold,
                "species_id": species_id
            }
    
    return {
        "enable": False,
        "score_threshold": 0.4,
        "species_id": "9606"
    }

def filter_network_by_string(
    tf_network: nx.DiGraph,
    string_data: Dict[str, Dict[str, float]],
    score_threshold: float = 0.0
) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """
    ネットワーク内の全エッジをSTRING-dbデータで評価
    
    Args:
        tf_network: 元のネットワーク
        string_data: STRING-dbからの相互作用データ
        score_threshold: 最小スコア閾値
    """
    # 新しいネットワークを作成
    filtered_network = nx.DiGraph()
    filtered_network.add_nodes_from(tf_network.nodes())
    
    # ネットワーク情報
    st.write(f"Original network:")
    st.write(f"Nodes: {tf_network.number_of_nodes()}")
    st.write(f"Edges: {tf_network.number_of_edges()}")
    
    # 元のエッジと相互作用情報を表示
    original_edges = list(tf_network.edges())
 #   st.write("\nEvaluating original network edges:")
    for i, (source, target) in enumerate(original_edges, 1):
     #   st.write(f"\nEdge {i}: {source} -> {target}")
        
        # STRING-dbでの相互作用を確認
        score = None
        if source in string_data and target in string_data[source]:
            score = string_data[source][target]
        elif target in string_data and source in string_data[target]:
            score = string_data[target][source]
        
        if score is not None:
          #  st.write(f"  STRING-db interaction found with score: {score}")
            if score >= score_threshold:
                filtered_network.add_edge(source, target, weight=score)
             #   st.write(f"  Edge kept (score >= threshold)")
          #  else:
           #     st.write(f"  Edge removed (score < threshold)")
       # else:
        #    st.write(f"  No STRING-db interaction found")
    
    # 削除されたエッジを特定
    kept_edges = list(filtered_network.edges())
    removed_edges = [(u, v) for u, v in original_edges if (u, v) not in kept_edges]
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    # 結果の要約
    st.write("\nFiltering Results:")
    st.write(f"Edges in filtered network: {filtered_network.number_of_edges()}")
    st.write(f"Edges removed: {len(removed_edges)}")
    st.write(f"Nodes in filtered network: {filtered_network.number_of_nodes()}")
    st.write(f"Nodes removed: {len(isolated_nodes)}")
    
    # 保持されたエッジの詳細
#    if kept_edges:
      #  st.write("\nKept edges (sorted by STRING-db score):")
 #       edge_scores = []
  #      for u, v in kept_edges:
   #         score = filtered_network[u][v]['weight']
    #        edge_scores.append((u, v, score))
     #   
      #  for u, v, score in sorted(edge_scores, key=lambda x: x[2], reverse=True):
       #     st.write(f"{u} -> {v}: {score}")
    
    # 削除されたエッジの詳細
  #  if removed_edges:
     #   st.write("\nRemoved edges:")
   #     for u, v in removed_edges:
    #        score = None
     #       if u in string_data and v in string_data[u]:
      #          score = string_data[u][v]
       #     elif v in string_data and u in string_data[v]:
        #        score = string_data[v][u]
            
      #      if score is not None:
      #          st.write(f"{u} -> {v}: score {score} (below threshold)")
      #      else:
      #          st.write(f"{u} -> {v}: no STRING-db interaction")
    
    return filtered_network, removed_edges

def show_string_db_stats(st, original_network, filtered_network, removed_edges, string_data):
    """
    STRING-dbフィルタリングの統計情報を表示
    
    Args:
        st: Streamlitインスタンス
        original_network: 元のネットワーク
        filtered_network: フィルタリング後のネットワーク
        removed_edges: 除去されたエッジのリスト
        string_data: STRING-dbの相互作用データ
    """
    st.subheader("STRING-db Filtering Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Original Edges", original_network.number_of_edges())
    
    with col2:
        st.metric("Filtered Edges", filtered_network.number_of_edges())
    
    with col3:
        st.metric("Removed Edges", len(removed_edges))
    
    if removed_edges:
        with st.expander("View Removed Edges", expanded=False):
            removed_df = pd.DataFrame(removed_edges, columns=['Source', 'Target'])
            st.dataframe(removed_df)
            
            # CSVダウンロードボタンの追加
            csv = removed_df.to_csv(index=False)
            st.download_button(
                "Download removed edges as CSV",
                csv,
                "string_db_removed_edges.csv",
                "text/csv",
                key='download-string-removed-edges'
            )

    # STRING-db相互作用の統計
    if string_data:
        st.subheader("STRING-db Interaction Statistics")
        
        # 相互作用スコアの収集
        scores = []
        for protein1, interactions in string_data.items():
            scores.extend(interactions.values())
        
        if scores:
            stats_dict = {
                "Mean Score": np.mean(scores),
                "Median Score": np.median(scores),
                "Min Score": np.min(scores),
                "Max Score": np.max(scores),
                "Std Dev": np.std(scores)
            }
            
            col1, col2 = st.columns(2)
            
            # 統計情報の表示
            with col1:
                for key, value in stats_dict.items():
                    st.metric(key, f"{value:.3f}")
            
            # スコア分布のヒストグラム
            with col2:
                fig, ax = plt.subplots()
                ax.hist(scores, bins=20, edgecolor='black')
                ax.set_xlabel("STRING Interaction Score")
                ax.set_ylabel("Frequency")
                ax.set_title("Distribution of STRING Interaction Scores")
                plt.tight_layout()
                st.pyplot(fig)
                plt.close()


@st.cache_data
def load_trrust_data(species: str) -> Dict[str, List[Tuple[str, str, str]]]:
    """TRRUSTデータの読み込み
    
    Args:
        species: 'mouse' または 'human'
    
    Returns:
        Dict[str, List[Tuple[str, str, str]]]: 転写因子をキーとし、
        (ターゲット遺伝子, 相互作用タイプ, PMID)のリストを値とする辞書
    """
    trrust_data = defaultdict(list)
    
    # ファイルパスの決定
    if species == 'mouse':
        filepath = 'db/TRRUST/trrust_rawdata.mouse.tsv'
    elif species == 'human':
        filepath = 'db/TRRUST/trrust_rawdata.human.tsv'
    else:
        st.warning("Invalid species. Choose 'mouse' or 'human'.")
        return trrust_data
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                # 改行文字や余分な空白を除去
                line = line.strip()
                if not line:  # 空行をスキップ
                    continue
                
                parts = line.split('\t')
                if len(parts) < 4:  # 不完全な行をスキップ
                    continue
                
                tf, target, interaction_type, pmid = parts
                trrust_data[tf].append((target, interaction_type, pmid))
        
        # デバッグ出力
        st.write(f"TRRUST data loaded. Total TFs: {len(trrust_data)}")
    except FileNotFoundError:
        st.error(f"TRRUST data file not found: {filepath}")
    except Exception as e:
        st.error(f"Error reading TRRUST data: {e}")
    
    return trrust_data

def add_trrust_options(st, species="human"):
    """TRRUST関連のオプションをUIに追加"""
    st.subheader("TRRUST Database Integration")
    
    with st.expander("TRRUST Settings", expanded=True):
        # TRRUSTフィルタリングの有効化
        enable_trrust = st.checkbox(
            "Enable TRRUST Filtering", 
            help="Filter network edges based on TRRUST interaction data"
        )
        
        if enable_trrust:
            col1, col2 = st.columns(2)
            
            with col1:
                # 相互作用タイプの選択
                interaction_types = st.multiselect(
                    "Select Interaction Types",
                    options=["Activation", "Repression", "Unknown"],
                    default=["Activation", "Repression", "Unknown"],
                    help="Choose interaction types to include"
                )
            
            with col2:
                # 最小PMID数の設定
                min_pmid_count = st.slider(
                    "Minimum PMID Count",
                    min_value=1,
                    max_value=10,
                    value=1,
                    help="Minimum number of unique publications supporting the interaction"
                )
            
            # 生物種の表示
            st.info(f"Using {species.capitalize()} TRRUST data")
            
            return {
                "enable": enable_trrust,
                "interaction_types": interaction_types,
                "min_pmid_count": min_pmid_count
            }
        
    return {
        "enable": False,
        "interaction_types": [],
        "min_pmid_count": 1
    }



def filter_network_by_trrust(
    tf_network: nx.DiGraph, 
    trrust_data: Dict[str, List[Tuple[str, str, str]]],
    interaction_types: List[str],
    min_pmid_count: int
) -> Tuple[nx.DiGraph, List[Tuple[str, str]], List[Dict]]:
    """TRRUSTデータに基づいてネットワークをフィルタリング"""
    filtered_network = tf_network.copy()
    removed_edges = []
    trrust_edge_details = []  # TRRUSTデータのエッジ詳細を保存するリスト
    
    for source, target in list(tf_network.edges()):
        # 大文字小文字の違いを考慮するため、大文字・小文字両方でチェック
        source_variations = [source, source.lower(), source.upper(), 
                             source.capitalize(), source.title()]
        target_variations = [target, target.lower(), target.upper(), 
                             target.capitalize(), target.title()]
        
        # いずれかの大文字小文字のバリエーションで相互作用を探索
        trrust_interactions = []
        for src_var in source_variations:
            if src_var in trrust_data:
                trrust_interactions.extend([
                    (t, it, pmid) for (t, it, pmid) in trrust_data[src_var] 
                    if t in target_variations
                ])
        
        if trrust_interactions:
            # 選択された相互作用タイプと最小PMID数を満たすか確認
            valid_interactions = [
                (t, it, pmid) for (t, it, pmid) in trrust_interactions
                if it in interaction_types
            ]
            
            if valid_interactions:
                # 重複を除いたPMIDの数を計算
                unique_pmids = set()
                for _, _, pmid_str in valid_interactions:
                    pmids = pmid_str.split(';')
                    unique_pmids.update(pmids)
                
                # PMIDの数が最小要件を満たす場合
                if len(unique_pmids) >= min_pmid_count:
                    # エッジの詳細情報を保存
                    edge_detail = {
                        'source': source,
                        'target': target,
                        'interactions': valid_interactions,
                        'unique_pmids': list(unique_pmids),
                        'pmid_count': len(unique_pmids)
                    }
                    trrust_edge_details.append(edge_detail)
                else:
                    filtered_network.remove_edge(source, target)
                    removed_edges.append((source, target))
            else:
                filtered_network.remove_edge(source, target)
                removed_edges.append((source, target))
        else:
            filtered_network.remove_edge(source, target)
            removed_edges.append((source, target))
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    # TRRUSTデータのエッジ詳細をDataFrameで表示
    if trrust_edge_details:
        trrust_edges_df = pd.DataFrame(trrust_edge_details)
        st.write("### TRRUST Network Edges")
        st.write("Edges found in TRRUST database:")
        st.dataframe(trrust_edges_df)
    
    
    return filtered_network, removed_edges


def show_trrust_stats(st, original_network, filtered_network, removed_edges, trrust_data):
    """TRRUSTフィルタリングの統計情報を表示"""
    st.subheader("TRRUST Filtering Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Original Edges", original_network.number_of_edges())
        
    with col2:
        st.metric("Filtered Edges", filtered_network.number_of_edges())
        
    with col3:
        st.metric("Removed Edges", len(removed_edges))
    
    if removed_edges:
        with st.expander("View Removed Edges", expanded=False):
            removed_df = pd.DataFrame(removed_edges, columns=['Source', 'Target'])
            st.dataframe(removed_df)
            
            # CSVダウンロードボタンの追加
            csv = removed_df.to_csv(index=False)
            st.download_button(
                "Download removed edges as CSV",
                csv,
                "trrust_removed_edges.csv",
                "text/csv",
                key='download-trrust-removed-edges'
            )
def format_tf_name_for_chip_atlas(tf_name: str, species: str) -> str:
    """ChIP-Atlas用にTF名をフォーマット
    
    Args:
        tf_name: オリジナルのTF名
        species: 生物種 ("human" or "mouse")
        
    Returns:
        str: フォーマットされたTF名
    """
    # マウスの場合、最初の文字を大文字に、残りを小文字に
    if species == "mouse":
        if tf_name.startswith("Zfp") or tf_name.startswith("ZFP"):
            # Zfpの場合は特別処理
            return f"{tf_name[0].upper()}{tf_name[1:].lower()}"
        return f"{tf_name[0].upper()}{tf_name[1:].lower()}"
    
    # ヒトの場合、すべて大文字に
    return tf_name.upper()

def get_chip_atlas_data(tf: str, species: str = "mouse", distance: int = 1000) -> pd.DataFrame:
    """ChIP-Atlas APIからターゲット遺伝子データを取得"""
    genome = "hg38" if species == "human" else "mm10"
    distance_kb = distance // 1000
    
    url = f"https://chip-atlas.dbcls.jp/data/{genome}/target/{tf}.{distance_kb}.tsv"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
    #    st.write(f"Request URL: {url}")
    #    st.write(f"Response status: {response.status_code}")
        
        if response.text:
            from io import StringIO
            df = pd.read_csv(StringIO(response.text), sep='\t')
            
            # デバッグ情報を表示
          #  st.write("Column names:", df.columns.tolist())
          #  st.write("First few rows:")
          #  st.write(df.head())
            
            # カラム名を標準化
            column_mapping = {
                'Gene ID': 'gene_id',
                'Gene Symbol': 'target_gene',
                'Distance': 'distance',
                'Average': 'binding_score'
            }
            df = df.rename(columns=column_mapping)
            
            return df
        else:
            st.warning(f"No data found for {tf}")
            return pd.DataFrame()
            
    except requests.exceptions.RequestException as e:
        st.warning(f"Error retrieving data for {tf}: {str(e)}")
        return pd.DataFrame()
    except pd.errors.EmptyDataError:
        st.warning(f"No data in response for {tf}")
        return pd.DataFrame()

def filter_network_by_chip_atlas(tf_network: nx.DiGraph, 
                               chip_data: Dict[str, pd.DataFrame], 
                               min_score: float = 0) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """ChIP-Atlasデータに基づいてネットワークをフィルタリング"""
    filtered_network = tf_network.copy()
    removed_edges = []
    
    for source, target in list(tf_network.edges()):
        if source in chip_data:
            df = chip_data[source]
            if not df.empty:
                # カラム名の存在確認
                if 'Target_genes' not in df.columns:
                    st.warning(f"Missing 'Target_genes' column in data for {source}")
                    st.write("Available columns:", df.columns.tolist())
                    continue
                    
                if target in df['Target_genes'].values:
                    score = df[df['Target_genes'] == target].iloc[0,1]
                    if score < min_score:
                        filtered_network.remove_edge(source, target)
                        removed_edges.append((source, target))
                else:
                    filtered_network.remove_edge(source, target)
                    removed_edges.append((source, target))
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    return filtered_network, removed_edges

def get_network_chip_data(network_tfs: List[str], species: str = "mouse", 
                         distance: int = 1000) -> Dict[str, pd.DataFrame]:
    """ネットワーク内の全TFのChIP-Atlasデータを取得"""
    chip_data = {}
    
    progress_text = st.empty()
    progress_bar = st.progress(0)
    
    for i, tf in enumerate(network_tfs):
        progress_text.write(f"Fetching data for {tf} ({i+1}/{len(network_tfs)})")
        
        df = get_chip_atlas_data(tf, species, distance)
        if not df.empty:
            chip_data[tf] = df
            
        progress_bar.progress((i + 1) / len(network_tfs))
        time.sleep(1)  # APIの負荷を考慮
    
    progress_text.empty()
    progress_bar.empty()
    
    return chip_data


def get_edge_chip_scores(tf_network: nx.DiGraph, 
                        chip_data: Dict[str, pd.DataFrame]) -> Dict[Tuple[str, str], float]:
    """ネットワークの各エッジのChIP-Atlasスコアを取得"""
    edge_scores = {}
    
    for source, target in tf_network.edges():
        if source in chip_data:
            df = chip_data[source]
            if not df.empty and target in df['target_gene'].values:
                score = df[df['target_gene'] == target]['binding_score'].iloc[0]
                edge_scores[(source, target)] = score
            else:
                edge_scores[(source, target)] = 0.0
    
    return edge_scores

def determine_species(regulon_file) -> str:
    """regulonファイルの内容から生物種を判定"""
    content = regulon_file.getvalue().decode('utf-8').splitlines()
    
    is_mostly_uppercase = all(
        len(line.split('\t')[0].split('(')[0].strip()) > 0 and 
        line.split('\t')[0].split('(')[0].strip().isupper() 
        for line in content[1:] if len(line.split('\t')) >= 2
    )
    
    return "human" if is_mostly_uppercase else "mouse"


def add_chip_atlas_options(st, species="human"):
    """ChIP-Atlas関連のオプションをUIに追加"""
    st.subheader("ChIP-Atlas Integration")
    
    with st.expander("ChIP-Atlas Settings", expanded=True):
        # ChIP-Atlasフィルタリングの有効化
        enable_chip_atlas = st.checkbox(
            "Enable ChIP-Atlas Filtering", 
            help="Filter network edges based on ChIP-Atlas binding data"
        )
        
        if enable_chip_atlas:
            col1, col2 = st.columns(2)
            
            with col1:
                # TSSからの距離の選択
                tss_distance = st.selectbox(
                    "Distance from TSS",
                    options=[1000, 5000, 10000],
                    format_func=lambda x: f"± {x//1000} kb",
                    help="Select the distance from Transcription Start Site",
                    index=1
                )
                
            with col2:
                # 最小バインディングスコアの設定
                min_binding_score = st.slider(
                    "Minimum Average Binding Score",
                    min_value=0,
                    max_value=1000,
                    value=13,
                    step=1,
                    help="Filter edges based on ChIP-Atlas binding score. -10*Log10[MACS2 Q-value]. -10*log10(0.05) ≈ 13)."
                )
            
            # 生物種の表示
            st.info(f"Using {species.capitalize()} genome ({('hg38' if species == 'human' else 'mm10')})")
            
            return {
                "enable": enable_chip_atlas,
                "distance": tss_distance,
                "min_score": min_binding_score
            }
        
    return {
        "enable": False,
        "distance": 5000,
        "min_score": 50
    }

def show_chip_atlas_stats(st, original_network, filtered_network, removed_edges, chip_data):
    """ChIP-Atlasフィルタリングの統計情報を表示"""
    st.subheader("ChIP-Atlas Filtering Results")
    
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.metric("Original Edges", original_network.number_of_edges())
        
    with col2:
        st.metric("Filtered Edges", filtered_network.number_of_edges())
        
    with col3:
        st.metric("Removed Edges", len(removed_edges))
    
    if removed_edges:
        with st.expander("View Removed Edges", expanded=False):
            removed_df = pd.DataFrame(removed_edges, columns=['Source', 'Target'])
            st.dataframe(removed_df)
            
            # CSVダウンロードボタンの追加
            csv = removed_df.to_csv(index=False)
            st.download_button(
                "Download removed edges as CSV",
                csv,
                "removed_edges.csv",
                "text/csv",
                key='download-removed-edges'
            )


@st.cache_data
def load_tfs_list(file) -> Set[str]:
    """
    Regulonデータの大文字/小文字の状況に基づいて適切なTFリストを読み込む
    
    Args:
        file: アップロードされたregulonファイル
    
    Returns:
        Set[str]: TFのリスト
    """
    # ファイルの内容を読み込み
    content = file.getvalue().decode('utf-8').splitlines()
    
    # 最初の行をスキップしてデータを処理
    # 遺伝子名の大文字/小文字パターンを確認
    is_mostly_uppercase = all(
        len(line.split('\t')[0].split('(')[0].strip()) > 0 and 
        line.split('\t')[0].split('(')[0].strip().isupper() 
        for line in content[1:] if len(line.split('\t')) >= 2
    )
    
    # TFリストを読み込む
    try:
        if is_mostly_uppercase:
            # ヒト(hg38)のTFリストを読み込む
            filepath = 'db/SCENIC/allTFs_hg38.txt'
        else:
            # マウス(mm)のTFリストを読み込む
            filepath = 'db/SCENIC/allTFs_mm.txt'
        
        with open(filepath, 'r') as f:
            return {line.strip() for line in f if line.strip()}
    
    except FileNotFoundError:
        st.warning(f"TFリストファイル {filepath} が見つかりませんでした。")
        return set()

@st.cache_data
def load_regulon_data(file) -> Tuple[Dict[str, Set[str]], Set[str]]:
    """
    regulonデータの読み込み
    Returns:
        Tuple[Dict[str, Set[str]], Set[str]]: (regulon_data, all_genes)
        - regulon_data: TFとその制御遺伝子の辞書
        - all_genes: 全ての遺伝子（TFと標的遺伝子）のセット
    """
    regulon_data = defaultdict(set)
    all_genes = set()
    
    content = file.getvalue().decode('utf-8').splitlines()
    for line in content[1:]:  # Skip header
        row = line.split('\t')
        if len(row) >= 2:
            tf = row[0].split('(')[0].strip()
            genes = {gene.strip() for gene in row[1].split(',')}
            regulon_data[tf] = genes
            all_genes.add(tf)
            all_genes.update(genes)
    
    return regulon_data, all_genes

@st.cache_data
def load_csi_data(file) -> Dict[str, Dict[str, float]]:
    """CSIデータの読み込み"""
    csi_data = defaultdict(dict)
    content = file.getvalue().decode('utf-8').splitlines()
    for line in content[1:]:  # Skip header
        row = line.split(',')
        if len(row) >= 3:
            # 1列目と2列目から転写因子名を抽出 (数字以降を除去)
            tf1 = row[0].split('(')[0].split()[0].strip()
            tf2 = row[1].split('(')[0].split()[0].strip()
            try:
                csi = float(row[2])
                csi_data[tf1][tf2] = csi
            except ValueError:
                continue
    return csi_data


def get_connected_tfs(regulon_data: Dict[str, Set[str]], 
                     csi_data: Dict[str, Dict[str, float]], 
                     source_tf: str,
                     network_mode: bool,
                     csi_threshold: float = None,
                     max_tfs: int = None) -> List[Tuple[str, float, str]]:
    """TFに接続する全てのTFとそのCSI値を取得"""
    connected_tfs = []
    
    # 追加のデバッグ出力
    print(f"\nget_connected_tfs for {source_tf}")
    print(f"network_mode: {network_mode}")
    print(f"csi_threshold: {csi_threshold}")
    print(f"source_tf in regulon_data: {source_tf in regulon_data}")
    print(f"Source TF downstream genes: {regulon_data.get(source_tf, 'No data')}")
    
    # 上流TFの取得
    for tf, genes in regulon_data.items():
        if source_tf in genes:
            weight = csi_data[tf].get(source_tf, 0.0)
            if csi_threshold is None or weight >= csi_threshold:
                connected_tfs.append((tf, weight, 'upstream'))
    
    # network_modeの場合は下流も取得
    if network_mode and source_tf in regulon_data:
        for target in regulon_data[source_tf]:
            if target in regulon_data:  # TFの場合のみ
                weight = csi_data[source_tf].get(target, 0.0)
                print(f"Potential downstream TF: {target}, CSI: {weight}")
                if csi_threshold is None or weight >= csi_threshold:
                    connected_tfs.append((target, weight, 'downstream'))
    
    # CSIの値でソートして上位を選択
    connected_tfs.sort(key=lambda x: x[1], reverse=True)
    if max_tfs is not None:
        connected_tfs = connected_tfs[:max_tfs]
    
    print(f"Connected TFs: {connected_tfs}")
    return connected_tfs


def calculate_weighted_centrality(tf_network: nx.DiGraph) -> Tuple[Dict, Dict, Dict, Dict]:
    """重み付きネットワークの各種中心性指標を計算"""
    try:
        # 重み付き次数中心性
        weighted_degree = {}
        for node in tf_network.nodes():
            in_weight = sum(d['weight'] for u, v, d in tf_network.in_edges(node, data=True))
            out_weight = sum(d['weight'] for u, v, d in tf_network.out_edges(node, data=True))
            weighted_degree[node] = (in_weight + out_weight) / (2 * tf_network.number_of_nodes() - 2)
        
        # 重み付きbetweenness centrality
        betweenness_centrality = nx.betweenness_centrality(tf_network, weight='weight')
        
        # 重み付きPageRank
        pagerank = nx.pagerank(tf_network, weight='weight')
        
        # 重み付き固有ベクトル中心性
        eigenvector_centrality = nx.eigenvector_centrality(tf_network, weight='weight')
        
    except Exception as e:
        print(f"Error in centrality calculation: {str(e)}")
        # フォールバック：重みなし中心性
        weighted_degree = nx.degree_centrality(tf_network)
        betweenness_centrality = nx.betweenness_centrality(tf_network)
        pagerank = nx.pagerank(tf_network)
        eigenvector_centrality = pagerank
    
    return weighted_degree, betweenness_centrality, pagerank, eigenvector_centrality


def identify_hubs_weighted(tf_network: nx.DiGraph, weight_threshold: float = 0.9) -> Set[str]:
    """重み付きネットワークでのハブの同定"""
    # 重み付き次数の計算
    weighted_degrees = {}
    for node in tf_network.nodes():
        in_weight = sum(d['weight'] for u, v, d in tf_network.in_edges(node, data=True))
        out_weight = sum(d['weight'] for u, v, d in tf_network.out_edges(node, data=True))
        weighted_degrees[node] = in_weight + out_weight
    
    threshold = np.percentile(list(weighted_degrees.values()), weight_threshold*100)
    hubs = {node for node, degree in weighted_degrees.items() if degree >= threshold}
    return hubs


def identify_bottlenecks_weighted(tf_network: nx.DiGraph, betweenness_threshold: float = 0.9) -> Set[str]:
    """重み付きネットワークでのボトルネックの同定"""
    betweenness = nx.betweenness_centrality(tf_network, weight='weight')
    threshold = np.percentile(list(betweenness.values()), betweenness_threshold*100)
    bottlenecks = {node for node, bc in betweenness.items() if bc >= threshold}
    return bottlenecks


def filter_network_by_csi(tf_network: nx.DiGraph, csi_threshold: float = 0.0) -> nx.DiGraph:
    """CSIの閾値に基づいてネットワークをフィルタリング"""
    filtered_network = nx.DiGraph()
    
    for u, v, data in tf_network.edges(data=True):
        if data['weight'] >= csi_threshold:
            filtered_network.add_edge(u, v, **data)
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    return filtered_network


def analyze_csi_distribution(tf_network: nx.DiGraph):
    """CSIの分布を分析"""
    edge_weights = [d['weight'] for (u, v, d) in tf_network.edges(data=True)]
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    
    # CSIの分布
    ax1.hist(edge_weights, bins=30, alpha=0.7)
    ax1.set_title('CSI Distribution')
    ax1.set_xlabel('CSI')
    ax1.set_ylabel('Frequency')
    
    # Box plot
    ax2.boxplot(edge_weights)
    ax2.set_title('CSI Box Plot')
    ax2.set_ylabel('CSI')
    
    plt.tight_layout()
    st.pyplot(fig)
    
    # 基本統計量を返す
    stats_dict = {
        "Mean": np.mean(edge_weights),
        "Median": np.median(edge_weights),
        "Std": np.std(edge_weights),
        "Min": np.min(edge_weights),
        "Max": np.max(edge_weights)
    }
    if len(edge_weights) > 0:
        stats_dict["Skewness"] = stats.skew(edge_weights)
    
    return stats_dict

def get_edge_color(source_node, node_categories, color_scheme, target, hubs=None):
   """エッジの色を取得する関数
   
   Args:
       source_node (str): エッジの起点となるノード
       node_categories (dict): カテゴリーごとのノードのセット
       color_scheme (dict): カテゴリーごとの色の辞書
       target (str): ターゲット遺伝子の名前
       hubs (set, optional): hubノードのセット. Defaults to None.
   
   Returns:
       str: 16進数のカラーコード
   """    
   # hub nodeのチェックを最初に行う
   if hubs and source_node in hubs:
       return color_scheme['hub_node']
       
   if source_node == target:
       return color_scheme['target']
   elif source_node in node_categories['upstream_1']:
       return color_scheme['upstream_1']
   elif source_node in node_categories['upstream_2']:
       return color_scheme['upstream_2']
   elif source_node in node_categories['downstream_1']:
       return color_scheme['downstream_1']
   elif source_node in node_categories['downstream_2']:
       return color_scheme['downstream_2']
   
   return color_scheme['target_genes']  # デフォルト色

def build_expanded_tf_network(
    regulon_data: Dict[str, Set[str]], 
    csi_data: Dict[str, Dict[str, float]], 
    all_tfs: Set[str],
    target_genes: List[str],
    max_upstream: int,
    max_downstream: int,
    network_mode: bool,
    csi_threshold: float = None,
    max_tfs_per_level: int = None
) -> Tuple[nx.DiGraph, Dict[str, Set[str]]]:
    """複数のターゲット遺伝子に対する拡張ネットワークを構築"""
    tf_network = nx.DiGraph()
    regulators_by_level = {
        'upstream_1': set(),
        'upstream_2': set(),
        'downstream_1': set(),
        'downstream_2': set(),
        'additional_tfs': set(),
        'target_genes': set(target_genes)  # ターゲット遺伝子を記録
    }
    
    def get_connected_tfs_for_level(current_tf: str, direction: str) -> List[Tuple[str, float, str]]:
        """特定の方向の接続TFを取得"""
        connected_tfs = []
        
        # 上流と下流を探索
        for tf, genes in regulon_data.items():
            # 上流の探索
            if direction in ['upstream', 'both']:
                if current_tf in genes and tf in all_tfs:
                    weight = csi_data[tf].get(current_tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'upstream'))
            
            # 下流の探索
            if direction in ['downstream', 'both']:
                if current_tf in regulon_data and tf in regulon_data[current_tf] and tf in all_tfs:
                    weight = csi_data[current_tf].get(tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'downstream'))
        
        # CSIの値でソート
        connected_tfs.sort(key=lambda x: x[1], reverse=True)
        if max_tfs_per_level is not None:
            connected_tfs = connected_tfs[:max_tfs_per_level]
        
        return connected_tfs

    # 各ターゲット遺伝子に対して処理
    for target_gene in target_genes:
        # 上流の探索
        current_level_tfs = {target_gene}
        for level in range(1, max_upstream + 1):
            next_level_tfs = set()
            level_key = f'upstream_{level}'
            
            for current_tf in current_level_tfs:
                direction = 'upstream' if not network_mode else 'both'
                connected_tfs = get_connected_tfs_for_level(current_tf, direction)
                
                for tf, weight, conn_direction in connected_tfs:
                    if tf not in all_tfs:
                        regulators_by_level['additional_tfs'].add(tf)
                    else:
                        regulators_by_level[level_key].add(tf)
                    
                    if conn_direction == 'upstream':
                        tf_network.add_edge(tf, current_tf, weight=weight)
                    else:
                        tf_network.add_edge(current_tf, tf, weight=weight)
                    
                    next_level_tfs.add(tf)
            
            current_level_tfs = next_level_tfs
        
        # 下流の探索
        if target_gene in regulon_data:
            current_level_tfs = {target_gene}
            for level in range(1, max_downstream + 1):
                next_level_tfs = set()
                level_key = f'downstream_{level}'
                
                for current_tf in current_level_tfs:
                    direction = 'downstream' if not network_mode else 'both'
                    connected_tfs = get_connected_tfs_for_level(current_tf, direction)
                    
                    for tf, weight, conn_direction in connected_tfs:
                        if tf not in all_tfs:
                            regulators_by_level['additional_tfs'].add(tf)
                        else:
                            regulators_by_level[level_key].add(tf)
                        
                        if conn_direction == 'downstream':
                            tf_network.add_edge(current_tf, tf, weight=weight)
                        else:
                            tf_network.add_edge(tf, current_tf, weight=weight)
                        
                        next_level_tfs.add(tf)
                
                current_level_tfs = next_level_tfs
    
    return tf_network, regulators_by_level

def cleanup_network(network: nx.DiGraph, target_genes: List[str]) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """
    いずれかのターゲット遺伝子に接続されていないサブネットワークを削除
    
    Args:
        network: 元のネットワーク
        target_genes: ターゲット遺伝子のリスト
    
    Returns:
        Tuple[nx.DiGraph, List[Tuple[str, str]]]: 
            - クリーニング後のネットワーク
            - 除外されたエッジのリスト
    """
    cleaned_network = network.copy()
    
    # 有向グラフを無向グラフに変換して到達可能性をチェック
    undirected = network.to_undirected()
    
    # いずれかのターゲット遺伝子から到達可能なノードを取得
    reachable_nodes = set()
    for target_gene in target_genes:
        reachable_nodes.update(nx.node_connected_component(undirected, target_gene))
    
    # 到達不可能なノードを特定
    unreachable_nodes = set(network.nodes()) - reachable_nodes
    
    # 除外されるエッジを記録
    removed_edges = []
    for u, v in network.edges():
        if u in unreachable_nodes or v in unreachable_nodes:
            removed_edges.append((u, v))
    
    # 到達不可能なノードを削除
    cleaned_network.remove_nodes_from(unreachable_nodes)
    
    # 除外情報の表示
    if removed_edges:
        st.markdown("### Removing disconnected subnetworks:")
        st.write(f"Nodes removed: {len(unreachable_nodes)}")
        st.write(f"Edges removed: {len(removed_edges)}")
        
        removed_df = pd.DataFrame(removed_edges, columns=['Source', 'Target'])
        
        with st.expander("View removed disconnected edges", expanded=False):
            st.write("Edges in disconnected subnetworks:")
            st.dataframe(removed_df)
            
            csv = removed_df.to_csv(index=False)
            st.download_button(
                "Download removed edges as CSV",
                csv,
                "disconnected_edges.csv",
                "text/csv",
                key='download-disconnected-edges'
            )
    
    return cleaned_network, removed_edges


def get_node_color(node: str, 
                  regulators_by_level: Dict[str, Set[str]], 
                  color_map: Dict[str, str], 
                  target_genes: List[str], 
                  hubs: Optional[Set[str]] = None) -> str:
    """ノードの色を決定する
    
    Args:
        node: 対象ノード
        regulators_by_level: レベルごとの制御因子
        color_map: 色マップ
        target_genes: ターゲット遺伝子のリスト
        hubs: ハブノードのセット（オプション）
    
    Returns:
        str: 色を表す16進数コード
    """
    # ターゲット遺伝子の場合
    if node in target_genes:
        if len(target_genes) == 1:
            return color_map['target']
        else:
            idx = target_genes.index(node)
            return color_map[f'target_{idx}']
    
    # ハブノードの場合
    if hubs and node in hubs:
        return color_map['hub_node']
    
    # レベルに応じた色分け
    if node in regulators_by_level['upstream_1']:
        return color_map['upstream_1']
    elif node in regulators_by_level['upstream_2']:
        return color_map['upstream_2']
    elif node in regulators_by_level['downstream_1']:
        return color_map['downstream_1']
    elif node in regulators_by_level['downstream_2']:
        return color_map['downstream_2']
    
    # それ以外の遺伝子
    return color_map['target_genes']


def visualize_static_network(tf_network: nx.DiGraph,
                           regulators_by_level: Dict[str, Set[str]],
                           target_genes: List[str],
                           font_size: int = 12,
                           node_size_factor: float = 1.0,
                           edge_width_factor: float = 1.0,
                           hub_centric: bool = False,
                           hubs: Optional[Set[str]] = None) -> None:
    """静的なネットワーク図を作成（複数ターゲット対応）"""
    plt.clf()
    fig = plt.figure(figsize=(15, 12))
    
    # レイアウトの設定
    if hub_centric and hubs:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
        center_x = np.mean([coord[0] for coord in pos.values()])
        center_y = np.mean([coord[1] for coord in pos.values()])
        
        hub_set = set(hubs) | set(target_genes)
        hub_radius = 0.3
        hub_list = list(hub_set)
        
        for i, hub in enumerate(hub_list):
            angle = 2 * np.pi * i / len(hub_list)
            pos[hub] = np.array([
                center_x + hub_radius * np.cos(angle),
                center_y + hub_radius * np.sin(angle)
            ])
    else:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
    
    # カラーマップの定義（複数ターゲット用）
    color_map = {
        'target': '#ff0000',          # Red (単一ターゲット用)
        'target_0': '#ff0000',        # Red (最初のターゲット用)
        'target_1': '#0000ff',        # Blue (2番目のターゲット用)
        'upstream_1': '#add8e6',      # Light blue
        'upstream_2': '#90ee90',      # Light green
        'downstream_1': '#ffff00',    # Yellow
        'downstream_2': '#ffc0cb',    # Pink
        'target_genes': '#d3d3d3',    # Light gray
        'hub_node': '#FFA500'         # Orange
    }
    
    # ノードの色とサイズの設定
    node_colors = []
    node_sizes = []
    for node in tf_network.nodes():
        size = 1000 * node_size_factor
        color = get_node_color(node, regulators_by_level, color_map, 
                             target_genes, hubs if hub_centric else None)
        
        if node in target_genes or (hub_centric and hubs and node in hubs):
            size *= 2
        elif node in regulators_by_level['upstream_1'] or \
             node in regulators_by_level['downstream_1']:
            size *= 1.5
            
        node_colors.append(color)
        node_sizes.append(size)
    
    # エッジの描画設定
    edge_weights = [tf_network[u][v]['weight'] for u, v in tf_network.edges()]
    max_weight = max(edge_weights) if edge_weights else 1.0
    edge_widths = [(w / max_weight * 3) * edge_width_factor for w in edge_weights]
    
    # エッジの色を設定（ターゲット遺伝子からのエッジは対応する色に）
    edge_colors = []
    for u, v in tf_network.edges():
        if u in target_genes:
            edge_colors.append(color_map[f'target_{target_genes.index(u)}'])
        elif v in target_genes:
            edge_colors.append(color_map[f'target_{target_genes.index(v)}'])
        else:
            edge_colors.append('gray')
    
    # ネットワークの描画
    nx.draw_networkx_nodes(tf_network, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.7)
    nx.draw_networkx_edges(tf_network, pos, width=edge_widths, 
                          edge_color=edge_colors, arrows=True, 
                          arrowsize=20, alpha=0.4,
                          node_size=1000,
                          min_source_margin=10,
                          min_target_margin=15)
    nx.draw_networkx_labels(tf_network, pos, font_size=font_size, 
                          font_weight='bold')
    
    # 凡例の作成
    legend_elements = []
    
    # ターゲット遺伝子の凡例
    for i, target in enumerate(target_genes):
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      label=f'Target {i+1} ({target})', 
                      markerfacecolor=color_map[f'target_{i}'], markersize=15)
        )
    
    if hub_centric and hubs:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      label='Hub Nodes', 
                      markerfacecolor=color_map['hub_node'], markersize=15)
        )
    
    legend_elements.extend([
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Upstream (1st)', 
                  markerfacecolor=color_map['upstream_1'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Upstream (2nd)', 
                  markerfacecolor=color_map['upstream_2'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Downstream (1st)', 
                  markerfacecolor=color_map['downstream_1'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Downstream (2nd)', 
                  markerfacecolor=color_map['downstream_2'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Other Genes', 
                  markerfacecolor=color_map['target_genes'], markersize=15)
    ])
    
    plt.legend(handles=legend_elements, loc='center left', 
              bbox_to_anchor=(1, 0.5))
    plt.title(f"Gene Regulatory Network for {' and '.join(target_genes)}")
    plt.axis('off')
    st.pyplot(fig)

def create_interactive_network(tf_network: nx.DiGraph,
                             regulators_by_level: Dict[str, Set[str]],
                             target_genes: List[str],
                             font_size: int = 12,
                             node_size_factor: float = 1.0,
                             edge_width_factor: float = 1.0,
                             height: int = 600,
                             layout_type: str = "static_spring",
                             layout_seed: int = 42,
                             hub_centric: bool = False,
                             hubs: Optional[Set[str]] = None) -> None:
    """Pyvisを使用してインタラクティブなネットワーク図を作成
    
    Args:
        tf_network (nx.DiGraph): 可視化するネットワーク
        regulators_by_level (Dict[str, Set[str]]): 遺伝子の階層レベル
        target_genes (List[str]): ターゲット遺伝子のリスト
        font_size (int, optional): フォントサイズ. Defaults to 12.
        node_size_factor (float, optional): ノードサイズの倍率. Defaults to 1.0.
        edge_width_factor (float, optional): エッジ幅の倍率. Defaults to 1.0.
        height (int, optional): 描画の高さ. Defaults to 600.
        layout_type (str, optional): レイアウトのタイプ. 
            選択肢: "static_spring", "barnes_hut", "hierarchical", "force_atlas2", "repulsion". 
            Defaults to "static_spring".
        layout_seed (int, optional): レイアウトのランダムシード. Defaults to 42.
        hub_centric (bool, optional): ハブ中心のレイアウトを使用するかどうか. Defaults to False.
        hubs (Optional[Set[str]], optional): ハブ遺伝子のセット. Defaults to None.
    """
    height_px = f"{height}px"
    
    net = Network(height=height_px, width="100%", bgcolor="#ffffff", 
                 font_color="black", directed=True)
    
    # レイアウトオプションの定義
    layout_options = {
        "static_spring": {
            "physics": {
                "enabled": False
            },
            "interaction": {
                "dragNodes": True,
                "dragView": True,
                "zoomView": True,
                "navigationButtons": True,
                "hover": True
            },
            "layout": {
                "randomSeed": layout_seed
            }
        },
        "barnes_hut": {
            "physics": {
                "enabled": True,
                "stabilization": {
                    "enabled": True,
                    "iterations": 100,
                    "updateInterval": 50,
                    "fit": True
                },
                "barnesHut": {
                    "gravitationalConstant": -2000,
                    "centralGravity": 0.3,
                    "springLength": 150,
                    "springConstant": 0.04,
                    "damping": 0.09,
                    "avoidOverlap": 1
                },
                "minVelocity": 0.75,
                "maxVelocity": 30
            }
        },
        "hierarchical": {
            "physics": {
                "enabled": True,
                "hierarchicalRepulsion": {
                    "centralGravity": 0.0,
                    "springLength": 100,
                    "springConstant": 0.01,
                    "nodeDistance": 120
                }
            },
            "layout": {
                "hierarchical": {
                    "enabled": True,
                    "direction": "UD",
                    "sortMethod": "directed",
                    "levelSeparation": 150,
                    "nodeSpacing": 200
                }
            }
        },
        "force_atlas2": {
            "physics": {
                "enabled": True,
                "forceAtlas2Based": {
                    "gravitationalConstant": -50,
                    "centralGravity": 0.01,
                    "springLength": 100,
                    "springConstant": 0.08,
                    "damping": 0.4,
                    "avoidOverlap": 1
                }
            }
        },
        "repulsion": {
            "physics": {
                "enabled": True,
                "repulsion": {
                    "nodeDistance": 200,
                    "centralGravity": 0.2,
                    "springLength": 200,
                    "springConstant": 0.05,
                    "damping": 0.09
                }
            }
        }
    }
    
    # 共通のオプションを追加
    base_options = {
        "interaction": {
            "dragNodes": True,
            "dragView": True,
            "zoomView": True,
            "navigationButtons": True,
            "hover": True
        }
    }

    # レイアウトの選択
    if layout_type not in layout_options:
        raise ValueError(f"Invalid layout_type. Choose from {list(layout_options.keys())}")
    
    # 選択されたレイアウトのオプションを適用
    selected_layout = layout_options[layout_type].copy()
    selected_layout.update(base_options)
    net.options = selected_layout

    # カラーマップの定義
    color_map = {
        'target': '#ff0000',          # Red for single target
        'target_0': '#ff0000',        # Red for first target
        'target_1': '#0000ff',        # Blue for second target
        'upstream_1': '#add8e6',      # Light blue
        'upstream_2': '#90ee90',      # Light green
        'downstream_1': '#ffff00',    # Yellow
        'downstream_2': '#ffc0cb',    # Pink
        'target_genes': '#d3d3d3',    # Light gray
        'hub_node': '#FFA500'         # Orange
    }
    
    # 初期レイアウトの計算（静的レイアウトの場合）
    if layout_type == "static_spring":
        pos = nx.spring_layout(tf_network, k=1, iterations=50, seed=layout_seed)
        if hub_centric and hubs:
            center_x = np.mean([coord[0] for coord in pos.values()])
            center_y = np.mean([coord[1] for coord in pos.values()])
            
            hub_set = set(hubs) | set(target_genes)
            for i, hub in enumerate(hub_set):
                angle = 2 * np.pi * i / len(hub_set)
                pos[hub] = np.array([
                    center_x + 0.3 * np.cos(angle),
                    center_y + 0.3 * np.sin(angle)
                ])
    
    # ノードの追加
    for node in tf_network.nodes():
        size = 20 * node_size_factor
        
        # ノードの色とサイズの決定
        if node in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(node)}']
            size *= 2
        elif hub_centric and hubs and node in hubs:
            color = color_map['hub_node']
            size *= 2
        elif node in regulators_by_level['upstream_1']:
            color = color_map['upstream_1']
            size *= 1.5
        elif node in regulators_by_level['upstream_2']:
            color = color_map['upstream_2']
        elif node in regulators_by_level['downstream_1']:
            color = color_map['downstream_1']
            size *= 1.5
        elif node in regulators_by_level['downstream_2']:
            color = color_map['downstream_2']
        else:
            color = color_map['target_genes']
        
        # ノードの追加（レイアウトに応じて位置を設定）
        if layout_type == "static_spring":
            x, y = pos[node]
            x = x * 500
            y = y * 500
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size}, x=x, y=y, physics=False)
        else:
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size})
    
    # エッジの追加
    edge_weights = [d['weight'] for (u, v, d) in tf_network.edges(data=True)]
    max_weight = max(edge_weights) if edge_weights else 1.0
    
    for u, v, d in tf_network.edges(data=True):
        width = 1 + (d['weight'] / max_weight * 5) * edge_width_factor
        
        # エッジの色の決定
        if u in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(u)}']
        elif v in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(v)}']
        else:
            color = 'gray'
            
        net.add_edge(u, v, width=width, color=color, title=f"CSI: {d['weight']:.3f}")

    # JavaScriptの追加（barnes_hutレイアウトの場合）
    if layout_type == "barnes_hut":
        net.html = net.html.replace('</head>',
            '''
            <script>
            window.addEventListener('load', function() {
                let stabilizationTimeout;
                
                function disablePhysics() {
                    if (network) {
                        network.setOptions({
                            physics: {
                                enabled: false,
                                stabilization: {
                                    enabled: false
                                },
                                barnesHut: {
                                    gravitationalConstant: 0,
                                    centralGravity: 0,
                                    springLength: 0,
                                    springConstant: 0,
                                    damping: 1
                                },
                                minVelocity: 0,
                                maxVelocity: 0
                            }
                        });
                    }
                }

                network.on("stabilizationIterationsDone", function() {
                    clearTimeout(stabilizationTimeout);
                    disablePhysics();
                });

                stabilizationTimeout = setTimeout(disablePhysics, 5000);

                network.on("dragStart", disablePhysics);
                network.on("dragEnd", disablePhysics);
                network.on("drag", disablePhysics);
                network.on("zoom", disablePhysics);
                network.on("click", disablePhysics);
            });
            </script>
            </head>
            ''')

    # ネットワークの保存と表示
    with tempfile.NamedTemporaryFile(delete=False, suffix='.html', mode='w', encoding='utf-8') as f:
        net.save_graph(f.name)
        components.html(open(f.name, 'r', encoding='utf-8').read(), height=height, scrolling=True)
    os.unlink(f.name)


def save_interactive_network(tf_network: nx.DiGraph,
                            regulators_by_level: Dict[str, Set[str]],
                            target_gene: str,
                            output_path: str,
                            font_size: int = 12,
                            node_size_factor: float = 1.0,
                            edge_width_factor: float = 1.0,
                            height: int = 600,
                            layout_type: str = "static_spring",
                            layout_seed: int = 42) -> None:
    """インタラクティブなネットワーク図を保存
    
    Args:
        output_path (str): 保存するファイルのパス（.htmlである必要がある）
    """
    height_px = f"{height}px"
    
    net = Network(height=height_px, width="100%", bgcolor="#ffffff", 
                 font_color="black", directed=True)
    
    # レイアウトオプションの定義
    layout_options = {
        "static_spring": {
            "physics": {
                "enabled": False
            },
            "interaction": {
                "dragNodes": True,
                "dragView": True,
                "zoomView": True,
                "navigationButtons": True,
                "hover": True
            },
            "layout": {
                "randomSeed": layout_seed
            }
        },
        "barnes_hut": {
            "physics": {
                "enabled": True,
                "stabilization": {
                    "enabled": True,
                    "iterations": 100,
                    "updateInterval": 50,
                    "fit": True
                },
                "barnesHut": {
                    "gravitationalConstant": -2000,
                    "centralGravity": 0.3,
                    "springLength": 150,
                    "springConstant": 0.04,
                    "damping": 0.09,
                    "avoidOverlap": 1
                },
                "minVelocity": 0.75,
                "maxVelocity": 30
            }
        },
        "hierarchical": {
            "physics": {
                "enabled": True,
                "hierarchicalRepulsion": {
                    "centralGravity": 0.0,
                    "springLength": 100,
                    "springConstant": 0.01,
                    "nodeDistance": 120
                }
            },
            "layout": {
                "hierarchical": {
                    "enabled": True,
                    "direction": "UD",
                    "sortMethod": "directed",
                    "levelSeparation": 150,
                    "nodeSpacing": 200
                }
            }
        },
        "force_atlas2": {
            "physics": {
                "enabled": True,
                "forceAtlas2Based": {
                    "gravitationalConstant": -50,
                    "centralGravity": 0.01,
                    "springLength": 100,
                    "springConstant": 0.08,
                    "damping": 0.4,
                    "avoidOverlap": 1
                }
            }
        },
        "repulsion": {
            "physics": {
                "enabled": True,
                "repulsion": {
                    "nodeDistance": 200,
                    "centralGravity": 0.2,
                    "springLength": 200,
                    "springConstant": 0.05,
                    "damping": 0.09
                }
            }
        }
    }
    
    # 共通のオプションを追加
    base_options = {
        "interaction": {
            "dragNodes": True,
            "dragView": True,
            "zoomView": True,
            "navigationButtons": True,
            "hover": True
        }
    }

    # レイアウトの選択
    if layout_type not in layout_options:
        raise ValueError(f"Invalid layout_type. Choose from {list(layout_options.keys())}")
    
    # 選択されたレイアウトのオプションを適用
    selected_layout = layout_options[layout_type].copy()
    selected_layout.update(base_options)
    net.options = selected_layout

    color_map = {
        'target': '#ff0000',  # Red
        'upstream_1': '#add8e6',  # Light blue
        'upstream_2': '#90ee90',  # Light green
        'downstream_1': '#ffff00',  # Yellow
        'downstream_2': '#ffc0cb',  # Pink
        'target_genes': '#d3d3d3'  # Light gray
    }
    
    # 初期レイアウトの計算（静的レイアウトの場合）
    if layout_type == "static_spring":
        pos = nx.spring_layout(tf_network, k=1, iterations=50, seed=layout_seed)
    
    # ノードの追加
    for node in tf_network.nodes():
        size = 20 * node_size_factor
        color = color_map['target_genes']
        
        if node == target_gene:
            color = color_map['target']
            size *= 2
        elif node in regulators_by_level['upstream_1']:
            color = color_map['upstream_1']
            size *= 1.5
        elif node in regulators_by_level['upstream_2']:
            color = color_map['upstream_2']
        elif node in regulators_by_level['downstream_1']:
            color = color_map['downstream_1']
            size *= 1.5
        elif node in regulators_by_level['downstream_2']:
            color = color_map['downstream_2']
        
        # 静的レイアウトの場合は位置を設定
        if layout_type == "static_spring":
            # スケーリングと位置の調整
            x, y = pos[node]
            x = x * 500  # スケーリング係数
            y = y * 500
            
            net.add_node(node, label=node, color=color, size=size, 
                        font={'size': font_size},
                        x=x, y=y, physics=False)
        else:
            # 他のレイアウトの場合は通常のノード追加
            net.add_node(node, label=node, color=color, size=size, 
                        font={'size': font_size})
    
    # エッジの追加
    edge_weights = [d['weight'] for (u, v, d) in tf_network.edges(data=True)]
    max_weight = max(edge_weights) if edge_weights else 1.0
    
    for u, v, d in tf_network.edges(data=True):
        width = 1 + (d['weight'] / max_weight * 5) * edge_width_factor
        net.add_edge(u, v, width=width, title=f"CSI: {d['weight']:.3f}")

    # 静的レイアウトの場合は追加のJavaScriptを挿入しない
    if layout_type == "barnes_hut":
        net.html = net.html.replace('</head>',
            '''
            <script>
            window.addEventListener('load', function() {
                let stabilizationTimeout;
                
                function disablePhysics() {
                    if (network) {
                        network.setOptions({
                            physics: {
                                enabled: false,
                                stabilization: {
                                    enabled: false
                                },
                                barnesHut: {
                                    gravitationalConstant: 0,
                                    centralGravity: 0,
                                    springLength: 0,
                                    springConstant: 0,
                                    damping: 1
                                },
                                minVelocity: 0,
                                maxVelocity: 0
                            }
                        });
                    }
                }

                // 安定化完了時に物理演算を無効化
                network.on("stabilizationIterationsDone", function() {
                    clearTimeout(stabilizationTimeout);
                    disablePhysics();
                });

                // バックアップとして一定時間後に強制的に無効化
                stabilizationTimeout = setTimeout(disablePhysics, 5000);

                // 各種イベントで物理演算が再開しないようにする
                network.on("dragStart", disablePhysics);
                network.on("dragEnd", disablePhysics);
                network.on("drag", disablePhysics);
                network.on("zoom", disablePhysics);
                network.on("click", disablePhysics);
            });
            </script>
            </head>
            ''')

    # ネットワークの保存（新しいコード）
    net.save_graph(output_path)

def save_static_network(tf_network: nx.DiGraph,
                        regulators_by_level: Dict[str, Set[str]],
                        target_genes: List[str],
                        output_path: str,
                        font_size: int = 12,
                        node_size_factor: float = 1.0,
                        edge_width_factor: float = 1.0,
                        hub_centric: bool = False,
                        hubs: Optional[Set[str]] = None,
                        format: str = 'png') -> None:
    """静的なネットワーク図を保存（複数ターゲット対応）"""
    # サポートされる形式の検証
    if format.lower() not in ['png', 'pdf']:
        raise ValueError("Format must be either 'png' or 'pdf'")
    
    if not output_path.lower().endswith(f'.{format}'):
        output_path = f"{output_path}.{format}"

    plt.clf()
    fig = plt.figure(figsize=(15, 12))
    
    # レイアウトの設定
    if hub_centric and hubs:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
        center_x = np.mean([coord[0] for coord in pos.values()])
        center_y = np.mean([coord[1] for coord in pos.values()])
        
        hub_set = set(hubs) | set(target_genes)
        hub_radius = 0.3
        hub_list = list(hub_set)
        
        for i, hub in enumerate(hub_list):
            angle = 2 * np.pi * i / len(hub_list)
            pos[hub] = np.array([
                center_x + hub_radius * np.cos(angle),
                center_y + hub_radius * np.sin(angle)
            ])
    else:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
    
    # カラーマップ
    color_map = {
        'target': '#ff0000',          # Red for single target
        'target_0': '#ff0000',        # Red for first target
        'target_1': '#0000ff',        # Blue for second target
        'upstream_1': '#add8e6',      # Light blue
        'upstream_2': '#90ee90',      # Light green
        'downstream_1': '#ffff00',    # Yellow
        'downstream_2': '#ffc0cb',    # Pink
        'target_genes': '#d3d3d3',    # Light gray
        'hub_node': '#FFA500'         # Orange
    }
    
    # ノードの色とサイズの設定
    node_colors = []
    node_sizes = []
    for node in tf_network.nodes():
        size = 1000 * node_size_factor
        if node in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(node)}']
            size *= 2
        elif hub_centric and hubs and node in hubs:
            color = color_map['hub_node']
            size *= 2
        elif node in regulators_by_level['upstream_1']:
            color = color_map['upstream_1']
            size *= 1.5
        elif node in regulators_by_level['upstream_2']:
            color = color_map['upstream_2']
        elif node in regulators_by_level['downstream_1']:
            color = color_map['downstream_1']
            size *= 1.5
        elif node in regulators_by_level['downstream_2']:
            color = color_map['downstream_2']
        else:
            color = color_map['target_genes']
            
        node_colors.append(color)
        node_sizes.append(size)
    
    # エッジの設定
    edge_weights = [tf_network[u][v]['weight'] for u, v in tf_network.edges()]
    max_weight = max(edge_weights) if edge_weights else 1.0
    edge_widths = [(w / max_weight * 3) * edge_width_factor for w in edge_weights]
    
    # エッジの色を設定
    edge_colors = []
    for u, v in tf_network.edges():
        if u in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(u)}']
        elif v in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(v)}']
        else:
            color = 'gray'
        edge_colors.append(color)
    
    # ネットワークの描画
    nx.draw_networkx_nodes(tf_network, pos, node_color=node_colors, 
                          node_size=node_sizes, alpha=0.7)
    nx.draw_networkx_edges(tf_network, pos, width=edge_widths, 
                          edge_color=edge_colors, arrows=True, 
                          arrowsize=20, alpha=0.4,
                          node_size=1000,
                          min_source_margin=10,
                          min_target_margin=15)
    nx.draw_networkx_labels(tf_network, pos, font_size=font_size, 
                          font_weight='bold')
    
    # 凡例の作成
    legend_elements = []
    
    # ターゲット遺伝子の凡例
    for i, target in enumerate(target_genes):
        if len(target_genes) == 1:
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          label=f'Target ({target})', 
                          markerfacecolor=color_map['target'], markersize=15)
            )
        else:
            legend_elements.append(
                plt.Line2D([0], [0], marker='o', color='w', 
                          label=f'Target {i+1} ({target})', 
                          markerfacecolor=color_map[f'target_{i}'], markersize=15)
            )
    
    if hub_centric and hubs:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      label='Hub Nodes', 
                      markerfacecolor=color_map['hub_node'], markersize=15)
        )
    
    legend_elements.extend([
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Upstream (1st)', 
                  markerfacecolor=color_map['upstream_1'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Upstream (2nd)', 
                  markerfacecolor=color_map['upstream_2'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Downstream (1st)', 
                  markerfacecolor=color_map['downstream_1'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Downstream (2nd)', 
                  markerfacecolor=color_map['downstream_2'], markersize=15),
        plt.Line2D([0], [0], marker='o', color='w', 
                  label='Other Genes', 
                  markerfacecolor=color_map['target_genes'], markersize=15)
    ])
    
    plt.legend(handles=legend_elements, loc='center left', 
              bbox_to_anchor=(1, 0.5))
    plt.title(f"Gene Regulatory Network for {' and '.join(target_genes)}")
    plt.axis('off')
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', 
                dpi=300 if format == 'png' else 600)
    plt.close(fig)



def save_interactive_network(tf_network: nx.DiGraph,
                            regulators_by_level: Dict[str, Set[str]],
                            target_genes: List[str],
                            output_path: str,
                            font_size: int = 12,
                            node_size_factor: float = 1.0,
                            edge_width_factor: float = 1.0,
                            height: int = 600,
                            layout_type: str = "static_spring",
                            layout_seed: int = 42,
                            hub_centric: bool = False,
                            hubs: Optional[Set[str]] = None) -> None:
    """インタラクティブなネットワーク図を保存（複数ターゲット対応）"""
    net = Network(height=f"{height}px", width="100%", 
                 bgcolor="#ffffff", font_color="black", directed=True)
    
# レイアウトオプションの定義と適用
    layout_options = {
        "static_spring": {
            "physics": {
                "enabled": False
            },
            "interaction": {
                "dragNodes": True,
                "dragView": True,
                "zoomView": True,
                "navigationButtons": True,
                "hover": True
            },
            "layout": {
                "randomSeed": layout_seed
            }
        },
        "barnes_hut": {
            "physics": {
                "enabled": True,
                "stabilization": {
                    "enabled": True,
                    "iterations": 100,
                    "updateInterval": 50,
                    "fit": True
                },
                "barnesHut": {
                    "gravitationalConstant": -2000,
                    "centralGravity": 0.3,
                    "springLength": 150,
                    "springConstant": 0.04,
                    "damping": 0.09,
                    "avoidOverlap": 1
                },
                "minVelocity": 0.75,
                "maxVelocity": 30
            }
        },
        "hierarchical": {
            "physics": {
                "enabled": True,
                "hierarchicalRepulsion": {
                    "centralGravity": 0.0,
                    "springLength": 100,
                    "springConstant": 0.01,
                    "nodeDistance": 120
                }
            },
            "layout": {
                "hierarchical": {
                    "enabled": True,
                    "direction": "UD",
                    "sortMethod": "directed",
                    "levelSeparation": 150,
                    "nodeSpacing": 200
                }
            }
        },
        "force_atlas2": {
            "physics": {
                "enabled": True,
                "forceAtlas2Based": {
                    "gravitationalConstant": -50,
                    "centralGravity": 0.01,
                    "springLength": 100,
                    "springConstant": 0.08,
                    "damping": 0.4,
                    "avoidOverlap": 1
                }
            }
        },
        "repulsion": {
            "physics": {
                "enabled": True,
                "repulsion": {
                    "nodeDistance": 200,
                    "centralGravity": 0.2,
                    "springLength": 200,
                    "springConstant": 0.05,
                    "damping": 0.09
                }
            }
        }
    }
    
    # 共通のオプションを追加
    base_options = {
        "interaction": {
            "dragNodes": True,
            "dragView": True,
            "zoomView": True,
            "navigationButtons": True,
            "hover": True
        }
    }
    
    if layout_type not in layout_options:
        layout_type = "static_spring"
    net.options = layout_options[layout_type]
    
    # カラーマップ
    color_map = {
        'target': '#ff0000',          # Red for single target
        'target_0': '#ff0000',        # Red for first target
        'target_1': '#0000ff',        # Blue for second target
        'upstream_1': '#add8e6',      # Light blue
        'upstream_2': '#90ee90',      # Light green
        'downstream_1': '#ffff00',    # Yellow
        'downstream_2': '#ffc0cb',    # Pink
        'target_genes': '#d3d3d3',    # Light gray
        'hub_node': '#FFA500'         # Orange
    }
    
    # レイアウト計算
    if layout_type == "static_spring":
        pos = nx.spring_layout(tf_network, k=1, iterations=50, seed=layout_seed)
        if hub_centric and hubs:
            center_x = np.mean([coord[0] for coord in pos.values()])
            center_y = np.mean([coord[1] for coord in pos.values()])
            
            hub_set = set(hubs) | set(target_genes)
            for i, hub in enumerate(hub_set):
                angle = 2 * np.pi * i / len(hub_set)
                pos[hub] = np.array([
                    center_x + 0.3 * np.cos(angle),
                    center_y + 0.3 * np.sin(angle)
                ])
    
    # ノードの追加
    for node in tf_network.nodes():
        size = 20 * node_size_factor
        
        if node in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(node)}']
            size *= 2
        elif hub_centric and hubs and node in hubs:
            color = color_map['hub_node']
            size *= 2
        elif node in regulators_by_level['upstream_1']:
            color = color_map['upstream_1']
            size *= 1.5
        elif node in regulators_by_level['upstream_2']:
            color = color_map['upstream_2']
        elif node in regulators_by_level['downstream_1']:
            color = color_map['downstream_1']
            size *= 1.5
        elif node in regulators_by_level['downstream_2']:
            color = color_map['downstream_2']
        else:
            color = color_map['target_genes']
        
        if layout_type == "static_spring":
            x, y = pos[node]
            x = x * 500
            y = y * 500
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size}, x=x, y=y, physics=False)
        else:
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size})
    
    # エッジの追加
    edge_weights = [d['weight'] for (u, v, d) in tf_network.edges(data=True)]
    max_weight = max(edge_weights) if edge_weights else 1.0
    
    for u, v, d in tf_network.edges(data=True):
        width = 1 + (d['weight'] / max_weight * 5) * edge_width_factor
        
        if u in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(u)}']
        elif v in target_genes:
            if len(target_genes) == 1:
                color = color_map['target']
            else:
                color = color_map[f'target_{target_genes.index(v)}']
        else:
            color = 'gray'
            
        net.add_edge(u, v, width=width, color=color, title=f"CSI: {d['weight']:.3f}")

   # barnes_hutレイアウトの場合はJavaScriptを追加
    if layout_type == "barnes_hut":
        net.html = net.html.replace('</head>',
            '''
            <script>
            window.addEventListener('load', function() {
                let stabilizationTimeout;
                
                function disablePhysics() {
                    if (network) {
                        network.setOptions({
                            physics: {
                                enabled: false,
                                stabilization: {
                                    enabled: false
                                },
                                barnesHut: {
                                    gravitationalConstant: 0,
                                    centralGravity: 0,
                                    springLength: 0,
                                    springConstant: 0,
                                    damping: 1
                                },
                                minVelocity: 0,
                                maxVelocity: 0
                            }
                        });
                    }
                }

                network.on("stabilizationIterationsDone", function() {
                    clearTimeout(stabilizationTimeout);
                    disablePhysics();
                });

                stabilizationTimeout = setTimeout(disablePhysics, 5000);

                network.on("dragStart", disablePhysics);
                network.on("dragEnd", disablePhysics);
                network.on("drag", disablePhysics);
                network.on("zoom", disablePhysics);
                network.on("click", disablePhysics);
            });
            </script>
            </head>
            ''')

    net.save_graph(output_path)


def add_save_buttons_to_main(st, filtered_network, regulators_by_level, target_genes: List[str], 
                             font_size, node_size, edge_width, viz_height, layout_type, 
                             static_format, hub_centric=False, hubs=None):
    """Save buttons for network visualizations (複数ターゲット対応)"""
    import tempfile
    import os
    import base64

    st.subheader("Save Visualizations")
    col1, col2 = st.columns(2)
    
    # ファイル名のためのターゲット遺伝子名を結合
    target_genes_str = "_and_".join(target_genes)
    
    with col1:
        st.write("**Static Network**")
        
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=f'.{static_format}') as tmp_file:
                save_static_network(
                    tf_network=filtered_network,
                    regulators_by_level=regulators_by_level,
                    target_genes=target_genes,  # 変更: target_gene -> target_genes
                    output_path=tmp_file.name,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    hub_centric=hub_centric,
                    hubs=hubs,
                    format=static_format
                )
                tmp_file_path = tmp_file.name
            
            with open(tmp_file_path, "rb") as file:
                base64_file = base64.b64encode(file.read()).decode('utf-8')
            
            href = f'<a href="data:application/octet-stream;base64,{base64_file}" download="{target_genes_str}_network.{static_format}">Download Static Network ({static_format.upper()})</a>'
            st.markdown(href, unsafe_allow_html=True)
        
        except Exception as e:
            st.error(f"Error creating static network file: {e}")
        finally:
            if 'tmp_file_path' in locals() and os.path.exists(tmp_file_path):
                os.unlink(tmp_file_path)
    
    with col2:
        st.write("**Interactive Network**")
        
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix='.html') as tmp_file:
                save_interactive_network(
                    tf_network=filtered_network,
                    regulators_by_level=regulators_by_level,
                    target_genes=target_genes,  # 変更: target_gene -> target_genes
                    output_path=tmp_file.name,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    height=viz_height,
                    layout_type=layout_type,
                    hub_centric=hub_centric,
                    hubs=hubs
                )
                tmp_file_path = tmp_file.name
            
            with open(tmp_file_path, "rb") as file:
                base64_file = base64.b64encode(file.read()).decode('utf-8')
            
            href = f'<a href="data:text/html;base64,{base64_file}" download="{target_genes_str}_network.html">Download Interactive Network (HTML)</a>'
            st.markdown(href, unsafe_allow_html=True)
        
        except Exception as e:
            st.error(f"Error creating interactive network file: {e}")
        finally:
            if 'tmp_file_path' in locals() and os.path.exists(tmp_file_path):
                os.unlink(tmp_file_path)

# ハブの計算部分を修正
def get_hub_nodes(method: str, network: nx.DiGraph, 
                  centrality_data: Dict[str, Dict], 
                  n: int = 5, percentile: float = 0.9) -> Set[str]:
    """
    指定された手法でハブノードを特定

    Args:
        method: ハブ特定手法
        network: ネットワーク
        centrality_data: 中心性指標のデータ
            (weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector)
        n: トップN（中心性ベースの場合）
        percentile: パーセンタイル閾値（従来のハブ手法の場合）
    """
    weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector = centrality_data
    
    if method == "top_degree":
        return set(sorted(weighted_degree.items(), key=lambda x: x[1], reverse=True)[:n])
    elif method == "top_pagerank":
        return set(sorted(weighted_pagerank.items(), key=lambda x: x[1], reverse=True)[:n])
    elif method == "top_betweenness":
        return set(sorted(weighted_between.items(), key=lambda x: x[1], reverse=True)[:n])
    elif method == "top_eigenvector":
        return set(sorted(weighted_eigenvector.items(), key=lambda x: x[1], reverse=True)[:n])
    elif method == "degree_hub":
        return identify_hubs_weighted(network, weight_threshold=percentile/100)
    elif method == "bottleneck_hub":
        return identify_bottlenecks_weighted(network, betweenness_threshold=percentile/100)
    elif method == "common_hub":
        hubs = identify_hubs_weighted(network, weight_threshold=percentile/100)
        bottlenecks = identify_bottlenecks_weighted(network, betweenness_threshold=percentile/100)
        return hubs & bottlenecks
    
    return set()


def main():
    st.title("SCENIC TF Network Analysis")
    
    # メインエリアでのファイルアップロード
    st.header("Data Upload")
    regulon_file = st.file_uploader("Upload Regulon Data (e.g., regulons.seurat_object.txt)", type=['txt', 'tsv'])
    csi_file = st.file_uploader("Upload CSI Data (e.g., csi_results.csv)", type=['csv'])
    
    if regulon_file and csi_file:
        st.header("Analysis Settings")
        
        # データの読み込み
        regulon_data, all_genes = load_regulon_data(regulon_file)
        csi_data = load_csi_data(csi_file)
        all_tfs = load_tfs_list(regulon_file)
        species = determine_species(regulon_file)

        with st.form("input_groups and batch"):     
            # ターゲット遺伝子選択
            target_genes = st.multiselect(
                "Select Target Genes (max 2)",
                sorted(list(all_genes)),
                max_selections=2,
                help="Select up to two target genes for analysis"
            )
            
            # 基本設定
            col1, col2, col3 = st.columns(3)
            with col1:
                max_upstream = st.number_input("Max Upstream Level", 
                                             min_value=1, max_value=3, value=2)
            with col2:
                max_downstream = st.number_input("Max Downstream Level", 
                                               min_value=1, max_value=3, value=2)
            
            # ネットワークモードの設定
            network_mode = st.checkbox("Enable Network Mode", value=False,
                help="If enabled, search both upstream and downstream connections for each TF")
            
            # フィルタリング設定
            st.subheader("Filtering Settings")
            with st.expander("CSI Filtering", expanded=True):
                filtering_method = st.radio("Filtering Method",
                    options=["CSI Threshold", "Max TFs"],
                    help="Choose how to filter TFs")
                
                if filtering_method == "CSI Threshold":
                    csi_threshold = st.slider("CSI Threshold", 0.0, 1.0, 0.0)
                    max_tfs = None
                else:
                    max_tfs = st.slider("Max TFs per level", 1, 50, 20)
                    csi_threshold = None

            # TRRUST設定
            trrust_options = add_trrust_options(st, species)

            # ChIP-Atlas設定
            chip_atlas_options = add_chip_atlas_options(st, species)

            # STRING-db設定
            string_db_options = add_string_db_options(st, species)

            # 中心性とハブの設定
            st.subheader("Hub Selection for Visualization")
            with st.expander("Hub Settings", expanded=True):
                hub_methods = {
                    "Top Weighted Degree Centrality": "top_degree",
                    "Top Weighted PageRank": "top_pagerank",
                    "Top Weighted Betweenness Centrality": "top_betweenness",
                    "Top Weighted Eigenvector Centrality": "top_eigenvector",
                    "Weighted Degree-based Hubs": "degree_hub",
                    "Weighted Bottleneck Nodes": "bottleneck_hub",
                    "Common Weighted Hubs": "common_hub"
                }
                
                col1, col2 = st.columns(2)
                with col1:
                    selected_hub_method = st.selectbox(
                        "Select hub method for visualization",
                        options=list(hub_methods.keys()),
                        index=list(hub_methods.keys()).index("Common Weighted Hubs"),
                        help="Choose method to identify hub nodes for visualization"
                    )
                    
                    top_n = st.number_input(
                        "Number of top nodes (for centrality-based methods)",
                        min_value=1,
                        max_value=20,
                        value=5
                    )
                
                with col2:
                    percentile_threshold = st.slider(
                        "Percentile threshold (for hub-based methods)",
                        min_value=80,
                        max_value=99,
                        value=90,
                        help="Nodes above this percentile will be considered as hubs"
                    )
                    
                    hub_centric = st.checkbox("Center hub nodes in visualization?", 
                                            value=True,
                                            help="If enabled, hub nodes will be placed at the center of the network")
            
            # 可視化設定
            st.subheader("Visualization Settings")
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                font_size = st.slider("Font Size", 10, 30, 12)
            with col2:
                node_size = st.slider("Node Size", 0.5, 2.0, 1.0)
            with col3:
                edge_width = st.slider("Edge Width", 0.2, 2.0, 1.0)
            with col4:
                viz_height = st.slider("Graph Height", 400, 1000, 800)
            
            layout_type = st.selectbox(
                "Interactive network layout type",
                ["static_spring", "barnes_hut", "hierarchical", "force_atlas2", "repulsion"],
                help="Choose the network layout algorithm"
            )
            
            static_format = st.selectbox("Static network format to download", ["pdf", 'png'])

            submitted = st.form_submit_button("Submit")
        
        st.markdown("#### To change options, you need to press submit everytime.")
        
        if len(target_genes) == 0 and submitted:
            st.warning("Please select at least one target gene.")
            return
        
        # 解析開始
        if st.button("Start Analysis", type="primary"):
            with st.spinner('Analyzing network...'):
                # ネットワーク構築
                tf_network, regulators_by_level = build_expanded_tf_network(
                    regulon_data, csi_data, all_tfs,
                    target_genes, max_upstream, max_downstream,
                    network_mode, csi_threshold, max_tfs
                )
    
                # 各種フィルタリングの適用
                filtered_network = tf_network.copy()
                
                if trrust_options["enable"]:
                    trrust_data = load_trrust_data(species)
                    filtered_network, removed_edges = filter_network_by_trrust(
                        filtered_network, trrust_data,
                        trrust_options["interaction_types"],
                        trrust_options["min_pmid_count"]
                    )
                    show_trrust_stats(st, tf_network, filtered_network, removed_edges, trrust_data)
                
                if chip_atlas_options["enable"]:
                    network_tfs = [node for node in filtered_network.nodes() if node in regulon_data]
                    chip_data = get_network_chip_data(network_tfs, species, chip_atlas_options["distance"])
                    filtered_network, removed_edges = filter_network_by_chip_atlas(
                        filtered_network, chip_data,
                        chip_atlas_options["min_score"]
                    )
                    show_chip_atlas_stats(st, tf_network, filtered_network, removed_edges, chip_data)
                
                if string_db_options["enable"]:
                    string_data = get_string_interactions(
                        list(filtered_network.nodes()),
                        string_db_options["species_id"],
                        string_db_options["score_threshold"]
                    )
                    if string_data:
                        filtered_network, removed_edges = filter_network_by_string(
                            filtered_network, string_data,
                            string_db_options["score_threshold"]
                        )
                        show_string_db_stats(st, tf_network, filtered_network, removed_edges, string_data)
                
                # ターゲット遺伝子に接続されていないサブネットワークを削除
                filtered_network, disconnected_edges = cleanup_network(filtered_network, target_genes)
                
                # ネットワーク解析
                if filtered_network.number_of_nodes() >= 2:
                    st.subheader("Network Analysis Results")
                    
                    # 中心性の計算
                    weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector = \
                        calculate_weighted_centrality(filtered_network)
                    
                    # ハブの特定
                    hub_method = hub_methods[selected_hub_method]
                    hubs = get_hub_nodes(
                        method=hub_method,
                        network=filtered_network,
                        centrality_data=(weighted_degree, weighted_between, 
                                       weighted_pagerank, weighted_eigenvector),
                        n=top_n,
                        percentile=percentile_threshold
                    )
                    
                    # 可視化
                    st.subheader("Network Visualization")
                    st.write("Static Network:")
                    visualize_static_network(
                        filtered_network, regulators_by_level, target_genes,
                        font_size, node_size, edge_width,
                        hub_centric=hub_centric, hubs=hubs
                    )
                    
                    st.write("Interactive Network:")
                    create_interactive_network(
                        filtered_network, regulators_by_level, target_genes,
                        font_size, node_size, edge_width, viz_height,
                        layout_type, hub_centric=hub_centric, hubs=hubs
                    )
                    
                    # 保存ボタン
                    add_save_buttons_to_main(
                        st, filtered_network, regulators_by_level, target_genes,
                        font_size, node_size, edge_width, viz_height,
                        layout_type, static_format,
                        hub_centric=hub_centric, hubs=hubs
                    )
                    
                else:
                    st.warning("Not enough nodes for analysis after filtering.")
    else:
        st.info("Please upload both Regulon and CSI data files to start the analysis.")

if __name__ == "__main__":
    main()