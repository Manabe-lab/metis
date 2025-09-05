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
import time
from io import StringIO
import hashlib


def visualize_single_distribution(data: Dict[str, Dict[str, float]], 
                                threshold: float = 0.1,
                                title: str = "Score Distribution"):
    """
    Visualize distribution of a single score type (CSI or Adjacencies)
    
    Args:
        data: Score dictionary (nested dict format)
        threshold: Threshold value to mark on the plot
        title: Title for the plot
    """
    st.subheader(title)
    
    # Collect all scores
    scores = []
    
    # Collect all TFs
    all_tfs = set(data.keys())
    
    for tf1 in all_tfs:
        for tf2 in all_tfs:
            if tf1 == tf2:
                continue
                
            score = data.get(tf1, {}).get(tf2, 0.0)
            
            if score > 0.0:  # Only consider relationships with non-zero score
                scores.append(score)
    
    if not scores:
        st.warning("No score data available for visualization.")
        return
    
    # Basic statistics
    stats_dict = {
        "Count": len(scores),
        "Mean": np.mean(scores),
        "Median": np.median(scores),
        "Min": np.min(scores),
        "Max": np.max(scores),
        "Std Dev": np.std(scores)
    }
    
    col1, col2 = st.columns(2)
    
    with col1:
        for key, value in stats_dict.items():
            if key == "Count":
                st.metric(key, f"{value}")
            else:
                st.metric(key, f"{value:.4f}")
    
    with col2:
        fig, ax = plt.subplots()
        ax.hist(scores, bins=20, edgecolor='black', alpha=0.7)
        ax.axvline(x=threshold, color='red', linestyle='--', label=f'Threshold: {threshold}')
        ax.set_xlabel(f"{title}")
        ax.set_ylabel("Frequency")
        ax.legend()
        st.pyplot(fig)
        plt.close(fig)

# 1. ランクベース変換関数（どのデータソースにも適用可能）
@st.cache_data
def convert_to_rank_based(edge_weights: Dict[str, Dict[str, float]], 
                         use_ties: bool = True) -> Dict[str, Dict[str, float]]:
    """
    エッジウェイトデータをランクベースに変換する
    
    Args:
        edge_weights: 元のエッジウェイトデータ（CSIまたはAdjacencies）
        use_ties: 同じスコアに同じランクを割り当てる場合はTrue、連続ランクの場合はFalse
    
    Returns:
        Dict[str, Dict[str, float]]: ランクベースに変換されたエッジウェイト
    """
    # すべてのエッジとスコアを収集
    all_edges = []
    for source, targets in edge_weights.items():
        for target, score in targets.items():
            if score > 0:  # 正のスコアのみ考慮
                all_edges.append((source, target, score))
    
    if not all_edges:
        return edge_weights
    
    # スコアでソート（降順）
    all_edges.sort(key=lambda x: x[2], reverse=True)
    
    # 結果格納用の辞書
    ranked_weights = defaultdict(dict)
    
    if use_ties:
        # 同じスコアには同じランクを割り当て
        unique_scores = sorted(set(edge[2] for edge in all_edges), reverse=True)
        max_rank = len(unique_scores) - 1
        
        # スコアごとのランクを計算
        score_to_rank = {score: i/max_rank if max_rank > 0 else 0.5 
                       for i, score in enumerate(unique_scores)}
        
        # 各エッジにランクを割り当て
        for source, target, score in all_edges:
            normalized_score = score_to_rank[score]
            ranked_weights[source][target] = normalized_score
    else:
        # 連続ランク
        total_edges = len(all_edges)
        
        for i, (source, target, _) in enumerate(all_edges):
            # ランクを0-1範囲に正規化（最高ランク=1、最低ランク=0）
            normalized_score = 1.0 - (i / (total_edges - 1) if total_edges > 1 else 0)
            ranked_weights[source][target] = normalized_score
    
    return ranked_weights

# 1. adjacenciesデータロード関数
@st.cache_data
def load_adjacencies_data(file) -> Dict[str, Dict[str, float]]:
    """
    Load regulatory relationships and importance scores from SCENIC adjacencies.tsv
    
    Returns:
        Dict[str, Dict[str, float]]: Dictionary with TF as key and another dictionary 
                                     with target genes and importance scores as values
    """


    # ファイルの内容を一度だけ読み込み
    file_content = file.getvalue()
    
    # 安定したハッシュを計算（デバッグ用にログ出力）
    content_hash = hashlib.md5(file_content).hexdigest()

    adjacencies_data = defaultdict(dict)
    content = file_content.decode('utf-8').splitlines()
    
    # Get header line to identify column indices
    header = content[0].split('\t')
    try:
        tf_idx = header.index('TF')
        target_idx = header.index('target')
        importance_idx = header.index('importance')
    except ValueError:
        # Fall back to default column positions if headers don't match
        tf_idx, target_idx, importance_idx = 0, 1, 7  # Common positions in SCENIC adjacencies files
        st.warning("Could not find expected column headers in adjacencies file. Using default positions.")
    
    for line in content[1:]:  # Skip header
        row = line.split('\t')
        if len(row) > max(tf_idx, target_idx, importance_idx):
            tf = row[tf_idx].strip()
            target = row[target_idx].strip()
            try:
                importance = float(row[importance_idx])
                adjacencies_data[tf][target] = importance
            except (ValueError, IndexError):
                continue
    
    return adjacencies_data

# 2. importanceスコア正規化関数
@st.cache_data
def normalize_adjacency_importance(adjacencies_data: Dict[str, Dict[str, float]]) -> Dict[str, Dict[str, float]]:
    """
    Normalize adjacencies importance scores:
    1. Apply log10 transformation to handle extreme ranges
    2. Scale to 0-1 range for compatibility with thresholding
    """
    # Collect all importance scores
    all_scores = []
    for tf, targets in adjacencies_data.items():
        all_scores.extend(targets.values())
    
    if not all_scores:
        return adjacencies_data
    
    # Find minimum non-zero value
    positive_scores = [s for s in all_scores if s > 0]
    if not positive_scores:
        return adjacencies_data
    
    min_positive = min(positive_scores)
    offset = min_positive / 10  # Small offset to handle zeros
    
    # Step 1: Apply log10 transformation
    log_transformed = defaultdict(dict)
    for tf, targets in adjacencies_data.items():
        for target, score in targets.items():
            log_transformed[tf][target] = np.log10(score + offset) if score > 0 else np.log10(offset)
    
    # Collect log-transformed values
    log_values = []
    for tf, targets in log_transformed.items():
        log_values.extend(targets.values())
    
    # Step 2: Min-max normalization to 0-1 range
    min_log = min(log_values)
    max_log = max(log_values)
    range_log = max_log - min_log
    
    # Handle case where all values are identical
    if range_log == 0:
        return adjacencies_data
    
    # Normalize to 0-1 range
    normalized_data = defaultdict(dict)
    for tf, targets in log_transformed.items():
        for target, log_score in targets.items():
            normalized_score = (log_score - min_log) / range_log
            normalized_data[tf][target] = normalized_score
    
    return normalized_data


# 両方向のフィルタリングをサポートする関数
def filter_sequentially(
    csi_data: Dict[str, Dict[str, float]],
    importance_data: Dict[str, Dict[str, float]],
    primary_threshold: float = 0.3,
    secondary_threshold: float = 0.1,
    primary_is_csi: bool = True
) -> Dict[str, Dict[str, float]]:
    """
    Two-step sequential filtering with configurable order
    
    Args:
        csi_data: CSI scores dictionary
        importance_data: Normalized importance scores dictionary
        primary_threshold: Threshold for primary filter
        secondary_threshold: Threshold for secondary filter
        primary_is_csi: If True, filter by CSI first, then by importance
                        If False, filter by importance first, then by CSI
    
    Returns:
        Dict[str, Dict[str, float]]: Filtered scores (uses values from the secondary filter)
    """
    filtered_data = defaultdict(dict)
    
    # Apply filters in the specified order
    if primary_is_csi:
        primary_data = csi_data
        secondary_data = importance_data
    else:
        primary_data = importance_data
        secondary_data = csi_data
    
    # Iterate through all TFs in both datasets
    all_tfs = set(list(csi_data.keys()) + list(importance_data.keys()))
    
    for tf1 in all_tfs:
        for tf2 in all_tfs:
            if tf1 == tf2:
                continue
                
            # First apply primary filter
            primary_score = primary_data.get(tf1, {}).get(tf2, 0.0)
            if primary_score >= primary_threshold:
                # Then apply secondary filter
                secondary_score = secondary_data.get(tf1, {}).get(tf2, 0.0)
                if secondary_score >= secondary_threshold:
                    # Use secondary score as the final value
                    filtered_data[tf1][tf2] = secondary_score
    
    return filtered_data

# 4. ネットワークフィルタリング関数
def filter_network_by_weights(tf_network: nx.DiGraph, weight_threshold: float = 0.0) -> nx.DiGraph:
    """Edge weightsの閾値に基づいてネットワークをフィルタリング"""
    filtered_network = nx.DiGraph()
    
    for u, v, data in tf_network.edges(data=True):
        if data['weight'] >= weight_threshold:
            filtered_network.add_edge(u, v, **data)
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    return filtered_network

# 5. スコア分布可視化関数
def visualize_score_distributions(csi_data: Dict[str, Dict[str, float]], 
                                importance_data: Dict[str, Dict[str, float]],
                                csi_threshold: float = 0.3,
                                importance_threshold: float = 0.1):
    """
    Visualize the distributions of CSI and normalized importance scores
    
    Args:
        csi_data: CSI scores dictionary
        importance_data: Normalized importance scores dictionary
        csi_threshold: Current CSI threshold for reference
        importance_threshold: Current importance threshold for reference
    """
    st.subheader("CSI and Importance Score Distributions")
    
    # Collect all scores
    csi_scores = []
    importance_scores = []
    
    # Pairs of scores for scatter plot
    paired_scores = []
    
    # Collect all TFs
    all_tfs = set(list(csi_data.keys()) + list(importance_data.keys()))
    
    for tf1 in all_tfs:
        for tf2 in all_tfs:
            if tf1 == tf2:
                continue
                
            csi_score = csi_data.get(tf1, {}).get(tf2, 0.0)
            imp_score = importance_data.get(tf1, {}).get(tf2, 0.0)
            
            if csi_score > 0.0:  # Only consider relationships with non-zero CSI
                csi_scores.append(csi_score)
                
            if imp_score > 0.0:  # Only consider relationships with non-zero importance
                importance_scores.append(imp_score)
                
            # If both scores exist, add to paired data
            if csi_score > 0.0 and imp_score > 0.0:
                paired_scores.append((csi_score, imp_score))
    
    # Skip if no scores found
    if not csi_scores and not importance_scores:
        st.warning("No score data available for visualization.")
        return
    
    # Basic statistics
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("### CSI Score Statistics")
        if csi_scores:
            stats_dict = {
                "Count": len(csi_scores),
                "Mean": np.mean(csi_scores),
                "Median": np.median(csi_scores),
                "Min": np.min(csi_scores),
                "Max": np.max(csi_scores),
                "Std Dev": np.std(csi_scores)
            }
            
            for key, value in stats_dict.items():
                if key == "Count":
                    st.metric(key, f"{value}")
                else:
                    st.metric(key, f"{value:.4f}")
        else:
            st.write("No CSI scores available.")
    
    with col2:
        st.write("### Normalized Importance Score Statistics")
        if importance_scores:
            stats_dict = {
                "Count": len(importance_scores),
                "Mean": np.mean(importance_scores),
                "Median": np.median(importance_scores),
                "Min": np.min(importance_scores),
                "Max": np.max(importance_scores),
                "Std Dev": np.std(importance_scores)
            }
            
            for key, value in stats_dict.items():
                if key == "Count":
                    st.metric(key, f"{value}")
                else:
                    st.metric(key, f"{value:.4f}")
        else:
            st.write("No importance scores available.")
    
    # Distribution plots
    st.write("### Score Distributions")
    col1, col2 = st.columns(2)
    
    with col1:
        if csi_scores:
            fig, ax = plt.subplots()
            ax.hist(csi_scores, bins=20, edgecolor='black', alpha=0.7)
            ax.axvline(x=csi_threshold, color='red', linestyle='--', label=f'Threshold: {csi_threshold}')
            ax.set_xlabel("CSI Score")
            ax.set_ylabel("Frequency")
            ax.set_title("Distribution of CSI Scores")
            ax.legend()
            st.pyplot(fig)
            plt.close(fig)
    
    with col2:
        if importance_scores:
            fig, ax = plt.subplots()
            ax.hist(importance_scores, bins=20, edgecolor='black', alpha=0.7)
            ax.axvline(x=importance_threshold, color='red', linestyle='--', label=f'Threshold: {importance_threshold}')
            ax.set_xlabel("Normalized Importance Score")
            ax.set_ylabel("Frequency")
            ax.set_title("Distribution of Normalized Importance Scores")
            ax.legend()
            st.pyplot(fig)
            plt.close(fig)
    
    # Scatter plot of paired scores
    if paired_scores:
        st.write("### Relationship between CSI and Importance Scores")
        
        # Convert to numpy array for easier handling
        paired_array = np.array(paired_scores)
        
        fig, ax = plt.subplots(figsize=(10, 6))
        scatter = ax.scatter(paired_array[:, 0], paired_array[:, 1], alpha=0.3, s=5)
        
        # Draw quadrant lines for thresholds
        ax.axvline(x=csi_threshold, color='red', linestyle='--', alpha=0.7, label=f'CSI Threshold: {csi_threshold}')
        ax.axhline(y=importance_threshold, color='blue', linestyle='--', alpha=0.7, label=f'Importance Threshold: {importance_threshold}')
        
        # Label the quadrants
        ax.text(0.99, 0.99, "High CSI, High Importance", ha='right', va='top', transform=ax.transAxes, fontsize=10)
        ax.text(0.01, 0.99, "Low CSI, High Importance", ha='left', va='top', transform=ax.transAxes, fontsize=10)
        ax.text(0.99, 0.01, "High CSI, Low Importance", ha='right', va='bottom', transform=ax.transAxes, fontsize=10)
        ax.text(0.01, 0.01, "Low CSI, Low Importance", ha='left', va='bottom', transform=ax.transAxes, fontsize=10)
        
        # Calculate counts in each quadrant
        q1 = sum(1 for x, y in paired_scores if x >= csi_threshold and y >= importance_threshold)
        q2 = sum(1 for x, y in paired_scores if x < csi_threshold and y >= importance_threshold)
        q3 = sum(1 for x, y in paired_scores if x >= csi_threshold and y < importance_threshold)
        q4 = sum(1 for x, y in paired_scores if x < csi_threshold and y < importance_threshold)
        
        # Calculate percentages
        total = len(paired_scores)
        q1_pct = q1 / total * 100
        q2_pct = q2 / total * 100
        q3_pct = q3 / total * 100
        q4_pct = q4 / total * 100
        
        # Add count annotations to quadrants
        ax.text(0.99, 0.93, f"{q1} pairs ({q1_pct:.1f}%)", ha='right', transform=ax.transAxes, fontsize=9)
        ax.text(0.01, 0.93, f"{q2} pairs ({q2_pct:.1f}%)", ha='left', transform=ax.transAxes, fontsize=9)
        ax.text(0.99, 0.07, f"{q3} pairs ({q3_pct:.1f}%)", ha='right', transform=ax.transAxes, fontsize=9)
        ax.text(0.01, 0.07, f"{q4} pairs ({q4_pct:.1f}%)", ha='left', transform=ax.transAxes, fontsize=9)
        
        ax.set_xlabel("CSI Score")
        ax.set_ylabel("Normalized Importance Score")
        ax.set_title("CSI vs Normalized Importance Scores")
        ax.legend()
        ax.grid(True, alpha=0.3)
        
        st.pyplot(fig)
        plt.close(fig)
        
        # Calculate correlation
        if len(paired_scores) > 1:
            correlation = np.corrcoef(paired_array[:, 0], paired_array[:, 1])[0, 1]
            st.info(f"Correlation between CSI and Importance scores: {correlation:.4f}")
    
    st.markdown("---")

def extend_network_with_ppi(
    tf_network: nx.DiGraph, 
    target_gene: str, 
    regulon_data: Dict[str, Set[str]], 
    string_data: Dict[str, Dict[str, float]],
    ppi_score_threshold: float = 0.4
) -> Tuple[nx.DiGraph, List[Tuple[str, str]], List[str]]:
    extended_network = tf_network.copy()
    added_edges = []
    ppi_tfs = []

    # Step 1: ネットワークの既存ノードを取得
    existing_nodes = set(extended_network.nodes())
  #  st.write(f"\n### Step 1: Current network nodes")
  #  st.write(f"Number of nodes in network: {len(existing_nodes)}")
    
    # Step 2: target geneのdirect targetのうち、転写因子のみを特定
    target_direct_tfs = set()
    if target_gene in regulon_data:
        target_regulon = regulon_data[target_gene]
        target_direct_tfs = {gene for gene in target_regulon 
                           if gene in regulon_data and gene != target_gene}
    
  #  st.write(f"\n### Step 2: Direct target TFs of {target_gene}")
  #  st.write(f"Number of direct target TFs: {len(target_direct_tfs)}")
  #  st.write("Direct target TFs:", sorted(list(target_direct_tfs)))
    
    # Step 3: 選択基準に基づいてPPIパートナーを選択
    ppi_candidates = []
    for tf in regulon_data.keys():  # 全転写因子をチェック
        # 基準1: 既存ネットワークに含まれていない
        if tf in existing_nodes:
            continue
        
        # 基準2: target geneとPPI相互作用がある
        if (tf not in string_data or 
            target_gene not in string_data[tf] or 
            string_data[tf][target_gene] < ppi_score_threshold):
            continue
        
        # 基準3: target geneのdirect target TFを制御している
        tf_targets = regulon_data[tf]
        shared_tf_targets = target_direct_tfs & tf_targets
        
        if shared_tf_targets:
            ppi_candidates.append({
                'tf': tf,
                'ppi_score': string_data[tf][target_gene],
                'controlled_tfs': shared_tf_targets
            })

    
    # Step 4: ネットワークの拡張
    for candidate in ppi_candidates:
        ppi_tf = candidate['tf']
        ppi_score = candidate['ppi_score']
        controlled_tfs = candidate['controlled_tfs']
        
        # PPIノードを追加
        extended_network.add_node(ppi_tf)
        ppi_tfs.append(ppi_tf)
        
        # 1. target_geneとのPPI相互作用（無向エッジ）
        extended_network.add_edge(target_gene, ppi_tf, weight=ppi_score, type='ppi')
        added_edges.append((target_gene, ppi_tf))

        for target_tf in controlled_tfs:
            if target_tf in extended_network:
                extended_network.add_edge(ppi_tf, target_tf, weight=ppi_score, type='regulation')
                added_edges.append((ppi_tf, target_tf))
    st.markdown("---")
    st.write("### Network Extension Summary")
    # 結果の表示
    if ppi_tfs:
        st.write(f"Added {len(ppi_tfs)} PPI transcription factors")
        st.write(f"Added {len(added_edges)} new edges")
        
        edges_df = pd.DataFrame(added_edges, columns=['Source', 'Target'])
        edges_df['Type'] = ['PPI' if extended_network[s][t]['type'] == 'ppi' else 'Regulation'
                           for s, t in added_edges]
        edges_df['Weight'] = [extended_network[s][t]['weight'] for s, t in added_edges]
        
        st.write("\n**Added Edges Details:**")
        st.dataframe(edges_df)
    else:
        st.markdown("### No PPI TF are found.")
    
    return extended_network, added_edges, ppi_tfs


def cleanup_network(network: nx.DiGraph, target_gene: str) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """
    target geneとの関連性を考慮してネットワークをクリーニング
    
    パスの条件:
    1. target geneへ到達するパス上のノード（上流）
    2. target geneから到達可能なパス上のノード（下流）
        - 直接的な下流
        - 上流制御因子の下流も含む
    
    例：target gene -> A -> B -> C
                          B -> D
    A,B,C,Dは全て保持（DはBを介してtarget geneの制御下にある）
    """
    cleaned_network = network.copy()
    nodes_to_keep = {target_gene}
    
    # 1. target geneへのパス上のノード（上流）を追加
    for node in network.nodes():
        if node != target_gene:
            try:
                path = nx.shortest_path(network, node, target_gene)
                nodes_to_keep.update(path)
            except nx.NetworkXNoPath:
                continue
    
    # 2. target geneからの制御下にあるノードを追加
    upstream_nodes = nodes_to_keep.copy()  # target geneへのパスを持つノード
    
    # 各上流ノードの下流をすべて追加
    for node in upstream_nodes:
        try:
            descendants = nx.descendants(network, node)
            nodes_to_keep.update(descendants)
        except nx.NetworkXError:
            continue
    
    # 削除するノードを特定
    nodes_to_remove = set(network.nodes()) - nodes_to_keep
    
    # 除外されるエッジを記録
    removed_edges = []
    for u, v in network.edges():
        if u in nodes_to_remove or v in nodes_to_remove:
            removed_edges.append((u, v))
    
    # ノードを削除
    cleaned_network.remove_nodes_from(nodes_to_remove)
    
    if removed_edges:
        st.markdown("### Removing disconnected subnetworks:")
        st.write(f"Nodes removed: {len(nodes_to_remove)}")
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

def get_string_interactions(
    proteins: List[str], 
    species: Optional[str] = None,
    score_threshold: float = 0.4
) -> Dict[str, Dict[str, float]]:
    """STRING-dbからPPIデータを取得する（大文字小文字を無視した完全一致）"""
    # 種の判定（渡されていない場合は大文字/小文字パターンから判定）
    if not species:
        is_mostly_uppercase = all(p.isupper() for p in proteins)
        species = "9606" if is_mostly_uppercase else "10090"  # ヒト/マウス
    
    # 種名の変換
    species_map = {
        "human": "9606", 
        "mouse": "10090",
        "9606": "9606", 
        "10090": "10090"
    }
    species_id = species_map.get(species, "9606")  # デフォルトはヒト
    
    # 種の情報を表示
    st.markdown("#### Searching TFs that interacts with the target gene and controls the targets of the target gene.")
    species_name = "Human" if species_id == "9606" else "Mouse"
    st.write(f"Fetching {species_name} protein interaction data from STRING-db")
    
    # シンボルを大文字小文字を無視してノーマライズ
    normalized_proteins = [p.strip().upper() for p in proteins]
 #   st.write(f"Protein names: {normalized_proteins}")
    
    base_url = "https://string-db.org/api"
    output_format = "tsv-no-header"
    method = "network"
    
    params = {
        "identifiers": "%0d".join(normalized_proteins),
        "species": species_id,
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
                    except (ValueError, IndexError):
                        continue
        
        # デバッグ情報
    #    st.write(f"Processed {processed_lines} interactions from STRING-db")
    #    st.write("Raw interactions:")
    #    for interaction in all_interactions[:10]:  # 最初の10個を表示
    #        st.write(f"{interaction[0]} - {interaction[1]}: {interaction[2]}")
        
        # 相互作用の追加
        for protein1, protein2, score in all_interactions:
            # 入力タンパク質リストとの大文字小文字を無視した完全一致チェック
            if (protein1.upper() in normalized_proteins and 
                protein2.upper() in normalized_proteins):
                
                # オリジナルの大文字小文字を保持
                orig_p1 = next(p for p in proteins if p.upper() == protein1.upper())
                orig_p2 = next(p for p in proteins if p.upper() == protein2.upper())
                
                interactions[orig_p1][orig_p2] = score
                interactions[orig_p2][orig_p1] = score  # 双方向の相互作用を追加
                found_interactions += 1
        
        st.write(f"Found {found_interactions} relevant interactions between input proteins")
        
        if not interactions:
            st.warning(f"No {species_name} protein interactions found in STRING-db for the input proteins.")
        else:
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
    
    with st.expander("STRING-db Settings", expanded=False):
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
    """STRING-dbデータに基づいてネットワークをフィルタリング (PPIエッジを除く)"""
    filtered_network = nx.DiGraph()
    filtered_network.add_nodes_from(tf_network.nodes())
    
    # ネットワーク情報
    st.write(f"Original network:")
    st.write(f"Nodes: {tf_network.number_of_nodes()}")
    st.write(f"Edges: {tf_network.number_of_edges()}")
    
    # 元のエッジと相互作用情報を表示
    original_edges = list(tf_network.edges(data=True))
    kept_edges = []
    removed_edges = []
    
    for source, target, data in original_edges:
        # PPIエッジはスキップしてそのまま保持
        if data.get('type') == 'ppi':
            filtered_network.add_edge(source, target, **data)
            kept_edges.append((source, target))
            continue
            
        # STRING-dbでの相互作用を確認
        score = None
        if source in string_data and target in string_data[source]:
            score = string_data[source][target]
        elif target in string_data and source in string_data[target]:
            score = string_data[target][source]
        
        if score is not None and score >= score_threshold:
            filtered_network.add_edge(source, target, **data)
            kept_edges.append((source, target))
        else:
            removed_edges.append((source, target))
    
    # 孤立したノードを削除
    isolated_nodes = list(nx.isolates(filtered_network))
    filtered_network.remove_nodes_from(isolated_nodes)
    
    # 結果の要約
    st.write("\nFiltering Results:")
    st.write(f"Edges in filtered network: {filtered_network.number_of_edges()}")
    st.write(f"Edges removed: {len(removed_edges)}")
    st.write(f"Nodes in filtered network: {filtered_network.number_of_nodes()}")
    st.write(f"Nodes removed: {len(isolated_nodes)}")
        
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
    
    with st.expander("TRRUST Settings", expanded=False):
        # TRRUSTフィルタリングの有効化
        enable_trrust = st.checkbox(
            "Enable TRRUST Filtering", 
            help="Filter network edges based on TRRUST interaction data"
        )
        
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
    """TRRUSTデータに基づいてネットワークをフィルタリング (PPIエッジを除く)"""
    filtered_network = tf_network.copy()
    removed_edges = []
    trrust_edge_details = []  # TRRUSTデータのエッジ詳細を保存するリスト
    
    for source, target in list(tf_network.edges()):
        # PPIエッジはスキップ
        if tf_network[source][target].get('type') == 'ppi':
            continue
            
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
    """ChIP-Atlas データを取得。ローカルファイルを優先し、なければAPIからダウンロード"""
    
    genome = "hg38" if species == "human" else "mm10"
    distance_kb = distance // 1000
    
    # ローカルファイルのパスを構築
    local_dir = f"db/ChIP-Atlas/{genome}/{distance_kb}kb"
    local_file = f"{local_dir}/{tf}.{distance_kb}.tsv"
    
    try:
        # まずローカルファイルを確認
        if os.path.exists(local_file):
            try:
                df = pd.read_csv(local_file, sep='\t')
               # st.write(f"Using local ChIP-Atlas data for {tf}")
                
                # カラム名を標準化
                column_mapping = {
                    'Gene ID': 'gene_id',
                    'Gene Symbol': 'target_gene',
                    'Distance': 'distance',
                    'Average': 'binding_score'
                }
                df = df.rename(columns=column_mapping)
                
                return df
            
            except pd.errors.EmptyDataError:
                st.warning(f"Local file for {tf} is empty, trying ChIP-Atlas API")
            except Exception as e:
                st.warning(f"Error reading local file for {tf}: {str(e)}, trying ChIP-Atlas API")
        
        # ローカルファイルがない場合はAPIから取得
        url = f"https://chip-atlas.dbcls.jp/data/{genome}/target/{tf}.{distance_kb}.tsv"
        response = requests.get(url)
        response.raise_for_status()
        
        if response.text:
            # ディレクトリが存在しない場合は作成
            os.makedirs(local_dir, exist_ok=True)
            
            # データを保存
            with open(local_file, 'w', encoding='utf-8') as f:
                f.write(response.text)
                
            # データをDataFrameに変換
            df = pd.read_csv(StringIO(response.text), sep='\t')
            
            # カラム名を標準化
            column_mapping = {
                'Gene ID': 'gene_id',
                'Gene Symbol': 'target_gene',
                'Distance': 'distance',
                'Average': 'binding_score'
            }
            df = df.rename(columns=column_mapping)
            
           # st.write(f"Downloaded and saved ChIP-Atlas data for {tf}")
            return df
        else:
            st.warning(f"No data found for {tf}")
            return pd.DataFrame()
            
    except requests.exceptions.RequestException as e:
        st.warning(f"Error retrieving data for {tf}: {str(e)}")
        return pd.DataFrame(['No data'])
    except pd.errors.EmptyDataError:
        st.warning(f"No data in response for {tf}")
        return pd.DataFrame()
    except Exception as e:
        st.error(f"Unexpected error for {tf}: {str(e)}")
        return pd.DataFrame()


def filter_network_by_chip_atlas(
    tf_network: nx.DiGraph, 
    chip_data: Dict[str, pd.DataFrame], 
    min_score: float = 0
) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """ChIP-Atlasデータに基づいてネットワークをフィルタリング (PPIエッジを除く)"""
    filtered_network = tf_network.copy()
    removed_edges = []
    
    for source, target in list(tf_network.edges()):
        # PPIエッジはスキップ
        if tf_network[source][target].get('type') == 'ppi':
            continue
            
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

    No_TF = []
    
    for i, tf in enumerate(network_tfs):
        progress_text.write(f"Fetching data for {tf} ({i+1}/{len(network_tfs)})")
        
        df = get_chip_atlas_data(tf, species, distance)
        if not df.empty:
            if df.iloc[0,0] == "No data":
                No_TF.append(tf)
            else:
                chip_data[tf] = df
            
        progress_bar.progress((i + 1) / len(network_tfs))
        time.sleep(1)  # APIの負荷を考慮
    
    progress_text.empty()
    progress_bar.empty()
    
    if len(No_TF) > 0:
        st.markdown(f"##### No data available in ChIp-Atlas for { ', '.join(No_TF)}")
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
    
    with st.expander("ChIP-Atlas Settings", expanded=False):
        # ChIP-Atlasフィルタリングの有効化
        enable_chip_atlas = st.checkbox(
            "Enable ChIP-Atlas Filtering", 
            help="Filter network edges based on ChIP-Atlas binding data"
        )
        

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
                value=20,
                step=1,
                help="Filter edges based on ChIP-Atlas binding score. -10*Log10[MACS2 Q-value]. -10*log10(0.05) ≈ 13, q=0.01 -> 20)."
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

    
def build_expanded_tf_network(
    regulon_data: Dict[str, Set[str]], 
    csi_data: Dict[str, Dict[str, float]], 
    all_tfs: Set[str],  # 全TFのリスト
    target_gene: str,
    max_upstream: int,
    max_downstream: int,
    network_mode: bool,
    use_equal_weights: bool = False,  # 新しいパラメータ
    csi_threshold: float = None,
    max_tfs_per_level: int = None
) -> Tuple[nx.DiGraph, Dict[str, Set[str]]]:
    """拡張された双方向のTFネットワークを構築"""
    tf_network = nx.DiGraph()

    regulators_by_level = {
        **{f'upstream_{i}': set() for i in range(1, 4)},
        **{f'downstream_{i}': set() for i in range(1, 4)},
        'additional_tfs': set()
    }
    
    def get_connected_tfs_for_level(current_tf: str, direction: str) -> List[Tuple[str, float, str]]:
        """特定の方向の接続TFを取得"""
        connected_tfs = []
        
        for tf, genes in regulon_data.items():
            # 上流と下流を探索
            if direction in ['upstream', 'both']:
                if current_tf in genes and tf in all_tfs:
                    # CSIデータがない、または等しい重みを使用する場合は1.0を使用
                    weight = 1.0 if use_equal_weights else csi_data[tf].get(current_tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'upstream'))
            
            # 下流の探索
            if direction in ['downstream', 'both']:
                if current_tf in regulon_data and tf in regulon_data[current_tf] and tf in all_tfs:
                    # CSIデータがない、または等しい重みを使用する場合は1.0を使用
                    weight = 1.0 if use_equal_weights else csi_data[current_tf].get(tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'downstream'))
        
        # CSIの値でソート（等しい重みの場合はランダムな順序になる）
        connected_tfs.sort(key=lambda x: x[1], reverse=True)
        if max_tfs_per_level is not None:
            connected_tfs = connected_tfs[:max_tfs_per_level]
        
        return connected_tfs

    # 上流の探索
    current_level_tfs = {target_gene}
    for level in range(1, max_upstream + 1):
        next_level_tfs = set()
        level_key = f'upstream_{level}'
        
        for current_tf in current_level_tfs:
            direction = 'upstream' if not network_mode else 'both'
            connected_tfs = get_connected_tfs_for_level(current_tf, direction)
            
            for tf, weight, conn_direction in connected_tfs:
                # TFの分類と接続の追加
                if tf not in all_tfs:
                    regulators_by_level['additional_tfs'].add(tf)
                else:
                    regulators_by_level[level_key].add(tf)
                
                # エッジの方向を接続の方向に基づいて追加
                if conn_direction == 'upstream':
                    tf_network.add_edge(tf, current_tf, weight=weight)
                else:
                    tf_network.add_edge(current_tf, tf, weight=weight)
                
                next_level_tfs.add(tf)
        
        current_level_tfs = next_level_tfs

    # 下流の探索
    if target_gene in regulon_data and max_downstream > 0:
        current_level_tfs = {target_gene}
        for level in range(1, max_downstream + 1):
            next_level_tfs = set()
            level_key = f'downstream_{level}'
            
            for current_tf in current_level_tfs:
                direction = 'downstream' if not network_mode else 'both'
                connected_tfs = get_connected_tfs_for_level(current_tf, direction)
                
                for tf, weight, conn_direction in connected_tfs:
                    # TFの分類と接続の追加
                    if tf not in all_tfs:
                        regulators_by_level['additional_tfs'].add(tf)
                    else:
                        regulators_by_level[level_key].add(tf)
                    
                    # エッジの方向を接続の方向に基づいて追加
                    if conn_direction == 'downstream':
                        tf_network.add_edge(current_tf, tf, weight=weight)
                    else:
                        tf_network.add_edge(tf, current_tf, weight=weight)
                    
                    next_level_tfs.add(tf)
            
            current_level_tfs = next_level_tfs
    
    return tf_network, regulators_by_level


def calculate_weighted_centrality(tf_network: nx.DiGraph, 
                                  use_equal_weights: bool = False) -> Tuple[Dict, Dict, Dict, Dict]:
    """重み付きネットワークの各種中心性指標を計算"""
    try:
        # ネットワーク内にエッジの重みがない場合は1.0に設定
        if use_equal_weights:
            for u, v, data in tf_network.edges(data=True):
                if 'weight' not in data:
                    data['weight'] = 1.0
        
        # これ以降は既存のコードと同じ
        weighted_degree = {}
        for node in tf_network.nodes():
            in_weight = sum(d['weight'] for u, v, d in tf_network.in_edges(node, data=True))
            out_weight = sum(d['weight'] for u, v, d in tf_network.out_edges(node, data=True))
            weighted_degree[node] = (in_weight + out_weight) / (2 * tf_network.number_of_nodes() - 2)
        
        # 他の中心性計算も同様に修正可能
        betweenness_centrality = nx.betweenness_centrality(tf_network, weight='weight')
        pagerank = nx.pagerank(tf_network, weight='weight')
        eigenvector_centrality = nx.eigenvector_centrality(tf_network, weight='weight')
        
    except Exception as e:
        # フォールバック
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

def get_edge_color(source_node, node_categories, color_scheme, target, ppi_tfs=None, hubs=None):
   """エッジの色を取得する関数
   
   Args:
       source_node (str): エッジの起点となるノード
       node_categories (dict): カテゴリーごとのノードのセット
       color_scheme (dict): カテゴリーごとの色の辞書
       target (str): ターゲット遺伝子の名前
       ppi_tfs (list, optional): PPIトランスクリプション因子のリスト. Defaults to None.
       hubs (set, optional): hubノードのセット. Defaults to None.
   
   Returns:
       str: 16進数のカラーコード
   """    
   # ppi_tfsがNoneの場合は空のリストに
   ppi_tfs = ppi_tfs or []

   # PPIトランスクリプション因子の特別な処理
   if source_node in ppi_tfs:
       return color_scheme['target']
   
   # ターゲット遺伝子からのエッジは常にターゲット色を優先
   if source_node == target:
       return color_scheme['target']
   
   # hub nodeのチェック（ターゲット遺伝子以外）
   if hubs and source_node in hubs and source_node != target:
       return color_scheme['hub_node']
       
   if source_node in node_categories['upstream_1']:
       return color_scheme['upstream_1']
   elif source_node in node_categories['upstream_2']:
       return color_scheme['upstream_2']
   elif source_node in node_categories['downstream_1']:
       return color_scheme['downstream_1']
   elif source_node in node_categories['downstream_2']:
       return color_scheme['downstream_2']
   
   return color_scheme['target_genes']  # デフォルト色


def adjust_node_positions(pos: Dict[str, np.ndarray], min_dist: float = 0.2) -> Dict[str, np.ndarray]:
    """
    ノードの位置を調整して重なりを防ぐ
    
    Args:
        pos: 元のノード位置の辞書
        min_dist: ノード間の最小距離
    
    Returns:
        Dict[str, np.ndarray]: 調整後のノード位置
    """
    adjusted_pos = pos.copy()
    nodes = list(pos.keys())
    
    # 全てのノードペアについて
    for i, node1 in enumerate(nodes):
        for node2 in nodes[i+1:]:
            pos1 = adjusted_pos[node1]
            pos2 = adjusted_pos[node2]
            
            # ノード間の距離を計算
            dist = np.linalg.norm(pos1 - pos2)
            
            if dist < min_dist:
                # 距離が最小距離未満の場合、ノードを反発させる
                direction = (pos1 - pos2) / dist
                adjustment = direction * (min_dist - dist) / 2
                
                adjusted_pos[node1] = pos1 + adjustment
                adjusted_pos[node2] = pos2 - adjustment
    
    return adjusted_pos

def get_network_layout(tf_network: nx.DiGraph, 
                      layout_type: str,
                      hub_centric: bool = False,
                      hubs: Optional[Set[str]] = None,
                      target_gene: str = None,
                      scale: float = 3.0,      # スケール係数
                      k: float = 2.0,          # ノード間の距離係数
                      min_dist: float = 0.2,    # 最小ノード間距離
                      iterations: int = 100     # 反復回数
                      ) -> Dict[str, np.ndarray]:
    """
    改良されたネットワークレイアウトを生成する関数
    
    Args:
        tf_network: ネットワークグラフ
        layout_type: レイアウトタイプ
        hub_centric: ハブ中心のレイアウトにするかどうか
        hubs: ハブノードのセット
        target_gene: ターゲット遺伝子
        scale: レイアウトのスケール係数
        k: ノード間の距離係数
        min_dist: 最小ノード間距離
        iterations: レイアウトアルゴリズムの反復回数
    """
    
    # 基本レイアウトの生成
    if layout_type == "spring":
        pos = nx.spring_layout(tf_network, k=k, scale=scale, iterations=iterations)
    elif layout_type == "kamada_kawai":
        try:
            # スケールとweightパラメータを追加
            pos = nx.kamada_kawai_layout(tf_network, scale=scale, weight='weight')
            # 位置を更に調整してノードの重なりを防ぐ
            pos = adjust_node_positions(pos, min_dist=min_dist)
        except:
            pos = nx.spring_layout(tf_network, k=k, scale=scale)
    elif layout_type == "fruchterman_reingold":
        pos = nx.fruchterman_reingold_layout(tf_network, k=k, scale=scale, iterations=iterations)
    elif layout_type == "spectral":
        try:
            pos = nx.spectral_layout(tf_network, scale=scale)
            # スペクトラルレイアウトの後にFruchterman-Reingoldで微調整
            pos = nx.fruchterman_reingold_layout(tf_network, pos=pos, k=k/2, iterations=50)
            pos = adjust_node_positions(pos, min_dist=min_dist)
        except:
            pos = nx.spring_layout(tf_network, k=k, scale=scale)
    elif layout_type == "shell":
        pos = nx.shell_layout(tf_network, scale=scale)
    else:
        pos = nx.spring_layout(tf_network, k=k, scale=scale)

    # Hub-centricレイアウトの適用
    if hub_centric and hubs:
        # 中心点の計算
        center_x = sum(pos[node][0] for node in tf_network.nodes()) / len(tf_network.nodes())
        center_y = sum(pos[node][1] for node in tf_network.nodes()) / len(tf_network.nodes())
        
        # ハブノードとターゲット遺伝子を中心に配置
        hub_nodes = set(hubs) | {target_gene}
        n_hubs = len(hub_nodes)
        
        # 中央付近に円形に配置
        for i, hub in enumerate(hub_nodes):
            angle = 2 * np.pi * i / n_hubs
            r = 0.3 * scale  # スケールに応じた半径
            pos[hub] = np.array([
                center_x + r * np.cos(angle),
                center_y + r * np.sin(angle)
            ])
        
        # 他のノードを外側に配置
        for node in tf_network.nodes():
            if node not in hub_nodes:
                current_pos = pos[node]
                direction = current_pos - np.array([center_x, center_y])
                if np.linalg.norm(direction) > 0:
                    direction = direction / np.linalg.norm(direction)
                    pos[node] = np.array([center_x, center_y]) + direction * scale
        
        # 最終的な位置調整
        pos = adjust_node_positions(pos, min_dist=min_dist)
    
    return pos


def add_visualization_settings(st):
    """
    Adds visualization settings UI elements to the Streamlit app
    
    Returns:
        tuple: (static_layout, interactive_layout, layout_scale, min_node_dist, hub_centric, color_hubs)
    """
    st.subheader("Visualization Settings")
    
    # Layout settings
    col1, col2 = st.columns(2)
    with col1:
        static_layout = st.selectbox(
            "Static network layout type",
            ["spring", "kamada_kawai", "fruchterman_reingold", "spectral", "shell"],
            help="Choose the layout algorithm for static visualization"
        )
    with col2:    
        interactive_layout = st.selectbox(
            "Interactive network layout type",
            ["static_spring", "barnes_hut", "hierarchical", "force_atlas2", "repulsion"],
            help="Choose the layout algorithm for interactive visualization"
        )
    
    with st.expander("Static Layout Options", expanded=False):
        col1, col2 = st.columns(2)
        with col1:
            layout_scale = st.slider(
                "Layout scale",
                min_value=1.0,
                max_value=5.0,
                value=3.0,
                help="Adjust the overall scale of the network layout"
            )
        with col2:    
            min_node_dist = st.slider(
                "Minimum node distance",
                min_value=0.1,
                max_value=1.0,
                value=0.2,
                help="Minimum distance between nodes to prevent overlap"
            )

    # Hub visualization options
    col1, col2 = st.columns(2)
    with col1:
        hub_centric = st.checkbox(
            "Center hub nodes in visualization?", 
            value=False,
            help="If enabled, hub nodes will be placed at the center of the network"
        )
    with col2:
        color_hubs = st.checkbox(
            "Highlight hub nodes with color?",
            value=True,
            help="If enabled, hub nodes will be colored differently regardless of layout"
        )

    # Cleanup network オプションを追加
    with st.expander("Network Cleanup Options", expanded=False):
        enable_cleanup = st.checkbox(
            "Enable network cleanup",
            value=True,
            help="Remove nodes and edges that are not connected to the target gene"
        )
    
    return static_layout, interactive_layout, layout_scale, min_node_dist, hub_centric, color_hubs, enable_cleanup
    

def visualize_static_network(tf_network: nx.DiGraph,
                           regulators_by_level: Dict[str, Set[str]],
                           target_gene: str,
                           ppi_tfs: Optional[List[str]] = None,
                           font_size: int = 12,
                           node_size_factor: float = 1.0,
                           edge_width_factor: float = 1.0,
                           hub_centric: bool = False,
                           color_hubs: bool = True,
                           hubs: Optional[Set[str]] = None,
                           node_alpha=0.7,
                           edge_alpha=0.6,
                           layout_type: str = "spring",
                           layout_scale: float = 3.0,
                           min_node_dist: float = 0.2) -> None:
    plt.clf()  # ここでのみクリア
    # Add this before both visualization functions
    color_map = {
        'target': '#ff0000',      # Red
        'upstream_1': '#add8e6',  # Light blue
        'upstream_2': '#90ee90',  # Light green
        'downstream_1': '#ffff00', # Yellow
        'downstream_2': '#ffc0cb', # Pink
        'target_genes': '#d3d3d3', # Light gray
        'hub_node': '#FFA500',    # Orange
        'ppi_node': '#800080',    # Purple
        'ppi_edge': '#8B4513',    # Saddle Brown
        'regular_edge': '#add8e6'  # Light blue
    }

    fig = plt.figure(figsize=(15, 12))
    
    # Get layout
    pos = get_network_layout(
        tf_network=tf_network,
        layout_type=layout_type,
        hub_centric=hub_centric,
        hubs=hubs,
        target_gene=target_gene,
        scale=layout_scale,
        min_dist=min_node_dist
    )
    
    # ppi_tfsがNoneの場合の処理
    ppi_tfs = ppi_tfs or []
    
    # ノードの色とサイズの設定
    node_colors = []
    node_sizes = []
    
    for node in tf_network.nodes():
        size = 1000 * node_size_factor
        color = color_map['target_genes']
        
        if node == target_gene:
            color = color_map['target']
            size *= 2
        elif node in ppi_tfs:
            color = color_map['ppi_node']
            size *= 1.5
        elif hubs and node in hubs and (hub_centric or color_hubs):
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
            
        node_colors.append(color)
        node_sizes.append(size)
     
    # ノードの描画
    nx.draw_networkx_nodes(tf_network, pos, 
                          node_color=node_colors,
                          node_size=node_sizes,
                          alpha=node_alpha)
    
    # エッジの分類
    regular_edges = []
    ppi_edges = []
    regulation_edges = []
    
    for u, v, data in tf_network.edges(data=True):
        edge_type = data.get('type', 'normal')
        if edge_type == 'ppi':
            ppi_edges.append((u, v))
        elif edge_type == 'regulation':
            regulation_edges.append((u, v))
        else:
            regular_edges.append((u, v))

    # エッジの色の設定
    edge_colors = [
        get_edge_color(source, regulators_by_level, color_map, target_gene, ppi_tfs, hubs) 
        for source, _ in tf_network.edges()
    ]
    
    # 通常のエッジを描画
    if regular_edges:
        nx.draw_networkx_edges(tf_network, pos,
                             edgelist=regular_edges,
                             width=2.0 * edge_width_factor,
                             edge_color=edge_colors,
                             arrows=True,
                             arrowsize=20,
                             alpha=edge_alpha,
                             min_source_margin=15,
                             min_target_margin=20)
    
    # PPIエッジを描画（無向）
    if ppi_edges:
        nx.draw_networkx_edges(tf_network, pos,
                             edgelist=ppi_edges,
                             width=2.5 * edge_width_factor,
                             edge_color=edge_colors,
                             arrows=False,
                             alpha=edge_alpha,
                             style='dashed',
                             min_source_margin=15,
                             min_target_margin=15)
    
    # 制御エッジを描画（有向）
    if regulation_edges:
        nx.draw_networkx_edges(tf_network, pos,
                             edgelist=regulation_edges,
                             width=2.5 * edge_width_factor,
                             edge_color=edge_colors,
                             arrows=True,
                             arrowsize=20,
                             alpha=edge_alpha,
                             style='solid',
                             min_source_margin=15,
                             min_target_margin=20)
    
    # ラベルの描画
    nx.draw_networkx_labels(tf_network, pos, font_size=font_size, 
                          font_weight='bold')
    
    # 凡例の作成
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', 
                  label=f'Target ({target_gene})', 
                  markerfacecolor=color_map['target'], markersize=15),
    ]
    
    if hub_centric and hubs:
        legend_elements.append(
            plt.Line2D([0], [0], marker='o', color='w', 
                      label='Hub Nodes', 
                      markerfacecolor=color_map['hub_node'], markersize=15)
        )

    if ppi_tfs:
        legend_elements.extend([
            plt.Line2D([0], [0], marker='o', color='w', 
                      label='PPI Transcription Factors', 
                      markerfacecolor=color_map['ppi_node'], markersize=15),
            plt.Line2D([0], [0], color=color_map['ppi_edge'], 
                      label='PPI Interactions', 
                      linewidth=3),
            plt.Line2D([0], [0], color=color_map['ppi_edge'],
                      label='PPI TF Regulations',
                      linewidth=3, marker='>', markersize=10)
        ])
    
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
                  label='Target Genes', 
                  markerfacecolor=color_map['target_genes'], markersize=15)
    ])

    plt.legend(handles=legend_elements, loc='center left', 
              bbox_to_anchor=(1, 0.5))
    plt.title(f"Gene Regulatory Network Centered on {target_gene}")
    plt.axis('off')
    st.pyplot(fig)

def create_interactive_network(tf_network: nx.DiGraph,
                             regulators_by_level: Dict[str, Set[str]],
                             target_gene: str,
                             ppi_tfs: Optional[List[str]] = None,
                             font_size: int = 12,
                             node_size_factor: float = 1.0,
                             edge_width_factor: float = 1.0,
                             height: int = 600,
                             layout_type: str = "static_spring",
                             layout_seed: int = 42,
                             hub_centric: bool = False,
                             color_hubs: bool = True,
                             hubs: Optional[Set[str]] = None) -> None:
 
    height_px = f"{height}px"
    net = Network(height=height_px, width="100%", bgcolor="#ffffff", 
                 font_color="black", directed=True)
                 
    # ppi_tfsがNoneの場合の処理
    ppi_tfs = ppi_tfs or []

    # カラーマップの定義
    color_map = {
        'target': '#ff0000',      # Red
        'upstream_1': '#add8e6',  # Light blue
        'upstream_2': '#90ee90',  # Light green
        'downstream_1': '#ffff00', # Yellow
        'downstream_2': '#ffc0cb', # Pink
        'target_genes': '#d3d3d3', # Light gray
        'hub_node': '#FFA500',    # Orange
        'ppi_node': '#800080',    # Purple
        'ppi_edge': '#8B4513',    # Saddle Brown
        'regular_edge': '#add8e6'  # Light blue
    }
    
    # レイアウトオプションの定義
    layout_options = {
    "static_spring": {"physics": {"enabled": False}},
    "barnes_hut": {
        "physics": {
            "enabled": True,
            "barnesHut": {
                "gravitationalConstant": -500,
                "centralGravity": 0.8,
                "springLength": 80,
                "springConstant": 0.01,
                "damping": 0.5,
                "avoidOverlap": 1
            },
            "maxVelocity": 15,
            "minVelocity": 1,
            "timestep": 0.2,
            "stabilization": {
                "enabled": True,
                "iterations": 200,
                "updateInterval": 10,
                "onlyDynamicEdges": False,
                "fit": True
            }
        }
    },
"hierarchical": {
    "layout": {
        "hierarchical": {
            "enabled": True,
            "direction": "UD",
            "sortMethod": "hubsize",    # hubsizeに変更
            "levelSeparation": 75,      # 適度な縦間隔
            "nodeSpacing": 100,         # 水平間隔は維持
            "treeSpacing": 100,         # 水平間隔は維持
            "blockShifting": True,
            "edgeMinimization": True,
            "parentCentralization": True,
            "improvedLayout": True,
            "levelBias": 0.0,           # バイアスを解除
        }
    },
    "edges": {
        "smooth": {
            "enabled": True,
            "type": "cubicBezier",      # エッジの描画方法を変更
            "roundness": 0.5
        }
    },
    "physics": {
        "enabled": True,
        "hierarchicalRepulsion": {
            "centralGravity": 0.0,      # 中心への引力を無効化
            "springLength": 100,
            "springConstant": 0.01,
            "nodeDistance": 80,
            "damping": 0.3
        },
        "solver": "hierarchicalRepulsion",
        "stabilization": {
            "enabled": True,
            "iterations": 1000,         # 反復回数を増やす
            "updateInterval": 10,
            "fit": True
        }
    }
},
    "force_atlas2": {
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -50,
                "centralGravity": 0.01,
                "springLength": 100,
                "springConstant": 0.08
            }
        }
    },
    "repulsion": {
        "physics": {
            "repulsion": {
                "nodeDistance": 200,
                "centralGravity": 0.2,
                "springLength": 200,
                "springConstant": 0.05
            }
        }
    }
    }

    # Hub-centricレイアウトの計算
    if hub_centric and hubs:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
        center_x = sum(pos[node][0] for node in tf_network.nodes()) / len(tf_network.nodes())
        center_y = sum(pos[node][1] for node in tf_network.nodes()) / len(tf_network.nodes())
        
        # ハブ遺伝子とターゲット遺伝子を中央に配置
        hub_nodes = set(hubs) | {target_gene}
        n_hubs = len(hub_nodes)
        
        # 中央付近に円形に配置
        for i, hub in enumerate(hub_nodes):
            angle = 2 * np.pi * i / n_hubs
            r = 0.3  # 中心からの距離
            pos[hub] = np.array([
                center_x + r * np.cos(angle),
                center_y + r * np.sin(angle)
            ])
        
        # 他のノードは外側に配置
        for node in tf_network.nodes():
            if node not in hub_nodes:
                current_pos = pos[node]
                # 中心からの方向ベクトルを計算
                direction = current_pos - np.array([center_x, center_y])
                # 正規化して一定距離に配置
                if np.linalg.norm(direction) > 0:
                    direction = direction / np.linalg.norm(direction)
                    pos[node] = np.array([center_x, center_y]) + direction * 1.0
    else:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)

    # レイアウトの選択と適用
    if layout_type in layout_options:
        net.options = layout_options[layout_type]


    # ノードの追加
    for node in tf_network.nodes():
        size = 20 * node_size_factor
        color = color_map['target_genes']
        
        if node == target_gene:
            color = color_map['target']
            size *= 2
        elif node in ppi_tfs:
            color = color_map['ppi_node']
            size *= 1.5
        elif hubs and node in hubs and (hub_centric or color_hubs):
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
                    
        if layout_type == "static_spring":
            x, y = pos[node]
            x = x * 500
            y = y * 500
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size}, x=x, y=y, physics=False)
        else:
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size})

    # エッジの分類
    regular_edges = []
    ppi_edges = []
    regulation_edges = []
    
    for u, v, data in tf_network.edges(data=True):
        edge_type = data.get('type', 'normal')
        if edge_type == 'ppi':
            ppi_edges.append((u, v))
        elif edge_type == 'regulation':
            regulation_edges.append((u, v))
        else:
            regular_edges.append((u, v))

    # 1. 通常のエッジ
    for u, v in regular_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.0 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_gene, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}},
                    physics=True)

    # 2. PPIエッジ（無向）
    for u, v in ppi_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        net.add_edge(u, v, width=width, color=color_map['ppi_edge'],
                    arrows={'to': {'enabled': False}},
                    physics=True)

    # 3. 制御エッジ（有向）
    for u, v in regulation_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_gene, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}},
                    physics=True)

    # JavaScriptによる制御を修正
    javascript_code = '''
<script>
window.addEventListener('load', function() {
    let stabilizationTimeout;
    let isStabilized = false;
    let forceStopTimeout;

    function disablePhysics() {
        if (!isStabilized && typeof network !== 'undefined') {
            network.setOptions({physics: {enabled: false}});
            network.stabilize(0);  // 強制的に安定化
            isStabilized = true;
        }
    }

    // 短いタイムアウト (5秒)
    stabilizationTimeout = setTimeout(disablePhysics, 5000);

    // 絶対的なタイムアウト (8秒)
    forceStopTimeout = setTimeout(function() {
        if (!isStabilized) {
            console.log("Force stopping physics after absolute timeout");
            disablePhysics();
        }
    }, 8000);

    if (typeof network !== 'undefined') {
        network.on("stabilizationIterationsDone", function() {
            clearTimeout(stabilizationTimeout);
            clearTimeout(forceStopTimeout);
            disablePhysics();
        });

        network.on("stabilizationProgress", function(params) {
            let progress = Math.round((params.iterations / params.total) * 100);
            if (progress >= 80) {  // 80%で停止
                clearTimeout(stabilizationTimeout);
                clearTimeout(forceStopTimeout);
                disablePhysics();
            }
        });

        // 任意のユーザーインタラクションで停止
        ['click', 'dragStart', 'zoom'].forEach(function(event) {
            network.on(event, function() {
                if (!isStabilized) {
                    clearTimeout(stabilizationTimeout);
                    clearTimeout(forceStopTimeout);
                    disablePhysics();
                }
            });
        });
    }
});
</script>
    '''

    # 物理演算の制御（barnes_hutの場合）
    if layout_type == "barnes_hut":
        net.html = net.html.replace('</head>',
            '''
            <script>
            window.addEventListener('load', function() {
                // 30秒後に強制的に物理シミュレーションを停止
                setTimeout(function() {
                    if (network) {
                        network.setOptions({physics: {enabled: false}});
                        console.log("Physics simulation stopped by timeout");
                    }
                }, 30000);

                // 安定化後の処理
                network.on("stabilizationIterationsDone", function() {
                    console.log("Stabilization finished");
                    network.setOptions({physics: {enabled: false}});
                });

                // プログレスバーの処理
                network.on("stabilizationProgress", function(params) {
                    console.log(params.iterations + " / " + params.total);
                });
            });
            </script>
            </head>
            ''')

    # ネットワークの保存と表示
    # generate the full HTML and write it directly
    html_content = net.generate_html()

    # show it in‑page
    components.html(html_content, height=height, scrolling=True)


def save_static_network(tf_network: nx.DiGraph,
                       regulators_by_level: Dict[str, Set[str]],
                       target_gene: str,
                       output_path: str,
                       ppi_tfs: Optional[List[str]] = None,
                       font_size: int = 12,
                       node_size_factor: float = 1.0,
                       edge_width_factor: float = 1.0,
                       hub_centric: bool = False,
                       color_hubs: bool = True,
                       hubs: Optional[Set[str]] = None,
                       format: str = 'png',
                       node_alpha: float = 0.7,
                       edge_alpha: float = 0.6,
                       layout_type: str = "spring",
                       layout_scale: float = 3.0,
                       min_node_dist: float = 0.2) -> None:

    if format.lower() not in ['png', 'pdf']:
        raise ValueError("Format must be either 'png' or 'pdf'")
    
    if not output_path.lower().endswith(f'.{format}'):
        output_path = f"{output_path}.{format}"
    
    plt.clf()
    fig = plt.figure(figsize=(15, 12))
    
    # format パラメータを除外して他のパラメータを渡す
    visualize_static_network(
        tf_network=tf_network,
        regulators_by_level=regulators_by_level,
        target_gene=target_gene,
        ppi_tfs=ppi_tfs,
        font_size=font_size,
        node_size_factor=node_size_factor,
        edge_width_factor=edge_width_factor,
        hub_centric=hub_centric,
        color_hubs=color_hubs,
        hubs=hubs,
        node_alpha=node_alpha,
        edge_alpha=edge_alpha,
        layout_type=layout_type,
        layout_scale=layout_scale,
        min_node_dist=min_node_dist
    )
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', 
                dpi=300 if format == 'png' else 600)
    plt.close(fig)




def save_interactive_network(tf_network: nx.DiGraph,
                            regulators_by_level: Dict[str, Set[str]],
                            target_gene: str,
                            output_path: str,
                            ppi_tfs: Optional[List[str]] = None,
                            font_size: int = 12,
                            node_size_factor: float = 1.0,
                            edge_width_factor: float = 1.0,
                            height: int = 600,
                            layout_type: str = "static_spring",
                            layout_seed: int = 42,
                            hub_centric: bool = False,
                            hubs: Optional[Set[str]] = None,
                            color_hubs: bool = True) -> None:
    height_px = f"{height}px"
    
    net = Network(height=height_px, width="100%", bgcolor="#ffffff", 
                 font_color="black", directed=True)

    # ppi_tfsがNoneの場合の処理
    ppi_tfs = ppi_tfs or []
    
    # Add this before both visualization functions
    color_map = {
        'target': '#ff0000',      # Red
        'upstream_1': '#add8e6',  # Light blue
        'upstream_2': '#90ee90',  # Light green
        'downstream_1': '#ffff00', # Yellow
        'downstream_2': '#ffc0cb', # Pink
        'target_genes': '#d3d3d3', # Light gray
        'hub_node': '#FFA500',    # Orange
        'ppi_node': '#800080',    # Purple
        'ppi_edge': '#8B4513',    # Saddle Brown
        'regular_edge': '#add8e6'  # Light blue
    }
            
    # レイアウトオプションの定義
    layout_options = {
    "static_spring": {"physics": {"enabled": False}},
    "barnes_hut": {
        "physics": {
            "enabled": True,
            "barnesHut": {
                "gravitationalConstant": -500,
                "centralGravity": 0.8,
                "springLength": 80,
                "springConstant": 0.01,
                "damping": 0.5,
                "avoidOverlap": 1
            },
            "maxVelocity": 15,
            "minVelocity": 1,
            "timestep": 0.2,
            "stabilization": {
                "enabled": True,
                "iterations": 200,
                "updateInterval": 10,
                "onlyDynamicEdges": False,
                "fit": True
            }
        }
    },
"hierarchical": {
    "layout": {
        "hierarchical": {
            "enabled": True,
            "direction": "UD",
            "sortMethod": "hubsize",    # hubsizeに変更
            "levelSeparation": 75,      # 適度な縦間隔
            "nodeSpacing": 100,         # 水平間隔は維持
            "treeSpacing": 100,         # 水平間隔は維持
            "blockShifting": True,
            "edgeMinimization": True,
            "parentCentralization": True,
            "improvedLayout": True,
            "levelBias": 0.0,           # バイアスを解除
        }
    },
    "edges": {
        "smooth": {
            "enabled": True,
            "type": "cubicBezier",      # エッジの描画方法を変更
            "roundness": 0.5
        }
    },
    "physics": {
        "enabled": True,
        "hierarchicalRepulsion": {
            "centralGravity": 0.0,      # 中心への引力を無効化
            "springLength": 100,
            "springConstant": 0.01,
            "nodeDistance": 80,
            "damping": 0.3
        },
        "solver": "hierarchicalRepulsion",
        "stabilization": {
            "enabled": True,
            "iterations": 1000,         # 反復回数を増やす
            "updateInterval": 10,
            "fit": True
        }
    }
},
    "force_atlas2": {
        "physics": {
            "forceAtlas2Based": {
                "gravitationalConstant": -50,
                "centralGravity": 0.01,
                "springLength": 100,
                "springConstant": 0.08
            }
        }
    },
    "repulsion": {
        "physics": {
            "repulsion": {
                "nodeDistance": 200,
                "centralGravity": 0.2,
                "springLength": 200,
                "springConstant": 0.05
            }
        }
    }
    }
    # Hub-centricレイアウトの計算
    if hub_centric and hubs:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)
        center_x = sum(pos[node][0] for node in tf_network.nodes()) / len(tf_network.nodes())
        center_y = sum(pos[node][1] for node in tf_network.nodes()) / len(tf_network.nodes())
        
        # ハブ遺伝子とターゲット遺伝子を中央に配置
        hub_nodes = set(hubs) | {target_gene}
        n_hubs = len(hub_nodes)
        
        # 中央付近に円形に配置
        for i, hub in enumerate(hub_nodes):
            angle = 2 * np.pi * i / n_hubs
            r = 0.3  # 中心からの距離
            pos[hub] = np.array([
                center_x + r * np.cos(angle),
                center_y + r * np.sin(angle)
            ])
        
        # 他のノードは外側に配置
        for node in tf_network.nodes():
            if node not in hub_nodes:
                current_pos = pos[node]
                # 中心からの方向ベクトルを計算
                direction = current_pos - np.array([center_x, center_y])
                # 正規化して一定距離に配置
                if np.linalg.norm(direction) > 0:
                    direction = direction / np.linalg.norm(direction)
                    pos[node] = np.array([center_x, center_y]) + direction * 1.0
    else:
        pos = nx.spring_layout(tf_network, k=2, iterations=50)

    # レイアウトの選択と適用
    if layout_type in layout_options:
        net.options = layout_options[layout_type]
    
    # ノードの追加
    for node in tf_network.nodes():
        size = 20 * node_size_factor
        color = color_map['target_genes']
        
        if node == target_gene:
            color = color_map['target']
            size *= 2
        elif node in ppi_tfs:
            color = color_map['ppi_node']
            size *= 1.5
        elif hubs and node in hubs and (hub_centric or color_hubs):
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
        
        if layout_type == "static_spring":
            x, y = pos[node]
            x = x * 500
            y = y * 500
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size}, x=x, y=y, physics=False)
        else:
            net.add_node(node, label=node, color=color, size=size,
                        font={'size': font_size})

    # エッジの分類
    regular_edges = []
    ppi_edges = []
    regulation_edges = []
    
    for u, v, data in tf_network.edges(data=True):
        edge_type = data.get('type', 'normal')
        if edge_type == 'ppi':
            ppi_edges.append((u, v))
        elif edge_type == 'regulation':
            regulation_edges.append((u, v))
        else:
            regular_edges.append((u, v))

    # 1. 通常のエッジ
    for u, v in regular_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.0 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_gene, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
            arrows={'to': {'enabled': True, 'scaleFactor': 1}},
            physics=True)

    # 2. PPIエッジ（無向）
    for u, v in ppi_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        net.add_edge(u, v, width=width, color=color_map['ppi_edge'],
                    arrows={'to': {'enabled': False}},
                    physics=True)

    # 3. 制御エッジ（有向）
    for u, v in regulation_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_gene, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}},
                    physics=True)

    # JavaScriptによる制御を修正
    javascript_code = '''
<script>
window.addEventListener('load', function() {
    let stabilizationTimeout;
    let isStabilized = false;
    let forceStopTimeout;

    function disablePhysics() {
        if (!isStabilized && typeof network !== 'undefined') {
            network.setOptions({physics: {enabled: false}});
            network.stabilize(0);  // 強制的に安定化
            isStabilized = true;
        }
    }

    // 短いタイムアウト (5秒)
    stabilizationTimeout = setTimeout(disablePhysics, 5000);

    // 絶対的なタイムアウト (8秒)
    forceStopTimeout = setTimeout(function() {
        if (!isStabilized) {
            console.log("Force stopping physics after absolute timeout");
            disablePhysics();
        }
    }, 8000);

    if (typeof network !== 'undefined') {
        network.on("stabilizationIterationsDone", function() {
            clearTimeout(stabilizationTimeout);
            clearTimeout(forceStopTimeout);
            disablePhysics();
        });

        network.on("stabilizationProgress", function(params) {
            let progress = Math.round((params.iterations / params.total) * 100);
            if (progress >= 80) {  // 80%で停止
                clearTimeout(stabilizationTimeout);
                clearTimeout(forceStopTimeout);
                disablePhysics();
            }
        });

        // 任意のユーザーインタラクションで停止
        ['click', 'dragStart', 'zoom'].forEach(function(event) {
            network.on(event, function() {
                if (!isStabilized) {
                    clearTimeout(stabilizationTimeout);
                    clearTimeout(forceStopTimeout);
                    disablePhysics();
                }
            });
        });
    }
});
</script>
    '''

    # 物理演算の制御（barnes_hutの場合）
    if layout_type == "barnes_hut":
        net.html = net.html.replace('</head>',
            '''
            <script>
            window.addEventListener('load', function() {
                // 30秒後に強制的に物理シミュレーションを停止
                setTimeout(function() {
                    if (network) {
                        network.setOptions({physics: {enabled: false}});
                        console.log("Physics simulation stopped by timeout");
                    }
                }, 30000);

                // 安定化後の処理
                network.on("stabilizationIterationsDone", function() {
                    console.log("Stabilization finished");
                    network.setOptions({physics: {enabled: false}});
                });

                // プログレスバーの処理
                network.on("stabilizationProgress", function(params) {
                    console.log(params.iterations + " / " + params.total);
                });
            });
            </script>
            </head>
            ''')

    # ネットワークをファイルに保存
    #net.save_graph(output_path)
    # write the pyvis HTML straight to your chosen path
    html_content = net.generate_html()
    with open(output_path, 'w', encoding='utf-8') as out:
        out.write(html_content)


def add_save_buttons_to_main(st, filtered_network, regulators_by_level, target_gene, 
                             font_size, node_size, edge_width, viz_height, 
                             interactive_layout, static_layout, static_format, 
                             ppi_tfs=None, hub_centric=False, color_hubs=True, hubs=None, 
                             node_alpha=0.7, edge_alpha=0.6):

    """Save buttons for network visualizations"""
    import os
    import base64

    st.subheader("Save Visualizations")
    col1, col2 = st.columns(2)
    
    with col1:
        st.write("**Static Network**")
        
        try:
            with tempfile.NamedTemporaryFile(delete=False, suffix=f'.{static_format}') as tmp_file:
                save_static_network(
                    tf_network=filtered_network,
                    regulators_by_level=regulators_by_level,
                    target_gene=target_gene,
                    output_path=tmp_file.name,
                    ppi_tfs=ppi_tfs,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    hub_centric=hub_centric,
                    hubs=hubs,
                    format=static_format,
                    node_alpha=node_alpha,
                    edge_alpha=edge_alpha,
                    layout_type=static_layout
                )
                tmp_file_path = tmp_file.name
            
            with open(tmp_file_path, "rb") as file:
                base64_file = base64.b64encode(file.read()).decode('utf-8')
            
            href = f'<a href="data:application/octet-stream;base64,{base64_file}" download="{target_gene}_network.{static_format}">Download Static Network ({static_format.upper()})</a>'
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
                    target_gene=target_gene,
                    output_path=tmp_file.name,
                    ppi_tfs=ppi_tfs,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    height=viz_height,
                    layout_type=interactive_layout,
                    hub_centric=hub_centric,
                    color_hubs=color_hubs,
                    hubs=hubs
                )
                tmp_file_path = tmp_file.name
            
            with open(tmp_file_path, "rb") as file:
                base64_file = base64.b64encode(file.read()).decode('utf-8')
            
            href = f'<a href="data:text/html;base64,{base64_file}" download="{target_gene}_network.html">Download Interactive Network (HTML)</a>'
            st.markdown(href, unsafe_allow_html=True)
        
        except Exception as e:
            st.error(f"Error creating interactive network file: {e}")
        finally:
            if 'tmp_file_path' in locals() and os.path.exists(tmp_file_path):
                os.unlink(tmp_file_path)


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
    
    # Initialize session state if needed
    if 'network_analyzed' not in st.session_state:
        st.session_state.network_analyzed = False
    if 'filtered_network' not in st.session_state:
        st.session_state.filtered_network = None
    if 'regulators_by_level' not in st.session_state:
        st.session_state.regulators_by_level = None
    if 'ppi_tfs' not in st.session_state:
        st.session_state.ppi_tfs = []
    if 'hub_nodes_for_viz' not in st.session_state:
        st.session_state.hub_nodes_for_viz = None

    with st.expander("Cache Control when previous data remain"):
        if st.button("Clear Cache and Reset"):
            st.cache_data.clear()
            for key in list(st.session_state.keys()):
                if key.startswith("adjacencies_") or key in ["normalized_adjacencies_data", "raw_adjacencies_data"]:
                    del st.session_state[key]
            st.rerun()

    st.markdown("##### Select edge weight data source:") #スペースを入れないと次の選択がうまくできない!!!!!!!!!!!!!!!!!!!!!!!!

    # エッジに使用するデータの選択
    edge_weight_source = st.radio(
        "Data type:",
        ["CSI Data", "Adjacencies Data", "Equal Weights (No edge weight data)"],
        key="edge_data_source_radio",
        help="""CSI:2つのレギュロン間の活性パターンの類似性がどれだけ特異的か 高いCSIを共有するレギュロンは下流の遺伝子を共同で制御している可能性が高い
        　　Adjacencies:TFとその標的との間の調節関係の強さを示す。高い値は、より強い調節関係を示唆しTFの発現変化が標的遺伝子の発現に大きな影響を与える可能性　より信頼性の高い遺伝子制御関係を表す
        """
    )


    st.markdown("##### Upload Regulon Data")
    regulon_file = st.file_uploader(
        "regulons.seurat_object.filtered_by_pySCENIC.txt or regulons.seurat_object.filtered_by_AUCell_exploreThresholds.txt or regulons.seurat_object.txt", 
        type=['txt', 'tsv'],
        help ='regulons.seurat_object.filtered_by_pySCENIC.txtはpySCENICでフィルタリングされたもの、regulons.seurat_object.filtered_by_AUCell_exploreThresholds.txtではregulonが活性化している細胞が10以下と判断されたregulonは除かれている'
    )

    # 選択に応じたファイルアップロード
    edge_weights_data = None
    csi_data = None
    adjacencies_data = None
    csi_file = None
    adjacencies_file = None
    use_equal_weights = False

    if edge_weight_source == "CSI Data":
        csi_file = st.file_uploader(
            "**Upload CSI Data (e.g., csi_results.csv)**", 
            type=['csv']
        )
        if csi_file:
            csi_data = load_csi_data(csi_file)
            edge_weights_data = csi_data
            st.write("CSI: connection specificity index. Please calculate using SCENIC CSI app.")


    elif edge_weight_source == "Adjacencies Data":
        adjacencies_file = st.file_uploader(
            "**Upload Adjacencies Data (e.g., adjacencies.tsv)**", 
            type=['tsv']
        )

        if adjacencies_file:
            # ファイルが変更されたかどうかをチェック
            file_content = adjacencies_file.getvalue()
            content_hash = hashlib.md5(file_content).hexdigest()
            
            if 'adjacencies_hash' not in st.session_state or st.session_state.adjacencies_hash != content_hash:
                with st.spinner('Processing adjacencies file (this may take a while)...'):
                    # ロードと正規化の両方をセッションステートに保存
                    raw_adjacencies_data = load_adjacencies_data(adjacencies_file)
                    
                    # 正規化処理のための追加のスピナー
                    with st.spinner('Normalizing adjacencies importance scores...'):
                        st.session_state.normalized_adjacencies_data = normalize_adjacency_importance(raw_adjacencies_data)
                    
                    # 元のデータも必要な場合はこちらも保存
                    st.session_state.raw_adjacencies_data = raw_adjacencies_data
                    st.session_state.adjacencies_hash = content_hash
                    st.info(f"Processed and normalized adjacencies data with hash: {content_hash[:8]}")
            
            # キャッシュされたデータを使用
            adjacencies_data = st.session_state.normalized_adjacencies_data
            edge_weights_data = adjacencies_data
            st.info("Adjacencies importance scores have been log10-normalized and scaled to 0-1 range.")


    else:  # Equal Weights
        use_equal_weights = True
        edge_weights_data = defaultdict(dict)
        st.info("Using equal weights for all edges (no weight data needed).")


    # フィルタリングの選択（データがロードされた場合のみ表示）
    if (edge_weight_source in ["CSI Data", "Adjacencies Data"] and edge_weights_data is not None) or use_equal_weights:
        
        use_rank_based = st.checkbox(
            "Use rank-based weighting", 
            value=False,
            help="Convert all edge weights to ranks (1=highest, 0=lowest)"
        )

        if use_rank_based:
            # ランクベース変換のサブオプション
            rank_ties = st.radio(
                "Rank calculation method:",
                ["Same score = Same rank", "All ranks unique (continuous)"],
                help="Choose how to handle identical scores"
            )
            
            use_ties = rank_ties == "Same score = Same rank"
            
            # ハッシュキーの計算（入力データと設定に基づく）
            edge_data_hash = hashlib.md5(str(edge_weights_data).encode()).hexdigest()
            rank_config_hash = f"{edge_data_hash}_{use_ties}"
            
            # セッションステートにキャッシュキーが存在するか確認
            if ('rank_based_data' not in st.session_state or 
                st.session_state.get('rank_config_hash', '') != rank_config_hash):
                
                # キャッシュがない場合は変換を実行
                with st.spinner('Converting edge weights to rank-based...'):
                    # キャッシュ付き関数を使用して変換
                    ranked_data = convert_to_rank_based(edge_weights_data, use_ties=use_ties)
                    
                    # セッションステートに保存
                    st.session_state.rank_based_data = ranked_data
                    st.session_state.rank_config_hash = rank_config_hash
                    
                if use_ties:
                    st.info("Edge weights converted to rank-based with identical scores getting the same rank.")
                else:
                    st.info("Edge weights converted to continuous rank-based (all weights are unique).")
            else:
                # キャッシュから読み込み
                if use_ties:
                    st.info("Using cached rank-based weights (same scores = same rank).")
                else:
                    st.info("Using cached continuous rank-based weights.")
            
            # セッションステートから変換済みデータを使用
            edge_weights_data = st.session_state.rank_based_data
            
            if edge_weight_source == "CSI Data":
                csi_data = edge_weights_data
            elif edge_weight_source == "Adjacencies Data":
                adjacencies_data = edge_weights_data

        st.subheader("Edge Filtering Options")
        
        # フィルタリングを行うかどうか
        apply_filtering = st.checkbox(
            "Apply filtering to edges", 
            value=False,
            help="Filter network edges based on weight or count"
        )
        
        if apply_filtering:
            # フィルタリング設定
            with st.expander("Filtering Settings", expanded=True):
                filtering_method = st.radio(
                    "Filtering Method",
                    options=["Weight Threshold", "Max TFs"],
                    key="filtering_method_radio",
                    help="Choose how to filter TFs"
                )
                
                if filtering_method == "Weight Threshold":
                    weight_threshold = st.slider(
                        f"{edge_weight_source.split()[0]} Threshold", 
                        min_value= 0.0,
                        max_value= 1.0,
                        value = 0.1,
                        help="Higher values keep only stronger connections (0-1 range)"
                    )

                    fig, ax = plt.subplots()
                    ax.hist(edge_weights_data, bins=20, edgecolor='black', alpha=0.7)     
                    if edge_weight_source == "CSI Data":
                        visualize_single_distribution(
                            data=edge_weights_data,
                            threshold=weight_threshold,
                            title="Distribution of CSI Scores"
                        )
                    else:
                        visualize_single_distribution(
                            data=edge_weights_data,
                            threshold=weight_threshold,
                            title="Distribution of Adjacencies Scores"
                        )
                else:
                    max_tfs = st.slider("Max TFs per level", 1, 50, 20)
                    weight_threshold = None
            
            # 追加のフィルタリングオプション（シーケンシャルフィルタリング）
            apply_secondary_filter = st.checkbox(
                "Apply secondary filtering", 
                value=False,
                help="Apply additional filtering using another data source"
            )
            
            if apply_secondary_filter:
                with st.expander("Secondary Filtering", expanded=True):
                    # 現在のデータソースに基づいて、もう一方を選択
                    if edge_weight_source == "CSI Data":
                        secondary_file = st.file_uploader(
                            "**Upload Adjacencies Data for Secondary Filtering**", 
                            type=['tsv']
                        )

                        if secondary_file:
                            # ファイルが変更されたかどうかをチェック
                            secondary_content = secondary_file.getvalue()
                            secondary_hash = hashlib.md5(secondary_content).hexdigest()
                            
                            if 'secondary_hash' not in st.session_state or st.session_state.secondary_hash != secondary_hash:
                                with st.spinner('Processing secondary adjacencies file...'):
                                    # ロードと正規化の両方をセッションステートに保存
                                    raw_secondary_data = load_adjacencies_data(secondary_file)
                                    
                                    with st.spinner('Normalizing secondary importance scores...'):
                                        st.session_state.normalized_secondary_data = normalize_adjacency_importance(raw_secondary_data)
                                    
                                    # 元のデータも保存
                                    st.session_state.raw_secondary_data = raw_secondary_data
                                    st.session_state.secondary_hash = secondary_hash
                                    st.info("Adjacencies importance scores have been log10-normalized and scaled to 0-1 range.")
                            
                            # キャッシュされたデータを使用
                            secondary_data = st.session_state.normalized_secondary_data
                            
                            secondary_threshold = st.slider(
                                "Adjacencies Threshold", 
                                0.0, 1.0, 0.1,
                                help="Keep only relationships with Adjacencies score above this threshold"
                            )
                            
                            # フィルタリング順序
                            filter_order = st.radio(
                                "Filtering Order",
                                ["Filter by CSI, then by Adjacencies", "Filter by Adjacencies, then by CSI"],
                                key="filter_order_radio"
                            )
                            
                            if secondary_file:
                                # Visualize score distributions
                                if filter_order == "Filter by CSI, then by Adjacencies":
                                    primary_is_csi = True
                                    primary_threshold = weight_threshold if filtering_method == "Weight Threshold" else 0.0
                                    visualize_score_distributions(
                                        csi_data,
                                        secondary_data,
                                        primary_threshold,
                                        secondary_threshold
                                    )
                                    
                                    # シーケンシャルフィルタリングを適用
                                    edge_weights_data = filter_sequentially(
                                        csi_data, 
                                        secondary_data,
                                        primary_threshold,
                                        secondary_threshold,
                                        primary_is_csi=True
                                    )
                                else:
                                    primary_is_csi = False
                                    primary_threshold = secondary_threshold
                                    visualize_score_distributions(
                                        csi_data,
                                        secondary_data,
                                        weight_threshold if filtering_method == "Weight Threshold" else 0.0,
                                        secondary_threshold
                                    )
                                    
                                    # シーケンシャルフィルタリングを適用
                                    edge_weights_data = filter_sequentially(
                                        csi_data, 
                                        secondary_data,
                                        secondary_threshold,
                                        weight_threshold if filtering_method == "Weight Threshold" else 0.0,
                                        primary_is_csi=False
                                    )
                                
                                st.info(f"Applied sequential filtering according to the selected order.")
                                
                    # Adjacenciesデータを主に使用している場合
                    elif edge_weight_source == "Adjacencies Data":
                        secondary_file = st.file_uploader(
                            "**Upload CSI Data for Secondary Filtering**", 
                            type=['csv']
                        )
                        if secondary_file:
                            secondary_data = load_csi_data(secondary_file)
                            st.write("CSI: connection specificity index.")
                            
                            secondary_threshold = st.slider(
                                "CSI Threshold", 
                                0.0, 1.0, 0.3,
                                help="Keep only relationships with CSI score above this threshold"
                            )
                            
                            # フィルタリング順序
                            filter_order = st.radio(
                                "Filtering Order",
                                ["Filter by Adjacencies, then by CSI", "Filter by CSI, then by Adjacencies"],
                                key="filter_order_radio"
                            )
                            
                            if secondary_file:
                                # Visualize score distributions
                                if filter_order == "Filter by Adjacencies, then by CSI":
                                    primary_is_csi = False
                                    primary_threshold = weight_threshold if filtering_method == "Weight Threshold" else 0.0
                                    visualize_score_distributions(
                                        secondary_data,
                                        adjacencies_data,
                                        secondary_threshold,
                                        primary_threshold
                                    )
                                    
                                    # シーケンシャルフィルタリングを適用
                                    edge_weights_data = filter_sequentially(
                                        secondary_data, 
                                        adjacencies_data,
                                        secondary_threshold,
                                        primary_threshold,
                                        primary_is_csi=True
                                    )
                                else:
                                    primary_is_csi = True
                                    primary_threshold = secondary_threshold
                                    visualize_score_distributions(
                                        secondary_data,
                                        adjacencies_data,
                                        secondary_threshold,
                                        weight_threshold if filtering_method == "Weight Threshold" else 0.0
                                    )
                                    
                                    # シーケンシャルフィルタリングを適用
                                    edge_weights_data = filter_sequentially(
                                        secondary_data, 
                                        adjacencies_data,
                                        secondary_threshold,
                                        weight_threshold if filtering_method == "Weight Threshold" else 0.0,
                                        primary_is_csi=False
                                    )
                                
                                st.info(f"Applied sequential filtering according to the selected order.")
        else:
            # フィルタリングなしの場合のデフォルト値
            weight_threshold = 0.0
            max_tfs = None


    # 2. データ読み込み条件の修正
    # 必要なファイルがアップロードされたことを確認する条件を修正
    required_files_uploaded = regulon_file and (
        (edge_weight_source == "CSI Data" and csi_file) or
        (edge_weight_source == "Adjacencies Data" and adjacencies_file) or
        (edge_weight_source == "Sequential (CSI + Adjacencies)" and csi_file and adjacencies_file) or
        use_equal_weights
    )

    
    if required_files_uploaded:
        # データの読み込み
        regulon_data, all_genes = load_regulon_data(regulon_file)
        
        # TFリストを自動的に読み込み
        all_tfs = load_tfs_list(regulon_file)
        # 生物種の判定
        species = determine_species(regulon_file)

        # Get visualization settings
        static_layout, interactive_layout, layout_scale, min_node_dist, hub_centric, color_hubs, enable_cleanup = add_visualization_settings(st)

        # その他の可視化パラメータ
        with st.expander("Additional Visualization Parameters", expanded=False):
            col1, col2, col3 = st.columns(3)
            with col1:
                font_size = st.slider("Font Size", 10, 30, 12)
            with col2:
                node_size = st.slider("Node Size", 0.5, 2.0, 1.0)
            with col3:
                node_alpha = st.slider("Node Alpha", 0.1, 1.0, 0.7)
            
            col4, col5, col6 = st.columns(3)
            with col4:
                edge_width = st.slider("Edge Width", 0.2, 2.0, 1.0)
            with col5:
                edge_alpha = st.slider("Edge Alpha", 0.1, 1.0, 0.6)
            with col6:
                viz_height = st.slider("Graph Height", 400, 1000, 800)
            
            static_format = st.selectbox("Static network format to download", ["pdf", 'png'])

        with st.form("input_groups and batch"):     
            st.subheader("Basic Settings")
            col1, col2, col3 = st.columns(3)
            with col1:
                target_gene = st.selectbox("Select TF of Interest/Target Gene", 
                                         sorted(list(all_genes)))
            with col2:
                max_upstream = st.number_input("Max Upstream Level", 
                                             min_value=0, max_value=4, value=2)
            with col3:
                max_downstream = st.number_input("Max Downstream Level", 
                                               min_value=0, max_value=4, value=2)
            
            # ネットワークモードの設定
            network_mode = st.checkbox(
                "Enable Network Mode", 
                value=False,
                help="If enabled, search both upstream and downstream connections for each TF"
            )
            
            # フィルタリング設定 - Only show if not using sequential filtering
            if edge_weight_source != "Sequential (CSI + Adjacencies)":
                st.subheader("Filtering Settings")
                with st.expander("Edge Weight Filtering", expanded=True):
                    filtering_method = st.radio(
                        "Filtering Method",
                        options=["Weight Threshold", "Max TFs"],
                        help="Choose how to filter TFs"
                    )
                    
                    if filtering_method == "Weight Threshold":
                        weight_threshold = st.slider(
                            "Weight Threshold", 
                            0.0, 1.0, 0.0, 
                            help="Higher values keep only stronger connections (0-1 range)"
                        )
                        max_tfs = None
                    else:
                        max_tfs = st.slider("Max TFs per level", 1, 50, 20)
                        weight_threshold = None
            else:
                # For sequential filtering, set these values
                weight_threshold = 0.0  # No additional filtering needed since we already filtered
                max_tfs = None

            # TRRUST設定
            trrust_options = add_trrust_options(st, species)

            # ChIP-Atlas設定
            chip_atlas_options = add_chip_atlas_options(st, species)

            # STRING-db設定
            string_db_options = add_string_db_options(st, species)

            # PPI Extension settings
            st.subheader("PPI Extension")
            st.write("Add TF interacts with the target gene.")
            with st.expander("PPI Extension Settings", expanded=False):
                ppi_options = {
                    "enable": st.checkbox(
                        "Extend network with STRING-db PPI data", 
                        help='Target geneと相互作用し、target geneの標的遺伝子をregulon（標的遺伝子）に含む転写因子を追加する'
                    ),
                    "score_threshold": st.slider(
                        "PPI Score Threshold", 
                        min_value=0.0, 
                        max_value=1.0, 
                        value=0.4, 
                        step=0.1
                    )
                }

            # Hub Selection settings
            st.subheader("Hub Selection for Visualization")
            with st.expander("Hub Settings", expanded=False):
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
                    

            submitted = st.form_submit_button("Submit")
            
            if submitted:
                st.session_state.network_analyzed = False  # Reset analysis state when form is submitted
        
        st.markdown("#### To change options, you need to press submit everytime.")
        
        # Initial visualization settings if not already analyzed
        if not st.session_state.network_analyzed:
            static_layout = "spring"
            interactive_layout = "static_spring"
            font_size = 12
            node_size = 1.0
            node_alpha = 0.7
            edge_width = 1.0
            edge_alpha = 0.6
            viz_height = 800
            static_format = "pdf"

        # Analysis and Visualization
        if st.button("Start Analysis", type="primary") or (st.session_state.network_analyzed and 'static_layout' in locals()):
            if not st.session_state.network_analyzed:
                # Only perform network analysis if not already done
                with st.spinner('Analyzing network...'):
                    # Build network
                    tf_network, regulators_by_level = build_expanded_tf_network(
                        regulon_data, edge_weights_data, all_tfs,
                        target_gene, max_upstream, max_downstream,
                        network_mode, use_equal_weights, weight_threshold, max_tfs
                    )
        
                    # Apply filters - use the renamed function
                    filtered_network = filter_network_by_weights(tf_network, weight_threshold)

                    # For sequential filtering, display some statistics
                    if edge_weight_source == "Sequential (CSI + Adjacencies)" and csi_data and adjacencies_data:
                        st.subheader("Sequential Filtering Statistics")
                        
                        # Calculate statistics
                        total_possible_relationships = 0
                        total_csi_relationships = 0
                        total_final_relationships = 0
                        
                        all_tfs = set(list(csi_data.keys()) + list(adjacencies_data.keys()))
                        for tf1 in all_tfs:
                            for tf2 in all_tfs:
                                if tf1 != tf2:
                                    total_possible_relationships += 1
                                    
                                    # Count relationships passing CSI threshold
                                    csi_score = csi_data.get(tf1, {}).get(tf2, 0.0)
                                    if csi_score >= csi_threshold:
                                        total_csi_relationships += 1
                                        
                                        # Count relationships passing both filters
                                        importance = adjacencies_data.get(tf1, {}).get(tf2, 0.0)
                                        if importance >= importance_threshold:
                                            total_final_relationships += 1
                        
                        # Display statistics
                        col1, col2, col3 = st.columns(3)
                        with col1:
                            st.metric("Total Possible Relationships", total_possible_relationships)
                        with col2:
                            st.metric("Passed CSI Filter", total_csi_relationships)
                            st.metric("% Passed CSI", f"{total_csi_relationships / total_possible_relationships * 100:.1f}%")
                        with col3:
                            st.metric("Passed Both Filters", total_final_relationships)
                            st.metric("% Passed Both", f"{total_final_relationships / total_possible_relationships * 100:.1f}%")
                        
                        st.info(f"Sequential filtering reduced relationships by {100 - (total_final_relationships / total_possible_relationships * 100):.1f}%")

                    # Initialize ppi_tfs as empty list
                    ppi_tfs = []
                    if ppi_options["enable"]:
                        tf_list = list(regulon_data.keys())
                        string_data = get_string_interactions(
                            tf_list,
                            species,
                            string_db_options["score_threshold"]
                        )

                        if string_data:
                            filtered_network, ppi_edges, ppi_tfs = extend_network_with_ppi(
                                filtered_network,
                                target_gene, 
                                regulon_data, 
                                string_data,
                                ppi_options["score_threshold"]
                            )
                    
                    # Apply additional filters if enabled
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

                    # Cleanup networkを条件付きで実行
                    if enable_cleanup:
                        filtered_network, disconnected_edges = cleanup_network(filtered_network, target_gene)
                                        
                    # Network analysis if enough nodes
                    if filtered_network.number_of_nodes() >= 2:
                        st.markdown("---")
                        st.subheader("Network Analysis Results")
                        
                        # Calculate centrality measures
                        weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector = \
                            calculate_weighted_centrality(filtered_network, use_equal_weights)
                        
                        # Display centrality analysis
                        st.write("### Centrality Analysis")
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("**Top Weighted Degree Centrality Nodes:**")
                            for node, value in sorted(weighted_degree.items(), 
                                                    key=lambda x: x[1], 
                                                    reverse=True)[:5]:
                                st.write(f"{node}: {value:.4f}")
                            
                            st.write("\n**Top Weighted PageRank Nodes:**")
                            for node, value in sorted(weighted_pagerank.items(), 
                                                    key=lambda x: x[1], 
                                                    reverse=True)[:5]:
                                st.write(f"{node}: {value:.4f}")
                        
                        with col2:
                            st.write("**Top Weighted Betweenness Centrality Nodes:**")
                            for node, value in sorted(weighted_between.items(), 
                                                    key=lambda x: x[1], 
                                                    reverse=True)[:5]:
                                st.write(f"{node}: {value:.4f}")
                            
                            st.write("\n**Top Weighted Eigenvector Centrality Nodes:**")
                            for node, value in sorted(weighted_eigenvector.items(), 
                                                    key=lambda x: x[1], 
                                                    reverse=True)[:5]:
                                st.write(f"{node}: {value:.4f}")
 

                        # Hub analysis
                        st.write("\n### Hub Analysis")
                        # Use the same percentile threshold as visualization
                        threshold = percentile_threshold / 100
                        weighted_hubs = identify_hubs_weighted(filtered_network, weight_threshold=threshold)
                        weighted_bottlenecks = identify_bottlenecks_weighted(filtered_network, 
                                                                           betweenness_threshold=threshold)
                        weighted_common_hubs = weighted_hubs & weighted_bottlenecks
                        
                        col1, col2 = st.columns(2)
                        with col1:
                            st.write(f"**Weighted Degree-based Hubs (>{percentile_threshold}th percentile):**")
                            st.write(', '.join(sorted(weighted_hubs)))
                            
                            st.write(f"\n**Weighted Bottleneck Nodes (>{percentile_threshold}th percentile):**")
                            st.write(', '.join(sorted(weighted_bottlenecks)))
                        
                        with col2:
                            st.write(f"**Common Weighted Hubs (>{percentile_threshold}th percentile):**")
                            st.write(', '.join(sorted(weighted_common_hubs)))

                        
                        # Identify hub nodes for visualization
                        hub_method = hub_methods[selected_hub_method]
                        centrality_data = (weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector)
                        hub_nodes_for_viz = get_hub_nodes(
                            hub_method, filtered_network, centrality_data,
                            n=top_n, percentile=percentile_threshold
                        )
                        
                        # Display selected hubs for visualization
                        st.write("\n### Selected Hubs for Visualization")
                        st.write(f"**Method: {selected_hub_method}**")
                        st.write(', '.join(sorted(hub_nodes_for_viz)))
                        
                        # Store results in session state
                        st.session_state.filtered_network = filtered_network
                        st.session_state.regulators_by_level = regulators_by_level
                        st.session_state.ppi_tfs = ppi_tfs
                        st.session_state.hub_nodes_for_viz = hub_nodes_for_viz
                        st.session_state.network_analyzed = True
                    else:
                        st.warning("Not enough nodes for analysis after filtering.")
                        return

            # Visualization code (runs on both initial analysis and layout changes)
            if st.session_state.filtered_network.number_of_edges() > 0:
                st.subheader("Network Visualization")
                
                # Visualization calls
                visualize_static_network(
                    tf_network=st.session_state.filtered_network,
                    regulators_by_level=st.session_state.regulators_by_level,
                    target_gene=target_gene,
                    ppi_tfs=st.session_state.ppi_tfs,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    hub_centric=hub_centric,
                    color_hubs=color_hubs,
                    hubs=st.session_state.hub_nodes_for_viz,
                    node_alpha=node_alpha,
                    edge_alpha=edge_alpha,
                    layout_type=static_layout,
                    layout_scale=layout_scale,
                    min_node_dist=min_node_dist
                )
                
                create_interactive_network(
                    st.session_state.filtered_network,
                    st.session_state.regulators_by_level,
                    target_gene,
                    st.session_state.ppi_tfs,
                    font_size,
                    node_size,
                    edge_width,
                    viz_height,
                    layout_type=interactive_layout,
                    hub_centric=hub_centric,
                    color_hubs=color_hubs,
                    hubs=st.session_state.hub_nodes_for_viz
                )
                
                # Save buttons
                add_save_buttons_to_main(
                    st=st,
                    filtered_network=st.session_state.filtered_network,
                    regulators_by_level=st.session_state.regulators_by_level,
                    target_gene=target_gene,
                    font_size=font_size,
                    node_size=node_size,
                    edge_width=edge_width,
                    viz_height=viz_height,
                    interactive_layout=interactive_layout,
                    static_layout=static_layout,
                    static_format=static_format,
                    ppi_tfs=st.session_state.ppi_tfs,
                    hub_centric=hub_centric,
                    hubs=st.session_state.hub_nodes_for_viz,
                    node_alpha=node_alpha,
                    edge_alpha=edge_alpha
                )
            else:
                st.warning("No edges remain after filtering. Try adjusting the parameters.")
    else:
        if use_equal_weights:
            if not regulon_file:
                st.info("Please upload Regulon data file to start the analysis.")
        else:
            if edge_weight_source == "Sequential (CSI + Adjacencies)":
                if not (csi_file and adjacencies_file):
                    st.info("Please upload both CSI and Adjacencies data files to use sequential filtering.")
            elif edge_weight_source == "CSI Data" and not csi_file:
                st.info("Please upload CSI data file to start the analysis.")
            elif edge_weight_source == "Adjacencies Data" and not adjacencies_file:
                st.info("Please upload Adjacencies data file to start the analysis.")
            
            if not regulon_file:
                st.info("Please also upload Regulon data file.")


if __name__ == "__main__":
    main()