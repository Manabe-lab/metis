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

def cleanup_network_multi(network: nx.DiGraph, target_genes: List[str]) -> Tuple[nx.DiGraph, List[Tuple[str, str]]]:
    """
    複数のターゲット遺伝子に対してネットワークをクリーニング
    各ターゲット遺伝子に対して独立してパスの確認を行う
    
    パスの条件:
    1. いずれかのターゲット遺伝子へ到達するパス上のノード（上流）
    2. いずれかのターゲット遺伝子からの到達可能なパス上のノード（下流）
    """
    cleaned_network = network.copy()
    nodes_to_keep = set(target_genes)  # 全ターゲット遺伝子は保持
    
    # 1. 各ターゲット遺伝子について個別にパスを確認
    for target_gene in target_genes:
        if target_gene not in network.nodes():
            continue
        # 上流のパスを確認
        for node in network.nodes():
            if node != target_gene:
                try:
                    # ノードからターゲットへのパスが存在するか確認
                    if nx.has_path(network, node, target_gene):
                        path = nx.shortest_path(network, node, target_gene)
                        nodes_to_keep.update(path)
                except nx.NetworkXError:
                    # エラーが発生した場合はスキップ
                    continue
        # 下流のパスを確認
        # target_geneからの到達可能なノードを全て追加
        try:
            descendants = nx.descendants(network, target_gene)
            nodes_to_keep.update(descendants)
        except nx.NetworkXError:
            continue
    
    # 2. 保持するノードからの下流ノードも追加
    current_nodes = nodes_to_keep.copy()
    for node in current_nodes:
        try:
            descendants = nx.descendants(network, node)
            nodes_to_keep.update(descendants)
        except nx.NetworkXError:
            continue
    
    # 3. 削除するノードを特定
    nodes_to_remove = set(network.nodes()) - nodes_to_keep
    
    # 4. 除外されるエッジを記録
    removed_edges = []
    for u, v in network.edges():
        if u in nodes_to_remove or v in nodes_to_remove:
            removed_edges.append((u, v))
    
    # 5. ノードを削除
    cleaned_network.remove_nodes_from(nodes_to_remove)
    
    # 6. 結果を表示
    if removed_edges:
        st.markdown("### Removing disconnected subnetworks:")
        st.write(f"Nodes removed: {len(nodes_to_remove)}")
        st.write(f"Nodes kept: {len(nodes_to_keep)}")
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


def create_interactive_multi_target_network(
    tf_network: nx.DiGraph,
    regulators_by_level: Dict[str, Set[str]],
    target_genes: List[str],
    ppi_tfs: Optional[List[str]] = None,
    font_size: int = 12,
    node_size_factor: float = 1.0,
    edge_width_factor: float = 1.0,
    height: int = 600,
    layout_type: str = "static_spring",
    layout_seed: int = 42,
    hub_centric: bool = False,
    color_hubs: bool = True,
    hubs: Optional[Set[str]] = None
) -> None:
    height_px = f"{height}px"
    net = Network(height=height_px, width="100%", bgcolor="#ffffff", 
                 font_color="black", directed=True)
                 
    # ppi_tfsがNoneの場合の処理
    ppi_tfs = ppi_tfs or []

    # カラーマップの定義
    target_colors = [
        '#ff0000',  # Red
        '#ff6b6b',  # Light red
        '#ff9999',  # Lighter red
        '#ffb366',  # Orange
        '#ffcc00',  # Yellow
    ]
    
    color_map = {
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

    # ターゲット遺伝子の色を追加
    for i, target in enumerate(target_genes):
        color_idx = min(i, len(target_colors) - 1)
        color_map[f'target_{target}'] = target_colors[color_idx]
    
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
                    "sortMethod": "hubsize",
                    "levelSeparation": 75,
                    "nodeSpacing": 100,
                    "treeSpacing": 100,
                    "blockShifting": True,
                    "edgeMinimization": True,
                    "parentCentralization": True,
                    "improvedLayout": True,
                    "levelBias": 0.0,
                }
            },
            "edges": {
                "smooth": {
                    "enabled": True,
                    "type": "cubicBezier",
                    "roundness": 0.5
                }
            },
            "physics": {
                "enabled": True,
                "hierarchicalRepulsion": {
                    "centralGravity": 0.0,
                    "springLength": 100,
                    "springConstant": 0.01,
                    "nodeDistance": 80,
                    "damping": 0.3
                },
                "solver": "hierarchicalRepulsion",
                "stabilization": {
                    "enabled": True,
                    "iterations": 1000,
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
        hub_nodes = set(hubs) | set(target_genes)
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
                direction = current_pos - np.array([center_x, center_y])
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
        
        if node in target_genes:
            color = color_map[f'target_{node}']
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

    # エッジの追加
    for u, v in regular_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.0 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}})

    for u, v in ppi_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        net.add_edge(u, v, width=width, color=color_map['ppi_edge'],
                    arrows={'to': {'enabled': False}},
                    dashes=True)

    for u, v in regulation_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}})

    # JavaScriptの制御コード
    stabilization_script = '''
    <script>
    window.addEventListener('load', function() {
        let stabilizationTimeout;
        let isStabilized = false;
        let forceStopTimeout;

        function disablePhysics() {
            if (!isStabilized && typeof network !== 'undefined') {
                network.setOptions({physics: {enabled: false}});
                network.stabilize(0);
                isStabilized = true;
            }
        }

        stabilizationTimeout = setTimeout(disablePhysics, 5000);
        forceStopTimeout = setTimeout(function() {
            if (!isStabilized) {
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
                if (progress >= 80) {
                    clearTimeout(stabilizationTimeout);
                    clearTimeout(forceStopTimeout);
                    disablePhysics();
                }
            });

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

    # ネットワークの保存と表示
    with tempfile.NamedTemporaryFile(delete=False, suffix='.html', mode='w', encoding='utf-8') as f:
        net.save_graph(f.name)
        content = open(f.name, 'r', encoding='utf-8').read()
        content = content.replace('</head>', f'{stabilization_script}</head>')
        components.html(content, height=height, scrolling=True)
    os.unlink(f.name)

def build_multi_target_network(
    regulon_data: Dict[str, Set[str]], 
    csi_data: Dict[str, Dict[str, float]], 
    all_tfs: Set[str],
    target_genes: List[str],
    max_upstream: int,
    max_downstream: int,
    network_mode: bool,
    use_equal_weights: bool = False,
    csi_threshold: float = None,
    max_tfs_per_level: int = None
) -> Tuple[nx.DiGraph, Dict[str, Set[str]]]:
    """複数のターゲット遺伝子に対してネットワークを構築"""
    combined_network = nx.DiGraph()
    combined_regulators = {
        'target_genes': set(target_genes),
        'upstream_1': set(),
        'upstream_2': set(),
        'downstream_1': set(),
        'downstream_2': set(),
        'additional_tfs': set()
    }
    
    for target_gene in target_genes:
        network, regulators = build_expanded_tf_network(
            regulon_data=regulon_data,
            csi_data=csi_data,
            all_tfs=all_tfs,
            target_gene=target_gene,
            max_upstream=max_upstream,
            max_downstream=max_downstream,
            network_mode=network_mode,
            use_equal_weights=use_equal_weights,
            csi_threshold=csi_threshold,
            max_tfs_per_level=max_tfs_per_level
        )
        
        # ネットワークの統合
        combined_network.add_nodes_from(network.nodes())
        for u, v, data in network.edges(data=True):
            if combined_network.has_edge(u, v):
                # エッジが既に存在する場合は重みを最大値で更新
                existing_weight = combined_network[u][v]['weight']
                new_weight = data['weight']
                combined_network[u][v]['weight'] = max(existing_weight, new_weight)
            else:
                # 新しいエッジを追加
                combined_network.add_edge(u, v, **data)
        
        # regulatorの統合
        for level, tfs in regulators.items():
            if level != 'target_genes':  # target_genes は既に設定済み
                combined_regulators[level].update(tfs)
    
    return combined_network, combined_regulators

def build_expanded_tf_network(
    regulon_data: Dict[str, Set[str]], 
    csi_data: Dict[str, Dict[str, float]], 
    all_tfs: Set[str],
    target_gene: str,
    max_upstream: int,
    max_downstream: int,
    network_mode: bool,
    use_equal_weights: bool = False,
    csi_threshold: float = None,
    max_tfs_per_level: int = None
) -> Tuple[nx.DiGraph, Dict[str, Set[str]]]:
    """単一ターゲット遺伝子のネットワーク構築（既存の関数を修正）"""
    tf_network = nx.DiGraph()
    regulators_by_level = {
        'target_genes': {target_gene},
        'upstream_1': set(),
        'upstream_2': set(),
        'downstream_1': set(),
        'downstream_2': set(),
        'additional_tfs': set()
    }
    
    def get_connected_tfs_for_level(current_tf: str, direction: str) -> List[Tuple[str, float, str]]:
        connected_tfs = []
        
        for tf, genes in regulon_data.items():
            if direction in ['upstream', 'both']:
                if current_tf in genes and tf in all_tfs:
                    weight = 1.0 if use_equal_weights else csi_data[tf].get(current_tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'upstream'))
            
            if direction in ['downstream', 'both']:
                if current_tf in regulon_data and tf in regulon_data[current_tf] and tf in all_tfs:
                    weight = 1.0 if use_equal_weights else csi_data[current_tf].get(tf, 0.0)
                    if csi_threshold is None or weight >= csi_threshold:
                        connected_tfs.append((tf, weight, 'downstream'))
        
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
    if target_gene in regulon_data and max_downstream > 0:
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

def get_network_layout_multi(
    tf_network: nx.DiGraph,
    target_genes: List[str],
    layout_type: str = "spring",
    hub_centric: bool = False,
    hubs: Optional[Set[str]] = None,
    scale: float = 3.0,
    k: float = 2.0,
    min_dist: float = 0.2,
    iterations: int = 100
) -> Dict[str, np.ndarray]:
    """複数のターゲット遺伝子に対応したレイアウトの生成
    
    Args:
        tf_network: ネットワークグラフ
        target_genes: ターゲット遺伝子のリスト
        layout_type: レイアウトタイプ
        hub_centric: ハブノードを中心に配置するかどうか
        hubs: ハブノードのセット
        scale: レイアウトのスケール
        k: ノード間の距離係数
        min_dist: 最小ノード間距離
        iterations: レイアウトアルゴリズムの反復回数
        
    Returns:
        Dict[str, np.ndarray]: ノードの位置情報
    """
    # 基本レイアウトの生成
    if layout_type == "spring":
        pos = nx.spring_layout(tf_network, k=k, scale=scale, iterations=iterations)
    elif layout_type == "kamada_kawai":
        try:
            pos = nx.kamada_kawai_layout(tf_network, scale=scale, weight='weight')
            pos = adjust_node_positions(pos, min_dist=min_dist)
        except:
            pos = nx.spring_layout(tf_network, k=k, scale=scale)
    elif layout_type == "fruchterman_reingold":
        pos = nx.fruchterman_reingold_layout(tf_network, k=k, scale=scale, iterations=iterations)
    elif layout_type == "spectral":
        try:
            pos = nx.spectral_layout(tf_network, scale=scale)
            pos = nx.fruchterman_reingold_layout(tf_network, pos=pos, k=k/2, iterations=50)
            pos = adjust_node_positions(pos, min_dist=min_dist)
        except:
            pos = nx.spring_layout(tf_network, k=k, scale=scale)
    elif layout_type == "shell":
        pos = nx.shell_layout(tf_network, scale=scale)
    else:
        pos = nx.spring_layout(tf_network, k=k, scale=scale)

    # Hub-centricレイアウトの適用
    if hub_centric and (hubs or target_genes):
        # 中心点の計算
        center_x = sum(pos[node][0] for node in tf_network.nodes()) / len(tf_network.nodes())
        center_y = sum(pos[node][1] for node in tf_network.nodes()) / len(tf_network.nodes())
        
        # ハブノードとターゲット遺伝子を中心に配置
        central_nodes = set(hubs or set()) | set(target_genes)
        n_central = len(central_nodes)
        
        if n_central > 0:
            # 中央付近に円形に配置
            for i, node in enumerate(central_nodes):
                angle = 2 * np.pi * i / n_central
                r = 0.3 * scale  # スケールに応じた半径
                pos[node] = np.array([
                    center_x + r * np.cos(angle),
                    center_y + r * np.sin(angle)
                ])
            
            # 他のノードを外側に配置
            for node in tf_network.nodes():
                if node not in central_nodes:
                    current_pos = pos[node]
                    direction = current_pos - np.array([center_x, center_y])
                    if np.linalg.norm(direction) > 0:
                        direction = direction / np.linalg.norm(direction)
                        pos[node] = np.array([center_x, center_y]) + direction * scale
            
            # 最終的な位置調整
            pos = adjust_node_positions(pos, min_dist=min_dist)
    
    return pos

def extend_multi_target_network_with_ppi(
    tf_network: nx.DiGraph, 
    target_genes: List[str], 
    regulon_data: Dict[str, Set[str]], 
    string_data: Dict[str, Dict[str, float]],
    ppi_score_threshold: float = 0.4
) -> Tuple[nx.DiGraph, List[Tuple[str, str]], List[str]]:
    """複数のターゲット遺伝子に対してPPIネットワークを拡張"""
    extended_network = tf_network.copy()
    added_edges = []
    ppi_tfs = []
    
    for target_gene in target_genes:
        network, edges, tfs = extend_network_with_ppi(
            tf_network=extended_network,
            target_gene=target_gene,
            regulon_data=regulon_data,
            string_data=string_data,
            ppi_score_threshold=ppi_score_threshold
        )
        
        extended_network = network
        added_edges.extend(edges)
        ppi_tfs.extend(tfs)
    
    # 重複を除去
    ppi_tfs = list(set(ppi_tfs))
    
    return extended_network, added_edges, ppi_tfs



def visualize_multi_target_network(
   tf_network: nx.DiGraph,
   regulators_by_level: Dict[str, Set[str]],
   target_genes: List[str],
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
   min_node_dist: float = 0.2,
   height: int=800
) -> None:
   # ターゲット遺伝子用の色
   target_colors = [
       '#ff0000',  # Red
       '#ff6b6b',  # Light red
       '#ff9999',  # Lighter red
       '#ffb366',  # Orange
       '#ffcc00',  # Yellow
   ]
   
   # 基本的な色マップ
   color_map = {
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
   
   # ターゲット遺伝子の色を追加
   for i, target in enumerate(target_genes):
       color_idx = min(i, len(target_colors) - 1)
       color_map[f'target_{target}'] = target_colors[color_idx]

   plt.clf()
 #  fig = plt.figure(figsize=(15, 12))
   width = height * (15/12)  # 16:9のアスペクト比
   fig = plt.figure(figsize=(width/96, height/96))  # matplotlib用にインチに変換  
   pos = get_network_layout_multi(
       tf_network=tf_network,
       layout_type=layout_type,
       hub_centric=hub_centric,
       hubs=hubs,
       target_genes=target_genes,
       scale=layout_scale,
       min_dist=min_node_dist
   )
   
   # ノードの色とサイズの設定
   node_colors = []
   node_sizes = []
   
   for node in tf_network.nodes():
       size = 1000 * node_size_factor
       color = color_map['target_genes']
       
       if node in target_genes:
           color = color_map[f'target_{node}']
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

   # エッジごとに個別に描画
   if regular_edges:
       for u, v in regular_edges:
           edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
           nx.draw_networkx_edges(tf_network, pos,
                                edgelist=[(u, v)],
                                width=2.0 * edge_width_factor,
                                edge_color=[edge_color],
                                arrows=True,
                                arrowsize=20,
                                alpha=edge_alpha,
                                 min_source_margin=15,
                                 min_target_margin=20)

   if ppi_edges:
       for u, v in ppi_edges:
           nx.draw_networkx_edges(tf_network, pos,
                                edgelist=[(u, v)],
                                width=2.5 * edge_width_factor,
                                edge_color=[color_map['ppi_edge']],
                                arrows=False,
                                alpha=edge_alpha,
                                style='dashed',
                                 min_source_margin=15,
                                 min_target_margin=15)

   if regulation_edges:
       for u, v in regulation_edges:
           edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
           nx.draw_networkx_edges(tf_network, pos,
                                edgelist=[(u, v)],
                                width=2.5 * edge_width_factor,
                                edge_color=[edge_color],
                                arrows=True,
                                arrowsize=20,
                                alpha=edge_alpha,
                                min_source_margin=15,
                                min_target_margin=20)
   
   # ラベルの描画
   nx.draw_networkx_labels(tf_network, pos, font_size=font_size, 
                         font_weight='bold')
   
   # 凡例の作成
   legend_elements = []
   
   # ターゲット遺伝子の凡例
   for i, target_gene in enumerate(target_genes):
       color_idx = min(i, len(target_colors) - 1)
       legend_elements.append(
           plt.Line2D([0], [0], marker='o', color='w', 
                     label=f'Target: {target_gene}', 
                     markerfacecolor=target_colors[color_idx], markersize=15)
       )
   
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
                     linewidth=3, linestyle='--'),
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
                 label='Other Genes', 
                 markerfacecolor=color_map['target_genes'], markersize=15)
   ])

   plt.legend(handles=legend_elements, loc='center left', 
             bbox_to_anchor=(1, 0.5))
   plt.title(f"Gene Regulatory Network for Multiple Targets")
   plt.axis('off')
   st.pyplot(fig)


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

def get_edge_color(source_node, node_categories, color_scheme, target_genes, ppi_tfs=None, hubs=None):
    """エッジの色を取得する関数
    
    Args:
        source_node (str): エッジの起点となるノード
        node_categories (dict): カテゴリーごとのノードのセット
        color_scheme (dict): カテゴリーごとの色の辞書
        target_genes (list): ターゲット遺伝子のリスト
        ppi_tfs (list, optional): PPIトランスクリプション因子のリスト
        hubs (set, optional): hubノードのセット
    """
    
    # ppi_tfsがNoneの場合は空のリストに
    ppi_tfs = ppi_tfs or []
    
    # PPIトランスクリプション因子の特別な処理
    if source_node in ppi_tfs:
        return color_scheme['ppi_edge']
    
    # ターゲット遺伝子からのエッジは優先
    if source_node in target_genes:
        color_idx = f'target_{source_node}'
        return color_scheme.get(color_idx, color_scheme['target_genes'])  # 'target' を 'target_genes' に変更
    
    # hub nodeのチェック（ターゲット遺伝子以外）
    if hubs and source_node in hubs and source_node not in target_genes:
        return color_scheme['hub_node']
    
    # 他のノードカテゴリーのチェック
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

def save_static_network(tf_network: nx.DiGraph,
                       regulators_by_level: Dict[str, Set[str]],
                       target_genes: List[str],
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
                       layout_scale: float = 3.0,  # 追加
                       min_node_dist: float = 0.2,  # 追加
                       height: int = 800) -> None:
    
    # 以下は以前のコードと同じ
    if format.lower() not in ['png', 'pdf']:
        raise ValueError("Format must be either 'png' or 'pdf'")
    
    if not output_path.lower().endswith(f'.{format}'):
        output_path = f"{output_path}.{format}"
    
    plt.clf()
    #fig = plt.figure(figsize=(15, 12))
    width = height * (15/12)  # 16:9のアスペクト比
    fig = plt.figure(figsize=(width/96, height/96))  # matplotlib用にインチに変換  

    
    visualize_multi_target_network(
        tf_network=tf_network,
        regulators_by_level=regulators_by_level,
        target_genes=target_genes,
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
        min_node_dist=min_node_dist,
        height=height
    )
    
    plt.tight_layout()
    plt.savefig(output_path, bbox_inches='tight', 
                dpi=300 if format == 'png' else 600)
    plt.close(fig)


def save_interactive_network(tf_network: nx.DiGraph,
                            regulators_by_level: Dict[str, Set[str]],
                            target_genes: List[str],
                            output_path: str,
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
    target_colors = [
        '#ff0000',  # Red
        '#ff6b6b',  # Light red
        '#ff9999',  # Lighter red
        '#ffb366',  # Orange
        '#ffcc00',  # Yellow
    ]
    
    color_map = {
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

    # ターゲット遺伝子の色を追加
    for i, target in enumerate(target_genes):
        color_idx = min(i, len(target_colors) - 1)
        color_map[f'target_{target}'] = target_colors[color_idx]
    
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
                    "sortMethod": "hubsize",
                    "levelSeparation": 75,
                    "nodeSpacing": 100,
                    "treeSpacing": 100,
                    "blockShifting": True,
                    "edgeMinimization": True,
                    "parentCentralization": True,
                    "improvedLayout": True,
                    "levelBias": 0.0,
                }
            },
            "edges": {
                "smooth": {
                    "enabled": True,
                    "type": "cubicBezier",
                    "roundness": 0.5
                }
            },
            "physics": {
                "enabled": True,
                "hierarchicalRepulsion": {
                    "centralGravity": 0.0,
                    "springLength": 100,
                    "springConstant": 0.01,
                    "nodeDistance": 80,
                    "damping": 0.3
                },
                "solver": "hierarchicalRepulsion",
                "stabilization": {
                    "enabled": True,
                    "iterations": 1000,
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
        hub_nodes = set(hubs) | set(target_genes)
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
                direction = current_pos - np.array([center_x, center_y])
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
        
        if node in target_genes:
            color = color_map[f'target_{node}']
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

    # エッジの追加
    for u, v in regular_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.0 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}})

    for u, v in ppi_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        net.add_edge(u, v, width=width, color=color_map['ppi_edge'],
                    arrows={'to': {'enabled': False}},
                    dashes=True)

    for u, v in regulation_edges:
        weight = tf_network[u][v].get('weight', 1.0)
        width = (2.5 * weight * edge_width_factor)
        edge_color = get_edge_color(u, regulators_by_level, color_map, target_genes, ppi_tfs, hubs)
        net.add_edge(u, v, width=width, color=edge_color,
                    arrows={'to': {'enabled': True, 'scaleFactor': 1}})

    # JavaScriptの制御コード
    stabilization_script = '''
    <script>
    window.addEventListener('load', function() {
        let stabilizationTimeout;
        let isStabilized = false;
        let forceStopTimeout;

        function disablePhysics() {
            if (!isStabilized && typeof network !== 'undefined') {
                network.setOptions({physics: {enabled: false}});
                network.stabilize(0);
                isStabilized = true;
            }
        }

        stabilizationTimeout = setTimeout(disablePhysics, 5000);
        forceStopTimeout = setTimeout(function() {
            if (!isStabilized) {
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
                if (progress >= 80) {
                    clearTimeout(stabilizationTimeout);
                    clearTimeout(forceStopTimeout);
                    disablePhysics();
                }
            });

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

    with tempfile.NamedTemporaryFile(delete=False, suffix='.html', mode='w', encoding='utf-8') as f:
        net.save_graph(f.name)
        content = open(f.name, 'r', encoding='utf-8').read()
        content = content.replace('</head>', f'{stabilization_script}</head>')
        with open(output_path, 'w', encoding='utf-8') as out_file:
            out_file.write(content)
        os.unlink(f.name)

def add_save_buttons_to_main(st, filtered_network, regulators_by_level, target_genes, 
                             font_size, node_size, edge_width, viz_height, 
                             interactive_layout, static_layout, static_format, 
                             ppi_tfs=None, hub_centric=False, color_hubs=True, hubs=None, 
                             node_alpha=0.7, edge_alpha=0.6,
                             layout_scale=3.0,  # 追加
                             min_node_dist=0.2):  # 追加

    """Save buttons for network visualizations"""
    import tempfile
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
                    target_genes=target_genes,
                    output_path=tmp_file.name,
                    ppi_tfs=ppi_tfs,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    hub_centric=hub_centric,
                    color_hubs=color_hubs,
                    hubs=hubs,
                    format=static_format,
                    node_alpha=node_alpha,
                    edge_alpha=edge_alpha,
                    layout_type=static_layout,
                    layout_scale=layout_scale,  # 追加
                    min_node_dist=min_node_dist,
                    height=viz_height
                )
                tmp_file_path = tmp_file.name
            
            with open(tmp_file_path, "rb") as file:
                base64_file = base64.b64encode(file.read()).decode('utf-8')
            
            # ファイル名を全ターゲット遺伝子を含むものに変更
            target_names = '_'.join(target_genes)
            href = f'<a href="data:application/octet-stream;base64,{base64_file}" download="network_{target_names}.{static_format}">Download Static Network ({static_format.upper()})</a>'
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
                    target_genes=target_genes,  # 修正: target_gene -> target_genes
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
            
            # ファイル名を全ターゲット遺伝子を含むものに変更
            target_names = '_'.join(target_genes)
            href = f'<a href="data:text/html;base64,{base64_file}" download="network_{target_names}.html">Download Interactive Network (HTML)</a>'
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
    st.title("Multi-target SCENIC TF Network Analysis")
    
    # Initialize session state
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
    
    st.header("Data Upload")

    # フィルタリング設定
    use_equal_weights = st.checkbox(
        "Do not use CSI file", 
        value=False, 
        help="Ignore CSI values and treat all edges with equal weight"
    )
    
    regulon_file = st.file_uploader(
        "**Upload Regulon Data (e.g., regulons.seurat_object.txt)**", 
        type=['txt', 'tsv']
    )
    
    csi_file = False
    if not use_equal_weights:
        csi_file = st.file_uploader(
            "**Upload CSI Data (e.g., csi_results.csv)**", 
            type=['csv']
        )
        st.write("CSI: connection specificity index")
    
    if regulon_file and (csi_file or use_equal_weights):
        # データの読み込み
        regulon_data, all_genes = load_regulon_data(regulon_file)
        if not use_equal_weights:
            csi_data = load_csi_data(csi_file)
        else:
            csi_data = defaultdict(dict)
        
        # TFリストを自動的に読み込み
        all_tfs = load_tfs_list(regulon_file)
        # 生物種の判定
        species = determine_species(regulon_file)

        # Get visualization settings
        static_layout, interactive_layout, layout_scale, min_node_dist, hub_centric, color_hubs, enable_cleanup = add_visualization_settings(st)

        # 可視化パラメータ
        with st.expander("Visualization Parameters", expanded=False):
            col1, col2, col3 = st.columns(3)
            with col1:
                font_size = st.slider("Font Size", 8, 30, 12)
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
            
            static_format = st.selectbox("Static network format", ["pdf", 'png'])

        with st.form("network_settings"):     
            st.subheader("Basic Settings")
            
            # Change to multiselect for target genes
            target_genes = st.multiselect(
                "Select TFs of Interest/Target Genes",
                sorted(list(all_genes)),
                default=[sorted(list(all_genes))[0]],
                help="Select one or more target genes (max 5)",
                max_selections=5  # 最大5つまで選択可能
            )
            
            col1, col2 = st.columns(2)
            with col1:
                max_upstream = st.number_input(
                    "Max Upstream Level", 
                    min_value=0, 
                    max_value=4, 
                    value=2
                )
            with col2:
                max_downstream = st.number_input(
                    "Max Downstream Level", 
                    min_value=0, 
                    max_value=4, 
                    value=2
                )
            
            # ネットワークモード
            network_mode = st.checkbox(
                "Enable Network Mode", 
                value=False,
                help="If enabled, search both upstream and downstream connections for each TF"
            )
            
            # フィルタリング設定
            st.subheader("Filtering Settings")
            with st.expander("CSI Filtering", expanded=True):
                filtering_method = st.radio(
                    "Filtering Method",
                    options=["CSI Threshold", "Max TFs"],
                    help="Choose how to filter TFs"
                )
                
                if filtering_method == "CSI Threshold":
                    csi_threshold = st.slider("CSI Threshold", 0.0, 1.0, 0.0)
                    max_tfs = None
                else:
                    max_tfs = st.slider("Max TFs per level", 1, 50, 20)
                    csi_threshold = None

            # その他のオプション
            trrust_options = add_trrust_options(st, species)
            chip_atlas_options = add_chip_atlas_options(st, species)
            string_db_options = add_string_db_options(st, species)

            # PPI Extension
            st.subheader("PPI Extension")
            with st.expander("PPI Extension Settings", expanded=False):
                ppi_options = {
                    "enable": st.checkbox(
                        "Extend network with STRING-db PPI data", 
                        help='Add TFs that interact with target genes and regulate their targets'
                    ),
                    "score_threshold": st.slider(
                        "PPI Score Threshold", 
                        min_value=0.0, 
                        max_value=1.0, 
                        value=0.4, 
                        step=0.1
                    )
                }

            # Hub Selection
            st.subheader("Hub Selection")
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
                        "Select hub method",
                        options=list(hub_methods.keys()),
                        index=list(hub_methods.keys()).index("Common Weighted Hubs")
                    )
                    
                    top_n = st.number_input(
                        "Number of top nodes",
                        min_value=1,
                        max_value=20,
                        value=5
                    )
                
                with col2:
                    percentile_threshold = st.slider(
                        "Percentile threshold",
                        min_value=80,
                        max_value=99,
                        value=90
                    )

            submitted = st.form_submit_button("Submit")
            
            if submitted:
                st.session_state.network_analyzed = False

        if len(target_genes) == 0:
            st.warning("Please select at least one target gene.")
            return

        # Analysis and Visualization
        if st.button("Start Analysis", type="primary") or (st.session_state.network_analyzed and 'static_layout' in locals()):
            if not st.session_state.network_analyzed:
                with st.spinner('Analyzing network...'):
                    # Build multi-target network
                    tf_network, regulators_by_level = build_multi_target_network(
                        regulon_data, csi_data, all_tfs,
                        target_genes, max_upstream, max_downstream,
                        network_mode, use_equal_weights, csi_threshold, max_tfs
                    )
        
                    # Apply filters
                    filtered_network = filter_network_by_csi(tf_network, csi_threshold)

                    # PPI extension
                    ppi_tfs = []
                    if ppi_options["enable"]:
                        tf_list = list(regulon_data.keys())
                        string_data = get_string_interactions(
                            tf_list,
                            species,
                            string_db_options["score_threshold"]
                        )

                        if string_data:
                            filtered_network, ppi_edges, ppi_tfs = extend_multi_target_network_with_ppi(
                                filtered_network,
                                target_genes, 
                                regulon_data, 
                                string_data,
                                ppi_options["score_threshold"]
                            )
                    
                    # Additional filters
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

                    # Cleanup network
                    if enable_cleanup:
                        filtered_network, disconnected_edges = cleanup_network_multi(filtered_network, target_genes)
 

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
                        
                        # この部分を追加
                        hub_methods = {
                            "Top Weighted Degree Centrality": "top_degree",
                            "Top Weighted PageRank": "top_pagerank",
                            "Top Weighted Betweenness Centrality": "top_betweenness",
                            "Top Weighted Eigenvector Centrality": "top_eigenvector",
                            "Weighted Degree-based Hubs": "degree_hub",
                            "Weighted Bottleneck Nodes": "bottleneck_hub",
                            "Common Weighted Hubs": "common_hub"
                        }

                        # Identify hub nodes for visualization
                        hub_method = hub_methods[selected_hub_method]
                        centrality_data = (weighted_degree, weighted_between, weighted_pagerank, weighted_eigenvector)
                        hub_nodes_for_viz = get_hub_nodes(
                            hub_method, filtered_network, centrality_data,
                            n=top_n, percentile=percentile_threshold
                        )
                        
                        # Display identified hub nodes
                        st.write(f"\n### Identified Hub Nodes ({selected_hub_method})")
                        if hub_nodes_for_viz:
                            # Create a DataFrame for hub nodes with their centrality metrics
                            hub_data = []
                            for node in hub_nodes_for_viz:
                                hub_data.append({
                                    'Node': node,
                                    'Degree Centrality': weighted_degree.get(node, 0),
                                    'Betweenness Centrality': weighted_between.get(node, 0),
                                    'PageRank': weighted_pagerank.get(node, 0),
                                    'Eigenvector Centrality': weighted_eigenvector.get(node, 0)
                                })
                            
                            hub_df = pd.DataFrame(hub_data)
                            hub_df = hub_df.round(4)  # Round to 4 decimal places
                            st.dataframe(hub_df, hide_index=True)
                            
                            # Display hub statistics
                            st.write("\n**Hub Statistics:**")
                            col1, col2 = st.columns(2)
                            with col1:
                                st.write(f"Total number of hub nodes: {len(hub_nodes_for_viz)}")
                                st.write(f"Percentage of network: {(len(hub_nodes_for_viz) / filtered_network.number_of_nodes() * 100):.1f}%")
                            with col2:
                                # Calculate average connectivity for hub nodes
                                hub_connectivity = sum(filtered_network.degree(node) for node in hub_nodes_for_viz) / len(hub_nodes_for_viz)
                                st.write(f"Average hub connectivity: {hub_connectivity:.2f}")
                                st.write(f"Network nodes: {filtered_network.number_of_nodes()}")
                        else:
                            st.warning("No hub nodes identified with current settings.")

                        # Store results in session state
                        st.session_state.filtered_network = filtered_network
                        st.session_state.regulators_by_level = regulators_by_level
                        st.session_state.ppi_tfs = ppi_tfs
                        st.session_state.hub_nodes_for_viz = hub_nodes_for_viz
                        st.session_state.network_analyzed = True
                    else:
                        st.warning("Not enough nodes for analysis after filtering.")
                        return

            # Visualization
            if st.session_state.filtered_network.number_of_edges() > 0:
                st.subheader("Network Visualization")
                
                # Static network visualization
                visualize_multi_target_network(
                    tf_network=st.session_state.filtered_network,
                    regulators_by_level=st.session_state.regulators_by_level,
                    target_genes=target_genes,
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
                    min_node_dist=min_node_dist,
                    height=viz_height
                )
                
                # Interactive network visualization
                create_interactive_multi_target_network(
                    tf_network=st.session_state.filtered_network,
                    regulators_by_level=st.session_state.regulators_by_level,
                    target_genes=target_genes,
                    ppi_tfs=st.session_state.ppi_tfs,
                    font_size=font_size,
                    node_size_factor=node_size,
                    edge_width_factor=edge_width,
                    height=viz_height,
                    layout_type=interactive_layout,
                    hub_centric=hub_centric,
                    color_hubs=color_hubs,
                    hubs=st.session_state.hub_nodes_for_viz
                )
                
                add_save_buttons_to_main(
                    st=st,
                    filtered_network=st.session_state.filtered_network,
                    regulators_by_level=st.session_state.regulators_by_level,
                    target_genes=target_genes,
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
                    edge_alpha=edge_alpha,
                    layout_scale=layout_scale,  # 追加
                    min_node_dist=min_node_dist,
                )
            else:
                st.warning("No edges remain after filtering. Try adjusting the parameters.")
    else:
        st.info("Please upload both Regulon and CSI data files to start the analysis.")

if __name__ == "__main__":
    main()
