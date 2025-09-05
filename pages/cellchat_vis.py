import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import networkx as nx
import scipy.sparse
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import warnings
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from networkx.algorithms.centrality import current_flow_betweenness_centrality, current_flow_closeness_centrality
from matplotlib.lines import Line2D
import scanpy as sc
import re
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

def netAnalysis_signalingRole_scatter(
    net,
    signaling=None,
    color_use=None,
    slot_name="netP",
    group=None,
    weight_MinMax=None,
    dot_size=(2, 6),
    point_shape=('o', 's', '^', 'D', 'v', '<', '>', 'p', 'X'),
    label_size=3,
    dot_alpha=0.6,
    x_measure="outdeg",
    y_measure="indeg",
    xlabel="Outgoing interaction strength",
    ylabel="Incoming interaction strength",
    title=None,
    font_size=10,
    font_size_title=10,
    do_label=True,
    show_legend=True,
    show_axes=True,
    width=8,
    height=5,
    xlim=None,
    ylim=None,
    use_count=False,
    sources_use=None,
    targets_use=None,
    sorted_order=None
):
    """
    2D visualization of dominant senders (sources) and receivers (targets)
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    from matplotlib.lines import Line2D
    import traceback
    
    # 単純化のために変換関数
    def transform_size(size):
        return size**2 * 10
    
    try:
        # 中心性情報の存在チェック
        if slot_name not in net or 'centr' not in net[slot_name]:
            raise ValueError("Please run netAnalysis_computeCentrality to compute network centrality scores")
        
        centr = net[slot_name]['centr']
        
        # 細胞タイプの取得
        if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            all_cell_types = net['net']['dimnames'][0]
        else:
            raise ValueError("Cell type information not found in network data")
        
        # Apply filtering if specified - for signaling role, show selected cells in full analysis
        if sources_use is not None or targets_use is not None:
            if sources_use is not None and targets_use is not None:
                # Keep union of sources and targets
                filtered_cells = list(set(sources_use) | set(targets_use))
            elif sources_use is not None:
                # Keep only sources
                filtered_cells = sources_use
            elif targets_use is not None:
                # Keep only targets
                filtered_cells = targets_use
            else:
                filtered_cells = all_cell_types
            
            # Filter to existing cells only
            cell_types = [ct for ct in all_cell_types if ct in filtered_cells]
        else:
            cell_types = all_cell_types
        
        n_cells = len(cell_types)
        
        # Create mapping from original indices to filtered indices
        cell_type_indices = {ct: i for i, ct in enumerate(all_cell_types)}
        filtered_indices = [cell_type_indices[ct] for ct in cell_types]
        
        # 中心性値の取得 - 最初は全細胞数で初期化
        n_all_cells = len(all_cell_types)
        outgoing_cells_all = np.zeros(n_all_cells)
        incoming_cells_all = np.zeros(n_all_cells)
        
        # 集計モード
        if signaling is None:
            print("Signaling role analysis on the aggregated cell-cell communication network")
            if isinstance(centr, dict) and "aggregate" in centr:
                outgoing_cells_all = np.array(centr["aggregate"][x_measure])
                incoming_cells_all = np.array(centr["aggregate"][y_measure])
            else:
                for path_idx in centr:
                    if x_measure in centr[path_idx] and y_measure in centr[path_idx]:
                        try:
                            out_vals = np.array(centr[path_idx][x_measure])
                            in_vals = np.array(centr[path_idx][y_measure])
                            # 長さが一致しない場合の処理
                            if len(out_vals) == n_all_cells:
                                outgoing_cells_all += out_vals
                            else:
                                print(f"Warning: Pathway {path_idx} has {len(out_vals)} cells, expected {n_all_cells}")
                            if len(in_vals) == n_all_cells:
                                incoming_cells_all += in_vals
                        except Exception as e:
                            print(f"Error processing pathway {path_idx}: {e}")
            
            # フィルタリング後の細胞のみを抽出
            outgoing_cells = outgoing_cells_all[filtered_indices]
            incoming_cells = incoming_cells_all[filtered_indices]
        
        # 特定パスウェイモード
        else:
            print(f"Signaling role analysis on pathway: {signaling}")
            
            # 文字列を配列に変換
            if isinstance(signaling, str):
                signaling = [signaling]
            
            # パスウェイ情報を取得
            if 'pathways' in net[slot_name]:
                pathways = net[slot_name]['pathways']
                if isinstance(pathways, np.ndarray):
                    pathways = pathways.tolist()
                
                # マッチングパスウェイの検索
                valid_pathways = []
                for s in signaling:
                    if s in pathways:
                        valid_pathways.append(s)
                
                if valid_pathways:
                    print(f"Found valid pathways: {valid_pathways}")
                    
                    # パスウェイごとに中心性を集計
                    for path in valid_pathways:
                        path_idx = pathways.index(path)
                        print(f"Processing pathway: {path} (index: {path_idx})")
                        
                        # 中心性値へのアクセス（辞書形式）
                        if isinstance(centr, dict):
                            for key in centr:
                                if (str(key) == str(path) or 
                                    key == path or 
                                    (isinstance(key, int) and key == path_idx)):
                                    print(f"Found centrality for key: {key}")
                                    if x_measure in centr[key] and y_measure in centr[key]:
                                        # Filter centrality values to selected cells only
                                        full_outgoing = np.array(centr[key][x_measure])
                                        full_incoming = np.array(centr[key][y_measure])
                                        
                                        # 全細胞の値を初期化
                                        temp_outgoing = np.zeros(n_all_cells)
                                        temp_incoming = np.zeros(n_all_cells)
                                        
                                        # 値を設定
                                        if len(full_outgoing) == n_all_cells:
                                            temp_outgoing = full_outgoing
                                        else:
                                            print(f"Warning: Pathway {key} has {len(full_outgoing)} cells, expected {n_all_cells}")
                                        
                                        if len(full_incoming) == n_all_cells:
                                            temp_incoming = full_incoming
                                        
                                        # フィルタリング後の細胞のみを加算
                                        outgoing_cells_all[filtered_indices] += temp_outgoing[filtered_indices]
                                        incoming_cells_all[filtered_indices] += temp_incoming[filtered_indices]
                        
                        # 中心性値へのアクセス（リスト形式）
                        elif path_idx < len(centr):
                            if isinstance(centr[path_idx], dict):
                                if x_measure in centr[path_idx] and y_measure in centr[path_idx]:
                                    # Filter centrality values to selected cells only
                                    full_outgoing = np.array(centr[path_idx][x_measure])
                                    full_incoming = np.array(centr[path_idx][y_measure])
                                    
                                    # 配列長チェックとフィルタリング
                                    if len(full_outgoing) == n_all_cells:
                                        outgoing_cells += full_outgoing[filtered_indices]
                                    elif len(full_outgoing) == len(filtered_indices):
                                        outgoing_cells += full_outgoing
                                    else:
                                        print(f"Warning: Outgoing array length {len(full_outgoing)} doesn't match expected {n_all_cells} or {len(filtered_indices)}")
                                    
                                    if len(full_incoming) == n_all_cells:
                                        incoming_cells += full_incoming[filtered_indices]
                                    elif len(full_incoming) == len(filtered_indices):
                                        incoming_cells += full_incoming
                                    else:
                                        print(f"Warning: Incoming array length {len(full_incoming)} doesn't match expected {n_all_cells} or {len(filtered_indices)}")
                else:
                    print(f"Warning: No valid pathways found for: {signaling}")
                    # 空のプロットを回避するため、何らかの値を追加
                    # centeridyの全エントリを走査して見つかったものを使用
                    if isinstance(centr, dict):
                        for key in centr:
                            if x_measure in centr[key] and y_measure in centr[key]:
                                full_outgoing = np.array(centr[key][x_measure])
                                full_incoming = np.array(centr[key][y_measure])
                                
                                # 配列長チェックとフィルタリング（fallback用）
                                if len(full_outgoing) == n_all_cells:
                                    outgoing_cells += full_outgoing[filtered_indices] * 0.0001
                                elif len(full_outgoing) == len(filtered_indices):
                                    outgoing_cells += full_outgoing * 0.0001
                                
                                if len(full_incoming) == n_all_cells:
                                    incoming_cells += full_incoming[filtered_indices] * 0.0001
                                elif len(full_incoming) == len(filtered_indices):
                                    incoming_cells += full_incoming * 0.0001
            
            # 直接パスウェイ指定の場合
            else:
                for path in signaling:
                    if path in centr:
                        if x_measure in centr[path] and y_measure in centr[path]:
                            full_outgoing = np.array(centr[path][x_measure])
                            full_incoming = np.array(centr[path][y_measure])
                            
                            # 配列長チェックとフィルタリング（直接パスウェイ指定用）
                            if len(full_outgoing) == n_all_cells:
                                outgoing_cells += full_outgoing[filtered_indices]
                            elif len(full_outgoing) == len(filtered_indices):
                                outgoing_cells += full_outgoing
                            else:
                                print(f"Warning: Path {path} outgoing array length {len(full_outgoing)} doesn't match expected {n_all_cells} or {len(filtered_indices)}")
                            
                            if len(full_incoming) == n_all_cells:
                                incoming_cells += full_incoming[filtered_indices]
                            elif len(full_incoming) == len(filtered_indices):
                                incoming_cells += full_incoming
                            else:
                                print(f"Warning: Path {path} incoming array length {len(full_incoming)} doesn't match expected {n_all_cells} or {len(filtered_indices)}")
        
        # リンク数の計算
        if use_count:
            if 'network' in net and 'count_matrix' in net['network']:
                count_matrix = net['network']['count_matrix']
            elif 'net' in net and 'count' in net['net']:
                count_matrix = net['net']['count']
            else:
                prob = net['net']['prob']
                count_matrix = np.sum(prob, axis=2) > 0
        else:
            prob = net['net']['prob']
            count_matrix = np.sum(prob, axis=2) > 0
        
        # DataFrameの場合はnumpy配列に変換
        if hasattr(count_matrix, 'values'):
            count_matrix = count_matrix.values

        num_link = np.zeros(n_cells)
        for i, cell_idx in enumerate(filtered_indices):
            num_link[i] = np.sum(count_matrix[cell_idx, :]) + np.sum(count_matrix[:, cell_idx]) - count_matrix[cell_idx, cell_idx]
        
        # データフレームの作成
        df = pd.DataFrame({
            'x': outgoing_cells,
            'y': incoming_cells,
            'labels': cell_types,
            'Count': num_link
        })
        
        # Apply sorted_order if provided
        if sorted_order is not None:
            # Get valid cell types from sorted_order (those that are in the data)
            valid_sorted = [cell for cell in sorted_order if cell in df['labels'].values]
            
            # Get remaining cells (those not in sorted_order)
            remaining_cells = [cell for cell in df['labels'].values if cell not in valid_sorted]
            
            # Create new order combining valid_sorted and remaining_cells
            new_order = valid_sorted + remaining_cells
            
            # Reorder dataframe to match sorted_order
            df['sort_key'] = df['labels'].map({cell: i for i, cell in enumerate(new_order)})
            df = df.sort_values('sort_key').drop('sort_key', axis=1).reset_index(drop=True)
            
            # Update cell_types to match the new order
            cell_types = new_order
        
        print(f"Data ranges - X: {min(outgoing_cells)}-{max(outgoing_cells)}, Y: {min(incoming_cells)}-{max(incoming_cells)}")
        
        # グループ情報の追加
        if group is not None:
            if isinstance(group, dict):
                df['Group'] = df['labels'].map(lambda x: group.get(x, "Others"))
            else:
                df['Group'] = group
        
        # カラーマッピングの設定
        if color_use is None:
            cmap = plt.cm.tab10
            node_colors = [cmap(i % 10) for i in range(n_cells)]
            color_map_dict = {ct: node_colors[i] for i, ct in enumerate(cell_types)}
        elif isinstance(color_use, str):
            try:
                cmap = plt.cm.get_cmap(color_use)
                node_colors = [cmap(i / max(1, n_cells - 1)) for i in range(n_cells)]
                color_map_dict = {ct: node_colors[i] for i, ct in enumerate(cell_types)}
            except:
                cmap = plt.cm.tab10
                node_colors = [cmap(i % 10) for i in range(n_cells)]
                color_map_dict = {ct: node_colors[i] for i, ct in enumerate(cell_types)}
        elif isinstance(color_use, dict):
            color_map_dict = color_use
        else:
            color_map_dict = {ct: color_use[i % len(color_use)] for i, ct in enumerate(cell_types)}
        
        # マーカー設定（group情報がある場合）
        marker_dict = {}
        if group is not None:
            unique_groups = sorted(df['Group'].unique())
            marker_dict = {g: point_shape[i % len(point_shape)] for i, g in enumerate(unique_groups)}
        
        # ドットサイズの正規化
        sizes = df['Count'].values.astype(float)
        if weight_MinMax is not None:
            min_w, max_w = weight_MinMax
            sizes = np.clip(sizes, min_w, max_w)
        size_min, size_max = np.min(sizes), np.max(sizes)
        if size_max > size_min:
            sizes_norm = (sizes - size_min) / (size_max - size_min)
        else:
            sizes_norm = np.ones_like(sizes) * 0.5
            
        dot_size_min, dot_size_max = dot_size
        sizes_plot = dot_size_min + (dot_size_max - dot_size_min) * sizes_norm
        
        # 単一シグナルモードのデータ調整（必要に応じて）
        if signaling is not None:
            # 軸範囲を調整して中心位置を合わせる必要がある場合
            # outgoing_cells と incoming_cells が非常に小さい値の場合に調整
            
            # 0以外の最小値を取得
            nonzero_x = outgoing_cells[outgoing_cells > 0]
            nonzero_y = incoming_cells[incoming_cells > 0]
            
            if len(nonzero_x) > 0 and len(nonzero_y) > 0:
                min_nonzero_x = np.min(nonzero_x)
                min_nonzero_y = np.min(nonzero_y)
                
                # 非常に小さい値（10^-7以下など）の場合はスケーリング
                if max(outgoing_cells) < 1e-7 or max(incoming_cells) < 1e-7:
                    scale_factor = 1e5
                    outgoing_cells = outgoing_cells * scale_factor
                    incoming_cells = incoming_cells * scale_factor
                    df['x'] = outgoing_cells
                    df['y'] = incoming_cells
                    
                    print(f"Applied scaling factor {scale_factor} to pathway data")
                    print(f"New data ranges - X: {min(outgoing_cells)}-{max(outgoing_cells)}, Y: {min(incoming_cells)}-{max(incoming_cells)}")
        
        # 単一プロット形式に変更（GridSpecを使わない）
        fig, ax = plt.subplots(figsize=(width, height))
        
        # 散布図のプロット
        if group is not None:
            for g in df['Group'].unique():
                group_df = df[df['Group'] == g]
                for idx in group_df.index:
                    row = group_df.loc[idx]
                    marker = marker_dict.get(g, 'o')
                    color = color_map_dict.get(row['labels'], 'gray')
                    ax.scatter(
                        row['x'], row['y'],
                        s=sizes_plot[idx] ** 2 * 10,
                        c=[color],
                        marker=marker,
                        alpha=dot_alpha,
                        edgecolor='black',
                        linewidth=0.5
                    )
        else:
            for idx, row in df.iterrows():
                color = color_map_dict.get(row['labels'], 'gray')
                ax.scatter(
                    row['x'], row['y'],
                    s=transform_size(sizes_plot[idx]),
                    c=[color],
                    marker='o',
                    alpha=dot_alpha,
                    edgecolor='black',
                    linewidth=0.5
                )
        
        # ラベル（細胞名は黒で、少し上に配置）
        if do_label:
            for idx, row in df.iterrows():
                # 単一シグナルモードでラベル位置を調整
                y_offset = 0.4
                if signaling is not None and max(incoming_cells) < 0.2:  # 小さい値の範囲の場合
                    y_offset = max(incoming_cells) * 0.02  # より小さいオフセット
                
                ax.text(
                    row['x'], row['y'] + y_offset,  # オフセット調整
                    row['labels'],
                    fontsize=label_size,
                    ha='center', va='bottom',
                    color='black',
                    fontweight='bold',
                    bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', pad=0)
                )
        
        # タイトルと軸ラベルの設定
        if title:
            # タイトルをプロットの上部に配置（余白を増やす）
            plt.subplots_adjust(top=0.85)  # タイトル用の余白を確保
            ax.set_title(title, fontsize=font_size_title, pad=20)
        
        ax.set_xlabel(xlabel, fontsize=font_size)
        ax.set_ylabel(ylabel, fontsize=font_size)
        
        # 軸の範囲設定
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
            
        # 単一シグナルモードで軸範囲を自動調整
        if signaling is not None:
            # 0と最大値の間に適切なマージンを設ける
            x_max = max(outgoing_cells)
            y_max = max(incoming_cells)
            
            # 最小値が0よりかなり大きい場合は、最小値も考慮
            x_min = min(outgoing_cells)
            y_min = min(incoming_cells)
            
            # マージンを計算（データ範囲の5%）
            x_margin = (x_max - x_min) * 0.05
            y_margin = (y_max - y_min) * 0.05
            
            # 最小値が0からほとんど変わらない場合は0から開始
            if x_min < 0.001:
                ax.set_xlim(-x_margin, x_max + x_margin)
            else:
                ax.set_xlim(x_min - x_margin, x_max + x_margin)
                
            if y_min < 0.001:
                ax.set_ylim(-y_margin, y_max + y_margin)
            else:
                ax.set_ylim(y_min - y_margin, y_max + y_margin)
        
        # グリッドとスタイル
        ax.grid(True, linestyle='--', alpha=0.3)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        
        # 凡例の作成と配置
        if show_legend:
            # 細胞タイプ凡例要素
            cell_legend_elements = []
            for ct in cell_types:
                color = color_map_dict.get(ct, 'gray')
                cell_legend_elements.append(
                    Line2D([0], [0], marker='o', color='w', label=ct,
                         markerfacecolor=color, markersize=8, alpha=dot_alpha)
                )
            
            # グループ凡例要素
            if group is not None:
                for g in sorted(marker_dict.keys()):
                    marker = marker_dict.get(g, 'o')
                    cell_legend_elements.append(
                        Line2D([0], [0], marker=marker, color='w', label=f"Group: {g}",
                              markerfacecolor='gray', markersize=8)
                    )
            
            # サイズ凡例要素
            size_legend_elements = []
            if size_max > size_min:
                min_count = df['Count'].min()
                max_count = df['Count'].max()
                mid_count = (min_count + max_count) / 2
                
                # 凡例用のサイズ計算
                min_size = dot_size_min
                mid_size = (dot_size_min + dot_size_max)/2
                max_size = dot_size_max
                
                # サイズ凡例要素を作成
                for count, size in zip([min_count, mid_count, max_count], [min_size, mid_size, max_size]):
                    adjusted_size = np.sqrt(transform_size(size))
                    size_legend_elements.append(
                        Line2D([0], [0], marker='o', color='w', label=f"Count: {count:.1f}",
                              markerfacecolor='gray', markersize=adjusted_size, alpha=0.7)
                    )
            
            # 細胞タイプ凡例をプロットの右側に配置
            if len(cell_legend_elements) > 0:
                first_legend = ax.legend(
                    handles=cell_legend_elements,
                    loc='center left',
                    bbox_to_anchor=(1.05, 0.5),
                    fontsize=font_size-1,
                    title="Cell Types"
                )
                ax.add_artist(first_legend)
            
            # サイズ凡例をプロットの下部に配置
            if len(size_legend_elements) > 0:
                second_legend = ax.legend(
                    handles=size_legend_elements,
                    loc='upper center',
                    bbox_to_anchor=(0.5, -0.15),
                    ncol=len(size_legend_elements),
                    fontsize=font_size-1,
                    title="Circle Size",
                    frameon=True
                )
        
        # レイアウト調整
        plt.tight_layout()

        # 余白を追加して凡例が見切れないようにする
        fig.subplots_adjust(right=0.75, bottom=0.2)  # この行を追加
        
        return fig
        
    except Exception as e:
        print(f"Error in netAnalysis_signalingRole_scatter: {str(e)}")
        print(traceback.format_exc())
        
        # エラー時にはシンプルな図を返す
        fig, ax = plt.subplots(figsize=(width, height))
        ax.text(0.5, 0.5, f"Error plotting data: {str(e)}", 
               ha='center', va='center', wrap=True)
        ax.axis('off')
        return fig


def compute_flow_betweenness(G):
    """
    Calculate flow betweenness centrality by calling R's sna::flowbet function directly via rpy2
    This function ensures precise matching with R's implementation of CellChat
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import numpy2ri
        import numpy as np
        
        # Activate automatic conversion between R and NumPy arrays
        numpy2ri.activate()
        
        # Import the required R packages
        sna = importr('sna')
        
        # 元のグラフから直接隣接行列を作成 (修正部分)
        adj_matrix = nx.to_numpy_array(G, weight='weight')
        
        # Create the R matrix
        r_adj_matrix = ro.r.matrix(adj_matrix, nrow=adj_matrix.shape[0], ncol=adj_matrix.shape[1])
        
        # Call sna::flowbet directly
        flowbet_result = sna.flowbet(r_adj_matrix)
        
        # Convert the result back to a Python dictionary
        flowbet_dict = {node: float(flowbet_result[i]) for i, node in enumerate(G.nodes())}
        
        return flowbet_dict
        
    except Exception as e:
        print(f"Error calling R's flowbet function: {str(e)}")
        traceback_info = traceback.format_exc()
        print(f"Traceback: {traceback_info}")
        
        # Improved fallback that attempts to match R's behavior better
        # Instead of using NetworkX's betweenness, we create zeros to be consistent with the info centrality fallback
        st.warning("R's flowbet error")
        return {node: 0.0 for node in G.nodes()}


def compute_information_centrality(G):
    """
    Calculate information centrality by calling R's sna::infocent function directly via rpy2
    This function ensures precise matching with R's implementation
    """
    try:
        import rpy2.robjects as ro
        from rpy2.robjects.packages import importr
        from rpy2.robjects import numpy2ri
        import numpy as np
        import traceback
        
        # Activate automatic conversion between R and NumPy arrays
        numpy2ri.activate()
        
        # Import the required R packages
        sna = importr('sna')
        base = importr('base')

        # 元のグラフから直接隣接行列を作成（修正部分）
        adj_matrix = nx.to_numpy_array(G, weight='weight')
        
        # Create the R matrix
        r_adj_matrix = ro.r.matrix(adj_matrix, nrow=adj_matrix.shape[0], ncol=adj_matrix.shape[1])
        
        # R CellChatと同じパラメータを使用
        infocent_result = sna.infocent(r_adj_matrix, diag=True, rescale=True, cmode="lower")
        
        # G_u を G に変更（修正部分）
        infocent_dict = {node: float(infocent_result[i]) for i, node in enumerate(G.nodes())}
        
        return infocent_dict
            
    except Exception as e:
        print(f"Error calling R's infocent function: {str(e)}")
        traceback_info = traceback.format_exc()
        print(f"Traceback: {traceback_info}")
        
        # Fallback to zeros
        return {node: 0.0 for node in G.nodes()}

def netAnalysis_signalingRole_network(
    net, 
    signaling=None, 
    font_size=10, 
    width=10, 
    height=4, 
    color_heatmap="OrBr", 
    sorted_order=None,
    show_value=False,
    cmap_name="tab10",
    hide_color_bar=False,
    color_use=None  # 新しいパラメータ追加
):
    """
    シグナリングネットワークの役割分析をヒートマップで可視化する関数
    
    Parameters
    ----------
    net : dict
        ネットワークオブジェクト。
    signaling : str, optional
        シグナル経路名。Noneの場合は集計モード。
    font_size : int, optional
        ラベル等のフォントサイズ
    width : int, optional
        プロットの幅（インチ単位）
    height : int, optional
        プロットの高さ（インチ単位）
    color_heatmap : str, optional
        ヒートマップのカラースキーム
    sorted_order : list, optional
        セルの表示順序を指定するリスト
    show_value : bool, optional
        ヒートマップの値を表示するかどうか
    cmap_name : str, optional
        カラーマップ名
    hide_color_bar : bool, optional
        上部の細胞カラーバーを非表示にするかどうか
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        ヒートマップ描画結果の図オブジェクト
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import matplotlib.gridspec as gridspec
    from scipy.cluster.hierarchy import linkage, leaves_list
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    import matplotlib.patches as patches
    
    # --- 1. 集計モードかパスウェイ別モードかを判定 ---
    if signaling == "Aggregate":
        if 'net' not in net or 'centr' not in net['net']:
            raise ValueError("Aggregate network centrality data not found.")
        centrality = net['net']['centr']
        title_prefix = "Aggregated "
        
        # 集計モードの場合の隣接行列を取得
        adj = np.array(net['net']['prob'])
        if adj.ndim == 3:
            adj = np.sum(adj, axis=2)
    else:
        if 'netP' not in net:
            raise ValueError("Network object must contain 'netP' data for pathway analysis")
        if 'pathways' not in net['netP']:
            raise ValueError("Pathways not found in network object")
        
        pathways = net['netP']['pathways']
        if isinstance(pathways, np.ndarray):
            pathways = pathways.tolist()
            
        if signaling not in pathways:
            raise ValueError(f"Pathway '{signaling}' not found in the dataset")
            
        pathway_idx = pathways.index(signaling)
        
        if 'centr' not in net['netP']:
            raise ValueError("Centrality data not found in network object")
            
        centrality = net['netP']['centr'][pathway_idx]
        title_prefix = f"{signaling} "
        
        # パスウェイ別モードの隣接行列を取得
        adj_all = np.array(net['netP']['prob'])
        if adj_all.ndim == 3:
            adj = adj_all[:, :, pathway_idx]
        else:
            adj = adj_all
    
    # --- 2. 細胞タイプ情報の取得 ---
    if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
        cell_types = list(net['net']['dimnames'][0])
    else:
        if 'outdeg' in centrality:
            cell_types = [f"Cell{i+1}" for i in range(len(centrality['outdeg']))]
        else:
            raise ValueError("Cannot determine cell types from network object")
    
    # --- 3. 中心性指標の取得 ---
    measures = ["outdeg", "indeg", "flowbet", "info"]
    measure_names = ["Sender", "Receiver", "Mediator", "Influencer"]
    
    # データの配列を事前に準備
    data_array = np.zeros((len(measure_names), len(cell_types)), dtype=float)
    
    # 中心性指標を配列に挿入
    for i, m in enumerate(measures):
        if m in centrality and isinstance(centrality[m], (list, np.ndarray)):
            # 数値型に変換して確実にfloat型にする
            data_array[i, :] = np.array([float(val) if val is not None else 0.0 
                                        for val in centrality[m]], dtype=float)
    
    # --- 4. 表示順の処理 ---
    # 元の細胞順序を保持
    original_cell_types = cell_types.copy()
    
    # 表示順序が指定されている場合は並べ替え
    if sorted_order is not None:
        # sorted_orderのうち、cell_typesに実際に存在する要素のみを使用
        valid_sorted = [cell for cell in sorted_order if cell in cell_types]
        
        # sorted_orderに含まれていない残りの細胞タイプ
        remaining_cells = [cell for cell in cell_types if cell not in valid_sorted]
        
        # 最終的な表示順序を作成
        display_cell_types = valid_sorted + remaining_cells
        
        # インデックスマッピングを作成
        original_indices = {cell: i for i, cell in enumerate(cell_types)}
        
        # 表示データの並べ替え
        display_data = np.zeros_like(data_array)
        for i, cell in enumerate(display_cell_types):
            if cell in original_indices:
                orig_idx = original_indices[cell]
                display_data[:, i] = data_array[:, orig_idx]
    else:
        # 並べ替えがない場合は元のデータをそのまま使用
        display_cell_types = cell_types
        display_data = data_array
    
    # --- 5. 各行を0-1に正規化 ---
    for i in range(display_data.shape[0]):
        max_val = np.max(display_data[i, :])
        if max_val > 0:
            display_data[i, :] = display_data[i, :] / max_val
    
    # --- 6. ヒートマップの作成 ---
    fig = plt.figure(figsize=(width, height))
    
    # グリッドスペック設定 - カラーバーの有無に応じて調整
    if hide_color_bar:
        # カラーバーを非表示にする場合、単純なグリッドレイアウト
        gs = gridspec.GridSpec(1, 2, width_ratios=[10, 1], wspace=0.05)
        ax_heatmap = plt.subplot(gs[0])
        cbar_ax = plt.subplot(gs[1])
    else:
        # カラーバー付きレイアウト
        gs = gridspec.GridSpec(2, 2, width_ratios=[10, 1], height_ratios=[1, 10],
                              wspace=0.05, hspace=0)  # hspaceを0に設定
        ax_heatmap = plt.subplot(gs[1, 0])
        cbar_ax = plt.subplot(gs[1, 1])
        ax_colors = plt.subplot(gs[0, 0])
    
    # カラーマップ設定
    if color_heatmap == "OrBr":
        colors = [(1, 1, 1), (0.996, 0.941, 0.851), (0.992, 0.8, 0.6), 
                 (0.992, 0.6, 0.4), (0.8, 0.4, 0.2), (0.6, 0.2, 0.1)]
        cmap = mcolors.LinearSegmentedColormap.from_list("custom_OrBr", colors)
    else:
        cmap = plt.get_cmap(color_heatmap)
    
    # ヒートマップ描画
    im = ax_heatmap.imshow(display_data, cmap=cmap, aspect='auto', vmin=0, vmax=1)
 
    if show_value: 
        # 数値をセルに表示
        for i in range(display_data.shape[0]):
            for j in range(display_data.shape[1]):
                # 小数点3桁で表示
                text = ax_heatmap.text(j, i, f"{display_data[i, j]:.3f}",
                                   ha="center", va="center", color="black",
                                   fontsize=font_size)
    
    # 軸ラベルの設定
    ax_heatmap.set_xticks(np.arange(len(display_cell_types)))
    ax_heatmap.set_yticks(np.arange(len(measure_names)))
    ax_heatmap.set_xticklabels(display_cell_types, rotation=90, fontsize=font_size)
    ax_heatmap.set_yticklabels(measure_names, fontsize=font_size)
    
    # 境界線を追加してセルを区分
    for edge, spine in ax_heatmap.spines.items():
        spine.set_visible(True)
        spine.set_color('black')
        spine.set_linewidth(0.5)
    
    # カラーバーの追加
    cbar = plt.colorbar(im, cax=cbar_ax, shrink=0.5)
    cbar.set_label('Importance', fontsize=font_size)
    cbar.ax.tick_params(labelsize=font_size-2)
    
    # 細胞タイプのカラーバーを追加（非表示オプションでなければ）

    # 細胞タイプのカラーバー部分の修正
    if not hide_color_bar:
        # 色の割り当て
        if color_use is not None and isinstance(color_use, dict):
            # 提供された色マッピングを使用
            original_colors = {}
            for cell in original_cell_types:
                original_colors[cell] = color_use.get(cell, (0.7, 0.7, 0.7, 1.0))
        else:
            # 従来通りカラーマップから生成
            try:
                cmap = cm.get_cmap(cmap_name)
            except:
                cmap = cm.get_cmap('tab10')
            
            original_colors = {}
            for i, cell in enumerate(original_cell_types):
                idx = i % 20 if cmap.N > 20 else i % cmap.N
                color = cmap(idx / max(1, min(20, cmap.N) - 1))
                original_colors[cell] = color
        
        # 表示順序に従って色のリストを作成（各細胞は元の色を維持）
        cell_colors = []
        for cell in display_cell_types:
            if cell in original_colors:
                cell_colors.append(original_colors[cell])
            else:
                # 元のリストにない細胞（通常ないはず）
                cell_colors.append((0.7, 0.7, 0.7, 1.0))
            
        # カラーバーをImageで描画（これが最も正確に余白なしで表示できる）
        color_array = np.array([cell_colors])
        ax_colors.imshow(color_array, aspect='auto', interpolation='nearest')
        ax_colors.set_xticks([])
        ax_colors.set_yticks([])
        ax_colors.axis('off')
    
    # タイトル
    plt.suptitle(f"{title_prefix}Signaling Role Analysis", fontsize=font_size+4, y=0.98)
    
    plt.tight_layout()
    plt.subplots_adjust(top=0.9, bottom=0.1)
    
    return fig

# R版のcolorRamp3関数をPythonで実装
def colorRamp3(breaks, colors, alpha=1.0):
    """
    RのcolorRamp3関数に相当する機能を実装
    
    Parameters
    ----------
    breaks : list or np.ndarray
        色の変化点となる値のリスト
    colors : list
        breaksに対応する色のリスト
    alpha : float, optional
        透明度（0.0〜1.0）
    
    Returns
    -------
    function
        数値を受け取り色を返す関数
    """
    import numpy as np
    import matplotlib.colors as mcolors
    
    breaks = np.array(breaks)
    
    # 色をRGBA形式に変換
    rgba_colors = []
    for color in colors:
        # 色が既にRGBAタプルの場合
        if isinstance(color, tuple) and len(color) in [3, 4]:
            if len(color) == 3:
                rgba_colors.append(color + (alpha,))
            else:
                rgba_colors.append(color)
        # 色が文字列の場合
        else:
            rgba = mcolors.to_rgba(color, alpha=alpha)
            rgba_colors.append(rgba)
    
    # 色を返す関数
    def color_mapper(x):
        x = np.asarray(x)
        
        # スカラー値の場合
        if np.isscalar(x) or x.size == 1:
            if np.isnan(x):
                return (1, 1, 1, 0)  # NaNは透明に
                
            # 最小値未満または最大値超過の処理
            if x <= breaks[0]:
                return rgba_colors[0]
            elif x >= breaks[-1]:
                return rgba_colors[-1]
                
            # 適切な区間を見つける
            for i in range(len(breaks) - 1):
                if breaks[i] <= x <= breaks[i + 1]:
                    # 線形補間
                    t = (x - breaks[i]) / (breaks[i + 1] - breaks[i])
                    r = rgba_colors[i][0] * (1 - t) + rgba_colors[i + 1][0] * t
                    g = rgba_colors[i][1] * (1 - t) + rgba_colors[i + 1][1] * t
                    b = rgba_colors[i][2] * (1 - t) + rgba_colors[i + 1][2] * t
                    a = rgba_colors[i][3] * (1 - t) + rgba_colors[i + 1][3] * t
                    return (r, g, b, a)
        
        # 配列の場合は要素ごとに処理
        else:
            result = []
            for val in x.flatten():
                result.append(color_mapper(val))
            return np.array(result).reshape(x.shape + (4,))
    
    return color_mapper

def netVisual_circle_individual(
    net, 
    vertex_weight=None, 
    weight_scale=True, 
    edge_weight_max=8, 
    edge_width_max=6,
    arrow=True, 
    vertex_size_max=5, 
    ncols=4, 
    measure="weight",
    pathway_name=None,
    title_name=None,
    cmap_name="tab10",
    sorted_order=None,
    alpha_edge=0.6,
    color_use=None  # 細胞名と色のマッピング辞書を追加
):
    """
    各細胞タイプごとのサークルプロットを複数サブプロットで描画する。
    ノード色は color_use によるマッピングまたは cmap_name から割り当て、
    エッジは始点ノード色の半透明で描画する。

    Parameters
    ----------
    net : dict or DataFrame or np.ndarray
        Interaction matrix or dict containing interaction matrices.
        For a dict, if 'net' key is present, it should contain either a 'weight' or 'count' key.
    vertex_weight : list or np.ndarray, optional
        Weights for each vertex/node (e.g., cell counts) used to scale node sizes.
    weight_scale : bool, optional
        Whether to scale edge width by weights (default: True).
    edge_weight_max : float, optional
        Maximum edge weight for scaling. If None, computed automatically.
    arrow : bool, optional
        Whether to add arrowheads to edges.
    vertex_size_max : float, optional
        Maximum vertex size.
    ncols : int, optional
        Number of columns in the subplot grid.
    measure : str, optional
        "weight" or "count". Determines which interaction matrix to use.
    pathway_name : str または list, optional
        特定のパスウェイ名または複数のパスウェイ名のリスト。None または "Aggregate" の場合は集計ネットワークを使用。
    title_name : str or None
        Overall title for the entire figure (optional).
    cmap_name : str
        matplotlib colormap name (default "tab10").
    sorted_order : list or None
        細胞タイプの並び順。指定があれば、cell_types とネットワーク行列を並び替える。
    alpha_edge : float, optional
        エッジの透明度（default: 0.6）。
    color_use : dict, optional
        細胞名をキー、色を値とする辞書。指定があればそれを優先して使用する。

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing individual circle plots.
    """
    import numpy as np
    import matplotlib.pyplot as plt
    import networkx as nx
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import math

    # パスウェイ名を常にリストとして扱う
    if pathway_name is not None and not isinstance(pathway_name, list):
        pathway_name = [pathway_name]

    net_array = None
    cell_types = None

    # パスウェイ指定あり、かつ "Aggregate" でない場合
    if pathway_name is not None and pathway_name[0] != "Aggregate" and pathway_name[0] is not None:
        # 複数パスウェイの場合の処理
        if isinstance(net, dict) and 'netP' in net and 'pathways' in net['netP'] and 'prob' in net['netP']:
            pathways = net['netP']['pathways']
            if isinstance(pathways, np.ndarray):
                pathways = pathways.tolist()
            
            # パスウェイインデックスを取得
            pathway_indices = []
            valid_pathways = []
            
            for p in pathway_name:
                try:
                    idx = pathways.index(p)
                    pathway_indices.append(idx)
                    valid_pathways.append(p)
                except ValueError:
                    print(f"Warning: Pathway '{p}' not found in dataset")
            
            if not pathway_indices:
                raise ValueError(f"None of the specified pathways were found.")
            
            # 複数パスウェイの確率行列を集計
            prob_arrays = []
            for idx in pathway_indices:
                prob_arrays.append(net['netP']['prob'][:, :, idx])
            
            # 確率行列を合計
            net_array = np.sum(prob_arrays, axis=0)
            
            # セルタイプの取得
            if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                cell_types = net['net']['dimnames'][0]
            
            # タイトル設定（指定がない場合）
            if title_name is None:
                if len(valid_pathways) == 1:
                    title_name = f"Cell-Cell Interaction: {valid_pathways[0]}"
                else:
                    if len(valid_pathways) <= 3:
                        title_name = f"Combined Pathways: {', '.join(valid_pathways)}"
                    else:
                        title_name = f"Combined {len(valid_pathways)} Pathways"

    # 1. net が dict なら measure に応じた行列を取得
    if net_array is None and isinstance(net, dict) and 'net' in net:
        if measure == "weight" and 'weight' in net['net']:
            net_array = net['net']['weight']
            if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                cell_types = net['net']['dimnames'][0]
        elif measure == "count" and 'count' in net['net']:
            net_array = net['net']['count']
            if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                cell_types = net['net']['dimnames'][0]

    if net_array is None and isinstance(net, dict) and 'network' in net:
        if measure == "weight" and 'strength_matrix' in net['network']:
            net_array = net['network']['strength_matrix']
            if hasattr(net_array, 'index'):
                cell_types = list(net_array.index)
        elif measure == "count" and 'count_matrix' in net['network']:
            net_array = net['network']['count_matrix']
            if hasattr(net_array, 'index'):
                cell_types = list(net_array.index)

    if net_array is None and isinstance(net, dict):
        if measure == "weight" and 'weight' in net:
            net_array = net['weight']
            if hasattr(net_array, 'index'):
                cell_types = list(net_array.index)
        elif measure == "count" and 'count' in net:
            net_array = net['count']
            if hasattr(net_array, 'index'):
                cell_types = list(net_array.index)

    if net_array is None and not isinstance(net, dict):
        net_array = net

    if net_array is None:
        raise ValueError("ネットワークマトリックスデータが見つかりません。")

    # 2. cell_types 未定義なら自動生成
    if cell_types is None:
        if hasattr(net_array, 'index'):  # DataFrame の場合
            cell_types = list(net_array.index)
        else:
            cell_types = [f"Cell{i+1}" for i in range(net_array.shape[0])]

    if not isinstance(cell_types, list):
        cell_types = list(cell_types)

    # DataFrame → values
    if hasattr(net_array, 'values'):
        net_array = net_array.values

    net_array = np.array(net_array)

    # --- 並び替え (sorted_order) があれば適用 ---
    if sorted_order is not None:
        reorder_idx = []
        for cell in sorted_order:
            if cell in cell_types:
                reorder_idx.append(cell_types.index(cell))
        if reorder_idx:
            cell_types = [cell_types[i] for i in reorder_idx]
            net_array = net_array[np.ix_(reorder_idx, reorder_idx)]

    n_cells = len(cell_types)

    # 3. edge_weight_max 自動計算
    if edge_weight_max is None:
        edge_weight_max = np.max(net_array)
    
    # 4. ノードの色設定 (color_use が指定されていればそれを優先)
    node_colors_dict = {}
    if color_use is not None and isinstance(color_use, dict):
        node_colors_dict = color_use.copy()
    else:
        cmap = cm.get_cmap(cmap_name)
        for i, cell in enumerate(cell_types):
            node_colors_dict[cell] = cmap(i / max(1, n_cells - 1))
    
    # 5. subplot レイアウト
    nrows = int(np.ceil(n_cells / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 4*nrows))
    if nrows * ncols == 1:
        axes = np.array([axes])
    else:
        axes = axes.flatten()

    # 6. 各サブプロットごとに描画
    for i, cell_type in enumerate(cell_types):
        if i < len(axes):
            # 各セルタイプについて、i 行のみ取り出す
            individual_net = np.zeros_like(net_array)
            individual_net[i, :] = net_array[i, :]

            # NetworkX グラフ作成
            G = nx.DiGraph()
            for j, ct in enumerate(cell_types):
                G.add_node(ct)
            for j, source in enumerate(cell_types):
                if j == i:  # i行のみ
                    for k, target in enumerate(cell_types):
                        w = individual_net[j, k]
                        if w > 0:
                            G.add_edge(source, target, weight=w)
            
            # ノード配置: 0時（12時）から時計回りに配置
            pos = {}
            radius = 5.0
            for k, ct in enumerate(cell_types):
                angle = math.pi/2 - (2 * math.pi * k / n_cells)
                x = radius * math.cos(angle)
                y = radius * math.sin(angle)
                pos[ct] = (x, y)

            # ノードサイズ
            if vertex_weight is not None and len(vertex_weight) == len(cell_types):
                vw = np.array(vertex_weight)
                vmax = vw.max() if vw.max() > 0 else 1
                node_sizes = vw / vmax * vertex_size_max * 100 + 100
            else:
                node_sizes = [300] * n_cells

            # サブプロット用のノード色リストを生成
            node_colors_plot = []
            for ct in cell_types:
                if ct in node_colors_dict:
                    node_colors_plot.append(node_colors_dict[ct])
                else:
                    node_colors_plot.append((0.7, 0.7, 0.7, 1.0))
                    print(f"Warning: No color defined for cell type '{ct}', using default gray")

            # ノード描画
            nx.draw_networkx_nodes(
                G, pos,
                node_color=node_colors_plot,
                node_size=node_sizes,
                alpha=0.9,
                ax=axes[i]
            )
            nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=axes[i])

            # エッジ描画 (始点ノード色の半透明)
            edges = list(G.edges(data=True))
            edge_colors = []
            widths = []
            for (u, v, data) in edges:
                w = data['weight']
                if weight_scale and edge_weight_max > 0:
                    w_scaled = 0.3 + (w / edge_weight_max) * edge_width_max
                else:
                    w_scaled = 0.3 + edge_width_max * w
                widths.append(w_scaled)

                # 始点 u の色（color_use あるいは cmap から取得）を半透明に
                if u in node_colors_dict:
                    c = node_colors_dict[u]
                else:
                    c = (0.7, 0.7, 0.7, 1.0)
                if len(c) == 4:
                    edge_colors.append((c[0], c[1], c[2], alpha_edge))
                else:
                    c_rgba = mcolors.to_rgba(c)
                    edge_colors.append((c_rgba[0], c_rgba[1], c_rgba[2], alpha_edge))

            if edges:
                if arrow:
                    nx.draw_networkx_edges(
                        G, pos,
                        edgelist=edges,
                        width=widths,
                        edge_color=edge_colors,
                        arrows=True,
                        arrowstyle='-|>',
                        arrowsize=15,
                        connectionstyle='arc3,rad=0.1',
                        ax=axes[i]
                    )
                else:
                    nx.draw_networkx_edges(
                        G, pos,
                        edgelist=edges,
                        width=widths,
                        edge_color=edge_colors,
                        arrows=False,
                        connectionstyle='arc3,rad=0.1',
                        ax=axes[i]
                    )
            
            axes[i].set_title(cell_type, fontsize=12)
            axes[i].axis('off')

    # unused subplot を削除
    for i in range(n_cells, len(axes)):
        fig.delaxes(axes[i])

    if title_name:
        fig.suptitle(title_name, fontsize=14)

    plt.tight_layout()
    return fig



def netVisual_circle(
    net, 
    measure="weight",
    pathway_name=None,
    vertex_weight=None, 
    weight_scale=True, 
    edge_weight_max=None, 
    arrow=True, 
    vertex_size_max=5,
    sorted_order=None,
    title_name=None,
    cmap_name="tab10",
    alpha_edge=0.6,
    edge_width_max=8,
    color_use=None,  # 細胞名と色のマッピング辞書
    sources_use=None,  # 新しいパラメータ: 使用するソース細胞タイプのリスト
    targets_use=None  # 新しいパラメータ: 使用するターゲット細胞タイプのリスト
):
    """
    1枚の円形ネットワーク図を生成する関数。
    """
    import networkx as nx
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import math


    # パスウェイ名を常にリストとして扱う
    if pathway_name is not None and not isinstance(pathway_name, list):
        pathway_name = [pathway_name]

    net_array = None
    cell_types = None

    # 1) net が dict の場合、measure に応じた行列を探す
    if isinstance(net, dict):
        # パスウェイ指定あり、かつ "Aggregate" でない場合
        if pathway_name is not None and pathway_name[0] != "Aggregate" and pathway_name[0] is not None:
            # 複数パスウェイの場合の処理
            if 'netP' in net and 'pathways' in net['netP'] and 'prob' in net['netP']:
                pathways = net['netP']['pathways']
                if isinstance(pathways, np.ndarray):
                    pathways = pathways.tolist()
                
                # パスウェイインデックスを取得
                pathway_indices = []
                valid_pathways = []
                
                for p in pathway_name:
                    try:
                        idx = pathways.index(p)
                        pathway_indices.append(idx)
                        valid_pathways.append(p)
                    except ValueError:
                        print(f"Warning: Pathway '{p}' not found in dataset")
                
                if not pathway_indices:
                    raise ValueError(f"None of the specified pathways were found.")
                
                # 複数パスウェイの確率行列を集計
                prob_arrays = []
                for idx in pathway_indices:
                    prob_arrays.append(net['netP']['prob'][:, :, idx])
                
                # 確率行列を合計
                net_array = np.sum(prob_arrays, axis=0)
                
                # セルタイプの取得
                if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                    cell_types = net['net']['dimnames'][0]
                
                # タイトル設定（指定がない場合）
                if title_name is None:
                    if len(valid_pathways) == 1:
                        title_name = f"Cell-Cell Interaction: {valid_pathways[0]}"
                    else:
                        if len(valid_pathways) <= 3:
                            title_name = f"Combined Pathways: {', '.join(valid_pathways)}"
                        else:
                            # 長すぎる場合は省略
                            title_name = f"Combined {len(valid_pathways)} Pathways"
        
        # Aggregate または パスウェイ指定なしの場合の処理
        if net_array is None:
            # 1-a) net['net'] に weight or count がある場合
            if 'net' in net:
                if measure == "weight" and 'weight' in net['net']:
                    net_array = net['net']['weight']
                    if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                        cell_types = net['net']['dimnames'][0]
                elif measure == "count" and 'count' in net['net']:
                    net_array = net['net']['count']
                    if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                        cell_types = net['net']['dimnames'][0]

            # 1-b) aggregated (net['network']) に strength_matrix / count_matrix がある場合
            if net_array is None and 'network' in net:
                if measure == "weight" and 'strength_matrix' in net['network']:
                    net_array = net['network']['strength_matrix']
                    if hasattr(net_array, 'index'):
                        cell_types = list(net_array.index)
                elif measure == "count" and 'count_matrix' in net['network']:
                    net_array = net['network']['count_matrix']
                    if hasattr(net_array, 'index'):
                        cell_types = list(net_array.index)

            # 1-c) net 直下に weight / count がある場合
            if net_array is None:
                if measure == "weight" and 'weight' in net:
                    net_array = net['weight']
                    if hasattr(net_array, 'index'):
                        cell_types = list(net_array.index)
                elif measure == "count" and 'count' in net:
                    net_array = net['count']
                    if hasattr(net_array, 'index'):
                        cell_types = list(net_array.index)

    # 2) net が dict でない場合、net 自体を行列とみなす
    if net_array is None and not isinstance(net, dict):
        net_array = net

    if net_array is None:
        raise ValueError("ネットワークマトリックスデータが見つかりません。")

    # 3) cell_types が未定義なら自動生成
    if cell_types is None:
        if hasattr(net_array, 'index'):  # DataFrame の場合
            cell_types = list(net_array.index)
        else:
            cell_types = [f"Cell{i+1}" for i in range(net_array.shape[0])]

    # list 化
    if not isinstance(cell_types, list):
        cell_types = list(cell_types)

    # DataFrame なら values に変換
    if hasattr(net_array, 'values'):
        net_array = net_array.values

    # numpy 配列化
    net_array = np.array(net_array)
    
    # Filter by sources_use and targets_use if specified
    if sources_use is not None or targets_use is not None:
        # Create a DataFrame for easier filtering
        df_matrix = pd.DataFrame(net_array, index=cell_types, columns=cell_types)
        
        # Determine which cell types to keep based on filtering
        if sources_use is not None and targets_use is not None:
            # Keep only the specified sources and targets
            all_cells = list(set(sources_use) | set(targets_use))
            cells_to_keep = [c for c in cell_types if c in all_cells]
        elif sources_use is not None:
            # Keep only specified sources (for all targets)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        elif targets_use is not None:
            # Keep only specified targets (for all sources)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        else:
            cells_to_keep = cell_types
        
        # Filter the DataFrame and ensure consistency
        if sources_use is not None:
            # Set non-source rows to 0
            for idx in df_matrix.index:
                if idx not in sources_use:
                    df_matrix.loc[idx, :] = 0
        
        if targets_use is not None:
            # Set non-target columns to 0
            for col in df_matrix.columns:
                if col not in targets_use:
                    df_matrix.loc[:, col] = 0
        
        # Keep only cells that have connections after filtering
        if sources_use is not None and targets_use is not None:
            # When both are specified, keep only the union of sources and targets
            cells_with_connections = cells_to_keep
        else:
            # Find cells that still have connections
            row_sums = df_matrix.sum(axis=1)
            col_sums = df_matrix.sum(axis=0)
            cells_with_connections = list(set(
                list(df_matrix.index[row_sums > 0]) + 
                list(df_matrix.columns[col_sums > 0])
            ))
            # Preserve original order
            cells_with_connections = [c for c in cell_types if c in cells_with_connections]
        
        if cells_with_connections:
            df_matrix = df_matrix.loc[cells_with_connections, cells_with_connections]
            net_array = df_matrix.values
            cell_types = cells_with_connections
    
    # ===== 色のマッピング設定 =====
    # 提供された色マッピングと固定の色マッピングをマージ
    node_colors_dict = {}
    
    # 提供された色マッピングがあればそれを優先
    if color_use is not None and isinstance(color_use, dict):
        for cell, color in color_use.items():
            node_colors_dict[cell] = color
    
    # ユーザーが新しいcmap_nameを指定した場合のフォールバック
    if cmap_name != "tab10" and color_use is None:
        try:
            cmap = cm.get_cmap(cmap_name)
            for i, cell in enumerate(set(cell_types)):
                if cell not in node_colors_dict:
                    idx = i % min(20, getattr(cmap, 'N', 10))
                    node_colors_dict[cell] = cmap(idx / max(1, min(20, getattr(cmap, 'N', 10)) - 1))
        except Exception as e:
            print(f"Warning: Failed to use colormap {cmap_name}: {e}")

    # 4) 孤立ノードの除外
    total_weights = [net_array[i, :].sum() + net_array[:, i].sum() for i in range(len(cell_types))]
    connected_indices = [i for i, w in enumerate(total_weights) if w > 0]
    if not connected_indices:
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.text(0.5, 0.5, "該当する細胞間相互作用がありません", ha='center', va='center', fontsize=14)
        ax.axis('off')
        return fig

    # フィルタリング
    cell_types = [cell_types[i] for i in connected_indices]
    net_array = net_array[np.ix_(connected_indices, connected_indices)]

    # 5) sorted_order があれば並び替え
    if sorted_order is not None:
        reorder_idx = []
        for cell in sorted_order:
            if cell in cell_types:
                reorder_idx.append(cell_types.index(cell))
        if reorder_idx:
            cell_types = [cell_types[i] for i in reorder_idx]
            net_array = net_array[np.ix_(reorder_idx, reorder_idx)]

    # 6) グラフ構築
    G = nx.DiGraph()
    for i, cell in enumerate(cell_types):
        G.add_node(cell)
    for i, source in enumerate(cell_types):
        for j, target in enumerate(cell_types):
            w = net_array[i, j]
            if w > 0:
                G.add_edge(source, target, weight=w)

    # 7) ===== 重要な変更: 明示的に0時から時計回りに位置を設定 =====
    pos = {}
    radius = 5.0
    n_cells = len(cell_types)
    
    # ノードを0時から時計回りに配置
    for i, cell in enumerate(cell_types):
        # 角度計算: 3π/2（270度、0時位置）から時計回りに増加
        angle = math.pi/2 - (2 * math.pi * i / n_cells)  # ここを修正
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        pos[cell] = (x, y)

    # 8) edge_weight_max が指定されていない場合、最大値を自動取得
    if edge_weight_max is None:
        edge_weight_max = net_array.max() if net_array.size > 0 else 1

    # 9) カラー設定 (ノード)
    node_colors = []
    for cell in cell_types:
        if cell in node_colors_dict:
            node_colors.append(node_colors_dict[cell])
        else:
            # デフォルト色はグレー
            node_colors.append((0.7, 0.7, 0.7, 1.0))
            # デバッグ情報
            print(f"Warning: No color defined for cell type '{cell}', using default gray")

    # 10) ノードサイズ計算
    if vertex_weight is not None and len(vertex_weight) == len(cell_types):
        vw = np.array(vertex_weight)
        vmax = vw.max() if vw.max() > 0 else 1
        node_sizes = vw / vmax * vertex_size_max * 200 + 200
    else:
        node_sizes = [300] * len(cell_types)

    # 11) 描画
    fig, ax = plt.subplots(figsize=(8, 8))
    # ノード描画
    nx.draw_networkx_nodes(
        G, pos,
        node_color=node_colors,
        node_size=node_sizes,
        alpha=0.9,
        ax=ax
    )
    # ラベル描画
    nx.draw_networkx_labels(G, pos, font_size=8, font_weight='bold', ax=ax)

    # 12) エッジ描画
    edges = list(G.edges(data=True))
    edge_colors = []
    widths = []
    if edges:
        for (u, v, data) in edges:
            w = data['weight']
            if weight_scale and edge_weight_max > 0:
                w_scaled = 0.3 + (w / edge_weight_max) * edge_width_max
            else:
                w_scaled = 0.3 + edge_width_max * w
            widths.append(w_scaled)
            
            # エッジの色 = 始点ノード u の色
            if u in node_colors_dict:
                c = node_colors_dict[u]
            else:
                c = (0.7, 0.7, 0.7, 1.0)  # デフォルト色
                
            if len(c) == 4:
                edge_colors.append((c[0], c[1], c[2], alpha_edge))
            else:
                c_rgba = mcolors.to_rgba(c)
                edge_colors.append((c_rgba[0], c_rgba[1], c_rgba[2], alpha_edge))

    if arrow:
        nx.draw_networkx_edges(
            G, pos,
            edgelist=edges,
            width=widths,
            edge_color=edge_colors,
            arrows=True,
            arrowstyle='-|>',
            arrowsize=15,
            connectionstyle='arc3,rad=0.1',
            ax=ax
        )
    else:
        nx.draw_networkx_edges(
            G, pos,
            edgelist=edges,
            width=widths,
            edge_color=edge_colors,
            arrows=False,
            connectionstyle='arc3,rad=0.1',
            ax=ax
        )

    if title_name:
        ax.set_title(title_name, fontsize=12, pad=20)
    ax.axis('off')
    ax.set_aspect('equal')
    fig.tight_layout()
    return fig

def netVisual_aggregate(object, signaling, signaling_name=None, color_use=None, thresh=0.05, vertex_receiver=None, sources_use=None, targets_use=None, idents_use=None, top=1, remove_isolate=False,
                           vertex_weight=1, vertex_weight_max=None, vertex_size_max=None,
                           weight_scale=True, edge_weight_max=None, edge_width_max=8,
                           layout="circle",
                           pt_title=12, title_space=6, vertex_label_cex=0.8,
                           sample_use=None, alpha_image=0.15, point_size=1.5,
                           group=None, cell_order=None, small_gap=1, big_gap=10, scale=False, reduce=-1, show_legend=False, legend_pos_x=20, legend_pos_y=20):
    """
    集計されたネットワークをグラフ形式で描画
    
    Parameters
    ----------
    object : dict
        CellChat結果のdictionary
    signaling : str
        シグナリングパスウェイ名
    signaling_name : str, optional
        代替シグナルパスウェイ名（表示用）
    color_use : dict or list, optional
        各細胞グループの色
    ... (他のパラメータ) ...
    
    Returns
    -------
    gg : matplotlib.figure.Figure or plotly.graph_objects.Figure
        プロット図
    """
    layout = layout.lower()  # 小文字に統一して処理を簡略化
    
    if signaling_name is None:
        signaling_name = signaling
    
    # シグナリングパスウェイ情報を検索
    pairLR = None  # 必要に応じて実装
    
    # ネットワーク情報を取得
    net = object.get('net', {})
    netP = object.get('netP', {})
    
    # 確率とp値を取得
    prob = net.get('prob', np.array([]))
    pval = net.get('pval', np.array([]))
    
    if prob.size == 0 or pval.size == 0:
        raise ValueError("No probability or p-value data available")
    
    # 閾値でフィルタリング
    prob_filtered = prob.copy()
    prob_filtered[pval > thresh] = 0
    
    # シグナリングまたはL-Rペアのフィルタリング（必要に応じて）
    
    # レイアウトによって処理分岐
    if layout == "hierarchy":
        # ヒエラルキーレイアウトの処理（必要に応じて実装）
        gg = None  # 実装が必要
        
    elif layout == "circle":
        # サークルレイアウトの処理
        prob_sum = np.sum(prob_filtered, axis=2)
        if edge_weight_max is None:
            edge_weight_max = np.max(prob_sum)

        gg = netVisual_circle(prob_sum, sources_use=sources_use, targets_use=targets_use, 
                             idents_use=idents_use, remove_isolate=remove_isolate, 
                             top=top, color_use=color_use, 
                             vertex_weight=vertex_weight, vertex_weight_max=vertex_weight_max, 
                             vertex_size_max=vertex_size_max, weight_scale=weight_scale, 
                             edge_weight_max=edge_weight_max, edge_width_max=8,
                             title_name=f"{signaling_name} signaling pathway network", 
                             vertex_label_cex=vertex_label_cex)
    
    elif layout == "chord":
        # コードダイアグラムの処理（必要に応じて実装）
        gg = None  # 実装が必要
        
    elif layout == "spatial":
        # 空間プロットの処理（必要に応じて実装）
        gg = None  # 実装が必要
    
    else:
        raise ValueError(f"Unsupported layout: {layout}")
    
    return gg


def netVisual_chord(net, group_colors=None, title="Cell-Cell Communication Network", figsize=(10, 10), sources_use=None, targets_use=None):
    """
    Create a chord diagram for the cell-cell communication network using matplotlib
    
    Parameters
    ----------
    net : dict
        Network object containing communication data
    group_colors : dict or list, optional
        Colors for each cell group
    title : str, optional
        Title for the plot
    figsize : tuple, optional
        Figure size (width, height) in inches
    sources_use : list or None
        List of source cell types to use. If None, use all
    targets_use : list or None
        List of target cell types to use. If None, use all
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import networkx as nx
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import math
    
    # Extract the appropriate matrix and cell types
    weight_matrix = None
    cell_types = None
    
    # Check various possible data structures
    # 1. Check for 'net' with 'prob'
    if isinstance(net, dict) and 'net' in net:
        if 'prob' in net['net']:
            # Sum across all interactions if it's a 3D array
            prob = net['net']['prob']
            if prob.ndim == 3:
                weight_matrix = np.sum(prob, axis=2)
            else:
                weight_matrix = prob
            
            # Get cell types if available
            if 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
                cell_types = net['net']['dimnames'][0]
    
    # 2. Check for 'network' with 'strength_matrix'
    if weight_matrix is None and isinstance(net, dict) and 'network' in net:
        if 'strength_matrix' in net['network']:
            weight_matrix = net['network']['strength_matrix']
            # Try to get cell types
            if hasattr(weight_matrix, 'index'):
                cell_types = weight_matrix.index
    
    # 3. Check if net itself has these keys directly
    if weight_matrix is None and isinstance(net, dict):
        for key in ['strength_matrix', 'weight', 'count']:
            if key in net:
                weight_matrix = net[key]
                if hasattr(weight_matrix, 'index'):
                    cell_types = weight_matrix.index
                break
    
    # 4. If net is the matrix directly
    if weight_matrix is None and not isinstance(net, dict):
        weight_matrix = net
    
    # If still no matrix found, raise error
    if weight_matrix is None:
        raise ValueError("Weight matrix not found in network object")
    
    # Generate cell types if not found
    if cell_types is None:
        if hasattr(weight_matrix, 'index'):  # It's a DataFrame
            cell_types = weight_matrix.index
        else:
            # Create generic names
            cell_types = [f"Cell{i+1}" for i in range(weight_matrix.shape[0])]
    
    # Convert to list if not already
    if not isinstance(cell_types, list):
        cell_types = list(cell_types)
    
    # Filter by sources_use and targets_use if specified
    if sources_use is not None or targets_use is not None:
        # Convert to DataFrame for easier filtering
        if not hasattr(weight_matrix, 'index'):
            df_matrix = pd.DataFrame(weight_matrix, index=cell_types, columns=cell_types)
        else:
            df_matrix = weight_matrix.copy()
        
        # Determine which cell types to keep based on filtering
        if sources_use is not None and targets_use is not None:
            # Keep only the specified sources and targets
            all_cells = list(set(sources_use) | set(targets_use))
            cells_to_keep = [c for c in cell_types if c in all_cells]
        elif sources_use is not None:
            # Keep only specified sources (for all targets)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        elif targets_use is not None:
            # Keep only specified targets (for all sources)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        else:
            cells_to_keep = cell_types
        
        # Filter the DataFrame and ensure consistency
        if sources_use is not None:
            # Set non-source rows to 0
            for idx in df_matrix.index:
                if idx not in sources_use:
                    df_matrix.loc[idx, :] = 0
        
        if targets_use is not None:
            # Set non-target columns to 0
            for col in df_matrix.columns:
                if col not in targets_use:
                    df_matrix.loc[:, col] = 0
        
        # Keep only cells that have connections after filtering
        if sources_use is not None and targets_use is not None:
            # When both are specified, keep only the union of sources and targets
            cells_with_connections = cells_to_keep
        else:
            # Find cells that still have connections
            row_sums = df_matrix.sum(axis=1)
            col_sums = df_matrix.sum(axis=0)
            cells_with_connections = list(set(
                list(df_matrix.index[row_sums > 0]) + 
                list(df_matrix.columns[col_sums > 0])
            ))
            # Preserve original order
            cells_with_connections = [c for c in cell_types if c in cells_with_connections]
        
        if cells_with_connections:
            df_matrix = df_matrix.loc[cells_with_connections, cells_with_connections]
            weight_matrix = df_matrix.values
            cell_types = cells_with_connections
    
    # If weight_matrix is a DataFrame, convert to numpy array
    if hasattr(weight_matrix, 'values'):
        weight_matrix = weight_matrix.values
    
    # Create directed graph
    G = nx.DiGraph()
    
    # Add nodes
    for i, cell in enumerate(cell_types):
        G.add_node(cell)
    
    # Add edges with weights
    for i, source in enumerate(cell_types):
        for j, target in enumerate(cell_types):
            if i != j and weight_matrix[i, j] > 0:
                G.add_edge(source, target, weight=float(weight_matrix[i, j]))
    
    # Remove isolated nodes if any
    G.remove_nodes_from(list(nx.isolates(G)))
    
    # If no nodes left (all isolated), return empty figure with message
    if len(G) == 0:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No communication links found above threshold", 
                ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.title(title, fontsize=16)
        return fig
    
    # Update cell_types to only include connected nodes
    cell_types = list(G.nodes())
    n_cells = len(cell_types)
    
    # Generate colors if not provided
    if group_colors is None:
        cmap = plt.cm.tab20
        group_colors = {cell: cmap(i % 20) for i, cell in enumerate(cell_types)}
    elif isinstance(group_colors, list):
        group_colors = {cell: group_colors[i % len(group_colors)] for i, cell in enumerate(cell_types)}
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Calculate node positions in a circle
    pos = {}
    radius = 5
    for i, node in enumerate(cell_types):
        theta = 2 * math.pi * i / n_cells
        x = radius * math.cos(theta)
        y = radius * math.sin(theta)
        pos[node] = (x, y)
    
    # Calculate edge widths based on weights
    edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
    max_weight = max(edge_weights) if edge_weights else 1
    
    # Scale widths for visibility
    scaled_widths = [1 + (w / max_weight) * 5 for w in edge_weights]
    
    # Draw nodes
    for node in G.nodes():
        x, y = pos[node]
        node_color = group_colors.get(node, (0.7, 0.7, 0.7, 1.0))
        
        # Draw node as circle
        circle = plt.Circle((x, y), 0.5, facecolor=node_color, edgecolor='black', zorder=10)
        ax.add_patch(circle)
        
        # Add node label
        offset = 0.2  # Label offset from node
        label_x = x * (1 + offset/radius)
        label_y = y * (1 + offset/radius)
        ax.text(label_x, label_y, node, ha='center', va='center', fontsize=12, 
                fontweight='bold', zorder=15)
    
    # Draw edges as curved arrows
    for i, (source, target) in enumerate(G.edges()):
        sx, sy = pos[source]
        tx, ty = pos[target]
        
        # Calculate control points for curved path
        # Find midpoint
        mx, my = (sx + tx) / 2, (sy + ty) / 2
        
        # Calculate perpendicular offset based on edge weight
        dx, dy = tx - sx, ty - sy
        length = math.sqrt(dx**2 + dy**2)
        
        # Make curves stronger for more visual separation
        curve_strength = 0.3 + (scaled_widths[i] - 1) * 0.1
        
        # Calculate control point by moving perpendicular to the line
        nx, ny = -dy/length, dx/length  # Normalized perpendicular vector
        cx = mx + nx * curve_strength * radius
        cy = my + ny * curve_strength * radius
        
        # Create curved path
        path = Path([(sx, sy), (cx, cy), (tx, ty)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
        
        # Arrow width and style based on weight
        width = scaled_widths[i] * 0.1
        
        # Use the source node color for the edge
        edge_color = group_colors.get(source, (0.7, 0.7, 0.7, 0.7))
        alpha = 0.7  # Transparency
        
        # Create patch with arrow
        arrow_patch = FancyArrowPatch(
            path=path,
            arrowstyle=f'-|>,head_length=10,head_width={10 + 5*width}',
            connectionstyle=f'arc3,rad=0.3',
            linewidth=width*3,
            color=edge_color,
            alpha=alpha,
            zorder=5
        )
        ax.add_patch(arrow_patch)
    
    # Add title
    plt.title(title, fontsize=16)
    
    # Create legend
    legend_elements = []
    for node, color in group_colors.items():
        if node in G:  # Only include connected nodes
            legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                             markerfacecolor=color, markersize=10, label=node))
    
    # Add legend outside the main plot area
    ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
    
    # Set limits and remove axes
    margin = 1.5
    ax.set_xlim(-radius*margin, radius*margin)
    ax.set_ylim(-radius*margin, radius*margin)
    ax.set_aspect('equal')
    ax.axis('off')
    
    plt.tight_layout()
    return fig

def netVisual_chord_by_pathway(net, pathway_name, group_colors=None, 
                                         title_prefix="Signaling Pathway: ", figsize=(10, 10)):
    """
    Create a chord diagram for a specific pathway using matplotlib
    
    Parameters
    ----------
    net : dict
        Dictionary containing interaction data with 'netP' and pathway information
    pathway_name : str
        Name of the pathway to visualize
    group_colors : dict or list, optional
        Colors for each cell group
    title_prefix : str, optional
        Prefix for the plot titles
    figsize : tuple, optional
        Figure size (width, height) in inches
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import networkx as nx
    from matplotlib.patches import FancyArrowPatch
    from matplotlib.path import Path
    import matplotlib.patches as patches
    import math
    
    try:
        # Find the pathway index in the netP list
        if 'netP' not in net or 'pathways' not in net['netP']:
            raise ValueError("No pathway information found in the provided network data")
        
        pathways = net['netP']['pathways']
        if isinstance(pathways, np.ndarray):
            pathways = pathways.tolist()
        
        if pathway_name not in pathways:
            raise ValueError(f"Pathway '{pathway_name}' not found in the dataset")
        
        pathway_idx = pathways.index(pathway_name)
        
        # Extract interaction matrix for this pathway
        if 'prob' not in net['netP']:
            raise ValueError("Probability information missing for pathways")
        
        pathway_matrix = net['netP']['prob'][:, :, pathway_idx]
        
        # Get cell types
        if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 1:
            cell_types = net['net']['dimnames'][0]
        else:
            cell_types = [f"Cell{i+1}" for i in range(pathway_matrix.shape[0])]
        
        # Convert cell_types to list if it's not already
        if not isinstance(cell_types, list):
            cell_types = list(cell_types)
        
        # Create directed graph
        G = nx.DiGraph()
        
        # Add nodes
        for i, cell in enumerate(cell_types):
            G.add_node(cell)
        
        # Add edges with weights
        for i, source in enumerate(cell_types):
            for j, target in enumerate(cell_types):
                if i != j and pathway_matrix[i, j] > 0:
                    G.add_edge(source, target, weight=float(pathway_matrix[i, j]))
        
        # Remove isolated nodes if any
        G.remove_nodes_from(list(nx.isolates(G)))
        
        # If no nodes left (all isolated), return empty figure with message
        if len(G) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, f"No communication links found for '{pathway_name}' pathway", 
                   ha='center', va='center', fontsize=14)
            ax.axis('off')
            plt.title(f"{title_prefix}{pathway_name}", fontsize=16)
            return fig
        
        # Update cell_types to only include connected nodes
        cell_types = list(G.nodes())
        n_cells = len(cell_types)
        
        # Generate colors if not provided
        if group_colors is None:
            cmap = plt.cm.tab20
            group_colors = {cell: cmap(i % 20) for i, cell in enumerate(cell_types)}
        elif isinstance(group_colors, list):
            group_colors = {cell: group_colors[i % len(group_colors)] for i, cell in enumerate(cell_types)}
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Calculate node positions in a circle
        pos = {}
        radius = 5
        for i, node in enumerate(cell_types):
            theta = 2 * math.pi * i / n_cells
            x = radius * math.cos(theta)
            y = radius * math.sin(theta)
            pos[node] = (x, y)
        
        # Calculate edge widths based on weights
        edge_weights = [G[u][v]['weight'] for u, v in G.edges()]
        max_weight = max(edge_weights) if edge_weights else 1
        
        # Scale widths for visibility
        scaled_widths = [1 + (w / max_weight) * 5 for w in edge_weights]
        
        # Draw nodes
        for node in G.nodes():
            x, y = pos[node]
            node_color = group_colors.get(node, (0.7, 0.7, 0.7, 1.0))
            
            # Draw node as circle
            circle = plt.Circle((x, y), 0.5, facecolor=node_color, edgecolor='black', zorder=10)
            ax.add_patch(circle)
            
            # Add node label
            offset = 0.2  # Label offset from node
            label_x = x * (1 + offset/radius)
            label_y = y * (1 + offset/radius)
            ax.text(label_x, label_y, node, ha='center', va='center', fontsize=12, 
                   fontweight='bold', zorder=15)
        
        # Draw edges as curved arrows
        for i, (source, target) in enumerate(G.edges()):
            sx, sy = pos[source]
            tx, ty = pos[target]
            
            # Calculate control points for curved path
            # Find midpoint
            mx, my = (sx + tx) / 2, (sy + ty) / 2
            
            # Calculate perpendicular offset based on edge weight
            dx, dy = tx - sx, ty - sy
            length = math.sqrt(dx**2 + dy**2)
            
            # Make curves stronger for more visual separation
            curve_strength = 0.3 + (scaled_widths[i] - 1) * 0.1
            
            # Calculate control point by moving perpendicular to the line
            nx, ny = -dy/length, dx/length  # Normalized perpendicular vector
            cx = mx + nx * curve_strength * radius
            cy = my + ny * curve_strength * radius
            
            # Create curved path
            path = Path([(sx, sy), (cx, cy), (tx, ty)], [Path.MOVETO, Path.CURVE3, Path.CURVE3])
            
            # Arrow width and style based on weight
            width = scaled_widths[i] * 0.1
            
            # Use the source node color for the edge
            edge_color = group_colors.get(source, (0.7, 0.7, 0.7, 0.7))
            alpha = 0.7  # Transparency
            
            # Create patch with arrow
            arrow_patch = FancyArrowPatch(
                path=path,
                arrowstyle=f'-|>,head_length=10,head_width={10 + 5*width}',
                connectionstyle=f'arc3,rad=0.3',
                linewidth=width*3,
                color=edge_color,
                alpha=alpha,
                zorder=5
            )
            ax.add_patch(arrow_patch)
        
        # Add title
        plt.title(f"{title_prefix}{pathway_name}", fontsize=16)
        
        # Create legend
        legend_elements = []
        for node, color in group_colors.items():
            if node in G:  # Only include connected nodes
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', 
                                                markerfacecolor=color, markersize=10, label=node))
        
        # Add legend outside the main plot area
        ax.legend(handles=legend_elements, loc='upper left', bbox_to_anchor=(1, 1))
        
        # Set limits and remove axes
        margin = 1.5
        ax.set_xlim(-radius*margin, radius*margin)
        ax.set_ylim(-radius*margin, radius*margin)
        ax.set_aspect('equal')
        ax.axis('off')
        
        plt.tight_layout()
        return fig
        
    except Exception as e:
        import traceback
        print(f"Error generating chord diagram: {str(e)}")
        print(traceback.format_exc())
        
        # Return an empty figure with error message
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"Error: {str(e)}", 
               ha='center', va='center', fontsize=14, color='red')
        ax.axis('off')
        plt.title(f"{title_prefix}{pathway_name}", fontsize=16)
        return fig

def netVisual_heatmap(
    net, 
    measure="weight", 
    color_heatmap="Oranges", 
    font_size=14, 
    font_size_title=20, 
    sorted_order=None,
    title=None, 
    annot=True,
    color_use=None,  # 細胞名と色のマッピング辞書
    show_color_bar=False,  # カラーバーを表示するか
    vmin=None,  # 新しいパラメータ: カラースケールの最小値
    vmax=None,
    figx=12,
    figy=10,  # 新しいパラメータ: カラースケールの最大値
    sources_use=None,  # 新しいパラメータ: 使用するソース細胞タイプのリスト
    targets_use=None  # 新しいパラメータ: 使用するターゲット細胞タイプのリスト
):
    """
    Create a heatmap for cell-cell communication
    Similar to netVisual_heatmap in R CellChat package
    
    Parameters
    ----------
    net : dict
        Network object containing communication data
    measure : str, optional
        Measure to use for heatmap, options are 'weight' or 'count'
    color_heatmap : str or list, optional
        Color scheme for the heatmap. Can be a matplotlib colormap name or a list of colors for diverging data
    font_size : int, optional
        Font size for axis labels and text
    font_size_title : int, optional
        Font size for the plot title
    sorted_order : list or None
        細胞タイプの並び順。指定があれば順序を変更
    title : str or None
        図のタイトル。Noneの場合は自動生成
    annot : bool
        ヒートマップにデータ値を表示するか
    color_use : dict or None
        細胞名と色のマッピング辞書。カラーバー表示に使用
    show_color_bar : bool
        カラーバーを軸に表示するか
    sources_use : list or None
        使用するソース細胞タイプのリスト。Noneの場合は全て使用
    targets_use : list or None
        使用するターゲット細胞タイプのリスト。Noneの場合は全て使用
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """    
    def mix_colors(color1, color2, ratio):
        """2つの色を指定された比率で混合する"""
        # 色をRGBAに変換
        c1 = mcolors.to_rgba(color1)
        c2 = mcolors.to_rgba(color2)
        
        # 混合
        mixed = tuple(c1[i] * (1 - ratio) + c2[i] * ratio for i in range(4))
        return mixed
    
    # Extract the appropriate matrix based on different possible structures
    matrix = None
    cell_types = None
    
    # Check if net contains 'network' with 'strength_matrix'/'count_matrix' structure
    if matrix is None and 'network' in net:
        if measure == "weight" and 'strength_matrix' in net['network']:
            matrix = net['network']['strength_matrix']
            if title is None:
                title = "Interaction weights/strength"
        elif measure == "count" and 'count_matrix' in net['network']:
            matrix = net['network']['count_matrix']
            if title is None:
                title = "Number of interactions"
    
    # If no cell types were found earlier, try to extract from matrix
    if cell_types is None:
        if hasattr(matrix, 'index'):  # It's a DataFrame
            cell_types = matrix.index
        else:
            # No cell type info available, create generic names
            cell_types = [f"Cell{i+1}" for i in range(matrix.shape[0])]
    
    # Convert to DataFrame if it's not already
    if not isinstance(matrix, pd.DataFrame):
        matrix = pd.DataFrame(matrix, index=cell_types, columns=cell_types)
    
    # STEP 1: Filter the dataframe by sources (rows) and targets (columns)
    # Both sources_use and targets_use are always specified together when filtering is enabled
    if sources_use is not None and targets_use is not None:
        # Filter sources (rows) and targets (columns)
        valid_sources = [s for s in sources_use if s in matrix.index]
        valid_targets = [t for t in targets_use if t in matrix.columns]
        
        # Filter matrix by sources (rows) and targets (columns)
        if valid_sources and valid_targets:
            matrix = matrix.loc[valid_sources, valid_targets].copy()
        else:
            # No valid cells found, return empty heatmap
            fig = plt.figure(figsize=(figx, figy))
            ax = fig.add_subplot(111)
            ax.text(0.5, 0.5, 'No valid sources or targets found.\nPlease check your cell type selections.', 
                    horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
                    fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.axis('off')
            return fig
    
    # STEP 2: Apply sorted_order to both rows and columns separately
    if sorted_order is not None:
        current_sources = list(matrix.index)
        current_targets = list(matrix.columns)
        
        # Reorder sources (rows)
        ordered_sources = [c for c in sorted_order if c in current_sources]
        remaining_sources = [c for c in current_sources if c not in ordered_sources]
        final_source_order = ordered_sources + remaining_sources
        
        # Reorder targets (columns)
        ordered_targets = [c for c in sorted_order if c in current_targets]
        remaining_targets = [c for c in current_targets if c not in ordered_targets]
        final_target_order = ordered_targets + remaining_targets
        
        # Apply reordering
        if final_source_order and final_target_order:
            matrix = matrix.reindex(index=final_source_order, columns=final_target_order)
    
    # Check if matrix is empty or invalid after filtering and ordering
    if matrix is None or matrix.empty or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        fig = plt.figure(figsize=(figx, figy))
        ax = fig.add_subplot(111)
        ax.text(0.5, 0.5, 'No data available after filtering/ordering.\nPlease check your cell type selections.', 
                horizontalalignment='center', verticalalignment='center', transform=ax.transAxes,
                fontsize=12, bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgray"))
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        return fig
    
    # レイアウト設定 - 値のcolorbarは常に表示
    fig = plt.figure(figsize=(figx, figy))
    
    if show_color_bar and color_use is not None:
        # 細胞タイプのカラーバー付きレイアウト
        main_width_ratio = max(6, figx * 1.2)  # メインプロットの幅比率（figxに比例）
        main_height_ratio = max(6, figy * 1.2)  # メインプロットの高さ比率（figyに比例）
        
        gs = plt.GridSpec(6, 4, 
                         width_ratios=[0.6, main_width_ratio, 0.3, 0.6], 
                         height_ratios=[0.6, main_height_ratio * 0.25, main_height_ratio * 0.5, main_height_ratio * 0.25, 0.3, 0.1],
                         wspace=0.01, hspace=0.01)
        
        # 各領域のサブプロット
        ax_corner = fig.add_subplot(gs[0, 0])
        ax_corner.set_visible(False)
        
        # 上部のカラーバー（X軸の細胞タイプ）
        ax_top = fig.add_subplot(gs[0, 1])
        
        # 左側のカラーバー（Y軸の細胞タイプ）
        ax_left = fig.add_subplot(gs[1:4, 0])  # ヒートマップと同じ高さ
        
        # メインのヒートマップ
        ax_heatmap = fig.add_subplot(gs[1:4, 1])  # 3行分を使用
        
        # 右側のカラーバー（値の凡例用）- heatmapの縦サイズの半分程度
        cbar_ax = fig.add_subplot(gs[2, 3])  # 中央の1行のみ使用
        
    else:
        # シンプルレイアウト（値のcolorbarのみ）
        main_width_ratio = max(6, figx * 1.2)
        main_height_ratio = max(6, figy * 1.2)
        
        gs = plt.GridSpec(4, 3, 
                         width_ratios=[main_width_ratio, 0.2, 0.4], 
                         height_ratios=[0.2, main_height_ratio * 0.25, main_height_ratio * 0.5, main_height_ratio * 0.25],
                         wspace=0.05, hspace=0.02)
        
        # メインのヒートマップ
        ax_heatmap = fig.add_subplot(gs[1:4, 0])  # 3行分を使用
        
        # 右側のカラーバー（値の凡例用）- heatmapの縦サイズの半分程度
        cbar_ax = fig.add_subplot(gs[2, 2])  # 中央の1行のみ使用
        
        # cell-type colorbarは無し
        ax_top = None
        ax_left = None

    fmt = ".2f" if measure == "weight" else ".0f"

    # カラーマップの設定
    # 正負両方の値がある場合は発散カラーマップを使用
    if np.min(matrix.values) < 0:
        # 発散カラーマップ用の設定
        if isinstance(color_heatmap, str):
            if color_heatmap == "Oranges" or color_heatmap == "Reds":
                # デフォルトカラーの場合、青赤の発散カラーマップに変更
                color_heatmap = ['#2166ac', '#f7f7f7', '#b2182b']
        
        # カラーマップの作成
        if isinstance(color_heatmap, list) and len(color_heatmap) >= 3:
            # カスタムカラーマップの作成
            neg_color = color_heatmap[0]
            mid_color = color_heatmap[1]
            pos_color = color_heatmap[2]
            n_bins = 100
            colors = []
            # 負の部分のカラーグラデーション
            for i in range(n_bins // 2):
                r = i / (n_bins // 2)
                colors.append(mix_colors(neg_color, mid_color, r))
            # 正の部分のカラーグラデーション
            for i in range(n_bins // 2):
                r = i / (n_bins // 2)
                colors.append(mix_colors(mid_color, pos_color, r))
            cmap = LinearSegmentedColormap.from_list('custom_diverging', colors, N=n_bins)
        else:
            # デフォルトの発散カラーマップ
            cmap = 'RdBu_r'
        
        # ゼロを中心にした境界値
        vmin = np.min(matrix.values)
        vmax = np.max(matrix.values)
        center = 0

        # seabornのheatmapにcbar_axパラメータを直接渡す
        # 常にcbar_axを使用（cbar_axは常に定義されている）
        sns.heatmap(
            matrix,
            cmap=cmap,
            center=center,
            vmin=vmin,
            vmax=vmax,
            square=False,
            linewidths=0.5,
            cbar_ax=cbar_ax,  # 常にcbar_axを使用
            annot=annot,
            fmt=fmt,
            ax=ax_heatmap
        )
    else:
        # 正の値のみの場合は通常のカラーマップを使用
        # 常にcbar_axを使用（cbar_axは常に定義されている）
        sns.heatmap(
            matrix,
            cmap=color_heatmap,
            square=False,
            linewidths=0.5,
            cbar_ax=cbar_ax,  # 常にcbar_axを使用
            annot=annot,
            fmt=fmt,
            ax=ax_heatmap,
            vmin=vmin,
            vmax=vmax,
        )
    
    # 細胞タイプのカラーバーを描画（オプション）
    if show_color_bar and color_use is not None:
        # ヒートマップの座標位置を取得
        heatmap_pos = ax_heatmap.get_position()
        heatmap_width = heatmap_pos.width
        heatmap_height = heatmap_pos.height
        
        # 1. 上部のカラーバー作成（X軸の細胞タイプ）
        if ax_top is not None:
            # カラーバーの位置をヒートマップに合わせる
            top_pos = ax_top.get_position()
            ax_top.set_position([heatmap_pos.x0, top_pos.y0, heatmap_width, top_pos.height])
            
            # 色配列の作成
            cell_colors_x = []
            for cell in matrix.columns:
                if cell in color_use:
                    cell_colors_x.append(color_use[cell])
                else:
                    # デフォルト色（グレー）
                    cell_colors_x.append((0.7, 0.7, 0.7, 1.0))
            
            # 細胞タイプのカラーバーを描画
            color_array_x = np.array([cell_colors_x])
            ax_top.imshow(color_array_x, aspect='auto', interpolation='nearest')
            ax_top.set_xticks([])
            ax_top.set_yticks([])
            ax_top.axis('off')
        
        # 2. 左側のカラーバー作成（Y軸の細胞タイプ）
        if ax_left is not None:
            # カラーバーの位置をヒートマップに合わせる
            left_pos = ax_left.get_position()
            ax_left.set_position([left_pos.x0, heatmap_pos.y0, left_pos.width, heatmap_height])
            
            # 色配列の作成（Y軸は転置して使用）
            cell_colors_y = []
            for cell in matrix.index:
                if cell in color_use:
                    cell_colors_y.append(color_use[cell])
                else:
                    # デフォルト色（グレー）
                    cell_colors_y.append((0.7, 0.7, 0.7, 1.0))
            
            # Y軸用には色配列を転置して表示（縦方向）
            color_array_y = np.array([cell_colors_y]).T  # 転置して縦方向に
            ax_left.imshow(color_array_y, aspect='auto', interpolation='nearest')
            ax_left.set_xticks([])
            ax_left.set_yticks([])
            ax_left.axis('off')
    
    # 軸ラベルのフォーマット - すべて黒色に
    ax_heatmap.set_xlabel('')
    ax_heatmap.set_ylabel('')
    ax_heatmap.tick_params(axis='both', labelsize=font_size, colors='black')
    
    # すべてのラベルを黒色に設定
    for tick_label in ax_heatmap.get_xticklabels():
        tick_label.set_color('black')
    
    for tick_label in ax_heatmap.get_yticklabels():
        tick_label.set_color('black')
    
    # タイトル設定
    ax_heatmap.set_title(title, fontsize=font_size_title)
    
    # ヒートマップのアスペクト比を自動調整（matplotlibに任せる）
    # square=Falseに設定したので、軸のサイズに合わせて自動調整される
    ax_heatmap.set_aspect('auto')
    
    # X軸ラベルの回転
    plt.xticks(rotation=45, ha='right')
    plt.yticks(rotation=0)

    ax_heatmap.set_xlabel('Targets (Receiver)', fontsize=font_size+2)
    ax_heatmap.set_ylabel('Sources (Sender)', fontsize=font_size+2)
    
    # Apply tight_layout before final size adjustment
    plt.tight_layout()
    
    # Ensure the figure size is correctly set
    try:
        current_size = fig.get_size_inches()
        if abs(current_size[0] - figx) > 0.1 or abs(current_size[1] - figy) > 0.1:
            fig.set_size_inches(figx, figy, forward=True)
    except Exception as e:
        pass  # Silently continue if size setting fails
    
    return fig


def mix_colors(color1, color2, ratio):
    """
    2つの色を指定された比率で混合する
    
    Parameters
    ----------
    color1 : str
        最初の色（HTMLカラーコードまたは名前）
    color2 : str
        2番目の色（HTMLカラーコードまたは名前）
    ratio : float
        混合比率（0: 100% color1, 1: 100% color2）
    
    Returns
    -------
    tuple
        混合色のRGBAタプル
    """
    import matplotlib.colors as mcolors
    
    # 色をRGBAに変換
    c1 = mcolors.to_rgba(color1)
    c2 = mcolors.to_rgba(color2)
    
    # 混合
    mixed = tuple(c1[i] * (1 - ratio) + c2[i] * ratio for i in range(4))
    return mixed



def plotGeneExpression(cellchat, signaling=None, adata=None, features=None, group_by=None, 
                      color_use=None, threshold=0.05, title=None, width=10, height=8):
    """
    Plot the expression of genes involved in a particular signaling pathway
    Similar to plotGeneExpression in R CellChat package
    
    Parameters
    ----------
    cellchat : dict
        CellChat results object
    signaling : str, optional
        Name of the signaling pathway to visualize
    adata : AnnData, optional
        AnnData object with gene expression data
    features : list, optional
        List of genes to visualize (if not using pathway)
    group_by : str, optional
        Column in adata.obs for grouping cells
    color_use : dict, optional
        Colors for each cell group
    threshold : float, optional
        Significance threshold for pathways
    title : str, optional
        Title for the plot
    width : int, optional
        Width of the plot in inches
    height : int, optional
        Height of the plot in inches
    
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object
    """
    try:
        import scanpy as sc
    except ImportError:
        raise ImportError("scanpy is required for plotGeneExpression. Please install it with 'pip install scanpy'.")
    
    # Get the AnnData object from cellchat if not provided
    if adata is None:
        if 'adata' in cellchat:
            adata = cellchat['adata']
        else:
            raise ValueError("AnnData object not found in cellchat and not provided")
    
    # Get group_by information
    if group_by is None:
        # Try to infer from cellchat object
        if hasattr(cellchat, 'groupby'):
            group_by = cellchat.groupby
        elif isinstance(cellchat, dict) and 'groupby' in cellchat:
            group_by = cellchat['groupby']
        else:
            # Try to find a categorical column in adata.obs
            cat_cols = [col for col in adata.obs.columns if pd.api.types.is_categorical_dtype(adata.obs[col])]
            if cat_cols:
                group_by = cat_cols[0]
            else:
                raise ValueError("No grouping variable found; please provide group_by")
    
    # Check if the group_by column exists
    if group_by not in adata.obs.columns:
        raise ValueError(f"Column '{group_by}' not found in adata.obs")
    
    # Get genes to visualize
    genes_to_plot = features
    
    # If signaling pathway is provided, find relevant genes
    if signaling is not None and genes_to_plot is None:
        # Find pathway in cellchat object
        if 'netP' not in cellchat or 'pathways' not in cellchat['netP']:
            raise ValueError("Pathway information not found in cellchat object")
        
        pathways = cellchat['netP']['pathways']
        if signaling not in pathways:
            raise ValueError(f"Pathway '{signaling}' not found in the dataset")
        
        # Get pathway index
        pathway_idx = np.where(pathways == signaling)[0][0]
        
        # Get significant interactions for this pathway
        if 'pval' in cellchat['netP'] and 'prob' in cellchat['netP']:
            pathway_pval = cellchat['netP']['pval'][:, :, pathway_idx]
            pathway_prob = cellchat['netP']['prob'][:, :, pathway_idx]
            
            # Filter by threshold
            sig_idx = np.where((pathway_pval <= threshold) & (pathway_prob > 0))
            
            if len(sig_idx[0]) == 0:
                raise ValueError(f"No significant interactions found for pathway {signaling}")
            
            # Get corresponding LR pairs
            if 'results' in cellchat and not cellchat['results'].empty:
                # Filter results for this pathway and significant interactions
                if 'interaction_name' in cellchat['results'].columns and 'ligand' in cellchat['results'].columns and 'receptor' in cellchat['results'].columns:
                    # Get pathway interactions
                    pathway_interactions = cellchat['results'][
                        cellchat['results']['interaction_name'].str.contains(signaling, case=False, na=False)]
                    
                    # Extract ligands and receptors
                    ligands = pathway_interactions['ligand'].unique()
                    receptors = pathway_interactions['receptor'].unique()
                    
                    # Combine and clean
                    genes_to_plot = list(set(list(ligands) + list(receptors)))
                    genes_to_plot = [g for g in genes_to_plot if g in adata.var_names]
                else:
                    raise ValueError("Required columns not found in results")
            else:
                # Fallback if results are not available
                raise ValueError("Detailed results not found in cellchat object")
    
    if genes_to_plot is None or len(genes_to_plot) == 0:
        raise ValueError("No genes found to plot")
    
    # Filter genes by availability in the AnnData object
    genes_to_plot = [g for g in genes_to_plot if g in adata.var_names]
    
    if len(genes_to_plot) == 0:
        raise ValueError("None of the selected genes found in the AnnData object")
    
    # Calculate mean expression by group
    adata_subset = adata[:, genes_to_plot].copy()
    
    if len(adata_subset) == 0:
        raise ValueError("No data found for the selected genes")
    
    # Use scanpy's function to calculate mean expression by group
    sc.pp.normalize_total(adata_subset, target_sum=1e4)
    sc.pp.log1p(adata_subset)
    
    # Calculate mean expression per group
    group_means = pd.DataFrame(
        {gene: np.zeros(len(adata_subset.obs[group_by].unique())) for gene in genes_to_plot},
        index=sorted(adata_subset.obs[group_by].unique())
    )
    
    for group in group_means.index:
        cells = adata_subset[adata_subset.obs[group_by] == group].X
        # If it's a sparse matrix, convert to dense
        if scipy.sparse.issparse(cells):
            cells = cells.toarray()
        group_means.loc[group] = cells.mean(axis=0)
    
    # Scale expression to 0-1 range for visualization
    scaled_means = group_means.copy()
    for col in scaled_means.columns:
        min_val = scaled_means[col].min()
        max_val = scaled_means[col].max()
        if max_val > min_val:
            scaled_means[col] = (scaled_means[col] - min_val) / (max_val - min_val)
    
    # Determine layout
    n_genes = len(genes_to_plot)
    n_cols = min(5, n_genes)
    n_rows = int(np.ceil(n_genes / n_cols))
    
    # Create figure and subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(width, height))
    
    # Handle case of single gene
    if n_genes == 1:
        axes = np.array([axes])
    
    # Flatten axes for easier iteration
    if n_rows > 1 and n_cols > 1:
        axes = axes.flatten()
    
    # Plot each gene
    for i, gene in enumerate(genes_to_plot):
        if i < len(axes):
            # Get expression data for this gene
            gene_data = pd.DataFrame({'Expression': group_means[gene]})
            gene_data['Group'] = gene_data.index
            
            # Create violin plot
            if len(gene_data) > 1:
                sns.barplot(x='Group', y='Expression', data=gene_data, 
                           palette=color_use, ax=axes[i])
                axes[i].set_title(gene, fontsize=10)
                axes[i].set_ylabel('Mean Expression', fontsize=8)
                axes[i].set_xlabel('')
                
                # Rotate x-axis labels if there are many groups
                if len(gene_data) > 5:
                    axes[i].set_xticklabels(axes[i].get_xticklabels(), rotation=45, ha='right')
            else:
                axes[i].text(0.5, 0.5, f"Insufficient data for {gene}", 
                             ha='center', va='center', fontsize=10)
                axes[i].axis('off')
    
    # Remove unused subplots
    for i in range(n_genes, len(axes)):
        fig.delaxes(axes[i])
    
    # Add overall title if provided or signaling is specified
    if title:
        plt.suptitle(title, fontsize=14)
    elif signaling:
        plt.suptitle(f"Gene Expression in {signaling} Signaling Pathway", fontsize=14)
    
    plt.tight_layout()
    
    return fig


def create_chord_diagram_pycirclize(
    net, 
    pathway_name=None, 
    figsize=(10, 10), 
    cmap_name="tab10", 
    sorted_order=None, 
    measure="weight",
    color_use=None,
    space=5,  # 細胞間の空白 (度)
    alpha_edge=0.7,  # リンクの透明度
    show_edge_border=True,  # リンクの黒い輪郭を表示するかどうか
    edge_border_width=0.5,
    sources_use=None,  # 使用するソース細胞タイプのリスト
    targets_use=None  # 使用するターゲット細胞タイプのリスト
):
    """
    Create chord diagram using pycirclize
    
    Parameters
    ----------
    net : dict
        ネットワークオブジェクト
    pathway_name : str or list, optional
        パスウェイ名（またはパスウェイ名のリスト）
    figsize : tuple, optional
        図のサイズ
    cmap_name : str, optional
        カラーマップ名
    sorted_order : list, optional
        細胞の順序
    measure : str, optional
        測定値 ("weight" or "count")
    color_use : dict, optional
        細胞名と色のマッピング辞書
    space : float, optional
        細胞間の空白 (度)
    alpha_edge : float, optional
        リンクの透明度 (0.0〜1.0)
    show_edge_border : bool, optional
        リンクの黒い輪郭を表示するかどうか
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        図オブジェクト
    """
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.colors as mcolors
    import traceback

    # パスウェイ名を整理
    # もし文字列なら配列に変換、Noneなら["Aggregate"]に設定
    if pathway_name is None:
        pathway_name = ["Aggregate"]
    elif isinstance(pathway_name, str):
        pathway_name = [pathway_name]
    
    # 行列とセルタイプを抽出
    if pathway_name[0] == "Aggregate":
        # 集計モード
        if 'network' not in net:
            raise ValueError("Aggregated network information is not available.")
        if measure == "weight":
            if 'strength_matrix' not in net['network']:
                raise ValueError("Strength matrix is not available.")
            matrix = net['network']['strength_matrix']
        else:
            if 'count_matrix' not in net['network']:
                raise ValueError("Count matrix is not available.")
            matrix = net['network']['count_matrix']
        if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            cell_types = net['net']['dimnames'][0]
        else:
            cell_types = [f"Cell{i+1}" for i in range(matrix.shape[0])]
        title = "Overall Cell-Cell Interactions (Aggregate)"
    else:
        # パスウェイモード
        pathways = net['netP']['pathways']
        if isinstance(pathways, np.ndarray):
            pathways = pathways.tolist()
        
        # パスウェイインデックスを取得
        pathway_indices = []
        valid_pathways = []
        
        for p in pathway_name:
            try:
                idx = pathways.index(p)
                pathway_indices.append(idx)
                valid_pathways.append(p)
            except ValueError:
                print(f"Warning: Pathway '{p}' not found in dataset")
        
        if not pathway_indices:
            raise ValueError(f"None of the specified pathways were found.")
        
        # 確率行列の集計
        if 'prob' in net['netP']:
            prob_arrays = []
            for idx in pathway_indices:
                prob_arrays.append(net['netP']['prob'][:, :, idx])
            
            matrix = np.sum(prob_arrays, axis=0)
        else:
            raise ValueError("Probability matrix not found in netP")
        
        if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            cell_types = net['net']['dimnames'][0]
        else:
            cell_types = [f"Cell{i+1}" for i in range(matrix.shape[0])]
        
        # タイトル
        if len(valid_pathways) == 1:
            title = f"Cell-Cell Interaction: {valid_pathways[0]}"
        else:
            if len(valid_pathways) <= 3:
                title = f"Combined Pathways: {', '.join(valid_pathways)}"
            else:
                title = f"Combined {len(valid_pathways)} Pathways"
    
    # リスト化
    if not isinstance(cell_types, list):
        cell_types = list(cell_types)
    
    # 配列に変換
    if not isinstance(matrix, np.ndarray):
        matrix = np.array(matrix)
    
    # Filter by sources_use and targets_use if specified
    if sources_use is not None or targets_use is not None:
        # Create a DataFrame for easier filtering
        df_matrix = pd.DataFrame(matrix, index=cell_types, columns=cell_types)
        
        # Determine which cell types to keep based on filtering
        if sources_use is not None and targets_use is not None:
            # Keep only the specified sources and targets
            all_cells = list(set(sources_use) | set(targets_use))
            cells_to_keep = [c for c in cell_types if c in all_cells]
        elif sources_use is not None:
            # Keep only specified sources (for all targets)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        elif targets_use is not None:
            # Keep only specified targets (for all sources)
            cells_to_keep = cell_types  # Keep all cells but filter matrix later
        else:
            cells_to_keep = cell_types
        
        # Filter the DataFrame and ensure consistency
        if sources_use is not None:
            # Set non-source rows to 0
            for idx in df_matrix.index:
                if idx not in sources_use:
                    df_matrix.loc[idx, :] = 0
        
        if targets_use is not None:
            # Set non-target columns to 0
            for col in df_matrix.columns:
                if col not in targets_use:
                    df_matrix.loc[:, col] = 0
        
        # Keep only cells that have connections after filtering
        if sources_use is not None and targets_use is not None:
            # When both are specified, keep only the union of sources and targets
            cells_with_connections = cells_to_keep
        else:
            # Find cells that still have connections
            row_sums = df_matrix.sum(axis=1)
            col_sums = df_matrix.sum(axis=0)
            cells_with_connections = list(set(
                list(df_matrix.index[row_sums > 0]) + 
                list(df_matrix.columns[col_sums > 0])
            ))
            # Preserve original order
            cells_with_connections = [c for c in cell_types if c in cells_with_connections]
        
        if cells_with_connections:
            df_matrix = df_matrix.loc[cells_with_connections, cells_with_connections]
            matrix = df_matrix.values
            cell_types = cells_with_connections
    
    # 孤立ノードのフィルタリング
    total_weights = [matrix[i, :].sum() + matrix[:, i].sum() for i in range(len(cell_types))]
    connected_indices = [i for i, w in enumerate(total_weights) if w > 0]
    
    if not connected_indices:
        # 接続がない場合、空の図を返す
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"'{title}' に接続が見つかりません", ha='center', va='center', fontsize=14)
        ax.axis('off')
        plt.title(title, fontsize=16)
        return fig
    
    # 接続のある細胞のみを保持
    connected_cell_types = [cell_types[i] for i in connected_indices]
    connected_matrix = matrix[np.ix_(connected_indices, connected_indices)]
    
    # 並び替え処理：色の順番は完全にcell orderの順番に合わせる
    if sources_use is not None and targets_use is not None:
        # source/targetに含まれる細胞の和集合を取得
        relevant_cells = set(sources_use) | set(targets_use)
        
        if sorted_order is not None:
            # sorted_orderの中で、source/targetに含まれ、かつ接続のある細胞のみを順番通りに取得
            final_cell_types = [cell for cell in sorted_order 
                              if cell in relevant_cells and cell in connected_cell_types]
            print(f"[CHORD DEBUG] sorted_order: {sorted_order}")
            print(f"[CHORD DEBUG] sources_use: {sources_use}")
            print(f"[CHORD DEBUG] targets_use: {targets_use}")
            print(f"[CHORD DEBUG] relevant_cells (union): {relevant_cells}")
            print(f"[CHORD DEBUG] final_cell_types (sorted by cell order): {final_cell_types}")
        else:
            # sorted_orderがない場合、connected_cell_typesの順序を使用
            final_cell_types = [cell for cell in connected_cell_types if cell in relevant_cells]
            print(f"[CHORD DEBUG] Without sorted_order - final_cell_types: {final_cell_types}")
    elif sorted_order is not None:
        # 通常のsorted_order処理
        valid_cells = [cell for cell in sorted_order if cell in connected_cell_types]
        remaining_cells = [cell for cell in connected_cell_types if cell not in valid_cells]
        final_cell_types = valid_cells + remaining_cells
        print(f"[CHORD DEBUG] Normal sorted_order - final_cell_types: {final_cell_types}")
    else:
        final_cell_types = connected_cell_types
        print(f"[CHORD DEBUG] No ordering - final_cell_types: {final_cell_types}")
    
    # インデックスマッピングを作成して行列を並び替え
    if final_cell_types != connected_cell_types:
        idx_map = {cell: i for i, cell in enumerate(connected_cell_types)}
        reorder_idx = [idx_map[cell] for cell in final_cell_types]
        final_matrix = connected_matrix[np.ix_(reorder_idx, reorder_idx)]
    else:
        final_matrix = connected_matrix
    
    # 色のマッピング設定（並び替え後の順序で色を生成）
    # 常に新しい色を生成（既存のcolor_useを無視してテスト）
    print(f"[CHORD DEBUG] Original color_use: {color_use}")
    
    # カラーマップから色を生成
    try:
        cmap = cm.get_cmap(cmap_name)
    except:
        cmap = cm.get_cmap('tab10')
    
    # 色のマッピングを作成（並び替え後の順序で色を割り当て）
    new_color_use = {}
    for i, cell in enumerate(final_cell_types):
        n_colors = min(cmap.N, 20) if hasattr(cmap, 'N') else 10
        new_color_use[cell] = cmap(i % n_colors / max(1, n_colors - 1))
    
    print(f"[CHORD DEBUG] Color assignment order: {final_cell_types}")
    print(f"[CHORD DEBUG] New colors: {[(cell, new_color_use[cell]) for cell in final_cell_types]}")
    
    # テスト用に新しい色を使用
    color_use = new_color_use
    
    # DataFrameを作成
    matrix_df = pd.DataFrame(final_matrix, index=final_cell_types, columns=final_cell_types)
    
    try:
        # フィギュアサイズを設定
        plt.figure(figsize=figsize)
        
        # pycirclizeをインポート
        from pycirclize import Circos
        from pycirclize.parser import Matrix as PyMatrix
        
        # カスタム色辞書を作成 - セクター用
        sector_colors = {}
        for cell in final_cell_types:
            if cell in color_use:
                sector_colors[cell] = color_use[cell]
            else:
                sector_colors[cell] = (0.7, 0.7, 0.7, 1.0)  # デフォルトはグレー
        
        print(f"[CHORD DEBUG] Final sector_colors: {sector_colors}")
        print(f"[CHORD DEBUG] First 3 sector colors:")
        for i, (cell, color) in enumerate(list(sector_colors.items())[:3]):
            print(f"  {i+1}. {cell}: {color}")
        
        # chord_diagramのlink_kws設定
        link_kws = {
            'direction': 1,
            'alpha': alpha_edge  # リンクの透明度を設定
        }
        
        # リンクの黒い輪郭を表示するかどうか
        if show_edge_border:
            link_kws['ec'] = "black"
            link_kws['lw'] = edge_border_width
        else:
            link_kws['ec'] = "none"
            link_kws['lw'] = 0
        
        # sector_colorをmatplotlibで扱える形式に変換
        converted_sector_colors = {}
        for cell, color in sector_colors.items():
            if isinstance(color, tuple):
                if len(color) == 3:
                    # RGBを16進数カラーコードに変換
                    converted_sector_colors[cell] = mcolors.rgb2hex(color)
                elif len(color) == 4:
                    # RGBAを16進数カラーコードに変換
                    converted_sector_colors[cell] = mcolors.rgb2hex(color[:3])
            else:
                # そのまま使用
                converted_sector_colors[cell] = color
        
        # マトリックスをpycirclizeのMatrix形式に変換
        py_matrix = PyMatrix(matrix_df)
        
        # Chord Diagramを作成（spaceパラメータを設定）
        chord_diagram = Circos.chord_diagram(
            py_matrix,
            space=space,  # 細胞間の空白を設定
            endspace=True,  # 最後の細胞の後にも空白を入れる
            cmap=converted_sector_colors,  # colormap - セクター色を指定
            r_lim=(97, 100),  # 半径の範囲
            link_kws=link_kws,  # リンクのプロパティ（alpha、輪郭の設定を含む）
        )
        
        # プロット
        chord_diagram.plotfig()
        
        # タイトル設定
        plt.title(title, fontsize=16, pad=20)
        plt.tight_layout()
        
        # 図を返す
        return plt.gcf()
        
    except Exception as e:
        print(f"コードダイアグラム生成エラー: {str(e)}")
        print(traceback.format_exc())
        
        # エラー時は空の図を返す
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"コードダイアグラム生成エラー: {str(e)}", 
                ha='center', va='center', fontsize=12, wrap=True)
        ax.axis('off')
        plt.title(title, fontsize=16)
        return fig



def netAnalysis_signalingRole_heatmap(
    net, 
    signaling=None, 
    pattern="outgoing",
    title=None,
    color_heatmap=None,
    thresh=0.05,
    width=10, 
    height=6,
    font_size=10, 
    cluster_rows=True, 
    cluster_cols=True,
    display_numbers=False,
    slot_name="netP",
    sorted_order=None,
    cmap_name=None,
    sources_use=None,
    targets_use=None
):
    """
    Complete redesign of the heatmap function using direct matplotlib control
    to ensure perfect alignment between heatmap and bar plots.
    
    Parameters
    ----------
    net : dict
        CellChat object with network information
    signaling : list or None
        Specific signaling pathways to analyze. If None, analyze all pathways.
    pattern : str
        "outgoing" or "incoming" signaling analysis
    title : str or None
        Title for the plot
    color_heatmap : str or list
        Color mapping for the heatmap
    thresh : float
        Threshold for significance filtering
    width : float
        Width of the plot
    height : float
        Height of the plot
    font_size : float
        Base font size for the plot
    cluster_rows : bool
        Whether to cluster rows
    cluster_cols : bool
        Whether to cluster columns
    display_numbers : bool
        Whether to display the values in the heatmap cells
    slot_name : str
        Slot name where network information is stored
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure with the heatmap
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import matplotlib.gridspec as gridspec
    from scipy.cluster.hierarchy import linkage, leaves_list
    import matplotlib.colors as mcolors
    import matplotlib.cm as cm
    import traceback
    
    try:
        # Check if necessary data exists
        if slot_name not in net:
            raise ValueError(f"Slot '{slot_name}' not found in net object")
        
        # Get cell types
        if 'net' in net and 'dimnames' in net['net'] and len(net['net']['dimnames']) > 0:
            all_cell_types = list(net['net']['dimnames'][0])  # 明示的にリストに変換
        else:
            raise ValueError("Cell type information not found in network data")
        
        # Apply filtering for contribution analysis - pattern-specific filtering
        if sources_use is not None or targets_use is not None:
            if pattern == "outgoing":
                # For outgoing (sources sending), filter by sources_use only
                if sources_use is not None:
                    cell_types = [ct for ct in all_cell_types if ct in sources_use]
                else:
                    cell_types = all_cell_types
            elif pattern == "incoming":
                # For incoming (targets receiving), filter by targets_use only  
                if targets_use is not None:
                    cell_types = [ct for ct in all_cell_types if ct in targets_use]
                else:
                    cell_types = all_cell_types
            else:
                cell_types = all_cell_types
        else:
            cell_types = all_cell_types
        
        # Create mapping from original indices to filtered indices
        cell_type_indices = {ct: i for i, ct in enumerate(all_cell_types)}
        filtered_indices = [cell_type_indices[ct] for ct in cell_types]
        
        # Get pathway information
        if signaling is None:
            if 'pathways' in net[slot_name]:
                pathway_names = list(net[slot_name]['pathways'])
                # 文字列のリストに変換
                pathway_names = [str(p) for p in pathway_names]
            else:
                raise ValueError("No pathway information found in the network object")
        else:
            # Filter to specified pathways
            if 'pathways' in net[slot_name]:
                all_paths = list(net[slot_name]['pathways'])
                # 文字列に変換
                all_paths = [str(p) for p in all_paths]
                signaling = [str(s) for s in signaling]
                
                pathway_names = [p for p in signaling if p in all_paths]
                if not pathway_names:
                    raise ValueError("None of the specified pathways found in the network")
            else:
                pathway_names = [str(s) for s in signaling]
        
        # Prepare data for the heatmap
        if 'prob' in net[slot_name]:
            prob = net[slot_name]['prob']
        elif 'prob' in net['net']:
            # Use general network probability
            prob = net['net']['prob']
        else:
            raise ValueError("Probability matrix not found in the network object")
        
        # Ensure prob is numpy array
        prob = np.array(prob, dtype=float)
        
        # Filter by pathway if signaling is specified
        pathway_indices = []
        if signaling is not None and 'pathways' in net[slot_name]:
            all_paths = [str(p) for p in net[slot_name]['pathways']]
            
            for p in pathway_names:
                if p in all_paths:
                    idx = all_paths.index(p)
                    if prob.ndim == 3 and idx < prob.shape[2]:
                        pathway_indices.append(idx)
            
            if not pathway_indices:
                raise ValueError("None of the specified pathways found or valid")
            
            # Filter the prob array by pathway indices
            if prob.ndim == 3:
                prob = prob[:, :, pathway_indices]
                # Update pathway names to match filtered probabilities
                pathway_names = [all_paths[i] for i in pathway_indices]
        
        # Apply thresholding if pval is available
        if 'pval' in net[slot_name]:
            pval = net[slot_name]['pval']
            # Apply the same filtering as for prob
            if signaling is not None and 'pathways' in net[slot_name] and pval.ndim == 3:
                pval = pval[:, :, pathway_indices]
            
            # Apply threshold
            prob_filtered = prob.copy()
            prob_filtered[pval > thresh] = 0
        else:
            prob_filtered = prob
        
        # Calculate signaling contribution matrix based on pattern
        if prob_filtered.ndim == 3:
            # 3D case - multiple pathways in prob
            num_pathways = prob_filtered.shape[2]
            
            # Create matrix for heatmap [filtered_cell_types × pathways]
            contribution = np.zeros((len(cell_types), num_pathways), dtype=float)
            
            if pattern == "outgoing":
                # For outgoing, sum over each row (sender) for filtered cells only
                for i, cell_idx in enumerate(filtered_indices):
                    for j in range(num_pathways):
                        contribution[i, j] = np.sum(prob_filtered[cell_idx, :, j])
            else:  # incoming
                # For incoming, sum over each column (receiver) for filtered cells only
                for i, cell_idx in enumerate(filtered_indices):
                    for j in range(num_pathways):
                        contribution[i, j] = np.sum(prob_filtered[:, cell_idx, j])
            
            # Ensure pathway_names is the correct length
            if len(pathway_names) > num_pathways:
                pathway_names = pathway_names[:num_pathways]
            elif len(pathway_names) < num_pathways:
                pathway_names = pathway_names + [f"Pathway_{i}" for i in range(len(pathway_names), num_pathways)]
        else:
            # 2D case - only one pathway in prob
            
            # Create matrix for heatmap [filtered_cell_types × 1]
            if pattern == "outgoing":
                # For outgoing, sum over each row (sender) for filtered cells only
                contribution = np.array([np.sum(prob_filtered[cell_idx, :]) for cell_idx in filtered_indices]).reshape(-1, 1)
            else:  # incoming
                # For incoming, sum over each column (receiver) for filtered cells only
                contribution = np.array([np.sum(prob_filtered[:, cell_idx]) for cell_idx in filtered_indices]).reshape(-1, 1)
            
            # If only one pathway, use its name
            if len(pathway_names) == 1:
                pathway_names = pathway_names
            else:
                pathway_names = ["Pathway"]
        
        # cell_types should now match the contribution matrix rows exactly
        # (no adjustment needed since we filtered properly above)
        
        # Create pandas DataFrame for better handling
        df_contribution = pd.DataFrame(
            contribution, 
            index=cell_types,
            columns=pathway_names
        )
        
        # Convert to float explicitly
        df_contribution = df_contribution.astype(float)
        
        # Check for invalid values
        if np.any(np.isnan(df_contribution.values)) or np.any(np.isinf(df_contribution.values)):
            df_contribution = df_contribution.replace([np.inf, -np.inf], 0).fillna(0)

        # Apply sorted_order if provided
        if sorted_order is not None:
            # Get valid cell types from sorted_order (those that are in the data)
            valid_sorted = [cell for cell in sorted_order if cell in df_contribution.index]
            
            # Get remaining cells (those not in sorted_order)
            remaining_cells = [cell for cell in df_contribution.index if cell not in valid_sorted]
            
            # Create new order combining valid_sorted and remaining_cells
            new_order = valid_sorted + remaining_cells
            
            # Reindex dataframe
            df_contribution = df_contribution.reindex(new_order)
        
        # Row-scale the values (normalize each row)
        df_norm = pd.DataFrame(index=df_contribution.index, columns=df_contribution.columns, dtype=float)
        for i, row in enumerate(df_contribution.index):
            row_values = df_contribution.loc[row, :].values
            row_max = np.max(row_values)
            if row_max > 0:
                df_norm.loc[row, :] = row_values / row_max
            else:
                df_norm.loc[row, :] = 0.0
        
        # Calculate row and column sums for bar plots
        row_sums = df_contribution.sum(axis=1)  # Cell type totals
        col_sums = df_contribution.sum(axis=0)  # Pathway totals
        
        # Apply clustering if requested
        if cluster_rows or cluster_cols:
            try:
                if cluster_rows and len(df_norm) > 1:
                    # Use scipy's linkage for hierarchical clustering
                    row_data = df_norm.values.astype(float)
                    row_linkage = linkage(row_data, method='average')
                    row_order = leaves_list(row_linkage)
                    # Reorder rows
                    df_norm = df_norm.iloc[row_order, :]
                    row_sums = row_sums.iloc[row_order]
                
                if cluster_cols and df_norm.shape[1] > 1:
                    # Cluster columns (pathways)
                    col_data = df_norm.values.T.astype(float)
                    col_linkage = linkage(col_data, method='average')
                    col_order = leaves_list(col_linkage)
                    # Reorder columns
                    df_norm = df_norm.iloc[:, col_order]
                    col_sums = col_sums.iloc[col_order]
            except Exception as e:
                print(f"Clustering error: {str(e)}")
                # Continue without clustering
        
        # Now we have our final data frame with correct order
        # Store the final cell and pathway orders
        final_cell_order = df_norm.index.tolist()
        final_pathway_order = df_norm.columns.tolist()
        
        # Transpose the data for our visualization
        df_t = df_norm.T
        
        # Prepare colors for heatmap
        if color_heatmap is None:
            # Default to YlOrRd colormap for heatmap
            cmap_name = "YlOrRd"
        else:
            cmap_name = color_heatmap
            
        try:
            cmap = plt.get_cmap(cmap_name)
        except:
            print(f"Warning: Colormap '{cmap_name}' not found, using YlOrRd instead.")
            cmap = plt.get_cmap("YlOrRd")
            
        cmap = plt.get_cmap(cmap_name)
        
        # ==== Create figure with direct matplotlib control ====
        fig = plt.figure(figsize=(width, height))
        
        # Create layout with precise control
        # Main heatmap takes most space, with margins for bars
        heatmap_left = 0.15   # Left margin for y-labels
        heatmap_bottom = 0.15 # Bottom margin for x-labels
        heatmap_width = 0.65  # Width of heatmap
        heatmap_height = 0.65 # Height of heatmap
        
        # Top bar takes same width as heatmap, small height
        top_bar_height = 0.1
        
        # Right bar takes same height as heatmap, small width
        right_bar_width = 0.1
        
        # Create the axes with precise positioning
        ax_heatmap = fig.add_axes([heatmap_left, heatmap_bottom, heatmap_width, heatmap_height])
        ax_top = fig.add_axes([heatmap_left, heatmap_bottom + heatmap_height + 0.01, heatmap_width, top_bar_height])
        ax_right = fig.add_axes([heatmap_left + heatmap_width + 0.01, heatmap_bottom, right_bar_width, heatmap_height])
        
        # Add colorbar
        cbar_ax = fig.add_axes([heatmap_left + heatmap_width + right_bar_width + 0.05, heatmap_bottom, 0.02, heatmap_height])
        
        # Axes for labels/titles
        title_ax = fig.add_axes([0, 0.9, 1, 0.1])
        title_ax.axis('off')
        
        # Draw the main heatmap
        # Convert dataframe to numpy array for direct plotting
        data = df_t.values
        
        # Create a normalized colormap
        norm = mcolors.Normalize(vmin=0, vmax=1)
        
        # Draw heatmap as a pcolormesh
        x_positions = np.arange(len(final_cell_order) + 1)
        y_positions = np.arange(len(final_pathway_order) + 1)
        
        mesh = ax_heatmap.pcolormesh(x_positions, y_positions, data, cmap=cmap, norm=norm)
        
        # Add text annotations if requested
        if display_numbers:
            for i in range(len(final_pathway_order)):
                for j in range(len(final_cell_order)):
                    ax_heatmap.text(j + 0.5, i + 0.5, f"{data[i, j]:.2f}", 
                                   ha="center", va="center", 
                                   fontsize=font_size-2)
        
        # Set labels on the axes
        ax_heatmap.set_xticks(np.arange(len(final_cell_order)) + 0.5)
        ax_heatmap.set_yticks(np.arange(len(final_pathway_order)) + 0.5)
        
        ax_heatmap.set_xticklabels(final_cell_order, rotation=45, ha='right', fontsize=font_size-1)
        ax_heatmap.set_yticklabels(final_pathway_order, fontsize=font_size-1)
        
        # Add a colorbar
        cb = plt.colorbar(mesh, cax=cbar_ax)
        cb.set_label("Relative Signaling Strength", fontsize=font_size-1)
        
        # Draw the top bar (cell types)
        cell_sums = row_sums.reindex(final_cell_order).values
        bars_top = ax_top.bar(np.arange(len(final_cell_order)) + 0.5, cell_sums, 
                             width=0.8, align='center', color='skyblue')
        
        # Make the top axis match the heatmap exactly
        ax_top.set_xlim(0, len(final_cell_order))
        ax_top.set_xticks([])
        ax_top.set_title(f"{pattern.capitalize()} Signaling Strength", fontsize=font_size)
        ax_top.spines['top'].set_visible(False)
        ax_top.spines['right'].set_visible(False)
        
        # Draw the right bar (pathways)
        pathway_sums = col_sums.reindex(final_pathway_order).values
        bars_right = ax_right.barh(np.arange(len(final_pathway_order)) + 0.5, pathway_sums, 
                                 height=0.8, align='center', color='lightgrey')
        
        # Make the right axis match the heatmap exactly
        ax_right.set_ylim(0, len(final_pathway_order))
        ax_right.set_yticks([])
        ax_right.set_title("Pathway\nStrength", fontsize=font_size, rotation=-90, 
                         x=1.4, y=0.5, ha='center', va='center')
        ax_right.spines['top'].set_visible(False)
        ax_right.spines['right'].set_visible(False)
        
        # Add overall title
        if title:
            title_text = title
        else:
            pattern_title = "Outgoing" if pattern == "outgoing" else "Incoming"
            signaling_title = f"Signaling Pathways: {', '.join(final_pathway_order[:3])}" + (f" + {len(final_pathway_order)-3} more" if len(final_pathway_order) > 3 else "")
            title_text = f"{pattern_title} Signaling Contribution - {signaling_title}"
            
        title_ax.text(0.5, 0.5, title_text, ha='center', va='center', fontsize=font_size+2)
        
        return fig
    
    except Exception as e:
        print(f"Error in netAnalysis_signalingRole_heatmap_transpose: {str(e)}")
        print(traceback.format_exc())
        
        # Create an error figure
        fig, ax = plt.subplots(figsize=(width, height))
        ax.text(0.5, 0.5, f"Error generating heatmap:\n{str(e)}", 
               ha='center', va='center', fontsize=12, wrap=True)
        ax.axis('off')
        return fig



def plot_network_centrality(network_summary, figsize=(12, 8)):
    """
    各細胞タイプのネットワーク中心性を描画
    
    Parameters
    ----------
    network_summary : dict
        ネットワーク集計結果
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    try:
        if 'network_centrality' not in network_summary or network_summary['network_centrality'].empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なネットワーク中心性データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Network Centrality Measures')
            return fig
            
        centrality_df = network_summary['network_centrality']
        
        if len(centrality_df) == 0:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なネットワーク中心性データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Network Centrality Measures')
            return fig
        
        # 使用可能な列を確認
        available_cols = ['cell_type']
        expected_cols = {
            'degree': ['degree'],
            'in_degree': ['in_degree', 'indeg'],
            'out_degree': ['out_degree', 'outdeg'],
            'betweenness': ['betweenness'],
            'eigenvector': ['eigenvector', 'eigen'],
            'hub': ['hub'],
            'authority': ['authority'],
            'page_rank': ['page_rank'],
            'flowbet': ['flowbet'],
            'info': ['info']
        }
        
        # データフレーム内の列に基づいて使用可能な中心性指標を決定
        col_mapping = {}
        for display_name, possible_cols in expected_cols.items():
            for col in possible_cols:
                if col in centrality_df.columns:
                    col_mapping[col] = display_name
                    available_cols.append(col)
                    break
        
        if len(available_cols) <= 1:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "中心性指標が見つかりません", ha='center', va='center')
            ax.axis('off')
            plt.title('Network Centrality Measures')
            return fig
        
        # データを整形
        centrality_df_melt = pd.melt(
            centrality_df, 
            id_vars=['cell_type'], 
            value_vars=[col for col in available_cols if col != 'cell_type'],
            var_name='centrality_type',
            value_name='value'
        )
        
        # 列名をマッピングで置換
        centrality_df_melt['centrality_type'] = centrality_df_melt['centrality_type'].map(
            lambda x: col_mapping.get(x, x)
        )
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # グループごとの棒グラフを描画
        sns.barplot(
            data=centrality_df_melt,
            x='cell_type',
            y='value',
            hue='centrality_type',
            ax=ax
        )
        
        plt.title('Network Centrality Measures for Cell Types')
        plt.xticks(rotation=45, ha='right')
        plt.legend(title='Centrality Type', bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        
        return fig
    except Exception as e:
        import traceback
        print(f"ネットワーク中心性プロット作成エラー: {str(e)}")
        print(traceback.format_exc())
        
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"ネットワーク中心性プロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Network Centrality Measures')
        return fig

def plot_lr_contribution(network_summary, top_n=20, figsize=(12, 8)):
    """
    リガンド-レセプターペアの寄与度を棒グラフで描画（直接的な解決策）
    
    Parameters
    ----------
    network_summary : dict
        ネットワーク集計結果
    top_n : int
        表示する上位ペアの数
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import traceback
    
    try:
        # 寄与度データを取得
        if 'lr_contribution' not in network_summary or np.all(network_summary['lr_contribution'] == 0):
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なLR寄与度データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Contributions')
            return fig
        
        # 結果テーブルからリガンド-レセプター情報を取得
        # これは必ず存在する想定
        if 'results' not in network_summary or network_summary['results'].empty:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "結果テーブルが見つかりません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Contributions')
            return fig
        
        # LR寄与度とインタラクション名のマッピング
        lr_values = network_summary['lr_contribution']
        interaction_names = []
        
        # dimnames[2]からインタラクション名を取得
        if 'net' in network_summary and 'dimnames' in network_summary['net'] and len(network_summary['net']['dimnames']) > 2:
            interaction_names = list(network_summary['net']['dimnames'][2])
        
        # この時点でinteraction_namesが空の場合はダミー名を生成
        if not interaction_names:
            interaction_names = [f"LR_{i}" for i in range(len(lr_values))]
        
        # 各インタラクション名に対応するリガンド-レセプターペア名を見つける
        # 結果テーブルから直接マッピング
        results_df = network_summary['results']
        
        # 必要な列があるか確認
        required_cols = ['interaction_name', 'ligand', 'receptor']
        if not all(col in results_df.columns for col in required_cols):
            # 必要な列がない場合、interaction_name列だけで処理
            display_names = interaction_names
        else:
            # リガンド-レセプターの表示名を作成
            name_map = {}
            for _, row in results_df.iterrows():
                # interaction_nameをキーとして、ligand-receptor形式のテキストをvalueとする
                if pd.notna(row['interaction_name']) and pd.notna(row['ligand']) and pd.notna(row['receptor']):
                    name_map[row['interaction_name']] = f"{row['ligand']}-{row['receptor']}"
            
            # display_namesリストを作成（存在しないものはそのままinteraction_nameを使用）
            display_names = [name_map.get(name, name) for name in interaction_names]
        
        # 寄与度と表示名をまとめたデータフレームを作成
        data = []
        for i, (value, name) in enumerate(zip(lr_values, display_names)):
            if value > 0:  # 寄与度が正の場合のみ含める
                data.append({
                    'index': i,
                    'name': name,
                    'value': value
                })
        
        # データフレームに変換
        df = pd.DataFrame(data)
        
        # 寄与度で降順ソート（明示的にascending=Falseを指定）
        df = df.sort_values('value', ascending=False).head(top_n)
        
        # 図の高さを動的に調整
        height = max(figsize[1], 0.4 * len(df) + 3)
        adjusted_figsize = (figsize[0], height)
        
        # 図を作成
        fig, ax = plt.subplots(figsize=adjusted_figsize)
        
        # カラーマップ（viridisの場合、高い値ほど明るい色に）
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(df)))
        
        # 横向きの棒グラフで描画（値の大きい順）
        bars = ax.barh(
            range(len(df)),
            df['value'].values,
            color=colors
        )
        
        # Y軸のラベルを明示的に設定
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['name'].values)  # ここで名前を明示的に設定
        
        # 棒グラフに数値を表示
        for i, bar in enumerate(bars):
            width = bar.get_width()
            label_x_pos = width + 0.05
            ax.text(label_x_pos, i, f"{width:.2f}", va='center')
        
        # 軸ラベルなどの設定
        ax.set_xlabel('Contribution Strength')
        ax.set_title(f'Top {min(top_n, len(df))} Ligand-Receptor Contributions')
        
        # Y軸のラベル設定を調整（長い名前に対応）
        plt.yticks(fontsize=10)
        plt.tight_layout()
        
        # グリッド線の追加
        ax.grid(axis='x', linestyle='--', alpha=0.3)
        
        # スタイル調整
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_axisbelow(True)  # グリッド線を背景に
        
        return fig
    except Exception as e:
        print(f"LR寄与度プロット作成エラー: {str(e)}")
        print(traceback.format_exc())
        
        # エラーの場合は説明付きの空のプロットを返す
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"LR寄与度プロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Ligand-Receptor Contributions')
        return fig

def plot_lr_contribution(network_summary, top_n=20, figsize=(12, 8)):
    """
    リガンド-レセプターペアの寄与度を棒グラフで描画（改良版）
    
    Parameters
    ----------
    network_summary : dict
        ネットワーク集計結果
    top_n : int
        表示する上位ペアの数
    figsize : tuple
        図のサイズ
        
    Returns
    -------
    matplotlib.figure.Figure
        プロットのFigureオブジェクト
    """

    
    try:
        # 寄与度データを取得
        if 'lr_contribution' not in network_summary or np.all(network_summary['lr_contribution'] == 0):
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なLR寄与度データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Contributions')
            return fig
        
        # リガンド-レセプターペアの情報取得
        lr_pairs = {}
        
        # 1. 直接渡されたLR名辞書を優先的に使用
        if 'lr_names' in network_summary and isinstance(network_summary['lr_names'], dict):
            lr_pairs = network_summary['lr_names']
        
        # 2. resultsからLR情報を取得（バックアップ）
        elif 'results' in network_summary and not network_summary['results'].empty:
            results = network_summary['results']
            # 結果からLR情報を抽出
            if all(col in results.columns for col in ['interaction_name', 'ligand', 'receptor']):
                for idx, row in results.iterrows():
                    lr_pairs[row['interaction_name']] = f"{row['ligand']}-{row['receptor']}"
        
        # インタラクション名の取得
        interaction_names = []
        
        # 1. 直接渡されたインタラクション名リストを優先
        if 'interaction_names' in network_summary and isinstance(network_summary['interaction_names'], (list, np.ndarray)):
            interaction_names = network_summary['interaction_names']
        # 2. dimnames[2]からインタラクション名を取得（バックアップ）
        elif 'net' in network_summary and 'dimnames' in network_summary['net'] and len(network_summary['net']['dimnames']) > 2:
            interaction_names = network_summary['net']['dimnames'][2]
        
        # LR寄与度データを取得し、DataFrameを構築
        lr_contrib_data = []
        
        for i, val in enumerate(network_summary['lr_contribution']):
            if val > 0:
                # インタラクション名の取得
                interaction_name = interaction_names[i] if i < len(interaction_names) else f"Interaction_{i+1}"
                
                # LR名の取得
                # 1. lr_pairs辞書から取得
                if interaction_name in lr_pairs:
                    lr_name = lr_pairs[interaction_name]
                # 2. インタラクション名そのものを使用
                else:
                    lr_name = interaction_name
                    
                # pathway情報の取得（可能な場合）
                pathway_info = ""
                if 'results' in network_summary and not network_summary['results'].empty:
                    if all(col in network_summary['results'].columns for col in ['interaction_name', 'pathway_name']):
                        matching_rows = network_summary['results'][
                            network_summary['results']['interaction_name'] == interaction_name
                        ]
                        if not matching_rows.empty and 'pathway_name' in matching_rows.iloc[0]:
                            pathway_info = matching_rows.iloc[0]['pathway_name']
                
                # 寄与度データに追加
                lr_contrib_data.append({
                    'interaction_name': interaction_name,
                    'lr_name': lr_name,
                    'pathway': pathway_info,
                    'interaction_prob_normalized': val
                })
        
        if not lr_contrib_data:
            fig, ax = plt.subplots(figsize=figsize)
            ax.text(0.5, 0.5, "有効なLR寄与度データがありません", ha='center', va='center')
            ax.axis('off')
            plt.title('Ligand-Receptor Contributions')
            return fig
        
        # DataFrameを作成して降順ソート
        lr_contrib_df = pd.DataFrame(lr_contrib_data)
        lr_contrib_df = lr_contrib_df.sort_values('interaction_prob_normalized', ascending=True).tail(top_n)

        
        # 表示用ラベルの作成
        display_labels = lr_contrib_df['lr_name'].tolist()
        
        # 図の高さを動的に調整（ペア数に比例）
        height = max(figsize[1], 0.3 * len(display_labels) + 3)
        adjusted_figsize = (figsize[0], height)
        
        # 図を作成
        fig, ax = plt.subplots(figsize=adjusted_figsize)
        
        # カラーマップ（viridisの場合、高い値ほど明るい色に）
        colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(lr_contrib_df)))
        
        # 横向きの棒グラフで描画（値の大きい順に表示）
        bars = ax.barh(
            range(len(display_labels)),  # Y軸の位置を逆順に
            lr_contrib_df['interaction_prob_normalized'].values,
            color=colors
        )
        
        # Y軸のラベルを設定
        ax.set_yticks(range(len(display_labels)))
        ax.set_yticklabels(display_labels)
        
        # 棒グラフに数値を表示
        for i, bar in enumerate(bars):
            width = bar.get_width()
            label_x_pos = width + 0.05
            ax.text(label_x_pos, bar.get_y() + bar.get_height()/2, f"{width:.2f}", 
                    va='center', fontsize=14)
        
        # 軸ラベルなどの設定
        ax.set_xlabel('Contribution Strength', fontsize=18)
        ax.set_title(f'Top {min(top_n, len(lr_contrib_df))} Ligand-Receptor Contributions', fontsize=20)
        
        # Y軸のラベル設定を調整（長い名前に対応）
        plt.yticks(fontsize=18)
        plt.xticks(fontsize=16)
        plt.tight_layout()
        
        # スタイル調整
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.set_axisbelow(True)  # グリッド線を背景に
        
        return fig
    except Exception as e:
        print(f"LR寄与度プロット作成エラー: {str(e)}")
        print(traceback.format_exc())
        
        # エラーの場合は説明付きの空のプロットを返す
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f"LR寄与度プロット作成エラー: {str(e)}", ha='center', va='center', wrap=True)
        ax.axis('off')
        plt.title('Ligand-Receptor Contributions')
        return fig


def plotGeneExpression(result, signaling, cellchatdb=None, group_by=None, cmap_name='tab10'):
    """
    遺伝子発現をScanpyのビオリンプロットで可視化する関数（改良版）
    すべての遺伝子を1つのプロットに統合し、各ビオリンプロットを縦に並べる
    
    Parameters
    ----------
    result : dict
        CellChat解析結果
    signaling : str
        シグナリングパスウェイ名
    cellchatdb : dict, optional
        CellChatDBデータ（interaction, complex, cofactorなどを含む辞書）
    group_by : str, optional
        細胞をグループ化するためのobsの列名
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        プロット図
    """
    import scanpy as sc
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    import re
    from matplotlib.backends.backend_pdf import PdfPages
    
    # グループ化の設定
    if group_by is None:
        # 利用可能なカテゴリカル列を探す
        cat_cols = [col for col in result['adata'].obs.columns if pd.api.types.is_categorical_dtype(result['adata'].obs[col])]
        if cat_cols:
            group_by = cat_cols[0]  # 最初のカテゴリカル列を使用
        else:
            # カテゴリカル列がない場合
            group_by = result['adata'].obs.columns[0]  # 最初の列を使用
            print(f"Warning: Using non-categorical column '{group_by}' for grouping")
    
    print(f"Using column '{group_by}' for cell grouping")
    
    # パスウェイ関連の遺伝子を取得する方法を選択
    pathway_genes = []
    
    # 方法1: CellChatDBから直接取得（より信頼性が高い）
    if cellchatdb is not None and 'interaction' in cellchatdb:
        print(f"Searching for genes related to pathway '{signaling}' in CellChatDB...")
        
        # 部分一致も含めてパスウェイを検索
        pathway_pattern = re.compile(f".*{re.escape(signaling)}.*", re.IGNORECASE)
        matching_interactions = cellchatdb['interaction'][
            cellchatdb['interaction']['pathway_name'].apply(
                lambda x: bool(pathway_pattern.match(str(x))) if pd.notna(x) else False
            )
        ]
        
        if len(matching_interactions) == 0:
            print(f"Warning: No interactions found for pathway '{signaling}' in CellChatDB")
        else:
            print(f"Found {len(matching_interactions)} interactions for pathway pattern '{signaling}'")
            
            # リガンドとレセプターを抽出
            ligands = set()
            receptors = set()
            
            for _, row in matching_interactions.iterrows():
                if pd.notna(row['ligand']):
                    ligands.add(str(row['ligand']))
                if pd.notna(row['receptor']):
                    receptors.add(str(row['receptor']))
            
            # 複合体を分解 (例: "COL1A1_COL1A2" → ["COL1A1", "COL1A2"])
            all_genes = set()
            for gene in ligands.union(receptors):
                if '_' in gene:
                    all_genes.update(gene.split('_'))
                else:
                    all_genes.add(gene)
            
            pathway_genes = list(all_genes)
            print(f"Extracted {len(pathway_genes)} genes from CellChatDB: {', '.join(pathway_genes[:5])}...")
    
    # 方法2: 結果からパスウェイ関連遺伝子を抽出（バックアップ）
    if not pathway_genes and 'results' in result and isinstance(result['results'], pd.DataFrame):
        print(f"Searching for genes related to pathway '{signaling}' in results...")
        
        # 完全一致と部分一致の両方を試す
        matching_results = result['results'][
            result['results']['interaction_name'].str.contains(signaling, case=False, na=False)
        ]
        
        if 'pathway_name' in result['results'].columns:
            # pathway_name列でも検索
            pathway_matches = result['results'][
                result['results']['pathway_name'].str.contains(signaling, case=False, na=False)
            ]
            matching_results = pd.concat([matching_results, pathway_matches]).drop_duplicates()
        
        if len(matching_results) == 0:
            print(f"Warning: No results found for pathway '{signaling}'")
        else:
            print(f"Found {len(matching_results)} results for pathway pattern '{signaling}'")
            
            # リガンドとレセプターを抽出
            ligands = set()
            receptors = set()
            
            if 'ligand' in matching_results.columns and 'receptor' in matching_results.columns:
                for _, row in matching_results.iterrows():
                    if pd.notna(row['ligand']):
                        ligands.add(str(row['ligand']))
                    if pd.notna(row['receptor']):
                        receptors.add(str(row['receptor']))
                
                # 複合体を分解
                all_genes = set()
                for gene in ligands.union(receptors):
                    if '_' in gene:
                        all_genes.update(gene.split('_'))
                    else:
                        all_genes.add(gene)
                
                pathway_genes.extend(list(all_genes))
                print(f"Extracted {len(all_genes)} genes from results")
    
    # 方法3: 最後の手段 - パスウェイ名を直接遺伝子名として試す
    if not pathway_genes:
        print(f"Warning: No genes found for pathway '{signaling}'. Trying pathway name as gene...")
        pathway_genes = [signaling]
    
    # 大文字小文字を無視した遺伝子マッチング
    existing_genes = []
    
    # デバッグ: 全ての遺伝子名を小文字に変換して比較
    adata_genes_lower = {gene.lower(): gene for gene in result['adata'].var_names}
    
    for gene in pathway_genes:
        # 完全一致と部分一致の両方を試す
        if gene.lower() in adata_genes_lower:
            # 完全一致
            existing_genes.append(adata_genes_lower[gene.lower()])
        else:
            # 部分一致を試す
            matching_genes = [
                adata_genes_lower[g_lower] 
                for g_lower in adata_genes_lower
                if gene.lower() in g_lower or g_lower in gene.lower()
            ]
            existing_genes.extend(matching_genes[:5])  # 最大5つまで
    
    # 重複排除
    existing_genes = list(set(existing_genes))
    
    print(f"Found {len(existing_genes)} matching genes in the dataset")
    print(f"Matching genes: {', '.join(existing_genes[:10])}" + ("..." if len(existing_genes) > 10 else ""))
    
    if not existing_genes:
        raise ValueError(f"パスウェイ {signaling} の遺伝子がデータセットに見つかりません。")
    
    # 最大表示する遺伝子数（必要に応じて調整）
    max_genes = 20
    genes_to_plot = existing_genes[:max_genes]
    n_genes = len(genes_to_plot)
    
    # 全ての遺伝子に対して1つの図にプロット
    fig = plt.figure(figsize=(12, n_genes * 1.5))
    
    # グリッドスペックを設定（最後の行だけラベルを表示）
    gs = fig.add_gridspec(n_genes, 1, hspace=0.05)
    
    # カテゴリカルデータから一意なカテゴリを取得
    categories = list(result['adata'].obs[group_by].cat.categories)
    
    # プロットのループ
    for i, gene in enumerate(genes_to_plot):
        ax = fig.add_subplot(gs[i, 0])
        
        # 個別のビオリンプロットを作成
        sc.pl.violin(
            result['adata'], 
            keys=gene, 
            groupby=group_by, 
            ax=ax,
            show=False,
            palette=cmap_name
        )
        
        # 枠線を非表示に
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        
        # 最後の行以外はX軸のラベルを非表示
        if i < n_genes - 1:
            ax.set_xlabel('')
            ax.set_xticklabels([])
            ax.set_xticks([])
        else:
            # 最後の行はX軸のラベルを回転して表示
            ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
        
        # 各プロットには遺伝子名だけを表示
        ax.set_title('')
        
        # Y軸ラベルを遺伝子名に設定
        ax.set_ylabel(gene, rotation=0, ha='right', va='center', fontsize=10)
    
    # 全体のタイトルを設定
    plt.suptitle(f"{signaling} Pathway Gene Expression", y=1.01, fontsize=16)
    
    # 余白を調整
    plt.tight_layout()
    
    return fig



def netAnalysis_signalingChanges_scatter(
    merged_cellchat, 
    idents_use=None, 
    signaling_exclude=None, 
    signaling_include=None,
    color_use="teal", 
    use_color_by_pattern=False,  # 新しいパラメータ: パターンごとに色分けするかどうか
    color_palette="Set1",  # 新しいパラメータ: 使用するカラーパレット名
    title=None,
    font_size=12,
    dot_size=3,
    max_label=15,
    use_arrows=True,
    figsize=(8, 6)
):
    """
    Identify the signaling changes of specific cell populations
    
    Parameters
    ----------
    merged_cellchat : dict
        Merged CellChat object containing two datasets
    idents_use : str
        Cell population to analyze
    signaling_exclude : list or str, optional
        Signaling pathways to exclude
    signaling_include : list or str, optional
        Signaling pathways to include (takes precedence over exclude)
    color_use : str or dict, optional
        Color for the scatter points (when use_color_by_pattern=False)
    use_color_by_pattern : bool, optional
        Whether to use different colors for each pattern
    color_palette : str, optional
        Name of the color palette to use (e.g., "Set1", "Tab10", "Dark2")
    title : str, optional
        Title for the plot
    font_size : int, optional
        Font size for the text
    dot_size : int, optional
        Size of the dots
    max_label : int, optional
        Maximum number of pathways to label
    use_arrows : bool, optional
        Whether to use arrows to point to data points instead of text boxes
    figsize : tuple, optional
        Figure size (width, height) in inches
        
    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure with signaling changes scatter plot
    """
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Check if merged object has the right structure
    if not isinstance(merged_cellchat, dict) or 'group1' not in merged_cellchat or 'group2' not in merged_cellchat:
        raise ValueError("Merged CellChat object must contain 'group1' and 'group2' entries")
    
    # Extract data for each group
    result1 = merged_cellchat['group1']['result']
    result2 = merged_cellchat['group2']['result']
    group1_name = merged_cellchat['group1']['name']
    group2_name = merged_cellchat['group2']['name']
    
    # Check if netP exists in both results
    if 'netP' not in result1 or 'netP' not in result2:
        raise ValueError("Network pathway data not found in both datasets")
    
    # Get all pathways from both datasets
    pathways1 = result1['netP'].get('pathways', [])
    pathways2 = result2['netP'].get('pathways', [])
    
    # Convert to list if they are numpy arrays
    if isinstance(pathways1, np.ndarray):
        pathways1 = list(pathways1)
    if isinstance(pathways2, np.ndarray):
        pathways2 = list(pathways2)
    
    # Combine unique pathways
    all_pathways = list(set(pathways1 + pathways2))
    
    # Filter pathways if specified
    if signaling_include is not None:
        if isinstance(signaling_include, str):
            signaling_include = [signaling_include]
        all_pathways = [p for p in all_pathways if p in signaling_include]
    elif signaling_exclude is not None:
        if isinstance(signaling_exclude, str):
            signaling_exclude = [signaling_exclude]
        all_pathways = [p for p in all_pathways if p not in signaling_exclude]
    
    # Get cell types and check if idents_use is valid
    if 'net' not in result1 or 'dimnames' not in result1['net'] or len(result1['net']['dimnames']) == 0:
        raise ValueError("Cell type information not found in first dataset")
    cell_types1 = result1['net']['dimnames'][0]
    
    if 'net' not in result2 or 'dimnames' not in result2['net'] or len(result2['net']['dimnames']) == 0:
        raise ValueError("Cell type information not found in second dataset")
    cell_types2 = result2['net']['dimnames'][0]
    
    # Check if specified cell population exists in both datasets
    if idents_use is None:
        raise ValueError("Must specify a cell population to analyze")
    
    if idents_use not in cell_types1 or idents_use not in cell_types2:
        available_cells = list(set(cell_types1).intersection(set(cell_types2)))
        raise ValueError(
            f"Cell population '{idents_use}' not found in both datasets. Available cell populations: {', '.join(available_cells)}"
        )
    
    # Find the index of the specified cell type in each dataset
    cell_idx1 = list(cell_types1).index(idents_use)
    cell_idx2 = list(cell_types2).index(idents_use)
    
    # Calculate differential outgoing and incoming signaling
    diff_data = []
    
    for pathway in all_pathways:
        # Try to find the pathway in each dataset
        pathway_idx1 = pathways1.index(pathway) if pathway in pathways1 else -1
        pathway_idx2 = pathways2.index(pathway) if pathway in pathways2 else -1
        
        # Skip if pathway not found in either dataset
        if pathway_idx1 == -1 and pathway_idx2 == -1:
            continue
        
        # Get probability matrices
        if pathway_idx1 != -1 and 'prob' in result1['netP']:
            prob1 = result1['netP']['prob'][:, :, pathway_idx1]
        else:
            prob1 = np.zeros((len(cell_types1), len(cell_types1)))
        
        if pathway_idx2 != -1 and 'prob' in result2['netP']:
            prob2 = result2['netP']['prob'][:, :, pathway_idx2]
        else:
            prob2 = np.zeros((len(cell_types2), len(cell_types2)))
        
        # Calculate outgoing strength for specified cell
        outgoing1 = np.sum(prob1[cell_idx1, :])
        outgoing2 = np.sum(prob2[cell_idx2, :])
        
        # Calculate incoming strength for specified cell
        incoming1 = np.sum(prob1[:, cell_idx1])
        incoming2 = np.sum(prob2[:, cell_idx2])
        
        # Calculate differences (Group2 - Group1, e.g., CTRL - cKO)
        outgoing_diff = outgoing2 - outgoing1
        incoming_diff = incoming2 - incoming1
        
        # Add to data if there's a difference
        if outgoing_diff != 0 or incoming_diff != 0:
            # Determine the pattern (outgoing, incoming, or both)
            if outgoing_diff > 0 and incoming_diff > 0:
                pattern = "Incoming & Outgoing specific"
                shape = "diamond"
            elif outgoing_diff > 0:
                pattern = "Outgoing specific"
                shape = "triangle-up"
            elif incoming_diff > 0:
                pattern = "Incoming specific"
                shape = "square"
            else:
                pattern = "Other"
                shape = "circle"
            
            diff_data.append({
                "pathway": pathway,
                "outgoing_diff": outgoing_diff,
                "incoming_diff": incoming_diff,
                "pattern": pattern,
                "shape": shape
            })
    
    # Create DataFrame
    df = pd.DataFrame(diff_data)
    
    # To avoid issues with saving figures, explicitly set backend to Agg
    # This line might be important for some environments
    import matplotlib
    matplotlib.use('Agg')
    
    # Create figure with specified size and a conservative DPI
    fig, ax = plt.subplots(figsize=figsize, dpi=100)
    
    # Early return if no data
    if len(df) == 0:
        ax.text(0.5, 0.5, f"No signaling changes found for {idents_use}", 
               ha='center', va='center', fontsize=font_size+2)
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Differential outgoing interaction strength", fontsize=font_size)
        ax.set_ylabel("Differential incoming interaction strength", fontsize=font_size)
        if title is None:
            title = f"Signaling changes of {idents_use} ({group2_name} vs. {group1_name})"
        ax.set_title(title, fontsize=font_size+2)
        ax.axis('on')
        ax.grid(True, linestyle='--', alpha=0.3)
        return fig
    
    # Safety check for extreme values - clip to reasonable range if necessary
    max_val = df[['outgoing_diff', 'incoming_diff']].abs().max().max()
    if max_val > 1000:  # If values are extremely large
        # Log a warning
        print(f"Warning: Extreme values detected (max: {max_val}). Clipping to ±1000.")
        # Clip values to prevent scaling issues
        df['outgoing_diff'] = df['outgoing_diff'].clip(-1000, 1000)
        df['incoming_diff'] = df['incoming_diff'].clip(-1000, 1000)
    
    # Constrain dot size to reasonable values
    dot_size = min(max(1, dot_size), 10)  # Between 1 and 10
    
    # Calculate absolute values for magnitude
    df['outgoing_abs'] = df['outgoing_diff'].abs()
    df['incoming_abs'] = df['incoming_diff'].abs()
    df['magnitude'] = df['outgoing_abs'] + df['incoming_abs']
    
    # Plot data with appropriate markers
    markers = {'Incoming specific': 's', 
               'Outgoing specific': '^', 
               'Incoming & Outgoing specific': 'D',
               'Other': 'o'}
    
    scatter_plots = {}
    
    # Prepare colors for each pattern
    if use_color_by_pattern:
        # Get the selected color palette
        import matplotlib.cm as cm
        from matplotlib.colors import ListedColormap, to_rgba
        
        # Define pattern order for consistent color assignment
        pattern_order = ["Incoming specific", "Outgoing specific", "Incoming & Outgoing specific", "Other"]
        
        # Get available patterns in the data
        available_patterns = [p for p in pattern_order if p in df['pattern'].unique()]
        
        # Get colors from specified palette
        try:
            # Try to get the palette as a colormap
            cmap = plt.get_cmap(color_palette)
            colors = [cmap(i / max(1, len(available_patterns) - 1)) for i in range(len(available_patterns))]
        except:
            # Fallback to Set1 if palette not found
            cmap = plt.get_cmap('Set1')
            colors = [cmap(i / max(1, len(available_patterns) - 1)) for i in range(len(available_patterns))]
        
        # Create color mapping
        pattern_colors = {pattern: colors[i] for i, pattern in enumerate(available_patterns)}
        
        # Plot each pattern with its assigned color
        for pattern in available_patterns:
            subset = df[df['pattern'] == pattern]
            if len(subset) == 0:
                continue
                
            marker = markers[pattern]
            color = pattern_colors[pattern]
            
            scatter = ax.scatter(
                subset['outgoing_diff'], 
                subset['incoming_diff'],
                s=dot_size**2 * 10,  # Scale dot size
                c=[color],  # Use assigned pattern color
                marker=marker,
                alpha=0.8,
                label=pattern,
                edgecolors='gray',
                linewidths=0.5
            )
            scatter_plots[pattern] = scatter
    else:
        # Original single color plotting
        for pattern in df['pattern'].unique():
            subset = df[df['pattern'] == pattern]
            marker = markers[pattern]
            
            if pattern == "Incoming & Outgoing specific":
                subset_color = 'teal' if color_use == 'teal' else color_use
            else:
                subset_color = color_use
                
            scatter = ax.scatter(
                subset['outgoing_diff'], 
                subset['incoming_diff'],
                s=dot_size**2 * 10,  # Scale dot size
                c=subset_color,
                marker=marker,
                alpha=0.7,
                label=pattern,
                edgecolors='gray',
                linewidths=0.5
            )
            scatter_plots[pattern] = scatter
    
    # Add vertical and horizontal reference lines
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    ax.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    
    # Sort by magnitude to prioritize points with largest values for labeling
    df_sorted = df.sort_values('magnitude', ascending=False)
    
    # Limit number of labels to prevent overcrowding
    max_label = min(max_label, len(df))
    
    # Add pathway labels for top signaling changes
    if use_arrows and len(df_sorted) > 0:
        # Identify points to label (top max_label by magnitude)
        label_data = df_sorted.head(max_label).copy()
        
        if len(label_data) > 0:
            # Get the plot boundaries for text placement
            x_min, x_max = ax.get_xlim()
            y_min, y_max = ax.get_ylim()
            x_range = x_max - x_min
            y_range = y_max - y_min
            
            # For each label point, determine placement based on quadrant
            for i, row in label_data.iterrows():
                x, y = row['outgoing_diff'], row['incoming_diff']
                text = row['pathway']
                
                # Determine quadrant and set text position
                if x >= 0 and y >= 0:  # Top right
                    text_x = x + 0.05 * x_range
                    text_y = y + 0.05 * y_range
                    ha = 'left'
                    va = 'bottom'
                elif x < 0 and y >= 0:  # Top left
                    text_x = x - 0.05 * x_range
                    text_y = y + 0.05 * y_range
                    ha = 'right'
                    va = 'bottom'
                elif x < 0 and y < 0:  # Bottom left
                    text_x = x - 0.05 * x_range
                    text_y = y - 0.05 * y_range
                    ha = 'right'
                    va = 'top'
                else:  # Bottom right
                    text_x = x + 0.05 * x_range
                    text_y = y - 0.05 * y_range
                    ha = 'left'
                    va = 'top'
                
                # Add arrow line
                ax.plot([x, text_x], [y, text_y], 'k-', linewidth=0.5, alpha=0.5)
                
                # Add label text
                ax.text(text_x, text_y, text, fontsize=font_size-2,
                        ha=ha, va=va, alpha=0.9)
    else:
        # Traditional text display without arrows
        for i, row in df_sorted.head(max_label).iterrows():
            x_offset = 0.02 * ax.get_xlim()[1] if row['outgoing_diff'] >= 0 else -0.02 * ax.get_xlim()[0]
            y_offset = 0.02 * ax.get_ylim()[1] if row['incoming_diff'] >= 0 else -0.02 * ax.get_ylim()[0]
            
            ha = 'left' if row['outgoing_diff'] >= 0 else 'right'
            va = 'bottom' if row['incoming_diff'] >= 0 else 'top'
            
            ax.text(
                row['outgoing_diff'] + x_offset,
                row['incoming_diff'] + y_offset,
                row['pathway'],
                fontsize=font_size-2,
                ha=ha,
                va=va
            )
    
    # Set labels and title
    ax.set_xlabel("Differential outgoing interaction strength", fontsize=font_size)
    ax.set_ylabel("Differential incoming interaction strength", fontsize=font_size)
    
    if title is None:
        title = f"Signaling changes of {idents_use} ({group2_name} vs. {group1_name})"
    ax.set_title(title, fontsize=font_size+2)
    
    # Add legend
    ax.legend(title="Patterns", fontsize=font_size-2)
    
    # Add grid and style
    ax.grid(True, linestyle='--', alpha=0.3)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    # Auto-scale axes with a small margin
    ax.margins(0.15)
    
    # Apply tight layout to prevent clipping
    plt.tight_layout()
    
    return fig