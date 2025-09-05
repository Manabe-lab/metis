import streamlit as st
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
import numpy as np
import pandas as pd
import pickle
import os
import tempfile
import traceback

# アプリのタイトルとヘッダー
st.title("CellChat QS Converter")
st.markdown("CellChat QSファイルをPython用pickleファイルに変換するツール")

# R環境の初期化 (セッション開始時に一度だけ実行)
@st.cache_resource
def initialize_r_environment():
    pandas2ri.activate()
    numpy2ri.activate()
    try:
        base = importr('base')
        qs = importr('qs')
        # CellChatもインポート
        ro.r('library(CellChat)')
        return True
    except Exception as e:
        st.error(f"R環境の初期化に失敗しました: {str(e)}")
        st.error("必要なRパッケージ(CellChat, qs)がインストールされているか確認してください。")
        return False

# R環境を初期化
r_initialized = initialize_r_environment()

# QSファイルを処理する関数
def process_qs_file(uploaded_file):
    # 一時ファイルに保存
    with tempfile.NamedTemporaryFile(suffix='.qs', delete=False) as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        tmp_path = tmp_file.name
    
    try:
        # 特殊文字をエスケープ
        tmp_path_escaped = tmp_path.replace('\\', '\\\\')
        
        # 最も単純なアプローチ: まずRで変数を設定し、それを使用する
        ro.r(f'temp_file_path <- "{tmp_path_escaped}"')
        
        r_code = """
        library(CellChat)
        library(qs)

        # QSファイルを読み込む
        tryCatch({
            cellchat_obj <- qs::qread(temp_file_path)
            
            # 必要なデータを抽出
            extracted_data <- list()

            # 1. メタデータ
            if (!is.null(cellchat_obj@meta)) {
                extracted_data$meta <- as.data.frame(cellchat_obj@meta)
            }

            # 2. ネットワーク情報
            if (!is.null(cellchat_obj@net$prob)) {
                # 確率行列
                extracted_data$net_prob <- cellchat_obj@net$prob
                
                # 次元名を個別に保存
                if (!is.null(dimnames(cellchat_obj@net$prob)[[1]])) {
                    extracted_data$dim_rows <- dimnames(cellchat_obj@net$prob)[[1]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[2]])) {
                    extracted_data$dim_cols <- dimnames(cellchat_obj@net$prob)[[2]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[3]])) {
                    extracted_data$dim_paths <- dimnames(cellchat_obj@net$prob)[[3]]
                }
            }

            if (!is.null(cellchat_obj@net$pval)) {
                extracted_data$net_pval <- cellchat_obj@net$pval
            }

            if (!is.null(cellchat_obj@net$weight)) {
                extracted_data$net_weight <- as.matrix(cellchat_obj@net$weight)
            }
            
            if (!is.null(cellchat_obj@net$count)) {
                extracted_data$net_count <- as.matrix(cellchat_obj@net$count)
            }

            # 3. パスウェイ情報
            if (!is.null(cellchat_obj@netP$pathways)) {
                extracted_data$pathways <- cellchat_obj@netP$pathways
            }

            if (!is.null(cellchat_obj@netP$prob)) {
                extracted_data$netP_prob <- cellchat_obj@netP$prob
            }

            # 4. LR情報
            if (!is.null(cellchat_obj@LR$LRsig)) {
                extracted_data$lr_sig <- as.data.frame(cellchat_obj@LR$LRsig)
            }

            # 返り値
            extracted_data
        }, error = function(e) {
            list(error = paste("QSファイルの読み込みエラー:", e$message))
        })
        """
        
        # Rコードを実行
        result_r = ro.r(r_code)
        
        # エラーチェック
        if isinstance(result_r, ro.vectors.ListVector) and 'error' in result_r.names:
            error_msg = result_r.rx2('error')
            st.error(f"QSファイルの処理に失敗しました: {error_msg}")
            return None

        # 結果の変換
        cellchat_data = {}
        for key in result_r.names:
            try:
                # キーの型によって変換方法を変える
                if key in ["net_prob", "net_pval", "netP_prob"]:
                    # R配列からNumPy配列への変換
                    cellchat_data[key] = np.array(result_r.rx2(key))
                else:
                    # その他のオブジェクトはpandas2riで変換
                    r_value = result_r.rx2(key)
                    try:
                        cellchat_data[key] = pandas2ri.rpy2py(r_value)
                    except Exception as e2:
                        # 変換に失敗した場合、基本的な方法でデータを取得
                        if isinstance(r_value, ro.vectors.StrVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.FloatVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.IntVector):
                            cellchat_data[key] = list(r_value)
                        else:
                            cellchat_data[key] = str(r_value)
                st.success(f"{key}を正常に抽出しました")
            except Exception as e:
                st.warning(f"{key}の変換に失敗しました: {e}")
                st.warning(traceback.format_exc())

        # 次元名の整理
        if all(key in cellchat_data for key in ['dim_rows', 'dim_cols', 'dim_paths']):
            cellchat_data['net_prob_dims'] = {
                'rows': cellchat_data['dim_rows'],
                'cols': cellchat_data['dim_cols'],
                'pathways': cellchat_data['dim_paths']
            }
            # 個別の次元キーを削除（オプション）
            for key in ['dim_rows', 'dim_cols', 'dim_paths']:
                if key in cellchat_data:
                    del cellchat_data[key]

        # 結果を整形
        result = format_cellchat_for_visualization(cellchat_data)
        
        # 一時ファイルの削除
        os.unlink(tmp_path)
        
        return result
        
    except Exception as e:
        st.error(f"処理中にエラーが発生しました: {e}")
        st.error(traceback.format_exc())
        # 一時ファイルの削除を試みる
        try:
            os.unlink(tmp_path)
        except:
            pass
        return None

def format_cellchat_for_visualization(cellchat_data):
    """抽出したデータをcellchat_vis.pyのビジュアライゼーション用に整形する"""
    # cellchat.pyの出力形式に合わせてデータを構造化
    result = {}
    
    # 1. ネットワーク情報
    # 3次元配列のprobとpvalがあるか確認
    has_3d_prob = 'net_prob' in cellchat_data and isinstance(cellchat_data['net_prob'], np.ndarray) and cellchat_data['net_prob'].ndim == 3
    has_3d_pval = 'net_pval' in cellchat_data and isinstance(cellchat_data['net_pval'], np.ndarray) and cellchat_data['net_pval'].ndim == 3
    
    result['net'] = {}
    
    # probとpvalが3次元の場合のみ設定
    if has_3d_prob:
        result['net']["prob"] = cellchat_data['net_prob']
    else:
        # データがなければ空の3次元配列を作成
        if 'net_prob_dims' in cellchat_data and all(key in cellchat_data['net_prob_dims'] for key in ['rows', 'cols']):
            rows = len(cellchat_data['net_prob_dims']['rows'])
            cols = len(cellchat_data['net_prob_dims']['cols'])
            paths = len(cellchat_data['net_prob_dims'].get('pathways', []))
            # 少なくとも1つのパスウェイを保証
            paths = max(1, paths)
            result['net']["prob"] = np.zeros((rows, cols, paths))
        else:
            # 次元情報もない場合は最小サイズの配列
            result['net']["prob"] = np.zeros((1, 1, 1))
    
    if has_3d_pval:
        result['net']["pval"] = cellchat_data['net_pval']
    else:
        # probと同じサイズの配列を作成
        result['net']["pval"] = np.ones_like(result['net']["prob"])
    
    # 次元名の情報を追加
    if 'net_prob_dims' in cellchat_data:
        dims = cellchat_data['net_prob_dims']
        cell_types_rows = dims.get('rows', [])
        cell_types_cols = dims.get('cols', [])
        pathways = dims.get('pathways', [])
        
        # 文字列リストに変換
        cell_types_rows = [str(ct) for ct in cell_types_rows]
        cell_types_cols = [str(ct) for ct in cell_types_cols]
        pathways = [str(p) for p in pathways]
        
        # 最小限のサイズを保証
        if not cell_types_rows:
            cell_types_rows = [f"Cell_{i}" for i in range(result['net']["prob"].shape[0])]
        if not cell_types_cols:
            cell_types_cols = [f"Cell_{i}" for i in range(result['net']["prob"].shape[1])]
        if not pathways:
            pathways = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]
        
        result['net']['dimnames'] = [
            cell_types_rows,
            cell_types_cols,
            pathways
        ]
    else:
        # 次元名がない場合は配列のサイズに合わせて生成
        shape = result['net']["prob"].shape
        result['net']['dimnames'] = [
            [f"Cell_{i}" for i in range(shape[0])],
            [f"Cell_{i}" for i in range(shape[1])],
            [f"Pathway_{i}" for i in range(shape[2])]
        ]
    
    # 集計したweight行列を作成
    cell_types = result['net']['dimnames'][0]
    num_cells = len(cell_types)
    
    # 全パスウェイにわたる合計（集計ネットワーク）
    prob_sum = np.sum(result['net']["prob"], axis=2)
    
    # weight行列とcount行列をDataFrameに変換
    result['net']["weight"] = prob_sum
    result['net']["count"] = (prob_sum > 0).astype(int)
    
    # 集計ネットワークの中心性を計算（NetworkXを使用）
    import networkx as nx
    
    # 初期化
    net_centr = {
        'outdeg': np.zeros(num_cells),
        'indeg': np.zeros(num_cells),
        'outdeg_unweighted': np.zeros(num_cells),
        'indeg_unweighted': np.zeros(num_cells),
        'betweenness': np.zeros(num_cells),
        'page_rank': np.zeros(num_cells),
        'hub': np.zeros(num_cells),
        'authority': np.zeros(num_cells),
        'eigen': np.zeros(num_cells),
        'flowbet': np.zeros(num_cells),
        'info': np.zeros(num_cells)
    }
    
    # ネットワークが有効（エッジがある）場合のみ計算
    if np.any(prob_sum > 0):
        # グラフ作成
        G = nx.DiGraph()
        for i in range(num_cells):
            G.add_node(i, name=cell_types[i])
        
        # エッジ追加
        for i in range(num_cells):
            for j in range(num_cells):
                if prob_sum[i, j] > 0:
                    G.add_edge(i, j, weight=float(prob_sum[i, j]))
        
        # 出入次数（加重/非加重）
        for i in range(num_cells):
            # 加重次数
            net_centr['outdeg'][i] = sum(G[i][j]['weight'] for j in G.successors(i)) if i in G else 0
            net_centr['indeg'][i] = sum(G[j][i]['weight'] for j in G.predecessors(i)) if i in G else 0
            
            # 非加重次数
            net_centr['outdeg_unweighted'][i] = G.out_degree(i) if i in G else 0
            net_centr['indeg_unweighted'][i] = G.in_degree(i) if i in G else 0
        
        # 重み逆数グラフ（媒介中心性用）
        G_inv = G.copy()
        for u, v, d in G_inv.edges(data=True):
            if d['weight'] > 0:
                d['weight'] = 1.0 / d['weight']
        
        try:
            # 媒介中心性
            betweenness = nx.betweenness_centrality(G_inv, weight='weight')
            for i in range(num_cells):
                net_centr['betweenness'][i] = betweenness.get(i, 0)
                # flowbetはbetweennessで代用
                net_centr['flowbet'][i] = betweenness.get(i, 0)
        except:
            pass
        
        try:
            # PageRank
            page_rank = nx.pagerank(G, weight='weight')
            for i in range(num_cells):
                net_centr['page_rank'][i] = page_rank.get(i, 0)
        except:
            pass
        
        try:
            # HubとAuthority
            hub, authority = nx.hits(G, max_iter=1000)
            for i in range(num_cells):
                net_centr['hub'][i] = hub.get(i, 0)
                net_centr['authority'][i] = authority.get(i, 0)
        except:
            pass
        
        try:
            # 固有ベクトル中心性
            eigen = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
            for i in range(num_cells):
                net_centr['eigen'][i] = eigen.get(i, 0)
        except:
            pass
        
        # Information Centrality - カスタム計算（近似）
        try:
            # info中心性の簡易近似 (正規化したdegreeとbetweennessの組み合わせ)
            in_degree = np.array([net_centr['indeg'][i] for i in range(num_cells)])
            out_degree = np.array([net_centr['outdeg'][i] for i in range(num_cells)])
            between = np.array([net_centr['betweenness'][i] for i in range(num_cells)])
            
            # 正規化
            norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
            norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
            norm_between = between / (np.max(between) if np.max(between) > 0 else 1)
            
            # 情報中心性を合成
            info_score = (norm_in + norm_out + norm_between) / 3
            for i in range(num_cells):
                net_centr['info'][i] = info_score[i]
        except:
            pass
    
    # 中心性指標をnetに追加
    result['net']['centr'] = net_centr
    
    # 2. パスウェイ情報
    result['netP'] = {}
    
    # パスウェイ情報の追加
    if 'pathways' in cellchat_data and cellchat_data['pathways'] is not None:
        result['netP']['pathways'] = cellchat_data['pathways']
    elif 'net_prob_dims' in cellchat_data and 'pathways' in cellchat_data['net_prob_dims']:
        result['netP']['pathways'] = cellchat_data['net_prob_dims']['pathways']
    else:
        # それでも見つからない場合はnet['dimnames']から取得
        if 'dimnames' in result['net'] and len(result['net']['dimnames']) > 2:
            result['netP']['pathways'] = result['net']['dimnames'][2]
        else:
            result['netP']['pathways'] = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]
    
    # ネットワーク確率情報の追加
    if 'netP_prob' in cellchat_data and cellchat_data['netP_prob'] is not None:
        result['netP']['prob'] = cellchat_data['netP_prob']
    else:
        # netP_probがない場合、パスウェイごとの集計を行う
        prob = result['net']["prob"]
        num_pathways = len(result['netP']['pathways'])
        
        # 各パスウェイの総和を計算
        netP_prob = np.zeros((num_cells, num_cells, num_pathways))
        for k in range(min(prob.shape[2], num_pathways)):
            netP_prob[:, :, k] = prob[:, :, k]
        
        result['netP']['prob'] = netP_prob
    
    # netP用の中心性情報を計算
    netP_centr = {}
    
    # 各パスウェイについて中心性を計算
    for k in range(min(result['netP']['prob'].shape[2], len(result['netP']['pathways']))):
        # このパスウェイの初期中心性辞書
        netP_centr[k] = {
            'outdeg': np.zeros(num_cells),
            'indeg': np.zeros(num_cells),
            'outdeg_unweighted': np.zeros(num_cells),
            'indeg_unweighted': np.zeros(num_cells),
            'betweenness': np.zeros(num_cells),
            'page_rank': np.zeros(num_cells),
            'hub': np.zeros(num_cells),
            'authority': np.zeros(num_cells),
            'eigen': np.zeros(num_cells),
            'flowbet': np.zeros(num_cells),
            'info': np.zeros(num_cells)
        }
        
        # パスウェイ行列
        pathway_weight = result['netP']['prob'][:, :, k]
        
        # ネットワークが有効（エッジがある）場合のみ計算
        if np.any(pathway_weight > 0):
            # グラフ作成
            G_pathway = nx.DiGraph()
            for i in range(num_cells):
                G_pathway.add_node(i, name=cell_types[i])
            
            # エッジ追加
            for i in range(num_cells):
                for j in range(num_cells):
                    if pathway_weight[i, j] > 0:
                        G_pathway.add_edge(i, j, weight=float(pathway_weight[i, j]))
            
            # 出入次数（加重/非加重）
            for i in range(num_cells):
                # 加重次数
                netP_centr[k]['outdeg'][i] = sum(G_pathway[i][j]['weight'] for j in G_pathway.successors(i)) if i in G_pathway else 0
                netP_centr[k]['indeg'][i] = sum(G_pathway[j][i]['weight'] for j in G_pathway.predecessors(i)) if i in G_pathway else 0
                
                # 非加重次数
                netP_centr[k]['outdeg_unweighted'][i] = G_pathway.out_degree(i) if i in G_pathway else 0
                netP_centr[k]['indeg_unweighted'][i] = G_pathway.in_degree(i) if i in G_pathway else 0
            
            # 重み逆数グラフ（媒介中心性用）
            G_inv = G_pathway.copy()
            for u, v, d in G_inv.edges(data=True):
                if d['weight'] > 0:
                    d['weight'] = 1.0 / d['weight']
            
            try:
                # 媒介中心性
                betweenness = nx.betweenness_centrality(G_inv, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['betweenness'][i] = betweenness.get(i, 0)
                    # flowbetはbetweennessで代用
                    netP_centr[k]['flowbet'][i] = betweenness.get(i, 0)
            except:
                pass
            
            try:
                # PageRank
                page_rank = nx.pagerank(G_pathway, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['page_rank'][i] = page_rank.get(i, 0)
            except:
                pass
            
            try:
                # HubとAuthority
                hub, authority = nx.hits(G_pathway, max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['hub'][i] = hub.get(i, 0)
                    netP_centr[k]['authority'][i] = authority.get(i, 0)
            except:
                pass
            
            try:
                # 固有ベクトル中心性
                eigen = nx.eigenvector_centrality(G_pathway, weight='weight', max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['eigen'][i] = eigen.get(i, 0)
            except:
                pass
            
            # Information Centrality - カスタム計算
            try:
                # info中心性の簡易近似 (正規化したdegreeとbetweennessの組み合わせ)
                in_degree = np.array([netP_centr[k]['indeg'][i] for i in range(num_cells)])
                out_degree = np.array([netP_centr[k]['outdeg'][i] for i in range(num_cells)])
                between = np.array([netP_centr[k]['betweenness'][i] for i in range(num_cells)])
                
                # 正規化
                norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
                norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
                norm_between = between / (np.max(between) if np.max(between) > 0 else 1)
                
                # 情報中心性を合成
                info_score = (norm_in + norm_out + norm_between) / 3
                for i in range(num_cells):
                    netP_centr[k]['info'][i] = info_score[i]
            except:
                pass
    
    # 中心性情報を追加
    result['netP']['centr'] = netP_centr
    
    # 3. 結果DataFrame作成
    # 細胞タイプ名とパスウェイ名を取得
    cell_types_rows = result['net']['dimnames'][0]
    cell_types_cols = result['net']['dimnames'][1]
    pathways = result['net']['dimnames'][2]
    
    # prob, pval配列から相互作用データを抽出
    prob = result['net']["prob"]
    pval = result['net']["pval"]
    
    results_data = []
    for i in range(prob.shape[0]):
        for j in range(prob.shape[1]):
            for k in range(min(prob.shape[2], len(pathways))):
                if prob[i, j, k] > 0:
                    # 相互作用名
                    interaction_name = pathways[k] if k < len(pathways) else f"interaction_{k}"
                    
                    # リガンド-レセプター情報（あれば）
                    ligand = ""
                    receptor = ""
                    
                    # LR情報取得
                    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
                        lr_df = cellchat_data['lr_sig']
                        if interaction_name in lr_df.index:
                            if 'ligand' in lr_df.columns:
                                ligand = lr_df.loc[interaction_name, 'ligand']
                            if 'receptor' in lr_df.columns:
                                receptor = lr_df.loc[interaction_name, 'receptor']
                    
                    results_data.append({
                        'source': cell_types_rows[i] if i < len(cell_types_rows) else f"Cell_{i}",
                        'target': cell_types_cols[j] if j < len(cell_types_cols) else f"Cell_{j}",
                        'interaction_name': interaction_name,
                        'ligand': str(ligand),
                        'receptor': str(receptor),
                        'prob': float(prob[i, j, k]),
                        'pval': float(pval[i, j, k])
                    })
    
    result['results'] = pd.DataFrame(results_data)
    
    # 4. ネットワークサマリー作成
    # 細胞タイプごとの集計（probの合計）
    prob = result['net']["prob"]
    strength_matrix = np.sum(prob, axis=2)
    
    # DataFrameに変換
    strength_df = pd.DataFrame(strength_matrix, 
                              index=cell_types, 
                              columns=cell_types)
    
    # カウント行列（非ゼロの相互作用数）
    count_matrix = np.sum(prob > 0, axis=2)
    count_df = pd.DataFrame(count_matrix,
                           index=cell_types,
                           columns=cell_types)
    
    # 送受信総量
    outgoing = pd.Series(np.sum(strength_matrix, axis=1), index=cell_types)
    incoming = pd.Series(np.sum(strength_matrix, axis=0), index=cell_types)
    
    # LR寄与度の計算
    lr_contribution = np.zeros(prob.shape[2])
    for k in range(prob.shape[2]):
        lr_contribution[k] = np.sum(prob[:, :, k])
    
    # ネットワーク中心性データフレームを作成
    network_centrality_dict = {
        'cell_type': cell_types,
    }
    
    # 各中心性指標をディクショナリに追加
    for metric in ['outdeg', 'indeg', 'outdeg_unweighted', 'indeg_unweighted', 
                  'betweenness', 'page_rank', 'hub', 'authority', 'eigen', 'flowbet', 'info']:
        network_centrality_dict[metric] = net_centr[metric]
    
    network_centrality = pd.DataFrame(network_centrality_dict)
    
    # コア中心性指標をより意味のある名前に変更（シグナリング役割分析用）
    role_data = {
        'cell_type': cell_types,
        'sender': net_centr['outdeg'],        # 送信者としての役割
        'receiver': net_centr['indeg'],       # 受信者としての役割
        'mediator': net_centr['flowbet'],     # 仲介者としての役割
        'influencer': net_centr['info']       # 影響力としての役割
    }
    
    # 役割データフレームを作成
    signaling_role = pd.DataFrame(role_data)
    
    result['network'] = {
        'strength_matrix': strength_df,
        'count_matrix': count_df,
        'outgoing': outgoing,
        'incoming': incoming,
        'lr_contribution': lr_contribution,
        'network_centrality': network_centrality,
        'signaling_role': signaling_role      # シグナリング役割分析用
    }
    
    # 5. メタデータ情報
    if 'meta' in cellchat_data:
        result['meta'] = cellchat_data['meta']
    
    # 6. GroupBy情報（デフォルト設定）
    result['groupby'] = 'cell.ident'  # デフォルト値
    
    # 7. LR情報
    if 'lr_sig' in cellchat_data:
        result['lr_sig'] = cellchat_data['lr_sig']
    
    # 8. adata情報 - 疑似的なanndata構造を作成
    from types import SimpleNamespace
    
    # adataオブジェクトを作成
    adata = SimpleNamespace()
    
    # メタデータを追加
    if 'meta' in cellchat_data and isinstance(cellchat_data['meta'], pd.DataFrame):
        adata.obs = cellchat_data['meta']
    else:
        adata.obs = pd.DataFrame(index=cell_types)
    
    # 遺伝子名リストを抽出
    gene_names = []
    
    # LRsigからリガンドとレセプター情報を抽出
    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
        lr_df = cellchat_data['lr_sig']
        if 'ligand' in lr_df.columns:
            ligands = lr_df['ligand'].dropna().astype(str).unique()
            gene_names.extend(ligands)
        if 'receptor' in lr_df.columns:
            receptors = lr_df['receptor'].dropna().astype(str).unique()
            gene_names.extend(receptors)
    
    # 遺伝子名がない場合のフォールバック
    if not gene_names:
        gene_names = ["Gene1", "Gene2", "Gene3"]
    
    # 重複を削除
    gene_names = sorted(list(set(gene_names)))
    
    # AnnDataのvar属性とvar_namesを模倣
    adata.var_names = pd.Index(gene_names)
    adata.var = pd.DataFrame(index=adata.var_names)
    
    # AnnDataのshape属性を追加
    adata.shape = (len(adata.obs), len(adata.var_names))
    
    # X属性を追加
    adata.X = np.zeros(adata.shape)
    
    result['adata'] = adata
    
    return result

# メイン処理
if r_initialized:
    # ファイルアップローダー
    uploaded_file = st.file_uploader("CellChat QSファイルをアップロード", type=["qs"])
    
    if uploaded_file is not None:
        st.info(f"ファイル「{uploaded_file.name}」を処理しています...")
        
        with st.spinner("QSファイルを処理中..."):
            result = process_qs_file(uploaded_file)
        
        if result is not None:
            st.success("QSファイルの変換が完了しました！")
            
            # 結果の概要を表示
            st.subheader("変換結果の概要")
            
            # 細胞タイプ情報の表示
            if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 0:
                cell_types = result['net']['dimnames'][0]
                cell_count = len(cell_types)
                if cell_count > 0:
                    st.write(f"検出された細胞タイプ数: {cell_count}")
                    if cell_count <= 20:
                        st.write(f"細胞タイプ: {', '.join(cell_types)}")
                    else:
                        with st.expander("細胞タイプ一覧"):
                            st.write(", ".join(cell_types))
            
            # パスウェイ情報の表示
            if 'netP' in result and 'pathways' in result['netP'] and result['netP']['pathways'] is not None:
                pathways = result['netP']['pathways']
                pathway_count = len(pathways)
                if pathway_count > 0:
                    st.write(f"シグナル経路数: {pathway_count}")
                    if pathway_count <= 10:
                        st.write(f"シグナル経路: {', '.join(pathways)}")
                    else:
                        with st.expander("シグナル経路一覧"):
                            st.write(", ".join(pathways))
            
            # 相互作用情報の表示
            if 'results' in result and not result['results'].empty:
                st.write(f"抽出された相互作用数: {len(result['results'])}")
                with st.expander("相互作用データのサンプル（上位5件）"):
                    st.dataframe(result['results'].head())
            
            # 結果をPickleファイルに保存
            output_filename = os.path.splitext(uploaded_file.name)[0] + "_python.pkl"
            pickle_data = pickle.dumps(result)
            
            # ダウンロードボタン
            st.download_button(
                label="変換結果をダウンロード",
                data=pickle_data,
                file_name=output_filename,
                mime="application/octet-stream"
            )
            
            st.info("""
            ダウンロードしたpickleファイルは、cellchat.pyでCellChatの可視化に使用できます。
            Streamlitアプリにアップロードするか、コード内でpickle.load()を使用してロードしてください。
            """)
else:
    st.error("R環境の初期化に失敗したため、アプリを実行できません。必要なパッケージがインストールされているか確認してください。")