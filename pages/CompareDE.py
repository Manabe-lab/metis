import streamlit as st
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import numpy as np
from upsetplot import plot
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import io
from difflib import SequenceMatcher
import itertools

def process_filename(filenames):
    """
    ファイル名を2つずつ比較し、共通する部分を除去する
    """
    # 拡張子を除去
    names = [name.rsplit('.', 1)[0] for name in filenames]
    
    def compare_two_names(name1, name2):
        """2つの名前を比較し、それぞれの固有の部分を見つける"""
        # 名前を分割し、空の要素を除去
        parts1 = [p for p in name1.replace('-', '_').split('_') if p]
        parts2 = [p for p in name2.replace('-', '_').split('_') if p]
        
        # 共通部分を見つける
        common = set(parts1) & set(parts2)
        
        # 各名前の固有部分を抽出
        unique1 = [p for p in parts1 if p not in common]
        unique2 = [p for p in parts2 if p not in common]
        
        return unique1, unique2

    # 各ファイルの固有部分を保持
    processed = {f: set() for f in filenames}
    
    # 2つずつ比較
    for i in range(len(filenames)):
        for j in range(i + 1, len(filenames)):
            unique1, unique2 = compare_two_names(names[i], names[j])
            # 固有部分を蓄積
            processed[filenames[i]].update(unique1)
            processed[filenames[j]].update(unique2)
    
    # 結果をクリーンアップ
    result = {}
    for f in filenames:
        unique_parts = list(processed[f])
        if unique_parts:
            # 空の要素を除去して結合
            result[f] = '_'.join(filter(None, unique_parts)).strip('_')
        else:
            # 元の名前から空の要素を除去
            parts = [p for p in names[filenames.index(f)].replace('-', '_').split('_') if p]
            result[f] = '_'.join(parts)
    
    return result
    

@st.cache_data
def read_excel(file, index_col=None, header=0):
    return pd.read_excel(file, index_col=index_col, header=header)

@st.cache_data
def read_csv(file, index_col=None, sep=',', header=0):
    return pd.read_csv(file, index_col=index_col, header=header, sep=sep, engine='python')

def get_regulated_genes(df, lfc_col, pval_col, direction='up', 
                       lfc_threshold=1, pval_threshold=0.05):
    """指定した方向に発現変動する遺伝子を抽出する関数"""
    if direction == 'up':
        return set(df.index[
            (df[lfc_col] > lfc_threshold) & 
            (df[pval_col] < pval_threshold)
        ])
    else:  # down
        return set(df.index[
            (df[lfc_col] < -lfc_threshold) & 
            (df[pval_col] < pval_threshold)
        ])

def create_venn_diagram(datasets, names, title):
    """ベン図を作成する関数"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    if len(datasets) == 2:
        venn2(
            subsets=[datasets[0], datasets[1]],
            set_labels=names,
            ax=ax
        )
    elif len(datasets) == 3:
        venn3(
            subsets=[datasets[0], datasets[1], datasets[2]],
            set_labels=names,
            ax=ax
        )
    
    plt.title(title)
    return fig

def create_upset_plot(datasets, names, title):
    """
    UpSetプロット風の可視化を作成する関数
    """
    # 共通部分の組み合わせを見つける
    def get_combinations():
        n = len(datasets)
        combos = []
        sizes = []
        labels = []
        
        # 各組み合わせについて
        for i in range(1, n + 1):
            for combo in itertools.combinations(range(n), i):
                # その組み合わせに含まれるデータセットの共通部分
                common = set.intersection(*[datasets[j] for j in combo])
                if len(common) > 0:  # 共通部分が存在する場合のみ追加
                    combos.append(combo)
                    sizes.append(len(common))
                    labels.append(' & '.join([names[j] for j in combo]))
        
        return combos, sizes, labels

    combos, sizes, labels = get_combinations()
    
    # サイズでソート
    sorted_indices = np.argsort(sizes)[::-1]
    sizes = [sizes[i] for i in sorted_indices]
    labels = [labels[i] for i in sorted_indices]
    combos = [combos[i] for i in sorted_indices]

    # プロットの作成
    fig, (ax_sets, ax_matrix) = plt.subplots(2, 1, 
                                            figsize=(12, 8),
                                            gridspec_kw={'height_ratios': [3, 1]})
    
    # バーチャートの描画
    bars = ax_sets.bar(range(len(sizes)), sizes)
    ax_sets.set_xticks([])
    ax_sets.set_ylabel('Intersection Size')
    
    # 行列の描画
    matrix = np.zeros((len(names), len(combos)))
    for i, combo in enumerate(combos):
        for j in combo:
            matrix[j, i] = 1
    
    ax_matrix.imshow(matrix, cmap='binary', aspect='auto')
    ax_matrix.set_yticks(range(len(names)))
    ax_matrix.set_yticklabels(names)
    ax_matrix.set_xticks(range(len(labels)))
    ax_matrix.set_xticklabels(labels, rotation=45, ha='right')
    
    # レイアウトの調整
    plt.tight_layout()
    ax_sets.set_title(title)
    
    return fig


def create_2d_scatter(df1, df2, x_col, y_col, name1, name2, pval_col1, pval_col2,
                     highlight_significant=True, lfc_threshold=1, pval_threshold=0.05,
                     use_union=False):
    """2次元散布図を作成"""
    # インデックスの処理（共通 or 和集合）
    if use_union:
        # 和集合を使用し、欠損値を0で埋める
        all_index = df1.index.union(df2.index)
        df1 = df1.reindex(all_index).fillna(0)
        df2 = df2.reindex(all_index).fillna(0)
    else:
        # 共通部分のみ使用
        common_index = df1.index.intersection(df2.index)
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]

    fig = go.Figure()
    
    # ホバーテキストの作成
    hovertemplate = (
        "遺伝子: %{customdata}<br>" +
        f"{name1} {x_col}: %{{x:.3f}}<br>" +
        f"{name2} {y_col}: %{{y:.3f}}<br>" +
        "<extra></extra>"
    )

    if 'neg_log_padj' not in x_col:  # Log2FC の散布図の場合
        if highlight_significant:
            # データセット1で有意
            sig1 = (abs(df1[x_col]) > lfc_threshold) & (df1[pval_col1] < pval_threshold)
            # データセット2で有意
            sig2 = (abs(df2[y_col]) > lfc_threshold) & (df2[pval_col2] < pval_threshold)
            
            # 両方で有意
            both_sig = sig1 & sig2
            # どちらかのみで有意
            only_sig1 = sig1 & ~sig2
            only_sig2 = ~sig1 & sig2
            # 非有意
            not_sig = ~sig1 & ~sig2

            # 各カテゴリーのプロット
            fig.add_trace(go.Scatter(
                x=df1[x_col][both_sig],
                y=df2[y_col][both_sig],
                mode='markers',
                name='Both significant',
                customdata=df1.index[both_sig],
                hovertemplate=hovertemplate,
                marker=dict(color='red', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][only_sig1],
                y=df2[y_col][only_sig1],
                mode='markers',
                name=f'Only in {name1}',
                customdata=df1.index[only_sig1],
                hovertemplate=hovertemplate,
                marker=dict(color='blue', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][only_sig2],
                y=df2[y_col][only_sig2],
                mode='markers',
                name=f'Only in {name2}',
                customdata=df1.index[only_sig2],
                hovertemplate=hovertemplate,
                marker=dict(color='green', size=4)
            ))
            fig.add_trace(go.Scatter(
                x=df1[x_col][not_sig],
                y=df2[y_col][not_sig],
                mode='markers',
                name='Not significant',
                customdata=df1.index[not_sig],
                hovertemplate=hovertemplate,
                marker=dict(color='grey', size=2, opacity=0.5)
            ))
        else:
            fig.add_trace(go.Scatter(
                x=df1[x_col],
                y=df2[y_col],
                mode='markers',
                name='Genes',
                customdata=df1.index,
                hovertemplate=hovertemplate,
                marker=dict(
                    color='blue',
                    size=3,
                    opacity=0.6
                )
            ))
    else:  # P値の散布図の場合
        # P値による色分け
        low_p1 = df1[pval_col1] < pval_threshold
        low_p2 = df2[pval_col2] < pval_threshold
        
        # 両方で有意
        both_sig = low_p1 & low_p2
        # どちらかのみで有意
        only_sig1 = low_p1 & ~low_p2
        only_sig2 = ~low_p1 & low_p2
        # 非有意
        not_sig = ~low_p1 & ~low_p2
        
        # 各カテゴリーのプロット
        fig.add_trace(go.Scatter(
            x=df1[x_col][both_sig],
            y=df2[y_col][both_sig],
            mode='markers',
            name='Both significant',
            customdata=df1.index[both_sig],
            hovertemplate=hovertemplate,
            marker=dict(color='red', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][only_sig1],
            y=df2[y_col][only_sig1],
            mode='markers',
            name=f'Only in {name1}',
            customdata=df1.index[only_sig1],
            hovertemplate=hovertemplate,
            marker=dict(color='blue', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][only_sig2],
            y=df2[y_col][only_sig2],
            mode='markers',
            name=f'Only in {name2}',
            customdata=df1.index[only_sig2],
            hovertemplate=hovertemplate,
            marker=dict(color='green', size=4)
        ))
        fig.add_trace(go.Scatter(
            x=df1[x_col][not_sig],
            y=df2[y_col][not_sig],
            mode='markers',
            name='Not significant',
            customdata=df1.index[not_sig],
            hovertemplate=hovertemplate,
            marker=dict(color='grey', size=2, opacity=0.5)
        ))

    # 補助線の追加
    fig.add_hline(y=0, line_dash="dash", line_color="gray")
    fig.add_vline(x=0, line_dash="dash", line_color="gray")

    fig.update_layout(
        title=f"Comparison of {name1} vs {name2}",
        xaxis_title=f"{name1} - {x_col}",
        yaxis_title=f"{name2} - {y_col}",
        hovermode='closest',
        showlegend=True
    )

    return fig


def create_3d_scatter(df1, df2, df3, x_col, y_col, z_col, name1, name2, name3,
                     pval_col1, pval_col2, pval_col3,
                     highlight_significant=True, lfc_threshold=1, pval_threshold=0.05,
                     use_union=False):
    """3次元散布図を作成"""
    # インデックスの処理（共通 or 和集合）
    if use_union:
        # 和集合を使用し、欠損値を0で埋める
        all_index = df1.index.union(df2.index).union(df3.index)
        df1 = df1.reindex(all_index).fillna(0)
        df2 = df2.reindex(all_index).fillna(0)
        df3 = df3.reindex(all_index).fillna(0)
    else:
        # 共通部分のみ使用
        common_index = df1.index.intersection(df2.index).intersection(df3.index)
        df1 = df1.loc[common_index]
        df2 = df2.loc[common_index]
        df3 = df3.loc[common_index]

    fig = go.Figure()
    
    # ホバーテキストの作成
    hovertemplate = (
        "遺伝子: %{customdata}<br>" +
        f"{name1} {x_col}: %{{x:.3f}}<br>" +
        f"{name2} {y_col}: %{{y:.3f}}<br>" +
        f"{name3} {z_col}: %{{z:.3f}}<br>" +
        "<extra></extra>"
    )

    if highlight_significant:
        # 各データセットでの有意性
        sig1 = (abs(df1[x_col]) > lfc_threshold) & (df1[pval_col1] < pval_threshold)
        sig2 = (abs(df2[y_col]) > lfc_threshold) & (df2[pval_col2] < pval_threshold)
        sig3 = (abs(df3[z_col]) > lfc_threshold) & (df3[pval_col3] < pval_threshold)

        # 全てで有意
        all_sig = sig1 & sig2 & sig3
        # それ以外
        other = ~all_sig

        # プロット
        fig.add_trace(go.Scatter3d(
            x=df1[x_col][all_sig],
            y=df2[y_col][all_sig],
            z=df3[z_col][all_sig],
            mode='markers',
            name='Significant in all',
            customdata=df1.index[all_sig],
            hovertemplate=hovertemplate,
            marker=dict(size=4, color='red')
        ))
        fig.add_trace(go.Scatter3d(
            x=df1[x_col][other],
            y=df2[y_col][other],
            z=df3[z_col][other],
            mode='markers',
            name='Other',
            customdata=df1.index[other],
            hovertemplate=hovertemplate,
            marker=dict(size=2, color='grey', opacity=0.5)
        ))
    else:
        fig.add_trace(go.Scatter3d(
            x=df1[x_col],
            y=df2[y_col],
            z=df3[z_col],
            mode='markers',
            name='Genes',
            customdata=df1.index,
            hovertemplate=hovertemplate,
            marker=dict(size=3, color='blue', opacity=0.6)
        ))

    fig.update_layout(
        title="3D Comparison",
        scene=dict(
            xaxis_title=f"{name1} - {x_col}",
            yaxis_title=f"{name2} - {y_col}",
            zaxis_title=f"{name3} - {z_col}"
        ),
        margin=dict(l=0, r=0, b=0, t=30)
    )

    return fig


def main():
    st.title("DE results comparison")
    
    # ファイルアップロード
    uploaded_files = st.file_uploader(
        "DE解析結果ファイルを選択してください（複数可）", 
        accept_multiple_files=True,
        type=['csv', 'tsv', 'txt', 'xlsx', 'xls']
    )
    st.write("P-value, FCを持つファイルを２つ")
    
    if not uploaded_files:
        st.info("解析結果ファイルをアップロードしてください。")
        return
    
    # ファイル名の処理
    filename_mapping = process_filename([f.name for f in uploaded_files])
    
    # 閾値の設定
    col1, col2 = st.columns(2)
    with col1:
        lfc_threshold = st.number_input(
            "Log2 Fold Change threshold",
            value=1.0,
            step=0.1,
            min_value = 0.0
        )
    with col2:
        pval_threshold = st.number_input(
            "FDR threshold",
            value=0.05,
            step=0.01,
            min_value = 0.00,
            max_value = 1.00,
            format="%.3f"
        )
    
    # データセットの読み込みと処理
    datasets = {}
    for uploaded_file in uploaded_files:
        try:
            try:
                df = read_csv(uploaded_file, index_col=0, sep=None)
            except:
                df = read_excel(uploaded_file, index_col=0)
            numeric_cols = df.select_dtypes(include=[np.number]).columns
            datasets[filename_mapping[uploaded_file.name]] = {
                'data': df,
                'columns': numeric_cols
            }
        except Exception as e:
            st.error(f"エラー ({uploaded_file.name}): {e}")
    
    if len(datasets) < 2:
        st.warning("比較には少なくとも2つのデータセットが必要です。")
        return
    
    # データセット設定
    st.subheader("データセットとカラムの設定")
    
    dataset_settings = []
    for name, dataset in datasets.items():
        st.markdown(f"#### {name}")
        # カラム名のパターンマッチング
        pvalue = [i for i in dataset['columns'] if ('adj' in i.lower()) or  ('fdr' in i.lower()) or
         ('pvalue' in i.lower()) or ('p-val' in i.lower()) or  ('p val' in i.lower()) 
                 ]




        fc = [i for i in dataset['columns'] if ('log2fc' in i.lower()) or 
              ('fold change' in i.lower()) or ('log2foldchange' in i.lower()) or 
              ('logfc' in i.lower()) or ('coef' in i.lower())]
        
        lfc_col = st.selectbox(
            "Log2 Fold Change",
            fc,
            key=f"lfc_{name}"
        )

        pval_col = st.selectbox(
            "adjusted P",
            pvalue,
            key=f"pval_{name}"
        )
        use_in_analysis = True
        invert_sign = st.checkbox(
            "Log2 Fold Changeの符号を反転",
            key=f"invert_{name}",
            help="コントロールの定義が逆の場合にチェックしてください"
        )
        
        # データセットの設定を保存
        dataset_settings.append({
            'name': name,
            'data': dataset['data'].copy(),  # Create a copy to avoid modifying original
            'lfc_col': lfc_col,
            'pval_col': pval_col,
            'invert_sign': invert_sign,
            'use_in_analysis': use_in_analysis
        })
        
        # Invert sign if requested
        if invert_sign:
            dataset_settings[-1]['data'][lfc_col] = -dataset_settings[-1]['data'][lfc_col]

        st.markdown("---")
    
    # 解析タイプの選択
    analysis_type = st.radio(
        "解析タイプを選択：",
        ["発現変動遺伝子の比較", "散布図による比較"],
        horizontal=True
    )
    
    active_datasets = [ds for ds in dataset_settings if ds['use_in_analysis']]
    
    if analysis_type == "発現変動遺伝子の比較":
        # 発現方向性の選択
        direction = st.radio(
            "発現変動の方向性を選択：",
            ["Up-regulated", "Down-regulated", "Both"],
            horizontal=True
        )
        
        # プロットタイプの選択
        plot_type = st.radio(
            "プロットタイプを選択：",
            ["ベン図", "UpSetプロット"],
            horizontal=True
        )
        
        if direction in ["Up-regulated", "Down-regulated"]:
            dir_type = 'up' if direction == "Up-regulated" else 'down'
            significant_sets = []
            set_names = []
            
            for ds in active_datasets:
                sig_genes = get_regulated_genes(
                    ds['data'], 
                    ds['lfc_col'], 
                    ds['pval_col'],
                    dir_type,
                    lfc_threshold,
                    pval_threshold
                )
                significant_sets.append(sig_genes)
                set_names.append(ds['name'])
            
            title = f"{direction}"
            
            if plot_type == "ベン図" and len(significant_sets) <= 3:
                fig = create_venn_diagram(significant_sets, set_names, title)
                st.pyplot(fig)
            else:
                fig = create_upset_plot(significant_sets, set_names, title)
                st.pyplot(fig)
            
            # 共通遺伝子の表示
            common_genes = set.intersection(*significant_sets)
            st.write(f"共通して{direction}する遺伝子数: {len(common_genes)}")
            
            if st.checkbox(f"共通{direction}遺伝子リストを表示"):
                st.write(", ".join(sorted(list(common_genes))))
            
            # 共通遺伝子のエクスポート
            if len(common_genes) > 0:
                csv = pd.DataFrame(sorted(list(common_genes)), columns=['Gene_ID'])
                csv_data = csv.to_csv(index=False)
                st.download_button(
                    label=f"共通{direction}遺伝子リストをダウンロード",
                    data=csv_data,
                    file_name=f"common_{dir_type}_regulated_genes.csv",
                    mime="text/csv"
                )
        
        else:  # 両方別々に表示
            for dir_type, dir_name in [('up', 'Up-regulated'), ('down', 'Down-regulated')]:
                st.subheader(f"{dir_name}遺伝子の比較")
                significant_sets = []
                set_names = []
                
                for ds in active_datasets:
                    sig_genes = get_regulated_genes(
                        ds['data'], 
                        ds['lfc_col'], 
                        ds['pval_col'],
                        dir_type,
                        lfc_threshold,
                        pval_threshold
                    )
                    significant_sets.append(sig_genes)
                    set_names.append(ds['name'])
                
                title = f"{dir_name}遺伝子の比較"
                
                if plot_type == "ベン図" and len(significant_sets) <= 3:
                    fig = create_venn_diagram(significant_sets, set_names, title)
                    st.pyplot(fig)
                else:
                    fig = create_upset_plot(significant_sets, set_names, title)
                    st.pyplot(fig)
                
                # 共通遺伝子の表示
                common_genes = set.intersection(*significant_sets)
                st.write(f"共通して{dir_name}する遺伝子数: {len(common_genes)}")
                
                if st.checkbox(f"共通{dir_name}遺伝子リストを表示"):
                    st.write(", ".join(sorted(list(common_genes))))
                
                # 共通遺伝子のエクスポート
                if len(common_genes) > 0:
                    csv = pd.DataFrame(sorted(list(common_genes)), columns=['Gene_ID'])
                    csv_data = csv.to_csv(index=False)
                    st.download_button(
                        label=f"共通{dir_name}遺伝子リストをダウンロード",
                        data=csv_data,
                        file_name=f"common_{dir_type}_regulated_genes.csv",
                        mime="text/csv"
                    )
    
    else:  # 散布図による比較
        st.subheader("散布図による比較")
        
        plot_dimension = st.radio(
            "散布図の次元を選択：",
            ["2D散布図", "3D散布図"],
            horizontal=True
        )
        
        highlight = st.checkbox("有意な遺伝子をハイライト表示", value=True)
        
        if plot_dimension == "2D散布図" and len(active_datasets) >= 2:
            col1, col2 = st.columns(2)
            with col1:
                dataset1 = st.selectbox("1つ目のデータセット", [ds['name'] for ds in active_datasets], key='scatter1', index=0)
            with col2:
                dataset2 = st.selectbox("2つ目のデータセット", [ds['name'] for ds in active_datasets], key='scatter2', index=1)
            
            if dataset1 != dataset2:
                ds1 = next(ds for ds in active_datasets if ds['name'] == dataset1)
                ds2 = next(ds for ds in active_datasets if ds['name'] == dataset2)
                
                # 遺伝子セットの選択
                use_union = st.checkbox("全ての遺伝子を表示（存在しない遺伝子は0として扱う）", 
                                      value=False,
                                      help="チェックを外すと共通の遺伝子のみを表示します")
                
                # 遺伝子数の表示
                common_genes = len(ds1['data'].index.intersection(ds2['data'].index))
                total_genes = len(ds1['data'].index.union(ds2['data'].index))
                st.write(f"共通の遺伝子数: {common_genes}")
                st.write(f"総遺伝子数: {total_genes}")
                
                # 値の選択
                value_type = st.radio(
                    "比較する値を選択：",
                    ["Log2 Fold Change", "adjusted P-value (-log10)"],
                    horizontal=True
                )
                
                if value_type == "Log2 Fold Change":
                    fig = create_2d_scatter(
                        ds1['data'], ds2['data'],
                        ds1['lfc_col'], ds2['lfc_col'],
                        dataset1, dataset2,
                        ds1['pval_col'], ds2['pval_col'],
                        highlight,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                else:
                    # adjusted P-valueの場合は-log10に変換
                    ds1['data']['neg_log_padj'] = -np.log10(ds1['data'][ds1['pval_col']])
                    ds2['data']['neg_log_padj'] = -np.log10(ds2['data'][ds2['pval_col']])
                    fig = create_2d_scatter(
                        ds1['data'], ds2['data'],
                        'neg_log_padj', 'neg_log_padj',
                        dataset1, dataset2,
                        ds1['pval_col'], ds2['pval_col'],
                        False,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                
                st.plotly_chart(fig, use_container_width=True)

                # 相関係数の計算と表示
                if value_type == "Log2 Fold Change":
                    if use_union:
                        all_index = ds1['data'].index.union(ds2['data'].index)
                        df1 = ds1['data'].reindex(all_index).fillna(0)
                        df2 = ds2['data'].reindex(all_index).fillna(0)
                    else:
                        common_index = ds1['data'].index.intersection(ds2['data'].index)
                        df1 = ds1['data'].loc[common_index]
                        df2 = ds2['data'].loc[common_index]
                    
                    # NaNを除外して相関係数を計算
                    mask = ~np.isnan(df1[ds1['lfc_col']]) & ~np.isnan(df2[ds2['lfc_col']])
                    correlation = np.corrcoef(
                        df1[ds1['lfc_col']][mask].astype(float),
                        df2[ds2['lfc_col']][mask].astype(float)
                    )[0,1]
                    
                    # 相関係数が計算できた場合のみ表示
                    if not np.isnan(correlation):
                        st.write(f"Log2 Fold Change の相関係数: {correlation:.3f}")
                    else:
                        st.warning("Log2 Fold Change の相関係数を計算できませんでした")
                else:
                    if use_union:
                        all_index = ds1['data'].index.union(ds2['data'].index)
                        df1 = ds1['data'].reindex(all_index).fillna(1)  # P値の欠損は1として扱う
                        df2 = ds2['data'].reindex(all_index).fillna(1)
                    else:
                        common_index = ds1['data'].index.intersection(ds2['data'].index)
                        df1 = ds1['data'].loc[common_index]
                        df2 = ds2['data'].loc[common_index]
                    
                    df1['neg_log_padj'] = -np.log10(df1[ds1['pval_col']])
                    df2['neg_log_padj'] = -np.log10(df2[ds2['pval_col']])
                    
                    correlation = np.corrcoef(
                        df1['neg_log_padj'].values,
                        df2['neg_log_padj'].values
                    )[0,1]
                    st.write(f"-log10(adjusted P-value) の相関係数: {correlation:.3f}")
                
        elif plot_dimension == "3D散布図" and len(active_datasets) >= 3:
            col1, col2, col3 = st.columns(3)
            with col1:
                dataset1 = st.selectbox("X軸のデータセット", [ds['name'] for ds in active_datasets], key='scatter3d1')
            with col2:
                dataset2 = st.selectbox("Y軸のデータセット", [ds['name'] for ds in active_datasets], key='scatter3d2')
            with col3:
                dataset3 = st.selectbox("Z軸のデータセット", [ds['name'] for ds in active_datasets], key='scatter3d3')
            
            if len({dataset1, dataset2, dataset3}) == 3:
                ds1 = next(ds for ds in active_datasets if ds['name'] == dataset1)
                ds2 = next(ds for ds in active_datasets if ds['name'] == dataset2)
                ds3 = next(ds for ds in active_datasets if ds['name'] == dataset3)
                
                # 遺伝子セットの選択
                use_union = st.checkbox("全ての遺伝子を表示（存在しない遺伝子は0として扱う）", 
                                      value=False,
                                      help="チェックを外すと共通の遺伝子のみを表示します")
                
                # 遺伝子数の表示
                common_genes = len(set.intersection(
                    set(ds1['data'].index),
                    set(ds2['data'].index),
                    set(ds3['data'].index)
                ))
                total_genes = len(set.union(
                    set(ds1['data'].index),
                    set(ds2['data'].index),
                    set(ds3['data'].index)
                ))
                st.write(f"共通の遺伝子数: {common_genes}")
                st.write(f"総遺伝子数: {total_genes}")
                

                value_type = st.radio(
                    "比較する値を選択：",
                    ["Log2 Fold Change", "adjusted P-value (-log10)"],
                    horizontal=True
                )
                
                if value_type == "Log2 Fold Change":
                    fig = create_3d_scatter(
                        ds1['data'], ds2['data'], ds3['data'],
                        ds1['lfc_col'], ds2['lfc_col'], ds3['lfc_col'],
                        dataset1, dataset2, dataset3,
                        ds1['pval_col'], ds2['pval_col'], ds3['pval_col'],
                        highlight,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                else:
                    # adjusted P-valueの場合は-log10に変換
                    ds1['data']['neg_log_padj'] = -np.log10(ds1['data'][ds1['pval_col']])
                    ds2['data']['neg_log_padj'] = -np.log10(ds2['data'][ds2['pval_col']])
                    ds3['data']['neg_log_padj'] = -np.log10(ds3['data'][ds3['pval_col']])
                    fig = create_3d_scatter(
                        ds1['data'], ds2['data'], ds3['data'],
                        'neg_log_padj', 'neg_log_padj', 'neg_log_padj',
                        dataset1, dataset2, dataset3,
                        ds1['pval_col'], ds2['pval_col'], ds3['pval_col'],
                        False,
                        lfc_threshold,
                        pval_threshold,
                        use_union
                    )
                
                st.plotly_chart(fig, use_container_width=True)


if __name__ == "__main__":
    main()
