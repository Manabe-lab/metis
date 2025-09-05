import streamlit as st
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm
from statsmodels.othermod.betareg import BetaModel
import plotly.express as px
from helper_func import remove_sample_num
import os

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []    
    # 最も短い文字列の長さを取得
    min_length = min(len(s) for s in strings)
    # 共通の末尾部分の長さを見つける
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break            
    # 共通の末尾部分が見つからない場合は元のリストを返す
    if suffix_length == 0:
        return strings        
    # 共通の末尾部分を削除して新しいリストを作成
    return [s[:-suffix_length] for s in strings]


def read_file_with_options(file, file_extension):
    """ファイルを読み込む関数"""
    try:
        file_extension = file_extension.lower()
        
        # Excelファイルの場合
        if file_extension in ['xlsx', 'xls']:
            df = pd.read_excel(file, index_col=0)
        # その他のテキストファイル（csv, tsv, txt）の場合
        else:
            try:
                # まずタブ区切りで試す
                df = pd.read_csv(file, sep='\t', index_col=0)
            except:
                # 区切り文字を自動検出
                df = pd.read_csv(file, sep=None, engine='python', index_col=0)
        
        return df
        
    except Exception as e:
        st.error(f'ファイル読み込みエラー: {str(e)}')
        return None

def analyze_data(data, groups, df):
    """クラスターデータの解析を実行"""
    # 列合計を計算して比率に変換
    column_sums = data.sum(axis=0)
    proportions = data / column_sums
    
    # ロジット変換
    logit_transformed = np.log(proportions / (1 - proportions))
    
    # グループインデックスの取得
    unique_groups = np.unique(groups)
    if len(unique_groups) != 2:
        st.error("現在2グループの比較のみ対応しています")
        return None, None, None
        
    group1, group2 = unique_groups
    group1_idx = np.where(groups == group1)[0]
    group2_idx = np.where(groups == group2)[0]
    
    # グループごとにデータを分割
    group1_data = logit_transformed[:, group1_idx]
    group2_data = logit_transformed[:, group2_idx]
    
    # 各クラスターでt検定を実施
    t_stats = []
    p_values = []
    means_group1 = []
    means_group2 = []
    cell_types = df.index.tolist()  # ファイルの行名を使用
    
    for i in range(len(data)):
        t_stat, p_val = stats.ttest_ind(group1_data[i, :], group2_data[i, :])
        t_stats.append(t_stat)
        p_values.append(p_val)
        means_group1.append(np.mean(proportions[i, group1_idx] * 100))
        means_group2.append(np.mean(proportions[i, group2_idx] * 100))
    
    # 多重比較の補正
    rejected_holm, p_adjusted_holm, _, _ = multipletests(p_values, method='holm')
    rejected_bh, p_adjusted_bh, _, _ = multipletests(p_values, method='fdr_bh')
    rejected_by, p_adjusted_by, _, _ = multipletests(p_values, method='fdr_by')
    
    # 結果をデータフレームにまとめる
    results_t = pd.DataFrame({
        'Cell Type': cell_types,
        f'Mean % {group1}': means_group1,
        f'Mean % {group2}': means_group2,
        't-statistic': t_stats,
        'Original p-value': p_values,
        'Holm p-value': p_adjusted_holm,
        'BH FDR p-value': p_adjusted_bh,
        'BY FDR p-value': p_adjusted_by
    })
    

    # ベータ回帰用のデータ準備
    cell_types_beta = []
    group_values = []
    values = []
    
    for i in range(len(data)):
        cell_name = df.index[i]  # ファイルの行名を使用
        for j in range(data.shape[1]):
            cell_types_beta.append(cell_name)  # Cell type {i} の代わりに実際の行名を使用
            group_values.append(groups[j])
            values.append(proportions[i, j])
    
    beta_df = pd.DataFrame({
        'CellType': cell_types_beta,
        'Group': group_values,
        'Proportion': values
    })
    
    # ベータ回帰の結果を格納
    beta_results = []
    
    for cell_type in np.unique(cell_types_beta):
        subset = beta_df[beta_df['CellType'] == cell_type].copy()
        
        # デザイン行列の作成
        X = pd.get_dummies(subset['Group'], drop_first=True)
        X = sm.add_constant(X)
        X = np.asarray(X, dtype=float)
        
        y = np.asarray(subset['Proportion'], dtype=float)
        
        try:
            model = BetaModel(y, X)
            results_model = model.fit()
            
            mean_group1 = subset[subset['Group'] == group1]['Proportion'].mean() * 100
            mean_group2 = subset[subset['Group'] == group2]['Proportion'].mean() * 100
            
            beta_results.append({
                'Cell_Type': cell_type,
                f'Mean_{group1}': mean_group1,
                f'Mean_{group2}': mean_group2,
                'Estimate': results_model.params[1],
                'Std_Error': results_model.bse[1],
                'z_value': results_model.tvalues[1],
                'p_value': results_model.pvalues[1]
            })
        except Exception as e:
            st.error(f"Error in {cell_type}: {str(e)}")
    
    if beta_results:
        results_beta = pd.DataFrame(beta_results)
        
        # FDR補正
        _, results_beta['p_adjusted_BH'] = multipletests(results_beta['p_value'], method='fdr_bh')[:2]
        _, results_beta['p_adjusted_BY'] = multipletests(results_beta['p_value'], method='fdr_by')[:2]
    
    return results_t, results_beta, proportions

def main():
    st.title('Comparison of ratios')
    st.markdown('### logit-transformed T test & beta regression')
    uploaded_file = st.file_uploader(
        "File upload （txt, csv, tsv, excel）",
        type=['txt', 'csv', 'tsv', 'xlsx', 'xls']
    )
    st.write('File format:')
    st.markdown("""
    |  | Cell1 | Cell2 | Cell3 |
    | --- | --- | --- | --- |
    | Sample1 |  |  |  |
    | Sample2 | ||  |

    """)

    transout = st.checkbox('Transpose the data?')

    if uploaded_file is not None:
        try:
            file_extension = uploaded_file.name.split('.')[-1]
            df = read_file_with_options(uploaded_file, file_extension)
            
            if df is not None:
                if transout:
                    df = df.T
                    st.markdown("#### Data are transposed.")

                st.subheader("Uploaded data")
                st.dataframe(df)

                condition = [str(i) for i in df.columns.tolist()] #error防止
                group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
                group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字と_を除く
                group_condition = [x.replace('_', '.') for x in group_condition] #_を.に

                # グループ設定用のデータフレーム作成
                group_df = pd.DataFrame({
                    'Sample': df.columns,
                    'Group': group_condition
                })
                
                with st.form("Set_groups"):
                    edited_group_df = st.data_editor(
                        group_df,
                        column_config={
                            "Sample": st.column_config.TextColumn("Sample", disabled=True),
                            "Group": st.column_config.TextColumn("Group")
                        }
                    )
                    submitted = st.form_submit_button("Submit")
                
                if st.button("Run analysis"):
                    groups = edited_group_df['Group'].values
                    data = df.values
                    
                    results_t, results_beta, proportions = analyze_data(data, groups, df)  # dfを引数として渡す
                    
                    if results_t is not None:
                        group1, group2 = np.unique(groups)
                        
                        st.subheader("t検定の結果")
                        st.dataframe(results_t.style.format({
                            f'Mean % {group1}': "{:.4f}",
                            f'Mean % {group2}': "{:.4f}",
                            't-statistic': "{:.4f}",
                            'Original p-value': "{:.4f}",
                            'Holm p-value': "{:.4f}",
                            'BH FDR p-value': "{:.4f}",
                            'BY FDR p-value': "{:.4f}"
                        }))
                        
                        st.subheader("ベータ回帰の結果")
                        st.dataframe(results_beta.style.format({
                            f'Mean_{group1}': "{:.4f}",
                            f'Mean_{group2}': "{:.4f}",
                            'Estimate': "{:.4f}",
                            'Std_Error': "{:.4f}",
                            'z_value': "{:.4f}",
                            'p_value': "{:.4e}",
                            'p_adjusted_BH': "{:.4e}",
                            'p_adjusted_BY': "{:.4e}"
                        }))
                        
                        # 有意な結果の数を表示
                        st.subheader("有意な結果の数 (α = 0.05)")
                        
                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.write("t test:")
                            st.write(f"Original p-value: {sum(results_t['Original p-value'] < 0.05)}")
                            st.write(f"Holm p-value: {sum(results_t['Holm p-value'] < 0.05)}")
                            st.write(f"BH FDR: {sum(results_t['BH FDR p-value'] < 0.05)}")
                            st.write(f"BY FDR: {sum(results_t['BY FDR p-value'] < 0.05)}")
                        
                        with col2:
                            st.write("beta regression:")
                            st.write(f"Original p-value: {sum(results_beta['p_value'] < 0.05)}")
                            st.write(f"BH FDR: {sum(results_beta['p_adjusted_BH'] < 0.05)}")
                            st.write(f"BY FDR: {sum(results_beta['p_adjusted_BY'] < 0.05)}")
                        
                        # クラスターごとの可視化
                        st.subheader("クラスター比率の分布")
                        
                        # プロット用のデータ準備
                        plot_data = []
                        for i in range(len(data)):
                            cell_name = df.index[i]  # ファイルの行名を使用
                            for j in range(data.shape[1]):
                                plot_data.append({
                                    'Cluster': cell_name,  # f"Cell type {i}"の代わりに実際の行名を使用
                                    'Group': groups[j],
                                    'Proportion': proportions[i,j] * 100
                                })
                        plot_df = pd.DataFrame(plot_data)
                        
                        # 各クラスターについて箱ひげ図を作成
                        for cluster in plot_df['Cluster'].unique():
                            cluster_data = plot_df[plot_df['Cluster'] == cluster]
                            fig = px.box(
                                cluster_data,
                                x='Group',
                                y='Proportion',
                                title=f'{cluster} の比率分布 (%)',
                                points='all'  # 個々のデータ点も表示
                            )
                            
                            fig.update_layout(
                                yaxis_title='Proportion (%)',
                                showlegend=True,
                                boxmode='group'
                            )
                            
                            st.plotly_chart(fig, use_container_width=True)
                        
                        file_name_head = os.path.splitext(uploaded_file.name)[0]
                            

                        # 結果のダウンロード
                        st.subheader("結果のダウンロード")
                                                
                        # メモリ上のバッファにExcelファイルを保存
                        from io import BytesIO
                        buffer = BytesIO()
                        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
                            results_t.to_excel(writer, sheet_name='t_test_results', index=False)
                            results_beta.to_excel(writer, sheet_name='beta_regression_results', index=False)
                            
                        # ダウンロードボタン
                        st.download_button(
                            label="解析結果をダウンロード",
                            data=buffer.getvalue(),
                            file_name=f'{file_name_head}_results.xlsx',
                            mime='application/vnd.ms-excel'
                        )
                
        except Exception as e:
            st.error(f'エラーが発生しました: {str(e)}')

if __name__ == '__main__':
    main()