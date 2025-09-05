import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import zipfile
import re
from scipy import stats
from statsmodels.stats.multitest import fdrcorrection

# Function to load data
def load_data(file_path):
    return pd.read_csv(file_path, sep='\t', index_col=0)

# Function to get top regulons for each cell type
def get_top_regulons_per_cell_type(df, n):
    top_regulons = set()
    for column in df.columns:
        top_n = df[column].nlargest(n).index.tolist()
        top_regulons.update(top_n)
    return list(top_regulons)

# Function to parse user input for regulon names
def parse_regulon_input(input_string):
    # Split the input string by commas, spaces, and newlines
    regulons = re.split(r'[,\s\n]+', input_string)
    # Remove empty strings and strip whitespace
    regulons = [r.strip() for r in regulons if r.strip()]
    return regulons

# Function to create clustered heatmap using seaborn
def create_clustered_heatmap(data, is_zscore=False, figsize=(16, 12), show_all_labels=False, cluster_columns=True):
    if is_zscore:
        cmap = "RdBu_r"
        center = 0
    else:
        cmap = "YlOrRd"
        center = None

    # Create the clustermap
    clustergrid = sns.clustermap(data, cmap=cmap, center=center, 
                                 dendrogram_ratio=(0.2, 0.2),
                                 cbar_pos=None,  # Disable default colorbar
                                 figsize=figsize,
                                 col_cluster=cluster_columns)  # Add this line to control column clustering


    # Get the positions of the heatmap and dendrograms
    heatmap_pos = clustergrid.ax_heatmap.get_position()
    dend_row_pos = clustergrid.ax_row_dendrogram.get_position()
    dend_col_pos = clustergrid.ax_col_dendrogram.get_position()

    # Adjust the position of the heatmap
    clustergrid.ax_heatmap.set_position([heatmap_pos.x0, heatmap_pos.y0, 
                                         heatmap_pos.width * 0.9, heatmap_pos.height])

    # Adjust the position of the row dendrogram
    clustergrid.ax_row_dendrogram.set_position([dend_row_pos.x0, heatmap_pos.y0, 
                                                dend_row_pos.width, heatmap_pos.height])

    # Adjust the position of the column dendrogram
    clustergrid.ax_col_dendrogram.set_position([heatmap_pos.x0, dend_col_pos.y0, 
                                                heatmap_pos.width * 0.9, dend_col_pos.height])

    # Add a new axes for the colorbar
    cbar_ax = clustergrid.fig.add_axes([0.92, heatmap_pos.y0, 0.02, heatmap_pos.height])
    
    # Add colorbar to the new axes
    norm = plt.Normalize(vmin=data.min().min(), vmax=data.max().max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    clustergrid.fig.colorbar(sm, cax=cbar_ax)
    
    # Show all labels if requested
    if show_all_labels:
        yticks = np.arange(len(data.index)) + 0.5
        clustergrid.ax_heatmap.set_yticks(yticks)
        clustergrid.ax_heatmap.set_yticklabels(data.index, fontsize=8)
    else:
        # Show only every nth label to avoid overcrowding
        n = max(1, len(data.index) // 20)  # Show at most 20 labels
        yticks = np.arange(0, len(data.index), n) + 0.5
        clustergrid.ax_heatmap.set_yticks(yticks)
        clustergrid.ax_heatmap.set_yticklabels(data.index[::n], fontsize=8)

    # Adjust label positions to be centered
    clustergrid.ax_heatmap.yaxis.set_tick_params(pad=0)
    for label in clustergrid.ax_heatmap.get_yticklabels():
        label.set_verticalalignment('center')

    # Remove tick marks
    clustergrid.ax_heatmap.tick_params(axis='y', which='both', length=0)

    return clustergrid

# 修正: SCENICスタイルでサブセットからTop Nを取得する関数
def get_top_regulons_for_subsets(df, n, subsets):
    top_regulons = set()
    subset_data = df[subsets]
    for column in subset_data.columns:
        # 各細胞タイプ（列）ごとにTop Nを選択
        top_n = subset_data[column].nlargest(n).index.tolist()
        top_regulons.update(top_n)
    return list(top_regulons)


def compare_rss_test_vs_control(rss_data, test_groups, control_groups):
    test_data = rss_data[test_groups]
    control_data = rss_data[control_groups]
    
    # 各グループの平均を計算
    test_mean = test_data.mean(axis=1)
    control_mean = control_data.mean(axis=1)
    
    # 差分の計算
    diff = test_mean - control_mean
    
    # 記述統計
    descriptive_stats = {
        "test_mean": test_mean.mean(),
        "test_std": test_mean.std(),
        "control_mean": control_mean.mean(),
        "control_std": control_mean.std(),
        "mean_diff": diff.mean(),
        "median_diff": diff.median(),
        "std_diff": diff.std(),
        "max_diff": diff.max(),
        "min_diff": diff.min(),
    }
    
    results = {
        "descriptive_stats": descriptive_stats,
        "diff": diff,
        "test_mean": test_mean,
        "control_mean": control_mean
    }
    
    # 統計的検定を実行
    w_statistics = []
    w_pvalues = []
    ks_statistics = []
    ks_pvalues = []
    
    for regulon in rss_data.index:
        w_stat, w_pval = stats.ranksums(test_data.loc[regulon], control_data.loc[regulon])
        ks_stat, ks_pval = stats.ks_2samp(test_data.loc[regulon], control_data.loc[regulon])
        
        w_statistics.append(w_stat)
        w_pvalues.append(w_pval)
        ks_statistics.append(ks_stat)
        ks_pvalues.append(ks_pval)
    
    # FDR補正
    _, w_fdr = fdrcorrection(w_pvalues)
    _, ks_fdr = fdrcorrection(ks_pvalues)
    
    results["wilcoxon"] = {"statistic": w_statistics, "p_value": w_pvalues, "fdr": w_fdr}
    results["ks_test"] = {"statistic": ks_statistics, "p_value": ks_pvalues, "fdr": ks_fdr}
    
    # 効果量（Cohen's d）の計算
    d = (test_mean.mean() - control_mean.mean()) / np.sqrt((test_mean.std()**2 + control_mean.std()**2) / 2)
    results["effect_size"] = d
    
    # 視覚化
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 散布図
    ax1.scatter(control_mean, test_mean, alpha=0.5)
    ax1.set_xlabel("Control Mean RSS")
    ax1.set_ylabel("Test Mean RSS")
    ax1.set_title("Mean RSS: Test vs Control")
    max_val = max(control_mean.max(), test_mean.max())
    ax1.plot([0, max_val], [0, max_val], 'r--')
    
    # ボックスプロット（修正版）
    test_data_flat = test_data.values.flatten()
    control_data_flat = control_data.values.flatten()
    
    boxplot_data = [
        pd.DataFrame({'RSS': test_data_flat, 'Group': 'Test'}),
        pd.DataFrame({'RSS': control_data_flat, 'Group': 'Control'})
    ]
    boxplot_data = pd.concat(boxplot_data, ignore_index=True)
    
    sns.boxplot(x='Group', y='RSS', data=boxplot_data, ax=ax2)
    ax2.set_title("RSS Distribution: Test vs Control")
    ax2.set_ylabel("RSS")
    
    plt.tight_layout()
    
    results["plot"] = fig
    
    return results


# 修正: 統計結果をDataFrameとして作成
def create_statistical_results_df(comparison_results, rss_data):
    regulons = rss_data.index
    
    results_df = pd.DataFrame({
        'Regulon': regulons,
        'Test Mean': comparison_results['test_mean'],
        'Control Mean': comparison_results['control_mean'],
        'Difference': comparison_results['diff']
    })
    
    if 'wilcoxon' in comparison_results:
        results_df['Wilcoxon statistic'] = comparison_results['wilcoxon']['statistic']
        results_df['Wilcoxon p-value'] = comparison_results['wilcoxon']['p_value']
        results_df['Wilcoxon FDR'] = comparison_results['wilcoxon']['fdr']
        
    if 'ks_test' in comparison_results:
        results_df['KS statistic'] = comparison_results['ks_test']['statistic']
        results_df['KS p-value'] = comparison_results['ks_test']['p_value']
        results_df['KS FDR'] = comparison_results['ks_test']['fdr']
    
    return results_df

# Main app

def main():
    st.title('SCENIC GRN viewer')
    st.sidebar.title("Options")

    # File upload
    rss_file = st.file_uploader("Upload rss_regulon_by_cell_type_FULL_TABLE...(_filtered/_unfiltered)", type=['txt'], help="""_filteredには10細胞以上で活性化している信頼性の高いregulonが選別されており、通常はこちらを用いる。""")
    regulon_file = st.file_uploader("Upload scaled_regulon_activity_by_cell_type_FULL_TABLE...(_filtered/_unfiltered)", type=['txt'])

    if rss_file and regulon_file:
        # Load data
        rss_data = load_data(rss_file)
        regulon_data = load_data(regulon_file)

        # Sidebar controls
        st.sidebar.header('Parameters')
        matrix_selection = st.sidebar.radio('Matrix Selection', ['RSS', 'Regulon Activity'])
        if matrix_selection == 'Regulon Activity':
            fig_header = "RegulonActivity_"
        else:
            fig_header = "RSS_"
        selection_method = st.sidebar.radio("Regulon selection method", ["Top N", "Custom"])

        # RSS比較オプションをメインフォームに統合
        compare_rss = st.sidebar.checkbox("Compare RSS between cell types")

        # すべてのオプションを1つのフォームに統合
        with st.sidebar.form("analysis_options_form"):
            st.write("Analysis Options")
            
            # サブセット選択のチェックボックス
            use_subsets = st.checkbox("Use subsets for analysis")
            
            if use_subsets:
                all_groups = rss_data.columns.tolist()
                selected_subsets = st.multiselect("Available subsets", all_groups, default=all_groups)

            if selection_method == "Top N":
                top_n = st.number_input('Number of top regulons per cell type', min_value=1, max_value=200, value=10)
            else:
                custom_regulons = st.text_area("Enter regulon names (without '+', separated by space, comma, or newline)")

            show_all_labels = st.checkbox('Show all regulon names in heatmap', value=False, help="行名をすべて表示する")
            cluster_columns = st.checkbox('Cluster columns (x-axis)', value=True)

            plot_width = st.number_input('Plot width', min_value=4, max_value=24, value=14)
            plot_height = st.number_input('Plot height', min_value=4, max_value=50, value=12)
            download_format = st.radio('Heatmap download format', ['PNG', 'PDF'], index=1)
     
            if compare_rss:
                st.markdown("#### RSS comparison")
                cell_types = rss_data.columns.tolist()
                test_groups = st.multiselect("Select test cell types", cell_types, default=[cell_types[0]])
                control_groups = st.multiselect("Select control cell types", cell_types, default=[cell_types[1]] if len(cell_types) > 1 else [])

            submit_button = st.form_submit_button(label='Apply Settings')

        if submit_button:
            # 以下は既存のコードと同じ処理を続ける
            if use_subsets:
                if selected_subsets:
                    rss_data = rss_data[selected_subsets]
                    regulon_data = regulon_data[selected_subsets]
                else:
                    st.warning("Please select at least one subset for analysis.")
                    return
            
            if selection_method == "Top N":
                top_regulons = get_top_regulons_per_cell_type(rss_data, top_n)
            else:
                input_regulons = parse_regulon_input(custom_regulons)
                if matrix_selection == 'Regulon Activity':
                    top_regulons = [item for item in regulon_data.index if any(k in item for k in input_regulons)]
                else:
                    top_regulons = [item for item in rss_data.index if any(k in item for k in input_regulons)]
                st.sidebar.write(f"Applied {len(top_regulons)} valid regulons")

            # Filter top_regulons to only include those present in data
            if matrix_selection == 'Regulon Activity':    
                valid_regulons = [reg for reg in top_regulons if reg in regulon_data.index]
            else:
                valid_regulons = [reg for reg in top_regulons if reg in rss_data.index]
                
            if not valid_regulons:
                st.error("No valid regulons found. Please check your data or input.")
                return

            # RSS比較の処理
            comparison_results = None
            statistical_results_df = None

            if compare_rss:
                if len(test_groups) == 0 or len(control_groups) == 0:
                    st.warning("Please select at least one test and one control cell type for comparison.")
                else:
                    comparison_results = compare_rss_test_vs_control(rss_data, test_groups, control_groups)
                    statistical_results_df = create_statistical_results_df(comparison_results, rss_data)
                    
                    st.write("### Descriptive Statistics")
                    st.write(f"Test group mean RSS: {comparison_results['descriptive_stats']['test_mean']:.4f} ± {comparison_results['descriptive_stats']['test_std']:.4f}")
                    st.write(f"Control group mean RSS: {comparison_results['descriptive_stats']['control_mean']:.4f} ± {comparison_results['descriptive_stats']['control_std']:.4f}")
                    st.write(f"Mean difference: {comparison_results['descriptive_stats']['mean_diff']:.4f}")
                    st.write(f"Median difference: {comparison_results['descriptive_stats']['median_diff']:.4f}")
                    st.write(f"Effect size (Cohen's d): d = {comparison_results['effect_size']:.4f}")
                    
                    st.write("### Statistical Tests")
                    st.write("Detailed Wilcoxon rank-sum test and Kolmogorov-Smirnov test results (including FDR-corrected p-values) for each regulon are available in the downloadable ZIP file.")
                    
                    st.pyplot(comparison_results['plot'])

                    # Top regulons with largest differences
                    diff_series = comparison_results['diff']
                    top_diff = diff_series.abs().nlargest(10)
                    st.write("### Top 10 Differentially Expressed Regulons")
                    for regulon, diff in top_diff.items():
                        direction = "higher" if diff > 0 else "lower"
                        w_fdr = statistical_results_df.loc[statistical_results_df['Regulon'] == regulon, 'Wilcoxon FDR'].values[0]
                        ks_fdr = statistical_results_df.loc[statistical_results_df['Regulon'] == regulon, 'KS FDR'].values[0]
                        st.write(f"{regulon}: {abs(diff):.4f} ({direction} in Test group, Wilcoxon FDR: {w_fdr:.4e}, KS FDR: {ks_fdr:.4e})")



            # Prepare data for heatmap
            if matrix_selection == 'RSS':
                heatmap_data = rss_data.loc[valid_regulons]
                is_zscore = False
            else:
                heatmap_data = regulon_data.loc[valid_regulons]
                is_zscore = True

            # Create heatmap with dendrograms
            clustergrid = create_clustered_heatmap(heatmap_data, is_zscore, 
                                               figsize=(plot_width, plot_height),
                                               show_all_labels=show_all_labels,
                                               cluster_columns=cluster_columns)
            
            # Convert matplotlib figure to image for display
            buf = io.BytesIO()
            clustergrid.savefig(buf, format='png', dpi=300, bbox_inches='tight')
            buf.seek(0)
            st.image(buf, use_container_width=True)

            # Display matrix
            st.markdown('#### Matrix')
            if matrix_selection == 'RSS':
                display_matrix = rss_data
            else:
                display_matrix = regulon_data
            
            st.dataframe(display_matrix)

            # Download data as ZIP
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                # Add heatmap data
                heatmap_csv = heatmap_data.to_csv().encode()
                zip_file.writestr('heatmap_data.csv', heatmap_csv)
                
                # Add display matrix
                matrix_csv = display_matrix.to_csv().encode()
                zip_file.writestr('matrix_data.csv', matrix_csv)
                
                # Add heatmap image
                img_buf = io.BytesIO()
                clustergrid.savefig(img_buf, format=download_format.lower(), dpi=300, bbox_inches='tight')
                img_buf.seek(0)
                zip_file.writestr(f'{fig_header}heatmap.{download_format.lower()}', img_buf.getvalue())

                # Add statistical results if available
                if statistical_results_df is not None:
                    statistical_results_tsv = statistical_results_df.to_csv(sep='\t', index=False).encode()
                    zip_file.writestr('statistical_results.tsv', statistical_results_tsv)

            zip_buffer.seek(0)
            st.download_button(
                label="Download ZIP",
                data=zip_buffer.getvalue(),
                file_name="grn_analysis_data.zip",
                mime="application/zip"
            )
        else:
            st.info("Please set your analysis options and click 'Apply Settings' to update the analysis.")



if __name__ == '__main__':
    main()