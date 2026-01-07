import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import io
import zipfile
import re
from scipy import stats
from scipy.spatial.distance import jensenshannon
from statsmodels.stats.multitest import fdrcorrection
import scanpy as sc

# Function to load data
def load_data(file_path):
    return pd.read_csv(file_path, sep='\t', index_col=0)

# Accurately reproduce pySCENIC's RSS calculation
def regulon_specificity_scores(auc_mtx, cell_type_series):
    """
    Calculate regulon specificity scores (RSS) exactly as pySCENIC does.

    Parameters
    ----------
    auc_mtx : pd.DataFrame
        AUC matrix with cells as columns and regulons as rows
    cell_type_series : pd.Series
        Cell type annotations with cell IDs as index

    Returns
    -------
    pd.DataFrame
        RSS scores with cell types as rows and regulons as columns
    """
    def rss_score(aucs, labels):
        """
        Calculate RSS using Jensen-Shannon divergence.
        RSS = 1 - JS divergence
        """
        # Normalize to probability distributions
        aucs_norm = aucs / aucs.sum() if aucs.sum() > 0 else aucs
        labels_norm = labels / labels.sum() if labels.sum() > 0 else labels
        return 1.0 - jensenshannon(aucs_norm, labels_norm)

    # Align cell types with AUC matrix columns
    cell_types = cell_type_series.loc[auc_mtx.columns]
    unique_cell_types = cell_types.unique()

    # Initialize RSS matrix
    rss_matrix = pd.DataFrame(
        index=unique_cell_types,
        columns=auc_mtx.index,
        dtype=float
    )

    # Calculate RSS for each regulon and cell type
    for regulon in auc_mtx.index:
        aucs = auc_mtx.loc[regulon].values

        for cell_type in unique_cell_types:
            # Binary indicator: 1 if cell is of this type, 0 otherwise
            labels = (cell_types == cell_type).astype(float).values

            # Calculate RSS
            rss_matrix.loc[cell_type, regulon] = rss_score(aucs, labels)

    return rss_matrix

def calculate_scaled_regulon_activity(auc_mtx, cell_type_series):
    """
    Calculate scaled (Z-score normalized) regulon activity per cell type.

    Parameters
    ----------
    auc_mtx : pd.DataFrame
        AUC matrix with cells as columns and regulons as rows
    cell_type_series : pd.Series
        Cell type annotations with cell IDs as index

    Returns
    -------
    pd.DataFrame
        Z-score normalized activity with cell types as columns and regulons as rows
    """
    # Align cell types with AUC matrix columns
    cell_types = cell_type_series.loc[auc_mtx.columns]

    # Calculate mean activity per cell type
    activity = pd.DataFrame(
        index=auc_mtx.index,
        columns=cell_types.unique()
    )

    for cell_type in cell_types.unique():
        cell_mask = cell_types == cell_type
        activity[cell_type] = auc_mtx.loc[:, cell_mask].mean(axis=1)

    # Z-score normalization across cell types (row-wise)
    scaled_activity = activity.T  # Transpose to match output format
    scaled_activity = (scaled_activity - scaled_activity.mean(axis=0)) / scaled_activity.std(axis=0)

    return scaled_activity

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
                                 col_cluster=cluster_columns)  # Control column clustering


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

# Modified: Function to get Top N from subsets in SCENIC style
def get_top_regulons_for_subsets(df, n, subsets):
    top_regulons = set()
    subset_data = df[subsets]
    for column in subset_data.columns:
        # Select Top N for each cell type (column)
        top_n = subset_data[column].nlargest(n).index.tolist()
        top_regulons.update(top_n)
    return list(top_regulons)


def compare_rss_test_vs_control(rss_data, test_groups, control_groups):
    test_data = rss_data[test_groups]
    control_data = rss_data[control_groups]

    # Calculate mean for each group
    test_mean = test_data.mean(axis=1)
    control_mean = control_data.mean(axis=1)

    # Calculate difference
    diff = test_mean - control_mean

    # Descriptive statistics
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

    # Perform statistical tests
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

    # FDR correction
    _, w_fdr = fdrcorrection(w_pvalues)
    _, ks_fdr = fdrcorrection(ks_pvalues)

    results["wilcoxon"] = {"statistic": w_statistics, "p_value": w_pvalues, "fdr": w_fdr}
    results["ks_test"] = {"statistic": ks_statistics, "p_value": ks_pvalues, "fdr": ks_fdr}

    # Calculate effect size (Cohen's d)
    d = (test_mean.mean() - control_mean.mean()) / np.sqrt((test_mean.std()**2 + control_mean.std()**2) / 2)
    results["effect_size"] = d

    # Visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Scatter plot
    ax1.scatter(control_mean, test_mean, alpha=0.5)
    ax1.set_xlabel("Control Mean RSS")
    ax1.set_ylabel("Test Mean RSS")
    ax1.set_title("Mean RSS: Test vs Control")
    max_val = max(control_mean.max(), test_mean.max())
    ax1.plot([0, max_val], [0, max_val], 'r--')

    # Box plot (modified version)
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


# Modified: Create statistical results as DataFrame
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

    # Option: RSS/Activity recalculation
    st.header("Options")
    st.subheader("RSS/Activity Recalculation")
    recalculate_mode = st.checkbox(
        "Recalculate RSS/Activity with different clustering",
        value=False,
        help="Recalculate RSS and Activity using the original AUC matrix and different clustering information from h5ad file"
    )

    if recalculate_mode:
        st.info("📊 **Recalculation Mode**: Upload AUC matrix and h5ad file to recalculate RSS and activity with different cell type annotations")

        # Upload AUC file
        auc_file = st.file_uploader(
            "Upload AUC_per_cell.txt",
            type=['txt'],
            help="AUC_per_cell.txt file from SCENIC results directory"
        )

        # Upload h5ad file
        h5ad_file = st.file_uploader(
            "Upload h5ad file (with cell metadata)",
            type=['h5ad'],
            help="h5ad file containing cell metadata (clustering information)"
        )

        if auc_file and h5ad_file:
            # Load AUC matrix
            with st.spinner("Loading AUC matrix..."):
                auc_mtx = pd.read_csv(auc_file, sep='\t', index_col=0)
                st.success(f"✓ Loaded AUC matrix: {auc_mtx.shape[0]} regulons × {auc_mtx.shape[1]} cells")

            # Load h5ad file
            with st.spinner("Loading h5ad file..."):
                import tempfile
                import os

                # Save uploaded file temporarily
                with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
                    tmp.write(h5ad_file.read())
                    tmp_path = tmp.name

                try:
                    adata = sc.read_h5ad(tmp_path)
                    st.success(f"✓ Loaded h5ad: {adata.n_obs} cells")
                finally:
                    os.unlink(tmp_path)

            # Select clustering column
            categorical_cols = [col for col in adata.obs.columns
                              if adata.obs[col].dtype.name == 'category'
                              or (adata.obs[col].dtype == 'object' and adata.obs[col].nunique() < 100)]

            if not categorical_cols:
                st.error("❌ No categorical columns found in h5ad metadata!")
                return

            cluster_column = st.selectbox(
                "Select clustering column",
                categorical_cols,
                help="Select the column containing clustering information from metadata"
            )

            # Match cell IDs
            cell_type_series = adata.obs[cluster_column]

            # Match cell IDs between AUC matrix and h5ad
            common_cells = list(set(auc_mtx.columns) & set(cell_type_series.index))

            if len(common_cells) == 0:
                st.error("❌ No matching cell IDs between AUC matrix and h5ad file!")
                st.write("AUC cell IDs (first 5):", auc_mtx.columns[:5].tolist())
                st.write("h5ad cell IDs (first 5):", cell_type_series.index[:5].tolist())
                return

            st.info(f"✓ Found {len(common_cells)} matching cells ({len(common_cells)/len(auc_mtx.columns)*100:.1f}% of AUC matrix)")

            # Calculate RSS and Activity
            if st.button("🔄 Recalculate RSS and Activity", type="primary"):
                with st.spinner("Calculating RSS..."):
                    # Subset to common cells
                    auc_subset = auc_mtx[common_cells]
                    cell_types_subset = cell_type_series.loc[common_cells]

                    # Calculate RSS
                    rss_data = regulon_specificity_scores(auc_subset, cell_types_subset)
                    st.success(f"✓ Calculated RSS: {rss_data.shape}")

                with st.spinner("Calculating scaled activity..."):
                    # Calculate scaled activity
                    regulon_data = calculate_scaled_regulon_activity(auc_subset, cell_types_subset)
                    st.success(f"✓ Calculated scaled activity: {regulon_data.shape}")

                # Use data in recalculation mode
                st.session_state['recalculated_rss'] = rss_data
                st.session_state['recalculated_activity'] = regulon_data
                st.success("✅ Recalculation complete! Scroll down to visualize results.")
        else:
            st.warning("⚠️ Please upload both AUC_per_cell.txt and h5ad file")
            return

        # Use recalculated data
        if 'recalculated_rss' in st.session_state:
            rss_data = st.session_state['recalculated_rss']
            regulon_data = st.session_state['recalculated_activity']
        else:
            return

    else:
        # Normal mode: Load existing RSS/Activity files
        # File upload
        rss_file = st.file_uploader("Upload rss_regulon_by_cell_type_FULL_TABLE...(_filtered/_unfiltered)", type=['txt'], help="""_filtered contains reliable regulons activated in 10+ cells and is typically used.""")
        regulon_file = st.file_uploader("Upload scaled_regulon_activity_by_cell_type_FULL_TABLE...(_filtered/_unfiltered)", type=['txt'])

        if not (rss_file and regulon_file):
            st.info("👆 Upload SCENIC result files to begin")
            return

        # Load data
        rss_data = load_data(rss_file)
        regulon_data = load_data(regulon_file)

    # Visualization processing using rss_data and regulon_data (common to both modes)
    if True:  # When rss_data and regulon_data exist

        # Main panel controls
        st.subheader('Parameters')
        matrix_selection = st.radio('Matrix Selection', ['RSS', 'Regulon Activity'])
        if matrix_selection == 'Regulon Activity':
            fig_header = "RegulonActivity_"
        else:
            fig_header = "RSS_"
        selection_method = st.radio("Regulon selection method", ["Top N", "Custom"])

        # Integrate RSS comparison option into main form
        compare_rss = st.checkbox("Compare RSS between cell types")

        # Subset selection checkbox (placed outside the form)
        use_subsets = st.checkbox("Use subsets for analysis")

        # Integrate all options into one form
        with st.form("analysis_options_form"):
            st.write("Analysis Options")

            if use_subsets:
                all_groups = rss_data.columns.tolist()
                selected_subsets = st.multiselect("Available subsets", all_groups, default=all_groups)

            if selection_method == "Top N":
                top_n = st.number_input('Number of top regulons per cell type', min_value=1, max_value=200, value=10)
            else:
                custom_regulons = st.text_area("Enter regulon names (without '+', separated by space, comma, or newline)")

            show_all_labels = st.checkbox('Show all regulon names in heatmap', value=False, help="Display all row names")
            cluster_columns = st.checkbox('Cluster columns (x-axis)', value=True)

            plot_width = st.number_input('Plot width', min_value=3, max_value=24, value=14)
            plot_height = st.number_input('Plot height', min_value=4, max_value=50, value=12)
            download_format = st.radio('Heatmap download format', ['PNG', 'PDF'], index=1)

            if compare_rss:
                st.markdown("#### RSS comparison")
                cell_types = rss_data.columns.tolist()
                test_groups = st.multiselect("Select test cell types", cell_types, default=[cell_types[0]])
                control_groups = st.multiselect("Select control cell types", cell_types, default=[cell_types[1]] if len(cell_types) > 1 else [])

            submit_button = st.form_submit_button(label='Apply Settings')

        if submit_button:
            # Continue with the same processing as existing code
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
                st.write(f"Applied {len(top_regulons)} valid regulons")

            # Filter top_regulons to only include those present in data
            if matrix_selection == 'Regulon Activity':
                valid_regulons = [reg for reg in top_regulons if reg in regulon_data.index]
            else:
                valid_regulons = [reg for reg in top_regulons if reg in rss_data.index]

            if not valid_regulons:
                st.error("No valid regulons found. Please check your data or input.")
                return

            # RSS comparison processing
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