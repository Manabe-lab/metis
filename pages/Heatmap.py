import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import shutil
import os
import sys
import matplotlib.colors as mcolors
from helper_func import mk_temp_dir
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage as scipy_linkage
from scipy import stats
from sklearn.decomposition import NMF
from statsmodels.stats.multitest import multipletests

# Add font settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42  # Use TrueType fonts
plt.rcParams['ps.fonttype'] = 42   # Use TrueType fonts


st.set_page_config(page_title="Heatmap", page_icon="🌡")

# Functions for P-value overlay feature
def jaro_winkler_similarity(s1, s2):
    """
    Calculate Jaro-Winkler similarity between two strings
    """
    from jellyfish import jaro_winkler_similarity as jw
    return jw(s1, s2)

def find_nes_fdr_pairs(df_columns):
    """
    Detect NES columns and suggest the closest FDR/adj column for each NES column

    Parameters
    ----------
    df_columns : list
        List of dataframe column names

    Returns
    -------
    nes_fdr_pairs : dict
        Mapping of {nes_col: suggested_fdr_col}
    nes_columns : list
        List of detected NES columns
    fdr_columns : list
        List of detected FDR/adj columns
    """
    # Detect columns containing NES
    nes_columns = [col for col in df_columns if 'NES' in col or 'nes' in col.lower()]

    # Detect columns containing FDR or adj
    fdr_columns = [col for col in df_columns
                   if 'FDR' in col or 'fdr' in col.lower() or
                      'adj' in col.lower() or 'p.adj' in col or 'padj' in col.lower()]

    # Suggest the closest FDR column for each NES column
    nes_fdr_pairs = {}

    for nes_col in nes_columns:
        if not fdr_columns:
            nes_fdr_pairs[nes_col] = None
            continue

        # Calculate similarity with each FDR column
        similarities = {}
        for fdr_col in fdr_columns:
            sim = jaro_winkler_similarity(nes_col, fdr_col)
            similarities[fdr_col] = sim

        # Select the column with highest similarity
        best_fdr_col = max(similarities.keys(), key=lambda k: similarities[k])
        nes_fdr_pairs[nes_col] = best_fdr_col

    return nes_fdr_pairs, nes_columns, fdr_columns

def prepare_heatmap_with_pvalues(df, nes_fdr_mapping, nes_nominal_p_mapping=None):
    """
    Prepare heatmap data based on NES and FDR column mappings

    Parameters
    ----------
    df : pd.DataFrame
        Original dataframe
    nes_fdr_mapping : dict
        Mapping of {nes_col: fdr_col}
    nes_nominal_p_mapping : dict, optional
        Mapping of {nes_col: nominal_p_col} (for Stouffer analysis)

    Returns
    -------
    df_values : pd.DataFrame
        Value data for heatmap (NES columns only)
    df_pvalues : pd.DataFrame
        P-value data (corresponding FDR columns)
    df_nominal_pvalues : pd.DataFrame or None
        Nominal P-value data (for Stouffer analysis, only if specified)
    """
    # Extract only NES columns
    nes_cols = list(nes_fdr_mapping.keys())
    df_values = df[nes_cols].copy()

    # Convert value columns to numeric type
    for col in df_values.columns:
        df_values[col] = pd.to_numeric(df_values[col], errors='coerce')

    # Extract corresponding FDR columns and rename
    fdr_cols = [nes_fdr_mapping[col] for col in nes_cols]
    df_pvalues = df[fdr_cols].copy()
    df_pvalues.columns = nes_cols  # Change to same column names as NES columns

    # Convert P-value columns to numeric type (properly handle scientific notation strings)
    for col in df_pvalues.columns:
        df_pvalues[col] = pd.to_numeric(df_pvalues[col], errors='coerce')

    # Process nominal P-value columns
    df_nominal_pvalues = None
    if nes_nominal_p_mapping is not None:
        nominal_p_cols = [nes_nominal_p_mapping[col] for col in nes_cols]
        df_nominal_pvalues = df[nominal_p_cols].copy()
        df_nominal_pvalues.columns = nes_cols  # Change to same column names as NES columns

        # Convert nominal P-value columns to numeric type
        for col in df_nominal_pvalues.columns:
            df_nominal_pvalues[col] = pd.to_numeric(df_nominal_pvalues[col], errors='coerce')

    return df_values, df_pvalues, df_nominal_pvalues

def stouffer_combine_pvalues(pval_matrix, effect_matrix):
    """
    Combine multiple P-values using Stouffer's method (signed) and apply BH correction

    Parameters
    ----------
    pval_matrix : pd.DataFrame
        Nominal P-value DataFrame (rows=genes, columns=conditions)
    effect_matrix : pd.DataFrame
        Effect size (NES, logFC, etc.) DataFrame (rows=genes, columns=conditions)
        Sign indicates positive/negative direction

    Returns
    -------
    combined_pvalues : pd.Series
        Combined P-value for each gene (before BH correction)
    combined_fdr : pd.Series
        Combined FDR for each gene (after BH correction)
    mean_effects : pd.Series
        Mean effect (NES/logFC) for each gene
    """
    combined_pvalues = []
    mean_effects = []

    for idx in pval_matrix.index:
        row_pvals = pval_matrix.loc[idx].values
        row_effects = effect_matrix.loc[idx].values

        # Extract valid P-value and effect pairs
        valid_mask = (~np.isnan(row_pvals)) & (~np.isnan(row_effects)) & (row_pvals > 0) & (row_pvals < 1)
        valid_pvals = row_pvals[valid_mask]
        valid_effects = row_effects[valid_mask]

        if len(valid_pvals) == 0:
            combined_pvalues.append(np.nan)
            mean_effects.append(np.nan)
            continue

        # Convert each P-value to Z-score and multiply by the sign of effect
        z_scores = np.sign(valid_effects) * stats.norm.ppf(1 - valid_pvals)

        # Stouffer's combined Z-score
        combined_z = np.sum(z_scores) / np.sqrt(len(z_scores))

        # Convert to combined P-value (two-tailed test)
        combined_p = 2 * (1 - stats.norm.cdf(np.abs(combined_z)))
        combined_pvalues.append(combined_p)

        # Calculate mean effect (NES/logFC)
        mean_effects.append(np.mean(row_effects))

    # Apply BH correction for multiple testing
    combined_pvalues_series = pd.Series(combined_pvalues, index=pval_matrix.index)

    # Correct only non-NaN values
    valid_indices = ~combined_pvalues_series.isna()
    combined_fdr_series = pd.Series(np.nan, index=pval_matrix.index)

    if valid_indices.sum() > 0:
        valid_pvals = combined_pvalues_series[valid_indices].values
        _, corrected_pvals, _, _ = multipletests(valid_pvals, method='fdr_bh')
        combined_fdr_series[valid_indices] = corrected_pvals

    mean_effects_series = pd.Series(mean_effects, index=pval_matrix.index)

    return combined_pvalues_series, combined_fdr_series, mean_effects_series

def add_significance_dots(ax, data_matrix, pval_matrix, max_log_p=5,
                         p_thresholds=[0.001, 0.01, 0.05],
                         dot_color='white', dot_alpha=1.0, max_dot_size_percent=100,
                         dot_edge_color='black', dot_edge_width=0.1,
                         p_threshold_visualization=1.0, scale_method='log10'):
    """
    Add dots to each heatmap cell according to P-value

    Parameters
    ----------
    ax : matplotlib axes
        Axes of the heatmap (expected to be g.ax_heatmap)
    data_matrix : np.ndarray or pd.DataFrame
        Heatmap value data
    pval_matrix : np.ndarray or pd.DataFrame
        P-value data (same shape as data_matrix)
    max_log_p : float
        Maximum value of -log10(p) (for scaling)
    p_thresholds : list
        P-value thresholds to display in legend
    dot_color : str
        Dot fill color (default: 'white')
    dot_alpha : float
        Dot transparency
    max_dot_size_percent : int
        Maximum dot size as percentage of cell height (default: 100)
    dot_edge_color : str
        Dot edge (outline) color (default: 'black')
    dot_edge_width : float
        Dot edge (outline) width (default: 0.1)
    p_threshold_visualization : float
        Maximum P-value threshold for dot display (default: 1.0)
        P-values larger than this value will not show dots
    scale_method : str
        Scaling method: 'log10' or 'linear'

    Notes
    -----
    If P=0 values are detected, a warning will be displayed and they will be drawn as maximum-sized dots.
    """
    from matplotlib.lines import Line2D

    # Convert DataFrame to ndarray
    if isinstance(data_matrix, pd.DataFrame):
        data_matrix = data_matrix.values
    if isinstance(pval_matrix, pd.DataFrame):
        pval_matrix = pval_matrix.values

    # Detect and replace P=0 (treat P=0 as maximum significance)
    zero_p_mask = (pval_matrix == 0) & (~np.isnan(pval_matrix))
    zero_p_count = np.sum(zero_p_mask)
    if zero_p_count > 0:
        # Get minimum non-zero P-value in data
        non_zero_pvals = pval_matrix[(pval_matrix > 0) & (~np.isnan(pval_matrix))]
        if len(non_zero_pvals) > 0:
            min_nonzero_p = np.min(non_zero_pvals)
            replacement_p = min_nonzero_p  # Display with same size as minimum non-zero P-value
        else:
            replacement_p = 1e-10  # Fallback

        # Display warning
        st.warning(f"⚠️ {zero_p_count} cells with P-value=0 detected. These will be displayed with the same size as the minimum non-zero P-value ({replacement_p:.2e}).")
        # Replace P=0 with minimum non-zero P-value (display as same-sized dots)
        pval_matrix = pval_matrix.copy()  # Copy to avoid modifying original data
        pval_matrix[zero_p_mask] = replacement_p

    # Calculate cell height (in points)
    fig = ax.get_figure()
    bbox = ax.get_position()
    fig_height_inches = fig.get_figheight()
    ax_height_inches = bbox.height * fig_height_inches
    ax_height_points = ax_height_inches * 72  # 1 inch = 72 points

    n_rows = data_matrix.shape[0]
    cell_height_points = ax_height_points / n_rows

    # Limit maximum dot size to specified percentage of cell height
    max_dot_size = (max_dot_size_percent / 100.0) * cell_height_points

    # Get minimum P-value in actual data (used for both scaling methods)
    min_p_in_data = None
    # Extract valid P-values from numpy array
    valid_mask = (~np.isnan(pval_matrix)) & (pval_matrix < 1) & (pval_matrix > 0) & (pval_matrix <= p_threshold_visualization)
    valid_pvals = pval_matrix[valid_mask]
    if len(valid_pvals) > 0:
        min_p_in_data = np.min(valid_pvals)

    # For log10 scaling, dynamically calculate max_log_p based on minimum P-value in data
    dynamic_max_log_p = max_log_p  # Default value
    if scale_method == 'log10' and min_p_in_data is not None and min_p_in_data > 0:
        dynamic_max_log_p = -np.log10(min_p_in_data)

    # Draw dot in each cell
    for i in range(data_matrix.shape[0]):
        for j in range(data_matrix.shape[1]):
            p_val = pval_matrix[i, j]

            # Skip invalid P-values or P-values larger than threshold
            if np.isnan(p_val) or p_val >= 1 or p_val <= 0 or p_val > p_threshold_visualization:
                continue

            # Calculate value according to scaling method
            if scale_method == 'log10':
                # Calculate -log10(p) (dynamically scale based on minimum P-value in data)
                transformed_val = -np.log10(p_val)
                normalized = transformed_val / dynamic_max_log_p if dynamic_max_log_p > 0 else 0
                normalized = max(0.0, min(1.0, normalized))  # Clip to 0-1
            else:  # linear
                # Linear scaling: scale between p_threshold_visualization (minimum size) and min_p_in_data (maximum size)
                if min_p_in_data is not None and p_threshold_visualization > min_p_in_data:
                    # Smaller p_val results in larger normalized
                    normalized = (p_threshold_visualization - p_val) / (p_threshold_visualization - min_p_in_data)
                    # Clip to 0-1 range
                    normalized = max(0.0, min(1.0, normalized))
                else:
                    normalized = 1.0  # Fallback (all maximum size)

            # Continuous scaling (range 0.20 to 1.0)
            scaled_size = 0.20 + 0.80 * normalized
            final_size = max_dot_size * scaled_size  # Scale based on max_dot_size

            # Draw dot (adjust to seaborn heatmap coordinate system)
            ax.plot(j + 0.5, i + 0.5, 'o',
                   markersize=final_size,
                   markerfacecolor=dot_color,
                   markeredgecolor=dot_edge_color,
                   markeredgewidth=dot_edge_width,
                   alpha=dot_alpha)

    # Add legend
    legend_elements = []

    for p in sorted(p_thresholds):
        # Calculate value according to scaling method (same logic as heatmap)
        if scale_method == 'log10':
            transformed_val = -np.log10(p)
            normalized = transformed_val / dynamic_max_log_p if dynamic_max_log_p > 0 else 0
            normalized = max(0.0, min(1.0, normalized))  # Clip to 0-1
        else:  # linear
            # Same normalization as heatmap (scale between p_threshold_visualization and min_p_in_data)
            if min_p_in_data is not None and p_threshold_visualization > min_p_in_data:
                normalized = (p_threshold_visualization - p) / (p_threshold_visualization - min_p_in_data)
                # Clip to 0 if legend P-value is larger than threshold, clip to 1 if too small
                normalized = max(0.0, min(1.0, normalized))
            else:
                normalized = 1.0  # Fallback

        scaled_size = 0.20 + 0.80 * normalized
        # Same size calculation as heatmap
        legend_size = max_dot_size * scaled_size

        legend_elements.append(
            Line2D([0], [0], marker='o', color='w',
                   markerfacecolor=dot_color,
                   markeredgecolor=dot_edge_color,
                   markeredgewidth=dot_edge_width,
                   alpha=dot_alpha,
                   markersize=legend_size,
                   label=f'p={p}')
        )

    # Place legend (positioned outside upper right of figure)
    fig = ax.get_figure()
    p_legend = fig.legend(
        handles=legend_elements,
        title="P-value",
        loc='upper left',
        bbox_to_anchor=(1.02, 1.0),
        fontsize=9,
        frameon=True,
        title_fontsize=10,
        bbox_transform=fig.transFigure
    )

    return ax

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t', header = None).encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header= 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
    return df_c

if 'filename_add' not in globals(): # Keep previous data when starting over
 #   st.write('file name kept')
    filename_add = ""


#https://discuss.streamlit.io/t/dynamically-created-multiple-checkbox/18273/2
def checkbox_container(data):
#    st.header('Select columns or rows')
#    new_data = st.text_input('Enter country Code to add')
    cols = st.columns(5)
#    if cols[0].button('Add Coutry'):
#        dummy_data.append(new_data)
    if cols[0].button('Select All'):
        for i in data:
            st.session_state['dynamic_checkbox_' + i] = True
        st.rerun()
    if cols[1].button('UnSelect All'):
        for i in data:
            st.session_state['dynamic_checkbox_' + i] = False
        st.rerun()
    for i in data:
        st.checkbox(i, key='dynamic_checkbox_' + i)

def get_selected_checkboxes():
    return [i.replace('dynamic_checkbox_','') for i in st.session_state.keys() if i.startswith('dynamic_checkbox_') and st.session_state[i]]


def create_custom_cmap(colors, name='custom'):
    if len(colors) == 2:
        # For two colors, create a linear gradient
        return mcolors.LinearSegmentedColormap.from_list(name, colors, N=256)
    elif len(colors) == 3:
        # For three colors, create two segments
        cmap1 = mcolors.LinearSegmentedColormap.from_list("cmap1", colors[:2], N=128)
        cmap2 = mcolors.LinearSegmentedColormap.from_list("cmap2", colors[1:], N=128)
        # Combine the two segments
        newcolors = np.vstack((cmap1(np.linspace(0, 1, 128)),
                               cmap2(np.linspace(0, 1, 128))))
        return mcolors.ListedColormap(newcolors, name=name)
    else:
        raise ValueError("Only 2 or 3 colors are supported for custom colormap creation.")

# Categorize Matplotlib colormaps by category
colormap_categories = {
    'Sequential': [
        'viridis', 'plasma', 'inferno', 'magma', 'cividis',
        'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
        'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
        'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn'
    ],
    'Sequential (2)': [
        'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
        'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
        'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper'
    ],
    'Seaborn': [
        'rocket', 'mako', 'flare', 'crest', 'vlag', 'icefire'
    ],
    'Diverging': [
        'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
        'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic'
    ],
    'Cyclic': ['twilight', 'twilight_shifted', 'hsv'],
    'Qualitative': [
        'Pastel1', 'Pastel2', 'Paired', 'Accent', 'Dark2',
        'Set1', 'Set2', 'Set3', 'tab10', 'tab20', 'tab20b', 'tab20c'
    ],
    'Miscellaneous': [
        'flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
        'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
        'gist_rainbow', 'rainbow', 'jet', 'turbo', 'nipy_spectral',
        'gist_ncar'
    ]
}

def show_colormap(cmap_name):
    gradient = np.linspace(0, 1, 256)
    gradient = np.vstack((gradient, gradient))
    fig, ax = plt.subplots(figsize=(10, 1))
    ax.imshow(gradient, aspect='auto', cmap=plt.get_cmap(cmap_name))
    ax.set_axis_off()
    st.pyplot(fig)


@st.cache_data
def perform_clustering_computation(data, method, metric):
    st.write(method)
    st.write(metric)
    # Perform data normalization and preprocessing here (if needed)

    # Row clustering
    row_distances = pdist(data, metric)
    row_linkage = scipy_linkage(row_distances, method)

    # Column clustering
    col_distances = pdist(data.T, metric)
    col_linkage = scipy_linkage(col_distances, method)

    return row_linkage, col_linkage, row_distances, col_distances

def plot_clustermap(df, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c, xticklabels, yticklabels, annot, fmt, linewidths, linecolor, py_x_size, py_y_size, x_font_size, y_font_size):
    g = sns.clustermap(
        df,
        row_linkage=row_linkage,
        col_linkage=col_linkage,
        center=v_center,
        cmap=cmap,
        vmin=v_min,
        vmax=v_max,
        row_cluster=y_c,
        col_cluster=x_c,
        xticklabels=xticklabels,
        yticklabels=yticklabels,
        annot=annot,
        fmt=fmt,
        linewidths=linewidths,
        linecolor=linecolor,
        figsize=(py_x_size, py_y_size)
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=x_font_size)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
    return g

def clean_column_name(col):
    if isinstance(col, (int, float)):
        # For integers, convert to string as is
        if col.is_integer():
            return str(int(col))
        # For decimals, remove trailing zeros
        return f'{col:g}'
    elif isinstance(col, str):
        # For strings, try to convert to number if possible
        try:
            num = float(col)
            if num.is_integer():
                return str(int(num))
            return f'{num:g}'
        except ValueError:
            # If cannot convert to number, return as is
            return col
    else:
        return str(col)

st.markdown("### Heatmap")
st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")
# Save to temp directory
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir, res_dir = mk_temp_dir("Heatmap")
    st.session_state.temp_dir = temp_dir
else:
    temp_dir = st.session_state.temp_dir
    temp_dir, res_dir = mk_temp_dir("Heatmap", temp_dir)


use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Data in memory:")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
 #       st.write(st.session_state.uploaded_file_name)
        file_name_head = st.session_state.uploaded_file_name
        # Homer support
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list(): # When Row_name is included
            df = df.set_index('Row_name')
            df.index.name = "Gene"


if use_upload == 'Yes':
    input_file_type = st.radio(
        "Data format:",
        ('tsv','csv', 'excel'))
    #Gene_column = 0
    genome = st.checkbox("Genome occupancy data (eg., Homer's annotatePeaks.pl output)?")
    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])

    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, header = None, index_col = None)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', header=None, index_col = None)
        else:
            df = read_excel(uploaded_file, index_col = None, header = None)
        st.write("Uploaded data:")
        st.write(df.head(3))
        st.write('Data Dimension: '+str(df.shape))

        st.markdown("###### Data format should be genes as rows and sample as columns.")
        st.markdown("""
    |  | Sample1 | Sample2 |
    | --- | --- | --- |
    | Gene1 |  |  |
    | Gene2 | ||

    """)
        st.write("    ")

        transpose_df = st.checkbox('Transpose the data?')

        if transpose_df:
            df = df.T
        if not genome:
            df.columns = df.iloc[0,:].tolist()   # Determine columns after transpose
        else: # For homer annotate, add numbers to the beginning in order of appearance
            org_col = df.iloc[0,:].tolist()
            from collections import defaultdict
            # Dictionary to track column name occurrences
            column_counts = defaultdict(int)
            # Generate new column names
            new_columns = [org_col[0]]  # First column name as is
            for col in org_col[1:]:
                col = clean_column_name(col)
                #col = str(col).rstrip('.0')
                column_counts[col] += 1
                if column_counts[col] > 1:
                    new_col_name = f"{column_counts[col]}_{col}"
                else:
                    new_col_name = col
                new_columns.append(new_col_name)

            # Set new column names
            df.columns = new_columns

        df = df.drop(0, axis = 0) # Use first row as column name and remove
        content = df.columns.tolist()
        # Handle data where R's (0,0) is blank ==============================
        if isinstance(content[0], float) and np.isnan(content[0]):
            st.write("0,0isnan")
            content[0] = "Gene"  # Change to "Gene" instead of "NaN"
            df.columns = content
        #========================================================

        Gene_column = content[0]
        if "Annotation/Divergence" in content:
              # Convert column names
            search_word = '([^\ \(]+).*'

            for i in range(1, len(content)):
                match = re.search(search_word, content[i])
                if match:
                    content[i] = match.group(1).replace(' ', '_')
            df.columns = content # Change name temporarily
            df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel support

            pattern = "([^|]*)"
            repatter = re.compile(pattern)
            f_annotation = lambda x: repatter.match(x).group(1)
            df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
    #        df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
            # Remove everything before annotation/divergence
            df = df.loc[:,'Annotation/Divergence':]
            content = df.columns.tolist()
            content[0] = 'Gene'
            df.columns = content
            Gene_column = "Gene"
            st.write("Converted Annotation/Divergence to gene symbols.")

        elif "Gene" in content:
            Gene_column =  "Gene"
        else:
            Gene_column =  st.selectbox('Select gene name column', content)

        df = df.set_index(Gene_column)
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        st.session_state.uploaded_file_name = file_name_head
        st.session_state.df = df

    else:
        sys.exit(1)
   ##### file upload

if df is not None:
    # First check if P-value overlay mode
    st.markdown('---')
    st.markdown("##### P-value overlay mode?")
    st.markdown("###### For GSEA results or datasets with value columns (e.g., NES) and p-value columns (e.g., FDR)")
    use_pvalue_overlay_check = st.checkbox('Enable P-value overlay with significance dots', value=False)

    # Numeric conversion processing
    if use_pvalue_overlay_check:
        # In P-value overlay mode, apply numeric conversion only to selected columns later
        st.info("💡 In P-value overlay mode, only selected columns will be converted to numeric values.")
        # Do not delete rows containing NA (string columns may exist)
    else:
        # In normal mode, convert all columns to numeric
        nonzero = st.checkbox('Remove all zero genes?', value=True)
        if nonzero:
            df = df.loc[~(df==0).all(axis=1)] # Remove rows where all are 0
        df = df.dropna(how='any', axis=0)
        df = df.astype(float) # Strings may remain

    st.write(df.head(3))

    # Check and fix duplicate column names
    if len(df.columns) != len(set(df.columns)):
        duplicate_cols = [col for col in df.columns if df.columns.tolist().count(col) > 1]
        unique_duplicates = list(set(duplicate_cols))
        st.warning(f"⚠️ Duplicate column names detected: {unique_duplicates}")
        st.info("💡 Renaming duplicate columns to make them unique...")

        # Make column names unique
        cols = pd.Series(df.columns)
        for dup in unique_duplicates:
            dup_indices = [i for i, x in enumerate(df.columns) if x == dup]
            for idx, pos in enumerate(dup_indices):
                if idx > 0:  # Keep first occurrence as is, add numbers to subsequent ones
                    cols[pos] = f"{dup}_{idx}"
        df.columns = cols
        st.write("Renamed columns:")
        st.write(df.columns.tolist())

    # P-value overlay mode flow (checkbox already displayed above)
    df_pvalues = None
    df_nominal_pvalues = None
    use_stouffer = False
    if use_pvalue_overlay_check:
        st.markdown("---")
        st.markdown("##### Configure value and p-value column pairs")

        # Get column name list
        all_columns = df.columns.tolist()

        # Auto-detect NES and FDR columns (use as recommended candidates)
        nes_fdr_pairs, nes_candidates_auto, fdr_candidates_auto = find_nes_fdr_pairs(all_columns)

        if not nes_candidates_auto:
            st.info("💡 No columns containing 'NES' found. All columns will be available for selection.")
            nes_candidates_auto = []

        if not fdr_candidates_auto:
            st.info("💡 No columns containing 'FDR' or 'adj' found. All columns will be available for selection.")
            fdr_candidates_auto = []

        if nes_candidates_auto:
            st.info(f"🔍 Auto-detected {len(nes_candidates_auto)} value columns (containing 'NES') and {len(fdr_candidates_auto)} p-value columns (containing 'FDR/adj')")

        # Select number of pairs
        st.markdown("###### Step 1: How many column pairs?")

        # Calculate max_value (minimum 1)
        max_pairs = max(1, len(all_columns) // 2)

        # Set default value
        if 'num_pairs' not in st.session_state:
            default_pairs = len(nes_candidates_auto) if nes_candidates_auto else 1
            # Adjust not to exceed max_pairs
            st.session_state.num_pairs = min(default_pairs, max_pairs)

        # If existing num_pairs exceeds max_pairs, adjust
        if st.session_state.num_pairs > max_pairs:
            st.session_state.num_pairs = max_pairs

        num_pairs = st.number_input(
            "Number of value-pvalue column pairs:",
            min_value=1,
            max_value=max_pairs,
            value=st.session_state.num_pairs,
            help="How many value columns (e.g., NES) do you want to visualize?"
        )
        st.session_state.num_pairs = num_pairs

        st.markdown("---")
        st.markdown("###### Step 2: Configure each column pair")
        st.markdown("*Select a value column and its corresponding p-value column*")

        # Keep selections in session state
        if 'nes_fdr_mapping' not in st.session_state:
            st.session_state.nes_fdr_mapping = nes_fdr_pairs

        # Select value and P-value columns for each pair
        selected_pairs = {}
        already_selected_nes = []  # Track already selected value columns
        already_selected_fdr = []  # Track already selected P-value columns

        for idx in range(num_pairs):
            st.markdown(f"##### Pair {idx+1}")

            col1, col2 = st.columns(2)

            with col1:
                # Select value column
                st.markdown("**Value column** (e.g., NES)")

                # If NES candidates found, display them preferentially
                if nes_candidates_auto:
                    # Place NES candidates first, others after
                    other_columns = [col for col in all_columns if col not in nes_candidates_auto]
                    value_column_options = nes_candidates_auto + other_columns

                    # Set default value (select from unselected NES columns)
                    # List of NES candidates excluding already selected columns
                    available_nes_candidates = [col for col in nes_candidates_auto if col not in already_selected_nes]

                    if available_nes_candidates:
                        # If unselected NES columns exist, select first unselected column
                        default_nes = available_nes_candidates[0]
                        default_nes_index = value_column_options.index(default_nes)
                    elif idx < len(nes_candidates_auto):
                        # If all NES columns selected, select in order
                        default_nes = nes_candidates_auto[idx % len(nes_candidates_auto)]
                        default_nes_index = value_column_options.index(default_nes)
                    else:
                        default_nes_index = 0

                    st.caption(f"💡 {len(nes_candidates_auto)} NES-containing columns listed first")
                else:
                    # If no NES candidates, display all columns
                    value_column_options = all_columns
                    default_nes_index = min(idx, len(all_columns) - 1)

                selected_nes = st.selectbox(
                    "Value column:",
                    options=value_column_options,
                    index=default_nes_index,
                    key=f"nes_select_{idx}",
                    label_visibility="collapsed",
                    help="Select the column containing values (e.g., NES, log2FC, t-statistic)"
                )

                # Record selected column (to exclude in next pair)
                already_selected_nes.append(selected_nes)

            with col2:
                # Select P-value column
                st.markdown("**P-value column** (e.g., FDR)")

                # Extract columns containing P-value related strings as candidates (exclude already selected columns)
                pvalue_keywords = ['p', 'P', 'fdr', 'FDR', 'adj', 'padj', 'p.adj', 'qvalue', 'q.value', 'pval', 'p_val', 'p-val']
                pvalue_candidates = [col for col in all_columns
                                    if any(keyword in col for keyword in pvalue_keywords)
                                    and col not in already_selected_fdr]

                # If P-value candidates found, auto-select most similar one based on string distance
                if pvalue_candidates:
                    # Calculate string distance with selected value column
                    similarities = {}
                    for pval_col in pvalue_candidates:
                        sim = jaro_winkler_similarity(selected_nes, pval_col)
                        similarities[pval_col] = sim

                    # Recommend (auto-select) column with highest similarity
                    suggested_fdr = max(similarities.keys(), key=lambda k: similarities[k])

                    # Place recommended column first in candidate list
                    pvalue_options = [suggested_fdr] + [col for col in pvalue_candidates if col != suggested_fdr]
                    # Make remaining columns also selectable (unselected columns not in candidates)
                    other_cols = [col for col in all_columns
                                 if col not in pvalue_candidates and col not in already_selected_fdr]
                    pvalue_options.extend(other_cols)

                    default_fdr_index = 0  # Recommended column first (auto-selected)

                    st.caption(f"💡 Auto-selected: **{suggested_fdr}** (similarity: {similarities[suggested_fdr]:.2f})")
                else:
                    # If no P-value candidates, select from unselected columns
                    available_cols = [col for col in all_columns if col not in already_selected_fdr]
                    if available_cols:
                        pvalue_options = available_cols
                        default_fdr_index = 0
                    else:
                        # If all columns selected, select from all columns
                        pvalue_options = all_columns
                        default_fdr_index = 0

                selected_fdr = st.selectbox(
                    "P-value column:",
                    options=pvalue_options,
                    index=default_fdr_index,
                    key=f"fdr_select_{idx}",
                    label_visibility="collapsed",
                    help="Select the column containing p-values (e.g., FDR, adj.P.Val, padj)"
                )

                # Record selected P-value column (to exclude in next pair)
                already_selected_fdr.append(selected_fdr)

            selected_pairs[selected_nes] = selected_fdr

            # Display pair confirmation
            st.caption(f"✓ Value: **{selected_nes}** ↔ P-value: **{selected_fdr}**")

            if idx < num_pairs - 1:  # Separator line if not last element
                st.markdown("---")

        # Duplicate check
        value_columns = list(selected_pairs.keys())
        if len(value_columns) != len(set(value_columns)):
            st.error("⚠️ ERROR: Same value column selected multiple times! Each value column should be unique.")
            st.stop()

        # Save selections
        st.session_state.nes_fdr_mapping = selected_pairs

        # Display summary of selected pairs
        st.markdown("---")
        st.success(f"✓ {len(selected_pairs)} column pairs configured successfully")
        with st.expander("📋 View selected column pairs", expanded=False):
            for idx, (nes_col, fdr_col) in enumerate(selected_pairs.items(), 1):
                st.write(f"{idx}. Value: **{nes_col}** ↔ P-value: **{fdr_col}**")

        # Data preview
        with st.expander("📊 Preview selected columns", expanded=False):
            preview_cols = []
            for nes_col, fdr_col in selected_pairs.items():
                preview_cols.extend([nes_col, fdr_col])
            # Remove duplicates (preserve order)
            seen = set()
            preview_cols_unique = []
            for col in preview_cols:
                if col not in seen:
                    seen.add(col)
                    preview_cols_unique.append(col)

            if len(preview_cols) != len(preview_cols_unique):
                st.warning(f"⚠️ Some columns are selected multiple times. Showing unique columns only.")

            st.dataframe(df[preview_cols_unique].head(5))

        # Stouffer's method settings
        st.markdown("---")
        st.markdown("##### Stouffer's method for combining p-values")
        use_stouffer = st.checkbox(
            'Use Stouffer\'s method to combine multiple p-values',
            value=False,
            help='Combine p-values from multiple conditions and add a single Combined column.'
        )

        nes_nominal_p_mapping = None
        if use_stouffer:
            st.info("💡 Stouffer's method will combine nominal p-values across conditions and apply BH correction")
            st.markdown("###### Configure nominal p-value columns for Stouffer's method")
            st.markdown("*Select the nominal p-value column for each value column*")

            # Select nominal P-value columns
            selected_nominal_p_pairs = {}
            already_selected_nominal_p = []

            for idx, (nes_col, fdr_col) in enumerate(selected_pairs.items(), 1):
                st.markdown(f"**{idx}. {nes_col}**")

                # Extract columns containing nominal P-value related strings as candidates
                nominal_p_keywords = ['NOM', 'nom', 'nominal', 'Nominal', 'p.val', 'pval', 'PValue']
                nominal_p_candidates = [col for col in all_columns
                                        if any(keyword in col for keyword in nominal_p_keywords)
                                        and col not in already_selected_nominal_p]

                # If nominal P candidates found, auto-select most similar one based on string distance
                if nominal_p_candidates:
                    similarities = {}
                    for pval_col in nominal_p_candidates:
                        sim = jaro_winkler_similarity(nes_col, pval_col)
                        similarities[pval_col] = sim

                    suggested_nominal_p = max(similarities.keys(), key=lambda k: similarities[k])

                    # Place recommended column first in candidate list
                    nominal_p_options = [suggested_nominal_p] + [col for col in nominal_p_candidates if col != suggested_nominal_p]
                    # Make remaining columns also selectable
                    other_cols = [col for col in all_columns
                                 if col not in nominal_p_candidates and col not in already_selected_nominal_p]
                    nominal_p_options.extend(other_cols)

                    default_nominal_p_index = 0
                    st.caption(f"💡 Auto-selected: **{suggested_nominal_p}** (similarity: {similarities[suggested_nominal_p]:.2f})")
                else:
                    # If no nominal P candidates
                    available_cols = [col for col in all_columns if col not in already_selected_nominal_p]
                    if available_cols:
                        nominal_p_options = available_cols
                        default_nominal_p_index = 0
                    else:
                        nominal_p_options = all_columns
                        default_nominal_p_index = 0

                selected_nominal_p = st.selectbox(
                    f"Nominal p-value column for {nes_col}:",
                    options=nominal_p_options,
                    index=default_nominal_p_index,
                    key=f"nominal_p_select_{idx}",
                    help="Select the column containing nominal p-values (not FDR-corrected)"
                )

                already_selected_nominal_p.append(selected_nominal_p)
                selected_nominal_p_pairs[nes_col] = selected_nominal_p

                st.caption(f"✓ Value: **{nes_col}** ↔ Nominal P: **{selected_nominal_p}**")

                if idx < len(selected_pairs):
                    st.markdown("---")

            # Display summary of selected nominal P columns
            st.success(f"✓ {len(selected_nominal_p_pairs)} nominal p-value columns configured for Stouffer's method")
            with st.expander("📋 View nominal p-value column pairs", expanded=False):
                for idx, (nes_col, nominal_p_col) in enumerate(selected_nominal_p_pairs.items(), 1):
                    st.write(f"{idx}. Value: **{nes_col}** ↔ Nominal P: **{nominal_p_col}**")

            nes_nominal_p_mapping = selected_nominal_p_pairs

        # Prepare data
        try:
            df, df_pvalues, df_nominal_pvalues = prepare_heatmap_with_pvalues(
                df, st.session_state.nes_fdr_mapping, nes_nominal_p_mapping
            )

            # Convert selected columns to numeric
            st.write("Converting selected columns to numeric...")
            for col in df.columns:
                try:
                    df[col] = pd.to_numeric(df[col], errors='coerce')
                except:
                    st.warning(f"Could not convert column '{col}' to numeric")

            for col in df_pvalues.columns:
                try:
                    df_pvalues[col] = pd.to_numeric(df_pvalues[col], errors='coerce')
                except:
                    st.warning(f"Could not convert p-value column '{col}' to numeric")

            if df_nominal_pvalues is not None:
                for col in df_nominal_pvalues.columns:
                    try:
                        df_nominal_pvalues[col] = pd.to_numeric(df_nominal_pvalues[col], errors='coerce')
                    except:
                        st.warning(f"Could not convert nominal p-value column '{col}' to numeric")

            # Remove NA rows
            original_len = len(df)
            df = df.dropna(how='any', axis=0)
            df_pvalues = df_pvalues.loc[df.index]  # Align to same rows
            if df_nominal_pvalues is not None:
                df_nominal_pvalues = df_nominal_pvalues.loc[df.index]  # Align to same rows

            if len(df) < original_len:
                st.info(f"ℹ️ Removed {original_len - len(df)} rows with NA values after numeric conversion")

            st.success(f"✓ Prepared {len(selected_pairs)} value-pvalue pairs for heatmap with significance dots")
            st.write("Heatmap will show:", df.columns.tolist())
            st.write(f"Data shape: {df.shape[0]} genes × {df.shape[1]} conditions")

        except Exception as e:
            st.error(f"Error preparing data: {str(e)}")
            st.stop()

    st.markdown('---')
    st.markdown("##### Filter and transform data?")
    calc_z = False
    center0_z = False  # Set to True for Z-score
    howlog = 'No'
    Manip = st.checkbox('minip', label_visibility = 'collapsed')
    if Manip:
        f_inf = -float('inf')
        p_inf = float('inf')
        min_val = f_inf
        max_val = f_inf
        high_min_val = p_inf
        high_max_val = p_inf
        delta_val = 1
        fold_val = 1
        min_variance = 0
        top_n = p_inf
        more_filt = st.checkbox('Additional filtering options (e.g., minimal values, FC...)')
        if more_filt:
            min_val =  float(st.text_input("All values of each gene are larger than",  value=f_inf))
            max_val =  float(st.text_input("Max value is larger than",  value=f_inf))
            delta_val =  float(st.text_input("Delta (max - min value) >",  value=0))
            fold_val =  float(st.text_input("Fold (max / min) >",  value=1))
            min_variance =  float(st.text_input("Minmum variance across samples > (e.g., 0.3)",  value=0))
            high_min_val =  float(st.text_input("All values of each gene are smaller than or equal",  value=p_inf))
            high_max_val =  float(st.text_input("Min value is smaller than or equal",  value=p_inf))
            top_n =  float(st.text_input("Top n in mean",  value=p_inf))

        st.markdown("######   ")
        calc_div = st.checkbox('Divided by (e.g., 1000, 1000,000)?', value = False)
        if calc_div:
            div_unit =  int(st.text_input("Unit: ",  value=1))
            df = df/div_unit
        calc_log = st.checkbox('Log transformation?')
        if calc_log:
            howlog = st.radio('Method', ['log2+1', 'log2', 'loge+1', 'loge','log10+1','log10', 'asinh'])
        else:
            howlog ='No'

        st.markdown("######   ")

        if min_val != f_inf:
            df = df[df.apply(min, axis=1) > min_val]


        if max_val != f_inf:
            df =  df[df.apply(max, axis=1) > max_val] # If this is min, it will be deleted if even one 0 exists.

        if delta_val > 1:
            df = df[df.apply(max, axis=1) > df.apply(min, axis=1) + delta_val]


        if fold_val > 1:
            df = df[df.apply(max, axis=1) > df.apply(min, axis=1) * fold_val]

        if min_variance > 0:
            df = df.loc[(df.index[(df.var(axis=1) > min_variance)]),:]

        if high_min_val != p_inf:
            df = df[df.apply(max, axis=1) <= high_min_val]

        if high_max_val != p_inf:
            df =  df[df.apply(min, axis=1) <= high_max_val] # If this is min, it will be deleted if even one 0 exists.

        if top_n != p_inf:
            top_ix = df.mean(axis = 1).sort_values(ascending=False).head(10).index
            new_index = [x for x in df.index.to_list() if x in top_ix]
            df = df.loc[new_index,:]

        df = df.astype('float')

        if howlog == 'log2+1':
            df = np.log2(df+1)
        elif howlog == 'log2':
            df = np.log2(df)
        elif howlog == 'loge+1':
            df = np.log1p(df)
        elif howlog == 'loge':
            df = np.log(df)
        elif howlog == 'log10+1':
            df = np.log10(df+1)
        elif howlog == 'log10':
            df = np.log10(df)
        elif howlog == 'asinh':
            df = np.arcsinh(df)



    st.markdown("##### Data standardization (Z-score transformation) ?")
    calc_z = st.checkbox('Z-score?', label_visibility =  'collapsed')
    if calc_z:
        center0_z= True
        df_z = df.copy()
        m = df_z.mean(1)
        s = df_z.std(1)
        df_z = df_z.sub(m, axis=0).div(s, axis = 0)
        df_z = np.round(df_z, decimals=10)
        df_z = df_z.loc[~(df_z==0).all(axis=1)] # Remove rows where all are 0
        df_z = df_z.dropna(how='any', axis=0) # Error handling

        # Check for columns with near-zero variance after Z-score transformation
        # These columns will appear white/neutral because values are near 0
        col_stds = df_z.std(axis=0)
        low_var_cols = col_stds[col_stds < 0.1].index.tolist()
        if low_var_cols:
            st.warning(f"⚠️ Z-score warning: The following columns have very low variance after transformation and will appear white in the heatmap: {low_var_cols}")
            st.info("💡 The original data for these columns is very close to the mean value for each row, so they become nearly 0 after Z-score transformation.")

        df = df_z


    st.markdown('---')


    st.markdown("##### Use subset of genes (or rows)?")
    subset_gene = st.checkbox('Use subset of genes (rows)?',label_visibility =  'collapsed')
    if subset_gene:
        st.markdown("##### Genes (comma, semicolon, space, CR separated):")
        genes = st.text_input("genes", label_visibility='collapsed')
        keep_all = st.checkbox('Do not remove duplicated genes?')
        gene_list = []
        if len(genes) > 0:
            genes = genes.replace("'", "")
            genes = genes.replace('"', "")
            gene_list = genes.split(' ')  # First separate by space
            gene_list = list(filter(lambda a: a != '', gene_list))  # Remove only spaces

            # Split by each delimiter while preserving order
            if ',' in genes:
                gene_list = sum([x.split(',') for x in gene_list], [])
            if ';' in genes:
                gene_list = sum([x.split(';') for x in gene_list], [])
            if '\t' in genes:
                gene_list = sum([x.split('\t') for x in gene_list], [])
            if '\n' in genes:
                gene_list = sum([x.split('\n') for x in gene_list], [])

            # Remove duplicates while preserving order
            seen = set()
            ordered_unique_genes = []
            for gene in gene_list:
                gene_lower = gene.lower()
                if gene_lower not in seen:
                    seen.add(gene_lower)
                    ordered_unique_genes.append(gene)

            if keep_all:
                ordered_unique_genes = gene_list

            # Match with dataframe index (preserve order)
            df_index_set = set(df.index.tolist())
            gene_subset = []
            for gene in ordered_unique_genes:
                if any(x.lower() == gene.lower() for x in df_index_set):
                    matching_gene = next(x for x in df_index_set if x.lower() == gene.lower())
                    gene_subset.append(matching_gene)

            df = df.loc[gene_subset, :]

            # Apply same subset to P-value matrix
            if use_pvalue_overlay_check and df_pvalues is not None:
                df_pvalues = df_pvalues.loc[gene_subset, :]
            if use_stouffer and df_nominal_pvalues is not None:
                df_nominal_pvalues = df_nominal_pvalues.loc[gene_subset, :]

    st.markdown('---')

    df = df.astype(float)
    st.markdown("##### Cleaned data:")
 #   df.iloc[:3,:]
    st.write(df.head())
    st.write('Data Dimension: '+str(df.shape))
    st.markdown('---')

    with st.sidebar:
        only_minmax = st.checkbox('**Two colors?**')
        if only_minmax:# Make it easy to change to two colors here
            v_center = None
            c_cmap = 'viridis'
            cmap = 'viridis'
        else:
            c_cmap = 'bwr'
            cmap = 'bwr'

        v_min = None
        v_max = None
        if calc_z:
            v_center = 0
        else:
            v_center = None


        annot  = False
        fmt = None
        linewidths= 0
        change_c = st.checkbox('**Change color map?**')
        if change_c:
            create_c = st.checkbox('Create custom color map?')
            if not create_c:
                category = st.selectbox("Select colormap category", options=list(colormap_categories.keys()))
                c_cmap = st.selectbox("Select colormap", options=colormap_categories[category])
                inv_c = st.checkbox('Inverse color map?')
                if inv_c:
                    c_cmap = c_cmap + "_r"
                show_colormap(c_cmap)

#            from matplotlib.colors import ListedColormap
            if create_c:
            #    two_c = st.checkbox('Make two-color map?')
                if only_minmax:
                    min_c = st.color_picker('Min color:', '#ffffff')
                    max_c = st.color_picker('Max color:', '#EE4B2B')
                    colors = [min_c, max_c]
                else:
                    min_c = st.color_picker('Min color:', '#0096FF')
                    center_c = st.color_picker('Center color:', '#ffffff')
                    max_c = st.color_picker('Max color:', '#EE4B2B')
                    colors = [min_c, center_c, max_c]
                c_cmap = create_custom_cmap(colors, name="custom")
                show_colormap(c_cmap)
            st.write(f"Selected colormap: {c_cmap}")
            cmap = c_cmap

            st.markdown('---')

        # Plot range settings (always visible)
        # Check if data has both positive and negative values
        has_positive = (df.values > 0).any()
        has_negative = (df.values < 0).any()
        center0_default = center0_z or (has_positive and has_negative)
        center0 = st.checkbox('Set center as 0?', value = center0_default)
        if center0:
            v_center = 0
        else:
            v_center = None
        minmax = st.checkbox('Change min/center/max?')
        if minmax:
            st.write('min: ' + str(df.to_numpy().min()) + '    mean: ' + str(df.to_numpy().mean()) + "    max: " + str(df.to_numpy().max()) )
            if only_minmax:
                v_min =  float(st.text_input("Min value",  value=df.to_numpy().min()))
                v_max =  float(st.text_input("Max value",  value=df.to_numpy().max()))
                v_center = None
            else:
                v_min =  float(st.text_input("Min value",  value=df.to_numpy().min()))
                v_center = float(st.text_input("Center",  value=df.to_numpy().mean()))
                v_max =  float(st.text_input("Max value",  value=df.to_numpy().max()))

        st.markdown('---')

        annot_on = st.checkbox('Show value on each cell?')
        if annot_on:
            annot = True
            annot_digit = st.text_input("Number of decimal palces",  value=0)
            annot_digit = int(annot_digit)
    #        if annot_digit == 0:
    #            fmt = "d" #"d" doesn't work well
    #        else:
            fmt = "." + str(annot_digit) + 'f'
    #        st.write(fmt)


        grid_on = st.checkbox('Show grid lines?')
        if grid_on:
            linewidths = st.number_input("Grid line width",  value=0.1, step=0.001, format="%.3f")
            line_col = st.selectbox("Grid line color", ('white', 'black','red','blue','gray','yellow'), index=0)
        else:
            line_col = 'white'

        # Boundary line options for clustering
        show_boundary = st.checkbox('Show cluster boundary lines?', value=False)
        if show_boundary:
            boundary_line = st.number_input('Boundary line width', value=0.5, step=0.1, format="%.1f")
            boundary_color = st.selectbox("Boundary line color", ('black', 'white', 'red','blue','gray','yellow'), index=0)
        else:
            boundary_line = 0.0
            boundary_color = 'black'

        st.markdown('---')


    # Currently cannot change plot size
        py_x_size = 8
        py_y_size = 8

        st.markdown('#### Plot size')
        py_x_size = float(st.text_input("Plot x size:", value = 8))
        py_y_size = float(st.text_input("Plot y size:", value = 8))
        st.markdown('#### Font size')
    #    sns_font_scale = float(st.text_input("Font scale:", value = 1))
        x_font_size = float(st.text_input("Sample name (column) font size:", value = 12))
        y_font_size = float(st.text_input("Gene name (row) font size:", value = 12))

        xticklabels= "auto"
        yticklabels= "auto"
        x_all = st.checkbox("Show all sample (column) names?", value = False)
        y_all = st.checkbox("Show all gene (row) names? !!!Do not check this for a large number of genes!!!", value = False)
        if x_all:
            xticklabels=1

        if y_all:
            yticklabels=1

        # P-value dot display options (display only in P-value overlay mode)
        if use_pvalue_overlay_check and df_pvalues is not None:
            st.markdown('---')
            st.markdown('#### P-value dot options')

            # Dot color and transparency
            col1, col2 = st.columns(2)
            with col1:
                dot_color = st.color_picker('Dot fill color', '#FFFFFF')
            with col2:
                dot_alpha = st.slider('Dot transparency', min_value=0.0, max_value=1.0, value=1.0, step=0.1)

            # Dot edge (outline) settings
            st.markdown('##### Dot edge (outline) settings')
            col3, col4 = st.columns(2)
            with col3:
                dot_edge_color = st.color_picker('Edge color', '#000000')
            with col4:
                dot_edge_width = st.slider('Edge width', min_value=0.0, max_value=3.0, value=0.1, step=0.1)

            # Dot size scaling method
            st.markdown('##### Dot size scaling method')
            scale_method = st.radio(
                'Select scaling method for dot size based on p-value:',
                options=['log10', 'linear'],
                index=0,
                help='log10: -log10(p) transformation (emphasizes small p-values)\n'
                     'linear: 1-p linear scaling (proportional to significance)'
            )

            # P-value threshold setting for visualization
            st.markdown('##### Threshold for visualization')
            p_threshold_visualization = st.slider(
                'Show dots only for p-values ≤ threshold',
                min_value=0.0,
                max_value=1.0,
                value=0.05,
                step=0.01,
                format="%.2f",
                help='Set the maximum P-value to display dots.\n'
                     'Default (0.05) shows dots only for p≤0.05.\n'
                     'Setting to 1.0 will display all P-values.'
            )

            # P-value threshold settings (for legend)
            st.markdown('##### P-value thresholds for legend')
            custom_thresholds = st.checkbox('Customize p-value thresholds', value=False)

            if custom_thresholds:
                p_threshold_1 = st.number_input('Threshold 1:', min_value=0.0, max_value=1.0, value=0.001, format="%.4f")
                p_threshold_2 = st.number_input('Threshold 2:', min_value=0.0, max_value=1.0, value=0.01, format="%.4f")
                p_threshold_3 = st.number_input('Threshold 3:', min_value=0.0, max_value=1.0, value=0.05, format="%.4f")
                p_thresholds = sorted([p_threshold_1, p_threshold_2, p_threshold_3])
            else:
                p_thresholds = [0.001, 0.01, 0.05]

            # Maximum -log10(p) value setting
            max_log_p = st.slider('Max -log10(p) for scaling', min_value=2, max_value=10, value=5,
                                 help='Controls the scaling range of dot size.\n'
                                      'Larger value = size changes gradually up to smaller P-values\n'
                                      'Smaller value = reaches maximum size earlier\n'
                                      'Example: value 5 → max size at p≤0.00001, value 3 → max size at p≤0.001')

            st.info(f"💡 Dot size changes gradually from p={max(p_thresholds):.3f} to p={10**(-max_log_p):.1e}\n\n"
                   f"**Explanation**: Dots with -log10(p-value) ≥ {max_log_p} (i.e., p≤{10**(-max_log_p):.1e}) will be maximum size")

            # Specify maximum dot size as percentage of cell height
            max_dot_size_percent = st.slider('Max dot size (% of cell height)',
                                            min_value=50, max_value=100, value=100, step=5,
                                            help='Specify maximum dot size as percentage of cell height.')

    st.markdown('##### Clustering:')
    if calc_z:
        # Exclude NMF during Z-score transformation as NMF cannot handle negative values
        clustering_type = st.radio("Clustering:", ('Nonclustering','Hierarchical','k-means','x-means', 'g-means'), label_visibility='collapsed')
    else:
        clustering_type = st.radio("Clustering:", ('Nonclustering','Hierarchical','k-means','x-means', 'g-means', 'NMF'), label_visibility='collapsed')
    if clustering_type == 'k-means':
        from kneed import KneeLocator
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        from sklearn.preprocessing import StandardScaler
        st.markdown('##### k-means options:')
        elbow = st.checkbox("Draw elbow plot and determine K automaticllay?", value = False)
        if "k" not in st.session_state:
            st.session_state.k = 3
        else:
            k_number = st.session_state.k
        if elbow:
            if st.button('Generate Elbow Plot'):
                with st.spinner('Generating elbow plot...'):
                    try:
                        sse = []
                        K = range(1, 11)
                        for k in K:
                            kmeans = KMeans(n_clusters=k, init='k-means++', n_init='auto', random_state=42)
                            kmeans.fit(df)
                            sse.append(kmeans.inertia_)
                        fig, ax = plt.subplots()
                        ax.plot(K, sse, 'bx-')
                        ax.set_xlabel('k')
                        ax.set_ylabel('Sum of squared distances')
                        ax.set_title('Elbow Method For Optimal k')
                        st.pyplot(fig)

                        kl = KneeLocator(K, sse, curve="convex", direction="decreasing")
                        if kl.elbow:
                            st.write(f"Optimal K suggested by elbow method: {kl.elbow}")
                            k_number = kl.elbow
                            st.session_state.k = k_number
                        else:
                            st.write("Could not determine optimal K automatically. Please select K manually.")
                    except Exception as e:
                        st.error(f"An error occurred while generating the elbow plot: {str(e)}")

        else:
            k_number = int(st.number_input("K clusters:", value = 3))
            st.session_state.k = k_number
        st.write("K number: " + str(k_number))

    elif clustering_type == 'NMF':
        st.markdown('##### NMF options:')
        st.markdown('NMF works best with non-negative data. Negative values will be set to 0.')

        with st.expander("ℹ️ About NMF (Non-negative Matrix Factorization)"):
            st.markdown("""
            **NMF** decomposes your data matrix into two non-negative matrices:
            - **W matrix**: Gene loadings on components (genes × components)
            - **H matrix**: Component loadings on samples (components × samples)

            **When to use NMF:**
            - Best for count data (RNA-seq, ChIP-seq peaks)
            - Finding parts-based representations
            - Identifying additive biological processes
            - Topic modeling in gene expression

            **Parameter Guide:**
            - **Number of components**: Start with 2-5 for exploratory analysis. More components = finer granularity
            - **Initialization methods**:
              - `nndsvd`: Best for sparse data (recommended for most genomics data)
              - `nndsvda`: NNDSVD with zeros filled with small random values
              - `nndsvdar`: NNDSVD with zeros filled with the data average
              - `random`: Random initialization (less reproducible)
            - **Max iterations**: Usually 200 is sufficient; increase if not converged
            """)

        if "n_components" not in st.session_state:
            st.session_state.n_components = 3
        n_components = int(st.number_input("Number of components:", value = st.session_state.n_components, min_value=2, max_value=min(df.shape),
                                         help="Number of latent components to extract. Each component represents a pattern in your data."))
        st.session_state.n_components = n_components

        nmf_init = st.selectbox("Initialization method:", ('nndsvd', 'nndsvda', 'nndsvdar', 'random'), index=0,
                               help="Algorithm for matrix initialization. 'nndsvd' is recommended for sparse genomics data.")

        nmf_max_iter = int(st.number_input("Max iterations:", value = 200, min_value=100, max_value=1000,
                                          help="Maximum number of iterations. Increase if the algorithm doesn't converge."))

        st.write(f"Number of components: {n_components}")

#    save_cluster = st.checkbox("Save cluster info?", value = False)
    st.markdown('---')
    y_c = False
    x_c = False
    if clustering_type == 'Hierarchical':
        import fastcluster
        y_c = st.checkbox("Cluster rows (Y axis)?", value = True)
        x_c = st.checkbox("Cluster colums (X axis)?", value = False)
        st.markdown('---')
        method_type = st.radio("Clustering method:", ('average','weighted', 'ward', 'median','single','centroid'))
        metric_type = st.radio("Clustering metric:", ('euclidean', 'seuclidean', 'sqeuclidean', 'minkowski', 'correlation', 'mahalanobis', 'cityblock', 'jaccard', 'jensenshannon'))
        st.markdown('---')


    save_type = st.radio("Save heatmap as: (Preparing PDF may take a time.)", ('png','pdf'))

    show_cor = st.checkbox('Show correlation coeficient matrix?')
    st.write("Clustering: " + clustering_type)
#    if "clustering" not in st.session_state:
#        st.session_state.clustering = False
    make_plot = st.button('Make plot')
    if make_plot:
        make_plot = False

        # If both are False, also change clustering_type # This prevents it from going to X-means
        if show_cor:
            correlation_coefficients = df.corr()
            fig_c, ax_c = plt.subplots() # This format is required to avoid errors
            ax_c = sns.heatmap(correlation_coefficients, vmax=1, vmin=-1, cmap='seismic', square=True,
                annot=False, xticklabels=1, yticklabels=1)
            st.pyplot(fig_c)
            fig_c.savefig(res_dir + "/corrlation." + save_type, format=save_type)

        df_file_name = file_name_head + '.Data4Heatmap.tsv'
        df.to_csv(res_dir + "/" + df_file_name,sep= '\t')

        # Stouffer's method: Calculate combined P-value and add Combined column
        if use_stouffer and df_nominal_pvalues is not None:
            st.markdown("---")
            st.markdown("##### Stouffer's method results")
            st.info("🔬 Combining p-values across conditions using Stouffer's method with BH correction...")

            # Execute Stouffer analysis
            combined_p, combined_fdr, mean_effects = stouffer_combine_pvalues(df_nominal_pvalues, df)

            # Display results
            stouffer_results = pd.DataFrame({
                'Gene': df.index,
                'Mean_NES': mean_effects,
                'Combined_P': combined_p,
                'Combined_FDR': combined_fdr
            })
            stouffer_results = stouffer_results.sort_values('Combined_P')

            st.write("Top 10 genes by combined p-value:")
            st.dataframe(stouffer_results.head(10))

            # Add Combined column to rightmost of df and df_pvalues
            df['Combined'] = mean_effects
            df_pvalues['Combined'] = combined_fdr

            # Save results to file
            stouffer_file_name = file_name_head + '.Stouffer_results.tsv'
            stouffer_results.to_csv(res_dir + "/" + stouffer_file_name, sep='\t', index=False)
            st.success(f"✓ Stouffer's method results saved to: {stouffer_file_name}")
            st.write(f"📊 'Combined' column added to heatmap (mean NES as color, combined FDR as dot size)")

        if clustering_type =="Nonclustering":
            g = sns.clustermap(df, center = v_center, cmap = cmap,
                    vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                    xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths, linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)

            # Add P-value dots
            if use_pvalue_overlay_check and df_pvalues is not None:
                df_pvalues_to_plot = df_pvalues.copy()

                add_significance_dots(
                    ax=g.ax_heatmap,
                    data_matrix=df,
                    pval_matrix=df_pvalues_to_plot,
                    max_log_p=max_log_p,
                    p_thresholds=p_thresholds,
                    dot_color=dot_color,
                    dot_alpha=dot_alpha,
                    max_dot_size_percent=max_dot_size_percent,
                    dot_edge_color=dot_edge_color,
                    dot_edge_width=dot_edge_width,
                    p_threshold_visualization=p_threshold_visualization,
                    scale_method=scale_method
                )

            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'Hierarchical':
            with st.spinner('Performing hierarchical clustering...'):
                # Exclude Combined column when using Stouffer for clustering
                has_combined = use_stouffer and 'Combined' in df.columns
                if has_combined:
                    df_for_clustering = df.drop(columns=['Combined'])
                    df_pvalues_for_clustering = df_pvalues.drop(columns=['Combined']) if df_pvalues is not None else None
                    st.info("ℹ️ 'Combined' column is excluded from clustering calculation (displayed on right)")
                else:
                    df_for_clustering = df
                    df_pvalues_for_clustering = df_pvalues

                # Execute clustering computation (cacheable)
                row_linkage, col_linkage, row_distances, col_distances  = perform_clustering_computation(df_for_clustering.values, method_type, metric_type)

                # If Combined column exists: pre-sort only columns, let clustermap handle rows
                if has_combined:
                    # Get column order (only if x_c is True)
                    if x_c:
                        g_temp = sns.clustermap(df_for_clustering, col_linkage=col_linkage, col_cluster=True, row_cluster=False)
                        col_order = g_temp.dendrogram_col.reordered_ind
                        plt.close(g_temp.fig)
                        # Sort columns and add Combined column to rightmost
                        df_with_combined = df_for_clustering.iloc[:, col_order].copy()
                        df_with_combined['Combined'] = df['Combined']
                        # Same for P-values
                        if df_pvalues is not None:
                            df_pvalues_with_combined = df_pvalues_for_clustering.iloc[:, col_order].copy()
                            df_pvalues_with_combined['Combined'] = df_pvalues['Combined']
                        else:
                            df_pvalues_with_combined = None
                    else:
                        # No column clustering: add Combined as is
                        df_with_combined = df_for_clustering.copy()
                        df_with_combined['Combined'] = df['Combined']
                        if df_pvalues is not None:
                            df_pvalues_with_combined = df_pvalues_for_clustering.copy()
                            df_pvalues_with_combined['Combined'] = df_pvalues['Combined']
                        else:
                            df_pvalues_with_combined = None

                    # Cluster only rows (use row_linkage), don't cluster columns
                    g = sns.clustermap(df_with_combined, row_linkage=row_linkage, col_linkage=None,
                                      center=v_center, cmap=cmap, vmin=v_min, vmax=v_max,
                                      row_cluster=y_c, col_cluster=False,
                                      xticklabels=xticklabels, yticklabels=yticklabels,
                                      annot=annot, fmt=fmt, linewidths=linewidths, linecolor=line_col,
                                      figsize=(py_x_size, py_y_size))
                    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=x_font_size)
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)

                    # Add P-value dots (corresponding to row clustering order)
                    if use_pvalue_overlay_check and df_pvalues_with_combined is not None:
                        df_pvalues_ordered = df_pvalues_with_combined.copy()
                        # Sort P-value data in row clustering order
                        if y_c:
                            row_order = g.dendrogram_row.reordered_ind
                            df_pvalues_ordered = df_pvalues_ordered.iloc[row_order, :]

                        add_significance_dots(
                            ax=g.ax_heatmap,
                            data_matrix=df_with_combined,
                            pval_matrix=df_pvalues_ordered,
                            max_log_p=max_log_p,
                            p_thresholds=p_thresholds,
                            dot_color=dot_color,
                            dot_alpha=dot_alpha,
                            max_dot_size_percent=max_dot_size_percent,
                            dot_edge_color=dot_edge_color,
                            dot_edge_width=dot_edge_width,
                            p_threshold_visualization=p_threshold_visualization,
                            scale_method=scale_method
                        )
                else:
                    # If no Combined column: conventional processing
                    g = plot_clustermap(
                        df_for_clustering, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c,
                        xticklabels, yticklabels, annot, fmt, linewidths, line_col,
                        py_x_size, py_y_size, x_font_size, y_font_size
                    )

                    # Add P-value dots (corresponding to clustering order)
                    if use_pvalue_overlay_check and df_pvalues_for_clustering is not None:
                        df_pvalues_ordered = df_pvalues_for_clustering.copy()

                        # Sort P-value data in clustering order
                        if y_c:  # Row clustering enabled
                            row_order = g.dendrogram_row.reordered_ind
                            df_pvalues_ordered = df_pvalues_ordered.iloc[row_order, :]

                        if x_c:  # Column clustering enabled
                            col_order = g.dendrogram_col.reordered_ind
                            df_pvalues_ordered = df_pvalues_ordered.iloc[:, col_order]

                        add_significance_dots(
                            ax=g.ax_heatmap,
                            data_matrix=df_for_clustering,
                            pval_matrix=df_pvalues_ordered,
                            max_log_p=max_log_p,
                            p_thresholds=p_thresholds,
                            dot_color=dot_color,
                            dot_alpha=dot_alpha,
                            max_dot_size_percent=max_dot_size_percent,
                            dot_edge_color=dot_edge_color,
                            dot_edge_width=dot_edge_width,
                            p_threshold_visualization=p_threshold_visualization,
                            scale_method=scale_method
                        )

            st.pyplot(g)

        elif clustering_type == 'k-means':
            with st.spinner('Performing k-means clustering...'):
                # Exclude Combined column when using Stouffer for clustering
                has_combined = use_stouffer and 'Combined' in df.columns
                if has_combined:
                    df_for_clustering = df.drop(columns=['Combined'])
                    df_pvalues_for_clustering = df_pvalues.drop(columns=['Combined']) if df_pvalues is not None else None
                    st.info("ℹ️ 'Combined' column is excluded from clustering calculation (displayed on right)")
                else:
                    df_for_clustering = df
                    df_pvalues_for_clustering = df_pvalues

                kmeans = KMeans(n_clusters=int(k_number), init ='k-means++', n_init='auto', random_state=42)
                clusters = kmeans.fit_predict(df_for_clustering)
                df2 = df_for_clustering.copy()
                df2["cluster"] = clusters
                df3 = pd.DataFrame(df2['cluster'], index = df2.index)
                st.write(df3.head(3))

                cluster_file_name = file_name_head + '.k-' + str(k_number) + '.tsv'

                df3.sort_values('cluster').to_csv(res_dir + "/" + cluster_file_name,sep= '\t')
                df2 = df2.sort_values('cluster')
                df2_sorted = df2.copy(deep=True)
                df2 = df2.drop('cluster', axis =1)

                # Add Combined column to rightmost
                if has_combined:
                    df2['Combined'] = df.loc[df2.index, 'Combined']

                g = sns.clustermap(df2, center = v_center, cmap = cmap,
                vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                         figsize = (py_x_size,py_y_size))
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
                if show_boundary:
                    ax = g.ax_heatmap
                    cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                    for boundary in cluster_boundaries[:-1]:
                        ax.axhline(y=boundary, color=boundary_color, linewidth=boundary_line)
                    # Create df_boundary
                    df_boundary = pd.DataFrame({
                        'position': cluster_boundaries
                    }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                    # Calculate size (first row is same as position, thereafter is difference)
                    df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                    st.write(df_boundary)
                    df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_kmeans' + str(k_number) + '_ClusterSize.tsv',sep= '\t')

                # Add P-value dots (corresponding to clustering order)
                if use_pvalue_overlay_check and df_pvalues_for_clustering is not None:
                    # Since df2 is sorted by clustering, also sort P-value data in same order
                    df_pvalues_sorted = df_pvalues_for_clustering.loc[df2.index]
                    # Add Combined column to rightmost
                    if has_combined and df_pvalues is not None:
                        df_pvalues_sorted['Combined'] = df_pvalues.loc[df_pvalues_sorted.index, 'Combined']

                    add_significance_dots(
                        ax=g.ax_heatmap,
                        data_matrix=df2,
                        pval_matrix=df_pvalues_sorted,
                        max_log_p=max_log_p,
                        p_thresholds=p_thresholds,
                        dot_color=dot_color,
                        dot_alpha=dot_alpha,
                        max_dot_size_percent=max_dot_size_percent,
                        dot_edge_color=dot_edge_color,
                        dot_edge_width=dot_edge_width,
                        p_threshold_visualization=p_threshold_visualization,
                        scale_method=scale_method
                    )

            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'x-means':
            st.write("Calculating x-means...")
            from pyclustering.cluster import xmeans

            # Exclude Combined column when using Stouffer for clustering
            if use_stouffer and 'Combined' in df.columns:
                df_for_clustering = df.drop(columns=['Combined'])
                df_pvalues_for_clustering = df_pvalues.drop(columns=['Combined']) if df_pvalues is not None else None
                st.info("ℹ️ 'Combined' column is excluded from clustering calculation")
            else:
                df_for_clustering = df
                df_pvalues_for_clustering = df_pvalues

            initial_centers = xmeans.kmeans_plusplus_initializer(df_for_clustering, 2).initialize()
            xm = xmeans.xmeans(df_for_clustering, initial_centers=initial_centers, )
            xm.process()
            clusters = xm.predict(df_for_clustering)
            df2 = df_for_clustering.copy()
            df2["cluster"] = clusters
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head + '.xmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            df2 = df2.sort_values('cluster')
            df2_sorted = df2.copy(deep=True)
            df2 = df2.drop('cluster', axis =1)

            # Add Combined column to rightmost
            has_combined = use_stouffer and 'Combined' in df.columns
            if has_combined:
                df2['Combined'] = df.loc[df2.index, 'Combined']

            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
            if boundary_line > 0:
                ax = g.ax_heatmap
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                for boundary in cluster_boundaries[:-1]:
                    ax.axhline(y=boundary, color=line_col, linewidth=boundary_line)
                # Create df_boundary
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                # Calculate size (first row is same as position, thereafter is difference)
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_xmeans_ClusterSize.tsv',sep= '\t')

            # Add P-value dots (corresponding to clustering order)
            if use_pvalue_overlay_check and df_pvalues_for_clustering is not None:
                df_pvalues_sorted = df_pvalues_for_clustering.loc[df2.index]
                # Add Combined column to rightmost
                if has_combined and df_pvalues is not None:
                    df_pvalues_sorted['Combined'] = df_pvalues.loc[df_pvalues_sorted.index, 'Combined']

                add_significance_dots(
                    ax=g.ax_heatmap,
                    data_matrix=df2,
                    pval_matrix=df_pvalues_sorted,
                    max_log_p=max_log_p,
                    p_thresholds=p_thresholds,
                    dot_color=dot_color,
                    dot_alpha=dot_alpha,
                    max_dot_size_percent=max_dot_size_percent,
                    dot_edge_color=dot_edge_color,
                    dot_edge_width=dot_edge_width,
                    p_threshold_visualization=p_threshold_visualization,
                    scale_method=scale_method
                )

            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'g-means':
            from pyclustering.cluster import gmeans

            # Exclude Combined column when using Stouffer for clustering
            if use_stouffer and 'Combined' in df.columns:
                df_for_clustering = df.drop(columns=['Combined'])
                df_pvalues_for_clustering = df_pvalues.drop(columns=['Combined']) if df_pvalues is not None else None
                st.info("ℹ️ 'Combined' column is excluded from clustering calculation")
            else:
                df_for_clustering = df
                df_pvalues_for_clustering = df_pvalues

            ar = df_for_clustering.to_numpy()
            with st.spinner('This takes a long time...'):
                initial_centers = gmeans.kmeans_plusplus_initializer(ar, 2).initialize()
                gm = gmeans.gmeans(ar, initial_centers=initial_centers, )
                gm.process()
                clusters = gm.predict(ar)
            df2 = df_for_clustering.copy()
            df2["cluster"] = clusters
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head +  '.gmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            df2 = df2.sort_values('cluster')
            df2_sorted = df2.copy(deep=True)
            df2 = df2.drop('cluster', axis =1)

            # Add Combined column to rightmost
            has_combined = use_stouffer and 'Combined' in df.columns
            if has_combined:
                df2['Combined'] = df.loc[df2.index, 'Combined']

            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt, linewidths= linewidths,linecolor=line_col,
                     figsize = (py_x_size,py_y_size))
            g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
            g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)
            if boundary_line > 0:
                ax = g.ax_heatmap
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                for boundary in cluster_boundaries[:-1]:
                    ax.axhline(y=boundary, color=line_col, linewidth=boundary_line)
                # Create df_boundary
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                # Calculate size (first row is same as position, thereafter is difference)
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_gmeans_ClusterSize.tsv',sep= '\t')

            # Add P-value dots (corresponding to clustering order)
            if use_pvalue_overlay_check and df_pvalues_for_clustering is not None:
                df_pvalues_sorted = df_pvalues_for_clustering.loc[df2.index]
                # Add Combined column to rightmost
                if has_combined and df_pvalues is not None:
                    df_pvalues_sorted['Combined'] = df_pvalues.loc[df_pvalues_sorted.index, 'Combined']

                add_significance_dots(
                    ax=g.ax_heatmap,
                    data_matrix=df2,
                    pval_matrix=df_pvalues_sorted,
                    max_log_p=max_log_p,
                    p_thresholds=p_thresholds,
                    dot_color=dot_color,
                    dot_alpha=dot_alpha,
                    max_dot_size_percent=max_dot_size_percent,
                    dot_edge_color=dot_edge_color,
                    dot_edge_width=dot_edge_width,
                    p_threshold_visualization=p_threshold_visualization,
                    scale_method=scale_method
                )

            st.pyplot(g)

            st.markdown('---')

        elif clustering_type == 'NMF':
            with st.spinner('Performing NMF decomposition...'):
                # Exclude Combined column when using Stouffer for clustering
                if use_stouffer and 'Combined' in df.columns:
                    df_for_clustering = df.drop(columns=['Combined'])
                    df_pvalues_for_clustering = df_pvalues.drop(columns=['Combined']) if df_pvalues is not None else None
                    st.info("ℹ️ 'Combined' column is excluded from clustering calculation")
                else:
                    df_for_clustering = df
                    df_pvalues_for_clustering = df_pvalues

                # Ensure non-negative data for NMF
                df_nmf = df_for_clustering.copy()
                if df_nmf.min().min() < 0:
                    st.warning("Negative values detected. Setting them to 0 for NMF.")
                    df_nmf = df_nmf.clip(lower=0)

                # Perform NMF
                nmf_model = NMF(n_components=n_components, init=nmf_init, max_iter=nmf_max_iter, random_state=42)
                W = nmf_model.fit_transform(df_nmf)
                H = nmf_model.components_

                # Assign each gene to the component with the highest weight
                cluster_assignments = np.argmax(W, axis=1)

                # Add cluster information to dataframe
                df2 = df_for_clustering.copy()
                df2["cluster"] = cluster_assignments
                df3 = pd.DataFrame(df2['cluster'], index = df2.index)

                st.write("NMF Component Assignments:")
                st.dataframe(df3.head(10))

                # Save cluster assignments
                cluster_file_name = file_name_head + f'.NMF_{n_components}components.tsv'
                df3.sort_values('cluster').to_csv(res_dir + "/" + cluster_file_name, sep='\t')

                # Save W matrix (gene loadings on components)
                W_df = pd.DataFrame(W, index=df_for_clustering.index, columns=[f'Component_{i}' for i in range(n_components)])
                W_df.to_csv(res_dir + "/" + file_name_head + f'.NMF_W_matrix_{n_components}comp.tsv', sep='\t')

                # Save H matrix (component loadings on samples)
                H_df = pd.DataFrame(H, index=[f'Component_{i}' for i in range(n_components)], columns=df_for_clustering.columns)
                H_df.to_csv(res_dir + "/" + file_name_head + f'.NMF_H_matrix_{n_components}comp.tsv', sep='\t')

                # Sort genes by cluster for visualization
                df2 = df2.sort_values('cluster')
                df2_sorted = df2.copy(deep=True)
                df2 = df2.drop('cluster', axis=1)

                # Add Combined column to rightmost
                has_combined = use_stouffer and 'Combined' in df.columns
                if has_combined:
                    df2['Combined'] = df.loc[df2.index, 'Combined']

                # Create heatmap
                g = sns.clustermap(df2, center = v_center, cmap = cmap,
                    vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                    xticklabels=xticklabels, yticklabels=yticklabels, annot = annot, fmt = fmt,
                    linewidths= linewidths, linecolor=line_col,
                    figsize = (py_x_size, py_y_size))
                g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize = x_font_size)
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize = y_font_size)

                # Add cluster boundaries
                if show_boundary:
                    ax = g.ax_heatmap
                    cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                    for boundary in cluster_boundaries[:-1]:
                        ax.axhline(y=boundary, color=boundary_color, linewidth=boundary_line)

                    # Create boundary dataframe
                    df_boundary = pd.DataFrame({
                        'position': cluster_boundaries
                    }, index=pd.Index(range(len(cluster_boundaries)), name='component'))
                    df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                    st.write("Component sizes:")
                    st.write(df_boundary)
                    df_boundary.to_csv(res_dir + "/" + file_name_head + f'_NMF_{n_components}_ComponentSize.tsv', sep='\t')

                # Add P-value dots (corresponding to clustering order)
                if use_pvalue_overlay_check and df_pvalues_for_clustering is not None:
                    df_pvalues_sorted = df_pvalues_for_clustering.loc[df2.index]
                    # Add Combined column to rightmost
                    if has_combined and df_pvalues is not None:
                        df_pvalues_sorted['Combined'] = df_pvalues.loc[df_pvalues_sorted.index, 'Combined']

                    add_significance_dots(
                        ax=g.ax_heatmap,
                        data_matrix=df2,
                        pval_matrix=df_pvalues_sorted,
                        max_log_p=max_log_p,
                        p_thresholds=p_thresholds,
                        dot_color=dot_color,
                        dot_alpha=dot_alpha,
                        max_dot_size_percent=max_dot_size_percent,
                        dot_edge_color=dot_edge_color,
                        dot_edge_width=dot_edge_width,
                        p_threshold_visualization=p_threshold_visualization,
                        scale_method=scale_method
                    )

                st.pyplot(g)

                # Display reconstruction error
                reconstruction_error = nmf_model.reconstruction_err_
                st.write(f"NMF Reconstruction Error: {reconstruction_error:.4f}")

                # Display W matrix as heatmap
                with st.expander("View Gene-Component Matrix (W)"):
                    st.markdown("**W matrix**: Shows how much each gene contributes to each component")
                    # Adjust figure size based on W matrix size
                    w_fig_height = min(max(len(df.index) * 0.15, 5), 20)
                    fig_w, ax_w = plt.subplots(figsize=(n_components + 2, w_fig_height))

                    # Option to show only top genes
                    show_top_genes = st.checkbox("Show only top genes per component", value=True, key="w_matrix_top")
                    if show_top_genes:
                        top_n_genes = st.slider("Number of top genes per component:", min_value=5, max_value=50, value=20, key="w_matrix_topn")
                        # Get top N genes for each component
                        top_gene_indices = set()
                        for comp in range(n_components):
                            top_indices = np.argsort(W[:, comp])[-top_n_genes:]
                            top_gene_indices.update(top_indices)
                        top_gene_indices = sorted(list(top_gene_indices))
                        W_display = W_df.iloc[top_gene_indices, :]
                        st.write(f"Showing top {top_n_genes} genes per component ({len(top_gene_indices)} unique genes total)")
                    else:
                        W_display = W_df
                        st.write(f"Showing all {len(W_df)} genes")

                    sns.heatmap(W_display, cmap='YlOrRd', annot=False, xticklabels=1,
                               yticklabels=1 if len(W_display) <= 50 else 'auto', ax=ax_w, cbar_kws={'label': 'Weight'})
                    ax_w.set_title(f"Gene loadings on NMF components (W matrix)")
                    ax_w.set_xlabel("Components")
                    ax_w.set_ylabel("Genes")
                    st.pyplot(fig_w)

                    # Display top genes per component in text
                    st.markdown("**Top genes per component:**")
                    for comp in range(n_components):
                        top_gene_idx = np.argsort(W[:, comp])[-10:][::-1]
                        top_genes = df.index[top_gene_idx].tolist()
                        top_weights = W[top_gene_idx, comp]
                        st.write(f"Component_{comp}: {', '.join([f'{gene} ({weight:.3f})' for gene, weight in zip(top_genes[:5], top_weights[:5])])}")

                # Display H matrix as heatmap
                with st.expander("View Component-Sample Matrix (H)"):
                    st.markdown("**H matrix**: Shows how much each component is expressed in each sample")
                    fig_h, ax_h = plt.subplots(figsize=(py_x_size, n_components))
                    sns.heatmap(H_df, cmap='YlGnBu', annot=False, xticklabels=1, yticklabels=1, ax=ax_h, cbar_kws={'label': 'Weight'})
                    ax_h.set_title("NMF component expression across samples (H matrix)")
                    ax_h.set_xlabel("Samples")
                    ax_h.set_ylabel("Components")
                    st.pyplot(fig_h)

            st.markdown('---')


        if howlog == "No":
            logmethod = ""
        else:
            logmethod = "_" + howlog
        if calc_z:
            logmethod = logmethod + '.Z'
        if save_type == 'pdf':
            if clustering_type == 'k-means':
                file_name = file_name_head + logmethod +  '.k-' + str(k_number) + '.heatmap.pdf'
            elif clustering_type == 'NMF':
                file_name = file_name_head + logmethod + f'.NMF-{n_components}.heatmap.pdf'
            else:
                file_name = file_name_head + logmethod + '.heatmap.pdf'
        else:
            if clustering_type == 'k-means':
                file_name = file_name_head + logmethod +  '.k-' + str(k_number) + '.heatmap.png'
            elif clustering_type == 'NMF':
                file_name = file_name_head + logmethod + f'.NMF-{n_components}.heatmap.png'
            else:
                file_name = file_name_head + logmethod + '.heatmap.png'
        try:
            with st.spinner('Generating figure files...'):
                # Save both PNG and PDF
                import io

                # Generate filenames
                if clustering_type == 'k-means':
                    base_filename = file_name_head + logmethod + '.k-' + str(k_number) + '.heatmap'
                elif clustering_type == 'NMF':
                    base_filename = file_name_head + logmethod + f'.NMF-{n_components}.heatmap'
                else:
                    base_filename = file_name_head + logmethod + '.heatmap'

                png_filename = base_filename + '.png'
                pdf_filename = base_filename + '.pdf'

                # Save PNG
                g.savefig(res_dir + "/" + png_filename, format='png', dpi=300, bbox_inches='tight')
                # Save PDF
                g.savefig(res_dir + "/" + pdf_filename, format='pdf', bbox_inches='tight')

                # Store in session state for download buttons
                st.session_state['heatmap_fig'] = g
                st.session_state['heatmap_base_filename'] = base_filename

        except Exception as e:
            st.error(f"Error saving figure: {str(e)}")
            pass
        else:
            st.success("✓ Heatmap generated successfully!")

            # Download buttons for PNG and PDF
            col1, col2, col3 = st.columns(3)

            with col1:
                # PNG download
                png_buffer = io.BytesIO()
                g.savefig(png_buffer, format='png', dpi=300, bbox_inches='tight')
                png_buffer.seek(0)
                st.download_button(
                    label="📥 Download PNG",
                    data=png_buffer,
                    file_name=png_filename,
                    mime="image/png"
                )

            with col2:
                # PDF download
                pdf_buffer = io.BytesIO()
                g.savefig(pdf_buffer, format='pdf', bbox_inches='tight')
                pdf_buffer.seek(0)
                st.download_button(
                    label="📥 Download PDF",
                    data=pdf_buffer,
                    file_name=pdf_filename,
                    mime="application/pdf"
                )

            with col3:
                # ZIP download (all results)
                shutil.make_archive(temp_dir + "/Heatmap", format='zip', root_dir=res_dir)
                with open(temp_dir + "/Heatmap.zip", "rb") as fp:
                    st.download_button(
                        label="📦 Download All (ZIP)",
                        data=fp,
                        file_name=file_name_head + logmethod + '_' + clustering_type + ".Heatmap.zip",
                        mime="application/zip"
                    )
