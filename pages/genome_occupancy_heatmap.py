import streamlit as st
import pandas as pd
import numpy as np
import warnings
import seaborn as sns
import matplotlib.pyplot as plt
import re
import shutil
import os
import sys
import matplotlib.colors as mcolors
import io
from helper_func import mk_temp_dir
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage as scipy_linkage, fcluster

# pyclustering uses numpy.warnings, which has been removed in newer numpy
if not hasattr(np, 'warnings'):
    np.warnings = warnings

# Font settings
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['pdf.fonttype'] = 42  # Use TrueType fonts
plt.rcParams['ps.fonttype'] = 42   # Use TrueType fonts


st.set_page_config(page_title="Genome Occupancy Heatmap", page_icon="🌡")


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

if 'filename_add' not in globals():
    filename_add = ""


def checkbox_container(data):
    cols = st.columns(5)
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
        return mcolors.LinearSegmentedColormap.from_list(name, colors, N=256)
    elif len(colors) == 3:
        cmap1 = mcolors.LinearSegmentedColormap.from_list("cmap1", colors[:2], N=128)
        cmap2 = mcolors.LinearSegmentedColormap.from_list("cmap2", colors[1:], N=128)
        newcolors = np.vstack((cmap1(np.linspace(0, 1, 128)),
                               cmap2(np.linspace(0, 1, 128))))
        return mcolors.ListedColormap(newcolors, name=name)
    else:
        raise ValueError("Only 2 or 3 colors are supported for custom colormap creation.")

# Classify Matplotlib colormaps by category
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
    plt.close(fig)


@st.cache_data
def perform_clustering_computation(data, method, metric):
    st.write(method)
    st.write(metric)
    row_distances = pdist(data, metric)
    row_linkage = scipy_linkage(row_distances, method)
    col_distances = pdist(data.T, metric)
    col_linkage = scipy_linkage(col_distances, method)
    return row_linkage, col_linkage, row_distances, col_distances

def plot_clustermap(df, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c, xticklabels, yticklabels, py_x_size, py_y_size, x_font_size, y_font_size):
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
        figsize=(py_x_size, py_y_size)
    )
    clean_heatmap_xticks(g, x_font_size)
    if hide_gene_names:
        g.ax_heatmap.set_yticks([])
        g.ax_heatmap.set_ylabel('')
    else:
        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
    return g

def draw_sample_boundaries(g, df, color='black', linewidth=0.5, gap=0.0):
    """Draw vertical lines between samples in multi-sample data.
    Uses sample boundaries detected by parse_sample_columns.
    If gap > 0, draws white bands at boundaries to create gaps."""
    samples = parse_sample_columns(df)
    if samples is None or len(samples) <= 1:
        return
    ax = g.ax_heatmap
    # Cumulative column count per sample -> boundary positions
    cumulative = 0
    boundaries = []
    for sample_name, cols_info in samples.items():
        cumulative += len(cols_info)
        boundaries.append(cumulative)
    # Do not draw a line after the last sample
    for b in boundaries[:-1]:
        if gap > 0:
            ax.axvspan(b - gap / 2, b + gap / 2, color='white', zorder=2, linewidth=0)
        ax.axvline(x=b, color=color, linewidth=linewidth, zorder=3)


def draw_cluster_boundaries(g, boundaries, color='black', linewidth=0.5, gap=0.0):
    """Draw horizontal lines at cluster boundaries. If gap > 0, add white band gaps."""
    ax = g.ax_heatmap
    for b in boundaries:
        if gap > 0:
            ax.axhspan(b - gap / 2, b + gap / 2, color='white', zorder=2, linewidth=0)
        ax.axhline(y=b, color=color, linewidth=linewidth, zorder=3)


def clean_heatmap_xticks(g, x_font_size):
    """Remove |sample_name and N_ prefix from heatmap x-axis labels, showing position numbers only"""
    _prefix_re = re.compile(r'^\d+_(.+)$')
    labels = g.ax_heatmap.get_xmajorticklabels()
    new_labels = []
    for l in labels:
        t = l.get_text()
        if '|' in t:
            new_labels.append(t.split('|')[0])
        else:
            m = _prefix_re.match(t)
            new_labels.append(m.group(1) if m else t)
    g.ax_heatmap.set_xticklabels(new_labels, fontsize=x_font_size)


def clean_column_name(col):
    if isinstance(col, (int, float)):
        if col.is_integer():
            return str(int(col))
        return f'{col:g}'
    elif isinstance(col, str):
        try:
            num = float(col)
            if num.is_integer():
                return str(int(num))
            return f'{num:g}'
        except ValueError:
            return col
    else:
        return str(col)

def filter_transform_data(df, nonzero,
                          do_manip, div_unit, howlog,
                          min_val, max_val, delta_val, fold_val, min_variance,
                          high_min_val, high_max_val, top_n,
                          calc_z, gene_subset_tuple):
    """Cached filtering, transformation, Z-score, and gene subset"""
    if nonzero:
        df = df.loc[~(df == 0).all(axis=1)]
    df = df.dropna(how='any', axis=0)
    df = df.astype(float)

    if do_manip:
        if div_unit != 1:
            df = df / div_unit

        f_inf = -float('inf')
        p_inf = float('inf')

        if min_val != f_inf:
            df = df[df.min(axis=1) > min_val]
        if max_val != f_inf:
            df = df[df.max(axis=1) > max_val]
        if delta_val > 1:
            df = df[df.max(axis=1) > df.min(axis=1) + delta_val]
        if fold_val > 1:
            df = df[df.max(axis=1) > df.min(axis=1) * fold_val]
        if min_variance > 0:
            df = df.loc[df.index[df.var(axis=1) > min_variance], :]
        if high_min_val != p_inf:
            df = df[df.max(axis=1) <= high_min_val]
        if high_max_val != p_inf:
            df = df[df.min(axis=1) <= high_max_val]
        if top_n != p_inf:
            top_ix = df.mean(axis=1).sort_values(ascending=False).head(int(top_n)).index
            new_index = [x for x in df.index.to_list() if x in top_ix]
            df = df.loc[new_index, :]

        df = df.astype('float')

        if howlog == 'log2+1': df = np.log2(df + 1)
        elif howlog == 'log2': df = np.log2(df)
        elif howlog == 'loge+1': df = np.log1p(df)
        elif howlog == 'loge': df = np.log(df)
        elif howlog == 'log10+1': df = np.log10(df + 1)
        elif howlog == 'log10': df = np.log10(df)
        elif howlog == 'asinh': df = np.arcsinh(df)

    center0_z = False
    low_var_cols = []
    if calc_z:
        center0_z = True
        m = df.mean(1)
        s = df.std(1)
        df = df.sub(m, axis=0).div(s, axis=0)
        df = np.round(df, decimals=10)
        df = df.loc[~(df == 0).all(axis=1)]
        df = df.dropna(how='any', axis=0)
        col_stds = df.std(axis=0)
        low_var_cols = col_stds[col_stds < 0.1].index.tolist()

    if gene_subset_tuple is not None:
        df_index_set = set(df.index.tolist())
        gene_subset = []
        for gene in gene_subset_tuple:
            matches = [x for x in df_index_set if x.lower() == gene.lower()]
            if matches:
                gene_subset.append(matches[0])
        if gene_subset:
            df = df.loc[gene_subset, :]

    df = df.astype(float)
    return df, center0_z, low_var_cols


def sort_dataframe(df, sort_regions, sort_using):
    """deepTools plotHeatmap-style row sorting"""
    if sort_regions == 'no':
        return df
    agg_funcs = {
        'mean': np.mean, 'median': np.median,
        'max': np.max, 'min': np.min, 'sum': np.sum,
    }
    func = agg_funcs[sort_using]
    # Exclude tracking columns from sort calculation
    _exclude = {'_bed_pos', 'cluster', '_hc_cluster'}
    _data_cols = [c for c in df.columns if c not in _exclude]
    sort_vals = df[_data_cols].apply(func, axis=1)
    ascending = (sort_regions == 'ascend')
    return df.loc[sort_vals.sort_values(ascending=ascending).index]


def sort_within_clusters(df, cluster_col, sort_regions, sort_using):
    """Group by cluster order and sort rows within each cluster.
    In deepTools plotHeatmap, rows are globally sorted then grouped by cluster,
    effectively sorting within each cluster.
    Even when sort_regions='no', rows are grouped by cluster order."""
    parts = []
    for c in sorted(df[cluster_col].unique()):
        cluster_df = df[df[cluster_col] == c]
        if sort_regions != 'no':
            cluster_df = sort_dataframe(cluster_df, sort_regions, sort_using)
        parts.append(cluster_df)
    return pd.concat(parts)


def plot_average_profile(df, ax, x_font_size=12):
    """Draw average profile (line plot) above the heatmap"""
    means = df.mean(axis=0)
    x = np.arange(len(means))
    ax.plot(x, means.values, color='blue', linewidth=1.5)
    ax.fill_between(x, means.values, alpha=0.15, color='blue')
    # Display standard error
    sems = df.sem(axis=0)
    ax.fill_between(x, (means - sems).values, (means + sems).values,
                    alpha=0.25, color='blue')
    ax.set_xlim(0, len(means) - 1)
    ax.set_ylabel('Mean signal')
    ax.set_xticks([])
    ax.tick_params(axis='y', labelsize=x_font_size * 0.8)


def parse_sample_columns(df):
    """Group (position, col_name) by sample from column names.
    Supported formats:
      1. 'position|sample_name' (Homer annotatePeaks partial output)
      2. 'N_position' (app duplicate column rename: 2_-3000, 3_-3000, ...)
         Sample 1 has no prefix (-3000), Sample 2+ has N_ prefix
      3. Numeric only -> single sample
    """
    samples = {}

    # --- Format 1: position|sample_name ---
    for col in df.columns:
        s = str(col)
        if '|' in s:
            parts = s.split('|', 1)
            try:
                pos = float(parts[0])
                sample = parts[1]
            except ValueError:
                continue
            samples.setdefault(sample, []).append((pos, col))
    if samples:
        for s in samples:
            samples[s].sort(key=lambda x: x[0])
        return samples

    # --- Format 2: N_position (duplicate column name renaming) ---
    # Pattern: "N_-3000" or "N_3000" (N >= 2) vs "-3000" (Sample 1)
    import re as _re
    _prefix_pat = _re.compile(r'^(\d+)_(.+)$')
    has_prefix = False
    for col in df.columns:
        if _prefix_pat.match(str(col)):
            has_prefix = True
            break

    if has_prefix:
        sample_1_label = 'Sample_1'
        for col in df.columns:
            s = str(col)
            m = _prefix_pat.match(s)
            if m:
                sample_idx = int(m.group(1))
                pos_str = m.group(2)
                try:
                    pos = float(pos_str)
                except ValueError:
                    continue
                sample_label = f'Sample_{sample_idx}'
                samples.setdefault(sample_label, []).append((pos, col))
            else:
                try:
                    pos = float(s)
                except ValueError:
                    continue
                samples.setdefault(sample_1_label, []).append((pos, col))
        if samples:
            for s in samples:
                samples[s].sort(key=lambda x: x[0])
            # Sort by sample number order
            samples = dict(sorted(samples.items(),
                                  key=lambda kv: int(kv[0].split('_')[1])))
            return samples

    # --- Format 3: Numeric only -> single sample ---
    positions = []
    for col in df.columns:
        try:
            pos = float(str(col))
            positions.append((pos, col))
        except ValueError:
            continue
    if positions:
        positions.sort(key=lambda x: x[0])
        samples['Signal'] = positions
        return samples

    return None


def plot_cluster_profiles(df, cluster_labels, n_clusters, x_font_size=12, figsize_width=8):
    """Draw average profiles per cluster, plotted by sample.
    Samples are overlaid as different colored lines on a single graph."""
    samples = parse_sample_columns(df)
    if samples is None:
        return None

    palette = sns.color_palette("tab10", len(samples))
    fig, axes = plt.subplots(n_clusters, 1,
                             figsize=(figsize_width, 2.5 * n_clusters),
                             squeeze=False)

    for cl in range(n_clusters):
        ax = axes[cl, 0]
        mask = (cluster_labels == cl)
        cluster_df = df[mask]
        n_regions = mask.sum()

        for idx, (sample_name, cols_info) in enumerate(samples.items()):
            positions = [c[0] for c in cols_info]
            col_names = [c[1] for c in cols_info]
            means = cluster_df[col_names].mean(axis=0).values
            sems = cluster_df[col_names].sem(axis=0).values
            ax.plot(positions, means, color=palette[idx], linewidth=1.5, label=sample_name)
            ax.fill_between(positions, means - sems, means + sems,
                            alpha=0.2, color=palette[idx])

        ax.set_title(f"Cluster {cl} (n={n_regions})", fontsize=x_font_size)
        ax.set_ylabel('Mean signal', fontsize=x_font_size * 0.8)
        ax.tick_params(axis='both', labelsize=x_font_size * 0.7)
        if cl == 0:
            ax.legend(fontsize=x_font_size * 0.7, loc='upper right')

    axes[-1, 0].set_xlabel('Position (bp)', fontsize=x_font_size * 0.8)
    fig.tight_layout()
    return fig


def export_clustered_bed(bed_df, bed_positions, cluster_labels, res_dir, file_name_head, method_name):
    """Split BED file by cluster based on clustering results"""
    bed_rows = bed_df.iloc[bed_positions].reset_index(drop=True).copy()
    bed_rows['_cluster'] = cluster_labels
    all_beds = {}
    for c in sorted(bed_rows['_cluster'].unique()):
        c_bed = bed_rows[bed_rows['_cluster'] == c].drop('_cluster', axis=1)
        fname = f"{file_name_head}_{method_name}_cluster{int(c)}.bed"
        fpath = os.path.join(res_dir, fname)
        c_bed.to_csv(fpath, sep='\t', index=False, header=False)
        all_beds[int(c)] = {
            'path': fpath, 'count': len(c_bed), 'name': fname,
            'data': c_bed.to_csv(sep='\t', index=False, header=False).encode('utf-8')
        }
    return all_beds


st.markdown("### Genome Occupancy Heatmap")
st.markdown("Heatmap dedicated to Homer `annotatePeaks.pl` Genome Occupancy output")
st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")

if "go_temp_dir" not in st.session_state:
    st.session_state.go_temp_dir = True
    temp_dir, res_dir = mk_temp_dir("GenomeOccupancyHeatmap")
    st.session_state.go_temp_dir = temp_dir
else:
    temp_dir = st.session_state.go_temp_dir
    temp_dir, res_dir = mk_temp_dir("GenomeOccupancyHeatmap", temp_dir)


use_upload = 'Yes'
if 'go_df' in st.session_state:
    st.write("Data in memory:")
    st.write(st.session_state.go_df.head())
    if st.session_state.go_df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.go_df
        input_file_type = 'tsv'
        file_name_head = st.session_state.go_uploaded_file_name
        # Homer compatibility
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list():
            df = df.set_index('Row_name')
            df.index.name = "Gene"


if use_upload == 'Yes':
    input_file_type = st.radio(
        "Data format:",
        ('tsv','csv', 'excel'))
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

        # Genome occupancy: add sequential numbers to Homer annotate output columns
        org_col = df.iloc[0,:].tolist()
        from collections import defaultdict

        # Homer -ghist detection: pattern of repeating numeric column names
        _data_cols = [clean_column_name(c) for c in org_col[1:]]
        _homer_ghist = False
        _homer_n_samples = 0
        _homer_n_per_sample = 0
        try:
            _pos_values = [float(c) for c in _data_cols]
            for i in range(1, len(_pos_values)):
                if _pos_values[i] < _pos_values[i - 1]:
                    _homer_n_per_sample = i
                    break
            if _homer_n_per_sample > 0 and len(_pos_values) % _homer_n_per_sample == 0:
                _homer_n_samples = len(_pos_values) // _homer_n_per_sample
                if _homer_n_samples >= 2:
                    _first = _pos_values[:_homer_n_per_sample]
                    if all(_pos_values[s * _homer_n_per_sample:(s + 1) * _homer_n_per_sample] == _first
                           for s in range(1, _homer_n_samples)):
                        _homer_ghist = True
        except (ValueError, TypeError):
            pass

        if _homer_ghist:
            # Convert column names to position|sample_name format
            _default_names = [f"Sample_{i + 1}" for i in range(_homer_n_samples)]
            st.session_state['go_homer_ghist'] = True
            st.session_state['go_homer_n_samples'] = _homer_n_samples
            if 'go_homer_sample_names' not in st.session_state:
                st.session_state['go_homer_sample_names'] = _default_names
            _snames = st.session_state['go_homer_sample_names']
            new_columns = [clean_column_name(org_col[0])]
            for s_idx in range(_homer_n_samples):
                for p_idx in range(_homer_n_per_sample):
                    pos = _data_cols[s_idx * _homer_n_per_sample + p_idx]
                    new_columns.append(f"{pos}|{_snames[s_idx]}")
            st.info(f"Homer -ghist pattern detected: {_homer_n_samples} samples × {_homer_n_per_sample} positions")
        else:
            st.session_state['go_homer_ghist'] = False
            column_counts = defaultdict(int)
            new_columns = [org_col[0]]
            for col in org_col[1:]:
                col = clean_column_name(col)
                column_counts[col] += 1
                if column_counts[col] > 1:
                    new_col_name = f"{column_counts[col]}_{col}"
                else:
                    new_col_name = col
                new_columns.append(new_col_name)

        df.columns = new_columns

        df = df.drop(0, axis = 0)
        content = df.columns.tolist()
        if isinstance(content[0], float) and np.isnan(content[0]):
            st.write("0,0isnan")
            content[0] = "Gene"
            df.columns = content

        Gene_column = content[0]
        if "Annotation/Divergence" in content:
            search_word = '([^\ \(]+).*'

            for i in range(1, len(content)):
                match = re.search(search_word, content[i])
                if match:
                    content[i] = match.group(1).replace(' ', '_')
            df.columns = content
            df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str)

            pattern = "([^|]*)"
            repatter = re.compile(pattern)
            f_annotation = lambda x: repatter.match(x).group(1)
            df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
            df = df.loc[:,'Annotation/Divergence':]
            content = df.columns.tolist()
            content[0] = 'Gene'
            df.columns = content
            Gene_column = "Gene"
            st.write("Converted Annotation/Divergence to gene symbols.")

        elif "Gene" in content:
            Gene_column = "Gene"
        else:
            Gene_column = st.selectbox('Select gene name column', content)

        df = df.set_index(Gene_column)
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        st.session_state.go_uploaded_file_name = file_name_head
        st.session_state.go_df = df

    else:
        sys.exit(1)

if df is not None:
    st.markdown('---')
    nonzero = st.checkbox('Remove all zero genes?', value=True)

    st.write(df.head(3))

    # Check and fix duplicate column names
    if len(df.columns) != len(set(df.columns)):
        duplicate_cols = [col for col in df.columns if df.columns.tolist().count(col) > 1]
        unique_duplicates = list(set(duplicate_cols))
        st.warning(f"⚠️ Duplicate column names detected: {unique_duplicates}")
        st.info("💡 Renaming duplicate columns to make them unique...")

        cols = pd.Series(df.columns)
        for dup in unique_duplicates:
            dup_indices = [i for i, x in enumerate(df.columns) if x == dup]
            for idx, pos in enumerate(dup_indices):
                if idx > 0:
                    cols[pos] = f"{dup}_{idx}"
        df.columns = cols
        st.write("Renamed columns:")
        st.write(df.columns.tolist())

    # --- BED file upload (optional) ---
    st.markdown("##### BED file (optional)")
    st.markdown("Upload the original peak BED file to export cluster-specific BED files.")
    _bed_uploader = st.file_uploader("BED file", type=['bed', 'txt', 'tsv'], key='go_bed_uploader')
    has_bed = False
    if _bed_uploader is not None:
        _bed_df = pd.read_csv(_bed_uploader, sep='\t', header=None, comment='#')
        if len(_bed_df) == len(df):
            st.session_state['go_bed_df'] = _bed_df
            has_bed = True
            st.success(f"BED file loaded: {len(_bed_df)} regions (matches data rows)")
        else:
            st.error(f"BED row count ({len(_bed_df)}) does not match data row count ({len(df)}).")
            if 'go_bed_df' in st.session_state:
                del st.session_state['go_bed_df']
    elif 'go_bed_df' in st.session_state:
        _bed_df = st.session_state['go_bed_df']
        if len(_bed_df) == len(df):
            has_bed = True
            st.info(f"BED file in memory: {len(_bed_df)} regions")
        else:
            del st.session_state['go_bed_df']

    # Save pre-filter index for BED position tracking
    if has_bed:
        st.session_state['go_pre_filter_index'] = df.index.tolist()

    st.markdown('---')
    st.markdown("##### Filter and transform data?")
    calc_z = False
    center0_z = False
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
        else:
            div_unit = 1
        calc_log = st.checkbox('Log transformation?')
        if calc_log:
            howlog = st.radio('Method', ['log2+1', 'log2', 'loge+1', 'loge','log10+1','log10', 'asinh'])
        else:
            howlog ='No'


    st.markdown("##### Data standardization (Z-score transformation) ?")
    calc_z = st.checkbox('Z-score?', label_visibility =  'collapsed')

    st.markdown('---')

    st.markdown("##### Use subset of genes (or rows)?")
    subset_gene = st.checkbox('Use subset of genes (rows)?',label_visibility =  'collapsed')
    _gene_subset_tuple = None
    if subset_gene:
        st.markdown("##### Genes (comma, semicolon, space, CR separated):")
        genes = st.text_input("genes", label_visibility='collapsed')
        keep_all = st.checkbox('Do not remove duplicated genes?')
        gene_list = []
        if len(genes) > 0:
            genes = genes.replace("'", "")
            genes = genes.replace('"', "")
            gene_list = genes.split(' ')
            gene_list = list(filter(lambda a: a != '', gene_list))

            if ',' in genes:
                gene_list = sum([x.split(',') for x in gene_list], [])
            if ';' in genes:
                gene_list = sum([x.split(';') for x in gene_list], [])
            if '\t' in genes:
                gene_list = sum([x.split('\t') for x in gene_list], [])
            if '\n' in genes:
                gene_list = sum([x.split('\n') for x in gene_list], [])

            seen = set()
            ordered_unique_genes = []
            for gene in gene_list:
                gene_lower = gene.lower()
                if gene_lower not in seen:
                    seen.add(gene_lower)
                    ordered_unique_genes.append(gene)

            if keep_all:
                ordered_unique_genes = gene_list

            _gene_subset_tuple = tuple(ordered_unique_genes)

    st.markdown('---')

    # === Cached preprocessing (session_state) ===
    if not Manip:
        div_unit = 1
        min_val = -float('inf')
        max_val = -float('inf')
        high_min_val = float('inf')
        high_max_val = float('inf')
        delta_val = 1
        fold_val = 1
        min_variance = 0
        top_n = float('inf')

    _filter_key = (
        df.shape, tuple(df.columns),
        nonzero,
        Manip, div_unit, howlog,
        min_val, max_val, delta_val, fold_val, min_variance,
        high_min_val, high_max_val, top_n,
        calc_z, _gene_subset_tuple,
    )
    if (st.session_state.get('_go_filter_key') == _filter_key
            and '_go_filtered_df' in st.session_state):
        df = st.session_state['_go_filtered_df']
        center0_z = st.session_state['_go_center0_z']
        _low_var_cols = st.session_state['_go_low_var_cols']
    else:
        df, center0_z, _low_var_cols = filter_transform_data(
            df, nonzero,
            Manip, div_unit, howlog,
            min_val, max_val, delta_val, fold_val, min_variance,
            high_min_val, high_max_val, top_n,
            calc_z, _gene_subset_tuple,
        )
        st.session_state['_go_filtered_df'] = df
        st.session_state['_go_center0_z'] = center0_z
        st.session_state['_go_low_var_cols'] = _low_var_cols
        st.session_state['_go_filter_key'] = _filter_key

    # --- BED position tracking after filtering ---
    _bed_positions = None
    if has_bed:
        _orig_idx = st.session_state.get('go_pre_filter_index', [])
        bed_df = st.session_state.get('go_bed_df')
        if bed_df is not None and len(bed_df) == len(_orig_idx):
            # Consume-from-left matching: map filtered row names to original positions
            _remaining = list(range(len(_orig_idx)))
            _bed_positions = []
            for name in df.index:
                matched = False
                for i, pos in enumerate(_remaining):
                    if _orig_idx[pos] == name:
                        _bed_positions.append(pos)
                        _remaining.pop(i)
                        matched = True
                        break
                if not matched:
                    _bed_positions = None
                    has_bed = False
                    st.warning("BED position tracking failed: index mismatch.")
                    break
            if _bed_positions is not None and len(_bed_positions) != len(df):
                _bed_positions = None
                has_bed = False
        else:
            has_bed = False

    st.markdown("##### Cleaned data:")
    st.write(df.head())
    st.write('Data Dimension: '+str(df.shape))

    # Homer -ghist: sample name input UI
    if st.session_state.get('go_homer_ghist'):
        _n_samp = st.session_state['go_homer_n_samples']
        _cur_names = st.session_state.get('go_homer_sample_names',
                                           [f"Sample_{i+1}" for i in range(_n_samp)])
        with st.expander(f"Sample names ({_n_samp} samples)", expanded=False):
            _names_text = st.text_area(
                "One sample name per line:",
                value='\n'.join(_cur_names),
                height=min(200, 30 * _n_samp),
                key="go_sample_names_input")
            _new_names = [n.strip() for n in _names_text.strip().split('\n') if n.strip()]
            if len(_new_names) == _n_samp and _new_names != _cur_names:
                # Update column names
                _rename_map = {}
                for col in df.columns:
                    s = str(col)
                    if '|' in s:
                        pos_part, old_name = s.split('|', 1)
                        if old_name in _cur_names:
                            new_name = _new_names[_cur_names.index(old_name)]
                            _rename_map[col] = f"{pos_part}|{new_name}"
                if _rename_map:
                    df = df.rename(columns=_rename_map)
                    st.session_state['go_homer_sample_names'] = _new_names
                    st.session_state.go_df = df
                    st.rerun()
            elif len(_new_names) != _n_samp and _names_text.strip():
                st.warning(f"Expected {_n_samp} names, got {len(_new_names)}")

    st.markdown('---')

    with st.sidebar:
        only_minmax = st.checkbox('**Two colors?**', value=True)
        if only_minmax:
            v_center = None
            c_cmap = 'Reds'
            cmap = 'Reds'
        else:
            c_cmap = 'bwr'
            cmap = 'bwr'

        v_min = None
        v_max = None
        if not only_minmax and calc_z:
            v_center = 0
        else:
            v_center = None

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

            if create_c:
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

        if not only_minmax:
            center0 = st.checkbox('Set center as 0?', value = center0_z)
            if center0:
                v_center = 0
            else:
                v_center = None
        # deepTools compatible zMin/zMax (percentile input)
        _mat_vals = df.to_numpy().flatten()
        _mat_vals = _mat_vals[~np.isnan(_mat_vals)]

        minmax = st.checkbox('**zMin / zMax**' if only_minmax else '**zMin / center / zMax**', value=True)
        if minmax:
            _zmin_pct = float(st.number_input("zMin percentile", value=1.0, min_value=0.0, max_value=100.0, step=0.5, key="go_zmin_pct"))
            _zmax_pct = float(st.number_input("zMax percentile", value=98.0, min_value=0.0, max_value=100.0, step=0.5, key="go_zmax_pct"))
            v_min = float(np.percentile(_mat_vals, _zmin_pct))
            v_max = float(np.percentile(_mat_vals, _zmax_pct))
            st.caption(f"zMin = {v_min:.4g} ({_zmin_pct}‰) / zMax = {v_max:.4g} ({_zmax_pct}‰)")
            st.caption(f"actual min: {float(_mat_vals.min()):.4g} / max: {float(_mat_vals.max()):.4g}")
            if only_minmax:
                v_center = None
            else:
                v_center = float(st.text_input("Center", value=f"{float(_mat_vals.mean()):.4g}"))

        st.markdown('---')

        show_boundary = st.checkbox('Show cluster boundary lines?', value=False)
        if show_boundary:
            boundary_line = st.number_input('Boundary line width', value=0.5, step=0.1, format="%.1f")
            boundary_color = st.selectbox("Boundary line color", ('black', 'white', 'red','blue','gray','yellow'), index=0)
            cluster_gap = st.number_input('Cluster gap width (rows)', value=0.0, min_value=0.0, step=0.5, format="%.1f", key="go_cl_gap",
                                           help="Draw white bands on both sides of cluster boundaries to create gaps (in row units). Gaps may not align correctly with dendrograms.")
        else:
            boundary_line = 0.0
            boundary_color = 'black'
            cluster_gap = 0.0

        # Inter-sample boundary lines (vertical)
        show_sample_boundary = st.checkbox('Show sample boundary lines?', value=False,
                                            help="Show vertical lines between samples in multi-sample data")
        if show_sample_boundary:
            sample_boundary_width = st.number_input('Sample boundary line width', value=0.5, step=0.1, format="%.1f", key="go_sb_width")
            sample_boundary_color = st.selectbox("Sample boundary line color",
                ('black', 'white', 'red', 'blue', 'gray', 'yellow', 'green', 'orange'),
                index=0, key="go_sb_color")
            sample_boundary_gap = st.number_input('Gap width (columns)', value=0.0, min_value=0.0, step=0.5, format="%.1f", key="go_sb_gap",
                                                   help="Draw white bands on both sides of boundary lines to create gaps (in column units)")
        else:
            sample_boundary_width = 0.0
            sample_boundary_color = 'black'
            sample_boundary_gap = 0.0

        st.markdown('---')

        py_x_size = 8
        py_y_size = 8

        st.markdown('#### Plot size')
        py_x_size = float(st.text_input("Plot x size:", value = 8))
        py_y_size = float(st.text_input("Plot y size:", value = 8))
        st.markdown('#### Font size')
        x_font_size = float(st.text_input("Sample name (column) font size:", value = 12))

        hide_gene_names = st.checkbox("Hide gene names and tick marks?", value=True,
                                       help="Hide gene names and Y-axis tick marks")
        if hide_gene_names:
            yticklabels = False
            y_font_size = 12
        else:
            y_font_size = float(st.text_input("Gene name (row) font size:", value = 12))
            yticklabels = "auto"
            y_all = st.checkbox("Show all gene (row) names? !!!Do not check this for a large number of genes!!!", value = False)
            if y_all:
                yticklabels = 1

        xticklabels = "auto"
        x_all = st.checkbox("Show all sample (column) names?", value = False)
        if x_all:
            xticklabels = 1

    st.markdown('##### Sort rows:')
    sort_regions = st.selectbox("Sort regions:",
        ('no', 'descend', 'ascend'),
        index=0, key="go_sort_regions",
        help="no: keep original order, descend: descending, ascend: ascending")
    if sort_regions != 'no':
        sort_using = st.selectbox("Sort using:",
            ('mean', 'median', 'max', 'min', 'sum'),
            index=0, key="go_sort_using",
            help="Method to aggregate values for each row")
    else:
        sort_using = 'mean'

    st.markdown('---')

    st.markdown('##### Average profile plot:')
    show_profile = st.checkbox("Show average profile plot above heatmap?", value=False, key="go_show_profile")
    if show_profile:
        profile_height = float(st.number_input("Profile plot height (cm):", value=4.0, min_value=1.0, max_value=20.0, step=0.5, key="go_profile_h"))
    st.markdown('---')

    st.markdown('##### Clustering:')
    clustering_type = st.radio("Clustering:", ('Nonclustering','Hierarchical','k-means','x-means', 'g-means (slow)'), label_visibility='collapsed')

    # --- Common clustering: z-score normalization option ---
    clustering_zscore = False
    if clustering_type != 'Nonclustering':
        if calc_z:
            clustering_zscore = True  # Already z-scored, so automatically treated as enabled
        else:
            _zs_default = False
            clustering_zscore = st.checkbox(
                "Row-wise z-score normalization for clustering",
                value=_zs_default, key="go_clustering_zscore",
                help="Apply z-score normalization only for clustering distance calculation. "
                     "Heatmap display values remain as original data. "
                     "Mitigates the problem of outliers being isolated due to differences in signal intensity.")
        st.markdown('---')

    show_cluster_profile = False
    if clustering_type != 'Nonclustering':
        show_cluster_profile = st.checkbox(
            "Show per-cluster average profile plot?",
            value=False, key="go_cluster_profile",
            help="After clustering, show per-cluster average profiles by sample")

    if clustering_type == 'k-means':
        from kneed import KneeLocator
        from sklearn.cluster import KMeans
        from sklearn.metrics import silhouette_score
        from sklearn.preprocessing import StandardScaler
        st.markdown('##### k-means options:')
        elbow = st.checkbox("Draw elbow plot and determine K automaticllay?", value = False)
        if "go_k" not in st.session_state:
            st.session_state.go_k = 3
        else:
            k_number = st.session_state.go_k
        if elbow:
            if st.button('Generate Elbow Plot'):
                with st.spinner('Generating elbow plot...'):
                    try:
                        # Calculate elbow plot with z-score data as well
                        if clustering_zscore and not calc_z:
                            from scipy.stats import zscore
                            _elbow_data = zscore(df.values, axis=1, nan_policy='omit')
                            _elbow_data = np.nan_to_num(_elbow_data, nan=0.0)
                        else:
                            _elbow_data = df.values
                        sse = []
                        K = range(1, 11)
                        for k in K:
                            kmeans = KMeans(n_clusters=k, init='k-means++', n_init='auto', random_state=42)
                            kmeans.fit(_elbow_data)
                            sse.append(kmeans.inertia_)
                        fig, ax = plt.subplots()
                        ax.plot(K, sse, 'bx-')
                        ax.set_xlabel('k')
                        ax.set_ylabel('Sum of squared distances')
                        ax.set_title('Elbow Method For Optimal k')
                        st.pyplot(fig)
                        plt.close(fig)

                        kl = KneeLocator(K, sse, curve="convex", direction="decreasing")
                        if kl.elbow:
                            st.write(f"Optimal K suggested by elbow method: {kl.elbow}")
                            k_number = kl.elbow
                            st.session_state.go_k = k_number
                        else:
                            st.write("Could not determine optimal K automatically. Please select K manually.")
                    except Exception as e:
                        st.error(f"An error occurred while generating the elbow plot: {str(e)}")
        else:
            k_number = int(st.number_input("K clusters:", value = 3))
            st.session_state.go_k = k_number
        st.write("K number: " + str(k_number))

    st.markdown('---')
    y_c = False
    x_c = False
    if clustering_type == 'Hierarchical':
        import fastcluster
        y_c = st.checkbox("Cluster rows (Y axis)?", value = True)
        x_c = st.checkbox("Cluster colums (X axis)?", value = False)
        st.markdown('---')
        method_type = st.radio("Clustering method:", ('ward', 'average','weighted', 'median','single','centroid'))
        if method_type == 'ward':
            metric_type = 'euclidean'
            st.info("Ward linkage requires euclidean distance (fixed)")
        else:
            metric_type = st.radio("Clustering metric:", ('euclidean', 'seuclidean', 'sqeuclidean', 'minkowski', 'correlation', 'mahalanobis', 'cityblock', 'jaccard', 'jensenshannon'))
        st.markdown('---')
        cut_dendrogram = st.checkbox("Cut dendrogram into clusters?", value=False)
        if cut_dendrogram:
            n_clusters_cut = int(st.number_input("Number of clusters:", value=3, min_value=2, max_value=50, key="go_hc_nclust"))
        st.markdown('---')


    save_type = st.radio("Save heatmap as: (Preparing PDF may take a time.)", ('png','pdf'))

    show_cor = st.checkbox('Show correlation coeficient matrix?')
    st.write("Clustering: " + clustering_type)
    make_plot = st.button('Make plot')
    _btn_clicked = make_plot
    if make_plot:

        if show_cor:
            correlation_coefficients = df.corr()
            fig_c, ax_c = plt.subplots()
            ax_c = sns.heatmap(correlation_coefficients, vmax=1, vmin=-1, cmap='seismic', square=True,
                annot=False, xticklabels=1, yticklabels=1)
            st.pyplot(fig_c)
            fig_c.savefig(res_dir + "/corrlation." + save_type, format=save_type)
            plt.close(fig_c)

        df_file_name = file_name_head + '.Data4Heatmap.tsv'
        df.to_csv(res_dir + "/" + df_file_name,sep= '\t')

        if clustering_type =="Nonclustering":
            # Nonclustering: global sort
            if has_bed and _bed_positions is not None:
                _df_tmp = df.copy()
                _df_tmp['_bed_pos'] = _bed_positions
                _df_tmp = sort_dataframe(_df_tmp, sort_regions, sort_using)
                _sorted_bed_pos = _df_tmp['_bed_pos'].values.tolist()
                df_sorted = _df_tmp.drop('_bed_pos', axis=1)
            else:
                df_sorted = sort_dataframe(df, sort_regions, sort_using)

            if show_profile:
                fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                plot_average_profile(df_sorted, ax_profile, x_font_size)
                st.pyplot(fig_profile)
                plt.close(fig_profile)

            g = sns.clustermap(df_sorted, center = v_center, cmap = cmap,
                    vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                    xticklabels=xticklabels, yticklabels=yticklabels,
                     figsize = (py_x_size,py_y_size))
            clean_heatmap_xticks(g, x_font_size)
            if hide_gene_names:
                g.ax_heatmap.set_yticks([])
                g.ax_heatmap.set_ylabel('')
            else:
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
            if show_sample_boundary:
                draw_sample_boundaries(g, df_sorted, sample_boundary_color, sample_boundary_width, sample_boundary_gap)

            st.pyplot(g)

            # BED export (nonclustering = single group)
            if has_bed and _bed_positions is not None:
                bed_df = st.session_state['go_bed_df']
                _labels = np.zeros(len(_sorted_bed_pos), dtype=int)
                st.session_state['go_bed_export'] = export_clustered_bed(
                    bed_df, _sorted_bed_pos, _labels, res_dir, file_name_head, 'sorted')

            st.markdown('---')

        elif clustering_type == 'Hierarchical':
            with st.spinner('Performing hierarchical clustering...'):
                # z-score normalization: applied only for distance calculation (display uses original data)
                # If calc_z is already done, z-score transform is already applied, so skip
                if clustering_zscore and not calc_z:
                    from scipy.stats import zscore
                    data_for_clustering = zscore(df.values, axis=1, nan_policy='omit')
                    data_for_clustering = np.nan_to_num(data_for_clustering, nan=0.0)
                    st.info("Row-wise z-score normalization applied for clustering distance calculation.")
                else:
                    data_for_clustering = df.values

                # Linkage cache: skip recalculation if same data, method, and metric
                import hashlib
                _hc_data_hash = hashlib.md5(data_for_clustering.tobytes()).hexdigest()
                _hc_cache_key = (_hc_data_hash, method_type, metric_type)
                _hc_prev = st.session_state.get('go_hc_cache_key')
                if _hc_prev == _hc_cache_key and 'go_hc_linkage' in st.session_state:
                    row_linkage, col_linkage, row_distances, col_distances = st.session_state['go_hc_linkage']
                    st.info("Linkage cache hit — skipped distance/linkage recalculation.")
                else:
                    row_linkage, col_linkage, row_distances, col_distances = perform_clustering_computation(data_for_clustering, method_type, metric_type)
                    st.session_state['go_hc_linkage'] = (row_linkage, col_linkage, row_distances, col_distances)
                    st.session_state['go_hc_cache_key'] = _hc_cache_key

                if cut_dendrogram and y_c:
                    # Assign clusters using fcluster
                    cluster_labels = fcluster(row_linkage, t=n_clusters_cut, criterion='maxclust')
                    df_clusters = pd.DataFrame({'cluster': cluster_labels}, index=df.index)

                    # Get row order from dendrogram
                    from scipy.cluster.hierarchy import dendrogram as scipy_dendrogram
                    dendro = scipy_dendrogram(row_linkage, no_plot=True)
                    row_order = dendro['leaves']

                    # Reorder rows by dendrogram order
                    df_ordered = df.iloc[row_order]
                    clusters_ordered = cluster_labels[row_order]

                    # Renumber clusters by their first appearance order in dendrogram
                    seen = {}
                    new_label = 0
                    remapped = np.zeros_like(clusters_ordered)
                    for i, c in enumerate(clusters_ordered):
                        if c not in seen:
                            seen[c] = new_label
                            new_label += 1
                        remapped[i] = seen[c]

                    # Save cluster information
                    df_cluster_export = df_clusters.copy()
                    df_cluster_export.to_csv(
                        res_dir + "/" + file_name_head + f'.hierarchical_k{n_clusters_cut}.tsv', sep='\t')

                    # Cluster sizes
                    cluster_sizes = pd.Series(remapped).value_counts().sort_index()
                    df_boundary_info = pd.DataFrame({
                        'cluster': cluster_sizes.index,
                        'size': cluster_sizes.values,
                        'cumulative': cluster_sizes.cumsum().values,
                    })
                    st.write(f"Hierarchical clustering: {n_clusters_cut} clusters")
                    st.write(df_boundary_info)

                    # Palette for cluster color sidebar
                    _hc_palette = sns.color_palette("tab10", n_clusters_cut)

                    if sort_regions != 'no':
                        # --- With sorting: sort_within_clusters + color sidebar ---
                        df_ordered_with_cluster = df_ordered.copy()
                        df_ordered_with_cluster['_hc_cluster'] = remapped

                        if has_bed and _bed_positions is not None:
                            _hc_bed_pos = [_bed_positions[i] for i in row_order]
                            df_ordered_with_cluster['_bed_pos'] = _hc_bed_pos

                        df_plot = sort_within_clusters(
                            df_ordered_with_cluster, '_hc_cluster', sort_regions, sort_using)
                        plot_clusters = df_plot['_hc_cluster'].values

                        # BED export
                        if has_bed and '_bed_pos' in df_plot.columns:
                            _hc_final_bed_pos = df_plot['_bed_pos'].values.tolist()
                            df_plot = df_plot.drop(['_hc_cluster', '_bed_pos'], axis=1)
                            bed_df = st.session_state['go_bed_df']
                            st.session_state['go_bed_export'] = export_clustered_bed(
                                bed_df, _hc_final_bed_pos, plot_clusters.astype(int),
                                res_dir, file_name_head, f'hc_k{n_clusters_cut}')
                        else:
                            df_plot = df_plot.drop('_hc_cluster', axis=1)

                        # Color sidebar
                        row_colors = pd.Series(
                            [_hc_palette[int(c)] for c in plot_clusters],
                            index=df_plot.index, name='Cluster')

                        if show_profile:
                            fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                            plot_average_profile(df_plot, ax_profile, x_font_size)
                            st.pyplot(fig_profile)
                            plt.close(fig_profile)

                        g = sns.clustermap(df_plot, center=v_center, cmap=cmap,
                                vmin=v_min, vmax=v_max, row_cluster=False, col_cluster=x_c,
                                col_linkage=col_linkage if x_c else None,
                                row_colors=row_colors,
                                xticklabels=xticklabels, yticklabels=yticklabels,
                                figsize=(py_x_size, py_y_size))

                        # Draw boundary lines
                        if show_boundary:
                            _hc_bounds = np.cumsum(np.bincount(plot_clusters.astype(int)))
                            draw_cluster_boundaries(g, _hc_bounds[:-1], boundary_color, boundary_line, cluster_gap)

                        # per-cluster profile (with sorting)
                        _hc_cluster_labels_for_profile = plot_clusters.astype(int)
                        _hc_df_for_profile = df_plot

                    else:
                        # --- Without sorting: dendrogram display + color sidebar ---
                        # Use row_linkage order as-is -> dendrogram display possible
                        # Color sidebar created in original data row order (clustermap will reorder)
                        row_colors = pd.Series(
                            [_hc_palette[int(seen[c])] for c in cluster_labels],
                            index=df.index, name='Cluster')

                        # BED export (dendrogram order)
                        if has_bed and _bed_positions is not None:
                            _hc_bed_pos = [_bed_positions[i] for i in row_order]
                            bed_df = st.session_state['go_bed_df']
                            st.session_state['go_bed_export'] = export_clustered_bed(
                                bed_df, _hc_bed_pos, remapped.astype(int),
                                res_dir, file_name_head, f'hc_k{n_clusters_cut}')

                        if show_profile:
                            fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                            plot_average_profile(df, ax_profile, x_font_size)
                            st.pyplot(fig_profile)
                            plt.close(fig_profile)

                        g = sns.clustermap(df, center=v_center, cmap=cmap,
                                vmin=v_min, vmax=v_max, row_cluster=True, col_cluster=x_c,
                                row_linkage=row_linkage,
                                col_linkage=col_linkage if x_c else None,
                                row_colors=row_colors,
                                xticklabels=xticklabels, yticklabels=yticklabels,
                                figsize=(py_x_size, py_y_size))

                        if show_boundary and cluster_gap > 0:
                            st.warning("Cluster gap with dendrogram: gap positions may not align with dendrogram branches. Enable Sort regions for accurate display.")

                        # per-cluster profile (without sorting): use data in dendrogram order
                        _hc_cluster_labels_for_profile = remapped.astype(int)
                        _hc_df_for_profile = df_ordered

                    clean_heatmap_xticks(g, x_font_size)
                    if hide_gene_names:
                        g.ax_heatmap.set_yticks([])
                        g.ax_heatmap.set_ylabel('')
                    else:
                        g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)

                else:
                    if show_profile:
                        fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                        plot_average_profile(df, ax_profile, x_font_size)
                        st.pyplot(fig_profile)
                        plt.close(fig_profile)

                    g = plot_clustermap(
                        df, row_linkage, col_linkage, v_center, cmap, v_min, v_max, y_c, x_c,
                        xticklabels, yticklabels,
                        py_x_size, py_y_size, x_font_size, y_font_size
                    )

            if show_sample_boundary:
                draw_sample_boundaries(g, df, sample_boundary_color, sample_boundary_width, sample_boundary_gap)
            st.pyplot(g)

            # Per-cluster average profile plot (Hierarchical)
            if cut_dendrogram and y_c and show_cluster_profile:
                _cp_fig = plot_cluster_profiles(
                    _hc_df_for_profile, _hc_cluster_labels_for_profile,
                    n_clusters_cut, x_font_size, py_x_size)
                if _cp_fig:
                    st.pyplot(_cp_fig)
                    plt.close(_cp_fig)

        elif clustering_type == 'k-means':
            with st.spinner('Performing k-means clustering...'):
                # z-score normalization: applied only for clustering
                if clustering_zscore and not calc_z:
                    from scipy.stats import zscore
                    _km_data = zscore(df.values, axis=1, nan_policy='omit')
                    _km_data = np.nan_to_num(_km_data, nan=0.0)
                    st.info("Row-wise z-score normalization applied for k-means.")
                else:
                    _km_data = df.values
                kmeans = KMeans(n_clusters=int(k_number), init ='k-means++', n_init='auto', random_state=42)
                clusters = kmeans.fit_predict(_km_data)
                df2 = df.copy()
                df2["cluster"] = clusters
                if has_bed and _bed_positions is not None:
                    df2['_bed_pos'] = _bed_positions
                df3 = pd.DataFrame(df2['cluster'], index = df2.index)
                st.write(df3.head(3))

                cluster_file_name = file_name_head + '.k-' + str(k_number) + '.tsv'

                df3.sort_values('cluster').to_csv(res_dir + "/" + cluster_file_name,sep= '\t')
                # Sort within clusters
                df2 = sort_within_clusters(df2, 'cluster', sort_regions, sort_using)
                df2_sorted = df2.copy(deep=True)

                # Save cluster labels before dropping
                _km_cluster_labels_for_profile = df2_sorted['cluster'].values.astype(int)

                # Extract BED positions and export clustered BED
                if has_bed and '_bed_pos' in df2.columns:
                    _km_bed_pos = df2['_bed_pos'].values.tolist()
                    _km_clusters = _km_cluster_labels_for_profile
                    df2 = df2.drop(['cluster', '_bed_pos'], axis=1)
                    bed_df = st.session_state['go_bed_df']
                    st.session_state['go_bed_export'] = export_clustered_bed(
                        bed_df, _km_bed_pos, _km_clusters,
                        res_dir, file_name_head, f'kmeans_k{k_number}')
                else:
                    df2 = df2.drop('cluster', axis =1)

                if show_profile:
                    fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                    plot_average_profile(df2, ax_profile, x_font_size)
                    st.pyplot(fig_profile)
                    plt.close(fig_profile)

                g = sns.clustermap(df2, center = v_center, cmap = cmap,
                vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
                xticklabels=xticklabels, yticklabels=yticklabels,
                         figsize = (py_x_size,py_y_size))
                clean_heatmap_xticks(g, x_font_size)
                if hide_gene_names:
                    g.ax_heatmap.set_yticks([])
                    g.ax_heatmap.set_ylabel('')
                else:
                    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
                if show_boundary:
                    cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                    draw_cluster_boundaries(g, cluster_boundaries[:-1], boundary_color, boundary_line, cluster_gap)
                    df_boundary = pd.DataFrame({
                        'position': cluster_boundaries
                    }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                    df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                    st.write(df_boundary)
                    df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_kmeans' + str(k_number) + '_ClusterSize.tsv',sep= '\t')

            if show_sample_boundary:
                draw_sample_boundaries(g, df2, sample_boundary_color, sample_boundary_width, sample_boundary_gap)
            st.pyplot(g)

            # Per-cluster average profile plot (k-means)
            if show_cluster_profile:
                _cp_fig = plot_cluster_profiles(
                    df2, _km_cluster_labels_for_profile,
                    int(k_number), x_font_size, py_x_size)
                if _cp_fig:
                    st.pyplot(_cp_fig)
                    plt.close(_cp_fig)

            st.markdown('---')

        elif clustering_type == 'x-means':
            st.write("Calculating x-means...")
            from pyclustering.cluster import xmeans

            # z-score normalization: applied only for clustering
            if clustering_zscore and not calc_z:
                from scipy.stats import zscore
                _xm_data = zscore(df.values, axis=1, nan_policy='omit')
                _xm_data = np.nan_to_num(_xm_data, nan=0.0)
                _xm_data = pd.DataFrame(_xm_data, index=df.index, columns=df.columns)
                st.info("Row-wise z-score normalization applied for x-means.")
            else:
                _xm_data = df
            initial_centers = xmeans.kmeans_plusplus_initializer(_xm_data, 2).initialize()
            xm = xmeans.xmeans(_xm_data, initial_centers=initial_centers, )
            xm.process()
            clusters = xm.predict(_xm_data)
            df2 = df.copy()
            df2["cluster"] = clusters
            if has_bed and _bed_positions is not None:
                df2['_bed_pos'] = _bed_positions
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head + '.xmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            # Sort within clusters
            df2 = sort_within_clusters(df2, 'cluster', sort_regions, sort_using)
            df2_sorted = df2.copy(deep=True)

            # Save cluster labels before dropping
            _xm_cluster_labels_for_profile = df2_sorted['cluster'].values.astype(int)
            _xm_n_clusters = len(np.unique(_xm_cluster_labels_for_profile))

            # Extract BED positions and export clustered BED
            if has_bed and '_bed_pos' in df2.columns:
                _xm_bed_pos = df2['_bed_pos'].values.tolist()
                _xm_clusters = _xm_cluster_labels_for_profile
                df2 = df2.drop(['cluster', '_bed_pos'], axis=1)
                bed_df = st.session_state['go_bed_df']
                st.session_state['go_bed_export'] = export_clustered_bed(
                    bed_df, _xm_bed_pos, _xm_clusters,
                    res_dir, file_name_head, 'xmeans')
            else:
                df2 = df2.drop('cluster', axis =1)

            if show_profile:
                fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                plot_average_profile(df2, ax_profile, x_font_size)
                st.pyplot(fig_profile)
                plt.close(fig_profile)

            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels,
                     figsize = (py_x_size,py_y_size))
            clean_heatmap_xticks(g, x_font_size)
            if hide_gene_names:
                g.ax_heatmap.set_yticks([])
                g.ax_heatmap.set_ylabel('')
            else:
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
            if show_boundary:
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                draw_cluster_boundaries(g, cluster_boundaries[:-1], boundary_color, boundary_line, cluster_gap)
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_xmeans_ClusterSize.tsv',sep= '\t')

            if show_sample_boundary:
                draw_sample_boundaries(g, df2, sample_boundary_color, sample_boundary_width, sample_boundary_gap)
            st.pyplot(g)

            # Per-cluster average profile plot (x-means)
            if show_cluster_profile:
                _cp_fig = plot_cluster_profiles(
                    df2, _xm_cluster_labels_for_profile,
                    _xm_n_clusters, x_font_size, py_x_size)
                if _cp_fig:
                    st.pyplot(_cp_fig)
                    plt.close(_cp_fig)

            st.markdown('---')

        elif clustering_type == 'g-means (slow)':
            from pyclustering.cluster import gmeans

            # z-score normalization: applied only for clustering
            if clustering_zscore and not calc_z:
                from scipy.stats import zscore
                ar = zscore(df.values, axis=1, nan_policy='omit')
                ar = np.nan_to_num(ar, nan=0.0)
                st.info("Row-wise z-score normalization applied for g-means.")
            else:
                ar = df.to_numpy()
            with st.spinner('This takes a long time...'):
                initial_centers = gmeans.kmeans_plusplus_initializer(ar, 2).initialize()
                gm = gmeans.gmeans(ar, initial_centers=initial_centers, )
                gm.process()
                clusters = gm.predict(ar)
            df2 = df.copy()
            df2["cluster"] = clusters
            if has_bed and _bed_positions is not None:
                df2['_bed_pos'] = _bed_positions
            df3 = pd.DataFrame(df2['cluster'], index = df2.index)
            st.dataframe(df3)

            cluster_file_name = file_name_head +  '.gmeans.tsv'

            df3.sort_values('cluster').to_csv(res_dir  + "/" + cluster_file_name,sep= '\t')
            # Sort within clusters
            df2 = sort_within_clusters(df2, 'cluster', sort_regions, sort_using)
            df2_sorted = df2.copy(deep=True)

            # Save cluster labels before dropping
            _gm_cluster_labels_for_profile = df2_sorted['cluster'].values.astype(int)
            _gm_n_clusters = len(np.unique(_gm_cluster_labels_for_profile))

            # Extract BED positions and export clustered BED
            if has_bed and '_bed_pos' in df2.columns:
                _gm_bed_pos = df2['_bed_pos'].values.tolist()
                _gm_clusters = _gm_cluster_labels_for_profile
                df2 = df2.drop(['cluster', '_bed_pos'], axis=1)
                bed_df = st.session_state['go_bed_df']
                st.session_state['go_bed_export'] = export_clustered_bed(
                    bed_df, _gm_bed_pos, _gm_clusters,
                    res_dir, file_name_head, 'gmeans')
            else:
                df2 = df2.drop('cluster', axis =1)

            if show_profile:
                fig_profile, ax_profile = plt.subplots(figsize=(py_x_size, profile_height * 0.39))
                plot_average_profile(df2, ax_profile, x_font_size)
                st.pyplot(fig_profile)
                plt.close(fig_profile)

            g = sns.clustermap(df2, center = v_center, cmap = cmap,
            vmin= v_min, vmax = v_max, row_cluster= False, col_cluster = False,
            xticklabels=xticklabels, yticklabels=yticklabels,
                     figsize = (py_x_size,py_y_size))
            clean_heatmap_xticks(g, x_font_size)
            if hide_gene_names:
                g.ax_heatmap.set_yticks([])
                g.ax_heatmap.set_ylabel('')
            else:
                g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=y_font_size)
            if show_boundary:
                cluster_boundaries = np.cumsum(df2_sorted.groupby("cluster").size())
                draw_cluster_boundaries(g, cluster_boundaries[:-1], boundary_color, boundary_line, cluster_gap)
                df_boundary = pd.DataFrame({
                    'position': cluster_boundaries
                }, index=pd.Index(range(len(cluster_boundaries)), name='cluster'))
                df_boundary['size'] = df_boundary['position'].diff().fillna(df_boundary['position'].iloc[0].astype(int))
                st.write(df_boundary)
                df_boundary.to_csv(res_dir  + "/" +  file_name_head + '_gmeans_ClusterSize.tsv',sep= '\t')

            if show_sample_boundary:
                draw_sample_boundaries(g, df2, sample_boundary_color, sample_boundary_width, sample_boundary_gap)
            st.pyplot(g)

            # Per-cluster average profile plot (g-means)
            if show_cluster_profile:
                _cp_fig = plot_cluster_profiles(
                    df2, _gm_cluster_labels_for_profile,
                    _gm_n_clusters, x_font_size, py_x_size)
                if _cp_fig:
                    st.pyplot(_cp_fig)
                    plt.close(_cp_fig)

            st.markdown('---')


        if howlog == "No":
            logmethod = ""
        else:
            logmethod = "_" + howlog
        if calc_z:
            logmethod = logmethod + '.Z'
        if clustering_type == 'k-means':
            base_filename = file_name_head + logmethod + '.k-' + str(k_number) + '.heatmap'
        elif clustering_type == 'Hierarchical' and cut_dendrogram and y_c:
            base_filename = file_name_head + logmethod + '.hc-' + str(n_clusters_cut) + '.heatmap'
        else:
            base_filename = file_name_head + logmethod + '.heatmap'

        png_filename = base_filename + '.png'
        pdf_filename = base_filename + '.pdf'

        try:
            with st.spinner('Generating figure...'):
                _png_buf = io.BytesIO()
                g.savefig(_png_buf, format='png', dpi=300, bbox_inches='tight')
                st.session_state['go_heatmap_png_bytes'] = _png_buf.getvalue()
                st.session_state['go_heatmap_png_filename'] = png_filename
                st.session_state['go_heatmap_pdf_filename'] = pdf_filename
                st.session_state['go_heatmap_base_filename'] = base_filename
                st.session_state['go_heatmap_fig'] = g
                st.session_state['go_heatmap_pdf_bytes'] = None

        except Exception as e:
            st.error(f"Error saving figure: {str(e)}")
            pass
        else:
            st.success("✓ Heatmap generated successfully!")

        plt.close('all')

    # Display cached heatmap
    if 'go_heatmap_png_bytes' in st.session_state and st.session_state['go_heatmap_png_bytes']:
        if not _btn_clicked:
            st.image(st.session_state['go_heatmap_png_bytes'])
            st.success("Previous plot result (click Make plot to update after changing options)")

        col1, col2, col3 = st.columns(3)
        with col1:
            st.download_button(
                label="📥 Download PNG",
                data=st.session_state['go_heatmap_png_bytes'],
                file_name=st.session_state.get('go_heatmap_png_filename', 'heatmap.png'),
                mime="image/png",
                key="go_cached_dl_png"
            )
        with col2:
            if st.session_state.get('go_heatmap_pdf_bytes'):
                st.download_button(
                    label="📥 Download PDF",
                    data=st.session_state['go_heatmap_pdf_bytes'],
                    file_name=st.session_state.get('go_heatmap_pdf_filename', 'heatmap.pdf'),
                    mime="application/pdf",
                    key="go_cached_dl_pdf"
                )
            else:
                if st.button("📄 Generate PDF", key="go_gen_pdf"):
                    _fig = st.session_state.get('go_heatmap_fig')
                    if _fig is not None:
                        with st.spinner("Generating PDF..."):
                            _pdf_buf = io.BytesIO()
                            _fig.savefig(_pdf_buf, format='pdf', bbox_inches='tight')
                            st.session_state['go_heatmap_pdf_bytes'] = _pdf_buf.getvalue()
                            plt.close(_fig)
                            del st.session_state['go_heatmap_fig']
                        st.rerun()
        with col3:
            if st.button("📦 Generate ZIP", key="go_gen_zip"):
                import zipfile as _zf
                _zip_buf = io.BytesIO()
                with _zf.ZipFile(_zip_buf, 'w', _zf.ZIP_DEFLATED) as zf:
                    _base = st.session_state.get('go_heatmap_base_filename', 'heatmap')
                    zf.writestr(_base + '.png', st.session_state['go_heatmap_png_bytes'])
                    if st.session_state.get('go_heatmap_pdf_bytes'):
                        zf.writestr(_base + '.pdf', st.session_state['go_heatmap_pdf_bytes'])
                    _tsv_path = os.path.join(res_dir, file_name_head + '.Data4Heatmap.tsv')
                    if os.path.exists(_tsv_path):
                        zf.write(_tsv_path, os.path.basename(_tsv_path))
                st.session_state['go_heatmap_zip_bytes'] = _zip_buf.getvalue()
                st.rerun()

            if st.session_state.get('go_heatmap_zip_bytes'):
                _base = st.session_state.get('go_heatmap_base_filename', 'heatmap')
                st.download_button(
                    label="📦 Download ZIP",
                    data=st.session_state['go_heatmap_zip_bytes'],
                    file_name=_base + '.zip',
                    mime="application/zip",
                    key="go_cached_dl_zip"
                )

        # --- Clustered BED file downloads ---
        if 'go_bed_export' in st.session_state and st.session_state['go_bed_export']:
            st.markdown("##### Clustered BED files")
            _bed_export = st.session_state['go_bed_export']
            _bed_cols = st.columns(min(len(_bed_export), 4))
            for i, (c, info) in enumerate(sorted(_bed_export.items())):
                with _bed_cols[i % len(_bed_cols)]:
                    st.download_button(
                        label=f"Cluster {c} ({info['count']} regions)",
                        data=info['data'],
                        file_name=info['name'],
                        mime="text/tab-separated-values",
                        key=f"go_bed_dl_c{c}"
                    )
            # ZIP all BED files together
            if len(_bed_export) > 1:
                import zipfile as _zf_bed
                _bed_zip_buf = io.BytesIO()
                with _zf_bed.ZipFile(_bed_zip_buf, 'w', _zf_bed.ZIP_DEFLATED) as zf:
                    for c, info in sorted(_bed_export.items()):
                        zf.writestr(info['name'], info['data'])
                st.download_button(
                    label="📦 Download all BED (ZIP)",
                    data=_bed_zip_buf.getvalue(),
                    file_name=f"{file_name_head}_clustered_beds.zip",
                    mime="application/zip",
                    key="go_bed_dl_zip"
                )
