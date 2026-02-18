import streamlit as st
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import re
import shutil
import os
import plotly.express as px
import sys
from helper_func import mk_temp_dir, remove_after_space, remove_sample_num, remove_underscore, select_underscore
import plotly.graph_objects as go
# Add imports for statistical analysis
from scipy import stats
import itertools
from statsmodels.stats.multitest import multipletests
import math
from streamlit_sortables import sort_items

st.set_page_config(page_title="Box_Violin_plot", page_icon="◫")

st.sidebar.title("Options")


# Function to set plot style
def set_plot_style():
    plt.style.use('default')  # Reset to default style
    plt.rcParams['figure.facecolor'] = 'white'  # Set figure background to white
    plt.rcParams['axes.facecolor'] = 'white'    # Set plot area background to white
    sns.set_style("white")                      # Set seaborn style to white background

# Call before generating graphs
set_plot_style()

# Extended statistical test function: Added Student's t-test
def perform_statistical_test(data1, data2, test_type):
    """
    Perform statistical test between two groups

    Args:
        data1, data2: Data for the two groups
        test_type: Type of test ('t', 'student_t', 'u', or 'wilcoxon')

    Returns:
        p_value: p-value from the test
    """
    # Convert to numeric and remove NaN values
    data1 = pd.to_numeric(data1, errors='coerce').dropna()
    data2 = pd.to_numeric(data2, errors='coerce').dropna()

    if len(data1) == 0 or len(data2) == 0:
        return float('nan')

    try:
        if test_type == 't':
            # Perform t-test (Welch's t-test, does not assume equal variances)
            _, p_value = stats.ttest_ind(data1, data2, equal_var=False)
        elif test_type == 'student_t':
            # Perform Student's t-test (assumes equal variances)
            _, p_value = stats.ttest_ind(data1, data2, equal_var=True)
        elif test_type == 'u':
            # Perform Mann-Whitney U test
            _, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        elif test_type == 'wilcoxon':
            # Check if lengths are equal for paired test
            if len(data1) == len(data2):
                # Perform Wilcoxon signed-rank test
                _, p_value = stats.wilcoxon(data1, data2)
            else:
                # Fall back to Mann-Whitney U test for unpaired data
                _, p_value = stats.mannwhitneyu(data1, data2, alternative='two-sided')
        else:
            return float('nan')

        return p_value
    except:
        return float('nan')

# Multiple comparison correction function extended: Added Bonferroni, Holm-Bonferroni, and Benjamini-Yekutieli methods
def perform_multiple_tests(df, gene, test_type, groups, adjust_method=None):
    """
    Perform all pairwise statistical tests for a gene across groups

    Args:
        df: DataFrame with the data
        gene: Gene name to test
        test_type: Type of test ('t', 'student_t', 'u', or 'wilcoxon')
        groups: List of group names
        adjust_method: Method for p-value adjustment
                      (None, 'bh', 'bonferroni', 'holm', 'by', 'tukey')

    Returns:
        results: Dictionary with group pairs as keys and p-values as values
    """
    results = {}
    p_values = []
    group_pairs = []

    # Generate all possible pairs of groups
    for group1, group2 in itertools.combinations(groups, 2):
        group1_data = df.T[df.T['Group'] == group1][gene]
        group2_data = df.T[df.T['Group'] == group2][gene]

        p_value = perform_statistical_test(group1_data, group2_data, test_type)
        group_pairs.append((group1, group2))
        p_values.append(p_value)

    # Apply multiple comparison correction if specified
    if adjust_method == 'bh' and len(p_values) > 1:
        # Benjamini-Hochberg procedure
        try:
            adjusted_p = multipletests(p_values, method='fdr_bh')[1]
        except:
            adjusted_p = p_values
    elif adjust_method == 'bonferroni' and len(p_values) > 1:
        # Bonferroni correction
        try:
            adjusted_p = multipletests(p_values, method='bonferroni')[1]
        except:
            adjusted_p = [min(1.0, p * len(p_values)) for p in p_values]
    elif adjust_method == 'holm' and len(p_values) > 1:
        # Holm-Bonferroni method
        try:
            adjusted_p = multipletests(p_values, method='holm')[1]
        except:
            adjusted_p = p_values
    elif adjust_method == 'by' and len(p_values) > 1:
        # Benjamini-Yekutieli method
        try:
            adjusted_p = multipletests(p_values, method='fdr_by')[1]
        except:
            adjusted_p = p_values
    elif adjust_method == 'tukey' and len(p_values) > 1:
        # Tukey-Kramer method
        try:
            # Use a simple Bonferroni-like correction as approximation
            adjusted_p = [min(1.0, p * len(p_values)) for p in p_values]
        except:
            adjusted_p = p_values
    else:
        adjusted_p = p_values

    # Store results
    for i, (group1, group2) in enumerate(group_pairs):
        results[(group1, group2)] = adjusted_p[i]

    return results

def logit_transform(df, gene):
    """
    Apply logit transformation to ratio data

    Args:
        df: DataFrame with the data
        gene: Gene name to transform

    Returns:
        transformed_data: Transformed data series
    """
    # Convert to numeric and handle NaN
    data = pd.to_numeric(df[gene], errors='coerce')

    # Apply small offset to 0 and 1 values to avoid infinity
    epsilon = 1e-6
    data = data.apply(lambda x: epsilon if x <= 0 else (1-epsilon if x >= 1 else x))

    # Apply logit transformation: log(x/(1-x))
    transformed = np.log(data / (1 - data))

    return transformed

# Further improved function to add statistical annotations
def add_stat_annotations(fig, test_results, y_pos, groups, bar_height=15, x_offset=0.1, p_value_font_size=10, show_all_p=False):
    """
    Add statistical annotation bars and p-values to plot

    Args:
        fig: Plotly figure object
        test_results: Dictionary with group pairs as keys and p-values as values
        y_pos: Y position for the annotation
        groups: List of groups in the order they appear on the X-axis
        bar_height: Height of the annotation bar
        x_offset: Horizontal offset for the annotation bar
        p_value_font_size: Font size for p-value text
        show_all_p: Whether to show all p-values, including non-significant ones

    Returns:
        fig: Updated figure with annotations
        next_y_pos: Next Y position for annotations
    """
    max_y_pos = y_pos

    # Use the provided groups list to maintain X-axis order
    groups_list = list(groups)

    # Calculate offsets to shift bar height for each group pair
    group_pairs = []
    for (group1, group2), p_value in test_results.items():
        if not np.isnan(p_value):
            if p_value < 0.05 or show_all_p:
                # Store group pair along with index
                x1_pos = groups_list.index(group1)
                x2_pos = groups_list.index(group2)
                span = abs(x2_pos - x1_pos)
                group_pairs.append((group1, group2, span))

    # Sort by span (place shorter spans below)
    group_pairs.sort(key=lambda x: x[2], reverse=False)

    # Calculate offset according to number of bars (ensure sufficient spacing between bars)
    gap = bar_height * 1.2  # Ensure sufficient spacing between bars
    bar_offsets = {}
    current_offset = 0
    for i, (group1, group2, _) in enumerate(group_pairs):
        bar_offsets[(group1, group2)] = current_offset
        current_offset += gap

    # Add significance bars and p-value annotations
    for (group1, group2), p_value in test_results.items():
        if not np.isnan(p_value):
            # Skip non-significant cases (depending on option)
            if p_value >= 0.05 and not show_all_p:
                continue

            # Convert group names to numeric positions
            x1_pos = groups_list.index(group1)
            x2_pos = groups_list.index(group2)

            # Apply position offset for this bar
            current_offset = bar_offsets.get((group1, group2), 0)
            current_y_pos = y_pos + current_offset

            # Add horizontal bar
            fig.add_shape(
                type="line",
                x0=x1_pos, x1=x1_pos,
                y0=current_y_pos, y1=current_y_pos + bar_height,
                line=dict(color="black", width=1.5),
                xref="x", yref="y"
            )

            fig.add_shape(
                type="line",
                x0=x1_pos, x1=x2_pos,
                y0=current_y_pos + bar_height, y1=current_y_pos + bar_height,
                line=dict(color="black", width=1.5),
                xref="x", yref="y"
            )

            fig.add_shape(
                type="line",
                x0=x2_pos, x1=x2_pos,
                y0=current_y_pos, y1=current_y_pos + bar_height,
                line=dict(color="black", width=1.5),
                xref="x", yref="y"
            )

            # Display p-value as numbers only
            p_text = f"{p_value:.4g}"
            fig.add_annotation(
                x=(x1_pos + x2_pos) / 2,
                y=current_y_pos + bar_height,  # Display above bar
                text=p_text,
                showarrow=False,
                yshift=0,  # Slight upward shift
                font=dict(size=p_value_font_size),
                xanchor="center",
                yanchor="bottom"  # Changed to bottom alignment
            )

            max_y_pos = max(max_y_pos, current_y_pos + bar_height + 20)

    return fig, max_y_pos


def create_plot_with_points(df_show, gene, plot_type='box', color_split=False, points='all', c_choice=None,
                          py_x_size=600, py_y_size=400,
                          legend_font_size=12, x_font_size=12, y_font_size=12, gene_font_size=16,
                          jitter=0.3, y_min=None, y_max=None, show_box=True, marker_size=6,
                          stat_test=None, adjust_method=None, apply_logit=False, p_value_font_size=10,
                          show_all_p=False, stat_bar_height=15):

    # Apply logit transformation if specified
    plot_data = df_show.copy()
    if apply_logit and gene in plot_data.index:
        for col in plot_data.columns:
            if pd.to_numeric(plot_data.loc[gene, col], errors='coerce') is not None:
                try:
                    # Get column data and transform it
                    values = pd.to_numeric(plot_data.T[gene], errors='coerce')
                    # Apply logit transform only if values are between 0 and 1
                    if (values >= 0).all() and (values <= 1).all():
                        plot_data.loc[gene, :] = logit_transform(plot_data.T, gene)
                except:
                    pass

    gene_title = f"<i>{gene}</i>"

    if color_split:
        fig = go.Figure()

        # Use sorted order for traces if available, otherwise use data order
        color_group_order = list(df_show.T['Color_group'].unique())
        if st.session_state.get('sorted_groups') is not None:
            color_group_order = st.session_state.sorted_groups

        for idx, color_group in enumerate(color_group_order):
            group_data = df_show.T[df_show.T['Color_group'] == color_group]

            if plot_type == 'box':
                fig.add_trace(go.Box(
                    y=group_data[gene],
                    x=group_data['Condition'],
                    name=color_group,
                    marker_color=c_choice[idx % len(c_choice)],
                    boxpoints=points,
                    jitter=jitter,
                    pointpos=0,
                    marker_size=marker_size,
                    showlegend=True
                ))
            else:  # violin
                fig.add_trace(go.Violin(
                    y=group_data[gene],
                    x=group_data['Condition'],
                    name=color_group,
                    marker_color=c_choice[idx % len(c_choice)],
                    points=points,
                    jitter=jitter,
                    pointpos=0,
                    marker_size=marker_size,
                    showlegend=True,
                    box_visible=show_box,
                    meanline_visible=True
                ))
    else:
        fig = go.Figure()

        # Use sorted order for traces if available, otherwise use data order
        group_order = list(df_show.T['Group'].unique())
        if st.session_state.get('sorted_groups') is not None:
            group_order = st.session_state.sorted_groups

        for idx, group in enumerate(group_order):
            group_data = df_show.T[df_show.T['Group'] == group]

            if plot_type == 'box':
                fig.add_trace(go.Box(
                    y=group_data[gene],
                    name=group,
                    marker_color=c_choice[idx % len(c_choice)],
                    boxpoints=points,
                    jitter=jitter,
                    pointpos=0,
                    marker_size=marker_size,
                    showlegend=True
                ))
            else:  # violin
                fig.add_trace(go.Violin(
                    y=group_data[gene],
                    name=group,
                    marker_color=c_choice[idx % len(c_choice)],
                    points=points,
                    jitter=jitter,
                    pointpos=0,
                    marker_size=marker_size,
                    showlegend=True,
                    box_visible=show_box,
                    meanline_visible=True
                ))

    fig.update_layout(
        plot_bgcolor='rgb(255, 255, 255)',
        paper_bgcolor='rgb(255, 255, 255)',
        width=py_x_size,
        height=py_y_size,
        legend_font=dict(size=legend_font_size, family='Arial'),
        showlegend=True,
        font=dict(family='Arial')
    )

    # Set Y-axis range (only if specified)
    y_axis_range = {}
    if y_min is not None:
        y_axis_range['min'] = y_min
    if y_max is not None:
        y_axis_range['max'] = y_max

    # Set X-axis category order
    if not color_split:
        category_order = list(df_show.T['Group'].unique())
        # Use manually sorted order if available
        if st.session_state.get('sorted_groups') is not None:
            category_order = st.session_state.sorted_groups

        fig.update_xaxes(tickfont_size=x_font_size, title='', tickfont_family='Arial',
                        categoryorder='array', categoryarray=category_order)
    else:
        category_order = list(df_show.T['Color_group'].unique())
        # Use manually sorted order if available
        if st.session_state.get('sorted_groups') is not None:
            category_order = st.session_state.sorted_groups

        fig.update_xaxes(tickfont_size=x_font_size, title='', tickfont_family='Arial')

    fig.update_yaxes(
        tickfont_size=y_font_size,
        title=gene_title,
        titlefont_size=gene_font_size,
        range=[y_axis_range.get('min'), y_axis_range.get('max')],
        tickfont_family='Arial',
        titlefont_family='Arial'
    )

    # Add statistical annotations if test is specified
    if stat_test in ['t', 'student_t', 'u', 'wilcoxon']:
        # Use sorted order for groups if available
        if color_split:
            # For color-split, we need to perform tests between colors within each condition
            groups = list(df_show.T['Color_group'].unique())
        else:
            # Otherwise, perform tests between groups
            groups = list(df_show.T['Group'].unique())
        if st.session_state.get('sorted_groups') is not None:
            groups = st.session_state.sorted_groups

        # Perform statistical tests
        test_results = perform_multiple_tests(df_show, gene, stat_test, groups, adjust_method)

        # Get the maximum y value from the data for positioning annotations
        y_values = pd.to_numeric(df_show.T[gene], errors='coerce').dropna().values
        if len(y_values) > 0:
            if y_max is not None:
                max_y = y_max
            else:
                max_y = np.max(y_values)

            # Move bar position closer to box (1.15 -> 1.05)
            annotation_start_y = max_y * 1.05

            # Allow margin for multiple bars by considering bar height
            # Automatically adjust plot height considering number of bars
            significant_pairs = sum(1 for _, p in test_results.items() if p < 0.05 or show_all_p)
            extra_height = significant_pairs * stat_bar_height * 1.5  # Total height of bars and spaces

            fig, _ = add_stat_annotations(fig, test_results, annotation_start_y, groups,
                                         bar_height=stat_bar_height,
                                         p_value_font_size=p_value_font_size,
                                         show_all_p=show_all_p)

            # Adjust plot height to accommodate annotations
            fig.update_layout(height=py_y_size + extra_height)

    return fig



# Further improved statistical bar display function for combined plots - align bar position with boxes
def add_combined_stat_annotations(fig, test_results, idx, gene, df_show, y_max, groups, stat_bar_height, p_value_font_size, show_all_p):
    """
    Add statistical annotations to a combined plot for a specific gene index

    Args:
        fig: Plotly figure object
        test_results: Dictionary with group pairs as keys and p-values as values
        idx: Index of the gene plot
        gene: Current gene name
        df_show: DataFrame with the data
        y_max: Maximum y value for this gene
        groups: List of groups in the order they appear on the X-axis
        stat_bar_height: Height of the annotation bar (relative to data scale)
        p_value_font_size: Font size for p-value text
        show_all_p: Whether to show all p-values, including non-significant ones

    Returns:
        fig: Updated figure with annotations
    """
    # Use the provided groups list to maintain X-axis order
    groups_list = list(groups)

    # Calculate offsets to shift bar height for each group pair
    group_pairs = []
    for (group1, group2), p_value in test_results.items():
        if not np.isnan(p_value):
            if p_value < 0.05 or show_all_p:
                group_pairs.append((group1, group2))

    # Calculate offset according to number of bars (set spacing between bars)
    gap = stat_bar_height * 0.6  # Set spacing between bars
    bar_offsets = {}
    current_offset = 0
    for i, (group1, group2) in enumerate(group_pairs):
        bar_offsets[(group1, group2)] = current_offset
        current_offset += gap

    # Calculate relative bar height for each gene's value range
    y_values = pd.to_numeric(df_show.T[gene], errors='coerce').dropna().values
    if len(y_values) > 0:
        data_range = y_max - min(y_values)
        rel_bar_height = data_range * 0.1  # Set bar height to 10% of data range
        rel_gap = rel_bar_height * 0.6
    else:
        rel_bar_height = stat_bar_height
        rel_gap = gap

    for (group1, group2), p_value in test_results.items():
        if not np.isnan(p_value):
            # Skip non-significant cases (depending on option)
            if p_value >= 0.05 and not show_all_p:
                continue

            # Apply position offset for this bar
            current_offset = bar_offsets.get((group1, group2), 0) * rel_gap / gap

            # Place bar just above the maximum value of data
            current_y_pos = y_max * 1.02 + current_offset

            # Get group information
            # Use sorted order for positions if available
            position_groups = list(df_show.T['Group'].unique())
            if st.session_state.get('sorted_groups') is not None:
                position_groups = st.session_state.sorted_groups
            if group1 in position_groups and group2 in position_groups:
                # Get position information of groups within subplot
                group_positions = {}
                for i, g in enumerate(position_groups):
                    group_positions[g] = i / (len(position_groups) - 1) if len(position_groups) > 1 else 0.5

                # Calculate relative position of each group
                x1_pos = group_positions.get(group1, 0)
                x2_pos = group_positions.get(group2, 1)
            else:
                # Use default positions
                x1_pos = 0
                x2_pos = 1

            # Add horizontal bars - align to box positions
            fig.add_shape(
                type="line",
                x0=x1_pos, x1=x1_pos,
                y0=current_y_pos, y1=current_y_pos + rel_bar_height,
                line=dict(color="black", width=1.5),
                xref=f'x{idx+1} domain' if idx > 0 else 'x domain',
                yref=f'y{idx+1}' if idx > 0 else 'y'
            )

            fig.add_shape(
                type="line",
                x0=x1_pos, x1=x2_pos,
                y0=current_y_pos + rel_bar_height, y1=current_y_pos + rel_bar_height,
                line=dict(color="black", width=1.5),
                xref=f'x{idx+1} domain' if idx > 0 else 'x domain',
                yref=f'y{idx+1}' if idx > 0 else 'y'
            )

            fig.add_shape(
                type="line",
                x0=x2_pos, x1=x2_pos,
                y0=current_y_pos, y1=current_y_pos + rel_bar_height,
                line=dict(color="black", width=1.5),
                xref=f'x{idx+1} domain' if idx > 0 else 'x domain',
                yref=f'y{idx+1}' if idx > 0 else 'y'
            )

            # Display p-value text at center of bar
            p_text = f"{p_value:.4g}"
            fig.add_annotation(
                x=(x1_pos + x2_pos) / 2,
                y=current_y_pos + rel_bar_height,
                text=p_text,
                showarrow=False,
                yshift=5,
                font=dict(size=p_value_font_size),
                xref=f'x{idx+1} domain' if idx > 0 else 'x domain',
                yref=f'y{idx+1}' if idx > 0 else 'y',
                xanchor="center",
                yanchor="bottom"
            )

    return fig

def create_combined_plots(df_show, genes, plot_type='box', color_split=False, points='all', c_choice=None,
                        py_x_size=None, py_y_size=400,
                        legend_font_size=12, x_font_size=12, y_font_size=12, gene_font_size=16,
                        jitter=0.3, y_min=None, y_max=None, show_box=True, marker_size=6,
                        stat_test=None, adjust_method=None, apply_logit=False, p_value_font_size=10,
                        show_all_p=False, stat_bar_height=15):

    # Apply logit transformation if specified
    plot_data = df_show.copy()
    if apply_logit:
        for gene in genes:
            if gene in plot_data.index:
                try:
                    # Get column data and transform it
                    values = pd.to_numeric(plot_data.T[gene], errors='coerce')
                    # Apply logit transform only if values are between
                    if (values >= 0).all() and (values <= 1).all():
                        plot_data.loc[gene, :] = logit_transform(plot_data.T, gene)
                except:
                    pass

    # Calculate value range for each gene
    gene_y_ranges = {}
    for gene in genes:
        values = pd.to_numeric(df_show.T[gene].values, errors='coerce')
        values = values[~np.isnan(values)]
        if len(values) > 0:
            gene_min = float(min(values))
            gene_max = float(max(values))
            range_padding = (gene_max - gene_min) * 0.2  # Add 20% padding
            gene_y_ranges[gene] = (gene_min - range_padding, gene_max + range_padding)
        else:
            gene_y_ranges[gene] = (None, None)

    # When using common Y scale
    if y_min is None or y_max is None:
        all_values = []
        for gene in genes:
            values = pd.to_numeric(df_show.T[gene].values, errors='coerce')
            values = values[~np.isnan(values)]
            all_values.extend(values)

        if all_values:
            global_min = float(min(all_values))
            global_max = float(max(all_values))

            range_padding = (global_max - global_min) * 0.05
            if y_min is None:
                y_min = global_min - range_padding
            if y_max is None:
                y_max = global_max + range_padding

    if py_x_size is None:
        py_x_size = 400 * len(genes)

    fig = go.Figure()
    plot_width = 1.0 / len(genes)

    # Store statistical bar information for each gene
    stat_annotations_info = []

    # Use sorted order for traces if available, otherwise use data order
    if color_split:
        group_order_for_traces = list(df_show.T['Color_group'].unique())
    else:
        group_order_for_traces = list(df_show.T['Group'].unique())
    if st.session_state.get('sorted_groups') is not None:
        group_order_for_traces = st.session_state.sorted_groups

    for idx, gene in enumerate(genes):
        x_domain = [idx * plot_width, (idx + 1) * plot_width - 0.02]

        # Get value range for this gene
        gene_range = gene_y_ranges[gene]

        if color_split:
            for color_idx, color_group in enumerate(group_order_for_traces):
                group_data = df_show.T[df_show.T['Color_group'] == color_group]
                group_data[gene] = pd.to_numeric(group_data[gene], errors='coerce')

                if plot_type == 'box':
                    fig.add_trace(go.Box(
                        y=group_data[gene],
                        x=group_data['Condition'],
                        name=color_group,
                        marker_color=c_choice[color_idx % len(c_choice)],
                        boxpoints=points,
                        jitter=jitter,
                        pointpos=0,
                        marker_size=marker_size,
                        showlegend=(idx == 0),
                        xaxis=f'x{idx+1}' if idx > 0 else 'x',
                        yaxis=f'y{idx+1}' if idx > 0 else 'y'
                    ))
                else:  # violin
                    fig.add_trace(go.Violin(
                        y=group_data[gene],
                        x=group_data['Condition'],
                        name=color_group,
                        marker_color=c_choice[color_idx % len(c_choice)],
                        points=points,
                        jitter=jitter,
                        pointpos=0,
                        marker_size=marker_size,
                        showlegend=(idx == 0),
                        xaxis=f'x{idx+1}' if idx > 0 else 'x',
                        yaxis=f'y{idx+1}' if idx > 0 else 'y',
                        box_visible=show_box,
                        meanline_visible=True
                    ))
        else:
            for group_idx, group in enumerate(group_order_for_traces):
                group_data = df_show.T[df_show.T['Group'] == group]
                group_data[gene] = pd.to_numeric(group_data[gene], errors='coerce')

                if plot_type == 'box':
                    fig.add_trace(go.Box(
                        y=group_data[gene],
                        name=group,
                        marker_color=c_choice[group_idx % len(c_choice)],
                        boxpoints=points,
                        jitter=jitter,
                        pointpos=0,
                        marker_size=marker_size,
                        showlegend=(idx == 0),
                        xaxis=f'x{idx+1}' if idx > 0 else 'x',
                        yaxis=f'y{idx+1}' if idx > 0 else 'y'
                    ))
                else:  # violin
                    fig.add_trace(go.Violin(
                        y=group_data[gene],
                        name=group,
                        marker_color=c_choice[group_idx % len(c_choice)],
                        points=points,
                        jitter=jitter,
                        pointpos=0,
                        marker_size=marker_size,
                        showlegend=(idx == 0),
                        xaxis=f'x{idx+1}' if idx > 0 else 'x',
                        yaxis=f'y{idx+1}' if idx > 0 else 'y',
                        box_visible=show_box,
                        meanline_visible=True
                    ))

        # Set Y-axis range for each subplot - use individual ranges
        # Allow margin at top to accommodate statistical bars
        y_axis_min = gene_range[0] if y_min is None else y_min
        y_axis_max = gene_range[1] if y_max is None else y_max

        # Set X-axis category order
        xaxis_settings = {
            'domain': x_domain,
            'showticklabels': True,
            'tickfont_size': x_font_size,
            'title': '',
            'tickangle': 45,
            'tickfont_family': 'Arial'
        }
        if not color_split:
            category_order = list(df_show.T['Group'].unique())
            # Use manually sorted order if available
            if st.session_state.get('sorted_groups') is not None:
                category_order = st.session_state.sorted_groups

            xaxis_settings['categoryorder'] = 'array'
            xaxis_settings['categoryarray'] = category_order
        else:
            category_order = list(df_show.T['Color_group'].unique())
            # Use manually sorted order if available
            if st.session_state.get('sorted_groups') is not None:
                category_order = st.session_state.sorted_groups

        # Calculate Y-axis range (allow space for statistical bars)
        y_range = y_axis_max - y_axis_min
        y_upper = y_axis_max + y_range * 0.25  # Add 25% of range to top

        fig.update_layout(**{
            f'xaxis{idx+1 if idx > 0 else ""}': xaxis_settings,
            f'yaxis{idx+1 if idx > 0 else ""}': {
                'title': '',
                'titlefont_size': gene_font_size,
                'tickfont_size': y_font_size,
                'range': [y_axis_min, y_upper],
                'tickfont_family': 'Arial'
            }
        })

        # Store statistical test information
        if stat_test in ['t', 'student_t', 'u', 'wilcoxon']:
            # Use sorted order for groups if available
            if color_split:
                groups = list(df_show.T['Color_group'].unique())
            else:
                groups = list(df_show.T['Group'].unique())
            if st.session_state.get('sorted_groups') is not None:
                groups = st.session_state.sorted_groups

            test_results = perform_multiple_tests(df_show, gene, stat_test, groups, adjust_method)

            # Calculate actual maximum value for each gene
            gene_max_values = pd.to_numeric(df_show.T[gene], errors='coerce').dropna().values
            if len(gene_max_values) > 0:
                gene_max_actual = float(max(gene_max_values))
            else:
                gene_max_actual = gene_range[1] if y_max is None else y_max

            stat_annotations_info.append((idx, gene, test_results, gene_max_actual, groups))

    # Update overall plot layout
    annotations = []
    for idx, gene in enumerate(genes):
        x_pos = (idx * plot_width) + (plot_width / 2)
        annotations.append(dict(
            x=x_pos,
            y=1.1,
            xref='paper',
            yref='paper',
            text=f'<i>{gene}</i>',
            showarrow=False,
            font=dict(size=gene_font_size, family='Arial'),
            xanchor='center',
            yanchor='bottom'
        ))

    fig.update_layout(
        plot_bgcolor='rgb(255, 255, 255)',
        paper_bgcolor='rgb(255, 255, 255)',
        width=py_x_size,
        height=py_y_size * 1.5,  # Ensure sufficient height for statistical bars
        legend_font=dict(size=legend_font_size, family='Arial'),
        showlegend=True,
        annotations=annotations,
        font=dict(family='Arial')
    )

    # After all graph drawing is complete, add statistical bars
    for idx, gene, test_results, gene_max, groups in stat_annotations_info:
        fig = add_combined_stat_annotations(
            fig, test_results, idx, gene, df_show, gene_max, groups,
            stat_bar_height, p_value_font_size, show_all_p
        )

    return fig


@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t', header = None).encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header = 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep, engine='python')
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

# Initialize session state
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir, res_dir = mk_temp_dir("Boxplot")
    st.session_state.temp_dir = temp_dir
else:
    temp_dir = st.session_state.temp_dir
    temp_dir, res_dir = mk_temp_dir("Boxplot", temp_dir)

if "set_group" not in st.session_state:
    st.session_state.set_group = False

if 'user_input' not in st.session_state:
    st.session_state.user_input = ''

if 'sorted_groups' not in st.session_state:
    st.session_state.sorted_groups = None

st.markdown("## Box/Violin plot generator")

# Data structure selection
st.markdown("### Data Structure")
data_structure = st.radio(
    "Select data structure:",
    ("Matrix (rows=genes, cols=samples)", "Long format (rows=observations)"),
    horizontal=True,
)
is_long_format = data_structure.startswith("Long")

# Show example data format
with st.expander("Example data format", expanded=False):
    if is_long_format:
        st.markdown("**Long format example:**")
        example_long = pd.DataFrame({
            'Gene': ['GeneA', 'GeneA', 'GeneA', 'GeneB', 'GeneB', 'GeneB'],
            'CellType': ['HSC', 'HSC', 'MPP', 'HSC', 'HSC', 'MPP'],
            'Value': [0.85, 0.82, 0.91, 0.73, 0.75, 0.88]
        })
        st.dataframe(example_long, hide_index=True)
        st.markdown("""
        - **Category column** (optional): X-axis grouping (e.g., Gene, CellType)
        - **Value column** (required): Numeric data
        - **Group column** (optional): Color grouping axis (e.g., CellType, Condition)

        **ISP result example:**
        - 1 gene, multiple cell types -> Category: perturbed_gene or (None), Group: cell_type
        - Multiple genes, multiple cell types -> Category: perturbed_gene, Group: cell_type
        """)
    else:
        st.markdown("**Matrix format example:**")
        example_matrix = pd.DataFrame({
            'Gene': ['GeneA', 'GeneB', 'GeneC'],
            'Control_1': [10.5, 5.2, 20.1],
            'Control_2': [12.3, 6.1, 22.5],
            'Treatment_1': [8.7, 4.8, 18.9],
            'Treatment_2': [9.1, 5.0, 19.2]
        })
        st.dataframe(example_matrix, hide_index=True)
        st.markdown("""
        - **Rows**: Gene names
        - **Columns**: Sample names (Group_Number)
        - **Values**: Expression levels, etc.
        """)

use_upload = 'Yes'
if 'df' in st.session_state:
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name

if use_upload == 'Yes':
    input_file_type = st.radio(
            "File format:",
            ('auto', 'tsv','csv', 'excel'))

    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        if input_file_type == "auto":
            try:
                df = read_csv(uploaded_file, header=None, index_col=None, sep=None)
            except:
                df = read_excel(uploaded_file, index_col=None, header=None)

        elif input_file_type == "csv":
            df = read_csv(uploaded_file, header=None, index_col=None)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep='\t', header=None, index_col=None)
        else:
            df = read_excel(uploaded_file, index_col=None, header=None)

        st.write('Data Dimension: ' + str(df.shape))

        # ========================================
        # Long format handling
        # ========================================
        if is_long_format:
            # For Long format: first row is header
            df.columns = df.iloc[0, :].tolist()
            df = df.drop(0, axis=0).reset_index(drop=True)

            st.write("**Preview (Long format):**")
            st.write(df.head(5))

            content = df.columns.tolist()

            st.markdown("#### Column Selection")

            # Value column (numeric data) - Required
            numeric_cols = []
            for c in content:
                try:
                    pd.to_numeric(df[c], errors='raise')
                    numeric_cols.append(c)
                except:
                    pass

            # Try to find a good default for value column
            default_value_cols = ['Value', 'value', 'Expression', 'expression', 'cosine_similarity', 'Score', 'score']
            default_val_idx = 0
            for dv in default_value_cols:
                if dv in numeric_cols:
                    default_val_idx = numeric_cols.index(dv)
                    break

            value_col = st.selectbox(
                "Value column (required)",
                numeric_cols if numeric_cols else content,
                index=default_val_idx,
                help="Numeric column to plot",
                key="long_value_col"
            )

            col1, col2 = st.columns(2)

            with col1:
                # Category column (optional) - X-axis grouping
                category_options = ["(None)"] + content
                default_cat_cols = ['Gene', 'gene', 'perturbed_gene', 'Category', 'Name']
                default_cat_idx = 0  # Default to "(None)"
                for dc in default_cat_cols:
                    if dc in content:
                        default_cat_idx = category_options.index(dc)
                        break

                category_col = st.selectbox(
                    "Category column (optional)",
                    category_options,
                    index=default_cat_idx,
                    help="X-axis grouping (e.g., Gene, CellType)",
                    key="long_category_col"
                )

            with col2:
                # Group column (optional) - Color grouping
                group_options = ["(None)"] + content
                default_group_cols = ['Group', 'group', 'Condition', 'condition', 'cell_type', 'CellType', 'Treatment', 'treatment']
                default_grp_idx = 0  # Default to "(None)"
                for dg in default_group_cols:
                    if dg in content:
                        default_grp_idx = group_options.index(dg)
                        break

                group_col = st.selectbox(
                    "Group column (optional)",
                    group_options,
                    index=default_grp_idx,
                    help="Color grouping axis (e.g., CellType, Condition)",
                    key="long_group_col"
                )

            # Handle "(None)" selections
            if category_col == "(None)":
                category_col = None
            if group_col == "(None)":
                group_col = None

            # Validate: at least one of category or group must be specified
            if category_col is None and group_col is None:
                st.error("Please specify at least one of Category or Group.")
                st.stop()

            # Convert Long format to Matrix format for compatibility
            st.markdown("#### Data Conversion")

            # Ensure value column is numeric
            df[value_col] = pd.to_numeric(df[value_col], errors='coerce')

            # Determine the effective category for X-axis
            if category_col is None and group_col is not None:
                # Use group as the X-axis (e.g., cell_type comparison)
                effective_category = group_col
                effective_group = None
                st.info(f"X-axis: comparing by {group_col}")
            elif category_col is not None and group_col is None:
                # Use category as X-axis, no color grouping
                effective_category = category_col
                effective_group = None
                st.info(f"X-axis: comparing by {category_col}")
            else:
                # Both specified
                effective_category = category_col
                effective_group = group_col
                st.info(f"X-axis: {category_col}, Color: {group_col}")

            # Get unique categories
            unique_categories = df[effective_category].unique().tolist()

            if effective_group is not None:
                # With group column for color
                unique_groups = df[effective_group].unique().tolist()
                st.write(f"Categories: {len(unique_categories)}, Groups: {len(unique_groups)}")
                # Create sample names: Group_Index
                df['_sample_name'] = df[effective_group].astype(str) + '_' + df.groupby(effective_group).cumcount().add(1).astype(str)
            else:
                # Without group column
                st.write(f"Categories: {len(unique_categories)}")
                # Create sample names: Category_Index
                df['_sample_name'] = df[effective_category].astype(str) + '_' + df.groupby(effective_category).cumcount().add(1).astype(str)

            # Pivot to matrix format
            if len(unique_categories) > 1:
                # Multiple categories - traditional matrix
                df_matrix = df.pivot_table(
                    index=effective_category,
                    columns='_sample_name',
                    values=value_col,
                    aggfunc='first'
                )
            else:
                # Single category - create matrix with one row
                single_cat = unique_categories[0]
                df_matrix = pd.DataFrame(
                    [df[value_col].tolist()],
                    index=[single_cat],
                    columns=df['_sample_name'].tolist()
                )

            df = df_matrix
            st.success(f"Converted to Matrix format: {df.shape[0]} categories x {df.shape[1]} samples")

            file_name_head = os.path.splitext(uploaded_file.name)[0]
            st.session_state.df = df
            st.session_state.uploaded_file_name = file_name_head
            st.session_state.is_long_format = True
            st.session_state.long_format_info = {
                'category_col': category_col,
                'value_col': value_col,
                'group_col': group_col,
                'effective_category': effective_category,
                'effective_group': effective_group
            }

        # ========================================
        # Matrix format handling (original logic)
        # ========================================
        else:
            transpose_df = st.checkbox('Transpose the data?')
            if transpose_df:
                df = df.T

            df.columns = df.iloc[0,:].tolist()
            df = df.drop(0, axis=0)
            content = df.columns.tolist()
            Gene_column = content[0]

            if "Annotation/Divergence" in content:
                search_word = '([^\ \(]+).*'
                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str)

                df = df.loc[:,'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                Gene_column = "Gene"
                st.write("Found Annotation/Divergence column.")
            elif "Gene" in content:
                Gene_column = "Gene"
            else:
                Gene_column = st.selectbox(
                    'Select gene column',
                    content)

            df = df.set_index(Gene_column)

            pattern = "^([^|]*)"
            repatter = re.compile(pattern)
            f_annotation = lambda x: repatter.match(x).group(1) if repatter.match(x) else x

            original_index = list(df.index)
            new_index = [f_annotation(x) for x in original_index]
            df.index = new_index

            if original_index != new_index:
                st.warning("Removed content after '|' from gene symbols.")

            file_name_head = os.path.splitext(uploaded_file.name)[0]
            st.session_state.df = df
            st.session_state.uploaded_file_name = file_name_head
            st.session_state.is_long_format = False
    else:
        sys.exit(1)

if df is not None:
    st.write(df.head(3))

    st.markdown('---')
    st.markdown("##### Filter and transform data?")
    calc_z = False
    center0_z = False
    howlog = 'No'
    Manip = st.checkbox('', label_visibility = 'collapsed')
    if Manip:
        f_inf = -float('inf')
        p_inf = float('inf')

        st.markdown("######   ")
        calc_div = st.checkbox('Divided by (e.g., 1000, 1000,000)?', value = False)
        if calc_div:
            div_unit =  int(st.text_input("Unit: ",  value=1))
            df = df/div_unit
        calc_log = st.checkbox('Log transformation?')
        if calc_log:
            howlog = st.radio('Method', ['No', 'log2+1', 'log2', 'loge+1', 'loge','log10+1','log10', 'asinh'])

        st.markdown("######   ")
        calc_z = st.checkbox('Z-score?')

        st.markdown("######   ")

        st.write(df.tail())
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

        if calc_z:
            center0_z= True
            df_z = df.copy()
            m = df_z.mean(1)
            s = df_z.std(1)
            df_z = df_z.sub(m, axis=0).div(s, axis = 0)
            df_z = np.round(df_z, decimals=10)
            df_z = df_z.loc[~(df_z==0).all(axis=1)]
            df_z = df_z.dropna(how='any', axis=0)
            df = df_z

    # Add the logit transformation option
    st.markdown("##### Data type transformation")
    apply_logit = st.checkbox('Apply logit transformation? (for ratio data between 0-1)')

    st.markdown('---')

    condition = [str(i) for i in df.columns.tolist()]
    group_condition = [remove_after_space(x) for x in condition]
    group_condition = [remove_sample_num(x) for x in group_condition]

    reduced_condition = [remove_underscore(x) for x in group_condition]
    color_condition = [select_underscore(x) for x in group_condition]

    st.markdown("##### Cleaned data:")
    st.write(df.head(3))
    st.markdown('---')

    with st.sidebar:
        st.markdown("### Statistical analysis")
        perform_stats = st.checkbox('Perform statistical tests?')

        if perform_stats:
            stat_test = st.radio(
                "Test type:",
                ('t', 'student_t', 'u', 'wilcoxon'),
                format_func=lambda x: {'t': 't-test (Welch)',
                                       'student_t': "Student's t-test",
                                       'u': 'Mann-Whitney U test',
                                       'wilcoxon': 'Wilcoxon test'}[x],
                                       index=1

            )


            adjust_method = st.radio(
                "P-value adjustment:",
                (None, 'bh', 'bonferroni', 'holm', 'by', 'tukey'),
                format_func=lambda x: {None: 'None',
                                      'bh': 'Benjamini-Hochberg (FDR)',
                                      'bonferroni': 'Bonferroni',
                                      'holm': 'Holm-Bonferroni',
                                      'by': 'Benjamini-Yekutieli',
                                      'tukey': 'Tukey-Kramer'}[x]
            )

            # Add p-value font size option
            p_value_font_size = st.number_input("P-value font size:", min_value=3, max_value=24, value=8)

            # Add options for p-value display
            show_all_p = st.checkbox("Show all p-values (including non-significant)", value=False)

            # Add bar height option
            stat_bar_height = st.number_input("Statistical bar height:", min_value=0.1, max_value=20.0, value=10.0)
            st.markdown("---")
        else:
            stat_test = None
            adjust_method = None
            p_value_font_size = 10
            show_all_p = False
            stat_bar_height = 20

        st.markdown("### Plot type")
        plot_type = st.radio("Plot type:", ('Individual plots', 'Combined plot'), index=0, label_visibility="collapsed")

        st.markdown("### Plot style")
        plot_style = st.radio("Plot style:", ('Box', 'Violin'), index=0)
        plot_style = plot_style.lower()  # 'box' or 'violin'

        # Show box option for violin plot
        show_box = True
        if plot_style == 'violin':
            show_box = st.checkbox('Show box in violin plot?', value=True)

        st.markdown("### Data point settings")
        show_point = st.checkbox('Show data points?')
        if show_point:
            points = 'all'
            jitter = st.number_input("Jitter", min_value=0.0, max_value=0.5, value=0.2)
            marker_size = st.number_input("Marker size", min_value=1, max_value=20, value=6)
        else:
            points = 'outliers'
            jitter = 0.2
            marker_size = 6

        st.markdown("### Y-axis range")
        y_min = st.number_input("Y-axis min", value=None)
        y_max = st.number_input("Y-axis max", value=None)

        st.markdown("### Plot color")
        color_choice = st.selectbox('Color sequence:',
            ('Plotly','D3','G10','T10','Alphabet','Dark24','Light24','Set1','Pastel1','Dark2','Set2',
            'Pastel2','Set3','Antique','Bold','Pastel','Prism','Safe','Vivid'), key = 'Plotly')
        c_choice = getattr(px.colors.qualitative, color_choice)
        show_color = st.checkbox('Show color sequence?')
        if show_color:
            color_fig = plt.figure()
            color_fig = px.colors.qualitative.swatches()
            st.plotly_chart(color_fig)

        st.markdown('#### Plot size')
        py_x_size = float(st.text_input("Plot x size:", value = 600))
        py_y_size = float(st.text_input("Plot y size:", value = 400))
        st.markdown('#### Font size')
        st.markdown('##### Y-axis')
        gene_font_size = float(st.text_input("Gene name font size:", value = 16))
        y_font_size = float(st.text_input("Y-axis ticks font size:", value = 12))
        st.markdown('##### X-axis')
        x_font_size = float(st.text_input("Group name font size:", value = 12))
        st.markdown('##### Legend')
        legend_font_size = float(st.text_input("Legend font size:", value = 12))

    save_type = st.radio("Save plot as: (Preparing PDF may take a time.)", ('pdf','png','html'))
    color_split = st.checkbox("Grouped by 'condition x group_color'?")

    df_e = pd.DataFrame(index = condition, columns = ['Group', 'Condition', 'Color_group'])
    df_e["Group"] = group_condition
    df_e["Condition"] = reduced_condition
    df_e["Color_group"] = color_condition
    df_show = df.copy(deep=True)

    with st.form("Group"):
        st.write('Set group and color:')
        df_e = st.data_editor(df_e)
        submitted = st.form_submit_button("Set groups")
        st.session_state.set_group = True
        df_show.loc['Group',:] = df_e['Group'].tolist()
        df_show.loc['Color_group',:] = df_e['Color_group'].tolist()
        df_show.loc['Condition',:] = df_e['Condition'].tolist()

    if submitted or st.session_state.set_group:
        # Add checkbox for changing group order
        change_group_order = st.checkbox("Change group order?", value=False)

        if change_group_order:
            # Get current groups based on color_split setting
            if color_split:
                current_groups = list(df_show.T['Color_group'].unique())
            else:
                current_groups = list(df_show.T['Group'].unique())

            # Initialize sorted_groups with current groups if not set
            if st.session_state.sorted_groups is None:
                st.session_state.sorted_groups = current_groups.copy()

            st.markdown("##### Sort groups (drag and drop):")
            with st.form("Group_Sorter"):
                sorted_groups = sort_items(st.session_state.sorted_groups)
                submitted_sort = st.form_submit_button("Apply group order")
                if submitted_sort:
                    st.session_state.sorted_groups = sorted_groups
                    st.rerun()
        else:
            # Reset sorted_groups if checkbox is not checked
            st.session_state.sorted_groups = None

        # For Long format, automatically use all categories
        # For Matrix format, allow user to specify genes
        is_long = st.session_state.get('is_long_format', False)

        if is_long:
            # Long format: use all categories from the data
            gene_subset = df.index.tolist()
            st.info(f"Long format: {len(gene_subset)} categories")

            if st.button('Make plot'):
                if plot_type == 'Individual plots':
                    for i in gene_subset:
                        gene_name = "___" + i + "___"
                        gene_title = "<i>" + i + "</i>"
                        st.markdown(gene_name)

                        fig = create_plot_with_points(df_show, i,
                                                    plot_type=plot_style,
                                                    color_split=color_split,
                                                    points=points,
                                                    c_choice=c_choice,
                                                    py_x_size=py_x_size,
                                                    py_y_size=py_y_size,
                                                    legend_font_size=legend_font_size,
                                                    x_font_size=x_font_size,
                                                    y_font_size=y_font_size,
                                                    gene_font_size=gene_font_size,
                                                    jitter=jitter,
                                                    y_min=y_min,
                                                    y_max=y_max,
                                                    show_box=show_box,
                                                    marker_size=marker_size,
                                                    stat_test=stat_test,
                                                    adjust_method=adjust_method,
                                                    apply_logit=apply_logit,
                                                    p_value_font_size=p_value_font_size,
                                                    show_all_p=show_all_p,
                                                    stat_bar_height=stat_bar_height)

                        st.plotly_chart(fig, theme=None)
                        if save_type == 'html':
                            fig.write_html(res_dir + "/" + i + "." + save_type)
                        else:
                            fig.write_image(res_dir + "/" + i + "." + save_type, format=save_type)

                else:  # Combined plot
                    combined_x_size = py_x_size * len(gene_subset)

                    fig = create_combined_plots(df_show, gene_subset,
                                              plot_type=plot_style,
                                              color_split=color_split,
                                              points=points,
                                              c_choice=c_choice,
                                              py_x_size=combined_x_size,
                                              py_y_size=py_y_size,
                                              legend_font_size=legend_font_size,
                                              x_font_size=x_font_size,
                                              y_font_size=y_font_size,
                                              gene_font_size=gene_font_size,
                                              jitter=jitter,
                                              y_min=y_min,
                                              y_max=y_max,
                                              show_box=show_box,
                                              marker_size=marker_size,
                                              stat_test=stat_test,
                                              adjust_method=adjust_method,
                                              apply_logit=apply_logit,
                                              p_value_font_size=p_value_font_size,
                                              show_all_p=show_all_p,
                                              stat_bar_height=stat_bar_height)

                    st.plotly_chart(fig, theme=None)

                    if save_type == 'html':
                        fig.write_html(res_dir + "/combined_plot." + save_type)
                    else:
                        fig.write_image(res_dir + "/combined_plot." + save_type, format=save_type)

                # Zip file creation and download button
                shutil.make_archive(temp_dir + "/Boxplot", format='zip', root_dir=res_dir)
                with open(temp_dir + "/Boxplot.zip", "rb") as fp:
                    btn = st.download_button(
                        label="Download plots?",
                        data=fp,
                        file_name=file_name_head + ".Boxplot.zip",
                        mime="zip"
                    )

        else:
            # Matrix format: user specifies genes
            st.markdown("##### Genes (comma, space, CR separated):")
            genes = st.text_input("genes", label_visibility='collapsed', value=st.session_state.user_input)
            st.session_state.user_input = genes
            gene_list = []
            if len(genes) > 0:
                gene_list = genes.split(' ')
                gene_list = list(filter(lambda a: a != '', gene_list))
                if ',' in genes:
                    gene_list = sum([x.split(',') for x in gene_list], [])
                if '\t' in genes:
                    gene_list = sum([x.split('\t') for x in gene_list], [])
                if '\n' in genes:
                    gene_list = sum([x.split('\n') for x in gene_list], [])

                df_index_set = list(set(df.index.tolist()))
                df_index_set_lower = [x.lower() for x in df_index_set]
                gene_subset = [df_index_set[df_index_set_lower.index(x.lower())]
                              for x in gene_list if (x.lower() in df_index_set_lower)]
                gene_subset = sorted(set(gene_subset), key=gene_subset.index)

                if st.button('Make plot', key='matrix_make_plot'):
                    if plot_type == 'Individual plots':
                        for i in gene_subset:
                            gene_name = "___" + i + "___"
                            gene_title = "<i>" + i + "</i>"
                            st.markdown(gene_name)

                            fig = create_plot_with_points(df_show, i,
                                                        plot_type=plot_style,
                                                        color_split=color_split,
                                                        points=points,
                                                        c_choice=c_choice,
                                                        py_x_size=py_x_size,
                                                        py_y_size=py_y_size,
                                                        legend_font_size=legend_font_size,
                                                        x_font_size=x_font_size,
                                                        y_font_size=y_font_size,
                                                        gene_font_size=gene_font_size,
                                                        jitter=jitter,
                                                        y_min=y_min,
                                                        y_max=y_max,
                                                        show_box=show_box,
                                                        marker_size=marker_size,
                                                        stat_test=stat_test,
                                                        adjust_method=adjust_method,
                                                        apply_logit=apply_logit,
                                                        p_value_font_size=p_value_font_size,
                                                        show_all_p=show_all_p,
                                                        stat_bar_height=stat_bar_height)

                            st.plotly_chart(fig, theme=None)
                            if save_type == 'html':
                                fig.write_html(res_dir + "/" + i + "." + save_type)
                            else:
                                fig.write_image(res_dir + "/" + i + "." + save_type, format=save_type)

                    else:  # Combined plot
                        combined_x_size = py_x_size * len(gene_subset)

                        fig = create_combined_plots(df_show, gene_subset,
                                                  plot_type=plot_style,
                                                  color_split=color_split,
                                                  points=points,
                                                  c_choice=c_choice,
                                                  py_x_size=combined_x_size,
                                                  py_y_size=py_y_size,
                                                  legend_font_size=legend_font_size,
                                                  x_font_size=x_font_size,
                                                  y_font_size=y_font_size,
                                                  gene_font_size=gene_font_size,
                                                  jitter=jitter,
                                                  y_min=y_min,
                                                  y_max=y_max,
                                                  show_box=show_box,
                                                  marker_size=marker_size,
                                                  stat_test=stat_test,
                                                  adjust_method=adjust_method,
                                                  apply_logit=apply_logit,
                                                  p_value_font_size=p_value_font_size,
                                                  show_all_p=show_all_p,
                                                  stat_bar_height=stat_bar_height)

                        st.plotly_chart(fig, theme=None)

                        if save_type == 'html':
                            fig.write_html(res_dir + "/combined_plot." + save_type)
                        else:
                            fig.write_image(res_dir + "/combined_plot." + save_type, format=save_type)

                    # Zip file creation and download button
                    shutil.make_archive(temp_dir + "/Boxplot", format='zip', root_dir=res_dir)
                    with open(temp_dir + "/Boxplot.zip", "rb") as fp:
                        btn = st.download_button(
                            label="Download plots?",
                            data=fp,
                            file_name=file_name_head + ".Boxplot.zip",
                            mime="zip"
                        )
