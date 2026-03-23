"""
COMPASS Postprocessing - Differential Metabolic Activity Analysis
Equivalent to compassR functionality implemented in Python
"""

import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from io import BytesIO
import zipfile
import os
from pathlib import Path


def load_subsystem_annotations():
    """Load subsystem annotations from Human-GEM complete mapping"""
    # Primary source: Complete Human1 subsystem mapping (147 subsystems, 12971 reactions)
    # Use relative path from script location
    script_dir = os.path.dirname(os.path.abspath(__file__))
    primary_path = os.path.join(script_dir, '..', 'data', 'human1_subsystem_mapping.csv')

    if os.path.exists(primary_path):
        try:
            df = pd.read_csv(primary_path)
            # Ensure reaction_name column exists
            if 'reaction_name' not in df.columns:
                df['reaction_name'] = ''
            return df
        except Exception as e:
            st.warning(f"Error loading subsystem mapping: {e}")
            pass

    # Fallback: COMPASS core reactions (limited coverage)
    _compass_base = Path.home() / "anaconda3" / "envs" / "shiny" / "lib" / "python3.11" / "site-packages" / "compass" / "Resources" / "Metabolic Models"
    model_paths = [
        str(_compass_base / "Human1" / "core_reactions_md.csv"),
        str(_compass_base / "RECON2_mat" / "model" / "core_reactions_md.csv"),
    ]

    for path in model_paths:
        if os.path.exists(path):
            try:
                df = pd.read_csv(path)
                # Create MAR ID to subsystem mapping
                annot = df[['ID', 'SUBSYSTEM', 'NAME']].copy()
                annot.columns = ['reaction_id', 'subsystem', 'reaction_name']
                annot = annot.dropna(subset=['subsystem'])
                return annot
            except Exception as e:
                continue
    return None


@st.cache_data
def get_subsystem_mapping(_cache_key=None):
    """Get cached subsystem mapping"""
    return load_subsystem_annotations()


def clear_subsystem_cache():
    """Clear the subsystem mapping cache"""
    get_subsystem_mapping.clear()

st.set_page_config(page_title="COMPASS Postprocess", page_icon="📊", layout="wide")

st.title("COMPASS Postprocessing")
st.markdown("""
Differential analysis tool for COMPASS results (Python implementation equivalent to compassR)

### Features
- **Wilcoxon rank-sum test**: Compare reaction activity between groups
- **Effect size**: Choose between Cohen's d (parametric) or Rank-biserial correlation (non-parametric)
- **Multiple testing correction**: Benjamini-Hochberg method
- **Volcano plot**: Visualization of results
- **Subsystem aggregation**: Pathway-level analysis (when annotations are available)
""")

# ========================================
# Step 1: Load Data
# ========================================
st.header("Step 1: Load data")

upload_method = st.radio(
    "Data loading method",
    ["compassR format (3 files)", "reactions.tsv only"],
    help="compassR format: cell_metadata.csv, reactions.tsv, linear_gene_expression_matrix.tsv"
)

if upload_method == "compassR format (3 files)":
    col1, col2, col3 = st.columns(3)

    with col1:
        cell_metadata_file = st.file_uploader(
            "cell_metadata.csv",
            type=['csv'],
            key="cell_metadata"
        )

    with col2:
        reactions_file = st.file_uploader(
            "reactions.tsv",
            type=['tsv', 'txt'],
            key="reactions"
        )

    with col3:
        expression_file = st.file_uploader(
            "linear_gene_expression_matrix.tsv (optional)",
            type=['tsv', 'txt'],
            key="expression"
        )

    if cell_metadata_file and reactions_file:
        # Load data
        cell_metadata = pd.read_csv(cell_metadata_file, index_col=0)
        reactions = pd.read_csv(reactions_file, sep='\t', index_col=0)

        # Store in session state
        st.session_state.compass_cell_metadata = cell_metadata
        st.session_state.compass_reactions = reactions

        if expression_file:
            expression = pd.read_csv(expression_file, sep='\t', index_col=0)
            st.session_state.compass_expression = expression

        st.success(f"Data loaded: {len(reactions)} reactions x {reactions.shape[1]} cells")

else:
    reactions_file = st.file_uploader(
        "reactions.tsv (COMPASS scores)",
        type=['tsv', 'txt'],
        key="reactions_only"
    )

    if reactions_file:
        reactions = pd.read_csv(reactions_file, sep='\t', index_col=0)
        st.session_state.compass_reactions = reactions

        # Create dummy metadata from column names
        cell_ids = reactions.columns.tolist()
        cell_metadata = pd.DataFrame({'cell_id': cell_ids}, index=cell_ids)

        # Try to infer groups from cell names
        if any('CTRL' in c or 'HFD' in c for c in cell_ids):
            cell_metadata['group'] = ['CTRL' if 'CTRL' in c else 'HFD' for c in cell_ids]

        st.session_state.compass_cell_metadata = cell_metadata
        st.success(f"Data loaded: {len(reactions)} reactions x {reactions.shape[1]} cells")

# ========================================
# Step 2: Configure Analysis
# ========================================
if 'compass_reactions' in st.session_state:
    st.header("Step 2: Analysis settings")

    reactions = st.session_state.compass_reactions
    cell_metadata = st.session_state.compass_cell_metadata

    # Show data preview
    with st.expander("Data preview"):
        st.write("**Cell metadata:**")
        st.dataframe(cell_metadata.head())
        st.write(f"Columns: {list(cell_metadata.columns)}")

        st.write("**Reactions (top 10):**")
        st.dataframe(reactions.iloc[:10, :5])

    # Select group column
    st.subheader("Group settings")

    group_col = st.selectbox(
        "Select group column",
        options=cell_metadata.columns.tolist(),
        index=0 if 'cond' not in cell_metadata.columns else cell_metadata.columns.tolist().index('cond'),
        help="Column containing the groups to compare"
    )

    # Get unique groups
    unique_groups = cell_metadata[group_col].unique().tolist()

    col1, col2 = st.columns(2)
    with col1:
        group1 = st.selectbox(
            "Group 1 (control)",
            options=unique_groups,
            index=0
        )
    with col2:
        group2 = st.selectbox(
            "Group 2 (treatment)",
            options=[g for g in unique_groups if g != group1],
            index=0
        )

    # Show group sizes
    g1_cells = cell_metadata[cell_metadata[group_col] == group1].index.tolist()
    g2_cells = cell_metadata[cell_metadata[group_col] == group2].index.tolist()

    # Align with reactions columns
    g1_cells = [c for c in g1_cells if c in reactions.columns]
    g2_cells = [c for c in g2_cells if c in reactions.columns]

    col1, col2 = st.columns(2)
    with col1:
        st.metric(f"{group1} cells", len(g1_cells))
    with col2:
        st.metric(f"{group2} cells", len(g2_cells))

    # Analysis parameters
    st.subheader("Analysis parameters")

    # Effect size selector
    effect_size_method = st.radio(
        "Effect size",
        ["Cohen's d", "Rank-biserial correlation"],
        horizontal=True,
        help="Cohen's d: parametric (assumes normality). Rank-biserial: non-parametric (derived from Mann-Whitney U, range -1 to +1)"
    )
    is_rank_biserial = (effect_size_method == "Rank-biserial correlation")
    effect_col = 'rank_biserial' if is_rank_biserial else 'cohens_d'
    effect_label = "Rank-biserial r" if is_rank_biserial else "Cohen's d"

    col1, col2 = st.columns(2)
    with col1:
        padj_threshold = st.number_input(
            "Adjusted p-value threshold",
            min_value=0.001,
            max_value=0.1,
            value=0.05,
            step=0.01
        )
    with col2:
        if is_rank_biserial:
            effect_threshold = st.number_input(
                f"{effect_label} threshold",
                min_value=0.0,
                max_value=1.0,
                value=0.3,
                step=0.05,
                key="compass_effect_rb"
            )
        else:
            effect_threshold = st.number_input(
                f"{effect_label} threshold",
                min_value=0.0,
                max_value=2.0,
                value=0.5,
                step=0.1,
                key="compass_effect_cd"
            )

    # ========================================
    # Step 3: Run Analysis
    # ========================================
    if st.button("Run Differential Analysis", type="primary"):
        st.header("Step 3: Running analysis")

        progress_bar = st.progress(0)
        status_text = st.empty()

        results = []
        n_reactions = len(reactions)

        for i, reaction in enumerate(reactions.index):
            if i % 100 == 0:
                progress_bar.progress(i / n_reactions)
                status_text.text(f"Processing {i}/{n_reactions} reactions...")

            g1_values = reactions.loc[reaction, g1_cells].values.astype(float)
            g2_values = reactions.loc[reaction, g2_cells].values.astype(float)

            # Skip if no variance
            if np.std(g1_values) == 0 and np.std(g2_values) == 0:
                continue

            # Mann-Whitney U test
            try:
                stat, pval = mannwhitneyu(g2_values, g1_values, alternative='two-sided')
            except:
                continue

            # Cohen's d
            pooled_std = np.sqrt(
                ((len(g2_values)-1)*np.std(g2_values, ddof=1)**2 +
                 (len(g1_values)-1)*np.std(g1_values, ddof=1)**2) /
                (len(g2_values) + len(g1_values) - 2)
            )
            if pooled_std > 0:
                cohens_d = (np.mean(g2_values) - np.mean(g1_values)) / pooled_std
            else:
                cohens_d = 0

            # Rank-biserial correlation from U statistic
            n1, n2 = len(g1_values), len(g2_values)
            rank_biserial = (2 * stat) / (n1 * n2) - 1

            # Log2 fold change
            mean_g1 = np.mean(g1_values)
            mean_g2 = np.mean(g2_values)
            log2fc = np.log2((mean_g2 + 0.001) / (mean_g1 + 0.001))

            results.append({
                'reaction': reaction,
                'pvalue': pval,
                'cohens_d': cohens_d,
                'rank_biserial': rank_biserial,
                'log2FC': log2fc,
                f'mean_{group1}': mean_g1,
                f'mean_{group2}': mean_g2
            })

        progress_bar.progress(100)
        status_text.text("Computing adjusted p-values...")

        results_df = pd.DataFrame(results)

        # Multiple testing correction
        results_df['padj'] = multipletests(results_df['pvalue'], method='fdr_bh')[1]

        # Sort by adjusted p-value
        results_df = results_df.sort_values('padj')

        st.session_state.compass_results = results_df
        st.session_state.compass_group1 = group1
        st.session_state.compass_group2 = group2

        status_text.text("Analysis complete!")
        st.success("Analysis complete")

    # ========================================
    # Step 4: Results
    # ========================================
    if 'compass_results' in st.session_state:
        st.header("Step 4: Results")

        results_df = st.session_state.compass_results
        group1 = st.session_state.compass_group1
        group2 = st.session_state.compass_group2

        # Add dynamic significance flags based on selected effect size
        results_df['effect_size'] = results_df[effect_col]
        results_df['significant'] = (results_df['padj'] < padj_threshold) & (np.abs(results_df['effect_size']) > effect_threshold)
        results_df['direction'] = np.where(results_df['effect_size'] > 0, f'Up in {group2}', f'Down in {group2}')

        # Summary
        sig_up = ((results_df['padj'] < padj_threshold) & (results_df['effect_size'] > effect_threshold)).sum()
        sig_down = ((results_df['padj'] < padj_threshold) & (results_df['effect_size'] < -effect_threshold)).sum()

        col1, col2, col3, col4 = st.columns(4)
        with col1:
            st.metric("Total reactions", len(results_df))
        with col2:
            st.metric(f"Significant (padj<{padj_threshold})", (results_df['padj'] < padj_threshold).sum())
        with col3:
            st.metric(f"Up in {group2}", sig_up)
        with col4:
            st.metric(f"Down in {group2}", sig_down)

        # ========================================
        # Prepare subsystem annotation (shared across tabs)
        # ========================================
        subsystem_annot = get_subsystem_mapping(_cache_key="v4_with_names")

        results_with_subsystem = None
        if subsystem_annot is not None:
            results_with_subsystem = results_df.copy()
            results_with_subsystem['mar_id'] = results_with_subsystem['reaction'].str.replace('_pos|_neg', '', regex=True)
            results_with_subsystem = results_with_subsystem.merge(
                subsystem_annot,
                left_on='mar_id',
                right_on='reaction_id',
                how='left'
            )
            results_with_subsystem['subsystem'] = results_with_subsystem['subsystem'].fillna('Unknown')
            st.session_state.compass_results_annotated = results_with_subsystem

        # ========================================
        # Fragment functions for each tab
        # ========================================
        @st.fragment
        def render_volcano_tab():
            from matplotlib.patches import Patch
            fig, ax = plt.subplots(figsize=(10, 8))
            v_colors = np.where(
                (results_df['padj'] < padj_threshold) & (np.abs(results_df['effect_size']) > effect_threshold),
                np.where(results_df['effect_size'] > 0, 'red', 'blue'),
                'gray'
            )
            ax.scatter(results_df['effect_size'], -np.log10(results_df['padj']),
                       c=v_colors, alpha=0.5, s=10)
            ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--', alpha=0.5)
            ax.axvline(effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(-effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel(f"{effect_label} ({group2} - {group1})", fontsize=12)
            ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
            ax.set_title(f"COMPASS Differential Analysis: {group2} vs {group1}", fontsize=14)
            legend_elements = [
                Patch(facecolor='red', alpha=0.5, label=f'Up in {group2}'),
                Patch(facecolor='blue', alpha=0.5, label=f'Down in {group2}'),
                Patch(facecolor='gray', alpha=0.5, label='Not significant')
            ]
            ax.legend(handles=legend_elements, loc='lower right')
            plt.tight_layout()
            st.pyplot(fig)
            # Download buttons
            dl1, dl2 = st.columns(2)
            with dl1:
                buf_png = BytesIO()
                fig.savefig(buf_png, format='png', dpi=150, bbox_inches='tight')
                buf_png.seek(0)
                st.download_button("Volcano Plot (PNG)", data=buf_png.getvalue(),
                                   file_name=f"compass_volcano_{group2}_vs_{group1}.png",
                                   mime="image/png", key="volcano_main_png")
            with dl2:
                buf_pdf = BytesIO()
                fig.savefig(buf_pdf, format='pdf', bbox_inches='tight')
                buf_pdf.seek(0)
                st.download_button("Volcano Plot (PDF)", data=buf_pdf.getvalue(),
                                   file_name=f"compass_volcano_{group2}_vs_{group1}.pdf",
                                   mime="application/pdf", key="volcano_main_pdf")
            plt.close()

        @st.fragment
        def render_top_reactions_tab():
            st.subheader("Top Differential Reactions")
            top_source = results_with_subsystem if results_with_subsystem is not None else results_df
            top_display_cols = ['reaction', 'effect_size', 'log2FC', 'pvalue', 'padj']
            if results_with_subsystem is not None:
                top_display_cols = ['reaction', 'reaction_name', 'subsystem', 'effect_size', 'log2FC', 'pvalue', 'padj']
                # Ensure effect_size column exists in annotated df
                if 'effect_size' not in top_source.columns:
                    top_source = top_source.copy()
                    top_source['effect_size'] = top_source[effect_col]

            sub_tab1, sub_tab2 = st.tabs([f"Up in {group2}", f"Down in {group2}"])
            with sub_tab1:
                up_df = top_source[top_source['effect_size'] > 0].head(50).copy()
                up_df['effect_size'] = up_df['effect_size'].apply(lambda x: f"{x:.3f}")
                up_df['log2FC'] = up_df['log2FC'].apply(lambda x: f"{x:.3f}")
                up_df['pvalue'] = up_df['pvalue'].apply(lambda x: f"{x:.2e}")
                up_df['padj'] = up_df['padj'].apply(lambda x: f"{x:.2e}")
                show_cols = [c for c in top_display_cols if c in up_df.columns]
                st.dataframe(up_df[show_cols].rename(columns={'effect_size': effect_label}))
            with sub_tab2:
                down_df = top_source[top_source['effect_size'] < 0].head(50).copy()
                down_df['effect_size'] = down_df['effect_size'].apply(lambda x: f"{x:.3f}")
                down_df['log2FC'] = down_df['log2FC'].apply(lambda x: f"{x:.3f}")
                down_df['pvalue'] = down_df['pvalue'].apply(lambda x: f"{x:.2e}")
                down_df['padj'] = down_df['padj'].apply(lambda x: f"{x:.2e}")
                show_cols = [c for c in top_display_cols if c in down_df.columns]
                st.dataframe(down_df[show_cols].rename(columns={'effect_size': effect_label}))

        @st.fragment
        def render_bar_plot_tab():
            if results_with_subsystem is None:
                st.warning("Subsystem annotation file not found")
                return
            st.subheader("Up/Down-regulated Reactions per Subsystem")
            sig_annotated = results_with_subsystem[
                (results_with_subsystem['padj'] < padj_threshold) &
                (np.abs(results_with_subsystem['effect_size']) > effect_threshold)
            ].copy()
            if len(sig_annotated) == 0:
                st.info("No significant reactions found")
                return
            sig_annotated['direction'] = np.where(sig_annotated['effect_size'] > 0, 'Up', 'Down')
            subsystem_direction = sig_annotated.groupby(['subsystem', 'direction']).size().unstack(fill_value=0)
            subsystem_direction['total'] = subsystem_direction.sum(axis=1)
            subsystem_direction = subsystem_direction.sort_values('total', ascending=True)
            top_n = st.slider("Number of subsystems to display", min_value=10, max_value=50, value=20, key="bar_top_n")
            plot_data = subsystem_direction.tail(top_n).drop(columns=['total'])
            fig, ax = plt.subplots(figsize=(10, max(6, top_n * 0.3)))
            y_pos = np.arange(len(plot_data))
            bar_height = 0.8
            if 'Up' in plot_data.columns:
                ax.barh(y_pos, plot_data['Up'], bar_height, label=f'Up in {group2}', color='red', alpha=0.7)
            if 'Down' in plot_data.columns:
                ax.barh(y_pos, -plot_data['Down'], bar_height, label=f'Down in {group2}', color='blue', alpha=0.7)
            ax.set_yticks(y_pos)
            ax.set_yticklabels(plot_data.index, fontsize=9)
            ax.set_xlabel('Number of Reactions')
            ax.set_title(f'Differential Reactions per Subsystem\n{group2} vs {group1}')
            ax.axvline(0, color='black', linewidth=0.5)
            ax.legend(loc='lower right')
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()
            with st.expander("Subsystem summary table"):
                display_df = subsystem_direction.drop(columns=['total']).sort_values(
                    by=['Up', 'Down'] if 'Up' in subsystem_direction.columns else ['Down'],
                    ascending=False
                )
                st.dataframe(display_df)

        @st.fragment
        def render_dot_plot_tab():
            if results_with_subsystem is None:
                st.warning("Subsystem annotation file not found")
                return
            st.subheader("Up- and Down-Regulated Reactions Cross Pathway Boundaries")
            st.caption(f"Plot each reaction's {effect_label} by subsystem (compassR style)")
            min_reactions_for_plot = st.number_input(
                "Minimum number of reactions (subsystem display condition)",
                min_value=1, max_value=50, value=5,
                help="Only display subsystems containing at least this many reactions",
                key="dot_min_reactions"
            )
            subsystem_counts_for_dot = results_with_subsystem['subsystem'].value_counts()
            subsystems_to_plot = subsystem_counts_for_dot[subsystem_counts_for_dot >= min_reactions_for_plot].index.tolist()
            if len(subsystems_to_plot) == 0:
                st.warning(f"No subsystems meet the minimum reaction count of {min_reactions_for_plot}")
                return
            dot_plot_data = results_with_subsystem[
                results_with_subsystem['subsystem'].isin(subsystems_to_plot)
            ].copy()
            subsystem_mean_d = dot_plot_data.groupby('subsystem')['effect_size'].mean().sort_values()
            sorted_subsystems = subsystem_mean_d.index.tolist()
            max_subsystems = st.slider(
                "Number of subsystems to display",
                min_value=10, max_value=min(100, len(sorted_subsystems)),
                value=min(50, len(sorted_subsystems)),
                key="cross_pathway_max_subsystems"
            )
            sorted_subsystems = sorted_subsystems[:max_subsystems]
            dot_plot_data = dot_plot_data[dot_plot_data['subsystem'].isin(sorted_subsystems)]
            fig, ax = plt.subplots(figsize=(12, max(8, len(sorted_subsystems) * 0.25)))
            subsystem_to_y = {s: i for i, s in enumerate(sorted_subsystems)}
            dot_plot_data['y_pos'] = dot_plot_data['subsystem'].map(subsystem_to_y)
            dot_plot_data['dot_sig'] = (
                (dot_plot_data['padj'] < padj_threshold) &
                (np.abs(dot_plot_data['effect_size']) > effect_threshold)
            )
            d_colors = []
            d_alphas = []
            d_sizes = []
            for _, row in dot_plot_data.iterrows():
                d_colors.append('red' if row['effect_size'] > 0 else 'blue')
                if row['dot_sig']:
                    d_alphas.append(0.6)
                    d_sizes.append(40)
                else:
                    d_alphas.append(0.2)
                    d_sizes.append(20)
            for i, (_, row) in enumerate(dot_plot_data.iterrows()):
                ax.scatter(
                    row['effect_size'], row['y_pos'],
                    c=d_colors[i], alpha=d_alphas[i], s=d_sizes[i],
                    edgecolors=d_colors[i], linewidths=0.3
                )
            for y in range(len(sorted_subsystems)):
                ax.axhline(y, color='gray', linewidth=0.3, alpha=0.3)
            ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
            ax.set_yticks(range(len(sorted_subsystems)))
            ax.set_yticklabels(sorted_subsystems, fontsize=8)
            ax.set_xlabel(effect_label, fontsize=12)
            ax.set_title(f"Up- and Down-Regulated Reactions Cross Pathway Boundaries\n{group2} vs {group1}", fontsize=12)
            from matplotlib.lines import Line2D
            legend_elements = [
                Line2D([0], [0], marker='o', color='w', markerfacecolor='gray',
                       markersize=6, alpha=0.3, label='insignificant'),
                Line2D([0], [0], marker='o', color='w', markerfacecolor='black',
                       markersize=8, alpha=0.8, label=f'BH-adjusted p-value < {padj_threshold}')
            ]
            ax.legend(handles=legend_elements, loc='lower center', bbox_to_anchor=(0.5, -0.15),
                     ncol=2, frameon=False, fontsize=10)
            plt.tight_layout()
            plt.subplots_adjust(bottom=0.15)
            st.pyplot(fig)
            plt.close()
            with st.expander("Cross-pathway summary"):
                sig_up_dot = len(dot_plot_data[(dot_plot_data['dot_sig']) & (dot_plot_data['effect_size'] > 0)])
                sig_down_dot = len(dot_plot_data[(dot_plot_data['dot_sig']) & (dot_plot_data['effect_size'] < 0)])
                st.write(f"**Significantly up-regulated:** {sig_up_dot} reactions")
                st.write(f"**Significantly down-regulated:** {sig_down_dot} reactions")

        @st.fragment
        def render_subsystem_volcano_tab():
            if results_with_subsystem is None:
                st.warning("Subsystem annotation file not found")
                sig_results = results_df[(results_df['padj'] < padj_threshold) & (np.abs(results_df['effect_size']) > effect_threshold)]
                if len(sig_results) > 0:
                    sig_results = sig_results.copy()
                    sig_results['mar_prefix'] = sig_results['reaction'].str.extract(r'(MAR\d{3})')
                    prefix_summary = sig_results.groupby(['mar_prefix', 'direction']).size().unstack(fill_value=0)
                    st.dataframe(prefix_summary)
                return
            from matplotlib.patches import Patch as Patch2
            st.subheader("Subsystem-specific Volcano Plot")
            if st.button("Reload subsystem mapping", help="Clear cache and reload"):
                clear_subsystem_cache()
                st.rerun()
            with st.expander("Debug: Subsystem mapping info", expanded=False):
                matched = results_with_subsystem['subsystem'].notna().sum()
                st.write(f"Annotation file: {len(subsystem_annot)} reactions")
                st.write(f"Sample annotation IDs: {subsystem_annot['reaction_id'].head(5).tolist()}")
                st.write(f"Sample result MAR IDs: {results_with_subsystem['mar_id'].head(5).tolist()}")
                st.write(f"Matched: {matched}/{len(results_with_subsystem)}")
                unmatched = results_with_subsystem[results_with_subsystem['subsystem'] == 'Unknown']['mar_id'].head(10).tolist()
                st.write(f"Unmatched examples: {unmatched}")
            subsystem_counts = results_with_subsystem['subsystem'].value_counts()
            st.write(f"**{len(subsystem_counts)} subsystems found**")
            subsystem_reaction_counts = results_with_subsystem['subsystem'].value_counts()
            available_subsystems = subsystem_reaction_counts[subsystem_reaction_counts >= 5].index.tolist()
            selected_subsystem = st.selectbox(
                "Select subsystem",
                options=available_subsystems,
                index=0 if available_subsystems else None,
                help="Only subsystems with 5 or more reactions are shown"
            )
            if selected_subsystem:
                subsystem_data = results_with_subsystem[results_with_subsystem['subsystem'] == selected_subsystem]
                fig, ax = plt.subplots(figsize=(10, 8))
                sv_colors = np.where(
                    (subsystem_data['padj'] < padj_threshold) & (np.abs(subsystem_data['effect_size']) > effect_threshold),
                    np.where(subsystem_data['effect_size'] > 0, 'red', 'blue'),
                    'gray'
                )
                ax.scatter(subsystem_data['effect_size'], -np.log10(subsystem_data['padj']),
                           c=sv_colors, alpha=0.6, s=30)
                ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--', alpha=0.5)
                ax.axvline(effect_threshold, color='gray', linestyle='--', alpha=0.5)
                ax.axvline(-effect_threshold, color='gray', linestyle='--', alpha=0.5)
                ax.set_xlabel(f"{effect_label} ({group2} - {group1})", fontsize=12)
                ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
                ax.set_title(f"{selected_subsystem}\n({len(subsystem_data)} reactions)", fontsize=14)
                top_sig = subsystem_data[
                    (subsystem_data['padj'] < padj_threshold) &
                    (np.abs(subsystem_data['effect_size']) > effect_threshold)
                ].nsmallest(5, 'padj')
                for _, row in top_sig.iterrows():
                    label = row['reaction_name'] if pd.notna(row.get('reaction_name')) else row['reaction']
                    label = label[:30] + '...' if len(str(label)) > 30 else label
                    ax.annotate(label, (row['effect_size'], -np.log10(row['padj'])),
                               fontsize=8, alpha=0.8)
                legend_elements = [
                    Patch2(facecolor='red', alpha=0.6, label=f'Up in {group2}'),
                    Patch2(facecolor='blue', alpha=0.6, label=f'Down in {group2}'),
                    Patch2(facecolor='gray', alpha=0.6, label='Not significant')
                ]
                ax.legend(handles=legend_elements, loc='lower right')
                plt.tight_layout()
                # Save figure to buffers before closing
                fig_buf = BytesIO()
                fig.savefig(fig_buf, format='png', dpi=150, bbox_inches='tight')
                fig_buf.seek(0)
                fig_buf_pdf = BytesIO()
                fig.savefig(fig_buf_pdf, format='pdf', bbox_inches='tight')
                fig_buf_pdf.seek(0)
                st.pyplot(fig)
                plt.close()

                with st.expander(f"Reaction list for {selected_subsystem}"):
                    display_cols = ['reaction', 'reaction_name', 'effect_size', 'log2FC', 'padj']
                    display_cols = [c for c in display_cols if c in subsystem_data.columns]
                    st.dataframe(subsystem_data[display_cols].rename(columns={'effect_size': effect_label}).sort_values('padj').head(50))

                # Download for this subsystem
                safe_name = selected_subsystem.replace('/', '_').replace(' ', '_')
                dl_col1, dl_col2, dl_col3 = st.columns(3)
                with dl_col1:
                    sub_csv = BytesIO()
                    sub_dl_cols = ['reaction', 'reaction_name', 'effect_size', 'log2FC', 'pvalue', 'padj']
                    sub_dl_cols = [c for c in sub_dl_cols if c in subsystem_data.columns]
                    subsystem_data[sub_dl_cols].sort_values('padj').to_csv(sub_csv, sep='\t', index=False)
                    sub_csv.seek(0)
                    st.download_button(
                        label=f"TSV: {selected_subsystem}",
                        data=sub_csv.getvalue(),
                        file_name=f"compass_{safe_name}_{group2}_vs_{group1}.tsv",
                        mime="text/tab-separated-values",
                        key="dl_subsystem_tsv"
                    )
                with dl_col2:
                    st.download_button(
                        label=f"PNG: {selected_subsystem}",
                        data=fig_buf.getvalue(),
                        file_name=f"compass_{safe_name}_{group2}_vs_{group1}.png",
                        mime="image/png",
                        key="dl_subsystem_png"
                    )
                with dl_col3:
                    st.download_button(
                        label=f"PDF: {selected_subsystem}",
                        data=fig_buf_pdf.getvalue(),
                        file_name=f"compass_{safe_name}_{group2}_vs_{group1}.pdf",
                        mime="application/pdf",
                        key="dl_subsystem_pdf"
                    )

        @st.fragment
        def render_download_tab():
            st.subheader("Download")
            col1, col2 = st.columns(2)
            with col1:
                csv_buffer = BytesIO()
                results_df.to_csv(csv_buffer, sep='\t', index=False)
                csv_buffer.seek(0)
                st.download_button(
                    label="All results (TSV)",
                    data=csv_buffer.getvalue(),
                    file_name=f"compass_differential_{group1}_vs_{group2}.tsv",
                    mime="text/tab-separated-values"
                )
            with col2:
                sig_csv = BytesIO()
                sig_only = results_df[results_df['significant']]
                sig_only.to_csv(sig_csv, sep='\t', index=False)
                sig_csv.seek(0)
                st.download_button(
                    label="Significant results only (TSV)",
                    data=sig_csv.getvalue(),
                    file_name=f"compass_significant_{group1}_vs_{group2}.tsv",
                    mime="text/tab-separated-values"
                )
            from matplotlib.patches import Patch as Patch3
            fig, ax = plt.subplots(figsize=(10, 8))
            dl_colors = np.where(
                (results_df['padj'] < padj_threshold) & (np.abs(results_df['effect_size']) > effect_threshold),
                np.where(results_df['effect_size'] > 0, 'red', 'blue'),
                'gray'
            )
            ax.scatter(results_df['effect_size'], -np.log10(results_df['padj']),
                       c=dl_colors, alpha=0.5, s=10)
            ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--', alpha=0.5)
            ax.axvline(effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(-effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel(f"{effect_label} ({group2} - {group1})", fontsize=12)
            ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
            ax.set_title(f"COMPASS Differential Analysis: {group2} vs {group1}", fontsize=14)
            legend_elements = [
                Patch3(facecolor='red', alpha=0.5, label=f'Up in {group2}'),
                Patch3(facecolor='blue', alpha=0.5, label=f'Down in {group2}'),
                Patch3(facecolor='gray', alpha=0.5, label='Not significant')
            ]
            ax.legend(handles=legend_elements, loc='lower right')
            plt.tight_layout()
            fig_buffer = BytesIO()
            fig.savefig(fig_buffer, format='png', dpi=150, bbox_inches='tight')
            fig_buffer.seek(0)
            fig_buffer_pdf = BytesIO()
            fig.savefig(fig_buffer_pdf, format='pdf', bbox_inches='tight')
            fig_buffer_pdf.seek(0)
            plt.close()
            dl_col1, dl_col2 = st.columns(2)
            with dl_col1:
                st.download_button(
                    label="Volcano Plot (PNG)",
                    data=fig_buffer.getvalue(),
                    file_name=f"compass_volcano_{group1}_vs_{group2}.png",
                    mime="image/png"
                )
            with dl_col2:
                st.download_button(
                    label="Volcano Plot (PDF)",
                    data=fig_buffer_pdf.getvalue(),
                    file_name=f"compass_volcano_{group1}_vs_{group2}.pdf",
                    mime="application/pdf"
                )

        @st.fragment
        def render_pathway_tab():
            from matplotlib.patches import Patch as PatchPW
            if subsystem_annot is None:
                st.warning("Subsystem annotation file not found")
                return

            st.subheader("Pathway-level Differential Analysis")
            st.caption("Aggregate reaction scores per cell within each subsystem, then perform Wilcoxon test at the pathway level")

            pw_col1, pw_col2, pw_col3 = st.columns(3)
            with pw_col1:
                agg_method = st.radio("Aggregation method", ["mean", "median"], horizontal=True, key="compass_pw_agg")
            with pw_col2:
                min_rxns = st.number_input("Minimum number of reactions", min_value=1, max_value=20, value=3, key="compass_pw_min_rxns")
            with pw_col3:
                split_direction = st.checkbox("Split by reaction direction", value=False, key="compass_pw_split",
                    help="Split pathways by _pos (forward) / _neg (reverse) reaction suffix")

            agg_func = np.mean if agg_method == "mean" else np.median

            # Map reactions to subsystems
            rxn_to_sub = dict(zip(subsystem_annot['reaction_id'], subsystem_annot['subsystem']))
            rxn_series = reactions.index.to_series()
            mar_ids = rxn_series.str.replace(r'_pos|_neg', '', regex=True)
            rxn_subsystems = mar_ids.map(rxn_to_sub)

            # Determine per-reaction direction from _pos/_neg suffix
            if split_direction:
                rxn_pos_neg = {}
                for rxn in reactions.index:
                    if rxn.endswith('_pos'):
                        rxn_pos_neg[rxn] = 'pos'
                    elif rxn.endswith('_neg'):
                        rxn_pos_neg[rxn] = 'neg'
                    else:
                        rxn_pos_neg[rxn] = 'pos'  # default
                n_pos = sum(1 for v in rxn_pos_neg.values() if v == 'pos')
                n_neg = sum(1 for v in rxn_pos_neg.values() if v == 'neg')
                st.caption(f"Reaction direction: _pos={n_pos}, _neg={n_neg} reactions")

            # Aggregate per subsystem per cell
            unique_subs = rxn_subsystems.dropna().unique()
            pathway_scores = {}
            n_rxns = {}
            if split_direction:
                for sub in unique_subs:
                    sub_rxn_idx = reactions.index[rxn_subsystems == sub]
                    pos_rxns = [r for r in sub_rxn_idx if rxn_pos_neg.get(r) == 'pos']
                    neg_rxns = [r for r in sub_rxn_idx if rxn_pos_neg.get(r) == 'neg']
                    if len(pos_rxns) >= min_rxns:
                        pathway_scores[f"{sub} (forward)"] = agg_func(reactions.loc[pos_rxns].values, axis=0)
                        n_rxns[f"{sub} (forward)"] = len(pos_rxns)
                    if len(neg_rxns) >= min_rxns:
                        pathway_scores[f"{sub} (reverse)"] = agg_func(reactions.loc[neg_rxns].values, axis=0)
                        n_rxns[f"{sub} (reverse)"] = len(neg_rxns)
            else:
                for sub in unique_subs:
                    sub_rxn_idx = reactions.index[rxn_subsystems == sub]
                    if len(sub_rxn_idx) >= min_rxns:
                        pathway_scores[sub] = agg_func(reactions.loc[sub_rxn_idx].values, axis=0)
                        n_rxns[sub] = len(sub_rxn_idx)

            if len(pathway_scores) == 0:
                st.warning("No pathways meet the specified criteria")
                return

            pathway_df = pd.DataFrame(pathway_scores, index=reactions.columns)

            # Wilcoxon + effect size per pathway
            pw_results = []
            for pathway in pathway_df.columns:
                g1_v = pathway_df.loc[g1_cells, pathway].values.astype(float)
                g2_v = pathway_df.loc[g2_cells, pathway].values.astype(float)

                if np.std(g1_v) == 0 and np.std(g2_v) == 0:
                    continue
                try:
                    stat, pval = mannwhitneyu(g2_v, g1_v, alternative='two-sided')
                except:
                    continue

                pooled_std = np.sqrt(
                    ((len(g2_v)-1)*np.std(g2_v, ddof=1)**2 +
                     (len(g1_v)-1)*np.std(g1_v, ddof=1)**2) /
                    (len(g2_v) + len(g1_v) - 2)
                )
                cd = (np.mean(g2_v) - np.mean(g1_v)) / pooled_std if pooled_std > 0 else 0
                n1, n2 = len(g1_v), len(g2_v)
                rb = (2 * stat) / (n1 * n2) - 1
                es = rb if is_rank_biserial else cd

                pw_results.append({
                    'pathway': pathway,
                    'n_reactions': n_rxns[pathway],
                    'pvalue': pval,
                    'effect_size': es,
                    'log2FC': np.log2((np.mean(g2_v) + 0.001) / (np.mean(g1_v) + 0.001)),
                    f'mean_{group1}': np.mean(g1_v),
                    f'mean_{group2}': np.mean(g2_v),
                })

            if len(pw_results) == 0:
                st.warning("No analysis results")
                return

            pw_df = pd.DataFrame(pw_results)
            pw_df['padj'] = multipletests(pw_df['pvalue'], method='fdr_bh')[1]
            pw_df = pw_df.sort_values('padj')
            pw_df['significant'] = (pw_df['padj'] < padj_threshold) & (np.abs(pw_df['effect_size']) > effect_threshold)

            # Summary
            sig_up_pw = ((pw_df['padj'] < padj_threshold) & (pw_df['effect_size'] > effect_threshold)).sum()
            sig_down_pw = ((pw_df['padj'] < padj_threshold) & (pw_df['effect_size'] < -effect_threshold)).sum()

            c1, c2, c3 = st.columns(3)
            with c1:
                st.metric("Total pathways", len(pw_df))
            with c2:
                st.metric(f"Up in {group2}", sig_up_pw)
            with c3:
                st.metric(f"Down in {group2}", sig_down_pw)

            # Volcano plot
            fig, ax = plt.subplots(figsize=(10, 8))
            pw_colors = np.where(
                pw_df['significant'],
                np.where(pw_df['effect_size'] > 0, 'red', 'blue'),
                'gray'
            )
            ax.scatter(pw_df['effect_size'], -np.log10(pw_df['padj']),
                       c=pw_colors, alpha=0.6, s=50)
            ax.axhline(-np.log10(padj_threshold), color='red', linestyle='--', alpha=0.5)
            ax.axvline(effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.axvline(-effect_threshold, color='gray', linestyle='--', alpha=0.5)
            ax.set_xlabel(f"{effect_label} ({group2} - {group1})", fontsize=12)
            ax.set_ylabel("-log10(adjusted p-value)", fontsize=12)
            ax.set_title(f"Pathway-level Analysis ({agg_method})\n{group2} vs {group1}", fontsize=14)

            for _, row in pw_df[pw_df['significant']].iterrows():
                label = row['pathway'][:35] + '...' if len(row['pathway']) > 35 else row['pathway']
                ax.annotate(label, (row['effect_size'], -np.log10(row['padj'])),
                           fontsize=7, alpha=0.8)

            legend_elements = [
                PatchPW(facecolor='red', alpha=0.6, label=f'Up in {group2}'),
                PatchPW(facecolor='blue', alpha=0.6, label=f'Down in {group2}'),
                PatchPW(facecolor='gray', alpha=0.6, label='Not significant')
            ]
            ax.legend(handles=legend_elements, loc='lower right')
            plt.tight_layout()
            st.pyplot(fig)
            # Download buttons for pathway volcano
            pw_dl1, pw_dl2 = st.columns(2)
            with pw_dl1:
                pw_png_buf = BytesIO()
                fig.savefig(pw_png_buf, format='png', dpi=150, bbox_inches='tight')
                pw_png_buf.seek(0)
                st.download_button("Pathway Volcano (PNG)", data=pw_png_buf.getvalue(),
                                   file_name=f"compass_pathway_volcano_{group2}_vs_{group1}.png",
                                   mime="image/png", key="dl_pw_volcano_png")
            with pw_dl2:
                pw_pdf_buf = BytesIO()
                fig.savefig(pw_pdf_buf, format='pdf', bbox_inches='tight')
                pw_pdf_buf.seek(0)
                st.download_button("Pathway Volcano (PDF)", data=pw_pdf_buf.getvalue(),
                                   file_name=f"compass_pathway_volcano_{group2}_vs_{group1}.pdf",
                                   mime="application/pdf", key="dl_pw_volcano_pdf")
            plt.close()

            # Results table
            with st.expander("Pathway results table", expanded=True):
                disp = pw_df[['pathway', 'n_reactions', 'effect_size', 'log2FC', 'pvalue', 'padj']].copy()
                disp['effect_size'] = disp['effect_size'].apply(lambda x: f"{x:.3f}")
                disp['log2FC'] = disp['log2FC'].apply(lambda x: f"{x:.3f}")
                disp['pvalue'] = disp['pvalue'].apply(lambda x: f"{x:.2e}")
                disp['padj'] = disp['padj'].apply(lambda x: f"{x:.2e}")
                st.dataframe(disp.rename(columns={'effect_size': effect_label}))

            # Download
            pw_csv = BytesIO()
            pw_df.to_csv(pw_csv, sep='\t', index=False)
            pw_csv.seek(0)
            st.download_button(
                label="Pathway results (TSV)",
                data=pw_csv.getvalue(),
                file_name=f"compass_pathway_{group2}_vs_{group1}.tsv",
                mime="text/tab-separated-values",
                key="dl_compass_pathway"
            )

        # ========================================
        # Tabs - each calls its fragment function
        # ========================================
        tab_names = ["Volcano Plot", "Top Reactions", "Bar Plot", "Dot Plot", "Subsystem Volcano", "Pathway Analysis", "Download"]
        tab_volcano, tab_top, tab_bar, tab_dot, tab_sub_volcano, tab_pathway, tab_download = st.tabs(tab_names)

        with tab_volcano:
            render_volcano_tab()
        with tab_top:
            render_top_reactions_tab()
        with tab_bar:
            render_bar_plot_tab()
        with tab_dot:
            render_dot_plot_tab()
        with tab_sub_volcano:
            render_subsystem_volcano_tab()
        with tab_pathway:
            render_pathway_tab()
        with tab_download:
            render_download_tab()

        st.caption("""
        For detailed pathway information, visit [Metabolic Atlas](https://metabolicatlas.org/).
        """)

else:
    st.info("Upload data files to begin analysis")

# Footer
st.markdown("---")
st.markdown("""
**Reference:** Wagner, A., et al. (2021). Metabolic modeling of single Th17 cells reveals
regulators of autoimmunity. Cell, 184(16), 4168-4185.

**Note:** This tool provides equivalent functionality to compassR (deprecated) using Python.
""")
