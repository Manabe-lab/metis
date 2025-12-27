import streamlit as st
import decoupler as dc
import pandas as pd
import numpy as np
import sys
import os
import re
from jellyfish import jaro_winkler_similarity
import plotly.express as px
from helper_func import clear_old_directories
from helper_func import clear_old_files
from helper_func import check_species_index
import time
import shutil
import matplotlib.pyplot as plt
from io import StringIO

# Custom component temporarily disabled due to module loading issues in dynamic navigation
# import streamlit.components.v1 as components
#
# import os
#
# # Initialize component as a global variable
# if 'pathway_summary' not in st.session_state:
#     # Get the directory of the current script
#     current_dir = os.path.dirname(os.path.abspath(__file__))
#     # Build absolute path to components directory
#     components_dir = os.path.join(os.path.dirname(current_dir), 'components', 'pathway_summary')
#
#     # Define custom component - use components directly
#     st.session_state.pathway_summary = components.declare_component(
#         "pathway_summary",
#         path=components_dir,
#         url=None  # Explicitly set url to None
#     )
#
# # Get from session state
# pathway_summary = st.session_state.pathway_summary

import os

st.set_page_config(page_title="decoupleR_pathway_analysis", page_icon="⛖",
    layout="wide")

st.sidebar.title("Options")


def get_tf_databases(species):
    """Return TF database definitions"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"

    return [
        {'name': 'DoRothEA A', 'path': 'db', 'file': f'dorothea.{species}.tsv', 'filter': {'confidence': ['A']},
         'source': 'tf', 'target': 'target', 'type': 'dorothea'},
        {'name': 'DoRothEA B', 'path': 'db', 'file': f'dorothea.{species}.tsv', 'filter': {'confidence': ['A', 'B']},
         'source': 'tf', 'target': 'target', 'type': 'dorothea'},
 #       {'name': 'MSigDB TFT', 'path': msigdb_dir, 'file': 'c3.tft.v2023.2.Hs.symbols.gmt' if species=='human' else 'c3.tft.v2023.2.Hs.symbols.Mouse.gmt',
 #        'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'CollecTRI', 'path': 'db', 'file': f'TRI.{species}.tsv',
         'source': 'source', 'target': 'target', 'type': 'collecTRI'},
        {'name': 'ChEA 2022', 'path': enrichr_dir, 'file': 'ChEA_2022.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'ENCODE and ChEA', 'path': enrichr_dir, 'file': 'ENCODE_and_ChEA_Consensus_TFs_from_ChIP-X.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'ARCHS4 TF Coexp', 'path': enrichr_dir, 'file': 'ARCHS4_TFs_Coexp.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'TRRUST 2019', 'path': enrichr_dir, 'file': 'TRRUST_Transcription_Factors_2019.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
        {'name': 'TRANSFAC and JASPAR PWMs', 'path': enrichr_dir, 'file': 'TRANSFAC_and_JASPAR_PWMs.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
    ]

def get_pathway_databases(species):
    """Return pathway database definitions"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"

    return [
            {'name': 'MSigDB Hallmark', 'path': msigdb_dir, 'file': 'h.all.v2023.2.Hs.symbols.gmt' if species=='human' else 'mh.all.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB CP', 'path': msigdb_dir, 'file': 'c2.cp.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB GOBP', 'path': msigdb_dir, 'file': 'c5.go.bp.v2023.2.Hs.symbols.gmt' if species=='human' else 'm5.go.bp.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB Reactome', 'path': msigdb_dir, 'file': 'c2.cp.reactome.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.reactome.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr Reactome 2022', 'path': enrichr_dir, 'file': 'Reactome_2022.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr GOBP 2023', 'path': enrichr_dir, 'file': 'GO_Biological_Process_2023.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'MSigDB KEGG legacy', 'path': msigdb_dir, 'file': 'c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt' if species=='human' else 'm2.cp.kegg_legacy.v2023.2.Hs.symbols2Mouse.txt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'Enrichr KEGG 2021', 'path': enrichr_dir, 'file': 'KEGG_2021_Human.gmt' if species=='human' else 'KEGG_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
        ]

def get_celltype_databases(species):
    """Return cell type database definitions"""
    msigdb_dir = "db/mSigDB_mouse" if species == "mouse" else "db/mSigDB"
    enrichr_dir = "db/enrichr_gmt_mouse" if species == "mouse" else "db/enrichr_gmt"

    return [
            {'name': 'MSigDB cell type signature', 'path': msigdb_dir, 'file': 'c8.all.v2023.2.Hs.symbols.gmt' if species=='human' else 'm8.all.v2023.2.Mm.symbols.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'CellMarker 2024', 'path': enrichr_dir, 'file': 'CellMarker_2024.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
            {'name': 'CellMarker augmented 2021', 'path': enrichr_dir, 'file': 'CellMarker_Augmented_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'Tabula Sapiens', 'path': enrichr_dir, 'file': 'Tabula_Sapiens.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'Tabula Muris', 'path': enrichr_dir, 'file': 'Tabula_Muris.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'},
         {'name': 'PanglaoDB augmented', 'path': enrichr_dir, 'file': 'PanglaoDB_Augmented_2021.gmt',
         'source': 'source', 'target': 'target', 'type': 'gmt'}
        ]

def run_summary_analysis(databases, gene_list, n_background):
    """
    Common analysis execution function
    """
    all_results = {}
    progress_bar = st.progress(0)

    # Create 2-column layout
    col1, col2 = st.columns(2)

    # Run analysis for each database
    for idx, db in enumerate(databases):
        try:
            filepath = os.path.join(db['path'], db['file'])
            source = db['source']
            target = db['target']

            if db.get('type') == 'dorothea':
                net = pd.read_csv(filepath, sep='\t')
                net = net[net['confidence'].isin(db['filter']['confidence'])]

            elif db.get('type') == 'CollecTRI':
                net = pd.read_csv(filepath, sep='\t')

            else:
                net = dc.read_gmt(filepath)

            # Run analysis using cached function
            ora_res = run_single_database_analysis(
                filepath,
                gene_list,
                source,
                target,
                n_background,
                db_type=db.get('type'),
                filter_dict=db.get('filter')
            )

            if ora_res is not None:
                # Save results
                all_results[db['name']] = ora_res
                progress_bar.progress((idx + 1) / len(databases))

                # Display plots (split into 2 columns)
                with col1 if idx % 2 == 0 else col2:
                    if len(ora_res) > 0:
                        display_results(db, ora_res)
            else:
                st.warning(f"Database file not found: {filepath}")

        except Exception as e:
            st.error(f"Error processing {db['name']}: {str(e)}")

    progress_bar.empty()
    return all_results

def display_results(db, ora_res):
    """
    Common results display function
    """
    with st.expander(f"📊 {db['name']}", expanded=True):
        # Create plot
        fig, ax = plt.subplots(figsize=(6, 4))

        # Get top 5 and sort by -log10(FDR)
        plot_data = ora_res.head(5).copy()
        plot_data['neg_log10_fdr'] = -np.log10(plot_data['FDR p-value'])
        plot_data = plot_data.sort_values('neg_log10_fdr', ascending=True)

        # Create bar plot and color
        y_pos = range(len(plot_data))
        scores = plot_data['neg_log10_fdr']
        cmap = plt.cm.Reds
        colors = cmap(scores / scores.max())
        bars = ax.barh(y_pos, scores)

        for bar, color in zip(bars, colors):
            bar.set_color(color)

        ax.set_yticks(y_pos)
        ax.set_yticklabels(plot_data['Term'].str[:50])
        ax.set_xlabel('-log10(FDR)')
        plt.tight_layout()
        st.pyplot(fig)
        plt.close()

        # Detailed view button
        if st.button(f"Analyze {db['name']} in detail"):
            display_detailed_analysis(ora_res)

def display_detailed_analysis(ora_res):
    """
    Common detailed analysis display function
    """
    try:
        fig = create_enrichment_plot(ora_res, '')
        st.pyplot(fig)
    except Exception as e:
        st.error(f"Error creating enrichment plot: {str(e)}")

    try:
        fig_dotplot = create_dotplot(ora_res)
        st.pyplot(fig_dotplot)
    except Exception as e:
        st.error(f"Error creating dot plot: {str(e)}")

    # Display results table
    st.markdown("### Results Table")
    # Direct display instead of .style (avoid jinja2 errors)
    display_df = ora_res[['Term', 'FDR p-value', 'Features']].copy()
    display_df['FDR p-value'] = display_df['FDR p-value'].apply(lambda x: f'{x:.2e}')
    st.dataframe(display_df)

@st.cache_data
def run_single_database_analysis(filepath, gene_list, source, target, n_background, db_type=None, filter_dict=None):
    """Cache analysis for a single database"""
    if not os.path.exists(filepath):
        return None

    try:
        if db_type == 'dorothea':
            net = pd.read_csv(filepath, sep='\t')
            if filter_dict:
                for key, values in filter_dict.items():
                    net = net[net[key].isin(values)]

        elif db_type == "collecTRI":
            net = pd.read_csv(filepath, sep='\t')
        else:
            net = dc.read_gmt(filepath)

        ora_res = dc.get_ora_df(gene_list, net, source=source, target=target,
                               n_background=n_background, verbose=False)

        return ora_res.sort_values('FDR p-value', ascending=True)
    except Exception as e:
        st.error(f"Error in run_single_database_analysis: {str(e)}")
        return None

@st.cache_data
def create_enrichment_plot(ora_res_detail, title):
    """Generate and cache enrichment plot"""
    enr_pvals = ora_res_detail[['FDR p-value']]
    enr_pvals.index = ora_res_detail['Term']
    enr_pvals.values[enr_pvals.values == 0] = np.min(enr_pvals.values[enr_pvals.values != 0])
    enr_pvals = -np.log10(enr_pvals)

    g = dc.plot_barplot(enr_pvals.T, 'FDR p-value', vertical=True, top=15,
                       figsize=(8, 6), return_fig=True)
    g.gca().invert_yaxis()
    return g

@st.cache_data
def create_dotplot(ora_res_detail):
    """Generate and cache dot plot"""
    ora_res_detail = ora_res_detail.copy()
    ora_res_detail['count'] = ora_res_detail['Features'].str.split(';').str.len()

    max_count = ora_res_detail['count'].max()
    scale = 4 / max_count if max_count <= 10 else (2 / max_count if max_count <= 50 else 1 / max_count)

    dotplot = dc.plot_dotplot(
        ora_res_detail.sort_values('Combined score', ascending=False).head(15),
        x='Combined score',
        y='Term',
        s='count',
        c='FDR p-value',
        scale=scale,
        figsize=(8, 6),
        return_fig=True
    )
    return dotplot

@st.cache_data
def run_pathway_analysis(gene_list, net, source, target, n_background):
    """Run and cache pathway analysis"""
    ora_res_detail = dc.get_ora_df(gene_list, net, source=source, target=target,
                                  n_background=n_background, verbose=False)
    return ora_res_detail.sort_values('FDR p-value', ascending=True)


# Function to set plot style
def set_plot_style():
    plt.style.use('default')  # Reset default style
    plt.rcParams['figure.facecolor'] = 'white'  # Set figure background color to white
    plt.rcParams['axes.facecolor'] = 'white'    # Set plot area background color to white
#    sns.set_style("white")                      # Set seaborn style to white background

# Call before graph generation
set_plot_style()


def delete_file(file_list):
    for i in file_list:
        if os.path.exists(i):
            os.remove(i)

@st.cache_data
def read_xl(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header=0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep, engine='python')
    return df_c


#@st.cache_data
#def read_csv(file, index_col=None, sep=',', header = 0):
#    df_xl = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
#    return df_xl

@st.cache_data
def run_method(method, mat, net, source, target, weight, verbose=True, min_n = 0):
    norm = None
    corr = None
    if method == 'ulm':
        score, pvalue = dc.run_ulm(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'consensus':
        score, pvalue = dc.run_consensus(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'mlm':
        score, pvalue = dc.run_mlm(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    elif method == 'wsum_norm':
        score, norm, corr, pvalue = dc.run_wsum(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
        score = norm # Use wsum z-score
    elif method == 'viper':
        score, pvalue = dc.run_viper(mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
    if method == 'wsum_norm':
        return score, pvalue, norm, corr
    else:
        return score, pvalue

@st.cache_data
def run_GSEA_df(mat, stat='stat', net = 'net', source='source', target='target', min_n=5, seed=42):
    GSEA_res = dc.get_gsea_df(mat, stat=stat, net = net, source=source, target=target, min_n=min_n, seed=seed)
    return GSEA_res

@st.cache_data
def calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch):
    orig_len = len(df)
    df = df[np.isfinite(df[P_column]) & pd.notnull(df[P_column])]     # Remove rows where FC or p is NA
    df = df[np.isfinite(df[FC_column]) & pd.notnull(df[FC_column])]    # Remove rows where FC or p is NA
    if len(df) < orig_len:
        st.warning("The P or FC columns contain inf or NA")
    inv_parameter = 1
    if inv_switch:
        inv_parameter = -1
    # Check if p=0
    p_0 = (df.loc[:,P_column] == 0)
    if not any(p_0):
        # Create score
        if rank_metric == 'sign(LFC) x -log10(P)':
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
        else:
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
     # When p=0
    else:
        st.write("p=0 data are:")
        st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])
        # Read as 0e0, LogFC should also be 0
        # Identify genes with FC=0 and p=0
        problematic_mask = (df[FC_column] == 0) & (df[P_column] == 0)
        if any(problematic_mask):
            st.warning(f"Found {sum(problematic_mask)} genes with FC=0 and p=0. These will be excluded from analysis.")
            excluded_genes = df.loc[problematic_mask, Gene_column].tolist()
            st.write("Excluded genes:", ", ".join(excluded_genes[:10]), "..." if len(excluded_genes) > 10 else "")
            df = df[~problematic_mask]
            p_0 = (df.loc[:,P_column] == 0) # 0 with FC>0
            if any(p_0):
                st.write("Remaining p=0 data are:")
                st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])

        if rank_metric == 'sign(LFC) x -log10(P)':
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
        else:
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
        # Seurat "MAST" around 318?
        if input_file_type == 'Seurat':
            # max_score = np.log10(1e-324) # 1e-324 == 0 is TRUE, calculating log10 gives inf
            max_score = -324
            st.write("\nMax score: "+str(max_score))
        else:
            # max_score = np.log10(1e-324) # 1e-324 == 0 is TRUE, same in python, 1e-324 + 1e-323 is also calculated
            max_score = -324
            st.write("\nMax score: "+str(max_score))
        # Add FC value for ranking
        df.loc[(p_0 & (df.loc[:,FC_column]>0)),'score'] = max_score * -1 + df.loc[:,FC_column]  * inv_parameter  # Must enclose conditions in parentheses!!!
        df.loc[(p_0 & (df.loc[:,FC_column]<0)),'score'] = max_score + df.loc[:,FC_column] * inv_parameter
        st.write('Ranking score are -log10(P-values)')
    return df['score'].to_frame() # Convert to DF before returning



def set_back_func():
    input_file_type = st.radio(
        "Data format of the file containing all gene names (e.g., DESeq2 output):",
        ('tsv','csv', 'excel'))
    uploaded_file = st.file_uploader("Upload a file containing gene names (e.g., gene list, DESeq2, Homer)", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, header = None, index_col = None)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', header=None, index_col = None)
        else:
            df = read_excel(uploaded_file, index_col = None, header = None)

        # If data is 1 column and doesn't start with Gene
        if df.shape[1] == 1:
            bk_genes = df.iloc[:,1].values
            if bk_genes[0] == "Gene" or bk_genes[0] == "GENE":
                bk_genes = bk_genes[1:]

        else:
            df.columns = df.iloc[0,:].tolist()  # Determine columns after transpose to avoid issues
            df = df.drop(0, axis = 0) # Use first row as column names and remove it

            st.write(df.head())
            content = df.columns.tolist()
            Gene_column = content[0]
            if "Annotation/Divergence" in content:
                  # Convert column names
                search_word = '([^\ \(]*)\ \(.*'

                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content # Temporarily change names
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel compatible

                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
        #       df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
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
                Gene_column =  st.selectbox(
                'Select gene column',
                content)

            bk_genes = df[Gene_column].values
        n_background = len(set(df[Gene_column].to_list()))
        st.write('Background gene number: ' + str(n_background))
        return n_background

#####
# --- Initialising SessionState ---
if "decouplercalc" not in st.session_state:
      st.session_state.decouplercalc = False
if "ORA" not in st.session_state:
      st.session_state.ORA = False

# Save in temp directory
# --- Initialising SessionState ---
if "dc_temp_dir" not in st.session_state:
    st.session_state.dc_temp_dir = True
    dc_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
    os.mkdir(dc_temp_dir)
    st.session_state.dc_temp_dir = dc_temp_dir
else:
    dc_temp_dir = st.session_state.dc_temp_dir

# Control progeny loading

if "net" not in st.session_state:
    st.session_state.net = None

if "num_progeny" not in st.session_state:
    st.session_state.num_progeny = "500"

if 'up_list' not in st.session_state:
    st.session_state.up_list = None
if 'down_list' not in st.session_state:
    st.session_state.down_list = None
if 'P_column' not in st.session_state:
    st.session_state.P_column = None
if 'rank_source' not in st.session_state:
    st.session_state.rank_source = None
if 'Loading_column' not in st.session_state:
    st.session_state.Loading_column = None

#############################


st.markdown('''## Signal pathway and TF activity inference using decoupleR
##### For rank mode, upload DEG result file, PCA loading result file, or rank file
###### e.g.) FD-unshrunk DESeq2 file or PCA loadings file
###### decoupler originally uses stat values from DESeq2 or t value from limma. PCA loadings can be used to identify pathways associated with specific PCs.
---
''')

use_upload = "No"

Analysis_mode = st.radio(
    "##### Analysis mode:",
    ('Rank','Over-representation'), key='Rank')
st.markdown("---")


if Analysis_mode == "Rank":
    input_file_type = st.radio(
    "Data format of DESeq2, PCA loading result or rank file:",
    ('tsv','csv', 'excel', 'rank'), index = 0)
    st.write("DESeq2/PCA loading results: tsv, csv, excel; rank file: rank")

    rnk_input = False

    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx','rnk','rank'])
    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file, index_col = 0)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t', index_col = 0)
        elif input_file_type == 'rank':
            rnk_input = True
            df = read_csv(uploaded_file, sep = '\t', index_col = 0, header = None)
        else:
            df = read_xl(uploaded_file, index_col = 0)

        down_file_name = os.path.basename(uploaded_file.name)

        df.iloc[0:3,:]

else: # ORA
    ORA_mode = st.radio(
    "##### How to provide gene list:",
    ('DEG/PCA loading file','Cluster info file', "Input genes"), key='DEG result file')
    n_background = None
    set_back = False

    if ORA_mode == 'Cluster info file':
        cluster_file_type = st.radio(
            "Data from:",
            ('auto', 'tsv', 'csv', 'excel'))
        cluster_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
        if cluster_file is not None:
            if cluster_file_type == 'auto':
                try:
                    df_cluster = read_csv(cluster_file, sep = None)
                except:# excel
                    df_cluster = read_excel(cluster_file)
            if cluster_file_type == "csv":
                df_cluster = read_csv(cluster_file)
            elif cluster_file_type == 'tsv':
                df_cluster = read_csv(cluster_file, sep = '\t')
            elif cluster_file_type == 'excel':
                df_cluster = read_xl(cluster_file)


            st.write("Preview of uploaded data:")
            st.dataframe(df_cluster.head())


            if "cluster_submitted" not in st.session_state:
            	st.session_state.cluster_submitted = False


            with st.form("Set column"):
                # Allow user to select column names
                gene_column = st.selectbox("Select the column containing gene names:", df_cluster.columns, index=df_cluster.columns.get_loc("Gene") if "Gene" in df_cluster.columns else 0)
                cluster_column = st.selectbox("Select the column containing cluster information:", df_cluster.columns, index=df_cluster.columns.get_loc("Cluster") if "Cluster" in df_cluster.columns else 0)

                # Get unique clusters
                clusters = sorted(df_cluster[cluster_column].unique())
                column_submitted = st.form_submit_button("Set gene and cluster columns")

            with st.form("Choose cluster"):
                selected_cluster = st.multiselect("Cluster:", clusters)
                cluster_submitted = st.form_submit_button("Select clusters")
                st.session_state.cluster_submitted = True


            if st.session_state.cluster_submitted:
#                cluster_genes = df_cluster[df_cluster[cluster_column] == selected_cluster][gene_column].tolist()
                cluster_genes = df_cluster[df_cluster[cluster_column].isin(selected_cluster)][gene_column].tolist()
                gene_list = list(set(cluster_genes))
                st.write(f"{len(gene_list)} genes in cluster {selected_cluster}:")
                st.write('Gene list: ' + ' '.join(gene_list[:10]))
                st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")

                set_back = st.checkbox('Set background genes?', value=False)
                if set_back:
                    n_background = set_back_func()

    elif ORA_mode == 'Input genes':
        st.markdown("##### Or input genes (comma, semicolon, space, CR separated):")
        genes = st.text_input("genes",label_visibility = 'collapsed')
        gene_list = []
        if len(genes) > 0:
#            if ',' not in genes:
#                gene_list = genes.split(' ')
#            else:
#                genes =  ''.join(genes.split()) # Remove spaces
#                gene_list = genes.split(',')
            raw_genes = re.split(r'[,;\s]+', genes)
            # Remove spaces from each gene name and filter out empty strings
            gene_list = [re.sub(r'\s', '', gene) for gene in raw_genes if gene.strip()]

            genes = list(filter(lambda x: x != "", genes)) # Remove spaces
            gene_list =sorted(set(gene_list), key=gene_list.index)
            st.write(gene_list[:3])
            st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")

            set_back = st.checkbox('Set background genes?', value=False)
            if set_back:
                n_background = set_back_func()


    else: # use DEG file

        if "up" not in st.session_state:
            st.session_state.up = None
        if "down" not in st.session_state:
            st.session_state.down = None

        use_upload = 'Yes' # Yes when deseq2 is not available
        df_res = None
        if 'deseq2' in st.session_state:
            st.write("There is a deseq2 result in the cache. If you use it, do not upload a new file.")
            if st.session_state.deseq2 is not None:
                use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
            if use_upload == "No" and "df_res" not in st.session_state: # Create df_res
                if 'df_res' not in st.session_state:
                    df_res = st.session_state.deseq2
                    df_res['Gene'] = df_res.index
                    if  "deseq2_uploaded_file_name" in st.session_state:
                        file_name_head = st.session_state.deseq2_uploaded_file_name
                    else:
                        file_name_head = "res"
                    input_file_type = 'tsv'
                    if "Row_name" in df_res.columns.to_list(): # When Row_name is included
                        df_res = df_res.set_index('Row_name')
                        df_res.index.name = "Gene"
                    st.session_state.df_res = df_res # Record in df_res
            elif "df_res" in st.session_state:
                df_res = st.session_state.df_res
            else:
                st.write("something is wrong...")



        if use_upload == 'Yes': # When using DEG results df_res
    #        st.session_state.deseq2 = None
            input_file_type = st.radio(
                "Data from:",
                ('auto', 'tsv', 'csv', 'excel'))
            uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
            if uploaded_file is not None:
                if input_file_type == 'auto':
                    try:
                        df_res = read_csv(uploaded_file, sep = None)
                    except:# excel
                        df_res = read_excel(uploaded_file)


                if input_file_type == "csv":
                    df_res = read_csv(uploaded_file)
                elif input_file_type == 'tsv':
                    df_res = read_csv(uploaded_file, sep = '\t')
                elif input_file_type == 'excel':
                    df_res = read_xl(uploaded_file)
            #    st.session_state.deseq2 = df_res
                st.session_state.df_res = df_res
                file_name_head = os.path.splitext(uploaded_file.name)[0]
                st.session_state.deseq2_uploaded_file_name = file_name_head
            #    if 'seurat_res' not in st.session_state: # True when Seurat processing is done
            #        st.session_state.seurat_res = False

            else:
                sys.exit(1)
           ##### file upload

        if df_res is not None:
            if 'gene_list' not in st.session_state:
                st.session_state.gene_list = None

            st.write(df_res.head(3))


#            if use_upload == 'Yes': # When Seurat is clicked
    #        if 'seurat_res' not in st.session_state: # True when Seurat processing is done
    #            st.session_state.seurat_res = False
    #        elif not st.session_state.seurat_res: # When Seurat processing is not done yet
            seurat = st.checkbox('Seurat results?', value=False)
            if seurat:
#                df_res.columns = ['col', 'p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene']
#                df_res = df_res.drop("col", axis = 1)

            #    st.write(df_res.head())

                # Identify unique clusters
                clusters = df_res['cluster'].unique()
                # clusters = list(set(df_res['cluster']))
            #    st.write(clusters)
            #    st.write(f"Unique clusters: {clusters}")

                # Create list to hold DataFrames for each cluster
                cluster_dfs = []

                for cluster in clusters:
                    # Filter data for this cluster
                    cluster_data = df_res[df_res['cluster'] == cluster].copy()

            #        st.write(f"Processing cluster {cluster}, shape: {cluster_data.shape}")
                    # Set 'gene' column as index
                    cluster_data = cluster_data.set_index('gene')
                    # Drop 'cluster' column
                    cluster_data = cluster_data.drop('cluster', axis=1)
                    # Change column names to add cluster number at the beginning
                    new_columns = [f'{cluster}_{col}' for col in cluster_data.columns]
                    cluster_data.columns = new_columns
                    # Add cluster DataFrame to list
                    cluster_dfs.append(cluster_data)

                # Merge all cluster DataFrames
                result_df = pd.concat(cluster_dfs, axis=1)

#                # Rearrange columns to make 'gene' the first column
#                cols = ['gene'] + [col for col in result_df.columns if col != 'gene']
#                result_df = result_df[cols]
                df_res = result_df
                df_res["Gene"] = df_res.index.to_list()
                st.write(df_res.head(3))
                # st.session_state.seurat_res = True  This would prevent it from working when returning to the beginning

            content = df_res.columns.tolist()
            p_patterns = ['p.val', 'pvalue', 'p-val', 'p val', 'p_val', 'pval']
            pvalue = [i for i in content if any(p in i.lower() for p in p_patterns)] # and 'adj.pval' not in i.lower()]
            fc_patterns = ['log2fc', 'fold change', 'log2foldchange', 'coef', 'logfc']
            fc = [i for i in content if any(pattern in i.lower() for pattern in fc_patterns)]

            # Auto-detect PCA loadings
            loading_patterns = ['pc', 'loadings', 'component']
            loadings_cols = [i for i in content if any(p in i.lower() for p in loading_patterns)]

            # If no P-value or FDR but has PC columns, treat as PCA loadings
            is_pca_loading = (len(pvalue) == 0 and len(loadings_cols) > 0)

            if is_pca_loading:
                st.info("🔍 PCA loadings file detected (no p-value columns found, PC columns present)")
                ora_mode_auto = 'PCA loadings'

                # If not found, use all numeric columns as candidates
                if not loadings_cols:
                    loadings_cols = df_res.select_dtypes(include=[np.number]).columns.tolist()

                Loading_column = st.selectbox('Select PCA loading column', loadings_cols)
                file_name_add = Loading_column
                P_column = None  # Don't use P_column
                FC_column = None
            else:
                st.info("📊 DEG result file detected (p-value columns found)")
                ora_mode_auto = 'DEG'

                gene = [i for i in content if (i not in pvalue) and (i not in fc)]
                P_column = st.selectbox(
                    'Select adjusted P-value column',
                    pvalue)
                # Jaro-Winkler distance method
                # JW_dist = [Levenshtein.jaro_winkler(P_column, x) for x in fc]

                JW_dist = [jaro_winkler_similarity(P_column, x) for x in fc]
                try:
                    FC_column = st.selectbox(
                        'Select FC column',
                        fc, index = JW_dist.index(max(JW_dist)))
                except:
                    FC_column = st.selectbox(
                        'Select FC column', fc)

                file_name_add = P_column
                file_name_add = file_name_add.replace('.adj.pvalue','')
                file_name_add = file_name_add.replace(' adj. p-value','')

            file_name = file_name_add

            set_gene = st.checkbox("Specify gene column?", value=False)

            if set_gene:
                if is_pca_loading:
                    gene_for_selection = [i for i in content if i not in loadings_cols]
                else:
                    gene_for_selection = gene
                Gene_column =  st.selectbox(
                'Select gene column',
                gene_for_selection)
            else:

                if "Annotation/Divergence" in content:
                    pattern = "([^|]*)"
                    Gene_column = 'Annotation/Divergence'
                    df_res.loc[:,'Annotation/Divergence'] = df_res.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
                    st.write("Converted Annotation/Divergence to gene symbols.")
                elif "Gene" in content:
                    Gene_column =  "Gene"
                else:
                    if is_pca_loading:
                        gene_for_selection = [i for i in content if i not in loadings_cols]
                    else:
                        gene_for_selection = gene
                    Gene_column =  st.selectbox(
                    'Select gene column',
                    gene_for_selection)

            if use_upload == 'Yes':
                # Set index to Gene
                df_res.index = df_res[Gene_column].tolist()

            if 'df_res' in st.session_state:
                st.session_state.df_res = df_res

            with st.form("Basic settings:"):
                if is_pca_loading:
                    # PCA loadings mode: top N selection only
                    st.info("📊 PCA loadings mode: Select top N genes by loading values")
                    up_or_down = st.radio("Positive or negative loadings:", ('Positive','Negative'))
                    top_n = st.number_input('Number of top genes', min_value =1, step = 1, value=50)
                    set_back = st.checkbox('Set background genes?', value=False)
                    st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the analysis are a better background.")
                    ssubmitted_basic = st.form_submit_button("Change the parameters")
                    up=None
                    down=None
                    df_thre  = None

                    # PCA loadings mode processing
                    df_thre = df_res.copy(deep=True)
                    df_thre = df_thre.dropna(subset=[Loading_column])  # Remove NA

                    # Sort by absolute value
                    df_thre = df_thre.sort_values(Loading_column, ascending=False, key=abs)

                    # Separate by Positive or Negative loadings
                    if up_or_down == 'Positive':
                        up = df_thre[df_thre[Loading_column] > 0].index.to_list()
                        down = []
                    else:
                        down = df_thre[df_thre[Loading_column] < 0].index.to_list()
                        up = []

                else:
                    # DEG mode: original behavior
                    p_or_top = st.radio(
                        "P threshold or top/bottom:",
                        ('P_threshold','top', 'both'),
                        help="P_threshold: Select all genes meeting thresholds | top: Select top N genes (ignores P-value) | both: Apply thresholds then select top N | Note: FC threshold (default=0) applies to all modes when set >0"
                    )
                    up_or_down = st.radio("Up or down genes:", ('Up','Down'))
                    sort_val = sort_val = st.radio("Top based on:", ('FC','P value'), index = 0,
                        help="This setting is only used when 'top' or 'both' mode is selected")
                    p_thre = st.number_input('Threshold for adjusted P', min_value =0.000, max_value=1.000,
                    step =0.002, value=0.050)
                    FC_thre = st.number_input('Threshold for log FC', min_value =0.0, step =0.1, value=0.0)

                    top_n = st.number_input('Number of top genes', min_value =1, step = 1, value=50)

                    st.write("Top genes include only those that meet the adjusted P-value threshold.")
                    set_back = st.checkbox('Set background genes?', value=False)
                    st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")
                    ssubmitted_basic = st.form_submit_button("Change the parameters")
                    up=None
                    down=None
                    df_thre  = None

                    if p_or_top == 'top': # Force p_thre to 1 when only top is selected
                        p_thre = 1


                    # Add before creating df_thre, after the "Basic settings:" form
                    # Check and exclude genes with FC=0 and p=0
                    problematic_mask = (df_res[FC_column] == 0) & (df_res[P_column] == 0)
                    if any(problematic_mask):
                        st.warning(f"Found {sum(problematic_mask)} genes with FC=0 and p=0 in the DEG results.")
                        excluded_genes = df_res.loc[problematic_mask, Gene_column].tolist()
                        st.write("These genes will be excluded from ORA analysis:", ', '.join(excluded_genes[:10]),
                                 "..." if len(excluded_genes) > 10 else "")

                        # Exclude from df_res
                        df_res = df_res[~problematic_mask]

                    df_thre = df_res.copy(deep=True)
                    df_thre = df_thre[df_thre[P_column] < p_thre]
                    df_thre = df_thre[abs(df_thre[FC_column]) > FC_thre]
                    if sort_val == "P value":
                        df_thre =df_thre.sort_values(P_column, ascending=True)
                        up = df_thre[df_thre[FC_column]>0].index.to_list()
                        down = df_thre[df_thre[FC_column]<0].index.to_list()
                    else:
                        df_thre = df_thre.sort_values(FC_column, ascending=False)
                        up = df_thre[df_thre[FC_column]>0].index.to_list()
                        df_thre = df_thre.sort_values(FC_column, ascending=True)
                        down = df_thre[df_thre[FC_column]<0].index.to_list()
                # Top N selection (DEG mode only)
                if not is_pca_loading and p_or_top in ["top",  "both"]:
                    up = up[:top_n]
                    down = down[:top_n]

                # Top N selection for PCA loadings mode (already filtered above, just limit here)
                if is_pca_loading:
                    up = up[:top_n]
                    down = down[:top_n]

                up = sorted(set(up), key=up.index) # Remove duplicates
                down = sorted(set(down), key=down.index)

                st.session_state.up = up
                st.session_state.down = down
            try:
                if len(df_thre) > 0:
                    # For PCA loadings mode, handle Positive/Negative
                    if is_pca_loading:
                        if up_or_down == "Positive":
                            st.write(','.join(st.session_state.up))
                            st.write("Number of genes: " + str(len(up)))
                            gene_list = st.session_state.up
                        else:
                            st.write(','.join(st.session_state.down))
                            st.write("Number of genes: " + str(len(down)))
                            gene_list = st.session_state.down
                    else:
                        # DEG mode
                        if up_or_down == "Up":
                            st.write(','.join(st.session_state.up))
                            st.write("Number of genes: " + str(len(up)))
                            gene_list = st.session_state.up
                        else:
                            st.write(','.join(st.session_state.down))
                            st.write("Number of genes: " + str(len(down)))
                            gene_list = st.session_state.down
                    st.session_state.gene_list = gene_list
            except:
                pass



        if set_back:
            # bk_genes = df_res[Gene_column].values
            n_background = len(set(df_res[Gene_column].to_list()))
            st.write('Background gene number: ' + str(n_background))

        gene_list = st.session_state.gene_list

        # Generate file name
        if is_pca_loading:
            gene_list_file = Loading_column + "." + up_or_down + "-" + str(top_n) + '.txt'
        else:
            if p_or_top == 'top':
                gene_list_file = P_column + "." + up_or_down + "-" + str(top_n) + '.txt'
            else:
                gene_list_file = P_column + "." + up_or_down + ".p" + str(p_thre) + '.txt'

        tsv_string = '\n'.join(gene_list)


        st.download_button(
            label="Download the gene list",
            data=tsv_string,
            file_name=gene_list_file,
            mime="text/tsv"
        )


if 'df' in locals()  or 'gene_list' in locals() or 'df_res' in locals(): # When df or genes are entered


    if Analysis_mode == "Rank":
        generated_rnk = False
        use_stat = False
        rank_calc = False
        rank_metric = None
        if not rnk_input:
            # Copy index to Gene column
            df['Gene'] = df.index
            # Remove index name
            df.index.name = None
            content = df.columns.tolist()

            # Auto-detection logic
            p_patterns = ['p.val', 'pvalue', 'p-val', 'p val', 'p_val', 'pval']
            pvalue = [i for i in content if any(p in i.lower() for p in p_patterns) and 'adj.pval' not in i.lower()]
            loading_patterns = ['pc', 'loadings', 'component']
            loadings_cols = [i for i in content if any(p in i.lower() for p in loading_patterns)]

            # If no P-value or FDR but has PC columns, treat as PCA loadings
            is_pca_loading = (len(pvalue) == 0 and len(loadings_cols) > 0)

            if is_pca_loading:
                st.info("🔍 PCA loadings file detected (no p-value columns found, PC columns present)")
                rank_source = 'PCA loadings'
            else:
                st.info("📊 DEG result file detected (p-value columns found)")
                rank_source = 'P values (DEA results)'

            if rank_source == 'P values (DEA results)':
                rank_metric = st.radio(
                    "Ranking metric:",
                    ('sign(LFC) x -log10(P)', 'LFC x -log10(p)', "DESeq2 stat/limma t"), index = 0)
                    # calculate stat value
        if rank_metric ==  "DESeq2 stat/limma t" and not rnk_input and rank_source == 'P values (DEA results)':
            statvalue = [i for i in content if ('stat' in i) or ('t' in i)]
            stat_column = st.selectbox('Select stat column', statvalue)

        elif not rnk_input and rank_source == 'P values (DEA results)':
            generated_rnk = True # When rnk file is created within decoupler
            st.write("Select pvalue and logFC")
            # pvalue is already defined by auto-detection
            fc = [i for i in content if ('log2FC' in i) or ('Fold Change' in i) or ('log2FoldChange' in i) or ('coef' in i) or ('logFC' in i)]
            gene = [i for i in content if (i not in pvalue) and (i not in fc)]
            P_column = st.selectbox('Select P-value column', pvalue)
            stat_column = re.match(r'([^\.]+)', P_column).group(1) # Change name
            # Jaro-Winkler distance method
            JW_dist = [jaro_winkler_similarity(P_column, x) for x in fc]
            try:
                FC_column = st.selectbox(
                    'Select FC column',
                    fc, index = JW_dist.index(max(JW_dist)))
            except:
                FC_column = st.selectbox(
                    'Select FC column', fc)

            if "Gene" in content:
                Gene_column =  "Gene"
            elif "Symbol" in content:
                Gene_column =  "Symbol"
            else:
                Gene_column =  st.selectbox(
                'Select gene column',
                gene)

            inv_switch = st.checkbox('Invert the sign')

            df_sub = df[[P_column, FC_column]]
            df_score = calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch)
            stat_column = 'score'
            mat = df_score.transpose()
            rank_calc = True
            st.write(mat.iloc[:,:10])
            df = df_score

            # Delete files if P_column is changed
            if st.session_state.P_column != P_column:
                shutil.rmtree(dc_temp_dir)
                os.mkdir(dc_temp_dir)
                st.session_state.P_column = P_column

            # Save to session state (for later reference)
            st.session_state.rank_source = 'P values (DEA results)'

        elif not rnk_input and rank_source == 'PCA loadings':
            # PCA loadings mode
            generated_rnk = True
            st.write("Select PCA loadings column")

            # Loading column patterns (PC1, PC2, etc.)
            loading_patterns = ['pc', 'loadings', 'component']
            loadings_cols = [i for i in content if any(p in i.lower() for p in loading_patterns)]

            # If not found, use all numeric columns as candidates
            if not loadings_cols:
                loadings_cols = df.select_dtypes(include=[np.number]).columns.tolist()

            Loading_column = st.selectbox('Select PCA loading column', loadings_cols)
            stat_column = Loading_column

            # Select Gene column
            gene = [i for i in content if i not in loadings_cols]
            if "Gene" in content:
                Gene_column = "Gene"
            elif "Symbol" in content:
                Gene_column = "Symbol"
            else:
                Gene_column = st.selectbox('Select gene column', gene)

            inv_switch = st.checkbox('Invert the sign')

            # Convert loadings values to rank file
            df_score = df[[Gene_column, Loading_column]].copy()
            df_score = df_score.dropna()  # Remove NA
            df_score.columns = ['Gene', 'score']
            df_score = df_score.set_index('Gene')

            # Invert if needed
            if inv_switch:
                df_score['score'] = -1 * df_score['score']

            # Sort by score (descending)
            df_score = df_score.sort_values('score', ascending=False)
            stat_column = 'score'
            mat = df_score.transpose()
            rank_calc = True
            st.write(mat.iloc[:,:10])
            df = df_score

            # Delete files if Loading_column is changed
            if st.session_state.P_column != Loading_column:
                shutil.rmtree(dc_temp_dir)
                os.mkdir(dc_temp_dir)
                st.session_state.P_column = Loading_column  # Save to P_column even in PCA loadings mode

            # Save to session state (for later reference)
            st.session_state.rank_source = 'PCA loadings'
            # Unify to P_column (contains Loading_column name in PCA loadings case)

        else: # rank file
            stat_column = 'Rank'
            df.columns = ['Rank']

        if not rank_calc:

            if list(df.index.duplicated()).count(True) > 0:
                st.markdown("#### There are duplicated genes.")
                st.write('Duplicated genes:' +  ', '.join(list(df[df.index.duplicated()].index)))
                st.write("The first instances will be kept.")
                st.markdown("---")
                df = df[~df.index.duplicated(keep='first')]


            try:
                mat = df[[stat_column]].T
                st.write(mat.iloc[:,:10])
            except:
                st.markdown("#### Error. Prerank file?")
                sys.exit(1)


#--------------------------------------
    st.markdown("---")
    if Analysis_mode == 'Rank':
        species = st.radio("Species:", ('mouse','human'), index = check_species_index(df.index.to_list()[:50]))
        path = st.radio(
        "Pathway:",
        ('PROGENy', 'CollecTRI', 'DoRothEA', 'mSigDB', 'Enrichr', 'Homemade', 'Your own GMT file'))
        st.write("""
        PROGENy: signaling pathway\n
        CollecTRI: TF targets\n
        DoRothEA: TF targets (subset of TRI)
        """)
    else:
        species = st.radio("Species:", ('mouse','human'), index = check_species_index(gene_list[:50]))

        path = st.radio(
            "Pathway:",
            ( 'Pathway_summary','TF_summary', 'Celltype_summary', 'CollecTRI', 'DoRothEA', 'mSigDB', 'Enrichr', 'Homemade', 'Your own GMT file'),
            index = 0)
        st.write("""
        PROGENy: signaling pathway\n
        CollecTRI: TF targets\n
        DoRothEA: TF targets (subset of TRI)
        """)

    if path == 'PROGENy':
        num_progeny = st.radio("Number of top genes:", ('500','2000','5000', 'all'), index = 0,
            help="Number of footprint genes used for calculating each pathway's activity. 500 (high confidence) ~ all (all genes). Default 500 recommended.")
        if st.button('Load PROGENy db') or (num_progeny != st.session_state.num_progeny):
            net = pd.read_csv('./db/progeny.' + species + "." + num_progeny + '.tsv', sep = '\t')
            source='source'
            target='target'
            weight = 'weight'
            st.session_state.net = net
            st.session_state.num_progeny = num_progeny

        elif st.session_state.net is None:
            st.stop()
        else:
            net = st.session_state.net
            source='source'
            target='target'
            weight = 'weight'

    elif path == 'CollecTRI':
        net = pd.read_csv('./db/TRI.' + species + '.tsv', sep = '\t')
        source='source'
        target='target'
        weight = 'weight'
    elif path == 'DoRothEA':
        conf = st.selectbox("Minimal confidence level:", ("A","B","C","D"), index = 1)
        st.write("B is recommended.")
        conf_levels = ['A']
        if conf == 'B':
            conf_levels = ['A','B']
        elif conf == 'C':
            conf_levels = ['A', 'B','C']
        elif conf == 'D':
            conf_levels = ['A', 'B','C', 'D']
        net = pd.read_csv('./db/dorothea.' + species + '.tsv', sep = '\t')
        net = net[net['confidence'].isin( conf_levels)]
        # Rename 'mor' to 'weight' before passing to decoupler
        net = net.rename(columns={'mor': 'weight'})
        source='tf'
        target='target'
        weight = 'weight'


    elif path in ['TF_summary', 'Pathway_summary', 'Celltype_summary']:
        source = 'source'
        target = 'target'
        weight = None

        if path == 'TF_summary':
            st.write("### Transcription Factor Activity Summary")
            databases = get_tf_databases(species)  # TF database definitions
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()

        elif path == 'Pathway_summary':
            st.write("### Pathway Enrichment Summary")
            databases = get_pathway_databases(species)  # Pathway database definitions
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()

        elif path == 'Celltype_summary':
            st.write("### Cell Type Summary")
            databases = get_celltype_databases(species)
            all_results = run_summary_analysis(databases, gene_list, n_background)
            st.stop()



    elif path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
        source='source'
        target='target'
        weight = None
        if path == 'mSigDB':
            if species == 'mouse':
                dir_path = "db/mSigDB_mouse"
            else:
                dir_path = "db/mSigDB"
        elif path == 'Enrichr':
            if species == 'mouse':
                dir_path = "db/enrichr_gmt_mouse"
            else:
                dir_path = "db/enrichr_gmt"
        elif path == 'Homemade':
            if species == 'mouse':
                dir_path = "db/custom_gmt_mouse"
            else:
                dir_path = "db/custom_gmt"

        files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
        files_file.sort()
        key_index = None
        net_list = []
        if path == 'mSigDB':
            key_index = [item for item in files_file if "h.all" in item]
        else:
            key_index = files_file[0]
        GO_name = st.multiselect('Select gene sets',files_file, default = key_index)
        if len(GO_name) > 0:
            for i in GO_name:
                net_list.append(dc.read_gmt(dir_path + "/" + i))
            net = pd.concat(net_list, ignore_index = True)

    else:
        uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt','gmt'])
        if uploaded_gmt is not None:
            source='source'
            target='target'
            weight = None
            GO_name = uploaded_gmt.name
            path = uploaded_gmt.name.replace('.gmt','')
            stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
            s = stringio.read()
            # Spaces cause errors
            s = s.replace(' ', '_')

            with open('temp.gmt', mode='w') as f:
                f.write(s)
            net = dc.read_gmt('temp.gmt')
    #        os.remove("temp.gmt")

        else:
            st.stop()

    if "net" in locals():
        st.write('Pathway data:')
        st.write(net.head(2))
    #-------------------------------------
        if Analysis_mode == 'Rank':
            st.markdown("##### For CollecTRI TF enrichment: ULM, PROGENy pathway enrichment: ULM or MLM")
            st.write("For other weightless datasets (e.g., gmt gene sets), ulm treats each gene set independently while mlm accounts for all of them at the same time during fitting. If your gene sets contain many common targets, it is better to use ulm, as mlm may not be able to run at all.")
            method_index = 2
            if path ==  'CollecTRI' or path == 'DoRothEA':
                method_index = 0
            method = st.radio("Analytical model:", ("ulm", "consensus", "mlm",  "wsum_norm", 'viper', 'GSEA'), index = method_index)
            st.write('https://decoupler-py.readthedocs.io/en/latest/api.html#methods')
            st.write("Consensus score is calculated from ulm, mlm and wsum_norm.")


            min_n = 5
            min_n = int(st.number_input('Minimum number of targets per TF/set', min_value =0, step =1, value=5))
            st.write("Gene sets that have fewer genes in the uploaded data than this number will be ignored.")



            if method != 'GSEA':
                if 'dc_res' not in st.session_state:
                    st.session_state.dc_res = None
                    st.session_state.dc_score = None

                if st.button('Run analysis'): # or not st.session_state.decouplercalc:
                    if method == 'wsum_norm':
                        score, pvalue, norm, corr  = run_method(method= 'wsum_norm', mat=mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)
                        score = norm # Use wsum z-score
                    else:
                        score, pvalue = run_method(method= method, mat=mat, net=net, source=source, target=target, weight=weight, verbose=True, min_n = min_n)

                    res = score.T
                    res['p-value'] = pvalue.iloc[0,:]
                    res['adj.p_value'] = dc.p_adjust_fdr(pvalue.iloc[0,:])
                    res.columns = ['Score','p-value','adj.p_value']

                    res = res.sort_values('Score', ascending=False)
                    st.dataframe(res.head())
                    st.session_state.dc_res = res
                    st.session_state.dc_score = score

                st.markdown("##### You must run the analysis again if the parameters are changed.")
                st.markdown("---")

                if st.session_state.dc_res is not None:
                    res = st.session_state.dc_res
                    score = st.session_state.dc_score

                    with st.sidebar:
                        st.markdown("#### Barplot parameters")
                        if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                            if path == 'Homemade':
                                GO_file_name = '.'.join(GO_name)
                            else:
                                GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                            bar_name_org =  GO_file_name + '_barplot'
                        elif path == 'DoRothEA':
                            bar_name_org =  path + '-' + conf + '_barplot'
                        else:
                            bar_name_org =  path + '_barplot'
                        bar_name_head = st.text_input("Barplot file name: ", value = bar_name_org)
                        bar_name = bar_name_head + ".pdf"

                        # Filtering method selection
                        bar_filter_method = st.radio(
                            "Filter by:",
                            ('Top N', 'Adjusted P-value threshold'),
                            index=0,
                            help="Choose how to filter terms for the barplot"
                        )

                        if bar_filter_method == 'Top N':
                            bar_top = int(st.number_input('How many top terms to show', min_value=1, step=1, value=15))
                            bar_adjp_threshold = None
                        else:
                            bar_adjp_threshold = st.number_input(
                                'Adjusted P-value threshold',
                                min_value=0.0,
                                max_value=1.0,
                                step=0.01,
                                value=0.05,
                                help="Show terms with adj.p_value <= threshold"
                            )
                            bar_top = None

                        bar_vertical = st.checkbox("Vertical plot?", value=True)
                        bar_vc = st.checkbox("Change center value?", value=False)
                        if bar_vc:
                            bar_v_center = st.number_input('Center value', value=0.0)
                        else:
                            bar_v_center = None
                        bar_x_size = st.number_input('X size', min_value =1, value=8)
                        bar_y_size = st.number_input('Y size', min_value =1, value=6)

                    bar_plot_name = dc_temp_dir + "/" + bar_name

                    # Filter score based on adj.p_value if threshold is set
                    if bar_adjp_threshold is not None and 'adj.p_value' in res.columns:
                        significant_terms = res[res['adj.p_value'] <= bar_adjp_threshold].index.tolist()
                        if len(significant_terms) == 0:
                            st.warning(f"No terms with adj.p_value <= {bar_adjp_threshold}. Showing top 15 instead.")
                            bar_top = 15
                            score_filtered = score
                        else:
                            score_filtered = score[[col for col in score.columns if col in significant_terms]]
                            st.info(f"Found {len(significant_terms)} terms with adj.p_value <= {bar_adjp_threshold}")
                    else:
                        score_filtered = score

                    fig_bar = dc.plot_barplot(score_filtered, stat_column, top=bar_top, vertical=bar_vertical, figsize = (bar_x_size, bar_y_size),
                        vcenter = bar_v_center,  return_fig=True) # save = dc_temp_dir + "/" + bar_name,
                    fig_bar.gca().invert_yaxis()
                    fig_bar.savefig(dc_temp_dir + "/" + bar_name,bbox_inches='tight')
                    st.pyplot(fig_bar)
                #        dc.plot_barplot(score, stat_column, top=bar_top, vertical=bar_vertical, figsize = (bar_x_size, bar_y_size),
                #            vcenter = bar_v_center, save = dc_temp_dir + "/" + bar_name)


                    # Create plot with DECOUPLER (with FDR coloring)
                    fig = dc.plot_barplot(score_filtered, stat_column, top=bar_top, vertical=bar_vertical,
                                         figsize=(bar_x_size, bar_y_size), vcenter=bar_v_center, return_fig=True)

                    # Remove existing colorbar
                    # Get all axes
                    axes = fig.axes
                    # Delete last axes (colorbar)
                    if len(axes) > 1:  # If there are axes other than the main plot
                        fig.delaxes(axes[-1])

                    # Get main axes
                    ax = axes[0]

                    # Get pathway names displayed in bar plot (in display order)
                    displayed_pathways = [label.get_text() for label in ax.get_yticklabels()]

                    # Get adj.p_value in display order
                    log_padj = -np.log10(res.loc[displayed_pathways, 'adj.p_value'])

                    # Create colormap
                    cmap = plt.cm.Reds

                    # Change bar colors
                    for i, bar in enumerate(ax.containers[0]):
                        bar.set_color(cmap(log_padj.iloc[i] / log_padj.max()))

                    # Add new colorbar
                    sm = plt.cm.ScalarMappable(cmap=cmap)
                    sm.set_array(log_padj)
                    cbar = fig.colorbar(sm, ax=ax)
                    cbar.set_label('-log10(adj P-value)')

                    ax.invert_yaxis()

                    # Save and display plot
                    plt.tight_layout()
                    fig.savefig(dc_temp_dir + "/FDR_" + bar_name, bbox_inches='tight')
                    st.pyplot(fig)

                    if path == 'PROGENy' or path == 'CollecTRI':
                        tf_list = score.columns.to_list()
                        tf = st.selectbox('Select pathway to visualize', tf_list)
     #                   st.set_option('deprecation.showPyplotGlobalUse', False)

                        with st.sidebar:
                            st.markdown("#### Pathway plot parameters")
                            tf_name_head = st.text_input("TF file name: ", value = tf + '_targets')
                            tf_name = tf_name_head + ".pdf"
                            tf_top = int(st.number_input('How many top TFs to show', min_value =1, step =1, value=15))
                            tf_x_size = st.number_input('TF X size', min_value =1, value=8)
                            tf_y_size = st.number_input('TF Y size', min_value =1, value=6)


                        fig_targets = dc.plot_targets(df, stat=stat_column, source_name=tf, net=net, top=tf_top, figsize = (tf_x_size, tf_y_size),
                        save = dc_temp_dir + "/" + tf_name, return_fig=True)
                        st.pyplot(fig_targets)
                        if path == 'CollecTRI':
                            st.write("Weight shows positive/negative targets of the TF")

                    if generated_rnk:
                        # P_column contains Loading_column name in PCA loadings mode
                        add_head = st.session_state.P_column + '.' # Add name of rank file
                    else:
                        add_head = ""
                    if path == 'mSigDB' or path == 'Enrichr' or path == "Homemade":
                        if path == 'Homemade':
                            GO_file_name = ''.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method + "." + GO_file_name + '.tsv'
                        zip_name = down_file_name+ "." + add_head +  method + "." + GO_file_name
                    else:
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method+ "." + path + '.tsv'
                        zip_name = down_file_name + "." +add_head +  method + "." + path

                    if st.button('Prepare enrichment files to download'):
                        res.to_csv(out_file_name, sep = '\t')
                        shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)

                        with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                            btn = st.download_button(
                                label="Download Results of enrichment",
                            data=fp,
                            file_name=zip_name + ".zip",
                            mime = "zip"
                            )

            else: # GSEA
                if 'dc_gsea' not in st.session_state:
                    st.session_state.dc_gsea = None

                if st.button('Run GSEA analysis'):
                    GSEA_res = run_GSEA_df(mat.T, stat_column, net = net, source=source, target=target, min_n=min_n, seed=42)
                    st.session_state.dc_gsea = GSEA_res

                st.markdown("##### You must run the analysis again if the parameters are changed.")
                st.markdown("---")

                if st.session_state.dc_gsea is not None:
                    GSEA_res = st.session_state.dc_gsea

                    GO_terms = GSEA_res["Term"].to_list()
                    GSEA_res = GSEA_res.sort_values("FDR p-value", ascending = True)
                    nes_show = st.selectbox('Select term to visualize', GO_terms)

                    with st.sidebar:
                        st.markdown("#### GSEA enrichment plot parameters")
                        gsea_title_size = st.number_input('GSEA: title font size', min_value =1, value=14)
                        gsea_legend_size = st.number_input('GSEA: legend font size', min_value =1, value=12)
                        gsea_legend_x = st.number_input('GSEA: legend X pos', min_value =0.00, max_value = 1.00, value=0.2)
                        gsea_legend_y = st.number_input('GSEA: legend y pos', min_value =0.00, max_value = 1.00, value=0.5)
                        st.write('The left bottom is (0, 0). (0-1)')
                        gsea_x_size = st.number_input('GSEA: X size', min_value =1, value=5)
                        gsea_y_size = st.number_input('GSEA: Y size', min_value =1, value=5)

                    fig, d =  dc.plot_running_score(mat.T, stat_column, net = net, source=source, target=target,
                        set_name=nes_show, cmap='RdBu_r', figsize=(gsea_x_size, gsea_y_size), dpi=100, return_fig=True, save=None)
                    # Returns tuple object of fig, ax; d is gene name

                    nes = GSEA_res[GSEA_res['Term']==nes_show]['NES'].iloc[-1]
                    fdr = GSEA_res[GSEA_res['Term']==nes_show]['FDR p-value'].iloc[-1]

                    s = "NES: " + str(nes) + "\nFDR:" + str(fdr)
                    # plt.text(len(mat.T)/10, 18, s, fontsize=14)
                    plt.figtext(gsea_legend_x, gsea_legend_y, s, fontsize = gsea_legend_size) # Written in plot coordinates, from bottom left. 0-1
                    # Extract axis and modify title
                    gsea_title = nes_show.replace("_", " ")
                    fig.axes[0].set_title(gsea_title, wrap=True, fontsize= gsea_title_size)
                    fig.axes[0].set_ylabel("Enrichment Score")

                    st.pyplot(fig)
                    fig.savefig(dc_temp_dir + "/"+ nes_show + '.pdf', format='pdf')

                    st.dataframe(GSEA_res)
                    if generated_rnk:
                        # P_column contains Loading_column name in PCA loadings mode
                        add_head = st.session_state.P_column + '.' # Add name of rank file
                    else:
                        add_head = ""
                    if path == 'mSigDB' or path == 'Enrichr' or path == "Homemade":
                        if path == 'Homemade':
                            GO_file_name = ''.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method + "." + GO_file_name + '.tsv'
                        zip_name = down_file_name+ "." + add_head +  method + "." + GO_file_name
                    else:
                        out_file_name = dc_temp_dir + "/" + down_file_name + "." + method+ "." + path + '.tsv'
                        zip_name = down_file_name + "." +add_head +  method + "." + path

                    if st.button('Prepare GSEA files to download'):
                        GSEA_res.to_csv(out_file_name, sep = '\t')
                        shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)
                        with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                            btn = st.download_button(
                                label="Download GSEA Results",
                            data=fp,
                            file_name=zip_name + ".zip",
                            mime = "zip"
                            )


        else: # ORA

            if "ORA" not in st.session_state:
                st.session_state.ORA = False
            if "ORA_res" not in st.session_state:
                st.session_state.ORA_res = None

            if st.button('Run ORA analysis') or not st.session_state.ORA: # Don't recalculate if button is not pressed.
                try:
                    ORA_res = dc.get_ora_df(gene_list, net, source=source, target=target, n_background=n_background, verbose=False)
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.markdown("#### Error. Possible wrong choice of gene sets or species.")
                    sys.exit(1)
                ORA_res = ORA_res.sort_values('FDR p-value', ascending = True)
                st.session_state.ORA = True
                st.session_state.ORA_res = ORA_res

            st.markdown("##### You must run the analysis again if the parameters are changed.")
            st.markdown("---")

            if st.session_state.ORA_res is not None and st.session_state.ORA:
                ORA_res = st.session_state.ORA_res

                # Set 0s to min p-value
                enr_pvals = ORA_res[['FDR p-value']]
                enr_pvals.index = ORA_res['Term']
                enr_pvals.values[enr_pvals.values == 0] = np.min(enr_pvals.values[enr_pvals.values != 0])
                enr_pvals = -np.log10(enr_pvals)

                with st.sidebar:
                    st.markdown("#### ORA barplot parameters")
                    if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                        if path == 'Homemade':
                            GO_file_name = '-'.join(GO_name)
                        else:
                            GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                        bar_name_org =  GO_file_name + '_barplot'
                    elif path == 'DoRothEA':
                        bar_name_org =  path + '-' + conf + '_barplot'
                    else:
                        bar_name_org =  path + '_barplot'

                    if ORA_mode == 'Cluster info file' and st.session_state.cluster_submitted:
                        bar_name_org = bar_name_org + "." + "_".join(map(str, selected_cluster))
                    bar_name_head = st.text_input("ORA: Barplot file name: ", value = bar_name_org)
                    bar_name = bar_name_head + ".pdf"

                    # ORA Filtering method selection
                    ora_bar_filter_method = st.radio(
                        "ORA: Filter by:",
                        ('Top N', 'FDR threshold'),
                        index=0,
                        help="Choose how to filter terms for the ORA barplot"
                    )

                    if ora_bar_filter_method == 'Top N':
                        bar_top = int(st.number_input('ORA: How many top terms to show', min_value=1, step=1, value=15))
                        ora_fdr_threshold = None
                    else:
                        ora_fdr_threshold = st.number_input(
                            'ORA: FDR threshold',
                            min_value=0.0,
                            max_value=1.0,
                            step=0.01,
                            value=0.05,
                            help="Show terms with FDR p-value <= threshold"
                        )
                        bar_top = None

                    bar_vertical = st.checkbox("ORA: Vertical plot?", value=True)
                    bar_vc = st.checkbox("ORA: Change center value?", value=False)
                    if bar_vc:
                        bar_v_center = st.number_input('ORA: Center value', value=0.0)
                    else:
                        bar_v_center = None
                    bar_x_size = st.number_input('ORA: X size', min_value =1, value=8)
                    bar_y_size = st.number_input('ORA: Y size', min_value =1, value=6)

                # Filter ORA results by FDR if threshold is set
                if ora_fdr_threshold is not None:
                    significant_terms = ORA_res[ORA_res['FDR p-value'] <= ora_fdr_threshold]['Term'].tolist()
                    if len(significant_terms) == 0:
                        st.warning(f"ORA: No terms with FDR p-value <= {ora_fdr_threshold}. Showing top 15 instead.")
                        bar_top = 15
                        enr_pvals_filtered = enr_pvals
                    else:
                        enr_pvals_filtered = enr_pvals.loc[enr_pvals.index.isin(significant_terms)]
                        st.info(f"ORA: Found {len(significant_terms)} terms with FDR p-value <= {ora_fdr_threshold}")
                else:
                    enr_pvals_filtered = enr_pvals

                try:
                    fig_bar2 = dc.plot_barplot(enr_pvals_filtered.T, 'FDR p-value', vertical=bar_vertical, top=bar_top,
                        figsize = (bar_x_size, bar_y_size), vcenter = bar_v_center, return_fig=True)
                    fig_bar2.gca().invert_yaxis()
                    st.pyplot(fig_bar2)
                    st.markdown("#### Activity = -log10(adjP)")
                    st.markdown("###### -log10(0.05) = 1.301")
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.write("Probably little difference in FDR.")
                    st.markdown("#### The following graph likely has no use!")
                    st.markdown("##### -log10(0.05) = 1.301")
                    st.write(enr_pvals_filtered)
                    vmn = enr_pvals_filtered['FDR p-value'].min()
                    vmx = enr_pvals_filtered['FDR p-value'].max()
                    vc = enr_pvals_filtered['FDR p-value'].mean()
                    fig_bar3 = dc.plot_barplot(enr_pvals_filtered.T, 'FDR p-value', vertical=bar_vertical, top=bar_top,
                        figsize = (bar_x_size, bar_y_size), vmin = vmn, vmax=vmx, vcenter =vc, return_fig=True)
                    fig_bar3.gca().invert_yaxis()
                    st.pyplot(fig_bar3)
                    st.markdown("##### Activity = -log10(adjP)")




                # Log-transform
                enr_pvals = -np.log10(enr_pvals)
                # Add count
                ORA_res['count'] = ORA_res['Features'].str.split(';').str.len()
                st.dataframe(ORA_res)

                # Calculate scale based on max count value
                max_count = ORA_res['count'].max()

                # As a rule of thumb, use larger scale for small counts, smaller scale for large counts
                if max_count <= 10:
                    scale = 4 / max_count  # Larger scale for small counts
                elif max_count <= 50:
                    scale = 2 / max_count
                else:
                    scale = 1 / max_count   # Smaller scale for large counts


                # dot plot
                try:
                    fig_dotplot2 = dc.plot_dotplot(
                        ORA_res.sort_values('Combined score', ascending=False).head(bar_top),
                        x='Combined score',
                        y='Term',
                        s='count',
                        c='FDR p-value',
                        scale = scale,
                        figsize = (bar_x_size, bar_y_size),  return_fig=True
                    )
                    st.pyplot(fig_dotplot2)
                    st.write("Combined score = -log10(P)")
                except Exception as e:
                    st.error(f"Error: {str(e)}")
                    st.write("Cannot generate the dot plot")


                if use_upload == "Yes": # Add to file_name
                    file_name_head = os.path.splitext(uploaded_file.name)[0]
                    file_name_add = file_name_head[:12] + "__"
                else:
                    file_name_add = ""
                if path == 'mSigDB' or path == 'Enrichr' or path == 'Homemade':
                    if path == 'Homemade':
                        GO_file_name = '-'.join(GO_name)
                    else:
                        GO_file_name = '-'.join([''.join(i.replace('.gmt','').replace('txt','').split('.')[:3]) for i in GO_name])
                    out_file_name = dc_temp_dir + "/ORA_" + file_name_add + GO_file_name + '.tsv'
                    zip_name = "ORA_" +file_name_add +  GO_file_name
                else:
                    out_file_name = dc_temp_dir + "/ORA_" + file_name_add + path + '.tsv'
                    zip_name = "ORA_" + file_name_add + path

                if ORA_mode == 'Cluster info file' and st.session_state.cluster_submitted:
                    out_file_name =  os.path.splitext(out_file_name)[0] + "." + "_".join(map(str, selected_cluster)) + '.tsv'

                st.write(out_file_name)
                st.write(zip_name)


                if st.button('Prepare ORA files to download'):
                    # Save figures
                    bar_name = "ORA_barplot.pdf"
                    try:
                        fig_bar3.savefig(dc_temp_dir + "/" + bar_name, bbox_inches='tight')
                    except:
                        st.warning("Bar plot not available to save")
                    try:
                        fig_dotplot2.savefig(dc_temp_dir + "/dot_" + bar_name, bbox_inches='tight')
                    except:
                        st.warning("Dot plot not available to save")
                    ORA_res.to_csv(out_file_name, sep = '\t')
                    shutil.make_archive('temp' + "/" + zip_name, format='zip',root_dir= dc_temp_dir)
                    with open('temp' + "/" + zip_name + '.zip', "rb") as fp:
                        btn = st.download_button(
                            label="Download Results of ORA",
                        data=fp,
                        file_name=zip_name + ".zip",
                        mime = "zip",
                        on_click = delete_file([out_file_name, dc_temp_dir + "/" + bar_name])# Delete files when downloaded
                        )

            else:
                st.write("Click 'Run ORA analysis' to proceed.")
