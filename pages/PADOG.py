# PADOG (Pathway Analysis with Down-weighting of Overlapping Genes)
# Requires BiocManager::install("PADOG")
# Uses rpy2==3.5.1

import streamlit as st
import csv
import re
import os
import numpy as np
import pandas as pd
import shutil
import time
import sys
from pathlib import Path
from io import StringIO
from helper_func import clear_old_directories, clear_old_files, make_r_compatible_column_names, remove_after_space, remove_sample_num, check_species_index
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import StrVector
import plotly.express as px


def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []
    min_length = min(len(s) for s in strings)
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    if suffix_length == 0:
        return strings
    return [s[:-suffix_length] for s in strings]


def check_excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep|Oct|Dec|Feb|Nov)')
    index_name = dfx.index.values
    k = 0
    for i in dfx.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.markdown("#### There are Excel-autoconverted gene names")
                st.write("Gene names are not converted.")
                k = 1
            st.write(i)


def read_gmt(file_path):
    """Read a GMT file and return a dict {term_name: [gene1, gene2, ...]}"""
    gene_sets = dict()
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 3 and parts[0].strip():
                gene_sets[parts[0]] = [g for g in parts[2:] if g.strip()]
    return gene_sets


st.set_page_config(page_title="PADOG Pathway Analysis", page_icon="📃")


@st.cache_data
def convert_df(df):
    return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header=0):
    df_xl = pd.read_excel(file, index_col=index_col, header=header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col=index_col, header=0, sep=sep, engine='python')
    return df_c


# --- Temp directory setup ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)


# ===== Main UI =====
st.markdown("## PADOG Pathway Analysis")
st.markdown("""PADOG (Pathway Analysis with Down-weighting of Overlapping Genes) is a method
that down-weights the influence of genes overlapping across multiple pathways.
It evaluates pathway significance using a permutation-based approach.""")

st.markdown("---")

# --- Input mode selection ---
input_mode = st.radio("Input data type:",
    ("Raw count matrix (TMM + logCPM normalization)",
     "Pre-normalized expression matrix (logCPM/voom)"),
    help="For raw counts, TMM normalization followed by logCPM transformation is performed internally. Pre-normalized data such as limma voom output can also be used.")

df = None
file_name_head = "PADOG"

use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Available data")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes', 'No'), index=1)
    if use_upload == "No":
        df = st.session_state.df
        file_name_head = st.session_state.uploaded_file_name
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:, 8:]
        if "Row_name" in df.columns.to_list():
            df = df.set_index('Row_name')
            df.index.name = "Gene"

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio("",
        ('auto', 'Homer', 'tsv', 'csv', 'excel'), index=0, label_visibility='collapsed')
    uploaded_file = st.file_uploader("Choose a file", type=['txt', 'tsv', 'csv', 'xls', 'xlsx'])

    if uploaded_file is not None:
        if file_type == 'auto':
            try:
                df = read_csv(uploaded_file, sep=None)
                content = df.columns.tolist()

                if "Annotation/Divergence" in content:
                    search_word = '([^\ \(]*)\ \(.*'
                    for i in range(1, len(content)):
                        match = re.search(search_word, content[i])
                        if match:
                            content[i] = match.group(1).replace(' ', '_')
                    df.columns = content
                    df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str)
                    pattern = "([^|]*)"
                    repatter = re.compile(pattern)
                    f_annotation = lambda x: repatter.match(x).group(1)
                    df.loc[:, 'Annotation/Divergence'] = df.loc[:, 'Annotation/Divergence'].apply(f_annotation)
                    df = df.loc[:, 'Annotation/Divergence':]
                    st.write("Converted Annotation/Divergence to gene symbols.")

                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content

            except:
                df = read_excel(uploaded_file)
                content = df.columns.tolist()
                if "Annotation/Divergence" in content:
                    search_word = '([^\ \(]*)\ \(.*'
                    for i in range(1, len(content)):
                        match = re.search(search_word, content[i])
                        if match:
                            content[i] = match.group(1).replace(' ', '_')
                    df.columns = content
                    df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str)
                    pattern = "([^|]*)"
                    repatter = re.compile(pattern)
                    f_annotation = lambda x: repatter.match(x).group(1)
                    df.loc[:, 'Annotation/Divergence'] = df.loc[:, 'Annotation/Divergence'].apply(f_annotation)
                    df = df.loc[:, 'Annotation/Divergence':]
                    content = df.columns.tolist()
                    content[0] = 'Gene'
                    df.columns = content
                    st.write("Converted Annotation/Divergence to gene symbols.")
                else:
                    colnames = df.columns.tolist()
                    colnames[0] = 'Gene'
                    df.columns = colnames

        elif file_type != 'excel':
            if file_type == 'csv':
                df = read_csv(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep='\t')
            st.write("Original:")
            st.write(df.head())
            if file_type == 'Homer':
                df = df.iloc[:, 7:]
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(colnames)):
                    match = re.search(search_word, colnames[i])
                    if match:
                        colnames[i] = match.group(1).replace(' ', '_')
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                try:
                    df.iloc[:, 0] = df.iloc[:, 0].apply(f_annotation)
                    df.columns = colnames
                except:
                    st.markdown("### File format error. Non-Homer file?")
            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames
        else:  # excel
            df = read_excel(uploaded_file)
            content = df.columns.tolist()
            if "Annotation/Divergence" in content:
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str)
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:, 'Annotation/Divergence'] = df.loc[:, 'Annotation/Divergence'].apply(f_annotation)
                df = df.loc[:, 'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                st.write("Converted Annotation/Divergence to gene symbols.")
            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames

        df = df.set_index('Gene')
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        df = make_r_compatible_column_names(df)

    else:
        st.stop()


if df is not None:
    st.write('Gene number: ' + str(len(df)))
    st.write(df.head())

    # --- Data preprocessing ---
    is_raw_count = (input_mode == "Raw count matrix (TMM + logCPM normalization)")

    if is_raw_count:
        df = df.astype(float)
        if not float.is_integer(df.iloc[:, 0].sum() * 1000):
            st.markdown("# It is likely that your data are normalized. Please upload unnormalized raw count data, or select 'Pre-normalized' mode.")
        df = df.round(0)

    df = df.loc[~(df == 0).all(axis=1)]
    st.write("All zero count genes are removed.")

    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio("NaN:",
            ('remove NaN containing genes', 'convert to 0'), key='remove_nan')
        if convert_nan == "convert to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')

    check_excel_autoconversion(df)

    if len(df.index.values) != len(set(df.index.values)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene name of the second occurrence has _2 at the end.")
        lis = df.index.values
        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]

    df = make_r_compatible_column_names(df)
    df.columns = df.columns.str.replace('[^A-Za-z0-9]+', '_')
    df.columns = df.columns.str.replace('-', '_')

    st.write(df.head())

    # ===== Group assignment =====
    st.markdown("---")
    st.markdown("### Group assignment")

    condition = [str(i) for i in df.columns.tolist()]
    group_condition = remove_common_suffix(condition)
    group_condition = [remove_sample_num(x) for x in group_condition]

    df_e = pd.DataFrame(group_condition, index=condition, columns=["Group"])

    with st.form("input_groups"):
        st.write('Set groups:')
        edited_df_e = st.data_editor(df_e)
        condition = edited_df_e.iloc[:, 0].tolist()

        # Paired design
        paired_design = st.checkbox('Paired design', value=False,
            help="Check this if samples are paired (e.g., before/after from the same individual)")
        if paired_design:
            block_values = [str(i) for i in range(len(condition))]
            df_block = pd.DataFrame(block_values, index=[str(i) for i in df.columns.tolist()], columns=["Block"])
            st.write("Set block IDs (same ID = paired samples):")
            edited_block = st.data_editor(df_block)
            block_ids = edited_block.iloc[:, 0].tolist()
        else:
            block_ids = None

        submitted = st.form_submit_button("Submit")

    st.write('Group: ' + ' '.join(condition))

    # --- Select 2 groups for comparison ---
    unique_groups = list(dict.fromkeys(condition))

    if len(unique_groups) < 2:
        st.error("At least 2 groups are required for PADOG analysis.")
        st.stop()

    if len(unique_groups) == 2:
        control_group = unique_groups[0]
        disease_group = unique_groups[1]
        st.write(f"Comparison: **{disease_group}** (d) vs **{control_group}** (c)")
    else:
        st.markdown("##### Select 2 groups for comparison (PADOG is a 2-group method):")
        control_group = st.selectbox("Control group (c):", unique_groups, index=0)
        remaining = [g for g in unique_groups if g != control_group]
        disease_group = st.selectbox("Treatment group (d):", remaining, index=0)
        st.write(f"Comparison: **{disease_group}** (d) vs **{control_group}** (c)")

    # Filter to selected 2 groups
    selected_mask = [g in (control_group, disease_group) for g in condition]
    df_subset = df.loc[:, selected_mask]
    condition_subset = [g for g in condition if g in (control_group, disease_group)]

    # Map to "c" and "d"
    padog_group = ["c" if g == control_group else "d" for g in condition_subset]

    # Filter block_ids too if paired
    if paired_design and block_ids is not None:
        block_ids_subset = [b for b, m in zip(block_ids, selected_mask) if m]
    else:
        block_ids_subset = None

    # Validate sample counts
    n_control = padog_group.count("c")
    n_disease = padog_group.count("d")
    st.write(f"Control samples: {n_control}, Treatment samples: {n_disease}")

    if n_control < 2 or n_disease < 2:
        st.error("Each group must have at least 2 samples.")
        st.stop()

    # ===== Sidebar: Gene set selection + PADOG parameters =====
    with st.sidebar:
        st.markdown("## Gene Set Selection")

        # Species detection
        species_idx = check_species_index(df.index.to_list()[:100])
        species = st.radio("Species:", ('mouse', 'human'), index=species_idx)

        db = st.radio("Gene Set Database:", ('mSigDB', 'Homemade', 'Your own GMT file'))

        gene_sets = None

        if db == 'mSigDB':
            if species == 'mouse':
                dir_path = str(Path(__file__).resolve().parent.parent / "db" / "mSigDB_mouse")
            else:
                dir_path = str(Path(__file__).resolve().parent.parent / "db" / "mSigDB")
        elif db == 'Homemade':
            if species == 'mouse':
                dir_path = str(Path(__file__).resolve().parent.parent / "db" / "custum_gmt_mouse")
            else:
                dir_path = str(Path(__file__).resolve().parent.parent / "db" / "custum_gmt")
        else:
            dir_path = None

        if db in ('mSigDB', 'Homemade') and dir_path is not None:
            files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
            files_file.sort()
            key_index = len(files_file) - 1 if db == 'mSigDB' else 0
            GO_name = st.multiselect('Select gene set(s)', files_file,
                                     default=files_file[key_index] if files_file else None)

            if GO_name:
                gene_sets = dict()
                for gmt_name in GO_name:
                    gmt_path = os.path.join(dir_path, gmt_name)
                    gs = read_gmt(gmt_path)
                    gene_sets.update(gs)
        else:
            uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt', 'gmt'])
            if uploaded_gmt is not None:
                stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
                s = stringio.read()
                t = s.split('\n')
                gmt = [x.split('\t') for x in t]
                gene_sets = dict()
                for i in gmt:
                    if len(i) >= 3 and i[0].strip():
                        gene_sets[i[0]] = [g for g in i[2:] if g.strip()]

        if gene_sets is not None:
            st.write(f"Gene sets loaded: {len(gene_sets)}")
            all_gs_genes = set(g for genes in gene_sets.values() for g in genes)
            overlap = len(set(df_subset.index) & all_gs_genes)
            st.write(f"Gene overlap: {overlap} / {len(df_subset)}")
            if overlap <= 10:
                st.error("Gene overlap is too low (<= 10). PADOG requires > 10 overlapping genes. Check species/gene ID format.")

        st.markdown("---")
        st.markdown("## PADOG Parameters")

        NI = st.number_input("Number of permutations (NI)", value=1000, min_value=50, max_value=10000, step=100,
            help="Number of permutations. Higher values yield more accurate results but increase computation time.")

        Nmin = st.number_input("Minimum gene set size (Nmin)", value=3, min_value=1, max_value=50, step=1,
            help="Pathways with fewer genes than this threshold will be excluded.")

        dseed = st.number_input("Random seed", value=1, min_value=0, step=1,
            help="Seed for reproducibility. 0 = no seed.")

        use_parallel = st.checkbox("Parallel processing", value=False,
            help="Parallel computation using multiple cores. Faster but uses more memory.")

    # ===== Run PADOG =====
    st.markdown("---")

    if gene_sets is None or len(gene_sets) == 0:
        st.warning("Please select a gene set (sidebar).")
        st.stop()

    if st.button('Run PADOG', type='primary'):
        with st.spinner('Running PADOG analysis... (permutation-based: may take several minutes)'):
            try:
                # Save expression data to temp file
                df_subset.to_csv(temp_dir + "/expr.tsv", sep='\t')

                # Assign group vector to R
                r_group = StrVector(padog_group)
                ro.r.assign('padog_group', r_group)

                # Assign block vector if paired
                if paired_design and block_ids_subset is not None:
                    r_block = StrVector(block_ids_subset)
                    ro.r.assign('block_vector', r_block)

                # Prepare expression matrix in R
                if is_raw_count:
                    ro.r(f'''
                    library(edgeR)
                    rawdata <- read.csv('{temp_dir}/expr.tsv', sep='\\t', row.names=1)
                    y <- DGEList(counts=rawdata)
                    y <- calcNormFactors(y)
                    logCPM <- cpm(y, log=TRUE, prior.count=2)
                    ''')
                    st.info("TMM normalization + logCPM transformation applied.")
                else:
                    ro.r(f'''
                    logCPM <- as.matrix(read.csv('{temp_dir}/expr.tsv', sep='\\t', row.names=1))
                    ''')
                    st.info("Pre-normalized data used as-is.")

                # Prepare gene sets as R named list
                ro.r('gslist <- list()')
                for term, genes in gene_sets.items():
                    r_genes = StrVector([g for g in genes if g.strip()])
                    ro.r.assign('tmp_genes', r_genes)
                    ro.r.assign('tmp_term', term)
                    ro.r('gslist[[tmp_term]] <- tmp_genes')

                # gs.names=NULL lets PADOG use names(gslist) automatically

                # Set organism
                organism = "hsa" if species == "human" else "mmu"

                # Build PADOG call
                paired_str = "TRUE" if paired_design else "FALSE"
                parallel_str = "TRUE" if use_parallel else "FALSE"
                dseed_str = str(dseed) if dseed > 0 else "NULL"
                block_str = "block = block_vector," if (paired_design and block_ids_subset is not None) else ""

                ro.r(f'''
                library(PADOG)

                result <- padog(
                    esetm = logCPM,
                    group = padog_group,
                    paired = {paired_str},
                    {block_str}
                    gslist = gslist,
                    organism = "{organism}",
                    annotation = NULL,
                    gs.names = NULL,
                    NI = {NI},
                    plots = FALSE,
                    Nmin = {Nmin},
                    verbose = TRUE,
                    parallel = {parallel_str},
                    dseed = {dseed_str}
                )

                # Fill Name from rownames if NA
                if (all(is.na(result$Name))) {{
                    result$Name <- rownames(result)
                }}

                # Add BH-corrected p-values
                result$Ppadog_adj <- p.adjust(result$Ppadog, method="BH")
                result$PmeanAbsT_adj <- p.adjust(result$PmeanAbsT, method="BH")

                # Sort by Ppadog
                result <- result[order(result$Ppadog),]
                ''')

                # Convert R result to pandas
                with localconverter(ro.default_converter + pandas2ri.converter):
                    result_df = ro.conversion.rpy2py(ro.r('result'))

                # PADOG result already has Name, ID, Size, etc. columns
                # rownames are pathway IDs - keep as index or drop if duplicate
                if result_df.index.name is None:
                    result_df.index.name = 'rowname'
                if 'ID' in result_df.columns:
                    result_df = result_df.reset_index(drop=True)
                else:
                    result_df = result_df.reset_index()
                    result_df = result_df.rename(columns={result_df.columns[0]: 'ID'})

                # If Name column is all NaN, fill from ID or index
                if 'Name' in result_df.columns and result_df['Name'].isna().all():
                    if 'ID' in result_df.columns:
                        result_df['Name'] = result_df['ID']
                    else:
                        result_df['Name'] = result_df.index

                # Add gene list for each pathway
                # Get overlap between GMT-derived gene_sets and expression data genes
                expr_genes = set(df_subset.index.astype(str))
                genes_col = []
                name_col_tmp = 'Name' if 'Name' in result_df.columns else 'ID'
                id_col_tmp = 'ID' if 'ID' in result_df.columns else None
                for _, row in result_df.iterrows():
                    pathway_name = str(row[name_col_tmp]) if pd.notna(row[name_col_tmp]) else ""
                    pathway_id = str(row[id_col_tmp]) if id_col_tmp and pd.notna(row[id_col_tmp]) else ""
                    # Search gene set by pathway name or ID from GMT dictionary
                    gs_genes = gene_sets.get(pathway_name, gene_sets.get(pathway_id, []))
                    # Keep only genes present in expression data
                    matched = sorted(set(gs_genes) & expr_genes)
                    genes_col.append(";".join(matched) if matched else "")
                result_df["Genes"] = genes_col

                st.session_state.padog_result = result_df
                st.session_state.padog_gene_sets = gene_sets  # Save for aPEAR integration
                st.success("PADOG analysis completed.")

            except Exception as e:
                st.error(f"PADOG error: {str(e)}")
                st.stop()

    # ===== Results display =====
    if 'padog_result' in st.session_state and st.session_state.padog_result is not None:
        result_df = st.session_state.padog_result

        st.markdown("### PADOG Results")
        st.write(f"Total gene sets analyzed: {len(result_df)}")

        p_threshold = st.number_input("Adjusted P-value threshold (Ppadog_adj):",
            value=0.05, min_value=0.0, max_value=1.0, step=0.01)

        if 'Ppadog_adj' in result_df.columns:
            sig_results = result_df[result_df['Ppadog_adj'] < p_threshold]
        elif 'Ppadog' in result_df.columns:
            sig_results = result_df[result_df['Ppadog'] < p_threshold]
        else:
            sig_results = result_df

        st.write(f"Significant pathways: **{len(sig_results)}**")

        st.dataframe(result_df)

        # Bar plot of significant pathways
        p_col = 'Ppadog_adj' if 'Ppadog_adj' in result_df.columns else 'Ppadog'
        name_col = 'Name' if 'Name' in result_df.columns else 'ID'

        if len(sig_results) > 0:
            top_n = min(30, len(sig_results))
            plot_df = sig_results.head(top_n).copy()
            plot_df['-log10(p)'] = -np.log10(plot_df[p_col].astype(float).clip(lower=1e-300))
            plot_df[name_col] = plot_df[name_col].astype(str)

            fig = px.bar(plot_df, x='-log10(p)', y=name_col,
                         orientation='h',
                         title=f'Top {top_n} Significant Pathways (PADOG, {p_col} < {p_threshold})',
                         text=plot_df['-log10(p)'].round(2))
            fig.update_layout(yaxis={'categoryorder': 'total ascending', 'type': 'category'},
                              height=max(400, top_n * 30))
            fig.update_traces(textposition='outside')
            st.plotly_chart(fig, use_container_width=True)

        # Download
        tsv = result_df.to_csv(index=False, sep='\t')
        st.download_button(
            label="Download results as TSV",
            data=tsv,
            file_name=f"{file_name_head}_PADOG_results.tsv",
            mime="text/tab-separated-values"
        )
