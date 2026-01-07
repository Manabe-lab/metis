#!!!!!!!!!!!!!! pip install rpy2==3.5.1  Newer versions cause errors

# Basically calculate using global variables.
# Variables assigned from Python are global variables


import streamlit as st
import csv
import re
import os
import numpy as np
import pandas as pd
import shutil
from PIL import Image
from helper_func import clear_old_directories, clear_old_files, make_r_compatible_column_names, remove_after_space, remove_sample_num
import time
import sys
from rpy2 import robjects
import io
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from itertools import combinations
from rpy2.robjects import StrVector
from rpy2.robjects.lib import grdevices
import matplotlib.pyplot as plt
import pdf2image

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []
    # Get the length of the shortest string
    min_length = min(len(s) for s in strings)
    # Find the length of the common suffix
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break
    # If no common suffix is found, return the original list
    if suffix_length == 0:
        return strings
    # Create a new list with the common suffix removed
    return [s[:-suffix_length] for s in strings]

def format_comparison(comparison_str):
    # Remove parts containing group_factor
    comparison_str = re.sub(r'([-]?\d*\*)?group_factor\w+\s*', '', comparison_str)

    # Parse coefficient string
    coefficients = comparison_str.split()
    pos_group = None
    neg_group = None

    for coef in coefficients:
        if coef.startswith('1*'):
            pos_group = coef[2:]  # Remove '1*'
        elif coef.startswith('-1*'):
            neg_group = coef[3:]  # Remove '-1*'

    if pos_group and neg_group:
        return f"{pos_group} vs. {neg_group}"
    elif pos_group:
        return f"{pos_group} vs. Control"
    elif neg_group:
        return f"Control vs. {neg_group}"
    else:
        return "Comparison information not available"  # Message when parsing fails


def create_integer_contrasts(groups):
    # Get the number of groups
    n = len(groups)

    # Generate all possible combinations of 2 groups
    group_pairs = list(combinations(groups, 2))

    # List to store contrasts
    contrasts = []

    for pair in group_pairs:
        contrast = np.zeros(n, dtype=int)  # Add dtype=int to make integer type
        first_index = groups.index(pair[0])
        second_index = groups.index(pair[1])

        # Set the first group to -1, the second group to +1
        contrast[first_index] = -1
        contrast[second_index] = 1

        contrasts.append(contrast)

    # Create contrast matrix
    contrast_matrix = np.array(contrasts, dtype=int)  # Add dtype=int

    # Convert to DataFrame
    contrast_df = pd.DataFrame(contrast_matrix,
                               columns=groups,
                               index=[f"{pair[1]}_vs_{pair[0]}" for pair in group_pairs])

    return contrast_df

def capture_r_output_as_dataframe(r_code):
    # Settings to convert R data frames to Python pandas DataFrames
    pandas2ri.activate()

    # Execute R code and get result
    result = ro.r(r_code)

    # If result is a data.frame, convert to pandas DataFrame
    if isinstance(result, ro.vectors.DataFrame):
        df = pandas2ri.rpy2py(result)
        return df
    else:
        # If not a data frame, return as string
        return str(result)

def capture_r_output(r_code):
    # Settings for capturing standard output
    old_stdout = sys.stdout
    sys.stdout = io.StringIO()

    try:
        # Execute R code
        ro.r(r_code)
        # Get captured output
        output = sys.stdout.getvalue()
    finally:
        # Restore standard output
        sys.stdout = old_stdout

    return output


# March-1 Sept-1 support
def excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.write("There are Excel-autoconverted gene names")
                k = 1
            autoconvert_flag = True
            st.write("Converting " + i)
            if x.group(2) == "Mar":
                index_name[j] = "March" + x.group(1)
            elif x.group(2) == "Sep":
                index_name[j] = "Sept" + x.group(1)
        j += 1
    dfx.index = index_name
    return(dfx)


def check_excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep|Oct|Dec|Feb|Nov)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.markdown("#### There are Excel-autoconverted gene names")
                st.write("Gene names are not converted.")
                k = 1
            autoconvert_flag = True
            st.write(i)

st.set_page_config(page_title="Calculate limma-voom", page_icon="📃")


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c

@st.cache_data
def calc_barplot(data, ylabel):
    fig, ax = plt.subplots()
    ax = sns.barplot(data=data)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    ax.set_ylabel(ylabel, fontsize = 14)
    return fig


st.sidebar.title("Options")
# Save in temp directory
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    # Delete old directories and files
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


st.markdown("### Use raw count data")

use_sf = False # Use size factor

use_upload = 'Yes'
if 'df' in st.session_state:
    st.write("Available data")
    st.write(st.session_state.df.head())
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        # Homer support
        if "Transcript/RepeatID" in df.columns[0]:
            df = df.iloc[:,8:]
            st.write(df.head())
        if "Row_name" in df.columns.to_list(): # When Row_name is included
            df = df.set_index('Row_name')
            df.index.name = "Gene"

if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('auto', 'Homer','tsv','csv','excel'), index = 0, label_visibility = 'collapsed')
    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])
    use_sf = st.checkbox('Upload Size Factors (Optional)')

    if uploaded_file is not None:

        if file_type == 'auto':
            try:
                df = read_csv(uploaded_file, sep = None)
                st.write("Uploaded file:")
                st.write(df.head())

                content = df.columns.tolist()

                if "Annotation/Divergence" in content:
                     # Convert column names
                    search_word = '([^\ \(]*)\ \(.*'

                    for i in range(1, len(content)):
                        match = re.search(search_word, content[i])
                        if match:
                            content[i] = match.group(1).replace(' ', '_')
                    df.columns = content # Temporarily change names
                    df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel support
                    pattern = "([^|]*)"
                    repatter = re.compile(pattern)
                    f_annotation = lambda x: repatter.match(x).group(1)
                    df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                    # Exclude columns before annotation/divergence
                    df = df.loc[:,'Annotation/Divergence':]
                    st.write("Converted Annotation/Divergence to gene symbols.")
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content

            except:# excel
                df = read_excel(uploaded_file)
                content = df.columns.tolist()
                if "Annotation/Divergence" in content:
                     # Convert column names
                    search_word = '([^\ \(]*)\ \(.*'

                    for i in range(1, len(content)):
                        match = re.search(search_word, content[i])
                        if match:
                            content[i] = match.group(1).replace(' ', '_')
                    df.columns = content # Temporarily change names
                    df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel support
                    pattern = "([^|]*)"
                    repatter = re.compile(pattern)
                    f_annotation = lambda x: repatter.match(x).group(1)
                    df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                    # Exclude columns before annotation/divergence
                    df = df.loc[:,'Annotation/Divergence':]
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
                df = read_csv(uploaded_file, sep = '\t')
            st.write("Original:")
            st.write(df.head())
            if file_type == 'Homer':
                df = df.iloc[:,7:]
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                # Convert column names
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(colnames)):
                    match = re.search(search_word, colnames[i])
                    if match:
                        colnames[i] = match.group(1).replace(' ', '_')
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                try:
                    df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)
                    df.columns = colnames
                except:
                    st.markdown("### File format error. Non-Homer file?")

            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames
        else: # excel
            df = read_excel(uploaded_file)
            content = df.columns.tolist()
            if "Annotation/Divergence" in content:
                 # Convert column names
                search_word = '([^\ \(]*)\ \(.*'

                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content # Temporarily change names
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel support
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                # Exclude columns before annotation/divergence
                df = df.loc[:,'Annotation/Divergence':]
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

        # Remove special characters
        df = make_r_compatible_column_names(df)

    else:
        sys.exit(1)

if df is not None:
    st.write('Original gene number:  ' + str(len(df)))
    st.write(df.head())

    # Convert to float - typo correction
    df = df.astype(float)

    if not float.is_integer(df.iloc[:,0].sum()*1000):
        st.markdown("# It is likely that your data are normalized. Please upload unnormalized raw count data.")

    df = df.round(0)

    df = df.loc[~(df==0).all(axis=1)] # Remove rows where all values are 0

########## Excel support?
    st.write("All zero count genes are removed.")
    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio( "NaN:",
        ('remove Nan containing genes', 'conver to 0' ), key='remove Nan containing genes')
        if convert_nan == "conver to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')
############
    check_excel_autoconversion(df)


    if len(df.index.values) != len(set(df.index.values)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene name of the second occurrence has _2 at the end.")
        lis = df.index.values
        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]

    # Remove special characters
    df = make_r_compatible_column_names(df)

    st.write(df.head())
    total_count = pd.DataFrame(df.sum()[1:])
    total_count.columns= ['Total counts']
    large_var = False
    if max(total_count['Total counts']) > min(total_count['Total counts']) * 2:
        large_var = True
        st.markdown("### Large difference (>2x) in counts")
        import matplotlib.pyplot as plt
        import seaborn as sns
        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']


        f1 = calc_barplot(df_sum.T, ylabel = "Total counts")
        st.pyplot(f1)

        f2 = calc_barplot(np.log1p(df), ylabel = "ln(x+1)")
        st.pyplot(f2)


    with st.sidebar:
        st.markdown("##### Filter out weakly-expressed genes:")
        independentFiltering = st.checkbox('Yes', value= True)
        st.markdown("""Filtering low-expression genes improves FDR calculation.
        If many genes are filtered out, uncheck this for GSEA to rank all genes.""")

        min_threshold = 0
        max_threshold = 0
        sample_threshold = 0

        # edgeR::filterByExpr option
        use_edgeR_filter = st.checkbox("Use edgeR::filterByExpr()", value=False,
                                      help="CPM-based statistical filtering (recommended). Performs appropriate low-expression gene filtering considering group composition. The simple filtering below will be disabled.")
        if use_edgeR_filter:
            st.markdown("##### filterByExpr parameters:")
            filter_min_count = st.number_input("min.count", value=10, min_value=0,
                                              help="Minimum count per sample (default: 10)")
            filter_min_total_count = st.number_input("min.total.count", value=15, min_value=0,
                                                    help="Minimum total count across all samples (default: 15)")
            filter_min_prop = st.number_input("min.prop", value=0.7, min_value=0.0, max_value=1.0, step=0.1,
                                             help="Proportion of samples in smallest group (default: 0.7)")
            filter_large_n = st.number_input("large.n", value=10, min_value=0,
                                            help="Group size considered 'large' (default: 10)")
        else:
            filter_min_count = None
            filter_min_total_count = None
            filter_min_prop = None
            filter_large_n = None

            # Simple filtering options (only shown when filterByExpr is not used)
            if independentFiltering:
                st.markdown("##### Filter the genes > counts in all samples:")
                min_threshold = st.number_input("count minimum", value = 0, label_visibility = 'collapsed')
                min_threshold = int(min_threshold)
                st.markdown("##### Filter the genes > counts in at least one sample:")
                max_threshold = st.number_input("count max", value = 0, label_visibility = 'collapsed')
                max_threshold = int(max_threshold)

            if large_var:
                st.markdown("##### Filter the samples <= counts:")
                sample_threshold = st.number_input("Minimum total cout", value = 0, label_visibility = 'collapsed')

        st.markdown("---")
        st.markdown("##### eBayes options:")
        use_robust_ebayes = st.checkbox("Use robust eBayes", value=True,
                                       help="Use eBayes(robust=TRUE). Performs robust estimation against outliers (recommended)")

        st.markdown("---")
        st.markdown("##### Differential expression test method:")
        test_method = st.radio(
            "Select test method:",
            ["eBayes + topTable (standard)", "TREAT + topTreat (minimum fold-change)"],
            index=0,
            help="""**eBayes + topTable (standard)**: Tests for log2FC=0. Detects even small changes (suitable for GSEA ranking).

**TREAT + topTreat (minimum fold-change)**: Tests for changes exceeding minimum fold-change threshold. Focuses on biologically meaningful genes (suitable for candidate gene selection/validation).

Example: With lfc=0.585, only genes with log2FC>0.585 (1.5-fold or more change) are considered significant. Small changes (log2FC=0.05, etc.) are excluded.

Common thresholds: 0.585 (1.5-fold, standard), 1.0 (2-fold, stringent)

Reference: McCarthy & Smyth (2009) Bioinformatics"""
        )

        use_treat = (test_method == "TREAT + topTreat (minimum fold-change)")

        if use_treat:
            treat_lfc = st.number_input(
                "Log2 fold-change threshold",
                value=0.585,
                min_value=0.0,
                step=0.1,
                help="Minimum log2FC threshold. 0.585=1.5-fold change (standard), 1.0=2-fold change, 1.5=3-fold change. Common thresholds: 0.585 (standard), 1.0 (stringent)"
            )
        else:
            treat_lfc = None

        st.markdown("---")
        st.markdown("##### Voom transformation method:")
        which_voom = st.radio(
            "Select voom method:",
            ['voom', 'voomWithQualityWeights', 'voomLmFit'],
            index=2,
            help="""**voom**: Standard method. Only gene-level precision weighting.

**voomWithQualityWeights**: For datasets with low-quality samples. Sample-level quality weights downweight low-quality samples. Iterative calculation for convergence.

**voomLmFit**: For sparse data (scRNA-seq, etc.). Runs voom and lmFit simultaneously. Efficient. Can use sample.weights option for sample quality weighting."""
        )

        # sample.weights option for voomLmFit
        use_sample_weights = False
        if which_voom == 'voomLmFit':
            use_sample_weights = st.checkbox(
                "Use sample quality weights",
                value=False,
                help="Use sample.weights=TRUE. Performs sample-level quality weighting (similar to voomWithQualityWeights). Recommended if low-quality samples are present."
            )


    # Skip simple filtering when using filterByExpr
    if not use_edgeR_filter:
        if min_threshold > 0:
            df = df[df.apply(min, axis=1) > min_threshold]
        if max_threshold > 0:
            df = df[df.apply(max, axis=1) > max_threshold]

        st.write('Filtered gene number:  ' + str(len(df)))

        if any(df.sum() <= sample_threshold): # Remove columns with count 0
            st.markdown('#### There are the samples that have counts <= ' + str(sample_threshold))
            st.write(", ".join(df.columns[df.sum() <= sample_threshold].to_list()))
            st.write('They are removed. Now data are:')
            df = df.drop(df.columns[df.sum() <= sample_threshold].to_list(), axis = 1)
            st.write(df.head())
    else:
        st.write('Gene number before filterByExpr:  ' + str(len(df)))

    # Remove special characters
    df.columns = df.columns.str.replace('[^A-Za-z0-9]+', '_')
    df.columns = df.columns.str.replace('-', '_')

    condition = [str(i) for i in df.columns.tolist()] # Error prevention
    group_condition = remove_common_suffix(condition) # Remove common suffix
    group_condition = [remove_sample_num(x) for x in group_condition] # Remove trailing numbers and underscores


    df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])
    df_b = pd.DataFrame(condition, index =condition, columns = ["Batch"])

    batch_in = st.checkbox('Setting batch?')
    with st.form("input_groups and batch"):
        st.write('Set groups:')
        edited_df_e = st.data_editor(df_e)

        condition = edited_df_e.iloc[:,0].tolist()

        if batch_in:
            st.write('Set batch:')
            edited_df_b = st.data_editor(df_b)

        if batch_in:
            batch = edited_df_b.iloc[:,0].tolist()
            st.write('Batch: ' + '  '.join(batch))
        else:
            batch = ["No batch"] # When there is no batch

        submitted = st.form_submit_button("Submit")
    st.write('Group: ' + ' '.join(condition))

    if (len(condition) != len(df.columns)):
            st.write("The number of group name does not match the data.")


    for i in df.columns.values:
        a = df.select_dtypes(exclude='number')
        if not a.empty:
            st.write("There is a non-nemeric value in ")
            st.write(a)


    st.markdown("---")



    if st.button('Run limma-voom'):
        with st.spinner('Calculating limma-voom...'):
            df_path = temp_dir + "/df.tsv"

            r_condition = StrVector(condition)
            ro.r.assign('condition', r_condition)

            # Pass batch information to R if available
            if batch_in and batch != ["No batch"]:
                r_batch = StrVector(batch)
                ro.r.assign('batch', r_batch)
                has_batch = True
            else:
                has_batch = False

            df.to_csv(df_path, sep = '\t')
            ro.r(f'''
            library(limma)
            library(edgeR)

            # Read data
            rawdata <- read.csv('{df_path}', sep = '\t', row.names = 1)

            # Create DGEList object
            y <- DGEList(counts=rawdata)

            # Normalize library sizes
            y <- calcNormFactors(y)

            # Create design matrix
            group_factor <- factor(condition)
            ''')

            # Build design matrix formula based on batch
            if has_batch:
                ro.r('''
                batch_factor <- factor(batch)
                design <- model.matrix(~0 + group_factor + batch_factor)
                # Clean up column names
                colnames(design) <- gsub("group_factor", "", colnames(design))
                colnames(design) <- gsub("batch_factor", "batch_", colnames(design))
                ''')
                st.info("ℹ️ Design matrix includes batch adjustment: ~0 + group + batch")
            else:
                ro.r('''
                design <- model.matrix(~0 + group_factor)
                colnames(design) <- levels(group_factor)
                ''')
                st.info("ℹ️ Design matrix: ~0 + group (no batch adjustment)")

            # Apply filterByExpr if enabled
            if use_edgeR_filter:
                st.markdown("#### Applying edgeR::filterByExpr filtering...")
                original_genes = len(df)
                ro.r(f'''
                # Apply filterByExpr with group and parameters
                keep <- filterByExpr(y, design=design,
                                    min.count={filter_min_count},
                                    min.total.count={filter_min_total_count},
                                    min.prop={filter_min_prop},
                                    large.n={filter_large_n})
                y <- y[keep, , keep.lib.sizes=FALSE]
                cat('Genes retained after filterByExpr:', sum(keep), 'out of', length(keep), '\\n')
                ''')
                genes_kept = int(ro.r('sum(keep)')[0])
                st.write(f"**filterByExpr**: {genes_kept} / {original_genes} genes retained")
                st.write(f"- Parameters: min.count={filter_min_count}, min.total.count={filter_min_total_count}, min.prop={filter_min_prop}, large.n={filter_large_n}")


            # Execute different R code based on selected method
            if which_voom != 'voomLmFit':
                ro.r(f'''
                # Apply voom transformation and quality weights
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                v <- {which_voom}(y, design, plot=TRUE)
                dev.off()

                # Fit linear model
                fit <- lmFit(v, design)
                ''')
            else:  # edgeR::voomLmFit
                # Apply sample.weights option
                sample_weights_param = "TRUE" if use_sample_weights else "FALSE"
                ro.r(f'''
                # Perform voom transformation and linear model fitting at once
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                fit <- voomLmFit(y, design, plot=TRUE, sample.weights={sample_weights_param})
                dev.off()
                ''')
                if use_sample_weights:
                    st.info("ℹ️ voomLmFit with sample.weights=TRUE (sample quality weighting enabled)")

            # Display the voom plot
            try:
                with open(ro.r('png_file')[0], "rb") as file:
                    st.image(file.read(), caption='Voom Plot', use_container_width=True)
            except Exception as e:
                st.error(f"Failed to display Voom plot: {str(e)}")

            # Create contrast matrix using makeContrasts with group names
            unique_groups = list(dict.fromkeys(condition))
            group_pairs = list(combinations(unique_groups, 2))

            res = dict()
            for pair in group_pairs:
                comparison_name = f"{pair[1]}_vs_{pair[0]}"

                # Use makeContrasts with group names - handles batch columns automatically
                contrast_formula = f"{pair[1]} - {pair[0]}"

                if use_treat:
                    # Use TREAT for testing against minimum fold-change
                    ro.r(f'''
                    contrast <- makeContrasts({contrast_formula}, levels=design)
                    fit2 <- contrasts.fit(fit, contrast)
                    fit2 <- treat(fit2, lfc={treat_lfc}, robust={str(use_robust_ebayes).upper()})
                    top_table <- topTreat(fit2, sort.by="p", n=Inf)
                    ''')
                    method_label = f"TREAT (lfc={treat_lfc}, robust={use_robust_ebayes})"
                else:
                    # Standard eBayes approach
                    ro.r(f'''
                    contrast <- makeContrasts({contrast_formula}, levels=design)
                    fit2 <- contrasts.fit(fit, contrast)
                    fit2 <- eBayes(fit2, robust={str(use_robust_ebayes).upper()})
                    top_table <- topTable(fit2, sort.by="P", n=Inf)
                    ''')
                    method_label = f"eBayes (robust={use_robust_ebayes})"

                # Convert R data frame to pandas DataFrame
                with localconverter(ro.default_converter + pandas2ri.converter):
                    top_table_df = ro.conversion.rpy2py(ro.r('top_table'))

                top_table_df = top_table_df.reset_index()
                top_table_df = top_table_df.rename(columns={'index': 'gene'})

                st.write(f"**{comparison_name}** ({method_label})")
                st.write(top_table_df.head())

                # Generate and display MA plot
                ro.r(f'''
                png_file <- tempfile(fileext = ".png")
                png(png_file, width = 800, height = 600)
                limma::plotMA(fit2, main=paste("MA Plot -", "{comparison_name}"))
                dev.off()
                ''')

                try:
                    with open(ro.r('png_file')[0], "rb") as file:
                        st.image(file.read(), caption=f'MA Plot - {comparison_name}', use_container_width=True)
                except Exception as e:
                    st.error(f"Failed to display MA plot for {comparison_name}: {str(e)}")

                res[comparison_name] = top_table_df

            # Merge all results
            new_dfs = []
            for key, df in res.items():
                # Rename columns to include the comparison name, but keep 'gene' as is
                df = df.rename(columns={col: f"{key}.{col}" for col in df.columns if col != 'gene'})
                new_dfs.append(df)

            # Use 'gene' column as the key for merging
            merged_df = new_dfs[0]
            for df in new_dfs[1:]:
                merged_df = pd.merge(merged_df, df, on='gene', how='outer')

            # Set 'gene' as the index after merging
            merged_df = merged_df.set_index('gene')

            st.markdown("### Merged results (all comparisons)")
            st.write(f"Total genes: {len(merged_df)}")
            st.write(merged_df.head(20))

            # Convert DataFrame to TSV
            tsv = merged_df.to_csv(index=True, sep='\t')

            # Create download button
            st.download_button(
                label="Download data as TSV",
                data=tsv,
                file_name=file_name_head + ".limma-voom.tsv",
                mime="text/tab-separated-values"
            )
