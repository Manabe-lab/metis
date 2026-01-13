import streamlit as st
import csv
import re
import os
import io
import numpy as np
import pandas as pd
import shutil
from helper_func import clear_old_directories
from helper_func import clear_old_files
import time
import sys
import PyWGCNA

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import pymupdf
from PIL import Image
import traceback  # Added at the beginning of the code

import subprocess
from typing import Dict, Any


def run_r_script(df, temp_dir, MEDissThres, deepSplit, soft_power, min_module_size):
    """
    Run the original R WGCNA code using subprocess
    """
    # Save input data
    temp_input = os.path.join(temp_dir, "input_data.csv")
    df.to_csv(temp_input)

    # Create R script content
    r_script = f'''
        library(WGCNA)
        library(flashClust)

        options(stringsAsFactors = FALSE)
        enableWGCNAThreads(nThreads = 10)

        cat("Reading input data...\\n")
        data <- read.csv("{temp_input}", row.names=1)
        print(head(data))
        cat("Initial dimensions (samples x genes):", dim(data), "\\n")

        original_genes = rownames(data)
        original_samples = colnames(data)
        n_genes = ncol(data)
        n_samples = nrow(data)

        input_matrix = t(data)
        cat("Analysis matrix dimensions (genes x samples):", dim(input_matrix), "\\n")

        powers = c(c(1:10), seq(from = 12, to=20, by=2))
        sft = pickSoftThreshold(
            input_matrix,
            powerVector = powers,
            networkType = "signed hybrid",
            verbose = 5
        )

        pdf("{temp_dir}/topology_analysis.pdf", width=12, height=6)
        par(mfrow = c(1,2))
        cex1 = 0.8

        plot(sft$fitIndices[, 1],
             -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             xlab = "Soft Threshold (power)",
             ylab = "Scale Free Topology Model Fit, signed R^2",
             main = paste("Scale independence")
        )
        text(sft$fitIndices[, 1],
             -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
             labels = powers, cex = cex1, col = "red"
        )
        abline(h = 0.80, col = "red")

        plot(sft$fitIndices[, 1],
             sft$fitIndices[, 5],
             xlab = "Soft Threshold (power)",
             ylab = "Mean Connectivity",
             type = "n",
             main = paste("Mean connectivity")
        )
        text(sft$fitIndices[, 1],
             sft$fitIndices[, 5],
             labels = powers,
             cex = cex1, col = "red")
        dev.off()

        power = {soft_power}
        cat("Using user-specified power:", power, "\\n")

        net = blockwiseModules(
            input_matrix,
            power = power,
            maxBlockSize = 5000,
            networkType = "signed hybrid",
            TOMType = "signed",
            minModuleSize = {min_module_size},
            deepSplit = {deepSplit},
            mergeCutHeight = {MEDissThres},
            numericLabels = TRUE,
            pamRespectsDendro = FALSE,
            saveTOMs = FALSE,
            verbose = 10
        )

        moduleLabels = net$colors
        moduleColors = labels2colors(moduleLabels)
        cat("Module detection results:\\n")
        cat("Number of detected modules:", length(unique(moduleColors)), "\\n")
        print(table(moduleColors))

        tryCatch({{
            pdf("{temp_dir}/dendrogram.pdf", width=12, height=6)
            plotDendroAndColors(
                dendro = net$dendrograms[[1]],
                colors = moduleColors,
                main = "Gene dendrogram and module colors",
                dendroLabels = FALSE,
                hang = 0.03,
                addGuide = TRUE,
                guideHang = 0.05
            )
            dev.off()
        }}, error = function(e) {{
            cat("Error in dendrogram plotting:", e$message, "\\n")
        }})

        adjacency = adjacency(input_matrix, power = power, type = "signed hybrid")
        TOM = TOMsimilarity(adjacency, TOMType = "signed")
        write.csv(TOM, "{temp_dir}/TOM.csv")

        results = data.frame(
            gene = original_genes,
            color = moduleColors,
            label = moduleLabels
        )
        write.csv(results, "{temp_dir}/module_colors.csv", row.names = FALSE)

        MEs = moduleEigengenes(input_matrix, colors = moduleColors)$eigengenes
        datME = orderMEs(MEs)
        write.csv(MEs, "{temp_dir}/module_eigengenes.csv")
        write.csv(datME, "{temp_dir}/datME.csv")

        sample_info = data.frame(sample = original_samples)
        write.csv(sample_info, "{temp_dir}/sample_order.csv", row.names = FALSE)

        # Explicitly close all open devices and quit R
        graphics.off()
        quit(save = "no")
    '''

    # Write R script to file
    script_path = os.path.join(temp_dir, "wgcna_analysis.R")
    with open(script_path, 'w') as f:
        f.write(r_script)

    try:
        result = subprocess.run(
            ['conda', 'run', '-n', 'shiny', 'Rscript', script_path],
            check=True,
            capture_output=True,
            text=True
        )
        
        # Get module information from gene-color.tsv
        gene_colors = pd.read_csv(f"{temp_dir}/module_colors.csv")

        # Load TOM and get it as numpy array
        tom = pd.read_csv(f"{temp_dir}/TOM.csv", index_col=0)
        tom_array = tom.values
        
        # Create PyWGCNA object
        pyWGCNA_df = PyWGCNA.WGCNA(
            name=file_name_head, 
            species=species, 
            geneExp=df.T, 
            outputPath=temp_dir + "/",
            save=True,
            MEDissThres=MEDissThres
        )

        # Set TOM
        pyWGCNA_df.TOM = tom_array

        # Set module color information
        pyWGCNA_df.datExpr.var['moduleColors'] = gene_colors['color'].values

        return pyWGCNA_df

    except subprocess.CalledProcessError as e:
        print(f"Error running R script: {e}")
        print("R output:", e.output)
        raise
    finally:
        # After completion, check for remaining R processes and terminate them
        import psutil
        for proc in psutil.process_iter(['pid', 'name']):
            try:
                if proc.name() in ['R', 'Rscript']:
                    proc.kill()
            except (psutil.NoSuchProcess, psutil.AccessDenied):
                continue


def check_and_kill_r_processes():
    """Check and stop existing R processes"""
    import psutil

    r_processes = []
    for proc in psutil.process_iter(['pid', 'name', 'cmdline']):
        try:
            # Search for R-related processes
            if proc.info['name'] in ['R', 'Rscript'] or \
               (proc.info['cmdline'] and any('rpy2' in cmd.lower() for cmd in proc.info['cmdline'])):
                r_processes.append(proc)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue

    if r_processes:
        st.warning(f"Found {len(r_processes)} running R processes. Stopping them before proceeding.")
        for proc in r_processes:
            try:
                proc.terminate()
                proc.wait(timeout=3)  # Wait 3 seconds
            except psutil.TimeoutExpired:
                proc.kill()  # Force termination
            except Exception as e:
                st.error(f"Error killing process {proc.pid}: {e}")
        return True
    return False

def remove_files_in_directory(directory):
    # Check if directory exists
    if not os.path.exists(directory):
        print(f"Directory {directory} does not exist.")
        return

    # Loop through all items in the directory
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        try:
            if os.path.isfile(file_path) or os.path.islink(file_path):
                # Delete if file or symbolic link
                os.unlink(file_path)
            elif os.path.isdir(file_path):
                # Delete recursively if directory
                shutil.rmtree(file_path)
        except Exception as e:
            print(f'Failed to delete {file_path}. Reason: {e}')

def validate_metadata(df_e):
    problems = []
    for col in df_e.columns:
        if df_e[col].isnull().any():
            problems.append(f"Column '{col}' contains missing values.")
        if len(df_e[col].unique()) == 1:
            problems.append(f"Column '{col}' contains only a single value.")
        if len(df_e[col].unique()) == len(df_e):
            problems.append(f"All values in column '{col}' are different. There may be a grouping issue.")
    return problems

def preprocess_metadata(df_e):
    # Check the length of each column
    column_lengths = df_e.apply(len)

    if column_lengths.nunique() != 1:
        st.warning("Metadata column lengths do not match. Adjusting data.")
        # Adjust other columns to match the longest column length
        max_length = column_lengths.max()
        for col in df_e.columns:
            if len(df_e[col]) < max_length:
                # Add '_dummy' values for the missing rows
                df_e[col] = df_e[col].append(pd.Series(['_dummy'] * (max_length - len(df_e[col]))))
        st.write("Metadata are modified.")
        st.write(df_e)

    return df_e

def convert_pdf_to_images(pdf_path):
    pdf = pymupdf.open(pdf_path)
    page = pdf[0]  # Get the first page
    pix = page.get_pixmap()
    img = Image.frombytes("RGB", [pix.width, pix.height], pix.samples)
    return img


def auto_color_assignment(df, column_name):
    unique_values = df[column_name].unique()
    n_unique = len(unique_values)

    # When number of categories is 2 (binary data)
    if n_unique == 2:
        return dict(zip(unique_values, ['blue', 'red']))

    # When number of categories is 3 to 10
    elif 3 <= n_unique <= 10:
        palette = sns.color_palette("husl", n_unique).as_hex()
        return dict(zip(unique_values, palette))

    # When number of categories exceeds 10
    else:
        # Use cyclic colormap for many categories
        cmap = plt.cm.get_cmap("tab20")  # 20-color colormap
        colors = [mcolors.rgb2hex(cmap(i % 20)) for i in range(n_unique)]
        return dict(zip(unique_values, colors))


@st.cache_data
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


st.set_page_config(page_title="WGCNA", page_icon="📃")


@st.cache_data
def read_csv(file, index_col=None, sep=',', header = 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
    return df_c

@st.cache_data
def read_csv2(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_pywgcna(temp_file_path):
   with open(temp_file_path, "wb") as f:
       f.write(uploaded_file.getbuffer())
   pyWGCNA_df = PyWGCNA.readWGCNA(temp_file_path)
   return pyWGCNA_df

def remove_sample_num(i):
    i = re.sub('[_-][^-_]*\d+$', '', i)
    i = re.sub('\d+$', '', i)
    return i

def remove_after_space(i):
    m = re.match(r'([^\ ]+)(\ )+.+',i)
    if m is not None:
        return m.group(1)
    else:
        return i

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
    os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/figures'
    st.session_state.res_dir = res_dir
    os.mkdir(res_dir)

else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/figures'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)


st.markdown("# WGCNA")
st.markdown("## Data should be vst-normalized.")
st.markdown("### R WGCNA would produce more modules than pyWGCNA.")
st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")

with st.sidebar:
    st.markdown("##### Filter out weakly-expressed genes:")
    independentFiltering = st.checkbox('Yes', value= False)
    min_threshold = 0
    max_threshold = 0
    if independentFiltering:
        st.markdown("##### Filter the genes > counts in all samples:")
        min_threshold = st.number_input("count minimum", value = 0.0, label_visibility = 'collapsed')
        min_threshold = float(min_threshold)
        st.markdown("##### Filter the genes > counts in at least one sample:")
        max_threshold = st.number_input("count max", value = 0.0, label_visibility = 'collapsed')
        max_threshold = float(max_threshold)
    MEDissThres = st.number_input("##### MEDissThres:", value = 0.20, min_value=0.05, max_value=1.00, label_visibility="visible",
        help='''Represents the dissimilarity between modules: Module dissimilarity = 1 - correlation coefficient between Module Eigengenes (ME)
Example: MEDissThres = 0.2 corresponds to a correlation coefficient of 0.8
Merging decision:
If the dissimilarity between two modules is less than this threshold, those modules will be merged
Lowering the threshold (e.g., 0.15) means only highly similar modules will be merged
Raising the threshold (e.g., 0.25) means more modules will be merged
''')
    st.write("Module Eigengene Dissimilarity Threshold. Lower values produce more modules.")

    # Add R WGCNA options
    use_R = st.checkbox("Use R WGCNA", value=True,
        help="Use original R implementation of WGCNA via rpy2. This may provide more stable results.")
    if use_R:
        st.markdown("##### R WGCNA Parameters:")
        deepSplit = st.number_input("##### cutree deepSplit:", value = 2, min_value=0, max_value=4,
        help="Higher values tend to produce more, smaller modules. Default: 2, integer.")
        soft_power = st.number_input("Soft Power:", value=6, min_value=1, max_value=30,
            help="Power for soft-thresholding. Default is 6. Higher values make the network more sparse.")
        min_module_size = st.number_input("Minimum Module Size:", value=20, min_value=5, max_value=100,
            help="Minimum number of genes in a module. Default is 20.")

    vis_net = st.checkbox("Visualize module as newtork", value = False)
    PPI = st.checkbox("Idenitfy STRING PPI pairs in each module", value = False)
    top_n = st.number_input("##### Number of top hub genes", value = 20)

    

species = st.radio("Species:", ('mouse','human'))
if species == 'mouse':
    species = 'mus musculus'
    STRING_species = 10090
else:
    species = 'homo sapiens'
    STRING_species = 9606

if 'skip_first' not in st.session_state:
    st.session_state.skip_first = False

if 'df' not in st.session_state:
    st.session_state.df = None

use_upload = 'Yes'
if 'df' in st.session_state or st.session_state.skip_first:
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        if "Row_name" in df.columns.to_list(): # When it contains Row_name
            df = df.set_index('Row_name')
            df.index.name = "Gene"


if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('Homer','tsv','csv','excel','Upload PyWGCNA object and skip the first part'), index = 1, label_visibility = 'collapsed')

    if file_type == 'Upload PyWGCNA object and skip the first part':
        skip_first = True
        df = None
    else:
        skip_first = False
    st.session_state.skip_first = skip_first


    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx', 'p'])
    if uploaded_file is not None and not st.session_state.skip_first:
        if file_type is not 'excel':
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
            df = read_excel(uploaded_file, index_col = 0)
            content = df.columns.tolist()
            if "Annotation/Divergence" in content:
                 # Convert column names
                search_word = '([^\ \(]*)\ \(.*'

                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content # Temporarily rename columns
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # Excel support
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                # Remove everything before annotation/divergence
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

    elif uploaded_file is not None and st.session_state.skip_first:
        st.write("Reading PyWGCNA object...")
        temp_file_path = os.path.join(temp_dir, uploaded_file.name)
    #    with open(temp_file_path, "wb") as f:
    #        f.write(uploaded_file.getbuffer())
    #    pyWGCNA_df = read_pywgcna(temp_file_path)
        pyWGCNA_df = read_pywgcna(temp_file_path)
        st.write("Metadata:")
        st.write(pyWGCNA_df.datExpr.obs)
        pyWGCNA_df.outputPath=temp_dir + "/"
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        condition = [str(i) for i in pyWGCNA_df.datExpr.obs.index.tolist()[:]] # Error prevention
        group_condition = [remove_after_space(x) for x in condition] # Remove text after space
        group_condition = [remove_sample_num(x) for x in group_condition] # Remove trailing numbers

    else:
        sys.exit(1)


if df is not None and not st.session_state.skip_first:

############ Convert hyphen in sample names to underscore (causes R errors)
    if "-" in "".join(df.columns.values):
        st.write("Minus in sample name will be converted to _.")
        new_columns = [x.replace('-','_') for x in df.columns.values]
        df.columns = new_columns
############

    # Handle cases where name starts with a number
    # Change the first character
    numericstart = False
    colnames = df.columns.to_list()
    for i in range(len(colnames)):
        if re.search('^\d+', colnames[i]) is not None:
            numericstart = True
            colnames[i] = "X" + colnames[i]
    if numericstart:
        df.columns = colnames
        st.write("Some sample names start with numbers. They will be converted to X...")

    condition = [str(i) for i in df.columns.tolist()[:]] # Error prevention
    group_condition = [remove_after_space(x) for x in condition] # Remove text after space
    group_condition = [remove_sample_num(x) for x in group_condition] # Remove trailing numbers

    st.write('Original gene number:  ' + str(len(df)))

    # Convert to float with rounding
    df = df.astype(float)

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
############ Convert hyphen in sample names to underscore (causes R errors)
    if "-" in "".join(df.columns.values):
        st.write("Minus in sample name will be converted to _.")
        new_columns = [x.replace('-','_') for x in df.columns.values]
        df.columns = new_columns
############

    if min_threshold > 0:
        df = df[df.apply(min, axis=1) > min_threshold]
    if max_threshold > 0:
        df = df[df.apply(max, axis=1) > max_threshold]

    st.write('Filtered gene number:  ' + str(len(df)))

    st.write(df.head())

if df is not None or (st.session_state.skip_first and pyWGCNA_df):

    if 'condition_col' not in st.session_state:
        st.session_state.condition_col = []

    with st.form("add columns"):
        st.markdown("##### Add conditions, such as genotype, time, age (comma, space, CR separated):")
        genes = st.text_input("genes",label_visibility = 'collapsed')
        gene_list = []
        if len(genes) > 0:
            gene_list = genes.split(' ') # First split by spaces
            gene_list = list(filter(lambda a: a != '', gene_list)) # Remove empty strings
            if ',' in genes:
                gene_list = sum([x.split(',') for x in gene_list],[]) # Flatten with sum(x, [])
            if '\t' in genes:
                gene_list = sum([x.split('\t') for x in gene_list],[])
            if '\n' in genes:
                gene_list = sum([x.split('\n') for x in gene_list],[])
            gene_list = [a for a in gene_list if a != ''] # Remove empty elements
        condition_col = sum([gene_list], [] )
        st.session_state.condition_col = condition_col
        submitted_group = st.form_submit_button("Submit")

    if submitted_group or len(st.session_state.condition_col) > 0:
        with st.form("Edit conditions"):
            df_e = pd.DataFrame(index = condition, columns = condition_col)
            for i in df_e.columns.values:
                df_e[i] = group_condition
            st.write('Set conditions:')
        #    edited_df_e = st.experimental_data_editor(df_e)
            df_e = st.data_editor(df_e)
            submitted = st.form_submit_button("Set")

        condition = df_e.iloc[:,0].tolist()

        for i in df_e.columns.values:
            st.write(' '.join(df_e.loc[:,i].tolist()))

        if not skip_first:
            if (len(condition) != len(df.columns)):
                st.write("The number of group name does not match the data.")

        #    df_condition = pd.DataFrame(condition)
        #    df_batch = pd.DataFrame(batch)

        # Handle erroneous conversions like 1-Mar
            check_excel_autoconversion(df)

            if len(df.index.values) != len(set(df.index.values)):
                st.markdown("#### There are duplicated rows. Converting the names...")
                st.write("The gene name of the second occurrence has _2 at the end.")
                lis = df.index.values
                df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]

        st.markdown("""
    --------------------------------------------------------------------------
            """)



        if st.button('Run WGCNA'):
            # Metadata preprocessing (existing code)
            for col in df_e.select_dtypes(include=['object']).columns:
                df_e[col] = df_e[col].str.replace('_', '.')
            st.write('Using this metadata:')
            st.write(df_e)
            df_e = preprocess_metadata(df_e)
            metadata_problems = validate_metadata(df_e)
            if metadata_problems:
                st.warning("The following issues were found in the metadata:")
                for problem in metadata_problems:
                    st.write(f"- {problem}")
                st.warning("These issues may affect the WGCNA analysis.")

            # Directory preparation (existing code)
            remove_files_in_directory(res_dir)
            if PPI:
                if not os.path.exists(res_dir + "/PPI"):
                    os.mkdir(res_dir + "/PPI")
            if not os.path.exists(res_dir + "/hub"):
                os.mkdir(res_dir + "/hub")

            if not st.session_state.skip_first:
                if use_R:
                    try:
                        # Check and kill any existing R processes
                        if check_and_kill_r_processes():
                            time.sleep(2)
                        
                        # Run the R script
                        run_r_script(
                            df=df,
                            temp_dir=temp_dir,
                            MEDissThres=MEDissThres,
                            deepSplit=deepSplit,
                            soft_power=soft_power,
                            min_module_size=min_module_size
                        )

                        # Copy PDF files generated by R WGCNA to res_dir
                        pdf_files = ['topology_analysis.pdf', 'dendrogram.pdf']
                        for pdf_file in pdf_files:
                            src = os.path.join(temp_dir, pdf_file)
                            dst = os.path.join(res_dir, pdf_file)
                            if os.path.exists(src):
                                shutil.copy2(src, dst)

                        # Display topology analysis plot
                        try:
                            st.write("### Network topology analysis:")
                            img = convert_pdf_to_images(f"{res_dir}/topology_analysis.pdf")
                            st.image(img, use_container_width=True)
                        except Exception as e:
                            st.error(f"Error displaying topology analysis plot: {str(e)}")

                        # Display dendrogram
                        try:
                            st.write("### Module dendrogram:")
                            img = convert_pdf_to_images(f"{res_dir}/dendrogram.pdf")
                            st.image(img, use_container_width=True)
                        except Exception as e:
                            st.error(f"Error displaying dendrogram: {str(e)}")
                        
                        # Load results and create PyWGCNA object
                        module_colors = pd.read_csv(f"{temp_dir}/module_colors.csv")
                        module_eigengenes = pd.read_csv(f"{temp_dir}/module_eigengenes.csv", index_col=0)
                        datME = pd.read_csv(f"{temp_dir}/datME.csv", index_col=0)
                        sample_order = pd.read_csv(f"{temp_dir}/sample_order.csv")
                        tom = pd.read_csv(f"{temp_dir}/TOM.csv", index_col=0)
                        
                        # Create PyWGCNA object
                        pyWGCNA_df = PyWGCNA.WGCNA(
                            name=file_name_head, 
                            species=species, 
                            geneExp=df.T, 
                            outputPath=temp_dir + "/",
                            save=True,
                            MEDissThres=MEDissThres
                        )

                        # Set R WGCNA results
                        pyWGCNA_df.datExpr.var['moduleColors'] = module_colors['color'].values
                        pyWGCNA_df.TOM = tom.values
                        pyWGCNA_df.MEs = module_eigengenes
                        pyWGCNA_df.datME = datME

                    except Exception as e:
                        st.error(f"Error in R WGCNA: {str(e)}")
                        st.error(traceback.format_exc())
                        st.stop()

                else:
                    # Existing PyWGCNA implementation
                    pyWGCNA_df = PyWGCNA.WGCNA(name=file_name_head, 
                                          species=species, 
                                          geneExp=df.T, 
                                          outputPath=temp_dir + "/",
                                          save=True,
                                          MEDissThres=MEDissThres)
                    pyWGCNA_df.preprocess()

                    
                    try:
                        # Convert PDF to image (existing code)
                        st.write(res_dir + "/sample_clustering_cleaning.pdf")
                        img = convert_pdf_to_images(res_dir + "/sample_clustering_cleaning.pdf")
                        st.image(img, use_container_width=True)
                    except Exception as e:
                        st.error(f"An error occurred: {str(e)}")

                    pyWGCNA_df.findModules()
                    st.write("Done finding modules")

            # Common processing from here on (existing code)
            pyWGCNA_df.updateSampleInfo(df_e.astype('object'))


            # Assign colors to each column of metadata
            for column in pyWGCNA_df.datExpr.obs.columns:
                color_dict = auto_color_assignment(pyWGCNA_df.datExpr.obs, column)
                pyWGCNA_df.setMetadataColor(column, color_dict)

            # Confirm color settings
            for col in pyWGCNA_df.datExpr.obs.columns:
                print(f"Color mapping for {col}:")
                for category, color in pyWGCNA_df.metadataColors[col].items():
                    print(f"  {category}: {color}")
            pyWGCNA_df.datExpr.var['gene_name'] = pyWGCNA_df.datExpr.var.index
            pyWGCNA_df.geneExpr.var['gene_name'] = pyWGCNA_df.geneExpr.var.index

            pyWGCNA_df.analyseWGCNA()

            pyWGCNA_df.saveWGCNA()

            pyWGCNA_df.datExpr.var.to_csv(res_dir + '/' + file_name_head +'_gene-color.tsv', sep = '\t')

            try:
                # Convert PDF to image
                img = convert_pdf_to_images(res_dir + "/module-traitRelationships.pdf")
                # Display image
                st.image(img, use_container_width=True)
            except Exception as e:
                st.error(f"An error occurred: {str(e)}")

       #     module_names = pyWGCNA_df.moduleTraitCor.index.tolist()
       #     module_names = [name[2:] for name in module_names]
            module_names = pyWGCNA_df.datExpr.var.moduleColors.unique().tolist()
            for i in module_names:
                st.markdown(f'#### Module: {i}')
                try:
                    # Convert PDF to image
                    img = convert_pdf_to_images(res_dir + "/module_heatmap_eigengene_" + i + ".pdf")
                    # Display image
                    st.image(img, use_container_width=True)
                    img = convert_pdf_to_images(res_dir + "/module_barplot_eigengene_" + i + ".pdf")
                    # Display image
                    st.image(img, use_container_width=True)

                except Exception as e:
                    st.error(f"An error occurred: {str(e)}")

                if vis_net:
                    try:
                        pyWGCNA_df.CoexpressionModulePlot(modules=[i], numGenes=10, numConnections=100, minTOM=0)
                    except Exception as e:
                        st.warning(f"Unable to create network visualization for module {i}. Error: {str(e)}")
                        st.write("Continuing with the rest of the analysis...")

                hub = pyWGCNA_df.top_n_hub_genes(moduleName=i, n=top_n)
                st.write(hub)
                hub.to_csv(res_dir + '/hub/' + i +'_hub.tsv', sep = '\t')

                if PPI:
                    genes = pyWGCNA_df.datExpr.var[pyWGCNA_df.datExpr.var.moduleColors == i]
                    genes = genes.index.astype(str).tolist()
                    PPI_pairs = pyWGCNA_df.request_PPI(genes=genes, species=STRING_species)
                    PPI_pairs.to_csv(res_dir + "/PPI/" +  i +'_PPI.pairs.tsv', sep = '\t')

            if vis_net:
                try:
                    pyWGCNA_df.CoexpressionModulePlot(modules=module_names, numGenes=100, numConnections=1000, minTOM=0, file_name="all")
                except Exception as e:
                    st.warning(f"Unable to create network visualization for all modules. Error: {str(e)}")
                    st.write("Continuing with the rest of the analysis...")


            # Create data for heatmap
            file_list = os.listdir(res_dir)
            top_file = [x for x in file_list if 'rss_Top' in x]
            if len(top_file) > 0:
                top_cont = pd.read_csv(f'{res_dir}/{top_file[0]}', sep = '\t', index_col=0)
                top_rss = list(set(top_cont.to_numpy().flatten()))
                rss = pd.read_csv(f'{res_dir}/rss_cellType.tsv', sep = '\t', index_col = 0)
                rss_df = rss.T.loc[top_rss,:]
                genes = [x.replace("(+)","") for x in rss_df.index.to_list()]
                rss_df.index =genes
                rss_df = rss_df[sample_list ]
                zscore = pd.read_csv(f'{res_dir}/regulon_zscore_results.csv', sep = ',')
                zscore = zscore.set_index('regulon')
                zscore_df = pd.DataFrame(index=genes, columns = sample_list)
                for i in genes:
                    zscore_df.loc[i] = [zscore.loc[i][zscore.loc[i]['cell_type'] == x].values[0][1] for x in sample_list]
                zscore_df.to_csv(f'{res_dir}/{file_name_head}.Top.Z.heatmap.tsv', sep = '\t')
                rss_df.to_csv(f'{res_dir}/{file_name_head}.Top.rss.heatmap.tsv', sep = '\t')

            shutil.move(temp_dir + '/' + file_name_head + '.p', res_dir + "/" + file_name_head + '.p') # Move object to figures

            shutil.make_archive("res", format='zip',root_dir= res_dir)

            if pyWGCNA_df is not None:
                with open("res.zip", "rb") as fp:
                    btn = st.download_button(
                        label="Download Results",
                    data=fp,
                    file_name=file_name_head + "_WGCNA.zip",
                    mime = "zip"
                    )
                try:
                    os.remove(file_name_head + "_WGCNA.zip")
                    shutil.rmtree(temp_dir)
                    os.mkdir(temp_dir)
                except:
                    pass


# All-zero data should be removed before sending data


# Adjust filename when ref is specified?
