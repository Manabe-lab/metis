import streamlit as st
import scanpy as sc
import numpy as np
import pandas as pd
import os
import pickle
import re
import time
import matplotlib.pyplot as plt
import seaborn as sns
import zipfile
import io
import concurrent.futures
from joblib import Parallel, delayed
from streamlit_sortables import sort_items

from pages.cellchat import (
    cellchat_analysis, get_cellchatdb_from_r, debug_cellchatdb, check_species_index, optimize_matrix_operations,
    compute_hill_outer_vectorized,precompute_complex_mapping, extractGene,identify_overexpressed_genes,
    calculate_mean_expression_optimized, computeExpr_LR, computeExpr_coreceptor, precompute_gene_expressions, 
    computeExpr_agonist, computeExpr_antagonist, hill_function,
    netAnalysis_computeCentrality, create_cell_color_mapping)

from pages.cellchat_vis import (
    netVisual_circle, netVisual_circle_individual, netVisual_chord_by_pathway,
    plotGeneExpression, plot_lr_contribution,
    netVisual_heatmap, create_chord_diagram_pycirclize, netAnalysis_signalingRole_scatter,
    netAnalysis_signalingChanges_scatter)

# Configure logging first
import logging
import traceback
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger("CellChat")

# Create temp directory for saving files
if "cellchat_temp_dir" not in st.session_state:
    st.session_state.cellchat_temp_dir = True
    cellchat_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    os.mkdir(cellchat_temp_dir)
    st.session_state.cellchat_temp_dir = cellchat_temp_dir
else:
    cellchat_temp_dir = st.session_state.cellchat_temp_dir

@st.cache_data
def adata_read(uploaded_file):
    adata = sc.read_h5ad(uploaded_file)
    return adata

def run_cellchat_for_group(
    adata,
    group_name,
    gene_use,
    complex_input,
    cofactor_input,
    resource,
    celltype_key,
    data_layer,
    min_cells,
    expr_prop,
    n_perms,
    n_cpus,
    features,
    r_patcher,
    population_size,
    fun_type,
    trim

):
    """
    Encapsulate the CellChat analysis for a single group in a function.
    Raises RuntimeError if an error occurs.
    """

    st.write(f"funtype {fun_type}")

    result = cellchat_analysis(
        adata=adata,
        groupby=celltype_key,
#        db=cellchatdb,
        gene_use=gene_use,
        complex_input=complex_input,
        cofactor_input=cofactor_input,
        resource=resource,
        use_layer=data_layer,
        min_cells=min_cells,
        min_cells_expr=min_cells_expr,
        expr_prop=expr_prop,
        pseudocount=1.0,
        trim_threshold=0.05,
        k=0.5,
        n=1,
        fun_type=fun_type,
        raw_use=True,
        nboot=n_perms,
        seed=1,  # Match R default
        n_jobs=n_cpus,
        key_added="cellchat_res",
        trim=trim,
        apply_pval_filter=True,
        features=features,
        r_patcher=r_patcher,
        population_size=population_size,
        do_de=do_de
    )
    # If an error occurred, raise an Exception (instead of using st.error/st.stop here).
    if 'error' in result:
        raise RuntimeError(
            f"Error in {group_name} analysis: {result['error']}\n"
            f"{result.get('traceback','')}"
        )
    return result

# Function to sanitize filenames
def sanitize_filename(filename, max_length=20):
    """Sanitize and limit the length of a filename"""
    filename = re.sub(r'[\\/*?:"<>|]', "_", filename)
    if len(filename) > max_length:
        base, ext = os.path.splitext(filename)
        filename = base[:max_length] + ext
    return filename

# Function to infer species from gene names
def check_species_index(gene_list):
    """Infer species (human or mouse) from gene symbols case pattern"""
    if not gene_list:
        return 0  # Default to human if empty
    
    sample_genes = gene_list[:500] if len(gene_list) > 500 else gene_list
    uppercase_ratios = []
    for gene in sample_genes:
        if not gene or not isinstance(gene, str):
            continue
        uppercase_count = sum(1 for char in gene if char.isupper())
        ratio = uppercase_count / len(gene) if len(gene) > 0 else 0
        uppercase_ratios.append(ratio)
    
    avg_uppercase_ratio = sum(uppercase_ratios) / len(uppercase_ratios) if uppercase_ratios else 0
    return 1 if avg_uppercase_ratio > 0.5 else 0  # 0.5+ = human, less = mouse

# Function to split AnnData by a grouping variable
def split_anndata(adata, split_key, group1, group2):
    """Split AnnData object into two based on a grouping variable"""
    # Create masks for each group
    mask1 = adata.obs[split_key] == group1
    mask2 = adata.obs[split_key] == group2
    
    # Create subsets
    adata1 = adata[mask1].copy()
    adata2 = adata[mask2].copy()
    
    return adata1, adata2

# Function to save analysis results
def save_cellchat_results(result1, result2, group1_name, group2_name, uploaded_filename=None, selected_types=None):
    """Save the CellChat results to a pickle file"""
    try:
        # Create a filename based on the input
        if uploaded_filename:
            base_filename = sanitize_filename(uploaded_filename, 20)
        else:
            base_filename = "cellchat_comparison"
        
        # Create a signal types string if provided
        if selected_types:
            signal_str = "_".join([sig_type[:3] for sig_type in selected_types])
            filename = f"{base_filename}_{group1_name}_vs_{group2_name}_{signal_str}.pkl"
        else:
            filename = f"{base_filename}_{group1_name}_vs_{group2_name}.pkl"
        
        # Create a full results object
        combined_results = {
            'group1': {
                'name': group1_name,
                'result': result1
            },
            'group2': {
                'name': group2_name,
                'result': result2
            },
            'upload_filename': uploaded_filename,
            'selected_types': selected_types,
            'created_at': time.time()
        }
        
        # Serialize to bytes
        output = pickle.dumps(combined_results)
        
        return output, filename
    except Exception as e:
        st.error(f"Error saving results: {str(e)}")
        import traceback
        st.error(traceback.format_exc())
        return None, None

# Function to filter result by selected pathways
def filter_result_by_pathways(result, selected_pathways):
    """Filter CellChat result to include only selected pathways"""
    filtered_result = result.copy()
    
    # Check if netP exists
    if 'netP' not in result or 'pathways' not in result['netP']:
        return filtered_result
    
    # Get indices of selected pathways that exist in this result
    all_pathways = list(result['netP']['pathways'])
    pathway_indices = [i for i, p in enumerate(all_pathways) if p in selected_pathways]
    
    # If no selected pathways exist in this result, return empty matrices
    if not pathway_indices:
        # Create empty matrices with same structure
        if 'network' in filtered_result and 'strength_matrix' in filtered_result['network']:
            cell_names = filtered_result['network']['strength_matrix'].index
            filtered_result['network']['strength_matrix'] = pd.DataFrame(
                0, index=cell_names, columns=cell_names)
            filtered_result['network']['count_matrix'] = pd.DataFrame(
                0, index=cell_names, columns=cell_names)
        return filtered_result
    
    # Filter netP probability matrix (dim: sources x targets x pathways)
    if 'prob' in result['netP']:
        filtered_prob = result['netP']['prob'][:, :, pathway_indices]
        
        # Sum across selected pathways to get aggregated network
        aggregated_prob = np.sum(filtered_prob, axis=2)
        
        # Update network matrices
        if 'network' in filtered_result:
            cell_names = result['net']['dimnames'][0]
            
            # Create strength matrix from aggregated probabilities
            filtered_result['network']['strength_matrix'] = pd.DataFrame(
                aggregated_prob,
                index=cell_names,
                columns=cell_names
            )
            
            # Create count matrix (binary: 1 if interaction exists, 0 otherwise)
            filtered_result['network']['count_matrix'] = pd.DataFrame(
                (aggregated_prob > 0).astype(int),
                index=cell_names,
                columns=cell_names
            )
    
    return filtered_result

# Function to save multiple files as ZIP
def create_zip_from_files(file_paths, zip_filename="cellchat_files.zip"):
    """Create a ZIP file from multiple file paths"""
    try:
        # Create an in-memory ZIP file
        zip_buffer = io.BytesIO()
        
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for file_path in file_paths:
                if os.path.exists(file_path):
                    file_name = os.path.basename(file_path)
                    zip_file.write(file_path, file_name)
                else:
                    st.warning(f"File not found: {file_path}")
        
        # Reset buffer position
        zip_buffer.seek(0)
        
        return zip_buffer
    except Exception as e:
        st.error(f"Error creating ZIP file: {str(e)}")
        import traceback
        st.error(traceback.format_exc())
        return None

# Function to find first matching item in a list or return default
def find_first_index_or_default(lst, elements, default=0):
    """Find the first matching element in a list or return default index"""
    for element in elements:
        try:
            return lst.index(element)
        except ValueError:
            continue
    return default


def identify_union_overexpressed_genes(adata1, adata2, gene_use,
                            group_by=None, 
                            min_cells=10,
                            expr_prop=0,
                            thresh_p=0.05,
                            use_union_genes="union"):


    st.write(f"original data1: {adata1.shape[0]} cells x {adata1.shape[1]} genes")

    # 1. Check cell counts in cell groups
   # cell_counts = adata1.obs[group_by].value_counts()
  #  valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
 #   if len(valid_groups) < len(cell_counts):
 #       st.warning(f"{len(cell_counts) - len(valid_groups)} cell groups with fewer than {min_cells} cells have been excluded")                   
    # Keep only cells in valid groups
  #  adata1_filtered = adata1[adata1.obs[group_by].isin(valid_groups)].copy()
    valid_signaling_genes = [g for g in gene_use if g in adata1.var_names]
    # Keep only signaling genes (equivalent to R's subsetData)
    adata1_filtered = adata1[:, valid_signaling_genes].copy()
    st.write(f"LR gene data of data1: {adata1_filtered.shape[0]} cells x {adata1_filtered.shape[1]} genes")

   # cell_counts = adata2.obs[group_by].value_counts()
   # valid_groups = cell_counts[cell_counts >= min_cells].index.tolist()
   # if len(valid_groups) < len(cell_counts):
   #     st.warning(f"{len(cell_counts) - len(valid_groups)} cell groups with fewer than {min_cells} cells have been excluded")                   
    # Keep only cells in valid groups
   # adata2_filtered = adata2[adata2.obs[group_by].isin(valid_groups)].copy()
    valid_signaling_genes = [g for g in gene_use if g in adata2.var_names]
    # Keep only signaling genes (equivalent to R's subsetData)
    adata2_filtered = adata2[:, valid_signaling_genes].copy()
    st.write(f"LR gene data of data2: {adata2_filtered.shape[0]} cells x {adata2_filtered.shape[1]} genes")


    # Identify overexpressed genes in group 1
    result1 = identify_overexpressed_genes(
        adata1_filtered,
        group_by=group_by,
        thresh_pct=expr_prop,
        thresh_p=thresh_p,
        only_pos=True,
        do_de=do_de,
        do_fast=True,
        min_cells=min_cells
    )
    
    # Identify overexpressed genes in group 2
    result2 = identify_overexpressed_genes(
        adata2_filtered,
        group_by=group_by,
        thresh_pct=expr_prop,
        thresh_p=thresh_p,
        only_pos=True,
        do_de=do_de,
        do_fast=True,
        min_cells=min_cells
    )
    
    # Create union of overexpressed genes from both groups
    genes1 = set(result1['features'])
    genes2 = set(result2['features'])
    if use_union_genes=="union":
        union_genes = list(genes1.union(genes2))
    elif use_union_genes=="intersection":
        union_genes = list(genes1 & genes2)
    else:
        union_genes = None
    
    return union_genes, genes1, genes2



if __name__ == "__main__":
    # Set page config
    st.set_page_config(page_title="CellChat Comparison Analysis", page_icon="💬", layout="wide")

    # Main title
    st.title("CellChat Comparison Analysis")

    # Initialize and manage analysis_mode
    if 'completed_analysis' not in st.session_state:
        st.session_state.completed_analysis = False

    # Dynamically set default value for mode selection
    default_mode_index = 1 if st.session_state.get('completed_analysis', False) else 0

    analysis_mode = st.radio(
        "Choose analysis mode:",
        ["Upload H5AD file for new analysis", "Load saved results for visualization"],
        index=default_mode_index
    )

    # Mode 1: Upload new H5AD file for analysis
    if analysis_mode == "Upload H5AD file for new analysis":
        st.header("Data Settings")

        uploaded_file = st.file_uploader("Upload H5AD file", type=['h5ad'], 
                                       help="AnnData file for CellChat analysis")

        if uploaded_file is not None:
            # Try to read the file
            try:
#                with st.spinner('Reading AnnData file...'):
#                    adata = sc.read_h5ad(uploaded_file)
                
                adata = adata_read(uploaded_file)
                
                st.success(f"Successfully loaded file: {uploaded_file.name}")
                
                # Display summary info
                st.write(f"Cells: {adata.n_obs}, Genes: {adata.n_vars}")
                
                # Show available categorical variables for cell type
                categorical_cols = [col for col in adata.obs.columns if pd.api.types.is_categorical_dtype(adata.obs[col])]
                other_cols = [col for col in adata.obs.columns if col not in categorical_cols]
                
         
                # Select cell type annotation column
                default_celltype_idx = find_first_index_or_default(
                    categorical_cols, 
                    ['cell.ident', 'cell_type', 'celltype', 'seurat_clusters'],
                    0
                )
                celltype_key = st.selectbox(
                    "Select cell type annotation column:",
                    categorical_cols,
                    index=min(default_celltype_idx, len(categorical_cols)-1) if categorical_cols else 0,
                    help="Column to use for cell type annotations"
                )
                cell_list = sorted(adata.obs[celltype_key].cat.categories.to_list() if hasattr(adata.obs[celltype_key], 'cat') else sorted(adata.obs[celltype_key].unique()))
                st.write(", ".join(cell_list))


                # Select group comparison column
                default_group_idx = find_first_index_or_default(
                    categorical_cols, 
                    ['orig.ident', 'stimulus', 'condition', 'treatment', 'sample'],
                    0
                )
                splitby_key = st.selectbox(
                    "Select grouping variable for comparison:",
                    categorical_cols,
                    index=min(default_group_idx, len(categorical_cols)-1) if categorical_cols else 0,
                    help="Column to split the dataset for comparison"
                )

                # Get available groups from the selected grouping variable
                if splitby_key:
                    available_groups = adata.obs[splitby_key].cat.categories.tolist() if hasattr(adata.obs[splitby_key], 'cat') else sorted(adata.obs[splitby_key].unique())
                    
                    # If there are exactly 2 groups, select them by default
                    if len(available_groups) == 2:
                        group1 = available_groups[0]
                        group2 = available_groups[1]
                    else:
                        # Otherwise let the user select
                        group1 = st.selectbox("Select first group", available_groups, index=0)
                        remaining_groups = [g for g in available_groups if g != group1]
                        group2 = st.selectbox("Select second group", remaining_groups, index=0)
                    
                    # Display cell counts for selected groups
                    g1_count = (adata.obs[splitby_key] == group1).sum()
                    g2_count = (adata.obs[splitby_key] == group2).sum()
                    st.write(f"Group 1 ({group1}): {g1_count} cells")
                    st.write(f"Group 2 ({group2}): {g2_count} cells")
                    
                    # Store the selected groups in session state
                    st.session_state.group1 = group1
                    st.session_state.group2 = group2
                    st.session_state.celltype_key = celltype_key
                    st.session_state.splitby_key = splitby_key

                
                # Available data layers
                available_layers = list(adata.layers.keys()) if hasattr(adata, 'layers') else []
                st.write(f"Available data layers: {', '.join(available_layers) if available_layers else 'No additional layers'}")
                
                data_layer = st.selectbox(
                    "Select data layer:",
                    ['X (default)'] + available_layers,
                    index=0,
                    help="Data layer to use for analysis"
                )
                
                # Store the selected layer in session state
                st.session_state.data_layer = None if data_layer == 'X (default)' else data_layer
                
                # Infer species from gene names
                species = st.radio("Species:", ('mouse', 'human'), 
                                  index=check_species_index(adata.var.index.to_list()[:50]))

                # Store species in session state
                st.session_state.species = species
                
                # Store AnnData in session state
                st.session_state.adata = adata
                
            except Exception as e:
                st.error(f"Error loading file: {str(e)}")

        # Main panel for analysis options and results
        if 'adata' in st.session_state:
            adata = st.session_state.adata
            
            # Analysis settings form
            with st.form("analysis_settings"):
                st.header("Analysis Settings") 
                col1, col2 = st.columns(2)            
                with col1:
                    species = st.radio("Species:", ('mouse', 'human'), 
                                      index=check_species_index(adata.var.index.to_list()[:50]))
                    
                    # Define available_layers within the form scope
                    available_layers = list(adata.layers.keys()) if hasattr(adata, 'layers') and adata is not None else []
                    data_layer = st.selectbox("Using layer (should be normalized data):", 
                                             ['X (default)'] + available_layers,
                                             index=0)
                    
                    selected_types = st.multiselect(
                        "Signaling types to analyze (can choose multiple):",
                        ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling"],
                        default=["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"]
                    )

                    st.markdown("---")
                    st.markdown("##### Do not change the following parameters unless necessary.")
                    n_cpus = st.slider("Num CPUs:", 
                                  min_value=1, max_value=os.cpu_count(), step=1, value=os.cpu_count()-2)
                   # r_patcher = st.checkbox("Fine-tune to match R calc results", value=False,
                   #     help="Absorbs some differences from R calculation results. Weight approximates even without checking. Count increases if not checked, but trend matches. Count threshold criteria are somewhat arbitrary.")
                    r_patcher = False

                    n_perms = st.slider("Permutation number:", 
                                min_value=20, max_value=500, step=20, value=100)     


                    min_cells = st.number_input('Min cells in a cell-type:', 
                                               min_value=0, max_value=100, step=1, value=10,
                                               help="Filtering is off in SCALA's cellchat.")

                    expr_prop = st.number_input('Min fraction of expressing cells:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.0)
                    st.markdown('---')
                    population_size = st.checkbox("population.size", value=False, help="Use when you want to evaluate communication between larger clusters more strongly (with higher probability).")



                with col2:

                    do_de = st.checkbox("Filter pathways by differential expression in at least one cell type.",
                                        value=True, help="Filter for DEGs in at least one cluster. If off, filter by number of expressing cells.")

                    thresh_p = st.number_input('Threshold p value for overexpressed genes:', 
                                               min_value=0.0, max_value=0.9, step=0.01, value=0.05,
                                               help = "Threshold for selecting changed genes by Wilcoxon test")

                    # Indent with nested columns
                    subcol1, subcol2 = st.columns([0.1, 0.9])
                    with subcol2:
                        st.write("If DEG filtering is off:")
                        min_cells_expr = st.number_input('Gene filtering for min expressing cells:', 
                                                   min_value=0, max_value=100, step=1, value=10,
                                                   help="Filter by number of expressing cells among all cells when DEG filtering is off. Must be quite low or filtering becomes stricter.")
                   
                    st.markdown("---")
                    fun_type = st.radio("Method to calculate average expression per cell group", ['triMean','truncatedMean'], index = 0, help="The number of predicted ligand-receptor pairs depends on how average gene expression per cell group is calculated. trimean generates fewer interactions than other methods. CellChat is known to perform well in predicting stronger interactions, which is very useful for narrowing down interactions for experimental validation. trimean approximates a 25% trimmed mean, meaning that if the proportion of expressing cells in a group is less than 25%, the average gene expression becomes zero. With truncatedMean, you can set 5% or 10% trimmed means, etc.")

                    st.markdown("*To reduce stringency for pathway selection, use truncatedMean.*")
                    trim = st.number_input("Trimming for truncated mean:", min_value=0.00, max_value=0.25, value=0.10)

                    st.markdown("---")                    
                    apply_pval_filter = st.checkbox("Filter interaction pathways by P value", value=True)

                    trim_threshold = st.number_input('Filtering P threshold:', 
                                               min_value=0.01, max_value=0.20, step=0.01, value=0.05)
                    st.markdown('---')
      
                    use_union_genes = st.radio(
                        "Use individual overexpressed genes, or union or intersection from both groups",
                        ['individual','union','intersection'],
                        index=0,
                        help="Identifies overexpressed genes individually or from both groups and uses their union/intersect for analysis. In CellChat tutorial, data is filtered individually by DEG etc., and LR analysis is performed on those genes. Therefore, some genes may not be included in the analysis. Also, LR pairs that are not significant in individual analyses will have an interaction score of 0."
                    )
                
                submit_button = st.form_submit_button("Run Comparison Analysis")              
            
            if submit_button:
                try:
                    # Fetch variables from session state
                    group1 = st.session_state.group1
                    group2 = st.session_state.group2
                    celltype_key = st.session_state.celltype_key
                    splitby_key = st.session_state.splitby_key
                    data_layer = st.session_state.data_layer
                    species = st.session_state.species
                    
                    # Split AnnData into two groups
                    with st.spinner(f"Splitting data into {group1} and {group2} groups..."):
                        adata1, adata2 = split_anndata(adata, splitby_key, group1, group2)
                        
                        # Check if either group has too few cells
                        if adata1.n_obs < min_cells:
                            st.error(f"Group '{group1}' has fewer than {min_cells} cells.")
                            st.stop()
                        if adata2.n_obs < min_cells:
                            st.error(f"Group '{group2}' has fewer than {min_cells} cells.")
                            st.stop()
                        
                        st.success(f"Successfully split data: {group1} ({adata1.n_obs} cells) and {group2} ({adata2.n_obs} cells)")
                    
                    # Get CellChatDB
                    with st.spinner(f"Getting {species} CellChatDB..."):
                        cellchatdb = get_cellchatdb_from_r(species=species)
                   #     cellchatdb = debug_cellchatdb(cellchatdb)
                        st.success(f"{species} CellChatDB successfully obtained")
                        
                        # Filter CellChatDB to selected signaling types
                        interaction = cellchatdb['interaction']
                        interaction_filtered = interaction[interaction['annotation'].isin(selected_types)]
                        cellchatdb['interaction'] = interaction_filtered
                        
                        st.write(f"Using {len(cellchatdb['interaction'])} ligand-receptor pairs")

                    # Expand database here
                    gene_use, resource, complex_input, cofactor_input, gene_info = extractGene(cellchatdb)

                    # Validate database
                    if resource.empty:
                        raise ValueError("The interaction information in the database is empty. Please verify that you have a valid CellChatDB.")
                        
                    # Check for required columns
                    for required_col in ['ligand', 'receptor']:
                        if required_col not in resource.columns:
                            raise ValueError(f"Resource dataframe requires '{required_col}' column")

                    st.write(f"LR genes in db: {len(gene_use)}")
                    
                    # Set up progress tracking and results container
                    st.write("Running cellchat analysis...")
                    progress = st.progress(0)
                    results_container = st.container()
                    
                    # Create status expanders for each group
                    group1_status = st.expander(f"Status for {group1} analysis", expanded=True)
                    group2_status = st.expander(f"Status for {group2} analysis", expanded=True)

                    if use_union_genes != "individual":
                        # First do preprocessing
                        with st.spinner("Identifying overexpressed genes from both groups..."):
                            union_genes, genes1, genes2 = identify_union_overexpressed_genes(
                                adata1, adata2, gene_use,
                                group_by=celltype_key, 
                                min_cells=min_cells,
                                expr_prop=expr_prop,
                                thresh_p=thresh_p,
                                use_union_genes=use_union_genes
                            )

                    if use_union_genes=="union":
                        st.write(f"Union of overexpressed genes: {len(union_genes)}")
                    elif use_union_genes=="intersection":
                        st.write(f"Intersection of overexpressed genes: {len(union_genes)}")

                    # Somewhere in your Streamlit code, after you split adata into adata1, adata2, etc.
                    progress.progress(0.1)
                    results_container = st.container()

                    # Prepare the arguments for the run_cellchat_for_group function
                    args_group1 = (
                        adata1,
                        group1,
                        gene_use,
                        complex_input,
                        cofactor_input,
                        resource,        
                        celltype_key,
                        data_layer,
                        min_cells,
                        expr_prop,
                        n_perms,
                        n_cpus,
                        union_genes if use_union_genes in ["union", "intersection"] else None,
                        r_patcher,
                        population_size,
                        fun_type,
                        trim
                    )
                    args_group2 = (
                        adata2,
                        group2,
                        gene_use,
                        complex_input,
                        cofactor_input,
                        resource,
                        celltype_key,
                        data_layer,
                        min_cells,
                        expr_prop,
                        n_perms,
                        n_cpus,
                        union_genes if use_union_genes in ["union", "intersection"] else None,
                        r_patcher,
                        population_size,
                        fun_type,
                        trim
                    )

                    # Sequential execution instead of parallel to avoid rpy2 issues
                    try:
                        progress.progress(0.25)
                        st.info(f"Starting CellChat analysis for {group1}...")
                        result1 = run_cellchat_for_group(*args_group1)
                        st.success(f"CellChat analysis for {group1} completed!")
                        
                        progress.progress(0.75)
                        st.info(f"Starting CellChat analysis for {group2}...")
                        result2 = run_cellchat_for_group(*args_group2)
                        st.success(f"CellChat analysis for {group2} completed!")
                        
                    except Exception as e:
                        st.error(f"Error during CellChat analysis: {str(e)}")
                        st.stop()



                    # Now that both analyses are done, we can do the rest of your logic
                    progress.progress(80)
                    st.success("Both CellChat analyses have finished!")

                    # Store results in session state for visualization
                    st.session_state.result1 = result1
                    st.session_state.result2 = result2
                    st.session_state.group1_name = group1
                    st.session_state.group2_name = group2
                    
                    # Combine the results into a single comparison object for visualization
                    combined_results = {
                        'group1': {
                            'name': group1,
                            'result': result1
                        },
                        'group2': {
                            'name': group2,
                            'result': result2
                        },
                        'celltype_key': celltype_key,
                        'splitby_key': splitby_key
                    }
                    
                    st.session_state.combined_results = combined_results
                    
                    progress.progress(100)
                    
                    # Show success message
                    st.success("Comparison analysis completed successfully! Proceed to the visualization tab.")
                    
                    # Add download buttons for saving results
                    st.markdown("#### 💾 Save Analysis Results")
                    
                    # Save individual results as PKL
                    col1, col2 = st.columns(2)
                    
                    # Prepare data for download
                    serialized_results, filename = save_cellchat_results(
                        result1, result2, group1, group2,
                        uploaded_filename=uploaded_file.name,
                        selected_types=selected_types
                    )
                    
                    if serialized_results:
                        st.download_button(
                            label="Download Combined Results",
                            data=serialized_results,
                            file_name=filename,
                            mime="application/octet-stream",
                            help="Download both group analyses as a single file for later import"
                        )

                    # When analysis succeeds
                    st.session_state.completed_analysis = True
                    st.session_state.combined_results = combined_results
                    
                    # Success message and guide to next steps
                    st.success("Comparison analysis completed successfully!")
                    st.info("Please select 'Load saved results for visualization' from the mode selector above to view the results.")

                    
                except Exception as e:
                    st.error(f"Error during comparison analysis: {str(e)}")
                    st.exception(e)
    

    # Mode 2: Load saved results for visualization
    else:  # analysis_mode == "Load saved results for visualization"
        st.header("Load Saved CellChat Results")
        
        load_mode = st.radio(
            "Choose how to load results:",
            ["Upload combined results file", "Upload two separate result files"],
            index=0
        )
        
        if load_mode == "Upload combined results file":
            uploaded_results = st.file_uploader("Upload saved CellChat results file", type=['pkl'], 
                                             help="Upload a previously saved CellChat results file (.pkl)")
            
            if uploaded_results is not None:
                try:
                    with st.spinner('Loading saved CellChat results...'):
                        # Read the pickle file
                        bytes_data = uploaded_results.read()
                        loaded_results = pickle.loads(bytes_data)
                        
                        # Validate the loaded results
                        if isinstance(loaded_results, dict) and 'group1' in loaded_results and 'group2' in loaded_results:
                            st.success(f"Successfully loaded results: {uploaded_results.name}")
                            
                            # Extract necessary information
                            group1_name = loaded_results['group1']['name']
                            group2_name = loaded_results['group2']['name']
                            result1 = loaded_results['group1']['result']
                            result2 = loaded_results['group2']['result']
                            
                            # Store in session state
                            st.session_state.result1 = result1
                            st.session_state.result2 = result2
                            st.session_state.group1_name = group1_name
                            st.session_state.group2_name = group2_name
                            
                            # Store the combined results
                            st.session_state.combined_results = {
                                'group1': {
                                    'name': group1_name,
                                    'result': result1
                                },
                                'group2': {
                                    'name': group2_name,
                                    'result': result2
                                }
                            }
                            
                            # Display summary information
                            # ... (existing summary display code)
                        else:
                            st.error("The uploaded file does not contain valid CellChat results. Please upload a properly formatted results file.")
                
                except Exception as e:
                    st.error(f"Error loading results: {str(e)}")
                    import traceback
                    st.error(traceback.format_exc())
        
        else:  # "Upload two separate results files"
            # Create two columns for file uploaders
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("### Group 1")
                group1_name = st.text_input("Group 1 name:", value="Group1")
                uploaded_file1 = st.file_uploader("Upload Group 1 results file", type=['pkl'], 
                                                help="Upload CellChat results file for Group 1 (.pkl)")
            
            with col2:
                st.markdown("### Group 2")
                group2_name = st.text_input("Group 2 name:", value="Group2")
                uploaded_file2 = st.file_uploader("Upload Group 2 results file", type=['pkl'], 
                                                help="Upload CellChat results file for Group 2 (.pkl)")
            
            if uploaded_file1 is not None and uploaded_file2 is not None:
                try:
                    with st.spinner('Loading separate CellChat results...'):
                        # Read the first file
                        bytes_data1 = uploaded_file1.read()
                        result1 = pickle.loads(bytes_data1)
                        
                        # Read the second file
                        bytes_data2 = uploaded_file2.read()
                        result2 = pickle.loads(bytes_data2)
                        
                        # Validate the results (basic validation for CellChat results)
                        if isinstance(result1, dict) and 'results' in result1 and isinstance(result2, dict) and 'results' in result2:
                            st.success("Successfully loaded both result files")
                            
                            # Store in session state
                            st.session_state.result1 = result1
                            st.session_state.result2 = result2
                            st.session_state.group1_name = group1_name
                            st.session_state.group2_name = group2_name
                            
                            # Store the combined results
                            st.session_state.combined_results = {
                                'group1': {
                                    'name': group1_name,
                                    'result': result1
                                },
                                'group2': {
                                    'name': group2_name,
                                    'result': result2
                                }
                            }
                            
                            # Display summary information
                            st.subheader("Loaded Results Summary")
                            
                            col1, col2 = st.columns(2)
                            
                            with col1:
                                st.write(f"### Group 1: {group1_name}")
                                if 'results' in result1 and not result1['results'].empty:
                                    st.write(f"Significant interactions: {len(result1['results'])}")
                                    
                                    # Display cell types
                                    if 'net' in result1 and 'dimnames' in result1['net'] and len(result1['net']['dimnames']) > 0:
                                        cell_types = result1['net']['dimnames'][0]
                                        st.write(f"Cell types: {', '.join(cell_types)}")
                                    
                                    # Display pathways
                                    if 'netP' in result1 and 'pathways' in result1['netP']:
                                        st.write(f"Pathways: {len(result1['netP']['pathways'])}")
                            
                            with col2:
                                st.write(f"### Group 2: {group2_name}")
                                if 'results' in result2 and not result2['results'].empty:
                                    st.write(f"Significant interactions: {len(result2['results'])}")
                                    
                                    # Display cell types
                                    if 'net' in result2 and 'dimnames' in result2['net'] and len(result2['net']['dimnames']) > 0:
                                        cell_types = result2['net']['dimnames'][0]
                                        st.write(f"Cell types: {', '.join(cell_types)}")
                                    
                                    # Display pathways
                                    if 'netP' in result2 and 'pathways' in result2['netP']:
                                        st.write(f"Pathways: {len(result2['netP']['pathways'])}")
                                        
                            # Add button to continue to visualization
                            if st.button("Proceed to Visualization"):
                                print("rerun")
                        else:
                            st.error("One or both of the uploaded files do not contain valid CellChat results. Please upload properly formatted results files.")
                
                except Exception as e:
                    st.error(f"Error loading results: {str(e)}")
                    import traceback
                    st.error(traceback.format_exc())




    # Show visualization options if results are available
    if 'combined_results' in st.session_state:
        st.header("Visualizations")
        st.sidebar.title("Options")
        st.markdown("##### Options are displayed at the bottom of the left side panel")        
        # Create tabs for different visualization types
        tabs = st.tabs([
            "Comparison Overview", 
            "Circle Plots", 
            "Heatmaps", 
            "Pathway Analysis",
            "LR Contribution",
            "Cell Signaling Roles",  # New tab
            "Signaling Changes"      # New tab
        ])
        
        combined_results = st.session_state.combined_results
        result1 = combined_results['group1']['result']
        result2 = combined_results['group2']['result']
        group1_name = combined_results['group1']['name']
        group2_name = combined_results['group2']['name']


        with st.sidebar:
            colormap_options = [
                "YlOrRd", "OrRd", "YlOrBr", "Oranges", "Reds", "RdPu", "Purples", 
                "PuRd", "Blues", "Greens", "YlGn", "YlGnBu", "GnBu", "Greys", "binary",
                "viridis", "plasma", "inferno", "magma", "cividis", 
                "viridis_r", "plasma_r", "inferno_r", "magma_r", "cividis_r"
            ]
            color_heatmap = st.selectbox(
                "Heatmap color:",
                colormap_options,
                index=2
            )
            cmap_name = st.sidebar.selectbox(
                "Select colormap for other plots:",
                ["tab10","Set1", "Set2", "Set3", "tab20", "Paired", "Dark2",
                "tab20b", "tab20c","Pastel1",
                "Pastel2",  "Accent", "viridis", "plasma", "inferno", "magma"],
                index=0
            )
            # Get cell types from results for filtering options
            if result1 and 'net' in result1 and 'dimnames' in result1['net']:
                cell_list = list(result1['net']['dimnames'][0])
            elif result2 and 'net' in result2 and 'dimnames' in result2['net']:
                cell_list = list(result2['net']['dimnames'][0])
            else:
                cell_list = []
            
            # Add after colormap selection section
            if cell_list:
                # If colormap changed or first time
                if 'cell_color_map' not in st.session_state or st.session_state.get('current_cmap', '') != cmap_name:
                    st.session_state.cell_color_map = create_cell_color_mapping(cell_list, cmap_name)
                    st.session_state.current_cmap = cmap_name

            sort_cell = st.checkbox("Change cell order?")
            if sort_cell and cell_list:
                with st.form("Sorter"):
                    sorted_order = sort_items(cell_list.copy())
                    submitted_sort = st.form_submit_button("Done sorting")
            else:
                sorted_order = None
            
            # Add source and target filtering options
            if cell_list:
                st.markdown("#### Cell Type Filtering")
                filter_cells = st.checkbox("Filter source/target cells?")
                if filter_cells:
                    sources_use = st.multiselect(
                        "Select source cell types:",
                        options=cell_list,
                        default=cell_list,
                        help="Select which cell types to use as sources (signal senders)"
                    )
                    targets_use = st.multiselect(
                        "Select target cell types:",
                        options=cell_list,
                        default=cell_list,
                        help="Select which cell types to use as targets (signal receivers)"
                    )
                else:
                    sources_use = None
                    targets_use = None
            else:
                sources_use = None
                targets_use = None
        
        # Tab 1: Comparison Overview
        with tabs[0]:
            st.subheader("Comparison Overview")
            
            # Display basic statistics
            col1, col2 = st.columns(2)
            
            # Function to get basic stats from a result
            def get_result_stats(result):
                stats = {}
                if 'results' in result and not result['results'].empty:
                    stats['total_interactions'] = len(result['results'])
                    stats['significant_interactions'] = len(result['results'][result['results']['pval'] <= 0.05])
                else:
                    stats['total_interactions'] = 0
                    stats['significant_interactions'] = 0
                
                if 'netP' in result and 'pathways' in result['netP']:
                    stats['pathways'] = len(result['netP']['pathways'])
                else:
                    stats['pathways'] = 0
                    
                return stats
            
            # Get stats for both groups
            stats1 = get_result_stats(result1)
            stats2 = get_result_stats(result2)
            
            with col1:
                st.metric(f"{group1_name} - Total Interactions", stats1['total_interactions'])
                st.metric(f"{group1_name} - Significant Interactions", stats1['significant_interactions'])
                st.metric(f"{group1_name} - Pathways", stats1['pathways'])
                
            with col2:
                st.metric(f"{group2_name} - Total Interactions", stats2['total_interactions'])
                st.metric(f"{group2_name} - Significant Interactions", stats2['significant_interactions'])
                st.metric(f"{group2_name} - Pathways", stats2['pathways'])
            
            # Compare pathways between groups
            st.subheader("Pathway Comparison")
            
            pathways1 = set(result1['netP']['pathways']) if 'netP' in result1 and 'pathways' in result1['netP'] else set()
            pathways2 = set(result2['netP']['pathways']) if 'netP' in result2 and 'pathways' in result2['netP'] else set()
            
            common_pathways = pathways1.intersection(pathways2)
            unique_pathways1 = pathways1 - pathways2
            unique_pathways2 = pathways2 - pathways1
            
            st.write(f"Common pathways: {len(common_pathways)}")
            st.write(f"Unique to {group1_name}: {len(unique_pathways1)}")
            st.write(f"Unique to {group2_name}: {len(unique_pathways2)}")
            
            # Display pathway lists in expanders
            with st.expander("View common pathways"):
                st.write(", ".join(sorted(common_pathways)))
            
            with st.expander(f"View pathways unique to {group1_name}"):
                st.write(", ".join(sorted(unique_pathways1)))
            
            with st.expander(f"View pathways unique to {group2_name}"):
                st.write(", ".join(sorted(unique_pathways2)))
            
            # Compare ligand-receptor pairs
            st.subheader("Ligand-Receptor Pair Comparison")
            
            # Function to extract LR pairs
            def get_lr_pairs(result):
                if 'results' not in result or result['results'].empty:
                    return set()
                
                pairs = set()
                if 'ligand' in result['results'].columns and 'receptor' in result['results'].columns:
                    for _, row in result['results'].iterrows():
                        if pd.notna(row['ligand']) and pd.notna(row['receptor']):
                            pairs.add(f"{row['ligand']}-{row['receptor']}")
                
                return pairs
            
            lr_pairs1 = get_lr_pairs(result1)
            lr_pairs2 = get_lr_pairs(result2)
            
            common_lr = lr_pairs1.intersection(lr_pairs2)
            unique_lr1 = lr_pairs1 - lr_pairs2
            unique_lr2 = lr_pairs2 - lr_pairs1
            
            st.write(f"Common L-R pairs: {len(common_lr)}")
            st.write(f"Unique to {group1_name}: {len(unique_lr1)}")
            st.write(f"Unique to {group2_name}: {len(unique_lr2)}")
            
            # Display top LR pairs for each group
            st.subheader("Top Interactions by Group")
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.write(f"### Top interactions in {group1_name}")
                if 'results' in result1 and not result1['results'].empty:
                    top_lr1 = result1['results'].sort_values('prob', ascending=False).head(10)
                    st.dataframe(top_lr1[['source', 'target', 'ligand', 'receptor', 'prob', 'pval']])
                else:
                    st.info("No interaction data available")
            
            with col2:
                st.write(f"### Top interactions in {group2_name}")
                if 'results' in result2 and not result2['results'].empty:
                    top_lr2 = result2['results'].sort_values('prob', ascending=False).head(10)
                    st.dataframe(top_lr2[['source', 'target', 'ligand', 'receptor', 'prob', 'pval']])
                else:
                    st.info("No interaction data available")
        
        # Tab 2: Circle Plots
        with tabs[1]:
            st.subheader("Circle Plots Comparison")
            
            # Choose aggregate or specific pathway
            if ('netP' in result1 and 'pathways' in result1['netP']) or ('netP' in result2 and 'pathways' in result2['netP']):

                col1, col2 = st.columns(2)

                with col1:

                    # Options: "Aggregate" or multiple selection of individual pathways
                    option_type = st.radio(
                        "Pathway selection:",
                        ["Aggregate", "Specific pathways"],
                        horizontal=True
                    )
                    
                    if option_type == "Aggregate":
                        selected_pathway = "Aggregate"
                        measure_key = "weight"
                    else:
                        # Get all pathways from both groups (union)
                        pathways1 = set(result1['netP']['pathways']) if 'netP' in result1 and 'pathways' in result1['netP'] else set()
                        pathways2 = set(result2['netP']['pathways']) if 'netP' in result2 and 'pathways' in result2['netP'] else set()
                        common_pathways = sorted(pathways1.union(pathways2))
                        
                        if common_pathways:
                            selected_pathway = st.multiselect(
                                "Select pathway:",
                                sorted(common_pathways),
                                default=[common_pathways[0]] if common_pathways else []
                            )
                            measure_key = "weight"  # Always use weight for specific pathways
                        else:
                            st.warning("No pathways found in either group")
                            selected_pathway = None
                        
                        measure_key = "weight"

                    circle_type = st.radio(
                        "Plot type:",
                        ["Total network", "Individual cell-type"],
                        horizontal=True
                    )


                with col2:
                    # Display customization options
                    edge_width_max = st.slider("Max edge width:", min_value=1, max_value=20, value=8, step=1)
                    circle_alpha_edge = st.slider("Edge transparency:", min_value=0.1, max_value=1.0, value=0.6, step=0.1)
                    vertex_size_max= st.slider("Max node size:", min_value=1, max_value=15, value=6, step=1)
                    show_vertex = st.checkbox("Node size reflects the number of cells", value=False)

                

            else:
                st.warning("Pathway information not available")
                option_type = "Aggregate"
                selected_pathway = "Aggregate"
                measure_key = "weight"
            
            # Generate circle plots side by side
            if selected_pathway is not None:
                if st.button("Generate circle plots"):

                    plotcol1, plotcol2 = st.columns(2)
                    circle_plot_files = []

                    with plotcol1:
                        st.write(f"### {group1_name}")
                        try:
                            # Title for plot
                            if selected_pathway == "Aggregate":
                                title = f"{group1_name} - Aggregate network"
                            else:
                                title = f"{group1_name} - {'.'.join(selected_pathway)} pathway"
                            
                            with st.spinner("Generating circle plot..."):
                                if circle_type == "Total network":
                                    fig1 = netVisual_circle(
                                        net=result1,
                                        title_name=title,
                                        pathway_name=None if selected_pathway == "Aggregate" else selected_pathway,
                                        measure=measure_key,
                                        cmap_name=cmap_name,
                                        sorted_order=sorted_order,
                                        arrow=True,
                                        edge_width_max=edge_width_max,
                                        alpha_edge=circle_alpha_edge,
                                     #   vertex_weight=vertex_weight, # Not calculated yet
                                        vertex_size_max=vertex_size_max,
                                        color_use=st.session_state.get('cell_color_map', None),  # Add color mapping
                                        sources_use=sources_use,
                                        targets_use=targets_use
                                    )
                                else:
                                    fig1 = netVisual_circle_individual(
                                        net=result1,
                                        pathway_name=None if selected_pathway == "Aggregate" else selected_pathway,
                                        measure=measure_key,
                                        cmap_name=cmap_name,
                                        sorted_order=sorted_order,
                                        arrow=True,
                                        edge_width_max=edge_width_max,
                                        alpha_edge=circle_alpha_edge,
                                     #   vertex_weight=vertex_weight, # Not calculated yet
                                        vertex_size_max=vertex_size_max,
                                        color_use=st.session_state.get('cell_color_map', None)  # Add color mapping
                                    )
                            
                            st.pyplot(fig1)
                            
                            # Save and offer download
                            pdf_path = f"{cellchat_temp_dir}/circle_{group1_name}_{selected_pathway}.pdf"
                            fig1.savefig(pdf_path, bbox_inches='tight')
                            circle_plot_files.append(pdf_path)
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label=f"Download {group1_name} PDF",
                                    data=pdf_bytes,
                                    file_name=f'cellchat_circle_{group1_name}_{selected_pathway}.pdf',
                                    mime='application/octet-stream'
                                )
                        except Exception as e:
                            st.error(f"Error generating circle plot for {group1_name}: {str(e)}")
                    
                    with plotcol2:
                        st.write(f"### {group2_name}")
                        try:
                            # Title for plot
                            if selected_pathway == "Aggregate":
                                title = f"{group1_name} - Aggregate network"
                            else:
                                title = f"{group1_name} - {'.'.join(selected_pathway)} pathway"
                            
                            with st.spinner("Generating circle plot..."):
                                if circle_type == "Total network":
                                    fig2 = netVisual_circle(
                                        net=result2,
                                        title_name=title,
                                        pathway_name=None if selected_pathway == "Aggregate" else selected_pathway,
                                        measure=measure_key,
                                        cmap_name=cmap_name,
                                        sorted_order=sorted_order,
                                        arrow=True,
                                        edge_width_max=edge_width_max,
                                        alpha_edge=circle_alpha_edge,
                                     #   vertex_weight=vertex_weight, # Not calculated yet
                                        vertex_size_max=vertex_size_max,
                                        color_use=st.session_state.get('cell_color_map', None),  # Add color mapping
                                        sources_use=sources_use,
                                        targets_use=targets_use
                                    )
                                else:
                                    fig2 = netVisual_circle_individual(
                                        net=result2,
                                        pathway_name=None if selected_pathway == "Aggregate" else selected_pathway,
                                        measure=measure_key,
                                        title_name=title,
                                        cmap_name=cmap_name,
                                        sorted_order=sorted_order,
                                        arrow=True,
                                        edge_width_max=edge_width_max,
                                        alpha_edge=circle_alpha_edge,
                                      #  vertex_weight=vertex_weight,
                                        vertex_size_max=vertex_size_max,
                                        color_use=st.session_state.get('cell_color_map', None)  # Add color mapping
                                    )
                            
                            st.pyplot(fig2)
                            
                            # Save and offer download
                            pdf_path = f"{cellchat_temp_dir}/circle_{group2_name}_{selected_pathway}.pdf"
                            fig2.savefig(pdf_path, bbox_inches='tight')
                            circle_plot_files.append(pdf_path)
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label=f"Download {group2_name} PDF",
                                    data=pdf_bytes,
                                    file_name=f'cellchat_circle_{group2_name}_{selected_pathway}.pdf',
                                    mime='application/octet-stream'
                                )
                        except Exception as e:
                            st.error(f"Error generating circle plot for {group2_name}: {str(e)}")
                    
                    # Download all plots as ZIP
                    if len(circle_plot_files) > 0:
                        try:
                            zip_buffer = create_zip_from_files(circle_plot_files)
                            if zip_buffer:
                                st.download_button(
                                    label="Download all circle plots as ZIP",
                                    data=zip_buffer,
                                    file_name=f'cellchat_circle_plots_{group1_name}_vs_{group2_name}.zip',
                                    mime='application/zip'
                                )
                        except Exception as e:
                            st.error(f"Error creating ZIP file: {str(e)}")
        
        # Tab 3: Heatmaps
        with tabs[2]:
            st.subheader("Interaction Heatmaps Comparison")
            st.info("""
            This tab shows heatmaps of cell-cell interaction strength or counts for both groups.
            The differential heatmap shows the difference between groups (Group 2 - Group 1).
            Red cells indicate higher values in Group 1, blue cells indicate higher values in Group 2.
            """)
            
            # Heatmap options
            col1, col2 = st.columns(2)
            
            with col1:
                
                heatmap_type = st.radio(
                    "Heatmap type:",
                    ["Interaction strength", "Interaction number"],
                    horizontal=True
                )
                
                heatmap_annot = st.checkbox("Show values?", value=True)
                if 'net' in result1 and 'dimnames' in result1['net'] and 'net' in result2 and 'dimnames' in result2['net']:
                    cell_types1 = set(result1['net']['dimnames'][0])
                    cell_types2 = set(result2['net']['dimnames'][0])
                    common_cell_types = sorted(cell_types1.intersection(cell_types2))
                    
                    if sort_cell and sorted_order:
                        # Filter sorted_order to include only common cell types
                        use_sorted_order = [ct for ct in sorted_order if ct in common_cell_types]
                    else:
                        use_sorted_order = common_cell_types
                else:
                    use_sorted_order = None

                reverse_comparison = st.checkbox(
                    "Reverse comparison?",
                    value=False,
                    help="When checked, the differential heatmap will show Control - Test instead of Test - Control"
                )
                
                # Pathway selection for heatmap
                st.markdown("#### Pathway Selection")
                if ('netP' in result1 and 'pathways' in result1['netP']) or ('netP' in result2 and 'pathways' in result2['netP']):
                    pathways1 = set(result1['netP']['pathways']) if 'netP' in result1 and 'pathways' in result1['netP'] else set()
                    pathways2 = set(result2['netP']['pathways']) if 'netP' in result2 and 'pathways' in result2['netP'] else set()
                    common_pathways = sorted(pathways1.union(pathways2))
                    
                    pathway_option = st.radio(
                        "Heatmap scope:",
                        ["All pathways", "Specific pathways"],
                        horizontal=True,
                        help="Choose to show aggregate of all pathways or filter by specific pathways"
                    )
                    
                    if pathway_option == "Specific pathways" and common_pathways:
                        selected_pathways = st.multiselect(
                            "Select pathways:",
                            options=common_pathways,
                            default=common_pathways[:min(3, len(common_pathways))],
                            help="Select specific pathways to include in the heatmap"
                        )
                    else:
                        selected_pathways = None
                else:
                    st.warning("Pathway information not available")
                    pathway_option = "All pathways"
                    selected_pathways = None

            with col2:
                heatmap_title = st.slider("Title font size:", min_value=0, max_value=30, value=20)
                heatmap_font = st.slider("Font size:", min_value=0, max_value=30, value=14)
                heatmap_x = st.slider("Fig width:", min_value=1.0, max_value=20.0, value=10.0)
                heatmap_y = st.slider("Fig height:", min_value=1.0, max_value=20.0, value=8.0)
                    
            # Generate heatmap button
            if st.button("Generate heatmaps"):
                heatmap_files = []
                
                # Get data matrices and calculate maximum values
                with st.spinner("Computing unified scale..."):
                    # Get matrices based on selected heatmap type
                    if heatmap_type == "Interaction strength":
                        matrix1 = result1['network']['strength_matrix'] 
                        matrix2 = result2['network']['strength_matrix']
                    else:  # "Interaction number"
                        matrix1 = result1['network']['count_matrix']
                        matrix2 = result2['network']['count_matrix']
                    
                    # Align to common cell types
                    if use_sorted_order:
                        matrix1 = matrix1.reindex(index=use_sorted_order, columns=use_sorted_order)
                        matrix2 = matrix2.reindex(index=use_sorted_order, columns=use_sorted_order)
                    
                    # Calculate global maximum value
                    global_max = max(matrix1.values.max(), matrix2.values.max())
                    st.write(f"Max value: {global_max:.2f}")
                
                # Generate heatmap using modified function
                with st.spinner("Generating heatmaps with unified scale..."):
                    # Configure measure
                    measure = "weight" if heatmap_type == "Interaction strength" else "count"
                    
                    # Pathway filtering: create filtered result if pathways are selected
                    if selected_pathways:
                        # Filter result1 by selected pathways
                        filtered_result1 = filter_result_by_pathways(result1, selected_pathways)
                        filtered_result2 = filter_result_by_pathways(result2, selected_pathways)
                        
                        # Re-calculate global max for filtered data
                        if heatmap_type == "Interaction strength":
                            matrix1_filtered = filtered_result1['network']['strength_matrix'] 
                            matrix2_filtered = filtered_result2['network']['strength_matrix']
                        else:
                            matrix1_filtered = filtered_result1['network']['count_matrix']
                            matrix2_filtered = filtered_result2['network']['count_matrix']
                        
                        global_max_filtered = max(matrix1_filtered.values.max(), matrix2_filtered.values.max())
                        
                        # Use filtered results and max
                        result1_to_use = filtered_result1
                        result2_to_use = filtered_result2
                        max_to_use = global_max_filtered
                        
                        pathway_str = ', '.join(selected_pathways[:3])  # Show first 3 pathways in title
                        if len(selected_pathways) > 3:
                            pathway_str += f" (+{len(selected_pathways)-3} more)"
                        title1 = f"{group1_name} - {heatmap_type} ({pathway_str})"
                        title2 = f"{group2_name} - {heatmap_type} ({pathway_str})"
                    else:
                        # Use original results
                        result1_to_use = result1
                        result2_to_use = result2
                        max_to_use = global_max
                        title1 = f"{group1_name} - {heatmap_type}"
                        title2 = f"{group2_name} - {heatmap_type}"
                    
                    # Generate first heatmap
                    fig_heatmap1 = netVisual_heatmap(
                        result1_to_use, 
                        measure=measure, 
                        color_heatmap=color_heatmap, 
                        font_size=heatmap_font, 
                        font_size_title=heatmap_title,
                        sorted_order=use_sorted_order,
                        title=title1,
                        annot=heatmap_annot,
                        vmin=0,            # Explicitly set minimum value
                        vmax=max_to_use,    # Set common maximum value
                        figx=heatmap_x,
                        figy=heatmap_y,
                        sources_use=sources_use if filter_cells else None,
                        targets_use=targets_use if filter_cells else None
                    )
                    
                    # Generate second heatmap
                    fig_heatmap2 = netVisual_heatmap(
                        result2_to_use, 
                        measure=measure, 
                        color_heatmap=color_heatmap, 
                        font_size=heatmap_font, 
                        font_size_title=heatmap_title,
                        sorted_order=use_sorted_order,
                        title=title2,
                        annot=heatmap_annot,
                        vmin=0,            # Explicitly set minimum value
                        vmax=max_to_use,    # Set common maximum value
                        figx=heatmap_x,
                        figy=heatmap_y,
                        sources_use=sources_use if filter_cells else None,
                        targets_use=targets_use if filter_cells else None
                    )
                
                # Display each heatmap
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"### {group1_name}")
                    st.pyplot(fig_heatmap1)
                    
                    # Save PDF and download
                    # Create pathway suffix for filename
                    if selected_pathways:
                        pathway_filename = '_'.join(selected_pathways[:2])  # Use first 2 pathways
                        if len(selected_pathways) > 2:
                            pathway_filename += f"_plus{len(selected_pathways)-2}"
                        filename_suffix = f"_{pathway_filename}"
                    else:
                        filename_suffix = ""
                    
                    pdf_path = f"{cellchat_temp_dir}/heatmap_{group1_name}{filename_suffix}.pdf"
                    fig_heatmap1.savefig(pdf_path, bbox_inches='tight')
                    heatmap_files.append(pdf_path)
                    
                    with open(pdf_path, "rb") as pdf_file:
                        pdf_bytes = pdf_file.read()
                        st.download_button(
                            label=f"Download {group1_name} heatmap",
                            data=pdf_bytes,
                            file_name=f'cellchat_{group1_name}_heatmap{filename_suffix}.pdf',
                            mime='application/octet-stream'
                        )
                
                with col2:
                    st.write(f"### {group2_name}")
                    st.pyplot(fig_heatmap2)
                    
                    # Save PDF and download
                    pdf_path = f"{cellchat_temp_dir}/heatmap_{group2_name}{filename_suffix}.pdf"
                    fig_heatmap2.savefig(pdf_path, bbox_inches='tight')
                    heatmap_files.append(pdf_path)
                    
                    with open(pdf_path, "rb") as pdf_file:
                        pdf_bytes = pdf_file.read()
                        st.download_button(
                            label=f"Download {group2_name} heatmap",
                            data=pdf_bytes,
                            file_name=f'cellchat_{group2_name}_heatmap{filename_suffix}.pdf',
                            mime='application/octet-stream'
                        )
                
                # Generate difference heatmap
                if 'network' in result1 and 'strength_matrix' in result1['network'] and 'network' in result2 and 'strength_matrix' in result2['network']:
                    st.subheader("Differential Heatmap")
                    
                    try:
                        with st.spinner("Generating differential heatmap..."):
                            # Use filtered or original matrices based on pathway selection
                            if selected_pathways:
                                # Use the already filtered matrices
                                filtered_matrix1 = matrix1_filtered.copy()
                                filtered_matrix2 = matrix2_filtered.copy()
                            else:
                                # Use original matrices
                                filtered_matrix1 = matrix1.copy()
                                filtered_matrix2 = matrix2.copy()
                            
                            # Apply source/target filtering
                            if filter_cells and sources_use and targets_use:
                                # Keep only cells specified in sources_use/targets_use
                                valid_sources = [s for s in sources_use if s in filtered_matrix1.index and s in filtered_matrix1.columns]
                                valid_targets = [t for t in targets_use if t in filtered_matrix1.columns and t in filtered_matrix1.index]
                                
                                if valid_sources and valid_targets:
                                    filtered_matrix1 = filtered_matrix1.loc[valid_sources, valid_targets]
                                    filtered_matrix2 = filtered_matrix2.loc[valid_sources, valid_targets]
                            
                            # Differential heatmap title with pathway info
                            if selected_pathways:
                                pathway_str = ', '.join(selected_pathways[:2])  # Show first 2 pathways in diff title
                                if len(selected_pathways) > 2:
                                    pathway_str += f" (+{len(selected_pathways)-2} more)"
                                pathway_suffix = f" ({pathway_str})"
                            else:
                                pathway_suffix = ""
                            
                            # Calculate difference matrix based on comparison direction
                            if reverse_comparison:
                                diff_matrix = filtered_matrix1 - filtered_matrix2
                                title = f"Differential {heatmap_type}: {group1_name} - {group2_name}{pathway_suffix}"
                            else:
                                diff_matrix = filtered_matrix2 - filtered_matrix1
                                title = f"Differential {heatmap_type}: {group2_name} - {group1_name}{pathway_suffix}"
                            
                            # Get symmetric boundary values
                            max_abs_val = max(abs(diff_matrix.min().min()), abs(diff_matrix.max().max()))
                            vmin = -max_abs_val
                            vmax = max_abs_val
                            
                            # Create figure
                            fig_diff, ax = plt.subplots(figsize=(10, 8))
                            fmt = ".2f" if measure == "weight" else ".0f"
                            
                            # Plot heatmap with diverging colormap
                            sns.heatmap(
                                diff_matrix,
                                cmap="RdBu_r",  # Diverging colormap: blue for positive, red for negative
                                center=0,
                                vmin=vmin,
                                vmax=vmax,
                                square=True,
                                linewidths=0.5,
                                annot=heatmap_annot,
                                fmt=fmt,
                                ax=ax
                            )
                            
                            ax.set_title(title, fontsize=heatmap_title)
                            ax.set_xlabel('Targets (Receiver)', fontsize=heatmap_font+2)
                            ax.set_ylabel('Sources (Sender)', fontsize=heatmap_font+2)
                            plt.xticks(rotation=45, ha='right',fontsize=heatmap_font)
                            plt.yticks(rotation=0, fontsize=heatmap_font)
                            fig_diff.tight_layout()

                            # Adjust colorbar after tight_layout
                            cbar_ax = fig_diff.axes[-1]
                            current_pos = cbar_ax.get_position()
                            new_pos = [current_pos.x0, current_pos.y0, 0.01, 0.2]  # Change width to 0.02
                            cbar_ax.set_position(new_pos)
                            
                            # Display plot
                            st.pyplot(fig_diff)
                            
                            # Save PDF and download
                            pdf_path = f"{cellchat_temp_dir}/diff_heatmap{filename_suffix}.pdf"
                            fig_diff.savefig(pdf_path, bbox_inches='tight')
                            heatmap_files.append(pdf_path)
                            
                            with open(pdf_path, "rb") as pdf_file:
                                pdf_bytes = pdf_file.read()
                                st.download_button(
                                    label="Download differential heatmap",
                                    data=pdf_bytes,
                                    file_name=f'cellchat_diff_heatmap{filename_suffix}.pdf',
                                    mime='application/octet-stream'
                                )
                    except Exception as e:
                        st.error(f"Error occurred while generating difference heatmap: {str(e)}")
                        st.exception(e)
                
                # Download all heatmaps as ZIP
                if len(heatmap_files) > 0:
                    try:
                        zip_buffer = create_zip_from_files(heatmap_files)
                        
                        if zip_buffer:
                            st.download_button(
                                label="Download all heatmaps as ZIP",
                                data=zip_buffer,
                                file_name=f'cellchat_all_heatmaps_{group1_name}_vs_{group2_name}.zip',
                                mime='application/zip'
                            )
                    except Exception as e:
                        st.error(f"Error occurred while creating ZIP file: {str(e)}")

        # Tab 4: Pathway Analysis
        with tabs[3]:
            st.subheader("Pathway Analysis Comparison")
            
            def get_pathway_contribution(result):
                if 'netP' not in result:
                    return pd.DataFrame()
                
                if 'pathways' not in result['netP']:
                    return pd.DataFrame()
                
                # Check if pathways is empty (different ways to check depending on type)
                pathways = result['netP']['pathways']
                if isinstance(pathways, (list, np.ndarray)):
                    if len(pathways) == 0:
                        return pd.DataFrame()
                elif isinstance(pathways, pd.Series):
                    if pathways.empty:
                        return pd.DataFrame()
                elif pathways is None:
                    return pd.DataFrame()
                
                # Get pathway strength
                pathway_strength = []
                for i, pathway in enumerate(pathways):
                    # Sum the strength for this pathway across all cell pairs
                    strength = np.sum(result['netP']['prob'][:, :, i])
                    pathway_strength.append({"pathway": pathway, "strength": strength})
                
                # Convert to DataFrame and sort
                df = pd.DataFrame(pathway_strength)
                if not df.empty:
                    df = df.sort_values("strength", ascending=False)
                
                return df
            
            # Get pathway contributions for both groups
            pathway_contrib1 = get_pathway_contribution(result1)
            pathway_contrib2 = get_pathway_contribution(result2)
            
            # Display pathway contributions
            if not pathway_contrib1.empty and not pathway_contrib2.empty:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"### Top Pathways in {group1_name}")
                    st.dataframe(pathway_contrib1.head(10))
                    
                    # Plot top pathways
                    fig1, ax1 = plt.subplots(figsize=(10, 6))
                    top_n = min(10, len(pathway_contrib1))
                    
                    if top_n > 0:
                        top_pathways1 = pathway_contrib1.head(top_n)
                        bars1 = ax1.barh(top_pathways1['pathway'], top_pathways1['strength'], color='skyblue')
                        ax1.set_xlabel('Interaction Strength')
                        ax1.set_title(f'Top Pathways in {group1_name}')
                        plt.tight_layout()
                        
                        st.pyplot(fig1)
                
                with col2:
                    st.write(f"### Top Pathways in {group2_name}")
                    st.dataframe(pathway_contrib2.head(10))
                    
                    # Plot top pathways
                    fig2, ax2 = plt.subplots(figsize=(10, 6))
                    top_n = min(10, len(pathway_contrib2))
                    
                    if top_n > 0:
                        top_pathways2 = pathway_contrib2.head(top_n)
                        bars2 = ax2.barh(top_pathways2['pathway'], top_pathways2['strength'], color='lightcoral')
                        ax2.set_xlabel('Interaction Strength')
                        ax2.set_title(f'Top Pathways in {group2_name}')
                        plt.tight_layout()
                        
                        st.pyplot(fig2)
                
                # Comparative pathway analysis
                st.subheader("Comparative Pathway Analysis")
                
                # Merge the pathway data
                merged_pathway = pd.merge(
                    pathway_contrib1, pathway_contrib2,
                    on='pathway', how='outer', suffixes=('_' + group1_name, '_' + group2_name)
                )
                
                # Fill NaN values with 0
                merged_pathway = merged_pathway.fillna(0)
                
                # Calculate fold change and sort
                strength1_col = f"strength_{group1_name}"
                strength2_col = f"strength_{group2_name}"
                
                # Add a small pseudocount to avoid division by zero
                epsilon = 1e-10
                merged_pathway['log2_fold_change'] = np.log2((merged_pathway[strength2_col] + epsilon) / (merged_pathway[strength1_col] + epsilon))
                
                # Sort by absolute fold change
                merged_pathway['abs_log2_fold_change'] = np.abs(merged_pathway['log2_fold_change'])
                merged_pathway = merged_pathway.sort_values('abs_log2_fold_change', ascending=False)
                
                # Display the comparative table
                st.dataframe(merged_pathway[[
                    'pathway', strength1_col, strength2_col, 'log2_fold_change'
                ]].head(15))
                
                # Plot the fold changes
                try:
                    fig_fc, ax_fc = plt.subplots(figsize=(12, 8))
                    top_n_fc = min(15, len(merged_pathway))
                    
                    if top_n_fc > 0:
                        top_pathways_fc = merged_pathway.head(top_n_fc)
                        
                        # Create bar colors based on fold change direction
                        bar_colors = ['red' if x < 0 else 'blue' for x in top_pathways_fc['log2_fold_change']]
                        
                        bars_fc = ax_fc.barh(
                            top_pathways_fc['pathway'],
                            top_pathways_fc['log2_fold_change'],
                            color=bar_colors
                        )
                        
                        # Add a vertical line at x=0
                        ax_fc.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
                        
                        ax_fc.set_xlabel('Log2 Fold Change')
                        ax_fc.set_title(f'Pathway Fold Change: {group2_name} vs {group1_name}')
                        
                        # Add legends
                        import matplotlib.patches as mpatches
                        red_patch = mpatches.Patch(color='red', label=f'Higher in {group1_name}')
                        blue_patch = mpatches.Patch(color='blue', label=f'Higher in {group2_name}')
                        ax_fc.legend(handles=[red_patch, blue_patch])
                        
                        plt.tight_layout()
                        st.pyplot(fig_fc)
                        
                        # Save the fold change plot
                        fc_pdf_path = f"{cellchat_temp_dir}/pathway_fold_change.pdf"
                        fig_fc.savefig(fc_pdf_path, bbox_inches='tight')
                        
                        # Add download button
                        with open(fc_pdf_path, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                            st.download_button(
                                label="Download Pathway Fold Change Plot",
                                data=pdf_bytes,
                                file_name=f'cellchat_pathway_fold_change_{group1_name}_vs_{group2_name}.pdf',
                                mime='application/octet-stream'
                            )
                except Exception as e:
                    st.error(f"Error generating fold change plot: {str(e)}")
            else:
                st.warning("Pathway information not available for one or both groups")
        
        # Tab 5: LR Contribution
        with tabs[4]:
            st.subheader("Ligand-Receptor Contribution Analysis")
            
            # Function to get top LR pairs
            def get_top_lr_pairs(result, top_n=10):
                if 'results' not in result or result['results'].empty:
                    return pd.DataFrame()
                
                # Copy and sort by probability
                lr_df = result['results'].copy()
                lr_df = lr_df.sort_values('prob', ascending=False)
                
                # Take top N
                top_lr = lr_df.head(top_n)
                
                return top_lr
            
            # Get top LR pairs for both groups
            top_lr1 = get_top_lr_pairs(result1, top_n=15)
            top_lr2 = get_top_lr_pairs(result2, top_n=15)
            
            if not top_lr1.empty and not top_lr2.empty:
                # Display top LR pairs
                col1, col2 = st.columns(2)
                
                with col1:
                    st.write(f"### Top L-R Pairs in {group1_name}")
                    
                    # Create a more readable format for display
                    display_cols1 = top_lr1[['source', 'target', 'ligand', 'receptor', 'prob']].copy()
                    display_cols1['L-R Pair'] = display_cols1['ligand'] + ' → ' + display_cols1['receptor']
                    display_cols1['Communication'] = display_cols1['source'] + ' → ' + display_cols1['target']
                    
                    st.dataframe(display_cols1[['Communication', 'L-R Pair', 'prob']])
                    
                    # Use the lr_contribution visualization
                    try:
                        from pages.cellchat_vis import plot_lr_contribution
                        
                        # Prepare network data like in cellchat.py
                        network_data1 = result1['network'].copy()
                        
                        # Add interaction names and LR pairs information
                        if 'results' in result1 and not result1['results'].empty:
                            results_df = result1['results']
                            interaction_names = results_df['interaction_name'].tolist()
                            lr_pairs = {}
                            for _, row in results_df.iterrows():
                                lr_pairs[row['interaction_name']] = f"{row['ligand']}-{row['receptor']}"
                            
                            network_data1['lr_names'] = lr_pairs
                            network_data1['interaction_names'] = interaction_names
                        
                        fig_lr1 = plot_lr_contribution(
                            network_data1,
                            top_n=10,
                            figsize=(10, 8)
                        )
                        st.pyplot(fig_lr1)
                        
                        # Save figure
                        lr_pdf_path1 = f"{cellchat_temp_dir}/lr_contribution_{group1_name}.pdf"
                        fig_lr1.savefig(lr_pdf_path1, bbox_inches='tight')
                        
                        with open(lr_pdf_path1, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                            st.download_button(
                                label=f"Download LR Contribution ({group1_name})",
                                data=pdf_bytes,
                                file_name=f'cellchat_lr_contribution_{group1_name}.pdf',
                                mime='application/octet-stream'
                            )
                    except Exception as e:
                        st.error(f"Error generating LR contribution plot for {group1_name}: {str(e)}")
                
                with col2:
                    st.write(f"### Top L-R Pairs in {group2_name}")
                    
                    # Create a more readable format for display
                    display_cols2 = top_lr2[['source', 'target', 'ligand', 'receptor', 'prob']].copy()
                    display_cols2['L-R Pair'] = display_cols2['ligand'] + ' → ' + display_cols2['receptor']
                    display_cols2['Communication'] = display_cols2['source'] + ' → ' + display_cols2['target']
                    
                    st.dataframe(display_cols2[['Communication', 'L-R Pair', 'prob']])
                    
                    # Use the lr_contribution visualization
                    try:
                        from pages.cellchat_vis import plot_lr_contribution
                        
                        # Prepare network data like in cellchat.py
                        network_data2 = result2['network'].copy()
                        
                        # Add interaction names and LR pairs information
                        if 'results' in result2 and not result2['results'].empty:
                            results_df = result2['results']
                            interaction_names = results_df['interaction_name'].tolist()
                            lr_pairs = {}
                            for _, row in results_df.iterrows():
                                lr_pairs[row['interaction_name']] = f"{row['ligand']}-{row['receptor']}"
                            
                            network_data2['lr_names'] = lr_pairs
                            network_data2['interaction_names'] = interaction_names
                        
                        fig_lr2 = plot_lr_contribution(
                            network_data2,
                            top_n=10,
                            figsize=(10, 8)
                        )
                        st.pyplot(fig_lr2)
                        
                        # Save figure
                        lr_pdf_path2 = f"{cellchat_temp_dir}/lr_contribution_{group2_name}.pdf"
                        fig_lr2.savefig(lr_pdf_path2, bbox_inches='tight')
                        
                        with open(lr_pdf_path2, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                            st.download_button(
                                label=f"Download LR Contribution ({group2_name})",
                                data=pdf_bytes,
                                file_name=f'cellchat_lr_contribution_{group2_name}.pdf',
                                mime='application/octet-stream'
                            )
                    except Exception as e:
                        st.error(f"Error generating LR contribution plot for {group2_name}: {str(e)}")
                
                # Comparative LR analysis
                st.subheader("Common L-R Pairs Between Groups")
                
                # Find common LR pairs
                lr_pairs1 = set(top_lr1['ligand'] + '_' + top_lr1['receptor'])
                lr_pairs2 = set(top_lr2['ligand'] + '_' + top_lr2['receptor'])
                common_pairs = lr_pairs1.intersection(lr_pairs2)
                
                if common_pairs:
                    # Get data for common pairs
                    common_data = []
                    
                    for pair in common_pairs:
                        ligand, receptor = pair.split('_', 1)
                        
                        # Get data from group 1
                        row1 = top_lr1[(top_lr1['ligand'] == ligand) & (top_lr1['receptor'] == receptor)]
                        if not row1.empty:
                            row1 = row1.iloc[0]
                            prob1 = row1['prob']
                            source1 = row1['source']
                            target1 = row1['target']
                        else:
                            # Find in the full results
                            row1 = result1['results'][(result1['results']['ligand'] == ligand) & 
                                                     (result1['results']['receptor'] == receptor)]
                            if not row1.empty:
                                row1 = row1.sort_values('prob', ascending=False).iloc[0]
                                prob1 = row1['prob']
                                source1 = row1['source']
                                target1 = row1['target']
                            else:
                                prob1 = 0
                                source1 = ""
                                target1 = ""
                        
                        # Get data from group 2
                        row2 = top_lr2[(top_lr2['ligand'] == ligand) & (top_lr2['receptor'] == receptor)]
                        if not row2.empty:
                            row2 = row2.iloc[0]
                            prob2 = row2['prob']
                            source2 = row2['source']
                            target2 = row2['target']
                        else:
                            # Find in the full results
                            row2 = result2['results'][(result2['results']['ligand'] == ligand) & 
                                                     (result2['results']['receptor'] == receptor)]
                            if not row2.empty:
                                row2 = row2.sort_values('prob', ascending=False).iloc[0]
                                prob2 = row2['prob']
                                source2 = row2['source']
                                target2 = row2['target']
                            else:
                                prob2 = 0
                                source2 = ""
                                target2 = ""
                        
                        # Add to common data
                        common_data.append({
                            'ligand': ligand,
                            'receptor': receptor,
                            'pair': f"{ligand} → {receptor}",
                            f'prob_{group1_name}': prob1,
                            f'comm_{group1_name}': f"{source1} → {target1}" if source1 and target1 else "-",



                            f'prob_{group2_name}': prob2,
                            f'comm_{group2_name}': f"{source2} → {target2}" if source2 and target2 else "-",
                            'log2_fold_change': np.log2((prob2 + 1e-10) / (prob1 + 1e-10))
                        })
                    
                    # Create DataFrame and display
                    common_df = pd.DataFrame(common_data)
                    common_df = common_df.sort_values('log2_fold_change', key=abs, ascending=False)
                    
                    st.dataframe(common_df[[
                        'pair', 
                        f'comm_{group1_name}', f'prob_{group1_name}',
                        f'comm_{group2_name}', f'prob_{group2_name}',
                        'log2_fold_change'
                    ]])
                    
                    # Plot log2 fold changes
                    try:
                        fig_lr_fc, ax_lr_fc = plt.subplots(figsize=(12, len(common_df) * 0.5 + 2))
                        
                        # Create bar colors based on fold change direction
                        bar_colors = ['red' if x < 0 else 'blue' for x in common_df['log2_fold_change']]
                        
                        bars_lr_fc = ax_lr_fc.barh(
                            common_df['pair'],
                            common_df['log2_fold_change'],
                            color=bar_colors
                        )
                        
                        # Add a vertical line at x=0
                        ax_lr_fc.axvline(x=0, color='black', linestyle='-', linewidth=0.5)
                        
                        ax_lr_fc.set_xlabel('Log2 Fold Change')
                        ax_lr_fc.set_title(f'L-R Pair Fold Change: {group2_name} vs {group1_name}')
                        
                        # Add legends
                        import matplotlib.patches as mpatches
                        red_patch = mpatches.Patch(color='red', label=f'Higher in {group1_name}')
                        blue_patch = mpatches.Patch(color='blue', label=f'Higher in {group2_name}')
                        ax_lr_fc.legend(handles=[red_patch, blue_patch])
                        
                        plt.tight_layout()
                        st.pyplot(fig_lr_fc)
                        
                        # Save the fold change plot
                        lr_fc_pdf_path = f"{cellchat_temp_dir}/lr_fold_change.pdf"
                        fig_lr_fc.savefig(lr_fc_pdf_path, bbox_inches='tight')
                        
                        # Add download button
                        with open(lr_fc_pdf_path, "rb") as pdf_file:
                            pdf_bytes = pdf_file.read()
                            st.download_button(
                                label="Download LR Fold Change Plot",
                                data=pdf_bytes,
                                file_name=f'cellchat_lr_fold_change_{group1_name}_vs_{group2_name}.pdf',
                                mime='application/octet-stream'
                            )
                    except Exception as e:
                        st.error(f"Error generating LR fold change plot: {str(e)}")
                else:
                    st.info("No common L-R pairs found in the top interactions of both groups.")
                
                # Create ZIP of all generated files
                try:
                    all_files = []
                    # Add all PDF files in the temp directory
                    for file in os.listdir(cellchat_temp_dir):
                        if file.endswith('.pdf'):
                            all_files.append(os.path.join(cellchat_temp_dir, file))
                    
                    if all_files:
                        zip_buffer = create_zip_from_files(all_files)
                        if zip_buffer:
                            st.download_button(
                                label="Download all figures as ZIP",
                                data=zip_buffer,
                                file_name=f'cellchat_all_figures_{group1_name}_vs_{group2_name}.zip',
                                mime='application/zip'
                            )
                except Exception as e:
                    st.error(f"Error creating ZIP for all figures: {str(e)}")
            else:
                st.warning("LR interaction data not available for one or both groups")


        with tabs[5]:
            st.markdown("#### Cell Type Signaling Role Comparison")
            
            st.info("""
            In this tab, we compare changes in signaling roles of each cell type between two datasets.
            The left-right position (X-axis) in the scatter plot represents outgoing signal strength, and the up-down position (Y-axis) represents incoming signal strength.
            The size of the circles represents the number of interactions (links).
            Cell types positioned in the upper right have more central roles in the network.
            """)
            
            # Standardize dot size based on link count
            col1, col2 = st.columns(2)
            
            with col1:
                dot_size_min = st.slider("Minimum dot size:", 1, 10, 2)
                dot_size_max = st.slider("Maximum dot size:", dot_size_min+1, 20, 6)
                font_size = st.slider("Font size:", 6, 14, 10)
                do_label = st.checkbox("Show cell type labels", value=True)
            
            with col2:
                role_x = st.slider("Fig width:", min_value=1, max_value=20, value=8)
                role_y = st.slider("Fig height:", min_value=1, max_value=20, value=5)
            
            if st.button("Generate Cell Signaling Role Plots"):
                with st.spinner("Generating plots..."):
                    try:
                        # Get network link counts for each dataset
                        from pages.cellchat_vis import netAnalysis_signalingRole_scatter
                        
                        # Get link counts from both datasets and apply unified scale
                        count_matrix1 = result1['network']['count_matrix']
                        count_matrix2 = result2['network']['count_matrix']
                        
                        # Calculate link count: row sum + column sum - diagonal
                        if hasattr(count_matrix1, 'values'):
                            num_link1 = np.sum(count_matrix1.values, axis=0) + np.sum(count_matrix1.values, axis=1) - np.diag(count_matrix1.values)
                            num_link2 = np.sum(count_matrix2.values, axis=0) + np.sum(count_matrix2.values, axis=1) - np.diag(count_matrix2.values)
                        else:
                            num_link1 = np.sum(count_matrix1, axis=0) + np.sum(count_matrix1, axis=1) - np.diag(count_matrix1)
                            num_link2 = np.sum(count_matrix2, axis=0) + np.sum(count_matrix2, axis=1) - np.diag(count_matrix2)
                        
                        num_link = np.concatenate([num_link1, num_link2])
                    #    st.write(num_link)
                        weight_MinMax = [np.min(num_link), np.max(num_link)]
                    #    st.write(weight_MinMax)

                        # Plot for group 1
                        fig1 = netAnalysis_signalingRole_scatter(
                            result1,
                            title=f"Signaling Roles - {group1_name}",
                            weight_MinMax=weight_MinMax,
                            dot_size=(dot_size_min, dot_size_max),
                            font_size=font_size,
                            do_label=do_label,
                            color_use=st.session_state.get('cell_color_map', None),
                            width=role_x,
                            height=role_y,
                            use_count=True,
                            sources_use=sources_use if filter_cells else None,
                            targets_use=targets_use if filter_cells else None
                        )

                                     
                        # Plot for group 2
                        fig2 = netAnalysis_signalingRole_scatter(
                            result2,
                            title=f"Signaling Roles - {group2_name}",
                            weight_MinMax=weight_MinMax,
                            dot_size=(dot_size_min, dot_size_max),
                            font_size=font_size,
                            do_label=do_label,
                            color_use=st.session_state.get('cell_color_map', None),
                            width=role_x,
                            height=role_y,
                            use_count=True,
                            sources_use=sources_use if filter_cells else None,
                            targets_use=targets_use if filter_cells else None
                        )
                        # Align min-max for x, y
                        ax1 = fig1.axes[0]  # Get Axes object from fig1
                        ax2 = fig2.axes[0]  # Get Axes object from fig2

                        xlim1 = ax1.get_xlim()
                        xlim2 = ax2.get_xlim()
                        ylim1 = ax1.get_ylim()
                        ylim2 = ax2.get_ylim()
                        global_xmin = min(xlim1[0], xlim2[0])
                        global_xmax = max(xlim1[1], xlim2[1])
                        global_ymin = min(ylim1[0], ylim2[0])
                        global_ymax = max(ylim1[1], ylim2[1])

                        ax1.set_xlim((global_xmin, global_xmax))
                        ax2.set_xlim((global_xmin, global_xmax))
                        ax1.set_ylim((global_ymin, global_ymax))
                        ax2.set_ylim((global_ymin, global_ymax))

                        col1, col2 = st.columns(2)
                        
                        with col1:
                            st.pyplot(fig1)
                            
                            # Save PDF and download
                            pdf_path1 = f"{cellchat_temp_dir}/signaling_roles_{group1_name}.pdf"
                            fig1.savefig(pdf_path1, bbox_inches='tight')
                            
                            with open(pdf_path1, "rb") as pdf_file:
                                st.download_button(
                                    label=f"Download {group1_name} PDF",
                                    data=pdf_file.read(),
                                    file_name=f'cellchat_signaling_roles_{group1_name}.pdf',
                                    mime='application/octet-stream'
                                )
                        with col2:
                            st.pyplot(fig2)                       
                            # Save PDF and download
                            pdf_path2 = f"{cellchat_temp_dir}/signaling_roles_{group2_name}.pdf"
                            fig2.savefig(pdf_path2, bbox_inches='tight')
                            
                            with open(pdf_path2, "rb") as pdf_file:
                                st.download_button(
                                    label=f"Download {group2_name} PDF",
                                    data=pdf_file.read(),
                                    file_name=f'cellchat_signaling_roles_{group2_name}.pdf',
                                    mime='application/octet-stream'
                                )
                        
                        # Download all plots as ZIP
                        zip_files = [pdf_path1, pdf_path2]
                        if len(zip_files) > 0:
                            try:
                                zip_buffer = create_zip_from_files(zip_files)
                                if zip_buffer:
                                    st.download_button(
                                        label="Download all as ZIP",
                                        data=zip_buffer,
                                        file_name=f'cellchat_signaling_roles_{group1_name}_vs_{group2_name}.zip',
                                        mime='application/zip'
                                    )
                            except Exception as e:
                                st.error(f"Error creating ZIP file: {str(e)}")
                        
                    except Exception as e:
                        st.error(f"Error generating signaling role plots: {str(e)}")
                        st.exception(e)

        # New tab - Cell-specific signaling changes
        # Update only the 'Signaling Changes' tab section


        with tabs[6]:
            st.markdown("#### Cell-Specific Signaling Changes Analysis")

            st.info("""
            This tab analyzes signaling pathway changes between two datasets, focusing on a specific cell type.

            - X-axis: Difference in outgoing signal strength (Group 2 - Group 1)
            - Y-axis: Difference in incoming signal strength (Group 2 - Group 1)
            - Shape:
              - Square (□): Incoming-specific changes
              - Triangle (△): Outgoing-specific changes
              - Diamond (◇): Both incoming and outgoing specific changes

            Pathways in the upper-right region represent increased signals in the second dataset.
            """)

            # Extract common cell types
            cell_types1 = result1['net']['dimnames'][0]
            cell_types2 = result2['net']['dimnames'][0]
            common_cells = sorted(set(cell_types1).intersection(set(cell_types2)))
            
            if not common_cells:
                st.warning("No common cell types found between the two datasets")
            else:
                # Cell type selection
                idents_use = st.selectbox(
                    "Select cell type to analyze:",
                    options=common_cells,
                    index=0
                )
                
                # Pathway selection/exclusion options
                all_pathways1 = set(result1['netP']['pathways']) if 'netP' in result1 and 'pathways' in result1['netP'] else set()
                all_pathways2 = set(result2['netP']['pathways']) if 'netP' in result2 and 'pathways' in result2['netP'] else set()
                all_pathways = sorted(all_pathways1.union(all_pathways2))
                
                if not all_pathways:
                    st.warning("Pathway information not found")
                else:
                    # Pathway selection method
                    pathway_selection = st.radio(
                        "Pathway selection method:",
                        ["Show all", "Select specific pathways", "Exclude specific pathways"],
                        horizontal=True
                    )
                    
                    signaling_include = None
                    signaling_exclude = None
                    
                    if pathway_selection == "Select specific pathways":
                        signaling_include = st.multiselect(
                            "Select pathways to include:",
                            options=all_pathways,
                            default=all_pathways[:min(5, len(all_pathways))]
                        )
                    elif pathway_selection == "Exclude specific pathways":
                        default_exclude = [p for p in ["MIF"] if p in all_pathways]
                        signaling_exclude = st.multiselect(
                            "Select pathways to exclude:",
                            options=all_pathways,
                            default=default_exclude
                        )
                    
                    # Plot settings
                    col1, col2, col3 = st.columns(3)
                    
                    with col1:
                        dot_size = st.slider("Dot size:", 2, 10, 3)
                        font_size = st.slider("Font size:", 8, 16, 12)
                        scatter_x = st.slider("Scatter width:", min_value=1, max_value=20, value=6)
                        scatter_y = st.slider("Scatter height:", min_value=1, max_value=20, value=6)
                                                    
                    with col2:
                        # Add color settings
                        color_mode = st.radio(
                            "Color mode:",
                            ["Single color", "Color by pattern"],
                            horizontal=True,
                            help="Single color: Use same color for all dots. Color by pattern: Use different colors for each pattern (Incoming/Outgoing/Both/Other)",
                            index=1
                        )
                        
                        if color_mode == "Single color":
                            use_color_by_pattern = False
                            color_use = st.selectbox(
                                "Dot color:",
                                ["teal", "blue", "green", "red", "purple", "orange"],
                                index=0
                            )
                            color_palette = None
                        else:
                            use_color_by_pattern = True
                            color_use = None
                            color_palette = st.selectbox(
                                "Color palette:",
                                ["Set1", "Set2", "Set3", "tab10", "tab20", "Dark2", "Pastel1", "Pastel2", "Paired", "Accent"],
                                index=9,
                                help="Color palette used for coloring by pattern"
                            )
                                                # Label settings
                        max_label = st.slider("Number of labels to display:", 5, 50, 20)

                        use_arrows = st.checkbox("Show labels with arrows", value=True, 
                                                help="When ON, labels are displayed using arrows. When OFF, labels are displayed directly next to points.")

                        reverse_groups = st.checkbox("Reverse group order", value=False,
                                                 help="When checked, the difference between group 1 and group 2 is reversed during analysis.")


                    
                    # Run button
                    if st.button("Generate Signaling Changes Plot"):
                        with st.spinner("Analyzing signaling changes..."):
                            try:
                                # Create merged CellChat object
                                merged_cellchat = {
                                    'group1': {
                                        'name': group1_name,
                                        'result': result1
                                    },
                                    'group2': {
                                        'name': group2_name,
                                        'result': result2
                                    }
                                }
                                
                                # If reverse switch is on, swap group order
                                if reverse_groups:
                                    merged_cellchat = {
                                        'group1': {
                                            'name': group2_name,
                                            'result': result2
                                        },
                                        'group2': {
                                            'name': group1_name,
                                            'result': result1
                                        }
                                    }
                                
                                # Run analysis
                                fig = netAnalysis_signalingChanges_scatter(
                                    merged_cellchat,
                                    idents_use=idents_use,
                                    signaling_exclude=signaling_exclude,
                                    signaling_include=signaling_include,
                                    color_use=color_use,
                                    use_color_by_pattern=use_color_by_pattern,  # New parameter
                                    color_palette=color_palette,  # New parameter
                                    font_size=font_size,
                                    dot_size=dot_size,
                                    max_label=max_label,
                                    use_arrows=use_arrows,
                                    figsize=(scatter_x, scatter_y)
                                )
                                
                                # Display plot
                                st.pyplot(fig)
                                
                                # Save PDF and download
                                pdf_path = f"{cellchat_temp_dir}/signaling_changes_{idents_use}.pdf"
                                fig.savefig(pdf_path, bbox_inches='tight')
                                
                                with open(pdf_path, "rb") as pdf_file:
                                    st.download_button(
                                        label=f"Download PDF",
                                        data=pdf_file.read(),
                                        file_name=f'cellchat_signaling_changes_{idents_use}_{group1_name}_vs_{group2_name}.pdf',
                                        mime='application/octet-stream'
                                    )
                            
                            except Exception as e:
                                st.error(f"Error analyzing signaling changes: {str(e)}")
                                st.exception(e)

        # Add download section at the bottom (always visible when results are available)
        st.markdown("---")
        st.markdown("#### 💾 Save Analysis Results")
        st.info("Download your analysis results for future use or sharing.")
        
        # Prepare data for download
        serialized_results, filename = save_cellchat_results(
            result1, result2, group1_name, group2_name,
            uploaded_filename=None,  # We don't have the original filename in visualization mode
            selected_types=None
        )
        
        if serialized_results:
            col1, col2 = st.columns(2)
            with col1:
                st.download_button(
                    label="📥 Download Combined Results",
                    data=serialized_results,
                    file_name=filename,
                    mime="application/octet-stream",
                    help="Download both group analyses as a single file for later import",
                    use_container_width=True
                )
            with col2:
                st.info(f"**File name:** {filename}\n\n**Contains:** Both group analyses and metadata")
        
    else:
        st.info("To visualize a CellChat results, please upload the pkl file.")
        
        # Add example download option
        st.subheader("No saved results?")
        st.write("If you don't have saved results yet, you can:")
        st.write("1. Use the 'Upload H5AD file for new analysis' option to generate new results")
        st.write("2. After analysis, use the 'Download Combined Results' button to save your results for future use")

