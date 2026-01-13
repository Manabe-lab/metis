import streamlit as st
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, numpy2ri
from rpy2.robjects.packages import importr
import numpy as np
import pandas as pd
import pickle
import os
import tempfile
import traceback

# App title and header
st.title("CellChat QS Converter")
st.markdown("A tool to convert CellChat QS files to Python pickle files")

# Initialize R environment (runs only once at session start)
@st.cache_resource
def initialize_r_environment():
    pandas2ri.activate()
    numpy2ri.activate()
    try:
        base = importr('base')
        qs = importr('qs')
        # Also import CellChat
        ro.r('library(CellChat)')
        return True
    except Exception as e:
        st.error(f"Failed to initialize R environment: {str(e)}")
        st.error("Please verify that the required R packages (CellChat, qs) are installed.")
        return False

# Initialize R environment
r_initialized = initialize_r_environment()

# Function to process QS files
def process_qs_file(uploaded_file):
    # Save to temporary file
    with tempfile.NamedTemporaryFile(suffix='.qs', delete=False) as tmp_file:
        tmp_file.write(uploaded_file.getvalue())
        tmp_path = tmp_file.name
    
    try:
        # Escape special characters
        tmp_path_escaped = tmp_path.replace('\\', '\\\\')

        # Simplest approach: first set a variable in R, then use it
        ro.r(f'temp_file_path <- "{tmp_path_escaped}"')
        
        r_code = """
        library(CellChat)
        library(qs)

        # Load QS file
        tryCatch({
            cellchat_obj <- qs::qread(temp_file_path)

            # Extract required data
            extracted_data <- list()

            # 1. Metadata
            if (!is.null(cellchat_obj@meta)) {
                extracted_data$meta <- as.data.frame(cellchat_obj@meta)
            }

            # 2. Network information
            if (!is.null(cellchat_obj@net$prob)) {
                # Probability matrix
                extracted_data$net_prob <- cellchat_obj@net$prob

                # Save dimension names individually
                if (!is.null(dimnames(cellchat_obj@net$prob)[[1]])) {
                    extracted_data$dim_rows <- dimnames(cellchat_obj@net$prob)[[1]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[2]])) {
                    extracted_data$dim_cols <- dimnames(cellchat_obj@net$prob)[[2]]
                }
                if (!is.null(dimnames(cellchat_obj@net$prob)[[3]])) {
                    extracted_data$dim_paths <- dimnames(cellchat_obj@net$prob)[[3]]
                }
            }

            if (!is.null(cellchat_obj@net$pval)) {
                extracted_data$net_pval <- cellchat_obj@net$pval
            }

            if (!is.null(cellchat_obj@net$weight)) {
                extracted_data$net_weight <- as.matrix(cellchat_obj@net$weight)
            }

            if (!is.null(cellchat_obj@net$count)) {
                extracted_data$net_count <- as.matrix(cellchat_obj@net$count)
            }

            # 3. Pathway information
            if (!is.null(cellchat_obj@netP$pathways)) {
                extracted_data$pathways <- cellchat_obj@netP$pathways
            }

            if (!is.null(cellchat_obj@netP$prob)) {
                extracted_data$netP_prob <- cellchat_obj@netP$prob
            }

            # 4. LR information
            if (!is.null(cellchat_obj@LR$LRsig)) {
                extracted_data$lr_sig <- as.data.frame(cellchat_obj@LR$LRsig)
            }

            # Return value
            extracted_data
        }, error = function(e) {
            list(error = paste("Error reading QS file:", e$message))
        })
        """
        
        # Execute R code
        result_r = ro.r(r_code)

        # Error check
        if isinstance(result_r, ro.vectors.ListVector) and 'error' in result_r.names:
            error_msg = result_r.rx2('error')
            st.error(f"Failed to process QS file: {error_msg}")
            return None

        # Convert results
        cellchat_data = {}
        for key in result_r.names:
            try:
                # Change conversion method based on key type
                if key in ["net_prob", "net_pval", "netP_prob"]:
                    # Convert R array to NumPy array
                    cellchat_data[key] = np.array(result_r.rx2(key))
                else:
                    # Convert other objects using pandas2ri
                    r_value = result_r.rx2(key)
                    try:
                        cellchat_data[key] = pandas2ri.rpy2py(r_value)
                    except Exception as e2:
                        # If conversion fails, get data using basic methods
                        if isinstance(r_value, ro.vectors.StrVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.FloatVector):
                            cellchat_data[key] = list(r_value)
                        elif isinstance(r_value, ro.vectors.IntVector):
                            cellchat_data[key] = list(r_value)
                        else:
                            cellchat_data[key] = str(r_value)
                st.success(f"Successfully extracted {key}")
            except Exception as e:
                st.warning(f"Failed to convert {key}: {e}")
                st.warning(traceback.format_exc())

        # Organize dimension names
        if all(key in cellchat_data for key in ['dim_rows', 'dim_cols', 'dim_paths']):
            cellchat_data['net_prob_dims'] = {
                'rows': cellchat_data['dim_rows'],
                'cols': cellchat_data['dim_cols'],
                'pathways': cellchat_data['dim_paths']
            }
            # Remove individual dimension keys (optional)
            for key in ['dim_rows', 'dim_cols', 'dim_paths']:
                if key in cellchat_data:
                    del cellchat_data[key]

        # Format results
        result = format_cellchat_for_visualization(cellchat_data)

        # Delete temporary file
        os.unlink(tmp_path)

        return result

    except Exception as e:
        st.error(f"An error occurred during processing: {e}")
        st.error(traceback.format_exc())
        # Try to delete temporary file
        try:
            os.unlink(tmp_path)
        except:
            pass
        return None

def format_cellchat_for_visualization(cellchat_data):
    """Format extracted data for visualization in cellchat_vis.py"""
    # Structure data to match cellchat.py output format
    result = {}

    # 1. Network information
    # Check if 3D arrays for prob and pval exist
    has_3d_prob = 'net_prob' in cellchat_data and isinstance(cellchat_data['net_prob'], np.ndarray) and cellchat_data['net_prob'].ndim == 3
    has_3d_pval = 'net_pval' in cellchat_data and isinstance(cellchat_data['net_pval'], np.ndarray) and cellchat_data['net_pval'].ndim == 3
    
    result['net'] = {}

    # Set only if prob and pval are 3-dimensional
    if has_3d_prob:
        result['net']["prob"] = cellchat_data['net_prob']
    else:
        # Create empty 3D array if no data
        if 'net_prob_dims' in cellchat_data and all(key in cellchat_data['net_prob_dims'] for key in ['rows', 'cols']):
            rows = len(cellchat_data['net_prob_dims']['rows'])
            cols = len(cellchat_data['net_prob_dims']['cols'])
            paths = len(cellchat_data['net_prob_dims'].get('pathways', []))
            # Ensure at least one pathway
            paths = max(1, paths)
            result['net']["prob"] = np.zeros((rows, cols, paths))
        else:
            # Create minimum size array if no dimension information
            result['net']["prob"] = np.zeros((1, 1, 1))
    
    if has_3d_pval:
        result['net']["pval"] = cellchat_data['net_pval']
    else:
        # Create array with same size as prob
        result['net']["pval"] = np.ones_like(result['net']["prob"])

    # Add dimension name information
    if 'net_prob_dims' in cellchat_data:
        dims = cellchat_data['net_prob_dims']
        cell_types_rows = dims.get('rows', [])
        cell_types_cols = dims.get('cols', [])
        pathways = dims.get('pathways', [])

        # Convert to string list
        cell_types_rows = [str(ct) for ct in cell_types_rows]
        cell_types_cols = [str(ct) for ct in cell_types_cols]
        pathways = [str(p) for p in pathways]

        # Ensure minimum size
        if not cell_types_rows:
            cell_types_rows = [f"Cell_{i}" for i in range(result['net']["prob"].shape[0])]
        if not cell_types_cols:
            cell_types_cols = [f"Cell_{i}" for i in range(result['net']["prob"].shape[1])]
        if not pathways:
            pathways = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]
        
        result['net']['dimnames'] = [
            cell_types_rows,
            cell_types_cols,
            pathways
        ]
    else:
        # Generate dimension names based on array size if not available
        shape = result['net']["prob"].shape
        result['net']['dimnames'] = [
            [f"Cell_{i}" for i in range(shape[0])],
            [f"Cell_{i}" for i in range(shape[1])],
            [f"Pathway_{i}" for i in range(shape[2])]
        ]

    # Create aggregated weight matrix
    cell_types = result['net']['dimnames'][0]
    num_cells = len(cell_types)

    # Sum across all pathways (aggregated network)
    prob_sum = np.sum(result['net']["prob"], axis=2)

    # Convert weight matrix and count matrix to DataFrame
    result['net']["weight"] = prob_sum
    result['net']["count"] = (prob_sum > 0).astype(int)

    # Calculate centrality of aggregated network (using NetworkX)
    import networkx as nx

    # Initialize
    net_centr = {
        'outdeg': np.zeros(num_cells),
        'indeg': np.zeros(num_cells),
        'outdeg_unweighted': np.zeros(num_cells),
        'indeg_unweighted': np.zeros(num_cells),
        'betweenness': np.zeros(num_cells),
        'page_rank': np.zeros(num_cells),
        'hub': np.zeros(num_cells),
        'authority': np.zeros(num_cells),
        'eigen': np.zeros(num_cells),
        'flowbet': np.zeros(num_cells),
        'info': np.zeros(num_cells)
    }

    # Calculate only if network is valid (has edges)
    if np.any(prob_sum > 0):
        # Create graph
        G = nx.DiGraph()
        for i in range(num_cells):
            G.add_node(i, name=cell_types[i])

        # Add edges
        for i in range(num_cells):
            for j in range(num_cells):
                if prob_sum[i, j] > 0:
                    G.add_edge(i, j, weight=float(prob_sum[i, j]))

        # In/out degree (weighted/unweighted)
        for i in range(num_cells):
            # Weighted degree
            net_centr['outdeg'][i] = sum(G[i][j]['weight'] for j in G.successors(i)) if i in G else 0
            net_centr['indeg'][i] = sum(G[j][i]['weight'] for j in G.predecessors(i)) if i in G else 0

            # Unweighted degree
            net_centr['outdeg_unweighted'][i] = G.out_degree(i) if i in G else 0
            net_centr['indeg_unweighted'][i] = G.in_degree(i) if i in G else 0

        # Inverse weight graph (for betweenness centrality)
        G_inv = G.copy()
        for u, v, d in G_inv.edges(data=True):
            if d['weight'] > 0:
                d['weight'] = 1.0 / d['weight']

        try:
            # Betweenness centrality
            betweenness = nx.betweenness_centrality(G_inv, weight='weight')
            for i in range(num_cells):
                net_centr['betweenness'][i] = betweenness.get(i, 0)
                # Use betweenness as substitute for flowbet
                net_centr['flowbet'][i] = betweenness.get(i, 0)
        except:
            pass

        try:
            # PageRank
            page_rank = nx.pagerank(G, weight='weight')
            for i in range(num_cells):
                net_centr['page_rank'][i] = page_rank.get(i, 0)
        except:
            pass

        try:
            # Hub and Authority
            hub, authority = nx.hits(G, max_iter=1000)
            for i in range(num_cells):
                net_centr['hub'][i] = hub.get(i, 0)
                net_centr['authority'][i] = authority.get(i, 0)
        except:
            pass

        try:
            # Eigenvector centrality
            eigen = nx.eigenvector_centrality(G, weight='weight', max_iter=1000)
            for i in range(num_cells):
                net_centr['eigen'][i] = eigen.get(i, 0)
        except:
            pass

        # Information Centrality - custom calculation (approximation)
        try:
            # Simple approximation of info centrality (combination of normalized degree and betweenness)
            in_degree = np.array([net_centr['indeg'][i] for i in range(num_cells)])
            out_degree = np.array([net_centr['outdeg'][i] for i in range(num_cells)])
            between = np.array([net_centr['betweenness'][i] for i in range(num_cells)])

            # Normalize
            norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
            norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
            norm_between = between / (np.max(between) if np.max(between) > 0 else 1)

            # Synthesize information centrality
            info_score = (norm_in + norm_out + norm_between) / 3
            for i in range(num_cells):
                net_centr['info'][i] = info_score[i]
        except:
            pass

    # Add centrality metrics to net
    result['net']['centr'] = net_centr

    # 2. Pathway information
    result['netP'] = {}

    # Add pathway information
    if 'pathways' in cellchat_data and cellchat_data['pathways'] is not None:
        result['netP']['pathways'] = cellchat_data['pathways']
    elif 'net_prob_dims' in cellchat_data and 'pathways' in cellchat_data['net_prob_dims']:
        result['netP']['pathways'] = cellchat_data['net_prob_dims']['pathways']
    else:
        # If still not found, get from net['dimnames']
        if 'dimnames' in result['net'] and len(result['net']['dimnames']) > 2:
            result['netP']['pathways'] = result['net']['dimnames'][2]
        else:
            result['netP']['pathways'] = [f"Pathway_{i}" for i in range(result['net']["prob"].shape[2])]

    # Add network probability information
    if 'netP_prob' in cellchat_data and cellchat_data['netP_prob'] is not None:
        result['netP']['prob'] = cellchat_data['netP_prob']
    else:
        # If netP_prob is not available, aggregate by pathway
        prob = result['net']["prob"]
        num_pathways = len(result['netP']['pathways'])

        # Calculate sum for each pathway
        netP_prob = np.zeros((num_cells, num_cells, num_pathways))
        for k in range(min(prob.shape[2], num_pathways)):
            netP_prob[:, :, k] = prob[:, :, k]

        result['netP']['prob'] = netP_prob

    # Calculate centrality information for netP
    netP_centr = {}

    # Calculate centrality for each pathway
    for k in range(min(result['netP']['prob'].shape[2], len(result['netP']['pathways']))):
        # Initial centrality dictionary for this pathway
        netP_centr[k] = {
            'outdeg': np.zeros(num_cells),
            'indeg': np.zeros(num_cells),
            'outdeg_unweighted': np.zeros(num_cells),
            'indeg_unweighted': np.zeros(num_cells),
            'betweenness': np.zeros(num_cells),
            'page_rank': np.zeros(num_cells),
            'hub': np.zeros(num_cells),
            'authority': np.zeros(num_cells),
            'eigen': np.zeros(num_cells),
            'flowbet': np.zeros(num_cells),
            'info': np.zeros(num_cells)
        }

        # Pathway matrix
        pathway_weight = result['netP']['prob'][:, :, k]

        # Calculate only if network is valid (has edges)
        if np.any(pathway_weight > 0):
            # Create graph
            G_pathway = nx.DiGraph()
            for i in range(num_cells):
                G_pathway.add_node(i, name=cell_types[i])

            # Add edges
            for i in range(num_cells):
                for j in range(num_cells):
                    if pathway_weight[i, j] > 0:
                        G_pathway.add_edge(i, j, weight=float(pathway_weight[i, j]))

            # In/out degree (weighted/unweighted)
            for i in range(num_cells):
                # Weighted degree
                netP_centr[k]['outdeg'][i] = sum(G_pathway[i][j]['weight'] for j in G_pathway.successors(i)) if i in G_pathway else 0
                netP_centr[k]['indeg'][i] = sum(G_pathway[j][i]['weight'] for j in G_pathway.predecessors(i)) if i in G_pathway else 0

                # Unweighted degree
                netP_centr[k]['outdeg_unweighted'][i] = G_pathway.out_degree(i) if i in G_pathway else 0
                netP_centr[k]['indeg_unweighted'][i] = G_pathway.in_degree(i) if i in G_pathway else 0

            # Inverse weight graph (for betweenness centrality)
            G_inv = G_pathway.copy()
            for u, v, d in G_inv.edges(data=True):
                if d['weight'] > 0:
                    d['weight'] = 1.0 / d['weight']

            try:
                # Betweenness centrality
                betweenness = nx.betweenness_centrality(G_inv, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['betweenness'][i] = betweenness.get(i, 0)
                    # Use betweenness as substitute for flowbet
                    netP_centr[k]['flowbet'][i] = betweenness.get(i, 0)
            except:
                pass

            try:
                # PageRank
                page_rank = nx.pagerank(G_pathway, weight='weight')
                for i in range(num_cells):
                    netP_centr[k]['page_rank'][i] = page_rank.get(i, 0)
            except:
                pass

            try:
                # Hub and Authority
                hub, authority = nx.hits(G_pathway, max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['hub'][i] = hub.get(i, 0)
                    netP_centr[k]['authority'][i] = authority.get(i, 0)
            except:
                pass

            try:
                # Eigenvector centrality
                eigen = nx.eigenvector_centrality(G_pathway, weight='weight', max_iter=1000)
                for i in range(num_cells):
                    netP_centr[k]['eigen'][i] = eigen.get(i, 0)
            except:
                pass

            # Information Centrality - custom calculation
            try:
                # Simple approximation of info centrality (combination of normalized degree and betweenness)
                in_degree = np.array([netP_centr[k]['indeg'][i] for i in range(num_cells)])
                out_degree = np.array([netP_centr[k]['outdeg'][i] for i in range(num_cells)])
                between = np.array([netP_centr[k]['betweenness'][i] for i in range(num_cells)])

                # Normalize
                norm_in = in_degree / (np.max(in_degree) if np.max(in_degree) > 0 else 1)
                norm_out = out_degree / (np.max(out_degree) if np.max(out_degree) > 0 else 1)
                norm_between = between / (np.max(between) if np.max(between) > 0 else 1)

                # Synthesize information centrality
                info_score = (norm_in + norm_out + norm_between) / 3
                for i in range(num_cells):
                    netP_centr[k]['info'][i] = info_score[i]
            except:
                pass

    # Add centrality information
    result['netP']['centr'] = netP_centr

    # 3. Create results DataFrame
    # Get cell type names and pathway names
    cell_types_rows = result['net']['dimnames'][0]
    cell_types_cols = result['net']['dimnames'][1]
    pathways = result['net']['dimnames'][2]

    # Extract interaction data from prob, pval arrays
    prob = result['net']["prob"]
    pval = result['net']["pval"]
    
    results_data = []
    for i in range(prob.shape[0]):
        for j in range(prob.shape[1]):
            for k in range(min(prob.shape[2], len(pathways))):
                if prob[i, j, k] > 0:
                    # Interaction name
                    interaction_name = pathways[k] if k < len(pathways) else f"interaction_{k}"

                    # Ligand-receptor information (if available)
                    ligand = ""
                    receptor = ""

                    # Get LR information
                    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
                        lr_df = cellchat_data['lr_sig']
                        if interaction_name in lr_df.index:
                            if 'ligand' in lr_df.columns:
                                ligand = lr_df.loc[interaction_name, 'ligand']
                            if 'receptor' in lr_df.columns:
                                receptor = lr_df.loc[interaction_name, 'receptor']
                    
                    results_data.append({
                        'source': cell_types_rows[i] if i < len(cell_types_rows) else f"Cell_{i}",
                        'target': cell_types_cols[j] if j < len(cell_types_cols) else f"Cell_{j}",
                        'interaction_name': interaction_name,
                        'ligand': str(ligand),
                        'receptor': str(receptor),
                        'prob': float(prob[i, j, k]),
                        'pval': float(pval[i, j, k])
                    })
    
    result['results'] = pd.DataFrame(results_data)

    # 4. Create network summary
    # Aggregate by cell type (sum of prob)
    prob = result['net']["prob"]
    strength_matrix = np.sum(prob, axis=2)

    # Convert to DataFrame
    strength_df = pd.DataFrame(strength_matrix,
                              index=cell_types,
                              columns=cell_types)

    # Count matrix (number of non-zero interactions)
    count_matrix = np.sum(prob > 0, axis=2)
    count_df = pd.DataFrame(count_matrix,
                           index=cell_types,
                           columns=cell_types)

    # Total outgoing and incoming
    outgoing = pd.Series(np.sum(strength_matrix, axis=1), index=cell_types)
    incoming = pd.Series(np.sum(strength_matrix, axis=0), index=cell_types)

    # Calculate LR contribution
    lr_contribution = np.zeros(prob.shape[2])
    for k in range(prob.shape[2]):
        lr_contribution[k] = np.sum(prob[:, :, k])

    # Create network centrality dataframe
    network_centrality_dict = {
        'cell_type': cell_types,
    }

    # Add each centrality metric to dictionary
    for metric in ['outdeg', 'indeg', 'outdeg_unweighted', 'indeg_unweighted',
                  'betweenness', 'page_rank', 'hub', 'authority', 'eigen', 'flowbet', 'info']:
        network_centrality_dict[metric] = net_centr[metric]

    network_centrality = pd.DataFrame(network_centrality_dict)

    # Rename core centrality metrics to more meaningful names (for signaling role analysis)
    role_data = {
        'cell_type': cell_types,
        'sender': net_centr['outdeg'],        # Role as sender
        'receiver': net_centr['indeg'],       # Role as receiver
        'mediator': net_centr['flowbet'],     # Role as mediator
        'influencer': net_centr['info']       # Role as influencer
    }

    # Create role dataframe
    signaling_role = pd.DataFrame(role_data)

    result['network'] = {
        'strength_matrix': strength_df,
        'count_matrix': count_df,
        'outgoing': outgoing,
        'incoming': incoming,
        'lr_contribution': lr_contribution,
        'network_centrality': network_centrality,
        'signaling_role': signaling_role      # For signaling role analysis
    }

    # 5. Metadata information
    if 'meta' in cellchat_data:
        result['meta'] = cellchat_data['meta']

    # 6. GroupBy information (default setting)
    result['groupby'] = 'cell.ident'  # Default value

    # 7. LR information
    if 'lr_sig' in cellchat_data:
        result['lr_sig'] = cellchat_data['lr_sig']

    # 8. adata information - create pseudo anndata structure
    from types import SimpleNamespace

    # Create adata object
    adata = SimpleNamespace()

    # Add metadata
    if 'meta' in cellchat_data and isinstance(cellchat_data['meta'], pd.DataFrame):
        adata.obs = cellchat_data['meta']
    else:
        adata.obs = pd.DataFrame(index=cell_types)

    # Extract gene name list
    gene_names = []

    # Extract ligand and receptor information from LRsig
    if 'lr_sig' in cellchat_data and isinstance(cellchat_data['lr_sig'], pd.DataFrame):
        lr_df = cellchat_data['lr_sig']
        if 'ligand' in lr_df.columns:
            ligands = lr_df['ligand'].dropna().astype(str).unique()
            gene_names.extend(ligands)
        if 'receptor' in lr_df.columns:
            receptors = lr_df['receptor'].dropna().astype(str).unique()
            gene_names.extend(receptors)

    # Fallback if no gene names
    if not gene_names:
        gene_names = ["Gene1", "Gene2", "Gene3"]

    # Remove duplicates
    gene_names = sorted(list(set(gene_names)))

    # Mimic AnnData var attribute and var_names
    adata.var_names = pd.Index(gene_names)
    adata.var = pd.DataFrame(index=adata.var_names)

    # Add AnnData shape attribute
    adata.shape = (len(adata.obs), len(adata.var_names))

    # Add X attribute
    adata.X = np.zeros(adata.shape)
    
    result['adata'] = adata
    
    return result

# Main processing
if r_initialized:
    # File uploader
    uploaded_file = st.file_uploader("Upload CellChat QS file", type=["qs"])

    if uploaded_file is not None:
        st.info(f"Processing file '{uploaded_file.name}'...")

        with st.spinner("Processing QS file..."):
            result = process_qs_file(uploaded_file)

        if result is not None:
            st.success("QS file conversion completed!")

            # Display summary of results
            st.subheader("Conversion Results Summary")

            # Display cell type information
            if 'net' in result and 'dimnames' in result['net'] and len(result['net']['dimnames']) > 0:
                cell_types = result['net']['dimnames'][0]
                cell_count = len(cell_types)
                if cell_count > 0:
                    st.write(f"Number of detected cell types: {cell_count}")
                    if cell_count <= 20:
                        st.write(f"Cell types: {', '.join(cell_types)}")
                    else:
                        with st.expander("Cell type list"):
                            st.write(", ".join(cell_types))

            # Display pathway information
            if 'netP' in result and 'pathways' in result['netP'] and result['netP']['pathways'] is not None:
                pathways = result['netP']['pathways']
                pathway_count = len(pathways)
                if pathway_count > 0:
                    st.write(f"Number of signaling pathways: {pathway_count}")
                    if pathway_count <= 10:
                        st.write(f"Signaling pathways: {', '.join(pathways)}")
                    else:
                        with st.expander("Signaling pathway list"):
                            st.write(", ".join(pathways))

            # Display interaction information
            if 'results' in result and not result['results'].empty:
                st.write(f"Number of extracted interactions: {len(result['results'])}")
                with st.expander("Sample interaction data (top 5)"):
                    st.dataframe(result['results'].head())

            # Save results to Pickle file
            output_filename = os.path.splitext(uploaded_file.name)[0] + "_python.pkl"
            pickle_data = pickle.dumps(result)

            # Download button
            st.download_button(
                label="Download conversion results",
                data=pickle_data,
                file_name=output_filename,
                mime="application/octet-stream"
            )

            st.info("""
            The downloaded pickle file can be used for CellChat visualization in cellchat.py.
            Upload it to the Streamlit app or load it in your code using pickle.load().
            """)
else:
    st.error("Cannot run the app because R environment initialization failed. Please verify that the required packages are installed.")