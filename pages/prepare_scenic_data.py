import streamlit as st
import pandas as pd
import numpy as np
import scanpy as sc
from scipy import stats
# from pyscenic.rss import regulon_specificity_scores  # Has numpy deprecation issue
from scipy.spatial.distance import jensenshannon
import warnings
import zipfile
import io
import plotly.express as px
warnings.filterwarnings('ignore')

# Cache functions for data loading
@st.cache_data
def load_auc_data(file_contents, filename):
    """Load AUC data with caching"""
    return pd.read_csv(file_contents, sep='\t', index_col=0)

@st.cache_data
def load_metadata_tsv(file_contents, filename):
    """Load metadata from TSV with caching"""
    return pd.read_csv(file_contents, sep='\t', index_col=0)

@st.cache_data
def load_metadata_h5ad(file_contents, filename):
    """Load metadata from h5ad with caching"""
    adata = sc.read_h5ad(file_contents)
    return adata.obs

@st.cache_data
def filter_regulons_by_cells(auc_df, min_cells=10):
    """Filter regulons based on number of cells where they are active (AUC > 0)"""
    # Count cells with AUC > 0 for each regulon
    active_cells_per_regulon = (auc_df > 0).sum(axis=1)
    # Keep regulons active in more than min_cells
    filtered_regulons = active_cells_per_regulon[active_cells_per_regulon > min_cells].index
    return filtered_regulons, active_cells_per_regulon

@st.cache_data
def calculate_mean_auc_by_celltype(auc_df, metadata, identity_column):
    """Calculate mean AUC per cell type with caching"""
    # Ensure cell order matches
    common_cells = list(set(auc_df.columns) & set(metadata.index))
    if len(common_cells) == 0:
        return None, None, "No matching cells found between AUC data and metadata!"
    
    auc_df_filtered = auc_df[common_cells]
    meta_df_filtered = metadata.loc[common_cells]
    
    # Calculate mean AUC per cell type
    cell_types = meta_df_filtered[identity_column].unique()
    mean_auc_by_celltype = pd.DataFrame(
        index=auc_df_filtered.index,
        columns=cell_types
    )
    
    for cell_type in cell_types:
        cells_of_type = meta_df_filtered[meta_df_filtered[identity_column] == cell_type].index
        cells_of_type = list(set(cells_of_type) & set(auc_df_filtered.columns))
        if len(cells_of_type) > 0:
            mean_auc_by_celltype[cell_type] = auc_df_filtered[cells_of_type].mean(axis=1)
    
    return mean_auc_by_celltype, meta_df_filtered, None

@st.cache_data
def calculate_z_scores(mean_auc_by_celltype):
    """Calculate Z-scores with caching"""
    # Row-wise normalization across cell types
    z_scores = mean_auc_by_celltype.apply(lambda row: (row - row.mean()) / row.std(), axis=1)
    z_scores = z_scores.replace([np.inf, -np.inf], np.nan).fillna(0)
    return z_scores

def regulon_specificity_scores(auc_mtx, cell_type_series):
    """
    Calculates the Regulon Specificity Scores (RSS). [doi: 10.1016/j.celrep.2018.10.045]
    Re-implemented to avoid numpy deprecation issue in pySCENIC.
    
    :param auc_mtx: The dataframe with the AUC values for all cells and regulons (n_cells x n_regulons).
    :param cell_type_series: A pandas Series object with cell identifiers as index and cell type labels as values.
    :return: A pandas dataframe with the RSS values (cell type x regulon).
    """
    cell_types = list(cell_type_series.unique())
    n_types = len(cell_types)
    regulons = list(auc_mtx.columns)
    n_regulons = len(regulons)
    rss_values = np.empty(shape=(n_types, n_regulons), dtype=np.float32)
    
    def rss(aucs, labels):
        # jensenshannon function provides distance which is the sqrt of the JS divergence.
        return 1.0 - jensenshannon(aucs / aucs.sum(), labels / labels.sum())
    
    for cidx, regulon_name in enumerate(regulons):
        for ridx, cell_type in enumerate(cell_types):
            rss_values[ridx, cidx] = rss(
                auc_mtx[regulon_name], (cell_type_series == cell_type).astype(int)
            )
    
    return pd.DataFrame(data=rss_values, index=cell_types, columns=regulons)

@st.cache_data
def calculate_rss(auc_df, metadata, identity_column):
    """Calculate RSS with caching"""
    # Ensure cell order matches
    common_cells = list(set(auc_df.columns) & set(metadata.index))
    auc_df_filtered = auc_df[common_cells]
    meta_df_filtered = metadata.loc[common_cells]
    
    cell_type_series = meta_df_filtered[identity_column]
    # regulon_specificity_scores returns cell types as index, regulons as columns
    # We need to transpose it to match scaled activity format (regulons as rows)
    rss_df = regulon_specificity_scores(auc_df_filtered.T, cell_type_series)
    return rss_df.T  # Transpose to have regulons as rows, cell types as columns

st.set_page_config(page_title="SCENIC Data Preparation", layout="wide")

st.title("SCENIC Data Preparation for SCENICviewer")
st.markdown("This app prepares AUC data from pySCENIC for visualization in SCENICviewer")

# Initialize session state
if 'auc_data' not in st.session_state:
    st.session_state.auc_data = None
if 'metadata' not in st.session_state:
    st.session_state.metadata = None
if 'selected_identity' not in st.session_state:
    st.session_state.selected_identity = None

# File upload section
st.header("1. Data Upload")
col1, col2 = st.columns(2)

with col1:
    st.subheader("AUC Data")
    auc_file = st.file_uploader("Upload AUC_per_cell.txt", type=['txt', 'tsv'])
    if auc_file is not None:
        try:
            st.session_state.auc_data = load_auc_data(auc_file, auc_file.name)
            st.success(f"AUC data loaded: {st.session_state.auc_data.shape[0]} regulons × {st.session_state.auc_data.shape[1]} cells")
            with st.expander("Preview AUC data"):
                st.dataframe(st.session_state.auc_data.iloc[:5, :5])
        except Exception as e:
            st.error(f"Error loading AUC data: {e}")

with col2:
    st.subheader("Metadata")
    metadata_option = st.radio("Choose metadata source:", ["meta.tsv file", "h5ad file"])
    
    if metadata_option == "meta.tsv file":
        meta_file = st.file_uploader("Upload meta.tsv", type=['tsv', 'txt'])
        if meta_file is not None:
            try:
                st.session_state.metadata = load_metadata_tsv(meta_file, meta_file.name)
                st.success(f"Metadata loaded: {st.session_state.metadata.shape[0]} cells × {st.session_state.metadata.shape[1]} features")
                with st.expander("Preview metadata"):
                    st.dataframe(st.session_state.metadata.head())
            except Exception as e:
                st.error(f"Error loading metadata: {e}")
    
    else:  # h5ad file
        h5ad_file = st.file_uploader("Upload h5ad file", type=['h5ad'])
        if h5ad_file is not None:
            try:
                st.session_state.metadata = load_metadata_h5ad(h5ad_file, h5ad_file.name)
                st.success(f"Metadata loaded from h5ad: {st.session_state.metadata.shape[0]} cells × {st.session_state.metadata.shape[1]} features")
                with st.expander("Preview metadata"):
                    st.dataframe(st.session_state.metadata.head())
            except Exception as e:
                st.error(f"Error loading h5ad file: {e}")

# Cell identity selection
if st.session_state.metadata is not None:
    st.header("2. Select Cell Identity")
    
    # Show available columns
    identity_columns = st.session_state.metadata.columns.tolist()
    st.session_state.selected_identity = st.selectbox(
        "Select identity column for cell type classification:",
        identity_columns,
        help="Choose the column that contains cell type annotations"
    )
    
    if st.session_state.selected_identity:
        # Get unique types and counts for the selected identity
        selected_column = st.session_state.metadata[st.session_state.selected_identity]
        unique_types = selected_column.unique()
        
        # Filter out NaN values if any
        unique_types = [x for x in unique_types if pd.notna(x)]
        
        st.info(f"Found {len(unique_types)} unique cell types in '{st.session_state.selected_identity}': {', '.join(map(str, unique_types[:5]))}{' ...' if len(unique_types) > 5 else ''}")
        
        # Show cell type distribution
        show_distribution = st.checkbox("Show cell type distribution", key=f"dist_{st.session_state.selected_identity}")
        if show_distribution:
            cell_counts = selected_column.value_counts().dropna()
            
            # Display as table and chart
            col1, col2 = st.columns([1, 2])
            with col1:
                st.dataframe(cell_counts.to_frame("Count"))
            with col2:
                # Use plotly for better bar chart visualization
                fig = px.bar(
                    x=cell_counts.index, 
                    y=cell_counts.values,
                    labels={'x': 'Cell Type', 'y': 'Count'},
                    title=f'Cell Type Distribution ({st.session_state.selected_identity})'
                )
                fig.update_layout(xaxis_tickangle=45)
                st.plotly_chart(fig, use_container_width=True)

# Calculate Z-scores and RSS
if st.session_state.auc_data is not None and st.session_state.metadata is not None and st.session_state.selected_identity:
    st.header("3. Calculate Z-scores and RSS")
    
    # Clear previous results if identity changed
    if 'last_identity' not in st.session_state:
        st.session_state.last_identity = None
    
    if st.session_state.last_identity != st.session_state.selected_identity:
        # Clear previous results when identity changes
        for key in ['z_scores', 'rss', 'mean_auc_by_celltype', 'z_scores_unfiltered', 'rss_unfiltered', 'mean_auc_by_celltype_unfiltered']:
            if key in st.session_state:
                del st.session_state[key]
        st.session_state.last_identity = st.session_state.selected_identity
    
    if st.button("Calculate", type="primary"):
        with st.spinner("Calculating Z-scores and RSS..."):
            try:
                # Get filtered regulons
                st.info("Identifying regulons active in >10 cells...")
                filtered_regulons, active_cells_count = filter_regulons_by_cells(st.session_state.auc_data, min_cells=10)
                st.info(f"Found {len(filtered_regulons)} filtered regulons (out of {len(st.session_state.auc_data.index)} total)")
                
                # Calculate mean AUC per cell type using cached function
                st.info("Calculating mean AUC per cell type...")
                mean_auc_by_celltype, meta_df_filtered, error = calculate_mean_auc_by_celltype(
                    st.session_state.auc_data, 
                    st.session_state.metadata, 
                    st.session_state.selected_identity
                )
                
                if error:
                    st.error(error)
                    st.stop()
                
                # Always calculate both unfiltered and filtered versions
                
                # Calculate unfiltered first
                st.info("Calculating Z-scores (unfiltered)...")
                z_scores_unfiltered = calculate_z_scores(mean_auc_by_celltype)
                
                st.info("Calculating RSS (unfiltered)...")
                rss_unfiltered = calculate_rss(
                    st.session_state.auc_data,
                    st.session_state.metadata,
                    st.session_state.selected_identity
                )
                
                st.session_state.z_scores_unfiltered = z_scores_unfiltered
                st.session_state.rss_unfiltered = rss_unfiltered
                st.session_state.mean_auc_by_celltype_unfiltered = mean_auc_by_celltype
                
                # Then calculate filtered
                mean_auc_filtered = mean_auc_by_celltype.loc[filtered_regulons]
                
                st.info("Calculating Z-scores (filtered)...")
                z_scores_filtered = calculate_z_scores(mean_auc_filtered)
                
                st.info("Calculating RSS (filtered)...")
                # For RSS, we need to filter the AUC data first
                auc_data_filtered = st.session_state.auc_data.loc[filtered_regulons]
                rss_filtered = calculate_rss(
                    auc_data_filtered,
                    st.session_state.metadata,
                    st.session_state.selected_identity
                )
                
                # Store filtered as main results (for display)
                st.session_state.z_scores = z_scores_filtered
                st.session_state.rss = rss_filtered
                st.session_state.mean_auc_by_celltype = mean_auc_filtered
                
                st.success("Calculations completed!")
                
                # Display results (filtered versions)
                col1, col2, col3 = st.columns(3)
                with col1:
                    st.subheader("Mean AUC by Cell Type (Filtered)")
                    st.dataframe(st.session_state.mean_auc_by_celltype.iloc[:10, :5])
                
                with col2:
                    st.subheader("Z-scores (Filtered)")
                    st.dataframe(st.session_state.z_scores.iloc[:10, :5])
                
                with col3:
                    st.subheader("RSS (Filtered)")
                    st.dataframe(st.session_state.rss.iloc[:10, :5])
                
            except Exception as e:
                st.error(f"Error during calculation: {e}")
                st.exception(e)

# Export results
if 'z_scores' in st.session_state and 'rss' in st.session_state:
    st.header("4. Export Results")
    
    # Create file names with identity
    identity_name = st.session_state.selected_identity.replace(' ', '_').replace('/', '_').replace('\\', '_')
    
    # Show filtering status
    if 'z_scores_unfiltered' in st.session_state and 'z_scores' in st.session_state:
        col1, col2 = st.columns(2)
        with col1:
            st.metric("Filtered regulons", st.session_state.z_scores.shape[0])
        with col2:
            st.metric("Unfiltered regulons", st.session_state.z_scores_unfiltered.shape[0])
    
    if st.button("Download All Results as ZIP", type="primary"):
        try:
            # Prepare data with proper column names
            scaled_activity = st.session_state.z_scores.copy()
            scaled_activity.columns = [str(col) for col in scaled_activity.columns]
            
            rss_export = st.session_state.rss.copy()
            rss_export.columns = [str(col) for col in rss_export.columns]
            
            mean_auc_export = st.session_state.mean_auc_by_celltype.copy()
            mean_auc_export.columns = [str(col) for col in mean_auc_export.columns]
            
            # Create file names with identity
            scaled_filename = f"scaled_regulon_activity_by_{identity_name}_FULL_TABLE_filtered.txt"
            rss_filename = f"rss_regulon_by_{identity_name}_FULL_TABLE_filtered.txt"
            mean_auc_filename = f"mean_auc_by_{identity_name}.txt"
            
            # Create ZIP file in memory
            zip_buffer = io.BytesIO()
            with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                # Add filtered versions (main files for SCENICviewer)
                zip_file.writestr(scaled_filename, scaled_activity.to_csv(sep='\t'))
                zip_file.writestr(rss_filename, rss_export.to_csv(sep='\t'))
                zip_file.writestr(mean_auc_filename, mean_auc_export.to_csv(sep='\t'))
                
                # Add unfiltered versions if available
                if 'z_scores_unfiltered' in st.session_state:
                    scaled_unfiltered = st.session_state.z_scores_unfiltered.copy()
                    scaled_unfiltered.columns = [str(col) for col in scaled_unfiltered.columns]
                    
                    rss_unfiltered = st.session_state.rss_unfiltered.copy()
                    rss_unfiltered.columns = [str(col) for col in rss_unfiltered.columns]
                    
                    mean_auc_unfiltered = st.session_state.mean_auc_by_celltype_unfiltered.copy()
                    mean_auc_unfiltered.columns = [str(col) for col in mean_auc_unfiltered.columns]
                    
                    # Create unfiltered file names
                    scaled_unfiltered_filename = f"scaled_regulon_activity_by_{identity_name}_FULL_TABLE_unfiltered.txt"
                    rss_unfiltered_filename = f"rss_regulon_by_{identity_name}_FULL_TABLE_unfiltered.txt"
                    mean_auc_unfiltered_filename = f"mean_auc_by_{identity_name}_unfiltered.txt"
                    
                    # Add unfiltered files
                    zip_file.writestr(scaled_unfiltered_filename, scaled_unfiltered.to_csv(sep='\t'))
                    zip_file.writestr(rss_unfiltered_filename, rss_unfiltered.to_csv(sep='\t'))
                    zip_file.writestr(mean_auc_unfiltered_filename, mean_auc_unfiltered.to_csv(sep='\t'))
                
                # Add metadata file with analysis info
                has_unfiltered = 'z_scores_unfiltered' in st.session_state
                metadata_info = f"""Analysis Information
===================
Identity column used: {st.session_state.selected_identity}
Number of filtered regulons: {scaled_activity.shape[0]}
{"Number of unfiltered regulons: " + str(st.session_state.z_scores_unfiltered.shape[0]) if has_unfiltered else ""}
Number of cell types: {scaled_activity.shape[1]}
Cell types: {', '.join(scaled_activity.columns)}

Filtering criteria: Regulons active in >10 cells

Files included (filtered - for SCENICviewer.py):
- {scaled_filename}: Z-scores for regulon activity
- {rss_filename}: Regulon Specificity Scores
- {mean_auc_filename}: Mean AUC values by cell type

{"Files included (unfiltered):" if has_unfiltered else ""}
{f"- {scaled_unfiltered_filename}: Z-scores for regulon activity" if has_unfiltered else ""}
{f"- {rss_unfiltered_filename}: Regulon Specificity Scores" if has_unfiltered else ""}
{f"- {mean_auc_unfiltered_filename}: Mean AUC values by cell type" if has_unfiltered else ""}

Processing method:
- Z-scores calculated using row-wise normalization (same as server.R tutorial)
- RSS calculated using re-implemented regulon_specificity_scores function
- Filtering applied after RSS calculation (same as server.R)
"""
                zip_file.writestr("analysis_info.txt", metadata_info)
            
            zip_buffer.seek(0)
            
            # Download button for ZIP file
            st.download_button(
                label="📦 Download ZIP file",
                data=zip_buffer.getvalue(),
                file_name=f"scenic_results_by_{identity_name}.zip",
                mime="application/zip"
            )
            
            st.success(f"ZIP file prepared with {scaled_activity.shape[0]} regulons × {scaled_activity.shape[1]} cell types")
            
            # Show file contents
            with st.expander("📁 Files included in ZIP"):
                st.write("**Filtered files (for SCENICviewer.py):**")
                st.write(f"- {scaled_filename}")
                st.write(f"- {rss_filename}")
                st.write(f"- {mean_auc_filename}")
                
                if has_unfiltered:
                    st.write("\n**Unfiltered files:**")
                    st.write(f"- {scaled_unfiltered_filename}")
                    st.write(f"- {rss_unfiltered_filename}")
                    st.write(f"- {mean_auc_unfiltered_filename}")
                
                st.write("\n**Metadata:**")
                st.write("- analysis_info.txt")
                
        except Exception as e:
            st.error(f"Error creating ZIP file: {e}")
            st.exception(e)

# Notes and information
with st.expander("ℹ️ About this tool"):
    st.markdown("""
    ### Data Processing Steps:
    
    1. **AUC Data**: Load the AUC matrix from pySCENIC (regulons × cells)
    2. **Metadata**: Load cell annotations from meta.tsv or h5ad file
    3. **Z-score Calculation** (same as server.R tutorial):
       - Calculate mean AUC for each regulon in each cell type
       - Row-wise normalization: Z-score = (value - row_mean) / row_std
       - Each regulon is normalized across cell types (mean=0, std=1)
    4. **RSS Calculation**: 
       - Uses pySCENIC's regulon_specificity_scores function
       - Measures how specific each regulon is to each cell type
    
    ### Output Files:
    - **scaled_regulon_activity_by_cell_type_FULL_TABLE_filtered.txt**: Z-scores for regulon activity
    - **rss_regulon_by_cell_type_FULL_TABLE_filtered.txt**: Regulon Specificity Scores
    
    These files are compatible with SCENICviewer.py for visualization.
    """)

# Validation notes
with st.expander("⚠️ Important Notes"):
    st.markdown("""
    ### Compatibility with server.R:
    
    1. **Z-score Calculation**: 
       - Uses the same row-wise normalization as server.R tutorial
       - Each regulon is normalized across cell types: (value - mean) / std
    
    2. **RSS Calculation**:
       - Uses pySCENIC's built-in function, which should match the standard pipeline
    
    3. **Data Format**:
       - Ensure your AUC data has regulons as rows and cells as columns
       - Cell barcodes must match between AUC data and metadata
    
    ### Recommendations:
    - Verify that your AUC data was generated using standard pySCENIC steps
    - Check that cell type annotations are biologically meaningful
    - Consider filtering out regulons with very low activity before analysis
    """)