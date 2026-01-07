import scanpy as sc
import anndata as ad
import streamlit as st
import pandas as pd
import numpy as np
import os
import hashlib
import celltypist
from celltypist import models
import matplotlib.pyplot as plt
import io
import tempfile
import time

st.set_page_config(page_title="CellTypist", page_icon="🔬")

# Set custom model directory (relative to streamlit directory)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(SCRIPT_DIR)  # /home/ichiro/streamlit
MODEL_DIR = os.path.join(BASE_DIR, "db", "celltypist", "data", "models")

st.markdown("## CellTypist - Automated Cell Type Annotation")
st.markdown("""
[CellTypist](https://www.celltypist.org/) is a tool for automated cell type annotation
using machine learning classifiers trained on large-scale single-cell RNA-seq data.
""")

# Get available models from local directory
@st.cache_data
def get_model_list():
    """Get list of available CellTypist models from local directory"""
    model_info = []
    if os.path.exists(MODEL_DIR):
        all_models = [f for f in os.listdir(MODEL_DIR) if f.endswith('.pkl')]
        for m in sorted(all_models):
            name = m.replace('.pkl', '')
            # Categorize by species
            name_lower = name.lower()
            if 'Mouse' in name or 'mouse' in name_lower:
                species = 'Mouse'
            elif ('Human' in name or 'human' in name_lower or
                  'COVID' in name or 'PBMC' in name or 'Immune_All' in name or
                  'cHSPCs' in name or  # Human circulating HSPCs
                  name.startswith('Cells_') or  # Human cell atlases
                  name.startswith('Nuclei_') or  # Human nuclei atlases
                  name.startswith('Healthy_') or  # Healthy human tissues
                  'Fetal' in name or 'Developing' in name):  # Human developmental
                species = 'Human'
            elif 'Macaque' in name or 'Pig' in name:
                species = 'Other (Primate/Pig)'
            else:
                species = 'Other'
            model_info.append({'model': m, 'name': name, 'species': species})
    return pd.DataFrame(model_info)

def get_file_hash(file_bytes):
    """Get MD5 hash of file content"""
    return hashlib.md5(file_bytes).hexdigest()

@st.cache_data
def read_h5ad_cached(_file_bytes, file_hash):
    """Read h5ad file from bytes with caching based on content hash"""
    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_file:
        tmp_file.write(_file_bytes)
        tmp_path = tmp_file.name
    adata = sc.read_h5ad(tmp_path)
    os.unlink(tmp_path)
    return adata

@st.cache_resource
def load_model(model_name, convert_to_mouse=False):
    """Load CellTypist model from local directory, optionally convert to mouse"""
    model_path = os.path.join(MODEL_DIR, model_name)
    model = models.Model.load(model=model_path)
    if convert_to_mouse:
        model.convert()  # Convert human gene names to mouse orthologs
    return model

def run_celltypist(_adata, model_name, mode='best match', majority_voting=True, convert_to_mouse=False):
    """Run CellTypist annotation"""
    model = load_model(model_name, convert_to_mouse=convert_to_mouse)
    predictions = celltypist.annotate(
        _adata,
        model=model,
        mode=mode,
        majority_voting=majority_voting
    )
    return predictions

# Sidebar for model selection
st.sidebar.header("Model Selection")

model_df = get_model_list()

# Species filter
species_filter = st.sidebar.multiselect(
    "Filter by species",
    options=['Human', 'Mouse', 'Other (Primate/Pig)', 'Other'],
    default=['Human']
)

filtered_models = model_df[model_df['species'].isin(species_filter)]

# Set default to Immune_All_Low if available
default_model = 'Immune_All_Low.pkl'
model_options = filtered_models['model'].tolist()
if default_model in model_options:
    default_index = model_options.index(default_model)
else:
    default_index = 0

selected_model = st.sidebar.selectbox(
    "Select model",
    options=model_options,
    index=default_index,
    format_func=lambda x: x.replace('.pkl', '')
)

if selected_model:
    st.sidebar.markdown(f"**Selected:** {selected_model.replace('.pkl', '')}")

# Check if selected model is Human (can be converted to Mouse)
model_info = model_df[model_df['model'] == selected_model]
is_human_model = len(model_info) > 0 and model_info.iloc[0]['species'] == 'Human'

# Main content
st.markdown("### Upload Data")

uploaded_file = st.file_uploader(
    "Upload h5ad file (AnnData format)",
    type=['h5ad'],
    help="Upload a preprocessed single-cell RNA-seq dataset in h5ad format"
)

if uploaded_file is not None:
    with st.spinner("Loading data..."):
        # Read h5ad with caching based on file content hash
        file_bytes = uploaded_file.getvalue()
        file_hash = get_file_hash(file_bytes)
        adata = read_h5ad_cached(file_bytes, file_hash)

    st.success(f"Loaded: {adata.n_obs} cells x {adata.n_vars} genes")

    # Detect species from gene names
    def detect_species_from_genes(adata):
        """Detect species by checking gene name patterns"""
        genes = adata.var_names[:100].tolist()  # Check first 100 genes
        # Mouse genes: Title case (e.g., Gapdh, Actb)
        # Human genes: UPPERCASE (e.g., GAPDH, ACTB)
        mouse_markers = ['Gapdh', 'Actb', 'Ptprc', 'Cd3e', 'Cd4', 'Cd8a']
        human_markers = ['GAPDH', 'ACTB', 'PTPRC', 'CD3E', 'CD4', 'CD8A']

        mouse_count = sum(1 for g in mouse_markers if g in genes)
        human_count = sum(1 for g in human_markers if g in genes)

        if mouse_count > human_count:
            return 'Mouse'
        elif human_count > mouse_count:
            return 'Human'
        else:
            # Check case pattern
            upper_count = sum(1 for g in genes if g.isupper())
            if upper_count > len(genes) / 2:
                return 'Human'
            else:
                return 'Mouse'

    detected_species = detect_species_from_genes(adata)

    # Show basic info
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("**Cell metadata columns:**")
        st.write(list(adata.obs.columns))
    with col2:
        st.markdown("**Available embeddings:**")
        st.write(list(adata.obsm.keys()) if adata.obsm else "None")

    st.info(f"**Detected species from gene names:** {detected_species}")

    # CellTypist options
    st.markdown("### Annotation Settings")

    # Model conversion option (Human -> Mouse)
    if is_human_model:
        # Auto-check if mouse data detected
        default_convert = (detected_species == 'Mouse')
        convert_to_mouse = st.checkbox(
            "Convert Human model to Mouse",
            value=default_convert,
            help="Convert human gene names to mouse orthologs. Recommended when analyzing mouse data with human-trained models."
        )
        if convert_to_mouse:
            st.warning("⚠️ Cross-species results should be interpreted with caution.")
    else:
        convert_to_mouse = False
        if detected_species == 'Mouse':
            st.success("✓ Mouse model selected for mouse data")

    # PCA selection for majority voting
    pca_keys = [k for k in adata.obsm.keys() if 'pca' in k.lower()]

    if len(pca_keys) > 0:
        col1, col2, col3 = st.columns(3)
        with col1:
            selected_pca = st.selectbox(
                "Select PCA for majority voting",
                options=pca_keys,
                help="Select which PCA embedding to use for neighborhood calculation"
            )
            # Show PCA dimensions
            n_pcs = adata.obsm[selected_pca].shape[1]
            st.caption(f"Dimensions: {n_pcs} PCs")
            if n_pcs < 50:
                st.warning(f"Warning: {n_pcs} PCs < 50. Majority voting may fail.")
    else:
        selected_pca = None
        st.warning("No PCA found in data. Majority voting will be disabled.")

    col1, col2 = st.columns(2)
    with col1:
        mode = st.selectbox(
            "Prediction mode",
            options=['best match', 'prob match'],
            help="'best match': assigns the most likely cell type. 'prob match': assigns cell type based on probability threshold."
        )
    with col2:
        majority_voting_disabled = selected_pca is None
        majority_voting = st.checkbox(
            "Majority voting",
            value=not majority_voting_disabled,
            disabled=majority_voting_disabled,
            help="Refine predictions by local majority voting within cell neighborhoods (requires PCA with >= 50 components)"
        )

    # Create a unique key for this analysis
    analysis_key = f"{file_hash}_{selected_model}_{mode}_{majority_voting}_{convert_to_mouse}"

    # Run annotation
    if st.button("Run CellTypist", type="primary"):
        # Copy selected PCA to X_pca if needed
        if majority_voting and selected_pca and selected_pca != 'X_pca':
            adata.obsm['X_pca'] = adata.obsm[selected_pca].copy()
            st.info(f"Using {selected_pca} as X_pca for majority voting")
        with st.spinner(f"Running CellTypist with {selected_model}..."):
            start_time = time.time()

            try:
                predictions = run_celltypist(adata, selected_model, mode, majority_voting, convert_to_mouse)
                elapsed = time.time() - start_time

                # Get results
                adata_result = predictions.to_adata()

                # Store in session state with analysis key
                st.session_state['celltypist_key'] = analysis_key
                st.session_state['celltypist_result'] = adata_result
                st.session_state['celltypist_predictions'] = predictions
                st.session_state['celltypist_elapsed'] = elapsed
                st.session_state['celltypist_convert'] = convert_to_mouse
                st.session_state['celltypist_model'] = selected_model
                st.session_state['celltypist_majority'] = majority_voting

            except Exception as e:
                st.error(f"Error during annotation: {str(e)}")
                st.exception(e)

    # Display results if available in session state
    if ('celltypist_key' in st.session_state and
        st.session_state.get('celltypist_key') == analysis_key and
        'celltypist_result' in st.session_state):

        adata_result = st.session_state['celltypist_result']
        predictions = st.session_state['celltypist_predictions']
        elapsed = st.session_state.get('celltypist_elapsed', 0)
        convert_to_mouse = st.session_state.get('celltypist_convert', False)
        majority_voting = st.session_state.get('celltypist_majority', True)

        if convert_to_mouse:
            st.success(f"Annotation completed in {elapsed:.1f} seconds (Human→Mouse converted)")
        else:
            st.success(f"Annotation completed in {elapsed:.1f} seconds")

        # Show prediction summary
        st.markdown("### Results Summary")

        # Cell type distribution
        if majority_voting:
            pred_col = 'majority_voting'
        else:
            pred_col = 'predicted_labels'

        cell_counts = adata_result.obs[pred_col].value_counts()

        col1, col2 = st.columns([1, 2])
        with col1:
            st.markdown("**Cell type counts:**")
            st.dataframe(cell_counts.reset_index().rename(
                columns={'index': 'Cell Type', pred_col: 'Count'}
            ))

        with col2:
            # Bar plot
            fig, ax = plt.subplots(figsize=(10, max(4, len(cell_counts) * 0.3)))
            cell_counts.sort_values().plot(kind='barh', ax=ax)
            ax.set_xlabel('Number of cells')
            ax.set_ylabel('Cell type')
            ax.set_title('Cell Type Distribution')
            plt.tight_layout()
            st.pyplot(fig)
            plt.close()

        # UMAP visualization
        st.markdown("### UMAP Visualization")

        # Find all UMAP-like embeddings
        umap_keys = [k for k in adata_result.obsm.keys() if 'umap' in k.lower()]

        if len(umap_keys) > 0:
            selected_umap = st.selectbox(
                "Select UMAP embedding",
                options=umap_keys,
                index=0 if 'X_umap' not in umap_keys else umap_keys.index('X_umap')
            )

            # Copy selected UMAP to X_umap for scanpy plotting
            if selected_umap != 'X_umap':
                adata_result.obsm['X_umap'] = adata_result.obsm[selected_umap].copy()

            fig, axes = plt.subplots(1, 2, figsize=(16, 6))

            # Original clusters if available
            if 'leiden' in adata_result.obs:
                sc.pl.umap(adata_result, color='leiden', ax=axes[0], show=False, title='Original Clusters')
            elif 'louvain' in adata_result.obs:
                sc.pl.umap(adata_result, color='louvain', ax=axes[0], show=False, title='Original Clusters')
            elif 'seurat_clusters' in adata_result.obs:
                sc.pl.umap(adata_result, color='seurat_clusters', ax=axes[0], show=False, title='Original Clusters')
            else:
                axes[0].text(0.5, 0.5, 'No cluster info', ha='center', va='center', transform=axes[0].transAxes)
                axes[0].set_title('Original Clusters')

            # CellTypist predictions
            sc.pl.umap(adata_result, color=pred_col, ax=axes[1], show=False,
                      title=f'CellTypist ({selected_model.replace(".pkl", "")})', legend_loc='on data', legend_fontsize=6)

            plt.tight_layout()
            st.pyplot(fig)
            plt.close()
        else:
            st.warning("No UMAP embedding found (X_umap, X_*.umap, etc.). Run sc.tl.umap() first for visualization.")

        # Decision matrix (probability)
        st.markdown("### Prediction Probabilities")

        prob_matrix = predictions.probability_matrix

        with st.expander("Show probability matrix (top cells)"):
            st.dataframe(prob_matrix.head(100))

        # Download results
        st.markdown("### Download Results")

        col1, col2, col3 = st.columns(3)

        with col1:
            # Download predictions as TSV
            # Build column list without duplicates
            cols = [pred_col]
            if pred_col != 'predicted_labels' and 'predicted_labels' in adata_result.obs.columns:
                cols.append('predicted_labels')
            if 'conf_score' in adata_result.obs.columns:
                cols.append('conf_score')
            pred_df = adata_result.obs[cols]
            pred_df.index.name = 'cell_barcode'  # Add index name for SCALA compatibility
            tsv = pred_df.to_csv(sep='\t').encode('utf-8')
            st.download_button(
                label="Download predictions (TSV)",
                data=tsv,
                file_name="celltypist_predictions.tsv",
                mime="text/tab-separated-values"
            )

        with col2:
            # Download probability matrix
            prob_tsv = prob_matrix.to_csv(sep='\t').encode('utf-8')
            st.download_button(
                label="Download probabilities (TSV)",
                data=prob_tsv,
                file_name="celltypist_probabilities.tsv",
                mime="text/tab-separated-values"
            )

        with col3:
            # Download annotated h5ad
            with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp_out:
                adata_result.write_h5ad(tmp_out.name)
                with open(tmp_out.name, 'rb') as f:
                    h5ad_data = f.read()
                os.unlink(tmp_out.name)

            st.download_button(
                label="Download annotated h5ad",
                data=h5ad_data,
                file_name="celltypist_annotated.h5ad",
                mime="application/octet-stream"
            )

else:
    st.info("Please upload an h5ad file to begin annotation.")

    # Show model information
    st.markdown("### Available Models")
    st.markdown(f"**Total models:** {len(model_df)}")

    with st.expander("Show all models"):
        for species in ['Human', 'Mouse', 'Other (Primate/Pig)', 'Other']:
            species_models = model_df[model_df['species'] == species]
            if len(species_models) > 0:
                st.markdown(f"**{species}** ({len(species_models)} models)")
                for _, row in species_models.iterrows():
                    st.markdown(f"- {row['name']}")

# Footer
st.markdown("---")
st.markdown("""
**Reference:** Domínguez Conde, C., et al. (2022). Cross-tissue immune cell analysis reveals tissue-specific
features in humans. Science, 376(6594), eabl5197.
""")
