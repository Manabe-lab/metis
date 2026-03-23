import scanpy as sc
import anndata as ad
import numpy as np
import streamlit as st
import pandas as pd
import time
import os
from helper_func import clear_old_directories
from helper_func import clear_old_files
import matplotlib.pyplot as plt
import io
import seaborn as sns
import SEACells
import tempfile
import shutil
import time
import cupyx.scipy.sparse


# st.set_option('deprecation.showPyplotGlobalUse', False)

st.set_page_config(page_title="Metacell aggregation by SEACells", page_icon="💬")
st.markdown("### Metacell generation by SEACells")
st.write("Stores the sum of raw counts from cells in each metacell into adata.X / adata.raw.X.")
st.write("The CPM normalization option (ON by default) applies CPM normalization to adata.X.")

@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    return adata

def generate_group_name(row):
    # Use "all" when cell_type is None
    cell_type_val = "all" if cell_type is None else row[cell_type]

    if sample_id == 'single_sample':
        key = cell_type_val
    else:
        key = (row[sample_id], cell_type_val)
    if key not in group_counter:
        group_counter[key] = 0
    group_counter[key] += 1
    if sample_id == 'single_sample':
        return f"{cell_type_val}_{group_counter[key]}"
    else:
        return f"{cell_type_val}_{row[sample_id]}_{group_counter[key]}"

def calc_df(adata, adata_df):
    df = pd.DataFrame(adata.raw.X.T)
    df.columns = adata_df['SEACell_group'].to_list()
    df.index = adata.var.index.to_list()
    df = df.sort_index(axis=1)
    return df

if "seacell_temp_dir" not in st.session_state or not os.path.exists('temp/seacell'):
    st.session_state.seacell_temp_dir = True
    seacell_temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
    os.mkdir(seacell_temp_dir)
    os.mkdir(seacell_temp_dir + '/seacell')
    st.session_state.seacell_temp_dir = seacell_temp_dir
else:
    seacell_temp_dir = st.session_state.seacell_temp_dir


uploaded_file = st.file_uploader("Upload a h5ad file", type=['h5ad'])


#---------------

if uploaded_file is not None:
    adata = read_h5ad(uploaded_file)
    st.write("Uploaded data")
    st.write(adata)

    # ========================================
    # Raw Counts Layer Selection
    # ========================================
    st.subheader("📊 Raw Counts Layer Selection")
    st.markdown("SEACells uses **raw counts** to aggregate metacells.")

    # Check current .X statistics
    if hasattr(adata.X, 'toarray'):
        x_sample = adata.X[:100, :100].toarray()
    else:
        x_sample = adata.X[:100, :100]
    x_max = float(np.max(x_sample))
    x_has_decimals = not np.allclose(x_sample, np.round(x_sample))

    # Check .raw.X if available
    raw_available = hasattr(adata, 'raw') and adata.raw is not None
    if raw_available:
        if hasattr(adata.raw.X, 'toarray'):
            raw_sample = adata.raw.X[:100, :100].toarray()
        else:
            raw_sample = adata.raw.X[:100, :100]
        raw_max = float(np.max(raw_sample))
        raw_has_decimals = not np.allclose(raw_sample, np.round(raw_sample))

    # Build layer options
    layer_options = []
    if raw_available:
        layer_options.append("adata.raw.X (default)")
    layer_options.append("adata.X")
    if adata.layers:
        for layer_name in adata.layers.keys():
            layer_options.append(f"adata.layers['{layer_name}']")

    # Show layer info
    with st.expander("Layer statistics", expanded=False):
        if raw_available:
            st.write(f"**adata.raw.X**: max={raw_max:.2f}, has_decimals={raw_has_decimals}")
        st.write(f"**adata.X**: max={x_max:.2f}, has_decimals={x_has_decimals}")
        if adata.layers:
            for layer_name in adata.layers.keys():
                layer_data = adata.layers[layer_name]
                if hasattr(layer_data, 'toarray'):
                    layer_sample = layer_data[:100, :100].toarray()
                else:
                    layer_sample = layer_data[:100, :100]
                layer_max = float(np.max(layer_sample))
                layer_decimals = not np.allclose(layer_sample, np.round(layer_sample))
                st.write(f"**adata.layers['{layer_name}']**: max={layer_max:.2f}, has_decimals={layer_decimals}")

    selected_raw_layer = st.selectbox(
        "Select raw counts source:",
        options=layer_options,
        index=0,
        help="Select the layer containing raw counts for SEACells metacell aggregation"
    )

    # Apply layer selection - copy to adata.raw if needed
    if selected_raw_layer == "adata.X":
        if x_has_decimals:
            st.warning("⚠️ adata.X contains decimal values. This data may already be normalized.")
        # Copy X to raw
        adata.raw = adata.copy()
        st.info("✓ Copied adata.X to adata.raw")
    elif "adata.layers['" in selected_raw_layer:
        layer_name = selected_raw_layer.split("'")[1]
        # Create raw from layer
        adata_raw = adata.copy()
        adata_raw.X = adata.layers[layer_name].copy()
        adata.raw = adata_raw
        st.info(f"✓ Copied adata.layers['{layer_name}'] to adata.raw")
    elif selected_raw_layer == "adata.raw.X (default)":
        if not raw_available:
            st.error("adata.raw does not exist. Please select a different layer.")
            st.stop()
        st.info("✓ Using adata.raw.X")

    st.write("Snippet of the count matrix (raw):")
    st.write(adata.raw.X[:5,:5])
    st.write(f"The number of cells: {len(adata)}")

    meta_data = adata.obs.columns.to_list()

    meta_data = list(filter(lambda x: x not in ['nCount_RNA', 'nFeature_RNA','percent.mt', 'Cell_id'], meta_data))
    obsm_keys = [x.replace('X_','') for x in list(adata.obsm)]
    if len(obsm_keys) == 0:
        st.markdown("#### No dimensionality reduced data!!! Please upload a collect data set. At least the reduced data that was used to determine clusters (e.g., PCA) are needed.")
        time.sleep(10)
        st.rerun()
    umap_keys = list(filter(lambda x: 'umap' in x.lower(), obsm_keys))
    umap_key = umap_keys[0] if umap_keys else None
    if umap_key is not None:
        umap_index = umap_keys.index(umap_key)
    else:
        umap_index = 0

    reduction_vis = umap_key
    reduction_use = obsm_keys[0]
    reduction_subset = 'X_' + reduction_use + '_sub'


    split_sel = st.checkbox(
        'Split by samples?',
        value=True,
        help="Split cells by sample for processing. When OFF, all cells are processed together, and cells from different samples may be mixed in the same metacell."
    )
    st.caption("⚠️ If you want to aggregate metacells per sample, make sure to enable this option.")

    treat_as_single_type = st.checkbox(
        'Treat all cells as single type',
        value=False,
        help="When enabled, treats all cells as the same type without distinguishing cell types."
    )

    wgcna = st.checkbox("Generate single matrix file for downsteram analysis (e.g., WGCNA)?", value = True)

    normalize_cpm = st.checkbox(
        "Apply CPM normalization to adata.X",
        value=True,
        help="Converts adata.X from raw count sum to CPM normalization (target_sum=1,000,000). "
             "Recommended for downstream analyses such as COMPASS, scFEA, etc. adata.raw.X retains the raw count sum."
    )

    sample_id = "single_sample"

    with st.form("Basic settings:"):
        if split_sel:
            sample_id = st.selectbox("Split by:", meta_data)

        if treat_as_single_type:
            # When treating cells as single type, use dummy cell_type column
            cell_type = None
        else:
            cell_type = st.selectbox("Identity of cell types:", meta_data)

        reduction_use = st.selectbox("Reduction used for clustering (e.g., PCA, MNN)", obsm_keys)
        reduction_subset = 'X_' + reduction_use + '_sub'

        n_dim = st.number_input('Number of dimensions used for clustering', min_value =1, value=30)

        obsm_keys = [x.replace('X_','') for x in list(adata.obsm)]

        reduction_vis = st.selectbox("Reducition to visualize (e.g. umap):",umap_keys, index = umap_index)

        fracion_n_SEACs = st.number_input('n_SEACells: approximate the number of cells to generate', min_value =0, value=round(len(adata)/75))
        st.write("As a rule of thumb, choosing one metacell for every 75 single-cells.")
        st.write("This is the total number across all samples.")
        n_waypoint_eigs = st.number_input('n_waypoint_eigs:',  value=10)

        preserve_small_samples = st.checkbox(
            "Preserve small samples (minimum 1 metacell)",
            value=False,
            help="Creates at least 1 metacell even for samples with few cells. When OFF, samples with n_SEACells=0 are skipped."
        )

        use_gpu = st.checkbox(
            "Use GPU (cupy)",
            value=True,
            help="Use GPU for acceleration. Disable if CUBLAS errors occur."
        )

        submitted_basic = st.form_submit_button("Set the parameters and run")


    if submitted_basic:
        fig, ax = plt.subplots()

        fig, ax = plt.subplots(figsize=(4, 4))
        if cell_type is not None:
            sc.pl.scatter(adata, basis=reduction_vis, color=cell_type, frameon=True, ax=ax, show=False)
        else:
            sc.pl.scatter(adata, basis=reduction_vis, frameon=True, ax=ax, show=False)
        st.pyplot(fig)

        sc.pp.highly_variable_genes(adata, n_top_genes=1500)


        adata.obsm[reduction_subset] = adata.obsm['X_' + reduction_use][:,:n_dim]



        # Create a dictionary to store split AnnData objects
        adata_dict = {}

        if split_sel:

            # Get unique values in 'orig.ident'
            unique_idents = adata.obs[sample_id].unique()

            for ident in unique_idents:
                adata_dict[ident] = adata[adata.obs[sample_id] == ident].copy()

        else:
            adata_dict["all"] = adata


        dict_keys = list(adata_dict.keys())
        SEACell_adata = {}

        total_cell_num = len(adata)

        for i in dict_keys:
            st.markdown("### " + i)
            st.write(f"{len(adata_dict[i])} cells.")
            if len(adata_dict[i]) == 0:
                st.markdown("#### No cells. Skipe this sample.")
                continue

            # Recalculate fracion_n_SEACs based on total sample count

            ## 75 cell:bin
#            n_cells = round(len(adata_dict[i].obs) / fracion_n_SEACs)
            n_SEACells_raw = round(len(adata_dict[i]) / total_cell_num * fracion_n_SEACs)

            if preserve_small_samples:
                n_SEACells = max(1, n_SEACells_raw)
                if n_SEACells_raw == 0:
                    st.warning(f"⚠️ Few cells, setting n_SEACells=1 (calculated value: {n_SEACells_raw})")
            else:
                n_SEACells = n_SEACells_raw
                if n_SEACells == 0:
                    st.warning(f"⚠️ Skipping this sample because n_SEACells=0")
                    continue

            st.write(f"nSEACells for this sample: {n_SEACells}")

            ## Core parameters
            build_kernel_on =reduction_subset  # key in adata.obsm to use for computing metacells
                                      # This would be replaced by 'X_svd' for ATAC data

            ## adataditional parameters
             # Number of eigenvalues to consider when initializing metacells
            model = SEACells.core.SEACells(adata_dict[i],
                          build_kernel_on=build_kernel_on,
                          n_SEACells=n_SEACells,
                          n_waypoint_eigs=n_waypoint_eigs,
                          convergence_epsilon = 1e-5,
                          use_gpu = use_gpu)


            try:
                model.construct_kernel_matrix()
            except Exception as e:
                st.error(f"Error during kernel matrix construction: {e}")
                if "CUBLAS" in str(e) or "cupy" in str(e).lower() or "cuda" in str(e).lower():
                    st.warning("⚠️ GPU error occurred. Please uncheck 'Use GPU (cupy)' and re-run in CPU mode.")
                    st.stop()
                raise
            M = model.kernel_matrix
            st.write(f"M shape: {M.shape}")

            # Initialize archetypes
            try:
                model.initialize_archetypes()
            except Exception as e:
                st.error(f"Error during archetype initialization: {e}")
                st.warning(f"Skipping this sample ({len(adata_dict[i])} cells). Please set `n_waypoint_eigs` (current: {n_waypoint_eigs}) to a smaller value and re-run.")
                del adata_dict[i]
                continue

            # Clear any existing matplotlib figures
            plt.clf()
            plt.close('all')

            # Prepare a buffer to save the plot
            buf = io.BytesIO()

            # Plot the initialization
            st.write("Initial points. Make sure they are spread across phenotypic space.")
            SEACells.plot.plot_initialization(
                adata_dict[i],
                model,
                plot_basis='X_' + reduction_vis,
                save_as=buf,  # Save to our buffer
                show=False    # Don't show the plot (which would close it)
            )

            # Reset buffer position
            buf.seek(0)
            # Display the image in Streamlit
            st.image(buf)
            # Clear any remaining matplotlib state
            plt.clf()
            plt.close('all')
            try:
                model.fit(min_iter=10, max_iter=50)
            except Exception as e:
                st.error(f"Error during model fitting: {e}")
                if "CUBLAS" in str(e) or "cupy" in str(e).lower() or "cuda" in str(e).lower():
                    st.warning("⚠️ GPU error occurred. Please uncheck 'Use GPU (cupy)' and re-run in CPU mode.")
                    st.stop()
                else:
                    st.write('Please add "import cupyx.scipy.sparse" and "import cupy as cp" to core.py of SEACells.')


            # Check for convergence
            model.plot_convergence(save_as=buf,  show=False)
            buf.seek(0)
            st.image(buf)
            plt.clf()
            plt.close('all')


            try:
                model.fit(min_iter=10, max_iter=50)
            except Exception as e:
                st.error(f"Error during final model fitting: {e}")
                if "CUBLAS" in str(e) or "cupy" in str(e).lower() or "cuda" in str(e).lower():
                    st.warning("⚠️ GPU error occurred. Please uncheck 'Use GPU (cupy)' and re-run in CPU mode.")
                    st.stop()
                raise


            adata_dict[i].obs[['SEACell']].head()

            SEACell_adata[i] = SEACells.core.summarize_by_SEACell(adata_dict[i], SEACells_label='SEACell', summarize_layer='raw')

            SEACells.plot.plot_2D(adata_dict[i], key='X_'+reduction_vis, colour_metacells=True)

            SEACells.plot.plot_SEACell_sizes(adata_dict[i], bins=5, save_as=buf,  show=False)
            buf.seek(0)
            st.image(buf)
            plt.clf()
            plt.close('all')


            # Compute celltype purity only when cell_type is specified
            if cell_type is not None:
                SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata_dict[i], cell_type)

                # Create a new figure with 3 subplots in a row
                fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

                # Plot 1: Celltype Purity
                sns.boxplot(data=SEACell_purity, y=cell_type + '_purity', ax=ax1)
                ax1.set_title('Celltype Purity')
                sns.despine(ax=ax1)
            else:
                # Only two plots when cell_type is not specified
                fig, (ax2, ax3) = plt.subplots(1, 2, figsize=(8, 4))

            # Plot 2: Compactness
            compactness = SEACells.evaluate.compactness(adata_dict[i], 'X_'+reduction_use)
            sns.boxplot(data=compactness, y='compactness', ax=ax2)
            ax2.set_title('Compactness')
            sns.despine(ax=ax2)

            # Plot 3: Separation
            separation = SEACells.evaluate.separation(adata_dict[i], 'X_'+reduction_use, nth_nbr=1)
            sns.boxplot(data=separation, y='separation', ax=ax3)
            ax3.set_title('Separation')
            sns.despine(ax=ax3)

            # Adjust the layout and show the plot
            plt.tight_layout()
            st.pyplot(fig)
        if split_sel:
        # Merge AnnData objects
            adata_merged = ad.concat(
                adata_dict,
                join='outer',
                merge='same',
                label='sample'
            )
            # Rename SEAcell names: {sample_id}_{cell_type}_SEACell-{n}
            # Convert to string to avoid categorical column assignment issues
            seacell_series = adata_merged.obs['SEACell'].astype(str)
            sample_series = adata_merged.obs[sample_id].astype(str)

            # Get cell_type (empty string if None)
            if cell_type is not None:
                celltype_series = adata_merged.obs[cell_type].astype(str)
            else:
                celltype_series = pd.Series([''] * len(seacell_series), index=seacell_series.index)

            new_seacell_names = []
            for i in range(len(seacell_series)):
                sample_val = sample_series.iloc[i]
                celltype_val = celltype_series.iloc[i]
                seacell_val = seacell_series.iloc[i]

                # Extract the SEACell-{n} part
                if 'SEACell-' in seacell_val:
                    seacell_num = seacell_val.split('SEACell-')[-1]
                    seacell_suffix = f"SEACell-{seacell_num}"
                else:
                    seacell_suffix = seacell_val

                # Build name: {sample_id}_{cell_type}_SEACell-{n}
                if cell_type is not None and celltype_val:
                    new_name = f"{sample_val}_{celltype_val}_{seacell_suffix}"
                else:
                    new_name = f"{sample_val}_{seacell_suffix}"

                new_seacell_names.append(new_name)

            adata_merged.obs['SEACell'] = new_seacell_names

        else:
            adata_merged = adata_dict["all"]

            # Without split: {cell_type}_SEACell-{n} (no sample_id)
            seacell_series = adata_merged.obs['SEACell'].astype(str)

            if cell_type is not None:
                celltype_series = adata_merged.obs[cell_type].astype(str)

                new_seacell_names = []
                for i in range(len(seacell_series)):
                    celltype_val = celltype_series.iloc[i]
                    seacell_val = seacell_series.iloc[i]

                    # Extract the SEACell-{n} part
                    if 'SEACell-' in seacell_val:
                        seacell_num = seacell_val.split('SEACell-')[-1]
                        seacell_suffix = f"SEACell-{seacell_num}"
                    else:
                        seacell_suffix = seacell_val

                    # Build name: {cell_type}_SEACell-{n}
                    new_name = f"{celltype_val}_{seacell_suffix}"
                    new_seacell_names.append(new_name)

                adata_merged.obs['SEACell'] = new_seacell_names


        del adata_merged.obsm[reduction_subset]

    #    SEA2Cell_ad = SEACells.core.summarize_by_SEACell(adata_merged, SEACells_label='SEACell', celltype_label=cell_type)
        # SEA3Cell_ad.X contains the summarized data (sum)
        # Results are the same as create_seurat_adata.
    #    SEA2Cell_ad.raw = ad.AnnData(
    #        X=SEA2Cell_ad.X.copy(),
    #        var=SEA2Cell_ad.var.copy(),
    #        obs=SEA2Cell_ad.obs.copy()
    #    ) # Copy adata.X to adata.raw.X



        def create_seacell_adata(adata_merged, seacell_key='SEACell'):
            # Group by SEACell
            seacell_groups = adata_merged.obs.groupby(seacell_key)

            # Separate categorical and numeric columns
            cat_cols = adata_merged.obs.select_dtypes(include=['category', 'object']).columns
            num_cols = adata_merged.obs.select_dtypes(include=[np.number]).columns

            # Process categorical columns
            new_obs = pd.DataFrame()
            if len(cat_cols) > 0:
                cat_data = adata_merged.obs[cat_cols]
                new_obs = seacell_groups.apply(lambda x: cat_data.loc[x.index].mode().iloc[0])

            # Count number of cells
            new_obs['n_cells'] = seacell_groups.size()

            # Process raw.X data (raw count sum)
            if adata_merged.raw is not None:
                raw_X = np.zeros((len(new_obs), adata_merged.raw.n_vars))
                for i, (seacell, group) in enumerate(seacell_groups):
                    raw_X[i] = adata_merged.raw[group.index, :].X.sum(axis=0)
                raw_var = adata_merged.raw.var
            else:
                raw_X = None
                raw_var = None

            # Store raw count sum in .X (using the same gene set as raw.X)
            if raw_X is not None:
                new_X = raw_X.copy()
                use_var = raw_var
            else:
                # Fallback to .X sum when raw.X is not available
                new_X = np.zeros((len(new_obs), adata_merged.n_vars))
                for i, (seacell, group) in enumerate(seacell_groups):
                    new_X[i] = adata_merged[group.index, :].X.sum(axis=0)
                use_var = adata_merged.var

            # Calculate mean UMAP coordinates per SEACell
            umap_coords = pd.DataFrame(adata_merged.obsm['X_' + reduction_vis], index=adata_merged.obs.index)
            umap_coords[seacell_key] = adata_merged.obs[seacell_key]
            new_umap = umap_coords.groupby(seacell_key).mean().values
            new_obsm = {'X_' + reduction_vis: new_umap}

            # Create new AnnData object
            adata_seacell = sc.AnnData(X=new_X, obs=new_obs, var=use_var, obsm=new_obsm)

            # Save raw count sum in raw.X
            if raw_X is not None:
                adata_seacell.raw = sc.AnnData(X=raw_X, var=raw_var)

            return adata_seacell

        # Create new AnnData object using the function
        adata_seacell = create_seacell_adata(adata_merged)

        # CPM normalization
        if normalize_cpm:
            # Display raw count sum statistics
            if hasattr(adata_seacell.X, 'toarray'):
                total_counts = np.array(adata_seacell.X.sum(axis=1)).flatten()
            else:
                total_counts = np.array(adata_seacell.X.sum(axis=1)).flatten()

            st.markdown("#### CPM Normalization")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("Total counts (min)", f"{total_counts.min():.0f}")
            with col2:
                st.metric("Total counts (mean)", f"{total_counts.mean():.0f}")
            with col3:
                st.metric("Total counts (max)", f"{total_counts.max():.0f}")

            # Apply CPM normalization (no log transform)
            sc.pp.normalize_total(adata_seacell, target_sum=1e6)
            st.success("✅ CPM normalization applied (target_sum=1,000,000, no log transform)")
            st.caption("adata.X: CPM normalized / adata.raw.X: raw count sum")
        else:
            st.info("ℹ️ CPM normalization is disabled. adata.X contains raw count sums.")
            st.caption("⚠️ Downstream analyses (COMPASS, scFEA, etc.) require separate normalization.")

        fig, ax = plt.subplots(figsize=(4, 4))
        if cell_type is not None:
            sc.pl.scatter(adata_seacell, basis=reduction_vis, color=cell_type, frameon=True, ax=ax, show=False)
        else:
            sc.pl.scatter(adata_seacell, basis=reduction_vis, frameon=True, ax=ax, show=False)
        st.pyplot(fig)

        st.markdown("---")
#        st.write("adata_with_SEACells")
#        st.write(SEA2Cell_ad)
        st.write("SEACell.summarized")
        st.write(adata_seacell)

        # Generate group names
        adata_df = adata_seacell.obs
        group_dict = {}
        group_counter = {}

        adata_df['SEACell_group'] = adata_df.apply(generate_group_name, axis=1)
        st.write(adata_df['SEACell_group'])
        df = calc_df(adata_seacell, adata_df)

        file_name_head = os.path.splitext(uploaded_file.name)[0]
        save_dir_name = seacell_temp_dir + "/seacell/"

        if cell_type is not None:
            categories = adata_seacell.obs[cell_type].cat.categories.tolist()
            df_dict = {category: df.filter(regex=f'^{category}') for category in categories}

            # Display results
            for category, df_split in df_dict.items():
                st.write(f"\n{category}:")
                st.write(df_split.head())
                df_split.to_csv(os.path.join(save_dir_name, file_name_head + "_" + f"{category}_{cell_type}.SEAcells.tsv"), sep='\t', index=True)
        else:
            # Save all data as a single file when cell_type is not specified
            st.write("\nAll cells (single type):")
            st.write(df.head())
            df.to_csv(os.path.join(save_dir_name, file_name_head + "_all.SEAcells.tsv"), sep='\t', index=True)


        # Display results
        print(df)

        adata_seacell.write_h5ad(save_dir_name +file_name_head + "_SEAcells.summarized.h5ad", compression="gzip")

    #    SEA2Cell_ad.write_h5ad(save_dir_name + file_name_head + "_with_SEACells.h5ad", compression="gzip")

        if wgcna:
            df.to_csv(os.path.join(save_dir_name, file_name_head + "_SEACells.all.tsv"), sep='\t', index=True )

        shutil.make_archive(seacell_temp_dir + "/seacell", format='zip',root_dir= seacell_temp_dir + "/seacell/")

        if normalize_cpm:
            st.info("""
            **📝 About Metacell Data**

            - **adata.X**: raw count sum → **CPM normalized** (can be used directly for COMPASS, scFEA, etc.)
            - **adata.raw.X**: raw count sum (not normalized)
            - **TSV files**: raw count sum (based on raw.X)
            """)
        else:
            st.info("""
            **📝 About Metacell Data**

            The count for each metacell is the sum of raw counts from all cells within that metacell.

            ⚠️ **Note**: Normalization is required when using metacell data for downstream analyses.
            """)

        with open(seacell_temp_dir + "/seacell" + '.zip', "rb") as fp:
            btn = st.download_button(
                label="Download Results",
                data=fp,
                file_name=file_name_head + ".zip",
                mime="zip"
            )

