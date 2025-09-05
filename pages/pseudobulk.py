import scanpy as sc
import anndata as ad
import streamlit as st
import pandas as pd
import re
import os
from helper_func import clear_old_directories
from helper_func import clear_old_files
import matplotlib.pyplot as plt
import io
import seaborn as sns
import tempfile
import decoupler as dc
import shutil
import time
import scipy
import pickle
import numpy as np

#st.set_option('deprecation.showPyplotGlobalUse', False)

st.set_page_config(page_title="Pseudobulk by decoupleR", page_icon="💬")

# Add this new function to clean data
def clean_counts_layer(adata):
    """
    Clean NaN and inf values only from adata.layers['counts']
    """
    # Clean only layers['counts'] if it exists
    if 'counts' in adata.layers:
        if scipy.sparse.issparse(adata.layers['counts']):
            counts_data = adata.layers['counts'].data
            # Count NaN and inf values before cleaning
            nan_count = np.sum(np.isnan(counts_data))
            inf_count = np.sum(np.isinf(counts_data))
            
            # Replace NaN and inf with 0
            counts_data[np.isnan(counts_data)] = 0
            counts_data[np.isinf(counts_data)] = 0
            adata.layers['counts'].data = counts_data
            
            return nan_count, inf_count
        else:
            # Count NaN and inf values before cleaning
            nan_count = np.sum(np.isnan(adata.layers['counts']))
            inf_count = np.sum(np.isinf(adata.layers['counts']))
            
            adata.layers['counts'] = np.nan_to_num(adata.layers['counts'], nan=0.0, posinf=0.0, neginf=0.0)
            
            return nan_count, inf_count
    
    return 0, 0

@st.cache_data
def convert_df(df):
    return df.to_csv(index=True, sep='\t').encode('utf-8')


@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    if adata.raw.X is not None:
        adata.layers['counts'] = adata.raw.X # layersにコピー
    else:
        st.markdown("### No adata.raw.X data. Use adata.X?")
    return adata

@st.cache_data
def calc_pseudobulk1(_adata, groups_col, sample_col, min_cells, min_counts, mode="sum", cache_key=None, skip_checks=False):
    if mode == "mean_then_logit":
        # まず平均を計算
        pdata = dc.get_pseudobulk(
            _adata,
            groups_col=groups_col,
            sample_col=sample_col,
            layer='counts',
            mode='mean',
            min_cells=0,
            min_counts=0,
            skip_checks=skip_checks
        )
        # logit変換を適用
        # 0と1の値を調整
        eps = 1e-6
        X = pdata.X.copy()
        X = np.clip(X, eps, 1-eps)
        pdata.X = scipy.special.logit(X)
        return pdata
    else:
        return dc.get_pseudobulk(
            _adata,
            groups_col=groups_col,
            sample_col=sample_col,
            layer='counts',
            mode=mode,
            min_cells=0,
            min_counts=0,
            skip_checks=skip_checks
        )

@st.cache_data
def calc_pseudobulk2(_adata, groups_col, sample_col, min_cells, min_counts, mode="sum", cache_key=None, skip_checks=False):
    # まず通常のpseudolulk計算
    if mode == "mean_then_logit":
        pdata = dc.get_pseudobulk(
            _adata,
            groups_col=groups_col,
            sample_col=sample_col,
            layer='counts',
            mode='mean',  # まず平均を計算
            min_cells=min_cells,
            min_counts=min_counts,
            skip_checks=skip_checks
        )
        # logit変換を後から適用
        eps = 1e-6
        X = pdata.X.copy()
        X = np.clip(X, eps, 1-eps)
        pdata.X = scipy.special.logit(X)
        return pdata
    else:
        return dc.get_pseudobulk(
            _adata,
            groups_col=groups_col,
            sample_col=sample_col,
            layer='counts',
            mode=mode,
            min_cells=min_cells,
            min_counts=min_counts,
            skip_checks=skip_checks
        )

#@st.cache_data
def calc_df():
    df = pd.DataFrame(pdata.raw.X.T)
    df.columns = pdata.obs.index.to_list()
    df.index = pdata.var.index.to_list()
    df = df.sort_index(axis=1)
    return df

def download_h5ad():
    file_name = file_name_head + "_SEAcells.h5ad"
    # Create a temporary file
    with tempfile.NamedTemporaryFile(delete=False, suffix='.h5ad') as tmp:
        # Write the AnnData object to the temporary file

        pdata.write_h5ad(tmp.name, compression="gzip")

        # Read the temporary file into a BytesIO object
        with open(tmp.name, "rb") as f:
            buffer = io.BytesIO(f.read())

    # Remove the temporary file
    os.unlink(tmp.name)

def reshape_and_save_pseudobulk_data(pdata, sample_col, groups_col, save_dir):
    # Ensure the save directory exists
    os.makedirs(save_dir, exist_ok=True)
    
    # Extract data from pdata.X
    data = pdata.X.toarray() if scipy.sparse.issparse(pdata.X) else pdata.X
    
    # Get gene names
    gene_names = pdata.var_names
    
    # Create a dictionary to store DataFrames for each cell type (group)
    df_dict = {}
    
    # Get unique group values (cell types)
    unique_groups = pdata.obs[groups_col].unique()
    
    for group in unique_groups:
        # Clean group name
        cleaned_group = str(group).replace(' ', '_')
        # Filter data for the current group
        group_mask = pdata.obs[groups_col] == group
        group_data = data[group_mask]
        group_samples = pdata.obs.loc[group_mask, sample_col]
        
        # Create DataFrame for the current group
        df = pd.DataFrame(group_data, columns=gene_names, index=group_samples)
        
        # Group by sample_col and sum the values
        df = df.groupby(sample_col).sum().T

        # Clean column names by replacing spaces with underscores
        df.columns = [f"{cleaned_group}_{str(col).replace(' ', '_')}" for col in df.columns]
        
        # Store the DataFrame in the dictionary
        df_dict[group] = df
        
        # Save the DataFrame as TSV
        file_name = f"{group}_pseudobulk.tsv"
        file_path = os.path.join(save_dir, file_name)
        df.to_csv(file_path, sep='\t', index=True)
        
    return df_dict

def clean_column_names(adata: ad.AnnData) -> (ad.AnnData, dict, bool):
    """
    Clean column names in adata.obs to ensure they are valid, preserving dots and case.
    
    Parameters:
    -----------
    adata : anndata.AnnData
        The AnnData object whose obs column names need to be cleaned.
    
    Returns:
    --------
    anndata.AnnData
        The AnnData object with cleaned obs column names.
    dict
        A dictionary mapping original column names to cleaned column names (only for changed names).
    bool
        A flag indicating whether any changes were made to the column names.
    """
    def clean_name(name: str) -> str:
        # Remove any non-alphanumeric characters (except underscores and dots)
        # and replace spaces with underscores
        cleaned = re.sub(r'[^\w\s.]', '', name)
        cleaned = re.sub(r'\s+', '_', cleaned)
        
        # Ensure the name starts with a letter or underscore
        if not cleaned[0].isalpha() and cleaned[0] != '_':
            cleaned = '_' + cleaned
        
        return cleaned

    # Create a mapping of original to cleaned names
    name_mapping = {col: clean_name(col) for col in adata.obs.columns}
    
    # Check for duplicate cleaned names and modify if necessary
    seen_names = {}
    for original, cleaned in name_mapping.items():
        if cleaned in seen_names:
            counter = 1
            while f"{cleaned}_{counter}" in seen_names:
                counter += 1
            name_mapping[original] = f"{cleaned}_{counter}"
        seen_names[name_mapping[original]] = original

    # Keep only the changed mappings
    changed_mapping = {orig: cleaned for orig, cleaned in name_mapping.items() if orig != cleaned}

    # Check if any changes were made
    changes_made = bool(changed_mapping)

    # Rename the columns in adata.obs only if changes were made
    if changes_made:
        adata.obs.rename(columns=name_mapping, inplace=True)
    
    return adata, changed_mapping, changes_made

if "pseudobulk_temp_dir" not in st.session_state or type(st.session_state.pseudobulk_temp_dir) == bool:
    pseudobulk_temp_dir = os.path.join("temp", str(round(time.time())))
    # Create main temp directory if it doesn't exist
    os.makedirs("temp", exist_ok=True)
    # Clean up old directories
    clear_old_directories("temp")
    clear_old_files("temp")
    # Create new directories
    os.makedirs(pseudobulk_temp_dir, exist_ok=True)
    os.makedirs(os.path.join(pseudobulk_temp_dir, "df"), exist_ok=True)
    st.session_state.pseudobulk_temp_dir = pseudobulk_temp_dir
else:
    pseudobulk_temp_dir = st.session_state.pseudobulk_temp_dir


if "pdata" not in st.session_state:
    st.session_state.pdata = None

if "mete_edited" not in st.session_state:
    st.session_state.mete_edited = False


uploaded_file = st.file_uploader("Upload a h5ad file", type=['h5ad'])
use_X = st.checkbox('Use adata.X instead of adata.raw.X? [Counts are usually in adata.raw.X.]', value =False)
st.markdown("##### __通常はadata.raw.Xを用いsummarize methodはsum。__")
st.write("Normalized dataを用いて平均値を求める場合は、adata.Xを用い、summarize methodはmean。")


#---------------

if uploaded_file is not None:
    adata = read_h5ad(uploaded_file)
    # adata.raw.Xのデータをadata.Xへ
    if not use_X:
        adata.X = adata.raw.X
    else:  # use_X = Trueの場合
        adata.layers['counts'] = adata.X

    nan_count, inf_count = clean_counts_layer(adata)
    if nan_count > 0 or inf_count > 0:
        st.info(f"📊 Counts layerでNaN値 {nan_count}個、無限大値 {inf_count}個を0に置換しました。")

    st.markdown("__Uploaded data__")
    temp_df = pd.DataFrame(
    adata.X[:5,:8].toarray() if scipy.sparse.issparse(adata.X) else adata.X[:5,:8],
    index=adata.obs_names[:5],
    columns=adata.var_names[:8]
    )
    st.dataframe(temp_df) 
    st.write("Countsは基本的に整数 整数でない場合はデータ構造に注意")
 #   st.write(adata.X[:5,:5])
    st.markdown("__adata.obs (metadata)__")
    st.write(adata.obs.head())


    # Clean column names
    adata, name_mapping, changes_made = clean_column_names(adata)
    if changes_made:
        # Print mapping of original to cleaned names
        st.write("\nIdentity names have special characters. Mapping of original to cleaned names:")
        for original, cleaned in name_mapping.items():
            st.write(f"{original} -> {cleaned}")


    meta_data = adata.obs.columns.to_list()


    meta_data = list(filter(lambda x: x not in ['nCount_RNA', 'nFeature_RNA','percent.mt', 'Cell_id'], meta_data))

    meta_df = adata.obs[meta_data]

    one_cell = st.checkbox('Treat all data as one cell-type (Split only by samples)?', value =False)


    with st.form("Basic settings:"):
        sample_col = st.selectbox("Sample:", meta_data)

        if not one_cell:
            groups_col = st.selectbox("Cell type:", meta_data)

        
        edit_sample = st.checkbox("Edit sample names?", value=False)
        edit_cell = st.checkbox("Edit cell names?", value=False)

        assemble_all = st.checkbox("Generate a matrix file with all data?", value=False)
        save_h5ad = st.checkbox("Save h5ad file?", value = False)

        st.markdown("#### Quality threshold")
        min_cells = st.number_input('Minimum number of cells in a subgroup:', min_value =0,  value=10)

        min_counts= st.number_input('Minimum total counts in a subgroup:', min_value =0,  value=100)

        method = st.radio(
            "Summarzie method:",
            ["sum", "mean", "median", "mean_then_logit"],
            index=0,
        )
        st.write("mean_then_logit for SCENIC AUC")


        cache_key = time.time()

        submitted_basic = st.form_submit_button("Change the parameters")

    edit_meta = False
    if edit_sample or edit_cell:
        edit_meta = True
        with st.form("Edit metadata"):
            
            meta_df_sample = pd.DataFrame(meta_df[sample_col].cat.categories.tolist(), index = meta_df[sample_col].cat.categories.tolist(), columns = ['new sample name']).copy()
            meta_df_cell = pd.DataFrame(meta_df[groups_col].cat.categories.tolist(), index = meta_df[groups_col].cat.categories.tolist(), columns = ['new cell name']).copy()

            st.write('Editable:')

            if edit_sample:
                edited_df_sample = st.data_editor(meta_df_sample)
            if edit_cell:
                edited_df_cell = st.data_editor(meta_df_cell)

            finish_edit = st.form_submit_button("Finish")

            if finish_edit or st.session_state.mete_edited:
                if edit_sample:
                    category_mapping = edited_df_sample['new sample name'].to_dict()
                    adata.obs[sample_col] = adata.obs[sample_col].map(category_mapping)
                if edit_cell:
                    category_mapping = edited_df_cell['new cell name'].to_dict()
                    adata.obs[groups_col] = adata.obs[groups_col].map(category_mapping)

                st.write(adata.obs.head())
                st.session_state.mete_edited = True

    if (submitted_basic or st.session_state.pdata) and not (edit_meta and not st.session_state.mete_edited):

        if one_cell:
            adata.obs['One_cell'] = 'one_cell_type'
            groups_col = 'One_cell'
            st.write("One_cell is added.")
            cache_key = time.time()
        # Get pseudo-bulk profile for all
        skip_checks_flag = False

        try:
            pdata =  calc_pseudobulk1(adata, groups_col, sample_col, min_cells, min_counts, mode=method, cache_key=cache_key)
        except ValueError:
            st.warning("countsが少数を含みます。データが正しいかどうか確認を。少数のまま計算します。")
            pdata = calc_pseudobulk1(adata, groups_col, sample_col, min_cells, min_counts, mode=method, cache_key=cache_key, skip_checks=True)
            skip_checks_flag = True

        if one_cell:
            pdata.obs['One_cell'] = 'one_cell_type' #pdataにも追加する
        st.write(pdata.obs.head(3))

        st.markdown("#### All data without QC")

        dc.plot_psbulk_samples(pdata, groupby=[sample_col, groups_col], figsize=(12, 4))
        fig = plt.gcf() 
        st.pyplot(fig)  # 明示的に図をst.pyplot()に渡す
        plt.close(fig)  # メモリを解放するために図を閉じる



        # Get pseudo-bulk profile
        pdata =  calc_pseudobulk2(adata, groups_col, sample_col, min_cells, min_counts, mode=method, cache_key=cache_key,  skip_checks=skip_checks_flag)
        if one_cell:
            pdata.obs['One_cell'] = 'one_cell_type' #pdataにも追加する

        st.markdown('#### After QC')

        dc.plot_psbulk_samples(pdata, groupby=[sample_col, groups_col], figsize=(12, 4))
        fig = plt.gcf() 
        st.pyplot(fig)
        plt.close(fig) 

        pdata.raw = ad.AnnData(
                X=pdata.X.copy(),
                var=pdata.var.copy(),
                obs=pdata.obs.copy()
            )


        st.session_state.pdata = True

        df = calc_df()


        st.write(df.head())


        st.markdown("---")
        file_name_head = os.path.splitext(uploaded_file.name)[0]

        file_name_head = file_name_head + sample_col + '_' + groups_col + ".pseudobulk"

        if st.button('Prepare files to download'):


            reshaped_df_dict = reshape_and_save_pseudobulk_data(pdata, sample_col, groups_col, pseudobulk_temp_dir + "/df")

            
            st.success("Data reshaping complete!")
            
            # Display preview for each sample
            for sample, df in reshaped_df_dict.items():
                st.subheader(f"Sample: {sample}")
                st.write("Shape:", df.shape)
                st.write("Preview:")
                st.dataframe(df.head())
                st.write("---")


            if save_h5ad:
                pdata.write_h5ad(pseudobulk_temp_dir + "/df/" + file_name_head + ".h5ad", compression="gzip")


       #     with open("temp/pseudo.pkl", "wb") as tf:
       #         pickle.dump(reshaped_df_dict, tf)

            if assemble_all:
                dict_keys = list(reshaped_df_dict.keys())
                df_all = reshaped_df_dict[dict_keys[0]]
                # 元々のカラム名を保持
                df_all.columns = [col for col in df_all.columns]
                
                for x in range(len(dict_keys)-1):
                    df_temp = reshaped_df_dict[dict_keys[x+1]]
                    # 元々のカラム名を保持
                    df_temp.columns = [col for col in df_temp.columns]
                    df_all = pd.concat([df_all, df_temp], axis=1)
            
                df_all.to_csv(pseudobulk_temp_dir + "/df/" + file_name_head + ".all.tsv", sep='\t', index=True)


            shutil.make_archive(pseudobulk_temp_dir + "/pseudobulk", format='zip',root_dir= pseudobulk_temp_dir + "/df")

            with open(pseudobulk_temp_dir + "/pseudobulk" + '.zip', "rb") as fp:
                btn = st.download_button(
                    label="Download Results",
                data=fp,
                file_name=file_name_head + ".zip",
                mime = "zip"
                )