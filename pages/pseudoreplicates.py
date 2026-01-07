import scanpy as sc
import anndata as ad
import numpy as np
import streamlit as st
import pandas as pd
import os
import matplotlib.pyplot as plt
import io
import tempfile
from sklearn.model_selection import KFold
import shutil
from anndata import Raw
import time
from helper_func import clear_old_directories
from helper_func import clear_old_files

# st.set_option('deprecation.showPyplotGlobalUse', False)

st.set_page_config(page_title="Random pseudo-replicates", page_icon="💬")

st.markdown("### Generate random pseudo-replicates")

@st.cache_data
def read_h5ad(file):
	adata = sc.read_h5ad(file)
	return adata


def process_adata(adata, max_splits=3, sample_col=None, cell_type=None, dim_use='X_umap'):
    # "sample"and"cell_type"ofsetmimatchwasewithGroupize
    if sample_col is not None:
        adata.obs['group'] = (
            adata.obs[cell_type].astype(str) + '_' +
            adata.obs[sample_col].astype(str)
        )
    else:
        adata.obs['group'] = adata.obs[cell_type].astype(str)

    new_obs = []
    new_X = []
    new_raw_X = []
    new_umap = []

    for group in adata.obs['group'].unique():
        group_adata = adata[adata.obs['group'] == group]
        n_samples = len(group_adata)

        # GroupofSamplenumtobaseduiten_splitstheadjustarrange
        n_splits = min(max_splits, n_samples)

        if n_splits == 1:
            # Sampleis1tsushikanotplacematch、divratiosezutosoofmamauseuse
            split_adata = group_adata
            count_sum = np.sum(split_adata.raw.X.toarray(), axis=0)
            umap_mean = np.mean(split_adata.obsm[dim_use], axis=0)
            new_obs_entry = {
                cell_type: split_adata.obs[cell_type].iloc[0],
                'group': group,
                'split': 0,
                'n_cells': len(split_adata)
            }
            if sample_col is not None:
                new_obs_entry[sample_col] = split_adata.obs[sample_col].iloc[0]
            new_obs.append(new_obs_entry)
            new_X.append(count_sum)
        #    new_raw_X.append(count_sum)  # sourceofkauntoData
            new_umap.append(umap_mean)
        else:
            # ndivratio
            kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
            for i, (_, split_idx) in enumerate(kf.split(group_adata.obs_names)):
                split_adata = group_adata[group_adata.obs_names[split_idx]]
                # kauntoDataofsumtheCalculation
                count_sum = np.sum(split_adata.raw.X.toarray(), axis=0)
                # UMAPofflatavgvaltheCalculation
                umap_mean = np.mean(split_adata.obsm[dim_use], axis=0)
                # newshiiDatatheaddadd
                new_obs_entry = {
                    cell_type: split_adata.obs[cell_type].iloc[0],
                    'group': group,
                    'split': i,
                    'n_cells': len(split_adata)
                }
                if sample_col is not None:
                    new_obs_entry[sample_col] = split_adata.obs[sample_col].iloc[0]
                new_obs.append(new_obs_entry)
                new_X.append(count_sum)
           #     new_raw_X.append(count_sum)  # sourceofkauntoData
                new_umap.append(umap_mean)

    # newshiiAnnDataobujiekutothemakebecome
    new_obs = pd.DataFrame(new_obs)
    new_X = np.vstack(new_X)
#    new_X = np.array(new_X) #korewithisumakuikanot
 #   st.write(new_X)
 #   new_raw_X = np.array(new_raw_X)
    new_umap = np.array(new_umap)
    new_adata = ad.AnnData(X=new_X, obs=new_obs, var=adata.var)
    new_adata.raw = ad.AnnData(X=new_X, var=adata.var)
    new_adata.obsm[dim_use] = new_umap

    return new_adata


def calc_df(adata):
	df = pd.DataFrame(adata.raw.X.T)
	combined_list = [f"{group}_{split}" for group, split in zip(adata.obs['group'], adata.obs['split'])]
	df.columns = combined_list
	df.index = adata.var.index.to_list()
	df = df.sort_index(axis=1)
	return df

uploaded_file = st.file_uploader("Upload a h5ad file", type=['h5ad'])

if "pseudobulk_temp_dir" not in st.session_state or type(st.session_state.pseudobulk_temp_dir) == bool:
	pseudobulk_temp_dir = "temp/" + str(round(time.time()))
	if not os.path.exists('temp'):
		os.mkdir('temp')
	else:
		clear_old_directories("temp")
		clear_old_files("temp")
	os.mkdir(pseudobulk_temp_dir)
	os.mkdir(pseudobulk_temp_dir + "/df")
	st.session_state.pseudobulk_temp_dir = pseudobulk_temp_dir
else:
	pseudobulk_temp_dir = st.session_state.pseudobulk_temp_dir

#---------------

if uploaded_file is not None:
	adata = read_h5ad(uploaded_file)
	st.write("Uploaded data")
	st.write(adata)
	st.write(adata.X[:5,:5])



	meta_data = adata.obs.columns.to_list()


	meta_data = list(filter(lambda x: x not in ['nCount_RNA', 'nFeature_RNA','percent.mt', 'Cell_id'], meta_data))
	obsm_keys = [x.replace('X_','') for x in list(adata.obsm)]
	umap_keys = list(filter(lambda x: 'umap' in x.lower(), obsm_keys))
	umap_key = umap_keys[0] if umap_keys else None
	if umap_key is not None:
		umap_index = umap_keys.index(umap_key)
	else:
		umap_index = 0

	reduction_vis = umap_key
	reduction_use = obsm_keys[0]
	reduction_subset = 'X_' + reduction_use + '_sub'


	split_sel = st.checkbox('Split by samples?', value= False)

	with st.form("Basic settings:"):
		if split_sel:
			sample_col = st.selectbox("Split by:", meta_data)
		else:
			sample_col = None

		cell_type = st.selectbox("Identity of cell types:", meta_data)

		max_splits = st.number_input('How many divisions in each cell-type?', value=3)

		obsm_keys = [x.replace('X_','') for x in list(adata.obsm)]

		reduction_vis = st.selectbox("Reduciton to show (e.g. umap):",umap_keys, index = umap_index)

		submitted_basic = st.form_submit_button("Change the parameters and go")




	if submitted_basic:
		fig, ax = plt.subplots()

		fig, ax = plt.subplots(figsize=(4, 4))
		sc.pl.scatter(adata, basis=reduction_vis, color=cell_type, frameon=True, ax = ax, show = False)
		st.pyplot(fig)

		adata_sub = process_adata(adata,
							  max_splits=max_splits,
							  sample_col=sample_col,
							  cell_type=cell_type,
							  dim_use=f'X_{reduction_vis}')
		adata_sub.raw = adata_sub.copy()

		sc.pp.normalize_per_cell(adata_sub, counts_per_cell_after=1e4)
		sc.pp.log1p(adata_sub)

		st.markdown("#### Pseudoreplicates")

		fig, ax = plt.subplots(figsize=(4, 4))
		sc.pl.scatter(adata_sub, basis=reduction_vis, color=cell_type,  size=50, ax = ax, show = False)
		st.pyplot(fig)
		st.write(adata_sub)
		file_name_head = os.path.splitext(uploaded_file.name)[0]
		adata_sub.write_h5ad(pseudobulk_temp_dir + "/df/" + file_name_head + ".pseudoreplicates.h5ad", compression="gzip")

		df = calc_df(adata_sub)
		st.write(df.head())

		categories = adata_sub.obs[cell_type].cat.categories.tolist()
		# eachkategori-tobaseduiteDatafure-muthedivratio

		df_dict = {category: df.filter(regex=f'^{category}') for category in categories}



		# ResulttheDisplay
		for category, df_split in df_dict.items():
			st.write(f"\n{category}:")
			st.write(df_split.head())
			df_split.to_csv(os.path.join(pseudobulk_temp_dir, "df",file_name_head + "_" + f"{category}_{cell_type}.pseudoreplicates.tsv"), sep='\t', index=True)

		st.markdown("---")


		shutil.make_archive(pseudobulk_temp_dir + "/pseudoreplicates", format='zip',root_dir= pseudobulk_temp_dir + "/df")

		with open(pseudobulk_temp_dir + "/pseudoreplicates" + '.zip', "rb") as fp:
			btn = st.download_button(
				label="Download Results",
			data=fp,
			file_name=file_name_head + ".zip",
			mime = "zip"
			)
