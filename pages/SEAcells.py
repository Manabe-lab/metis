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
st.write("各metacellに含まれる細胞のraw countの合計がcountとなる。")
st.write("metacellのデータを用いた解析ではcountデータとして扱い、normalizationが必要。")

@st.cache_data
def read_h5ad(file):
    adata = sc.read_h5ad(file)
    return adata

def generate_group_name(row):
    if sample_id == 'single_sample':
        key = row[cell_type]
    else:
        key = (row[sample_id], row[cell_type])
    if key not in group_counter:
        group_counter[key] = 0
    group_counter[key] += 1
    if sample_id == 'single_sample':
        return f"{row[cell_type]}_{group_counter[key]}"
    else:
        return f"{row[cell_type]}_{row[sample_id]}_{group_counter[key]}"

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
    st.write("Snippet of the count matrix, adata.raw.X")
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


    split_sel = st.checkbox('Split by samples?', value= False)
    wgcna = st.checkbox("Generate single matrix file for downsteram analysis (e.g., WGCNA)?", value = True)

    sample_id = "single_sample"

    with st.form("Basic settings:"):
        if split_sel:
            sample_id = st.selectbox("Split by:", meta_data)
        cell_type = st.selectbox("Identity of cell types:", meta_data)

        reduction_use = st.selectbox("Reduciton to calculate クラスタリングに用いた次元圧縮 (eg. PCA, mnn)", obsm_keys)
        reduction_subset = 'X_' + reduction_use + '_sub'

        n_dim = st.number_input('Number of dimension used to identify clusters クラスタリングに用いた次元数', min_value =1, value=30)

        obsm_keys = [x.replace('X_','') for x in list(adata.obsm)]

        reduction_vis = st.selectbox("Reducition to visualize (e.g. umap):",umap_keys, index = umap_index)

        fracion_n_SEACs = st.number_input('n_SEACells: approximate the number of cells to generate', min_value =0, value=round(len(adata)/75))
        st.write("As a rule of thumb, choosing one metacell for every 75 single-cells.")
        st.write("サンプル全体での数.")
        n_waypoint_eigs = st.number_input('n_waypoint_eigs:',  value=10)

        submitted_basic = st.form_submit_button("Set the parameters and run")


    if submitted_basic:
        fig, ax = plt.subplots()

        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.scatter(adata, basis=reduction_vis, color=cell_type, frameon=True, ax = ax, show = False)
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

            # total sample数でfracion_n_SEACsは計算しなおす

            ## 75 cell:bin
#            n_cells = round(len(adata_dict[i].obs) / fracion_n_SEACs)
            n_SEACells = round(len(adata_dict[i]) / total_cell_num * fracion_n_SEACs)
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
                          use_gpu = True            )


            model.construct_kernel_matrix()
            M = model.kernel_matrix
            st.write(f"M shape: {M.shape}")

            # Initialize archetypes
            try:
                model.initialize_archetypes()
            except Exception as e:
                st.write(f"Error: {e}")
                st.write('多分このサンプルのnSEAcellsが少なすぎます。このサンプルを除いて続行します。')
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
                st.write(f"Error: {e}")
                st.write('Please add "import cupyx.scipy.sparse" and "import cupy as cp" to core.py of SEACells.')


            # Check for convergence
            model.plot_convergence(save_as=buf,  show=False)
            buf.seek(0)
            st.image(buf)
            plt.clf()
            plt.close('all')


            model.fit(min_iter=10, max_iter=50)


            adata_dict[i].obs[['SEACell']].head()

            SEACell_adata[i] = SEACells.core.summarize_by_SEACell(adata_dict[i], SEACells_label='SEACell', summarize_layer='raw')

            SEACells.plot.plot_2D(adata_dict[i], key='X_'+reduction_vis, colour_metacells=True)

            SEACells.plot.plot_SEACell_sizes(adata_dict[i], bins=5, save_as=buf,  show=False)
            buf.seek(0)
            st.image(buf)
            plt.clf()
            plt.close('all')


            SEACell_purity = SEACells.evaluate.compute_celltype_purity(adata_dict[i], cell_type)

            # Create a new figure with 3 subplots in a row
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12, 4))

            # Plot 1: Celltype Purity
            sns.boxplot(data=SEACell_purity, y=cell_type + '_purity', ax=ax1)
            ax1.set_title('Celltype Purity')
            sns.despine(ax=ax1)

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
        # AnnData オブジェクトをマージ
            adata_merged = ad.concat(
                adata_dict,
                join='outer',
                merge='same',
                label='sample'
            )
            # SEAcellの名前を変更
            for i in range(len(adata_merged.obs['SEACell'])):
                if adata_merged.obs[sample_id].iloc[i] not in adata_merged.obs['SEACell'].iloc[i]:
                    adata_merged.obs['SEACell'].iloc[i] = adata_merged.obs[sample_id].iloc[i] + '_' + adata_merged.obs['SEACell'].iloc[i]

        else:
            adata_merged = adata_dict["all"]


        del adata_merged.obsm[reduction_subset]

    #    SEA2Cell_ad = SEACells.core.summarize_by_SEACell(adata_merged, SEACells_label='SEACell', celltype_label=cell_type)
        # SEA3Cell_ad.Xにsummarizeしたデータが入る sum
        # create_seurat_adataと結果は同じになる。
    #    SEA2Cell_ad.raw = ad.AnnData(
    #        X=SEA2Cell_ad.X.copy(),
    #        var=SEA2Cell_ad.var.copy(),
    #        obs=SEA2Cell_ad.obs.copy()
    #    ) # adata.Xをadta.raw.Xにコピー



        def create_seacell_adata(adata_merged, seacell_key='SEACell'):
            # SEACellごとにグループ化
            seacell_groups = adata_merged.obs.groupby(seacell_key)
            
            # カテゴリカル列と数値列を分離
            cat_cols = adata.obs.select_dtypes(include=['category', 'object']).columns
            num_cols = adata.obs.select_dtypes(include=[np.number]).columns

            # カテゴリカル列の処理
            new_obs = pd.DataFrame()
            if len(cat_cols) > 0:
                cat_data = adata.obs[cat_cols]
                new_obs = seacell_groups.apply(lambda x: cat_data.loc[x.index].mode().iloc[0])

            # 細胞数のカウント
            new_obs['n_cells'] = seacell_groups.size()

            # 新しいXマトリックスを作成（細胞カウントの合計）
            new_X = np.zeros((len(new_obs), adata_merged.n_vars))
            for i, (seacell, group) in enumerate(seacell_groups):
                new_X[i] = adata_merged[group.index, :].X.sum(axis=0)

            # SEACellごとのUMAP座標の平均を計算
            umap_coords = pd.DataFrame(adata_merged.obsm['X_' + reduction_vis], index=adata_merged.obs.index)
            umap_coords[seacell_key] = adata_merged.obs[seacell_key]
            new_umap = umap_coords.groupby(seacell_key).mean().values
            new_obsm = {'X_' + reduction_vis: new_umap}

            # 新しいAnnDataオブジェクトを作成
            adata_seacell = sc.AnnData(X=new_X, obs=new_obs, var=adata_merged.var, obsm=new_obsm)

            # raw.Xのデータを処理
            if adata_merged.raw is not None:
                raw_X = np.zeros((len(new_obs), adata_merged.raw.n_vars))
                for i, (seacell, group) in enumerate(seacell_groups):
                    raw_X[i] = adata_merged.raw[group.index, :].X.sum(axis=0)
                adata_seacell.raw = sc.AnnData(X=raw_X, var=adata_merged.raw.var)
            return adata_seacell

        # 関数を使用して新しいAnnDataオブジェクトを作成
        adata_seacell = create_seacell_adata(adata_merged)


        fig, ax = plt.subplots(figsize=(4, 4))
        sc.pl.scatter(adata_seacell, basis=reduction_vis, color=cell_type, frameon=True, ax = ax, show = False)
        st.pyplot(fig)

        st.markdown("---")
#        st.write("adata_with_SEACells")
#        st.write(SEA2Cell_ad)
        st.write("SEACell.summarized")
        st.write(adata_seacell)

        # グループ名の生成
        adata_df = adata_seacell.obs
        group_dict = {}
        group_counter = {}

        adata_df['SEACell_group'] = adata_df.apply(generate_group_name, axis=1)
        st.write(adata_df['SEACell_group'])
        df = calc_df(adata_seacell, adata_df)
        categories = adata_seacell.obs[cell_type].cat.categories.tolist()

        df_dict = {category: df.filter(regex=f'^{category}') for category in categories}

        file_name_head = os.path.splitext(uploaded_file.name)[0]

        save_dir_name = seacell_temp_dir + "/seacell/"


        # 結果を表示
        for category, df_split in df_dict.items():
            st.write(f"\n{category}:")
            st.write(df_split.head())
            df_split.to_csv(os.path.join(save_dir_name, file_name_head + "_" + f"{category}_{cell_type}.SEAcells.tsv"), sep='\t', index=True)


        # 結果の表示
        print(df)

        adata_seacell.write_h5ad(save_dir_name +file_name_head + "_SEAcells.summarized.h5ad", compression="gzip")

    #    SEA2Cell_ad.write_h5ad(save_dir_name + file_name_head + "_with_SEACells.h5ad", compression="gzip")

        if wgcna:
            df.to_csv(os.path.join(save_dir_name, file_name_head + "_SEACells.all.tsv"), sep='\t', index=True )

        shutil.make_archive(seacell_temp_dir + "/seacell", format='zip',root_dir= seacell_temp_dir + "/seacell/")


        with open(seacell_temp_dir + "/seacell" + '.zip', "rb") as fp:
            btn = st.download_button(
                label="Download Results",
            data=fp,
            file_name=file_name_head + ".zip",
            mime = "zip"
            )

