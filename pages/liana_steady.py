import liana as li
# needed for visualization and toy data
import scanpy as sc

from liana.method import singlecellsignalr, connectome, cellphonedb, natmi, logfc, cellchat, geometric_mean, scseqcomm
from liana.mt import rank_aggregate
from functools import lru_cache
import itertools
from plotnine import scale_y_discrete

import numpy as np
import streamlit as st
import pandas as pd
import numpy as np
import os
from helper_func import clear_old_directories
from helper_func import clear_old_files
from helper_func import check_species_index
import time
import matplotlib.pyplot as plt

import datetime

# we import plotnine
import plotnine as p9
from plotnine import scale_color_cmap

st.set_page_config(page_title="Liana_steady-state", page_icon="💬")


st.sidebar.title("Options")

@st.cache_resource
def load_geneinfo():
	geneinfo_human = pd.read_csv("db/nichenetr.db/geneinfo_human.tsv", sep='\t')
	geneinfo_2022 = pd.read_csv("db/nichenetr.db/geneinfo_2022.tsv", sep='\t')
	return geneinfo_human, geneinfo_2022

@st.cache_data
def read_h5ad(file):
	adata = sc.read_h5ad(file)
	return adata

@st.cache_data
def read_map_df():
	map_df = li.rs.get_hcop_orthologs(url='https://ftp.ebi.ac.uk/pub/databases/genenames/hcop/human_mouse_hcop_fifteen_column.txt.gz',
			  columns=['human_symbol', 'mouse_symbol'],
			   # NOTE: HCOP integrates multiple resource, so we can filter out mappings in at least 3 of them for confidence
			   min_evidence=3
			   )
	# rename the columns to source and target, respectively for the original organism and the target organism
	map_df = map_df.rename(columns={'human_symbol':'source', 'mouse_symbol':'target'})
	return map_df

#@st.cache_data
#def mouse_conversion(resource, map_df):
#	mouse = li.rs.translate_resource(resource,
#								 map_df=map_df,
#								 columns=['ligand', 'receptor'],
#								 replace=True,
#								 # Here, we will be harsher and only keep mappings that don't map to more than 1 mouse gene
#								 one_to_many=1
#								 )
#	return mouse


def find_first_index_or_zero(lst, elements):
	for element in elements:
		try:
			return lst.index(element)
		except ValueError:
			continue
	return 0

def has_nan_like(lst):
	nan_values = {'nan', 'NaN', 'NAN', 'n/a', 'N/A', '', ' '}
	return any(str(item).strip().lower() in nan_values for item in lst)

@st.cache_data
def convert_human_to_mouse_symbols(symbols, version=1): # nichenetrのRスクリプトをClaude3.5で変換
	if not isinstance(symbols, (list, pd.Series)):
		raise ValueError("symbols should be a list or pandas Series of human gene symbols")
	if version == 1:
		geneinfo = geneinfo_human
	elif version == 2:
		geneinfo = geneinfo_2022
	else:
		raise ValueError("version must be 1 or 2")
	unambiguous_mouse_genes = (
		geneinfo.dropna()
		.groupby('symbol_mouse').size()
		.reset_index(name='count')
		.query('count < 2')['symbol_mouse']
		.tolist()
	)
	ambiguous_mouse_genes = (
		geneinfo.dropna()
		.groupby('symbol_mouse').size()
		.reset_index(name='count')
		.query('count >= 2')['symbol_mouse']
		.tolist()
	)
	geneinfo_ambiguous_solved = geneinfo[
		(geneinfo['symbol_mouse'].isin(ambiguous_mouse_genes)) &
		(geneinfo['symbol'] == geneinfo['symbol_mouse'].str.upper())
	]
	geneinfo = pd.concat([
		geneinfo[geneinfo['symbol_mouse'].isin(unambiguous_mouse_genes)],
		geneinfo_ambiguous_solved
	]).dropna()
	humansymbol2mousesymbol = dict(zip(geneinfo['symbol'], geneinfo['symbol_mouse']))
	converted_symbols = [humansymbol2mousesymbol.get(symbol, np.nan) for symbol in symbols]
	return converted_symbols

#@st.cache_data
#def mouse_conversion(resource):
#	st.write("Converting to mouse genes...This takes a while...")
#	nonconversion = pd.DataFrame(columns=["ligand", "receptor", "mouse_ligand","mouse_receptor"])
#	mouse = pd.DataFrame(columns=["ligand", "receptor"])
#	to_drop = []
#	for i in range(len(resource)):
#		row = resource.iloc[i,:].to_list()
#		new_row = convert_human_to_mouse_symbols(row)
#		if has_nan_like(new_row):
#			row.extend(new_row)
#			nonconversion.loc[len(nonconversion)] = row
#			to_drop.append(i)
#		else:
#			mouse.loc[len(mouse)] = new_row
#	mouse = mouse.drop_duplicates()
#	return mouse, nonconversion

@lru_cache(maxsize=None)
def cached_convert_human_to_mouse_symbols(symbol):
	result = convert_human_to_mouse_symbols([symbol])
	return result[0] if result else np.nan


@st.cache_data
def mouse_conversion(resource):
 #   st.write("Converting to mouse genes... This may take a moment...")

	def convert_symbol(symbol):
		if pd.isna(symbol):
			return np.nan
		if '_' in symbol:
			parts = symbol.split('_')
			converted = [cached_convert_human_to_mouse_symbols(part) for part in parts]
			return '_'.join(converted) if all(not pd.isna(part) for part in converted) else np.nan
		return cached_convert_human_to_mouse_symbols(symbol)

	# Vectorize the convert_symbol function
	vectorized_convert = np.vectorize(convert_symbol, otypes=[object])

	# Convert all symbols at once
	mouse_ligands = vectorized_convert(resource['ligand'])
	mouse_receptors = vectorized_convert(resource['receptor'])

	# Create a DataFrame with both human and mouse symbols
	result = pd.DataFrame({
		'ligand': resource['ligand'],
		'receptor': resource['receptor'],
		'mouse_ligand': mouse_ligands,
		'mouse_receptor': mouse_receptors
	})

	# Apply has_nan_like to each row
	has_nan = result[['mouse_ligand', 'mouse_receptor']].isna().any(axis=1)

	# Split into converted and non-converted
	nonconversion = result[has_nan]
	mouse = result[~has_nan][['mouse_ligand', 'mouse_receptor']].rename(
		columns={'mouse_ligand': 'ligand', 'mouse_receptor': 'receptor'}
	)

	# Drop duplicates from mouse DataFrame
	mouse = mouse.drop_duplicates()

	return mouse, nonconversion
#
#def mouse_conversion(resource):
#	st.write("Converting to mouse genes... This takes a while...")
#	# Convert all symbols at once
#	ligands = resource['ligand'].tolist()
#	receptors = resource['receptor'].tolist()
#	mouse_ligands = convert_human_to_mouse_symbols(ligands)
#	mouse_receptors = convert_human_to_mouse_symbols(receptors)
#	# Create a DataFrame with both human and mouse symbols
#	result = pd.DataFrame({
#		'ligand': ligands,
#		'receptor': receptors,
#		'mouse_ligand': mouse_ligands,
#		'mouse_receptor': mouse_receptors
#	})
#	# Apply has_nan_like to each row
#	has_nan = result.apply(lambda row: has_nan_like(row['mouse_ligand':'mouse_receptor']), axis=1)
#	# Split into converted and non-converted
#	nonconversion = result[has_nan]
#	mouse = result[~has_nan][['mouse_ligand', 'mouse_receptor']].rename(
#		columns={'mouse_ligand': 'ligand', 'mouse_receptor': 'receptor'}
#	)
#	# Drop duplicates from mouse DataFrame
#	mouse = mouse.drop_duplicates()
#	return mouse, nonconversion

@st.cache_data
def gnerate_interaction_matrix(df, LR_count_threshold):
	def calculate_interaction_metrics(df, source_cell, target_cell, LR_count_threshold):
		# Filter the data for the specific source and target cell types
			filtered_df = df[(df['source'] == source_cell) & (df['target'] == target_cell)]
			# Calculate the total interaction strength
			interaction_strength = filtered_df[magnitude_col].sum()
			# Count the number of unique ligand-receptor pairs
			lr_count = len(filtered_df[filtered_df[specificity_col]<LR_count_threshold])
			return interaction_strength, lr_count

	# Get unique cell types
	cell_types = set(df['source'].unique()) | set(df['target'].unique())

	# Calculate interaction metrics for all combinations
	results = []
	for source, target in itertools.product(cell_types, repeat=2):
		strength, count = calculate_interaction_metrics(df, source, target, LR_count_threshold)
		results.append({
			'source': source,
			'target': target,
			'interaction_strength': strength,
			'lr_pair_count': count
			})

	# Convert results to a DataFrame
	results_df = pd.DataFrame(results)

	# Sort results by interaction strength in descending order
	results_df = results_df.sort_values('interaction_strength', ascending=False)
	return results_df

@st.cache_data
def calc_lr(_adata,rankaggregate_options, func_name):
	function_dict[func_name](adata,**rankaggregate_options)
	return adata


if "liana_res" not in st.session_state:
	st.session_state.liana_res = []

#if "liana_basic_change" not in st.session_state:
#	st.session_state.liana_basic_change = False

if "liana_temp_dir" not in st.session_state:
	st.session_state.liana_temp_dir = True
	liana_temp_dir = "temp/" + str(round(time.time()))
	if not os.path.exists('temp'):
		os.mkdir('temp')
	else:
		clear_old_directories("temp")
		clear_old_files("temp")
	os.mkdir(liana_temp_dir)
	st.session_state.liana_temp_dir = liana_temp_dir
else:
	liana_temp_dir = st.session_state.liana_temp_dir

st.markdown("### L-R interaction calculation using LIANA+ in one condition")

function_dict = {
'singlecellsignalr': singlecellsignalr,
'connectome': connectome,
'cellphonedb': cellphonedb,
'natmi': natmi,
'logfc': logfc,
'cellchat': cellchat,
'geometric_mean': geometric_mean,
'rank_aggregate': rank_aggregate,
'scseqcomm': scseqcomm
}

method_list = ['CellPhoneDB','Connectome','log2FC','NATMI','SingleCellSignalR','Rank_Aggregate','GeometricMean','scSeqComm','CellChat']
func_list = ['cellphonedb','connectome','logfc','natmi','singlecellsignalr','rank_aggregate','geometric_mean','scseqcomm','cellchat']
lr_means_list = ['lr_means','expr_prod','None','expr_prod','lrscore','magnitude_rank','lr_gmeans','inter_score','lr_probs']
pval_col = ['cellphone_pvals','scaled_weight','lr_logfc','spec_weight','None','specificity_rank','gmean_pvals','None','cellchat_pvals']


#============ 新しいファイルをアップロードしたときは、cacheをclearする

def get_file_identifier(file):
	if file is not None:
		return f"{file.name}_{file.size}"
	return None


uploaded_file = st.file_uploader("Upload a h5ad file", type=['h5ad'])

if uploaded_file  is not None:
	current_file_id = get_file_identifier(uploaded_file)

	if 'last_file_id' not in st.session_state:
		st.session_state.last_file_id = None

	if current_file_id != st.session_state.last_file_id:
		st.cache_data.clear()
		st.cache_resource.clear()
		st.session_state.last_file_id = current_file_id
		st.success("新しいファイルが検出されました。キャッシュをクリアしました。")

#---------------

	adata = read_h5ad(uploaded_file)
	st.write("Uploaded data")
	st.write(adata.X[:3,:3])

	meta = adata.obs.columns.to_list()

	for i in ['nFeature_RNA','nCount_RNA','percent.mt', 'Cell_id']:
		try:
			meta.remove(i)
		except:
			pass
	submitted_basic = False


	groupby = st.selectbox("cell identity:", meta, index = find_first_index_or_zero(meta, ["cell.ident",
	"seurat_clusters"]))
	cell_list = sorted(adata.obs[groupby].cat.categories.to_list())
	st.write(", ".join(cell_list))
	with st.form("Basic settings:"):
		species = st.radio("Species:", ('mouse','human'), index = check_species_index(adata.var.index.to_list()[:50]))
		method = st.selectbox(
		"Method:",
		('Rank_Aggregate','CellPhoneDB','CellChat','Connectome','log2FC','NATMI','SingleCellSignalR','Geometric Mean','scSeqComm'), index = 0)

		func_name = func_list[method_list.index(method)]

		if species == "mouse":
			db_default = 1
		else:
			db_default = 0
		db = st.selectbox("database:", ('consensus', 'mouseconsensus', 'cellphonedb','cellchatdb', 'cellchat_secreted_signaling', 'nichenet','baccin2019', 'cellcall', 'cellinker', 'celltalkdb', 'connectomedb2020',  'embrace', 'guide2pharma', 'hpmr', 'icellnet', 'italk', 'kirouac2010', 'lrdb','ramilowski2015'), index = db_default)
		st.write("LR DBs are in human symbols except mouseconsensus. For mouse, they will be converted.")

		expr_prop = st.number_input('Threshold for min fraction of cells expressing the gene', min_value =0.0, max_value=0.9,
		step =0.01, value=0.1)
		remove_complex = st.checkbox("Remove ligand/receptor complexes (e.g., L17A_IL17F, IL17RA_IL17RC)?")
		st.write("Note that some LR pairs have only comlexes.")
		submitted_basic = st.form_submit_button("Change settings")
	if submitted_basic:
		st.session_state.liana_res = []

	if db == 'cellchat_secreted_signaling': #cellchat secreted signalingを追加
		if species == 'mouse':
			resource = pd.read_csv("db/cellchat_secreted_signaling_mouse.tsv", sep = '\t')
		else:
			resource = pd.read_csv("db/Cellchat_secretory.tsv", sep = '\t')
	elif db == 'nichenet':
		if species == 'mouse':
			resource = pd.read_csv("db/nichenet.tsv", sep = '\t')
		else:
			resorce = pd.read_csv("db/nichenet_human.tsv", sep = '\t')
	elif db == 'cellchatdb': #cellchat secreted signalingを追加
		if species == 'mouse':
			resource = pd.read_csv("db/cellchat_mouse.tsv", sep = '\t')
		else:
			resource = pd.read_csv("db/cellchat.tsv", sep = '\t')
	elif db != 'mouseconsensus' and species == 'mouse':
		geneinfo_human, geneinfo_2022 = load_geneinfo()
#		map_df = read_map_df()
		resource = li.rs.select_resource(db)
		mouse, unconverted = mouse_conversion(resource)
		if len(unconverted) >0:
			st.write("Uncorverted LR")
			st.write(unconverted)
		st.write("Original LR numbers: " + str(len(resource)))
		st.write("Mouse LR numbers: " + str(len(mouse)))
		resource = mouse
#		st.write(mouse.head(100))
	else:
		resource = li.rs.select_resource(db)
	
	st.write("LR numbers: " + str(len(resource)))

	resource = resource.astype({'ligand': str, 'receptor': str}) #型を修正

	if remove_complex:
		resource = resource[~resource['ligand'].str.contains('_') & ~resource['receptor'].str.contains('_')]
		st.write("Filtered LR numbers: " + str(len(resource)))

	st.write(resource.head(5))

	if db == "consensus":
		res_name = "liana" + "_res"
	else:
		res_name = db + "_res"

	st.write(" ")

	if st.button('Run analysis') or (len(st.session_state.liana_res) > 0):
		#optionを変えたときは待つようにする

		if method in ['Rank_Aggregate','CellChat','CellPhoneDB','Connectome','log2FC','NATMI','SingleCellSignalR','Geometric Mean','scSeqComm']:
			rankaggregate_options = {'groupby': groupby, 'resource': resource, "expr_prop": expr_prop,
			"use_raw": False, "key_added": res_name, "n_jobs": 8}

			try: #既存データを除く
				del adata.uns[res_name]
			except:
				pass

			adata = calc_lr(adata, rankaggregate_options, func_name)

#			li.mt.rank_aggregate(adata,**rankaggregate_options)
			st.write(adata.uns[res_name].head())
			st.session_state.liana_res = adata.uns[res_name]
			st.session_state.liana_resart = False

		specificity_col = pval_col[method_list.index(method)]
		magnitude_col = lr_means_list[method_list.index(method)]
		if specificity_col == "None":
			specificity_col = magnitude_col
			st.write("No res for LR specificity by this method.")
		if magnitude_col == "None":
			magnitude_col = specificity_col
			st.write("No res for LR magnitude by this method.")
		st.write("Magnitude score: " + magnitude_col)
		st.write("Specificity score: " + specificity_col)

		dot_inverse_colour = False
		if method =="Rank_Aggregate":
			dot_inverse_colour = True


		if st.button('Show the whole result table'):
			st.write(st.session_state.liana_res)

		ligand_list = sorted(list(set(st.session_state.liana_res['ligand_complex'])))
		receptor_list = sorted(list(set(st.session_state.liana_res['receptor_complex'])))

		with st.form("LR cell selection:"):
			source_labels = st.multiselect('Select cells expressing ligands',cell_list, default = cell_list[0])
			target_labels = st.multiselect('Select cells expressing receptors',cell_list, default = cell_list[0])
			submitted = st.form_submit_button("Set sourc and target cells")
		st.write(" ")

		with st.sidebar:
			st.markdown("#### Options for dotplot")
			dot_type = st.radio("Dotplot of:", ('top n','specific ligands', "specific receptors"), index = 0)
			dot_top_n = None
			dot_source = None
			dot_target = None
			if dot_type == 'top n':
				dot_top_n = st.number_input('Number of top interactions to show', min_value =1, step = 1, value=20)
				dot_down_name = "dotplot_top_" + str(dot_top_n) + '.pdf'
			elif dot_type == 'specific ligands':
				dot_source = st.multiselect('Select ligands to show', ligand_list, default = ligand_list[0])
				dot_down_name = "dotplot_" + "".join(dot_source) + '.pdf'
			else:
				dot_target = st.multiselect('Select receptors to show', receptor_list, default = receptor_list[0])
				dot_down_name = "dotplot_" + "".join(dot_target) + '.pdf'
			dot_filter = st.number_input("Filter p val <", value = 0.05, step = 0.01, min_value = 0.0, max_value = 1.0)
			dot_ranking = st.selectbox("Dotplot order by", ['magnitude (LR expression)', 'specificity (p-val)'], index =0)
			if dot_ranking == 'specificity (p-val)':
				dot_ranking_col = specificity_col
				if method in ["NATMI","Connectome"]:
					dot_orderby_ascending = False
				else:
					dot_orderby_ascending = True
			else:
				dot_ranking_col = magnitude_col
				dot_orderby_ascending = False
			palette = ['viridis','magma','plasma','inferno','turbo','black_yellow','white_green',
			'blue_yellow', 'cool_warm','green_purple']
			dot_color = st.selectbox("Color:", palette, index = 0)
			dot_x_size = st.number_input("Dotplot x size:", value = 8, step = 1, min_value = 1)
			dot_y_size = st.number_input("Dotplot y size:", value = 7, step = 1, min_value = 1)
			st.markdown("---")
			st.markdown("#### Options for tileplot")
			tile_type = st.radio("Tileplot of:", ('top n','specific ligands', "specific receptors"), index = 0)
			tile_top_n = None
			tile_source = None
			tile_target = None
			if tile_type == 'top n':
				tile_top_n = st.number_input('N_top to show', min_value =1, step = 1, value=20)
				tile_down_name = "tileplot_top_" + str(tile_top_n) + '.pdf'
			elif tile_type == 'specific ligands':
				tile_source = st.multiselect('Select ligands', ligand_list, default = ligand_list[0])
				tile_down_name = "tileplot_" + "".join(tile_source) + '.pdf'
			else:
				tile_target = st.multiselect('Select receptors', receptor_list, default = receptor_list[0])
				tile_down_name = "tileplot_" + "".join(tile_target) + '.pdf'
			tile_x_size = st.number_input("tileplot x size:", value = 8, step = 1, min_value = 1)
			tile_y_size = st.number_input("tileplot y size:", value = 7, step = 1, min_value = 1)
			st.markdown("---")
			st.markdown("#### Options for circosplot")
			use_strength_var = st.radio("Interaction strength:", ('Sum of interaction magnitude','Count of interactions above specificity threshold'))
			if use_strength_var == 'Sum of interaction magnitude':
				use_strength = "interaction_strength"
			else:
				use_strength = "lr_pair_count"
			LR_count_threshold = st.number_input("Threshold for specificity (P-val)", value = 0.05, step = 0.01, min_value = 0.0)
			circos_size = st.number_input("circosplot size:", value = 8, step = 1, min_value = 1)

			highlight_from = []
			highlight_to = []
			with st.form("Circoid from to:"):
				highlight_from = st.multiselect('Highlight from',cell_list, default = None)
				highlight_to = st.multiselect('Highlight to',cell_list, default = None)
				submitted = st.form_submit_button("Highlight from/to cell-types")

		if st.button('Draw dotplot'):

			if method in ["NATMI","Connectome"]: #NATMIのspecificity weightは高いほうがよい。
				dot_inverse_size = False
			else:
				dot_inverse_size = True

			try:


				dotplot = li.pl.dotplot(#adata = adata,
					liana_res = st.session_state.liana_res,
					colour=magnitude_col,
					inverse_colour=dot_inverse_colour,
					size=specificity_col,
					inverse_size=dot_inverse_size, # we inverse sign since we want small p-values to have large sizes
					source_labels=source_labels,
					target_labels=target_labels,
					orderby=dot_ranking_col,
					orderby_ascending=dot_orderby_ascending,
					top_n=dot_top_n,
					ligand_complex = dot_source,
					receptor_complex = dot_target,
					figure_size=(dot_x_size, dot_y_size),
					# finally, since cpdbv2 suggests using a filter to FPs
					# we filter the pvals column to <= 0.05
					filter_fun=lambda x: x[specificity_col] <= dot_filter,
					#uns_key=res_name, # uns_key to use, default is 'liana_res',
					 )


				if palette.index(dot_color) < 5:
					dotplot = dotplot + scale_color_cmap(cmap_name=dot_color)
				elif palette.index(dot_color) < 7:
					color_split = dot_color.split('_')
					dotplot = dotplot + p9.scale_color_gradient(low=color_split[0], high=color_split[1])
				elif dot_color == 'cool_warm':
					dotplot = dotplot + p9.scale_color_gradient(low="#4575B4", high="#D73027")
				elif dot_color == "blue_yellow":
					dotplot = dotplot + p9.scale_color_gradient(low="#352A86", high="#FFFF00")
				elif dot_color == "green_purple":
					dotplot = dotplot + p9.scale_color_gradient(low="#00CC00", high="#CC00CC")

				# ligand等を設定する場合はY axisを逆転する　完全に逆転するというより大きさ順になる印象
#				if (dot_source is not None) or (dot_target is not None):
				# Get unique values of the 'interaction' column
				interactions = dotplot.data['interaction'].unique()
				# Reverse the order
				reversed_interactions = list(reversed(interactions))
				dotplot = dotplot + scale_y_discrete(limits=reversed_interactions)
#				dotplot = dotplot + scale_y_discrete(reverse=True)　#　これは動かない
				fig = dotplot.draw()
				st.pyplot(fig)
				if method not in ["NATMI","Connectome"]:
					st.write("Dot size is -log10 transformed")
				pdf_path = liana_temp_dir + 'dotplot.pdf'
				dotplot.save(pdf_path, bbox_inches='tight')
				with open(pdf_path, "rb") as pdf_file:
					PDFbyte = pdf_file.read()

					st.download_button(label="Download dotplot",
									data=PDFbyte,
									file_name=dot_down_name,
									mime='application/octet-stream')
			except:
				st.write("Nothing to draw. You may change Filter p val < on the side panel.")
		st.write(" ")
		if method != "Rank_Aggregate":
			if st.button('Draw tileplot'):
				tileplot = li.pl.tileplot(#adata = adata,
					# NOTE: fill & label need to exist for both
					# ligand_ and receptor_ columns
					liana_res = st.session_state.liana_res,
					fill=magnitude_col,
					label='props',
					label_fun=lambda x: f'{x:.2f}',
					top_n=tile_top_n,
					orderby=specificity_col,
					orderby_ascending=True,
					source_labels=source_labels,
					target_labels=target_labels,
					ligand_complex = tile_source,
					receptor_complex = tile_target,
					#uns_key=res_name, # NOTE: default is 'liana_res'
					source_title='Ligand',
					target_title='Receptor',
					figure_size=(tile_x_size, tile_y_size),
							 )

				fig = tileplot.draw()
				st.pyplot(fig)
				pdf_path = liana_temp_dir + 'tileplot.pdf'
				tileplot.save(pdf_path, bbox_inches='tight')
				with open(pdf_path, "rb") as pdf_file:
					PDFbyte = pdf_file.read()
					st.download_button(label="Download tileplot",
									data=PDFbyte,
									file_name=dot_down_name,
									mime='application/octet-stream')


		st.write(" ")

		if st.button('Draw circolsplot'):

			interaction_df = gnerate_interaction_matrix(st.session_state.liana_res, LR_count_threshold)
	#		st.write(interaction_df)

			import pycirclize
			from pycirclize import Circos
			from pycirclize.parser import Matrix

			interaction_matrix = interaction_df.pivot(index='source',columns='target',values = use_strength)

			st.write(interaction_matrix)

			matrix_data=interaction_matrix.values
			row_names = list(interaction_matrix.index)
			col_names = row_names

			# Normalize the data
			matrix_df_normalized = interaction_matrix / interaction_matrix.max().max() * 100


			# Define link_kws handler function to customize each link property
			def link_kws_handler(from_label: str, to_label: str):
				if len(highlight_from) > 0 or len(highlight_to) > 0:
					if from_label in highlight_from or to_label in highlight_to:
						return dict(alpha=0.5, zorder=1.0)
					else:
						 return dict(alpha=0.1, zorder=1.0)
				else:
					return dict(alpha=0.5, zorder=1.0)


#			fig, ax = plt.subplots(figsize=(circos_size, circos_size), dpi=150)


			# Initialize from matrix (Can also directly load tsv matrix file)
			circos = Circos.initialize_from_matrix(
				matrix_df_normalized,
				space=3,
				cmap="tab10",
				r_lim=(93, 100),
				ticks_interval=500,
				label_kws=dict(r=94, size=12, color="white"),
				link_kws=dict(direction=1, ec="black", lw=0.5),
				link_kws_handler=link_kws_handler
			)

			#print(matrix_df)
			fig = circos.plotfig()
			fig.set_size_inches(circos_size, circos_size)

			st.pyplot(fig)
			pdf_path = liana_temp_dir + 'circosplot.pdf'
			plt.savefig(pdf_path, bbox_inches='tight')
			plt.close()
			with open(pdf_path, "rb") as pdf_file:
				PDFbyte = pdf_file.read()
				st.download_button(label="Download circosplot",
								data=PDFbyte,
								file_name=method + "_circos.pdf",
								mime='application/octet-stream')


		st.markdown("---")
		if st.button('Rerun calculation'):
			st.session_state.liana_res = []
			st.rerun()



