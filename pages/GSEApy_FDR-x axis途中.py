import streamlit as st
import pandas as pd
import csv
import re
import os
import numpy as np
import matplotlib.pyplot as plt
import gseapy as gp
from io import StringIO
import pickle
import plotly.express as px
from gseapy import gseaplot
from gseapy import dotplot
import shutil
import glob
import networkx as nx
from gseapy import enrichment_map
import glob
from PIL import Image
from natsort import natsorted
from helper_func import clear_old_directories
from helper_func import clear_old_files
import time
import seaborn as sns
from matplotlib.cm import ScalarMappable


st.set_option('deprecation.showPyplotGlobalUse', False)

st.set_page_config(page_title="GSEApy analysis.", page_icon="√")


def file_name_check(file_name): #\は変換がうまく行かない
	file_name = re.sub(r'[/|:|?|.|"|<|>|\ |\\]', '-', file_name)
	return file_name


@st.cache_data
def read_excel(file, index_col=None, header = 0):
	df_xl = pd.read_excel(file, index_col = index_col, header = header)
	return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=',', header = 0):
	df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
	return df_c

@st.cache_data
def convert_df(df):
	return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def get_library_name():
	return gp.get_library_name()

@st.cache_data
def calc_enrich(gene_list, gene_sets, background=None):
	enr = gp.enrich(gene_list=gene_list,
					 gene_sets=gene_sets, # kegg is a dict object
					 background=background, # or "hsapiens_gene_ensembl", or int, or text file, or a list of genes
					 outdir=None,
					 verbose=True)
	return enr

@st.cache_data
def calc_enrichr(gene_list,gene_sets):
	enr = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
					 gene_sets=GO_name,
					 outdir=None, # don't write to disk
					)
	return enr

@st.cache_data
def calc_prerank(rnk, gene_sets, min_size, max_size, permutation_num, seed):
	pre_res = gp.prerank(rnk=rnk, # or rnk = rnk,
					gene_sets=GO,
					threads=20,
					min_size=min_size,
					max_size=max_size,
					permutation_num=permutation_num, # reduce number to speed up testing
					outdir=None, # don't write to disk
					seed=seed,
					verbose=False,
				  )
	return pre_res

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

def add_GO_term(i):
	new_cont = [i, i]
	new_cont.extend(GO[i])
	return new_cont


# --- Initialising SessionState ---
if "gseacalc" not in st.session_state:
	  st.session_state.gseacalc = False
if "gsea" not in st.session_state:
	  st.session_state.gsea = False

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #	st.write('file name kept')
	filename_add = ""

if "gsea_temp_dir" not in st.session_state:
	st.session_state.gsea_temp_dir = True
	gsea_temp_dir = "temp/" + str(round(time.time()))
	if not os.path.exists('temp'):
		os.mkdir('temp')
	else:
		clear_old_directories("temp")
		clear_old_files("temp")
	os.mkdir(gsea_temp_dir)
	st.session_state.gsea_temp_dir = gsea_temp_dir
else:
	gsea_temp_dir = st.session_state.gsea_temp_dir




st.markdown("### GSEA and overrepresentation test by GSEApy")

Analysis_mode = st.radio(
	"Analysis mode:",
	('Prerank','Over-representation'), key='Prerank')

if Analysis_mode == "Over-representation":
	Test_mode = st.radio("Test mode:", ('Hypergeometric test','Enrichr web'), key='Hypergeometric test')
	st.markdown("##### Genes (comma, space, CR separated):")
	genes = st.text_input("genes",label_visibility = 'collapsed')
	gene_list = []
	if len(genes) > 0:
		if ',' not in genes:
			gene_list = genes.split(' ')
		else:
			genes =  ''.join(genes.split()) #空白を除く
			gene_list = genes.split(',')
		genes = list(filter(lambda x: x != "", genes)) #空白を除く
		st.write(gene_list[:3])

	if Test_mode == 'Enrichr web':
		test_name = 'Enrichr_Web'
		GO_list = get_library_name()

		GO_name = st.selectbox('Select gene set', GO_list)

		if st.button('Run Enrichr web analysis'):

			enr = calc_enrichr(gene_list=gene_list, gene_sets=GO_name)
#			enr = gp.enrichr(gene_list=gene_list, # or "./tests/data/gene_list.txt",
#					 gene_sets=GO_name,
#					 outdir=None, # don't write to disk
#					)
			st.dataframe(enr.results)

	else:
		test_name = 'GSEApy_Hypergeometric'
		species = st.radio("Species:", ('mouse','human'))
		db = st.radio("DB:", ('mSigDB','Enrichr', 'Your own GMT file'))

		if db == 'mSigDB':
			if species == 'mouse':
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "mSigDB_mouse")
			else:
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "mSigDB")
		elif db == 'Enrichr':
			if species == 'mouse':
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "enrichr_database_mouse")
			else:
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "enrichr_database")
		else:
			GO = None
			uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt','gmt'])
			if uploaded_gmt is not None:
				GO_name = uploaded_gmt.name
				stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
				s = stringio.read()
				t = s.split('\n')
				gmt =[x.split('\t') for x in t]
				GO = dict()
				for i in gmt:
					GO[i[0]] = i[2:]

		if db != "Your own GMT file":
			files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
			files_file.sort()
			key_index = 0
			GO_name = st.multiselect('Select gene set',files_file, default = files_file[-1])
#			GO_name = st.selectbox('Select gene set',files_file, index=key_index)
			if db == 'mSigDB': #GMT fileからdictへ変換
				GO = dict()
				for i in GO_name:
					GO_file = dir_path + '/' + i
					with open(GO_file) as f:
						s = f.read()
					t = s.split('\n')
					gmt =[x.split('\t') for x in t]
					GO_dic = dict()
					for i in gmt:
						GO_dic[i[0]] = i[2:]
					GO = GO | GO_dic

			else:
				GO = dict()
				for i in GO_name:
					with open(dir_path + '/' + i, 'rb') as handle:
						GO_dic = pickle.load(handle)
					GO = GO | GO_dic


		set_back = st.checkbox('Set background genes?', value=False)
		st.markdown("By default, all genes in the gene sets are used as background. However, all the expressed genes are a better background. To do this, set background genes.")
		if set_back:
			input_file_type = st.radio(
		"Data format:",
		('tsv','csv', 'excel'))
			uploaded_file = st.file_uploader("Upload a file containing gene names (e.g., gene list, DESeq2, Homer)", type=['txt','tsv','csv','xls','xlsx'])
			if uploaded_file is not None:
				if input_file_type == "csv":
					df = read_csv(uploaded_file, header = None, index_col = None)
				elif input_file_type == 'tsv':
					df = read_csv(uploaded_file, sep = '\t', header=None, index_col = None)
				else:
					df = read_excel(uploaded_file, index_col = None, header = None)

				# もし1列のデータで最初にGeneがないとき
				if df.shape[1] == 1:
					bk_genes = df.iloc[:,1].values
					if bk_genes[0] == "Gene" or bk_genes[0] == "GENE":
						bk_genes = bk_genes[1:]

				else:
					df.columns = df.iloc[0,:].tolist()	# transposeすると狂うので、transposeした後にcolumnsを決める
					df = df.drop(0, axis = 0) # 1行目を列名にして除く

					st.write(df.head())
					content = df.columns.tolist()
					Gene_column = content[0]
					if "Annotation/Divergence" in content:
						  # colnamesの変換
						search_word = '([^\ \(]*)\ \(.*'

						for i in range(1, len(content)):
							match = re.search(search_word, content[i])
							if match:
								content[i] = match.group(1).replace(' ', '_')
						df.columns = content # 一旦名前を変更
						df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # excel 対応

						pattern = "([^|]*)"
						repatter = re.compile(pattern)
						f_annotation = lambda x: repatter.match(x).group(1)
						df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
				#		df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
						# annotation/divergence以前を除く
						df = df.loc[:,'Annotation/Divergence':]
						content = df.columns.tolist()
						content[0] = 'Gene'
						df.columns = content
						Gene_column = "Gene"
						st.write("Converted Annotation/Divergence to gene symbols.")

					elif "Gene" in content:
						Gene_column =  "Gene"
					else:
						Gene_column =  st.selectbox(
						'Select gene column',
						content)

					bk_genes = df[Gene_column].values
				st.write(bk_genes[:3])


		if (st.button('Run local hypergeometric test') or st.session_state.gseacalc) and GO:
			st.session_state.gseacalc = True
			if set_back:
				enr = calc_enrich(gene_list=gene_list, gene_sets=GO, background=bk_genes)


			else:
				enr = calc_enrich(gene_list=gene_list, gene_sets=GO, background=None)

			enr.results = enr.results.sort_values('Adjusted P-value', ascending=True)

			st.dataframe(enr.results)
			p_thre = 0.05
			num_bar = 10
			py_x_size = 600
			py_y_size = 400
			with st.sidebar:
				p_thre = st.number_input('Visuallzation threshold for adj. P', min_value =0.0, step = 0.01, value=0.05)
				num_bar = int(st.number_input('Number of terms to visualize', min_value =1, step = 1, value=10))

				py_x_size = int(st.number_input("Plot x size:", value = 600, step = 100, min_value = 100))
				py_y_size = int(st.number_input("Plot y size:", value = 400, step = 100, min_value = 100))

			sig = enr.results.loc[(enr.results['Adjusted P-value']<p_thre),:]

			if len(sig) == 0:
				st.markdown("### Nothing passed adjuste P < " + str(p_thre))
			else:
				sig['-log10adjP'] = -np.log10(sig['Adjusted P-value'])
				if len(sig) < num_bar:
					chart_len = len(sig)
				else:
					chart_len = num_bar
				sig_sub = sig.iloc[:chart_len,:]
				sig_sub = sig_sub.sort_values('Combined Score', ascending = True)
				sig_sub = sig_sub.sort_values('-log10adjP', ascending = True)

				fig = px.bar(sig_sub, y="Term", x="-log10adjP", color ='Combined Score', orientation='h',  width=py_x_size, height=py_y_size)

				st.plotly_chart(fig,use_container_width=True)


			csv = convert_df(enr.results)
			GO_name = "_".join(GO_name).replace('.gmt','')
			st.download_button(
				"Press to Download",
				csv,
				GO_name + "." + test_name + '.tsv',
				"text/csv",
				key='download-csv'
			)

# prerank
else:
	uploaded_rnk = st.file_uploader("Upload rnk file", type=['rnk','rank','txt'])
	if uploaded_rnk is not None:
		rnk_name = uploaded_rnk.name
		rnk = read_csv(uploaded_rnk, index_col=0, sep='\t', header = None)

		species = st.radio("Species:", ('mouse','human'))
		db = st.radio("DB:", ('mSigDB','Enrichr', 'Your own GMT file'))

		if db == 'mSigDB':
			if species == 'mouse':
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "mSigDB_mouse")
			else:
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "mSigDB")
		elif db == 'Enrichr':
			if species == 'mouse':
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "enrichr_database_mouse")
			else:
				dir_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "enrichr_database")
		else:
			uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt','gmt'])
			if uploaded_gmt is not None:
				GO_name = uploaded_gmt.name
				stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
				s = stringio.read()
				t = s.split('\n')
				gmt =[x.split('\t') for x in t]
				GO = dict()
				for i in gmt:
					GO[i[0]] = i[2:]

		if db != "Your own GMT file":
			files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
			files_file.sort()
			key_index = 0
#			if db == 'mSigDB':
#				key_index=len(files_file) - 1

			GO_name = st.multiselect('Select gene set',files_file, default = files_file[-1])
#			GO_name = st.selectbox('Select gene set',files_file, index=key_index)
			if db == 'mSigDB': #GMT fileからdictへ変換
				GO = dict()
				for i in GO_name:
					GO_file = dir_path + '/' + i
					with open(GO_file) as f:
						s = f.read()
					t = s.split('\n')
					gmt =[x.split('\t') for x in t]
					GO_dic = dict()
					for i in gmt:
						GO_dic[i[0]] = i[2:]
					GO = GO | GO_dic

			else:
				GO = dict()
				for i in GO_name:
					with open(dir_path + '/' + i, 'rb') as handle:
						GO_dic = pickle.load(handle)
					GO = GO | GO_dic

		# termの_をスペースへ　　enrichment graphでタイトルがラップされないため
		for i in list(GO.keys()):
			GO[i.replace("_"," ")] = GO.pop(i)

		min_size = 15
		max_size = 500
		seed = 6
		permutation_num = 1000
		with st.sidebar:
			st.markdown("#### GSEA parameters")
			min_size = int(st.number_input('Minimum number of genes in a gene set', min_value =0, step =1, value=15))
			max_size = int(st.number_input('Maximum number of genes in a gene set', min_value =0, step =100, value=500))
			permutation_num = int(st.number_input('Number of permutation', min_value =100, step =100, value=1000))
			seed = int(st.number_input('Seed for permutation (0: time stamp)', min_value =0, step =1, value=6))
			if seed == 0:
			 	import time
			 	seed = int(time.time())
			rev_rnk = st.checkbox("Reverse rank order?", value = False)
			if rev_rnk:
			 	rnk[1] = -rnk[1]
			 	rev_rnk_name = st.text_input("New rank name:", rnk_name + "_rev")
			 	rnk_name = rev_rnk_name

		GO_name_dir = "_".join(GO_name).replace('.gmt','')
#		GO_name_dir  = GO_name.replace(".gmt","")
		rnk_name_dir = rnk_name.replace(".rnk",'')
		rnk_name_dir = rnk_name_dir.replace(".rank",'')


		gsea_dir = gsea_temp_dir + "/" + rnk_name_dir + "_" + GO_name_dir # dir作成
		st.write(gsea_dir)
	#	if os.path.exists(gsea_dir) and not st.session_state.gsea:
	#		shutil.rmtree(gsea_dir)
	#		st.write("GSEA_dir exists")
		if not os.path.exists(gsea_dir):
			st.write("GSEA_dir does not exists")
			if not os.path.exists('temp'):
				os.mkdir('temp')
			os.mkdir(gsea_dir)
			os.mkdir(gsea_dir + '/upregulated_enrichment')
			os.mkdir(gsea_dir + '/downregulated_enrichment')

		if st.button('Run prerank GSEA test') and GO: # 状態にかかわらず、ボタンが押されたら計算する。
			st.session_state.gsea = True

			if list(rnk.index.duplicated()).count(True) > 0:
				st.markdown("#### There are duplicated genes in the rank file.")
				st.write('Dupliated genes:' +  ', '.join(list(rnk[rnk.index.duplicated()].index)))
				st.write("The first instances will be kept.")
				st.markdown("---")
				rnk = rnk[~rnk.index.duplicated(keep='first')]


			pre_res = calc_prerank(rnk=rnk, gene_sets=GO, min_size=min_size, max_size=max_size,
							permutation_num=permutation_num, seed=seed)

			st.session_state.pre_res = pre_res
			terms = list(pre_res.results.keys())
			gsea_res = pd.DataFrame(columns = ['SIZE','ES','NES','NOM p-val','FDR q-val','FWER p-val','TAG %','GENE %', 'LEADING GENES'])
			gsea_res.index.name = "NAME"
			for i in terms:
				cont = pre_res.results[i]
				gsea_res.loc[i] = [len(GO[i]), cont['es'],cont['nes'],cont['pval'],cont['fdr'],
				cont['fwerp'],cont['tag %'],cont['gene %'],cont['lead_genes']]
			gsea_res = gsea_res.sort_values('FDR q-val', ascending=True)
		#	gsea_res.to_csv('temp.tsv', sep = '\t')
			up_res = gsea_res.loc[gsea_res['ES'] > 0]
			up_res =  up_res.sort_values('FDR q-val', ascending=True)
			down_res = gsea_res.loc[gsea_res['ES'] < 0]
			down_res =  down_res.sort_values('FDR q-val', ascending=True)
			st.session_state.gsea_res = gsea_res
			st.session_state.up_res = up_res
			st.session_state.down_res = down_res  #必要なデータは状態変数に入れておく
			st.session_state.terms = terms

		if st.session_state.gsea: # 一度計算した後は、ここからスタート
			pre_res = st.session_state.pre_res
			gsea_res = st.session_state.gsea_res
			up_res = st.session_state.up_res
			down_res = st.session_state.down_res

			terms = st.session_state.terms
			geneset_num = len(terms)


			st.markdown("##### Upregulated gene sets")
			st.write(str(len(up_res)) + " / " + str(geneset_num))
			st.write(str([up_res.loc[x,'FDR q-val'] < 0.25 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 25%')
			st.write(str([up_res.loc[x,'FDR q-val'] < 0.1 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 10%')
			st.write(str([up_res.loc[x,'FDR q-val'] < 0.05 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 5%')

			st.dataframe(up_res)

			st.markdown("##### Downregulated gene sets")
			st.write(str(len(down_res)) + " / " + str(geneset_num))


			st.write(str([down_res.loc[x,'FDR q-val'] < 0.25 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 25%')
			st.write(str([down_res.loc[x,'FDR q-val'] < 0.1 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 10%')
			st.write(str([down_res.loc[x,'FDR q-val'] < 0.05 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 5%')

			st.dataframe(down_res)

			term_thre = st.number_input('Terms with q <', min_value =0.0, value=0.1)

			select_term = gsea_res.loc[gsea_res['FDR q-val'] < term_thre].index.values
			enrichment_term = st.selectbox('Select term to draw enrichment graph', select_term)
			es_x_size = 5
			es_y_size = 5
			with st.sidebar:
				st.markdown("#### Enrichment plot parameters")
				es_x_size = int(st.number_input("Enrichment plot x size:", value = 5, step = 1, min_value = 2))
				es_y_size = int(st.number_input("Ebrichment plot y size:", value = 5, step = 1, min_value = 2))

			if enrichment_term:

				g = pre_res.plot(terms=enrichment_term, figsize=(es_x_size,es_y_size))
				st.pyplot(g)

				gseaplot(rank_metric=pre_res.ranking, term=enrichment_term,
					ofname=gsea_dir + '/' + file_name_check(enrichment_term) + '.pdf', #スペースを_へ
					**pre_res.results[enrichment_term], figsize=(es_x_size,es_y_size))
			else:
				st.markdown("#### No terms with FDR < " + str(term_thre) + " for enrichment plot.")

			res2d = pd.DataFrame(gsea_res.loc[:,["FDR q-val","LEADING GENES","NES"]])
			res2d['Term'] = res2d.index.values
			res2d.columns = ['FDR q-val', 'Lead_genes', 'NES', 'Term']

#			dot_fdr  = 0.25
#			dot_size = 6
#			dot_x_size = 8
#			dot_y_size = 8
			with st.sidebar:
				st.markdown("#### Dotplot parameters")
				dot_fdr  = st.number_input("Dotplot FDR threshold:", value = 0.1, min_value = 0.001)
				dot_num  = st.number_input("Max number of terms to show:", value = 12, min_value = 1)
				dot_pos = st.selectbox('Show terms:', ('Both',"Upregulated","Downregulated"), index = 0)
				dot_size = int(st.number_input("Dot size:", value = 5, step = 1, min_value = 1))
				dot_x_size = int(st.number_input("Dotplot x size:", value = 8, step = 1, min_value = 2))
				dot_y_size = int(st.number_input("Dotplot y size:", value = 8, step = 1, min_value = 2))

			st.write(" ")

			res2d_dot = res2d.copy(deep = True)
			if dot_pos == "Upregulated":
				res2d_dot = res2d_dot.loc[res2d_dot['NES'] > 0]
			elif dot_pos == "Downregulated":
				res2d_dot = res2d_dot.loc[res2d_dot['NES'] < 0]

			res2d_dot = res2d_dot.sort_values('FDR q-val', ascending = True)

			if len(res2d_dot) > dot_num:
				res2d_dot = res2d_dot.iloc[0:dot_num]

			if len(res2d_dot.loc[res2d_dot['FDR q-val'] < dot_fdr]) > 1:
				try:
					p = dotplot(res2d_dot, column="FDR q-val",
						x = 'NES',
					 title=GO_name_dir,
					 cmap=plt.cm.viridis,
					 size=dot_size, # adjust dot size
					 figsize=(dot_x_size, dot_y_size), cutoff=dot_fdr, show_ring=False)

					st.pyplot(p.figure)
					dotplot(res2d_dot, column="FDR q-val",
						x='NES',
					 title=GO_name_dir,
					 cmap=plt.cm.viridis,
					 size=dot_size, # adjust dot size
					 figsize=(dot_x_size,dot_y_size), cutoff=dot_fdr, show_ring=False, ofname=gsea_dir + '/' + GO_name_dir + '.dotplot.pdf',)
				except:
					st.write("Error in generating dotplot.")
			else:
				st.markdown("#### No or only one term with FDR < " + str(dot_fdr) + " for dot plot.")


			with st.sidebar:
				st.markdown("#### Barplot parameters")
				bar_vis_thr  = st.number_input("Barplot FDR threshold:", value = 0.05, min_value = 0.001)
				bar_num  = st.number_input("Max number of terms for barplot", value = 12, min_value = 1)
				bar_type = st.selectbox('Barplot type:', ('Single plot','Separete up/down plots'), index = 0)
				bar_FDR = st.checkbox('FDR as X-axis?', value=False)
				bar_cmap = st.selectbox('Barplot color map:', ('Default', 'Greys_r', 'Blues_r', 'BrBG_r', 'BuGn_r', 'BuPu_r','GnBu_r', 'Greens_r','OrRd_r', 'Oranges_r', 'autumn', 'binary_r', 'bone',  'gist_gray' 'magma_r' 'viridis'), index = 0)
				bar_x_size = int(st.number_input("Barplot x size:", value = 8, step = 1, min_value = 2))
				bar_y_size = int(st.number_input("Barplot y size:", value = 8, step = 1, min_value = 2))

			st.write(" ")

			#ーーーーーーーーーーbarplot
			gsea_nes = gsea_res.copy(deep=True)
			gsea_nes["Gene set"] = gsea_nes.index.to_list() # indexはsnsで使えない様子
			gsea_nes = gsea_nes.sort_values('NES', ascending=False)
			gsea_nes = gsea_nes.loc[gsea_nes['FDR q-val']<bar_vis_thr]
			up_nes = gsea_nes.loc[gsea_nes['NES'] > 0]
			down_nes = gsea_nes.loc[gsea_nes['NES'] < 0]

			if bar_type == "Separete up/down plots":
				if bar_FDR:
					NES_max = max(gsea_nes["NES"])
					NES_min = min(gsea_nes["NES"])
					gsea_nes['-log10FDR'] = np.log10(1/gsea_nes['FDR q-val'])
					gsea_nes = gsea_nes.sort_values('FDR q-val', ascending=True)
					up_nes = gsea_nes.loc[gsea_nes['NES'] > 0]
					down_nes = gsea_nes.loc[gsea_nes['NES'] < 0]
				if len(up_nes) >0:
					fig, ax = plt.subplots(figsize=(bar_x_size, bar_y_size))
					if len(up_nes) > bar_num:
						up_nes = up_nes.iloc[0:bar_num]
					if bar_cmap == "Default":
						bar_cmap_use = "OrRd_r"
					else:
						bar_cmap_use = bar_cmap
					if bar_FDR:
						if "_r" in bar_cmap_use:
							bar_cmap_use = bar_cmap_use.replace('_r',"")
						else:
							bar_cmap_use = bar_cmap_use + "_r"
					pal = sns.color_palette(bar_cmap_use, as_cmap=True)
					if bar_FDR:
						ax =  sns.barplot(data= up_nes,  y='Gene set', x="-log10FDR", palette=pal(up_nes["NES"]/NES_max  ))
					else:
						ax =  sns.barplot(data= up_nes,  y='Gene set', x="NES",  palette=pal(up_nes["FDR q-val"]/bar_vis_thr  ))
					ax.set_ylabel('')
					ax.spines['top'].set_visible(False)
					ax.spines['right'].set_visible(False)
					ax.yaxis.set_ticks_position('none')
					if bar_cmap:
						sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,NES_max))
					else:
						sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,bar_vis_thr))
					sm.set_array([])

					cax = plt.axes([1,0.2, 0.04, 0.6])
					cbar  = plt.colorbar(sm,  cax=cax)
					if bar_cmap:
						cbar.set_label('NES', rotation=270,labelpad=25)
					else:
						cbar.set_label('Adj P-value', rotation=270,labelpad=25)
						cax.invert_yaxis()
					st.pyplot(fig)
					fig.savefig(gsea_dir + '/' + "Barplot_up.pdf", format='pdf', bbox_inches='tight')

				if len(down_nes) >0:
					fig, ax = plt.subplots(figsize=(bar_x_size, bar_y_size))
					if len(down_nes) > bar_num:
						down_nes = down_nes.iloc[0:bar_num]
					if bar_cmap == "Default":
						bar_cmap_use = "Blues_r"
					else:
						bar_cmap_use = bar_cmap
					pal = sns.color_palette(bar_cmap_use, as_cmap=True)
					if bar_FDR:
						ax =  sns.barplot(data= down_nes,  y='Gene set', x="-log10FDR", palette=pal(-down_nes["NES"]/NES_min  ))
					else:
						ax =  sns.barplot(data= down_nes,  y='Gene set', x="NES",  palette=pal(down_nes["FDR q-val"]/bar_vis_thr  ))
					ax.set_ylabel('')
					ax.spines['top'].set_visible(False)
					ax.spines['left'].set_visible(False)
					ax.yaxis.set_ticks_position('none')
					if bar_cmap:
						sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0, NES_min))
					else:
						sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,bar_vis_thr))
					sm.set_array([])

					cax = plt.axes([1,0.2, 0.04, 0.6])
					cbar  = plt.colorbar(sm,  cax=cax)
					if bar_cmap:
						cbar.set_label('NES', rotation=270,labelpad=25)
						cax.invert_yaxis()
					else:
						cbar.set_label('Adj P-value', rotation=270,labelpad=25)
						cax.invert_yaxis()
					st.pyplot(fig)
					fig.savefig(gsea_dir + '/' + "Barplot_down.pdf", format='pdf', bbox_inches='tight')


			else:
				if len(gsea_nes) >0:
					fig, ax = plt.subplots(figsize=(bar_x_size, bar_y_size))
					if len(gsea_nes) > bar_num: #NESの絶対値でtop をとる
						gsea_nes['abs NES'] = np.abs(gsea_nes['NES'] )
						gsea_nes = gsea_nes.sort_values('abs NES', ascending=False)
						gsea_nes = gsea_nes.iloc[0:bar_num]
						gsea_nes = gsea_nes.sort_values('NES', ascending=False)
					if bar_cmap == "Default":
						bar_cmap_use = "viridis"
					else:
						bar_cmap_use = bar_cmap
					pal = sns.color_palette(bar_cmap_use, as_cmap=True)
					ax =  sns.barplot(data= gsea_nes,  y='Gene set', x="NES",  palette=pal(gsea_nes["FDR q-val"]/bar_vis_thr  ))
					ax.set_ylabel('')
					ax.spines['top'].set_visible(False)
					ax.spines['right'].set_visible(False)
					ax.yaxis.set_ticks_position('none')
					sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,bar_vis_thr))
					sm.set_array([])

					cax = plt.axes([1,0.2, 0.04, 0.6])
					cbar  = plt.colorbar(sm,  cax=cax)
					cbar.set_label('Adj P-value', rotation=270,labelpad=25)
					cax.invert_yaxis()
					st.pyplot(fig)
					fig.savefig(gsea_dir + '/' + "Barplot.pdf", format='pdf', bbox_inches='tight')




			# return two dataframe
			with st.sidebar:
				st.markdown("#### Enrichment map parameters")
				map_fdr  = st.number_input("Enrichment map FDR threshold:", value = 0.1, min_value = 0.0)
				map_num  = st.number_input("Number of top terms to show:", value = 10, min_value = 1)
				map_pos = st.selectbox('Show terms of:', ('Both',"Upregulated","Downregulated"), index = 0)
				map_node  = st.number_input("Node size:", value = 800, min_value = 100)
				node_font = st.number_input("Node label size:", value = 10, min_value = 1)
				node_cmap = st.selectbox('Node color map:', ('Accent', 'Blues', 'BrBG', 'BuGn', 'BuPu', 'CMRmap', 'Dark2', 'GnBu', 'Greens', 'Greys', 'OrRd', 'Oranges', 'PRGn', 'Paired', 'Pastel1', 'Pastel2', 'PiYG', 'PuBu', 'PuBuGn', 'PuOr', 'PuRd', 'Purples', 'RdBu', 'RdGy', 'RdPu', 'RdYlBu', 'RdYlGn', 'Reds', 'Set1', 'Set2', 'Set3', 'Spectral', 'Wistia', 'YlGn', 'YlGnBu', 'YlOrBr', 'YlOrRd', 'afmhot', 'autumn', 'binary', 'bone', 'brg', 'bwr', 'cividis', 'cool', 'coolwarm', 'copper', 'cubehelix', 'flag', 'gist_earth', 'gist_gray', 'gist_heat', 'gist_ncar', 'gist_rainbow', 'gist_stern', 'gist_yarg', 'gnuplot', 'gnuplot2', 'gray', 'hot', 'hsv', 'inferno', 'jet', 'magma', 'nipy_spectral', 'ocean', 'pink', 'plasma', 'prism', 'rainbow', 'seismic', 'spring', 'summer', 'tab10', 'tab20', 'tab20b', 'tab20c', 'terrain', 'turbo', 'twilight', 'twilight_shifted', 'viridis', 'winter'), index = 11)
				c_r = st.checkbox("Reverse color order?", value=False)
				if c_r:
					node_cmap = node_cmap + "_r"
				map_x_size = int(st.number_input("Enrichment plot x size:", value = 8, step = 1, min_value = 2))
				map_y_size = int(st.number_input("Ebrichment plot y size:", value = 8, step = 1, min_value = 2))


			# return two dataframe
			try:
				res2d_map = res2d.copy(deep = True)
				if map_pos == "Upregulated":
					res2d_map = res2d_map.loc[res2d_map['NES'] > 0]
				elif map_pos == "Downregulated":
					res2d_map = res2d_map.loc[res2d_map['NES'] < 0]
				nodes, edges = enrichment_map(res2d_map, column = 'FDR q-val', cutoff = map_fdr, top_term = map_num)
				#nodes, edges = enrichment_map(res2d.loc[res2d['FDR q-val'] < map_fdr])
				# build graph
				G = nx.from_pandas_edgelist(edges,
								source='src_idx',
								target='targ_idx',
								edge_attr=['jaccard_coef', 'overlap_coef', 'overlap_genes'])
				fig, ax = plt.subplots(figsize=(map_x_size, map_y_size))


				# 以下、間違いを修正済み
				# Gで使用されているtermはnodesの一部
				fig, ax = plt.subplots(figsize=(map_x_size, map_y_size))

				# init node cooridnates
				pos=nx.layout.spiral_layout(G)
				#node_size = nx.get_node_attributes()
				# draw node
				nx.draw_networkx_nodes(G,
									   pos=pos,
									   cmap=node_cmap,
									   node_color=nodes.loc[list(G),'NES'],
									   node_size=nodes.loc[list(G),'Hits_ratio'] * map_node)
				# draw node label
				nx.draw_networkx_labels(G,
										pos=pos,font_size = node_font,
										labels=nodes.loc[list(G),'Term'])
				# draw edge
				edge_weight = nx.get_edge_attributes(G, 'jaccard_coef').values()
				nx.draw_networkx_edges(G,
									   pos=pos,
									   width=list(map(lambda x: x*10, edge_weight)),
									   edge_color='#CDDBD4')

				plt.gca().spines['right'].set_visible(False)
				plt.gca().spines['top'].set_visible(False)
				plt.gca().spines['bottom'].set_visible(False)
				plt.gca().spines['left'].set_visible(False)

				st.pyplot(fig)

				st.markdown('##### Network of gene sets that share leading edge genes.\nNode color: NES, Edge weigth: Jaccard coefficient')

				fig.savefig(gsea_dir + '/' + "enrichment_map.pdf", format='pdf', bbox_inches='tight')
			except:
				st.markdown("#### No terms for network.")



			if st.button('Prepare result files to download'):
				progress_text = "Operation in progress. Please wait."
				my_bar = st.progress(0, text=progress_text)
				basic_stat = "Upregulated gene sets:\n" + str(len(up_res)) + " / " + str(geneset_num) + '\n' + str([up_res.loc[x,'FDR q-val'] < 0.25 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 25%\n' + str([up_res.loc[x,'FDR q-val'] < 0.1 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 10%\n' +str([up_res.loc[x,'FDR q-val'] < 0.05 for x in up_res.index.values ].count(True)) + ' gene sets are significant at FDR < 5%\n' + "\nDownregulated gene sets:\n" + str(len(down_res)) + " / " + str(geneset_num) + '\n' + str([down_res.loc[x,'FDR q-val'] < 0.25 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 25%\n' +str([down_res.loc[x,'FDR q-val'] < 0.1 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 10%\n' +str([down_res.loc[x,'FDR q-val'] < 0.05 for x in down_res.index.values ].count(True)) + ' gene sets are significant at FDR < 5%\n'
				with open(gsea_dir  + '/basic_stat.txt', mode='w') as f:
					f.write(basic_stat)
				num_up_fig = 16
				if num_up_fig > len(up_res):
					num_up_fig = len(up_res)
				num_down_fig = 16
				if num_down_fig > len(down_res):
					num_down_fig = len(down_res)
				percent_count = 0.5 / (num_up_fig + num_down_fig + 18)
				percent_complete = 0
				for i in range(num_up_fig):
					enrichment_term = up_res.index.values[i]
					gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, ofname=gsea_dir + '/upregulated_enrichment/' + file_name_check(enrichment_term) + '.pdf',
						**pre_res.results[enrichment_term], figsize=(es_x_size,es_y_size))
					percent_complete = percent_complete + percent_count
					my_bar.progress(percent_complete, text=progress_text)
				for i in range(num_up_fig):
					enrichment_term = up_res.index.values[i]
					gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, ofname=gsea_dir + '/downregulated_enrichment/' + file_name_check(enrichment_term)  + '.pdf',
						**pre_res.results[enrichment_term], figsize=(es_x_size,es_y_size))
					percent_complete = percent_complete + percent_count
					my_bar.progress(percent_complete, text=progress_text)
				if not os.path.exists(gsea_dir + '/up_png'):
					os.mkdir(gsea_dir + '/up_png')
				if not os.path.exists(gsea_dir + '/down_png'):
					os.mkdir(gsea_dir + '/down_png')

				for i in range(num_up_fig):
					enrichment_term = up_res.index.values[i]
					gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, ofname=gsea_dir +'/up_png/' + str(i) + file_name_check(enrichment_term)  + '.png',
						**pre_res.results[enrichment_term], figsize=(es_x_size,es_y_size))
					percent_complete += percent_count
					my_bar.progress(percent_complete, text=progress_text)
				for i in range(num_down_fig):
					enrichment_term = down_res.index.values[i]
					gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, ofname=gsea_dir + '/down_png/' + str(i) + file_name_check(enrichment_term)  + '.png',
						**pre_res.results[enrichment_term], figsize=(es_x_size,es_y_size))
					percent_complete += percent_count
					my_bar.progress(percent_complete, text=progress_text)


				files = glob.glob(gsea_dir +'/up_png/*.png')
				# タイル状に pm × pm 枚配置
				pm = 4
				d = []
				for i in natsorted(files):
					img = Image.open(i)
					img = np.asarray(img)
					#img = cv2.resize(img, (300, 300), cv2.INTER_LANCZOS4)
					d.append(img)
				fig, ax = plt.subplots(pm, pm, figsize=(16, 16))
				fig.subplots_adjust(hspace=0, wspace=0)
				less_fig = False
				for i in range(pm):
					for j in range(pm):
						try:
							ax[i, j].xaxis.set_major_locator(plt.NullLocator())
							ax[i, j].yaxis.set_major_locator(plt.NullLocator())
							ax[i, j].imshow(d[pm*i+j], cmap="bone")
						except:
							less_fig = True
						ax[i, j].axis('off')
				plt.show()
				st.markdown("#### Upregulated")
				if less_fig:
					st.write("Less than 16 images.")
				st.pyplot(fig)
				fig.savefig(gsea_dir + '/upregulated_enrichment.png', bbox_inches='tight')


				files = glob.glob(gsea_dir +'/down_png/*.png')
				# タイル状に pm × pm 枚配置
				pm = 4
				d = []
				for i in natsorted(files):
					img = Image.open(i)
					img = np.asarray(img)
					#img = cv2.resize(img, (300, 300), cv2.INTER_LANCZOS4)
					d.append(img)
				fig, ax = plt.subplots(pm, pm, figsize=(16, 16))
				fig.subplots_adjust(hspace=0, wspace=0)
				less_fig = False
				for i in range(pm):
					for j in range(pm):
						try:
							ax[i, j].xaxis.set_major_locator(plt.NullLocator())
							ax[i, j].yaxis.set_major_locator(plt.NullLocator())
							ax[i, j].imshow(d[pm*i+j], cmap="bone")
						except:
							less_fig = True
						ax[i, j].axis('off')
				st.markdown("#### Downregulated")
				if less_fig:
					st.write("Less than 16 images.")
				st.pyplot(fig)
				fig.savefig(gsea_dir + '/downregulated_enrichment.png', bbox_inches='tight')


				gsea_res.to_csv(gsea_dir + '/' + 'gsea_report.tsv', sep = '\t')
				up_res.to_csv(gsea_dir + '/' + 'gsea_report_for_na_pos.tsv', sep = '\t')
				down_res.to_csv(gsea_dir + '/' + 'gsea_report_for_na_neg.tsv', sep = '\t')

				st.markdown("##### For large gene sets, it may take some time.")

				if not os.path.exists(gsea_dir + '/edb'):
					os.mkdir(gsea_dir + '/edb')
				rnk.to_csv(gsea_dir +'/edb/' + rnk_name, sep = '\t',header=False)
				# gmt fileの作成
#				gene_sets_gmt = []
#				for i in list(GO.keys()):
#					new_cont = [i, i]
#					new_cont.extend(GO[i])
#					gene_sets_gmt.append(new_cont)
				gene_sets_gmt = [add_GO_term(i) for i in list(GO.keys())] #時間がかかるので内包表記にする
				gmt_str = ""
				percent_complete += percent_count
				my_bar.progress(percent_complete, text=progress_text)
				for i in gene_sets_gmt:
					gmt_str = gmt_str + '\t'.join(i) + '\n'

				percent_complete += percent_count
				my_bar.progress(percent_complete, text=progress_text)
				with open(gsea_dir +'/edb/gene_set.gmt', 'w') as f:
					f.writelines(gmt_str)

				shutil.rmtree(gsea_dir +'/up_png/')
				shutil.rmtree(gsea_dir +'/down_png/')

				down_name = rnk_name_dir + "_" + GO_name_dir
				shutil.make_archive(gsea_temp_dir + "/" + down_name, format='zip',root_dir= gsea_dir)

				my_bar.empty()

				with open(gsea_temp_dir + "/" + down_name + '.zip', "rb") as fp:
					btn = st.download_button(
						label="Download Results",
					data=fp,
					file_name=down_name + ".zip",
					mime = "zip"
					)
