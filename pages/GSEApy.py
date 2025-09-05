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
from helper_func import check_species_index
import time
import seaborn as sns
from matplotlib.cm import ScalarMappable
import sys
from jellyfish import jaro_winkler_similarity
import subprocess
import json


# st.set_option('deprecation.showPyplotGlobalUse', False)
st.set_page_config(page_title="GSEApy analysis.", page_icon="√")

# ボタンスタイルを定義（黄色背景に黒テキスト）
st.markdown("""
<style>
    div.stButton > button {
        background-color: #FFD700;  /* 明るい黄色 */
        color: #000000;  /* 黒テキスト */
        font-weight: 800;  /* テキストをより太く */
        border: none;
        padding: 0.5rem 1rem;
        border-radius: 0.3rem;
        font-size: 1.1rem;  /* テキストサイズを大きく */
    }
    div.stButton > button:hover {
        background-color: #FFC000;  /* ホバー時は少し暗い黄色 */
        box-shadow: 0 4px 8px rgba(0,0,0,0.2);
    }
</style>
""", unsafe_allow_html=True)


# プロットスタイルを設定する関数
def set_plot_style():
    plt.style.use('default')  # デフォルトスタイルをリセット
    plt.rcParams['figure.facecolor'] = 'white'  # 図の背景色を白に
    plt.rcParams['axes.facecolor'] = 'white'    # プロット領域の背景色を白に
    sns.set_style("white")                      # seabornのスタイルを白背景に
    
# グラフ生成前に呼び出す
set_plot_style()


def file_name_check(file_name): #\は変換がうまく行かない
    file_name = re.sub(r'[/|:|?|.|"|<|>|\ |\\]', '-', file_name)
    return file_name


#March-1 Sept-1対応
def excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.write("There are Excel-autoconverted gene names")
                k = 1
            autoconvert_flag = True
            st.write("Converting " + i)
            if x.group(2) == "Mar":
                index_name[j] = "March" + x.group(1)
            elif x.group(2) == "Sep":
                index_name[j] = "Sept" + x.group(1)
        j += 1
    dfx.index = index_name
    return(dfx)


@st.cache_data
def read_xl(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv2(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c


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
    # 実際の遺伝子セットを事前にフィルタリング
    filtered_gene_sets = {}
    small_sets = {}
    
    # 遺伝子セットの総数を表示
    st.write(f"Total number of gene sets: {len(gene_sets)}")
#    st.write("Gene sets dictionary keys:")
#    st.write(list(gene_sets.keys()))
    
    for term, genes in gene_sets.items():
        # 名前を文字列に変換し、空の場合は代替名を使用
        term_str = str(term) if term else "Unnamed Set"
        
        # 遺伝子リストの妥当性チェック
        if not isinstance(genes, (list, set, tuple)) or len(genes) == 0:
            small_sets[term_str] = 0
            continue
            
        # min_sizeより小さい遺伝子セットをチェック
        if len(genes) < min_size:
            small_sets[term_str] = len(genes)
            continue
            
        # 有効な遺伝子セットを追加
        filtered_gene_sets[term] = genes
    
    # 小さな遺伝子セットの報告をexpandableに変更
    if small_sets:
        with st.expander(f"⚠️ {len(small_sets)} gene sets have fewer genes than min_size ({min_size}) or are empty (click to view details). You may change the threshold on the side panel."):
            st.write("The following gene sets were excluded from the analysis:")
            for term, count in sorted(small_sets.items(), key=lambda x: x[1]):
                st.info(f"- {term}: {count} genes")
    
    # フィルタリングされた遺伝子セットがなければエラー
    if not filtered_gene_sets:
        st.error("No gene sets remain after filtering. Please reduce min_size or use different gene sets.")
        return None
        
    # 1つしか遺伝子セットがない場合の警告
    if len(filtered_gene_sets) == 1:
        st.warning("Only one gene set remains after filtering. This may cause errors in GSEApy.")
        st.info("Consider using fgsea (R) method instead or adding more gene sets.")
        return None
    
    # フィルタリングされた遺伝子セットを使用
    try:
        # GSEApyが1つの遺伝子セットでエラーを起こす可能性があるための対策
        if len(filtered_gene_sets) == 1 and 'prerank_method' in st.session_state and st.session_state.prerank_method == 'GSEApy':
            # ダミーの遺伝子セットを追加（実際の分析には使用されない）
            dummy_term = "DUMMY_GENESET_FOR_COMPATIBILITY"
            filtered_gene_sets[dummy_term] = list(rnk.index[:min_size])  # ランキングの最初のmin_size個の遺伝子
            st.info(f"Added a dummy gene set '{dummy_term}' to prevent errors. This will not affect your results.")
            
        pre_res = gp.prerank(rnk=rnk, 
                        gene_sets=filtered_gene_sets,
                        threads=20,
                        min_size=min_size,
                        max_size=max_size,
                        permutation_num=permutation_num,
                        outdir=None,
                        seed=seed,
                        verbose=False,
                      )
        return pre_res
    except Exception as e:
        st.error(f"Error in GSEApy: {str(e)}")
        st.info("If using a single gene set, try reducing min_size or switching to 'fgsea (R)' method.")
        return None

def fix_term_for_filename(term):
    """ファイル名として安全な形式に変換する"""
    # スペースをドットに変換、その他の特殊文字も処理
    return re.sub(r'[^a-zA-Z0-9_]', '.', term)



@st.cache_data
# calc_fgsea関数を一時的にキャッシュから外して、クラス定義を適切に行います
# @st.cache_data アノテーションを削除して新しいバージョンを使います

def calc_fgsea(rnk, gene_sets, min_size, max_size):
    # 一時ファイルを作成
    rnk_file = f"{gsea_temp_dir}/rnk_for_fgsea.tsv"
    gmt_file = f"{gsea_temp_dir}/gene_sets_for_fgsea.gmt"
    result_file = f"{gsea_temp_dir}/fgsea_results.json"
    res_dir = f"{gsea_temp_dir}/res_data"
    
    # rankファイルを保存
    rnk.to_csv(rnk_file, sep='\t', header=False)
    
    # GMTファイルを作成
    with open(gmt_file, 'w') as f:
        for term, genes in gene_sets.items():
            joined_genes = "\t".join(genes)
            f.write(f"{term}\t{term}\t{joined_genes}\n")
    
    # RESデータを格納するディレクトリを作成
    if not os.path.exists(res_dir):
        os.makedirs(res_dir)
    
    # Rスクリプトを作成（RESデータの計算も含む）
    r_script = f'''
    library(fgsea)
    library(data.table)

    # 遺伝子ランクを読み込む
    ranks <- read.table("{rnk_file}", header=FALSE, sep="\t")
    ranks <- setNames(ranks$V2, ranks$V1)

    # 遺伝子セットを読み込む
    pathways <- gmtPathways("{gmt_file}")

    # fgseaを実行
    fgseaRes <- fgsea(pathways=pathways, 
                       stats=ranks,
                       minSize={min_size},
                       maxSize={max_size})

    # 安全なファイル名を生成する関数
    makeSafeFilename <- function(name) {{
      # スペースや特殊文字をドットに置き換え
      safe_name <- gsub("[^a-zA-Z0-9_]", ".", name)
      return(safe_name)
    }}

    # RESデータを計算して保存する関数
    calculateAndSaveRES <- function(pathway, ranks, pathways, output_dir) {{
        # 遺伝子セットを取得
        pathway_genes <- pathways[[pathway]]
        
        # ランク付けされた遺伝子リスト
        gene_list <- names(ranks)
        gene_ranks <- ranks
        
        # ランキングに基づいてソート
        sorted_idx <- order(gene_ranks, decreasing=TRUE)
        gene_list <- gene_list[sorted_idx]
        gene_ranks <- gene_ranks[sorted_idx]
        
        # ヒットの配列を作成（遺伝子セットに含まれる場合は1、そうでない場合は0）
        N <- length(gene_list)
        hits <- integer(N)
        for (i in 1:N) {{
            if (gene_list[i] %in% pathway_genes) {{
                hits[i] <- 1
            }}
        }}
        
        # 遺伝子セット内の遺伝子数
        NH <- sum(hits)
        
        # RESを計算
        RES <- numeric(N)
        running_sum <- 0
        
        for (i in 1:N) {{
            if (hits[i] == 1) {{
                # 遺伝子セット内の遺伝子
                running_sum <- running_sum + 1/NH
            }} else {{
                # 遺伝子セット外の遺伝子
                running_sum <- running_sum - 1/(N-NH)
            }}
            RES[i] <- running_sum
        }}
        
        # 安全なファイル名を生成
        safe_pathway <- makeSafeFilename(pathway)
        
        # 結果を保存
        res_file <- file.path(output_dir, paste0(safe_pathway, "_res.tsv"))
        write.table(RES, file=res_file, sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        
        # ヒット位置も保存
        hits_file <- file.path(output_dir, paste0(safe_pathway, "_hits.tsv"))
        hit_indices <- which(hits == 1)
        write.table(hit_indices, file=hits_file, sep="\\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
        
        # デバッグ用にファイル名を出力
        cat("Saved RES and hits for", pathway, "as", res_file, "and", hits_file, "\\n")
        
        return(RES)
    }}

    # 各パスウェイのRESを計算して保存
    if (!dir.exists("{res_dir}")) {{
        dir.create("{res_dir}", recursive=TRUE)
    }}

    for (pathway in names(pathways)) {{
        if (pathway %in% fgseaRes$pathway) {{
            tryCatch({{
                calculateAndSaveRES(pathway, ranks, pathways, "{res_dir}")
            }}, error = function(e) {{
                cat("Error calculating RES for", pathway, ":", conditionMessage(e), "\\n")
            }})
        }}
    }}

    # Leading edgeを文字列に変換
    fgseaRes$leadingEdge <- lapply(fgseaRes$leadingEdge, function(x) paste(x, collapse=","))

    # データフレームを行ごとにリストに変換
    result_list <- split(fgseaRes, seq(nrow(fgseaRes)))

    # pathwayを名前として使用
    names(result_list) <- fgseaRes$pathway

    # JSONとして保存
    library(jsonlite)
    write_json(result_list, "{result_file}", pretty=TRUE)
    '''
    
    r_script_file = f"{gsea_temp_dir}/run_fgsea.R"
    with open(r_script_file, 'w') as f:
        f.write(r_script)
    
    # Rスクリプトを実行
    try:
        result = subprocess.run(['Rscript', r_script_file], check=True, capture_output=True, text=True)
        st.write("R script executed successfully")
    except subprocess.CalledProcessError as e:
        st.error(f"Error executing R script: {e}")
        st.error(f"stderr: {e.stderr}")
        return None
    
    # fgseaの結果を読み込む
    try:
        with open(result_file, 'r') as f:
            fgsea_json = json.load(f)
        
        st.write("Successfully loaded fgsea results")
        
        if not fgsea_json:
            st.error("No results found in fgsea output")
            return None
            
    except Exception as e:
        st.error(f"Error reading fgsea results: {e}")
        return None
    
    # RESデータを読み込んでfgsea_jsonに追加
    for term in fgsea_json:
        # ファイル名のために用語を安全な形式に変換
        safe_term = fix_term_for_filename(term)
        res_file = f"{res_dir}/{safe_term}_res.tsv"
        hits_file = f"{res_dir}/{safe_term}_hits.tsv"
        
        try:
            if os.path.exists(res_file):
                # RESデータを読み込む
                res_data = np.loadtxt(res_file)
                # データがリストなら辞書に変換
                if isinstance(fgsea_json[term], list) and len(fgsea_json[term]) > 0:
                    # リストの最初の要素を取得し辞書に変換
                    data_item = fgsea_json[term][0]
                    fgsea_json[term] = data_item
                
                # fgsea_json[term]が辞書でない場合は辞書に変換
                if not isinstance(fgsea_json[term], dict):
                    fgsea_json[term] = {'ES': 0, 'NES': 0, 'pval': 1, 'padj': 1, 'leadingEdge': ''}
                
                # RESデータを追加
                fgsea_json[term]['res_data'] = res_data.tolist()
            
            if os.path.exists(hits_file):
                # ヒット位置を読み込む
                hit_indices = np.loadtxt(hits_file, dtype=int)
                
                # fgsea_json[term]が辞書でない場合は対応
                if not isinstance(fgsea_json[term], dict):
                    fgsea_json[term] = {'ES': 0, 'NES': 0, 'pval': 1, 'padj': 1, 'leadingEdge': ''}
                
                # ヒット位置を追加
                fgsea_json[term]['hit_indices'] = hit_indices.tolist()
                
        except Exception as e:
            st.write(f"Warning: Could not load RES data for {term}: {e}")
    
    return FgseaResults(fgsea_json, rnk, gene_sets)


# 正しく定義されたFgseaResultsクラス
class FgseaResults:
    def __init__(self, fgsea_json, rnk_data, gene_sets=None):
        self.results = {}
        self.ranking = rnk_data
        self.gene_sets = gene_sets if gene_sets is not None else {}
        
        for term, data in fgsea_json.items():
            # データを抽出
            if isinstance(data, list) and len(data) > 0:
                data_item = data[0]
            elif isinstance(data, dict):
                data_item = data
            else:
                continue
                
            # 結果エントリを作成
            self.results[term] = {
                'es': data_item.get('ES', 0),
                'nes': data_item.get('NES', 0),
                'pval': data_item.get('pval', 1),
                'fdr': data_item.get('padj', 1),
                'fwerp': data_item.get('pval', 1),
                'tag %': 0,
                'gene %': 0,
                'lead_genes': data_item.get('leadingEdge', ''),
                # RESとヒット位置
                'res_data': data_item.get('res_data', []),
                'hit_indices': data_item.get('hit_indices', [])
            }
    
    def plot(self, terms, **kwargs):
        """
        GSEApyのgseaplot関数を使用してプロット
        """
        # gseaplotをインポート
        from gseapy.plot import gseaplot
        
        if isinstance(terms, str):
            terms = [terms]
        
        fig, axes = plt.subplots(len(terms), 1, figsize=kwargs.get('figsize', (6, 6*len(terms))))
        if len(terms) == 1:
            axes = [axes]
        
        for i, term in enumerate(terms):
            if term in self.results:
                result = self.results[term]
                
                # ランクメトリックを取得
                if hasattr(self.ranking, 'iloc'):
                    rank_metric = self.ranking.iloc[:, 0].values
                elif hasattr(self.ranking, 'values'):
                    rank_metric = self.ranking.values
                else:
                    rank_metric = np.array(self.ranking).flatten()
                
                # ヒット位置を取得
                hits = []
                if 'hit_indices' in result and result['hit_indices']:
                    hits = result['hit_indices']
                    # 数値型に変換
                    hits = [int(x) for x in hits]
                
                # RESデータを取得
                res = []
                if 'res_data' in result and result['res_data']:
                    res = result['res_data']
                    # 数値型に変換
                    res = np.array(res, dtype=float)
                
                # データの長さの一貫性をチェック
                if len(res) != len(rank_metric):
                    print(f"Warning: RES length ({len(res)}) doesn't match rank_metric length ({len(rank_metric)})")
                    # 長さを合わせる
                    if len(res) > 0:
                        from scipy.interpolate import interp1d
                        x_old = np.linspace(0, 1, len(res))
                        x_new = np.linspace(0, 1, len(rank_metric))
                        f = interp1d(x_old, res, bounds_error=False, fill_value=(res[0], res[-1]))
                        res = f(x_new)
                
                # GSEApyのgseaplot関数を使用
                try:
                    # 統計値はfgseaから取得
                    gseaplot(
                        rank_metric=rank_metric,    # ランキングメトリック
                        term=term,                  # 遺伝子セット名
                        hits=hits,                  # ヒット位置 
                        nes=result['nes'],          # NES
                        pval=result['pval'],        # p値
                        fdr=result['fdr'],          # FDR
                        RES=res,                    # RESデータ
                        ax=axes[i],                 # プロット先の軸
                        ofname=None                 # ファイル出力なし
                    )
                except Exception as e:
                    print(f"Error with gseaplot: {e}")
                    axes[i].text(0.5, 0.5, f"Error plotting: {str(e)}", 
                               ha='center', va='center', transform=axes[i].transAxes)
            else:
                axes[i].text(0.5, 0.5, f"Term '{term}' not found in results", 
                           ha='center', va='center', transform=axes[i].transAxes)
        
        plt.tight_layout()
        return fig
    

    
    def plot_with_r(self, term, output_file, figsize=(6, 6), aspect_ratio=0.8):
        """
        Generate an enrichment plot for a specific term using R's fgsea package
        
        Parameters:
        term (str): The gene set term to plot
        output_file (str): Path to save the PDF output
        figsize (tuple): Figure size in inches (width, height)
        aspect_ratio (float): Height to width ratio for the plot (default: 0.6)
        
        Returns:
        bool: True if successful, False otherwise
        """
        if term not in self.results:
            print(f"Term '{term}' not found in results")
            return False
        
        # Get result data for this term
        result = self.results[term]
        
        # Create a temporary directory for R script and data
        temp_dir = os.path.dirname(output_file)
        r_script_file = os.path.join(temp_dir, f"plot_script_{fix_term_for_filename(term)}.R")
        rnk_file = os.path.join(temp_dir, f"ranking_{fix_term_for_filename(term)}.tsv")
        
        # Save ranking data to a file
        if hasattr(self.ranking, 'to_csv'):
            self.ranking.to_csv(rnk_file, sep='\t', header=False)
        else:
            with open(rnk_file, 'w') as f:
                for gene, value in zip(self.ranking.index, self.ranking.iloc[:, 0]):
                    f.write(f"{gene}\t{value}\n")
        
        # Get lead genes (if available)
        lead_genes = result.get('lead_genes', '')
        if lead_genes and isinstance(lead_genes, str):
            lead_genes_list = lead_genes.split(',')
        else:
            # If lead genes are not available, use an empty list
            lead_genes_list = []
        
        # Adjust dimensions based on aspect ratio
        pdf_width = figsize[0]
        pdf_height = figsize[0] * aspect_ratio  # Set height based on width
        
        # Create R script for plotting
        r_script = f"""
        library(fgsea)
        library(ggplot2)
        
        # Set output dimensions
        pdf("{output_file}", width={pdf_width}, height={pdf_height})
        
        # Load ranking data
        ranks <- read.table("{rnk_file}", header=FALSE, sep="\t")
        names(ranks) <- c("gene", "rank")
        ranks <- setNames(ranks$rank, ranks$gene)
        
        # Get pathway genes from the gene_sets if available, otherwise use lead genes
        pathways <- list()
        """
        
        # Use gene_sets if available, otherwise use lead genes
        if term in self.gene_sets:
            genes_str = ','.join([f'"{gene}"' for gene in self.gene_sets[term]])
            r_script += f'pathways[["{term}"]] <- c({genes_str})\n'
        else:
            genes_str = ','.join([f'"{gene}"' for gene in lead_genes_list])
            if genes_str:
                r_script += f'pathways[["{term}"]] <- c({genes_str})\n'
            else:
                r_script += f'pathways[["{term}"]] <- c()\n'
        
        r_script += f"""
        # Get statistics from pre-computed results
        nes <- {result['nes']}
        pval <- {result['pval']}
        padj <- {result['fdr']}
        
        # Create basic plot with statistics
        p <- plotEnrichment(pathways[["{term}"]], ranks) +
          labs(title="{term}",
               subtitle=paste("NES =", round(nes, 2), 
                             "  P-value =", signif(pval, 3),
                             "  FDR =", signif(padj, 3))) +
          theme_minimal()
        
        # Print the plot to PDF
        print(p)
        dev.off()
        
        # Create PNG for Streamlit display
        png_file <- gsub("\\\\.pdf$", ".png", "{output_file}")
        png(png_file, width={pdf_width*100}, height={pdf_height*100}, res=100)
        print(p)
        dev.off()
        """
        
        # Write R script to file
        with open(r_script_file, 'w') as f:
            f.write(r_script)
        
        # Run R script
        try:
            result = subprocess.run(['Rscript', r_script_file], check=True, capture_output=True, text=True)
            print(f"Successfully generated plot for {term}")
            return True
        except subprocess.CalledProcessError as e:
            print(f"Error executing R script for {term}: {e}")
            print(f"R stderr: {e.stderr}")
            return False



def add_GO_term(i):
    new_cont = [i, i]
    new_cont.extend(GO[i])
    return new_cont


@st.cache_data
def calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch):
    orig_len = len(df)
    df = df[np.isfinite(df[P_column]) & pd.notnull(df[P_column])]     # FCやpがNAのものを除く
    df = df[np.isfinite(df[FC_column]) & pd.notnull(df[FC_column])] 
    if len(df) < orig_len:
        st.warning("The P or FC columns contain inf or NA")
    # p=0がないか、みる
    inv_parameter = 1
    if inv_switch:
        inv_parameter = -1
    p_0 = (df.loc[:,P_column] == 0)
    if not any(p_0):
        #scoreを作る
        if rank_metric == 'sign(LFC) x -log10(P)':
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
        else:
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
     # p=0があるとき
    else:
        st.write("p=0 data are:")
        st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])
        # 0e0がとして読まれる　LogFCも0のはず
        if any((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)):
            st.write("And FC=0. Probably original 0.00E+00 means 1.")
            st.write(df.loc[((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)), [Gene_column,FC_column,P_column]])
            st.write("Convert 0 to 1.")
            df.loc[((df.loc[:,FC_column] == 0) & (df.loc[:,P_column] == 0)), [FC_column,P_column]] = [1,1]
            p_0 = (df.loc[:,P_column] == 0) # FC>0の0
            if any(p_0):
                st.write("Remaining p=0 data are:")
                st.write(df.loc[(df.loc[:,P_column] == 0), (Gene_column, FC_column, P_column)])

        if rank_metric == 'sign(LFC) x -log10(P)':
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * np.sign(x[FC_column]) * inv_parameter, axis =1)
        else:
            df.loc[:, 'score'] = df.apply(lambda x: -1 * np.log10(x[P_column]) * x[FC_column] * inv_parameter, axis =1)
        #Seurat "MAST"だと318あたり？
        if input_file_type == 'Seurat':
            #max_score = np.log10(1e-324) # 1e-324 == 0でTRUEになる log10を計算するとinf
            max_score = -324
            st.write("\nMax score: "+str(max_score))
        else:
            #max_score = np.log10(1e-324) # 1e-324 == 0でTRUEになる pythonでも同じ　1e-324 + 1e-323でも計算される
            max_score = -324
            st.write("\nMax score: "+str(max_score))
        # 順位付けのためにFCの値を足す
        df.loc[(p_0 & (df.loc[:,FC_column]>0)),'score'] = max_score * -1 + df.loc[:,FC_column] * inv_parameter #条件を括弧で囲むこと！！！
        df.loc[(p_0 & (df.loc[:,FC_column]<0)),'score'] = max_score + df.loc[:,FC_column] * inv_parameter
        st.write('Ranking score are -log10(P-values)')
    return df['score'].to_frame() #DFに変換してから返す

st.sidebar.title("Options")
# --- Initialising SessionState ---
if "gseacalc" not in st.session_state:
      st.session_state.gseacalc = False
if "gsea" not in st.session_state:
      st.session_state.gsea = False

if 'filename_add' not in globals(): #最初からやり直しになるときに以前のデータを保持
 #    st.write('file name kept')
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

if "rnk" not in st.session_state:
    st.session_state.rnk = None

if "rnk_name" not in st.session_state:
    st.session_state.rnk_name = None

if "pre_res" not in st.session_state:
    st.session_state.pre_res = None

st.markdown("### GSEA and overrepresentation test by GSEApy")

Analysis_mode = st.radio(
    "Analysis mode:",
    ('Prerank','Over-representation', 'Gene set score'), key='Prerank')
st.markdown("Prerank: for GSEA")
st.markdown("Over-representation for a group of genes (e.g., DEG)")
st.markdown("Gene set score for calculation of gene set score in each sample using ssGSEA or GSVA")

if Analysis_mode == "Over-representation":
    st.markdown("### Over-representation mode")
    Test_mode = st.radio("Test mode:", ('Hypergeometric test','Enrichr web'), key='Hypergeometric test')
    st.markdown("##### Genes (comma, space, CR separated):")
    genes = st.text_input("genes",label_visibility = 'collapsed')
    gene_list = []
    gsea_dir = gsea_temp_dir + "/ORA"  # dir作成
    if len(genes) > 0:
        genes = genes.replace("'","")
        genes = genes.replace('"',"")
        gene_list = genes.split(' ') #まず空白で分離
        gene_list = list(filter(lambda a: a != '', gene_list)) #空白のみを除く
        if ',' in genes:
            gene_list = sum([x.split(',') for x in gene_list],[]) #sumで平坦化 sum(x, [])
        if '\t' in genes:
            gene_list = sum([x.split('\t') for x in gene_list],[])
        if '\n' in genes:
            gene_list = sum([x.split('\n') for x in gene_list],[])
        gene_list =sorted(set(gene_list), key=gene_list.index)
        st.write(gene_list[:3])

    if Test_mode == 'Enrichr web':
        test_name = 'Enrichr_Web'
        GO_list = get_library_name()

        GO_name = st.selectbox('Select gene set', GO_list)

        if st.button('Run Enrichr web analysis'):

            enr = calc_enrichr(gene_list=gene_list, gene_sets=GO_name)

            st.dataframe(enr.results)

    else:
        test_name = 'GSEApy_Hypergeometric'
        species = st.radio("Species:", ('mouse','human'))
        db = st.radio("DB:", ('mSigDB','Enrichr', 'Homemade', 'Your own GMT file'))

        if db == 'mSigDB':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/mSigDB_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/mSigDB"
        elif db == 'Enrichr':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/enrichr_database_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/enrichr_database"
        elif db == 'Homemade':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt"
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
            else:
                st.stop()

        if db != "Your own GMT file":
            files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
            files_file.sort()
            key_index = 0
            if db == 'mSigDB':
                key_index=len(files_file) - 1
            GO_name = st.multiselect('Select gene set',files_file, default = files_file[key_index])
#            GO_name = st.selectbox('Select gene set',files_file, index=key_index)
            if db == 'mSigDB' or db == 'Homemade': #GMT fileからdictへ変換
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
        st.markdown("By default, all genes in the gene sets are used as background. However, all genes in the DEG analysis are a better background. To do this, define background genes.")
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
                    df.columns = df.iloc[0,:].tolist()    # transposeすると狂うので、transposeした後にcolumnsを決める
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
                #        df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
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

            if not os.path.exists(gsea_dir):
                st.write("GSEA_dir does not exists")
                if not os.path.exists('temp'):
                    os.mkdir('temp')
                os.mkdir(gsea_dir)


            enr.results = enr.results.sort_values('Adjusted P-value', ascending=True)

            st.dataframe(enr.results)
            p_thre = 0.05
            num_bar = 10
            py_x_size = 600
            py_y_size = 400
            with st.sidebar:
                p_thre = st.number_input('Visuallzation threshold for adj. P', min_value =0.0, step = 0.01, value=0.05)
                num_bar = int(st.number_input('Number of terms to visualize', min_value =1, step = 1, value=10))
                ORA_bar_cmap = st.selectbox('Barplot color map:', ('viridis_r', 'Blues', 'Greys', 'BrBG', 'BuGn', 'BuPu','GnBu', 'Greens','OrRd', 'Oranges', 'autumn', 'binary', 'bone',  'gist_gray', 'magma'), index = 0)
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

#                fig = px.bar(sig_sub, y="Term", x="-log10adjP", color ='Combined Score', orientation='h',  width=py_x_size, height=py_y_size)
#                st.plotly_chart(fig,use_container_width=True)

                pal = sns.color_palette(ORA_bar_cmap, as_cmap=True)
                fig, ax = plt.subplots(figsize=(py_x_size/72, py_y_size/72))
                ax =  sns.barplot(data= sig_sub,  y='Term', x="-log10adjP", palette=pal(sig_sub["Combined Score"]/max(sig_sub["Combined Score"]) ))
                ax.set_ylabel('')
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.yaxis.set_ticks_position('none')
                ax.invert_yaxis()
                sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,max(sig_sub["Combined Score"])))
                sm.set_array([])


                cax = plt.axes([1,0.2, 0.04, 0.6])
                cbar  = plt.colorbar(sm,  cax=cax)
                cbar.set_label('Combined score', rotation=270,labelpad=25)
#                cax.invert_yaxis()
                st.pyplot(fig)

#            csv = convert_df(enr.results)
                if st.button('Prepare files to download'):
                    GO_name = "_".join(GO_name).replace('.gmt','')
                    enr.results.to_csv(gsea_dir + "/" + GO_name + "." + test_name + '.tsv', sep = '\t'  )
                    fig.savefig(gsea_dir +  "/" + GO_name + "." + test_name + '.pdf',  format='pdf', bbox_inches='tight')
                    down_name = 'ORA_' + ''.join(GO_name)
                    shutil.make_archive(gsea_temp_dir + "/" + down_name, format='zip',root_dir= gsea_dir)

                    with open(gsea_temp_dir + "/" + down_name + '.zip', "rb") as fp:
                        btn = st.download_button(
                            label="Download Results",
                        data=fp,
                        file_name=down_name + ".zip",
                        mime = "zip"
                        )


# prerank
elif Analysis_mode == "Prerank":
    st.markdown("### Prerank-GSEA mode")

    # fgseaオプションを追加
    prerank_method = st.radio(
        "GSEA method:",
        ('GSEApy', 'fgsea (R)'), key='prerank_method')

    rnk_type = st.radio("Data format:",
    ('Differential expression analysis (DEA) file', 'rank'), index = 0)

    if rnk_type != 'rank':
        input_file_type = st.radio(
             "Data format of DEA file:",
             ('tsv','csv', 'excel', 'Seurat'), index = 0)
        st.write("Seurat: Seurat/SCALA's marker gene analysis output file")

        uploaded_file = st.file_uploader("Upload rnk file", type=['tsv','csv','txt','xls','lxsx'])
        if uploaded_file is not None:
            if input_file_type == "csv":
                df = read_csv(uploaded_file, index_col = 0)
            elif input_file_type == 'tsv':
                df = read_csv(uploaded_file, sep = '\t', index_col = 0)
            elif input_file_type == 'Seurat':
                # Read the Seurat format file
                df_orig = read_csv(uploaded_file, sep = '\t')
                
                # Check if we have required columns
                required_cols = ['p_val', 'avg_log2FC', 'gene', 'cluster']
                missing_cols = [col for col in required_cols if col not in df_orig.columns]
                if missing_cols:
                    st.error(f"Missing required columns in Seurat file: {', '.join(missing_cols)}")
                    st.stop()
                
                # Get the list of unique clusters
                clusters = df_orig['cluster'].unique()
                
                # Ask the user to select a cluster
                selected_cluster = st.selectbox('Select cluster to analyze:', clusters)
                
                # Filter data for the selected cluster
                df_cluster = df_orig[df_orig['cluster'] == selected_cluster].copy()
                
                # Convert the Seurat format to the expected format for GSEA
                # Set gene as index
                df = df_cluster.set_index('gene')
                
                # Rename columns to include cluster name for clarity
                df = df.rename(columns={
                    'p_val': f'{selected_cluster}_p_val',
                    'avg_log2FC': f'{selected_cluster}_avg_log2FC',
                    'p_val_adj': f'{selected_cluster}_p_val_adj',
                    'pct.1': f'{selected_cluster}_pct.1',
                    'pct.2': f'{selected_cluster}_pct.2'
                })
                
                # Drop the cluster column as it's redundant now
                if 'cluster' in df.columns:
                    df = df.drop('cluster', axis=1)
                
                st.write(f"Analyzing cluster: {selected_cluster}")
                st.write(f"Number of genes: {len(df)}")
            else:
                df = read_xl(uploaded_file, index_col = 0)

            down_file_name = os.path.basename(uploaded_file.name)

            df.iloc[0:3,:]
            rank_metric = st.radio(
                "Ranking metric:",
                ('sign(LFC) x -log10(P)', 'LFC x -log10(p)', "DESeq2 stat"), index = 0)
                # calculate stat value
            # indexをGeneカラムにコピー
            df['Gene'] = df.index
            # indexの名前を除く
            df.index.name = None
            content = df.columns.tolist()
            st.write("Select pvalue and logFC")
            p_patterns = ['p.value', 'pvalue', 'p-val', 'p val', 'p_val', 'pval', 'p_val_adj']
            pvalue = [i for i in content if any(p in i.lower() for p in p_patterns) and 'adj.pval' not in i.lower()]
            fc_patterns = ['log2fc', 'fold change', 'log2foldchange', 'coef', 'logfc', 'avg_log2fc']
            fc = [i for i in content if any(pattern in i.lower() for pattern in fc_patterns)]
            gene = [i for i in content if (i not in pvalue) and (i not in fc)]
            P_column = st.selectbox('Select P-value column', pvalue)
            stat_column = re.match(r'([^\.]+)', P_column).group(1) #名前を変更する
            # ジャロ・ウィンクラー距離法
            JW_dist = [jaro_winkler_similarity(P_column, x) for x in fc]
            try:
                FC_column = st.selectbox(
                    'Select FC column',
                    fc, index = JW_dist.index(max(JW_dist)))
            except:
                FC_column = st.selectbox(
                    'Select FC column', fc)

            if "Gene" in content:
                Gene_column =  "Gene"
            elif "Symbol" in content:
                Gene_column =  "Symbol"
            else:
                Gene_column =  st.selectbox(
                'Select gene column',
                gene)

            inv_switch = st.checkbox('Invert the sign')


            df_sub = df[[P_column, FC_column]]
            rnk = calc_rank(df, P_column, FC_column, rank_metric, Gene_column, inv_switch)
            st.write(rnk.head())
            rnk_name = P_column
            st.session_state.rnk = rnk
            st.session_state.rnk_name = rnk_name

    else:# rank file

        uploaded_rnk = st.file_uploader("Upload rnk file", type=['rnk','rank','txt'])
        if uploaded_rnk is not None:
            rnk_name = uploaded_rnk.name
            rnk = read_csv(uploaded_rnk, index_col=0, sep='\t', header = None)
            st.session_state.rnk = rnk
            st.session_state.rnk_name = rnk_name

    if st.session_state.rnk is not None and st.session_state.rnk_name is not None:
        rnk = st.session_state.rnk
        rnk_name = st.session_state.rnk_name
        species = st.radio("Species:", ('mouse','human'), index = check_species_index(rnk.index.to_list()[:50]))
        db = st.radio("DB:", ('mSigDB','Enrichr', 'Homemade', 'Your own GMT file'))

        if db == 'mSigDB':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/mSigDB_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/mSigDB"
        elif db == 'Enrichr':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/enrichr_database_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/enrichr_database"
        elif db == 'Homemade':
            if species == 'mouse':
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt_mouse"
            else:
                dir_path = "/home/cellxgene/streamlit/db/custum_gmt"

        else:
            uploaded_gmt = st.file_uploader("Upload GMT file", type=['txt','gmt'])
            if uploaded_gmt is not None:
                GO_name = [uploaded_gmt.name]
                stringio = StringIO(uploaded_gmt.getvalue().decode("utf-8"))
                s = stringio.read()
                t = s.split('\n')
                gmt =[x.split('\t') for x in t]
                GO = dict()
                for i in gmt:
                    GO[i[0]] = i[2:]
            else:
                st.stop()



        if db != "Your own GMT file":
            files_file = [f for f in os.listdir(dir_path) if os.path.isfile(os.path.join(dir_path, f))]
            files_file.sort()
            key_index = 0
            if db == 'mSigDB':
                key_index=len(files_file) - 1

            GO_name = st.multiselect('Select gene set',files_file, default = files_file[key_index])

            if db == 'mSigDB' or db == 'Homemade': #GMT fileからdictへ変換
                GO = dict()
                for i in GO_name:
                    GO_file = dir_path + '/' + i
                    with open(GO_file) as f:
                        s = f.read()
                    t = s.split('\n')
                    gmt =[x.split('\t') for x in t]
                    GO_dic = dict()
                    for i in gmt:
                        if len(i) >= 3 and i[0].strip():  # キーが空でなく、少なくとも3要素あることを確認
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
            min_size = int(st.number_input('Minimum number of genes in a gene set', min_value =1, step =1, value=15))
            max_size = int(st.number_input('Maximum number of genes in a gene set', min_value =1, step =100, value=500))
            
            # GSEApyの場合のみpermutation_numを表示
            if prerank_method == 'GSEApy':
                permutation_num = int(st.number_input('Number of permutation', min_value =100, step =100, value=1000))
                seed = int(st.number_input('Seed for permutation (0: time stamp)', min_value =0, step =1, value=6))
                if seed == 0:
                    import time
                    seed = int(time.time())
            
            rev_rnk = st.checkbox("Reverse rank order?", value = False)
            # 修正後のコード:
            if rev_rnk:
                # DataFrameの構造を確認して適切に処理
                if isinstance(rnk, pd.DataFrame):
                    # DataFrameの最初の列（スコア値）を反転
                    col_name = rnk.columns[0]
                    rnk[col_name] = -rnk[col_name]
                    st.write(f"Reversed values in column: {col_name}")
                else:
                    # DataFrameでない場合（念のため）
                    st.error("Error: rank data is not in expected format")
                    
                rev_rnk_name = st.text_input("New rank name:", rnk_name + "_rev")
                rnk_name = rev_rnk_name
                
                # 反転後のデータを表示
                st.write("Preview of reversed ranking:")
                st.write(rnk.head())

      #  st.write(GO_name)
        GO_name_dir = "_".join(GO_name).replace('.gmt','')
        rnk_name_dir = rnk_name.replace(".rnk",'')
        rnk_name_dir = rnk_name_dir.replace(".rank",'')


        gsea_dir = gsea_temp_dir + "/" + rnk_name_dir + "_" + GO_name_dir # dir作成
     #   st.write(gsea_dir)
    #    if os.path.exists(gsea_dir) and not st.session_state.gsea:
    #        shutil.rmtree(gsea_dir)
    #        st.write("GSEA_dir exists")
        if not os.path.exists(gsea_dir):
            st.write("GSEA_dir does not exists")
            if not os.path.exists('temp'):
                os.mkdir('temp')
            os.mkdir(gsea_dir)
            os.mkdir(gsea_dir + '/upregulated_enrichment')
            os.mkdir(gsea_dir + '/downregulated_enrichment')


        if st.button('Run prerank GSEA test', type = 'primary') and GO: # 状態にかかわらず、ボタンが押されたら計算する。
            st.session_state.pre_res = None

            if list(rnk.index.duplicated()).count(True) > 0:
                st.markdown("#### There are duplicated genes in the rank file.")
                st.write('Dupliated genes:' +  ', '.join(list(rnk[rnk.index.duplicated()].index)))
                st.write("The instance of the largest absolute value will be kept.")
                st.markdown("---")
                # absoluteで最大値を取る
                rnk['abs'] = np.abs(rnk.iloc[:,0])
                rnk = rnk.sort_values('abs', ascending=False)
                rnk = rnk[~rnk.index.duplicated(keep='first')]
                rnk.drop('abs', axis =1, inplace=True)
                st.write('Updated rank')
                st.write(rnk.head())

            # 選択されたメソッドに応じて処理を分岐
            if prerank_method == 'GSEApy':
                pre_res = calc_prerank(rnk=rnk, gene_sets=GO, min_size=min_size, max_size=max_size,
                                permutation_num=permutation_num, seed=seed)
            else:  # fgsea (R)
                pre_res = calc_fgsea(rnk=rnk, gene_sets=GO, min_size=min_size, max_size=max_size)
                if pre_res is None:
                    st.error("fgsea analysis failed. Make sure R and the fgsea package are installed correctly.")
                    st.stop()
                # ここでNoneチェックを追加
            if pre_res is None:
                st.error("Analysis failed. Please check your inputs or try a different method.")
                st.stop()  # 処理を中断して以降のコードを実行しない

            st.session_state.pre_res = pre_res
            terms = list(pre_res.results.keys())
            gsea_res = pd.DataFrame(columns = ['SIZE','ES','NES','NOM p-val','FDR q-val','FWER p-val','TAG %','GENE %', 'LEADING GENES'])
            gsea_res.index.name = "NAME"
            for i in terms:
                cont = pre_res.results[i]
                gsea_res.loc[i] = [len(GO[i]), cont['es'],cont['nes'],cont['pval'],cont['fdr'],
                cont['fwerp'],cont['tag %'],cont['gene %'],cont['lead_genes']]
            gsea_res = gsea_res.sort_values('FDR q-val', ascending=True)
        #    gsea_res.to_csv('temp.tsv', sep = '\t')
            up_res = gsea_res.loc[gsea_res['ES'] > 0]
            up_res =  up_res.sort_values('FDR q-val', ascending=True)
            down_res = gsea_res.loc[gsea_res['ES'] < 0]
            down_res =  down_res.sort_values('FDR q-val', ascending=True)
            st.session_state.gsea_res = gsea_res
            st.session_state.up_res = up_res
            st.session_state.down_res = down_res  #必要なデータは状態変数に入れておく
            st.session_state.terms = terms

        if st.session_state.pre_res is not None: # 一度計算した後は、ここからスタート
            st.markdown("##### You must run the analysis again if the parameters are changed.")
            st.markdown("---")
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

                # fgseaメソッドが選択されている場合のみアスペクト比のスライダーを表示
                if prerank_method == 'fgsea (R)':
                    aspect_ratio = st.slider(
                        "Plot aspect ratio (height/width):", 
                        min_value=0.3, 
                        max_value=1.0, 
                        value=0.6, 
                        step=0.1,
                        help="Adjust this to make the plot wider (lower values) or taller (higher values)"
                    )
                else:
                    aspect_ratio = 1.0  # デフォルト値

            if enrichment_term:
                # Check if we're using fgsea (different plotting approach needed)
                if prerank_method == 'fgsea (R)':
                    # For fgsea results, use R to create the plot
                    output_file = gsea_dir + '/' + file_name_check(enrichment_term) + '.pdf'
                    success = pre_res.plot_with_r(enrichment_term, output_file, 
                                                  figsize=(es_x_size, es_y_size), 
                                                  aspect_ratio=aspect_ratio)
                    
                    # Display the plot in Streamlit
                    png_file = output_file.replace('.pdf', '.png')
                    if os.path.exists(png_file):
                        img = plt.imread(png_file)
                        fig, ax = plt.subplots(figsize=(es_x_size, es_y_size * aspect_ratio))
                        ax.imshow(img)
                        ax.axis('off')
                        st.pyplot(fig)
                    else:
                        st.error(f"Failed to generate plot for {enrichment_term}")
                else:
                    # For GSEApy, we can use the built-in plot and gseaplot functions
                    g = pre_res.plot(terms=enrichment_term, figsize=(es_x_size, es_y_size))
                    st.pyplot(g)
                    
                    gseaplot(rank_metric=pre_res.ranking, term=enrichment_term,
                        ofname=gsea_dir + '/' + file_name_check(enrichment_term) + '.pdf',
                        **pre_res.results[enrichment_term], figsize=(es_x_size, es_y_size))
            else:
                st.markdown("#### No terms with FDR < " + str(term_thre) + " for enrichment plot.")


            res2d = pd.DataFrame(gsea_res.loc[:,["FDR q-val","LEADING GENES","NES"]])
            res2d['Term'] = res2d.index.values
            res2d.columns = ['FDR q-val', 'Lead_genes', 'NES', 'Term']

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
                if bar_type == 'Separete up/down plots':
                    Y_reverse = st.checkbox("Reverse Y order in downregulated plot?", value=False)
                #dot_size = int(st.number_input("Dot size:", value = 5, step = 1, min_value = 1))
                bar_cmap = st.selectbox('Barplot color map:', ('Default', 'Greys_r', 'Blues_r', 'BrBG_r', 'BuGn_r', 'BuPu_r','GnBu_r', 'Greens_r','OrRd_r', 'Oranges_r', 'autumn', 'binary_r', 'bone',  'gist_gray', 'magma_r', 'viridis'), index = 0)
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
                if len(up_nes) >0:
                    fig, ax = plt.subplots(figsize=(bar_x_size, bar_y_size))
                    if len(up_nes) > bar_num:
                        up_nes = up_nes.iloc[0:bar_num]
                    if bar_cmap == "Default":
                        bar_cmap_use = "OrRd_r"
                    else:
                        bar_cmap_use = bar_cmap
                    pal = sns.color_palette(bar_cmap_use, as_cmap=True)
                    ax =  sns.barplot(data= up_nes,  y='Gene set', x="NES",  palette=pal(up_nes["FDR q-val"]/bar_vis_thr  ))
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
                    fig.savefig(gsea_dir + '/' + "Barplot_up.pdf", format='pdf', bbox_inches='tight')

                if len(down_nes) >0:
                    fig, ax = plt.subplots(figsize=(bar_x_size, bar_y_size))
                    if len(down_nes) > bar_num:
                        down_nes = down_nes.iloc[0:bar_num]
                        if Y_reverse:
                            down_nes = down_nes.sort_values('NES', ascending=True)
                    if bar_cmap == "Default":
                        bar_cmap_use = "Blues_r"
                    else:
                        bar_cmap_use = bar_cmap
                    pal = sns.color_palette(bar_cmap_use, as_cmap=True)
                    ax =  sns.barplot(data= down_nes,  y='Gene set', x="NES",  palette=pal(down_nes["FDR q-val"]/bar_vis_thr  ))
                    ax.set_ylabel('')
                    ax.spines['top'].set_visible(False)
                    ax.spines['left'].set_visible(False)
                    ax.yaxis.set_ticks_position('none')
                    sm = ScalarMappable(cmap=pal, norm=plt.Normalize(0,bar_vis_thr))
                    sm.set_array([])

                    cax = plt.axes([1,0.2, 0.04, 0.6])
                    cbar  = plt.colorbar(sm,  cax=cax)
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
                    if prerank_method == 'fgsea (R)':
                        # Save plot for fgsea results using R
                        output_file = gsea_dir + '/upregulated_enrichment/' + file_name_check(enrichment_term) + '.pdf'
                        pre_res.plot_with_r(enrichment_term, output_file, 
                                            figsize=(es_x_size, es_y_size),
                                            aspect_ratio=aspect_ratio)
                    else:
                        # Save plot for GSEApy results
                        gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, 
                                ofname=gsea_dir + '/upregulated_enrichment/' + file_name_check(enrichment_term) + '.pdf',
                                **pre_res.results[enrichment_term], figsize=(es_x_size, es_y_size))
                    percent_complete = percent_complete + percent_count
                    my_bar.progress(percent_complete, text=progress_text)

                for i in range(num_down_fig):
                    enrichment_term = down_res.index.values[i]
                    if prerank_method == 'fgsea (R)':
                        # Save plot for fgsea results using R
                        output_file = gsea_dir + '/downregulated_enrichment/' + file_name_check(enrichment_term) + '.pdf'
                        pre_res.plot_with_r(enrichment_term, output_file, 
                                            figsize=(es_x_size, es_y_size),
                                            aspect_ratio=aspect_ratio)
                    else:
                        # Save plot for GSEApy results
                        gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, 
                                ofname=gsea_dir + '/downregulated_enrichment/' + file_name_check(enrichment_term) + '.pdf',
                                **pre_res.results[enrichment_term], figsize=(es_x_size, es_y_size))
                    percent_complete = percent_complete + percent_count
                    my_bar.progress(percent_complete, text=progress_text)

                if not os.path.exists(gsea_dir + '/up_png'):
                    os.mkdir(gsea_dir + '/up_png')
                if not os.path.exists(gsea_dir + '/down_png'):
                    os.mkdir(gsea_dir + '/down_png')

                for i in range(num_up_fig):
                    enrichment_term = up_res.index.values[i]
                    if prerank_method == 'fgsea (R)':
                        # Save plot for fgsea results using R
                        output_file = gsea_dir + '/up_png/' + str(i) + file_name_check(enrichment_term) + '.png'
                        pdf_file = output_file.replace('.png', '.pdf')
                        pre_res.plot_with_r(enrichment_term, pdf_file, 
                                            figsize=(es_x_size, es_y_size),
                                            aspect_ratio=aspect_ratio)
                        # The PNG is created automatically in the plot_with_r method
                    else:
                        # Save plot for GSEApy results
                        gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, 
                                ofname=gsea_dir + '/up_png/' + str(i) + file_name_check(enrichment_term) + '.png',
                                **pre_res.results[enrichment_term], figsize=(es_x_size, es_y_size))
                    percent_complete += percent_count
                    my_bar.progress(percent_complete, text=progress_text)

                for i in range(num_down_fig):
                    enrichment_term = down_res.index.values[i]
                    if prerank_method == 'fgsea (R)':
                        # Save plot for fgsea results using R
                        output_file = gsea_dir + '/down_png/' + str(i) + file_name_check(enrichment_term) + '.png'
                        pdf_file = output_file.replace('.png', '.pdf')
                        pre_res.plot_with_r(enrichment_term, pdf_file, 
                                            figsize=(es_x_size, es_y_size),
                                            aspect_ratio=aspect_ratio)
                        # The PNG is created automatically in the plot_with_r method
                    else:
                        # Save plot for GSEApy results
                        gseaplot(rank_metric=pre_res.ranking, term=enrichment_term, 
                                ofname=gsea_dir + '/down_png/' + str(i) + file_name_check(enrichment_term) + '.png',
                                **pre_res.results[enrichment_term], figsize=(es_x_size, es_y_size))
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
#                gene_sets_gmt = []
#                for i in list(GO.keys()):
#                    new_cont = [i, i]
#                    new_cont.extend(GO[i])
#                    gene_sets_gmt.append(new_cont)
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