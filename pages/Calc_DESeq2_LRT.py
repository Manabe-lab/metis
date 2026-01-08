#!!!!!!!!!!!!!! pip install rpy2==3.5.1  新しいバージョンはエラーが出る

# 基本的にglobal変数で計算する。
# pythonからassgnされるのはglobal変数


import streamlit as st
import rpy2
import csv
import re
import os
import numpy as np
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector
import pyper
import shutil
from PIL import Image
import itertools
from helper_func import clear_old_directories, clear_old_files, remove_after_space, remove_sample_num
import time
import sys
from collections import Counter

def remove_common_suffix(strings):
    if not strings or len(strings) == 0:
        return []    
    # 最も短い文字列の長さを取得
    min_length = min(len(s) for s in strings)
    # 共通の末尾部分の長さを見つける
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break            
    # 共通の末尾部分が見つからない場合は元のリストを返す
    if suffix_length == 0:
        return strings        
    # 共通の末尾部分を削除して新しいリストを作成
    return [s[:-suffix_length] for s in strings]


def rename_duplicates(df):
    """
    Rename duplicate indices by adding _2, _3, etc. to subsequent occurrences   
    Args:
        df: pandas DataFrame
    Returns:
        DataFrame with renamed indices
    """
    # Get current index values
    lis = df.index.values
    
    # Count occurrences of each value
    counts = Counter()
    new_indices = []
    
    for x in lis:
        counts[x] += 1
        if counts[x] == 1:
            new_indices.append(x)
        else:
            new_indices.append(f"{x}_{counts[x]}")
    
    # Check if there were any duplicates
    if len(lis) != len(set(lis)):
        st.markdown("#### There are duplicated rows. Converting the names...")
        st.write("The gene names of subsequent occurrences have _2, _3, etc. at the end.")
        
        # Display which names were changed
        for name, count in counts.items():
            if count > 1:
                st.write(f"'{name}' appears {count} times → {name}, " + 
                        ", ".join([f"{name}_{i}" for i in range(2, count + 1)]))
    
    # Set new index
    df.index = new_indices
    return df

@st.cache_data
def check_excel_autoconversion(dfx):
    p = re.compile(r'(\d+)\-(Mar|Sep|Oct|Dec|Feb|Nov)')
    index_name = dfx.index.values
    j = 0
    k = 0
    for i in df.index.values:
        x = p.match(i)
        if x:
            if k == 0:
                st.markdown("#### There are Excel-autoconverted gene names")
                st.write("Gene names are not converted.")
                k = 1
            autoconvert_flag = True
            st.write(i)


r = pyper.R(use_pandas=True)
f = ro.r("source('pages/deseq2_func.R')") # full pathが必要

st.set_page_config(page_title="DESeq2-LRT", page_icon="📃")

@st.cache_data
def read_csv(file, index_col=None, sep=',', header = 0):
    df_c = pd.read_csv(file, index_col = index_col, header = header, sep = sep)
    return df_c


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def calc_barplot(data, ylabel):
    fig, ax = plt.subplots()
    ax = sns.barplot(data=data)
    ax.set_xticklabels(ax.get_xticklabels(),rotation = 90)
    ax.set_ylabel(ylabel, fontsize = 14)
    return fig


# temp内に保存する
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    #古いdirecotryとファイルを削除する
    temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
    os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/res_LRT'
    st.session_state.res_dir = res_dir
    os.mkdir(res_dir)

else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/res_LRT'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)

st.sidebar.title("Options")
with st.sidebar:
    st.markdown("#### Test medhod:")
    test_method = st.radio("Test method:",
                         ["DESeq2-LRT", "limma eBayes", "edgeR-QL", "Beta Regression",
                          "Generalized Linear Model (GLM)", "Generalized Additive Model (GAM)", "maSigPro"],
                         index=0, label_visibility = 'collapsed')

    st.markdown("##### limma eBayes: Linear model + moderated t/F-test (not LRT)")
    st.markdown("##### edgeR-QL: GLM + quasi-likelihood F-test")
    st.markdown("##### limma with logit transformation and beta regression are for proportion data (0-1)")

    if test_method == 'DESeq2-LRT':
        st.markdown("---")
        st.markdown("### DESeq2 Polynomial Options:")
        add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（時間経過に伴う非線形変化）を検出するために多項式の項を加えます")

        
        if add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2次多項式：1つの方向変化（増加→減少など）を検出。3次多項式：2つの方向変化（増加→減少→増加など）を検出")
            
            # 多項式実装方法の選択
            poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()関数：直交多項式も使用可能で推奨。I()関数：単純なべき乗項の使用")
            
            # poly()関数を使う場合の追加オプション
            use_poly_function = poly_implementation.startswith("poly()")
            
            if use_poly_function:
                # 多項式タイプの選択オプション
                poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="直交多項式（Orthogonal）：共線性を避けるために推奨。Raw：解釈しやすい係数を得るが共線性の問題がある場合あり")
                use_raw = poly_type.startswith("Raw")
            else:
                use_raw = False  # I()関数使用時は関係ないのでFalse
                            
            polynomial_degree = 2 if polynomial_term.startswith('2:') else 3
            
        else:
            polynomial_variable = None
            polynomial_degree = 1
            use_raw = False
            use_poly_function = True  # デフォルト値


    if test_method == 'limma eBayes':
        limma_data = st.radio("Data type:",
            ["RNA-seq count", "Non-count data", "0-1 data (proportion, AUC etc) to logit transformation"],
            index = 1)

        if limma_data == "RNA-seq count":
            apply_logit = False
            limma_count = True
        elif limma_data == "Non-count data":
            apply_logit = False
            limma_count = False
        else:
            apply_logit = True
            limma_count = False

        # RNA-seq count特有のオプション
        if limma_count:
            st.markdown("---")
            st.markdown("### RNA-seq Count Data Options:")

            # filterByExpr option
            use_edgeR_filter = st.checkbox("Use edgeR::filterByExpr()", value=False,
                                          help="CPMベースの統計的フィルタリング（推奨）。グループ構成を考慮した適切な低発現遺伝子除外を実行します。")
            if use_edgeR_filter:
                st.markdown("##### filterByExpr parameters:")
                filter_min_count = st.number_input("min.count", value=10, min_value=0,
                                                  help="各サンプルの最小カウント数（デフォルト: 10）")
                filter_min_total_count = st.number_input("min.total.count", value=15, min_value=0,
                                                        help="全サンプル合計の最小カウント数（デフォルト: 15）")
                filter_min_prop = st.number_input("min.prop", value=0.7, min_value=0.0, max_value=1.0, step=0.1,
                                                 help="最小グループのサンプル数の割合（デフォルト: 0.7）")
                filter_large_n = st.number_input("large.n", value=10, min_value=0,
                                                help="'大きい'と見なされるグループサイズ（デフォルト: 10）")
            else:
                filter_min_count = None
                filter_min_total_count = None
                filter_min_prop = None
                filter_large_n = None

            st.markdown("---")
            st.markdown("##### Voom transformation method:")
            which_voom = st.radio(
                "Select voom method:",
                ['voom', 'voomWithQualityWeights', 'voomLmFit'],
                index=2,
                help="""**voom**: 標準的な方法。遺伝子レベルの精度重み付けのみ。

**voomWithQualityWeights**: 低品質サンプルを含むデータセット向け。サンプルレベルのquality weightsで低品質サンプルをダウンウェイト。反復計算で収束。

**voomLmFit**: スパースデータ（scRNA-seq等）向け。voomとlmFitを同時実行。効率的。sample.weightsオプションでサンプル品質重み付け可能。"""
            )

            # voomLmFit用のsample.weightsオプション
            use_sample_weights = False
            if which_voom == 'voomLmFit':
                use_sample_weights = st.checkbox(
                    "Use sample quality weights",
                    value=False,
                    help="sample.weights=TRUEを使用。サンプルレベルの品質重み付け。voomWithQualityWeightsと同様の効果を効率的に実現。"
                )

            st.markdown("---")
            st.markdown("##### eBayes options:")
            use_robust_ebayes = st.checkbox("Use robust eBayes", value=True,
                                           help="eBayes(robust=TRUE) を使用。外れ値に頑健な推定を実行（推奨）")

            st.markdown("---")
            st.markdown("##### Repeated measures / Donor correlation:")
            use_duplicate_correlation = st.checkbox(
                "Use duplicateCorrelation for repeated measures",
                value=False,
                help="""重複測定・ドナー相関の処理。同一ドナー内の複数サブタイプ測定、time-courseデータ等で使用。

例：各ドナーから複数のサブタイプ（例: Arterial, Capillary）の擬似バルクを作成した場合、同一ドナー内のサブタイプ間相関をduplicateCorrelation(block=donor)で推定し、lmFit(..., block=donor, correlation=ρ)で考慮。

これにより、サブタイプ横断の単一モデルでドナー間の生物学的変動を適切に扱えます。

注意: Blocking変数の選択は、メタデータ入力後に表示されます。"""
            )
        else:
            # Non-count or logit transformation
            use_edgeR_filter = False
            filter_min_count = None
            filter_min_total_count = None
            filter_min_prop = None
            filter_large_n = None
            which_voom = 'voom'
            use_sample_weights = False
            use_duplicate_correlation = False
            block_variable = None

            # eBayes options for non-count data
            st.markdown("---")
            st.markdown("### Options for Non-count / Logit-transformed Data:")
            st.markdown("##### eBayes options:")
            use_robust_ebayes = st.checkbox("Use robust eBayes", value=True,
                                           help="eBayes(robust=TRUE) を使用。外れ値に頑健な推定を実行（推奨）")

        st.markdown("---")
        st.markdown("### Batch Correction (SVA):")
        use_sva = st.checkbox('Use SVA (Surrogate Variable Analysis) for batch correction?', value=False,
                             help="SVAを使用して未知のバッチ効果や潜在変数を推定し、デザイン行列に追加します。")

        if use_sva:
            sva_n_sv = st.number_input(
                'Number of surrogate variables to estimate',
                min_value=1,
                max_value=10,
                value=2,
                help="推定するsurrogate variableの数。デフォルトは2。多すぎると生物学的シグナルも除去してしまう可能性があります。"
            )
        else:
            sva_n_sv = 2  # デフォルト値

        st.markdown("---")
        st.markdown("### limma Polynomial Options:")
        limma_add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（時間経過に伴う非線形変化）を検出するために多項式の項を加えます")

        if limma_add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            limma_polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2次多項式：1つの方向変化（増加→減少など）を検出。3次多項式：2つの方向変化（増加→減少→増加など）を検出")
            
            # 多項式実装方法の選択
            limma_poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()関数：直交多項式も使用可能で推奨。I()関数：単純なべき乗項の使用")
            
            # poly()関数を使う場合の追加オプション
            limma_use_poly_function = limma_poly_implementation.startswith("poly()")
            
            if limma_use_poly_function:
                # 多項式タイプの選択オプション
                limma_poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="直交多項式（Orthogonal）：共線性を避けるために推奨。Raw：解釈しやすい係数を得るが共線性の問題がある場合あり")
                limma_use_raw = limma_poly_type.startswith("Raw")
            else:
                limma_use_raw = False  # I()関数使用時は関係ないのでFalse
                            
            limma_polynomial_degree = 2 if limma_polynomial_term.startswith('2:') else 3
            
        else:
            limma_polynomial_variable = None
            limma_polynomial_degree = 1
            limma_use_raw = False
            limma_use_poly_function = True  # デフォルト値

    # edgeR-QL特有のオプション（RNA-seq count data only）
    if test_method == 'edgeR-QL':
        st.markdown("---")
        st.markdown("### edgeR-QL Options (RNA-seq count data):")

        # filterByExpr option
        use_edgeR_filter = st.checkbox("Use edgeR::filterByExpr()", value=True,
                                      help="CPMベースの統計的フィルタリング（推奨）。グループ構成を考慮した適切な低発現遺伝子除外を実行します。")
        if use_edgeR_filter:
            st.markdown("##### filterByExpr parameters:")
            filter_min_count = st.number_input("min.count", value=10, min_value=0,
                                              help="各サンプルの最小カウント数（デフォルト: 10）")
            filter_min_total_count = st.number_input("min.total.count", value=15, min_value=0,
                                                    help="全サンプル合計の最小カウント数（デフォルト: 15）")
            filter_min_prop = st.number_input("min.prop", value=0.7, min_value=0.0, max_value=1.0, step=0.1,
                                             help="最小グループのサンプル数の割合（デフォルト: 0.7）")
            filter_large_n = st.number_input("large.n", value=10, min_value=0,
                                            help="'大きい'と見なされるグループサイズ（デフォルト: 10）")
        else:
            filter_min_count = None
            filter_min_total_count = None
            filter_min_prop = None
            filter_large_n = None

        st.markdown("---")
        st.markdown("##### Quasi-likelihood test method:")
        edger_test_method = st.radio(
            "Select test method:",
            ["glmQLFTest (standard)", "glmTreat (minimum fold-change)"],
            index=0,
            help="""**glmQLFTest (standard)**: log2FC=0を検定。微小な変化も検出（GSEA用ランキングに適）。

**glmTreat (minimum fold-change)**: 最小fold-change閾値を超える変化を検定。生物学的に意味のある遺伝子に絞り込み（候補遺伝子選定・バリデーションに適）。

例: lfc=0.585設定時、log2FC>0.585（1.5倍以上変化）のみ有意判定。微小変化（log2FC=0.05等）は除外される。

一般的な基準: 0.585（1.5倍・標準的）、1.0（2倍・厳しめ）

参考: McCarthy & Chen (2012) Nucleic Acids Res"""
        )

        use_glmtreat = (edger_test_method == "glmTreat (minimum fold-change)")

        if use_glmtreat:
            treat_lfc = st.number_input(
                "Log2 fold-change threshold",
                value=0.585,
                min_value=0.0,
                step=0.1,
                help="最小log2FC閾値。0.585=1.5倍変化（標準的）、1.0=2倍変化、1.5=3倍変化。一般的な基準: 0.585（標準）、1.0（厳しめ）"
            )
        else:
            treat_lfc = None

        st.markdown("---")
        st.markdown("##### Robust estimation:")
        use_robust_edger = st.checkbox("Use robust estimation", value=True,
                                       help="estimateDisp()とglmQLFit()でrobust=TRUEを使用。外れ値に頑健な推定を実行（推奨）")

    # ベータ回帰特有のオプション
    if test_method == 'Beta Regression':
        st.markdown("### Beta Regression Options:")
        epsilon = st.number_input("Epsilon for boundary adjustment (0-1 data)", 
                                min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        add_higher = st.checkbox("Add polynomial terms?", value=False, help="非線形的な変化を捉えるため多項項を加える")
        beta_polynomial_degree = 1
        if add_higher:
            polynomial_term = st.radio("Degree", ['2:Quadratic term','3:Cubic term'], index = 0, help = "2次の項を加えるとU-shaped/inverted U-shaped patternsが捉えられる。3次の項を追加することで、より複雑な発現パターンを捉える。例えば：急上昇後に横ばいになり、その後低下するパターン、波形パターン。")
            st.markdown("#### The first item in full model will be used for polynominal term.")
            if polynomial_term == "2:Quadratic term":
                beta_polynomial_degree = 2
            else:
                beta_polynomial_degree = 3

        
    # GLM特有のオプション
    if test_method == 'Generalized Linear Model (GLM)':
        st.markdown("### GLM Options:")
        # 分布ファミリーの選択オプション
        glm_dist_family = st.radio("Probability distribution", 
                              ["Beta (0-1)", 
                               "Gaussian", 
                               "Poisson", 
                               "Negative Binomial"],
                              index=0,
                              help="データタイプに応じた確率分布を選択してください。Beta: 0-1の値（割合など）, Gaussian: 連続値データ, Poisson: シンプルなカウントデータ, Negative Binomial: 過分散のあるカウントデータ（RNA-seq, scRNA-seq等）")

        if glm_dist_family == "Beta (0-1)":
            glm_dist_short = "beta"
        elif glm_dist_family == "Gaussian":
            glm_dist_short = "gaussian"
        elif glm_dist_family == "Poisson":
            glm_dist_short = "poisson"
        elif glm_dist_family == "Negative Binomial":
            glm_dist_short = "nb"

        glm_epsilon = 0
        glm_nb_theta = None

        # リンク関数の選択
        if glm_dist_family == "Beta (0-1)":
            glm_link = st.radio("Link function for Beta distribution",
                               ["logit", "probit", "cloglog"],
                               index=0,
                               help="logit: 最も一般的、probit: 正規分布ベース、cloglog: 極値分布ベース")
            glm_epsilon = st.number_input("Epsilon for boundary adjustment", 
                                   min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        elif glm_dist_family == "Gaussian":
            glm_link = st.radio("Link function for Gaussian distribution",
                               ["identity", "log", "inverse"],
                               index=0,
                               help="identity: 線形関係、log: 正の値のみ、inverse: 逆数変換")
        elif glm_dist_family == "Poisson":
            glm_link = st.radio("Link function for Poisson distribution",
                               ["log", "identity", "sqrt"],
                               index=0,
                               help="log: 最も一般的（カウントデータ）、identity: 線形、sqrt: 平方根変換")
        elif glm_dist_family == "Negative Binomial":
            glm_link = st.radio("Link function for Negative Binomial distribution",
                               ["log", "identity", "sqrt"],
                               index=0,
                               help="log: 最も一般的（過分散カウントデータ）、identity: 線形、sqrt: 平方根変換")
            glm_nb_theta = st.number_input("過分散パラメータ (theta)", 
                                     min_value=0.1, max_value=100.0, value=10.0,
                                     help="値が小さいほど過分散が大きいことを意味します。RNA-seqでは通常5-10程度。scRNA:0.5-3 (10x:1-2)")

        # Polynomial options for GLM
        st.markdown("---")
        st.markdown("### GLM Polynomial Options:")
        glm_add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（時間経過に伴う非線形変化）を検出するために多項式の項を加えます")

        if glm_add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            glm_polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2次多項式：1つの方向変化（増加→減少など）を検出。3次多項式：2つの方向変化（増加→減少→増加など）を検出")
            
            # 多項式実装方法の選択
            glm_poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()関数：直交多項式も使用可能で推奨。I()関数：単純なべき乗項の使用")
            
            # poly()関数を使う場合の追加オプション
            glm_use_poly_function = glm_poly_implementation.startswith("poly()")
            
            if glm_use_poly_function:
                # 多項式タイプの選択オプション
                glm_poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="直交多項式（Orthogonal）：共線性を避けるために推奨。Raw：解釈しやすい係数を得るが共線性の問題がある場合あり")
                glm_use_raw = glm_poly_type.startswith("Raw")
            else:
                glm_use_raw = False  # I()関数使用時は関係ないのでFalse
                            
            glm_polynomial_degree = 2 if glm_polynomial_term.startswith('2:') else 3
            
        else:
            glm_polynomial_variable = None
            glm_polynomial_degree = 1
            glm_use_raw = False
            glm_use_poly_function = True  # デフォルト値

    # GAM特有のオプション
    if test_method == 'Generalized Additive Model (GAM)':
        st.markdown("### GAM Options:")
        # 分布ファミリーの選択オプションを追加
        dist_family = st.radio("Probability distribution", 
                              ["Beta (0-1)", 
                               "Gaussian", 
                               "Poisson", 
                               "Negative Binomial"],
                              index=0,
                              help="データタイプに応じた確率分布を選択してください。Beta: 0-1の値（割合など）, Gaussian: 連続値データ, Poisson: シンプルなカウントデータ, Negative Binomial: 過分散のあるカウントデータ（RNA-seq, scRNA-seq等）")

        if dist_family == "Beta (0-1)":
            dist_short = "beta"
        elif dist_family == "Gaussian":
            dist_short = "gaussian"
        elif dist_family == "Poisson":
            dist_short = "poisson"
        elif dist_family == "Negative Binomial":
            dist_short = "nb"

        st.write(dist_short)

        epsilon = 0
        nb_theta = None


        # 分布に応じたパラメータ設定
        if dist_family == "Beta (0-1)":
            epsilon = st.number_input("Epsilon for boundary adjustment", 
                                   min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        elif dist_family == "Negative Binomial":
            nb_theta = st.number_input("過分散パラメータ (theta)", 
                                     min_value=0.1, max_value=100.0, value=10.0,
                                     help="値が小さいほど過分散が大きいことを意味します。RNA-seqでは通常5-10程度。scRNA:0.5-3 (10x:1-2)")


        gam_k = st.slider("Spline basis dimension (k)", min_value=3, max_value=20, value=4, help = "非線形関係をモデル化するたスムージング関数の柔軟性（複雑さ）を制御するパラメータ。k の値が大きいほど、モデルはより複雑な非線形パターンを捉える。「時間点の数 + 0〜1」")
        gam_method = st.radio("Smoothing parameter estimation method", 
                            ["REML", "GCV.Cp", "ML"], index=0,
                            help="REML (Restricted Maximum Likelihood):分散パラメータと平滑化パラメータを同時に推定。バイアスが少なく、より信頼性の高い推定。小さなサンプルサイズでも比較的安定。   ML (Maximum Likelihood):最尤推定法。モデル比較（AIC, BICなど）に適す。小さなサンプルサイズではバイアスが生じる可能性。  GCV.Cp (Generalized Cross Validation / Mallows' Cp):クロスバリデーションに基づく。モデルの予測能力を最適化。予測が目的の場合に有用。")

        selected_spline = st.radio("Spline type",
                        ['Thin Plate Regression Splines (tp)', 'Cubic Regression Splines (cr)', 'Cubic Smoothing Splines (cs)'],
                        index =1, help="tp:最も汎用的なスプラインタイプ。一般的にデフォルト。  cr:3次多項式の区分的な組み合わせ。時間点が少ない場合に適しているケースがある。  cs:観測データ点にノットを配置。非常に柔軟で、データ点間を滑らかに補間。データ点が少ない場合に有用なことがある")
        # 選択に基づいて変数を設定
        if selected_spline == 'Thin Plate Regression Splines (tp)':
            spline_type = "tp"
        elif selected_spline == 'Cubic Regression Splines (cr)':
            spline_type = "cr"
        elif selected_spline == 'Cubic Smoothing Splines (cs)':
            spline_type = "cs"

    #    beta_norm = st.checkbox("Normalization by maximum value of time variable", value = False,
    #        help='最大値による正規化.時間を含むときに収束しない場合に試みる') #効果がなさそう

#        beta_normalization = "TRUE" if beta_norm else "FALSE"
        beta_normalization = "FALSE"
    if test_method in ['Beta Regression', 'Generalized Linear Model (GLM)', 'Generalized Additive Model (GAM)']:
        n_cores = st.slider("Parallel cores", min_value=1, 
                           max_value=os.cpu_count()-1, 
                           value=max(1, os.cpu_count()//2-4))


    if test_method == 'maSigPro':
        st.markdown("### maSigPro Options:")
        
        # データタイプの選択
        data_type = st.radio(
            "Data type:",
            ["RNA-seq count data (GLM)", "qPCR/continuous data (Gaussian)", "0-1 data (logit transformation)"],
            index=0
        )
        
        # データタイプに応じたパラメータ設定
        if data_type == "0-1 data (logit transformation)":
            st.markdown("##### Boundary adjustment for 0-1 data:")
            epsilon = st.number_input("Epsilon", 
                                    min_value=1e-8, 
                                    max_value=0.01, 
                                    value=1e-6,
                                    format="%.8f")
        elif data_type == "qPCR/continuous data (Gaussian)":
            st.markdown("##### qPCR/Continuous data options:")
            log_transform = st.checkbox("Log2 transform data", value=True, 
                                      help="qPCRデータの場合、通常log2変換を行います（ΔCt値など）")
            normalization = st.checkbox("Z-score normalization across samples", value=False,
                                       help="qPCRでは通常不要")
        
        # 共通のパラメータ
        degree = st.slider("Polynomial degree", min_value=1, max_value=3, value=2)
        rsq = st.number_input("R-squared cutoff", min_value=0.1, max_value=0.9, value=0.7, step=0.05)
        q_value = st.number_input("Q-value (FDR)", min_value=0.001, max_value=0.5, value=0.05, step=0.01)
        
        # クラスタリングオプション
        st.markdown("##### Clustering options:")
        cluster_method = st.radio("Clustering method", 
                                ["hclust", "kmeans", "clara"], index=0)
        k = st.slider("Number of clusters", min_value=2, max_value=15, value=9)
        
        # 可視化オプション
        st.markdown("##### Visualization:")
        plot_top_n = st.slider("Number of genes to plot", min_value=5, max_value=50, value=20)

    st.markdown("---")

st.markdown("### Time-course and ANOVA-like Differential Expression Analysis")
st.markdown("### maSigPro for time-course test")
st.markdown("""
##### Statistical Methods:
- **DESeq2-LRT**: Negative binomial GLM + Likelihood Ratio Test (RNA-seq count data)
- **edgeR-QL**: Negative binomial GLM + Quasi-likelihood F-test (RNA-seq count data)
- **limma eBayes**: Linear model + moderated t/F-test (voom-transformed counts or continuous data)
- **Beta Regression**: Beta distribution GLM (proportion data 0-1)
- **GAM**: Generalized Additive Model (flexible non-linear relationships)
""")
st.markdown("##### DESeq2-LRT, beta regression can use polynomial terms for time-course analysis")
st.markdown("##### limma uses empirical Bayes moderated statistics, NOT likelihood ratio test")
st.markdown("##### See left sidebar for options")
st.write(" ")
use_upload = 'Yes'
if 'df' in st.session_state:
    if st.session_state.df is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.df
        input_file_type = 'tsv'
        file_name_head = st.session_state.uploaded_file_name
        if "Row_name" in df.columns.to_list(): # Row_nameを含むとき
            df = df.set_index('Row_name')
            df.index.name = "Gene"



if use_upload == 'Yes':
    st.markdown("##### Data format:")
    file_type = st.radio(
        "",    ('Homer','tsv','csv','excel'), index = 1, label_visibility = 'collapsed')


    uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])
    if uploaded_file is not None:
        if file_type is not 'excel':
            if file_type == 'csv':
                df = read_csv(uploaded_file)
            else:
                df = read_csv(uploaded_file, sep = '\t')
            st.write("Original:")
            st.write(df.head())
            if file_type == 'Homer':
                df = df.iloc[:,7:]
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                # colnamesの変換
                search_word = '([^\ \(]*)\ \(.*'
                for i in range(1, len(colnames)):
                    match = re.search(search_word, colnames[i])
                    if match:
                        colnames[i] = match.group(1).replace(' ', '_')
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                try:
                    df.iloc[:,0] = df.iloc[:,0].apply(f_annotation)
                    df.columns = colnames
                except:
                    st.markdown("### File format error. Non-Homer file?")

            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames
        else: # excel
            df = read_excel(uploaded_file, index_col = 0)
            content = df.columns.tolist()
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
                # annotation/divergence以前を除く
                df = df.loc[:,'Annotation/Divergence':]
                content = df.columns.tolist()
                content[0] = 'Gene'
                df.columns = content
                st.write("Converted Annotation/Divergence to gene symbols.")
            else:
                colnames = df.columns.tolist()
                colnames[0] = 'Gene'
                df.columns = colnames

        df = df.set_index('Gene')
        file_name_head = os.path.splitext(uploaded_file.name)[0]
    else:
        sys.exit(1)

if df is not None:

############ sample名のRで使えない文字を修正
    def make_valid_r_names(names):
        """Rで有効な変数名に変換する関数"""
        valid_names = []
        changes_made = False
        
        for name in names:
            original_name = name
            
            # 1. 特殊文字を置換
            name = re.sub(r'[ ]+', '.', name)  # スペースを.に
            name = re.sub(r'[-]+', '_', name)  # ハイフンを_に
            name = re.sub(r'[^\w.]', '_', name)  # 英数字・ピリオド・アンダースコア以外を_に
            
            # 2. 先頭が数字の場合はXを付ける
            if re.match(r'^\d', name):
                name = 'X' + name
            
            # 3. 先頭がピリオドで数字が続く場合はXを付ける
            if re.match(r'^\.\d', name):
                name = 'X' + name
            
            # 4. 予約語チェック（基本的なもの）
            r_reserved = ['if', 'else', 'repeat', 'while', 'function', 'for', 'in', 'next', 'break', 
                         'TRUE', 'FALSE', 'NULL', 'Inf', 'NaN', 'NA', 'NA_integer_', 'NA_real_', 
                         'NA_complex_', 'NA_character_']
            if name in r_reserved:
                name = name + '_'
            
            if original_name != name:
                changes_made = True
            
            valid_names.append(name)
        
        return valid_names, changes_made
    
    # サンプル名を修正
    new_columns, changes_made = make_valid_r_names(df.columns.tolist())
    if changes_made:
        st.write("Sample names have been converted to be R-compatible:")
        for old, new in zip(df.columns.tolist(), new_columns):
            if old != new:
                st.write(f"  '{old}' → '{new}'")
        df.columns = new_columns
############


    st.write('Original gene number:  ' + str(len(df)))

    # floatに変換 誤射悟入
    df = df.astype(float)

    if not float.is_integer(df.iloc[:,0].sum()*1000):
        if test_method == "DESeq2-LRT":
            st.markdown("## It is likely that your data are normalized. Please upload unnormalized raw count data.")

    if test_method == "DESeq2-LRT": #DESeq2は整数化
        df = df.round(0)

    df = df.loc[~(df==0).all(axis=1)] #すべて0のrowを除く

########## excel対応?
    st.write("All zero count genes are removed.")
    if df.isnull().values.sum() > 0:
        st.write("There are " + str(df.isnull().values.sum()) + " NaN in :")
        st.write(df[df.isnull().any(axis=1)])
        convert_nan = st.radio( "NaN:",
        ('remove Nan containing genes', 'conver to 0' ), key='remove Nan containing genes')
        if convert_nan == "conver to 0":
            df = df.fillna(0)
        else:
            df = df.dropna(how='any')
############ sample名に-がある場合はunderscoreへ Rでエラーになる
    if "-" in "".join(df.columns.values):
        st.write("Minus in sample name will be converted to _.")
        new_columns = [x.replace('-','_') for x in df.columns.values]
        df.columns = new_columns
############


  #  st.write(df.head())
    total_count = pd.DataFrame(df.sum()[1:])
    total_count.columns= ['Total counts']
    large_var = False
    if max(total_count['Total counts']) > min(total_count['Total counts']) * 2:
        large_var = True
        st.markdown("### Large difference (>2x) in counts")
        st.write(f"Minimum total counts: {min(total_count['Total counts'])}")
        st.markdown("##### Low count samples can be filtered on the side panel.")
        import matplotlib.pyplot as plt
        import seaborn as sns
        df_sum = pd.DataFrame(df.sum())
        df_sum.columns = ['Counts']


        f1 = calc_barplot(df_sum.T, ylabel = "Total counts")
        st.pyplot(f1)

        f2 = calc_barplot(np.log1p(df), ylabel = "ln(x+1)")
        st.pyplot(f2)

    with st.sidebar:

        sample_threshold = 0

        if large_var:
            st.markdown("### Filter the samples <= counts:")
            sample_threshold = st.number_input("Minimum total cout", value = 0, label_visibility = 'collapsed')
            st.markdown("---")

        if test_method == 'DESeq2-LRT':
            st.markdown("### Filter out weakly-expressed genes before multiple test correction:",help = "independentFiltering default:TRUE 平均正規化カウントに基づいて遺伝子をフィルタリングし、多重検定補正の負担を減らすことで統計的検出力を向上させる")
            independentFiltering = st.checkbox('Yes', value=True)
            st.markdown("---")
        else:
            # For other methods, independentFiltering is not applicable
            independentFiltering = False

        st.markdown("### Furhter filtering of genes")
        st.markdown("""#### 低発現遺伝子の除外はFDRの計算を改善する""")

        st.markdown("#### Filter the genes > counts in all samples:")
        min_threshold = st.number_input("count minimum", value = 0, label_visibility = 'collapsed')
        min_threshold = int(min_threshold)
        st.markdown("#### Filter the genes > counts in at least one sample:")
        max_threshold = st.number_input("count max", value = 0, label_visibility = 'collapsed')
        max_threshold = int(max_threshold)
        st.markdown("---")

        if test_method == 'DESeq2-LRT':
            st.markdown("### Batch correction:")
            sva = st.checkbox('SVA batch removal?')
            sva_calc = True
            if sva:
                sva_calc = st.checkbox('Calculate only 2 surrogate variables? Deselect if want to calculate up to the recommended number.', value = True)
                st.markdown("---")

    if any(df.sum() <= sample_threshold): # count 0の列を除く
        st.markdown('#### There are the samples that have counts <= ' + str(sample_threshold))
        st.write(", ".join(df.columns[df.sum() <= sample_threshold].to_list()))
        st.markdown('##### They are removed. Now data are:')
        df = df.drop(df.columns[df.sum() <= sample_threshold].to_list(), axis = 1)

    st.write(df.head())

    if min_threshold > 0:
        df = df[df.apply(min, axis=1) > min_threshold]
    if max_threshold > 0:
        df = df[df.apply(max, axis=1) > max_threshold]

    st.markdown(f'#### Filtered gene number: {str(len(df))}')


    condition = [str(i) for i in df.columns.tolist()[:]] #error防止
    group_condition = remove_common_suffix(condition) #末尾の共通要素を除く
  #  group_condition = [remove_after_space(x) for x in condition] #スペース以降を除く
    group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字を除く


    st.markdown("##### Add conditions other than group, such as genotype (comma, space, CR separated):")
    genes = st.text_input("genes",label_visibility = 'collapsed')
    gene_list = []
    if len(genes) > 0:
        gene_list = genes.split(' ') #まず空白で分離
        gene_list = list(filter(lambda a: a != '', gene_list)) #空白のみを除く
        if ',' in genes:
            gene_list = sum([x.split(',') for x in gene_list],[]) #sumで平坦化 sum(x, [])
        if '\t' in genes:
            gene_list = sum([x.split('\t') for x in gene_list],[])
        if '\n' in genes:
            gene_list = sum([x.split('\n') for x in gene_list],[])
        gene_list = [a for a in gene_list if a != ''] #空を除く
    condition_col = sum([['Group'], gene_list], [] )

    with st.form("input_groups and batch"):
        df_e = pd.DataFrame(index = condition, columns = condition_col)
        for i in df_e.columns.values:
            df_e[i] = group_condition
        st.write('Set conditions:')
    #    edited_df_e = st.experimental_data_editor(df_e)
        df_e = st.data_editor(df_e)
        submitted = st.form_submit_button("Submit")

    condition = df_e.iloc[:,0].tolist()

    for i in df_e.columns.values:
        st.write(' '.join(df_e.loc[:,i].tolist()))

    # 各変数の型を選択するためのセクション
    st.write('Select variable types:')
    var_types = {}
    cols = st.columns(len(condition_col))
    for i, col in enumerate(condition_col):
        with cols[i]:
            var_types[col] = st.radio(
                f"{col}",
                options=["categorical", "continuous"],
                index=0,  # デフォルトはカテゴリカル
                key=f"type_{col}"
            )
   
# modelを作るための要素をリストにする
    comb = [':'.join(x) for x in  list(itertools.combinations(condition_col, 2))]
#ここで ':'.joint(x)とすると、このあと、モデルを作るときに:以降が除かれる
    selections = selections = sum([condition_col, comb],[])

    # duplicateCorrelation用のblocking変数選択（limma eBayes + RNA-seq countの場合のみ）
    block_variable = None
    if test_method == 'limma eBayes' and limma_count and use_duplicate_correlation:
        st.markdown("---")
        st.markdown("### Select Blocking Variable for duplicateCorrelation:")
        st.markdown("""
**重要**: 重複測定の単位（ドナーID、個体ID等）を示す変数を選択してください。

**使用方法**:
1. 上記のメタデータ入力で、ドナーID列を追加済みであることを確認
   - 例: "Add conditions other than group"に"Donor"を追加
   - メタデータテーブルで各サンプルにドナーIDを入力
2. 下記でその列名を選択

**注意**:
- この変数はdesign matrixには含めず、blockingのみに使用されます
- 同一ドナーから複数サンプル（例: 異なるsubtype）がある場合に使用
- Full modelには他の変数（batch, subtype, sex等）を含めてください
        """)
        block_variable = st.selectbox(
            "Blocking variable (column name from metadata):",
            options=["None"] + list(df_e.columns),
            help="メタデータテーブルで入力したドナーID列の名前を選択。例: Donor, donor_id, patient_id, individual"
        )
        if block_variable == "None":
            block_variable = None
            st.warning("⚠️ Blocking variable not selected. duplicateCorrelation will not be used.")
        else:
            st.info(f"ℹ️ Will use duplicateCorrelation with blocking variable: **{block_variable}**")
            st.write(f"Selected column values: {', '.join(df_e[block_variable].unique())}")
        st.markdown("---")

    null_model = st.checkbox("Null model as reduced model?", value = False, help="帰無モデルをreduced modelにする。つまりfull modelに設定した要因のいずれか／すべてに関連した発現変動を検出。")

    st.markdown("##### Select conditions for full model:")
    full = st.multiselect('fullmodel',selections, label_visibility = 'collapsed')

    # 最初の変数を時間変数として保存
    time_var = None
    if len(full) > 0:
        time_var = full[0]
    #  st.write(time_var)

    if not null_model:
        st.markdown("##### Select conditions for reduced model:")
        reduced = st.multiselect('reducedmodel',selections, label_visibility = 'collapsed')
    else:
        reduced = []
    # 多項式適用のロジック（DESeq2-LRT、limma、GLMでpolynomialがTrueのとき）
    polynomial_enabled = False
    if test_method == 'DESeq2-LRT' and add_polynomial and len(full) > 0 and time_var is not None:
        polynomial_enabled = True
        poly_degree = polynomial_degree
        poly_use_raw = use_raw
        poly_use_poly_function = use_poly_function
    elif test_method == 'limma eBayes' and limma_add_polynomial and len(full) > 0 and time_var is not None:
        polynomial_enabled = True
        poly_degree = limma_polynomial_degree
        poly_use_raw = limma_use_raw
        poly_use_poly_function = limma_use_poly_function
    elif test_method == 'Generalized Linear Model (GLM)' and glm_add_polynomial and len(full) > 0 and time_var is not None:
        polynomial_enabled = True
        poly_degree = glm_polynomial_degree
        poly_use_raw = glm_use_raw
        poly_use_poly_function = glm_use_poly_function
    
    if polynomial_enabled:
        
        if poly_use_poly_function:
            # poly()関数を使った実装
            if poly_degree == 2:
                # raw=TRUEが指定されている場合は追加
                raw_param = ", raw=TRUE" if poly_use_raw else ""
                full[0] = f"poly({time_var}, degree=2{raw_param})"
            else:  # 3次多項式
                raw_param = ", raw=TRUE" if poly_use_raw else ""
                full[0] = f"poly({time_var}, degree=3{raw_param})"
            
            st.markdown(f"##### Using {'raw' if poly_use_raw else 'orthogonal'} polynomial term for {time_var}: {full[0]}")
        else:
            # I()関数を使った実装（明示的なべき乗）
            if poly_degree == 2:
                # 元の時間変数を置き換え、2次の項を追加
                new_terms = [time_var, f"I({time_var}^2)"]
                # fullの最初の要素を置き換え
                full[0] = new_terms[0]
                # 2次の項を挿入
                full.insert(1, new_terms[1])
                
                st.markdown(f"##### Using explicit powers for {time_var}: {time_var} + I({time_var}^2)")
            else:  # 3次多項式
                # 元の時間変数を置き換え、2次と3次の項を追加
                new_terms = [time_var, f"I({time_var}^2)", f"I({time_var}^3)"]
                # fullの最初の要素を置き換え
                full[0] = new_terms[0]
                # 2次と3次の項を挿入
                full.insert(1, new_terms[1])
                full.insert(2, new_terms[2])
                
                st.markdown(f"##### Using explicit powers for {time_var}: {time_var} + I({time_var}^2) + I({time_var}^3)")


    full = [x.replace(':','\:') for x in full] # :のままだと文字列が消失する
    reduced = [x.replace(':','\:') for x in reduced]

    if len(reduced) > 0 and not null_model:
        null_model = False

    full_model = "~ " + " + ".join(full)
    if null_model:
        reduced_model = "~ 1"
    elif len(reduced) == 0: #reducedを指定していないときはnull modelにする
        reduced_model = "~ 1"
        st.markdown("#### Null model is uses as reduced model.")
    else:
        reduced_model = "~ " + " + ".join(reduced)
    st.markdown("##### Full model:  " + full_model)
    st.markdown("##### Reduced model:  " + reduced_model)
    st.markdown("""
timeは第1項にする。\n
Full modelとReduced modelとの違いが検定される。
例えばgenotypと時系列のデータのときにgenptypeは関係なく、時系列変化をする遺伝子を検出する場合は
~ time + genotype と ~ genotype の比較となる。\n
もし、genotype特異的で時間で変化する遺伝子を検出する場合は
~ time + genotype + genotype\:time と ~ time + genotype の比較となる。\n
\n
Reduced modelにnull modelを設定するとFull modelの要因で変化する遺伝子を検出する。\n
例えばWTの細胞の時系列データだけの場合、timeをFull modelにNull modelをReduced modelにする。\n
        """)

    if (len(condition) != len(df.columns)):
            st.write("The number of group name does not match the data.")

#    df_condition = pd.DataFrame(condition)
#    df_batch = pd.DataFrame(batch)

# 1-Marなどの誤変換への対応
    check_excel_autoconversion(df)

    if len(df.index.values) != len(set(df.index.values)):
#        st.markdown("#### There are duplicated rows. Converting the names...")
#        st.write("The gene name of the second occurrence has _2 at the end.")
#        lis = df.index.values
#        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]
        df = rename_duplicates(df)


    # 多項式次数が1より大きい場合の時間変数チェック
    if (test_method == 'Beta Regression' and 'beta_polynomial_degree' in locals() and beta_polynomial_degree > 1) or (test_method == 'DESeq2-LRT' and 'add_polynomial' in locals() and add_polynomial) or (test_method == 'limma eBayes' and 'limma_add_polynomial' in locals() and limma_add_polynomial) or (test_method == 'Generalized Linear Model (GLM)' and 'glm_add_polynomial' in locals() and glm_add_polynomial):
        if len(full) > 0:
          #  st.write("Using polynomial")
            try:
                coldata_file = os.path.join(temp_dir, 'coldata.tsv')
                df_e.to_csv(coldata_file, sep='\t', index=False)
                coldata = pd.read_table(coldata_file)
                if time_var in coldata.columns:
                    # 時間変数の型チェック
                    time_col = coldata[time_var]
                    is_numeric = pd.api.types.is_numeric_dtype(time_col)
                    
                    if not is_numeric:
                        # 数値に変換可能かチェック
                        try:
                            # 数字だけを抽出する正規表現パターン
                            numeric_values = time_col.str.extract(r'(\d+\.?\d*)')[0].astype(float)
                            st.info(f"情報: 時間変数 '{time_var}' は文字列ですが、数値として抽出できます。解析時に自動的に変換されます。")
                            coldata[time_var] = numeric_values #coldataを変えておく
                            
                            # ユニークポイント数チェック
                            unique_points = len(numeric_values.unique())
                            
                            # 多項式次数のチェック
                            current_poly_degree = None
                            if test_method == 'Beta Regression' and 'beta_polynomial_degree' in locals():
                                current_poly_degree = beta_polynomial_degree
                            elif test_method == 'DESeq2-LRT' and 'polynomial_degree' in locals():
                                current_poly_degree = polynomial_degree
                            elif test_method == 'limma eBayes' and 'limma_polynomial_degree' in locals():
                                current_poly_degree = limma_polynomial_degree
                            elif test_method == 'Generalized Linear Model (GLM)' and 'glm_polynomial_degree' in locals():
                                current_poly_degree = glm_polynomial_degree
                            
                            if current_poly_degree and current_poly_degree >= unique_points:
                                st.error(f"エラー: 多項式の次数（{current_poly_degree}）がユニークなタイムポイント数（{unique_points}）以上です。")
                                st.error(f"あなたのデータで使用可能な最大次数: {unique_points - 1}")
                                st.error(f"タイムポイント: {sorted(numeric_values.unique())}")
                                st.stop()
                                
                        except:
                            current_poly_degree = None
                            if test_method == 'Beta Regression' and 'beta_polynomial_degree' in locals():
                                current_poly_degree = beta_polynomial_degree
                            elif test_method == 'DESeq2-LRT' and 'polynomial_degree' in locals():
                                current_poly_degree = polynomial_degree
                            elif test_method == 'limma eBayes' and 'limma_polynomial_degree' in locals():
                                current_poly_degree = limma_polynomial_degree
                            elif test_method == 'Generalized Linear Model (GLM)' and 'glm_polynomial_degree' in locals():
                                current_poly_degree = glm_polynomial_degree
                                
                            if current_poly_degree:
                                st.warning(f"警告: 時間変数 '{time_var}' は数値ではありません。多項式モデル（次数{current_poly_degree}）を使用するには数値が必要です。モデルが収束しない可能性があります。")
                else:
                    st.warning(f"警告: 時間変数 '{time_var}' が実験デザインファイルに見つかりません。")
            except Exception as e:
                st.warning(f"実験デザインファイルの読み込み中にエラーが発生しました: {str(e)}")

    st.markdown("""
--------------------------------------------------------------------------
        """)

    # Collinearity check before running analysis (for limma eBayes)
    collinear_pairs = []
    if test_method == 'limma eBayes' and df_e is not None:
        import re
        from itertools import combinations
        # Extract variable names from formulas
        full_vars = re.findall(r'\b[a-zA-Z_]\w*\b', full_model.replace('~', '').replace(':', '_').replace('*', ' '))
        full_vars = [v for v in full_vars if v in df_e.columns]

        if len(full_vars) >= 2:
            for var1, var2 in combinations(full_vars, 2):
                if df_e[var1].dtype == 'object' and df_e[var2].dtype == 'object':
                    contingency = df_e.groupby([var1, var2]).size().unstack(fill_value=0)
                    row_sums = (contingency > 0).sum(axis=1)
                    col_sums = (contingency > 0).sum(axis=0)
                    var1_predicts_var2 = (row_sums == 1).all()
                    var2_predicts_var1 = (col_sums == 1).all()
                    if var1_predicts_var2 or var2_predicts_var1:
                        collinear_pairs.append((var1, var2, contingency))

    # Show collinearity warning and bypass option
    skip_collinearity_check = False
    if collinear_pairs:
        st.warning("⚠️ **Potential collinearity detected!**")
        for var1, var2, contingency in collinear_pairs:
            with st.expander(f"📊 {var1} vs {var2}"):
                st.dataframe(contingency)
        # Use session_state to persist checkbox state across reruns
        if 'skip_collinearity' not in st.session_state:
            st.session_state.skip_collinearity = False
        skip_collinearity_check = st.checkbox(
            "⚡ Proceed anyway (let R handle collinearity)",
            value=st.session_state.skip_collinearity,
            key="skip_collinearity"
        )

    if st.button('Run analysis'):
        #まずRのdfに変換
        if test_method == 'DESeq2-LRT':
    #        ro.r.assign('cts',cts) # エラーが出るのでファイルに一旦保存する
            r.assign('df',df)
            r.assign('df_e',df_e)

            # 変数タイプ情報をR側に渡す
            continuous_vars = [col for col in condition_col if var_types[col] == "continuous"]
            r_continuous_vars = ro.StrVector(continuous_vars)
            ro.r.assign('continuous_vars', r_continuous_vars)

            pyper_df_path = "saveRDS(df, '" + temp_dir + "/pyper_df.RDS')"
            r(pyper_df_path)
            pyper_df_e_path = "saveRDS(df_e, '" + temp_dir + "/pyper_df_e.RDS')"
            r(pyper_df_e_path)
            read_pyper_df = "cts <- readRDS('" + temp_dir + "/pyper_df.RDS')"
            read_pyper_df_e = "coldata <- readRDS('" + temp_dir + "/pyper_df_e.RDS')"
            ro.r(read_pyper_df)
            ro.r(read_pyper_df_e)
            #まずベクターに変換
            r_condition =  ro.StrVector(condition)
            ro.r.assign('condition', r_condition)
            full_model = full_model.replace('\:',':')
            reduced_model = reduced_model.replace('\:',':')
            ro.r.assign('full_model', full_model)
            ro.r.assign('reduced_model', reduced_model)
            ro.r.assign('sva',sva)
            ro.r.assign('sva_calc', sva_calc)
            ro.r.assign('independentFiltering', independentFiltering)
            ro.r.assign('res_dir', res_dir)
            ro.r.assign('temp_dir', temp_dir)

            if 'add_polynomial' in locals() and add_polynomial:
                # 変数をRに渡す
                if len(full) > 0:
                    polynomial_variable = full[0]
                    ro.r.assign('add_polynomial', True)
                    ro.r.assign('polynomial_degree', polynomial_degree)
                    ro.r.assign('polynomial_variable', polynomial_variable)
                    ro.r.assign('use_raw', use_raw)
                    ro.r.assign('use_poly_function', use_poly_function)
                    #counts_file = os.path.join(temp_dir, 'counts.tsv')
                    #df.to_csv(counts_file, sep='\t')

                else:
                    st.error("Cannot use polynomial terms: full model is empty")
                    ro.r.assign('add_polynomial', False)
            else:
                ro.r.assign('add_polynomial', False)

            with st.spinner('Calculating DESeq2...'):
                ro.r('calc_dds_LRT()')


            image = Image.open(res_dir + '/DispersionEstimates.png')
            st.image(image, caption='Despersion Estimates')

            res_df = pd.read_csv(res_dir + '/DESeq2_LRT_res.tsv', sep = '\t', index_col= 0)

            st.warning("""
            ⚠️ **LRT結果のlog2FoldChangeについての重要な注意:**  
            log2FoldChange列はfullモデルの最後の変数の係数を表示しています。
            これはLRTが検定している効果とは一致しない可能性があります。
            - **LRTが検定**: モデル間の改善（複数の効果を同時に検定する場合もある）
            - **log2FCが表示**: 最後の係数の効果のみ

            - baseMean: 全サンプルの平均正規化カウント
            - log2FoldChange: デザイン式の最後の変数の効果サイズ（LRTの検定対象とは異なる可能性！）
            - lfcSE: log2FoldChangeの標準誤差
            - stat: LRT検定統計量（カイ二乗値）
            - pvalue: LRTのP値（fullモデル vs reducedモデルの検定）
            - padj: 多重検定補正後のP値

            """)

            # 日本語の説明ファイルを作成
            with open(os.path.join(res_dir, '必ずお読みください.txt'), 'w', encoding='utf-8') as f:
                f.write("""
            DESeq2 LRT結果についての重要な説明
            =====================================

            DESeq2_LRT_res.tsvの各列の意味：
            - baseMean: 全サンプルの平均正規化カウント ✓
            - log2FoldChange: デザイン式の最後の変数の効果サイズ（LRTの検定対象とは異なる可能性！）
            - lfcSE: log2FoldChangeの標準誤差
            - stat: LRT検定統計量（カイ二乗値） ✓
            - pvalue: LRTのP値（fullモデル vs reducedモデルの検定） ✓
            - padj: 多重検定補正後のP値 ✓

            【注意】log2FoldChangeはLRTが検定している効果とは異なる場合があります。

            例：
            - time効果を検定しているが、genotypeが最後の変数の場合
              → LRTはtime効果を検定、log2FCはgenotype効果を表示

            正しく解釈するために：
            1. どのモデルが比較されたか確認（full vs reduced）
            2. その差分がLRTの検定対象
            3. log2FCは参考値であり、誤解を招く可能性がある

            【推奨】
            - 方向性（上昇/下降）が必要な場合は、別途Wald testを実行
            - または、検定したい変数を最後に配置してモデルを再構築
            """)

            st.write("FDR < 0.05: " + str(len(res_df.loc[(res_df['padj']<0.05),])))

            res_df= res_df.loc[(res_df['padj']<0.1),'padj']
            st.dataframe(res_df)
            if sva:
                st.markdown("#### =======SVA=======")
                with st.spinner('Preparing SVAseq...'):
                    sva_n = ro.r("sv_n <- calc_sva_n()")
                st.write("Recommended number of SVA covariates: " + str(int(sva_n[0])))
                with st.spinner('Calculating SVAseq...'):
                    ro.r("calc_svseq_LRT()")

            if sva:
                st.session_state.deseq2lrt = read_csv(res_dir + "/SVA_LRT_res.tsv", sep = '\t', index_col=0)
            else:
                st.session_state.deseq2lrt = read_csv(res_dir + "/DESeq2_LRT_res.tsv", sep = '\t', index_col=0)


            file_name = file_name_head + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")

            shutil.make_archive("res", format='zip',root_dir= res_dir)

        elif test_method == 'limma eBayes':
            # Check if collinearity was detected and user hasn't bypassed
            if 'collinear_pairs' in dir() and collinear_pairs and not skip_collinearity_check:
                st.error("⚠️ Collinearity detected. Check the bypass option above to proceed.")
                st.stop()

            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            # Save input to files for R import
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            voom_plot_path = os.path.join(res_dir, 'voom_plot.png')
        #    if os.path.exists(voom_plot_path):
        #        os.remove(voom_plot_path)

            ro.r.assign('temp_dir', temp_dir)

            # SVA関連の変数をRに渡す
            if 'use_sva' in locals() and use_sva:
                ro.r.assign('use_sva', True)
                ro.r.assign('sva_n_sv', sva_n_sv)
            else:
                ro.r.assign('use_sva', False)
                ro.r.assign('sva_n_sv', 2)

            # 変数タイプ情報をR側に渡す
            continuous_vars = [col for col in condition_col if var_types[col] == "continuous"]
            r_continuous_vars = ro.StrVector(continuous_vars)
            ro.r.assign('continuous_vars', r_continuous_vars)

            # polynomial関連の変数をRに渡す
            if 'limma_add_polynomial' in locals() and limma_add_polynomial and len(full) > 0:
                polynomial_variable = full[0]
                ro.r.assign('add_polynomial', True)
                ro.r.assign('polynomial_degree', limma_polynomial_degree)
                ro.r.assign('polynomial_variable', polynomial_variable)
                ro.r.assign('use_raw', limma_use_raw)
                ro.r.assign('use_poly_function', limma_use_poly_function)
            else:
                ro.r.assign('add_polynomial', False)

            if apply_logit:
                # For logit-transformed data, use this R code
                robust_param = "TRUE" if use_robust_ebayes else "FALSE"
                r_code = f"""
                sink()
                sink(paste0(temp_dir, "/limma_output.txt"))
                library(limma)
                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')

                # For logit transformed data, apply the transformation in R instead
                eps <- 1e-6
                counts <- pmax(counts, eps)
                counts <- pmin(counts, 1-eps)
                counts <- log(counts/(1-counts))

                # Create design matrices
                design_full <- model.matrix(as.formula('{full_model}'), data=coldata)
                design_reduced <- model.matrix(as.formula('{reduced_model}'), data=coldata)

                svobj <- NULL  # SVA not implemented for logit-transformed data yet

                # Identify coefficients specific to the full model
                add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))

                # For logit-transformed data, skip voom and directly fit with limma
                fit_full <- lmFit(counts, design_full)
                fit_full <- eBayes(fit_full, robust={robust_param})
                cat('Used eBayes with robust={robust_param} for logit-transformed data\\n')

                if (length(add_coefs) == 1) {{
                  res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                }} else {{
                  cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                  colnames(cm) <- add_coefs
                  for (i in 1:length(add_coefs)) {{
                    cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                  }}

                  fit_contrast <- contrasts.fit(fit_full, cm)
                  fit_contrast <- eBayes(fit_contrast, robust={robust_param})

                  res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                }}
                output_file <- paste0('{res_dir}/', if(use_sva && !is.null(svobj)) 'sva{sva_n_sv}_' else '', 'limma_res.tsv')
                write.table(res, file=output_file, sep='\t', quote=FALSE, col.names=NA)
                cat('Results written to:', output_file, '\n')
                cat('SVA was', if(use_sva && !is.null(svobj)) 'applied' else 'NOT applied', '\n')
                sink()
                """
            elif limma_count:
                # For RNA-seq count data with advanced options
                log_file = f"{res_dir}/limma_debug.log"

                # Prepare block variable if duplicateCorrelation is used
                block_vector = None
                if use_duplicate_correlation and block_variable:
                    block_vector = df_e[block_variable].values
                    ro.r.assign('block_vector', StrVector(block_vector))
                    st.info(f"ℹ️ Using duplicateCorrelation with blocking variable: {block_variable}")

                # Build R code with conditional blocks
                r_code_parts = [f"""
                # Disable X11 graphics device
                options(bitmapType='cairo')
                pdf(NULL)

                # ログファイルにリダイレクト
                sink('{log_file}', append=FALSE, split=TRUE)
                library(edgeR)
                library(limma)
                tryCatch({{
                    counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
                    counts <- round(counts)
                    coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')
                    y <- DGEList(counts=counts)
                """]

                # filterByExpr with custom parameters or no filtering
                if use_edgeR_filter:
                    r_code_parts.append(f"""
                    # Apply filterByExpr with custom parameters
                    design_temp <- model.matrix(as.formula('{full_model}'), data=coldata)
                    keep <- filterByExpr(y, design=design_temp,
                                        min.count={filter_min_count},
                                        min.total.count={filter_min_total_count},
                                        min.prop={filter_min_prop},
                                        large.n={filter_large_n})
                    y <- y[keep, , keep.lib.sizes=FALSE]
                    cat('Genes retained after filterByExpr:', sum(keep), 'out of', length(keep), '\\n')
                    """)
                else:
                    r_code_parts.append("""
                    # No filtering applied (user disabled filterByExpr)
                    cat('No gene filtering applied\\n')
                    """)

                r_code_parts.append("""
                    # Normalize
                    y <- calcNormFactors(y)

                    # Create initial design matrices
                    design_full_base <- model.matrix(as.formula('""" + full_model + """'), data=coldata)
                    design_reduced_base <- model.matrix(as.formula('""" + reduced_model + """'), data=coldata)
                    """)

                # SVA処理を追加
                if use_sva:
                    r_code_parts.append(f"""
                    # SVA for batch correction
                    cat('\\n===== SVA ANALYSIS =====\\n')
                    suppressMessages(library(sva))

                    # Prepare data for SVA
                    # svaseq expects normalized count-like data (not log-transformed)
                    # Use TMM-normalized CPM values (cpm() uses normalization factors from calcNormFactors)
                    # This accounts for library size differences and TMM normalization
                    dat_for_sva <- cpm(y, log=FALSE)

                    # Estimate surrogate variables using TMM-normalized CPM data
                    cat('Estimating', {sva_n_sv}, 'surrogate variables...\\n')
                    cat('Input data range: min=', min(dat_for_sva), ' max=', max(dat_for_sva), '\\n')
                    svobj <- tryCatch({{
                        svaseq(dat_for_sva, design_full_base, design_reduced_base, n.sv={sva_n_sv})
                    }}, error = function(e) {{
                        cat('SVA estimation failed:', e$message, '\\n')
                        cat('Continuing without SVA\\n')
                        return(NULL)
                    }})

                    if (!is.null(svobj)) {{
                        cat('Successfully estimated', ncol(svobj$sv), 'surrogate variables\\n')

                        # Visualize SVs
                        png(paste0('{res_dir}/sva_plots.png'), width=1200, height=600)
                        par(mfrow=c(1, min(ncol(svobj$sv), 3)))
                        for (i in 1:min(ncol(svobj$sv), 3)) {{
                            # Get condition factor for coloring
                            cond_col <- 1
                            for (cname in colnames(coldata)) {{
                                if (length(unique(coldata[[cname]])) > 1 && length(unique(coldata[[cname]])) < nrow(coldata)) {{
                                    cond_col <- as.numeric(factor(coldata[[cname]]))
                                    break
                                }}
                            }}
                            plot(svobj$sv[,i], pch=19, col=cond_col,
                                 main=paste("SV", i),
                                 xlab="Sample", ylab="Surrogate Variable Value")
                        }}
                        dev.off()

                        # Add SVs to design matrices
                        design_full <- cbind(design_full_base, svobj$sv)
                        design_reduced <- cbind(design_reduced_base, svobj$sv)

                        sv_names <- paste0("SV", 1:ncol(svobj$sv))
                        colnames(design_full) <- c(colnames(design_full_base), sv_names)
                        colnames(design_reduced) <- c(colnames(design_reduced_base), sv_names)

                        cat('Added', ncol(svobj$sv), 'SVs to design matrices\\n')
                        cat('Design matrix columns (with SVs):', colnames(design_full), '\\n')
                    }} else {{
                        # SVA failed, use original designs
                        design_full <- design_full_base
                        design_reduced <- design_reduced_base
                        cat('Using original design matrices without SVs\\n')
                    }}
                    """)
                else:
                    r_code_parts.append("""
                    # No SVA - use base design matrices
                    design_full <- design_full_base
                    design_reduced <- design_reduced_base
                    svobj <- NULL  # Explicitly set to NULL when SVA is not used
                    """)

                r_code_parts.append("""
                    # Check design matrix rank
                    cat('\\n===== DESIGN MATRIX CHECK =====\\n')
                    cat('Full model formula:', '""" + full_model + """', '\\n')
                    cat('Full design matrix rank:', qr(design_full)$rank, '\\n')
                    cat('Full design matrix ncol:', ncol(design_full), '\\n')
                    cat('Full design matrix columns:', colnames(design_full), '\\n')

                    if (qr(design_full)$rank < ncol(design_full)) {
                        stop('Design matrix is not of full rank. Coefficients not estimable: ',
                             paste(colnames(design_full)[qr(design_full)$pivot[(qr(design_full)$rank+1):ncol(design_full)]], collapse=', '),
                             '\\nThis indicates collinearity between variables. Please check your experimental design.')
                    }

                    add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))
                """)

                # Voom transformation
                if which_voom == 'voom':
                    r_code_parts.append(f"""
                    png('{res_dir}/voom_plot.png')
                    v <- voom(y, design_full, plot=TRUE)
                    dev.off()
                    """)
                elif which_voom == 'voomWithQualityWeights':
                    r_code_parts.append(f"""
                    png('{res_dir}/voom_plot.png')
                    v <- voomWithQualityWeights(y, design_full, plot=TRUE)
                    dev.off()
                    cat('Sample quality weights:', v$sample.weights, '\\n')
                    """)
                elif which_voom == 'voomLmFit':
                    sample_weights_param = "TRUE" if use_sample_weights else "FALSE"
                    r_code_parts.append(f"""
                    png('{res_dir}/voom_plot.png')
                    # voomLmFit combines voom and lmFit
                    """)

                # duplicateCorrelation and lmFit
                if use_duplicate_correlation and block_variable:
                    if which_voom == 'voomLmFit':
                        r_code_parts.append(f"""
                    fit_full <- voomLmFit(y, design_full, plot=TRUE, sample.weights={sample_weights_param}, block=block_vector)
                    dev.off()
                    cat('Used voomLmFit with sample.weights={sample_weights_param} and block variable\\n')
                    """)
                    else:
                        r_code_parts.append(f"""
                    # Estimate correlation for repeated measures
                    corfit <- duplicateCorrelation(v, design_full, block=block_vector)
                    cat('Consensus correlation:', corfit$consensus.correlation, '\\n')

                    # Re-run voom with correlation
                    png('{res_dir}/voom_plot_with_correlation.png')
                    v <- {which_voom}(y, design_full, plot=TRUE, block=block_vector, correlation=corfit$consensus)
                    dev.off()

                    # Re-estimate correlation with updated voom weights
                    corfit <- duplicateCorrelation(v, design_full, block=block_vector)
                    cat('Updated consensus correlation:', corfit$consensus.correlation, '\\n')

                    # Fit with correlation
                    fit_full <- lmFit(v, design_full, block=block_vector, correlation=corfit$consensus)
                    """)
                else:
                    # No duplicateCorrelation
                    if which_voom == 'voomLmFit':
                        r_code_parts.append(f"""
                    fit_full <- voomLmFit(y, design_full, plot=TRUE, sample.weights={sample_weights_param})
                    dev.off()
                    cat('Used voomLmFit with sample.weights={sample_weights_param}\\n')
                    """)
                    else:
                        r_code_parts.append("""
                    fit_full <- lmFit(v, design_full)
                    """)

                # eBayes
                robust_param = "TRUE" if use_robust_ebayes else "FALSE"
                r_code_parts.append(f"""
                    fit_full <- eBayes(fit_full, robust={robust_param})
                    cat('Used eBayes with robust={robust_param}\\n')

                    # Extract results
                    if (length(add_coefs) == 1) {{
                      res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                    }} else {{
                      cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                      colnames(cm) <- add_coefs
                      for (i in 1:length(add_coefs)) {{
                        cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                      }}

                      fit_contrast <- contrasts.fit(fit_full, cm)
                      fit_contrast <- eBayes(fit_contrast, robust={robust_param})
                      res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                    }}
                    output_file <- paste0('{res_dir}/', if(use_sva && !is.null(svobj)) 'sva{sva_n_sv}_' else '', 'limma_res.tsv')
                    write.table(res, file=output_file, sep='\\t', quote=FALSE, col.names=NA)
                    cat('Results written to:', output_file, '\\n')
                    cat('SVA was', if(use_sva && !is.null(svobj)) 'applied' else 'NOT applied', '\\n')

                    }}, error = function(e) {{
                        cat("\\n===== ERROR =====\\n")
                        cat("Error in limma analysis:", conditionMessage(e), "\\n")
                        cat("Error occurred at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
                        traceback()
                    }}, finally = {{
                        cat("\\n===== ANALYSIS COMPLETE =====\\n")
                        cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
                        sink()
                    }})
                """)

                r_code = "".join(r_code_parts)

            else:
                # Non-count data
                robust_param = "TRUE" if use_robust_ebayes else "FALSE"

                sva_code = ""
                if use_sva:
                    sva_code = f"""
                # SVA for batch correction
                cat('\\n===== SVA ANALYSIS =====\\n')
                suppressMessages(library(sva))

                # Estimate surrogate variables using raw data
                cat('Estimating', {sva_n_sv}, 'surrogate variables...\\n')
                svobj <- tryCatch({{
                    sva(as.matrix(counts), design_full_base, design_reduced_base, n.sv={sva_n_sv})
                }}, error = function(e) {{
                    cat('SVA estimation failed:', e$message, '\\n')
                    cat('Continuing without SVA\\n')
                    return(NULL)
                }})

                if (!is.null(svobj)) {{
                    cat('Successfully estimated', ncol(svobj$sv), 'surrogate variables\\n')

                    # Visualize SVs
                    png(paste0('{res_dir}/sva_plots.png'), width=1200, height=600)
                    par(mfrow=c(1, min(ncol(svobj$sv), 3)))
                    for (i in 1:min(ncol(svobj$sv), 3)) {{
                        cond_col <- 1
                        for (cname in colnames(coldata)) {{
                            if (length(unique(coldata[[cname]])) > 1 && length(unique(coldata[[cname]])) < nrow(coldata)) {{
                                cond_col <- as.numeric(factor(coldata[[cname]]))
                                break
                            }}
                        }}
                        plot(svobj$sv[,i], pch=19, col=cond_col,
                             main=paste("SV", i),
                             xlab="Sample", ylab="Surrogate Variable Value")
                    }}
                    dev.off()

                    # Add SVs to design matrices
                    design_full <- cbind(design_full_base, svobj$sv)
                    design_reduced <- cbind(design_reduced_base, svobj$sv)

                    sv_names <- paste0("SV", 1:ncol(svobj$sv))
                    colnames(design_full) <- c(colnames(design_full_base), sv_names)
                    colnames(design_reduced) <- c(colnames(design_reduced_base), sv_names)

                    cat('Added', ncol(svobj$sv), 'SVs to design matrices\\n')
                }} else {{
                    design_full <- design_full_base
                    design_reduced <- design_reduced_base
                }}
                """
                else:
                    sva_code = """
                design_full <- design_full_base
                design_reduced <- design_reduced_base
                svobj <- NULL  # Explicitly set to NULL when SVA is not used
                """

                r_code = f"""
                # Disable X11 graphics device
                options(bitmapType='cairo')
                pdf(NULL)

                sink()
                sink(paste0(temp_dir, "/limma_output.txt"))
                library(limma)
                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')

                # Create base design matrices
                design_full_base <- model.matrix(as.formula('{full_model}'), data=coldata)
                design_reduced_base <- model.matrix(as.formula('{reduced_model}'), data=coldata)

                {sva_code}

                add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))

                # 非カウントデータなのでvoomはスキップ
                fit_full <- lmFit(counts, design_full)
                fit_full <- eBayes(fit_full, robust={robust_param})
                cat('Used eBayes with robust={robust_param} for non-count data\\n')

                if (length(add_coefs) == 1) {{
                  res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                }} else {{
                  cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                  colnames(cm) <- add_coefs
                  for (i in 1:length(add_coefs)) {{
                    cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                  }}

                  fit_contrast <- contrasts.fit(fit_full, cm)
                  fit_contrast <- eBayes(fit_contrast, robust={robust_param})

                  res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                }}
                output_file <- paste0('{res_dir}/', if(use_sva && !is.null(svobj)) 'sva{sva_n_sv}_' else '', 'limma_res.tsv')
                write.table(res, file=output_file, sep='\t', quote=FALSE, col.names=NA)
                cat('Results written to:', output_file, '\n')
                cat('SVA was', if(use_sva && !is.null(svobj)) 'applied' else 'NOT applied', '\n')
                sink()
                """

            # Save R code for debugging
            r_code_file = os.path.join(res_dir, 'generated_r_code.R')
            with open(r_code_file, 'w') as f:
                f.write(r_code)
            st.info(f"ℹ️ Generated R code saved to: {r_code_file}")

            ro.r(r_code)

            # Check for log file and display any errors
            # Try both possible log file locations
            log_file = os.path.join(res_dir, 'limma_debug.log')  # RNA-seq count
            if not os.path.exists(log_file):
                log_file = os.path.join(temp_dir, 'limma_output.txt')  # Non-count

            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                if 'ERROR' in log_content or 'Error' in log_content:
                    st.error("Error occurred during limma analysis:")
                    st.code(log_content)

                # Always show log in an expander for debugging
                with st.expander("📋 View analysis log"):
                    st.code(log_content)

            # Display SVA plots if they exist AND SVA was actually applied (check log)
            sva_plot_path = f"{res_dir}/sva_plots.png"
            sva_actually_applied = False
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_text = f.read()
                    sva_actually_applied = "===== SVA ANALYSIS =====" in log_text and "Successfully estimated" in log_text

            if sva_actually_applied and os.path.exists(sva_plot_path):
                st.success("✅ SVA analysis was performed!")
                st.image(sva_plot_path, caption='Surrogate Variables')
            elif use_sva and not sva_actually_applied:
                st.warning("⚠️ SVA was enabled but not applied. Check the log file for details.")

            # Display voom plot if it exists (only for RNA-seq count data)
            voom_plot_path = f"{res_dir}/voom_plot.png"
            if os.path.exists(voom_plot_path):
                st.image(voom_plot_path, caption='Voom mean-variance trend')

            # Check if result file exists before reading
            # Try SVA result file first if SVA was enabled
            sva_filename = f'sva{sva_n_sv}_limma_res.tsv'
            if use_sva and os.path.exists(os.path.join(res_dir, sva_filename)):
                limma_res_file = os.path.join(res_dir, sva_filename)
                st.info(f"📊 Reading SVA-corrected results (n_sv={sva_n_sv})")
            else:
                limma_res_file = os.path.join(res_dir, 'limma_res.tsv')
                if use_sva:
                    st.warning(f"⚠️ SVA was enabled but {sva_filename} not found. Using standard results.")

            if not os.path.exists(limma_res_file):
                st.error(f"Result file not found: {limma_res_file}")
                st.info("Checking log file for errors...")
                if os.path.exists(log_file):
                    with open(log_file, 'r') as f:
                        st.code(f.read())
                st.stop()

            res_df = pd.read_csv(limma_res_file, sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)

            sva_prefix = f"sva{sva_n_sv}_" if (use_sva and os.path.exists(os.path.join(res_dir, sva_filename))) else ""
            file_name = file_name_head + "_" + sva_prefix + "limma_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")

            shutil.make_archive("res", format='zip',root_dir= res_dir)


        elif test_method == 'edgeR-QL':
            # edgeR-QL for LRT-like analysis (RNA-seq count data)
            st.info(f"ℹ️ Starting edgeR-QL analysis with {len(df)} genes")

            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            ro.r.assign('temp_dir', temp_dir)

            log_file = f"{res_dir}/edger_debug.log"
            robust_param = "TRUE" if use_robust_edger else "FALSE"

            # Build R code
            r_code_parts = [f"""
            sink('{log_file}', append=FALSE, split=TRUE)
            library(edgeR)
            library(limma)
            tryCatch({{
                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\\t')
                counts <- round(counts)
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\\t')
                y <- DGEList(counts=counts)
            """]

            # filterByExpr
            if use_edgeR_filter:
                r_code_parts.append(f"""
                # Apply filterByExpr with custom parameters
                design_temp <- model.matrix(as.formula('{full_model}'), data=coldata)
                keep <- filterByExpr(y, design=design_temp,
                                    min.count={filter_min_count},
                                    min.total.count={filter_min_total_count},
                                    min.prop={filter_min_prop},
                                    large.n={filter_large_n})
                y <- y[keep, , keep.lib.sizes=FALSE]
                cat('Genes retained after filterByExpr:', sum(keep), 'out of', length(keep), '\\n')
                """)
            else:
                r_code_parts.append("""
                # No filtering applied
                cat('No gene filtering applied\\n')
                """)

            r_code_parts.append(f"""
                # Normalize and estimate dispersion
                y <- calcNormFactors(y, method="TMM")
                design_full <- model.matrix(as.formula('{full_model}'), data=coldata)
                design_reduced <- model.matrix(as.formula('{reduced_model}'), data=coldata)

                # Estimate dispersion and fit QL model
                y <- estimateDisp(y, design_full, robust={robust_param})
                fit <- glmQLFit(y, design_full, robust={robust_param})
                cat('Used robust={robust_param} for edgeR-QL\\n')

                # Identify coefficients to test (LRT-like: full vs reduced)
                add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))
                cat('Testing coefficients:', paste(add_coefs, collapse=", "), '\\n')
                cat('Number of coefficients being tested:', length(add_coefs), '\\n')

                # Determine test type based on number of coefficients
                if (length(add_coefs) == 1) {{
                    cat('Single coefficient test: will use t-test (logFC direction available)\\n')
                }} else {{
                    cat('Multiple coefficient test: will use F-test (no logFC direction, total effect only)\\n')
                }}

                # Perform QL test
            """)

            if use_glmtreat:
                r_code_parts.append(f"""
                # glmTreat for minimum fold-change testing
                qlf <- glmTreat(fit, coef=add_coefs, lfc={treat_lfc})
                cat('Used glmTreat with lfc={treat_lfc}\\n')
                """)
            else:
                r_code_parts.append("""
                # glmQLFTest for standard testing
                qlf <- glmQLFTest(fit, coef=add_coefs)
                cat('Used glmQLFTest (standard)\\n')
                """)

            r_code_parts.append(f"""
                # Extract results
                res <- topTags(qlf, n=Inf, sort.by="PValue")$table

                # Add gene column name for compatibility with GSEA.py
                res_out <- cbind(gene=rownames(res), res)
                write.table(res_out, file='{res_dir}/edger_res.tsv', sep='\\t', quote=FALSE, row.names=FALSE)

                # Generate BCV plot
                png('{res_dir}/bcv_plot.png')
                plotBCV(y)
                dev.off()

                # Generate MD plot
                png('{res_dir}/md_plot.png')
                plotMD(qlf)
                abline(h=c(-1, 1), col="blue", lty=2)
                dev.off()

                }}, error = function(e) {{
                    cat("\\n===== ERROR =====\\n")
                    cat("Error in edgeR-QL analysis:", conditionMessage(e), "\\n")
                    traceback()
                }}, finally = {{
                    cat("\\n===== ANALYSIS COMPLETE =====\\n")
                    sink()
                }})
            """)

            r_code = "".join(r_code_parts)
            ro.r(r_code)

            # Display filtering and test type information from log
            try:
                log_content = open(log_file, 'r').read()

                # Filtering results
                if use_edgeR_filter:
                    if 'Genes retained after filterByExpr:' in log_content:
                        for line in log_content.split('\n'):
                            if 'Genes retained after filterByExpr:' in line:
                                st.info(f"ℹ️ {line}")
                                break

                # Test type information
                if 'Number of coefficients being tested:' in log_content:
                    for line in log_content.split('\n'):
                        if 'Testing coefficients:' in line:
                            st.write(f"**{line.strip()}**")
                        if 'Number of coefficients being tested:' in line:
                            st.write(f"**{line.strip()}**")
                        if 'Single coefficient test:' in line:
                            st.success("✓ Single coefficient test: logFC direction is available (t-test)")
                        if 'Multiple coefficient test:' in line:
                            st.warning("⚠️ Multiple coefficient test: F-test only (total effect, no logFC direction)")
            except:
                pass

            # Display BCV plot
            bcv_plot_path = f"{res_dir}/bcv_plot.png"
            if os.path.exists(bcv_plot_path):
                st.image(bcv_plot_path, caption='BCV (Biological Coefficient of Variation) Plot')

            # Display MD plot
            md_plot_path = f"{res_dir}/md_plot.png"
            if os.path.exists(md_plot_path):
                st.image(md_plot_path, caption='MD (Mean-Difference) Plot')

            # Display results
            # Read with 'gene' column as index (compatible with GSEA.py)
            res_df = pd.read_csv(os.path.join(res_dir, 'edger_res.tsv'), sep='\t', index_col='gene')
            st.markdown("### edgeR-QL Results:")
            st.write(f"**Total genes analyzed**: {len(res_df)}")
            st.write(f"**Significant (FDR<0.05)**: {(res_df['FDR']<0.05).sum()}")
            st.dataframe(res_df)

            file_name = file_name_head + "_edgeR-QL_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")
            shutil.make_archive("res", format='zip', root_dir=res_dir)


        elif test_method == 'Beta Regression':
            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            # ファイルの保存と設定は同じ
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            ro.r.assign('temp_dir', temp_dir)


            r_code = f"""
            sink()
            sink(paste0(temp_dir, "/beta_output.txt"))
            library(betareg)
            library(lmtest)
            library(parallel)

            counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
            coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')

            # 0-1の境界の調整
            eps <- {epsilon}
            counts <- pmax(pmin(counts, 1-eps), eps)

            # 並列処理クラスターの設定
            n_cores <- {n_cores}
            cl <- makeCluster(n_cores)

            # 並列処理に必要なパッケージをクラスターに読み込み
            clusterEvalQ(cl, {{
              library(betareg)
              library(lmtest)
            }})

            time_var <- "{full[0]}"
            cat("time_var")
            cat(time_var)

            # 時間変数の確認と変換
            if(time_var %in% colnames(coldata)) {{
              cat("Time variable exists in coldata. Values:", "\\n")
              print(coldata[[time_var]])
              
              # 時間変数が数値でない場合は変換
              if(!is.numeric(coldata[[time_var]])) {{
                cat("Converting time variable to numeric\\n")
                # 数値抽出と変換
                coldata[[time_var]] <- as.numeric(gsub("[^0-9.]", "", as.character(coldata[[time_var]])))
                cat("After conversion:", "\\n")
                print(coldata[[time_var]])
              }}
            }} else {{
              cat("WARNING: Time variable not found in coldata!\\n")
            }}

            coldata[[time_var]] <- coldata[[time_var]] / max(coldata[[time_var]]) #最大値によるnormalization

            # 多項式の次数に基づく項を構築
            polynomial_terms <- ""
            if ({beta_polynomial_degree} >= 2) {{
              polynomial_terms <- paste0(polynomial_terms, " + I(", time_var, "^2)")
            }}
            if ({beta_polynomial_degree} >= 3) {{
              polynomial_terms <- paste0(polynomial_terms, " + I(", time_var, "^3)")
            }}

            # 変数をクラスターに送信
            clusterExport(cl, c("counts", "coldata", "eps", "time_var", "polynomial_terms"))

            # 処理開始メッセージ
            cat("Starting parallel beta regression on", n_cores, "cores for", nrow(counts), "genes\\n")

            # テストモデル実行
            test_model_result <- tryCatch({{
              # テスト用のデータ
              test_gene_data <- data.frame(y=as.numeric(counts[1,]), coldata)
              
              # フォーミュラを構築
              full_formula <- paste("{full_model.replace('~', '')}", polynomial_terms)
              
              # モデル適合を試行
              test_fit <- betareg(as.formula(paste("y ~", full_formula)), data=test_gene_data)
              "success"
            }}, error=function(e) {{
              # エラーメッセージを返す
              return(conditionMessage(e))
            }})

            # エラーカウンター初期化
            error_counter <- 0
            error_message <- ""

            # 並列処理関数
            process_gene <- function(i) {{
              gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
              
              full_formula <- paste("{full_model.replace('~', '')}", polynomial_terms)
              reduced_formula <- "{reduced_model.replace('~', '')}"
              
              tryCatch({{
                # フルモデルとリデュースドモデルのフィッティング
                full_fit <- betareg(as.formula(paste("y ~", full_formula)), data=gene_data)
                reduced_fit <- betareg(as.formula(paste("y ~", reduced_formula)), data=gene_data)
                
                # 尤度比検定
                lr_test <- lrtest(reduced_fit, full_fit)
                
                # 結果を返す
                c(statistic = lr_test$Chisq[2],
                  df = lr_test$Df[2],
                  p_value = lr_test$`Pr(>Chisq)`[2],
                  logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
              }}, error=function(e) {{
                # エラーが発生した場合はNAを返す
                if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
              }})
            }}

            # 並列処理の実行
            system.time(
              results_list <- parLapply(cl, 1:nrow(counts), process_gene)
            )

            # クラスターの終了
            stopCluster(cl)

            # 結果をマトリックスに変換
            results_matrix <- do.call(rbind, results_list)
            rownames(results_matrix) <- rownames(counts)

            # NULLの数をカウント
            na_count <- sum(is.na(results_matrix[, "statistic"]))
            total_genes <- nrow(results_matrix)
            na_percent <- round(100 * na_count / total_genes, 2)

            # 結果ファイルにモデル情報を追加
            cat("\\n### モデル収束情報 ###\\n", file='{res_dir}/model_convergence_info.txt')
            cat("多項式次数:", {beta_polynomial_degree}, "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            cat("フルモデル式:", paste("y ~", paste("{full_model.replace('~', '')}", polynomial_terms)), "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            cat("縮小モデル式:", paste("y ~", "{reduced_model.replace('~', '')}"), "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)

            if (na_count > 0) {{
              cat("警告: ", na_count, " 個の遺伝子 (", na_percent, "%) でモデルが収束しませんでした。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              
              if (na_count == total_genes) {{
                cat("すべての遺伝子でモデルが収束しませんでした。多項式次数を下げる、時間変数をスケーリング(現在最大値でスケーリングしてある)、不均等の時間を均等化することを検討してください。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
                cat("テストモデルのエラー: ", test_model_result, "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              }} else {{
                cat(total_genes - na_count, " 個の遺伝子 (", 100 - na_percent, "%) で正常に解析できました。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              }}
            }} else {{
              cat("すべての遺伝子で正常にモデルが収束しました。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            }}

            # 多重検定補正
            results_matrix <- cbind(results_matrix, 
                                  adj.P.Val = p.adjust(results_matrix[, "p_value"], method="BH"))

            # 結果の保存
            res_df <- as.data.frame(results_matrix)
            res_df <- res_df[order(res_df$p_value), ]
            write.table(res_df, file='{res_dir}/betareg_res.tsv', sep='\\t', quote=FALSE, col.names=NA)
            sink()
            """

        elif test_method == 'Generalized Linear Model (GLM)':
            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            # ファイルの保存と設定
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            ro.r.assign('temp_dir', temp_dir)
            
            # 変数タイプ情報をR側に渡す
            continuous_vars = [col for col in condition_col if var_types[col] == "continuous"]
            r_continuous_vars = ro.StrVector(continuous_vars)
            ro.r.assign('continuous_vars', r_continuous_vars)
            
            # polynomial関連の変数をRに渡す
            if 'glm_add_polynomial' in locals() and glm_add_polynomial and len(full) > 0:
                polynomial_variable = full[0]
                ro.r.assign('add_polynomial', True)
                ro.r.assign('polynomial_degree', glm_polynomial_degree)
                ro.r.assign('polynomial_variable', polynomial_variable)
                ro.r.assign('use_raw', glm_use_raw)
                ro.r.assign('use_poly_function', glm_use_poly_function)
            else:
                ro.r.assign('add_polynomial', False)

            r_code = f"""
            sink()
            sink(paste0(temp_dir, "/glm_output.txt"))
            library(parallel)

            counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
            coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')

            # 連続変数を処理
            if (exists("continuous_vars") && length(continuous_vars) > 0) {{
                cat("Processing continuous variables...\\n")
                for (col_name in continuous_vars) {{
                    if (col_name %in% colnames(coldata)) {{
                        cat(paste0("Treating '", col_name, "' as continuous variable\\n"))
                        if (!is.numeric(coldata[[col_name]])) {{
                            original_values <- coldata[[col_name]]
                            numeric_values <- as.numeric(gsub("[^0-9.]", "", as.character(original_values)))
                            if (any(is.na(numeric_values))) {{
                                warning(paste0("Cannot convert '", col_name, "' to numeric. Using as factor."))
                                coldata[[col_name]] <- factor(coldata[[col_name]])
                            }} else {{
                                coldata[[col_name]] <- numeric_values
                                cat(paste0("Converted '", col_name, "' to numeric values: ", 
                                         paste(head(numeric_values), collapse=", "), "...\\n"))
                            }}
                        }}
                    }}
                }}
            }}

            # 残りの変数を因子型に変換（連続変数以外）
            for (i in c(1:dim(coldata)[2])) {{
                col_name <- colnames(coldata)[i]
                if (!exists("continuous_vars") || !(col_name %in% continuous_vars)) {{
                    cat(paste0("Treating '", col_name, "' as categorical variable\\n"))
                    coldata[,i] <- factor(coldata[,i])
                }}
            }}

            # 分布ファミリーに応じたデータ前処理
            if ("{glm_dist_short}" == "beta") {{
                # 0-1の境界の調整
                eps <- {glm_epsilon}
                counts <- pmax(pmin(counts, 1-eps), eps)
            }} else if ("{glm_dist_short}" == "gaussian") {{
                # ガウスの場合は特に前処理不要
            }} else if ("{glm_dist_short}" == "poisson" || "{glm_dist_short}" == "nb") {{
                # ポアソンとNBはカウントデータを想定（整数化）
                counts <- round(counts)
            }}

            # 分布ファミリーの設定関数
            get_family <- function() {{
                if ("{glm_dist_short}" == "beta") {{
                    # Betaの場合はbetaregパッケージを使用
                    library(betareg)
                    if ("{glm_link}" == "logit") {{
                        return(list(family = "beta", link = "logit"))
                    }} else if ("{glm_link}" == "probit") {{
                        return(list(family = "beta", link = "probit"))
                    }} else if ("{glm_link}" == "cloglog") {{
                        return(list(family = "beta", link = "cloglog"))
                    }}
                }} else if ("{glm_dist_short}" == "gaussian") {{
                    if ("{glm_link}" == "identity") {{
                        return(gaussian(link = "identity"))
                    }} else if ("{glm_link}" == "log") {{
                        return(gaussian(link = "log"))
                    }} else if ("{glm_link}" == "inverse") {{
                        return(gaussian(link = "inverse"))
                    }}
                }} else if ("{glm_dist_short}" == "poisson") {{
                    if ("{glm_link}" == "log") {{
                        return(poisson(link = "log"))
                    }} else if ("{glm_link}" == "identity") {{
                        return(poisson(link = "identity"))
                    }} else if ("{glm_link}" == "sqrt") {{
                        return(poisson(link = "sqrt"))
                    }}
                }} else if ("{glm_dist_short}" == "nb") {{
                    library(MASS)
                    if ("{glm_link}" == "log") {{
                        return(negative.binomial(theta = {glm_nb_theta if glm_nb_theta else 1}, link = "log"))
                    }} else if ("{glm_link}" == "identity") {{
                        return(negative.binomial(theta = {glm_nb_theta if glm_nb_theta else 1}, link = "identity"))
                    }} else if ("{glm_link}" == "sqrt") {{
                        return(negative.binomial(theta = {glm_nb_theta if glm_nb_theta else 1}, link = "sqrt"))
                    }}
                }}
            }}

            # 並列クラスターの設定
            n_cores <- {n_cores}
            cl <- makeCluster(n_cores)

            # クラスターに必要なパッケージとデータを送信
            clusterEvalQ(cl, {{
                library(MASS)
                if ("{glm_dist_short}" == "beta") {{
                    library(betareg)
                    library(lmtest)
                }}
            }})

            clusterExport(cl, c("counts", "coldata", "get_family"))

            # 処理関数
            process_gene <- function(i) {{
                gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
                family_to_use <- get_family()
                
                result <- tryCatch({{
                    if ("{glm_dist_short}" == "beta") {{
                        # Beta回帰の場合はbetaregパッケージを使用
                        library(betareg)
                        full_fit <- betareg(as.formula("{full_model}"), 
                                          data=gene_data, 
                                          link=family_to_use$link)
                        reduced_fit <- betareg(as.formula("{reduced_model}"), 
                                             data=gene_data, 
                                             link=family_to_use$link)
                        
                        # 尤度比検定
                        library(lmtest)
                        lr_test <- lrtest(reduced_fit, full_fit)
                        
                        # 結果を返す
                        c(statistic = lr_test$Chisq[2],
                          df = lr_test$Df[2],
                          p_value = lr_test$`Pr(>Chisq)`[2],
                          logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
                    }} else {{
                        # 通常のGLMの場合
                        full_fit <- glm(as.formula("{full_model}"), family=family_to_use, data=gene_data)
                        reduced_fit <- glm(as.formula("{reduced_model}"), family=family_to_use, data=gene_data)
                        
                        # 尤度比検定
                        anova_result <- anova(reduced_fit, full_fit, test="LRT")
                        
                        # 結果を返す
                        c(statistic = anova_result$Deviance[2],
                          df = anova_result$Df[2],
                          p_value = anova_result$`Pr(>Chi)`[2],
                          logLik_diff = logLik(full_fit) - logLik(reduced_fit))
                    }}
                }}, error=function(e) {{
                    # エラーの場合はNAを返す
                    if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                    c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
                }})
                
                return(result)
            }}

            # 並列処理の実行
            cat("Starting parallel GLM regression on", n_cores, "cores for", nrow(counts), "genes\\n")
            cat("Using distribution family:", "{glm_dist_family}", "\\n")
            system.time(
                results_list <- parLapply(cl, 1:nrow(counts), process_gene)
            )

            # クラスターの終了
            stopCluster(cl)

            # 結果をデータフレームに変換
            results_matrix <- do.call(rbind, results_list)
            results <- as.data.frame(results_matrix)
            rownames(results) <- rownames(counts)

            # NAの数をカウント
            na_count <- sum(is.na(results$statistic))
            if (na_count > 0) {{
                cat("警告:", na_count, "個の遺伝子でモデルが収束しませんでした (", 
                    round(100 * na_count / nrow(results), 2), "%)\\n")
            }}

            # 多重検定補正
            results$adj.P.Val <- p.adjust(results$p_value, method="BH")

            # モデル情報の保存
            cat("\\n### GLMモデル情報 ###\\n", file='{res_dir}/glm_model_info.txt')
            cat("分布ファミリー:", "{glm_dist_family}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("リンク関数:", "{glm_link}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("フルGLMモデル式:", "{full_model}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("縮小GLMモデル式:", "{reduced_model}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)

            # 結果を保存
            write.table(results[order(results$p_value), ], 
                        file='{res_dir}/glm_{glm_dist_short}_{glm_link}_res.tsv', 
                        sep='\\t', quote=FALSE, col.names=NA)

            cat("GLM regression analysis completed\\n")
            sink()
            """

        elif test_method == 'Generalized Additive Model (GAM)':
            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            ro.r.assign('beta_normalization', beta_normalization)
            ro.r.assign('spline_type', spline_type)
            if dist_short == "nb":
                ro.r.assign('nb_theta', nb_theta)
            ro.r.assign('dist_short', dist_short)
            # Save input to files for R import
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)
            ro.r.assign('temp_dir', temp_dir)


            # 更新されたGAM Rコード
            r_code = f"""
                sink()
                sink(paste0({temp_dir}, "/GAM_output.txt"))
                library(mgcv)
                library(lmtest)
                library(parallel)

                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\\t')
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\\t')

                # 分布ファミリーに応じたデータ前処理
                if ("{dist_short}" == "beta") {{
                    # 0-1の境界の調整
                    eps <- {epsilon}
                    counts <- pmax(pmin(counts, 1-eps), eps)
                }} else if ("{dist_short}" == "gaussian") {{
                    # ガウスの場合は特に前処理不要
                }} else if ("{dist_short}" == "poisson" || "{dist_short}" == "nb") {{
                    # ポアソンとNBはカウントデータを想定（整数化）
                    counts <- round(counts)
                }}

                # 時間変数の特定とフォーミュラの作成
                time_var <- "{full[0]}"
                cat("Time variable:", time_var, "\\n")

                # 時間変数の確認と変換
                if(time_var %in% colnames(coldata)) {{
                  cat("Time variable exists in coldata. Values:", "\\n")
                  print(coldata[[time_var]])
                  
                  # 時間変数が数値でない場合は変換
                  if(!is.numeric(coldata[[time_var]])) {{
                    cat("Converting time variable to numeric\\n")
                    coldata[[time_var]] <- as.numeric(gsub("[^0-9.]", "", as.character(coldata[[time_var]])))
                    cat("After conversion:", "\\n")
                    print(coldata[[time_var]])
                  }}
                }} else {{
                  cat("WARNING: Time variable not found in coldata!\\n")
                }}

                # normalization（ベータ分布以外でも使える）
                if ({beta_normalization} == "TRUE"){{
                    coldata[[time_var]] <- coldata[[time_var]] / max(coldata[[time_var]])
                    cat("Time is normalized by max value.")
                }}

                # GAMモデル式の作成
                gam_full_formula <- "{full_model.replace('~', '')}"
                gam_reduced_formula <- "{reduced_model.replace('~', '')}"

                # 平滑化項の追加 (時間変数に対して)
                if(time_var %in% colnames(coldata)) {{
                  if(length(unique(coldata[[time_var]])) >= 3) {{  # 少なくとも3つの異なる値が必要
                    # フルモデルにのみ平滑化項を追加
                    if(grepl(time_var, gam_full_formula)) {{
                      gam_full_formula <- gsub(
                        paste0("\\\\b", time_var, "\\\\b"), 
                        paste0("s(", time_var, ", k={gam_k}, bs='{spline_type}')"), 
                        gam_full_formula
                      )
                    }} else {{
                      # time_varが明示的に含まれていない場合は追加
                      gam_full_formula <- paste(gam_full_formula, "+", paste0("s(", time_var, ", k={gam_k}, bs='{spline_type}')"))
                    }}
                    cat("Full GAM formula:", gam_full_formula, "\\n")
                    cat("Reduced GAM formula:", gam_reduced_formula, "\\n")
                  }} else {{
                    cat("Not enough unique time points for smoothing, using linear terms\\n")
                  }}
                }}

                cat("{dist_short}")

                # 分布ファミリーの設定関数
                get_family <- function() {{
                    if ("{dist_short}" == "beta") {{
                        return(betar())
                    }} else if ("{dist_short}" == "gaussian") {{
                        return(gaussian())
                    }} else if ("{dist_short}" == "poisson") {{
                        return(poisson())
                    }} else if ("{dist_short}" == "nb") {{
                        library(mgcv)
                        # mgcvのnb関数はthetaパラメータを受け取る
                        return(negbin(theta = {nb_theta}))
                    }}
                }}

                # 並列クラスターの設定
                n_cores <- {n_cores}
                cl <- makeCluster(n_cores)

                # クラスターに必要なパッケージとデータを送信
                clusterEvalQ(cl, {{
                  library(mgcv)
                  library(lmtest)
                }})

                clusterExport(cl, c("counts", "coldata", "gam_full_formula", 
                                    "gam_reduced_formula", "get_family"))

                # 処理関数
                process_gene <- function(i) {{
                  gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
                  family_to_use <- get_family()
                  
                  result <- tryCatch({{
                    # フルモデルと縮小モデルのフィッティング
                    full_fit <- gam(as.formula(paste("y ~", gam_full_formula)), 
                                    family=family_to_use, data=gene_data, method="{gam_method}")
                    
                    reduced_fit <- gam(as.formula(paste("y ~", gam_reduced_formula)), 
                                       family=family_to_use, data=gene_data, method="{gam_method}")
                    
                    # 尤度比検定
                    lr_test <- lrtest(reduced_fit, full_fit)
                    
                    # 結果を返す
                    c(statistic = lr_test$Chisq[2],
                      df = lr_test$Df[2],
                      p_value = lr_test$`Pr(>Chisq)`[2],
                      logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
                  }}, error=function(e) {{
                    # エラーの場合はNAを返す
                    if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                    c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
                  }})
                  
                  return(result)
                }}

                # 並列処理の実行
                cat("Starting parallel GAM regression on", n_cores, "cores for", nrow(counts), "genes\\n")
                cat("Using distribution family:", "{dist_family}", "\\n")
                system.time(
                  results_list <- parLapply(cl, 1:nrow(counts), process_gene)
                )

                # クラスターの終了
                stopCluster(cl)

                # 結果をデータフレームに変換
                results_matrix <- do.call(rbind, results_list)
                results <- as.data.frame(results_matrix)
                rownames(results) <- rownames(counts)

                # モデル情報の保存
                cat("\\n### モデル情報 ###\\n", file='{res_dir}/gam_model_info.txt')
                cat("分布ファミリー:", "{dist_family}", "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("GAM平滑化パラメータ k:", {gam_k}, "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("推定方法:", "{gam_method}", "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("フルGAMモデル式:", paste("y ~", gam_full_formula), "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("縮小GAMモデル式:", paste("y ~", gam_reduced_formula), "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)

                # NAの数をカウント
                na_count <- sum(is.na(results$statistic))
                if (na_count > 0) {{
                  cat("警告:", na_count, "個の遺伝子でモデルが収束しませんでした (", 
                      round(100 * na_count / nrow(results), 2), "%)  k, spline, 不均等時系列の正規化等を検討", 
                      file='{res_dir}/gam_model_info.txt', append=TRUE)
                }}

                # 多重検定補正
                results$adj.P.Val <- p.adjust(results$p_value, method="BH")

                # 結果を保存
                write.table(results[order(results$p_value), ], 
                            file='{res_dir}/gam_{dist_short}_{spline_type}_res.tsv', 
                            sep='\\t', quote=FALSE, col.names=NA)

                cat("GAM regression analysis completed\\n")
                sink()
            """



        elif test_method == 'maSigPro':
            # Restore : in interaction terms (was escaped as \: for display)
            full_model = full_model.replace('\\:', ':')
            reduced_model = reduced_model.replace('\\:', ':')

            # ファイル保存設定
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)


            # 時間情報を含む適切なedesignデータフレームを作成するRコード
            r_code = f"""
            library(maSigPro)
            
            # データの読み込み
            cat("Loading data...\\n")
            counts <- read.table("{counts_file}", header=TRUE, row.names=1, sep="\\t")
            coldata <- read.table("{coldata_file}", header=TRUE, sep="\\t")
            print(coldata)

            # maSigPro用の適切なデザインデータフレームを作成
            cat("Creating proper design matrix for maSigPro...\\n")
            
            # グループ情報を取得（"Group"列）
            time_col <- as.character(coldata${full[0]})
            

            # 時間情報を抽出（例："0w", "1w", "4w"から数値に変換）
            time_values <- as.numeric(gsub("[^0-9.]", "", time_col))
            cat("time_values")
            cat(time_values)
            
            # レプリケート情報を作成
            # 同じ時間値を持つサンプルにユニークな番号を割り当て
            replicates <- numeric(length(time_values))
            for (t in unique(time_values)) {{
                idx <- which(time_values == t)
                replicates[idx] <- 1:length(idx)
            }}
            
            # maSigPro用の正しいedesignデータフレームを作成
            edesign <- data.frame(
                Time = time_values,
                Replicate = replicates
            )
            rownames(edesign) <- colnames(counts)
            
      #      # 他の実験条件があれば追加
      #      if (ncol(coldata) > 1) {{
      #          for (i in 2:ncol(coldata)) {{
      #              col_name <- colnames(coldata)[i]
      #              edesign[[col_name]] <- coldata[[i]]
      #          }}
      #      }}

            # Add Group column (all 1s for single condition)
            edesign$Group <- rep(1, nrow(edesign)) #tutorialに合わせて全部1にする
            
            # データ型の確認
            cat("Time values:", paste(time_values, collapse=", "), "\\n")
            cat("Time values are numeric:", is.numeric(edesign$Time), "\\n")
            cat("edesign")
            print(edesign)
            
            # 前処理
            """
            
            # データタイプに応じた処理の追加
            if data_type == "0-1 data (logit transformation)":
                r_code += f"""
            # 0-1データの処理
            eps <- {epsilon}
            counts <- pmax(pmin(counts, 1-eps), eps)
            counts <- log(counts/(1-counts))
            use_counts_param <- FALSE
            cat("Applied logit transformation\\n")
            """
            elif data_type == "qPCR/continuous data (Gaussian)":
                r_code += f"""
            # qPCR/連続データの処理
            use_counts_param <- FALSE
            cat("Using Gaussian model for continuous data\\n")
            
            # データの前処理
            original_counts <- counts
            """
                
                # log変換オプション
                if log_transform:
                    r_code += """
            # Log2変換を適用
            counts <- log2(counts + 1)  # +1を加えて0値を避ける
            cat("Applied log2(x+1) transformation\\n")
            """
                
                # 正規化オプション
                if normalization:
                    r_code += """
            # Z-score正規化（サンプル間）
            counts <- t(scale(t(counts)))
            cat("Applied z-score normalization across samples\\n")
            """
                
                r_code += """
            cat("qPCR data preprocessing completed\\n")
            """
            else:  # RNA-seq count data
                r_code += """
            # RNA-seqデータの処理
            use_counts_param <- TRUE
            cat("Using GLM for count data\\n")
            """
            
            # デザイン行列を使用した分析
            r_code += f"""
            # デザイン行列を使用して分析
            cat("Running maSigPro analysis...\\n")
            
            # 指定された次数でデザイン行列を作成
            design <- make.design.matrix(edesign, degree={degree})
            
            # 回帰分析の実行
            cat("Running p.vector...\\n")
            fit <- p.vector(counts, design$edesign, Q={q_value}, MT.adjust="none", counts=use_counts_param)
           # fit <- p.vector(counts, design$edesign, Q={q_value}, MT.adjust="BH", counts=use_counts_param)
            
            # 有意な遺伝子数の確認
            sig_count <- sum(fit$p < {q_value}, na.rm=TRUE)
            cat("Genes with p <", {q_value}, ":", sig_count, "\\n")

            # After running p.vector() and finding no significant genes

            
            # 有意な遺伝子がある場合のみ次のステップへ
            if (sig_count > 0) {{
                cat("Running T.fit...\\n")
                tstep <- T.fit(fit, step.method="backward", alfa={q_value})
                
                cat("Getting significant genes...\\n")
                sigs <- get.siggenes(tstep, rsq={rsq}, vars="each")
                
                # 結果を保存
                if (!is.null(sigs) && !is.null(sigs$sig.genes) && !is.null(sigs$sig.genes$sig.profiles) && nrow(sigs$sig.genes$sig.profiles) > 0) {{
                    cat("Found", nrow(sigs$sig.genes$sig.profiles), "significant genes\\n")
                    
                    # プロファイルと係数を保存
                    write.table(sigs$sig.genes$sig.profiles, file="{res_dir}/maSigPro_sig_profiles.tsv", sep="\\t", quote=FALSE)
                    write.table(sigs$coefficients, file="{res_dir}/maSigPro_coefficients.tsv", sep="\\t", quote=FALSE)
                    
                    # 要約情報の保存
                    cat("maSigPro Analysis Results\\n",
                        "------------------------\\n",
                        "Total genes analyzed: ", nrow(counts), "\\n",
                        "Significant genes (p <", {q_value}, "): ", sig_count, "\\n",
                        "Significant genes (rsq >", {rsq}, "): ", nrow(sigs$sig.genes$sig.profiles), "\\n",
                        file="{res_dir}/maSigPro_summary.txt")
                }} else {{
                    cat("No genes passed R-squared threshold\\n")
                    cat("No genes passed R-squared threshold of", {rsq}, "\\n", file="{res_dir}/maSigPro_summary.txt")
                }}
            }} else {{
                cat("No significant genes found\\n")
                cat("No significant genes found at Q-value", {q_value}, "\\n", file="{res_dir}/maSigPro_summary.txt")
            }}
            """
            
            # Rコードを実行
            with st.spinner('Calculating maSigPro... This may take a while.'):
                try:
                    # Rコードをデバッグ用に保存
                    with open(os.path.join(temp_dir, 'debug_maSigPro.R'), 'w') as f:
                        f.write(r_code)
                    
                    # Rコード実行
                    ro.r(r_code)
                    
                    # 結果用の空データフレームを作成（エラー回避用）
                    res_df = pd.DataFrame()
                    
                    # 結果の表示
                    st.markdown("### maSigPro Analysis Results")
                    
                    # 結果ファイルの確認と表示
                    summary_file = os.path.join(res_dir, 'maSigPro_summary.txt')
                    if os.path.exists(summary_file):
                        with open(summary_file, 'r') as f:
                            summary = f.read()
                        st.text(summary)
                    
                    # 有意な遺伝子の結果を表示
                    sig_profiles_file = os.path.join(res_dir, 'maSigPro_sig_profiles.tsv')
                    if os.path.exists(sig_profiles_file):
                        res_df = pd.read_csv(sig_profiles_file, sep='\t', index_col=0)
                        if not res_df.empty:
                            st.write("### Top significant genes:")
                            st.dataframe(res_df.head(10))
                            
                            # 係数の表示
                            coef_file = os.path.join(res_dir, 'maSigPro_coefficients.tsv')
                            if os.path.exists(coef_file):
                                st.write("### Regression coefficients:")
                                coef = pd.read_csv(coef_file, sep='\t', index_col=0)
                                st.dataframe(coef.head(10))
                            
                            # 結果のダウンロード
                            data_type_short = ""
                            if data_type == "RNA-seq count data (GLM)":
                                data_type_short = "RNAseq"
                            elif data_type == "qPCR/continuous data (Gaussian)":
                                data_type_short = "qPCR"
                            elif data_type == "0-1 data (logit transformation)":
                                data_type_short = "logit"
                            
                            file_name = file_name_head + f"_maSigPro_{data_type_short}"
                            shutil.make_archive("res", format='zip', root_dir=res_dir)
                    else:
                        st.info("No significant genes were found. Try adjusting the Q-value or R-squared threshold.")
                
                except Exception as e:
                    st.error(f"Error executing R code: {str(e)}")
                    # エラー時にも空のDataFrameを作成
                    res_df = pd.DataFrame()
            
     
            
            # 結果のZIP生成
            data_type_short = ""
            if data_type == "RNA-seq count data (GLM)":
                data_type_short = "RNAseq"
            elif data_type == "qPCR/continuous data (Gaussian)":
                data_type_short = "qPCR"
            elif data_type == "0-1 data (logit transformation)":
                data_type_short = "logit"
            
            file_name = file_name_head + f"_maSigPro_{data_type_short}"
            shutil.make_archive("res", format='zip', root_dir=res_dir)

        # 結果の表示と保存
        if test_method == 'Beta Regression':
            ro.r(r_code)
            res_df = pd.read_csv(os.path.join(res_dir, 'betareg_res.tsv'), sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)

            # モデル収束情報の確認と表示
            convergence_file = os.path.join(res_dir, 'model_convergence_info.txt')
            if os.path.exists(convergence_file):
                with open(convergence_file, 'r') as f:
                    convergence_info = f.read()
        
            # 収束に問題があるかチェック
            if "警告" in convergence_info:
                st.warning(convergence_info)
            else:
                st.success(convergence_info)
                
            file_name = file_name_head + "_betareg_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")
            shutil.make_archive("res", format='zip',root_dir= res_dir)

        elif test_method == 'Generalized Linear Model (GLM)':
            ro.r(r_code)
            res_df = pd.read_csv(os.path.join(res_dir, f'glm_{glm_dist_short}_{glm_link}_res.tsv'), sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)
            
            # モデル情報の表示
            model_info_file = os.path.join(res_dir, 'glm_model_info.txt')
            if os.path.exists(model_info_file):
                with open(model_info_file, 'r') as f:
                    model_info = f.read()
                st.text(model_info)
            
            file_name = file_name_head + f"_glm_{glm_dist_short}_{glm_link}_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")
            shutil.make_archive("res", format='zip',root_dir= res_dir)

        elif test_method == 'Generalized Additive Model (GAM)':
            ro.r(r_code)
            res_df = pd.read_csv(os.path.join(res_dir, f'gam_{dist_short}_{spline_type}_res.tsv'), sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)
            
            file_name = file_name_head + "_gam_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")
            shutil.make_archive("res", format='zip',root_dir= res_dir)


            # モデル収束情報の確認と表示
            convergence_file = os.path.join(res_dir, 'gam_model_info.txt')
            if os.path.exists(convergence_file):
                with open(convergence_file, 'r') as f:
                    convergence_info = f.read()
        
            # 収束に問題があるかチェック
            if "警告" in convergence_info:
                st.warning(convergence_info)
            else:
                st.success(convergence_info)

        if res_df is not None:
            with open("res.zip", "rb") as fp:
                btn = st.download_button(
                    label="Download Results",
                data=fp,
                file_name=file_name + "_DESeq2-LRT.zip",
                mime = "zip"
                )
            try:
                os.remove(file_name + "_DESeq2-LRT.zip")
                shutil.rmtree(temp_dir)
                os.mkdir(temp_dir)
            except:
                pass


#　データを送る前にすべてゼロのデータは除くべき


# refが指定されているときはファイル名を調整する?
