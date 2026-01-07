#!!!!!!!!!!!!!! pip install rpy2==3.5.1  newshiiba-jiyonisErrorisoutru

# basemainaltoglobalchangenumwithCalculationdo。
# pythonfromassgnsareruofisglobalchangenum


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
    # mostmoshortitextcharcoloflongsathegetget
    min_length = min(len(s) for s in strings)
    # copassofendtailpartdivoflongsatheviewtsukeru
    suffix_length = 0
    for i in range(1, min_length + 1):
        suffix = strings[0][-i:]
        if all(s.endswith(suffix) for s in strings):
            suffix_length = i
        else:
            break            
    # copassofendtailpartdivisviewtsufromnotplacematchissourceofrisutothereturnsu
    if suffix_length == 0:
        return strings        
    # copassofendtailpartdivthedeleteremoveshitenewshiirisutothemakebecome
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
f = ro.r("source('pages/deseq2_func.R')") # full pathisrequired

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


# tempintoSavedo
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    #oldidirecotryandFilethedeleteremovedo
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
                         ["DESeq2-LRT", "limma eBayes", "Beta Regression", 
                          "Generalized Linear Model (GLM)", "Generalized Additive Model (GAM)", "maSigPro"], 
                         index=0, label_visibility = 'collapsed')

    st.markdown("##### limma eBayes with logit transformation and beta regression are for proportion data.")

    if test_method == 'DESeq2-LRT':
        st.markdown("---")
        st.markdown("### DESeq2 Polynomial Options:")
        add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（timebetweenpassovertoaccompunonlineshapechangeize）thecheckoutdofortomanyitemformatofitemtheaddemasu")

        
        if add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2nextmanyitemformat：1tsuofwaydirchangeize（increaseadd→reducefewetc）thecheckout。3nextmanyitemformat：2tsuofwaydirchangeize（increaseadd→reducefew→increaseaddetc）thecheckout")
            
            # manyitemformatrealinstallwaymethodofSelect
            poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()relnum：directcrossmanyitemformatmouseusepossiblewithinferrec。I()relnum：simplepurenabekimultitemofuseuse")
            
            # poly()relnumtheuseuplacematchofaddaddOption
            use_poly_function = poly_implementation.startswith("poly()")
            
            if use_poly_function:
                # manyitemformattaipuofSelectOption
                poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="directcrossmanyitemformat（Orthogonal）：colinenaturetheavoidkerufortoinferrec。Raw：solveinterpshiyasuirelatenumthegetruiscolinenatureofquesttopicisexistplacematchari")
                use_raw = poly_type.startswith("Raw")
            else:
                use_raw = False  # I()relnumuseusetimeisrelrelatenotofwithFalse
                            
            polynomial_degree = 2 if polynomial_term.startswith('2:') else 3
            
        else:
            polynomial_variable = None
            polynomial_degree = 1
            use_raw = False
            use_poly_function = True  # defuorutoval


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
            
        st.markdown("---")
        st.markdown("### limma Polynomial Options:")
        limma_add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（timebetweenpassovertoaccompunonlineshapechangeize）thecheckoutdofortomanyitemformatofitemtheaddemasu")

        if limma_add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            limma_polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2nextmanyitemformat：1tsuofwaydirchangeize（increaseadd→reducefewetc）thecheckout。3nextmanyitemformat：2tsuofwaydirchangeize（increaseadd→reducefew→increaseaddetc）thecheckout")
            
            # manyitemformatrealinstallwaymethodofSelect
            limma_poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()relnum：directcrossmanyitemformatmouseusepossiblewithinferrec。I()relnum：simplepurenabekimultitemofuseuse")
            
            # poly()relnumtheuseuplacematchofaddaddOption
            limma_use_poly_function = limma_poly_implementation.startswith("poly()")
            
            if limma_use_poly_function:
                # manyitemformattaipuofSelectOption
                limma_poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="directcrossmanyitemformat（Orthogonal）：colinenaturetheavoidkerufortoinferrec。Raw：solveinterpshiyasuirelatenumthegetruiscolinenatureofquesttopicisexistplacematchari")
                limma_use_raw = limma_poly_type.startswith("Raw")
            else:
                limma_use_raw = False  # I()relnumuseusetimeisrelrelatenotofwithFalse
                            
            limma_polynomial_degree = 2 if limma_polynomial_term.startswith('2:') else 3
            
        else:
            limma_polynomial_variable = None
            limma_polynomial_degree = 1
            limma_use_raw = False
            limma_use_poly_function = True  # defuorutoval
            
    # be-tatimereturnspechaveofOption
    if test_method == 'Beta Regression':
        st.markdown("### Beta Regression Options:")
        epsilon = st.number_input("Epsilon for boundary adjustment (0-1 data)", 
                                min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        add_higher = st.checkbox("Add polynomial terms?", value=False, help="nonlineshapealnachangeizethecaptureeruformanyitemitemtheadderu")
        beta_polynomial_degree = 1
        if add_higher:
            polynomial_term = st.radio("Degree", ['2:Quadratic term','3:Cubic term'], index = 0, help = "2nextofitemtheadderuandU-shaped/inverted U-shaped patternsiscaptureerareru。3nextofitemtheaddadddokoandwith、thanmultimiscnaExpressionpata-nthecaptureeru。Exampleeba：急upascaftertohorizbaitonari、soofafterlowbelowdopata-n、波shapepata-n。")
            st.markdown("#### The first item in full model will be used for polynominal term.")
            if polynomial_term == "2:Quadratic term":
                beta_polynomial_degree = 2
            else:
                beta_polynomial_degree = 3

        
    # GLMspechaveofOption
    if test_method == 'Generalized Linear Model (GLM)':
        st.markdown("### GLM Options:")
        # divdistfuamiri-ofSelectOption
        glm_dist_family = st.radio("Probability distribution", 
                              ["Beta (0-1)", 
                               "Gaussian", 
                               "Poisson", 
                               "Negative Binomial"],
                              index=0,
                              help="DatataiputorespondjitacertainratedivdisttheSelectshitekudasai。Beta: 0-1ofval（ratiomatchetc）, Gaussian: connectcontinuevalData, Poisson: shinpurunakauntoData, Negative Binomial: overdivscatterofexistkauntoData（RNA-seq, scRNA-seqeq）")

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

        # rinkurelnumofSelect
        if glm_dist_family == "Beta (0-1)":
            glm_link = st.radio("Link function for Beta distribution",
                               ["logit", "probit", "cloglog"],
                               index=0,
                               help="logit: mostmoonegenal、probit: correctruledivdistbe-su、cloglog: 極valdivdistbe-su")
            glm_epsilon = st.number_input("Epsilon for boundary adjustment", 
                                   min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        elif glm_dist_family == "Gaussian":
            glm_link = st.radio("Link function for Gaussian distribution",
                               ["identity", "log", "inverse"],
                               index=0,
                               help="identity: lineshaperelrelate、log: correctofvalofmi、inverse: invnumchangechange")
        elif glm_dist_family == "Poisson":
            glm_link = st.radio("Link function for Poisson distribution",
                               ["log", "identity", "sqrt"],
                               index=0,
                               help="log: mostmoonegenal（kauntoData）、identity: lineshape、sqrt: flatwayrootchangechange")
        elif glm_dist_family == "Negative Binomial":
            glm_link = st.radio("Link function for Negative Binomial distribution",
                               ["log", "identity", "sqrt"],
                               index=0,
                               help="log: mostmoonegenal（overdivscatterkauntoData）、identity: lineshape、sqrt: flatwayrootchangechange")
            glm_nb_theta = st.number_input("overdivscatterParameter (theta)", 
                                     min_value=0.1, max_value=100.0, value=10.0,
                                     help="valissmallsaihodooverdivscatterisbigkiikoandthemeanmeanshimasu。RNA-seqwithispassnormal5-10extentdegree。scRNA:0.5-3 (10x:1-2)")

        # Polynomial options for GLM
        st.markdown("---")
        st.markdown("### GLM Polynomial Options:")
        glm_add_polynomial = st.checkbox("Add polynomial terms for time variable?", value=False, 
                                   help="Non-linear changes over time（timebetweenpassovertoaccompunonlineshapechangeize）thecheckoutdofortomanyitemformatofitemtheaddemasu")

        if glm_add_polynomial:
            st.markdown("#### Using first variable from full model as time variable, which must be numeric.")

            glm_polynomial_term = st.radio("Polynomial degree", 
                                    ['2: Quadratic: U-shaped pattern', 
                                     '3: Cubic: S-shaped pattern'], 
                                    index=0, 
                                    help="2nextmanyitemformat：1tsuofwaydirchangeize（increaseadd→reducefewetc）thecheckout。3nextmanyitemformat：2tsuofwaydirchangeize（increaseadd→reducefew→increaseaddetc）thecheckout")
            
            # manyitemformatrealinstallwaymethodofSelect
            glm_poly_implementation = st.radio("Implementation method",
                                        ["poly() function", "I() function (explicit powers)"],
                                        index=0,
                                        help="poly()relnum：directcrossmanyitemformatmouseusepossiblewithinferrec。I()relnum：simplepurenabekimultitemofuseuse")
            
            # poly()relnumtheuseuplacematchofaddaddOption
            glm_use_poly_function = glm_poly_implementation.startswith("poly()")
            
            if glm_use_poly_function:
                # manyitemformattaipuofSelectOption
                glm_poly_type = st.radio("Polynomial type",
                                   ["Orthogonal", "Raw"],
                                   index=0,
                                   help="directcrossmanyitemformat（Orthogonal）：colinenaturetheavoidkerufortoinferrec。Raw：solveinterpshiyasuirelatenumthegetruiscolinenatureofquesttopicisexistplacematchari")
                glm_use_raw = glm_poly_type.startswith("Raw")
            else:
                glm_use_raw = False  # I()relnumuseusetimeisrelrelatenotofwithFalse
                            
            glm_polynomial_degree = 2 if glm_polynomial_term.startswith('2:') else 3
            
        else:
            glm_polynomial_variable = None
            glm_polynomial_degree = 1
            glm_use_raw = False
            glm_use_poly_function = True  # defuorutoval

    # GAMspechaveofOption
    if test_method == 'Generalized Additive Model (GAM)':
        st.markdown("### GAM Options:")
        # divdistfuamiri-ofSelectOptiontheaddadd
        dist_family = st.radio("Probability distribution", 
                              ["Beta (0-1)", 
                               "Gaussian", 
                               "Poisson", 
                               "Negative Binomial"],
                              index=0,
                              help="DatataiputorespondjitacertainratedivdisttheSelectshitekudasai。Beta: 0-1ofval（ratiomatchetc）, Gaussian: connectcontinuevalData, Poisson: shinpurunakauntoData, Negative Binomial: overdivscatterofexistkauntoData（RNA-seq, scRNA-seqeq）")

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


        # divdisttorespondjitaParameterSettings
        if dist_family == "Beta (0-1)":
            epsilon = st.number_input("Epsilon for boundary adjustment", 
                                   min_value=0.0000001, max_value=0.01, value=0.000001, format="%.7f")
        elif dist_family == "Negative Binomial":
            nb_theta = st.number_input("overdivscatterParameter (theta)", 
                                     min_value=0.1, max_value=100.0, value=10.0,
                                     help="valissmallsaihodooverdivscatterisbigkiikoandthemeanmeanshimasu。RNA-seqwithispassnormal5-10extentdegree。scRNA:0.5-3 (10x:1-2)")


        gam_k = st.slider("Spline basis dimension (k)", min_value=3, max_value=20, value=4, help = "nonlineshaperelrelatethemoderuizedotasumu-jingurelnumofflexsoftnature（multimiscsa）thecontrolctrldoParameter。k ofvalisbigkiihodo、moderuisthanmultimiscnanonlineshapepata-nthecaptureeru。「timebetweenpointofnum + 0〜1」")
        gam_method = st.radio("Smoothing parameter estimation method", 
                            ["REML", "GCV.Cp", "ML"], index=0,
                            help="REML (Restricted Maximum Likelihood):divscatterParameterandflatsmoothizeParameterthesametimetoEstimation。baiasuisfewnaku、thantrustrelynatureofhighiEstimation。smallsanaSamplesaizuwithmoComparisonalsafeset。   ML (Maximum Likelihood):mostlikelyEstimationmethod。moderuComparison（AIC, BICetc）tofitsu。smallsanaSamplesaizuwithisbaiasuisgenjirupossiblenature。  GCV.Cp (Generalized Cross Validation / Mallows' Cp):kurosubaride-shiyontobaseduku。moderuofadvmeasablepowerthemostfitize。advmeasisidxalofplacematchtohaveuse。")

        selected_spline = st.radio("Spline type",
                        ['Thin Plate Regression Splines (tp)', 'Cubic Regression Splines (cr)', 'Cubic Smoothing Splines (cs)'],
                        index =1, help="tp:mostmogenericusealnasupuraintaipu。onegenaltodefuoruto。  cr:3nextmanyitemformatofareadivalnasetmimatchwase。timebetweenpointisfewnotplacematchtofitshiteexistke-suisexist。  cs:observemeasDatapointtonotothedistplace。nonnormaltoflexsoftwith、Datapointbetweenthesmoothrakatosuppbetween。Datapointisfewnotplacematchtohaveusenakoandisexist")
        # SelecttobaseduitechangenumtheSettings
        if selected_spline == 'Thin Plate Regression Splines (tp)':
            spline_type = "tp"
        elif selected_spline == 'Cubic Regression Splines (cr)':
            spline_type = "cr"
        elif selected_spline == 'Cubic Smoothing Splines (cs)':
            spline_type = "cs"

    #    beta_norm = st.checkbox("Normalization by maximum value of time variable", value = False,
    #        help='maximumvaltoyoruNormalization.timebetweentheincludemuandkitocollbundleshinotplacematchtotrymiru') #effectresultisnasasou

#        beta_normalization = "TRUE" if beta_norm else "FALSE"
        beta_normalization = "FALSE"
    if test_method in ['Beta Regression', 'Generalized Linear Model (GLM)', 'Generalized Additive Model (GAM)']:
        n_cores = st.slider("Parallel cores", min_value=1, 
                           max_value=os.cpu_count()-1, 
                           value=max(1, os.cpu_count()//2-4))


    if test_method == 'maSigPro':
        st.markdown("### maSigPro Options:")
        
        # DatataipuofSelect
        data_type = st.radio(
            "Data type:",
            ["RNA-seq count data (GLM)", "qPCR/continuous data (Gaussian)", "0-1 data (logit transformation)"],
            index=0
        )
        
        # DatataiputorespondjitaParameterSettings
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
                                      help="qPCRDataofplacematch、passnormallog2changechangetherowimasu（ΔCtvaletc）")
            normalization = st.checkbox("Z-score normalization across samples", value=False,
                                       help="qPCRwithispassnormalnotneed")
        
        # copassofParameter
        degree = st.slider("Polynomial degree", min_value=1, max_value=3, value=2)
        rsq = st.number_input("R-squared cutoff", min_value=0.1, max_value=0.9, value=0.7, step=0.05)
        q_value = st.number_input("Q-value (FDR)", min_value=0.001, max_value=0.5, value=0.05, step=0.01)
        
        # kurasutaringuOption
        st.markdown("##### Clustering options:")
        cluster_method = st.radio("Clustering method", 
                                ["hclust", "kmeans", "clara"], index=0)
        k = st.slider("Number of clusters", min_value=2, max_value=15, value=9)
        
        # VisualizationOption
        st.markdown("##### Visualization:")
        plot_top_n = st.slider("Number of genes to plot", min_value=5, max_value=50, value=20)

    st.markdown("---")

st.markdown("### DESeq2 likelihood-ratio test (LRT), limma eBayes, betareg, GAM for time-course and ANOVA-like test")
st.markdown("### maSigPro for time-course test")
st.markdown("##### DESeq2-LRT, beta regression can use polynomial terms that help time-course analysis")
st.markdown("##### limma eBayes and GAM can be used with both count and non-count data, including AUC")
st.markdown("##### beta regression, GAM with beta regression and limma with logit transformation are for proportion (0-1) data")
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
        if "Row_name" in df.columns.to_list(): # Row_nametheincludemuandki
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
                # colnamesofchangechange
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
                 # colnamesofchangechange
                search_word = '([^\ \(]*)\ \(.*'

                for i in range(1, len(content)):
                    match = re.search(search_word, content[i])
                    if match:
                        content[i] = match.group(1).replace(' ', '_')
                df.columns = content # oneoncenamebeforethechangefurther
                df['Annotation/Divergence'] = df['Annotation/Divergence'].astype(str) # excel pairrespond
                pattern = "([^|]*)"
                repatter = re.compile(pattern)
                f_annotation = lambda x: repatter.match(x).group(1)
                df.loc[:,'Annotation/Divergence'] = df.loc[:,'Annotation/Divergence'].apply(f_annotation)
                # annotation/divergencebeforetheremoveku
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

############ samplenameofRwithuseenottextcharthemodcorrect
    def make_valid_r_names(names):
        """Rwithvalidnachangenumnametochangechangedorelnum"""
        valid_names = []
        changes_made = False
        
        for name in names:
            original_name = name
            
            # 1. specspecialtextchartheplacechange
            name = re.sub(r'[ ]+', '.', name)  # supe-suthe.to
            name = re.sub(r'[-]+', '_', name)  # haifunthe_to
            name = re.sub(r'[^\w.]', '_', name)  # 英numchar,piriodo,anda-sukoaoroutthe_to
            
            # 2. firstheadisnumcharofplacematchisXtheattachkeru
            if re.match(r'^\d', name):
                name = 'X' + name
            
            # 3. firstheadispiriodowithnumchariscontinuekuplacematchisXtheattachkeru
            if re.match(r'^\.\d', name):
                name = 'X' + name
            
            # 4. advaboutwordchieku（basemainalnamoof）
            r_reserved = ['if', 'else', 'repeat', 'while', 'function', 'for', 'in', 'next', 'break', 
                         'TRUE', 'FALSE', 'NULL', 'Inf', 'NaN', 'NA', 'NA_integer_', 'NA_real_', 
                         'NA_complex_', 'NA_character_']
            if name in r_reserved:
                name = name + '_'
            
            if original_name != name:
                changes_made = True
            
            valid_names.append(name)
        
        return valid_names, changes_made
    
    # Samplenamethemodcorrect
    new_columns, changes_made = make_valid_r_names(df.columns.tolist())
    if changes_made:
        st.write("Sample names have been converted to be R-compatible:")
        for old, new in zip(df.columns.tolist(), new_columns):
            if old != new:
                st.write(f"  '{old}' → '{new}'")
        df.columns = new_columns
############


    st.write('Original gene number:  ' + str(len(df)))

    # floattochangechange errcast悟in
    df = df.astype(float)

    if not float.is_integer(df.iloc[:,0].sum()*1000):
        if test_method == "DESeq2-LRT":
            st.markdown("## It is likely that your data are normalized. Please upload unnormalized raw count data.")

    if test_method == "DESeq2-LRT": #DESeq2isarrangenumize
        df = df.round(0)

    df = df.loc[~(df==0).all(axis=1)] #all0ofrowtheremoveku

########## excelpairrespond?
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
############ samplenameto-isexistplacematchisunderscorehe RwithErrortobecome
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

        if test_method == 'DESeq2-LRT' or (test_method == 'limma eBayes' and limma_count):
            st.markdown("### Filter out weakly-expressed genes before multiple test correction:",help = "independentFiltering default:TRUE flatavgNormalizationkauntotobaseduiteGenetheFilteringshi、manyweightchecksetsuppcorrectofneg担thereducerasukoandwithStatisticalalcheckOutputthedirupsaseru")
            independentFiltering = st.checkbox('Yes', value=True)
            st.markdown("---")

        st.markdown("### Furhter filtering of genes")
        st.markdown("""#### lowExpressionGeneofremoveoutisFDRofCalculationthereform善do""")

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

    if any(df.sum() <= sample_threshold): # count 0ofcoltheremoveku
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


    condition = [str(i) for i in df.columns.tolist()[:]] #errorpreventstop
    group_condition = remove_common_suffix(condition) #endtailofcopassneedelemtheremoveku
  #  group_condition = [remove_after_space(x) for x in condition] #supe-suordesctheremoveku
    group_condition = [remove_sample_num(x) for x in group_condition] #endtailofnumchartheremoveku


    st.markdown("##### Add conditions other than group, such as genotype (comma, space, CR separated):")
    genes = st.text_input("genes",label_visibility = 'collapsed')
    gene_list = []
    if len(genes) > 0:
        gene_list = genes.split(' ') #mazuemptywhitewithdivsep
        gene_list = list(filter(lambda a: a != '', gene_list)) #emptywhiteofmitheremoveku
        if ',' in genes:
            gene_list = sum([x.split(',') for x in gene_list],[]) #sumwithflatflatize sum(x, [])
        if '\t' in genes:
            gene_list = sum([x.split('\t') for x in gene_list],[])
        if '\n' in genes:
            gene_list = sum([x.split('\n') for x in gene_list],[])
        gene_list = [a for a in gene_list if a != ''] #emptytheremoveku
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

    # eachchangenumoftypetheSelectdoforofsekushiyon
    st.write('Select variable types:')
    var_types = {}
    cols = st.columns(len(condition_col))
    for i, col in enumerate(condition_col):
        with cols[i]:
            var_types[col] = st.radio(
                f"{col}",
                options=["categorical", "continuous"],
                index=0,  # defuorutoiskategorikaru
                key=f"type_{col}"
            )
   
# modelthemakeruforofneedelemtherisutotodo
    comb = [':'.join(x) for x in  list(itertools.combinations(condition_col, 2))]
#kokowith ':'.joint(x)anddoand、koofaand、moderuthemakeruandkito:ordescisremovekareru
    selections = selections = sum([condition_col, comb],[])

    null_model = st.checkbox("Null model as reduced model?", value = False, help="returnnomoderuthereduced modeltodo。tsumarifull modeltoSettingsshitaneedcauseofizureka／alltorelconnectshitaExpressionchangemovethecheckout。")

    st.markdown("##### Select conditions for full model:")
    full = st.multiselect('fullmodel',selections, label_visibility = 'collapsed')

    # mostfirstofchangenumthetimebetweenchangenumandshiteSave
    time_var = None
    if len(full) > 0:
        time_var = full[0]
    #  st.write(time_var)

    if not null_model:
        st.markdown("##### Select conditions for reduced model:")
        reduced = st.multiselect('reducedmodel',selections, label_visibility = 'collapsed')
    else:
        reduced = []
    # manyitemformatfituseofrojiku（DESeq2-LRT、limma、GLMwithpolynomialisTrueofandki）
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
            # poly()relnumtheusetarealinstall
            if poly_degree == 2:
                # raw=TRUEispointsetsareteexistplacematchisaddadd
                raw_param = ", raw=TRUE" if poly_use_raw else ""
                full[0] = f"poly({time_var}, degree=2{raw_param})"
            else:  # 3nextmanyitemformat
                raw_param = ", raw=TRUE" if poly_use_raw else ""
                full[0] = f"poly({time_var}, degree=3{raw_param})"
            
            st.markdown(f"##### Using {'raw' if poly_use_raw else 'orthogonal'} polynomial term for {time_var}: {full[0]}")
        else:
            # I()relnumtheusetarealinstall（clearshowalnabekimult）
            if poly_degree == 2:
                # sourceoftimebetweenchangenumtheplacekichangee、2nextofitemtheaddadd
                new_terms = [time_var, f"I({time_var}^2)"]
                # fullofmostfirstofneedelemtheplacekichangee
                full[0] = new_terms[0]
                # 2nextofitemthe挿in
                full.insert(1, new_terms[1])
                
                st.markdown(f"##### Using explicit powers for {time_var}: {time_var} + I({time_var}^2)")
            else:  # 3nextmanyitemformat
                # sourceoftimebetweenchangenumtheplacekichangee、2nextand3nextofitemtheaddadd
                new_terms = [time_var, f"I({time_var}^2)", f"I({time_var}^3)"]
                # fullofmostfirstofneedelemtheplacekichangee
                full[0] = new_terms[0]
                # 2nextand3nextofitemthe挿in
                full.insert(1, new_terms[1])
                full.insert(2, new_terms[2])
                
                st.markdown(f"##### Using explicit powers for {time_var}: {time_var} + I({time_var}^2) + I({time_var}^3)")


    full = [x.replace(':','\:') for x in full] # :ofmamadaandtextcharcolisdeletefaildo
    reduced = [x.replace(':','\:') for x in reduced]

    if len(reduced) > 0 and not null_model:
        null_model = False

    full_model = "~ " + " + ".join(full)
    if null_model:
        reduced_model = "~ 1"
    elif len(reduced) == 0: #reducedthepointsetshiteinotandkiisnull modeltodo
        reduced_model = "~ 1"
        st.markdown("#### Null model is uses as reduced model.")
    else:
        reduced_model = "~ " + " + ".join(reduced)
    st.markdown("##### Full model:  " + full_model)
    st.markdown("##### Reduced model:  " + reduced_model)
    st.markdown("""
Full modelandReduced modelandofdiffiischecksetsareru。
ExampleebagenotypandtimesyscolofDataofandkitogenptypeisrelrelatenaku、timesyscolchangeizethedoGenethecheckoutdoplacematchis
~ genotype + time and ~ genotype ofComparisonandbecome。\n
moshi、genotypespecdiffalwithtimebetweenwithchangeizedoGenethecheckoutdoplacematchis
~ genotype + time + genotype\:time and ~ genotype + time ofComparisonandbecome。\n
\n
Reduced modeltonull modeltheSettingsdoandFull modelofneedcausewithchangeizedoGenethecheckoutdo。\n
ExampleebaWTofCelloftimesyscolDatadakeofplacematch、timetheFull modeltoNull modeltheReduced modeltodo。
        """)

    if (len(condition) != len(df.columns)):
            st.write("The number of group name does not match the data.")

#    df_condition = pd.DataFrame(condition)
#    df_batch = pd.DataFrame(batch)

# 1-Maretcoferrchangechangeheofpairrespond
    check_excel_autoconversion(df)

    if len(df.index.values) != len(set(df.index.values)):
#        st.markdown("#### There are duplicated rows. Converting the names...")
#        st.write("The gene name of the second occurrence has _2 at the end.")
#        lis = df.index.values
#        df.index = [x + ['', '_2'][x in lis[0:i]] for i, x in enumerate(lis)]
        df = rename_duplicates(df)


    # manyitemformatnextnumis1thanbigkiiplacematchoftimebetweenchangenumchieku
    if (test_method == 'Beta Regression' and 'beta_polynomial_degree' in locals() and beta_polynomial_degree > 1) or (test_method == 'DESeq2-LRT' and 'add_polynomial' in locals() and add_polynomial) or (test_method == 'limma eBayes' and 'limma_add_polynomial' in locals() and limma_add_polynomial) or (test_method == 'Generalized Linear Model (GLM)' and 'glm_add_polynomial' in locals() and glm_add_polynomial):
        if len(full) > 0:
          #  st.write("Using polynomial")
            try:
                coldata_file = os.path.join(temp_dir, 'coldata.tsv')
                df_e.to_csv(coldata_file, sep='\t', index=False)
                coldata = pd.read_table(coldata_file)
                if time_var in coldata.columns:
                    # timebetweenchangenumoftypechieku
                    time_col = coldata[time_var]
                    is_numeric = pd.api.types.is_numeric_dtype(time_col)
                    
                    if not is_numeric:
                        # numvaltochangechangepossiblekachieku
                        try:
                            # numchardaketheextractoutdocorrectruletablepresentpata-n
                            numeric_values = time_col.str.extract(r'(\d+\.?\d*)')[0].astype(float)
                            st.info(f"infoinfo: timebetweenchangenum '{time_var}' istextcharcolwithsuis、numvalandshiteextractoutwithkimasu。Analysistimetoselfmovealtochangechangesaremasu。")
                            coldata[time_var] = numeric_values #coldatathechangeeteoku
                            
                            # yuni-kupointonumchieku
                            unique_points = len(numeric_values.unique())
                            
                            # manyitemformatnextnumofchieku
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
                                st.error(f"Error: manyitemformatofnextnum（{current_poly_degree}）isyuni-kunataimupointonum（{unique_points}）orupwithsu。")
                                st.error(f"anataofDatawithuseusepossiblenamaximumnextnum: {unique_points - 1}")
                                st.error(f"taimupointo: {sorted(numeric_values.unique())}")
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
                                st.warning(f"Warning: timebetweenchangenum '{time_var}' isnumvalwithisarimasen。manyitemformatmoderu（nextnum{current_poly_degree}）theuseusedotoisnumvalisrequiredwithsu。moderuiscollbundleshinotpossiblenatureisarimasu。")
                else:
                    st.warning(f"Warning: timebetweenchangenum '{time_var}' isrealtestdezainFiletoviewtsukarimasen。")
            except Exception as e:
                st.warning(f"realtestdezainFileofLoadingmidtoErrorisoccurgenshimashita: {str(e)}")

    st.markdown("""
--------------------------------------------------------------------------
        """)
    if st.button('Run analysis'):
        #mazuRofdftochangechange
        if test_method == 'DESeq2-LRT':
    #        ro.r.assign('cts',cts) # ErrorisoutruofwithFiletooneonceSavedo
            r.assign('df',df)
            r.assign('df_e',df_e)

            # changenumtaipuinfoinfotheRsidetotransfersu
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
            #mazubekuta-tochangechange
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
                # changenumtheRtotransfersu
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
            # Save input to files for R import
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            voom_plot_path = os.path.join(res_dir, 'voom_plot.png')
        #    if os.path.exists(voom_plot_path):
        #        os.remove(voom_plot_path)

            ro.r.assign('temp_dir', temp_dir)
            
            # changenumtaipuinfoinfotheRsidetotransfersu
            continuous_vars = [col for col in condition_col if var_types[col] == "continuous"]
            r_continuous_vars = ro.StrVector(continuous_vars)
            ro.r.assign('continuous_vars', r_continuous_vars)
            
            # polynomialrelconnectofchangenumtheRtotransfersu
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
                
                # Identify coefficients specific to the full model
                add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))
                
                # For logit-transformed data, skip voom and directly fit with limma
                fit_full <- lmFit(counts, design_full)
                fit_full <- eBayes(fit_full)
                
                if (length(add_coefs) == 1) {{
                  res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                }} else {{
                  cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                  colnames(cm) <- add_coefs
                  for (i in 1:length(add_coefs)) {{
                    cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                  }}
                  
                  fit_contrast <- contrasts.fit(fit_full, cm)
                  fit_contrast <- eBayes(fit_contrast)
                  
                  res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                }}
                write.table(res, file='{res_dir}/limma_res.tsv', sep='\t', quote=FALSE, col.names=NA)
                sink()
                """
            elif limma_count:
                # For regular count data, use the original approach
                ro.r.assign('independentFiltering', independentFiltering)

                log_file = f"{res_dir}/limma_debug.log"
                ro.r.assign('log_file', log_file)

                r_code = f"""

                # roguFiletoridairekuto
                sink('{log_file}', append=FALSE, split=TRUE)
                library(edgeR)
                library(limma)
                tryCatch({{
                    counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')

                    # weightneed: kauntoDataofarrangenumchangechangetheConfirm
                    counts <- round(counts)  # realnumthekauntoDataandshiteuseuplacematchisarrangenumtochangechange

                    coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')
                    y <- DGEList(counts=counts)


                    # Filtering（Option）
                    if (independentFiltering) {{
                        keep <- filterByExpr(y, design=model.matrix(as.formula('{full_model}'), data=coldata))
                        y <- y[keep, ]
                    }}

                    # Normalization
                    y <- calcNormFactors(y)

                    design_full <- model.matrix(as.formula('{full_model}'), data=coldata)
                    design_reduced <- model.matrix(as.formula('{reduced_model}'), data=coldata)
                    
                    add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))
                    
                    png('{res_dir}/voom_plot.png')
                    v <- voom(y, design_full, plot=TRUE)
                    dev.off()

                    # voomResultConfirm
                    print("Post-voom data:")
                    print(dim(v$E))
                    print(range(v$E))
                    
                    fit_full <- lmFit(v, design_full)
                    fit_full <- eBayes(fit_full)
                    
                    if (length(add_coefs) == 1) {{
                      res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                    }} else {{
                      cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                      colnames(cm) <- add_coefs
                      for (i in 1:length(add_coefs)) {{
                        cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                      }}
                      
                      fit_contrast <- contrasts.fit(fit_full, cm)
                      fit_contrast <- eBayes(fit_contrast)
                      
                      res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                    }}
                    write.table(res, file='{res_dir}/limma_res.tsv', sep='\t', quote=FALSE,  col.names=NA)

                    }}, error = function(e) {{
                        cat("\\n===== ERROR =====\\n")
                        cat("Error in limma analysis:", conditionMessage(e), "\\n")
                        cat("Error occurred at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
                        traceback()
                    }}, finally = {{
                        cat("\\n===== ANALYSIS COMPLETE =====\\n")
                        cat("Time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n")
                        # sinkthesolveremove
                        sink()
                    }})

                """

                st.write(f"{res_dir}/voom_plot.png")
                st.image(f"{res_dir}/voom_plot.png", caption='Voom mean-variance trend')

            else:
                r_code = f"""
                sink()
                sink(paste0(temp_dir, "/limma_output.txt"))
                library(limma)
                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\t')
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\t')

                design_full <- model.matrix(as.formula('{full_model}'), data=coldata)
                design_reduced <- model.matrix(as.formula('{reduced_model}'), data=coldata)

                add_coefs <- setdiff(colnames(design_full), colnames(design_reduced))

                # nonkauntoDatanaofwithvoomissukipu
                fit_full <- lmFit(counts, design_full)
                fit_full <- eBayes(fit_full)

                if (length(add_coefs) == 1) {{
                  res <- topTable(fit_full, coef=add_coefs, number=Inf, adjust.method='fdr')
                }} else {{
                  cm <- matrix(0, nrow=ncol(design_full), ncol=length(add_coefs))
                  colnames(cm) <- add_coefs
                  for (i in 1:length(add_coefs)) {{
                    cm[which(colnames(design_full) == add_coefs[i]), i] <- 1
                  }}
                  
                  fit_contrast <- contrasts.fit(fit_full, cm)
                  fit_contrast <- eBayes(fit_contrast)
                  
                  res <- topTable(fit_contrast, number=Inf, sort.by='F', adjust.method='fdr')
                }}
                write.table(res, file='{res_dir}/limma_res.tsv', sep='\t', quote=FALSE, col.names=NA)
                sink()
                """

            ro.r(r_code)
            res_df = pd.read_csv(os.path.join(res_dir, 'limma_res.tsv'), sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)

            file_name = file_name_head + "_limma_" + full_model.replace(" + ", "_").replace(" ", "") + "_vs_" + reduced_model.replace(" + ", "_").replace(" ", "")

            shutil.make_archive("res", format='zip',root_dir= res_dir)


        elif test_method == 'Beta Regression':

            # FileofSaveandSettingsissame
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

            # 0-1ofboundworldofadjustarrange
            eps <- {epsilon}
            counts <- pmax(pmin(counts, 1-eps), eps)

            # aligncolprocprocClusterofSettings
            n_cores <- {n_cores}
            cl <- makeCluster(n_cores)

            # aligncolprocproctorequirednapake-jitheClustertoLoading
            clusterEvalQ(cl, {{
              library(betareg)
              library(lmtest)
            }})

            time_var <- "{full[0]}"
            cat("time_var")
            cat(time_var)

            # timebetweenchangenumofConfirmandchangechange
            if(time_var %in% colnames(coldata)) {{
              cat("Time variable exists in coldata. Values:", "\\n")
              print(coldata[[time_var]])
              
              # timebetweenchangenumisnumvalwithnotplacematchischangechange
              if(!is.numeric(coldata[[time_var]])) {{
                cat("Converting time variable to numeric\\n")
                # numvalextractoutandchangechange
                coldata[[time_var]] <- as.numeric(gsub("[^0-9.]", "", as.character(coldata[[time_var]])))
                cat("After conversion:", "\\n")
                print(coldata[[time_var]])
              }}
            }} else {{
              cat("WARNING: Time variable not found in coldata!\\n")
            }}

            coldata[[time_var]] <- coldata[[time_var]] / max(coldata[[time_var]]) #maximumvaltoyorunormalization

            # manyitemformatofnextnumtobasedukuitemthestructbuild
            polynomial_terms <- ""
            if ({beta_polynomial_degree} >= 2) {{
              polynomial_terms <- paste0(polynomial_terms, " + I(", time_var, "^2)")
            }}
            if ({beta_polynomial_degree} >= 3) {{
              polynomial_terms <- paste0(polynomial_terms, " + I(", time_var, "^3)")
            }}

            # changenumtheClustertosendtrust
            clusterExport(cl, c("counts", "coldata", "eps", "time_var", "polynomial_terms"))

            # procprocstartmese-ji
            cat("Starting parallel beta regression on", n_cores, "cores for", nrow(counts), "genes\\n")

            # tesutomoderuRun
            test_model_result <- tryCatch({{
              # tesutouseofData
              test_gene_data <- data.frame(y=as.numeric(counts[1,]), coldata)
              
              # fuo-miyurathestructbuild
              full_formula <- paste("{full_model.replace('~', '')}", polynomial_terms)
              
              # moderufitmatchthetryrow
              test_fit <- betareg(as.formula(paste("y ~", full_formula)), data=test_gene_data)
              "success"
            }}, error=function(e) {{
              # Errormese-jithereturnsu
              return(conditionMessage(e))
            }})

            # Errorkaunta-initialize
            error_counter <- 0
            error_message <- ""

            # aligncolprocprocrelnum
            process_gene <- function(i) {{
              gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
              
              full_formula <- paste("{full_model.replace('~', '')}", polynomial_terms)
              reduced_formula <- "{reduced_model.replace('~', '')}"
              
              tryCatch({{
                # furumoderuandrideyu-sudomoderuoffuiteingu
                full_fit <- betareg(as.formula(paste("y ~", full_formula)), data=gene_data)
                reduced_fit <- betareg(as.formula(paste("y ~", reduced_formula)), data=gene_data)
                
                # likelydegreeratiocheckset
                lr_test <- lrtest(reduced_fit, full_fit)
                
                # Resultthereturnsu
                c(statistic = lr_test$Chisq[2],
                  df = lr_test$Df[2],
                  p_value = lr_test$`Pr(>Chisq)`[2],
                  logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
              }}, error=function(e) {{
                # ErrorisoccurgenshitaplacematchisNAthereturnsu
                if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
              }})
            }}

            # aligncolprocprocofRun
            system.time(
              results_list <- parLapply(cl, 1:nrow(counts), process_gene)
            )

            # Clusterofend
            stopCluster(cl)

            # Resultthematorikusutochangechange
            results_matrix <- do.call(rbind, results_list)
            rownames(results_matrix) <- rownames(counts)

            # NULLofnumthekaunto
            na_count <- sum(is.na(results_matrix[, "statistic"]))
            total_genes <- nrow(results_matrix)
            na_percent <- round(100 * na_count / total_genes, 2)

            # ResultFiletomoderuinfoinfotheaddadd
            cat("\\n### moderucollbundleinfoinfo ###\\n", file='{res_dir}/model_convergence_info.txt')
            cat("manyitemformatnextnum:", {beta_polynomial_degree}, "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            cat("furumoderuformat:", paste("y ~", paste("{full_model.replace('~', '')}", polynomial_terms)), "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            cat("shrinksmallmoderuformat:", paste("y ~", "{reduced_model.replace('~', '')}"), "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)

            if (na_count > 0) {{
              cat("Warning: ", na_count, " pieceofGene (", na_percent, "%) withmoderuiscollbundleshimasenwithshita。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              
              if (na_count == total_genes) {{
                cat("allofGenewithmoderuiscollbundleshimasenwithshita。manyitemformatnextnumthebelowgeru、timebetweenchangenumthesuke-ringu(currentmaximumvalwithsuke-ringushiteexist)、notavgeqoftimebetweentheavgeqizedokoandthecheckdiscussshitekudasai。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
                cat("tesutomoderuofError: ", test_model_result, "\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              }} else {{
                cat(total_genes - na_count, " pieceofGene (", 100 - na_percent, "%) withcorrectnormaltoAnalysiswithkimashita。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
              }}
            }} else {{
              cat("allofGenewithcorrectnormaltomoderuiscollbundleshimashita。\\n", file='{res_dir}/model_convergence_info.txt', append=TRUE)
            }}

            # manyweightchecksetsuppcorrect
            results_matrix <- cbind(results_matrix, 
                                  adj.P.Val = p.adjust(results_matrix[, "p_value"], method="BH"))

            # ResultofSave
            res_df <- as.data.frame(results_matrix)
            res_df <- res_df[order(res_df$p_value), ]
            write.table(res_df, file='{res_dir}/betareg_res.tsv', sep='\\t', quote=FALSE, col.names=NA)
            sink()
            """

        elif test_method == 'Generalized Linear Model (GLM)':
            # FileofSaveandSettings
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)

            ro.r.assign('temp_dir', temp_dir)
            
            # changenumtaipuinfoinfotheRsidetotransfersu
            continuous_vars = [col for col in condition_col if var_types[col] == "continuous"]
            r_continuous_vars = ro.StrVector(continuous_vars)
            ro.r.assign('continuous_vars', r_continuous_vars)
            
            # polynomialrelconnectofchangenumtheRtotransfersu
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

            # connectcontinuechangenumtheprocproc
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

            # remainriofchangenumthecausechildtypetochangechange（connectcontinuechangenumorout）
            for (i in c(1:dim(coldata)[2])) {{
                col_name <- colnames(coldata)[i]
                if (!exists("continuous_vars") || !(col_name %in% continuous_vars)) {{
                    cat(paste0("Treating '", col_name, "' as categorical variable\\n"))
                    coldata[,i] <- factor(coldata[,i])
                }}
            }}

            # divdistfuamiri-torespondjitaDataPreprocessing
            if ("{glm_dist_short}" == "beta") {{
                # 0-1ofboundworldofadjustarrange
                eps <- {glm_epsilon}
                counts <- pmax(pmin(counts, 1-eps), eps)
            }} else if ("{glm_dist_short}" == "gaussian") {{
                # gausuofplacematchisspectoPreprocessingnotneed
            }} else if ("{glm_dist_short}" == "poisson" || "{glm_dist_short}" == "nb") {{
                # poasonandNBiskauntoDatathe想set（arrangenumize）
                counts <- round(counts)
            }}

            # divdistfuamiri-ofSettingsrelnum
            get_family <- function() {{
                if ("{glm_dist_short}" == "beta") {{
                    # Betaofplacematchisbetaregpake-jitheuseuse
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

            # aligncolClusterofSettings
            n_cores <- {n_cores}
            cl <- makeCluster(n_cores)

            # Clustertorequirednapake-jiandDatathesendtrust
            clusterEvalQ(cl, {{
                library(MASS)
                if ("{glm_dist_short}" == "beta") {{
                    library(betareg)
                    library(lmtest)
                }}
            }})

            clusterExport(cl, c("counts", "coldata", "get_family"))

            # procprocrelnum
            process_gene <- function(i) {{
                gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
                family_to_use <- get_family()
                
                result <- tryCatch({{
                    if ("{glm_dist_short}" == "beta") {{
                        # Betatimereturnofplacematchisbetaregpake-jitheuseuse
                        library(betareg)
                        full_fit <- betareg(as.formula("{full_model}"), 
                                          data=gene_data, 
                                          link=family_to_use$link)
                        reduced_fit <- betareg(as.formula("{reduced_model}"), 
                                             data=gene_data, 
                                             link=family_to_use$link)
                        
                        # likelydegreeratiocheckset
                        library(lmtest)
                        lr_test <- lrtest(reduced_fit, full_fit)
                        
                        # Resultthereturnsu
                        c(statistic = lr_test$Chisq[2],
                          df = lr_test$Df[2],
                          p_value = lr_test$`Pr(>Chisq)`[2],
                          logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
                    }} else {{
                        # passnormalofGLMofplacematch
                        full_fit <- glm(as.formula("{full_model}"), family=family_to_use, data=gene_data)
                        reduced_fit <- glm(as.formula("{reduced_model}"), family=family_to_use, data=gene_data)
                        
                        # likelydegreeratiocheckset
                        anova_result <- anova(reduced_fit, full_fit, test="LRT")
                        
                        # Resultthereturnsu
                        c(statistic = anova_result$Deviance[2],
                          df = anova_result$Df[2],
                          p_value = anova_result$`Pr(>Chi)`[2],
                          logLik_diff = logLik(full_fit) - logLik(reduced_fit))
                    }}
                }}, error=function(e) {{
                    # ErrorofplacematchisNAthereturnsu
                    if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                    c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
                }})
                
                return(result)
            }}

            # aligncolprocprocofRun
            cat("Starting parallel GLM regression on", n_cores, "cores for", nrow(counts), "genes\\n")
            cat("Using distribution family:", "{glm_dist_family}", "\\n")
            system.time(
                results_list <- parLapply(cl, 1:nrow(counts), process_gene)
            )

            # Clusterofend
            stopCluster(cl)

            # ResulttheDatafure-mutochangechange
            results_matrix <- do.call(rbind, results_list)
            results <- as.data.frame(results_matrix)
            rownames(results) <- rownames(counts)

            # NAofnumthekaunto
            na_count <- sum(is.na(results$statistic))
            if (na_count > 0) {{
                cat("Warning:", na_count, "pieceofGenewithmoderuiscollbundleshimasenwithshita (", 
                    round(100 * na_count / nrow(results), 2), "%)\\n")
            }}

            # manyweightchecksetsuppcorrect
            results$adj.P.Val <- p.adjust(results$p_value, method="BH")

            # moderuinfoinfoofSave
            cat("\\n### GLMmoderuinfoinfo ###\\n", file='{res_dir}/glm_model_info.txt')
            cat("divdistfuamiri-:", "{glm_dist_family}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("rinkurelnum:", "{glm_link}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("furuGLMmoderuformat:", "{full_model}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)
            cat("shrinksmallGLMmoderuformat:", "{reduced_model}", "\\n", file='{res_dir}/glm_model_info.txt', append=TRUE)

            # ResulttheSave
            write.table(results[order(results$p_value), ], 
                        file='{res_dir}/glm_{glm_dist_short}_{glm_link}_res.tsv', 
                        sep='\\t', quote=FALSE, col.names=NA)

            cat("GLM regression analysis completed\\n")
            sink()
            """

        elif test_method == 'Generalized Additive Model (GAM)':
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


            # furthernewsaretaGAM Rko-do
            r_code = f"""
                sink()
                sink(paste0({temp_dir}, "/GAM_output.txt"))
                library(mgcv)
                library(lmtest)
                library(parallel)

                counts <- read.table('{counts_file}', header=TRUE, row.names=1, sep='\\t')
                coldata <- read.table('{coldata_file}', header=TRUE, sep='\\t')

                # divdistfuamiri-torespondjitaDataPreprocessing
                if ("{dist_short}" == "beta") {{
                    # 0-1ofboundworldofadjustarrange
                    eps <- {epsilon}
                    counts <- pmax(pmin(counts, 1-eps), eps)
                }} else if ("{dist_short}" == "gaussian") {{
                    # gausuofplacematchisspectoPreprocessingnotneed
                }} else if ("{dist_short}" == "poisson" || "{dist_short}" == "nb") {{
                    # poasonandNBiskauntoDatathe想set（arrangenumize）
                    counts <- round(counts)
                }}

                # timebetweenchangenumofspecsetandfuo-miyuraofmakebecome
                time_var <- "{full[0]}"
                cat("Time variable:", time_var, "\\n")

                # timebetweenchangenumofConfirmandchangechange
                if(time_var %in% colnames(coldata)) {{
                  cat("Time variable exists in coldata. Values:", "\\n")
                  print(coldata[[time_var]])
                  
                  # timebetweenchangenumisnumvalwithnotplacematchischangechange
                  if(!is.numeric(coldata[[time_var]])) {{
                    cat("Converting time variable to numeric\\n")
                    coldata[[time_var]] <- as.numeric(gsub("[^0-9.]", "", as.character(coldata[[time_var]])))
                    cat("After conversion:", "\\n")
                    print(coldata[[time_var]])
                  }}
                }} else {{
                  cat("WARNING: Time variable not found in coldata!\\n")
                }}

                # normalization（be-tadivdistoroutwithmouseeru）
                if ({beta_normalization} == "TRUE"){{
                    coldata[[time_var]] <- coldata[[time_var]] / max(coldata[[time_var]])
                    cat("Time is normalized by max value.")
                }}

                # GAMmoderuformatofmakebecome
                gam_full_formula <- "{full_model.replace('~', '')}"
                gam_reduced_formula <- "{reduced_model.replace('~', '')}"

                # flatsmoothizeitemofaddadd (timebetweenchangenumtopairshite)
                if(time_var %in% colnames(coldata)) {{
                  if(length(unique(coldata[[time_var]])) >= 3) {{  # fewnakuandmo3tsuofdiffbecomevalisrequired
                    # furumoderutoofmiflatsmoothizeitemtheaddadd
                    if(grepl(time_var, gam_full_formula)) {{
                      gam_full_formula <- gsub(
                        paste0("\\\\b", time_var, "\\\\b"), 
                        paste0("s(", time_var, ", k={gam_k}, bs='{spline_type}')"), 
                        gam_full_formula
                      )
                    }} else {{
                      # time_varisclearshowaltoincludemareteinotplacematchisaddadd
                      gam_full_formula <- paste(gam_full_formula, "+", paste0("s(", time_var, ", k={gam_k}, bs='{spline_type}')"))
                    }}
                    cat("Full GAM formula:", gam_full_formula, "\\n")
                    cat("Reduced GAM formula:", gam_reduced_formula, "\\n")
                  }} else {{
                    cat("Not enough unique time points for smoothing, using linear terms\\n")
                  }}
                }}

                cat("{dist_short}")

                # divdistfuamiri-ofSettingsrelnum
                get_family <- function() {{
                    if ("{dist_short}" == "beta") {{
                        return(betar())
                    }} else if ("{dist_short}" == "gaussian") {{
                        return(gaussian())
                    }} else if ("{dist_short}" == "poisson") {{
                        return(poisson())
                    }} else if ("{dist_short}" == "nb") {{
                        library(mgcv)
                        # mgcvofnbrelnumisthetaParametertherecvkegetru
                        return(negbin(theta = {nb_theta}))
                    }}
                }}

                # aligncolClusterofSettings
                n_cores <- {n_cores}
                cl <- makeCluster(n_cores)

                # Clustertorequirednapake-jiandDatathesendtrust
                clusterEvalQ(cl, {{
                  library(mgcv)
                  library(lmtest)
                }})

                clusterExport(cl, c("counts", "coldata", "gam_full_formula", 
                                    "gam_reduced_formula", "get_family"))

                # procprocrelnum
                process_gene <- function(i) {{
                  gene_data <- data.frame(y=as.numeric(counts[i,]), coldata)
                  family_to_use <- get_family()
                  
                  result <- tryCatch({{
                    # furumoderuandshrinksmallmoderuoffuiteingu
                    full_fit <- gam(as.formula(paste("y ~", gam_full_formula)), 
                                    family=family_to_use, data=gene_data, method="{gam_method}")
                    
                    reduced_fit <- gam(as.formula(paste("y ~", gam_reduced_formula)), 
                                       family=family_to_use, data=gene_data, method="{gam_method}")
                    
                    # likelydegreeratiocheckset
                    lr_test <- lrtest(reduced_fit, full_fit)
                    
                    # Resultthereturnsu
                    c(statistic = lr_test$Chisq[2],
                      df = lr_test$Df[2],
                      p_value = lr_test$`Pr(>Chisq)`[2],
                      logLik_diff = lr_test$LogLik[2] - lr_test$LogLik[1])
                  }}, error=function(e) {{
                    # ErrorofplacematchisNAthereturnsu
                    if (i <= 5) cat("Gene", i, "Error:", conditionMessage(e), "\\n")
                    c(statistic = NA, df = NA, p_value = NA, logLik_diff = NA)
                  }})
                  
                  return(result)
                }}

                # aligncolprocprocofRun
                cat("Starting parallel GAM regression on", n_cores, "cores for", nrow(counts), "genes\\n")
                cat("Using distribution family:", "{dist_family}", "\\n")
                system.time(
                  results_list <- parLapply(cl, 1:nrow(counts), process_gene)
                )

                # Clusterofend
                stopCluster(cl)

                # ResulttheDatafure-mutochangechange
                results_matrix <- do.call(rbind, results_list)
                results <- as.data.frame(results_matrix)
                rownames(results) <- rownames(counts)

                # moderuinfoinfoofSave
                cat("\\n### moderuinfoinfo ###\\n", file='{res_dir}/gam_model_info.txt')
                cat("divdistfuamiri-:", "{dist_family}", "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("GAMflatsmoothizeParameter k:", {gam_k}, "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("Estimationwaymethod:", "{gam_method}", "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("furuGAMmoderuformat:", paste("y ~", gam_full_formula), "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)
                cat("shrinksmallGAMmoderuformat:", paste("y ~", gam_reduced_formula), "\\n", file='{res_dir}/gam_model_info.txt', append=TRUE)

                # NAofnumthekaunto
                na_count <- sum(is.na(results$statistic))
                if (na_count > 0) {{
                  cat("Warning:", na_count, "pieceofGenewithmoderuiscollbundleshimasenwithshita (", 
                      round(100 * na_count / nrow(results), 2), "%)  k, spline, notavgeqtimesyscolofNormalizationeqthecheckdiscuss", 
                      file='{res_dir}/gam_model_info.txt', append=TRUE)
                }}

                # manyweightchecksetsuppcorrect
                results$adj.P.Val <- p.adjust(results$p_value, method="BH")

                # ResulttheSave
                write.table(results[order(results$p_value), ], 
                            file='{res_dir}/gam_{dist_short}_{spline_type}_res.tsv', 
                            sep='\\t', quote=FALSE, col.names=NA)

                cat("GAM regression analysis completed\\n")
                sink()
            """



        elif test_method == 'maSigPro':
            # FileSaveSettings
            counts_file = os.path.join(temp_dir, 'counts.tsv')
            df.to_csv(counts_file, sep='\t')
            coldata_file = os.path.join(temp_dir, 'coldata.tsv')
            df_e.to_csv(coldata_file, sep='\t', index=False)


            # timebetweeninfoinfotheincludemufitcutnaedesignDatafure-muthemakebecomedoRko-do
            r_code = f"""
            library(maSigPro)
            
            # DataofLoading
            cat("Loading data...\\n")
            counts <- read.table("{counts_file}", header=TRUE, row.names=1, sep="\\t")
            coldata <- read.table("{coldata_file}", header=TRUE, sep="\\t")
            print(coldata)

            # maSigProuseoffitcutnadezainDatafure-muthemakebecome
            cat("Creating proper design matrix for maSigPro...\\n")
            
            # Groupinfoinfothegetget（"Group"col）
            time_col <- as.character(coldata${full[0]})
            

            # timebetweeninfoinfotheextractout（Example："0w", "1w", "4w"fromnumvaltochangechange）
            time_values <- as.numeric(gsub("[^0-9.]", "", time_col))
            cat("time_values")
            cat(time_values)
            
            # repurike-toinfoinfothemakebecome
            # sametimebetweenvaltheholdtsuSampletoyuni-kunanumnumtheratioricurrentte
            replicates <- numeric(length(time_values))
            for (t in unique(time_values)) {{
                idx <- which(time_values == t)
                replicates[idx] <- 1:length(idx)
            }}
            
            # maSigProuseofcorrectshiiedesignDatafure-muthemakebecome
            edesign <- data.frame(
                Time = time_values,
                Replicate = replicates
            )
            rownames(edesign) <- colnames(counts)
            
      #      # otherofrealtestConditionisarebaaddadd
      #      if (ncol(coldata) > 1) {{
      #          for (i in 2:ncol(coldata)) {{
      #              col_name <- colnames(coldata)[i]
      #              edesign[[col_name]] <- coldata[[i]]
      #          }}
      #      }}

            # Add Group column (all 1s for single condition)
            edesign$Group <- rep(1, nrow(edesign)) #tutorialtomatchwaseteallpart1todo
            
            # DatatypeofConfirm
            cat("Time values:", paste(time_values, collapse=", "), "\\n")
            cat("Time values are numeric:", is.numeric(edesign$Time), "\\n")
            cat("edesign")
            print(edesign)
            
            # Preprocessing
            """
            
            # Datataiputorespondjitaprocprocofaddadd
            if data_type == "0-1 data (logit transformation)":
                r_code += f"""
            # 0-1Dataofprocproc
            eps <- {epsilon}
            counts <- pmax(pmin(counts, 1-eps), eps)
            counts <- log(counts/(1-counts))
            use_counts_param <- FALSE
            cat("Applied logit transformation\\n")
            """
            elif data_type == "qPCR/continuous data (Gaussian)":
                r_code += f"""
            # qPCR/connectcontinueDataofprocproc
            use_counts_param <- FALSE
            cat("Using Gaussian model for continuous data\\n")
            
            # DataofPreprocessing
            original_counts <- counts
            """
                
                # logchangechangeOption
                if log_transform:
                    r_code += """
            # Log2changechangethefituse
            counts <- log2(counts + 1)  # +1theaddete0valtheavoidkeru
            cat("Applied log2(x+1) transformation\\n")
            """
                
                # NormalizationOption
                if normalization:
                    r_code += """
            # Z-scoreNormalization（Samplebetween）
            counts <- t(scale(t(counts)))
            cat("Applied z-score normalization across samples\\n")
            """
                
                r_code += """
            cat("qPCR data preprocessing completed\\n")
            """
            else:  # RNA-seq count data
                r_code += """
            # RNA-seqDataofprocproc
            use_counts_param <- TRUE
            cat("Using GLM for count data\\n")
            """
            
            # dezainrowcoltheuseuseshitadivanalyze
            r_code += f"""
            # dezainrowcoltheuseuseshitedivanalyze
            cat("Running maSigPro analysis...\\n")
            
            # pointsetsaretanextnumwithdezainrowcolthemakebecome
            design <- make.design.matrix(edesign, degree={degree})
            
            # timereturndivanalyzeofRun
            cat("Running p.vector...\\n")
            fit <- p.vector(counts, design$edesign, Q={q_value}, MT.adjust="none", counts=use_counts_param)
           # fit <- p.vector(counts, design$edesign, Q={q_value}, MT.adjust="BH", counts=use_counts_param)
            
            # SignificantnaGenenumofConfirm
            sig_count <- sum(fit$p < {q_value}, na.rm=TRUE)
            cat("Genes with p <", {q_value}, ":", sig_count, "\\n")

            # After running p.vector() and finding no significant genes

            
            # SignificantnaGeneisexistplacematchofminextofsutepuhe
            if (sig_count > 0) {{
                cat("Running T.fit...\\n")
                tstep <- T.fit(fit, step.method="backward", alfa={q_value})
                
                cat("Getting significant genes...\\n")
                sigs <- get.siggenes(tstep, rsq={rsq}, vars="each")
                
                # ResulttheSave
                if (!is.null(sigs) && !is.null(sigs$sig.genes) && !is.null(sigs$sig.genes$sig.profiles) && nrow(sigs$sig.genes$sig.profiles) > 0) {{
                    cat("Found", nrow(sigs$sig.genes$sig.profiles), "significant genes\\n")
                    
                    # puroFileandrelatenumtheSave
                    write.table(sigs$sig.genes$sig.profiles, file="{res_dir}/maSigPro_sig_profiles.tsv", sep="\\t", quote=FALSE)
                    write.table(sigs$coefficients, file="{res_dir}/maSigPro_coefficients.tsv", sep="\\t", quote=FALSE)
                    
                    # needaboutinfoinfoofSave
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
            
            # Rko-dotheRun
            with st.spinner('Calculating maSigPro... This may take a while.'):
                try:
                    # Rko-dothedebaguusetoSave
                    with open(os.path.join(temp_dir, 'debug_maSigPro.R'), 'w') as f:
                        f.write(r_code)
                    
                    # Rko-doRun
                    ro.r(r_code)
                    
                    # ResultuseofemptyDatafure-muthemakebecome（Errortimeavoiduse）
                    res_df = pd.DataFrame()
                    
                    # ResultofDisplay
                    st.markdown("### maSigPro Analysis Results")
                    
                    # ResultFileofConfirmandDisplay
                    summary_file = os.path.join(res_dir, 'maSigPro_summary.txt')
                    if os.path.exists(summary_file):
                        with open(summary_file, 'r') as f:
                            summary = f.read()
                        st.text(summary)
                    
                    # SignificantnaGeneofResulttheDisplay
                    sig_profiles_file = os.path.join(res_dir, 'maSigPro_sig_profiles.tsv')
                    if os.path.exists(sig_profiles_file):
                        res_df = pd.read_csv(sig_profiles_file, sep='\t', index_col=0)
                        if not res_df.empty:
                            st.write("### Top significant genes:")
                            st.dataframe(res_df.head(10))
                            
                            # relatenumofDisplay
                            coef_file = os.path.join(res_dir, 'maSigPro_coefficients.tsv')
                            if os.path.exists(coef_file):
                                st.write("### Regression coefficients:")
                                coef = pd.read_csv(coef_file, sep='\t', index_col=0)
                                st.dataframe(coef.head(10))
                            
                            # ResultofDownload
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
                    # ErrortimetomoemptyofDataFramethemakebecome
                    res_df = pd.DataFrame()
            
     
            
            # ResultofZIPgenbecome
            data_type_short = ""
            if data_type == "RNA-seq count data (GLM)":
                data_type_short = "RNAseq"
            elif data_type == "qPCR/continuous data (Gaussian)":
                data_type_short = "qPCR"
            elif data_type == "0-1 data (logit transformation)":
                data_type_short = "logit"
            
            file_name = file_name_head + f"_maSigPro_{data_type_short}"
            shutil.make_archive("res", format='zip', root_dir=res_dir)

        # ResultofDisplayandSave
        if test_method == 'Beta Regression':
            ro.r(r_code)
            res_df = pd.read_csv(os.path.join(res_dir, 'betareg_res.tsv'), sep='\t', index_col=0)
            st.write(f"Significant (FDR<0.05): {(res_df['adj.P.Val']<0.05).sum()}")
            st.dataframe(res_df)

            # moderucollbundleinfoinfoofConfirmandDisplay
            convergence_file = os.path.join(res_dir, 'model_convergence_info.txt')
            if os.path.exists(convergence_file):
                with open(convergence_file, 'r') as f:
                    convergence_info = f.read()
        
            # collbundletoquesttopicisexistkachieku
            if "Warning" in convergence_info:
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
            
            # moderuinfoinfoofDisplay
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


            # moderucollbundleinfoinfoofConfirmandDisplay
            convergence_file = os.path.join(res_dir, 'gam_model_info.txt')
            if os.path.exists(convergence_file):
                with open(convergence_file, 'r') as f:
                    convergence_info = f.read()
        
            # collbundletoquesttopicisexistkachieku
            if "Warning" in convergence_info:
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


#　DatathesendrubeforetoallzeroofDataisremovekubeki


# refispointsetsareteexistandkiisFilenametheadjustarrangedo?
