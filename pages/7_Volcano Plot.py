
import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import shutil
import plotly
import os
import re
import Levenshtein
import sys
from helper_func import mk_temp_dir
from pdf2jpg import pdf2jpg
import math

st.set_page_config(page_title="Volcano_Plot", page_icon="🌋")
st.markdown("### Volcano plot")

@st.cache_data
def convert_df(df):
   return df.to_csv(index=False, sep='\t', header = None).encode('utf-8')

@st.cache_data
def read_xl(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_xl = pd.read_csv(file, index_col = index_col, header = 0, sep = sep)
    return df_xl

def hex2rgb(hex_color, c_alpha):
    rgb = 'rgba' + plotly.colors.convert_colors_to_same_type(hex_color, colortype='rgb')[0][0][3:-1] + ", " + str(c_alpha)  + ")"
    return rgb

def pyhex2rgb(hex_color, c_alpha):
    t = plotly.colors.convert_colors_to_same_type(hex_color, colortype='rgb')[0][0][4:-1].split(', ')
    a = [float(x)/255 for x in t]
    a.append(c_alpha)
    return np.array([a])


# Function to read gene lists from files
@st.cache_data
def read_gene_list(file_path):
    try:
        with open(file_path, 'r') as f:
            # Read lines and strip whitespace
            genes = [line.strip() for line in f.readlines() if line.strip()]
        return genes
    except FileNotFoundError:
        st.warning(f"File not found: {file_path}")
        return []


st.sidebar.title("Options")
st.markdown("##### Options are displayed at the bottom of the left side panel")
with st.sidebar:
    goi = None
    use_goi = st.checkbox("Specify gene functional class to label?", value=False)
    if use_goi:
        species = st.radio("Species", ['mouse','human'], index=0)
        selected_class = st.selectbox("Functional class",["TF","Epigenetic regulator","CytokineGF", "LR-pair Ligand", "LR-pair Receptor", "Secreted factor",], index = 0)
        if species == 'mouse':
            # Define the functional classes and their corresponding file paths
            functional_classes = {
                "TF": "Mus_musculus_TF_llst.txt",
                "Epigenetic regulator": "Mus_musculus_epigenetic_factors_list.txt",
                "CytokineGF": "cytokine-genes-mouse_symbols.txt", 
                "LR-pair Ligand": "Mouse_ligand.txt",
                "LR-pair Receptor": "Mouse_receptor.txt",
                "Secreted Factor": "secreted.mouse.list.txt"
            }
        else:
            functional_classes = {
                "TF": "Homo_sapiens_TF_list.txt",
                "Epigenetic regulator": "Human_epigenetic_factors_list.txt",
                "CytokineGF": "cytokine-genes-mouse_symbols.txt", 
                "LR-pair Ligand": "Human_ligand.txt",
                "LR-pair Receptor": "Human_receptor.txt",
                "Secreted Factor": "secreted.mouse.list.txt"
            }
            st.write("Human secreted factor list is not available yet")


        # Read the gene list for the selected class
        file_path = "db/functionalclass/" + functional_classes[selected_class]
        goi = read_gene_list(file_path)
# temp内に保存する
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir, res_dir = mk_temp_dir("Volcano")
    st.session_state.temp_dir = temp_dir
else:
    temp_dir = st.session_state.temp_dir
    temp_dir, res_dir = mk_temp_dir("Volcano", temp_dir)


use_upload = 'Yes'
if 'deseq2' in st.session_state:
    if st.session_state.deseq2 is not None:
        use_upload = st.radio("Upload new file?", ('Yes','No'), index = 1)
    if use_upload == "No":
        df = st.session_state.deseq2
        df['Gene'] = df.index
        if  "deseq2_uploaded_file_name" in st.session_state:
            file_name_head = st.session_state.deseq2_uploaded_file_name
        else:
            file_name_head = "res"
        input_file_type = 'tsv'
        if "Row_name" in df.columns.to_list(): # Row_nameを含むとき
            df = df.set_index('Row_name')
            df.index.name = "Gene"

if "inv" not in st.session_state:
    st.session_state.inv = False #inverseのスイッチ invされている間はTrue

input_file_type = st.radio(
    "Data format:",
    ('tsv','csv', 'excel'))


if use_upload == 'Yes':
    uploaded_file = st.file_uploader(" ", type=['txt','tsv','csv','xls','xlsx'])
    if uploaded_file is not None:
        if input_file_type == "csv":
            df = read_csv(uploaded_file)
        elif input_file_type == 'tsv':
            df = read_csv(uploaded_file, sep = '\t')
        else:
            df = read_xl(uploaded_file)
        st.session_state.deseq2 = df
        file_name_head = os.path.splitext(uploaded_file.name)[0]
        st.session_state.deseq2_uploaded_file_name = file_name_head

    else:
        sys.exit(1)
   ##### file upload

if df is not None:

    df.iloc[0:3,:]

    content = df.columns.tolist()
    pvalue = [i for i in content if ('pvalue' in i) or ('p-val' in i) or ('p val' in i) or ('padj' in i) or ('p_val' in i)  or ('adj' in i)]
    fc = [i for i in content if ('log2FC' in i) or ('Fold Change' in i) or ('log2FoldChange' in i) or ('logFC' in i) or ('coef' in i)]
    gene = [i for i in content if (i not in pvalue) and (i not in fc)]
    P_column = st.selectbox(
        'Select adjusted P-value column',
        pvalue)
    # ジャロ・ウィンクラー距離法
    JW_dist = [Levenshtein.jaro_winkler(P_column, x) for x in fc]
    try:
        FC_column = st.selectbox(
            'Select FC column',
            fc, index = JW_dist.index(max(JW_dist)))
    except:
        FC_column = st.selectbox(
            'Select FC column', fc)

    file_name_add = P_column
    file_name_add = file_name_add.replace('.adj.pvalue','')
    file_name_add = file_name_add.replace(' adj. p-value','')
    file_name = file_name_add

    # Modified section for the inversion handling
    # Store original data when loading new file or initializing
    if "original_fc_data" not in st.session_state or use_upload == 'Yes':
        st.session_state.original_fc_data = df[FC_column].copy()
        st.session_state.inv = False  # Reset inversion state for new file

    inv_switch = st.checkbox('Flip reference (Reverse comparison)?', value=st.session_state.inv)
    if inv_switch:
        file_name = file_name + "_ReferenceFlipped"

    # Handle sign inversion
    if inv_switch != st.session_state.inv:  # Only update when switch changes
        if inv_switch:
            df[FC_column] = -st.session_state.original_fc_data
        else:
            df[FC_column] = st.session_state.original_fc_data.copy()  # Use copy to avoid reference issues
        st.session_state.inv = inv_switch

 #   st.write(st.session_state.inv)

    df.iloc[0:3,:]
 #   st.write(st.session_state.inv)



    set_gene = st.checkbox("Specify gene column?", value=False)

    if set_gene:
        Gene_column =  st.selectbox(
        'Select gene column',
        gene)
    else:

        if "Annotation/Divergence" in content:
            pattern = "([^|]*)"
            Gene_column = 'Annotation/Divergence'
            df.loc[:,'Annotation/Divergence'] = df.apply(lambda x: re.sub(r'([^|]*).*', r'\1', x['Annotation/Divergence']), axis=1)
            st.write("Converted Annotation/Divergence to gene symbols.")
        elif "Gene" in content:
            Gene_column =  "Gene"
        else:
            Gene_column =  st.selectbox(
            'Select gene column',
            gene)



    if use_upload == 'Yes':
        # indexをGeneにする
        df.index = df[Gene_column].tolist()

    st.markdown("#### Use R EnhancedVolcano?")
    enhanced =  st.checkbox('Use R EnhancedVolcano?', label_visibility = 'collapsed')
    if not enhanced:
        transform = st.radio(
            "Y axis transformation:",
            ('-Log10','-Log2', 'As is'), index =0)

        if (P_column is not None) & (FC_column is not None):
            if transform == "-Log10":
                df["Y_data"] = -np.log10(df[P_column])
                ylabel = '-Log10 P'

            elif transform == "-Log2":
                df["Y_data"] = -np.log2(df[P_column])
                ylabel = '-Log2 P'
            else:
                df["Y_data"] = df[P_column]
                ylabel = 'P'

# pが0のときの対応
    zero_adj = "No"
    if min(df[P_column]) == 0:
        st.markdown("##### There are zeros in P values.")
        zero_adj = st.radio(
            "Assing a minimal value (transformed miniman value + 1) to zeros for visualization?",
            ('Yes','No'))

    if zero_adj == "Yes":
        if (P_column is not None) & (FC_column is not None):
            # maxにするとinfが結果になるため、まずnanに変換した後、np.nanmaxを使う
            df.replace([np.inf, -np.inf], np.nan, inplace=True)
            min_val = np.nanmax(df["Y_data"])
            if transform == "-Log10" or "-Log2":
                df.loc[(df[P_column]==0),"Y_data"] = min_val + 1
                st.markdown("###### Minimal P values have been converted to :" + str(min_val + 1))
            else:
                df["Y_data"] = df[P_column]
                ylabel = 'P'

    blue_c = '#D3D3D3'
    up_c = '#ef553b'
    down_c = '#636efa'
#    red_c = "#ef553b"
    Pt_size = 4
    Pt_alpha = 1
    blue_alpha = 0.8
    up_alpha = 1
    down_alpha = 1

    if enhanced:
        with st.sidebar:
            st.markdown("#### EnhancedVolcano option")
            graph_format = st.radio('Output format', ['png', 'pdf'])
            num_label = str(int(st.number_input('Max number of genes to label',min_value=0, value = 50)))
            EV_X_size = str(st.number_input("Plot X size", min_value=1.0, value=10.0))
            EV_Y_size = str(st.number_input("Plot Y size", min_value=1.0, value=11.0))
            change_limit = st.checkbox('Change axis limits?', value = False)
            if change_limit:
                EV_X_min = str(st.number_input("Min of X axis ", value=-5.0))
                EV_X_max = str(st.number_input("Max of X axis ", value=5.0))
                EV_Y_min = str(st.number_input("Min of Y axis ",  value=0.0))
                EV_Y_max = str(st.number_input("Max of Y axis ", value=5.0))
                EV_limits = ", xlim = c(" + EV_X_min + ", " + EV_X_max + "), ylim = c(" + EV_Y_min + ", " + EV_Y_max + ")"
            dConnectors = st.checkbox('Draw connectors for gele labels?', value = True)
            if dConnectors:
                wConnectors = str(st.number_input("Connector width", min_value=0.10, value=0.40))
            st.markdown("##### Gene point attributes")
            pSize = str(st.number_input("Point size", min_value=0.5, value=2.0))
            cAlpha = str(st.number_input("Point transparency (alpha)", min_value=0.0, value=1.0, max_value =1.0))
            lSize = str(st.number_input("Gene label size", min_value=0.5, value=5.0))
            EV_color = st.checkbox('Change point colors?', value = False)
            if EV_color:
                bg_gray = st.checkbox('Set all backgroud gray?', value = False)
                if not bg_gray:
                    B_c = st.color_picker('Positive point color', '#ee0000')
                    P_c = st.color_picker('Only P positive color', '#4169e1')
                    FC_c = st.color_picker('Only FC positive color', '#228b22')
                    NS_c = st.color_picker('Background point color', '#4d4d4d')
                else:
                    B_c = st.color_picker('Positive point color', '#ee0000')
                    NS_c = st.color_picker('Background point color', '#D3D3D3')
                    P_c = NS_c
                    FC_c = NS_c
                color_codes = "col = c('" + NS_c + "', '"  +  FC_c  + "', '"  +  P_c + "', '"  + B_c + "')"
                st.write(color_codes)
            st.markdown("##### Legend attributes")
            lLabSize = str(st.number_input("Legend font size", min_value=1.0, value=12.0))
            lIconSize = str(st.number_input("Legend icon size", min_value=1.0, value=6.0))


    else:
        with st.sidebar:

            st.markdown("##### Modify point size, color, transparency?")
            Mod_point = st.checkbox('mod point', label_visibility = 'collapsed')
            if Mod_point:
                Pt_size = float(st.text_input("Point size:", value = 4))
                st.markdown("##### NS point color:")
                blue_c = st.color_picker('', '#D3D3D3', label_visibility="collapsed")
                st.write("ex: gray #808080  light gray #D3D3D3  blue gray #7393B3  sage green #8A9A5B")
                blue_alpha = st.slider('NS point transparency (alpha) 0-1:', 0.0, 1.0, 0.8, 0.05)
                st.markdown("##### Significantly upregulated point color:")
                up_c = st.color_picker('Significant point color', '#ef553b', label_visibility="collapsed")
                st.write("ex: burgandy #800020  bordeaux #4C0013")
                up_alpha = st.slider('Upregulated point transparency (alpha) 0-1:', 0.0, 1.0, 1.0, 0.05)
                st.markdown("##### Significant point color:")
                down_c = st.color_picker('Significantly downregulated point color', '#636efa', label_visibility="collapsed")
                st.write("ex: cobalt blue #0047AB  blue #0000FF")
                down_alpha = st.slider('Downregulated point transparency (alpha) 0-1:', 0.0, 1.0, 1.0, 0.05)
                st.markdown("---")
            st.markdown("##### Modify plot ranges?")

            Mod_range = st.checkbox('range', label_visibility = 'collapsed')
            x_min = -float("inf")
            x_max = float("inf")
            y_max = float("inf")
            x_axis_mod = False
            y_axis_mod = False
            if Mod_range:
                x_min = float(st.text_input("X-min:", value = math.floor(min(df[FC_column]))-1))
                x_max = float(st.text_input("X-max:", value = math.floor(max(df[FC_column]))+1))
                y_max = float(st.text_input("Y-max:", value = float("inf")))
                x_axis_mod = True
                if  y_max != float("inf"):
                    y_axis_mod = True
                st.markdown("---")

            blue_rgb =  hex2rgb(blue_c, blue_alpha)
            up_rgb = hex2rgb(up_c, up_alpha)
            down_rgb = hex2rgb(down_c, down_alpha)

            st.markdown("##### Figure Y size:")
            Y_size = int(st.text_input("", value = 700, label_visibility = 'collapsed'))

    with st.sidebar:
        st.markdown("##### Coloring thresholds:")
        P_threshold = float(st.text_input("P threshold:", value = '0.05'))
        FC_threshold = float(st.text_input("LFC threshold:", value = '1'))
    L_P_threshold = P_threshold
    L_FC_threshold = FC_threshold

   # if not enhanced:
    with st.sidebar:
        st.markdown("##### Gene Name Labeling thresholds:")
        L_P_threshold = float(st.text_input("Labeling P threshold:", value = '0.05'))
        L_FC_threshold = float(st.text_input("Labeling LFC threshold:", value = '1'))
        st.markdown("---")

    df['Color2'] = '#D3D3D3'
    if not enhanced:
        with st.sidebar:
            Font_size = 12
            Font_color = "#000000"
            Font_alpha = 1
            Font_position = 'top center'
            Font_italic = False  # デフォルトは非イタリック
            Mod_label = st.checkbox('Modify gene label attributes?')
            if Mod_label:
                Font_size = int(st.text_input("Gene label size:", value = 12))
                Font_alpha = st.slider('Gene label transparency (alpha) 0-1:', 0.0, 1.0, 1.0, 0.05)
                Font_color = st.color_picker('Gene label color', '#000000')
                Font_italic = st.checkbox('Italic font?', value=False)  # イタリックのオプションを追加
                Font_position = st.radio("Label position:", ('top center','top left', 'top right'))
                st.markdown("---")
            Font_rgb = hex2rgb(Font_color, Font_alpha)
            # フォントファミリーの設定
     #       Font_family = "Arial Italic, Times New Roman Italic" if Font_italic else "Arial, Times New Roman"


         #   df['Color2'] = blue_rgb
            df['Color2'][(df[P_column] < P_threshold) & (df[FC_column] > FC_threshold)] = up_rgb
            df['Color2'][(df[P_column] < P_threshold) & (df[FC_column] < -FC_threshold)] = down_rgb

            # colorに合わせてデータの順番をソート
            df2 = df.copy()
         #   st.write("Before sort:" + str(len(df2)))
            up_rows = df2.loc[df2['Color2'] == up_rgb]
          #  st.write(up_rows.head())
            df2.drop(up_rows.index, inplace=True)
            df2 = pd.concat([df2, up_rows])
         #   st.write("After up:" + str(len(df2)))
            down_rows = df2.loc[df2['Color2'] == down_rgb]
            df2.drop(down_rows.index, inplace=True)
            df2 = pd.concat([df2, down_rows])
          #  st.write("After sort:" + str(len(df2)))

            df = df2



            P_line = 0
            FC_line = 0
            PFC_width = 0.4
            PFC_color = '#000000'
            PFC_lines = st.checkbox('Plot P and FC threshold lines?')
            if  PFC_lines:
                P_line = float(st.text_input("P line (set 0 to remove):", value = str(P_threshold)))
                FC_line = float(st.text_input("FC line (set 0 to remove):", value = str(FC_threshold)))
                PFC_width = float(st.text_input("Line width:", value = 0.4))
                PFC_color = st.color_picker('Axis font color', '#000000')
                if transform == "-Log10":
                    P_line = -np.log10(P_line)
                elif transoform == "-Log2":
                    P_line = -np.log2(P_line)
                st.markdown("---")


    text_label = df[Gene_column].copy()
 #   text_label[(df[P_column] >= L_P_threshold) | (np.abs(df[FC_column]) <= L_FC_threshold)] = None

    text_label[((L_P_threshold > 0) & (df[P_column] >= L_P_threshold)) | 
           ((L_FC_threshold > 0) & (np.abs(df[FC_column]) <= L_FC_threshold))] = None


    # adj_pがNAのものが多い
    text_label[ np.isnan(df[P_column])] = None
    gene_list = [x for x in text_label if x is not None]

    if use_goi and goi is not None: # goiで制限
        gene_list = list(set(gene_list) & set(goi))
        # Update text_label to match gene_list
        text_label_filtered = df[Gene_column].copy()
        text_label_filtered[:] = None  # Reset all labels to None
        # Only set labels for genes in gene_list
        for gene in gene_list:
            if gene in df.index:
                text_label_filtered.loc[gene] = gene
        text_label = text_label_filtered

    with st.sidebar:
        Sp_genes = st.checkbox('Specify the genes to label?')
        if Sp_genes:
            st.markdown("##### Genes (comma or space  or tab separated):")
            genes = st.text_input("",label_visibility = 'collapsed')

            Sp_point_size = st.number_input("Point size multiplier for specified genes:", min_value=1.0, value=1.2, step=0.1)
            Sp_sig_point_size = st.number_input("Additional size multiplier for significant specified genes:", min_value=1.0, value=1.4, step=0.1)
            st.markdown("*Final size for significant specified genes will be: Base size × First multiplier × Second multiplier*")


            # 指定された遺伝子のリスト作成部分を修正
            gene_list_sp = []
            if len(genes) > 0:
                text_label = df[Gene_column].copy(deep=True)
                genes = genes.replace("'","")
                genes = genes.replace('"',"")
                gene_list_sp = genes.split(' ')
                gene_list_sp = list(filter(lambda a: a != '', gene_list_sp))
                if ',' in genes:
                    gene_list_sp = sum([x.split(',') for x in gene_list_sp],[])
                if '\t' in genes:
                    gene_list_sp = sum([x.split('\t') for x in gene_list_sp],[])
                if '\n' in genes:
                    gene_list_sp = sum([x.split('\n') for x in gene_list_sp],[])

                gene_set_lower = [str(x).lower() for x in set(gene_list_sp)]
                text_label_dat = []
                
                # Create point size array
                point_sizes = [Pt_size] * len(df)  # Default size for all points
                
                for idx, x in enumerate(text_label):
                    if x != None:
                        if str(x).lower() in gene_set_lower:
                            text_label_dat.append(x)
                            # Base multiplier for specified genes
                            point_sizes[idx] = Pt_size * Sp_point_size
                            # Additional multiplier for significant specified genes
                            if df.iloc[idx][P_column] < P_threshold and abs(df.iloc[idx][FC_column]) > FC_threshold:
                                point_sizes[idx] = point_sizes[idx] * Sp_sig_point_size
                        else:
                            text_label_dat.append(None)
                    else:
                        text_label_dat.append(None)
                        
                text_label.iloc[:] = text_label_dat
                gene_list_all = [x for x in text_label if x is not None]  # ここで gene_list_all を定義
                df['Point_Size'] = point_sizes

                text_label_color = text_label.copy() #色付け用
                text_label_color[(df[P_column] >= P_threshold) | (np.abs(df[FC_column]) <= FC_threshold)] = None
                gene_list_color = [x for x in text_label_color if x is not None]


                text_label[(df[P_column] >= L_P_threshold) | (np.abs(df[FC_column]) <= L_FC_threshold)] = None
                # adj_pがNAのものが多い
                text_label[ np.isnan(df[P_column])] = None
                gene_list = [x for x in text_label if x is not None]

                Sp_attribute = st.checkbox('Change the point color of the specified genes?')
                Sp_alpha = 1
                Sp_color = "#ef553b"

                if Sp_attribute:
                    Sp_all = st.checkbox('Color all specified genes irrespective of FDR/FC thresholds?')
                    if Sp_all:
                        gene_list_color = gene_list_all
                    Sp_c = st.color_picker('Color:', '#4C0013', label_visibility="collapsed")
                    Sp_alpha = st.slider('NS point transparency:', 0.0, 1.0, 1.0, 0.05)
                    Sp_sig_alpha = 1.0  # 有意な遺伝子用の透明度を1.0（完全に不透明）に固定

                    # 非有意な指定遺伝子用のRGBA色
                    Sp_rgb_ns = hex2rgb(Sp_c, Sp_alpha)
                    # 有意な指定遺伝子用のRGBA色
                    Sp_rgb_sig = hex2rgb(Sp_c, Sp_sig_alpha)

                    for i in gene_list_color:
                        # 有意性に基づいて適切な透明度を選択
                        is_significant = (df.loc[i, P_column] < P_threshold) and (abs(df.loc[i, FC_column]) > FC_threshold)
                        df.loc[df[Gene_column] == i, 'Color2'] = Sp_rgb_sig if is_significant else Sp_rgb_ns

                    df2 = df.copy()
                    # 有意な指定遺伝子を最前面に
                    sig_rows = df2.loc[(df2['Color2'] == Sp_rgb_sig)]
                    df2.drop(sig_rows.index, inplace=True)
                    df2 = pd.concat([df2, sig_rows])
                    
                    # 非有意な指定遺伝子を次に
                    ns_rows = df2.loc[(df2['Color2'] == Sp_rgb_ns)]
                    df2.drop(ns_rows.index, inplace=True)
                    df2 = pd.concat([df2, ns_rows])

                    text_label = df2[Gene_column].copy()
                    text_label_dat = []
                    for x in text_label:
                        if x != None:
                            if str(x).lower() in gene_set_lower:
                                text_label_dat.append(x)
                            else:
                                text_label_dat.append(None)
                        else:
                            text_label_dat.append(None)
                    text_label.iloc[:] = text_label_dat
                    text_label[(df[P_column] >= L_P_threshold) | (np.abs(df[FC_column]) <= L_FC_threshold)] = None
                    # adj_pがNAのものが多い
                    text_label[ np.isnan(df[P_column])] = None
                    gene_list = [x for x in text_label if x is not None]
                    df = df2.copy()

            st.markdown("---")

    if enhanced:
        with st.sidebar:
            Hl_gene_option = st.checkbox('Highlight specific genes?')
            if Hl_gene_option:
                st.markdown("##### Genes (comma or space  or tab separated):")
                Hl_genes = st.text_input("highlight genes",label_visibility = 'collapsed')
                Hl_gene_list_sp = []
                if len(Hl_genes) > 0:
                    Hl_genes = Hl_genes.replace("'","")
                    Hl_genes = Hl_genes.replace('"',"")
                    st.write(Hl_genes)
                    Hl_gene_list_sp = Hl_genes.split(' ') #まず空白で分離
                    Hl_gene_list_sp = list(filter(lambda a: a != '', Hl_gene_list_sp)) #空白のみを除く
                    if ',' in Hl_genes:
                        Hl_gene_list_sp = sum([x.split(',') for x in Hl_gene_list_sp],[]) #sumで平坦化 sum(x, [])
                    if '\t' in Hl_genes:
                        Hl_gene_list_sp = sum([x.split('\t') for x in Hl_gene_list_sp],[])
                    if '\n' in Hl_genes:
                        Hl_gene_list_sp = sum([x.split('\n') for x in Hl_gene_list_sp],[])
              #      st.write(text_label)
         #           text_label[-text_label.isin(gene_list)] = None
                    Hl_gene_set_lower = [x.lower() for x in set(Hl_gene_list_sp)] #大文字小文字関係なくサーチ
                    Hl_text_label_dat = []
                    for x in df[Gene_column]:
                        if x != None:
                            if x.lower() in Hl_gene_set_lower:
                                Hl_text_label_dat.append(x)
                            else:
                                Hl_text_label_dat.append(None)
                        else:
                            Hl_text_label_dat.append(None)
                    Hl_gene_list = [x for x in Hl_text_label_dat if x is not None]
                    Hl_pSize = str(st.number_input("Highlighted point size", min_value=0.5, value=float(pSize) + 4.0))



    if not enhanced:
        with st.sidebar:
            Axis_size = 16
            Tick_size = 12
            Axis_color = '#000000'
            Tick_color = '#000000'
            Mod_label = st.checkbox('Modify legend and tick fonts?')
            if Mod_label:
                Axis_size = int(st.text_input("Axis font size:", value = 16))
                Tick_size = int(st.text_input("Tick font size:", value = 12))
                Axis_color = st.color_picker('Axis font color', '#000000')
                Tick_color = st.color_picker('Tick font color', '#000000')


            repel = st.checkbox('Make a figure with gene name position adusted?')
            if repel:
                st.markdown("#### This may take a long time for many genes. Advised for <30 genes.")
                py_x_size = float(st.text_input("Plot x size:", value = 14))
                py_y_size = float(st.text_input("Plot y size:", value = 14))


    st.markdown("---")
    st.write("Number of genes to label: " + str(len(gene_list)))

    if st.button('Make plot'):
        if not enhanced:

            fig = go.Figure()
            data_o = go.Scatter(
             x=df[FC_column],
             y=df['Y_data'],
             name='test',
             hovertext=list(df[Gene_column]),
            mode='markers',
            marker_color=df['Color2'],
            marker_size=df['Point_Size'] if 'Point_Size' in df.columns else Pt_size
            )


            fig.add_trace(data_o)
            if FC_line != 0:
                fig.add_vline(x=FC_line, line_width=PFC_width, line_color = PFC_color)
                fig.add_vline(x=-FC_line, line_width=PFC_width, line_color = PFC_color)
            if P_line != 0:
                fig.add_hline(y=P_line, line_width=PFC_width, line_color = PFC_color)

            fig.update_layout(height=Y_size) #サイズの設定
            fig.update_layout(xaxis=dict(title=FC_column, titlefont = dict(size = Axis_size, color=Axis_color), tickfont=dict(size=Tick_size, color=Tick_color)),yaxis = dict(title =ylabel , titlefont = dict(size = Axis_size, color=Axis_color), tickfont=dict(size=Tick_size, color=Tick_color)))
            if x_axis_mod:
                fig.update_layout(xaxis_range = [x_min,x_max])
            if y_axis_mod:
                fig.update_layout(yaxis_range = [0,y_max])

            st.plotly_chart(fig)


            fig3 = go.Figure()

            # フォントファミリーの設定（Arialのみに変更）とstyle-italicを追加
            if Font_italic:
                text_label = [f'<i>{label}</i>' if label is not None else None for label in text_label]

            data = go.Scatter(
                x=df[FC_column],
                y=df['Y_data'],
                name='test',
                hovertext=list(df[Gene_column]),
                text=text_label,
                textfont=dict(
                    size=Font_size, 
                    color=Font_rgb, 
                    family="Arial"
                ),
                mode='markers+text',
                marker_color=df['Color2'],
                marker_size=df['Point_Size'] if 'Point_Size' in df.columns else Pt_size
            )


            fig3.add_trace(data)
            fig3.update_traces(textposition=Font_position)
            if FC_line != 0:
                fig3.add_vline(x=FC_line, line_width=PFC_width, line_color = PFC_color)
                fig3.add_vline(x=-FC_line, line_width=PFC_width, line_color = PFC_color)
            if P_line != 0:
                fig3.add_hline(y=P_line, line_width=PFC_width, line_color = PFC_color)


            fig3.update_layout(height=Y_size) #サイズの設定
            fig3.update_layout(xaxis=dict(title=FC_column, titlefont = dict(size = Axis_size, color=Axis_color), tickfont=dict(size=Tick_size, color=Tick_color)),yaxis = dict(title =ylabel , titlefont = dict(size = Axis_size, color=Axis_color), tickfont=dict(size=Tick_size, color=Tick_color)))
            if x_axis_mod:
                fig3.update_layout(xaxis_range = [x_min,x_max])
            if y_axis_mod:
                fig3.update_layout(yaxis_range = [0,y_max])

            st.plotly_chart(fig3)

            fig.write_image(res_dir + "/" +  file_name + "_Volcano_no_label.pdf")
            fig3.write_image(res_dir + "/" +  file_name + "_Volcano.pdf")

            if repel:
                import matplotlib.pyplot as plt
                from adjustText import adjust_text
                # NSとSigのデータを作る
                NS_df = df.copy()
                S_df = df.copy()
                NS_df = NS_df[ (NS_df[P_column]>=P_threshold) | (np.abs(NS_df[FC_column]) <= FC_threshold) ]
                up_df = S_df[ (S_df[P_column]<P_threshold) & (df[FC_column] > FC_threshold)]
                down_df = S_df[ (S_df[P_column]<P_threshold) & (df[FC_column] < -FC_threshold)]
        #        st.write(up_df)
        #        st.write(down_df)
                if len(gene_list) > 70:
                    st.markdown("### !!! Too many genes to label !!!")
                    st.write(str(len(gene_list)) + " genes.")
                else:
                    xns = NS_df[FC_column]
                    yns = NS_df['Y_data']
                    xup = up_df[FC_column]
                    yup = up_df['Y_data']
                    xdown = down_df[FC_column]
                    ydown = down_df['Y_data']

                    fig4 = plt.figure(figsize=(py_x_size, py_y_size))
                    py_blue_c = pyhex2rgb(blue_c, blue_alpha)
                    py_up_c = pyhex2rgb(up_c, up_alpha)
                    py_down_c = pyhex2rgb(down_c, down_alpha)
                    plt.scatter(xns, yns, c=py_blue_c, edgecolor=(1,1,1,0),  s=Pt_size*1.5)
                    plt.scatter(xup, yup, c=py_up_c, edgecolor=(1,1,1,0),  s=Pt_size*1.5)
                    plt.scatter(xdown, ydown, c=py_down_c, edgecolor=(1,1,1,0),  s=Pt_size*1.5)
                    texts = []
                    if len(gene_list) >0:
                        xl = []
                        yl = []
                        for i in gene_list:
                            xl.append(df.loc[i, FC_column])
                            yl.append(df.loc[i, 'Y_data'])
                        for x, y, l in zip(xl, yl, gene_list):
                            texts.append(plt.text(x, y, l, size=Font_size*0.6))
                        adjust_text(texts, arrowprops=dict(arrowstyle="-", color='k', lw=0.4))
                    plt.gca().spines['right'].set_visible(False)
                    plt.gca().spines['top'].set_visible(False)
                    plt.xlabel(FC_column)
                    plt.ylabel(ylabel)
                    if x_axis_mod:
                        plt.xlim([x_min,x_max])
                    if y_axis_mod:
                        plt.ylim([0,y_max])

                    if FC_line != 0:
                        plt.axvline(x=FC_line, ymin=0, ymax=1, lw=PFC_width, color = PFC_color )
                        plt.axvline(x=-FC_line, ymin=0, ymax=1, lw=PFC_width, color = PFC_color)
                    if P_line != 0:
                        plt.axhline(y=P_line, xmin=0, xmax=1, lw=PFC_width, color = PFC_color)
                    st.pyplot(fig4)

                    fig4.savefig(res_dir + "/" +  file_name + "_Volcano_arrow_labels.pdf")

        if enhanced:
            if use_upload == 'Yes' or "pyper" not in st.session_state:
                st.session_state.pyper = False
            if not st.session_state.pyper:
                st.session_state.pyper = True
                df.to_csv('pyper_temp.csv')
            import pyper
            r = pyper.R(use_pandas='True')
            r("library(EnhancedVolcano)")
            st.session_state.pyper = True
#                r.assign("p_df", df)
            r('options(ggrepel.max.overlaps = Inf)')
            r("p_df <- read.csv('pyper_temp.csv', row.names =1, check.names = FALSE)")
            r.assign('graph_format', graph_format)
            r.assign('EV_Y_size', str(EV_Y_size))
            r.assign('EV_X_size', str(EV_X_size))

            if len(gene_list) > int(num_label):
                df=df.sort_values(P_column, ascending=True)
                gene_list = df.index.tolist()[:int(num_label)]

            select_genes =  ", " +  "selectLab = c(" + '"' + '", "'.join(gene_list) +  '")'
#            st.write(select_genes)


            r_enhanced_command_head = "p <- EnhancedVolcano(p_df, lab = rownames(p_df), x = '" + FC_column + "', y = '" + P_column + "'" + select_genes
            EV_options = ""

            xlab =  ", xlab = bquote(~Log[2]~ 'fold change')"
            pCutoff = ", pCutoff = " + str(P_threshold)
            FCCutoff =  ", FCcutoff = " + str(FC_threshold)
            if Hl_gene_option:
                pt_size_list =  [Hl_pSize if x in Hl_gene_list else pSize for x in df.index.to_list()]
                pt_size_list_str = "c("+ ", ".join(pt_size_list)+ ")"
#                r.assign('pt_size_list', pt_size_list)
                pointSize = ", pointSize = " + pt_size_list_str
            else:
                pointSize = ", pointSize = " + pSize
            labSize  = ", labSize = " + lSize
            colAlpha =  ", colAlpha = " +  cAlpha
            legendLabSize =  ", legendLabSize = " + lLabSize
            legendIconSize =  ", legendIconSize = " + lIconSize
            if dConnectors:
                drawConnectors =  ", drawConnectors = TRUE"
            else:
                drawConnectors =  ", drawConnectors = FALSE"
            widthConnectors =  ", widthConnectors = " + wConnectors
            if EV_color:
                EV_options = EV_options + ", " + color_codes

            if change_limit:
                EV_options = EV_options + EV_limits
 
            r_enhanced_command = r_enhanced_command_head + xlab + pCutoff + FCCutoff + pointSize + labSize + colAlpha + legendLabSize + legendIconSize + drawConnectors + widthConnectors + EV_options + ")"
            r(r_enhanced_command)

            # まずPNGを保存して表示
            png_command = "ggsave(filename = '" + res_dir + "/" + file_name + "_EnhancedVolcano.png', plot = p, device = 'png', width = " + EV_X_size + ", height = " + EV_Y_size + ", dpi = 300)"
            r(png_command)
            st.image(res_dir + "/" + file_name + "_EnhancedVolcano.png")

            # PDFも保存
            pdf_command = "ggsave(filename = '" + res_dir + "/" + file_name + "_EnhancedVolcano.pdf', plot = p, device = 'pdf', width = " + EV_X_size + ", height = " + EV_Y_size + ")"
            r(pdf_command)

        # zipファイルの作成と保存
        shutil.make_archive(temp_dir + "/Volcano", format='zip', root_dir=res_dir)
        with open(temp_dir + "/Volcano.zip", "rb") as fp:
            btn = st.download_button(
                label="Download plots?",
                data=fp,
                file_name=file_name_head + "_" + file_name + "_Volcano.zip",
                mime="zip"
            )

        # 一時ディレクトリの掃除
        shutil.rmtree(temp_dir)
        os.mkdir(temp_dir)
