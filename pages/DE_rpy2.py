import streamlit as st
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri, Formula
from rpy2.robjects.packages import importr
from rpy2.robjects import StrVector
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import io
import os
import time
from helper_func import clear_old_directories, clear_old_files, make_r_compatible_column_names, remove_after_space, remove_sample_num
import re
from itertools import combinations
import zipfile
import sys
from collections import Counter

pandas2ri.activate()

# import R libraries
DESeq2 = importr('DESeq2')
edgeR = importr('edgeR')
Limma = importr('limma')
stats = importr('stats')
base = importr('base')

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
#    return(dfx)

st.set_page_config(page_title="Compare DESeq2, edgeR, limma-voom", page_icon="📃")


@st.cache_data
def convert_df(df):
   return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def read_excel(file, index_col=None, header = 0):
    df_xl = pd.read_excel(file, index_col = index_col, header = header)
    return df_xl

@st.cache_data
def read_csv(file, index_col=None, sep=','):
    df_c = pd.read_csv(file, index_col = index_col, header = 0, sep = sep, engine='python')
    return df_c

def create_venn_diagram(set1, set2, set3, labels):
    plt.figure(figsize=(10, 10))
    venn3([set(set1), set(set2), set(set3)], set_labels=labels)
    plt.title("Venn Diagram of significant genes")
    return plt

class DE_rpy2:
    def __init__(self, count_matrix_path, design_matrix_path, design_formula=None, gene_column='id', temp_dir='temp'):
        self.count_matrix_path = count_matrix_path
        self.design_matrix_path = design_matrix_path
        self.temp_dir = temp_dir
        
        if not os.path.exists(self.temp_dir):
            os.makedirs(self.temp_dir)

        self.count_matrix = pd.read_csv(count_matrix_path, sep='\t', index_col=0)
        self.design_matrix = pd.read_csv(design_matrix_path, sep='\t', index_col=0)

        self.gene_ids = self.count_matrix.index
        self.gene_column = gene_column

        if design_formula is None:
            condition = self.design_matrix.columns[0]
            self.design_formula = Formula('~ ' + condition)
        else:
            self.design_formula = Formula(design_formula)

        self.deseq2_dds = None
        self.edger_fit = None
        self.limma_fit = None

    def run_deseq2(self):
        ro.r(f'''
        count_data <- as.matrix(read.csv("{self.count_matrix_path}", sep="\t", row.names=1, check.names=FALSE))
        count_data <- round(count_data) #四捨五入する
        coldata <- read.csv("{self.design_matrix_path}", sep="\t", row.names=1, check.names=FALSE)
        stopifnot(all(colnames(count_data) == rownames(coldata)))
        
        dds <- DESeqDataSetFromMatrix(countData = count_data,
                                      colData = coldata,
                                      design = {self.design_formula.r_repr()})
        dds <- DESeq(dds)
        
        saveRDS(dds, file="{os.path.join(self.temp_dir, 'deseq2_dds.rds')}")
        ''')
        
        self.deseq2_dds = ro.r(f'readRDS("{os.path.join(self.temp_dir, "deseq2_dds.rds")}")')

    def get_deseq2_results(self, group1, group2, threshold=0.05):
        ro.r(f'''
        dds <- readRDS("{os.path.join(self.temp_dir, 'deseq2_dds.rds')}")
        res <- results(dds, contrast=c("condition", "{group2}", "{group1}"))
        write.csv(as.data.frame(res), file="{os.path.join(self.temp_dir, f'deseq2_results_{group1}_vs_{group2}.csv')}")
        ''')

        result = pd.read_csv(os.path.join(self.temp_dir, f'deseq2_results_{group1}_vs_{group2}.csv'), index_col=0)
        result[self.gene_column] = self.gene_ids
        label = result.index[result['padj'] < threshold].tolist()

        return result, label

    def run_edger(self):
        ro.r(f'''
        count_data <- as.matrix(read.csv("{self.count_matrix_path}", sep="\t", row.names=1, check.names=FALSE))
        coldata <- read.csv("{self.design_matrix_path}", sep="\t", row.names=1, check.names=FALSE)
        stopifnot(all(colnames(count_data) == rownames(coldata)))

        y <- DGEList(counts=count_data, group=coldata$condition)
        y <- calcNormFactors(y)
        condition <- factor(coldata$condition)
        design <- model.matrix(~0 + condition, data=y$samples)
        print('design in R')
        print(design)

        y <- estimateDisp(y, design, robust =TRUE)
        fit <- glmQLFit(y, design)
        
        saveRDS(list(y=y, fit=fit, design=design), file="{os.path.join(self.temp_dir, 'edger_fit.rds')}")
        ''')
        
        self.edger_fit = ro.r(f'readRDS("{os.path.join(self.temp_dir, "edger_fit.rds")}")')

    def get_edger_results(self, group1, group2, threshold=0.05):
        ro.r(f'''
        edger_data <- readRDS("{os.path.join(self.temp_dir, 'edger_fit.rds')}")
        y <- edger_data$y
        fit <- edger_data$fit
        design <- edger_data$design
        
        contrast <- makeContrasts(contrasts=paste0("condition", "{group2}", "-", "condition", "{group1}"), levels=design)
        print(contrast)
        qlf <- glmQLFTest(fit, contrast=contrast)
        res <- topTags(qlf, n=Inf)
        write.csv(as.data.frame(res), file="{os.path.join(self.temp_dir, f'edger_results_{group1}_vs_{group2}.csv')}")
        ''')

        result = pd.read_csv(os.path.join(self.temp_dir, f'edger_results_{group1}_vs_{group2}.csv'), index_col=0)
        label = result.index[result['FDR'] < threshold].tolist()

        return result, label

    def run_limma(self):
        ro.r(f'''
        library(edgeR)
        library(limma)
        count_data <- as.matrix(read.csv("{self.count_matrix_path}", sep="\t", row.names=1, check.names=FALSE))
        coldata <- read.csv("{self.design_matrix_path}", sep="\t", row.names=1, check.names=FALSE)
        stopifnot(all(colnames(count_data) == rownames(coldata)))
        
        y <- DGEList(counts=count_data, group=coldata$condition)
        y <- calcNormFactors(y)
        condition <- factor(coldata$condition)
        l_design <- model.matrix(~0 + condition, data=y$samples)
        ''')

        # 選択された方法に基づいて異なるRコードを実行
        if which_voom != 'voomLmFit':
            ro.r(f'''
            # voom変換と品質重みの適用
            l_v <- {which_voom}(y, l_design)
            l_fit <- lmFit(l_v, l_design)
            l_fit <- eBayes(l_fit)
            saveRDS(list(fit=l_fit, design=l_design), file="{os.path.join(self.temp_dir, 'limma_fit.rds')}")
            ''')
        else:  # edgeR::voomLmFit
        #    st.write("using voomLmFit")
            ro.r(f'''
            l_fit <- voomLmFit(y, l_design)
            l_fit <- eBayes(l_fit)
            saveRDS(list(fit=l_fit, design=l_design), file="{os.path.join(self.temp_dir, 'limma_fit.rds')}")
            ''')
        
        self.limma_fit = ro.r(f'readRDS("{os.path.join(self.temp_dir, "limma_fit.rds")}")')

    def get_limma_results(self, group1, group2, threshold=0.05):
        ro.r(f'''
        limma_data <- readRDS("{os.path.join(self.temp_dir, 'limma_fit.rds')}")
        fit <- limma_data$fit
        design <- limma_data$design
        
        contrast <- makeContrasts(paste0("condition", "{group2}", "-", "condition", "{group1}"), levels=design)
        fit2 <- contrasts.fit(fit, contrast)
        fit2 <- eBayes(fit2)
        res <- topTable(fit2, number=Inf)
        write.csv(as.data.frame(res), file="{os.path.join(self.temp_dir, f'limma_results_{group1}_vs_{group2}.csv')}")
        ''')

        result = pd.read_csv(os.path.join(self.temp_dir, f'limma_results_{group1}_vs_{group2}.csv'), index_col=0)
        label = result.index[result['adj.P.Val'] < threshold].tolist()

        return result, label


# Set up temporary directory
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir = "temp/" + str(round(time.time()))
    if not os.path.exists('temp'):
        os.mkdir('temp')
    else:
        clear_old_directories("temp")
        clear_old_files("temp")
        if not os.path.exists(temp_dir):
            os.mkdir(temp_dir)
    st.session_state.temp_dir = temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)
else:
    temp_dir = st.session_state.temp_dir
    res_dir = temp_dir + '/res'
    st.session_state.res_dir = res_dir
    if not os.path.exists(temp_dir):
        os.mkdir(temp_dir)
        os.mkdir(res_dir)
    if not os.path.exists(res_dir):
        os.mkdir(res_dir)

st.sidebar.title("Options")
st.markdown('### Compare DESeq2, edgeR, limma-voom')
st.markdown("#### Upload raw count data")
st.markdown("##### Data format:")
file_type = st.radio(
    "",    ('auto', 'Homer','tsv','csv','excel'), index = 0, label_visibility = 'collapsed')
uploaded_file = st.file_uploader("Choose a file", type=['txt','tsv', 'csv', 'xls','xlsx'])

if uploaded_file is not None:

    if file_type == 'auto':
        try:
            df = read_csv(uploaded_file, sep = None)
            st.write("Uploaded file:")
            st.write(df.head())

            content = df.columns.tolist()
#                Gene_column = content[0]

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
              #  st.write(df.head())
                # annotation/divergence以前を除く
                df = df.loc[:,'Annotation/Divergence':]
              #  st.write(df.head())
                st.write("Converted Annotation/Divergence to gene symbols.")
            df = make_r_compatible_column_names(df)
            content = df.columns.tolist()
            content[0] = 'Gene'
            df.columns = content

     #       df.set_index("Gene", inplace = True)

        except:# excel
            df = read_excel(uploaded_file)
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

    elif file_type != 'excel':
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
        df = read_excel(uploaded_file)
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

    #特殊文字を除く
    df.columns = df.columns.str.replace('[^A-Za-z0-9]+', '_')
    df.columns = df.columns.str.replace('-', '_')

else:
    sys.exit(1)

if uploaded_file is not None:
    if len(df.index.values) != len(set(df.index.values)):
        df = rename_duplicates(df)
    #特殊文字を除く
    df.columns = df.columns.str.replace('[^A-Za-z0-9]+', '_')
    df.columns = df.columns.str.replace('-', '_')

    with st.sidebar:
        st.markdown("##### Filter out weakly-expressed genes:")
        independentFiltering = st.checkbox('Yes', value= True)
        st.markdown("""低発現遺伝子の除外はFDRの計算を改善する。
        filtered outされた遺伝子が多い場合、GSEA用に全遺伝子をランキングするためにはチェックを外す。""")
        min_threshold = 0
        max_threshold = 0
        if independentFiltering:

            st.markdown("##### Filter the genes > counts in all samples:")
            min_threshold = st.number_input("count minimum", value = 0, label_visibility = 'collapsed')
            min_threshold = int(min_threshold)
            st.markdown("##### Filter the genes > counts in at least one sample:")
            max_threshold = st.number_input("count max", value = 0, label_visibility = 'collapsed')
            max_threshold = int(max_threshold)

        which_voom = st.radio("Voom method: ", ['voom', 'voomWithQualityWeights', 'voomLmFit'] , index=2)
        st.write("voomWithQualityWeights for a dataset containing low quality data")
        st.write("voomLmFit for sparse data (e.g., scRNA)")

    sample_threshold = 0


    if min_threshold > 0:
        df = df[df.apply(min, axis=1) > min_threshold]
    if max_threshold > 0:
        df = df[df.apply(max, axis=1) > max_threshold]

    st.write('Filtered gene number:  ' + str(len(df)))

    if any(df.sum() <= sample_threshold): # count 0の列を除く
        st.markdown('#### There are the samples that have counts <= ' + str(sample_threshold))
        st.write(", ".join(df.columns[df.sum() <= sample_threshold].to_list()))
        st.write('They are removed. Now data are:')
        df = df.drop(df.columns[df.sum() <= sample_threshold].to_list(), axis = 1)
        st.write(df.head())

    condition = [str(i) for i in df.columns.tolist()] #error防止
    group_condition = [remove_after_space(x) for x in condition] #スペース以降を除く
    group_condition = [remove_sample_num(x) for x in group_condition] #末尾の数字を除く


    st.write(df.head())

    df_e = pd.DataFrame(group_condition, index = condition, columns = ["Group"])

    with st.form("input_groups and batch"):
        st.write('Set groups:')
    #    edited_df_e = st.experimental_data_editor(df_e)
        edited_df_e = st.data_editor(df_e)

        condition = edited_df_e.iloc[:,0].tolist()
        

        submitted = st.form_submit_button("Submit")
    st.write('Group: ' + ' '.join(condition))


    if (len(condition) != len(df.columns)):
            st.write("The number of group name does not match the data.")


    design_matrix = pd.DataFrame({'condition': group_condition}, index=df.columns)

    r_condition = StrVector(condition)
    ro.r.assign('condition', r_condition)

    count_matrix_path = os.path.join(temp_dir, "count_matrix.tsv")
    design_matrix_path = os.path.join(temp_dir, "design_matrix.tsv")

    def is_numeric(x):
        return isinstance(x, (int, float)) and not pd.isna(x)
    numeric_row = df.applymap(lambda x: is_numeric(x))
    non_numeric_rows = df.apply(lambda row: row.apply(lambda x: not is_numeric(x)).any(), axis=1)
    # 少なくとも1つの非数値要素を含む行を選択
    if len(numeric_row) < len(df):
        st.write("There are nonumberic rows.")
        st.write(df_with_non_numeric)
        df = df[numeric_row.all(axis=1)]
        st.write(f"Remaining genes: {len(df)}")

    df.to_csv(count_matrix_path, sep='\t')
    design_matrix.to_csv(design_matrix_path, sep='\t')

    # Add a slider for threshold selection
    threshold = st.slider("##### Select p-value threshold for Venn diagram", 0.0, 0.2, 0.05, 0.01)

    if st.button('Run Differential Expression Analysis'):
        with st.spinner('Calculating...'):
            de_analysis = DE_rpy2(count_matrix_path, design_matrix_path, temp_dir=temp_dir)
            
            # Run analysis on all data
            de_analysis.run_deseq2()
            de_analysis.run_edger()
            de_analysis.run_limma()
            
            # Get all combinations of conditions
            conditions = design_matrix['condition'].unique()
            condition_pairs = list(combinations(conditions, 2))
            
            # Create a ZipFile object
            zip_path = os.path.join(temp_dir, 'results.zip')
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for group1, group2 in condition_pairs:
                    st.write(f"### Analysis for {group1} vs {group2}")
                    
                    deseq2_result, deseq2_label = de_analysis.get_deseq2_results(group1, group2, threshold)
                    st.write("DESeq2 Results:")
                    st.write(deseq2_result.head())
                    deseq2_result.to_csv(os.path.join(temp_dir, f'DESeq2.{group1}_vs_{group2}.tsv'), sep='\t')
                    zipf.write(os.path.join(temp_dir, f'DESeq2.{group1}_vs_{group2}.tsv'), f'DESeq2.{group1}_vs_{group2}.tsv')

                    edger_result, edger_label = de_analysis.get_edger_results(group1, group2, threshold)
                    st.write("edgeR Results:")
                    st.write(edger_result.head())
                    edger_result.to_csv(os.path.join(temp_dir, f'edgeR.{group1}_vs_{group2}.tsv'), sep='\t')
                    zipf.write(os.path.join(temp_dir, f'edgeR.{group1}_vs_{group2}.tsv'), f'edgeR.{group1}_vs_{group2}.tsv')

                    limma_result, limma_label = de_analysis.get_limma_results(group1, group2, threshold)
                    st.write("limma Results:")
                    st.write(limma_result.head())
                    limma_result.to_csv(os.path.join(temp_dir, f'limma.{group1}_vs_{group2}.tsv'), sep='\t')
                    zipf.write(os.path.join(temp_dir, f'limma.{group1}_vs_{group2}.tsv'), f'limma.{group1}_vs_{group2}.tsv')

                    # Create and display Venn diagram
                    st.write(f'DEG: DEseq2, {len(deseq2_label)}; edgeR, {len(edger_label)}; limma, {len(limma_label)}')
                    fig = create_venn_diagram(deseq2_label, edger_label, limma_label, ['DESeq2', 'edgeR', 'limma'])
                    st.pyplot(fig)

            st.success("Analysis completed for all condition pairs!")
            
            # Provide download link for the zip file
            with open(zip_path, "rb") as fp:
                btn = st.download_button(
                    label="Download All Results",
                    data=fp,
                    file_name="de_analysis_results.zip",
                    mime="application/zip"
                )