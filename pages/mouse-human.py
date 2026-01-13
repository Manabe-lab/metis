import streamlit as st
import csv
import re
import pickle
import sys
import os
import pandas as pd
import itertools
import shutil
import numpy as np

def load_m2h():
    import pickle
    data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "mouse2human.dic")
    with open(data_path, mode='rb') as f:
        data = pickle.load(f)
    return data

def load_h2m():
    import pickle
    data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "human2mouse.dic")
    with open(data_path, mode='rb') as f:
        data = pickle.load(f)
    return data


def mouse_human_conversion(dic, x):
    try:
        y = dic[x]
        return y
    except:
        return None

@st.cache_data
def convert_human_to_mouse_symbols(symbols, version=1): # R script from nichenetr converted with Claude3.5
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

@st.cache_data
def convert_mouse_to_human_symbols(symbols, version=1):
    if not isinstance(symbols, (list, pd.Series)):
        raise ValueError("symbols should be a list or pandas Series of mouse gene symbols")
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
    mousesymbol2humansymbol = dict(zip(geneinfo['symbol_mouse'], geneinfo['symbol']))
    converted_symbols = [mousesymbol2humansymbol.get(symbol, np.nan) for symbol in symbols]
    return converted_symbols


if os.path.exists('res'):
    shutil.rmtree("res")
os.mkdir("res") # Clear res first

# File type selection
file_type = st.radio(
    "Input file type",
    ('Gene list', 'Expression matrix'), index=0)

species = st.radio(
    "Species from",
    ('Mouse','Human'), index =1)

method = st.radio(
    "Conversion method",
    ('in-house','Nicehnetr v1', 'Nichenetr v2', 'Nicehnetr v1 (corrected)', 'Nichenetr v2 (corrected)', 'Consensus', 'Consensus (corrected)'), index = 3)

st.markdown("##### Methods:")
st.write("in-house: one to all orthologs mapping")
st.write("Nichnet: one to one")
st.write("v1: older symbols, v2: 2022 version. v1 may be better.")
#st.write("The mouth ortholog of CCL2 is set to Ccl2 in v1. CCL13 is also mapped to Ccl2. CCL2 is often mapped to Ccl12.")
st.write("Nicehnetr v1 (corrected): Ccl2→CCL2, Ccl3→CCL3, Cxcl1→CXCL1 corrected")
st.write("Nicehnetr v2 (corrected): v2 with correction")
st.write("Consensus: HomoloGene + Ensembl Compara consensus mapping")
st.write("Consensus (corrected): with correction")

# Help for creating consensus database
with st.expander("📚 How to create consensus database"):
    st.markdown("""
    ### How to create consensus ortholog mapping

    **1. Data sources**
    - **HomoloGene**: Conservative ortholog database provided by NCBI
    - **Ensembl Compara**: Phylogenetic ortholog database from Ensembl

    **2. Creation process**
    - Extract high-confidence mouse-human ortholog pairs from HomoloGene
    - Extract ortholog pairs with confidence score ≥ 75 from Ensembl Compara
    - Keep only pairs that match in both databases (intersection)
    - Keep only 1:1 mappings

    **3. Quality control**
    - Remove duplicate mappings
    - Standardize gene symbols
    - Remove non-standard gene names

    """)
st.markdown("---")

if file_type == 'Gene list':
    st.markdown("##### Upload gene list file or enter genes manually:")
    uploaded_file = st.file_uploader("Choose a file", type=['txt','csv','tsv'])
    
    st.markdown("##### Or Genes (comma, space, CR separated):")
    genes = st.text_input("genes",label_visibility = 'collapsed')
    gene_list = []
    
    if len(genes) > 0:
        genes = genes.replace("'","")
        genes = genes.replace('"',"")
        gene_list = genes.split(' ') # First split by spaces
        gene_list = list(filter(lambda a: a != '', gene_list)) # Remove only empty strings
        if ',' in genes:
            gene_list = sum([x.split(',') for x in gene_list],[]) # Flatten with sum(x, [])
        if '\t' in genes:
            gene_list = sum([x.split('\t') for x in gene_list],[])
        if '\n' in genes:
            gene_list = sum([x.split('\n') for x in gene_list],[])
    
    st.write("Gene symbols should be separated by new line.")
    
    if uploaded_file is not None:
        df = pd.read_csv(uploaded_file, header= None)
        data = df.iloc[:,0].tolist()
        uploaded = True
        is_matrix = False
    
    elif len(genes) >0:
        df = pd.DataFrame(gene_list)
        data = gene_list
        uploaded = False
        is_matrix = False

else:  # Expression matrix
    st.markdown("##### Upload expression matrix (CSV/TSV format):")
    st.write("Expected format: 1st column = gene names (from row 2), other columns = sample data")
    uploaded_file = st.file_uploader("Choose an expression matrix file", type=['csv','tsv'])
    
    if uploaded_file is not None:
        # Load with header
        if uploaded_file.name.endswith('.tsv'):
            df = pd.read_csv(uploaded_file, sep='\t')
        else:
            df = pd.read_csv(uploaded_file)
        
        # Extract gene names from 1st column (from row 2)
        data = df.iloc[:, 0].tolist()
        uploaded = True
        is_matrix = True
        
        st.write("Expression matrix preview:")
        st.write(df.head())



if st.button('Run conversion'):

    if file_type == 'Gene list':
        st.write(data[:3])
        df.columns = [species]
        df = df.set_index(species)
    else:
        st.write(f"Converting {len(data)} genes from expression matrix...")

    if species == 'Mouse':
        to_species = 'Human'
    else:
        to_species = 'Mouse'

    if method == 'in-house':

        if species == 'Mouse':
            m2h_dic = load_m2h()
        else:
            h2m_dic = load_h2m()

        df[to_species] = ''

        converted_list=[]


        for i in range(len(data)):
            try:
                if species == "Mouse":
                    res = mouse_human_conversion(m2h_dic, data[i])
                else:
                    res = mouse_human_conversion(h2m_dic, data[i])
                if res != None and res != set():
                    converted_list.append(res)
                    df.loc[data[i],to_species] = ', '.join(res)
                print(res)
            except:
                pass
        flat_list = list(itertools.chain.from_iterable(converted_list))
        d = "\n".join(flat_list)

    else:
        geneinfo_human = pd.read_csv("db/nichenetr.db/geneinfo_human.tsv", sep = '\t')
        geneinfo_2022 = pd.read_csv("db/nichenetr.db/geneinfo_2022.tsv", sep = '\t')


        if species == 'Mouse':
            if method == 'Nicehnetr v1':
                converted_genes = convert_mouse_to_human_symbols(data, version=1)
            elif method == 'Nicehnetr v1 (corrected)':
                # Use corrected version of NichenetR v1
                db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "nichenetr.db", "nichenetr_geneinfo_human_corrected.csv")
                geneinfo_corrected = pd.read_csv(db_path)
                unambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count < 2')['symbol_mouse']
                    .tolist()
                )
                ambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count >= 2')['symbol_mouse']
                    .tolist()
                )
                geneinfo_ambiguous_solved = geneinfo_corrected[
                    (geneinfo_corrected['symbol_mouse'].isin(ambiguous_mouse_genes)) &
                    (geneinfo_corrected['symbol'] == geneinfo_corrected['symbol_mouse'].str.upper())
                ]
                geneinfo_processed = pd.concat([
                    geneinfo_corrected[geneinfo_corrected['symbol_mouse'].isin(unambiguous_mouse_genes)],
                    geneinfo_ambiguous_solved
                ]).dropna()
                mousesymbol2humansymbol = dict(zip(geneinfo_processed['symbol_mouse'], geneinfo_processed['symbol']))
                converted_genes = [mousesymbol2humansymbol.get(symbol, np.nan) for symbol in data]
            elif method == 'Consensus':
                # Use consensus table
                data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "consensus_orthologs_one_to_one.csv")
                consensus_df = pd.read_csv(data_path)
                mapping_dict = dict(zip(consensus_df['mouse_symbol'], consensus_df['human_symbol']))
                converted_genes = [mapping_dict.get(symbol, np.nan) for symbol in data]
            elif method == 'Consensus (corrected)':
                # Use corrected version of consensus table
                data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "consensus_orthologs_one_to_one_corrected.csv")
                consensus_df = pd.read_csv(data_path)
                mapping_dict = dict(zip(consensus_df['mouse_symbol'], consensus_df['human_symbol']))
                converted_genes = [mapping_dict.get(symbol, np.nan) for symbol in data]
            elif method == 'Nichenetr v2 (corrected)':
                # Use corrected version of NichenetR v2
                db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "nichenetr.db", "nichenetr_geneinfo_2022_corrected.csv")
                geneinfo_corrected = pd.read_csv(db_path)
                unambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count < 2')['symbol_mouse']
                    .tolist()
                )
                ambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count >= 2')['symbol_mouse']
                    .tolist()
                )
                geneinfo_ambiguous_solved = geneinfo_corrected[
                    (geneinfo_corrected['symbol_mouse'].isin(ambiguous_mouse_genes)) &
                    (geneinfo_corrected['symbol'] == geneinfo_corrected['symbol_mouse'].str.upper())
                ]
                geneinfo_processed = pd.concat([
                    geneinfo_corrected[geneinfo_corrected['symbol_mouse'].isin(unambiguous_mouse_genes)],
                    geneinfo_ambiguous_solved
                ]).dropna()
                mousesymbol2humansymbol = dict(zip(geneinfo_processed['symbol_mouse'], geneinfo_processed['symbol']))
                converted_genes = [mousesymbol2humansymbol.get(symbol, np.nan) for symbol in data]
            else:
                converted_genes = convert_mouse_to_human_symbols(data, version=2)

        else:
            if method == 'Nicehnetr v1':
                converted_genes = convert_human_to_mouse_symbols(data, version=1)
            elif method == 'Nicehnetr v1 (corrected)':
                # Use corrected version of NichenetR v1
                db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "nichenetr.db", "nichenetr_geneinfo_human_corrected.csv")
                geneinfo_corrected = pd.read_csv(db_path)
                unambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count < 2')['symbol_mouse']
                    .tolist()
                )
                ambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count >= 2')['symbol_mouse']
                    .tolist()
                )
                geneinfo_ambiguous_solved = geneinfo_corrected[
                    (geneinfo_corrected['symbol_mouse'].isin(ambiguous_mouse_genes)) &
                    (geneinfo_corrected['symbol'] == geneinfo_corrected['symbol_mouse'].str.upper())
                ]
                geneinfo_processed = pd.concat([
                    geneinfo_corrected[geneinfo_corrected['symbol_mouse'].isin(unambiguous_mouse_genes)],
                    geneinfo_ambiguous_solved
                ]).dropna()
                humansymbol2mousesymbol = dict(zip(geneinfo_processed['symbol'], geneinfo_processed['symbol_mouse']))
                converted_genes = [humansymbol2mousesymbol.get(symbol, np.nan) for symbol in data]
            elif method == 'Consensus':
                # Use consensus table
                data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "consensus_orthologs_one_to_one.csv")
                consensus_df = pd.read_csv(data_path)
                mapping_dict = dict(zip(consensus_df['human_symbol'], consensus_df['mouse_symbol']))
                converted_genes = [mapping_dict.get(symbol, np.nan) for symbol in data]
            elif method == 'Consensus (corrected)':
                # Use corrected version of consensus table
                data_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data", "consensus_orthologs_one_to_one_corrected.csv")
                consensus_df = pd.read_csv(data_path)
                mapping_dict = dict(zip(consensus_df['human_symbol'], consensus_df['mouse_symbol']))
                converted_genes = [mapping_dict.get(symbol, np.nan) for symbol in data]
            elif method == 'Nichenetr v2 (corrected)':
                # Use corrected version of NichenetR v2
                db_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), "db", "nichenetr.db", "nichenetr_geneinfo_2022_corrected.csv")
                geneinfo_corrected = pd.read_csv(db_path)
                unambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count < 2')['symbol_mouse']
                    .tolist()
                )
                ambiguous_mouse_genes = (
                    geneinfo_corrected.dropna()
                    .groupby('symbol_mouse').size()
                    .reset_index(name='count')
                    .query('count >= 2')['symbol_mouse']
                    .tolist()
                )
                geneinfo_ambiguous_solved = geneinfo_corrected[
                    (geneinfo_corrected['symbol_mouse'].isin(ambiguous_mouse_genes)) &
                    (geneinfo_corrected['symbol'] == geneinfo_corrected['symbol_mouse'].str.upper())
                ]
                geneinfo_processed = pd.concat([
                    geneinfo_corrected[geneinfo_corrected['symbol_mouse'].isin(unambiguous_mouse_genes)],
                    geneinfo_ambiguous_solved
                ]).dropna()
                humansymbol2mousesymbol = dict(zip(geneinfo_processed['symbol'], geneinfo_processed['symbol_mouse']))
                converted_genes = [humansymbol2mousesymbol.get(symbol, np.nan) for symbol in data]
            else:
                converted_genes = convert_human_to_mouse_symbols(data, version=2)

        if file_type == 'Gene list':
            #converted_list = [x for x in converted_genes if not pd.isna(x)]
            df[to_species] = converted_genes
            converted_genes = [x for x in converted_genes if pd.isnull(x) == False] # Remove NaN
            converted_genes = sorted(set(converted_genes), key=converted_genes.index) # Remove duplicates
            d = "\n".join(converted_genes) # Remove None and convert to string
        else:
            # For expression matrix: change gene names in 1st column
            df_converted = df.copy()
            df_converted.iloc[:, 0] = converted_genes
            
            # Conversion statistics
            total_genes = len(converted_genes)
            successful_conversions = sum(1 for x in converted_genes if pd.notna(x))
            
            st.write(f"Conversion completed: {successful_conversions}/{total_genes} genes converted")
            
            # Identify genes that could not be converted
            failed_genes = [data[i] for i, x in enumerate(converted_genes) if pd.isna(x)]
            if len(failed_genes) > 0:
                st.write(f"**Genes removed (no conversion found): {len(failed_genes)}**")
                if len(failed_genes) <= 20:
                    st.write(", ".join(failed_genes))
                else:
                    st.write(", ".join(failed_genes[:20]) + f"... (and {len(failed_genes)-20} more)")
            
            # Delete rows that could not be converted (rows with None)
            mask = pd.notna(converted_genes)
            df_converted = df_converted[mask].reset_index(drop=True)
            converted_genes_clean = [x for x in converted_genes if pd.notna(x)]
            
            # Handle duplicate gene names: calculate mean for many-to-one mappings
            if len(converted_genes_clean) > 0:
                df_converted.iloc[:, 0] = converted_genes_clean
                
                # If duplicate gene names exist, calculate mean of numeric columns
                numeric_cols = df_converted.select_dtypes(include=[np.number]).columns
                if len(numeric_cols) > 0:
                    # Group by 1st column (gene name) and calculate mean of numeric columns
                    gene_col = df_converted.columns[0]
                    agg_dict = {col: 'mean' for col in numeric_cols}
                    # Get first value for non-numeric columns
                    non_numeric_cols = [col for col in df_converted.columns if col not in numeric_cols and col != gene_col]
                    for col in non_numeric_cols:
                        agg_dict[col] = 'first'
                    
                    # Identify duplicate genes (before merge)
                    duplicate_genes = df_converted[df_converted.duplicated(gene_col, keep=False)][gene_col].unique()
                    
                    df_converted = df_converted.groupby(gene_col, as_index=False).agg(agg_dict)
                    
                    duplicates_removed = len(converted_genes_clean) - len(df_converted)
                    if duplicates_removed > 0:
                        st.write(f"**Merged {duplicates_removed} duplicate genes by averaging expression values:**")
                        if len(duplicate_genes) <= 50:
                            st.write(", ".join(duplicate_genes))
                        else:
                            st.write(", ".join(duplicate_genes[:50]) + f"... (and {len(duplicate_genes)-50} more)")
                
                # Final gene list
                final_genes = df_converted.iloc[:, 0].tolist()
                d = "\n".join(final_genes)
            else:
                d = ""

    if file_type == 'Gene list':
        st.write(df)
    else:
        st.write("Converted expression matrix preview:")
        st.write(df_converted.head())


    if uploaded:
        base_name = os.path.splitext(uploaded_file.name)[0]
        if file_type == 'Expression matrix':
            OutDataName = base_name + '.' + to_species + '_genes.txt'
            OutDfName = base_name + '.' + species + '-' + to_species + '_matrix.tsv'
        else:
            OutDataName = base_name + '.' + to_species + '.txt'
            OutDfName = base_name + '.' + species + '-' + to_species + '.tsv'
        out_zip_name = base_name + ".zip"
    else:
        OutDataName = to_species + '.txt'
        OutDfName = to_species + '.tsv'
        out_zip_name = to_species + ".zip"

#    st.write(d)
    with open('res/' + OutDataName, mode = 'w') as f:
        f.write(d)

    if file_type == 'Gene list':
        df.to_csv('res/' + OutDfName, sep = '\t')
    else:
        df_converted.to_csv('res/' + OutDfName, sep='\t', index=False)


    shutil.make_archive("res", format='zip',root_dir= 'res')


    with open("res.zip", "rb") as fp:
        btn = st.download_button(
            label="Download Results",
        data=fp,
        file_name=out_zip_name,
        mime = "zip"
        )



