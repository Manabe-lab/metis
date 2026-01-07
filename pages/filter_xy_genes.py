import streamlit as st
import pandas as pd
import pickle
import sys
import os
import re

st.set_page_config(page_title="Filter Special Gene Categories", page_icon="🧬")

st.title("🧬 Filter Special Gene Categories")
st.markdown("""
Filter or extract genes from special categories:
- **Sex chromosomes**: X and Y chromosome genes
- **Mitochondrial genes**: chrM/chrMT genes
- **Ribosomal genes**: Ribosomal protein genes (Rpl*, Rps*, RPL*, RPS*)
- **Pseudogenes**: Annotated pseudogenes (gene_type contains 'pseudogene')
- **Retrogenes**: Retrotransposed genes (tag contains 'retrogene')

**Supported genomes:**
- **Mouse**: mm10 (10x), mm39 (10x), mm10 (HOMER), mm39 (HOMER)
- **Human**: hg38 (10x 2020A), hg38 (10x 2024A), hg38 (HOMER)
""")

@st.cache_data
def load_xy_gene_dict():
    """Load all X/Y chromosome gene reference dictionaries"""
    # Use relative paths from streamlit root (where metis.py is located)
    possible_data_dirs = [
        './data',
        '../data',
        'data',
        os.path.join(os.path.dirname(__file__), '../data'),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')
    ]

    data_dir = None
    for dir_path in possible_data_dirs:
        if os.path.exists(dir_path):
            data_dir = dir_path
            break

    if data_dir is None:
        st.error(f"❌ Data directory not found. Please ensure X/Y gene reference files are in the data directory.")
        return {}

    xy_genes = {}
    missing_files = []

    # Load 10x references
    for genome in ['mm10_10x', 'mm39_10x', 'hg38_2020_10x', 'hg38_2024_10x']:
        pkl_file = os.path.join(data_dir, f'{genome}_xy_genes.pkl')
        if os.path.exists(pkl_file):
            try:
                with open(pkl_file, 'rb') as f:
                    xy_genes[genome] = pickle.load(f)
            except Exception as e:
                st.warning(f"Failed to load {pkl_file}: {e}")
        else:
            missing_files.append(f'{genome}_xy_genes.pkl')

    # Load HOMER references
    for genome in ['mm10_homer', 'mm39_homer', 'hg38_homer']:
        pkl_file = os.path.join(data_dir, f'{genome}_xy_genes.pkl')
        if os.path.exists(pkl_file):
            try:
                with open(pkl_file, 'rb') as f:
                    xy_genes[genome] = pickle.load(f)
            except Exception as e:
                st.warning(f"Failed to load {pkl_file}: {e}")
        else:
            missing_files.append(f'{genome}_xy_genes.pkl')

    # Show informative message about missing files
    if missing_files:
        with st.expander(f"⚠️ {len(missing_files)} reference file(s) not found (click to see details)"):
            st.write(f"**Data directory**: `{data_dir}`")
            st.write("**Missing files**:")
            for f in missing_files:
                st.write(f"- {f}")
            st.info("These references will not be available. To add them, place the pkl files in the data directory.")

    return xy_genes

@st.cache_data
def load_special_gene_dict():
    """Load all special gene reference dictionaries (mitochondrial, ribosomal, pseudogene)"""
    possible_data_dirs = [
        './data',
        '../data',
        'data',
        os.path.join(os.path.dirname(__file__), '../data'),
        os.path.join(os.path.dirname(os.path.abspath(__file__)), '../data')
    ]

    data_dir = None
    for dir_path in possible_data_dirs:
        if os.path.exists(dir_path):
            data_dir = dir_path
            break

    if data_dir is None:
        return {}

    special_genes = {}
    missing_files = []

    # Load 10x references only (HOMER not yet supported for special genes)
    for genome in ['mm10_10x', 'mm39_10x', 'hg38_2020_10x', 'hg38_2024_10x']:
        pkl_file = os.path.join(data_dir, f'{genome}_special_genes.pkl')
        if os.path.exists(pkl_file):
            try:
                with open(pkl_file, 'rb') as f:
                    special_genes[genome] = pickle.load(f)
            except Exception as e:
                st.warning(f"Failed to load {pkl_file}: {e}")
        else:
            missing_files.append(f'{genome}_special_genes.pkl')

    if missing_files and len(special_genes) == 0:
        with st.expander(f"⚠️ Special gene reference files not found"):
            st.write(f"**Data directory**: `{data_dir}`")
            st.info("Special gene filtering (MT, ribosomal, pseudogene) will not be available. Run `scripts/extract_special_genes.py` to generate these files.")

    return special_genes

@st.cache_data
def convert_df_to_csv(df):
    """Convert dataframe to CSV for download"""
    return df.to_csv(index=False).encode('utf-8')

@st.cache_data
def convert_list_to_txt(gene_list):
    """Convert gene list to text for download"""
    return '\n'.join(gene_list).encode('utf-8')

# Load reference data
xy_gene_dict = load_xy_gene_dict()
special_gene_dict = load_special_gene_dict()

if not xy_gene_dict:
    st.error("❌ **No X/Y chromosome gene reference files found**")
    st.markdown("""
### How to set up reference files:

1. **Generate reference files** by running these scripts on your local machine:
   ```bash
   python3 scripts/extract_xy_genes.py
   python3 scripts/extract_xy_genes_homer.py
   ```

2. **Upload pkl files** to the server's data directory:
   - `mm10_10x_xy_genes.pkl`
   - `mm39_10x_xy_genes.pkl`
   - `hg38_2020_10x_xy_genes.pkl`
   - `hg38_2024_10x_xy_genes.pkl`
   - `mm10_homer_xy_genes.pkl`
   - `mm39_homer_xy_genes.pkl`
   - `hg38_homer_xy_genes.pkl`

3. **Data directory location**: Check the expandable warning above for the expected path.
    """)
    st.stop()

# UI - Genome selection
st.markdown("### 1. Select Reference Genome")

# Step 1: Species selection
species = st.radio(
    "Select species:",
    ["Human", "Mouse"],
    index=1,  # Default to Mouse
    horizontal=True
)

# Step 2: Genome reference selection based on species
genome_options = {
    'Human': {
        '10x Genomics': {
            'hg38_2020_10x': 'hg38 (2020-A reference)',
            'hg38_2024_10x': 'hg38 (2024-A reference)'
        },
        'HOMER': {
            'hg38_homer': 'hg38/hg19'
        }
    },
    'Mouse': {
        '10x Genomics': {
            'mm10_10x': 'mm10',
            'mm39_10x': 'mm39'
        },
        'HOMER': {
            'mm10_homer': 'mm10',
            'mm39_homer': 'mm39'
        }
    }
}

# Get available options for selected species
species_options = genome_options[species]

# Create columns for annotation source selection
col1, col2 = st.columns(2)

with col1:
    annotation_source = st.radio(
        "Select annotation source:",
        list(species_options.keys()),
        help="10x Genomics: Official Cell Ranger references\nHOMER: HOMER genome annotations"
    )

with col2:
    # Get available genomes for selected source
    available_refs = {k: v for k, v in species_options[annotation_source].items() if k in xy_gene_dict}

    if not available_refs:
        st.error(f"⚠️ No {annotation_source} references available for {species}")
        st.stop()

    genome_choice = st.selectbox(
        f"Select {annotation_source} version:",
        list(available_refs.keys()),
        format_func=lambda x: available_refs[x]
    )

# Display reference info
if genome_choice:
    ref_data = xy_gene_dict[genome_choice]

    # Create display name
    source_display = "10x Genomics" if "10x" in genome_choice else "HOMER"
    version_display = available_refs[genome_choice]

    info_text = f"📊 **{species} - {source_display} - {version_display}**\n"
    info_text += f"- X chromosome: {len(ref_data['X'])} genes\n"
    info_text += f"- Y chromosome: {len(ref_data['Y'])} genes"

    # Add special gene info if available (only for 10x references)
    if genome_choice in special_gene_dict:
        special_ref = special_gene_dict[genome_choice]
        info_text += f"\n- Mitochondrial: {len(special_ref.get('mitochondrial', []))} genes"
        info_text += f"\n- Ribosomal: {len(special_ref.get('ribosomal', []))} genes"
        info_text += f"\n- Pseudogenes: {len(special_ref.get('pseudogene', []))} genes"
        if 'retrogene' in special_ref:
            info_text += f"\n- Retrogenes: {len(special_ref['retrogene'])} genes"
        if 'lncrna' in special_ref:
            info_text += f"\n- lncRNA: {len(special_ref['lncrna'])} genes"
        if 'small_rna' in special_ref:
            info_text += f"\n- Small RNA: {len(special_ref['small_rna'])} genes"

    st.info(info_text)

# UI - Input method
st.markdown("### 2. Input Gene List")

input_method = st.radio(
    "How to provide gene list:",
    ('Gene list file', 'Data table file (filter rows)', 'Manual input')
)

gene_list = []
filtered_df = None  # For data table filtering mode
base_name = None  # For storing original filename

if input_method == 'Gene list file':
    st.markdown("Upload a file containing gene symbols (one per line, or as a column in CSV/TSV/Excel)")

    uploaded_file = st.file_uploader(
        "Choose a file",
        type=['txt', 'csv', 'tsv', 'xls', 'xlsx']
    )

    if uploaded_file is not None:
        file_ext = uploaded_file.name.split('.')[-1].lower()

        if file_ext == 'txt':
            gene_list = uploaded_file.read().decode('utf-8').strip().split('\n')
            gene_list = [g.strip() for g in gene_list if g.strip()]

        elif file_ext == 'csv':
            df = pd.read_csv(uploaded_file)
            st.write("**Preview of uploaded file:**")
            st.dataframe(df.head())

            gene_column = st.selectbox('Select gene column:', df.columns.tolist())
            gene_list = df[gene_column].dropna().astype(str).tolist()

        elif file_ext == 'tsv':
            df = pd.read_csv(uploaded_file, sep='\t')
            st.write("**Preview of uploaded file:**")
            st.dataframe(df.head())

            gene_column = st.selectbox('Select gene column:', df.columns.tolist())
            gene_list = df[gene_column].dropna().astype(str).tolist()

        elif file_ext in ['xls', 'xlsx']:
            df = pd.read_excel(uploaded_file)
            st.write("**Preview of uploaded file:**")
            st.dataframe(df.head())

            gene_column = st.selectbox('Select gene column:', df.columns.tolist())
            gene_list = df[gene_column].dropna().astype(str).tolist()

elif input_method == 'Data table file (filter rows)':
    st.markdown("Upload a data table file (1st column = gene names). Rows will be filtered based on X/Y chromosome genes.")

    # File type selection
    file_type = st.radio(
        "Data format:",
        ('auto', 'tsv', 'csv', 'excel'),
        horizontal=True
    )

    uploaded_file = st.file_uploader(
        "Choose a file",
        type=['txt', 'csv', 'tsv', 'xls', 'xlsx'],
        key='data_table_upload'
    )

    if uploaded_file is not None:
        # Store original filename
        original_filename = uploaded_file.name
        base_name = os.path.splitext(original_filename)[0]

        # Read file with auto detection
        if file_type == 'auto':
            try:
                df = pd.read_csv(uploaded_file, sep=None, engine='python', index_col=0)
            except:  # Excel
                df = pd.read_excel(uploaded_file, index_col=0)
        elif file_type == 'csv':
            df = pd.read_csv(uploaded_file, index_col=0)
        elif file_type == 'tsv':
            df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
        elif file_type == 'excel':
            df = pd.read_excel(uploaded_file, index_col=0)

        st.write("**Preview of uploaded data:**")
        st.dataframe(df.head())

        # Get gene list from index (1st column)
        gene_list = df.index.astype(str).tolist()
        filtered_df = df  # Store for later filtering

        st.info(f"✓ Loaded {len(gene_list)} genes from data table (index column)")

else:  # Manual input
    manual_input = st.text_area(
        'Enter gene symbols (separated by newline, tab, comma, or space):',
        height=200,
        placeholder='Xist\nDdx3y\nKdm5d\n...\n\nOr: Xist, Ddx3y, Kdm5d\nOr: Xist\tDdx3y\tKdm5d',
        help="Supported delimiters: newline (\\n), tab (\\t), comma (,), space"
    )
    if manual_input:
        # Split by multiple delimiters: newline, tab, comma, space
        gene_list = re.split(r'[\n\t,\s]+', manual_input.strip())
        gene_list = [g.strip() for g in gene_list if g.strip()]

        st.info(f"✓ Parsed {len(gene_list)} genes from input")

# UI - Filter operation
st.markdown("### 3. Select Gene Categories to Filter")

# Gene category selection
gene_categories = st.multiselect(
    'Select gene categories:',
    ['X chromosome', 'Y chromosome', 'Mitochondrial', 'Ribosomal', 'Pseudogene', 'Retrogene', 'lncRNA', 'Small RNA', 'Custom gene list'],
    default=['X chromosome', 'Y chromosome'],
    help="Select which gene categories to filter. Only categories available for your selected genome will be active. 'Custom gene list' allows you to input your own gene names to filter."
)

# Check if special genes are available for selected genome
special_available = genome_choice in special_gene_dict if 'genome_choice' in locals() else False

# Disable special categories if not available
if not special_available and any(cat in gene_categories for cat in ['Mitochondrial', 'Ribosomal', 'Pseudogene', 'Retrogene', 'lncRNA', 'Small RNA']):
    st.warning("⚠️ Special gene categories (Mitochondrial, Ribosomal, Pseudogene, Retrogene, lncRNA, Small RNA) are only available for 10x Genomics references. Please select a 10x genome or run `scripts/extract_special_genes.py` to generate reference files.")

# Custom gene list input
custom_filter_genes = []
if 'Custom gene list' in gene_categories:
    st.markdown("#### Custom Gene List for Filtering")
    st.markdown("Enter gene names to filter (matching is **case-insensitive**):")
    custom_genes_input = st.text_area(
        'Custom gene names:',
        height=150,
        placeholder='Enter gene names separated by:\n- Newline\n- Comma (,)\n- Semicolon (;)\n- Tab\n- Space\n\nExample: Gapdh, Actb; Rpl13a\nOr one per line:\nGapdh\nActb\nRpl13a',
        help="Supported delimiters: newline, comma, semicolon, tab, space. Matching is case-insensitive.",
        key='custom_filter_genes_input'
    )
    if custom_genes_input:
        # Parse with multiple delimiters: newline, comma, semicolon, tab, space
        custom_filter_genes = re.split(r'[\n\t,;\s]+', custom_genes_input.strip())
        custom_filter_genes = [g.strip() for g in custom_filter_genes if g.strip()]
        st.info(f"✓ Parsed {len(custom_filter_genes)} custom genes for filtering")

st.markdown("### 4. Select Operation")

operation = st.radio(
    "What would you like to do?",
    ('Remove selected genes (exclude)', 'Extract selected genes only', 'Separate by category')
)

# Process genes
if gene_list and st.button('🔬 Process Genes'):

    st.markdown("---")
    st.markdown("### Results")

    # Get gene reference data
    ref_data = xy_gene_dict[genome_choice]

    # Create case-insensitive lookup dictionaries for all categories
    category_lookups = {}

    # X/Y chromosomes
    if 'X chromosome' in gene_categories:
        category_lookups['X chromosome'] = {g.lower(): g for g in ref_data['X']}
    if 'Y chromosome' in gene_categories:
        category_lookups['Y chromosome'] = {g.lower(): g for g in ref_data['Y']}

    # Special genes (if available)
    if genome_choice in special_gene_dict:
        special_ref = special_gene_dict[genome_choice]
        if 'Mitochondrial' in gene_categories and 'mitochondrial' in special_ref:
            category_lookups['Mitochondrial'] = {g.lower(): g for g in special_ref['mitochondrial']}
        if 'Ribosomal' in gene_categories and 'ribosomal' in special_ref:
            category_lookups['Ribosomal'] = {g.lower(): g for g in special_ref['ribosomal']}
        if 'Pseudogene' in gene_categories and 'pseudogene' in special_ref:
            category_lookups['Pseudogene'] = {g.lower(): g for g in special_ref['pseudogene']}
        if 'Retrogene' in gene_categories and 'retrogene' in special_ref:
            category_lookups['Retrogene'] = {g.lower(): g for g in special_ref['retrogene']}
        if 'lncRNA' in gene_categories and 'lncrna' in special_ref:
            category_lookups['lncRNA'] = {g.lower(): g for g in special_ref['lncrna']}
        if 'Small RNA' in gene_categories and 'small_rna' in special_ref:
            category_lookups['Small RNA'] = {g.lower(): g for g in special_ref['small_rna']}

    # Custom gene list (case-insensitive matching)
    if 'Custom gene list' in gene_categories and custom_filter_genes:
        category_lookups['Custom gene list'] = {g.lower(): g for g in custom_filter_genes}

    # Classify input genes (case-insensitive matching)
    input_genes_set = set(gene_list)

    # Find matches for each category
    category_genes = {cat: set() for cat in gene_categories}

    for gene in input_genes_set:
        gene_lower = gene.lower()
        for category, lookup in category_lookups.items():
            if gene_lower in lookup:
                category_genes[category].add(gene)

    # Calculate genes found in any selected category
    special_found = set()
    for cat_set in category_genes.values():
        special_found.update(cat_set)

    # Remaining genes (not in any selected category)
    other_genes = input_genes_set - special_found

    # Legacy variables for backward compatibility
    x_found = category_genes.get('X chromosome', set())
    y_found = category_genes.get('Y chromosome', set())
    xy_found = x_found | y_found
    autosomal = other_genes

    # Display results based on operation
    if operation == 'Remove selected genes (exclude)':
        st.success(f"✅ Removed {len(special_found)} genes from selected categories")
        st.info(f"📋 Remaining genes: {len(other_genes)}")

        if len(special_found) > 0:
            st.markdown("**Removed genes:**")
            # Determine category for each removed gene
            gene_categories_list = []
            for g in sorted(list(special_found)):
                cats = [cat for cat, genes in category_genes.items() if g in genes]
                gene_categories_list.append(', '.join(cats) if cats else 'Unknown')

            removed_df = pd.DataFrame({
                'Gene': sorted(list(special_found)),
                'Category': gene_categories_list
            })
            st.dataframe(removed_df)

        # For data table mode, filter rows
        if input_method == 'Data table file (filter rows)' and filtered_df is not None:
            # Filter DataFrame to keep only other genes (not in selected categories)
            result_table = filtered_df[filtered_df.index.isin(other_genes)]
            st.markdown("**Filtered data table (remaining genes):**")
            st.dataframe(result_table.head(20))
            st.info(f"Showing first 20 rows of {len(result_table)} total rows")

            # Download filtered data table
            csv = result_table.to_csv(index=True, sep='\t').encode('utf-8')
            download_filename = f"{base_name}_filtered_remaining.tsv" if base_name else "filtered_remaining_data.tsv"
            st.download_button(
                label="📥 Download filtered data table (TSV)",
                data=csv,
                file_name=download_filename,
                mime='text/tab-separated-values'
            )
        else:
            # For gene list mode
            result_df = pd.DataFrame({'Gene': sorted(list(other_genes))})
            st.markdown("**Filtered gene list (remaining genes):**")
            st.dataframe(result_df)

            # Download button
            tsv = result_df.to_csv(index=False, sep='\t').encode('utf-8')
            st.download_button(
                label="📥 Download filtered genes (TSV)",
                data=tsv,
                file_name='filtered_remaining_genes.tsv',
                mime='text/tab-separated-values'
            )

    elif operation == 'Extract selected genes only':
        st.success(f"✅ Found {len(special_found)} genes in selected categories")

        # For data table mode, filter rows
        if input_method == 'Data table file (filter rows)' and filtered_df is not None:
            # Filter DataFrame to keep only genes in selected categories
            result_table = filtered_df[filtered_df.index.isin(special_found)]
            st.markdown("**Filtered data table (selected categories only):**")
            st.dataframe(result_table.head(20))
            st.info(f"Showing first 20 rows of {len(result_table)} total rows")

            # Download filtered data table
            csv = result_table.to_csv(index=True, sep='\t').encode('utf-8')
            download_filename = f"{base_name}_filtered_selected.tsv" if base_name else "selected_categories_data.tsv"
            st.download_button(
                label="📥 Download selected genes data table (TSV)",
                data=csv,
                file_name=download_filename,
                mime='text/tab-separated-values'
            )
        else:
            # For gene list mode
            # Determine category for each gene
            gene_categories_list = []
            for g in sorted(list(special_found)):
                cats = [cat for cat, genes in category_genes.items() if g in genes]
                gene_categories_list.append(', '.join(cats) if cats else 'Unknown')

            result_df = pd.DataFrame({
                'Gene': sorted(list(special_found)),
                'Category': gene_categories_list
            })
            st.dataframe(result_df)

            # Download button
            tsv = result_df.to_csv(index=False, sep='\t').encode('utf-8')
            st.download_button(
                label="📥 Download selected genes (TSV)",
                data=tsv,
                file_name='selected_category_genes.tsv',
                mime='text/tab-separated-values'
            )

    else:  # Separate by category
        # Create metrics for all selected categories
        cols = st.columns(min(len(gene_categories) + 1, 5))  # Max 5 columns

        for idx, category in enumerate(gene_categories):
            with cols[idx % len(cols)]:
                st.metric(category, len(category_genes[category]))

        # Add "Other" category
        if len(gene_categories) < len(cols):
            with cols[len(gene_categories)]:
                st.metric("Other genes", len(other_genes))

        # Display each category in tabs
        tab_names = list(gene_categories) + ["Other genes"]
        tabs = st.tabs(tab_names)

        # Dynamically create tabs for each category
        for idx, tab in enumerate(tabs):
            with tab:
                if idx < len(gene_categories):
                    # Category tabs
                    category = gene_categories[idx]
                    category_set = category_genes[category]

                    if len(category_set) > 0:
                        if input_method == 'Data table file (filter rows)' and filtered_df is not None:
                            cat_table = filtered_df[filtered_df.index.isin(category_set)]
                            st.dataframe(cat_table.head(20))
                            st.info(f"Showing first 20 rows of {len(cat_table)} total rows")

                            tsv = cat_table.to_csv(index=True, sep='\t').encode('utf-8')
                            safe_cat_name = category.replace(' ', '_').lower()
                            download_filename = f"{base_name}_filtered_{safe_cat_name}.tsv" if base_name else f"{safe_cat_name}_data.tsv"
                            st.download_button(
                                label=f"📥 Download {category} data (TSV)",
                                data=tsv,
                                file_name=download_filename,
                                mime='text/tab-separated-values',
                                key=f'download_{safe_cat_name}'
                            )
                        else:
                            cat_df = pd.DataFrame({'Gene': sorted(list(category_set))})
                            st.dataframe(cat_df)

                            tsv = cat_df.to_csv(index=False, sep='\t').encode('utf-8')
                            safe_cat_name = category.replace(' ', '_').lower()
                            st.download_button(
                                label=f"📥 Download {category} genes (TSV)",
                                data=tsv,
                                file_name=f'{safe_cat_name}_genes.tsv',
                                mime='text/tab-separated-values',
                                key=f'download_{safe_cat_name}'
                            )
                    else:
                        st.info(f"No {category} genes found")
                else:
                    # "Other genes" tab
                    if len(other_genes) > 0:
                        if input_method == 'Data table file (filter rows)' and filtered_df is not None:
                            other_table = filtered_df[filtered_df.index.isin(other_genes)]
                            st.dataframe(other_table.head(20))
                            st.info(f"Showing first 20 rows of {len(other_table)} total rows")

                            tsv = other_table.to_csv(index=True, sep='\t').encode('utf-8')
                            download_filename = f"{base_name}_filtered_other.tsv" if base_name else "other_genes_data.tsv"
                            st.download_button(
                                label="📥 Download other genes data (TSV)",
                                data=tsv,
                                file_name=download_filename,
                                mime='text/tab-separated-values',
                                key='download_other'
                            )
                        else:
                            other_df = pd.DataFrame({'Gene': sorted(list(other_genes))})
                            st.dataframe(other_df)

                            tsv = other_df.to_csv(index=False, sep='\t').encode('utf-8')
                            st.download_button(
                                label="📥 Download other genes (TSV)",
                                data=tsv,
                                file_name='other_genes.tsv',
                                mime='text/tab-separated-values',
                                key='download_other'
                            )
                    else:
                        st.info("No other genes found")

    # Summary statistics
    st.markdown("---")
    st.markdown("### Summary Statistics")

    summary_data = {'Category': ['Total input genes'], 'Count': [len(input_genes_set)], 'Percentage': [100.0]}

    # Add each selected category
    for category in gene_categories:
        count = len(category_genes[category])
        pct = count/len(input_genes_set)*100 if len(input_genes_set) > 0 else 0
        summary_data['Category'].append(category)
        summary_data['Count'].append(count)
        summary_data['Percentage'].append(pct)

    # Add total selected and other genes
    summary_data['Category'].extend(['Total selected', 'Other genes'])
    summary_data['Count'].extend([len(special_found), len(other_genes)])
    summary_data['Percentage'].extend([
        len(special_found)/len(input_genes_set)*100 if len(input_genes_set) > 0 else 0,
        len(other_genes)/len(input_genes_set)*100 if len(input_genes_set) > 0 else 0
    ])

    summary_df = pd.DataFrame(summary_data)
    summary_df['Percentage'] = summary_df['Percentage'].apply(lambda x: f'{x:.2f}%')
    st.dataframe(summary_df)

# Footer
st.markdown("---")
st.markdown("""
**Note:** Gene symbol matching is **case-insensitive** (e.g., "eif2s3y", "Eif2s3y", "EIF2S3Y" are all treated as the same gene).

**Custom gene list:** You can input your own gene names to filter. Supported delimiters include:
- Newline
- Comma (,)
- Semicolon (;)
- Tab
- Space

**Reference sources:**
- 10x Genomics: Official Cell Ranger reference packages
- HOMER: HOMER genome annotation files
""")
