#!/usr/bin/env python3
"""
Microarray Gene Name Filter
A tool for extracting, filtering, and aggregating gene names from microarray data
"""

import streamlit as st
import pandas as pd
import numpy as np
import io
import re
from collections import defaultdict, Counter

class GeneInfoExtractor:
    """Class for learning position patterns and extracting gene information"""

    def __init__(self):
        self.position_patterns = {}
        self.separator_pattern = None
        self.learned = False

    def detect_separator(self, sample_texts):
        """Automatically detect the delimiter used in sample texts"""
        separators = {
            ' // ': r'\s*//\s*',
            '//': r'//',
            ' / ': r'\s*/\s*',
            '/': r'/',
            '\t': r'\t',
            ', ': r',\s*',
            ',': r',',
            ' ': r'\s+',
            '|': r'\|',
            ';': r';'
        }

        # Count occurrences of each delimiter
        counts = {}
        for name, pattern in separators.items():
            total_count = 0
            for text in sample_texts:
                if pd.notna(text) and text != "---":
                    count = len(re.findall(pattern, str(text)))
                    total_count += count
            if total_count > 0:
                counts[name] = total_count

        # Return the most frequent delimiter
        if counts:
            most_common = max(counts, key=counts.get)
            return most_common, separators[most_common]
        else:
            return ' // ', r'\s*//\s*'  # Default

    def classify_element(self, element):
        """Classify the type of element"""
        element = str(element).strip()

        if not element or element == "---":
            return "empty"

        # RefSeq ID
        if any(element.startswith(prefix) for prefix in ["NM_", "NR_", "XM_", "XR_", "NP_", "XP_"]):
            return "refseq"

        # Ensembl ID
        if element.startswith("ENSMUST") or element.startswith("ENSMUSG") or element.startswith("ENSMUS") or \
           element.startswith("ENST") or element.startswith("ENSG") or element.startswith("ENS"):
            return "ensembl"

        # GenBank clone ID
        if element.startswith("BC") and len(element) > 2 and element[2:].replace(".", "").isdigit():
            return "genbank"

        if element.startswith("AK") and len(element) > 2 and element[2:].isdigit():
            return "genbank"

        # Entrez Gene ID (numbers only)
        if element.isdigit():
            return "entrez"

        # Chromosome location information
        if "|" in element or "cM" in element or element.startswith("chr"):
            return "location"

        # Description detection
        words = element.split()
        if len(words) >= 3:
            # Long description
            if not re.match(r'^[A-Za-z0-9]+[-]?[A-Za-z0-9]+$', element):
                return "description"

        # Description containing specific keywords
        skip_keywords = ["uncharacterized", "predicted", "hypothetical", "family", "domain",
                        "containing", "protein", "member", "homolog", "like", "antigen", "receptor"]
        if any(keyword in element.lower() for keyword in skip_keywords):
            if len(element.split()) > 2 or len(element) > 30:
                return "description"

        # Likely a gene symbol
        if re.match(r'^([A-Za-z][A-Za-z0-9_-]*\d*[A-Za-z0-9]*|[0-9]+[A-Za-z][A-Za-z0-9]*[A-Za-z]+)$', element):
            return "gene_symbol"

        # Other (probably description)
        return "description"

    def learn_pattern(self, sample_data, separator=None):
        """Learn position patterns from sample data"""
        # Detect delimiter
        if separator is None:
            sample_texts = [str(text) for text in sample_data if pd.notna(text) and text != "---"]
            sep_name, self.separator_pattern = self.detect_separator(sample_texts[:20])
            st.info(f"🔍 Detected delimiter: '{sep_name}'")
        else:
            self.separator_pattern = separator

        # Count element type occurrence frequency at each position
        position_counts = defaultdict(Counter)
        max_positions = 0

        for text in sample_data:
            if pd.isna(text) or text == "---" or text == "":
                continue

            text = str(text)
            # For multiple entries (/// delimiter), use only the first entry
            if " /// " in text:
                text = text.split(" /// ")[0]

            # Split by delimiter
            parts = re.split(self.separator_pattern, text)
            max_positions = max(max_positions, len(parts))

            for pos, part in enumerate(parts):
                element_type = self.classify_element(part)
                if element_type != "empty":
                    position_counts[pos][element_type] += 1

        # Determine the most frequent element type at each position
        self.position_patterns = {}

        st.write("📊 Position pattern analysis results:")
        pattern_summary = []

        for pos in range(max_positions):
            if pos in position_counts:
                # Get the most frequent type
                most_common = position_counts[pos].most_common()
                if most_common:
                    top_type = most_common[0][0]
                    top_count = most_common[0][1]
                    total = sum(position_counts[pos].values())
                    confidence = top_count / total * 100

                    # Only adopt as pattern if confidence is 40% or higher
                    if confidence >= 40:
                        self.position_patterns[pos] = top_type
                        pattern_summary.append({
                            "Position": pos + 1,
                            "Estimated Type": self.get_type_label(top_type),
                            "Confidence": f"{confidence:.1f}%",
                            "Sample Count": total
                        })

        if pattern_summary:
            pattern_df = pd.DataFrame(pattern_summary)
            st.dataframe(pattern_df, use_container_width=True)

        self.learned = True
        return self.position_patterns

    def get_type_label(self, type_key):
        """Convert type key to label"""
        labels = {
            "gene_symbol": "Gene Symbol",
            "refseq": "RefSeq ID",
            "ensembl": "Ensembl ID",
            "genbank": "GenBank ID",
            "entrez": "Entrez ID",
            "description": "Description",
            "location": "Location",
            "empty": "Empty"
        }
        return labels.get(type_key, type_key)

    def extract_with_pattern(self, text):
        """Extract information using learned patterns"""
        result = {
            'gene_symbol': '',
            'refseq_id': '',
            'ensembl_id': '',
            'genbank_id': '',
            'entrez_id': '',
            'description': '',
            'location': ''
        }

        if pd.isna(text) or text == "---" or text == "":
            return result

        text = str(text)

        # For multiple entries (/// delimiter), use only the first entry
        if " /// " in text:
            text = text.split(" /// ")[0]

        # Split by delimiter
        parts = re.split(self.separator_pattern, text)

        # Collect values for each type
        symbols = []
        refseqs = []
        ensembls = []
        genbanks = []
        entrezs = []
        descriptions = []
        locations = []

        for pos, part in enumerate(parts):
            part = part.strip()
            if not part or part == "---":
                continue

            # Prioritize if there is a learned pattern
            if pos in self.position_patterns:
                expected_type = self.position_patterns[pos]
                actual_type = self.classify_element(part)

                # When expected type matches actual type
                if expected_type == actual_type:
                    if expected_type == "gene_symbol" and not symbols:
                        symbols.append(part)
                    elif expected_type == "refseq":
                        refseqs.append(part)
                    elif expected_type == "ensembl":
                        ensembls.append(part)
                    elif expected_type == "genbank":
                        genbanks.append(part)
                    elif expected_type == "entrez":
                        entrezs.append(part)
                    elif expected_type == "description":
                        descriptions.append(part)
                    elif expected_type == "location":
                        locations.append(part)
                # Even if types don't match, save based on actual type
                else:
                    if actual_type == "gene_symbol" and not symbols:
                        symbols.append(part)
                    elif actual_type == "refseq":
                        refseqs.append(part)
                    elif actual_type == "ensembl":
                        ensembls.append(part)
                    elif actual_type == "genbank":
                        genbanks.append(part)
                    elif actual_type == "entrez":
                        entrezs.append(part)
                    elif actual_type == "description":
                        descriptions.append(part)
                    elif actual_type == "location":
                        locations.append(part)
            else:
                # Classify positions without patterns by actual type
                actual_type = self.classify_element(part)
                if actual_type == "gene_symbol" and not symbols:
                    symbols.append(part)
                elif actual_type == "refseq":
                    refseqs.append(part)
                elif actual_type == "ensembl":
                    ensembls.append(part)
                elif actual_type == "genbank":
                    genbanks.append(part)
                elif actual_type == "entrez":
                    entrezs.append(part)
                elif actual_type == "description":
                    descriptions.append(part)
                elif actual_type == "location":
                    locations.append(part)

        # Compile results
        if symbols:
            result['gene_symbol'] = symbols[0]
        if refseqs:
            result['refseq_id'] = "; ".join(refseqs)
        if ensembls:
            result['ensembl_id'] = "; ".join(ensembls)
        if genbanks:
            result['genbank_id'] = "; ".join(genbanks)
        if entrezs:
            result['entrez_id'] = "; ".join(entrezs)
        if descriptions:
            result['description'] = " | ".join(descriptions)
        if locations:
            result['location'] = "; ".join(locations)

        return result

@st.cache_data
def convert_df(df):
    """Convert DataFrame to TSV format"""
    return df.to_csv(index=True, sep='\t').encode('utf-8')

@st.cache_data
def load_and_parse_file(file_content, file_name):
    """Load and parse file (with caching)"""
    delimiter = '\t' if '\t' in file_content.split('\n')[0] else ','
    df = pd.read_csv(
        io.StringIO(file_content),
        sep=delimiter,
        engine='python'
    )
    return df

@st.cache_data
def extract_gene_info(df_dict, gene_column, sample_size):
    """Extract gene information (with caching)"""
    # Restore DataFrame from dict
    df = pd.DataFrame(df_dict)

    extractor = GeneInfoExtractor()

    # Learn from sample data
    sample_size_actual = min(sample_size, len(df))
    sample_data = df[gene_column].head(sample_size_actual).dropna()
    extractor.learn_pattern(sample_data)

    if not extractor.learned:
        return None, None

    # Efficient batch processing
    batch_size = 1000
    results = []

    for i in range(0, len(df), batch_size):
        batch = df[gene_column].iloc[i:i+batch_size]
        batch_results = batch.apply(extractor.extract_with_pattern)
        results.extend(batch_results.tolist())

    # Expand extracted information into individual columns
    info_df = pd.DataFrame(results)

    # Remove empty columns
    cols_to_keep = []
    for col in info_df.columns:
        if info_df[col].str.len().gt(0).any():
            cols_to_keep.append(col)

    return info_df, cols_to_keep

def main():
    st.set_page_config(page_title="Microarray Gene Name Filter", page_icon="🧬", layout="wide")
    st.title("🧬 Microarray Gene Name Filter")
    st.markdown("---")

    st.markdown("""
    ### Features
    - Extract gene information from microarray data
    - Select gene name column and place it in the first column
    - Option to remove rows where gene name is None
    - Detection and aggregation of duplicate genes (Mean/Max)
    """)

    # File upload
    uploaded_file = st.file_uploader(
        "Select a TSV/CSV file",
        type=['tsv', 'txt', 'csv'],
        help="Microarray data file"
    )

    if uploaded_file is not None:
        try:
            # Load file (with caching)
            content = uploaded_file.getvalue().decode('utf-8')
            df = load_and_parse_file(content, uploaded_file.name)

            st.success(f"✅ File loaded successfully: {len(df):,} rows x {len(df.columns)} columns")

            # Data preview
            st.subheader("📊 Data Preview")

            preview_rows = st.slider(
                "Number of preview rows",
                min_value=5,
                max_value=min(100, len(df)),
                value=10
            )

            st.dataframe(df.head(preview_rows))

            # Gene column selection
            st.subheader("🔍 Gene Information Column Selection")

            # Find default column
            default_col = None
            for col in df.columns:
                if 'gene' in col.lower() or 'symbol' in col.lower() or 'description' in col.lower():
                    default_col = col
                    break

            gene_column = st.selectbox(
                "Select the column containing gene information",
                options=df.columns.tolist(),
                index=df.columns.tolist().index(default_col) if default_col else 0,
                help="Usually the 'gene_assignment' column contains gene information"
            )

            # Display sample of selected column
            if gene_column:
                st.info("Sample of selected column (first 5 rows):")
                sample_data = df[gene_column].head(5).to_frame()
                st.dataframe(sample_data)

                # Learning sample size settings
                st.subheader("⚙️ Learning Settings")

                col1, col2 = st.columns(2)
                with col1:
                    sample_size = st.number_input(
                        "Number of samples for learning",
                        min_value=20,
                        max_value=min(500, len(df)),
                        value=min(100, len(df)),
                        help="More samples improve accuracy but increase processing time"
                    )

                with col2:
                    st.info(f"""
                    📌 Recommended settings:
                    - Small data (<1000 rows): 50-100
                    - Medium data (1000-10000 rows): 100-200
                    - Large data (>10000 rows): 200-500
                    """)

                # Process execution button
                if st.button("🚀 Extract Gene Information", type="primary"):
                    st.subheader("📈 Learning Position Patterns and Extraction")

                    with st.spinner("Extracting gene information..."):
                        # Convert DataFrame to dict for caching
                        df_dict = df.to_dict('list')
                        info_df, cols_to_keep = extract_gene_info(df_dict, gene_column, sample_size)

                    if info_df is None:
                        st.error("Pattern learning failed")
                    else:
                        # Save to session state
                        st.session_state['info_df'] = info_df
                        st.session_state['cols_to_keep'] = cols_to_keep
                        st.session_state['original_df'] = df.copy()
                        st.success("✨ Extraction complete!")

                # If extraction results exist
                if 'info_df' in st.session_state and 'cols_to_keep' in st.session_state:
                    info_df = st.session_state['info_df']
                    cols_to_keep = st.session_state['cols_to_keep']
                    df = st.session_state['original_df'].copy()

                    # Display statistics
                    st.subheader("📈 Extraction Results Statistics")

                    cols = st.columns(len(cols_to_keep))
                    for i, col_name in enumerate(cols_to_keep):
                        with cols[i]:
                            count = info_df[col_name].str.len().gt(0).sum()
                            percentage = count / len(info_df) * 100

                            label = {
                                'gene_symbol': "Gene Symbol",
                                'refseq_id': "RefSeq ID",
                                'ensembl_id': "Ensembl ID",
                                'genbank_id': "GenBank ID",
                                'entrez_id': "Entrez ID",
                                'description': "Description",
                                'location': "Location"
                            }.get(col_name, col_name)

                            st.metric(
                                label,
                                f"{count:,} rows",
                                f"{percentage:.1f}%"
                            )

                    # Select which gene name column to use
                    st.subheader("🎯 Select Gene Name to Place in First Column")

                    # Display available gene name types
                    available_types = []
                    for col_name in cols_to_keep:
                        if col_name in ['gene_symbol', 'refseq_id', 'ensembl_id', 'genbank_id', 'entrez_id']:
                            count = info_df[col_name].str.len().gt(0).sum()
                            if count > 0:
                                available_types.append(col_name)

                    if available_types:
                        selected_gene_col = st.selectbox(
                            "Select gene name type to place in first column",
                            options=available_types,
                            format_func=lambda x: {
                                'gene_symbol': "Gene Symbol",
                                'refseq_id': "RefSeq ID",
                                'ensembl_id': "Ensembl ID",
                                'genbank_id': "GenBank ID",
                                'entrez_id': "Entrez ID"
                            }.get(x, x),
                            help="The selected gene name will become the row index"
                        )

                        # Option to remove None rows
                        remove_none = st.checkbox(
                            "Remove rows with empty gene names",
                            value=True,
                            help="Removes rows where the selected gene name type value is empty"
                        )

                        if st.button("✅ Place Gene Name in First Column", type="primary"):
                            # Get selected gene name column
                            gene_names = info_df[selected_gene_col].copy()

                            # Prepare other extracted information columns
                            other_cols = [col for col in cols_to_keep if col != selected_gene_col]
                            other_info = info_df[other_cols].copy()

                            # Combine with original data
                            result_df = pd.concat([other_info, df], axis=1)
                            result_df.index = gene_names
                            result_df.index.name = "Gene"

                            # Remove None rows
                            if remove_none:
                                before_len = len(result_df)
                                result_df = result_df[result_df.index.str.len() > 0]
                                removed_count = before_len - len(result_df)
                                st.info(f"Removed {removed_count} rows with empty gene names")

                            # Save to session state
                            st.session_state['result_df'] = result_df
                            st.session_state['selected_gene_col'] = selected_gene_col

                            st.success("✨ Gene name placed in first column!")

                        # If results exist
                        if 'result_df' in st.session_state:
                            result_df = st.session_state['result_df']

                            st.subheader("📋 Processing Results Preview")
                            st.write(f"Data length: {len(result_df):,} rows")
                            st.dataframe(result_df.head(preview_rows))

                            # Detect duplicate genes
                            st.subheader("🔍 Duplicate Gene Detection")

                            dup_d = result_df.loc[result_df.index.duplicated(keep=False), :].sort_index()
                            dup_count = len(set(dup_d.index))

                            if dup_count > 0:
                                st.warning(f"⚠️ Number of duplicate genes: {dup_count}")
                                st.dataframe(dup_d)

                                # Aggregate options
                                st.subheader("📊 Duplicate Gene Aggregation")

                                agg_method = st.radio(
                                    "Select aggregation method:",
                                    ('Mean', 'Max'),
                                    help="Mean: Average value, Max: Maximum value"
                                )

                                if st.button("🔄 Aggregate Duplicate Genes", type="primary"):
                                    dup_gene = set(result_df.index[result_df.index.duplicated(keep=False)])
                                    df_nodup = result_df[~result_df.index.duplicated(keep='first')]
                                    grouping = result_df.groupby(level=0)

                                    if agg_method == "Mean":
                                        df_agg = grouping.mean(numeric_only=True)
                                    else:
                                        df_agg = grouping.max(numeric_only=True)

                                    for gene in dup_gene:
                                        # Update numeric columns only
                                        df_nodup.loc[gene, df_agg.columns] = df_agg.loc[gene, df_agg.columns].values

                                    st.success(f"✨ Aggregated using {agg_method} method")
                                    st.write("Duplicate genes after aggregation:")
                                    st.dataframe(df_nodup.loc[list(dup_gene), :].sort_index())

                                    # Update session state
                                    st.session_state['result_df'] = df_nodup
                                    st.session_state['aggregated'] = True

                                    st.rerun()
                            else:
                                st.success("✅ No duplicate genes found")

                            # Download
                            st.markdown("---")
                            st.subheader("💾 File Download")

                            st.write(f"Final data length: {len(result_df):,} rows")

                            original_name = uploaded_file.name.rsplit('.', 1)[0]
                            output_filename = f"{original_name}_filtered.tsv"

                            csv_data = convert_df(result_df)
                            st.download_button(
                                label="📥 Download Processed File (TSV)",
                                data=csv_data,
                                file_name=output_filename,
                                mime="text/tab-separated-values"
                            )
                    else:
                        st.warning("No gene name information was extracted")

        except Exception as e:
            st.error(f"❌ An error occurred: {str(e)}")
            st.markdown("**Debug information:**")
            st.code(str(e))

    # Footer
    st.markdown("---")
    st.markdown("""
    <div style='text-align: center; color: gray; font-size: 0.8em;'>
    Microarray Gene Name Filter v1.0 |
    Extract, filter, and aggregate gene names from microarray data
    </div>
    """, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
