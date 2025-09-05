import streamlit as st
import pandas as pd
import numpy as np
import io
import zipfile
import re
import sys
import matplotlib.pyplot as plt

def determine_file_type(filename):
    if filename.endswith('.bed'):
        return 'bed'
    elif 'narrowPeak' in filename:
        return 'narrowPeak'
    elif 'broadPeak' in filename:
        return 'broadPeak'
    elif filename.endswith('.txt'):
        return 'homer'
    return 'unknown'

def get_column_names(file_type):
    if file_type == 'bed':
        return ['chrom', 'start', 'end', 'name', 'score', 'strand']
    elif file_type == 'narrowPeak':
        return ['chrom', 'start', 'end', 'name', 'score', 'strand', 'FC', '-log10p', '-log10q', 'summit_pos']
    elif file_type == 'broadPeak':
        return ['chrom', 'start', 'end', 'name', 'score', 'strand', 'FC', '-log10p', '-log10q']
    else:
        return None  # For homer or unknown file types

def has_header(first_line, file_type):
    if file_type == 'homer':
        return True
    elif file_type in ['bed', 'narrowPeak', 'broadPeak']:
        fields = first_line.split('\t')
        # MACS3 and BED files typically start with chromosome names
        return not (fields[0].lower().startswith('chr') or fields[0].isdigit())
    return False

def identify_columns(df):
    numeric_columns = df.select_dtypes(include=[int, float]).columns
    if len(numeric_columns) >= 2:
        return numeric_columns[0], numeric_columns[1]
    else:
        raise ValueError("Unable to identify numeric columns for start and end positions")

def calculate_fragment_length(row, start_col, end_col):
    return int(row[end_col]) - int(row[start_col])

def process_file(file, file_type, use_fragment_filter, fragment_min, fragment_max, use_value_filter, value_column, value_min, value_max):
    file_type = determine_file_type(file.name)
    first_line = file.readline().decode('utf-8').strip()
    file.seek(0)  # Reset file pointer to the beginning

    header_present = has_header(first_line, file_type)
    columns = get_column_names(file_type)

#    st.write(f"Debug: File type: {file_type}, Header present: {header_present}")

    df = pd.read_csv(file, sep='\t', header=None, names=columns)

#    st.write(f"Debug: Columns in dataframe: {df.columns.tolist()}")
#    st.write(f"Debug: First few rows of data:")
    st.write(df.head(3))
#    st.write(f"Debug: Data types: {df.dtypes}")

    # Convert all columns to numeric where possible
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='ignore')

    # Identify the correct column names for start and end
    start_col, end_col = 'start', 'end'

    df['fragment_length'] = df.apply(lambda row: calculate_fragment_length(row, start_col, end_col), axis=1)

    if use_fragment_filter:
        initial_rows = len(df)
        df = df[(df['fragment_length'] >= fragment_min) & (df['fragment_length'] <= fragment_max)]
        st.write(f"Rows after fragment filter: {len(df)} (removed {initial_rows - len(df)} rows)")

    if use_value_filter and value_column in df.columns:
        initial_rows = len(df)
        df[value_column] = pd.to_numeric(df[value_column], errors='coerce')
        if value_min is not None:
            df = df[df[value_column] >= value_min]
        if value_max is not None:
            df = df[df[value_column] <= value_max]
        st.write(f"Rows after value filter: {len(df)} (removed {initial_rows - len(df)} rows)")
    elif use_value_filter:
        st.write(f"Value filter not applied. Column '{value_column}' not found in dataframe.")

    return df, file_type

def get_common_columns(dataframes):
    if not dataframes:
        return []
    common_columns = set(dataframes[0].columns)
    for df in dataframes[1:]:
        common_columns = common_columns.intersection(set(df.columns))
    return list(common_columns)

def get_filterable_columns(dataframes):
    common_columns = get_common_columns(dataframes)
    numeric_columns = [col for col in common_columns if all(pd.api.types.is_numeric_dtype(df[col]) for df in dataframes)]
    # Exclude 'chrom', 'start', 'end', and 'strand' columns
    return [col for col in numeric_columns if col.lower() not in ['chrom', 'start', 'end', 'strand']]

def calculate_initial_values(dataframes, column):
    all_values = pd.concat([df[column] for df in dataframes])
    q1 = all_values.quantile(0.25)
    q3 = all_values.quantile(0.75)
    iqr = q3 - q1
    lower_bound = max(all_values.min(), q1 - 1.5 * iqr)
    upper_bound = min(all_values.max(), q3 + 1.5 * iqr)
    return lower_bound, upper_bound

def plot_histogram_multiple(dataframes, column, title, file_names):
    fig, ax = plt.subplots(figsize=(10, 6))
    for df, name in zip(dataframes, file_names):
        ax.hist(df[column].dropna(), bins=50, alpha=0.5, label=name)
    ax.set_title(title)
    ax.set_xlabel(column)
    ax.set_ylabel('Frequency')
    ax.legend()
    return fig

st.title('Bed, macs and Homer peak file filter for peak size and other parameters')

uploaded_files = st.file_uploader("Upload your files", accept_multiple_files=True, type=['bed', 'txt', 'narrowPeak', 'broadPeak'])

if uploaded_files:
    st.header("Filtering Options")

    # Process all files without filters first
    processed_dfs = []
    file_types = []
    for file in uploaded_files:
        df, file_type = process_file(file, None, False, 0, sys.maxsize, False, None, None, None)
        processed_dfs.append(df)
        file_types.append(file_type)

    filterable_columns = get_filterable_columns(processed_dfs)

    use_fragment_filter = st.checkbox('##### Apply peak length filter')
    if use_fragment_filter:
        fragment_min, fragment_max = calculate_initial_values(processed_dfs, 'fragment_length')
        col1, col2 = st.columns(2)
        with col1:
            fragment_min = st.number_input('Minimum fragment length (bp)', value=int(fragment_min), step=1)
        with col2:
            fragment_max = st.number_input('Maximum fragment length (bp)', value=int(fragment_max), step=1)
    else:
        fragment_min, fragment_max = 0, sys.maxsize

    use_value_filter = st.checkbox('##### Apply additional filter')
    if use_value_filter and filterable_columns:
        value_column = st.selectbox('Select value column for filtering', filterable_columns)
        lower_bound, upper_bound = calculate_initial_values(processed_dfs, value_column)
        col1, col2 = st.columns(2)
        with col1:
            value_min = st.number_input('Minimum value', value=float(lower_bound), format="%.2f")
        with col2:
            value_max = st.number_input('Maximum value', value=float(upper_bound), format="%.2f")
    else:
        value_column, value_min, value_max = None, None, None

    # Re-process all files with filters
    filtered_dfs = []
    for file, file_type in zip(uploaded_files, file_types):
        df, _ = process_file(file, file_type, use_fragment_filter, fragment_min, fragment_max,
                             use_value_filter, value_column, value_min, value_max)
        filtered_dfs.append(df)

    # Display summary of all files
    st.header("Summary of Processed Files")
    summary_data = []
    for file, df in zip(uploaded_files, filtered_dfs):
        summary_data.append({
            "File Name": file.name,
            "File Type": determine_file_type(file.name),
            "Rows": len(df),
            "Columns": len(df.columns)
        })
    st.table(pd.DataFrame(summary_data))

    # Plot histograms
    if use_fragment_filter:
        st.subheader("Fragment Length Distribution")
        fig = plot_histogram_multiple(filtered_dfs, 'fragment_length', 'Fragment Length Distribution', [f.name for f in uploaded_files])
        st.pyplot(fig)

    if use_value_filter and value_column:
        st.subheader(f"{value_column} Distribution")
        fig = plot_histogram_multiple(filtered_dfs, value_column, f'{value_column} Distribution', [f.name for f in uploaded_files])
        st.pyplot(fig)

    # Download filtered results
    if st.button('Prepare filtered results for download'):
        output = io.BytesIO()
        with zipfile.ZipFile(output, 'w') as zf:
            for file, df in zip(uploaded_files, filtered_dfs):
                csv_buffer = io.StringIO()
                df.to_csv(csv_buffer, sep='\t', index=False, header=True)
                zf.writestr(f'filtered_{file.name}', csv_buffer.getvalue())

        output.seek(0)
        st.download_button(
            label="Download ZIP",
            data=output,
            file_name="filtered_files.zip",
            mime="application/zip"
        )

else:
    st.write("Please upload one or more files to begin.")