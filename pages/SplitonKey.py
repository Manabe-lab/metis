import streamlit as st
import pandas as pd
import zipfile
from io import BytesIO

# Function to load data
def load_data(file):
    file_extension = file.name.split('.')[-1].lower()
    if file_extension in ['txt', 'tsv']:
        return pd.read_csv(file, sep="\t")
    elif file_extension == 'csv':
        return pd.read_csv(file)
    else:
        st.error(f"Unsupported file format: {file_extension}")
        return None

# Function to split data based on selected keys
def split_data(data, clusters, data_key, cluster_key, cluster_column):
    grouped_data = {cluster: data[data[data_key].isin(clusters[clusters[cluster_column] == cluster][cluster_key])] 
                    for cluster in clusters[cluster_column].unique()}
    return grouped_data

# Function to create a ZIP file from split data
def create_zip(split_data_dict, original_file_name, key_name, file_extension):
    zip_buffer = BytesIO()
    with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zf:
        for cluster, df in split_data_dict.items():
            file_name = f"{original_file_name}_{key_name}_{cluster}.{file_extension}"
            if file_extension in ['txt', 'tsv']:
                file_data = df.to_csv(index=False, sep="\t")
            else:
                file_data = df.to_csv(index=False)
            zf.writestr(file_name, file_data)
    return zip_buffer

# Streamlit app
st.title('File Splitter by Gene Name and Cluster')

# Upload files
data_file = st.file_uploader("Upload the data file", type=["txt", "tsv", "csv"])
cluster_file = st.file_uploader("Upload the cluster file", type=["txt", "tsv", "csv"])

if data_file and cluster_file:
    # Load data
    data = load_data(data_file)
    clusters = load_data(cluster_file)

    if data is not None and clusters is not None:
        # Display data
        st.write("Data File:")
        st.dataframe(data.head())
        st.write("Cluster File:")
        st.dataframe(clusters.head())

        # Select keys for splitting
        data_key = st.selectbox("Select the column with gene names in the data file", data.columns)
        cluster_key = st.selectbox("Select the column with gene names in the cluster file", clusters.columns)
        cluster_column = st.selectbox("Select the column with cluster information in the cluster file", clusters.columns)

        # Split button
        if st.button("Split"):
            # Split data
            split_data_dict = split_data(data, clusters, data_key, cluster_key, cluster_column)

            # Display split data
            for cluster, df in split_data_dict.items():
                st.write(f"Cluster: {cluster}")
                st.dataframe(df.head())

            # Create ZIP file
            original_file_name = ''.join(data_file.name.rsplit('.', 1)[0])[:25]  # Get the base name of the original file without extension

            file_extension = data_file.name.split('.')[-1].lower()
            zip_buffer = create_zip(split_data_dict, original_file_name, cluster_column, file_extension)

            # Download ZIP file
            st.download_button(
                label="Download ZIP",
                data=zip_buffer.getvalue(),
                file_name="split_data.zip",
                mime="application/zip"
            )
