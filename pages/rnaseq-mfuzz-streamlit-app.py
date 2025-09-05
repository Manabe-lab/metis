import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import io
import zipfile
from sklearn.preprocessing import StandardScaler
import tempfile
import os
import re
from PIL import Image
import psutil

# Monkey patch for pandas2ri to use items() instead of iteritems()
from rpy2.robjects import conversion
@conversion.py2rpy.register(pd.DataFrame)
def py2rpy_pandasdataframe(obj):
    from rpy2.robjects.vectors import DataFrame, ListVector, StrVector
    from rpy2.rinterface import SexpVector

    if obj.index.name:
        obj = obj.reset_index()
    columns = StrVector(obj.columns)
    items = [conversion.py2rpy(obj.iloc[:, i]) for i in range(obj.shape[1])]
    return DataFrame(ListVector(zip(columns, items)))

def get_memory_usage():
    process = psutil.Process(os.getpid())
    return process.memory_info().rss / 1024 ** 2  # メモリ使用量をMB単位で返す

# Function to install R packages
def install_r_packages(packages):
    utils = importr('utils')
    for package in packages:
        if not robjects.r(f'require({package}, quietly = TRUE)')[0]:
            st.warning(f"Installing R package: {package}")
            utils.install_packages(package, repos="http://cran.us.r-project.org")

# Install required R packages
required_packages = ['e1071', 'Mfuzz', 'ggplot2']
install_r_packages(required_packages)

# Initialize R
robjects.r('options(warn=-1)')
robjects.r('suppressMessages(library(Mfuzz))')
robjects.r('suppressMessages(library(ggplot2))')

# Enable automatic conversion between pandas and R dataframes
pandas2ri.activate()

@st.cache_data
def load_data(file):
    if file.name.endswith('.tsv'):
        df = pd.read_csv(file, sep='\t', index_col=0)
    elif file.name.endswith('.csv'):
        df = pd.read_csv(file, index_col=0)
    elif file.name.endswith('.xlsx'):
        df = pd.read_excel(file, index_col=0)
    else:
        st.error("Unsupported file format. Please upload a TSV, CSV, or Excel file.")
        return None
    return df

def create_time_input_table(sample_names):

    data_df = pd.DataFrame({'Sample': sample_names, 'Time': [0.0] * len(sample_names)})
    if 'time_data' not in st.session_state:
        st.session_state.time_data = data_df   
    
    edited_df = st.data_editor(
        data_df,
        num_rows="fixed",
        hide_index=True,
        column_config={
            "Sample": st.column_config.TextColumn("Sample", disabled=True),
            "Time": st.column_config.NumberColumn("Time", min_value=0.0, format="%.2f")
        },
        key="time_input"
    )
    
    if not edited_df.equals(st.session_state.time_data):
        st.session_state.time_data = edited_df
    
    return st.session_state.time_data

def preprocess_data(data, remove_low_expression, log_transform, normalize, low_expression_threshold):
#    if remove_low_expression:
#        data = data[data.mean(axis=1) > low_expression_threshold]

    if log_transform:
        data = np.log2(data + 1)

    if normalize:
        # Gene-wise standardization
        data = (data.T - data.T.mean()) / data.T.std()
        data = data.T

    return data

def clean_column_name(name):
    # Remove special characters and spaces, replace with underscores
    cleaned = re.sub(r'[^\w\s]', '', str(name))
    cleaned = re.sub(r'\s+', '_', cleaned)
    # Ensure the name starts with a letter or underscore (R requirement)
    if not cleaned[0].isalpha() and cleaned[0] != '_':
        cleaned = 'X' + cleaned
    return cleaned

def prepare_data_for_r(data):
    # Create a copy of the dataframe
    cleaned_data = data.copy()
    
    # Clean column names
    cleaned_data.columns = [clean_column_name(col) for col in cleaned_data.columns]
    
    # Clean index names
    cleaned_data.index = [clean_column_name(idx) for idx in cleaned_data.index]
    
    # Create mappings for later reference
    column_mapping = dict(zip(cleaned_data.columns, data.columns))
    index_mapping = dict(zip(cleaned_data.index, data.index))
    
    return cleaned_data, column_mapping, index_mapping

def calculate_optimal_m_python(data):
    N, D = data.shape
    m = 1 + (1418/N + 22.05) * D**(-2) + (12.33/N + 0.243) * D**(-0.0406 * np.log(N) - 0.1134)
    return m

def calculate_optimal_m_r(data):
    r_code = """
    function(data) {
        library(Mfuzz)
        eset <- ExpressionSet(assayData=as.matrix(data))
        m <- mestimate(eset)
        return(m)
    }
    """
    r_func = robjects.r(r_code)
    
    # Convert pandas DataFrame to R matrix
    r_matrix = robjects.r.matrix(robjects.FloatVector(data.values.flatten()), nrow=data.shape[0], ncol=data.shape[1])
    
    optimal_m = r_func(r_matrix)[0]
    return optimal_m

def calculate_optimal_m(data):
    m_python = calculate_optimal_m_python(data)
    m_r = calculate_optimal_m_r(data)
    return m_python, m_r

def plot_mfuzz_results_python(data, time_points, clusters, membership):
    n_clusters = len(set(clusters))
    fig, axes = plt.subplots(n_clusters, 1, figsize=(10, 5*n_clusters))

    if n_clusters == 1:
        axes = [axes]

    for i in range(n_clusters):
        cluster_genes = clusters[clusters == i+1].index
        cluster_data = data.loc[cluster_genes]

        # Calculate mean and standard deviation for the cluster
        mean = cluster_data.mean()
        std = cluster_data.std()

        # Plot mean line
        axes[i].plot(time_points, mean, color='blue', linewidth=2, label='Mean')
        
        # Plot confidence interval
        axes[i].fill_between(time_points, mean-std, mean+std, alpha=0.3, color='blue', label='±1 SD')

        # Plot individual gene trajectories (first 50 genes)
        for gene in cluster_genes[:50]:
            axes[i].plot(time_points, cluster_data.loc[gene], alpha=0.1, color='gray')

        axes[i].set_title(f"Cluster {i+1}")
        axes[i].set_xlabel("Time")
        axes[i].set_ylabel("Expression")
        axes[i].legend()

    plt.tight_layout()
    return fig

def run_mfuzz(data, time_points, c, m, remove_low_expression, log_transform, normalize, low_expression_threshold):
    # Preprocess data
    data_preprocessed = preprocess_data(data, remove_low_expression, log_transform, normalize, low_expression_threshold)
    
    # Prepare data for R
    cleaned_data, column_mapping, index_mapping = prepare_data_for_r(data_preprocessed)
    
    # Create temporary files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as temp_data_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as temp_time_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as temp_result_file, \
         tempfile.NamedTemporaryFile(mode='w', suffix='.png', delete=False) as temp_plot_file:

        # Save cleaned data and time points to temporary files
        cleaned_data.to_csv(temp_data_file.name)
        pd.DataFrame({'time': time_points}).to_csv(temp_time_file.name, index=False)

        # R code to run Mfuzz and generate plot
        r_code = f"""
        library(Mfuzz)
        data <- read.csv("{temp_data_file.name}", row.names=1)
        time <- read.csv("{temp_time_file.name}")$time
        
        cat("Data dimensions:", dim(data)[1], "x", dim(data)[2], "\\n")
        cat("First 5 rows and 5 columns of the data:\\n")
        print(head(data))
        cat("Time points:", time, "\\n")
        
        eset <- ExpressionSet(assayData=as.matrix(data))
        pData(eset)$time <- time
        result <- mfuzz(eset, c={c}, m={m})
        
        # Generate plot
        png("{temp_plot_file.name}", width=1500, height=1500, res=150)
        mfuzz.plot(eset, cl=result, mfrow=c(ceiling({c}/2), 2), time.labels=time, new.window=FALSE)
        dev.off()
        
        write.csv(result$membership, "{temp_result_file.name}")
        write.csv(result$cluster, "{temp_result_file.name}.cluster")
        
        # Display cluster information
        cat("Cluster sizes:\\n")
        print(table(result$cluster))
        cat("\\nFirst few rows of membership matrix:\\n")
        print(head(result$membership))
        """

        # Run R code and capture output
        r_output = robjects.r(r_code)
        
        # Display R output
        st.subheader("R Mfuzz Results:")
        st.text(str(r_output))

        # Read results
        try:
            membership = pd.read_csv(temp_result_file.name, index_col=0)
            clusters = pd.read_csv(f"{temp_result_file.name}.cluster", index_col=0, header=None, names=['Cluster'])
            
            # Remove any non-numeric cluster assignments and convert to integer
            clusters = clusters[pd.to_numeric(clusters['Cluster'], errors='coerce').notnull()]
            clusters['Cluster'] = clusters['Cluster'].astype(int)
            
            # Convert clusters to a Series
            clusters = clusters['Cluster']

            # Ensure membership and clusters have the same index
            common_index = membership.index.intersection(clusters.index)
            membership = membership.loc[common_index]
            clusters = clusters.loc[common_index]

            # Debug information
            st.write("Membership shape:", membership.shape)
            st.write("Clusters shape:", clusters.shape)
            st.write("Clusters head:", clusters.head())
            st.write("Clusters value counts:", clusters.value_counts())
            st.write("Clusters dtype:", clusters.dtype)

            # Map results back to original index names
            membership.index = [index_mapping.get(idx, idx) for idx in membership.index]
            clusters.index = [index_mapping.get(idx, idx) for idx in clusters.index]

        except Exception as e:
            st.error(f"Error reading Mfuzz results: {str(e)}")
            return None, None, None, None

    # Clean up temporary files
    os.unlink(temp_data_file.name)
    os.unlink(temp_time_file.name)
    os.unlink(temp_result_file.name)
    os.unlink(f"{temp_result_file.name}.cluster")

    return membership, clusters, data_preprocessed, temp_plot_file.name

def plot_mfuzz_results_python(data, time_points, clusters, membership):
    n_clusters = len(set(clusters))
    fig, axes = plt.subplots(n_clusters, 1, figsize=(10, 5*n_clusters))

    if n_clusters == 1:
        axes = [axes]

    unique_time_points = sorted(set(time_points))
    
    for i in range(n_clusters):
        cluster_genes = clusters[clusters == i+1].index
        cluster_data = data.loc[cluster_genes]

        # Calculate mean and standard deviation for the cluster
        mean = cluster_data.mean()
        std = cluster_data.std()

        # Plot individual gene trajectories (first 50 genes)
        for gene in cluster_genes:
            axes[i].plot(unique_time_points, cluster_data.loc[gene], alpha=0.1, color='lightgray')

        # Plot mean line
        axes[i].plot(unique_time_points, mean, color='blue', linewidth=2, label='Mean')
        
        # Plot confidence interval
        axes[i].fill_between(unique_time_points, mean-std, mean+std, alpha=0.3, color='blue', label='±1 SD')

        axes[i].set_title(f"Cluster {i+1}")
        axes[i].set_xlabel("Time")
        axes[i].set_ylabel("Expression")
        axes[i].legend()

        # Set x-axis ticks to match the actual time points
        axes[i].set_xticks(unique_time_points)
        axes[i].set_xticklabels(unique_time_points)

    plt.tight_layout()
    return fig

def main():
    st.title("Mfuzz analysis")
    st.markdown("#### Please upload normallzed data. Replicated data will be averaged for each time point.")
#    st.sidebar.write(f"Current memory usage: {get_memory_usage():.2f} MB")

    uploaded_file = st.file_uploader("Choose a TSV, CSV, or Excel file", type=["tsv", "csv", "xlsx"])

    if uploaded_file is not None:
        data = load_data(uploaded_file)
        if data is not None:
            st.write("Data preview:")
            st.write(data.head())
            st.info("Note: For R compatibility, column names and index values will be cleaned by removing special characters and spaces. Original names will be restored in the final results.")

            st.subheader("Data Preprocessing Options")

            st.write(f"Original data shape: {data.shape}")
            remove_low_expression = st.checkbox("Remove low expression genes", value=True)
            low_expression_threshold = 1.0
            if remove_low_expression:
                low_expression_threshold = st.number_input("Low expression threshold", min_value=0.0, value=1.0, step=0.1)
                data = data[data.mean(axis=1) > low_expression_threshold]
                st.write(f"Data shape after removing low expression genes: {data.shape}")

            if 'submitted' not in st.session_state:
                st.session_state.submitted = False

            with st.form("input_time"):
                st.subheader("Enter time points for each sample:")
                time_df = create_time_input_table(data.columns)

                st.session_state.submitted = True

                submitted = st.form_submit_button("Submit time")

            if st.session_state.submitted:
                time_points = time_df['Time'].tolist()
                unique_time_points = sorted(set(time_points))

                # Average data for samples with the same time point
                data_averaged = pd.DataFrame(index=data.index)
                for time in unique_time_points:
                    cols = [col for col, t in zip(data.columns, time_points) if t == time]
                    data_averaged[time] = data[cols].mean(axis=1)

                st.write(f"Data shape after averaging replicates: {data_averaged.shape}")

                log_transform = st.checkbox("Apply log2(x+1) transformation", value=True)
                normalize = st.checkbox("Normalize data (gene-wise standardization)", value=True)

                if normalize:
                    st.info("Gene-wise standardization: Each gene's expression values will be centered to have mean 0 and scaled to have standard deviation 1 across all time points.")

                c = st.number_input("Number of clusters (c)", min_value=2, max_value=20, value=6)
                
                optimal_m_python, optimal_m_r = calculate_optimal_m(data_averaged)
                st.write(f"Calculated optimal Fuzzification parameter (m):")
                st.write(f"Python implementation: {optimal_m_python:.2f}")
                st.write(f"R implementation (using mestimate): {optimal_m_r:.2f}")
                
                st.info("""
                The optimal m value is calculated using the mestimate function from the Mfuzz package in R.
                The Python implementation is an approximation of this function.
                """)
                
                if abs(optimal_m_python - optimal_m_r) > 0.01:
                    st.warning("The Python and R implementations produced different results. This might be due to slight differences in the calculation method or data handling.")
                
                use_python_m = st.checkbox("Use Python-calculated m value", value=True)
                optimal_m = optimal_m_python if use_python_m else optimal_m_r
                
                if optimal_m > 3.0:
                    st.warning(f"The calculated optimal m value ({optimal_m:.2f}) is unusually high. This might indicate issues with your data or preprocessing. Consider using a default value or investigating your data further.")
                    use_default_m = st.checkbox("Use default m value (2.0)", value=True)
                    if use_default_m:
                        optimal_m = 2.0
                
                m = st.number_input("Fuzzification parameter (m)", min_value=1.1, max_value=10.0, value=float(optimal_m), step=0.1)
                
                if m > 3.0:
                    st.warning("The selected m value is higher than the usual range (1.1 to 3.0). This might lead to very fuzzy clustering results.")

                plot_option = st.radio("Select plot type", ["R plot", "Python plot"])


        if st.session_state.submitted and st.button("Run Mfuzz Analysis"):
            time_points = time_df['Time'].tolist()

            with st.spinner("Running Mfuzz analysis..."):
                membership, clusters, data_preprocessed, r_plot_file = run_mfuzz(data_averaged, unique_time_points, c, m, remove_low_expression, log_transform, normalize, low_expression_threshold)

            if membership is not None and clusters is not None:
                st.subheader("Mfuzz Results")
                st.write("Membership shape:", membership.shape)
                st.write("Clusters shape:", clusters.shape)
                st.write("Preprocessed data shape:", data_preprocessed.shape)
                st.write("Clusters head:", clusters.head())
                st.write("Clusters value counts:", clusters.value_counts())
                st.write("Clusters dtype:", clusters.dtype)



                if plot_option == "R plot":
                    st.subheader("R-generated Mfuzz Plot")
                    try:
                        r_plot_image = Image.open(r_plot_file)
                        st.image(r_plot_image, use_column_width=True)
                    except Exception as e:
                        st.error(f"Failed to display R plot: {str(e)}")
                        st.info("You can still download the plot from the results ZIP file.")
                else:
                    st.subheader("Python-generated Mfuzz Plot")
                    fig = plot_mfuzz_results_python(data_preprocessed, time_points, clusters, membership.values)
                    st.pyplot(fig)


                # Prepare results for download
                cluster_info = pd.DataFrame({
                    'Gene': data_preprocessed.index,
                    'Cluster': clusters,
                    'Membership': membership.values.max(axis=1)
                })


                # Create a ZIP file with results
                zip_buffer = io.BytesIO()
                with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
                    # Save cluster information
                    cluster_info_buffer = io.StringIO()
                    cluster_info.to_csv(cluster_info_buffer, index=False)
                    zip_file.writestr("cluster_info.csv", cluster_info_buffer.getvalue())

                    # Save plot
                    plot_buffer = io.BytesIO()
                    if plot_option == "R plot":
                        with open(r_plot_file, 'rb') as f:
                            plot_buffer.write(f.read())
                        zip_file.writestr("mfuzz_R_plot.png", plot_buffer.getvalue())
                    else:
                        fig.savefig(plot_buffer, format='pdf')
                        zip_file.writestr("mfuzz_plot.pdf", plot_buffer.getvalue())

                    # Save preprocessed data
                    data_buffer = io.StringIO()
                    data_preprocessed.to_csv(data_buffer)
                    zip_file.writestr("preprocessed_data.csv", data_buffer.getvalue())

                # Offer ZIP file for download
                st.download_button(
                    label="Download Results",
                    data=zip_buffer.getvalue(),
                    file_name="mfuzz_results.zip",
                    mime="application/zip"
                )
            else:
                st.error("Mfuzz analysis failed. Please check your input data and parameters.")

if __name__ == "__main__":
    main()