import streamlit as st
import pandas as pd
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.vectors import StrVector
from rpy2.robjects.packages import importr
import io
import base64
import tempfile
import os
from scipy import stats

# Import R packages
pandas2ri.activate()
impulsede2 = importr('ImpulseDE2')
base = importr('base')

@st.cache_data
def load_data(uploaded_file):
    if uploaded_file.name.endswith('.csv'):
        df = pd.read_csv(uploaded_file, index_col=0)
    elif uploaded_file.name.endswith('.tsv'):
        df = pd.read_csv(uploaded_file, sep='\t', index_col=0)
    elif uploaded_file.name.endswith('.txt'):
        df =  pd.read_csv(uploaded_file, index_col = 0, sep=None, engine='python') # Auto detect
    elif uploaded_file.name.endswith(('.xls', '.xlsx')):
        df = pd.read_excel(uploaded_file, index_col=0)
    else:
        st.error("Unsupported file format. Please upload a CSV, TSV, or Excel file.")
        return None
    return df

def filter_genes(count_data, threshold):
    """
    Filtering function to keep only genes with values above the specified threshold
    """
    filtered_data = count_data[(count_data >= threshold).all(axis=1)]
    return filtered_data

def run_impulsede2(count_data, annotation, case_ctrl, identify_transients, confounders):
    try:
        # Create temporary files
        with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as temp_annotation:
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.csv') as temp_counts:
                # Save data to temporary files
                annotation.to_csv(temp_annotation.name, index=False)
                count_data.to_csv(temp_counts.name)

                # R code to run ImpulseDE2
                r_code = f"""
                library(ImpulseDE2)

                # Read data
                dfAnnotation <- read.csv("{temp_annotation.name}")
                matCountData <- as.matrix(read.csv("{temp_counts.name}", row.names=1))

                # Run ImpulseDE2
                objectImpulseDE2 <- runImpulseDE2(
                  matCountData           = matCountData,
                  dfAnnotation           = dfAnnotation,
                  boolCaseCtrl           = {"TRUE" if case_ctrl else "FALSE"},
                  vecConfounders         = {confounders if confounders else "NULL"},
                  boolIdentifyTransients = {"TRUE" if identify_transients else "FALSE"},
                  scaNProc               = 10 )

                # Save results
                write.csv(objectImpulseDE2$dfImpulseDE2Results, "impulsede2_results.csv", row.names = TRUE)
                """

                # Execute R code
                robjects.r(r_code)

        # Read results
        results_df = pd.read_csv("impulsede2_results.csv", index_col=0)

        # Delete temporary files
        os.unlink(temp_annotation.name)
        os.unlink(temp_counts.name)
        os.unlink("impulsede2_results.csv")

        return results_df
    except Exception as e:
        st.error(f"An error occurred while running ImpulseDE2: {str(e)}")
        return None

def plot_genes(results_df, counts_df, annotation_df, n_top_ids=5):
    # Sort genes by adjusted p-value and select top genes
    top_genes = results_df.sort_values('padj').head(n_top_ids)

    # Create a line plot for the top genes
    fig = go.Figure()
    for gene in top_genes.index:
        if gene in counts_df.index:
            y_values = counts_df.loc[gene].values
            x_values = annotation_df['Time'].values
            fig.add_trace(go.Scatter(x=x_values, y=y_values, mode='lines+markers', name=gene))

    fig.update_layout(
        title="Top Genes Expression",
        xaxis_title="Time",
        yaxis_title="Expression",
        legend_title="Genes",
        hovermode="x unified"
    )

    return fig

def plot_heatmap(results_df, counts_df, annotation_df):
    # Select top 50 differentially expressed genes
    top_genes = results_df.sort_values('padj').head(50).index
    expression_data = counts_df.loc[top_genes]

    # Calculate Z-scores for each gene
    z_scores = expression_data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)

    # Sort columns by time
    sorted_times = annotation_df.sort_values('Time')
    z_scores = z_scores[sorted_times['Sample']]  # Use 'Sample' column to reorder

    fig = go.Figure(data=go.Heatmap(
        z=z_scores.values,
        x=sorted_times['Time'].values,
        y=z_scores.index,
        colorscale='RdBu_r'
    ))

    fig.update_layout(
        title="Heatmap of Top 50 Differentially Expressed Genes (Z-scores)",
        xaxis_title="Time",
        yaxis_title="Genes",
        height=800
    )

    return fig

def main():
    st.title("ImpulseDE2 time-course data analysis")

    uploaded_file = st.file_uploader("Upload count data file (CSV, TSV, Excel)", type=["csv", "tsv", "txt", "xls", "xlsx"])

    if uploaded_file is not None:
        count_data = load_data(uploaded_file)
        if count_data is not None:
            st.write("Count data preview:")
            st.write(count_data.head())

            # Gene filtering option
            st.subheader("Gene Filtering")
            filter_threshold = st.number_input("Enter threshold for gene filtering", min_value=0, value=0, step=1)

            if st.button("Apply Filter"):
                original_gene_count = count_data.shape[0]
                filtered_count_data = filter_genes(count_data, filter_threshold)
                filtered_gene_count = filtered_count_data.shape[0]

                st.write(f"Original number of genes: {original_gene_count}")
                st.write(f"Number of genes after filtering: {filtered_gene_count}")

                count_data = filtered_count_data  # Update the count_data with filtered data

            # Case-control analysis option
            case_ctrl = st.checkbox("Case-control differential expression analysis")

            # Identify transients option
            identify_transients = st.checkbox("Identify transiently regulated genes")

            # Confounders input
            confounders = st.text_input("Confounders (comma-separated, leave blank if none)", "")
            confounders = [c.strip() for c in confounders.split(',')] if confounders else None

            # Metadata input
            st.subheader("Time and condition")
            n_samples = count_data.shape[1]

            # Initialize metadata with default values
            metadata = pd.DataFrame({
                "Sample": count_data.columns,
                "Time": range(1, n_samples + 1),
                "Condition": ["case"] * n_samples
            })

            # Add confounders to metadata if specified
            if confounders:
                for confounder in confounders:
                    metadata[confounder] = [''] * n_samples

            if 'inpulse_submitted' not in st.session_state:
                st.session_state.inpulse_submitted = False

            with st.form("input_time and condition"):

                # Create an editable dataframe for metadata
                edited_metadata = st.data_editor(
                    metadata,
                    column_config={
                        "Sample": st.column_config.TextColumn("Sample", disabled=True),
                        "Time": st.column_config.NumberColumn("Time", min_value=0, format="%d"),
                        "Condition": st.column_config.SelectboxColumn("Condition", options=["case", "control"]) if case_ctrl else None,
                    },
                    hide_index=True,
                    num_rows="fixed"
                )
                st.session_state.inpulse_submitted = True
                submitted = st.form_submit_button("Submit")

            if st.session_state.inpulse_submitted:

                if st.button("Run ImpulseDE2"):
                    if not case_ctrl:
                        edited_metadata['Condition'] = 'case'  # Set all to 'case' if not case-control

                    results_df = run_impulsede2(count_data, edited_metadata, case_ctrl, identify_transients, confounders)

                    if results_df is not None:
                        st.subheader("Results")
                        st.write(results_df)

                        # Plot genes
                        st.subheader("Gene Expression Plot")
                        gene_plot = plot_genes(results_df, count_data, edited_metadata)
                        st.plotly_chart(gene_plot)

                        # Plot heatmap
                        st.subheader("Heatmap")
                        heatmap = plot_heatmap(results_df, count_data, edited_metadata)
                        st.plotly_chart(heatmap)

                        # Download results
                        csv = results_df.to_csv(index=True)
                        b64 = base64.b64encode(csv.encode()).decode()
                        href = f'<a href="data:file/csv;base64,{b64}" download="impulsede2_results.csv">Download results as CSV</a>'
                        st.markdown(href, unsafe_allow_html=True)

                        st.markdown("""
#### Main Result Columns
- **Gene**: Gene ID
- **p**: P-value for differential expression
- **padj**: Benjamini-Hochberg FDR-adjusted P-value
- **loglik_full**: Log-likelihood of full model
- **loglik_red**: Log-likelihood of reduced model
- **df_full**: Degrees of freedom of full model
- **df_red**: Degrees of freedom of reduced model
- **mean**: Estimated mean parameter of constant model for the first batch
- **allZero**: Whether the gene had no non-zero observations (if TRUE, fitting and DE analysis are skipped, and entry is NA)

#### Additional Columns for Case-only DE Analysis
- **converge_impulse**: Convergence status of impulse model fit (full model)
- **converge_const**: Convergence status of constant model fit (reduced model)

#### Additional Columns for Case-control DE Analysis
- **converge_combined**: Convergence status of impulse model fit combining case and control samples (reduced model)
- **converge_case**: Convergence status of impulse model fit for case condition samples (full model 1/2)
- **converge_control**: Convergence status of impulse model fit for control condition samples (full model 2/2)

#### Transient Expression Identification (when boolIdentifyTransients = TRUE)
- **converge_sigmoid**: Convergence status of sigmoid model fit for case condition samples
- **impulseTOsigmoid_p**: P-value of log-likelihood ratio test between impulse and sigmoid model fits for case condition samples
- **impulseTOsigmoid_padj**: Benjamini-Hochberg corrected P-value of the above
- **sigmoidTOconst_p**: P-value of log-likelihood ratio test between sigmoid and constant model fits for case condition samples
- **sigmoidTOconst_padj**: Benjamini-Hochberg corrected P-value of the above
- **isTransient**: Whether the gene is transiently activated or deactivated and differentially expressed
- **isMonotonous**: Whether the gene is differentially expressed but not transiently (monotonic increase or decrease in expression level)
""")

    else:
        st.info("Please upload a count data file to begin.")

if __name__ == "__main__":
    main()
