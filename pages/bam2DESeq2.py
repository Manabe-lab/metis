import streamlit as st
import pandas as pd
import numpy as np
import subprocess
import os
import glob
import zipfile
import tempfile
from scipy.stats.mstats import gmean

# Set up the Streamlit app
st.title("CUT&RUN greenlist and peak counts for DESeq2")
st.markdown("### Counts are normalized for greenlist")

# Function to run shell commands
def run_command(command):
    st.text(f"Executing command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    return output.decode("utf-8"), error.decode("utf-8")

# Function to index BAM files
def index_bam_file(bam_file):
    if os.path.exists(bam_file + '.bai'):
        st.info(f"Index file for {bam_file} already exists. Skipping indexing.")
        return True
    command = f"samtools index {bam_file}"
    output, error = run_command(command)
    if error:
        st.error(f"Error indexing {bam_file}: {error}")
        return False
    return True

# Function to calculate counts using greenlist
def calculate_greenlist_counts(bam_files, greenlist_bed, output_prefix):
    # Index BAM files
    for bam_file in bam_files:
        if not index_bam_file(bam_file):
            return False

    bam_files_str = " ".join(bam_files)
    command = f"multiBamSummary BED-file --BED {greenlist_bed} --smartLabels -e --centerReads -o {output_prefix}.npz -b {bam_files_str} --outRawCounts {output_prefix}_output"
    output, error = run_command(command)
    
    # Combine output and error, as multiBamSummary seems to use stderr for all output
    full_output = output + error
    st.text("multiBamSummary output:")
    st.text(full_output)
    
    if "Number of bins found:" not in full_output:
        st.error("Error in multiBamSummary. Expected output not found.")
        return False
    
    # Check if output file exists
    if not os.path.exists(f"{output_prefix}_output"):
        st.error(f"Output file {output_prefix}_output was not created.")
        return False
    
    # Process the output to create the glist_quant.tsv file
    st.text(f"Processing output file: {output_prefix}_output")
    try:
        with open(f"{output_prefix}_output", 'r') as f:
            with open(f"{output_prefix}_glist_quant.tsv", 'w') as out:
                header = f.readline().strip().split('\t')
                sample_names = header[3:]
                out.write("Chr_Start_End\t" + "\t".join(sample_names) + "\n")
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        new_line = f"{parts[0]}_{parts[1]}_{parts[2]}\t" + "\t".join(parts[3:]) + "\n"
                        out.write(new_line)
    except Exception as e:
        st.error(f"Error processing output file: {str(e)}")
        return False
    
    st.success("Greenlist counts calculated successfully.")
    return True

# Function to calculate counts for peak files
def calculate_peak_counts(bam_files, peak_file, output_prefix):
    # Index BAM files
    for bam_file in bam_files:
        if not index_bam_file(bam_file):
            return False

    bam_files_str = " ".join(bam_files)
    command = f"multiBamSummary BED-file --BED {peak_file} --smartLabels -e --centerReads -o {output_prefix}.npz -b {bam_files_str} --outRawCounts {output_prefix}_output"
    output, error = run_command(command)
    
    # Combine output and error, as multiBamSummary seems to use stderr for all output
    full_output = output + error
    st.text("multiBamSummary output:")
    st.text(full_output)
    
    if "Number of bins found:" not in full_output:
        st.error("Error in multiBamSummary. Expected output not found.")
        return False
    
    # Check if output file exists
    if not os.path.exists(f"{output_prefix}_output"):
        st.error(f"Output file {output_prefix}_output was not created.")
        return False
    
    # Process the output to create the peak_quant.tsv file
    st.text(f"Processing output file: {output_prefix}_output")
    try:
        with open(f"{output_prefix}_output", 'r') as f:
            with open(f"{output_prefix}_peak_quant.tsv", 'w') as out:
                header = f.readline().strip().split('\t')
                sample_names = header[3:]
                out.write("Chr_Start_End\t" + "\t".join(sample_names) + "\n")
                for line in f:
                    if not line.startswith('#'):
                        parts = line.strip().split('\t')
                        new_line = f"{parts[0]}_{parts[1]}_{parts[2]}\t" + "\t".join(parts[3:]) + "\n"
                        out.write(new_line)
    except Exception as e:
        st.error(f"Error processing output file: {str(e)}")
        return False
    
    st.success("Peak counts calculated successfully.")
    return True


def calculate_normalization_factors(input_file, output_file):
    try:
        # Read the input file
        counts = pd.read_csv(input_file, sep='\t', index_col=0)
        
        # Calculate size factors using DESeq2 method
        def estimate_size_factors(counts):
            # Calculate geometric means of rows
            geo_means = gmean(counts, axis=1)
            
            # Calculate ratios of each sample to the geometric mean
            ratios = counts.div(geo_means, axis=0)
            
            # Calculate size factors as median of ratios, excluding zeros and NaNs
            size_factors = ratios.replace(0, np.nan).median(skipna=True)
            
            return size_factors
        
        sf = estimate_size_factors(counts)
        
        # Create DataFrame with size factors and normalizers
        sf_df = pd.DataFrame({
            'sf': sf,
            'normalizer': 1 / sf
        })
        
        # Write the results to the output file
        sf_df.to_csv(output_file, sep='\t', index=True, index_label='sample')
        
        # Display the results
        st.subheader("Normalization Factors")
        st.dataframe(sf_df.style.format({
            'sf': '{:.6f}',
            'normalizer': '{:.6f}'
        }))
        
        st.success("Normalization factors calculated successfully.")
        return True
    except Exception as e:
        st.error(f"Error calculating normalization factors: {str(e)}")
        return False


# File upload for BAM and BAI files
uploaded_files = st.file_uploader("Upload BAM files and their indices (optional)", type=["bam", "bai"], accept_multiple_files=True)

# File upload for peak file
uploaded_peak_file = st.file_uploader("Upload peak file (BED/broadPeak/narrowPeak format)", type=["bed", "broadPeak", "narrowPeak"])

# Genome selection
genome = st.selectbox("Select genome", ["hg38", "mm39", "mm10"])

# Classify uploaded files
bam_files = []
bai_files = []
if uploaded_files:
    for file in uploaded_files:
        if file.name.endswith('.bam'):
            bam_files.append(file)
        elif file.name.endswith('.bai'):
            bai_files.append(file)

# Display table of uploaded files and sample names
if bam_files:
    st.subheader("Uploaded Files")
    file_data = {
        'Sample Name': [os.path.splitext(file.name)[0] for file in bam_files],
        'BAM File': [file.name for file in bam_files],
        'Index File': [
            next((bai_file.name for bai_file in bai_files if bai_file.name == f"{os.path.splitext(bam_file.name)[0]}.bai" or 
                                                              bai_file.name == f"{bam_file.name}.bai"), "Not provided")
            for bam_file in bam_files
        ]
    }
    st.table(pd.DataFrame(file_data))

if st.button("Run Analysis"):
    if not bam_files:
        st.error("Please upload BAM files.")
    elif not uploaded_peak_file:
        st.error("Please upload a peak file.")
    else:
        # Create temporary directory using tmpfs
        with tempfile.TemporaryDirectory(dir="/dev/shm") as temp_dir:
            st.info(f"Using temporary directory: {temp_dir}")
            
            # Save uploaded BAM files
            saved_bam_files = []
            for file in bam_files:
                file_path = os.path.join(temp_dir, file.name)
                with open(file_path, "wb") as f:
                    f.write(file.getbuffer())
                saved_bam_files.append(file_path)
            
            # Save uploaded BAI files
            for file in bai_files:
                file_path = os.path.join(temp_dir, file.name)
                with open(file_path, "wb") as f:
                    f.write(file.getbuffer())
            
            # Save uploaded peak file
            peak_file_path = os.path.join(temp_dir, uploaded_peak_file.name)
            with open(peak_file_path, "wb") as f:
                f.write(uploaded_peak_file.getbuffer())
            
            # Calculate greenlist counts
            greenlist_bed = f"db/blacklist/{genome}_CUTnRUN_greenlist.v1.bed"  # Assume this file exists
            st.info(f"Using greenlist file: {greenlist_bed}")
            if calculate_greenlist_counts(saved_bam_files, greenlist_bed, os.path.join(temp_dir, "greenlist")):
                # Calculate counts for peak file
                if calculate_peak_counts(saved_bam_files, peak_file_path, os.path.join(temp_dir, "peaks")):
                    # Calculate normalization factors
                    if calculate_normalization_factors(os.path.join(temp_dir, "greenlist_glist_quant.tsv"), os.path.join(temp_dir, "glist_sizeFactors.tsv")):
                        # Create a zip file containing the results
                        zip_filename = f"{os.path.splitext(uploaded_peak_file.name)[0]}_results.zip"
                        with zipfile.ZipFile(zip_filename, 'w') as zipf:
                            zipf.write(os.path.join(temp_dir, "peaks_peak_quant.tsv"), "peak_counts.tsv")
                            zipf.write(os.path.join(temp_dir, "greenlist_glist_quant.tsv"), "greenlist_counts.tsv")
                            zipf.write(os.path.join(temp_dir, "glist_sizeFactors.tsv"), "size_factors.tsv")
                        
                        # Offer download
                        with open(zip_filename, "rb") as f:
                            st.download_button(
                                label="Download Results",
                                data=f,
                                file_name=zip_filename,
                                mime="application/zip"
                            )
                        
                        # Clean up
                        os.remove(zip_filename)
                    else:
                        st.error("Error in calculating normalization factors. Please check the logs above.")
                else:
                    st.error("Error in calculating peak counts. Please check the logs above.")
            else:
                st.error("Error in calculating greenlist counts. Please check the logs above.")

        st.success("Analysis completed and temporary files cleaned up.")

st.info("Use SF column of size_factors.tsv as Size Factors.")