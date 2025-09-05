import streamlit as st
import subprocess
import tempfile
import os
import zipfile
import gzip
import shutil
import pandas as pd
import sys
import re
import collections

class Bedgraph_element:
    def __init__(self, chrom, start, end, signal):
        assert(isinstance(chrom, str))
        assert(isinstance(start, int))
        assert(isinstance(end, int))
        assert(isinstance(signal, float) or isinstance(signal, int))
        self.chrom = chrom
        self.start = start
        self.end = end
        self.signal = float(signal)

    def __str__(self):
        return '\t'.join([str(x) for x in [self.chrom, self.start, self.end, self.signal]])

    def len(self):
        return self.end - self.start

    def get_signal(self):
        return self.signal

    def set_signal(self, signal):
        self.signal = signal

def file_gzipped(f):
    return re.search('.gz$', f) is not None

def open_bedgraph(bedgraph_file):
    if file_gzipped(bedgraph_file):
        return gzip.open(bedgraph_file, 'rt')
    else:
        return open(bedgraph_file, 'r')

def parse_line(line):
    line = line.rstrip().split('\t')
    chrom = line[0]
    start = int(line[1])
    end = int(line[2])
    signal = float(line[3])
    return Bedgraph_element(chrom, start, end, signal)

def save_uploaded_file(uploaded_file, destination):
    if uploaded_file.name.endswith('.gz'):
        with gzip.open(uploaded_file, 'rb') as f_in:
            with open(destination, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
    else:
        with open(destination, "wb") as f:
            f.write(uploaded_file.getvalue())

def remove_suffix(text):
    suffixes = ['bedgraph', 'bg']
    lower_text = text.lower()
    for suffix in suffixes:
        if lower_text.endswith(suffix):
            return text[:-len(suffix)]
    return text

def normalize_bedgraph(input_file, output_file, normalization_method, target_value, read_length=None):
    total_signal = 0
    bases = 0

    with open_bedgraph(input_file) as b:
        for line in b:
            element = parse_line(line)
            total_signal += (element.len() * element.get_signal())
            bases += element.len()

    current_avg_signal = float(total_signal) / bases
    st.info(f'Current average signal per base: {current_avg_signal}')

    if normalization_method == 'to_mean_signal':
        normalization_factor = target_value / current_avg_signal
    elif normalization_method == 'to_number_reads':
        if read_length is None or read_length <= 0:
            st.error("Invalid read length. Please provide a positive integer.")
            return False
        number_reads = round(total_signal / read_length)
        st.info(f'Estimated number of reads: {number_reads}')
        normalization_factor = float(target_value) / number_reads
    else:
        st.error('Invalid normalization method')
        return False

    st.info(f'Normalization factor: {normalization_factor}')

    with open_bedgraph(input_file) as infile, open(output_file, 'w') as outfile:
        for line in infile:
            element = parse_line(line)
            element.set_signal(element.get_signal() * normalization_factor)
            outfile.write(str(element) + '\n')

    return True

st.title("Peak calling with SEACR v1.4 for CUT&RUN")

# Multiple target file uploads
target_files = st.file_uploader("##### Target data bedgraph file(s)", type=["bedgraph", "gz"], accept_multiple_files=True)

st.markdown("---")
# Single control file upload
control_file = st.file_uploader("##### Control (IgG) data bedgraph file", type=["bedgraph", "gz"])
st.markdown("***Or***")
# Numeric threshold as an alternative to control file
numeric_threshold = st.number_input("###### Numeric threshold (alternative to control file)", 
    min_value=0.0, max_value=1.0, value=0.01, step=0.01,
    help = "A numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.")
st.markdown("---")

# Add normalization options
normalization_method = st.radio("##### Normalization of bedgraphs", ["None", "To mean signal", "To number of reads"],
    help = 'Normalize bedgraphs. May modestly affect peak calling.'
    )

if normalization_method == "To mean signal":
    target_mean_signal = st.number_input("Target mean signal", min_value=0.0, value=1.0, step=0.1)
elif normalization_method == "To number of reads":
    target_number_reads = st.number_input("Target number of reads", min_value=1, value=10000000, step=1000000)
    read_length = st.number_input("Read length", min_value=1, value=100, step=1)
    st.info("If you're unsure about the read length, you can try common values like 50, 75, or 100.")

# Normalization option for SEACR
seacr_normalization = st.radio("##### SEACR Normalization", ["norm", "non"], help = 
'"norm" denotes normalization of control to target data, "non" skips this behavior. "norm" is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).')

# Mode selection
mode = st.radio("##### Mode", ["relaxed", "stringent"], 
    help = '"relaxed" uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the "relaxed" mode described in the text, whereas "stringent" uses the peak of the curve, and corresponds to "stringent" mode.')

# Peak extension constant
peak_extension = st.number_input("Peak extension constant", min_value=0.0, max_value=1.0, value=0.1, step=0.1, help = 
"A numeric peak extension factor (-e) during peak merging, where all peaks will be extended by e*(mean(peak length)).")

# Remove peaks overlapping IgG peaks
remove_overlapping = st.radio("##### Remove peaks overlapping IgG peaks", ["yes", "no"], help = 
'Peaks overlapping robust IgG peaks get filtered out or not.')

# Run SEACR button
if st.button("Run SEACR"):
    if target_files:
        # Create a temporary directory for SEACR input and output
        with tempfile.TemporaryDirectory(dir="/dev/shm") as tmpdirname:
            # Save control file if provided
            if control_file is not None:
                control_path = os.path.join(tmpdirname, "control.bedgraph")
                save_uploaded_file(control_file, control_path)
                
                # Normalize control file if needed
                if normalization_method != "None":
                    normalized_control_path = os.path.join(tmpdirname, "normalized_control.bedgraph")
                    if normalization_method == "To mean signal":
                        success = normalize_bedgraph(control_path, normalized_control_path, 'to_mean_signal', target_mean_signal)
                    elif normalization_method == "To number of reads":
                        success = normalize_bedgraph(control_path, normalized_control_path, 'to_number_reads', target_number_reads, read_length)
                    
                    if success:
                        control_path = normalized_control_path
                        st.success("Control file normalized successfully.")
                    else:
                        st.error("Failed to normalize control file. Using original file.")
            
            # Process each target file
            for i, target_file in enumerate(target_files):
                output_prefix = remove_suffix(os.path.splitext(target_file.name)[0])
                target_path = os.path.join(tmpdirname, f"target_{i}.bedgraph")
                save_uploaded_file(target_file, target_path)
                
                # Normalize target file if needed
                if normalization_method != "None":
                    normalized_target_path = os.path.join(tmpdirname, f"normalized_target_{i}.bedgraph")
                    if normalization_method == "To mean signal":
                        success = normalize_bedgraph(target_path, normalized_target_path, 'to_mean_signal', target_mean_signal)
                    elif normalization_method == "To number of reads":
                        success = normalize_bedgraph(target_path, normalized_target_path, 'to_number_reads', target_number_reads, read_length)
                    
                    if success:
                        target_path = normalized_target_path
                        st.success(f"Target file {target_file.name} normalized successfully.")
                    else:
                        st.warning(f"Failed to normalize target file {target_file.name}. Using original file.")

                # Prepare command
                cmd = [
                    "bash", "SEACR_1.4.sh",
                    "-b", target_path,
                    "-c", control_path if control_file else str(numeric_threshold),
                    "-n", seacr_normalization,
                    "-m", mode,
                    "-o", os.path.join(tmpdirname, f"{output_prefix}"),
                    "-e", str(peak_extension),
                    "-r", remove_overlapping
                ]
                
                # Run SEACR
                try:
                    result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                    st.success(f"SEACR analysis completed successfully for {target_file.name}!")
                    st.text(result.stdout)
                except subprocess.CalledProcessError as e:
                    st.error(f"Error running SEACR for {target_file.name}: {e.stderr}")
            
            # Create a zip file containing all results
            zip_path = os.path.join(tmpdirname, "seacr_results.zip")
            with zipfile.ZipFile(zip_path, 'w') as zipf:
                for file in os.listdir(tmpdirname):
                    if file.endswith(f".{mode}.bed"):
                        zipf.write(os.path.join(tmpdirname, file), file)
            
            # Provide download link for the zip file
            with open(zip_path, "rb") as f:
                st.download_button(
                    label="Download All Results (ZIP)",
                    data=f,
                    file_name="seacr_results.zip",
                    mime="application/zip"
                )

st.markdown("""
---
Output bed file:

Field 1: Chromosome

Field 2: Start coordinate

Field 3: End coordinate

Field 4: Total signal contained within denoted coordinates

Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates

Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal

Meers, MP, Tenenbaum, D and Henikoff S (2019). Peak calling by sparse enrichment analysis 
for CUT&RUN chromatin profiling. *Epigenetics & Chromatin* 2019 **12:42**.
""")