import streamlit as st
import subprocess
import os
import tempfile
import zipfile
import io
import sys

def run_command(command, description):
    st.subheader(description)
    cmd_str = " ".join(command)
    st.code(cmd_str, language="bash")
    try:
        result = subprocess.run(command, capture_output=True, text=True, check=True)
        st.subheader("MACS3 Output:")
        if result.stdout:
            st.text(result.stdout)
        if result.stderr:
            st.text("Log information:")
            st.text(result.stderr)
        return result.stdout, result.stderr
    except subprocess.CalledProcessError as e:
        st.error(f"Error occurred while running {description}")
        st.text("Error output:")
        st.text(e.stderr)
        st.text("Please check your MACS3 installation and Python environment.")
        st.text("Recommended Python version for MACS3: 3.7 - 3.10")
        st.text(f"Current Python version: {sys.version}")
        return None, e.stderr

def create_zip(directory):
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
        for root, _, files in os.walk(directory):
            for file in files:
                if not file.endswith('.bam'):  # BAMファイルを除外
                    file_path = os.path.join(root, file)
                    arcname = os.path.relpath(file_path, directory)
                    zip_file.write(file_path, arcname)
    return zip_buffer

st.title("MACS3 Peak Calling")

# Main panel
col1, col2 = st.columns(2)

species = st.selectbox("Select species", ["Mouse", "Human"])
genome_size = "mm" if species == "Mouse" else "hs"

target_files = st.file_uploader("#### Upload target paired-end bam file(s)", type="bam", accept_multiple_files=True)
control_file = st.file_uploader("#### Upload control paired-end bam file", type="bam")
is_broad = st.checkbox("##### Enable broad peak calling (--broad)", help="Use this for broad peaks like H3K27ac")

st.write("---")

# Sidebar for parameter settings
st.sidebar.header("MACS3 callpeak Parameters")

# Set recommended parameters based on broad peak calling
if is_broad:
    default_qvalue = 0.1
    default_broad_cutoff = 0.1
else:
    default_qvalue = 0.05
    default_broad_cutoff = None

# Common parameters
name = st.sidebar.text_input("Experiment name (-n)", value="MACS3_output")
qvalue = st.sidebar.number_input("q-value cutoff (-q)", value=default_qvalue, format="%.3f", step=0.001)
pvalue = st.sidebar.number_input("p-value cutoff (-p)", value=None, format="%.3f", step=0.001)
#gsize = st.sidebar.text_input("Effective genome size (-g)", value=genome_size)

# Advanced parameters
st.sidebar.subheader("Advanced Parameters")
tsize = st.sidebar.number_input("Tag size (-s)", value=None)
bw = st.sidebar.number_input("Band width (--bw)", value=300)
mfold_lower = st.sidebar.number_input("Lower mfold bound (--mfold)", value=5)
mfold_upper = st.sidebar.number_input("Upper mfold bound (--mfold)", value=50)
nolambda = st.sidebar.checkbox("Use fixed background lambda (--nolambda)")
slocal = st.sidebar.number_input("Small local region (-slocal)", value=1000)
llocal = st.sidebar.number_input("Large local region (-llocal)", value=10000)

# Parameters specific to broad peak calling
if is_broad:
    broad_cutoff = st.sidebar.number_input("Broad cutoff (--broad-cutoff)", value=default_broad_cutoff, format="%.3f", step=0.001)

# Other options
keep_dup = st.sidebar.selectbox("How to handle duplicates (--keep-dup)", ["1", "auto", "all"])
call_summits = st.sidebar.checkbox("Call summits (--call-summits)")
cutoff_analysis = st.sidebar.checkbox("Perform cutoff analysis (--cutoff-analysis)")

if st.button("Run MACS3 callpeak"):
    if target_files and control_file:
        with tempfile.TemporaryDirectory() as temp_dir:
            # Save uploaded files
            target_paths = []
            for file in target_files:
                path = os.path.join(temp_dir, file.name)
                with open(path, "wb") as f:
                    f.write(file.getbuffer())
                target_paths.append(path)
            
            control_path = os.path.join(temp_dir, control_file.name)
            with open(control_path, "wb") as f:
                f.write(control_file.getbuffer())
            
            # Get the name of the first target file without extension
            first_target_name = os.path.splitext(os.path.basename(target_paths[0]))[0]
            
            # Update the output name to include the first target file name
            output_name = f"{first_target_name}_{name}"
            
            # Construct MACS3 callpeak command
            cmd = [
                "macs3", "callpeak",
                "-t", *target_paths,
                "-c", control_path,
                "-f", "BAMPE",
                "-g", gsize,
                "-n", output_name,
                "-q", str(qvalue),
                "--outdir", temp_dir
            ]

            if is_broad:
                cmd.extend(["--broad", "--broad-cutoff", str(broad_cutoff)])
            if pvalue:
                cmd.extend(["-p", str(pvalue)])
            if tsize:
                cmd.extend(["-s", str(tsize)])
            if nolambda:
                cmd.append("--nolambda")
            cmd.extend(["--bw", str(bw)])
            cmd.extend(["--mfold", str(mfold_lower), str(mfold_upper)])
            cmd.extend(["--slocal", str(slocal)])
            cmd.extend(["--llocal", str(llocal)])
            cmd.extend(["--keep-dup", keep_dup])
            if call_summits and not is_broad:
                cmd.append("--call-summits")
            if cutoff_analysis:
                cmd.append("--cutoff-analysis")

            # Run MACS3 callpeak
            stdout, stderr = run_command(cmd, "Running MACS3 callpeak")

            if stderr and "Done!" in stderr:  # MACS3の実行が成功した場合
                # Create ZIP file for download (BAMファイルは除外される)
                st.subheader("Download Results")
                zip_buffer = create_zip(temp_dir)
                st.download_button(
                    label="Download All Results (ZIP)",
                    data=zip_buffer.getvalue(),
                    file_name=f"{output_name}_macs3_results.zip",
                    mime="application/zip"
                )
                
                # Display summary of results
                st.subheader("Summary of Results")
                peak_file = os.path.join(temp_dir, f"{output_name}_peaks.xls")
                if os.path.exists(peak_file):
                    with open(peak_file, 'r') as f:
                        peak_count = sum(1 for line in f if not line.startswith('#')) - 1  # Subtract header
                    st.write(f"Number of peaks detected: {peak_count}")
                else:
                    st.write("Peak file not found. Please check the output.")
            else:
                st.warning("MACS3 callpeak may have failed. Please check the output above and your settings.")

    else:
        st.error("Please upload both target and control BAMPE files.")