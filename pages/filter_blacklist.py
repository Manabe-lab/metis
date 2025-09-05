import streamlit as st
import pandas as pd
import os
import subprocess
import tempfile
import zipfile
import io
from concurrent.futures import ThreadPoolExecutor, as_completed

def sort_bed(input_file, output_file):
    command = f'sort -k1,1 -k2,2n {input_file} > {output_file}'
    subprocess.run(command, shell=True, check=True)

def remove_blacklist(input_file, blacklist_file, output_file):
    command = f'bedtools intersect -v -a {input_file} -b {blacklist_file} > {output_file}'
    subprocess.run(command, shell=True, check=True)

def process_file(file, genome, blacklist_file):
    with tempfile.TemporaryDirectory() as tmpdir:
        try:
            input_path = os.path.join(tmpdir, file.name)
            with open(input_path, 'wb') as f:
                f.write(file.getvalue())
            
            root, ext = os.path.splitext(file.name)
            bed_filename = f"{root}.filtered.bed"
            peak_filename = f"{root}.filtered.txt"
            
            if ext == '.txt':
                # Convert peak to BED
                df = pd.read_csv(input_path, sep='\t', comment='#', header=None)
                if len(df.columns) >= 7:  # Not GroSeq
                    bed_df = pd.DataFrame({
                        'chr': df[1],
                        'start': df[2].astype(int) - 1,
                        'end': df[3],
                        'name': df[0],
                        'score': df[7],
                        'strand': df[4],
                        'tag_count': df[5]
                    })
                else:  # GroSeq
                    bed_df = pd.DataFrame({
                        'chr': df[1],
                        'start': df[2].astype(int) - 1,
                        'end': df[3],
                        'name': df[0],
                        'score': df[5],
                        'strand': df[4]
                    })
                bed_path = os.path.join(tmpdir, f"{root}.bed")
                bed_df.to_csv(bed_path, sep='\t', header=False, index=False)
            else:
                bed_path = input_path
            
            # Sort BED file
            sorted_bed_path = os.path.join(tmpdir, f"{root}.sorted.bed")
            sort_bed(bed_path, sorted_bed_path)

            # Remove blacklist regions
            if blacklist_file == 'ATAC':
                blacklist_file = f"{genome}_peaks.narrowPeak"
            elif blacklist_file == 'CUT&RUN':
                blacklist_file = f"{genome}_CUTnRUN_blacklist.v1.bed"
            elif blacklist_file == 'ENCODEv1':
                blacklist_file = f"{genome}-blacklist.v1.bed"
            elif blacklist_file == 'ENCODEv2':
                blacklist_file = f"{genome}-blacklist.v2.bed"
            elif blacklist_file == 'ATAC+ENCODEv2':
                blacklist_file = f"{genome}.ENCODE.ATAC.bed"
            
            blacklist_path = f"./db/blacklist/{blacklist_file}"
            filtered_bed_path = os.path.join(tmpdir, bed_filename)
            remove_blacklist(sorted_bed_path, blacklist_path, filtered_bed_path)
            
            # Convert filtered BED back to peak if necessary
            if ext == '.txt':
                filtered_df = pd.read_csv(filtered_bed_path, sep='\t', header=None)
                peak_df = pd.DataFrame({
                    'PeakID': filtered_df[3],
                    'chr': filtered_df[0],
                    'start': filtered_df[1].astype(int) + 1,
                    'end': filtered_df[2],
                    'strand': filtered_df[5],
                    'findPeaks Score': filtered_df[4]
                })
                if len(filtered_df.columns) > 6:
                    peak_df['Normalized Tag Count'] = filtered_df[6]
                
                output_path = os.path.join(tmpdir, peak_filename)
                peak_df.to_csv(output_path, sep='\t', index=False)
            else:
                output_path = filtered_bed_path

            # Check if output file exists
            if not os.path.exists(output_path):
                raise FileNotFoundError(f"Output file was not created: {output_path}")

            # Copy file outside of temporary directory
            permanent_output_path = os.path.join(os.getcwd(), os.path.basename(output_path))
            with open(output_path, 'rb') as src, open(permanent_output_path, 'wb') as dst:
                dst.write(src.read())

            return permanent_output_path, 'txt' if ext == '.txt' else 'bed'
        except Exception as e:
            st.error(f"Error processing file {file.name}: {str(e)}")
            return None, None

st.title('Blacklist Filter App')
st.write("ChIP-seq: ENCODE")
st.write("ATAC: https://github.com/buenrostrolab/mitoblacklist/tree/master/peaks")
st.write("maxATAC: Genome-scale transcription-factor binding prediction from ATAC-seq with deep neural networks. use included (i) blacklisted regions from ENCODE data, (ii) centromeres, telomeres, and annotated gaps available from UCSC table browser for hg38, (iii) regions ≥1kb with ≥ 90% sequence identity to chrM, and (iv) regions with low mappability on chr21.")
st.write('CUT&RUN: The CUT&RUN greenlist: genomic regions of consistent noise are effective normalizing factors for quantitative epigenome mapping. Bioinformatics 25:bbad538')

genome = st.radio('Select genome', ['hg38', 'mm10', 'mm39'], index=1)
blacklist_file = st.radio('Select blacklist', ['ENCODEv1', 'ENCODEv2','CUT&RUN', 'ATAC', 'ATAC+ENCODEv2'], index=1)

uploaded_files = st.file_uploader("Choose peak/bed files", accept_multiple_files=True, type=['txt','bed'])

# リストをグローバルスコープで定義
processed_files = []

if uploaded_files and st.button('Process Files'):
    progress_bar = st.progress(0)
    status_text = st.empty()
    
    processed_files = []
    with ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_file, file, genome, blacklist_file) for file in uploaded_files]
        for i, future in enumerate(as_completed(futures)):
            processed_file, ext = future.result()
            if processed_file:
                processed_files.append((processed_file, ext))
            progress = (i + 1) / len(uploaded_files)
            progress_bar.progress(progress)
            status_text.text(f'Processed {i+1}/{len(uploaded_files)} files')
    
    if processed_files:
        # Create ZIP file
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, 'w', zipfile.ZIP_DEFLATED) as zip_file:
            for file_path, ext in processed_files:
                if os.path.exists(file_path):
                    file_name = os.path.basename(file_path)
                    with open(file_path, 'rb') as f:
                        zip_file.writestr(file_name, f.read())
                else:
                    st.warning(f"File not found: {file_path}")
        
        # Offer ZIP file for download
        st.download_button(
            label="Download processed files",
            data=zip_buffer.getvalue(),
            file_name="processed_files.zip",
            mime="application/zip"
        )
    else:
        st.error("No files were successfully processed.")

# Clean up temporary files after processing
for file_path, _ in processed_files:
    if os.path.exists(file_path):
        os.remove(file_path)
