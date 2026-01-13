import streamlit as st
import pandas as pd
import os
from pybedtools import BedTool
import tempfile
import shutil
import io
from zipfile import ZipFile

@st.cache_data
def detect_bed_format(file_path):
    """
    Function to detect the format of a BED file
    """
    with open(file_path, 'r') as f:
        # Skip header
        line = f.readline()
        while line.startswith('#') or is_header_line(line):
            line = f.readline()

        fields = line.strip().split('\t')

        # Check characteristics of MACS format
        if len(fields) >= 4:
            try:
                # Standard BED format (chr, start, end)
                if fields[0].startswith('chr') and fields[1].isdigit() and fields[2].isdigit():
                    return 'standard'
                # MACS format (name, chr, start, end)
                elif fields[1].startswith('chr') and fields[2].isdigit() and fields[3].isdigit():
                    return 'macs'
            except:
                pass
    return 'unknown'

@st.cache_data
def convert_to_standard_bed(input_path, output_path, format_type):
    """
    Function to convert various BED formats to standard BED format
    """
    if format_type == 'macs':
        # Convert MACS format to standard BED format
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                if line.startswith('#') or is_header_line(line):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    # Extract and write only chr, start, end
                    bed_line = f"{fields[1]}\t{fields[2]}\t{fields[3]}\n"
                    fout.write(bed_line)
        return True
    elif format_type == 'standard':
        # If already in standard format, copy as is
        shutil.copy(input_path, output_path)
        return True
    return False

@st.cache_data
def is_header_line(line):
    """
    Function to determine if a line is a header (column names)
    """
    fields = line.strip().split('\t')
    if len(fields) >= 3:
        common_headers = ['chr', 'chrom', 'chromosome', 'start', 'end', 'begin', 'stop', 'name', 'peak']
        first_fields_lower = [field.lower() for field in fields[:3]]
        return any(header in first_fields_lower for header in common_headers)
    return False

@st.cache_data
def convert_bed_to_fasta(bed_file, genome, output_file):
    """
    Function to convert a BED file to a FASTA file
    """
    bed = BedTool(bed_file)
    bed.sequence(fi=f'db/genome/{genome}.fa', fo=output_file)

def show_file_preview(file_path, format_type=None):
    """
    Function to display a file preview
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()[:10]  # Display only the first 10 lines

        st.write("File preview:")
        for i, line in enumerate(lines, 1):
            if i == 1 and is_header_line(line):
                st.markdown(f"**Header line:** `{line.strip()}`")
            else:
                st.text(f"Line {i}: {line.strip()}")
                
        if format_type:
            st.info(f"Detected format: {format_type}")
            
        if len(lines) == 10:
            st.text("...")
    except Exception as e:
        st.error(f"An error occurred while displaying the preview: {str(e)}")

def main():
    st.title('BED to FASTA Converter')

    # Genome selection
    genome = st.selectbox(
        'Genome:',
        ('mm10', 'mm39', 'hg38')
    )

    # Multiple file upload
    uploaded_files = st.file_uploader(
        "Upload bed files", 
        type=['bed', 'broadPeak', 'narrowPeak'],
        accept_multiple_files=True
    )
    st.write("can handle MACS outputs")

    if uploaded_files:
        # Progress bar to display processing status
        progress_bar = st.progress(0)
        status_text = st.empty()

        # In-memory buffer to store converted files
        zip_buffer = io.BytesIO()

        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # Create ZIP file
                with ZipFile(zip_buffer, 'w') as zip_file:
                    # Process each file
                    for i, uploaded_file in enumerate(uploaded_files):
                        status_text.text(f'Processing: {uploaded_file.name}')

                        # Update progress status
                        progress = (i + 1) / len(uploaded_files)
                        progress_bar.progress(progress)

                        # Save input file
                        input_path = os.path.join(temp_dir, uploaded_file.name)
                        with open(input_path, 'wb') as f:
                            f.write(uploaded_file.getbuffer())

                        # Detect format
                        format_type = detect_bed_format(input_path)
                        st.info(f"Format type of {uploaded_file.name}: {format_type}")

                        # Convert to standard BED format
                        standard_bed_path = os.path.join(temp_dir, f"standard_{uploaded_file.name}")
                        if not convert_to_standard_bed(input_path, standard_bed_path, format_type):
                            st.error(f"{uploaded_file.name}: Unsupported format.")
                            show_file_preview(input_path, format_type)
                            continue

                        # Set output filename (change .bed to .fa)
                        output_filename = os.path.splitext(uploaded_file.name)[0] + '.fa'
                        output_path = os.path.join(temp_dir, output_filename)

                        try:
                            # Convert standard format BED file to FASTA
                            convert_bed_to_fasta(standard_bed_path, genome, output_path)

                            # Add converted file to ZIP
                            zip_file.write(output_path, output_filename)
                            
                        except Exception as e:
                            st.error(f"""
                            An error occurred:

                            File: {uploaded_file.name}
                            {str(e)}
                            """)

                            # Display file preview
                            show_file_preview(standard_bed_path, format_type)
                            continue

                # Display when all processing is complete
                status_text.text('Done!')

                # ZIP file download button
                zip_buffer.seek(0)
                st.download_button(
                    label="Download fasta files",
                    data=zip_buffer,
                    file_name="converted_fasta_files.zip",
                    mime="application/zip"
                )

                # Display processing results summary
                st.write("### Results")
                for file in uploaded_files:
                    st.write(f"- {file.name} → {os.path.splitext(file.name)[0]}.fa")

        except Exception as e:
            st.error(f"An unexpected error occurred: {str(e)}")

if __name__ == '__main__':
    main()