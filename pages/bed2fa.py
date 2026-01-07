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
    BEDFileoffuo-matothecheckoutdorelnum
    """
    with open(file_path, 'r') as f:
        # heda-thesukipu
        line = f.readline()
        while line.startswith('#') or is_header_line(line):
            line = f.readline()
        
        fields = line.strip().split('\t')
        
        # MACSfuo-matoofspecfeaturetheConfirm
        if len(fields) >= 4:
            try:
                # passnormalofBEDfuo-mato（chr, start, end）ofplacematch
                if fields[0].startswith('chr') and fields[1].isdigit() and fields[2].isdigit():
                    return 'standard'
                # MACSfuo-matoofplacematch（name, chr, start, end）
                elif fields[1].startswith('chr') and fields[2].isdigit() and fields[3].isdigit():
                    return 'macs'
            except:
                pass
    return 'unknown'

@st.cache_data
def convert_to_standard_bed(input_path, output_path, format_type):
    """
    verschiedene BEDfuo-matothemarklevelalnaBEDfuo-matotochangechangedorelnum
    """
    if format_type == 'macs':
        # MACSfuo-matothemarklevelBEDfuo-matotochangechange
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                if line.startswith('#') or is_header_line(line):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    # chr, start, endofmitheextractoutshitewritekioutshi
                    bed_line = f"{fields[1]}\t{fields[2]}\t{fields[3]}\n"
                    fout.write(bed_line)
        return True
    elif format_type == 'standard':
        # alreadytomarklevelfuo-matoofplacematchissoofmamakopi-
        shutil.copy(input_path, output_path)
        return True
    return False

@st.cache_data
def is_header_line(line):
    """
    rowisheda-（karamuname）whetherthejudgesetdorelnum
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
    BEDFiletheFastaFiletochangechangedorelnum
    """
    bed = BedTool(bed_file)
    bed.sequence(fi=f'db/genome/{genome}.fa', fo=output_file)

def show_file_preview(file_path, format_type=None):
    """
    Fileofpurebiyu-theDisplaydorelnum
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()[:10]  # mostfirstof10rowofmiDisplay
            
        st.write("File preview:")
        for i, line in enumerate(lines, 1):
            if i == 1 and is_header_line(line):
                st.markdown(f"**📋 heda-row:** `{line.strip()}`")
            else:
                st.text(f"row {i}: {line.strip()}")
                
        if format_type:
            st.info(f"Detected format: {format_type}")
            
        if len(lines) == 10:
            st.text("...")
    except Exception as e:
        st.error(f"purebiyu-ofDisplaymidtoErrorisoccurgenshimashita: {str(e)}")

def main():
    st.title('BED to FASTA Converter')

    # genomeofSelect
    genome = st.selectbox(
        'Genome:',
        ('mm10', 'mm39', 'hg38')
    )

    # multinumFileUpload
    uploaded_files = st.file_uploader(
        "Upload bed files", 
        type=['bed', 'broadPeak', 'narrowPeak'],
        accept_multiple_files=True
    )
    st.write("can handle MACS outputs")

    if uploaded_files:
        # progprogstatestatetheDisplaydopuroguresuba-
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # changechangedonemiFilethe格storedomemoriupofbafua
        zip_buffer = io.BytesIO()
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # ZIPFileofmakebecome
                with ZipFile(zip_buffer, 'w') as zip_file:
                    # eachFileofprocproc
                    for i, uploaded_file in enumerate(uploaded_files):
                        status_text.text(f'Processing: {uploaded_file.name}')
                        
                        # progprogstatestateoffurthernew
                        progress = (i + 1) / len(uploaded_files)
                        progress_bar.progress(progress)

                        # InputFileofSave
                        input_path = os.path.join(temp_dir, uploaded_file.name)
                        with open(input_path, 'wb') as f:
                            f.write(uploaded_file.getbuffer())

                        # fuo-matoofcheckout
                        format_type = detect_bed_format(input_path)
                        st.info(f"Format type of {uploaded_file.name}: {format_type}")

                        # marklevelBEDfuo-matoheofchangechange
                        standard_bed_path = os.path.join(temp_dir, f"standard_{uploaded_file.name}")
                        if not convert_to_standard_bed(input_path, standard_bed_path, format_type):
                            st.error(f"{uploaded_file.name}: sapo-tosareteinotfuo-matowithsu。")
                            show_file_preview(input_path, format_type)
                            continue

                        # OutputFilenameofSettings（.bedthe.fatochangefurther）
                        output_filename = os.path.splitext(uploaded_file.name)[0] + '.fa'
                        output_path = os.path.join(temp_dir, output_filename)

                        try:
                            # marklevelfuo-matoofBEDFiletheFastatochangechange
                            convert_bed_to_fasta(standard_bed_path, genome, output_path)
                            
                            # changechangeshitaFiletheZIPtoaddadd
                            zip_file.write(output_path, output_filename)
                            
                        except Exception as e:
                            st.error(f"""
                            Errorisoccurgenshimashita:
                            
                            File: {uploaded_file.name}
                            {str(e)}
                            """)
                            
                            # Filepurebiyu-theDisplay
                            show_file_preview(standard_bed_path, format_type)
                            continue
                
                # allteofprocprocisCompleteshitaraDisplay
                status_text.text('Done！')
                
                # ZIPFileofDownloadbotan
                zip_buffer.seek(0)
                st.download_button(
                    label="Download fasta files",
                    data=zip_buffer,
                    file_name="converted_fasta_files.zip",
                    mime="application/zip"
                )

                # procprocResultofsamari-Display
                st.write("### Results")
                for file in uploaded_files:
                    st.write(f"- {file.name} → {os.path.splitext(file.name)[0]}.fa")

        except Exception as e:
            st.error(f"advperiodsenuErrorisoccurgenshimashita: {str(e)}")

if __name__ == '__main__':
    main()