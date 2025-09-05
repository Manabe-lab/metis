import streamlit as st
from streamlit.runtime.uploaded_file_manager import UploadedFile
import pandas as pd
import subprocess
import os
import zipfile
import io
import sys
import csv
import re
import shutil
import time
from typing import List, Tuple


# /dev/shm ディレクトリのパス
shm_dir = '/dev/shm'

# streamlit_temp ディレクトリのパス
streamlit_temp_dir = os.path.join(shm_dir, 'streamlit_temp')

# streamlit_temp ディレクトリが存在しない場合、作成する
if not os.path.exists(streamlit_temp_dir):
    try:
        os.makedirs(streamlit_temp_dir)
        print(f"Created directory: {streamlit_temp_dir}")
    except OSError as e:
        print(f"Error creating directory {streamlit_temp_dir}: {e}")
else:
    print(f"Directory already exists: {streamlit_temp_dir}")

# 作成したディレクトリのパスを使用する
TMPFS_MOUNT = streamlit_temp_dir

# 新しい関数: 古いファイルを削除
def clean_old_files(directory: str, days: int = 2):
    current_time = time.time()
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            file_age = current_time - os.path.getmtime(file_path)
            if file_age > (days * 24 * 60 * 60):  # days to seconds
                try:
                    os.remove(file_path)
                    st.text(f"Removed old file: {file_path}")
                except Exception as e:
                    st.error(f"Error removing {file_path}: {e}")

# アプリケーションの開始時に古いファイルを削除
clean_old_files(TMPFS_MOUNT)


# ブラックリストファイルのパスを定義
BLACKLIST_FILES = {
    'mm10': {
        'ChIP-seq': 'db/blacklist/mm10-blacklist.v2.bed',
        'CUT&RUN-seq': 'db/blacklist/mm10_CUTnRUN_blacklist.v1.bed',
        'ATAC-seq': 'db/blacklist/mm10_peaks.narrowPeak'
    },
    'mm39': {
        'ChIP-seq': 'db/blacklist/mm39-blacklist.v2.bed',
        'CUT&RUN-seq': 'db/blacklist/mm39_CUTnRUN_blacklist.v1.bed',
         'ATAC-seq': 'db/blacklist/mm39_peaks.narrowPeak'
            },
    'hg38': {
        'ChIP-seq': 'db/blacklist/hg38-blacklist.v2.bed',
        'CUT&RUN-seq': 'db/blacklist/hg38_CUTnRUN_blacklist.v1.bed',
        'ATAC-seq': 'db/blacklist/hg38_peaks.narrowPeak'
    }
}


def run_command(command: str) -> None:
    st.text(f"Running command: {command}")
    try:
        result = subprocess.run(command, shell=True, check=True, capture_output=True, text=True)
        if result.stderr:
            st.warning(f"Command stderr: {result.stderr[:4000]}")
    except subprocess.CalledProcessError as e:
        st.error(f"Command failed with error: {e.stderr[:1000]}")
        raise

def filter_blacklist(file_path: str, genome: str, seq_type: str) -> str:
    if not os.path.exists(file_path):
        st.error(f"Input file not found: {file_path}")
        raise FileNotFoundError(f"Input file not found: {file_path}")

    output_file = os.path.join(TMPFS_MOUNT, os.path.basename(file_path).replace('.bed', '.blacklist-filtered.bed'))

    # Sort the BED file
    sorted_bed = os.path.join(TMPFS_MOUNT, f"sorted_{os.path.basename(file_path)}")
    sort_command = f"sort -k1,1 -k2,2n {file_path} > {sorted_bed}"
    run_command(sort_command)

    # Remove blacklisted regions
    blacklist_file = BLACKLIST_FILES[genome][seq_type]
    if not os.path.exists(blacklist_file):
        st.error(f"Blacklist file not found: {blacklist_file}")
        raise FileNotFoundError(f"Blacklist file not found: {blacklist_file}")

    filter_command = f"bedtools intersect -v -a {sorted_bed} -b {blacklist_file} > {output_file}"
    run_command(filter_command)

    with open(output_file, 'r') as f:
        content = f.read()
        st.text(f"Contents of output file {output_file}:")
        st.text(content[:200])

    if not content:
        st.warning(f"The filtered file {output_file} is empty. Using the original file.")
        shutil.copy(file_path, output_file)
    else:
        st.success(f"Filtered file saved as: {output_file}")

    return output_file

def run_annotate_peaks(file_path: str, genome: str) -> str:
    output_file = os.path.join(TMPFS_MOUNT, os.path.basename(file_path).replace('.bed', '.annotated.bed'))
    command = f"annotatePeaks.pl {file_path} {genome} > {output_file}"
    run_command(command)
    return output_file

def filter_annotate_file(file_path: str, search_column: str, search_words: List[str]) -> str:
    output_file = os.path.join(TMPFS_MOUNT, os.path.basename(file_path).replace('.annotated.bed', '.filtered.annotated.bed'))
    
    try:
        with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
            reader = csv.reader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t', quoting=csv.QUOTE_NONE)
            
            header = next(reader)
            writer.writerow(header)

            if search_column not in header:
                st.error(f"Column '{search_column}' not found in the file header.")
                return file_path

            lookup_index = header.index(search_column)

            for row in reader:
                if len(row) > lookup_index:
                    annotation = str(row[lookup_index])
                    if not any(word.lower() in annotation.lower() for word in search_words):
                        writer.writerow(row)
                else:
                    st.warning(f"Skipping a row with insufficient columns: {row}")

    except Exception as e:
        st.error(f"An error occurred while filtering the file {file_path}: {str(e)}")
        return file_path

    return output_file

def bed_to_pos(bed_file: str) -> str:
    pos_file = bed_file.replace('.bed', '.txt')
    command = f"bed2pos.pl {bed_file} > {pos_file}"
    run_command(command)
    return pos_file

@st.cache_data
def get_annotation_options(file_path: str) -> List[str]:
    df = pd.read_csv(file_path, sep='\t')
    if 'Detailed Annotation' not in df.columns:
        st.error(f"'Detailed Annotation' column not found in {file_path}")
        return []
    annotations = df['Detailed Annotation'].astype(str).unique().tolist()
    cleaned_annotations = []
    for ann in annotations:
        if isinstance(ann, str):
            if '|' in ann:
                parts = ann.split('|')
                if len(parts) >= 3:
                    cleaned_annotations.append(parts[1])
            elif '(' in ann:
                cleaned_annotations.append(ann.split(' (')[0])
            elif ann.startswith('CpG-'):
                cleaned_annotations.append('CpG')
            else:
                cleaned_annotations.append(ann)
    return sorted(set(cleaned_annotations))

@st.cache_data
def process_files(uploaded_files: List[UploadedFile], genome: str, seq_type: str, process_level: int) -> List[str]:
    processed_files = []
    for uploaded_file in uploaded_files:
        original_file_name = uploaded_file.name
        original_file_path = os.path.join(TMPFS_MOUNT, original_file_name)

        with open(original_file_path, "wb") as f:
            f.write(uploaded_file.getbuffer())
        st.text(f"Saved original file: {original_file_path}")

        if original_file_name.endswith('.txt'):
            bed_file = os.path.join(TMPFS_MOUNT, original_file_name.replace('.txt', '.bed'))
            pos2bed_command = f"pos2bed.pl {original_file_path} > {bed_file}"
            run_command(pos2bed_command)
            input_file = bed_file
        elif original_file_name.endswith('Peak'):
            input_file = original_file_path + ".bed"
            os.rename(original_file_path, input_file)
        else:
            input_file = original_file_path

        # レベル1: アノテーション
        if process_level >= 1:
            annotated_file = run_annotate_peaks(input_file, genome)
            if os.path.getsize(annotated_file) > 0:
                processed_files.append(annotated_file)
                if original_file_name.endswith('.txt'):
                    processed_files.append(bed_to_pos(annotated_file))
            else:
                st.warning(f"The annotated file {annotated_file} is empty.")

        # レベル2: ブラックリストフィルタリング
        if process_level >= 2:
            blacklist_filtered_file = filter_blacklist(input_file, genome, seq_type)
            if os.path.getsize(blacklist_filtered_file) > 0:
                processed_files.append(blacklist_filtered_file)
                if original_file_name.endswith('.txt'):
                    processed_files.append(bed_to_pos(blacklist_filtered_file))
                
                blacklist_filtered_annotated_file = run_annotate_peaks(blacklist_filtered_file, genome)
                if os.path.getsize(blacklist_filtered_annotated_file) > 0:
                    processed_files.append(blacklist_filtered_annotated_file)
                    if original_file_name.endswith('.txt'):
                        processed_files.append(bed_to_pos(blacklist_filtered_annotated_file))
                else:
                    st.warning(f"The annotated file {blacklist_filtered_annotated_file} is empty.")
            else:
                st.warning(f"The blacklist filtered file {blacklist_filtered_file} is empty.")

        # レベル3は処理レベル2のファイルを使用するため、ここでは特に処理を行いません。

    return processed_files




def is_index_row(row: List[str]) -> bool:
    if len(row) < 4:  # BEDファイルは少なくとも3列（染色体、開始、終了）必要
        return False
    if row[0].startswith('PeakID'):
        return True
    else:
        return False


def remove_index_from_bed(file_path: str) -> str:
    output_file = os.path.join(TMPFS_MOUNT, f"no_header.{os.path.basename(file_path)}")
    with open(file_path, 'r') as infile, open(output_file, 'w') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t', quoting=csv.QUOTE_MINIMAL)
        
        first_row = next(reader, None)
        if first_row:
            if not is_index_row(first_row):
                # 最初の行がインデックス行でない場合、書き込む
                writer.writerow(first_row)
            
            # 残りの行を処理
            for row in reader:
                if not is_index_row(row):
                    writer.writerow(row)
    
    return output_file

def create_zip(files: List[str]) -> io.BytesIO:
    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "a", zipfile.ZIP_DEFLATED, False) as zip_file:
        for file in files:
            if file.endswith('.bed'):
                no_index_file = remove_index_from_bed(file)
                zip_file.write(file, os.path.basename(file))
                zip_file.write(no_index_file, os.path.basename(no_index_file))
                os.remove(no_index_file)  # 一時ファイルを削除
            else:
                zip_file.write(file, os.path.basename(file))
    return zip_buffer


def filter_annotations(files: List[str], annotations_to_remove: List[str]) -> List[str]:
    filtered_files = []
    for file in files:
        if not os.path.exists(file):
            st.error(f"File not found: {file}")
            continue
        if os.path.getsize(file) > 0:
            output_file = filter_annotate_file(file, "Detailed Annotation", annotations_to_remove)
            filtered_files.append(output_file)
            if file.endswith('.txt'):
                filtered_files.append(bed_to_pos(output_file))
        else:
            st.warning(f"Skipping empty file: {file}")
    return filtered_files

# Streamlitアプリケーションのメイン部分
st.title('Bed/Peak file annotator and filter')
st.markdown('#### peak annotation, filtering of blacklisted peaks and specific class of peaks (e.g., simple repeats)')

if 'processed_files' not in st.session_state:
    st.session_state.processed_files = []
if 'annotation_options' not in st.session_state:
    st.session_state.annotation_options = []
if 'selected_annotations' not in st.session_state:
    st.session_state.selected_annotations = []

def update_selected_annotations():
    st.session_state.selected_annotations = st.session_state.multiselect_annotations

if st.button("Clear cache for new uploads?"):
    st.cache_data.clear()
    st.session_state.processed_files = []
    st.session_state.annotation_options = []
    st.session_state.selected_annotations = []
    # Clean up temporary files
    if st.session_state.processed_files:
        for file in st.session_state.processed_files:
            if os.path.exists(file):
                os.remove(file)
        st.session_state.processed_files = []
        st.success('Cache and temporary files cleaned up')


uploaded_files = st.file_uploader("Choose bed, macs' Peak, or Homer peak txt files", accept_multiple_files=True)
st.write("Homer peak file must be .txt")

if uploaded_files or st.session_state.annotation_options:
    st.write(f"Uploaded files: {[file.name for file in uploaded_files]}")
    
    genome = st.selectbox('Select Genome', ['mm10', 'mm39', 'hg38'])
    seq_type = st.radio("Select sequencing type", ('ChIP-seq', 'CUT&RUN-seq', 'ATAC-seq'), index=1)

    process_level = st.radio("Select processing level:",
                             ["1. Annotation only",
                              "2. Filter blacklist",
                              "3. Further filtering"],
                             index=0)

    if st.button('Process Files') or st.session_state.annotation_options:
        try:
            st.session_state.processed_files = process_files(uploaded_files, genome, seq_type, int(process_level[0]))
            
            st.success('Processing complete!')

            if int(process_level[0]) == 3:
                if st.session_state.processed_files:
                    last_processed_file = [f for f in st.session_state.processed_files if f.endswith('.annotated.bed')][-1]
                    
                    if not st.session_state.annotation_options:
                        st.session_state.annotation_options = get_annotation_options(last_processed_file)
                    
                    st.multiselect(
                        '##### Select annotations to filter out:',
                        st.session_state.annotation_options,
                        default=st.session_state.selected_annotations,
                        key='multiselect_annotations',
                        on_change=update_selected_annotations
                    )

                    if st.button('Apply Annotation Filter'):
                        try:
                            with st.spinner('Applying annotation filter...'):
                                annotated_files = [f for f in st.session_state.processed_files if f.endswith('.annotated.bed')]
                                filtered_files = filter_annotations(annotated_files, st.session_state.selected_annotations)
                                st.session_state.processed_files.extend(filtered_files)

                            st.success('Annotation filtering complete!')
                            
                            for file in filtered_files:
                                st.text(f"Contents of filtered file ({file}):")
                                with open(file, 'r') as f:
                                    st.text(f.read()[:200])
                            
                            if filtered_files:
                                zip_buffer = create_zip(filtered_files)
                                st.download_button(
                                    label="Download filtered files",
                                    data=zip_buffer.getvalue(),
                                    file_name="filtered_files.zip",
                                    mime="application/zip"
                                )
                            else:
                                st.warning("No files were processed after filtering.")
                        except Exception as e:
                            st.error(f"An error occurred during filtering: {str(e)}")
                            st.text(f"Error details: {sys.exc_info()}")
                else:
                    st.warning("No files were processed. Please check the previous steps for any errors.")
                st.markdown("---")
            
            # Prepare download button for all processed files
            if st.session_state.processed_files:
                zip_buffer = create_zip(st.session_state.processed_files)
                st.download_button(
                    label="Download all processed files",
                    data=zip_buffer.getvalue(),
                    file_name="all_processed_files.zip",
                    mime="application/zip"
                )
            else:
                st.warning("No files were processed. All input files were empty.")
        except Exception as e:
            st.error(f"An error occurred during processing: {str(e)}")
            st.text(f"Error details: {sys.exc_info()}")

else:
    st.error('Please upload files first.')

st.markdown('---')
# Clean up temporary files
if st.session_state.processed_files:
    if st.button('Clean up temporary files'):
        for file in st.session_state.processed_files:
            if os.path.exists(file):
                os.remove(file)
        st.session_state.processed_files = []
        st.success('Temporary files cleaned up')

st.text(f"Files in {TMPFS_MOUNT}: {os.listdir(TMPFS_MOUNT)}")