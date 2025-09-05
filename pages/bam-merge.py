import streamlit as st
import os
import tempfile
import subprocess
import re
import multiprocessing
import zipfile
from functools import reduce

ROOT_DIR = '/home/lab/sftp-data/METIS_data'

def get_directory_structure(rootdir):
    dir_structure = {}
    rootdir = rootdir.rstrip(os.sep)
    start = rootdir.rfind(os.sep) + 1
    for path, dirs, files in os.walk(rootdir):
        folders = path[start:].split(os.sep)
        subdir = dict.fromkeys([f for f in files if f.endswith('.bam')])
        parent = reduce(dict.get, folders[:-1], dir_structure)
        parent[folders[-1]] = subdir
    return dir_structure

def render_directory_tree(directory, path='', depth=0):
    for key, value in directory.items():
        full_path = os.path.join(path, key)
        indent = '  ' * depth
        if isinstance(value, dict):
            st.markdown(f"{indent}📁 **{key}**")
            if st.checkbox(f"Open {key}", key=f"open_{full_path}"):
                render_directory_tree(value, full_path, depth + 1)
        elif key.endswith('.bam'):
            st.checkbox(f"{indent}📄 {key}", key=full_path)

def get_selected_files(directory, path='', selected_files=None):
    if selected_files is None:
        selected_files = []
    for key, value in directory.items():
        full_path = os.path.join(path, key)
        if isinstance(value, dict):
            get_selected_files(value, full_path, selected_files)
        elif key.endswith('.bam') and st.session_state.get(full_path, False):
            selected_files.append(full_path)
    return selected_files

@st.cache_data
def get_cached_directory_structure(rootdir):
    return get_directory_structure(rootdir)

def merge_and_sort_bams(bam_files, output_file, threads):
    merge_command = f"sambamba merge -t {threads} {output_file} {' '.join(bam_files)}"
    subprocess.run(merge_command, shell=True, check=True)
    
    sort_command = f"sambamba sort -t {threads} -o {output_file}.sorted {output_file}"
    subprocess.run(sort_command, shell=True, check=True)
    
    os.rename(f"{output_file}.sorted", output_file)
    
    index_command = f"sambamba index -t {threads} {output_file}"
    subprocess.run(index_command, shell=True, check=True)

def validate_bam_file(file_path):
    try:
        result = subprocess.run(["sambamba", "view", "-H", file_path], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False

def generate_output_filename(input_filenames):
    basenames = [os.path.splitext(os.path.basename(filename))[0] for filename in input_filenames]
    
    common_prefix = os.path.commonprefix(basenames)
    common_suffix = os.path.commonprefix([name[::-1] for name in basenames])[::-1]
    
    different_parts = []
    for basename in basenames:
        different_part = basename[len(common_prefix):len(basename)-len(common_suffix) if common_suffix else None]
        different_part = re.sub(r'[^\w\-]', '_', different_part)
        if different_part:
            different_parts.append(different_part)
    
    if different_parts:
        middle_part = '_'.join(different_parts)
        result = f"{common_prefix}{middle_part}{common_suffix}_merged.bam"
    else:
        result = f"{common_prefix}{common_suffix}_merged.bam"
    
    result = re.sub(r'_{2,}', '_', result)
    
    return result

def create_zip_file(bam_file, bai_file, zip_file):
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        zipf.write(bam_file, os.path.basename(bam_file))
        zipf.write(bai_file, os.path.basename(bai_file))

def main():
    st.title("BAM File Merger and Sorter")

    st.write("Select BAM files from the directory tree:")
    directory_structure = get_cached_directory_structure(ROOT_DIR)
    render_directory_tree(directory_structure)

    selected_files = get_selected_files(directory_structure)
    
    if selected_files:
        st.write(f"Selected {len(selected_files)} files")

        max_threads = multiprocessing.cpu_count()
        threads = st.slider("Select number of threads", min_value=1, max_value=max_threads, value=max_threads//2)

        if st.button("Merge and Sort BAM files"):
            progress_bar = st.progress(0)
            status_text = st.empty()

            with tempfile.TemporaryDirectory() as temp_dir:
                temp_bam_files = []
                original_filenames = []

                for i, file in enumerate(selected_files):
                    if validate_bam_file(file):
                        temp_bam_files.append(file)
                        original_filenames.append(os.path.basename(file))
                    else:
                        st.error(f"Invalid BAM file: {os.path.basename(file)}")
                        return

                    progress_bar.progress((i + 1) / len(selected_files) * 0.5)
                    status_text.text(f"Validating files: {i+1}/{len(selected_files)}")

                if not temp_bam_files:
                    st.error("No valid BAM files were selected.")
                    return

                output_filename = generate_output_filename(original_filenames)
                output_file = os.path.join(os.path.dirname(temp_bam_files[0]), output_filename)

                status_text.text("Merging and sorting BAM files...")
                merge_and_sort_bams(temp_bam_files, output_file, threads)
                progress_bar.progress(0.75)

                status_text.text("Creating ZIP file for download...")
                zip_filename = output_filename.replace('.bam', '.zip')
                zip_path = os.path.join(temp_dir, zip_filename)
                create_zip_file(output_file, output_file + ".bai", zip_path)

                progress_bar.progress(1.0)
                status_text.text("Processing complete. Preparing download...")

                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download merged BAM and index (ZIP)",
                        data=f,
                        file_name=zip_filename,
                        mime="application/zip"
                    )

                st.success(f"BAM files have been merged and sorted successfully. The output file is saved as {output_file}")

if __name__ == "__main__":
    main()