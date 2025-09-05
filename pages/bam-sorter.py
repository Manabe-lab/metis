import streamlit as st
import os
import tempfile
import subprocess
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

def sort_bam(input_file, output_file, threads):
    sort_command = f"sambamba sort -t {threads} -o {output_file} {input_file}"
    subprocess.run(sort_command, shell=True, check=True)

def validate_bam_file(file_path):
    try:
        result = subprocess.run(["sambamba", "view", "-H", file_path], check=True, capture_output=True)
        return True
    except subprocess.CalledProcessError:
        return False

def create_zip_file(files, zip_file):
    with zipfile.ZipFile(zip_file, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for file in files:
            zipf.write(file, os.path.basename(file))

def main():
    st.title("BAM File Sorter for METIS Data")

    st.write("Select BAM files from the directory tree:")
    directory_structure = get_cached_directory_structure(ROOT_DIR)
    render_directory_tree(directory_structure)

    selected_files = get_selected_files(directory_structure)
    
    if selected_files:
        st.write(f"Selected {len(selected_files)} files")

        max_threads = multiprocessing.cpu_count()
        threads = st.slider("Select number of threads", min_value=1, max_value=max_threads, value=max_threads//2)

        if st.button("Sort Selected BAM files"):
            progress_bar = st.progress(0)
            status_text = st.empty()

            with tempfile.TemporaryDirectory() as temp_dir:
                sorted_files = []

                for i, file in enumerate(selected_files):
                    if validate_bam_file(file):
                        output_bam = os.path.join(os.path.dirname(file), os.path.basename(file).replace('.bam', '.sorted.bam'))
                        status_text.text(f"Sorting file: {os.path.basename(file)}")
                        sort_bam(file, output_bam, threads)
                        sorted_files.append(output_bam)
                    else:
                        st.error(f"Invalid BAM file: {os.path.basename(file)}")

                    progress_bar.progress((i + 1) / len(selected_files))

                if not sorted_files:
                    st.error("No valid BAM files were processed.")
                    return

                status_text.text("Creating ZIP file for download...")
                zip_filename = "sorted_bam_files.zip"
                zip_path = os.path.join(temp_dir, zip_filename)
                create_zip_file(sorted_files, zip_path)

                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download sorted BAM files (ZIP)",
                        data=f,
                        file_name=zip_filename,
                        mime="application/zip"
                    )

                status_text.text("Processing complete.")
                progress_bar.progress(1.0)
                st.success(f"Sorted {len(sorted_files)} BAM files successfully.")

if __name__ == "__main__":
    main()