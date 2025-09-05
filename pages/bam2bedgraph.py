import streamlit as st
import subprocess
import os
import tempfile
import zipfile
import shutil
from helper_func import mk_temp_dir

# temp内に保存する
# --- Initialising SessionState ---
if "temp_dir" not in st.session_state:
    st.session_state.temp_dir = True
    temp_dir, res_dir = mk_temp_dir("bedgraph")
    st.session_state.temp_dir = temp_dir
else:
    temp_dir = st.session_state.temp_dir
    temp_dir, res_dir = mk_temp_dir("bedgraph", temp_dir)

GENOME_FILES = {
    "mm10": "/home/cellxgene/streamlit/db/genome_size/mm10.chrom.sizes",
    "mm39": "./db/genome_size/mm39.chrom.sizes",
    "hg38": "./db/genome_size/hg38.chrom.sizes"
}

def process_bam_file(bam_file, sample_name, genome_file):
    # Copy the BAM file to the temp directory
    temp_bam = os.path.join(temp_dir, f"{sample_name}.bam")
    shutil.copy(bam_file, temp_bam)

    # Process commands
    commands = [
        f"sambamba sort -n -t 16 -o {sample_name}_sorted_by_name.bam {sample_name}.bam",
        f"bedtools bamtobed -bedpe -i {sample_name}_sorted_by_name.bam > {sample_name}.bed",
        f"awk '$1==$4 && $6-$2 < 1000 {{print $0}}' {sample_name}.bed > {sample_name}.clean.bed",
        f"cut -f 1,2,6 {sample_name}.clean.bed | sort -k1,1 -k2,2n -k3,3n > {sample_name}.fragments.bed",
        f"bedtools genomecov -bg -i {sample_name}.fragments.bed -g {genome_file} > {sample_name}.PE.bedgraph",
        f"gzip -f {sample_name}.PE.bedgraph"
    ]

    for cmd in commands:
        subprocess.run(cmd, shell=True, check=True, cwd=temp_dir)

    # Return the path of the final compressed output file
    return os.path.join(temp_dir, f"{sample_name}.PE.bedgraph.gz")

def main():
    st.title("BAM File Converter")

    uploaded_files = st.file_uploader("Upload BAM files", type="bam", accept_multiple_files=True)
    genome_option = st.radio("Select genome", ("mm10", "mm39", "hg38"), index=0)

    if uploaded_files:
        if st.button("Process Files"):
            genome_file = GENOME_FILES[genome_option]
            if not os.path.exists(genome_file):
                st.error(f"Genome file not found: {genome_file}")
                return

            with tempfile.TemporaryDirectory() as temp_dir:
                output_files = []
                total_files = len(uploaded_files)

                # プログレスバーの初期化
                progress_bar = st.progress(0)
                status_text = st.empty()

                for i, uploaded_file in enumerate(uploaded_files):
                    sample_name = os.path.splitext(uploaded_file.name)[0]
                    
                    # ステータス更新
                    status_text.text(f"Processing {sample_name}... ({i+1}/{total_files})")
                    
                    bam_path = os.path.join(temp_dir, uploaded_file.name)
                    with open(bam_path, "wb") as f:
                        f.write(uploaded_file.getvalue())

                    output_file = process_bam_file(bam_path, sample_name, genome_file)
                    output_files.append(output_file)

                    # プログレスバー更新
                    progress_bar.progress((i + 1) / total_files)

                # ZIPファイルの作成
                status_text.text("Creating ZIP file...")
                zip_path = os.path.join(temp_dir, "results.zip")
                with zipfile.ZipFile(zip_path, "w") as zip_file:
                    for file in output_files:
                        zip_file.write(file, os.path.basename(file))

                # 完了メッセージ
                status_text.text("Processing complete!")

                # ダウンロードボタンの表示
                with open(zip_path, "rb") as f:
                    st.download_button(
                        label="Download Results",
                        data=f.read(),
                        file_name="results.zip",
                        mime="application/zip"
                    )


if __name__ == "__main__":
    main()