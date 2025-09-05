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
    BEDファイルのフォーマットを検出する関数
    """
    with open(file_path, 'r') as f:
        # ヘッダーをスキップ
        line = f.readline()
        while line.startswith('#') or is_header_line(line):
            line = f.readline()
        
        fields = line.strip().split('\t')
        
        # MACSフォーマットの特徴を確認
        if len(fields) >= 4:
            try:
                # 通常のBEDフォーマット（chr, start, end）の場合
                if fields[0].startswith('chr') and fields[1].isdigit() and fields[2].isdigit():
                    return 'standard'
                # MACSフォーマットの場合（name, chr, start, end）
                elif fields[1].startswith('chr') and fields[2].isdigit() and fields[3].isdigit():
                    return 'macs'
            except:
                pass
    return 'unknown'

@st.cache_data
def convert_to_standard_bed(input_path, output_path, format_type):
    """
    verschiedene BEDフォーマットを標準的なBEDフォーマットに変換する関数
    """
    if format_type == 'macs':
        # MACSフォーマットを標準BEDフォーマットに変換
        with open(input_path, 'r') as fin, open(output_path, 'w') as fout:
            for line in fin:
                if line.startswith('#') or is_header_line(line):
                    continue
                fields = line.strip().split('\t')
                if len(fields) >= 4:
                    # chr, start, endのみを抽出して書き出し
                    bed_line = f"{fields[1]}\t{fields[2]}\t{fields[3]}\n"
                    fout.write(bed_line)
        return True
    elif format_type == 'standard':
        # 既に標準フォーマットの場合はそのままコピー
        shutil.copy(input_path, output_path)
        return True
    return False

@st.cache_data
def is_header_line(line):
    """
    行がヘッダー（カラム名）かどうかを判定する関数
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
    BEDファイルをFastaファイルに変換する関数
    """
    bed = BedTool(bed_file)
    bed.sequence(fi=f'db/genome/{genome}.fa', fo=output_file)

def show_file_preview(file_path, format_type=None):
    """
    ファイルのプレビューを表示する関数
    """
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()[:10]  # 最初の10行のみ表示
            
        st.write("File preview:")
        for i, line in enumerate(lines, 1):
            if i == 1 and is_header_line(line):
                st.markdown(f"**📋 ヘッダー行:** `{line.strip()}`")
            else:
                st.text(f"行 {i}: {line.strip()}")
                
        if format_type:
            st.info(f"Detected format: {format_type}")
            
        if len(lines) == 10:
            st.text("...")
    except Exception as e:
        st.error(f"プレビューの表示中にエラーが発生しました: {str(e)}")

def main():
    st.title('BED to FASTA Converter')

    # genomeの選択
    genome = st.selectbox(
        'Genome:',
        ('mm10', 'mm39', 'hg38')
    )

    # 複数ファイルアップロード
    uploaded_files = st.file_uploader(
        "Upload bed files", 
        type=['bed', 'broadPeak', 'narrowPeak'],
        accept_multiple_files=True
    )
    st.write("can handle MACS outputs")

    if uploaded_files:
        # 進捗状況を表示するプログレスバー
        progress_bar = st.progress(0)
        status_text = st.empty()
        
        # 変換済みファイルを格納するメモリ上のバッファ
        zip_buffer = io.BytesIO()
        
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                # ZIPファイルの作成
                with ZipFile(zip_buffer, 'w') as zip_file:
                    # 各ファイルの処理
                    for i, uploaded_file in enumerate(uploaded_files):
                        status_text.text(f'Processing: {uploaded_file.name}')
                        
                        # 進捗状況の更新
                        progress = (i + 1) / len(uploaded_files)
                        progress_bar.progress(progress)

                        # 入力ファイルの保存
                        input_path = os.path.join(temp_dir, uploaded_file.name)
                        with open(input_path, 'wb') as f:
                            f.write(uploaded_file.getbuffer())

                        # フォーマットの検出
                        format_type = detect_bed_format(input_path)
                        st.info(f"Format type of {uploaded_file.name}: {format_type}")

                        # 標準BEDフォーマットへの変換
                        standard_bed_path = os.path.join(temp_dir, f"standard_{uploaded_file.name}")
                        if not convert_to_standard_bed(input_path, standard_bed_path, format_type):
                            st.error(f"{uploaded_file.name}: サポートされていないフォーマットです。")
                            show_file_preview(input_path, format_type)
                            continue

                        # 出力ファイル名の設定（.bedを.faに変更）
                        output_filename = os.path.splitext(uploaded_file.name)[0] + '.fa'
                        output_path = os.path.join(temp_dir, output_filename)

                        try:
                            # 標準フォーマットのBEDファイルをFastaに変換
                            convert_bed_to_fasta(standard_bed_path, genome, output_path)
                            
                            # 変換したファイルをZIPに追加
                            zip_file.write(output_path, output_filename)
                            
                        except Exception as e:
                            st.error(f"""
                            エラーが発生しました:
                            
                            ファイル: {uploaded_file.name}
                            {str(e)}
                            """)
                            
                            # ファイルプレビューを表示
                            show_file_preview(standard_bed_path, format_type)
                            continue
                
                # 全ての処理が完了したら表示
                status_text.text('Done！')
                
                # ZIPファイルのダウンロードボタン
                zip_buffer.seek(0)
                st.download_button(
                    label="Download fasta files",
                    data=zip_buffer,
                    file_name="converted_fasta_files.zip",
                    mime="application/zip"
                )

                # 処理結果のサマリー表示
                st.write("### Results")
                for file in uploaded_files:
                    st.write(f"- {file.name} → {os.path.splitext(file.name)[0]}.fa")

        except Exception as e:
            st.error(f"予期せぬエラーが発生しました: {str(e)}")

if __name__ == '__main__':
    main()