import streamlit as st
import os
import sys
import time
from datetime import datetime
import asyncio
import aiofiles
import shutil
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import threading
from typing import Optional
import tempfile

# Streamlitの設定
st.set_page_config(
    page_title="File Uploader",
    page_icon="📁",
    layout="wide",
    initial_sidebar_state="expanded",
)

# 定数設定
CHUNK_SIZE = 1024 * 1024  # 1MB
TIMEOUT = 300  # 5分
MAX_RETRIES = 3
RETRY_DELAY = 5  # 5秒

def list_dirs(dir_path):
    try:
        return [name for name in os.listdir(dir_path) if os.path.isdir(os.path.join(dir_path, name))]
    except Exception as e:
        st.error(f"Error listing directories: {str(e)}")
        return []

async def save_uploaded_file_with_retry(uploaded_file, destination_path, processed_size, total_size, total_progress):
    for attempt in range(MAX_RETRIES):
        try:
            # 一時ファイルを使用して保存
            with tempfile.NamedTemporaryFile(delete=False) as temp_file:
                temp_path = temp_file.name
                
            async with aiofiles.open(temp_path, "wb") as temp_out:
                buffer = uploaded_file.getbuffer()
                position = 0
                
                while position < len(buffer):
                    chunk = buffer[position:position + CHUNK_SIZE]
                    await temp_out.write(chunk)
                    position += len(chunk)
                    
                    # 全体の進捗を更新
                    total_processed = processed_size + position
                    total_progress.progress(min(1.0, total_processed / total_size))

            # 保存先ディレクトリの作成とファイルの移動
            os.makedirs(os.path.dirname(destination_path), exist_ok=True)
            shutil.move(temp_path, destination_path)
            return True

        except Exception as e:
            if attempt < MAX_RETRIES - 1:
                await asyncio.sleep(RETRY_DELAY)
                continue
            else:
                st.error(f"Failed to save file after {MAX_RETRIES} attempts: {str(e)}")
                return False
        
        finally:
            if 'temp_path' in locals() and os.path.exists(temp_path):
                try:
                    os.unlink(temp_path)
                except:
                    pass

async def process_upload_queue(upload_queue, base_path, total_size, total_progress):
    processed_size = 0
    while True:
        try:
            upload_info = await upload_queue.get()
            if upload_info is None:  # 終了シグナル
                break
                
            uploaded_file = upload_info
            relative_path = uploaded_file.name.replace('\\', '/')
            if '..' in relative_path:
                st.error(f"Invalid path detected in {relative_path}")
                continue
                
            destination_path = os.path.join(base_path, relative_path)
            success = await save_uploaded_file_with_retry(
                uploaded_file, 
                destination_path,
                processed_size,
                total_size,
                total_progress
            )
            
            if success:
                processed_size += uploaded_file.size
            
        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
        finally:
            upload_queue.task_done()

def process_uploaded_files(uploaded_files, base_path):
    async def main():
        upload_queue = asyncio.Queue()
        
        # 全体のサイズを計算
        total_size = sum(file.size for file in uploaded_files)
        st.write(f"Total size to upload: {total_size / (1024 * 1024):.2f} MB")
        
        # 全体の進捗バー
        progress_container = st.empty()
        total_progress = progress_container.progress(0)
        start_time = datetime.now()
        
        # ワーカータスクの作成
        workers = [asyncio.create_task(process_upload_queue(upload_queue, base_path, total_size, total_progress))
                  for _ in range(3)]  # 同時処理数を3に制限
        
        # ファイルをキューに追加
        for file in uploaded_files:
            await upload_queue.put(file)
            
        # 終了シグナルの送信
        for _ in workers:
            await upload_queue.put(None)
            
        # すべてのワーカーの完了を待つ
        await asyncio.gather(*workers)
        
        # 完了時の表示
        end_time = datetime.now()
        duration = (end_time - start_time).total_seconds()
        progress_container.empty()
        st.success(f"🎉 Upload completed! ({len(uploaded_files)} files, {total_size / (1024 * 1024):.2f} MB in {duration:.1f} seconds)")
    
    asyncio.run(main())

# セッション状態の初期化
if 'key_count' not in st.session_state:
   st.session_state['key_count'] = 1

if 'SCALAorcellxgene' not in st.session_state:
    st.session_state['SCALAorcellxgene'] = "SCALA"

st.markdown("### Upload files for SCALA or cellxgene\n")
SCALAorcellxgene = st.radio(
    "Upload a file for SCALA or cellxgene or CHIP-seq",
    ["SCALA", "cellxgene", "Bulk"], 
    label_visibility='collapsed',
    index=0
)

# ベースディレクトリの設定
if 'dir_path' not in st.session_state:
    st.session_state['dir_path'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders'
    st.session_state['home_dir'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders'

# アップロード先の切り替え
if SCALAorcellxgene != st.session_state['SCALAorcellxgene']:
    if SCALAorcellxgene == "SCALA":
        st.session_state['home_dir'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders'
    elif SCALAorcellxgene == "cellxgene":
        st.session_state['home_dir'] = '/home/lab/sftp-data/SCALA-data/SCALA-download/cellxgene-data/Uploaded'
    elif SCALAorcellxgene == "Bulk":
        st.session_state['home_dir'] = '/home/lab/sftp-data/METIS_data'

    st.session_state['dir_path'] = st.session_state['home_dir']
    st.session_state['SCALAorcellxgene'] = SCALAorcellxgene
    st.rerun()

# ディレクトリ選択UI
st.markdown("### Select the folder to store the file\n")
st.markdown("##### Now selected directory is " + st.session_state.dir_path.replace(st.session_state['home_dir'], ''))
dirs = list_dirs(st.session_state['dir_path'])
dirs.sort()
dirs.insert(0, '')
dirs.insert(1, '..')
st.write("##### Select a directory:")
choice = st.selectbox('Select a directory', dirs, key=st.session_state['key_count'], label_visibility='collapsed')

if (choice != '') and (choice != '..'):
    st.session_state['dir_path'] = os.path.join(st.session_state['dir_path'], choice)
    st.session_state['key_count'] += 2
elif (choice == '..') and (len(st.session_state['dir_path']) > len(st.session_state['home_dir'])):
    st.session_state['dir_path'] = os.path.dirname(st.session_state['dir_path'])
    st.session_state['key_count'] += 2

# 新規ディレクトリ作成UI
st.markdown("###### To make a new directory in " + st.session_state.dir_path.replace(st.session_state['home_dir'], '') + ':')
new_dir_name = st.text_input('Type a name for the new directory', value="", key=st.session_state['key_count'] + 1)

if st.button("Make a new directory"):
    if new_dir_name:
        st.session_state['key_count'] += 2
        new_dir_path = os.path.join(st.session_state['dir_path'], new_dir_name)
        if not os.path.exists(new_dir_path):
            try:
                os.makedirs(new_dir_path)
                st.success(f'New directory `{new_dir_name}` has been created.')
                st.session_state["dir_path"] = new_dir_path
            except Exception as e:
                st.error(f"Error creating directory: {str(e)}")

# ファイルアップロードUI
st.write("")
target_dir = st.session_state["dir_path"].replace(st.session_state['home_dir'], "")
if st.button(f'Select this folder: {target_dir}', type='primary') or st.session_state["dir_path"] != st.session_state['home_dir']:
    if st.session_state["dir_path"] == st.session_state['home_dir']:
        sys.exit(1)
    st.markdown("#### Uploaded files will be in " + st.session_state['SCALAorcellxgene'] + ": " + 
                st.session_state['dir_path'].replace(st.session_state['home_dir'], ""))

    st.markdown('---\n')

    # ファイルアップローダー
    uploaded_files = st.file_uploader(
        "Choose files or folders to upload", 
        type=None, 
        accept_multiple_files=True,
        help="You can drag and drop multiple files or folders here. Folder upload is supported in modern browsers."
    )

    if uploaded_files:
        try:
            process_uploaded_files(uploaded_files, st.session_state["dir_path"])
        except Exception as e:
            st.error(f'Upload failed: {str(e)}')