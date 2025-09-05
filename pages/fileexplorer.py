import streamlit as st
import os
import shutil
from datetime import datetime
import asyncio
import aiofiles
from concurrent.futures import ThreadPoolExecutor
import grp
import pwd
import traceback

executor = ThreadPoolExecutor(max_workers=4)

def get_parent_path(path):
    return os.path.dirname(path)

async def async_listdir(path):
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(executor, os.listdir, path)

async def async_is_dir(path):
    loop = asyncio.get_event_loop()
    return await loop.run_in_executor(executor, os.path.isdir, path)

async def get_contents(path, sort_by='name', page=1, items_per_page=50):
    try:
        all_items = await async_listdir(path)

        dirs = []
        files = []

        for item in all_items:
            item_path = os.path.join(path, item)
            if await async_is_dir(item_path):
                dirs.append(item)
            else:
                files.append(item)

        if sort_by == 'name':
            dirs.sort()
            files.sort()

        total_items = len(dirs) + len(files)
        start_index = (page - 1) * items_per_page
        end_index = start_index + items_per_page

        paginated_dirs = dirs[max(0, start_index):end_index]
        if len(paginated_dirs) < items_per_page:
            remaining = items_per_page - len(paginated_dirs)
            paginated_files = files[:remaining]
        else:
            paginated_files = []

    except PermissionError:
        st.error(f"Access denied: {path}")
        return [], [], 0

    return paginated_dirs, paginated_files, total_items

def is_subpath(path, parent):
    return os.path.abspath(path).startswith(os.path.abspath(parent))

async def get_file_content(file_path):
    try:
        async with aiofiles.open(file_path, "rb") as file:
            return await file.read()
    except Exception as e:
        st.error(f"Error reading file: {e}")
        return None

async def modify_item(item_path, item_type):
    action = st.radio("Select an action", ["Rename", "Delete"])

    if action == "Rename":
        new_name = st.text_input("Enter new name")
        if st.button("Execute Rename"):
            try:
                new_path = os.path.join(os.path.dirname(item_path), new_name)
                await asyncio.get_event_loop().run_in_executor(executor, os.rename, item_path, new_path)
                st.success("Rename completed.")
                st.session_state.selected_file = new_path
                st.rerun()
            except PermissionError:
                st.error("You do not have permission to rename.")

    elif action == "Delete":
        st.warning("Warning: This action cannot be undone.")
        if st.button("Execute Delete", type="primary"):
            try:
                if item_type == "file":
                    await asyncio.get_event_loop().run_in_executor(executor, os.remove, item_path)
                else:
                    await asyncio.get_event_loop().run_in_executor(executor, shutil.rmtree, item_path)
                st.success("Delete completed.")
                st.session_state.pop('selected_file', None)
                st.rerun()
            except PermissionError:
                st.error("You do not have permission to delete.")


async def create_directory_with_permissions(path):
    loop = asyncio.get_event_loop()
    try:
        # ディレクトリを作成（lambda関数を使用）
        await loop.run_in_executor(executor, lambda: os.makedirs(path, exist_ok=True))

        # ユーザー情報を取得
        uid = os.getuid()
        user = pwd.getpwuid(uid).pw_name

        # グループ情報を取得（例: 'sftp'グループ）
        try:
            gid = grp.getgrnam('sftp').gr_gid
        except KeyError:
            # 'sftp'グループが存在しない場合、ユーザーのプライマリグループを使用
            gid = pwd.getpwuid(uid).pw_gid

        # 所有者とグループを設定
        await loop.run_in_executor(executor, os.chown, path, uid, gid)

        # パーミッションを設定 (770 = rwxrwx---)
        await loop.run_in_executor(executor, os.chmod, path, 0o770)

        return f"Directory created successfully. Owner: {user}, Group: {grp.getgrgid(gid).gr_name}"
    except PermissionError as e:
        return f"Permission denied: {e}\n{traceback.format_exc()}"
    except Exception as e:
        return f"Error creating directory: {e}\n{traceback.format_exc()}"

async def main():
    st.title("File Explorer")

    home_dirs = {
        "SCALA": "/home/lab/sftp-data/SCALA-data/SCALA-download/Personal_folders",
        "cellxgene": "/home/lab/sftp-data/cellxgene-upload",
        "bulk": "/home/lab/sftp-data/METIS_data"
    }

    # Initialize selected_home if not present
    if 'selected_home' not in st.session_state:
        st.session_state.selected_home = list(home_dirs.values())[0]
        st.session_state.current_path = st.session_state.selected_home

    # Get the current selection
    selected_home_key = st.selectbox(
        "Select home directory", 
        list(home_dirs.keys()), 
        index=list(home_dirs.keys()).index(
            [k for k, v in home_dirs.items() if v == st.session_state.selected_home][0]
        )
    )


    # Update paths when home directory changes
    new_home = home_dirs[selected_home_key]
    if new_home != st.session_state.selected_home:
        st.session_state.selected_home = new_home
        st.session_state.current_path = new_home
        st.session_state.page = 1
        st.session_state.pop('selected_file', None)
        st.rerun()

    st.write(f"Current path: {st.session_state.current_path}")

    parent_path = get_parent_path(st.session_state.current_path)
    if parent_path != st.session_state.current_path and is_subpath(parent_path, st.session_state.selected_home):
        if st.button("Go to parent directory"):
            st.session_state.current_path = parent_path
            st.session_state.page = 1
            st.session_state.pop('selected_file', None)
            st.rerun()

    sort_method = st.selectbox("Sort by", ["Name"], index=0)
    sort_by = 'name'

    if 'page' not in st.session_state:
        st.session_state.page = 1
    items_per_page = 50

    with st.spinner("Loading directory contents..."):
        dirs, files, total_items = await get_contents(st.session_state.current_path, sort_by, st.session_state.page, items_per_page)

    st.write("Directories:")
    for dir in dirs:
        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown(f"📁 {dir}")
        with col2:
            if st.button("Open", key=f"open_{dir}"):
                st.session_state.current_path = os.path.join(st.session_state.current_path, dir)
                st.session_state.page = 1
                st.session_state.pop('selected_file', None)
                st.rerun()

    st.write("Files:")
    for file in files:
        file_path = os.path.join(st.session_state.current_path, file)
        is_selected = 'selected_file' in st.session_state and st.session_state.selected_file == file_path

        col1, col2 = st.columns([3, 1])
        with col1:
            if is_selected:
                st.markdown(f"**:blue[📄 {file}]**")
            else:
                st.markdown(f"📄 {file}")
        with col2:
            if st.button("Modify/Download", key=f"modify_download_{file}"):
                st.session_state.selected_file = file_path
                st.rerun()

    if 'selected_file' in st.session_state:
        st.write(f"Selected file: **:blue[{os.path.basename(st.session_state.selected_file)}]**")
        action = st.radio("Choose action", ["Modify", "Download"])

        if action == "Modify":
            await modify_item(st.session_state.selected_file, "file")
        elif action == "Download":
            file_content = await get_file_content(st.session_state.selected_file)
            if file_content is not None:
                st.download_button(
                    label="Download File",
                    data=file_content,
                    file_name=os.path.basename(st.session_state.selected_file),
                    mime="application/octet-stream"
                )

    # Pagination controls
    total_pages = (total_items + items_per_page - 1) // items_per_page
    col1, col2, col3 = st.columns([1, 3, 1])
    with col1:
        if st.button("Previous", disabled=(st.session_state.page == 1)):
            st.session_state.page -= 1
            st.rerun()
    with col2:
        st.write(f"Page {st.session_state.page} of {total_pages}")
    with col3:
        if st.button("Next", disabled=(st.session_state.page == total_pages)):
            st.session_state.page += 1
            st.rerun()

    st.write("Upload a file:")
    uploaded_file = st.file_uploader("Choose a file to upload", type="all")
    if uploaded_file is not None:
        file_path = os.path.join(st.session_state.current_path, uploaded_file.name)
        async with aiofiles.open(file_path, "wb") as f:
            await f.write(uploaded_file.getbuffer())
        st.success(f"File {uploaded_file.name} uploaded successfully.")
        st.rerun()

    new_dir = st.text_input("New directory name")
    if st.button("Create directory"):
        if new_dir:
            new_dir_path = os.path.join(st.session_state.current_path, new_dir)
            result = await create_directory_with_permissions(new_dir_path)
            if "successfully" in result:
                st.success(result)
            else:
                st.error(result)

            # ディレクトリ作成の詳細情報を表示
            st.write("Debug Information:")
            st.write(f"Attempted to create directory: {new_dir_path}")
            st.write(f"Current working directory: {os.getcwd()}")
            st.write(f"Current user: {pwd.getpwuid(os.getuid()).pw_name}")
            st.write(f"Current user's groups: {[g.gr_name for g in grp.getgrall() if pwd.getpwuid(os.getuid()).pw_name in g.gr_mem]}")

            try:
                parent_dir = os.path.dirname(new_dir_path)
                parent_stat = os.stat(parent_dir)
                st.write(f"Parent directory permissions: {oct(parent_stat.st_mode)[-3:]}")
                st.write(f"Parent directory owner: {pwd.getpwuid(parent_stat.st_uid).pw_name}")
                st.write(f"Parent directory group: {grp.getgrgid(parent_stat.st_gid).gr_name}")
            except Exception as e:
                st.write(f"Error getting parent directory info: {e}")

    # 現在のディレクトリの権限情報を表示
    current_path = st.session_state.current_path
    try:
        stat_info = os.stat(current_path)
        owner = pwd.getpwuid(stat_info.st_uid).pw_name
        try:
            group = grp.getgrgid(stat_info.st_gid).gr_name
        except KeyError:
            group = f"Unknown ({stat_info.st_gid})"
        permissions = oct(stat_info.st_mode)[-3:]  # 最後の3桁（8進数）を取得

        st.write(f"Current directory: {current_path}")
        st.write(f"Owner: {owner}")
        st.write(f"Group: {group}")
        st.write(f"Permissions: {permissions}")
    except Exception as e:
        st.error(f"Error getting current directory info: {e}")

if __name__ == "__main__":
    asyncio.run(main())