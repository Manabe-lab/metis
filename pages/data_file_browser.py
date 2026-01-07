"""
Data File Browser
Multiple directory file browser with different access permissions
"""

import streamlit as st
from streamlit_file_browser import st_file_browser
import os
import shutil

st.set_page_config(page_title="Data File Browser", page_icon="📁", layout="wide")

st.title("📁 Data File Browser")

st.markdown("""
sa-ba-upofDataFiletheviewbrowse,管procwithkimasu。deirekutoritoyoteakusesuauthlimitisdiffnarimasu。
""")

# Define directories with their permissions
DIRECTORIES = {
    "cellxgene-data (Read Only)": {
        "path": "/home/lab/cellxgene-data",
        "description": "cellxgenewith公opensareteexistData。viewbrowseandDownloadofmipossible。",
        "can_download": True,
        "can_upload": False,
        "can_new_folder": False,
        "can_delete": False,
    },
    "cellxgene Upload": {
        "path": "/home/lab/sftp-data/cellxgene_upload",
        "description": "cellxgeneuseofDataUploadfuoruda。FileofUpload、fuorudamakebecome、deleteremoveispossible。",
        "can_download": True,
        "can_upload": True,
        "can_new_folder": True,
        "can_delete": True,
    },
    "SCALA Personal Folders": {
        "path": "/home/lab/sftp-data/SCALA-data/SCALA-upload/Personal_folders",
        "description": "SCALAuseofpiece人fuoruda。FileofUpload、fuorudamakebecome、deleteremoveispossible。",
        "can_download": True,
        "can_upload": True,
        "can_new_folder": True,
        "can_delete": True,
    },
}

# Directory selection
selected_dir = st.selectbox(
    "deirekutoritheSelect:",
    options=list(DIRECTORIES.keys()),
    key="selected_directory"
)

dir_config = DIRECTORIES[selected_dir]
base_path = dir_config["path"]

# Show directory info
st.info(f"📂 **{selected_dir}**\n\n{dir_config['description']}")

# Show permissions
col1, col2, col3, col4 = st.columns(4)
with col1:
    if dir_config["can_download"]:
        st.success("✓ Download")
    else:
        st.error("✗ Download")
with col2:
    if dir_config["can_upload"]:
        st.success("✓ Upload")
    else:
        st.error("✗ Upload")
with col3:
    if dir_config["can_new_folder"]:
        st.success("✓ New Folder")
    else:
        st.error("✗ New Folder")
with col4:
    if dir_config["can_delete"]:
        st.success("✓ Delete")
    else:
        st.error("✗ Delete")

st.divider()

# Check if directory exists
if not os.path.exists(base_path):
    st.error(f"❌ deirekutoriisexistatshimasen: {base_path}")
    st.stop()

# ========================================
# File Browser (read-only mode)
# ========================================
st.subheader("📂 File Browser")

try:
    event = st_file_browser(
        path=base_path,
        key=f"browser_{selected_dir.replace(' ', '_').replace('(', '').replace(')', '')}",
        show_download_file=False,  # We'll handle this ourselves
        show_upload_file=False,    # We'll handle this ourselves
        show_new_folder=False,     # We'll handle this ourselves
        show_delete_file=False,    # We'll handle this ourselves
        show_choose_file=False,
        use_static_file_server=False,
        show_preview=False,
    )
except FileNotFoundError as e:
    st.error(f"❌ FileakusesuError: 壊retashinborikurinkuisarimasu。")
    event = None
except Exception as e:
    st.error(f"❌ Error: {str(e)}")
    event = None

st.divider()

# ========================================
# File Operations
# ========================================
st.subheader("📋 File Operations")

# Get selected file info
selected_file_path = None
selected_file_name = None
selected_is_dir = False

if event is not None:
    event_type = event.get('type', '')
    # Handle both file and folder selection events
    if event_type in ('SELECT_FILE', 'SELECT_FOLDER', 'DOUBLE_CLICK'):
        target = event.get('target', {})
        if target:
            relative_path = target.get('path', '')
            selected_file_name = target.get('name', '')
            selected_file_path = os.path.join(base_path, relative_path)
            selected_is_dir = os.path.isdir(selected_file_path) if selected_file_path else False

# Show selected file
if selected_file_path and os.path.exists(selected_file_path):
    file_size = os.path.getsize(selected_file_path) if not selected_is_dir else 0
    size_str = f"{file_size / (1024*1024):.1f} MB" if file_size > 0 else "Directory"

    st.success(f"Selectmid: **{selected_file_name}** ({size_str})")

    col1, col2 = st.columns(2)

    # Download button (lazy loading to avoid slow file reads)
    with col1:
        if dir_config["can_download"] and not selected_is_dir:
            if file_size > 500 * 1024 * 1024:
                st.warning("⚠️ bigkinaFilewithsu。Downloadtotimebetweeniskakarimasu。")

            # Two-step download to avoid loading file on every selection
            if st.button("📥 Downloadlevelprep", key="prepare_download_btn"):
                st.session_state.ready_to_download = selected_file_path

            if st.session_state.get('ready_to_download') == selected_file_path:
                with open(selected_file_path, 'rb') as f:
                    st.download_button(
                        label="📥 DownloadRun",
                        data=f,
                        file_name=selected_file_name,
                        key="download_btn"
                    )

    # Delete button
    with col2:
        if dir_config["can_delete"]:
            if st.button("🗑️ deleteremove", key="delete_btn", type="secondary"):
                st.session_state.confirm_delete = selected_file_path

            if st.session_state.get('confirm_delete') == selected_file_path:
                st.warning(f"maincurrenttodeleteremoveshimasuka？ **{selected_file_name}**")
                col_yes, col_no = st.columns(2)
                with col_yes:
                    if st.button("isi、deleteremovedo", key="confirm_yes", type="primary"):
                        try:
                            if selected_is_dir:
                                shutil.rmtree(selected_file_path)
                            else:
                                os.remove(selected_file_path)
                            st.success(f"✓ deleteremoveshimashita: {selected_file_name}")
                            st.session_state.confirm_delete = None
                            st.rerun()
                        except Exception as e:
                            st.error(f"deleteremoveError: {e}")
                with col_no:
                    if st.button("kiyanseru", key="confirm_no"):
                        st.session_state.confirm_delete = None
                        st.rerun()
else:
    st.info("👆 FiletheSelectshitekudasai")

st.divider()

# ========================================
# Upload Section
# ========================================
if dir_config["can_upload"]:
    st.subheader("📤 Upload")

    # Determine upload directory
    if selected_file_path:
        if selected_is_dir:
            # Folder selected - use it as upload destination
            upload_dir = selected_file_path
        else:
            # File selected - use its parent directory
            upload_dir = os.path.dirname(selected_file_path)

        # Show relative path (or "root" if it's the base path)
        if upload_dir == base_path:
            st.info(f"Uploadfirst: **ru-todeirekutori**")
        else:
            st.info(f"Uploadfirst: **{os.path.relpath(upload_dir, base_path)}/**")
    else:
        upload_dir = base_path
        st.info(f"Uploadfirst: **ru-todeirekutori**（File/fuorudatheSelectdoandchangefurtherwithkimasu）")

    uploaded_file = st.file_uploader(
        "FiletheSelect",
        key="file_uploader",
        help="UploaddoFiletheSelectshitekudasai"
    )

    if uploaded_file is not None:
        dest_path = os.path.join(upload_dir, uploaded_file.name)

        if os.path.exists(dest_path):
            st.warning(f"⚠️ samenameFileisexistatshimasu: {uploaded_file.name}")
            overwrite = st.checkbox("upwritekido", key="overwrite_check")
        else:
            overwrite = True

        if st.button("📤 UploadRun", key="upload_btn", type="primary"):
            if overwrite:
                try:
                    with open(dest_path, 'wb') as f:
                        f.write(uploaded_file.getbuffer())
                    st.success(f"✓ UploadComplete: {uploaded_file.name}")
                    st.rerun()
                except Exception as e:
                    st.error(f"UploadError: {e}")
            else:
                st.error("upwritekitheallowcanshitekudasai")

# ========================================
# New Folder Section
# ========================================
if dir_config["can_new_folder"]:
    st.subheader("📁 New Folder")

    # Determine parent directory
    if selected_file_path:
        if selected_is_dir:
            # Folder selected - create inside it
            parent_dir = selected_file_path
        else:
            # File selected - create in its parent directory
            parent_dir = os.path.dirname(selected_file_path)

        # Show relative path (or "root" if it's the base path)
        if parent_dir == base_path:
            st.info(f"makebecomefirst: **ru-todeirekutori**")
        else:
            st.info(f"makebecomefirst: **{os.path.relpath(parent_dir, base_path)}/** in")
    else:
        parent_dir = base_path
        st.info(f"makebecomefirst: **ru-todeirekutori**（File/fuorudatheSelectdoandchangefurtherwithkimasu）")

    new_folder_name = st.text_input("fuorudaname", key="new_folder_name")

    if st.button("📁 fuorudamakebecome", key="create_folder_btn"):
        if new_folder_name:
            # Sanitize folder name
            safe_name = "".join(c for c in new_folder_name if c.isalnum() or c in ('_', '-', '.', ' '))
            new_folder_path = os.path.join(parent_dir, safe_name)

            if os.path.exists(new_folder_path):
                st.error(f"samenameoffuorudaisalreadytoexistatshimasu: {safe_name}")
            else:
                try:
                    os.makedirs(new_folder_path)
                    st.success(f"✓ fuorudathemakebecomeshimashita: {safe_name}")
                    st.rerun()
                except Exception as e:
                    st.error(f"fuorudamakebecomeError: {e}")
        else:
            st.error("fuorudanametheInputshitekudasai")

# Show event details (for debugging)
if event is not None:
    with st.expander("Event details", expanded=False):
        st.json(event)
