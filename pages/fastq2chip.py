"""
FASTQ -> ChIP-seq / CUT&RUN / ATAC-seq BAM

Data sources:
  1. SFTP server: Navigate run directories to select FASTQ files
  2. Azenta OSS: Download from Alibaba Cloud OSS via ossutil
  3. Upload: Upload local FASTQ files

Pipeline:
  fastp -> STAR (EndToEnd) -> sort/dedup -> HOMER -> [CUT&RUN: fragment bedGraph]
"""

import streamlit as st
import paramiko
import os
import json
import time
import subprocess
import stat as stat_module
import re
from pathlib import Path
from datetime import datetime
from helper_func import clear_old_directories

# Relax Streamlit static file serving limits (for large files like BAM)
from streamlit.web.server import app_static_file_handler as _asfh
_asfh.MAX_APP_STATIC_FILE_SIZE = 10 * 1024 * 1024 * 1024  # 10 GB
_asfh.SAFE_APP_STATIC_FILE_EXTENSIONS = (
    *_asfh.SAFE_APP_STATIC_FILE_EXTENSIONS,
    ".bam", ".bai", ".bedgraph", ".bw", ".txt", ".tsv", ".csv",
)

# =============================================================================
#  Constants
# =============================================================================

_STREAMLIT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS_DIR = os.path.join(_STREAMLIT_DIR, "temp", "fastq2chip_jobs")
WORKER_SCRIPT = os.path.join(_STREAMLIT_DIR, "fastq2chip_worker.py")
PYTHON_PATH = "python3"

GENOME_OPTIONS = ["mm10", "hg38"]
ASSAY_TYPES = ["chipseq", "cutandrun", "atacseq"]
ASSAY_LABELS = {
    "chipseq": "ChIP-seq",
    "cutandrun": "CUT&RUN",
    "atacseq": "ATAC-seq",
}

STATUS_COLORS = {
    "queued": "gray",
    "downloading": "blue",
    "trimming": "orange",
    "mapping": "orange",
    "dedup": "violet",
    "homer": "violet",
    "bedgraph": "violet",
    "finalizing": "violet",
    "done": "green",
    "error": "red",
}

os.makedirs(JOBS_DIR, exist_ok=True)

# =============================================================================
#  SFTP Connection (Aging)
# =============================================================================

def get_sftp_connection():
    """Get SFTP connection (cached in session_state)"""
    if "sftp_ssh" in st.session_state and st.session_state.sftp_ssh is not None:
        try:
            transport = st.session_state.sftp_ssh.get_transport()
            if transport is not None and transport.is_active():
                return st.session_state.sftp_ssh, st.session_state.sftp_client
        except Exception:
            pass

    try:
        cred = st.secrets["aging_sftp"]
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(cred["host"], username=cred["username"], password=cred["password"])
        sftp = ssh.open_sftp()
        st.session_state.sftp_ssh = ssh
        st.session_state.sftp_client = sftp
        return ssh, sftp
    except Exception as e:
        st.error(f"SFTP connection failed: {e}")
        return None, None


def list_remote_dir(sftp, path):
    """List the contents of a remote directory"""
    try:
        items = sftp.listdir_attr(path)
        dirs = []
        files = []
        for item in items:
            if stat_module.S_ISDIR(item.st_mode):
                dirs.append(item.filename)
            else:
                files.append((item.filename, item.st_size))
        return sorted(dirs), sorted(files, key=lambda x: x[0])
    except Exception as e:
        st.error(f"Cannot list directory: {e}")
        return [], []


def scan_directory_tree(sftp, path):
    """Scan the run directory tree focusing on FASTQ files"""
    projects = []
    try:
        top_items = sftp.listdir_attr(path)
        top_dirs = {i.filename: i for i in top_items if stat_module.S_ISDIR(i.st_mode)}

        if "SampleSheet" in top_dirs:
            base = f"{path}/SampleSheet"
        else:
            base = path

        proj_items = sftp.listdir_attr(base)
        proj_dirs = sorted(
            [i for i in proj_items if stat_module.S_ISDIR(i.st_mode)],
            key=lambda x: x.filename,
        )

        for pd in proj_dirs:
            proj_name = pd.filename
            proj_path = f"{base}/{proj_name}"
            samples = []

            other_path = f"{proj_path}/Other"
            try:
                sample_items = sftp.listdir_attr(other_path)
            except Exception:
                other_path = proj_path
                try:
                    sample_items = sftp.listdir_attr(other_path)
                except Exception:
                    continue

            sample_dirs = sorted(
                [i for i in sample_items if stat_module.S_ISDIR(i.st_mode)],
                key=lambda x: x.filename,
            )

            for sd in sample_dirs:
                sample_path = f"{other_path}/{sd.filename}"
                try:
                    files = sftp.listdir_attr(sample_path)
                except Exception:
                    continue
                fastqs = [
                    (f.filename, f.st_size)
                    for f in files
                    if f.filename.endswith((".fastq.gz", ".fq.gz"))
                ]
                if fastqs:
                    samples.append({"name": sd.filename, "files": sorted(fastqs)})

            if samples:
                projects.append({"project": proj_name, "samples": samples})

    except Exception:
        pass

    return projects


def format_project_tree(projects):
    """Format the project tree for text display"""
    lines = []
    for proj in projects:
        total_size = sum(s for sample in proj["samples"] for _, s in sample["files"])
        n_samples = len(proj["samples"])
        lines.append(f"📁 {proj['project']}/  ({n_samples} samples, {format_size(total_size)})")
        for sample in proj["samples"]:
            lines.append(f"    📁 {sample['name']}/")
            for fname, sz in sample["files"]:
                lines.append(f"        🧬 {fname}  ({format_size(sz)})")
    return "\n".join(lines)


def scan_samples(sftp, project_dir):
    """Detect samples within a project directory"""
    samples = []

    # Pattern 1: Samples inside Other/ subdirectory
    other_dir = f"{project_dir}/Other"
    try:
        sample_dirs = sftp.listdir_attr(other_dir)
        for sd in sorted(sample_dirs, key=lambda x: x.filename):
            if not stat_module.S_ISDIR(sd.st_mode):
                continue
            sample_path = f"{other_dir}/{sd.filename}"
            fastq = []
            try:
                for f in sftp.listdir_attr(sample_path):
                    if f.filename.endswith((".fastq.gz", ".fq.gz")):
                        fastq.append((f.filename, f.st_size))
            except:
                continue
            if fastq:
                r1 = [f for f in fastq if "_R1" in f[0] or "_1." in f[0]]
                r2 = [f for f in fastq if "_R2" in f[0] or "_2." in f[0]]
                samples.append({
                    "name": sd.filename,
                    "remote_dir": sample_path,
                    "files": fastq,
                    "r1": r1[0] if r1 else None,
                    "r2": r2[0] if r2 else None,
                })
    except:
        pass

    # Pattern 2: fastq.gz files directly under project
    if not samples:
        try:
            for f in sftp.listdir_attr(project_dir):
                if f.filename.endswith((".fastq.gz", ".fq.gz")):
                    samples.append({
                        "name": f.filename,
                        "remote_dir": project_dir,
                        "files": [(f.filename, f.st_size)],
                        "r1": (f.filename, f.st_size) if "_R1" in f.filename or "_1." in f.filename else None,
                        "r2": (f.filename, f.st_size) if "_R2" in f.filename or "_2." in f.filename else None,
                    })
        except:
            pass

    return samples


# =============================================================================
#  Azenta OSS Utilities
# =============================================================================

def list_oss_files(oss_path, access_key_id, access_key_secret, endpoint):
    """List files in an OSS path using ossutil"""
    cmd = [
        "ossutil", "ls", oss_path,
        "-i", access_key_id,
        "-k", access_key_secret,
        "-e", endpoint,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
    if result.returncode != 0:
        return None, result.stderr

    files = []
    for line in result.stdout.strip().split("\n"):
        if not line.strip() or line.startswith("Last") or line.startswith("Object"):
            continue
        parts = line.split()
        if len(parts) >= 8 and parts[-1].startswith("oss://"):
            obj_path = parts[-1]
            fname = os.path.basename(obj_path.rstrip("/"))
            try:
                size = int(parts[4])
            except (ValueError, IndexError):
                size = 0
            if fname and size > 0:
                files.append({"name": fname, "path": obj_path, "size": size})
    return files, None


def parse_oss_samples(files):
    """Detect paired-end samples from OSS file listing"""
    fastq_files = [f for f in files if f["name"].endswith((".fastq.gz", ".fq.gz"))]

    paired = {}
    for f in fastq_files:
        name = f["name"]
        m = re.match(r'(.+)_(?:L\d+_)?1\.(f(?:ast)?q\.gz)', name)
        if m:
            key = m.group(1)
            if key not in paired:
                paired[key] = {"name": key, "r1": None, "r2": None}
            paired[key]["r1"] = f
            continue
        m = re.match(r'(.+)_(?:L\d+_)?2\.(f(?:ast)?q\.gz)', name)
        if m:
            key = m.group(1)
            if key not in paired:
                paired[key] = {"name": key, "r1": None, "r2": None}
            paired[key]["r2"] = f
            continue
        m = re.match(r'(.+)_R1(.+)', name)
        if m:
            key = m.group(1)
            if key not in paired:
                paired[key] = {"name": key, "r1": None, "r2": None}
            paired[key]["r1"] = f
            continue
        m = re.match(r'(.+)_R2(.+)', name)
        if m:
            key = m.group(1)
            if key not in paired:
                paired[key] = {"name": key, "r1": None, "r2": None}
            paired[key]["r2"] = f
            continue

    return list(paired.values())


# =============================================================================
#  Common Utilities
# =============================================================================

def format_size(size_bytes):
    """Format file size into a human-readable string"""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


LAB_NGS_BASE = os.environ.get("NGS_DATA_DIR", os.path.expanduser("~/ngs_data"))


def _build_dir_tree(path, prefix="", max_depth=3, depth=0):
    """Generate a tree string of directories only"""
    if depth >= max_depth:
        return []
    try:
        entries = sorted([
            d for d in os.listdir(path)
            if os.path.isdir(os.path.join(path, d)) and not d.startswith(".")
        ])
    except Exception:
        return []
    lines = []
    for i, name in enumerate(entries):
        is_last = (i == len(entries) - 1)
        connector = "└── " if is_last else "├── "
        lines.append(f"{prefix}{connector}{name}/")
        child_prefix = prefix + ("    " if is_last else "│   ")
        lines.extend(_build_dir_tree(
            os.path.join(path, name), child_prefix, max_depth, depth + 1,
        ))
    return lines


def pipeline_config_ui(key_prefix, default_save_to_lab=False):
    """Pipeline configuration UI"""
    col1, col2 = st.columns(2)
    with col1:
        genome = st.selectbox("Genome", GENOME_OPTIONS, index=0, key=f"{key_prefix}_genome")
    with col2:
        assay_type = st.selectbox(
            "Assay type", ASSAY_TYPES,
            format_func=lambda x: ASSAY_LABELS[x],
            index=0,
            key=f"{key_prefix}_assay_type",
        )

    if assay_type == "cutandrun":
        st.info("CUT&RUN: A fragment bedGraph for SEACR will also be generated.")
    elif assay_type == "atacseq":
        st.info("ATAC-seq: chrM reads will be removed.")

    # --- Save directory selection ---
    st.markdown("**Save directory**")
    lab_ngs_available = os.path.isdir(LAB_NGS_BASE)

    if lab_ngs_available:
        save_to_lab = st.checkbox(
            "Save to lab_NGS (shared storage)",
            value=default_save_to_lab,
            key=f"{key_prefix}_save_to_lab",
        )
    else:
        save_to_lab = False

    if save_to_lab:
        # --- Directory browser ---
        browse_key = f"{key_prefix}_lab_path"
        if browse_key not in st.session_state:
            st.session_state[browse_key] = LAB_NGS_BASE

        current_path = st.session_state[browse_key]

        # Display current path
        rel = os.path.relpath(current_path, LAB_NGS_BASE)
        if rel == ".":
            st.markdown("📁 **lab_NGS/**")
        else:
            st.markdown(f"📁 **lab_NGS/{rel}/**")

        # Folder tree display
        if current_path != LAB_NGS_BASE:
            tree_lines = _build_dir_tree(current_path)
            if tree_lines:
                rel_name = os.path.basename(current_path)
                tree_text = f"{rel_name}/\n" + "\n".join(tree_lines)
                with st.expander("📂 Folder structure", expanded=True):
                    st.code(tree_text, language="text")
            else:
                st.caption("(no subdirectories)")

        # Subdirectory listing
        try:
            subdirs = sorted([
                d for d in os.listdir(current_path)
                if os.path.isdir(os.path.join(current_path, d)) and not d.startswith(".")
            ])
        except Exception:
            subdirs = []

        col_nav1, col_nav2, col_nav3 = st.columns([4, 1, 1])
        with col_nav1:
            if subdirs:
                selected_sub = st.selectbox(
                    "Folder", subdirs,
                    key=f"{key_prefix}_subdir_sel",
                )
            else:
                selected_sub = None
                st.caption("(empty)")
        with col_nav2:
            if selected_sub and st.button("📂 Open", key=f"{key_prefix}_open_dir"):
                st.session_state[browse_key] = os.path.join(current_path, selected_sub)
                st.rerun()
        with col_nav3:
            if current_path != LAB_NGS_BASE:
                if st.button("⬆ Back", key=f"{key_prefix}_up_dir"):
                    st.session_state[browse_key] = os.path.dirname(current_path)
                    st.rerun()

        # Create new folder
        col_new1, col_new2 = st.columns([4, 1])
        with col_new1:
            new_folder_name = st.text_input(
                "New folder",
                key=f"{key_prefix}_new_folder_name",
                placeholder="e.g. 2026-02_CUT-RUN",
                label_visibility="collapsed",
            )
        with col_new2:
            if st.button("📁+ Create", key=f"{key_prefix}_create_folder"):
                if new_folder_name:
                    new_path = os.path.join(current_path, new_folder_name)
                    os.makedirs(new_path, exist_ok=True)
                    st.session_state[browse_key] = new_path
                    st.rerun()

        local_dir = current_path
        st.success(f"Save to: `{local_dir}`")
    else:
        local_dir = st.text_input(
            "Local directory",
            value=os.path.join(str(Path.home()), "NGS", "metis") + "/",
            key=f"{key_prefix}_local_dir",
        )

    include_homer_tags = st.checkbox(
        "Include Homer tag directories in results",
        value=False,
        key=f"{key_prefix}_include_homer_tags",
        help="Homer tag directories are large. "
             "Enable to include them in collected results.",
    )

    cleanup = st.checkbox(
        "Cleanup after completion (delete FASTQ and intermediate files)",
        value=True,
        key=f"{key_prefix}_cleanup",
    )

    return genome, assay_type, local_dir, save_to_lab, include_homer_tags, cleanup


# =============================================================================
#  Job Management
# =============================================================================

def load_job_status(job_dir):
    status_path = os.path.join(job_dir, "status.json")
    if os.path.exists(status_path):
        with open(status_path) as f:
            return json.load(f)
    return {"status": "unknown", "step": "", "progress_pct": 0}


def load_job_config(job_dir):
    config_path = os.path.join(job_dir, "config.json")
    if os.path.exists(config_path):
        with open(config_path) as f:
            return json.load(f)
    return {}


def get_all_jobs():
    if not os.path.isdir(JOBS_DIR):
        return []
    jobs = []
    for name in os.listdir(JOBS_DIR):
        job_dir = os.path.join(JOBS_DIR, name)
        if os.path.isdir(job_dir):
            status = load_job_status(job_dir)
            config = load_job_config(job_dir)
            jobs.append({"name": name, "dir": job_dir, "status": status, "config": config})
    jobs.sort(key=lambda x: x["name"], reverse=True)
    return jobs


def submit_job(config):
    job_id = str(int(time.time()))
    job_dir = os.path.join(JOBS_DIR, job_id)
    os.makedirs(job_dir, exist_ok=True)

    config["job_id"] = job_id
    config["created_at"] = datetime.now().isoformat(timespec='seconds')
    with open(os.path.join(job_dir, "config.json"), "w") as f:
        json.dump(config, f, indent=2, ensure_ascii=False)

    with open(os.path.join(job_dir, "status.json"), "w") as f:
        json.dump({
            "status": "queued", "step": "Waiting to start",
            "progress_pct": 0, "started_at": None, "finished_at": None, "error": None,
        }, f, indent=2)

    log_file = open(os.path.join(job_dir, "nohup.out"), "w")
    env = os.environ.copy()
    subprocess.Popen(
        [PYTHON_PATH, WORKER_SCRIPT, job_dir],
        stdout=log_file, stderr=subprocess.STDOUT, start_new_session=True,
        env=env,
    )
    return job_id


# =============================================================================
#  UI
# =============================================================================

st.title("FASTQ to ChIP/CUT&RUN/ATAC BAM")

tab_sftp, tab_oss, tab_upload, tab_monitor = st.tabs(["Aging SFTP", "Azenta OSS", "Upload", "Job Monitor"])

# =============================================================================
#  Tab 1: Aging SFTP
# =============================================================================

with tab_sftp:
    st.subheader("1. Select FASTQ files from Aging server")
    st.caption("Navigate: Nova-seq → {run} → SampleSheet → **{project}** → Scan")
    st.caption("Each project within the SampleSheet (e.g., Aoki_ChIP-seq) should be submitted as a separate job.")

    remote_dir = st.text_input(
        "Remote directory path",
        value=st.session_state.get("chip_sftp_remote_dir", "/home/Aging/Nova-seq"),
        key="chip_remote_dir_input",
    )

    if st.button("Browse", type="primary", key="chip_sftp_browse"):
        st.session_state.chip_sftp_remote_dir = remote_dir

    if "chip_sftp_remote_dir" in st.session_state:
        ssh, sftp = get_sftp_connection()
        if sftp is not None:
            current_dir = st.session_state.chip_sftp_remote_dir
            st.caption(f"Current: `{current_dir}`")
            dirs, files = list_remote_dir(sftp, current_dir)

            if dirs:
                selected_run = st.selectbox(
                    "Run directory",
                    ["-- select --"] + dirs,
                    key="chip_subdir_nav",
                )

                if selected_run != "-- select --":
                    prev = st.session_state.get("_chip_prev_subdir_nav")
                    if prev != selected_run:
                        st.session_state._chip_prev_subdir_nav = selected_run
                        target = f"{current_dir}/{selected_run}"
                        with st.spinner("Scanning directory tree..."):
                            projects = scan_directory_tree(sftp, target)
                            st.session_state.chip_sftp_dir_tree = (selected_run, target, projects)

                if "chip_sftp_dir_tree" in st.session_state and len(st.session_state.chip_sftp_dir_tree) == 3:
                    tree_name, tree_base, tree_projects = st.session_state.chip_sftp_dir_tree
                    if tree_projects:
                        with st.expander(f"**{tree_name}/ project structure**", expanded=True):
                            st.code(format_project_tree(tree_projects), language="text")

                        proj_names = [p["project"] for p in tree_projects]
                        selected_proj = st.selectbox(
                            "Project to process",
                            ["-- select --"] + proj_names,
                            key="chip_proj_nav",
                        )
                        if selected_proj != "-- select --":
                            proj_path = f"{tree_base}/SampleSheet/{selected_proj}"
                            try:
                                sftp.stat(proj_path)
                            except Exception:
                                proj_path = f"{tree_base}/{selected_proj}"

                            prev_proj = st.session_state.get("_chip_prev_proj_nav")
                            if prev_proj != selected_proj:
                                st.session_state._chip_prev_proj_nav = selected_proj
                                st.session_state.chip_sftp_remote_dir = proj_path
                                with st.spinner("Scanning samples..."):
                                    samples = scan_samples(sftp, proj_path)
                                    st.session_state.chip_sftp_samples = samples

                            if "chip_sftp_samples" in st.session_state and st.session_state.chip_sftp_samples:
                                samples = st.session_state.chip_sftp_samples
                                st.success(f"**{len(samples)} samples found**")
                                for s in samples:
                                    paired = "paired" if s['r1'] and s['r2'] else "INCOMPLETE"
                                    st.text(f"  {s['name']}: {paired}  R1={format_size(s['r1'][1]) if s['r1'] else '?'} R2={format_size(s['r2'][1]) if s['r2'] else '?'}")
                                total_size = sum(f[1] for s in samples for f in s["files"])
                                st.info(f"Total download size: **{format_size(total_size)}**")
                            elif "chip_sftp_samples" in st.session_state:
                                st.warning("No FASTQ samples found in this project.")

                    elif tree_name:
                        st.warning(f"No FASTQ samples found in {tree_name}/.")

    st.markdown("---")

    st.subheader("2. Pipeline Configuration")
    sftp_genome, sftp_assay_type, sftp_local_dir, sftp_save_to_lab, sftp_include_homer, sftp_cleanup = pipeline_config_ui("chip_sftp", default_save_to_lab=True)

    st.markdown("---")

    st.subheader("3. Submit")

    has_samples = "chip_sftp_samples" in st.session_state and len(st.session_state.get("chip_sftp_samples", [])) > 0

    if has_samples:
        samples = st.session_state.chip_sftp_samples
        paired_samples = [s for s in samples if s["r1"] and s["r2"]]
        project_name = os.path.basename(st.session_state.get("chip_sftp_remote_dir", "").rstrip("/"))
        st.markdown(f"**Project: {project_name}** — {len(paired_samples)} paired-end samples / {ASSAY_LABELS[sftp_assay_type]} / {sftp_genome}")

        if not paired_samples:
            st.error("No paired-end samples found.")
        elif st.button("Submit Job", type="primary", key="chip_sftp_submit"):
            cred = st.secrets["aging_sftp"]
            project_name = os.path.basename(st.session_state.chip_sftp_remote_dir.rstrip("/"))
            sftp_samples = []
            for s in paired_samples:
                sftp_samples.append({
                    "name": s["name"],
                    "remote_dir": s["remote_dir"],
                    "r1": s["r1"][0],
                    "r2": s["r2"][0],
                })
            config = {
                "data_source": "sftp",
                "project_name": project_name,
                "sftp_host": cred["host"],
                "sftp_username": cred["username"],
                "sftp_password": cred["password"],
                "sftp_samples": sftp_samples,
                "local_dir": sftp_local_dir,
                "save_to_lab": sftp_save_to_lab,
                "include_homer_tags": sftp_include_homer,
                "genome": sftp_genome,
                "assay_type": sftp_assay_type,
                "cleanup": sftp_cleanup,
            }
            job_id = submit_job(config)
            st.success(f"Job submitted! ID: {job_id}")
            st.info("Processing continues even if you close the browser. Check progress in the Job Monitor tab.")
    else:
        st.warning("First connect to the SFTP server and scan for FASTQ files.")


# =============================================================================
#  Tab 2: Azenta OSS
# =============================================================================

with tab_oss:
    st.subheader("1. Azenta/GENEWIZ OSS Settings")
    st.caption("Enter the OSS information from the GENEWIZ/Azenta email.")

    oss_path = st.text_input(
        "OSS Path (Add to OSS)",
        value=st.session_state.get("chip_oss_path", ""),
        key="chip_oss_path_input",
        placeholder="oss://gwjapan/2026.2/60-1286302729/N2602837_60-1286302729_20260213152423",
    )

    col_oss1, col_oss2 = st.columns(2)
    with col_oss1:
        oss_key_id = st.text_input("AccessKeyId", key="chip_oss_key_id")
    with col_oss2:
        oss_key_secret = st.text_input("AccessKeySecret", type="password", key="chip_oss_key_secret")

    oss_region = st.selectbox(
        "Region",
        ["oss-ap-northeast-1 (Tokyo)", "oss-ap-southeast-1 (Singapore)", "oss-cn-shanghai"],
        index=0,
        key="chip_oss_region",
    )
    oss_endpoint = oss_region.split(" ")[0] + ".aliyuncs.com"

    if st.button("List files", type="primary", key="chip_oss_list"):
        if not oss_key_id or not oss_key_secret:
            st.error("Please enter AccessKeyId and AccessKeySecret.")
        else:
            with st.spinner("Listing OSS files..."):
                path = oss_path.rstrip("/") + "/"
                files, err = list_oss_files(path, oss_key_id, oss_key_secret, oss_endpoint)
                if err:
                    st.error(f"OSS error: {err}")
                elif files:
                    st.session_state.chip_oss_files = files
                    st.session_state.chip_oss_path = oss_path
                else:
                    st.warning("No files found at this path.")

    if "chip_oss_files" in st.session_state and st.session_state.chip_oss_files:
        files = st.session_state.chip_oss_files
        fastq_files = [f for f in files if f["name"].endswith((".fastq.gz", ".fq.gz"))]
        other_files = [f for f in files if not f["name"].endswith((".fastq.gz", ".fq.gz"))]

        st.success(f"**{len(fastq_files)} FASTQ files** / {len(other_files)} other files")

        samples = parse_oss_samples(files)
        st.session_state.chip_oss_samples = samples

        if samples:
            for s in samples:
                paired = "paired" if s["r1"] and s["r2"] else "INCOMPLETE"
                r1_size = format_size(s["r1"]["size"]) if s["r1"] else "?"
                r2_size = format_size(s["r2"]["size"]) if s["r2"] else "?"
                st.text(f"  {s['name']}: {paired}  R1={r1_size} R2={r2_size}")

            total_size = sum(
                f["size"] for s in samples for f in [s["r1"], s["r2"]] if f
            )
            st.info(f"Total download size: **{format_size(total_size)}**")
        else:
            st.warning("Could not pair FASTQ files.")

        if other_files:
            with st.expander(f"Other files ({len(other_files)})"):
                for f in other_files:
                    st.text(f"  {f['name']} ({format_size(f['size'])})")

    st.markdown("---")

    st.subheader("2. Pipeline Configuration")
    oss_genome, oss_assay_type, oss_local_dir, oss_save_to_lab, oss_include_homer, oss_cleanup = pipeline_config_ui("chip_oss")

    st.markdown("---")

    st.subheader("3. Submit")

    has_oss_samples = "chip_oss_samples" in st.session_state and len(st.session_state.get("chip_oss_samples", [])) > 0

    if has_oss_samples:
        samples = st.session_state.chip_oss_samples
        paired_samples = [s for s in samples if s["r1"] and s["r2"]]
        oss_project = os.path.basename(st.session_state.get("chip_oss_path", "").rstrip("/"))
        st.markdown(f"**Project: {oss_project}** — {len(paired_samples)} paired-end samples / {ASSAY_LABELS[oss_assay_type]} / {oss_genome}")

        if not paired_samples:
            st.error("No paired-end samples found.")
        elif not oss_key_id or not oss_key_secret:
            st.error("Please enter AccessKeyId and AccessKeySecret.")
        elif st.button("Submit Job", type="primary", key="chip_oss_submit"):
            oss_samples_config = []
            for s in paired_samples:
                oss_samples_config.append({
                    "name": s["name"],
                    "r1_name": s["r1"]["name"],
                    "r1_path": s["r1"]["path"],
                    "r2_name": s["r2"]["name"],
                    "r2_path": s["r2"]["path"],
                })
            config = {
                "data_source": "oss",
                "project_name": oss_project,
                "oss_path": oss_path.rstrip("/") + "/",
                "oss_access_key_id": oss_key_id,
                "oss_access_key_secret": oss_key_secret,
                "oss_endpoint": oss_endpoint,
                "oss_samples": oss_samples_config,
                "local_dir": oss_local_dir,
                "save_to_lab": oss_save_to_lab,
                "include_homer_tags": oss_include_homer,
                "genome": oss_genome,
                "assay_type": oss_assay_type,
                "cleanup": oss_cleanup,
            }
            job_id = submit_job(config)
            st.success(f"Job submitted! ID: {job_id}")
            st.info("Processing continues even if you close the browser. Check progress in the Job Monitor tab.")
    else:
        st.warning("First retrieve the OSS file list.")


# =============================================================================
#  Tab 3: Upload
# =============================================================================

with tab_upload:
    st.subheader("1. Upload FASTQ files")
    st.caption("Upload local FASTQ files to run the pipeline.")

    uploaded_files = st.file_uploader(
        "FASTQ files (.fastq.gz / .fq.gz)",
        type=["gz"],
        accept_multiple_files=True,
        key="chip_upload_fastq",
    )

    if uploaded_files:
        fastq_files = [f for f in uploaded_files if f.name.endswith((".fastq.gz", ".fq.gz"))]
        other_files = [f for f in uploaded_files if not f.name.endswith((".fastq.gz", ".fq.gz"))]

        if other_files:
            st.warning(f"{len(other_files)} non-FASTQ files ignored: {', '.join(f.name for f in other_files)}")

        if fastq_files:
            paired = {}
            for f in fastq_files:
                name = f.name
                m = re.match(r'(.+?)_(?:L\d+_)?(?:R?1)\.(f(?:ast)?q\.gz)', name)
                if m:
                    key = m.group(1)
                    if key not in paired:
                        paired[key] = {"name": key, "r1": None, "r2": None}
                    paired[key]["r1"] = f
                    continue
                m = re.match(r'(.+?)_(?:L\d+_)?(?:R?2)\.(f(?:ast)?q\.gz)', name)
                if m:
                    key = m.group(1)
                    if key not in paired:
                        paired[key] = {"name": key, "r1": None, "r2": None}
                    paired[key]["r2"] = f
                    continue
                m = re.match(r'(.+?)_R1(.+)', name)
                if m:
                    key = m.group(1)
                    if key not in paired:
                        paired[key] = {"name": key, "r1": None, "r2": None}
                    paired[key]["r1"] = f
                    continue
                m = re.match(r'(.+?)_R2(.+)', name)
                if m:
                    key = m.group(1)
                    if key not in paired:
                        paired[key] = {"name": key, "r1": None, "r2": None}
                    paired[key]["r2"] = f
                    continue

            samples = list(paired.values())
            paired_samples = [s for s in samples if s["r1"] and s["r2"]]

            st.success(f"**{len(fastq_files)} files uploaded → {len(paired_samples)} paired-end samples**")
            for s in samples:
                status = "paired" if s["r1"] and s["r2"] else "INCOMPLETE"
                r1_name = s["r1"].name if s["r1"] else "?"
                r2_name = s["r2"].name if s["r2"] else "?"
                st.text(f"  {s['name']}: {status}  R1={r1_name}  R2={r2_name}")

            total_size = sum(f.size for f in fastq_files)
            st.info(f"Total size: **{format_size(total_size)}**")

            st.session_state.chip_upload_paired_samples = paired_samples
            st.session_state.chip_upload_fastq_files = fastq_files
        else:
            st.warning("No FASTQ files (.fastq.gz / .fq.gz) found.")

    st.markdown("---")

    st.subheader("2. Pipeline Configuration")
    up_genome, up_assay_type, up_local_dir, up_save_to_lab, up_include_homer, up_cleanup = pipeline_config_ui("chip_upload")

    st.markdown("---")

    st.subheader("3. Submit")

    has_upload_samples = len(st.session_state.get("chip_upload_paired_samples", [])) > 0

    if has_upload_samples:
        paired_samples = st.session_state.chip_upload_paired_samples
        st.markdown(f"**{len(paired_samples)} paired-end samples** / {ASSAY_LABELS[up_assay_type]} / {up_genome}")

        if st.button("Submit Job", type="primary", key="chip_upload_submit"):
            save_dir = up_local_dir
            os.makedirs(save_dir, exist_ok=True)

            with st.spinner("Saving uploaded files..."):
                for f in st.session_state.chip_upload_fastq_files:
                    local_path = os.path.join(save_dir, f.name)
                    with open(local_path, "wb") as out:
                        out.write(f.getbuffer())

            config = {
                "data_source": "local",
                "project_name": "upload",
                "local_dir": save_dir,
                "save_to_lab": up_save_to_lab,
                "include_homer_tags": up_include_homer,
                "genome": up_genome,
                "assay_type": up_assay_type,
                "cleanup": up_cleanup,
            }
            job_id = submit_job(config)
            st.success(f"Job submitted! ID: {job_id}")
            st.info("Processing continues even if you close the browser. Check progress in the Job Monitor tab.")
    else:
        st.warning("First upload FASTQ files.")


# =============================================================================
#  Tab 4: Job Monitor
# =============================================================================

with tab_monitor:
    st.subheader("Job Monitor")

    clear_old_directories(JOBS_DIR, date=7)

    auto_refresh = st.checkbox("Auto-refresh (10 sec)", value=False, key="chip_auto_refresh")

    if st.button("Refresh", key="chip_refresh"):
        st.rerun()

    @st.fragment(run_every=10 if auto_refresh else None)
    def _chip_monitor_fragment():
        jobs = get_all_jobs()

        if not jobs:
            st.info("No jobs found.")
        else:
            for job in jobs:
                s = job["status"]
                c = job["config"]
                status_str = s.get("status", "unknown")
                color = STATUS_COLORS.get(status_str, "gray")

                col_id, col_status, col_info = st.columns([2, 2, 4])
                with col_id:
                    project = c.get("project_name", "")
                    source = c.get("data_source", "sftp").upper()
                    label = f"**{project}** [{source}]" if project else f"**Job {job['name']}**"
                    st.markdown(label)
                with col_status:
                    st.markdown(f":{color}[{status_str.upper()}]")
                with col_info:
                    genome_str = c.get("genome", "?")
                    assay_str = ASSAY_LABELS.get(c.get("assay_type", ""), c.get("assay_type", "?"))
                    n_samples = len(c.get("sftp_samples", c.get("oss_samples", [])))
                    st.markdown(f"{genome_str} / {assay_str} / {n_samples} samples")

                progress = s.get("progress_pct", 0) / 100
                st.progress(progress, text=s.get("step", ""))

                with st.expander("Details"):
                    st.markdown(f"**Created:** {c.get('created_at', '?')}")
                    st.markdown(f"**Started:** {s.get('started_at', 'Not yet')}")
                    st.markdown(f"**Finished:** {s.get('finished_at', 'Not yet')}")
                    st.markdown(f"**Data source:** {c.get('data_source', 'sftp')}")
                    st.markdown(f"**Assay type:** {ASSAY_LABELS.get(c.get('assay_type', ''), '?')}")
                    st.markdown(f"**Local dir:** `{c.get('local_dir', '?')}`")

                    if s.get("error"):
                        st.error(f"Error: {s['error']}")

                    log_path = os.path.join(job["dir"], "pipeline.log")
                    if os.path.exists(log_path):
                        if st.button("Show log", key=f"chip_log_{job['name']}"):
                            with open(log_path) as f:
                                lines = f.readlines()
                                tail = lines[-50:] if len(lines) > 50 else lines
                                st.code("".join(tail), language="text")

                    if status_str == "done":
                        import zipfile
                        import io
                        import streamlit.components.v1 as components
                        import glob as glob_mod

                        save_to_lab = s.get("save_to_lab", False)

                        # Resolve file paths: filename -> full_path
                        file_map = {}
                        fastp_dir = ""

                        if save_to_lab:
                            base_dir = s.get("results_dir", c.get("local_dir", ""))
                            if base_dir and os.path.isdir(base_dir):
                                st.info(f"Results saved to: `{base_dir}`")
                                bam_dir = os.path.join(base_dir, "trimmed", "STAR_out", "bam")
                                seacr_dir = os.path.join(base_dir, "trimmed", "STAR_out", "seacr_input")
                                star_out_dir = os.path.join(base_dir, "trimmed", "STAR_out")
                                fastp_dir = os.path.join(base_dir, "fastp_reports")
                                if os.path.isdir(bam_dir):
                                    for f in os.listdir(bam_dir):
                                        if f.endswith((".dedup.bam", ".dedup.bam.bai", "_metrics.txt")):
                                            file_map[f] = os.path.join(bam_dir, f)
                                if os.path.isdir(seacr_dir):
                                    for f in os.listdir(seacr_dir):
                                        if f.endswith(".bedgraph"):
                                            file_map[f] = os.path.join(seacr_dir, f)
                                if os.path.isdir(star_out_dir):
                                    for f in os.listdir(star_out_dir):
                                        if f.endswith("Log.final.out"):
                                            file_map[f] = os.path.join(star_out_dir, f)
                        else:
                            results_dir = os.path.join(job["dir"], "results")
                            if os.path.isdir(results_dir):
                                for f in os.listdir(results_dir):
                                    fpath = os.path.join(results_dir, f)
                                    if os.path.isfile(fpath):
                                        file_map[f] = fpath
                                fastp_dir = os.path.join(results_dir, "fastp_reports")

                        if file_map:
                            # Categorize files
                            bam_files = sorted([f for f in file_map if f.endswith(".dedup.bam")])
                            bai_files = sorted([f for f in file_map if f.endswith(".dedup.bam.bai")])
                            metrics_files = sorted([f for f in file_map if f.endswith("_metrics.txt")])
                            bedgraph_files = sorted([f for f in file_map if f.endswith(".bedgraph")])
                            log_files = sorted([f for f in file_map if f.endswith("Log.final.out")])

                            # Create symlinks in static directory
                            static_dir = os.path.join(
                                os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "static", f"chip_{job['name']}",
                            )
                            os.makedirs(static_dir, exist_ok=True)
                            for fname, src_path in file_map.items():
                                dst = os.path.join(static_dir, fname)
                                if os.path.isfile(src_path) and not os.path.exists(dst):
                                    os.symlink(src_path, dst)

                            def static_url(fname):
                                return f"/app/static/chip_{job['name']}/{fname}"

                            def download_all_button(label, file_list, btn_key):
                                """Generate a button that downloads files sequentially via JS"""
                                files_js = json.dumps([
                                    {"name": f, "url": static_url(f)}
                                    for f in file_list
                                ])
                                total_size = sum(
                                    os.path.getsize(file_map[f])
                                    for f in file_list if f in file_map
                                )
                                components.html(f"""
                                <button id="btn_{btn_key}" onclick="downloadAll_{btn_key}()" style="
                                    background-color:#ff4b4b; color:white; border:none;
                                    padding:8px 16px; border-radius:6px; cursor:pointer;
                                    font-size:14px; font-weight:500;">
                                    {label} ({format_size(total_size)})
                                </button>
                                <span id="status_{btn_key}" style="margin-left:10px; font-size:13px; color:#666;"></span>
                                <script>
                                async function downloadAll_{btn_key}() {{
                                    const btn = document.getElementById('btn_{btn_key}');
                                    if (btn.disabled) return;
                                    btn.disabled = true;
                                    btn.style.opacity = '0.5';
                                    btn.style.cursor = 'not-allowed';
                                    const files = {files_js};
                                    const status = document.getElementById('status_{btn_key}');
                                    status.textContent = 'If the browser asks for permission to download multiple files, click "Allow".';
                                    await new Promise(r => setTimeout(r, 1000));
                                    for (let i = 0; i < files.length; i++) {{
                                        status.textContent = `Downloading ${{i+1}}/${{files.length}}: ${{files[i].name}}`;
                                        const a = document.createElement('a');
                                        a.href = files[i].url;
                                        a.download = files[i].name;
                                        document.body.appendChild(a);
                                        a.click();
                                        document.body.removeChild(a);
                                        if (i < files.length - 1) {{
                                            await new Promise(r => setTimeout(r, 3000));
                                        }}
                                    }}
                                    status.textContent = `${{files.length}} files queued. Large files (BAM etc.) may take several minutes to complete.`;
                                    btn.disabled = false;
                                    btn.style.opacity = '1';
                                    btn.style.cursor = 'pointer';
                                }}
                                </script>
                                """, height=45)

                            # --- BAM + BAI ---
                            if bam_files:
                                bam_bai = []
                                for f in bam_files:
                                    bam_bai.append(f)
                                    bai = f + ".bai"
                                    if bai in bai_files:
                                        bam_bai.append(bai)
                                st.markdown(f"**BAM + index** ({len(bam_files)} samples)")
                                download_all_button(
                                    f"Download All BAM + BAI",
                                    bam_bai, f"bam_{job['name']}",
                                )
                                with st.expander("Individual downloads"):
                                    for f in bam_files:
                                        size = format_size(os.path.getsize(file_map[f])) if f in file_map else ""
                                        bai_name = f + ".bai"
                                        bai_link = f" | [.bai]({static_url(bai_name)})" if bai_name in bai_files else ""
                                        st.markdown(
                                            f"[{f}]({static_url(f)}) ({size}){bai_link}",
                                            unsafe_allow_html=False,
                                        )

                            # --- bedGraph files ---
                            if bedgraph_files:
                                st.markdown(f"**Fragment bedGraph** ({len(bedgraph_files)} files)")
                                download_all_button(
                                    f"Download All bedGraph",
                                    bedgraph_files, f"bg_{job['name']}",
                                )

                            # --- Dedup metrics (zipped) ---
                            if metrics_files:
                                st.markdown(f"**Dedup metrics** ({len(metrics_files)} files)")
                                buf = io.BytesIO()
                                with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                                    for fname in metrics_files:
                                        zf.write(file_map[fname], fname)
                                st.download_button(
                                    f"dedup_metrics.zip ({format_size(buf.tell())})",
                                    data=buf.getvalue(),
                                    file_name="dedup_metrics.zip",
                                    key=f"chip_dl_{job['name']}_metrics_zip",
                                )

                            # --- STAR logs (zipped) ---
                            if log_files:
                                st.markdown(f"**STAR logs** ({len(log_files)} files)")
                                buf = io.BytesIO()
                                with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                                    for fname in log_files:
                                        zf.write(file_map[fname], fname)
                                st.download_button(
                                    f"star_logs.zip ({format_size(buf.tell())})",
                                    data=buf.getvalue(),
                                    file_name="star_logs.zip",
                                    key=f"chip_dl_{job['name']}_logs_zip",
                                )

                            # --- fastp reports (HTML only, zipped) ---
                            if fastp_dir and os.path.isdir(fastp_dir):
                                html_files = sorted(
                                    f for f in os.listdir(fastp_dir) if f.endswith(".html")
                                )
                                if html_files:
                                    st.markdown(f"**fastp reports** ({len(html_files)} HTML files)")
                                    buf = io.BytesIO()
                                    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
                                        for hf in html_files:
                                            zf.write(os.path.join(fastp_dir, hf), hf)
                                    st.download_button(
                                        f"fastp_reports.zip ({format_size(buf.tell())})",
                                        data=buf.getvalue(),
                                        file_name="fastp_reports.zip",
                                        key=f"chip_dl_{job['name']}_fastp_zip",
                                    )

                st.markdown("---")

    _chip_monitor_fragment()
