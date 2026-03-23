"""
FASTQ -> RNA-seq count table

Data sources:
  1. SFTP server: Navigate run directories to select FASTQ files
  2. Azenta OSS: Download from Alibaba Cloud OSS via ossutil
  3. Upload: Upload local FASTQ files

Runs STAR/Salmon pipeline in background to generate count tables.
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

# =============================================================================
#  Constants
# =============================================================================

_STREAMLIT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
JOBS_DIR = os.path.join(_STREAMLIT_DIR, "temp", "fastq2rna_jobs")
WORKER_SCRIPT = os.path.join(_STREAMLIT_DIR, "fastq2rna_worker.py")
PYTHON_PATH = "python3"

GENOME_OPTIONS = ["mm10", "hg38", "mm39"]
COUNT_TOOLS = ["salmon", "featurecounts", "homer"]

STATUS_COLORS = {
    "queued": "gray",
    "downloading": "blue",
    "trimming": "orange",
    "mapping": "orange",
    "counting": "violet",
    "postprocessing": "violet",
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
    """Scan the run directory tree focusing on FASTQ files

    Returns: list of dicts:
        {"project": str, "samples": [{"name": str, "files": [(name, size)]}]}
    """
    projects = []
    try:
        # Look for SampleSheet/
        top_items = sftp.listdir_attr(path)
        top_dirs = {i.filename: i for i in top_items if stat_module.S_ISDIR(i.st_mode)}

        # If SampleSheet exists, explore inside it; otherwise explore directly
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

            # Look for Other/ subdirectory
            other_path = f"{proj_path}/Other"
            try:
                sample_items = sftp.listdir_attr(other_path)
            except Exception:
                # If Other/ doesn't exist, look directly
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
            fastq_info = "  ".join(
                f"{fname} ({format_size(sz)})" for fname, sz in sample["files"]
            )
            lines.append(f"    📁 {sample['name']}/")
            for fname, sz in sample["files"]:
                lines.append(f"        🧬 {fname}  ({format_size(sz)})")
    return "\n".join(lines)


def scan_samples(sftp, project_dir):
    """Detect samples within a project directory

    Structure: R1/R2 fastq.gz files in {project}/Other/{sample}/
    Also handles fastq.gz files directly under {project}/
    """
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
        # Line format: date time tz_offset tz_name  size  StorageClass  ETAG  oss://path
        # Example: 2026-02-13 17:10:11 +0900 JST   1494546648   Standard   ETAG   oss://...
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

    # Pairing: _1.fq.gz / _2.fq.gz or _R1 / _R2
    paired = {}
    for f in fastq_files:
        name = f["name"]
        # _1.fq.gz / _1.fastq.gz pattern
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
        # _R1 / _R2 pattern
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
    """Pipeline configuration UI (common section)"""
    col1, col2 = st.columns(2)
    with col1:
        genome = st.selectbox("Genome", GENOME_OPTIONS, index=0, key=f"{key_prefix}_genome")
        count_tool = st.selectbox("Count tool", COUNT_TOOLS, index=0, key=f"{key_prefix}_count_tool")
    with col2:
        annotation = "refseq"
        strandedness = 0
        convert_bam = False
        make_bigwig = False
        if count_tool != "salmon":
            default_ann_idx = 1 if count_tool == "homer" else 0
            annotation = st.selectbox(
                "Annotation", ["gencode", "refseq"],
                index=default_ann_idx,
                key=f"{key_prefix}_annotation",
                help="Unify STAR index and GTF. "
                     "RefSeq is the default for Homer.")
        if count_tool == "featurecounts":
            strandedness = st.selectbox(
                "Strandedness", [0, 1, 2],
                format_func=lambda x: {0: "0: unstranded", 1: "1: stranded", 2: "2: reverse"}[x],
                key=f"{key_prefix}_strandedness",
            )
        if count_tool != "salmon":
            convert_bam = st.checkbox("Convert to BAM", key=f"{key_prefix}_bam")
            if convert_bam:
                make_bigwig = st.checkbox("Generate bigWig", key=f"{key_prefix}_bigwig")

    if count_tool == "salmon":
        st.info("Salmon: No STAR needed. fastp -> Salmon quant -> gene-level count table.")

    # --- Save directory selection ---
    st.markdown("**Save directory**")
    lab_ngs_available = os.path.isdir(LAB_NGS_BASE)

    if lab_ngs_available:
        save_to_lab = st.checkbox(
            "Save to shared storage",
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

        rel = os.path.relpath(current_path, LAB_NGS_BASE)
        if rel == ".":
            st.markdown("📁 **shared_storage/**")
        else:
            st.markdown(f"📁 **shared_storage/{rel}/**")

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
                placeholder="e.g. 2026-02_RNA-seq",
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

    cleanup = st.checkbox(
        "Cleanup after completion (delete FASTQ and intermediate files)",
        value=True,
        key=f"{key_prefix}_cleanup",
    )

    return genome, count_tool, annotation, strandedness, convert_bam, make_bigwig, local_dir, save_to_lab, cleanup


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

st.title("FASTQ to RNA-seq Count Table")

tab_sftp, tab_oss, tab_upload, tab_monitor = st.tabs(["Aging SFTP", "Azenta OSS", "Upload", "Job Monitor"])

# =============================================================================
#  Tab 1: Aging SFTP
# =============================================================================

with tab_sftp:
    st.subheader("1. Select FASTQ files from Aging server")
    st.caption("Navigate: Nova-seq → {run} → SampleSheet → **{project}** → Scan")
    st.caption("Each project within SampleSheet (e.g. Aoki_RNA-seq, Kim_RNA-seq) should be submitted as a separate job.")

    # SFTP browser
    remote_dir = st.text_input(
        "Remote directory path",
        value=st.session_state.get("sftp_remote_dir", "/data/sequencing"),
        key="remote_dir_input",
    )

    if st.button("Browse", type="primary"):
        st.session_state.sftp_remote_dir = remote_dir

    if "sftp_remote_dir" in st.session_state:
        ssh, sftp = get_sftp_connection()
        if sftp is not None:
            current_dir = st.session_state.sftp_remote_dir
            st.caption(f"Current: `{current_dir}`")
            dirs, files = list_remote_dir(sftp, current_dir)

            # Subdirectory listing (run selection)
            if dirs:
                selected_run = st.selectbox(
                    "Run directory",
                    ["-- select --"] + dirs,
                    key="subdir_nav",
                )

                # Auto-scan tree when selection changes
                if selected_run != "-- select --":
                    prev = st.session_state.get("_prev_subdir_nav")
                    if prev != selected_run:
                        st.session_state._prev_subdir_nav = selected_run
                        target = f"{current_dir}/{selected_run}"
                        with st.spinner("Scanning directory tree..."):
                            projects = scan_directory_tree(sftp, target)
                            st.session_state.sftp_dir_tree = (selected_run, target, projects)

                # Directory tree display + project selection
                if "sftp_dir_tree" in st.session_state and len(st.session_state.sftp_dir_tree) == 3:
                    tree_name, tree_base, tree_projects = st.session_state.sftp_dir_tree
                    if tree_projects:
                        with st.expander(f"**{tree_name}/ project structure**", expanded=True):
                            st.code(format_project_tree(tree_projects), language="text")

                        # Project selection below the tree
                        proj_names = [p["project"] for p in tree_projects]
                        selected_proj = st.selectbox(
                            "Project to process",
                            ["-- select --"] + proj_names,
                            key="proj_nav",
                        )
                        if selected_proj != "-- select --":
                            # Build path via SampleSheet/
                            proj_path = f"{tree_base}/SampleSheet/{selected_proj}"
                            try:
                                sftp.stat(proj_path)
                            except Exception:
                                proj_path = f"{tree_base}/{selected_proj}"

                            # Auto-scan when project is selected
                            prev_proj = st.session_state.get("_prev_proj_nav")
                            if prev_proj != selected_proj:
                                st.session_state._prev_proj_nav = selected_proj
                                st.session_state.sftp_remote_dir = proj_path
                                with st.spinner("Scanning samples..."):
                                    samples = scan_samples(sftp, proj_path)
                                    st.session_state.sftp_samples = samples

                            # Display samples
                            if "sftp_samples" in st.session_state and st.session_state.sftp_samples:
                                samples = st.session_state.sftp_samples
                                st.success(f"**{len(samples)} samples found**")
                                for s in samples:
                                    paired = "paired" if s['r1'] and s['r2'] else "INCOMPLETE"
                                    st.text(f"  {s['name']}: {paired}  R1={format_size(s['r1'][1]) if s['r1'] else '?'} R2={format_size(s['r2'][1]) if s['r2'] else '?'}")
                                total_size = sum(f[1] for s in samples for f in s["files"])
                                st.info(f"Total download size: **{format_size(total_size)}**")
                            elif "sftp_samples" in st.session_state:
                                st.warning("No FASTQ samples found in this project.")

                    elif tree_name:
                        st.warning(f"No FASTQ samples found in {tree_name}/.")

    st.markdown("---")

    # Pipeline configuration
    st.subheader("2. Pipeline Configuration")
    sftp_genome, sftp_count_tool, sftp_annotation, sftp_strandedness, sftp_convert_bam, sftp_make_bigwig, sftp_local_dir, sftp_save_to_lab, sftp_cleanup = pipeline_config_ui("sftp", default_save_to_lab=True)

    st.markdown("---")

    # Submit
    st.subheader("3. Submit")

    has_samples = "sftp_samples" in st.session_state and len(st.session_state.get("sftp_samples", [])) > 0

    if has_samples:
        samples = st.session_state.sftp_samples
        paired_samples = [s for s in samples if s["r1"] and s["r2"]]
        project_name = os.path.basename(st.session_state.get("sftp_remote_dir", "").rstrip("/"))
        st.markdown(f"**Project: {project_name}** — {len(paired_samples)} paired-end samples / {sftp_count_tool} / {sftp_genome}")

        if not paired_samples:
            st.error("No paired-end samples found.")
        elif st.button("Submit Job", type="primary", key="sftp_submit"):
            cred = st.secrets["aging_sftp"]
            project_name = os.path.basename(st.session_state.sftp_remote_dir.rstrip("/"))
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
                "genome": sftp_genome,
                "count_tool": sftp_count_tool,
                "annotation": sftp_annotation,
                "strandedness": sftp_strandedness,
                "convert_bam": sftp_convert_bam,
                "make_bigwig": sftp_make_bigwig,
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
    st.caption("Enter the OSS information provided in the GENEWIZ email.")

    oss_path = st.text_input(
        "OSS Path (Add to OSS)",
        value=st.session_state.get("oss_path", ""),
        key="oss_path_input",
        placeholder="oss://gwjapan/2026.2/60-1286302729/N2602837_60-1286302729_20260213152423",
    )

    col_oss1, col_oss2 = st.columns(2)
    with col_oss1:
        oss_key_id = st.text_input("AccessKeyId", key="oss_key_id")
    with col_oss2:
        oss_key_secret = st.text_input("AccessKeySecret", type="password", key="oss_key_secret")

    oss_region = st.selectbox(
        "Region",
        ["oss-ap-northeast-1 (Tokyo)", "oss-ap-southeast-1 (Singapore)", "oss-cn-shanghai"],
        index=0,
        key="oss_region",
    )
    oss_endpoint = oss_region.split(" ")[0] + ".aliyuncs.com"

    # OSS file listing
    if st.button("List files", type="primary", key="oss_list"):
        if not oss_key_id or not oss_key_secret:
            st.error("Please enter AccessKeyId and AccessKeySecret.")
        else:
            with st.spinner("Listing OSS files..."):
                # Append / to the end of the path
                path = oss_path.rstrip("/") + "/"
                files, err = list_oss_files(path, oss_key_id, oss_key_secret, oss_endpoint)
                if err:
                    st.error(f"OSS error: {err}")
                elif files:
                    st.session_state.oss_files = files
                    st.session_state.oss_path = oss_path
                else:
                    st.warning("No files found at this path.")

    # File listing display and sample detection
    if "oss_files" in st.session_state and st.session_state.oss_files:
        files = st.session_state.oss_files
        fastq_files = [f for f in files if f["name"].endswith((".fastq.gz", ".fq.gz"))]
        other_files = [f for f in files if not f["name"].endswith((".fastq.gz", ".fq.gz"))]

        st.success(f"**{len(fastq_files)} FASTQ files** / {len(other_files)} other files")

        # Sample detection
        samples = parse_oss_samples(files)
        st.session_state.oss_samples = samples

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

        # Other files
        if other_files:
            with st.expander(f"Other files ({len(other_files)})"):
                for f in other_files:
                    st.text(f"  {f['name']} ({format_size(f['size'])})")

    st.markdown("---")

    # Pipeline configuration
    st.subheader("2. Pipeline Configuration")
    oss_genome, oss_count_tool, oss_annotation, oss_strandedness, oss_convert_bam, oss_make_bigwig, oss_local_dir, oss_save_to_lab, oss_cleanup = pipeline_config_ui("oss")

    st.markdown("---")

    # Submit
    st.subheader("3. Submit")

    has_oss_samples = "oss_samples" in st.session_state and len(st.session_state.get("oss_samples", [])) > 0

    if has_oss_samples:
        samples = st.session_state.oss_samples
        paired_samples = [s for s in samples if s["r1"] and s["r2"]]
        # Infer project name from OSS path
        oss_project = os.path.basename(st.session_state.get("oss_path", "").rstrip("/"))
        st.markdown(f"**Project: {oss_project}** — {len(paired_samples)} paired-end samples / {oss_count_tool} / {oss_genome}")

        if not paired_samples:
            st.error("No paired-end samples found.")
        elif not oss_key_id or not oss_key_secret:
            st.error("Please enter AccessKeyId and AccessKeySecret.")
        elif st.button("Submit Job", type="primary", key="oss_submit"):
            # OSS sample information
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
                "genome": oss_genome,
                "count_tool": oss_count_tool,
                "annotation": oss_annotation,
                "strandedness": oss_strandedness,
                "convert_bam": oss_convert_bam,
                "make_bigwig": oss_make_bigwig,
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
        key="upload_fastq",
    )

    if uploaded_files:
        # Pairing detection
        fastq_files = [f for f in uploaded_files if f.name.endswith((".fastq.gz", ".fq.gz"))]
        other_files = [f for f in uploaded_files if not f.name.endswith((".fastq.gz", ".fq.gz"))]

        if other_files:
            st.warning(f"{len(other_files)} non-FASTQ files ignored: {', '.join(f.name for f in other_files)}")

        if fastq_files:
            # Pairing
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

            st.session_state.upload_paired_samples = paired_samples
            st.session_state.upload_fastq_files = fastq_files
        else:
            st.warning("No FASTQ files (.fastq.gz / .fq.gz) found.")

    st.markdown("---")

    # Pipeline configuration
    st.subheader("2. Pipeline Configuration")
    up_genome, up_count_tool, up_annotation, up_strandedness, up_convert_bam, up_make_bigwig, up_local_dir, up_save_to_lab, up_cleanup = pipeline_config_ui("upload")

    st.markdown("---")

    # Submit
    st.subheader("3. Submit")

    has_upload_samples = len(st.session_state.get("upload_paired_samples", [])) > 0

    if has_upload_samples:
        paired_samples = st.session_state.upload_paired_samples
        st.markdown(f"**{len(paired_samples)} paired-end samples** / {up_count_tool} / {up_genome}")

        if st.button("Submit Job", type="primary", key="upload_submit"):
            # Save files locally
            save_dir = up_local_dir
            os.makedirs(save_dir, exist_ok=True)

            with st.spinner("Saving uploaded files..."):
                for f in st.session_state.upload_fastq_files:
                    local_path = os.path.join(save_dir, f.name)
                    with open(local_path, "wb") as out:
                        out.write(f.getbuffer())

            config = {
                "data_source": "local",
                "project_name": "upload",
                "local_dir": save_dir,
                "save_to_lab": up_save_to_lab,
                "genome": up_genome,
                "count_tool": up_count_tool,
                "annotation": up_annotation,
                "strandedness": up_strandedness,
                "convert_bam": up_convert_bam,
                "make_bigwig": up_make_bigwig,
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

    auto_refresh = st.checkbox("Auto-refresh (10 sec)", value=False)

    if st.button("Refresh"):
        st.rerun()

    @st.fragment(run_every=10 if auto_refresh else None)
    def _rna_monitor_fragment():
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
                    tool_str = c.get("count_tool", "?")
                    n_samples = len(c.get("sftp_samples", c.get("oss_samples", [])))
                    st.markdown(f"{genome_str} / {tool_str} / {n_samples} samples")

                progress = s.get("progress_pct", 0) / 100
                st.progress(progress, text=s.get("step", ""))

                with st.expander("Details"):
                    st.markdown(f"**Created:** {c.get('created_at', '?')}")
                    st.markdown(f"**Started:** {s.get('started_at', 'Not yet')}")
                    st.markdown(f"**Finished:** {s.get('finished_at', 'Not yet')}")
                    st.markdown(f"**Data source:** {c.get('data_source', 'sftp')}")
                    st.markdown(f"**Local dir:** `{c.get('local_dir', '?')}`")

                    if s.get("error"):
                        st.error(f"Error: {s['error']}")

                    log_path = os.path.join(job["dir"], "pipeline.log")
                    if os.path.exists(log_path):
                        if st.button("Show log", key=f"log_{job['name']}"):
                            with open(log_path) as f:
                                lines = f.readlines()
                                tail = lines[-50:] if len(lines) > 50 else lines
                                st.code("".join(tail), language="text")

                    if status_str == "done":
                        import zipfile
                        import io

                        save_to_lab = s.get("save_to_lab", False)

                        # Resolve file paths: filename -> full_path
                        file_map = {}
                        fastp_dir = ""

                        if save_to_lab:
                            base_dir = s.get("results_dir", c.get("local_dir", ""))
                            if base_dir and os.path.isdir(base_dir):
                                st.info(f"Results saved to: `{base_dir}`")
                                trimmed_dir = os.path.join(base_dir, "trimmed")
                                star_out_dir = os.path.join(trimmed_dir, "STAR_out")
                                salmon_out_dir = os.path.join(trimmed_dir, "salmon_out")
                                fastp_dir = os.path.join(base_dir, "fastp_reports")
                                # Count table
                                for ct_dir, ct_name in [
                                    (salmon_out_dir, "counts_salmon.tsv"),
                                    (star_out_dir, "counts.txt"),
                                    (os.path.join(star_out_dir, "Homer"), "Raw.txt"),
                                ]:
                                    ct_path = os.path.join(ct_dir, ct_name)
                                    if os.path.isfile(ct_path):
                                        file_map[ct_name] = ct_path
                                # BAM
                                bam_dir = os.path.join(star_out_dir, "bam")
                                if os.path.isdir(bam_dir):
                                    for f in os.listdir(bam_dir):
                                        if f.endswith((".sorted.bam", ".sorted.bam.bai")):
                                            file_map[f] = os.path.join(bam_dir, f)
                                # bigWig
                                bigwig_dir = os.path.join(star_out_dir, "bigwig")
                                if os.path.isdir(bigwig_dir):
                                    for f in os.listdir(bigwig_dir):
                                        if f.endswith(".bw"):
                                            file_map[f] = os.path.join(bigwig_dir, f)
                                # STAR logs
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
                            st.markdown("**Results:**")
                            # Count table (small - direct download)
                            count_files = [f for f in sorted(file_map)
                                           if f in ("counts_salmon.tsv", "counts.txt", "Raw.txt")]
                            for fname in count_files:
                                with open(file_map[fname], "rb") as fh:
                                    st.download_button(
                                        f"📊 {fname}",
                                        data=fh.read(),
                                        file_name=fname,
                                        key=f"dl_{job['name']}_{fname}",
                                    )

                            # BAM files
                            bam_files = sorted([f for f in file_map if f.endswith(".sorted.bam")])
                            if bam_files:
                                st.markdown(f"**BAM** ({len(bam_files)} files)")
                                for fname in bam_files:
                                    with open(file_map[fname], "rb") as fh:
                                        st.download_button(
                                            f"Download {fname} ({format_size(os.path.getsize(file_map[fname]))})",
                                            data=fh.read(),
                                            file_name=fname,
                                            key=f"dl_{job['name']}_{fname}",
                                        )

                            # bigWig files
                            bw_files = sorted([f for f in file_map if f.endswith(".bw")])
                            if bw_files:
                                st.markdown(f"**bigWig** ({len(bw_files)} files)")
                                for fname in bw_files:
                                    with open(file_map[fname], "rb") as fh:
                                        st.download_button(
                                            f"Download {fname}",
                                            data=fh.read(),
                                            file_name=fname,
                                            key=f"dl_{job['name']}_{fname}",
                                        )

                            # STAR logs (zipped)
                            log_files = sorted([f for f in file_map if f.endswith("Log.final.out")])
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
                                    key=f"dl_{job['name']}_logs_zip",
                                )

                            # fastp reports (HTML only, zipped)
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
                                        key=f"dl_{job['name']}_fastp_zip",
                                    )

                st.markdown("---")

    _rna_monitor_fragment()
