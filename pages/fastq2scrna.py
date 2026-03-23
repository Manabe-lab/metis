"""
FASTQ -> scRNA-seq (Cell Ranger)

Data sources:
  1. Local NovaSeq (Aging SFTP): BCL run -> mkfastq -> count/multi
  2. Azenta (OSS/SFTP/Upload): FASTQ -> rename -> count/multi

Runs the Cell Ranger pipeline in the background and
saves results to the SCALA personal folder in 10X format.
"""

import streamlit as st
import paramiko
import os
from dotenv import load_dotenv

load_dotenv()
import json
import time
import shutil
import subprocess
import stat as stat_module
import re
import pandas as pd
from pathlib import Path
from datetime import datetime

# =============================================================================
#  Constants
# =============================================================================

_STREAMLIT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_HOME = str(Path.home())
JOBS_DIR = os.path.join(_STREAMLIT_DIR, "temp", "fastq2scrna_jobs")
WORKER_SCRIPT = os.path.join(_STREAMLIT_DIR, "fastq2scrna_worker.py")
PYTHON_PATH = "python3"
CELLRANGER_PATH = os.environ.get("CELLRANGER_PATH", "cellranger")

_CELLRANGER_REF_DIR = os.environ.get(
    "CELLRANGER_REF_DIR", os.path.expanduser("~/cellranger_refs")
)
_ALL_REFERENCES = {
    "GRCm39-2024-A": os.path.join(_CELLRANGER_REF_DIR, "refdata-gex-GRCm39-2024-A"),
    "mm10-2020A": os.path.join(_CELLRANGER_REF_DIR, "mm10-2020A"),
    "GRCh38-2024-A": os.path.join(_CELLRANGER_REF_DIR, "refdata-gex-GRCh38-2024-A"),
    "GRCh38-2020-A": os.path.join(_CELLRANGER_REF_DIR, "refdata-gex-GRCh38-2020-A"),
    "mm10-2020A+EGFP": os.path.join(_CELLRANGER_REF_DIR, "mm10-2020AplusEGFP"),
    "mm10+LncFAO": os.path.join(_CELLRANGER_REF_DIR, "mm10plusLncFAO"),
    "mm10+LncFAO-EGFP": os.path.join(_CELLRANGER_REF_DIR, "mm10plusLncFAO-EGFP"),
    "mm10-2020A+EGFP (new)": os.path.join(_CELLRANGER_REF_DIR, "build_EGFP_ref/mm10-2020A-plusEGFP"),
    "GRCm39-2024-A+EGFP": os.path.join(_CELLRANGER_REF_DIR, "build_EGFP_ref/mm39-2024A-plusEGFP"),
}
# Show only references that exist on disk
REFERENCES = {k: v for k, v in _ALL_REFERENCES.items() if os.path.isdir(v)}

SCALA_BASE = os.environ.get(
    "SCALA_DATA_DIR",
    os.path.expanduser("~/scala_data/Personal_folders"),
)

# Master hashtag database
HASHTAG_DB_PATH = os.path.join(_STREAMLIT_DIR, "data", "totalseq_hashtag_barcodes.csv")


@st.cache_data
def load_hashtag_db():
    """Load TotalSeq hashtag master data"""
    if os.path.isfile(HASHTAG_DB_PATH):
        return pd.read_csv(HASHTAG_DB_PATH)
    return pd.DataFrame()

# =============================================================================
#  FASTQ prefix -> library type auto-classification
# =============================================================================

_LIBRARY_PATTERNS = {
    "Gene Expression": re.compile(
        r'(?i)(^|[_\-])(GEX|GE|RNA|GeneExp(ression)?|Transcriptome)([_\-]|$)'
    ),
    "Antibody Capture": re.compile(
        r'(?i)(^|[_\-])(SP|ADT|AB|Antibody|Surface[\s_\-]?Protein|CITE|Protein|TotalSeq|FeatureBarcod(e|ing))([_\-]|$)'
    ),
    "Multiplexing Capture": re.compile(
        r'(?i)(^|[_\-])(CMO|HTO|Hashtag|Hash|Multiplex|CellPlex)([_\-]|$)'
    ),
}


def classify_fastq_prefix(prefix):
    """Infer library type from a FASTQ prefix. Returns None if unknown."""
    for lib_type, pattern in _LIBRARY_PATTERNS.items():
        if pattern.search(prefix):
            return lib_type
    return None


def classify_all_prefixes(prefixes, llm_classification=None):
    """Classify all prefixes. LLM results take priority; otherwise use pattern matching.
    Unclassified prefixes default to Gene Expression."""
    result = {}
    for p in prefixes:
        if llm_classification and p in llm_classification:
            result[p] = llm_classification[p]
        else:
            cls = classify_fastq_prefix(p)
            result[p] = cls if cls else "Gene Expression"
    return result


STATUS_COLORS = {
    "queued": "gray",
    "downloading": "blue",
    "download_complete": "green",
    "renaming": "orange",
    "mkfastq": "orange",
    "cellranger": "violet",
    "copying_results": "violet",
    "done": "green",
    "error": "red",
}

os.makedirs(JOBS_DIR, exist_ok=True)

# =============================================================================
#  SFTP connection (Aging) -- same pattern as fastq2rna.py
# =============================================================================

def get_sftp_connection():
    """Get SFTP connection (cached in session_state)"""
    if "scrna_sftp_ssh" in st.session_state and st.session_state.scrna_sftp_ssh is not None:
        try:
            transport = st.session_state.scrna_sftp_ssh.get_transport()
            if transport is not None and transport.is_active():
                return st.session_state.scrna_sftp_ssh, st.session_state.scrna_sftp_client
        except Exception:
            pass

    try:
        cred = st.secrets["aging_sftp"]
        ssh = paramiko.SSHClient()
        ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        ssh.connect(cred["host"], username=cred["username"], password=cred["password"])
        sftp = ssh.open_sftp()
        st.session_state.scrna_sftp_ssh = ssh
        st.session_state.scrna_sftp_client = sftp
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


# =============================================================================
#  Azenta OSS utilities -- same as fastq2rna.py
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


# =============================================================================
#  Common utilities
# =============================================================================

def format_size(size_bytes):
    """Format file size into a human-readable string"""
    for unit in ["B", "KB", "MB", "GB", "TB"]:
        if size_bytes < 1024:
            return f"{size_bytes:.1f} {unit}"
        size_bytes /= 1024
    return f"{size_bytes:.1f} PB"


def _build_dir_tree(path, prefix="", max_depth=3, depth=0):
    """Generate a tree string showing directories only"""
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


# =============================================================================
#  SCALA output directory browser
# =============================================================================

def scala_dir_browser(key_prefix):
    """Directory browser for selecting a SCALA personal folder"""
    st.markdown("**SCALA output directory**")

    scala_available = os.path.isdir(SCALA_BASE)
    if not scala_available:
        st.warning(f"SCALA directory not accessible: `{SCALA_BASE}`")
        scala_dest = st.text_input(
            "SCALA destination path",
            value=SCALA_BASE + "/",
            key=f"{key_prefix}_scala_dest_text",
        )
        return scala_dest

    browse_key = f"{key_prefix}_scala_path"
    if browse_key not in st.session_state:
        st.session_state[browse_key] = SCALA_BASE

    current_path = st.session_state[browse_key]

    rel = os.path.relpath(current_path, SCALA_BASE)
    if rel == ".":
        st.markdown("📁 **Personal_folders/**")
    else:
        st.markdown(f"📁 **Personal_folders/{rel}/**")

    # Folder tree display
    if current_path != SCALA_BASE:
        tree_lines = _build_dir_tree(current_path)
        if tree_lines:
            rel_name = os.path.basename(current_path)
            tree_text = f"{rel_name}/\n" + "\n".join(tree_lines)
            with st.expander("Folder structure", expanded=True):
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
                key=f"{key_prefix}_scala_subdir_sel",
            )
        else:
            selected_sub = None
            st.caption("(empty)")
    with col_nav2:
        if selected_sub and st.button("Open", key=f"{key_prefix}_scala_open_dir"):
            st.session_state[browse_key] = os.path.join(current_path, selected_sub)
            st.rerun()
    with col_nav3:
        if current_path != SCALA_BASE:
            if st.button("Back", key=f"{key_prefix}_scala_up_dir"):
                st.session_state[browse_key] = os.path.dirname(current_path)
                st.rerun()

    # Create new folder
    col_new1, col_new2 = st.columns([4, 1])
    with col_new1:
        new_folder_name = st.text_input(
            "New folder",
            key=f"{key_prefix}_scala_new_folder",
            placeholder="e.g. project_scRNA",
            label_visibility="collapsed",
        )
    with col_new2:
        if st.button("+ Create", key=f"{key_prefix}_scala_create_folder"):
            if new_folder_name:
                new_path = os.path.join(current_path, new_folder_name)
                os.makedirs(new_path, exist_ok=True)
                st.session_state[browse_key] = new_path
                st.rerun()

    st.success(f"SCALA output: `{current_path}`")
    return current_path


# =============================================================================
#  Sample information extraction via LLM
# =============================================================================

def extract_sample_info_with_llm(file_content, file_name="", extraction_type="multi",
                                  fastq_filenames=None):
    """Extract sample information from an uploaded file using an LLM.

    extraction_type:
        "multi" -- hashtag, libraries, samples (for multi mode)
        "mkfastq" -- GEX/CMO sample sheet (for mkfastq)
        "local" -- integrated mkfastq + multi extraction (for Local NovaSeq)
    fastq_filenames: list of FASTQ filenames (available for Local Folder, etc.)
    """
    try:
        import anthropic
    except ImportError:
        st.error("The anthropic package is not installed. Run `pip install anthropic`.")
        return None

    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        st.error("ANTHROPIC_API_KEY is not set.")
        return None

    if extraction_type == "mkfastq":
        prompt = f"""Extract the sample sheet information needed for cellranger mkfastq (10x Chromium scRNA-seq) from the following file content.

File name: {file_name}

File content:
```
{file_content}
```

This file typically contains the following information:
- Sample names (e.g., PB, LSK, SampleGEX, etc.)
- Dual Index TT Set A information (e.g., H12 -> SI-TT-H12) -- for Gene Expression libraries
- Dual Index NT/NN Set A information (e.g., H10 -> SI-NT-H10) -- for Cell Surface Protein / Hashtag libraries
- Run numbers or sheet distinctions

Return in the following JSON format:

```json
{{
  "gex_sample_sheet": [
    {{
      "lane": "*",
      "sample": "sample_name",
      "index": "SI-TT-XXX (10x Dual Index Kit TT Set A)"
    }}
  ],
  "cmo_sample_sheet": [
    {{
      "lane": "*",
      "sample": "Hashtag-sample_name",
      "index": "SI-NT-XXX (10x Dual Index Kit NT/NN Set A)"
    }}
  ]
}}
```

Notes:
- Indexes must include the full prefix "SI-TT-" or "SI-NT-".
  Example: "H12" -> "SI-TT-H12", "A1" -> "SI-TT-A1"
- If there is no CMO / Hashtag library, set cmo_sample_sheet to an empty list [].
- If there are multiple runs, combine all samples and return them together.
- Return JSON only. No explanatory text.
"""
    elif extraction_type == "local":
        prompt = f"""Extract **all** information needed for 10x Chromium scRNA-seq analysis from the following file content.
This file contains sample information from a BCL run executed on a Local NovaSeq.

File name: {file_name}

File content:
```
{file_content}
```

This file typically contains the following information:
- Sample names (e.g., PB, LSK) and their corresponding Dual Index (TT/NT)
- Hashtag numbers (1-24) and corresponding barcode sequences (15bp)
- Sample-to-Hashtag mapping
- GEX libraries use Dual Index TT Set A; CMO/Hashtag libraries use Dual Index NT Set A

**IMPORTANT: Rows sharing the same Dual Index TT (e.g., SI-TT-H12) belong to the same GEM well.
Group samples by GEM well.**

Return in the following JSON format:

```json
{{
  "mode": "multi or count",
  "gex_sample_sheet": [
    {{
      "lane": "*",
      "sample": "GEX sample name (e.g., PB, LSK, WT)",
      "index": "SI-TT-XXX"
    }}
  ],
  "cmo_sample_sheet": [
    {{
      "lane": "*",
      "sample": "CMO sample name (e.g., PB, LSK, WT)",
      "index": "SI-NT-XXX"
    }}
  ],
  "hashtag_numbers": [1, 2, 3, 4],
  "species": "mouse or human",
  "totalseq_series": "B or C (null if unknown)",
  "hashtags_with_sequence": [
    {{
      "hashtag_number": 1,
      "sequence": "ACCCACCAGTAAGAC",
      "sample_label": "corresponding label"
    }}
  ],
  "samples": [
    {{
      "sample_id": "WT_1",
      "hashtag_numbers": [1],
      "group": "WT_PB",
      "gex_sample": "PB"
    }},
    {{
      "sample_id": "WT_2",
      "hashtag_numbers": [2],
      "group": "WT_PB",
      "gex_sample": "PB"
    }},
    {{
      "sample_id": "HFpEF_1",
      "hashtag_numbers": [5],
      "group": "HFpEF_PB",
      "gex_sample": "PB"
    }}
  ]
}}
```

Rules:
- **mode**: "multi" if CMO/Hashtag library (NT index) exists, otherwise "count"
- **gex_sample_sheet / cmo_sample_sheet**: **Only one row per physical GEM well**. Never duplicate rows with the same index.
  One GEM well may contain samples from multiple biological conditions, but the mkfastq sample sheet is per GEM well.
  Example: Even if the PB well (SI-TT-H12) mixes WT and HFpEF mice, gex_sample_sheet has only one "PB" row.
  Sample names should come from the Excel sample label column (e.g., "PB", "LSK", "WT", "KO_Iso", "AAV_MSC").
  **Do not append index numbers to sample names** (e.g., "PB" not "PB_H12").
  Only when multiple conditions share the same Dual Index in one well, use the well name instead of condition names.
  If each GEM well has a different index, use the respective sample labels as sample names.
  "TT" and "NT" are index set names, not sample names.
- **hashtag_numbers**: Hashtag numbers (1-24). Not product code numbers.
  Example: TotalSeq-B0301 = Mouse Hashtag **1**, B0302 = **2**, ..., B0308 = **8**
  Example: TotalSeq-B0251 = Human Hashtag **1**, ..., B0258 = **8**
- **species**: mouse or human. Infer from sample names and context.
- **hashtags_with_sequence**: Only if barcode sequences are listed in the file.
- **samples**: **1 hashtag = 1 sample is the principle.** Describe each individual (each hashtag) as a separate sample.
  Do not combine multiple hashtags into 1 sample.
  - sample_id: name identifying the individual (e.g., WT_1, WT_2, HFpEF_1)
  - hashtag_numbers: hashtag number(s) for this sample (usually 1)
  - group: condition group label + GEM well name (e.g., WT_PB, HFpEF_LSK). Used as parent directory for output.
  - gex_sample: **Must exactly match a sample name in gex_sample_sheet** (e.g., "PB", "LSK"). GEM well name, not condition name.
  - Samples sharing the same TT index have the same gex_sample.
- If no CMO / Hashtag library: cmo_sample_sheet is an empty list [], hashtag-related fields are also empty.
- Return JSON only. No explanatory text.
"""
    elif extraction_type == "count":
        # count mode: extract sample names, species, and library_classification
        _fq_info = ""
        if fastq_filenames:
            _fq_prefixes = set()
            for _fn in fastq_filenames:
                _m = re.match(r'^(.+?)_S\d+_L\d{3}', _fn)
                if _m:
                    _fq_prefixes.add(_m.group(1))
            _prefix_types = classify_all_prefixes(_fq_prefixes)
            _fq_info = f"""

FASTQ filenames:
```
{chr(10).join(fastq_filenames)}
```

Known abbreviation patterns:
- Gene Expression: GEX, GE, RNA, GeneExpression, Transcriptome
- Antibody Capture (Surface Protein/ADT): SP, ADT, AB, SurfaceProtein, CITE
- Multiplexing Capture (CMO/Hashtag): CMO, HTO, Hashtag, Multiplex

Pre-classification by pattern matching:
{json.dumps(_prefix_types, ensure_ascii=False, indent=2)}
"""

        prompt = f"""Extract the information needed for scRNA-seq Cell Ranger count analysis from the following file content.

File name: {file_name}

File content:
```
{file_content}
```
{_fq_info}
Return in the following JSON format. Use null for missing information:

```json
{{
  "sample_names": ["sample_name_1", "sample_name_2"],
  "species": "mouse or human",
  "library_classification": {{
    "FASTQ_prefix_1": "Gene Expression / Antibody Capture / Multiplexing Capture"
  }}
}}
```

Notes:
- sample_names: FASTQ prefix names to specify with cellranger count --sample
- species: mouse or human. null if undetermined.
- library_classification: library type for each FASTQ prefix (see patterns above)
- Return JSON only. No explanatory text.
"""
    else:
        # Add FASTQ filename information
        _fq_info = ""
        if fastq_filenames:
            # Extract prefixes from FASTQ filenames and classify by pattern matching
            _fq_prefixes = set()
            for _fn in fastq_filenames:
                _m = re.match(r'^(.+?)_S\d+_L\d{3}', _fn)
                if _m:
                    _fq_prefixes.add(_m.group(1))
            _prefix_types = classify_all_prefixes(_fq_prefixes)
            _fq_gex_pfx = sorted([p for p, t in _prefix_types.items() if t == "Gene Expression"])
            _fq_adt_pfx = sorted([p for p, t in _prefix_types.items() if t == "Antibody Capture"])
            _fq_cmo_pfx = sorted([p for p, t in _prefix_types.items() if t == "Multiplexing Capture"])

            # If SP/ADT is used for hashtags, reclassify as CMO (when no dedicated CMO prefix exists)
            if _fq_adt_pfx and not _fq_cmo_pfx:
                _fq_cmo_pfx = _fq_adt_pfx
                _fq_adt_pfx = []

            _fq_info = f"""

FASTQ filenames (actual files in the data directory):
```
{chr(10).join(fastq_filenames)}
```

Known abbreviation patterns:
- Gene Expression: GEX, GE, RNA, GeneExpression, Transcriptome
- Antibody Capture (Surface Protein/ADT): SP, ADT, AB, SurfaceProtein, CITE
- Multiplexing Capture (CMO/Hashtag): CMO, HTO, Hashtag, Multiplex

Classify each FASTQ prefix based on the sample sheet context and the patterns above.

GEX FASTQ prefixes (= GEM well names): {', '.join(_fq_gex_pfx)}
ADT FASTQ prefixes (Antibody Capture): {', '.join(_fq_adt_pfx) if _fq_adt_pfx else '(none)'}
CMO FASTQ prefixes (Multiplexing Capture): {', '.join(_fq_cmo_pfx)}

Determine the fastq_id for each GEX / ADT / CMO library from the filenames.
Example: "Control_S0_L001_I1_001.fastq.gz" -> GEX fastq_id = "Control"
    "Control-Hashtag_S0_L001_I1_001.fastq.gz" -> CMO fastq_id = "Control-Hashtag"
    "MSC_AAVshRNA_SP_S0_L001_R1_001.fastq.gz" -> ADT fastq_id = "MSC_AAVshRNA_SP"
If there are multiple GEX/CMO pairs, multiple GEM wells exist.

**IMPORTANT**: Always use the GEX FASTQ prefix names for the samples' gex_sample field.
Infer the mapping between Excel genotype/condition names and FASTQ prefixes from context.
Example: genotype "CKO" -> GEX prefix "CCL3-CKO", genotype "WT" -> GEX prefix "Control"
"""

        prompt = f"""Extract the information needed for scRNA-seq Cell Ranger multi analysis from the following file content.

File name: {file_name}

File content:
```
{file_content}
```
{_fq_info}
This file typically contains the following information:
- Hashtag numbers used (e.g., Hashtag 5, 6, 7, 8)
- TotalSeq series (A, B, C) and species (human, mouse)
- Sample names and which Hashtag corresponds to each
- Dual Index information (for GEX/CMO)
- Barcode sequences (if included)

Return in the following JSON format. Use null for missing information:

```json
{{
  "hashtag_numbers": [5, 6, 7, 8],
  "species": "mouse or human",
  "totalseq_series": "B or C (default B if unknown)",
  "hashtags_with_sequence": [
    {{
      "hashtag_number": 5,
      "sequence": "barcode sequence 15bp (only if listed in the file)",
      "sample_label": "sample name/label corresponding to this Hashtag (if available)"
    }}
  ],
  "samples": [
    {{
      "sample_id": "individual name (e.g., WT_1, WT_2, HFpEF_1)",
      "hashtag_numbers": [5],
      "group": "condition group label (e.g., WT, HFpEF, KO)",
      "gex_sample": "corresponding gex_sample_sheet sample name (which GEM well)",
      "cmo_ids": null
    }}
  ],
  "libraries": {{
    "gex": {{
      "fastq_id": "GEX FASTQ sample ID (must match gex_sample_sheet sample name)",
      "description": "Gene Expression"
    }},
    "adt": {{
      "fastq_id": "ADT FASTQ sample ID (only if Antibody Capture exists, otherwise null)",
      "description": "Antibody Capture"
    }},
    "cmo": {{
      "fastq_id": "CMO FASTQ sample ID (must match cmo_sample_sheet sample name)",
      "description": "Multiplexing Capture"
    }}
  }},
  "library_classification": {{
    "FASTQ_prefix_1": "Gene Expression / Antibody Capture / Multiplexing Capture",
    "FASTQ_prefix_2": "classification result"
  }},
  "gex_sample_sheet": [
    {{
      "lane": "*",
      "sample": "GEX sample name",
      "index": "SI-TT-XXX"
    }}
  ],
  "cmo_sample_sheet": [
    {{
      "lane": "*",
      "sample": "CMO sample name",
      "index": "SI-NT-XXX"
    }}
  ]
}}
```

Notes:
- hashtag_numbers: **Hashtag numbers** (1-24) as a list. **Not product code numbers**.
  Example: TotalSeq-B0255 = Hashtag **5**, TotalSeq-B0301 = Ms.Hashtag **1**
  Product code to Hashtag number mapping:
    - Human B/C: B0251=1, B0252=2, ..., B0255=5, B0256=6, B0257=7, B0258=8, ..., B0260=10
    - Mouse B/C: B0301=1, B0302=2, ..., B0305=5, ..., B0310=10, ..., B0319=19, B0324=24
- species: mouse or human. null if undetermined.
- hashtags_with_sequence: Only if barcode sequences are listed in the file. Otherwise empty list.
- samples: **1 hashtag = 1 sample is the principle.** Describe each individual (each hashtag) as a separate sample.
  Do not combine multiple hashtags into 1 sample.
  - sample_id: name identifying the individual (e.g., WT_1, WT_2, HFpEF_1)
  - group: condition group label (e.g., WT, HFpEF, KO). Used as parent directory for output.
  - Samples sharing the same Dual Index (SI-TT-XX) belong to the same GEM well group.
  - hashtag_numbers in samples should also use Hashtag numbers (1-24).
  - **gex_sample**: Must specify the gex_sample_sheet sample name to indicate which GEM well the sample belongs to.
    Example: If gex_sample_sheet has sample="LSK_HFD0w" and this sample is from that GEM well, gex_sample="LSK_HFD0w"
  - When there are multiple GEM wells, the mapping of each sample to its well is critical.
- Indexes must use the full form with "SI-TT-" or "SI-NT-" prefix.
- **libraries fastq_id**: Must match the sample name in gex_sample_sheet / cmo_sample_sheet.
  "TT" and "NT" are index set names, not fastq_ids.
  Example: If gex_sample_sheet has sample="WT", then libraries.gex.fastq_id should also be "WT".
  For multiple GEM wells, use the first well's name for the libraries fastq_id.
- **Built-in CMO ID**: Some experiments use Cell Ranger built-in CMO barcodes
  (CMO301-CMO312) instead of TotalSeq hashtags. In this case:
  - Do not include these numbers in hashtag_numbers (hashtag_numbers is for TotalSeq hashtag numbers only)
  - Specify IDs directly in the sample's cmo_ids: "cmo_ids": ["CMO305", "CMO306", ...]
  - Example: 4 GEM wells, 2 using TotalSeq Hashtag 5-8 and 2 using CMO305-CMO311
    -> hashtag_numbers: [5, 6, 7, 8]; TotalSeq well samples use hashtag_numbers;
       CMO well samples use cmo_ids: ["CMO305"] directly.
- cmo_ids: null for normal samples using TotalSeq hashtags (auto-resolved from hashtag_numbers).
  Only specify directly for samples using Built-in CMO IDs.
- Return JSON only. No explanatory text.
"""

    try:
        client = anthropic.Anthropic(api_key=api_key)
        max_tok = 4000 if extraction_type in ("local", "multi") else 2000
        response = client.messages.create(
            model="claude-sonnet-4-6",
            max_tokens=max_tok,
            messages=[{"role": "user", "content": prompt}],
        )
        text = response.content[0].text

        # Extract JSON portion
        json_match = re.search(r'```json\s*(.*?)\s*```', text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group(1))
        # Also try JSON without ``` delimiters
        json_match = re.search(r'\{.*\}', text, re.DOTALL)
        if json_match:
            return json.loads(json_match.group(0))
        st.error("Could not extract JSON from the LLM response.")
        return None
    except Exception as e:
        st.error(f"LLM API error: {e}")
        return None


def resolve_hashtags_from_llm(llm_result):
    """Resolve hashtag information from LLM extraction results using the master DB

    Returns: list of dicts with id, name, read, pattern, sequence, feature_type
    """
    db = load_hashtag_db()
    if db.empty:
        return _fallback_hashtags_from_llm(llm_result)

    numbers = llm_result.get("hashtag_numbers", [])
    species = llm_result.get("species", "mouse")
    series = llm_result.get("totalseq_series", "B")

    if not numbers:
        return _fallback_hashtags_from_llm(llm_result)

    # Strategy 1: Filter by hashtag_number from master DB
    filtered = db[
        (db["species"] == species) &
        (db["series"] == series) &
        (db["hashtag_number"].isin(numbers))
    ]

    if filtered.empty:
        # Filter by species only
        filtered = db[
            (db["species"] == species) &
            (db["hashtag_number"].isin(numbers))
        ]
        if not filtered.empty:
            filtered = filtered.drop_duplicates(subset=["hashtag_number"], keep="first")

    # Strategy 2: Handle case where LLM returned product code numbers (255 -> B0255)
    if filtered.empty and numbers:
        # Build IDs from product code numbers and match
        product_ids = [f"HTO_{series}{str(n).zfill(4)}" for n in numbers]
        filtered = db[db["id"].isin(product_ids)]
        if filtered.empty:
            # Try with species only
            for s in ["B", "C"]:
                product_ids = [f"HTO_{s}{str(n).zfill(4)}" for n in numbers]
                filtered = db[db["id"].isin(product_ids)]
                if not filtered.empty:
                    break

    # Strategy 3: Match by sequence if LLM directly extracted barcode sequences
    if filtered.empty:
        seqs = [h.get("sequence") for h in llm_result.get("hashtags_with_sequence", []) if h.get("sequence")]
        if seqs:
            filtered = db[db["sequence"].isin(seqs)]

    if filtered.empty:
        st.warning(f"No matching hashtags found in master DB (species={species}, series={series}, numbers={numbers})")
        return _fallback_hashtags_from_llm(llm_result)

    hashtags = []
    for _, row in filtered.iterrows():
        hashtags.append({
            "id": row["id"],
            "name": row["name"],
            "read": row["read"],
            "pattern": row["pattern"],
            "sequence": row["sequence"],
            "feature_type": row["feature_type"],
        })

    return hashtags


def _fallback_hashtags_from_llm(llm_result):
    """Use sequence information directly extracted by LLM when master DB is unavailable"""
    hashtags = []
    for h in llm_result.get("hashtags_with_sequence", []):
        if h.get("sequence"):
            num = h.get("hashtag_number", "?")
            hashtags.append({
                "id": f"Hashtag_{num}",
                "name": h.get("sample_label", f"Hashtag_{num}"),
                "read": "R2",
                "pattern": "5PNNNNNNNNNN(BC)",
                "sequence": h["sequence"],
                "feature_type": "Multiplexing Capture",
            })
    return hashtags


def resolve_samples_from_llm(llm_result, hashtags):
    """Map LLM sample information to hashtag IDs"""
    # Build hashtag_number -> ID mapping table
    # Supports both master DB hashtag_number (5,6,7,8) and product code numbers (255,256,...)
    db = load_hashtag_db()
    num_to_id = {}
    for h in hashtags:
        hid = h["id"]
        # Get hashtag_number from master DB
        db_row = db[db["id"] == hid]
        if not db_row.empty:
            ht_num = int(db_row.iloc[0]["hashtag_number"])
            num_to_id[ht_num] = hid
            # Also map product code numbers (e.g. 255 from HTO_B0255)
            code_match = re.search(r'(\d+)$', hid)
            if code_match:
                num_to_id[int(code_match.group(1))] = hid
        # Also match by number contained in the ID
        code_match = re.search(r'(\d+)$', hid)
        if code_match:
            num_to_id[int(code_match.group(1))] = hid

    samples = []
    for s in llm_result.get("samples", []):
        # Use directly if Built-in CMO ID is specified
        direct_cmo = s.get("cmo_ids")
        if direct_cmo and isinstance(direct_cmo, list) and all(
            isinstance(c, str) and re.match(r'^CMO\d+$', c) for c in direct_cmo
        ):
            cmo_ids = direct_cmo
        elif direct_cmo and isinstance(direct_cmo, str) and re.match(r'^CMO\d+$', direct_cmo):
            cmo_ids = [direct_cmo]
        else:
            # Normal hashtag_numbers -> hashtag ID resolution
            ht_numbers = s.get("hashtag_numbers", [])
            cmo_ids = []
            for num in ht_numbers:
                if num in num_to_id:
                    cmo_ids.append(num_to_id[num])
                else:
                    # Fallback: search the hashtag list directly
                    found = False
                    for h in hashtags:
                        if str(num) in h["id"]:
                            cmo_ids.append(h["id"])
                            found = True
                            break
                    if not found:
                        cmo_ids.append(f"Hashtag_{num}")
        # Sanitize sample_id for Cell Ranger compatibility (ASCII only, spaces -> _)
        sid = s["sample_id"]
        sid = re.sub(r'[^\x20-\x7E]', '', sid)  # Remove non-ASCII characters
        sid = sid.strip().replace(" ", "_").replace("(", "").replace(")", "")
        sid = re.sub(r'_+', '_', sid).strip("_")  # Collapse consecutive underscores
        if not sid:
            sid = f"Sample_{len(samples)+1}"
        entry = {
            "sample_id": sid,
            "cmo_ids": cmo_ids,
        }
        # Retain GEM well association info
        if s.get("gex_sample"):
            entry["gex_sample"] = s["gex_sample"]
        # Retain group info (for output folder structure)
        if s.get("group"):
            grp_name = s["group"]
            grp_name = re.sub(r'[^\x20-\x7E]', '', grp_name)
            grp_name = grp_name.strip().replace(" ", "_").replace("(", "").replace(")", "")
            grp_name = re.sub(r'_+', '_', grp_name).strip("_")
            entry["group"] = grp_name
        samples.append(entry)
    return samples


def read_uploaded_file_content(uploaded_file):
    """Read the contents of an uploaded file as a string"""
    name = uploaded_file.name.lower()
    if name.endswith((".xlsx", ".xls")):
        try:
            xls = pd.ExcelFile(uploaded_file)
            parts = []
            for sheet in xls.sheet_names:
                df = pd.read_excel(xls, sheet_name=sheet, header=None)
                parts.append(f"=== Sheet: {sheet} ===\n{df.to_string(index=False)}")
            uploaded_file.seek(0)
            return "\n\n".join(parts)
        except Exception as e:
            st.error(f"Excel read error: {e}")
            return None
    else:
        # Text-based files (csv, tsv, txt, etc.)
        try:
            content = uploaded_file.read().decode("utf-8")
            uploaded_file.seek(0)
            return content
        except UnicodeDecodeError:
            try:
                uploaded_file.seek(0)
                content = uploaded_file.read().decode("shift_jis")
                uploaded_file.seek(0)
                return content
            except Exception:
                st.error("Could not detect the file encoding.")
                return None


# =============================================================================
#  FASTQ rename utility (-> Cell Ranger format)
# =============================================================================
#
# FASTQ filename format required by Cell Ranger:
#   {SampleName}_S{Number}_L00{Lane}_R{Read}_001.fastq.gz
#
#   - SampleName: arbitrary (must match the fastq_id in cellranger multi)
#   - S{Number}: sample number (e.g. S1; Cell Ranger ignores the value)
#   - L00{Lane}: lane number (optional in v4.0+; safe to standardize as L001)
#   - R{Read}: R1 = Read 1, R2 = Read 2 (required; I1/I2 are index reads)
#   - 001: file chunk number (always 001)
#   - Extension: .fastq.gz
#
# Files already in Cell Ranger format match _S\d+_L\d{3}_R[12]_001\.fastq\.gz.

CELLRANGER_FASTQ_RE = re.compile(
    r'^.+_S\d+_L\d{3}_[RI][12]_001\.fastq\.gz$'
)

RENAME_PATTERNS = [
    {
        "name": "Azenta: {SAMPLE}_L{N}_{READ}.fq.gz",
        "pattern": r'^(.+?)_L\d+_(\d+)\.fq\.gz$',
        "replacement": lambda m: f"{m.group(1)}_S1_L001_R{m.group(2)}_001.fastq.gz",
    },
    {
        "name": "Azenta: {SAMPLE}_{READ}.fq.gz",
        "pattern": r'^(.+?)_(\d+)\.fq\.gz$',
        "replacement": lambda m: f"{m.group(1)}_S1_L001_R{m.group(2)}_001.fastq.gz",
    },
    {
        "name": "{SAMPLE}_L{N}_{READ}.fastq.gz (lane → L001)",
        "pattern": r'^(.+?)_L\d+_(\d+)\.fastq\.gz$',
        "replacement": lambda m: f"{m.group(1)}_S1_L001_R{m.group(2)}_001.fastq.gz",
    },
]


def is_cellranger_format(filenames):
    """Check if all files are in Cell Ranger format"""
    return all(CELLRANGER_FASTQ_RE.match(fn) for fn in filenames)


def detect_rename_mapping(filenames):
    """Detect applicable rename patterns from a list of filenames and return a mapping"""
    # No renaming needed if already in Cell Ranger format
    if is_cellranger_format(filenames):
        return None, "already Cell Ranger format"

    for pat_info in RENAME_PATTERNS:
        mapping = {}
        matched_all = True
        for fn in filenames:
            m = re.match(pat_info["pattern"], fn)
            if m:
                mapping[fn] = pat_info["replacement"](m)
            else:
                matched_all = False
                break
        if matched_all and mapping:
            return mapping, pat_info["name"]
    return None, None


def apply_custom_rename(filenames, pattern, replacement):
    """Generate a rename mapping using a user-specified regex"""
    mapping = {}
    try:
        regex = re.compile(pattern)
    except re.error as e:
        st.error(f"Invalid regex: {e}")
        return None
    for fn in filenames:
        m = regex.match(fn)
        if m:
            try:
                new_name = m.expand(replacement)
                mapping[fn] = new_name
            except Exception as e:
                st.error(f"Replacement error for {fn}: {e}")
                return None
        else:
            st.warning(f"Pattern did not match: {fn}")
            return None
    return mapping


def rename_with_llm(filenames):
    """Generate FASTQ file rename mapping using LLM

    Fallback for filenames that don't match regex patterns.
    Cell Ranger format: {SampleName}_S1_L001_R{1or2}_001.fastq.gz
    """
    # No LLM needed if already in Cell Ranger format
    if all(CELLRANGER_FASTQ_RE.match(f) for f in filenames):
        st.info("Files are already in Cell Ranger format. No renaming needed.")
        mapping = {f: f for f in filenames}
        return mapping, "Already in Cell Ranger format"

    try:
        import anthropic
    except ImportError:
        st.error("The anthropic package is required. Run `pip install anthropic`.")
        return None

    api_key = os.getenv("ANTHROPIC_API_KEY")
    if not api_key:
        st.error("ANTHROPIC_API_KEY is not set.")
        return None

    file_list = "\n".join(filenames)

    prompt = f"""Rename the following FASTQ filenames to the format required by Cell Ranger.

## Cell Ranger FASTQ filename format (required)

```
{{SampleName}}_S1_L001_R{{Read}}_001.fastq.gz
```

Field descriptions:
- **SampleName**: Sample name. Preserve the sample-identifying portion from the original filename.
  This name will be used as the fastq_id in cellranger multi, so it should be meaningful.
- **S1**: Sample number. Always S1 is fine.
- **L001**: Lane number. Always standardize to L001.
- **R1 / R2**: Read 1 or Read 2. Determine from the original filename.
  - "1", "R1", "_1." -> R1
  - "2", "R2", "_2." -> R2
- **001**: Chunk number. Always 001.
- Extension must be **.fastq.gz** (convert .fq.gz -> .fastq.gz)

## Files to rename

```
{file_list}
```

## Output format

Return in the following JSON format. Return **JSON only**:

```json
{{
  "mapping": {{
    "original_filename_1": "new_filename_1",
    "original_filename_2": "new_filename_2"
  }},
  "detected_pattern": "one-line description of the detected filename pattern"
}}
```

Notes:
- The SampleName portion of paired R1 and R2 must match exactly.
- R1/R2 of the same sample must use the same SampleName.
- Return mappings for all files.
"""

    try:
        client = anthropic.Anthropic(api_key=api_key)
        response = client.messages.create(
            model="claude-sonnet-4-6",
            max_tokens=2000,
            messages=[{"role": "user", "content": prompt}],
        )
        text = response.content[0].text
        if not text:
            st.error("Empty response from LLM. Please try again.")
            return None

        json_match = re.search(r'```json\s*(.*?)\s*```', text, re.DOTALL)
        if json_match:
            result = json.loads(json_match.group(1))
        else:
            json_match = re.search(r'\{.*\}', text, re.DOTALL)
            if json_match:
                result = json.loads(json_match.group(0))
            else:
                st.error("Could not extract JSON from the LLM response.")
                return None

        mapping = result.get("mapping", {})
        pattern_desc = result.get("detected_pattern", "")

        # Validation: check if all renamed files are in Cell Ranger format
        invalid = [v for v in mapping.values() if not CELLRANGER_FASTQ_RE.match(v)]
        if invalid:
            st.warning(f"Some filenames generated by LLM do not match Cell Ranger format: {invalid[:3]}")

        return mapping, pattern_desc

    except Exception as e:
        st.error(f"LLM rename error: {e}")
        return None


# =============================================================================
#  multi.csv generation
# =============================================================================

def build_multi_csv(config):
    """Generate multi.csv content from multi mode configuration"""
    lines = []

    # [gene-expression]
    lines.append("[gene-expression]")
    lines.append(f"reference,{config['reference']}")
    if config.get("hashtags"):
        lines.append(f"cmo-set,{config['feature_ref_path']}")
    if config.get("create_bam"):
        lines.append("create-bam,true")
    lines.append("")

    # [libraries]
    lines.append("[libraries]")
    lines.append("fastq_id,fastqs,feature_types")
    for lib in config.get("libraries", []):
        lines.append(f"{lib['fastq_id']},{lib['fastqs']},{lib['feature_types']}")
    lines.append("")

    # [samples]
    if config.get("samples"):
        lines.append("[samples]")
        lines.append("sample_id,cmo_ids")
        for s in config["samples"]:
            cmo_str = s["cmo_ids"] if isinstance(s["cmo_ids"], str) else "|".join(s["cmo_ids"])
            lines.append(f"{s['sample_id']},{cmo_str}")

    return "\n".join(lines)


def build_feature_ref_csv(hashtags):
    """Generate feature reference CSV from a hashtag list"""
    lines = ["id,name,read,pattern,sequence,feature_type"]
    for h in hashtags:
        lines.append(f"{h['id']},{h['name']},{h['read']},{h['pattern']},{h['sequence']},{h['feature_type']}")
    return "\n".join(lines)


# =============================================================================
#  Job management -- same pattern as fastq2rna.py
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
#  Common settings UI (shared by count / multi)
# =============================================================================

def common_cellranger_settings(key_prefix):
    """Return common Cell Ranger settings"""
    col1, col2 = st.columns(2)
    with col1:
        reference_name = st.selectbox(
            "Reference genome",
            list(REFERENCES.keys()),
            key=f"{key_prefix}_reference",
        )
    with col2:
        create_bam = st.checkbox("Create BAM (--create-bam)", value=True, key=f"{key_prefix}_create_bam")
        also_run_count = st.checkbox(
            "Also run cellranger count for GEX (combined BAM)",
            value=False,
            key=f"{key_prefix}_also_run_count",
            help="Additionally run cellranger count on GEX FASTQ in multi mode to generate "
                 "a combined BAM (possorted_genome_bam.bam) containing all cells before demultiplexing."
        )
        run_annotate = st.checkbox(
            "Run cell type annotation (cellranger annotate)",
            value=False,
            key=f"{key_prefix}_run_annotate",
            help="Run automatic cell type annotation after cellranger count/multi completes (beta). "
                 "Requires a 10x Genomics Cloud Access Token.",
        )

    reference = REFERENCES[reference_name]

    if also_run_count:
        st.info(
            "**BAM files in multiplexing mode:**\n"
            "- Cell Ranger multi -> per-sample BAM (`sample_alignments.bam`): individual BAM per sample after demultiplexing\n"
            "- Cell Ranger count -> combined BAM (`possorted_genome_bam.bam`): combined BAM of all cells\n\n"
            "**per-sample BAM**: For use when Cell Ranger demultiplexing results are used directly\n"
            "**combined BAM**: Needed when re-demultiplexing with SCALA etc., or running DropletQC / velocyto "
            "on all cells"
        )

    col3, col4 = st.columns(2)
    with col3:
        localcores = st.number_input(
            "Local cores (--localcores)",
            min_value=1, max_value=256,
            value=min(64, os.cpu_count() or 32),
            key=f"{key_prefix}_localcores",
        )
    with col4:
        try:
            with open('/proc/meminfo') as _mf:
                _mi = {l.split(':')[0]: int(l.split()[1]) for l in _mf if len(l.split()) >= 2}
            _avail_gb = max(1, _mi.get('MemAvailable', 0) // (1024 * 1024))
        except Exception:
            _avail_gb = max(1, os.sysconf('SC_PAGE_SIZE') * os.sysconf('SC_PHYS_PAGES') // (1024 ** 3))
        localmem = st.number_input(
            f"Local memory GB (--localmem) [available: {_avail_gb} GB]",
            min_value=32, max_value=512,
            value=max(32, min(128, _avail_gb)),
            key=f"{key_prefix}_localmem",
        )

    return reference, reference_name, create_bam, also_run_count, int(localcores), int(localmem), run_annotate


# =============================================================================
#  count mode settings UI
# =============================================================================

def count_mode_ui(key_prefix, fastq_dir_default="", sample_default="",
                  fastq_filenames=None):
    """Settings UI for count mode"""
    st.markdown("#### Cell Ranger count settings")

    # --- LLM extraction from sample sheet (optional) ---
    llm_key = f"{key_prefix}_count_llm_result"
    st.markdown("**Sample information (optional: Excel/Text upload + LLM extraction)**")
    st.caption("Upload a sample sheet to automatically infer sample names and library types.")

    count_uploaded = st.file_uploader(
        "Sample info file",
        type=["xlsx", "xls", "csv", "tsv", "txt"],
        key=f"{key_prefix}_count_sample_info_file",
    )

    if count_uploaded is not None:
        if st.button("Extract with LLM", key=f"{key_prefix}_count_extract_llm"):
            try:
                file_name = count_uploaded.name
                if file_name.endswith((".xlsx", ".xls")):
                    import openpyxl
                    xls = pd.ExcelFile(count_uploaded)
                    sheets_text = []
                    for sheet_name in xls.sheet_names:
                        df = xls.parse(sheet_name)
                        sheets_text.append(f"--- Sheet: {sheet_name} ---\n{df.to_string()}")
                    file_content = "\n\n".join(sheets_text)
                else:
                    file_content = count_uploaded.read().decode("utf-8", errors="replace")

                with st.spinner("Extracting sample info with LLM..."):
                    result = extract_sample_info_with_llm(
                        file_content, file_name,
                        extraction_type="count",
                        fastq_filenames=fastq_filenames,
                    )
                if result:
                    st.session_state[llm_key] = result
                    st.success("LLM extraction complete")
                    st.rerun()
            except Exception as e:
                st.error(f"File read error: {e}")

    # Apply LLM results
    llm_count_data = st.session_state.get(llm_key, {})
    if llm_count_data:
        with st.expander("LLM extraction result", expanded=False):
            st.json(llm_count_data)
        # Auto-set sample names
        llm_sample_names = llm_count_data.get("sample_names", [])
        if llm_sample_names and not sample_default:
            sample_default = ", ".join(llm_sample_names)
        # Display library classification
        lib_cls = llm_count_data.get("library_classification", {})
        if lib_cls:
            st.info("Library classification: " + ", ".join(f"{k} = {v}" for k, v in lib_cls.items()))

    sample_name = st.text_input(
        "Sample name (FASTQ sample prefix)",
        value=sample_default,
        key=f"{key_prefix}_count_sample_name",
        placeholder="e.g. SampleX",
    )
    fastqs_dir = st.text_input(
        "FASTQ directory",
        value=fastq_dir_default,
        key=f"{key_prefix}_count_fastqs_dir",
    )

    return sample_name, fastqs_dir


# =============================================================================
#  multi mode settings UI
# =============================================================================

def multi_mode_ui(key_prefix, fastq_dir_default="", multi_gem_well=False,
                   fastq_filenames=None):
    """Settings UI for multi mode -- hashtag, libraries, samples.
    When multi_gem_well=True, Libraries/Samples are hidden because they are edited in the GEM well groups section.
    """
    st.markdown("#### Cell Ranger multi settings")

    # --- Sample information extraction via LLM ---
    st.markdown("**Sample information (Excel/Text upload + LLM extraction)**")
    st.caption("Upload a sample sheet to infer hashtag numbers and automatically look up barcode sequences from the master DB.")

    info_file = st.file_uploader(
        "Sample info file (Excel, CSV, TSV, TXT)",
        type=["xlsx", "xls", "csv", "tsv", "txt"],
        key=f"{key_prefix}_info_file",
    )

    llm_key = f"{key_prefix}_llm_result"
    resolved_ht_key = f"{key_prefix}_resolved_hashtags"
    resolved_samples_key = f"{key_prefix}_resolved_samples"

    if info_file is not None:
        if st.button("Extract with LLM", type="primary", key=f"{key_prefix}_extract_llm"):
            content = read_uploaded_file_content(info_file)
            if content:
                with st.spinner("Extracting sample info with LLM..."):
                    result = extract_sample_info_with_llm(
                        content, info_file.name,
                        fastq_filenames=fastq_filenames,
                    )
                    if result:
                        st.session_state[llm_key] = result
                        # Look up hashtags from master DB
                        resolved = resolve_hashtags_from_llm(result)
                        st.session_state[resolved_ht_key] = resolved
                        # Also look up sample information
                        resolved_samples = resolve_samples_from_llm(result, resolved)
                        st.session_state[resolved_samples_key] = resolved_samples

                        if resolved:
                            # Directly update the Hashtag CSV widget
                            ht_csv = "id,name,read,pattern,sequence,feature_type\n"
                            for h in resolved:
                                ht_csv += f"{h['id']},{h['name']},{h['read']},{h['pattern']},{h['sequence']},{h['feature_type']}\n"
                            st.session_state[f"{key_prefix}_hashtag_csv"] = ht_csv
                            species = result.get("species", "?")
                            series = result.get("totalseq_series", "?")
                            nums = result.get("hashtag_numbers", [])
                            st.success(f"Extraction complete! TotalSeq-{series} ({species}), Hashtag {nums} -> {len(resolved)} entries resolved from master DB")
                        else:
                            st.warning("Could not identify hashtag numbers. Please enter them manually.")

    # --- Hashtag (Feature Reference) ---
    st.markdown("**Hashtag / Feature Reference**")
    st.caption("Only the hashtags being used will be included in the feature reference.")

    # Set default values from LLM + master DB lookup results
    resolved_hashtags = st.session_state.get(resolved_ht_key, [])

    if resolved_hashtags:
        default_ht_text = "id,name,read,pattern,sequence,feature_type\n"
        for h in resolved_hashtags:
            default_ht_text += f"{h['id']},{h['name']},{h['read']},{h['pattern']},{h['sequence']},{h['feature_type']}\n"
    else:
        default_ht_text = "id,name,read,pattern,sequence,feature_type\n"

    hashtag_text = st.text_area(
        "Hashtag CSV (feature reference) -- editable",
        value=default_ht_text,
        height=150,
        key=f"{key_prefix}_hashtag_csv",
    )

    # Parse
    hashtags = []
    if hashtag_text.strip():
        lines = [l.strip() for l in hashtag_text.strip().split("\n") if l.strip()]
        if len(lines) > 1:
            header = lines[0].lower()
            if "id" in header and "sequence" in header:
                for line in lines[1:]:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 6:
                        hashtags.append({
                            "id": parts[0], "name": parts[1], "read": parts[2],
                            "pattern": parts[3], "sequence": parts[4], "feature_type": parts[5],
                        })
    if hashtags:
        st.success(f"{len(hashtags)} hashtags defined (only those in use)")

    # --- Libraries / Samples ---
    # Multiple GEM wells -> editable in the GEM well groups section, so skip here
    if multi_gem_well:
        st.caption("Libraries / Samples can be edited in the **GEM well groups** section below.")
        return hashtags, [], []

    # --- Auto-detect sample prefixes from FASTQ filenames ---
    llm_data = st.session_state.get(llm_key, {})
    llm_libs = llm_data.get("libraries", {})

    _gex_prefixes = []
    _adt_prefixes = []
    _cmo_prefixes = []
    if fastq_filenames:
        _prefixes = set()
        for fn in fastq_filenames:
            m = re.match(r'^(.+?)_S\d+_L\d{3}', fn)
            if m:
                _prefixes.add(m.group(1))
        _llm_lib_cls = llm_data.get("library_classification", {})
        _prefix_types = classify_all_prefixes(_prefixes, _llm_lib_cls)
        _gex_prefixes = sorted([p for p, t in _prefix_types.items() if t == "Gene Expression"])
        _adt_prefixes = sorted([p for p, t in _prefix_types.items() if t == "Antibody Capture"])
        _cmo_prefixes = sorted([p for p, t in _prefix_types.items() if t == "Multiplexing Capture"])

        # Context-dependent SP/ADT classification: when no CMO prefix exists but ADT does,
        # promote ADT -> Multiplexing Capture if hashtag info (cmo_ids in samples) is present
        if _adt_prefixes and not _cmo_prefixes:
            _llm_samples = llm_data.get("samples", [])
            _has_cmo_ids = any(s.get("cmo_ids") for s in _llm_samples) if _llm_samples else False
            _has_hashtags = bool(hashtags) or bool(llm_data.get("hashtag_numbers"))
            if _has_cmo_ids or _has_hashtags:
                _cmo_prefixes = _adt_prefixes
                _adt_prefixes = []

    # Multiple GEX prefixes -> multi-GEM-well mode
    _is_multi_gem = len(_gex_prefixes) > 1

    if _is_multi_gem:
        # ================================================================
        # Multi-GEM-well UI
        # ================================================================
        st.markdown("**GEM wells** (auto-detected from FASTQ filenames)")
        st.caption(f"{len(_gex_prefixes)} GEM wells detected. `cellranger multi` will be run for each GEM well.")

        # Auto-match GEX <-> ADT / CMO
        _gem_wells = []
        for gp in _gex_prefixes:
            # Match CMO prefix by name (e.g. Control -> Control-Hashtag)
            matched_cmo = ""
            for cp in _cmo_prefixes:
                if cp.startswith(gp):
                    matched_cmo = cp
                    break
            # Match ADT prefix by name (e.g. MSC -> MSC_SP)
            matched_adt = ""
            for ap in _adt_prefixes:
                if ap.startswith(gp) or gp.startswith(ap.rsplit("_", 1)[0] if "_" in ap else ap):
                    matched_adt = ap
                    break
            _gem_wells.append({"gex": gp, "cmo": matched_cmo, "adt": matched_adt})

        # LLM / resolved samples
        resolved_samples = st.session_state.get(resolved_samples_key, [])
        llm_samples_raw = llm_data.get("samples", [])
        all_samples = resolved_samples if resolved_samples else llm_samples_raw

        # Debug: verify LLM results (always visible)
        with st.expander("Debug: LLM extraction status", expanded=not all_samples):
            st.write(f"**llm_key**: `{llm_key}`")
            st.write(f"**llm_data keys**: `{list(llm_data.keys()) if llm_data else '(empty)'}`")
            st.write(f"**resolved_samples**: {len(resolved_samples)} entries")
            st.write(f"**llm_samples_raw**: {len(llm_samples_raw)} entries")
            st.write(f"**all_samples**: {len(all_samples)} entries")
            if llm_data:
                st.json(llm_data)
            else:
                st.info("No LLM results yet. Upload a sample info file and click Extract with LLM.")

        # Check for mismatch between LLM gex_sample and FASTQ prefixes
        if all_samples:
            _llm_gex_names = sorted(set(
                s.get("gex_sample", "") for s in all_samples if s.get("gex_sample")
            ))
            _fastq_well_names = [w["gex"] for w in _gem_wells]
            _exact_match = set(_llm_gex_names) == set(_fastq_well_names)
            if not _exact_match and _llm_gex_names:
                st.warning(
                    f"GEM well names in the sample sheet do not match FASTQ prefixes.\n\n"
                    f"- **Sample sheet**: {', '.join(_llm_gex_names)}\n"
                    f"- **FASTQ prefix**: {', '.join(_fastq_well_names)}\n\n"
                    f"Attempting mapping by partial match or order. Please correct manually if incorrect."
                )

        all_multi_configs = []

        # Pre-set default values in session_state (value= is only effective on first render in forms)
        for wi, well in enumerate(_gem_wells):
            # Match LLM samples by gex_sample
            well_samples = [s for s in all_samples
                            if s.get("gex_sample", "") == well["gex"]]
            # Fallback: partial substring match
            if not well_samples and all_samples:
                well_samples = [s for s in all_samples
                                if s.get("gex_sample", "") and (
                                    s["gex_sample"].lower() in well["gex"].lower()
                                    or well["gex"].lower() in s["gex_sample"].lower()
                                )]
            # Fallback 2: order-based mapping
            if not well_samples and all_samples:
                _unique_gex = list(dict.fromkeys(
                    s.get("gex_sample", "") for s in all_samples if s.get("gex_sample")
                ))
                if len(_unique_gex) == len(_gem_wells) and wi < len(_unique_gex):
                    _target = _unique_gex[wi]
                    well_samples = [s for s in all_samples
                                    if s.get("gex_sample", "") == _target]
            _n_key = f"{key_prefix}_gem{wi}_n_samples"
            _desired_n = max(len(well_samples), 2)
            if well_samples and st.session_state.get(_n_key, 2) < _desired_n:
                st.session_state[_n_key] = _desired_n
            for si, ws in enumerate(well_samples):
                _sid_key = f"{key_prefix}_gem{wi}_sid_{si}"
                _cmo_key = f"{key_prefix}_gem{wi}_cmo_{si}"
                if not st.session_state.get(_sid_key):
                    st.session_state[_sid_key] = ws["sample_id"]
                if not st.session_state.get(_cmo_key):
                    cids = ws.get("cmo_ids", [])
                    st.session_state[_cmo_key] = cids if isinstance(cids, str) else ", ".join(cids) if cids else ""

        with st.form(key=f"{key_prefix}_gem_wells_form"):
            for wi, well in enumerate(_gem_wells):
                st.markdown(f"##### GEM well {wi+1}: {well['gex']}")
                col_g1, col_g2, col_g3 = st.columns(3)
                with col_g1:
                    w_gex = st.text_input(
                        "GEX fastq_id", value=well["gex"],
                        key=f"{key_prefix}_gem{wi}_gex_id")
                with col_g2:
                    w_adt = st.text_input(
                        "ADT fastq_id", value=well.get("adt", ""),
                        key=f"{key_prefix}_gem{wi}_adt_id")
                with col_g3:
                    w_cmo = st.text_input(
                        "CMO fastq_id", value=well["cmo"],
                        key=f"{key_prefix}_gem{wi}_cmo_id")

                st.caption(f"FASTQ dirs: {fastq_dir_default}")

                _n_key = f"{key_prefix}_gem{wi}_n_samples"
                well_samples = [s for s in all_samples
                                if s.get("gex_sample", "") == well["gex"]]
                if not well_samples and all_samples:
                    well_samples = [s for s in all_samples
                                    if s.get("gex_sample", "") and (
                                        s["gex_sample"].lower() in well["gex"].lower()
                                        or well["gex"].lower() in s["gex_sample"].lower()
                                    )]
                if not well_samples and all_samples:
                    _unique_gex = list(dict.fromkeys(
                        s.get("gex_sample", "") for s in all_samples if s.get("gex_sample")
                    ))
                    if len(_unique_gex) == len(_gem_wells) and wi < len(_unique_gex):
                        _target = _unique_gex[wi]
                        well_samples = [s for s in all_samples
                                        if s.get("gex_sample", "") == _target]

                st.caption("Sample ID / CMO IDs are editable (they will be used as output folder names).")
                n_ws = st.number_input(
                    "Samples in this well",
                    min_value=1, max_value=50,
                    value=max(len(well_samples), 2),
                    key=_n_key,
                )
                for si in range(int(n_ws)):
                    col_s1, col_s2 = st.columns(2)
                    d_sid = well_samples[si]["sample_id"] if si < len(well_samples) else ""
                    d_cmo = ""
                    if si < len(well_samples):
                        cids = well_samples[si].get("cmo_ids", [])
                        d_cmo = cids if isinstance(cids, str) else ", ".join(cids) if cids else ""
                    with col_s1:
                        st.text_input(
                            f"Sample {si+1} ID", value=d_sid,
                            key=f"{key_prefix}_gem{wi}_sid_{si}")
                    with col_s2:
                        st.text_input(
                            f"Sample {si+1} CMO IDs", value=d_cmo,
                            key=f"{key_prefix}_gem{wi}_cmo_{si}",
                            placeholder="e.g. HTO_B0305, HTO_B0306")
                if wi < len(_gem_wells) - 1:
                    st.divider()

            st.form_submit_button("Apply", type="primary")

        # After form submission, read values from session_state to build multi_configs
        for wi, well in enumerate(_gem_wells):
            w_gex = st.session_state.get(f"{key_prefix}_gem{wi}_gex_id", well["gex"])
            w_adt = st.session_state.get(f"{key_prefix}_gem{wi}_adt_id", well.get("adt", ""))
            w_cmo = st.session_state.get(f"{key_prefix}_gem{wi}_cmo_id", well["cmo"])
            _n_key = f"{key_prefix}_gem{wi}_n_samples"
            n_ws = int(st.session_state.get(_n_key, 2))

            w_samples = []
            for si in range(n_ws):
                sid = st.session_state.get(f"{key_prefix}_gem{wi}_sid_{si}", "")
                cmo = st.session_state.get(f"{key_prefix}_gem{wi}_cmo_{si}", "")
                if sid and cmo:
                    cmo_list = [c.strip() for c in cmo.split(",")]
                    w_samples.append({"sample_id": sid, "cmo_ids": cmo_list})

            if w_samples:
                st.success(f"GEM well {wi+1} ({w_gex}): {len(w_samples)} samples")

            w_libs = [
                {"fastq_id": w_gex, "fastqs": fastq_dir_default, "feature_types": "Gene Expression"},
            ]
            if w_adt:
                w_libs.append(
                    {"fastq_id": w_adt, "fastqs": fastq_dir_default, "feature_types": "Antibody Capture"},
                )
            if w_cmo:
                w_libs.append(
                    {"fastq_id": w_cmo, "fastqs": fastq_dir_default, "feature_types": "Multiplexing Capture"},
                )

                # multi_config entry
                all_multi_configs.append({
                    "id": w_gex,
                    "gex_fastq_id": w_gex,
                    "adt_fastq_id": w_adt,
                    "cmo_fastq_id": w_cmo,
                    "libraries": w_libs,
                    "samples": w_samples,
                    "hashtags": hashtags,
                })

        # Return multi_configs (processed on the submit side)
        return hashtags, all_multi_configs, []

    # ================================================================
    # Single GEM well UI (legacy)
    # ================================================================
    st.markdown("**Libraries**")

    _auto_gex_id = llm_libs.get("gex", {}).get("fastq_id", "")
    _auto_adt_id = llm_libs.get("adt", {}).get("fastq_id", "") if llm_libs.get("adt") else ""
    _auto_cmo_id = llm_libs.get("cmo", {}).get("fastq_id", "")
    if not _auto_gex_id and _gex_prefixes:
        _auto_gex_id = _gex_prefixes[0]
    if not _auto_adt_id and _adt_prefixes:
        _auto_adt_id = _adt_prefixes[0]
    if not _auto_cmo_id and _cmo_prefixes:
        _auto_cmo_id = _cmo_prefixes[0]
    # Set directly in session_state
    if _auto_gex_id and not st.session_state.get(f"{key_prefix}_gex_fastq_id"):
        st.session_state[f"{key_prefix}_gex_fastq_id"] = _auto_gex_id
    if _auto_adt_id and not st.session_state.get(f"{key_prefix}_adt_fastq_id"):
        st.session_state[f"{key_prefix}_adt_fastq_id"] = _auto_adt_id
    if _auto_cmo_id and not st.session_state.get(f"{key_prefix}_cmo_fastq_id"):
        st.session_state[f"{key_prefix}_cmo_fastq_id"] = _auto_cmo_id

    col_lib1, col_lib2, col_lib3 = st.columns(3)
    with col_lib1:
        gex_fastq_id = st.text_input(
            "GEX fastq_id",
            value=_auto_gex_id,
            key=f"{key_prefix}_gex_fastq_id",
            placeholder="e.g. GEX_sample",
        )
        gex_fastqs = st.text_input(
            "GEX FASTQ directory",
            value=fastq_dir_default,
            key=f"{key_prefix}_gex_fastqs",
        )
    with col_lib2:
        adt_fastq_id = st.text_input(
            "ADT fastq_id",
            value=_auto_adt_id,
            key=f"{key_prefix}_adt_fastq_id",
            placeholder="e.g. ADT_sample (optional)",
        )
        adt_fastqs = st.text_input(
            "ADT FASTQ directory",
            value=fastq_dir_default if _auto_adt_id else "",
            key=f"{key_prefix}_adt_fastqs",
        )
    with col_lib3:
        cmo_fastq_id = st.text_input(
            "CMO fastq_id",
            value=_auto_cmo_id,
            key=f"{key_prefix}_cmo_fastq_id",
            placeholder="e.g. CMO_sample",
        )
        cmo_fastqs = st.text_input(
            "CMO FASTQ directory",
            value=fastq_dir_default,
            key=f"{key_prefix}_cmo_fastqs",
        )

    # --- Samples (sample_id ↔ cmo_ids) ---
    st.markdown("**Samples (sample_id → cmo_ids)**")

    resolved_samples = st.session_state.get(resolved_samples_key, [])
    llm_samples_raw = llm_data.get("samples", [])
    default_samples = resolved_samples if resolved_samples else llm_samples_raw
    default_n_samples = max(len(default_samples), 2)
    n_samples = st.number_input(
        "Number of samples",
        min_value=1, max_value=50,
        value=default_n_samples,
        key=f"{key_prefix}_n_samples",
    )

    samples = []
    for i in range(int(n_samples)):
        col_s1, col_s2 = st.columns(2)
        default_sid = default_samples[i]["sample_id"] if i < len(default_samples) else ""
        default_cmo = ""
        if i < len(default_samples):
            cids = default_samples[i].get("cmo_ids", [])
            default_cmo = cids if isinstance(cids, str) else ", ".join(cids) if cids else ""
        with col_s1:
            sid = st.text_input(f"Sample {i+1} ID", value=default_sid, key=f"{key_prefix}_sid_{i}")
        with col_s2:
            cmo = st.text_input(f"Sample {i+1} CMO IDs", value=default_cmo, key=f"{key_prefix}_cmo_{i}",
                                placeholder="e.g. HTO_B0305, HTO_B0306")
        if sid and cmo:
            cmo_list = [c.strip() for c in cmo.split(",")]
            samples.append({"sample_id": sid, "cmo_ids": cmo_list})

    if samples:
        st.success(f"{len(samples)} samples configured")

    # --- multi.csv preview ---
    if gex_fastq_id and gex_fastqs:
        libraries = [
            {"fastq_id": gex_fastq_id, "fastqs": gex_fastqs, "feature_types": "Gene Expression"},
        ]
        if adt_fastq_id and adt_fastqs:
            libraries.append(
                {"fastq_id": adt_fastq_id, "fastqs": adt_fastqs, "feature_types": "Antibody Capture"},
            )
        if cmo_fastq_id and cmo_fastqs:
            libraries.append(
                {"fastq_id": cmo_fastq_id, "fastqs": cmo_fastqs, "feature_types": "Multiplexing Capture"},
            )

        preview_config = {
            "reference": "(selected reference)",
            "create_bam": True,
            "hashtags": hashtags,
            "feature_ref_path": "feature_reference.csv",
            "libraries": libraries,
            "samples": samples,
        }
        preview_csv = build_multi_csv(preview_config)
        with st.expander("multi.csv preview", expanded=False):
            st.code(preview_csv, language="text")

        return hashtags, libraries, samples

    return hashtags, [], samples


# =============================================================================
#  UI
# =============================================================================

st.title("FASTQ to scRNA (Cell Ranger)")

tab_local, tab_azenta, tab_monitor = st.tabs(["Local NovaSeq", "Azenta / Dropbox / Upload", "Job Monitor"])

# =============================================================================
#  Tab 1: Local NovaSeq
# =============================================================================

with tab_local:
    local_input_source = st.radio(
        "Input source",
        ["SFTP (BCL run)", "Local BCL run"],
        index=0,
        key="scrna_local_input_source",
        horizontal=True,
    )

    if local_input_source == "SFTP (BCL run)":
        st.subheader("1. Select BCL run from Aging server")
        st.caption("Navigate: Nova-seq -> select a {run} directory")

        # SFTP browser
        remote_dir = st.text_input(
            "Remote directory path",
            value=st.session_state.get("scrna_sftp_remote_dir", "/data/sequencing"),
            key="scrna_remote_dir_input",
        )

        if st.button("Browse", type="primary", key="scrna_sftp_browse"):
            st.session_state.scrna_sftp_remote_dir = remote_dir

        if "scrna_sftp_remote_dir" in st.session_state:
            ssh, sftp = get_sftp_connection()
            if sftp is not None:
                current_dir = st.session_state.scrna_sftp_remote_dir
                st.caption(f"Current: `{current_dir}`")
                dirs, files = list_remote_dir(sftp, current_dir)

                if dirs:
                    selected_run = st.selectbox(
                        "Run directory",
                        ["-- select --"] + dirs,
                        key="scrna_subdir_nav",
                    )

                    if selected_run != "-- select --":
                        run_path = f"{current_dir}/{selected_run}"
                        st.session_state.scrna_selected_run = run_path
                        st.success(f"Selected: `{run_path}`")

    else:  # Local BCL run
        st.subheader("1. Select local BCL run directory")
        st.caption("Specify a local BCL run directory directly (no SFTP download required).")

        local_bcl_path = st.text_input(
            "BCL run directory path",
            value=st.session_state.get("scrna_local_bcl_path", ""),
            key="scrna_local_bcl_path",
            placeholder="/data/sequencing/run_directory/...",
        )

        if local_bcl_path:
            if os.path.isdir(local_bcl_path):
                # Basic validation of BCL run directory
                has_rta = os.path.exists(os.path.join(local_bcl_path, "RTAComplete.txt"))
                has_runinfo = os.path.exists(os.path.join(local_bcl_path, "RunInfo.xml"))
                if has_rta or has_runinfo:
                    st.success(f"BCL run directory: `{local_bcl_path}`")
                    st.session_state.scrna_selected_run = local_bcl_path
                else:
                    st.warning(f"RTAComplete.txt / RunInfo.xml not found. Please verify this is a BCL run directory.")
                    # Warning only -- user can proceed if they confirm it's correct
                    st.session_state.scrna_selected_run = local_bcl_path
            else:
                st.error(f"Directory does not exist: `{local_bcl_path}`")

    st.markdown("---")

    # Extract sample sheet information (done before mode selection)
    st.subheader("2. Sample Sheet & Mode")

    # Automatic sample sheet extraction via LLM
    st.caption("Automatically extract sample sheet information from an Excel file using LLM and infer the mode (multi/count).")
    mkfastq_info_file = st.file_uploader(
        "Sample info file (Excel, CSV, TSV, TXT)",
        type=["xlsx", "xls", "csv", "tsv", "txt"],
        key="local_mkfastq_info_file",
    )

    mkfastq_llm_key = "local_mkfastq_llm_result"
    if mkfastq_info_file is not None:
        if st.button("Extract with LLM", type="primary", key="local_mkfastq_extract_llm"):
            content = read_uploaded_file_content(mkfastq_info_file)
            if content:
                with st.spinner("Extracting sample info with LLM (integrated mkfastq + multi)..."):
                    result = extract_sample_info_with_llm(content, mkfastq_info_file.name, extraction_type="local")
                    if result:
                        st.session_state[mkfastq_llm_key] = result
                        # Clear old GEM well widget keys (prevent stale values from persisting)
                        keys_to_clear = [k for k in st.session_state if k.startswith("gem_")]
                        for k in keys_to_clear:
                            del st.session_state[k]
                        # Also clear local_gem_groups
                        st.session_state.pop("local_gem_groups", None)
                        # Mode estimation
                        detected_mode_val = result.get("mode", "")
                        has_cmo = bool(result.get("cmo_sample_sheet"))
                        if detected_mode_val == "multi" or has_cmo:
                            detected_mode = "mkfastq + multi"
                        else:
                            detected_mode = "mkfastq + count"
                        st.session_state["local_detected_mode"] = detected_mode
                        # Also directly update the Mode selectbox
                        st.session_state["local_mode"] = detected_mode

                        # Directly update the GEX/CMO sample sheet widgets
                        gex_text = "Lane,Sample,Index\n"
                        for row in result.get("gex_sample_sheet", []):
                            gex_text += f"{row.get('lane','*')},{row.get('sample','')},{row.get('index','')}\n"
                        st.session_state["local_gex_sheet_text"] = gex_text

                        if has_cmo:
                            cmo_text = "Lane,Sample,Index\n"
                            for row in result.get("cmo_sample_sheet", []):
                                cmo_text += f"{row.get('lane','*')},{row.get('sample','')},{row.get('index','')}\n"
                            st.session_state["local_cmo_sheet_text"] = cmo_text

                        # Also save multi info (hashtag, samples)
                        if has_cmo and result.get("hashtag_numbers"):
                            # Save with the same keys for multi_mode_ui
                            st.session_state["local_llm_result"] = result
                            resolved = resolve_hashtags_from_llm(result)
                            st.session_state["local_resolved_hashtags"] = resolved
                            resolved_samples = resolve_samples_from_llm(result, resolved)
                            st.session_state["local_resolved_samples"] = resolved_samples
                            # Directly update the Hashtag CSV widget
                            if resolved:
                                ht_csv = "id,name,read,pattern,sequence,feature_type\n"
                                for h in resolved:
                                    ht_csv += f"{h['id']},{h['name']},{h['read']},{h['pattern']},{h['sequence']},{h['feature_type']}\n"
                                st.session_state["local_hashtag_csv"] = ht_csv
                            # Directly update the Samples widget
                            if resolved_samples:
                                st.session_state["local_n_samples"] = len(resolved_samples)
                                for i, s in enumerate(resolved_samples):
                                    st.session_state[f"local_sid_{i}"] = s["sample_id"]
                                    cids = s.get("cmo_ids", [])
                                    cmo_str = cids if isinstance(cids, str) else ", ".join(cids)
                                    st.session_state[f"local_cmo_{i}"] = cmo_str
                            ht_nums = result.get("hashtag_numbers", [])
                            species = result.get("species", "?")
                            st.success(f"Extraction complete! -> **{detected_mode}** | Hashtag {ht_nums} ({species}), {len(resolved)} entries resolved")
                        else:
                            st.success(f"Extraction complete! -> **{detected_mode}**")

    # Mode selection (using default values estimated from LLM results)
    mode_options = ["mkfastq + multi", "mkfastq + count", "mkfastq only"]
    detected = st.session_state.get("local_detected_mode", None)
    default_idx = mode_options.index(detected) if detected in mode_options else 0
    local_mode = st.selectbox(
        "Mode",
        mode_options,
        index=default_idx,
        key="local_mode",
    )

    # mkfastq settings
    st.markdown("#### mkfastq settings")

    # Build default values from LLM results
    mkfastq_llm_data = st.session_state.get(mkfastq_llm_key, {})
    gex_default = "Lane,Sample,Index\n"
    cmo_default = "Lane,Sample,Index\n"
    if mkfastq_llm_data.get("gex_sample_sheet"):
        for row in mkfastq_llm_data["gex_sample_sheet"]:
            gex_default += f"{row.get('lane','*')},{row.get('sample','')},{row.get('index','')}\n"
    else:
        gex_default += "*,SampleGEX,SI-TT-H12\n"
    if mkfastq_llm_data.get("cmo_sample_sheet"):
        for row in mkfastq_llm_data["cmo_sample_sheet"]:
            cmo_default += f"{row.get('lane','*')},{row.get('sample','')},{row.get('index','')}\n"
    else:
        cmo_default += "*,SampleCMO,SI-NT-H10\n"

    st.caption("GEX Sample Sheet CSV (text input or file upload)")
    gex_sheet_method = st.radio(
        "GEX Sample Sheet input method",
        ["Text input", "File upload"],
        key="local_gex_sheet_method",
        horizontal=True,
    )

    if gex_sheet_method == "Text input":
        gex_sample_sheet = st.text_area(
            "GEX Sample Sheet CSV",
            value=gex_default,
            height=100,
            key="local_gex_sheet_text",
        )
    else:
        gex_sheet_file = st.file_uploader("GEX Sample Sheet CSV", type=["csv"], key="local_gex_sheet_file")
        gex_sample_sheet = ""
        if gex_sheet_file:
            gex_sample_sheet = gex_sheet_file.read().decode("utf-8")
            gex_sheet_file.seek(0)

    cmo_sample_sheet = ""
    if "multi" in local_mode:
        st.caption("CMO Sample Sheet (hashtag library index; required for multi mode)")
        cmo_sheet_method = st.radio(
            "CMO Sample Sheet input method",
            ["Text input", "File upload"],
            key="local_cmo_sheet_method",
            horizontal=True,
        )
        if cmo_sheet_method == "Text input":
            cmo_sample_sheet = st.text_area(
                "CMO Sample Sheet CSV",
                value=cmo_default,
                height=100,
                key="local_cmo_sheet_text",
            )
        else:
            cmo_sheet_file = st.file_uploader("CMO Sample Sheet CSV", type=["csv"], key="local_cmo_sheet_file")
            if cmo_sheet_file:
                cmo_sample_sheet = cmo_sheet_file.read().decode("utf-8")
                cmo_sheet_file.seek(0)

    mkfastq_cores = st.number_input(
        "mkfastq --localcores",
        min_value=1, max_value=256,
        value=min(128, os.cpu_count() or 64),
        key="local_mkfastq_cores",
    )

    # count/multi settings
    if local_mode != "mkfastq only":
        st.markdown("---")
        st.subheader("3. Cell Ranger Settings")
        local_ref, local_ref_name, local_create_bam, local_also_run_count, local_cores, local_mem, local_run_annotate = common_cellranger_settings("local")

        if "count" in local_mode:
            # mkfastq + count: detect multiple samples from GEX sample sheet
            gex_lines = [l.strip() for l in gex_sample_sheet.strip().split("\n") if l.strip()]
            gex_samples_from_sheet = []
            if len(gex_lines) > 1:
                for line in gex_lines[1:]:
                    parts = [p.strip() for p in line.split(",")]
                    if len(parts) >= 2 and parts[1]:
                        gex_samples_from_sheet.append(parts[1])

            if len(gex_samples_from_sheet) > 1:
                st.markdown("#### Cell Ranger count settings (multiple samples)")
                st.caption(f"{len(gex_samples_from_sheet)} samples detected in the GEX sample sheet. `cellranger count` will be run for each sample.")

                local_count_configs = []
                for i, sname in enumerate(gex_samples_from_sheet):
                    col1, col2 = st.columns([3, 3])
                    with col1:
                        edited_name = st.text_input(
                            f"Sample {i+1} name",
                            value=sname,
                            key=f"local_count_sample_{i}",
                        )
                    with col2:
                        edited_id = st.text_input(
                            f"Sample {i+1} output ID",
                            value=sname,
                            key=f"local_count_id_{i}",
                        )
                    local_count_configs.append({
                        "sample_name": edited_name,
                        "id": edited_id,
                    })
                st.success(f"{len(local_count_configs)} count runs will be executed")
            else:
                local_count_sample, local_count_fastqs = count_mode_ui("local")
                local_count_configs = None
        elif "multi" in local_mode:
            # Pre-determine the number of GEM wells
            llm_data = st.session_state.get("local_llm_result", {})
            gex_sheets = llm_data.get("gex_sample_sheet", [])
            is_multi_gem = len(gex_sheets) > 1

            local_hashtags, local_libraries, local_samples = multi_mode_ui("local", multi_gem_well=is_multi_gem)

            # GEM well group display -- when LLM extracted from sample sheet
            cmo_sheets = llm_data.get("cmo_sample_sheet", [])
            llm_samples = llm_data.get("samples", [])

            if len(gex_sheets) > 1 and llm_samples:
                st.markdown("---")
                st.markdown("#### GEM well groups (auto-detected from sample sheet)")
                st.caption("Samples sharing the same Dual Index are processed in one GEM well. `cellranger multi` is run for each GEM well.")

                # Group by GEX index
                gem_groups = {}
                for gs in gex_sheets:
                    idx = gs["index"]
                    if idx not in gem_groups:
                        gem_groups[idx] = {
                            "gex_index": idx,
                            "gex_sample": gs["sample"],
                            "cmo_index": "",
                            "cmo_sample": "",
                            "sample_ids": [],
                        }

                # Match CMO index
                for cs in cmo_sheets:
                    # Match CMO to GEX: use entries at the same position
                    pass  # Match by sample below

                # unique GEX indices
                unique_gex = list(gem_groups.keys())
                # CMO indices (same order)
                unique_cmo = list(dict.fromkeys(cs["index"] for cs in cmo_sheets))
                for i, gex_idx in enumerate(unique_gex):
                    if i < len(unique_cmo):
                        gem_groups[gex_idx]["cmo_index"] = unique_cmo[i]
                        # Get CMO sample name
                        for cs in cmo_sheets:
                            if cs["index"] == unique_cmo[i]:
                                gem_groups[gex_idx]["cmo_sample"] = cs["sample"]
                                break

                # Distribute LLM samples to groups
                resolved_ht = st.session_state.get("local_resolved_hashtags", [])
                resolved_smp = st.session_state.get("local_resolved_samples", [])
                if not resolved_smp:
                    resolved_smp = resolve_samples_from_llm(llm_data, resolved_ht)

                # GEX sample name -> GEX index mapping
                gex_sample_to_idx = {gs["sample"]: gs["index"] for gs in gex_sheets}

                # Strategy 1: Distribute by gex_sample if LLM returned it
                assigned_by_gex_sample = False
                for s in resolved_smp:
                    gs = s.get("gex_sample", "")
                    if gs and gs in gex_sample_to_idx:
                        idx = gex_sample_to_idx[gs]
                        gem_groups[idx]["sample_ids"].append(s)
                        assigned_by_gex_sample = True
                    elif gs:
                        # Also try partial match
                        matched = False
                        for name, idx in gex_sample_to_idx.items():
                            if gs in name or name in gs:
                                gem_groups[idx]["sample_ids"].append(s)
                                assigned_by_gex_sample = True
                                matched = True
                                break
                        if not matched:
                            assigned_by_gex_sample = False
                            break
                    else:
                        assigned_by_gex_sample = False
                        break

                # Strategy 2: If no gex_sample, use 1:1 mapping when gex_sample_sheet count matches samples count
                if not assigned_by_gex_sample:
                    # Reset
                    for idx in gem_groups:
                        gem_groups[idx]["sample_ids"] = []
                    if len(gex_sheets) == len(resolved_smp):
                        for i, gs in enumerate(gex_sheets):
                            idx = gs["index"]
                            gem_groups[idx]["sample_ids"].append(resolved_smp[i])
                    elif len(unique_gex) > 1:
                        # Fallback: display as unassigned sample list (prompt user for manual adjustment)
                        st.warning("Could not automatically determine sample-to-GEM-well assignment. Please verify manually.")
                        # Put all samples in the first well for now
                        for s in resolved_smp:
                            gem_groups[unique_gex[0]]["sample_ids"].append(s)

                # Per-GEM-well hashtag filtering
                # Retain only the cmo_ids used in each well as the hashtag list
                all_hashtag_ids = {h["id"] for h in local_hashtags}
                for gex_idx, grp in gem_groups.items():
                    used_cmo_ids = set()
                    for s in grp["sample_ids"]:
                        cids = s.get("cmo_ids", [])
                        if isinstance(cids, list):
                            used_cmo_ids.update(cids)
                        else:
                            used_cmo_ids.add(cids)
                    # Filter to only the hashtags used in this well
                    well_custom_ht = [h for h in local_hashtags if h["id"] in used_cmo_ids]
                    # Set hashtags to empty for wells using Built-in CMO IDs (CMO301-CMO312)
                    has_builtin_cmo = any(
                        re.match(r'^CMO\d+$', cid) for cid in used_cmo_ids
                    )
                    grp["hashtags"] = well_custom_ht
                    grp["uses_builtin_cmo"] = has_builtin_cmo and not well_custom_ht

                # UI display
                for gex_idx, grp in gem_groups.items():
                    grp_id = gex_idx.replace("SI-TT-", "")
                    n_smp = len(grp["sample_ids"])
                    n_ht = len(grp.get("hashtags", []))
                    uses_builtin = grp.get("uses_builtin_cmo", False)
                    ht_label = "Built-in CMO" if uses_builtin else f"{n_ht} hashtags"
                    with st.expander(f"GEM well {grp_id}: {grp['gex_sample']} ({n_smp} samples, {ht_label})", expanded=True):
                        st.text(f"GEX: {grp['gex_sample']} ({gex_idx})")
                        st.text(f"CMO: {grp['cmo_sample']} ({grp['cmo_index']})")
                        if uses_builtin:
                            # Display built-in CMO IDs
                            builtin_ids = set()
                            for s in grp["sample_ids"]:
                                cids = s.get("cmo_ids", [])
                                for cid in (cids if isinstance(cids, list) else [cids]):
                                    if re.match(r'^CMO\d+$', cid):
                                        builtin_ids.add(cid)
                            st.text(f"CMO IDs: {', '.join(sorted(builtin_ids))} (Cell Ranger built-in)")
                        elif grp.get("hashtags"):
                            ht_names = [f"{h['id']} ({h['name']})" for h in grp["hashtags"]]
                            st.text(f"Hashtags: {', '.join(ht_names)}")
                        if grp["sample_ids"]:
                            for si, s in enumerate(grp["sample_ids"]):
                                cmo_default = ", ".join(s["cmo_ids"]) if isinstance(s.get("cmo_ids"), list) else str(s.get("cmo_ids", ""))
                                col_sid, col_grp, col_cmo = st.columns([2, 2, 3])
                                with col_sid:
                                    new_sid = st.text_input(
                                        f"Sample {si+1} ID",
                                        value=s["sample_id"],
                                        key=f"gem_{gex_idx}_sid_{si}",
                                    )
                                with col_grp:
                                    new_grp = st.text_input(
                                        f"Sample {si+1} Group",
                                        value=s.get("group", ""),
                                        key=f"gem_{gex_idx}_grp_{si}",
                                    )
                                with col_cmo:
                                    new_cmo = st.text_input(
                                        f"Sample {si+1} CMO IDs",
                                        value=cmo_default,
                                        key=f"gem_{gex_idx}_cmo_{si}",
                                    )
                                # Apply user edits
                                s["sample_id"] = new_sid
                                s["group"] = new_grp
                                s["cmo_ids"] = [c.strip() for c in new_cmo.split(",") if c.strip()]

                st.session_state["local_gem_groups"] = gem_groups

    # Output destination
    st.markdown("---")
    st.subheader("4. Output" if local_mode != "mkfastq only" else "3. Output")

    local_work_dir = st.text_input(
        "Work directory",
        value=str(Path.home() / "analysis"),
        key="local_work_dir",
        placeholder="/data/analysis/project_scRNA",
    )

    if local_mode != "mkfastq only":
        local_scala_dest = scala_dir_browser("local")

    # Submit
    st.markdown("---")
    step_num = "5" if local_mode != "mkfastq only" else "4"
    st.subheader(f"{step_num}. Submit")

    has_run = "scrna_selected_run" in st.session_state
    if not has_run:
        st.warning("First select a run directory.")
    elif not gex_sample_sheet.strip():
        st.warning("Please enter the GEX Sample Sheet.")
    elif st.button("Submit Job", type="primary", key="local_submit"):
        config = {
            "data_source": "local_bcl" if local_input_source == "Local BCL run" else "sftp",
            "mode": local_mode.split(" + ")[-1] if local_mode != "mkfastq only" else "mkfastq",
            "work_dir": local_work_dir,
            "mkfastq": {
                "run_dir": st.session_state.scrna_selected_run,
                "gex_sample_sheet": gex_sample_sheet,
                "cmo_sample_sheet": cmo_sample_sheet,
                "localcores": int(mkfastq_cores),
            },
        }

        if local_input_source == "SFTP (BCL run)":
            cred = st.secrets["aging_sftp"]
            config["sftp_config"] = {
                "host": cred["host"],
                "username": cred["username"],
                "password": cred["password"],
            }

        if local_mode != "mkfastq only":
            config["reference"] = local_ref
            config["reference_name"] = local_ref_name
            config["create_bam"] = local_create_bam
            config["also_run_count"] = local_also_run_count
            config["run_annotate"] = local_run_annotate
            config["localcores"] = local_cores
            config["localmem"] = local_mem
            config["scala_dest"] = local_scala_dest

            if "count" in local_mode:
                if local_count_configs and len(local_count_configs) > 1:
                    # Multiple sample count
                    config["count"] = {
                        "count_configs": local_count_configs,
                    }
                else:
                    config["count"] = {
                        "sample_name": local_count_sample,
                        "fastqs_dir": local_count_fastqs,
                    }
            elif "multi" in local_mode:
                gem_groups = st.session_state.get("local_gem_groups")
                if gem_groups and len(gem_groups) > 1:
                    # Multiple GEM wells -> build multi_configs (per-well hashtags)
                    multi_configs = []
                    for gex_idx, grp in gem_groups.items():
                        # GEM well ID: use GEX sample name (e.g. EC_8w, E17_EC)
                        well_id = grp["gex_sample"]
                        # Use per-GEM-well hashtags (filtered)
                        well_hashtags = grp.get("hashtags", local_hashtags)
                        mc = {
                            "id": well_id,
                            "hashtags": well_hashtags,
                            "libraries": [
                                {"fastq_id": grp["gex_sample"], "fastqs": "(auto)", "feature_types": "Gene Expression"},
                                {"fastq_id": grp["cmo_sample"], "fastqs": "(auto)", "feature_types": "Multiplexing Capture"},
                            ],
                            "samples": grp["sample_ids"],
                        }
                        multi_configs.append(mc)
                    config["multi"] = {"multi_configs": multi_configs}
                else:
                    # Single GEM well
                    config["multi"] = {
                        "hashtags": local_hashtags,
                        "libraries": local_libraries,
                        "samples": local_samples,
                    }

        job_id = submit_job(config)
        st.success(f"Job submitted! ID: {job_id}")
        st.info("Processing continues even if you close the browser. Check progress in the Job Monitor tab.")


# =============================================================================
#  Tab 2: Azenta
# =============================================================================

with tab_azenta:
    st.subheader("1. FASTQ Acquisition")

    azenta_source = st.radio(
        "Data source",
        ["Azenta OSS", "SFTP", "Upload", "Dropbox", "Google Drive", "Local Folder"],
        key="azenta_source",
        horizontal=True,
    )

    fastq_filenames = []
    azenta_data_config = {}

    if azenta_source == "Azenta OSS":
        st.caption("Enter the OSS information from the GENEWIZ/Azenta email.")
        oss_path = st.text_input(
            "OSS Path",
            key="scrna_oss_path",
            placeholder="oss://gwjapan/2026.2/...",
        )
        col_oss1, col_oss2 = st.columns(2)
        with col_oss1:
            oss_key_id = st.text_input("AccessKeyId", key="scrna_oss_key_id")
        with col_oss2:
            oss_key_secret = st.text_input("AccessKeySecret", type="password", key="scrna_oss_key_secret")

        oss_region = st.selectbox(
            "Region",
            ["oss-ap-northeast-1 (Tokyo)", "oss-ap-southeast-1 (Singapore)", "oss-cn-shanghai"],
            key="scrna_oss_region",
        )
        oss_endpoint = oss_region.split(" ")[0] + ".aliyuncs.com"

        if st.button("List files", type="primary", key="scrna_oss_list"):
            if not oss_key_id or not oss_key_secret:
                st.error("Please enter AccessKeyId and AccessKeySecret.")
            else:
                with st.spinner("Listing OSS files..."):
                    path = oss_path.rstrip("/") + "/"
                    files, err = list_oss_files(path, oss_key_id, oss_key_secret, oss_endpoint)
                    if err:
                        st.error(f"OSS error: {err}")
                    elif files:
                        st.session_state.scrna_oss_files = files
                    else:
                        st.warning("No files found.")

        if "scrna_oss_files" in st.session_state:
            files = st.session_state.scrna_oss_files
            fastq_files = [f for f in files if f["name"].endswith((".fastq.gz", ".fq.gz"))]
            st.success(f"**{len(fastq_files)} FASTQ files found**")
            for f in fastq_files:
                st.text(f"  {f['name']} ({format_size(f['size'])})")
            total_size = sum(f["size"] for f in fastq_files)
            st.info(f"Total: **{format_size(total_size)}**")
            fastq_filenames = [f["name"] for f in fastq_files]

            azenta_data_config = {
                "data_source": "oss",
                "oss_path": oss_path.rstrip("/") + "/",
                "oss_access_key_id": oss_key_id,
                "oss_access_key_secret": oss_key_secret,
                "oss_endpoint": oss_endpoint,
                "oss_files": [{"name": f["name"], "path": f["path"], "size": f["size"]} for f in fastq_files],
            }

    elif azenta_source == "SFTP":
        st.caption("Download FASTQ files from an SFTP server.")
        sftp_host = st.text_input("SFTP Host", key="scrna_az_sftp_host")
        col_s1, col_s2 = st.columns(2)
        with col_s1:
            sftp_user = st.text_input("Username", key="scrna_az_sftp_user")
        with col_s2:
            sftp_pass = st.text_input("Password", type="password", key="scrna_az_sftp_pass")
        sftp_dir = st.text_input("Remote directory", key="scrna_az_sftp_dir", placeholder="/data/fastq/")

        if st.button("Browse SFTP", type="primary", key="scrna_az_sftp_browse"):
            if sftp_host and sftp_user and sftp_pass and sftp_dir:
                try:
                    ssh = paramiko.SSHClient()
                    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
                    ssh.connect(sftp_host, username=sftp_user, password=sftp_pass)
                    sftp_client = ssh.open_sftp()
                    items = sftp_client.listdir_attr(sftp_dir)
                    fq_files = []
                    for item in items:
                        if not stat_module.S_ISDIR(item.st_mode):
                            if item.filename.endswith((".fastq.gz", ".fq.gz")):
                                fq_files.append({"name": item.filename, "size": item.st_size,
                                                 "path": f"{sftp_dir.rstrip('/')}/{item.filename}"})
                    st.session_state.scrna_az_sftp_files = fq_files
                    sftp_client.close()
                    ssh.close()
                except Exception as e:
                    st.error(f"SFTP error: {e}")

        if "scrna_az_sftp_files" in st.session_state:
            fq_files = st.session_state.scrna_az_sftp_files
            st.success(f"**{len(fq_files)} FASTQ files found**")
            for f in fq_files:
                st.text(f"  {f['name']} ({format_size(f['size'])})")
            fastq_filenames = [f["name"] for f in fq_files]

            azenta_data_config = {
                "data_source": "sftp",
                "sftp_config": {
                    "host": sftp_host, "username": sftp_user, "password": sftp_pass,
                    "remote_dir": sftp_dir,
                },
                "sftp_files": fq_files,
            }

    elif azenta_source == "Upload":
        st.caption("Upload local FASTQ files.")
        uploaded_fqs = st.file_uploader(
            "FASTQ files (.fastq.gz / .fq.gz)",
            type=["gz"],
            accept_multiple_files=True,
            key="scrna_upload_fastq",
        )
        if uploaded_fqs:
            fq_files = [f for f in uploaded_fqs if f.name.endswith((".fastq.gz", ".fq.gz"))]
            st.success(f"**{len(fq_files)} FASTQ files uploaded**")
            for f in fq_files:
                st.text(f"  {f.name} ({format_size(f.size)})")
            fastq_filenames = [f.name for f in fq_files]
            st.session_state.scrna_upload_files = fq_files

            azenta_data_config = {
                "data_source": "upload",
                "upload_files": [f.name for f in fq_files],
            }

    elif azenta_source == "Dropbox":
        st.caption("Download a ZIP from a Dropbox shared link, extract, and configure parameters.")

        # --- Auto-restore completed download jobs ---
        if "scrna_dropbox_downloaded" not in st.session_state:
            for _jname in sorted(os.listdir(JOBS_DIR), reverse=True) if os.path.isdir(JOBS_DIR) else []:
                _jdir = os.path.join(JOBS_DIR, _jname)
                _sp = os.path.join(_jdir, "status.json")
                _cp = os.path.join(_jdir, "config.json")
                if os.path.isfile(_sp) and os.path.isfile(_cp):
                    try:
                        with open(_sp) as _f:
                            _js = json.load(_f)
                        with open(_cp) as _f:
                            _jc = json.load(_f)
                        if (_jc.get("data_source") == "dropbox"
                                and _js.get("status") == "download_complete"
                                and _js.get("fastq_dir")
                                and os.path.isdir(_js["fastq_dir"])):
                            st.session_state.scrna_dropbox_downloaded = {
                                "fastq_dir": _js["fastq_dir"],
                                "files": _js.get("extracted_files", []),
                            }
                            break
                    except Exception:
                        pass

        dropbox_url = st.text_input(
            "Dropbox shared link URL",
            placeholder="https://www.dropbox.com/scl/fo/...",
            key="scrna_dropbox_url",
        )
        fastq_filenames = []

        # Show file list even without URL if already downloaded
        dl_info = st.session_state.get("scrna_dropbox_downloaded")

        if dl_info:
            # === Phase 2: Download complete -> display file list ===
            files = dl_info["files"]
            fastq_dir = dl_info["fastq_dir"]
            st.success(f"**{len(files)} FASTQ files** extracted to `{fastq_dir}`")
            with st.expander("Extracted files", expanded=False):
                for f in files:
                    st.text(f"  {f['name']} ({format_size(f['size'])})")
            fastq_filenames = [f["name"] for f in files]
            azenta_data_config = {
                "data_source": "dropbox",
                "dropbox_url": dropbox_url,
                "fastq_dir": fastq_dir,
            }
            if st.button("Reset (re-download)", key="scrna_dropbox_reset"):
                st.session_state.pop("scrna_dropbox_downloaded", None)
                st.session_state.pop("scrna_dropbox_job_id", None)
                st.rerun()

        elif dropbox_url:
            if "dropbox.com" not in dropbox_url:
                st.error("This is not a Dropbox URL.")
            else:
                # === Phase 1: Download ===
                dl_job_id = st.session_state.get("scrna_dropbox_job_id")

                if dl_job_id:
                    # Download progress display
                    job_dir = os.path.join(JOBS_DIR, dl_job_id)
                    status_path = os.path.join(job_dir, "status.json")
                    if os.path.exists(status_path):
                        with open(status_path) as _f:
                            dl_status = json.load(_f)
                        dl_st = dl_status.get("status", "unknown")
                        if dl_st == "download_complete":
                            st.session_state.scrna_dropbox_downloaded = {
                                "fastq_dir": dl_status["fastq_dir"],
                                "files": dl_status["extracted_files"],
                            }
                            st.rerun()
                        elif dl_st == "error":
                            st.error(f"Download failed: {dl_status.get('error')}")
                            if st.button("Retry", key="scrna_dropbox_retry"):
                                st.session_state.pop("scrna_dropbox_job_id", None)
                                st.rerun()
                        else:
                            pct = dl_status.get("progress_pct", 0)
                            step = dl_status.get("step", "Downloading...")
                            st.progress(pct / 100, text=step)
                            st.info("Downloading... Refresh the page to check progress.")

                            @st.fragment(run_every=5)
                            def _dropbox_dl_auto_refresh():
                                """Auto-refresh during download"""
                                _sp = os.path.join(JOBS_DIR, dl_job_id, "status.json")
                                if os.path.exists(_sp):
                                    with open(_sp) as __f:
                                        _s = json.load(__f)
                                    if _s.get("status") in ("download_complete", "error"):
                                        st.rerun()
                            _dropbox_dl_auto_refresh()
                    else:
                        st.warning("Job status not found. Waiting...")
                else:
                    # Download start UI
                    dropbox_work_dir = st.text_input(
                        "Work directory (download destination)",
                        value=str(Path.home() / "analysis"),
                        key="dropbox_work_dir",
                        placeholder="/data/analysis/project_scRNA",
                    )
                    if st.button("Download & Extract", type="primary", key="scrna_dropbox_dl"):
                        if not dropbox_work_dir:
                            st.error("Please specify a work directory.")
                        else:
                            dl_config = {
                                "data_source": "dropbox",
                                "dropbox_url": dropbox_url,
                                "work_dir": dropbox_work_dir,
                                "phase": "download_only",
                            }
                            job_id = submit_job(dl_config)
                            st.session_state.scrna_dropbox_job_id = job_id
                            st.rerun()

    elif azenta_source == "Google Drive":
        st.caption("Download files from a Google Drive shared folder, then configure parameters.")

        # --- Auto-restore completed download jobs ---
        if "scrna_gdrive_downloaded" not in st.session_state:
            for _jname in sorted(os.listdir(JOBS_DIR), reverse=True) if os.path.isdir(JOBS_DIR) else []:
                _jdir = os.path.join(JOBS_DIR, _jname)
                _sp = os.path.join(_jdir, "status.json")
                _cp = os.path.join(_jdir, "config.json")
                if os.path.isfile(_sp) and os.path.isfile(_cp):
                    try:
                        with open(_sp) as _f:
                            _js = json.load(_f)
                        with open(_cp) as _f:
                            _jc = json.load(_f)
                        if (_jc.get("data_source") == "gdrive"
                                and _js.get("status") == "download_complete"
                                and _js.get("fastq_dir")
                                and os.path.isdir(_js["fastq_dir"])):
                            st.session_state.scrna_gdrive_downloaded = {
                                "fastq_dir": _js["fastq_dir"],
                                "files": _js.get("extracted_files", []),
                            }
                            break
                    except Exception:
                        pass

        gdrive_url = st.text_input(
            "Google Drive shared folder URL",
            placeholder="https://drive.google.com/drive/folders/...",
            key="scrna_gdrive_url",
        )
        fastq_filenames = []

        dl_info_gd = st.session_state.get("scrna_gdrive_downloaded")

        if dl_info_gd:
            # === Phase 2: Download complete -> display file list ===
            files = dl_info_gd["files"]
            fastq_dir = dl_info_gd["fastq_dir"]
            st.success(f"**{len(files)} FASTQ files** downloaded to `{fastq_dir}`")
            with st.expander("Downloaded files", expanded=False):
                for f in files:
                    st.text(f"  {f['name']} ({format_size(f['size'])})")
            fastq_filenames = [f["name"] for f in files]
            azenta_data_config = {
                "data_source": "gdrive",
                "gdrive_url": gdrive_url,
                "fastq_dir": fastq_dir,
            }
            if st.button("Reset (re-download)", key="scrna_gdrive_reset"):
                st.session_state.pop("scrna_gdrive_downloaded", None)
                st.session_state.pop("scrna_gdrive_job_id", None)
                st.rerun()

        elif gdrive_url:
            if "drive.google.com" not in gdrive_url:
                st.error("This is not a Google Drive URL.")
            else:
                # === Phase 1: Download ===
                dl_job_id_gd = st.session_state.get("scrna_gdrive_job_id")

                if dl_job_id_gd:
                    job_dir = os.path.join(JOBS_DIR, dl_job_id_gd)
                    status_path = os.path.join(job_dir, "status.json")
                    if os.path.exists(status_path):
                        with open(status_path) as _f:
                            dl_status = json.load(_f)
                        dl_st = dl_status.get("status", "unknown")
                        if dl_st == "download_complete":
                            st.session_state.scrna_gdrive_downloaded = {
                                "fastq_dir": dl_status["fastq_dir"],
                                "files": dl_status["extracted_files"],
                            }
                            st.rerun()
                        elif dl_st == "error":
                            st.error(f"Download failed: {dl_status.get('error')}")
                            if st.button("Retry", key="scrna_gdrive_retry"):
                                st.session_state.pop("scrna_gdrive_job_id", None)
                                st.rerun()
                        else:
                            pct = dl_status.get("progress_pct", 0)
                            step = dl_status.get("step", "Downloading...")
                            st.progress(pct / 100, text=step)
                            st.info("Downloading... Refresh the page to check progress.")

                            @st.fragment(run_every=5)
                            def _gdrive_dl_auto_refresh():
                                _sp = os.path.join(JOBS_DIR, dl_job_id_gd, "status.json")
                                if os.path.exists(_sp):
                                    with open(_sp) as __f:
                                        _s = json.load(__f)
                                    if _s.get("status") in ("download_complete", "error"):
                                        st.rerun()
                            _gdrive_dl_auto_refresh()
                    else:
                        st.warning("Job status not found. Waiting...")
                else:
                    gdrive_work_dir = st.text_input(
                        "Work directory (download destination)",
                        value=str(Path.home() / "analysis"),
                        key="gdrive_work_dir",
                        placeholder="/data/analysis/project_scRNA",
                    )
                    if st.button("Download", type="primary", key="scrna_gdrive_dl"):
                        if not gdrive_work_dir:
                            st.error("Please specify a work directory.")
                        else:
                            dl_config = {
                                "data_source": "gdrive",
                                "gdrive_url": gdrive_url,
                                "work_dir": gdrive_work_dir,
                                "phase": "download_only",
                            }
                            job_id = submit_job(dl_config)
                            st.session_state.scrna_gdrive_job_id = job_id
                            st.rerun()

    elif azenta_source == "Local Folder":
        st.caption("Specify a local FASTQ directory directly (no download required).")

        local_folder_path = st.text_input(
            "FASTQ directory path",
            value=st.session_state.get("scrna_local_folder_path", ""),
            key="scrna_local_folder_path",
            placeholder="/data/sequencing/project/run_dir",
        )

        if local_folder_path and os.path.isdir(local_folder_path):
            # Collect FASTQ files via recursive search
            _fq_found = {}  # subdir -> [filenames]
            for root, dirs, files in os.walk(local_folder_path, followlinks=True):
                fqs = sorted([f for f in files if f.endswith((".fastq.gz", ".fq.gz"))])
                if fqs:
                    _fq_found[root] = fqs

            if _fq_found:
                _all_dirs = sorted(_fq_found.keys())
                _total_fq = sum(len(v) for v in _fq_found.values())
                # Collect unique filenames from all directories (used for GEM well detection etc.)
                _all_fq_names = set()
                for _fqs in _fq_found.values():
                    _all_fq_names.update(_fqs)
                fastq_filenames = sorted(_all_fq_names)

                if len(_all_dirs) == 1:
                    st.success(f"**{_total_fq} FASTQ files** in `{_all_dirs[0]}`")
                else:
                    st.success(f"**{_total_fq} FASTQ files** in {len(_all_dirs)} subdirectories (all used)")
                    for d in _all_dirs:
                        rel = os.path.relpath(d, local_folder_path)
                        st.caption(f"  `{rel}/` — {len(_fq_found[d])} files")

                with st.expander("File list", expanded=False):
                    for d in _all_dirs:
                        if len(_all_dirs) > 1:
                            rel = os.path.relpath(d, local_folder_path)
                            st.markdown(f"**{rel}/**")
                        for f in _fq_found[d]:
                            sz = os.path.getsize(os.path.join(d, f))
                            st.text(f"  {f} ({format_size(sz)})")

                # cellranger recursively searches --fastqs directories, so
                # exclude nested directories and pass only top-level ones
                _top_dirs = []
                for d in _all_dirs:
                    if not any(d.startswith(other + os.sep) for other in _all_dirs if other != d):
                        _top_dirs.append(d)
                azenta_data_config = {
                    "data_source": "local_folder",
                    "fastq_dirs": _top_dirs,
                }
            else:
                st.warning("No FASTQ files found under the specified directory.")
        elif local_folder_path:
            st.error(f"Directory does not exist: `{local_folder_path}`")

    # --- FASTQ Rename ---
    st.markdown("---")
    st.subheader("2. FASTQ Rename → Cell Ranger format")

    if fastq_filenames:
        # Check if already in Cell Ranger format
        already_cr = is_cellranger_format(fastq_filenames)

        if already_cr:
            st.success("FASTQ files are already in Cell Ranger format. No renaming needed.")
            st.session_state.scrna_rename_mapping = None
        else:
            # Auto-detect
            auto_mapping, auto_pattern_name = detect_rename_mapping(fastq_filenames)

            method_options = ["Auto-detect", "LLM (auto-infer)", "Custom regex", "No rename needed"]
            default_idx = 0

            rename_method = st.radio(
                "Rename method",
                method_options,
                index=default_idx,
                key="azenta_rename_method",
                horizontal=True,
            )

            rename_mapping = None

            if rename_method == "Auto-detect":
                if auto_mapping:
                    st.success(f"Pattern detected: **{auto_pattern_name}**")
                    rename_mapping = auto_mapping
                else:
                    st.warning("Auto-detection failed. Please use LLM or Custom regex.")

            elif rename_method == "LLM (auto-infer)":
                st.caption(
                    "Cell Ranger format: `{SampleName}_S1_L001_R{1or2}_001.fastq.gz`"
                )
                llm_rename_key = "scrna_llm_rename_result"
                if st.button("Infer rename with LLM", type="primary", key="azenta_llm_rename_btn"):
                    with st.spinner("Analyzing filename pattern with LLM..."):
                        result = rename_with_llm(fastq_filenames)
                        if result:
                            mapping, pattern_desc = result
                            st.session_state[llm_rename_key] = (mapping, pattern_desc)

                if llm_rename_key in st.session_state:
                    mapping, pattern_desc = st.session_state[llm_rename_key]
                    if pattern_desc:
                        st.info(f"Detected: {pattern_desc}")
                    rename_mapping = mapping

            elif rename_method == "Custom regex":
                custom_pattern = st.text_input(
                    "Regex pattern",
                    value=r'^(.+?)_L?(\d+)_?(\d+)\.fq\.gz$',
                    key="azenta_custom_pattern",
                )
                custom_replacement = st.text_input(
                    "Replacement",
                    value=r'\1_S1_L001_R\3_001.fastq.gz',
                    key="azenta_custom_replacement",
                )
                if custom_pattern and custom_replacement:
                    rename_mapping = apply_custom_rename(fastq_filenames, custom_pattern, custom_replacement)

            # Rename preview
            if rename_mapping:
                st.markdown("**Rename preview:**")
                rename_df = pd.DataFrame([
                    {"Before": k, "After": v} for k, v in rename_mapping.items()
                ])
                st.dataframe(rename_df, use_container_width=True, hide_index=True)

                # Validation
                invalid = [v for v in rename_mapping.values() if not CELLRANGER_FASTQ_RE.match(v)]
                if invalid:
                    st.error(f"Some files do not match Cell Ranger format: {invalid[:3]}")
                else:
                    st.success(f"All {len(rename_mapping)} files will be converted to Cell Ranger format.")

                st.session_state.scrna_rename_mapping = rename_mapping
            elif rename_method == "No rename needed":
                st.session_state.scrna_rename_mapping = None
                st.info("Continuing without renaming.")
            else:
                st.session_state.pop("scrna_rename_mapping", None)
    else:
        st.info("First acquire FASTQ files.")

    # --- Mode selection ---
    st.markdown("---")

    dropbox_downloaded = (azenta_source == "Dropbox"
                          and "scrna_dropbox_downloaded" in st.session_state)
    gdrive_downloaded = (azenta_source == "Google Drive"
                         and "scrna_gdrive_downloaded" in st.session_state)

    if (azenta_source == "Dropbox" and not dropbox_downloaded) or \
       (azenta_source == "Google Drive" and not gdrive_downloaded):
        st.subheader("3. Settings")
        _src_label = "Dropbox" if azenta_source == "Dropbox" else "Google Drive"
        st.info(f"You can configure the mode and samples after downloading from {_src_label}.")
        az_ref, az_ref_name, az_create_bam, az_also_run_count, az_cores, az_mem, az_run_annotate = common_cellranger_settings("azenta")
        azenta_mode = "auto"
        az_count_sample = az_count_fastqs = ""
        az_hashtags = az_libraries = az_samples = []
    else:
        st.subheader("3. Cell Ranger Mode & Settings")
        azenta_mode = st.selectbox(
            "Mode",
            ["multi", "count"],
            key="azenta_mode",
        )
        az_ref, az_ref_name, az_create_bam, az_also_run_count, az_cores, az_mem, az_run_annotate = common_cellranger_settings("azenta")
        # Default FASTQ directory from Local Folder / Dropbox / Google Drive
        _fq_dir_default = ""
        if azenta_source == "Local Folder" and azenta_data_config.get("fastq_dirs"):
            _fq_dir_default = ",".join(azenta_data_config["fastq_dirs"])
        elif dropbox_downloaded:
            _fq_dir_default = st.session_state["scrna_dropbox_downloaded"]["fastq_dir"]
        elif gdrive_downloaded:
            _fq_dir_default = st.session_state["scrna_gdrive_downloaded"]["fastq_dir"]
        # Set directly in session_state (widget value is only effective on first render)
        if _fq_dir_default:
            if not st.session_state.get("azenta_gex_fastqs"):
                st.session_state["azenta_gex_fastqs"] = _fq_dir_default
            if not st.session_state.get("azenta_cmo_fastqs"):
                st.session_state["azenta_cmo_fastqs"] = _fq_dir_default

        if azenta_mode == "count":
            _count_sample_default = ""
            # Auto-infer sample names from FASTQ filenames
            _detect_files = []
            if dropbox_downloaded:
                _detect_files = [f["name"] for f in st.session_state["scrna_dropbox_downloaded"]["files"]]
            elif gdrive_downloaded:
                _detect_files = [f["name"] for f in st.session_state["scrna_gdrive_downloaded"]["files"]]
            elif fastq_filenames:
                _detect_files = fastq_filenames
            if _detect_files:
                _sample_names = sorted(set(
                    re.match(r'^(.+?)_S\d+_L\d{3}', fn).group(1)
                    for fn in _detect_files
                    if re.match(r'^(.+?)_S\d+_L\d{3}', fn)
                    and not re.match(r'^(.+?)_S\d+_L\d{3}', fn).group(1).lower().startswith("undetermined")
                ))
                if _sample_names:
                    _count_sample_default = ", ".join(_sample_names)
            az_count_sample, az_count_fastqs = count_mode_ui(
                "azenta", fastq_dir_default=_fq_dir_default,
                sample_default=_count_sample_default,
                fastq_filenames=fastq_filenames if fastq_filenames else _detect_files if _detect_files else None,
            )
        else:
            az_hashtags, az_libraries, az_samples = multi_mode_ui(
                "azenta", fastq_dir_default=_fq_dir_default,
                fastq_filenames=fastq_filenames if fastq_filenames else None,
            )

    # Output destination
    st.markdown("---")
    st.subheader("4. Output")

    # When Dropbox / Google Drive download is complete, default work_dir to the download destination
    _default_work_dir = str(Path.home() / "analysis")
    if dropbox_downloaded:
        _dl_fdir = st.session_state["scrna_dropbox_downloaded"]["fastq_dir"]
        _default_work_dir = os.path.dirname(_dl_fdir)
    elif gdrive_downloaded:
        _dl_fdir = st.session_state["scrna_gdrive_downloaded"]["fastq_dir"]
        _default_work_dir = os.path.dirname(_dl_fdir)

    az_work_dir = st.text_input(
        "Work directory",
        value=_default_work_dir,
        key="azenta_work_dir",
        placeholder="/data/analysis/project_scRNA",
    )

    az_scala_dest = scala_dir_browser("azenta")

    # Submit
    st.markdown("---")
    st.subheader("5. Submit")

    can_submit = bool(azenta_data_config)
    if azenta_source == "Dropbox" and not dropbox_downloaded:
        can_submit = False
    if azenta_source == "Google Drive" and not gdrive_downloaded:
        can_submit = False

    # --- multi.csv preview ---
    if can_submit and azenta_mode == "multi":
        _preview_fq_dirs = ""
        if azenta_source == "Local Folder" and azenta_data_config.get("fastq_dirs"):
            _preview_fq_dirs = ",".join(azenta_data_config["fastq_dirs"])
        elif dropbox_downloaded:
            _preview_fq_dirs = st.session_state["scrna_dropbox_downloaded"]["fastq_dir"]
        elif gdrive_downloaded:
            _preview_fq_dirs = st.session_state["scrna_gdrive_downloaded"]["fastq_dir"]

        # Branch based on multi_configs or single
        if az_libraries and isinstance(az_libraries[0], dict) and "id" in az_libraries[0]:
            _preview_configs = az_libraries  # multi_configs
        else:
            _preview_configs = [{
                "id": "run",
                "hashtags": az_hashtags,
                "libraries": az_libraries,
                "samples": az_samples,
            }]

        with st.expander("Preview: Cell Ranger multi config", expanded=True):
            for pc in _preview_configs:
                well_id = pc.get("id", "run")
                st.markdown(f"**`cellranger multi --id={well_id}`**")

                csv_lines = ["[gene-expression]", f"reference,{az_ref}"]
                pc_hashtags = pc.get("hashtags", az_hashtags)
                if pc_hashtags:
                    csv_lines.append(f"cmo-set,(auto-generated)")
                csv_lines += ["", "[libraries]", "fastq_id,fastqs,feature_types"]
                for lib in pc.get("libraries", []):
                    fqs = lib.get("fastqs", _preview_fq_dirs) or _preview_fq_dirs
                    fq_id = lib["fastq_id"]
                    for fq_dir in fqs.split(","):
                        fq_dir = fq_dir.strip()
                        if not fq_dir:
                            continue
                        # Only show directories that contain files matching the fastq_id
                        if os.path.isdir(fq_dir):
                            has_match = any(
                                f.startswith(fq_id + "_S") and f.endswith((".fastq.gz", ".fq.gz"))
                                for f in os.listdir(fq_dir)
                            )
                            if not has_match:
                                continue
                        csv_lines.append(f"{fq_id},{fq_dir},{lib['feature_types']}")
                pc_samples = pc.get("samples", [])
                if pc_samples:
                    csv_lines += ["", "[samples]", "sample_id,cmo_ids"]
                    for s in pc_samples:
                        cmo_str = s["cmo_ids"] if isinstance(s["cmo_ids"], str) else "|".join(s["cmo_ids"])
                        csv_lines.append(f"{s['sample_id']},{cmo_str}")

                st.code("\n".join(csv_lines), language="text")

    if not can_submit:
        if azenta_source == "Dropbox":
            st.warning("First download and extract from Dropbox.")
        elif azenta_source == "Google Drive":
            st.warning("First download from Google Drive.")
        else:
            st.warning("First acquire FASTQ files.")
    elif st.button("Submit Job", type="primary", key="azenta_submit"):
        config = {
            **azenta_data_config,
            "mode": azenta_mode,
            "work_dir": az_work_dir,
            "reference": az_ref,
            "reference_name": az_ref_name,
            "create_bam": az_create_bam,
            "also_run_count": az_also_run_count,
            "run_annotate": az_run_annotate,
            "localcores": az_cores,
            "localmem": az_mem,
            "scala_dest": az_scala_dest,
        }

        # Rename settings
        rename_mapping = st.session_state.get("scrna_rename_mapping")
        if rename_mapping:
            config["rename"] = {
                "enabled": True,
                "mapping": rename_mapping,
            }

        if azenta_mode == "count":
            # Support multiple samples via comma-separated values
            _samples = [s.strip() for s in az_count_sample.split(",") if s.strip()]
            if len(_samples) > 1:
                config["count"] = {
                    "count_configs": [
                        {"sample_name": s, "id": s, "fastqs_dir": az_count_fastqs}
                        for s in _samples
                    ],
                }
            else:
                config["count"] = {
                    "sample_name": _samples[0] if _samples else az_count_sample,
                    "fastqs_dir": az_count_fastqs,
                }
        elif azenta_mode == "multi":
            # Determine the return value type from multi_mode_ui:
            # - multi-GEM-well: az_libraries = list of multi_config dicts (each with "id" key)
            # - single GEM well: az_libraries = list of library dicts (each with "fastq_id" key)
            if az_libraries and isinstance(az_libraries[0], dict) and "id" in az_libraries[0]:
                # Multi-GEM-well: az_libraries is multi_configs
                config["multi"] = {"multi_configs": az_libraries}
            else:
                config["multi"] = {
                    "hashtags": az_hashtags,
                    "libraries": az_libraries,
                    "samples": az_samples,
                }

        # For Upload source, save files to work_dir
        if azenta_source == "Upload" and "scrna_upload_files" in st.session_state:
            upload_dir = os.path.join(az_work_dir, "fastq")
            os.makedirs(upload_dir, exist_ok=True)
            with st.spinner("Saving uploaded files..."):
                for f in st.session_state.scrna_upload_files:
                    local_path = os.path.join(upload_dir, f.name)
                    with open(local_path, "wb") as out:
                        out.write(f.getbuffer())
            config["upload_dir"] = upload_dir

        job_id = submit_job(config)
        st.success(f"Job submitted! ID: {job_id}")
        st.info("Processing continues even if you close the browser. Check progress in the Job Monitor tab.")


# =============================================================================
#  Tab 3: Job Monitor
# =============================================================================

with tab_monitor:
    st.subheader("Job Monitor")

    # Delete job directories and work_dirs older than 7 days
    _cleanup_days = 7
    _now = time.time()
    _old_dt = _now - _cleanup_days * 86400
    if os.path.isdir(JOBS_DIR):
        for _jname in os.listdir(JOBS_DIR):
            _jdir = os.path.join(JOBS_DIR, _jname)
            if not os.path.isdir(_jdir):
                continue
            if os.path.getmtime(_jdir) >= _old_dt:
                continue
            # Also delete work_dir
            _cp = os.path.join(_jdir, "config.json")
            if os.path.isfile(_cp):
                try:
                    with open(_cp) as _f:
                        _jc = json.load(_f)
                    _wd = _jc.get("work_dir", "")
                    if _wd and os.path.isdir(_wd) and _wd != "/":
                        shutil.rmtree(_wd, ignore_errors=True)
                except Exception:
                    pass
            shutil.rmtree(_jdir, ignore_errors=True)

    auto_refresh = st.checkbox("Auto-refresh (10 sec)", value=False, key="scrna_auto_refresh")

    if st.button("Refresh", key="scrna_refresh"):
        st.rerun()

    @st.fragment(run_every=10 if auto_refresh else None)
    def _scrna_monitor_fragment():
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
                    source = c.get("data_source", "?").upper()
                    mode = c.get("mode", "?")
                    label = f"**Job {job['name']}** [{source}]"
                    st.markdown(label)
                with col_status:
                    st.markdown(f":{color}[{status_str.upper()}]")
                with col_info:
                    ref_name = c.get("reference_name", "?")
                    st.markdown(f"{mode} / {ref_name}")

                progress = s.get("progress_pct", 0) / 100
                st.progress(progress, text=s.get("step", ""))

                with st.expander("Details"):
                    st.markdown(f"**Created:** {c.get('created_at', '?')}")
                    st.markdown(f"**Started:** {s.get('started_at', 'Not yet')}")
                    st.markdown(f"**Finished:** {s.get('finished_at', 'Not yet')}")
                    st.markdown(f"**Mode:** {c.get('mode', '?')}")
                    st.markdown(f"**Work dir:** `{c.get('work_dir', '?')}`")
                    if c.get("scala_dest"):
                        st.markdown(f"**SCALA dest:** `{c['scala_dest']}`")

                    if s.get("error"):
                        st.error(f"Error: {s['error']}")

                    # web_summary.html download
                    if status_str == "done":
                        summaries = s.get("web_summaries", [])
                        if summaries:
                            st.markdown("**Web Summaries:**")
                            for ws in summaries:
                                if os.path.isfile(ws):
                                    with open(ws, "rb") as fh:
                                        fname = os.path.basename(ws)
                                        st.download_button(
                                            f"Download {fname}",
                                            data=fh.read(),
                                            file_name=fname,
                                            mime="text/html",
                                            key=f"dl_{job['name']}_{fname}",
                                        )
                        # cell_types.csv download (cellranger annotate)
                        work_dir = c.get("work_dir", "")
                        if work_dir and c.get("run_annotate"):
                            import glob as glob_mod
                            ct_files = glob_mod.glob(os.path.join(work_dir, "*_annotate", "outs", "cell_types", "cell_types.csv"))
                            if ct_files:
                                st.markdown("**Cell Type Annotation:**")
                                for ct in ct_files:
                                    annotate_id = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(ct))))
                                    run_id = annotate_id.replace("_annotate", "")
                                    with open(ct, "rb") as fh:
                                        st.download_button(
                                            f"Download cell_types.csv ({run_id})",
                                            data=fh.read(),
                                            file_name=f"{run_id}_cell_types.csv",
                                            mime="text/csv",
                                            key=f"dl_{job['name']}_{run_id}_celltypes",
                                        )
                        if s.get("scala_copied"):
                            st.success(f"Results copied to SCALA: `{c.get('scala_dest', '')}`")

                    log_path = os.path.join(job["dir"], "pipeline.log")
                    if os.path.exists(log_path):
                        if st.button("Show log", key=f"scrna_log_{job['name']}"):
                            with open(log_path) as f:
                                lines = f.readlines()
                                tail = lines[-80:] if len(lines) > 80 else lines
                                st.code("".join(tail), language="text")

                st.markdown("---")

    _scrna_monitor_fragment()
