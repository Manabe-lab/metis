#!/usr/bin/env python3
"""
METIS GitHub Upload Processing Script

This script processes Python files for GitHub upload by:
1. Translating Japanese comments and UI text to English
2. Converting absolute paths to relative paths
3. Excluding MISC files from processing

Usage:
    python process_files_for_upload.py

Workflow:
1. Run this script from /home/ichiro/metis-github/
2. Review changes in the pages/ directory
3. Commit and push to GitHub:
   git add .
   git commit -m "Update files with translations"
   git push origin main
"""

import os
import re
import shutil
from datetime import datetime

# Configuration
SOURCE_DIR = "/home/ichiro/streamlit/pages"
DEST_DIR = "/home/ichiro/metis-github/pages"

# MISC files - excluded from upload
MISC_FILES = [
    "merge_excel.py",
    "union.py",
    "SplitonKey.py",
    "paperqa2.py",
    "grants.py",
    "tts_generator.py",
    "form408_to_ms_json.py"
]

# Files to exclude (backups, copies, etc.)
EXCLUDE_PATTERNS = [
    "copy",
    "backup",
    ".bak",
    ".org",
    "org_this_works",
    "batch_becomes_smaller",
    "_old",
    "test_",
]

# Path replacements (absolute to relative)
PATH_REPLACEMENTS = [
    ("/home/cellxgene/streamlit/", ""),
    ("/home/ichiro/streamlit/", ""),
    ('"/home/cellxgene/streamlit/', '"'),
    ('"/home/ichiro/streamlit/', '"'),
]

def should_exclude(filename):
    """Check if file should be excluded from processing"""
    if filename in MISC_FILES:
        return True
    filename_lower = filename.lower()
    for pattern in EXCLUDE_PATTERNS:
        if pattern.lower() in filename_lower:
            return True
    return False

def convert_paths(content):
    """Convert absolute paths to relative paths"""
    for old_path, new_path in PATH_REPLACEMENTS:
        content = content.replace(old_path, new_path)
    return content

def get_files_to_process():
    """Get list of files that need processing (newer in source or new)"""
    files_to_process = []

    for filename in os.listdir(SOURCE_DIR):
        if not filename.endswith('.py'):
            continue
        if should_exclude(filename):
            continue

        src_path = os.path.join(SOURCE_DIR, filename)
        dst_path = os.path.join(DEST_DIR, filename)

        src_mtime = os.path.getmtime(src_path)

        if os.path.exists(dst_path):
            dst_mtime = os.path.getmtime(dst_path)
            if src_mtime > dst_mtime:
                files_to_process.append((filename, 'updated'))
        else:
            files_to_process.append((filename, 'new'))

    return sorted(files_to_process)

def process_file(filename, translate_func=None):
    """Process a single file"""
    src_path = os.path.join(SOURCE_DIR, filename)
    dst_path = os.path.join(DEST_DIR, filename)

    with open(src_path, 'r', encoding='utf-8') as f:
        content = f.read()

    # Convert absolute paths to relative
    content = convert_paths(content)

    # Apply translation if provided
    if translate_func:
        content = translate_func(content)

    # Write to destination
    with open(dst_path, 'w', encoding='utf-8') as f:
        f.write(content)

    return True

if __name__ == "__main__":
    print("METIS GitHub Upload Processor")
    print("=" * 50)

    files = get_files_to_process()
    print(f"\nFiles to process: {len(files)}")

    for filename, status in files:
        print(f"  [{status}] {filename}")

    print("\nNote: Japanese to English translation must be done manually")
    print("or by calling this script with a translation function.")
