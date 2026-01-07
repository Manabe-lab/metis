# METIS GitHub Repository Update Log - 2026-01-07

## Summary
Updated metis-github repository to sync with local streamlit version, removed unused files, and translated Japanese UI text to English.

---

## 1. Removed Unused Files (29 files)

Files not referenced by `metis.py` were removed from the repository:

### CellChat Variants
- cellchat_CPU_original.py
- cellchat_gpu.py
- cellchat_GPU.py
- cellchat_optimized.py
- cellchat_vis.py

### RNA Velocity (old versions)
- velocity_analysis.py (replaced with new version from local)
- velocity_filter.py (replaced with new version from local)
- velocity_visualization.py (replaced with new version from local)
- veloviz.py (replaced with new version from local)
- TFvelo.py (replaced with new version from local)

### Dynamo (old versions)
- dynamo_analysis.py (replaced with new version from local)
- dynamo_visualization.py (replaced with new version from local)
- dynamo_perturbation.py (replaced with dynamo_perturbation_v2.py)
- dynamo_perturbation_visualization.py
- dynamo_lap_correct.py (replaced with new version from local)

### CellRank (old versions)
- cellrank_analysis.py (replaced with new version from local)
- cellrank_visualization.py (replaced with new version from local)
- deepvelo_analysis.py (replaced with new version from local)

### Others
- add_metadata.py (replaced with new version from local)
- celltypist.py (replaced with new version from local)
- fip.py
- KStest.py (replaced with new version from local)
- mouse-human-updated.py
- paperqa2_cleaned.py
- pca_statistics.py (replaced with new version from local)
- pseudotime_gene_expression.py (replaced with new version from local)
- sccoda_analysis.py (replaced with new version from local)
- SCENIC_network_multi.py
- selenium_helpers.py

---

## 2. Updated metis.py Structure

Synced with local version, adding new sections and pages:

### New Sections Added
- **SCENIC** (separate section from scRNA-seq)
- **RNA velocity** (new dedicated section)
- **scRNA file operation** (expanded)

### New Pages Added (25 files from local)
- add_metadata.py
- cellrank_analysis.py
- cellrank_visualization.py
- celltypist.py
- data_file_browser.py
- deepvelo_analysis.py
- download.py
- dynamo_analysis.py
- dynamo_lap_correct.py
- dynamo_perturbation_v2.py
- dynamo_visualization.py
- filebrowser.py
- fileexplorer.py
- filter_xy_genes.py
- ftp.py
- KStest.py
- microarray_gene_filter.py
- pca_statistics.py
- pseudotime_gene_expression.py
- sccoda_analysis.py
- TFvelo.py
- velocity_analysis.py
- velocity_filter.py
- velocity_visualization.py
- veloviz.py

### Excluded from GitHub (local only)
- grants.py
- paperqa2.py
- tts_generator.py
- form408_to_ms_json.py

---

## 3. Japanese to English Translation

Translated Japanese UI text to English across 66+ files.

### Major Files Translated
| File | Japanese Lines (Before) | Japanese Lines (After) |
|------|------------------------|----------------------|
| pca_statistics.py | 1105 | 0 |
| cellchat.py | 889 | 0 |
| Calc_DESeq2_LRT.py | 320 | ~64 |
| SCENIC_multinetwork.py | 281 | ~48 |
| KStest.py | 209 | ~43 |
| cellchat_permutation.py | 207 | ~40 |
| sccoda_analysis.py | 183 | ~21 |
| And 60+ more files... | | |

### Translation Coverage
- **Before**: ~2500+ lines with Japanese text
- **After**: ~126 lines with some remaining kanji characters
- **Coverage**: ~95% translated

### Remaining Kanji
Some single kanji characters remain that require context for proper translation. These are mostly in specialized scientific/statistical terminology.

---

## 4. Git Commits

### Commit 1: Remove unused files
```
Remove 29 unused files not referenced by metis.py
```

### Commit 2: Sync with local
```
Sync metis.py with local version and add missing page files
```

### Commit 3: Translation
```
Translate Japanese UI text to English in all page files
```

---

## Repository Status After Update

- **Total pages**: 80+ Python files
- **Language**: Primarily English UI
- **Structure**: Organized by analysis category (Data Normalization, Gene conversion, Data Analysis, etc.)
- **New features**: RNA velocity section, SCENIC section, expanded scRNA tools

---

## Notes

- Debug files with Japanese names were removed earlier in the session
- `.gitignore` was updated to exclude `__pycache__` directories
- All changes pushed to: https://github.com/Manabe-lab/metis
