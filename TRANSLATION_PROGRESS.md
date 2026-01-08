# Translation Progress - 2026-01-07

## Current Status: PAUSED

The previous character-by-character translation created garbled text. The approach has been changed to:
1. Copy original files from local (/home/ichiro/streamlit/pages/) that have proper Japanese
2. Perform proper translation

---

## Completed Steps

### 1. Copied files from local (with proper Japanese)
The following 36 files were copied from local to replace the garbled versions:

**First batch (top garbled files):**
- cellchat.py
- dynamo_visualization.py
- dynamo_analysis.py
- sccoda_analysis.py
- dynamo_perturbation_v2.py
- dynamo_lap_correct.py
- pseudotime_gene_expression.py
- data_file_browser.py
- CompareDE.py
- microarray_gene_filter.py
- cellrank_visualization.py
- cellchatR2py.py
- velocity_visualization.py
- TFvelo.py
- KStest.py
- Calc_DESeq2_LRT.py
- add_metadata.py
- velocity_filter.py
- cellrank_analysis.py
- velocity_analysis.py

**Second batch:**
- 7_Volcano Plot.py
- bed2fa.py
- cellchat_permutation.py
- DE_rpy2.py
- liana_comparison.py
- liana_steady.py
- Memento2D.py
- Memento.py
- pca_statistics.py
- PyWGCNA_comparison.py
- PyWGCNA.py
- SCENIC_multinetwork.py
- SEAcells.py
- veloviz.py
- venn_upset.py
- WGCNA_plot.py

---

## Remaining Work

### TODO: Translate Japanese to English
All 36 files now have proper Japanese text and need to be translated to English.

**Translation approach should be:**
1. Use context-aware translation (not character-by-character)
2. Preserve Python code structure
3. Only translate strings within quotes (help text, labels, comments, docstrings)
4. Test files after translation to ensure no syntax errors

### Files with most Japanese text (priority order):
1. cellchat.py - ~1000+ lines
2. pca_statistics.py - ~1100+ lines
3. Calc_DESeq2_LRT.py - ~300+ lines
4. SCENIC_multinetwork.py - ~280 lines
5. KStest.py - ~200 lines
6. cellchat_permutation.py - ~200 lines
7. Others...

---

## Commands to Resume

```bash
# Check current git status
cd /home/ichiro/metis-github/pages
git status

# Check which files have Japanese text
grep -l '[ぁ-んァ-ン一-龥]' *.py | wc -l

# Count Japanese lines per file
for f in *.py; do count=$(grep -c '[ぁ-んァ-ン一-龥]' "$f" 2>/dev/null || echo 0); if [ "$count" -gt 0 ]; then echo "$count $f"; fi; done | sort -rn

# After translation, commit:
git add -A
git commit -m "Translate Japanese UI text to English (proper translation)"
git push
```

---

## Note
Do NOT use character-by-character sed translation. Use proper context-aware translation for Japanese → English.
