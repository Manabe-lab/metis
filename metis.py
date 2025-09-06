#!/usr/bin/env python3
"""
METIS - Molecular Exploration and Transcriptomic Investigation Suite

Copyright (C) 2024 METIS Development Team

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

"""

import streamlit as st

st.set_page_config(
    page_title="METIS",
    page_icon="🧝🏻‍♀️",
    layout="wide"
)

# Home page function
def home():
    # Display logo image
    import os
    logo_path = os.path.join(os.path.dirname(__file__), "Metis_favicon.png")
    if os.path.exists(logo_path):
        st.image(logo_path, width=150)
    
    st.markdown(
    """
    #### ***M***olecular ***E***xploration and ***T***ranscriptomic ***I***nvestigation ***S***uite

    #### RNA-seq workflow
    1. Count Data Normalization → PCA, Heatmap (QC and data distribution check)
    2. DESeq2 → Volcano plot
    3. Gene set analysis: DESeq2 → GSEApy, decoupleR
        - Alternatively, over-representation analysis using DEG lists (GSEApy, decoupleR)

    #### Usage Examples
    - __Creating data for Heatmap:__
    Count Data Normalization. Or extract rlog data using "Homer to Data" from DESeq2 results.
    Filter by expression levels or changes using "Filter Log Z-score". Basic filtering can also be done in Heatmap.

    - __Updating gene names from Homer etc:__
    Use "Update Gene Symbol". Since many gene names are updated, it's better to update depending on downstream analysis methods.

    - __Converting gene names to Ensembl ID:__
    Use "Homer to Ensembl". For example, DIANE uses Ensembl IDs, so convert using this tool.

    - __Data file organization:__
    Remove unnecessary columns/rows. Extract only needed columns. Change row/column names, reorder columns, transpose rows and columns using "Manipulate Data Table".


    #### Data Normalization/Manipulation
    - __Count Data Normalization:__
    RNA-seq count data normalization, QC PCA, heatmap, box plot and other preprocessing 

    - __Homer or DEseq2 to Data:__
    Convert Homer output files to data-only files
    Extract only rlog data from DESeq2

    - __Count to TPM:__
    Convert raw count data to TPM

    - __Manipulate Data Table:__
    Delete/extract rows/columns from data table. Change column order, row names, column names

    - __Merge Data Table:__
    Merge data files

    - __Filter and Transform Data:__
    Data filtering. Data XY transposition.

    - __Spread Sheet:__
    Spreadsheet editor

    #### Gene name conversion
    - __Update Gene Symbol:__
    Update gene symbols in Homer output or general data files (first column is Symbol)

    - __Gene Symbol to EnsemblID:__
    Convert gene symbols to Ensembl IDs in Homer output or general data files (for DIANE file creation)

    - __EnsemblID to Gene Symbol__
    Handle duplicate genes as well

    - __Mouse Human Symbol Conversion__

    - __GMT Mouse Human Conversion:__
    Convert gene symbols in gene set files (gmt files) used in GSEA. Format validation for gmt files is also available.


    #### Data Analysis
    - __DESeq2:__
    DESeq2, limma-eBayes, Beta regression, GLM group comparisons

    - __DESeq2-LRT etc:__
    ANOVA-like, time series analysis etc. limma-eBayes for non-count data analysis. Beta regression for 0-1 data

    - __edgeR:__
    edgeR

    - __limma:__
    limma voom calculation

    - __DE method comparison:__
    Compare DEG calculation methods (DESeq2, edgeR, limma) results

    - __Permutation test:__
    Freedman-Lane test including batch effects

    - __Make ranking file for GSEA:__
    Create GSEA ranking files

    - __Batch removal by Combat-seq:__
    Batch removal using Combat-seq

    - __impulseDE2:__
    Time course analysis using impulseDE2

    - __Mfuzz clustering:__
    Time course data clustering using Mfuzz

    - __DEA result comparison:__
    Compare significantly different genes from DEA results

    - __Compare ratios:__
    Statistical analysis of ratio data using t-test and β regression on logit transformed data

    #### Data Visualization
    - __PCA:__
    PCA, UMAP, tSNE, MDS

    - __Volcano plot__

    - __Heatmap__

    - __Box/Violin plot__

    - __Venn and Upset Plot__ Venn diagram and UpSet plot

    #### Pathway Analysis
    - __decoupleR:__
    Signal pathway, TF activity analysis, and GSEA using mSigDB. Can create publication quality enrichment plots

    - __GSEApy:__
    GSEA

    - __PPI analysis:__
    Protein-protein interaction analysis using STRING-db

    #### WGCNA
    - __WGCNA:__
    WGCNA analysis using WGCNApy or R WGCNA

    - __WGCNA network plot:__
    Network visualization of WGCNA modules

    - __WGCNA hub UMAP:__
    UMAP visualization of WGCNA module genes

    - __WGCNA objects comparison:__
    Display relationships between WGCNA modules

    - __Generate gmt from cluster info:__
    Create gmt files from cluster results or WGCNA module genes

    #### scRNA-seq
    - __Pseudobulk:__
    Create pseudobulk data from anndata (h5ad)

    - __Metacells by SEACells:__
    Create metacells using SEACells

    - __Random pseudo-replicates:__
    Create pseudo-replicates by random cell splitting

    - __SCENIC heatmap:__
    Heatmap visualization of SCENIC gene regulatory networks

    - __SCENIC CSI:__
    Calculate connection specificity index (CSI) of SCENIC regulons

    - __SCENIC network analysis:__
    Network visualization of SCENIC regulons centered on transcription factors or their targets

    - __SCENIC multinetwork analysis:__
    Network visualization of SCENIC regulons centered on multiple transcription factors and targets

    #### Cell communication
    - __LIANA LR analysis:__
    Ligand-receptor analysis using LIANA+

    - __CellChat:__
    CellChat analysis from h5ad files - faithful Python implementation of R version

    - __CellChat comparison:__
    Two-condition comparison from h5ad files containing both conditions

    - __CellChat permutation test:__
    Statistical testing of two-condition comparison using permutation test

    - __CellChat R qs to python:__
    Convert SCALA CellChat analysis result qs files for use in metis

    #### scRNA file operation
    - __Download public data for SCALA and cellxgene:__
    Download files analyzable in SCALA/cellxgene. Can download from GEO and other sources if download links are available

    - __File explorer:__
    Browse, delete, and download files accessible by SCALA

    - __File uploader:__
    Upload files accessible by SCALA/cellxgene

    - __SCALA file browser:__
    Browse files accessible by SCALA

    #### ChIP-seq
    - __Sort BAM file:__
    Sort BAM files in METIS_data directory

    - __Merge BAM files:__
    Merge BAM files in METIS_data directory

    - __Convert bam to bedGraph for SEACR:__
    Convert BAM files to bedGraph format for SEACR

    - __Peak calling with SEACR:__

    - __Peak calling with macs3:__

    - __Annotating and filtering peaks:__
    Annotation and annotation-based filtering of peak files and bed files

    - __Bed file filter for length score:__
    Filtering by peak length and score. Peak length distribution is also displayed

    - __Bam to DESeq2:__
    Normalize CUT&RUN peak counts based on greenlist counts for DESeq2 analysis

    - __Filter out blacklist:__
    Filter peaks overlapping with blacklist regions

    - __Convert bed to fasta:__
    Can also handle MACS peak files

    """
    )

# Navigation setup - maintaining section divisions
pg = st.navigation({
    "Home": [
        st.Page(home, title="Home", icon="🏛"),
    ],
    "Data Normalization/Manipulation 📝": [
        st.Page("pages/normalization.py", title="Count Data Normalization"),
        st.Page("pages/count2tpm.py", title="Count to TPM"),
        st.Page("pages/1_Homer to Data.py", title="Homer or DESeq2 to Data"),
        st.Page("pages/4_Manipulate Data Table.py", title="Manipulate Data Table"),
        st.Page("pages/merge_data.py", title="Merge Data Files"),
        st.Page("pages/5_Filter Log Z-score.py", title="Filter and Transform Data"),
        st.Page("pages/SpreadSheet.py", title="Spread Sheet"),
    ],
    "Gene name conversion 🔁": [
        st.Page("pages/2_Update Gene Symbol.py", title="Update Gene Symbol"),
        st.Page("pages/3_Homer to Ensembl.py", title="Gene Symbol to Ensembl ID"),
        st.Page("pages/Ensembl2Symbol.py", title="Ensembl ID to Gene Symbol"),
        st.Page("pages/mouse-human.py", title="Mouse Human Symbol Conversion"),
        st.Page("pages/gmt-mouse-human.py", title="GMT Mouse Human Conversion"),
    ],
    "Data Analysis 🧮": [
        st.Page("pages/6_Calc_DESeq2.py", title="DESeq2 etc"),
        st.Page("pages/Calc_DESeq2_LRT.py", title="DESeq2-LRT etc"),
        st.Page("pages/edgeR.py", title="edgeR"),
        st.Page("pages/limma.py", title="limma"),
        st.Page("pages/DE_rpy2.py", title="DE method comparison"),
        st.Page("pages/permutation_test.py", title="Permutation test"),
        st.Page("pages/8_rnkgene.py", title="Make ranking file for GSEA"),
        st.Page("pages/combat-seq.py", title="Batch removal by Combat-seq"),
        st.Page("pages/impulsede2-streamlit-app.py", title="impulseDE2"),
        st.Page("pages/rnaseq-mfuzz-streamlit-app.py", title="Mfuzz clustering"),
        st.Page("pages/CompareDE.py", title="DEA result comparison"),
        st.Page("pages/analyze_freq.py", title="Compare ratios"),
    ],
    "Data Visualization 🌋": [
        st.Page("pages/pca.py", title="PCA"),
        st.Page("pages/7_Volcano Plot.py", title="Volcano plot"),
        st.Page("pages/Heatmap.py", title="Heatmap"),
        st.Page("pages/Boxplot.py", title="Box_Violin plot"),
        st.Page("pages/venn_upset.py", title="Venn_Upset Plot"),
    ],
    "Pathway Analysis 🔀": [
        st.Page("pages/decoupler.py", title="decoupleR"),
        st.Page("pages/GSEApy.py", title="GSEApy"),
        st.Page("pages/StringDB.py", title="PPI analysis"),
    ],
    "WGCNA 🥅": [
        st.Page("pages/PyWGCNA.py", title="WGCNA"),
        st.Page("pages/WGCNA_plot.py", title="WGCNA network plot"),
        st.Page("pages/WGCNAumap.py", title="WGCNA hub UMAP"),
        st.Page("pages/PyWGCNA_comparison.py", title="WGCNA objects comparison"),
        st.Page("pages/cluster2gmt.py", title="Generate gmt from cluster info"),
    ],
    "scRNA-seq 🎡": [
        st.Page("pages/pseudobulk.py", title="Pseudobulk"),
        st.Page("pages/SEAcells.py", title="Metacells by SEACells"),
        st.Page("pages/pseudoreplicates.py", title="Random pseudo-replicates"),
        st.Page("pages/Memento.py", title="memento DE analysis"),
        st.Page("pages/Memento2D.py", title="memento 2D analysis"),
        st.Page("pages/SCENICviewer.py", title="SCENIC heatmap"),
        st.Page("pages/prepare_scenic_data.py", title="Prepare regulon data for heatmap"),
        st.Page("pages/SCENICcsi.py", title="SCENIC CSI"),
        st.Page("pages/SCENIC_network.py", title="SCENIC network analysis"),
        st.Page("pages/SCENIC_multinetwork.py", title="SCENIC multinetwork analysis"),
    ],
    "Cell communication 💬": [
        st.Page("pages/liana_steady.py", title="LIANA LR analysis"),
        st.Page("pages/liana_comparison.py", title="LIANA comparison"),
        st.Page("pages/cellchat.py", title="CellChat"),
        st.Page("pages/cellchat_comparison.py", title="CellChat comparison"),
        st.Page("pages/cellchat_permutation.py", title="CellChat permutation test"),
        st.Page("pages/cellchatR2py.py", title="CellChat R qs to python"),
    ],
    "ChIP-seq 🧬": [
        st.Page("pages/bam-sorter.py", title="Sort BAM file"),
        st.Page("pages/bam-merge.py", title="Merge BAM files"),
        st.Page("pages/bam2bedgraph.py", title="Bam to bedGraph for SEACR"),
        st.Page("pages/seacr-streamlit-app.py", title="SEACR peak calling"),
        st.Page("pages/macs3-peak-calling-app.py", title="Macs3 peak calling"),
        st.Page("pages/bed-peak-filter.py", title="Annotate_Filter peaks"),
        st.Page("pages/peak-filter.py", title="Bed length_score filter"),
        st.Page("pages/bam2DESeq2.py", title="Bam to DESeq2"),
        st.Page("pages/filter_blacklist.py", title="Blacklist filter"),
        st.Page("pages/bed2fa.py", title="Bed to fasta"),
        st.Page("pages/cut-and-run-denoise.py", title="Denoise bedgraph bigwig"),
    ],
    "MISC Ⓜ": [
        st.Page("pages/merge_excel.py", title="Merge excel files"),
        st.Page("pages/union.py", title="Remove duplicates"),
        st.Page("pages/SplitonKey.py", title="Split data file on key in another file"),
    ],
})

# Run navigation
pg.run()