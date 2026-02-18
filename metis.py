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
import streamlit.components.v1 as components

st.set_page_config(
    page_title="METIS",
    page_icon="Metis_favicon.png",
    layout="wide"
)

# Close all sections while keeping Streamlit internal state in sync
components.html("""
<script>
function closeAllSectionsProper() {
    const doc = window.parent.document;
    const headers = doc.querySelectorAll('[data-testid="stNavSectionHeader"]');

    let closedCount = 0;

    headers.forEach((header, idx) => {
        // Check if next element is visible
        const nextSibling = header.nextElementSibling;
        const isExpanded = nextSibling &&
                          (nextSibling.style.display !== 'none' &&
                           nextSibling.offsetHeight > 0);

        if (isExpanded) {
            console.log(`Closing section ${idx} via click`);
            // Properly call Streamlit handler
            header.click();
            closedCount++;
        }
    });

    console.log(`Closed ${closedCount} sections`);
    return closedCount;
}

// Execute only on initial load
if (!sessionStorage.getItem('sections-properly-closed')) {
    window.addEventListener('load', () => {
        // Wait for Streamlit initialization
        setTimeout(() => {
            const closed = closeAllSectionsProper();
            if (closed > 0) {
                sessionStorage.setItem('sections-properly-closed', 'true');
            }
        }, 1000);

        // Re-execute just in case
        setTimeout(() => {
            closeAllSectionsProper();
        }, 2000);
    });
}
</script>
""", height=0)


# Home page function
def home():
    # Display logo image
    import os
    logo_path = os.path.join(os.path.dirname(__file__), "Metis_logo.png")
    if os.path.exists(logo_path):
        st.image(logo_path, width=150)

    st.markdown(
    """
    #### ***M***olecular ***E***xploration and ***T***ranscriptomic ***I***nvestigation ***S***uite

    #### Data Normalization/Manipulation
    - __Count Data Normalization:__
    RNA-seq count data normalization, QC PCA, heatmap, box plot and other preprocessing

    - __Homer or DEseq2 to Data:__
    Convert Homer output files to data-only files.
    Extract only rlog data from DESeq2.

    - __Count to TPM:__
    Convert raw count data to TPM

    - __Manipulate Data Table:__
    Delete/extract rows and columns from data table. Reorder columns, rename rows and columns.

    - __Merge Data Table:__
    Merge data files

    - __Filter and Transform Data:__
    Data filtering. Data XY transposition.

    - __Spread Sheet:__
    Spreadsheet editor

    - __Batch removal by Combat-seq:__
    Batch removal using Combat-seq

    - __Filter X/Y Chromosome Genes:__
    Filter genes on X/Y chromosomes (remove/extract/classify). Supports 10x Genomics and HOMER annotation.
    Multiple delimiters (newline, tab, comma, space) for gene list input.

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

    - __Microarray Gene Name Filter:__
    Extract and filter gene information from microarray data. Select gene name column to place in first column. Supports duplicate gene aggregation.


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

    - __impulseDE2:__
    Time course analysis using impulseDE2

    - __Mfuzz clustering:__
    Time course data clustering using Mfuzz

    - __DEA result comparison:__
    Compare significantly different genes from DEA results

    - __Compare ratios:__
    Statistical analysis of ratio data using t-test and beta regression on logit transformed data

    - __KS Test Distribution Analysis:__
    Two-group distribution comparison using Kolmogorov-Smirnov test and Anderson-Darling test

    #### Data Visualization
    - __PCA:__
    PCA, UMAP, tSNE, MDS

    - __PCA statistics:__
    Statistical analysis of PCA principal component scores and covariates. OLS (standard/robust SE), LMM (linear mixed model), Freedman-Lane permutation test for multiple comparison adjusted p-values. Estimated marginal means (EMM) calculation supported.

    - __Volcano plot__

    - __Heatmap__

    - __Box/Violin plot__

    - __Venn and Upset Plot__ Venn diagram and UpSet plot

    #### Pathway Analysis
    - __decoupleR:__
    Signal pathway, TF activity analysis, and GSEA using mSigDB. Can create publication quality enrichment plots

    - __GSEApy:__
    GSEA

    - __aPEAR Enrichment Network:__
    Network visualization of pathway enrichment results. Auto-detection from GSEApy, decoupleR, PADOG output TSV. Clustering and labeling based on gene overlap between pathways

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
    - __Dimension reduction:__
    Dimension reduction for scRNA-seq data (h5ad) (UMAP, t-SNE, Diffusion Map, PHATE, PaCMAP, TriMap, etc.)

    - __Pseudobulk:__
    Create pseudobulk data from anndata (h5ad)

    - __Metacells by SEACells:__
    Create metacells using SEACells

    - __Random pseudo-replicates:__
    Create pseudo-replicates by random cell splitting

    - __Pseudobulk DEA (RUV):__
    Pseudobulk from h5ad -> RUV correction (RUVg/RUV2/RUVIII) -> Differential expression analysis with edgeR/DESeq2

    - __memento DE analysis:__
    Differential expression analysis using memento

    - __memento 2D analysis:__
    memento 2D analysis

    - __scCODA compositional analysis:__
    Differential analysis of cell type composition using scCODA hierarchical Bayesian model for synthetic data

    - __Add metadata to h5ad:__
    Add metadata from another h5ad or TSV file to h5ad file. Match by cell barcode

    - __CellTypist annotation:__
    Automatic cell type annotation using machine learning. Supports 60 pre-trained models (Human/Mouse tissues)

    - __COMPASS metabolic analysis:__
    Single-cell metabolic flux analysis. Estimates metabolic reaction activity from gene expression using Flux Balance Analysis (FBA)

    - __COMPASS postprocess:__
    Differential analysis of COMPASS results (Wilcoxon test, Cohen's d, volcano plot, subsystem aggregation)

    - __scFEA metabolic analysis:__
    Fast metabolic flux estimation using Graph Neural Network. Significantly faster than COMPASS with no license required

    - __scFEA postprocess:__
    Differential analysis of scFEA results (Wilcoxon test, Cohen's d, volcano plot, supermodule aggregation). Supports both Flux/Balance modes

    - __Geneformer perturbation:__
    In silico perturbation analysis using Transformer foundation model. Predicts effects of gene knockout/overexpression

    #### SCENIC
    - __SCENIC heatmap:__
    Heatmap visualization of SCENIC gene regulatory networks

    - __Prepare regulon data for heatmap:__
    Preprocess SCENIC regulon data for heatmap

    - __SCENIC CSI:__
    Calculate connection specificity index (CSI) of SCENIC regulons

    - __SCENIC network analysis:__
    Network visualization of SCENIC regulons centered on transcription factors or their targets

    - __SCENIC multinetwork analysis:__
    Network visualization of SCENIC regulons centered on multiple transcription factors and targets

    - __SCENIC network analysis 2:__
    SCENIC network analysis with Community Detection, Feed-Forward Loop detection, AUCell differential analysis, Feedback Loop detection, and TF family annotation

    #### RNA velocity
    - __Data filtering:__
    Filter and merge loom files generated by velocyto based on cells in h5ad file

    - __scVelo analysis:__
    RNA velocity analysis using scVelo

    - __scVelo visualization:__
    Visualization of scVelo analysis results (velocity stream, grid, phase portraits, etc.)

    - __CellRank analysis:__
    Cell fate and lineage analysis using CellRank

    - __CellRank visualization:__
    Visualization of CellRank analysis results (terminal states, fate probabilities, gene trends, etc.)

    - __DeepVelo analysis:__
    Deep learning-based RNA velocity estimation

    - __UniTVelo analysis:__
    Unified-time RNA velocity analysis. Uses scVelo analysis results as input. Selectable Unified-time/Independent mode

    - __noSpliceVelo analysis:__
    RNA velocity analysis without spliced/unspliced separation. Can analyze with standard scRNA-seq count data only

    - __Pseudotime gene expression:__
    Visualization of gene expression trends along pseudotime and cluster density

    - __Dynamo analysis:__
    Advanced RNA velocity analysis and vector field reconstruction using Dynamo

    - __Dynamo visualization:__
    Visualization of Dynamo analysis results (streamline, vector field, topography, geometric features, etc.)

    - __Dynamo perturbation:__
    Gene perturbation simulation analysis

    - __Dynamo LAP:__
    Least Action Path analysis for computing optimal transition paths between cell states

    - __TFvelo analysis:__
    Transcription factor velocity analysis

    - __TFvelo to Dynamo:__
    Convert TFvelo results to Dynamo format for Vector Field construction and perturbation analysis

    - __VeloViz:__
    2D embedding considering RNA velocity information. Unlike traditional UMAP/tSNE, differentiation direction is more clearly visualized

    #### Cell communication
    - __LIANA LR analysis:__
    Ligand-receptor analysis using LIANA+

    - __LIANA comparison:__
    LIANA condition comparison

    - __CellChat:__
    CellChat analysis from h5ad files - faithful Python implementation of R version

    - __CellChat comparison:__
    Two-condition comparison from h5ad files containing both conditions

    - __CellChat permutation test:__
    Statistical testing of two-condition comparison using permutation test

    - __CellChat R qs to python:__
    Convert SCALA CellChat analysis result qs files for use in METIS

    #### scRNA file operation
    - __Download public data for SCALA and cellxgene:__
    Download files analyzable in SCALA/cellxgene. Can download from GEO and other sources if download links are available

    - __File explorer:__
    Browse, delete, and download files accessible by SCALA

    - __File uploader:__
    Upload files accessible by SCALA/cellxgene

    - __SCALA file browser:__
    Browse files accessible by SCALA

    - __Data file browser:__
    Browse and manage data files

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
    Filtering by peak length and score.
    \tPeak length distribution is also displayed

    - __Bam to DESeq2:__
    Normalize CUT&RUN peak counts based on greenlist counts for DESeq2 analysis

    - __Filter out blacklist:__
    Filter peaks overlapping with blacklist regions

    - __Convert bed to fasta:__
    Can also handle MACS peak files

    - __Denoise bedgraph bigwig:__
    Denoise CUT&RUN data (bedgraph/bigwig format)

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
        st.Page("pages/combat-seq.py", title="Batch removal by Combat-seq"),
        st.Page("pages/filter_xy_genes.py", title="Filter X/Y Chromosome Genes"),
    ],
    "Gene name conversion 🔁": [
        st.Page("pages/2_Update Gene Symbol.py", title="Update Gene Symbol"),
        st.Page("pages/3_Homer to Ensembl.py", title="Gene Symbol to Ensembl ID"),
        st.Page("pages/Ensembl2Symbol.py", title="Ensembl ID to Gene Symbol"),
        st.Page("pages/mouse-human.py", title="Mouse Human Symbol Conversion"),
        st.Page("pages/gmt-mouse-human.py", title="GMT Mouse Human Conversion"),
        st.Page("pages/microarray_gene_filter.py", title="Microarray Gene Name Filter"),
    ],
    "Data Analysis 🧮": [
        st.Page("pages/6_Calc_DESeq2.py", title="DESeq2 etc"),
        st.Page("pages/Calc_DESeq2_LRT.py", title="DESeq2-LRT etc"),
        st.Page("pages/edgeR.py", title="edgeR"),
        st.Page("pages/limma.py", title="limma"),
        st.Page("pages/DE_rpy2.py", title="DE method comparison"),
        st.Page("pages/permutation_test.py", title="Permutation test"),
        st.Page("pages/8_rnkgene.py", title="Make ranking file for GSEA"),
        st.Page("pages/impulsede2-streamlit-app.py", title="impulseDE2"),
        st.Page("pages/rnaseq-mfuzz-streamlit-app.py", title="Mfuzz clustering"),
        st.Page("pages/CompareDE.py", title="DEA result comparison"),
        st.Page("pages/analyze_freq.py", title="Compare ratios"),
        st.Page("pages/KStest.py", title="KS Test Distribution Analysis"),
    ],
    "Data Visualization 🌋": [
        st.Page("pages/pca.py", title="PCA"),
        st.Page("pages/pca_statistics.py", title="PCA statistics"),
        st.Page("pages/7_Volcano Plot.py", title="Volcano plot"),
        st.Page("pages/Heatmap.py", title="Heatmap"),
        st.Page("pages/Boxplot.py", title="Box_Violin plot"),
        st.Page("pages/venn_upset.py", title="Venn_Upset Plot"),
    ],
    "Pathway Analysis 🔀": [
        st.Page("pages/decoupler.py", title="decoupleR"),
        st.Page("pages/GSEApy.py", title="GSEApy"),
        st.Page("pages/aPEAR.py", title="aPEAR Enrichment Network"),
        st.Page("pages/PADOG.py", title="PADOG"),
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
        st.Page("pages/dimension_reduction.py", title="Dimension reduction"),
        st.Page("pages/pseudobulk.py", title="Pseudobulk"),
        st.Page("pages/SEAcells.py", title="Metacells by SEACells"),
        st.Page("pages/pseudoreplicates.py", title="Random pseudo-replicates"),
        st.Page("pages/pseudobulk_DEA.py", title="Pseudobulk DEA (RUV)"),
        st.Page("pages/Memento.py", title="memento DE analysis"),
        st.Page("pages/Memento2D.py", title="memento 2D analysis"),
        st.Page("pages/sccoda_analysis.py", title="scCODA compositional analysis"),
        st.Page("pages/add_metadata.py", title="Add metadata to h5ad"),
        st.Page("pages/celltypist.py", title="CellTypist annotation"),
        st.Page("pages/compass_analysis.py", title="COMPASS metabolic analysis"),
        st.Page("pages/compass_postprocess.py", title="COMPASS postprocess"),
        st.Page("pages/scFEA_analysis.py", title="scFEA metabolic analysis"),
        st.Page("pages/scFEA_postprocess.py", title="scFEA postprocess"),
        st.Page("pages/geneformer_analysis.py", title="Geneformer perturbation"),
    ],
    "SCENIC 🎭": [
        st.Page("pages/SCENICviewer.py", title="SCENIC heatmap"),
        st.Page("pages/prepare_scenic_data.py", title="Prepare regulon data for heatmap"),
        st.Page("pages/SCENICcsi.py", title="SCENIC CSI"),
        st.Page("pages/SCENIC_network.py", title="SCENIC network analysis"),
        st.Page("pages/SCENIC_multinetwork.py", title="SCENIC multinetwork analysis"),
        st.Page("pages/SCENIC_network2.py", title="SCENIC network analysis 2"),
    ],
    "RNA velocity 🚀": [
        st.Page("pages/velocity_filter.py", title="Data filtering"),
        st.Page("pages/velocity_analysis.py", title="scVelo analysis"),
        st.Page("pages/velocity_visualization.py", title="scVelo visualization"),
        st.Page("pages/cellrank_analysis.py", title="CellRank analysis"),
        st.Page("pages/cellrank_visualization.py", title="CellRank visualization"),
        st.Page("pages/deepvelo_analysis.py", title="DeepVelo analysis"),
        st.Page("pages/unitvelo_analysis.py", title="UniTVelo analysis"),
        st.Page("pages/nosplicevelo_analysis.py", title="noSpliceVelo analysis"),
        st.Page("pages/pseudotime_gene_expression.py", title="Pseudotime gene expression"),
        st.Page("pages/dynamo_analysis.py", title="Dynamo analysis"),
        st.Page("pages/dynamo_visualization.py", title="Dynamo visualization"),
        st.Page("pages/dynamo_perturbation_v2.py", title="Dynamo perturbation"),
        st.Page("pages/dynamo_lap_correct.py", title="Dynamo LAP"),
        st.Page("pages/TFvelo.py", title="TFvelo analysis"),
        st.Page("pages/tfvelo_to_dynamo.py", title="TFvelo to Dynamo"),
        st.Page("pages/veloviz.py", title="VeloViz"),
    ],
    "Cell communication 💬": [
        st.Page("pages/liana_steady.py", title="LIANA LR analysis"),
        st.Page("pages/liana_comparison.py", title="LIANA comparison"),
        st.Page("pages/cellchat.py", title="CellChat"),
        st.Page("pages/cellchat_comparison.py", title="CellChat comparison"),
        st.Page("pages/cellchat_permutation.py", title="CellChat permutation test"),
        st.Page("pages/cellchatR2py.py", title="CellChat R qs to python"),
    ],
    "scRNA file operation 🗄": [
        st.Page("pages/download.py", title="Download public data for SCALA and cellxgene"),
        st.Page("pages/fileexplorer.py", title="File explorer"),
        st.Page("pages/ftp.py", title="File uploader"),
        st.Page("pages/filebrowser.py", title="SCALA file browser"),
        st.Page("pages/data_file_browser.py", title="Data file browser"),
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
        st.Page("pages/paperqa2.py", title="PaperQA2"),
        st.Page("pages/grants.py", title="Grants eval"),
        st.Page("pages/tts_generator.py", title="TTS Generator"),
        st.Page("pages/form408_to_ms_json.py", title="MS JSON converter"),
    ],
})

# Run navigation
pg.run()
