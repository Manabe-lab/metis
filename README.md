# 🧝🏻‍♀️ METIS

## **M**olecular **E**xploration and **T**ranscriptomic **I**nvestigation **S**uite

METIS is a comprehensive bioinformatics web application suite built with Streamlit for RNA-seq analysis, single-cell RNA-seq analysis, ChIP-seq analysis, and various data manipulation tasks. This toolkit provides an intuitive interface for researchers to perform complex genomic analyses without requiring extensive programming knowledge.

## 🚀 Quick Start

### Installation

1. Clone this repository:
```bash
git clone https://github.com/yourusername/metis.git
cd metis
```

2. Install required packages:
```bash
pip install -r requirements.txt
```

3. Run the main application:
```bash
streamlit run metis.py
```

## 📋 Required Packages

### Core Dependencies
- **streamlit** - Web application framework
- **pandas** - Data manipulation and analysis
- **numpy** - Numerical computing
- **matplotlib** - Plotting library
- **seaborn** - Statistical data visualization
- **plotly** - Interactive plotting
- **scipy** - Scientific computing

### Bioinformatics Libraries
- **scanpy** - Single-cell analysis
- **anndata** - Annotated data structures
- **pyscenic** - Gene regulatory network analysis
- **decoupler** - Pathway activity inference
- **GSEApy** - Gene Set Enrichment Analysis
- **PyWGCNA** - Weighted Gene Co-expression Network Analysis
- **liana** - Ligand-receptor analysis
- **cellchat** - Cell communication analysis

### R Integration
- **rpy2** - R interface for Python
- **R packages**: DESeq2, edgeR, limma, WGCNA, Combat-seq, impulseDE2, Mfuzz

### File Processing
- **openpyxl** - Excel file handling
- **h5py** - HDF5 file format
- **pysam** - SAM/BAM file processing
- **pybedtools** - BED file manipulation

## 📊 Available Applications Overview

METIS includes **67 specialized applications** organized into 9 functional categories:

| Category | Apps | Primary Functions |
|----------|------|-------------------|
| **📝 Data Normalization/Manipulation** | 7 apps | Raw data processing, format conversion, table manipulation |
| **🔁 Gene Name Conversion** | 5 apps | ID mapping, symbol updates, cross-species conversion |
| **🧮 Data Analysis** | 12 apps | Differential expression, statistical testing, batch correction |
| **🌋 Data Visualization** | 5 apps | PCA, volcano plots, heatmaps, box plots, Venn diagrams |
| **🔀 Pathway Analysis** | 3 apps | Gene set enrichment, pathway activity, protein interactions |
| **🥅 WGCNA** | 5 apps | Co-expression networks, module detection, hub analysis |
| **🎡 scRNA-seq** | 11 apps | Single-cell analysis, pseudobulk, metacells, SCENIC |
| **💬 Cell Communication** | 6 apps | Ligand-receptor analysis, CellChat, intercellular signaling |
| **🧬 ChIP-seq** | 11 apps | Peak calling, annotation, BAM processing, quality control |
| **Ⓜ Miscellaneous** | 3 apps | File operations, data cleaning, utility functions |

### 🎯 Key Application Highlights

**Most Popular Tools:**
- **DESeq2** - Industry-standard differential expression analysis
- **GSEApy** - Comprehensive gene set enrichment analysis
- **CellChat** - Advanced cell-cell communication analysis
- **SCENIC** - Gene regulatory network reconstruction
- **PyWGCNA** - Co-expression network analysis

**Data Processing Pipeline:**
1. **Import & Normalize** → Count Data Normalization, TPM conversion
2. **Quality Control** → PCA, heatmaps, filtering
3. **Statistical Analysis** → DESeq2, edgeR, limma comparisons
4. **Pathway Analysis** → GSEA, decoupleR enrichment
5. **Visualization** → Publication-ready plots and figures

**Specialized Workflows:**
- **Single-cell RNA-seq**: Complete pipeline from raw counts to cell communication
- **ChIP-seq/CUT&RUN**: Peak calling to differential binding analysis
- **Time-course**: impulseDE2 and Mfuzz clustering analysis
- **Multi-omics**: Integration tools for diverse data types

## 📝 Data Normalization/Manipulation

| Tool | Features & Parameters | Input/Output |
|------|---------------------|--------------|
| **Count Data Normalization** | **Analysis:** RNA-seq count normalization and QC<br>**Key Features:** Multiple methods (TPM, FPKM, CPM), PCA, heatmaps, batch effect detection | Input: Raw count matrix<br>Output: Normalized data, QC plots |
| **Homer or DESeq2 to Data** | **Analysis:** Extract clean data from Homer/DESeq2 results<br>**Key Features:** rlog/vst extraction, automated formatting, multi-comparison support | Input: Homer/DESeq2 files<br>Output: Clean data matrices |
| **Count to TPM** | **Analysis:** Convert raw counts to TPM normalization<br>**Key Features:** Built-in transcript lengths, custom length files, batch processing | Input: Count matrix + lengths<br>Output: TPM normalized data |
| **Manipulate Data Table** | **Analysis:** Data table manipulation and restructuring<br>**Key Features:** Row/column operations, regex filtering, missing value handling | Input: Any tabular data<br>Output: Processed data table |
| **Merge Data Files** | **Analysis:** Combine multiple data files<br>**Key Features:** Flexible join types, automatic alignment, duplicate handling | Input: Multiple data files<br>Output: Merged dataset |
| **Filter and Transform Data** | **Analysis:** Data filtering and mathematical transformations<br>**Key Features:** Expression thresholds, variance filtering, log/Z-score transforms | Input: Expression matrix<br>Output: Filtered/transformed data |
| **Spread Sheet** | **Analysis:** Interactive data editing interface<br>**Key Features:** Excel-like editing, real-time validation, multiple formats | Input: Various formats<br>Output: Edited data files |

## 🔁 Gene Name Conversion

| Tool | Features & Parameters | Input/Output |
|------|---------------------|--------------|
| **Update Gene Symbol** | **Analysis:** Update gene symbols to latest HGNC/MGI standards<br>**Key Features:** Batch processing, deprecated symbol detection, synonym resolution | Input: Homer files, TSV with symbols<br>Output: Updated symbol files |
| **Gene Symbol to Ensembl ID** | **Analysis:** Convert gene symbols to Ensembl IDs<br>**Key Features:** Human/mouse support, one-to-many handling, canonical transcript selection | Input: Symbol lists, Homer files<br>Output: Ensembl ID mappings |
| **Ensembl ID to Gene Symbol** | **Analysis:** Convert Ensembl IDs to gene symbols<br>**Key Features:** Duplicate handling, biotype filtering, cross-reference validation | Input: Ensembl ID lists<br>Output: Symbol mappings |
| **Mouse Human Symbol Conversion** | **Analysis:** Cross-species ortholog mapping<br>**Key Features:** Multiple databases, confidence scores, one-to-many handling | Input: Mouse/Human symbols<br>Output: Cross-species mappings |
| **GMT Mouse Human Conversion** | **Analysis:** Convert gene sets between mouse and human<br>**Key Features:** GMT file validation, batch processing, coverage statistics | Input: GMT files<br>Output: Converted gene sets |

## 🧮 Data Analysis

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **DESeq2** | **Analysis Types:**<br>• Pairwise comparisons (treatment vs control)<br>• Multi-group comparisons (one-way ANOVA)<br>• Multi-factorial designs (treatment + batch effects)<br>• Time course analysis<br>• Dose-response studies<br>**Test Methods:**<br>• Wald test for single comparisons<br>• GLM-based modeling with covariates<br>• Beta regression for proportion data<br>• FDR and Bonferroni correction | Input: Count matrix + metadata<br>Output: DE results, MA plots |
| **DESeq2-LRT** | **Analysis Types:**<br>• Multi-group time course experiments<br>• Complex experimental designs with multiple factors<br>• Nested model comparisons<br>• Non-linear time trajectory analysis<br>**Test Methods:**<br>• Likelihood Ratio Test (LRT) for model comparison<br>• ANOVA-like F-test for multiple groups<br>• Chi-square test for nested models<br>• Custom contrast matrix testing<br>• Joint significance testing across conditions | Input: Count data + time/group info<br>Output: LRT results, time profiles |
| **edgeR** | **Analysis Types:**<br>• Two-group comparisons (classical exact test)<br>• Multi-group designs with GLM framework<br>• Paired sample analysis<br>• Time course and factorial experiments<br>**Test Methods:**<br>• Fisher's exact test for small samples<br>• Quasi-likelihood F-test for complex designs<br>• Likelihood ratio test for nested models<br>• TMM normalization with robust dispersion estimation<br>• Empirical Bayes moderated statistics | Input: Count matrix + design<br>Output: DE genes, dispersion plots |
| **limma** | **Analysis Types:**<br>• Microarray differential expression<br>• RNA-seq analysis via voom transformation<br>• Multi-factorial experimental designs<br>• Paired and repeated measures analysis<br>• Meta-analysis across studies<br>**Test Methods:**<br>• Moderated t-test with empirical Bayes<br>• F-test for multiple contrasts<br>• Linear mixed models for correlated samples<br>• Voom-limma pipeline for count data<br>• Robust regression options | Input: Expression/count matrix<br>Output: Moderated statistics, p-values |
| **DE method comparison** | **Analysis:** Compare DESeq2, edgeR, limma results<br>**Key Features:** Concordance analysis, Venn diagrams, ROC curves | Input: Multiple DE results<br>Output: Comparison reports |
| **Permutation test** | **Analysis:** Freedman-Lane permutation testing<br>**Key Features:** Batch-aware randomization, bootstrap confidence intervals | Input: Expression + design matrix<br>Output: Permutation p-values |
| **Make ranking file for GSEA** | **Analysis:** Create GSEA-compatible ranking files<br>**Key Features:** Multiple ranking metrics, tie handling, batch processing | Input: DE results<br>Output: Ranked gene lists (.rnk) |
| **Batch removal by Combat-seq** | **Analysis:** sva Combat-seq batch correction for counts<br>**Key Features:** Biological variation preservation, before/after visualization | Input: Count matrix + batch info<br>Output: Batch-corrected counts |
| **impulseDE2** | **Analysis:** R impulseDE2 time course analysis<br>**Key Features:** Impulse model fitting, pattern classification, trajectory clustering | Input: Time course data<br>Output: Impulse model results |
| **Mfuzz clustering** | **Analysis:** R Mfuzz soft clustering for time series<br>**Key Features:** Fuzzy membership values, optimal cluster determination | Input: Time series expression<br>Output: Cluster assignments |
| **DEA result comparison** | **Analysis:** Compare differential expression results across studies<br>**Key Features:** Meta-analysis, effect size comparisons, forest plots | Input: Multiple DE result sets<br>Output: Meta-analysis results |
| **Compare ratios** | **Analysis:** Statistical analysis of ratio/proportion data<br>**Key Features:** Beta regression, logit transformation, confidence intervals | Input: Ratio/proportion data<br>Output: Statistical test results |

## 🌋 Data Visualization

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **PCA** | **Analysis:** PCA, UMAP, tSNE, MDS dimensionality reduction<br>**Key Features:** Interactive 2D/3D plotting, batch effect visualization, loading plots | Input: Expression matrix<br>Output: PCA plots, coordinates |
| **Volcano plot** | **Analysis:** Interactive volcano plots for DE results<br>**Key Features:** Plotly-based, gene labeling, customizable thresholds | Input: DE results<br>Output: Interactive volcano plots |
| **Heatmap** | **Analysis:** Hierarchical clustering heatmaps<br>**Key Features:** Multiple distance metrics, row/column annotations, dendrograms | Input: Expression matrix<br>Output: Clustered heatmaps |
| **Box/Violin plot** | **Analysis:** Box, violin, and strip plots with statistics<br>**Key Features:** Group comparisons, outlier detection, statistical testing | Input: Expression + metadata<br>Output: Statistical plots |
| **Venn/Upset Plot** | **Analysis:** Multi-set intersection analysis<br>**Key Features:** Interactive UpSet plots, set size statistics, overlap testing | Input: Gene/feature lists<br>Output: Intersection visualizations |

## 🔀 Pathway Analysis

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **decoupleR** | **Analysis:** Python decoupleR pathway and TF activity inference<br>**Key Features:** mSigDB integration, multiple statistical methods, publication plots | Input: Expression matrix + gene sets<br>Output: Activity scores, plots |
| **GSEApy** | **Analysis:** GSEApy gene set enrichment analysis<br>**Key Features:** Pre-ranked and standard GSEA, multiple databases, leading edge genes | Input: Gene lists/expression + gene sets<br>Output: Enrichment results, plots |
| **PPI analysis** | **Analysis:** STRING-db protein-protein interaction networks<br>**Key Features:** Confidence filtering, network visualization, Cytoscape export | Input: Gene/protein lists<br>Output: Interaction networks, annotations |

## 🥅 WGCNA (Weighted Gene Co-expression Network Analysis)

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **WGCNA** | **Analysis:** PyWGCNA and R-WGCNA co-expression network analysis<br>**Key Features:** Module detection, trait correlation, hub gene identification | Input: Expression matrix + traits<br>Output: Modules, networks, correlations |
| **WGCNA network plot** | **Analysis:** WGCNA module network visualization<br>**Key Features:** Interactive networks, dendrograms, trait-module heatmaps | Input: WGCNA results<br>Output: Network visualizations |
| **WGCNA hub UMAP** | **Analysis:** UMAP visualization of WGCNA module genes<br>**Key Features:** Hub gene highlighting, module-specific plots, gene annotations | Input: WGCNA modules + expression<br>Output: UMAP plots, hub gene lists |
| **WGCNA objects comparison** | **Analysis:** Compare WGCNA modules between conditions/studies<br>**Key Features:** Module preservation, consensus modules, statistical testing | Input: Multiple WGCNA objects<br>Output: Preservation statistics, plots |
| **Generate gmt from cluster info** | **Analysis:** Convert clustering results to GMT gene sets<br>**Key Features:** WGCNA module conversion, format validation, batch processing | Input: Cluster assignments<br>Output: GMT gene set files |

## 🎡 Single-cell RNA-seq (scRNA-seq)

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Pseudobulk** | **Analysis:** Create pseudobulk data from single-cell AnnData<br>**Key Features:** Cell type aggregation, multiple methods, quality metrics | Input: AnnData (h5ad)<br>Output: Pseudobulk expression matrices |
| **Metacells by SEACells** | **Analysis:** SEACells metacell construction algorithm<br>**Key Features:** Cell similarity grouping, size balancing, parameter optimization | Input: AnnData with embeddings<br>Output: Metacell assignments, objects |
| **Random pseudo-replicates** | **Analysis:** Create pseudo-replicates by random cell splitting<br>**Key Features:** Statistical power enhancement, bootstrap approaches, balance optimization | Input: Single-cell data<br>Output: Pseudo-replicated datasets |
| **memento DE analysis** | **Analysis:** Memento single-cell differential expression<br>**Key Features:** Multi-condition comparisons, batch handling, effect size estimation | Input: AnnData + conditions<br>Output: DE results, statistics |
| **memento 2D analysis** | **Analysis:** Memento 2D trajectory and spatial analysis<br>**Key Features:** Gradient-based testing, pseudotime integration, spatial patterns | Input: Spatial/trajectory data<br>Output: 2D analysis results |
| **SCENIC heatmap** | **Analysis:** Visualize SCENIC regulon activity patterns<br>**Key Features:** Hierarchical clustering, cell type annotation, interactive plots | Input: SCENIC regulon activities<br>Output: Activity heatmaps |
| **Prepare regulon data for heatmap** | **Analysis:** Prepare SCENIC results for heatmap visualization<br>**Key Features:** Activity score calculation, normalization, format conversion | Input: Raw SCENIC results<br>Output: Processed regulon data |
| **SCENIC CSI** | **Analysis:** Calculate Connection Specificity Index for SCENIC regulons<br>**Key Features:** Connectivity metrics, statistical testing, comparative analysis | Input: SCENIC network data<br>Output: CSI scores, statistics |
| **SCENIC network analysis** | **Analysis:** SCENIC gene regulatory network visualization<br>**Key Features:** TF-target networks, centrality measures, Cytoscape export | Input: SCENIC regulons<br>Output: Network visualizations |
| **SCENIC multinetwork analysis** | **Analysis:** Multi-TF regulatory network integration<br>**Key Features:** Comparative analysis, regulatory cascades, cross-condition comparison | Input: Multiple regulon sets<br>Output: Integrated network analysis |

## 💬 Cell Communication

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **LIANA LR analysis** | **Analysis:** LIANA+ ligand-receptor interaction analysis<br>**Key Features:** Multiple databases, statistical methods, cell type-specific analysis | Input: AnnData with cell types<br>Output: L-R interaction results |
| **LIANA comparison** | **Analysis:** Compare LIANA results between conditions<br>**Key Features:** Statistical testing, effect size calculations, visualization | Input: Multiple LIANA results<br>Output: Comparative analysis |
| **CellChat** | **Analysis:** Python CellChat cell-cell communication analysis<br>**Key Features:** Signaling pathway database, network centrality, pattern recognition | Input: AnnData (h5ad)<br>Output: Communication networks, plots |
| **CellChat comparison** | **Analysis:** Compare CellChat results between conditions<br>**Key Features:** Network topology changes, pathway comparisons, effect quantification | Input: Dual-condition h5ad<br>Output: Differential communication |
| **CellChat permutation test** | **Analysis:** Statistical testing of CellChat comparisons<br>**Key Features:** Permutation-based testing, bootstrap intervals, multiple correction | Input: CellChat comparison results<br>Output: Statistical significance |
| **CellChat R qs to python** | **Analysis:** Convert SCALA CellChat results to Python format<br>**Key Features:** Cross-platform compatibility, format conversion, data validation | Input: R CellChat .qs files<br>Output: Python-compatible data |


## 🧬 ChIP-seq Analysis

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Sort BAM file** | **Analysis:** SAMtools coordinate-based BAM sorting<br>**Key Features:** Memory optimization, multi-threading, automatic indexing | Input: Unsorted BAM files<br>Output: Sorted, indexed BAM |
| **Merge BAM files** | **Analysis:** SAMtools BAM file merging<br>**Key Features:** Header compatibility, quality preservation, duplicate marking | Input: Multiple BAM files<br>Output: Merged BAM file |
| **Bam to bedGraph for SEACR** | **Analysis:** Convert BAM to SEACR-compatible bedGraph<br>**Key Features:** Fragment extension, normalization options, strand-specific processing | Input: BAM files<br>Output: bedGraph files |
| **SEACR peak calling** | **Analysis:** SEACR CUT&RUN peak calling<br>**Key Features:** Stringent/relaxed modes, control integration, FDR thresholding | Input: bedGraph + control<br>Output: Peak files |
| **Macs3 peak calling** | **Analysis:** MACS3 model-based peak calling<br>**Key Features:** Broad/narrow modes, statistical testing, parameter optimization | Input: BAM/bedGraph files<br>Output: Peak calls, summits |
| **Annotating and filtering peaks** | **Analysis:** Peak annotation and genomic feature filtering<br>**Key Features:** Distance-based filtering, functional annotation, statistical enrichment | Input: Peak files (BED)<br>Output: Annotated peaks |
| **Bed length/score filter** | **Analysis:** Filter peaks by length and score thresholds<br>**Key Features:** Distribution visualization, custom criteria, quality metrics | Input: Peak/BED files<br>Output: Filtered peaks, plots |
| **Bam to DESeq2** | **Analysis:** Create count matrices for CUT&RUN DESeq2 analysis<br>**Key Features:** Greenlist normalization, quality integration, statistical prep | Input: BAM + peak files<br>Output: Count matrices |
| **Blacklist filter** | **Analysis:** Remove peaks overlapping blacklist regions<br>**Key Features:** Genome blacklist integration, overlap statistics, quality assessment | Input: Peak/BED files<br>Output: Filtered peak files |
| **Bed to fasta** | **Analysis:** Extract sequences from BED coordinates<br>**Key Features:** Multiple genomes, custom flanking, format validation | Input: BED files + genome<br>Output: FASTA sequences |
| **Denoise bedgraph bigwig** | **Analysis:** Signal smoothing and noise reduction<br>**Key Features:** Multiple algorithms, quality assessment, format conversion | Input: bedGraph/bigWig<br>Output: Denoised signal files |

## Ⓜ Miscellaneous Tools

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Merge excel files** | **Analysis:** Combine multiple Excel files<br>**Key Features:** Multi-sheet processing, header alignment, data type preservation | Input: Multiple Excel files<br>Output: Merged Excel/CSV |
| **Remove duplicates** | **Analysis:** Remove duplicate entries from datasets<br>**Key Features:** Multiple criteria, configurable rules, statistics reporting | Input: Any tabular data<br>Output: Deduplicated datasets |
| **Split data file on key** | **Analysis:** Split data based on key values<br>**Key Features:** Custom criteria, metadata preservation, multiple output formats | Input: Data file + key file<br>Output: Split data files |

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

METIS integrates many excellent bioinformatics tools and libraries. We thank all the developers and researchers who created these foundational tools that make METIS possible.