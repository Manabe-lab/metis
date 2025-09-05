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


## 📝 Data Normalization/Manipulation

| Tool | Features & Parameters | Input/Output |
|------|---------------------|--------------|
| **Count Data Normalization** | • Multiple normalization methods (TPM, FPKM, CPM, log2)<br>• QC metrics calculation<br>• PCA analysis with customizable components<br>• Interactive heatmaps with clustering options<br>• Box plots for sample comparison<br>• Batch effect visualization | Input: Raw count matrix<br>Output: Normalized data, QC plots |
| **Homer or DESeq2 to Data** | • Extract data from Homer differential analysis results<br>• Convert DESeq2 results to clean matrices<br>• Extract rlog/vst transformed data<br>• Automated column filtering and renaming<br>• Support for multiple comparison extractions | Input: Homer/DESeq2 files<br>Output: Clean data matrices |
| **Count to TPM** | • Transcript length-based normalization<br>• Built-in transcript length databases<br>• Custom length file support<br>• Batch processing capabilities<br>• Quality control metrics | Input: Count matrix + lengths<br>Output: TPM normalized data |
| **Manipulate Data Table** | • Advanced row/column operations<br>• Regular expression filtering<br>• Data type conversions<br>• Missing value handling<br>• Custom sorting and indexing<br>• Batch rename operations | Input: Any tabular data<br>Output: Processed data table |
| **Merge Data Files** | • Multi-file merging with flexible join types<br>• Automatic column alignment<br>• Duplicate handling strategies<br>• Missing data interpolation<br>• Memory-efficient processing | Input: Multiple data files<br>Output: Merged dataset |
| **Filter and Transform Data** | • Expression threshold filtering<br>• Variance-based feature selection<br>• Log transformations<br>• Z-score standardization<br>• Data transposition<br>• Outlier detection and removal | Input: Expression matrix<br>Output: Filtered/transformed data |
| **Spread Sheet** | • Excel-like interface for data editing<br>• Real-time data validation<br>• Formula calculations<br>• Import/export multiple formats<br>• Undo/redo functionality | Input: Various formats<br>Output: Edited data files |

## 🔁 Gene Name Conversion

| Tool | Description |
|------|-------------|
| **Update Gene Symbol** | Update gene symbols in Homer output or general data files (first column is Symbol) |
| **Gene Symbol to Ensembl ID** | Convert gene symbols to Ensembl IDs in Homer output or general data files (for DIANE file creation) |
| **Ensembl ID to Gene Symbol** | Handle duplicate genes as well |
| **Mouse Human Symbol Conversion** | Convert between mouse and human gene symbols |
| **GMT Mouse Human Conversion** | Convert gene symbols in gene set files (gmt files) used in GSEA. Format validation for gmt files is also available |

## 🧮 Data Analysis

| Tool | Description |
|------|-------------|
| **DESeq2** | DESeq2, limma-eBayes, Beta regression, GLM group comparisons |
| **DESeq2-LRT** | ANOVA-like, time series analysis etc. limma-eBayes for non-count data analysis. Beta regression for 0-1 data |
| **edgeR** | edgeR differential expression analysis |
| **limma** | limma voom calculation |
| **DE method comparison** | Compare DEG calculation methods (DESeq2, edgeR, limma) results |
| **Permutation test** | Freedman-Lane test including batch effects |
| **Make ranking file for GSEA** | Create GSEA ranking files |
| **Batch removal by Combat-seq** | Batch removal using Combat-seq |
| **impulseDE2** | Time course analysis using impulseDE2 |
| **Mfuzz clustering** | Time course data clustering using Mfuzz |
| **DEA result comparison** | Compare significantly different genes from DEA results |
| **Compare ratios** | Statistical analysis of ratio data using t-test and β regression on logit transformed data |

## 🌋 Data Visualization

| Tool | Description |
|------|-------------|
| **PCA** | PCA, UMAP, tSNE, MDS |
| **Volcano plot** | Interactive volcano plots for differential expression results |
| **Heatmap** | Customizable heatmaps with clustering |
| **Box/Violin plot** | Box plots and violin plots for expression visualization |
| **Venn/Upset Plot** | Venn diagram and UpSet plot for set comparisons |

## 🔀 Pathway Analysis

| Tool | Description |
|------|-------------|
| **decoupleR** | Signal pathway, TF activity analysis, and GSEA using mSigDB. Can create publication quality enrichment plots |
| **GSEApy** | Gene Set Enrichment Analysis |
| **PPI analysis** | Protein-protein interaction analysis using STRING-db |

## 🥅 WGCNA (Weighted Gene Co-expression Network Analysis)

| Tool | Description |
|------|-------------|
| **WGCNA** | WGCNA analysis using WGCNApy or R WGCNA |
| **WGCNA network plot** | Network visualization of WGCNA modules |
| **WGCNA hub UMAP** | UMAP visualization of WGCNA module genes |
| **WGCNA objects comparison** | Display relationships between WGCNA modules |
| **Generate gmt from cluster info** | Create gmt files from cluster results or WGCNA module genes |

## 🎡 Single-cell RNA-seq (scRNA-seq)

| Tool | Description |
|------|-------------|
| **Pseudobulk** | Create pseudobulk data from anndata (h5ad) |
| **Metacells by SEACells** | Create metacells using SEACells |
| **Random pseudo-replicates** | Create pseudo-replicates by random cell splitting |
| **memento DE analysis** | Differential expression analysis using memento |
| **memento 2D analysis** | Two-dimensional analysis using memento |
| **SCENIC heatmap** | Heatmap visualization of SCENIC gene regulatory networks |
| **Prepare regulon data for heatmap** | Prepare SCENIC regulon data for visualization |
| **SCENIC CSI** | Calculate connection specificity index (CSI) of SCENIC regulons |
| **SCENIC network analysis** | Network visualization of SCENIC regulons centered on transcription factors or their targets |
| **SCENIC multinetwork analysis** | Network visualization of SCENIC regulons centered on multiple transcription factors and targets |

## 💬 Cell Communication

| Tool | Description |
|------|-------------|
| **LIANA LR analysis** | Ligand-receptor analysis using LIANA+ |
| **LIANA comparison** | Compare LIANA results between conditions |
| **CellChat** | CellChat analysis from h5ad files - faithful Python implementation of R version |
| **CellChat comparison** | Two-condition comparison from h5ad files containing both conditions |
| **CellChat permutation test** | Statistical testing of two-condition comparison using permutation test |
| **CellChat R qs to python** | Convert SCALA CellChat analysis result qs files for use in metis |

## 🗄 scRNA File Operation

| Tool | Description |
|------|-------------|
| **Download public data** | Download files analyzable in SCALA/cellxgene. Can download from GEO and other sources if download links are available |
| **File explorer** | Browse, delete, and download files accessible by SCALA |
| **File uploader** | Upload files accessible by SCALA/cellxgene |
| **SCALA file browser** | Browse files accessible by SCALA |

## 🧬 ChIP-seq Analysis

| Tool | Description |
|------|-------------|
| **Sort BAM file** | Sort BAM files in METIS_data directory |
| **Merge BAM files** | Merge BAM files in METIS_data directory |
| **Bam to bedGraph for SEACR** | Convert BAM files to bedGraph format for SEACR |
| **SEACR peak calling** | Peak calling using SEACR algorithm |
| **Macs3 peak calling** | Peak calling using MACS3 |
| **Annotating and filtering peaks** | Annotation and annotation-based filtering of peak files and bed files |
| **Bed length/score filter** | Filtering by peak length and score. Peak length distribution is also displayed |
| **Bam to DESeq2** | Normalize CUT&RUN peak counts based on greenlist counts for DESeq2 analysis |
| **Blacklist filter** | Filter peaks overlapping with blacklist regions |
| **Bed to fasta** | Can also handle MACS peak files |
| **Denoise bedgraph bigwig** | Denoise bedGraph and bigWig files |

## Ⓜ Miscellaneous Tools

| Tool | Description |
|------|-------------|
| **Merge excel files** | Combine multiple Excel files |
| **Remove duplicates** | Remove duplicate entries from datasets |
| **Split data file on key** | Split data files based on key values in another file |

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

METIS integrates many excellent bioinformatics tools and libraries. We thank all the developers and researchers who created these foundational tools that make METIS possible.