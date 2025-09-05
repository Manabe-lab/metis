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

| Tool | Features & Parameters | Input/Output |
|------|---------------------|--------------|
| **Update Gene Symbol** | • Automatic symbol updating to latest HGNC/MGI standards<br>• Batch processing of multiple files<br>• Deprecated symbol detection and replacement<br>• Synonym resolution<br>• Quality control reports<br>• Custom symbol mapping support | Input: Homer files, TSV with symbols<br>Output: Updated symbol files |
| **Gene Symbol to Ensembl ID** | • Human and mouse symbol conversion<br>• Multiple database sources (Ensembl, NCBI)<br>• Version-specific mappings<br>• One-to-many relationship handling<br>• Canonical transcript selection<br>• Missing symbol reporting | Input: Symbol lists, Homer files<br>Output: Ensembl ID mappings |
| **Ensembl ID to Gene Symbol** | • Reverse ID conversion with validation<br>• Duplicate gene handling strategies<br>• Biotype filtering (protein-coding, lncRNA, etc.)<br>• Version compatibility checking<br>• Bulk conversion processing<br>• Cross-reference validation | Input: Ensembl ID lists<br>Output: Symbol mappings |
| **Mouse Human Symbol Conversion** | • Ortholog mapping using multiple databases<br>• One-to-one and one-to-many conversions<br>• Confidence score reporting<br>• Evolutionary relationship validation<br>• Batch processing capabilities<br>• Custom ortholog definitions | Input: Mouse/Human symbols<br>Output: Cross-species mappings |
| **GMT Mouse Human Conversion** | • Gene set file format conversion<br>• Species-specific pathway translation<br>• GMT file validation and repair<br>• Batch processing multiple gene sets<br>• Coverage statistics reporting<br>• Custom gene set creation | Input: GMT files<br>Output: Converted gene sets |

## 🧮 Data Analysis

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **DESeq2** | • Multiple experimental design support<br>• GLM-based differential analysis<br>• Beta regression for proportion data<br>• Batch effect correction<br>• Multiple testing correction (FDR, Bonferroni)<br>• Interactive parameter tuning<br>• Results filtering and export | Input: Count matrix + metadata<br>Output: DE results, MA plots |
| **DESeq2-LRT** | • Likelihood ratio tests for complex designs<br>• Time course analysis capabilities<br>• ANOVA-like multi-group comparisons<br>• Model comparison statistics<br>• Trajectory analysis<br>• Custom contrast definitions | Input: Count data + time/group info<br>Output: LRT results, time profiles |
| **edgeR** | • TMM normalization<br>• Exact tests for small samples<br>• Quasi-likelihood F-tests<br>• Robust dispersion estimation<br>• Multi-dimensional scaling plots<br>• FDR control options | Input: Count matrix<br>Output: DE genes, dispersion plots |
| **limma** | • Linear modeling framework<br>• Empirical Bayes moderation<br>• Voom transformation for RNA-seq<br>• Duplicate correlation handling<br>• Contrast matrix definitions<br>• Gene set testing integration | Input: Expression matrix<br>Output: Moderated t-statistics |
| **DE method comparison** | • Side-by-side method evaluation<br>• Concordance analysis<br>• Venn diagram generation<br>• ROC curve analysis<br>• Benchmark statistics<br>• Method recommendation engine | Input: Multiple DE results<br>Output: Comparison reports |
| **Permutation test** | • Freedman-Lane permutation procedure<br>• Batch-aware randomization<br>• Multiple correction strategies<br>• Bootstrap confidence intervals<br>• Non-parametric hypothesis testing<br>• Custom test statistics | Input: Expression + design matrix<br>Output: Permutation p-values |
| **Make ranking file for GSEA** | • Multiple ranking metrics (t-stat, log2FC, p-value)<br>• Custom gene filtering<br>• Tie handling strategies<br>• Format validation<br>• Batch processing multiple comparisons | Input: DE results<br>Output: Ranked gene lists (.rnk) |
| **Batch removal by Combat-seq** | • Count-based batch correction<br>• Biological variation preservation<br>• Multiple batch variables<br>• Before/after visualization<br>• Quality assessment metrics | Input: Count matrix + batch info<br>Output: Batch-corrected counts |
| **impulseDE2** | • Impulse model fitting<br>• Time course DE detection<br>• Pattern classification<br>• Trajectory clustering<br>• Pseudotime analysis<br>• Model selection criteria | Input: Time course data<br>Output: Impulse model results |
| **Mfuzz clustering** | • Soft clustering algorithm<br>• Fuzzy membership values<br>• Optimal cluster number determination<br>• Time profile visualization<br>• Cluster validation metrics<br>• Gene annotation integration | Input: Time series expression<br>Output: Cluster assignments |
| **DEA result comparison** | • Multi-study result integration<br>• Meta-analysis capabilities<br>• Effect size comparisons<br>• Heterogeneity assessment<br>• Publication bias testing<br>• Forest plot generation | Input: Multiple DE result sets<br>Output: Meta-analysis results |
| **Compare ratios** | • Logit transformation<br>• Beta regression modeling<br>• Proportion data analysis<br>• Confidence interval estimation<br>• Model diagnostics<br>• Effect size calculations | Input: Ratio/proportion data<br>Output: Statistical test results |

## 🌋 Data Visualization

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **PCA** | • Multiple dimensionality reduction methods<br>• Interactive 2D/3D plotting<br>• Variance explained calculations<br>• Sample clustering analysis<br>• Batch effect visualization<br>• Custom color/shape mappings<br>• Loading plot generation | Input: Expression matrix<br>Output: PCA plots, coordinates |
| **Volcano plot** | • Interactive plotly-based visualization<br>• Customizable significance thresholds<br>• Gene label highlighting<br>• Color scheme options<br>• Export to multiple formats<br>• Zoom and pan functionality<br>• Statistical annotation | Input: DE results<br>Output: Interactive volcano plots |
| **Heatmap** | • Hierarchical clustering options<br>• Multiple distance metrics<br>• Row/column annotations<br>• Color palette customization<br>• Dendrogram display control<br>• Missing value handling<br>• Batch processing capabilities | Input: Expression matrix<br>Output: Clustered heatmaps |
| **Box/Violin plot** | • Multiple plot types (box, violin, strip)<br>• Statistical significance testing<br>• Group comparison options<br>• Outlier detection and highlighting<br>• Custom color schemes<br>• Batch processing<br>• Export options | Input: Expression + metadata<br>Output: Statistical plots |
| **Venn/Upset Plot** | • Multi-set intersection analysis<br>• Interactive UpSet plots<br>• Set size statistics<br>• Intersection export<br>• Custom set labeling<br>• Batch processing<br>• Statistical testing of overlaps | Input: Gene/feature lists<br>Output: Intersection visualizations |

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