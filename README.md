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

## 🧬 RNA-seq Workflow

### Standard Analysis Pipeline
1. **Count Data Normalization** → PCA, Heatmap (QC and data distribution check)
2. **DESeq2** → Volcano plot
3. **Gene set analysis**: DESeq2 → GSEApy, decoupleR
   - Alternatively, over-representation analysis using DEG lists (GSEApy, decoupleR)

### Usage Examples
- **Creating data for Heatmap**: Count Data Normalization. Or extract rlog data using "Homer to Data" from DESeq2 results. Filter by expression levels or changes using "Filter Log Z-score". Basic filtering can also be done in Heatmap.
- **Updating gene names from Homer etc**: Use "Update Gene Symbol". Since many gene names are updated, it's better to update depending on downstream analysis methods.
- **Converting gene names to Ensembl ID**: Use "Homer to Ensembl". For example, DIANE uses Ensembl IDs, so convert using this tool.
- **Data file organization**: Remove unnecessary columns/rows. Extract only needed columns. Change row/column names, reorder columns, transpose rows and columns using "Manipulate Data Table".

## 📝 Data Normalization/Manipulation

| Tool | Description |
|------|-------------|
| **Count Data Normalization** | RNA-seq count data normalization, QC PCA, heatmap, box plot and other preprocessing |
| **Homer or DESeq2 to Data** | Convert Homer output files to data-only files; Extract only rlog data from DESeq2 |
| **Count to TPM** | Convert raw count data to TPM |
| **Manipulate Data Table** | Delete/extract rows/columns from data table. Change column order, row names, column names |
| **Merge Data Files** | Merge data files |
| **Filter and Transform Data** | Data filtering. Data XY transposition |
| **Spread Sheet** | Spreadsheet editor |

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

## 🐛 Troubleshooting

### Common Issues
1. **Memory errors**: Some analyses require substantial memory, especially for large single-cell datasets
2. **R package dependencies**: Ensure all required R packages are installed
3. **File format compatibility**: Check that input files match expected formats

### Performance Tips
- Use filtering options to reduce dataset size before analysis
- Consider using pseudobulk analysis for very large single-cell datasets
- Utilize the batch processing features for multiple samples

## 🤝 Contributing

We welcome contributions! Please feel free to submit issues, feature requests, or pull requests.

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 📞 Support

If you encounter any issues or have questions, please:
1. Check the troubleshooting section above
2. Search existing issues on GitHub
3. Create a new issue with detailed information about your problem

## 🙏 Acknowledgments

METIS integrates many excellent bioinformatics tools and libraries. We thank all the developers and researchers who created these foundational tools that make METIS possible.