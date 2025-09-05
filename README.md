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

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **decoupleR** | • Multi-omics pathway activity inference<br>• TF activity estimation from expression<br>• Multiple statistical methods (mlm, norm, gsva)<br>• mSigDB integration<br>• Custom gene set support<br>• Publication-quality enrichment plots<br>• Batch processing capabilities | Input: Expression matrix + gene sets<br>Output: Activity scores, plots |
| **GSEApy** | • Pre-ranked and standard GSEA<br>• Multiple gene set databases<br>• Enrichment plot generation<br>• Leading edge gene identification<br>• Multiple testing correction<br>• Custom ranking metrics<br>• Batch processing support | Input: Gene lists/expression + gene sets<br>Output: Enrichment results, plots |
| **PPI analysis** | • STRING-db protein interaction networks<br>• Confidence score filtering<br>• Network visualization<br>• Functional annotation mapping<br>• Cluster detection algorithms<br>• Export to Cytoscape format<br>• Batch gene list processing | Input: Gene/protein lists<br>Output: Interaction networks, annotations |

## 🥅 WGCNA (Weighted Gene Co-expression Network Analysis)

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **WGCNA** | • Scale-free network construction<br>• Soft thresholding power optimization<br>• Module detection algorithms<br>• Trait correlation analysis<br>• Hub gene identification<br>• Both PyWGCNA and R-WGCNA support<br>• Parameter optimization guidance | Input: Expression matrix + traits<br>Output: Modules, networks, correlations |
| **WGCNA network plot** | • Interactive module network visualization<br>• Hierarchical clustering dendrograms<br>• Module-trait relationship heatmaps<br>• Gene connectivity plots<br>• Customizable layout algorithms<br>• Export to multiple formats | Input: WGCNA results<br>Output: Network visualizations |
| **WGCNA hub UMAP** | • UMAP dimensionality reduction<br>• Hub gene highlighting<br>• Module-specific visualizations<br>• Interactive plotting<br>• Gene annotation overlay<br>• Batch processing multiple modules | Input: WGCNA modules + expression<br>Output: UMAP plots, hub gene lists |
| **WGCNA objects comparison** | • Cross-study module preservation<br>• Module similarity metrics<br>• Consensus module identification<br>• Statistical significance testing<br>• Visualization of comparisons<br>• Batch comparison capabilities | Input: Multiple WGCNA objects<br>Output: Preservation statistics, plots |
| **Generate gmt from cluster info** | • Convert cluster results to GMT format<br>• WGCNA module to gene set conversion<br>• Custom gene set creation<br>• Multiple clustering method support<br>• Batch processing<br>• Format validation | Input: Cluster assignments<br>Output: GMT gene set files |

## 🎡 Single-cell RNA-seq (scRNA-seq)

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Pseudobulk** | • Cell type-specific aggregation<br>• Multiple aggregation methods (sum, mean, median)<br>• Minimum cell count filtering<br>• Metadata preservation<br>• Quality control metrics<br>• Batch processing<br>• Export to bulk analysis formats | Input: AnnData (h5ad)<br>Output: Pseudobulk expression matrices |
| **Metacells by SEACells** | • Metacell construction algorithm<br>• Cell similarity-based grouping<br>• Size-balanced metacells<br>• Parameter optimization<br>• Quality assessment metrics<br>• Visualization tools<br>• Integration with downstream analysis | Input: AnnData with embeddings<br>Output: Metacell assignments, objects |
| **Random pseudo-replicates** | • Statistical power enhancement<br>• Random cell sampling strategies<br>• Replicate balance optimization<br>• Bootstrap-based approaches<br>• Multiple splitting methods<br>• Quality control assessment | Input: Single-cell data<br>Output: Pseudo-replicated datasets |
| **memento DE analysis** | • Single-cell differential expression<br>• Multi-condition comparisons<br>• Batch effect handling<br>• Statistical modeling<br>• Multiple testing correction<br>• Effect size estimation | Input: AnnData + conditions<br>Output: DE results, statistics |
| **memento 2D analysis** | • Two-dimensional trajectory analysis<br>• Spatial expression patterns<br>• Gradient-based testing<br>• Pseudotime integration<br>• Multi-modal analysis<br>• Visualization tools | Input: Spatial/trajectory data<br>Output: 2D analysis results |
| **SCENIC heatmap** | • Regulon activity heatmaps<br>• Hierarchical clustering<br>• Cell type annotation<br>• Interactive visualization<br>• Export capabilities<br>• Batch processing | Input: SCENIC regulon activities<br>Output: Activity heatmaps |
| **Prepare regulon data for heatmap** | • Data preprocessing for visualization<br>• Activity score calculation<br>• Normalization options<br>• Filtering parameters<br>• Format conversion<br>• Quality control | Input: Raw SCENIC results<br>Output: Processed regulon data |
| **SCENIC CSI** | • Connection Specificity Index calculation<br>• Regulon connectivity metrics<br>• Statistical significance testing<br>• Comparative analysis<br>• Visualization tools<br>• Batch processing | Input: SCENIC network data<br>Output: CSI scores, statistics |
| **SCENIC network analysis** | • TF-target network visualization<br>• Interactive network plots<br>• Centrality measures<br>• Community detection<br>• Export to Cytoscape<br>• Custom layout algorithms | Input: SCENIC regulons<br>Output: Network visualizations |
| **SCENIC multinetwork analysis** | • Multi-TF network integration<br>• Comparative network analysis<br>• Regulatory cascade identification<br>• Cross-condition comparisons<br>• Advanced visualization<br>• Statistical testing | Input: Multiple regulon sets<br>Output: Integrated network analysis |

## 💬 Cell Communication

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **LIANA LR analysis** | • Multiple ligand-receptor databases<br>• Statistical method integration<br>• Cell type-specific analysis<br>• Confidence scoring<br>• Batch processing<br>• Custom L-R pair support<br>• Visualization tools | Input: AnnData with cell types<br>Output: L-R interaction results |
| **LIANA comparison** | • Cross-condition comparisons<br>• Statistical significance testing<br>• Effect size calculations<br>• Visualization of differences<br>• Multiple comparison correction<br>• Export capabilities | Input: Multiple LIANA results<br>Output: Comparative analysis |
| **CellChat** | • Comprehensive cell communication analysis<br>• Signaling pathway database integration<br>• Network centrality analysis<br>• Pattern recognition<br>• Statistical testing<br>• Python implementation of R CellChat<br>• Batch processing support | Input: AnnData (h5ad)<br>Output: Communication networks, plots |
| **CellChat comparison** | • Two-condition differential analysis<br>• Network topology changes<br>• Pathway-specific comparisons<br>• Statistical significance testing<br>• Visualization of differences<br>• Effect size quantification | Input: Dual-condition h5ad<br>Output: Differential communication |
| **CellChat permutation test** | • Permutation-based significance testing<br>• Multiple testing correction<br>• Bootstrap confidence intervals<br>• Custom test statistics<br>• Batch processing<br>• Result validation | Input: CellChat comparison results<br>Output: Statistical significance |
| **CellChat R qs to python** | • Cross-platform compatibility<br>• SCALA result integration<br>• Format conversion<br>• Data validation<br>• Batch conversion<br>• Quality control | Input: R CellChat .qs files<br>Output: Python-compatible data |


## 🧬 ChIP-seq Analysis

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Sort BAM file** | • Coordinate-based sorting<br>• Memory optimization<br>• Multi-threading support<br>• Quality control metrics<br>• Batch processing<br>• Index generation | Input: Unsorted BAM files<br>Output: Sorted, indexed BAM |
| **Merge BAM files** | • Multiple file merging<br>• Header compatibility checking<br>• Quality score preservation<br>• Memory-efficient processing<br>• Batch operations<br>• Duplicate marking | Input: Multiple BAM files<br>Output: Merged BAM file |
| **Bam to bedGraph for SEACR** | • SEACR-compatible format conversion<br>• Normalization options<br>• Fragment extension<br>• Strand-specific processing<br>• Quality filtering<br>• Batch conversion | Input: BAM files<br>Output: bedGraph files |
| **SEACR peak calling** | • CUT&RUN optimized algorithm<br>• Stringent/relaxed modes<br>• Control sample integration<br>• FDR-based thresholding<br>• Batch processing<br>• Quality metrics | Input: bedGraph + control<br>Output: Peak files |
| **Macs3 peak calling** | • Model-based peak calling<br>• Multiple experimental designs<br>• Statistical significance testing<br>• Multiple output formats<br>• Parameter optimization<br>• Broad/narrow peak modes | Input: BAM/bedGraph files<br>Output: Peak calls, summits |
| **Annotating and filtering peaks** | • Genomic feature annotation<br>• Distance-based filtering<br>• Functional annotation<br>• Custom annotation databases<br>• Statistical enrichment<br>• Batch processing | Input: Peak files (BED)<br>Output: Annotated peaks |
| **Bed length/score filter** | • Size-based filtering<br>• Score threshold application<br>• Distribution visualization<br>• Quality control metrics<br>• Batch processing<br>• Custom criteria support | Input: Peak/BED files<br>Output: Filtered peaks, plots |
| **Bam to DESeq2** | • CUT&RUN specific normalization<br>• Greenlist region counting<br>• DESeq2-compatible output<br>• Quality control integration<br>• Batch processing<br>• Statistical modeling prep | Input: BAM + peak files<br>Output: Count matrices |
| **Blacklist filter** | • Genome blacklist integration<br>• Artifact region removal<br>• Custom blacklist support<br>• Overlap statistics<br>• Batch processing<br>• Quality assessment | Input: Peak/BED files<br>Output: Filtered peak files |
| **Bed to fasta** | • Sequence extraction<br>• Multiple genome support<br>• Custom flanking regions<br>• Batch processing<br>• Quality control<br>• Format validation | Input: BED files + genome<br>Output: FASTA sequences |
| **Denoise bedgraph bigwig** | • Signal smoothing algorithms<br>• Noise reduction techniques<br>• Multiple filtering methods<br>• Quality assessment<br>• Batch processing<br>• Format conversion | Input: bedGraph/bigWig<br>Output: Denoised signal files |

## Ⓜ Miscellaneous Tools

| Tool | Features & Parameters | Input/Output |
|------|---------------------|-------------|
| **Merge excel files** | • Multi-sheet processing<br>• Flexible merging strategies<br>• Header alignment<br>• Data type preservation<br>• Batch processing<br>• Quality control checks | Input: Multiple Excel files<br>Output: Merged Excel/CSV |
| **Remove duplicates** | • Multiple deduplication criteria<br>• Configurable matching rules<br>• Statistics reporting<br>• Memory-efficient processing<br>• Batch operations<br>• Quality assessment | Input: Any tabular data<br>Output: Deduplicated datasets |
| **Split data file on key** | • Key-based file splitting<br>• Custom splitting criteria<br>• Metadata preservation<br>• Batch processing<br>• Quality control<br>• Multiple output formats | Input: Data file + key file<br>Output: Split data files |

## 📜 License

This project is licensed under the MIT License - see the LICENSE file for details.

## 🙏 Acknowledgments

METIS integrates many excellent bioinformatics tools and libraries. We thank all the developers and researchers who created these foundational tools that make METIS possible.