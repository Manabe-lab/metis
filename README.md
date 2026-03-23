<div align="center">
<img src="Metis_logo_github.png" alt="METIS Logo" width="200">
</div>

## **M**olecular **E**xploration and **T**ranscriptomic **I**nvestigation **S**uite

METIS is a comprehensive bioinformatics web application suite built with [Streamlit](https://streamlit.io/) for RNA-seq, single-cell RNA-seq, ChIP-seq/CUT&RUN/ATAC-seq analyses, and various data manipulation tasks. It provides an intuitive interface for researchers to perform complex genomic analyses without requiring extensive programming knowledge.

## Quick Start

### Installation

1. Clone this repository:
```bash
git clone https://github.com/Manabe-lab/metis.git
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

## Available Applications (125 tools in 17 categories)

### NGS Pipeline (3 apps)
| Tool | Description |
|------|-------------|
| FASTQ to RNA-seq count table | Automated pipeline from FASTQ to gene count matrix |
| FASTQ to ChIP/CUT&RUN/ATAC BAM | Alignment pipeline for epigenomic experiments |
| FASTQ to scRNA (Cell Ranger) | Single-cell RNA-seq preprocessing with Cell Ranger |

### Data Normalization / Manipulation (9 apps)
| Tool | Description |
|------|-------------|
| Count Data Normalization | RNA-seq count normalization with multiple methods |
| Count to TPM | Convert raw counts to TPM normalization |
| Homer or DESeq2 to Data | Extract clean data from Homer/DESeq2 results |
| Manipulate Data Table | Data table manipulation and restructuring |
| Merge Data Files | Combine multiple data files |
| Filter and Transform Data | Data filtering and mathematical transformations |
| Spread Sheet | Interactive data editing interface |
| Batch removal by Combat-seq | sva Combat-seq batch correction for count data |
| Filter X/Y Chromosome Genes | Remove sex chromosome genes from expression data |

### Gene Name Conversion (6 apps)
| Tool | Description |
|------|-------------|
| Update Gene Symbol | Update gene symbols to latest HGNC/MGI standards |
| Gene Symbol to Ensembl ID | Convert gene symbols to Ensembl IDs |
| Ensembl ID to Gene Symbol | Convert Ensembl IDs to gene symbols |
| Mouse Human Symbol Conversion | Cross-species ortholog mapping |
| GMT Mouse Human Conversion | Convert gene sets between mouse and human |
| Microarray Gene Name Filter | Filter and map microarray probe names to gene symbols |

### Data Analysis (12 apps)
| Tool | Description |
|------|-------------|
| DESeq2 | Differential expression analysis with Wald test and GLM |
| DESeq2-LRT | Time course and multi-group analysis with Likelihood Ratio Test |
| edgeR | Differential expression with exact test and quasi-likelihood |
| limma | Linear modeling for microarray and RNA-seq (voom) |
| DE method comparison | Compare multiple DE methods on the same dataset |
| Permutation test | Non-parametric permutation-based differential testing |
| Make ranking file for GSEA | Create GSEA-compatible ranking files from DE results |
| impulseDE2 | Impulse model for time-course differential expression |
| Mfuzz clustering | Soft clustering of time-course expression profiles |
| DEA result comparison | Compare and integrate differential expression results |
| Compare ratios | Compare cell-type or condition ratios |
| KS Test Distribution Analysis | Kolmogorov-Smirnov test for distribution comparisons |

### Data Visualization (6 apps)
| Tool | Description |
|------|-------------|
| PCA | PCA, UMAP, tSNE, MDS dimensionality reduction |
| PCA statistics | Statistical analysis of PCA components and covariates |
| Volcano plot | Interactive volcano plots for DE results |
| Heatmap | Hierarchical clustering heatmaps |
| Box/Violin plot | Box, violin, and strip plots with statistics |
| Venn/Upset Plot | Multi-set intersection analysis |

### Pathway Analysis (5 apps)
| Tool | Description |
|------|-------------|
| decoupleR | Pathway and TF activity inference |
| GSEApy | Gene Set Enrichment Analysis |
| aPEAR | Enrichment network visualization with LLM-powered cluster naming |
| PADOG | Pathway Analysis with Down-weighting of Overlapping Genes |
| PPI analysis | STRING-db protein-protein interaction networks |

### WGCNA (5 apps)
| Tool | Description |
|------|-------------|
| WGCNA | PyWGCNA/R-WGCNA co-expression network analysis |
| WGCNA network plot | Module network visualization |
| WGCNA hub UMAP | UMAP visualization of module genes |
| WGCNA objects comparison | Compare modules between conditions/studies |
| Generate gmt from cluster info | Convert clustering results to GMT gene sets |

### scRNA-seq Manipulation (8 apps)
| Tool | Description |
|------|-------------|
| Seurat to h5ad conversion | Convert between Seurat and AnnData formats |
| Dimension reduction | PCA, UMAP, tSNE for single-cell data |
| Clustering | Graph-based clustering (Leiden, Louvain) |
| Batch integration | Harmony, scVI, BBKNN integration methods |
| Subset / Merge h5ad | Subset or merge AnnData objects |
| Cell cycle scoring | Cell cycle phase assignment |
| Edit metadata | Modify cell metadata and annotations |
| CellTypist annotation | Automated cell type annotation |

### scRNA-seq Cell Aggregation (4 apps)
| Tool | Description |
|------|-------------|
| Pseudobulk | Create pseudobulk data from single-cell AnnData |
| Metacells by SEACells | SEACells metacell construction |
| Random pseudo-replicates | Create pseudo-replicates by random cell splitting |
| Pseudobulk DEA (RUV) | Pseudobulk differential expression with RUV normalization |

### scRNA-seq Analysis (8 apps)
| Tool | Description |
|------|-------------|
| Feature & Expression plots | UMAP, violin, dot plots for gene expression |
| Marker gene analysis | Find and visualize marker genes per cluster |
| memento DE analysis | Memento single-cell differential expression |
| memento 2D analysis | Memento 2D trajectory analysis |
| scCODA compositional analysis | Compositional analysis of cell-type proportions |
| Milo differential abundance | Differential abundance testing with neighborhoods |
| CINEMA-OT perturbation analysis | Optimal transport-based perturbation analysis |
| Geneformer perturbation | Transformer-based in silico perturbation |

### scRNA-seq Pathway (4 apps)
| Tool | Description |
|------|-------------|
| COMPASS metabolic analysis | Single-cell metabolic flux analysis |
| COMPASS postprocess | Visualize and analyze COMPASS results |
| scFEA metabolic analysis | scFEA flux estimation |
| scFEA postprocess | Visualize and analyze scFEA results |

### SCENIC (6 apps)
| Tool | Description |
|------|-------------|
| SCENIC heatmap | Visualize regulon activity patterns |
| Prepare regulon data for heatmap | Process raw SCENIC results |
| SCENIC CSI | Connection Specificity Index for regulons |
| SCENIC network analysis | Gene regulatory network visualization |
| SCENIC multinetwork analysis | Multi-TF regulatory network integration |
| SCENIC network analysis 2 | Advanced network with TRRUST/ChIP-Atlas/STRING integration |

### RNA Velocity (16 apps)
| Tool | Description |
|------|-------------|
| Data filtering | Pre-filtering for velocity analysis |
| scVelo analysis | RNA velocity estimation with scVelo |
| scVelo visualization | Velocity stream and embedding plots |
| CellRank analysis | Fate probability and lineage analysis |
| CellRank visualization | Visualize CellRank results |
| DeepVelo analysis | Deep learning-based velocity estimation |
| UniTVelo analysis | Unified time-informed velocity |
| noSpliceVelo analysis | Velocity without splicing information |
| Pseudotime gene expression | Gene expression along pseudotime |
| Dynamo analysis | Differential geometry-based velocity analysis |
| Dynamo visualization | Visualize Dynamo results |
| Dynamo perturbation | In silico perturbation with vector fields |
| Dynamo LAP | Least Action Path analysis |
| TFvelo analysis | Transcription factor-informed velocity |
| TFvelo to Dynamo | Convert TFvelo results to Dynamo format |
| VeloViz | Velocity-informed 2D embedding |

### Cell Communication (6 apps)
| Tool | Description |
|------|-------------|
| LIANA LR analysis | LIANA+ ligand-receptor interaction analysis |
| LIANA comparison | Compare LIANA results between conditions |
| CellChat | Python CellChat cell-cell communication analysis |
| CellChat comparison | Compare CellChat results between conditions |
| CellChat permutation test | Statistical testing of CellChat comparisons |
| CellChat R qs to python | Convert R CellChat results to Python format |

### BAM / BED Operation (4 apps)
| Tool | Description |
|------|-------------|
| Sort BAM file | SAMtools coordinate-based BAM sorting |
| Merge BAM files | SAMtools BAM file merging |
| BAM to bigWig/bedgraph | Convert BAM to bigWig or bedGraph format |
| Bed to fasta | Extract sequences from BED coordinates |

### Peak Calling (14 apps)
| Tool | Description |
|------|-------------|
| SEACR peak caller | CUT&RUN sparse enrichment peak calling |
| GoPeaks / Summit | GoPeaks broad and narrow peak calling |
| LanceOtron peak caller | Deep learning-based peak calling |
| Macs3 peak caller | Model-based peak calling |
| Homer peak caller | Homer findPeaks peak calling |
| Omnipeak peak caller | Omnipeak consensus peak calling |
| Blacklist filter | Remove peaks in blacklist regions |
| Annotate/Filter peaks | Peak annotation and genomic feature filtering |
| Bed length/score filter | Filter peaks by length and score |
| Bed peak merger | Merge overlapping peak regions |
| Peak comparison | Compare peak sets between samples |
| BAM QC Dashboard | BAM file quality metrics dashboard |
| IDR Analysis | Irreproducibility Discovery Rate analysis |
| TOBIAS Footprinting | Transcription factor footprinting |

### Downstream Analysis (8 apps)
| Tool | Description |
|------|-------------|
| Peak counter for histogram/heatmap | Homer quantification for visualization |
| Heatmap for genome occupancy | Genome-wide signal heatmaps |
| Peak counter for differential analysis | Count reads in peaks for DE |
| Differential binding (no replicate) | Differential binding without replicates |
| Differential binding (with replicates) | Differential binding with biological replicates |
| Homer Motif Analysis | De novo and known motif analysis |
| STREME Motif Analysis | MEME Suite motif discovery |
| ChromHMM Chromatin State | Chromatin state segmentation |

## Required Packages

### Core Dependencies
- **streamlit** - Web application framework
- **pandas**, **numpy** - Data manipulation and numerical computing
- **matplotlib**, **seaborn**, **plotly** - Visualization
- **scipy**, **scikit-learn** - Scientific computing and machine learning

### Bioinformatics Libraries
- **scanpy**, **anndata** - Single-cell analysis
- **scvelo**, **dynamo-release**, **cellrank** - RNA velocity
- **pyscenic** - Gene regulatory network analysis
- **decoupler**, **GSEApy** - Pathway activity inference
- **PyWGCNA** - Weighted Gene Co-expression Network Analysis
- **liana** - Ligand-receptor analysis
- **cellchat** - Cell communication analysis

### R Integration
- **rpy2** - R interface for Python
- **R packages**: DESeq2, edgeR, limma, WGCNA, Combat-seq, impulseDE2, Mfuzz, PADOG

### File Processing
- **openpyxl** - Excel file handling
- **h5py** - HDF5 file format
- **pysam** - SAM/BAM file processing
- **pybedtools** - BED file manipulation

### Specialized Tools
- **pyscenic** - SCENIC gene regulatory networks
- **compass** - Metabolic flux analysis
- **celltypist** - Automated cell type annotation
- **geneformer** - Transformer-based perturbation

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

METIS integrates many excellent bioinformatics tools and libraries. We thank all the developers and researchers who created these foundational tools that make METIS possible.
