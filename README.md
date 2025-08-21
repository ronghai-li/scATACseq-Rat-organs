# scATACseq-Rat-organs
Code and pipelines for reproducibility of 'Single-nucleus multiple-organ chromatin accessibility landscape in the adult rat' study

The code includes preprocessing, quality control, integration of scATAC-seq and scRNA-seq data, visualization, and downstream analyses such as motif enrichment, differential accessibility and across species analysis.

# Repository Structure
| File                          | Description                                           |
|-------------------------------|-------------------------------------------------------|
| **01_Functions.R**            | Custom R functions for preprocessing & QC             |
| **02_Preprocessing.R**        | Initial data preprocessing (ArchR input, filtering)   |
| **03_Annotation.R**           | Cell type annotation pipeline                         |
| **0401_Data_Convert.R**       | Data format conversion scripts                        |
| **0402_Plotting.py**          | Python plotting scripts                               |
| **05_Peak_Calling.R**         | Peak calling pipeline                                 |
| **0601_across_organs_analysis.R** | Comparative analysis across rat organs            |
| **0602_Plotting.py**          | Visualization of cross-organ results                  |
| **07_across_species_analysis.R** | Comparative analysis across species                |
| **metadata_AR.tsv**           | Metadata file (sample info, organ)                    |
| **LICENSE**                   | Open-source license (MIT)                             |
| **README.md**                 | Project documentation (this file)                     |
| **.gitignore**                | Ignore unnecessary files                              |

# Requirements

## R Environment
- **R version**: ≥ 4.3.0  

| Package         | Source        |
|-----------------|---------------|
| ArchR           | Bioconductor  |
| patchwork       | CRAN          |
| ggplot2         | CRAN          |
| dplyr           | CRAN          |
| ggrepel         | CRAN          |
| BiocGenerics    | Bioconductor  |
| Seurat          | CRAN/Bioconductor |
| ggpubr          | CRAN          |
| cowplot         | CRAN          |
| tidyr           | CRAN          |
| mclust          | CRAN          |
| gtools          | CRAN          |
| magrittr        | CRAN          |
| tidyverse       | CRAN          |
| magick          | CRAN          |
| circlize        | Bioconductor  |

---

## Python Environment
- **Python version**: ≥ 3.10  

| Package     | Install via |
|-------------|-------------|
| numpy       | pip/conda   |
| pandas      | pip/conda   |
| scanpy      | pip/conda   |
| matplotlib  | pip/conda   |
| seaborn     | pip/conda   |
| plotly      | pip/conda   |
| anndata     | pip/conda   |

# Setup

## R dependencies (install from CRAN/Bioconductor)
```r
install.packages(c("patchwork","ggplot2","dplyr","ggrepel","ggpubr",
                   "cowplot","tidyr","mclust","gtools","magrittr",
                   "tidyverse","magick"))
BiocManager::install(c("ArchR","Seurat","BiocGenerics","circlize"))
```

## Python dependencies (pip or conda)
```bash
pip install numpy pandas scanpy matplotlib seaborn plotly anndata
```

# Input / Output format
Input: raw fragments (tsv.gz), metadata (metadata_AR.tsv)

Output: ArchRProject objects, peak matrices, UMAP/TSNE plots, PDF figures
