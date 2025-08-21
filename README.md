# scATACseq-Rat-organs
Code and pipelines for reproducibility of 'Single-nucleus multiple-organ chromatin accessibility landscape in the adult rat' study

The code includes preprocessing, quality control, integration of scATAC-seq and scRNA-seq data, visualization, and downstream analyses such as motif enrichment, differential accessibility and across species analysis.

.
├── 01_Functions.R              # Custom R functions for preprocessing & QC
├── 02_Preprocessing.R          # Initial data preprocessing (ArchR input, filtering)
├── 03_Annotation.R             # Cell type annotation pipeline
├── 0401_Data_Convert.R         # Data format conversion scripts
├── 0402_Plotting.py            # Python plotting scripts
├── 05_Peak_Calling.R           # Peak calling pipeline
├── 0601_across_organs_analysis.R   # Comparative analysis across rat organs
├── 0602_Plotting.py            # Visualization of cross-organ results
├── 07_across_species_analysis.R    # Comparative analysis across species
├── metadata_AR.tsv             # Metadata file (sample info, organ)
├── LICENSE                     # Open-source license (MIT)
├── README.md                   # This file
└── .gitignore                  # Ignore unnecessary files
