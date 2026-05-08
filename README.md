# proj_rapeseed
*Large-scale multi-omics unveils host–microbiome interactions driving root development and nitrogen acquisition*

This repository provides the data and scripts required to reproduce the analyses in the paper.

## Contents
- **Original material list**
- **Cleaned multi-omics datasets** for downstream analyses
- **Kinship matrices** derived from each omics dataset
- **Input data and plotting scripts** for all main and supplementary figures
  - Each folder contains one or two R scripts. By reading the comments at the beginning of each script, users can quickly understand its main purpose. The input files required for running the script can be easily identified from the data-reading commands within the script, allowing users to rerun the analyses and reproduce the initial draft figures. All input files are either located within the same folder as the script or in the 02.CleanOriData folder.
- **R scripts** for multi-omics prediction and association analyses:
  - Genomic prediction of 203 ASVs  
  - Genomic prediction of 13 ionomic  
  - GWAS for 203 ASVs  
  - eGWAS for 17,006 genes  
  - TWAS for all 203 ASVs  

## Reproducibility
All scripts are thoroughly annotated with detailed usage instructions to ensure reproducibility.

## Citation
Li, N., Li, G., Huang, X. *et al.*  **Large-scale multi-omics unveils host–microbiome interactions driving root development and nitrogen acquisition.**  *Nature Plants* **12**, 319–336 (2026).  https://doi.org/10.1038/s41477-025-02210-7
