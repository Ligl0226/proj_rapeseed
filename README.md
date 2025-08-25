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
> Li N, Li G, Wang D, Ma L, Huang X, Bai Z, et al.  
> **Large-scale multi-omics analyses identified root-microbiome associations underlying plant nitrogen nutrition** [Internet].  
> 2024 [cited 2025 Aug 25]. Available from: [http://biorxiv.org/lookup/doi/10.1101/2024.02.05.578621](http://biorxiv.org/lookup/doi/10.1101/2024.02.05.578621)
