# Structural Analysis

## Dependencies

### Python

The following python modules are needed to reproduce the analysis

```
-----
numpy               1.22.4
pandas              1.5.3
session_info        1.0.0
structuremap        0.0.8
tqdm                4.65.0
-----
Python 3.9.12 (main, Apr  5 2022, 01:53:17) [Clang 12.0.0 ]
macOS-10.16-x86_64-i386-64bit
-----
```

### R

The following R packages are needed to reproduce the analysis

```
-----
 [1] ggforce_0.3.4        openxlsx_4.2.5       cowplot_1.1.1        pbapply_1.7-0       
 [5] ggpubr_0.4.0         forcats_0.5.1        stringr_1.4.0        dplyr_1.0.9         
 [9] purrr_0.3.4          readr_2.1.2          tidyr_1.2.0          tibble_3.1.8        
[13] ggplot2_3.3.6        tidyverse_1.3.2      seqinr_4.2-16        here_1.0.1          
[17] org.Hs.eg.db_3.15.0  org.Mm.eg.db_3.15.0  AnnotationDbi_1.58.0 IRanges_2.30.0      
[21] S4Vectors_0.34.0     Biobase_2.56.0       BiocGenerics_0.42.0 
-----
R version 4.2.1 (2022-06-23)
Platform: x86_64-apple-darwin17.0 (64-bit)
Running under: macOS Ventura 13.4
-----
```

### Summary

The jupyter notebooks `1_download_topological_domains.ipynb`, `2_retrieve_annotated_alphafold.ipynb` and `3_preprocess_sift_scores.ipynb` are used to generate the needed inputs for the `main.Rmd` notebook, which uses all the structural, domain and conservation information to explore the properties of glycosites with different levels of microheterogeneity. The `main.Rmd` uses the preprocessed PTM tables provided in the manuscript as Supplementary Files.

To execute `3_preprocess_sift_scores.ipynb`, precomputed SIFT scores for a given organism need to be retrieved from https://sift.bii.a-star.edu.sg/sift4g/ and collapsed to a single file. For  `1_download_topological_domains.ipynb` and `2_retrieve_annotated_alphafold.ipynb` the only requirement is a stable internet connection. 


