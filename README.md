#



![Repository summarising all analysis related to the DQGlyco paper.](DQGlyco_logopng.png)


## Background

Protein glycosylation is a highly diverse post-translational modification, modulating key cellular processes such as cell signaling, adhesion and cell-cell interactions. Its deregulation has been associated with various pathologies, including cancer and neurological diseases. Methods capable of quantifying glycosylation dynamics are essential to start unraveling the biological functions of protein glycosylation. Here we present Deep Quantitative Glycoprofiling (DQGlyco), a method that combines high-throughput sample preparation, high-sensitivity detection, and precise multiplexed quantification of protein glycosylation. We used  DQGlyco to profile the mouse brain glycoproteome, in which we identify close to 180,000 unique N-glycopeptides. We observed extensive heterogeneity of glycoforms and determined their functional and structural preferences. We used our quantitative approach to characterize glycosites tissue-specificity and demonstrated that the presence of a defined gut microbiota resulted in extensive remodeling of the brain glycoproteome when compared to that of germ-free animals, exemplifying how the gut microbiome may affect brain protein functions. All results and data are available at https://apps.embl.de/glycoapp/.

---

## Repository structure

The analysis scripts can be found in analysis, which further is divided in structural analysis and data analysis. 


### 1. Data analysis

The data analysis starts with psm/protein files which are the output of an MSFragger search (https://msfragger.nesvilab.org0) and can be retrieved from the **PRIDE repository** of the study (PXD042237 (Username: reviewer_pxd042237@ebi.ac.uk,  Password: yuVC5zjW)). 

**In [*analysis/data_analysis/01_data_processing.Rmd*](analysis/data_analysis/01_data_processing.Rmd), these psm files are processed and saved in a standardised way.** 

---

All further analyses start from the processed files and are done per dataset (eg labelfree, tissue data, solubility data etc).

---


### Figures and Tables

In [*analysis/data_analysis/figures_tables_supplements*](analysis/data_analysis/figures_tables_supplements) code used to generate all figures, tables and supplementary data is deposited. Code either starts from the processed results or results of the analyses described in [*analysis/data_analysis*](analysis/data_analysis). Please reach out in case of questions. Please note that figure colour, legends, axis title etc were modified in Inkscape to be included in the final manuscripts, so slight differences arise from that.

### 2. Structural analysis

The structural analysis code can be found in [*analysis/structural_analysis*](analysis/structural_analysis) and has a separate README with additional explanation.





