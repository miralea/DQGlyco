#



![Repository summarising all analysis related to the DQGlyco paper.](DQGlyco_logopng.png)


## Background

Protein glycosylation, a diverse post-translational modification, regulates key cellular processes such as signaling, adhesion, and cell-cell interactions. Its dysregulation is linked to diseases like cancer and neurological disorders. To explore the biological roles of glycosylation, we present Deep Quantitative Glycoprofiling (DQGlyco), a method combining high-throughput sample preparation, sensitive detection, and precise multiplexed quantification. Using DQGlyco, we profiled the mouse brain glycoproteome, identifying nearly 180,000 unique N-glycopeptides — 25 times more than previous studies. This revealed extensive glycoform heterogeneity and their structural preferences. We applied DQGlyco to quantify glycopeptide regulation in human cells treated with a fucosylation inhibitor and to study tissue-specific glycosylation patterns in mice. Additionally, we showed that the presence of a defined gut microbiota induces significant remodeling of the mouse brain glycoproteome, exemplifying how the gut microbiome may affect brain protein functions. Finally, we introduce a strategy for the systematic characterization of glycoforms biophysical properties, offering insights into glycosylation functionality. Overall, DQGlyco’s in-depth profiling uncovered previously unappreciated complexity in glycosylation regulation. All results and data are available at https://apps.embl.de/glycoapp/.

---

## Repository structure

The analysis scripts can be found in analysis, which further is divided in structural analysis and data analysis. 


### 1. Data analysis

The data analysis starts with psm/protein files which are the output of an MSFragger search (https://msfragger.nesvilab.org0) and can be retrieved from the **PRIDE repository** of the study (PXD042237 (Username: reviewer_pxd042237@ebi.ac.uk,  Password: yuVC5zjW)). 

**In [*analysis/data_analysis/01_data_processing.Rmd*](analysis/data_analysis/01_data_processing.Rmd), these label free psm files are processed and saved in a standardised way. The psm files of the TMT labelled data is processed the same way in the respective scripts and result objects used for analysis and plotting are provided in data.7z (https://oc.embl.de/index.php/s/FiJI0NKFgVB9uDM). ** 

---

All further analyses start from the processed files and are done per dataset (eg labelfree, tissue data, solubility data etc). All results, external datasets etc are deposited in data.7z (https://oc.embl.de/index.php/s/FiJI0NKFgVB9uDM) and result files are additionally available at *https://apps.embl.de/glycoapp/* in the download section or as supplementary material of the manuscript.

---


### Figures and Tables

In [*analysis/data_analysis/figures_tables_supplements*](analysis/data_analysis/figures_tables_supplements) code used to generate all figures, tables and supplementary data is deposited. Code either starts from the processed results or results of the analyses described in [*analysis/data_analysis*](analysis/data_analysis). Please reach out in case of questions (mira.burtscher(at)embl.de). Please note that figure colour, legends, axis title etc were modified in Inkscape to be included in the final manuscripts, so slight differences arise from that.

### 2. Structural analysis

The structural analysis code can be found in [*analysis/structural_analysis*](analysis/structural_analysis) and has a separate README with additional explanation.





