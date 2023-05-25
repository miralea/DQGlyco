#



![Repository summarising all analysis related to the DQGlyco paper.](DQGlyco_logopng.png)


## Background

Protein glycosylation is a highly diverse post-translational modification, modulating key cellular processes such as cell signaling, adhesion and cell-cell interactions. Its deregulation has been associated with various pathologies, including cancer and neurological diseases. Methods capable of quantifying glycosylation dynamics are essential to start unraveling the  biological functions of protein glycosylation. Here we present, Deep Quantitative Glycoprofiling method, (DQGlyco), that combines high-throughput sample preparation, high-sensitivity detection, and precise multiplexed quantification of protein glycosylation. We used  DQGlyco to profile the mouse brain glycoproteome, in which we identify 158,972 and 15,056 unique N- and O-glycopeptides localized on 3,199 and 2,365 glycoproteins, respectively, constituting a 25-fold improvement at the glycopeptide level when compared to previous studies.  We observed extensive heterogeneity of glycoforms and determined their functional and structural preferences. The presence of a defined gut microbiota resulted in extensive remodeling of the brain glycoproteome when compared to that of germ-free animals, exemplifying how the gut microbiome may affect brain protein functions.

---

## Repository structure

The analysis scripts can be found in analysis, which further is divided in structural analysis and data analysis. By using the provided structure in [*data/*](data/) and following the instructions below, all analyses can be repeated.


### Data analysis

The data analysis starts with psm/protein files which are the output of an MSFragger search (https://msfragger.nesvilab.org0) and can be retrieved from the **PRIDE repository** of the study (PXD042237 (Username: reviewer_pxd042237@ebi.ac.uk,  Password: yuVC5zjW)). In total 9 MSFragger outputfiles are used in this study saved to **data/psm_files**.

**In [*analysis/data_analysis/01_data_processing.Rmd*](analysis/data_analysis/01_data_processing.Rmd), these psm files are processed and saved in a standardised way ins data/processed.** 

---

All further analyses start from the processed files saved to **data/processed**.

- [*analysis/data_analysis/02_enrichment_analyses.Rmd*](analysis/data_analysis/02_enrichment_analyses.Rmd) contains code for all enrichment analyses performed in this study.

- [*analysis/data_analysis/03_comparison_to_other_datasets.Rmd*](analysis/data_analysis/03_comparison_to_other_datasets.Rmd) contains code for all comparisons to other studies performed in this study. Tables from other studies have to be downloaded as indicated in the script (in data/comparison_to_datasets). For this, Mus musculus glycosylation data from uniprot has to be provided (in data/uniprot).

- [*analysis/data_analysis/04_sitecorrelation.Rmd*](analysis/data_analysis/04_sitecorrelation.Rmd) contains code for all site-centric analyses performed in this study (correlation of glycosylation profiles, comparison to phosphorylation data). For this, Mus musculus phosphorylation data from uniprot has to be provided (in data/uniprot).

- [*analysis/data_analysis/05_motif_analysis.Rmd*](analysis/data_analysis/05_motif_analysis.Rmd) contains code for all motif analyses performed in this study. For this a Mus musculus fasta file has to be provided.

- [*analysis/data_analysis/06_quantitativeprofiling_microbiomegroups.Rmd*](analysis/data_analysis/06_quantitativeprofiling_microbiomegroups.Rmd) contains code for the analysis of the quantitative glycoproteomics and proteomics data used in this study. For this, Homo sapiens glycosyltransferases data from uniprot has to be provided (in data/uniprot).

### Figures and Tables

In [*analysis/data_analysis/figures_tables_supplements*](analysis/data_analysis/figures_tables_supplements) code used to generate all figures, tables and supplementary data is deposited. Code either starts from the processed results in **data/processed** or results of the analyses described in [*analysis/data_analysis*](analysis/data_analysis). Please reach out in case of questions. Please note that figure colour, legends, axis title etc were modified in Inkscape to be included in the final manuscripts, so slight differences arise from that.





