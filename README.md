# DQGlyco
Repository summarising all analysis related to the DQGlyco paper.

## Background

Protein glycosylation is a highly diverse post-translational modification, modulating key cellular processes such as cell signaling, adhesion and cell-cell interactions. Its deregulation has been associated with various pathologies, including cancer and neurological diseases. Methods capable of quantifying glycosylation dynamics are essential to start unraveling the  biological functions of protein glycosylation. Here we present, Deep Quantitative Glycoprofiling method, (DQGlyco), that combines high-throughput sample preparation, high-sensitivity detection, and precise multiplexed quantification of protein glycosylation. We used  DQGlyco to profile the mouse brain glycoproteome, in which we identify 158,972 and 15,056 unique N- and O-glycopeptides localized on 3,199 and 2,365 glycoproteins, respectively, constituting a 25-fold improvement at the glycopeptide level when compared to previous studies.  We observed extensive heterogeneity of glycoforms and determined their functional and structural preferences. The presence of a defined gut microbiota resulted in extensive remodeling of the brain glycoproteome when compared to that of germ-free animals, exemplifying how the gut microbiome may affect brain protein functions.

## Repository structure

The analysis scripts can be found in analysis, which further is divided in structural analysis and data analysis.

### Data analysis

The data analysis starts with psm/protein files which are the output of an MSFragger search (https://msfragger.nesvilab.org0) and can be retrieved from the PRIDE repository of the study (link). In total 9 MSFragger outputfiles are used in this study.

- In *analysis/data_analysis/01_data_processing.Rmd*, these psm files are processed in a standardized way and in *figures_tables_supplements/supp_tables.Rmd* parts of them are used to generate the supplementary tables provided with the paper.
- *analysis/data_analysis/02_enrichment_analyses.Rmd* contains code for all enrichment analyses performed in this study.
- *analysis/data_analysis/03_comparison_to_other_datasets.Rmd* contains code for all comparisons to other studies performed in this study. To run this code, tables from other studies have to be downloaded as indicated in the script. For this, Mus musculus glycosylation data from uniprot has to provided (in data/uniprot).
- *analysis/data_analysis/04_sitecorrelation.Rmd* contains code for all site-centric analyses performed in this study (correlation of glycosylation profiles, comparison to phosphorylation data). For this, Mus musculus phosphorylation data from uniprot has to provided (in data/uniprot).
- *analysis/data_analysis/05_motif_analysis.Rmd* contains code for all motif analyses performed in this study. For this a Mus musculus fasta file has to be provided
- *analysis/data_analysis/06_quantitativeprofiling_microbiomegroups.Rmd* contains code for the analysis of the quantitative glycoproteomics and proteomics data used in this study. For this, Homo sapiens glycosyltransferases data from uniprot has to provided (in data/uniprot).






