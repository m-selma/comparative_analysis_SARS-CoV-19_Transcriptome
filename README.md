# Comparative Analysis of SARS-CoV-19 Transcriptome: classification approaches using features extracted from single-cell vs. bulk RNA-seq data

## Problem statement and strategy
SARS-CoV-19 is a multi-system disease characterized by an interplay between immunological and inflammatory cascades. Attempts to describe the heterogeneous host response to the virus will enable a precision approach to therapy. 

Aim: to leverage transcriptomic data analysis and supervised ML techniques to understand the whole-blood transcriptomic host response to SARS-CoV-2. 

### Part 1: scRNA-seq data
- Feature Extraction: use scRNA-seq data from SARS-CoV-19 patients vs. healthy controls to identify differentially expressed genes (DEGs) between the groups. 
- Classification: use supervised binary classification to differentiate between SARS-CoV-19 patients and healthy controls in bulk datasets using identified DEGs as features.
- Gene Ontology Analysis: identify biological pathways associated with the identified SARS-CoV-19 biomarker genes.

### Part 2: Bulk Data
- Feature Extraction: use bulk RNA-seq data from SARS-CoV-19 patients vs. healthy controls to identify DEGs. 
- Classification: repeat classification process on the same datasets with DEGs from bulk data as features.

--> compare classification performance using features selected from sc vs. bulk data.


## Datasets
All datasets used are publicly available.
### For part 1:
- [feature extraction] - scRNA-seq: GSE166992, gene expression data of PBMN cells of 5 SARS-CoV-19 patients and 3 healthy controls.
- [classification] - bulk RNA-seq: GSE152641 (62 SARS-CoV-19/24 healthy), GSE152075 (430 SARS-CoV-19/54 healthy), GSE152418 (17 SARS-CoV-19/17 healthy).
### For part 2:
- [feature extraction] - bulk RNA-seq: GSE196822 (37 SARS-CoV-19/9 healthy).
- [classification] - bulk RNA-seq: GSE152641 (62 SARS-CoV-19/24 healthy), GSE152075 (430 SARS-CoV-19/54 healthy), GSE152418 (17 SARS-CoV-19/17 healthy).


