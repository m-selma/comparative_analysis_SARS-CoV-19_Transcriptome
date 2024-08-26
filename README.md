# comparative_analysis_SARS-CoV-19_Transcriptome

Problem statement: SARS-CoV-19 is a multi-systems disease characterized by an interplay between immunological and inflammatory cascades. Attempts to describe the heterogeneous host response to the virus will enable a precision approach to therapy. 

Aim: to leverage transcriptomic data analysis and supervised ML techniques to understand the whole-blood transcriptomic host response to SARS-CoV-2. 

### Part 1: scRNA-seq data
- Feature Extraction: use scRNA-seq data from SARS-CoV-19 patients vs. healthy controls to identify differentially expressed genes (DEGs) between the groups. 
- Classification: use supervised binary classification to differentiate between SARS-CoV-19 patients and healthy controls in bulk datasets using identified DEGs as features.
- Gene Ontology Analysis: identify biological pathways associated with the identified SARS-CoV-19 biomarker genes.

### Part 2: Bulk Data
- Feature Extraction: use bulk RNA-seq data from SARS-CoV-19 patients vs. healthy controls to identify DEGs. 
- Classification: repeat classification process on the same datasets with DEGs from bulk data as features.

--> compare classification performance using features selected from sc vs. bulk data.


### Datasets
The datasets used are publicly available.
### For part 1:
- [feature extraction] - scRNA-seq: gene expression data of PBMN cells of 5 SARS-CoV-19 patients and 3 healthy controls (see /sc_data).
- [classification] - bulk RNA-seq: GSE152641 (62 SARS-CoV-19/24 healthy), GSE152075 (430 SARS-CoV-19/54 healthy), GSE152418 (17 SARS-CoV-19/17 healthy).
### For part 2:
- [feature extraction] - bulk RNA-seq: GSE196822 (37 SARS-CoV-19/9 healthy).
- [classification] - bulk RNA-seq: GSE152641 (62 SARS-CoV-19/24 healthy), GSE152075 (430 SARS-CoV-19/54 healthy), GSE152418 (17 SARS-CoV-19/17 healthy).

1. GSE152075: RNA-sequencing profiles of nasopharyngeal swabs from 430 individuals with SARS-CoV-2 and 54 negative controls. (Platform: Illumina NextSeq 500)
2. GSE152418: RNAseq analysis of PBMCs in a group of 17 COVID-19 subjects and 17 healthy controls. (Platform: Illumina NovaSeq 6000)

