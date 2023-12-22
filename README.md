# PATH40040
Final Project

The code can be read in four parts, these can be seen as sectioned off in the R script. Data are taken from https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018. \
*Comments throughout the script can serve as a guide. 

Section 1: \
The first section involves reading in the data. The following documents are used in our analysis: data_clinical_patient.txt, data_mrna_seq_v2_rsem.txt, and data_mutations.txt. 

Section 2: \
The second section of the code prepares and filters the data. 
We format the data so that it is in an acceptable format to become a DESeq data set object. We locate the ERBB2 gene from the cna dataset and match this to the rnaseq data from patients. We then create a vector which contains the metadata indicating if a given patient has ERBB2 positive or negative breast cancer. 

Section 3: \
The third section performs differential expression analysis using the DESeq2 package. 

Section 4: \
The fourth section performs principal component analysis on the count data 

Section 5: \
The fifth section performs pathway enrichment using the enrichKEGG function. 
