# PATH40040 Final Project Code

# Morgan Gallagher, 23205638, 20/12/2023

library(DESeq2)
library(BiocManager)
library(clusterProfiler)
library(dplyr)
library(ggrepel)

#read in data ------------------------------------------------------------------

#set file path
path <- "/Users/morgan/Desktop/UCD/Bio Princ/Final Project/"
file_name <- "brca_tcga_pan_can_atlas_2018.tar.gz"

#untar folder
untar(paste0(path, file_name))

# change directory 
setwd(paste(getwd() , "/brca_tcga_pan_can_atlas_2018", sep = ""))

#read patient info data
clinical = read.delim("data_clinical_patient.txt")

#read rnaseq data
rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

#read mutation data
data_mutations = read.delim("data_mutations.txt")

#read copy number data
cna = read.delim('data_cna.txt')

#data preparation --------------------------------------------------------------

#create vector if row is pseudogene or not
keep = !duplicated(rnaseq[,1])

#first gene is pseudogene, change to FALSE
keep[1] = keep[1] == FALSE

#keep TRUE rows in rnaseq
rnaseq = rnaseq[keep,]

# set rownames of rnaseq to hugo symbols
rownames(rnaseq)  = rnaseq[,1]

# find ERBB2 row in cna data
erbb2_indx = which(cna[,1] == 'ERBB2')

# match patients in rnaseq to patients in cna.
rna_cna_id = which(is.element(colnames(rnaseq[,-c(1,2)]), 
                              colnames(cna[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = 
  sum(is.element(colnames(rnaseq[,2+rna_cna_id]), colnames(cna[,-c(1,2)]))) 

#create empty metadata matrix
meta_erbb2 = matrix(0,length(rna_cna_id),1)

#indicate of each patient has amplified erbb2 or not
for (i in 1:length(rna_cna_id)){
  # access the colnames of i, patient ids
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(cna[erbb2_indx,col_cna]>0)
  
}

# rename column in metadata
colnames(meta_erbb2) = 'ERBB2Amp'

# transform rnaseq data to int for analysis
rna_cna_sub = round(rna_cna_sub)

#impute zeros for NA and negative values
rna_cna_sub[is.na(rna_cna_sub)] = 0  
rna_cna_sub[rna_cna_sub<0] = 0

#differential expression analysis ----------------------------------------------

#create DESSeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = rna_cna_sub,
                                  colData = meta_erbb2,
                                  design = ~ ERBB2Amp)

#set threshold for group (number of people with more than 10 counts)
smallestGroupSize <- 3

#only keep genes where 3 or more people have true values
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize

#only keep the keep genes in the DESeqDataSet
dds <- dds[keep,]

#Normalize data
dds <- DESeq(dds)

#store results of diff ex analysis
res <- results(dds)

#summary of results
summary(res)

#rename rows as genes
rownames(res) = rnaseq[keep,1]

#get indices of which adjusted p values are less than 0.05, stat sig
signif = which(res$padj<0.05)

#store results that are only statistically significant
deg = res[signif,]

summary(deg)

#sort genes by log2 fold
sorted_res <- deg[order(abs(deg$log2FoldChange), 
                                decreasing = TRUE), ]

#store top 10 genes by their log2 fold
top_10_genes <- head(sorted_res, 10)

#print top 10 diff ex genes
print(top_10_genes)

#upreg and downreg genes
upregulated_genes = sum(deg$log2FoldChange > 0)
downregulated_genes = sum(deg$log2FoldChange < 0)

cat("Number of upregulated genes:", upregulated_genes, "\n")
cat("Number of downregulated genes:", downregulated_genes, "\n")

# Convert DESeqResults object to data frame
res_df <- as.data.frame(res)

# Identify the top 10 differentially expressed genes with padj < 0.05
top_genes <- head(res_df[order(-abs(res_df$log2FoldChange)), ], 10)

# Create MA plot with labels for top 10 genes with padj < 0.05
ma_plot <- ggplot(res_df, aes(x = log2(baseMean), y = log2FoldChange)) +
  geom_point(aes(color = ifelse(padj < 0.05, "Significant", "Not Significant")), alpha = 0.7) +
  scale_color_manual(values = c("black", "red"), guide = guide_legend(title = NULL)) +  
  labs(title = "MA Plot",
       x = "Mean Expression (log2)",
       y = "Log2 Fold Change") +
        theme_minimal()

# Display the MA plot
print(ma_plot)

#pca ---------------------------------------------------------------------------

# transform data using vst
rld <- vst(dds, blind = FALSE)

# perform pca
pc <- prcomp(assay(rld))

# create a data frame for plotting
pc_data <- as.data.frame(pc$rotation[, c(1, 2)])
pc_data$Group <- factor(meta_erbb2)

# plot pca
ggplot(pc_data, aes(PC1, PC2, color = `Group`)) +
  geom_point(size = 2) +
  labs(title = "PCA Plot",
       x = "Principal Component 1",
       y = "Principal Component 2",
       color = "Group") +
  theme_minimal()

  
#pathway enrichment analysis ---------------------------------------------------

#obtain the entrez id's
entrez_ids = rnaseq[keep,2]

#id's that were significant in diff seq
entrez_all = entrez_ids[signif]

#id's that were upregulated
entrez_up = entrez_all[signif[deg[,2]>0.]]

#id's that were downregulated
entrez_down = entrez_all[signif[deg[,2]<0.]]

#all pathways
all_paths =   enrichKEGG(gene = entrez_all, 
                             organism = 'hsa', pvalueCutoff = 0.05)


sorted_paths <- all_paths[order(-all_paths$Count), ]
top_pathways <- head(sorted_paths, 5)
print(top_pathways)

#upregulated pathways
up_paths =   enrichKEGG(gene = entrez_up, 
                            organism = 'hsa', pvalueCutoff = 0.05)

sorted_up_paths <- up_paths[order(-up_paths$Count), ]
head(sorted_up_paths)

head(up_paths)

#downregulated pathways
down_paths =   enrichKEGG(gene = entrez_down, 
                              organism = 'hsa', pvalueCutoff = 0.05)

sorted_down_paths <- down_paths[order(-down_paths$Count), ]
head(sorted_down_paths)

head(down_paths)
