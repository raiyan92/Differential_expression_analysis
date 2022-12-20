#importing libraries
library(tidyverse)
library(forestmangr)
library(DESeq2)
library(pheatmap)
library(viridis)
#fetching read count file
counts = read.csv("D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\raw_counts.csv")
#fetching metadata file
coldata = read.csv("D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\sample1.csv")
#check whether columns in count matches rows in metadata
all(colnames(counts) %in% rownames(coldata))
#creating deseq dataset to perform normalization and achieve DEGs
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ infection)
dds
#pre-filter to remove any genes without avg 1 count per sample # keeping rows that have at least 484 reads total

keep <- rowSums(counts(dds)) >= 484
dds <- dds[keep,]

dds
# set the factor level
dds$infection <- relevel(dds$infection, ref = "negative")
#perorm deseq2
dds <- DESeq(dds, parallel=TRUE, quiet = FALSE) 
res <- results(dds)
#dispersion plot
plotDispEsts(dds)
#MA plot
plotMA(res)
#filtering result by p and padj value
sum(res$pvalue <0.05, na.rm = TRUE)
sum(res$padj <0.1, na.rm = TRUE)
resSig <- res[which(res$padj < 0.1 ), ]  
resSig <- as.data.frame(resSig)
resSig <- rownames_to_column(resSig)
resSig <- dplyr::arrange(resSig, padj)
norm_counts <- as.data.frame(counts(dds, normalized = TRUE)) 
sig_genes <- resSig$rowname
norm_sig <- norm_counts[sig_genes,]

#write csv
write.csv(rsults(dds), "D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\posneg_results_pre-filtered.csv")
write.csv(resSig, "D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\posneg_significant_genes.csv")
write.csv(as.data.frame(counts(dds, normalized = TRUE)), "D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\posneg_norm_counts_pre-filtered.csv")
write.csv(norm_sig, "D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\posneg_sig_gene_counts.csv")
#fetching top 50 genes
#Figure 1A: heatmap of top 50 most significant genes. 
df <- coldata[, c("Age", "infection")]  #add columns as desired to add to metadata. 
colnames(df) <- c("Age", "infection")
de <- (norm_sig[1:50,] + 1)
fc <- de / rowMeans(de)
l2fc <- log2(fc)
#creates dendrogram for samples and genes. if cluster_cols=FALSE it does not cluster
print(pheatmap(l2fc, annotation_col=df, cluster_cols=TRUE))
png(filename="D:\\pendrive\\M.Sc bioinformatics\\Semester 3\\project\\1A_heatmap.png", width=90, height=20, units="in", res=300)