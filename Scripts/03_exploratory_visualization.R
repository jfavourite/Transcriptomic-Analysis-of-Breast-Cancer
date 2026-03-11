# Module 03: Data Visualization & Quality Control
# Purpose: PCA, Volcano, and Heatmap generation

library(DESeq2)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)

dds <- readRDS("data/dds_analyzed.rds")
res <- readRDS("data/res_raw.rds")

# 1. Variance Stabilizing Transformation (VST) for PCA/Heatmap
vsd <- vst(dds, blind = FALSE)

# 2. PCA Plot (Checking for batch effects/clustering)
pca_plot <- plotPCA(vsd, intgroup="condition") + 
  theme_minimal() + 
  ggtitle("PCA: Clustering of Breast Cancer Tumor vs Normal Samples")
ggsave("results/figures/PCA_plot.png", pca_plot)

# 3. Volcano Plot Construction
res_df <- as.data.frame(res)
res_df$status <- "Not Significant"
res_df$status[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated (Tumor)"
res_df$status[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated (Normal)"

volcano <- ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=status)) +
  geom_point(alpha=0.4, size=1.5) +
  scale_color_manual(values=c("blue", "grey", "red")) +
  geom_vline(xintercept=c(-1, 1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_classic() +
  labs(title="Volcano Plot: Breast Cancer Somatic Expression Profile")
ggsave("results/figures/Volcano_plot.png", volcano)

# 4. Heatmap of Top 50 Genes
top50 <- head(order(res$padj), 50)

vsd_mat <- assay(vsd)

top_genes <- head(order(res$padj), 50)
plot_mat <- vsd_mat[top_genes, ]

plot_mat <- t(apply(plot_mat, 1, scale))
colnames(plot_mat) <- colnames(vsd)

annotation_col <- data.frame(Group = colData(dds)$condition)
rownames(annotation_col) <- colnames(vsd)


ensembl_ids <- gsub("\\..*|\\|.*", "", rownames(plot_mat))

gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys = ensembl_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")
gene_symbols[is.na(gene_symbols)] <- ensembl_ids[is.na(gene_symbols)]


pheatmap(plot_mat, 
         main = "TUMOR vs Adjacent Normal of WGS BREAST CANCER samples: Top 50 Genes",
         labels_row = gene_symbols, # Use the new symbols here
         annotation_col = annotation_col,
         show_colnames = TRUE, 
         show_rownames = TRUE, 
         cluster_cols = FALSE,
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(150),
         fontsize_row = 10)


#Saving the top 50 genes in a csv file
res_final <- as.data.frame(res)

res_final$symbol <- mapIds(org.Hs.eg.db, 
                           keys = gsub("\\..*", "", rownames(res_final)), 
                           column = "SYMBOL", 
                           keytype = "ENSEMBL", 
                           multiVals = "first")

res_final <- res_final[order(res_final$padj), ]


write.csv(head(res_final, 50), "WGS_of_tumor_and_adjacent_breast_cancer_samples_Top50_Genes.csv")

write.csv(res_final, "WGS_of_tumor_and_adjacent_breast_cancer_samples.csv")
