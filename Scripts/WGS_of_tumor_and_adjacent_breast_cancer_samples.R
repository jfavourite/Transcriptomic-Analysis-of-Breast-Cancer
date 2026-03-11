# Installing and loading the necessary packages
BiocManager::install(c("tximport", "DESeq2", "GenomicFeatures", "txdbmaker","AnnotationDbi"))

library(tximport)
library(DESeq2)
library(GenomicFeatures)
library(txdbmaker)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)

#Creating transcript for gene mapping

txdb <- makeTxDbFromGFF("reference/gencode.v46.annotation.gtf")

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(txdb,
                                 keys = k,
                                 columns = "GENEID",
                                 keytype = "TXNAME")

tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

#Arranging quant files

files <- c(
  "quant/tumor1/quant.sf",
  "quant/tumor2/quant.sf",
  "quant/tumor3/quant.sf",
  "quant/normal1/quant.sf",
  "quant/normal2/quant.sf",
  "quant/normal3/quant.sf"
)


names(files) <- c("T1","T2","T3","N1","N2","N3")

#Importing with tximport
txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                ignoreAfterBar = TRUE)

#Creating a metadata
condition <- factor(c("Tumor","Tumor","Tumor",
                      "Normal","Normal","Normal"))

coldata <- data.frame(row.names = names(files),
                      condition)

#Creating a DESeq object

dds <- DESeqDataSetFromTximport(txi,
                                colData = coldata,
                                design = ~ condition)

#Filtering low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

#Running the differential analysis

dds <- DESeq(dds)

res <- results(dds) 

#PCA plot to view normal and tumor

vsd <- vst(dds, blind = FALSE)

plotPCA(vsd, intgroup = "condition")
#Broader and more explanatory PCA plot
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

#Volcano plot

res_df <- as.data.frame(res)
res_df$diffexpressed <- "Not Significant"


res_df$diffexpressed[res_df$log2FoldChange > 1 & res_df$padj < 0.05] <- "Up-regulated (Tumor)"
res_df$diffexpressed[res_df$log2FoldChange < -1 & res_df$padj < 0.05] <- "Down-regulated (Normal)"

ggplot(res_df, aes(x=log2FoldChange, y=-log10(padj), color=diffexpressed)) +
  geom_point(alpha=0.5, size=1.5) +
  theme_minimal() +
  scale_color_manual(values=c("blue", "grey", "red")) + # Custom colors
  geom_vline(xintercept=c(-1, 1), col="black", linetype="dashed") + # Fold Change lines
  geom_hline(yintercept=-log10(0.05), col="black", linetype="dashed") + # Significance line
  labs(title="Tumor vs Normal WGS of breast cancer samples",
       x="Log2 Fold Change", 
       y="-Log10 Adjusted P-value",
       color="Gene Status")
# Creating a heatmap for the top 50 genes

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

#GO ENRICHMENT ANALYSIS(FUCTIONAL ANALYSIS)

clean_ids <- gsub("\\..*", "", rownames(res))
sig_mask <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
sig_ensembl <- clean_ids[sig_mask]


geneIDs <- bitr(sig_ensembl,
                fromType = "ENSEMBL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# Run Enrichment with "Readable" Gene Symbols
ego <- enrichGO(gene          = geneIDs$ENTREZID,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",        
                pAdjustMethod = "BH",       
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05,       
                readable      = TRUE)        

# Save outputs
write.csv(as.data.frame(ego), "GO_enrichment_results.csv")

# Dotplot
# showCategory=20 ensures the top 20 processes are seen

dotplot(ego, showCategory=20) + 
  scale_color_gradient(low="red", high="blue") + 
  ggtitle("GO Biological Processes: Tumor vs adjacent Normal og WGS breast cancer samples") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"))

# KEGG enrichment
ekegg <- enrichKEGG(gene         = geneIDs$ENTREZID,
                    organism     = 'hsa',         
                    pvalueCutoff = 0.05,
                    pAdjustMethod = "BH")

# Visualize the results
dotplot(ekegg, showCategory = 15) + 
  ggtitle("KEGG Pathway Enrichment: Tumor vs Adjacent Normal; WGS breast cancer samples")


#GSEA analysis

# Preparimg the Gene List

geneList <- res$log2FoldChange

#  Assigning the CLEANED IDs to the list
names(geneList) <- clean_ids

# Sorting in decreasing order
geneList <- sort(geneList, decreasing = TRUE)

# Run GSEA
gsea_res <- gseGO(geneList      = geneList,
                  OrgDb         = org.Hs.eg.db,
                  ont           = "BP",        
                  keyType       = "ENSEMBL",  
                  minGSSize     = 10,
                  maxGSSize     = 500,
                  pvalueCutoff  = 0.05,
                  verbose       = FALSE)

dotplot(gsea_res, showCategory=10, split=".sign") + 
  facet_grid(.~.sign) +
  ggtitle("GSEA Results: WGSof breast cancer samples(Tumor vs adjacent normal) Pathway Dynamics")

# Save the full GSEA table
write.csv(as.data.frame(gsea_res), "GSEA_Full_Results.csv", row.names = FALSE)

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

