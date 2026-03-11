# Module 04: Functional Enrichment & Pathway Analysis
# Purpose: GO/KEGG/GSEA to identify biological mechanisms

library(clusterProfiler)
library(org.Hs.eg.db)

res <- readRDS("data/res_raw.rds")
clean_ids <- gsub("\\..*", "", rownames(res))

sig_mask <- which(res$padj < 0.05 & abs(res$log2FoldChange) > 1)
sig_ensembl <- clean_ids[sig_mask]


geneIDs <- bitr(sig_ensembl,
                fromType = "ENSEMBL",
                toType   = "ENTREZID",
                OrgDb    = org.Hs.eg.db)

# Run GO Enrichment with "Readable" Gene Symbols
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
library(ggplot2)
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


#GSEA 

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