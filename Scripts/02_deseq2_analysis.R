# Module 02: Differential Expression Analysis
# Purpose: Statistical modeling of Tumor vs Normal counts

library(DESeq2)

txi <- readRDS("data/txi_summarized.rds")

# 1. Experimental Design
coldata <- data.frame(
  row.names = colnames(txi$counts),
  condition = factor(c(rep("Tumor", 3), rep("Normal", 3)),
                     levels = c("Normal", "Tumor"))
)

# 2. DESeq2 Workflow
dds <- DESeqDataSetFromTximport(txi, colData = coldata, design = ~ condition)

# 3. Pre-filtering (Reducing noise/low-count genes)
# Rationale: Removes genes with insufficient data to provide statistical power
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# 4. Run Differential Expression
dds <- DESeq(dds)
res <- results(dds, alpha = 0.05) # Setting significance threshold at 5%

saveRDS(dds, "data/dds_analyzed.rds")
saveRDS(res, "data/res_raw.rds")
