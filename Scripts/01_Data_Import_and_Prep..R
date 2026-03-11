# Module 01: Data Retrieval and Transcript-to-Gene Mapping
# Purpose: Map GENCODE transcripts to Gene IDs and import Salmon quant.sf files

library(tximport)
library(GenomicFeatures)
library(txdbmaker)
library(AnnotationDbi)

# 1. Create TxDb from GENCODE GTF
# Technical Note: Using GENCODE v46 (Latest) ensures high annotation accuracy
txdb <- makeTxDbFromGFF("reference/gencode.v46.annotation.gtf")

k <- keys(txdb, keytype = "TXNAME")

tx2gene <- AnnotationDbi::select(txdb,
                                 keys = k,
                                 columns = "GENEID",
                                 keytype = "TXNAME")

tx2gene <- tx2gene[, c("TXNAME", "GENEID")]

# 2. Locate Quant Files
# Note: Ensure paths are relative for portability across different HPC environments
sample_ids <- c("T1","T2","T3","N1","N2","N3")

files <- file.path("quant", c("tumor1", "tumor2", "tumor3", "normal1", "normal2", "normal3"), "quant.sf")

names(files) <- sample_ids

# 3. Import with tximport (Gene-level summarization)
txi <- tximport(files, type = "salmon", 
                tx2gene = tx2gene,
                ignoreAfterBar = TRUE)

# Save intermediate object for next script
saveRDS(txi, "data/txi_summarized.rds")
