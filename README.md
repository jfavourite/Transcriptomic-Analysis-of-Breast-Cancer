# Transcriptomic-Analysis-of-Breast-Cancer
## Tumor vs Adjacent Normal Tissue (RNA-seq, Homo sapiens)
This project presents a full end-to-end RNA-seq analysis workflow comparing:
3 Breast Tumor samples
3 Adjacent Normal Breast Tissue samples
Organism: Homo sapiens

The analysis covers:

Raw data acquisition from SRA
FASTQ generation
Transcript-level quantification
Gene-level summarization
Differential expression analysis
Functional enrichment (GO, KEGG)
Gene Set Enrichment Analysis (GSEA)

This repository demonstrates computational reproducibility, statistical rigor, and biological interpretation suitable for high-level research.

## 1️⃣ Data Acquisition (Command Line Workflow)

All command-line processing was performed in WSL (Ubuntu).

Step 1: Download raw sequencing data
Using NCBI SRA Toolkit:
prefetch SRRXXXXXXX
Then FASTQ generation:
#### fasterq-dump SRRXXXXXXX --split-files -O raw -e 3
This generated paired-end FASTQ files:
#### sample_1.fastq
#### sample_2.fastq

## 2️⃣ Transcriptome Reference Preparation
Reference annotation:
GENCODE v46
#### gencode.v46.transcripts.fa
#### gencode.v46.annotation.gtf
Salmon index built using:
#### salmon index \
 #### -t gencode.v46.transcripts.fa \
 #### -i salmon_index \
 #### -k 31
 
## 3️⃣ Transcript Quantification
Transcript-level quantification was performed using: Salmon v1.10.3
For each sample:
salmon quant \
  -i salmon_index \
  -l A \
  -1 sample_1.fastq \
  -2 sample_2.fastq \
  -p 3 \
  -o sample_output

Output: quant.sf
This was done for all 6 samples.

## 4️⃣ Gene-Level Summarization

Transcript-to-gene mapping was created from:
gencode.v46.annotation.gtf
Using R:
tximport
GenomicFeatures
Transcripts were aggregated to gene-level counts.

## 5️⃣ Differential Expression Analysis

Differential expression performed using:
DESeq2
Design formula: design = ~ condition
Filtering: Genes with total counts ≤ 10 removed
Statistical thresholds:
FDR < 0.05
|log2FoldChange| > 1

## 6️⃣ Results
### Global Structure (PCA)
PC1 explains 54% of total variance
PC2 explains 22% of total variance
Tumor and normal samples cluster distinctly along PC1, indicating strong tumor-driven transcriptional reprogramming.
!|
Differential Expression
3008 genes significant at FDR < 0.05
2640 genes with |log2FC| > 1
This indicates extensive transcriptomic remodeling in tumor tissue.

Biological themes:
Extracellular matrix remodeling
Immune activation
Chemokine signaling

## 7️⃣ Functional Enrichment

Functional enrichment performed using:
clusterProfiler
#### GO Biological Processes

Top enriched:
Leukocyte cell-cell adhesion
Regulation of T cell activation
Leukocyte proliferation
Lymphocyte differentiation
Indicates strong immune microenvironment involvement.
!(

#### KEGG Pathways
Top enriched:
Systemic lupus erythematosus
Cell adhesion molecule interaction
Neutrophil extracellular trap formation
ECM-receptor interaction
Integrin signaling
These support:
Tumor–stroma interaction
Immune activation
Adhesion dysregulation

## 8️⃣ Gene Set Enrichment Analysis (GSEA)

Ranked gene list analysis confirmed coordinated enrichment of immune-related and adhesion-related biological processes, reinforcing over-representation findings.
!

## 9️⃣ Biological Interpretation

Tumor samples exhibit:
Strong extracellular matrix remodeling 
Chemokine-driven immune signaling 
Cell adhesion pathway dysregulation
Immune infiltration signatures
The high variance explained by PC1 (54%) suggests tumor status is the dominant driver of gene expression variation.
Overall, results are consistent with:
Tumor microenvironment remodeling
Immune engagement
Stromal activation in breast cancer

## 🔬 Technical Stack
Salmon
tximport
DESeq2
clusterProfiler
org.Hs.eg.db
pheatmap
ggplot2

All analyses were performed in R.
