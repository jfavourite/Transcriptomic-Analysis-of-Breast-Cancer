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
#### salmon quant \
 #### -i salmon_index \
 #### -l A \
 #### -1 sample_1.fastq \
 #### -2 sample_2.fastq \
 #### -p 3 \
#### -o sample_output

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
<img width="431" height="236" alt="plot pca of the 6 files" src="https://github.com/user-attachments/assets/34dc9bba-dc1b-4605-8ce3-dd76e558a8c9" />
### Differential Expression
3008 genes significant at FDR < 0.05
2640 genes with |log2FC| > 1
This indicates extensive transcriptomic remodeling in tumor tissue.
<img width="431" height="236" alt="Volcano plot WGS" src="https://github.com/user-attachments/assets/f67e2e64-646f-4bc4-8dfb-df5b8297bfaf" />
<img width="1030" height="563" alt="Heatmap" src="https://github.com/user-attachments/assets/30717f1b-f2d6-4c0d-9766-bdeb59470a72" />

### Biological themes:
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
<img width="961" height="525" alt="GO enrichment" src="https://github.com/user-attachments/assets/be1b87bf-d6c0-4c67-b90e-1fbfd2074c99" />

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
<img width="824" height="450" alt="KEGG ENRICHMENT PLOT" src="https://github.com/user-attachments/assets/b9d5ef3a-f660-4f5c-9985-b6954c30376c" />

## 8️⃣ Gene Set Enrichment Analysis (GSEA)

Ranked gene list analysis confirmed coordinated enrichment of immune-related and adhesion-related biological processes, reinforcing over-representation findings.
<img width="1373" height="750" alt="GSEA results" src="https://github.com/user-attachments/assets/bd4f8955-7971-4265-a079-b4e9369181d9" />

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
