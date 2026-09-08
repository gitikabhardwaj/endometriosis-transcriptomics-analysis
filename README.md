# Precision Medicine for Endometriosis

## Overview

This project analyzes publicly available RNA-seq data to investigate transcriptomic differences associated with endometriosis.

Using the NCBI GEO dataset GSE153740, I developed a reproducible workflow to compare mid-secretory eutopic endometrial samples from women with endometriosis and controls (n=4 per group).

The analysis includes sample-level quality assessment, differential gene expression using DESeq2, gene annotation, and ranked gene-set enrichment analysis using MSigDB Hallmark pathways.

Given the small cohort, results are interpreted as exploratory rather than as clinically validated biomarkers.

## Research Questions

- Do endometriosis and control samples show distinct global transcriptomic patterns?
- Which genes differ in expression between endometriosis and control samples?
- Which biological pathways show coordinated enrichment across the ranked transcriptome?
- How should gene- and pathway-level findings be interpreted given sample-level heterogeneity and limited cohort size?

## Analysis Workflow

1. Selected and evaluated human endometriosis RNA-seq dataset GSE153740
2. Reviewed sample metadata and experimental design
3. Assessed sequencing depth and filtered low-count genes
4. Applied variance-stabilizing transformation for sample-level PCA
5. Performed differential expression analysis with DESeq2
6. Annotated gene-level differential-expression results
7. Ranked genes using DESeq2 Wald statistics
8. Performed Hallmark gene-set enrichment analysis with fgsea and MSigDB
9. Generated PCA, volcano, and pathway-enrichment visualizations

## Tools

- R
- DESeq2
- fgsea
- msigdbr / MSigDB Hallmark gene sets
- Python
- pandas
- NumPy
- matplotlib / seaborn
- Jupyter Notebook
- Git / GitHub

## Results Summary

### Sample-level variation

Variance-stabilizing transformation followed by PCA showed substantial sample-level heterogeneity. PC1 and PC2 explained 68.6% and 19.0% of the total variance, respectively.

Endometriosis and control samples did not separate cleanly in PCA space, suggesting that disease status was not the dominant source of global transcriptomic variation in this small cohort.

![VST PCA](results/figures/vst_pca.png)

### Differential expression

Differential expression analysis was performed using DESeq2 on NCBI-generated raw gene counts.

Genes were retained if they had at least 10 counts in at least 4 samples. At an adjusted p-value threshold of 0.05:

- 12 genes were differentially expressed
- 7 showed higher expression in endometriosis
- 5 showed lower expression in endometriosis

Leading signals included S100P, ORM1, MT1X, MT1F, STC1, LRRN1, ENC1, CBLN1, and LYPD6.

![Volcano plot](results/figures/volcano_plot.png)

### Pathway enrichment

Ranked gene-set enrichment analysis was performed using the DESeq2 Wald statistic and MSigDB Hallmark gene sets.

Thirteen Hallmark pathways were significant at FDR < 0.05.

Positive enrichment toward endometriosis included TNF-alpha/NF-kB signaling, late estrogen response, hypoxia, and xenobiotic metabolism.

Negative enrichment included E2F targets, G2M checkpoint, mitotic spindle, MYC targets, DNA repair, Notch signaling, and other pathways.

![Hallmark enrichment](results/figures/hallmark_enrichment.png)

### Interpretation

The cohort does not show clear global transcriptomic separation by disease status, and the small sample size limits statistical power and generalizability.

However, the analysis identifies a focused set of gene-level differences and coordinated pathway-level signals that may warrant further investigation in larger independent cohorts.

These findings should be considered exploratory and should not be interpreted as clinically validated biomarkers.
## Repository Structure

```text
data/        Dataset documentation and metadata
notebooks/   Exploratory analyses and notebooks
src/         Reusable analysis scripts
results/     Generated figures and analysis outputs
