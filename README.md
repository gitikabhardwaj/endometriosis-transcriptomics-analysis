# Precision Medicine for Endometriosis

## Overview

This project explores publicly available transcriptomic datasets to investigate molecular signatures associated with endometriosis.

The goal is to build a reproducible bioinformatics workflow for transcriptomic data processing, exploratory analysis, differential gene expression, pathway analysis, and eventual evaluation of machine learning approaches for biomarker discovery.

This project is currently in progress.

## Research Questions

This project aims to explore questions such as:

- Which genes show differential expression between endometriosis and control samples?
- Which biological pathways are enriched among differentially expressed genes?
- Are there transcriptomic patterns that distinguish disease and control samples?
- Can machine learning approaches help identify candidate molecular signatures associated with endometriosis?

## Planned Workflow

1. Identify and acquire relevant public datasets from NCBI GEO
2. Inspect sample metadata and experimental design
3. Perform data preprocessing and quality control
4. Conduct exploratory data analysis
5. Perform differential gene expression analysis
6. Conduct pathway and gene-set enrichment analysis
7. Explore machine learning approaches for transcriptomic classification
8. Generate reproducible visualizations and documentation

## Tools

Planned tools include:

- Python
- R
- pandas
- NumPy
- matplotlib / seaborn
- scikit-learn
- Bioconductor packages
- Jupyter Notebook

Additional tools will be added as the analysis develops.
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
