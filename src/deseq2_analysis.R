library(DESeq2)

counts_file <- "data/raw/GSE153740_raw_counts_GRCh38.p13_NCBI.tsv.gz"

counts_df <- read.delim(
  counts_file,
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

dim(counts_df)
head(counts_df)
sample_metadata <- data.frame(
  row.names = colnames(counts_df),
  sample = c("En1", "En2", "En3", "En6", "F10", "F14", "F15", "F9"),
  group = factor(
    c(
      "Endometriosis",
      "Endometriosis",
      "Endometriosis",
      "Endometriosis",
      "Control",
      "Control",
      "Control",
      "Control"
    ),
    levels = c("Control", "Endometriosis")
  )
)

sample_metadata

dds <- DESeqDataSetFromMatrix(
  countData = counts_df,
  colData = sample_metadata,
  design = ~ group
)

dds
keep <- rowSums(counts(dds) >= 10) >= 4

dds_filtered <- dds[keep, ]

cat("Genes before filtering:", nrow(dds), "\n")
cat("Genes after filtering:", nrow(dds_filtered), "\n")
cat("Genes removed:", nrow(dds) - nrow(dds_filtered), "\n")
dds_filtered <- DESeq(dds_filtered)

results_deseq <- results(
  dds_filtered,
  contrast = c("group", "Endometriosis", "Control")
)

summary(results_deseq)
results_df <- as.data.frame(results_deseq)

results_df$GeneID <- rownames(results_df)

results_df <- results_df[
  order(results_df$padj),
]

head(results_df, 20)
dir.create("results/tables", recursive = TRUE, showWarnings = FALSE)

write.csv(
  results_df,
  "results/tables/deseq2_results.csv",
  row.names = FALSE
)
annotation_file <- "data/Human.GRCh38.p13.annot.tsv.gz"

annotation <- read.delim(
  annotation_file,
  header = TRUE,
  check.names = FALSE
)

# Make sure GeneID has the same type in both tables
results_df$GeneID <- as.character(results_df$GeneID)
annotation$GeneID <- as.character(annotation$GeneID)

# Add selected gene annotation fields to DESeq2 results
annotated_results <- merge(
  results_df,
  annotation[, c(
    "GeneID",
    "Symbol",
    "Description",
    "GeneType",
    "EnsemblGeneID"
  )],
  by = "GeneID",
  all.x = TRUE,
  sort = FALSE
)

# Re-sort by adjusted p-value
annotated_results <- annotated_results[
  order(annotated_results$padj),
]

# Show the top 20 results
print(
  head(
    annotated_results[, c(
      "GeneID",
      "Symbol",
      "Description",
      "GeneType",
      "baseMean",
      "log2FoldChange",
      "pvalue",
      "padj"
    )],
    20
  )
)

# Save annotated results
write.csv(
  annotated_results,
  "results/tables/deseq2_results_annotated.csv",
  row.names = FALSE
)
sig_005 <- subset(
  annotated_results,
  !is.na(padj) & padj < 0.05
)

sig_010 <- subset(
  annotated_results,
  !is.na(padj) & padj < 0.10
)

cat("Genes with padj < 0.05:", nrow(sig_005), "\n")
cat("Genes with padj < 0.10:", nrow(sig_010), "\n")

cat(
  "Upregulated at padj < 0.05:",
  sum(sig_005$log2FoldChange > 0),
  "\n"
)

cat(
  "Downregulated at padj < 0.05:",
  sum(sig_005$log2FoldChange < 0),
  "\n"
)
# Volcano plot

volcano_df <- annotated_results

volcano_df$significance <- "Not significant"
volcano_df$significance[
  !is.na(volcano_df$padj) &
  volcano_df$padj < 0.05 &
  volcano_df$log2FoldChange > 0
] <- "Upregulated"

volcano_df$significance[
  !is.na(volcano_df$padj) &
  volcano_df$padj < 0.05 &
  volcano_df$log2FoldChange < 0
] <- "Downregulated"

volcano_df$minus_log10_padj <- -log10(volcano_df$padj)

dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

png(
  "results/figures/volcano_plot.png",
  width = 1600,
  height = 1200,
  res = 180
)

plot(
  volcano_df$log2FoldChange,
  volcano_df$minus_log10_padj,
  pch = 16,
  col = ifelse(
    volcano_df$significance == "Upregulated",
    "firebrick",
    ifelse(
      volcano_df$significance == "Downregulated",
      "steelblue",
      "grey75"
    )
  ),
  xlab = "Log2 Fold Change (Endometriosis vs Control)",
  ylab = "-Log10 Adjusted P-value",
  main = "Differential Expression in Endometriosis"
)

abline(
  h = -log10(0.05),
  lty = 2,
  col = "darkgrey"
)

abline(
  v = 0,
  lty = 2,
  col = "darkgrey"
)
# Label genes significant at FDR < 0.05
label_df <- volcano_df[
  !is.na(volcano_df$padj) &
  volcano_df$padj < 0.05 &
  !is.na(volcano_df$Symbol),
]

text(
  x = label_df$log2FoldChange,
  y = label_df$minus_log10_padj,
  labels = label_df$Symbol,
  pos = 3,
  cex = 0.75,
  offset = 0.5
)
dev.off()
# Create ranked gene list for downstream pathway analysis

ranked_results <- annotated_results[
  !is.na(annotated_results$stat) &
  !is.na(annotated_results$Symbol) &
  annotated_results$Symbol != "",
]

ranked_results <- ranked_results[
  order(ranked_results$stat, decreasing = TRUE),
]

ranked_gene_list <- ranked_results[, c(
  "GeneID",
  "Symbol",
  "stat",
  "log2FoldChange",
  "pvalue",
  "padj"
)]

write.csv(
  ranked_gene_list,
  "results/tables/ranked_gene_list.csv",
  row.names = FALSE
)

cat("Genes available for pathway ranking:", nrow(ranked_gene_list), "\n")

print(head(ranked_gene_list, 10))
print(tail(ranked_gene_list, 10))
