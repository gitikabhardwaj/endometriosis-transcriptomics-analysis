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
