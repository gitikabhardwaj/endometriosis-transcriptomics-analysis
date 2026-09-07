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
