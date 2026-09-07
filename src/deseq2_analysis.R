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
