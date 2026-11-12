## run_replication.R
##
## Replication of the GSE153740 discovery analysis in GSE234354 (Mortlock et al. cohort).
##
## GSE234354 is the bulk RNA-seq from 206 women (143 endometriosis cases, 63 controls).
## The GEO release carries cycle stage, age, pod and batch, but NOT endometriosis
## case/control status. This script is written to run the moment those labels arrive.
##
## Two label sources are supported:
##   (a) a CSV from the dataset authors: columns `gsm` (or `key`) and `group`
##   (b) any file with per-sample disease status keyed on GSM ID or biopsy name
##
## Usage:
##   Rscript src/run_replication.R
##
## Prerequisites:
##   - data/external/GSE234354/counts.tsv            (Entrez GeneID x GSM columns)
##   - data/external/GSE234354/metadata.csv          (written by this script from GEO)
##   - data/external/GSE234354/disease_labels.csv    (from the authors; optional)
##   - results/tables/deseq2_results.csv             (discovery DESeq2 output)

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(org.Hs.eg.db)
})

## ------------------------------------------------------------------ config ----
counts_file  <- "data/external/GSE234354/counts.tsv"
matrix_files <- c("data/external/GSE234354/GSE234354-GPL20301_series_matrix.txt",
                  "data/external/GSE234354/GSE234354-GPL11154_series_matrix.txt")
label_file   <- "data/external/GSE234354/disease_labels.csv"   # from the authors
disc_file    <- "results/tables/deseq2_results.csv"
phase_keep   <- "6"          # 6 = mid secretory (days 20-23) per the series description
out_dir      <- "results/replication"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ------------------------------------------------- metadata from GEO matrices ----
read_meta <- function(f) {
  if (!file.exists(f)) return(NULL)
  lines <- readLines(f)
  get <- function(pat) lapply(grep(pat, lines, value = TRUE),
                              function(l) gsub('"', "", strsplit(l, "\t")[[1]][-1]))
  m <- data.frame(sample = get("^!Sample_geo_accession")[[1]],
                  title  = get("^!Sample_title")[[1]],
                  stringsAsFactors = FALSE)
  stopifnot(nrow(m) > 0, ncol(m) > 0)
  for (v in get("^!Sample_characteristics_ch1")) {
    key <- sub(":.*", "", v[1])
    if (length(v) == nrow(m)) m[[key]] <- sub("^[^:]*: ?", "", v)
  }
  desc <- get("^!Sample_description")
  if (length(desc) == nrow(m))
    m$key <- sub(".*Sample ([^ ]+) from.*", "\\1", desc)
  else
    m$key <- sub("^endometrial biopsy, ", "", m$title)
  m
}

metas <- lapply(matrix_files, read_meta)
metas <- Filter(Negate(is.null), metas)
if (!length(metas)) stop("no series matrix files found in data/external/GSE234354/")

meta <- Reduce(function(a, b) {
  cols <- union(names(a), names(b))
  for (cl in setdiff(cols, names(a))) a[[cl]] <- NA
  for (cl in setdiff(cols, names(b))) b[[cl]] <- NA
  rbind(a[cols], b[cols])
}, metas)

meta$stage  <- suppressWarnings(as.integer(meta[["cycle stage"]]))
meta$woman  <- sub("_2$", "", meta$key)
meta$batch  <- meta[["batch"]]
message("GEO metadata rows: ", nrow(meta))

## ------------------------------------------------ attach disease labels if present ----
read_labels <- function(f) {
  if (!file.exists(f)) return(NULL)
  lab <- fread(f)
  cand <- intersect(c("gsm", "sample", "key", "title", "biopsy"), names(lab))
  if (!length(cand)) stop("label file needs one of: gsm, sample, key, title, biopsy")
  id <- cand[1]
  gcol <- intersect(c("group", "disease", "endometriosis", "status", "condition"),
                    names(lab))
  if (!length(gcol)) stop("label file needs one of: group, disease, status, condition")
  out <- data.frame(id = as.character(lab[[id]]),
                    group_raw = as.character(lab[[gcol[1]]]),
                    stringsAsFactors = FALSE)
  out
}

lab <- read_labels(label_file)
if (is.null(lab)) {
  message("\nNo disease label file at ", label_file,
          "\n -> writing metadata only; the disease test cannot run yet.")
  fwrite(meta, file.path(out_dir, "GSE234354_metadata.csv"))
} else {
  key2gsm <- unique(meta[, c("key", "sample")])
  lab$sample <- lab$id
  hit <- lab$id %in% meta$sample
  if (!all(hit)) {
    m <- merge(lab[!hit, "id", drop = FALSE], key2gsm, by.x = "id", by.y = "key", all.x = TRUE)
    lab$sample[!hit] <- m$sample[match(lab$id[!hit], m$id)]
  }
  meta <- merge(meta, lab[, c("sample", "group_raw")], by = "sample", all.x = TRUE)
  meta$group <- ifelse(grepl("endometrio|case|yes|1", meta$group_raw, ignore.case = TRUE),
                       "endometriosis",
                ifelse(grepl("control|no|0", meta$group_raw, ignore.case = TRUE),
                       "control", NA))
  message("group counts over all labelled samples:")
  print(table(meta$group, useNA = "ifany"))
  fwrite(meta, file.path(out_dir, "GSE234354_metadata.csv"))
}

## ------------------------------------------------------------- counts matrix ----
cts <- fread(counts_file)
gene_ids <- as.character(cts[[1]])
cts <- as.matrix(cts[, -1, with = FALSE])
rownames(cts) <- gene_ids
stopifnot(all(abs(cts - round(cts)) < 1e-6))    # must be integer counts
message("counts: ", nrow(cts), " genes x ", ncol(cts), " samples")

## ------------------------------------------------------------- disease test ----
if (!is.null(lab) && sum(!is.na(meta$group)) >= 6) {
  sel <- meta[meta$sample %in% colnames(cts) & meta$stage == as.integer(phase_keep) &
              !is.na(meta$group), ]
  sel <- sel[!duplicated(sel$woman), ]          # one biopsy per woman
  sel <- sel[order(sel$sample), ]
  message("\nmid-secretory samples used:")
  print(table(sel$group))

  if (nlevels(factor(sel$group)) < 2 || min(table(sel$group)) < 3) {
    message("-> too few samples per group for a disease test; stopping here.")
  } else {
    m <- round(cts[, sel$sample, drop = FALSE])
    m <- m[rowSums(m >= 10) >= 4, ]             # same filter as discovery

    coldata <- data.frame(row.names = sel$sample, group = factor(sel$group,
                             levels = c("control", "endometriosis")))
    if (length(unique(na.omit(sel$batch))) > 1) {
      coldata$batch <- factor(sel$batch)
      chk <- table(coldata$batch, coldata$group)
      print(chk)
      # only adjust for batch if it is not confounded with group
      design <- if (any(chk == 0)) ~ group else ~ batch + group
    } else {
      design <- ~ group
    }
    message("design: ", paste(deparse(design), collapse = ""))

    dds <- DESeq(DESeqDataSetFromMatrix(m, coldata, design))
    res <- as.data.frame(results(dds, contrast = c("group", "endometriosis", "control")))
    res$GeneID <- rownames(res)
    res$gene <- mapIds(org.Hs.eg.db, res$GeneID, "SYMBOL", "ENTREZID")
    res <- res[order(res$padj), ]
    fwrite(res, file.path(out_dir, "GSE234354_midsecretory_deseq2.csv"))
    message("significant at FDR < 0.05: ", sum(res$padj < 0.05, na.rm = TRUE))

    ## ---- concordance against the discovery gene set ----
    d <- fread(disc_file)
    d[, gene := mapIds(org.Hs.eg.db, as.character(GeneID), "SYMBOL", "ENTREZID")]
    disc <- d[!is.na(padj) & padj < 0.05 & !is.na(gene),
              .(gene, disc_lfc = log2FoldChange)][!duplicated(gene)]
    conc <- merge(disc, res[, c("gene", "log2FoldChange", "stat", "pvalue", "padj")],
                  by = "gene", all.x = TRUE)
    conc$same_sign <- sign(conc$disc_lfc) == sign(conc$log2FoldChange)
    conc <- conc[order(-abs(conc$disc_lfc)), ]
    print(conc)
    cat("same direction:", sum(conc$same_sign, na.rm = TRUE), "of",
        sum(!is.na(conc$log2FoldChange)), "| not testable:",
        sum(is.na(conc$log2FoldChange)), "\n")
    fwrite(conc, file.path(out_dir, "GSE234354_gene_concordance.csv"))

    ## ---- Hallmark GSEA on the Wald statistic ----
    ranks <- res$stat
    names(ranks) <- res$gene
    ranks <- ranks[!is.na(ranks)]
    ranks <- sort(ranks, decreasing = TRUE)
    hm <- msigdbr(species = "Homo sapiens", category = "H")
    sets <- split(hm$gene_symbol, hm$gs_name)
    set.seed(1)
    fg <- fgsea(sets, ranks, minSize = 15, maxSize = 500)
    fg <- fg[order(fg$padj), ]
    print(fg[, .(pathway, NES, padj)])
    fwrite(fg, file.path(out_dir, "GSE234354_hallmark_fgsea.csv"))

    disc_fg_path <- "results/tables/hallmark_fgsea_results.csv"
    if (file.exists(disc_fg_path)) {
      df <- fread(disc_fg_path)
      if (all(c("pathway", "NES", "padj") %in% names(df))) {
        pm <- merge(df[padj < 0.05, .(pathway, disc_NES = NES)],
                    fg[, .(pathway, NES, padj)], by = "pathway")
        pm$same_sign <- sign(pm$disc_NES) == sign(pm$NES)
        print(pm)
        cat("pathways same direction:", sum(pm$same_sign), "of", nrow(pm), "\n")
        fwrite(pm, file.path(out_dir, "GSE234354_pathway_concordance.csv"))
      }
    }
    writeLines(capture.output(sessionInfo()),
               file.path(out_dir, "sessionInfo_GSE234354.txt"))
  }
}

## ------------------------------------------------- secretory-stage drift check ----
## Runs whether or not disease labels exist. Uses all stage 5-7 biopsies.
sel2 <- meta[meta$sample %in% colnames(cts) & meta$stage %in% 5:7, ]
sel2 <- sel2[!duplicated(sel2$woman), ]
if (nrow(sel2) >= 12 && length(unique(sel2$stage)) == 3) {
  m2 <- round(cts[, sel2$sample, drop = FALSE])
  m2 <- m2[rowSums(m2 >= 10) >= 4, ]
  cd2 <- data.frame(row.names = sel2$sample, stage = factor(sel2$stage))
  v <- assay(varianceStabilizingTransformation(
        DESeqDataSetFromMatrix(m2, cd2, ~ stage), blind = TRUE))
  rownames(v) <- mapIds(org.Hs.eg.db, rownames(v), "SYMBOL", "ENTREZID")
  v <- v[!is.na(rownames(v)) & !duplicated(rownames(v)), ]

  genes <- intersect(c("S100P","ORM1","MT1X","MT1F","STC1","LRRN1","ENC1",
                       "CBLN1","LYPD6","PTK7"), rownames(v))
  tab <- sapply(genes, function(g) tapply(v[g, ], cd2$stage, median))
  out <- data.frame(gene = genes, t(tab), delta_7_minus_5 = tab["7", ] - tab["5", ])
  rownames(out) <- NULL
  print(out)
  fwrite(out, file.path(out_dir, "stage_check_stage5-7.csv"))

  stage_tab <- sapply(genes, function(g) {
    k <- kruskal.test(v[g, ] ~ cd2$stage)$p.value
    k
  })
  print(round(sort(stage_tab), 4))
  fwrite(data.frame(gene = names(stage_tab), kruskal_p = as.numeric(stage_tab)),
         file.path(out_dir, "stage_check_kruskal.csv"))
} else {
  message("stage-drift check skipped: ", nrow(sel2), " usable biopsies")
}

message("\nDone. Outputs in ", out_dir)
