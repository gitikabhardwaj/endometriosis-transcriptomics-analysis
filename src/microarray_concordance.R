## microarray_concordance.R
##
## Replication of the GSE153740 discovery gene set in GSE6364
## (Affymetrix U133 Plus 2.0, mid-secretory eutopic endometrium, 9 vs 8).
##
## Result, reproducible from this script:
##   1 gene at FDR < 0.05 array-wide
##   4 of 11 testable candidates concordant in direction
##   effect sizes attenuated 5-10x relative to the discovery cohort
##
## Usage:
##   Rscript src/microarray_concordance.R
##
## Prerequisites:
##   - results/tables/deseq2_results.csv   (discovery DESeq2 output; GeneID is Entrez)
##   - internet access on first run (GEOquery downloads GSE6364 and GPL570)
##
## Notes on the series matrix:
##   GSE6364's `characteristics_ch1` encodes both cycle phase and disease status in a
##   single string, e.g. "Mid Secretory Phase Endometriosis". The expr matrix is on the
##   LINEAR scale (median ~194, max ~150000), so it must be log2(x + 1) transformed
##   before limma. Do NOT run normalizeBetweenArrays: the series is already RMA-processed.

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(limma)
  library(org.Hs.eg.db)
})

## ------------------------------------------------------------------ config ----
acc        <- "GSE6364"
phase_pat  <- "Mid Secretory"          # mid-secretory samples only
case_pat   <- "Endometriosis"
ctrl_pat   <- "Normal"
disc_file  <- "results/tables/deseq2_results.csv"
out_dir    <- "results/replication"
species    <- "Homo sapiens"
expect_n   <- 17                       # 9 endometriosis + 8 control
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

## ----------------------------------------------------- discovery gene set ----
disc_raw <- fread(disc_file)
required <- c("GeneID", "log2FoldChange", "padj")
if (!all(required %in% names(disc_raw)))
  stop("discovery table needs columns: ", paste(required, collapse = ", "))

disc_raw[, gene := mapIds(org.Hs.eg.db, as.character(GeneID), "SYMBOL", "ENTREZID")]
disc <- disc_raw[!is.na(padj) & padj < 0.05 & !is.na(gene),
                 .(gene, disc_lfc = log2FoldChange)]
disc <- disc[!duplicated(gene)]
message("discovery genes at FDR < 0.05 with a symbol: ", nrow(disc))
print(disc)
print(table(direction = ifelse(disc$disc_lfc > 0, "higher", "lower")))

## ------------------------------------------------------------ load GSE6364 ----
eset <- getGEO(acc, GSEMatrix = TRUE, getGPL = TRUE, destdir = "data/external")[[1]]
ch   <- pData(eset)[["characteristics_ch1"]]
if (is.null(ch)) stop("no characteristics_ch1 field in ", acc)

keep <- grepl(phase_pat, ch, ignore.case = TRUE)
grp  <- ifelse(grepl(case_pat, ch, ignore.case = TRUE), "endometriosis",
        ifelse(grepl(ctrl_pat, ch, ignore.case = TRUE), "control", NA))
message("\nphase x group table:")
print(table(keep = keep, grp = grp, useNA = "ifany"))

sel <- which(keep)
if (length(sel) != expect_n)
  warning("expected ", expect_n, " mid-secretory samples, found ", length(sel))

ex_raw <- exprs(eset)
stopifnot(ncol(ex_raw) == length(ch))

## ---- linear -> log2, then subset to the selected samples ----
message("\nraw intensity range: ", paste(round(range(ex_raw, na.rm = TRUE), 1),
                                        collapse = " - "))
ex  <- log2(ex_raw[, sel, drop = FALSE] + 1)
grp <- factor(grp[sel], levels = c("control", "endometriosis"))
stopifnot(ncol(ex) == length(grp), ncol(ex) == expect_n)

message("samples used: ", ncol(ex))
print(table(grp))
message("log2 intensity range: ", paste(round(range(ex, na.rm = TRUE), 2),
                                        collapse = " - "))

## ---------------------------------------------------- probes -> symbols ----
fd  <- fData(eset)
sym <- grep("gene.?symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
if (is.na(sym)) stop("no Gene Symbol column in fData; columns: ",
                     paste(colnames(fd), collapse = ", "))
message("\nsymbol column: ", sym)

map <- data.frame(probe  = rownames(fd),
                  symbol = fd[[sym]],
                  stringsAsFactors = FALSE)
map <- map[!is.na(map$symbol) & map$symbol != "" & map$probe %in% rownames(ex), ]
message("probes with a symbol: ", nrow(map))

ex <- ex[map$probe, , drop = FALSE]
rownames(ex) <- map$symbol
message("duplicate symbols collapsed: ", sum(duplicated(rownames(ex))))
ex <- ex[!duplicated(rownames(ex)), , drop = FALSE]
stopifnot(ncol(ex) == length(grp))
message("genes entering limma: ", nrow(ex))

## ------------------------------------------------------------------ limma ----
design <- model.matrix(~ grp)
fit    <- eBayes(lmFit(ex, design))
tt     <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
tt$gene <- rownames(tt)

message("\nlimma logFC distribution:")
print(summary(tt$logFC))
message("genes at FDR < 0.05 array-wide: ", sum(tt$adj.P.Val < 0.05))
fwrite(tt[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val")],
       file.path(out_dir, paste0(acc, "_midsecretory_limma.csv")))

## --------------------------------------------------------- concordance ----
conc <- merge(disc, tt[, c("gene", "logFC", "P.Value", "adj.P.Val")],
              by = "gene", all.x = TRUE)
conc$same_sign  <- sign(conc$disc_lfc) == sign(conc$logFC)
conc$attenuation <- abs(conc$disc_lfc) / abs(conc$logFC)
conc <- conc[order(-abs(conc$disc_lfc)), ]

message("\n--- candidate concordance ---")
print(conc)

testable <- sum(!is.na(conc$logFC))
n_conc   <- sum(conc$same_sign, na.rm = TRUE)
cat("\nsame direction:", n_conc, "of", testable, "\n")
cat("not on array:", sum(is.na(conc$logFC)), "\n")
cat("median effect attenuation (discovery / array):",
    round(median(conc$attenuation, na.rm = TRUE), 1), "x\n")

## binomial check against chance
if (testable > 0) {
  bt <- binom.test(n_conc, testable, 0.5)
  cat("binomial test vs chance:", round(bt$p.value, 3), "\n")
}

fwrite(conc, file.path(out_dir, paste0(acc, "_midsecretory_concordance.csv")))

## --------------------------------------------------------------- figure ----
png(file.path("results/figures", paste0(acc, "_concordance.png")),
    width = 1400, height = 900, res = 150)
op <- par(mar = c(4, 8, 3, 1))
o <- order(conc$disc_lfc)
bp <- barplot(conc$disc_lfc[o], horiz = TRUE, las = 1,
              col = ifelse(conc$disc_lfc[o] > 0, "#B2182B", "#2166AC"),
              xlab = "log2 fold change (endometriosis vs control)",
              main = paste0(acc, ": discovery effect sizes"))
barplot(conc$logFC[o], horiz = TRUE, las = 1, add = TRUE,
        col = ifelse(conc$logFC[o] > 0, "#EF8A62", "#67A9CF"),
        border = NA, density = 30)
axis(2, at = bp, labels = conc$gene[o], las = 1, tick = FALSE, line = -0.5)
abline(v = 0, lty = 2)
legend("bottomright", fill = c("#B2182B", "#2168AC"),
       legend = c("discovery log2FC", "GSE6364 log2FC"), bty = "n", cex = 0.8)
par(op); dev.off()

writeLines(capture.output(sessionInfo()),
           file.path(out_dir, "sessionInfo_GSE6364.txt"))
message("\nDone. Outputs in ", out_dir)
