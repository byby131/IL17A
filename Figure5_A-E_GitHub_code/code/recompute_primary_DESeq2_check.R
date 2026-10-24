#!/usr/bin/env Rscript
# Transparent independent reimplementation of the archived primary DESeq2 model.
# Never overwrites the frozen source data or the published statistics.
# The original run's executable source was not present in the archived delivery.

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
})

root <- normalizePath(".", mustWork = TRUE)
counts_path <- file.path(root, "data/shared/original_counts_input.csv")
metadata_path <- file.path(root, "data/shared/sample_metadata.csv")
frozen_path <- file.path(root, "data/A_volcano/all_tested_gene_statistics.csv")
stopifnot(all(file.exists(c(counts_path, metadata_path, frozen_path))))

counts <- fread(counts_path)
gene_ids <- as.character(counts[[1L]])
stopifnot(!anyDuplicated(gene_ids), !anyNA(gene_ids))
mat <- as.matrix(counts[, -1L, with = FALSE])
storage.mode(mat) <- "integer"
rownames(mat) <- gene_ids
metadata <- fread(metadata_path)
stopifnot(all(metadata$count_column %in% colnames(mat)), nrow(metadata) == 6L)
mat <- mat[, metadata$count_column, drop = FALSE]
metadata[, condition := factor(condition, levels = c("vehicle", "IL17A"))]
stopifnot(all(!is.na(metadata$condition)))

keep <- rowSums(mat >= 10L) >= 3L
stopifnot(sum(keep) == 17887L)
dds <- DESeqDataSetFromMatrix(
  countData = mat[keep, , drop = FALSE],
  colData = data.frame(condition = metadata$condition, row.names = metadata$count_column),
  design = ~ condition
)
dds <- DESeq(dds, test = "Wald", fitType = "parametric", sfType = "ratio", quiet = TRUE)
res <- results(dds, contrast = c("condition", "IL17A", "vehicle"),
               alpha = 0.05, independentFiltering = FALSE)

recomputed <- data.table(
  gene_id = rownames(res), baseMean = res$baseMean,
  log2FC = res$log2FoldChange, lfcSE = res$lfcSE,
  stat = res$stat, pvalue = res$pvalue, padj = res$padj
)
frozen <- fread(frozen_path)
setkey(recomputed, gene_id); setkey(frozen, gene_id)
stopifnot(identical(recomputed$gene_id, frozen$gene_id))
fields <- c("baseMean", "log2FC", "lfcSE", "stat", "pvalue", "padj")
comparison <- rbindlist(lapply(fields, function(field) {
  a <- recomputed[[field]]; b <- frozen[[field]]
  ok <- is.finite(a) & is.finite(b)
  data.table(field = field, finite_pairs = sum(ok),
             max_abs_diff = if (any(ok)) max(abs(a[ok] - b[ok])) else NA_real_,
             na_discordant = sum(is.na(a) != is.na(b)))
}))

out <- file.path(root, "recomputed_analysis")
dir.create(out, showWarnings = FALSE)
fwrite(recomputed, file.path(out, "DESeq2_independent_recompute.csv"))
fwrite(comparison, file.path(out, "comparison_to_frozen_statistics.csv"))
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))
print(comparison)
message("Independent recomputation saved under recomputed_analysis/; archived statistics remain untouched.")
