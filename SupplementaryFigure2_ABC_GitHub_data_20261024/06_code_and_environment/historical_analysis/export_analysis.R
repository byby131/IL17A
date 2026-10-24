if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
invisible(Sys.setlocale("LC_COLLATE", "C"))
.libPaths(c("E:/codex-转录组/adsc_il17_rnaseq/.Rlibrary", .libPaths()))
suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(fgsea)
})

root <- "E:/codex-转录组"
source_csv <- file.path(root, "ADSC_IL17A_timecourse_RNAseq_EXPRESSION_MATRIX_COUNTS.CSV")
prior <- file.path(root, "HFF-1/round3_20260914/cache")
out <- file.path(root, "analysis/adsc_timecourse_counts")
dir.create(out, recursive=TRUE, showWarnings=FALSE)

fit <- readRDS(file.path(prior, "sender_fit.rds"))
dds <- readRDS(file.path(prior, "sender_dds.rds"))
ann <- read.csv(file.path(prior, "sender_gene_annotation.csv"), check.names=FALSE)
input_audit <- jsonlite::fromJSON(file.path(prior, "sender_input_audit.json"))
stopifnot(identical(deparse(design(dds)), "~replicate + time"))
stopifnot(identical(levels(colData(dds)$time), c("0","3","12","24")))
stopifnot(all(table(colData(dds)$replicate, colData(dds)$time) == 1))
stopifnot(nrow(dds) == 14425L, input_audit$exact_versioned_matches == 232117L)

write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE, na="")
meta <- data.frame(
  sample_id=colnames(dds),
  time_hours=as.integer(as.character(colData(dds)$time)),
  replicate=as.character(colData(dds)$replicate),
  condition="IL-17A 50 ng/mL",
  baseline=as.character(colData(dds)$time)=="0",
  stringsAsFactors=FALSE
)
write_csv(meta, "sample_metadata.csv")

raw <- counts(dds, normalized=FALSE)
norm <- counts(dds, normalized=TRUE)
qc <- data.frame(
  sample_id=colnames(dds),
  time_hours=meta$time_hours,
  replicate=meta$replicate,
  filtered_rounded_gene_count_sum=colSums(raw),
  size_factor=sizeFactors(dds),
  genes_count_gt_0=colSums(raw > 0),
  genes_count_ge_10=colSums(raw >= 10),
  normalized_count_sum=colSums(norm),
  stringsAsFactors=FALSE
)
write_csv(qc, "sample_QC.csv")

vsd <- vst(dds, blind=FALSE)
v <- assay(vsd)
pc <- prcomp(t(v), center=TRUE, scale.=FALSE)
pca <- cbind(meta, PC1=pc$x[,1], PC2=pc$x[,2], PC3=pc$x[,3])
attr(pca, "variance") <- NULL
write_csv(pca, "PCA_VST_coordinates.csv")
write_csv(data.frame(PC=paste0("PC", seq_along(pc$sdev)), variance_explained=pc$sdev^2/sum(pc$sdev^2)), "PCA_variance.csv")
cormat <- cor(v, method="pearson")
write.csv(cormat, file.path(out, "sample_correlation_VST.csv"), row.names=TRUE)

mean_by_time <- sapply(levels(colData(dds)$time), function(tt) rowMeans(norm[,colData(dds)$time==tt,drop=FALSE]))
colnames(mean_by_time) <- paste0("mean_normalized_count_", colnames(mean_by_time), "h")
ann_idx <- match(rownames(dds), ann$gene_id)
annotation <- data.frame(gene_id=rownames(dds), gene_name=ann$gene_name[ann_idx], gene_type=ann$gene_type[ann_idx], stringsAsFactors=FALSE)

contrast_names <- c("3_vs_0","12_vs_0","24_vs_0")
de_tables <- list()
for (nm in contrast_names) {
  rr <- fit$results[[nm]]
  stopifnot(identical(rr$gene_id, rownames(dds)))
  rr <- cbind(rr[,c("gene_id","gene_name","gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj","converged")], mean_by_time)
  rr$direction_FDR_0.05 <- ifelse(!is.na(rr$padj) & rr$padj < 0.05 & rr$log2FoldChange > 0, "up",
                                  ifelse(!is.na(rr$padj) & rr$padj < 0.05 & rr$log2FoldChange < 0, "down", "not_significant"))
  rr <- rr[order(rr$padj, -abs(rr$stat), na.last=TRUE),]
  de_tables[[nm]] <- rr
  write_csv(rr, paste0("DESeq2_", nm, ".csv"))
}

deg_summary <- do.call(rbind, lapply(contrast_names, function(nm) {
  x <- de_tables[[nm]]
  data.frame(contrast=nm, tested=nrow(x), up_FDR_lt_0.05=sum(x$direction_FDR_0.05=="up"), down_FDR_lt_0.05=sum(x$direction_FDR_0.05=="down"), not_significant=sum(x$direction_FDR_0.05=="not_significant"), median_abs_log2FC=median(abs(x$log2FoldChange),na.rm=TRUE))
}))
write_csv(deg_summary, "DEG_summary.csv")

top_genes <- do.call(rbind, lapply(contrast_names, function(nm) {
  x <- de_tables[[nm]]
  up <- head(x[x$direction_FDR_0.05=="up",], 15)
  down <- head(x[x$direction_FDR_0.05=="down",], 15)
  z <- rbind(up,down)
  z$contrast <- nm
  z$rank_within_direction <- ave(seq_len(nrow(z)), z$direction_FDR_0.05, FUN=seq_along)
  z[,c("contrast","direction_FDR_0.05","rank_within_direction","gene_id","gene_name","gene_type","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj",colnames(mean_by_time))]
}))
write_csv(top_genes, "top15_up_down_each_contrast.csv")

# Omnibus likelihood-ratio test for any time effect, keeping replicate blocking.
dds_lrt <- DESeq(dds, test="LRT", reduced=~replicate, sfType="ratio", minReplicatesForReplace=Inf, parallel=FALSE)
lrt <- results(dds_lrt, independentFiltering=FALSE)
lrt <- data.frame(annotation, baseMean=lrt$baseMean, deviance=lrt$stat, pvalue=lrt$pvalue, padj=lrt$padj, mean_by_time, check.names=FALSE)
lrt <- lrt[order(lrt$padj, -lrt$deviance, na.last=TRUE),]
write_csv(lrt, "DESeq2_LRT_any_time_effect.csv")

# Unique Ensembl gene -> Entrez mapping for rank-based Hallmark GSEA.
ens <- sub("[.][0-9]+$", "", rownames(dds))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), keytype="ENSEMBL", columns=c("ENTREZID","SYMBOL")))
map <- map[!is.na(map$ENTREZID),]
ambiguous <- names(which(table(map$ENSEMBL) > 1))
map_unique <- map[!map$ENSEMBL %in% ambiguous & !duplicated(map$ENSEMBL),]
gmt <- fgsea::gmtPathways(file.path(root, "adsc_il17_rnaseq/references/h.all.v2025.1.Hs.entrez.gmt"))
set.seed(20260920)
gsea_all <- list()
rank_audit <- list()
for (nm in contrast_names) {
  x <- fit$results[[nm]]
  r <- data.frame(ENSEMBL=sub("[.][0-9]+$", "", x$gene_id), stat=x$stat, stringsAsFactors=FALSE)
  r <- merge(r, map_unique[,c("ENSEMBL","ENTREZID","SYMBOL")], by="ENSEMBL")
  r <- r[is.finite(r$stat),]
  r <- r[order(abs(r$stat), decreasing=TRUE),]
  r <- r[!duplicated(r$ENTREZID),]
  ranks <- r$stat; names(ranks) <- r$ENTREZID; ranks <- sort(ranks, decreasing=TRUE)
  fg <- fgseaMultilevel(pathways=gmt, stats=ranks, minSize=15, maxSize=500, eps=0)
  fg$contrast <- nm
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
  fg <- as.data.frame(fg)
  gsea_all[[nm]] <- fg
  rank_audit[[nm]] <- data.frame(contrast=nm, input_tested_genes=nrow(x), valid_Wald_statistics=sum(is.finite(x$stat)), unique_Entrez_ranked=length(ranks), excluded_ambiguous_Ensembl=length(ambiguous), stringsAsFactors=FALSE)
}
gsea <- do.call(rbind, gsea_all)
gsea$padj_global_all_147 <- p.adjust(gsea$pval, method="BH")
gsea <- gsea[order(gsea$contrast, gsea$padj, -abs(gsea$NES)),]
write_csv(gsea, "Hallmark_GSEA_all.csv")
write_csv(gsea[gsea$padj<0.05,], "Hallmark_GSEA_significant_per_contrast.csv")
write_csv(do.call(rbind,rank_audit), "Hallmark_GSEA_rank_audit.csv")

summary <- list(
  source_sha256=input_audit$raw_sha256,
  reference=input_audit$reference,
  transcript_rows=input_audit$transcripts,
  exact_versioned_transcript_matches=input_audit$exact_versioned_matches,
  aggregated_genes=input_audit$gene_count,
  tested_genes=nrow(dds),
  design="~ replicate + time",
  reference_time="0 h untreated/pre-IL-17A baseline",
  contrasts=c("3 h vs 0 h","12 h vs 0 h","24 h vs 0 h"),
  DEG_summary=deg_summary,
  LRT_FDR_lt_0.05=sum(lrt$padj<0.05,na.rm=TRUE),
  PCA_variance=as.list((pc$sdev^2/sum(pc$sdev^2))[1:3]),
  correlation_range_off_diagonal=range(cormat[upper.tri(cormat)]),
  Hallmark_significant_per_contrast=as.list(table(gsea$contrast[gsea$padj<0.05])),
  Hallmark_tests_total=nrow(gsea),
  Hallmark_significant_global_all=as.list(table(gsea$contrast[gsea$padj_global_all_147<0.05])),
  limitations=input_audit$limitations
)
jsonlite::write_json(summary, file.path(out,"analysis_summary.json"), pretty=TRUE, auto_unbox=TRUE, digits=10)

writeLines(c(
  "# ADSC IL-17A time-course RNA-seq analysis",
  "",
  "Input: Salmon estimated transcript counts; TPM was not used.",
  "Annotation: all 232,117 versioned transcript IDs exactly match GENCODE v36 CHR (GRCh38.p13); counts were aggregated to 60,660 genes before rounding, with sample totals conserved.",
  "Filter: rounded gene counts >=10 in at least 3 samples; 14,425 genes tested.",
  "Design: DESeq2 ~ replicate + time. Replicate IDs 1-3 block the three complete biological experiments. Time reference is 0 h untreated/pre-IL-17A baseline.",
  "Contrasts: 3 h, 12 h, and 24 h versus 0 h. Positive log2FoldChange means higher at the indicated time point.",
  "DEG rule: BH FDR <0.05, no log2FC cutoff; independent filtering disabled; no LFC shrinkage.",
  "Omnibus test: DESeq2 LRT compares ~ replicate + time against ~ replicate.",
  "Hallmark GSEA: Wald statistic ranking; MSigDB 2025.1 Hs Hallmark Entrez GMT; unique unambiguous mappings; gene sets 15-500; fgseaMultilevel; 49 sets per contrast, with per-contrast BH FDR plus a joint 147-test BH column.",
  "",
  "Limitations:",
  "- Original quant.sf, tx2gene and effective transcript lengths were unavailable. Counts-only gene aggregation lacks sample-specific effective-length offsets; inference is provisional.",
  "- 0 h is the pre-treatment baseline. There is no time-matched untreated control, so IL-17A exposure and culture-time effects cannot be completely separated.",
  "- Estimated transcript counts were aggregated before rounding. TPM was not read or used.",
  "- Gene annotation was reconstructed using the exact matching GENCODE v36 reference; the original quantification index metadata was not supplied."
), file.path(out,"methods_and_limitations.md"), useBytes=TRUE)

writeLines(capture.output(sessionInfo()), file.path(out,"sessionInfo.txt"), useBytes=TRUE)
cat(jsonlite::toJSON(summary, pretty=TRUE, auto_unbox=TRUE, digits=5))
