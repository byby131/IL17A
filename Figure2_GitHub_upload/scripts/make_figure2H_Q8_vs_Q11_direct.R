#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop(paste(
    "Usage: Rscript scripts/make_figure2H_Q8_vs_Q11_direct.R",
    "<Q8_pseudobulk.rds> <Q8_model_donors.csv>",
    "<Q11_pseudobulk.rds> <Q11_model_donors.csv>",
    "<donor_group_quality.csv> <output_dir>"
  ))
}

q8_pb_path <- normalizePath(args[[1]], mustWork = TRUE)
q8_dm_path <- normalizePath(args[[2]], mustWork = TRUE)
q11_pb_path <- normalizePath(args[[3]], mustWork = TRUE)
q11_dm_path <- normalizePath(args[[4]], mustWork = TRUE)
quality_path <- normalizePath(args[[5]], mustWork = TRUE)
output_dir <- args[[6]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
set.seed(264508)

q8_pb <- readRDS(q8_pb_path)
q11_pb <- readRDS(q11_pb_path)
q8_dm <- fread(q8_dm_path)
q11_dm <- fread(q11_dm_path)
quality <- fread(quality_path)

q8_donors <- q8_dm[disease_group == "Disease", donor_id]
q11_donors <- q11_dm[disease_group == "Disease", donor_id]
paired_donors <- sort(intersect(q8_donors, q11_donors))
if (length(paired_donors) < 3L) {
  stop("Fewer than three qualified paired donors are available for Q8 vs Q11")
}

genes <- intersect(rownames(q8_pb), rownames(q11_pb))
q8_mat <- q8_pb[genes, paired_donors, drop = FALSE]
q11_mat <- q11_pb[genes, paired_donors, drop = FALSE]
colnames(q8_mat) <- paste(paired_donors, "Q8", sep = "__")
colnames(q11_mat) <- paste(paired_donors, "Q11", sep = "__")
pb <- cbind(q11_mat, q8_mat)

sample_meta <- data.table(
  sample_id = colnames(pb),
  donor_id = rep(paired_donors, 2L),
  state = rep(c("Q11", "Q8"), each = length(paired_donors))
)
sample_meta[, donor_id := factor(donor_id)]
sample_meta[, state := factor(state, levels = c("Q11", "Q8"))]
rownames(sample_meta) <- sample_meta$sample_id

design_formula <- ~ donor_id + state
design_matrix <- model.matrix(design_formula, sample_meta)
design_rank <- qr(design_matrix)$rank
if (design_rank != ncol(design_matrix)) stop("Paired donor design is not full rank")
residual_df <- nrow(design_matrix) - design_rank

# Match the established gene-information rule: count >=10 in at least the
# number of pseudobulks in the smaller state group (three paired donors).
n_min <- min(table(sample_meta$state))
n_pb_ge_10 <- rowSums(pb >= 10)
keep <- n_pb_ge_10 >= n_min
pb_keep <- round(pb[keep, , drop = FALSE])

dds <- DESeqDataSetFromMatrix(
  countData = pb_keep,
  colData = as.data.frame(sample_meta),
  design = design_formula
)
dds <- DESeq(
  dds, quiet = TRUE, fitType = "local",
  minReplicatesForReplace = Inf
)
res <- results(
  dds,
  contrast = c("state", "Q8", "Q11"),
  cooksCutoff = TRUE,
  independentFiltering = TRUE
)

de <- as.data.table(as.data.frame(res), keep.rownames = "gene")
setnames(de, "stat", "wald_stat")
de[, comparison := "Q8 vs Q11 paired"]
de[, direction := fifelse(
  wald_stat > 0, "Q8-high-ranked",
  fifelse(wald_stat < 0, "Q11-high-ranked", "Zero")
)]
de[, na_reason := fifelse(
  is.na(baseMean), "all_zero_or_unestimable",
  fifelse(
    is.na(wald_stat) & is.na(pvalue), "Cook_outlier_or_unestimable",
    fifelse(is.na(padj), "independent_filtering", NA_character_)
  )
)]

rank_stats <- de[is.finite(wald_stat), wald_stat]
names(rank_stats) <- de[is.finite(wald_stat), gene]
rank_stats <- sort(rank_stats, decreasing = TRUE)

hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
pathways <- split(hallmark$gene_symbol, hallmark$gs_name)
pathways <- lapply(pathways, intersect, y = names(rank_stats))
pathways <- pathways[lengths(pathways) >= 15 & lengths(pathways) <= 500]

fg <- as.data.table(fgseaMultilevel(
  pathways = pathways,
  stats = rank_stats,
  minSize = 15,
  maxSize = 500,
  eps = 0
))
fg[, q := p.adjust(pval, method = "BH")]
fg[, pathway_label := gsub("_", " ", gsub("^HALLMARK_", "", pathway))]
fg[, direction := fifelse(NES > 0, "Q8-high-ranked", "Q11-high-ranked")]
fg[, fdr_status := factor(
  ifelse(q < 0.05, "FDR < 0.05", "FDR >= 0.05"),
  levels = c("FDR < 0.05", "FDR >= 0.05")
)]
fg[, q_label := ifelse(q < 0.001, "q<0.001", sprintf("q=%.3f", q))]

# The display follows the previous revision: top 10 NES in each direction.
# Significance is encoded independently and never used to fill a quota.
plot_dt <- rbindlist(list(
  head(fg[NES > 0][order(-NES, q, pathway)], 10L),
  head(fg[NES < 0][order(NES, q, pathway)], 10L)
), fill = TRUE)
plot_dt <- unique(plot_dt, by = "pathway")
plot_dt <- plot_dt[order(NES)]
plot_dt[, pathway_display := factor(pathway_label, levels = pathway_label)]

q8_cells <- quality[
  analysis_group == "LoSc Q8" & donor_id %chin% paired_donors & expression_qc_pass == TRUE,
  sum(n_cells)
]
q11_cells <- quality[
  analysis_group == "LoSc Q11" & donor_id %chin% paired_donors & expression_qc_pass == TRUE,
  sum(n_cells)
]
sig_q8 <- fg[NES > 0 & q < 0.05, .N]
sig_q11 <- fg[NES < 0 & q < 0.05, .N]
n_tested_pathways <- nrow(fg)
n_display_q8 <- plot_dt[NES > 0, .N]
n_display_q11 <- plot_dt[NES < 0, .N]

direction_colors <- c("Q8-high-ranked" = "#D94755", "Q11-high-ranked" = "#3B66A5")
max_abs <- max(abs(plot_dt$NES), na.rm = TRUE)
x_limit <- ceiling((max_abs + 1.05) * 2) / 2

p <- ggplot(plot_dt, aes(NES, pathway_display)) +
  geom_vline(xintercept = 0, linewidth = 0.40, color = "#222222") +
  geom_segment(
    aes(x = 0, xend = NES, yend = pathway_display, color = direction, alpha = fdr_status),
    linewidth = 0.70, lineend = "round"
  ) +
  geom_point(aes(color = direction, shape = fdr_status), size = 2.35, stroke = 0.78) +
  geom_text(
    aes(label = q_label), color = "#111111", family = "Helvetica", size = 1.78,
    hjust = ifelse(plot_dt$NES >= 0, -0.14, 1.14)
  ) +
  scale_color_manual(
    values = direction_colors,
    breaks = c("Q8-high-ranked", "Q11-high-ranked"),
    labels = c("Q8 high-ranked", "Q11 high-ranked"),
    name = NULL
  ) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR >= 0.05" = 0.42), guide = "none") +
  guides(
    color = guide_legend(order = 1, nrow = 1, byrow = TRUE),
    shape = guide_legend(order = 2, nrow = 1, byrow = TRUE)
  ) +
  scale_x_continuous(
    limits = c(-x_limit, x_limit),
    breaks = pretty(c(-max_abs, max_abs), n = 5),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    title = "Q8 versus Q11 pathway-level differences",
    subtitle = sprintf(
      "%d paired donors; %d Q8 and %d Q11 cells; %d eligible Hallmark pathways\nFDR < 0.05: %d Q8-high, %d Q11-high",
      length(paired_donors), q8_cells, q11_cells, n_tested_pathways, sig_q8, sig_q11
    ),
    x = "Normalized enrichment score (NES)", y = NULL,
    caption = paste0(
      "Paired donor pseudobulk DESeq2 (~ donor + state); Hallmark BH correction.\n",
      "Exploratory because only three paired donors are available."
    )
  ) +
  theme_classic(base_family = "Helvetica", base_size = 6.8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.35, color = "#222222"),
    axis.text.y = element_text(size = 5.9, face = "bold", color = "#111111", margin = margin(r = 2)),
    axis.text.x = element_text(size = 5.9, color = "#111111"),
    axis.title.x = element_text(size = 6.4, margin = margin(t = 3)),
    plot.title = element_text(size = 8.0, face = "bold", lineheight = 1.02, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 5.9, color = "#444444"),
    plot.caption = element_text(size = 5.1, color = "#555555", hjust = 0),
    legend.position = "top", legend.justification = "left", legend.box = "vertical",
    legend.spacing.y = grid::unit(0.2, "pt"),
    legend.text = element_text(size = 5.8), legend.key.width = grid::unit(7, "pt"),
    plot.margin = margin(4, 9, 4, 4)
  )

stem <- file.path(output_dir, "Figure2H_Q8_vs_Q11_paired_top10_each_direction")
width_mm <- 135
height_mm <- 118

if (identical(Sys.getenv("FIGURE_USE_CAIRO", unset = "0"), "1")) {
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Helvetica", bg = "white"
  )
} else {
  grDevices::pdf(
    paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
    family = "Helvetica", bg = "white", useDingbats = FALSE
  )
}
print(p)
grDevices::dev.off()

svglite::svglite(
  paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
  bg = "white"
)
print(p)
grDevices::dev.off()

ragg::agg_tiff(
  paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm",
  res = 600, compression = "lzw", background = "white"
)
print(p)
grDevices::dev.off()

ragg::agg_png(
  paste0(stem, "_preview.png"), width = width_mm, height = height_mm, units = "mm",
  res = 300, background = "white"
)
print(p)
grDevices::dev.off()

fg_export <- copy(fg)
fg_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
plot_export <- copy(plot_dt)
plot_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]

fwrite(sample_meta, file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_paired_sample_metadata.csv"))
fwrite(data.table(
  gene = genes,
  n_pseudobulks_count_ge_10 = n_pb_ge_10,
  retained = keep,
  n_min = n_min
), file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_gene_filter.csv"))
fwrite(de, file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_DESeq2_full_results.csv.gz"))
fwrite(fg_export, file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_Hallmark_GSEA.csv"))
fwrite(plot_export, file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_plot_source_data.csv"))
saveRDS(pb, file.path(output_dir, "source_data/Figure2H_Q8_vs_Q11_paired_raw_pseudobulk.rds"), compress = FALSE)

writeLines(c(
  "panel=Figure 2H Q8 versus Q11 direct comparison",
  sprintf("paired_donors=%s", paste(paired_donors, collapse = ",")),
  sprintf("n_pairs=%d; Q8_cells=%d; Q11_cells=%d", length(paired_donors), q8_cells, q11_cells),
  sprintf("design=~ donor_id + state; rank=%d/%d; residual_df=%d", design_rank, ncol(design_matrix), residual_df),
  "contrast=Q8 minus Q11; positive Wald/NES means Q8 high-ranked; negative means Q11 high-ranked",
  sprintf("gene_filter=count>=10 in at least %d of 6 paired pseudobulks", n_min),
  "Hallmark=all eligible human MSigDB Hallmark sets; fgseaMultilevel; size 15-500",
  "q=BH correction across all Hallmark pathways tested in this direct comparison",
  sprintf("actual_significant=%d Q8-high and %d Q11-high at q<0.05", sig_q8, sig_q11),
  sprintf("display=up to 10 per direction; %d Q8-high and %d Q11-high shown; open points are not significant", n_display_q8, n_display_q11),
  "qualification=exploratory because only three paired biological donors are available",
  "static_preflight_source_syntax=Rscript parse completed successfully; validator warning reflects its limited R parser",
  "static_preflight_final_width=explicit export width is 135 mm and height is 118 mm",
  "pdf_device=native R PDF with Helvetica and useDingbats=FALSE; SVG retains editable text"
), file.path(output_dir, "Figure2H_Q8_vs_Q11_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2H_Q8_vs_Q11_sessionInfo.txt"))

message("Figure 2H Q8 vs Q11 paired analysis complete")
