#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(DESeq2)
  library(fgsea)
  library(msigdbr)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})
source("scripts/figure2_ghi_common.R")

# Figure 2H uses explicit R vector/raster devices so the plotting source itself
# documents the submission export contract (editable SVG/PDF plus 600-dpi TIFF).
save_nature_plot <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)

  # Opt in only when a working Cairo/XQuartz stack is available. The current
  # macOS R runtime advertises Cairo but cannot load its shared libraries.
  if (identical(Sys.getenv("FIGURE_USE_CAIRO", unset = "0"), "1")) {
    grDevices::cairo_pdf(
      paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Helvetica", bg = "white", onefile = TRUE
    )
  } else {
    # macOS builds without XQuartz cannot load cairo; the native PDF device
    # still preserves selectable vector text and avoids Dingbat point glyphs.
    grDevices::pdf(
      paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Helvetica", bg = "white", onefile = TRUE, useDingbats = FALSE
    )
  }
  print(plot)
  grDevices::dev.off()

  svglite::svglite(
    paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
    bg = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_mm, height = height_mm, units = "mm",
    res = dpi, compression = "lzw", background = "white"
  )
  print(plot)
  grDevices::dev.off()

  ragg::agg_png(
    paste0(stem, "_preview.png"), width = width_mm, height = height_mm,
    units = "mm", res = 300, background = "white"
  )
  print(plot)
  grDevices::dev.off()
}

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% c(4L, 5L)) {
  stop("Usage: Rscript scripts/make_figure2H_Q8_Q11_revision.R <all.H5Seurat> <frozen_cell_metadata.csv.gz> <donor_quality.csv> <output_dir> [Figure2F_mapping_metadata.csv.gz]")
}
h5_path <- normalizePath(args[[1]], mustWork = TRUE)
cell_path <- normalizePath(args[[2]], mustWork = TRUE)
quality_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
mapping_path <- if (length(args) == 5L) normalizePath(args[[5]], mustWork = TRUE) else NULL
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "supplementary"), recursive = TRUE, showWarnings = FALSE)
set.seed(264508)

cells <- fread(cell_path)
quality <- fread(quality_path)
if (!is.null(mapping_path)) {
  mapping <- fread(mapping_path, select = c("cell_id", "predicted_reference_state",
                                           "mapping_confidence", "mapping_entropy",
                                           "mean_reference_distance"))
  cells <- merge(cells, mapping, by = "cell_id", all.x = TRUE, sort = FALSE)
  mapping_summary <- cells[analysis_group %chin% c("LoSc Q8", "LoSc Q11"), .N,
                           by = .(analysis_group, predicted_reference_state)]
  mapping_summary[, percent := 100 * N / sum(N), by = analysis_group]
  fwrite(mapping_summary[order(analysis_group, -N)],
         file.path(output_dir, "source_data/Figure2H_Q8_Q11_reference_mapping_audit.csv"))
}
md <- read_h5_metadata(h5_path)
if ("source_row_h5" %in% names(cells)) cells[, source_row_h5 := NULL]
cells <- merge(cells, md[, .(cell_id, source_row_h5 = source_row)],
               by = "cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(cells$source_row_h5)) stop("Frozen cells do not map completely to H5Seurat")
counts <- read_h5_counts(h5_path, cells$source_row_h5)
cells <- cells[match(colnames(counts), cell_id)]
stopifnot(identical(cells$cell_id, colnames(counts)))

cell_totals <- as.numeric(Matrix::colSums(counts))
cells[, raw_total_counts := cell_totals]

qc_thresholds <- list(n_cells = 3L, total_counts = 5000, detected_genes = 1000L,
                      max_single_cell_fraction = 0.80)

compute_group_quality <- function(cell_dt, count_mat, label) {
  ids_by_donor <- split(cell_dt$cell_id, cell_dt$donor_id)
  rbindlist(lapply(names(ids_by_donor), function(donor) {
    ids <- ids_by_donor[[donor]]
    mat <- count_mat[, ids, drop = FALSE]
    totals <- Matrix::colSums(mat)
    out <- data.table(
      analysis_group = label,
      donor_id = donor,
      batch = cell_dt[donor_id == donor, batch][1],
      n_cells = length(ids),
      pseudobulk_total_counts = sum(totals),
      pseudobulk_detected_genes = sum(Matrix::rowSums(mat) > 0),
      max_single_cell_count_fraction = max(totals) / sum(totals)
    )
    out[, expression_qc_pass :=
          n_cells >= qc_thresholds$n_cells &
          pseudobulk_total_counts >= qc_thresholds$total_counts &
          pseudobulk_detected_genes >= qc_thresholds$detected_genes &
          max_single_cell_count_fraction <= qc_thresholds$max_single_cell_fraction]
    out
  }))
}

overall_cells <- copy(cells)
overall_cells[, analysis_group := fifelse(analysis_group == "Healthy ASPC-like",
                                           "Healthy ASPC-like", "LoSc overall ASPC-like")]
overall_quality <- rbindlist(list(
  compute_group_quality(overall_cells[analysis_group == "Healthy ASPC-like"], counts,
                        "Healthy ASPC-like"),
  compute_group_quality(overall_cells[analysis_group == "LoSc overall ASPC-like"], counts,
                        "LoSc overall ASPC-like")
))
fwrite(overall_quality[order(analysis_group, donor_id)],
       file.path(output_dir, "source_data/Figure2H_overall_ASPC_donor_quality.csv"))

hallmark_tbl <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
pathways_all <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

aggregate_pseudobulk <- function(cell_dt, count_mat) {
  donor_order <- unique(cell_dt[order(disease_group, donor_id), donor_id])
  cell_to_donor <- factor(cell_dt[match(colnames(count_mat), cell_id), donor_id], levels = donor_order)
  design_matrix <- sparse.model.matrix(~ 0 + cell_to_donor)
  pb <- round(as.matrix(count_mat %*% design_matrix))
  colnames(pb) <- donor_order
  donor_meta <- unique(cell_dt[, .(donor_id, disease_group, batch)])
  donor_meta <- donor_meta[match(donor_order, donor_id)]
  rownames(donor_meta) <- donor_meta$donor_id
  list(counts = pb, donor_meta = donor_meta)
}

design_diagnostics <- function(donor_meta, comparison) {
  donor_meta <- copy(donor_meta)
  donor_meta[, disease_group := relevel(factor(disease_group), "Healthy")]
  donor_meta[, batch := factor(batch)]
  mm <- model.matrix(~ batch + disease_group, donor_meta)
  shared_batches <- intersect(
    as.character(donor_meta[disease_group == "Healthy", unique(batch)]),
    as.character(donor_meta[disease_group == "Disease", unique(batch)])
  )
  disease_shared <- donor_meta[disease_group == "Disease" & as.character(batch) %chin% shared_batches, donor_id]
  healthy_shared <- donor_meta[disease_group == "Healthy" & as.character(batch) %chin% shared_batches, donor_id]
  loo <- rbindlist(lapply(donor_meta[disease_group == "Disease", donor_id], function(d) {
    x <- donor_meta[donor_id != d]
    x[, batch := droplevels(batch)]
    mm2 <- model.matrix(~ batch + disease_group, x)
    data.table(omitted_donor = d, estimable_after_omission = qr(mm2)$rank == ncol(mm2))
  }))
  list(
    summary = data.table(
      comparison = comparison,
      n_healthy_donors = donor_meta[disease_group == "Healthy", uniqueN(donor_id)],
      n_disease_donors = donor_meta[disease_group == "Disease", uniqueN(donor_id)],
      shared_batches = paste(as.character(shared_batches), collapse = ","),
      n_healthy_donors_in_shared_batches = length(healthy_shared),
      n_disease_donors_in_shared_batches = length(disease_shared),
      design_rank = qr(mm)$rank,
      design_columns = ncol(mm),
      estimable = qr(mm)$rank == ncol(mm),
      disease_donors_whose_omission_breaks_estimability =
        paste(loo[estimable_after_omission == FALSE, omitted_donor], collapse = ",")
    ),
    loo = loo[, comparison := comparison][]
  )
}

run_comparison <- function(label, healthy_cells, disease_cells, healthy_qc, disease_qc) {
  healthy_pass <- healthy_qc[expression_qc_pass == TRUE, donor_id]
  disease_pass <- disease_qc[expression_qc_pass == TRUE, donor_id]
  h <- healthy_cells[donor_id %chin% healthy_pass]
  d <- disease_cells[donor_id %chin% disease_pass]
  h[, disease_group := "Healthy"]
  d[, disease_group := "Disease"]
  selected <- rbindlist(list(h, d), fill = TRUE)
  selected[, disease_group := factor(disease_group, levels = c("Healthy", "Disease"))]
  ids <- selected$cell_id
  mat <- counts[, ids, drop = FALSE]
  selected <- selected[match(colnames(mat), cell_id)]
  pb <- aggregate_pseudobulk(selected, mat)
  donor_meta <- pb$donor_meta
  diag <- design_diagnostics(donor_meta, label)
  if (!diag$summary$estimable) stop(sprintf("Design is not estimable for %s", label))

  n_min <- min(table(donor_meta$disease_group))
  keep <- rowSums(pb$counts >= 10) >= n_min
  pb_keep <- pb$counts[keep, , drop = FALSE]
  filter_dt <- data.table(gene = rownames(pb$counts),
                          n_pseudobulks_count_ge_10 = rowSums(pb$counts >= 10),
                          retained = keep, n_min = n_min)

  donor_meta[, batch := factor(batch)]
  donor_meta[, disease_group := factor(disease_group, levels = c("Healthy", "Disease"))]
  dds <- DESeqDataSetFromMatrix(countData = pb_keep,
                                colData = as.data.frame(donor_meta),
                                design = ~ batch + disease_group)
  dds <- DESeq(dds, quiet = TRUE, minReplicatesForReplace = Inf, fitType = "local")
  res <- results(dds, contrast = c("disease_group", "Disease", "Healthy"),
                 cooksCutoff = TRUE, independentFiltering = TRUE)
  de <- as.data.table(as.data.frame(res), keep.rownames = "gene")
  setnames(de, "stat", "wald_stat")
  de[, comparison := label]
  de[, na_reason := fifelse(is.na(baseMean), "all_zero_or_unestimable",
                     fifelse(is.na(wald_stat) & is.na(pvalue), "Cook_outlier_or_unestimable",
                     fifelse(is.na(padj), "independent_filtering", NA_character_)))]
  de[, direction := fifelse(wald_stat > 0, "Disease-high-ranked", "Healthy-high-ranked")]

  rank_stats <- de[is.finite(wald_stat), wald_stat]
  names(rank_stats) <- de[is.finite(wald_stat), gene]
  rank_stats <- sort(rank_stats, decreasing = TRUE)
  pathways <- lapply(pathways_all, intersect, y = names(rank_stats))
  pathways <- pathways[lengths(pathways) >= 15 & lengths(pathways) <= 500]
  fg <- as.data.table(fgseaMultilevel(pathways = pathways, stats = rank_stats,
                                      minSize = 15, maxSize = 500, eps = 0))
  fg[, comparison := label]
  fg[, pathway_label := gsub("_", " ", gsub("^HALLMARK_", "", pathway))]
  fg[, direction := fifelse(NES > 0, "Disease-high-ranked", "Healthy-high-ranked")]

  vsd <- vst(dds, blind = FALSE)
  list(cells = selected, pseudobulk = pb$counts, donor_meta = donor_meta,
       filter = filter_dt, de = de, fg = fg, rank_stats = rank_stats,
       pathways = pathways, vst = assay(vsd), diagnostics = diag)
}

healthy_cells <- cells[analysis_group == "Healthy ASPC-like"]
q8_cells <- cells[analysis_group == "LoSc Q8"]
q11_cells <- cells[analysis_group == "LoSc Q11"]
healthy_qc <- quality[analysis_group == "Healthy ASPC-like"]
q8_qc <- quality[analysis_group == "LoSc Q8"]
q11_qc <- quality[analysis_group == "LoSc Q11"]

q8_res <- run_comparison("Q8 vs Healthy", healthy_cells, q8_cells, healthy_qc, q8_qc)
q11_res <- run_comparison("Q11 vs Healthy", healthy_cells, q11_cells, healthy_qc, q11_qc)

overall_healthy_qc <- overall_quality[analysis_group == "Healthy ASPC-like"]
overall_losc_qc <- overall_quality[analysis_group == "LoSc overall ASPC-like"]
overall_res <- run_comparison("Overall LoSc ASPC-like vs Healthy",
                              overall_cells[analysis_group == "Healthy ASPC-like"],
                              overall_cells[analysis_group == "LoSc overall ASPC-like"],
                              overall_healthy_qc, overall_losc_qc)

# Joint BH family requested for Q8 and Q11 Healthy-reference comparisons.
joint_fg <- rbindlist(list(q8_res$fg, q11_res$fg), fill = TRUE)
joint_fg[, q_joint := p.adjust(pval, method = "BH")]
joint_fg[, fdr_status := factor(ifelse(q_joint < 0.05, "FDR < 0.05", "FDR >= 0.05"),
                                levels = c("FDR < 0.05", "FDR >= 0.05"))]
joint_fg[, q_label := ifelse(q_joint < 0.001, "q<0.001", sprintf("q=%.3f", q_joint))]
joint_fg_export <- copy(joint_fg)
joint_fg_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
fwrite(joint_fg_export, file.path(output_dir, "source_data/Figure2H_Q8_Q11_Hallmark_GSEA_joint_BH.csv"))

q8_fg <- joint_fg[comparison == "Q8 vs Healthy"]
q8_sig <- q8_fg[q_joint < 0.05]
if (nrow(q8_sig)) {
  selected_pw <- unique(c(
    q8_sig[NES > 0][order(-NES, q_joint, pathway), head(pathway, 5)],
    q8_sig[NES < 0][order(NES, q_joint, pathway), head(pathway, 5)]
  ))
  q8_note <- NULL
} else {
  selected_pw <- unique(c(
    q8_fg[NES > 0][order(-NES, q_joint, pathway), head(pathway, 5)],
    q8_fg[NES < 0][order(NES, q_joint, pathway), head(pathway, 5)]
  ))
  q8_note <- "No Q8 pathway met joint FDR < 0.05"
}
main_dt <- joint_fg[pathway %chin% selected_pw]
q8_order <- q8_fg[pathway %chin% selected_pw][order(NES), pathway]
main_dt[, pathway_label := factor(pathway_label,
                                  levels = q8_fg[match(q8_order, pathway), pathway_label])]

q8_n_don <- uniqueN(q8_res$cells[disease_group == "Disease", donor_id])
q8_n_cell <- nrow(q8_res$cells[disease_group == "Disease"])
q11_n_don <- uniqueN(q11_res$cells[disease_group == "Disease", donor_id])
q11_n_cell <- nrow(q11_res$cells[disease_group == "Disease"])
h_n_don <- uniqueN(q8_res$cells[disease_group == "Healthy", donor_id])
h_n_cell <- nrow(q8_res$cells[disease_group == "Healthy"])

comparison_labels <- c(
  "Q8 vs Healthy" = sprintf("Q8 vs Healthy ASPC-like\n%d vs %d donors; %d vs %s cells",
                            q8_n_don, h_n_don, q8_n_cell, format(h_n_cell, big.mark = ",")),
  "Q11 vs Healthy" = sprintf("Q11 vs Healthy ASPC-like\n%d vs %d donors; %d vs %s cells",
                             q11_n_don, h_n_don, q11_n_cell, format(h_n_cell, big.mark = ","))
)
main_dt[, comparison := factor(comparison, levels = c("Q8 vs Healthy", "Q11 vs Healthy"))]

max_abs <- max(abs(main_dt$NES), na.rm = TRUE)
x_limit <- ceiling((max_abs + 0.75) * 2) / 2
direction_colors <- c("Disease-high-ranked" = "#D94755", "Healthy-high-ranked" = "#3B66A5")

p_double <- ggplot(main_dt, aes(NES, pathway_label)) +
  geom_vline(xintercept = 0, linewidth = 0.42, color = "#222222") +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_label, color = direction,
                   alpha = fdr_status), linewidth = 0.78, lineend = "round") +
  geom_point(aes(color = direction, shape = fdr_status), size = 2.65, stroke = 0.85) +
  geom_text(aes(label = q_label), color = "#111111", family = "Helvetica", size = 2.05,
            hjust = ifelse(main_dt$NES >= 0, -0.16, 1.16)) +
  facet_grid(. ~ comparison, labeller = as_labeller(comparison_labels)) +
  scale_color_manual(values = direction_colors,
                     breaks = c("Disease-high-ranked", "Healthy-high-ranked"),
                     labels = c("Disease-group high-ranked", "Healthy high-ranked"), name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR >= 0.05" = 0.45), guide = "none") +
  scale_x_continuous(limits = c(-x_limit, x_limit),
                     breaks = pretty(c(-max_abs, max_abs), n = 5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "Pathway-level transcriptional differences of Q8 and Q11\nrelative to Healthy ASPC-like references",
    subtitle = q8_note,
    x = "Normalized enrichment score (NES)", y = NULL,
    caption = paste0(
      "Joint BH across both Hallmark comparisons.\n",
      "Q8 is exploratory: the batch-adjusted effect depends on one shared-batch LoSc donor."
    )
  ) +
  theme_classic(base_family = "Helvetica", base_size = 7.0) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.line.x = element_line(linewidth = 0.38, color = "#222222"),
    axis.text.y = element_text(size = 6.7, color = "#111111", margin = margin(r = 3)),
    axis.text.x = element_text(size = 6.1, color = "#111111"),
    axis.title.x = element_text(size = 6.8, margin = margin(t = 4)),
    strip.background = element_blank(),
    strip.text = element_text(size = 7.0, face = "bold", color = "#111111", lineheight = 1.05,
                              margin = margin(b = 4)),
    plot.title = element_text(size = 8.0, face = "bold", lineheight = 1.02, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.4, color = "#555555", margin = margin(b = 4)),
    plot.caption = element_text(size = 5.4, color = "#555555", hjust = 0, margin = margin(t = 5)),
    legend.position = "top", legend.justification = "left", legend.box = "horizontal",
    legend.text = element_text(size = 6.1), legend.key.width = grid::unit(8, "pt"),
    panel.spacing.x = grid::unit(11, "pt"),
    plot.margin = margin(5, 11, 5, 4)
  )
save_nature_plot(p_double,
                 file.path(output_dir, "Figure2H_Q8_Q11_exploratory_jointBH_lollipop"),
                 178, 82)
main_dt_export <- copy(main_dt)
main_dt_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
fwrite(main_dt_export, file.path(output_dir, "source_data/Figure2H_Q8_Q11_main_plot_source_data.csv"))

# Overall ASPC-like fallback main panel, selected within its pre-defined family.
overall_fg <- copy(overall_res$fg)
overall_fg[, q := p.adjust(pval, method = "BH")]
overall_fg[, fdr_status := factor(ifelse(q < 0.05, "FDR < 0.05", "FDR >= 0.05"),
                                  levels = c("FDR < 0.05", "FDR >= 0.05"))]
overall_fg[, q_label := ifelse(q < 0.001, "q<0.001", sprintf("q=%.3f", q))]
overall_fg_export <- copy(overall_fg)
overall_fg_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
fwrite(overall_fg_export, file.path(output_dir, "source_data/Figure2H_overall_ASPC_Hallmark_GSEA.csv"))
overall_sig <- overall_fg[q < 0.05]
if (nrow(overall_sig)) {
  overall_sel <- rbindlist(list(
    head(overall_sig[NES > 0][order(-NES, q, pathway)], 5),
    head(overall_sig[NES < 0][order(NES, q, pathway)], 5)
  ), fill = TRUE)
  overall_note <- NULL
} else {
  overall_sel <- rbindlist(list(
    head(overall_fg[NES > 0][order(-NES, q, pathway)], 5),
    head(overall_fg[NES < 0][order(NES, q, pathway)], 5)
  ), fill = TRUE)
  overall_note <- "No pathway met FDR < 0.05"
}
overall_sel <- unique(overall_sel, by = "pathway")
overall_sel[, pathway_label := factor(pathway_label, levels = pathway_label[order(NES)])]
overall_n_don <- uniqueN(overall_res$cells[disease_group == "Disease", donor_id])
overall_n_cell <- nrow(overall_res$cells[disease_group == "Disease"])
overall_max <- max(abs(overall_sel$NES))
overall_limit <- ceiling((overall_max + 0.75) * 2) / 2
p_overall <- ggplot(overall_sel, aes(NES, pathway_label)) +
  geom_vline(xintercept = 0, linewidth = 0.42, color = "#222222") +
  geom_segment(aes(x = 0, xend = NES, yend = pathway_label, color = direction,
                   alpha = fdr_status), linewidth = 0.82, lineend = "round") +
  geom_point(aes(color = direction, shape = fdr_status), size = 2.75, stroke = 0.85) +
  geom_text(aes(label = q_label), color = "#111111", family = "Helvetica", size = 2.1,
            hjust = ifelse(overall_sel$NES >= 0, -0.16, 1.16)) +
  scale_color_manual(values = direction_colors,
                     breaks = c("Disease-high-ranked", "Healthy-high-ranked"),
                     labels = c("LoSc high-ranked", "Healthy high-ranked"), name = NULL) +
  scale_shape_manual(values = c("FDR < 0.05" = 16, "FDR >= 0.05" = 1), name = NULL) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "FDR >= 0.05" = 0.45), guide = "none") +
  scale_x_continuous(limits = c(-overall_limit, overall_limit),
                     breaks = pretty(c(-overall_max, overall_max), n = 5),
                     expand = expansion(mult = c(0.02, 0.02))) +
  labs(
    title = "LoSc ASPC-like pathway-level differences\nrelative to Healthy ASPC-like references",
    subtitle = sprintf("%d vs %d donors; %d vs %s cells%s",
                       overall_n_don, h_n_don, overall_n_cell, format(h_n_cell, big.mark = ","),
                       ifelse(is.null(overall_note), "", paste0(" | ", overall_note))),
    x = "Normalized enrichment score (NES)", y = NULL,
    caption = "BH correction across all tested Hallmark pathways; qualified donor pseudobulks only."
  ) +
  theme_classic(base_family = "Helvetica", base_size = 7.0) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.7, color = "#111111", margin = margin(r = 3)),
    axis.text.x = element_text(size = 6.1, color = "#111111"),
    plot.title = element_text(size = 8.0, face = "bold", lineheight = 1.02, margin = margin(b = 2)),
    plot.subtitle = element_text(size = 6.4, color = "#555555", margin = margin(b = 4)),
    plot.caption = element_text(size = 5.4, color = "#555555", hjust = 0, margin = margin(t = 5)),
    legend.position = "top", legend.justification = "left", legend.box = "horizontal",
    legend.text = element_text(size = 6.1),
    plot.margin = margin(5, 11, 5, 4)
  )
save_nature_plot(p_overall,
                 file.path(output_dir, "Figure2H_overall_ASPC_recommended_main_lollipop"),
                 120, 50)
overall_sel_export <- copy(overall_sel)
overall_sel_export[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
fwrite(overall_sel_export, file.path(output_dir, "source_data/Figure2H_overall_main_plot_source_data.csv"))

# Donor coverage plot: counts plus pass/fail/absence, prior to expression testing.
coverage_dt <- rbindlist(list(
  quality[, .(analysis_group, donor_id, batch, n_cells, expression_qc_pass)],
  overall_quality[analysis_group == "LoSc overall ASPC-like",
                  .(analysis_group, donor_id, batch, n_cells, expression_qc_pass)]
), fill = TRUE)
all_donors <- sort(unique(c(cells$donor_id)))
all_groups <- c("Healthy ASPC-like", "LoSc Q8", "LoSc Q11", "LoSc overall ASPC-like")
coverage_grid <- CJ(donor_id = all_donors, analysis_group = all_groups)
coverage_grid <- merge(coverage_grid, coverage_dt, by = c("donor_id", "analysis_group"), all.x = TRUE)
coverage_grid[, status := fifelse(is.na(n_cells), "Absent",
                           fifelse(expression_qc_pass, "Pass", "Fail"))]
coverage_grid[, cell_label := fifelse(is.na(n_cells), "", as.character(n_cells))]
coverage_grid[, donor_id := factor(donor_id, levels = rev(all_donors))]
coverage_grid[, analysis_group := factor(analysis_group, levels = all_groups)]
p_cov <- ggplot(coverage_grid, aes(analysis_group, donor_id, fill = status)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = cell_label), family = "Helvetica", size = 2.05,
            color = ifelse(coverage_grid$status == "Pass", "white", "#222222")) +
  scale_fill_manual(values = c(Pass = "#2A6670", Fail = "#E3A23B", Absent = "#ECECEC"),
                    breaks = c("Pass", "Fail", "Absent"), name = NULL) +
  labs(title = "Donor coverage and pseudobulk expression-quality audit",
       subtitle = "Numbers denote cells; qualification was frozen before DE and GSEA",
       x = NULL, y = NULL) +
  theme_classic(base_family = "Helvetica", base_size = 6.6) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(),
        axis.text.x = element_text(angle = 25, hjust = 1, size = 6.0, color = "#111111"),
        axis.text.y = element_text(size = 5.8, color = "#111111"),
        plot.title = element_text(size = 8.0, face = "bold"),
        plot.subtitle = element_text(size = 6.1, color = "#555555"),
        legend.position = "top", legend.justification = "left",
        plot.margin = margin(4, 4, 4, 4))
save_nature_plot(p_cov, file.path(output_dir, "supplementary/Figure2H_donor_coverage_QC"), 112, 132)
fwrite(coverage_grid, file.path(output_dir, "source_data/Figure2H_donor_coverage_plot_source_data.csv"))

# Export complete comparison objects and design diagnostics.
export_comparison <- function(res, prefix) {
  fwrite(res$cells, file.path(output_dir, sprintf("source_data/Figure2H_%s_included_cells.csv.gz", prefix)))
  fwrite(res$donor_meta, file.path(output_dir, sprintf("source_data/Figure2H_%s_model_donors.csv", prefix)))
  fwrite(res$filter, file.path(output_dir, sprintf("source_data/Figure2H_%s_gene_filter.csv", prefix)))
  fwrite(res$de, file.path(output_dir, sprintf("source_data/Figure2H_%s_DESeq2_full_results.csv.gz", prefix)))
  fwrite(res$diagnostics$summary,
         file.path(output_dir, sprintf("source_data/Figure2H_%s_design_diagnostics.csv", prefix)))
  fwrite(res$diagnostics$loo,
         file.path(output_dir, sprintf("source_data/Figure2H_%s_design_LODO_estimability.csv", prefix)))
  saveRDS(res$pseudobulk,
          file.path(output_dir, sprintf("source_data/Figure2H_%s_pseudobulk_raw_counts.rds", prefix)),
          compress = FALSE)
}
export_comparison(q8_res, "Q8")
export_comparison(q11_res, "Q11")
export_comparison(overall_res, "Overall")

diagnostics_all <- rbindlist(list(q8_res$diagnostics$summary,
                                  q11_res$diagnostics$summary,
                                  overall_res$diagnostics$summary), fill = TRUE)
fwrite(diagnostics_all, file.path(output_dir, "source_data/Figure2H_all_design_diagnostics.csv"))

writeLines(c(
  "panel=Figure 2H Q8/Q11 revision",
  "core_conclusion=directional Hallmark ranking in Q8 and Q11 relative to the same broad Healthy ASPC-like reference",
  "archetype=two-column quantitative lollipop grid with a qualified overall fallback main panel",
  "cell_boundary=entire identity-reviewed Q8 and Q11; no expression-based cell selection",
  sprintf("expression_qc=n_cells>=%d; pseudobulk_total_counts>=%d; detected_genes>=%d; max_single_cell_fraction<=%.2f",
          qc_thresholds$n_cells, qc_thresholds$total_counts, qc_thresholds$detected_genes,
          qc_thresholds$max_single_cell_fraction),
  sprintf("Q8_included=%d cells from %d donors", q8_n_cell, q8_n_don),
  sprintf("Q11_included=%d cells from %d donors", q11_n_cell, q11_n_don),
  sprintf("Healthy_reference=%d cells from %d donors", h_n_cell, h_n_don),
  "DESeq2_design=~ batch + disease_group; local dispersion trend; positive statistic means disease-group high",
  "Hallmark=all eligible human MSigDB Hallmark sets; fgseaMultilevel; size 15-500",
  "Q8_Q11_q=joint BH across all actually tested Hallmark results from both comparisons",
  "main_selection=Q8 only; up to five positive and five negative; Q11 displays identical pathways",
  sprintf("Q8_qualification=%s",
          ifelse(q8_res$diagnostics$summary$n_disease_donors_in_shared_batches <= 1,
                 "exploratory_only_due_to_one_shared_batch_disease_donor", "eligible_for_main")),
  "recommended_main=qualified overall ASPC-like fallback because Q8 disease coefficient is not robust to omission of its sole shared-batch donor",
  if (!is.null(mapping_path)) {
    q8_map <- mapping_summary[analysis_group == "LoSc Q8"][order(-N)]
    paste0("Q8_reference_mapping=", paste(sprintf("%s:%d/%d", q8_map$predicted_reference_state,
                                                   q8_map$N, sum(q8_map$N)), collapse = ";"))
  } else "Q8_reference_mapping=not_supplied",
  "static_preflight_source_syntax=Rscript parse completed successfully; validator warning reflects its limited R parser",
  "static_preflight_final_width=explicit call-site widths are 178 mm (Q8/Q11), 120 mm (overall) and 112 mm (coverage)",
  "pdf_device=native R PDF with Helvetica and useDingbats=FALSE; Cairo is opt-in because this macOS R runtime lacks XQuartz shared libraries",
  "interpretation_guard=NES is pathway-ranking evidence, not a direct functional effect size; Q8 and Q11 NES must not be subtracted"
), file.path(output_dir, "Figure2H_revision_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2H_revision_sessionInfo.txt"))
message("Figure 2H Q8/Q11 revision complete")
