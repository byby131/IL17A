#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ComplexHeatmap)
  library(circlize)
  library(svglite)
  library(ragg)
  library(grid)
})

options(stringsAsFactors = FALSE)

# ComplexHeatmap measures text on multiple grid devices. Register an Arial-named
# metric family for the base PDF/PostScript metric databases while the actual
# PDF/SVG/ragg exports continue to use the installed Arial TrueType font.
if (!("Arial" %in% names(grDevices::pdfFonts()))) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("sans")[[1L]])
}
if (!("Arial" %in% names(grDevices::postscriptFonts()))) {
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts("sans")[[1L]])
}
if (capabilities("aqua") && !("Arial" %in% names(grDevices::quartzFonts()))) {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

project_dir <- normalizePath(".", mustWork = TRUE)
counts_path <- file.path(project_dir, "data/private/ADSC_IL17A_Dose_RNAseq_expression_matrix_counts.csv")
annotation_path <- file.path(project_dir, "work/ADSC_IL17A_Dose_RNAseq/Ensembl_transcript_gene_mapping_current.csv")
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_0_vs_50_DEG_heatmaps")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(counts_path)) stop("Counts matrix is missing.")
if (!file.exists(annotation_path)) stop("Ensembl transcript-to-gene annotation is missing.")

counts_dt <- data.table::fread(counts_path, check.names = FALSE)
annotation_dt <- data.table::fread(annotation_path)

required_annotation_columns <- c(
  "ensembl_transcript_id", "ensembl_gene_id", "external_gene_name", "transcript_biotype"
)
if (!all(required_annotation_columns %in% names(annotation_dt))) {
  stop("Annotation file does not contain the required Ensembl fields.")
}

feature_id_versioned <- as.character(counts_dt[[1L]])
feature_id <- sub("[.][0-9]+$", "", feature_id_versioned)
sample_columns <- names(counts_dt)[-1L]
selected_samples <- c("nc_1", "nc_2", "nc_3", "s50_1", "s50_2", "s50_3")
if (!all(selected_samples %in% sample_columns)) stop("Required 0 or 50 ng/ml samples are missing.")

annotation_dt <- annotation_dt[!is.na(ensembl_transcript_id) & nzchar(ensembl_transcript_id)]
setkey(annotation_dt, ensembl_transcript_id)
matched_annotation <- annotation_dt[J(feature_id)]

mapped <- !is.na(matched_annotation$ensembl_gene_id) & nzchar(matched_annotation$ensembl_gene_id)
if (sum(mapped) < 10000L) stop("Insufficient transcript-to-gene annotation coverage.")

transcript_counts <- as.matrix(counts_dt[mapped, ..selected_samples])
storage.mode(transcript_counts) <- "double"
if (any(!is.finite(transcript_counts)) || any(transcript_counts < 0)) {
  stop("Selected estimated-count values must be finite and non-negative.")
}

gene_id_per_transcript <- matched_annotation$ensembl_gene_id[mapped]
gene_counts <- rowsum(transcript_counts, group = gene_id_per_transcript, reorder = FALSE)

gene_annotation <- unique(
  matched_annotation[mapped, .(ensembl_gene_id, external_gene_name)]
)
gene_annotation[, has_symbol := !is.na(external_gene_name) & nzchar(external_gene_name)]
setorder(gene_annotation, ensembl_gene_id, -has_symbol)
gene_annotation <- gene_annotation[!duplicated(ensembl_gene_id)]

sample_metadata <- data.table(
  sample = selected_samples,
  display = c("0-1", "0-2", "0-3", "50-1", "50-2", "50-3"),
  group = factor(
    c(rep("0 ng/ml", 3L), rep("50 ng/ml", 3L)),
    levels = c("0 ng/ml", "50 ng/ml")
  ),
  replicate = rep(seq_len(3L), 2L)
)

y <- edgeR::DGEList(
  counts = gene_counts[, selected_samples, drop = FALSE],
  group = sample_metadata$group
)
design <- model.matrix(~ group, data = sample_metadata)
keep_gene <- edgeR::filterByExpr(y, design = design)
if (sum(keep_gene) < 1000L) stop("Too few genes pass expression filtering.")

y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- edgeR::calcNormFactors(y, method = "TMM")
y <- edgeR::estimateDisp(y, design = design, robust = TRUE)
fit <- edgeR::glmQLFit(y, design = design, robust = TRUE)
qlf <- edgeR::glmQLFTest(fit, coef = "group50 ng/ml")

de_table <- data.table::as.data.table(
  edgeR::topTags(qlf, n = Inf, sort.by = "PValue")$table,
  keep.rownames = "ensembl_gene_id"
)
de_table <- merge(
  de_table,
  gene_annotation[, .(ensembl_gene_id, external_gene_name)],
  by = "ensembl_gene_id",
  all.x = TRUE,
  sort = FALSE
)
de_table[, gene_label := fifelse(
  !is.na(external_gene_name) & nzchar(external_gene_name),
  external_gene_name,
  ensembl_gene_id
)]
de_table[, direction := fifelse(
  FDR < 0.05 & logFC >= 1, "Higher in 50 ng/ml",
  fifelse(FDR < 0.05 & logFC <= -1, "Higher in 0 ng/ml", "Not selected")
)]

eligible <- de_table[direction != "Not selected"]
eligible[, abs_logFC := abs(logFC)]
setorder(eligible, direction, FDR, -abs_logFC)
if (eligible[direction == "Higher in 50 ng/ml", .N] < 25L ||
    eligible[direction == "Higher in 0 ng/ml", .N] < 25L) {
  stop("Fewer than 25 significant genes are available in one direction.")
}

select_balanced <- function(n_high_50, n_high_0) {
  up <- copy(eligible[direction == "Higher in 50 ng/ml"])
  down <- copy(eligible[direction == "Higher in 0 ng/ml"])
  setorder(up, FDR, -abs_logFC)
  setorder(down, FDR, -abs_logFC)
  rbind(up[seq_len(n_high_50)], down[seq_len(n_high_0)])
}

selected_top50 <- select_balanced(25L, 25L)
selected_top25 <- select_balanced(13L, 12L)

log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)

make_matrix_bundle <- function(selected_dt) {
  ids <- selected_dt$ensembl_gene_id
  mat_logcpm <- log_cpm[ids, selected_samples, drop = FALSE]
  mat_z <- t(scale(t(mat_logcpm)))
  mat_z[!is.finite(mat_z)] <- 0
  mat_z <- pmax(pmin(mat_z, 2), -2)

  # Display-only row-wise min-max scaling reproduces the reference heatmap
  # full red/blue contrast with six samples. Differential statistics and all
  # clustering distances continue to use the original normalized expression.
  row_min <- apply(mat_logcpm, 1L, min)
  row_max <- apply(mat_logcpm, 1L, max)
  row_range <- row_max - row_min
  row_range[row_range == 0 | !is.finite(row_range)] <- 1
  mat_display <- sweep(mat_logcpm, 1L, row_min, FUN = "-")
  mat_display <- sweep(mat_display, 1L, row_range, FUN = "/")
  mat_display <- (mat_display * 4) - 2
  mat_display[!is.finite(mat_display)] <- 0

  display_labels <- make.unique(selected_dt$gene_label, sep = " · ")
  rownames(mat_logcpm) <- display_labels
  rownames(mat_z) <- display_labels
  rownames(mat_display) <- display_labels
  colnames(mat_logcpm) <- sample_metadata$display
  colnames(mat_z) <- sample_metadata$display
  colnames(mat_display) <- sample_metadata$display

  row_split <- factor(
    selected_dt$direction,
    levels = c("Higher in 50 ng/ml", "Higher in 0 ng/ml"),
    labels = c("50 ng/ml high", "0 ng/ml high")
  )
  list(logcpm = mat_logcpm, z = mat_z, display = mat_display, row_split = row_split)
}

top50_bundle <- make_matrix_bundle(selected_top50)
top25_bundle <- make_matrix_bundle(selected_top25)

heatmap_colors <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#0000FF", "#FFFFFF", "#FF0000")
)

make_heatmap <- function(bundle, top_n) {
  cluster_rows_from_original <- function(m) {
    original <- bundle$z[rownames(m), , drop = FALSE]
    stats::hclust(stats::dist(original), method = "complete")
  }
  cluster_columns_from_original <- function(m) {
    sample_labels <- rownames(m)
    original <- t(bundle$z[, sample_labels, drop = FALSE])
    stats::hclust(stats::dist(original), method = "complete")
  }

  ComplexHeatmap::Heatmap(
    bundle$display,
    name = "Scaled expression",
    col = heatmap_colors,
    na_col = "#F7F7F7",
    rect_gp = grid::gpar(col = NA),
    border = TRUE,
    border_gp = grid::gpar(col = "#303030", lwd = 0.7),
    cluster_rows = cluster_rows_from_original,
    cluster_row_slices = FALSE,
    row_dend_side = "left",
    row_dend_width = grid::unit(15, "mm"),
    row_dend_gp = grid::gpar(col = "black", lwd = 0.7),
    row_split = bundle$row_split,
    row_gap = grid::unit(1.4, "mm"),
    row_title_side = "right",
    row_title_rot = 90,
    row_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "plain"),
    row_names_side = "right",
    row_names_gp = grid::gpar(fontfamily = "Arial", fontsize = 10, col = "black"),
    row_names_max_width = grid::unit(34, "mm"),
    cluster_columns = cluster_columns_from_original,
    cluster_column_slices = FALSE,
    column_dend_side = "top",
    column_dend_height = grid::unit(12, "mm"),
    column_dend_gp = grid::gpar(col = "black", lwd = 0.7),
    column_split = sample_metadata$group,
    column_gap = grid::unit(3, "mm"),
    column_title = c("0 ng/ml", "50 ng/ml"),
    column_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold"),
    column_names_gp = grid::gpar(fontfamily = "Arial", fontsize = 10, col = "black"),
    column_names_rot = 0,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(
      title = NULL,
      at = c(2, 1, 0, -1, -2),
      labels = c("2", "1", "0", "-1", "-2"),
      direction = "horizontal",
      legend_width = grid::unit(42, "mm"),
      labels_gp = grid::gpar(fontfamily = "Arial", fontsize = 10),
      border = "#404040"
    ),
    use_raster = FALSE
  )
}

draw_heatmap <- function(bundle, top_n) {
  ht <- make_heatmap(bundle, top_n)
  ComplexHeatmap::draw(
    ht,
    heatmap_legend_side = "top",
    annotation_legend_side = "top",
    merge_legends = TRUE,
    column_title = sprintf("Top %d DEGs: 0 vs 50 ng/ml", top_n),
    column_title_gp = grid::gpar(fontfamily = "Arial", fontsize = 10, fontface = "bold"),
    padding = grid::unit(c(2, 2, 2, 2), "mm")
  )
}

export_heatmap <- function(bundle, top_n, width_mm, height_mm) {
  prefix <- file.path(
    output_dir,
    sprintf("ADSC_IL17A_0_vs_50_Top%d_DEG_heatmap_contrast_scaled", top_n)
  )

  # Cairo PDF remains available for portable vector output. Quartz is preferred on
  # macOS so the installed Arial font is retained as editable text.
  if (isTRUE(getOption("deg.heatmap.use_cairo_pdf", FALSE))) {
    grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_mm / 25.4,
                         height = height_mm / 25.4, family = "Arial")
  } else if (capabilities("aqua")) {
    grDevices::quartz(type = "pdf", file = paste0(prefix, ".pdf"),
                      width = width_mm / 25.4, height = height_mm / 25.4,
                      family = "Arial")
  } else {
    grDevices::pdf(paste0(prefix, ".pdf"), width = width_mm / 25.4,
                   height = height_mm / 25.4, family = "Arial", useDingbats = FALSE)
  }
  draw_heatmap(bundle, top_n)
  dev.off()

  svglite::svglite(paste0(prefix, ".svg"), width = width_mm / 25.4,
                   height = height_mm / 25.4, system_fonts = list(sans = "Arial"))
  draw_heatmap(bundle, top_n)
  dev.off()

  ragg::agg_tiff(paste0(prefix, ".tiff"), width = width_mm / 25.4,
                 height = height_mm / 25.4, units = "in", res = 600,
                 compression = "lzw")
  draw_heatmap(bundle, top_n)
  dev.off()

  ragg::agg_png(paste0(prefix, "_preview.png"), width = width_mm / 25.4,
                height = height_mm / 25.4, units = "in", res = 300,
                background = "white")
  draw_heatmap(bundle, top_n)
  dev.off()
}

export_heatmap(top50_bundle, 50L, width_mm = 118, height_mm = 220)
export_heatmap(top25_bundle, 25L, width_mm = 118, height_mm = 142)

data.table::fwrite(de_table, file.path(source_dir, "ADSC_IL17A_0_vs_50_full_DEG_statistics.csv"))
data.table::fwrite(sample_metadata, file.path(source_dir, "ADSC_IL17A_0_vs_50_sample_metadata.csv"))
data.table::fwrite(
  data.table(
    sample = selected_samples,
    library_size = y$samples$lib.size,
    TMM_norm_factor = y$samples$norm.factors,
    group = sample_metadata$group
  ),
  file.path(source_dir, "ADSC_IL17A_0_vs_50_normalization_factors.csv")
)

write_selected_source <- function(selected_dt, bundle, top_n) {
  selected_out <- copy(selected_dt)
  selected_out[, heatmap_rank := seq_len(.N)]
  data.table::fwrite(
    selected_out,
    file.path(source_dir, sprintf("ADSC_IL17A_0_vs_50_Top%d_selected_DEGs.csv", top_n))
  )

  z_dt <- data.table::as.data.table(bundle$z, keep.rownames = "gene_label")
  logcpm_dt <- data.table::as.data.table(bundle$logcpm, keep.rownames = "gene_label")
  data.table::fwrite(
    z_dt,
    file.path(source_dir, sprintf("ADSC_IL17A_0_vs_50_Top%d_row_Zscores.csv", top_n))
  )
  display_dt <- data.table::as.data.table(bundle$display, keep.rownames = "gene_label")
  data.table::fwrite(
    display_dt,
    file.path(source_dir, sprintf("ADSC_IL17A_0_vs_50_Top%d_display_scaled.csv", top_n))
  )
  data.table::fwrite(
    logcpm_dt,
    file.path(source_dir, sprintf("ADSC_IL17A_0_vs_50_Top%d_logCPM.csv", top_n))
  )
}

write_selected_source(selected_top50, top50_bundle, 50L)
write_selected_source(selected_top25, top25_bundle, 25L)

qa_lines <- c(
  "Figure: ADSC IL-17A 0 ng/ml versus 50 ng/ml Top DEG heatmaps",
  "Role: comparison / transcriptional response discovery",
  "Input unit: transcript-level fractional estimated counts",
  sprintf("Input transcripts: %s", format(nrow(counts_dt), big.mark = ",")),
  sprintf("Mapped transcripts: %s", format(sum(mapped), big.mark = ",")),
  sprintf("Aggregated Ensembl genes: %s", format(nrow(gene_counts), big.mark = ",")),
  sprintf("Genes retained by edgeR filterByExpr: %s", format(sum(keep_gene), big.mark = ",")),
  "Groups: 0 ng/ml (nc_1, nc_2, nc_3) versus 50 ng/ml (s50_1, s50_2, s50_3)",
  "Differential analysis: edgeR TMM normalization, robust dispersion estimation, quasi-likelihood GLM",
  "Contrast: 50 ng/ml minus 0 ng/ml; positive logFC means higher in 50 ng/ml",
  "DEG eligibility: Benjamini-Hochberg FDR < 0.05 and absolute log2 fold change >= 1",
  sprintf("Eligible DEGs: %d total; %d higher at 50 ng/ml; %d higher at 0 ng/ml",
          eligible[, .N], eligible[direction == "Higher in 50 ng/ml", .N],
          eligible[direction == "Higher in 0 ng/ml", .N]),
  "Top 50 selection: 25 genes per direction, ordered by FDR then absolute logFC",
  "Top 25 selection: 13 genes higher at 50 ng/ml and 12 genes higher at 0 ng/ml; nested within the Top 50 rule",
  "Heatmap colors: display-only row-wise min-max scaling of normalized log2-CPM to [-2, 2]",
  "Rows are clustered within direction and columns within dose group using the original row Z-scores, not the display-scaled matrix",
  "No sample was excluded. Only the requested 0 and 50 ng/ml groups entered this comparison.",
  "Replicate status was inferred from filenames; biological-versus-technical replicate status was not provided.",
  "Ensembl transcript-to-gene mapping was obtained from the official Ensembl BioMart service on 2026-09-11.",
  "Font: Arial 10 pt; reference-matched pure red = row maximum and pure blue = row minimum",
  "The display stretch does not alter gene selection, differential statistics, normalized expression values, or clustering input.",
  "Top 50 size: 118 x 220 mm; Top 25 size: 118 x 142 mm",
  "Exports: editable-text PDF/SVG, 600-dpi TIFF, and 300-dpi PNG previews"
)
writeLines(qa_lines, file.path(output_dir, "ADSC_IL17A_0_vs_50_DEG_heatmaps_QA_notes.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "ADSC_IL17A_0_vs_50_DEG_heatmaps_sessionInfo.txt"))

message("DEG heatmaps written to: ", output_dir)
