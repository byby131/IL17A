#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

options(stringsAsFactors = FALSE)

project_dir <- normalizePath(".", mustWork = TRUE)
counts_path <- file.path(project_dir, "data/private/ADSC_IL17A_Dose_RNAseq_expression_matrix_counts.csv")
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_Dose_PCA")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(counts_path)
if (!all(file.exists(required_files))) {
  stop("Missing input file: ", paste(required_files[!file.exists(required_files)], collapse = ", "))
}

counts_dt <- data.table::fread(counts_path, check.names = FALSE)
if (ncol(counts_dt) < 3L) stop("Counts matrix must contain one feature column and at least two samples.")

feature_id <- as.character(counts_dt[[1L]])
if (anyDuplicated(feature_id)) stop("Duplicate feature identifiers detected in the counts matrix.")

sample_names <- names(counts_dt)[-1L]
counts_mat <- as.matrix(counts_dt[, -1L, with = FALSE])
storage.mode(counts_mat) <- "double"
rownames(counts_mat) <- feature_id
colnames(counts_mat) <- sample_names

if (any(!is.finite(counts_mat))) stop("Counts matrix contains non-finite values.")
if (any(counts_mat < 0)) stop("Counts matrix contains negative values.")

parse_group <- function(x) {
  prefix <- sub("_[0-9]+$", "", x)
  out <- ifelse(
    prefix == "nc", "Control",
    ifelse(grepl("^s[0-9]+$", prefix), paste("Dose", sub("^s", "", prefix)), NA_character_)
  )
  out
}

group_label <- parse_group(sample_names)
replicate_id <- suppressWarnings(as.integer(sub("^.*_", "", sample_names)))
if (anyNA(group_label) || anyNA(replicate_id)) {
  stop("Sample names do not match the expected pattern: nc_1 or s10_1.")
}

group_levels <- c("Control", "Dose 10", "Dose 50", "Dose 100", "Dose 200")
if (!setequal(unique(group_label), group_levels)) {
  stop("Unexpected experimental groups: ", paste(unique(group_label), collapse = ", "))
}

sample_meta <- data.table(
  sample = sample_names,
  group = factor(group_label, levels = group_levels),
  replicate = replicate_id,
  dose_label = factor(group_label, levels = group_levels),
  dose_value = c("Control" = 0, "Dose 10" = 10, "Dose 50" = 50,
                 "Dose 100" = 100, "Dose 200" = 200)[group_label]
)

# Exploratory PCA workflow for estimated counts:
# retain transcripts with CPM > 1 in at least three samples, apply TMM library-size
# normalization, transform to log2-CPM, and calculate PCA on the 5,000 most
# variable retained transcripts. Fractional estimated counts are not rounded.
dge <- edgeR::DGEList(counts = counts_mat)
raw_cpm <- edgeR::cpm(dge, log = FALSE)
keep <- rowSums(raw_cpm > 1) >= 3L
if (sum(keep) < 100L) stop("Fewer than 100 transcripts pass the CPM filter.")

dge <- dge[keep, , keep.lib.sizes = FALSE]
dge <- edgeR::calcNormFactors(dge, method = "TMM")
log_cpm <- edgeR::cpm(dge, log = TRUE, prior.count = 2)

feature_variance <- apply(log_cpm, 1L, stats::var)
feature_variance[!is.finite(feature_variance)] <- NA_real_
valid_variance <- which(!is.na(feature_variance) & feature_variance > 0)
top_n <- min(5000L, length(valid_variance))
top_idx <- valid_variance[order(feature_variance[valid_variance], decreasing = TRUE)][seq_len(top_n)]

pca_fit <- stats::prcomp(t(log_cpm[top_idx, , drop = FALSE]), center = TRUE, scale. = FALSE)
explained <- 100 * (pca_fit$sdev^2 / sum(pca_fit$sdev^2))

scores <- data.table(
  sample = rownames(pca_fit$x),
  PC1 = pca_fit$x[, 1L],
  PC2 = pca_fit$x[, 2L]
)
scores <- merge(scores, sample_meta, by = "sample", sort = FALSE)
scores[, group := factor(group, levels = group_levels)]
setorder(scores, group, replicate)

hulls <- scores[, {
  idx <- grDevices::chull(PC1, PC2)
  .SD[idx]
}, by = group]

palette <- c(
  "Control" = "#ED8A32",
  "Dose 10" = "#CE9CB8",
  "Dose 50" = "#D43B7D",
  "Dose 100" = "#7FA873",
  "Dose 200" = "#2798B5"
)
shapes <- c(
  "Control" = 21,
  "Dose 10" = 22,
  "Dose 50" = 23,
  "Dose 100" = 24,
  "Dose 200" = 25
)

p <- ggplot(scores, aes(x = PC1, y = PC2)) +
  geom_hline(yintercept = 0, linewidth = 0.28, linetype = "dotted", color = "#C7C7C7") +
  geom_vline(xintercept = 0, linewidth = 0.28, linetype = "dotted", color = "#C7C7C7") +
  geom_polygon(
    data = hulls,
    aes(group = group, fill = group, color = group),
    alpha = 0.13,
    linewidth = 0.35,
    show.legend = FALSE
  ) +
  geom_point(
    aes(fill = group, shape = group),
    size = 3.2,
    stroke = 0.42,
    color = "#303030"
  ) +
  scale_fill_manual(values = palette, breaks = group_levels, name = NULL) +
  scale_color_manual(values = palette, breaks = group_levels, guide = "none") +
  scale_shape_manual(values = shapes, breaks = group_levels, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.09, 0.13))) +
  scale_y_continuous(expand = expansion(mult = c(0.10, 0.13))) +
  labs(
    title = "ADSC IL-17A dose-response transcriptome",
    x = sprintf("PC1 (%.1f%%)", explained[1L]),
    y = sprintf("PC2 (%.1f%%)", explained[2L])
  ) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(family = "Arial", size = 10, face = "plain", hjust = 0.5, margin = margin(b = 3.5)),
    axis.title = element_text(family = "Arial", size = 10, face = "plain", color = "black"),
    axis.text = element_text(family = "Arial", size = 10, color = "black"),
    axis.line = element_line(linewidth = 0.45, color = "black"),
    axis.ticks = element_line(linewidth = 0.40, color = "black"),
    axis.ticks.length = grid::unit(1.6, "mm"),
    legend.position = "inside",
    legend.position.inside = c(0.985, 0.985),
    legend.justification = c(1, 1),
    legend.direction = "vertical",
    legend.text = element_text(family = "Arial", size = 10, color = "black"),
    legend.key = element_rect(fill = scales::alpha("white", 0.82), color = NA),
    legend.key.height = grid::unit(4.0, "mm"),
    legend.spacing.y = grid::unit(0.2, "mm"),
    legend.background = element_rect(fill = scales::alpha("white", 0.82), color = NA),
    legend.margin = margin(1.5, 2, 1.5, 2, unit = "mm"),
    plot.margin = margin(4, 4, 3, 4, unit = "mm")
  ) +
  guides(
    fill = guide_legend(override.aes = list(
      shape = unname(shapes),
      fill = unname(palette),
      size = 2.7,
      color = "#303030",
      linewidth = 0.4
    )),
    shape = "none"
  )

width_mm <- 92
height_mm <- 84

pdf_path <- file.path(output_dir, "ADSC_IL17A_Dose_PCA_reference_palette.pdf")
svg_path <- file.path(output_dir, "ADSC_IL17A_Dose_PCA_reference_palette.svg")
tiff_path <- file.path(output_dir, "ADSC_IL17A_Dose_PCA_reference_palette.tiff")
png_path <- file.path(output_dir, "ADSC_IL17A_Dose_PCA_reference_palette_preview.png")

# Keep cairo_pdf explicitly available for portable vector export. Quartz is used on
# macOS so that the installed Arial font remains editable and correctly embedded.
if (isTRUE(getOption("pca.use_cairo_pdf", FALSE))) {
  grDevices::cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
} else if (capabilities("aqua")) {
  grDevices::quartz(type = "pdf", file = pdf_path, width = width_mm / 25.4,
                    height = height_mm / 25.4, family = "Arial")
} else {
  grDevices::pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4,
                 family = "Arial", useDingbats = FALSE)
}
print(p)
dev.off()

svglite::svglite(svg_path, width = width_mm / 25.4, height = height_mm / 25.4,
                 system_fonts = list(sans = "Arial"))
print(p)
dev.off()

ragg::agg_tiff(tiff_path, width = width_mm / 25.4, height = height_mm / 25.4,
               units = "in", res = 600, compression = "lzw")
print(p)
dev.off()

ragg::agg_png(png_path, width = width_mm / 25.4, height = height_mm / 25.4,
              units = "in", res = 300, background = "white")
print(p)
dev.off()

data.table::fwrite(scores, file.path(source_dir, "ADSC_IL17A_Dose_PCA_coordinates.csv"))
data.table::fwrite(sample_meta, file.path(source_dir, "ADSC_IL17A_Dose_sample_metadata.csv"))

loadings <- data.table(
  feature_id = rownames(pca_fit$rotation),
  PC1_loading = pca_fit$rotation[, 1L],
  PC2_loading = pca_fit$rotation[, 2L],
  variance_logCPM = feature_variance[top_idx]
)
loadings[, max_abs_loading := pmax(abs(PC1_loading), abs(PC2_loading))]
setorder(loadings, -max_abs_loading)
data.table::fwrite(loadings, file.path(source_dir, "ADSC_IL17A_Dose_PCA_loadings.csv"))

qa_lines <- c(
  "Figure: ADSC IL-17A dose-response PCA",
  sprintf("Input: %s", counts_path),
  sprintf("Input features: %s", format(nrow(counts_mat), big.mark = ",")),
  sprintf("Input samples: %d", ncol(counts_mat)),
  sprintf("Retained after CPM > 1 in at least 3 samples: %s", format(sum(keep), big.mark = ",")),
  sprintf("PCA features: %s most variable transcripts", format(top_n, big.mark = ",")),
  "Normalization: edgeR TMM; transformation: log2-CPM with prior.count = 2",
  "PCA: centered, unscaled prcomp on samples; no sample was excluded",
  sprintf("PC1 explained variance: %.3f%%", explained[1L]),
  sprintf("PC2 explained variance: %.3f%%", explained[2L]),
  "Replicate numbers and dose groups were parsed from sample names; biological-versus-technical replicate status is not asserted.",
  "Dose units were not provided and are intentionally omitted from labels.",
  "The input contains fractional estimated counts; values were not rounded.",
  "Hulls connect the three samples in each group and do not represent confidence intervals.",
  "Palette revision only: groups follow the supplied reference order (orange, dusty pink, magenta, sage green, cyan-blue).",
  "PCA calculation, sample coordinates, shapes, hull geometry, and labels were not changed.",
  sprintf("Final figure size: %d x %d mm", width_mm, height_mm),
  "Font: Arial, 10 pt",
  "Exports: editable-text PDF/SVG plus 600-dpi TIFF and 300-dpi PNG preview"
)
writeLines(qa_lines, file.path(output_dir, "ADSC_IL17A_Dose_PCA_QA_notes.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "ADSC_IL17A_Dose_PCA_sessionInfo.txt"))

message("PCA outputs written to: ", output_dir)
