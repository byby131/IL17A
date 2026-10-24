#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/make_figure2B_new_strategy.R <metadata_umap.csv.gz> <output_dir>")
}

metadata_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_family <- "Helvetica"
set.seed(264508)

group_order <- c(
  "Keratinocytes",
  "Fibroblasts",
  "ADSC/ASPC-like stromal cells",
  "Endothelial cells",
  "Pericytes/VSMCs",
  "T cells",
  "Myeloid cells",
  "Other immune cells",
  "Other skin cells"
)

# Palette reconstructed from the visual grammar of Steele et al., Nature
# Immunology (2025), Fig. 3. Assignments are fixed across the revised Figure 2.
group_palette <- c(
  "Keratinocytes" = "#FFC17F",
  "Fibroblasts" = "#9A7FBD",
  "ADSC/ASPC-like stromal cells" = "#6FB3A8",
  "Endothelial cells" = "#67A8CD",
  "Pericytes/VSMCs" = "#3581B7",
  "T cells" = "#50AA4B",
  "Myeloid cells" = "#F36569",
  "Other immune cells" = "#704BA3",
  "Other skin cells" = "#CF9F88"
)

cluster_to_group <- data.table(
  seurat_clusters = as.character(0:42),
  global_group = c(
    "Keratinocytes", "Keratinocytes", "Keratinocytes", "Fibroblasts",
    "Endothelial cells", "Myeloid cells", "Endothelial cells", "ADSC/ASPC-like stromal cells",
    "Other immune cells", "T cells", "Keratinocytes", "Pericytes/VSMCs",
    "Pericytes/VSMCs", "Fibroblasts", "Keratinocytes", "Keratinocytes",
    "Pericytes/VSMCs", "Endothelial cells", "Pericytes/VSMCs", "Fibroblasts",
    "T cells", "Keratinocytes", "Keratinocytes", "Endothelial cells",
    "Myeloid cells", "Keratinocytes", "Endothelial cells", "Endothelial cells",
    "Keratinocytes", "Other skin cells", "Pericytes/VSMCs", "Fibroblasts",
    "Other skin cells", "Endothelial cells", "Other skin cells", "Other immune cells",
    "Other skin cells", "Keratinocytes", "Other immune cells", "Keratinocytes",
    "Keratinocytes", "Fibroblasts", "Keratinocytes"
  )
)

group_definition <- data.table(
  global_group = group_order,
  included_public_labels = c(
    "Keratinocytes",
    "Fibroblasts excluding cluster 7",
    "Fibroblast cluster 7; PI16/DPP4/CFD-enriched",
    "Endothelial Cells + Lymphatic Endothelial Cells",
    "Pericytes + Smooth Muscle Cells",
    "T Cells",
    "Macrophage + Dendritic Cells",
    "B Cells + NK Cells + Mast Cells",
    "Eccrine Cells + Melanocytes + marker-corrected Schwann cluster 32"
  )
)

md <- fread(metadata_path)
required_columns <- c("cell_id", "seurat_clusters", "library_id", "health", "UMAP_1", "UMAP_2")
missing_columns <- setdiff(required_columns, names(md))
if (length(missing_columns) > 0L) {
  stop("Missing metadata columns: ", paste(missing_columns, collapse = ", "))
}

md[, seurat_clusters := as.character(seurat_clusters)]
md[, input_row := .I]
md <- merge(md, cluster_to_group, by = "seurat_clusters", all.x = TRUE, sort = FALSE)
if (anyNA(md$global_group)) {
  stop("At least one cluster was not assigned to a revised global group")
}
setorder(md, input_row)
md[, input_row := NULL]
md[, global_group := factor(global_group, levels = group_order)]

group_counts <- md[, .(
  n_cells = .N,
  n_donors = uniqueN(library_id),
  healthy_cells = sum(health == "Healthy"),
  losc_cells = sum(health != "Healthy")
), by = global_group]
group_counts[, fraction_percent := 100 * n_cells / sum(n_cells)]

# Draw the largest classes first so small populations remain visible.
draw_order <- group_counts[order(-n_cells), as.character(global_group)]
md[, draw_rank := match(as.character(global_group), draw_order)]
setorder(md, draw_rank)

label_positions <- data.table(
  global_group = factor(group_order, levels = group_order),
  label = c(
    "Keratinocytes",
    "Fibroblasts",
    "ADSC/ASPC-like\nstromal cells",
    "Endothelial cells",
    "Pericytes/VSMCs",
    "T cells",
    "Myeloid cells",
    "Other immune cells",
    "Other skin cells"
  ),
  UMAP_1 = c(-3.4, 7.7, 7.15, 4.1, 9.55, -10.7, -12.45, -9.0, 0.0),
  UMAP_2 = c(2.55, 7.85, 10.55, -10.25, -0.75, -3.9, 4.2, -6.55, 6.9)
)

x_range <- range(md$UMAP_1)
y_range <- range(md$UMAP_2)
x_span <- diff(x_range)
y_span <- diff(y_range)
axis_x0 <- x_range[1] + 0.018 * x_span
axis_y0 <- y_range[1] + 0.022 * y_span
axis_len_x <- 0.090 * x_span
axis_len_y <- 0.090 * y_span

p_umap <- ggplot(md, aes(UMAP_1, UMAP_2, color = global_group)) +
  geom_point(size = 0.19, alpha = 0.95, stroke = 0) +
  geom_text(
    data = label_positions,
    aes(label = label),
    color = "#111111",
    family = "Arial",
    fontface = "bold",
    size = 8.5 / ggplot2::.pt,
    lineheight = 0.95,
    show.legend = FALSE
  ) +
  annotate(
    "segment",
    x = axis_x0, xend = axis_x0 + axis_len_x,
    y = axis_y0, yend = axis_y0,
    linewidth = 0.34, color = "#22272A",
    arrow = arrow(length = unit(1.15, "mm"), type = "closed")
  ) +
  annotate(
    "segment",
    x = axis_x0, xend = axis_x0,
    y = axis_y0, yend = axis_y0 + axis_len_y,
    linewidth = 0.34, color = "#22272A",
    arrow = arrow(length = unit(1.15, "mm"), type = "closed")
  ) +
  annotate(
    "text",
    x = axis_x0 + axis_len_x * 0.52,
    y = axis_y0 - 0.035 * y_span,
    label = "UMAP1", family = base_family, size = 2.20, color = "#111111"
  ) +
  annotate(
    "text",
    x = axis_x0 - 0.035 * x_span,
    y = axis_y0 + axis_len_y * 0.52,
    label = "UMAP2", angle = 90, family = base_family, size = 2.20, color = "#111111"
  ) +
  scale_color_manual(values = group_palette, drop = FALSE) +
  coord_equal(
    xlim = c(x_range[1] - 0.7, x_range[2] + 0.7),
    ylim = c(y_range[1] - 0.8, y_range[2] + 0.8),
    clip = "off"
  ) +
  labs(title = "Global skin atlas") +
  theme_void(base_family = base_family) +
  theme(
    plot.title = element_text(
      size = 8.2, face = "plain", color = "#111111",
      hjust = 0.5, margin = margin(b = 2.2)
    ),
    plot.margin = margin(5.0, 4.5, 4.0, 5.0),
    legend.position = "none"
  ) +
  annotate(
    "text",
    x = -Inf, y = Inf, label = "B",
    hjust = -0.15, vjust = 1.15,
    family = base_family, fontface = "bold", size = 3.2, color = "black"
  )

final_figure_width_mm <- 183
width_mm <- 105
height_mm <- 92
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4
stem <- file.path(output_dir, "Figure2B_global_UMAP_new_strategy")

svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
print(p_umap)
dev.off()

if (identical(Sys.info()[["sysname"]], "Darwin")) {
  quartz(
    type = "pdf", file = paste0(stem, ".pdf"),
    width = width_in, height = height_in, family = base_family, bg = "white"
  )
} else {
  cairo_pdf(
    paste0(stem, ".pdf"), width = width_in, height = height_in,
    family = base_family, bg = "white"
  )
}
print(p_umap)
dev.off()

agg_tiff(
  paste0(stem, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white", compression = "lzw"
)
print(p_umap)
dev.off()

agg_png(
  paste0(stem, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 300, background = "white"
)
print(p_umap)
dev.off()

fwrite(md[, .(cell_id, seurat_clusters, library_id, health, UMAP_1, UMAP_2, global_group)],
       file.path(output_dir, "Figure2B_global_UMAP_source_data.csv.gz"))
fwrite(group_counts, file.path(output_dir, "Figure2B_global_group_counts.csv"))
fwrite(group_definition, file.path(output_dir, "Figure2B_global_group_definition.csv"))
fwrite(
  data.table(global_group = group_order, hex = unname(group_palette)),
  file.path(output_dir, "Figure2B_global_palette.csv")
)

qa_lines <- c(
  "panel=Figure 2B",
  "core_conclusion=GSE264508 resolves nine manuscript-level global skin compartments, including an ADSC/ASPC-like stromal compartment.",
  "backend=R only",
  sprintf("dimensions_mm=%s x %s", width_mm, height_mm),
  sprintf("final_composite_width_mm=%s", final_figure_width_mm),
  sprintf("cells=%s", nrow(md)),
  sprintf("donors=%s", uniqueN(md$library_id)),
  "sampling=all cells shown; no display downsampling",
  "main_panel_scope=one global UMAP only",
  "supplementary_scope=Control/LoSc split UMAP, global marker dot plot, and per-sample UMAP",
  "annotation=43 published clusters merged into nine prespecified major classes; cluster 7 retained as PI16/DPP4/CFD-enriched ADSC/ASPC-like stroma",
  "statistics=no hypothesis test is shown in Figure 2B",
  "vector=SVG and PDF; all cell points and text remain vector elements",
  "palette=exact colors selected from the supplied UMAP reference card; non-fluorescent blue/teal/green/coral/purple/brown-orange family",
  "label_style=all UMAP group labels are Arial 8.5 pt bold black text without boxes or background fill",
  "integrity=no cells excluded and no coordinates altered",
  "preflight_source_syntax=Rscript parse passed",
  "preflight_width_note=105 mm is the individual panel width within a planned 183 mm two-column composite"
)
writeLines(qa_lines, file.path(output_dir, "Figure2B_global_UMAP_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2B_global_UMAP_sessionInfo.txt"))
