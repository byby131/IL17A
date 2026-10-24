#!/usr/bin/env Rscript
# Reconstructed plotting implementation for the current Figure 5A–E.
# It consumes the archived panel source tables without rerunning statistics.
# The original plotting source was not present in the archived delivery.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

root <- normalizePath(".", mustWork = TRUE)
stopifnot(file.exists(file.path(root, "data/shared/original_counts_input.csv")))
out <- file.path(root, "reproduced_panels")
dir.create(out, showWarnings = FALSE)
samples <- c("NC_1", "NC_2", "NC_3", "IL17_1", "IL17_2", "IL17_3")
sample_labels <- c("Vehicle 1", "Vehicle 2", "Vehicle 3", "IL-17A 1", "IL-17A 2", "IL-17A 3")

save_panel <- function(p, name, width_mm, height_mm) {
  w <- width_mm / 25.4; h <- height_mm / 25.4
  if (capabilities("aqua")) {
    grDevices::quartz(type = "pdf", file = file.path(out, paste0(name, ".pdf")),
                      width = w, height = h, family = "Arial")
  } else {
    grDevices::cairo_pdf(file.path(out, paste0(name, ".pdf")), width = w, height = h,
                         family = "Arial")
  }
  print(p); grDevices::dev.off()
  svglite::svglite(file.path(out, paste0(name, ".svg")), width = w, height = h,
                   system_fonts = list(sans = "Arial"))
  print(p); grDevices::dev.off()
  ragg::agg_png(file.path(out, paste0(name, ".png")), width = w, height = h,
                units = "in", res = 300, background = "white")
  print(p); grDevices::dev.off()
  ragg::agg_tiff(file.path(out, paste0(name, ".tiff")), width = w, height = h,
                 units = "in", res = 600, compression = "lzw", background = "white")
  print(p); grDevices::dev.off()
}

theme_pub <- theme_classic(base_family = "Arial", base_size = 9) +
  theme(plot.title = element_text(face = "bold", size = 10),
        axis.title = element_text(face = "bold", size = 9),
        axis.text = element_text(size = 8),
        legend.position = "none",
        plot.margin = margin(5, 5, 5, 5))

# A: all 17,882 plotted genes; the archived x/y clipping is display-only.
vol <- fread("data/A_volcano/plotted_points.csv")
labs <- fread("data/A_volcano/top_gene_labels.csv")
stopifnot(nrow(vol) == 17882L, all(labs$gene_id %in% vol$gene_id))
vol[, point_color := as.character(point_color)]
vlabels <- vol[match(labs$gene_id, vol$gene_id)]
pA <- ggplot(vol, aes(plot_x_log2FC, plot_y_minus_log10_FDR)) +
  geom_point(aes(colour = point_color), size = 0.33, alpha = 0.72) +
  scale_colour_identity() +
  geom_vline(xintercept = 0, colour = "#333333", linewidth = 0.25) +
  ggrepel::geom_text_repel(data = vlabels, aes(label = gene_id), family = "Arial",
                           fontface = "bold", size = 2.5, max.overlaps = Inf,
                           min.segment.length = 0, seed = 20261025,
                           box.padding = 0.12, point.padding = 0.05) +
  scale_x_continuous(limits = c(-8, 8)) +
  scale_y_continuous(limits = c(0, 66), breaks = c(0, 20, 40, 60)) +
  labs(x = "Log2 FC", y = "−Log10 FDR") + theme_pub
save_panel(pA, "Figure5A_volcano_frozen_replot", 80, 66)

read_heatmap <- function(panel) {
  mat <- fread(file.path("data", panel, "heatmap_VST_rowZ.csv"))
  stopifnot(identical(names(mat), c("gene_id", samples)))
  long <- melt(mat, id.vars = "gene_id", measure.vars = samples,
               variable.name = "sample", value.name = "z")
  long[, gene_row := match(gene_id, mat$gene_id)]
  long[, y := nrow(mat) + 1L - gene_row]
  long[, x := match(as.character(sample), samples)]
  list(mat = mat, long = long)
}

z_scale <- function(limit) {
  scale_fill_gradient2(low = "#4D78A7", mid = "#F7F7F7", high = "#C95B6C",
                       midpoint = 0, limits = c(-limit, limit), oob = scales::squish,
                       name = "VST row Z")
}

# B: all significant genes in the archived order; no re-selection/re-clustering.
b <- read_heatmap("B_DEG_heatmap")
b_ann <- fread("data/B_DEG_heatmap/row_annotations.csv")
stopifnot(nrow(b$mat) == nrow(b_ann), identical(b$mat$gene_id, b_ann$gene_id))
up_n <- sum(b_ann$direction == "up")
pB <- ggplot(b$long, aes(x, y, fill = z)) +
  geom_raster() + z_scale(2.5) +
  geom_hline(yintercept = nrow(b$mat) - up_n + 0.5, linewidth = 0.35, colour = "white") +
  scale_x_continuous(breaks = 1:6, labels = sample_labels, expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.5, 6.5), ylim = c(0.5, nrow(b$mat) + 0.5)) +
  labs(x = NULL, y = NULL, title = sprintf("Up (%d)     Down (%d)", up_n, nrow(b$mat) - up_n)) +
  theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                    axis.text.y = element_blank(), axis.ticks.y = element_blank(),
                    legend.position = "bottom", legend.title = element_text(size = 7))
save_panel(pB, "Figure5B_all_DEGs_frozen_replot", 82, 58)

fdr_label <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x == 0, "0*", formatC(x, format = "e", digits = 1)))
}

small_heatmap <- function(panel, title, with_programs = FALSE) {
  h <- read_heatmap(panel)
  ann <- fread(file.path("data", panel, "row_annotations.csv"))
  stopifnot(nrow(h$mat) == nrow(ann), identical(h$mat$gene_id, ann$gene_id))
  ann[, y := nrow(ann) + 1L - row_order]
  ann[, fdr_text := fdr_label(padj)]
  p <- ggplot(h$long, aes(x, y, fill = z)) + geom_tile(width = 0.94, height = 0.94) +
    z_scale(2) +
    geom_text(data = ann, aes(x = 7.2, y = y, label = fdr_text),
              inherit.aes = FALSE, hjust = 0.5, family = "Arial", size = 2.15) +
    scale_x_continuous(breaks = 1:6, labels = sample_labels, limits = c(0.5, if (with_programs) 9.5 else 8.0),
                       expand = c(0, 0)) +
    scale_y_continuous(breaks = ann$y, labels = ann$gene_id, expand = expansion(add = 0.45)) +
    labs(x = NULL, y = NULL, title = title) +
    theme_pub + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7),
                      axis.text.y = element_text(size = 7, face = "bold"),
                      axis.ticks.y = element_blank(), legend.position = "bottom",
                      legend.title = element_text(size = 7))
  if (with_programs) {
    key <- fread("data/C_functional_program_heatmap/program_color_key.csv")
    programs <- rbindlist(lapply(seq_len(nrow(ann)), function(i) {
      codes <- strsplit(ann$all_contributing_programs[i], ";", fixed = TRUE)[[1L]]
      data.table(x = 8.0 + match(codes, key$program_code) * 0.25,
                 y = ann$y[i], color = key$color[match(codes, key$program_code)])
    }))
    p <- p + geom_tile(data = programs, aes(x = x, y = y), inherit.aes = FALSE,
                       fill = programs$color, width = 0.18, height = 0.70)
  }
  p
}

pC <- small_heatmap("C_functional_program_heatmap", "Functional programs", TRUE)
save_panel(pC, "Figure5C_programs_frozen_replot", 116, 115)
pD <- small_heatmap("D_adipogenic_TF_heatmap", "Adipogenic transcription factors")
save_panel(pD, "Figure5D_TFs_frozen_replot", 103, 81)

# E: use the archived 20 selected pathways and their stored plotting lengths.
ora <- fread("data/E_KEGG/all_20_displayed_pathways.csv")
stopifnot(nrow(ora) == 20L, all(table(ora$direction) == 10L))
draw_ora <- function(direction, fill) {
  wanted_direction <- direction
  d <- copy(ora[ora[["direction"]] == wanted_direction])
  stopifnot(nrow(d) == 10L)
  d[, y := 11L - display_order_within_direction]
  ggplot(d) +
    geom_rect(aes(xmin = 0, xmax = bar_length_minus_log10_FDR,
                  ymin = y - 0.32, ymax = y + 0.32), fill = fill) +
    geom_text(aes(x = 0.06, y = y, label = pathway_name), hjust = 0,
              family = "Arial", fontface = "bold", size = 2.15) +
    scale_x_continuous(limits = c(0, 8), breaks = seq(0, 8, 2), expand = c(0, 0)) +
    scale_y_continuous(breaks = 1:10, labels = rep("", 10), expand = expansion(add = 0.55)) +
    labs(x = "−log10(FDR)", y = NULL, title = if (direction == "up") "KEGG Up" else "KEGG Down") +
    theme_pub + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
}
pE <- draw_ora("up", "#C76B7B") + draw_ora("down", "#7EA4BF") +
  plot_layout(ncol = 2)
save_panel(pE, "Figure5E_KEGG_frozen_replot", 165, 90)

message("Replots saved to: ", out)
