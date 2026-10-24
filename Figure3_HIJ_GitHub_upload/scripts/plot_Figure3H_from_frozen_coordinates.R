#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(data.table); library(ggplot2); library(svglite); library(ragg)})
project_dir <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_Dose_PCA")
scores <- data.table::fread(file.path(output_dir, "source_data/ADSC_IL17A_Dose_PCA_coordinates.csv"))
group_levels <- c("Control", "Dose 10", "Dose 50", "Dose 100", "Dose 200")
scores[, group := factor(group, levels = group_levels)]
explained <- c(9.528, 7.993) # Frozen PCA variance percentages, not recomputed from coordinates.
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

