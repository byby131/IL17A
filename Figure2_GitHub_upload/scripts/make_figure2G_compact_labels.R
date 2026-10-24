#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(svglite)
  library(ragg)
  library(scales)
})
suppressWarnings(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"))

# Figure contract:
# - Claim: donor-balanced marker expression supports the stromal identities of
#   Healthy references and LoSc Q1-Q12 clusters without consuming excessive
#   width in row labels or legends.
# - Archetype: compact quantitative grid with a side statistics table.
# - Data/statistics: layout-only revision of the frozen Figure 2G source data;
#   no values, rows, genes, clusters or identity calls are altered.
# - Export: 183 x 118 mm, Arial, editable SVG/PDF and 600-dpi TIFF.

input_dir <- "outputs/Figure2G_strict"
output_dir <- "outputs/Figure2G_compact_labels"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)

dot <- fread(file.path(input_dir, "Figure2G_dotplot_values_donor_balanced.csv"))
clusters <- fread(file.path(input_dir, "Figure2G_query_cluster_summary.csv"))

marker_modules <- list(
  `General stromal` = c("DCN", "LUM", "COL1A1"),
  `ASPC-like` = c("PDGFRA", "CD34", "PI16", "DPP4", "CFD", "APOD"),
  `Papillary FB` = c("COL13A1", "APCDD1", "WIF1"),
  `ECM remodeling` = c("POSTN", "THBS2", "SFRP4", "FN1"),
  `Pericyte` = c("RGS5", "CSPG4", "PDGFRB"),
  `VSMC` = c("MYH11", "CNN1", "TAGLN"),
  `Other-lineage checks` = c("PTPRC", "PECAM1", "EPCAM")
)
marker_genes <- unique(unlist(marker_modules, use.names = FALSE))

ref_full <- c(
  "Homeostatic/adipogenic ASPC",
  "Papillary/conventional fibroblast",
  "ECM-remodeling fibroblast",
  "Pericyte",
  "VSMC"
)
ref_short <- c(
  "Homeostatic/adipogenic ASPC" = "ASPC",
  "Papillary/conventional fibroblast" = "Papillary FB",
  "ECM-remodeling fibroblast" = "ECM-FB",
  "Pericyte" = "Pericyte",
  "VSMC" = "VSMC"
)
query_order <- paste0("Q", 1:12)
query_short <- setNames(query_order, query_order)
query_short[c("Q8", "Q11")] <- c("Q8 · ASPC-like", "Q11 · ASPC-like")
display_map <- c(ref_short, query_short)
top_to_bottom <- c(unname(ref_short[ref_full]), unname(query_short[query_order]))
row_levels <- rev(top_to_bottom)

plot_dt <- copy(dot)
plot_dt[, row_label := unname(display_map[expression_group])]
if (anyNA(plot_dt$row_label)) stop("Unmapped Figure 2G row labels")
plot_dt[, row_label := factor(row_label, levels = row_levels)]
plot_dt[, gene := factor(gene, levels = marker_genes)]
cap <- quantile(plot_dt$color_value, 0.98, na.rm = TRUE)
plot_dt[, color_plot := pmin(color_value, cap)]

base_dot <- ggplot(plot_dt, aes(gene, row_label)) +
  geom_point(
    aes(size = pct_detected_donor_balanced, fill = color_plot),
    shape = 21, colour = "#2A2A2A", stroke = 0.18
  ) +
  annotate(
    "segment", x = 0.5, xend = length(marker_genes) + 0.5,
    y = length(query_order) + 0.5, yend = length(query_order) + 0.5,
    linewidth = 0.45, colour = "#333333"
  ) +
  scale_size_continuous(
    range = c(0.15, 3.15), limits = c(0, 100), breaks = c(20, 50, 80),
    name = "Detected cells (%)"
  ) +
  scale_fill_stepsn(
    colours = c("#FFFFFF", "#F4D2C3", "#D4775B", "#8D251E"),
    limits = c(0, cap), n.breaks = 6, oob = squish,
    name = "Mean expression (log1p)"
  ) +
  scale_x_discrete(position = "top", expand = expansion(add = 0.30)) +
  scale_y_discrete(expand = expansion(add = 0.48)) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "#222222", fill = NA, linewidth = 0.40),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90, hjust = 0, vjust = 0.5, size = 10,
      colour = "#111111", margin = margin(b = 2)
    ),
    axis.text.y = element_text(size = 10, colour = "#111111", margin = margin(r = 2)),
    legend.title = element_text(size = 10, colour = "#111111"),
    legend.text = element_text(size = 10, colour = "#111111"),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.box.just = "left",
    legend.spacing.x = unit(3, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(1, 1, 0, 1)
  ) +
  guides(
    size = guide_legend(
      order = 1, title.position = "top", title.hjust = 0,
      direction = "horizontal", nrow = 1,
      override.aes = list(fill = "white", colour = "#2A2A2A")
    ),
    fill = guide_coloursteps(
      order = 2, title.position = "top", title.hjust = 0,
      direction = "horizontal", barwidth = unit(34, "mm"),
      barheight = unit(3.0, "mm"), ticks = FALSE
    )
  )

p_dot <- base_dot + theme(legend.position = "none")
legend_device <- ragg::agg_capture(width = 2000, height = 700, res = 300)
legend_grob <- cowplot::get_legend(base_dot)
invisible(legend_device())
p_legend <- wrap_elements(full = legend_grob)

ref_stats <- unique(dot[expression_group %chin% ref_full,
  .(expression_group, n_cells, n_donors)])
ref_stats[, `:=`(
  display_group = unname(ref_short[expression_group]),
  confidence = "Ref.", uncertain_percent = 0,
  disease_deviated_percent = 0
)]
query_stats <- clusters[, .(
  expression_group = query_cluster,
  n_cells, n_donors,
  display_group = unname(query_short[query_cluster]),
  confidence = fifelse(identity_confidence == "Moderate", "Mod.", identity_confidence),
  uncertain_percent, disease_deviated_percent
)]
side_dt <- rbindlist(list(ref_stats, query_stats), use.names = TRUE, fill = TRUE)
side_dt[, `:=`(n_cells = as.numeric(n_cells), n_donors = as.numeric(n_donors))]
side_dt[, display_group := factor(display_group, levels = row_levels)]

side_num <- melt(
  side_dt,
  id.vars = "display_group",
  measure.vars = c("n_cells", "n_donors", "uncertain_percent", "disease_deviated_percent"),
  variable.name = "metric", value.name = "value"
)
side_num[, label := fifelse(
  metric %chin% c("uncertain_percent", "disease_deviated_percent"),
  sprintf("%.1f", value), as.character(round(value))
)]
side_long <- rbindlist(list(
  side_num[, .(display_group, metric = as.character(metric), label)],
  side_dt[, .(display_group, metric = "confidence", label = confidence)]
))
side_long[, metric := factor(
  metric,
  levels = c("n_cells", "n_donors", "confidence", "uncertain_percent", "disease_deviated_percent"),
  labels = c("Cells", "Donors", "Conf.", "Unc.%", "Dev.%")
)]

p_side <- ggplot(side_long, aes(metric, display_group, label = label)) +
  geom_tile(fill = "white", colour = "#B8B8B8", linewidth = 0.23) +
  geom_text(family = "Arial", size = 10 / ggplot2::.pt, colour = "#111111") +
  scale_x_discrete(position = "top", expand = c(0, 0)) +
  scale_y_discrete(expand = expansion(add = 0.48)) +
  labs(x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    axis.text.x = element_text(
      size = 10, angle = 90, hjust = 0, vjust = 0.5,
      colour = "#111111", margin = margin(b = 2)
    ),
    plot.margin = margin(1, 1, 0, 0)
  )

caption_text <- paste0(
  "Figure 2G | Marker expression supporting stromal identities and reference-mapping outliers. ",
  "Dot size indicates the donor-balanced percentage of detected cells and colour indicates log1p-transformed ",
  "donor-balanced mean expression. Healthy reference abbreviations: ASPC, homeostatic/adipogenic ASPC; ",
  "Papillary FB, papillary/conventional fibroblast; ECM-FB, ECM-remodeling fibroblast. ",
  "Disease-cluster identities: Q1, VSMC-like; Q2, VSMC-like; Q3, Pericyte-like; Q4, Pericyte-like; ",
  "Q5, VSMC-like; Q6, provisional VSMC-like; Q7, provisional papillary/conventional fibroblast-like; ",
  "Q8, ASPC-like; Q9, provisional ECM-remodeling fibroblast-like; Q10, ECM-remodeling fibroblast-like; ",
  "Q11, ASPC-like; Q12, unresolved/possible non-stromal admixture. Cells, number of cells; Donors, number of donors; ",
  "Conf., identity confidence; Unc.%, uncertain fraction; Dev.%, disease-deviated fraction."
)

p_grid <- p_dot + p_side + plot_layout(widths = c(5.0, 1.75))
p_bottom <- p_legend + plot_spacer() + plot_layout(widths = c(5.0, 1.75))
p_final <- p_grid / p_bottom +
  plot_layout(heights = c(1, 0.18)) +
  plot_annotation(
    title = "Marker expression supporting stromal identities and reference-mapping outliers",
    theme = theme(
      text = element_text(family = "Arial", colour = "#111111"),
      plot.title = element_text(size = 10, face = "bold", hjust = 0, margin = margin(b = 1)),
      plot.margin = margin(1, 1, 1, 1)
    )
  )

stem <- file.path(output_dir, "Figure2G_compact_stromal_identity_dotplot")
width_mm <- 183
height_mm <- 105

svglite::svglite(paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(p_final)
dev.off()

if (isTRUE(getOption("figure2g.use_cairo_pdf", FALSE))) {
  cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
            family = "Arial", onefile = TRUE)
} else if (capabilities("aqua")) {
  quartz(type = "pdf", file = paste0(stem, ".pdf"), width = width_mm / 25.4,
         height = height_mm / 25.4, family = "Arial")
} else {
  pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial", useDingbats = FALSE)
}
print(p_final)
dev.off()

agg_tiff(paste0(stem, ".tiff"), width = width_mm / 25.4, height = height_mm / 25.4,
         units = "in", res = 600, compression = "lzw", background = "white")
print(p_final)
dev.off()

agg_png(paste0(stem, "_preview.png"), width = width_mm / 25.4, height = height_mm / 25.4,
        units = "in", res = 300, background = "white")
print(p_final)
dev.off()

fwrite(plot_dt, file.path(output_dir, "source_data", "Figure2G_compact_dotplot_source_data.csv"))
fwrite(side_long, file.path(output_dir, "source_data", "Figure2G_compact_side_table_source_data.csv"))
writeLines(caption_text, file.path(output_dir, "Figure2G_compact_figure_legend.txt"))
writeLines(c(
  "Figure 2G compact-label revision QA notes",
  "Reuse level=style-only inheritance from the frozen Figure 2G result.",
  "Data changes=none; all 17 rows, 25 genes and side-table values retained.",
  "Healthy row labels=ASPC, Papillary FB, ECM-FB, Pericyte, VSMC.",
  "Disease row labels=Q1-Q12, with Q8 and Q11 retained as ASPC-like priorities.",
  "Gene labels=Arial 10 pt, vertical 90 degrees.",
  "Side headers=Cells, Donors, Conf., Unc.%, Dev.%.",
  "Legends=two horizontal guides placed side by side below the dot matrix.",
  "Full identities and abbreviation definitions=accompanying Figure legend text file.",
  "Final size=183 x 105 mm; margins minimized; editable PDF/SVG and 600-dpi TIFF."
), file.path(output_dir, "Figure2G_compact_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2G_compact_sessionInfo.txt"))

message("Done: ", normalizePath(output_dir))
