#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

options(stringsAsFactors = FALSE)

if (!"Arial" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("sans")[[1L]])
}
if (!"Arial" %in% names(grDevices::postscriptFonts())) {
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts("sans")[[1L]])
}
if (capabilities("aqua") && !"Arial" %in% names(grDevices::quartzFonts())) {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

project_dir <- normalizePath(".", mustWork = TRUE)
kegg_source <- file.path(
  project_dir,
  "outputs/ADSC_IL17A_0_vs_50_KEGG/source_data/ADSC_IL17A_0_vs_50_KEGG_all_tested_pathways.csv"
)
go_mf_source <- file.path(
  project_dir,
  "outputs/ADSC_IL17A_0_vs_50_GO_all_ontologies/source_data/ADSC_IL17A_0_vs_50_GO_MF_all_tested_terms.csv"
)
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_0_vs_50_KEGG_GO_MF_74x104mm")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

for (path in c(kegg_source, go_mf_source)) {
  if (!file.exists(path)) stop("Required source table is missing: ", basename(path))
}

width_mm <- 74
height_mm <- 104

save_figure <- function(plot, prefix) {
  if (capabilities("aqua")) {
    quartz(
      type = "pdf", file = paste0(prefix, ".pdf"),
      width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial"
    )
  } else {
    cairo_pdf(
      paste0(prefix, ".pdf"),
      width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial"
    )
  }
  print(plot)
  dev.off()

  svglite(
    paste0(prefix, ".svg"),
    width = width_mm / 25.4, height = height_mm / 25.4,
    system_fonts = list(sans = "Arial")
  )
  print(plot)
  dev.off()

  agg_tiff(
    paste0(prefix, ".tiff"),
    width = width_mm / 25.4, height = height_mm / 25.4,
    units = "in", res = 600, compression = "lzw",
    background = "white"
  )
  print(plot)
  dev.off()

  agg_png(
    paste0(prefix, ".png"),
    width = width_mm / 25.4, height = height_mm / 25.4,
    units = "in", res = 600, background = "white"
  )
  print(plot)
  dev.off()
}

wrap_label <- function(x, width) {
  vapply(x, function(value) paste(strwrap(value, width = width), collapse = "\n"), character(1L))
}

# -----------------------------------------------------------------------------
# KEGG: top 10 BH-adjusted significant pathways in each direction.
# -----------------------------------------------------------------------------
kegg_all <- fread(kegg_source, check.names = FALSE)
required_kegg <- c("Description", "p.adjust", "Count", "Direction")
if (!all(required_kegg %in% names(kegg_all))) {
  stop("The KEGG source table is missing required columns.")
}
kegg_all <- kegg_all[
  is.finite(p.adjust) & p.adjust > 0 &
    Direction %chin% c("Up in 50 ng/ml", "Down in 50 ng/ml")
]
kegg_all[, minus_log10_adjusted_p := -log10(p.adjust)]

select_top10 <- function(x, direction) {
  out <- copy(x[Direction == direction & p.adjust < 0.05])
  setorder(out, p.adjust, -Count, Description)
  if (nrow(out) < 10L) stop("Fewer than 10 significant KEGG pathways for ", direction)
  out[seq_len(10L)]
}

kegg_up <- select_top10(kegg_all, "Up in 50 ng/ml")
kegg_down <- select_top10(kegg_all, "Down in 50 ng/ml")
kegg_xmax <- max(c(kegg_up$minus_log10_adjusted_p, kegg_down$minus_log10_adjusted_p)) * 1.03

compact_kegg_label <- function(x) {
  replacements <- c(
    "IL-17 signaling pathway" = "IL-17 signaling",
    "Viral protein interaction with cytokine and cytokine receptor" = "Viral protein-cytokine axis",
    "Cytoskeleton in muscle cells" = "Cytoskeleton-muscle",
    "Pathways in cancer" = "Cancer pathways",
    "Human papillomavirus infection" = "HPV infection",
    "Proteoglycans in cancer" = "Cancer proteoglycans",
    "PI3K-Akt signaling pathway" = "PI3K-Akt signaling",
    "Arrhythmogenic right ventricular cardiomyopathy" = "Arrhythmogenic RV CM"
  )
  out <- x
  matched <- match(out, names(replacements))
  out[!is.na(matched)] <- unname(replacements[matched[!is.na(matched)]])
  out
}

prepare_kegg <- function(x) {
  x <- copy(x)
  x[, pathway_label := compact_kegg_label(Description)]
  x[, pathway_label := factor(pathway_label, levels = rev(pathway_label))]
  score_range <- range(x$minus_log10_adjusted_p, finite = TRUE)
  switch_point <- score_range[1L] + 0.44 * diff(score_range)
  x[, label_color := ifelse(minus_log10_adjusted_p >= switch_point, "white", "black")]
  x
}

kegg_up_plot <- prepare_kegg(kegg_up)
kegg_down_plot <- prepare_kegg(kegg_down)

make_kegg_panel <- function(plot_dt, direction_label, colors, show_x, show_title) {
  score_range <- range(plot_dt$minus_log10_adjusted_p, finite = TRUE)
  if (diff(score_range) == 0) score_range <- score_range + c(-0.1, 0.1)

  ggplot(plot_dt, aes(x = minus_log10_adjusted_p, y = pathway_label)) +
    geom_col(
      aes(fill = minus_log10_adjusted_p),
      width = 0.78, color = NA
    ) +
    geom_text(
      aes(x = 0.055, label = pathway_label, color = label_color),
      hjust = 0, vjust = 0.5,
      family = "Arial", fontface = "bold",
      size = 10 / ggplot2::.pt, lineheight = 0.82,
      show.legend = FALSE
    ) +
    scale_fill_gradientn(
      colors = colors,
      limits = score_range,
      breaks = scales::pretty_breaks(n = 3)(score_range),
      name = "-log10 Padj",
      guide = guide_colorbar(
        title.position = "top", title.hjust = 0,
        direction = "horizontal",
        barwidth = unit(18, "mm"), barheight = unit(2.4, "mm"),
        ticks.colour = "#333333", frame.colour = "#555555"
      )
    ) +
    scale_color_identity() +
    scale_x_continuous(
      limits = c(0, kegg_xmax),
      expand = expansion(mult = c(0, 0)),
      breaks = scales::pretty_breaks(n = 4)
    ) +
    labs(
      title = if (show_title) "KEGG pathways" else NULL,
      # The horizontal colour scale below each panel already names this metric.
      # Keep the lower-panel tick labels, but avoid a duplicated axis title.
      x = NULL,
      y = direction_label
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
      text = element_text(family = "Arial", face = "bold", color = "black", size = 10),
      plot.title = element_text(
        family = "Arial", size = 10, face = "bold", hjust = 0.5,
        margin = margin(b = 0.8, unit = "mm")
      ),
      axis.title = element_text(family = "Arial", size = 10, face = "bold", color = "black"),
      axis.title.x = element_text(margin = margin(t = 0.7, unit = "mm")),
      axis.title.y = element_text(margin = margin(r = 0.8, unit = "mm")),
      axis.text.x = if (show_x) {
        element_text(family = "Arial", size = 10, face = "bold", color = "black")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.30, color = "black")
      } else {
        element_blank()
      },
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "#3F3F3F", linewidth = 0.38),
      panel.grid.major.x = element_line(
        color = "#CFCFCF", linewidth = 0.25, linetype = "dashed"
      ),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.box = "horizontal",
      legend.box.just = "center",
      legend.title = element_text(
        family = "Arial", size = 10, face = "bold", color = "black",
        lineheight = 0.90
      ),
      legend.text = element_text(family = "Arial", size = 10, face = "bold", color = "black"),
      legend.margin = margin(0.1, 0, 0, 0, unit = "mm"),
      legend.box.spacing = unit(0.2, "mm"),
      plot.margin = margin(
        t = if (show_title) 0.8 else 0.4,
        r = 0.8, b = if (show_x) 0.5 else 0.3, l = 0.8,
        unit = "mm"
      )
    )
}

kegg_red <- c("#FDE8E5", "#F8AAA1", "#ED6A5A", "#C92D3A", "#760016")
kegg_blue <- c("#E6F1FA", "#A9CDE8", "#639FCF", "#2768AE", "#08306B")

kegg_plot <- make_kegg_panel(
  kegg_up_plot, "Up in 50 ng/ml", kegg_red, FALSE, TRUE
) / make_kegg_panel(
  kegg_down_plot, "Down in 50 ng/ml", kegg_blue, TRUE, FALSE
) + plot_layout(heights = c(1, 1))

kegg_prefix <- file.path(
  output_dir, "ADSC_IL17A_0_vs_50_KEGG_top10_up_down_74x104mm"
)
save_figure(kegg_plot, kegg_prefix)

# -----------------------------------------------------------------------------
# GO Molecular Function: top 10 BH-adjusted significant terms per direction.
# -----------------------------------------------------------------------------
go_all <- fread(go_mf_source, check.names = FALSE)
required_go <- c("Description", "p.adjust", "GeneRatio", "overlap", "Direction")
if (!all(required_go %in% names(go_all))) {
  stop("The GO Molecular Function source table is missing required columns.")
}
go_all <- go_all[
  is.finite(p.adjust) & p.adjust > 0 & is.finite(GeneRatio) & GeneRatio > 0 &
    Direction %chin% c("Up in 50 ng/ml", "Down in 50 ng/ml")
]

select_go_top10 <- function(x, direction) {
  out <- copy(x[Direction == direction & p.adjust < 0.05])
  setorder(out, p.adjust, -overlap, Description)
  if (nrow(out) < 10L) stop("Fewer than 10 significant GO-MF terms for ", direction)
  out[seq_len(10L)]
}

go_up <- select_go_top10(go_all, "Up in 50 ng/ml")
go_down <- select_go_top10(go_all, "Down in 50 ng/ml")

compact_go_label <- function(x) {
  x <- sub("^positive regulation of ", "pos. reg. of ", x)
  x <- sub("^negative regulation of ", "neg. reg. of ", x)
  x <- sub("^regulation of ", "reg. of ", x)
  x <- gsub("extracellular matrix", "ECM", x, fixed = TRUE)
  x <- gsub("structural constituent of ribosome", "ribosomal structural constituent", x, fixed = TRUE)
  x <- gsub("ATP-dependent protein folding chaperone", "ATP-dependent folding chaperone", x, fixed = TRUE)
  x <- gsub("signaling receptor activator activity", "receptor activator activity", x, fixed = TRUE)
  x <- gsub("signaling receptor regulator activity", "receptor regulator activity", x, fixed = TRUE)
  x
}

prepare_go <- function(x) {
  x <- copy(x)
  x[, process_label := compact_go_label(Description)]
  x[, plot_y := rev(seq_len(.N))]
  x
}

go_up_plot <- prepare_go(go_up)
go_down_plot <- prepare_go(go_down)
go_ratio_range <- range(c(go_up$GeneRatio, go_down$GeneRatio), finite = TRUE)
go_p_range <- range(c(go_up$p.adjust, go_down$p.adjust), finite = TRUE)
go_p_breaks <- 10^seq(log10(go_p_range[1L]), log10(go_p_range[2L]), length.out = 3L)
go_colors <- c("#173CCB", "#6A28C2", "#C3159A", "#FF1744")

make_go_panel <- function(plot_dt, title_text) {
  ggplot(plot_dt, aes(y = plot_y)) +
    annotate(
      "rect", xmin = 0.965, xmax = 1.035,
      ymin = 0.5, ymax = nrow(plot_dt) + 0.5,
      fill = NA, color = "#555555", linewidth = 0.36
    ) +
    geom_text(
      aes(x = 0.92, label = process_label),
      hjust = 1, family = "Arial", fontface = "bold",
      size = 9 / ggplot2::.pt, color = "black", lineheight = 0.88
    ) +
    geom_point(
      aes(x = 1, size = GeneRatio, color = p.adjust),
      alpha = 1, stroke = 0
    ) +
    scale_size_continuous(
      name = "GeneRatio", range = c(1.8, 4.0), limits = go_ratio_range,
      breaks = c(0.02, 0.05, 0.08),
      guide = guide_legend(
        order = 1, title.position = "top", title.hjust = 0,
        nrow = 1, override.aes = list(color = "black", alpha = 1)
      )
    ) +
    scale_color_gradientn(
      name = "Adjusted P",
      colors = go_colors, limits = go_p_range, trans = "log10",
      breaks = go_p_breaks,
      labels = scales::label_scientific(digits = 1),
      guide = guide_colorbar(
        order = 2, title.position = "top", title.hjust = 0,
        direction = "horizontal",
        barwidth = unit(19, "mm"), barheight = unit(2.6, "mm"),
        frame.colour = "#555555", ticks.colour = "#555555"
      )
    ) +
    scale_x_continuous(limits = c(-1.26, 1.12), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, nrow(plot_dt) + 0.5), expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_classic(base_family = "Arial", base_size = 9) +
    theme(
      text = element_text(family = "Arial", face = "bold", color = "black", size = 9),
      plot.title = element_text(
        family = "Arial", size = 10, face = "bold", hjust = 0.5,
        margin = margin(b = 0.5, unit = "mm")
      ),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_blank(),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.box = "vertical",
      legend.box.just = "center",
      legend.title = element_text(family = "Arial", size = 9, face = "bold", color = "black"),
      legend.text = element_text(family = "Arial", size = 9, face = "bold", color = "black"),
      legend.key.height = unit(2.8, "mm"),
      legend.spacing.x = unit(0.7, "mm"),
      legend.spacing.y = unit(0.4, "mm"),
      legend.margin = margin(0.2, 0, 0, 0, unit = "mm"),
      plot.margin = margin(0.3, 0.8, 0.2, 0.8, unit = "mm")
    )
}

go_plot <- make_go_panel(
  go_up_plot, "Upregulated in 50 ng/ml"
) / make_go_panel(
  go_down_plot, "Downregulated in 50 ng/ml"
) +
  plot_layout(heights = c(1, 1), guides = "collect") +
  plot_annotation(
    title = "GO Molecular Function",
    theme = theme(
      plot.title = element_text(
        family = "Arial", size = 10, face = "bold", hjust = 0.5,
        color = "black", margin = margin(b = 0.3, unit = "mm")
      )
    )
  ) & theme(legend.position = "bottom", legend.box = "vertical")

go_prefix <- file.path(
  output_dir, "ADSC_IL17A_0_vs_50_GO_MF_top10_up_down_74x104mm"
)
save_figure(go_plot, go_prefix)

fwrite(
  rbindlist(list(kegg_up, kegg_down), use.names = TRUE, fill = TRUE),
  file.path(source_dir, "ADSC_IL17A_0_vs_50_KEGG_top10_each_direction_74x104mm.csv")
)
fwrite(
  rbindlist(list(go_up, go_down), use.names = TRUE, fill = TRUE),
  file.path(source_dir, "ADSC_IL17A_0_vs_50_GO_MF_top10_each_direction_74x104mm.csv")
)

qa_notes <- c(
  "Core conclusion: the panels summarize the most significant up- and downregulated pathway enrichments after 50 ng/ml IL-17A treatment.",
  "Archetype: compact quantitative grid designed for a fixed A4 manuscript panel slot.",
  "Backend: R only (ggplot2 and patchwork).",
  "Source reuse level: structural adaptation of the existing KEGG and GO-MF enrichment outputs; enrichment statistics were not recomputed or altered.",
  "KEGG selection: BH-adjusted P < 0.05; top 10 per direction ordered by adjusted P, then overlap count and pathway name.",
  sprintf("KEGG availability: %d significant upregulated and %d significant downregulated pathways before top-10 display selection.",
          kegg_all[Direction == "Up in 50 ng/ml" & p.adjust < 0.05, .N],
          kegg_all[Direction == "Down in 50 ng/ml" & p.adjust < 0.05, .N]),
  "GO-MF selection: BH-adjusted P < 0.05; top 10 per direction ordered by adjusted P, then overlap count and term name.",
  sprintf("GO-MF availability: %d significant upregulated and %d significant downregulated terms before top-10 display selection.",
          go_all[Direction == "Up in 50 ng/ml" & p.adjust < 0.05, .N],
          go_all[Direction == "Down in 50 ng/ml" & p.adjust < 0.05, .N]),
  "Both figures use a final physical size of exactly 74 mm wide by 104 mm high.",
  "KEGG typography: Arial Bold 10 pt for pathway labels, axes, titles, and color scales.",
  "GO-MF typography: Arial Bold 9 pt for term labels and legends; 10 pt for panel and figure titles.",
  "KEGG panels are stacked vertically with red and blue adjusted-P gradients for up- and downregulated pathways.",
  "GO-MF panels are stacked vertically with GeneRatio encoded by dot size and BH-adjusted P encoded by color.",
  "Exports: editable-text PDF and SVG plus 600 dpi TIFF and PNG."
)
writeLines(
  qa_notes,
  file.path(output_dir, "ADSC_IL17A_0_vs_50_KEGG_GO_MF_74x104mm_QA_notes.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "ADSC_IL17A_0_vs_50_KEGG_GO_MF_74x104mm_sessionInfo.txt")
)

message("Reformatted KEGG and GO-MF figures written to: ", output_dir)
