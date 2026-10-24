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

project_dir <- normalizePath(".", mustWork = TRUE)
source_file <- file.path(
  project_dir,
  "outputs/ADSC_IL17A_0_vs_50_KEGG/source_data/ADSC_IL17A_0_vs_50_KEGG_all_tested_pathways.csv"
)
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_0_vs_50_KEGG_bubble_uncompressed_74x125mm")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(source_file)) stop("KEGG source table is missing.")

width_mm <- 74
# The 104 mm draft forced 20 pathways, two titles, axes and two scales into an
# overly compressed vertical space. The selected uncompressed version keeps
# the 74 mm column width and extends the panel height to 125 mm.
height_mm <- 125

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

# Figure contract:
# - Claim: 50 ng/ml IL-17A produces distinct up- and downregulated KEGG programs.
# - Evidence: the same 10 most significant pathways per direction are shown in
#   two visual encodings, without changing enrichment statistics.
# - Archetype: compact quantitative comparison.
# - Export: editable PDF/SVG plus 600 dpi TIFF/PNG at exactly 74 x 104 mm.

kegg <- fread(source_file, check.names = FALSE)
required <- c(
  "Description", "p.adjust", "Count", "Direction", "RichFactor",
  "GeneRatio", "gene_symbols"
)
if (!all(required %in% names(kegg))) {
  stop("The KEGG source table is missing required columns.")
}

kegg <- kegg[
  is.finite(p.adjust) & p.adjust > 0 &
    is.finite(Count) & Count > 0 &
    is.finite(RichFactor) & RichFactor > 0 &
    Direction %chin% c("Up in 50 ng/ml", "Down in 50 ng/ml")
]
kegg[, minus_log10_adjusted_p := -log10(p.adjust)]

select_top10 <- function(dt, direction) {
  out <- copy(dt[Direction == direction & p.adjust < 0.05])
  setorder(out, p.adjust, -Count, Description)
  if (nrow(out) < 10L) stop("Fewer than 10 significant pathways for ", direction)
  out[seq_len(10L)]
}

selected <- rbindlist(list(
  select_top10(kegg, "Up in 50 ng/ml"),
  select_top10(kegg, "Down in 50 ng/ml")
), use.names = TRUE)

compact_label <- function(x) {
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
  idx <- match(out, names(replacements))
  out[!is.na(idx)] <- unname(replacements[idx[!is.na(idx)]])
  out
}

selected[, plot_label := compact_label(Description)]
selected[, direction_short := fifelse(
  Direction == "Up in 50 ng/ml", "Up in 50 ng/ml", "Down in 50 ng/ml"
)]

fwrite(
  selected,
  file.path(source_dir, "ADSC_IL17A_0_vs_50_KEGG_top10_each_direction_used.csv")
)

font_family <- "Arial"
font_size <- 10
text_mm <- font_size / ggplot2::.pt

# -----------------------------------------------------------------------------
# Version 1: reference-image style count circles plus horizontal bars.
# -----------------------------------------------------------------------------
bar_xmax <- ceiling(max(selected$minus_log10_adjusted_p) * 1.03)

prepare_bar <- function(direction) {
  x <- copy(selected[Direction == direction])
  x[, row_id := rev(seq_len(.N))]
  x
}

bar_panel <- function(dt, strip_fill, bar_fill, show_title, show_x) {
  strip_text <- unique(dt$direction_short)

  ggplot(dt) +
    annotate(
      "rect", xmin = -1.46, xmax = -1.00, ymin = 0.5, ymax = 10.5,
      fill = strip_fill, color = NA
    ) +
    annotate(
      "text", x = -1.23, y = 5.5, label = strip_text, angle = 90,
      family = font_family, fontface = "bold", size = text_mm,
      color = "white"
    ) +
    geom_col(
      aes(x = minus_log10_adjusted_p, y = row_id),
      width = 0.73, fill = bar_fill, color = NA
    ) +
    geom_point(
      aes(x = -0.55, y = row_id),
      shape = 21, size = 4.8, stroke = 0.42,
      fill = bar_fill, color = "#2F2F2F"
    ) +
    geom_text(
      aes(x = -0.55, y = row_id, label = Count),
      family = font_family, fontface = "bold", size = 8.7 / ggplot2::.pt,
      color = "#202020"
    ) +
    geom_text(
      aes(x = 0.08, y = row_id, label = plot_label),
      hjust = 0, family = font_family, fontface = "bold",
      size = text_mm, color = "#171717", lineheight = 0.88
    ) +
    scale_x_continuous(
      limits = c(-1.50, bar_xmax),
      breaks = seq(0, bar_xmax, by = 2),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_y_continuous(limits = c(0.5, 10.5), expand = c(0, 0)) +
    labs(
      title = if (show_title) "KEGG pathways" else NULL,
      x = if (show_x) "-log10(adjusted P)" else NULL,
      y = NULL
    ) +
    coord_cartesian(clip = "on") +
    theme_classic(base_family = font_family, base_size = font_size) +
    theme(
      text = element_text(
        family = font_family, face = "bold", size = font_size, color = "black"
      ),
      plot.title = element_text(
        family = font_family, face = "bold", size = font_size,
        hjust = 0.5, margin = margin(b = 1.0, unit = "mm")
      ),
      axis.title.x = element_text(
        family = font_family, face = "bold", size = font_size,
        margin = margin(t = 1.0, unit = "mm")
      ),
      axis.text.x = if (show_x) {
        element_text(family = font_family, face = "bold", size = font_size)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.35, color = "black")
      } else {
        element_blank()
      },
      axis.line.x = if (show_x) {
        element_line(linewidth = 0.35, color = "black")
      } else {
        element_blank()
      },
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(
        t = if (show_title) 0.7 else 0.2,
        r = 0.5, b = if (show_x) 0.8 else 0.2, l = 0.5,
        unit = "mm"
      )
    )
}

bar_plot <- bar_panel(
  prepare_bar("Up in 50 ng/ml"),
  strip_fill = "#C95E68", bar_fill = "#E8A4AA",
  show_title = TRUE, show_x = FALSE
) / bar_panel(
  prepare_bar("Down in 50 ng/ml"),
  strip_fill = "#4A7FA8", bar_fill = "#A9C8DE",
  show_title = FALSE, show_x = TRUE
) + plot_layout(heights = c(1, 1))

bar_prefix <- file.path(
  output_dir,
  "ADSC_IL17A_0_vs_50_KEGG_top10_up_down_count_bar_74x104mm"
)
# The user selected the bubble version for final delivery. Keep the bar code as
# an archived alternative, but do not render it in the final workflow.
export_bar_version <- FALSE
if (export_bar_version) save_figure(bar_plot, bar_prefix)

# -----------------------------------------------------------------------------
# Version 2: reference-image style enrichment bubble plot.
# -----------------------------------------------------------------------------
prepare_bubble <- function(direction) {
  x <- copy(selected[Direction == direction])
  x[, plot_label := factor(plot_label, levels = rev(plot_label))]
  x
}

common_rich_range <- range(selected$RichFactor, finite = TRUE)
rich_pad <- max(diff(common_rich_range) * 0.18, 0.040)
common_sig_range <- range(selected$minus_log10_adjusted_p, finite = TRUE)
# Reference palette requested by the user: lower -log10(Padj) values are
# bright pink and progressively stronger enrichment shifts through magenta
# and violet to deep blue. Only the colour mapping changes; data are untouched.
reference_sig_colors <- c("#F21865", "#C51BC4", "#7336D4", "#244FD6")

bubble_panel <- function(dt, direction_title, show_title, show_x) {
  ggplot(
    dt,
    aes(
      x = RichFactor, y = plot_label,
      size = Count, fill = minus_log10_adjusted_p
    )
  ) +
    geom_point(shape = 21, color = "#454545", stroke = 0.28) +
    scale_x_continuous(
      limits = c(common_rich_range[1] - rich_pad, common_rich_range[2] + rich_pad),
      breaks = scales::pretty_breaks(n = 3),
      expand = expansion(mult = c(0, 0))
    ) +
    scale_size_continuous(
      range = c(1.7, 4.8), breaks = c(15, 30, 45),
      limits = c(min(selected$Count), max(selected$Count)),
      name = "Gene count",
      guide = guide_legend(
        title.position = "top", nrow = 1, byrow = TRUE, order = 1,
        override.aes = list(fill = "white", color = "black")
      )
    ) +
    scale_fill_gradientn(
      colors = reference_sig_colors,
      limits = common_sig_range,
      breaks = scales::pretty_breaks(n = 3)(common_sig_range),
      name = "-log10 Padj",
      guide = guide_colorbar(
        direction = "horizontal", title.position = "top", title.hjust = 0,
        barwidth = unit(20, "mm"), barheight = unit(2.4, "mm"),
        frame.colour = "#555555", ticks.colour = "#555555", order = 2
      )
    ) +
    labs(
      title = NULL,
      subtitle = direction_title,
      x = if (show_x) "Rich factor" else NULL,
      y = NULL
    ) +
    theme_classic(base_family = font_family, base_size = font_size) +
    theme(
      text = element_text(
        family = font_family, face = "bold", size = font_size, color = "black"
      ),
      plot.subtitle = element_text(
        family = font_family, face = "bold", size = font_size,
        hjust = 0.5, margin = margin(b = 0.7, unit = "mm")
      ),
      plot.title.position = "plot",
      axis.title.x = element_text(
        family = font_family, face = "bold", size = font_size,
        margin = margin(t = 0.7, unit = "mm")
      ),
      axis.text.x = if (show_x) {
        element_text(family = font_family, face = "bold", size = font_size)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.32, color = "black")
      } else {
        element_blank()
      },
      axis.line.x = if (show_x) {
        element_line(linewidth = 0.32, color = "black")
      } else {
        element_blank()
      },
      axis.text.y = element_text(
        family = font_family, face = "bold", size = font_size,
        color = "black", margin = margin(r = 1.2, unit = "mm")
      ),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.border = element_rect(
        fill = NA, color = "#333333", linewidth = 0.42
      ),
      panel.grid.major.x = element_line(color = "#D6D6D6", linewidth = 0.25),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.margin = margin(
        t = 0.4,
        r = 0.6, b = if (show_x) 0.7 else 0.5, l = 0.5,
        unit = "mm"
      )
    )
}

# Use a single 20-row panel rather than two short facets. This preserves the
# requested 10 pt type while giving every pathway a full, non-overlapping row.
bubble_dt <- copy(selected)
bubble_dt[, row_id := c(20:11, 10:1)]
bubble_y_labels <- setNames(bubble_dt$plot_label, bubble_dt$row_id)

bubble_body <- ggplot(
  bubble_dt,
  aes(
    x = RichFactor, y = row_id,
    size = Count, fill = minus_log10_adjusted_p
  )
) +
  annotate(
    "rect", xmin = 0.407, xmax = 0.438, ymin = 10.5, ymax = 20.5,
    fill = "#C95E68", color = NA
  ) +
  annotate(
    "rect", xmin = 0.407, xmax = 0.438, ymin = 0.5, ymax = 10.5,
    fill = "#4A7FA8", color = NA
  ) +
  annotate(
    "text", x = 0.4225, y = 15.5, label = "UP", angle = 90,
    family = font_family, fontface = "bold", size = text_mm, color = "white"
  ) +
  annotate(
    "text", x = 0.4225, y = 5.5, label = "DOWN", angle = 90,
    family = font_family, fontface = "bold", size = text_mm, color = "white"
  ) +
  geom_hline(yintercept = 10.5, color = "white", linewidth = 0.7) +
  geom_point(shape = 21, color = "#454545", stroke = 0.28) +
  scale_x_continuous(
    # Extra left padding prevents large low-RichFactor bubbles from being
    # clipped by the full panel border; data coordinates are unchanged.
    limits = c(0.10, 0.44), breaks = c(0.15, 0.25, 0.35),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(
    limits = c(0.5, 20.5), breaks = 1:20,
    labels = unname(bubble_y_labels[as.character(1:20)]), expand = c(0, 0)
  ) +
  scale_size_continuous(
    range = c(1.7, 4.0), limits = c(min(selected$Count), max(selected$Count)),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colors = reference_sig_colors,
    limits = common_sig_range, guide = "none"
  ) +
  labs(title = "KEGG enrichment", x = "Rich factor", y = NULL) +
  theme_classic(base_family = font_family, base_size = font_size) +
  theme(
    text = element_text(
      family = font_family, face = "bold", size = font_size, color = "black"
    ),
    plot.title = element_text(
      family = font_family, face = "bold", size = font_size,
      hjust = 0.5, margin = margin(b = 0.7, unit = "mm")
    ),
    plot.title.position = "plot",
    axis.title.x = element_text(
      family = font_family, face = "bold", size = font_size,
      margin = margin(t = 0.7, unit = "mm")
    ),
    axis.text.x = element_text(
      family = font_family, face = "bold", size = font_size, color = "black"
    ),
    axis.text.y = element_text(
      family = font_family, face = "bold", size = font_size,
      color = "black", margin = margin(r = 1.0, unit = "mm")
    ),
    axis.ticks.x = element_line(linewidth = 0.32, color = "black"),
    axis.ticks.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    panel.border = element_rect(
      fill = NA, color = "#333333", linewidth = 0.42
    ),
    panel.grid.major.x = element_line(color = "#D6D6D6", linewidth = 0.25),
    panel.grid.minor = element_blank(),
    plot.margin = margin(0.5, 0.5, 0.3, 0.5, unit = "mm")
  )

# Final selected layout: UP and DOWN are separate, independently framed
# panels, vertically arranged with identical Rich factor and color/size scales.
bubble_body <- bubble_panel(
  prepare_bubble("Up in 50 ng/ml"),
  direction_title = "Up in 50 ng/ml", show_title = FALSE, show_x = FALSE
) / bubble_panel(
  prepare_bubble("Down in 50 ng/ml"),
  direction_title = "Down in 50 ng/ml", show_title = FALSE, show_x = TRUE
) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "KEGG enrichment",
    theme = theme(
      plot.title = element_text(
        family = font_family, face = "bold", size = font_size,
        hjust = 0.5, margin = margin(b = 0.4, unit = "mm")
      )
    )
  )

# Compact custom legend: keeps 10 pt text while limiting the legend to 16% of
# the 104 mm panel height. Rectangular gradient tiles remain vector objects.
count_legend <- data.table(
  x = c(7, 25, 43), y = 3.8,
  Count = c(15, 30, 45), label = c("15", "30", "45")
)
gradient_legend <- data.table(
  x = seq(66, 96, length.out = 80), y = 3.8,
  significance = seq(common_sig_range[1], common_sig_range[2], length.out = 80)
)
gradient_breaks <- scales::pretty_breaks(n = 3)(common_sig_range)
gradient_breaks <- gradient_breaks[
  gradient_breaks >= common_sig_range[1] & gradient_breaks <= common_sig_range[2]
]
gradient_break_x <- 66 +
  (gradient_breaks - common_sig_range[1]) /
  diff(common_sig_range) * 30

bubble_legend <- ggplot() +
  geom_point(
    data = count_legend,
    aes(x = x, y = y, size = Count),
    shape = 21, fill = "white", color = "black", stroke = 0.38
  ) +
  geom_text(
    data = count_legend,
    aes(x = x + 5.3, y = y, label = label),
    hjust = 0, family = font_family, fontface = "bold", size = text_mm
  ) +
  geom_tile(
    data = gradient_legend,
    aes(x = x, y = y, fill = significance),
    width = 30 / 80 * 1.05, height = 1.55
  ) +
  annotate(
    "rect", xmin = 65.8, xmax = 96.2, ymin = 3.0, ymax = 4.6,
    fill = NA, color = "#555555", linewidth = 0.34
  ) +
  annotate(
    "text", x = 2, y = 8.0, label = "Gene count",
    hjust = 0, family = font_family, fontface = "bold", size = text_mm
  ) +
  annotate(
    "text", x = 65.8, y = 8.0, label = "-log10 Padj",
    hjust = 0, family = font_family, fontface = "bold", size = text_mm
  ) +
  geom_text(
    data = data.table(x = gradient_break_x, y = 1.25, label = format(gradient_breaks)),
    aes(x = x, y = y, label = label),
    family = font_family, fontface = "bold", size = text_mm
  ) +
  scale_size_continuous(
    range = c(2.0, 4.5), limits = c(min(selected$Count), max(selected$Count)),
    guide = "none"
  ) +
  scale_fill_gradientn(
    colors = reference_sig_colors,
    limits = common_sig_range, guide = "none"
  ) +
  coord_cartesian(xlim = c(0, 100), ylim = c(0, 10), clip = "off") +
  theme_void(base_family = font_family, base_size = font_size) +
  theme(plot.margin = margin(0, 0.4, 0, 0.4, unit = "mm"))

bubble_plot <- wrap_elements(full = bubble_body) /
  wrap_elements(full = bubble_legend) +
  plot_layout(heights = c(0.91, 0.09))

bubble_prefix <- file.path(
  output_dir,
  "ADSC_IL17A_0_vs_50_KEGG_top10_up_down_bubble_reference_palette_74x125mm"
)
save_figure(bubble_plot, bubble_prefix)

qa_lines <- c(
  "ADSC IL-17A KEGG two-style figure QA",
  paste0("Canvas: ", width_mm, " x ", height_mm, " mm"),
  "Font: Arial Bold 10 pt throughout both panels and legends.",
  "Selection: p.adjust < 0.05; top 10 per direction ordered by adjusted P, Count, and name.",
  "No TPM expression threshold was applied to KEGG input-gene selection.",
  "KEGG input follows the standard edgeR DEG rule: FDR < 0.05, split by the sign of logFC.",
  paste0("Up pathways shown: ", sum(selected$Direction == "Up in 50 ng/ml")),
  paste0("Down pathways shown: ", sum(selected$Direction == "Down in 50 ng/ml")),
  "Bar version: bar length = -log10(adjusted P); numbered circle = gene Count.",
  "Bubble version: x = RichFactor; bubble size = Count; color = -log10(adjusted P).",
  "Bubble color palette: vivid pink -> magenta -> violet -> deep blue, matching the supplied reference image.",
  "The two bubble plotting areas use a complete thin dark outer frame.",
  "Full pathway names and unmodified enrichment statistics are retained in source_data.",
  "No enrichment statistic was recomputed or altered for styling."
)
writeLines(qa_lines, file.path(output_dir, "ADSC_IL17A_KEGG_two_styles_QA_notes.txt"))

sink(file.path(output_dir, "ADSC_IL17A_KEGG_two_styles_sessionInfo.txt"))
sessionInfo()
sink()

message("KEGG two-style figures written to: ", output_dir)
