suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
input_csv <- if (length(args) >= 1L) args[[1]] else
  "outputs/Figure2H_disease_up_two_layouts/source_data/Figure2H_disease_up_selected_pathways.csv"
output_dir <- if (length(args) >= 2L) args[[2]] else
  "outputs/Figure2H_compact_horizontal_bar"

input_csv <- normalizePath(input_csv, mustWork = TRUE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)

# Figure contract
# Core conclusion: compactly compare the strongest disease-up Hallmark programs
# in Q8 and Q11 and the best-supported Q8-vs-Q11 directional differences.
# Evidence: Q8 and Q11 show the eight disease-high pathways selected previously;
# the direct panel shows four Q8-high and four Q11-high pathways.
# Archetype: compact quantitative small-multiple grid.
# Integrity: NES and q values are preserved exactly; no donor dots or error bars
# are added because the GSEA table contains one enrichment estimate per pathway.
# Export: 183 x 52 mm, editable PDF/SVG, 600-dpi TIFF, and PNG preview.

d <- fread(input_csv)
required <- c(
  "comparison_display", "direction_display", "pathway_label", "NES",
  "q_value", "selection_rank", "selection_rule"
)
missing_cols <- setdiff(required, names(d))
if (length(missing_cols)) {
  stop("Input table is missing: ", paste(missing_cols, collapse = ", "))
}

expected_counts <- c(
  "Q8 vs Healthy" = 8L,
  "Q11 vs Healthy" = 8L,
  "Q8 vs Q11" = 8L
)
observed_counts <- table(d$comparison_display)
if (!all(names(expected_counts) %in% names(observed_counts)) ||
    any(observed_counts[names(expected_counts)] != expected_counts)) {
  stop("Expected exactly eight pathways in each comparison panel")
}

palette <- data.table(
  comparison_display = c("Q8 vs Healthy", "Q11 vs Healthy", "Q8 vs Q11"),
  strip_fill = c("#EBC6DA", "#FFF2A8", "#C9E8EF"),
  bar_fill = c("#E4B9D0", "#F8E996", "#B7DEE7"),
  bar_line = c("#555555", "#555555", "#555555"),
  strip_text = c("#9E456F", "#8B7410", "#277886")
)

format_q <- function(q) {
  ifelse(
    q < 0.001, "q<0.001",
    ifelse(q < 0.01, "q<0.01", ifelse(q < 0.05, "q<0.05", ""))
  )
}

wrap_label <- function(x, width = 17L) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

d[, `:=`(
  q_label = format_q(q_value),
  fdr_significant = q_value < 0.05
)]

# Standard pathway abbreviations are used only for display so that the 52-mm
# high layout remains legible at 8 pt. Full Hallmark names remain in source data.
compact_labels <- c(
  "PI3K AKT MTOR SIGNALING" = "PI3K-AKT-MTOR",
  "UNFOLDED PROTEIN RESPONSE" = "UPR",
  "EPITHELIAL MESENCHYMAL TRANSITION" = "EMT",
  "OXIDATIVE PHOSPHORYLATION" = "OXPHOS",
  "INTERFERON ALPHA RESPONSE" = "IFN-ALPHA RESPONSE",
  "TNFA SIGNALING VIA NFKB" = "TNFA/NFKB SIGNALING",
  "CHOLESTEROL HOMEOSTASIS" = "CHOLESTEROL\nHOMEOSTASIS"
)
d[, pathway_display := pathway_label]
d[pathway_label %in% names(compact_labels),
  pathway_display := unname(compact_labels[pathway_label])]

panel_titles <- c(
  "Q8 vs Healthy" = "Q8 disease-up",
  "Q11 vs Healthy" = "Q11 disease-up",
  "Q8 vs Q11" = "Q8-Q11"
)

font_family <- "Arial"
font_pt <- 8
text_size <- font_pt / ggplot2::.pt

prepare_panel <- function(comp) {
  z <- copy(d[comparison_display == comp])
  if (comp == "Q8 vs Q11") {
    order_labels <- c(
      z[direction_display == "Q8-high"][order(selection_rank), pathway_display],
      z[direction_display == "Q11-high"][order(selection_rank), pathway_display]
    )
    z[, pathway_plot := pathway_display]
  } else {
    order_labels <- z[order(selection_rank), pathway_display]
    z[, pathway_plot := pathway_display]
  }
  z[, pathway_plot := factor(pathway_plot, levels = rev(order_labels))]
  z[]
}

make_panel <- function(comp) {
  z <- prepare_panel(comp)
  pal <- palette[comparison_display == comp]
  z[, strip_label := unname(panel_titles[[comp]])]

  if (comp == "Q8 vs Healthy") {
    x_limits <- c(0, 1.72)
    x_breaks <- c(0, 0.75, 1.50)
  } else if (comp == "Q11 vs Healthy") {
    x_limits <- c(0, 2.10)
    x_breaks <- c(0, 1.0, 2.0)
  } else {
    x_limits <- c(-2.15, 2.48)
    x_breaks <- c(-2, 0, 2)
  }

  span <- diff(x_limits)
  z[, q_x := NES - sign(NES) * 0.045 * span]
  z[, q_hjust := fifelse(NES >= 0, 1, 0)]
  z[, bar_alpha := fifelse(fdr_significant, 0.92, 0.48)]

  p <- ggplot(z, aes(x = NES, y = pathway_plot)) +
    geom_vline(xintercept = 0, colour = "#202020", linewidth = 0.42) +
    geom_col(
      aes(alpha = bar_alpha), width = 0.52,
      fill = pal$bar_fill, colour = pal$bar_line,
      linewidth = 0.12
    ) +
    geom_text(
      aes(x = q_x, label = q_label, hjust = q_hjust),
      family = font_family, fontface = "plain", size = text_size,
      colour = "#111111", lineheight = 0.92
    ) +
    facet_wrap(~strip_label, ncol = 1) +
    scale_alpha_identity() +
    scale_x_continuous(
      breaks = x_breaks,
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = x_limits, clip = "off") +
    labs(x = NULL, y = NULL) +
    theme_classic(base_family = font_family, base_size = font_pt) +
    theme(
      text = element_text(
        family = font_family, face = "plain", size = font_pt,
        colour = "#111111"
      ),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.45, colour = "#202020"),
      axis.ticks.x = element_line(linewidth = 0.40, colour = "#202020"),
      axis.text.x = element_text(size = font_pt, face = "plain", margin = margin(t = 2)),
      axis.text.y = element_text(
        size = font_pt, face = "plain", lineheight = 0.88,
        margin = margin(r = 2.5)
      ),
      strip.background = element_rect(fill = pal$strip_fill, colour = NA),
      strip.text = element_text(
        size = font_pt, face = "plain", colour = pal$strip_text,
        lineheight = 0.88, margin = margin(t = 2.2, r = 2, b = 2.2, l = 2)
      ),
      panel.grid = element_blank(),
      panel.spacing = unit(0, "pt"),
      plot.margin = margin(2, 2, 2, 2)
    )

  if (comp == "Q8 vs Q11") {
    p <- p +
      geom_hline(yintercept = 4.5, linewidth = 0.34, linetype = 2, colour = "#777777")
  }
  p
}

p_q8 <- make_panel("Q8 vs Healthy")
p_q11 <- make_panel("Q11 vs Healthy")
p_direct <- make_panel("Q8 vs Q11")

figure <- (p_q8 | p_q11 | p_direct) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(
    caption = "Normalized enrichment score (NES)",
    theme = theme(
      plot.background = element_rect(fill = "white", colour = NA),
      plot.caption = element_text(
        family = font_family, face = "plain", size = font_pt,
        colour = "#111111", hjust = 0.5, margin = margin(t = 2)
      ),
      plot.margin = margin(2, 3, 2, 3)
    )
  )

base <- file.path(output_dir, "Figure2H_compact_horizontal_bar")
width_mm <- 183
height_mm <- 52
width_in <- width_mm / 25.4
height_in <- height_mm / 25.4

# cairo_pdf is the preferred cross-platform device. Set
# options(figure2h.use_cairo_pdf = TRUE) on systems with a working Cairo build.
# This macOS R build uses Quartz by default because its Cairo X11 libraries are
# unavailable; Quartz keeps all paths and text vector-editable.
if (isTRUE(getOption("figure2h.use_cairo_pdf", FALSE))) {
  grDevices::cairo_pdf(
    filename = paste0(base, ".pdf"), width = width_in, height = height_in,
    family = font_family, onefile = TRUE
  )
} else {
  quartz(
    type = "pdf", file = paste0(base, ".pdf"), width = width_in,
    height = height_in, family = font_family
  )
}
print(figure)
dev.off()

svglite::svglite(paste0(base, ".svg"), width = width_in, height = height_in,
                 system_fonts = list(Arial = "Arial"))
print(figure)
dev.off()

ragg::agg_tiff(
  paste0(base, ".tiff"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white"
)
print(figure)
dev.off()

ragg::agg_png(
  paste0(base, "_preview.png"), width = width_in, height = height_in,
  units = "in", res = 600, background = "white"
)
print(figure)
dev.off()

source_out <- copy(d)
for (col in names(source_out)) {
  if (is.list(source_out[[col]])) {
    source_out[, (col) := vapply(get(col), paste, collapse = ";", FUN.VALUE = character(1))]
  }
}
fwrite(
  source_out,
  file.path(output_dir, "source_data", "Figure2H_compact_horizontal_bar_source_data.csv")
)
fwrite(
  palette,
  file.path(output_dir, "source_data", "Figure2H_compact_horizontal_bar_palette.csv")
)

qa_lines <- c(
  "Figure 2H compact horizontal bar chart QA",
  "Canvas: 183 x 52 mm; height reduced by one third from the 78-mm draft.",
  "Typeface: Arial Regular, 8 pt for strip, axis, pathway, and q-value text.",
  "Q8 and Q11: disease-up pathways only; top 8 retained from the approved selection.",
  "Q8 vs Q11: four Q8-high and four Q11-high pathways retained from the approved selection.",
  "NES and q values were not recalculated or altered.",
  "Standard display abbreviations (EMT, OXPHOS, UPR, etc.) preserve legibility at 52-mm height; full names remain in source data.",
  "CHOLESTEROL HOMEOSTASIS is displayed on two lines to balance panel widths and inter-panel spacing.",
  "Only significant q values are displayed: q<0.001, q<0.01, or q<0.05; non-significant q labels are hidden.",
  "Exact unrounded q values remain available in the source-data CSV.",
  "Bar outlines use an ultra-fine 0.12-mm light-black (#555555) stroke.",
  "No donor dots or error bars were added because they are not present in the GSEA result table.",
  "Outputs: editable PDF/SVG plus 600-dpi TIFF and PNG preview."
)
writeLines(qa_lines, file.path(output_dir, "Figure2H_compact_horizontal_bar_QA.txt"))
capture.output(sessionInfo(), file = file.path(output_dir, "Figure2H_compact_horizontal_bar_sessionInfo.txt"))

message("Wrote Figure2H compact horizontal bar outputs to: ", normalizePath(output_dir))
