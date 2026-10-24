#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
healthy_path <- if (length(args) >= 1L) args[[1]] else
  "outputs/Figure2H_Q8Q11_revision/delivery/source_data/Figure2H_Q8_Q11_Hallmark_GSEA_joint_BH.csv"
direct_path <- if (length(args) >= 2L) args[[2]] else
  "outputs/Figure2H_Q8Q11_revision/Q8_vs_Q11_direct/delivery/source_data/Figure2H_Q8_vs_Q11_Hallmark_GSEA.csv"
output_dir <- if (length(args) >= 3L) args[[3]] else
  "outputs/Figure2H_disease_up_two_layouts"

healthy_path <- normalizePath(healthy_path, mustWork = TRUE)
direct_path <- normalizePath(direct_path, mustWork = TRUE)
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

# Figure contract
# Core conclusion: show the strongest disease-up Hallmark programs in Q8 and Q11,
# and the four best-supported pathways in each direction of the paired Q8-vs-Q11 comparison.
# Archetype: quantitative small-multiple grid.
# Selection: Q8/Q11 NES > 0, top 8 by joint q; direct comparison top 4 by q
# independently for NES > 0 (Q8-high) and NES < 0 (Q11-high).
# Statistics: joint BH across 92 Q8/Q11-vs-Healthy tests; direct BH across 17
# eligible pathways from three paired donors.
# Export: two layouts, editable SVG/PDF plus 600-dpi TIFF and PNG preview.

healthy <- fread(healthy_path)
direct <- fread(direct_path)

required_healthy <- c("comparison", "pathway", "pathway_label", "NES", "pval", "q_joint")
required_direct <- c("pathway", "pathway_label", "NES", "pval", "q")
if (length(setdiff(required_healthy, names(healthy)))) {
  stop("Healthy-reference GSEA table is missing: ",
       paste(setdiff(required_healthy, names(healthy)), collapse = ", "))
}
if (length(setdiff(required_direct, names(direct)))) {
  stop("Direct Q8-vs-Q11 GSEA table is missing: ",
       paste(setdiff(required_direct, names(direct)), collapse = ", "))
}

wrap_label <- function(x, width = 25L) {
  vapply(x, function(z) paste(strwrap(z, width = width), collapse = "\n"), character(1))
}

format_q <- function(q) {
  ifelse(q < 0.001, "q<0.001", sprintf("q=%.3f", q))
}

select_disease_up <- function(comp, display, n = 8L) {
  z <- copy(healthy[comparison == comp & is.finite(NES) & NES > 0 & is.finite(q_joint)])
  setorder(z, q_joint, -NES, pathway)
  z <- head(z, n)
  if (nrow(z) < n) stop("Fewer than ", n, " disease-up pathways for ", comp)
  z[, `:=`(
    comparison_display = display,
    direction_display = "Disease-high",
    q_value = q_joint,
    selection_rank = seq_len(.N),
    selection_rule = paste0("NES > 0; top ", n, " by joint q, then NES")
  )]
  z[]
}

select_direct <- function(n = 4L) {
  pos <- copy(direct[is.finite(NES) & NES > 0 & is.finite(q)])
  neg <- copy(direct[is.finite(NES) & NES < 0 & is.finite(q)])
  setorder(pos, q, -NES, pathway)
  setorder(neg, q, NES, pathway)
  pos <- head(pos, n)
  neg <- head(neg, n)
  if (nrow(pos) < n || nrow(neg) < n) {
    stop("Fewer than ", n, " eligible pathways in one direct-comparison direction")
  }
  pos[, `:=`(
    comparison_display = "Q8 vs Q11",
    direction_display = "Q8-high",
    q_value = q,
    selection_rank = seq_len(.N),
    selection_rule = paste0("NES > 0; top ", n, " by direct-comparison q, then NES")
  )]
  neg[, `:=`(
    comparison_display = "Q8 vs Q11",
    direction_display = "Q11-high",
    q_value = q,
    selection_rank = seq_len(.N),
    selection_rule = paste0("NES < 0; top ", n, " by direct-comparison q, then absolute NES")
  )]
  rbindlist(list(pos, neg), fill = TRUE)
}

plot_dt <- rbindlist(list(
  select_disease_up("Q8 vs Healthy", "Q8 vs Healthy", 8L),
  select_disease_up("Q11 vs Healthy", "Q11 vs Healthy", 8L),
  select_direct(4L)
), fill = TRUE)

plot_dt[, `:=`(
  fdr_significant = q_value < 0.05,
  fdr_status = fifelse(q_value < 0.05, "q < 0.05", "q >= 0.05"),
  q_label = format_q(q_value)
)]

for (col in c("leadingEdge")) {
  if (col %in% names(plot_dt) && is.list(plot_dt[[col]])) {
    plot_dt[, (col) := vapply(get(col), paste, collapse = ";", FUN.VALUE = character(1))]
  }
}

horizontal_palette <- data.table(
  comparison_display = c("Q8 vs Healthy", "Q11 vs Healthy", "Q8 vs Q11"),
  strip_fill = c("#EBC6DA", "#FFF2A8", "#C9E8EF"),
  mark_colour = c("#A94F7A", "#A98312", "#2F8797")
)
vertical_palette <- data.table(
  comparison_display = c("Q8 vs Healthy", "Q11 vs Healthy", "Q8 vs Q11"),
  strip_fill = c("#A8D18D", "#F5A675", "#F06C70"),
  mark_colour = c("#527F38", "#C85E29", "#A92F43")
)
fwrite(horizontal_palette, file.path(source_dir, "Figure2H_horizontal_palette.csv"))
fwrite(vertical_palette, file.path(source_dir, "Figure2H_vertical_palette.csv"))

font_family <- "Arial"
font_pt <- 8
geom_text_size <- font_pt / ggplot2::.pt

prepare_panel_data <- function(comp, wrap_width) {
  d <- copy(plot_dt[comparison_display == comp])
  if (comp == "Q8 vs Q11") {
    display_order <- c(
      d[direction_display == "Q8-high"][order(selection_rank), pathway_label],
      d[direction_display == "Q11-high"][order(selection_rank), pathway_label]
    )
  } else {
    display_order <- d[order(selection_rank), pathway_label]
  }
  d[, pathway_plot := wrap_label(pathway_label, width = wrap_width)]
  display_order_wrapped <- wrap_label(display_order, width = wrap_width)
  d[, pathway_plot := factor(pathway_plot, levels = rev(display_order_wrapped))]

  xmin_data <- min(c(0, d$NES), na.rm = TRUE)
  xmax_data <- max(c(0, d$NES), na.rm = TRUE)
  data_span <- xmax_data - xmin_data
  if (!is.finite(data_span) || data_span <= 0) data_span <- 1
  d[, label_x := NES + fifelse(NES >= 0, 0.045, 0.060) * data_span]
  d[, label_hjust := 0]
  d[, label_vjust := fifelse(NES >= 0, 0.5, -0.55)]

  if (xmin_data < 0) {
    x_limits <- c(xmin_data - 0.28 * data_span, xmax_data + 0.28 * data_span)
  } else {
    x_limits <- c(0, xmax_data + 0.42 * data_span)
  }
  list(data = d, limits = x_limits)
}

make_panel <- function(comp, palette_dt, wrap_width = 24L,
                       show_x_title = TRUE, vertical_layout = FALSE) {
  prep <- prepare_panel_data(comp, wrap_width)
  d <- prep$data
  palette_row <- palette_dt[comparison_display == comp]
  strip_fill <- palette_row$strip_fill
  mark_colour <- palette_row$mark_colour

  strip_text <- switch(
    comp,
    "Q8 vs Healthy" = "Q8 vs Healthy ASPC-like\nDisease-up pathways (top 8)",
    "Q11 vs Healthy" = "Q11 vs Healthy ASPC-like\nDisease-up pathways (top 8)",
    "Q8 vs Q11" = "Q8 vs Q11 (paired donors)\nTop 4 Q8-high and top 4 Q11-high"
  )
  d[, strip_label := strip_text]

  p <- ggplot(d, aes(NES, pathway_plot)) +
    geom_vline(xintercept = 0, linewidth = 0.42, colour = "#222222") +
    geom_segment(
      aes(x = 0, xend = NES, yend = pathway_plot),
      linewidth = 0.9, colour = mark_colour, lineend = "round"
    ) +
    geom_point(
      data = d[fdr_significant == TRUE],
      shape = 21, size = 2.9, stroke = 0.75,
      colour = mark_colour, fill = mark_colour
    ) +
    geom_point(
      data = d[fdr_significant == FALSE],
      shape = 21, size = 2.9, stroke = 0.85,
      colour = mark_colour, fill = "white"
    ) +
    geom_text(
      aes(x = label_x, label = q_label, hjust = label_hjust, vjust = label_vjust),
      family = font_family, fontface = "bold", size = geom_text_size,
      colour = "#111111", lineheight = 0.95
    ) +
    facet_wrap(~strip_label) +
    scale_x_continuous(
      breaks = pretty(prep$limits, n = if (vertical_layout) 6 else 4),
      expand = expansion(mult = 0)
    ) +
    coord_cartesian(xlim = prep$limits, clip = "off") +
    labs(x = if (show_x_title) "Normalized enrichment score (NES)" else NULL, y = NULL) +
    theme_classic(base_family = font_family, base_size = font_pt) +
    theme(
      text = element_text(family = font_family, face = "bold", size = font_pt, colour = "#111111"),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.x = element_line(linewidth = 0.45, colour = "#222222"),
      axis.ticks.x = element_line(linewidth = 0.4, colour = "#222222"),
      axis.text = element_text(size = font_pt, face = "bold", colour = "#111111"),
      axis.text.y = element_text(lineheight = 0.9, margin = margin(r = 3)),
      axis.title.x = element_text(size = font_pt, face = "bold", margin = margin(t = 4)),
      strip.background = element_rect(fill = strip_fill, colour = NA),
      strip.text = element_text(
        size = font_pt, face = "bold", colour = "#111111",
        lineheight = 0.92, margin = margin(t = 4, r = 4, b = 4, l = 4)
      ),
      panel.grid = element_blank(),
      plot.margin = if (vertical_layout) margin(4, 12, 4, 4) else margin(4, 14, 4, 4)
    )

  if (comp == "Q8 vs Q11") {
    p <- p +
      geom_hline(yintercept = 4.5, linewidth = 0.36, linetype = 2, colour = "#777777") +
      labs(caption = "Left: Q11-high     Right: Q8-high") +
      theme(plot.caption = element_text(
        size = font_pt, face = "bold", hjust = 0.5,
        colour = "#333333", margin = margin(t = 3)
      ))
  }
  p
}

horizontal <- (
  make_panel("Q8 vs Healthy", horizontal_palette, wrap_width = 19L,
             show_x_title = TRUE, vertical_layout = FALSE) |
    make_panel("Q11 vs Healthy", horizontal_palette, wrap_width = 19L,
               show_x_title = TRUE, vertical_layout = FALSE) |
    make_panel("Q8 vs Q11", horizontal_palette, wrap_width = 19L,
               show_x_title = TRUE, vertical_layout = FALSE)
) +
  plot_layout(widths = c(1, 1, 1)) +
  plot_annotation(
    title = "Hallmark pathways enriched in disease-associated ASPC states",
    caption = paste0(
      "Q8/Q11 panels show disease-up pathways only. Filled points: q < 0.05; open points: q >= 0.05.\n",
      "Q8 vs Q11 is exploratory (3 paired donors)."
    ),
    theme = theme(
      text = element_text(family = font_family, face = "bold", size = font_pt, colour = "#111111"),
      plot.title = element_text(size = font_pt, face = "bold", hjust = 0),
      plot.caption = element_text(size = font_pt, face = "bold", colour = "#444444", hjust = 0)
    )
  )

vertical <- (
  make_panel("Q8 vs Healthy", vertical_palette, wrap_width = 34L,
             show_x_title = FALSE, vertical_layout = TRUE) /
    make_panel("Q11 vs Healthy", vertical_palette, wrap_width = 34L,
               show_x_title = FALSE, vertical_layout = TRUE) /
    make_panel("Q8 vs Q11", vertical_palette, wrap_width = 34L,
               show_x_title = TRUE, vertical_layout = TRUE)
) +
  plot_layout(heights = c(1, 1, 1.08)) +
  plot_annotation(
    title = "Hallmark pathways enriched in disease-associated ASPC states",
    caption = paste0(
      "Q8/Q11 panels show disease-up pathways only. Filled points: q < 0.05; open points: q >= 0.05.\n",
      "Q8 vs Q11 is exploratory (3 paired donors)."
    ),
    theme = theme(
      text = element_text(family = font_family, face = "bold", size = font_pt, colour = "#111111"),
      plot.title = element_text(size = font_pt, face = "bold", hjust = 0),
      plot.caption = element_text(size = font_pt, face = "bold", colour = "#444444", hjust = 0)
    )
  )

open_pdf <- function(path, width_mm, height_mm) {
  if (capabilities("aqua")) {
    grDevices::quartz(
      type = "pdf", file = path,
      width = width_mm / 25.4, height = height_mm / 25.4,
      family = font_family, bg = "white"
    )
  } else if (capabilities("cairo")) {
    grDevices::cairo_pdf(
      path, width = width_mm / 25.4, height = height_mm / 25.4,
      family = font_family, bg = "white"
    )
  } else {
    grDevices::pdf(
      path, width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Helvetica", bg = "white", useDingbats = FALSE
    )
  }
}

save_plot <- function(plot, stem, width_mm, height_mm) {
  open_pdf(paste0(stem, ".pdf"), width_mm, height_mm)
  print(plot)
  dev.off()

  svglite::svglite(
    paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4,
    bg = "white", system_fonts = list(sans = font_family)
  )
  print(plot)
  dev.off()

  ragg::agg_tiff(
    paste0(stem, ".tiff"), width = width_mm, height = height_mm,
    units = "mm", res = 600, compression = "lzw", background = "white"
  )
  print(plot)
  dev.off()

  ragg::agg_png(
    paste0(stem, "_preview.png"), width = width_mm, height = height_mm,
    units = "mm", res = 300, background = "white"
  )
  print(plot)
  dev.off()
}

horizontal_stem <- file.path(output_dir, "Figure2H_disease_up_horizontal_pink_yellow_blue")
vertical_stem <- file.path(output_dir, "Figure2H_disease_up_vertical_green_orange_red")
save_plot(horizontal, horizontal_stem, width_mm = 270, height_mm = 118)
save_plot(vertical, vertical_stem, width_mm = 183, height_mm = 225)

export_dt <- copy(plot_dt)
export_dt[, comparison_order := match(
  comparison_display, c("Q8 vs Healthy", "Q11 vs Healthy", "Q8 vs Q11")
)]
export_dt[, direction_order := match(
  direction_display, c("Disease-high", "Q8-high", "Q11-high")
)]
setorder(export_dt, comparison_order, direction_order, selection_rank)
export_dt[, c("comparison_order", "direction_order") := NULL]
fwrite(export_dt, file.path(source_dir, "Figure2H_disease_up_selected_pathways.csv"))

selection_summary <- export_dt[, .(
  n_displayed = .N,
  n_q_lt_0_05 = sum(q_value < 0.05),
  minimum_q = min(q_value),
  maximum_q = max(q_value),
  selection_rule = selection_rule[1]
), by = .(comparison_display, direction_display)]
fwrite(selection_summary, file.path(source_dir, "Figure2H_disease_up_selection_summary.csv"))

writeLines(c(
  "Figure 2H disease-up pathway revision QA notes",
  "Q8 vs Healthy: only NES > 0 pathways retained; top 8 by joint BH q, then NES.",
  "Q11 vs Healthy: only NES > 0 pathways retained; top 8 by joint BH q, then NES.",
  "Q8 vs Q11: top 4 by q selected independently for NES > 0 (Q8-high) and NES < 0 (Q11-high).",
  "Q8/Q11 Healthy-reference q values use joint BH correction across 92 tested pathway-comparison results.",
  "Direct Q8-vs-Q11 q values use BH correction across 17 eligible Hallmark pathways; only 3 paired donors are available.",
  "No NES, p value or q value was changed; no simulated values were used.",
  "Filled endpoints indicate q < 0.05; open endpoints indicate q >= 0.05; exact q values are printed.",
  "Version 1 palette/layout: Q8 pink, Q11 yellow, Q8-vs-Q11 blue; panels run left to right.",
  "Version 2 palette/layout: Q8 green, Q11 orange, Q8-vs-Q11 red; panels run top to bottom.",
  "All visible text is Arial Bold 8 pt; horizontal canvas is 270 x 118 mm to preserve the requested physical font size.",
  "Vertical canvas is 183 x 225 mm and fits a Nature-style double-column width.",
  "PDF/SVG are vector; TIFF is 600 dpi."
), file.path(output_dir, "Figure2H_disease_up_two_layouts_QA_notes.txt"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2H_disease_up_sessionInfo.txt"))
message("Figure 2H disease-up two-layout revision complete: ", normalizePath(output_dir))
