#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
})

out_dir <- "outputs/CellChat_Q8_Q11_vertical_combined"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)

source_levels <- c("Naive CD4 T", "Memory CD4 T", "Th17-like", "Treg", "Activated CD4 T")
source_labels <- c("Naive\nCD4 T", "Memory\nCD4 T", "Th17-like", "Treg", "Activated\nCD4 T")

read_panel <- function(model) {
  base <- file.path(sprintf("outputs/CellChat_%s_donor_significance", model), "source_data")
  scores <- fread(file.path(base, "CellChat_donor_level_total_strength.csv"))
  stats <- fread(file.path(base, sprintf("CellChat_%s_vs_Healthy_donor_level_statistics.csv", model)))
  scores[, model := model]
  stats[, model := model]
  list(scores = scores, stats = stats)
}

q11 <- read_panel("Q11")
q8 <- read_panel("Q8")
all_scores <- rbindlist(list(q11$scores, q8$scores), use.names = TRUE, fill = TRUE)
all_stats <- rbindlist(list(q11$stats, q8$stats), use.names = TRUE, fill = TRUE)
fwrite(all_scores, file.path(out_dir, "source_data", "CellChat_Q8_Q11_donor_level_total_strength.csv"))
fwrite(all_stats, file.path(out_dir, "source_data", "CellChat_Q8_Q11_donor_level_statistics.csv"))

theme_panel <- function() {
  theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.30, colour = "#222222"),
      axis.text = element_text(colour = "#111111"),
      axis.text.x = element_text(size = 6.4, face = "bold", lineheight = 0.92),
      axis.text.y = element_text(size = 6.5),
      axis.title.y = element_text(size = 6.8, margin = margin(r = 4)),
      legend.position = "top",
      legend.justification = "left",
      legend.title = element_blank(),
      legend.text = element_text(size = 6.7),
      legend.key.height = unit(3.5, "mm"),
      legend.key.width = unit(4.5, "mm"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 8, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = 6.7, colour = "#4D4D4D", hjust = 0),
      plot.margin = margin(3, 5, 2, 4)
    )
}

make_panel <- function(model, scores, stats, disease_colour) {
  conditions <- c("Healthy", model)
  pal <- c(Healthy = "#3569A8")
  pal[model] <- disease_colour

  summ <- scores[, .(
    mean_strength = mean(donor_total_strength),
    sem = sd(donor_total_strength) / sqrt(.N)
  ), by = .(source, condition)]
  summ[, source_index := match(source, source_levels)]
  summ[, x := source_index + fifelse(condition == "Healthy", -0.18, 0.18)]
  scores <- copy(scores)
  scores[, source_index := match(source, source_levels)]
  scores[, x := source_index + fifelse(condition == "Healthy", -0.18, 0.18)]
  setorder(scores, source_index, condition, donor_id)
  scores[, x_jitter := x + if (.N == 1L) 0 else seq(-0.050, 0.050, length.out = .N),
         by = .(source, condition)]

  br <- merge(stats[, .(source, q_two_sided, significance)],
              scores[, .(y_max = max(donor_total_strength)), by = source], by = "source")
  br[, source_index := match(source, source_levels)]
  panel_max <- max(scores$donor_total_strength)
  br[, y := y_max * 1.15 + panel_max * 0.028]
  br[, label := paste0(significance, "  q=", formatC(q_two_sided, format = "f", digits = 3))]
  y_limit <- max(br$y) * 1.14

  ggplot() +
    geom_col(data = summ, aes(x, mean_strength, fill = condition), width = 0.30,
             colour = "#222222", linewidth = 0.25) +
    geom_errorbar(data = summ,
                  aes(x, ymin = pmax(0, mean_strength - sem), ymax = mean_strength + sem),
                  width = 0.09, linewidth = 0.35, colour = "#222222") +
    geom_point(data = scores, aes(x_jitter, donor_total_strength, fill = condition),
               shape = 21, size = 1.35, stroke = 0.28, colour = "white", alpha = 0.94) +
    geom_segment(data = br,
                 aes(x = source_index - 0.18, xend = source_index + 0.18, y = y, yend = y),
                 linewidth = 0.35, colour = "#222222") +
    geom_segment(data = br,
                 aes(x = source_index - 0.18, xend = source_index - 0.18,
                     y = y - 0.012 * y_limit, yend = y),
                 linewidth = 0.35, colour = "#222222") +
    geom_segment(data = br,
                 aes(x = source_index + 0.18, xend = source_index + 0.18,
                     y = y - 0.012 * y_limit, yend = y),
                 linewidth = 0.35, colour = "#222222") +
    geom_text(data = br, aes(source_index, y + 0.028 * y_limit, label = label),
              family = "Helvetica", size = 2.15, colour = "#111111") +
    scale_fill_manual(values = pal, breaks = conditions) +
    scale_x_continuous(breaks = seq_along(source_levels), labels = source_labels) +
    scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0))) +
    labs(
      title = sprintf("%s versus Healthy ASPC", model),
      subtitle = "Mean +/- SEM; points represent biological donors",
      x = NULL,
      y = "CellChat-informed total\ninteraction strength"
    ) +
    theme_panel()
}

p_q11 <- make_panel("Q11", q11$scores, q11$stats, "#D8902F")
p_q8 <- make_panel("Q8", q8$scores, q8$stats, "#D85A4F")

combined <- (p_q11 / p_q8) +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Donor-level CD4-to-ASPC interaction strength",
    caption = paste0(
      "Exact two-sided donor-label permutation tests; BH correction across five CD4 states within each comparison.\n",
      "Panel-specific fixed LR universes and y-axis scales; absolute Q8 and Q11 strengths are not directly comparable."
    ),
    theme = theme(
      text = element_text(family = "Helvetica", colour = "#111111"),
      plot.title = element_text(size = 9, face = "bold", hjust = 0, margin = margin(b = 2)),
      plot.caption = element_text(size = 5.6, colour = "#555555", hjust = 0, lineheight = 0.95,
                                  margin = margin(t = 3))
    )
  )

stem <- file.path(out_dir, "CellChat_Q11_Q8_vs_Healthy_vertical_combined")
width_mm <- 122
height_mm <- 168

svglite::svglite(paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(combined)
dev.off()

if (isTRUE(getOption("cellchat.use_cairo_pdf", FALSE))) {
  cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4,
            height = height_mm / 25.4, family = "Helvetica", onefile = TRUE)
} else if (capabilities("aqua")) {
  quartz(type = "pdf", file = paste0(stem, ".pdf"),
         width = width_mm / 25.4, height = height_mm / 25.4,
         family = "Helvetica")
} else {
  pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Helvetica", useDingbats = FALSE)
}
print(combined)
dev.off()

agg_tiff(paste0(stem, ".tiff"), width = width_mm / 25.4, height = height_mm / 25.4,
         units = "in", res = 600, compression = "lzw", background = "white")
print(combined)
dev.off()

agg_png(paste0(stem, "_preview.png"), width = width_mm / 25.4, height = height_mm / 25.4,
        units = "in", res = 300, background = "white")
print(combined)
dev.off()

writeLines(c(
  "CellChat Q11/Q8 vertical combined bar plot QA notes",
  "Core conclusion: show donor-level CD4-to-ASPC total interaction strength for Q11 and Q8 versus Healthy in one compact vertical figure.",
  "Layout: Q11 above Q8; 122 x 168 mm narrow portrait format.",
  "Data/statistics: unchanged from the final standalone Q11 and Q8 analyses.",
  "Bars: mean; error bars: SEM; points: biological donors.",
  "Tests: exact two-sided donor-label permutation; BH correction across five CD4 states within each comparison.",
  "Panel-specific LR universes and y-axis scales are retained; cross-panel absolute heights are not directly comparable.",
  "No observations or CD4 states were excluded during figure assembly."
), file.path(out_dir, "CellChat_Q11_Q8_vertical_combined_QA_notes.txt"))

writeLines(capture.output(sessionInfo()), file.path(out_dir, "CellChat_Q11_Q8_vertical_combined_sessionInfo.txt"))
message("Done: ", normalizePath(out_dir))
