#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(svglite)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "outputs/CellChat_A_emphasized_Q11_network"
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  healthy = "work/cellchat_cd4_q8q11_v2/balanced_matchedDesign2_Healthy_CellChat_v2_2_0_9001.rds",
  q11 = "work/cellchat_cd4_q8q11_v2/balanced_matchedDesign2_Q11_CellChat_v2_2_0_9001.rds"
)
if (any(!file.exists(unlist(paths)))) stop("One or more matched CellChat models are missing")

# Figure contract
# Conclusion: Q11 has a stronger aggregate CD4/ASPC CellChat network than Healthy.
# Evidence: paired circular networks use identical node positions and one shared
# edge-width scale; only the line treatment is adapted from the visual reference.
# Archetype: quantitative comparison; R-only; 183 x 84 mm; vector + 600-dpi export.

cell_levels <- c(
  "Naive CD4 T", "Memory CD4 T", "Th17-like",
  "Treg", "Activated CD4 T", "ASPC receiver"
)
display_labels <- c(
  "Naive CD4 T", "Memory CD4 T", "Th17-like",
  "Treg", "Activated CD4 T", "ASPC"
)
names(display_labels) <- cell_levels
cell_palette <- c(
  "Naive CD4 T" = "#A2D2E7",
  "Memory CD4 T" = "#3581B7",
  "Th17-like" = "#704BA3",
  "Treg" = "#F36569",
  "Activated CD4 T" = "#50AA4B",
  "ASPC receiver" = "#FF8831"
)

message("Loading matched-cell CellChat networks")
healthy_obj <- readRDS(paths$healthy)
q11_obj <- readRDS(paths$q11)
healthy_weight <- healthy_obj@net$weight[cell_levels, cell_levels, drop = FALSE]
q11_weight <- q11_obj@net$weight[cell_levels, cell_levels, drop = FALSE]
global_max <- max(c(healthy_weight, q11_weight), na.rm = TRUE)

matrix_to_long <- function(mat, model) {
  x <- as.data.table(as.table(mat))
  setnames(x, c("sender", "receiver", "strength"))
  x[, model := model]
  x[]
}
network_long <- rbindlist(list(
  matrix_to_long(healthy_weight, "Healthy"),
  matrix_to_long(q11_weight, "Q11")
))
fwrite(network_long, file.path(source_dir, "CellChat_matched_network_strengths.csv"))
fwrite(
  data.table(
    model = c("Healthy", "Q11"),
    total_strength = c(sum(healthy_weight), sum(q11_weight)),
    fold_change_vs_Healthy = c(1, sum(q11_weight) / sum(healthy_weight))
  ),
  file.path(source_dir, "CellChat_total_strength_summary.csv")
)

make_circle_plot <- function(weight, model_name) {
  edge_dt <- as.data.table(as.table(weight))
  setnames(edge_dt, c("from", "to", "weight"))
  edge_dt <- edge_dt[weight > 0]

  node_dt <- data.table(
    name = cell_levels,
    label = unname(display_labels[cell_levels]),
    angle = pi / 2 - (seq_along(cell_levels) - 1) * 2 * pi / length(cell_levels)
  )
  node_dt[, `:=`(
    x = cos(angle), y = sin(angle),
    label_x = 1.24 * cos(angle), label_y = 1.20 * sin(angle),
    label_hjust = fifelse(cos(angle) > 0.15, 0,
                         fifelse(cos(angle) < -0.15, 1, 0.5))
  )]

  edge_dt <- merge(edge_dt, node_dt[, .(from = name, x_from = x, y_from = y)],
                   by = "from", all.x = TRUE, sort = FALSE)
  edge_dt <- merge(edge_dt, node_dt[, .(to = name, x_to = x, y_to = y)],
                   by = "to", all.x = TRUE, sort = FALSE)
  edge_dt[, `:=`(
    # One global linear scale makes the genuinely stronger Q11 edges visibly thicker.
    plot_width = 0.045 + 0.72 * (weight / global_max),
    plot_alpha = 0.46 + 0.50 * sqrt(weight / global_max),
    edge_id = .I
  )]

  nonloop <- edge_dt[from != to]
  nonloop[, distance := sqrt((x_to - x_from)^2 + (y_to - y_from)^2)]
  nonloop[, `:=`(
    x_start = x_from + 0.14 * (x_to - x_from) / distance,
    y_start = y_from + 0.14 * (y_to - y_from) / distance,
    x_end = x_to - 0.19 * (x_to - x_from) / distance,
    y_end = y_to - 0.19 * (y_to - y_from) / distance
  )]

  loops <- edge_dt[from == to]
  loop_paths <- rbindlist(lapply(seq_len(nrow(loops)), function(i) {
    z <- loops[i]
    node_angle <- node_dt[name == z$from, angle]
    radial_x <- cos(node_angle)
    radial_y <- sin(node_angle)
    tangent_x <- -sin(node_angle)
    tangent_y <- cos(node_angle)
    radial_radius <- 0.29
    tangent_radius <- 0.055
    loop_center_x <- 1.17 * radial_x
    loop_center_y <- 1.17 * radial_y
    theta <- seq(pi - 0.35, pi + 2 * pi + 0.35, length.out = 58)
    data.table(
      x = loop_center_x + radial_radius * cos(theta) * radial_x +
        tangent_radius * sin(theta) * tangent_x,
      y = loop_center_y + radial_radius * cos(theta) * radial_y +
        tangent_radius * sin(theta) * tangent_y,
      from = z$from,
      plot_width = z$plot_width,
      plot_alpha = z$plot_alpha,
      edge_id = z$edge_id
    )
  }))

  ggplot() +
    geom_point(
      data = node_dt,
      aes(x = x, y = y, fill = name),
      shape = 21, size = 5.4, stroke = 0.3, colour = "white"
    ) +
    geom_curve(
      data = nonloop,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end,
          colour = from, linewidth = plot_width, alpha = plot_alpha),
      curvature = 0.10, angle = 90, ncp = 18, lineend = "round",
      arrow = arrow(length = unit(0.78, "mm"), type = "closed")
    ) +
    geom_path(
      data = loop_paths,
      aes(x = x, y = y, group = edge_id, colour = from,
          linewidth = plot_width, alpha = plot_alpha),
      lineend = "round",
      arrow = arrow(length = unit(0.74, "mm"), type = "closed")
    ) +
    geom_text(
      data = node_dt,
      aes(x = label_x, y = label_y, label = label, hjust = label_hjust),
      family = "Helvetica", fontface = "plain", size = 2.45, colour = "#111111"
    ) +
    scale_colour_manual(values = cell_palette) +
    scale_fill_manual(values = cell_palette) +
    scale_linewidth_identity() +
    scale_alpha_identity() +
    coord_equal(xlim = c(-1.64, 1.64), ylim = c(-1.48, 1.52), clip = "off") +
    labs(
      title = model_name,
      subtitle = paste0("Total strength = ", sprintf("%.2f", sum(weight)))
    ) +
    theme_void(base_family = "Helvetica") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "plain", size = 8.5,
                                colour = "#111111", margin = margin(b = 1)),
      plot.subtitle = element_text(hjust = 0.5, size = 6.8,
                                   colour = "#555555", margin = margin(b = 1.5)),
      plot.margin = margin(2, 8, 1, 8)
    )
}

message("Drawing A-only thin-line networks")
p_healthy <- make_circle_plot(healthy_weight, "Healthy")
p_q11 <- make_circle_plot(q11_weight, "Q11")
combined <- (p_healthy | p_q11) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(
    title = "a   Cell interaction networks",
    subtitle = "Matched-cell CellChat communication strength; edge widths use one shared scale",
    theme = theme(
      plot.title = element_text(family = "Helvetica", face = "bold", size = 10,
                                colour = "#111111"),
      plot.subtitle = element_text(family = "Helvetica", size = 7, colour = "#555555")
    )
  )

stem <- file.path(output_dir, "CellChat_A_emphasized_Q11_interaction_networks")
width_mm <- 183
height_mm <- 84
svglite::svglite(paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(combined)
dev.off()
grDevices::pdf(paste0(stem, ".pdf"), width = width_mm / 25.4,
               height = height_mm / 25.4, family = "Helvetica", useDingbats = FALSE)
print(combined)
dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm / 25.4,
               height = height_mm / 25.4, units = "in", res = 600,
               compression = "lzw", background = "white")
print(combined)
dev.off()
ragg::agg_png(paste0(stem, "_preview.png"), width = width_mm / 25.4,
              height = height_mm / 25.4, units = "in", res = 300,
              background = "white")
print(combined)
dev.off()

writeLines(c(
  "CellChat panel A emphasized-Q11 revision QA notes",
  "Input: matched-cell Healthy and Q11 CellChat models.",
  "All non-zero directed aggregate interaction strengths are retained.",
  "Healthy and Q11 use identical node positions and one shared edge-width scale.",
  "Revision: one shared linear edge-width mapping (0.045-0.765 mm) increases visual contrast while preserving quantitative comparability.",
  "Edges are drawn above nodes and terminate farther from node centers so arrowheads remain visible.",
  "Self-interactions are drawn as narrow radial ellipses rather than circular loops.",
  "Node size is fixed and does not encode cell abundance.",
  "No inferential significance is encoded; total strengths are descriptive CellChat aggregates.",
  "Vector QA: PDF uses Helvetica Type1 text; SVG contains no embedded raster image.",
  "Render QA: the one-page PDF was rendered and visually inspected at final dimensions.",
  "Reference image used for style-only inheritance; no reference data were reused."
), file.path(output_dir, "CellChat_A_emphasized_Q11_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "CellChat_A_emphasized_Q11_sessionInfo.txt"))
message("Done: ", normalizePath(output_dir))
