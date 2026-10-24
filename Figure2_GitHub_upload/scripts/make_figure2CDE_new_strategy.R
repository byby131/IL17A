#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(hdf5r)
  library(Matrix)
  library(Seurat)
  library(svglite)
  library(ragg)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/make_figure2CDE_new_strategy.R <Tcell_secondary.csv.gz> <all.H5Seurat> <output_dir>")
}

tcell_path <- normalizePath(args[[1]], mustWork = TRUE)
h5_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

base_family <- "Helvetica"
set.seed(264508)

state_order <- c(
  "Naive CD4 T", "Memory CD4 T", "Th17-like", "Treg", "Activated CD4 T"
)
state_palette <- c(
  "Naive CD4 T" = "#A2D2E7",
  "Memory CD4 T" = "#3581B7",
  "Th17-like" = "#F36569",
  "Treg" = "#704BA3",
  "Activated CD4 T" = "#FF8831"
)

cluster_annotation <- data.table(
  tcell_res_1 = as.character(0:9),
  final_state = c(
    "Th17-like", "Treg", "Excluded: CD8/cytotoxic-like", "Memory CD4 T",
    "Activated CD4 T", "Excluded: CD8/gamma-delta-like", "Naive CD4 T",
    "Activated CD4 T", "Naive CD4 T", "Activated CD4 T"
  ),
  disposition = c(
    "retained", "retained", "excluded", "retained", "retained",
    "excluded", "retained", "retained", "retained", "retained"
  ),
  evidence = c(
    "CCR6/KLRB1/RORA/RORC/CCL20; low CD8/gamma-delta program",
    "FOXP3/IL2RA/CTLA4/TIGIT/IKZF2/TNFRSF4",
    "CD8A/CD8B/CCL5/NKG7 cytotoxic program",
    "IL7R/LTB with low naive, CD8 and gamma-delta programs",
    "CD69/FOS/JUN/ICOS/HLA-DRA activation program",
    "CD8A/CD8B plus TRDC/TRGC and cytotoxic program",
    "CCR7/TCF7/LEF1/SELL/MAL",
    "CD69/FOS/JUN activation program; donor-supported",
    "CCR7/TCF7/LEF1/SELL; merged with cluster 6; no cycling evidence",
    "FOS/JUN/CD69 activation program; small but donor-supported subset"
  )
)

dt <- fread(tcell_path)
dt[, tcell_res_1 := as.character(tcell_res_1)]
dt <- merge(dt, cluster_annotation, by = "tcell_res_1", all.x = TRUE, sort = FALSE)
if (anyNA(dt$final_state)) stop("Unannotated T-cell clusters remain")

excluded <- dt[disposition == "excluded"]
cd4 <- dt[disposition == "retained"]
cd4[, final_state := factor(final_state, levels = state_order)]
if (nrow(cd4) == 0L) stop("No retained CD4 cells")

# Re-embed only the evidence-supported CD4 compartment using the published
# Harmony latent space. This changes visualization geometry, not labels.
h5 <- H5File$new(h5_path, "r")
on.exit(h5$close_all(), add = TRUE)
harmony <- h5[["reductions"]][["harmony"]][["cell.embeddings"]]$read()
harmony_cd4 <- harmony[cd4$source_row, 1:30, drop = FALSE]
rownames(harmony_cd4) <- cd4$cell_id
colnames(harmony_cd4) <- paste0("HARMONY_", seq_len(ncol(harmony_cd4)))

dummy <- sparseMatrix(
  i = rep.int(1L, nrow(cd4)), j = seq_len(nrow(cd4)), x = rep.int(1, nrow(cd4)),
  dims = c(1L, nrow(cd4)), dimnames = list("placeholder", cd4$cell_id)
)
obj <- CreateSeuratObject(counts = dummy, min.cells = 0, min.features = as.integer(0))
obj[["harmony"]] <- CreateDimReducObject(
  embeddings = harmony_cd4, key = "HARMONY_", assay = "RNA"
)
obj <- RunUMAP(
  obj, reduction = "harmony", dims = 1:30,
  n.neighbors = 25, min.dist = 0.30, metric = "cosine",
  seed.use = 264508, reduction.name = "cd4_umap", reduction.key = "CD4UMAP_",
  verbose = FALSE
)
umap <- as.data.table(Embeddings(obj, "cd4_umap"), keep.rownames = "cell_id")
setnames(umap, c("cell_id", "CD4_UMAP_1", "CD4_UMAP_2"))
cd4 <- merge(cd4, umap, by = "cell_id", all.x = TRUE, sort = FALSE)
cd4[, final_state := factor(final_state, levels = state_order)]

state_summary <- cd4[, .(
  n_cells = .N,
  n_donors = uniqueN(library_id),
  healthy_cells = sum(health == "Healthy"),
  losc_cells = sum(health == "LS")
), by = final_state]
state_summary[, fraction_percent := 100 * n_cells / sum(n_cells)]

label_positions <- cd4[, {
  center_x <- median(CD4_UMAP_1)
  center_y <- median(CD4_UMAP_2)
  d2 <- (CD4_UMAP_1 - center_x)^2 + (CD4_UMAP_2 - center_y)^2
  .(CD4_UMAP_1 = CD4_UMAP_1[which.min(d2)], CD4_UMAP_2 = CD4_UMAP_2[which.min(d2)])
}, by = final_state]
label_positions[, label := c(
  "Naive", "Memory", "Th17-like", "Treg", "Activated"
)[match(as.character(final_state), state_order)]]

draw_order <- state_summary[order(-n_cells), as.character(final_state)]
cd4[, draw_rank := match(as.character(final_state), draw_order)]
setorder(cd4, draw_rank)

x_range <- range(cd4$CD4_UMAP_1)
y_range <- range(cd4$CD4_UMAP_2)
x_span <- diff(x_range)
y_span <- diff(y_range)
axis_x0 <- x_range[1] + 0.025 * x_span
axis_y0 <- y_range[1] + 0.030 * y_span

p_c <- ggplot(cd4, aes(CD4_UMAP_1, CD4_UMAP_2, color = final_state)) +
  geom_point(size = 0.38, alpha = 0.95, stroke = 0) +
  geom_text_repel(
    data = label_positions,
    aes(label = label),
    color = "#111111",
    family = "Arial", fontface = "bold", size = 8.5 / ggplot2::.pt,
    box.padding = 0.28, point.padding = 0.18,
    min.segment.length = Inf, segment.color = NA,
    seed = 264508, max.overlaps = Inf,
    show.legend = FALSE
  ) +
  annotate("segment", x = axis_x0, xend = axis_x0 + 0.10 * x_span,
           y = axis_y0, yend = axis_y0, linewidth = 0.34, color = "#22272A",
           arrow = arrow(length = unit(1.15, "mm"), type = "closed")) +
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + 0.10 * y_span, linewidth = 0.34,
           color = "#22272A", arrow = arrow(length = unit(1.15, "mm"), type = "closed")) +
  annotate("text", x = axis_x0 + 0.055 * x_span, y = axis_y0 - 0.040 * y_span,
           label = "UMAP1", family = base_family, size = 2.20, color = "#111111") +
  annotate("text", x = axis_x0 - 0.045 * x_span, y = axis_y0 + 0.055 * y_span,
           label = "UMAP2", angle = 90, family = base_family, size = 2.20, color = "#111111") +
  scale_color_manual(values = state_palette, drop = FALSE) +
  coord_equal(
    xlim = x_range + c(-0.04, 0.04) * x_span,
    ylim = y_range + c(-0.04, 0.04) * y_span,
    expand = FALSE, clip = "off"
  ) +
  labs(title = "CD4 T-cell states") +
  theme_void(base_family = base_family) +
  theme(
    plot.title = element_text(size = 8.2, color = "#111111", hjust = 0.5, margin = margin(b = 2)),
    plot.margin = margin(5.5, 5.0, 4.0, 5.0), legend.position = "none"
  ) +
  annotate("text", x = -Inf, y = Inf, label = "C", hjust = -0.15, vjust = 1.15,
           family = base_family, fontface = "bold", size = 3.2, color = "black")

make_dot_data <- function(input, states, genes, modules) {
  available <- intersect(genes, names(input))
  observed <- rbindlist(lapply(available, function(gene) {
    input[as.character(final_state) %chin% states, .(
      n_cells = .N,
      n_expressing = sum(get(gene) > 0),
      pct_exp = 100 * mean(get(gene) > 0),
      avg_exp = mean(get(gene))
    ), by = .(state = as.character(final_state))][, marker := gene]
  }))
  grid <- CJ(state = states, marker = genes, unique = TRUE)
  out <- merge(grid, observed, by = c("state", "marker"), all.x = TRUE)
  out[, available := marker %chin% available]
  out[available == FALSE, `:=`(n_cells = NA_integer_, n_expressing = NA_integer_, pct_exp = NA_real_, avg_exp = NA_real_)]
  out[, avg_exp_scaled := {
    if (all(is.na(avg_exp)) || sd(avg_exp, na.rm = TRUE) == 0 || is.na(sd(avg_exp, na.rm = TRUE))) {
      rep(0, .N)
    } else {
      pmax(-1.5, pmin(1.5, as.numeric(scale(avg_exp))))
    }
  }, by = marker]
  out[, module := modules[marker]]
  out
}

focused_genes <- c(
  "RORC", "CCR6", "KLRB1", "CCL20",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
  "TRDC", "CD8A", "CD8B", "NKG7"
)
focused_modules <- c(
  setNames(rep("Th17-like markers", 4), focused_genes[1:4]),
  setNames(rep("Treg identity", 5), focused_genes[5:9]),
  setNames(rep("Lineage-mixture indicators", 4), focused_genes[10:13])
)
focused_states <- c("Th17-like", "Treg")
# Scale each gene across all retained CD4 states before selecting the two rows.
# With only two rows, row-wise z scores collapse to artificial +/-0.707 values
# and exaggerate small differences between the candidate and Treg populations.
dot_d_all <- make_dot_data(cd4, state_order, focused_genes, focused_modules)
dot_d <- dot_d_all[state %chin% focused_states]
dot_d[, state_display := fifelse(state == "Th17-like", "Th17-like", "Treg")]
dot_d[, state_y := fifelse(state == "Th17-like", 0.915, 0.86)]
dot_d[, marker := factor(marker, levels = focused_genes)]
dot_d[, module := factor(
  module,
  levels = c("Th17-like markers", "Treg identity", "Lineage-mixture indicators")
)]
dot_d[, marker_x := match(as.character(marker), focused_genes)]

omitted_unavailable_th17 <- c("IL23R", "IL17A", "IL17F")

marker_audit <- dcast(
  copy(dot_d), marker + module + available ~ state,
  value.var = c("pct_exp", "avg_exp")
)
marker_audit[, marker := factor(marker, levels = focused_genes)]
setorder(marker_audit, marker)
marker_audit[, marker := as.character(marker)]

p_d <- ggplot(dot_d, aes(marker_x, state_y)) +
  annotate(
    "rect", xmin = 0.45, xmax = 13.55, ymin = 0.82, ymax = 0.95,
    fill = NA, color = "#303638", linewidth = 0.34
  ) +
  annotate(
    "rect", xmin = 0.55, xmax = 4.45, ymin = 0.96, ymax = 0.985,
    fill = "#E36F62", color = NA
  ) +
  annotate(
    "rect", xmin = 4.55, xmax = 9.45, ymin = 0.96, ymax = 0.985,
    fill = "#684079", color = NA
  ) +
  annotate(
    "rect", xmin = 9.55, xmax = 13.45, ymin = 0.96, ymax = 0.985,
    fill = "#AEB5B7", color = NA
  ) +
  annotate("text", x = 2.5, y = 1.025, label = "Th17-like",
           family = base_family, fontface = "plain", size = 1.8, color = "#D95A83") +
  annotate("text", x = 7.0, y = 1.025, label = "Treg",
           family = base_family, fontface = "plain", size = 1.8, color = "#684079") +
  annotate("text", x = 11.5, y = 1.025, label = "Lineage-mixture",
           family = base_family, fontface = "plain", size = 1.8, color = "#52686B") +
  geom_point(
    shape = 21, size = 0.30, fill = "white", color = "#7E8587", stroke = 0.16
  ) +
  geom_point(
    data = dot_d[available == TRUE & !is.na(pct_exp) & pct_exp > 0],
    aes(size = pct_exp, fill = avg_exp_scaled),
    shape = 21, color = "#303638", stroke = 0.22
  ) +
  scale_fill_gradientn(
    colors = c("#EEF3F4", "#F5C0B4", "#E45E54", "#7D1739"),
    values = rescale(c(-1.5, -0.4, 0.55, 1.5)), limits = c(-1.5, 1.5),
    breaks = c(-1.5, 0, 1.5), labels = c("-1.5", "0", "1.5"),
    name = "Mean expression\n(scaled)"
  ) +
  scale_size(
    range = c(0.42, 2.45), limits = c(0, 50), breaks = c(10, 20, 30, 40),
    name = "Fraction of cells\nin group (%)"
  ) +
  scale_x_continuous(
    breaks = seq_along(focused_genes), labels = focused_genes,
    limits = c(0.45, length(focused_genes) + 0.55), expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = c(0.86, 0.915), labels = c("Treg", "Th17-like"),
    limits = c(0.82, 1.05), expand = c(0, 0)
  ) +
  guides(
    size = guide_legend(
      order = 1, title.position = "top", direction = "horizontal", nrow = 1,
      override.aes = list(fill = "#E77762")
    ),
    fill = guide_colorbar(
      order = 2, title.position = "top", direction = "horizontal",
      barwidth = unit(12, "mm"), barheight = unit(1.5, "mm"), ticks = FALSE
    )
  ) +
  coord_cartesian(clip = "off") +
  labs(tag = "D", x = NULL, y = NULL) +
  theme_classic(base_family = base_family, base_size = 7) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 90, hjust = 1, vjust = 0.5, size = 5.0,
      face = "italic", color = "#111111", margin = margin(t = 1.0)
    ),
    axis.text.y = element_text(size = 5.2, face = "plain", color = "#111111"),
    plot.tag = element_text(size = 8.0, face = "bold", color = "black"),
    plot.tag.position = c(0.005, 0.995),
    legend.position = "right", legend.box = "vertical",
    legend.spacing.y = unit(0.5, "mm"),
    legend.margin = margin(0, 0, 0, 0.8),
    legend.title = element_text(size = 5.0, lineheight = 0.92),
    legend.text = element_text(size = 5.0),
    legend.key.height = unit(2.4, "mm"),
    legend.key.width = unit(3.2, "mm"),
    plot.margin = margin(4.2, 1.5, 1.5, 4.5)
  )

# Complete CD4-state marker panel for Supplementary Figure 2.
supp_genes <- c(
  "CCR7", "TCF7", "LEF1", "SELL", "IL7R", "LTB", "MAL",
  "CCR6", "KLRB1", "RORA", "RORC", "CCL20",
  "FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2",
  "CD69", "FOS", "JUN", "ICOS", "HLA-DRA", "CD8A", "TRDC", "NKG7"
)
supp_modules <- c(
  setNames(rep("Naive", 4), supp_genes[1:4]),
  setNames(rep("Memory", 3), supp_genes[5:7]),
  setNames(rep("Th17-like", 5), supp_genes[8:12]),
  setNames(rep("Treg", 5), supp_genes[13:17]),
  setNames(rep("Activated", 5), supp_genes[18:22]),
  setNames(rep("Exclusion", 3), supp_genes[23:25])
)
dot_s <- make_dot_data(cd4, state_order, supp_genes, supp_modules)
dot_s[, state := factor(state, levels = rev(state_order))]
dot_s[, marker := factor(marker, levels = supp_genes)]

p_supp <- ggplot(dot_s, aes(marker, state)) +
  geom_point(shape = 21, size = 1.9, fill = "#F1F2F2", color = "#D3D7D8", stroke = 0.2) +
  geom_point(aes(size = pct_exp, fill = avg_exp_scaled), shape = 21, color = "#34393B", stroke = 0.22) +
  scale_fill_gradientn(
    colors = c("#F4F2F7", "#B9DDF3", "#4C9CCF", "#283E8D"),
    values = rescale(c(-1.5, -0.5, 0.5, 1.5)), limits = c(-1.5, 1.5),
    name = "Mean expression\n(row-scaled)"
  ) +
  scale_size(range = c(0.5, 3.9), limits = c(0, 100), breaks = c(10, 30, 50), name = "% expressed") +
  labs(title = "CD4 T-cell state marker evidence", x = NULL, y = NULL) +
  theme_minimal(base_family = base_family, base_size = 7) +
  theme(
    panel.grid.major = element_line(color = "#E5E8E9", linewidth = 0.22), panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, color = "#34393B", linewidth = 0.35),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6.0, face = "italic", color = "#202426"),
    axis.text.y = element_text(size = 6.2, face = "bold", color = "#202426"),
    plot.title = element_text(size = 7.3, hjust = 0.5, margin = margin(b = 3)),
    legend.position = "right", legend.title = element_text(size = 5.8), legend.text = element_text(size = 5.5),
    plot.margin = margin(5, 5, 4, 5)
  )

# Donor is the inferential unit. Ratios require no pseudocount because every
# donor has at least one retained Treg cell; this is checked below.
donor_counts <- cd4[, .N, by = .(library_id, health, onset, final_state)]
donor_wide <- dcast(donor_counts, library_id + health + onset ~ final_state, value.var = "N", fill = 0)
donor_wide[, total_CD4 := rowSums(.SD), .SDcols = state_order]
if (any(donor_wide[["Treg"]] == 0)) stop("At least one donor has zero Treg cells; ratio rule must be revised")
donor_wide[, `:=`(
  Th17_pct = 100 * `Th17-like` / total_CD4,
  Treg_pct = 100 * Treg / total_CD4,
  Th17_Treg_ratio = `Th17-like` / Treg
)]

metric_labels <- c(
  "Th17_pct" = "Th17-like / CD4 T (%)",
  "Treg_pct" = "Treg / CD4 T (%)",
  "Th17_Treg_ratio" = "Th17-like / Treg ratio"
)
donor_long <- melt(
  donor_wide,
  id.vars = c("library_id", "health", "onset", state_order, "total_CD4"),
  measure.vars = names(metric_labels), variable.name = "metric", value.name = "value"
)
donor_long[, condition := factor(fifelse(health == "Healthy", "Healthy", "LoSc"), levels = c("Healthy", "LoSc"))]
donor_long[, metric := factor(metric, levels = names(metric_labels), labels = unname(metric_labels))]

cliffs_delta <- function(x, y) mean(outer(x, y, FUN = function(a, b) sign(a - b)))
stats <- rbindlist(lapply(names(metric_labels), function(metric_name) {
  h <- donor_long[as.character(metric) == metric_labels[[metric_name]] & condition == "Healthy", value]
  l <- donor_long[as.character(metric) == metric_labels[[metric_name]] & condition == "LoSc", value]
  wt <- wilcox.test(l, h, exact = FALSE, conf.int = TRUE)
  data.table(
    metric = metric_labels[[metric_name]],
    healthy_n = length(h), losc_n = length(l),
    healthy_median = median(h), losc_median = median(l),
    median_difference_LoSc_minus_Healthy = median(l) - median(h),
    cliffs_delta_LoSc_vs_Healthy = cliffs_delta(l, h),
    wilcoxon_W = unname(wt$statistic), p_value = wt$p.value,
    method = "two-sided Wilcoxon rank-sum; donor-level"
  )
}))
stats[, metric := factor(metric, levels = unname(metric_labels))]
stats[, annotation := sprintf("P = %s; \u0394med = %.2f", format.pval(p_value, digits = 2, eps = 1e-4), median_difference_LoSc_minus_Healthy)]

condition_palette <- c("Healthy" = "#3FAF63", "LoSc" = "#D85A5A")
setorder(donor_long, metric, condition, library_id)
donor_long[, point_x := as.numeric(condition) + seq(-0.105, 0.105, length.out = .N), by = .(metric, condition)]

p_e <- ggplot(donor_long, aes(as.numeric(condition), value, color = condition, group = condition)) +
  geom_boxplot(
    width = 0.48, outlier.shape = NA, fill = "white", linewidth = 0.42,
    color = "#24282A", fatten = 1.3
  ) +
  geom_point(aes(x = point_x, shape = onset), size = 1.65, alpha = 0.92, stroke = 0.32) +
  geom_text(
    data = stats, aes(x = 1.5, y = Inf, label = annotation), inherit.aes = FALSE,
    vjust = 1.35, family = base_family, size = 2.05, color = "#303538"
  ) +
  facet_wrap(~metric, scales = "free_y", nrow = 1) +
  scale_color_manual(values = condition_palette, guide = "none") +
  scale_shape_manual(values = c("Adult" = 16, "Peds" = 21), name = NULL) +
  scale_x_continuous(breaks = 1:2, labels = c("Healthy", "LoSc"), limits = c(0.55, 2.45)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.24))) +
  labs(x = NULL, y = NULL) +
  theme_classic(base_family = base_family, base_size = 7) +
  theme(
    axis.line = element_line(linewidth = 0.38, color = "#282D2F"),
    axis.ticks = element_line(linewidth = 0.34, color = "#282D2F"),
    axis.text.x = element_text(size = 6.2, face = "bold", color = "#202426", margin = margin(t = 2)),
    axis.text.y = element_text(size = 5.8, color = "#202426"),
    strip.background = element_blank(), strip.text = element_text(size = 6.5, face = "bold", margin = margin(b = 2.5)),
    panel.spacing = unit(5.5, "mm"), legend.position = "top",
    legend.text = element_text(size = 5.8), legend.key.width = unit(4, "mm"),
    plot.margin = margin(6, 4.5, 4, 7)
  ) +
  geom_text(
    data = data.table(
      metric = factor(metric_labels[[1]], levels = unname(metric_labels)),
      x = -Inf, y = Inf, label = "E"
    ),
    aes(x = x, y = y, label = label), inherit.aes = FALSE,
    hjust = -0.85, vjust = 1.10, family = base_family,
    fontface = "bold", size = 3.2, color = "black"
  )

save_plot <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot)
  dev.off()
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    quartz(type = "pdf", file = paste0(stem, ".pdf"), width = width_in, height = height_in,
           family = base_family, bg = "white")
  } else {
    cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
              family = base_family, bg = "white")
  }
  print(plot)
  dev.off()
  agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in,
           units = "in", res = 600, background = "white", compression = "lzw")
  print(plot)
  dev.off()
  agg_png(paste0(stem, "_preview.png"), width = width_in, height = height_in,
          units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
}

save_plot(p_c, file.path(output_dir, "Figure2C_CD4_state_UMAP"), 94, 84)
save_plot(p_d, file.path(output_dir, "Figure2D_Th17_Treg_marker_dotplot"), 110, 25)
save_plot(p_e, file.path(output_dir, "Figure2E_donor_abundance"), 178, 68)
save_plot(p_supp, file.path(output_dir, "Supplementary_Figure2_CD4_marker_dotplot"), 178, 82)

fwrite(
  cd4[, .(cell_id, source_row, library_id, health, onset, tcell_res_1, final_state, CD4_UMAP_1, CD4_UMAP_2)],
  file.path(output_dir, "Figure2C_CD4_cell_metadata.csv.gz")
)
fwrite(cluster_annotation, file.path(output_dir, "Figure2C_CD4_cluster_annotation.csv"))
fwrite(state_summary, file.path(output_dir, "Figure2C_CD4_state_summary.csv"))
fwrite(data.table(state = state_order, hex = unname(state_palette)), file.path(output_dir, "Figure2C_CD4_state_palette.csv"))
fwrite(
  excluded[, .(cell_id, source_row, library_id, health, onset, tcell_res_1, final_state)],
  file.path(output_dir, "Figure2C_excluded_Tcell_metadata.csv.gz")
)
fwrite(dot_d, file.path(output_dir, "Figure2D_marker_dotplot_source_data.csv"))
fwrite(data.table(marker = focused_genes, module = unname(focused_modules[focused_genes]), available = focused_genes %chin% names(cd4)),
       file.path(output_dir, "Figure2D_marker_panel.csv"))
fwrite(marker_audit, file.path(output_dir, "Figure2D_marker_interpretation_audit.csv"))
fwrite(
  data.table(
    marker = c("IL23R", "IL17A", "IL17F"),
    available_in_processed_matrix = c("IL23R", "IL17A", "IL17F") %chin% names(cd4),
    figure_action = "omitted from Figure 2D because unavailable in the processed matrix; no imputation"
  ),
  file.path(output_dir, "Figure2D_missing_classical_Th17_marker_audit.csv")
)
fwrite(donor_long, file.path(output_dir, "Figure2E_donor_abundance_source_data.csv"))
fwrite(stats, file.path(output_dir, "Figure2E_donor_abundance_statistics.csv"))
fwrite(dot_s, file.path(output_dir, "Supplementary_Figure2_CD4_marker_dotplot_source_data.csv"))
fwrite(data.table(marker = supp_genes, module = unname(supp_modules[supp_genes])),
       file.path(output_dir, "Supplementary_Figure2_CD4_marker_panel.csv"))

qa <- c(
  "panels=Figure 2C, Figure 2D, Figure 2E, Supplementary CD4 marker dot plot",
  "backend=R only",
  sprintf("input_T_cells=%s", nrow(dt)),
  sprintf("retained_CD4_cells=%s", nrow(cd4)),
  sprintf("excluded_CD8_or_gamma_delta_like_cells=%s", nrow(excluded)),
  sprintf("retained_donors=%s", uniqueN(cd4$library_id)),
  "state_definition=resolution-1 T-cell communities annotated by multi-gene modules and donor support",
  "cycling_decision=no donor-supported MKI67/TOP2A coexpressing community; no Cycling state forced",
  "CD4_visualization=retained CD4 cells re-embedded from author Harmony dimensions 1:30; labels not changed by UMAP",
  "Figure2D_encoding=dot size percent expressed; color mean normalized expression scaled per gene across all five retained CD4 states before focusing on the two displayed states",
  "Figure2D_display_label=Th17-like only; ambiguous memory-state wording removed",
  "Figure2D_layout=single framed compact horizontal two-row dot plot with near-touching rows, marker-family labels and color bars above the matrix border, and right-side legends; 110 x 25 mm",
  "Figure2D_identity_interpretation=Th17-like denotes a transcriptomic phenotype supported by multi-gene Th17-associated evidence; it is not equivalent to a proven IL-17-secreting effector population",
  "Figure2D_lineage_interpretation=TRDC,CD8A,CD8B,NKG7 are mixture/context indicators; no retained cell is excluded because of positivity for any single marker",
  paste0("processed_matrix_missing_markers=", paste(omitted_unavailable_th17, collapse = ",")),
  "missing_marker_rendering=IL23R,IL17A,IL17F omitted from Figure 2D per revision request and recorded in a separate audit table; values not imputed",
  "Figure2E_inferential_unit=donor",
  "Figure2E_denominator=all five retained CD4 states per donor",
  "Figure2E_ratio=raw Th17-like cell count divided by raw Treg cell count; no pseudocount required",
  "Figure2E_test=two-sided Wilcoxon rank-sum; all donor points displayed",
  "effect_size=LoSc-minus-Healthy median difference plus Cliff's delta in source statistics",
  "vector=SVG and PDF; raster=600 dpi TIFF; preview=300 dpi PNG",
  "integrity=no expression values or state proportions edited to enforce an expected direction"
)
writeLines(qa, file.path(output_dir, "Figure2CDE_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2CDE_sessionInfo.txt"))

message("Figure 2C-E complete: ", nrow(cd4), " CD4 cells from ", uniqueN(cd4$library_id), " donors")
