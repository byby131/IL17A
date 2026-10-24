#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(RANN)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/map_and_plot_figure2F.R <stromal_reference_query_expression.csv.gz> <all.H5Seurat> <output_dir>")
}

input_path <- normalizePath(args[[1]], mustWork = TRUE)
h5_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
base_family <- "Helvetica"
set.seed(264508)

state_order <- c(
  "Homeostatic/adipogenic ASPC",
  "Papillary/conventional fibroblast",
  "ECM-remodeling fibroblast",
  "Pericyte",
  "VSMC"
)
state_palette <- c(
  # Darkened derivatives of the supplied reference hues: lower luminance and
  # restrained chroma improve legibility without changing category identity.
  "Homeostatic/adipogenic ASPC" = "#1F6F65",
  "Papillary/conventional fibroblast" = "#644090",
  "ECM-remodeling fibroblast" = "#B85C38",
  "Pericyte" = "#357F9F",
  "VSMC" = "#184A73",
  "Uncertain" = "#D98A2B",
  "Disease-deviated" = "#A64D79"
)
mapped_query_color <- "#9B9B9B"

# Resolution 0.2 is intentionally merged to five biologically interpretable
# Healthy states. The merge is fixed before LoSc mapping.
healthy_cluster_map <- data.table(
  healthy_res_0.2 = as.character(0:6),
  reference_state = c(
    "Papillary/conventional fibroblast",
    "Homeostatic/adipogenic ASPC",
    "VSMC",
    "Pericyte",
    "ECM-remodeling fibroblast",
    "VSMC",
    "ECM-remodeling fibroblast"
  ),
  evidence = c(
    "COL14A1/WIF1/APCDD1/SLPI with DCN/LUM/COL1A1",
    "PDGFRA/CD34/C7/CFD/APOD/CXCL12 with PPARG/CEBPA",
    "ACTA2/TAGLN/MYH11/MCAM contractile mural state",
    "RGS5/PDGFRB with lower MYH11/CNN1 than VSMC",
    "POSTN/FN1/THBS2 extracellular-matrix program",
    "ACTA2/TAGLN/MYH11/CNN1 high contractile mural state",
    "SFRP4-high extracellular-matrix fibroblast state"
  )
)

st <- fread(input_path)
st[, healthy_res_0.2 := as.character(healthy_res_0.2)]
ref <- st[health == "Healthy"]
query <- st[health == "LS"]
ref <- merge(ref, healthy_cluster_map, by = "healthy_res_0.2", all.x = TRUE, sort = FALSE)
if (anyNA(ref$reference_state)) stop("Healthy reference clusters remain unannotated")
ref[, reference_state := factor(reference_state, levels = state_order)]

h5 <- H5File$new(h5_path, "r")
on.exit(h5$close_all())
harmony_all <- h5[["reductions"]][["harmony"]][["cell.embeddings"]]$read()
ref_latent <- harmony_all[ref$source_row, 1:30, drop = FALSE]
query_latent <- harmony_all[query$source_row, 1:30, drop = FALSE]

l2_normalize <- function(x) {
  norm <- sqrt(rowSums(x^2))
  norm[norm == 0] <- 1
  x / norm
}
ref_norm <- l2_normalize(ref_latent)
query_norm <- l2_normalize(query_latent)
k_map <- 30L

summarize_neighbors <- function(nn, labels, levels) {
  neighbor_labels <- matrix(labels[nn$nn.idx], nrow = nrow(nn$nn.idx), ncol = ncol(nn$nn.idx))
  votes <- t(apply(neighbor_labels, 1, function(x) tabulate(match(x, levels), nbins = length(levels))))
  probs <- votes / rowSums(votes)
  prediction <- levels[max.col(probs, ties.method = "first")]
  max_probability <- apply(probs, 1, max)
  entropy <- -rowSums(ifelse(probs > 0, probs * log(probs), 0)) / log(length(levels))
  list(
    prediction = prediction,
    max_probability = max_probability,
    entropy = entropy,
    mean_distance = rowMeans(nn$nn.dists),
    probabilities = probs
  )
}

# Leave-one-donor-out Healthy calibration estimates in-reference uncertainty
# without allowing a donor to map to itself.
calibration <- rbindlist(lapply(sort(unique(ref$library_id)), function(donor) {
  train <- ref$library_id != donor
  test <- !train
  nn <- RANN::nn2(ref_norm[train, , drop = FALSE], ref_norm[test, , drop = FALSE], k = k_map)
  s <- summarize_neighbors(nn, as.character(ref$reference_state[train]), state_order)
  data.table(
    cell_id = ref$cell_id[test], library_id = ref$library_id[test],
    truth = as.character(ref$reference_state[test]), prediction = s$prediction,
    max_probability = s$max_probability, entropy = s$entropy,
    mean_distance = s$mean_distance
  )
}))
calibration[, correct := truth == prediction]

# Prespecified 95% in-reference envelope. Out-of-manifold distance has priority;
# in-manifold low-confidence cells are labelled uncertain.
confidence_cutoff <- unname(quantile(calibration$max_probability, 0.05, na.rm = TRUE))
entropy_cutoff <- unname(quantile(calibration$entropy, 0.95, na.rm = TRUE))
distance_cutoff <- unname(quantile(calibration$mean_distance, 0.95, na.rm = TRUE))

nn_query <- RANN::nn2(ref_norm, query_norm, k = k_map)
query_map <- summarize_neighbors(nn_query, as.character(ref$reference_state), state_order)
query[, `:=`(
  predicted_reference_state = query_map$prediction,
  mapping_confidence = query_map$max_probability,
  mapping_entropy = query_map$entropy,
  mean_reference_distance = query_map$mean_distance
)]
for (j in seq_along(state_order)) {
  query[[paste0("prob_", make.names(state_order[[j]]))]] <- query_map$probabilities[, j]
}
query[, mapping_status := fifelse(
  mean_reference_distance > distance_cutoff,
  "Disease-deviated",
  fifelse(mapping_confidence < confidence_cutoff | mapping_entropy > entropy_cutoff, "Uncertain", "Mapped")
)]
query[, mapping_label := fifelse(mapping_status == "Mapped", predicted_reference_state, mapping_status)]

# Recluster all non-confident query cells as a separate disease-deviation
# compartment; this label is retained for downstream 2H profiling.
nonmapped <- query$mapping_status != "Mapped"
query[, disease_deviation_subcluster := NA_character_]
if (sum(nonmapped) >= 50L) {
  k_dev <- min(30L, max(10L, floor(sqrt(sum(nonmapped)))))
  dev_dummy <- sparseMatrix(
    i = rep.int(1L, sum(nonmapped)), j = seq_len(sum(nonmapped)), x = rep.int(1, sum(nonmapped)),
    dims = c(1L, sum(nonmapped)), dimnames = list("placeholder", query$cell_id[nonmapped])
  )
  dev_obj <- CreateSeuratObject(counts = dev_dummy, min.cells = 0, min.features = as.integer(0))
  dev_embed <- query_latent[nonmapped, , drop = FALSE]
  rownames(dev_embed) <- query$cell_id[nonmapped]
  colnames(dev_embed) <- paste0("HARMONY_", seq_len(ncol(dev_embed)))
  dev_obj[["harmony"]] <- CreateDimReducObject(embeddings = dev_embed, key = "HARMONY_", assay = "RNA")
  dev_obj <- FindNeighbors(dev_obj, reduction = "harmony", dims = 1:30, k.param = k_dev,
                           graph.name = c("dev_nn", "dev_snn"), verbose = FALSE)
  dev_obj <- FindClusters(dev_obj, graph.name = "dev_snn", resolution = 0.4,
                          algorithm = 1, random.seed = 264508,
                          cluster.name = "disease_deviation_subcluster", verbose = FALSE)
  dev_cluster <- as.character(dev_obj$disease_deviation_subcluster)
  names(dev_cluster) <- colnames(dev_obj)
  query[nonmapped, disease_deviation_subcluster := dev_cluster[cell_id]]
}

# Reference-query neighborhood mixing is a descriptive QC only.
all_norm <- rbind(ref_norm, query_norm)
all_origin <- c(rep("Healthy reference", nrow(ref)), rep("LoSc query", nrow(query)))
mix_nn <- RANN::nn2(all_norm, all_norm, k = 31L)
mix_index <- mix_nn$nn.idx[, -1, drop = FALSE]
mix_origin <- matrix(all_origin[mix_index], nrow = nrow(mix_index), ncol = ncol(mix_index))
mix_fraction <- rowMeans(mix_origin != all_origin)
ref[, opposite_origin_neighbor_fraction := mix_fraction[seq_len(nrow(ref))]]
query[, opposite_origin_neighbor_fraction := mix_fraction[nrow(ref) + seq_len(nrow(query))]]

# Threshold sensitivity is evaluated over fixed Healthy quantile envelopes.
sensitivity_quantiles <- c(0.90, 0.95, 0.975)
sensitivity <- rbindlist(lapply(sensitivity_quantiles, function(q) {
  conf_q <- unname(quantile(calibration$max_probability, 1 - q, na.rm = TRUE))
  ent_q <- unname(quantile(calibration$entropy, q, na.rm = TRUE))
  dist_q <- unname(quantile(calibration$mean_distance, q, na.rm = TRUE))
  status <- fifelse(
    query$mean_reference_distance > dist_q, "Disease-deviated",
    fifelse(query$mapping_confidence < conf_q | query$mapping_entropy > ent_q, "Uncertain", "Mapped")
  )
  data.table(
    envelope = paste0(q * 100, "%"), quantile = q,
    confidence_cutoff = conf_q, entropy_cutoff = ent_q, distance_cutoff = dist_q,
    status = c("Mapped", "Uncertain", "Disease-deviated"),
    n_cells = c(sum(status == "Mapped"), sum(status == "Uncertain"), sum(status == "Disease-deviated"))
  )[, fraction_percent := 100 * n_cells / sum(n_cells)]
}))

reference_counts <- ref[, .(n_cells = .N, n_donors = uniqueN(library_id)), by = reference_state]
query_counts <- query[, .(n_cells = .N, n_donors = uniqueN(library_id)), by = .(mapping_status, mapping_label)]
donor_contribution <- query[, .N, by = .(library_id, onset, mapping_label)]
donor_contribution[, fraction_percent := 100 * N / sum(N), by = library_id]

label_positions <- ref[, {
  x0 <- median(reference_UMAP_1)
  y0 <- median(reference_UMAP_2)
  d2 <- (reference_UMAP_1 - x0)^2 + (reference_UMAP_2 - y0)^2
  .(reference_UMAP_1 = reference_UMAP_1[which.min(d2)], reference_UMAP_2 = reference_UMAP_2[which.min(d2)])
}, by = reference_state]
label_positions[, label := c("Homeostatic/\nadipogenic ASPC", "Papillary/\nconventional FB",
                             "ECM-remodeling FB", "Pericyte", "VSMC")[
                               match(as.character(reference_state), state_order)]]

mapping_levels <- c(state_order, "Uncertain", "Disease-deviated")
query[, mapping_label := factor(mapping_label, levels = mapping_levels)]

status_labels <- query[, .(n_cells = .N), by = mapping_status]
status_labels[, fraction_percent := 100 * n_cells / nrow(query)]
status_labels <- status_labels[mapping_status %chin% c("Uncertain", "Disease-deviated")]
status_labels[, label := sprintf("%s\nn = %s (%.1f%%)", mapping_status, n_cells, fraction_percent)]

x_limits <- range(c(ref$reference_UMAP_1, query$reference_UMAP_1))
y_limits <- range(c(ref$reference_UMAP_2, query$reference_UMAP_2))
x_span <- diff(x_limits)
y_span <- diff(y_limits)
# Place the long ADSC label just above the dense core so 12-pt text does not
# obscure the cells while remaining visually tied to the population.
label_positions[reference_state == "Homeostatic/adipogenic ASPC", `:=`(
  label = "Homeostatic/\nadipogenic ASPC",
  reference_UMAP_1 = median(ref[reference_state == "Homeostatic/adipogenic ASPC", reference_UMAP_1]),
  reference_UMAP_2 = 1.25
)]
label_positions[reference_state == "Papillary/conventional fibroblast", `:=`(
  reference_UMAP_1 = 6.25,
  reference_UMAP_2 = 5.80
)]
label_positions[reference_state == "ECM-remodeling fibroblast", `:=`(
  reference_UMAP_1 = -0.45,
  reference_UMAP_2 = 4.35
)]
# Reserve a lower-left gutter so the compact UMAP axes remain outside the
# cellular point cloud while staying aligned across both panels.
axis_x0 <- x_limits[1] - 0.125 * x_span
axis_y0 <- y_limits[1] - 0.115 * y_span
panel_x_limits <- c(x_limits[1] - 0.175 * x_span, x_limits[2] + 0.020 * x_span)
panel_y_limits <- c(y_limits[1] - 0.175 * y_span, y_limits[2] + 0.020 * y_span)

umap_axes <- list(
  annotate("segment", x = axis_x0, xend = axis_x0 + 0.090 * x_span,
           y = axis_y0, yend = axis_y0, linewidth = 0.30, color = "#252A2C",
           arrow = arrow(length = unit(1.0, "mm"), type = "closed")),
  annotate("segment", x = axis_x0, xend = axis_x0,
           y = axis_y0, yend = axis_y0 + 0.090 * y_span, linewidth = 0.30,
           color = "#252A2C", arrow = arrow(length = unit(1.0, "mm"), type = "closed")),
  annotate("text", x = axis_x0 + 0.048 * x_span, y = axis_y0 - 0.038 * y_span,
           label = "UMAP1", family = base_family, size = 2.05, color = "#111111"),
  annotate("text", x = axis_x0 - 0.038 * x_span, y = axis_y0 + 0.048 * y_span,
           label = "UMAP2", angle = 90, family = base_family, size = 2.05, color = "#111111")
)

p_ref <- ggplot(ref, aes(reference_UMAP_1, reference_UMAP_2, color = reference_state)) +
  geom_point(size = 0.27, alpha = 1, stroke = 0) +
  geom_text(
    data = label_positions, aes(label = label), color = "#111111",
    family = "Arial", fontface = "bold", size = 8.5 / ggplot2::.pt,
    lineheight = 0.92, show.legend = FALSE
  ) +
  umap_axes +
  scale_color_manual(values = state_palette[mapping_levels], limits = mapping_levels,
                     drop = FALSE, name = NULL) +
  coord_equal(xlim = panel_x_limits, ylim = panel_y_limits,
              expand = FALSE, clip = "off") +
  labs(title = "Healthy stromal reference") +
  theme_void(base_family = base_family) +
  theme(plot.title = element_text(size = 7.6, color = "#111111", hjust = 0.5, margin = margin(b = 2)),
        legend.position = "none", plot.margin = margin(5, 3, 2, 4))

p_query_base <- ggplot() +
  geom_point(data = ref, aes(reference_UMAP_1, reference_UMAP_2),
             size = 0.12, color = "#E3E3E3", alpha = 0.52, stroke = 0) +
  geom_point(
    data = query[mapping_status == "Mapped"],
    aes(reference_UMAP_1, reference_UMAP_2),
    size = 0.24, color = mapped_query_color, alpha = 0.76, stroke = 0
  ) +
  geom_point(
    data = query[mapping_status == "Uncertain"],
    aes(reference_UMAP_1, reference_UMAP_2),
    shape = 21, size = 0.47, fill = state_palette[["Uncertain"]],
    color = "#87520C", stroke = 0.13, alpha = 1
  ) +
  geom_point(
    data = query[mapping_status == "Disease-deviated"],
    aes(reference_UMAP_1, reference_UMAP_2),
    shape = 21, size = 0.52, fill = state_palette[["Disease-deviated"]],
    color = "#672D4B", stroke = 0.13, alpha = 1
  ) +
  umap_axes +
  coord_equal(xlim = panel_x_limits, ylim = panel_y_limits,
              expand = FALSE, clip = "off") +
  labs(title = "LoSc query projected onto Healthy reference") +
  theme_void(base_family = base_family) +
  theme(plot.title = element_text(size = 7.6, color = "#111111", hjust = 0.5, margin = margin(b = 2)),
        legend.position = "none", plot.margin = margin(5, 4, 2, 3))

p_query_boxed <- p_query_base +
  annotate(
    "label",
    x = x_limits[1] + 0.56 * x_span, y = y_limits[2] - 0.035 * y_span,
    label = status_labels[mapping_status == "Disease-deviated", label],
    family = "Arial", fontface = "bold", size = 2.00, lineheight = 0.93,
    fill = state_palette[["Disease-deviated"]], color = "white",
    label.padding = unit(0.10, "lines"), label.r = unit(0.05, "lines"),
    linewidth = 0.18
  ) +
  annotate(
    "label",
    x = x_limits[1] + 0.83 * x_span, y = y_limits[2] - 0.035 * y_span,
    label = status_labels[mapping_status == "Uncertain", label],
    family = "Arial", fontface = "bold", size = 2.00, lineheight = 0.93,
    fill = state_palette[["Uncertain"]], color = "#2B230F",
    label.padding = unit(0.10, "lines"), label.r = unit(0.05, "lines"),
    linewidth = 0.18
  )

p_query_status_text <- p_query_base +
  annotate(
    "text",
    x = x_limits[1] + 0.56 * x_span, y = y_limits[2] - 0.035 * y_span,
    label = status_labels[mapping_status == "Disease-deviated", label],
    family = "Arial", fontface = "bold", size = 2.35, lineheight = 0.93,
    color = state_palette[["Disease-deviated"]]
  ) +
  annotate(
    "text",
    x = x_limits[1] + 0.83 * x_span, y = y_limits[2] - 0.035 * y_span,
    label = status_labels[mapping_status == "Uncertain", label],
    family = "Arial", fontface = "bold", size = 2.35, lineheight = 0.93,
    color = state_palette[["Uncertain"]]
  )

p_f <- (p_ref | p_query_boxed) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = list("F"), theme = theme(
    plot.tag = element_text(family = base_family, face = "bold", size = 9),
    plot.tag.position = c(0.003, 0.995), plot.margin = margin(1, 1, 1, 1)
  ))

p_f_status_text <- (p_ref | p_query_status_text) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = list("F"), theme = theme(
    plot.tag = element_text(family = base_family, face = "bold", size = 9),
    plot.tag.position = c(0.003, 0.995), plot.margin = margin(1, 1, 1, 1)
  ))

# Supplementary mapping QC.
qc_ref <- calibration[, .(
  origin = "Healthy LODO", max_probability, entropy, mean_distance,
  opposite_origin_neighbor_fraction = ref$opposite_origin_neighbor_fraction[match(cell_id, ref$cell_id)]
)]
qc_query <- query[, .(
  origin = "LoSc query", max_probability = mapping_confidence,
  entropy = mapping_entropy, mean_distance = mean_reference_distance,
  opposite_origin_neighbor_fraction
)]
qc_long <- melt(
  rbindlist(list(qc_ref, qc_query)), id.vars = "origin",
  measure.vars = c("max_probability", "entropy", "mean_distance", "opposite_origin_neighbor_fraction"),
  variable.name = "metric", value.name = "value"
)
metric_names <- c(
  max_probability = "Mapping confidence", entropy = "Normalized entropy",
  mean_distance = "Mean reference distance", opposite_origin_neighbor_fraction = "Opposite-origin neighbor fraction"
)
qc_long[, metric := factor(metric, levels = names(metric_names), labels = unname(metric_names))]
qc_origin_palette <- c("Healthy LODO" = "#3FAF63", "LoSc query" = "#D85A5A")
p_qc <- ggplot(qc_long, aes(value, color = origin, fill = origin)) +
  geom_density(alpha = 0.13, linewidth = 0.48, adjust = 1.1) +
  facet_wrap(~metric, scales = "free", nrow = 1) +
  scale_color_manual(values = qc_origin_palette) +
  scale_fill_manual(values = qc_origin_palette) +
  labs(x = NULL, y = "Density", color = NULL, fill = NULL,
       title = "Healthy-calibrated mapping quality") +
  theme_classic(base_family = base_family, base_size = 7) +
  theme(strip.background = element_blank(), strip.text = element_text(size = 6.2, face = "bold"),
        axis.text = element_text(size = 5.7), axis.title.y = element_text(size = 6),
        legend.position = "top", legend.text = element_text(size = 5.8),
        plot.title = element_text(size = 7.2, hjust = 0.5), panel.spacing = unit(4, "mm"))

donor_contribution[, mapping_label := factor(mapping_label, levels = mapping_levels)]
donor_order <- donor_contribution[, .(uncertain_fraction = sum(fraction_percent[mapping_label %in% c("Uncertain", "Disease-deviated")])), by = library_id][order(uncertain_fraction), library_id]
donor_contribution[, library_id := factor(library_id, levels = donor_order)]
p_donor <- ggplot(donor_contribution, aes(library_id, fraction_percent, fill = mapping_label)) +
  geom_col(width = 0.78) +
  scale_fill_manual(values = state_palette[mapping_levels], drop = FALSE, name = NULL) +
  labs(x = "LoSc donor", y = "Mapped composition (%)", title = "Query-state contribution by donor") +
  theme_classic(base_family = base_family, base_size = 7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 5.4),
        axis.text.y = element_text(size = 5.6), axis.title = element_text(size = 6.2),
        legend.position = "right", legend.text = element_text(size = 5.2),
        plot.title = element_text(size = 7.2, hjust = 0.5))

sensitivity[, status := factor(status, levels = c("Mapped", "Uncertain", "Disease-deviated"))]
p_sensitivity <- ggplot(sensitivity, aes(quantile * 100, fraction_percent, color = status, group = status)) +
  geom_line(linewidth = 0.65) + geom_point(size = 1.7) +
  scale_color_manual(values = c(
    "Mapped" = mapped_query_color,
    "Uncertain" = state_palette[["Uncertain"]],
    "Disease-deviated" = state_palette[["Disease-deviated"]]
  ), name = NULL) +
  scale_x_continuous(breaks = sensitivity_quantiles * 100) +
  labs(x = "Healthy calibration envelope (%)", y = "LoSc query cells (%)",
       title = "Mapping-threshold sensitivity") +
  theme_classic(base_family = base_family, base_size = 7) +
  theme(axis.text = element_text(size = 5.8), axis.title = element_text(size = 6.2),
        legend.position = "top", legend.text = element_text(size = 5.5),
        plot.title = element_text(size = 7.2, hjust = 0.5))

save_plot <- function(plot, stem, width_mm, height_mm) {
  width_in <- width_mm / 25.4
  height_in <- height_mm / 25.4
  svglite(paste0(stem, ".svg"), width = width_in, height = height_in, bg = "white")
  print(plot); dev.off()
  if (identical(Sys.info()[["sysname"]], "Darwin")) {
    quartz(type = "pdf", file = paste0(stem, ".pdf"), width = width_in, height = height_in,
           family = base_family, bg = "white")
  } else {
    cairo_pdf(paste0(stem, ".pdf"), width = width_in, height = height_in,
              family = base_family, bg = "white")
  }
  print(plot); dev.off()
  agg_tiff(paste0(stem, ".tiff"), width = width_in, height = height_in, units = "in",
           res = 600, background = "white", compression = "lzw")
  print(plot); dev.off()
  agg_png(paste0(stem, "_preview.png"), width = width_in, height = height_in, units = "in",
          res = 300, background = "white")
  print(plot); dev.off()
}

save_plot(p_f, file.path(output_dir, "Figure2F_healthy_reference_LoSc_mapping"), 183, 88)
save_plot(p_f_status_text, file.path(output_dir, "Figure2F_healthy_reference_LoSc_mapping_status_text_only"), 183, 88)
save_plot(p_qc, file.path(output_dir, "Supplementary_Figure2F_mapping_QC"), 183, 67)
save_plot(p_donor, file.path(output_dir, "Supplementary_Figure2F_donor_contribution"), 150, 80)
save_plot(p_sensitivity, file.path(output_dir, "Supplementary_Figure2F_threshold_sensitivity"), 92, 72)

fwrite(ref, file.path(output_dir, "Figure2F_healthy_reference_cell_metadata.csv.gz"))
fwrite(query, file.path(output_dir, "Figure2F_LoSc_query_mapping_metadata.csv.gz"))
fwrite(healthy_cluster_map, file.path(output_dir, "Figure2F_healthy_cluster_annotation.csv"))
fwrite(reference_counts, file.path(output_dir, "Figure2F_reference_state_counts.csv"))
fwrite(query_counts, file.path(output_dir, "Figure2F_query_mapping_counts.csv"))
fwrite(calibration, file.path(output_dir, "Supplementary_Figure2F_LODO_calibration.csv.gz"))
fwrite(qc_long, file.path(output_dir, "Supplementary_Figure2F_mapping_QC_source_data.csv.gz"))
fwrite(donor_contribution, file.path(output_dir, "Supplementary_Figure2F_donor_contribution_source_data.csv"))
fwrite(sensitivity, file.path(output_dir, "Supplementary_Figure2F_threshold_sensitivity_source_data.csv"))
fwrite(data.table(label = names(state_palette), hex = unname(state_palette)), file.path(output_dir, "Figure2F_palette.csv"))

thresholds <- c(
  sprintf("k_neighbors=%s", k_map),
  sprintf("confidence_cutoff_Healthy_5th_percentile=%.6f", confidence_cutoff),
  sprintf("entropy_cutoff_Healthy_95th_percentile=%.6f", entropy_cutoff),
  sprintf("distance_cutoff_Healthy_95th_percentile=%.6f", distance_cutoff),
  sprintf("Healthy_LODO_accuracy=%.6f", mean(calibration$correct)),
  sprintf("mapped_query_cells=%s", sum(query$mapping_status == "Mapped")),
  sprintf("uncertain_query_cells=%s", sum(query$mapping_status == "Uncertain")),
  sprintf("disease_deviated_query_cells=%s", sum(query$mapping_status == "Disease-deviated"))
)
writeLines(thresholds, file.path(output_dir, "Figure2F_mapping_thresholds.txt"))
writeLines(c(
  "panel=Figure 2F",
  "backend=R only",
  "core_conclusion=LoSc stromal cells are projected onto an equally scaled Healthy-only reference, with the 420 out-of-manifold disease-deviated cells highlighted as the visual focus and 243 uncertain cells shown separately.",
  sprintf("healthy_reference_cells=%s", nrow(ref)),
  sprintf("losc_query_cells=%s", nrow(query)),
  "reference_states=defined exclusively from Healthy resolution-0.2 clusters and multi-gene marker evidence",
  "mapping=30-nearest-neighbor voting in L2-normalized author Harmony dimensions 1:30",
  "calibration=leave-one-Healthy-donor-out; 95% in-reference envelope",
  "priority_rule=distance outside envelope -> Disease-deviated; otherwise low confidence/high entropy -> Uncertain",
  "reclustering=all non-mapped LoSc cells reclustered at resolution 0.4 in the same latent space",
  "query_visualization=out-of-sample transform into Healthy-only UMAP model",
  "panel_geometry=Healthy and LoSc panels use identical coordinate limits, equal aspect ratio, and equal layout widths; the redundant two-row bottom legend is removed",
  "status_emphasis=Disease-deviated and Uncertain cells are overplotted with larger deep-berry and amber symbols and direct n/percent labels",
  "status_annotation_versions=boxed color callouts and color-matched bold text without background are both exported",
  "query_color_rule=all ordinary mapped LoSc cells are grey; Uncertain and Disease-deviated use distinct amber and deep-berry colors absent from the five Healthy states",
  "label_style=Healthy state labels are Arial 8.5 pt bold black text without boxes or background fill; the ADSC label is manually positioned above its dense cell core",
  "axis_placement=compact UMAP1/UMAP2 arrows and labels sit in a dedicated lower-left gutter outside the cellular point cloud in both panels",
  "interpretation_guard=disease-deviated denotes cells outside the Healthy reference envelope; it is not presented as a donor-level differential-abundance significance test",
  "integrity=LoSc cells do not define Healthy states, UMAP geometry, or mapping thresholds",
  "vector=SVG and PDF; raster=600 dpi TIFF; preview=300 dpi PNG"
), file.path(output_dir, "Figure2F_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2F_sessionInfo.txt"))
message("Figure 2F complete: ", paste(thresholds[5:8], collapse = "; "))
