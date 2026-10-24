#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(CellChat)
  library(ggplot2)
  library(svglite)
  library(ragg)
})

source("scripts/figure2_ghi_common.R")

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else "outputs/CellChat_Q8_donor_significance"
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  h5 = "data/GSE264508/GSE264508_data.all.H5Seurat",
  cd4 = "outputs/Figure2C_CD4_cell_metadata.csv.gz",
  receiver = "outputs/Figure2H_Q8Q11_revision/audit/Figure2H_Q8_Q11_frozen_cell_metadata.csv.gz",
  pooled_lr = "outputs/CellChat_CD4_Q8Q11_vs_Healthy/source_data/CellChat_significant_CD4_to_ASPC_LR_interactions.csv.gz"
)
missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths)) stop("Missing inputs: ", paste(missing_paths, collapse = ", "))

set.seed(264508)
future::plan("sequential")
source_levels <- c("Naive CD4 T", "Memory CD4 T", "Th17-like", "Treg", "Activated CD4 T")
condition_levels <- c("Healthy", "Q8")
condition_palette <- c(Healthy = "#3569A8", Q8 = "#D85A4F")
min_source_cells <- 3L
min_receiver_cells <- 3L

# Figure contract:
# Core conclusion: test whether donor-level CD4-to-ASPC communication scores are
# higher in Q8 than in Healthy for each of five CD4 states.
# Archetype: one-panel quantitative comparison.
# Statistics: exact two-sided label-permutation test across biological donors;
# BH correction across the five CD4-state comparisons; mean +/- SEM shown.
# Output: 183 x 72 mm; editable SVG/PDF, 600-dpi TIFF, PNG preview and source CSV.

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.4, colour = "#222222"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#222222"),
      axis.text = element_text(colour = "#111111"),
      axis.title = element_text(colour = "#111111"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      legend.title = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = base_size + 1.2, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.1, colour = "#4D4D4D"),
      plot.caption = element_text(size = base_size - 1, colour = "#5A5A5A", hjust = 0),
      plot.margin = margin(4, 6, 3, 5)
    )
}

message("Reading frozen annotations")
cd4 <- fread(paths$cd4)[final_state %chin% source_levels, .(
  cell_id,
  source_row_h5 = as.integer(source_row),
  donor_id = as.character(library_id),
  condition = fifelse(health == "Healthy", "Healthy", "Q8"),
  cell_group = as.character(final_state),
  receiver_state = NA_character_
)]
receiver <- fread(paths$receiver)[
  qc_status == "Retained" & analysis_group %chin% c("Healthy ASPC-like", "LoSc Q8"),
  .(
    cell_id,
    source_row_h5 = as.integer(source_row_h5),
    donor_id = as.character(donor_id),
    condition = fifelse(analysis_group == "Healthy ASPC-like", "Healthy", "Q8"),
    cell_group = "ASPC receiver",
    receiver_state = as.character(analysis_group)
  )
]
meta <- unique(rbindlist(list(cd4, receiver), use.names = TRUE), by = "cell_id")
setorder(meta, source_row_h5)

message("Fixing a common ligand-receptor universe from the matched pooled CellChat analysis")
lr_pooled <- fread(paths$pooled_lr)[
  analysis == "cell-number matched" & model %chin% c("Healthy", "Q8") & source %chin% source_levels
]
lr_universe <- unique(lr_pooled[, .(source, interaction_name)])
db_interaction <- as.data.table(CellChatDB.human$interaction, keep.rownames = "db_row")
lr_universe <- merge(lr_universe, db_interaction, by = "interaction_name", all.x = TRUE)
if (anyNA(lr_universe$ligand) || anyNA(lr_universe$receptor)) stop("LR universe did not map to CellChatDB")
fwrite(lr_universe, file.path(source_dir, "CellChat_fixed_LR_universe_by_CD4_state.csv"))

complex_input <- CellChatDB.human$complex
needed_genes <- unique(c(lr_universe$ligand, lr_universe$receptor))
complex_names <- intersect(needed_genes, rownames(complex_input))
if (length(complex_names)) {
  complex_genes <- unlist(complex_input[complex_names, grep("^subunit", names(complex_input)), drop = FALSE],
                          use.names = FALSE)
  needed_genes <- unique(c(setdiff(needed_genes, complex_names), complex_genes[complex_genes != ""]))
}

message("Reading expression counts for donor-level scoring")
counts <- read_h5_counts(paths$h5, meta$source_row_h5, genes = needed_genes)
meta <- meta[match(colnames(counts), cell_id)]
stopifnot(identical(meta$cell_id, colnames(counts)))
libsize_all_genes <- read_h5_metadata(paths$h5)[source_row %in% meta$source_row_h5, .(
  source_row,
  nCount_RNA = as.numeric(nCount_RNA)
)]
libsize_all_genes <- libsize_all_genes[match(meta$source_row_h5, source_row), nCount_RNA]
fallback <- as.numeric(Matrix::colSums(counts))
libsize_all_genes[!is.finite(libsize_all_genes) | libsize_all_genes <= 0] <- fallback[
  !is.finite(libsize_all_genes) | libsize_all_genes <= 0
]
logdata <- counts %*% Diagonal(x = 1e4 / libsize_all_genes)
dimnames(logdata) <- dimnames(counts)
logdata@x <- log1p(logdata@x)
global_max <- max(logdata@x)
if (!is.finite(global_max) || global_max <= 0) stop("Invalid global expression scale")
logdata <- logdata / global_max

tri_mean <- function(x) {
  q <- as.numeric(stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE, type = 7))
  (q[[1]] + 2 * q[[2]] + q[[3]]) / 4
}

score_one_donor <- function(condition_oi, source_oi, donor_oi) {
  sender_cells <- meta[
    condition == condition_oi & cell_group == source_oi & donor_id == donor_oi,
    cell_id
  ]
  receiver_label <- if (condition_oi == "Healthy") "Healthy ASPC-like" else "LoSc Q8"
  receiver_cells <- meta[
    condition == condition_oi & receiver_state == receiver_label & donor_id == donor_oi,
    cell_id
  ]
  if (length(sender_cells) < min_source_cells || length(receiver_cells) < min_receiver_cells) return(NULL)

  avg_sender <- apply(as.matrix(logdata[, sender_cells, drop = FALSE]), 1, tri_mean)
  avg_receiver <- apply(as.matrix(logdata[, receiver_cells, drop = FALSE]), 1, tri_mean)
  avg <- cbind(sender = avg_sender, receiver = avg_receiver)
  rownames(avg) <- rownames(logdata)
  pair_lr <- as.data.frame(lr_universe[source == source_oi])
  rownames(pair_lr) <- pair_lr$interaction_name
  ligand_avg <- computeExpr_LR(pair_lr$ligand, avg, complex_input)
  receptor_avg <- computeExpr_LR(pair_lr$receptor, avg, complex_input)
  lr_product <- ligand_avg[, 1L] * receptor_avg[, 2L]
  probability <- lr_product / (0.5 + lr_product)
  detail <- data.table(
    condition = condition_oi,
    source = source_oi,
    donor_id = donor_oi,
    interaction_name = pair_lr$interaction_name,
    pathway_name = pair_lr$pathway_name,
    ligand = pair_lr$ligand,
    receptor = pair_lr$receptor,
    ligand_average = ligand_avg[, 1L],
    receptor_average = receptor_avg[, 2L],
    communication_probability = as.numeric(probability),
    n_source_cells = length(sender_cells),
    n_receiver_cells = length(receiver_cells)
  )
  list(
    summary = data.table(
      condition = condition_oi,
      source = source_oi,
      donor_id = donor_oi,
      donor_total_strength = sum(probability),
      n_lr_pairs = nrow(pair_lr),
      n_source_cells = length(sender_cells),
      n_receiver_cells = length(receiver_cells)
    ),
    detail = detail
  )
}

message("Computing donor-level CellChat-informed total interaction strength")
res <- list()
idx <- 0L
for (condition_oi in condition_levels) {
  receiver_label <- if (condition_oi == "Healthy") "Healthy ASPC-like" else "LoSc Q8"
  receiver_donors <- unique(meta[condition == condition_oi & receiver_state == receiver_label, donor_id])
  for (source_oi in source_levels) {
    source_donors <- unique(meta[condition == condition_oi & cell_group == source_oi, donor_id])
    for (donor_oi in intersect(receiver_donors, source_donors)) {
      z <- score_one_donor(condition_oi, source_oi, donor_oi)
      if (!is.null(z)) {
        idx <- idx + 1L
        res[[idx]] <- z
      }
    }
  }
}
donor_scores <- rbindlist(lapply(res, function(x) x[["summary"]]))
lr_scores <- rbindlist(lapply(res, function(x) x[["detail"]]))
fwrite(donor_scores, file.path(source_dir, "CellChat_donor_level_total_strength.csv"))
fwrite(lr_scores, file.path(source_dir, "CellChat_donor_level_LR_scores.csv.gz"))

exact_permutation <- function(healthy, q8) {
  values <- c(healthy, q8)
  n_q8 <- length(q8)
  n_healthy <- length(healthy)
  observed <- mean(q8) - mean(healthy)
  combos <- combn(length(values), n_q8)
  q8_means <- colMeans(matrix(values[combos], nrow = n_q8))
  healthy_means <- (sum(values) - n_q8 * q8_means) / n_healthy
  null_difference <- q8_means - healthy_means
  tolerance <- 1e-12
  list(
    difference = observed,
    p_two_sided = mean(abs(null_difference) >= abs(observed) - tolerance),
    p_q8_higher = mean(null_difference >= observed - tolerance),
    n_permutations = length(null_difference)
  )
}

message("Testing Q8 versus Healthy across biological donors")
stats <- rbindlist(lapply(source_levels, function(source_oi) {
  h <- donor_scores[source == source_oi & condition == "Healthy", donor_total_strength]
  q <- donor_scores[source == source_oi & condition == "Q8", donor_total_strength]
  test <- exact_permutation(h, q)
  cliffs_delta <- mean(outer(q, h, `>`)) - mean(outer(q, h, `<`))
  data.table(
    source = source_oi,
    n_healthy = length(h),
    n_q8 = length(q),
    healthy_mean = mean(h),
    q8_mean = mean(q),
    mean_difference = test$difference,
    fold_change = (mean(q) + 1e-8) / (mean(h) + 1e-8),
    p_two_sided = test$p_two_sided,
    p_q8_higher = test$p_q8_higher,
    n_exact_permutations = test$n_permutations,
    cliffs_delta = cliffs_delta
  )
}))
stats[, q_two_sided := p.adjust(p_two_sided, method = "BH")]
stats[, q_q8_higher := p.adjust(p_q8_higher, method = "BH")]
stats[, significance := fcase(
  q_two_sided < 0.001, "***",
  q_two_sided < 0.01, "**",
  q_two_sided < 0.05, "*",
  default = "ns"
)]
fwrite(stats, file.path(source_dir, "CellChat_Q8_vs_Healthy_donor_level_statistics.csv"))

summary_plot <- donor_scores[, .(
  mean_strength = mean(donor_total_strength),
  sem = sd(donor_total_strength) / sqrt(.N),
  n = .N
), by = .(source, condition)]
summary_plot[, source_index := match(source, source_levels)]
summary_plot[, x := source_index + fifelse(condition == "Healthy", -0.19, 0.19)]
donor_scores[, source_index := match(source, source_levels)]
donor_scores[, x := source_index + fifelse(condition == "Healthy", -0.19, 0.19)]
setorder(donor_scores, source_index, condition, donor_id)
donor_scores[, x_jitter := x + if (.N == 1L) 0 else seq(-0.055, 0.055, length.out = .N),
             by = .(source, condition)]

brackets <- merge(
  stats,
  donor_scores[, .(y_max = max(donor_total_strength)), by = source],
  by = "source"
)
brackets[, source_index := match(source, source_levels)]
brackets[, y := y_max * 1.18 + max(donor_scores$donor_total_strength) * 0.025]
brackets[, label := paste0(significance, "  q=", formatC(q_two_sided, format = "f", digits = 3))]
y_limit <- max(brackets$y) * 1.15

p <- ggplot() +
  geom_col(
    data = summary_plot,
    aes(x, mean_strength, fill = condition),
    width = 0.32,
    colour = "#222222",
    linewidth = 0.3
  ) +
  geom_errorbar(
    data = summary_plot,
    aes(x, ymin = pmax(0, mean_strength - sem), ymax = mean_strength + sem),
    width = 0.10,
    linewidth = 0.4,
    colour = "#222222"
  ) +
  geom_point(
    data = donor_scores,
    aes(x_jitter, donor_total_strength, fill = condition),
    shape = 21,
    size = 1.65,
    stroke = 0.35,
    colour = "white",
    alpha = 0.92
  ) +
  geom_segment(
    data = brackets,
    aes(x = source_index - 0.19, xend = source_index + 0.19, y = y, yend = y),
    linewidth = 0.4,
    colour = "#222222"
  ) +
  geom_segment(
    data = brackets,
    aes(x = source_index - 0.19, xend = source_index - 0.19,
        y = y - 0.012 * y_limit, yend = y),
    linewidth = 0.4,
    colour = "#222222"
  ) +
  geom_segment(
    data = brackets,
    aes(x = source_index + 0.19, xend = source_index + 0.19,
        y = y - 0.012 * y_limit, yend = y),
    linewidth = 0.4,
    colour = "#222222"
  ) +
  geom_text(
    data = brackets,
    aes(source_index, y + 0.018 * y_limit, label = label),
    family = "Helvetica",
    size = 2.35,
    colour = "#111111"
  ) +
  scale_fill_manual(values = condition_palette, breaks = condition_levels) +
  scale_x_continuous(breaks = seq_along(source_levels), labels = source_levels) +
  scale_y_continuous(limits = c(0, y_limit), expand = expansion(mult = c(0, 0))) +
  labs(
    title = "Donor-level CD4-to-ASPC interaction strength",
    subtitle = "Q8 versus Healthy; bars show mean +/- SEM and points denote biological donors",
    x = NULL,
    y = "CellChat-informed total interaction strength",
    caption = paste0(
      "Fixed LR universe from matched pooled CellChat; donors required >=", min_source_cells,
      " source and >=", min_receiver_cells,
      " receiver cells. Exact two-sided label-permutation test; BH correction across five CD4 states."
    )
  ) +
  theme_nature(7) +
  theme(
    legend.position = "top",
    legend.justification = "left",
    axis.text.x = element_text(angle = 22, hjust = 1, face = "bold"),
    plot.caption.position = "plot"
  )

stem <- file.path(output_dir, "CellChat_Q8_vs_Healthy_donor_level_total_strength")
width_mm <- 183
height_mm <- 78
svglite::svglite(paste0(stem, ".svg"), width = width_mm / 25.4, height = height_mm / 25.4)
print(p)
grDevices::dev.off()
if (isTRUE(getOption("cellchat.use_cairo_pdf", FALSE))) {
  grDevices::cairo_pdf(
    filename = paste0(stem, ".pdf"), width = width_mm / 25.4,
    height = height_mm / 25.4, family = "Helvetica", onefile = TRUE
  )
} else {
  quartz(
    type = "pdf", file = paste0(stem, ".pdf"), width = width_mm / 25.4,
    height = height_mm / 25.4, family = "Helvetica"
  )
}
print(p)
grDevices::dev.off()
ragg::agg_tiff(paste0(stem, ".tiff"), width = width_mm / 25.4,
               height = height_mm / 25.4, units = "in", res = 600,
               compression = "lzw", background = "white")
print(p)
grDevices::dev.off()
ragg::agg_png(paste0(stem, "_preview.png"), width = width_mm / 25.4,
              height = height_mm / 25.4, units = "in", res = 300,
              background = "white")
print(p)
grDevices::dev.off()

writeLines(c(
  "CellChat donor-level significance bar plot QA notes",
  "Core claim: donor-level CellChat-informed CD4-to-ASPC total strength differs between Q8 and Healthy.",
  paste0("Eligibility: >=", min_source_cells, " source cells and >=", min_receiver_cells, " receiver cells per donor."),
  "LR universe: union of significant interactions in the matched pooled Healthy and Q8 CellChat models, fixed before donor scoring.",
  "Score: sum of CellChat mass-action probabilities using donor-level triMean expression and one shared expression scale.",
  "Center/spread: mean +/- SEM; points are biological donors.",
  "Test: exact two-sided permutation of condition labels; BH correction across five CD4-state comparisons.",
  "The donor-level score is CellChat-informed and follows its mass-action transformation, but it is not a direct output of compareInteractions()."
), file.path(output_dir, "CellChat_Q8_donor_significance_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "CellChat_Q8_donor_significance_sessionInfo.txt"))
message("Done: ", normalizePath(output_dir))
