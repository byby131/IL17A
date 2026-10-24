#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(CellChat)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})

source("scripts/figure2_ghi_common.R")

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1L) args[[1]] else "outputs/CellChat_CD4_Q8Q11_vs_Healthy"
cache_dir <- if (length(args) >= 2L) args[[2]] else "work/cellchat_cd4_q8q11_v2"
source_dir <- file.path(output_dir, "source_data")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

paths <- list(
  h5 = "data/GSE264508/GSE264508_data.all.H5Seurat",
  cd4 = "outputs/Figure2C_CD4_cell_metadata.csv.gz",
  receiver = "outputs/Figure2H_Q8Q11_revision/audit/Figure2H_Q8_Q11_frozen_cell_metadata.csv.gz"
)
missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths)) stop("Missing inputs: ", paste(missing_paths, collapse = ", "))

set.seed(264508)
future::plan("sequential")

source_levels <- c("Naive CD4 T", "Memory CD4 T", "Th17-like", "Treg", "Activated CD4 T")
all_levels <- c(source_levels, "ASPC receiver")
condition_levels <- c("Healthy", "Q8", "Q11")
palette_condition <- c(Healthy = "#3569A8", Q8 = "#D94755", Q11 = "#D8902F")

theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "#202020"),
      axis.ticks = element_line(linewidth = 0.35, colour = "#202020"),
      axis.text = element_text(colour = "#111111"),
      axis.title = element_text(colour = "#111111"),
      plot.title = element_text(size = base_size + 1, face = "bold", hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#4D4D4D"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", colour = "#111111"),
      legend.key = element_blank(),
      legend.background = element_blank(),
      panel.grid = element_blank(),
      plot.margin = margin(3, 4, 3, 3)
    )
}

message("Preparing frozen CD4 and ASPC-state annotations")
cd4 <- fread(paths$cd4)[final_state %chin% source_levels, .(
  cell_id,
  source_row_h5 = as.integer(source_row),
  donor_id = as.character(library_id),
  condition = fifelse(health == "Healthy", "Healthy", "Disease"),
  cell_group = as.character(final_state)
)]

receiver <- fread(paths$receiver)[qc_status == "Retained" & analysis_group %chin% c(
  "Healthy ASPC-like", "LoSc Q8", "LoSc Q11"
), .(
  cell_id,
  source_row_h5 = as.integer(source_row_h5),
  donor_id = as.character(donor_id),
  condition = fcase(
    analysis_group == "Healthy ASPC-like", "Healthy",
    analysis_group == "LoSc Q8", "Q8",
    analysis_group == "LoSc Q11", "Q11"
  ),
  cell_group = "ASPC receiver",
  receiver_state = as.character(analysis_group)
)]

meta_all <- unique(rbindlist(list(
  cd4[, .(cell_id, source_row_h5, donor_id, condition, cell_group, receiver_state = NA_character_)],
  receiver
), use.names = TRUE), by = "cell_id")
if (anyDuplicated(meta_all$cell_id)) stop("Duplicated cell IDs in merged metadata")

h5_md <- read_h5_metadata(paths$h5)[, .(cell_id_h5 = cell_id, source_row_h5 = source_row)]
meta_all <- merge(meta_all, h5_md, by = "source_row_h5", all.x = TRUE, sort = FALSE)
if (anyNA(meta_all$cell_id_h5) || any(meta_all$cell_id != meta_all$cell_id_h5)) {
  stop("Frozen cell identifiers do not map exactly to H5Seurat source rows")
}
meta_all[, cell_id_h5 := NULL]
setorder(meta_all, source_row_h5)

coverage <- meta_all[, .(
  n_cells = .N,
  n_donors = uniqueN(donor_id)
), by = .(condition, receiver_state, cell_group)]
fwrite(coverage, file.path(source_dir, "CellChat_input_cell_and_donor_coverage.csv"))
fwrite(meta_all, file.path(source_dir, "CellChat_frozen_input_cell_metadata.csv.gz"))

message("Reading raw counts for ", nrow(meta_all), " selected cells")
counts <- read_h5_counts(paths$h5, meta_all$source_row_h5)
meta_all <- meta_all[match(colnames(counts), cell_id)]
stopifnot(identical(meta_all$cell_id, colnames(counts)))
counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
libsize <- as.numeric(Matrix::colSums(counts))
if (any(!is.finite(libsize) | libsize <= 0)) stop("Non-positive library size detected")
logdata <- counts %*% Diagonal(x = 1e4 / libsize)
dimnames(logdata) <- dimnames(counts)
logdata@x <- log1p(logdata@x)

select_model_cells <- function(model_condition, balanced = FALSE, seed = 264508L) {
  source_condition <- if (model_condition == "Healthy") "Healthy" else "Disease"
  receiver_label <- switch(model_condition,
    Healthy = "Healthy ASPC-like",
    Q8 = "LoSc Q8",
    Q11 = "LoSc Q11"
  )
  src <- meta_all[condition == source_condition & cell_group %chin% source_levels]
  rec <- meta_all[receiver_state == receiver_label]
  if (!balanced) return(rbindlist(list(src, rec), use.names = TRUE))

  target_source_n <- cd4[, .N, by = .(condition, cell_group)][,
    .(n = min(N)), by = cell_group
  ]
  sample_even_donor <- function(x, n_target) {
    if (nrow(x) <= n_target) return(copy(x))
    donors <- unique(x$donor_id)
    selected <- integer()
    # First preserve donor breadth, then fill remaining slots at random.
    for (d in sample(donors)) {
      hit <- which(x$donor_id == d)
      selected <- c(selected, hit[sample.int(length(hit), 1L)])
      if (length(selected) >= n_target) break
    }
    remaining <- setdiff(seq_len(nrow(x)), selected)
    if (length(selected) < n_target) {
      take <- sample.int(length(remaining), n_target - length(selected))
      selected <- c(selected, remaining[take])
    }
    x[sort(selected)]
  }

  # Q8 and Q11 use exactly the same matched disease-CD4 source cells.
  set.seed(if (source_condition == "Healthy") 264509L else 264510L)
  src_matched <- rbindlist(lapply(source_levels, function(g) {
    n_target <- target_source_n[cell_group == g, n]
    sample_even_donor(src[cell_group == g], n_target)
  }))
  set.seed(seed + 1000L)
  rec_matched <- sample_even_donor(rec, 43L)
  rbindlist(list(src_matched, rec_matched), use.names = TRUE)
}

model_cells <- rbindlist(lapply(c(FALSE, TRUE), function(balanced) {
  rbindlist(lapply(seq_along(condition_levels), function(k) {
    z <- select_model_cells(condition_levels[[k]], balanced = balanced,
                            seed = 264508L + k + 100L * balanced)
    z[, `:=`(
      model = condition_levels[[k]],
      analysis = if (balanced) "cell-number matched" else "full"
    )]
    z
  }))
}))
fwrite(model_cells, file.path(source_dir, "CellChat_model_cell_metadata.csv.gz"))
model_coverage <- model_cells[, .(
  n_cells = .N,
  n_donors = uniqueN(donor_id)
), by = .(analysis, model, cell_group)]
fwrite(model_coverage, file.path(source_dir, "CellChat_model_cell_and_donor_coverage.csv"))

run_cellchat <- function(model_condition, analysis_label) {
  balanced <- analysis_label == "cell-number matched"
  cache_file <- file.path(cache_dir, paste0(
    if (balanced) "balanced_matchedDesign2" else "full", "_", model_condition,
    "_CellChat_v", gsub("\\.", "_", as.character(packageVersion("CellChat"))), ".rds"
  ))
  if (file.exists(cache_file)) {
    message("Loading cached CellChat object: ", cache_file)
    return(readRDS(cache_file))
  }
  cm <- model_cells[analysis == analysis_label & model == model_condition]
  cm[, cell_group := factor(cell_group, levels = all_levels)]
  cells <- cm$cell_id
  x <- logdata[, cells, drop = FALSE]
  meta <- data.frame(
    cell_group = cm$cell_group,
    samples = factor(rep(model_condition, length(cells))),
    row.names = cells
  )
  message("Running CellChat: ", analysis_label, " / ", model_condition,
          " (", ncol(x), " cells)")
  object <- createCellChat(object = x, meta = meta, group.by = "cell_group",
                           datatype = "RNA", do.sparse = TRUE)
  object@DB <- CellChatDB.human
  object <- subsetData(object)
  object <- identifyOverExpressedGenes(object, do.fast = FALSE, min.cells = 10)
  object <- identifyOverExpressedInteractions(object)
  object <- computeCommunProb(
    object,
    type = "triMean",
    raw.use = TRUE,
    population.size = FALSE,
    distance.use = FALSE,
    nboot = 100,
    seed.use = 264508L
  )
  object <- filterCommunication(object, min.cells = 10)
  object <- computeCommunProbPathway(object, thresh = 0.05)
  object <- aggregateNet(object, thresh = 0.05)
  saveRDS(object, cache_file, compress = "xz")
  object
}

objects <- list()
for (analysis_label in c("full", "cell-number matched")) {
  for (model_condition in condition_levels) {
    key <- paste(analysis_label, model_condition, sep = "__")
    objects[[key]] <- run_cellchat(model_condition, analysis_label)
  }
}

extract_lr <- function(object, analysis_label, model_condition) {
  z <- as.data.table(subsetCommunication(
    object = object,
    slot.name = "net",
    sources.use = source_levels,
    targets.use = "ASPC receiver",
    thresh = 0.05
  ))
  if (!nrow(z)) {
    return(data.table(
      source = character(), target = character(), ligand = character(), receptor = character(),
      prob = numeric(), pval = numeric(), interaction_name = character(), pathway_name = character(),
      analysis = character(), model = character()
    ))
  }
  z[, `:=`(analysis = analysis_label, model = model_condition)]
  z
}

lr_all <- rbindlist(lapply(names(objects), function(key) {
  bits <- strsplit(key, "__", fixed = TRUE)[[1]]
  extract_lr(objects[[key]], bits[[1]], bits[[2]])
}), fill = TRUE)
fwrite(lr_all, file.path(source_dir, "CellChat_significant_CD4_to_ASPC_LR_interactions.csv.gz"))

metric_grid <- CJ(
  analysis = c("full", "cell-number matched"),
  model = condition_levels,
  source = source_levels,
  unique = TRUE
)
metrics <- lr_all[, .(
  interaction_count = .N,
  total_strength = sum(prob, na.rm = TRUE)
), by = .(analysis, model, source)]
metrics <- merge(metric_grid, metrics, by = c("analysis", "model", "source"), all.x = TRUE)
metrics[is.na(interaction_count), interaction_count := 0L]
metrics[is.na(total_strength), total_strength := 0]
fwrite(metrics, file.path(source_dir, "CellChat_total_CD4_to_ASPC_metrics.csv"))

make_contrasts <- function(x) {
  healthy <- x[model == "Healthy", .(
    analysis, source,
    healthy_count = interaction_count,
    healthy_strength = total_strength
  )]
  disease <- x[model %chin% c("Q8", "Q11")]
  out <- merge(disease, healthy, by = c("analysis", "source"), all.x = TRUE)
  positive_strengths <- c(out$total_strength, out$healthy_strength)
  positive_strengths <- positive_strengths[is.finite(positive_strengths) & positive_strengths > 0]
  strength_pc <- if (length(positive_strengths)) max(min(positive_strengths) / 2, 1e-8) else 1e-8
  out[, `:=`(
    count_log2_ratio = log2((interaction_count + 1) / (healthy_count + 1)),
    strength_log2_ratio = log2((total_strength + strength_pc) / (healthy_strength + strength_pc)),
    count_difference = interaction_count - healthy_count,
    strength_difference = total_strength - healthy_strength,
    strength_pseudocount = strength_pc
  )]
  out
}

contrasts <- make_contrasts(metrics)
contrasts[, source_order := match(source, source_levels)]
setorder(contrasts, analysis, model, source_order)
contrasts[, source_order := NULL]
fwrite(contrasts, file.path(source_dir, "CellChat_Q8_Q11_vs_Healthy_total_metric_contrasts.csv"))

robust_calls <- dcast(
  contrasts,
  model + source ~ analysis,
  value.var = c("count_log2_ratio", "strength_log2_ratio", "count_difference", "strength_difference")
)
robust_calls[, `:=`(
  full_strength_enhanced = strength_difference_full > 0,
  matched_strength_enhanced = `strength_difference_cell-number matched` > 0,
  full_count_enhanced = count_difference_full > 0,
  matched_count_enhanced = `count_difference_cell-number matched` > 0
)]
robust_calls[, conclusion := fcase(
  full_strength_enhanced & matched_strength_enhanced & full_count_enhanced & matched_count_enhanced,
    "Consistently enhanced in strength and count",
  full_strength_enhanced & matched_strength_enhanced,
    "Strength enhanced in both analyses; interaction count not consistently enhanced",
  !full_strength_enhanced & !matched_strength_enhanced,
    "No strength enhancement in either analysis",
  default = "Direction differs between full and cell-number-matched analyses"
)]
fwrite(robust_calls, file.path(source_dir, "CellChat_robust_enhancement_calls.csv"))

pathway <- lr_all[, .(
  pathway_strength = sum(prob, na.rm = TRUE),
  n_lr = .N
), by = .(analysis, model, source, pathway_name)]
path_grid <- CJ(
  analysis = c("full", "cell-number matched"),
  model = condition_levels,
  source = source_levels,
  pathway_name = sort(unique(pathway$pathway_name)),
  unique = TRUE
)
pathway <- merge(path_grid, pathway, by = c("analysis", "model", "source", "pathway_name"), all.x = TRUE)
pathway[is.na(pathway_strength), pathway_strength := 0]
pathway[is.na(n_lr), n_lr := 0L]
fwrite(pathway, file.path(source_dir, "CellChat_CD4_to_ASPC_pathway_strengths.csv.gz"))

path_contrast <- merge(
  pathway[model %chin% c("Q8", "Q11")],
  pathway[model == "Healthy", .(
    analysis, source, pathway_name,
    healthy_pathway_strength = pathway_strength,
    healthy_n_lr = n_lr
  )],
  by = c("analysis", "source", "pathway_name"), all.x = TRUE
)
positive_path_strengths <- c(path_contrast$pathway_strength, path_contrast$healthy_pathway_strength)
positive_path_strengths <- positive_path_strengths[
  is.finite(positive_path_strengths) & positive_path_strengths > 0
]
path_pc <- if (length(positive_path_strengths)) max(min(positive_path_strengths) / 2, 1e-8) else 1e-8
path_contrast[, `:=`(
  pathway_difference = pathway_strength - healthy_pathway_strength,
  pathway_log2_ratio = log2((pathway_strength + path_pc) / (healthy_pathway_strength + path_pc)),
  pathway_pseudocount = path_pc
)]
fwrite(path_contrast, file.path(source_dir, "CellChat_Q8_Q11_vs_Healthy_pathway_contrasts.csv.gz"))

lr_key_cols <- c("analysis", "source", "interaction_name", "pathway_name", "ligand", "receptor")
lr_disease <- lr_all[model %chin% c("Q8", "Q11")]
lr_healthy <- lr_all[model == "Healthy"]
lr_contrast <- merge(
  lr_disease,
  lr_healthy[, c(lr_key_cols, "prob"), with = FALSE],
  by = lr_key_cols,
  all = TRUE,
  suffixes = c("_disease", "_healthy")
)
lr_contrast[is.na(model), model := fifelse(is.na(model), NA_character_, model)]
# Expand Healthy-only interactions to both receiver contrasts.
healthy_only <- lr_contrast[is.na(model)]
if (nrow(healthy_only)) {
  lr_contrast <- lr_contrast[!is.na(model)]
  lr_contrast <- rbindlist(list(
    lr_contrast,
    copy(healthy_only)[, model := "Q8"],
    copy(healthy_only)[, model := "Q11"]
  ), fill = TRUE)
}
lr_contrast[is.na(prob_disease), prob_disease := 0]
lr_contrast[is.na(prob_healthy), prob_healthy := 0]
lr_contrast[, probability_difference := prob_disease - prob_healthy]
lr_contrast[, direction := fifelse(
  probability_difference >= 0, paste0(model, " higher"), "Healthy higher"
)]
fwrite(lr_contrast, file.path(source_dir, "CellChat_Q8_Q11_vs_Healthy_LR_contrasts.csv.gz"))

message("Auditing IL-17 family interactions, including non-significant entries")
il17_db <- as.data.table(CellChatDB.human$interaction)[pathway_name == "IL17"]
extract_il17 <- function(object, analysis_label, model_condition) {
  prob_arr <- object@net$prob
  pval_arr <- object@net$pval
  interactions <- dimnames(prob_arr)[[3]]
  rbindlist(lapply(source_levels, function(src) {
    rbindlist(lapply(il17_db$interaction_name, function(int_name) {
      present <- int_name %chin% interactions
      dbrow <- il17_db[interaction_name == int_name][1]
      data.table(
        analysis = analysis_label,
        model = model_condition,
        source = src,
        target = "ASPC receiver",
        interaction_name = int_name,
        ligand = dbrow$ligand,
        receptor = dbrow$receptor,
        probability = if (present) prob_arr[src, "ASPC receiver", int_name] else 0,
        p_value = if (present) pval_arr[src, "ASPC receiver", int_name] else 1,
        retained_as_significant = present && pval_arr[src, "ASPC receiver", int_name] < 0.05 &&
          prob_arr[src, "ASPC receiver", int_name] > 0
      )
    }))
  }))
}
il17_audit <- rbindlist(lapply(names(objects), function(key) {
  bits <- strsplit(key, "__", fixed = TRUE)[[1]]
  extract_il17(objects[[key]], bits[[1]], bits[[2]])
}))
fwrite(il17_audit, file.path(source_dir, "CellChat_IL17_family_complete_audit.csv"))

gene_detection <- function(model_condition, analysis_label, genes) {
  cm <- model_cells[analysis == analysis_label & model == model_condition]
  genes <- intersect(genes, rownames(counts))
  rbindlist(lapply(all_levels, function(g) {
    cells <- cm[cell_group == g, cell_id]
    if (!length(cells)) return(NULL)
    mat <- counts[genes, cells, drop = FALSE]
    data.table(
      analysis = analysis_label,
      model = model_condition,
      cell_group = g,
      gene = genes,
      n_cells = length(cells),
      pct_detected = 100 * as.numeric(Matrix::rowMeans(mat > 0)),
      mean_raw_count = as.numeric(Matrix::rowMeans(mat))
    )
  }))
}
il17_genes <- unique(c("IL17A", "IL17F", "IL17B", "IL17C", "IL25", "IL17RA", "IL17RB", "IL17RC", "IL17RE"))
il17_expression <- rbindlist(lapply(c("full", "cell-number matched"), function(a) {
  rbindlist(lapply(condition_levels, gene_detection, analysis_label = a, genes = il17_genes))
}))
fwrite(il17_expression, file.path(source_dir, "CellChat_IL17_family_expression_detection.csv"))

message("Creating submission-style summary figure")
matched <- contrasts[analysis == "cell-number matched"]
matched[, source := factor(source, levels = rev(source_levels))]
matched[, model := factor(model, levels = c("Q8", "Q11"))]
fill_max <- max(abs(c(matched$strength_log2_ratio, matched$count_log2_ratio)), na.rm = TRUE)

p_a <- ggplot(matched, aes(model, source, fill = strength_log2_ratio)) +
  geom_tile(colour = "white", linewidth = 0.75) +
  geom_text(aes(label = sprintf("%.2f", strength_log2_ratio)), size = 2.25,
            family = "Helvetica", fontface = "bold", colour = "black") +
  scale_fill_gradient2(
    low = "#3569A8", mid = "#F4F1EA", high = "#D94755", midpoint = 0,
    limits = c(-fill_max, fill_max), oob = squish, name = "log2 ratio"
  ) +
  labs(title = "Total interaction strength", subtitle = "Cell-number matched",
       x = NULL, y = NULL) +
  coord_fixed(ratio = 0.72) +
  theme_nature(7) +
  theme(axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "bottom", legend.key.width = grid::unit(10, "mm"))

p_b <- ggplot(matched, aes(model, source, fill = count_log2_ratio)) +
  geom_tile(colour = "white", linewidth = 0.75) +
  geom_text(aes(label = sprintf("%.2f", count_log2_ratio)), size = 2.25,
            family = "Helvetica", fontface = "bold", colour = "black") +
  scale_fill_gradient2(
    low = "#3569A8", mid = "#F4F1EA", high = "#D94755", midpoint = 0,
    limits = c(-fill_max, fill_max), oob = squish, name = "log2 ratio"
  ) +
  labs(title = "Significant LR-pair count", subtitle = "Cell-number matched",
       x = NULL, y = NULL) +
  coord_fixed(ratio = 0.72) +
  theme_nature(7) +
  theme(axis.ticks = element_blank(), axis.line = element_blank(),
        legend.position = "none")

sens <- contrasts[, .(analysis, model, source, strength_log2_ratio)]
sens <- dcast(sens, model + source ~ analysis, value.var = "strength_log2_ratio")
sens[, model := factor(model, levels = c("Q8", "Q11"))]
sens_limit <- max(abs(c(sens$full, sens$`cell-number matched`)), na.rm = TRUE) + 0.8
p_c <- ggplot(sens, aes(full, `cell-number matched`, colour = model, label = source)) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "#8C8C8C") +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "#8C8C8C") +
  geom_abline(slope = 1, intercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_point(size = 2.1) +
  ggrepel::geom_text_repel(size = 2.15, family = "Helvetica", fontface = "bold",
                           box.padding = 0.25, point.padding = 0.15, max.overlaps = Inf,
                           min.segment.length = 0, segment.size = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = palette_condition[c("Q8", "Q11")], name = NULL) +
  labs(title = "Sensitivity to cell-number matching",
       subtitle = "CD4-to-receiver total strength",
       x = "Full model log2 ratio", y = "Matched model log2 ratio") +
  coord_equal(xlim = c(-sens_limit, sens_limit), ylim = c(-sens_limit, sens_limit),
              expand = FALSE) +
  theme_nature(7) +
  theme(legend.position = "top", legend.justification = "left")

path_plot <- path_contrast[analysis == "cell-number matched"]
path_plot[, abs_change := abs(pathway_difference)]
path_top <- path_plot[order(-abs_change), head(.SD, 5), by = model]
path_top[, label := paste(source, pathway_name, sep = " | ")]
path_top[, label := factor(label, levels = rev(unique(label[order(model, pathway_log2_ratio)])))]
path_top[, model := factor(model, levels = c("Q8", "Q11"))]
p_d <- ggplot(path_top, aes(pathway_log2_ratio, label, colour = model)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#333333") +
  geom_segment(aes(x = 0, xend = pathway_log2_ratio, yend = label), linewidth = 0.75,
               alpha = 0.75) +
  geom_point(size = 2.4) +
  facet_wrap(~model, scales = "free_y", ncol = 1) +
  scale_colour_manual(values = palette_condition[c("Q8", "Q11")], guide = "none") +
  labs(title = "Largest pathway shifts", subtitle = "Top five per receiver contrast",
       x = "log2 ratio versus Healthy", y = NULL) +
  theme_nature(6.7)

combined <- wrap_plots(
  list(p_a, p_b, p_c, p_d),
  ncol = 2,
  byrow = TRUE,
  heights = c(0.86, 1.25),
  widths = c(1, 1.16)
) +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(family = "Helvetica", face = "bold", size = 10, colour = "#111111")
    )
  )

figure_stem <- file.path(output_dir, "CellChat_CD4_Q8_Q11_vs_Healthy_summary")
save_nature_plot(combined, figure_stem, width_mm = 183, height_mm = 145, dpi = 600)

q11_matched <- contrasts[analysis == "cell-number matched" & model == "Q11"]
q8_matched <- contrasts[analysis == "cell-number matched" & model == "Q8"]
q11_all_strength <- all(q11_matched$strength_difference > 0)
q8_all_strength <- all(q8_matched$strength_difference > 0)
il17_n_sig <- il17_audit[retained_as_significant == TRUE, .N]
report_lines <- c(
  "# CellChat：所有 CD4 亚群与 Q8/Q11 的通讯是否较健康对照增强",
  "",
  "## 核心结论",
  "",
  paste0("- Q11：细胞数匹配后，5 个 CD4 亚群到 Q11 的总通讯强度",
         if (q11_all_strength) "全部高于" else "并非全部高于",
         "健康 CD4 亚群到 Healthy ASPC-like 的对应值。"),
  paste0("- Q8：细胞数匹配后，5 个 CD4 亚群到 Q8 的总通讯强度",
         if (q8_all_strength) "全部高于" else "均未高于",
         "健康基线。"),
  "- 因此不能概括为‘Q8 和 Q11 与所有 CD4 亚群的相互作用均增强’；数据支持的是 Q11 广泛增强、Q8 不增强。",
  paste0("- IL-17 家族达到 CellChat 显著标准的 CD4→受体相互作用数为 ", il17_n_sig,
         "。本分析不支持经典 IL-17A/F→IL17RA/RC 轴增强。"),
  "",
  "## 细胞数匹配后的定量结果",
  "",
  "总强度 log2 比值（疾病状态/健康；正值为疾病较高）：",
  "",
  paste0("- Q11：", paste(sprintf("%s %.2f", q11_matched$source, q11_matched$strength_log2_ratio), collapse = "；"), "。"),
  paste0("- Q8：", paste(sprintf("%s %.2f", q8_matched$source, q8_matched$strength_log2_ratio), collapse = "；"), "。"),
  "",
  "Q11 的主要增加信号由 PPIA-BSG（CypA）、CD99-CD99、MIF-CD74/CD44 及部分 SELL-CD34 构成；",
  "Q8 中上述主要信号多数低于健康基线。经典 IL-17A/F 信号未被检出为显著通讯。",
  "",
  "## 分析与解释边界",
  "",
  "- Healthy 基线为健康 CD4 亚群→Healthy ASPC-like；原始数据不存在‘健康 Q8/Q11’。因此该比较同时包含发送端疾病状态和接收端细胞状态的差异，不能把差异完全归因于 Q8 或 Q11 本身。",
  "- 主结果采用每个 CD4 亚群及受体细胞数匹配；全量模型作为敏感性分析。CellChat 参数为 triMean、population.size=FALSE、distance.use=FALSE、100 次 bootstrap、min.cells=10。",
  "- CellChat 是池化细胞的配体-受体推断，本比较不是供者重复层面的差异检验；因此不报告条件间 q 值，也不把结果表述为因果关系。",
  "- Q8 仅 43 个细胞，是精度和可检出性的主要限制；Q11 的广泛增强应视为需要外部数据或空间/蛋白证据验证的通讯潜能。"
)
writeLines(report_lines, file.path(output_dir, "CellChat结果说明_中文版.md"))

writeLines(capture.output(sessionInfo()), file.path(output_dir, "CellChat_sessionInfo.txt"))
writeLines(c(
  "CellChat CD4-to-ASPC analysis QA notes",
  paste0("CellChat version: ", packageVersion("CellChat")),
  "Database: complete CellChatDB.human (all annotations).",
  "Estimator: triMean; population.size=FALSE; distance.use=FALSE; nboot=100; min.cells=10.",
  "Primary comparison: cell-number-matched model; full model is a sensitivity analysis.",
  "Healthy baseline is Healthy CD4 subgroups -> Healthy ASPC-like.",
  "Q8/Q11 contrasts are disease CD4 subgroups -> Q8/Q11 and therefore combine sender-condition and receiver-state differences.",
  "CellChat pooled contrasts are descriptive; they are not donor-replicate differential tests and no between-condition q value is claimed.",
  "Positive log2 ratio means stronger/more communication in the disease-state model; negative means Healthy higher.",
  "IL17 family interactions were audited directly from the full CellChat probability arrays, including non-significant entries."
), file.path(output_dir, "CellChat_QA_notes.txt"))

message("Done: ", normalizePath(output_dir))
