#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(scales)
})
source("scripts/figure2_ghi_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript scripts/make_figure2G_strict.R <all.H5Seurat> <healthy_F.csv.gz> <query_F.csv.gz> <output_dir>")
}
h5_path <- normalizePath(args[[1]], mustWork = TRUE)
healthy_path <- normalizePath(args[[2]], mustWork = TRUE)
query_path <- normalizePath(args[[3]], mustWork = TRUE)
output_dir <- args[[4]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "supplementary"), recursive = TRUE, showWarnings = FALSE)
set.seed(264508)

ref_order <- c(
  "Homeostatic/adipogenic ASPC", "Papillary/conventional fibroblast",
  "ECM-remodeling fibroblast", "Pericyte", "VSMC"
)
candidate_labels <- c(
  "Homeostatic/adipogenic ASPC" = "ASPC-like",
  "Papillary/conventional fibroblast" = "Papillary/conventional FB-like",
  "ECM-remodeling fibroblast" = "ECM-remodeling FB-like",
  "Pericyte" = "Pericyte-like", "VSMC" = "VSMC-like"
)
marker_modules <- list(
  `General stromal` = c("DCN", "LUM", "COL1A1"),
  `ASPC-like` = c("PDGFRA", "CD34", "PI16", "DPP4", "CFD", "APOD"),
  `Papillary FB` = c("COL13A1", "APCDD1", "WIF1"),
  `ECM remodeling` = c("POSTN", "THBS2", "SFRP4", "FN1"),
  `Pericyte` = c("RGS5", "CSPG4", "PDGFRB"),
  `VSMC` = c("MYH11", "CNN1", "TAGLN"),
  `Other-lineage checks` = c("PTPRC", "PECAM1", "EPCAM")
)
marker_genes <- unique(unlist(marker_modules, use.names = FALSE))
module_annotation <- rbindlist(lapply(names(marker_modules), function(x) {
  data.table(module = x, gene = marker_modules[[x]])
}))

ref <- fread(healthy_path)
query <- fread(query_path)
md <- read_h5_metadata(h5_path)
md_keep <- md[, .(cell_id, source_row_h5 = source_row, donor_id = library_id,
                  sample_id = library_id, disease_group = fifelse(health == "LS", "LoSc", "Healthy"),
                  batch, activity, processing, sex, subtype,
                  nCount_RNA = as.numeric(nCount_RNA), nFeature_RNA = as.numeric(nFeature_RNA),
                  percent.mt = as.numeric(percent.mt), percent.ribo = as.numeric(percent.ribo))]
ref <- merge(ref, md_keep, by = "cell_id", all.x = TRUE, sort = FALSE)
query <- merge(query, md_keep, by = "cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(query$source_row_h5) || anyNA(ref$source_row_h5)) stop("Figure 2F cells do not map completely to H5Seurat")

# Query clustering uses all qualified LoSc stromal cells and the raw RNA count
# matrix. The graph is built in 30-dimensional PCA space, never from 2-D UMAP.
query_counts <- read_h5_counts(h5_path, query$source_row_h5)
saveRDS(query_counts, file.path(output_dir, "Figure2G_LoSc_query_raw_counts.rds"), compress = FALSE)
qmeta <- as.data.frame(query[, .(cell_id, donor_id, sample_id, disease_group, batch, onset,
                                 activity, processing, sex, subtype, nCount_RNA,
                                 nFeature_RNA, percent.mt, percent.ribo)])
rownames(qmeta) <- qmeta$cell_id
qobj <- CreateSeuratObject(counts = query_counts, meta.data = qmeta, min.cells = 0, min.features = 0)
qobj <- NormalizeData(qobj, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
qobj <- FindVariableFeatures(qobj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
qobj <- ScaleData(qobj, features = VariableFeatures(qobj), verbose = FALSE)
qobj <- RunPCA(qobj, features = VariableFeatures(qobj), npcs = 30, verbose = FALSE)
qobj <- FindNeighbors(qobj, reduction = "pca", dims = 1:30, k.param = 30,
                      graph.name = c("query_nn", "query_snn"), verbose = FALSE)
for (res in c(0.2, 0.4, 0.6)) {
  qobj <- FindClusters(qobj, graph.name = "query_snn", resolution = res,
                       algorithm = 1, random.seed = 264508,
                       cluster.name = paste0("query_res_", res), verbose = FALSE)
}
qobj <- RunUMAP(qobj, reduction = "pca", dims = 1:30, n.neighbors = 30,
                min.dist = 0.3, metric = "cosine", seed.use = 264508,
                reduction.name = "query_umap", verbose = FALSE)
qmeta_out <- qobj[[]]
qmd <- as.data.table(qmeta_out)
qmd[, cell_id := rownames(qmeta_out)]
qumap <- Embeddings(qobj, "query_umap")
qmd[, `:=`(query_UMAP_1 = qumap[cell_id, 1], query_UMAP_2 = qumap[cell_id, 2])]
qmd[, query_res_0.4 := as.character(query_res_0.4)]
cent <- qmd[, .(x = median(query_UMAP_1), y = median(query_UMAP_2)), by = query_res_0.4][order(x, y)]
cent[, query_cluster := paste0("Q", seq_len(.N))]
qmd <- merge(qmd, cent[, .(query_res_0.4, query_cluster)], by = "query_res_0.4", all.x = TRUE, sort = FALSE)
query <- merge(query, qmd[, .(cell_id, query_res_0.2 = `query_res_0.2`, query_res_0.4,
                              query_res_0.6 = `query_res_0.6`, query_cluster,
                              query_UMAP_1, query_UMAP_2)], by = "cell_id", all.x = TRUE, sort = FALSE)

# Raw-count donor-balanced dotplot values.
all_cells <- rbindlist(list(
  ref[, .(cell_id, source_row_h5, donor_id, nCount_RNA,
          expression_group = as.character(reference_state), reference_or_query = "Healthy reference")],
  query[, .(cell_id, source_row_h5, donor_id, nCount_RNA,
            expression_group = query_cluster, reference_or_query = "LoSc query")]
))
marker_counts <- read_h5_counts(h5_path, all_cells$source_row_h5, marker_genes)
gene_presence <- data.table(gene_requested = marker_genes,
                            gene_mapped = marker_genes %chin% rownames(marker_counts))
fwrite(gene_presence, file.path(output_dir, "Figure2G_gene_mapping_and_missing.csv"))
expr <- donor_balanced_expression(marker_counts, all_cells, "expression_group", marker_genes)

# Independent query-cluster identity is assigned from multi-gene similarity to
# the five Healthy reference profiles; Figure 2F mapping status is not used to
# form clusters or to select the winning identity.
identity_genes <- setdiff(marker_genes, c(marker_modules[["General stromal"]], marker_modules[["Other-lineage checks"]]))
profile <- dcast(expr$balanced[gene %chin% identity_genes], expression_group ~ gene,
                 value.var = "color_value", fill = 0)
profile_mat <- as.matrix(profile[, -1])
rownames(profile_mat) <- profile$expression_group
scale_by_gene <- scale(profile_mat)
scale_by_gene[!is.finite(scale_by_gene)] <- 0
ref_mat <- scale_by_gene[ref_order, , drop = FALSE]
q_names <- cent$query_cluster
q_mat <- scale_by_gene[q_names, , drop = FALSE]
corr <- t(vapply(seq_len(nrow(q_mat)), function(i) {
  vapply(seq_len(nrow(ref_mat)), function(j) suppressWarnings(cor(q_mat[i, ], ref_mat[j, ], method = "spearman")), numeric(1))
}, numeric(nrow(ref_mat))))
rownames(corr) <- q_names
colnames(corr) <- ref_order
corr[!is.finite(corr)] <- -1
identity <- rbindlist(lapply(q_names, function(qn) {
  x <- sort(corr[qn, ], decreasing = TRUE)
  best <- names(x)[1]
  margin <- unname(x[1] - x[2])
  confidence <- if (x[1] >= 0.60 && margin >= 0.10) "High" else if (x[1] >= 0.35 && margin >= 0.05) "Moderate" else "Low"
  data.table(query_cluster = qn, best_reference = best,
             identity_candidate = unname(candidate_labels[best]),
             identity_confidence = confidence,
             best_spearman = unname(x[1]), second_spearman = unname(x[2]), margin = margin)
}))

detect_wide <- dcast(expr$balanced[expression_group %chin% q_names], expression_group ~ gene,
                     value.var = "pct_detected_donor_balanced", fill = 0)
for (g in marker_genes) if (!g %in% names(detect_wide)) detect_wide[[g]] <- 0
other_genes <- marker_modules[["Other-lineage checks"]]
other_max <- apply(as.matrix(detect_wide[, .SD, .SDcols = other_genes]), 1, max)
detect_wide[, other_lineage_max_pct := other_max]
identity <- merge(identity, detect_wide[, .(query_cluster = expression_group, other_lineage_max_pct)],
                  by = "query_cluster", all.x = TRUE)
identity[other_lineage_max_pct >= 20, `:=`(
  identity_candidate = "Unresolved / possible non-stromal admixture",
  identity_confidence = "Low"
)]
identity[identity_confidence == "Low" & !grepl("^Unresolved", identity_candidate),
         identity_candidate := paste0("Provisional ", identity_candidate)]
identity[, annotation_reason := sprintf(
  "Highest multi-gene Healthy-profile Spearman r=%.2f (margin %.2f); other-lineage max detection %.1f%%",
  best_spearman, margin, other_lineage_max_pct
)]

query <- merge(query, identity[, .(query_cluster, identity_candidate, identity_confidence,
                                    annotation_reason)], by = "query_cluster", all.x = TRUE, sort = FALSE)
query[, display_group := paste0(query_cluster, " | ", identity_candidate)]
display_map <- unique(query[, .(expression_group = query_cluster, display_group)])
expr$balanced <- merge(expr$balanced, display_map, by = "expression_group", all.x = TRUE)
expr$balanced[is.na(display_group), display_group := expression_group]
expr$pooled <- merge(expr$pooled, display_map, by = "expression_group", all.x = TRUE)
expr$pooled[is.na(display_group), display_group := expression_group]

ref_audit <- ref[, .(
  reference_or_query = "Healthy reference", donor_id, sample_id,
  reference_label_original = as.character(reference_state),
  mapping_status_original = "Healthy reference", mapping_score_original = 1,
  deviation_flag_original = FALSE, query_cluster = NA_character_,
  identity_candidate = as.character(reference_state), identity_confidence = "Reference",
  annotation_reason = evidence, qc_status = "Retained", cell_id,
  batch, onset, activity, processing, sex, subtype, nCount_RNA, nFeature_RNA,
  percent.mt, percent.ribo
)]
query_audit <- query[, .(
  reference_or_query = "LoSc query", donor_id, sample_id,
  reference_label_original = predicted_reference_state,
  mapping_status_original = mapping_status,
  mapping_score_original = mapping_confidence,
  deviation_flag_original = mapping_status == "Disease-deviated",
  query_cluster, identity_candidate, identity_confidence, annotation_reason,
  qc_status = "Retained", cell_id, batch, onset, activity, processing, sex, subtype,
  nCount_RNA, nFeature_RNA, percent.mt, percent.ribo
)]
cell_audit <- rbindlist(list(ref_audit, query_audit), fill = TRUE)

cluster_summary <- query[, .(
  n_cells = .N, n_donors = uniqueN(donor_id), median_nCount_RNA = median(nCount_RNA),
  median_nFeature_RNA = median(nFeature_RNA), median_percent_mt = median(percent.mt),
  uncertain_n = sum(mapping_status == "Uncertain"),
  disease_deviated_n = sum(mapping_status == "Disease-deviated"),
  uncertain_percent = 100 * mean(mapping_status == "Uncertain"),
  disease_deviated_percent = 100 * mean(mapping_status == "Disease-deviated"),
  identity_candidate = first(identity_candidate), identity_confidence = first(identity_confidence),
  annotation_reason = first(annotation_reason)
), by = query_cluster][order(as.integer(sub("Q", "", query_cluster)))]
donor_cluster_counts <- query[, .N, by = .(donor_id, query_cluster, identity_candidate, mapping_status)]
mapping_audit <- rbindlist(list(
  ref[, .(reference_or_query = "Healthy reference", mapping_status = "Healthy reference",
          reference_label = as.character(reference_state), n_cells = .N, n_donors = uniqueN(donor_id)),
      by = reference_state][, reference_state := NULL],
  query[, .(n_cells = .N, n_donors = uniqueN(donor_id)),
        by = .(mapping_status, reference_label = predicted_reference_state)][
          , reference_or_query := "LoSc query"]
), fill = TRUE)

fwrite(mapping_audit, file.path(output_dir, "Fig2F_mapping_audit.csv"))
fwrite(cell_audit, file.path(output_dir, "Figure2G_cell_metadata.csv.gz"))
fwrite(cluster_summary, file.path(output_dir, "Figure2G_query_cluster_summary.csv"))
fwrite(identity, file.path(output_dir, "Figure2G_identity_evidence.csv"))
fwrite(donor_cluster_counts, file.path(output_dir, "Figure2G_donor_cluster_counts.csv"))
fwrite(expr$by_donor, file.path(output_dir, "Figure2G_expression_by_donor.csv.gz"))
fwrite(expr$balanced, file.path(output_dir, "Figure2G_dotplot_values_donor_balanced.csv"))
fwrite(expr$pooled, file.path(output_dir, "Figure2G_dotplot_values_pooled_cells.csv"))
fwrite(query[grepl("ASPC-like", identity_candidate), .(
  cell_id, donor_id, sample_id, query_cluster, identity_candidate, identity_confidence,
  mapping_status_original = mapping_status, predicted_reference_state,
  mapping_confidence, mapping_entropy, mean_reference_distance
)], file.path(output_dir, "Figure2G_ASPC_candidate_cells.csv.gz"))
fwrite(module_annotation, file.path(output_dir, "Figure2G_marker_panel.csv"))

# Supplementary query UMAP.
cluster_palette <- setNames(hue_pal(h = c(15, 375), c = 100, l = 55)(length(q_names)), q_names)
label_pos <- query[, .(query_UMAP_1 = median(query_UMAP_1), query_UMAP_2 = median(query_UMAP_2)), by = query_cluster]
p_umap <- ggplot(query, aes(query_UMAP_1, query_UMAP_2, color = query_cluster)) +
  geom_point(size = 0.22, alpha = 0.9, stroke = 0) +
  geom_text(data = label_pos, aes(label = query_cluster), color = "black", fontface = "bold",
            family = "Helvetica", size = 2.7, show.legend = FALSE) +
  scale_color_manual(values = cluster_palette) + coord_equal() +
  labs(title = "LoSc stromal query clusters", x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_family = "Helvetica", base_size = 7) +
  theme(legend.position = "none", axis.line = element_line(linewidth = 0.35),
        axis.ticks = element_blank(), axis.text = element_blank(),
        plot.title = element_text(face = "bold", size = 8, hjust = 0.5))
save_nature_plot(p_umap, file.path(output_dir, "supplementary/Figure2G_query_cluster_UMAP"), 92, 82)

# Main donor-balanced dotplot. Reference and query blocks share the same genes.
row_ref <- rev(ref_order)
row_query <- rev(cluster_summary$query_cluster)
row_lookup <- c(setNames(ref_order, ref_order), setNames(unique(query$display_group), unique(query$query_cluster)))
plot_dt <- copy(expr$balanced)
plot_dt[, row_label := row_lookup[expression_group]]
row_levels <- c(row_lookup[row_query], row_ref)
plot_dt[, row_label := factor(row_label, levels = row_levels)]
plot_dt[, gene := factor(gene, levels = marker_genes)]
plot_dt[, block := fifelse(expression_group %chin% ref_order, "Healthy reference", "LoSc query clusters")]
cap <- quantile(plot_dt$color_value, 0.98, na.rm = TRUE)
plot_dt[, color_plot := pmin(color_value, cap)]

p_dot <- ggplot(plot_dt, aes(gene, row_label)) +
  geom_point(aes(size = pct_detected_donor_balanced, fill = color_plot), shape = 21,
             color = "#242424", stroke = 0.18) +
  annotate("segment", x = 0.5, xend = length(marker_genes) + 0.5,
           y = length(row_query) + 0.5, yend = length(row_query) + 0.5,
           linewidth = 0.55, color = "#333333") +
  scale_size_continuous(range = c(0.15, 3.1), limits = c(0, 100),
                        breaks = c(20, 50, 80), name = "Detected cells (%)") +
  scale_fill_stepsn(colors = c("#FFFFFF", "#F4D2C3", "#D4775B", "#8D251E"),
                    limits = c(0, cap), n.breaks = 6, oob = squish,
                    name = "Mean expression\nlog1p(donor mean)") +
  scale_x_discrete(position = "top", expand = expansion(add = 0.35)) +
  scale_y_discrete(expand = expansion(add = 0.35)) +
  labs(title = "Marker expression supporting stromal identities and annotation of reference-mapping outliers",
       x = NULL, y = NULL) +
  theme_classic(base_family = "Helvetica", base_size = 6.5) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(color = "#222222", fill = NA, linewidth = 0.45),
    axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 60, hjust = 0, vjust = 0, size = 5.7, color = "#111111"),
    axis.text.y = element_text(size = 5.9, color = "#111111"),
    plot.title = element_text(face = "bold", size = 8.2, hjust = 0, margin = margin(b = 4)),
    legend.position = "right", legend.title = element_text(size = 6.2),
    legend.text = element_text(size = 5.8), legend.key.height = unit(3.0, "mm"),
    plot.margin = margin(4, 4, 4, 4)
  )

side_dt <- rbindlist(list(
  ref[, .(display_group = as.character(reference_state), n_cells = .N, n_donors = uniqueN(donor_id),
          confidence = "Reference", uncertain_percent = 0, disease_deviated_percent = 0), by = reference_state][, reference_state := NULL],
  cluster_summary[, .(display_group = row_lookup[query_cluster], n_cells, n_donors,
                      confidence = identity_confidence, uncertain_percent, disease_deviated_percent)]
))
side_dt[, display_group := factor(display_group, levels = row_levels)]
side_long_num <- melt(side_dt, id.vars = "display_group",
                      measure.vars = c("n_cells", "n_donors", "uncertain_percent", "disease_deviated_percent"),
                      variable.name = "metric", value.name = "value")
side_long_num[, label := fifelse(metric %chin% c("uncertain_percent", "disease_deviated_percent"),
                                 sprintf("%.1f", value), as.character(round(value)))]
side_long <- rbindlist(list(
  side_long_num[, .(display_group, metric = as.character(metric), label)],
  side_dt[, .(display_group, metric = "confidence",
              label = fifelse(confidence == "Reference", "Ref.",
                              fifelse(confidence == "Moderate", "Mod.", confidence)))]
))
side_long[, metric := factor(metric,
  levels = c("n_cells", "n_donors", "confidence", "uncertain_percent", "disease_deviated_percent"),
  labels = c("Cells", "Donors", "Conf.", "Uncertain (%)", "Deviated (%)"))]
p_side <- ggplot(side_long, aes(metric, display_group, label = label)) +
  geom_tile(fill = "white", color = "#BDBDBD", linewidth = 0.25) +
  geom_text(family = "Helvetica", size = 1.8, color = "#111111") +
  scale_x_discrete(position = "top", expand = c(0, 0)) + scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL) + theme_void(base_family = "Helvetica") +
  theme(axis.text.x = element_text(size = 5.2, angle = 60, hjust = 0, vjust = 0, color = "#111111"),
        plot.margin = margin(23, 3, 4, 0))

p_main <- p_dot + p_side + plot_layout(widths = c(5.1, 1.2), guides = "collect") &
  theme(legend.position = "right")
save_nature_plot(p_main, file.path(output_dir, "Figure2G_stromal_identity_dotplot_strict"), 183, 112)

writeLines(c(
  "panel=Figure 2G strict rebuild",
  "query_clustering=all qualified LoSc stromal cells; raw RNA counts; LogNormalize; 2000 HVGs; PCA30; k30 SNN; Louvain resolution 0.4",
  "sensitivity=resolutions 0.2 and 0.6 retained in cell metadata",
  "identity=independent multi-gene similarity to five Healthy donor-balanced reference profiles",
  "dotplot_color=log1p(equal-weight mean of donor-level mean 10k-normalized counts)",
  "dotplot_size=equal-weight mean of donor-level percent detected",
  "donor_absent_from_group=missing, never zero",
  "status_guard=Uncertain and Disease-deviated remain mapping-status fields, not cell identities",
  sprintf("healthy_reference_cells=%d", nrow(ref)),
  sprintf("losc_query_cells=%d", nrow(query)),
  sprintf("query_clusters=%d", uniqueN(query$query_cluster)),
  sprintf("aspc_candidate_cells=%d", query[grepl("ASPC-like", identity_candidate), .N])
), file.path(output_dir, "Figure2G_QA_notes.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Figure2G_sessionInfo.txt"))
message("Figure 2G strict rebuild complete")
