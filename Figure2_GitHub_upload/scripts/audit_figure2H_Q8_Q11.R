#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
source("scripts/figure2_ghi_common.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/audit_figure2H_Q8_Q11.R <all.H5Seurat> <Figure2G_cell_metadata.csv.gz> <output_dir>")
}
h5_path <- normalizePath(args[[1]], mustWork = TRUE)
gmeta_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

gmeta <- fread(gmeta_path)
gmeta[, analysis_group := fifelse(
  reference_or_query == "Healthy reference" &
    reference_label_original == "Homeostatic/adipogenic ASPC",
  "Healthy ASPC-like",
  fifelse(reference_or_query == "LoSc query" & query_cluster == "Q8",
          "LoSc Q8", fifelse(reference_or_query == "LoSc query" & query_cluster == "Q11",
                             "LoSc Q11", NA_character_))
)]
cells <- gmeta[!is.na(analysis_group) & qc_status == "Retained"]

md <- read_h5_metadata(h5_path)
cells <- merge(cells, md[, .(cell_id, source_row_h5 = source_row)], by = "cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(cells$source_row_h5)) stop("Audit cells do not map completely to H5 source rows")
counts <- read_h5_counts(h5_path, cells$source_row_h5)
cells <- cells[match(colnames(counts), cell_id)]
stopifnot(identical(cells$cell_id, colnames(counts)))

cell_totals <- Matrix::colSums(counts)
cell_detected <- Matrix::colSums(counts > 0)
cells[, raw_total_counts := as.numeric(cell_totals)]
cells[, raw_detected_genes := as.integer(cell_detected)]

qc <- rbindlist(lapply(split(seq_len(nrow(cells)), interaction(cells$analysis_group, cells$donor_id,
                                                               drop = TRUE)), function(idx) {
  cm <- cells[idx]
  mat <- counts[, idx, drop = FALSE]
  totals <- Matrix::colSums(mat)
  data.table(
    analysis_group = cm$analysis_group[[1]],
    donor_id = cm$donor_id[[1]],
    batch = cm$batch[[1]],
    n_cells = nrow(cm),
    pseudobulk_total_counts = sum(totals),
    pseudobulk_detected_genes = sum(Matrix::rowSums(mat) > 0),
    max_single_cell_count_fraction = max(totals) / sum(totals),
    median_cell_counts = median(totals),
    median_cell_detected_genes = median(Matrix::colSums(mat > 0)),
    identity_confidence = paste(sort(unique(cm$identity_confidence)), collapse = ";"),
    mapping_statuses = paste(sort(unique(cm$mapping_status_original)), collapse = ";")
  )
}))

# These information thresholds are frozen before DE/GSEA inspection. They are
# used only to exclude donor-by-group pseudobulks that are too sparse or
# effectively represent a single cell; they are not tuned to pathway results.
qc[, pass_n_cells := n_cells >= 3L]
qc[, pass_total_counts := pseudobulk_total_counts >= 5000]
qc[, pass_detected_genes := pseudobulk_detected_genes >= 1000L]
qc[, pass_cell_dominance := max_single_cell_count_fraction <= 0.80]
qc[, expression_qc_pass := pass_n_cells & pass_total_counts & pass_detected_genes & pass_cell_dominance]
qc[, exclusion_reason := fifelse(expression_qc_pass, "Pass",
  paste(
    ifelse(!pass_n_cells, "n_cells<3", ""),
    ifelse(!pass_total_counts, "total_counts<5000", ""),
    ifelse(!pass_detected_genes, "detected_genes<1000", ""),
    ifelse(!pass_cell_dominance, "single_cell_fraction>0.80", ""),
    sep = ";"
  )
)]
qc[, exclusion_reason := gsub("(^;+|;+$)", "", gsub(";+", ";", exclusion_reason))]

coverage <- qc[, .(
  donors_nonempty = .N,
  donors_qc_pass = sum(expression_qc_pass),
  cells_total = sum(n_cells),
  cells_qc_pass_donors = sum(n_cells[expression_qc_pass]),
  batches_nonempty = paste(sort(unique(batch)), collapse = ","),
  batches_qc_pass = paste(sort(unique(batch[expression_qc_pass])), collapse = ",")
), by = analysis_group]

fwrite(cells, file.path(output_dir, "Figure2H_Q8_Q11_frozen_cell_metadata.csv.gz"))
fwrite(qc[order(analysis_group, donor_id)], file.path(output_dir, "Figure2H_donor_group_coverage_and_quality.csv"))
fwrite(coverage, file.path(output_dir, "Figure2H_group_coverage_summary.csv"))

cat("Coverage summary\n")
print(coverage)
cat("\nQualified donor-by-group pseudobulks\n")
print(qc[expression_qc_pass == TRUE, .(analysis_group, donor_id, batch, n_cells,
                              pseudobulk_total_counts, pseudobulk_detected_genes,
                              max_single_cell_count_fraction)][order(analysis_group, donor_id)])
