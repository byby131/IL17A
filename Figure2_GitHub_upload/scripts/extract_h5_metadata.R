#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript scripts/extract_h5_metadata.R <file.H5Seurat> <output_dir>")
}

h5_path <- normalizePath(args[[1]], mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_h5_column <- function(node) {
  if (inherits(node, "H5D")) return(node$read())
  if (inherits(node, "H5Group") && all(c("levels", "values") %in% names(node))) {
    lev <- node[["levels"]]$read()
    val <- as.integer(node[["values"]]$read())
    # H5Seurat stores R factor codes (1-based); retain missing values if present.
    out <- rep(NA_character_, length(val))
    ok <- !is.na(val) & val >= 1L & val <= length(lev)
    out[ok] <- lev[val[ok]]
    return(factor(out, levels = lev))
  }
  stop("Unsupported H5 metadata node")
}

h5 <- H5File$new(h5_path, "r")
on.exit(h5$close_all())

cell_ids <- h5[["cell.names"]]$read()
md_group <- h5[["meta.data"]]
md_names <- setdiff(names(md_group), "_index")
md <- data.table(cell_id = cell_ids)
for (nm in md_names) md[[nm]] <- read_h5_column(md_group[[nm]])

umap_matrix <- h5[["reductions"]][["umap"]][["cell.embeddings"]]$read()
if (nrow(umap_matrix) != nrow(md) || ncol(umap_matrix) != 2L) {
  stop("Unexpected UMAP dimensions")
}
md[, `:=`(UMAP_1 = umap_matrix[, 1], UMAP_2 = umap_matrix[, 2])]

fwrite(md, file.path(output_dir, "metadata_umap.csv.gz"))

fields <- intersect(
  c("celltype", "subtype", "health", "onset", "activity", "processing", "batch", "sex", "library_id", "seurat_clusters"),
  names(md)
)
for (field in fields) {
  tab <- md[, .(n_cells = .N), by = field][order(-n_cells)]
  fwrite(tab, file.path(output_dir, paste0("table_", field, ".csv")))
  cat("\n## ", field, "\n", sep = "")
  print(tab)
}

for (spec in list(
  celltype_by_health = c("celltype", "health"),
  subtype_by_health = c("subtype", "health"),
  donor_by_health = c("library_id", "health"),
  donor_by_onset = c("library_id", "onset"),
  celltype_by_donor = c("celltype", "library_id")
)) {
  nm <- names(spec)
}

cross_specs <- list(
  celltype_by_health = c("celltype", "health"),
  subtype_by_health = c("subtype", "health"),
  donor_by_health = c("library_id", "health"),
  donor_by_onset = c("library_id", "onset"),
  celltype_by_donor = c("celltype", "library_id")
)
for (nm in names(cross_specs)) {
  fields_i <- cross_specs[[nm]]
  if (all(fields_i %in% names(md))) {
    tab <- md[, .(n_cells = .N), by = fields_i]
    fwrite(tab, file.path(output_dir, paste0(nm, ".csv")))
  }
}

cat("\nmetadata_rows:", nrow(md), "\n")
cat("metadata_columns:", paste(names(md), collapse = ", "), "\n")
