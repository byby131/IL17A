suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
  library(Matrix)
})

read_h5_column <- function(node) {
  if (inherits(node, "H5D")) return(node$read())
  if (inherits(node, "H5Group") && all(c("levels", "values") %in% names(node))) {
    lev <- node[["levels"]]$read()
    val <- as.integer(node[["values"]]$read())
    out <- rep(NA_character_, length(val))
    ok <- !is.na(val) & val >= 1L & val <= length(lev)
    out[ok] <- lev[val[ok]]
    return(out)
  }
  stop("Unsupported H5 metadata node")
}

read_h5_metadata <- function(h5_path) {
  h5 <- H5File$new(h5_path, "r")
  on.exit(h5$close_all())
  ids <- h5[["cell.names"]]$read()
  grp <- h5[["meta.data"]]
  fields <- setdiff(names(grp), "_index")
  out <- data.table(cell_id = ids, source_row = seq_along(ids))
  for (nm in fields) out[[nm]] <- read_h5_column(grp[[nm]])
  out
}

# H5Seurat stores RNA matrices in CSC form with one compressed column per cell.
# This reader scans contiguous cell blocks and returns only requested source rows
# and, optionally, requested genes without materialising the full 60k-cell matrix.
read_h5_counts <- function(h5_path, source_rows, genes = NULL, chunk_cells = 1000L) {
  source_rows <- sort(unique(as.integer(source_rows)))
  if (!length(source_rows)) stop("No source rows supplied")
  h5 <- H5File$new(h5_path, "r")
  on.exit(h5$close_all())
  rna <- h5[["assays"]][["RNA"]]
  features_all <- rna[["features"]]$read()
  cell_ids_all <- h5[["cell.names"]]$read()
  if (is.null(genes)) {
    feature_keep <- seq_along(features_all)
  } else {
    feature_keep <- match(unique(genes), features_all)
    feature_keep <- feature_keep[!is.na(feature_keep)]
  }
  features <- features_all[feature_keep]
  feature_map <- integer(length(features_all))
  feature_map[feature_keep] <- seq_along(feature_keep)
  col_map <- integer(length(cell_ids_all))
  col_map[source_rows] <- seq_along(source_rows)

  grp <- rna[["counts"]]
  indptr <- as.numeric(grp[["indptr"]]$read())
  indices_ds <- grp[["indices"]]
  data_ds <- grp[["data"]]
  i_list <- vector("list", ceiling(length(cell_ids_all) / chunk_cells))
  j_list <- vector("list", length(i_list))
  x_list <- vector("list", length(i_list))
  out_index <- 0L

  for (cell_start in seq.int(1L, length(cell_ids_all), by = chunk_cells)) {
    cell_end <- min(length(cell_ids_all), cell_start + chunk_cells - 1L)
    selected_here <- source_rows[source_rows >= cell_start & source_rows <= cell_end]
    if (!length(selected_here)) next
    ptr <- indptr[cell_start:(cell_end + 1L)]
    n_by_cell <- diff(ptr)
    nz_start <- ptr[[1]] + 1
    nz_end <- ptr[[length(ptr)]]
    if (nz_end < nz_start) next
    feature_index_0 <- as.integer(indices_ds[nz_start:nz_end])
    value <- as.numeric(data_ds[nz_start:nz_end])
    local_cell <- rep.int(seq_len(cell_end - cell_start + 1L), n_by_cell)
    global_cell <- cell_start - 1L + local_cell
    row_new <- feature_map[feature_index_0 + 1L]
    col_new <- col_map[global_cell]
    keep <- row_new > 0L & col_new > 0L
    if (!any(keep)) next
    out_index <- out_index + 1L
    i_list[[out_index]] <- row_new[keep]
    j_list[[out_index]] <- col_new[keep]
    x_list[[out_index]] <- value[keep]
  }
  i_list <- i_list[seq_len(out_index)]
  j_list <- j_list[seq_len(out_index)]
  x_list <- x_list[seq_len(out_index)]
  sparseMatrix(
    i = unlist(i_list, use.names = FALSE),
    j = unlist(j_list, use.names = FALSE),
    x = unlist(x_list, use.names = FALSE),
    dims = c(length(features), length(source_rows)),
    dimnames = list(features, cell_ids_all[source_rows])
  )
}

save_nature_plot <- function(plot, stem, width_mm, height_mm, dpi = 600) {
  dir.create(dirname(stem), recursive = TRUE, showWarnings = FALSE)
  ggplot2::ggsave(paste0(stem, ".pdf"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = grDevices::pdf, bg = "white",
                  family = "Helvetica", useDingbats = FALSE)
  ggplot2::ggsave(paste0(stem, ".svg"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = svglite::svglite, bg = "white")
  ggplot2::ggsave(paste0(stem, ".tiff"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = ragg::agg_tiff, dpi = dpi,
                  compression = "lzw", bg = "white")
  ggplot2::ggsave(paste0(stem, "_preview.png"), plot, width = width_mm, height = height_mm,
                  units = "mm", device = ragg::agg_png, dpi = 300, bg = "white")
}

donor_balanced_expression <- function(counts, cell_meta, row_field, genes) {
  genes <- intersect(genes, rownames(counts))
  cells <- intersect(colnames(counts), cell_meta$cell_id)
  cm <- copy(cell_meta[match(cells, cell_id)])
  stopifnot(identical(cm$cell_id, cells))
  totals <- as.numeric(cm$nCount_RNA)
  totals[!is.finite(totals) | totals <= 0] <- as.numeric(Matrix::colSums(counts[, cells, drop = FALSE]))[!is.finite(totals) | totals <= 0]
  norm <- counts[genes, cells, drop = FALSE] %*% Diagonal(x = 10000 / totals)
  pooled <- rbindlist(lapply(seq_along(genes), function(k) {
    data.table(
      cell_id = cells,
      donor_id = cm$donor_id,
      expression_group = cm[[row_field]],
      gene = genes[[k]],
      norm_count = as.numeric(norm[k, ]),
      detected = as.numeric(counts[genes[[k]], cells] > 0)
    )
  }))
  by_donor <- pooled[, .(
    mean_expr = mean(norm_count),
    pct_detected = 100 * mean(detected)
  ), by = .(donor_id, expression_group, gene)]
  balanced <- by_donor[, .(
    n_donors = uniqueN(donor_id),
    mean_expr_donor_balanced = mean(mean_expr),
    pct_detected_donor_balanced = mean(pct_detected),
    color_value = log1p(mean(mean_expr))
  ), by = .(expression_group, gene)]
  pooled_group <- pooled[, .(
    n_cells = .N,
    mean_expr_pooled = mean(norm_count),
    pct_detected_pooled = 100 * mean(detected),
    color_value_pooled = log1p(mean(norm_count))
  ), by = .(expression_group, gene)]
  list(by_donor = by_donor, balanced = merge(balanced, pooled_group,
                                             by = c("expression_group", "gene"), all = TRUE),
       pooled = pooled_group)
}
