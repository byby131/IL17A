#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
  library(Matrix)
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/cluster_and_profile_tcells.R <all.H5Seurat> <metadata_umap.csv.gz> <output_dir>")
}

h5_path <- normalizePath(args[[1]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

set.seed(264508)
md <- fread(metadata_path)
md[, seurat_clusters := as.character(seurat_clusters)]
t_keep <- md$celltype == "T Cells" | md$seurat_clusters %chin% c("9", "20")
t_md <- copy(md[t_keep])
t_md[, source_row := which(t_keep)]

if (nrow(t_md) < 100L) stop("Too few T cells selected")

h5 <- H5File$new(h5_path, "r")
on.exit(h5$close_all())

harmony <- h5[["reductions"]][["harmony"]][["cell.embeddings"]]$read()
if (nrow(harmony) != nrow(md)) stop("Harmony embedding and metadata cell counts differ")
harmony_t <- harmony[t_keep, 1:30, drop = FALSE]
rownames(harmony_t) <- t_md$cell_id
colnames(harmony_t) <- paste0("HARMONY_", seq_len(ncol(harmony_t)))

dummy <- sparseMatrix(
  i = rep.int(1L, nrow(t_md)),
  j = seq_len(nrow(t_md)),
  x = rep.int(1, nrow(t_md)),
  dims = c(1L, nrow(t_md)),
  dimnames = list("placeholder", t_md$cell_id)
)
meta_df <- as.data.frame(t_md)
rownames(meta_df) <- meta_df$cell_id

obj <- CreateSeuratObject(counts = dummy, meta.data = meta_df, min.cells = 0, min.features = 0)
obj[["harmony"]] <- CreateDimReducObject(embeddings = harmony_t, key = "HARMONY_", assay = "RNA")
obj <- FindNeighbors(
  obj, reduction = "harmony", dims = 1:30,
  k.param = 25, prune.SNN = 1 / 15,
  graph.name = c("tcell_nn", "tcell_snn"), verbose = FALSE
)

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)
for (resolution in resolutions) {
  obj <- FindClusters(
    obj, graph.name = "tcell_snn", resolution = resolution,
    algorithm = 1, random.seed = 264508,
    cluster.name = paste0("tcell_res_", resolution), verbose = FALSE
  )
}

obj <- RunUMAP(
  obj, reduction = "harmony", dims = 1:30,
  n.neighbors = 25, min.dist = 0.28, metric = "cosine",
  seed.use = 264508, reduction.name = "tcell_umap",
  reduction.key = "tUMAP_", verbose = FALSE
)

out_md <- as.data.table(obj[[]], keep.rownames = "seurat_cell_id")
umap <- as.data.table(Embeddings(obj, "tcell_umap"), keep.rownames = "seurat_cell_id")
setnames(umap, c("seurat_cell_id", "tcell_UMAP_1", "tcell_UMAP_2"))
out_md <- merge(out_md, umap, by = "seurat_cell_id", all.x = TRUE, sort = FALSE)
out_md[, cell_id := seurat_cell_id]
out_md[, seurat_cell_id := NULL]
for (resolution in resolutions) {
  field <- paste0("tcell_res_", resolution)
  out_md[[field]] <- as.character(out_md[[field]])
}

marker_modules <- list(
  Pan_T = c("CD3D", "CD3E", "TRBC1", "LTB"),
  CD4_naive = c("CCR7", "LEF1", "TCF7", "SELL", "MAL", "LTB"),
  CD4_memory = c("IL7R", "LTB", "KLRB1", "CCR6", "RORA"),
  Th17_like = c("RORC", "RORA", "CCR6", "KLRB1", "CCL20", "IL23R", "IL17A", "IL17F"),
  Treg = c("FOXP3", "IL2RA", "CTLA4", "TIGIT", "IKZF2", "TNFRSF4"),
  Activated = c("CD69", "FOS", "JUN", "JUNB", "ICOS", "TNFRSF9", "HLA-DRA"),
  Cycling = c("MKI67", "TOP2A", "STMN1", "TUBA1B"),
  CD8_cytotoxic = c("CD8A", "CD8B", "CCL5", "NKG7", "GNLY", "PRF1"),
  Gamma_delta = c("TRDC", "TRGC1", "TRGC2"),
  NK_like = c("NKG7", "KLRD1", "GNLY", "PRF1")
)
marker_annotation <- unique(rbindlist(lapply(names(marker_modules), function(module) {
  data.table(module = module, marker = marker_modules[[module]])
})))

rna <- h5[["assays"]][["RNA"]]
features <- rna[["features"]]$read()
marker_annotation[, feature_index_0 := match(marker, features) - 1L]
marker_annotation[, present := !is.na(feature_index_0)]
fwrite(marker_annotation, file.path(output_dir, "Tcell_marker_presence.csv"))

present_markers <- unique(marker_annotation[present == TRUE, marker])
present_index <- match(present_markers, features) - 1L
names(present_index) <- present_markers

expr_dense <- matrix(
  0, nrow = nrow(t_md), ncol = length(present_markers),
  dimnames = list(t_md$cell_id, present_markers)
)
expr <- rna[["data"]]
indptr <- as.numeric(expr[["indptr"]]$read())
indices_ds <- expr[["indices"]]
data_ds <- expr[["data"]]

local_row_map <- integer(nrow(md))
local_row_map[t_keep] <- seq_len(nrow(t_md))
marker_col_map <- match(seq_along(features) - 1L, present_index)
chunk_cells <- 1000L

for (cell_start in seq.int(1L, nrow(md), by = chunk_cells)) {
  cell_end <- min(nrow(md), cell_start + chunk_cells - 1L)
  ptr <- indptr[cell_start:(cell_end + 1L)]
  n_by_cell <- diff(ptr)
  nz_start_1 <- ptr[[1]] + 1
  nz_end_1 <- ptr[[length(ptr)]]

  if (nz_end_1 >= nz_start_1) {
    feature_index_0 <- as.integer(indices_ds[nz_start_1:nz_end_1])
    value <- as.numeric(data_ds[nz_start_1:nz_end_1])
    cell_local <- rep.int(seq_len(cell_end - cell_start + 1L), n_by_cell)
    source_rows <- cell_start - 1L + cell_local
    marker_columns <- marker_col_map[feature_index_0 + 1L]
    keep <- t_keep[source_rows] & !is.na(marker_columns)

    if (any(keep)) {
      target_rows <- local_row_map[source_rows[keep]]
      expr_dense[cbind(target_rows, marker_columns[keep])] <- value[keep]
    }
  }

  if (cell_start == 1L || cell_end %% 10000L == 0L || cell_end == nrow(md)) {
    message("Scanned cells ", cell_start, "-", cell_end, " / ", nrow(md))
  }
}

expr_dt <- as.data.table(expr_dense, keep.rownames = "cell_id")
out_md <- merge(out_md, expr_dt, by = "cell_id", all.x = TRUE, sort = FALSE)

score_gene <- function(gene) {
  values <- out_md[[gene]]
  sd_value <- sd(values)
  if (is.na(sd_value) || sd_value == 0) return(rep(0, length(values)))
  pmax(-3, pmin(3, as.numeric(scale(values))))
}

module_definitions <- rbindlist(lapply(names(marker_modules), function(module) {
  genes <- intersect(marker_modules[[module]], present_markers)
  data.table(module = module, marker = marker_modules[[module]], present = marker_modules[[module]] %in% present_markers,
             used_for_score = marker_modules[[module]] %in% genes)
}))
fwrite(module_definitions, file.path(output_dir, "Tcell_module_definitions.csv"))

for (module in names(marker_modules)) {
  genes <- intersect(marker_modules[[module]], present_markers)
  score_matrix <- do.call(cbind, lapply(genes, score_gene))
  out_md[[paste0(module, "_score")]] <- rowMeans(score_matrix)
}

fwrite(out_md, file.path(output_dir, "Tcell_secondary_clustering_and_expression.csv.gz"))

resolution_summary <- rbindlist(lapply(resolutions, function(resolution) {
  field <- paste0("tcell_res_", resolution)
  out_md[, .(
    n_cells = .N,
    n_donors = uniqueN(library_id),
    healthy_cells = sum(health == "Healthy"),
    losc_cells = sum(health != "Healthy"),
    median_cells_per_donor = median(as.numeric(table(library_id))),
    CD4_naive_score = mean(CD4_naive_score),
    CD4_memory_score = mean(CD4_memory_score),
    Th17_like_score = mean(Th17_like_score),
    Treg_score = mean(Treg_score),
    Activated_score = mean(Activated_score),
    Cycling_score = mean(Cycling_score),
    CD8_cytotoxic_score = mean(CD8_cytotoxic_score),
    Gamma_delta_score = mean(Gamma_delta_score),
    NK_like_score = mean(NK_like_score)
  ), by = .(cluster = get(field))][, resolution := resolution]
}))
setcolorder(resolution_summary, c("resolution", "cluster", setdiff(names(resolution_summary), c("resolution", "cluster"))))
fwrite(resolution_summary, file.path(output_dir, "Tcell_resolution_module_summary.csv"))

marker_summary <- rbindlist(lapply(resolutions, function(resolution) {
  field <- paste0("tcell_res_", resolution)
  rbindlist(lapply(present_markers, function(marker) {
    out_md[, .(
      n_total = .N,
      n_expressing = sum(get(marker) > 0),
      pct_exp = 100 * mean(get(marker) > 0),
      avg_exp = mean(get(marker))
    ), by = .(cluster = get(field))][, `:=`(resolution = resolution, marker = marker)]
  }))
}))
marker_summary[, avg_exp_scaled := {
  scaled <- as.numeric(scale(avg_exp))
  if (all(is.na(scaled))) rep(0, .N) else pmax(-2.5, pmin(2.5, scaled))
}, by = .(resolution, marker)]
marker_summary <- merge(marker_summary, unique(marker_annotation[, .(marker, module)]), by = "marker", allow.cartesian = TRUE)
setcolorder(marker_summary, c("resolution", "cluster", "module", "marker", "n_total", "n_expressing", "pct_exp", "avg_exp", "avg_exp_scaled"))
fwrite(marker_summary, file.path(output_dir, "Tcell_cluster_marker_summary.csv"))

donor_cluster <- rbindlist(lapply(resolutions, function(resolution) {
  field <- paste0("tcell_res_", resolution)
  part <- out_md[, .N, by = .(library_id, health, onset, cluster = get(field))]
  part[, resolution := resolution]
  setcolorder(part, c("resolution", "library_id", "health", "onset", "cluster", "N"))
  part
}))
fwrite(donor_cluster, file.path(output_dir, "Tcell_donor_cluster_counts.csv"))

provenance <- c(
  sprintf("input_cells=%s", nrow(md)),
  sprintf("selected_T_cells=%s", nrow(out_md)),
  sprintf("selected_donors=%s", uniqueN(out_md$library_id)),
  "selection=published T Cells annotation or global clusters 9/20",
  "embedding=author Harmony dimensions 1:30",
  "neighbors=k25; prune.SNN=1/15",
  "resolutions=0.2,0.4,0.6,0.8,1.0",
  "umap=n.neighbors25; min.dist0.28; cosine; seed264508",
  "expression=processed RNA normalized data layer",
  "missing_required_markers=IL23R,IL17A,IL17F absent from processed 17,093-feature matrix",
  "integrity=no T cells removed before CD4/CD8/gamma-delta annotation"
)
writeLines(provenance, file.path(output_dir, "Tcell_clustering_provenance.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "Tcell_clustering_sessionInfo.txt"))

message("T-cell clustering and marker profiling complete: ", nrow(out_md), " cells")
