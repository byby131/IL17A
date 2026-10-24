#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(hdf5r)
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(uwot)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript scripts/prepare_healthy_stromal_reference.R <all.H5Seurat> <metadata_umap.csv.gz> <output_dir>")
}

h5_path <- normalizePath(args[[1]], mustWork = TRUE)
metadata_path <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(264508)

md <- fread(metadata_path)
md[, seurat_clusters := as.character(seurat_clusters)]
stromal_keep <- md$celltype %chin% c("Fibroblasts", "Pericytes", "Smooth Muscle Cells") &
  md$seurat_clusters != "32"
st <- copy(md[stromal_keep])
st[, source_row := which(stromal_keep)]
st[, reference_status := fifelse(health == "Healthy", "Healthy reference", "LoSc query")]

h5 <- H5File$new(h5_path, "r")
on.exit(h5$close_all())
harmony_all <- h5[["reductions"]][["harmony"]][["cell.embeddings"]]$read()
if (nrow(harmony_all) != nrow(md)) stop("Harmony embedding and metadata differ")
harmony_st <- harmony_all[stromal_keep, 1:30, drop = FALSE]
rownames(harmony_st) <- st$cell_id
colnames(harmony_st) <- paste0("HARMONY_", seq_len(ncol(harmony_st)))

ref_keep <- st$health == "Healthy"
ref_harmony <- harmony_st[ref_keep, , drop = FALSE]
query_harmony <- harmony_st[!ref_keep, , drop = FALSE]

dummy <- sparseMatrix(
  i = rep.int(1L, sum(ref_keep)), j = seq_len(sum(ref_keep)), x = rep.int(1, sum(ref_keep)),
  dims = c(1L, sum(ref_keep)), dimnames = list("placeholder", st$cell_id[ref_keep])
)
ref_meta <- as.data.frame(st[ref_keep])
rownames(ref_meta) <- ref_meta$cell_id
ref_obj <- CreateSeuratObject(counts = dummy, meta.data = ref_meta, min.cells = 0, min.features = 0)
ref_obj[["harmony"]] <- CreateDimReducObject(embeddings = ref_harmony, key = "HARMONY_", assay = "RNA")
ref_obj <- FindNeighbors(
  ref_obj, reduction = "harmony", dims = 1:30, k.param = 30, prune.SNN = 1 / 15,
  graph.name = c("healthy_stromal_nn", "healthy_stromal_snn"), verbose = FALSE
)
resolutions <- c(0.2, 0.4, 0.6, 0.8)
for (resolution in resolutions) {
  ref_obj <- FindClusters(
    ref_obj, graph.name = "healthy_stromal_snn", resolution = resolution,
    algorithm = 1, random.seed = 264508,
    cluster.name = paste0("healthy_res_", resolution), verbose = FALSE
  )
}

# The UMAP model is learned only from Healthy cells. LoSc coordinates are an
# out-of-sample transform and therefore cannot alter the Healthy geometry.
umap_fit <- uwot::umap(
  ref_harmony, n_neighbors = 30, min_dist = 0.30, metric = "cosine",
  n_components = 2, init = "spectral", ret_model = TRUE,
  n_threads = 1, verbose = TRUE, seed = 264508
)
ref_umap <- umap_fit$embedding
query_umap <- uwot::umap_transform(query_harmony, umap_fit, n_threads = 1, verbose = TRUE)

ref_out <- as.data.table(ref_obj[[]], keep.rownames = "cell_id")
for (resolution in resolutions) {
  field <- paste0("healthy_res_", resolution)
  ref_out[[field]] <- as.character(ref_out[[field]])
}
ref_out[, `:=`(reference_UMAP_1 = ref_umap[, 1], reference_UMAP_2 = ref_umap[, 2])]

query_out <- copy(st[!ref_keep])
query_out[, `:=`(
  reference_UMAP_1 = query_umap[, 1], reference_UMAP_2 = query_umap[, 2]
)]

metadata_fields <- c(
  "cell_id", "source_row", "library_id", "health", "onset", "celltype",
  "seurat_clusters", "reference_status"
)
ref_keep_fields <- c(metadata_fields, paste0("healthy_res_", resolutions), "reference_UMAP_1", "reference_UMAP_2")
ref_out <- ref_out[, ..ref_keep_fields]
query_out <- query_out[, c(metadata_fields, "reference_UMAP_1", "reference_UMAP_2"), with = FALSE]

marker_modules <- list(
  ADSC_ASPC_identity = c("PDGFRA", "CD34", "PI16", "DPP4", "CFD", "APOD", "C7", "CXCL12"),
  Fibroblast_identity = c("DCN", "LUM", "COL1A1", "COL3A1"),
  Papillary_fibroblast = c("COL18A1", "APCDD1", "WIF1", "PTGDS", "SLPI"),
  Reticular_fibroblast = c("C7", "APOD", "CFD", "CXCL12", "COL14A1"),
  Pericyte = c("RGS5", "PDGFRB", "CSPG4", "RBP1", "MCAM"),
  VSMC = c("ACTA2", "TAGLN", "MYH11", "CNN1"),
  Adipogenic_homeostatic = c("PPARG", "ADIRF", "CEBPA", "LPL", "FABP4", "PLIN1"),
  Inflammatory = c("NFKBIA", "NFKBIZ", "TNFAIP3", "IL6", "CXCL3", "CXCL8", "CSF3"),
  ECM_remodeling = c("POSTN", "SFRP4", "THBS2", "FN1", "COL11A1", "TNC", "MMP11", "TGFBI")
)
marker_annotation <- unique(rbindlist(lapply(names(marker_modules), function(module) {
  data.table(module = module, marker = marker_modules[[module]])
})))

rna <- h5[["assays"]][["RNA"]]
features <- rna[["features"]]$read()
marker_annotation[, feature_index_0 := match(marker, features) - 1L]
marker_annotation[, present := !is.na(feature_index_0)]
present_markers <- unique(marker_annotation[present == TRUE, marker])
present_index <- match(present_markers, features) - 1L

expr_dense <- matrix(
  0, nrow = nrow(st), ncol = length(present_markers),
  dimnames = list(st$cell_id, present_markers)
)
expr <- rna[["data"]]
indptr <- as.numeric(expr[["indptr"]]$read())
indices_ds <- expr[["indices"]]
data_ds <- expr[["data"]]
local_row_map <- integer(nrow(md))
local_row_map[stromal_keep] <- seq_len(nrow(st))
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
    keep <- stromal_keep[source_rows] & !is.na(marker_columns)
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
all_out <- rbindlist(list(ref_out, query_out), fill = TRUE)
all_out <- merge(all_out, expr_dt, by = "cell_id", all.x = TRUE, sort = FALSE)

# Scores are z-standardized per gene using Healthy reference cells only, then
# applied unchanged to both reference and query cells.
for (gene in present_markers) {
  mu <- mean(all_out[health == "Healthy", get(gene)])
  sigma <- sd(all_out[health == "Healthy", get(gene)])
  score_name <- paste0("z_", gene)
  if (is.na(sigma) || sigma == 0) {
    all_out[[score_name]] <- 0
  } else {
    all_out[[score_name]] <- pmax(-4, pmin(4, (all_out[[gene]] - mu) / sigma))
  }
}
for (module in names(marker_modules)) {
  genes <- intersect(marker_modules[[module]], present_markers)
  z_fields <- paste0("z_", genes)
  all_out[[paste0(module, "_score")]] <- rowMeans(as.matrix(all_out[, ..z_fields]))
}

ref_profile <- all_out[health == "Healthy"]
resolution_summary <- rbindlist(lapply(resolutions, function(resolution) {
  field <- paste0("healthy_res_", resolution)
  ref_profile[, .(
    n_cells = .N,
    n_donors = uniqueN(library_id),
    fibroblast_fraction = mean(celltype == "Fibroblasts"),
    pericyte_fraction = mean(celltype == "Pericytes"),
    vsmc_fraction = mean(celltype == "Smooth Muscle Cells"),
    ADSC_ASPC_identity_score = mean(ADSC_ASPC_identity_score),
    Fibroblast_identity_score = mean(Fibroblast_identity_score),
    Papillary_fibroblast_score = mean(Papillary_fibroblast_score),
    Reticular_fibroblast_score = mean(Reticular_fibroblast_score),
    Pericyte_score = mean(Pericyte_score),
    VSMC_score = mean(VSMC_score),
    Adipogenic_homeostatic_score = mean(Adipogenic_homeostatic_score),
    Inflammatory_score = mean(Inflammatory_score),
    ECM_remodeling_score = mean(ECM_remodeling_score)
  ), by = .(cluster = get(field))][, resolution := resolution]
}))
setcolorder(resolution_summary, c("resolution", "cluster", setdiff(names(resolution_summary), c("resolution", "cluster"))))

marker_summary <- rbindlist(lapply(resolutions, function(resolution) {
  field <- paste0("healthy_res_", resolution)
  rbindlist(lapply(present_markers, function(marker) {
    ref_profile[, .(
      n_cells = .N,
      n_expressing = sum(get(marker) > 0),
      pct_exp = 100 * mean(get(marker) > 0),
      avg_exp = mean(get(marker))
    ), by = .(cluster = get(field))][, `:=`(resolution = resolution, marker = marker)]
  }))
}))
marker_summary[, avg_exp_scaled := {
  z <- as.numeric(scale(avg_exp))
  if (all(is.na(z))) rep(0, .N) else pmax(-2.5, pmin(2.5, z))
}, by = .(resolution, marker)]

fwrite(all_out, file.path(output_dir, "stromal_reference_query_expression.csv.gz"))
fwrite(resolution_summary, file.path(output_dir, "healthy_stromal_resolution_summary.csv"))
fwrite(marker_summary, file.path(output_dir, "healthy_stromal_cluster_marker_summary.csv"))
fwrite(marker_annotation, file.path(output_dir, "stromal_marker_presence.csv"))
fwrite(
  rbindlist(lapply(names(marker_modules), function(module) {
    data.table(module = module, marker = marker_modules[[module]], present = marker_modules[[module]] %chin% present_markers)
  })),
  file.path(output_dir, "stromal_module_definitions.csv")
)

writeLines(c(
  sprintf("stromal_cells=%s", nrow(st)),
  sprintf("healthy_reference_cells=%s", sum(ref_keep)),
  sprintf("losc_query_cells=%s", sum(!ref_keep)),
  sprintf("healthy_reference_donors=%s", uniqueN(st[ref_keep, library_id])),
  sprintf("losc_query_donors=%s", uniqueN(st[!ref_keep, library_id])),
  "reference_definition=Healthy cells only",
  "latent_space=author Harmony dimensions 1:30",
  "healthy_clustering=k30 SNN; Louvain; resolutions 0.2,0.4,0.6,0.8",
  "reference_umap=Healthy-only uwot model; cosine; k30; min.dist0.30; seed264508",
  "query_umap=out-of-sample uwot transform; LoSc cells cannot change reference geometry",
  "expression=processed RNA normalized data layer",
  "score_calibration=per-gene mean/SD fitted on Healthy reference and applied unchanged to LoSc query"
), file.path(output_dir, "healthy_stromal_reference_provenance.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "healthy_stromal_reference_sessionInfo.txt"))
message("Healthy stromal reference prepared: ", sum(ref_keep), " reference + ", sum(!ref_keep), " query cells")
