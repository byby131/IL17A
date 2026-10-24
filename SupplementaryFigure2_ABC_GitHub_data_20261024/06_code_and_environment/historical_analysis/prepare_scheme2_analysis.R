if (.Platform$OS.type == "windows") invisible(Sys.setlocale("LC_CTYPE", "English_United States.utf8"))
invisible(Sys.setlocale("LC_COLLATE", "C"))
.libPaths(c("E:/codex-转录组/adsc_il17_rnaseq/.Rlibrary", .libPaths()))
suppressPackageStartupMessages({
  library(DESeq2)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(GO.db)
  library(fgsea)
})

root <- "E:/codex-转录组"
base_dir <- file.path(root, "analysis/figure4_scheme2")
out <- file.path(base_dir, "source_data")
dir.create(out, recursive=TRUE, showWarnings=FALSE)
write_csv <- function(x, name) write.csv(x, file.path(out, name), row.names=FALSE, na="")

cache <- file.path(root, "HFF-1/round3_20260914/cache")
dds <- readRDS(file.path(cache, "sender_dds.rds"))
ann <- read.csv(file.path(cache, "sender_gene_annotation.csv"), check.names=FALSE)
stopifnot(identical(deparse(design(dds)), "~replicate + time"))
stopifnot(identical(levels(colData(dds)$time), c("0","3","12","24")))
stopifnot(nrow(dds) == 14425L)

# Freeze the current canonical pathway sources from one MSigDB release.
ref_dir <- file.path(root, "adsc_il17_rnaseq/references")
msig_rds <- file.path(ref_dir, "msigdb.2025.1.Hs.rds")
kegg_gmt_file <- file.path(ref_dir, "c2.cp.kegg_medicus.v2025.1.Hs.entrez.gmt")
reactome_gmt_file <- file.path(ref_dir, "c2.cp.reactome.v2025.1.Hs.entrez.gmt")

if (!file.exists(kegg_gmt_file)) {
  msig <- readRDS(msig_rds)
  keep <- which(msig[["gs_collection"]] == "C2" & msig[["gs_subcollection"]] == "CP:KEGG_MEDICUS")
  kd <- unique(msig[keep, c("gs_name", "db_ncbi_gene")])
  kd <- kd[!is.na(kd$db_ncbi_gene) & nzchar(as.character(kd$db_ncbi_gene)), ]
  kp <- split(as.character(kd$db_ncbi_gene), kd$gs_name)
  kp <- lapply(kp, unique)
  lines <- vapply(names(kp), function(nm) paste(c(nm, "MSigDB_2025.1_Hs_CP_KEGG_MEDICUS", kp[[nm]]), collapse="\t"), character(1))
  writeLines(lines, kegg_gmt_file, useBytes=TRUE)
  rm(msig, kd, kp); invisible(gc())
}

kegg <- fgsea::gmtPathways(kegg_gmt_file)
reactome <- fgsea::gmtPathways(reactome_gmt_file)
canonical <- list(KEGG_MEDICUS=kegg, Reactome=reactome)

# Identical Ensembl-to-Entrez ranking map used for the prior Hallmark analysis.
ens <- sub("[.][0-9]+$", "", rownames(dds))
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), keytype="ENSEMBL", columns=c("ENTREZID","SYMBOL")))
map <- map[!is.na(map$ENTREZID),]
ambiguous_ens <- names(which(table(map$ENSEMBL) > 1))
map_unique <- map[!map$ENSEMBL %in% ambiguous_ens & !duplicated(map$ENSEMBL),]
entrez_to_symbol <- setNames(map_unique$SYMBOL, map_unique$ENTREZID)

contrast_def <- list(
  `3_vs_0`=c("time","3","0"),
  `12_vs_0`=c("time","12","0"),
  `24_vs_0`=c("time","24","0"),
  `12_vs_3`=c("time","12","3"),
  `24_vs_3`=c("time","24","3")
)

gsea_rows <- list()
rank_audit <- list()
de_results <- list()
counter <- 0L
for (cn in names(contrast_def)) {
  rr <- as.data.frame(results(dds, contrast=contrast_def[[cn]], independentFiltering=FALSE, alpha=0.05))
  rr$gene_id <- rownames(rr)
  de_results[[cn]] <- rr
  r <- data.frame(ENSEMBL=sub("[.][0-9]+$", "", rownames(dds)), stat=rr$stat, stringsAsFactors=FALSE)
  r <- merge(r, map_unique[,c("ENSEMBL","ENTREZID","SYMBOL")], by="ENSEMBL")
  r <- r[is.finite(r$stat),]
  r <- r[order(abs(r$stat), decreasing=TRUE),]
  r <- r[!duplicated(r$ENTREZID),]
  ranks <- r$stat; names(ranks) <- r$ENTREZID; ranks <- sort(ranks, decreasing=TRUE)
  rank_audit[[cn]] <- data.frame(contrast=cn, tested_genes=nrow(dds), unique_Entrez_ranked=length(ranks), stringsAsFactors=FALSE)
  for (lib in names(canonical)) {
    counter <- counter + 1L
    set.seed(20260921L + counter)
    fg <- fgseaMultilevel(pathways=canonical[[lib]], stats=ranks, minSize=15, maxSize=500, eps=0, nproc=1)
    fg <- as.data.frame(fg)
    fg$library <- lib
    fg$contrast <- cn
    fg$leadingEdge_Entrez <- vapply(fg$leadingEdge, paste, collapse=";", FUN.VALUE=character(1))
    fg$leadingEdge_symbols <- vapply(fg$leadingEdge, function(ids) {
      ss <- unname(entrez_to_symbol[as.character(ids)])
      paste(unique(ss[!is.na(ss) & nzchar(ss)]), collapse=";")
    }, FUN.VALUE=character(1))
    fg$leadingEdge <- NULL
    gsea_rows[[paste(cn,lib,sep="_")]] <- fg
  }
}

gsea <- do.call(rbind, gsea_rows)
gsea$padj_library_contrast <- gsea$padj
gsea$padj_combined_canonical_contrast <- ave(gsea$pval, gsea$contrast, FUN=function(x) p.adjust(x, method="BH"))
gsea <- gsea[order(gsea$contrast, gsea$padj_combined_canonical_contrast, gsea$padj_library_contrast, -abs(gsea$NES)),]
write_csv(gsea, "canonical_GSEA_KEGG_MEDICUS_Reactome_all_5_contrasts.csv")
write_csv(do.call(rbind, rank_audit), "canonical_GSEA_rank_audit.csv")

# Apply the same strict temporal classification to every canonical pathway.
metric_cols <- c("NES","pval","padj_library_contrast","padj_combined_canonical_contrast","size")
keys <- unique(gsea[,c("library","pathway")])
wide_parts <- lapply(names(contrast_def), function(cn) {
  z <- gsea[gsea$contrast == cn, c("library","pathway",metric_cols)]
  names(z)[-(1:2)] <- paste0(metric_cols, "_", cn)
  z
})
wide <- Reduce(function(x,y) merge(x,y,by=c("library","pathway"),all=TRUE), wide_parts)
sig_pos <- function(nes, q) is.finite(nes) & nes > 0 & is.finite(q) & q < 0.05
sig_neg <- function(nes, q) is.finite(nes) & nes < 0 & is.finite(q) & q < 0.05
qcol <- function(cn) wide[[paste0("padj_combined_canonical_contrast_",cn)]]
ncolv <- function(cn) wide[[paste0("NES_",cn)]]
wide$early_response <- sig_pos(ncolv("3_vs_0"), qcol("3_vs_0"))
wide$late_baseline_gain <- sig_pos(ncolv("12_vs_0"),qcol("12_vs_0")) | sig_pos(ncolv("24_vs_0"),qcol("24_vs_0"))
wide$late_vs_3_gain <- sig_pos(ncolv("12_vs_3"),qcol("12_vs_3")) | sig_pos(ncolv("24_vs_3"),qcol("24_vs_3"))
wide$late_other_time_not_negative <- !sig_neg(ncolv("12_vs_0"),qcol("12_vs_0")) & !sig_neg(ncolv("24_vs_0"),qcol("24_vs_0"))
wide$late_emergent <- !wide$early_response & wide$late_baseline_gain & wide$late_vs_3_gain & wide$late_other_time_not_negative
wide$temporal_class <- ifelse(wide$early_response,"Early response",ifelse(wide$late_emergent,"Late-emergent response","Other"))
wide$late_screen_FDR <- pmin(qcol("12_vs_3"), qcol("24_vs_3"), na.rm=TRUE)
wide$late_screen_NES <- pmax(ifelse(sig_pos(ncolv("12_vs_3"),qcol("12_vs_3")),ncolv("12_vs_3"),-Inf),
                             ifelse(sig_pos(ncolv("24_vs_3"),qcol("24_vs_3")),ncolv("24_vs_3"),-Inf))
wide <- wide[order(factor(wide$temporal_class,levels=c("Late-emergent response","Early response","Other")),wide$late_screen_FDR,-wide$late_screen_NES),]
write_csv(wide, "canonical_GSEA_temporal_classification_all.csv")
late_display <- head(wide[wide$late_emergent,], 8L)
write_csv(late_display, "canonical_GSEA_late_emergent_display.csv")
tgf_audit <- wide[grepl("TGF|TRANSFORMING_GROWTH_FACTOR|SMAD", wide$pathway, ignore.case=TRUE),]
write_csv(tgf_audit, "canonical_GSEA_TGF_SMAD_audit.csv")

# Multidatabase ORA for the objectively selected k=3 modules.
assign_file <- file.path(root, "analysis/figure4_abc/source_data/Figure4C/module_gene_assignments.csv")
assign <- read.csv(assign_file, check.names=FALSE)
gene_ids <- rownames(dds)
map_all <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=unique(ens), keytype="ENSEMBL", columns=c("ENTREZID","SYMBOL")))
map_all <- unique(map_all[!is.na(map_all$ENTREZID),c("ENSEMBL","ENTREZID","SYMBOL")])
ens_n <- table(map_all$ENSEMBL); entrez_n <- table(map_all$ENTREZID)
map_1to1 <- map_all[ens_n[map_all$ENSEMBL] == 1L & entrez_n[map_all$ENTREZID] == 1L,,drop=FALSE]
ens_to_entrez <- setNames(map_1to1$ENTREZID,map_1to1$ENSEMBL)
background_entrez <- unique(unname(ens_to_entrez[ens]))
background_entrez <- background_entrez[!is.na(background_entrez)]
module_entrez <- lapply(c("M1","M2","M3"), function(m) {
  ids <- assign$gene_id[assign$module == m]
  x <- unname(ens_to_entrez[sub("[.][0-9]+$","",ids)])
  unique(x[!is.na(x)])
})
names(module_entrez) <- c("M1","M2","M3")

hallmark <- fgsea::gmtPathways(file.path(ref_dir,"h.all.v2025.1.Hs.entrez.gmt"))
go_map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=background_entrez, keytype="ENTREZID", columns=c("GOALL","ONTOLOGYALL")))
go_map <- unique(go_map[!is.na(go_map$GOALL) & go_map$ONTOLOGYALL == "BP",c("ENTREZID","GOALL")])
go_bp <- lapply(split(go_map$ENTREZID,go_map$GOALL),unique)
go_terms <- AnnotationDbi::Term(GO.db::GOTERM[names(go_bp)])
go_labels <- setNames(unname(go_terms),names(go_terms))

ora_libraries <- list(Hallmark=hallmark,KEGG_MEDICUS=kegg,Reactome=reactome,GO_BP=go_bp)
ora_labels <- list(
  Hallmark=setNames(names(hallmark),names(hallmark)),
  KEGG_MEDICUS=setNames(names(kegg),names(kegg)),
  Reactome=setNames(names(reactome),names(reactome)),
  GO_BP=go_labels
)
for (lib in names(ora_libraries)) {
  pp <- lapply(ora_libraries[[lib]], function(x) intersect(unique(as.character(x)), background_entrez))
  pp <- pp[lengths(pp)>=15L & lengths(pp)<=500L]
  ora_libraries[[lib]] <- pp
}

ora_one <- function(module_name,genes,pathways,library_name,labels) {
  N <- length(background_entrez); genes <- intersect(unique(genes),background_entrez); n <- length(genes)
  rows <- lapply(names(pathways),function(term_id) {
    term_genes <- pathways[[term_id]]; K <- length(term_genes); overlap <- intersect(genes,term_genes); k <- length(overlap)
    p <- phyper(k-1L,K,N-K,n,lower.tail=FALSE)
    odds <- ((k+0.5)*(N-K-n+k+0.5))/((n-k+0.5)*(K-k+0.5))
    label <- labels[[term_id]]
    data.frame(module=module_name,library=library_name,term_id=term_id,term_name=label,
               background_mapped_genes=N,module_mapped_genes=n,term_background_genes=K,
               overlap_genes=k,odds_ratio=odds,enrichment_ratio=(k/n)/(K/N),pvalue=p,
               overlap_entrez=paste(overlap,collapse=";"),stringsAsFactors=FALSE)
  })
  do.call(rbind,rows)
}

ora <- do.call(rbind,lapply(names(module_entrez),function(m) {
  do.call(rbind,lapply(names(ora_libraries),function(lib) {
    ora_one(m,module_entrez[[m]],ora_libraries[[lib]],lib,ora_labels[[lib]])
  }))
}))
ora$padj_within_module_library <- ave(ora$pvalue,interaction(ora$module,ora$library,drop=TRUE),FUN=function(x)p.adjust(x,"BH"))
ora$padj_within_module_all_libraries <- ave(ora$pvalue,ora$module,FUN=function(x)p.adjust(x,"BH"))
ora$padj_global_all_modules <- p.adjust(ora$pvalue,"BH")
ora <- ora[order(ora$module,ora$padj_within_module_library,ora$padj_within_module_all_libraries,-ora$odds_ratio,-ora$overlap_genes),]
write_csv(ora,"module_ORA_all_Hallmark_KEGG_Reactome_GO_BP.csv")

select_top_terms <- function(x,max_terms=3L,jaccard_cutoff=0.5) {
  x <- x[x$padj_within_module_library < 0.05 & x$overlap_genes >= 3L,,drop=FALSE]
  if (nrow(x)==0L) return(x)
  x <- x[order(x$padj_within_module_library,x$padj_within_module_all_libraries,-x$odds_ratio,-x$overlap_genes,x$term_id),]
  chosen <- integer(0); hit_sets <- strsplit(x$overlap_entrez,";",fixed=TRUE)
  for (i in seq_len(nrow(x))) {
    if (!length(chosen)) chosen <- i else {
      jac <- vapply(chosen,function(j) length(intersect(hit_sets[[i]],hit_sets[[j]]))/length(union(hit_sets[[i]],hit_sets[[j]])),numeric(1))
      if (all(jac < jaccard_cutoff)) chosen <- c(chosen,i)
    }
    if (length(chosen)>=max_terms) break
  }
  x[chosen,,drop=FALSE]
}
top_terms <- do.call(rbind,lapply(split(ora,ora$module),select_top_terms))
if (nrow(top_terms)) top_terms$display_rank <- ave(seq_len(nrow(top_terms)),top_terms$module,FUN=seq_along)
write_csv(top_terms,"module_ORA_top_nonredundant_terms.csv")
module_tgf <- ora[grepl("TGF|TRANSFORMING GROWTH FACTOR|SMAD",ora$term_name,ignore.case=TRUE) | grepl("TGF|SMAD",ora$term_id,ignore.case=TRUE),]
write_csv(module_tgf,"module_ORA_TGF_SMAD_audit.csv")

# Effector kinetics are retained only because all five map to the objectively
# selected M2 and remain positively regulated at 3, 12 and 24 h vs 0 h.
effectors <- c("CXCL1","CXCL5","CXCL6","CSF3","NFKBIZ")
ann_idx <- match(rownames(dds),ann$gene_id)
gene_name <- ann$gene_name[ann_idx]
names(gene_name) <- rownames(dds)
eff_ids <- names(gene_name)[gene_name %in% effectors]
stopifnot(length(unique(gene_name[eff_ids])) == length(effectors))
norm <- counts(dds,normalized=TRUE)
v <- assay(vst(dds,blind=FALSE))
meta <- as.data.frame(colData(dds))
meta$sample_id <- rownames(meta)
kin <- do.call(rbind,lapply(eff_ids,function(gid) {
  module <- assign$module[match(gid,assign$gene_id)]
  stopifnot(length(module)==1L,!is.na(module),module=="M2")
  data.frame(gene_id=gid,gene_name=gene_name[[gid]],module=module,sample_id=colnames(dds),
             time_h=as.integer(as.character(meta[colnames(dds),"time"])),replicate=as.character(meta[colnames(dds),"replicate"]),
             normalized_count=as.numeric(norm[gid,]),log2_normalized_count=log2(as.numeric(norm[gid,])+1),
             VST_expression=as.numeric(v[gid,]),stringsAsFactors=FALSE)
}))
kin$gene_name <- factor(kin$gene_name,levels=effectors)
kin <- kin[order(kin$gene_name,kin$time_h,kin$replicate),]
kin$gene_name <- as.character(kin$gene_name)
write_csv(kin,"effector_gene_kinetics_samples.csv")

eff_stats <- do.call(rbind,lapply(eff_ids,function(gid) {
  ix <- match(gid,rownames(dds)); nm <- gene_name[[gid]]
  lrt <- read.csv(file.path(root,"analysis/adsc_timecourse_counts/DESeq2_LRT_any_time_effect.csv"),check.names=FALSE)
  li <- match(gid,lrt$gene_id)
  data.frame(gene_id=gid,gene_name=nm,module=assign$module[match(gid,assign$gene_id)],LRT_padj=lrt$padj[li],
             log2FC_3_vs_0=de_results[["3_vs_0"]]$log2FoldChange[ix],padj_3_vs_0=de_results[["3_vs_0"]]$padj[ix],
             log2FC_12_vs_0=de_results[["12_vs_0"]]$log2FoldChange[ix],padj_12_vs_0=de_results[["12_vs_0"]]$padj[ix],
             log2FC_24_vs_0=de_results[["24_vs_0"]]$log2FoldChange[ix],padj_24_vs_0=de_results[["24_vs_0"]]$padj[ix],stringsAsFactors=FALSE)
}))
write_csv(eff_stats,"effector_gene_statistics.csv")

summary <- list(
  design="~ replicate + time",
  tested_genes=nrow(dds),
  canonical_sources=c("MSigDB 2025.1 Hs C2:CP:KEGG_MEDICUS","MSigDB 2025.1 Hs C2:CP:REACTOME"),
  canonical_pathways_tested=as.list(table(gsea$library[gsea$contrast=="3_vs_0"])),
  strict_late_emergent_counts=as.list(table(wide$library[wide$late_emergent])),
  strict_late_emergent_total=sum(wide$late_emergent),
  TGF_SMAD_canonical_hits=nrow(tgf_audit),
  module_sizes=as.list(table(assign$module)),
  module_top_terms=if(nrow(top_terms)) split(top_terms[,c("library","term_name","padj_within_module_library")],top_terms$module) else list(),
  effector_genes=eff_stats,
  limitations=c("No time-matched untreated control at 3, 12 or 24 h.","Counts-only Salmon gene aggregation lacks sample-specific effective-length offsets.")
)
jsonlite::write_json(summary,file.path(out,"scheme2_analysis_summary.json"),pretty=TRUE,auto_unbox=TRUE,digits=10)
writeLines(c(
  "# Figure 4 scheme 2 analysis",
  "",
  "B: primary Hallmark GSEA shows only strict early-response programs (3 h vs 0 h NES > 0 and FDR < 0.05).",
  "Secondary canonical GSEA tests the complete MSigDB 2025.1 Hs KEGG_MEDICUS and Reactome libraries for all five contrasts with no pathway preselection.",
  "Strict late-emergent rules are identical to the frozen Hallmark rules and use BH FDR < 0.05 across both canonical libraries within each contrast.",
  "C: the previously selected k=3 modules are unchanged; Hallmark, KEGG_MEDICUS, Reactome and GO-BP ORA use the 14,425 tested genes as the declared background and one-to-one Entrez mappings for computation.",
  "D: CXCL1, CXCL5, CXCL6, CSF3 and NFKBIZ were retained because every gene maps to M2 and remains positively regulated at 3, 12 and 24 h versus 0 h.",
  "No pathway, module, threshold or gene was changed to force a TGF-beta result."
),file.path(out,"methods_scheme2.md"),useBytes=TRUE)

cat(jsonlite::toJSON(summary,pretty=TRUE,auto_unbox=TRUE,digits=5))
