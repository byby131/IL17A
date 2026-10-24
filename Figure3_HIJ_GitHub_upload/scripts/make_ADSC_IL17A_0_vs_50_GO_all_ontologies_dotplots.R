#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(GO.db)
  library(ggplot2)
  library(patchwork)
  library(svglite)
  library(ragg)
  library(grid)
})

options(stringsAsFactors = FALSE)

if (!"Arial" %in% names(grDevices::pdfFonts())) {
  grDevices::pdfFonts(Arial = grDevices::pdfFonts("sans")[[1L]])
}
if (!"Arial" %in% names(grDevices::postscriptFonts())) {
  grDevices::postscriptFonts(Arial = grDevices::postscriptFonts("sans")[[1L]])
}
if (capabilities("aqua") && !"Arial" %in% names(grDevices::quartzFonts())) {
  grDevices::quartzFonts(
    Arial = grDevices::quartzFont(c("Arial", "Arial Bold", "Arial Italic", "Arial Bold Italic"))
  )
}

project_dir <- normalizePath(".", mustWork = TRUE)
deg_path <- file.path(
  project_dir,
  "outputs/ADSC_IL17A_0_vs_50_DEG_heatmaps/source_data",
  "ADSC_IL17A_0_vs_50_full_DEG_statistics.csv"
)
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_0_vs_50_GO_all_ontologies")
source_dir <- file.path(output_dir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(deg_path)) stop("The 0 versus 50 ng/ml DEG statistics table is missing.")

deg <- data.table::fread(deg_path, check.names = FALSE)
required_columns <- c(
  "ensembl_gene_id", "external_gene_name", "gene_label", "logFC", "PValue", "FDR"
)
if (!all(required_columns %in% names(deg))) {
  stop("The DEG statistics table does not contain all required columns.")
}
if (any(!is.finite(deg$logFC)) || any(!is.finite(deg$FDR))) {
  stop("Non-finite logFC or FDR values are present in the DEG table.")
}

up_genes_all <- unique(deg[FDR < 0.05 & logFC > 0, ensembl_gene_id])
down_genes_all <- unique(deg[FDR < 0.05 & logFC < 0, ensembl_gene_id])
symbol_lookup <- setNames(deg$gene_label, deg$ensembl_gene_id)

# Query propagated GO annotations once, then construct a separate tested-gene
# universe and BH correction family for BP, MF, and CC.
go_all <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(deg$ensembl_gene_id),
  keytype = "ENSEMBL",
  columns = c("GOALL", "ONTOLOGYALL")
)
go_all <- unique(data.table::as.data.table(go_all)[
  !is.na(GOALL) & nzchar(GOALL) & ONTOLOGYALL %chin% c("BP", "MF", "CC"),
  .(ensembl_gene_id = ENSEMBL, GO_ID = GOALL, ontology_code = ONTOLOGYALL)
])
if (nrow(go_all) == 0L) stop("No propagated GO annotation was recovered.")

ontology_info <- data.table::data.table(
  ontology_code = c("BP", "MF", "CC"),
  ontology_name = c("Biological Process", "Molecular Function", "Cellular Component")
)

requested_ontology <- toupper(trimws(Sys.getenv("GO_ONTOLOGY", unset = "")))
if (nzchar(requested_ontology)) {
  if (!requested_ontology %chin% ontology_info$ontology_code) {
    stop("GO_ONTOLOGY must be one of BP, MF, or CC.")
  }
  ontology_info <- ontology_info[ontology_code == requested_ontology]
  go_all <- go_all[ontology_code == requested_ontology]
  invisible(gc())
}

shorten_label <- function(x) {
  x <- sub("^positive regulation of ", "pos. reg. of ", x)
  x <- sub("^negative regulation of ", "neg. reg. of ", x)
  x <- sub("^regulation of ", "reg. of ", x)
  x <- gsub("extracellular matrix", "ECM", x, fixed = TRUE)
  x <- gsub("transmembrane transporter activity", "transmembrane transport activity", x, fixed = TRUE)
  x <- gsub("protein-containing complex", "protein complex", x, fixed = TRUE)
  x
}

append_display_genes <- function(x, go_map, input_genes) {
  x[, geneID := vapply(
    GO_ID,
    function(go_id) {
      ids <- go_map[GO_ID == go_id & ensembl_gene_id %chin% input_genes, ensembl_gene_id]
      paste(unique(ids), collapse = "/")
    },
    character(1L)
  )]
  x[, gene_symbols := vapply(
    strsplit(geneID, "/", fixed = TRUE),
    function(ids) {
      symbols <- symbol_lookup[ids]
      symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
      paste(unique(symbols), collapse = "/")
    },
    character(1L)
  )]
  x
}

prepare_plot_data <- function(x) {
  x[, process_label := shorten_label(Description)]
  x[, process_label := as.character(process_label)]
  x[, plot_x := 1]
  x[, plot_y := rev(seq_len(.N))]
  x
}

save_vector_and_raster <- function(plot, prefix, width_mm, height_mm) {
  draw_plot <- function() print(plot)

  # Cairo PDF remains a portable vector branch; Quartz preserves installed
  # Arial as editable text on macOS.
  if (isTRUE(getOption("go.all.use_cairo_pdf", FALSE))) {
    grDevices::cairo_pdf(
      paste0(prefix, ".pdf"), width = width_mm / 25.4,
      height = height_mm / 25.4, family = "Arial"
    )
  } else if (capabilities("aqua")) {
    grDevices::quartz(
      type = "pdf", file = paste0(prefix, ".pdf"),
      width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial"
    )
  } else {
    grDevices::pdf(
      paste0(prefix, ".pdf"), width = width_mm / 25.4,
      height = height_mm / 25.4, family = "Arial", useDingbats = FALSE
    )
  }
  draw_plot()
  grDevices::dev.off()

  svglite::svglite(
    paste0(prefix, ".svg"), width = width_mm / 25.4,
    height = height_mm / 25.4, system_fonts = list(sans = "Arial")
  )
  draw_plot()
  grDevices::dev.off()

  ragg::agg_tiff(
    paste0(prefix, ".tiff"), width = width_mm / 25.4,
    height = height_mm / 25.4, units = "in", res = 600,
    compression = "lzw", background = "white"
  )
  draw_plot()
  grDevices::dev.off()

  ragg::agg_png(
    paste0(prefix, "_preview.png"), width = width_mm / 25.4,
    height = height_mm / 25.4, units = "in", res = 300,
    background = "white"
  )
  draw_plot()
  grDevices::dev.off()
}

analysis_summary <- list()
qa_all <- character()

for (i in seq_len(nrow(ontology_info))) {
  ont_code <- ontology_info$ontology_code[i]
  ont_name <- ontology_info$ontology_name[i]
  message("Processing GO ", ont_code, " (", ont_name, ") ...")

  go_map <- unique(go_all[ontology_code == ont_code, .(ensembl_gene_id, GO_ID)])
  if (nrow(go_map) == 0L) stop("No GO ", ont_code, " annotation was recovered.")

  universe_genes <- unique(go_map$ensembl_gene_id)
  universe_n <- length(universe_genes)
  term_sizes <- go_map[, .(term_size = .N), by = GO_ID]
  term_sizes <- term_sizes[term_size >= 10L & term_size <= 500L]

  run_ora <- function(input_genes, direction) {
    input_genes <- intersect(unique(input_genes), universe_genes)
    input_n <- length(input_genes)
    if (input_n < 10L) stop("Too few GO-annotated genes for ", ont_code, " ", direction, ".")

    overlaps <- go_map[
      ensembl_gene_id %chin% input_genes,
      .(overlap = .N),
      by = GO_ID
    ]
    result <- merge(term_sizes, overlaps, by = "GO_ID", all = FALSE)
    result <- result[overlap > 0L]
    result[, pvalue := stats::phyper(
      overlap - 1L,
      term_size,
      universe_n - term_size,
      input_n,
      lower.tail = FALSE
    )]
    result[, p.adjust := stats::p.adjust(pvalue, method = "BH")]
    result[, Description := AnnotationDbi::Term(GO.db::GOTERM[GO_ID])]
    result <- result[!is.na(Description) & nzchar(Description)]
    result[, GeneRatio := overlap / input_n]
    result[, GeneRatio_label := sprintf("%d/%d", overlap, input_n)]
    result[, BgRatio_label := sprintf("%d/%d", term_size, universe_n)]
    result[, FoldEnrichment := GeneRatio / (term_size / universe_n)]
    result[, Direction := direction]
    result[, Ontology := ont_name]
    result[, Ontology_code := ont_code]
    data.table::setorder(result, p.adjust, -overlap, Description)
    result
  }

  up_full <- run_ora(up_genes_all, "Up in 50 ng/ml")
  down_full <- run_ora(down_genes_all, "Down in 50 ng/ml")
  up_sig <- up_full[p.adjust < 0.05]
  down_sig <- down_full[p.adjust < 0.05]
  if (nrow(up_sig) == 0L || nrow(down_sig) == 0L) {
    stop("No BH-adjusted significant GO-", ont_code, " terms were detected for one direction.")
  }

  top_n <- 10L
  up_plot <- data.table::copy(up_sig[seq_len(min(top_n, .N))])
  down_plot <- data.table::copy(down_sig[seq_len(min(top_n, .N))])
  up_plot <- append_display_genes(up_plot, go_map, intersect(up_genes_all, universe_genes))
  down_plot <- append_display_genes(down_plot, go_map, intersect(down_genes_all, universe_genes))
  up_plot <- prepare_plot_data(up_plot)
  down_plot <- prepare_plot_data(down_plot)

  display_ratio_range <- range(c(up_plot$GeneRatio, down_plot$GeneRatio), finite = TRUE)
  display_p_range <- range(c(up_plot$p.adjust, down_plot$p.adjust), finite = TRUE)
  if (any(!is.finite(display_p_range)) || any(display_p_range <= 0)) {
    stop("Displayed adjusted P values must be finite and strictly positive for log10 scaling.")
  }
  if (diff(display_ratio_range) == 0) {
    display_ratio_range <- display_ratio_range + c(-1, 1) * max(display_ratio_range[1] * 0.05, 1e-4)
  }
  if (diff(log10(display_p_range)) == 0) {
    display_p_range <- display_p_range * c(0.8, 1.2)
  }
  display_p_breaks <- 10^seq(
    log10(display_p_range[1L]),
    log10(display_p_range[2L]),
    length.out = 3L
  )

  dot_colors <- c("#173CCB", "#6A28C2", "#C3159A", "#FF1744")

  make_go_panel <- function(plot_dt, title_text) {
    ggplot(plot_dt, aes(y = plot_y)) +
      annotate(
        "rect",
        xmin = 0.955,
        xmax = 1.045,
        ymin = 0.5,
        ymax = nrow(plot_dt) + 0.5,
        fill = NA,
        color = "#555555",
        linewidth = 0.45
      ) +
      geom_text(
        aes(x = 0.925, label = process_label),
        hjust = 1,
        family = "Arial",
        size = 10 / ggplot2::.pt,
        color = "black",
        lineheight = 0.92
      ) +
      geom_point(
        aes(x = 1, size = GeneRatio, color = p.adjust),
        alpha = 1,
        stroke = 0
      ) +
      scale_size_continuous(
        name = "GeneRatio",
        range = c(2.5, 5.8),
        limits = display_ratio_range,
        breaks = scales::pretty_breaks(n = 3),
        guide = guide_legend(
          order = 1,
          title.position = "top",
          title.hjust = 0,
          nrow = 1,
          override.aes = list(color = "black", alpha = 1)
        )
      ) +
      scale_color_gradientn(
        name = "P.adjust",
        colors = dot_colors,
        limits = display_p_range,
        trans = "log10",
        breaks = display_p_breaks,
        labels = scales::label_scientific(digits = 1),
        guide = guide_colorbar(
          order = 2,
          title.position = "top",
          title.hjust = 0,
          direction = "horizontal",
          barwidth = grid::unit(34, "mm"),
          barheight = grid::unit(3.2, "mm"),
          frame.colour = "#555555",
          ticks.colour = "#555555"
        )
      ) +
      scale_x_continuous(limits = c(-0.83, 1.12), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0.5, nrow(plot_dt) + 0.5), expand = c(0, 0)) +
      coord_cartesian(clip = "off") +
      labs(title = title_text, x = NULL, y = NULL) +
      theme_classic(base_family = "Arial", base_size = 10) +
      theme(
        plot.title = element_text(
          family = "Arial", size = 10, face = "bold", hjust = 0.5,
          color = "black", margin = margin(b = 2.0, unit = "mm")
        ),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom",
        legend.box = "horizontal",
        legend.box.just = "center",
        legend.title = element_text(family = "Arial", size = 10, face = "plain", color = "black"),
        legend.text = element_text(family = "Arial", size = 10, color = "black"),
        legend.key.height = grid::unit(4.0, "mm"),
        legend.spacing.x = grid::unit(2.5, "mm"),
        legend.margin = margin(1.0, 0, 0, 0, unit = "mm"),
        plot.margin = margin(2.0, 2.0, 1.0, 2.0, unit = "mm")
      )
  }

  p_up <- make_go_panel(up_plot, "Upregulated in 50 ng/ml")
  p_down <- make_go_panel(down_plot, "Downregulated in 50 ng/ml")
  combined <- p_up + p_down +
    patchwork::plot_layout(ncol = 2, widths = c(1, 1), guides = "collect") +
    patchwork::plot_annotation(
      title = paste0("GO ", ont_name),
      theme = theme(
        plot.title = element_text(
          family = "Arial", size = 10, face = "bold", hjust = 0.5,
          color = "black", margin = margin(b = 0.8, unit = "mm")
        )
      )
    ) &
    theme(legend.position = "bottom")

  fig_width_mm = 183
  width_mm <- fig_width_mm
  height_mm <- 112
  prefix <- file.path(
    output_dir,
    paste0("ADSC_IL17A_0_vs_50_GO_", ont_code, "_up_down_dotplot")
  )
  save_vector_and_raster(combined, prefix, width_mm, height_mm)

  full_terms <- data.table::rbindlist(list(up_full, down_full), use.names = TRUE, fill = TRUE)
  displayed <- data.table::rbindlist(list(up_plot, down_plot), use.names = TRUE, fill = TRUE)
  displayed[, process_label := as.character(process_label)]
  data.table::fwrite(
    full_terms,
    file.path(source_dir, paste0("ADSC_IL17A_0_vs_50_GO_", ont_code, "_all_tested_terms.csv"))
  )
  data.table::fwrite(
    displayed,
    file.path(source_dir, paste0("ADSC_IL17A_0_vs_50_GO_", ont_code, "_displayed_top10_each_direction.csv"))
  )

  analysis_summary[[ont_code]] <- data.table::data.table(
    ontology_code = ont_code,
    ontology_name = ont_name,
    tested_genes = nrow(deg),
    annotated_universe = universe_n,
    upregulated_FDR_lt_0.05 = length(up_genes_all),
    downregulated_FDR_lt_0.05 = length(down_genes_all),
    upregulated_annotated = length(intersect(up_genes_all, universe_genes)),
    downregulated_annotated = length(intersect(down_genes_all, universe_genes)),
    significant_up_terms = nrow(up_sig),
    significant_down_terms = nrow(down_sig),
    displayed_up_terms = nrow(up_plot),
    displayed_down_terms = nrow(down_plot)
  )

  qa_all <- c(
    qa_all,
    paste0("GO ", ont_code, " - ", ont_name),
    sprintf("  Annotated tested-gene universe: %s", format(universe_n, big.mark = ",")),
    sprintf(
      "  Significant upregulated genes: %s; annotated: %s",
      format(length(up_genes_all), big.mark = ","),
      format(length(intersect(up_genes_all, universe_genes)), big.mark = ",")
    ),
    sprintf(
      "  Significant downregulated genes: %s; annotated: %s",
      format(length(down_genes_all), big.mark = ","),
      format(length(intersect(down_genes_all, universe_genes)), big.mark = ",")
    ),
    sprintf("  Significant enriched up terms: %d", nrow(up_sig)),
    sprintf("  Significant enriched down terms: %d", nrow(down_sig)),
    sprintf("  Displayed: %d up and %d down terms", nrow(up_plot), nrow(down_plot)),
    ""
  )

  rm(
    go_map, universe_genes, term_sizes, up_full, down_full, up_sig, down_sig,
    up_plot, down_plot, full_terms, displayed, p_up, p_down, combined
  )
  invisible(gc())
}

summary_dt <- data.table::rbindlist(analysis_summary, use.names = TRUE, fill = TRUE)
summary_suffix <- if (nzchar(requested_ontology)) paste0("_", requested_ontology) else "_BP_MF_CC"
data.table::fwrite(
  summary_dt,
  file.path(source_dir, paste0("ADSC_IL17A_0_vs_50_GO", summary_suffix, "_analysis_summary.csv"))
)

qa_header <- c(
  "Figure set: GO BP, MF, and CC enrichment for ADSC IL-17A 0 versus 50 ng/ml",
  "Core conclusion: the three figures separately summarize biological-process, molecular-function, and cellular-component over-representation among genes altered by 50 ng/ml IL-17A.",
  "Archetype: quantitative comparison grid; each ontology is an independent figure with upregulated and downregulated panels.",
  "Input differential statistics: edgeR TMM robust quasi-likelihood analysis; contrast is 50 minus 0 ng/ml.",
  sprintf("Tested Ensembl genes: %s", format(nrow(deg), big.mark = ",")),
  "DEG definition: edgeR BH-FDR < 0.05; direction assigned by logFC sign; no absolute fold-change cutoff.",
  "GO annotation: org.Hs.eg.db GOALL mappings, including propagated ancestor terms.",
  "Over-representation test: one-sided hypergeometric test using the ontology-specific GO-annotated tested-gene universe.",
  "Eligible GO term size: 10 to 500 annotated background genes.",
  "Multiple-comparison correction: Benjamini-Hochberg separately for each ontology and expression direction.",
  "Displayed terms: ten smallest adjusted P values among terms with adjusted P < 0.05 in each direction.",
  "Dot size encodes GeneRatio; dot color encodes adjusted P, with lower adjusted P shown in blue and higher adjusted P in red.",
  "No GO term was manually removed, simplified, or reordered by biological preference; only display labels were abbreviated for fit.",
  "The supplied reference image was used for style-only inheritance.",
  "Final figure size: 183 x 112 mm per ontology.",
  "Typography: Arial 10 pt.",
  "Exports per ontology: editable-text PDF/SVG, 600-dpi TIFF, and 300-dpi PNG preview.",
  paste("org.Hs.eg.db version:", as.character(utils::packageVersion("org.Hs.eg.db"))),
  paste("GO.db version:", as.character(utils::packageVersion("GO.db"))),
  "",
  "Ontology-specific counts:"
)
writeLines(
  c(qa_header, qa_all),
  file.path(output_dir, paste0("ADSC_IL17A_0_vs_50_GO", summary_suffix, "_QA_notes.txt"))
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, paste0("ADSC_IL17A_0_vs_50_GO", summary_suffix, "_sessionInfo.txt"))
)

# Independent ontology runs are intentionally supported to limit memory use.
# Once all three summaries exist, consolidate them into one audit table and QA note.
per_ontology_summary_files <- file.path(
  source_dir,
  paste0("ADSC_IL17A_0_vs_50_GO_", c("BP", "MF", "CC"), "_analysis_summary.csv")
)
if (all(file.exists(per_ontology_summary_files))) {
  consolidated_summary <- data.table::rbindlist(
    lapply(per_ontology_summary_files, data.table::fread),
    use.names = TRUE,
    fill = TRUE
  )
  data.table::fwrite(
    consolidated_summary,
    file.path(source_dir, "ADSC_IL17A_0_vs_50_GO_BP_MF_CC_analysis_summary.csv")
  )
  consolidated_lines <- c(
    qa_header,
    "",
    apply(consolidated_summary, 1L, function(z) {
      sprintf(
        "%s (%s): universe=%s; up annotated=%s; down annotated=%s; significant terms up=%s, down=%s; displayed up=%s, down=%s.",
        z[["ontology_code"]], z[["ontology_name"]], z[["annotated_universe"]],
        z[["upregulated_annotated"]], z[["downregulated_annotated"]],
        z[["significant_up_terms"]], z[["significant_down_terms"]],
        z[["displayed_up_terms"]], z[["displayed_down_terms"]]
      )
    })
  )
  writeLines(
    consolidated_lines,
    file.path(output_dir, "ADSC_IL17A_0_vs_50_GO_BP_MF_CC_QA_notes.txt")
  )
}

message("GO BP/MF/CC outputs written to: ", output_dir)
