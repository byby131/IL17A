#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(KEGGREST)
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
output_dir <- file.path(project_dir, "outputs/ADSC_IL17A_0_vs_50_KEGG")
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

# Convert the tested Ensembl gene universe to Entrez identifiers. Unique Entrez
# IDs are used so one-to-many annotation mappings cannot duplicate enrichment input.
id_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(deg$ensembl_gene_id),
  keytype = "ENSEMBL",
  columns = c("ENTREZID", "SYMBOL")
)
id_map <- data.table::as.data.table(id_map)
id_map <- id_map[!is.na(ENTREZID) & nzchar(ENTREZID)]
id_map <- unique(id_map, by = c("ENSEMBL", "ENTREZID"))

deg_mapped <- merge(
  deg,
  id_map,
  by.x = "ensembl_gene_id",
  by.y = "ENSEMBL",
  all = FALSE,
  allow.cartesian = TRUE
)

tested_entrez <- unique(deg_mapped$ENTREZID)
up_entrez <- unique(deg_mapped[FDR < 0.05 & logFC > 0, ENTREZID])
down_entrez <- unique(deg_mapped[FDR < 0.05 & logFC < 0, ENTREZID])

if (length(up_entrez) < 10L || length(down_entrez) < 10L) {
  stop("Too few mapped significant genes are available for separate KEGG enrichment.")
}

run_kegg <- function(entrez_ids) {
  clusterProfiler::enrichKEGG(
    gene = entrez_ids,
    universe = tested_entrez,
    organism = "hsa",
    keyType = "ncbi-geneid",
    pvalueCutoff = 1,
    pAdjustMethod = "BH",
    qvalueCutoff = 1,
    minGSSize = 10,
    maxGSSize = 500,
    use_internal_data = FALSE
  )
}

up_result <- run_kegg(up_entrez)
down_result <- run_kegg(down_entrez)

up_full <- data.table::as.data.table(as.data.frame(up_result))
down_full <- data.table::as.data.table(as.data.frame(down_result))
if (nrow(up_full) == 0L || nrow(down_full) == 0L) {
  stop("KEGG enrichment returned no pathways for one direction.")
}

add_direction <- function(x, direction) {
  x <- data.table::copy(x)
  x[, Direction := direction]
  x[, minus_log10_adjusted_p := -log10(pmax(p.adjust, .Machine$double.xmin))]
  x
}

up_full <- add_direction(up_full, "Up in 50 ng/ml")
down_full <- add_direction(down_full, "Down in 50 ng/ml")

up_sig <- up_full[p.adjust < 0.05]
down_sig <- down_full[p.adjust < 0.05]
data.table::setorder(up_sig, p.adjust, -Count, Description)
data.table::setorder(down_sig, p.adjust, -Count, Description)

if (nrow(up_sig) == 0L || nrow(down_sig) == 0L) {
  stop("No BH-adjusted significant KEGG pathways were detected for one direction.")
}

top_n <- 10L
up_plot <- up_sig[seq_len(min(top_n, .N))]
down_plot <- down_sig[seq_len(min(top_n, .N))]

wrap_label <- function(x, width = 43L) {
  vapply(
    x,
    function(value) paste(strwrap(value, width = width), collapse = "\n"),
    character(1L)
  )
}

prepare_plot_data <- function(x) {
  x <- data.table::copy(x)
  x[, pathway_label := wrap_label(Description)]
  x[, pathway_label := factor(pathway_label, levels = rev(pathway_label))]
  score_range <- range(x$minus_log10_adjusted_p, finite = TRUE)
  text_switch <- score_range[1L] + 0.42 * diff(score_range)
  x[, label_color := ifelse(minus_log10_adjusted_p >= text_switch, "white", "black")]
  x
}

up_plot <- prepare_plot_data(up_plot)
down_plot <- prepare_plot_data(down_plot)
common_x_max <- max(c(up_plot$minus_log10_adjusted_p, down_plot$minus_log10_adjusted_p)) * 1.035

make_panel <- function(plot_dt, direction_label, gradient_colors, show_x, show_title) {
  score_range <- range(plot_dt$minus_log10_adjusted_p, finite = TRUE)
  if (diff(score_range) == 0) score_range <- score_range + c(-0.1, 0.1)

  ggplot(plot_dt, aes(x = minus_log10_adjusted_p, y = pathway_label)) +
    geom_col(
      aes(fill = minus_log10_adjusted_p),
      width = 0.76,
      color = NA
    ) +
    geom_text(
      aes(
        x = 0.065,
        label = pathway_label,
        color = label_color
      ),
      hjust = 0,
      vjust = 0.5,
      family = "Arial",
      size = 10 / ggplot2::.pt,
      lineheight = 0.90,
      show.legend = FALSE
    ) +
    scale_fill_gradientn(
      colors = gradient_colors,
      limits = score_range,
      name = "-Log10\n(adjusted P)",
      guide = guide_colorbar(
        title.position = "top",
        title.hjust = 0,
        barwidth = grid::unit(3.2, "mm"),
        barheight = grid::unit(19, "mm"),
        ticks.colour = "#333333",
        frame.colour = "#555555"
      )
    ) +
    scale_color_identity() +
    scale_x_continuous(
      limits = c(0, common_x_max),
      expand = expansion(mult = c(0, 0)),
      breaks = scales::pretty_breaks(n = 5)
    ) +
    labs(
      title = if (show_title) "KEGG pathways" else NULL,
      x = if (show_x) "-Log10(adjusted P value)" else NULL,
      y = direction_label
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(
        family = "Arial", size = 10, face = "bold", hjust = 0.5,
        margin = margin(b = 2.5, unit = "mm")
      ),
      axis.title.x = element_text(
        family = "Arial", size = 10, color = "black",
        margin = margin(t = 1.8, unit = "mm")
      ),
      axis.title.y = element_text(
        family = "Arial", size = 10, color = "black",
        margin = margin(r = 2.0, unit = "mm")
      ),
      axis.text.x = if (show_x) {
        element_text(family = "Arial", size = 10, color = "black")
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) {
        element_line(linewidth = 0.35, color = "black")
      } else {
        element_blank()
      },
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(fill = NA, color = "#353535", linewidth = 0.45),
      panel.grid.major.x = element_line(
        color = "#C9C9C9", linewidth = 0.32, linetype = "dashed"
      ),
      panel.grid.minor = element_blank(),
      legend.position = "right",
      legend.title = element_text(
        family = "Arial", size = 10, face = "plain", color = "black",
        lineheight = 0.95
      ),
      legend.text = element_text(family = "Arial", size = 10, color = "black"),
      legend.margin = margin(0, 0, 0, 1.5, unit = "mm"),
      plot.margin = margin(
        t = if (show_title) 2.5 else 1.0,
        r = 2.5,
        b = if (show_x) 2.5 else 0.8,
        l = 2.5,
        unit = "mm"
      )
    )
}

up_colors <- c("#FDE8E5", "#F8AAA1", "#ED6A5A", "#C92D3A", "#760016")
down_colors <- c("#E6F1FA", "#A9CDE8", "#639FCF", "#2768AE", "#08306B")

p_up <- make_panel(
  up_plot,
  direction_label = "Up in 50 ng/ml",
  gradient_colors = up_colors,
  show_x = FALSE,
  show_title = TRUE
)
p_down <- make_panel(
  down_plot,
  direction_label = "Down in 50 ng/ml",
  gradient_colors = down_colors,
  show_x = TRUE,
  show_title = FALSE
)

combined <- p_up / p_down +
  patchwork::plot_layout(
    heights = c(nrow(up_plot), nrow(down_plot)),
    guides = "keep"
  )

width_mm <- 183
height_mm <- 155
prefix <- file.path(output_dir, "ADSC_IL17A_0_vs_50_KEGG_up_down_barplot")

# Cairo PDF is retained as a portable vector branch. Quartz is preferred on
# macOS to preserve installed Arial as editable vector text.
if (isTRUE(getOption("kegg.bar.use_cairo_pdf", FALSE))) {
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
print(combined)
dev.off()

svglite::svglite(
  paste0(prefix, ".svg"), width = width_mm / 25.4,
  height = height_mm / 25.4, system_fonts = list(sans = "Arial")
)
print(combined)
dev.off()

ragg::agg_tiff(
  paste0(prefix, ".tiff"), width = width_mm / 25.4,
  height = height_mm / 25.4, units = "in", res = 600,
  compression = "lzw", background = "white"
)
print(combined)
dev.off()

ragg::agg_png(
  paste0(prefix, "_preview.png"), width = width_mm / 25.4,
  height = height_mm / 25.4, units = "in", res = 300,
  background = "white"
)
print(combined)
dev.off()

entrez_symbol <- unique(id_map[!is.na(SYMBOL) & nzchar(SYMBOL), .(ENTREZID, SYMBOL)])
symbol_map <- split(entrez_symbol$SYMBOL, entrez_symbol$ENTREZID)
append_symbols <- function(x) {
  x <- data.table::copy(x)
  x[, gene_symbols := vapply(
    strsplit(geneID, "/", fixed = TRUE),
    function(ids) {
      symbols <- unique(unlist(symbol_map[ids], use.names = FALSE))
      symbols <- symbols[!is.na(symbols) & nzchar(symbols)]
      paste(symbols, collapse = "/")
    },
    character(1L)
  )]
  x
}

up_full <- append_symbols(up_full)
down_full <- append_symbols(down_full)
displayed <- data.table::rbindlist(list(
  append_symbols(up_plot),
  append_symbols(down_plot)
), use.names = TRUE, fill = TRUE)
displayed[, pathway_label := as.character(pathway_label)]

data.table::fwrite(
  data.table::rbindlist(list(up_full, down_full), use.names = TRUE, fill = TRUE),
  file.path(source_dir, "ADSC_IL17A_0_vs_50_KEGG_all_tested_pathways.csv")
)
data.table::fwrite(
  displayed,
  file.path(source_dir, "ADSC_IL17A_0_vs_50_KEGG_displayed_top10_each_direction.csv")
)
data.table::fwrite(
  id_map,
  file.path(source_dir, "ADSC_IL17A_0_vs_50_Ensembl_Entrez_symbol_mapping.csv")
)
data.table::fwrite(
  deg_mapped[, .(
    ensembl_gene_id, ENTREZID, SYMBOL, gene_label, logFC, PValue, FDR
  )],
  file.path(source_dir, "ADSC_IL17A_0_vs_50_mapped_DEG_universe.csv")
)

kegg_info <- tryCatch(
  paste(capture.output(KEGGREST::keggInfo("kegg")), collapse = " | "),
  error = function(e) paste("KEGG release query unavailable:", conditionMessage(e))
)

qa_lines <- c(
  "Figure: KEGG pathway enrichment for ADSC IL-17A 0 versus 50 ng/ml",
  "Role: pathway-level interpretation of the transcriptional response to 50 ng/ml IL-17A.",
  "Input differential statistics: edgeR TMM robust quasi-likelihood analysis; contrast is 50 minus 0 ng/ml.",
  sprintf("Tested Ensembl genes: %s", format(nrow(deg), big.mark = ",")),
  sprintf("Unique mapped Entrez genes in enrichment universe: %s", format(length(tested_entrez), big.mark = ",")),
  "DEG definition for enrichment: edgeR BH-FDR < 0.05; direction assigned by the sign of logFC; no absolute fold-change cutoff.",
  sprintf("Mapped upregulated Entrez IDs: %s", format(length(up_entrez), big.mark = ",")),
  sprintf("Mapped downregulated Entrez IDs: %s", format(length(down_entrez), big.mark = ",")),
  "Enrichment: clusterProfiler enrichKEGG, organism hsa, live KEGG annotation, tested-gene universe, hypergeometric test.",
  "Multiple-comparison correction: Benjamini-Hochberg separately within the upregulated and downregulated enrichment analyses.",
  "Displayed pathways: the ten smallest adjusted P values among pathways with adjusted P < 0.05 in each direction.",
  sprintf("Significant upregulated KEGG pathways: %d", nrow(up_sig)),
  sprintf("Significant downregulated KEGG pathways: %d", nrow(down_sig)),
  "Bar length and color depth both encode -log10(adjusted P); red is upregulated and blue is downregulated.",
  "Pathway names are directly labeled within the bars, matching the supplied reference structure; dark bars use white labels and light bars use black labels.",
  "No pathway was manually removed, renamed, or reordered by biological preference.",
  "The supplied reference image was used for style-only inheritance; statistical logic was implemented from the RNA-seq data.",
  sprintf("Final figure size: %d x %d mm", width_mm, height_mm),
  "Typography: Arial 10 pt.",
  "Exports: editable-text PDF/SVG, 600-dpi TIFF, and 300-dpi PNG preview.",
  paste("KEGG database information:", kegg_info)
)
writeLines(
  qa_lines,
  file.path(output_dir, "ADSC_IL17A_0_vs_50_KEGG_QA_notes.txt")
)
capture.output(
  sessionInfo(),
  file = file.path(output_dir, "ADSC_IL17A_0_vs_50_KEGG_sessionInfo.txt")
)

message("KEGG outputs written to: ", output_dir)
