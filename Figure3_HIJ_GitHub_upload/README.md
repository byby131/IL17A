# Figure 3H–J analysis and plotting code

This folder contains the code and frozen result tables for the three computational panels in the current Figure 3. It does **not** contain panels A–G, or the unpublished transcript-level expression matrix.

| Panel | Analysis | Final plotting script | Frozen PDF |
|---|---|---|---|
| H | IL-17A dose-response PCA, 15 samples | `scripts/make_ADSC_IL17A_dose_PCA.R` | `outputs/ADSC_IL17A_Dose_PCA/ADSC_IL17A_Dose_PCA_reference_palette.pdf` |
| I | 0 versus 50 ng/mL GO Molecular Function, up/down | `scripts/reformat_ADSC_IL17A_KEGG_GO_MF_74x104mm.R` | `outputs/ADSC_IL17A_0_vs_50_KEGG_GO_MF_74x104mm/ADSC_IL17A_0_vs_50_GO_MF_top10_up_down_74x104mm.pdf` |
| J | 0 versus 50 ng/mL KEGG, up/down | `scripts/replot_ADSC_IL17A_KEGG_bar_bubble_74x104mm.R` | `outputs/ADSC_IL17A_0_vs_50_KEGG_bubble_uncompressed_74x125mm/ADSC_IL17A_0_vs_50_KEGG_top10_up_down_bubble_reference_palette_74x125mm.pdf` |

The script name for J retains an older `74x104mm` label; the code and final PDF use **74 × 125 mm**. Panel I is 74 × 104 mm. The PDF files are reference outputs, not substitutes for the analysis code.

## Code path

1. `make_ADSC_IL17A_dose_PCA.R` filters and TMM-normalizes transcript-level estimated counts, computes log2-CPM, selects the 5,000 most variable transcripts, and draws H. `plot_Figure3H_from_frozen_coordinates.R` redraws H from the saved coordinates without access to the expression matrix; it does **not** recompute PCA.
2. `make_ADSC_IL17A_0_vs_50_DEG_heatmaps.R` performs the shared gene-level edgeR differential analysis for I/J. It also generates historical Top-DEG heatmaps, which are *not* part of Figure 3H–J. The exported `full_DEG_statistics.csv` is the I/J input.
3. `make_ADSC_IL17A_0_vs_50_GO_all_ontologies_dotplots.R` generates the GO over-representation results. Only its MF result feeds I; BP and CC are historical byproducts, not Figure 3 panels here.
4. `make_ADSC_IL17A_0_vs_50_KEGG_barplot.R` generates the KEGG over-representation results and source tables for J. Its original bar plot is not the final J display.
5. The final I and J scripts above render from the saved full GO-MF and KEGG result tables, respectively. Both retain the actual selected top-ten pathways for each direction.

## Replot from the included frozen results

Run these commands from **this folder** with R and the listed packages installed:

```sh
Rscript scripts/plot_Figure3H_from_frozen_coordinates.R
Rscript scripts/reformat_ADSC_IL17A_KEGG_GO_MF_74x104mm.R
Rscript scripts/replot_ADSC_IL17A_KEGG_bar_bubble_74x104mm.R
```

The I script also writes a historical compact KEGG bar panel; the J script also writes a historical alternate bubble/bar style. Use only the PDF filenames in the table for the current Figure 3.

## Full analytical rerun

To recompute H and the I/J differential analysis, obtain the unpublished estimated-count matrix from the study data custodian and place it at:

`data/private/ADSC_IL17A_Dose_RNAseq_expression_matrix_counts.csv`

The private estimated-count matrix is not included in this package. The frozen transcript-to-gene mapping is supplied losslessly as `work/ADSC_IL17A_Dose_RNAseq/Ensembl_transcript_gene_mapping_current.csv.gz` so that every file fits GitHub's browser-upload limit. Before a full analytical rerun, restore the original CSV from this package's root with:

```sh
gunzip -k work/ADSC_IL17A_Dose_RNAseq/Ensembl_transcript_gene_mapping_current.csv.gz
```

On macOS, the `.gz` file can also be expanded in Finder. The SHA-256 of the restored CSV must be `68d52ec9145d3ec09c9b09135fc95c2a05c5b4d227ed2631cecbf501cae903f2`. Frozen-result replots do not require this step.

From the package root, execute in order:

```sh
Rscript scripts/make_ADSC_IL17A_dose_PCA.R
Rscript scripts/make_ADSC_IL17A_0_vs_50_DEG_heatmaps.R
Rscript scripts/make_ADSC_IL17A_0_vs_50_GO_all_ontologies_dotplots.R
Rscript scripts/make_ADSC_IL17A_0_vs_50_KEGG_barplot.R
Rscript scripts/reformat_ADSC_IL17A_KEGG_GO_MF_74x104mm.R
Rscript scripts/replot_ADSC_IL17A_KEGG_bar_bubble_74x104mm.R
```

Rerunning `enrichKEGG` queries a live KEGG annotation service, so exact terms or adjusted P values may differ from the saved results. The supplied full tested-pathway table is the frozen source for the published J plot.

## Statistical provenance

- H: 15 samples, five doses (0, 10, 50, 100, 200 ng/mL), three samples per dose; CPM > 1 in at least three samples, TMM log2-CPM with prior count 2, top 5,000 transcript variances, centered PCA without variance scaling. Convex hulls are not confidence intervals.
- I/J: six samples (three control and three 50 ng/mL). Transcript estimated counts were summed by mapped Ensembl gene ID. edgeR `filterByExpr` retained 14,069 genes; TMM, robust dispersion and quasi-likelihood GLM (`~ group`), then BH FDR < 0.05. There was no TPM or absolute log2FC threshold for the enrichment inputs.
- I: separate up/down GO-MF over-representation tests. The plot shows ten significant terms per direction; point size is GeneRatio and color is BH-adjusted P.
- J: separate up/down KEGG over-representation tests. The plot shows ten significant pathways per direction; x is rich factor, point size is overlap count, and color is −log10(BH-adjusted P).

Saved `sessionInfo.txt` files document the original software environment. `FILE_MANIFEST.csv` records SHA-256 hashes for the browser-upload version of the staged code and frozen inputs/results, including the compressed annotation. The packaged R scripts use the current working directory as project root; start them from this folder.
