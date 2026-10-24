# Figure 5A–E: code status and reproducibility

The current composite figure has panels A–K. The computational RNA-seq panels traceable to this analysis are **A–E**. Panels F–K contain microscopy, staining quantification, and an experimental schematic; their image-processing/statistics scripts were not found in the archived computational delivery and are not represented as our bioinformatics code here.

| Panel | Data and plotted quantity | Code in this folder |
|---|---|---|
| A | All 17,882 plotted genes, log2FC versus −log10 FDR | `code/plot_A_to_E_from_frozen_source.R` |
| B | All 1,823 DEGs, VST row-Z heatmap | Same |
| C | 24 genes from five functional programs, VST row-Z and gene FDR | Same; original selection rules under `metadata/` |
| D | 14 adipogenic transcription factors, VST row-Z and gene FDR | Same |
| E | Top 10 up- and top 10 down-input KEGG pathways, −log10 FDR | Same; all eligible ORA results under `data/E_KEGG/` |

## Important provenance distinction

The original *executable* scripts that rendered A–E and performed the archived KEGG/GSEA calculations were **not** present among the saved files. Do not label the new plotting script as the historical original. It is a transparent reconstruction that reads the frozen, audited source tables and creates new reference plots without changing any statistic, gene list, pathway selection or matrix value. Its visual layout is not a pixel-identical copy of the manuscript composite. For the precise manuscript artwork, retain the author-approved final figure file separately.

`code/recompute_primary_DESeq2_check.R` independently implements the documented primary model using the saved input counts and compares every modeled gene with the frozen table. It is not the historical original script and writes to `recomputed_analysis/`, leaving the published values unchanged. On the local test environment it retained 17,887 genes and reproduced the archived results with no NA discordance; maximum absolute differences were 1.9×10⁻⁹ for baseMean, 1.1×10⁻⁹ for log2FC, 5.1×10⁻⁸ for lfcSE, 2.1×10⁻⁶ for Wald statistic, 4.1×10⁻⁷ for P and 1.7×10⁻⁶ for BH-FDR. These tiny differences reflect numerical environment differences, not a change to the frozen figure data.

The existing `code/verify_source_data.py`, `metadata/analysis_config_used.R`, and `metadata/C_selection_rules.R` are original saved verification/configuration assets. They are not substitutes for the missing historical GSEA/KEGG scripts. The 2026-10-24 release manifest is preserved in `provenance/Figure5_A-E_GitHub_20261024_manifest.csv`. The top-level `manifest.csv` covers this expanded folder, including new code and its verification outputs.

## Run from this folder

```sh
Rscript code/plot_A_to_E_from_frozen_source.R
Rscript code/recompute_primary_DESeq2_check.R
python code/verify_source_data.py --root . --output validation/current_check.json
```

R requires `data.table`, `ggplot2`, `ggrepel`, `patchwork`, `svglite`, `ragg`, `DESeq2` and their dependencies; `scales` is used for display-only color saturation. The plotting script writes PDF, SVG, PNG and 600-dpi TIFF files under `reproduced_panels/`. The recomputation writes CSV and session information under `recomputed_analysis/`. Running it does **not** rerun GSEA or KEGG, which would require the exact frozen database snapshots and original executable code.

## Public-release check

This folder includes the six-sample count matrix because it was already present in the author-prepared Figure 5 GitHub data package. It does **not** include FASTQ/BAM. Before an actual public push, the authors should confirm that this count matrix and its licence may be released, then add the manuscript/data accession and repository licence. No upload has been performed here.
