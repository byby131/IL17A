# Supplementary Figure 3 E–G: HFF-1 RNA-seq figure source data

This package contains the quantified count input, frozen analytical results, exact figure source values and portable replotting code for **Supplementary Figure 3E–G**. It repackages existing analyses; it does not refit DESeq2 or rerun enrichment. Panels A–D and H–J are outside this package. `metadata/panel_mapping.csv` records the correspondence to the earlier figure assembly.

## Experiment and sample identifiers

HFF-1 cells received conditioned medium from unstimulated ADSCs (**Control CM_FB**, raw samples `FBcontrol_1`–`FBcontrol_3`) or from ADSCs exposed to 50 ng/mL IL-17A for **12 h** (**IL-17A CM_FB**, raw samples `ADSC12hCM_1`–`ADSC12hCM_3`). HFF-1 cells were exposed to conditioned medium for **24 h** before RNA collection. There are three independent biological replicates per group, six samples total; the HFF comparison is **unpaired**. Repeated suffixes across the two groups do not denote pairs. See `metadata/samples.csv`.

Positive log2 fold change, Wald statistic and enrichment direction refer to IL-17A CM_FB relative to Control CM_FB. This contrast concerns RNA from recipient HFF-1 cells, not RNA from ADSCs.

## Files to use

| Panel | Primary source files | Meaning |
|---|---|---|
| E | `source_data/Figure_S3E_PCA.csv`, `Figure_S3E_variance_explained.csv` | Six saved PCA coordinates and variance explained |
| F | `source_data/Figure_S3F_pathway_overview.csv` | All 24 displayed pathways in figure order; exact NES and global FDR |
| G | `source_data/Figure_S3G_sample_scores.csv` | 42 sample/program observations: original scores and heatmap row z-scores |
| G | `source_data/Figure_S3G_raw_scores_matrix.csv`, `Figure_S3G_row_z_matrix.csv` | The same observations in 7 × 6 matrix form |
| G | `source_data/Figure_S3G_GSEA.csv` | Seven original enrichment results, including nonsignificant and negative results |
| G audit | `source_data/Figure_S3G_score_comparison.csv` | Sample-level tests; these are not the GSEA dot significance values |

`data/input/Salmon_gene_counts.csv` is the **unaltered, unnormalized Salmon-derived gene count table**, including fractional estimated counts. It is quantified expression data, **not raw sequencing reads**. FASTQ files, Salmon quant.sf files and effective-length offsets are not included. Raw sequencing repository accessions are not available in these files and have not been invented.

`data/processed/` preserves all tested DESeq2 results, filter status for all input genes, normalized counts, VST values, the PCA gene list, complete mapped Wald ranking, mapping audit, and the full original GSEA result table. Complete families are retained so displayed selections and multiple-testing corrections can be audited. The GSEA table also includes the previously tested SSc signature because it belongs to the original global FDR family; that signature is not an additional panel in Figure 3E–G.

## Frozen analytical methods

- DESeq2 1.46.0, design `~ condition`, parametric dispersion fit. Salmon counts were rounded to integers; retain genes with rounded count ≥10 in at least 3 of 6 samples. No gene effective-length offsets were available. Wald tests, unshrunken log2 fold changes, default Cook handling, `independentFiltering=FALSE`; BH over estimable tests after filtering. Original missing values are retained.
- E: PCA from the 2,000 highest-variance genes in the saved VST matrix, `blind=FALSE`. PCA does not use group labels, but the VST is design-aware. PC1 explains 81.9% and PC2 6.8%. No PCA is recalculated by the replot script.
- GSEA: finite, converged, P-value-estimable Wald statistics with unique exact Ensembl-to-Entrez mapping; no DEG or fold-change selection. Ambiguous mappings are excluded. MSigDB 2025.1.Hs Hallmark and Reactome; mapped gene-set sizes 15–500; fgseaMultilevel 1.32.4, `gseaParam=1`, `scoreType="std"`, `eps=0`, `sampleSize=101`, `nPermSimple=10000`, `nproc=1`, seed 20260914. Annotation: org.Hs.eg.db 3.20.0. Original `FDR_global` is retained, including the original SSc test in its correction family. It is not recomputed over the displayed rows.
- F: significant Hallmark/Reactome pathways (global FDR <0.05), separately by NES sign, are grouped into connected components using full-set Jaccard similarity ≥0.5. Representative: lowest global FDR, then P value, then name. Select the top six representatives per sign per collection by FDR; the stored table retains their final plotted order. This is retrospective display-level redundancy reduction, not preregistration or new enrichment testing. Cluster membership is in `metadata/pathway_display_clusters.json`.
- G heatmap: seven fixed external Hallmark response programs; custom ssGSEA-compatible implementation on the 13,215 uniquely Entrez-mapped VST genes. Within-sample average ranks are cast to integers, weighted by exponent 0.25, with integrated hit-minus-miss walk and no across-sample range normalization. GSVA was not executed for these saved scores. Row z-scores use the mean and sample SD (ddof=1) across six samples. Display limits are −2 to +2; stored values are not clipped.
- G dots: original NES; red positive, blue negative; filled global FDR <0.05, open otherwise. These are **GSEA statistics**, not P values from sample-score comparisons. Sample-score inference used exact two-sided label permutation over 20 assignments and BH across all 25 tested signatures; the complete correction family is included. No signature was removed for a nonsignificant or opposite result.

## Reproduce the data checks and plots

Use Python 3.11 or later and install `requirements.txt`, then run from this folder:

```sh
python -m pip install -r requirements.txt
python scripts/validate_data.py
python scripts/replot_EFG.py
```

The replot script draws three independent PNG panels from the frozen source tables. It reproduces data encodings, not the exact Illustrator composite layout. All text is Arial Bold 10 pt at saved physical size; Arial Bold must be installed (or use `--font /path/to/arial-bold.ttf`). The script stops rather than substituting fonts. Third-party font files are not distributed. The scripts do not perform a counts-to-DESeq2/GSEA analysis and are not presented as a complete upstream analytical pipeline.

Complete third-party gene-set libraries are not redistributed in this figure-source package. The original MSigDB release and set identifiers are recorded above and in the tables; rerunning upstream enrichment or signature scoring requires the original reference resources. Replotting saved results does not require these libraries.

## Audit and upload

`metadata/file_provenance.csv` maps byte-identical copies to original project-relative names. `SHA256SUMS.txt` covers every distributed file except itself. `metadata/validation_report.json` records checks against frozen upstream results. Missing values are `NA` or empty cells, not zero. All CSV files are UTF-8. No licence, DOI, sequencing accession or public repository URL is asserted by this package.

Upload the **contents of this folder** to your intended GitHub repository, retaining the relative directory structure. The accompanying ZIP is an alternative transfer/archive format; uploading both the ZIP and its extracted contents is unnecessary. No repository upload has been performed.
