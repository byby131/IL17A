# Figure 5A-E: adipogenesis RNA-seq source data

This data package supports Figure 5A-E of the IL-17A/adipogenesis study. It contains the original input count matrix, saved differential-expression results, panel-specific expression matrices, and enrichment results. The release was assembled on 2026-10-24 from existing analysis outputs; no differential-expression model or enrichment analysis was rerun.

**Code update:** This expanded copy adds a plotting implementation and an independent DESeq2 recomputation check. Read [START_HERE_CODE.md](START_HERE_CODE.md) before using or uploading the code. The original executable figure/GSEA/KEGG scripts were not present in the saved delivery; newly written scripts are explicitly labelled as reconstructions and do not replace the frozen manuscript statistics.

## Start here

| Current panel | Data directory | Content |
|---|---|---|
| Figure 5A | [A_volcano](data/A_volcano) | 17,887 tested entries; 17,882 plotted points; 10 labels |
| Figure 5B | [B_DEG_heatmap](data/B_DEG_heatmap) | All 1,823 DEGs: 897 upregulated, 926 downregulated |
| Figure 5C | [C_functional_program_heatmap](data/C_functional_program_heatmap) | 24 genes representing five functional programs |
| Figure 5D | [D_adipogenic_TF_heatmap](data/D_adipogenic_TF_heatmap) | All 14 selected adipogenesis-related transcription factors |
| Figure 5E | [E_KEGG](data/E_KEGG) | Top 10 significant KEGG pathways from each DEG direction; complete ORA and background |

The historical source folder used E for the functional heatmap, F for the TF heatmap, and C for KEGG. This release uses the current figure numbering. [panel_mapping.csv](metadata/panel_mapping.csv) preserves the mapping. All 43 original data CSVs were copied byte-for-byte into the new directory structure. Historical annotations inside copied files retain their original wording.

## Experiment and samples

Human ADSCs were sampled at day 7 of adipogenic induction (3 induction days followed by 4 maintenance days). The comparison is MDI + IL-17A (50 ng/mL) versus MDI + vehicle, with three independent biological replicates per group. Samples are ordered NC_1, NC_2, NC_3, IL17_1, IL17_2, IL17_3. In the current figure NC denotes Vehicle. The reference condition is vehicle; positive log2FC indicates higher expression in IL-17A.

The saved metadata do not establish the donor count or treatment/replenishment schedule. This day-7 dataset is distinct from the ADSC 0/3/12/24-h time course and from the HFF-1 conditioned-medium dataset.

## Raw input, processed data, and source data

- [original_counts_input.csv](data/shared/original_counts_input.csv): complete original input, 77,775 rows and six sample columns. The first column contains the original gene identifiers; its header is preserved. This is an input count matrix, not raw sequencing reads.
- [sample_metadata.csv](data/shared/sample_metadata.csv): original sample-to-column mapping, group, replicate type, and available experimental metadata.
- [DESeq2_size_factors.csv](data/shared/DESeq2_size_factors.csv): fitted normalization factors.
- Each panel folder retains the input-count subsets, normalized counts, VST expression, plotted values, and relevant statistics.
- [supporting_analysis](supporting_analysis): full tested-gene normalized/VST matrices (including the five entries omitted from the volcano), filtering and mapping audits, complete GSEA results/ranking/leading edges, functional selection audits, and KEGG mapping/set-size audits. These are additional copies from the same archived run.
- [metadata](metadata): saved analysis configuration, selection rules, original R session information, and current panel/sample mappings.
- [provenance/file_origin.csv](provenance/file_origin.csv): source-to-release mapping and hashes for copied files. `user_source/` and `archived_run/` denote source categories, not directories required at runtime.

FASTQ, BAM, quant.sf, and original alignment/quantification references are not part of this package. The upstream quantification software and original annotation version for this count matrix were not recorded in the saved analysis. The package supports inspection of count-level and figure-level results; it is not a complete sequencing pipeline.

## Recorded methods

The original DESeq2 analysis used a negative-binomial model with design `~ condition`, Wald testing, ratio size factors, and unshrunk log2FC estimates. Entries with counts >=10 in at least three samples were retained (17,887). Independent filtering was disabled. Gene-level P values were adjusted by Benjamini-Hochberg (BH); the DEG definition was `padj < 0.05`, without an additional fold-change cutoff. Five Cook-excluded entries have NA P/FDR. Package versions and the exact saved settings are in `metadata/`.

**A:** x = log2FC, y = -log10(BH-FDR). Display x values are clipped to [-8,8] and display y to 60. For a recorded zero FDR, the stored display transformation uses the smallest positive normal float before clipping; this is not a replacement inferential P value.

**B-D:** cells show row Z scores of saved DESeq2 VST expression, centered over all six samples and divided by sample SD (ddof=1). The color saturation range is [-2.5,2.5] for B and [-2,2] for C/D. Unclipped values are retained.

**C:** programs were drawn from the complete Hallmark/Reactome GSEA (MSigDB 2025.1.Hs; 1,061 tested unique sets; effective set size 15-500), using finite Wald statistics from uniquely mapped tested genes, including non-DEGs. The saved analysis used fgseaMultilevel (`scoreType=std`, `gseaParam=1`, `eps=0`, `sampleSize=101`, `nPermSimple=10000`, seed 20260913). Program selection used global BH-FDR <0.05; ordering was FDR, descending |NES|, then set ID. Within each functional domain, Jaccard >=0.5 removed redundant sets; at most two programs were retained. Up to six true leading-edge genes per program were selected by descending |Wald statistic|, descending |log2FC|, detectable sample count, and gene ID, followed by cross-program deduplication. The heatmap has 24 genes. Its FDR column is gene-level DESeq2 BH-FDR, not program-level GSEA FDR. The P1-P5 color strips identify programs.

**D:** all 14 discussion-selected TFs are retained, including nonsignificant genes. A single star appended to a gene name means gene-level BH-FDR <0.05; star count does not encode further significance tiers.

**E:** KEGG over-representation analysis used a one-sided Fisher exact/hypergeometric upper-tail test. The background comprises 6,002 tested, uniquely mapped, KEGG-annotated genes. Mapped foreground sizes are 486 up and 476 down. All 329 eligible pathways (15-500 background genes), including zero overlaps, were tested in both directions; the 658 direction-specific tests share one BH family. The separate `all` foreground has its own 329-test family and is retained in the complete results, but is not plotted. Each direction displays the top 10 pathways with FDR <0.05, ordered by FDR, descending fold enrichment, then pathway ID. Bar length = -log10(FDR). Direction refers to the input DEG list. The saved KEGG reference manifest records annotation download dates and checksums; database source files are not bundled.

## Zero, missing values, stars, and identifiers

- `pvalue`: unadjusted test P value; `padj`: gene-level BH-adjusted P value (FDR).
- `0*` in panel C means numerical underflow. The stored numerical value is 0; it must not be interpreted as an exact zero probability or as a conventional significance star. Six displayed genes have this annotation; seven genes in the full result have zero P/FDR.
- `NA` denotes missing/unavailable values and must not be converted to 0. In sample metadata it can denote an unrecorded attribute; in statistical columns the accompanying `na_reason` explains the cause.
- Mixed gene-symbol/Ensembl identifiers are preserved. Ensembl version stripping and exact one-to-one Entrez mapping belong to the archived enrichment workflow; no aliases were added during packaging.
- `hit_genes` and `hit_entrez` are separate set lists. Do not pair them by string position. Use the explicit membership table and unique identifier map.
- The copied functional row annotations include historical TGF-response-set membership and other fields that are not displayed in Figure 5C. These are annotations, not new experimental results.
- CSV encoding is UTF-8, sometimes with BOM as in the originals. Import gene identifiers as text when using spreadsheet software.

See [DATA_DICTIONARY.md](DATA_DICTIONARY.md), [file_catalog.csv](metadata/file_catalog.csv), and [CSV_COLUMNS.md](metadata/CSV_COLUMNS.md) for file and field details.

## Validation

The original-folder audit passed 155 checks: 43 original SHA-256 hashes, lineage to archived analysis tables, raw-count subsets, sample order, normalization, VST row Z scores, color clipping, gene labels and FDR text, pathway order, bar lengths, and membership. Visible labels/order were checked against the user-supplied final Figure 5. This is a data/content verification, not a claim of pixel-level raster identity. See [AUDIT_REPORT_zh.md](AUDIT_REPORT_zh.md) and [original_source_check.json](validation/original_source_check.json).

To verify the distributed package:

```sh
python -m pip install -r code/requirements.txt
python code/verify_source_data.py --root . --output validation/current_check.json
```

The script checks hashes and saved data consistency. It does not run DESeq2, GSEA, or ORA. The original R environment is preserved separately from the Python verification environment. `.gitattributes` disables automatic newline conversion so the byte-level manifest remains valid.

## Publication metadata

This is a locally prepared release; no GitHub upload or public deposit has been made by this preparation step. Before public release, supply the dataset authors/contact, repository URL, associated manuscript identifier when available, and the author-approved data/code licences. No licence or accession number has been invented or granted by this package. Third-party annotations retain their source attribution and applicable terms. Record the final release tag/commit in the manuscript's availability statement.
