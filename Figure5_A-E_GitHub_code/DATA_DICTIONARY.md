# Data dictionary

## Shared matrix conventions

Rows are original gene identifiers; sample columns are NC_1, NC_2, NC_3, IL17_1, IL17_2, IL17_3. All values in copied CSVs retain the original serialized precision. Row/sample order is meaningful in plotting matrices. `NA` and blank fields represent missing/unrecorded values in context, never numerical zero.

| File/field | Meaning |
|---|---|
| raw_counts.csv; original_counts_input.csv | Original nonnegative integer counts; panel subsets versus complete input |
| normalized_counts.csv | Original counts divided by saved DESeq2 sample size factor |
| VST_expression.csv | Saved variance-stabilized expression, not original counts |
| heatmap_VST_rowZ.csv | Per-gene (VST - six-sample mean) / sample SD |
| heatmap_color_values_clipped.csv | Row Z after color-range saturation; not a new statistical result |
| gene_id | Original gene symbol or Ensembl identifier, kept as text |
| entrez_id | Unique mapped Entrez identifier, kept as text |
| sample_id; count_column | Analysis sample ID; matching column name in original input |
| condition | Internal model level: vehicle or IL17A |
| size_factor | Saved DESeq2 normalization factor |
| row_order; draw_order | One-based top-to-bottom heatmap row order; volcano point draw order |
| baseMean | Mean normalized count across the six samples |
| log2FC | Unshrunk log2(IL17A / vehicle) effect estimate |
| lfcSE | Standard error of log2FC |
| stat | DESeq2 Wald statistic |
| pvalue | Unadjusted P value |
| padj; DESeq2_padj | BH-adjusted gene-level P value; retained cross-check column from DESeq2 |
| log2FC_scale; shrinkage_method; log2FC_shrunk | Original effect-estimation metadata; no shrinkage, shrunk values unavailable |
| direction | up/down/not_significant/not_estimable for genes; enrichment-list or NES direction in enrichment tables |
| beta_converged; max_cooks; na_reason | Model diagnostics; maximum Cook distance; reason for unavailable statistics |
| sensitivity_abs_log2FC_ge_1 | Optional sensitivity flag; not the threshold used for Figure 5 |
| samples_count_ge_10; quality_flag | Count-detection support and saved quality annotation |
| plot_x_log2FC; plot_y_minus_log10_FDR | Display-clipped volcano coordinates |
| minus_log10_FDR | Unclipped display transformation with machine-normal floor for zero FDR |
| point_color; top_label; selection_source | Saved display color, label flag and original label-selection rule |
| primary_domain; domain | Functional domain assigned in the selection audit |
| all_contributing_programs | Semicolon-separated P1-P5 program membership; one gene may have multiple memberships |
| program_code; color; name | Program code, hex color, full original gene-set name |
| existing_readout | Previously discussed experimental-readout annotation |
| external_TGFB1_response_set_member | Historical membership annotation; not a measured TGF-beta response |
| observed_early_response | Unavailable historical field; the day-7 panel is not an early time-course analysis |

## GSEA / functional selection

| Field | Meaning |
|---|---|
| set_id; collection; all_source_set_ids | Unique gene-set ID, Hallmark/Reactome collection, source IDs merged by identical members |
| pval; software_padj | Original fgsea P value and software-adjusted result |
| ES; NES | Enrichment score; normalized enrichment score |
| size; leadingEdge | Effective tested set size; semicolon-separated leading-edge Entrez IDs |
| log2err | Original fgsea error estimate |
| FDR_global_HR | BH-adjusted value across all tested unique Hallmark + Reactome sets; used for selection |
| FDR_within_Hallmark; FDR_within_Reactome | Supplemental within-collection adjustments; not the primary selection threshold |
| rank_global | Rank by global FDR, ties assigned minimum rank |
| rank_in_program; selected; selection_reason; reason | Candidate order, display-selection flag, and saved selection/exclusion reason |
| comparison | IL17A_vs_vehicle in the leading-edge table |

## KEGG

| Field | Meaning |
|---|---|
| pathway_id; pathway_name; kegg_url | Human KEGG pathway identifier, saved name and source URL |
| direction | DEG foreground: up, down, or all (the all family is not shown in E) |
| Count | Number of foreground genes overlapping the pathway |
| foreground_size | Number of foreground genes represented in the KEGG background |
| pathway_size | Pathway members within the background |
| universe_size | Background gene count (6002) |
| GeneRatio; BgRatio | Count / foreground_size; pathway_size / universe_size |
| fold_enrichment | GeneRatio / BgRatio |
| odds_ratio | Fisher-test odds-ratio estimate |
| pvalue; FDR | One-sided ORA P value; BH-adjusted value for the stated correction family |
| correction_family | All eligible pathways jointly across up/down, or separate all-DEG family |
| hit_genes; hit_entrez | Independently sorted set lists; do not match positionally |
| category; subcategory | Saved KEGG hierarchy annotations |
| group; key | Historical plotting labels and direction/pathway key |
| display_order_within_direction | One-based pathway order within up/down panel |
| bar_length_minus_log10_FDR | Actual numerical horizontal bar length |
| effective_size; eligible | Size within background and eligibility under 15-500 rule |

## Provenance and file catalog

`metadata/file_catalog.csv` lists every scientific CSV payload with row/column counts and byte size. `metadata/CSV_COLUMNS.md` lists its exact column names. `provenance/file_origin.csv` gives copied-file origin and SHA-256. `provenance/original_manifest.csv` preserves the old directory names; the top-level `manifest.csv` describes this release. Reference manifests identify original database source files that are not bundled, rather than claiming they exist in this package.
