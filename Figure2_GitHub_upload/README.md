# Figure 2 analysis and plotting code (panels A–I)

This folder contains the R scripts and small, frozen analysis inputs corresponding to the current nine-panel Figure 2 layout shown in `reference/Figure2_current_layout.png`. Run commands from **this folder's root**. The panel letters here refer to the current manuscript layout; older script/output names retain their original development labels.

| Current panel | Content | Primary script(s) | Existing output stem |
|---|---|---|---|
| A | Global skin atlas | `make_figure2B_new_strategy.R` | `Figure2B_global_UMAP_new_strategy` |
| B | CD4 T-cell states | `cluster_and_profile_tcells.R`; `make_figure2CDE_new_strategy.R` | `Figure2C_CD4_state_UMAP` |
| C | Th17-like/Treg marker dot plot | `make_figure2CDE_new_strategy.R` | `Figure2D_Th17_Treg_marker_dotplot` |
| D | Donor-level CD4 proportions | `make_figure2CDE_new_strategy.R` | `Figure2E_donor_abundance` |
| E | Healthy stromal reference and LoSc mapping | `prepare_healthy_stromal_reference.R`; `map_and_plot_figure2F.R` | `Figure2F_healthy_reference_LoSc_mapping` |
| F | Stromal identity dot plot | `make_figure2G_strict.R`; `make_figure2G_compact_labels.R` | `Figure2G_compact_stromal_identity_dotplot` |
| G | Q8/Q11 Hallmark pathway NES bars | `make_figure2H_Q8_Q11_revision.R`; `make_figure2H_Q8_vs_Q11_direct.R`; `make_figure2H_disease_up_two_layouts.R`; `make_figure2H_compact_horizontal_bar.R` | `Figure2H_compact_horizontal_bar` |
| H | Donor-level CD4-to-ASPC CellChat-informed strength | `analyze_CellChat_CD4_Q8Q11_vs_Healthy.R`; Q11/Q8 donor scripts; `make_CellChat_Q8_Q11_vertical_combined_bar.R` | `CellChat_Q11_Q8_vs_Healthy_vertical_combined` |
| I | Matched-cell CellChat networks | `analyze_CellChat_CD4_Q8Q11_vs_Healthy.R`; `make_CellChat_A_thin_network.R` | `CellChat_A_emphasized_Q11_interaction_networks` |

`scripts/figure2_ghi_common.R` supplies shared H5Seurat readers and expression helpers. `scripts/audit_figure2H_Q8_Q11.R` creates the frozen cell selection and pseudobulk quality table used by G, H and I. `scripts/extract_h5_metadata.R` extracts the shared atlas metadata.

## Inputs and repository scope

The processed public GSE264508 H5Seurat file is **not included** because it is approximately 1.9 GB. Obtain the processed dataset separately and place it at `data/GSE264508/GSE264508_data.all.H5Seurat`. The metadata, smaller intermediate tables, selected pathway results, donor-level CellChat source data, and three matched-cell CellChat objects are included at their expected relative paths. `FILE_MANIFEST.csv` records SHA-256 digests for all copied files.

The 97 MB LoSc query count cache is also omitted; `make_figure2G_strict.R` regenerates it from the H5Seurat file. No clinical donor table is included in this Figure 2 code package. The included donor labels are those used in the public dataset's analysis tables.

For exact statistical reruns, use the package versions recorded in the included `*sessionInfo.txt` files. Key versions used in the saved results include R packages Seurat 5.5.0, DESeq2 1.46.0, fgsea 1.32.4, msigdbr 26.1.0 and CellChat 2.2.0.9001. Required plotting packages are loaded at the top of each script.

## Fast regeneration from included frozen inputs

These commands regenerate the core plot files without repeating upstream clustering or differential analyses. They do **not** require the 1.9 GB H5Seurat file, except where stated.

```bash
Rscript scripts/make_figure2B_new_strategy.R \
  work/gse264508_audit/metadata_umap.csv.gz outputs

Rscript scripts/make_figure2G_compact_labels.R
Rscript scripts/make_figure2H_compact_horizontal_bar.R
Rscript scripts/make_CellChat_Q8_Q11_vertical_combined_bar.R
Rscript scripts/make_CellChat_A_thin_network.R
```

Panels B–E and the full upstream analyses use the H5Seurat file. The complete command order is below.

## Full analysis order

```bash
Rscript scripts/extract_h5_metadata.R \
  data/GSE264508/GSE264508_data.all.H5Seurat work/gse264508_audit

Rscript scripts/make_figure2B_new_strategy.R \
  work/gse264508_audit/metadata_umap.csv.gz outputs

Rscript scripts/cluster_and_profile_tcells.R \
  data/GSE264508/GSE264508_data.all.H5Seurat \
  work/gse264508_audit/metadata_umap.csv.gz work/figure2_cd4_new

Rscript scripts/make_figure2CDE_new_strategy.R \
  work/figure2_cd4_new/Tcell_secondary_clustering_and_expression.csv.gz \
  data/GSE264508/GSE264508_data.all.H5Seurat outputs

Rscript scripts/prepare_healthy_stromal_reference.R \
  data/GSE264508/GSE264508_data.all.H5Seurat \
  work/gse264508_audit/metadata_umap.csv.gz work/figure2_stromal_new

Rscript scripts/map_and_plot_figure2F.R \
  work/figure2_stromal_new/stromal_reference_query_expression.csv.gz \
  data/GSE264508/GSE264508_data.all.H5Seurat outputs

Rscript scripts/make_figure2G_strict.R \
  data/GSE264508/GSE264508_data.all.H5Seurat \
  outputs/Figure2F_healthy_reference_cell_metadata.csv.gz \
  outputs/Figure2F_LoSc_query_mapping_metadata.csv.gz \
  outputs/Figure2G_strict

Rscript scripts/make_figure2G_compact_labels.R

Rscript scripts/audit_figure2H_Q8_Q11.R \
  data/GSE264508/GSE264508_data.all.H5Seurat \
  outputs/Figure2G_strict/Figure2G_cell_metadata.csv.gz \
  outputs/Figure2H_Q8Q11_revision/audit

Rscript scripts/make_figure2H_Q8_Q11_revision.R \
  data/GSE264508/GSE264508_data.all.H5Seurat \
  outputs/Figure2H_Q8Q11_revision/audit/Figure2H_Q8_Q11_frozen_cell_metadata.csv.gz \
  outputs/Figure2H_Q8Q11_revision/audit/Figure2H_donor_group_coverage_and_quality.csv \
  outputs/Figure2H_Q8Q11_revision/delivery \
  outputs/Figure2F_LoSc_query_mapping_metadata.csv.gz

Rscript scripts/make_figure2H_Q8_vs_Q11_direct.R \
  outputs/Figure2H_Q8Q11_revision/delivery/source_data/Figure2H_Q8_pseudobulk_raw_counts.rds \
  outputs/Figure2H_Q8Q11_revision/delivery/source_data/Figure2H_Q8_model_donors.csv \
  outputs/Figure2H_Q8Q11_revision/delivery/source_data/Figure2H_Q11_pseudobulk_raw_counts.rds \
  outputs/Figure2H_Q8Q11_revision/delivery/source_data/Figure2H_Q11_model_donors.csv \
  outputs/Figure2H_Q8Q11_revision/audit/Figure2H_donor_group_coverage_and_quality.csv \
  outputs/Figure2H_Q8Q11_revision/Q8_vs_Q11_direct/delivery

Rscript scripts/make_figure2H_disease_up_two_layouts.R
Rscript scripts/make_figure2H_compact_horizontal_bar.R

Rscript scripts/analyze_CellChat_CD4_Q8Q11_vs_Healthy.R
Rscript scripts/make_CellChat_Q11_donor_significance_bar.R
Rscript scripts/make_CellChat_Q8_donor_significance_bar.R
Rscript scripts/make_CellChat_Q8_Q11_vertical_combined_bar.R
Rscript scripts/make_CellChat_A_thin_network.R
```

Each script retains the original analysis settings and writes its outputs under `outputs/`. Running the full sequence can overwrite the included frozen tables; use a separate working copy if preservation of the submitted snapshot is needed.

## Statistical and layout notes

- G uses Hallmark GSEA on donor pseudobulk comparisons. Q8 and Q11 each show eight positive-NES pathways ranked by their original q values; a displayed pathway is **not necessarily significant**. The direct Q8–Q11 panel shows four pathways per direction. Only q < 0.05 values are labelled in the compact bar plot. Frozen full pathway and DESeq2 results are provided.
- H shows the mean and SEM of donor-level CellChat-informed interaction scores, with one point per biological donor. The two-sided exact donor-label permutation P values were corrected by BH across five CD4 states **within each Q11 or Q8 comparison**. Q8 and Q11 panels used different fixed ligand–receptor universes; their absolute heights should not be interpreted as a direct cross-panel test.
- I uses cell-number-matched pooled CellChat models and a shared edge-width scale. Healthy total strength 5.92 and Q11 total strength 12.03 are descriptive aggregates, not donor-level hypothesis tests.
- The supplied reference image shows the journal composite. Its G and H components were repositioned during editorial assembly: the original G plotting script exports the three comparisons horizontally, and the original H script exports Q11 above Q8. The source code and data reproduce the plotted values and component plots, while the exact composite placement and the red emphasis added to selected labels are not encoded as a separate assembly script in the available files.

The archived development drafts, unrelated figures, large H5Seurat file, local R package library, and temporary caches were deliberately excluded from this upload folder.
