# Data dictionary

| Files / field | Meaning |
|---|---|
| Counts `Name` | Versioned ENST transcript ID; rows are not genes |
| `S50_<time>_<replicate>` | Same sample IDs across matrices; time in hours; replicate identifies one complete experiment |
| Counts sample columns | Unnormalized Salmon estimated counts; fractions are expected |
| `gene_id` in sender gene files | Versioned ENSG ID after GENCODE v36 aggregation |
| `ensembl_ID`, `gene_id` in module matrices | Gene identifiers; follow mapping CSV for version conversion |
| `symbol`, `gene_name` | Gene symbol |
| `cluster` / `group` | Frozen iDEP Group 1-7; not a treatment or donor |
| Server heatmap sample columns | Frozen processed server-centered expression, not raw counts |
| `display_late_delta` | Mean server value at 12/24 h minus mean at 0/3 h, for row ordering |
| VST Z matrix sample columns | Gene-wise Z of DESeq2 VST across 12 samples |
| `matched_vst` / `dds_gene_id` | Whether/how a frozen server gene maps into the saved DESeq2 model |
| `module_score` | Per-sample mean gene Z over matched genes in that module |
| `n_genes` | 147 for G1 and 800 for G5; not sample count |
| `mean_module_score`, `sem`, `ci95_low`, `ci95_high` | Mean, standard error and 95% t interval across 3 experiments at that time |
| `fold_enrichment` | Current C horizontal bar value |
| `FDR_exported` | Original server BH-adjusted P, displayed/exported precision |
| `BH_FDR_reconstructed` / `raw_P_reconstructed` | Group5 numerical audit values; blank for Group1, where not reconstructed here |
| `displayed_global_pathway_size` | Whole database pathway size; not the filtered-background pathway size |
| `overlap_ensembl` | Comma-separated overlapping genes |
| D wide columns `0`, `3`, `12`, `24` | Mean server display expression at that time over 3 samples |
| D `heatmap_4_times_by_11_genes.csv` | Rows=time; columns=genes, exactly ordered for current D |
| `study_donor_label` | Study-local label donor_1; all 12 samples are from the same donor |

Blank reconstructed-statistic cells mean not supplied/reconstructed, not zero. Original CSV names/fields may retain historical Figure4C labels; current mapping is in PANEL_INDEX.csv. Original numeric exports retain their stored precision.
