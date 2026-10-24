# Data dictionary and current-panel naming

Original file contents are preserved byte-for-byte where copied. Current filenames relabel S3B→A, S3C→B and S3D→C; internal historical field names are retained. `PROVENANCE.csv` provides every original path.

| Field / file | Definition |
|---|---|
| Counts `Name` | Versioned ENST transcript identifier; not a gene row |
| `S50_<time>_<replicate>` | Sample columns; experimental replicate is matched across time |
| `gene_id` / `unversioned` / `ENTREZID` | Ensembl versioned / unversioned / Entrez gene identifiers; mapping is supplied |
| `sizeFactor` | Saved DESeq2 sample normalization factor |
| raw aggregated Salmon counts | Fractional estimates summed by gene, before rounding/normalization |
| model integer counts | Rounded saved DESeq2 input |
| normalized counts | Integer model counts divided by saved sample size factor |
| VST / gene Z | Variance-stabilized expression / gene-wise Z across 12 samples |
| `padj` in gene DE | Original genome-wide gene BH-FDR |
| `padj_library_contrast` in A GSEA | Authoritative complete-library-within-contrast GSEA BH-FDR |
| other GSEA `padj*` | Historical correction scopes; not interchangeable with A display FDR |
| `NES`, `ES`, `pval`, `size` | GSEA normalized score, enrichment score, nominal P, ranked-overlap set size |
| `in_Group5` | Frozen server Group5 membership; different from C's secretome S1/S2 |
| `total_detectable` | Number of testable mapped members for that pathway/contrast |
| `count`, `proportion` | Number in Up/Down/NS category and fraction of detectable denominator |
| `figure_column_order` | Four B pathways ordered left to right (1-4) |
| `gene_BH_FDR`, `direction` | Full-gene DE FDR and direction; no added fold-change threshold |
| `n_Group5`, `n_leading_edge`, `fraction_Group5` | Leading-edge G5 numerator, all leading-edge denominator and fraction |
| `pathway_BH_FDR` | Same library-wide GSEA FDR used by A |
| `secreted_module` / C `module` | S1 or S2; not Group1 or Group2 |
| `plotted_row_1based`, `leaf_order` | Final heatmap row and original dendrogram leaf coordinate |
| `module_score` | Mean gene Z over all member genes for one sample |
| `SD`, `shaded_lower`, `shaded_upper` | Sample SD across 3 experiments and mean±SD in C |
| FASTQ `read_mate` | Sequencing R1/R2 mate, not an experimental replicate |
| FASTQ `supplier_md5` | Original supplier checksum, not independently recomputed |

NA/blank means unavailable, not zero. Extreme nominal/FDR zero values may reflect stored numerical underflow. Original double precision is retained. In the boxplot data each point is a gene; do not interpret gene count as biological sample size.
