# Reproduction and scope QA

- Current Figure 5 computational panels: A–E. The experimental/schematic panels F–K are outside this bioinformatics code package.
- Frozen source-data verification: `PASS`, 205 checks. Counts: 77,775 input rows, 17,887 modeled genes, 17,882 volcano points, 897 upregulated and 926 downregulated DEGs, 24 C genes, 14 D genes, 20 E pathways and 6,002 KEGG-background genes.
- Primary DESeq2 model was independently reimplemented from the saved counts and documented configuration. It retained 17,887 genes, had no missing-value discordance with the frozen result, and produced only small cross-platform numeric differences (see `recomputed_analysis/comparison_to_frozen_statistics.csv`). This check did not replace archived P values or FDRs.
- The A–E plotting script parsed and rendered all five panels as PDF, SVG, PNG and 600-dpi TIFF from the frozen plotting tables. No simulated observations, reselected genes, re-ranked pathways or new inferential results were used.
- The reconstructed visual layout is **not** the historical manuscript artwork. In particular, labels and white space differ from the author-approved composite. Use the saved source data and code for audit/reuse, not the reconstruction as a drop-in figure replacement.
- Nature-figure static preflight found no failures. Its remaining warnings are: dynamic per-panel width values; the source contains the text `log10` in axis labels rather than an unguarded numeric log transform; and static syntax checking was followed by a successful R parse/runtime test.
- The original executable GSEA/KEGG/figure-assembly scripts were not among the preserved assets; do not describe these reconstructed scripts as those historical files.
