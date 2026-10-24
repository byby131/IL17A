# Supplementary Figure 2 A-C: GitHub data package

Prepared 2026-10-24 from the verified `SupplementaryFigure3_SourceData_20261009` archive for the author's final **Supplementary Figure 2**. The previous **S3B, S3C, S3D correspond to current S2A, S2B, S2C**. The old reserved fluorescence panel A is not part of the final figure. No new DESeq2 fitting, GSEA or clustering was performed.

## Upload contents

Extract the ZIP and upload this folder's contents, preserving the directory structure. `PANEL_INDEX.csv` identifies every plotted component. `PROVENANCE.csv` records source paths within the old package and copying/derivation operations. `MANIFEST_SHA256.csv` verifies all packaged files except itself. Nothing has been uploaded to a remote repository.

- `01_input_counts`: original unnormalized Salmon estimated transcript counts (232,117 transcripts x 12 samples), plus sample metadata/size factors. These counts are byte-identical to the input in the Figure4 GitHub package prepared on the same date. TPM is not used.
- `02_panel_source_data`: current A/B/C source tables, including every heatmap cell, gene point, composition count, leading-edge member and module trajectory observation.
- `03_complete_statistics`: all five DESeq2 contrasts, LRT, all 5,120 Hallmark/Reactome GSEA results, ranked backgrounds, pathway membership/selection records, and secretome eligibility/stability evidence. Historical duplicates are retained where they document an upstream stage.
- `04_frozen_inputs_and_models`: aggregated gene counts, annotation, saved DESeq2 objects, all-gene VST/normalized counts, original server Group1-7 membership and Group5 mappings.
- `05_reference_snapshots`: exact Hallmark/Reactome GMT and UniProt secreted-protein snapshots. The separately obtained GENCODE v36 GTF remains in the original archive; its verified hash is recorded here without duplicating that 44.5 MB reference file.
- `06_code_and_environment`: portable data validator plus unmodified historical analysis/plotting code and environment records. Historical code retains old paths and S3 labels and requires adaptation; it is provenance, not a tested one-command exporter of the current edited S2 layout.
- `07_provenance_and_QC`: current panel configuration, original integrity checks, screenshot checks and earlier numerical audits.
- `08_FASTQ_inventory`: 24 paired-end FASTQ filenames, sample mapping, sizes and supplier MD5 values. Reads are not included in this GitHub package.

## Experiment and raw data

Human ADSCs from **one donor, three independent experiments**. Each experiment includes 0, 3, 12 and 24 h; same-number experimental replicates are matched across time. 0 h is untreated / pre-IL-17A baseline; later samples were collected during 50 ng/mL IL-17A exposure. `Donor_1` is a study-local anonymized label. FASTQ R1/R2 means read mates, whereas plot R1/R2/R3 means independent experiments.

Salmon estimated counts can be fractional and are not raw sequencing reads. Exact versioned transcript IDs were mapped to reconstructed GENCODE v36 CHR GRCh38.p13, summed to 60,660 genes and rounded for DESeq2. Count >=10 in >=3 samples retained 14,425 genes. Design: `~ replicate + time`, 0 h reference; LRT reduced model: `~ replicate`. VST uses `blind=FALSE`. Original quant.sf, Salmon index/run records and sample-specific effective lengths are unavailable; the recovered FASTQ filenames do not by themselves prove the missing FASTQ-to-Salmon processing history.

## Current panel definitions

**A.** All 22 displayed pathways and 4 comparisons (3/12/24 vs 0, and 24 vs 3), 88 cells. Color is NES clipped to [-3,3]. Text and stars use **padj_library_contrast**, BH within the complete eligible library and comparison: 49 Hallmark and 975 Reactome tests per comparison, five comparisons in the full tables (5,120 rows). Do not substitute other historical `padj` fields. The omitted display contrast 12 vs 3 remains in complete results. All 88 final-figure text labels were checked exactly, including scientific notation and stars.

**B.** Four fixed pathway representatives: Hallmark TGF-beta, Reactome BMP, Reactome TGF-beta/SMADs, Reactome TGFB family. Detectable denominators are 51,24,47,141; G5 memberships are 4,4,3,7. Top bars use counts/proportions at 3,12,24 h vs 0. Up/Down requires original genome-wide gene BH-FDR<0.05 and log2FC sign; NS means not significant. Middle panels include all 526 gene observations at 12 and 24 h, partitioned by frozen Group5 membership, with 16 boxes. These points are genes, not extra biological replicates; shared pathway members can appear in more than one pathway. Box summaries use linear quartiles and 1.5-IQR whiskers; all scatter values remain available, including outliers. Bottom annotations use NES sign and G5 / leading-edge counts, with all 207 corresponding gene records provided. All 12 displayed pathway/contrast leading-edge results have FDR>=0.05 and are descriptive.

**C.** The full 500 dynamic secreted-protein-encoding genes, not Top500: frozen UniProt-reviewed human secreted annotation, genome-wide LRT FDR<0.05 and VST availability. There are 820 secreted-annotated tested genes; all 500 eligible genes enter clustering. Gene-wise Z=(VST-row mean)/sample SD across all 12 samples. Frozen Ward/Euclidean clustering on four time means gives S1=337 and S2=163. The heatmap uses its saved row order and a [-2.5,2.5] color range (visible ticks -2,0,2). Module scores average all member-gene Z values per sample. The curve and shaded band are mean +/- sample SD across **three independent experiments**, not SEM or 95% CI. All 6,000 gene/sample observations, raw aggregated estimates, integer model counts, normalized counts, VST and Z are supplied. CXCL8/IL6 are in S1; CSF3/CXCL1 are in S2.

## Verify

With Python, pandas and numpy available, run `python 06_code_and_environment/validate_source_data.py`. The validator checks SHA256, full-library BH values, the screenshot's 88 labels, all B composition and leading-edge counts, all plotted gene effects, and C's counts/normalization/VST-Z/module means/SD. It does not rerun Monte Carlo GSEA, refit DESeq2, recluster or read 36 GB of FASTQ.

No time-matched untreated controls were available; changes during IL-17A exposure cannot be fully separated from culture time. The original server clustering parameters and original quantification records remain incomplete. The UniProt release number was not recorded; the exact frozen file/hash is retained. No accession, DOI or license has been invented.

## FASTQ and repository scope

All 24 FASTQ files were found at their recorded G-drive locations on 2026-10-24, with file sizes matching the original inventory: **35,966,545,938 bytes**. Supplier MD5 values are copied, not independently recomputed. The GitHub manifest omits workstation-specific absolute paths; the separate local audit retains those paths.

Upload the extracted analysis/source files to GitHub. Individual FASTQ files exceed ordinary GitHub repository file limits; sequencing reads should be archived through the relevant sequencing-data repository workflow, with a real accession then linked from this repository. GitHub browser uploads allow individual files up to 25 MiB and ordinary Git blocks files over 100 MiB ([GitHub documentation](https://docs.github.com/en/repositories/working-with-files/managing-large-files/about-large-files-on-github), checked 2026-10-24). See the [NCBI sequencing/functional-genomics submission routes](https://www.ncbi.nlm.nih.gov/sra/docs/submit/) for FASTQ deposition. This package does not create a sequencing submission or fill missing repository metadata.
