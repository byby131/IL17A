# Figure 4 A-D: input counts and figure source data

Prepared 2026-10-24 for the user's current Figure 4 screenshot. All primary data were taken from the SHA256-verified `Figure4_submission_package_20260922`.

## Correct panel mapping

Current A-D are four components of the old package's composite **Figure4C**. The old package's A (PCA), B (Hallmark GSEA), and D (persistent inflammatory effectors) are different panels and are not included as current A-D data. See `PANEL_INDEX.csv` and `核对报告_中文.md`.

- A: frozen 5,000-gene server heatmap, 12 samples, Groups 1-7.
- B: Group 1 / Group 5 module trajectories, 24 sample values and 8 means with 95% t intervals.
- C: 4 Group 1 and 5 Group 5 KEGG enrichment bars; bar lengths are **fold enrichment**, not -log10(P/FDR).
- D: 11 Group 5 KEGG TGF-beta overlap genes, 12 sample values and 4 time-point means.

## What to upload

Upload the contents of this folder, preserving the directory structure. `MANIFEST_SHA256.csv` verifies all packaged files except itself. The adjacent ZIP and `.zip.sha256` are transport and integrity files; a `.sha256` file alone contains no expression data. No repository has been created or published by this preparation step.

`01_input_counts/ADSC_IL17A_timecourse_RNAseq_EXPRESSION_MATRIX_COUNTS.CSV` is the original available **analysis input**: 232,117 versioned transcript IDs x 12 samples of unnormalized Salmon estimated counts. TPM was not used. These estimated counts are not raw sequencing reads and need not be integers. FASTQ and original Salmon `quant.sf` files are not included or found in the inspected package/workspace.

`02_frozen_server_exports/` contains the original iDEP heatmap and enrichment exports. These are processed server outputs, not raw counts. `03_panel_source_data/` contains current-panel data. `04_processed_support/` contains aggregated gene counts, annotation, the saved DESeq2 object and fixed VST-derived values. The historical `sender` filename refers to this ADSC time course, not adipogenesis or HFF-1 RNA.

## Experimental design and authoritative metadata

One human ADSC donor; three independent experiments, each including 0, 3, 12 and 24 h. Same-number replicates across time belong to the same experiment. 0 h is untreated / pre-IL-17A baseline; subsequent samples were exposed to 50 ng/mL IL-17A. The model is `~ replicate + time`, with 0 h as reference. `donor_1` is a study-local label, not a personal or original donor identifier.

Use `01_input_counts/sample_metadata.csv`. The historical metadata is retained under `06_provenance/` for traceability; its generic IL-17A condition label at 0 h is superseded by the explicit untreated baseline here. The one-donor clarification follows the author's 2026-10-07 confirmation. Sample IDs and all expression values are unchanged.

## Methods and interpretation of values

A uses the frozen server-centered expression matrix. No reclustering was performed. Within each group, display rows are sorted by mean(12,24 h) minus mean(0,3 h). Group sizes are 333, 731, 811, 898, 1043, 618 and 566. The symmetric color limit is the 99th percentile of absolute values (2.3235826, displayed as 2.3).

B uses DESeq2 `vst(blind=FALSE)`, followed by gene-wise Z scoring over 12 samples (sample SD, denominator n-1), then averaging genes within each frozen module for each sample. Group 1 has 333 frozen members and 147 VST matches; Group 5 has 1,043 frozen members and 800 VST matches. Lines show the mean of three independent experiments; bands use mean +/- t(0.975, df=2) * SD/sqrt(3). Do not substitute A's server matrix for B's VST matrix.

C preserves the original iDEP data113 KEGG results, not a later MSigDB KEGG_LEGACY reanalysis. For Group 5 TGF-beta: query=229 KEGG-annotated genes, background=5487, pathway-in-background=87, overlap=11, one-sided hypergeometric P=0.0008926000825948, BH FDR=0.0428775059034482 (exported 0.0429), fold enrichment=3.0295136274657435. BH scope was 268 Group5-hit pathways. Only the available exported eight Group5 rows and their audit are included; a full 268-row result table and the complete iDEP database/background membership are not present in the old package. Group1 values remain as exported. See `05_KEGG_audit/` for the reconstruction and database identification.

D uses means of three replicates from the **server display matrix**, not VST or counts, to retain all 11 original overlap genes, including NOG and INHBE, which are absent from the rebuilt VST matrix. Display order: ID1, ID2, ID3, ACVR2B, ID4, SMAD9, PITX2, NOG, SMAD6, RGMA, INHBE. Symmetric color limit is 1.1226162 (displayed as 1.1). The four-time-by-eleven-gene CSV directly matches the current orientation.

## Reproduction scope and provenance

Run `python code/validate_source_data.py` after installing pandas and numpy. This validates hashes, dimensions, module scores/CI, gene means and C's bar values. It does not refit DESeq2 or rerun the historical server. The saved `sender_dds.rds` and original session information allow follow-up model work in a compatible R environment. `PROVENANCE.csv` maps each copied/derived table to the old package path. Copied files are byte-identical; derived tables only select, label, reorder, or transpose existing values.

The counts-to-gene reconstruction used exact versioned transcript matching to GENCODE v36 CHR GRCh38.p13, followed by rounding and count >=10 in at least 3 samples (14,425 retained genes). The subsequently obtained reference GTF is in the full old package, not this focused data subset; SHA256: `81359865bab048361cbeaf91c9441aac4c700b63978e6691d579d78bf2322d30`. The original quantification GTF/tx2gene and sample-specific effective-length information were unavailable. Frozen gene counts, annotation and the DESeq2 object are included; this is not a complete FASTQ-to-figure pipeline.

There were no time-matched untreated controls, so culture-time effects cannot be fully separated from IL-17A exposure. Frozen clusters and their trajectories derive from the same experiment. KEGG enrichment does not by itself demonstrate activation of the full canonical pathway. No accession, DOI or license has been invented or assigned.

## Visual correspondence

The original combined PDF and current screenshot share the same heatmap patterns, module trajectories, 9 KEGG terms/bar values and 11-gene ordering. The screenshot separates the old composite C into A-D and changes titles, colors, gridlines, legends and layout. Correspondence was assessed using visible features and frozen numeric source tables, not the original editable assembly file. This package delivers data, not a modified figure.
