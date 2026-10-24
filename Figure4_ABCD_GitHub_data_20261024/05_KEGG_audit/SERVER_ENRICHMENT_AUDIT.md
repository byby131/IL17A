# Figure 4C Group 5 server enrichment audit

## Conclusion

The original server export is internally coherent and reproducible as an iDEP KEGG over-representation analysis. The complete frozen Group 5 contains 1,043 unique genes. iDEP then restricted that list to 229 genes with at least one KEGG annotation and used 5,487 filtered expressed, KEGG-annotated genes as the effective background. For KEGG hsa04350, 11 Group 5 genes overlapped 87 background genes in the pathway (the export displays the global database size of 107). This reproduces raw P = 0.0008926000825948277, fold enrichment = 3.02951362746574, and exported BH FDR = 4.29e-02.

The exported FDR values across all eight Group 5 rows uniquely reproduce when 268 KEGG pathways are corrected within Group 5. iDEP only tested pathways hit by at least one Group 5 gene, then displayed FDR < 0.1 and at most 10 rows. TGF-beta signaling had raw-P rank 5; the monotone BH value propagated from rank 7 is 0.0428775059034482, rendered by iDEP as 4.29e-02.

## TGF-beta result

- Pathway: KEGG hsa04350, TGF-beta signaling pathway
- Complete Group 5 membership: 1,043
- Effective KEGG-annotated Group 5 query: 229
- Effective filtered KEGG-annotated background: 5,487
- Global pathway size displayed by iDEP: 107
- Effective pathway size in the filtered background: 87
- Overlap: 11
- Raw one-sided hypergeometric P: 0.0008926000825948277
- BH FDR over 268 Group 5-hit KEGG pathways: 4.29e-02
- Fold enrichment: 3.02951362746574
- Derived Fisher odds ratio: 3.44048768710768
- Overlap genes: ID1, ID2, ID3, ID4, SMAD6, SMAD9, PITX2, RGMA, INHBE, NOG, ACVR2B

## Why the later complete-module ORA differs

The later analysis did not repeat the server test. It used MSigDB 2025.1 Hs C2:CP:KEGG_LEGACY, restricted its universe to the 5,000 clustered genes, retained only terms of size 15-500 after restriction, and corrected 40 retained terms within Group 5 x KEGG_LEGACY. Its effective universe was 691 genes, its Group 5 query was 143 genes, and the restricted TGF-beta set contained 24 genes. It found 10 overlaps (the same original genes except RGMA), P = 0.0143826371677454, and FDR = 0.352412081046915.

Thus the difference is explained by database membership, mapping, universe, term-size filtering, and BH scope. It is not evidence that either arithmetic result is erroneous.

## Reproducibility assessment

The original iDEP analysis is suitable for the planned Figure 4C pathway panel if it is identified exactly as an iDEP data113 KEGG ORA using filtered expressed genes as background. The figure should retain the narrower claim that Group 5 is a late-rising module enriched for the KEGG TGF-beta signaling gene set. It should not claim activation of the entire canonical pathway.

The CSV does not contain the internal KEGG release string from the database `source` table. The reproducible database identifier supported by the server code and archive is iDEP `data113`; the human database file was last modified 2025-09-11 and the exported hsa04350 size was 107. Preserve the data113 archive identifier, the 11-gene overlap, and this audit with the source data.
