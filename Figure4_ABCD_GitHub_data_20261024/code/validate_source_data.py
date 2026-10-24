"""Validate the packaged source data without accessing the original drives.

Usage: python code/validate_source_data.py
Dependencies: pandas, numpy. Does not refit DESeq2 or rerun iDEP.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

def read(name):
    return pd.read_csv(ROOT / name)

def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def close(actual, expected, atol=1e-9, rtol=1e-10):
    np.testing.assert_allclose(actual, expected, atol=atol, rtol=rtol)

manifest = read('MANIFEST_SHA256.csv')
for row in manifest.itertuples(index=False):
    p = ROOT / row.path
    assert p.is_file(), row.path
    assert p.stat().st_size == row.bytes, row.path
    assert sha(p) == row.sha256, row.path
actual_files = {p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
assert actual_files == set(manifest.path) | {'MANIFEST_SHA256.csv'}, 'Unexpected or missing files'

samples = [f'S50_{t}_{r}' for t in (0, 3, 12, 24) for r in (1, 2, 3)]
counts = read('01_input_counts/ADSC_IL17A_timecourse_RNAseq_EXPRESSION_MATRIX_COUNTS.CSV')
assert counts.shape == (232117, 13)
assert counts.Name.is_unique and set(counts.columns) == {'Name', *samples}
assert np.isfinite(counts[samples].values).all() and (counts[samples] >= 0).all().all()
meta = read('01_input_counts/sample_metadata.csv')
assert meta.sample_id.is_unique and set(meta.sample_id) == set(samples)
assert meta.study_donor_label.nunique() == 1 and meta.experiment_id.nunique() == 3
assert meta.groupby('experiment_id').time_h.apply(lambda x: set(x) == {0, 3, 12, 24}).all()
gene_counts = read('04_processed_support/sender_gene_counts.csv.gz')
assert len(gene_counts) == 60660
close(counts[samples].sum().values, gene_counts[samples].sum().values, atol=1e-5)

heat = read('03_panel_source_data/A_heatmap/heatmap_5000_genes_display_order.csv')
server = read('02_frozen_server_exports/heatmap_data.csv')
assert heat.groupby('cluster').size().tolist() == [333,731,811,898,1043,618,566]
fields = ['ensembl_ID','symbol','cluster'] + samples
pd.testing.assert_frame_equal(heat[fields].sort_values('ensembl_ID').reset_index(drop=True), server[fields].sort_values('ensembl_ID').reset_index(drop=True))
late = heat[[f'S50_{t}_{r}' for t in (12,24) for r in (1,2,3)]].mean(axis=1)
early = heat[[f'S50_{t}_{r}' for t in (0,3) for r in (1,2,3)]].mean(axis=1)
close(heat.display_late_delta, late-early)
assert heat.cluster.is_monotonic_increasing
assert heat.groupby('cluster').display_late_delta.apply(lambda x:x.is_monotonic_decreasing).all()

z = read('04_processed_support/VST_gene_Z_matrix_all_matched_modules.csv')
close(z[samples].mean(axis=1), 0)
close(z[samples].std(axis=1,ddof=1), 1)
scores = read('03_panel_source_data/B_trajectories/Group1_Group5_sample_module_scores.csv')
means = read('03_panel_source_data/B_trajectories/Group1_Group5_mean_and_CI95.csv')
assert len(scores) == 24 and len(means) == 8
for row in scores.itertuples(index=False):
    matrix = z.loc[z.cluster == row.cluster]
    assert row.n_genes == len(matrix) == {1:147,5:800}[row.cluster]
    close(row.module_score, matrix[row.sample_id].mean())
for row in means.itertuples(index=False):
    v = scores.loc[(scores.cluster == row.cluster)&(scores.time_h == row.time_h),'module_score']
    assert len(v) == row.n_samples == 3
    m = v.mean(); se = v.std(ddof=1) / np.sqrt(3); halfwidth = 4.302652729696142 * se
    close([row.mean_module_score,row.sem,row.ci95_low,row.ci95_high], [m,se,m-halfwidth,m+halfwidth])

display = read('03_panel_source_data/C_KEGG/KEGG_displayed_9_pathways.csv')
enrich = read('02_frozen_server_exports/enrichment.csv')
audit = read('03_panel_source_data/C_KEGG/Group5_KEGG_audited_all8_exported_terms.csv')
assert display.groupby('group').size().to_dict() == {1:4,5:5}
for row in display.itertuples(index=False):
    if row.group == 1:
        original = enrich.loc[(enrich.group==1)&(enrich.Pathway==row.pathway)].iloc[0]
        close([row.fold_enrichment,row.FDR_exported,row.overlap_count],[original['Fold enriched'],original.FDR,original.nGenes])
    else:
        original = audit.loc[audit.pathway==row.pathway].iloc[0]
        close([row.fold_enrichment,row.FDR_exported,row.overlap_count,row.raw_P_reconstructed,row.BH_FDR_reconstructed],[original.fold_enrichment_reconstructed,original.server_BH_FDR_exported,original.overlap_count,original.raw_P_reconstructed,original.BH_FDR_reconstructed_from_268])
        close(row.fold_enrichment, (row.overlap_count/row.query_genes)/(row.effective_background_pathway_size/row.background_genes))

order = ['ID1','ID2','ID3','ACVR2B','ID4','SMAD9','PITX2','NOG','SMAD6','RGMA','INHBE']
d = read('03_panel_source_data/D_TGFbeta_genes/11_gene_timepoint_means_display_order.csv').set_index('symbol')
ds = read('03_panel_source_data/D_TGFbeta_genes/11_gene_server_expression_12_samples.csv').set_index('symbol').loc[order]
dt = read('03_panel_source_data/D_TGFbeta_genes/heatmap_4_times_by_11_genes.csv').set_index('time_h')
assert d.index.tolist() == dt.columns.tolist() == order
assert dt.index.tolist() == [0,3,12,24]
for t in (0,3,12,24):
    close(d[str(t)],ds[[f'S50_{t}_{r}' for r in (1,2,3)]].mean(axis=1))
    close(dt.loc[t],d[str(t)])
frozen_subset = server.set_index('symbol').loc[order]
close(ds[samples], frozen_subset[samples])
assert (frozen_subset.cluster==5).all()

print(json.dumps({'status':'PASS','manifest_files_verified':len(manifest),'counts_transcripts':len(counts),'samples':len(samples),'aggregated_genes':len(gene_counts),'A_genes':len(heat),'B_sample_values':len(scores),'B_time_means':len(means),'C_pathways':len(display),'D_genes':len(d),'checks':['SHA256','counts and sample metadata','transcript-to-gene count conservation','A original values and ordering','B scores and 95% t intervals','C original enrichment values','D original sample values and time means']},indent=2))
