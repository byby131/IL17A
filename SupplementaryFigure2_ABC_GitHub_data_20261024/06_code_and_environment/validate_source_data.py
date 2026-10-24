"""Validate current S2A-C files and statistics; no fitting or GSEA rerun.

Run: python 06_code_and_environment/validate_source_data.py
Requires numpy and pandas. The validation is read-only.
"""
from pathlib import Path
import hashlib,json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
A='02_panel_source_data/A_pathway_enrichment/'
B='02_panel_source_data/B_pathway_heterogeneity/'
C='02_panel_source_data/C_secreted_genes/'
SA='03_complete_statistics/A_GSEA/'
SC='03_complete_statistics/C_secretome/'

def read(name): return pd.read_csv(ROOT/name)
def close(a,b,atol=1e-10,rtol=1e-10):np.testing.assert_allclose(a,b,atol=atol,rtol=rtol,equal_nan=True)
def sha(p):
    with open(p,'rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def bh(p):
    p=np.asarray(p);order=np.argsort(p,kind='stable');q=p[order]*len(p)/np.arange(1,len(p)+1)
    q=np.minimum(1,np.minimum.accumulate(q[::-1])[::-1]);out=np.empty(len(p));out[order]=q;return out

manifest=read('MANIFEST_SHA256.csv')
for r in manifest.itertuples(index=False):
    p=ROOT/r.path
    assert p.is_file() and p.stat().st_size==r.bytes,r.path
    assert sha(p)==r.sha256,r.path
allfiles={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts}
assert allfiles==set(manifest.path)|{'MANIFEST_SHA256.csv'}
assert all((ROOT/p).is_file() for p in read('PANEL_INDEX.csv').primary_source_file)
assert all((ROOT/p).is_file() for p in read('PROVENANCE.csv').file)

cfg=json.loads((ROOT/'07_provenance_and_QC/current_panel_config.json').read_text())
samples=cfg['sample_order'];cons=cfg['A_contrast_order']
raw=read('01_input_counts/original_transcript_Salmon_counts.csv')
assert raw.shape==(232117,13) and raw.Name.is_unique
assert np.isfinite(raw[samples].values).all() and (raw[samples]>=0).all().all()
meta=read('01_input_counts/sample_metadata.csv').set_index('sample_id').loc[samples]
assert len(meta)==12 and meta.donor_code.nunique()==1 and meta.replicate.nunique()==3
assert meta.groupby('replicate').time.apply(lambda x:set(x)=={0,3,12,24}).all()
gene=read('04_frozen_inputs_and_models/sender_gene_counts.csv.gz').set_index('gene_id')
assert len(gene)==60660
close(raw[samples].sum(),gene[samples].sum(),atol=1e-5)

# A: recompute full-library BH from stored nominal P values.
full=read(SA+'A_GSEA_all_Hallmark_Reactome_5_contrasts.csv')
assert len(full)==5120 and not full.duplicated(['pathway','contrast']).any()
for (library,cn),d in full.groupby(['library','contrast']):
    assert len(d)=={'Hallmark':49,'Reactome':975}[library]
    close(d.padj_library_contrast,bh(d.pval),atol=1e-300,rtol=1e-9)
display=read(A+'A_display_GSEA_long_88_cells.csv')
assert len(display)==88
fullidx=full.set_index(['pathway','contrast'])
for r in display.itertuples(index=False):
    orig=fullidx.loc[(r.pathway,r.contrast)]
    close([r.NES,r.padj_library_contrast,r.pval,r.size],[orig.NES,orig.padj_library_contrast,orig.pval,orig['size']])
for filename,field in [('A_NES_matrix.csv','NES'),('A_padj_library_contrast_matrix.csv','padj_library_contrast'),('A_pval_matrix.csv','pval'),('A_size_matrix.csv','size')]:
    table=read(A+filename).set_index('pathway')
    assert table.index.tolist()==cfg['A_pathway_order'] and table.columns.tolist()==cons
    pivot=display.pivot(index='pathway',columns='contrast',values=field).loc[table.index,cons]
    close(table.values,pivot.values,atol=1e-300 if field in ['pval','padj_library_contrast'] else 1e-10)
labels=read('07_provenance_and_QC/A_all88_screenshot_label_checks.csv')
for r in labels.itertuples(index=False):
    q=fullidx.loc[(r.pathway,r.contrast),'padj_library_contrast']
    text=(f'{q:.0e}' if q<.001 else f'{q:.3f}')+('*' if q<.05 else '')
    assert text==r.figure_text==r.source_text

# B: every plotted effect comes from the original genome-wide DE table.
effects=read(B+'B_full_pathway_gene_effects_and_Group5_displayed4sets.csv')
assert len(effects)==789 and set(effects.pathway)==set(cfg['B_pathway_order'])
for cn,d in effects.groupby('contrast'):
    de=read('03_complete_statistics/DESeq2/DESeq2_'+cn+'.csv').set_index('gene_id')
    close(d.log2FoldChange,de.loc[d.gene_id,'log2FoldChange'])
    close(d.gene_BH_FDR,de.loc[d.gene_id,'padj'],atol=1e-300)
    close(d.pvalue,de.loc[d.gene_id,'pvalue'],atol=1e-300)
server=read('04_frozen_inputs_and_models/frozen_server_5000_gene_membership_and_display_values.csv')
g5=set(server.loc[server.cluster==5,'ensembl_ID'].str.replace(r'\.[0-9]+$','',regex=True))
assert (effects.in_Group5==effects.unversioned.isin(g5)).all()
composition=read(B+'B_whole_pathway_direction_composition_displayed4sets_all3times.csv')
for (path,cn),sub in effects.groupby(['pathway','contrast']):
    states=np.where(sub.gene_BH_FDR.isna()|sub.log2FoldChange.isna(),'Unavailable',np.where((sub.gene_BH_FDR<.05)&(sub.log2FoldChange>0),'Up',np.where((sub.gene_BH_FDR<.05)&(sub.log2FoldChange<0),'Down','Not significant')))
    assert (states==sub.direction).all()
    comp=composition[(composition.pathway==path)&(composition.contrast==cn)]
    for r in comp.itertuples(index=False):
        assert r.count==int((states==r.status).sum())
        assert r.total_detectable==len(sub)
        close(r.proportion,r.count/r.total_detectable)
    assert comp['count'].sum()==len(sub)
points=read(B+'B_boxplot_all526_gene_points.csv')
assert len(points)==526
expected=effects[effects.contrast.isin(['12_vs_0','24_vs_0'])]
pd.testing.assert_frame_equal(points.reset_index(drop=True),expected.reset_index(drop=True),check_exact=False,rtol=1e-12,atol=1e-14)
boxes=read(B+'B_boxplot_summary_16_boxes.csv');assert len(boxes)==16
for r in boxes.itertuples(index=False):
    vals=points.loc[(points.pathway==r.pathway)&(points.contrast==r.contrast)&(points.in_Group5==r.in_Group5),'log2FoldChange'].values
    q1,median,q3=np.percentile(vals,[25,50,75]);iqr=q3-q1
    assert len(vals)==r.n_genes
    close([r.q1,r.median,r.q3,r.whisker_low,r.whisker_high],[q1,median,q3,min(q1,vals[vals>=q1-1.5*iqr].min()),max(q3,vals[vals<=q3+1.5*iqr].max())])
le=read(B+'B_leading_edge_summary_displayed4sets.csv')
legenes=read(B+'B_leading_edge_genes_displayed4sets.csv')
assert len(le)==12 and len(legenes)==207
for r in le.itertuples(index=False):
    subset=legenes[(legenes.pathway==r.pathway)&(legenes.contrast==r.contrast)]
    assert len(subset)==r.n_leading_edge and subset.in_Group5.sum()==r.n_Group5
    fr=fullidx.loc[(r.pathway,r.contrast)]
    assert set(subset.ENTREZID.astype(int).astype(str))==set(str(fr.leadingEdge_Entrez).split(';'))
    close([r.NES,r.pathway_BH_FDR],[fr.NES,fr.padj_library_contrast])
    assert r.pathway_BH_FDR>=.05
bcheck=read('07_provenance_and_QC/B_screenshot_count_and_LE_checks.csv')
for r in bcheck.itertuples(index=False):
    cc=composition[(composition.pathway==r.pathway)&(composition.contrast==r.contrast)].set_index('status')
    assert '/'.join(str(int(cc.loc[s,'count'])) for s in ['Up','Down','Not significant'])==r.figure_counts
    ll=le[(le.pathway==r.pathway)&(le.contrast==r.contrast)].iloc[0]
    assert ('+' if ll.NES>0 else '-')+f'{int(ll.n_Group5)}/{int(ll.n_leading_edge)}'==r.figure_le
for path,expected_counts in zip(cfg['B_pathway_order'],[(4,51),(4,24),(3,47),(7,141)]):
    e=effects[(effects.pathway==path)&(effects.contrast=='24_vs_0')]
    assert (int(e.in_Group5.sum()),len(e))==expected_counts

# C: raw counts -> saved normalization, saved VST -> Z -> scores/SD.
hm=read(C+'C_all500genes_12samples_gene_wise_Z.csv')
assert len(hm)==500 and hm.gene_id.is_unique and hm.plotted_row_1based.tolist()==list(range(1,501))
assert hm.groupby('secreted_module').size().to_dict()=={'S1':337,'S2':163}
assert hm.iloc[:337].secreted_module.eq('S1').all() and hm.iloc[337:].secreted_module.eq('S2').all()
raw500=read(C+'C_all500genes_12samples_raw_aggregated_Salmon_estimated_counts.csv').set_index('gene_id').loc[hm.gene_id]
integers=read(C+'C_all500genes_12samples_model_integer_counts.csv').set_index('gene_id').loc[hm.gene_id]
norm=read(C+'C_all500genes_12samples_DESeq2_normalized_counts.csv').set_index('gene_id').loc[hm.gene_id]
vst=read(C+'C_all500genes_12samples_VST.csv').set_index('gene_id').loc[hm.gene_id]
close(raw500[samples],gene.loc[hm.gene_id,samples])
close(integers[samples],np.rint(raw500[samples]),atol=0,rtol=0)
close(norm[samples],integers[samples].to_numpy()/meta.sizeFactor.to_numpy())
allvst=read('04_frozen_inputs_and_models/VST_all_14425_genes.csv').set_index('gene_id')
close(vst[samples],allvst.loc[hm.gene_id,samples])
v=vst[samples].to_numpy();z=(v-v.mean(axis=1,keepdims=True))/v.std(axis=1,ddof=1,keepdims=True)
close(hm[samples],z)
elig=read(SC+'C_all_tested_genes_secreted_annotation_and_LRT.csv')
assert len(elig)==14425 and elig.secreted_annotation.sum()==820
eligible=elig.secreted_annotation & (elig.padj<.05) & elig.gene_id.isin(allvst.index)
assert set(elig.loc[eligible,'gene_id'])==set(hm.gene_id)
scores=read(C+'C_24_module_sample_points.csv');curves=read(C+'C_curve_mean_SD_and_shaded_bounds.csv')
assert len(scores)==24 and len(curves)==8
for r in scores.itertuples(index=False):
    sample=f'S50_{r.time_h}_{r.replicate}'
    members=hm[hm.secreted_module==r.module]
    assert len(members)==r.n_genes
    close(r.module_score,members[sample].mean())
for r in curves.itertuples(index=False):
    vals=scores[(scores.module==r.module)&(scores.time_h==r.time_h)].sort_values('replicate').module_score
    m=vals.mean();sd=vals.std(ddof=1)
    assert len(vals)==r.n_independent_experiments==3
    close([r.R1,r.R2,r.R3],vals)
    close([r.mean_module_score,r.SD,r.shaded_lower,r.shaded_upper],[m,sd,m-sd,m+sd])
for symbol,module in [('CXCL8','S1'),('IL6','S1'),('CSF3','S2'),('CXCL1','S2')]:
    assert hm.loc[hm.gene_name==symbol,'secreted_module'].tolist()==[module]
fq=read('08_FASTQ_inventory/FASTQ_24_files_manifest.csv')
assert len(fq)==24 and set(fq['sample'])==set(samples)
assert fq.groupby('sample').read_mate.apply(lambda x:set(x)=={'R1','R2'}).all()
assert not fq.reads_included_in_this_GitHub_package.any()

print(json.dumps({'status':'PASS','manifest_files_verified':len(manifest),'raw_transcripts':len(raw),'samples':12,'A_cells':len(display),'GSEA_full_library_tests':len(full),'B_composition_triplets':12,'B_boxplot_gene_points':len(points),'B_boxes':len(boxes),'B_leading_edge_gene_records':len(legenes),'C_genes':len(hm),'C_module_points':len(scores),'C_mean_SD_summaries':len(curves),'FASTQ_inventory_only_files':len(fq),'no_refitting_GSEA_or_clustering':True},indent=2))
