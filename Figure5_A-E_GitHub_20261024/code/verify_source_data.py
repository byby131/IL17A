"""Verify saved Figure 5 source data; never fit a model or run enrichment.

python code/verify_source_data.py --root . --output validation/current_check.json
Optional --legacy reads the original folder names; --analysis-root checks the
archived analysis tables. Numerical comparisons only check saved transformations.
"""
from pathlib import Path
import argparse, hashlib, json, re
import numpy as np
import pandas as pd

PANEL_DIRS = {
    'A': ('A_volcano', 'A_volcano'),
    'B': ('B_DEG_heatmap', 'B_DEG_heatmap'),
    'C': ('E_functional_program_heatmap', 'C_functional_program_heatmap'),
    'D': ('F_Adipogenic_TF_heatmap', 'D_adipogenic_TF_heatmap'),
    'E': ('C_KEGG', 'E_KEGG'),
}
SAMPLES = ['NC_1','NC_2','NC_3','IL17_1','IL17_2','IL17_3']
C_GENES = 'CD36 GLUL ALDH1A1 H2AZ1 MAOA CRYZ CXCL1 SOD2 CXCL6 NAMPT SERPINE1 TNFAIP6 CXCL5 CXCL8 CXCL12 CCL2 MMP2 COL16A1 MMP1 PLOD2 COL5A1 SERPINH1 MMP3 COL1A1'.split()
C_FDR = '6.8e-81 2.4e-42 1.5e-13 1.2e-5 3.9e-5 0.00039 0* 0* 0* 2.3e-179 1.3e-159 1.2e-92 0* 0* 0* 2.3e-72 1.9e-87 5e-76 6.5e-75 3.1e-41 2.5e-38 1.3e-19 3.7e-17 6e-17'.split()
D_GENES = 'PPARG CEBPA CEBPB CEBPD KLF6 EBF2 FOXO1 PATZ1 PPARD ATF4 SREBF1 KLF15 ZNF423 PRDM16'.split()
D_FDR = '3.89e-9 0.139 7.43e-7 1.3e-5 2.14e-30 0.0123 0.0385 0.00826 0.00109 1.04e-7 0.127 0.509 0.432 0.105'.split()
E_UP = ['Coronavirus disease','Viral protein interaction with cytokine and cytokine receptor','IL-17 signaling pathway','Cytokine-cytokine receptor interaction','TNF signaling pathway','Mineral absorption','Rheumatoid arthritis','Biosynthesis of amino acids','HIF-1 signaling pathway','Protein digestion and absorption']
E_DOWN = ['Regulation of actin cytoskeleton','Cytoskeleton in muscle cells','IgSF CAM signaling','Focal adhesion','Proteoglycans in cancer','Hypertrophic cardiomyopathy','Adherens junction','MicroRNAs in cancer','Motor proteins','Bacterial invasion of epithelial cells']

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p): return pd.read_csv(p, dtype=str, keep_default_na=False, encoding='utf-8-sig')
def nums(df): return df.replace({'NA':np.nan,'':np.nan}).astype(float)
def fmt(v,digits):
    if v==0: return '0*'
    return re.sub(r'e([+-])0+(\d+)',r'e\1\2',f'{v:.{digits}g}')

def verify(root, legacy=False, analysis=None):
    checks=[]
    def ok(label,condition,detail=''):
        checks.append(dict(check=label,passed=bool(condition),detail=detail))
        if not condition: raise AssertionError(label+': '+str(detail))
    def equal(label,a,b):
        pd.testing.assert_frame_equal(a.reset_index(drop=True),b.reset_index(drop=True),check_dtype=False)
        ok(label,True,f'{len(a)} rows; exact stored strings')
    def close(label,a,b):
        np.testing.assert_allclose(np.asarray(a,dtype=float),np.asarray(b,dtype=float),rtol=1e-11,atol=1e-12,equal_nan=True)
        ok(label,True,'rtol=1e-11; atol=1e-12 (display transformations only)')
    data=root if legacy else root/'data'
    folders={key:data/pair[0 if legacy else 1] for key,pair in PANEL_DIRS.items()}
    if legacy:
        mf=read(root/'manifest.csv')
        for r in mf.to_dict('records'): ok('Original SHA256: '+r['file'],sha(root/r['file'])==r['sha256'])
        ok('Original manifest covers every data CSV',set(mf.file)=={p.relative_to(root).as_posix() for p in root.rglob('*.csv') if p.name!='manifest.csv'})
    else:
        mf=read(root/'manifest.csv')
        for r in mf.to_dict('records'): ok('Release SHA256: '+r['path'],sha(root/r['path'])==r['sha256'])
        ok('Manifest covers release payload',set(mf.path)=={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and p.name!='manifest.csv' and p.relative_to(root).as_posix()!='validation/current_check.json' and '__pycache__' not in p.parts})
    meta=read(data/'shared/sample_metadata.csv')
    ok('Six samples in fixed plot order',meta.sample_id.tolist()==SAMPLES)
    ok('Day 7 independent comparison',meta.harvest_day.eq('7').all() and meta.replicate_type.eq('independent_biological').all() and meta.condition.tolist()==['vehicle']*3+['IL17A']*3)
    raw=read(data/'shared/original_counts_input.csv')
    raw=raw.rename(columns={raw.columns[0]:'gene_id',**dict(zip(meta.count_column,meta.sample_id))}).set_index('gene_id')
    ok('Complete input shape and unique IDs',raw.shape==(77775,6) and raw.index.is_unique)
    raw_n=nums(raw[SAMPLES]);ok('Counts nonnegative integers',np.isfinite(raw_n).all().all() and (raw_n>=0).all().all() and (raw_n%1==0).all().all())
    sf=read(data/'shared/DESeq2_size_factors.csv').set_index('sample_id').loc[SAMPLES,'size_factor'].astype(float)
    de=read(folders['A']/'all_tested_gene_statistics.csv').set_index('gene_id')
    ok('17887 tested genes, unique IDs',len(de)==17887 and de.index.is_unique)
    keep=raw_n.index[(raw_n>=10).sum(axis=1)>=3]
    ok('Saved tested genes match recorded count filter',set(keep)==set(de.index))
    padj=pd.to_numeric(de.padj,errors='coerce')
    ok('DEG direction totals',(de.direction.value_counts().to_dict()=={'not_significant':16059,'down':926,'up':897,'not_estimable':5}))
    ok('Five missing FDR, seven underflow zeros',padj.isna().sum()==5 and padj.eq(0).sum()==7)
    ok('NA FDR caused by Cook exclusions',de.loc[padj.isna(),'na_reason'].eq('cooks_outlier_pvalue_excluded').all())
    ok('DEG threshold is FDR only',set(de.index[padj<.05])==set(de.index[de.direction.isin(['up','down'])]))
    if analysis:
        equal('All statistics equal archived analysis',de.reset_index(),read(analysis/'tables/deseq2_all_tested.csv'))
        ok('Original input matches archived SHA',sha(data/'shared/original_counts_input.csv')=='4af7b20447ff0b8714a211ccd193c4d66a2bf1ff5b1cb30846c6886fbfae7913')
    points=read(folders['A']/'plotted_points.csv')
    ok('A: finite point count and membership',len(points)==17882 and set(points.gene_id)==set(de.index[padj.notna()]))
    equal('A: all plotted statistical fields',points[list(de.reset_index().columns)],de.loc[points.gene_id].reset_index())
    close('A: x display clipping',nums(points[['plot_x_log2FC']]).values[:,0],np.clip(points.log2FC.astype(float),-8,8))
    y=-np.log10(np.maximum(points.padj.astype(float),np.finfo(float).tiny))
    close('A: minus log10 FDR',points.minus_log10_FDR,y)
    close('A: y display clipping',points.plot_y_minus_log10_FDR,np.minimum(y,60))
    ok('A: point color matches direction',points.point_color.tolist()==points.direction.map({'up':'#FF0000','down':'#0000FF','not_significant':'#777777'}).tolist())
    labels=read(folders['A']/'top_gene_labels.csv')
    ok('A: image gene labels',set(labels.gene_id)==set('A2M SLC40A1 CIDEC MGP CD36 CXCL6 CXCL1 CXCL8 CXCL5 CXCL12'.split()))
    ok('A: top-label flags',set(points.loc[points.top_label.str.lower().eq('true'),'gene_id'])==set(labels.gene_id))
    for panel,folder in folders.items():
        for filename in ['raw_counts.csv','normalized_counts.csv','VST_expression.csv']:
            m=read(folder/filename).set_index('gene_id')
            ok(panel+': '+filename+' sample order',m.columns.tolist()==SAMPLES and m.index.is_unique)
            if filename=='raw_counts.csv':equal(panel+': raw count provenance',m,raw.loc[m.index,SAMPLES])
            elif filename=='normalized_counts.csv':close(panel+': normalization by saved size factors',nums(m),nums(raw.loc[m.index,SAMPLES]).div(sf,axis=1))
            if analysis and filename in ['normalized_counts.csv','VST_expression.csv']:
                target='normalized_counts.csv' if filename=='normalized_counts.csv' else 'vst_matrix.csv'
                ref=read(analysis/'tables'/target).set_index('gene_id').loc[m.index,SAMPLES]
                equal(panel+': '+filename+' archived values',m,ref)
        statfile='hit_gene_differential_statistics.csv' if panel=='E' else 'differential_statistics.csv'
        if panel!='A':
            stats=read(folder/statfile)
            equal(panel+': statistics match whole-genome table',stats,de.loc[stats.gene_id].reset_index())
        if panel not in ['B','C','D']: continue
        rows=read(folder/'row_annotations.csv');ids=rows.gene_id.tolist()
        ok(panel+': consecutive row order',rows.row_order.astype(int).tolist()==list(range(1,len(rows)+1)))
        v=read(folder/'VST_expression.csv').set_index('gene_id')
        z=read(folder/'heatmap_VST_rowZ.csv').set_index('gene_id')
        clip=read(folder/'heatmap_color_values_clipped.csv').set_index('gene_id')
        ok(panel+': all matrices share row order',v.index.tolist()==z.index.tolist()==clip.index.tolist()==ids)
        vn=nums(v); zn=nums(z)
        close(panel+': row Z from saved VST',zn,vn.sub(vn.mean(axis=1),axis=0).div(vn.std(axis=1,ddof=1),axis=0))
        limit=2.5 if panel=='B' else 2
        close(panel+': heatmap color clipping',nums(clip),zn.clip(-limit,limit))
        if panel=='B':
            ok('B: all 1823 DEGs, up then down',len(ids)==1823 and set(ids)==set(de.index[padj<.05]) and rows.direction.tolist()==['up']*897+['down']*926)
            if analysis:equal('B: exact archived row order',rows,read(analysis/'tables/panel_B_row_order.csv'))
        else:
            expected=C_GENES if panel=='C' else D_GENES
            ok(panel+': image gene order',ids==expected)
            fd=de.loc[ids,'padj'].astype(float)
            ok(panel+': all image FDR labels',[fmt(x,2 if panel=='C' else 3) for x in fd]==(C_FDR if panel=='C' else D_FDR))
            if panel=='D':ok('D: nine starred genes',list(fd.index[fd<.05])=='PPARG CEBPB CEBPD KLF6 EBF2 FOXO1 PATZ1 PPARD ATF4'.split())
            if panel=='C':
                ok('C: six underflow footnotes',list(fd.index[fd==0])=='CXCL1 SOD2 CXCL6 CXCL5 CXCL8 CXCL12'.split())
                if analysis:equal('C: all row and program annotations',rows.drop(columns='row_order'),read(analysis/'figure5_extension/tables/functional_heatmap_rows.csv'))
    pk=read(folders['C']/'program_color_key.csv')
    ok('C: program codes/colors',pk.program_code.tolist()==['P1','P2','P3','P4','P5'] and pk.color.tolist()==['#0072B2','#E69F00','#009E73','#CC79A7','#56B4E9'])
    k=read(folders['E']/'all_20_displayed_pathways.csv')
    ok('E: image pathway order',k[k.direction=='up'].pathway_name.tolist()==E_UP and k[k.direction=='down'].pathway_name.tolist()==E_DOWN)
    close('E: bar lengths are minus log10 FDR',k.bar_length_minus_log10_FDR,-np.log10(k.FDR.astype(float)))
    ora=read(folders['E']/'background_and_complete_results/KEGG_ORA_all_results.csv')
    ok('E: complete ORA families',ora.direction.value_counts().to_dict()=={'all':329,'up':329,'down':329})
    ok('E: joint direction correction family',ora.loc[ora.direction!='all','correction_family'].eq('up_and_down_joint_all_eligible_KEGG_pathways').all())
    ok('E: background size',len(read(folders['E']/'background_and_complete_results/KEGG_universe.csv'))==6002 and k.universe_size.eq('6002').all())
    for direction in ['up','down']:
        q=k[k.direction==direction]
        equal('E: '+direction+' split table',read(folders['E']/(direction+'_10_pathways.csv')),q)
        t=ora[ora.direction==direction].copy();t['f']=t.FDR.astype(float);t['fe']=t.fold_enrichment.astype(float)
        top=t[t.f<.05].sort_values(['f','fe','pathway_id'],ascending=[True,False,True]).head(10)
        equal('E: '+direction+' matches full-result top 10',q[list(ora.columns)],top[list(ora.columns)])
    members=read(folders['E']/'pathway_gene_membership.csv')
    ok('E: membership has no duplicate links',not members.duplicated(['direction','pathway_id','gene_id']).any())
    for r in k.to_dict('records'):
        genes=set(members.loc[(members.direction==r['direction'])&(members.pathway_id==r['pathway_id']),'gene_id'])
        ok('E: membership '+r['direction']+' '+r['pathway_id'],genes==set(r['hit_genes'].split(';')) and len(genes)==int(r['Count']))
    return dict(status='PASS',check_count=len(checks),checks=checks,
        scope='Saved CSVs, archived result lineage (when supplied), and manually transcribed visible labels/order in the user-provided Figure 5. No model fit or enrichment test rerun; no pixel-level identity claim.',
        totals=dict(input_rows=77775,tested_rows=17887,volcano_points=17882,up_DEG=897,down_DEG=926,functional_genes=24,TF_genes=14,KEGG_displayed=20,KEGG_background=6002))

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=Path(__file__).resolve().parents[1]);ap.add_argument('--legacy',action='store_true');ap.add_argument('--analysis-root',type=Path);ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    report=verify(args.root,args.legacy,args.analysis_root)
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True);args.output.write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='checks'},ensure_ascii=False,indent=2))
