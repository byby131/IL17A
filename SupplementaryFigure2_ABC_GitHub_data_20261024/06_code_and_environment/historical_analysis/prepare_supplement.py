from pathlib import Path
import sys, json, hashlib, shutil, itertools, platform
import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import linkage, fcluster, leaves_list
from sklearn.metrics import silhouette_score, adjusted_rand_score

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(BASE/'dependencies'))
import gseapy
SRC=BASE/'source_data'; SRC.mkdir(exist_ok=True)
SAMPLES=[f'S50_{h}_{r}' for h in (0,3,12,24) for r in (1,2,3)]
TIMES=[0,3,12,24]
CONTRASTS=['3_vs_0','12_vs_0','24_vs_0','24_vs_3']
def save(x,name): x.to_csv(SRC/name,index=False)
def hashfile(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def bh(x):
    x=np.asarray(x,float); out=np.full(x.shape,np.nan); ok=np.isfinite(x); a=x[ok]; ix=np.argsort(a); q=np.minimum.accumulate((a[ix]*len(a)/np.arange(1,len(a)+1))[::-1])[::-1]; b=np.empty(len(a)); b[ix]=np.minimum(q,1);out[ok]=b;return out
def gmt(p):
    return {a[0]:list(dict.fromkeys(a[2:])) for a in (s.rstrip('\n').split('\t') for s in p.open())}

input_paths=[ROOT/'analysis/adsc_timecourse_counts/DESeq2_LRT_any_time_effect.csv',ROOT/'analysis/figure4_final_frozen/source_data/Figure4B_Hallmark_GSEA_all_5_contrasts.csv',ROOT/'analysis/figure4_scheme2/source_data/canonical_GSEA_KEGG_MEDICUS_Reactome_all_5_contrasts.csv',ROOT/'HFF-1/round3_20260914/cache/sender_dds.rds',ROOT/'HFF-1/round3_20260914/cache/references/uniprot_human_reviewed_secreted.tsv']
for cn in CONTRASTS:input_paths.append(ROOT/f'analysis/figure4_abc/source_data/DESeq2_{cn}.csv')
refs=ROOT/'adsc_il17_rnaseq/references'
for nm in ['h.all.v2025.1.Hs.entrez.gmt','c2.cp.reactome.v2025.1.Hs.entrez.gmt']:input_paths.append(refs/nm)
manifest=[{'file':str(p.relative_to(ROOT)),'sha256':hashfile(p)} for p in input_paths]
save(pd.DataFrame(manifest),'reused_input_hashes.csv')
for p in [input_paths[0]]+input_paths[5:9]:shutil.copy2(p,SRC/p.name)
shutil.copy2(input_paths[4],SRC/'UniProt_reviewed_human_secreted_frozen.tsv')
for nm in ['h.all.v2025.1.Hs.entrez.gmt','c2.cp.reactome.v2025.1.Hs.entrez.gmt']:shutil.copy2(refs/nm,SRC/nm)
sets={**gmt(refs/'h.all.v2025.1.Hs.entrez.gmt'),**gmt(refs/'c2.cp.reactome.v2025.1.Hs.entrez.gmt')}
library=lambda p:'Hallmark' if p.startswith('HALLMARK_') else 'Reactome'
h=pd.read_csv(input_paths[1]);h['library']='Hallmark';h['padj_library_contrast']=h.padj
r=pd.read_csv(input_paths[2]);r=r[r.library.eq('Reactome')].copy()
g=pd.concat([h,r],ignore_index=True,sort=False)
save(g,'B_GSEA_all_Hallmark_Reactome_5_contrasts.csv')

# Display selector is based only on pathway names, never NES/FDR.
hall_focus=['HALLMARK_INFLAMMATORY_RESPONSE','HALLMARK_TNFA_SIGNALING_VIA_NFKB','HALLMARK_IL6_JAK_STAT3_SIGNALING','HALLMARK_TGF_BETA_SIGNALING','HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION']
react_exact=['REACTOME_SIGNALING_BY_BMP','REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX','REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS','REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING','REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING','REACTOME_INTERLEUKIN_6_SIGNALING','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX','REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS','REACTOME_COLLAGEN_FORMATION','REACTOME_COLLAGEN_DEGRADATION','REACTOME_NFKB_ACTIVATION_THROUGH_FADD_RIP_1_PATHWAY_MEDIATED_BY_CASPASE_8_AND_10']
react_inflammation=sorted(p for p in sets if p.startswith('REACTOME_') and any(k in p for k in ['INFLAMMA','NFKB_ACTIVATION','NFKB_SIGNAL','NF_KB_ACTIVATION']))
focus=list(dict.fromkeys(hall_focus+react_exact+react_inflammation))
display=[p for p in focus if p in set(g.pathway)]
save(g[g.pathway.isin(display)&g.contrast.isin(CONTRASTS)].copy(),'B_display_GSEA.csv')

v=pd.read_csv(SRC/'VST_all_14425_genes.csv').set_index('gene_id')
norm=pd.read_csv(SRC/'normalized_counts_all_14425_genes.csv').set_index('gene_id')
mapping=pd.read_csv(SRC/'Ensembl_Entrez_mapping_frozen_annotation.csv',dtype=str).dropna(subset=['ENTREZID'])
ambiguous=mapping.ENSEMBL.value_counts(); mapping=mapping[~mapping.ENSEMBL.isin(ambiguous[ambiguous>1].index)].drop_duplicates('ENSEMBL')
gene_map=pd.DataFrame({'gene_id':v.index,'ENSEMBL':v.index.str.replace(r'\.\d+$','',regex=True),'mean_normalized_count':norm.loc[v.index,SAMPLES].mean(axis=1).values})
gene_map=gene_map.merge(mapping,on='ENSEMBL').sort_values(['mean_normalized_count','gene_id'],ascending=[False,True]).drop_duplicates('ENTREZID').sort_values('ENTREZID')
save(gene_map,'C_ssGSEA_gene_mapping.csv')
X=v.loc[gene_map.gene_id,SAMPLES].copy();X.index=gene_map.ENTREZID
N=norm.loc[gene_map.gene_id,SAMPLES].copy();N.index=gene_map.ENTREZID
universe=set(X.index)
eligible={p:sorted(set(ids)&universe) for p,ids in sets.items() if 15<=len(set(ids)&universe)<=500}
audit=pd.DataFrame([{'library':library(p),'pathway':p,'database_size':len(ids),'ssGSEA_mapped_size':len(set(ids)&universe),'ssGSEA_eligible':p in eligible,'fgsea_tested':p in set(g.pathway),'focus':p in focus,'mapped_members':';'.join(sorted(set(ids)&universe))} for p,ids in sets.items()])
save(audit,'gene_set_size_and_membership_audit.csv')
save(pd.DataFrame([{'pathway':p,'present_in_database':p in sets,'present_in_frozen_GSEA':p in set(g.pathway),'displayed':p in display} for p in focus+['HALLMARK_BMP_SIGNALING','HALLMARK_EXTRACELLULAR_MATRIX_ORGANIZATION']]),'requested_program_availability.csv')
save(pd.DataFrame({'ENTREZID':sorted(universe)}),'C_ssGSEA_background_genes.csv')

print('Score complete libraries:',len(eligible),'sets,',len(X),'mapped genes',flush=True)
def score(data,weight,tag):
    # Explicit first ID column: GSEApy treats numeric-string indices as row numbers.
    table=data.copy(); table.index.name='Entrez'; table=table.reset_index()
    assert table.Entrez.is_unique and len(table.columns)==13
    res=gseapy.ssgsea(data=table,gene_sets=eligible,outdir=None,sample_norm_method='rank',correl_norm_type='rank',min_size=15,max_size=500,permutation_num=0,weight=weight,threads=2,no_plot=True,seed=20261008,verbose=False)
    q=res.res2d.rename(columns={'Name':'sample_id','Term':'pathway','ES':'ES','NES':'ssGSEA_normalized_score'}).copy()
    for c in ['ES','ssGSEA_normalized_score']:q[c]=q[c].astype(float)
    q['method']=tag;q['library']=q.pathway.map(library)
    q['time_h']=q.sample_id.str.split('_').str[1].astype(int);q['replicate']=q.sample_id.str.split('_').str[2].astype(int)
    save(q,f'C_ssGSEA_{tag}_all_samples.csv');return q
primary=score(X,0.25,'VST_rank_weight025')
alt_log=score(np.log2(N+1),0.25,'logcounts_rank_weight025')
alt_weight=score(X,0,'VST_rank_weight0')
panels=hall_focus[:4]+['REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX','REACTOME_SIGNALING_BY_BMP','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX']
save(pd.DataFrame({'pathway':panels,'eligible':[p in eligible for p in panels]}),'C_predefined_display_sets.csv')
save(primary[primary.pathway.isin(panels)],'C_display_sample_scores.csv')
statsrows=[];stability=[];summary=[]
tables={name:t.pivot(index='pathway',columns='sample_id',values='ssGSEA_normalized_score').reindex(columns=SAMPLES) for name,t in [('primary',primary),('logcounts',alt_log),('weight0',alt_weight)]}
for p,vals in tables['primary'].iterrows():
    a=vals.to_numpy().reshape(4,3)
    for it,t in enumerate(TIMES):summary.append({'pathway':p,'library':library(p),'time_h':t,'n':3,'mean':a[it].mean(),'sd':a[it].std(ddof=1)})
    for cn in CONTRASTS:
        late,early=map(int,cn.split('_vs_'));delta=a[TIMES.index(late)]-a[TIMES.index(early)]
        st,pv=stats.ttest_1samp(delta,0)
        row={'pathway':p,'library':library(p),'contrast':cn,'mean_paired_delta':delta.mean(),'sd_paired_delta':delta.std(ddof=1),'n_experiments':3,'n_positive':int((delta>0).sum()),'n_negative':int((delta<0).sum()),'t':st,'df':2,'pvalue':pv}
        for j,d in enumerate(delta):row[f'delta_R{j+1}']=d
        statsrows.append(row)
    for tag in ['logcounts','weight0']:
        b=tables[tag].loc[p].to_numpy().reshape(4,3)
        rs={'pathway':p,'library':library(p),'sensitivity':tag,'spearman_12_samples':stats.spearmanr(a.ravel(),b.ravel()).statistic}
        for cn in CONTRASTS:
            late,early=map(int,cn.split('_vs_'));i,j=TIMES.index(late),TIMES.index(early)
            rs[f'mean_direction_agrees_{cn}']=bool(np.sign((a[i]-a[j]).mean())==np.sign((b[i]-b[j]).mean()))
        stability.append(rs)
q=pd.DataFrame(statsrows);q['BH_FDR_library_contrast']=q.groupby(['library','contrast']).pvalue.transform(lambda x:bh(x))
save(q,'C_ssGSEA_paired_contrasts_all_sets.csv');save(pd.DataFrame(summary),'C_ssGSEA_time_summary.csv');save(pd.DataFrame(stability),'C_score_stability_all_sets.csv')

# Every eligible LRT-significant secreted-protein-encoding gene is included.
uni=pd.read_csv(input_paths[4],sep='\t').fillna('')
secreted={}
for _,row in uni.iterrows():
    for symbol in str(row['Gene Names (primary)']).split('; '):
        if symbol:secreted.setdefault(symbol,[]).append(row['Entry'])
lrt=pd.read_csv(input_paths[0]);lrt['secreted_annotation']=lrt.gene_name.isin(secreted)
lrt['UniProt_accessions']=lrt.gene_name.map(lambda x:';'.join(secreted.get(x,[])))
lrt['eligible_dynamic_secreted']=lrt.secreted_annotation & (lrt.padj<0.05) & lrt.gene_id.isin(v.index)
save(lrt,'D_all_tested_genes_secreted_annotation_and_LRT.csv')
d=lrt[lrt.eligible_dynamic_secreted].copy().sort_values('gene_id').reset_index(drop=True)
Y=v.loc[d.gene_id,SAMPLES].to_numpy(float);sd=Y.std(axis=1,ddof=1)
assert np.isfinite(Y).all() and (sd>0).all()
Z=(Y-Y.mean(axis=1,keepdims=True))/sd[:,None];T=Z.reshape(-1,4,3).mean(axis=2)
tree=linkage(T,method='ward',optimal_ordering=True)
select=[]
for k in range(2,min(8,len(T)-1)+1):
    labs=fcluster(tree,k,criterion='maxclust')
    select.append({'k':k,'actual_clusters':len(set(labs)),'silhouette':silhouette_score(T,labs),'minimum_cluster_size':int(pd.Series(labs).value_counts().min())})
ks=pd.DataFrame(select).sort_values(['silhouette','k'],ascending=[False,True]);save(ks,'D_clustering_k_selection.csv')
k=int(ks.iloc[0].k);lab=fcluster(tree,k,criterion='maxclust')
order=sorted(set(lab),key=lambda x:(int(T[lab==x].mean(axis=0).argmax()),tuple(T[lab==x].mean(axis=0))))
labelmap={old:f'S{i+1}' for i,old in enumerate(order)}
d['secreted_module']=[labelmap[x] for x in lab]
save(d,'D_complete_dynamic_secreted_membership.csv')
hm=d[['gene_id','gene_name','secreted_module','padj']].copy()
hm[SAMPLES]=Z
leaves=leaves_list(tree);rank={i:j for j,i in enumerate(leaves)}
hm['leaf_order']=[rank[i] for i in range(len(hm))]
hm=hm.sort_values(['secreted_module','leaf_order']);save(hm,'D_heatmap_all_eligible_genes_Z.csv')
save(pd.DataFrame(tree,columns=['child1','child2','distance','n_members']),'D_Ward_linkage.csv')
boot=[]
for leave in range(3):
    profile=Z.reshape(-1,4,3)[:,:,[r for r in range(3) if r!=leave]].mean(axis=2)
    alt=fcluster(linkage(profile,method='ward'),k,criterion='maxclust')
    boot.append({'omitted_experiment':leave+1,'fixed_k':k,'adjusted_rand_index':adjusted_rand_score(lab,alt)})
save(pd.DataFrame(boot),'D_cluster_leave_one_experiment_stability.csv')
gene_time=[];directions=[];modules=[]
for i,row in d.iterrows():
    a=Z[i].reshape(4,3)
    for it,t in enumerate(TIMES):gene_time.append({'gene_id':row.gene_id,'gene_name':row.gene_name,'module':row.secreted_module,'time_h':t,'mean_Z':a[it].mean(),'sd_Z':a[it].std(ddof=1),'R1':a[it,0],'R2':a[it,1],'R3':a[it,2]})
    for cn in CONTRASTS:
        late,early=map(int,cn.split('_vs_'));delta=a[TIMES.index(late)]-a[TIMES.index(early)]
        directions.append({'gene_id':row.gene_id,'gene_name':row.gene_name,'module':row.secreted_module,'contrast':cn,'mean_delta_Z':delta.mean(),'n_positive':int((delta>0).sum()),'n_negative':int((delta<0).sum())})
for m in sorted(d.secreted_module.unique()):
    a=Z[d.secreted_module.eq(m)].mean(axis=0).reshape(4,3)
    for it,t in enumerate(TIMES):
        for j in range(3):modules.append({'module':m,'n_genes':int(d.secreted_module.eq(m).sum()),'time_h':t,'replicate':j+1,'module_score':a[it,j]})
save(pd.DataFrame(gene_time),'D_gene_timepoint_mean_SD_and_replicates.csv');save(pd.DataFrame(directions),'D_gene_direction_consistency.csv');save(pd.DataFrame(modules),'D_module_sample_scores.csv')
save(lrt[lrt.gene_name.isin(['IL6','CXCL1','CXCL8','CSF3'])],'D_validated_marker_eligibility.csv')

tgf=[p for p in sets if any(x in p for x in ['TGF_BETA','TGFB_FAMILY','TGFBR'])]
bmp=[p for p in sets if 'BMP' in p]
overlaps=[]
id_symbol=gene_map.set_index('ENTREZID').SYMBOL.to_dict()
for a,b in itertools.product(tgf,bmp):
    A,B=set(sets[a]),set(sets[b]);common=A&B; tested=common&universe
    overlaps.append({'TGF_set':a,'BMP_set':b,'TGF_size':len(A),'BMP_size':len(B),'n_shared_full':len(common),'Jaccard_full':len(common)/len(A|B),'n_shared_tested':len(tested),'shared_Entrez':';'.join(sorted(common)),'shared_tested_symbols':';'.join(sorted({str(id_symbol.get(i,i)) for i in tested}))})
save(pd.DataFrame(overlaps),'TGF_BMP_all_set_overlap_audit.csv')
gsum=g.groupby(['library','contrast']).agg(tested_pathways=('pathway','size'),FDR05=('padj_library_contrast',lambda x:int((x<0.05).sum()))).reset_index()
save(gsum,'B_full_library_test_counts.csv')
config={'figure':'S3B-D','samples':SAMPLES,'donors':1,'independent_experiments':3,'no_DESeq2_refit':True,'no_fgsea_rerun':True,'gene_count':len(v),'LRT_FDR05':int((lrt.padj<.05).sum()),'ssGSEA_mapped_genes':len(X),'ssGSEA_sets':len(eligible),'GSEA_display_order':display,'ssGSEA_display_order':panels,'dynamic_secreted_genes':len(d),'secreted_annotated_tested':int(lrt.secreted_annotation.sum()),'secreted_k':k,'silhouette':float(ks.iloc[0].silhouette),'leave_one_out_ARI':[x['adjusted_rand_index'] for x in boot],'gseapy_version':gseapy.__version__,'python':platform.python_version(),'focus_selector':'gene-set names only; not significance','ssgsea_parameters':{'sample_norm_method':'rank','correl_norm_type':'rank','weight':0.25,'permutations':0,'min_size':15,'max_size':500},'secreted_annotation_sha256':hashfile(input_paths[4]),'original_inputs_unchanged':all(hashfile(ROOT/x['file'])==x['sha256'] for x in manifest)}
(BASE/'analysis_summary.json').write_text(json.dumps(config,indent=2),encoding='utf-8')
print(json.dumps(config,indent=2),flush=True)
