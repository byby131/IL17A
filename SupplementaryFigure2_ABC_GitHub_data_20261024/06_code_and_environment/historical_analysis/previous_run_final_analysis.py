from pathlib import Path
import json,shutil,hashlib
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

B=Path(__file__).resolve().parent;ROOT=B.parents[1];S=B/'source_data';PREV=ROOT/'analysis/figure4_S3_refinement_20261008';ORIG=ROOT/'analysis/supplementary_S3_20261008'
SAMPLES=[f'S50_{h}_{r}' for h in (0,3,12,24) for r in (1,2,3)];CONS=['3_vs_0','12_vs_0','24_vs_0']
def save(d,n):d.to_csv(S/n,index=False)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def gmt(p):return {x[0]:set(x[2:]) for x in (l.strip().split('\t') for l in p.open())}
def bh(p):
    p=np.asarray(p,float);idx=np.argsort(p);q=np.empty(len(p));q[idx]=np.minimum(1,np.minimum.accumulate((p[idx]*len(p)/np.arange(1,len(p)+1))[::-1])[::-1]);return q
for p in (PREV/'source_data').iterdir():
    if p.suffix in ['.csv','.tsv','.gmt','.txt']:shutil.copy2(p,S/p.name)
mem=pd.read_csv(S/'Figure4C_frozen_membership_mapping.csv');v=pd.read_csv(S/'VST_all_14425_genes.csv').set_index('gene_id')
U=set(v.index.str.replace(r'\.\d+$','',regex=True));assert len(U)==len(v)==14425
save(pd.DataFrame({'gene_id_unversioned':sorted(U)}),'ORA_common_detectable_universe.csv')
mods={i:set(mem.loc[mem.cluster==i,'gene_id_unversioned'])&U for i in range(1,8)}
ma=[]
for i in mods:ma.append({'group':i,'frozen_total_n':int((mem.cluster==i).sum()),'detectable_query_n':len(mods[i]),'excluded_not_detectable_n':int((mem.cluster==i).sum())-len(mods[i]),'common_universe_n':len(U)})
save(pd.DataFrame(ma),'ORA_module_mapping_audit.csv')
symbol=dict(zip(mem.gene_id_unversioned,mem.symbol));terms={};rows=[];ta=[]
for lib in ['GO_BP','Reactome']:
    database=gmt(S/f'{lib}.2025.1.Hs.ensembl.gmt')
    for name,allids in database.items():
        ids=allids&U;eligible=15<=len(ids)<=500;ta.append({'library':lib,'term_id':name,'database_size':len(allids),'detectable_size':len(ids),'tested':eligible})
        if not eligible:continue
        terms[(lib,name)]=ids
        for group,query in mods.items():
            hits=sorted(ids&query);k=len(hits);n=len(query);K=len(ids);N=len(U)
            rows.append({'group':group,'library':lib,'term_id':name,'N_background':N,'n_module':n,'K_pathway':K,'k_hits':k,'GeneRatio':k/n,'BgRatio':K/N,'fold_enrichment':(k/n)/(K/N),'pvalue':hypergeom.sf(k-1,N,K,n),'hit_ensembl':';'.join(hits),'hit_symbols':';'.join(str(symbol.get(g,g)) for g in hits)})
ora=pd.DataFrame(rows);ora['BH_FDR_global_7modules_2libraries']=bh(ora.pvalue)
ora['BH_FDR_within_module_library_sensitivity']=ora.groupby(['group','library']).pvalue.transform(lambda p:bh(p))
save(ora,'ORA_all7_GO_BP_Reactome_complete.csv');save(pd.DataFrame(ta),'ORA_term_eligibility_all.csv')
selection=[];display=[]
for group in range(1,8):
    sub=ora[ora.group==group].sort_values(['BH_FDR_global_7modules_2libraries','fold_enrichment','k_hits','library','term_id'],ascending=[True,False,False,True,True]);chosen=[];rank=0
    for _,r in sub.iterrows():
        reason='eligible';representative='';jac=np.nan
        if r.BH_FDR_global_7modules_2libraries>=.05:reason='global_FDR_not_significant'
        elif r.fold_enrichment<=1:reason='not_overrepresented'
        elif r.k_hits<3:reason='fewer_than_3_hits'
        else:
            ids=terms[(r.library,r.term_id)]
            for rr in chosen:
                other=terms[(rr['library'],rr['term_id'])];j=len(ids&other)/len(ids|other)
                if j>=.5:reason='redundant_Jaccard_ge_0.5';representative=rr['term_id'];jac=j;break
            if reason=='eligible':
                rank+=1;rr=r.to_dict();rr['nonredundant_rank']=rank;chosen.append(rr)
                if rank<=6:display.append(rr)
        selection.append({'group':group,'library':r.library,'term_id':r.term_id,'decision':reason,'represented_by':representative,'Jaccard':jac,'nonredundant_rank':rank if reason=='eligible' else np.nan,'main_display':reason=='eligible' and rank<=6 and group in [1,5]})
save(pd.DataFrame(selection),'ORA_selection_audit_all7.csv');disp=pd.DataFrame(display);save(disp,'ORA_top6_nonredundant_each_module.csv');save(disp[disp.group.isin([1,5])],'Figure4C_GO_Reactome_display.csv')
summary=ora.groupby(['group','library']).agg(tested=('term_id','size'),significant=('BH_FDR_global_7modules_2libraries',lambda x:(x<.05).sum())).reset_index();save(summary,'ORA_test_counts_summary.csv')

# Existing fixed TGF/BMP definitions, complete exact fgsea ranks and mapping.
sets={**gmt(S/'h.all.v2025.1.Hs.entrez.gmt'),**gmt(S/'c2.cp.reactome.v2025.1.Hs.entrez.gmt')}
ref=pd.read_csv(S/'C_fixed_TGF_BMP_database_membership.csv',dtype=str);names=sorted(ref.gs_name.unique())
gsea=pd.read_csv(S/'B_GSEA_all_Hallmark_Reactome_5_contrasts.csv');ranks=pd.read_csv(S/'B_ranked_gene_backgrounds_and_Wald_stats.csv',dtype={'ENTREZID':str})
group5=set(mem.loc[mem.cluster==5,'gene_id_unversioned']);canonicalrank=ranks[ranks.contrast=='24_vs_0'];crU=set(canonicalrank.ENTREZID)
crmap=dict(zip(canonicalrank.ENTREZID,canonicalrank.gene_id.str.replace(r'\.\d+$','',regex=True)))
elig=[]
for name in names:
    ids=sets[name]&crU;loc={e for e in ids if crmap[e] in group5}
    elig.append({'pathway':name,'library':'Hallmark' if name.startswith('HALLMARK_') else 'Reactome','tested_members':len(ids),'G5_members':len(loc),'G5_fraction':len(loc)/len(ids) if ids else np.nan,'eligible':15<=len(ids)<=500 and len(loc)>=3})
elig=pd.DataFrame(elig);selected=[];rules=[]
for _,r in elig.sort_values(['library','tested_members','pathway']).iterrows():
    reason='selected';other='';jac=np.nan
    if not r.eligible:reason='insufficient_G5_members_or_GSEA_size'
    elif r.library=='Reactome':
        ids=sets[r.pathway]&crU
        for p in selected:
            if not p.startswith('REACTOME'):continue
            b=sets[p]&crU;j=len(ids&b)/len(ids|b)
            if j>=.5:reason='redundant_full_set_Jaccard_ge_0.5';other=p;jac=j;break
    if reason=='selected':selected.append(r.pathway)
    rr=r.to_dict();rr.update({'decision':reason,'represented_by':other,'Jaccard':jac});rules.append(rr)
save(pd.DataFrame(rules),'S3C_representative_selection_audit.csv')
composition=[];gene_rows=[];leads=[];leadsgenes=[];g5membership=[]
for cn in CONS:
    rank=ranks[ranks.contrast==cn].copy();rank['unversioned']=rank.gene_id.str.replace(r'\.\d+$','',regex=True);rank['in_Group5']=rank.unversioned.isin(group5)
    d=pd.read_csv(S/f'DESeq2_{cn}.csv').set_index('gene_id')
    rank['pvalue']=rank.gene_id.map(d.pvalue);rank['gene_BH_FDR']=rank.gene_id.map(d.padj)
    rank['direction']=np.where(rank.gene_BH_FDR.isna(),'Unavailable',np.where((rank.gene_BH_FDR<.05)&(rank.log2FoldChange>0),'Up',np.where((rank.gene_BH_FDR<.05)&(rank.log2FoldChange<0),'Down','Not significant')))
    for name in names:
        members=rank[rank.ENTREZID.isin(sets[name])].copy();members['pathway']=name;gene_rows.append(members)
        N=len(members)
        for status in ['Up','Down','Not significant','Unavailable']:
            n=int((members.direction==status).sum());composition.append({'pathway':name,'contrast':cn,'status':status,'count':n,'total_detectable':N,'proportion':n/N if N else np.nan})
        for flag,scope in [(True,'Group5'),(False,'Outside_Group5')]:
            m=members[members.in_Group5==flag];g5membership.append({'pathway':name,'contrast':cn,'scope':scope,'n':len(m),'fraction_of_detectable_pathway':len(m)/N if N else np.nan,'median_log2FC':m.log2FoldChange.median(),'q25_log2FC':m.log2FoldChange.quantile(.25),'q75_log2FC':m.log2FoldChange.quantile(.75),'n_up_FDR05':int((m.direction=='Up').sum()),'n_down_FDR05':int((m.direction=='Down').sum()),'n_not_significant':int((m.direction=='Not significant').sum())})
        rr=gsea[(gsea.contrast==cn)&(gsea.pathway==name)]
        if not len(rr):
            leads.append({'pathway':name,'contrast':cn,'tested':False,'n_leading_edge':0,'interpretation':'not_tested_under_frozen_size_filter'});continue
        rr=rr.iloc[0];ids=set(str(rr.leadingEdge_Entrez).split(';'));l=members[members.ENTREZID.isin(ids)].copy();assert len(l)==len(ids)
        l['NES']=rr.NES;l['pathway_BH_FDR']=rr.padj_library_contrast;l['leading_edge_type']='positive-enrichment' if rr.NES>0 else 'negative-enrichment';leadsgenes.append(l)
        leads.append({'pathway':name,'contrast':cn,'tested':True,'NES':rr.NES,'pathway_BH_FDR':rr.padj_library_contrast,'n_leading_edge':len(l),'n_Group5':int(l.in_Group5.sum()),'fraction_Group5':l.in_Group5.mean(),'n_positive_log2FC':int((l.log2FoldChange>0).sum()),'n_negative_log2FC':int((l.log2FoldChange<0).sum()),'n_gene_significant_up':int((l.direction=='Up').sum()),'n_gene_significant_down':int((l.direction=='Down').sum()),'leading_edge_type':l.leading_edge_type.iloc[0],'interpretation':'significant_enrichment' if rr.padj_library_contrast<.05 else 'descriptive_only_pathway_FDR_ge_0.05'})
save(pd.DataFrame(composition),'S3C_whole_pathway_direction_composition_all12sets_all3times.csv');save(pd.concat(gene_rows),'S3C_full_pathway_gene_effects_and_Group5_all.csv');save(pd.DataFrame(g5membership),'S3C_Group5_vs_other_gene_effect_summaries.csv');save(pd.DataFrame(leads),'S3C_leading_edge_summary_all.csv');save(pd.concat(leadsgenes),'S3C_leading_edge_genes_all.csv')

# NOG: inspect original count availability and exact prefit filter, never impute.
raw=pd.read_csv(ROOT/'HFF-1/round3_20260914/cache/sender_gene_counts.csv.gz');raw['unversioned']=raw.iloc[:,0].str.replace(r'\.\d+$','',regex=True)
nog=raw[raw.unversioned=='ENSG00000183691'].copy();save(nog,'NOG_original_aggregated_counts.csv')
nog_qc={'rows_found':len(nog),'in_saved_DESeq2':any(x.startswith('ENSG00000183691.') for x in v.index),'reason':'not_found_in_aggregated_matrix'}
if len(nog)==1:
    vals=nog[SAMPLES].iloc[0].to_numpy(float);rounded=np.rint(vals);nog_qc.update({'raw_count_min':float(vals.min()),'raw_count_max':float(vals.max()),'n_samples_rounded_count_ge10':int((rounded>=10).sum()),'passes_original_ge10_in3samples_filter':bool((rounded>=10).sum()>=3),'reason':'fails_original_prefit_expression_filter' if (rounded>=10).sum()<3 else 'requires_further_mapping_audit'})
(B/'NOG_QC.json').write_text(json.dumps(nog_qc,indent=2))
# Original Figure4D is server display data (not formal VST),11 KEGG overlap genes.
fd=ROOT/'analysis/figure4C_group5_preview/source_data'
for nm in ['Figure4C_Group5_TGFbeta_overlap_membership_11.csv','Figure4C_Group5_TGFbeta_overlap_server_expression_12_samples.csv','Figure4C_Group5_TGFbeta_overlap_timepoint_means_wide.csv']:
    shutil.copy2(fd/nm,S/('Figure4D_original_'+nm))
dmem=pd.read_csv(fd/'Figure4C_Group5_TGFbeta_overlap_membership_11.csv');key='symbol' if 'symbol' in dmem.columns else 'SYMBOL';dgenes=set(dmem[key]);ds=[]
for cn in ['3_vs_0','12_vs_0','24_vs_0','12_vs_3','24_vs_3']:
    d=pd.read_csv(S/f'DESeq2_{cn}.csv');d=d[d.gene_name.isin(dgenes)].copy();d['contrast']=cn;ds.append(d)
save(pd.concat(ds),'Figure4D_original11_genes_available_DE_all5contrasts.csv')
save(pd.DataFrame([{'symbol':g,'available_in_fixed_DESeq2':g in set(pd.concat(ds).gene_name)} for g in sorted(dgenes)]),'Figure4D_original11_gene_DE_availability.csv')
cfg={'common_universe_n':len(U),'ORA_tests':len(ora),'ORA_primary_BH_scope':'all7modules_x_GO_BP_and_Reactome','main_Group1_terms':disp[disp.group==1].term_id.tolist(),'main_Group5_terms':disp[disp.group==5].term_id.tolist(),'S3C_selected_sets':selected,'S3C_all_defined_sets':names,'NOG':nog_qc,'original_model_and_membership_unchanged':True,'no_fitting_no_GSEA_rerun':True}
(B/'analysis_summary.json').write_text(json.dumps(cfg,indent=2),encoding='utf8')
print(json.dumps(cfg,indent=2))
