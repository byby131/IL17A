"""Layout-only reconstruction of frozen S3B-D; A4 top20% reserved for future fluorescence."""
from pathlib import Path
import sys,json,hashlib,shutil,importlib.util
import numpy as np,pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import fontManager,FontProperties
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch
B=Path(__file__).resolve().parent;ROOT=B.parents[1];PRE=ROOT/'analysis/Figure4_S3_submission_QC_20261009';S=B/'source_data';O=B/'figures';Q=B/'qa'
for p in [S,O,Q]:p.mkdir(parents=True,exist_ok=True)
files=['B_display_GSEA.csv','S3C_whole_pathway_direction_composition_all12sets_all3times.csv','S3C_full_pathway_gene_effects_and_Group5_all.csv','S3C_leading_edge_summary_all.csv','D_heatmap_all_eligible_genes_Z.csv','D_module_sample_scores.csv']
for n in files:
    if not(S/n).exists():shutil.copy2(PRE/'source_data'/n,S/n)
for n in ['analysis_summary.json','retained_B_D_config.json']:
    if not(B/n).exists():shutil.copy2(PRE/n,B/n)
if not (B/'qa_helpers/audit_panel_alignment.py').exists():shutil.copytree(PRE/'qa_helpers',B/'qa_helpers',dirs_exist_ok=True)
sys.path.insert(0,str(B/'qa_helpers'));from audit_panel_alignment import require_matplotlib_panel_alignment
f=Path(r'C:\Windows\Fonts\arialbd.ttf');assert f.exists();fontManager.addfont(str(f));assert FontProperties(fname=str(f)).get_name()=='Arial'
mpl.rcParams.update({'font.family':'Arial','font.size':10,'font.weight':'bold','axes.labelweight':'bold','axes.titleweight':'bold','axes.titlesize':10,'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','axes.linewidth':.5,'axes.spines.top':False,'axes.spines.right':False})
fig=plt.figure(figsize=(210/25.4,297/25.4));W,H=210,297
def ax(x,y,w,h):return fig.add_axes([x/W,1-(y+h)/H,w/W,h/H])
def text(x,y,s,**kw):return fig.text(x/W,1-y/H,s,va='top',**kw)
def clean(a):a.tick_params(length=2,pad=1);a.spines['top'].set_visible(False);a.spines['right'].set_visible(False)
BLUE='#2166AC';RED='#B2182B';ORANGE='#D66F4F';TEAL='#158EAD';GRAY='#AEB5BA';cmap=LinearSegmentedColormap.from_list('BWR',[BLUE,'white',RED])
cfg=json.loads((B/'analysis_summary.json').read_text());old=json.loads((B/'retained_B_D_config.json').read_text());samples=old['samples']
names={
'HALLMARK_INFLAMMATORY_RESPONSE':'H | Inflammatory response',
'HALLMARK_TNFA_SIGNALING_VIA_NFKB':'H | TNFA signaling via NF-kB',
'HALLMARK_IL6_JAK_STAT3_SIGNALING':'H | IL6 / JAK / STAT3 signaling',
'HALLMARK_TGF_BETA_SIGNALING':'H | TGF-beta signaling',
'HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION':'H | Epithelial-mesenchymal transition',
'REACTOME_SIGNALING_BY_BMP':'R | BMP signaling',
'REACTOME_SIGNALING_BY_TGFB_FAMILY_MEMBERS':'R | TGFB family signaling',
'REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX':'R | TGF-beta receptor complex',
'REACTOME_TGF_BETA_RECEPTOR_SIGNALING_ACTIVATES_SMADS':'R | TGF-beta receptor activates SMADs',
'REACTOME_DOWNREGULATION_OF_TGF_BETA_RECEPTOR_SIGNALING':'R | TGF-beta receptor downregulation',
'REACTOME_INTERLEUKIN_6_FAMILY_SIGNALING':'R | IL6 family signaling',
'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION':'R | Extracellular matrix organization',
'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX':'R | Extracellular matrix degradation',
'REACTOME_CELL_EXTRACELLULAR_MATRIX_INTERACTIONS':'R | Cell-extracellular matrix interactions',
'REACTOME_COLLAGEN_FORMATION':'R | Collagen formation',
'REACTOME_COLLAGEN_DEGRADATION':'R | Collagen degradation',
'REACTOME_ADORA2B_MEDIATED_ANTI_INFLAMMATORY_CYTOKINES_PRODUCTION':'R | ADORA2B anti-inflammatory cytokines',
'REACTOME_ANTI_INFLAMMATORY_RESPONSE_FAVOURING_LEISHMANIA_PARASITE_INFECTION':'R | Anti-inflammatory / Leishmania',
'REACTOME_FCERI_MEDIATED_NF_KB_ACTIVATION':'R | FCERI-mediated NF-kB activation',
'REACTOME_INFLAMMASOMES':'R | Inflammasomes',
'REACTOME_RIP_MEDIATED_NFKB_ACTIVATION_VIA_ZBP1':'R | RIP / ZBP1-mediated NF-kB activation',
'REACTOME_TRAF6_MEDIATED_NF_KB_ACTIVATION':'R | TRAF6-mediated NF-kB activation'}
# B: all22 preselected rows and every original contrast retained.
text(5,62,'B');text(11,62,'Time-resolved pathway enrichment')
g=pd.read_csv(S/'B_display_GSEA.csv');order=old['GSEA_display_order'];cons=['3_vs_0','12_vs_0','24_vs_0','24_vs_3']
nes=g.pivot(index='pathway',columns='contrast',values='NES').loc[order,cons];qs=g.pivot(index='pathway',columns='contrast',values='padj_library_contrast').loc[order,cons]
a=ax(101,69,104,87);im=a.imshow(nes,aspect='auto',cmap=cmap,vmin=-3,vmax=3,interpolation='nearest')
a.set_yticks(range(22),[names[p] for p in order]);a.set_xticks(range(4),['3 vs 0','12 vs 0','24 vs 0','24 vs 3']);a.tick_params(length=0,pad=2)
a.set_xticks(np.arange(-.5,4,1),minor=True);a.set_yticks(np.arange(-.5,22,1),minor=True);a.grid(which='minor',color='white',lw=.5);a.tick_params(which='minor',length=0)
for i in range(22):
    for j in range(4):
        q=qs.iloc[i,j];s=f'{q:.0e}' if q<.001 else f'{q:.3f}';s+='*' if q<.05 else ''
        a.text(j,i,s,ha='center',va='center',color='white' if abs(nes.iloc[i,j])>1.95 else '#17232D')
for sp in a.spines.values():sp.set_visible(False)
ca=ax(10,164,28,2);cb=fig.colorbar(im,cax=ca,orientation='horizontal',ticks=[-3,0,3]);cb.ax.tick_params(length=1,pad=1)
text(46,163,'NES (color), BH-FDR (text); * FDR < 0.05')
text(46,167,'H: Hallmark; R: Reactome; contrasts in hours')
# C: four pathways as columns. Each retains composition / all gene effects / all LE entries.
text(5,175,'C');text(11,175,'Pathway heterogeneity')
paths=cfg['S3C_selected_sets'];short=['Hallmark\nTGF-beta','Reactome\nBMP','Reactome\nTGF-beta / SMADs','Reactome\nTGFB family']
cp=pd.read_csv(S/'S3C_whole_pathway_direction_composition_all12sets_all3times.csv');ge=pd.read_csv(S/'S3C_full_pathway_gene_effects_and_Group5_all.csv');le=pd.read_csv(S/'S3C_leading_edge_summary_all.csv')
comp_axes=[];dist_axes=[];rng=np.random.default_rng(20261008)
for i,p in enumerate(paths):
    x=9+i*49;text(x+19,180.5,short[i],ha='center',linespacing=1.2)
    sub=cp[cp.pathway==p];count=sub.pivot(index='contrast',columns='status',values='count');N=int(sub.total_detectable.iloc[0])
    a=ax(x,192,23,12);comp_axes.append(a);base=np.zeros(3)
    for st,col in [('Up',ORANGE),('Down',BLUE),('Not significant','#D5D9DC')]:
        z=sub[sub.status==st].set_index('contrast').loc[cons[:3]];a.barh(range(3),z.proportion,left=base,height=.72,color=col,edgecolor='white',lw=.3);base+=z.proportion.to_numpy()
    a.set_yticks(range(3),['3','12','24']);a.invert_yaxis();a.set_xlim(0,1);a.set_xticks([0,1],['0','100']);a.spines['left'].set_visible(False);clean(a);a.tick_params(axis='y',length=0)
    for j,cn in enumerate(cons[:3]):a.text(1.10,j,'/'.join(str(int(count.loc[cn,s])) for s in ['Up','Down','Not significant']),va='center',clip_on=False)
    text(x+18,213,f'G5: {int(ge[(ge.pathway==p)&(ge.contrast=="24_vs_0")].in_Group5.sum())}/{N}',ha='center')
    m=ax(x+3,221,34,16);dist_axes.append(m);subg=ge[(ge.pathway==p)&ge.contrast.isin(cons[1:3])];arrays=[]
    for j,(cn,flag) in enumerate([('12_vs_0',True),('12_vs_0',False),('24_vs_0',True),('24_vs_0',False)]):
        vals=subg[(subg.contrast==cn)&(subg.in_Group5==flag)].log2FoldChange.to_numpy();arrays.append(vals)
        # x-only jitter displays existing observations; expression values never simulated.
        m.scatter(j+rng.uniform(-.13,.13,len(vals)),vals,s=4,c=TEAL if flag else GRAY,alpha=.8,edgecolors='none',zorder=3)
    box=m.boxplot(arrays,positions=range(4),widths=.55,showfliers=False,patch_artist=True,medianprops={'color':'#333','linewidth':.6},whiskerprops={'linewidth':.4},capprops={'linewidth':.4})
    for j,b in enumerate(box['boxes']):b.set(facecolor=TEAL if j%2==0 else GRAY,alpha=.20)
    m.axhline(0,color='#89939B',ls='--',lw=.4);m.set_ylim(-3,2.2);m.set_yticks([-2,0,2]);m.set_xlim(-.6,3.6);m.set_xticks([.5,2.5],['12 h','24 h']);clean(m)
    strings=[]
    for cn in cons[:3]:
        r=le[(le.pathway==p)&(le.contrast==cn)].iloc[0];strings.append(('+' if r.NES>0 else '−')+f'{int(r.n_Group5)}/{int(r.n_leading_edge)}')
    text(x+19,249,'  '.join(strings),ha='center')
text(7,209,'Up',color=ORANGE);text(14,209,'/');text(17,209,'Down',color=BLUE);text(28,209,'/');text(31,209,'NS',color='#66717A');text(39,209,': counts; bars = %; rows = hours')
text(113,175,'G5',color=TEAL);text(121,175,'/ other',color='#66717A');text(141,175,'log2FC vs 0 h')
text(7,244,'LE: 3 / 12 / 24 h; NES sign, G5/LE; all FDR ≥ 0.05')
# D: retain all500 genes /12 samples and the two individual-experiment trajectories.
text(5,256,'D');text(11,256,'Dynamic secreted-protein-encoding genes')
hm=pd.read_csv(S/'D_heatmap_all_eligible_genes_Z.csv');mods=pd.read_csv(S/'D_module_sample_scores.csv');assert len(hm)==500
d=ax(16,264,46,22);d.imshow(hm[samples],aspect='auto',cmap=cmap,vmin=-2.5,vmax=2.5,interpolation='nearest')
d.set_xticks(range(12),[f'R{r}' for t in range(4) for r in (1,2,3)],rotation=90);d.tick_params(length=0,pad=1)
d.set_yticks([168,418],['S1','S2']);d.axhline(336.5,color='white',lw=.6)
for j in [2.5,5.5,8.5]:d.axvline(j,color='white',lw=.5)
for j,t in enumerate([0,3,12,24]):text(16+(j*3+1.5)/12*46,260,f'{t} h',ha='center')
for sp in d.spines.values():sp.set_visible(False)
positions=sorted((np.flatnonzero(hm.gene_name.eq(n))[0],n) for n in ['IL6','CXCL8','CXCL1','CSF3'])
targets=[]
for k,(true,n) in enumerate(positions):
    yy=max(true,targets[-1]+85) if targets else true;targets.append(yy)
if targets[-1]>490:targets=[t-(targets[-1]-490) for t in targets]
for (true,n),yy in zip(positions,targets):
    d.plot([11.5,12.1,12.4],[true,true,yy],color='#56616C',lw=.4,clip_on=False);d.text(12.6,yy,n,va='center',clip_on=False)
d.set_xlim(-.5,11.5);d.set_ylim(499.5,-.5)
curves=[]
for i,mod in enumerate(['S1','S2']):
    a=ax(89+i*62,267,52,18);curves.append(a);dat=mods[mods.module==mod];color=ORANGE if i==0 else TEAL
    for r,mark in [(1,'o'),(2,'s'),(3,'^')]:
        z=dat[dat.replicate==r].sort_values('time_h');a.plot(z.time_h,z.module_score,c=color,alpha=.3,lw=.4);a.scatter(z.time_h,z.module_score,s=8,marker=mark,facecolor='white',edgecolor=color,lw=.5,zorder=3)
    su=dat.groupby('time_h').module_score.agg(['mean','std']).reindex([0,3,12,24]);a.fill_between(su.index,su['mean']-su['std'],su['mean']+su['std'],color=color,alpha=.16,lw=0);a.plot(su.index,su['mean'],c=color,lw=1)
    a.set_xlim(-1,25);a.set_ylim(-1.3,1.3);a.set_yticks([-1,0,1]);a.set_xticks([0,3,12,24]);clean(a);a.set_title(f'{mod} | {int(dat.n_genes.iloc[0])} genes',pad=3)
text(89,291,'Time (h); mean module Z ± SD')
ca=ax(65,294,17,1.5);cb=fig.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(-2.5,2.5),cmap=cmap),cax=ca,orientation='horizontal',ticks=[-2,0,2]);cb.ax.xaxis.set_ticks_position('top');cb.ax.tick_params(length=1,pad=1)
text(57,292,'Z')
panels=comp_axes+dist_axes+curves;ids=[f'c{i}' for i in range(4)]+[f'e{i}' for i in range(4)]+['s1','s2']
for a,k in zip(panels,ids):a.set_label(k)
require_matplotlib_panel_alignment(fig,axes=panels,panel_ids=ids,row_groups=[{'id':'composition','panels':ids[:4]},{'id':'effects','panels':ids[4:8]},{'id':'trajectories','panels':ids[8:]}],json_out=Q/'A4_alignment.json',strict=True)
# Validate exact font sizes and reserved area before exporting, with no tight-bbox cropping.
fig.canvas.draw();renderer=fig.canvas.get_renderer();bounds=[]
for t in fig.findobj(mpl.text.Text):
    if not t.get_visible() or not t.get_text():continue
    assert t.get_fontsize()==10 and t.get_fontweight()=='bold'
    bb=t.get_window_extent(renderer);top=297-bb.y1/fig.dpi*25.4;bottom=297-bb.y0/fig.dpi*25.4
    bounds.append({'text':t.get_text(),'top_mm':top,'bottom_mm':bottom,'left_mm':bb.x0/fig.dpi*25.4,'right_mm':bb.x1/fig.dpi*25.4})
pd.DataFrame(bounds).to_csv(Q/'text_bounds.csv',index=False)
prefix=O/'SupplementaryFigure3_A4_80percent'
fig.savefig(str(prefix)+'.pdf',dpi=600);fig.savefig(str(prefix)+'.svg',dpi=600);fig.savefig(str(prefix)+'.jpg',dpi=300);fig.savefig(str(prefix)+'.tiff',dpi=600,pil_kwargs={'compression':'tiff_lzw'})
pd.DataFrame([{'pathway':p,'display_label':names[p],'size':int(g[g.pathway==p]['size'].iloc[0])} for p in order]).to_csv(S/'B_display_label_key.csv',index=False)
(Q/'layout_contract.json').write_text(json.dumps({'page_mm':[210,297],'reserved_top_mm':59.4,'content_region_mm':[0,59.4,210,237.6],'font':'Arial Bold','font_pt':10,'analysis_changed':False,'B_rows':22,'C_sets':4,'D_genes':500,'D_samples':12},indent=2))
print('A4 rendered; all original data retained.')
