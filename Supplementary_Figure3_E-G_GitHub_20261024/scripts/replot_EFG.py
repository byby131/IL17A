"""Replot saved E–G values; no new DESeq2, PCA, GSEA or ssGSEA."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager, colors
from matplotlib.ft2font import FT2Font
from matplotlib.text import Text

ROOT=Path(__file__).resolve().parents[1]
parser=argparse.ArgumentParser()
parser.add_argument('--font',type=Path,help='Optional Arial Bold TTF path')
parser.add_argument('--qa-dir',type=Path,help='Optional local QA PDFs; requires audit_panel_alignment on PYTHONPATH')
args=parser.parse_args()
if args.font: font_manager.fontManager.addfont(str(args.font))
fontpath=args.font or Path(font_manager.findfont(font_manager.FontProperties(family='Arial',weight='bold'),fallback_to_default=False))
face=FT2Font(str(fontpath))
if face.family_name!='Arial' or 'Bold' not in face.style_name:
    raise RuntimeError('Actual Arial Bold font required; no substitution allowed.')
FP=font_manager.FontProperties(fname=str(fontpath),size=10)
plt.rcParams.update({'font.family':'Arial','font.weight':'bold','font.size':10,'axes.titlesize':10,'axes.titleweight':'bold','axes.labelsize':10,'axes.labelweight':'bold','xtick.labelsize':10,'ytick.labelsize':10,'legend.fontsize':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'svg.fonttype':'none'})
BLUE='#2867B2'; RED='#CE3B3B'
CMAP=colors.LinearSegmentedColormap.from_list('BlueWhitePink',['#63A5CA','#C1DAEB','#FFFFFF','#F0C3D5','#D45F92'],N=256)
def read(p): return pd.read_csv(ROOT/p)
def finish(fig,name):
    fig.canvas.draw()
    for t in fig.findobj(Text): t.set_fontproperties(FP)
    if args.qa_dir:
        from audit_panel_alignment import require_matplotlib_panel_alignment
        args.qa_dir.mkdir(parents=True,exist_ok=True)
        opts={}
        if name.startswith('Figure_S3G'):
            opts={'axes':fig.axes[:2],'panel_ids':['heatmap','NES'],'row_groups':[['heatmap','NES']], 'exemptions':[{'panels':['heatmap','NES'],'checks':['panel-width'],'reason':'Six sample columns versus one NES axis require different widths; row boundaries align.'}]}
        require_matplotlib_panel_alignment(fig,json_out=str(args.qa_dir/(name+'.alignment.json')),strict=True,**opts)
        fig.savefig(args.qa_dir/(name+'.pdf'))
    fig.savefig(ROOT/'preview'/name,dpi=300,facecolor='white')
    plt.close(fig)

pc=read('source_data/Figure_S3E_PCA.csv')
fig,ax=plt.subplots(figsize=(6.7,3.4))
fig.subplots_adjust(left=.14,right=.98,bottom=.30,top=.85)
for k,label,c in [(0,'Control CM_FB',BLUE),(3,'IL-17A CM_FB',RED)]:
    group=pc.iloc[k:k+3]
    ax.scatter(group.PC1,group.PC2,c=c,s=35,label=label)
    for i,row in enumerate(group.itertuples()):
        offset=(6,5) if k==0 else [(8,6),(-13,-3),(8,-9)][i]
        ax.annotate(str(i+1),(row.PC1,row.PC2),xytext=offset,textcoords='offset points')
ax.set(xlim=(-18,19),ylim=(-8,10),xlabel='PC1 (81.9%)',ylabel='PC2 (6.8%)',title='E  Sample-level PCA')
ax.legend(loc='upper center',bbox_to_anchor=(.5,-.27),ncol=2,frameon=False)
finish(fig,'Figure_S3E_PCA.png')

f=read('source_data/Figure_S3F_pathway_overview.csv')
fig,axs=plt.subplots(2,1,figsize=(7.5,8.2))
fig.subplots_adjust(left=.64,right=.98,bottom=.08,top=.93,hspace=.30)
for ax,collection in zip(axs,['Hallmark','Reactome']):
    sub=f[f.pathway.str.startswith(collection.upper()+'_')]
    y=np.arange(len(sub))
    ax.barh(y,sub.NES,color=[RED if n>0 else BLUE for n in sub.NES],height=.65)
    ax.set_yticks(y,sub.display_label); ax.invert_yaxis()
    ax.set(xlabel='NES',xlim=(-3,3),title=collection)
    ax.axvline(0,c='gray',lw=.6);ax.spines['left'].set_visible(False)
fig.suptitle('F  Pathway overview')
finish(fig,'Figure_S3F_pathway_overview.png')

scores=read('source_data/Figure_S3G_sample_scores.csv')
g=read('source_data/Figure_S3G_GSEA.csv')
order=read('metadata/panel_G_display_order.csv')
z=scores.pivot(index='signature',columns='sample',values='row_z').loc[order.signature,pc['sample']].to_numpy()
fig=plt.figure(figsize=(8,4.4))
ax=fig.add_axes([.43,.29,.25,.52]); dots=fig.add_axes([.76,.29,.19,.52])
im=ax.pcolormesh(np.arange(7)-.5,np.arange(8)-.5,z,cmap=CMAP,vmin=-2,vmax=2,edgecolors='white',linewidth=.45)
ax.set_ylim(6.5,-.5);ax.set_yticks(np.arange(7),order.display_label)
ax.set_xticks(np.arange(6),['1','2','3','1','2','3']);ax.tick_params(length=0)
ax.set_title('Score z')
for t,c in zip(ax.get_xticklabels(),[BLUE]*3+[RED]*3):t.set_color(c)
for spine in ax.spines.values():spine.set_visible(False)
for i,row in enumerate(g.itertuples()):
    c=RED if row.NES>0 else BLUE
    dots.hlines(i,0,row.NES,color=c,lw=.7)
    dots.scatter(row.NES,i,s=35,edgecolor=c,facecolor=c if row.FDR_global<.05 else 'white')
dots.set(ylim=(6.5,-.5),xlim=(-2,2.7),title='NES',xticks=[-2,0,2],yticks=[])
dots.axvline(0,c='gray',lw=.6);dots.spines['left'].set_visible(False)
fig.text(.5,.93,'G  Inflammatory programs',ha='center')
fig.text(.43,.18,'Control CM_FB',color=BLUE,ha='center')
fig.text(.73,.18,'IL-17A CM_FB',color=RED,ha='center')
cax=fig.add_axes([.44,.07,.22,.025]);fig.colorbar(im,cax=cax,orientation='horizontal',ticks=[-2,0,2])
finish(fig,'Figure_S3G_inflammatory_programs.png')
print('Three PNG panels written to preview/. Actual font: Arial Bold, 10 pt.')
