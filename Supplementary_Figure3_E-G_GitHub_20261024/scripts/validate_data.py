"""Validate frozen figure values; no analytical refitting or enrichment."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
def read(p): return pd.read_csv(ROOT/p)
checks=[]
def check(name, condition):
    assert bool(condition), name
    checks.append(name)

meta=read('metadata/samples.csv')
pc=read('source_data/Figure_S3E_PCA.csv')
check('Six unique samples, 3 per condition',len(pc)==6 and pc['sample'].is_unique and list(meta.groupby('condition').size())==[3,3])
check('PCA sample order matches metadata',pc['sample'].tolist()==meta['sample'].tolist())
de=read('data/processed/DESeq2_all_tested_genes.csv')
vst=read('data/processed/VST_expression.csv')
check('16713 tested/VST genes',len(de)==16713 and len(vst)==16713)
genes=read('data/processed/PCA_selected_genes.csv')
check('2000 selected PCA genes',len(genes)==2000 and genes.gene_id.is_unique and set(genes.gene_id)<=set(vst.gene_id))
variance=read('source_data/Figure_S3E_variance_explained.csv')
check('Variance fraction sums to one',np.isclose(variance.variance_fraction.sum(),1))
summary=json.loads((ROOT/'metadata/original_DESeq2_summary.json').read_text())
check('Unaltered Salmon input hash',hashlib.sha256((ROOT/'data/input/Salmon_gene_counts.csv').read_bytes()).hexdigest()==summary['source_sha256']['counts'])
full=read('data/processed/GSEA_all_C_E_results.csv').set_index('pathway')
f=read('source_data/Figure_S3F_pathway_overview.csv')
g=read('source_data/Figure_S3G_GSEA.csv')
for name,table,n in [('F',f,24),('G',g,7)]:
    check(name+' expected unique row count',len(table)==n and table.pathway.is_unique)
    for col in ['NES','FDR_global']:
        check(name+' '+col+' matches frozen full GSEA',np.allclose(table[col],full.loc[table.pathway,col],rtol=1e-12,atol=1e-14))
check('F all pathways retain original significant global FDR',f.FDR_global.lt(.05).all())
audit=json.loads((ROOT/'metadata/pathway_display_clusters.json').read_text())
check('F exact saved display selection and order',f.pathway.tolist()==audit['shown'])
scores=read('source_data/Figure_S3G_sample_scores.csv')
check('G exactly 42 unique signature/sample observations',len(scores)==42 and not scores.duplicated(['signature','sample']).any())
check('G program and sample identifiers consistent',set(scores.signature)==set(g.pathway) and set(scores['sample'])==set(pc['sample']))
for key,group in scores.groupby('signature'):
    expected=(group.ssGSEA_score-group.ssGSEA_score.mean())/group.ssGSEA_score.std(ddof=1)
    check('G row z matches sample SD: '+key,np.allclose(expected,group.row_z,rtol=1e-12,atol=1e-12))
for value,suffix in [('ssGSEA_score','raw_scores'),('row_z','row_z')]:
    matrix=read('source_data/Figure_S3G_'+suffix+'_matrix.csv').set_index('signature')
    expected=scores.pivot(index='signature',columns='sample',values=value).loc[g.pathway,pc['sample']]
    check('G '+suffix+' matrix matches long table',np.allclose(matrix.values,expected.values,rtol=1e-12,atol=1e-12))
ranking=read('data/processed/GSEA_rank_all_valid_genes.csv')
check('Mapped Wald ranking unique and descending',ranking.ENTREZID.is_unique and ranking.stat.is_monotonic_decreasing)
manifest=ROOT/'SHA256SUMS.txt'
if manifest.exists():
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest,filename=line.split('  ',1)
        check('SHA256 '+filename,hashlib.sha256((ROOT/filename).read_bytes()).hexdigest()==digest)
report={'status':'PASS','checks':checks,'sample_count':len(pc),'PCA_genes':len(genes),'panel_F_pathways':len(f),'panel_G_programs':len(g),'panel_G_sample_scores':len(scores),'ranked_genes':len(ranking),'full_GSEA_rows':len(full),'refitting_performed':False}
print(json.dumps(report,indent=2))
# Preserve the distributed validation report when a release manifest is present.
if not manifest.exists():
    (ROOT/'metadata/validation_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
