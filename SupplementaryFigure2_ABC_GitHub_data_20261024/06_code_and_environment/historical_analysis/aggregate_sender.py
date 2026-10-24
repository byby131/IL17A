from pathlib import Path
import sys,re,gzip,json,hashlib
import pandas as pd
import numpy as np
sys.stdout.reconfigure(encoding='utf-8')
O=Path(__file__).resolve().parent
f=O.parents[1]/'adsc_il17_rnaseq/early_timecourse/data/ADSC_IL17A_timecourse_RNAseq_EXPRESSION_MATRIX_COUNTS.CSV'
d=pd.read_csv(f,index_col=0);assert d.index.is_unique and np.isfinite(d.values).all() and (d.values>=0).all()
ann=[]
with gzip.open(O/'cache/references/gencode.v36.annotation.gtf.gz','rt') as g:
 for line in g:
  if line.startswith('#'):continue
  cols=line.rstrip().split('\t')
  if cols[2]!='transcript':continue
  a=dict(re.findall(r'(\w+) "([^"]+)"',cols[8]));ann.append([a['transcript_id'],a['gene_id'],a['gene_name'],a['gene_type']])
a=pd.DataFrame(ann,columns=['transcript_id','gene_id','gene_name','gene_type']).set_index('transcript_id')
assert a.index.is_unique
missing=d.index.difference(a.index)
if len(missing):
 print('Unmatched versioned IDs:',len(missing),missing[:12].tolist());raise RuntimeError('Exact reference compatibility not established')
a=a.loc[d.index];groups=a.gene_id
gene=d.groupby(groups,sort=True).sum();assert np.allclose(gene.sum().values,d.sum().values,rtol=1e-12)
gene.to_csv(O/'cache/sender_gene_counts.csv.gz')
ga=a.reset_index(drop=True).drop_duplicates();assert ga.gene_id.is_unique
ga.to_csv(O/'cache/sender_gene_annotation.csv',index=False)
meta=[]
for s in d.columns:
 m=re.fullmatch(r'S50_(0|3|12|24)_([123])',s);assert m
 meta.append({'sample':s,'time':m[1],'replicate':m[2]})
assert len(meta)==12
report={'raw_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'reference':'GENCODE v36 CHR GRCh38.p13','reference_sha256':hashlib.sha256((O/'cache/references/gencode.v36.annotation.gtf.gz').read_bytes()).hexdigest(),'transcripts':len(d),'exact_versioned_matches':len(a),'reference_transcripts':len(ann),'gene_count':len(gene),'count_conservation':True,'metadata':meta,'design':'~ replicate + time','reference_level':'0','contrasts':['3 vs 0','12 vs 0','24 vs 0'],'sender_contrast':'12 vs 0','limitations':['Salmon effective-length offsets unavailable; counts-only gene aggregation is provisional.','0 h is pre-treatment baseline; no time-matched untreated group to isolate culture-time effects.']}
(O/'cache/sender_input_audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print('Exact versioned transcript matches:',len(a),'/',len(d),'genes:',len(gene),'count totals conserved.')
