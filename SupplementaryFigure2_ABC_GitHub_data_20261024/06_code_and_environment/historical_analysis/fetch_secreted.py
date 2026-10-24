import urllib.request,urllib.parse,sys
from pathlib import Path
sys.stdout.reconfigure(encoding='utf-8')
O=Path(__file__).resolve().parent/'cache/references'
query=urllib.parse.urlencode({'query':'(organism_id:9606) AND (reviewed:true) AND (keyword:KW-0964)','format':'tsv','fields':'accession,gene_primary,protein_name,cc_subcellular_location'})
jobs={'uniprot_human_reviewed_secreted.tsv':'https://rest.uniprot.org/uniprotkb/stream?'+query,'nichenetr_master.zip':'https://codeload.github.com/saeyslab/nichenetr/zip/refs/heads/master'}
for name,u in jobs.items():
 try:
  with urllib.request.urlopen(urllib.request.Request(u,headers={'User-Agent':'Mozilla/5.0'}),timeout=45) as r:(O/name).write_bytes(r.read())
  print(name,(O/name).stat().st_size,flush=True)
 except Exception as e:print(name,e,flush=True)
