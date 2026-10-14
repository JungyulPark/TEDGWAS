from pathlib import Path
import gzip,json,urllib.request,time,sys
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TED_Data_Review_20261002';OUT.mkdir(parents=True,exist_ok=True)
M=OUT/'provenance';M.mkdir(exist_ok=True)
L=Path('C:/ProjectTEDGWAS')
paths=[*sorted((L/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes').glob('*')),*sorted(L.glob('finngen*.gz')),L/'TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz']
inventory=[]
for p in paths:
 if not p.is_file():continue
 op=gzip.open if p.suffix=='.gz' else open
 with op(p,'rt',encoding='utf-8-sig') as f:
  h=next(f).strip().split('\t');r=next(f).strip().split('\t')
 inventory.append({'file':str(p),'bytes':p.stat().st_size,'mtime_ns':p.stat().st_mtime_ns,'columns':h,'first_row':dict(zip(h,r)),'inspection':'header and first row only'})
(M/'local_file_inventory.json').write_text(json.dumps(inventory,ensure_ascii=False,indent=2),encoding='utf8')
urls={
 'BBJ_study':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018627',
 'BBJ_ancestry':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018627/ancestries',
 'UKB_study':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636',
 'UKB_ancestry':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636/ancestries',
 'FinnGen_R12_GRAVES_OPHT':'https://r12.finngen.fi/pheno/GRAVES_OPHT',
 'FinnGen_R12_E4_GRAVES_OPHT_STRICT':'https://r12.finngen.fi/pheno/E4_GRAVES_OPHT_STRICT',
 'FinnGen_R12_E4_GRAVES_STRICT':'https://r12.finngen.fi/pheno/E4_GRAVES_STRICT',
 'UKB_yaml':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038636/GCST90038636_buildGRCh37.tsv-meta.yaml',
 'BBJ_yaml':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018627/GCST90018627_buildGRCh37.tsv.gz-meta.yaml'
}
results=[]
for name,url in urls.items():
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'TED-local-data-review/1.0'})
  with urllib.request.urlopen(req,timeout=20) as resp:body=resp.read();ctype=resp.headers.get('Content-Type','')
  ext='.json' if 'json' in ctype else '.html' if 'html' in ctype else '.txt'
  dest=M/(name+ext);dest.write_bytes(body)
  results.append({'name':name,'url':url,'bytes':len(body),'saved':dest.name,'status':'retrieved'})
  print(name,len(body),flush=True)
 except Exception as e:
  results.append({'name':name,'url':url,'status':'unavailable','error':str(e)});print(name,type(e).__name__,str(e)[:180],flush=True)
(M/'retrieval_log.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf8')
