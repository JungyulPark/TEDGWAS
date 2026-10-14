from pathlib import Path
import urllib.request,json,time
ROOT=Path(__file__).resolve().parents[2];M=ROOT/'outputs/TED_Data_Review_20261002/provenance'
urls={
 'UKB_study_retry':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636',
 'FinnGen_R12_config':'https://r12.finngen.fi/static/resources/config.js',
 'BBJ_harmonized_directory':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90018001-GCST90019000/GCST90018627/harmonised/',
 'UKB_harmonized_directory':'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90038001-GCST90039000/GCST90038636/harmonised/',
 'UKB_paper_europepmc':'https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=EXT_ID:33959723%20AND%20SRC:MED&format=json',
 'FinnGen_R12_endpoint_definition':'https://r12.risteys.finngen.fi/endpoints/GRAVES_OPHT',
 'FinnGen_R12_strictGD_definition':'https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_STRICT',
 'FinnGen_R12_strictGO_definition':'https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_OPHT_STRICT'
}
logs=[]
for name,url in urls.items():
 try:
  with urllib.request.urlopen(urllib.request.Request(url,headers={'User-Agent':'Mozilla/5.0'}),timeout=25) as r:body=r.read();ctype=r.headers.get('content-type','')
  ext='.json' if 'json' in ctype else '.html' if 'html' in ctype else '.txt';dest=M/(name+ext);dest.write_bytes(body)
  logs.append({'name':name,'url':url,'bytes':len(body),'saved':dest.name,'status':'retrieved'});print(name,len(body),flush=True)
 except Exception as e:logs.append({'name':name,'url':url,'status':'unavailable','error':str(e)});print(name,str(e)[:160],flush=True)
(M/'retrieval_followup.json').write_text(json.dumps(logs,ensure_ascii=False,indent=2),encoding='utf8')
