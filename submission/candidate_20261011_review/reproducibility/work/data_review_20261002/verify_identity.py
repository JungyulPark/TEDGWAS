from pathlib import Path
import urllib.request,json,hashlib,re
ROOT=Path(__file__).resolve().parents[2];M=ROOT/'outputs/TED_Data_Review_20261002/provenance'
for code in ['GRAVES_OPHT','E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
 url='https://r12.finngen.fi/api/pheno/'+code
 try:
  with urllib.request.urlopen(url,timeout=25) as r:b=r.read()
  j=json.loads(b);(M/(code+'_api.json')).write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf8')
  print(code,json.dumps(j,ensure_ascii=False)[:3000],flush=True)
 except Exception as e:print(code,str(e)[:200],flush=True)
rows=[]
for cohort,acc in [('BBJ','GCST90018627'),('UKB','GCST90038636')]:
 p=Path('C:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes')/(acc+'_harmonised.tsv.gz')
 h=hashlib.md5()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 s=(M/(cohort+'_harmonized_yaml.txt')).read_text(encoding='utf8');exp=re.search(r'data_file_md5sum:\s*(\w+)',s).group(1)
 row={'cohort':cohort,'path':str(p),'local_md5':h.hexdigest(),'official_md5':exp,'exact_file_match':h.hexdigest()==exp};rows.append(row);print(row,flush=True)
(M/'file_identity.json').write_text(json.dumps(rows,indent=2),encoding='utf8')
