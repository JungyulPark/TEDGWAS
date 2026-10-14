from pathlib import Path
import urllib.request,json,sys,hashlib,datetime,shutil
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Extension_20261003'
O.mkdir(exist_ok=True);M=O/'Sources';M.mkdir(exist_ok=True)
plan=R/'outputs/TED_Allele_Integrated_20261003/Research_planning'
P=O/'Protocol';P.mkdir(exist_ok=True)
for p in plan.iterdir():
 if not (P/p.name).exists():shutil.copy2(p,P/p.name)
freeze=json.loads((P/'PROTOCOL_FREEZE.json').read_text())
assert hashlib.sha256((P/'PROTOCOL_v1.md').read_bytes()).hexdigest()==freeze['sha256']
urls={'openapi':'https://gtexportal.org/api/v2/openapi.json','download_metadata':'https://gtexportal.org/api/v2/dataset/openAccessFilesMetadata','cloud_objects':'https://storage.googleapis.com/storage/v1/b/adult-gtex/o?prefix=bulk-qtl/v8/single-tissue-cis-qtl/&maxResults=1000'}
for name,url in urls.items():
 dest=M/(name+'.json')
 if dest.exists():print('cached',name);continue
 try:
  req=urllib.request.Request(url,headers={'User-Agent':'TED-research-allele-audit/1.0'})
  with urllib.request.urlopen(req,timeout=40) as r:data=r.read();status=r.status
  value=json.loads(data);dest.write_bytes(data)
  (M/(name+'.provenance.json')).write_text(json.dumps({'url':url,'status':status,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':hashlib.sha256(data).hexdigest()},indent=2))
  print(name,status,len(data),list(value)[:8],flush=True)
 except Exception as e:print(name,type(e).__name__,str(e),flush=True)
