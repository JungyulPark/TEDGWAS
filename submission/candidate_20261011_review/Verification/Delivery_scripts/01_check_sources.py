from pathlib import Path
from pypdf import PdfReader
import hashlib,json,sys
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
reports=[]
for name in ['TED_Tissue_Integrated_20261003','TED_Tissue_Tables_Figures_20261003']:
 root=R/'outputs'/name;m=json.loads((root/'MANIFEST_SHA256.json').read_text(encoding='utf-8'));m=m.get('files',m)
 failures=[]
 for fn,meta in m.items():
  expected=meta['sha256'] if isinstance(meta,dict) else meta
  if not (root/fn).exists() or sha(root/fn)!=expected:failures.append(fn)
 assert not failures,failures
 reports.append({'package':name,'verified_files':len(m),'all_sha256_match':True})
main=R/'outputs/TED_Tissue_Integrated_20261003';p=PdfReader(main/'Supplementary_Material.pdf')
for i in [4,5,6]:
 t=p.pages[i].extract_text();assert f'Table S{i-3}.' in t
 print('Supplementary page',i+1,t[:100].replace('\n',' '))
for p in [main/'Tables.pdf']+list((main/'Figures').glob('*.pdf')):
 rd=PdfReader(p);print(p.name,len(rd.pages),'pages')
(W/'source_verification.json').write_text(json.dumps(reports,indent=2),encoding='utf-8')
print(json.dumps(reports,indent=2))
