from pathlib import Path
import hashlib,json
W=Path(__file__).resolve().parent
result={}
for d in (W/'render_release').iterdir():
 if not d.is_dir():continue
 old=W/('render_final' if d.name=='Manuscript' else 'render')/d.name
 pages=[]
 for p in sorted(d.glob('page-*.png'),key=lambda x:int(x.stem.split('-')[-1])):
  q=old/p.name
  identical=q.exists() and hashlib.sha256(p.read_bytes()).digest()==hashlib.sha256(q.read_bytes()).digest()
  pages.append({'page':p.name,'identical_to_reviewed':identical})
 result[d.name]=pages
(W/'render_comparison.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps({k:[x['page'] for x in v if not x['identical_to_reviewed']] for k,v in result.items()},indent=2))
