from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import json,hashlib,shutil
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Tables_Figures_20261003';A=R/'outputs/TED_Tissue_Extension_20261003'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for ext in ['pdf','png']:assert sha(O/f'Figures/Tissue_all_estimates.{ext}')==sha(A/f'Figures/Tissue_extension_forest.{ext}')
qa={'status':'PASS','new_tables':'Both rendered landscape pages inspected; all text and 20 estimates legible; no title border, clipping or overlap.','new_summary_figure':'Final PNG inspected, including availability symbols, four primary confidence intervals, raw/adjusted P values and all footnotes.','full_forest':'Previously verified figure copied byte-identically.','analysis':'No estimates changed; all display values checked against canonical tables.'}
(O/'Verification/visual_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8');shutil.copy2(W/'05_seal.py',O/'Verification/05_seal.py')
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.json')
manifest={p.relative_to(O).as_posix():sha(p) for p in files}
(O/'MANIFEST_SHA256.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
zf=R/'outputs/TED_Tissue_Tables_Figures_20261003.zip'
with ZipFile(zf,'w',ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,O.name+'/'+p.relative_to(O).as_posix())
with ZipFile(zf) as z:
 assert z.testzip() is None
 for n,h in manifest.items():assert hashlib.sha256(z.read(O.name+'/'+n)).hexdigest()==h
print(json.dumps({'zip':str(zf),'bytes':zf.stat().st_size,'files':len(manifest)+1,'checks':'PASS'},indent=2))
