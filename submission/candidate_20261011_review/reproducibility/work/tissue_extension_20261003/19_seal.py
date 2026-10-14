from pathlib import Path
from zipfile import ZipFile,ZIP_DEFLATED
import json,hashlib,shutil
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Integrated_20261003';V=O/'Verification'
qa={'status':'PASS','review':'Visual inspection of all modified DOCX pages (15 manuscript, 9 supplement, 1 checklist), final changed manuscript pages 2/5/6/12, the complete new forest plot, and 13 workbook range previews. No overlapping or clipped content observed. Other manuscript pages were pixel-identical after final wording changes.','main_tables_figures':'Byte-identical to the prior integrated package; not regenerated.','supplement':'Tables S1/S2 and embedded figures preserved; S3 replaced and inspected on page 7.','new_plot_role':'Author-review figure, not an additional numbered manuscript figure.'}
(V/'visual_QA.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
shutil.copy2(W/'19_seal.py',V/'Tissue_extension_scripts/19_seal.py')
files=sorted(p for p in O.rglob('*') if p.is_file() and p.name!='MANIFEST_SHA256.json')
manifest={str(p.relative_to(O)).replace('\\','/'):{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in files}
(O/'MANIFEST_SHA256.json').write_text(json.dumps({'package':'TED_Tissue_Integrated_20261003','files':manifest},indent=2),encoding='utf-8')
zpath=R/'outputs/TED_Tissue_Integrated_20261003.zip'
with ZipFile(zpath,'w',ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,O.name+'/'+p.relative_to(O).as_posix())
with ZipFile(zpath) as z:
 assert z.testzip() is None
 for n,meta in manifest.items():assert hashlib.sha256(z.read(O.name+'/'+n)).hexdigest()==meta['sha256'],n
summary={'zip':str(zpath),'zip_bytes':zpath.stat().st_size,'package_files':len(manifest)+1,'sha256':hashlib.sha256(zpath.read_bytes()).hexdigest(),'all_zip_contents_hash_verified':True}
(W/'delivery_checks.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
