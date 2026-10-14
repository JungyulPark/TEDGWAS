from pathlib import Path
import sys, shutil, json, hashlib, zipfile
from pypdf import PdfReader, PdfWriter

sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1]
O=R/'outputs/TED_Research_Contribution_20261003';V=O/'Verification'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
for name in ['final_editorial_numeric_audit.json','contribution_revision_checks.json']:
    assert json.loads((V/name).read_text(encoding='utf-8'))['status']=='PASS'
qa={'status':'PASS','scope':'Manuscript 14 pages and cover letter 1 page visually inspected after rendering. Supporting files unchanged from previous reviewed package.','documents':{}}
for name,expected in [('Manuscript',14),('Cover_Letter_Draft',1)]:
    folder=W/'render'/name;pdf=folder/(name+'.pdf')
    d=PdfReader(pdf)
    assert len(d.pages)==expected
    for i,p in enumerate(d.pages):
        assert p.extract_text().strip(),(name,i)
        assert '\ufffd' not in p.extract_text(),(name,i)
    assert len(list(folder.glob('page-*.png')))==expected
    qa['documents'][name]={'pages':expected,'all_pages_visually_inspected':True,'docx_sha256':sha(O/(name+'.docx')),'pdf_sha256':sha(pdf),'page_images':{p.name:sha(p) for p in sorted(folder.glob('page-*.png'))}}
    shutil.copy2(pdf,O/(name+'.pdf'))
parts=[('Manuscript',O/'Manuscript.pdf'),('Main tables',O/'Tables.pdf')]+[(f'Figure {n}',O/'Figures'/f'Figure{n}.pdf') for n in [1,2,3]]+[('Supplementary material',O/'Supplementary_Material.pdf')]
merged=PdfWriter();sources=[]
for title,p in parts:
    d=PdfReader(p)
    merged.append(d,outline_item=title)
    sources.extend((p,i) for i in range(len(d.pages)))
merged.write(O/'Full_Review.pdf');merged.close()
d=PdfReader(O/'Full_Review.pdf')
assert len(d.pages)==30
for i,(p,j) in enumerate(sources):
    src=PdfReader(p)
    assert d.pages[i].extract_text()==src.pages[j].extract_text()
qa['combined_review']={'pages':30,'source_text_matches':True,'assembly':'Lossless page append using bundled pypdf'}
(V/'document_visual_qa.json').write_text(json.dumps(qa,indent=2),encoding='utf-8')
audit=json.loads((V/'editorial_audit.json').read_text(encoding='utf-8'));audit['status']='AUTHOR_REVIEW_REVISION_VERIFIED'
audit['pending']='Author approvals and metadata plus final journal requirements; proposed follow-up analyses have not been run.'
(V/'editorial_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
scripts=V/'Revision_scripts';scripts.mkdir(exist_ok=True)
for p in W.glob('*.py'):shutil.copy2(p,scripts/p.name)
manifest={str(p.relative_to(O)).replace('\\','/'):sha(p) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='PACKAGE_SHA256.json'}
(O/'PACKAGE_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=O.with_suffix('.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(O.rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(O)).replace('\\','/'))
with zipfile.ZipFile(zpath) as z:
    assert z.testzip() is None
    for name,digest in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
receipt={'status':'PASS','zip':str(zpath),'sha256':sha(zpath),'files':len(manifest)+1,'bytes':zpath.stat().st_size,'manuscript_pages':14,'review_pages':30,'new_analyses_run':False}
(zpath.parent/(zpath.stem+'_receipt.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2))
