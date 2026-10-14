from pathlib import Path
import sys,shutil,json,hashlib,zipfile
W=Path(__file__).resolve().parent;R=W.parents[1]
sys.path.insert(0,str(R/'work/pydeps'))
import fitz
sys.stdout.reconfigure(encoding='utf-8')
O=R/'outputs/TED_Association_Manuscript_20261003';V=O/'Verification'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
qa={'status':'PASS','scope':'All affected document pages visually inspected; unchanged supporting files retain previous review and match SHA-256','documents':{}}
for name,expected in [('Manuscript',13),('Cover_Letter_Draft',1),('STROBE_MR_Checklist',1)]:
    folder=W/'render'/name;pdf=folder/(name+'.pdf')
    with fitz.open(pdf) as d:
        assert len(d)==expected,(name,len(d))
        for i,p in enumerate(d):
            assert p.get_text().strip(),(name,i)
            assert '\ufffd' not in p.get_text(),(name,i,'replacement glyph')
    # Count actual PDF pages; an earlier renderer pass left a stale PNG for
    # a removed page of the checklist. That image is not delivered or counted.
    imgs={f'page-{n}.png':sha(folder/f'page-{n}.png') for n in range(1,expected+1)}
    before=W/'render_before_polish'/name
    unchanged=[n for n,h in imgs.items() if (before/n).exists() and sha(before/n)==h]
    qa['documents'][name]={'pages':expected,'all_pages_visually_reviewed':True,'docx_sha256':sha(O/(name+'.docx')),'pdf_sha256':sha(pdf),'page_images':imgs,'unchanged_since_first_visual_pass':unchanged}
    shutil.copy2(pdf,O/(name+'.pdf'))

parts=[('Manuscript',O/'Manuscript.pdf'),('Main tables',O/'Tables.pdf')]+[(f'Figure {n}',O/'Figures'/f'Figure{n}.pdf') for n in [1,2,3]]+[('Supplementary material',O/'Supplementary_Material.pdf')]
merged=fitz.open();toc=[];page_sources=[]
for title,p in parts:
    with fitz.open(p) as d:
        toc.append([1,title,len(merged)+1]);merged.insert_pdf(d)
        page_sources.extend((p,i) for i in range(len(d)))
merged.set_toc(toc);merged.save(O/'Full_Review.pdf',garbage=4,deflate=True);merged.close()
with fitz.open(O/'Full_Review.pdf') as d:
    assert len(d)==29
    for i,(p,j) in enumerate(page_sources):
        with fitz.open(p) as src:
            assert d[i].get_text()==src[j].get_text()
            assert d[i].get_pixmap(matrix=fitz.Matrix(.6,.6)).samples==src[j].get_pixmap(matrix=fitz.Matrix(.6,.6)).samples
qa['full_review']={'pages':29,'source_text_and_pixels_match':True}
(V/'document_visual_qa.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')
audit=json.loads((V/'editorial_audit.json').read_text(encoding='utf-8'));audit['stage']='AUTHOR_REVIEW_PACKAGE_VERIFIED'
audit['qualification']='Focused editorial numerical and visual checks complete; author declarations, orbital methods metadata and target-journal requirements remain.'
(V/'editorial_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
scripts=V/'Editorial_scripts';scripts.mkdir(exist_ok=True)
for p in W.glob('*.py'):shutil.copy2(p,scripts/p.name)
(scripts/'README.md').write_text('Versioned editorial workflow; no raw statistical analysis was repeated. Run from work/association_manuscript_20261003 with the original source package and CSVs. Order: build_v1.py, refine_v2.py, render_v1.py, polish_v3.py, render_v1.py Manuscript STROBE_MR_Checklist, verify_v1.py, package_v1.py. The rendering helper uses the existing workspace portable LibreOffice because no soffice executable was present in the dependency bundle. Visual review must occur before packaging. The Korean README is maintained as an output document.\n',encoding='utf-8')
manifest={str(p.relative_to(O)).replace('\\','/'):sha(p) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='PACKAGE_SHA256.json'}
(O/'PACKAGE_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zpath=O.with_suffix('.zip')
with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
    for p in sorted(O.rglob('*')):
        if p.is_file():z.write(p,str(p.relative_to(O)).replace('\\','/'))
with zipfile.ZipFile(zpath) as z:
    assert z.testzip() is None
    for name,digest in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
receipt={'status':'PASS','zip':str(zpath),'sha256':sha(zpath),'files':len(manifest)+1,'bytes':zpath.stat().st_size,'full_review_pages':29,'manuscript_pages':13}
(zpath.parent/(zpath.stem+'_receipt.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2))
