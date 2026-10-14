from pathlib import Path
import sys,json,hashlib,shutil,zipfile,re
W=Path(__file__).resolve().parent; R=W.parents[1]
sys.path.insert(0,str(R/'work/pydeps'))
import fitz
sys.stdout.reconfigure(encoding='utf-8')
O=R/'outputs/TED_Manuscript_Reframed_20261002'; V=O/'Verification'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
names=['Manuscript','Tables','Supplementary_Material','Cover_Letter_Draft','STROBE_MR_Checklist']
expected=[13,5,8,1,1]
comparison=json.loads((W/'render_comparison.json').read_text())
qa={'status':'PASS','scope':'Layout review, not scientific-validity or author-approval certification','documents':{}}
for name,n in zip(names,expected):
 source=W/'render_release'/name/(name+'.pdf')
 with fitz.open(source) as doc:assert len(doc)==n
 shutil.copy2(source,O/(name+'.pdf'))
 qa['documents'][name]={'pages':n,'docx_sha256':sha(O/(name+'.docx')),'pdf_sha256':sha(source),'all_pages_visually_reviewed':True,'release_comparison':comparison[name],'page_images_sha256':{p.name:sha(p) for p in (W/'render_release'/name).glob('page-*.png')}}

merge_checks=[]
def merge(output,parts):
 merged=fitz.open();toc=[];source_pages=[]
 for title,path,pages in parts:
  with fitz.open(path) as doc:
   if pages is None:pages=list(range(len(doc)))
   toc.append([1,title,len(merged)+1])
   for i in pages:
    merged.insert_pdf(doc,from_page=i,to_page=i)
    source_pages.append((path,i))
 merged.set_toc(toc);merged.save(output,garbage=4,deflate=True);merged.close()
 with fitz.open(output) as result:
  assert len(result)==len(source_pages)
  for i,(path,j) in enumerate(source_pages):
   with fitz.open(path) as src:
    a=result[i];b=src[j]
    assert a.get_text()==b.get_text(),(output,i,'text')
    assert a.rect==b.rect,(output,i,'size')
    assert a.get_pixmap(matrix=fitz.Matrix(0.75,0.75)).samples==b.get_pixmap(matrix=fitz.Matrix(0.75,0.75)).samples,(output,i,'pixels')
 merge_checks.append({'file':str(output.relative_to(O)),'pages':len(source_pages),'text_and_render_equal_to_reviewed_sources':True})

for title,pages in [('Table1',[0]),('Table2',[1]),('Table3',[2,3,4])]:
 merge(O/'Tables'/(title+'.pdf'),[(title,O/'Tables.pdf',pages)])
figureparts=[(n,O/'Figures'/(n+'.pdf'),None) for n in ['Figure1','Figure2','Figure3','FigureS1','FigureS2']]
merge(O/'All_Figures.pdf',figureparts)
merge(O/'Full_Review.pdf', [('Manuscript',O/'Manuscript.pdf',None),('Main tables',O/'Tables.pdf',None)]+figureparts[:3]+[('Supplementary material',O/'Supplementary_Material.pdf',None)])
qa['merged_pdfs']=merge_checks
qa['figures']={n:{'visually_reviewed':True,'pdf_sha256':sha(p),'png_sha256':sha(p.with_suffix('.png'))} for n,p,_ in figureparts}
(V/'document_visual_qa.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf-8')

md=O/'Tables.md';txt=md.read_text(encoding='utf-8')
if 'In the formatted DOCX/PDF' not in txt:
 txt+='\nIn the formatted DOCX/PDF, bold P values denote nominal P < 0.05; screen and focused thresholds are stated separately.\n'
 md.write_text(txt,encoding='utf-8')

follow=R/'work/followup_validation_20261002'
for rel in ['plan_freeze.json','final_checks_v1/verification.json','alias_audit_v1/status.json','TSHR_v2/R_verification.csv','candidates_v1/R_verification.csv','instruments_v1/R_verification.csv','inputs/outcomes_manifest.json','inputs/eqtl_manifest.json']:
 p=follow/rel;assert p.exists(),p
 d=V/'Analysis_records'/rel;d.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,d)
scripts=V/'Scripts';scripts.mkdir(exist_ok=True)
for p in W.glob('*.py'):shutil.copy2(p,scripts/p.name)
(scripts/'README.md').write_text('These scripts retain workspace-relative analysis-source paths. Run them from their original work/manuscript_reframe_20261002 location with the original source files available. Build order: build_layout_v2.py (executes build_v1.py), editorial_v3.py, polish_v4.py, render_release_v3.py, compare_render_v1.py, package_v1.py. The pre-existing supplementary workbook is required by the layout rebuild; build_v1.py creates it during a full initial build. verify_v1.py records the source-to-workbook and Table 2 checks. This archive is an author-review package, not a self-contained raw-GWAS distribution.\n',encoding='utf-8')
audit=json.loads((V/'build_audit.json').read_text(encoding='utf-8'))
for path,digest in audit['source_hashes'].items():assert sha(R/path)==digest
audit['status']='PASS_NUMERICAL_STRUCTURE_AND_VISUAL_AUTHOR_REVIEW_PACKAGE'
audit['release_note']='Author declarations, orbital collection/sequencing metadata and journal-specific formatting remain to be confirmed.'
(V/'build_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
manuscript=(O/'Manuscript.md').read_text(encoding='utf-8');before=manuscript.split('## References')[0];used=[]
for group in re.findall(r'\[([0-9, ]+)\]',before):
 for n in map(int,re.findall(r'\d+',group)):
  if n not in used:used.append(n)
assert used==list(range(1,29))
readme=O/'READ_FIRST_KO.md';t=readme.read_text(encoding='utf-8')
if '## 검토 파일과 완료된 확인' not in t:
 t+='\n## 검토 파일과 완료된 확인\n\n- 먼저 읽기: Full_Review.pdf — 본문·주표·주그림·보충자료 29페이지 통합본.\n- 수정 가능 원고: Manuscript.docx. 초록 241단어, 본문 3,496단어, 참고문헌 28개.\n- 표: Tables.docx / Tables.pdf. Table 3은 A–C 세 패널로 묶었습니다.\n- 그림: All_Figures.pdf 및 Figures 폴더의 개별 PDF/PNG.\n- 보충자료: Supplementary_Material.docx / .pdf, Supplementary_Dataset.xlsx.\n- 투고문: Cover_Letter_Draft.docx — 저널 및 저자 선언 확인 전 초안.\n- 수치 확인: 원본과 상세 Excel의 16개 시트, 1,391,321개 셀을 비교했습니다. Table 2의 15개 비교와 참고문헌 등장 순서를 검사했습니다.\n- 문서 28페이지의 렌더링과 주그림·보충그림을 확인했습니다. 통합 PDF는 검토한 원본 PDF와 페이지별 텍스트·렌더링이 일치합니다.\n- Verification 폴더는 저자 검토용 기록입니다. 실제 저널 업로드 파일은 저널 요건에 맞춰 별도로 선택해야 합니다. 수치·레이아웃 검사는 인과성 또는 게재 가능성의 보증이 아닙니다.\n'
 readme.write_text(t,encoding='utf-8')
manifest={str(p.relative_to(O)).replace('\\','/'):sha(p) for p in sorted(O.rglob('*')) if p.is_file() and p.name!='PACKAGE_SHA256.json'}
(O/'PACKAGE_SHA256.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
zipout=O.with_suffix('.zip')
with zipfile.ZipFile(zipout,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(O.rglob('*')):
  if p.is_file():z.write(p,str(p.relative_to(O)).replace('\\','/'))
with zipfile.ZipFile(zipout) as z:
 assert z.testzip() is None
 for name,digest in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==digest
receipt={'zip':zipout.name,'sha256':sha(zipout),'bytes':zipout.stat().st_size,'status':'ARCHIVE_HASH_AND_CRC_PASS','files':len(manifest)+1}
(zipout.parent/(zipout.stem+'_receipt.json')).write_text(json.dumps(receipt,indent=2),encoding='utf-8')
print(json.dumps(receipt,indent=2),flush=True)
