from pathlib import Path
from pypdf import PdfReader,PdfWriter
import shutil,json,hashlib,sys,math
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];S=R/'outputs/TED_Tissue_Integrated_20261003';T=R/'outputs/TED_Tissue_Tables_Figures_20261003';O=R/'outputs/TED_Final_Review_20261011'
assert not O.exists(),'Keep previous deliveries immutable'
O.mkdir()
for p in S.rglob('*'):
 if not p.is_file():continue
 rel=p.relative_to(S)
 if rel.as_posix() in ['README_KO.md','MANIFEST_SHA256.json']:continue
 q=O/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
shutil.copy2(S/'MANIFEST_SHA256.json',O/'Verification/Integrated_source_manifest_20261003.json')
shutil.copy2(T/'MANIFEST_SHA256.json',O/'Verification/Tissue_display_source_manifest_20261003.json')
shutil.copy2(W/'source_verification.json',O/'Verification/Delivery_source_verification_20261011.json')
review=O/'Review_only';review.mkdir()
for n in ['Tissue_result_tables.docx','Tissue_result_tables.pdf','Tissue_results_review.pdf','PLACEMENT_KO.md']:
 shutil.copy2(T/n,review/n)
for ext in ['pdf','png']:shutil.copy2(T/'Figures'/('Tissue_primary_summary.'+ext),O/'Additional_analysis_figure'/('Tissue_primary_summary.'+ext))
# Assemble publication tables and figures without changing the content of any page.
writer=PdfWriter();index=[];expected=[]
def append(label,p,which=None):
 rd=PdfReader(p);ids=list(range(len(rd.pages))) if which is None else which
 start=len(writer.pages)
 for i in ids:
  pg=rd.pages[i];writer.add_page(pg)
  expected.append({'label':label,'source':str(p.relative_to(S)),'source_page':i+1,'text':pg.extract_text(),'size':[float(x) for x in pg.mediabox]})
 writer.add_outline_item(label,start);index.append({'item':label,'first_page':start+1,'pages':len(ids)})
for i in [1,2,3]:append(f'Table {i}',S/f'Tables/Table{i}.pdf')
for i in [1,2,3]:append(f'Figure {i}',S/f'Figures/Figure{i}.pdf')
for i,pg in [(1,4),(2,5),(3,6)]:append(f'Table S{i}',S/'Supplementary_Material.pdf',[pg])
for i in [1,2]:append(f'Figure S{i}',S/f'Figures/FigureS{i}.pdf')
writer.add_metadata({'/Title':'Final tables and figures for TSHR and IGF1R study','/Subject':'Scientific content finalized 2026-10-03; delivery verification 2026-10-11. Main Tables 1-3, Figures 1-3, Tables S1-S3, Figures S1-S2.'})
fp=O/'Final_Tables_Figures.pdf';writer.write(fp)
reader=PdfReader(fp);assert len(reader.pages)==13
for a,b in zip(reader.pages,expected):
 assert a.extract_text()==b['text'],b['label']
 assert all(math.isclose(float(x),y,rel_tol=0,abs_tol=1e-5) for x,y in zip(a.mediabox,b['size'])),b['label']
(O/'Verification/Tables_figures_contents.json').write_text(json.dumps(index,indent=2),encoding='utf-8')
print(json.dumps({'package':str(O),'combined_pdf_pages':len(reader.pages),'index':index},indent=2))
