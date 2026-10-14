from pathlib import Path
from pypdf import PdfReader
from PIL import Image,ImageDraw
import json,math,subprocess,sys
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];S=R/'outputs/TED_Tissue_Integrated_20261003';O=R/'outputs/TED_Final_Review_20261011'
expected=[];index=[]
def expect(label,p,which=None):
 rd=PdfReader(p);ids=list(range(len(rd.pages))) if which is None else which;start=len(expected)
 expected.extend((label,rd.pages[i]) for i in ids)
 index.append({'item':label,'first_page':start+1,'pages':len(ids)})
for i in [1,2,3]:expect(f'Table {i}',S/f'Tables/Table{i}.pdf')
for i in [1,2,3]:expect(f'Figure {i}',S/f'Figures/Figure{i}.pdf')
for i,pg in [(1,4),(2,5),(3,6)]:expect(f'Table S{i}',S/'Supplementary_Material.pdf',[pg])
for i in [1,2]:expect(f'Figure S{i}',S/f'Figures/FigureS{i}.pdf')
rd=PdfReader(O/'Final_Tables_Figures.pdf');assert len(rd.pages)==len(expected)==13
maxdiff=0
for a,(label,b) in zip(rd.pages,expected):
 assert a.extract_text()==b.extract_text(),label
 assert all(math.isclose(float(x),float(y),rel_tol=0,abs_tol=1e-5) for x,y in zip(a.mediabox,b.mediabox)),label
 maxdiff=max(maxdiff,max(abs(float(x)-float(y)) for x,y in zip(a.mediabox,b.mediabox)))
(O/'Verification/Tables_figures_contents.json').write_text(json.dumps(index,indent=2),encoding='utf-8')
(O/'Verification/PDF_assembly_checks.json').write_text(json.dumps({'status':'PASS','pages':13,'all_page_text_unchanged':True,'maximum_page_size_serialization_difference_points':maxdiff,'note':'PDF coordinate serialization rounds sub-micropoint precision; no layout change.'},indent=2),encoding='utf-8')
render=W/'render';render.mkdir(exist_ok=True)
exe='${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin/pdftoppm.exe'
subprocess.run([exe,'-png','-r','100',str(O/'Final_Tables_Figures.pdf'),str(render/'page')],check=True)
pages=sorted(render.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[-1]))
assert len(pages)==13
for k in range(0,len(pages),4):
 group=pages[k:k+4];can=Image.new('RGB',(1800,1450),'#dddddd');dr=ImageDraw.Draw(can)
 for j,p in enumerate(group):
  im=Image.open(p).convert('RGB');im.thumbnail((870,675))
  x=(j%2)*900+15;y=(j//2)*725+30
  can.paste(im,(x,y));dr.text((x,y-20),p.stem,fill='black')
 can.save(render/f'contact-{k//4+1}.jpg',quality=92)
print(json.dumps({'pages':len(pages),'max_size_difference_points':maxdiff,'index':index},indent=2))
