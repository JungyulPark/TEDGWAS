from pathlib import Path
from PIL import Image,ImageOps,ImageDraw
from pypdf import PdfReader
import json
W=Path(__file__).resolve().parent
Q=W/'qa';Q.mkdir(exist_ok=True)
report={}
for name,folder in [('Manuscript',W/'render_final/Manuscript'),('Supplementary_Material',W/'render/Supplementary_Material'),('STROBE_MR_Checklist',W/'render/STROBE_MR_Checklist')]:
 pages=sorted(folder.glob('page-*.png'),key=lambda p:int(p.stem.split('-')[-1]))
 for off in range(0,len(pages),4):
  selected=pages[off:off+4];canvas=Image.new('RGB',(1600,2180),'#dddddd');dr=ImageDraw.Draw(canvas)
  for i,p in enumerate(selected):
   im=Image.open(p).convert('RGB');im.thumbnail((780,1035))
   x=(i%2)*800+10;y=(i//2)*1090+35
   canvas.paste(im,(x,y));dr.text((x,y-23),f'{name} {p.stem}',fill='black')
  canvas.save(Q/f'{name}_{off//4+1}.jpg',quality=92)
 reader=PdfReader(folder/(name+'.pdf'))
 report[name]=[{'page':i+1,'characters':len(p.extract_text()),'start':p.extract_text()[:90],'end':p.extract_text()[-100:]} for i,p in enumerate(reader.pages)]
(Q/'pdf_page_text.json').write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(json.dumps({n:len(v) for n,v in report.items()}))
