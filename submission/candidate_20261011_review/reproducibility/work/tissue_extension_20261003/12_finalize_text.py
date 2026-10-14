from pathlib import Path
import sys,os,runpy,re
from docx import Document
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Integrated_20261003'
sys.stdout.reconfigure(encoding='utf-8')
md=O/'Manuscript.md'
t=md.read_text(encoding='utf-8')
assert t.count('Adipose *TSHR* associations passed correction')==1
t=t.replace('Adipose *TSHR* associations passed correction','Adipose *TSHR*–GD associations passed correction')
md.write_text(t,encoding='utf-8')
doc=Document(O/'Manuscript.docx')
p=next(p for p in doc.paragraphs if 'Adipose TSHR associations passed correction' in p.text)
assert any(' associations passed correction' in r.text for r in p.runs)
for r in p.runs:
 if ' associations passed correction' in r.text:r.text=r.text.replace(' associations passed correction','–GD associations passed correction')
doc.save(O/'Manuscript.docx')
os.environ['PATH']=str(R/'work/tools/lo/LibreOffice/program')+os.pathsep+'${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin'+os.pathsep+os.environ['PATH']
renderer='${LOCAL_USER_HOME}/.codex/plugins/cache/openai-primary-runtime/documents/26.905.11957/skills/documents/render_docx.py'
folder=W/'render_final/Manuscript';assert not folder.exists()
sys.argv=[renderer,str(O/'Manuscript.docx'),'--output_dir',str(folder),'--emit_pdf'];runpy.run_path(renderer,run_name='__main__')
print('Final manuscript rendered',flush=True)
