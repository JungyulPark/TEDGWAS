from pathlib import Path
import sys,re,ast,json,difflib
from docx import Document
from docx.shared import RGBColor
from docx.oxml.ns import qn
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1]
S=R/'outputs/TED_Manuscript_Reframed_20261002';O=R/'outputs/TED_Association_Manuscript_20261003'
t=(O/'Manuscript.md').read_text(encoding='utf-8')
t=t.replace('Colocalization support for *TSHR* varied by dataset and prior; *IGF1R* never reached 0.80.', 'Secondary colocalization gave *TSHR* shared-variant probabilities of 0.951 in BBJ and 0.264 in FinnGen GD, with prior sensitivity; *IGF1R* never reached 0.80.')
# Define the scale in the abstract without changing the estimates.
a,b=t.split('**Keywords:**',1)
a=a.replace('odds ratio [OR] 0.12','odds ratio [OR] per approximately one SD 0.12')
t=a+'**Keywords:**'+b
(O/'Manuscript.md').write_text(t,encoding='utf-8')
tree=ast.parse((W/'build_v1.py').read_text(encoding='utf-8'))
defs=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ['inline','rebuild']],type_ignores=[])
exec(compile(defs,str(W/'build_v1.py'),'exec'))
rebuild('Manuscript',t)
audit=json.loads((O/'Verification/editorial_audit.json').read_text(encoding='utf-8'))
abstract=t.split('## Abstract\n',1)[1].split('**Keywords:**',1)[0]
audit['abstract_words_including_labels']=len(abstract.replace('*','').split())
assert audit['abstract_words_including_labels']<=250
(O/'Verification/editorial_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
old=(S/'Manuscript.md').read_text(encoding='utf-8')
(O/'Verification/manuscript_changes.diff').write_text(''.join(difflib.unified_diff(old.splitlines(True),t.splitlines(True),fromfile='20261002/Manuscript.md',tofile='20261003/Manuscript.md')),encoding='utf-8')
print('Abstract:',audit['abstract_words_including_labels'])
