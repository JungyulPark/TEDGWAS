from pathlib import Path
import sys, os, subprocess, re, json, hashlib
from docx import Document

sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent; R=W.parents[1]
S=R/'outputs/TED_Association_Manuscript_20261003'
O=R/'outputs/TED_Research_Contribution_20261003'
V=O/'Verification'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()

# Reuse the previous source-to-prose numeric checker, with explicit version,
# citation-count and wording updates. No statistical model is rerun.
source=(R/'work/association_manuscript_20261003/verify_v1.py').read_text(encoding='utf-8')
source=source.replace("O=R/'outputs/TED_Association_Manuscript_20261003';S=R/'outputs/TED_Manuscript_Reframed_20261002'", "O=R/'outputs/TED_Research_Contribution_20261003';S=R/'outputs/TED_Association_Manuscript_20261003'")
source=source.replace('All 32 references', 'All 33 references').replace('range(1,33)', 'range(1,34)')
source=source.replace("'They do not show that an A-subunit shed by human orbital fibroblasts, acting alone, causes TED independently of GD.'", "'they do not show that an A-subunit shed by human orbital fibroblasts, acting alone, causes TED independently of GD.'")
source=source.replace('Original_package_verification/final_numeric_audit.json', 'Previous_revision/Original_package_verification/final_numeric_audit.json')
generated=W/'numeric_checks.py';generated.write_text(source,encoding='utf-8')
subprocess.run([sys.executable,str(generated)],check=True)

t=(O/'Manuscript.md').read_text(encoding='utf-8')
old=(S/'Manuscript.md').read_text(encoding='utf-8')
checks=[]
def check(label,value):
    assert value,label
    checks.append(label)
stripcit=lambda s:re.sub(r'\[[0-9, ]+\]','[REFERENCE]',s)
for start,end in [('## Abstract\n','## Introduction\n'),('## Methods\n','## Results\n'),('## Results\n','## Discussion\n'),('### Limitations\n','### Conclusions\n'),('### Conclusions\n','## Declarations\n'),('## Declarations\n','## References\n')]:
    get=lambda doc:doc.split(start,1)[1].split(end,1)[0]
    check('Text preserved apart from references '+start.strip(),stripcit(get(t))==stripcit(get(old)))
check('Five connected introduction paragraphs',len(t.split('## Introduction\n',1)[1].split('## Methods\n',1)[0].strip().split('\n\n'))==5)
check('Prior GD and ophthalmopathy findings explicitly acknowledged','Ji and colleagues integrated blood methylation, expression and protein associations' in t)
check('Shared public resources and lack of wholly independent replication acknowledged','agreement with earlier work is not wholly independent replication' in t)
check('No formal difference inferred from differing significance','not a formal test of different effects' in t)
check('Prior probability and decision threshold comparisons are accurately separated','5×10⁻⁵ and a PP.H4 threshold above 0.50' in t and '10⁻⁵ and at least 0.80' in t)
check('Verified 2026 primary source DOI present','doi:10.1186/s12916-026-04957-y.' in t)
check('No superiority or first-study assertion',not re.search(r'\b(first study|first to|unprecedented|superior to prior)\b',t,re.I))

# Trace original citations through reordering, rather than merely checking range.
audit=json.loads((V/'editorial_audit.json').read_text(encoding='utf-8'))
mapping={int(k):v for k,v in audit['reference_old_to_new'].items()}
refparse=lambda x:{int(n):s.strip() for n,s in re.findall(r'(?ms)^(\d+)\. (.*?)(?=^\d+\. |\Z)',x.split('## References\n',1)[1])}
oldrefs=refparse(old);newrefs=refparse(t)
for n,ref in oldrefs.items():check('Preserved reference identity '+str(n),newrefs[mapping[n]]==ref)
for start,end in [('## Methods\n','## Results\n'),('## Results\n','## Discussion\n')]:
    a=old.split(start,1)[1].split(end,1)[0]
    b=t.split(start,1)[1].split(end,1)[0]
    renumbered=re.sub(r'\[([0-9, ]+)\]',lambda m:'['+', '.join(str(mapping[int(n)]) for n in re.findall(r'\d+',m[1]))+']',a)
    check('Citation identity preserved '+start.strip(),renumbered==b)
check('Unchanged STROBE checklist still states eight discussion paragraphs',sha(S/'STROBE_MR_Checklist.docx')==sha(O/'STROBE_MR_Checklist.docx'))
report={'status':'PASS','scope':'Editorial framing, preservation and citation-identity checks','checks_count':len(checks),'checks':checks}
(V/'contribution_revision_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'editorial_checks':len(checks),'status':'PASS'}),flush=True)

subprocess.run([sys.executable,str(W/'render.py')],check=True)
