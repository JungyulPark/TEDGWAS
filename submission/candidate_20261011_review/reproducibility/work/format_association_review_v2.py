from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
p=root/'outputs/TED_Association_Review_20261003/Association_Review_KO.md'
t=p.read_text(encoding='utf-8')
t=re.sub(r'(?m)(^\|[^\n]*)\n\n(?=\|)',r'\1\n',t)
p.write_text(t,encoding='utf-8')
assert t.count('|---|')>=2
print('Review tables formatted; numerical results unchanged.')
