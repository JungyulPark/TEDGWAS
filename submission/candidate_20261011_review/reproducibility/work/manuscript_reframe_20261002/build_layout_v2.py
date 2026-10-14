"""Layout-only revision after visual inspection of version 1; frozen outputs retained."""
from pathlib import Path
W=Path(__file__).resolve().parent
src=(W/'build_v1.py').read_text(encoding='utf-8')
src=src.replace('from content_v1 import *','from content_v1 import *\nRESULTS = RESULTS.replace("approximately 0.978 in both ophthalmopathy endpoints", "0.978 and 0.979 in the two ophthalmopathy endpoints")')
src=src.replace('2,544 instrumented genes; 6,135 records','2,544 genes; 6,135 SNP–gene records')
src=src.replace('Two overlapping ophthalmopathy endpoints','Two overlapping GO definitions')
src=src.replace("ax.text(total*.68,threshold+.45", "ax.text(total*.02,threshold+.7")
src=src.replace('top=.94,bottom=.07','top=.90,bottom=.07')
src=src.replace("ax.set_xlabel('OR after omission (95% CI)');", "ax.set_xticks([.5,1,2,4] if g=='IGF1R' else [.1,1,10], ['0.5','1','2','4'] if g=='IGF1R' else ['0.1','1','10']);ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter());ax.set_xlabel('OR after omission (95% CI)');")
start=src.index('# A single consolidated results workbook')
end=src.index('# Plain academic DOCX')
src=src[:start]+"# Existing verified source-value workbook is unchanged in layout revision.\n"+src[end:]
exec(compile(src,str(W/'build_v1.py'),'exec'))
