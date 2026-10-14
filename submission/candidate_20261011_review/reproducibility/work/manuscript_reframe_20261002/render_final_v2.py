from pathlib import Path
p=Path(__file__).with_name('render_v1.py')
src=p.read_text(encoding='utf-8').replace("QA=W/'render'","QA=W/'render_final'")
exec(compile(src,str(p),'exec'))
