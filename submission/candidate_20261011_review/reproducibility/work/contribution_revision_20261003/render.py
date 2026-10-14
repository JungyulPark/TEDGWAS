from pathlib import Path
import sys, os, runpy, shutil
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1]
O=R/'outputs/TED_Research_Contribution_20261003'
# The bundled launcher resets PATH when starting another Python process.
# Set PATH after launch and invoke the unmodified skill renderer in-process.
# Windows has no bundled LibreOffice, so reuse the task-local portable copy.
lo=R/'work/tools/lo/LibreOffice/program'
os.environ['PATH']=str(lo)+os.pathsep+'${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin'+os.pathsep+os.environ['PATH']
assert Path(shutil.which('soffice.exe')).resolve()==(lo/'soffice.exe').resolve()
renderer=Path('${LOCAL_USER_HOME}/.codex/plugins/cache/openai-primary-runtime/documents/26.905.11957/skills/documents/render_docx.py')
for name in ['Manuscript','Cover_Letter_Draft']:
    folder=W/'render'/name
    assert not list(folder.glob('page-*.png')), 'Avoid stale rendered pages.'
    sys.argv=[str(renderer),str(O/(name+'.docx')),'--output_dir',str(folder),'--emit_pdf']
    runpy.run_path(str(renderer),run_name='__main__')
    print(name,len(list(folder.glob('page-*.png'))),'pages',flush=True)
