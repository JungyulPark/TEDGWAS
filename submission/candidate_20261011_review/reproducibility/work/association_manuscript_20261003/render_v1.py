from pathlib import Path
import sys,os,subprocess
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Association_Manuscript_20261003'
# Workspace dependency loader provides no soffice; use the existing portable
# workspace copy already used to render the source package, not desktop Office.
os.environ['PATH']=str(R/'work/tools/lo/LibreOffice/program')+os.pathsep+'${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin'+os.pathsep+os.environ['PATH']
renderer=Path('${LOCAL_USER_HOME}/.codex/plugins/cache/openai-primary-runtime/documents/26.905.11957/skills/documents/render_docx.py')
for name in sys.argv[1:] or ['Manuscript','Cover_Letter_Draft','STROBE_MR_Checklist']:
    p=subprocess.run([sys.executable,str(renderer),str(O/(name+'.docx')),'--output_dir',str(W/'render'/name),'--emit_pdf'],capture_output=True,text=True,encoding='utf-8',errors='replace')
    if p.returncode:print(p.stdout,p.stderr,flush=True);raise SystemExit(p.returncode)
    print(name,len(list((W/'render'/name).glob('page-*.png'))),'pages',flush=True)
