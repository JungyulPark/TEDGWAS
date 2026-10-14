from pathlib import Path
import os,sys,shutil,runpy
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Allele_Integrated_20261003'
lo=R/'work/tools/lo/LibreOffice/program'
os.environ['PATH']=str(lo)+os.pathsep+'${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin'+os.pathsep+os.environ['PATH']
assert Path(shutil.which('soffice.exe')).resolve()==(lo/'soffice.exe').resolve()
renderer=Path('${LOCAL_USER_HOME}/.codex/plugins/cache/openai-primary-runtime/documents/26.905.11957/skills/documents/render_docx.py')
for name in ['Manuscript','Supplementary_Material','STROBE_MR_Checklist']:
    folder=W/'render'/name
    assert not list(folder.glob('page-*.png')),'Avoid stale pages'
    sys.argv=[str(renderer),str(O/(name+'.docx')),'--output_dir',str(folder),'--emit_pdf']
    runpy.run_path(str(renderer),run_name='__main__')
    print(name,len(list(folder.glob('page-*.png'))),'pages',flush=True)
