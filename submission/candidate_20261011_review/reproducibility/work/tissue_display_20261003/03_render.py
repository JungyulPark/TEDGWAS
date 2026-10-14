from pathlib import Path
import runpy,sys,os
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Tables_Figures_20261003'
os.environ['PATH']=str(R/'work/tools/lo/LibreOffice/program')+os.pathsep+'${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/Library/bin'+os.pathsep+os.environ['PATH']
renderer='${LOCAL_USER_HOME}/.codex/plugins/cache/openai-primary-runtime/documents/26.905.11957/skills/documents/render_docx.py'
sys.argv=[renderer,str(O/'Tissue_result_tables.docx'),'--output_dir',str(W/'render_release'),'--emit_pdf'];runpy.run_path(renderer,run_name='__main__')
