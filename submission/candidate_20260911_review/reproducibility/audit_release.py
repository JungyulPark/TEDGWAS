"""Verify release records against delivered bytes; run after content/visual QA.

This verifies the identity of visually reviewed documents, not their appearance.
Rebuilding a DOCX invalidates its visual record until it is inspected again.
"""
from pathlib import Path
import hashlib,json,runpy

O=Path(__file__).resolve().parents[1]
P=O/'provenance'
master=O.parents[1]/'MANUSCRIPT_TED_TRAP_v5_MASTER.md'
def read(name):return json.loads((P/name).read_text(encoding='utf-8'))
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
errors=[]
def check(ok,label):
    if not ok:errors.append(label)
manifest=read('submission_manifest.json')
check(sha(master)==manifest['manuscript_sha256'],'Master SHA256')
check(hashlib.md5(master.read_bytes()).hexdigest()==manifest['manuscript_md5'],'Master MD5')
check(manifest['numeric_audit']==read('integrity_audit.json'),'Current numerical audit summary')
check(manifest['numeric_audit']['status']=='PASS','Numerical audit passed')
figure_audit=read('figure_numeric_audit.json')
check(manifest['figure_numeric_audit']==figure_audit,'Current figure audit summary')
check(figure_audit['status']=='PASS','Figure numerical audit passed')
for name,rec in manifest['files'].items():
    p=O/name
    check(p.exists() and sha(p)==rec['sha256'] and p.stat().st_size==rec['bytes'],'Manifest file: '+name)
qa=read('document_visual_qa.json')
check(qa['status']=='PASS','Visual review status')
check(qa['total_pages']==sum(r['page_count'] for r in qa['documents']),'Visual page total')
for rec in qa['documents']:
    check(sha(O/rec['document'])==rec['sha256'],'Reviewed DOCX identity: '+rec['document'])
    check(rec['all_pages_visually_reviewed'] and rec['page_count']==len(rec['page_images']),'All pages recorded: '+rec['document'])
figures=read('figure_verification.json')
pdf_text=read('figure_pdf_text_checks.json')
layout=read('figure_layout_checks.json')
check(figures['status']=='PASS','Figure visual review status')
check(figures['source_manifest_sha256']==sha(P/'clinical_figure_sources.json'),'Reviewed figure source identity')
check(figures['main_figures']==3 and figures['supplementary_figures']==1,'Complete figure set')
check({r['name'] for r in figures['figures']}=={'Figure1','Figure2','Figure3','FigureS1'},'Exact reviewed figure names')
for rec in figures['figures']:
    for ext in ['png','pdf']:
        check(sha(O/'figures'/(rec['name']+'.'+ext))==rec[ext+'_sha256'],'Reviewed figure identity: '+rec['name']+'.'+ext)
    check(min(rec['dpi'])>=299,'Figure resolution: '+rec['name'])
    check(rec['visual_review'].startswith('PASS'),'Figure reviewed: '+rec['name'])
    check(pdf_text[rec['name']]['all_pdf_text_black'],'Black PDF text: '+rec['name'])
    check(pdf_text[rec['name']]['pdf_sha256']==rec['pdf_sha256'],'Checked PDF text identity: '+rec['name'])
for name in ['Figure1','Figure2','Figure3']:
    check(layout[name]['all_text_black'],'Black plot labels: '+name)
check(layout['Figure3']['cell_labels_with_padding']==27,'All Figure3 matrix labels fit inside cells')
check(layout['Figure3']['strong_support_borders']==6 and layout['Figure3']['border_matches_cell_geometry'],'Six aligned strong-support borders')
wc=runpy.run_path(str(master.parent/'scripts/26_wordcount_main_text.py'))
sections=wc['split_sections'](master.read_text(encoding='utf-8'))
counts=read('word_counts.json')
check(counts['abstract_words']==wc['count_words'](sections['Abstract']),'Current abstract count')
check(counts['main_words_excluding_headings']==sum(wc['count_words'](sections[k]) for k in wc['MAIN']),'Current main count')
print(json.dumps({'status':'FAIL' if errors else 'PASS','manifest_files':len(manifest['files']),'reviewed_documents':len(qa['documents']),'reviewed_pages':qa['total_pages'],'reviewed_figures':len(figures['figures']),'errors':errors},indent=2))
raise SystemExit(bool(errors))
