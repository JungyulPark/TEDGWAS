from pathlib import Path
from pypdf import PdfReader,PdfWriter
from PIL import Image,ImageChops
import json,shutil,hashlib,sys
sys.stdout.reconfigure(encoding='utf-8')
W=Path(__file__).resolve().parent;R=W.parents[1];O=R/'outputs/TED_Tissue_Integrated_20261003';A=R/'outputs/TED_Tissue_Extension_20261003';S=R/'outputs/TED_Allele_Integrated_20261003';V=O/'Verification'
for name,folder in [('Manuscript',W/'render_release/Manuscript'),('Supplementary_Material',W/'render/Supplementary_Material'),('STROBE_MR_Checklist',W/'render/STROBE_MR_Checklist')]:shutil.copy2(folder/(name+'.pdf'),O/(name+'.pdf'))
writer=PdfWriter();reader=PdfReader(O/'Supplementary_Material.pdf');writer.add_page(reader.pages[6]);writer.write(O/'TableS3_Tissue_extension.pdf')
writer=PdfWriter();parts=[];start=0
for label,p in [('Manuscript',O/'Manuscript.pdf'),('Tables 1–3',O/'Tables.pdf')]+[(f'Figure {i}',O/f'Figures/Figure{i}.pdf') for i in [1,2,3]]+[('Supplementary material',O/'Supplementary_Material.pdf')]:
 reader=PdfReader(p)
 for page in reader.pages:writer.add_page(page)
 writer.add_outline_item(label,start);parts.append({'section':label,'first_page':start+1,'pages':len(reader.pages)});start+=len(reader.pages)
writer.add_metadata({'/Title':'TSHR and IGF1R — tissue extension author review','/Subject':'2026-10-03; updated tissue MR and Table S3; main tables and figures retained'})
writer.write(O/'Full_Review.pdf');(V/'full_review_contents.json').write_text(json.dumps(parts,indent=2),encoding='utf-8')
# Copy summary-statistic source extracts and provenance, but not reference individuals' genotype files.
dest=V/'Tissue_extension';dest.mkdir(exist_ok=True)
for p in A.rglob('*'):
 if not p.is_file() or p.suffix.lower() in ['.bed','.fam']:continue
 rel=p.relative_to(A);q=dest/rel;q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,q)
shutil.copytree(S/'Verification/Allele_audit',V/'Previous_allele_audit',dirs_exist_ok=True)
scripts=V/'Tissue_extension_scripts';scripts.mkdir(exist_ok=True)
for p in W.iterdir():
 if p.suffix in ['.py','.R','.mjs'] and p.is_file():shutil.copy2(p,scripts/p.name)
shutil.copy2(W/'blocks.json',scripts/'blocks.json')
(scripts/'REPRODUCE.md').write_text('''# Execution record

Scripts 01–06 retrieve the fixed sources, select instruments, align disease records, calculate the 50-cell grid and independently verify the 20 estimates in base R. 07 creates the review forest plot. Later scripts edit, render, validate and package the manuscript.

These are the executed workspace scripts with their original relative workspace assumptions, not a portable one-command installer. To reproduce, place them under work/tissue_extension_20261003 beneath the workspace root; restore source datasets at the recorded paths or explicitly adapt paths. Sources, protocol freezes, request metadata, hashes and full selected records are in the adjacent Tissue_extension folder. Python pandas/numpy, base R 4.3.3 and PLINK 1.9 were used; document/spreadsheet integration also used python-docx, pypdf, bundled LibreOffice and Artifact Tool. Original large GWAS files, cached GTEx compressed downloads and individual reference genotypes are not redistributed here.

The first extraction attempt with PLINK's 1 GB memory allocation failed; the preserved log identifies this resource limit. An 8 GB extraction then succeeded. It did not alter the selection criteria. Subsequent clumping and independent genotype-based selection agreed. Initial document assembly assertions were corrected before final output; final source, numeric and rendering checks are recorded in final_checks.json. R_sessionInfo.txt records the statistical verification runtime.

The workbook integration reuses the OPC merge routine from the previous allele extension after Artifact Tool authored all changed cells. The earlier helper is preserved as previous_merge_workbook.py. To run 13_prepare_merge.py unchanged, restore that helper to work/allele_integration_20261003/merge_workbook.py.

Completed analytical CSVs are immutable. To reproduce into a fresh run, change output directory names rather than overwriting this reviewed run.
''',encoding='utf-8')
shutil.copy2(W.parent/'allele_integration_20261003/merge_workbook.py',scripts/'previous_merge_workbook.py')
# Identify final-render pages requiring a second visual look; all other pages are pixel-identical to the inspected render.
changes=[]
for p in (W/'render_release/Manuscript').glob('page-*.png'):
 old=W/'render_final/Manuscript'/p.name
 if ImageChops.difference(Image.open(old).convert('RGB'),Image.open(p).convert('RGB')).getbbox():changes.append(int(p.stem.split('-')[-1]))
(W/'qa/final_changed_pages.json').write_text(json.dumps(sorted(changes)))
final=json.loads((V/'final_checks.json').read_text(encoding='utf-8'));summary=json.loads((V/'integration_summary.json').read_text(encoding='utf-8'))
summary.update(abstract_words=final['abstract_words'],main_words=final['main_words'],limitations_paragraphs=2)
(V/'integration_summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
print(json.dumps({'full_review':parts,'pages_requiring_final_visual_check':sorted(changes)},indent=2))
