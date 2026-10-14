from pathlib import Path
import csv, json, math, re, shutil, hashlib
from docx import Document
from docx.shared import Inches, Pt, RGBColor
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT

B = Path(__file__).resolve().parent
ROOT = B.parent.parent
SOURCE = ROOT / 'outputs/TED_Paper_Tables_Figures_20260930'
OUT = ROOT / 'outputs/TED_Paper_Tables_Figures_20260930_v2'
DATA = SOURCE / 'Detailed_Data'

def read_csv(name):
    with (DATA / name).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

summary = read_csv('Stage2_integrated_summary.csv')
genes = [r['gene'] for r in summary]
assert len(genes) == len(set(genes)) == 14 and genes[-1] == 'IGF1R'
full = {r['gene_symbol']: r for r in read_csv('Stage1_ALL_2544_genes.csv')}
primary = {(r['gene'], r['outcome']): r for r in read_csv('Stage2_14_genes_primary.csv')}
settings = read_csv('Stage2_all_504_colocalization_settings.csv')
low = {(r['gene'], r['outcome']): r for r in settings
       if r['scenario'] == 'original_reference' and math.isclose(float(r['p12']), 1e-6, rel_tol=1e-12)}
af = {(r['gene'], r['outcome']): r for r in settings
      if r['scenario'] == 'available_eqtlgen' and math.isclose(float(r['p12']), 1e-5, rel_tol=1e-12)}
assert len(primary) == len(low) == len(af) == 42
outcomes = [('BBJ_Graves', 'A  BioBank Japan Graves disease'),
            ('UKB_hyperthyroid', 'B  UK Biobank hyperthyroidism'),
            ('FinnGen_GO', 'C  FinnGen Graves ophthalmopathy')]

# Build a fresh delivery version without the superseded stand-alone Table 4.
OUT.mkdir(parents=True, exist_ok=True)
for sub in ['Tables', 'Figures', 'Detailed_Data']:
    (OUT / sub).mkdir(exist_ok=True)
for p in SOURCE.rglob('*'):
    if not p.is_file():
        continue
    rel = p.relative_to(SOURCE)
    if p.name in {'Table3.docx', 'Table3.pdf', 'Table4.docx', 'Table4.pdf',
                  'All_Tables.pdf', 'PACKAGE_SHA256.json', 'READ_ME_KO.md',
                  'Figure_Legends.docx', 'Figure_Legends.pdf', 'Figure_Legends.md'}:
        continue
    dest = OUT / rel
    dest.parent.mkdir(exist_ok=True, parents=True)
    shutil.copy2(p, dest)

SUP = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
def pformat(x):
    if x < .001:
        mant, exp = f'{x:.2e}'.split('e')
        return mant + '×10' + str(int(exp)).translate(SUP)
    return f'{x:.3g}'

def ppformat(x):
    return '<0.001' if x < .001 else f'{x:.3f}'

doc = Document()
sec = doc.sections[0]
# Preserve the A4 portrait page size and typography of the source tables.
sec.page_width, sec.page_height = Inches(8.27), Inches(11.69)
sec.left_margin = sec.right_margin = Inches(.65)
sec.top_margin = sec.bottom_margin = Inches(.62)
sec.footer_distance = Inches(.3)
for name in ['Normal', 'Title', 'Heading 1', 'Heading 2']:
    s = doc.styles[name]
    s.font.name = 'Times New Roman'
    s.font.color.rgb = RGBColor(0, 0, 0)
    for node in list(s._element.iter(qn('w:pBdr'))):
        node.getparent().remove(node)
doc.styles['Normal'].font.size = Pt(10.5)
doc.styles['Normal'].paragraph_format.space_after = Pt(5)
doc.styles['Title'].font.size = Pt(14)
doc.styles['Title'].font.bold = True
doc.styles['Title'].paragraph_format.space_after = Pt(7)
doc.styles['Heading 2'].font.size = Pt(11)
doc.styles['Heading 2'].font.bold = True
doc.styles['Heading 2'].paragraph_format.space_before = Pt(0)
doc.styles['Heading 2'].paragraph_format.space_after = Pt(7)
footer = sec.footer.paragraphs[0]
footer.alignment = WD_ALIGN_PARAGRAPH.RIGHT
field = OxmlElement('w:fldSimple')
field.set(qn('w:instr'), 'PAGE')
footer._p.append(field)

gene_pattern = '|'.join(re.escape(g) for g in sorted(genes, key=len, reverse=True))
italics = re.compile(r'\b(?:' + gene_pattern + r'|P)\b')
def rich(p, text, size=9.5, bold=False):
    pos = 0
    for m in italics.finditer(text):
        for part, italic in [(text[pos:m.start()], False), (m[0], True)]:
            if part:
                r = p.add_run(part); r.italic = italic; r.bold = bold; r.font.size = Pt(size)
        pos = m.end()
    if text[pos:]:
        r = p.add_run(text[pos:]); r.bold = bold; r.font.size = Pt(size)

def note(text):
    p = doc.add_paragraph()
    p.paragraph_format.space_before = Pt(0)
    p.paragraph_format.space_after = Pt(5)
    p.paragraph_format.line_spacing = 1.03
    rich(p, text)

def border(cell, side, style='single', size=4, color='B8B8B8'):
    pr = cell._tc.get_or_add_tcPr()
    borders = pr.find(qn('w:tcBorders'))
    if borders is None:
        borders = OxmlElement('w:tcBorders'); pr.append(borders)
    node = borders.find(qn('w:' + side))
    if node is None:
        node = OxmlElement('w:' + side); borders.append(node)
    for k, v in [('val', style), ('sz', str(size)), ('color', color)]:
        node.set(qn('w:' + k), v)

display = []
integrated = []
widths = [.95, .38, 1.67, 1.06, .65, .65, .69, .92]
assert abs(sum(widths) - 6.97) < 1e-12
for panel, (oc, label) in enumerate(outcomes):
    if panel:
        doc.add_page_break()
    title = 'Table 3 Integrated MR and colocalization results across 14 genes'
    if panel:
        title += ' continued'
    doc.add_paragraph(title, style='Title')
    doc.add_paragraph(label, style='Heading 2')
    t = doc.add_table(rows=2, cols=8)
    t.autofit = False
    tw = t._tbl.tblPr.find(qn('w:tblW'))
    tw.set(qn('w:type'), 'dxa'); tw.set(qn('w:w'), str(round(sum(widths) * 1440)))
    for c, w in zip(t.columns, widths):
        c.width = Inches(w)
    headers = ['Gene', 'IVs', 'OR (95% CI)', 'P value', 'PP.H2', 'PP.H3', 'PP.H4', 'PP.H4\nlower prior']
    for cell, text in zip(t.rows[1].cells, headers):
        cell.text = text
    t.cell(0, 0).merge(t.cell(1, 0)).text = 'Gene'
    t.cell(0, 1).merge(t.cell(0, 3)).text = 'Mendelian randomization'
    t.cell(0, 4).merge(t.cell(0, 7)).text = 'Colocalization posterior probabilities'
    for r in t.rows:
        repeat = OxmlElement('w:tblHeader'); r._tr.get_or_add_trPr().append(repeat)
    for gene in genes:
        mr, pp = full[gene], primary[(gene, oc)]
        lo, cohort = low[(gene, oc)], af[(gene, oc)]
        p_raw = mr['pvalue__' + oc]
        digits = 3 if gene == 'MAPKAPK5' else 2
        if p_raw:
            effect = float(mr['or__' + oc])
            lci, uci = float(mr['or_ci_lower__' + oc]), float(mr['or_ci_upper__' + oc])
            effect_text = f'{effect:.{digits}f} ({lci:.{digits}f}–{uci:.{digits}f})'
            n_iv = str(int(float(mr['n_iv__' + oc])))
            p_text = pformat(float(p_raw))
        else:
            n_iv, effect_text, p_text = 'NA', 'Not estimable', 'NA'
        vals = [gene + ('†' if gene == 'IGF1R' else ''), n_iv, effect_text, p_text,
                ppformat(float(pp['PP.H2'])), ppformat(float(pp['PP.H3'])),
                ppformat(float(pp['PP.H4'])), ppformat(float(lo['PP.H4']))]
        row = t.add_row()
        for cell, text in zip(row.cells, vals):
            cell.text = text
        display.append({'gene': gene, 'outcome': oc, 'cells': vals})
        record = {'gene': gene, 'outcome': oc, 'BBJ_discovery': gene != 'IGF1R',
                  'MHC': gene in ['HLA-A', 'HLA-DQA2', 'C4A', 'TUBB', 'PSMB8'],
                  'mr_method': ('Wald ratio' if n_iv == '1' else 'Multiplicative random-effects IVW') if p_raw else 'Not estimable'}
        for key in ['n_iv', 'beta', 'se', 'or', 'or_ci_lower', 'or_ci_upper', 'pvalue', 'status']:
            record[key] = mr[key + '__' + oc]
        record['n_overlap'] = pp['n_overlap']
        record['top_snp'] = pp['top_snp']
        for h in range(5):
            record['PP.H' + str(h)] = pp['PP.H' + str(h)]
        record['H4_lower_prior_1e-6'] = lo['PP.H4']
        record['H4_cohort_AF_1e-5'] = cohort['PP.H4']
        integrated.append(record)
    seen = set()
    for ri, row in enumerate(t.rows):
        nosplit = OxmlElement('w:cantSplit'); row._tr.get_or_add_trPr().append(nosplit)
        for ci, cell in enumerate(row.cells):
            if cell._tc in seen:
                continue
            seen.add(cell._tc)
            if ri >= 1 and not (ri == 1 and ci == 0):
                cell.width = Inches(widths[ci])
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            pr = cell._tc.get_or_add_tcPr()
            margins = OxmlElement('w:tcMar')
            for side, value in [('top', 85), ('bottom', 85), ('left', 45), ('right', 45)]:
                n = OxmlElement('w:' + side); n.set(qn('w:w'), str(value)); n.set(qn('w:type'), 'dxa'); margins.append(n)
            pr.append(margins)
            for side in ['top', 'bottom', 'left', 'right']:
                border(cell, side, 'nil')
            if ri == 0:
                border(cell, 'top', size=8, color='000000')
            border(cell, 'bottom', size=7 if ri == 1 or ri == len(t.rows)-1 else 3,
                   color='000000' if ri == 1 or ri == len(t.rows)-1 else 'D0D0D0')
            text = cell.text
            for p in cell.paragraphs:
                p.clear()
                p.paragraph_format.space_after = Pt(0)
                p.paragraph_format.space_before = Pt(0)
                p.paragraph_format.line_spacing = 1
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if ci == 0 else WD_ALIGN_PARAGRAPH.CENTER
                bold = ri < 2
                if ri >= 2:
                    gene = genes[ri-2]
                    if ci == 3 and full[gene]['pvalue__'+oc]:
                        bold = float(full[gene]['pvalue__'+oc]) < (.05/2544 if panel == 0 else .05)
                    if ci == 6:
                        bold = float(primary[(gene, oc)]['PP.H4']) >= .8
                    if ci == 7:
                        bold = float(low[(gene, oc)]['PP.H4']) >= .8
                rich(p, text, size=10 if ri < 2 else 10.5, bold=bold)
    spacer = doc.add_paragraph(); spacer.paragraph_format.space_after = Pt(0); spacer.paragraph_format.space_before = Pt(0)
    spacer.paragraph_format.line_spacing = Pt(3); spacer.add_run().font.size = Pt(3)
    threshold = ('Bold P values meet the BBJ discovery threshold P < 0.05/2,544 = 1.965×10⁻⁵.' if panel == 0
                 else 'Bold P values indicate nominal P < 0.05; the BBJ discovery threshold was not applied to this dataset.')
    note('The same 13 BBJ discovery genes and IGF1R† are shown in each panel. †Biologically selected comparator; not a BBJ discovery hit. Gene order follows BBJ discovery P values, with IGF1R last. ' + threshold)
    note('ORs are per approximately one SD higher reconstructed, genetically proxied blood expression. IVs, instrumental variants; CI, confidence interval. One IV: Wald ratio; multiple IVs: multiplicative random-effects inverse-variance weighting with under-dispersion correction. Two-sided P values; original reference-frequency estimates.')
    note('PP, posterior probability; H, hypothesis. Dataset 1 is disease and dataset 2 is blood expression. H2: expression only; H3: distinct causal variants; H4: a shared causal variant. H0 (neither trait) and H1 (disease only) account for the undisplayed probability. Values <0.001 are not zero.')
    note('Primary probabilities use p₁ = p₂ = 10⁻⁴ and p₁₂ = 10⁻⁵. The last column uses the lower shared-variant prior p₁₂ = 10⁻⁶. Bold PP.H4 values are ≥0.80. Full-precision probabilities, overlapping-variant counts, intermediate-prior and cohort-frequency results are retained in the accompanying detailed data.')
    if panel == 0:
        note('The five MHC genes HLA-A, HLA-DQA2, C4A, TUBB and PSMB8 were added to colocalization post hoc. Under the primary prior, only TSHR meets PP.H4 ≥0.80 in both BBJ and FinnGen; none meets both under the lower prior.')
    elif panel == 1:
        note('UKB hyperthyroidism is broader than GD; these disease definitions are not interchangeable replications. CTLA4 estimates depend mainly on rs13030124 in leave-one-out analysis. PP.H3 dominance is conditional on the single-causal-variant model and does not resolve multiple signals or exclude a gene effect.')
    else:
        note('FinnGen uses population controls and cannot isolate TED susceptibility within GD. VKORC1 MR was not estimable; NA is not a negative result. Regional colocalization was nevertheless estimable. CTLA4 estimates depend mainly on rs13030124 in leave-one-out analysis.')

doc.save(OUT / 'Tables/Table3.docx')
with (OUT / 'Detailed_Data/Table3_integrated_14genes.csv').open('w', encoding='utf-8-sig', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(integrated[0])); w.writeheader(); w.writerows(integrated)
(B / 'display_rows.json').write_text(json.dumps(display, ensure_ascii=False, indent=2), encoding='utf8')

# Correct the figure legend's cross-reference without changing figure images.
old_text = 'Table 4 gives the corresponding regional probabilities and sensitivity results.'
new_text = 'Table 3 integrates the corresponding MR estimates, regional probabilities and lower-prior results; full frequency-sensitivity results are supplied in the accompanying detailed data.'
legend = Document(SOURCE / 'Figure_Legends.docx')
changes = 0
for p in legend.paragraphs:
    for run in p.runs:
        if old_text in run.text:
            run.text = run.text.replace(old_text, new_text); changes += 1
assert changes == 1
legend.save(OUT / 'Figure_Legends.docx')
md = (SOURCE / 'Figure_Legends.md').read_text(encoding='utf8')
assert md.count(old_text) == 1
(OUT / 'Figure_Legends.md').write_text(md.replace(old_text, new_text), encoding='utf8')

readme = '''# Table 3 통합본을 반영한 표와 그림 묶음

기존 Table 3과 Table 4를 통합 Table 3으로 합쳤습니다. 본문 표는 Table 1–3, 보충표는 Table S1–S3입니다. Figures 1–3 및 S1–S3는 그대로 유지합니다.

## 통합 Table 3

- BBJ 발견 유전자 13개와 비교 대상 IGF1R을 포함한 14개 유전자를 표시합니다.
- BBJ GD, UKB 갑상선기능항진증, FinnGen Graves 안병증을 각 1쪽에 배치한 하나의 표입니다.
- 같은 행에서 MR 도구변수 수, OR·95% CI·P값, 기본 prior의 PP.H2·H3·H4 및 낮은 prior의 PP.H4를 비교할 수 있습니다.
- P값의 굵은 표시는 BBJ에서 발견 보정 기준, UKB·FinnGen에서 nominal P <0.05를 뜻합니다. 굵은 PP.H4는 ≥0.80을 뜻하며 P값과 다른 기준입니다.
- IGF1R은 발견 기준 미통과 비교 대상으로, VKORC1 FinnGen MR은 미추정으로 명시했습니다.

## 상세 수치 보존

기존 모든 상세 CSV를 그대로 유지했습니다. 새 Table3_integrated_14genes.csv는 42개 유전자·질병 조합의 MR, 전체 H0–H4, 공위치화 변이 수, 낮은 prior 및 cohort AF의 H4를 한 행에 연결합니다. 기존 Table 4의 cohort AF 열과 변이 수는 이 파일에서 확인할 수 있습니다. 504개 전체 민감도 설정도 유지했습니다. 수치를 재계산하거나 새로운 분석을 시행한 것이 아닙니다.

## 참조와 파일

Figure 3 범례의 Table 4 참조를 Table 3으로 수정했습니다. 그림 이미지 및 Tables 1–2·S1–S3는 원본과 동일합니다. All_Tables.pdf는 통합 표를 반영한 합본이며, 이 패키지에는 독립 Table 4가 없습니다.

이 패키지는 논문 본문·발표자료를 다시 빌드한 전체 제출본은 아닙니다. 기존 발표의 Tables 3/4 표기는 통합 Table 3을 뜻하도록 변경해야 하며, 원고에는 앞서 확장한 14개 유전자 분석 범위와 Table 3 및 Figure S3 인용의 동기화가 필요합니다.
'''
(OUT / 'READ_ME_KO.md').write_text(readme, encoding='utf8')
print(json.dumps({'out': str(OUT), 'rows': len(display), 'table_panels': 3, 'legend_reference_changes': changes}, ensure_ascii=False))
