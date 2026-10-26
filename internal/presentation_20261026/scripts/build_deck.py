"""Revised resident deck: same template, numbers and wording aligned with the final manuscript (f5efd35).

python3 build_deck.py <original.pptx> <figure_dir> <out.pptx>
Reuses original slides 1, 2, 4, 8, 11, 16 and 17 (edited); rebuilds the rest in the same style.
"""
import copy
import re
import sys
from pathlib import Path

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Emu, Inches, Pt

SRC, FIG, OUT = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3])
S1 = Path('/home/user/TEDGWAS/submission/candidate_20261011_review/Figures/FigureS1.png')
NAVY, BLUE, TEXT, GREY, RED = '1E2A5E', '2266AA', '1A1A1A', '666666', 'CC3333'
BOX, BOXLINE, PALE, RULE, CALL, CALLLINE = 'EEF4FF', 'B8C7DB', 'F7F8FA', 'DDDDDD', 'FFF8E8', 'B8860B'
FONT = 'Calibri'

prs = Presentation(SRC)
orig = list(prs.slides)
ORIG_NAMES = {sl.slide_id: {x.name for x in sl.shapes} for sl in orig}
LAYOUT = orig[1].slide_layout          # 'DEFAULT', no placeholders


# ---------------------------------------------------------------- text helpers
def runs(text):
    """Mini markup: **bold**, *italic*, ^superscript^, ~subscript~; \\* is a literal asterisk."""
    text = text.replace('\\*', '\x00').replace('\\~', '\x01')
    out, b, i, sup, sub = [], False, False, False, False
    for tok in re.split(r'(\*\*|\*|\^|~)', text):
        if tok == '**':
            b = not b
        elif tok == '*':
            i = not i
        elif tok == '^':
            sup = not sup
        elif tok == '~':
            sub = not sub
        elif tok:
            out.append((tok.replace('\x00', '*').replace('\x01', '~'), b, i, sup, sub))
    return out


def write_runs(p, text, size, color, bold=False, italic=False):
    for t, b, i, sup, sub in runs(text):
        r = p.add_run()
        r.text = t
        f = r.font
        f.size = Pt(size)
        f.bold = bold or b
        f.italic = italic or i
        f.color.rgb = RGBColor.from_string(color)
        f.name = FONT
        rpr = r._r.get_or_add_rPr()
        for tag in ('a:ea', 'a:cs'):
            e = rpr.find(qn(tag))
            if e is None:
                e = etree.SubElement(rpr, qn(tag))
            e.set('typeface', FONT)
        if sup:
            rpr.set('baseline', '30000')
        if sub:
            rpr.set('baseline', '-25000')


ALIGN = {'l': PP_ALIGN.LEFT, 'c': PP_ALIGN.CENTER, 'r': PP_ALIGN.RIGHT}


def fill_tf(tf, paras, size=14, color=TEXT, bold=False, align='l', after=0, before=0, line=None):
    tf.clear()
    first = True
    for item in paras:
        d = {'t': item} if isinstance(item, str) else dict(item)
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.alignment = ALIGN[d.get('align', align)]
        p.space_after = Pt(d.get('after', after))
        p.space_before = Pt(d.get('before', before))
        if d.get('line', line):
            p.line_spacing = d.get('line', line)
        ppr = p._p.get_or_add_pPr()
        if d.get('bullet'):
            lvl = d.get('level', 0)
            ppr.set('marL', str(int(Inches(0.2 + 0.22 * lvl))))
            ppr.set('indent', str(int(-Inches(0.2))))
            etree.SubElement(ppr, qn('a:buFont'), typeface='Arial')
            etree.SubElement(ppr, qn('a:buChar'), char='•' if lvl == 0 else '–')
        else:
            etree.SubElement(ppr, qn('a:buNone'))
        write_runs(p, d['t'], d.get('size', size), d.get('color', color), d.get('bold', bold), d.get('italic', False))


def tb(slide, x, y, w, h, paras, name, anchor='t', **kw):
    s = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    s.name = name
    tf = s.text_frame
    tf.word_wrap = True
    tf.auto_size = None
    for m in ('margin_left', 'margin_right', 'margin_top', 'margin_bottom'):
        setattr(tf, m, 0)
    tf.vertical_anchor = {'t': MSO_ANCHOR.TOP, 'm': MSO_ANCHOR.MIDDLE, 'b': MSO_ANCHOR.BOTTOM}[anchor]
    fill_tf(tf, paras, **kw)
    return s


def rect(slide, x, y, w, h, fill, line=None, name='box'):
    s = slide.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    s.name = name
    s.shadow.inherit = False
    s.fill.solid()
    s.fill.fore_color.rgb = RGBColor.from_string(fill)
    if line:
        s.line.color.rgb = RGBColor.from_string(line)
        s.line.width = Pt(0.75)
    else:
        s.line.fill.background()
    s.text_frame.text = ''
    return s


def retext(shape, paras, **kw):
    tf = shape.text_frame
    fill_tf(tf, paras, **kw)


def remove(shape):
    e = shape._element
    e.getparent().remove(e)


def by_name(slide, name):
    return [s for s in slide.shapes if s.name == name]


SLIDE_NO = [0]


def slide_number(slide, n):
    s = tb(slide, 9.28, 5.36, 0.5, 0.2, [''], 'Slide number', size=9, color=GREY, align='r')
    p = s.text_frame.paragraphs[0]._p
    for r in p.findall(qn('a:r')):
        p.remove(r)
    fld = etree.SubElement(p, qn('a:fld'), id='{3F1D2E4C-5B6A-4C7D-8E9F-%012d}' % n, type='slidenum')
    rpr = etree.SubElement(fld, qn('a:rPr'), lang='en-US', sz='900')
    sf = etree.SubElement(rpr, qn('a:solidFill'))
    etree.SubElement(sf, qn('a:srgbClr'), val=GREY)
    etree.SubElement(rpr, qn('a:latin'), typeface=FONT)
    etree.SubElement(fld, qn('a:t')).text = str(n)
    end = p.find(qn('a:endParaRPr'))
    if end is not None:
        p.remove(end)
        p.append(end)


def frame(slide, title):
    rect(slide, 0, 0, 10, 0.06, NAVY, name='Top bar')
    tb(slide, 0.60, 0.20, 8.94, 0.47, [title], 'Title', size=24, bold=True, color=TEXT)
    rect(slide, 0.60, 0.73, 8.80, 0.0104, RULE, name='Title rule')


def new_slide(title):
    s = prs.slides.add_slide(LAYOUT)
    frame(s, title)
    return s


def subhead(slide, text, y=0.92, size=16, w=8.8):
    return tb(slide, 0.60, y, w, 0.4, [text], 'Subheading', size=size, bold=True, color=BLUE)


def source(slide, text, y=5.30):
    return tb(slide, 0.60, y, 8.55, 0.22, [text], 'Source', size=8, color=GREY)


def pic(slide, path, x, y, name, descr, w=None, h=None):
    kw = {'width': Inches(w)} if w else {'height': Inches(h)}
    p = slide.shapes.add_picture(str(path), Inches(x), Inches(y), **kw)
    p.name = name
    p._element.nvPicPr.cNvPr.set('descr', descr)
    return p


def notes(slide, text):
    slide.notes_slide.notes_text_frame.text = text.strip()


# ---------------------------------------------------------------- tables
def table(slide, x, y, colw, rows, name, size=11, hsize=None, row_h=0.3, head_h=None, align=None, bold=(),
          fills=None):
    nr, nc = len(rows), len(rows[0])
    gf = slide.shapes.add_table(nr, nc, Inches(x), Inches(y), Inches(sum(colw)),
                                Inches((head_h or row_h) + row_h * (nr - 1)))
    gf.name = name
    t = gf.table
    tblPr = t._tbl.tblPr
    for c in list(tblPr):
        tblPr.remove(c)
    tblPr.set('firstRow', '1')
    tblPr.set('bandRow', '0')
    for j, w in enumerate(colw):
        t.columns[j].width = Inches(w)
    for i in range(nr):
        t.rows[i].height = Inches(head_h if (i == 0 and head_h) else row_h)
        for j in range(nc):
            cell = t.cell(i, j)
            head = i == 0
            al = 'l' if j == 0 or align is None else align[j]
            fill_tf(cell.text_frame, rows[i][j].split('\n'), size=(hsize or size) if head else size,
                    color='FFFFFF' if head else TEXT, bold=head or (i, j) in bold, align=al)
            tcPr = cell._tc.get_or_add_tcPr()
            for c in list(tcPr):
                tcPr.remove(c)
            for tag in ('a:lnL', 'a:lnR', 'a:lnT', 'a:lnB'):
                ln = etree.SubElement(tcPr, qn(tag), w=str(int(Pt(0.75))), cmpd='sng')
                sf = etree.SubElement(ln, qn('a:solidFill'))
                etree.SubElement(sf, qn('a:srgbClr'), val='7F8A99')
            fc = NAVY if head else (fills or {}).get((i, j), 'FFFFFF')
            sf = etree.SubElement(tcPr, qn('a:solidFill'))
            etree.SubElement(sf, qn('a:srgbClr'), val=fc)
            cell.margin_left = cell.margin_right = Inches(0.06)
            cell.margin_top = cell.margin_bottom = Inches(0.02)
            cell.vertical_anchor = MSO_ANCHOR.MIDDLE
    return gf


# ================================================================ slide 1: title (original 1)
s = orig[0]
sh = {x.name: x for x in s.shapes}
retext(sh['TSHR and IGF1R expression in Graves disease'] if 'TSHR and IGF1R expression in Graves disease' in sh
       else [x for x in s.shapes if x.has_text_frame and x.text_frame.text.startswith('TSHR and IGF1R')][0],
       ['*TSHR* and *IGF1R* in Graves disease', 'and thyroid eye disease'], size=30, bold=True, color=NAVY)
for x in s.shapes:
    if not x.has_text_frame:
        continue
    t = x.text_frame.text
    if t.startswith('Mendelian randomization'):
        retext(x, ['Mendelian randomization and colocalization'], size=20, color='555555')
    elif t.startswith('Genetic evidence linking'):
        retext(x, ['Genetic evidence from blood gene expression'], size=16.5, color=BLUE)
    elif t.startswith('Jungyul Park'):
        retext(x, ['Jungyul Park^1†^, Min-Seon Kim^2†^, Suk-Woo Yang^1\\*^'], size=14.25, bold=True, color=TEXT)
    elif t.startswith('¹ Ophthalmology'):
        retext(x, ['^1^ Department of Ophthalmology, Seoul St. Mary’s Hospital, College of Medicine, The Catholic University of Korea',
                   '^2^ Department of Ophthalmology, College of Medicine, The Catholic University of Korea'],
               size=11, color=GREY, after=2)
    elif t.startswith('* Co-corresponding'):
        x.width = Inches(6.5)
        retext(x, ['† Equal contribution (co-first authors)   \\* Corresponding author'], size=9, color=GREY)
    elif t.startswith('No financial'):
        x.left, x.width = Inches(0.60), Inches(6.0)
        retext(x, ['No conflicts of interest to disclose · Manuscript in preparation'], size=9, color=GREY)
notes(s, """
안녕하십니까. 가톨릭대학교 안과 김민선입니다.
오늘은 혈액 TSHR과 IGF1R 유전자 발현의 유전적 차이가 그레이브스병(GD), 갑상선기능항진증, 그리고 갑상선안병증(TED) 관련 진단과 어떻게 연관되는지를, 멘델리안 무작위화(MR)와 공위치 분석(colocalization)으로 비교한 연구를 말씀드리겠습니다.
먼저 강조드릴 점은, 이 연구는 치료 반응이 아니라 '질환 감수성', 즉 발병과의 연관성을 다룬다는 것입니다.
""")

# ================================================================ slide 2: background (original 2)
s = orig[1]
for x in list(s.shapes):
    if x.has_text_frame:
        t = x.text_frame.text
        if t.strip() == '2':
            remove(x)
        elif t.startswith('Genetic evidence for therapeutic'):
            x.top = Inches(0.92)
            retext(x, ['Genetic evidence for therapeutic targets'], size=16, bold=True, color=BLUE)
        elif t.startswith('Teprotumumab'):
            x.top, x.height, x.width = Inches(1.45), Inches(2.4), Inches(4.7)
            retext(x, [{'t': 'TED affects \\~40% of patients with Graves disease (GD)', 'bullet': True},
                       {'t': 'In orbital fibroblasts, TSHR forms a signaling complex with IGF-1R', 'bullet': True},
                       {'t': '**Teprotumumab**, an IGF-1R–blocking antibody, reduced proptosis in randomized trials in active TED',
                        'bullet': True},
                       {'t': 'Most approved drug mechanisms lack genetic support: an effective target need not be a '
                             'susceptibility gene', 'bullet': True}], size=13.5, color=TEXT, after=6)
        elif t.startswith('How strong'):
            x.top, x.height, x.width = Inches(3.95), Inches(0.9), Inches(4.7)
            retext(x, ['Is inherited variation in *IGF1R* expression, like *TSHR* variation, associated with GD or TED?'],
                   size=16, bold=True, color=NAVY)
for x in s.shapes:
    if x.shape_type == 13:
        x.left, x.top, x.width, x.height = Inches(5.45), Inches(1.0), Inches(3.95), Inches(3.85)
        x._element.nvPicPr.cNvPr.set('descr', 'Schematic of TED pathogenesis: thyroid-stimulating antibodies, T-cell help and '
                                              'IGF-1R on orbital fibroblasts, with teprotumumab blocking IGF-1R')
tb(s, 5.45, 4.92, 3.95, 0.3, ['Figure: add source. IGF-1R-stimulating autoantibodies (shown) remain debated.'],
   'Figure caption', size=8, color=GREY)
tb(s, 0.60, 5.02, 4.7, 0.42, ['Chin YH et al. Clin Endocrinol 2020; Tsui S et al. J Immunol 2008; Smith TJ et al. NEJM 2017; '
                              'Douglas RS et al. NEJM 2020; Nelson MR et al. Nat Genet 2015; Minikel EV et al. Nature 2024.'],
   'Source', size=8, color=GREY)
notes(s, """
갑상선안병증은 그레이브스병 환자의 약 40%에서 나타나며, 안구돌출, 복시, 시력 저하를 일으킬 수 있습니다.
안와 섬유모세포에서 TSHR은 IGF-1R과 신호 복합체를 이루고, IGF-1R 차단 항체인 teprotumumab은 활동성 TED 무작위 임상시험에서 안구돌출을 줄였습니다.
하지만 약이 듣는 표적이라고 해서 그 유전자가 발병 감수성 유전자라는 뜻은 아닙니다. 실제로 승인된 약물 기전의 상당수는 유전적 근거가 없습니다.
그래서 "IGF1R 발현의 유전적 차이도 TSHR처럼 GD나 TED와 연관되는가?"라는 질문을 하게 됩니다.
[발표 전 확인] 오른쪽 기전 그림의 출처를 슬라이드 하단 'Figure: add source' 자리에 기입하세요. 그림 속 IGF-1R 자극 자가항체는 아직 논란이 있는 개념입니다.
""")

# ================================================================ slide 3: prior studies and objectives (new)
s3 = new_slide('Background: prior studies and objectives')
rect(s3, 0.60, 0.95, 4.25, 3.15, PALE, RULE, name='Prior studies box')
tb(s3, 0.78, 1.07, 3.95, 2.95,
   [{'t': 'Previous genetic studies', 'size': 15, 'bold': True, 'color': BLUE, 'after': 6},
    {'t': '**Ji et al. 2025** (FinnGen R10): higher genetically proxied blood *TSHR* expression, lower odds of GD and strict ophthalmopathy',
     'bullet': True},
    {'t': '**Zhang et al. 2026** (FinnGen R12): shared-variant *TSHR* support with ophthalmopathy, weak with GD',
     'bullet': True},
    {'t': '**Ke et al. 2025**: proteome-wide MR prioritized other proteins', 'bullet': True},
    {'t': '**Gap:** none compared *IGF1R* with *TSHR*; all relied mainly on FinnGen', 'bullet': True, 'color': RED}],
   'Prior studies', size=13, after=7)
rect(s3, 5.05, 0.95, 4.35, 3.15, BOX, BOXLINE, name='Objectives box')
tb(s3, 5.23, 1.07, 4.0, 2.95,
   [{'t': 'Study objectives', 'size': 15, 'bold': True, 'color': BLUE, 'after': 6},
    '**1.** Compare genetic evidence linking blood *TSHR* and *IGF1R* expression with GD and ophthalmopathy endpoints, '
    'using MR and colocalization',
    '**2.** Place both genes in a screen of 2,544 potentially druggable genes, with *CTLA4* as an autoimmune comparator',
    '**3.** Test robustness to single instruments, allele frequencies and colocalization priors'],
   'Objectives', size=13, after=9)
rect(s3, 0.60, 4.25, 8.80, 0.55, CALL, CALLLINE, name='Key question box')
tb(s3, 0.78, 4.25, 8.5, 0.55, ['**Key question:** is inherited variation in blood *TSHR* and *IGF1R* expression associated with GD '
                                'and TED-related endpoints, and do the signals share one causal variant?'],
   'Key question', anchor='m', size=13, color=TEXT)
source(s3, 'Ji Q et al. Sci Rep 2025; Zhang H et al. BMC Med 2026; Ke C et al. ImmunoTargets Ther 2025.', y=5.30)
notes(s3, """
선행연구를 보면, Ji 등은 FinnGen release 10에서 혈액 TSHR 발현이 유전적으로 높을수록 GD와 strict 안병증의 odds가 낮다고 보고했습니다.
Zhang 등은 우리와 같은 FinnGen release 12 endpoint에서, TSHR이 안병증과는 공유 변이를 지지하지만 GD와는 지지가 약하다고 보고했습니다. FinnGen 자료에서는 우리 결과도 같은 패턴이지만, 같은 FinnGen R12와 eQTLGen 자료를 썼으므로 독립적인 재현은 아닙니다.
단백체 MR 연구는 다른 단백질들을 우선순위로 제시했습니다.
그러나 IGF1R을 TSHR과 직접 비교한 연구는 없었고, 대부분 FinnGen에만 의존했습니다.
그래서 본 연구는 첫째, 혈액 TSHR과 IGF1R 발현의 유전적 근거를 MR과 공위치 분석으로 비교하고, 둘째, 약물 표적 후보 유전자 2,544개 screen 안에서 CTLA4를 자가면역 비교 유전자로 두고 위치를 확인하며, 셋째, 단일 변이 의존성, 대립유전자 빈도, 공위치 사전확률에 대한 민감도를 확인했습니다.
""")

# ================================================================ slide 4: MR concept (original 4)
s = orig[3]
for x in list(s.shapes):
    if not x.has_text_frame:
        continue
    t = x.text_frame.text
    if t.strip() == '3':
        remove(x)
    elif t.startswith('Methods (1)'):
        retext(x, ['Methods (1): Two-sample Mendelian randomization'], size=24, bold=True, color=TEXT)
    elif t.startswith('Blood TSHR'):
        retext(x, ['Blood *TSHR* / *IGF1R* mRNA'], size=12, color=TEXT, align='c')
    elif t.strip() == 'Graves disease':
        retext(x, ['GD, hyperthyroidism, GO'], size=12, color=TEXT, align='c')
    elif t.strip() == 'SNP allele':
        retext(x, ['SNP allele'], size=12, color=TEXT, align='c')
    elif t.startswith('Observed associations'):
        retext(x, ['Observed associations:  ① SNP–expression relationship   ② SNP–disease relationship'], size=15.75,
               color=TEXT)
    elif t.startswith('A natural experiment'):
        x.top = Inches(0.92)
        retext(x, ['A natural experiment using genetic variation'], size=16, bold=True, color=BLUE)
    elif t.startswith('When these assumptions'):
        retext(x, ['If they hold, MR estimates the effect of lifelong, genetically proxied expression (Wald ratio or IVW)'],
               size=11.5, bold=True, color=TEXT)
    elif t.startswith('Assumptions'):
        retext(x, ['Assumptions: the SNP predicts expression, is unrelated to confounders, and affects disease only via expression'],
               size=11.5, color=TEXT)
for x in list(s.shapes):
    if x.shape_type == 13:
        remove(x)
tb(s, 6.30, 0.85, 3.1, 0.6, ['*OR* = exp(*β*~disease~ / *β*~expression~)'], 'MR formula', anchor='m', size=18, color=NAVY,
   align='c')
notes(s, """
멘델리안 무작위화는 수정(잉태) 시점에 무작위로 정해지는 유전 변이를 자연실험처럼 이용하는 방법입니다.
먼저 혈액에서 특정 유전자의 mRNA 발현량과 관련된 변이(SNP)를 찾고, 같은 변이가 질환 odds에 미치는 효과를 별도의 질환 자료에서 가져옵니다.
변이가 질환에 미치는 효과(log odds)를 변이가 발현에 미치는 효과로 나눈 뒤 지수로 바꾸면, 발현이 유전적으로 약 1 표준편차 높을 때의 odds ratio가 됩니다. 변이가 하나면 Wald ratio, 여러 개면 IVW로 합칩니다.
여기서 노출은 단백질이 아니라 혈액 mRNA 발현이고, 결과는 치료 반응이 아니라 진단입니다.
이 추정이 의미를 가지려면, 변이가 발현을 잘 예측하고, 교란요인과 무관하며, 발현을 거치지 않고는 질환에 영향을 주지 않아야 합니다.
그리고 MR P값이 아주 작더라도 그 자체가 발현이 질환을 매개한다는 증거는 아니어서, 다음에 설명할 공위치 분석을 함께 봅니다.
""")

# ================================================================ slide 5: datasets (new)
s5 = new_slide('Methods (2): Exposure and disease datasets')
subhead(s5, 'Summary statistics only; no individual-level genotypes')
rows = [['Dataset', 'Endpoint (label used here)', 'Ancestry', 'Cases / controls', 'Role'],
        ['eQTLGen', 'Blood cis-eQTL (exposure)', 'Predominantly European', 'Up to 31,684\nparticipants', 'Exposure'],
        ['BBJ', 'Graves disease (BBJ GD)', 'East Asian', '2,809 / 172,656', '**Discovery**'],
        ['FinnGen R12', 'E4_GRAVES_STRICT\n(FinnGen GD)', 'Finnish', '3,962 / 496,386', 'Same-disease comparison\n(post hoc)'],
        ['UK Biobank', 'Self-reported hyperthyroidism or\nthyrotoxicosis (UKB hyperthyroidism)', 'Not established as\nexclusively European', '3,731 / 480,867',
         'Broader-phenotype comparison'],
        ['FinnGen R12', 'GRAVES_OPHT (FinnGen GO)', 'Finnish', '858 / 499,490', 'TED-enriched comparison'],
        ['FinnGen R12', 'E4_GRAVES_OPHT_STRICT\n(FinnGen GO strict)', 'Finnish', '753 / 499,595',
         'Definition sensitivity\n(post hoc)']]
table(s5, 0.60, 1.30, [1.05, 2.45, 1.75, 1.45, 2.10], rows, 'Datasets table', size=11, hsize=11.5, row_h=0.42,
      head_h=0.36, align=['l', 'l', 'l', 'c', 'l'])
rect(s5, 0.60, 4.48, 8.80, 0.72, PALE, RULE, name='Dataset caveats box')
tb(s5, 0.75, 4.53, 8.5, 0.62,
   [{'t': 'The three FinnGen endpoints share participants: not independent replications, not meta-analyzed', 'bullet': True},
    {'t': 'Ophthalmopathy controls are all other participants, including GD without recorded eye disease, so GO estimates do not '
          'estimate TED risk among patients with GD', 'bullet': True}], 'Dataset caveats', size=11.5, after=2)
source(s5, 'GO, Graves ophthalmopathy; GD, Graves disease; BBJ, Biobank Japan; UKB, UK Biobank; eQTL, expression quantitative trait locus.')
notes(s5, """
노출 자료는 혈액 eQTL 연구인 eQTLGen으로, 최대 31,684명, 주로 유럽계입니다.
질환 자료는 다섯 개입니다. BBJ 그레이브스병이 발견(discovery) 자료이고, FinnGen GD는 같은 질환 비교, UK Biobank 갑상선기능항진증은 더 넓은 표현형 비교, FinnGen 안병증(GO)은 TED가 많이 포함된 비교, strict 안병증은 정의 민감도 분석입니다. FinnGen GD와 strict 안병증은 나중에(post hoc) 추가했습니다.
UK Biobank의 조상(ancestry)은 유럽계만으로 확인되지 않았고, FinnGen은 핀란드인 자료입니다.
주의할 점이 두 가지 있습니다. 첫째, FinnGen 세 endpoint는 참여자가 겹치기 때문에 독립적인 재현이 아니며, 메타분석하지 않았습니다.
둘째, 안병증 분석의 대조군은 나머지 모든 참여자, 즉 눈 증상이 기록되지 않은 GD 환자도 포함하기 때문에, GD 환자 중에서 TED가 생길 위험을 추정하는 것이 아닙니다.
한 사람이 모든 자료에 들어 있는 것이 아니라, 서로 다른 연구의 SNP별 요약통계를 결합하기 때문에 two-sample MR입니다.
""")

# ================================================================ slide 6: instruments (original 8)
s = orig[7]
for x in list(s.shapes):
    if not x.has_text_frame:
        continue
    t = x.text_frame.text
    if t.strip() == '5':
        remove(x)
    elif t.startswith('Methods (2)'):
        retext(x, ['Methods (3): Instruments and druggable-gene screen'], size=24, bold=True, color=TEXT)
    elif t.startswith('Selection of 2,544'):
        x.top = Inches(0.92)
        retext(x, ['Selecting expression instruments for 2,544 genes'], size=16, bold=True, color=BLUE)
    elif t.startswith('Druggable-gene list'):
        retext(x, ['Potentially druggable genome (Finan et al. 2017)'], size=12.75, color=GREY)
    elif t.startswith('31,684 participants'):
        retext(x, ['Up to 31,684 people, predominantly European'], size=12.75, color=GREY)
    elif t.startswith('① cis-eQTLs'):
        x.top = Inches(3.10)
        retext(x, ['① cis-eQTLs within 1 Mb of the gene, *P* < 5×10^−8^'], size=15, color=TEXT)
    elif t.startswith('② LD clumping'):
        x.top = Inches(3.42)
        retext(x, ['② LD clumping: r² < 0.001 within 10 Mb (1000 Genomes European)'], size=15, color=TEXT)
    elif t.startswith('2,544 genes /'):
        x.top = Inches(4.42)
        retext(x, ['2,544 genes / 6,135 SNP–gene instrument records'], size=17.25, bold=True, color=TEXT)
    elif t.startswith('↓'):
        x.top = Inches(4.08)
    elif t.startswith('Select SNPs'):
        retext(x, ['Select SNPs associated with expression of each gene'], size=15.75, bold=True, color=TEXT)
for x in s.shapes:
    if x.name == '직사각형 13' and not x.text_frame.text.strip():
        x.height = Inches(1.42)
tb(s, 0.80, 3.74, 8.36, 0.28, ['③ Instrument strength check: all F > 10 (minimum F = 29.7)'], 'Instrument strength', size=15,
   color=TEXT)
geo = {'직사각형 6': (0.60, 4.00), '직사각형 9': (5.40, 4.00), '약물 표적 후보 4,462개 유전자': (0.75, 3.70),
       'Finan et al.': (0.75, 3.70), 'eQTLGen 혈액 유전자 발현': (5.53, 3.75), '최대 31,684명, 주로 유럽계': (5.53, 3.75),
       '×': (4.72, 0.56), '직선 화살표 연결선 24': (2.60, 2.40), '직선 화살표 연결선 23': (5.00, 2.40)}
moved = set()
for x in s.shapes:
    for k, (left, width) in geo.items():
        if x.name.startswith(k):
            x.left, x.width = Inches(left), Inches(width)
            moved.add(k)
assert moved == set(geo), set(geo) - moved
rect(s, 0.60, 4.86, 8.80, 0.40, PALE, RULE, name='Threshold box')
tb(s, 0.75, 4.86, 8.6, 0.40, ['Screen threshold (BBJ): *P* < 0.05/2,544 = 1.965×10^−5^   ·   Focused FinnGen GD threshold: '
                              '*P* < 0.05/13 = 0.00385'], 'Thresholds', anchor='m', size=11.5, color=TEXT)
notes(s, """
어떤 SNP를 MR 도구로 썼는지 설명드리겠습니다.
기존 문헌의 약물 표적 후보 유전자 4,462개에서 출발해, eQTLGen 혈액 자료에서 유전자 주변 1 Mb 이내에서 발현과 유의하게 연관된 SNP(P < 5×10⁻⁸)를 골랐습니다.
서로 함께 유전되는 SNP는 LD clumping(r² < 0.001)으로 정리했고, 도구 강도를 나타내는 F 통계량은 모두 10을 넘었습니다(최소 29.7).
그 결과 2,544개 유전자, 6,135개 SNP–유전자 쌍을 얻었습니다.
BBJ screen 기준은 2,544개 전체에 대한 Bonferroni 보정, 즉 P < 1.965×10⁻⁵이고, BBJ에서 찾은 13개 유전자는 FinnGen GD에서 P < 0.05/13이면서 같은 방향인지로 다시 확인했습니다.
""")

# ================================================================ slide 7: colocalization (original 11)
s = orig[10]
for x in list(s.shapes):
    if x.has_text_frame:
        t = x.text_frame.text
        if x.name in ('PP.H4 ≥ 0.80', '본 연구의 강한 공위치 근거 기준') or x.name.startswith('PP.H4는') or \
                t.startswith('PP.H4 ≥ 0.80 :'):
            remove(x)
        elif t.startswith('Methods (3)'):
            retext(x, ['Methods (4): Colocalization'], size=24, bold=True, color=TEXT)
        elif t.startswith('Do expression'):
            x.top = Inches(0.92)
            retext(x, ['Do the expression and disease signals share the same causal variant?'], size=16, bold=True, color=BLUE)
        elif t.startswith('Bayesian analysis'):
            x.top, x.height = Inches(1.30), Inches(0.62)
            retext(x, ['An MR association can be highly significant even when two nearby, correlated variants drive the two signals. '
                       'coloc.abf (at most one causal variant per trait) gives posterior probabilities for five hypotheses.'],
                   size=13, color=TEXT)
    elif x.has_table:
        remove(x)
rows = [['Hypothesis', 'Interpretation in this study'],
        ['H0', 'No association with gene expression or disease risk'],
        ['H1', 'Association with disease risk only'],
        ['H2', 'Association with blood gene expression only'],
        ['H3', 'Both associated, but driven by **different causal variants**'],
        ['H4', 'Both associated and share the **same causal variant**']]
table(s, 0.60, 2.02, [1.20, 4.55], rows, 'Hypothesis table', size=13, hsize=13, row_h=0.42, head_h=0.40,
      bold={(i, 0) for i in range(1, 6)}, fills={(5, 0): 'FFE9A8', (5, 1): 'FFE9A8'})
rect(s, 6.55, 2.02, 2.85, 2.50, CALL, CALLLINE, name='Coloc settings box')
tb(s, 6.70, 2.12, 2.6, 2.30,
   [{'t': 'Settings', 'bold': True, 'color': NAVY, 'size': 14, 'after': 4},
    {'t': 'Priors: p1 = p2 = 10^−4^, p12 = 10^−5^', 'bullet': True},
    {'t': 'Sensitivity: p12 = 5×10^−6^ and 10^−6^', 'bullet': True},
    {'t': 'PP.H4 ≥ 0.80: descriptive reporting threshold', 'bullet': True},
    {'t': 'PP.H4 is a probability, not a *P* value', 'bullet': True}], 'Coloc settings', size=12, after=4)
source(s, 'Giambartolomei C et al. PLoS Genet 2014; Wallace C. PLoS Genet 2020.')
notes(s, """
공위치 분석이 왜 필요한지 말씀드리겠습니다.
가까이 있는 변이들은 함께 유전되기 때문에, 발현을 조절하는 변이와 질환과 연관된 변이가 서로 다르더라도 MR 연관성은 매우 유의하게 나올 수 있습니다.
coloc은 한 지역에서 형질마다 원인 변이가 최대 하나라고 가정하고, 다섯 가설의 사후확률을 계산합니다. H3는 서로 다른 변이, H4는 같은 하나의 변이를 공유한다는 가설입니다.
사전확률은 p1 = p2 = 10⁻⁴, p12 = 10⁻⁵로 두었고, p12를 5×10⁻⁶과 10⁻⁶로 낮춰 민감도를 봤습니다. PP.H4 0.80 이상을 보고 기준으로 썼습니다.
MR 연관이 주 결과이고 공위치 분석은 이를 해석하는 데 썼으며, 이 강조점은 사후(post hoc)에 정했습니다.
PP.H4는 P값이 아닙니다. PP.H4가 높다는 것은 이 모형에서 공유 변이가 지지된다는 뜻이지, "발현이 질병의 원인일 확률"을 뜻하지는 않습니다.
""")

# ================================================================ slide 8: screen results (new)
s8 = new_slide('Results (1): Druggable-gene screen in BBJ GD')
pic(s8, FIG / 'slide_fig1_manhattan.png', 0.45, 0.88, 'Figure 1B', 'Manhattan-style plot of BBJ GD MR P values for 2,234 genes; '
    '13 pass P < 0.05/2,544; IGF1R does not', h=4.3)
rect(s8, 6.62, 0.95, 2.78, 4.22, BOX, BOXLINE, name='Screen summary box')
tb(s8, 6.80, 1.05, 2.45, 4.02, anchor='m', paras=
   [{'t': '2,234', 'size': 28, 'bold': True, 'color': NAVY},
    {'t': 'genes estimable in BBJ', 'size': 13, 'after': 8},
    {'t': '13 genes', 'size': 28, 'bold': True, 'color': NAVY},
    {'t': 'passed *P* < 0.05/2,544', 'size': 13, 'after': 8},
    {'t': '*TSHR* and *CTLA4* passed; *IGF1R* did not (*P* = 0.0212)', 'size': 13, 'bold': True, 'color': RED, 'after': 8},
    {'t': 'Five MHC genes (linked region) and three 16p11.2 genes (one regional signal) cannot each count as independent '
          'targets', 'size': 11.5, 'color': GREY, 'after': 4},
    {'t': 'Only 14.6% of genes had 80% power for OR 1.5 at this threshold', 'size': 11.5, 'color': GREY}],
   name='Screen summary', size=13)
source(s8, 'Each point is one gene (strongest instrument position). Dashed line: screen threshold. Source: Figure 1B.')
notes(s8, """
BBJ에서 추정 가능한 2,234개 유전자 중 13개가 screen 기준을 통과했습니다.
여기에 TSHR과 CTLA4는 포함됐지만 IGF1R은 포함되지 않았습니다. IGF1R의 BBJ P값은 0.0212로 nominal 수준이었습니다.
13개 중 다섯 개는 MHC 영역, 세 개는 16번 염색체 16p11.2 영역 유전자로, 같은 지역 신호를 공유할 수 있어 각각을 독립적인 발견으로 보기는 어렵습니다.
참고로 이 엄격한 기준에서는 OR 1.5를 80% 검정력으로 찾을 수 있는 유전자가 14.6%뿐이어서, 기준을 통과하지 못했다고 효과가 없다는 뜻은 아닙니다.
그래프에서 점이 높을수록 연관 근거가 강하다는 뜻이지, 효과가 크거나 위험을 높이는 방향이라는 뜻은 아닙니다.
""")

# ================================================================ slide 9: overview figure (new)
s9 = new_slide('Results (2): *TSHR*, *IGF1R* and *CTLA4* across five datasets')
pic(s9, FIG / 'slide_fig2_overview.png', 0.6, 0.82, 'Figure 2', 'Forest plot of MR estimates and colocalization bars for TSHR, '
    'IGF1R and CTLA4 in five datasets', w=8.8)
source(s9, 'OR per \\~1 SD higher genetically proxied blood expression. Dashed line: OR = 1; dotted line: PP.H4 = 0.80. *TSHR* and '
           '*IGF1R* values on the next two slides; all values in Table 2. Source: Figure 2.', y=5.28)
notes(s9, """
세 유전자의 결과를 한 번에 보여드리는 그림입니다. 왼쪽이 MR 추정치, 오른쪽이 공위치 사후확률이며, 오른쪽의 파란 막대가 H4, 즉 공유 변이 확률이고 점선이 0.80입니다.
TSHR은 다섯 자료 모두에서 OR이 1보다 낮았습니다. 그러나 공위치 분석에서는 BBJ와 두 안병증 자료에서 H4가 높았던 반면, FinnGen GD와 UK Biobank에서는 주황색 H3, 즉 서로 다른 변이가 우세했습니다.
IGF1R은 다섯 자료 모두 OR이 1보다 컸지만 신뢰구간이 1에 가깝거나 걸쳐 있었고, 회색의 H2, 즉 발현 신호만 있는 패턴이 주로 나왔습니다.
CTLA4는 비교 유전자로 모든 자료에서 OR이 1보다 낮았습니다. BBJ(P = 5.45×10⁻¹⁵, PP.H4 0.201)와 FinnGen GD(PP.H4 0.254)에서는 주황색 H3가 우세했지만, P가 nominal 수준(0.0022–0.0118)이었던 UK Biobank와 두 안병증 자료에서는 PP.H4가 0.95 이상이었습니다. 연관의 강도와 공유 변이 근거는 서로 다른 질문에 답합니다.
TSHR과 IGF1R의 구체적인 숫자는 다음 두 슬라이드에서 보겠습니다.
""")

# ================================================================ slide 10: TSHR (new)
s10 = new_slide('Results (3): *TSHR* across five datasets')
subhead(s10, 'Lower odds in every dataset; distinct variants favored in UKB and FinnGen GD', size=16)
rows = [['Dataset', 'OR (95% CI)', '*P*', 'PP.H3', 'PP.H4', 'PP.H4, p12 = 10^−6^'],
        ['BBJ GD', '0.12 (0.07–0.21)', '1.09×10^−14^', '0.049', '0.951', '0.661'],
        ['FinnGen GD', '0.08 (0.05–0.12)', '3.19×10^−32^', '0.736', '0.264', '0.035'],
        ['UKB hyperthyroidism', '0.09 (0.06–0.14)', '8.77×10^−28^', '0.774', '0.226', '0.028'],
        ['FinnGen GO', '0.10 (0.04–0.24)', '2.82×10^−7^', '0.014', '0.986', '0.875'],
        ['FinnGen GO strict', '0.15 (0.06–0.39)', '9.02×10^−5^', '0.082', '0.891', '0.451']]
bold = {(1, 4), (4, 4), (5, 4), (4, 5), (2, 3), (3, 3)}
fills = {(2, 3): 'FBEBDD', (3, 3): 'FBEBDD', (1, 4): 'E3EEF7', (4, 4): 'E3EEF7', (5, 4): 'E3EEF7', (4, 5): 'E3EEF7'}
table(s10, 0.60, 1.38, [2.0, 1.65, 1.35, 1.05, 1.05, 1.7], rows, 'TSHR table', size=12.5, hsize=12, row_h=0.36,
      head_h=0.36, align=['l', 'c', 'c', 'c', 'c', 'c'], bold=bold, fills=fills)
tb(s10, 0.60, 3.70, 8.8, 1.5,
   [{'t': 'One instrument in every dataset, **rs179252**: its G allele raises blood *TSHR* by 0.105 SD (F = 176), with per-allele '
          'disease ORs of **0.77–0.82**; the per-SD ORs rescale this moderate effect', 'bullet': True},
    {'t': 'Shared variant favored (PP.H4 ≥ 0.80) in BBJ and both GO endpoints; **PP.H3 highest in UKB and FinnGen GD**, the two '
          'datasets with the most cases', 'bullet': True},
    {'t': 'At p12 = 10^−6^, BBJ and GO strict fell below 0.80; only FinnGen GO stayed above it', 'bullet': True}],
   'TSHR findings', size=13, after=5)
source(s10, 'OR per approximately 1 SD higher genetically proxied blood expression (Wald ratio). Blue: PP.H4 ≥ 0.80; orange: PP.H3 '
            'highest. Source: Table 2, Table S3.')
notes(s10, """
TSHR의 숫자입니다. 유전적으로 예측된 혈액 TSHR 발현이 높을수록 BBJ GD의 odds가 낮았고(OR 0.12, 95% CI 0.07–0.21, P = 1.09×10⁻¹⁴), FinnGen GD에서도 낮았습니다(OR 0.08, 95% CI 0.05–0.12, P = 3.19×10⁻³²). 갑상선기능항진증과 두 안병증 자료에서도 OR은 모두 1보다 낮았습니다.
OR 0.12가 매우 커 보이지만, 이것은 발현 1 표준편차당 값입니다. 다섯 자료 모두 rs179252 하나를 썼고, 이 변이의 G 대립유전자 하나당 질환 OR은 0.77–0.82 정도의 중간 크기 효과입니다. 이 효과를 발현 1 SD로 환산한 것이 표의 OR입니다.
공위치 분석에서는 BBJ, FinnGen GO, strict GO에서 PP.H4가 0.80을 넘었지만, MR이 가장 강했던 FinnGen GD(0.264)와 UK Biobank(0.226)에서는 서로 다른 변이(H3)가 우세했습니다. 이 두 자료가 환자 수가 가장 많습니다.
또 공유 변이의 사전확률을 10⁻⁶로 낮추면 BBJ는 0.951에서 0.661로, strict GO는 0.891에서 0.451로 낮아져, 이 두 자료의 지지는 사전확률에 민감합니다.
즉 연관성은 일관되지만, 같은 원인 변이가 발현과 질환을 함께 설명한다는 근거는 자료마다 달랐습니다.
""")

# ================================================================ slide 11: IGF1R (new)
s11 = new_slide('Results (4): *IGF1R* across five datasets')
subhead(s11, 'All ORs above 1, but nominal only in BBJ and UKB; FinnGen GD near the null', size=16)
rows = [['Dataset', 'IVs', 'OR (95% CI)', '*P*', 'PP.H2', 'PP.H4'],
        ['BBJ GD', '4', '1.56 (1.07–2.28)', '0.0212', '0.690', '0.073'],
        ['FinnGen GD', '3', '1.07 (0.86–1.33)', '0.539', '0.731', '0.010'],
        ['UKB hyperthyroidism', '4', '1.35 (1.07–1.70)', '0.0117', '0.400', '0.404'],
        ['FinnGen GO', '3', '1.41 (0.85–2.33)', '0.182', '0.623', '0.032'],
        ['FinnGen GO strict', '3', '1.52 (0.92–2.49)', '0.0996', '0.571', '0.038']]
fills = {(2, 2): 'E3EEF7', (2, 3): 'E3EEF7'}
table(s11, 0.60, 1.38, [2.15, 0.75, 2.0, 1.3, 1.3, 1.3], rows, 'IGF1R table', size=12.5, hsize=12, row_h=0.36, head_h=0.36,
      align=['l', 'c', 'c', 'c', 'c', 'c'], bold={(2, 2), (2, 3)}, fills=fills)
tb(s11, 0.60, 3.70, 8.8, 1.5,
   [{'t': 'None passed the screen threshold. **FinnGen GD**, the most precise estimate, had 80% power to detect an OR of 1.36 '
          'at *P* < 0.05', 'bullet': True},
    {'t': 'Weighted median: BBJ 1.56 (*P* = 0.0288), UKB 1.42 (*P* = 0.00244); in BBJ and UKB, two of four single-instrument '
          'omissions gave *P* > 0.05', 'bullet': True},
    {'t': 'Colocalization: **PP.H2 (expression signal only) highest in four of five datasets**; no PP.H4 ≥ 0.80', 'bullet': True}],
   'IGF1R findings', size=13, after=5)
source(s11, 'Multiplicative random-effects IVW. Source: Table 2, Figure S2; Supplementary Dataset (MR_estimators_original).')
notes(s11, """
IGF1R입니다. 다섯 자료 모두 OR은 1보다 컸지만, P < 0.05인 nominal 연관은 BBJ GD(OR 1.56, 95% CI 1.07–2.28, P = 0.0212)와 UK Biobank 갑상선기능항진증(OR 1.35, P = 0.0117)뿐이었고, screen 기준을 통과한 자료는 없었습니다.
가장 정밀한 추정치인 FinnGen GD는 OR 1.07(95% CI 0.86–1.33, P = 0.539)로 1에 가까웠습니다. 이 비교는 OR 1.36 정도의 효과를 80% 검정력으로 찾을 수 있었습니다. 즉 큰 효과는 어렵지만 작은 효과는 배제할 수 없습니다.
weighted median도 BBJ 1.56, UKB 1.42로 비슷했습니다. BBJ와 UKB에서 도구 변이를 하나씩 빼 보면 네 번 중 두 번은 P가 0.05를 넘었지만, 대부분 OR은 전체 추정치와 비슷했습니다. 유의성 소실은 정밀도 감소 때문일 수도 있습니다(UKB에서 rs2654980을 빼면 OR 1.35→1.16).
공위치 분석에서는 다섯 자료 중 네 개에서 H2, 즉 발현 신호만 있고 질환 신호는 약한 패턴이 가장 높았습니다.
따라서 혈액 발현을 통한 IGF1R의 유전적 근거는 제한적이지만, 유의하지 않다는 것이 유전적 역할이 없다는 뜻은 아닙니다.
""")

# ================================================================ slide 12: 13 BBJ discoveries (new)
s12 = new_slide('Results (5): Thirteen BBJ discoveries in FinnGen GD')
rows = [['Gene', 'Chr', 'BBJ OR', 'BBJ *P*', 'FinnGen GD OR', 'FinnGen GD *P*', 'Focused support', 'FinnGen GD PP.H4'],
        ['*HLA-A*', '6 (MHC)', '1.98', '2.60×10^−23^', '1.08', '0.648', 'No', '<0.001'],
        ['*HLA-DQA2*', '6 (MHC)', '2.39', '8.30×10^−17^', '1.90', '1.41×10^−17^', 'Yes', '<0.001'],
        ['*CTLA4*', '2', '0.18', '5.45×10^−15^', '0.17', '2.15×10^−6^', 'Yes', '0.254'],
        ['*TSHR*', '14', '0.12', '1.09×10^−14^', '0.08', '3.19×10^−32^', 'Yes', '0.264'],
        ['*C4A*', '6 (MHC)', '0.45', '2.23×10^−10^', '0.41', '4.97×10^−5^', 'Yes', '<0.001'],
        ['*HSD3B7*', '16p11.2', '0.29', '2.04×10^−7^', '0.90', '0.372', 'No', '0.011'],
        ['*TUBB*', '6 (MHC)', '1.49', '3.14×10^−7^', '3.28', '2.68×10^−7^', 'Yes', '<0.001'],
        ['*VKORC1*', '16p11.2', '0.13', '6.35×10^−7^', 'NA', 'NA', 'NA', '0.011'],
        ['*TNFSF14*', '19', '0.63', '1.48×10^−6^', '0.77', '0.00361', 'Yes', '0.303'],
        ['*PRSS36*', '16p11.2', '6.00', '2.02×10^−6^', '1.40', '0.0490', 'No', '0.063'],
        ['*MAPKAPK5*', '12', '0.977', '5.28×10^−6^', '0.990', '0.751', 'No', '0.041'],
        ['*PSMB8*', '6 (MHC)', '1.25', '6.86×10^−6^', '0.96', '0.512', 'No', '<0.001'],
        ['*IFNGR1*', '6', '2.10', '9.41×10^−6^', '1.62', '0.00330', 'Yes', '0.313']]
yes = {(i, 6): 'E3EEF7' for i in range(1, 14) if rows[i][6] == 'Yes'}
table(s12, 0.60, 0.86, [1.15, 0.95, 0.8, 1.2, 1.15, 1.2, 1.05, 1.3], rows, 'Discoveries table', size=10.5, hsize=10.5,
      row_h=0.235, head_h=0.42, align=['l', 'c', 'c', 'c', 'c', 'c', 'c', 'c'], fills=yes)
tb(s12, 0.60, 4.45, 8.8, 0.80,
   [{'t': 'Focused support (same direction, *P* < 0.05/13 in FinnGen GD): 7 of 13. *HLA-A*, the smallest BBJ *P*, was not supported',
     'bullet': True},
    {'t': 'No gene reached PP.H4 ≥ 0.80 in FinnGen GD; *HLA-DQA2* and *TUBB* < 0.001 in every dataset (MHC: H3 dominant; complex LD)',
     'bullet': True},
    {'t': '*CTLA4* FinnGen and UKB estimates depend on rs13030124 (omitting it: *P* = 0.323–0.972)', 'bullet': True}],
   'Discoveries findings', size=11.5, after=2)
source(s12, 'Listed by BBJ *P*. NA: no retained instrument, not a negative result. 95% CIs in Table 3.')
notes(s12, """
BBJ에서 찾은 13개 유전자를 FinnGen GD에서 다시 확인한 표입니다. 같은 방향이면서 P < 0.05/13을 만족한 유전자는 7개로, TSHR, CTLA4, HLA-DQA2, TUBB, C4A, TNFSF14, IFNGR1입니다.
BBJ에서 P값이 가장 작았던 HLA-A는 FinnGen GD에서 유의하지 않아, 발견 단계 P값이 가장 작다고 재현이 가장 잘 되는 것은 아니었습니다.
MHC 영역 유전자들은 공위치 분석에서 모두 H3, 즉 서로 다른 변이 패턴이 우세했고, HLA-DQA2와 TUBB는 모든 자료에서 PP.H4가 0.001 미만이었습니다. 다만 MHC는 LD가 복잡하고 여러 신호가 있을 수 있어, 원인 변이를 하나로 가정한 모형에서는 H3가 모형 불일치를 반영할 수 있습니다. FinnGen GD에서는 어떤 유전자도 PP.H4 0.80에 도달하지 못했습니다.
비교 유전자인 CTLA4는 FinnGen과 UKB 추정치가 rs13030124 하나에 크게 의존했습니다. 이 변이를 빼면 P = 0.323–0.972로 유의성이 사라집니다.
그래서 MHC 연관, 단일 변이 의존성, 16p11.2의 공통 지역 신호를 고려하면 이 유전자들을 각각 독립적인 원인 표적으로 셀 수는 없습니다.
""")

# ================================================================ slide 13: discussion TSHR (new)
s13 = new_slide('Discussion (1): *TSHR* and disease susceptibility')
subhead(s13, 'A consistent inverse association; expression-mediated causality is not established', size=16)
rect(s13, 0.60, 1.38, 4.3, 3.62, PALE, RULE, name='Direction box')
tb(s13, 0.78, 1.48, 3.98, 3.42,
   [{'t': 'Why lower odds?', 'size': 15, 'bold': True, 'color': BLUE, 'after': 6},
    {'t': 'rs179252-G is linked to **rs179247-G** (r = 0.978 European; 0.996 East Asian), the GD-protective allele associated '
          'with greater intrathymic *TSHR* transcription', 'bullet': True},
    {'t': 'Fits central (thymic) tolerance; the thymic step was **not tested**, and neither variant is a proven causal regulator',
     'bullet': True},
    {'t': 'Blood transcript, tissue expression and antibody-driven receptor activation are different exposures', 'bullet': True},
    {'t': 'Lifelong inherited differences, not advice to raise *TSHR* or a TSHR-directed treatment effect', 'bullet': True}],
   'Direction', size=13, after=6)
rect(s13, 5.10, 1.38, 4.3, 3.62, BOX, BOXLINE, name='Association box')
tb(s13, 5.28, 1.48, 3.98, 3.42,
   [{'t': 'Association is not colocalization', 'size': 15, 'bold': True, 'color': BLUE, 'after': 6},
    {'t': 'Very small MR *P* in UKB and FinnGen GD coexisted with distinct-variant patterns (also *CTLA4* in BBJ and FinnGen GD)',
     'bullet': True},
    {'t': 'Regional checks excluded a missing peak but did not identify the causal variant', 'bullet': True},
    {'t': 'GO controls include GD without recorded eye disease, and GO per-allele ORs resembled those for GD: '
          'these may mostly reflect **GD susceptibility**; higher GO PP.H4 is not TED-specific evidence', 'bullet': True}],
   'Association vs colocalization', size=13, after=6)
source(s13, 'Brand OJ et al. Hum Mol Genet 2009; Colobran R et al. Hum Mol Genet 2011; Marín-Sánchez A et al. Front Immunol 2019.')
notes(s13, """
첫 번째 고찰은 TSHR입니다.
TSHR 발현이 유전적으로 높을수록 GD odds가 낮다는 방향은, 자극 항체가 TSHR을 활성화해서 병이 생긴다는 기존 지식과 반대로 보일 수 있습니다.
하지만 rs179252의 G 대립유전자는, GD 위험이 낮은 쪽(보호 대립유전자)으로 보고된 rs179247-G와 강하게 연관되어 있고(유럽 r = 0.978, 동아시아 0.996), rs179247-G는 흉선에서 TSHR 전사가 더 많은 것과 관련이 있다고 보고되었습니다. 흉선에서 자가항원이 충분히 발현되어야 면역관용이 생긴다는 설명과 맞지만, 본 연구는 흉선 단계를 직접 검증하지 않았고, 두 변이 모두 원인 조절 변이로 증명되지는 않았습니다.
또한 혈액 mRNA 양, 조직 발현량, 자가항체에 의한 수용체 활성화는 서로 다른 노출이며, 이 OR은 평생의 유전적 차이를 뜻할 뿐, TSHR 발현을 높이라거나 TSHR 표적 치료의 효과를 뜻하지 않습니다.
그리고 UK Biobank와 FinnGen GD에서는 MR P값이 매우 작았지만 공위치 분석은 서로 다른 변이를 지지했습니다. 연관의 강도와 공유 변이 근거는 서로 대체할 수 없습니다.
안병증 자료의 대조군에는 눈 증상이 기록되지 않은 GD 환자도 포함되고, 안병증의 대립유전자당 OR도 GD와 비슷해, 이 결과는 주로 GD 감수성을 반영할 가능성이 있습니다. 안병증의 PP.H4가 더 높았던 것도 가장 작고 신호가 약한 두 자료에서 나온 것이어서, TED 특이적 효과를 뜻하지 않습니다.
""")

# ================================================================ slide 14: discussion IGF1R (new)
s14 = new_slide('Discussion (2): *IGF1R* findings and IGF-1R blockade')
subhead(s14, 'Limited blood-expression evidence for *IGF1R*; not an argument against IGF-1R blockade', size=16)
tb(s14, 0.60, 1.42, 8.8, 2.2,
   [{'t': 'Colocalization mostly favored an expression-only signal, consistent with weak regional disease signals', 'bullet': True},
    {'t': 'The near-null FinnGen GD estimate argues against large, but not modest, effects; non-significance is not evidence of '
          'no inherited role', 'bullet': True},
    {'t': 'Tissue effects remain open: single-instrument GTEx estimates, no orbital eQTL data', 'bullet': True},
    {'t': 'Significance for one gene but not the other is not a formal test of different effects', 'bullet': True}],
   'IGF1R interpretation', size=13.5, after=7)
rect(s14, 0.60, 3.30, 8.80, 1.15, CALL, CALLLINE, name='Treatment box')
tb(s14, 0.78, 3.38, 8.45, 1.0,
   [{'t': 'Disease susceptibility ≠ treatment response', 'size': 15, 'bold': True, 'color': NAVY, 'after': 4},
    {'t': 'Teprotumumab reduced proptosis in active TED and in long-duration, low-activity TED', 'bullet': True},
    {'t': 'TSHR–IGF-1R crosstalk concerns receptor signaling in established disease, not inherited susceptibility',
     'bullet': True}],
   'Treatment', size=13.5, after=4)
source(s14, 'Douglas RS et al. NEJM 2020; Douglas RS et al. J Clin Endocrinol Metab 2024; Krieger CC et al. Endocrinology 2019.')
notes(s14, """
두 번째 고찰은 IGF1R입니다.
혈액 발현으로 본 IGF1R의 유전적 근거는 제한적이었습니다. 공위치 분석은 주로 발현 신호만 있는 패턴이었고, 일부 도구 변이를 빼면 nominal 유의성이 사라졌습니다.
가장 정밀한 FinnGen GD 추정치가 1에 가까워 큰 효과는 어렵지만, 작은 효과는 배제할 수 없고, 유의하지 않다는 것이 유전적 역할이 없다는 뜻은 아닙니다. 조직별 효과도 아직 열려 있습니다.
또한 한 유전자는 유의하고 다른 유전자는 유의하지 않다는 것만으로 두 효과가 다르다고 통계적으로 검정한 것은 아닙니다.
가장 중요한 점은, 발병 감수성과 이미 생긴 질환에서의 치료 반응은 서로 다른 질문이라는 것입니다. Teprotumumab은 활동성 TED뿐 아니라 오래되고 활동성이 낮은 TED에서도 안구돌출을 줄였고, TSHR–IGF-1R 상호작용은 수용체 신호에 관한 것이지 유전 감수성에 관한 것이 아닙니다.
따라서 이번 IGF1R 결과는 IGF-1R 차단 치료에 반대하는 근거가 아닙니다.
""")

# ================================================================ slide 15: conclusions (original 16, rebuilt)
s = orig[15]
for x in list(s.shapes):
    if x.name not in ('Shape 0',) and not (x.has_text_frame and x.text_frame.text.startswith('Conclusions')) \
            and x.name != '직사각형 3':
        remove(x)
tb(s, 0.60, 0.90, 8.8, 2.35,
   [{'t': 'Conclusions', 'size': 16, 'bold': True, 'color': BLUE, 'after': 4},
    {'t': 'Higher genetically proxied blood *TSHR* expression was associated with lower odds of GD, hyperthyroidism and '
          'ophthalmopathy in all five datasets', 'bullet': True},
    {'t': 'Colocalization supported a shared *TSHR* variant in FinnGen GO and, prior-sensitively, in BBJ and GO strict, but favored '
          'distinct variants in UKB and FinnGen GD; expression-mediated causality is not established', 'bullet': True},
    {'t': '*IGF1R* ORs exceeded 1 but were nominal only in BBJ and UKB; FinnGen GD lay near the null. This does not argue against '
          'IGF-1R blockade', 'bullet': True},
    {'t': 'Susceptibility, not treatment response; GO analyses do not estimate TED risk among patients with GD', 'bullet': True}],
   'Conclusions', size=13, after=3)
rect(s, 0.60, 3.45, 8.80, 0.0104, RULE, name='Rule')
tb(s, 0.60, 3.52, 8.8, 1.2,
   [{'t': 'Limitations', 'size': 15, 'bold': True, 'color': TEXT, 'after': 3},
    {'t': 'Blood is neither thymus nor orbit; European eQTL effects applied to East Asian BBJ; one *TSHR* instrument', 'bullet': True},
    {'t': 'Three overlapping, population-controlled FinnGen endpoints; post hoc extensions and colocalization emphasis; '
          'single-causal-variant colocalization',
     'bullet': True}], 'Limitations', size=12.5, after=2)
rect(s, 0.60, 4.78, 8.80, 0.45, BOX, BOXLINE, name='Next box')
tb(s, 0.78, 4.78, 8.5, 0.45, ['**Next:** compare GD patients with and without TED; thymic and orbital regulatory data; '
                              'multi-signal colocalization'], 'Next', anchor='m', size=12.5, color=NAVY)
notes(s, """
결론입니다.
첫째, 유전적으로 예측된 혈액 TSHR 발현이 높을수록 다섯 자료 모두에서 GD, 갑상선기능항진증, 안병증의 odds가 낮았습니다.
둘째, 공위치 분석은 FinnGen GO에서 공유 변이를 지지했고, BBJ와 strict GO에서는 사전확률에 민감하게 지지했지만, UK Biobank와 FinnGen GD에서는 서로 다른 변이를 지지했습니다. 그래서 발현을 통한 인과성은 증명되지 않았습니다.
셋째, IGF1R은 모든 OR이 1보다 컸지만 nominal 연관은 BBJ와 UKB뿐이었고, FinnGen GD는 1에 가까웠습니다. 이것이 IGF-1R 차단 치료에 반대하는 근거는 아닙니다.
넷째, 이 연구는 치료 반응이 아니라 발병 감수성을 다루었고, 안병증 분석은 일반 인구 대조군을 써서 GD 환자 중 TED 감수성을 추정하지 못합니다.
제한점으로는 혈액이 흉선이나 안와 조직이 아니라는 점, 유럽인 eQTL을 동아시아인 BBJ에 적용한 점, TSHR 도구 변이가 하나뿐인 점, FinnGen 세 endpoint가 겹치고 일반 인구 대조군을 썼다는 점, FinnGen 확장과 공위치 분석의 강조가 post hoc이라는 점, 그리고 원인 변이를 하나로 가정한 공위치 모형이라는 점이 있습니다.
앞으로는 GD 환자 중 TED가 있는 사람과 없는 사람을 직접 비교하고, 흉선과 안와 조직의 조절 자료, 여러 신호를 허용하는 공위치 분석으로 확인해야 합니다. 감사합니다.
""")

# ================================================================ slide 16: thank you (original 17)
s = orig[16]
for x in s.shapes:
    if x.shape_type == 13:
        x.left, x.top, x.width, x.height = Inches(1.4), Inches(1.25), Inches(7.2), Inches(3.6)
        x._element.nvPicPr.cNvPr.set('descr', 'Cartoon of a thyroid, a chromosome detective and an eye linked by DNA')
    elif x.has_text_frame and x.text_frame.text.startswith('Oculoplastic'):
        x.left, x.top, x.width, x.height = Inches(1.0), Inches(4.86), Inches(8.0), Inches(0.36)
        retext(x, ['Oculoplastic Team, Seoul St. Mary’s Hospital'], size=18, bold=True, color=NAVY, align='c')
tb(s, 0.6, 5.24, 8.8, 0.3, ['Acknowledgement: Kyung-Hwa Shin (Pusan National University Hospital) provided the orbital '
                            'RNA-sequencing data.'], 'Acknowledgement', size=9.5, color=GREY, align='c')
notes(s, """
감사합니다. 안와 조직 RNA-sequencing 자료를 제공해 주신 부산대학교병원 진단검사의학과 신경화 선생님께 감사드립니다.
[발표 전 확인] 감사의 글 문구와 자료 사용 허락은 저자 확인 후 확정합니다.
""")

# ================================================================ backup divider and slides (new)
sd = prs.slides.add_slide(LAYOUT)
rect(sd, 0, 0, 10, 0.06, NAVY, name='Top bar')
tb(sd, 0.6, 2.2, 8.8, 0.7, ['Supplementary slides'], 'Divider title', size=32, bold=True, color=NAVY, align='c')
tb(sd, 0.6, 2.95, 8.8, 0.5, ['Colocalization across 24 genes · Leave-one-out · Tissue extension · Orbital RNA-seq'],
   'Divider list', size=15, color=GREY, align='c')
notes(sd, '질의응답용 보충 슬라이드입니다.')

b1 = new_slide('Supplement: Colocalization across 24 genes')
pic(b1, FIG / 'slide_fig3_heatmap.png', 0.55, 0.85, 'Figure 3A', 'Heatmap of shared-variant probability (PP.H4) for 24 genes '
    'in five datasets', w=8.9)
tb(b1, 0.60, 4.62, 8.8, 0.6,
   [{'t': 'Apart from *TSHR* and *CTLA4*, only *TNFSF14* and *IFNGR1* reached PP.H4 ≥ 0.80 (BBJ only); none in FinnGen GD',
     'bullet': True},
    {'t': 'Shading: PP.H4 (white 0 to blue 1); blank: < 0.001; box: ≥ 0.80 at p12 = 10^−5^; dot: still ≥ 0.80 at p12 = 10^−6^',
     'bullet': True}],
   'Heatmap notes', size=12, after=2)
source(b1, '*HLA-G* (right of divider): added because its MR result changed with the allele-frequency source. Values are posterior '
            'probabilities, not *P* values. Source: Figure 3A.')
notes(b1, """
24개 유전자의 공유 변이 확률(PP.H4)입니다. 테두리는 0.80 이상, 점은 사전확률을 10⁻⁶로 낮춰도 0.80 이상인 칸입니다.
TSHR과 CTLA4 외에는 TNFSF14와 IFNGR1만 BBJ에서 0.80을 넘었고, FinnGen GD에서는 어떤 유전자도 넘지 못했습니다.
24개 유전자 = BBJ 발견 13개 + IGF1R + FinnGen GD screen 추가 발견 9개 + HLA-G(대립유전자 빈도 출처에 따라 MR 결과가 바뀌어 추가). 이 24개만 공위치 분석을 했습니다.
""")

b2 = new_slide('Supplement: Single-instrument omission (leave-one-out)')
pic(b2, FIG / 'slide_figS2_loo.png', 0.55, 0.82, 'Figure S2', 'Leave-one-out MR estimates for IGF1R and CTLA4 in each dataset', w=8.9)
source(b2, 'Diamonds: all instruments; circles: one instrument omitted. Bold *P* values: nominal *P* < 0.05. BBJ *CTLA4* has one '
           'instrument. Leave-one-out reuses the same data; not an independent replication. Source: Figure S2.', y=5.20)
notes(b2, """
도구 변이를 하나씩 뺀 민감도 분석입니다.
IGF1R은 BBJ와 UKB에서 네 번 중 두 번 P가 0.05를 넘었고, UKB에서 rs2654980을 빼면 OR이 1.35에서 1.16으로 낮아졌습니다.
CTLA4는 rs13030124를 빼면 UKB와 FinnGen 세 endpoint 모두 유의성이 사라졌습니다(P = 0.323–0.972).
같은 자료를 다시 쓰는 분석이므로 독립적인 재현은 아닙니다.
""")

b3 = new_slide('Supplement: Tissue extension (GTEx v8)')
subhead(b3, '20 of 50 planned comparisons estimable; one instrument per gene–tissue combination')
rows = [['Gene / tissue', 'Instrument (F)', 'Result across five datasets'],
        ['*TSHR* / subcutaneous adipose', 'rs12101255 (38.5)', 'All ORs < 1 (0.31–0.51); 5 of 5 met family *P* threshold'],
        ['*TSHR* / visceral adipose', 'rs74064796 (32.1)', 'All ORs < 1 (0.39–0.78); 3 of 5 met family *P* threshold'],
        ['*IGF1R* / thyroid', 'rs1815009 (61.1)', 'No *P* < 0.05; primary FinnGen GD OR 0.76, other four above 1'],
        ['*IGF1R* / cultured skin fibroblasts', 'rs6598538 (33.8)', 'No *P* < 0.05; directions inconsistent'],
        ['Six other combinations', 'None', '30 comparisons not estimable (not negative results)']]
table(b3, 0.60, 1.30, [2.65, 1.75, 4.4], rows, 'Tissue table', size=12, hsize=12, row_h=0.4, head_h=0.36,
      align=['l', 'c', 'l'])
rect(b3, 0.60, 3.85, 8.80, 1.3, CALL, CALLLINE, name='Adipose box')
tb(b3, 0.78, 3.93, 8.45, 1.15,
   [{'t': 'The adipose instruments re-express the regional *TSHR* association', 'bold': True, 'color': NAVY, 'after': 3},
    {'t': 'rs12101255 and rs74064796 lie 50 bp apart in *TSHR* intron 1, about 15 kb from rs179252 (European r² 0.46 for rs12101255), and are '
          'blood *TSHR* eQTLs (*P* ≈ 10^−34^)', 'bullet': True},
    {'t': 'Same disease data and locus: no independent tissue confirmation; tissue colocalization not performed', 'bullet': True}],
   'Adipose', size=12.5, after=3)
source(b3, 'Primary family: 10 FinnGen GD comparisons (*P* < 0.005); secondary: 40 (*P* < 0.00125). GTEx adipose and fibroblasts '
            'are not orbital tissue. Source: Table S4.')
notes(b3, """
GTEx 조직 자료로 확장한 분석입니다. 계획한 50개 비교 중 20개만 추정 가능했고, 유전자–조직 조합마다 도구 변이가 하나뿐이었습니다.
지방조직 TSHR은 OR이 모두 1보다 낮았지만, 이 도구 변이들은 rs179252에서 약 15 kb 떨어진 TSHR intron 1에 있고 혈액 TSHR eQTL이기도 해서, 같은 지역의 연관을 다시 표현한 것이지 독립적인 조직 근거가 아닙니다.
갑상선과 배양 피부 섬유모세포(안와 조직 아님)의 IGF1R은 P < 0.05가 없었고 방향도 일관되지 않았습니다. 추정할 수 없었던 30개 비교는 음성 결과가 아닙니다.
""")

b4 = new_slide('Supplement: Descriptive orbital RNA-seq')
pic(b4, S1, 0.9, 0.9, 'Figure S1', 'Orbital transcript abundance (TPM) of TSHR, IGF1R and CTLA4 in four TED samples and one control', w=6.6)
rect(b4, 7.65, 0.95, 1.75, 3.0, PALE, RULE, name='Orbital note box')
tb(b4, 7.75, 1.05, 1.58, 2.85,
   [{'t': 'Descriptive only', 'bold': True, 'color': RED, 'size': 13, 'after': 4},
    {'t': '4 inactive TED (decompression) vs 1 control (nasal fat pad, different depot)', 'after': 4},
    {'t': 'No significance testing; not used to validate the genetic findings'}],
   'Orbital note', size=11.5)
tb(b4, 0.60, 4.25, 8.8, 0.9,
   [{'t': 'Mean orbital *TSHR*, *IGF1R* and *CTLA4* transcript levels in TED exceeded the single control; *IGF1R* was higher in '
          'only three of four TED samples', 'bullet': True},
    {'t': 'Data: Kyung-Hwa Shin, Pusan National University Hospital (IRB 2104-018-102)', 'bullet': True}],
   'Orbital findings', size=12.5, after=3)
source(b4, 'Points: individual samples; bars: TED means; axes differ by gene. TPM, transcripts per million; log~2~ fold change = '
            'log~2~[(mean TED TPM + 0.01)/(control TPM + 0.01)]. Source: Figure S1.')
notes(b4, """
안와 조직 RNA-sequencing 결과로, 비활동성 TED 환자 4명의 안와 지방과 대조군 1명의 코쪽 눈꺼풀 지방을 비교한 기술 통계입니다.
TED 평균 TSHR, IGF1R, CTLA4 발현이 대조군보다 높았지만, IGF1R은 네 명 중 세 명에서만 대조군보다 높았고, 대조군은 한 명이며 다른 부위의 지방입니다. 유의성 검정은 하지 않았고, 유전 분석 결과를 검증하는 자료로 쓰지 않았습니다.
""")

# ================================================================ order, delete, number
keep = [orig[0], orig[1], s3, orig[3], s5, orig[7], orig[10], s8, s9, s10, s11, s12, s13, s14, orig[15], orig[16],
        sd, b1, b2, b3, b4]
lst = prs.slides._sldIdLst
ids = {prs.part.related_part(x.rId): x for x in lst}
for sl in list(prs.slides):
    if sl not in keep:
        x = ids[sl.part]
        prs.part.drop_rel(x.rId)
        lst.remove(x)
for x in list(lst):
    lst.remove(x)
for sl in keep:
    lst.append(ids[sl.part])
for n, sl in enumerate(keep, 1):
    if n not in (1, 16, 17):
        slide_number(sl, n)
    el = sl._element
    timing = el.find(qn('p:timing'))
    if timing is not None:
        ids = {c.get('id') for c in el.iter(qn('p:cNvPr'))}
        if any(t.get('spid') not in ids for t in timing.iter(qn('p:spTgt'))):
            el.remove(timing)
HANGUL = re.compile('[\uac00-\ud7a3]')
for sl in (orig[0], orig[1], orig[3], orig[7], orig[10], orig[15], orig[16]):
    for x in sl.shapes:
        if x.name not in ORIG_NAMES[sl.slide_id]:
            continue
        if x.shape_type == 13:
            if not x._element.nvPicPr.cNvPr.get('descr'):
                x._element.nvPicPr.cNvPr.set('descr', 'The Catholic University of Korea logo')
            x.name = 'Picture'
        elif x.has_text_frame and x.text_frame.text.strip():
            t = x.text_frame.text.strip().split('\n')[0]
            x.name = 'Title' if abs(x.top - Inches(0.20)) < Inches(0.05) and x.left < Inches(1) else t[:40]
        elif x.name in ('Shape 0', '직사각형 10', '직사각형 12'):
            x.name = 'Top bar'
        elif x.name == '직사각형 3':
            x.name = 'Title rule'
        elif HANGUL.search(x.name) or x.name.startswith('Shape'):
            x.name = 'Connector' if x.shape_type == 9 else 'Box'
prs.core_properties.author = 'Jungyul Park; Min-Seon Kim; Suk-Woo Yang'
prs.core_properties.title = 'TSHR and IGF1R in Graves disease and thyroid eye disease'
prs.save(OUT)
print('saved', OUT, len(keep), 'slides')
