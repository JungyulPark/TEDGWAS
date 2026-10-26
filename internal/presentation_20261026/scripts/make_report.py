"""Change report: resident's original deck (v2, 17 slides) vs revised final deck (v4, 21 slides).

python3 make_report.py <thumb_dir> <out_dir>
Writes PRESENTATION_CHANGE_REPORT_20261026.html and .pdf (A4, Chromium via Playwright).
"""
import html
import re
import sys
from pathlib import Path

THUMBS, OUT = Path(sys.argv[1]), Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)
NAME = 'PRESENTATION_CHANGE_REPORT_20261026'


def md(s):
    """**bold**, *italic*, `code`, ^sup^ to HTML (text is escaped first)."""
    s = html.escape(s.replace('\\*', '\x00'), quote=False)
    s = re.sub(r'`([^`]+)`', r'<code>\1</code>', s)
    s = re.sub(r'\*\*(.+?)\*\*', r'<b>\1</b>', s)
    s = re.sub(r'(?<![\w*])\*([^*\s][^*]*?)\*(?![\w*])', r'<i>\1</i>', s)
    s = re.sub(r'\^([^^\s]+)\^', r'<sup>\1</sup>', s)
    return s.replace('\x00', '*')


def table(head, rows, cls='', widths=None):
    cg = ''.join(f'<col style="width:{w}">' for w in widths) if widths else ''
    th = ''.join(f'<th>{md(h)}</th>' for h in head)
    body = ''.join('<tr>' + ''.join(f'<td>{md(c)}</td>' for c in r) + '</tr>' for r in rows)
    return f'<table class="{cls}"><colgroup>{cg}</colgroup><thead><tr>{th}</tr></thead><tbody>{body}</tbody></table>'


def ul(items):
    return '<ul>' + ''.join(f'<li>{md(i)}</li>' for i in items) + '</ul>'


# ------------------------------------------------------------------ content
MAPPING = [
    ('1', '1', '수정', '제목을 원고 제목에 맞춤; 저자·소속·각주·원고 상태 수정'),
    ('2', '2', '수정', '근거 문장 4개로 재작성; 치료 서술 구체화; 출처 추가; 그림 축소와 캡션'),
    ('3', '3', '재작성', '선행연구(Ji, Zhang, Ke)와 연구 공백, 목적 3개, 핵심 질문'),
    ('4', '4', '수정', '노출(혈액 mRNA)·결과(GD, 갑상선기능항진증, GO) 표기; 공식을 텍스트로; 가정 문장'),
    ('5', '—', '삭제', '4번 복제 + 내용을 가리던 폭발 도형; 질문은 3번 Key question으로 옮김'),
    ('6', '5', '통합', 'eQTLGen/결과 자료 역할 표 → 데이터셋 표의 첫 행(노출)'),
    ('7', '5', '재작성', '결과 데이터셋 3개 → 5개; ancestry, FinnGen 중복, GO 대조군 주의'),
    ('8', '6', '수정', 'clumping 기준(r² < 0.001, 10 Mb), F 통계량 확인, 두 기준값; 상단 상자 여백'),
    ('9', '8', '재작성', '슬라이드용 새 Figure 1B; 그림을 가리던 설명 상자 제거; 요약 상자'),
    ('10', '8, 12', '통합', '유전자 범주 표 → 8번 요약 상자와 12번 표의 염색체 열'),
    ('11', '7', '수정', '오타 수정; 표를 가리던 설명 상자와 그 클릭 애니메이션 삭제; 표 재작성; 사전확률 상자'),
    ('12', '9, 10', '대체', '수정 전 Figure 2 → 새 Figure 2(5개 자료) + *TSHR* 숫자 표'),
    ('13', '12', '확대', 'MHC 5개 유전자(BBJ만) → BBJ 발견 13개 전체의 FinnGen GD 비교'),
    ('14', '13', '재작성', '방향 해석과 "연관 ≠ 공위치" 두 상자; 노트에서 지웠던 비공위치 내용 복원'),
    ('15', '11, 14', '분리', '표는 결과(11번, 5개 자료), 해석은 고찰(14번)'),
    ('16', '15', '재작성', '결과만 적은 결론 4개; 제한점 보완; 다음 단계'),
    ('17', '16', '수정', '감사의 글(신경화) 추가; 문구 위치 정리'),
    ('—', '17', '신규', '보충 슬라이드 구분'),
    ('—', '18', '신규', '24개 유전자 colocalization 히트맵(Figure 3A)'),
    ('—', '19', '신규', 'leave-one-out 민감도 분석(Figure S2)'),
    ('—', '20', '신규', '조직(GTEx v8) 확장 표(Table S4)'),
    ('—', '21', '신규', '안와 RNA-seq(Figure S1), 기술 통계만'),
]

MUST = [
    ('1. 저자·원고 상태',
     'Jungyul Park¹, Min-Sun Kim², **Kyung-Hwa Shin³\\***, Suk-Woo Yang¹\\*; "\\* Co-corresponding authors"; '
     '"Manuscript under review"',
     'Jungyul Park¹†, **Min-Seon** Kim²†, Suk-Woo Yang¹\\*; "† Equal contribution (co-first authors) \\* Corresponding author"; '
     '"Manuscript in preparation"; 감사 슬라이드에 신경화 선생님',
     '원고 저자란·Acknowledgements; 저자 결정(10월 19일); 아직 투고 전'),
    ('2. 결과 데이터셋 3 → 5',
     'BBJ, UKB, FinnGen GO 3개; UKB와 FinnGen 모두 "European"',
     'FinnGen GD(3,962 / 496,386)와 FinnGen GO strict(753 / 499,595) 추가; UKB는 "Not established as exclusively European", '
     'FinnGen은 Finnish; FinnGen 세 endpoint는 참여자가 겹쳐 독립 재현이 아님; GO 대조군에 안병증 기록이 없는 GD 포함',
     'Table 1; Methods'),
    ('3. *TSHR* colocalization',
     '"Strong colocalization support in BBJ and FinnGen"; FinnGen GD 없음; 14번 노트에서 비공위치·사전확률 내용을 "제거함"',
     'PP.H4 BBJ 0.951, FinnGen GO 0.986, GO strict 0.891 대 UKB 0.226, FinnGen GD 0.264(H3 우세); p12 = 10^−6^에서 BBJ 0.661, '
     'GO strict 0.451',
     'Table 2; Results; Abstract'),
    ('4. *IGF1R*',
     '표에 3개 자료(BBJ, UKB, FinnGen GO); "the FinnGen MR estimate was imprecise"',
     '5개 자료; FinnGen GD OR 1.07 (0.86–1.33), *P* = 0.539, PP.H4 0.010(가장 정밀, 1에 가까움, 검출 가능 OR 1.36); '
     'GO strict 1.52 (0.92–2.49), *P* = 0.0996; PP.H2가 5개 중 4개에서 최고; weighted median과 leave-one-out',
     'Table 2; Results; Figure S2'),
    ('5. *CTLA4*(노트)',
     '"BBJ 결합 기준은 통과하지 못했으며, Pvalue도 상당히 컸습니다"',
     'BBJ MR *P* = 5.45×10^−15^로 screen 통과, 기준에 못 미친 것은 PP.H4(0.201); UKB·GO의 *P*는 nominal이고 rs13030124에 의존'
     '(빼면 *P* = 0.323–0.972)',
     'Table 2; Table 3; Results'),
    ('6. 그림',
     '9번·12번에 수정 전 Figure 1·2(3개 자료, 그림 안 숫자, 옛 H4 막대); 설명 상자가 그림과 표를 가림',
     '워크북에서 슬라이드 크기로 다시 그림(Figure 1B, 2, 3A, S2); 가리는 도형 없음',
     'Supplementary_Dataset.xlsx; build_figures_20261014.py'),
    ('7. 치료 관련 서술',
     '노트 "TSHR 억제제 또한 갑상선안병증 치료에 효과가 있다고 하죠"; 슬라이드 "improves treatment outcomes"',
     '근거 없는 TSHR 억제제 문장 삭제; "reduced proptosis in randomized trials in active TED"; 그림 출처 기입 자리',
     'Introduction [9, 10]; Discussion [10, 44]'),
    ('8. 인과 표현',
     '"Assess the potential causal effects"; 노트 "인과관계를 … 나타내고자", "유전자 발현량이 질병의 원인일 가능성이 높다", '
     '"TSHR 발현량이 많을수록 GD가 발생할 확률은 낮아지는"; 노출을 "단백질"로 설명',
     '연관 표현으로 통일; "These data establish association, not expression-mediated causality"; 노출은 혈액 mRNA',
     'CLAUDE.md 해석 원칙; Discussion'),
    ('9. OR 0.12의 의미',
     '설명 없음(청중이 "88% 감소"로 오해하기 쉬움)',
     'rs179252 G allele: 혈액 *TSHR* 0.105 SD 증가(F = 176), allele당 OR 0.77–0.82; 1 SD당 OR은 이 중간 크기 효과의 환산',
     'Table 2 주석; Discussion'),
    ('10. 선행연구',
     '없음',
     'Ji 2025(FinnGen R10, strict ophthalmopathy), Zhang 2026(FinnGen R12, 같은 *TSHR* 패턴), Ke 2025(단백체 MR); 공백: *IGF1R* 비교 없음',
     'Introduction [18–20]'),
    ('11. 결론·제한점',
     '"Stronger genetic evidence for TSHR in BBJ GD and FinnGen ophthalmopathy"; 제한점 3줄',
     '결과만 4개: *TSHR* 역연관(5개 자료), colocalization(지지와 비지지 모두), *IGF1R*(nominal은 BBJ·UKB만), 감수성 ≠ 치료 반응; '
     '제한점: 혈액 ≠ 흉선·안와, 유럽 eQTL을 동아시아 BBJ에 적용, 단일 *TSHR* 도구변수, FinnGen 중복, post hoc 확장, 단일 원인변이 모형',
     'Conclusions; Limitations'),
]

NUMBERS = [
    ('결과 데이터셋 수', '3', '5', 'Table 1'),
    ('FinnGen GD 환자/대조군', '없음', '3,962 / 496,386', 'Table 1'),
    ('FinnGen GO strict 환자/대조군', '없음', '753 / 499,595', 'Table 1'),
    ('UKB ancestry', 'European', 'Not established as exclusively European', 'Table 1'),
    ('*TSHR* FinnGen GD OR (95% CI), *P*', '없음', '0.08 (0.05–0.12), 3.19×10^−32^', 'Table 2'),
    ('*TSHR* PP.H4, FinnGen GD / UKB', 'FinnGen GD 없음; UKB 0.226은 노트에만', '0.264 / 0.226 (둘 다 H3 우세)', 'Table 2'),
    ('*TSHR* PP.H4, p12 = 10^−6^ (BBJ / GO strict)', '노트에서 삭제', '0.661 / 0.451', 'Table 2'),
    ('*TSHR* allele당 OR', '없음', '0.77–0.82', 'Table 2 주석; Table S3'),
    ('*IGF1R* FinnGen GD', '없음', '1.07 (0.86–1.33), *P* = 0.539, PP.H4 0.010', 'Table 2'),
    ('*IGF1R* FinnGen GO strict', '없음', '1.52 (0.92–2.49), *P* = 0.0996, PP.H4 0.038', 'Table 2'),
    ('*IGF1R* FinnGen GD 검출 가능 OR(80% 검정력)', '없음', '1.36', 'Results'),
    ('*IGF1R* weighted median (BBJ / UKB)', '없음', '1.56 (*P* = 0.0288) / 1.42 (*P* = 0.00244)', 'Results'),
    ('*CTLA4* BBJ *P*', '"P value도 상당히 컸다"(노트)', '5.45×10^−15^ (screen 통과)', 'Table 2'),
    ('MHC 유전자 최대 PP.H4', '0.030 (*C4A*, 3개 자료)', '5개 자료로 보면 0.326 (*C4A*, GO strict); 모두 H3 우세', 'Colocalization 시트'),
    ('LD clumping 기준', '"approximately independent"', 'r² < 0.001, 10 Mb 창', 'Methods'),
    ('도구변수 강도', '없음', '모두 F > 10 (최소 29.7)', 'Methods'),
    ('FinnGen GD 확인 기준', '없음', '*P* < 0.05/13 = 0.00385, 같은 방향 (13개 중 7개 충족)', 'Methods; Table 3'),
    ('검정력', '없음', '80% 검정력으로 OR 1.5를 찾을 수 있는 유전자 14.6%', 'Results'),
]

KEPT = ('그대로 맞아 유지한 숫자: *TSHR* BBJ OR 0.12 (0.07–0.21), *P* = 1.09×10^−14^; *IGF1R* BBJ 1.56 (1.07–2.28), *P* = 0.0212, '
        'UKB 1.35 (1.07–1.70), *P* = 0.0117, FinnGen GO 1.41 (0.85–2.33), *P* = 0.182; MHC 5개 유전자의 BBJ OR과 *P*; '
        'BBJ·UKB·FinnGen GO 환자/대조군 수; 4,462 / 2,544 / 6,135 / 2,234 / 13; screen 기준 1.965×10^−5^.')

TABLES = [
    ('5', '데이터셋 표', '3행 → 노출 1행 + 결과 5행; endpoint 코드와 슬라이드 표기, ancestry, 역할; 셀 줄바꿈 정리'),
    ('7', '공위치 가설 표', '덱 표 서식으로 다시 만들고 H0 문장 정리; 가리던 설명 상자 대신 옆에 설정 상자'),
    ('10', '*TSHR* 표(신규)', '5개 자료 × OR (95% CI), *P*, PP.H3, PP.H4, PP.H4 (p12 = 10^−6^); 0.80 이상은 파란색, H3 우세는 주황색'),
    ('11', '*IGF1R* 표', '3행 → 5행; IVs와 PP.H2 열 추가; FinnGen GD 강조'),
    ('12', '13개 발견 유전자 표', 'MHC 5개(BBJ만) → 13개 × BBJ와 FinnGen GD OR·*P*, focused support, FinnGen GD PP.H4'),
    ('20', '조직 확장 표(신규)', '유전자–조직 조합별 도구변수와 5개 자료 결과; 30개 추정 불가 비교는 음성 결과가 아님'),
]

FIGURES = [
    ('8', 'Figure 1B (Manhattan)', '슬라이드 크기(6.0 × 4.3 in)로 다시 그림; 글자 10.5–12 pt; *TSHR*·*IGF1R*·*CTLA4* 굵게; 라벨 겹침 정리'),
    ('9', 'Figure 2 (MR + colocalization)', '5개 자료 × 3개 유전자; 그림 안 숫자 없음(표가 담당); H4를 0부터 쌓아 점선이 PP.H4 = 0.80'),
    ('18', 'Figure 3A (히트맵, 신규)', '가로형 5 × 24; 0.001 미만은 빈칸; 0.80 이상은 테두리; p12 = 10^−6^에서도 0.80 이상이면 점'),
    ('19', 'Figure S2 (leave-one-out, 신규)', '두 패널 가로 배치, 글자 10 pt; 굵은 *P*는 nominal *P* < 0.05'),
    ('21', 'Figure S1 (안와 RNA-seq)', '원고 그림 그대로 사용(개별 표본 자료가 있어야 다시 그릴 수 있음)'),
    ('4', 'MR 공식', '흐린 그림 → 슬라이드 텍스트(*OR* = exp(*β*disease / *β*expression))'),
]

FORMAT = [
    '유전자 기호와 *P*를 이탤릭으로 바꿈(원본은 이탤릭 0곳); 단백질 IGF-1R은 정체 유지',
    '쪽번호: 손으로 입력한 번호(3, 3, 3, 6, 5 … 11, 16, 17, 12) → 자동 번호',
    '내용을 가리던 도형 삭제: 5번 폭발 도형, 9번 설명 상자, 11번 설명 상자와 그 클릭 애니메이션, 14번 소제목 겹침',
    '13–15번의 맑은 고딕 → Calibri로 통일; 소제목 위치·크기 통일(위에서 0.92 in, 16 pt)',
    '6번 상단 두 상자를 본문 여백 안으로 정리; 2번 그림을 오른쪽 여백 안으로 축소',
    '그림 대체 텍스트(alt text) 추가; 옛 저자·인과 문구가 남아 있던 도형 이름 정리',
    '전공의 디자인(남색 상단 바, 회색 제목선, 파란 소제목, Calibri, 표 머리 색)은 그대로 유지',
]

NOTES = [
    '21장 모두 원고 기준으로 다시 씀(발표자 김민선, 말하는 순서대로)',
    '원본 노트의 오류 수정: TSHR 억제제 효과, *CTLA4* 설명, 인과 표현, "단백질" → 혈액 mRNA, "태어날 때" → 수정(잉태) 시점, '
    'log OR을 지수로 바꾼 값이 OR',
    '오타 수정: IGFR1, FinGen, 잇습니다, "맨델리안 무작위와", "그 안에서 bon"',
    '내부 파일명 삭제(Stage2_14_genes_primary.csv 등)',
    '질의응답 대비 내용 추가: 사전확률 민감도, *CTLA4*의 반대 패턴, MHC의 복잡한 LD, *HLA-G*를 넣은 이유, leave-one-out 해석',
    '"[발표 전 확인]" 메모 2곳: 2번 그림 출처, 감사 문구와 자료 사용 허락',
]

VERIFY = [
    '**1차 검토:** 원본 17장을 최종 원고·표·워크북·CLAUDE.md 고정 문구와 대조해 필수 수정 11개와 세부 항목을 정리함(review_checklist.md)',
    '**2차 독립 검증(v3 대상):** 숫자, 표현, 누락 항목, 그림, 레이아웃(2명) 검토 6개를 따로 돌리고, 검토자마다 반박 검증을 한 번 더 거침. '
    '70건 확인(중복 포함), 5건 기각. 모두 반영해 v4를 만듦',
    '확인된 지적의 성격: 표현의 정확성("the two largest datasets" → "환자 수가 가장 많은 두 자료", 14.6% 검정력 문장, '
    'r² 0.46은 rs12101255에만 해당, F > 10은 선택 기준이 아닌 확인), 빠진 단서(GO 결과는 주로 GD 감수성을 반영할 수 있음, '
    'colocalization 강조는 post hoc, MHC의 복잡한 LD), 레이아웃(여백, 줄바꿈, 소제목 위치), 삭제된 도형을 가리키던 클릭 애니메이션. '
    '**핵심 숫자(OR, CI, *P*, PP.H4)의 오류는 없었음**',
    '**형식 검사:** PowerPoint 구조 검사(validate.py) 통과; LibreOffice로 21장을 이미지로 렌더링해 확인; '
    '원고 저장소 audit PASS(357개 파일, 61개 문구 점검)',
]

OPEN = [
    '**2번 기전 그림 출처:** 슬라이드의 "Figure: add source" 자리에 출처를 넣어야 함(이 보고서에서는 그림을 가림)',
    '**감사 문구와 안와 자료 사용 허락:** 저자 확인 후 확정(PNUH IRB 2104-018-102)',
    '**원고 상태:** 발표 시점에 투고했다면 표지의 "Manuscript in preparation"을 고침',
    '**발표자 소속:** 노트 첫 줄은 "가톨릭대학교 안과 김민선입니다". 원고 소속 2에 서울성모병원이 없어서이며, 본인 소속에 맞게 고치면 됨',
    '**문서 메타데이터:** 원본에 있던 병원 문서보안(Fasoo) 추적 정보가 최종본에도 그대로 남아 있음. 외부에 공유하기 전에 기관 정책 확인',
    '**21번 Figure S1:** 원고 그림 그대로라 글꼴이 다른 슬라이드와 다름',
    '**발표 시간:** 본문 16장 기준. 시간이 정해지면 분량 조정',
]

# rows for the visual comparison: (original slides, final slide, note)
VISUAL = [
    ([1], 1, '저자·소속·각주·원고 상태'), ([2], 2, '근거 문장, 치료 서술, 출처'), ([3], 3, '선행연구와 목적'),
    ([4, 5], 4, 'MR 개념(5번 복제 슬라이드는 삭제)'), ([6, 7], 5, '데이터셋 5개'), ([8], 6, '도구변수와 기준값'),
    ([11], 7, 'colocalization 방법'), ([9, 10], 8, 'BBJ screen'), ([12], 9, '3개 유전자 개요(새 Figure 2)'),
    ([12, 14], 10, '*TSHR* 표'), ([15], 11, '*IGF1R* 표'), ([13, 10], 12, '13개 발견 유전자'),
    ([14], 13, '고찰: *TSHR*'), ([15], 14, '고찰: *IGF1R*와 치료'), ([16], 15, '결론·제한점'), ([17], 16, '감사'),
    ([], 17, '보충 슬라이드 구분'), ([], 18, '24개 유전자 colocalization'), ([], 19, 'leave-one-out'),
    ([], 20, '조직 확장'), ([], 21, '안와 RNA-seq'),
]


def thumb(kind, n):
    return f'<img src="{(THUMBS / f"{kind}-{n:02d}.jpg").as_uri()}" alt="{kind} slide {n}">'


vis_rows = []
for olist, f, note in VISUAL:
    left = ''.join(f'<figure>{thumb("orig", o)}<figcaption>원본 {o}</figcaption></figure>' for o in olist) or \
        '<div class="none">원본 없음(신규)</div>'
    vis_rows.append(f'<div class="vrow"><div class="vcol">{left}</div><div class="arrow">→</div>'
                    f'<div class="vcol"><figure>{thumb("final", f)}<figcaption>최종 {f} · {md(note)}</figcaption></figure></div></div>')

CSS = """
@page { size: A4; margin: 16mm 15mm 16mm 15mm; }
:root { --navy:#1E2A5E; --blue:#2266AA; --ink:#1A1A1A; --muted:#5B6472; --rule:#D5DAE2; --pale:#F4F6FA; --warm:#FFF8E8; }
* { box-sizing: border-box; }
body { font-family: 'Noto Sans KR', 'Carlito', sans-serif; color: var(--ink); font-size: 9.6pt; line-height: 1.55; margin: 0; }
h1 { font-size: 20pt; color: var(--navy); margin: 0 0 2mm; line-height: 1.25; letter-spacing: -0.2pt; }
.sub { font-size: 11pt; color: var(--blue); margin: 0 0 4mm; font-weight: 500; }
.meta { font-size: 8.4pt; color: var(--muted); border-top: 1.2pt solid var(--navy); border-bottom: 0.6pt solid var(--rule);
        padding: 2.2mm 0; margin-bottom: 5mm; }
.meta div { margin: 0.4mm 0; }
h2 { font-size: 12.5pt; color: var(--navy); margin: 7mm 0 2.5mm; padding-bottom: 1mm; border-bottom: 0.6pt solid var(--rule);
     break-after: avoid; }
h3 { font-size: 10.5pt; color: var(--blue); margin: 4mm 0 1.5mm; break-after: avoid; }
p { margin: 0 0 2mm; }
ul { margin: 0 0 2mm; padding-left: 5mm; }
li { margin: 0 0 1.2mm; }
code { font-family: 'DejaVu Sans Mono', 'Noto Sans KR', monospace; font-size: 8pt; background: var(--pale); padding: 0 1mm; border-radius: 1mm; }
table { width: 100%; border-collapse: collapse; font-size: 8.4pt; line-height: 1.45; margin: 1mm 0 3mm; table-layout: fixed; }
th { background: var(--navy); color: #fff; font-weight: 700; text-align: left; padding: 1.4mm 1.8mm; }
td { border-bottom: 0.5pt solid var(--rule); padding: 1.4mm 1.8mm; vertical-align: top; overflow-wrap: anywhere; }
tr { break-inside: avoid; }
tbody tr:nth-child(even) td { background: #FAFBFD; }
.tiles { display: grid; grid-template-columns: repeat(4, 1fr); gap: 3mm; margin: 2mm 0 4mm; }
.tile { background: var(--pale); border: 0.6pt solid var(--rule); border-radius: 1.5mm; padding: 2.5mm 3mm; }
.tile b { display: block; font-size: 17pt; color: var(--navy); line-height: 1.15; font-variant-numeric: tabular-nums; }
.tile span { font-size: 8pt; color: var(--muted); }
.box { background: var(--warm); border: 0.6pt solid #D9B55A; border-radius: 1.5mm; padding: 2.5mm 3.5mm; margin: 2mm 0 3mm; }
.small { font-size: 8.2pt; color: var(--muted); }
.vrow { display: grid; grid-template-columns: 1fr 6mm 1fr; align-items: center; gap: 2mm; padding: 2mm 0;
        border-bottom: 0.5pt solid var(--rule); break-inside: avoid; }
.vcol { display: flex; flex-direction: column; gap: 1.5mm; }
.vcol figure { margin: 0; }
.vcol img { width: 100%; border: 0.5pt solid #B8C0CC; display: block; }
.vcol figcaption { font-size: 7.6pt; color: var(--muted); margin-top: 0.6mm; }
.arrow { text-align: center; color: var(--blue); font-size: 13pt; }
.none { font-size: 8pt; color: var(--muted); text-align: center; border: 0.6pt dashed #B8C0CC; padding: 8mm 0; }
.pb { break-before: page; }
"""

body = f"""
<h1>발표자료 수정 보고서</h1>
<p class="sub">전공의 원본(v2, 17장)과 수정 최종본(v4, 21장)의 차이</p>
<div class="meta">
<div><b>발표 주제:</b> <i>TSHR</i> and <i>IGF1R</i> in Graves disease and thyroid eye disease (Mendelian randomization and colocalization)</div>
<div><b>원본:</b> <code>6f53284e-164749_2_______F_011__________v2.pptx</code> (파일명의 한글 부분은 업로드 과정에서 깨짐) · 17장</div>
<div><b>최종본:</b> <code>TSHR_IGF1R_MR_발표_최종.pptx</code>와 PDF · 21장(본문 16 + 보충 5), 한국어 발표자 노트</div>
<div><b>기준 자료:</b> 최종 원고 <code>MANUSCRIPT_TED_TRAP_v5_MASTER.md</code>(commit f5efd35), Tables 1–3, Tables S1–S4,
<code>Supplementary_Dataset.xlsx</code>, <code>CLAUDE.md</code>의 해석 원칙과 고정 문구</div>
<div><b>작성:</b> 2026년 10월 26일</div>
</div>

<h2>1. 한눈에 보기</h2>
<div class="tiles">
<div class="tile"><b>17 → 21</b><span>슬라이드(본문 16 + 보충 5)</span></div>
<div class="tile"><b>11</b><span>필수 수정 항목</span></div>
<div class="tile"><b>70</b><span>독립 검증에서 확인되어 반영한 지적(중복 포함)</span></div>
<div class="tile"><b>0</b><span>최종본 검증에서 나온 핵심 숫자 오류</span></div>
</div>
{ul([
    '원본은 10월 14일 이전 원고(결과 데이터셋 3개)를 바탕으로 만들어져, 화면의 숫자는 대체로 맞았지만 FinnGen GD와 strict 안병증, '
    '그리고 그에 달린 핵심 결과(*TSHR* 비공위치, *IGF1R* FinnGen GD의 null에 가까운 추정치)가 빠져 있었음. 저자 정보도 옛 판이었음',
    '원본 7장(1, 2, 4, 8, 11, 16, 17)은 디자인을 살려 수정했고, 나머지 10장은 통합·대체·삭제했으며, 14장을 새로 만듦',
    '표 6개를 새로 만들거나 다시 만들었고, 그림 4개를 워크북에서 슬라이드 크기로 다시 그림',
    '발표자 노트 21장을 모두 원고 기준으로 다시 씀',
    '전공의 디자인(색, 글꼴, 상단 바, 제목 위치)은 그대로 유지함',
])}

<h2>2. 슬라이드 대응표</h2>
{table(['원본', '최종', '처리', '주요 변경'], MAPPING, widths=['9%', '10%', '9%', '72%'])}

<h2>3. 필수 수정 11개</h2>
<p class="small">원본에서 무엇이 문제였고 최종본에서 어떻게 고쳤는지, 근거가 된 원고 위치와 함께 정리함.</p>
{table(['항목', '원본', '최종본', '근거'], MUST, widths=['13%', '32%', '40%', '15%'])}

<h2>4. 숫자 대조</h2>
<p class="small">새로 넣거나 고친 숫자. 모두 원고 표와 워크북에서 다시 확인함.</p>
{table(['항목', '원본', '최종본', '출처'], NUMBERS, widths=['30%', '22%', '33%', '15%'])}
<div class="box">{md(KEPT)}</div>

<h2>5. 표와 그림</h2>
<h3>표</h3>
{table(['최종', '표', '변경'], TABLES, widths=['8%', '22%', '70%'])}
<h3>그림</h3>
{table(['최종', '그림', '변경'], FIGURES, widths=['8%', '26%', '66%'])}
<p class="small">다시 그린 그림은 모두 <code>Supplementary_Dataset.xlsx</code>에서 원고 그림 스크립트
(<code>build_figures_20261014.py</code>)와 같은 데이터·논리로 그렸고, 글꼴(Carlito, Calibri와 같은 폭)과 크기만 슬라이드에 맞춤.</p>

<h2>6. 형식과 디자인</h2>
{ul(FORMAT)}

<h2>7. 발표자 노트</h2>
{ul(NOTES)}

<h2>8. 검증</h2>
{ul(VERIFY)}

<h2>9. 남은 확인 사항</h2>
{ul(OPEN)}

<h2 class="pb">부록. 슬라이드 비교(원본 → 최종)</h2>
<p class="small">왼쪽은 원본, 오른쪽은 최종본. 2번의 외부 기전 그림은 출처가 확인되지 않아 이 보고서에서는 회색 상자로 가림(실제 슬라이드에는 그대로 있음).</p>
{''.join(vis_rows)}
"""

page = f"""<!doctype html><html lang="ko"><head><meta charset="utf-8"><title>발표자료 수정 보고서</title>
<style>{CSS}</style></head><body>{body}</body></html>"""
(OUT / f'{NAME}.html').write_text(page, encoding='utf-8')

# ------------------------------------------------------------------ GitHub Markdown (same content, renders on GitHub)
def gfm(s):
    s = re.sub(r'\^([^^\s]+)\^', r'<sup>\1</sup>', s.replace('\\*', '\x00'))
    return s.replace('|', '\\|').replace('\x00', '\\*')


def gtable(head, rows):
    out = ['| ' + ' | '.join(gfm(h) for h in head) + ' |', '|' + '---|' * len(head)]
    return '\n'.join(out + ['| ' + ' | '.join(gfm(c) for c in r) + ' |' for r in rows])


def glist(items):
    return '\n'.join('- ' + gfm(i) for i in items)


def gthumb(kind, n):
    return f'<img src="thumbs/{kind}-{n:02d}.jpg" width="400" alt="{kind} slide {n}">'


vis_md = []
for olist, f, note in VISUAL:
    left = '<br>'.join(f'{gthumb("orig", o)}<br><sub>원본 {o}</sub>' for o in olist) or '원본 없음(신규)'
    vis_md.append(f'| {left} | {gthumb("final", f)}<br><sub>최종 {f} · {gfm(note)}</sub> |')

MD = f"""# 발표자료 수정 보고서

**전공의 원본(v2, 17장)과 수정 최종본(v4, 21장)의 차이**

[PDF로 보기]({NAME}.pdf)

- **발표 주제:** *TSHR* and *IGF1R* in Graves disease and thyroid eye disease (Mendelian randomization and colocalization)
- **원본:** `6f53284e-164749_2_______F_011__________v2.pptx` (파일명의 한글 부분은 업로드 과정에서 깨짐) · 17장
- **최종본:** `TSHR_IGF1R_MR_발표_최종.pptx`와 PDF · 21장(본문 16 + 보충 5), 한국어 발표자 노트
- **기준 자료:** 최종 원고 `MANUSCRIPT_TED_TRAP_v5_MASTER.md`(commit f5efd35), Tables 1–3, Tables S1–S4, `Supplementary_Dataset.xlsx`, `CLAUDE.md`의 해석 원칙과 고정 문구
- **작성:** 2026년 10월 26일
- 두 PPT 파일은 병원 문서보안 추적 정보가 들어 있어 이 공개 저장소에는 올리지 않음

## 1. 한눈에 보기

| 17 → 21 | 11 | 70 | 0 |
|---|---|---|---|
| 슬라이드(본문 16 + 보충 5) | 필수 수정 항목 | 독립 검증에서 확인되어 반영한 지적(중복 포함) | 최종본 검증에서 나온 핵심 숫자 오류 |

{glist([
    '원본은 10월 14일 이전 원고(결과 데이터셋 3개)를 바탕으로 만들어져, 화면의 숫자는 대체로 맞았지만 FinnGen GD와 strict 안병증, '
    '그리고 그에 달린 핵심 결과(*TSHR* 비공위치, *IGF1R* FinnGen GD의 null에 가까운 추정치)가 빠져 있었음. 저자 정보도 옛 판이었음',
    '원본 7장(1, 2, 4, 8, 11, 16, 17)은 디자인을 살려 수정했고, 나머지 10장은 통합·대체·삭제했으며, 14장을 새로 만듦',
    '표 6개를 새로 만들거나 다시 만들었고, 그림 4개를 워크북에서 슬라이드 크기로 다시 그림',
    '발표자 노트 21장을 모두 원고 기준으로 다시 씀',
    '전공의 디자인(색, 글꼴, 상단 바, 제목 위치)은 그대로 유지함',
])}

## 2. 슬라이드 대응표

{gtable(['원본', '최종', '처리', '주요 변경'], MAPPING)}

## 3. 필수 수정 11개

원본에서 무엇이 문제였고 최종본에서 어떻게 고쳤는지, 근거가 된 원고 위치와 함께 정리함.

{gtable(['항목', '원본', '최종본', '근거'], MUST)}

## 4. 숫자 대조

새로 넣거나 고친 숫자. 모두 원고 표와 워크북에서 다시 확인함.

{gtable(['항목', '원본', '최종본', '출처'], NUMBERS)}

> {gfm(KEPT)}

## 5. 표와 그림

### 표

{gtable(['최종', '표', '변경'], TABLES)}

### 그림

{gtable(['최종', '그림', '변경'], FIGURES)}

다시 그린 그림은 모두 `Supplementary_Dataset.xlsx`에서 원고 그림 스크립트(`build_figures_20261014.py`)와 같은 데이터·논리로 그렸고, 글꼴(Carlito, Calibri와 같은 폭)과 크기만 슬라이드에 맞춤. 슬라이드용 그림은 [`figures/`](figures/)에 있음.

## 6. 형식과 디자인

{glist(FORMAT)}

## 7. 발표자 노트

{glist(NOTES)}

## 8. 검증

{glist(VERIFY)}

검증 기록: [`review_checklist_20261025.md`](review_checklist_20261025.md), [`deck_review_findings_20261025.md`](deck_review_findings_20261025.md)

## 9. 남은 확인 사항

{glist(OPEN)}

## 부록. 슬라이드 비교(원본 → 최종)

왼쪽은 원본, 오른쪽은 최종본. 2번의 외부 기전 그림은 출처가 확인되지 않아 이 보고서에서는 회색 상자로 가림(실제 슬라이드에는 그대로 있음).

| 원본 | 최종 |
|---|---|
""" + '\n'.join(vis_md) + '\n'
(OUT / f'{NAME}.md').write_text(MD, encoding='utf-8')

import subprocess  # noqa: E402

subprocess.run(['node', str(Path(__file__).with_name('print_pdf.js')), str(OUT / f'{NAME}.html'), str(OUT / f'{NAME}.pdf')],
               check=True)
print('written', OUT / f'{NAME}.pdf')
