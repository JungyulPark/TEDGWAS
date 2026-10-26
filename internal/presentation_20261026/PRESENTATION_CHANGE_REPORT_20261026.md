# 발표자료 수정 보고서

**전공의 원본(v2, 17장)과 수정 최종본(v4, 21장)의 차이**

[PDF로 보기](PRESENTATION_CHANGE_REPORT_20261026.pdf)

- **발표 주제:** *TSHR* and *IGF1R* in Graves disease and thyroid eye disease (Mendelian randomization and colocalization)
- **원본:** `6f53284e-164749_2_______F_011__________v2.pptx` (파일명의 한글 부분은 업로드 과정에서 깨짐) · 17장
- **최종본:** `TSHR_IGF1R_MR_발표_최종.pptx`와 PDF · 21장(본문 16 + 보충 5), 한국어 발표자 노트
- **기준 자료:** 최종 원고 `MANUSCRIPT_TED_TRAP_v5_MASTER.md`(commit f5efd35), Tables 1–3, Tables S1–S4, `Supplementary_Dataset.xlsx`, `CLAUDE.md`의 해석 원칙과 고정 문구
- **버전:** v2 = 전공의 원본, v3 = 중간 수정본(독립 검증 대상), v4 = 최종본 · **작성:** 2026년 10월 26일
- 두 PPT 파일은 병원 문서보안 추적 정보가 들어 있어 이 공개 저장소에는 올리지 않음

## 1. 한눈에 보기

| 17 → 21 | 11 | 70 | 0 |
|---|---|---|---|
| 슬라이드(본문 16 + 보충 5) | 필수 수정 항목 | 중간 수정본(v3) 독립 검증에서 확인된 지적(중복 포함; 2건 외 반영) | 중간 수정본(v3) 독립 검증에서 나온 핵심 숫자(OR, CI, *P*, PP.H4) 오류 |

- 원본은 10월 14일 이전 원고(질환 데이터셋 3개)를 바탕으로 만들어져 화면의 숫자는 대체로 맞았지만, FinnGen GD와 FinnGen GO strict 두 자료와 여기서 나온 핵심 결과가 빠져 있었음: FinnGen GD에서 *TSHR* colocalization이 지지되지 않음(PP.H4 0.264, H3 우세), *IGF1R* FinnGen GD OR이 1에 가까움(1.07, *P* = 0.539). 저자 정보도 이전 버전이었음
- 원본 7장(1, 2, 4, 8, 11, 16, 17)은 원래 슬라이드를 고쳐 썼고(16번은 내용을 다시 씀), 나머지 10장의 내용은 새로 만든 본문 9장(최종 3, 5, 8–14번)으로 옮기거나 합치거나 나눴으며 5번은 삭제함. 보충 5장(최종 17–21번)을 더해 새로 만든 슬라이드는 14장. 장별 처리는 2절 대응표 참고
- 표 6개를 새로 만들거나 다시 만들었고, 그림 4개를 워크북에서 슬라이드 크기로 다시 그림
- 발표자 노트 21장을 모두 원고 기준으로 다시 씀
- 전공의 디자인(색, 글꼴, 상단 바, 제목 위치)은 그대로 유지함

## 2. 슬라이드 대응표

| 원본 | 최종 | 처리 | 주요 변경 |
|---|---|---|---|
| 1 | 1 | 수정 | 제목을 원고 제목에 맞춤; 저자·소속·각주·원고 상태 수정; "No financial interests to disclose" → "No conflicts of interest to disclose" |
| 2 | 2 | 수정 | 근거 문장 4개로 재작성; 치료 서술 구체화; 하단 참고문헌 추가; 그림 축소와 캡션 |
| 3 | 3 | 재작성 | 선행연구(Ji, Zhang, Ke)와 연구 공백, 목적 3개, 핵심 질문 |
| 4 | 4 | 수정 | 노출(혈액 mRNA)·질환(GD, 갑상선기능항진증, GO) 표기; 공식을 텍스트로; 가정 문장 |
| 5 | — | 삭제 | 4번 복제 + 내용을 가리던 폭발 도형. "*TSHR*이 *IGF1R*보다 더 강하게 연관되는가"라는 비교형 질문을 "두 유전자 발현이 GD·TED 관련 진단과 연관되는가, 같은 원인 변이를 공유하는가"로 바꿔 최종 3번 Key question에 넣음 |
| 6 | 5 | 통합 | eQTLGen/질환 자료 역할 표 → 데이터셋 표의 첫 행(노출); "31,684 participants" → "Up to 31,684" |
| 7 | 5 | 재작성 | 질환 데이터셋 3개 → 5개(FinnGen R12 GD·GO strict는 post hoc 표시); UKB 표현형 "Self-reported hyperthyroidism or thyrotoxicosis"; ancestry, FinnGen 참여자 중복, GO 대조군 주의; 소제목 "Summary statistics only" |
| 8 | 6 | 수정 | clumping 기준(r² < 0.001, 10 Mb), F 통계량 확인, 두 기준값; 상단 두 상자를 본문 여백 안으로 |
| 9 | 8 | 재작성 | 슬라이드용 새 Figure 1B; 그림을 가리던 Bonferroni 설명선과 클릭으로 나타나던 요약 상자 대신 그림 옆 요약 상자 |
| 10 | 8, 12 | 통합 | 유전자 범주 표 → 8번 요약 상자와 12번 표의 염색체 열 |
| 11 | 7 | 수정 | 오타 수정(casual → causal); 표를 가리던 설명선과 그 클릭 애니메이션 삭제; 표 재작성; 사전확률 상자; PP.H4 ≥ 0.80을 "강한 지지 기준"에서 "기술적 보고 기준"으로 바꾸고 단일 원인변이 가정과 post hoc 강조를 적음 |
| 12 | 9, 10 | 대체 | 수정 전 Figure 2 → 새 Figure 2(5개 자료) + *TSHR* 숫자 표 |
| 13 | 12 | 확대 | MHC 5개 유전자(BBJ만) → BBJ 발견 13개 전체의 FinnGen GD 비교 |
| 14 | 13, 10 | 재작성 | 방향 해석과 "연관 ≠ colocalization" 두 상자(13번); 원본 노트에서 "(제거함)"으로 표시됐던 UKB H3 우세는 13번에, 사전확률 민감도(BBJ 0.951 → 0.661)는 10번 표·노트에 복원 |
| 15 | 11, 14 | 분리 | 표는 결과(11번, 5개 자료), 해석은 고찰(14번). 14번에 추가: 한 유전자만 유의한 것은 차이 검정이 아님, 조직 효과 미해결(GTEx 단일 도구변수, 안와 eQTL 없음), 오래되고 활동성이 낮은 TED에서도 안구돌출 감소(Douglas 2024), TSHR–IGF-1R 상호작용은 감수성 문제가 아님 |
| 16 | 15 | 재작성 | 결과만 적은 결론 4개; 제한점 보완; 다음 단계 |
| 17 | 16 | 수정 | 감사의 글(신경화) 추가; 문구 위치 정리 |
| — | 17 | 신규 | 보충 슬라이드 구분 |
| — | 18 | 신규 | 24개 유전자 colocalization 히트맵(Figure 3A). 원본의 14개(BBJ 13 + *IGF1R*)에 FinnGen GD 전체 screen의 추가 발견 9개와 *HLA-G*를 더함 |
| — | 19 | 신규 | leave-one-out 민감도 분석(Figure S2) |
| — | 20 | 신규 | 조직(GTEx v8) 확장 표(Table S4) |
| — | 21 | 신규 | 안와 RNA-seq(Figure S1), 기술 통계만 |

## 3. 필수 수정 11개

원본에서 무엇이 문제였고 최종본에서 어떻게 고쳤는지, 근거가 된 원고 위치와 함께 정리함.

| 항목 | 원본 | 최종본 | 근거 |
|---|---|---|---|
| 1. 저자·원고 상태 | Jungyul Park¹, Min-Sun Kim², **Kyung-Hwa Shin³\***, Suk-Woo Yang¹\*; "\* Co-corresponding authors"; "Manuscript under review" | Jungyul Park¹†, **Min-Seon** Kim²†, Suk-Woo Yang¹\*; "† Equal contribution (co-first authors) \* Corresponding author"; "Manuscript in preparation"; 감사 슬라이드에 신경화 선생님 | 원고 저자란·감사의 글; 저자 결정(10월 19일); 아직 투고 전 |
| 2. 질환 데이터셋 3 → 5 | BBJ, UKB, FinnGen GO 3개; UKB와 FinnGen 모두 "European" | FinnGen R12 GD(3,962 / 496,386)와 GO strict(753 / 499,595)를 추가하고 post hoc으로 표시; UKB 표현형은 "Self-reported hyperthyroidism or thyrotoxicosis", ancestry는 "Not established as exclusively European"; FinnGen은 Finnish; FinnGen 세 endpoint는 참여자가 겹쳐 독립 재현이 아님; GO 대조군에 안병증 기록이 없는 GD 포함 | Table 1; Methods |
| 3. *TSHR* colocalization | "Strong colocalization support in BBJ and FinnGen under the default prior"(FinnGen은 GO만 해당); FinnGen GD 없음; 14번 노트에서 비공유 변이·사전확률 내용을 "(제거함)"으로 표시 | PP.H4 BBJ 0.951, FinnGen GO 0.986, GO strict 0.891 대 UKB 0.226, FinnGen GD 0.264(둘 다 H3 우세); p12 = 10<sup>−6</sup>에서 BBJ 0.661, GO strict 0.451 | Table 2; Results; Abstract |
| 4. *IGF1R* | 표에 3개 자료(BBJ, UKB, FinnGen GO); "the FinnGen MR estimate was imprecise" | 5개 자료; FinnGen GD OR 1.07 (0.86–1.33), *P* = 0.539(가장 정밀한 추정치로 1에 가까움; 80% 검정력의 검출 가능 OR 1.36), PP.H4 0.010; GO strict 1.52 (0.92–2.49), *P* = 0.0996; PP.H2가 5개 중 4개 자료에서 최고; weighted median과 leave-one-out | Table 2; Results; Figure S2 |
| 5. *CTLA4*(노트) | "BBJ 결합 기준은 통과하지 못했으며, Pvalue도 상당히 컸습니다" | BBJ MR *P* = 5.45×10<sup>−15</sup>로 screen 통과, 기준에 못 미친 것은 PP.H4(0.201); UKB와 두 GO 자료의 *P*는 nominal(0.0022–0.0118); UKB와 FinnGen 세 endpoint 추정치는 rs13030124에 의존(빼면 *P* = 0.323–0.972) | Table 2; Table 3; Results |
| 6. 그림 | 9번·12번에 수정 전 Figure 1·2(3개 자료, 그림 안 숫자, 옛 H4 막대); 설명 상자가 그림과 표를 가림 | 워크북에서 슬라이드 크기로 다시 그림(Figure 1B, 2, 3A, S2); 가리는 도형 없음 | 워크북(Supplementary Dataset); 원고 그림 스크립트 |
| 7. 치료 관련 서술 | 노트 "TSHR 억제제 또한 갑상선안병증 치료에 효과가 있다고 하죠"; 슬라이드 "improves treatment outcomes" | 근거 없는 TSHR 억제제 문장 삭제; "reduced proptosis in randomized trials in active TED"; 2번 기전 그림 아래에 출처 기입 자리("Figure: add source")만 만듦(출처는 미기입) | Introduction [9, 10]; Discussion [10, 44] |
| 8. 인과 표현 | "Assess the potential causal effects"; 노트 "인과관계를 … 나타내고자", "유전자 발현량이 질병의 원인일 가능성이 높다", "TSHR발현량이 많을수록 GD가 발생할 확률은 낮아지는"; 노출을 "단백질"로 설명 | 연관 표현으로 통일; 13번 소제목과 15번 결론에 "expression-mediated causality is not established"; 노출은 혈액 mRNA | CLAUDE.md 해석 원칙; Discussion |
| 9. OR 0.12의 의미 | 12번 캡션과 노트에 "약 1 SD당 값"이라는 설명만 있고, allele당 OR(0.77–0.82)과 그 환산 관계는 없음(청중이 "88% 감소"로 오해하기 쉬움) | rs179252 G allele: 혈액 *TSHR* 0.105 SD 증가(F = 176), allele당 OR 0.77–0.82; 1 SD당 OR은 이 중간 크기 효과를 환산한 값 | Table 2 주석; Discussion |
| 10. 선행연구 | *TSHR*·*IGF1R*에 대한 선행 MR·colocalization 연구(Ji, Zhang, Ke)와 연구 공백의 제시 없음(흉선 가설 등 다른 선행연구 언급만 있음) | Ji 2025(FinnGen R10, strict ophthalmopathy), Zhang 2026(FinnGen R12, 같은 *TSHR* 패턴), Ke 2025(단백체 MR); 연구 공백: *IGF1R*을 *TSHR*과 비교한 연구가 없고, 모두 주로 FinnGen에 의존함 | Introduction [18–20] |
| 11. 결론·제한점 | "Stronger genetic evidence for TSHR in BBJ GD and FinnGen ophthalmopathy"; 제한점 2줄과 다음 단계 1줄 | 결과만 4개: *TSHR* 역연관(5개 자료), colocalization(지지와 비지지 모두), *IGF1R*(nominal은 BBJ·UKB만), 감수성 ≠ 치료 반응(원본 제한점 "Population controls …"를 결론 4로 옮김); 제한점: 혈액 ≠ 흉선·안와, 유럽 eQTL을 동아시아 BBJ에 적용, 단일 *TSHR* 도구변수, FinnGen 중복, post hoc 확장, 단일 원인변이 모형 | Conclusions; Limitations |

## 4. 숫자 대조

새로 넣거나 고친 숫자 중 주요 항목(전체 목록은 아님). 모두 원고 표와 워크북에서 다시 확인함.

| 항목 | 원본 | 최종본 | 출처 |
|---|---|---|---|
| 질환 데이터셋 수 | 3 | 5 | Table 1 |
| FinnGen GD 환자/대조군 | 없음 | 3,962 / 496,386 | Table 1 |
| FinnGen GO strict 환자/대조군 | 없음 | 753 / 499,595 | Table 1 |
| UKB ancestry | European | Not established as exclusively European | Table 1 |
| *TSHR* FinnGen GD OR (95% CI), *P* | 없음 | 0.08 (0.05–0.12), 3.19×10<sup>−32</sup> | Table 2 |
| *TSHR* PP.H4, FinnGen GD / UKB | FinnGen GD 없음; UKB 0.226은 12번 그림(H4 열)과 노트에 있었으나 14번 고찰·16번 결론에는 빠짐 | 0.264 / 0.226 (둘 다 H3 우세) | Table 2 |
| *TSHR* PP.H4, p12 = 10<sup>−6</sup> (BBJ / GO strict) | BBJ 0.661은 14번 노트에 "(제거함)"으로 표시된 채 남아 있음(슬라이드에는 없음); GO strict는 원본에 없음 | 0.661 / 0.451 | Table 2 |
| *TSHR* allele당 OR | 없음 | 0.77–0.82 | Table 2 주석; Table S3 |
| *IGF1R* FinnGen GD | 없음 | 1.07 (0.86–1.33), *P* = 0.539, PP.H4 0.010 | Table 2 |
| *IGF1R* FinnGen GO strict | 없음 | 1.52 (0.92–2.49), *P* = 0.0996, PP.H4 0.038 | Table 2 |
| *IGF1R* FinnGen GD 검출 가능 OR(80% 검정력) | 없음 | 1.36 | Results |
| *IGF1R* weighted median (BBJ / UKB) | 없음 | 1.56 (*P* = 0.0288) / 1.42 (*P* = 0.00244) | Results |
| *CTLA4* BBJ *P* | 12번 그림에는 5.45×10<sup>−15</sup> (†)로 맞게 표시; 노트만 "Pvalue도 상당히 컸습니다"로 잘못 설명 | 5.45×10<sup>−15</sup> (screen 통과)를 노트에도 바르게 설명 | Table 2 |
| MHC 유전자 최대 PP.H4 | 0.030 (*C4A*, 3개 자료; 맞음) | 5개 자료로 보면 0.326 (*C4A*, GO strict); 모두 H3 우세 | Supplementary Dataset의 Colocalization 시트 |
| LD clumping 기준 | "approximately independent" | r² < 0.001, 10 Mb 창 | Methods |
| 도구변수 강도 | 없음 | 모두 F > 10 (최소 29.7) | Methods |
| FinnGen GD 확인 기준 | 없음 | *P* < 0.05/13 = 0.00385, 같은 방향 (13개 중 7개 충족) | Methods; Table 3 |
| 검정력 | 없음 | BBJ 추정 가능 유전자 2,234개 중 screen 기준에서 OR 1.5를 80% 검정력으로 검출할 수 있었던 유전자 14.6% | Results |

> 그대로 맞아 유지한 숫자: 원본 12번 그림의 *TSHR*·*IGF1R*·*CTLA4* OR (95% CI), *P*, PP.H4(BBJ·UKB·FinnGen GO; 모두 Table 2와 같음), 15번 *IGF1R* 표, 13번 MHC 5개 유전자의 BBJ OR·*P*와 최대 PP.H4 0.030, BBJ·UKB·FinnGen GO 환자/대조군 수, 약물 표적 후보 유전자 4,462개, 도구변수가 있는 유전자 2,544개(SNP–유전자 쌍 6,135개), BBJ에서 추정 가능한 유전자 2,234개, BBJ screen 통과 유전자 13개, screen 기준 *P* < 0.05/2,544 = 1.965×10<sup>−5</sup>.

## 5. 표와 그림

### 표

| 최종 | 표 | 변경 |
|---|---|---|
| 5 | 데이터셋 표 | 3행 → 노출 1행 + 질환 5행; endpoint 코드와 슬라이드 표기, ancestry, 역할; post hoc 표시; 셀 줄바꿈 정리 |
| 7 | colocalization 가설 표 | 덱 표 서식으로 다시 만들고 H0 문장 정리; 가리던 설명선 대신 옆에 설정 상자(Settings) |
| 10 | *TSHR* 표(신규) | 5개 자료 × OR (95% CI), *P*, PP.H3, PP.H4, PP.H4 (p12 = 10<sup>−6</sup>); 0.80 이상은 파란색, H3 우세는 주황색 |
| 11 | *IGF1R* 표 | 3행 → 5행; IVs와 PP.H2 열 추가; FinnGen GD 강조 |
| 12 | 13개 발견 유전자 표 | MHC 5개(BBJ만) → 13개 × BBJ와 FinnGen GD OR·*P*, focused support, FinnGen GD PP.H4 |
| 20 | 조직 확장 표(신규) | 유전자–조직 조합별 도구변수와 5개 자료 결과; 30개 추정 불가 비교는 음성 결과가 아님 |

### 그림

| 최종 | 그림 | 변경 |
|---|---|---|
| 4 | MR 공식 | 흐린 그림 → 슬라이드 텍스트(*OR* = exp(*β*<sub>disease</sub> / *β*<sub>expression</sub>)) |
| 8 | Figure 1B (Manhattan) | 원고 Figure 1의 B 패널만, 슬라이드 크기(6.0 × 4.3 in)로 다시 그림; 글자 10.5–12 pt; *TSHR*·*IGF1R*·*CTLA4* 굵게; 라벨 겹침 정리 |
| 9 | Figure 2 (MR + colocalization) | 5개 자료 × 3개 유전자; 숫자는 그림에서 빼고 표(10·11번)에 둠; H4 막대를 0부터 쌓아 점선이 곧 PP.H4 = 0.80 기준이 되게 함 |
| 18 | Figure 3A (히트맵, 신규) | 가로형(5개 자료 × 24개 유전자); PP.H4 0.001 미만은 빈칸, 0.80 이상은 테두리, p12 = 10<sup>−6</sup>에서도 0.80 이상이면 점 |
| 19 | Figure S2 (leave-one-out, 신규) | 두 패널 가로 배치; 행 이름과 *P* 10 pt, 축 11 pt, 패널 제목 12.5 pt; 굵은 *P* 값은 nominal *P* < 0.05 |
| 21 | Figure S1 (안와 RNA-seq) | 원고 그림 그대로 사용(개별 표본 자료가 있어야 다시 그릴 수 있음) |

다시 그린 그림은 모두 `Supplementary_Dataset.xlsx`에서 원고 그림 스크립트(`build_figures_20261014.py`)와 같은 데이터로 그렸고, 글꼴(Carlito, Calibri와 같은 폭), 크기, 배치를 슬라이드에 맞춤(Figure 1은 B 패널만, Figure 3은 A 패널을 가로로, Figure S2는 두 패널을 가로로). 슬라이드용 그림은 [`figures/`](figures/)에 있음.

## 6. 형식과 디자인

- 슬라이드 번호는 따로 적지 않으면 원본 기준
- 유전자 기호와 *P*를 이탤릭으로 바꿈(원본 슬라이드 텍스트에는 이탤릭이 없었음; 그림 이미지는 제외); 단백질 이름 IGF-1R은 관례대로 이탤릭 없이 둠
- 쪽번호: 손으로 입력한 번호(3, 3, 3, 6, 5 … 11, 16, 17, 12) → 자동 번호
- 내용을 가리던 도형 삭제: 5번 폭발 도형, 9번 Bonferroni 설명선과 클릭으로 나타나던 요약 상자, 11번 설명선; 14번은 두 줄로 넘친 소제목과 본문 글자 겹침을 해소
- 클릭 애니메이션을 모두 없앰(원본 9번 요약 상자, 11번 설명선, 12번 색 강조 상자); 노트의 "(클릭)" 표시도 지움
- 원본 4–7번 일부(특히 7번 데이터셋 표)와 13–15번에 섞여 있던 맑은 고딕 → Calibri로 통일; 소제목 위치·크기 통일(위에서 0.92 in, 16 pt)
- 최종 6번(원본 8번) 상단 두 상자를 본문 여백 안으로 정리; 2번 그림을 오른쪽 여백 안으로 축소
- 슬라이드 하단 참고문헌 표기 추가: 최종 2번(Chin, Tsui, Smith, Douglas 2020, Nelson, Minikel), 3번(Ji, Zhang, Ke), 7번(Giambartolomei, Wallace), 13번(Brand, Colobran, Marín-Sánchez), 14번(Douglas 2020·2024, Krieger); 6번 Finan에 연도 추가; 원본 13번 노트의 MHC sequencing consortium 1999 인용은 삭제
- 그림 대체 텍스트(alt text) 추가; 옛 저자·인과 문구가 남아 있던 도형 이름 정리
- 전공의 디자인(남색 상단 바, 회색 제목선, 파란 소제목, Calibri, 표 머리 색)은 그대로 유지

## 7. 발표자 노트

- 21장 모두 원고 기준으로 다시 씀(발표자 김민선, 말하는 순서대로)
- 원본 노트의 오류 수정: TSHR 억제제 효과, *CTLA4* 설명, 인과 표현, "단백질" → 혈액 mRNA, "태어날 때" → "수정(잉태) 시점", OR 설명은 "두 효과를 나눈 값(log OR)을 지수로 바꾸면 OR"로 풀어 씀
- 오타·끊긴 문장 수정: IGFR1 → IGF1R, FinGen → FinnGen, 잇습니다 → 있습니다, "맨델리안 무작위와" → "멘델리안 무작위화", 쓰다 만 "그 안에서 bon" 줄 삭제
- 내부 파일명 삭제(Stage2_14_genes_primary.csv 등)
- 질의응답 대비 내용 추가: 사전확률 민감도, *CTLA4*의 엇갈린 패턴(*P*가 작은 BBJ·FinnGen GD에서는 H3 우세, *P*가 nominal인 UKB와 두 GO 자료에서는 PP.H4 0.95 이상), MHC의 복잡한 LD, *HLA-G*를 넣은 이유, leave-one-out 해석
- "[발표 전 확인]" 메모 2곳: 2번 그림 출처, 감사 문구와 자료 사용 허락

## 8. 검증

- **1차 검토:** 원본 17장을 최종 원고·표·워크북·CLAUDE.md 고정 문구와 대조해 필수 수정 11개와 세부 항목을 정리함(review_checklist.md)
- **2차 독립 검증:** 대상은 중간 수정본 v3이며, 그 지적을 반영한 것이 최종본 v4임. 숫자·표현·누락 항목·그림·레이아웃(레이아웃만 2개) 분야로 나눈 검토 6개를 따로 진행했고, 각 검토의 지적은 반박 검증을 한 번 더 거침. 70건 확인(중복 포함), 5건 기각. 2건(문서 메타데이터, Figure S1 글꼴)은 9절 남은 확인 사항으로 넘기고 나머지를 모두 반영해 v4를 만듦
- 확인된 지적의 성격: 표현의 정확성("the two largest datasets" → "the two datasets with the most cases", 14.6% 검정력 문장, r² 0.46은 rs12101255에만 해당, F > 10은 선택 기준이 아닌 확인), 빠진 단서(GO 결과는 주로 GD 감수성을 반영할 수 있음, colocalization 강조는 post hoc, MHC의 복잡한 LD), 레이아웃(여백, 줄바꿈, 소제목 위치), 삭제된 도형을 가리키던 클릭 애니메이션. **핵심 숫자(OR, CI, *P*, PP.H4)의 오류는 없었음**
- **v4 형식 검사:** PowerPoint 구조 검사(validate.py) 통과; LibreOffice로 21장을 이미지로 렌더링해 확인
- **원고 패키지 무결성 점검(audit) 통과:** 패키지 파일 357개, 원고 고정 문구 61개 확인. 발표자료는 이 점검의 대상이 아님
- **이 보고서의 검증:** 원본 관련 서술, 최종본 관련 서술, 숫자·출처, 누락·일관성, 가독성 검토 5개와 반박 검증을 거쳐 확인된 지적을 반영함

검증 기록: [`review_checklist_20261025.md`](review_checklist_20261025.md), [`deck_review_findings_20261025.md`](deck_review_findings_20261025.md)

## 9. 남은 확인 사항

- **2번 기전 그림 출처:** 슬라이드의 "Figure: add source" 자리에 출처를 넣어야 함(이 보고서에서는 그림을 가림)
- **감사 문구와 안와 자료 사용 허락:** 저자 확인 후 확정(PNUH IRB 2104-018-102)
- **원고 상태:** 발표 시점에 투고했다면 표지의 "Manuscript in preparation"을 고침
- **발표자 소속:** 노트 첫 줄을 "가톨릭대학교 안과 김민선입니다"로 적은 것은 원고의 소속 2에 병원명 없이 "가톨릭대학교 의과대학 안과"만 있기 때문임. 실제 소속에 맞게 고치면 됨
- **문서 메타데이터:** 원본에 있던 병원 문서보안(Fasoo) 추적 정보가 최종본에도 그대로 남아 있음. 외부에 공유하기 전에 기관 정책 확인. 문서 속성의 슬라이드 수도 17장으로 남아 있음(PowerPoint에서 한 번 저장하면 갱신됨)
- **21번 Figure S1:** 원고 그림 그대로라 글꼴이 다른 슬라이드와 다름
- **발표 시간:** 본문 16장 기준. 시간이 정해지면 분량 조정

## 부록. 슬라이드 비교(원본 → 최종)

왼쪽은 원본, 오른쪽은 최종본. 2번의 외부 기전 그림은 출처가 확인되지 않아 이 보고서에서는 회색 상자로 가림(실제 슬라이드에는 그대로 있음).

| 원본 | 최종 |
|---|---|
| <img src="thumbs/orig-01.jpg" width="400" alt="orig slide 1"><br><sub>원본 1</sub> | <img src="thumbs/final-01.jpg" width="400" alt="final slide 1"><br><sub>최종 1 · 저자·소속·각주·원고 상태</sub> |
| <img src="thumbs/orig-02.jpg" width="400" alt="orig slide 2"><br><sub>원본 2</sub> | <img src="thumbs/final-02.jpg" width="400" alt="final slide 2"><br><sub>최종 2 · 근거 문장, 치료 서술, 출처</sub> |
| <img src="thumbs/orig-03.jpg" width="400" alt="orig slide 3"><br><sub>원본 3</sub> | <img src="thumbs/final-03.jpg" width="400" alt="final slide 3"><br><sub>최종 3 · 선행연구와 목적</sub> |
| <img src="thumbs/orig-04.jpg" width="400" alt="orig slide 4"><br><sub>원본 4</sub><br><img src="thumbs/orig-05.jpg" width="400" alt="orig slide 5"><br><sub>원본 5</sub> | <img src="thumbs/final-04.jpg" width="400" alt="final slide 4"><br><sub>최종 4 · MR 개념(5번 복제 슬라이드는 삭제)</sub> |
| <img src="thumbs/orig-06.jpg" width="400" alt="orig slide 6"><br><sub>원본 6</sub><br><img src="thumbs/orig-07.jpg" width="400" alt="orig slide 7"><br><sub>원본 7</sub> | <img src="thumbs/final-05.jpg" width="400" alt="final slide 5"><br><sub>최종 5 · 노출 자료와 질환 데이터셋 5개</sub> |
| <img src="thumbs/orig-08.jpg" width="400" alt="orig slide 8"><br><sub>원본 8</sub> | <img src="thumbs/final-06.jpg" width="400" alt="final slide 6"><br><sub>최종 6 · 도구변수와 기준값</sub> |
| <img src="thumbs/orig-11.jpg" width="400" alt="orig slide 11"><br><sub>원본 11</sub> | <img src="thumbs/final-07.jpg" width="400" alt="final slide 7"><br><sub>최종 7 · colocalization 방법</sub> |
| <img src="thumbs/orig-09.jpg" width="400" alt="orig slide 9"><br><sub>원본 9</sub><br><img src="thumbs/orig-10.jpg" width="400" alt="orig slide 10"><br><sub>원본 10</sub> | <img src="thumbs/final-08.jpg" width="400" alt="final slide 8"><br><sub>최종 8 · BBJ screen</sub> |
| <img src="thumbs/orig-12.jpg" width="400" alt="orig slide 12"><br><sub>원본 12</sub> | <img src="thumbs/final-09.jpg" width="400" alt="final slide 9"><br><sub>최종 9 · 3개 유전자 개요(새 Figure 2)</sub> |
| <img src="thumbs/orig-12.jpg" width="400" alt="orig slide 12"><br><sub>원본 12</sub><br><img src="thumbs/orig-14.jpg" width="400" alt="orig slide 14"><br><sub>원본 14</sub> | <img src="thumbs/final-10.jpg" width="400" alt="final slide 10"><br><sub>최종 10 · *TSHR* 표</sub> |
| <img src="thumbs/orig-15.jpg" width="400" alt="orig slide 15"><br><sub>원본 15</sub> | <img src="thumbs/final-11.jpg" width="400" alt="final slide 11"><br><sub>최종 11 · *IGF1R* 표</sub> |
| <img src="thumbs/orig-13.jpg" width="400" alt="orig slide 13"><br><sub>원본 13</sub><br><img src="thumbs/orig-10.jpg" width="400" alt="orig slide 10"><br><sub>원본 10</sub> | <img src="thumbs/final-12.jpg" width="400" alt="final slide 12"><br><sub>최종 12 · 13개 발견 유전자</sub> |
| <img src="thumbs/orig-14.jpg" width="400" alt="orig slide 14"><br><sub>원본 14</sub> | <img src="thumbs/final-13.jpg" width="400" alt="final slide 13"><br><sub>최종 13 · 고찰: *TSHR*</sub> |
| <img src="thumbs/orig-15.jpg" width="400" alt="orig slide 15"><br><sub>원본 15</sub> | <img src="thumbs/final-14.jpg" width="400" alt="final slide 14"><br><sub>최종 14 · 고찰: *IGF1R*와 치료</sub> |
| <img src="thumbs/orig-16.jpg" width="400" alt="orig slide 16"><br><sub>원본 16</sub> | <img src="thumbs/final-15.jpg" width="400" alt="final slide 15"><br><sub>최종 15 · 결론·제한점</sub> |
| <img src="thumbs/orig-17.jpg" width="400" alt="orig slide 17"><br><sub>원본 17</sub> | <img src="thumbs/final-16.jpg" width="400" alt="final slide 16"><br><sub>최종 16 · 감사</sub> |
| 원본 없음(신규) | <img src="thumbs/final-17.jpg" width="400" alt="final slide 17"><br><sub>최종 17 · 보충 슬라이드 구분</sub> |
| 원본 없음(신규) | <img src="thumbs/final-18.jpg" width="400" alt="final slide 18"><br><sub>최종 18 · 24개 유전자 colocalization</sub> |
| 원본 없음(신규) | <img src="thumbs/final-19.jpg" width="400" alt="final slide 19"><br><sub>최종 19 · leave-one-out</sub> |
| 원본 없음(신규) | <img src="thumbs/final-20.jpg" width="400" alt="final slide 20"><br><sub>최종 20 · 조직 확장</sub> |
| 원본 없음(신규) | <img src="thumbs/final-21.jpg" width="400" alt="final slide 21"><br><sub>최종 21 · 안와 RNA-seq</sub> |
