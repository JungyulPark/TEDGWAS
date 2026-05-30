# [Claude → Gemini] v5 Manuscript 업데이트 통지 (2026-05-27)

Claude verified canonical version이 기준입니다. Gemini 자율 작성본이 있으면 이 버전으로 reconcile하세요.
Canonical 파일: 07_manuscript/MANUSCRIPT_TED_TRAP_SUBMISSION.md

## 1. References — 27개로 확정 (전부 PubMed 검증, 등장순서 1~27)
- JEI Springer 스타일: 본문 [n], 목록 "저자 (연도) 제목. 잡지약어 권:페이지"
- 신규/확정 항목 중 주의:
  - #9 Ji 2025 (Sci Rep 15:37875, PMID 41162548)
  - #10 Li 2025 (Transl Vis Sci Technol 14:34, PMID 40576426) — IOVS-계열 안과 MR 선행
  - #26 GTEx Consortium 2020 (Science 369:1318-1330, PMID 32913098)
  - #27 Nelson 2015 (Nat Genet 47:856-860, PMID 26121088)
- 제외: Zhu/apeglm (본문 anchor 없음 → 인용 안 함)
- 데이터 인용: BBJ Graves + UKB hyperthyroid 둘 다 Sakaue 2021[#13], FinnGen=Kurki 2023[#14], eQTLGen=Vosa 2021[#12]
  → 공개 요약통계만 사용. UKB Application Number 불필요. "UK Biobank"를 title/abstract에 넣을 필요 없음.

## 2. Tables (main) — 3개로 확정
- Table 1: 데이터 소스 (변경 없음)
- Table 2: 강화됨 — backbone 3유전자 × 3 outcome, 통계 전체
  (N IV, OR 95%CI, IVW/Wald P, weighted median P, MR-Egger intercept P, Cochran Q P, coloc PP.H4, tissue)
  single-IV는 "NA (single IV)"로 정직 표기. CTLA4 UKB/FinnGen Cochran Q P 유의(0.001/0.049) 명시.
- Table 3: 신규 추가 (main 승격) — 13 Bonferroni hits 전체
  (N IV, OR 95%CI, BBJ P, coloc PP.H4 BBJ/FinnGen, 분류)
- 모든 값 출처: TaskD_03d_MR_all_outcomes_merged_v1, TaskE_01_coloc_results_v1, TaskE_02_candidate_coloc_v1 (Claude 검증)

## 3. Supplementary — 번호 유지 (S1~S7), 재번호 없음 ★중요★
- Table S1: **확장판으로 재정의** (중복 해소)
  Table 3 = BBJ discovery 요약(OR+coloc); Table S1 = cross-outcome 확장
  (각 hit의 BBJ/UKB/FinnGen β·P 3개 + coloc lead SNP + 분류)
- S2~S7: 변경 없음 (S2 coloc, S3 fine-mapping, S4 instruments, S5 chr16, S6 sensitivity, S7 STROBE-MR)
- ★ S 번호 재정렬 절대 하지 말 것. S1 삭제·재번호 아님. S1을 확장 표로 교체만.

## 4. Figures — 변경 없음
- Figure 1 = B&W 버전 사용 (Gemini R-version 폐기 유지)
- Figure 2 (volcano), Figure 3 (composite) — 기존 locked

## 5. Discussion — 인용 보강 완료 (14개 고유 인용, D2/D3/D4/D7/D8 각 2~3개)
- 내용/논리는 동일, 문헌 anchor만 추가. 본문 문장 자체는 friend-reviewed 버전 유지.

## 6. 검증 통과 (Claude 자동)
- 본문 인용 1~27 완전 오름차순, 누락·중복 0, 목록 27개 전부 anchor 있음
- Table/Figure/Supp cross-reference 전부 일치
- 본문 단어수 ≈ 4,452 (JEI 4,500 한도 — 여유 적음, docx 시 재확인)
