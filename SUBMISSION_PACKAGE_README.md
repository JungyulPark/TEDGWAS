# TED-TRAP v5 — JEI 제출 패키지 (최종)

**Target journal:** Journal of Endocrinological Investigation (JEI)
**Article type:** Original Article
**Master MD5:** `340a5a694a58916c0adb33ec1a173492` (Discussion compressed)
**Date assembled:** 2026-05-28

---

## 제출 파일 목록

### 1. 본문 원고
- **`MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx`** — 제출용 Word (pandoc 변환, 표 11개 보존)
- `MANUSCRIPT_TED_TRAP_v5_MASTER.md` — 마크다운 정본 (편집/버전관리용)

원고 구성:
- Title, Running head, Authors/affiliations, Abstract (262 words), Keywords
- Introduction → Methods → Results → Discussion → Declarations → References (27)
- Figure Legends (Fig 1–3)
- Tables 1–3 (본문 내 inline)
- Supplementary: Tables S1–S7 (실데이터 포함), Figure S1–S2 legends

### 2. Cover letter
- **`COVER_LETTER_JEI.docx`** — 제출용 (IGF1R framing이 Abstract와 정렬됨)
- `COVER_LETTER_JEI_aligned.md` — 마크다운 원본

### 3. 그림 파일 (Gemini 로컬에서 별도 첨부 필요 — 이 패키지에 미포함)
JEI 제출 시 아래 그림을 고해상도(≥300 DPI)로 별도 업로드:
- Figure 1 — Study design schematic
- Figure 2 — Druggable-gene-wide BBJ discovery volcano (2,235 genes)
- Figure 3 — Multi-layer evidence integration (A coloc / B MR forest / C tissue)
- Figure S1 — TSHR fine-mapping & LD
- Figure S2 — Candidate colocalization
> ⚠️ Figure 3B β축 = Table 2 값(−2.10/−2.44/−2.33)과 scale 일치 확인 완료 (Gemini, TaskD_03d 원본 대조).

### 4. 내부 문서 (제출하지 않음)
- `INTERNAL_external_GEO_research_note.md` — 외부 GEO 검토 기록 (옵션 A 근거, rebuttal 대비)

---

## 제출 전 최종 점검표

| 항목 | 상태 |
|---|---|
| Placeholder (TODO/빈칸) | ✅ 0개 (S3 PIP, S5 chr16 좌표/r² 채움 완료) |
| 수치 정합성 12개 감사 | ✅ 통과 (OR=exp(β), PP.H4 cross-table, 13 hits Bonferroni, 표본수) |
| #1/#3/#4 reviewer-proofing | ✅ 적용 (IGF1R effector는 해석으로, absence-of-coloc 방어문장) |
| Cover letter ↔ Abstract 정렬 | ✅ "most consistent with pharmacologic effector biology" |
| 외부 GEO 결정 | ✅ 옵션 A (3개 코호트 비재현 → 미포함, limitation 정직 유지) |
| Figure 3B β scale | ✅ Table 2와 일치 (Gemini 원본 대조) |
| 표 개수 (본문 3 + S1~S7 + S3 PIP sub) | ✅ docx에 11 table 보존 |
| 참고문헌 / 본문 인용 | ✅ 27개, [1]→[27] 오름차순 |

## ⚠️ 선생님 직접 결정 사항 (제출 전)

1. **단어 수** — 본문 실측 **4,637 words** (Intro 744 + Methods 964 + Results 1,278 + Discussion 1,651; Discussion compressed from 1,769).
   - JEI 권장 "preferably 4,500" (Methods 포함 여부는 저널이 명시 안 함).
   - Methods 포함 기준이면 ~137 초과(권장선 근접); Methods 제외 기준이면 3,673로 여유 통과. Discussion 문단 5·7 압축 + conclusion tissue 완화 완료.
   - 압축 후보: Discussion이 1,770으로 전체의 37%. "broader principle" 문단(~150) + "methodological features" 문단(~100)이 부연/중복이라 핵심 손실 없이 압축 가능.
2. **Abstract** — 현재 **262 words**. JEI abstract 한도(보통 250) 확인 필요. 12단어 내외 초과 가능성.

---

## 빌드 정보
- 변환: pandoc 3.1.3, markdown+pipe_tables → docx
- 검증: docx 표 11개 보존, 수치·편집 전부 반영 확인
