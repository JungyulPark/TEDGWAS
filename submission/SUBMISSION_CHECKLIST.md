> **Superseded by the 2026-09-05 candidate.** Use [candidate_20260905/README.md](candidate_20260905/README.md) and its revision report. The older files and build instructions below are historical. The older Word files, cover letter and figure copies were moved to `archives/submission_pre20260905/` on 2026-09-05 so that only one submission package remains under `submission/`; `scripts/29_build_submission_package.py` was retired with them and the builder is now `candidate_20260905/reproducibility/build_documents.py`.

# 제출 패키지 — TED-TRAP Paper 1

**Manuscript:** Genetic susceptibility and therapeutic target biology diverge at *TSHR* and
*IGF1R* in Graves disease and thyroid eye disease

**Authors:** Jungyul Park¹; Min-Seon Kim²; Kyung-Hwa Shin³\*; Suk-Woo Yang¹\* (\*corresponding)

**Target venue:** Endocrine Connections

Built from `MANUSCRIPT_TED_TRAP_v5_MASTER.md` (md5 `8ead93d32c35404a5bd30974e192b73c`) by
`python3 scripts/29_build_submission_package.py`. The master is the only source of truth;
every file below is derived from it and must be rebuilt, never hand-edited.

## 업로드할 파일

| # | 파일 | 내용 |
|---|---|---|
| 1 | `MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx` | 표제지 · 초록 · 본문 · 선언 · 참고문헌 26 · 그림 범례 · Tables 1–3 |
| 2 | `SUPPLEMENTARY_MATERIAL.docx` | Tables S1–S9 + Supplementary Figure 범례 |
| 3 | `COVER_LETTER_EndocrineConnections.docx` | 커버레터 |
| 4 | `figures/Figure1.png` | 연구 흐름도 + TSHR/IGF1R 비교 (흑백, CONSORT 방식) |
| 5 | `figures/Figure2.png` | BBJ discovery MR volcano |
| 6 | `figures/Figure3.png` | 3 outcome colocalization + MR forest (P값 표기) |
| 7 | `figures/FigureS1.png` | 안와조직 TPM (서술적) |
| 8 | `figures/FigureS2.png` | 후보 유전자 colocalization, 3 outcome |

`KSO_oral_abstract_KR_EN.docx`는 학회 초록으로, **저널에 업로드하지 않습니다.**
JEI용 커버레터는 삭제했습니다 — 구 저자명과 4층 근거 기준이라 재사용 불가. 필요하면
EC 버전에서 다시 만드십시오.

## 내용 확인 (자동 검사 통과)

- 본문 **4,959 단어** (한도 5,000) · Abstract **237 단어, 한 문단**
- 참고문헌 **26개**, 전부 인용, 범위 밖 없음
- 본문 표 3개 / 부록 표 9개 (S1–S9) / 그림 5개
- 근거층 **2개** — MR과 colocalization. Fine-mapping은 철회, 조직은 서술적 관찰
- `scripts/audit_paper1_integrity.py` 42개 검사 ALL PASS
- `scripts/29_build_submission_package.py` 패키지 검사 ALL PASS

## 제출 전 저자가 채워야 할 것

1. **Kyung-Hwa Shin (공동교신)** — 이메일, 전화번호, 학위(`MD` vs `MD, PhD`)
2. **Min-Seon Kim** — 소속에 병원명 (현재 "Department of Ophthalmology, College of
   Medicine, The Catholic University of Korea"로 병원 누락)
3. 저자 4인 **ORCID**, 공동교신저자 표기가 포털에서 허용되는지 확인
4. 포털이 생성한 PDF에서 표·그림·수식·줄번호 육안 확인
5. **GitHub 옛 커밋 `88eabcb`** — IRB 데이터가 서버에 잔존. 저장소 공유 전 처리
   (`../internal/GITHUB_SECURITY_CLEANUP.md`)

---

# Submission package — TED-TRAP Paper 1 (v5, final)

**Manuscript:** Genetic susceptibility and therapeutic target biology diverge at TSHR and
IGF1R in Graves disease and thyroid eye disease.

**Authors:** Jungyul Park¹; Min-Seon Kim²; Kyung-Hwa Shin³\*; Suk-Woo Yang¹\* (\*corresponding)

Master source of truth: `MANUSCRIPT_TED_TRAP_v5_MASTER.md`
(md5 `8ead93d32c35404a5bd30974e192b73c`, placeholders = 0, integrity audit ALL PASS).

**Target venue:** Endocrine Connections (1st choice) → Endocrine → JEI → Frontiers-Endo/BMC.
See `../internal/SUBMISSION_VENUE_CASCADE.md`.

## Files to upload

| Item | File | Status |
|---|---|---|
| Main manuscript (Word) | `MANUSCRIPT_TED_TRAP_v5_SUBMISSION.docx` | ✅ final — 13 tables embedded (Tables 1–3 + S1–S8), LibreOffice-validated |
| Cover letter (1st choice) | `COVER_LETTER_EndocrineConnections.docx` | ✅ **Endocrine Connections version** — sound-science framing, power/limitations, dual-use disclosure |
| Cover letter (fallback) | `COVER_LETTER_JEI.docx` | ⚠️ JEI version — **author list still says Yae-Eun Kang; fix before use** |
| Figure 1 | `figures/Figure1.png` | ✅ **claim figure** — screening funnel with reasons (A) + TSHR/IGF1R divergence (B); rebuilt 2026-08 |
| Figure 2 | `figures/Figure2.png` | ✅ BBJ MR volcano (2,235 genes, 13 hits, P = 1.965×10⁻⁵) |
| Figure 3 | `figures/Figure3.png` | ✅ composite A/B/C — **regenerated 2026-08** with corrected coloc legend |
| Figure S1 | `figures/FigureS1.png` | ✅ TSHR fine-mapping (rs179252 anchor) |
| Figure S2 | `figures/FigureS2.png` | ✅ candidate coloc — **regenerated 2026-08** with corrected legend |
| Conference abstract (not for journal) | `KSO_oral_abstract_KR_EN.docx` | ✅ KSO oral, KR + EN |

Tables 1–3 and Supplementary S1–S8 are embedded in the main `.docx` (no separate table files).

## Content inventory (verified)
- References 27; in-text citations [1]→[27], all cited.
- Main tables 3; Supplementary tables **S1–S8**; figure legends 5.
- OR = exp(β) verified for all 9 backbone estimates.
- Abstract **235 words, single paragraph** with inline structured headings (Endocrine Connections format).
- Main text **4,774 words** (Intro 655 / Methods 1,077 / Results 1,506 / Discussion 1,536) — under the 5,000-word
  Endocrine Connections limit. Recount with `python3 ../scripts/26_wordcount_main_text.py`.

## Substantive upgrades applied in the final round (2026-08)
1. **Table S8 — detectable effect range.** Design-based minimum detectable effect per gene/outcome
   (|β|min = (z₁₋α/₂ + z₀.₈₀) × SE) from the locked standard errors. Makes `robust_novel = 0`
   interpretable (median detectable OR 2.55 at the discovery threshold) and shows IGF1R's effect fell
   **below** its detection threshold in every outcome, so its sub-threshold status is expected rather
   than evidence against a genetic contribution. Script: `../scripts/22_power_analysis_null_screen.py`.
2. **Colocalization H1/H2 labels corrected (text + figures).** The pipeline passes dataset1 = outcome
   GWAS and dataset2 = *cis*-eQTL, so PP.H2 = "eQTL association only"; the manuscript and both figure
   builders had it inverted. 10 text occurrences fixed; Figures 3 and S2 rebuilt
   (`../scripts/23_rebuild_figures_coloc_label_fix.py`). No plotted value changed. See
   `../FIGURE_VERIFICATION.md`.
3. **FinnGen GO ascertainment overclaim removed.** GO cases are ascertained among Graves disease
   patients and compared with population controls, so that outcome re-measures GD susceptibility.
   TSHR is now framed as an anchor **for GD** whose signal is recovered in a TED-enriched Graves ophthalmopathy outcome,
   with a TED-specific effect explicitly not established (new Limitation "Sixth").
4. **New Limitations "Seventh" and "Eighth":** colocalization method limits (single-causal-variant
   assumption; one European allele-frequency vector used for both datasets including the East Asian
   discovery outcome; no colocalization power analysis) and the positive-control operating
   characteristic (CTLA4 colocalized only in the TED outcome, PP.H4 = 0.978 / rs1863800, and not in
   discovery, PP.H3 = 0.794 / rs231811 — one control does not establish filter sensitivity).
5. **Conceptual positioning made explicit:** the general "drug target ≠ susceptibility locus"
   principle is acknowledged as already established; the contribution is the worked
   within-disease demonstration.

## External review round 2 — applied 2026-08

| # | Item | Applied |
|---|---|---|
| 1 | Main text ≤ 5,000 words | ✅ 5,155 → **4,772** (Intro gaps paragraph, Results Methods-restatement, Discussion general-principle and TSHR paragraphs, Limitations duplication) |
| 2 | Abstract as a single paragraph | ✅ one paragraph, inline `Purpose:/Methods:/Results:/Conclusion:` headings, 235 words |
| 3 | IGF1R "nominal" precision | ✅ Abstract, Figure 1 and Figure 3 legends now say nominal in **BBJ and UKB but not FinnGen** (P = 0.18) |
| 4 | Soften the final Conclusion | ✅ Conclusion and the IGF1R Discussion paragraph now read "more compatible with … than with …", explicitly not establishing effector status nor excluding a modest genetic contribution |
| 5 | Figure 1 "identical instruments" error | ✅ **factual fix** — TSHR and IGF1R do not share instruments. Legend and the rendered Panel B subtitle now read "the same outcome hierarchy and analytic framework, with gene-specific *cis*-eQTL instruments"; Figure 1 regenerated |
| 6 | Complete author / corresponding-author info | ⛔ **not applied — requires the author's own data** (see to-do 1 above) |
| 7 | Justify "prespecified" | ✅ retained only for the outcome hierarchy (documented in `TrackA_MR/v5_upgrade/03_decision/TaskC_pre_analysis_plan_v1.md`); backbone genes now described as designated *a priori on biological and therapeutic grounds*, with a new Methods sentence stating that basis |
| 8 | Soften single-control tissue "significance" | ✅ Results and Figure 3 legend now report an *exploratory differential-expression signal* that "is not confirmatory because the analysis included only one control sample" |
| 9 | Stop calling FinnGen a "case series" | ✅ all 3 occurrences → "TED-enriched Graves ophthalmopathy outcome" |
| 10 | Update STROBE-MR Table S7 limitation row | ✅ item 18 now lists all ten limitations carried in the Discussion |

## External review round 3 — applied 2026-08

| # | Item | Applied |
|---|---|---|
| 1 | IGF1R summary overstated | ✅ **PP.H2 was being over-read as "no disease association" in 6 places.** IGF1R *does* have a nominal MR association (BBJ P = 0.021, UKB P = 0.012); PP.H2 concerns only the shared-variant hypothesis. Fixed in the Results paragraph, the Figure 3 legend, the Table 2 footnote (which sat directly under a row printing P = 0.021), the Figure 1 Panel B cell (rendered one row below "β = +0.45, P = 0.021"), the web review card, and `FIGURE_VERIFICATION.md` — the record that had seeded the wording. The web card now reads "IGF1R — effector-compatible profile / coloc PP.H2 — no shared-variant expression colocalization" |
| 2 | "every evidence layer" is not true | ✅ fine-mapping (SuSiE) was run only at *TSHR*. Methods now: "evaluated across the MR, colocalization, and orbital tissue layers … with additional fine-mapping performed for the single-instrument *TSHR* locus" |
| 3 | Power "excludes" too strong | ✅ only 35.6% of genes were powered for OR ≥ 2.0. Results and the Figure 1 legend and rendered verdict strip now say the null **constrains** additional large effects, particularly among well-powered genes, and **does not exclude** moderate ones — matching Table S8's own footnote |
| 4 | Corresponding-author details | ⛔ still outstanding (to-do 1) |
| 5 | Main text further trimmed | ✅ 4,886 → **4,772** words; Abstract 241 → **235** |
| 6 | Review page vs submission file | ✅ verified: the submission `.docx` opens on the running head → title → authors → Abstract, with **no** web chrome (checked all 16 zip parts: 0 hits for "Under review", "claudeusercontent", "Jump to", "at a glance"; no external hyperlink targets). The web page's masthead now states "Review copy — the journal file is the .docx, not this page" |
| 7 | **Stale duplicate master removed** | ✅ `submission/MANUSCRIPT_TED_TRAP_v5_MASTER.md` was a tracked mirror frozen at commit `35c8c22` (md5 `ade74638…`) — it still carried "identical instruments" and "case series". Deleted; the markdown source of truth is `../MANUSCRIPT_TED_TRAP_v5_MASTER.md`. The integrity audit now **fails** if a second master copy reappears |
| 8 | **Table S6 footnote credited *CTLA4* with support it does not have** | ✅ found by adversarial verification, not by the reviewer. The footnote read "Single-instrument loci (TSHR all outcomes; CTLA4 in BBJ) are supported by colocalization and fine-mapping" — but *CTLA4* was never fine-mapped (SuSiE was run only at *TSHR*) and its BBJ posterior favours **distinct** causal variants (PP.H3 = 0.794, PP.H4 = 0.206). Rewritten to state exactly what supports what |
| 9 | Table 2 footnote defined PP.H2 incorrectly | ✅ a sentence added earlier in this same round said "PP.H2 concerns the shared-variant hypothesis". PP.H2 is the *cis*-eQTL-only hypothesis (PP.H4 is the shared-variant one), which contradicted Methods, the Figure 3 and S2 legends and the Table S2 header. Now: "Colocalization tests only the shared-causal-variant hypothesis, so a PP.H2-dominant posterior does not contradict the nominal MR association reported above" |

## External review round 4 — applied 2026-08

| # | Item | Applied |
|---|---|---|
| 1 | **Methods outcome-hierarchy sentence was factually wrong** | ✅ it read "a prespecified hierarchy of decreasing sample size but increasing TED specificity", but the case counts run 2,809 (BBJ) → **3,731** (UKB) → 858 (FinnGen) — not monotonically decreasing — and TED specificity *dips* at the UKB hyperthyroidism step. Replaced with "a prespecified **functional** hierarchy", and the follow-on sentence now reads "Because adequately powered TED-specific GWAS are unavailable, this design prioritized discovery in Graves disease, assessed generalizability in a broader hyperthyroidism phenotype, and finally examined signal recovery in a smaller TED-enriched outcome" |
| 2 | Same claim in superseded cover-letter drafts | ✅ `COVER_LETTER.md` and `TrackA_MR/v5_upgrade/07_manuscript/COVER_LETTER.md` carried the identical sentence. Both corrected. The letter actually being submitted (`COVER_LETTER_EndocrineConnections.md`) never contained it |
| 3 | Corresponding-author details | ⛔ still outstanding (to-do 1) |

The integrity audit now asserts the three case counts are present and **fails** if the hierarchy is
described as monotone in sample size or TED specificity.

## Author to-do before clicking submit
1. **Author and corresponding-author details — BLOCKING, needs the author's own information.**
   The title page is still incomplete:
   - Kyung-Hwa Shin (co-corresponding): **no e-mail, no telephone**, and the degree is written as
     "MD" — confirm whether it should be "MD, PhD". Target format, matching the first
     corresponding author: `Kyung-Hwa Shin, MD[, PhD], Department of Laboratory Medicine, Pusan
     National University Hospital, Busan, Republic of Korea. E-mail: …; Tel: …`
   - Min-Seon Kim's affiliation ² reads "Department of Ophthalmology, College of Medicine, The
     Catholic University of Korea" with **no hospital named** — confirm the actual hospital.
   - ORCIDs for all four authors; confirm the portal accepts two corresponding authors; confirm all
     co-authors approve the author order.
   - Confirm IRB 2104-018-102 covers the orbital RNA-seq component and this investigator team.
2. **Cover letter** — use `COVER_LETTER_EndocrineConnections.docx` for Endocrine Connections. The JEI fallback letter still carries the superseded author name and must be corrected before reuse.
3. **Figures** — inspect each PNG at full size and confirm the journal's DPI/format requirement.
4. **Portal fields** — ORCID, author affiliations, ethics + informed-consent + data-availability
   statements (all already written in Declarations).
5. **Proof check** — after the portal builds the PDF, eyeball tables, equations, figures, and line
   numbers.
6. **GitHub security** — old commit `88eabcb` still resolvable server-side; recreate the repo or
   contact Support before the repo is shared. See `../internal/GITHUB_SECURITY_CLEANUP.md`.

## Not submitted (internal only)
`../internal/` — external-GEO decision, MVMR forensic verdict, literature landscape, venue cascade,
loop-engineering plan, security cleanup. Do not upload.
