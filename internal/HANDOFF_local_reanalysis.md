# 로컬 Claude Code 인계서 — TED-TRAP Paper 1 재분석 (Step 2–4)

_이 문서는 원본 데이터를 가진 로컬 머신의 Claude Code 세션이 그대로 읽고 실행할 수
있도록 자기완결적으로 작성했습니다. 이 대화의 맥락을 몰라도 됩니다._

---

## 0. 상황 요약 — 왜 재분석인가

TED-TRAP Paper 1(Endocrine Connections 제출 예정)이 외부 리뷰에서 **provenance 결함**
지적을 받았고, 저장소 원본 데이터로 검증한 결과 **5건이 사실로 확인**되어 **제출을 보류**했습니다.
전체 감사 기록: `internal/SUBMISSION_HOLD_provenance_audit.md`

| # | 결함 | 상태 |
|---|---|---|
| 1 | 도구변수 1개가 `P < 5e-8` 기준 위반 (`GAN`/`rs310019`) → 보고된 min F 14.2의 정체 | **원격에서 해결 완료** |
| 2 | FinnGen endpoint·증례수 미확정 (`s`가 placeholder) | **← Task 2, 로컬 필요** |
| 3 | Table S3 fine-mapping 수치가 원본과 불일치 | **← Task 3, 로컬 필요** |
| 4 | coloc이 p-value 기반 + 단일 EUR MAF, UKB coloc 부재 | **← Task 4, 로컬 필요** |

**절대 원칙:** 문장·표·그림을 먼저 고치지 마십시오. Task 2–4가 통과한 뒤에만
원고를 재생성합니다. 이 프로젝트의 지침은 `CLAUDE.md`에 있고, 핵심은
**정확하고 진실한게 생명 — 숫자를 기억으로 채우지 말 것**입니다.

---

## 1. 이미 끝난 것 (Step 1)

원격 세션에서 완료했습니다. `git pull` 하시면 받으실 수 있습니다.

**스크립트:** `scripts/27_verify_instrument_manifest.py`

```
python3 scripts/27_verify_instrument_manifest.py --write
```

**확인된 것.** 6,136개 instrument 중 **정확히 하나**가 선정 기준을 위반합니다.

```
GAN / rs310019   Z = -3.7684   P = 1.643e-04   F = Z^2 = 14.2008   N = 22,151
다음 최솟값       PPIB / rs895885            P = 4.986e-08   F = 29.724
```

`GAN`의 instrument는 그 SNP 하나뿐이라 애초에 2,545에 들어오면 안 됐습니다.
`TaskD_02_denominator_note.txt`가 "Genes with ≥1 valid cis-eQTL instrument
(**P<5e-8**, EUR-clumped): 2545"라고 적어 놓고 실제로는 위반 행을 포함하고 있으므로,
**D.2 단계의 필터에 버그가 있습니다.**

**정정된 수치.**

| 보고값 | 정정값 |
|---|---|
| 6,136 instruments | **6,135** |
| 2,545 MR-testable genes | **2,544** |
| 2,235 BBJ-estimable genes | **2,234** |
| min F = 14.2 | **29.72** |
| Bonferroni 0.05/2,545 = 1.9646e-05 | 0.05/2,544 = 1.9654e-05 |

**결론 영향 없음:** 두 threshold 모두 유효숫자 4자리에서 `1.965e-05`로 인쇄되고,
13개 hit 목록이 동일하며(C4A, CTLA4, HLA-A, HLA-DQA2, HSD3B7, IFNGR1, MAPKAPK5,
PRSS36, PSMB8, TNFSF14, TSHR, TUBB, VKORC1), backbone 3개 모두 영향 없습니다.
`GAN`의 BBJ MR은 β = 1.782, SE = 0.929, **P = 0.055**로 유의 근처도 아닙니다.

**산출물:**
- `TrackA_MR/v5_upgrade/04_druggable_mr/results/TaskD_02a_eqtl_instruments_snp_level_v2_verified.csv` (6,135행)
- `TrackA_MR/v5_upgrade/07_manuscript/supp_tables/TableS4_instruments_v2_verified.csv` (6,135행)

**로컬에서 추가로 해 주실 것 (Task 1b, 선택):**
D.2 필터가 왜 `GAN`을 통과시켰는지 원인을 찾아 주십시오. 같은 버그가 `CLUMP_FAILED`로
표시된 2개 유전자에도 영향을 줬을 수 있습니다(`TaskD_02_denominator_note.txt` 참조).
원인을 못 찾더라도 Task 2–4를 막지는 않습니다.

---

## 2. Task 2 — FinnGen 파일·endpoint·증례수 확정 (최우선)

### 문제

저장소 스크립트들이 **서로 다른 endpoint**를 참조하고, `858`이 **두 개의 다른 분모**에
붙어 있습니다.

| 파일 | endpoint | 분모 |
|---|---|---|
| `v5_upgrade/scripts/taskD_03_mr_3outcomes.R` | `finngen_R12_GRAVES_OPHT.gz` (broader) | `n_case = 858` |
| `v5_upgrade/scripts/taskE_01_coloc.R` | `finngen_R12_GRAVES_OPHT.gz` (broader) | `858 / 520387` |
| `v5_upgrade/scripts/taskE_02_coloc_candidates.R` | `finngen_R12_GRAVES_OPHT.gz` (broader) | `858 / 520387` |
| `scripts/03_mr_all_genes_full_sensitivity.R` | `finngen_R12_E4_GRAVES_OPHT` (**strict**) | `858 / 500348` |
| `scripts/06_coloc_all_loci.R` | `E4_GRAVES_OPHT` (**strict**) | `858 / 500348` |

그리고 PP.H4 = 0.986을 만든 스크립트에 이 주석이 **그대로 남아 있습니다**:

```r
# FinnGen s/N: CONFIRM control count. R12 GO ~858 cases. Total from file rows or meta.
# Placeholder — Gemini sets these from FinnGen R12 GRAVES_OPHT documentation:
FINNGEN_S <- 858 / 520387
FINNGEN_N <- 520387
```

`coloc.abf`에서 `type="cc"`일 때 `s`는 사후확률에 **직접** 들어갑니다. 미확인 값입니다.

외부 리뷰어는 FinnGen R12 공식 수치가 `E4_GRAVES_OPHT_STRICT` = **786 cases**,
`GRAVES_OPHT` = **894 cases**이며 858은 어느 쪽과도 맞지 않는다고 지적했습니다.
**이 주장을 검증 없이 받아들이지 마십시오 — 실제 파일과 Risteys로 직접 확인하십시오.**

### 해야 할 일

1. `c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz`가 실제로 존재하는지 확인하고,
   **헤더와 첫 몇 줄, 파일 크기, 다운로드 날짜**를 기록하십시오.
2. 같은 위치에 `E4_GRAVES_OPHT` 파일도 있는지 확인하십시오. 있다면 **어느 것이
   v5 결과를 만들었는지** 파일 mtime과 스크립트 실행 로그(`00_meta/analytical_log.md`)로
   대조하십시오.
3. FinnGen R12의 **공식 endpoint 문서**에서 해당 endpoint의 cases/controls를 확인하십시오
   (Risteys R12). 파일에 동봉된 manifest나 README가 있으면 그것을 우선하십시오.
4. `858`이 어디서 왔는지 추적하십시오. 후보: 다른 release(R11 등), phenotype 발생 건수와
   genotype-QC 후 GWAS case 수의 혼동, 또는 다른 논문에서 인용.

### 산출물

`internal/FINNGEN_endpoint_provenance.md`에 다음 표를 채워 주십시오. **추정 금지 —
확인된 것만 적고, 확인 못 한 항목은 "미확인"으로 두십시오.**

```
| 항목 | 값 | 출처 |
|---|---|---|
| 분석에 사용한 파일명 (전체 경로) |  | ls -la 출력 |
| 파일 크기 / mtime |  |  |
| endpoint code |  | 파일 헤더 또는 동봉 manifest |
| FinnGen release |  |  |
| cases |  | 공식 문서 |
| controls |  | 공식 문서 |
| total N |  | cases + controls |
| coloc 에 넣은 s |  | 재계산값 |
| 원고의 858 이 이 중 무엇과 일치하는가 |  |  |
```

### 통과 기준

- 파일 경로·endpoint·release·cases·controls가 **하나의 출처로** 확정될 것
- `s = cases / (cases + controls)`가 위 값에서 계산될 것
- 원고의 `858`과 일치하지 않으면 **일치하지 않는다고 명시**할 것

---

## 3. Task 3 — TSHR SuSiE 재실행

### 문제

원고의 Table S3 fine-mapping 서술이 저장소 원본과 **일치하지 않습니다.**

**원본** (`TrackA_MR/v5_upgrade/07_manuscript/supp_tables/TableS3_finemapping.csv`,
`TaskB_05_tshr_finemap_summary_v1.csv`, `WEEK1_MANIFEST.md`, `TaskC_decision_table_v1.csv`
— 넷 다 일치):

```
European (1000G), CS_1, SNPs=1, purity=1, rs179252, PIP=1, r2_EUR=1
European (1000G), CS_2, SNPs=1, purity=1, rs179248, PIP=1, r2_EUR=0.965
East Asian(1000G), CS_1, SNPs=1, purity=1, rs3783950, PIP=1
East Asian(1000G), CS_2, SNPs=1, purity=1, rs3783949, PIP=1
```

→ **EUR credible set 2개, 각각 SNP 1개, purity 1**

**원고 주장:**

| 원고 | 원본 |
|---|---|
| "a single 95% credible set" | 2개 |
| "purity = 0.993" | 1 |
| "log10 Bayes factor = 88.2" | 저장소에 없음 |
| 5개 변이 PIP 0.458 / 0.436 / 0.052 / 0.021 / 0.020 | 각 CS에 SNP 1개, PIP 1.0 |
| lead `rs11603529` | **저장소 어느 데이터 파일에도 없음** |
| `rs10137255` | **저장소 어느 데이터 파일에도 없음** |
| "~116-kb window" | 위 미확인 목록에서 유도됨 |

Table S3의 **clumping 표는 정확**합니다 (`TaskB_02_ld_clumping_sensitivity_v1.csv`와
행 단위로 일치). 문제는 PIP 블록입니다.

### 해야 할 일

`TrackA_MR/v5_upgrade/scripts/taskB_04_susie_finemap.R`를 **그대로 재실행**하고,
출력을 원본 그대로 저장하십시오. 스크립트를 고치지 말고, 먼저 무엇이 나오는지 보십시오.

재실행 시 반드시 함께 기록할 것:

- SuSiE run 날짜, `susieR` 버전, R 버전, random seed
- 입력 SNP 개수와 locus 좌표 (chr14 window, hg19)
- LD reference: 1000G Phase 3 어느 패널(EUR n=503 / EAS n=504), genome build
- **credible set 개수**와 각 CS의 membership
- **모든 SNP의 PIP** (상위 20개 이상)
- 각 CS의 **coverage**와 **purity**
- **rs179252의 실제 PIP**
- `susie_rss`의 **LD-vs-summary-statistics consistency 진단**
  (`susieR`의 `estimate_s_rss()` / `kriging_rss()` — reference LD와 summary stat이
  어긋나면 결과가 불안정합니다)

### 중요 — 범위 제한

eQTL만 fine-mapping한 것으로는 **질환과의 shared causal variant가 rs179252로 확정되지
않습니다.** 이를 주장하려면 최소한 coloc의 `SNP.PP.H4`와 H4 조건부 credible set이
필요하고(→ Task 4), 이상적으로는 outcome GWAS도 별도로 fine-mapping해야 합니다.

따라서 **현재 원고의 "fine-mapping resolved the shared causal variant to rs179252"는
철회 상태로 두십시오.**

### 산출물

- `TrackA_MR/v5_upgrade/02_tshr_finemap/results/TaskB_04_susie_rerun_<date>.csv` (전체 PIP)
- `internal/SUSIE_rerun_report.md` (위 항목 전부 + 원고 Table S3과의 대조표)

### 통과 기준

- credible set 개수, purity, rs179252 PIP가 **재실행 출력에서 직접** 나올 것
- `rs11603529` / `rs10137255` / `0.993` / `88.2`가 재실행 출력에 나타나는지 여부를 명시할 것
  (나타나지 않으면 원고 Table S3의 해당 블록은 **삭제 후 재작성** 대상)
- 살아남는 결론이 무엇인지 한 문장으로 적을 것
  (예상: "EUR CS 2개가 모두 단일 SNP이고 r²=0.965로 collinear → 독립 instrument 없음
  → 단일 IV locus-level anchor 취급은 타당")

---

## 4. Task 4 — 세 outcome colocalization 재분석 (가장 큰 작업)

### 문제

`TrackA_MR/v5_upgrade/scripts/taskE_01_coloc.R`의 실제 호출:

```r
coloc.abf(
  dataset1 = list(snp=m$SNP, pvalues=m$p_gwas,  type="cc",    s=s_val, N=n_val),
  dataset2 = list(snp=m$SNP, pvalues=m$p_eqtl,  type="quant", N=N_EQTL),
  MAF = m$MAF)     # <- g1000_eur_freq.frq, 두 dataset·두 outcome 모두에 동일 적용
```

세 가지 문제:

1. **beta/varbeta 미사용.** p-value 기반 ABF이므로 **MAF·N·s가 사후확률에 직접 반영**됩니다.
2. **단일 유럽 MAF 벡터**를 **동아시아 BBJ**에도 적용했습니다. PP.H4 = 0.951이 여기 의존합니다.
3. **UKB colocalization이 아예 없습니다.** `OUTCOMES` 리스트에 BBJ와 FinnGen뿐입니다
   (`TaskE_01_coloc_results_v1.csv` 확인 결과 outcome은 두 개).
   UKB는 eQTLGen과 인종이 가장 맞고 IGF1R MR이 가장 강한(P=0.012) outcome입니다.

부수적으로 `SNP.PP.H4`는 계산되지만(top SNP 선택에 사용) **보고되지 않습니다.**

### 해야 할 일

`taskE_01_coloc.R`을 다음 사양으로 개정해 `taskE_01b_coloc_v2.R`로 새로 만드십시오.
**기존 스크립트는 지우지 말고 남겨두십시오** (provenance 기록).

**(a) 입력 — 가능하면 beta/varbeta**

```r
# 우선순위 1: outcome 과 eQTL 모두 beta/se 가 있으면
dataset = list(snp=..., beta=..., varbeta=se^2, type="cc", s=..., N=...)
# 우선순위 2: beta 가 없으면 p-value + ancestry-matched MAF
```

- eQTLGen은 Z와 N을 주므로 β/SE를 유도할 수 있습니다
  (`taskD_03_mr_3outcomes.R`의 `prep_exposure()`가 이미 그렇게 합니다 — 그 로직 재사용).
- **BBJ에는 EAS MAF**를, UKB/FinnGen에는 EUR MAF를 쓰십시오.
  1000G Phase 3 EAS 패널이 `04_druggable_mr/data/`에 있는지 확인하고, 없으면
  BBJ GWAS 파일의 `effect_allele_frequency` 컬럼을 쓰십시오
  (`GCST90018627_harmonised.tsv.gz`에 있습니다).

**(b) outcome 3개 전부**

```r
OUTCOMES <- list(
  BBJ_Graves       = list(path=".../GCST90018627_harmonised.tsv.gz", s=2809/175465,  N=175465,  maf="EAS"),
  UKB_hyperthyroid = list(path=".../GCST90038636_harmonised.tsv.gz", s=3731/484598,  N=484598,  maf="EUR"),
  FinnGen_GO       = list(path="c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz",
                          s=<Task 2에서 확정>, N=<Task 2에서 확정>, maf="EUR")
)
```

> UKB의 `s`와 `N`은 `taskD_03_mr_3outcomes.R` L129의 `ukb_mu <- 3731/484598`에서
> 가져왔습니다. 이 값도 출처를 확인해 주십시오.

**(c) 보고 항목 — 전부**

유전자 × outcome 마다:
- `PP.H0 ~ PP.H4` 전부 (현재는 H0/H1도 저장은 되나 원고에 없음)
- **`SNP.PP.H4` 상위 변이와 그 값**
- **H4 조건부 95% credible set** (`res$results`를 `SNP.PP.H4` 내림차순 누적 0.95)
- 해당 locus에서 **outcome의 최소 P-value** (검출력 부족인지 진짜 null인지 구분용)
- 겹친 SNP 개수

**(d) prior sensitivity**

`p12 = 1e-6, 5e-6, 1e-5` 세 조건 전부 돌리고 PP.H4를 나란히 보고하십시오.
기본값 `p12=1e-5`는 조건에 따라 liberal할 수 있습니다 (Wallace 2020).

**(e) 대상 유전자**

backbone 3개 (`TSHR`, `IGF1R`, `CTLA4`) + 후보 5개
(`TNFSF14`, `IFNGR1`, `MAPKAPK5`, `HSD3B7`, `VKORC1`, `PRSS36`).

**(f) 가능하면 `coloc.susie`**

단일 인과변이 가정이 chr16p11.2처럼 복잡한 locus에서 무너집니다.
LD matrix를 만들 수 있으면 `coloc.susie`도 함께 돌려 비교하십시오.

### 산출물

- `TrackA_MR/v5_upgrade/05_coloc_candidates/TaskE_01b_coloc_v2_<date>.csv`
- `internal/COLOC_v2_report.md` — 위 항목 전부 + **v1 대비 PP.H4 변화표**

### 통과 기준

- 세 outcome 모두 결과가 있을 것
- `TSHR` BBJ PP.H4가 EAS MAF에서 얼마인지 명시 (v1의 0.951과 비교)
- `IGF1R` **UKB** PP.H4/PP.H2가 처음으로 보고될 것
- prior sensitivity 3조건에서 결론이 뒤집히는지 여부 명시

---

## 5. Task 5 — chr16p11.2 재검토 (Task 4 이후)

Table S5 각주가 스스로 `PRSS36`의 `rs78924645`가 나머지 둘과 **r² = 0.19 / 0.20**이라고
적어 놓고 "단일 LD 신호임을 확인한다"고 결론냅니다. **r²=0.19에서는 따라 나오지 않습니다.**

`rs4889606`(HSD3B7)–`rs34649473`(VKORC1)는 r²=0.90으로 하나의 신호가 맞지만,
`PRSS36`은 별개일 수 있습니다. Conditional analysis 또는 conditional colocalization으로
확인하십시오. 결과에 따라 **Figure 1 Panel A의 "3 — chr16p11.2 one LD block" 집계가
바뀝니다.**

---

## 6. 원고 표현 중 이미 확정된 수정 사항 (Task 2–4와 무관)

Task 결과와 상관없이 다음은 확정입니다. 다만 **원고 수정은 Task 2–4 통과 후에** 하십시오.

| 현재 | 수정 |
|---|---|
| "no robust novel target" | **"No novel candidate satisfied our stringent dual-phenotype colocalization criterion."** — `TNFSF14`는 BBJ coloc 0.994 + UKB β=−0.298, **P=0.0056** 방향 일치 복제. GD candidate requiring replication으로 남겨야 함 |
| IGF1R "no evidence of directional pleiotropy" | "no detectable evidence, with limited power" — instrument 3–4개뿐 |
| orbital RNA-seq를 TSHR의 근거층으로 사용 | descriptive illustration으로 강등, 초록·결론에서 `P = 0.032` 삭제 (control n=1) |
| "prespecified" | `TaskC_pre_analysis_plan_v1.md`에 등록·timestamp가 없으므로 **"defined analytical hierarchy"** 권장 |
| UKB LMM→log-odds 변환 미기술 | Methods에 명시: `β_logodds = β_LMM / (μ(1−μ))`, `μ = 3731/484598 ≈ 0.0077`, **β와 SE 모두 변환**, Z와 P는 보존 (`taskD_03_mr_3outcomes.R` L123–133) |
| Table S4 secondary-threshold 개수 | clumped인지 raw인지 명시 |

---

## 7. 로컬 세션에 붙여넣을 시작 프롬프트

```
이 저장소는 TED-TRAP Paper 1 (druggable-gene-wide MR + colocalization,
Graves disease / TED) 이다. 외부 리뷰에서 provenance 결함이 확인되어 제출을
보류한 상태다.

먼저 다음 두 문서를 읽어라:
  internal/SUBMISSION_HOLD_provenance_audit.md   (무엇이 왜 틀렸는지)
  internal/HANDOFF_local_reanalysis.md           (내가 할 일)
그리고 CLAUDE.md 의 프로젝트 규칙을 지켜라. 핵심은 숫자를 기억으로 채우지 말고
반드시 원본 파일에서 확인하는 것이다.

Step 1(instrument manifest)은 원격에서 끝났다. 나는 Task 2부터 시작한다.
Task 2 = FinnGen endpoint/증례수 확정. 이것부터 하나씩(하나씩) 진행하고,
각 Task 의 "통과 기준"을 충족했는지 명시적으로 보고해라.

원본 데이터는 c:/ProjectTEDGWAS/ 아래에 있다.
원고 문장·표·그림은 Task 2-4 가 전부 통과하기 전에는 절대 고치지 마라.
```

---

## 8. 되돌려 받아야 할 것

Task가 끝날 때마다 다음을 커밋하고 push해 주십시오. 원격 세션이 이어받아
원고·표·그림을 재생성합니다.

- `internal/FINNGEN_endpoint_provenance.md` (Task 2)
- `internal/SUSIE_rerun_report.md` + 재실행 PIP CSV (Task 3)
- `internal/COLOC_v2_report.md` + `TaskE_01b_coloc_v2_*.csv` (Task 4)
- 개정한 스크립트 (`taskE_01b_coloc_v2.R` 등) — **기존 스크립트는 남겨둘 것**

**커밋하지 말 것:** IRB 원자료(`data.txt`, FASTQ/BAM), TSHR-ATrap 특허 자료,
라이선스 데이터 원본(eQTLGen 전체 sumstats, FinnGen `.gz`, 1000G 패널).
`.gitignore`가 막고 있으니 약화시키지 마십시오.

---

## 9. 우선순위

```
Task 2 (FinnGen 확정)  ──┐
                         ├──> Task 4 (coloc 재분석) ──> Task 5 (chr16)
Task 3 (SuSiE 재실행)  ──┘
```

Task 2와 3은 서로 독립이라 병렬 가능합니다. **Task 4는 Task 2의 `s`/`N`이
확정되어야 시작할 수 있습니다.**

Task 4까지 통과하면 원고 재생성으로 넘어갑니다. 그 전에는 원고를 건드리지 않습니다.

---

# Task 6 (open) — eQTLGen allele-frequency sensitivity

The one item from external review round 7 that cannot be closed in the remote
container, because it needs the licence-bound eQTLGen allele-frequency file.

**The issue.** Exposure effect sizes were reconstructed from eQTLGen Z-scores and
per-SNP sample sizes using **1000 Genomes European allele frequencies (n = 503)**,
not the allele-frequency file eQTLGen distributes (n = 26,609). The manuscript now
states this and its consequence, but the check itself is missing.

**What is and is not affected** (this reasoning is now in Limitations, so the
re-analysis only has to confirm the magnitude):

- **Single-instrument loci** — β and SE both scale with the same factor, so β/SE,
  the direction and the *P* value are invariant. *TSHR* and the *CTLA4* discovery
  arm fall here, and their MR results cannot move.
- **Multi-instrument loci** — inverse-variance weights are 1/SE², so a per-SNP
  frequency error reweights the pooled estimate. *IGF1R* (4 instruments) is the
  case that matters.
- **Colocalization** — the v2 run uses beta and varbeta, both of which carry the
  frequency, so every posterior depends on it.
- **Table S8 power** — computed from the observed SEs, so it moves with them.

**What to run.**

1. Download the eQTLGen allele-frequency file
   (https://www.eqtlgen.org/cis-eqtls.html) and join it to the instrument manifest
   on rsID, checking allele orientation against the eQTLGen effect allele.
2. Re-derive β and SE for all 6,135 instruments with those frequencies.
3. Re-run: MR for the three backbone genes and the thirteen discovery hits across
   all three outcomes; colocalization for the nine genes in `taskE_01b_coloc_v2.R`;
   and Table S8.
4. Report, in `internal/EQTLGEN_MAF_SENSITIVITY.md`, a side-by-side table of β, SE,
   OR, *P* and PP.H4 under both frequency sources, plus the maximum absolute change
   in PP.H4 and whether any estimate crosses a reporting threshold.

**Acceptance.** If nothing crosses a threshold, add one sentence to Limitations
recording that the sensitivity analysis was performed and the conclusions were
unchanged, and the caveat can be shortened. If something does move, the affected
numbers are re-reported from the eQTLGen-frequency run, which is the better
primary analysis in any case.

**Note.** The instrument manifest to join against is
`TaskD_02a_eqtl_instruments_snp_level_v2_verified.csv` (6,135 rows), not the v1
file — v1 still contains the GAN row that failed the selection threshold.
