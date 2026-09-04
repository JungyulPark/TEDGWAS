# Task 2 — FinnGen endpoint / 파일 / 증례수 provenance 확정

**작성일:** 2026-09-04
**작업 머신:** `c:\ProjectTEDGWAS` (원본 데이터 보유)
**상태:** ✅ 통과 기준 3개 모두 충족

---

## 요약 (결론 먼저)

1. v5 colocalization이 사용한 파일은 **`c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz`** — FinnGen R12 **`GRAVES_OPHT`** (broader) endpoint. `taskE_01_coloc.R:20`에서 확정.
2. 그 endpoint의 공식 R12 GWAS 분석셋은 **858 cases / 499,490 controls / N = 500,348**.
   → **원고의 "858 cases"는 옳습니다.**
3. `taskE_01_coloc.R:76-77`의 **`N = 520,387`은 틀렸고 출처가 없습니다** (미치환 placeholder).
   올바른 값은 500,348. `s`는 0.00164877 → **0.00171481**.
4. 그러나 **PP.H4에 미치는 수치적 영향은 무시 가능** (절대차 ≈ 1.3×10⁻⁸). 이유는 §6.
   → **보고된 PP.H4 = 0.986은 이 결함만으로는 바뀌지 않습니다.**
5. 외부 리뷰어의 "786 / 894, 858은 어느 쪽과도 안 맞음" 주장은 **틀렸습니다.**
   786·894는 **Risteys** 수치(등록부/genotype-QC 단계)이고, 858은 **PheWeb = summary statistics
   실제 GWAS 분석셋** 수치입니다. 서로 다른 것을 셉니다. §5.
6. **별도의 실제 결함 2건**을 새로 발견했습니다 — `data_sources.csv`의 endpoint 라벨 오기와
   control 수 오기. §7.

---

## 1. 요구된 provenance 표

| 항목 | 값 | 출처 |
|---|---|---|
| 분석에 사용한 파일명 (전체 경로) | `c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz` | `TrackA_MR/v5_upgrade/scripts/taskE_01_coloc.R:20` (`FINNGEN <- ...`) |
| 파일 크기 | 798,476,004 bytes | `ls -la --time-style=full-iso` |
| 파일 mtime | 2026-04-23 13:13:38.928 +0900 | 동일 |
| 파일 MD5 | `1b26a77e1e4eb451e80ead00e4ec6f28` | `md5sum` (2026-09-04 실행) |
| endpoint code | `GRAVES_OPHT` ("Graves ophthalmopathy") | FinnGen R12 PheWeb `https://r12.finngen.fi/pheno/GRAVES_OPHT` — 페이지 헤더 "RELEASE 12 / Graves ophthalmopathy" |
| FinnGen release | R12 | 동일 (페이지 상단 "RELEASE 12"), 파일명과도 일치 |
| cases | **858** | R12 PheWeb 표시 "858 cases" |
| controls | **499,490** | R12 PheWeb 표시 "499490 controls" |
| total N | **500,348** | 858 + 499,490 = 500,348 |
| coloc 에 넣을 s | **0.001714810** | 858 / 500,348 (재계산) |
| 원고의 858이 무엇과 일치하는가 | **일치함** — R12 `GRAVES_OPHT` GWAS 분석셋 case 수 | §2·§3·§5 삼중 확인 |

**교차 검증:** FinnGen 공식 발표 "A total of 500,348 individuals, including 282,064 females and
218,284 males", "2,502 health endpoints", "> 21 M variants"
(https://www.finngen.fi/en/results-based-full-finngen-cohort-500000-participants-released) —
R12 전체 코호트 크기가 500,348로 독립 확인되며, 858 + 499,490과 정확히 일치합니다.

---

## 2. 파일 내부 증거 — 외부 문서에 의존하지 않는 s 유도

세 FinnGen `.gz` 파일 모두 `af_alt`, `af_alt_cases`, `af_alt_controls` 컬럼을 갖고 있습니다
(헤더: `#chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt af_alt_cases af_alt_controls`).

전체 대립유전자빈도는 case/control 빈도의 표본크기 가중평균이므로:

```
af_alt = s · af_alt_cases + (1 − s) · af_alt_controls
     ⇒  s = (af_alt − af_alt_controls) / (af_alt_cases − af_alt_controls)
```

각 파일 앞 400,000행에서 `af_alt_cases − af_alt_controls > 0.0005` 이고 `af_alt < 0.05`
(반올림 오차 대비 정보량이 큰 행)인 SNP만 사용해 SNP별 `s`를 계산했습니다.

| 파일 | 정보성 SNP 수 | median s | Q1 | Q3 |
|---|---|---|---|---|
| `finngen_R12_GRAVES_OPHT.gz` | 71,159 | **0.00171482** | 0.00171268 | 0.00171700 |
| `finngen_R12_E4_GRAVES_OPHT_STRICT.gz` | 72,687 | **0.00150495** | 0.00150306 | 0.00150683 |
| `finngen_R12_E4_GRAVES_STRICT.gz` | 43,766 | **0.00791852** | 0.00791307 | 0.00792393 |

IQR이 중앙값의 약 0.25%이고 n이 크므로 SEM ≈ 1×10⁻⁸ 수준입니다. 분포가 이렇게 좁다는 사실
자체가 `af_alt`가 실제로 해당 endpoint의 case/control 가중평균이라는 내적 검증입니다.

**N = 500,348을 대입하면 case 수가 정수로 떨어집니다:**

| 파일 | s × 500,348 | 반올림 |
|---|---|---|
| `GRAVES_OPHT` | 858.01 | **858** |
| `E4_GRAVES_OPHT_STRICT` | 752.98 | **753** |
| `E4_GRAVES_STRICT` | 3962.08 | **3962** |

---

## 3. 삼중 독립 확인

세 개의 서로 독립적인 출처가 동일한 값을 지시합니다.

| endpoint | ① 파일내부 유도 (§2) | ② `03b_04v3_finngen_local.R:51-53` | ③ R12 PheWeb (공식) |
|---|---|---|---|
| `GRAVES_OPHT` | 858 / 500,348 | `N=500348, cases=858, controls=499490` | **858 cases / 499490 controls** |
| `E4_GRAVES_OPHT_STRICT` | 753 / 500,348 | `N=500348, cases=753, controls=499595` | **753 cases / 499595 controls** |
| `E4_GRAVES_STRICT` | 3,962 / 500,348 | `N=500348, cases=3962, controls=496386` | **3962 cases / 496386 controls** |

①과 ②의 편차: −0.0014% / +0.0022% / −0.0020%. ③과는 정의상 완전 일치.

②는 v4 시절 스크립트 `TrackA_MR/scripts/03b_04v3_finngen_local.R`에 이미 올바른 메타데이터가
기록되어 있었다는 뜻입니다. **정답은 저장소 안에 이미 있었고, v5의 `taskE_01_coloc.R`이
그것을 참조하지 않고 placeholder를 새로 썼습니다.**

---

## 4. 실제 투입값 검정

측정 기준값: `GRAVES_OPHT` s = 0.00171483 (§2)

| 위치 | 투입한 s | 편차 | 판정 |
|---|---|---|---|
| `v5_upgrade/scripts/taskE_01_coloc.R:76` `858/520387` | 0.00164877 | **−3.852%** | ❌ **불일치** |
| `v5_upgrade/scripts/taskE_02_coloc_candidates.R:76` `858/520387` | 0.00164877 | **−3.852%** | ❌ **불일치** |
| `scripts/06_coloc_all_loci.R:139` `858/500348` | 0.00171481 | −0.001% | ✅ 일치 |
| `00_meta/data_sources.csv` `858 / 376419 controls` | 0.00227419 | **+32.619%** | ❌ **불일치** |

### 4.1 `520,387`의 출처 — 없음

저장소 전체에서 `520387`이 등장하는 곳은 네 군데뿐이며, 모두 **사용 기록**이지 **출처**가 아닙니다:

- `taskE_01_coloc.R:76,77`
- `taskE_02_coloc_candidates.R:76,77`
- `00_meta/analytical_log.md:142-144`

로그 원문:

```
## Note on FinnGen GO Coloc Parameters
- Total N used: 520,387
- Case count used: 858
- Computed s (case proportion): 858 / 520,387 = 0.001648
```

출처 인용이 없습니다. 스크립트의 주석은 이 값이 미확정임을 명시하고 있었습니다:

```r
# FinnGen s/N: CONFIRM control count. R12 GO ~858 cases. Total from file rows or meta.
# Placeholder — Gemini sets these from FinnGen R12 GRAVES_OPHT documentation:
FINNGEN_S <- 858 / 520387
FINNGEN_N <- 520387
```

`load_finngen()`도 `n_case=NA, n_ctrl=NA`를 반환하며 주석에 `# Gemini: fill s from FinnGen meta`
라고 되어 있습니다 (`taskE_01_coloc.R:64-68`). **치환이 끝내 이루어지지 않은 채 결과가 확정되었습니다.**

520,387이 어떤 실제 FinnGen 수치와 대응하는지는 **미확인**입니다. R12 코호트(500,348), R12
`GRAVES_OPHT` 총계(500,348), Risteys 등록부 수치 어느 것과도 일치하지 않습니다.

### 4.2 `scripts/06_coloc_all_loci.R`은 v5 결과의 출처가 아님

이 스크립트는 `858/500348`로 **올바른** 값을 갖고 있으나,
`finngen_file <- "TrackA_MR/data/finngen/finngen_R12_E4_GRAVES_OPHT.tsv.gz"` (L45)를 참조하는데
**`TrackA_MR/data/finngen/` 디렉토리가 존재하지 않습니다.** v3 시절 미실행 스크립트이며 v5
결과와 무관합니다.

---

## 5. 외부 리뷰어 주장 검증 — 리뷰어가 틀렸습니다

> 리뷰어 주장: "R12 공식 수치는 `E4_GRAVES_OPHT_STRICT` = 786, `GRAVES_OPHT` = 894이고
> 858은 어느 쪽과도 맞지 않는다"

**FinnGen은 같은 endpoint·같은 release에 대해 서로 다른 단계의 수치를 두 사이트에 게시합니다.**

| endpoint | Risteys R12 | PheWeb R12 (= summary stats) | 파일내부 유도 |
|---|---|---|---|
| `GRAVES_OPHT` | 894 ("Number of individuals", **FG R12 + FinRegistry** 등록부 합산) | **858 cases** | 858 |
| `E4_GRAVES_OPHT_STRICT` | 786 ("FinnGen genotyped cases after filtering based on genotype QC") | **753 cases** | 753 |

- 리뷰어의 **894**는 Risteys의 FinnGen+FinRegistry **등록부 개인 수**입니다. GWAS 분석셋이 아닙니다.
- 리뷰어의 **786**은 Risteys의 **genotype-QC 후 endpoint 정의 case 수**이며, 그마저도
  **다른 endpoint**(strict)의 것입니다. GWAS 분석셋(753)과도 다릅니다.
- 우리 분석이 쓴 **858**은 **summary statistics 파일이 실제로 만들어진 분석셋**의 case 수이고,
  이는 β·SE·p-value가 나온 바로 그 표본입니다. **coloc/MR에 넣어야 할 값은 이쪽이 맞습니다.**
- 결정적으로, 858은 **파일 내부 대립유전자빈도에서 독립적으로 재현**됩니다(§2). 외부 문서
  해석에 의존하지 않습니다.

리뷰어 수치를 받아들이려면 N이 각각 521,334 / 522,287이어야 하는데, 이는 R12 코호트
500,348보다 +4.2% / +4.4% 크며 파일 내용과 모순됩니다.

**→ 원고 본문·Table 1의 "858 cases"는 수정 대상이 아닙니다. 유지합니다.**

---

## 6. 수치적 영향 — PP.H4는 사실상 변하지 않습니다

설치된 `coloc` 5.2.3 소스에서 직접 확인:

```r
> coloc:::Var.data.cc
function (f, N, s)
{
    1/(2 * N * f * (1 - f) * s * (1 - s))
}
```

`process.dataset`은 `pvalues + MAF + N`이 주어지면 `approx.bf.p()`를 거쳐 이 함수를 씁니다.
**`N`과 `s`는 곱 `N·s·(1−s)`로만 들어갑니다.**

두 파라미터화 모두 동일한 case 수 858을 인코딩하고 있으므로:

| | N·s·(1−s) |
|---|---|
| 사용값 (520,387) | 856.5854 |
| 올바른 값 (500,348) | 856.5287 |
| 비 | 1.00006615 → **V 차이 0.0066%** |

**경험적 확인** — 합성 공유-인과변이 locus(1,500 SNP)에 실제 `coloc.abf`를 두 파라미터화로 실행:

| | used (520,387) | correct (500,348) | 절대차 |
|---|---|---|---|
| PP.H3 | 0.005389675 | 0.005389688 | +1.25×10⁻⁸ |
| **PP.H4** | **0.9946101916** | **0.9946101790** | **−1.25×10⁻⁸** |

재현 스크립트: `internal/scripts/verify_coloc_sN_impact.R`

**결론:** 이 결함은 **provenance 결함이지 수치 결함이 아닙니다.** 보고된 PP.H4
(TSHR FinnGen 0.986 등)는 이 오류만으로는 소수 3자리에서 변하지 않습니다.
다만 **기록·재현성 측면에서는 반드시 수정해야 합니다** — 논문이 명시하지 않은 미확인 값이
계산에 들어가 있었다는 사실 자체가 문제입니다.

⚠️ 이것은 Task 4의 **다른** 결함들(p-value 기반 ABF, EAS BBJ에 EUR MAF 적용, UKB coloc 부재)에
대해서는 **아무것도 말해주지 않습니다.** 그쪽은 여전히 실질적 영향이 있을 수 있으며 Task 4에서
별도로 측정해야 합니다.

---

## 7. 새로 발견된 결함 2건 — `00_meta/data_sources.csv`

`FinnGen_R12_GO_v1` 행에 두 개의 독립적인 오류가 있습니다.

| 필드 | 현재 기재 | 실제 | 판정 |
|---|---|---|---|
| `local_path` | `c:\ProjectTEDGWAS\finngen_R12_GRAVES_OPHT.gz` | 동일 | ✅ |
| `citation` | "FinnGen R12 endpoint **E4_GRAVES_OPHT_STRICT**" | 파일은 **`GRAVES_OPHT`** (broader) | ❌ **endpoint 라벨 오기** |
| `notes` | "858 cases / **376419** controls" | 858 / **499,490** | ❌ **control 수 오기 (−25%)** |

즉 이 CSV는 **broader 파일을 가리키면서 strict endpoint 이름을 붙이고, 어느 쪽과도 맞지 않는
control 수를 적어 두었습니다.** 858 + 376,419 = 377,277이며 이는 s를 32.6% 왜곡합니다.
이 값이 계산에 쓰인 흔적은 없으나(§4), 데이터 출처 표로서 논문 재현성 문서에 인용될 수 있으므로
수정 대상입니다.

---

## 8. 부수 확인 사항

- **MR은 영향 없음.** `taskD_03_mr_3outcomes_v1.R:50`은 `n_case = 858`을 메타데이터로만 보유하며,
  MR 추정은 β/SE로만 이루어집니다. s/N은 MR에 들어가지 않습니다.
- **파일 mtime(2026-04-23) < Task E.1 실행 시각(2026-05-26 16:33)** — 실행 후 파일이 교체된
  정황은 없습니다.
- **endpoint 선택 자체는 타당** — 원고는 FinnGen 결과를 "TED-specific sensitivity outcome"으로
  기술하며, broader `GRAVES_OPHT`(858 cases)를 쓴 것은 strict(753)보다 검정력이 높은 합리적 선택입니다.
  다만 **어느 endpoint를 썼는지 원고가 명시하지 않습니다.** Methods에 endpoint code를 적어야 합니다.

---

## 9. 통과 기준 대조

| 기준 | 결과 |
|---|---|
| 파일 경로·endpoint·release·cases·controls가 **하나의 출처로** 확정 | ✅ R12 PheWeb `GRAVES_OPHT` = 858 / 499,490, release 12. 추가로 파일내부 유도·v4 스크립트 기록이 독립 확인 |
| `s = cases/(cases+controls)`가 그 값에서 계산됨 | ✅ 858 / 500,348 = **0.001714810** |
| `858`과 불일치하면 명시 | ✅ **858은 일치합니다.** 불일치한 것은 **분모** — 520,387은 출처 없는 placeholder이며 올바른 값은 500,348 |

---

## 10. Task 4로 넘길 확정 파라미터

```r
FinnGen_GO = list(
  file     = "c:/ProjectTEDGWAS/finngen_R12_GRAVES_OPHT.gz",
  md5      = "1b26a77e1e4eb451e80ead00e4ec6f28",
  endpoint = "GRAVES_OPHT",      # FinnGen R12, "Graves ophthalmopathy"
  release  = "R12",
  cases    = 858,
  controls = 499490,
  N        = 500348,
  s        = 858/500348,         # 0.001714810
  maf      = "EUR"
)
```

참고 — 같은 파일군의 다른 두 endpoint (Task 4에서 민감도 분석에 쓸 경우):

```r
E4_GRAVES_OPHT_STRICT: cases=753,  controls=499595, N=500348, s=753/500348
E4_GRAVES_STRICT:      cases=3962, controls=496386, N=500348, s=3962/500348
```

---

## 11. 수정 대상 목록 (Task 4 통과 후 일괄 적용)

원고는 아직 건드리지 않았습니다. 아래는 **기록만** 해둡니다.

| # | 대상 | 수정 내용 | 근거 |
|---|---|---|---|
| 1 | `taskE_01_coloc.R:76-77` → 개정판 `taskE_01b_coloc_v2.R` | `520387` → `500348` | §1, §4 |
| 2 | `taskE_02_coloc_candidates.R:76-77` → 개정판 | 동일 | §4 |
| 3 | `00_meta/data_sources.csv` `FinnGen_R12_GO_v1` | citation의 endpoint를 `GRAVES_OPHT`로, notes의 controls를 `499490`으로 | §7 |
| 4 | `00_meta/analytical_log.md:141-144` | 정정 노트 추가 (기존 기록은 삭제하지 않음) | §4.1 |
| 5 | 원고 Methods | FinnGen endpoint code(`GRAVES_OPHT`)와 controls/N 명시 | §8 |
| 6 | 원고 본문·Table 1 "858 cases" | **변경 없음 — 옳습니다** | §5 |

기존 스크립트는 지우지 않았습니다. 개정은 Task 4에서 `taskE_01b_coloc_v2.R` 신규 파일로 수행합니다.

---

## 12. 재현 방법

```bash
# s 유도 (파일내부 증거)
bash internal/scripts/derive_finngen_s.sh

# coloc s/N 영향 검증
Rscript internal/scripts/verify_coloc_sN_impact.R
```

공식 수치 확인처:
- https://r12.finngen.fi/pheno/GRAVES_OPHT (RELEASE 12 / 858 cases / 499490 controls)
- https://r12.finngen.fi/pheno/E4_GRAVES_OPHT_STRICT (753 cases / 499595 controls)
- https://r12.finngen.fi/pheno/E4_GRAVES_STRICT (3962 cases / 496386 controls)
- https://www.finngen.fi/en/results-based-full-finngen-cohort-500000-participants-released (R12 = 500,348)
- https://r12.risteys.finngen.fi/endpoints/GRAVES_OPHT (894 = FG R12 + FinRegistry 등록부)
- https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_OPHT_STRICT (786 = genotype-QC 후)
