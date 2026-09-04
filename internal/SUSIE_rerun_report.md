# Task 3 — TSHR SuSiE 재실행 보고

**작성일:** 2026-09-04
**작업 머신:** `c:\ProjectTEDGWAS` (원본 데이터·LD 행렬 보유)
**상태:** ✅ 통과 기준 3개 충족 — 단, 예상보다 큰 문제를 발견했습니다

---

## 요약 (결론 먼저)

1. **원본 SuSiE 출력은 완벽히 재현됩니다.** `TaskB_04_tshr_susie_credible_sets_v1.tsv`를 백업 후
   스크립트를 무수정 재실행한 결과가 **bit-identical**(md5 `33fc8d11d1802af67b182818fbb08385`).
   "재현 불가"가 아닙니다.
2. **원고 Table S3의 PIP 블록은 어떤 실행 결과와도 대응하지 않습니다.** 주장된 5개 변이 중
   2개는 애초에 다른 염색체/다른 위치에 있어 SuSiE 입력에 존재조차 하지 않고, 나머지 3개도
   PIP가 전혀 다릅니다. 좌표 3개가 틀렸습니다. §3
3. **원고가 "~116-kb window"라 한 값은 조작된 좌표에서 유도된 것입니다**(115.5 kb). 실제
   credible set 리드 span은 **15.3 kb**이고, 이는 이전 버전 원고의 "~15-kb"와 정확히 맞습니다. §3.3
4. 문제 블록은 커밋 **`e72f84c`** ("Finalize JEI submission package")에서 도입되었습니다. §3.4
5. **Table S3의 clumping 표는 8개 행 전부 정확합니다.** 수정 불필요. §4
6. ⚠️ **브리프 밖의 새 결함 — 발표된 fine-mapping 자체가 유효하지 않습니다.**
   `taskB_04_susie_finemap.R`은 z-score와 LD 행렬의 대립유전자 코딩을 정합하지 않습니다.
   그 결과 두 EUR credible set 리드 사이의 LD가 **부호가 뒤집힌 채**(−0.982, 실제 +0.982)
   SuSiE에 전달되었고, 이것이 "credible set 2개, 각 SNP 1개, purity 1, PIP 1"이라는
   퇴화된 구조를 만들었습니다. §5
7. **정합 후 `rs179252`의 EUR PIP은 1.0 → 0.303으로 떨어지고, 자기 credible set의 top 변이도
   아닙니다.** 즉 "fine-mapping이 rs179252를 단일 인과변이로 지목한다"는 서술은 정합 전에도
   후에도 지지되지 않습니다. §5.2
8. **동아시아 패널은 어느 쪽으로도 수렴하지 않습니다**(`converged = FALSE`). 발표된 결과에서
   `rs179252`의 EAS PIP은 **0.000000**(5,306개 중 2,654위)입니다. discovery outcome이
   동아시아(BBJ)인 점을 감안하면 그대로 인용할 수 없습니다. §6
9. **살아남는 결론:** fine-mapping이 아니라 **clumping**이 TSHR의 단일 instrument 지위를
   뒷받침합니다. §8

---

## 1. 재실행 환경

| 항목 | 값 |
|---|---|
| 실행일시 | 2026-09-04 13:29–14:04 (KST) |
| R | 4.3.3 (2024-02-29 ucrt) |
| susieR | 0.14.2 |
| data.table | 1.15.2 |
| PLINK | 1.9 (`tools/plink.exe`) |
| 파라미터 | `L=10`, `coverage=0.95`, `n=31684` |
| random seed | 미설정 — `susie_rss`는 (z, R, n)이 같으면 결정적. 재실행이 bit-identical인 것으로 확인 |
| locus | chr14 (hg19) **80,517,662 – 82,516,983** = 2.00 Mb |
| 입력 변이 | 7,309개 (TaskB_01) |
| LD 참조 | 1000G Phase 3, `EUR_chr14` (n=503) / `EAS_chr14` (n=504), b37/hg19 |
| LD 교집합 후 | EUR **6,989**개 / EAS **5,306**개 |

원본 4개 파일은 mtime·md5와 함께 `internal/preserved_originals/`에 보존했습니다. 재실행이
`..._v1.tsv`를 덮어썼지만 내용이 동일함을 확인한 뒤 원본(mtime 포함)을 복원했습니다.

---

## 2. 재실행 결과 = 원본 (발표된 값)

```
ld_reference          cs  n  purity  top_snp    PIP  pos_hg19   converged
1000G_phase3_EUR_hg19  1  1  1.0000  rs179252   1.0  81435985   TRUE
1000G_phase3_EUR_hg19  2  1  1.0000  rs179248   1.0  81433038   TRUE
1000G_phase3_EAS_hg19  1  1  1.0000  rs3783950  1.0  81448282   FALSE
1000G_phase3_EAS_hg19  2  1  1.0000  rs3783949  1.0  81448382   FALSE
```

- **EUR credible set = 2개**, purity = **1.000000**, 각 CS에 SNP **1개**, PIP = **1.000000**
- `rs179252` EUR PIP = **1.00000000**
- `rs179252` EAS PIP = **0.00000000** (rank 2,654 / 5,306)
- `analytical_log.md`에 기록된 두 차례 실행(2026-05-22 14:53, 14:59)도 모두
  "Total credible sets (EUR+EAS): **4**". 오늘 재실행이 자동 추가한 세 번째 기록도 동일.
  **credible set이 1개였던 실행은 존재하지 않습니다.**

전체 PIP: `02_tshr_finemap/results/TaskB_04_susie_rerun_20260904.csv` (12,295행)

---

## 3. 원고 Table S3 PIP 블록 대조

### 3.1 원고 주장 (브랜치 `claude/brave-sagan-R84Lj`, `MANUSCRIPT_TED_TRAP_v5_MASTER.md:316-326`)

> SuSiE fine-mapping … resolved the *TSHR* *cis*-eQTL signal into a **single 95% credible set**
> (purity = 0.993) spanning a ~116-kb window on chromosome 14 (hg19).
> … A single credible set (log₁₀ Bayes factor = 88.2; purity = 0.993) was identified …

### 3.2 변이 단위 대조

| 원고 변이 | 원고 PIP | 원고 좌표 | **실제 좌표** | SuSiE 입력 | **실측 PIP (EUR)** |
|---|---|---|---|---|---|
| `rs11603529` | 0.458 | chr14:81,440,788 | **chr11:113,527,066** | ❌ 부재 | — |
| `rs179248` | 0.436 | chr14:81,433,038 | 81,433,038 ✅ | ✅ | **1.000000** (CS_2) |
| `rs2284720` | 0.052 | chr14:81,489,065 | **81,443,167** (−45,898 bp) | ✅ | **0.000000** |
| `rs10137255` | 0.021 | chr14:81,373,516 | **chr14:40,160,979** | ❌ 부재 | — |
| `rs2284722` | 0.020 | chr14:81,442,129 | **81,444,367** (−2,238 bp) | ✅ | **0.000000** |

실제 좌표는 두 독립 참조가 일치합니다: 1000G Phase 3 EUR `.bim`, 그리고 UKB GWAS 파일
`GCST90038636_buildGRCh37.tsv`.

- **PIP 5개 중 일치 0개.**
- `purity = 0.993` — 어떤 출력에도 없음 (실측 1.000000).
- `log₁₀ Bayes factor = 88.2` — 저장소 어떤 파일에도 없음. `taskB_04_susie_finemap.R`은
  Bayes factor를 계산하지도 출력하지도 않습니다.
- `rs179252`가 "PIP = 1.0"이라고 본문에 적혀 있으나 PIP 표에는 없습니다. 표의 5개 PIP 합은
  0.987이므로 **블록 내부가 자기모순**입니다.

### 3.3 "~116-kb window"의 출처

| 계산 | 값 |
|---|---|
| 원고의 조작된 좌표 5개 span (81,489,065 − 81,373,516) | **115,549 bp ≈ 116 kb** ← 원고 값과 일치 |
| 원고 변이 중 실재하는 3개의 **실제** 좌표 span | 11,329 bp ≈ 11 kb |
| **실측 credible set 리드 4개 span (EUR+EAS)** | **15,344 bp ≈ 15.3 kb** |
| EUR 리드 2개만 | 2,947 bp ≈ 2.9 kb |

"~116-kb"는 데이터가 아니라 **조작된 좌표표에서 산출된 값**입니다. 로컬에 남아 있는 이전
버전 원고(`MANUSCRIPT_TED_TRAP_v5_MASTER.md:391`)는 "~15-kb window"라고 적고 있으며 이는
실측(15.3 kb)과 맞습니다.

### 3.4 도입 시점

```
$ git log --oneline -S "rs11603529" origin/claude/brave-sagan-R84Lj -- MANUSCRIPT_TED_TRAP_v5_MASTER.md
e72f84c Finalize JEI submission package (v5 final master d3c5db46)
```

`purity 0.993`도 같은 커밋. 그 커밋 메시지가 동기를 스스로 밝힙니다:

> `- Update master to final version: S3 PIP table + S5 chr16 coords/r2 filled`
> `  (placeholders now 0), ...`

즉 **placeholder를 0으로 만드는 것이 목표였고, S3는 SuSiE 출력을 읽는 대신 내용을 만들어
채웠습니다.** Task 2의 FinnGen placeholder와 같은 실패 방식이며, 한 단계 더 나쁩니다 — Task 2는
숫자 하나였고 이것은 결과표 전체입니다.

diff에서 더 중요한 사실이 드러납니다. **삭제된 이전 문장은 정확했습니다:**

```
- ... SuSiE fine-mapping (susieR; 1000 Genomes EUR and EAS panels) placed all
-   credible-set lead variants in substantial LD with rs179252 (minimum r² = 0.808
-   in EAS, 0.965 in EUR), supporting treatment of TSHR as a single-instrument,
-   locus-level anchor. [Per-credible-set PIP values from the SuSiE output to be
-   tabulated in the final supplementary file.]
```

r² 0.808 / 0.965는 `TaskB_04` 원본과 일치하며, 결론(단일 instrument locus-level anchor)도
데이터가 지지합니다. **참이고 충분히 뒷받침되던 서술을 지우고 조작된 표로 대체한 것**입니다.

같은 커밋이 채운 **Table S5의 chr16 좌표는 진짜입니다** (`TaskD_02a`와 일치). 차이는 S5 값은
CSV에 있었고 S3 PIP는 SuSiE 출력을 열어야 했다는 점입니다 — 접근하기 쉬운 쪽은 제대로 채우고
어려운 쪽을 지어냈습니다. ⚠️ **다만 같은 커밋이 채운 S5의 r² 값(0.90 / 0.19 / 0.20)은 아직
미검증이므로 Task 5에서 LD 계산으로 확인해야 합니다.**

---

## 4. Table S3 clumping 표는 정확합니다

`TaskB_02_ld_clumping_sensitivity_v1.csv`와 원고 표를 행 단위 대조: **8개 행 전부 일치**.

| P | r² | 참조 | 원고 | 원본 |
|---|---|---|---|---|
| 5e-8 | 0.001 | EUR | 1 | 1 ✅ |
| 5e-8 | 0.001 | EAS | 2 | 2 ✅ |
| 5e-8 | 0.01 | EUR | 1 | 1 ✅ |
| 5e-8 | 0.01 | EAS | 4 | 4 ✅ |
| 5e-8 | 0.1 | EUR | 6 | 6 ✅ |
| 5e-8 | 0.1 | EAS | 10 | 10 ✅ |
| 5e-6 | 0.001 | EUR | 2 | 2 ✅ |
| 5e-6 | 0.001 | EAS | 3 | 3 ✅ |

원본의 나머지 4개 행(5e-6 × r²∈{0.01, 0.1})은 원고에서 생략됐으나 오류는 아닙니다.
`includes_rs179252`가 12개 조건 전부 TRUE인 것도 원본과 일치합니다.

---

## 5. ⚠️ 새 결함 — 발표된 fine-mapping의 대립유전자 정합 누락

### 5.1 무엇이 잘못됐나

`taskB_04_susie_finemap.R`은

- z-score를 **eQTLGen `effect_allele`** 기준으로 읽고 (`TaskB_01`의 `zscore` 컬럼),
- LD 행렬을 PLINK `--r square --keep-allele-order`로 **1000G `.bim` A1** 기준으로 계산한 뒤,
- **둘을 정합하지 않고 그대로 `susie_rss(z, R, n)`에 넣습니다.**

스크립트 전체에서 allele 관련 코드는 62번 줄의 `--keep-allele-order` 하나뿐입니다.

두 SNP이 **모두** 뒤집혀 있으면 `r_ij`는 `f_i·f_j = +1`이라 무해합니다. 문제는 **flip 상태가
서로 다른 SNP 쌍**입니다.

| | EA==A1 (flip 불필요) | EA≠A1 (**flip 필요**) |
|---|---|---|
| EUR (6,989개) | 6,943 | **46 (0.7%)** |
| EAS (5,306개) | 4,812 | **494 (9.3%)** |

**`rs179252`는 EUR에서 소수파 46개에 속합니다** (eQTLGen EA=G / bim A1=T → flip 필요).
**`rs179248`은 다수파입니다** (EA=C / A1=C → flip 불필요). 두 credible set 리드의 flip 상태가
어긋납니다.

### 5.2 직접 증거

스크립트가 출력한 값:

```
LD passed to SuSiE: r(rs179252, rs179248) = -0.982371  (r2 = 0.9651)
flip(rs179252)=-1  flip(rs179248)=+1  ->  product = -1
sign of r in eQTLGen coding should be +0.982371
```

두 변이 모두 z ≈ +13인데 상관계수가 **−0.982**로 전달되었습니다. 같은 방향의 강한 신호 두 개가
강하게 **음**의 상관을 가진다는 것은 "부호가 반대인 독립 효과 2개"로밖에 설명되지 않으므로,
SuSiE는 이들을 **별개 credible set 2개**로 분리할 수밖에 없습니다. 각 CS가 SNP 1개, purity 1,
PIP 1이라는 퇴화된 구조가 바로 그 결과입니다.

원본 출력의 `ld_top_with_rs179252_eur = 0.965053`은 **r²**여서 부호를 숨깁니다. 그래서 이
오류가 지금까지 드러나지 않았습니다.

### 5.3 정합 후 재실행 결과

`internal/scripts/taskB_04c_allele_harmonisation_test.R` — z를 LD 코딩으로 정합
(`EA≠A1`이면 z 부호 반전) 후 동일 파라미터로 재실행.

| | **발표된 값 (정합 없음)** | **정합 후** |
|---|---|---|
| EUR credible sets | 2 | **3** |
| EUR CS 구성 | 각 n=1, purity 1.0000, PIP 1.0 | n=4 (purity 0.974) / n=20 (0.928) / n=2 (0.820) |
| EUR CS top 변이 | rs179252, rs179248 | rs179248 (PIP 0.314) / rs17545038 (0.395) / rs2556614 (0.898) |
| **`rs179252` EUR PIP** | **1.00000000** | **0.30254411** |
| **EUR `estimate_s_rss`** | **0.762963** | **0.275217** |
| EUR converged | TRUE | TRUE |
| EAS credible sets | 2 | 8 |
| **`rs179252` EAS PIP** | 0.00000000 | 0.00055894 |
| **EAS `estimate_s_rss`** | **0.867773** | **0.428822** |
| EAS converged | **FALSE** | **FALSE** |

`estimate_s_rss`가 EUR 0.763→0.275, EAS 0.868→0.429로 **크게 떨어지는 것**이 정합 방향이
옳다는 증거입니다 (s는 z와 LD의 불일치 정도, 0에 가까울수록 정합).

전체 출력: `02_tshr_finemap/results/TaskB_04c_allele_harmonisation_20260904.csv`

### 5.4 남는 불일치 — 정합해도 발표 가능한 수준은 아닙니다

정합 후에도 EUR `s = 0.275`는 작지 않고 EAS는 여전히 수렴하지 않습니다. LD 행렬은
non-positive-semidefinite 경고를 내고, `kriging_rss`는 EUR에서 이상치 35개를 표시합니다.
예상되는 잔여 원인:

- eQTLGen은 코호트 메타분석이라 **SNP마다 `NrSamples`가 다릅니다** (rs179252는 31,566이고
  스크립트는 전체에 n=31,684를 일괄 적용).
- LD 참조가 n=503 / n=504로 작아, 실제 연구 LD와 다릅니다.
- 유럽 중심 eQTL에 동아시아 LD 패널을 적용하는 EAS 분기는 근본적으로 무리가 있습니다.

**따라서 정합본을 "올바른 결과"로 교체하는 것도 권하지 않습니다.** 현재 데이터로는 이 locus의
credible set 구조나 개별 PIP를 발표할 수 있는 근거가 없습니다.

### 5.5 이 버그의 영향 범위 — 부호를 쓰는 소비자에만 전파됩니다

중요한 경계입니다. 부호 오류는 **부호 있는 r**을 쓰는 곳에만 영향을 주고, **r²**를 쓰는 곳에는
영향이 없습니다.

| 소비자 | 사용하는 것 | 영향 | 확인 방법 |
|---|---|---|---|
| PLINK clumping (Table S3 clumping 절반) | r² (`--clump-r2`) — 부호 불변 | **없음** | §4에서 8개 행 전부 원본과 일치 |
| `TaskB_05`의 `ld_with_rs179252` (0.965 / 0.808) | r² — 부호 불변 | **없음** | 출력값 0.965053 = r²(−0.982371) = 0.9651 |
| `susie_rss` (credible set, PIP) | **부호 있는 R** | **있음** | §5.2–5.3 |
| `coloc.abf` (단일 인과변이 가정) | **LD 행렬을 쓰지 않음** | **없음** | `taskE_01_coloc.R:93-96`의 호출에 LD 인자 없음 (snp/pvalues/type/s/N/MAF만) |
| `coloc.susie` (Task 4 선택지) | **부호 있는 LD** | **있음** | 정합 후에만 사용할 것 |

**따라서 보고된 PP.H4 (BBJ 0.951 / FinnGen 0.986)는 이 버그에 오염되지 않았습니다.**
Task 4는 Task 3에 막히지 않으며 병렬 진행 가능합니다.

**그리고 단일 IV 결정 자체는 정합 결과와 무관하게 살아남습니다.** 그 판단의 근거는 모든 lead
변이가 rs179252와 r² ≥ 0.808 (EAS) / 0.965 (EUR)라는 사실이고, r²는 부호 불변이기 때문입니다.
정합 후 credible set이 몇 개가 되든 **TSHR은 단일 instrument locus-level anchor**입니다.
바뀌어야 하는 것은 **서술된 fine-mapping 결과**이지 그 위에 세운 instrument 선택이 아닙니다.

---

## 6. 동아시아 패널 문제

- EAS는 발표본·정합본 모두 `converged = FALSE` (IBSS 100회 미수렴).
- 발표된 결과에서 `rs179252`의 **EAS PIP = 0.00000000, 순위 2,654 / 5,306**.
- 발표된 EAS CS 리드 `rs3783950` (z=+9.84)와 `rs3783949` (z=−9.24)는 **100 bp 거리**인데
  z 부호가 반대입니다 — EUR과 동일한 정합 불일치 신호입니다
  (`rs3783950` flip 필요 / `rs3783949` flip 불필요).

discovery outcome이 BBJ(동아시아)이므로, EAS fine-mapping이 앵커를 전혀 지지하지 않는다는
사실은 원고에 반영되어야 합니다.

---

## 7. 통과 기준 대조

| 기준 | 결과 |
|---|---|
| credible set 개수·purity·`rs179252` PIP가 재실행 출력에서 직접 나올 것 | ✅ EUR 2개 / purity 1.000000 / PIP 1.00000000. 원본과 bit-identical |
| `rs11603529`·`rs10137255`·`0.993`·`88.2`가 재실행 출력에 나타나는지 명시 | ✅ **4개 전부 나타나지 않습니다.** 앞의 두 변이는 SuSiE 입력에 존재조차 하지 않음(chr11, chr14:40Mb). Table S3 PIP 블록은 삭제 후 재작성 대상 |
| 살아남는 결론을 한 문장으로 | ✅ §8 |

---

## 8. 살아남는 결론 (한 문장)

> **TSHR의 단일 instrument 지위는 fine-mapping이 아니라 LD clumping이 뒷받침한다** — P<5×10⁻⁸,
> r²<0.001의 1차 선정 기준에서 유럽 참조는 정확히 하나의 독립 instrument(rs179252)를 산출하며
> (12개 조건 전부에서 rs179252 유지, `TaskB_02`로 검증됨), 반면 이 locus의 SuSiE fine-mapping은
> 대립유전자 정합 누락으로 무효이고 정합 후에도 z–LD 불일치와 동아시아 미수렴이 남아
> credible set 구조나 PIP를 보고할 근거가 되지 못한다.

---

## 9. 원고 수정 대상 (Task 4 통과 후 일괄 적용, 지금은 기록만)

| # | 위치 | 조치 |
|---|---|---|
| 1 | `MANUSCRIPT_..._v5_MASTER.md:316-326` PIP 블록 전체 | **삭제.** 재현 불가이며 2개 변이는 다른 locus |
| 2 | "resolved … into a **single 95% credible set**" (Results 및 Table S3) | 철회. 실측은 2개이고 그마저 정합 오류의 산물 |
| 3 | "purity = 0.993", "log₁₀ Bayes factor = 88.2" | 삭제 — 근거 없음 |
| 4 | "~116-kb window" | 삭제. 실측 리드 span은 15.3 kb |
| 5 | Methods 'Fine-mapping' 문단 (`:72`) | SuSiE 근거를 철회하고 clumping+LD 근거로 재작성. §8 문장 사용 |
| 6 | Results `:108`, Discussion `:132`, `:142` fine-mapping 문장 | 동일하게 재작성 |
| 7 | Figure S1 (`:391`) 및 `FigureS1_TSHR_finemapping_LD.{png,pdf}` | credible set 기반 패널 재검토. "~15-kb"는 실측과 맞으므로 유지 가능 |
| 8 | Table S3 clumping 표 | **변경 없음 — 정확합니다** |
| 9 | Limitations | EAS 미수렴, z–LD 불일치, blood eQTL에 1000G LD 적용의 한계를 명시 |

원고·표·그림은 이 보고서 작성 시점까지 **전혀 수정하지 않았습니다.**
기존 스크립트도 삭제·수정하지 않았고, 개정은 `internal/scripts/` 아래 신규 파일로만 했습니다.

---

## 10. Task 4에 대한 함의

- Task 4의 coloc은 `coloc.abf`이며 SuSiE를 쓰지 않으므로 이 결함의 **직접** 영향은 없습니다.
- 다만 **동일한 대립유전자 정합 문제가 coloc에도 있는지 반드시 확인해야 합니다.**
  `taskE_01_coloc.R`은 p-value 기반이라 부호를 쓰지 않아 현재는 무해하지만, Task 4가 지시하는
  **beta/varbeta 전환 시 부호가 중요해집니다.** eQTLGen effect allele, 각 outcome GWAS의
  effect allele, MAF 참조의 A1 — 세 가지를 명시적으로 정합해야 합니다.
- 원고가 "rs179252가 shared causal variant"라고 할 때 근거는 coloc의 `SNP.PP.H4`이지 SuSiE가
  아닙니다. Task 4에서 `SNP.PP.H4`를 실제로 보고하면 이 주장을 fine-mapping 없이 지지할 수
  있는지 판정됩니다.

---

## 11. 재현 방법

```bash
# 1) 원본 재실행 (스크립트 무수정) — 원본 백업 후 실행할 것
Rscript TrackA_MR/v5_upgrade/scripts/taskB_04_susie_finemap.R

# 2) 전체 PIP + 진단 (estimate_s_rss, kriging_rss)
Rscript internal/scripts/taskB_04b_susie_rerun_diagnostics.R

# 3) 대립유전자 정합 검정 (결정적 증거)
Rscript internal/scripts/taskB_04c_allele_harmonisation_test.R
```

2)·3)은 각각 10–20분 소요됩니다 (`estimate_s_rss`의 6,989×6,989 고유분해).

보존 원본: `internal/preserved_originals/` (mtime·md5 유지)

| 파일 | md5 |
|---|---|
| `TaskB_02_ld_clumping_sensitivity_v1.csv` | `51c9d6474ce6e0ae1b932084c686394c` |
| `TaskB_03_tshr_cojo_independent_signals_v1.tsv` | `870c38d81f671ff765374bd24b0f1572` |
| `TaskB_04_tshr_susie_credible_sets_v1.tsv` | `33fc8d11d1802af67b182818fbb08385` |
| `TaskB_05_tshr_finemap_summary_v1.csv` | `3f327454582b1de398793672ff2db381` |
