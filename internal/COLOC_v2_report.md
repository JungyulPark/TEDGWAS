# Task 4 + Task 5 — Colocalization v2 및 chr16p11.2 재검토

**작성일:** 2026-09-04
**스크립트:** `TrackA_MR/v5_upgrade/scripts/taskE_01b_coloc_v2.R` (신규; v1 미삭제·미수정)
**출력:** `05_coloc_candidates/TaskE_01b_coloc_v2_20260904.csv` (81행 = 9 유전자 × 3 outcome × 3 prior)
**상태:** ✅ Task 4 통과 기준 4개 충족 / ✅ Task 5 해결

---

## 요약

1. **`TSHR`의 colocalization은 재분석에도 그대로입니다.** BBJ 0.951 → **0.951**, FinnGen 0.986 →
   **0.986** (소수 3자리 불변). beta/varbeta 전환, FinnGen N 정정, 인종 문제 제거를 모두 적용한
   뒤에도 변하지 않았습니다. **논문의 핵심 주장은 견고합니다.** §5
2. **UKB colocalization을 처음 산출했습니다.** 세 가지가 새로 드러났습니다:
   - `CTLA4` × UKB = **0.953** — 양성대조가 두 번째 outcome에서도 강하게 검증됨
   - `TSHR` × UKB = **0.226** (H3 = 0.774) — TSHR은 UKB에서 **공유 인과변이가 아닙니다**
   - `IGF1R` × UKB = **0.404** (H2 = 0.400) — IGF1R의 최고값. 0.80 미만이라 "강한 coloc 없음"은
     유지되나 **"H2-dominant"라는 서술은 UKB에서 성립하지 않습니다.** §6
3. **"novel target 없음" 결론은 강화됩니다.** `TNFSF14`(BBJ 0.994)와 `IFNGR1`(BBJ 0.989)가
   이제 **두 개의** 독립 replication outcome에서 실패합니다 (UKB 0.237 / 0.052, FinnGen
   0.017 / 0.020). §6.3
4. ⚠️ **prior 민감도에서 2건이 0.80 임계를 넘나듭니다** — `TSHR` × BBJ (p12=1e-6 에서 **0.661**)와
   `CTLA4` × UKB (**0.672**). `TSHR` × FinnGen은 0.875로 견고합니다. §7
5. ⚠️ **새 결함 — outcome 파일은 GRCh38, eQTLGen/1000G는 GRCh37입니다.** v1은 hg19 좌표로 만든
   윈도우를 hg38 위치에 적용했습니다(chr14 약 466 kb, chr2 약 865 kb 어긋남). v2는 rsID 매칭으로
   전환해 제거했고, **사후확률에 미친 영향은 미미**했습니다. §8
6. ⚠️ **새 결함 — UKB outcome의 인용이 틀렸습니다.** GCST90038636은 Sakaue 2021 Nat Genet이
   아니라 **Dönertaş 2021 Nat Aging (PMID 33959723)**입니다. §3.2
7. **Task 5 해결:** chr16p11.2 trio는 **단일 신호가 맞습니다.** GCTA-COJO 조건부 분석에서
   독립 신호 1개만 선택되고 세 변이 모두 조건부 P = 0.81–0.93으로 무너집니다. 다만 원고가
   근거로 든 **유럽 r²는 그 결론을 지지하지 않으며**, 동아시아 r²로 교체해야 합니다. §9

---

## 1. v1의 문제와 v2의 대응

| # | v1 | v2 |
|---|---|---|
| 1 | p-value 기반 ABF → MAF·N·s가 사후확률에 직접 반영 | **beta/varbeta** 사용. 각 연구의 자체 표준오차 사용 |
| 2 | 단일 1000G **EUR** MAF를 두 dataset·두 outcome 전부에 적용 (동아시아 BBJ 포함) | eQTLGen(유럽)에만 EUR MAF 사용해 β/SE 유도. **outcome 쪽은 MAF 불필요** — 문제가 패치가 아니라 소멸 |
| 3 | **UKB coloc 없음** | 세 outcome 전부 |
| 4 | `SNP.PP.H4` 계산하되 미보고 | top SNP·값·**H4 조건부 95% credible set** 보고 |
| 5 | FinnGen N = 520,387 (placeholder) | **N = 500,348, s = 858/500,348** (Task 2 확정) |
| 6 | prior 고정 (p12=1e-5) | **p12 = 1e-6 / 5e-6 / 1e-5** |
| 7 | hg19 윈도우를 hg38 위치에 적용 | **rsID 매칭** (§8) |
| 8 | UKB LMM 스케일 미처리 | **log-odds 재척도** (Lloyd-Jones 2018), MR 파이프라인과 동일 |

eQTLGen β/SE는 `taskD_03_mr_3outcomes_v1.R`의 `prep_exposure()` 로직을 재사용했습니다
(Zhu 2016): `se = 1/sqrt(2·eaf·(1−eaf)·(N_snp + Z²))`, `beta = Z·se`. **SNP별 `NrSamples`**를
사용하므로 Task 3에서 지적한 "전체에 31,684 일괄 적용" 문제도 여기서는 없습니다.
eQTLGen dataset은 `type="quant", sdY=1` — Zhu 근사가 표현형 SD 단위의 β를 주기 때문입니다.

---

## 2. 실행 환경

| 항목 | 값 |
|---|---|
| 실행 | 2026-09-04 |
| R / coloc | 4.3.3 / **5.2.3** |
| 유전자 (9) | TSHR, IGF1R, CTLA4 (backbone) + TNFSF14, IFNGR1, MAPKAPK5, HSD3B7, VKORC1, PRSS36 (후보) |
| 노출 | eQTLGen cis-eQTL, `type="quant"`, `sdY=1`, N=31,684 (β/SE는 SNP별 N) |
| 매칭 | **rsID** (빌드 무관), 대립유전자 정합 후 |

---

## 3. 확정된 outcome 파라미터

| outcome | cases | controls | N | s | 인종 | 척도 |
|---|---|---|---|---|---|---|
| BBJ Graves (GCST90018627) | 2,809 | 172,656 | 175,465 | 0.01600889 | East Asian (Japan) | log-odds |
| UKB hyperthyroid (GCST90038636) | 3,731 | 480,867 | 484,598 | 0.00769917 | NR (U.K.) | **linear LMM → 재척도** |
| FinnGen R12 GRAVES_OPHT | 858 | 499,490 | 500,348 | 0.00171481 | European (Finnish) | log-odds |

### 3.1 BBJ — 스크립트가 옳고 `data_sources.csv`가 틀렸습니다

GWAS Catalog GCST90018627: *"2,809 East Asian ancestry cases, 172,656 East Asian ancestry
controls"*, *"175465 East Asian (Japan)"*, PubMed 34594039 (Sakaue).
→ `taskE_01_coloc.R`의 `s=2809/175465, N=175465`는 **정확합니다.**

`00_meta/data_sources.csv`는 "2809 cases / **209644** controls"라고 적고 있는데, 209,644는
**다른 데이터셋**(OpenGWAS `bbj-a-123`, N=212,453, ncase=2,176)의 총계에서 2,809를 뺀 값입니다
(`TrackA_MR/results/graves_gwas.csv` 참조). FinnGen과 같은 유형의 오류입니다.
올바른 값은 **172,656**.

`TED_TRAP_Task.md:327`은 총계는 맞게(175,465 / 2,809 / 172,656) 적었으나 인종을 "European"으로
잘못 표기했습니다.

**파일 내부 교차검증:** BBJ의 `effect_allele_frequency`를 1000G chr14 패널과 대조한 결과
EAS와 r = **0.966**(평균절대차 0.028), EUR과 r = 0.827(0.109). 동아시아가 맞고 디스크의 파일도
올바른 파일입니다.

### 3.2 ⚠️ UKB — 수치는 맞고 **인용이 틀렸습니다**

GWAS Catalog GCST90038636: *"3,731 cases, 480,867 controls"*, *"484598 NR (U.K.)"*
→ `s = 3731/484598 = 0.00769917`, `N = 484,598` **확인**. 저장소 내 출처는
`TrackA_MR/results/hyper_gwas.csv`(OpenGWAS 메타데이터 덤프)이며 값이 일치합니다.

그러나 같은 페이지가 이 연구를 다음으로 명시합니다:

> **PubMed ID 33959723 / First author Dönertaş HM / Nat Aging / 2021-04-08**
> *"Common genetic associations between age-related diseases."*

원고는 이 outcome을 **참고문헌 [13] Sakaue et al. 2021 Nat Genet**으로 인용합니다
(본문 `:46`, Table 1 `:216`). `data_sources.csv`의 UKB 행도 pmid `34594039`, doi
`10.1038/s41588-021-00931-x`, citation "Sakaue et al. 2021 ... UKB hyperthyroid endpoint"로
적혀 있습니다. **모두 틀렸습니다.**

[13](Sakaue)은 **BBJ outcome에는 맞습니다.** UKB outcome에만 잘못 붙어 있습니다.

PubMed로 두 논문의 서지사항을 확인했습니다:

| | PMID | DOI | 서지 |
|---|---|---|---|
| Dönertaş HM, Fabian DK, Valenzuela MF, Partridge L, Thornton JM — *Common genetic associations between age-related diseases* | 33959723 | `10.1038/s43587-021-00051-5` | Nat Aging 2021;**1**(4):400–412 |
| Sakaue S, Kanai M, Tanigawa Y, et al — *A cross-population atlas of genetic associations for 220 human phenotypes* | 34594039 | `10.1038/s41588-021-00931-x` | Nat Genet 2021;**53**(10):1415–1424 |

Dönertaş 초록은 *"116 diseases in the **UK Biobank** data"*, 키워드에 `UK Biobank`·`GWAS`가
있습니다. Sakaue 초록은 *"220 … GWAS … in **BioBank Japan** (n = 179,000)"*입니다.

**혼동의 원인으로 보이는 것:** Sakaue 초록에 *"Meta-analyses with the UK Biobank and FinnGen
(n = 628,000)"*라는 문장이 있습니다. Sakaue 연구가 UKB를 포함한 것은 사실이나, **accession
GCST90038636 자체는 Dönertaş에 귀속**됩니다. 우리가 내려받아 쓴 파일은 그 accession입니다.

부수적으로, GWAS Catalog의 인종 라벨은 "NR (U.K.)"로 **미보고**입니다. 원고의 "European"은
UK Biobank 구성상 합리적 추정이나 카탈로그가 명시한 값은 아닙니다.

### 3.3 UKB 척도 재척도

`beta`가 BOLT-LMM 선형 척도이므로 `β_logodds = β_LMM / (μ(1−μ))`, SE도 동일하게 나눔
(μ = 3731/484598 = 0.0076992, 제수 = 0.00763989). Z와 P는 보존됩니다.
MR 파이프라인(`taskD_03_mr_3outcomes_v1.R`)과 동일한 처리입니다.

---

## 4. 결과 — 전체 (p12 = 1e-5)

| gene | outcome | n | PP.H2 | PP.H3 | PP.H4 | top SNP (SNP.PP.H4) | CS95 | outcome min P |
|---|---|---|---|---|---|---|---|---|
| **TSHR** | BBJ | 4,336 | 0.000 | 0.049 | **0.951** | rs179252 (0.549) | 3 | 1.83e-15 |
| **TSHR** | **UKB** | 5,667 | 0.000 | **0.774** | **0.226** | rs1023586 (0.489) | 3 | 1.40e-43 |
| **TSHR** | FinnGen | 6,579 | 0.000 | 0.014 | **0.986** | rs179252 (0.525) | 4 | 1.23e-07 |
| IGF1R | BBJ | 5,685 | 0.690 | 0.236 | 0.073 | rs2654980 (0.607) | 3 | 5.72e-04 |
| **IGF1R** | **UKB** | 7,394 | **0.400** | 0.196 | **0.404** | rs2654980 (0.766) | 3 | 9.00e-05 |
| IGF1R | FinnGen | 7,740 | 0.623 | 0.346 | 0.032 | rs2654980 (0.619) | 3 | 4.35e-04 |
| CTLA4 | BBJ | 3,014 | 0.000 | 0.799 | 0.201 | rs231811 (0.216) | 11 | 5.49e-18 |
| **CTLA4** | **UKB** | 4,283 | 0.000 | 0.047 | **0.953** | rs3087243 (0.268) | 11 | 9.00e-27 |
| CTLA4 | FinnGen | 4,924 | 0.000 | 0.022 | **0.978** | rs1863800 (0.170) | 12 | 2.11e-08 |
| TNFSF14 | BBJ | 4,790 | 0.001 | 0.005 | **0.994** | rs2291668 (1.000) | 1 | 9.51e-07 |
| TNFSF14 | UKB | 7,197 | 0.360 | 0.403 | 0.237 | rs2291668 (1.000) | 1 | 5.20e-06 |
| TNFSF14 | FinnGen | 7,703 | 0.626 | 0.357 | 0.017 | rs2291668 (1.000) | 1 | 1.03e-03 |
| IFNGR1 | BBJ | 3,973 | 0.005 | 0.006 | **0.989** | rs11754268 (0.956) | 1 | 9.41e-06 |
| IFNGR1 | UKB | 5,145 | 0.736 | 0.211 | 0.052 | rs11754268 (0.914) | 2 | 6.40e-04 |
| IFNGR1 | FinnGen | 5,891 | 0.641 | 0.339 | 0.020 | rs11754268 (0.931) | 2 | 1.92e-03 |
| MAPKAPK5 | BBJ | 2,077 | 0.000 | 1.000 | 0.000 | rs79271898 (0.385) | 4 | 5.27e-14 |
| MAPKAPK5 | UKB | 2,763 | 0.052 | 0.944 | 0.005 | rs79271898 (0.402) | 4 | 4.80e-07 |
| MAPKAPK5 | FinnGen | 3,294 | 0.663 | 0.311 | 0.026 | rs79271898 (0.385) | 4 | 4.44e-04 |
| HSD3B7 | BBJ | 1,630 | 0.000 | 0.364 | 0.636 | rs4889606 (0.966) | 1 | 1.15e-08 |
| HSD3B7 | UKB | 2,755 | 0.444 | 0.461 | 0.095 | rs4889606 (0.949) | 2 | 7.00e-05 |
| HSD3B7 | FinnGen | 3,104 | 0.733 | 0.240 | 0.026 | rs4889606 (0.948) | 2 | 5.11e-04 |
| VKORC1 | BBJ | 1,516 | 0.000 | 0.612 | 0.388 | rs34649473 (0.592) | 10 | 1.15e-08 |
| VKORC1 | UKB | 2,644 | 0.477 | 0.493 | 0.030 | rs34649473 (0.499) | 11 | 7.00e-05 |
| VKORC1 | FinnGen | 2,961 | 0.739 | 0.231 | 0.030 | rs2884737 (0.280) | 11 | 5.11e-04 |
| PRSS36 | BBJ | 1,479 | 0.000 | 0.832 | 0.168 | rs78924645 (0.554) | 2 | 1.15e-08 |
| PRSS36 | UKB | 2,577 | 0.436 | 0.449 | 0.115 | rs78924645 (0.529) | 2 | 7.00e-05 |
| PRSS36 | FinnGen | 2,821 | 0.744 | 0.219 | 0.037 | rs78924645 (0.550) | 2 | 5.11e-04 |

`outcome min P`는 해당 locus에서 outcome GWAS의 최소 P입니다 — 검정력 부족과 진짜 null을
구분하기 위한 값입니다.

---

## 5. v1 대비 변화 — 핵심 결과는 불변

| gene × outcome | v1 PP.H4 | v2 PP.H4 | 변화 |
|---|---|---|---|
| **TSHR × BBJ** | 0.951 | **0.951** | **불변** |
| **TSHR × FinnGen** | 0.986 | **0.986** | **불변** |
| CTLA4 × FinnGen | 0.978 | 0.978 | 불변 |
| CTLA4 × BBJ | 0.206 | 0.201 | −0.005 |
| TNFSF14 × BBJ | 0.994 | 0.994 | 불변 |
| IFNGR1 × BBJ | 0.989 | 0.989 | 불변 |
| IGF1R × BBJ | 0.043 | 0.073 | +0.030 |
| IGF1R × FinnGen | 0.036 | 0.032 | −0.004 |
| HSD3B7 × BBJ | 0.616 | 0.636 | +0.020 |
| VKORC1 × BBJ | 0.360 | 0.388 | +0.028 |
| PRSS36 × BBJ | 0.157 | 0.168 | +0.011 |
| 기타 FinnGen | — | — | 모두 ±0.004 이내 |

**모든 정성적 판정이 유지됩니다.** 최대 변화는 +0.030 (IGF1R × BBJ, 0.043 → 0.073)이며 어느
조합도 0.80 임계를 넘나들지 않습니다. Task 2에서 예측한 대로(`N·s(1−s)` 불변) FinnGen 값들이
사실상 그대로인 것도 확인됩니다.

---

## 6. UKB 추가로 드러난 것

### 6.1 `CTLA4` — 양성대조가 두 번째 outcome에서 검증됨

UKB PP.H4 = **0.953** (H3=0.047), FinnGen 0.978. 두 유럽 outcome 모두에서 강한 공유 인과변이.
BBJ에서는 H3=0.799로 교차인종 LD 복잡성이 확인됩니다. 양성대조로서의 역할이 강화됩니다.

### 6.2 `TSHR` — UKB에서는 공유 인과변이가 아닙니다

UKB PP.H4 = **0.226**, **PP.H3 = 0.774**. UKB의 해당 locus 최소 P는 **1.40e-43**으로 세 outcome
중 가장 강한 신호인데도, coloc은 **서로 다른 인과변이**(H3)를 지지합니다. top SNP도 rs179252가
아니라 rs1023586입니다.

이는 검정력 부족이 아닙니다 — 신호는 압도적으로 강합니다. 유럽 hyperthyroidism에서의 TSHR
연관은 eQTLGen 발현 신호와 **다른** 변이에서 오는 것으로 보입니다. 원고가 rs179252를
"the shared causal variant"로 서술할 때 이 반례를 명시해야 합니다.

### 6.3 `IGF1R` — 최고값이지만 여전히 0.80 미만

UKB PP.H4 = **0.404**, PP.H2 = **0.400** — 거의 균형. rs2654980의 `SNP.PP.H4` = 0.766.
BBJ(0.073) / FinnGen(0.032)의 5~12배입니다.

- **유지되는 것:** "IGF1R는 강한 발현 colocalization을 보이지 않는다" (0.404 < 0.80)
- **수정이 필요한 것:** 원고의 *"IGF1R is PP.H2-dominant (outcome association without a shared
  eQTL signal)"*. BBJ(H2=0.690)·FinnGen(H2=0.623)에서는 맞지만 **UKB에서는 H2 ≈ H4로
  불확정**입니다. 인종이 eQTLGen과 가장 잘 맞는 outcome에서 그렇다는 점이 중요합니다.

### 6.4 후보 유전자 — 결론 강화

| gene | BBJ | UKB | FinnGen |
|---|---|---|---|
| TNFSF14 | **0.994** | 0.237 | 0.017 |
| IFNGR1 | **0.989** | 0.052 | 0.020 |

v1에서는 "BBJ에서만 coloc, FinnGen에서 실패"였습니다. 이제 **두 개의 독립 replication
outcome에서 실패**합니다. FinnGen만으로는 검정력 부족 반론이 가능했지만, UKB는 FinnGen보다
훨씬 크고(N=484,598) 해당 locus 신호도 있습니다. **"no robust novel target" 결론이 강화됩니다.**

---

## 7. ⚠️ prior 민감도 — 2건이 임계를 넘나듭니다

| gene × outcome | p12=1e-6 | p12=5e-6 | p12=1e-5 | 판정 |
|---|---|---|---|---|
| **TSHR × BBJ** | **0.661** | 0.907 | 0.951 | ⚠️ 1e-6 에서 0.80 미만 |
| **CTLA4 × UKB** | **0.672** | 0.911 | 0.953 | ⚠️ 1e-6 에서 0.80 미만 |
| TSHR × FinnGen | 0.875 | 0.972 | 0.986 | 견고 |
| CTLA4 × FinnGen | 0.816 | 0.957 | 0.978 | 견고 |
| TNFSF14 × BBJ | 0.941 | 0.988 | 0.994 | 견고 |
| IFNGR1 × BBJ | 0.899 | 0.978 | 0.989 | 견고 |

Wallace 2020은 기본값 p12=1e-5가 조건에 따라 liberal할 수 있다고 지적합니다.
**`TSHR`의 BBJ colocalization(0.951)은 prior 의존적입니다** — 더 보수적인 prior에서는
0.661로 임계 미만입니다. 반면 **TED 특이 outcome인 FinnGen에서는 0.875로 견고**합니다.
원고는 이 사실을 밝혀야 하며, 다행히 논지에는 유리합니다: TED 관련 outcome 쪽이 더 견고합니다.

---

## 8. ⚠️ 새 결함 — 게놈 빌드 불일치

세 outcome 파일은 **GRCh38**, eQTLGen과 1000G 패널은 **GRCh37**입니다.

| SNP | BBJ / UKB / FinnGen | eQTLGen / 1000G |
|---|---|---|
| rs179252 | 14:**80,969,641** | 14:**81,435,985** |
| rs179248 | 14:80,966,694 | 14:81,433,038 |
| rs3087243 | 2:203,874,196 | 2:204,738,919 |

v1(`taskE_01_coloc.R:104-106`)은 hg19 유전자 좌표로 만든 ±1 Mb 윈도우를 **hg38 위치에** 적용해
outcome을 필터링했습니다. chr14에서 약 466 kb, chr2에서 약 865 kb 어긋납니다. 결과적으로
윈도우가 유전자를 비대칭적으로 감싸 상류 약 534 kb만 포함됩니다.

**v2의 대응:** outcome 필터링을 위치가 아니라 **rsID 소속**으로 전환했습니다. coloc의 분석
집합은 어차피 eQTLGen cis SNP과의 교집합이므로 위치 필터가 불필요하고, 이러면 빌드 문제가
소멸합니다.

**영향:** SNP 수는 늘었지만(예: TSHR×BBJ 3,756 → 4,336, ×FinnGen 5,560 → 6,579) **PP.H4는
사실상 불변**이었습니다(TSHR 세 outcome 모두 소수 3자리 동일). 회수된 SNP 대부분이 신호
영역 밖이었기 때문입니다. **실재하는 결함이나 발표된 수치를 바꾸지는 않습니다.**

부수 발견: `taskE_02_coloc_candidates.R`의 `TNFSF14 = list(chr=19L, start=6669934L,
end=6669934L)`는 start == end로, 유전자 구간이 아니라 단일 지점입니다. ±1 Mb 윈도우 덕에
분석에는 지장이 없었으나 좌표 정의 자체가 잘못돼 있습니다.

---

## 9. Task 5 — chr16p11.2

### 9.1 원고의 r²는 정확합니다

1000G Phase 3 **EUR** (n=503) 실측:

| 쌍 | 실측 r² | 원고 | 판정 |
|---|---|---|---|
| rs4889606 (HSD3B7) – rs34649473 (VKORC1) | **0.8957** | 0.90 | ✅ |
| rs4889606 – rs78924645 (PRSS36) | **0.1897** | 0.19 | ✅ |
| rs34649473 – rs78924645 | **0.2044** | 0.20 | ✅ |

좌표 3개(31,011,183 / 31,066,538 / 31,154,358)와 BBJ 주변부 P(2.04e-07 / 6.35e-07 / 2.02e-06)도
일치합니다. **S3 PIP 블록과 달리 Table S5의 값들은 진짜입니다.**

### 9.2 그러나 추론이 그 숫자와 모순됩니다

원고 각주:

> *Pairwise linkage disequilibrium in the 1000 Genomes Phase 3 **European** reference (n = 503)
> **confirms that the cluster does not comprise three independent signals**: … while rs78924645
> (PRSS36) is in **weaker LD with each (r² = 0.19 and 0.20**, respectively).*

**r² = 0.19–0.20은 단일 신호의 근거가 아니라 그 반대입니다.** 숫자를 제시한 뒤 반대되는 결론을
내리고 있습니다.

### 9.3 결론 자체는 옳습니다 — 근거를 바꿔야 합니다

**(a) 인종 정합 LD.** 후보들은 **BBJ(동아시아)**에서 발견됐으므로 EAS LD가 적절한 참조입니다.

| 쌍 | EUR r² | **EAS r²** |
|---|---|---|
| rs4889606 – rs34649473 | 0.8957 | **0.8537** |
| rs4889606 – rs78924645 | 0.1897 | **0.7612** |
| rs34649473 – rs78924645 | 0.2044 | **0.8630** |

동아시아에서는 셋 다 r² ≥ 0.76 — **하나의 블록이 맞습니다.**

**(b) 조건부 분석 (결정적).** BBJ chr16 구간에 GCTA-COJO stepwise selection (EAS 참조,
p < 5e-6): **독립 신호 1개** (`rs8050588`, chr16:30,937,799, P = 1.15e-08). 이에 조건화하면:

| SNP | 주변부 β / P | **조건부 β / P** |
|---|---|---|
| rs4889606 (HSD3B7) | −0.2411 / 2.04e-07 | −0.0042 / **0.928** |
| rs34649473 (VKORC1) | −0.2361 / 6.35e-07 | −0.0096 / **0.839** |
| rs78924645 (PRSS36) | +0.2543 / 2.02e-06 | +0.0132 / **0.806** |

**세 변이 모두 완전히 소멸합니다.** chr16p11.2는 단일 신호입니다.

**따라서 Figure 1 Panel A의 "3 — chr16p11.2 one LD block" 집계는 변경 불필요합니다.**
각주만 EAS r² + COJO 결과로 교체하면 됩니다. coloc 결과(§4)도 세 유전자 모두 어느 outcome에서도
0.80에 못 미쳐 후보 배제 결론과 정합적입니다.

---

## 10. Task 4 통과 기준 대조

| 기준 | 결과 |
|---|---|
| 세 outcome 모두 결과 존재 | ✅ 9 유전자 × 3 outcome × 3 prior = 81행 |
| `TSHR` BBJ PP.H4가 EAS MAF에서 얼마인지 명시 (v1의 0.951과 비교) | ✅ **0.951 — 불변.** beta/varbeta 전환으로 outcome 쪽 MAF 의존이 제거되어 "EUR MAF를 EAS에 적용" 문제 자체가 소멸 |
| `IGF1R` **UKB** PP.H4/PP.H2가 처음으로 보고됨 | ✅ **PP.H4 = 0.404, PP.H2 = 0.400** (rs2654980, SNP.PP.H4 = 0.766) |
| prior 민감도 3조건에서 결론이 뒤집히는지 명시 | ✅ **2건 뒤집힘** — TSHR × BBJ (0.661 @1e-6), CTLA4 × UKB (0.672 @1e-6). §7 |

`coloc.susie`는 **돌리지 않았습니다.** Task 3에서 확인된 대립유전자 정합 결함이 부호 있는 LD를
쓰는 모든 소비자에 상속되며, `coloc.susie`가 여기 해당합니다. 정합된 LD 파이프라인을 먼저
확립하지 않고 돌리면 Task 3과 같은 산물이 나옵니다.

---

## 11. 원고 수정 대상 (Task 2·3 항목과 함께 일괄 적용)

원고·표·그림은 **아직 수정하지 않았습니다.**

| # | 대상 | 조치 | 근거 |
|---|---|---|---|
| 1 | 본문 `:46`, Table 1 `:216` UKB 인용 [13] | **Dönertaş HM et al. (2021) Nat Aging, PMID 33959723**으로 교체. [13](Sakaue)은 BBJ에는 유지 | §3.2 |
| 2 | `data_sources.csv` UKB 행 pmid/doi/citation | 동일하게 정정 | §3.2 |
| 3 | `data_sources.csv` BBJ 행 controls 209,644 | **172,656**로 정정 | §3.1 |
| 4 | Table 2 / Figure 3A PP.H4 값 | **변경 없음** — 재분석에서 불변 | §5 |
| 5 | Table 2 / Figure 3A | **UKB 열 추가** (지금까지 공란) | §6 |
| 6 | "IGF1R is PP.H2-dominant" | UKB에서는 성립하지 않음을 명시 (H2 0.400 vs H4 0.404) | §6.3 |
| 7 | "rs179252 … the shared causal variant" | UKB에서 H3=0.774로 **다른** 인과변이를 지지함을 명시 | §6.2 |
| 8 | Methods colocalization 문단 | beta/varbeta, 인종별 처리, prior 민감도, `SNP.PP.H4` 보고 반영 | §1 |
| 9 | Results "no robust novel target" | UKB 추가로 **강화**됨을 반영 (두 replication outcome에서 실패) | §6.4 |
| 10 | Limitations | prior 민감도(TSHR×BBJ 0.661 @p12=1e-6), UKB 인종 라벨 "NR" 명시 | §7, §3.2 |
| 11 | Table S5 각주 | EUR r² 문장을 **EAS r² + COJO 조건부 결과**로 교체 | §9 |
| 12 | Figure 1 Panel A "chr16p11.2 one LD block" | **변경 없음 — 결론 유지** | §9.3 |
| 13 | `taskE_02_coloc_candidates.R` TNFSF14 좌표 | start == end 수정 | §8 |

---

## 12. 재현 방법

```bash
# 입력 준비 (사전 필터링; eQTLGen 3.9 GB 전체 로딩 회피)
bash internal/scripts/prepare_coloc_v2_inputs.sh

# coloc v2
Rscript TrackA_MR/v5_upgrade/scripts/taskE_01b_coloc_v2.R

# Task 5: chr16 LD + 조건부 분석
bash internal/scripts/chr16_ld_and_cojo.sh
```

출력:
- `05_coloc_candidates/TaskE_01b_coloc_v2_20260904.csv` (81행)
- `internal/COLOC_v2_report.md` (이 문서)
