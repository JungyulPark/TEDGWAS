# TED-TRAP 상세 실행 계획서
## Mendelian Randomization + Differential Pathway Network + Off-target Analysis

**프로젝트**: TED-TRAP (Thyroid Eye Disease — TSHR vs IGF-1R Pathway Analysis)
**연구자**: 박정열 (Park Jungyul, MD, PhD)
**작성일**: 2026-04-14
**작업 폴더**: `ProjectTEDGWAS/`
**예상 기간**: 4-5주 (분석) + 2주 (논문 초안)
**논문 1편 통합**: Track A (MR) + Track B (Network) + Track C (Off-target, Fig 4로 통합)

---

## 📁 폴더 구조 (첫날 생성)

```
ProjectTEDGWAS/
├── TrackA_MR/
│   ├── scripts/          # R 스크립트
│   ├── data/             # eQTL, GWAS 데이터
│   ├── results/          # MR 결과 테이블
│   └── figures/          # Forest plot, scatter plot
├── TrackB_Network/
│   ├── scripts/          # Python 스크립트
│   ├── data/             # CTD, DisGeNET, STRING 다운로드
│   │   ├── drug_targets/ # IGF1R pathway, TSHR pathway gene lists
│   │   ├── disease_genes/# TED/Graves disease genes
│   │   └── ppi/          # STRING PPI 데이터
│   ├── results/          # Intersection, enrichment 결과
│   └── figures/          # Venn, network, differential pathway
├── TrackC_Offtarget/
│   ├── scripts/
│   ├── data/             # Hearing loss genes, insulin pathway, GTEx
│   ├── results/
│   └── figures/          # Off-target network (→ 논문 Fig 4)
├── Manuscript/
│   ├── drafts/
│   ├── figures_final/
│   └── supplementary/
└── Literature/
    └── key_references/
```

---

## 파일 명명 규칙

`TED_{Track}_{Step}_{Description}_{YYYYMMDD}.{ext}`

예시:
- `TED_A_Step1_eQTL_extraction_20260414.R`
- `TED_B_Step1_IGF1R_pathway_genes_20260414.csv`
- `TED_C_Step2_hearing_loss_genes_20260416.csv`

---

## ════════════════════════════════════════
## PHASE 0: 사전 준비 (Day 0)
## ════════════════════════════════════════

### Task 0-1: Outcome GWAS 확인 ⚠️ 최우선

MR 전체가 이것에 달려 있으므로, 다른 모든 것보다 먼저 확인해야 함.

**실행 (R)**:
```r
library(ieugwasr)

# 1. OpenGWAS에서 Graves disease 검색
graves_search <- gwassearch("Graves")
print(graves_search[, c("id", "trait", "sample_size", "year")])

# 2. 후보 ID 정보 확인
gwas_info <- gwasinfo(id = c("ieu-a-1098"))
print(gwas_info[, c("id", "trait", "sample_size", "ncase", "ncontrol", "population")])

# 3. Hyperthyroidism도 확인
hyper_search <- gwassearch("hyperthyroidism")
print(hyper_search[, c("id", "trait", "sample_size")])
```

**FinnGen 확인**:
- https://r10.finngen.fi/ 접속
- "E4_GRAVES" 또는 "E4_THYROTOXICOSIS" 검색
- Case/Control 수, summary stats 다운로드 가능 여부 확인

**결정 기준**: N이 가장 크고 European population인 GWAS 선택
**기록**: 선택한 GWAS ID, N, population을 `TrackA_MR/data/outcome_gwas_info.txt`에 기록

**⚠️ 위험 요소**: 
- ieu-a-1098이 접근 불가거나 N이 너무 작을 수 있음
- FinnGen summary stats 다운로드에 시간 소요 (수 GB)
- Graves disease GWAS 자체가 few thousand 수준이면 MR power 부족 우려

---

### Task 0-2: eQTLGen 데이터 확인

**이미 보유한 파일 확인**:
```
c:\Projectbuild\CP3\data\2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz
```

이 파일을 `ProjectTEDGWAS/TrackA_MR/data/` 에 복사 또는 심볼릭 링크

**확인 사항**:
- 파일 크기 및 무결성 확인
- 헤더 구조 확인 (근시 때와 동일한지)
- TSHR (ENSG00000165409)가 이 파일에 존재하는지 빠르게 grep

```bash
zcat eQTLGen_file.gz | head -1
zcat eQTLGen_file.gz | grep "ENSG00000165409" | head -5
```

---

### Task 0-3: R/Python 환경 확인

**R 패키지**:
```r
# 필요 패키지 확인
library(TwoSampleMR)  # MR
library(ieugwasr)     # GWAS query
library(coloc)        # Colocalization
library(ggplot2)      # Plotting
library(dplyr)        # Data manipulation
library(forestplot)   # Forest plot (optional)
```

**Python 패키지**:
```python
# 필요 패키지
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib_venn import venn2, venn3  # Venn diagram
import networkx as nx                     # Network analysis
from scipy import stats
import requests                           # STRING API
```

---

## ════════════════════════════════════════
## PHASE 1: Track A — Mendelian Randomization (Week 1-2)
## ════════════════════════════════════════

### A-Step 1: Exposure eQTL 추출 (Day 1)

**목적**: 12개 TED 관련 유전자의 cis-eQTL을 eQTLGen에서 추출

**유전자 목록 (3개 모듈)**:

| # | Module | Gene | Ensembl ID | 역할 | 기대 |
|---|--------|------|------------|------|------|
| 1 | IGF-1R pathway | IGF1R | ENSG00000140443 | Teprotumumab target | 핵심 비교 대상 |
| 2 | IGF-1R pathway | IGF1 | ENSG00000017427 | IGF-1R ligand | 리간드 효과 |
| 3 | IGF-1R pathway | IRS1 | ENSG00000169047 | IGF-1R adaptor | Downstream |
| 4 | IGF-1R pathway | AKT1 | ENSG00000142208 | PI3K/AKT downstream | Downstream |
| 5 | TSHR pathway | TSHR | ENSG00000165409 | Primary autoantigen | 핵심 비교 대상 |
| 6 | TSHR pathway | ARRB1 | ENSG00000137486 | TSHR-IGF-1R scaffold | Bridge 유전자 |
| 7 | Inflammatory | IL6 | ENSG00000136244 | Tocilizumab target | 대안 치료 |
| 8 | Inflammatory | TNF | ENSG00000232810 | Pro-inflammatory | 보조 |
| 9 | HA synthesis | HAS2 | ENSG00000153446 | Hyaluronan 합성 | Effector |
| 10 | Adipogenesis | PPARG | ENSG00000132170 | Adipocyte 분화 | Effector |
| 11 | Fibrosis | TGFB1 | ENSG00000105329 | ECM remodeling | M-LIGHT 교훈 |
| 12 | Immune | CTLA4 | ENSG00000163599 | Immune checkpoint | GWAS hit |

**스크립트**: `TED_A_Step1_eQTL_extraction_20260414.R`

```r
# === eQTLGen에서 TED target gene eQTL 추출 ===

library(data.table)

# eQTLGen 파일 경로
eqtl_file <- "TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz"

# Target genes
target_genes <- data.frame(
  gene = c("IGF1R", "IGF1", "IRS1", "AKT1", "TSHR", "ARRB1",
           "IL6", "TNF", "HAS2", "PPARG", "TGFB1", "CTLA4"),
  ensembl = c("ENSG00000140443", "ENSG00000017427", "ENSG00000169047",
              "ENSG00000142208", "ENSG00000165409", "ENSG00000137486",
              "ENSG00000136244", "ENSG00000232810", "ENSG00000153446",
              "ENSG00000132170", "ENSG00000105329", "ENSG00000163599"),
  module = c(rep("IGF1R_pathway", 4), rep("TSHR_pathway", 2),
             rep("Inflammatory", 2), "HA_synthesis", "Adipogenesis",
             "Fibrosis", "Immune")
)

# 읽기 (대용량이므로 fread 사용)
cat("Reading eQTLGen data...\n")
eqtl <- fread(eqtl_file)
cat(sprintf("Total rows: %d\n", nrow(eqtl)))

# 각 유전자별 추출
for (i in 1:nrow(target_genes)) {
  gene_name <- target_genes$gene[i]
  ensembl_id <- target_genes$ensembl[i]
  
  gene_eqtl <- eqtl[Gene == ensembl_id | GeneSymbol == gene_name]
  
  cat(sprintf("\n%s (%s): %d cis-eQTLs total\n", gene_name, ensembl_id, nrow(gene_eqtl)))
  
  # Bonferroni significant만
  sig <- gene_eqtl[BonferroniP < 0.05]
  cat(sprintf("  Bonferroni significant: %d\n", nrow(sig)))
  
  # 저장
  outfile <- sprintf("TrackA_MR/data/eqtl_%s_all.csv", gene_name)
  fwrite(gene_eqtl, outfile)
  
  outfile_sig <- sprintf("TrackA_MR/data/eqtl_%s_sig.csv", gene_name)
  fwrite(sig, outfile_sig)
}

cat("\n=== Summary ===\n")
cat("Check: Do TSHR and IGF1R have sufficient IVs (>=3)?\n")
cat("If any gene has 0 significant eQTLs, exclude from MR\n")
```

**체크포인트**:
- [ ] TSHR에 충분한 IV가 있는가? (blood eQTL에서 TSHR 발현이 낮을 수 있음)
- [ ] IGF1R에 충분한 IV가 있는가?
- [ ] 각 유전자별 IV 수 기록 → `TrackA_MR/results/iv_summary.csv`

**⚠️ 핵심 위험**: TSHR은 주로 thyrocyte에서 발현되므로, blood eQTL에서 IV가 부족할 수 있음. 이 경우:
- p-value threshold를 5e-6으로 완화 (근시 때와 동일)
- GTEx thyroid tissue eQTL 대안 검토
- Limitation에 반드시 기술

---

### A-Step 2: MR 실행 — Primary Analysis (Day 2-3)

**스크립트**: `TED_A_Step2_MR_primary_20260415.R`

**방법 A: TwoSampleMR API 직접 사용** (OpenGWAS 접근 가능 시):

```r
library(TwoSampleMR)

# === 설정 ===
genes <- list(
  IGF1R = "eqtl-a-ENSG00000140443",
  TSHR  = "eqtl-a-ENSG00000165409",
  IL6   = "eqtl-a-ENSG00000136244",
  TGFB1 = "eqtl-a-ENSG00000105329",
  HAS2  = "eqtl-a-ENSG00000153446",
  PPARG = "eqtl-a-ENSG00000132170",
  ARRB1 = "eqtl-a-ENSG00000137486",
  IGF1  = "eqtl-a-ENSG00000017427",
  IRS1  = "eqtl-a-ENSG00000169047",
  AKT1  = "eqtl-a-ENSG00000142208",
  TNF   = "eqtl-a-ENSG00000232810",
  CTLA4 = "eqtl-a-ENSG00000163599"
)

outcome_id <- "ieu-a-1098"  # ← Task 0-1에서 확인한 ID로 수정

# === MR 루프 ===
all_results <- list()
all_steiger <- list()

for (gene_name in names(genes)) {
  cat(sprintf("\n========== %s ==========\n", gene_name))
  
  # 1) IV 추출
  exposure <- tryCatch(
    extract_instruments(
      outcomes = genes[[gene_name]],
      p1 = 5e-6,     # eQTL threshold (근시와 동일)
      clump = TRUE,
      r2 = 0.001,    # Clumping r2
      kb = 10000      # Clumping window
    ),
    error = function(e) { cat("  ERROR:", e$message, "\n"); NULL }
  )
  
  if (is.null(exposure) || nrow(exposure) == 0) {
    cat("  ❌ Insufficient IVs — SKIP\n")
    next
  }
  cat(sprintf("  IVs: %d\n", nrow(exposure)))
  
  # 2) Outcome extraction
  outcome <- extract_outcome_data(
    snps = exposure$SNP,
    outcomes = outcome_id
  )
  
  if (is.null(outcome) || nrow(outcome) == 0) {
    cat("  ❌ No outcome data — SKIP\n")
    next
  }
  
  # 3) Harmonize
  harmonised <- harmonise_data(exposure, outcome)
  cat(sprintf("  Harmonised SNPs: %d\n", nrow(harmonised[harmonised$mr_keep,])))
  
  # 4) MR
  mr_result <- mr(harmonised, method_list = c(
    "mr_ivw",                    # Primary
    "mr_egger_regression",       # Pleiotropy check
    "mr_weighted_median",        # Robustness
    "mr_weighted_mode"           # Robustness
  ))
  
  # 5) Steiger directionality
  steiger <- tryCatch(
    directionality_test(harmonised),
    error = function(e) { cat("  Steiger error:", e$message, "\n"); NULL }
  )
  
  # 6) Heterogeneity
  het <- mr_heterogeneity(harmonised)
  
  # 7) Pleiotropy
  pleio <- mr_pleiotropy_test(harmonised)
  
  # 8) Store results
  all_results[[gene_name]] <- list(
    mr = mr_result,
    steiger = steiger,
    het = het,
    pleio = pleio,
    n_iv = nrow(exposure),
    harmonised = harmonised
  )
  
  # Print key result
  ivw <- mr_result[mr_result$method == "Inverse variance weighted", ]
  if (nrow(ivw) > 0) {
    cat(sprintf("  IVW: β=%.3f, SE=%.3f, P=%.2e\n", ivw$b, ivw$se, ivw$pval))
  }
  
  Sys.sleep(2)  # API rate limit
}

# === 결과 정리 ===
summary_table <- do.call(rbind, lapply(names(all_results), function(g) {
  mr <- all_results[[g]]$mr
  ivw <- mr[mr$method == "Inverse variance weighted", ]
  data.frame(
    Gene = g,
    N_IV = all_results[[g]]$n_iv,
    Method = "IVW",
    Beta = round(ivw$b, 4),
    SE = round(ivw$se, 4),
    P = formatC(ivw$pval, format = "e", digits = 2),
    Significant = ifelse(ivw$pval < 0.05/12, "Yes (Bonferroni)", 
                         ifelse(ivw$pval < 0.05, "Yes (nominal)", "No"))
  )
}))

write.csv(summary_table, "TrackA_MR/results/MR_primary_summary.csv", row.names = FALSE)
print(summary_table)

# === Plots ===
for (gene_name in names(all_results)) {
  h <- all_results[[gene_name]]$harmonised
  
  # Scatter plot
  p <- mr_scatter_plot(all_results[[gene_name]]$mr, h)
  ggsave(sprintf("TrackA_MR/figures/scatter_%s.pdf", gene_name), p[[1]], width=8, height=6)
  
  # Forest plot (if >1 IV)
  if (all_results[[gene_name]]$n_iv > 1) {
    res_single <- mr_singlesnp(h)
    p_forest <- mr_forest_plot(res_single)
    ggsave(sprintf("TrackA_MR/figures/forest_%s.pdf", gene_name), p_forest[[1]], width=8, height=8)
  }
  
  # Funnel plot
  if (all_results[[gene_name]]$n_iv > 2) {
    res_single <- mr_singlesnp(h)
    p_funnel <- mr_funnel_plot(res_single)
    ggsave(sprintf("TrackA_MR/figures/funnel_%s.pdf", gene_name), p_funnel[[1]], width=8, height=6)
  }
}

# === 전체 결과 저장 ===
saveRDS(all_results, "TrackA_MR/results/MR_all_results.rds")
```

**체크포인트**:
- [ ] TSHR IVW p-value 기록
- [ ] IGF1R IVW p-value 기록
- [ ] 둘의 비교 해석 (4개 시나리오 중 어디에 해당하는지)
- [ ] Steiger directionality 모두 correct direction인지 확인
- [ ] MR-Egger intercept p > 0.05인지 (pleiotropy 없는지)
- [ ] Heterogeneity Q p-value 기록

**핵심 해석 시나리오**:

| TSHR MR | IGF1R MR | 해석 | 논문 메시지 |
|---------|----------|------|------------|
| Causal ✅ | Causal ✅ | 두 경로 모두 인과적 | 병합 치료 근거 |
| **Causal ✅** | **Null ❌** | **TSHR가 primary** | **TSHR-targeting 강력 지지 (최적 시나리오)** |
| Null ❌ | Causal ✅ | IGF-1R이 primary | Teprotumumab 정당화 |
| Null ❌ | Null ❌ | 미확인 | Limitation (power 부족 가능성) |

---

### A-Step 3: Colocalization (Day 3-4)

**목적**: MR에서 유의한 유전자에 대해, eQTL과 GWAS signal이 같은 causal variant를 공유하는지 확인

**스크립트**: `TED_A_Step3_coloc_20260416.R`

근시 `08_coloc_full_sumstats.R`을 수정하여 사용. 핵심 변경:
- gene_info: TSHR (chr14:81,041,426-81,226,811), IGF1R (chr15:99,192,200-99,507,759) 좌표
- outcome GWAS: Graves disease ID
- Window: gene ± 500kb

**⚠️ 주의**: 
- Coloc은 full summary statistics가 필요 → OpenGWAS에서 Graves GWAS full sumstats 접근 확인
- FinnGen은 full sumstats 제공하므로 대안
- eQTLGen full sumstats는 이미 보유

**결과 해석**:
- PP.H4 > 0.7: Strong colocalization (shared causal variant)
- PP.H3 > 0.7: Different causal variants (MR 결과 주의)
- PP.H4 > 0.5: Suggestive colocalization

---

### A-Step 4: FinnGen Replication (Day 4-5)

**목적**: Primary MR 결과를 독립 GWAS로 재현

**스크립트**: `TED_A_Step4_FinnGen_replication_20260417.R`

1. https://r10.finngen.fi/ 에서 Graves disease summary stats 다운로드
2. 로컬 MR 실행 (read_outcome_data 사용)
3. Primary 결과와 방향성, 유의성 비교

**체크포인트**:
- [ ] Primary에서 유의한 유전자가 FinnGen에서도 같은 방향인지
- [ ] Replication p < 0.05면 robust, 방향만 일치해도 supportive

---

### A-Step 5: MR Summary Figure (Day 5)

**Fig 2 설계**: 12개 유전자 MR 결과를 하나의 forest plot으로

좌측에 Gene name + Module, 우측에 OR [95% CI], P-value
색상: IGF-1R pathway (파란색), TSHR pathway (빨간색), Inflammatory (주황색), Others (회색)

**스크립트**: `TED_A_Step5_MR_figure_20260418.R`

---

## ════════════════════════════════════════
## PHASE 2: Track B — Differential Pathway Network (Week 2-3)
## ════════════════════════════════════════

### B-Step 1: Drug Target Gene Set 구축 (Day 1-2)

#### B-Step 1A: IGF-1R Pathway Genes (Teprotumumab coverage)

**데이터 소스 4개 (교차 확인)**:

| # | Source | 방법 | 기대 산출 |
|---|--------|------|----------|
| 1 | CTD | ctdbase.org → Genes → "IGF1R" → Interactions (Direct Evidence only) | IGF1R interacting genes |
| 2 | DrugBank | go.drugbank.com → "teprotumumab" (DB14910) → Targets + Pathways | Formal drug target |
| 3 | STRING | string-db.org → "IGF1R" → 1st shell (confidence ≥0.7) | PPI neighbors |
| 4 | Manual curation | 문헌 기반 (아래 목록) | Domain knowledge |

**수동 추가 목록 (문헌 확인 필수)**:
```
# IGF-1R direct signaling
IGF1R, IGF1, IGF2, IRS1, IRS2, SHC1, GRB2, SOS1

# PI3K/AKT arm
PIK3CA, PIK3CB, PIK3R1, AKT1, AKT2, MTOR, RPS6KB1

# RAS/MAPK arm
KRAS, BRAF, MAP2K1, MAP2K2, MAPK1, MAPK3

# JAK/STAT arm
JAK1, JAK2, STAT3, STAT5A

# TSHR crosstalk
ARRB1, ARRB2

# TED effectors (shared)
HAS1, HAS2, HAS3, PPARG, CEBPA
IL6, CXCL8, TNF, IL1B
```

**출력**: `TrackB_Network/data/drug_targets/IGF1R_pathway_genes.csv`
- 컬럼: Gene, Source (CTD/DrugBank/STRING/Manual), Evidence_Type, Reference

#### B-Step 1B: TSHR Pathway Genes (TSHR-targeting coverage)

**동일한 4개 소스**:

| # | Source | 방법 |
|---|--------|------|
| 1 | CTD | "TSHR" → Gene Interactions (Direct Evidence only) |
| 2 | DrugBank | TSHR-targeting drugs 검색 (K1-70은 아직 없을 수 있음) |
| 3 | STRING | "TSHR" → 1st shell (confidence ≥0.7) |
| 4 | Manual curation | 문헌 기반 |

**수동 추가 목록**:
```
# TSHR direct signaling
TSHR, GNAS, GNAQ, GNA11

# cAMP/PKA arm (IGF-1R independent!)
ADCY3, ADCY5, ADCY6, ADCY7, ADCY9
PRKACA, PRKACB
CREB1, CREB3, ATF1

# β-arrestin scaffold
ARRB1, ARRB2, GRK2 (ADRBK1), GRK5

# Gq/PLC arm
PLCB1, PLCB3
PRKCB, PRKCA

# TSHR-specific downstream
FOXO1, FOXO3

# TED effectors (shared)
HAS1, HAS2, HAS3, IL6, TNF
PPARG, CEBPA, FABP4, ADIPOQ
```

**출력**: `TrackB_Network/data/drug_targets/TSHR_pathway_genes.csv`

**⚠️ M-LIGHT 교훈 적용**:
- CTD에서 Direct Evidence만 사용 (Inference score 기반 제외)
- 각 유전자의 source를 명확히 기록 (reviewer 대비)
- STRING confidence threshold 0.7 일관 적용

---

### B-Step 2: TED Disease Gene Set 구축 (Day 2-3)

**데이터 소스 5개**:

| # | Source | 검색어 | 필터 |
|---|--------|--------|------|
| 1 | CTD | "Graves Ophthalmopathy" (D016644) → Genes | Direct Evidence only |
| 2 | CTD (확장) | "Graves Disease" (D006111) → Genes | Direct Evidence only |
| 3 | DisGeNET | "thyroid eye disease" OR "Graves ophthalmopathy" | Score > 0.1 |
| 4 | GWAS Catalog | Graves disease GWAS loci | Genome-wide significant |
| 5 | Manual curation | TED 핵심 유전자 (아래 목록) | 문헌 reference 필수 |

**수동 추가 (TED-specific 유전자, reference 포함)**:
```
# GWAS loci
TSHR, HLA-DRB1, CTLA4, CD40, PTPN22, ARID5B, TG, FOXP3

# Orbital fibroblast markers
CD34, THY1, SLIT2, ROBO1

# HA synthesis/degradation
HAS1, HAS2, HAS3, UGDH, HYAL1, HYAL2

# Adipogenesis
PPARG, CEBPA, CEBPB, FABP4, ADIPOQ, DLK1

# Fibrosis/ECM
TGFB1, CTGF (CCN2), COL1A1, COL3A1, FN1, MMP2, MMP9, LOX, TIMP1, ACTA2

# Inflammation
IL6, TNF, IL1B, IFNG, IL17A, IL4, IL13, CXCL8, CCL2, CCL5

# Immune cells
CD4, CD8A, FOXP3, CTLA4, CD40, CD40LG, ICOS, PDCD1

# Oxidative stress
SOD2, NOS2, HIF1A, NFE2L2

# Wnt pathway (TED orbital fibroblast)
CTNNB1, SFRP1, WNT5A

# IGF-1R/TSHR signaling (TED context)
IGF1R, IGF1, TSHR, ARRB1
```

**출력**: `TrackB_Network/data/disease_genes/TED_disease_genes.csv`
- 컬럼: Gene, Source, Evidence_Type (Direct/GWAS/Literature), Score (if available), Reference

**⚠️ 주의**:
- CTD "Graves Ophthalmopathy" 검색 시 유전자가 매우 적을 수 있음 → "Graves Disease"로 확장 필요
- DisGeNET에서 "thyroid eye disease"가 없을 수 있음 → "Graves ophthalmopathy", UMLS CUI 확인
- Manual curation 유전자는 intersection에 넣으면 순환논증 → Extension Layer로 분리 (M-LIGHT 교훈 #2)

---

### B-Step 3: Intersection + PPI Network (Day 3-4)

**스크립트**: `TED_B_Step3_intersection_PPI_20260417.py`

```
Step 3-1: Venn diagram
  - Set A = IGF1R pathway genes
  - Set B = TSHR pathway genes
  - Set D = TED disease genes
  
  - Intersection 1: A ∩ D = Teprotumumab-TED intersection
  - Intersection 2: B ∩ D = TSHR-targeting-TED intersection
  - Three zones: (A ∩ D) only, (B ∩ D) only, (A ∩ B ∩ D) shared

Step 3-2: STRING PPI (각 intersection에 대해)
  - STRING API 또는 웹에서 confidence ≥ 0.7
  - Network metrics: degree, betweenness centrality
  - Hub gene 식별

Step 3-3: 3-zone gene list 저장
  - IGF1R_only_genes.csv (Teprotumumab이 커버하지만 TSHR-targeting이 못 커버)
  - TSHR_only_genes.csv (TSHR-targeting이 커버하지만 Teprotumumab이 못 커버)
  - Shared_genes.csv (둘 다 커버)
```

**핵심 output**:
```
IGF-1R only zone → 이게 "Teprotumumab의 고유 coverage" = 부작용 원인 후보
TSHR only zone → 이게 "30% relapse의 원인" = IGF-1R 독립 경로
Shared zone → 두 약물 공통 coverage
```

---

### B-Step 4: KEGG/GO Enrichment (Day 4-5)

**각 zone별 enrichment 수행**:

| Zone | KEGG 기대 | GO BP 기대 |
|------|-----------|-----------|
| IGF-1R only | PI3K-Akt signaling, MAPK signaling | Cell proliferation, anti-apoptosis |
| TSHR only | cAMP signaling, Thyroid hormone signaling | GPCR signaling, adenylyl cyclase activity |
| Shared | Cytokine-cytokine interaction, ECM | Inflammation, HA synthesis, adipogenesis |

**도구**: 
- Python: gseapy (Enrichr wrapper) 또는 g:Profiler API
- 또는 웹: enrichr.maayanlab.cloud, g:profiler

**출력**: 각 zone의 top 20 KEGG + top 20 GO BP → `TrackB_Network/results/enrichment/`

---

### B-Step 5: Extension Layer (Day 5)

**M-LIGHT 교훈 적용**: 가설 유전자를 intersection에 넣으면 순환논증

**Extension Layer 구성**:
1. **β-arrestin scaffold module**: ARRB1, ARRB2, GRK2, GRK5
   - 두 drug network에서 이 모듈로의 edge 수 비교
   - TSHR network가 더 직접 연결? (TSHR → ARRB1은 direct)
   
2. **CD34+ fibrocyte differentiation module**: CD34, THY1, CXCR4, CXCL12
   - Teprotumumab은 CD34+ OF의 IGF-1R 감소시킴
   - TSHR-targeting은 CD34+ OF의 TSHR 활성화 자체를 차단

**Connectivity 비교**:
- IGF-1R network → Extension Layer: edge 수, shortest path
- TSHR network → Extension Layer: edge 수, shortest path

---

### B-Step 6: CMap Drug Repurposing (Day 5-6)

**Enrichr LINCS L1000 사용** (L1000CDS2와 혼동 금지 — M-LIGHT 교훈 #4):

1. Teprotumumab intersection gene list → Enrichr → LINCS L1000 Chem Pert Up/Down
2. TSHR intersection gene list → Enrichr → LINCS L1000 Chem Pert Up/Down
3. 비교: 같은 약물이 나오면 shared pathway, 다른 약물이면 differential

**출력**: `TrackB_Network/results/cmap/`

---

### B-Step 7: Differential Pathway Figure (Day 6-7)

**Fig 3 설계** (논문의 핵심 Figure):

```
┌─────────────────┐    ┌──────────────┐    ┌─────────────────┐
│  TEPROTUMUMAB    │    │   SHARED     │    │  TSHR-TARGETING  │
│  (IGF-1R block)  │    │              │    │  (TRAb removal)  │
├─────────────────┤    ├──────────────┤    ├─────────────────┤
│ IGF-1R only genes│←→ │ Common genes │←→ │ TSHR only genes  │
│                  │    │              │    │                  │
│ PI3K/AKT/mTOR   │    │ HA synthesis │    │ cAMP/PKA/CREB   │
│ RAS/MAPK        │    │ Adipogenesis │    │ Gq/PLC/PKC      │
│ JAK/STAT        │    │ Inflammation │    │ GPCR signaling   │
│                  │    │ ECM remodel  │    │                  │
├─────────────────┤    └──────────────┘    ├─────────────────┤
│ OFF-TARGET:      │                        │ OFF-TARGET:      │
│ ⚠ Cochlea       │                        │ ✅ None predicted│
│ ⚠ Pancreas      │                        │                  │
│ ⚠ CNS           │                        │                  │
└─────────────────┘                        └─────────────────┘
```

Cytoscape 또는 Python networkx로 제작. 색상 코딩:
- IGF-1R only: 파란색
- TSHR only: 빨간색  
- Shared: 보라색
- Off-target: 회색/경고색

---

## ════════════════════════════════════════
## PHASE 3: Track C — Off-target Pathway (Week 3, Fig 4용)
## ════════════════════════════════════════

### C-Step 1: IGF-1R Tissue Expression (Day 1)

**GTEx portal** (gtexportal.org):
1. IGF1R gene page → Multi-tissue expression plot 스크린샷/데이터 다운로드
2. 관심 조직: Nerve-Tibial (cochlea proxy), Pancreas, Brain 여러 영역, Thyroid, Adipose
3. Supplementary figure 용

**⚠️**: GTEx에 cochlea/inner ear 조직이 없음 → "Nerve" 계열을 proxy로 사용하되, limitation 명시

---

### C-Step 2: Hearing Loss Gene Set (Day 2)

**소스**:
1. DisGeNET: "sensorineural hearing loss" → gene list
2. OMIM: "deafness" → gene list
3. KEGG: 직접적인 "hearing" pathway는 없으므로, PI3K-Akt (hsa04151) + 관련 경로
4. 문헌: IGF-1 KO mouse cochlear phenotype genes (Gao & Nakagawa 2020, Okano et al.)
   - NTRN (Netrin1), GAP43, SOX2, ATOH1 등

**교차 분석**:
```
IGF-1R downstream genes ∩ Hearing loss genes = Cochlear off-target genes
```

**출력**: `TrackC_Offtarget/data/hearing_loss_genes.csv`, `cochlear_offtarget_intersection.csv`

---

### C-Step 3: Hyperglycemia/Insulin Pathway (Day 3)

**소스**:
1. KEGG: Insulin signaling pathway (hsa04910) → gene list 다운로드
2. KEGG: Type II diabetes mellitus (hsa04930)
3. IGF-1R/IR hybrid receptor 관련: INSR, IGF1R, IRS1, IRS2, IGFBP1-6

**교차 분석**:
```
IGF-1R downstream genes ∩ Insulin signaling genes = Pancreatic off-target genes
```

---

### C-Step 4: CNS Pathway (Day 3)

**소스**:
1. 문헌 기반: IGF-1의 CNS 기능 관련 유전자
   - Neuronal survival: BCL2, BAX, CASP3, CASP9
   - Synaptic plasticity: BDNF, ARC, CREB1
   - Myelination: MBP, PLP1, OLIG2
2. DisGeNET: "cognitive decline" 유전자

**교차 분석**:
```
IGF-1R downstream genes ∩ CNS function genes = CNS off-target genes
```

---

### C-Step 5: 통합 Off-target Network Figure (Day 4-5)

**Fig 4 설계**:

```
                    IGF-1R
                   /  |  \  \
                  /   |   \   \
            Orbit  Cochlea Pancreas CNS
            (✅)   (⚠️)    (⚠️)   (⚠️)
            
                    TSHR
                     |
                  Orbit only
                   (✅)
```

각 조직별 off-target gene list를 하위 네트워크로 표시:
- Cochlear off-target: IGF1R→PI3K→AKT→Netrin1/Gap43 (hair cell survival)
- Pancreatic off-target: IGF1R/INSR hybrid→IRS1/2→glucose homeostasis
- CNS off-target: IGF1R→AKT→CREB→neuronal survival

**핵심 메시지**: "Teprotumumab의 IGF-1R 차단은 4개 조직에 영향을 미치지만, TSHR 직접 차단은 orbital tissue에만 작용"

---

## ════════════════════════════════════════
## PHASE 4: 논문 작성 (Week 4-5)
## ════════════════════════════════════════

### 논문 구조

**Title**: "Mendelian randomization and differential pathway network analysis reveal distinct therapeutic mechanisms of IGF-1R versus TSHR inhibition in thyroid eye disease"

**Target**: Thyroid (IF 5.4) → JCEM (IF 5.8) → Front Endocrinol (IF 3.9)

#### Abstract (250 words)
- Background: Teprotumumab 성공 but 30% relapse + 부작용, TSHR vs IGF-1R 논쟁 미해결
- Methods: MR (12 genes, eQTLGen + Graves GWAS) + Differential pathway + Off-target
- Results: [결과에 따라 작성]
- Conclusion: [결과에 따라 작성]

#### Introduction (~800 words)
- P1: TED 개요, teprotumumab 승인
- P2: Teprotumumab의 한계 (30% relapse, hearing loss, hyperglycemia)
- P3: TSHR-IGF-1R 논쟁 (두 경로의 crosstalk, β-arrestin scaffold)
- P4: 질문 제기 + 연구 목적

#### Methods (~1500 words)
- 2.1 MR: Exposure (eQTLGen), Outcome (Graves GWAS), Methods (IVW, Egger, weighted median), Sensitivity (Steiger, coloc), Replication (FinnGen)
- 2.2 Drug pathway construction: IGF-1R pathway, TSHR pathway (sources, thresholds)
- 2.3 Disease gene set: CTD, DisGeNET, GWAS, manual curation
- 2.4 Differential pathway analysis: Intersection, PPI, enrichment
- 2.5 Off-target pathway: Hearing loss, insulin, CNS gene sets, intersection
- 2.6 Software and statistics

#### Results
- 3.1 MR identifies causal pathways (Fig 2, Table 1)
- 3.2 Differential pathway reveals distinct drug coverage (Fig 3, Table 2)
- 3.3 Off-target analysis explains teprotumumab adverse effects (Fig 4)
- 3.4 Extension layer: β-arrestin scaffold connectivity

#### Discussion (~1500 words)
- D1: 주요 발견 요약
- D2: TSHR 직접 차단의 이론적 장점 (or 결과에 따라 수정)
- D3: 30% relapse 설명 (TSHR only zone = IGF-1R 비의존 경로)
- D4: 부작용 기전 (off-target pathway)
- D5: TSHR-targeting therapy 전망 (K-1-70, cyclic peptide, decoy receptor)
- D6: Limitations (blood eQTL, TED-specific GWAS 부재, computational only, manual curation bias)
- D7: Conclusion

#### Figures
- **Fig 1**: Study design schematic (Graphical abstract 겸용)
- **Fig 2**: MR forest plot (12 genes, color-coded by module)
- **Fig 3**: Differential pathway network (3-zone: IGF-1R only / Shared / TSHR only)
- **Fig 4**: Off-target pathway network (IGF-1R → 4 tissues vs TSHR → 1 tissue)

#### Supplementary
- Table S1: Complete gene lists (IGF-1R pathway, TSHR pathway, TED disease genes)
- Table S2: MR full results (all methods, all genes)
- Table S3: Coloc results
- Table S4: KEGG/GO enrichment full results
- Table S5: CMap results
- Fig S1: Venn diagrams
- Fig S2: Individual MR scatter/forest/funnel plots
- Fig S3: GTEx tissue expression

---

## ════════════════════════════════════════
## 일별 실행 체크리스트
## ════════════════════════════════════════

### Week 1: 데이터 준비 + MR

| Day | Task | Output | 체크 |
|-----|------|--------|------|
| D1 | 폴더 생성 + GWAS 확인 (Task 0-1, 0-2) | outcome_gwas_info.txt | ☐ |
| D1 | eQTLGen에서 12개 유전자 eQTL 추출 (A-Step 1) | eqtl_*.csv × 12 | ☐ |
| D2 | MR primary 실행 (A-Step 2) | MR_primary_summary.csv | ☐ |
| D3 | MR 결과 확인 + plots | scatter/forest/funnel PDFs | ☐ |
| D3-4 | Coloc 실행 (A-Step 3) | coloc_results.csv | ☐ |
| D4-5 | FinnGen replication (A-Step 4) | FinnGen_MR_results.csv | ☐ |
| D5 | MR summary figure 제작 (A-Step 5) | Fig2_MR_forest.pdf | ☐ |

### Week 2: Network 구축

| Day | Task | Output | 체크 |
|-----|------|--------|------|
| D6 | IGF-1R pathway genes 수집 (B-Step 1A) | IGF1R_pathway_genes.csv | ☐ |
| D6 | TSHR pathway genes 수집 (B-Step 1B) | TSHR_pathway_genes.csv | ☐ |
| D7-8 | TED disease genes 수집 (B-Step 2) | TED_disease_genes.csv | ☐ |
| D8-9 | CTD Direct Evidence 필터 확인 | 필터된 gene list | ☐ |
| D9-10 | Intersection + PPI (B-Step 3) | 3-zone gene lists | ☐ |

### Week 3: Enrichment + Off-target + Figures

| Day | Task | Output | 체크 |
|-----|------|--------|------|
| D11 | KEGG/GO enrichment (B-Step 4) | enrichment results | ☐ |
| D11 | Extension Layer (B-Step 5) | connectivity 비교 | ☐ |
| D12 | CMap (B-Step 6) | CMap results | ☐ |
| D12-13 | GTEx + hearing loss/insulin genes (C-Step 1-3) | off-target gene lists | ☐ |
| D13-14 | Off-target network (C-Step 4-5) | Fig4 draft | ☐ |
| D14 | Differential pathway figure (B-Step 7) | Fig3 draft | ☐ |

### Week 4-5: 논문 작성

| Day | Task | Output | 체크 |
|-----|------|--------|------|
| D15-16 | Methods 작성 | Methods draft | ☐ |
| D17-18 | Results 작성 | Results draft | ☐ |
| D19-20 | Introduction + Discussion | Full draft | ☐ |
| D21 | Figure 최종 정리 | Fig 1-4 final | ☐ |
| D22-23 | Supplementary 정리 | Tables S1-S5, Figs S1-S3 | ☐ |
| D24-25 | 전체 교정 + reference 정리 | Submission-ready draft | ☐ |

---

## ════════════════════════════════════════
## 위험 관리 + 대응 전략
## ════════════════════════════════════════

| 위험 | 확률 | 영향 | 대응 |
|------|------|------|------|
| TSHR eQTL IV 부족 (blood에서 발현 낮음) | 높음 | 높음 | p-threshold 완화 (5e-6), GTEx thyroid eQTL 대안, limitation 기술 |
| Graves GWAS sample size 부족 | 중간 | 높음 | FinnGen 사용, 여러 GWAS 메타 |
| OpenGWAS API 접속 불안정 | 중간 | 중간 | 로컬 데이터로 전환 (FinnGen sumstats 다운로드) |
| CTD에서 "Graves Ophthalmopathy" 유전자 부족 | 높음 | 중간 | "Graves Disease"로 확장 + DisGeNET + manual curation |
| MR 결과가 TSHR/IGF1R 모두 null | 낮음 | 높음 | Power limitation 기술, Network 분석이 주축이 됨 |
| MR 결과가 IGF1R causal + TSHR null | 낮음 | 중간 | 객관적으로 보고, discussion에서 해석 (논문은 여전히 성립) |
| Reviewer: "실험 데이터 추가하라" | 높음 | 중간 | "Hypothesis-generating computational study" 프레이밍, TSHR-ATrap in vitro는 후속 연구로 |

---

## ════════════════════════════════════════
## Reviewer 예상 질문 + 대비
## ════════════════════════════════════════

| # | 예상 질문 | 대비 전략 |
|---|----------|----------|
| 1 | "Graves disease GWAS를 TED proxy로 쓴 근거는?" | TED-specific GWAS 부재 명시, Graves 환자의 40% TED 발생, TSHR/IGF-1R pathway는 TED와 Graves 공유 |
| 2 | "Blood eQTL로 orbital tissue를 대표할 수 있나?" | Limitation 인정, eQTLGen의 large N이 보상, tissue-specific 분석은 future work |
| 3 | "Manual curation의 bias는?" | Source를 모두 기록, sensitivity analysis (manual 제외해도 결과 유지되는지) |
| 4 | "TSHR-ATrap을 직접 검증하지 않았는데?" | 논문에서 TSHR-ATrap 직접 언급 최소화, "TSHR-targeting therapy" 일반화, K-1-70 등 기존 개발 약물 언급 |
| 5 | "Computational만으로 부작용 기전 설명 가능한가?" | 기존 in vitro/animal 데이터와 일치성 검토, IGF-1 KO mouse phenotype 인용 |
| 6 | "CMap 결과의 임상 의미는?" | Drug repurposing 가능성 시사, 직접 적용은 아님 |
| 7 | "Network hub가 반드시 중요한가?" | M-LIGHT 교훈: TGFB1은 hub 아니었지만 MR 인과관계 최강 → degree만으로 판단 금지 |
| 8 | "Single IV gene의 MR 결과 신뢰도는?" | Coloc으로 보완, Wald ratio는 report하되 해석 주의, FinnGen replication |

---

## ════════════════════════════════════════
## M-LIGHT 교훈 체크리스트 (매 단계 확인)
## ════════════════════════════════════════

- [ ] CTD Direct Evidence만 필터했는가? (Inference 제외)
- [ ] Extension Layer는 intersection 밖에 배치했는가? (순환논증 방지)
- [ ] Hub gene에 과도한 의미 부여하지 않았는가?
- [ ] CMap 출처를 정확히 기술했는가? (Enrichr LINCS L1000 vs L1000CDS2)
- [ ] Coloc을 수행했는가? (single IV 비판 대비)
- [ ] FinnGen replication을 했는가?
- [ ] Null 결과를 "부정"이 아닌 "미확인"으로 해석했는가?
- [ ] 원본 데이터와 결과 수치를 대조했는가?
- [ ] β값, degree, rank 등 수치가 정확한가?

---

## Claude와의 작업 흐름

1. **각 Step 시작 전**: Claude에게 해당 Step의 스크립트 작성 요청
2. **실행 후**: 결과 공유 → Claude가 해석 + 다음 Step 안내
3. **Figure 제작 시**: Claude에게 코드 요청 (Python matplotlib/seaborn 또는 R ggplot2)
4. **논문 작성 시**: Section별로 Claude에게 영문 초안 요청 → 검토 → 수정

---

*End of Execution Plan*
