# 재분석 표 읽는 법

이 표는 분석 확장 검토용이다. 기존 원고의 Table 번호나 Figure 번호를 바꾸지 않는다.

## 파일

| 파일 | 내용 |
|---|---|
| Table1_Datasets.csv | 질환 정의·출처·환자/대조군·자료 역할 |
| All_2544_genes_five_outcomes.csv | 2,544개 유전자 × 5개 질환, 추정 불가도 포함한 12,720행 |
| Table2_14_genes_MR_and_colocalization.csv | 기존 BBJ 발견 13개 + IGF1R, 총 70행의 MR·공위치화 통합표 |
| Table3_BBJ_discovery_followup_in_FinnGen_GD.csv | 14개 비교 유전자의 새 FinnGen GD 결과; 확인 분석 가족은 BBJ 발견 13개 |
| New_endpoint_exploratory_screen_hits.csv | 새 두 자료의 전체 스크린에서 P < 0.05/2,544를 만족한 결과. 사후 탐색이며 공위치화 미평가 유전자가 포함될 수 있음 |
| New_MR_all_frequency_scenarios.csv | 새 두 자료 × 2,544개 × 네 빈도/정렬 조건 |
| Colocalization_840_settings.csv | 기존 504개 + 새 336개 공위치화 조건. 독립 실험 840회라는 의미가 아님 |
| Leave_one_out_new_outcomes.csv | 14개 중 다중 IV 유전자의 새 두 자료에서 SNP 하나씩 제외한 MR |
| MR_instrument_coverage.csv | 선택·매칭·보존·제외 도구변수 수 |
| Colocalization_variant_coverage.csv | 전체 eQTL 및 질환 매칭·대립유전자 적합·빈도 가용 변이 수 |

## 주요 열

- `gene_symbol` / `gene`: 유전자 기호.
- `outcome`: 질환 GWAS 식별자. 치료 결과가 아니다. BBJ_Graves, UKB_hyperthyroid, FinnGen_GO는 기존 자료, E4_GRAVES_STRICT 및 E4_GRAVES_OPHT_STRICT는 이번 추가 자료.
- `n_iv`: 대립유전자 정렬 후 실제 MR에 남은 도구 SNP 수. 0 또는 누락이면 해당 조건에서 추정 불가이며 음성 결과가 아니다.
- `beta`, `se`: 유전적으로 예측한 발현 증가와 질환 log-odds의 연관 추정치 및 표준오차.
- `OR`, `CI_lower`, `CI_upper`: odds ratio와 95% 신뢰구간. 약 1 SD 혈액 발현 단위. 발현량을 실제로 치료로 변화시킨 효과가 아니다.
- `pvalue`: 양측 MR P값. 원 GWAS의 단일 SNP P값과 구분한다. 0.05 기준, 0.05/13 기준, 0.05/2,544 기준의 목적은 서로 다르다.
- `Q`, `Q_df`, `Q_pvalue`: 다중 IV의 Cochran 이질성 통계·자유도·P값. IV 하나에는 해당하지 않는다.
- `BBJ_discovery`: 기존 BBJ 2,544개 스크린의 보정 기준을 통과했는가.
- `same_direction_BBJ`: BBJ와 MR beta 부호가 같은가. 이것만으로 재현 또는 군간 차이를 입증하지 않는다.
- `focused_GD_corroboration`: 이번에 정한 13개 확인 가족에서 같은 방향이고 FinnGen GD P < 0.05/13인가. IGF1R는 이 가족의 발견 유전자가 아니다.
- `PP.H0`–`PP.H4`: 사후확률. 본 파일의 dataset1은 질환, dataset2는 발현이다. H0=둘 다 없음, H1=질환만, H2=발현만, H3=둘 다 있으나 서로 다른 변이, H4=공유변이. 다섯 값의 합은 1이다. H4는 P값이 아니다.
- `p12`: SNP 하나가 두 형질에 함께 작용한다는 사전확률. 결과의 H4와 구분한다.
- `n_overlap`: 해당 조건의 공위치화에 실제 사용한 영역 변이 수. MR IV 수와 다르다.
- `top_snp`: H4를 조건으로 했을 때 최대 SNP posterior를 갖는 변이. 검증된 유일한 원인 변이라는 뜻이 아니다.
- `H4_low_prior`: p12=10⁻⁶, 기존 EUR 빈도에서 PP.H4.
- `H4_eqtlgen`: p12=10⁻⁵, 가용 eQTLGen 빈도에서 PP.H4.
- `removed_snp`: leave-one-out에서 제외한 SNP. 원래 IV가 2개라면 남은 1개의 Wald ratio이다.

## 민감도 조건

| scenario | 의미 |
|---|---|
| original_reference | 기존 EUR 참조 빈도, 기본 분석 |
| paired_reference | 두 빈도를 모두 쓸 수 있는 동일 변이 집합에서 참조 빈도 |
| paired_eqtlgen | 위와 동일한 변이 집합에서 eQTLGen 빈도만 대체 |
| reharmonized_eqtlgen | MR: eQTLGen 빈도로 대립유전자를 다시 정렬하고 제외 기준을 다시 적용 |
| available_eqtlgen | 공위치화: eQTLGen 빈도가 가용한 적합 변이를 모두 사용 |

MR의 paired 조건은 기본 분석에서 보존된 SNP를 고정한다. 공위치화의 paired 조건은 두 빈도가 유효한 공통 영역 변이를 고정한다. 빈도만의 변화와 변이 가용성/정렬 변화가 섞이지 않게 구분한 것이다.

eQTL 좌표는 GRCh37, FinnGen 좌표는 GRCh38이다. 직접 위치를 서로 비교해 병합하지 않았으며 rsID·염색체·대립유전자를 대조했다. 하나의 rsID에 여러 대립유전자 기록이 있는 경우 적합성 확인 후 사용하며 P값으로 행을 고르지 않았다. 쉼표로 연결된 별칭 rsID는 기존 추출 규칙에 맞춰 제외했다.
