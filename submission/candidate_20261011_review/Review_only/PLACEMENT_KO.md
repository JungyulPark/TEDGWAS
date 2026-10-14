# 조직별 보강 결과의 표와 그림 및 논문 배치

권고는 **본문 Results의 독립된 보강분석 소절에서 핵심 수치를 보고하고, Table S3와 상세데이터로 뒷받침하는 것**입니다. 짧은 보조설명 한 문장으로만 처리하기에는 연구 질문과 직접 관련된 결과입니다. 다만 메인 Table 2나 Figure 2에 기존 혈액 결과와 섞어 넣거나, 새 주결과 표와 그림을 모두 추가할 정도의 독립적인 검증으로 제시하는 것은 권하지 않습니다.

## 어떤 답이 추가됐는가

TSHR의 낮은 GD odds 방향이 두 지방조직 eQTL에서도 유지됐고, FinnGen GD와 BBJ GD에서 보정 후 유의했습니다. 이는 혈액 표현형에만 한정된 관찰을 보강합니다. 그러나 피하지방과 내장지방의 안병증 결과가 모두 같은 수준의 근거를 보이지는 않았습니다.

IGF1R은 갑상선과 배양 섬유모세포에서 분석 가능했지만 10개 질병 비교 모두 P < 0.05에 도달하지 않았습니다. 따라서 제한적인 근거를 단순히 ‘혈액이라는 조직 선택 때문’이라고 단정하기 어려워졌습니다. 반대로 신뢰구간, 단일 도구변수, 비안와 조직이라는 조건 때문에 생물학적 역할이나 치료 효과를 부정할 수도 없습니다.

가용성 표와 그림의 전혈은 **GTEx 전혈(N = 670)**입니다. 기존 혈액 MR은 **eQTLGen(최대 31,684명)**을 사용했습니다. 따라서 GTEx 전혈에서 고정 임계값을 충족하는 도구변수가 없다는 표시는 기존 혈액 MR을 수행할 수 없었다는 뜻이 아닙니다.

## 권장 배치

| 위치 | 배치할 내용 | 판단 |
|---|---|---|
| 본문 Methods | 고정한 5개 조직, 도구변수 기준, 일차 10개·이차 40개 보정, 분석 불가 처리 | 필수. 이전 통합 원고에 반영됨 |
| 본문 Results | Tissue-specific extension 소절에서 네 일차 추정치, BBJ GD 확인, 안병증의 약한 결과, IGF1R 비유의 결과를 한 단락으로 보고 | 필수. 단순한 Discussion 보조설명으로만 남기지 않음 |
| 본문 Discussion | TSHR의 방향 일치가 혈액 밖으로 확장됐다는 의미와 IGF1R 해석의 변화 | 1개 논점으로 반영. 지방조직 매개 인과기전으로 단정하지 않음 |
| Table S3 | 열 개 조직–유전자 조합의 가용성 및 FinnGen GD 일차 결과 | 논문에 포함. 현재 S3 유지 |
| Supplementary Dataset | 전체 50개 계획 비교, 20개 추정치, 30개 분석 불가 사유, 도구변수·대립유전자·원 수치 | 하나의 기존 워크북 Tissue_extension 시트로 제공 |
| 메인 Table 1–3 | 기존 혈액 분석 및 발견 유전자 결과 | 새 조직 결과를 억지로 합치지 않고 현행 유지 |
| 메인 Figure 1–3 | 기존 연구 설계·혈액 MR·공위치화 결과 | 현행 유지 |
| 새 요약 그림·전체 forest plot | 가용성과 결과를 저자가 빠르게 검토하거나 발표에서 설명 | 검토·발표용 제공. 논문에서는 현재 S3와 상세데이터로 충분하므로 Figure S3를 의무적으로 추가하지 않음 |

이 배치는 보충표 3개·보충그림 2개·상세데이터 파일 1개라는 기존 구성을 유지합니다. 추가 검토표를 Table S4로 새로 번호 붙이지 않습니다. 이번 패키지는 시각화와 배치 판단 자료이며, 이전 통합 원고와 데이터의 수치 및 파일을 변경하지 않았습니다.

## 주결과와 구별하는 이유

1. 본래 혈액 분석 결과를 알고 수행한 후속 확장입니다. 조직 선택은 고정했지만 원 연구 전체의 사전등록 분석은 아닙니다. 이번에 부르는 ‘일차 비교’는 **조직 보강 내부의 FinnGen GD 비교**입니다.
2. 네 조합 모두 SNP 1개를 사용했습니다. MR-Egger, 이질성 Q, leave-one-out으로 조직별 결론을 추가 검증할 수 없습니다.
3. 조직별 공위치화는 전체 영역 자료를 확보하지 못해 시행하지 않았습니다. 같은 원인 변이가 발현과 질병을 함께 설명하는지 확인된 결과로 표시할 수 없습니다.
4. 기존 질병 GWAS를 재사용했고, 비슷한 유전체 영역의 변이를 이용합니다. 독립된 질병 코호트와 별개의 원인 변이로 재현한 결과로 세면 안 됩니다.
5. GTEx의 갑상선·배양 섬유모세포·피하/내장 지방·전혈은 안와조직이 아닙니다. GD 환자 내 TED 발생이나 치료 반응을 확인한 분석도 아닙니다.

이 이유들은 결과를 숨길 근거가 아니라, **보고는 충분히 하되 주장 수준을 보강분석에 맞춰야 하는 이유**입니다.

## 이번에 정리한 표와 그림

- Tissue_result_tables.docx / .pdf: A에서 10개 조합의 표본 수·도구변수·F·최소 eQTL P·분석 가능 여부, B에서 전체 20개 OR·95% CI·원 P·보정 P를 제시합니다. 보정 후 유의한 P만 굵게 표시했습니다.
- Table_S3_current.pdf: 논문에 이미 반영된 간결한 Table S3입니다. 새 상세 검토표와 역할을 구분합니다.
- Figures/Tissue_primary_summary.pdf / .png: A는 10개 조합의 가용성, B는 조직 보강의 일차 FinnGen GD 네 결과입니다. 색 있는 가용성 상자는 질병 유의성이 아니라 도구변수 확보를 뜻합니다. B의 채운 사각형과 별표는 다중검정 보정 후 유의함을 뜻합니다.
- Figures/Tissue_all_estimates.pdf / .png: 20개 추정치를 빠짐없이 보여주는 forest plot입니다. 피하지방·내장지방 안병증 결과의 차이도 포함했습니다.
- Source_data: 형식화 전 원 수치와 분석 불가 조합을 보존했습니다. 기존 Supplementary Dataset을 대체하는 새 제출 파일은 아닙니다.

## 본문에서 유지할 핵심 표현

> In the fixed-panel tissue extension, higher genetically proxied TSHR expression in subcutaneous and visceral adipose was associated with lower FinnGen GD odds, with family-corrected support also observed in BBJ GD. Ophthalmopathy support was less consistent for visceral adipose. Neither thyroid nor cultured-fibroblast IGF1R estimates reached P < 0.05 in any disease dataset. These single-instrument analyses extend the association evidence beyond blood but do not establish tissue-mediated causality.

수치가 포함된 상세 Results 단락은 이전 통합 원고의 Tissue-specific extension 소절에 이미 있습니다. 위 문장은 새 독립 결론을 추가하는 것이 아니라, 그 결과가 논문에서 갖는 역할을 요약합니다.
