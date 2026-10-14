from pathlib import Path
import json,shutil,hashlib,datetime,html
import pandas as pd
T=Path(__file__).resolve().parent;W=T.parent;O=W.parent/'outputs/TED_Reanalysis_20261002';V=O/'Verification';V.mkdir(exist_ok=True)
s=json.loads((O/'analysis_summary.json').read_text(encoding='utf8'));f=pd.read_csv(O/'Tables/Table2_14_genes_MR_and_colocalization.csv');mr=pd.read_csv(T/'mr_v1/MR_all_genes_all_scenarios.csv');coh=pd.read_csv(O/'Tables/Table1_Datasets.csv');allr=pd.read_csv(O/'Tables/All_2544_genes_five_outcomes.csv');col=pd.read_csv(O/'Tables/Colocalization_840_settings.csv')
hits=pd.read_csv(O/'Tables/New_endpoint_exploratory_screen_hits.csv')
def pf(v):return '추정 불가' if pd.isna(v) else f'{v:.3g}'
def ci(r):return '추정 불가' if pd.isna(r.OR) else f'{r.OR:.3f} ({r.CI_lower:.3f}–{r.CI_upper:.3f})'
labels=dict(zip(coh.outcome,coh.label))
key='''<h2>MR 연관성은 확인됐지만, 공유변이 근거는 일관되지 않습니다.</h2>
<p>새 FinnGen GD에서 <b>TSHR OR 0.080, P=3.19×10⁻³²</b>였습니다. 그러나 <b>공위치화 H4=0.264, H3=0.736</b>으로, 혈액 발현과 GD가 같은 원인 변이를 공유한다는 기준을 통과하지 못했습니다.</p>
<p>기존 BBJ 발견 13개 중 7개가 같은 MR 방향과 13개 보정 기준을 충족했습니다. 평가한 14개 유전자 중 새 GD에서 H4≥0.80인 유전자는 없었습니다. <b>IGF1R는 OR 1.07, P=0.539, H4=0.010</b>이었습니다.</p>
<p>따라서 ‘동일한 발현 기전이 GD에서 재현됐다’는 결론은 지지되지 않습니다. <b>질환별 MR 연관성과 공위치화 근거가 달라지는 패턴</b>을 중심으로 원고를 재정리해야 합니다.</p>'''
data={'cohorts':coh.to_dict('records'),'counts':s['counts'],'focus':json.loads(f.to_json(orient='records',double_precision=15)),'all':json.loads(allr.to_json(orient='records',double_precision=15)),'coloc':json.loads(col.to_json(orient='records',double_precision=15))}
template=(T/'review_template.html').read_text(encoding='utf8')
regional='''<div class="box"><h2>TSHR 영역에서 무엇이 다른가</h2><p>각 질환의 가장 강한 신호와 MR 도구 SNP가 완전히 같지는 않습니다. 이 그림은 불일치의 위치를 보여주는 사후 점검이며, LD 또는 여러 원인 신호를 분리한 분석은 아닙니다.</p><img class="figure" src="Figures/Figure_R4_TSHR_regional_patterns.png" alt="TSHR 영역의 혈액 발현과 다섯 질환 GWAS 신호"><a href="Figures/Figure_R4_TSHR_regional_patterns.pdf">영역별 신호 PDF</a></div>'''
template=template.replace('</section>\n<section id="focus"',regional+'</section>\n<section id="focus"',1)
(O/'START_HERE.html').write_text(template.replace('__KEY_FINDINGS__',key).replace('__DATA__',json.dumps(data,ensure_ascii=False,allow_nan=False).replace('</','<\\/')),encoding='utf8')

parts=['''# 다섯 질환 자료 재분석 — 결과 및 연구 방향 검토

2026-10-02. 이번에 FinnGen R12 GD strict와 GO strict의 SNP 수준 MR·공위치화·민감도 분석을 새로 시행했다. 기존 원고와 제출 파일은 수정하지 않았다. 추가 분석 계획은 새 결과를 계산하기 전에 기록했으며, 전체 연구 관점에서는 사후 확장 분석이다.

## 핵심 판단

**TSHR의 강한 MR 연관성은 새 GD 자료에서도 확인되었지만, 혈액 발현과 GD가 동일한 변이를 공유한다는 근거는 재현되지 않았다.** 따라서 이번 결과를 단순히 기존 결론의 강화로 쓰면 안 된다. MR 연관성과 발현 매개 인과해석을 구분하고, 표현형별 공위치화 차이와 prior 민감도를 논문의 중심 결과에 포함해야 한다.

IGF1R는 새 FinnGen GD에서도 MR P=0.539, H4=0.010으로, 이 혈액 eQTL 분석에서 질환 감수성 관련성을 확립하지 못했다. 이것은 무효과의 증명이나 IGF-1R 치료 효과의 반박이 아니다. TSHR와 IGF1R의 P값 크기를 직접 비교한 유전자 간 효과 차이 검정도 아니다.

## 1. 어떤 자료를 분석했나

| 자료 | 질환 | 환자 | 대조군 | 역할 |
|---|---|---:|---:|---|
| BBJ | GD | 2,809 | 172,656 | 기존 발견 |
| FinnGen E4_GRAVES_STRICT | GD strict | 3,962 | 496,386 | 이번 추가; BBJ와 더 가까운 질환 정의 |
| UKB | Hyperthyroidism or thyrotoxicosis | 3,731 | 480,867 | 기존; 더 넓은 표현형 |
| FinnGen GRAVES_OPHT | Graves ophthalmopathy | 858 | 499,490 | 기존 안병증 자료 |
| FinnGen E4_GRAVES_OPHT_STRICT | Graves ophthalmopathy strict | 753 | 499,595 | 이번 추가; 표현형 민감도 |

BBJ=BioBank Japan, UKB=UK Biobank, FinnGen=핀란드의 유전체·건강등록 연구 자원이다. 질환 GWAS 파일의 한 행은 한 유전변이이며, 환자 한 명이 아니다. 혈액 발현은 별도 eQTLGen 자료에서 유전적으로 예측했다.

FinnGen 세 자료는 참여자가 겹친다. 특히 GD strict 정의는 GO 관련 코드를 포함한다. 서로 독립적인 세 차례 재현으로 세거나, OR를 비교해 GD 환자 내 TED 특이성을 추론하면 안 된다. ‘strict’는 등록자료의 endpoint 명칭이며 임상 진단의 완벽함을 의미하지 않는다. 원래 FinnGen의 858명은 실제 GWAS 표본이며, Risteys의 별도 표현형 집계와 구분했다.

공식 표본 수·정의: [FinnGen GD strict](https://r12.finngen.fi/api/pheno/E4_GRAVES_STRICT), [GO strict](https://r12.finngen.fi/api/pheno/E4_GRAVES_OPHT_STRICT), [기존 GO](https://r12.finngen.fi/api/pheno/GRAVES_OPHT), [GD 정의](https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_STRICT). 응답 원문은 Verification에 보존했다.

## 2. 1단계 — MR 연관성

동일한 2,544개 유전자와 6,135개 SNP–유전자 도구 기록을 사용했다. 새 자료마다 2,480개를 추정할 수 있었고, 64개는 도구 매칭·정렬 후 추정 불가였다. 추정 불가를 음성 결과로 세지 않았다.

| 자료 | 추정 가능 | P < 0.05 | P < 0.05/2,544 |
|---|---:|---:|---:|
''']
for r in s['counts']:parts.append(f"| {labels[r['outcome']]} | {r['estimated']:,} | {r['nominal']} | {r['screen_bonferroni']} |\n")
parts.append('''
이 표의 각 열은 서로 다른 코호트를 같은 유의성 기준으로 기술한 것이다. 원래 발견은 BBJ에서 수행했고, 새 전체 스크린은 탐색 분석이다. 표본 수·표현형·도구 수가 달라 유의 유전자 수로 유전성을 비교할 수 없다.

기존 BBJ 발견 13개에 대한 이번 FinnGen GD 확인 기준은 **같은 방향 + P < 0.05/13=0.003846**이었다. HLA-DQA2, CTLA4, TSHR, C4A, TUBB, TNFSF14, IFNGR1의 7개가 통과했다. VKORC1은 추정 불가였고, 가족의 분모는 13으로 유지했다. 이 7개가 모두 공유변이까지 확인되었다는 의미는 아니다.

## 3. 2단계 — MR와 공유변이 근거를 함께 보기

| 유전자 | 새 FinnGen GD OR (95% CI) | MR P | H3 | H4 | 13개 확인 기준 |
|---|---|---:|---:|---:|---|
''')
for _,r in f[f.outcome.eq('E4_GRAVES_STRICT')].iterrows():
    status='비교 유전자' if r.gene_symbol=='IGF1R' else '통과' if r.focused_GD_corroboration else '추정 불가' if pd.isna(r.pvalue) else '미충족'
    parts.append(f'| {r.gene_symbol} | {ci(r)} | {pf(r.pvalue)} | {r["PP.H3"]:.3f} | {r["PP.H4"]:.3f} | {status} |\n')
parts.append('''
H4의 0.000 표시는 3자리 반올림이며 정확한 0이 아니다. 원시 정밀도는 CSV에 있다. **평가한 14개 중 새 FinnGen GD에서 H4≥0.80인 유전자는 없었다.** 나머지 2,530개는 공위치화 미평가이므로 ‘모든 유전자에서 공유변이 근거가 없다’고 확장할 수 없다.

H3가 높다는 결과는 단일 원인변이 모델 안에서 서로 다른 변이 설명의 상대적 지지가 높다는 뜻이다. 실제 복수 인과변이와 LD 구조를 분리해 확정한 결과가 아니다. 특히 MHC의 상관된 유전자·복잡한 LD에서는 유전자별 원인성을 단정하지 않는다.

## 4. TSHR의 해석이 어떻게 달라지는가

| 질환 자료 | OR (95% CI) | MR P | H4 기본 prior | H4 낮은 prior |
|---|---|---:|---:|---:|
''')
for _,r in f[f.gene_symbol.eq('TSHR')].iterrows():parts.append(f'| {labels[r.outcome]} | {ci(r)} | {pf(r.pvalue)} | {r["PP.H4"]:.3f} | {r.H4_low_prior:.3f} |\n')
parts.append('''
TSHR의 MR 방향은 모두 같지만, 공유변이 근거는 같지 않다. 새 GD는 H3=0.736/H4=0.264이고, UKB도 기존 H3 우세 패턴이었다. 기존 GO의 H4=0.986은 GO strict에서 0.891로 낮아졌고, GO strict는 낮은 prior에서 0.451로 기준 아래로 내려갔다. **표현형별 결과 차이는 관찰되지만, 이 차이가 생물학적 TED 특이성 때문인지 표본 정밀도·사례 정의·LD/복수 신호 때문인지는 이 분석만으로 분리할 수 없다.**

영역 점검에서 발현의 도구 SNP는 rs179252, 새 GD의 분석 영역 내 최대 |Z| 변이는 rs1023586이었다. 두 위치는 같은 TSHR 영역 안에 있다. Figure R4가 신호 분포를 보여준다. 서로 다른 최상위 SNP라는 사실만으로 서로 다른 원인 변이를 확정할 수는 없으며, 이 그림은 H3 결과를 설명하기 위한 사후 진단이다. 영역 그림의 P값은 beta/SE에서 계산한 정규 근사이고, FinnGen 원본의 보정 검정 P와 항상 동일하지 않을 수 있다.

**‘보호 방향’은 유전적으로 예측한 혈액 TSHR 발현 증가와 낮은 질환 odds의 연관을 뜻한다.** 실제 안와·갑상선 발현을 올리면 예방된다는 뜻이 아니다. 이번 결과는 흉선 면역관용 기전을 입증하지 않으며, 새 GD 공위치화 불일치 때문에 혈액 발현이 인과적으로 GD를 억제한다고 단정하기는 더 어렵다.

## 5. IGF1R와 CTLA4

IGF1R의 새 GD 결과는 OR 1.070 (0.861–1.330), P=0.539, H4=0.010이었다. GO strict는 OR 1.516 (0.924–2.486), P=0.0996, H4=0.038이었다. 따라서 ‘혈액 eQTL 기반 질환 감수성 근거를 확립하지 못했다’가 적절하다. GD의 CI가 비교적 좁아진 것은 정보가 늘었다는 뜻이지만, 작은 효과·다른 조직의 작용·치료반응의 부재를 증명하지 않는다. 통합 요약통계로 치료 기전을 직접 시험하지 않았다.

CTLA4의 새 GD MR는 OR 0.173, P=2.15×10⁻⁶이었지만 H4=0.254/H3=0.746이었다. 새 GO strict는 H4=0.979로 높았으며 낮은 prior에서도 0.820이었다. 그러나 MR는 rs13030124 의존성이 컸다. 새 GD에서 이 SNP를 제외하면 남은 SNP는 OR 0.570, P=0.323이고, 반대로 rs148849825를 제외하면 OR 0.154, P=2.17×10⁻²⁶이었다. GD 이질성 P=0.0275와 함께 해석해야 한다. GO strict에서도 같은 도구 의존성이 관찰됐다.

TSHR은 IV 하나여서 leave-one-out을 할 수 없다. 14개 중 다중 IV 유전자에서 시행한 새 두 자료의 leave-one-out **42행**을 모두 제공했다.

## 6. 민감도에서 유지된 것과 바뀐 것

eQTLGen 빈도로 대체해도 노출 Z의 최대 절대 변화는 1.42×10⁻¹⁴이었다. 동일 SNP를 고정한 단일 IV의 P 변화는 3.89×10⁻¹⁶ 이하였다. TSHR·IGF1R의 핵심 판정과 13개 확인 가족의 7개 통과 목록은 유지되었다.

공위치화의 같은 변이 집합에서 빈도만 바꾸었을 때 최대 |ΔH4|는 새 GD 0.000129, 새 GO strict 0.000215였다. 변이 가용성 변화까지 포함하면 각각 0.001914, 0.010975였다. 14개에서 H4=0.80을 넘나드는 유전자는 없었다.

**전체 MR 스크린은 완전히 불변이 아니었다.** 새 GD의 C4B는 P=1.10×10⁻⁵에서 1.48×10⁻⁴로, HLA-G는 P=0.107에서 1.72×10⁻⁵로 바뀌어 2,544개 보정 기준의 통과 여부가 달라졌다. 두 변화는 같은 SNP 집합에서 빈도에 따른 노출 스케일·가중치 변화만으로도 나타났다. 따라서 이 유전자를 강건한 새 발견이라고 보고해서는 안 된다. 전체 정렬을 다시 하면 각 새 자료에서 9개 유전자의 IV 수가 바뀌었다. 명목 유의성 변동은 GD 17개, GO strict 12개였다.

## 7. 탐색에서 새로 보인 유전자를 어떻게 다룰까

새 GD에서 2,544개 보정 기준을 통과한 13개는 BBJ의 13개와 동일한 목록이 아니다. 기존 BBJ 발견 목록 밖에서 APOM, NEU1, IGLV2-18, IGLV3-12, TNF, HIST1H4E, LTA, LTB, C4B가 나타났다. 이 9개에는 이번 고정 14개 공위치화 분석을 적용하지 않았다.

APOM과 NEU1는 동일 도구 rs3101018을 공유하고, IGLV2-18와 IGLV3-12는 rs5759410을 공유한다. 같은 극소 P값이 반복되는 것은 독립적인 네 유전자 발견을 의미하지 않는다. `Exploratory_hits_retained_instruments.csv`로 실제 SNP 구성을 제공한다. 후속 후보를 다룬다면 먼저 공유 SNP·유전자좌로 묶고, 별도의 사후 확장으로 전체 영역 공위치화를 평가해야 한다. 현재 단계에서 이 후보를 검증된 원인 유전자나 약물 표적으로 부를 수 없다.

## 8. 권장 논문 방향과 수정할 부분

연구 질문은 **‘약물표적이 될 수 있는 유전자의 혈액 발현 관련 변이가 GD·갑상선기능항진증·안병증에서 어떤 연관성과 공유변이 근거를 보이는가’**로 잡는 편이 자료에 맞는다. TSHR–IGF1R 대조는 임상적 동기와 주요 예시로 유지할 수 있지만, ‘TSHR susceptibility 대 IGF1R effector’라는 역할 분할을 이 자료가 입증한 것처럼 쓰면 안 된다.

- Methods: 추가한 두 endpoint와 겹치는 FinnGen 표본, 사후 확장, 13개 확인 가족 기준을 명시한다.
- Results: 7개 MR 확인과 14개 중 새 GD H4 통과 0개를 함께 제시한다. TSHR의 GD H3 우세, GO strict prior 민감도, 새 후보의 공유 SNP와 미평가 상태를 보여준다.
- Discussion: ‘MR 방향의 일관성’과 ‘혈액 발현의 인과매개 지지’가 다름을 설명한다. 통계적으로 유의한 한 SNP의 연결만으로 발현 기전을 확정하지 않는다. IGF1R 치료 효과에 대한 장황한 추측보다 현재 유전 자료가 답한 범위에 집중한다.
- Table: 자료표를 5개 endpoint로 확장하고, MR·H3·H4·낮은 prior·도구 수를 통합해 읽도록 한다. 전체 12,720행은 상세자료에 둔다.
- Figure: 검토용 R1(14개×5개 H4), R2(주요 MR forest), R3(TSHR prior), R4(TSHR 영역)를 제공했다. 기존 제출 Figure는 아직 교체하지 않았다.

가능한 결론 문장: *TSHR showed consistent MR associations across disease datasets, but colocalization support varied across phenotypes and prior assumptions. A shared causal relationship between blood TSHR expression and GD susceptibility was not consistently established. IGF1R lacked convergent MR and colocalization evidence in these blood-based analyses.*

제목 방향 예시: *Blood-expression Mendelian randomization and colocalization across Graves disease and ophthalmopathy phenotypes: a druggable-gene analysis*. 이는 검토 제안이며 최종 제목으로 반영하지 않았다.

## 9. 계산 검증과 남아 있는 범위

두 새 원자료는 각각 21,327,062행과 21,327,061행이었다. 유한 beta, 양의 SE, P·빈도 범위 및 대립유전자 누락 점검은 통과했다. 이는 INFO, 독립 LD, 모든 생물학적 가정의 검증을 대신하지 않는다. rsID가 없는 행과 쉼표 별칭은 기존 연결 규칙에서 제외됐다. 다중 대립유전자 행은 보존하고 적합성 확인 후 사용했으며, 최종 MR·공위치화 내 중복 rsID는 없었다.

기존 FinnGen MR 2,480개를 재현했다. 기존 FinnGen 공위치화 168개 조건의 차이는 최대 1.11×10⁻¹⁶이었다. 새 336개와 대조 168개를 합친 공위치화 504개 조건은 설치된 R coloc과 최대 5.53×10⁻¹⁴ 차이로 일치했다. 저장된 SNP별 입력에서 Python으로 독립 계산한 MR 12,402개도 최대 1.47×10⁻¹⁴ 차이로 일치했다.

R coloc은 일부 영역에서 약한 질환 신호 경고를 냈다. beta/SE 및 분산 입력을 대조했고 숫자 오류는 없었다. 약한 신호 자체는 음성 확정으로 바꾸지 않았으며 경고 원문을 보존했다. 최초 추출은 rsID 중복 금지 검사가 다중 대립유전자를 차단해 중단됐고, 원인을 확인해 버전 2에서 수정했다. 계획·중단·수정·실행 기록은 Verification에 있다.

이 단계는 계획한 두 FinnGen 확장 분석의 완료이지 전체 논문의 제출 준비 완료 판정은 아니다. 새 후보 9개 공위치화, 다중 신호 분리, GD 환자 내 TED 비교는 미시행이다. 기존 UKB의 상세 임상 코드·조상집단 선정 문헌 확인도 별도 남아 있다. **현재 결과를 원고에 반영하기 전에는 기존 ‘submit-ready’ 판단을 그대로 유지할 수 없다.**
''')
(O/'REPORT_KO.md').write_text(''.join(parts),encoding='utf8')
shutil.copy2(T/'DATA_DICTIONARY.md',O/'DATA_DICTIONARY.md')
for src in [T/'ANALYSIS_PLAN.md',T/'plan_freeze.json',T/'DEVIATIONS.md',T/'diagnostic_plan.json',T/'analysis.log',T/'inputs_v2/extraction_manifest.json',T/'mr_v1/verification.csv',T/'coloc_v1/R_coloc_verification.csv',T/'coloc_v1/R_warnings.txt',T/'coloc_v1/python_verification.json']:
    shutil.copy2(src,V/src.name)
for group in ['mr_v1','coloc_v1']:
    for name in ['run.log','sessionInfo.txt','R_sessionInfo.txt','warnings.txt']:
        p=T/group/name
        if p.exists():shutil.copy2(p,V/f'{group}_{name}')
for name in ['E4_GRAVES_STRICT_api.json','E4_GRAVES_OPHT_STRICT_api.json','GRAVES_OPHT_api.json','FinnGen_R12_strictGD_definition_text.txt','FinnGen_R12_strictGO_definition_text.txt']:
    shutil.copy2(W.parent/'outputs/TED_Data_Review_20261002/provenance'/name,V/name)
scripts=V/'scripts';scripts.mkdir(exist_ok=True)
for p in T.iterdir():
    if p.suffix in ['.py','.R','.html']:shutil.copy2(p,scripts/p.name)
verify={'status':'PASS for the specified analyses, not submission readiness','primary_new_rows':5088,'frequency_scenario_new_rows':20352,'five_outcome_primary_rows':12720,'focused_genes':14,'new_coloc_settings':336,'old_coloc_identity_settings':168,'R_coloc_vectors':504,'R_coloc_max_delta':float(pd.read_csv(T/'coloc_v1/R_coloc_verification.csv').set_index('test').loc['maximum_absolute_posterior_difference','value']),'independent_python_MR_checks':s['independent_python_MR_checks'],'independent_python_MR_max_delta':s['independent_python_MR_max_difference'],'leave_one_out_rows':42,'plan_hash_unchanged':hashlib.sha256((T/'ANALYSIS_PLAN.md').read_bytes()).hexdigest()==json.loads((T/'plan_freeze.json').read_text())['plan_sha256'],'limits':['Fixed 14-gene coloc only','Overlapping FinnGen endpoints','Single-causal-variant coloc model','No within-GD TED comparison','No clinical covariate subgroup analysis']}
assert verify['plan_hash_unchanged'] and len(mr[mr.outcome.ne('FinnGen_GO')])==20352
(V/'verification.json').write_text(json.dumps(verify,indent=2,ensure_ascii=False),encoding='utf8')
sources=[T/'inputs_v2/E4_GRAVES_STRICT_selected.tsv',T/'inputs_v2/E4_GRAVES_OPHT_STRICT_selected.tsv',T/'inputs_v2/eqtl_14genes.tsv',W/'maf_sensitivity_20260905/inputs/instruments_verified.csv',W/'maf_sensitivity_20260905/inputs/eur_freq_selected.tsv',W/'maf_sensitivity_20260905/inputs/eqtlgen_af_selected_v1.csv',W/'expanded_analysis_20260930/inputs/eur_freq.tsv',W/'expanded_analysis_20260930/inputs/eqtl_af.tsv']
(V/'input_hashes.json').write_text(json.dumps([{'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources],indent=2),encoding='utf8')
print(json.dumps(verify,indent=2,ensure_ascii=False))
