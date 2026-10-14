from pathlib import Path
import json,csv,html,math,hashlib,zipfile
import pandas as pd
ROOT=Path(__file__).resolve().parents[2];O=ROOT/'outputs/TED_Data_Review_20261002';M=O/'provenance'
A=json.loads((O/'analysis_summary.json').read_text(encoding='utf8'))
inv=json.loads((M/'local_file_inventory.json').read_text(encoding='utf8'))
mr=pd.read_csv(O/'All_2544_genes_3_outcomes.csv');sel=pd.read_csv(O/'Selected_14_gene_results.csv')
settings=pd.read_csv(ROOT/'work/expanded_analysis_20260930/results/Stage2_all_504_colocalization_settings.csv')
labels={'BBJ_Graves':'BBJ · Graves병','UKB_hyperthyroid':'UKB · 갑상선기능항진증','FinnGen_GO':'FinnGen · Graves 안병증'}
cohorts=[
 {'name':'BBJ','full':'BioBank Japan','accession':'GCST90018627','disease':'Graves병 (GD)','cases':2809,'controls':172656,'ancestry':'동아시아 · 일본','role':'현재 발견 분석','scale':'질환 log-OR','build':'사용 가공본 GRCh38; 원본 GRCh37','definition':'GD 환자와 연구의 비사례 대조군. 안병증 동반 여부를 나누는 비교가 아님.','source':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018627'},
 {'name':'UKB','full':'UK Biobank','accession':'GCST90038636','disease':'Hyperthyroidism or thyrotoxicosis','cases':3731,'controls':480867,'ancestry':'영국 모집; Catalog의 조상집단 표기는 NR','role':'더 넓은 질환 정의 비교','scale':'원자료 BOLT-LMM β → 근사 log-OR 변환','build':'사용 가공본 GRCh38; 원본 GRCh37','definition':'GD로 한정되지 않는 질환 명칭. 이번 검토에서 세부 진단 코드·조상집단 선정은 원논문 전문과 추가 대조가 필요함.','source':'https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636'},
 {'name':'FinnGen','full':'FinnGen · Finland','accession':'R12 GRAVES_OPHT','disease':'Graves ophthalmopathy','cases':858,'controls':499490,'ancestry':'핀란드','role':'안병증 표현형의 보조 평가','scale':'질환 log-OR','build':'GRCh38','definition':'안와·눈물기관 관련 조건과 갑상선중독증 진단을 조합한 레지스트리 정의. 대조군은 해당 안병증 사례가 아닌 사람들로, GD 환자로 제한되지 않음.','source':'https://r12.finngen.fi/api/pheno/GRAVES_OPHT'}
]
alternatives=[]
for code in ['E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT']:
 j=json.loads((M/(code+'_api.json')).read_text(encoding='utf8'))
 alternatives.append({'accession':code,'disease':j['phenostring'],'cases':j['num_cases'],'controls':j['num_controls'],'file':'C:/ProjectTEDGWAS/finngen_R12_'+code+'.gz','status':'로컬 파일 헤더 + R12 공식 정의/표본 수 확인. 이번에 MR·공위치화 재실행하지 않음.','source':'https://r12.finngen.fi/api/pheno/'+code})
column_defs={
 'rsid':'유전변이 식별자. 사람의 ID가 아님.','hm_rsid':'조화 처리된 rsID. 다른 효과값 열과 대립유전자 방향을 함께 확인.','rsids':'변이의 rsID. 일부 행은 없음; 여러 식별자가 포함될 수 있음.','variant_id':'원본 변이 식별자. 파일마다 좌표 표기 또는 rsID; 빈 값일 수도 있음.','hm_variant_id':'조화 처리된 변이 ID.',
 'chromosome':'염색체 번호.','hm_chrom':'조화 처리된 염색체 번호.','#chrom':'염색체 번호.','base_pair_location':'이 파일의 변이 좌표. GRCh37 원본과 GRCh38 가공본을 구분.','hm_pos':'조화 처리된 변이 좌표.','pos':'GRCh38 변이 좌표.',
 'effect_allele':'β가 표현하는 효과 대립유전자.','hm_effect_allele':'조화 처리 후 β가 표현하는 효과 대립유전자.','other_allele':'비교 대립유전자.','hm_other_allele':'조화 처리된 비교 대립유전자.','alt':'대체 대립유전자; FinnGen β의 기준.','ref':'참조 대립유전자.',
 'beta':'변이 한 사본 증가와 질환의 연관 효과. BBJ/FinnGen은 log-OR, UKB는 선형모형 척도. 유전자 발현의 MR β와 구별.','hm_beta':'조화 처리 후 β. 대립유전자가 바뀌면 부호도 바뀔 수 있음.','standard_error':'변이 효과 β의 표준오차(SE).','sebeta':'변이 효과 β의 표준오차(SE).','p_value':'해당 변이–질환 연관성 P값. 유전자 단위 MR P값과 다름.','pval':'해당 변이–질환 연관성 P값.','mlogp':'−log10(P). 값이 클수록 작은 P값.',
 'effect_allele_frequency':'효과 대립유전자 빈도. 환자 비율이 아님.','hm_effect_allele_frequency':'조화 처리된 효과 대립유전자 빈도.','af_alt':'전체 분석 표본의 대체 대립유전자 빈도.','af_alt_cases':'환자군의 대체 대립유전자 빈도. 개인별 자료가 아님.','af_alt_controls':'대조군의 대체 대립유전자 빈도.','nearest_genes':'위치상 가까운 유전자. 원인 유전자라는 판정이 아님.',
 'hm_code':'GWAS Catalog 조화 처리 상태 코드.','hm_coordinate_conversion':'좌표 변환 처리 표기.','odds_ratio':'원본 OR 열. UKB 가공본에서는 비어 있음.','hm_odds_ratio':'조화 처리된 OR 열; UKB에서는 비어 있음.','ci_lower':'원본 신뢰구간 하한; UKB에서는 비어 있음.','ci_upper':'원본 신뢰구간 상한; UKB에서는 비어 있음.','hm_ci_lower':'조화 처리된 신뢰구간 하한; UKB에서는 비어 있음.','hm_ci_upper':'조화 처리된 신뢰구간 상한; UKB에서는 비어 있음.',
 'INFO':'변이 imputation 품질 지표. UKB 원본에는 있고 현재 가공본에는 없음.','GENPOS':'유전지도 위치.','CHISQ_LINREG':'선형회귀 검정통계량.','P_LINREG':'선형회귀 P값.','CHISQ_BOLT_LMM_INF':'BOLT-LMM infinitesimal 모형 검정통계량.','P_BOLT_LMM_INF':'BOLT-LMM infinitesimal 모형 P값.','CHISQ_BOLT_LMM':'BOLT-LMM 검정통계량.',
 'Pvalue':'해당 SNP–유전자 발현 연관성 P값.','SNP':'발현 연관 변이 rsID.','SNPChr':'변이 염색체.','SNPPos':'GRCh37 변이 위치.','AssessedAllele':'Zscore가 표현하는 대립유전자.','OtherAllele':'비교 대립유전자.','Zscore':'SNP–혈액 유전자 발현 연관성의 방향과 통계량.','Gene':'Ensembl 유전자 ID.','GeneSymbol':'유전자명.','GeneChr':'유전자 염색체.','GenePos':'GRCh37 유전자 위치.','NrCohorts':'해당 연관성에 기여한 코호트 수.','NrSamples':'해당 SNP–유전자 연관성의 표본 수; 행마다 다름.','FDR':'발현 연관성의 FDR 보정값.','BonferroniP':'발현 연관성의 Bonferroni 보정 P값.'}
columns=[]
for f in inv:
 for k in f['columns']:columns.append({'file':Path(f['file']).name,'column':k,'meaning_ko':column_defs.get(k,'원본 문서와 추가 대조 필요'),'first_value':f['first_row'].get(k,'')})
pd.DataFrame(columns).to_csv(O/'Raw_file_column_dictionary.csv',index=False,encoding='utf-8-sig')
pd.DataFrame(cohorts).to_csv(O/'Cohort_comparison.csv',index=False,encoding='utf-8-sig');pd.DataFrame(alternatives).to_csv(O/'Additional_available_outcomes.csv',index=False,encoding='utf-8-sig')
def records(x):return x.astype(object).where(pd.notna(x),None).to_dict(orient='records')
payload={'analysis':A,'cohorts':cohorts,'alternatives':alternatives,'files':inv,'columns':columns,'genes':records(mr),'selected':records(sel),'low_prior':records(settings[(settings.scenario=='original_reference')&settings.p12.eq(1e-6)]),'labels':labels}
report='''# 연구 방향 재점검 — 1단계: 자료의 정체와 분석 가능 범위

2026-10-02. 기존 결론을 전제로 하지 않고 현재 사용한 3개 GWAS의 전체 행을 읽어 자료 품질을 점검했다. 기존 2,544개 유전자의 MR 결과는 다시 집계하고 β·SE로 P값, OR, CI의 산술 일치를 확인했다. SNP 수준 MR·공위치화를 새로 재실행하거나 원고를 수정한 단계는 아니다.

## 1. 가장 중요한 구분

BBJ, UKB, FinnGen은 연구 자원 이름이다. 이 연구가 다운로드한 것은 각각 특정 질환 한 가지의 GWAS 요약통계이다. 바이오뱅크 전체 자료나 환자별 진료기록이 아니다. 질환 GWAS의 한 행은 유전변이 한 개이며, eQTLGen의 한 행은 SNP–유전자 발현 연관성이다. 현재 파일에서 나이·성별·흡연·TRAb·CAS·안구돌출·복시·치료반응의 개인별 값을 볼 수 없고, 이를 보정한 임상 회귀나 임의 subgroup을 새로 만들 수 없다. 원 GWAS가 공변량을 보정했는지와 우리에게 그 개인별 변수가 있는지는 별개의 문제다.

| 자료 | 실제 질환 | 환자 | 대조군 | 현재 용도 |
|---|---|---:|---:|---|
'''
for c in cohorts:report+=f"| {c['name']} ({c['accession']}) | {c['disease']} | {c['cases']:,} | {c['controls']:,} | {c['role']} |\n"
report+='''
BBJ는 동아시아 GD 자료이다. UKB는 GD에 한정되지 않는 더 넓은 진단명이므로 동일 표현형의 반복 검증이라고 부르기 어렵다. UKB의 GWAS Catalog 조상집단 항목은 NR(not reported)이며, 영국 모집이라는 정보만으로 유럽계 선정 기준을 확정하지 않았다. 원논문 전문의 선정 기준을 대조하는 작업이 남아 있다. 이 접근 제한은 기존 수치가 틀렸다는 뜻이 아니다.

FinnGen 안병증 대조군은 GD 환자로 제한되어 있지 않다. 따라서 이 비교에는 GD 감수성과 안병증 과정이 함께 반영될 수 있으며, GD 환자 중 누가 TED를 갖게 되는지를 직접 분리하지 못한다. Risteys R12의 표현형 집계 894명과 실제 GWAS API의 858명은 서로 다른 집계이며, 사용한 GWAS에는 858/499,490을 적용한다. 이번에 공식 API에서 다시 확인했다.

출처: [BBJ GWAS Catalog API](https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90018627), [UKB GWAS Catalog API](https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636), [FinnGen R12 GWAS API](https://r12.finngen.fi/api/pheno/GRAVES_OPHT), [R12 표현형 정의](https://r12.risteys.finngen.fi/endpoints/GRAVES_OPHT).

## 2. 실제 파일의 항목

공통 핵심은 변이 ID, 염색체·위치, 효과/비교 대립유전자, β, SE, P값, 대립유전자 빈도다. FinnGen은 환자군·대조군 빈도와 nearest_genes를 추가 제공한다. nearest_genes는 위치상 가까운 유전자이고 원인 유전자 판정이 아니다. UKB 원본에는 INFO와 여러 BOLT-LMM 검정 열이 있으며, 사용 가공본에는 hm_ 접두사의 조화 처리 열이 있다. 원본 β와 hm_ 대립유전자를 무작정 섞으면 방향 오류가 생길 수 있다. UKB의 OR/CI 열은 비어 있으며 선형척도 β를 그대로 OR로 읽어서는 안 된다.

현재 분석의 UKB 변환은 β와 SE를 μ(1−μ)로 나누는 근사이며 μ=3,731/484,598이다. 실제 로지스틱 회귀를 새로 적합한 값과 동일하다고 단정하지 않는다. eQTLGen은 주로 혈액 발현의 Z 통계와 행별 N을 제공한다. 기존 연구는 이를 빈도로 보완해 β·SE를 재구성했다.

예를 들어 TSHR의 rs179252는 BBJ·FinnGen에서는 G 대립유전자에 음의 β, UKB에서는 T 대립유전자에 양의 β로 기록된다. 이는 방향 불일치의 증거가 아니다. 대립유전자를 맞추면 같은 방향이다. 원본 GRCh37 위치 81,435,985와 사용 가공본 GRCh38 위치 80,969,641도 구분해야 한다. 현재 세 질환 파일의 이 변이 위치는 모두 GRCh38과 일치하며, BBJ·UKB 가공본의 MD5는 공식 파일 메타데이터와 정확히 일치했다.

전체 열의 이름·한국어 설명·원본 예시는 `Raw_file_column_dictionary.csv` 및 대시보드의 ‘원본 항목’ 탭에 보존했다.

## 3. 이번 전체 파일 점검 결과

| 자료 | 실제 변이 행 수 | rsID 없는 행 | 핵심 β·SE·P/대립유전자 이상 행 | P < 5×10⁻⁸ 변이 행 |
|---|---:|---:|---:|---:|
'''
for r in A['raw_qc']:report+=f"| {r['cohort']} | {r['rows']:,} | {r['missing_rsid']:,} | {r['rows']-r['valid_core_rows']:,} | {r['p_lt_5e_8']:,} |\n"
report+='''
위 유의 변이 개수는 상관된 변이를 그대로 센 값으로 독립 유전자나 독립 유전자좌의 개수가 아니다. 서로 다른 조상집단·질환·표본 수·변이 집합이므로 이 개수로 어느 질환의 유전성이 더 강하다고 순위를 매길 수 없다. 핵심 이상 행 0은 P 범위, 유한 β, 양의 SE, 대립유전자 누락을 점검한 결과이며, 모든 유전학적 QC를 통과했다는 보증이 아니다. 이번에는 중복 변이 전수 확인·INFO 재필터·LD 독립성 검정을 수행하지 않았다. rsID 없는 행은 좌표 변이가 존재할 수 있지만 현재 rsID 기반 연결에서 손실될 수 있다.

공식 카탈로그에 기재된 SNP 수와 로컬 가공 파일 행 수는 처리 단계가 달라 서로 다를 수 있다. 분석 분모에는 실제 사용 파일의 집계를 명시해야 한다. 원본은 읽기만 했으며 새 보고서 폴더에만 산출물을 작성했다.

## 4. 기존 유전자 분석의 전체 양상 재집계

| 자료 | 적격 유전자 | MR 추정 가능 | 추정 불가 | P < 0.05 | 동일한 P < 0.05/2,544 기준 |
|---|---:|---:|---:|---:|---:|
'''
for r in A['mr_summary']:report+=f"| {labels[r['outcome']]} | {r['eligible_genes']:,} | {r['estimable_genes']:,} | {r['not_estimable_genes']:,} | {r['p_lt_0_05']} | {r['p_lt_0_05_over_2544_descriptive']} |\n"
report+='''
마지막 열의 UKB·FinnGen 값은 비교를 위해 이번에 재집계한 기술통계이다. 기존 연구의 발견 기준은 BBJ에 적용되었으며, 이를 사전에 정한 새로운 발견 결과로 바꾸어 보고하지 않는다. 추정 불가는 음성 결과가 아니다. OR·P값만으로 서로 다른 SNP 조합과 질환을 동일한 효과로 묶을 수 없다.

BBJ 발견 13개 중 UKB와 방향이 같은 유전자는 10/13개, FinnGen과 같은 유전자는 추정 가능한 12개 중 9개였다. 이는 기술적인 부호 비교이고 재현 성공률이나 군간 차이의 검정이 아니다. 전체 공통 추정 유전자에서는 방향 일치가 약 절반이어서 전체 2,544개의 결과가 일관되게 재현됐다는 주장은 부적절하다. 대다수의 작은 추정치는 잡음도 포함하므로 전체 부호 불일치를 생물학적 반대 효과로 해석해서도 안 된다.

해석해야 할 패턴은 여러 가지다. TSHR은 세 자료의 MR 방향이 같고 기본 prior에서 BBJ·FinnGen의 H4가 높지만 UKB H4는 낮고 BBJ prior 민감성이 있다. CTLA4는 유럽 자료의 H4는 높고 BBJ H4는 낮다. TNFSF14와 IFNGR1은 BBJ H4가 높지만 FinnGen에서는 낮다. MHC 5개는 강한 BBJ MR P값에 비해 공유변이 지지가 낮다. IGF1R은 BBJ·UKB의 다중검정 전 MR 연관성은 있지만 BBJ 발견 기준과 세 자료의 H4 기준을 충족하지 못한다. 이런 차이는 조직·조상집단·표현형·정밀도·도구 SNP 구성의 차이를 모두 고려해야 하며, 이번 재집계만으로 원인을 특정할 수 없다.

기존 14개 유전자만 공위치화를 했다는 선택 범위도 중요하다. 다른 2,530개 전체에 ‘공위치화 지지 없음’을 부여하면 안 된다. 대시보드에서는 분석하지 않은 항목을 미평가로 표시했다.

## 5. 혈액 발현 자료의 의미

'''+f"eQTLGen 도구 자료는 {A['instrument_records']:,}개의 SNP–유전자 기록, {A['instrument_unique_snps']:,}개 고유 SNP, {A['instrument_genes']:,}개 유전자다. 도구별 표본 수는 중앙값 {A['instrument_n_median']:,.0f}, 범위 {A['instrument_n_min']:,}–{A['instrument_n_max']:,}이다. 모든 도구가 31,684명에서 측정되었다고 쓰면 안 된다.\n\n"+'''
혈액 eQTL은 질환 GWAS에서 환자의 혈액을 직접 측정한 자료가 아니다. 별도 연구에서 SNP–발현 관계를 추정한 뒤 SNP–질환 관계와 결합한다. 따라서 IGF1R의 현재 결과는 ‘이번 혈액 eQTL 기반 분석으로 질환 감수성 관련성을 충분히 확립하지 못했다’는 뜻이다. 혈액에서의 검출 불가능, 안와에서만의 작용, 치료 효과의 기전을 입증하지 않는다. 소규모 안와 RNA-seq는 별도 기술통계 자료이며 이 간극을 단독으로 해결할 수 없다.

## 6. 지금 보유한 추가 자료와 다음 분석

| 추가 FinnGen R12 파일 | GWAS 환자 | 대조군 | 다음 분석에서의 가치 |
|---|---:|---:|---|
| E4_GRAVES_STRICT | 3,962 | 496,386 | BBJ GD와 더 가까운 질환 정의로 확인 |
| E4_GRAVES_OPHT_STRICT | 753 | 499,595 | 안병증 정의를 좁혔을 때 결과가 유지되는지 평가 |

두 파일 모두 로컬에 있고 공식 API의 표본 수·정의를 확인했다. 이번 단계에서는 헤더만 읽었으며 전체 파일 QC·SNP 수준 분석을 수행하지 않았다. strict는 데이터베이스 정의명이며 임상적 진단 정확성을 보증하는 형용사로 사용하지 않는다. GD strict는 정의에 안병증 endpoint도 포함하고, 안병증 strict는 진단 코드 조합으로 정의된다. 두 자료와 현재 안병증 자료는 같은 FinnGen R12에서 나오므로 표본 및 대조군 중복을 고려해야 한다. GD 파일과 안병증 파일의 β를 단순히 빼거나 OR를 비교한다고 GD 내 TED 분석이 되는 것은 아니다.

출처: [GD strict API](https://r12.finngen.fi/api/pheno/E4_GRAVES_STRICT), [GD strict 정의](https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_STRICT), [안병증 strict API](https://r12.finngen.fi/api/pheno/E4_GRAVES_OPHT_STRICT), [안병증 strict 정의](https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_OPHT_STRICT).

## 7. 논문 방향에 대한 잠정 권고

현재 자료에 가장 잘 맞는 중심 질문은 ‘혈액 발현과 GD 감수성을 연결하는 약물표적 후보는 무엇이고, 이 근거는 질환 정의와 인구집단이 달라져도 유지되는가?’이다. TSHR와 IGF1R는 임상적으로 중요한 비교 대상으로 두되, 전체 약물표적 유전자 스크리닝과 일관성 평가가 연구의 중심이 되어야 한다. ‘IGF1R의 effector biology를 입증했다’거나 ‘TED 특이적 발병 기전을 밝혔다’는 결론은 현재 자료 범위를 넘는다.

권하는 순서는 BBJ GD 발견 → FinnGen GD strict에서 동일 후보의 확인 → UKB의 더 넓은 갑상선기능항진증 및 FinnGen 안병증의 확장 비교이다. 이는 기존 자료를 본 뒤 수정하는 분석 계획이므로 후향적 확장임을 명시해야 한다. 새 결과를 보고 기준을 바꾸지 않도록, 후보군·주분석·검정 기준·민감도 분석·결측 처리 방식을 먼저 기록해야 한다.

다음 단계에서는 (1) 두 추가 FinnGen 파일의 QC, (2) 같은 도구 선정 원칙으로 MR을 실행하고 SNP별 제외 이유 확인, (3) 지정한 후보의 전체 영역 공위치화 및 prior/빈도 민감도, (4) 효과크기·CI·정밀도와 자료 간 차이를 함께 평가하는 순서가 적절하다. 조상집단과 표현형이 다른 자료를 곧바로 하나의 pooled OR로 합치거나, 유의/비유의만으로 효과 차이를 주장하지 않는다. TSHR가 최종 우선 후보로 남는지도 이 과정 후 판단해야 한다.

가능한 임시 제목은 “Genetically proxied blood expression of druggable genes in Graves disease: evidence across disease definitions and populations”이다. 추가 분석 전 제목과 결론을 확정하지 않는다. TED 환자 내 위험이나 치료반응이 주목적이라면 그 질문에 맞는 별도의 비교 자료가 필요하다.

## 8. 검증 범위와 남은 점검

완료: 세 사용 GWAS 전체 행 읽기, 핵심 열 범위·누락 점검, SHA-256 기록, BBJ·UKB 공식 MD5 대조, FinnGen R12 API 표본 수 대조, 기존 MR 결과 재집계 및 β/SE–P–OR–CI 산술 검증.

남음: UKB 원논문의 세부 phenotype·조상집단·공변량 선정 기준, 두 추가 파일 전체 QC, 모든 SNP의 새 harmonization 및 탈락 원인 확인, eQTL–GWAS 표본 중복 평가, 변경 분석 계획에 따른 새 MR/공위치화. 이번 점검이 기존 연구 전체의 재현을 대신하지 않는다.

이번 작업은 `outputs/TED_Data_Review_20261002`에 분리 보관했으며 기존 제출용 원고·발표·Table·Figure는 수정하지 않았다.
'''
(O/'REPORT_KO.md').write_text(report,encoding='utf8')
(O/'review_data.json').write_text(json.dumps(payload,ensure_ascii=False,allow_nan=False),encoding='utf8')
template=(Path(__file__).parent/'review_template.html').read_text(encoding='utf8')
(O/'START_HERE.html').write_text(template.replace('__DATA__',json.dumps(payload,ensure_ascii=False,allow_nan=False).replace('</','<\\/')),encoding='utf8')
print('Review report and interactive data explorer created.')
