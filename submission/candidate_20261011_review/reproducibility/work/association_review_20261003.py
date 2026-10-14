from pathlib import Path
import sys, json, hashlib, math
import pandas as pd
sys.stdout.reconfigure(encoding='utf-8')
R=Path(__file__).resolve().parents[1]
O=R/'outputs/TED_Association_Review_20261003';O.mkdir(exist_ok=True)
paths={
 'MR':R/'outputs/TED_Reanalysis_20261002/Tables/All_2544_genes_five_outcomes.csv',
 'selected':R/'outputs/TED_Followup_Validation_20261002/Tables/Table_F1_24_genes_five_outcomes.csv',
 'frequency':R/'outputs/TED_Reanalysis_20261002/Tables/New_MR_all_frequency_scenarios.csv',
 'workbook':R/'outputs/TED_Manuscript_Reframed_20261002/Supplementary_Dataset.xlsx'}
d=pd.read_csv(paths['MR']);e=pd.read_csv(paths['selected']);f=pd.read_csv(paths['frequency'])
order=['BBJ_Graves','E4_GRAVES_STRICT','UKB_hyperthyroid','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
labels=dict(zip(order,['BBJ GD','FinnGen GD','UKB hyperthyroidism','FinnGen ophthalmopathy','FinnGen strict ophthalmopathy']))
threshold=.05/2544;focused=.05/13
assert len(d)==2544*5 and len(e)==24*5
merged=e.merge(d,on=['gene_symbol','outcome'],suffixes=('_selected','_screen'),validate='one_to_one')
for c in ['OR','CI_lower','CI_upper','pvalue','beta','se']:
 for a,b in zip(merged[c+'_selected'],merged[c+'_screen']):
  assert (pd.isna(a) and pd.isna(b)) or math.isclose(a,b,rel_tol=1e-12,abs_tol=1e-15),(c,a,b)
bbj=d[(d.outcome==order[0]) & (d.pvalue<threshold)].sort_values('pvalue')
fg=d[(d.outcome==order[1]) & (d.pvalue<threshold)].sort_values('pvalue')
paired=bbj[['gene_symbol','beta','pvalue']].merge(d[d.outcome==order[1]][['gene_symbol','beta','pvalue']],on='gene_symbol',suffixes=('_BBJ','_FG'))
paired['same_direction']=paired.beta_BBJ*paired.beta_FG>0
paired['both_full_threshold']=paired.same_direction & (paired.pvalue_FG<threshold)
paired['focused_corroboration']=paired.same_direction & (paired.pvalue_FG<focused)
dual=set(paired.loc[paired.both_full_threshold,'gene_symbol']);focus=set(paired.loc[paired.focused_corroboration,'gene_symbol'])
assert dual=={'TSHR','CTLA4','HLA-DQA2','TUBB'}
assert focus=={'TSHR','CTLA4','HLA-DQA2','TUBB','C4A','TNFSF14','IFNGR1'}
loo=pd.read_excel(paths['workbook'],sheet_name='LOO_original')
igloo=loo[(loo.gene=='IGF1R') & (loo.excluded_SNP=='rs2654980')]
assert set(igloo.outcome)=={'BBJ_Graves','UKB_hyperthyroid','FinnGen_GO'}
ig=e[e.gene=='IGF1R'].set_index('outcome').loc[order].reset_index()
assert (ig.OR>1).all() and (ig.pvalue<.05).sum()==2 and (ig.pvalue<threshold).sum()==0
notes={
 'CTLA4':'Same-direction nominal associations in all five datasets; multi-IV results depend strongly on rs13030124. Colocalization varies across datasets.',
 'HLA-DQA2':'Both GD datasets and both ophthalmopathy definitions support MR association; UKB does not. Single IV for GD/GO; low PP.H4 throughout.',
 'TUBB':'Same-direction nominal associations in all five datasets. Heterogeneity in FinnGen GD (Q P=3.67e-5) and UKB; low PP.H4 throughout.',
 'C4A':'Supports focused GD corroboration, not the full FinnGen screen threshold. Both GO definitions significant; strong FinnGen GD instrument heterogeneity.',
 'TNFSF14':'Supports focused GD corroboration; nominal UKB association; neither GO definition significant.',
 'IFNGR1':'Supports focused GD corroboration; UKB and both GO definitions not nominally significant.'}
ig.to_csv(O/'IGF1R_five_datasets.csv',index=False,encoding='utf-8-sig')
paired.to_csv(O/'BBJ_discovery_GD_corroboration.csv',index=False,encoding='utf-8-sig')
e.to_csv(O/'Selected_24_genes_five_datasets.csv',index=False,encoding='utf-8-sig')
extra=e[e.gene.isin(set(fg.gene_symbol)-set(bbj.gene_symbol))]
assert len(extra)==45
extra.to_csv(O/'Additional_FinnGen_GD_hits_five_datasets.csv',index=False,encoding='utf-8-sig')
igloo.to_csv(O/'IGF1R_rs2654980_omission.csv',index=False,encoding='utf-8-sig')
diag=f[(f.gene_symbol.isin(list(notes)+['IGF1R','APOM','NEU1','C4B'])) & (f.scenario.isin(['original_reference','reharmonized_eqtlgen']))]
diag.to_csv(O/'Frequency_and_heterogeneity_diagnostics.csv',index=False,encoding='utf-8-sig')
def pval(x):
 if pd.isna(x):return 'NA'
 if x<.001:return f'{x:.2e}'
 return f'{x:.4g}'
def effect(q):return f'{q.OR:.2f} ({q.CI_lower:.2f}–{q.CI_upper:.2f}); P={pval(q.pvalue)}'
def row(g,o):return e[(e.gene==g)&(e.outcome==o)].iloc[0]
lines=['# 혈액 발현 도구변수와 질환의 연관성: 검토 요약',
 '작성일: 2026-10-03. 검증된 기존 분석 산출물을 다시 분류·대조했습니다. 새 환자 분석이나 새 유의성 검정을 시행하지 않았고, 현재 제출 검토 패키지는 변경하지 않았습니다.',
 '## 보고 방향',
 'MR 연관성의 방향·효과 추정치·정밀도와 질환 자료 간 일관성을 주요 결과로 보고하는 것이 적절합니다. 실제 측정한 혈중 발현과 질환의 관찰연구가 아니라 유전적 도구변수를 이용한 MR 결과입니다. 공위치화는 발현 매개 해석을 평가하는 핵심 결과로 함께 유지합니다. P값이 작다는 이유만으로 해당 유전자의 독립적인 인과 역할이나 치료 효과를 주장하지 않습니다.',
 'OR는 재구성한 혈액 발현 약 1 SD 증가에 대한 추정치입니다. OR 크기로 서로 다른 유전자의 생물학적 중요성을 순위화하지 않았습니다.',
 '## IGF1R',
 '| 질환 자료 | OR (95% CI); P | PP.H4 |', '|---|---|---|']
for _,q in ig.iterrows():lines.append(f'| {labels[q.outcome]} | {effect(q)} | {q["PP.H4"]:.3f} |')
lines += ['모든 점 추정치는 odds 증가 방향이지만 BBJ와 UKB만 보정 전 P<0.05입니다. 모든 결과가 원래의 전체 스크린 기준 P<0.05/2,544를 넘지 못합니다. 이를 강건하게 재현된 연관성 또는 효과 부재로 표현하지 않습니다. FinnGen GD OR 1.07의 95% CI는 0.86–1.33입니다. 세 FinnGen endpoint는 중복 참여자가 있어 세 번의 독립 재현이 아닙니다.',
 'rs2654980 제외 후: '+ '; '.join(f'{labels[q.outcome]} P={pval(q.pvalue)}' for _,q in igloo.iterrows())+'. 특정 SNP 제외 후 유의성 소실은 정밀도 감소로도 나타날 수 있으므로, 그 SNP가 반드시 편향을 일으켰거나 모든 효과를 설명한다고 단정하지 않습니다.',
 '권장 영문: IGF1R estimates were directionally consistent with higher disease odds, but nominal associations in BBJ GD and UKB hyperthyroidism did not meet the screen-wide significance threshold and were not corroborated at nominal significance in the FinnGen endpoints.',
 '## TSHR 외 BBJ 발견 유전자의 GD 비교',
 f'전체 스크린 기준은 P<{threshold:.8g}; BBJ 13개 발견에 대한 FinnGen GD 확인 기준은 P<{focused:.8g}입니다. 이 분류는 유전자 인과성 순위나 새로운 복합 점수가 아닙니다.',
 '| 유전자 | BBJ GD: OR (95% CI); P | FinnGen GD: OR (95% CI); P | 근거 범위 |', '|---|---|---|---|']
for g in notes:
 tier='두 GD 모두 전체 스크린 보정 기준 통과' if g in dual else 'BBJ 발견 + FinnGen GD 13개 확인 보정 기준 통과'
 lines.append(f'| {g} | {effect(row(g,order[0]))} | {effect(row(g,order[1]))} | {tier} |')
lines+=['### 해석에 필요한 단서']+[f'- {g}: {note}' for g,note in notes.items()]
lines+=['CTLA4와 TUBB는 다섯 자료 모두 같은 방향의 보정 전 유의성을 보입니다. HLA-DQA2는 UKB의 추정 방향이 반대이고 비유의입니다. 자료 간 유의/비유의 차이 자체를 효과 차이 검정으로 해석하지 않습니다.',
 '## 그 밖에 P값이 작은 결과',
 '- APOM과 NEU1: FinnGen GD P=1.43e-108, UKB P=3.62e-111. 같은 도구변수 rs3101018을 사용하므로 별개의 두 독립 발견으로 세면 안 됩니다. BBJ MR은 추정 불가이며 두 유전자 모두 강한 공위치화 지지가 없습니다.',
 '- TNF, HIST1H4E, LTA, LTB, IGLV2-18, IGLV3-12, C4B도 추가 FinnGen GD 스크린 결과입니다. IGLV 두 유전자는 같은 SNP를 공유합니다. LTA는 BBJ와 FinnGen GD의 점 추정 방향이 다릅니다. 이는 공식적인 질환/집단 차이 검정 결과가 아닙니다.',
 '- C4B는 다섯 자료에서 보정 전 유의하며 GO 연관성이 강하지만, eQTLGen 빈도 대체 후 FinnGen GD P가 1.10e-5에서 1.48e-4로 바뀌어 전체 스크린 기준을 통과하지 못합니다. 이질성도 동반됩니다.',
 '- HLA-A는 BBJ P=2.60e-23으로 가장 강하지만 FinnGen GD P=0.648이며 다른 자료에서도 보정 전 유의하지 않습니다. PSMB8, MAPKAPK5 역시 BBJ 발견을 다른 자료가 같은 수준으로 뒷받침하지 않습니다. VKORC1의 FinnGen MR은 추정 불가입니다.',
 '## 본문 배치 제안',
 'TSHR의 일관된 MR 연관성 → IGF1R의 증가 방향이나 제한적인 근거 → CTLA4·HLA-DQA2·TUBB 등 비교 유전자 → 공위치화·민감도 결과 순서가 적절합니다. 중심 질문은 TSHR/IGF1R에 두고 나머지 결과는 선택성 없이 전체 표에 보고합니다.',
 '권장 영문 핵심 문장: TSHR showed consistent MR associations with lower disease odds across the five datasets. IGF1R estimates were consistently above the null, but statistical support was limited to nominal associations in BBJ GD and UKB hyperthyroidism. CTLA4, HLA-DQA2 and TUBB also met the screen-wide threshold in both GD datasets. Colocalization and instrument diagnostics limited expression-mediated causal interpretation.',
 '## 자료와 해석 근거',
 '수치는 동봉 CSV와 검증된 2026-10-02 전체 분석 산출물에서 가져왔습니다. 24개 유전자·5개 질환 자료의 MR 수치 120행을 전체 2,544개 스크린 결과와 대조했습니다.',
 '[STROBE-MR 설명 및 해설](https://www.bmj.com/content/375/bmj.n2233): MR의 인과 해석에는 도구변수 가정 및 유전적 차이와 개입의 동등성에 대한 검토가 필요합니다. 이 문헌은 우리 연구의 수치 출처가 아니라 해석·보고 원칙의 근거입니다.']
(O/'Association_Review_KO.md').write_text('\n\n'.join(lines)+'\n',encoding='utf-8')
audit={'status':'PASS','existing_analysis_summary_only':True,'manuscript_package_modified':False,'screen_threshold':threshold,'focused_threshold':focused,'BBJ_discoveries':len(bbj),'FinnGen_GD_discoveries':len(fg),'both_GD_screen_significant':sorted(dual),'focused_GD_supported':sorted(focus),'selected_rows_crosschecked':len(merged),'IGF1R_nominal_count':2,'IGF1R_screen_significant_count':0,'source_sha256':{k:hashlib.sha256(p.read_bytes()).hexdigest() for k,p in paths.items()}}
(O/'verification.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(audit,ensure_ascii=False,indent=2))
