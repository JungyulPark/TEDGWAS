from coloc_core_v1 import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import json,shutil,hashlib,datetime
O=W.parent/'outputs/TED_Followup_Validation_20261002';TABLES=O/'Tables';FIG=O/'Figures';VERIFY=O/'Verification'
LABEL={'BBJ':'BBJ GD','UKB':'UKB hyperthyroidism','FinnGen_GO':'FinnGen GO','E4_GRAVES_STRICT':'FinnGen GD strict','E4_GRAVES_OPHT_STRICT':'FinnGen GO strict'}
ORDER=['BBJ','E4_GRAVES_STRICT','UKB','FinnGen_GO','E4_GRAVES_OPHT_STRICT']
COLORS=['#246991','#b65322','#606b7d','#497b45','#735a93']
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'text.color':'black','axes.labelcolor':'black','xtick.color':'black','ytick.color':'black','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.dpi':230})
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(fig,name):
    fig.savefig(FIG/f'{name}.png',bbox_inches='tight',facecolor='white');fig.savefig(FIG/f'{name}.pdf',bbox_inches='tight',facecolor='white');plt.close(fig)
def mdtable(df):
    return '| '+' | '.join(df.columns)+' |\n|'+'|'.join(['---']*len(df.columns))+'|\n'+'\n'.join('| '+' | '.join(map(str,row))+' |' for row in df.itertuples(index=False,name=None))
def pformat(p):
    if pd.isna(p):return '추정 불가'
    return f'{p:.3g}' if p>=.001 else f'{p:.2e}'
def ppformat(p):return f'{p:.3f}' if p>=.001 else f'{p:.2e}'
def main():
    assert json.loads((T/'final_checks_v1/verification.json').read_text())['PASS']
    alias=json.loads((T/'alias_audit_v1/status.json').read_text());assert alias['complete'] and alias['alias_hits']==0
    O.mkdir();TABLES.mkdir();FIG.mkdir();VERIFY.mkdir();(VERIFY/'Scripts').mkdir();(VERIFY/'Sources').mkdir()
    tshr=pd.read_csv(T/'TSHR_v2/TSHR_240_coloc_settings.csv');cand=pd.read_csv(T/'candidates_v1/Candidate_600_coloc_settings.csv')
    coverage=pd.read_csv(T/'final_checks_v1/marginal_signal_coverage.csv')
    # Complete combined evidence table: existing fixed 14 plus ten newly examined candidates.
    oldcol=pd.read_csv(W.parent/'outputs/TED_Reanalysis_20261002/Tables/Colocalization_840_settings.csv')
    newcol=cand.drop(columns=['window','p1_effective','p2_effective','p12_effective']).copy();newcol['outcome']=newcol.outcome.replace(ALIASES)
    allcol=pd.concat([oldcol,newcol],ignore_index=True);assert len(allcol)==1440 and not allcol.duplicated(['gene','outcome','scenario','p12']).any()
    allcol.to_csv(TABLES/'All_24_genes_1440_coloc_settings.csv',index=False)
    mr=pd.read_csv(W.parent/'outputs/TED_Reanalysis_20261002/Tables/All_2544_genes_five_outcomes.csv')
    primary=allcol[allcol.scenario.eq('original_reference')&allcol.p12.eq(1e-5)].copy()
    merged=primary.merge(mr,left_on=['gene','outcome'],right_on=['gene_symbol','outcome'],validate='one_to_one',how='left')
    assert len(merged)==120 and merged.gene_symbol.notna().all()
    lows=allcol[allcol.scenario.eq('original_reference')&allcol.p12.eq(1e-6)][['gene','outcome','PP.H4']].rename(columns={'PP.H4':'H4_low_prior'})
    newaf=allcol[allcol.scenario.eq('available_eqtlgen')&allcol.p12.eq(1e-5)][['gene','outcome','PP.H4']].rename(columns={'PP.H4':'H4_eqtlgen'})
    merged=merged.merge(lows,on=['gene','outcome']).merge(newaf,on=['gene','outcome'])
    merged['analysis_group']=np.where(merged.gene.isin(cand.gene),np.where(merged.gene.eq('HLA-G'),'frequency_sensitivity_only','new_exploratory_candidate'),'previous_fixed_14')
    merged['interpretation_note']='Single-causal-variant model; no proof of expression-mediated causality'
    weak=merged.gene.str.startswith('IGLV')&merged.outcome.eq('BBJ_Graves');merged.loc[weak,'interpretation_note']='Lead eQTL absent; weak expression and disease signals in overlap; uninformative comparison'
    merged.to_csv(TABLES/'Table_F1_24_genes_five_outcomes.csv',index=False)
    mappings={
      T/'TSHR_v2/TSHR_240_coloc_settings.csv':'Table_F2_TSHR_windows_priors_frequencies.csv',
      T/'TSHR_v2/regional_coverage_counts.csv':'TSHR_raw_region_coverage_counts.csv',
      T/'TSHR_v2/raw_GWAS_top10_by_scope.csv':'TSHR_raw_region_lead_variants.csv',
      T/'TSHR_v2/raw_GWAS_region_audit.csv':'TSHR_all_raw_regional_variants.csv',
      T/'candidates_v1/Candidate_primary_MR_coloc.csv':'Table_F3_new_candidates_MR_coloc.csv',
      T/'candidates_v1/Candidate_600_coloc_settings.csv':'New_candidates_600_coloc_settings.csv',
      T/'candidates_v1/coverage.csv':'New_candidates_variant_coverage.csv',
      T/'candidates_v1/shared_selected_instruments.csv':'New_candidates_shared_instruments.csv',
      T/'instruments_v1/full_MR.csv':'Table_F4_C4B_HLAG_frequency_MR.csv',
      T/'instruments_v1/leave_one_out.csv':'C4B_HLAG_leave_one_out_72.csv',
      T/'instruments_v1/SNP_estimates_weights.csv':'C4B_HLAG_SNP_weights_estimates.csv',
      T/'instruments_v1/exposure_frequency_scaling.csv':'C4B_HLAG_frequency_scaling.csv',
      T/'final_checks_v1/marginal_signal_coverage.csv':'All_11_genes_marginal_signal_coverage.csv',
      T/'inputs/eqtl_coverage.csv':'Source_eQTL_regions_sample_sizes.csv'}
    for src,name in mappings.items():shutil.copy2(src,TABLES/name)
    # Figure F1: predefined windows and prior dependence, same y scale.
    fig,(ax,bx)=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'width_ratios':[1.1,1]})
    windows=['full_source_cis','GenePos_pm500kb','GenePos_pm250kb','five_outcome_intersection']
    xlabels=['Full source\ncis region','Gene position\n±500 kb','Gene position\n±250 kb','Common SNPs\nacross 5 outcomes']
    for oc,color in zip(ORDER,COLORS):
        q=tshr[(tshr.outcome==oc)&tshr.scenario.eq('original_reference')&tshr.p12.eq(1e-5)].set_index('window').loc[windows]
        ax.plot(range(4),q['PP.H4'],marker='o',ms=5,lw=1.7,color=color,label=LABEL[oc])
        for i,val in enumerate(q['PP.H4']):
            off=(-15 if oc in ['BBJ','UKB','E4_GRAVES_OPHT_STRICT'] else 8)
            ax.annotate(f'{val:.3f}',(i,val),xytext=(0,off),textcoords='offset points',ha='center',fontsize=8.5,color='black')
        q=tshr[(tshr.outcome==oc)&tshr.scenario.eq('original_reference')&tshr.window.eq('full_source_cis')].sort_values('p12')
        bx.plot(q.p12,q['PP.H4'],marker='o',ms=5,lw=1.7,color=color,label=LABEL[oc])
    for a in [ax,bx]:
        a.axhline(.8,color='black',ls='--',lw=1);a.set_ylim(0,1.065);a.set_ylabel('PP.H4: probability of a shared variant');a.grid(axis='y',alpha=.18)
    ax.set_xticks(range(4),xlabels,fontsize=9);ax.set_title('A   TSHR results across predefined windows',loc='left',fontweight='bold',pad=12)
    bx.set_xscale('log');bx.set_xticks([1e-6,5e-6,1e-5],['10⁻⁶','5 × 10⁻⁶','10⁻⁵']);bx.minorticks_off();bx.set_xlabel('Prior probability per SNP of a shared association (p12)')
    bx.set_title('B   TSHR sensitivity to the shared-variant prior',loc='left',fontweight='bold',pad=12)
    handles,labels=ax.get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,frameon=False,bbox_to_anchor=(.5,.075),labelcolor='black')
    fig.suptitle('Follow-up Figure 1. TSHR colocalization remains dependent on the disease dataset',fontsize=15,fontweight='bold',y=.99)
    fig.subplots_adjust(bottom=.30,top=.86,wspace=.27)
    fig.text(.08,.025,'Reference allele frequencies. Dashed line: PP.H4 = 0.80; this is not a P-value threshold.\nOverlapping FinnGen endpoints are not independent replications.',fontsize=10)
    save(fig,'Figure_F1_TSHR_windows_and_priors')
    # Figure F2: all candidate posteriors, with an explicit coverage flag.
    primary=cand[cand.scenario.eq('original_reference')&cand.p12.eq(1e-5)]
    genes=['APOM','NEU1','C4B','HIST1H4E','LTA','LTB','TNF','HLA-G','IGLV2-18','IGLV3-12']
    matrix=primary.pivot(index='gene',columns='outcome',values='PP.H4').loc[genes,ORDER]
    fig,ax=plt.subplots(figsize=(12,8.3));im=ax.imshow(matrix,cmap='Blues',vmin=0,vmax=1,aspect='auto')
    for i,g in enumerate(genes):
        for j,oc in enumerate(ORDER):
            val=matrix.loc[g,oc];label=f'{val:.3f}' if val>=.001 else '<0.001'
            if g.startswith('IGLV') and oc=='BBJ':
                ax.add_patch(Rectangle((j-.5,i-.5),1,1,facecolor='#efefef',edgecolor='#cccccc',hatch='///',lw=.5));label='Limited†'
            ax.text(j,i,label,ha='center',va='center',color='black',fontsize=11)
    ax.set_xticks(range(5),[LABEL[x].replace('FinnGen ','FinnGen\n').replace('UKB ','UKB\n') for x in ORDER],fontsize=10);ax.set_yticks(range(10),[g+(' *' if g=='HLA-G' else '') for g in genes],fontstyle='italic')
    ax.set_xticks(np.arange(-.5,5,1),minor=True);ax.set_yticks(np.arange(-.5,10,1),minor=True);ax.grid(which='minor',color='white',linewidth=2);ax.tick_params(which='minor',bottom=False,left=False)
    for s in ax.spines.values():s.set_visible(False)
    cb=fig.colorbar(im,ax=ax,fraction=.035,pad=.04);cb.set_label('PP.H4 (shared-variant probability)');cb.ax.axhline(.8,color='black',ls='--',lw=1)
    ax.set_title('Follow-up Figure 2. New MR candidates did not reach PP.H4 ≥ 0.80',fontweight='bold',fontsize=15,pad=20)
    fig.subplots_adjust(left=.16,right=.93,top=.89,bottom=.19)
    fig.text(.16,.055,'Full source cis regions; reference frequencies; p12 = 10⁻⁵.\n* HLA-G was selected only in the allele-frequency sensitivity screen.\n† In BBJ, the lead eQTL is absent and both marginal signals are weak in the overlap.\nAPOM / NEU1 and IGLV2-18 / IGLV3-12 each share an MR instrument.',fontsize=10)
    save(fig,'Figure_F2_new_candidate_colocalization')
    # Figure F3: HLA-G weight change and the decisive leave-one-out comparison.
    snps=pd.read_csv(T/'instruments_v1/SNP_estimates_weights.csv');full=pd.read_csv(T/'instruments_v1/full_MR.csv');loo=pd.read_csv(T/'instruments_v1/leave_one_out.csv')
    s=snps[snps.gene.eq('HLA-G')&snps.outcome.eq('E4_GRAVES_STRICT')];a=s.pivot(index='SNP',columns='scenario',values='weight_fraction').sort_values('reharmonized_eqtlgen',ascending=False)
    fig,(ax,bx)=plt.subplots(1,2,figsize=(14,8),gridspec_kw={'width_ratios':[1.05,1]})
    y=np.arange(len(a));ax.barh(y-.17,100*a.original_reference,height=.32,color='#b2bdc8',label='1000 Genomes reference');ax.barh(y+.17,100*a.reharmonized_eqtlgen,height=.32,color='#3c718d',label='eQTLGen frequencies')
    ax.set_yticks(y,a.index,fontsize=9);ax.invert_yaxis();ax.set_xlabel('Contribution to the IVW estimate (%)');ax.set_xlim(0,67);ax.legend(loc='lower right',frameon=False,fontsize=9,labelcolor='black');ax.set_title('A   SNP contribution weights',loc='left',fontweight='bold',pad=15)
    a1=full[full.gene.eq('HLA-G')&full.outcome.eq('E4_GRAVES_STRICT')].copy();a1['condition']='All 16 instruments'
    a2=loo[loo.gene.eq('HLA-G')&loo.outcome.eq('E4_GRAVES_STRICT')&loo.removed_snp.eq('rs150254595')].copy();a2['condition']='Omit rs150254595'
    r=pd.concat([a1,a2],ignore_index=True)
    labels=[]
    for j,row in r.iterrows():
        color='#b2bdc8' if row.scenario=='original_reference' else '#3c718d';freq='Reference' if row.scenario=='original_reference' else 'eQTLGen'
        bx.errorbar(row.OR,j,xerr=[[row.OR-row.CI_lower],[row.CI_upper-row.OR]],fmt='o',color=color,capsize=4,ms=7)
        labels.append(row.condition+'\n'+freq+f'; P = {pformat(row.pvalue)}')
    bx.set_yticks(range(4),labels,fontsize=10);bx.invert_yaxis();bx.set_xscale('log');bx.set_xlim(.4,1.4);bx.set_xticks([.5,.7,1,1.3],['0.5','0.7','1.0','1.3']);bx.minorticks_off();bx.axvline(1,color='black',ls='--',lw=1);bx.set_xlabel('MR odds ratio (95% confidence interval)');bx.set_title('B   Removing the dominant SNP',loc='left',fontweight='bold',pad=15)
    fig.suptitle('Follow-up Figure 3. HLA-G sensitivity arises from changes in SNP weighting',fontsize=15,fontweight='bold',y=.98);fig.subplots_adjust(left=.12,right=.97,top=.87,bottom=.16,wspace=.65)
    fig.text(.12,.05,'FinnGen GD strict. Identical retained SNPs, disease estimates and expression Z-scores in both frequency analyses.\nThe same IVW estimator with an underdispersion floor was used throughout. Weights are not causal probabilities.',fontsize=10)
    save(fig,'Figure_F3_HLAG_weights_and_leave_one_out')
    # Figure F4: C4B two-instrument heterogeneity shown directly.
    fig,axes=plt.subplots(1,2,figsize=(13,5.5),sharex=True)
    for ax,sc,title in zip(axes,['original_reference','reharmonized_eqtlgen'],['A   1000 Genomes reference','B   eQTLGen frequencies']):
        one=full[full.gene.eq('C4B')&full.outcome.eq('E4_GRAVES_STRICT')&full.scenario.eq(sc)].iloc[0]
        qs=snps[snps.gene.eq('C4B')&snps.outcome.eq('E4_GRAVES_STRICT')&snps.scenario.eq(sc)].set_index('SNP')
        vals=[('Combined IVW',one.OR,one.CI_lower,one.CI_upper,one.pvalue)]
        for snp in ['rs3131383','rs9264533']:
            r=qs.loc[snp];vals.append((snp+' alone',r.OR_wald,r.CI_lower,r.CI_upper,r.P_wald))
        for j,(label,OR,lo,hi,p) in enumerate(vals):
            ax.errorbar(OR,j,xerr=[[OR-lo],[hi-OR]],fmt='D' if j==0 else 'o',color='#3c718d' if j==0 else '#6d7d89',capsize=4,ms=7)
            ax.text(2.65,j,'P = '+pformat(p),va='center',fontsize=9,color='black',clip_on=False)
        ax.set_yticks(range(3),[r[0] for r in vals]);ax.set_ylim(2.6,-.7);ax.set_xscale('log');ax.set_xlim(.95,2.55);ax.set_xticks([1,1.5,2,2.5],['1.0','1.5','2.0','2.5']);ax.minorticks_off();ax.axvline(1,color='black',ls='--',lw=1)
        ax.set_title(title+f'\nHeterogeneity P = {pformat(one.Q_pvalue)}',loc='left',fontweight='bold',fontsize=11,pad=12);ax.set_xlabel('MR odds ratio (95% confidence interval)')
    fig.suptitle('Follow-up Figure 4. C4B instruments agree in direction but differ in effect size',fontsize=14,fontweight='bold',y=.99)
    fig.subplots_adjust(left=.15,right=.84,wspace=1.25,top=.77,bottom=.22);fig.text(.15,.06,'FinnGen GD strict. Single-SNP P-values are unchanged by frequency rescaling.\nThe full IVW estimate reflects substantial between-instrument heterogeneity.',fontsize=10)
    save(fig,'Figure_F4_C4B_instrument_heterogeneity')
    # Sources and verification material; snapshots are distinguished from current retrieval.
    sourcebase=W.parent/'outputs/TED_Data_Review_20261002/provenance'
    for name in ['UKB_yaml.txt','UKB_study_retry.json','BBJ_study.json','E4_GRAVES_STRICT_api.json','E4_GRAVES_OPHT_STRICT_api.json','GRAVES_OPHT_api.json','FinnGen_R12_strictGD_definition_text.txt','FinnGen_R12_strictGO_definition_text.txt']:
        shutil.copy2(sourcebase/name,VERIFY/'Sources'/name)
    for p in T.glob('*.py'):shutil.copy2(p,VERIFY/'Scripts'/p.name)
    for p in T.glob('*.R'):shutil.copy2(p,VERIFY/'Scripts'/p.name)
    for name in ['PLAN.md','plan_freeze.json','frozen_windows.json','run.log']:shutil.copy2(T/name,VERIFY/name)
    for sub,names in {'inputs':['eqtl_manifest.json','outcomes_manifest.json'], 'TSHR_v2':['R_verification.csv','R_sessionInfo.txt','R_warnings.txt','baseline_identity.csv'],'candidates_v1':['R_verification.csv','R_sessionInfo.txt','R_warnings.txt'],'instruments_v1':['R_verification.csv','R_sessionInfo.txt'],'final_checks_v1':['verification.json','retained_instrument_source_identity.csv'],'alias_audit_v1':['summary.csv','hits.csv','status.json']}.items():
        dest=VERIFY/sub;dest.mkdir()
        for name in names:shutil.copy2(T/sub/name,dest/name)
    (VERIFY/'VERSION_NOTES.md').write_text('TSHR_v2 is canonical. TSHR_v1 matched raw coverage flags by rsID only, overcounting multiallelic raw rows; v2 uses rsID + chromosome + position + alleles. Coloc inputs and all 240 posterior vectors are unchanged. Earlier v1 is retained in the workspace for traceability. R locale startup warnings do not affect these numeric checks. R minimum-P warnings are preserved; weak marginal signals are flagged in Tables and the report.\n\nScripts reproduce within the original workspace layout and use read-only local source paths. The delivery is not a standalone raw-data archive.\n',encoding='utf8')
    # Compact tables for the readable Korean report.
    tt=tshr[tshr.scenario.eq('original_reference')&tshr.p12.eq(1e-5)].pivot(index='outcome',columns='window',values='PP.H4').loc[ORDER,windows]
    tt.index=[LABEL[x] for x in ORDER];tt=tt.map(lambda x:f'{x:.3f}');tt.columns=['전체 cis','±500 kb','±250 kb','5자료 공통 변이'];tt=tt.reset_index(names='자료')
    gd=merged[merged.outcome.eq('E4_GRAVES_STRICT')&merged.gene.isin(cand.gene)].copy().set_index('gene').loc[genes].reset_index()
    gt=pd.DataFrame({'유전자':gd.gene,'MR OR (95% CI)':[f'{r.OR:.3f} ({r.CI_lower:.3f}–{r.CI_upper:.3f})' for r in gd.itertuples()],'MR P':[pformat(x) for x in gd.pvalue],'PP.H3':[f'{x:.3f}' for x in gd['PP.H3']],'PP.H4':[ppformat(x) for x in gd['PP.H4']]})
    fgd=full[full.outcome.eq('E4_GRAVES_STRICT')].copy();ft=pd.DataFrame({'유전자':fgd.gene,'빈도':fgd.scenario.map({'original_reference':'1000G EUR','reharmonized_eqtlgen':'eQTLGen'}),'OR':[f'{x:.3f}' for x in fgd.OR],'MR P':[pformat(x) for x in fgd.pvalue],'이질성 P':[pformat(x) for x in fgd.Q_pvalue]})
    assert merged[merged.outcome.eq('E4_GRAVES_STRICT')]['PP.H4'].lt(.8).all()
    report=f'''# 추가 검증 결과 — TSHR부터 신규 후보까지

2026-10-02. 승인된 네 단계를 순서대로 수행했다. 기존 원고와 제출용 Table/Figure는 수정하지 않았고, 이번 결과는 연구 방향을 다시 정하기 위한 추가 검증본이다. 분석 선택은 전체 연구 관점에서 사후 분석이며, 이번 계산 전에 범위·후보·기준을 [계획서](Verification/PLAN.md)에 기록했다.

**핵심 결론: TSHR의 강한 MR 연관성은 유지되지만, 혈액 발현과 질환이 같은 변이를 공유한다는 근거는 자료별로 다르다. 이 차이는 이번에 검사한 변이 누락·영역 설정으로 설명되지 않았다. 신규 후보도 공유변이 기준을 통과하지 않았다.**

## 1. TSHR — 무엇을 새로 확인했나

- 전체 eQTL 127,341,798행에서 TSHR 7,309행을 재추출했고 기존 입력의 변이·대립유전자·Z·표본 수가 일치했다. 유의한 변이만 선택한 자료가 아니다.
- 각 질환 GWAS의 chr14:80,000,000–82,100,000(GRCh38)을 rsID 매칭과 무관하게 읽었다. 다섯 자료 모두 가장 강한 질환 연관 변이가 원래 공위치화 입력에 포함돼 있었다.
- 원본 cis 영역, 유전자 위치 ±500 kb, ±250 kb, 다섯 자료의 공통 변이를 비교했다. 위치 기준은 eQTL 원본 GenePos=81,516,989(GRCh37)이며 GRCh38과 직접 혼합하지 않았다.
- 같은 1000G 기준 빈도에서 FinnGen 세 자료는 각각 동일한 6,579개 변이를 사용했다. 다섯 자료 공통 집합은 4,030개이다. eQTLGen 가용 집합에서는 4,031개이다.

**기본 prior(p12=10⁻⁵)의 PP.H4:**

{mdtable(tt)}

새 FinnGen GD strict의 MR은 OR 0.080(95% CI 0.053–0.122), P=3.19×10⁻³²로 강하지만 H3=0.736, H4=0.264이다. 낮은 prior(p12=10⁻⁶)에서는 H4=0.035이다. BBJ는 기본 H4=0.951에서 낮은 prior에서 0.661로 내려간다. GO strict는 기본 0.891, 낮은 prior 0.451이다.

이는 **질환 연관성 자체와, 그 연관성을 혈액 TSHR 발현이 매개한다는 해석을 구분해야 한다**는 뜻이다. H3가 높아도 복수 원인변이나 LD를 조건부 분석으로 입증한 것은 아니다. 원본에 없는 eQTL 변이, 다중 신호, 조직 차이 등 모든 가능성을 제거한 검사도 아니다.

![TSHR 영역·prior 검증](Figures/Figure_F1_TSHR_windows_and_priors.png)

## 2. 신규 후보 10개 — MR 다음 단계까지 검증

새 발견 후보 9개(APOM, NEU1, IGLV2-18, IGLV3-12, TNF, HIST1H4E, LTA, LTB, C4B)에 빈도 민감도에서만 발견된 HLA-G를 추가했다. 전체 cis 통계에서 5질환 × 4빈도/가용성 조건 × 3prior, 총 600개 설정을 계산했다. **600개 모두 H4<0.8이었다.**

아래는 새 FinnGen GD strict의 기본 빈도·기본 prior 결과이다. HLA-G의 이 표 값은 원래 기준 빈도 결과라 유의하지 않다.

{mdtable(gt)}

기존 14개까지 합하면 이번에 공위치화를 평가한 **24개 유전자 중 새 FinnGen GD에서 H4≥0.8인 유전자는 없었다.** 나머지 2,520개는 공위치화 미평가이므로 전체 2,544개에 대한 음성 결론으로 확장하면 안 된다. 표 전체는 [24유전자 × 5자료 통합표](Tables/Table_F1_24_genes_five_outcomes.csv)에 있다.

APOM과 NEU1은 rs3101018을, IGLV2-18과 IGLV3-12는 rs5759410을 MR 도구로 공유한다. 이들을 서로 독립적인 두 발견으로 세거나 P값만으로 고유한 원인 유전자를 정할 수 없다. MHC의 복잡한 LD에서도 같은 주의가 필요하다.

**BBJ의 IGLV 두 유전자는 해석 제한을 별도로 표시했다.** 핵심 eQTL rs5759410이 없고, 겹치는 변이에서 가장 작은 eQTL P도 약 2.50×10⁻⁴ / 1.49×10⁻⁴이다. 따라서 낮은 H4를 강한 음성 증거로 읽으면 안 된다. UKB와 FinnGen에는 이 핵심 변이가 포함돼 있다. 누락된 주요 eQTL rsID는 원본 파일의 쉼표/세미콜론 별칭도 재검색했으나 복구되지 않았다. 이는 모든 좌표 기반 재매핑 가능성까지 배제한 것은 아니다.

![신규 후보 공위치화](Figures/Figure_F2_new_candidate_colocalization.png)

## 3. C4B·HLA-G — 빈도 민감성이 생기는 이유

새 GD와 GO strict 각각에서 동일한 SNP 집합·질환 효과·발현 Z가 유지됐다. 변화는 발현 β/SE 재구성과 SNP별 가중치에서 발생했으며, 단순히 도구가 추가·삭제된 결과가 아니다. 원본 빈도와 Z를 다시 추출한 값도 18개 유지 도구에서 정확히 일치했다.

{mdtable(ft)}

**HLA-G:** 새 GD에서 rs150254595의 가중치가 21.3%에서 58.9%로 커졌다. 이 변이를 제외하면 원래 빈도 P=0.820, eQTLGen 빈도 P=0.416으로 유의하지 않다. GO strict에서도 제외 후 P=0.821 / 0.929이다. 따라서 eQTLGen 조건에서만 발견된 HLA-G를 안정적인 새 표적으로 부르면 안 된다. 16개 도구가 있어도 효과가 고르게 지지되는 것은 아니다.

일부 희귀 변이의 빈도 차이가 크다. 예를 들어 rs141393161의 효과대립유전자 빈도는 1000G EUR 0.000994, eQTLGen 0.152528이다. 이 차이는 원본에 실제로 존재한다. 표본·희귀 변이 추정·MHC 변이 품질 중 무엇이 주원인인지까지 입증한 것은 아니며, 고유 원인 유전자 주장 전에는 별도 변이 품질 확인이 필요하다.

**C4B:** 두 SNP의 GD 효과 방향은 같지만 크기가 다르다. 기준 빈도에서 rs3131383 단독 OR=1.699, rs9264533 단독 OR=1.195이다. 이질성 P=7.79×10⁻⁷이며, eQTLGen 빈도에서 1.36×10⁻⁸이다. 전체 MR P는 1.10×10⁻⁵에서 1.48×10⁻⁴로 바뀌어 0.05/2,544 기준을 넘나든다. GO strict에서는 두 빈도 모두 스크린 기준을 통과하지만 H4는 기준 미달이다. 불리한 SNP를 제거한 결과를 주 분석으로 바꾸지 않았다.

![HLA-G 가중치와 leave-one-out](Figures/Figure_F3_HLAG_weights_and_leave_one_out.png)

![C4B의 도구별 이질성](Figures/Figure_F4_C4B_instrument_heterogeneity.png)

## 4. 자료 정의 — 원고에 바로잡아야 할 표현

| 자료 | 확인된 내용 | 원고에서의 표현 |
|---|---|---|
| BBJ | GD; 2,809 cases / 172,656 controls; East Asian | GD 발견 자료 |
| UKB GCST90038636 | Hyperthyroidism or thyrotoxicosis; 3,731 / 480,867; ancestry metadata NR | 자가 보고 갑상선기능항진증 비교 자료. 유럽계만 분석했다고 확정하지 않음 |
| FinnGen GD strict | 3,962 / 496,386; 정의에 GRAVES_OPHT 포함 | 추가 GD 분석. 안병증 없는 GD만의 분석이 아님 |
| FinnGen GO | 858 / 499,490 | 기존 안병증 자료, 일반 인구 대조군 |
| FinnGen GO strict | 753 / 499,595; H062&E059 코드 조합 | 동일 자원 내 표현형 민감도 분석 |

UKB 원 논문은 자가 보고 질환을 BOLT-LMM으로 분석했고, 성별·나이·BMI·평가센터·ethnicity·batch·20개 유전 주성분을 보정했다. 이는 우리가 개인자료에서 추가로 임상변수를 보정했다는 뜻은 아니다. 공개 메타데이터의 ancestry는 NR이며, Methods에 EUR-only 표본 제한을 확인하지 못했다. EUR LD reference 사용은 표본이 모두 유럽계라는 증거가 아니다. 원본 선형 β와 SE를 각각 μ(1−μ)로 나눈 근사 log-odds 변환(μ=3,731/484,598)을 유지했으며, 이는 로지스틱 회귀 재분석이 아니다. [원 논문 Methods, pp. 11–17](https://discovery.ucl.ac.uk/id/eprint/10152609/1/EMS122880.pdf), [GWAS Catalog](https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90038636).

원 논문의 후속 분석에서는 MHC를 제외했지만, 우리가 사용한 공개 GWAS 파일에는 해당 영역의 통계가 존재한다. 두 사실을 구분해야 한다.

FinnGen 자료들은 표본이 중첩되므로 독립적인 세 번의 재현이나 군간 효과 차이로 해석하지 않는다. GD 환자끼리 TED 유무를 비교한 분석도 아니다. GWAS API 표본 수를 사용했으며, Risteys의 별도 표현형 집계 숫자로 바꾸지 않았다. 이번 웹 도구의 Risteys 재접속은 실패하여 앞서 같은 날 저장한 공식 정의 원문을 확인했다. [GD 정의](https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_STRICT), [GO strict 정의](https://r12.risteys.finngen.fi/endpoints/E4_GRAVES_OPHT_STRICT), [저장 원문](Verification/Sources/FinnGen_R12_strictGD_definition_text.txt).

## 5. 논문 방향에 대한 판단

이제 중심 메시지는 “TSHR 발현의 인과효과가 여러 자료에서 확정되었다”가 아니라, **“MR 연관성과 공유변이 근거를 구분하면 TSHR 및 신규 후보에 대한 해석이 달라진다”**가 적절하다. TSHR의 OR<1은 유전적으로 예측한 혈액 발현 증가와 낮은 질환 odds의 연관 방향이다. 실제 발현을 치료로 높이면 예방된다는 뜻이나 흉선 면역관용 기전을 입증한 결과가 아니다.

IGF1R의 추가 GD 결과(OR 1.070, P=0.539, H4=0.010)는 이번 연구 조건에서 근거가 제한적이라는 해석을 유지한다. “혈액에서는 원래 확인할 수 없다”, “효과가 없다”, “치료 표적이 아니다” 중 어느 것도 입증하지 않았다.

원고를 개편할 때 필요한 조치는 다음과 같다.

1. **Methods/Table 1:** UKB의 정확한 출처·자가 보고 표현형·확인되지 않은 ancestry 범위·선형 척도 변환을 명시한다. FinnGen 중첩과 endpoint 정의를 구분한다.
2. **Results/Table:** 기존 14개와 새 10개를 선택 경로별로 나누되, 전부 MR OR·CI·P와 H3/H4를 함께 제시한다. 이번 Table F1은 검토용 전체 통합표이며 본문 표 배치는 다음 편집 단계에서 결정한다.
3. **Figures:** 이번 F1–F4를 추가 검증 그림으로 제공한다. 기존 논문의 Figure 번호·파일을 교체한 상태는 아니다. 본문에는 TSHR의 자료별 차이가 보이는 요약 그림, 상세 가중치/LOO는 보충자료가 적절하다.
4. **Discussion/Abstract:** 새 GD의 낮은 TSHR H4와 빈도에 민감한 후보를 숨기지 않고, 발현 매개 인과성과 치료 반응에 대한 단정을 낮춘다.

## 6. 검증 범위와 남은 사항

- R coloc 패키지로 TSHR 240개 + 신규 후보 600개 posterior 벡터를 대조했다. 최대 차이 1.13×10⁻¹³ 미만이다.
- 기존 TSHR 60개 설정을 재현했다. 이번 새 설정은 TSHR 영역 180개 + 신규 후보 600개이다.
- full MR 8개와 leave-one-out 72개를 R 패키지 및 별도 Python 공식으로 검증했다. 최대 차이 8.89×10⁻¹⁶ 미만이다.
- 모든 결과와 불리한 조건을 보존했다. PP.H4 기준을 결과에 맞춰 변경하거나 유전자를 선택적으로 제거하지 않았다.

이번 점검은 단일 인과변이 가정의 한계, 발현 조직 적합성, MHC에서의 고유 유전자 구분을 해결하지는 않는다. SuSiE·다중신호 공위치화는 시행하지 않았다. 원래 MR의 도구 선택/LD clumping을 독립된 최신 reference로 다시 검증한 작업도 아니다. 이 한계들을 풀어야만 가능한 강한 인과 주장을 지금 추가하지 않는다.

**분석 검증본은 준비되었지만, 기존 제출 원고가 이 결과까지 반영된 상태는 아니다.** 이 보고서는 제출 완료 판정이 아니다.

## 파일 안내

- [전체 24개 × 5자료 표](Tables/Table_F1_24_genes_five_outcomes.csv)
- [TSHR 모든 영역·prior·빈도 결과](Tables/Table_F2_TSHR_windows_priors_frequencies.csv)
- [신규 10개 기본 MR/공위치화](Tables/Table_F3_new_candidates_MR_coloc.csv)
- [C4B/HLA-G 빈도별 MR](Tables/Table_F4_C4B_HLAG_frequency_MR.csv)
- [leave-one-out 전체 72개](Tables/C4B_HLAG_leave_one_out_72.csv)
- [통계 열 설명](DATA_DICTIONARY_KO.md)
- [검증 기록](Verification/final_checks_v1/verification.json)

그림은 Figures 폴더에 PNG와 벡터 PDF로 각각 저장했다. 계산 코드는 Verification/Scripts에 보존했다. 출처 원본 대용량 파일은 포함하지 않는다.
'''
    (O/'REPORT_KO.md').write_text(report,encoding='utf8')
    dictionary='''# 상세 결과 열 설명

이 패키지는 추가 검증용이다. Table_F1은 24개 유전자 × 5질환=120행, 전체 공위치화 설정 파일은 24×5×4×3=1,440행이다. TSHR 별도 파일의 240행에는 기준 전체 cis 설정 60행이 중복 포함되므로 합산할 때 주의한다.

| 열 | 뜻 |
|---|---|
| gene / gene_symbol | 유전자 기호 |
| outcome | 질환 GWAS 데이터셋. 치료 결과가 아님 |
| n_iv | 실제 MR에 유지된 도구 SNP 수 |
| beta / se | MR log OR와 표준오차 |
| OR / CI_lower / CI_upper | OR와 95% CI; 근사 표준화 혈액 발현 단위 |
| pvalue | 양측 MR P값. Q_pvalue와 다름 |
| estimable | 도구 매칭 후 MR 추정 가능 여부. 불가는 음성이 아님 |
| P_lt_0_05_over_2544 | 전체 스크린 보정 기준 0.05/2,544의 충족 여부 |
| n_overlap | 공위치화에 들어간 영역 변이 수; n_iv가 아님 |
| PP.H0 | posterior probability of hypothesis 0: 두 형질 모두 영역 연관 없음 |
| PP.H1 | 이 코드의 dataset1인 질환만 영역 연관 |
| PP.H2 | dataset2인 혈액 발현만 영역 연관 |
| PP.H3 | 두 형질이 서로 다른 원인변이를 가진다는 단일 원인변이 모델 설명 |
| PP.H4 | 두 형질이 같은 원인변이를 공유한다는 설명 |
| p12 | 변이 하나가 두 형질 모두와 연관될 사전확률. P값이 아님 |
| top_snp | H4를 가정했을 때 최대 공동 ABF 변이. 원인 확정 변이가 아님 |
| H4_low_prior | p12=10^-6, 기준 빈도에서 H4 |
| H4_eqtlgen | p12=10^-5, eQTLGen 가용 변이에서 H4 |
| original_reference | 1000G EUR 빈도, 해당 빈도가 있는 변이 |
| paired_reference / paired_eqtlgen | 두 빈도가 모두 있는 동일 SNP 집합에서 빈도만 교체 |
| available_eqtlgen | eQTLGen 빈도가 있는 모든 호환 변이 |
| reharmonized_eqtlgen | MR에서 eQTLGen 빈도로 대립유전자를 다시 정렬한 조건 |
| Q / Q_df / Q_pvalue | SNP별 효과의 이질성 통계, 자유도, P값 |
| weight_fraction | bx²/sy²의 합을 1로 정규화한 IVW 기여도; 원인일 확률이 아님 |
| Q_contribution | 해당 SNP의 IVW 잔차 이질성 기여 |
| F_exposure | (발현 beta / 발현 SE)² |
| removed_snp | leave-one-out에서 제외한 SNP; none은 full model |
| eaf_ref / eaf_new | 효과대립유전자 빈도. minor allele frequency와 항상 같지 않음 |
| source_lead_status | 가장 강한 원본 eQTL이 매칭 집합에 포함되는지 |

H0–H4는 입력 형질 순서에 따라 H1/H2의 의미가 바뀐다. 여기서는 질환이 dataset1이다. 전체 영역, 단일 원인변이 가정, p1=p2=10^-4, p12=10^-5/5×10^-6/10^-6, 양적 발현 sdY=1, disease type=cc를 사용했다. beta/varbeta를 사용하므로 P값만을 입력하지 않았다. 표에서 0으로 표시된 매우 작은 수에는 부동소수점 underflow가 포함될 수 있다. 최대 signed Z를 검증했고 0 P값을 이유로 SNP를 제외하지 않았다.

일반 MR P<0.05, 스크린 P<0.05/2,544, PP.H4≥0.8은 서로 다른 판단 기준이다. H4 0.8은 이 연구의 공유변이 지지 기준이며, “유의확률 80%”나 “인과효과가 참일 확률 80%”를 뜻하지 않는다. [coloc 공식 설명](https://chr1swallace.github.io/coloc/articles/vignette.html).
'''
    (O/'DATA_DICTIONARY_KO.md').write_text(dictionary,encoding='utf8')
    print('Report and figures created:',O)
if __name__=='__main__':main()
