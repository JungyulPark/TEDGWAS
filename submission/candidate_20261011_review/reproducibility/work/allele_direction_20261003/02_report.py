from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, re, shutil, zipfile, sys
import pandas as pd

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TED_Allele_Direction_20261003'
T=OUT/'Tables'
result=json.loads((OUT/'AUDIT_RESULT.json').read_text())
assert result['all_pass'] and result['checks']==73
g=pd.read_csv(T/'Table_A3_both_variants_five_disease_GWAS.csv')
m=pd.read_csv(T/'Table_A4_MR_reproduction.csv')
ld=pd.read_csv(T/'Table_A2_signed_LD.csv')
eq=pd.read_csv(T/'Table_A1_eQTL_original_direction.csv')
checks=list({r['check']:r for r in pd.read_csv(OUT/'verification_checks.csv').to_dict('records')}.values())
def check(name,ok,detail=''):
    assert ok,(name,detail)
    checks[:]=[r for r in checks if r['check']!=name]
    checks.append(dict(check=name,pass_check=bool(ok),detail=str(detail)))

# Make the PLINK default A1 coding and haplotype/genotype distinction explicit.
hap=[]
for idx,row in ld.iterrows():
    pop=row.population
    pref=Path(f'C:/ProjectTEDGWAS/TrackA_MR/v5_upgrade/data/1000G_{pop}_phase3/{pop}_chr14')
    alleles={}
    for line in pref.with_suffix('.bim').open():
        f=line.split()
        if f[1] in ['rs179252','rs179247']:alleles[f[1]]=f[4:6]
    pr=pd.read_csv(OUT/f'LD_verification/{pop}_signed.ld',sep=r'\s+').iloc[0]
    oriented=pr.R*(1 if alleles['rs179252'][0]=='G' else -1)*(1 if alleles['rs179247'][0]=='G' else -1)
    check(f'{pop} PLINK default A1 coding reconciled',abs(oriented-row.signed_G_dosage_r)<5e-7)
    ld.loc[idx,'plink_A1_rs179252']=alleles['rs179252'][0]
    ld.loc[idx,'plink_A1_rs179247']=alleles['rs179247'][0]
    ld.loc[idx,'plink_raw_A1_correlation']=pr.R
    text=(OUT/f'LD_verification/{pop}_haplotype.log').read_text()
    for h,f in re.findall(r'^\s+(TG|GG|TA|GA)\s+([0-9.]+)\s+[0-9.]+\s*$',text,re.M):
        hap.append(dict(population=pop,rs179252_allele=h[0],rs179247_allele=h[1],haplotype=h,estimated_frequency=float(f),method='PLINK unphased two-locus maximum likelihood; log output rounded to six decimals'))
    ph=[x for x in hap if x['population']==pop]
    check(f'{pop} four haplotypes retained',len(ph)==4 and abs(sum(x['estimated_frequency'] for x in ph)-1)<2e-6)
    check(f'{pop} most common two haplotypes are GG and TA',set(x['haplotype'] for x in sorted(ph,key=lambda x:-x['estimated_frequency'])[:2])=={'GG','TA'})
ld.to_csv(T/'Table_A2_signed_LD.csv',index=False)
pd.DataFrame(hap).to_csv(T/'Table_A2b_inferred_haplotype_frequencies.csv',index=False)

check('All ten disease rows retained',len(g)==10 and g.groupby('outcome').size().eq(2).all())
check('All ten G-allele disease estimates lower odds',g.aligned_G_beta_logOR.lt(0).all())
check('Both G alleles increase blood TSHR',eq.AssessedAllele.eq('G').all() and eq.Zscore.gt(0).all())
plan=pd.read_csv(OUT/'Tissue_planned_results_50_cells.csv')
check('50 planned cells including ten primary cells',len(plan)==50 and plan.family.eq('primary').sum()==10 and plan.status.eq('not_started').all())
freeze=json.loads((OUT/'PROTOCOL_FREEZE.json').read_text())
check('Frozen protocol unchanged at reporting',hashlib.sha256((OUT/'PROTOCOL_v1.md').read_bytes()).hexdigest()==freeze['sha256'])

names={'BBJ':'BBJ: GD','E4_GRAVES_STRICT':'FinnGen: GD','UKB':'UKB: 갑상선기능항진증','FinnGen_GO':'FinnGen: 안병증','E4_GRAVES_OPHT_STRICT':'FinnGen: 엄격한 안병증 정의'}
rows=[]
for oc,name in names.items():
    a=g[g.outcome.eq(oc)&g.snp.eq('rs179252')].iloc[0]
    b=m[m.outcome.eq(oc)].iloc[0]
    rows.append(f"| {name} | {a.aligned_G_OR:.3f} ({a.aligned_G_CI_low:.3f}–{a.aligned_G_CI_high:.3f}) | {b.MR_OR:.3f} ({b.MR_CI_low:.3f}–{b.MR_CI_high:.3f}) | {b.MR_p:.3g} |")
e1=eq[eq.SNP.eq('rs179252')].iloc[0];e2=eq[eq.SNP.eq('rs179247')].iloc[0]
report=f'''# TSHR 대립유전자 방향 검증 — 완료 보고서

**판정: 혈액 TSHR 발현 증가 방향, 질환 odds 감소 방향, 기존 흉선 연구의 보호 대립유전자 방향이 일치합니다. 기존 주요 OR·P값은 재현됐습니다.**

검증일: 2026-10-03. 이 작업은 기존 연구의 사후 검증입니다. 후속 조직 분석 계획만 새로운 결과를 보기 전에 고정했으며, 원 연구 전체를 사전등록한 것으로 표현하지 않습니다.

## 1. 이번에 실제로 확인한 내용

원본 eQTLGen 압축 파일과 다섯 질환 GWAS 압축 파일에서 두 SNP의 행을 새로 추출했습니다. 기존 정리표만 읽고 결론 내리지 않았습니다. 1000 Genomes EUR 503명과 EAS 504명의 원본 BED 유전형에서 G 대립유전자 개수를 직접 계산하고 PLINK의 별도 계산과 대조했습니다. 혈액 발현의 beta/SE와 Wald 비율을 재계산했습니다.

| 확인 단계 | 관찰된 결과 | 의미 |
|---|---|---|
| 혈액 도구변수 rs179252 | G 대립유전자: Z = +{e1.Zscore:.4f}, P = {e1.Pvalue:.4g}, N = {int(e1.NrSamples):,} | G가 혈액 TSHR 발현 증가 방향 |
| 흉선 논문의 rs179247을 혈액 자료에서 직접 조회 | G 대립유전자: Z = +{e2.Zscore:.4f}, P = {e2.Pvalue:.4g}, N = {int(e2.NrSamples):,} | 같은 G 대립유전자가 혈액에서도 발현 증가 방향 |
| rs179252-G와 rs179247-G의 연결 | EUR: r = +{ld.iloc[0].signed_G_dosage_r:.6f}, r² = {ld.iloc[0].r_squared:.6f}; EAS: r = +{ld.iloc[1].signed_G_dosage_r:.6f}, r² = {ld.iloc[1].r_squared:.6f} | 두 인구집단 모두 G–G 및 T–A 조합이 우세 |
| 질환 GWAS | 두 SNP 모두 다섯 항목에서 G 대립유전자 OR < 1 | 대립유전자 정렬 후 질환 연관 방향도 일치 |

흉선 논문은 rs179247-G를 보호 대립유전자로 분류하고, 이 대립유전자의 상대 전사가 흉선에서 더 높다고 보고했습니다. 이 근거와 위 검증을 연결하면 **혈액에서 발현을 높이는 유전적 방향이 흉선 연구의 보호 방향과 일치한다**고 쓸 수 있습니다. [Marín-Sánchez et al., 2019](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2019.01695/full), [Colobran et al., 2011](https://pubmed.ncbi.nlm.nih.gov/21642385/)

두 SNP는 서로 강하게 연결돼 있으므로, 두 SNP의 결과를 별개의 독립적 재현으로 세면 안 됩니다. 참조자료의 LD로 추정한 조합이며 환자에서 직접 관찰한 흉선 발현이나 매개효과는 아닙니다.

## 2. 기존 MR 수치가 바뀌는가?

**바뀌지 않습니다.** 다음은 모두 같은 rs179252-G 방향으로 정렬한 값입니다. 첫 OR은 G 대립유전자 한 개당 질환 연관성이고, 두 번째 OR은 이를 혈액 발현 단위로 환산한 MR 결과입니다.

| 질환 자료 | SNP OR (95% CI), G 한 개당 | 기존 MR OR (95% CI), 재구성 발현 약 1 SD당 | 재계산 MR P |
|---|---|---|---|
{chr(10).join(rows)}

MR OR이 SNP OR보다 훨씬 작은 이유는 질환의 log OR을 작은 발현 효과계수(beta = {m.iloc[0].beta_exposure:.6f})로 나누기 때문입니다. 이를 TSHR를 실험적으로 증가시켰을 때의 치료 효과로 바꿔 해석할 수 없습니다. rs179252의 F 통계량은 {m.iloc[0].exposure_F:.2f}입니다.

BBJ/UKB/기존 안병증 결과의 기존 CI 계산은 1.96, 추가 FinnGen GD/엄격한 안병증 결과는 정규분포의 정확한 97.5% 분위수를 사용했습니다. 각 원래 계산 규칙으로 확인했으며 표시되는 CI나 해석에는 차이가 없습니다. 원본 GWAS P와 beta/SE에서 다시 계산한 MR P를 별도 열로 보존했습니다. 특히 UKB 원본 GWAS P를 재계산 MR P와 동일한 값으로 덮어쓰지 않았습니다.

## 3. 논문에서 강화할 부분과 유지할 부분

- **Results에 추가:** 효과 대립유전자 확인 결과를 짧게 보고하고 상세 검증표를 연결합니다.
- **Discussion을 보강:** 지금까지는 단순히 두 SNP가 높은 LD라는 점만 확인했습니다. 이제 G–G 연결과 혈액/문헌 간 방향 일치까지 보고할 수 있습니다. 흉선 면역관용 가설과 방향상 부합한다는 표현이 적절합니다.
- **Supplementary Methods에 추가:** 원본 행 추출, GRCh37/GRCh38 구분, G 대립유전자 기준 정렬, LD 두 방법 대조 및 해석 한계를 기록합니다.
- **기존 본문 Table 1–3와 Figure 1–3의 주요 수치:** 이번 검증으로 수정할 수치는 없습니다. 새로운 상세자료는 아래 검증표에 분리했습니다. 이전 제출 패키지를 덮어쓰지 않았습니다.

이번 결과가 흉선 발현을 통한 GD 예방 기전을 증명하거나, GD 환자 중 TED 발생 원인을 밝히거나, 갑상선/안와의 TSHR를 높이는 치료를 지지하는 것은 아닙니다. 보호 방향에 대한 **대립유전자 수준의 일관성**을 추가한 것입니다. 출판 논문의 같은 변이 정보와 우리 데이터를 연결한 검증으로, 새로운 흉선 실험을 수행한 것은 아닙니다.

## 4. 다음 조직 분석은 결과를 보기 전에 이렇게 고정했습니다

TSHR·IGF1R 모두 GTEx v8의 **갑상선, 배양 섬유모세포, 피하 지방, 내장 지방, 전혈**에서 조사합니다. 섬유모세포/지방은 안와의 해당 조직과 같다고 간주하지 않습니다. 전혈은 기존 eQTLGen과 비교하는 기준입니다.

| 판단 항목 | 고정한 기준 |
|---|---|
| 1차 비교 | 두 유전자 × 다섯 조직 × FinnGen GD = 10개; 유럽계 노출/질환 조합 |
| 1차 다중검정 기준 | P < 0.005 |
| 2차 비교 | 같은 10조합 × BBJ GD·UKB·두 FinnGen 안병증 정의 = 40개 |
| 2차 다중검정 기준 | P < 0.00125; BBJ는 인종 간 이전 가능성 검토 |
| 도구변수 | cis P < 5×10⁻⁸, F > 10, LD r² < 0.001; 기준을 낮춰 유의한 결과를 만들지 않음 |
| 분석 | SNP 1개는 Wald ratio, 2개 이상은 IVW; OR·95% CI·P·도구 수·강도를 모두 보고 |
| 해석 | 방향 일치/반대와 통계적 근거의 강도를 별도로 표기 |
| 분석 불가 | 도구 부족·자료 미공개·대립유전자 불명확 등을 이유와 함께 남김; 무효과로 표현하지 않음 |
| 공위치화 | 전체 영역 요약통계가 있을 때만 기존 단일 원인변이 모델로 보조 해석; 유의 SNP만으로 실행하지 않음 |

자료가 없는 조합도 삭제하지 않는 **50행 결과 틀**을 만들었습니다. 한 조직은 유의하고 다른 조직은 비유의하다는 이유만으로 조직 간 차이를 주장하지 않습니다. 반대 방향 결과도 같은 기준으로 보고합니다. FinnGen 자료들이 겹치므로 독립된 반복 연구로 합산하지 않습니다.

ImmuNexUT은 28개 면역세포 유형의 공개 범위와 이용 조건을 먼저 확인하는 별도 후속 단계입니다. 실행 전 구체적 세포 목록과 분석군·다중검정 기준을 추가로 고정해야 합니다. 이번에 면역세포 MR이나 다중신호 공위치화를 실행한 것은 아닙니다.

고정 시각(UTC): {freeze['frozen_at_utc']}. SHA-256: `{freeze['sha256']}`. 이후 변경은 원본 계획을 보존한 수정 기록으로 남깁니다. **조직별 MR 결과는 아직 생성하지 않았습니다.**

## 5. 파일 안내

- `Tables/Table_A1_eQTL_original_direction.csv`: 두 SNP의 원본 혈액 eQTL 방향.
- `Tables/Table_A2_signed_LD.csv`: EUR/EAS의 대립유전자 이름을 명시한 signed LD.
- `Tables/Table_A2b_inferred_haplotype_frequencies.csv`: 우세 조합과 드문 반대 조합을 모두 보존한 추정 빈도.
- `Tables/Table_A3_both_variants_five_disease_GWAS.csv`: 두 SNP × 다섯 질환, 원본 방향과 G 정렬 방향을 모두 보존.
- `Tables/Table_A4_MR_reproduction.csv`: 기존 MR beta/SE/OR/CI/P 재현.
- `Manuscript_insertions.md`: 현 원고 참고문헌 번호에 맞춘 Methods/Results/Discussion 삽입 문장.
- `PROTOCOL_v1.md`, `PROTOCOL_FREEZE.json`: 후속 분석 사전 계획 및 해시.
- `Tissue_planned_results_50_cells.csv`: 아직 시작하지 않은 조직 분석의 전체 결과 틀.
- `Source_extracts/`, `LD_verification/`, `SOURCE_MANIFEST.json`, `verification_checks.csv`: 원본 추출 행, 계산 로그, 해시와 검사 결과.

통계 원본에는 GRCh37의 rs179252 위치가 81,435,985, rs179247 위치가 81,432,546입니다. 질환 원본의 GRCh38 위치는 각각 80,969,641과 80,966,202입니다. 다른 build의 좌표를 같다고 비교하지 않았습니다. PLINK가 기본으로 세는 A1은 EUR rs179252에서 T라서 원시 상관계수는 음수입니다. 본 보고서는 두 SNP 모두 **G 개수**로 맞췄으며, 이때 두 집단 모두 양수입니다. PLINK의 haplotype 기반 r²와 genotype 상관의 제곱은 계산법이 달라 조금 다를 수 있어 혼용하지 않았습니다. [PLINK LD 문서](https://www.cog-genomics.org/plink/1.9/ld)
'''
(OUT/'READ_RESULTS_KO.md').write_text(report,encoding='utf-8')

insert='''# Manuscript integration text — allele-direction audit

Prepared against TED_Research_Contribution_20261003/Manuscript.md. Existing reference numbers 27 (Colobran 2011) and 28 (Marín-Sánchez 2019) remain correct. These are ready-to-integrate additions/replacements; the delivered manuscript DOCX/PDF has not been silently replaced. Merge the audit tables into the single detailed supplementary dataset at the next consolidated build, preserving Supplementary Tables S1–S3 and Figures S1–S2. A1–A4 here are audit-file identifiers, not new manuscript table numbers.

## Methods addition (post hoc audit, not a new primary hypothesis)

In a post hoc allele-direction audit, original eQTLGen and disease summary statistics were re-extracted for the TSHR instrument rs179252 and the thymic-expression-associated variant rs179247. Disease effects were aligned to the blood expression-increasing allele. Allele-labelled genotype correlations in 1000 Genomes European and East Asian reference samples were checked by direct genotype decoding and independently using PLINK. This analysis assessed directional consistency with published thymus findings, rather than testing thymic mediation.

## Results addition

Allele verification confirmed that rs179252-G increased blood TSHR expression and was associated with lower disease odds in all five outcomes. It was positively correlated with rs179247-G in European (r = 0.978; r² = 0.957; n = 503) and East Asian (r = 0.996; r² = 0.992; n = 504) reference samples. The latter allele also increased blood TSHR expression and corresponds to the protective allele preferentially transcribed in the published thymus study [28]. Recalculated MR estimates reproduced the existing results. These linked variants do not constitute independent genetic replications.

## Replacement for the current Discussion paragraph beginning “The lower-odds direction…”

The lower-odds direction for TSHR need not contradict the pathogenic role of receptor-stimulating antibodies. Blood transcript abundance, tissue-specific expression and receptor activation are different exposures. Allele verification showed that the blood expression-increasing rs179252-G allele was linked to rs179247-G in both European and East Asian reference samples. The latter is the GD-protective allele associated with greater thymic transcription in published studies [27, 28] and also increased blood expression in eQTLGen. This directional consistency strengthens the plausibility of a central-tolerance explanation, but does not demonstrate that thymic expression mediates the observed disease association. The variants remain correlated markers of a locus rather than proven causal regulators. The marked departure of the MR OR from 1 also reflects division by a small reconstructed expression coefficient. It should not be translated into the effect of experimentally increasing TSHR expression.

## Supplementary Methods addition

The audit used the original eQTLGen 2019-12-11 full cis-eQTL file and the same five disease GWAS sources as the main analysis. The GRCh37 positions of rs179252 and rs179247 were 14:81435985 and 14:81432546; their corresponding GRCh38 positions in disease data were 14:80969641 and 14:80966202. Both allele pairs were non-palindromic and matched between sources. UKB reported the opposite effect allele at both SNPs; its beta was reversed after the existing linear-to-log-odds conversion. No sign reversal was needed for BBJ or FinnGen. Reconstructed exposure coefficients and first-order Wald estimates were independently recalculated using the original reference-frequency scaling. Source GWAS P values and recalculated normal-approximation MR P values were retained separately.

European (n = 503) and East Asian (n = 504) reference samples had no missing genotypes at either variant. Correlations were explicitly oriented to G-allele counts at both SNPs, yielding r = 0.978322812 and 0.995803919, respectively. The corresponding genotype-based r² values were 0.957115524 and 0.991625445. PLINK's default A1 correlation in Europeans was negative because its A1 for rs179252 was T; reorientation to G resolved this apparent sign difference. Independent PLINK two-locus maximum-likelihood haplotype estimates supported G–G and T–A coupling. These are phase estimates from unphased reference genotypes, not directly measured phased haplotypes in the study participants. Haplotype-based and genotype-based r² were not interchanged. This audit does not resolve multiple causal signals, identify a causal variant, or establish thymic or orbital mediation.

## Integration checklist for the next manuscript build

- Add Methods/Results text and replace the relevant Discussion paragraph.
- Replace the earlier blanket statement that the matching direction is unestablished with the specific limitation that direction is consistent but mediation and causal identity are untested.
- Add A1–A4 audit tables as worksheets/sections in the existing single detailed supplementary dataset, with an updated data dictionary; do not proliferate separate numbered supplementary tables.
- No numeric changes to main Tables 1–3 or Figures 1–3 are required by this audit.
- Label the tissue protocol as a planned extension until the 50-cell grid is populated; do not write that replication has already been performed.
'''
(OUT/'Manuscript_insertions.md').write_text(insert,encoding='utf-8')

notes='''# Execution notes

The first audit attempt stopped before writing final results because an overly strict CI identity check assumed qnorm(0.975) for every historical result. The BBJ, UKB and original GO pipeline used 1.96; the later FinnGen GD and strict GO pipeline used qnorm(0.975). The audit was corrected to reproduce each documented historical convention. Beta, SE, OR and P calculations were unchanged. Freshly extracted raw rows were reused after this downstream failure; source files were SHA-256 hashed during recovery. The final successful audit retained both raw GWAS P and recomputed MR P.

No tissue-specific GTEx or ImmuNexUT association results were newly retrieved or analysed. No SuSiE or MVMR was run. The protocol is prospective only for this extension, not for the original manuscript or historical exploratory tissue analyses.

The supporting reference panel consists of unphased PLINK genotypes derived from 1000 Genomes phase 3. “Phase 3” is a project release label and does not mean that phased haplotypes are retained in BED files. Haplotype coupling reported here is inferred from those genotypes. The primary LD values are genotype correlations squared; the haplotype estimates are separately labelled.

Full-source extraction in the first successful input pass stopped after both target SNP rows were found: eQTLGen 2,397,210 rows, BBJ 9,971,383, UKB 7,751,092, FinnGen GO 16,383,732, FinnGen GD and strict GO 16,383,733 each. This is not an exhaustive duplicate-variant scan; earlier full-source and harmonization audits supply that context. Both exact allele pairs and source rows were independently compared with the prior extracts.

External sources consulted: eQTLGen cis-eQTL page (https://eqtlgen.org/cis-eqtls.html); Marín-Sánchez 2019 (https://doi.org/10.3389/fimmu.2019.01695); Colobran 2011 PubMed abstract (https://pubmed.ncbi.nlm.nih.gov/21642385/); PLINK LD and file-format documentation (https://www.cog-genomics.org/plink/1.9/ld, https://www.cog-genomics.org/plink/1.9/formats#bed). The eQTLGen README download could not be read through the web viewer; local raw column names and the official page's allele/Z-score definitions were retained as evidence. No TLS exception was applied in this task.
'''
(OUT/'EXECUTION_NOTES.md').write_text(notes,encoding='utf-8')
pd.DataFrame(checks).to_csv(OUT/'verification_checks.csv',index=False)
result['report_checks']=len(checks)
result['interpretation']='Directionally consistent with published thymic protective allele; no mediation established'
(OUT/'AUDIT_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
code=OUT/'Reproduce';code.mkdir(exist_ok=True)
for p in Path(__file__).parent.glob('*.py'):shutil.copy2(p,code/p.name)
(code/'README.md').write_text('''# Reproducing this audit

These are archived copies of the executed scripts. Their original working location is `<workspace>/work/allele_direction_20261003/`; copy them there before running. The workspace contains the prior canonical analysis inputs in `work/` and `outputs/`. Original large data are read from `C:/ProjectTEDGWAS`; those files are not bundled. Their SHA-256 hashes are in SOURCE_MANIFEST.json.

1. Preserve the delivered package before a rerun, since audit outputs use the same versioned directory.
2. Do not rerun freeze_protocol.py over the frozen plan. Its intentional existence guard prevents replacing the pre-analysis record.
3. Run 01_audit.py to extract the original raw records again, compute LD, and verify MR. The `--reuse-extracts` option is only for recovery using the exact already-generated raw extracts.
4. Run 02_report.py to build the interpretive report, final checks, and ZIP.

The original protocol-freeze and outcome files are provenance records. A future tissue analysis or amended protocol must use a new version and retain these files. Dependencies: Python with numpy and pandas; the configured PLINK 1.9 executable. The local run used the bundled Codex Python; no new packages were installed.
''',encoding='utf-8')
manifest=[]
for p in sorted(OUT.rglob('*')):
    if p.is_file() and p.name!='PACKAGE_MANIFEST.json':manifest.append(dict(path=p.relative_to(OUT).as_posix(),bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
(OUT/'PACKAGE_MANIFEST.json').write_text(json.dumps({'created_utc':datetime.now(timezone.utc).isoformat(),'files':manifest},indent=2),encoding='utf-8')
zip_path=OUT.with_suffix('.zip')
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.rglob('*')):
        if p.is_file():z.write(p,p.relative_to(OUT.parent))
with zipfile.ZipFile(zip_path) as z:
    assert z.testzip() is None
    for entry in manifest:assert hashlib.sha256(z.read(OUT.name+'/'+entry['path'])).hexdigest()==entry['sha256']
print(json.dumps({'report_checks':len(checks),'all_pass':all(x['pass_check'] for x in checks),'files':len(manifest)+1,'zip':str(zip_path),'zip_bytes':zip_path.stat().st_size,'zip_sha256':hashlib.sha256(zip_path.read_bytes()).hexdigest()},indent=2))
