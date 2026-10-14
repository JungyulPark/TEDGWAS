from pathlib import Path
from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor
import csv, gzip, hashlib, json, math, subprocess, sys
import numpy as np
import pandas as pd
from statistics import NormalDist
class norm:
    ppf = staticmethod(NormalDist().inv_cdf)
    sf = staticmethod(lambda x: 0.5 * math.erfc(x / math.sqrt(2)))

sys.stdout.reconfigure(encoding='utf-8')
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'outputs/TED_Allele_Direction_20261003'
RAW=OUT/'Source_extracts'; RAW.mkdir(exist_ok=True)
LD=OUT/'LD_verification'; LD.mkdir(exist_ok=True)
TABLE=OUT/'Tables'; TABLE.mkdir(exist_ok=True)
SRC=Path('C:/ProjectTEDGWAS')
SNPS=['rs179252','rs179247']
checks=[]; manifests=[]
def check(name,passed,detail=''):
    checks.append(dict(check=name,pass_check=bool(passed),detail=str(detail)))
    if not passed: raise AssertionError((name,detail))
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def record(p,role):
    return dict(path=str(p),role=role,bytes=p.stat().st_size,mtime_utc=datetime.fromtimestamp(p.stat().st_mtime,timezone.utc).isoformat(),sha256=sha(p))

freeze=json.loads((OUT/'PROTOCOL_FREEZE.json').read_text())
check('Protocol unchanged',sha(OUT/'PROTOCOL_v1.md')==freeze['sha256'])

# Two independent implementations: direct BED decoding and PLINK 1.9 reports.
(LD/'variants.txt').write_text('\n'.join(SNPS)+'\n',encoding='ascii')
(LD/'counted_alleles.txt').write_text('rs179252 G\nrs179247 G\n',encoding='ascii')
ldrows=[]; counts=[]
plink=SRC/'tools/plink.exe'
manifests.append(record(plink,'PLINK executable'))
for pop in ['EUR','EAS']:
    pref=SRC/f'TrackA_MR/v5_upgrade/data/1000G_{pop}_phase3/{pop}_chr14'
    n=sum(1 for line in pref.with_suffix('.fam').open() if line.strip())
    variants={}
    with pref.with_suffix('.bim').open() as f:
        for idx,line in enumerate(f):
            row=line.split()
            if row[1] in SNPS:variants[row[1]]=(idx,row)
    check(f'{pop} both variants present',len(variants)==2)
    width=(n+3)//4
    gdos={}
    with pref.with_suffix('.bed').open('rb') as f:
        check(f'{pop} BED header',f.read(3)==b'\x6c\x1b\x01')
        for snp,(idx,row) in variants.items():
            f.seek(3+idx*width);chunk=f.read(width)
            # BED 00=A1/A1, 01=missing, 10=heterozygous, 11=A2/A2.
            dose=np.array([{0:2,1:np.nan,2:1,3:0}[(chunk[i//4]>>(2*(i%4)))&3] for i in range(n)],dtype=float)
            check(f'{pop} {snp} G allele present','G' in row[4:6])
            if row[4]!='G': dose=2-dose
            gdos[snp]=dose
    valid=np.isfinite(gdos[SNPS[0]])&np.isfinite(gdos[SNPS[1]])
    x,y=gdos[SNPS[0]][valid],gdos[SNPS[1]][valid]
    r=float(np.corrcoef(x,y)[0,1]);r2=r*r
    for a in range(3):
        for b in range(3):counts.append(dict(population=pop,rs179252_G_copies=a,rs179247_G_copies=b,n=int(((x==a)&(y==b)).sum())))
    old=pd.read_csv(ROOT/f'work/local_evidence/major_review_20260912/ld/{pop}.ld',sep=r'\s+')
    old=old[old.SNP_A.eq('rs179247')&old.SNP_B.eq('rs179252')].iloc[0]
    check(f'{pop} old r2 reproduced',abs(r2-old.R2)<5e-7,{'new':r2,'old':old.R2})
    base=[str(plink),'--bfile',str(pref),'--extract',str(LD/'variants.txt'),'--keep-allele-order','--threads','2','--memory','1024']
    def run(extra,prefix):
        cmd=base+extra+['--out',str(LD/prefix)]
        p=subprocess.run(cmd,capture_output=True,text=True,encoding='utf-8',errors='replace')
        (LD/f'{prefix}.console.txt').write_text(p.stdout+p.stderr,encoding='utf-8')
        check(f'{prefix} PLINK execution',p.returncode==0,p.stderr)
    run(['--r','in-phase','with-freqs','--ld-window-r2','0'],f'{pop}_signed')
    run(['--ld',*SNPS],f'{pop}_haplotype')
    run(['--recode','A','--recode-allele',str(LD/'counted_alleles.txt')],f'{pop}_G_dosage')
    dr=pd.read_csv(LD/f'{pop}_G_dosage.raw',sep=r'\s+')
    for snp in SNPS: check(f'{pop} {snp} direct BED equals PLINK G dosage',np.array_equal(dr[snp+'_G'].to_numpy(),gdos[snp],equal_nan=True))
    pr=float(dr['rs179252_G'].corr(dr['rs179247_G']))
    check(f'{pop} signed r independent calculation',abs(r-pr)<1e-12)
    check(f'{pop} no missing pair genotypes',int(valid.sum())==n)
    # Retain actual G/G sign, not PLINK's potentially different A1/A1 sign.
    ldrows.append(dict(population=pop,n_samples=n,n_complete=int(valid.sum()),n_missing=int(n-valid.sum()),rs179252_counted_allele='G',rs179247_counted_allele='G',signed_G_dosage_r=r,r_squared=r2,rs179252_G_frequency=float(x.mean()/2),rs179247_G_frequency=float(y.mean()/2),coupling='G-G / T-A' if r>0 else 'G-A / T-G',phase_basis='inferred from unphased genotype LD; PLINK haplotype maximum likelihood available',rs179252_pos37=int(variants['rs179252'][1][3]),rs179247_pos37=int(variants['rs179247'][1][3])))
    for ext in ['.bed','.bim','.fam']:manifests.append(record(pref.with_suffix(ext),f'{pop} reference genotypes'))
    print('LD verified:',ldrows[-1],flush=True)
pd.DataFrame(ldrows).to_csv(TABLE/'Table_A2_signed_LD.csv',index=False)
pd.DataFrame(counts).to_csv(TABLE/'Reference_two_locus_genotype_counts.csv',index=False)

sources=[('eQTLGen',SRC/'TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz','SNP'),
('BBJ',SRC/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90018627_harmonised.tsv.gz','rsid'),
('UKB',SRC/'TrackA_MR/v5_upgrade/04_druggable_mr/data/outcomes/GCST90038636_harmonised.tsv.gz','hm_rsid'),
('FinnGen_GO',SRC/'finngen_R12_GRAVES_OPHT.gz','rsids'),
('E4_GRAVES_STRICT',SRC/'finngen_R12_E4_GRAVES_STRICT.gz','rsids'),
('E4_GRAVES_OPHT_STRICT',SRC/'finngen_R12_E4_GRAVES_OPHT_STRICT.gz','rsids')]

def extract(spec):
    name,src,idcol=spec
    dest=RAW/f'{name}_original_rows.tsv';found=[];n=0
    if '--reuse-extracts' in sys.argv and dest.exists():
        df=pd.read_csv(dest,sep='\t')
        if set(df[idcol])!=set(SNPS) or len(df)!=2:raise AssertionError('Invalid recovery extract')
        rec=record(src,name+' original source')
        rec.update(found=list(df[idcol]),extract_sha256=sha(dest),genome_build='GRCh37' if name=='eQTLGen' else 'GRCh38',extraction='freshly extracted in first attempt on 2026-10-03; reused after a downstream CI-convention check halted execution',rows_scanned_until_both_found=None)
        print(name,'verified recovered extract',flush=True)
        return name,df,rec
    with gzip.open(src,'rb') as f,dest.open('wb') as out:
        head=f.readline();out.write(head);headers=head.decode().strip().split('\t');ix=headers.index(idcol)
        gi=headers.index('Gene') if name=='eQTLGen' else None
        for line in f:
            n+=1
            if b'rs179252' not in line and b'rs179247' not in line:continue
            cells=line.decode().rstrip('\r\n').split('\t')
            ids=cells[ix].split(',')
            if not any(s in ids for s in SNPS):continue
            if gi is not None and cells[gi]!='ENSG00000165409':continue
            found.append(cells[ix]);out.write(line)
            if len(found)==2:break
    rec=record(src,name+' original source')
    rec.update(rows_scanned_until_both_found=n,found=found,extract_sha256=sha(dest),genome_build='GRCh37' if name=='eQTLGen' else 'GRCh38',extraction='fresh original compressed source, stopped after the two target records; not a complete source scan')
    print(name,'fresh extraction:',found,'after',n,'rows',flush=True)
    if set(found)!=set(SNPS):raise AssertionError(rec)
    return name,pd.read_csv(dest,sep='\t'),rec

with ThreadPoolExecutor(max_workers=3) as pool: extracted=list(pool.map(extract,sources))
data={name:d for name,d,rec in extracted};manifests.extend(rec for _,_,rec in extracted)
(OUT/'SOURCE_MANIFEST.json').write_text(json.dumps(manifests,indent=2),encoding='utf-8')
eq=data['eQTLGen'];cached=pd.read_csv(SRC/'TrackA_MR/data/eqtl_TSHR_all.csv')
for snp in SNPS:
    a=eq[eq.SNP.eq(snp)].iloc[0];b=cached[cached.SNP.eq(snp)].iloc[0]
    for col in ['AssessedAllele','OtherAllele','Zscore','NrSamples','SNPPos']:
        check(f'{snp} source versus gene extract {col}',a[col]==b[col])
eq.to_csv(TABLE/'Table_A1_eQTL_original_direction.csv',index=False)

exposure=pd.read_csv(ROOT/'work/reanalysis_20261002/mr_v1/exposure_comparison.csv')
ex=exposure[exposure.snp.eq('rs179252')&exposure.gene_symbol.eq('TSHR')].iloc[0]
iv=eq[eq.SNP.eq('rs179252')].iloc[0]
frfile=ROOT/'work/maf_sensitivity_20260905/inputs/eur_freq_selected.tsv'
fr=pd.read_csv(frfile,sep='\t');fr=fr[fr.SNP.eq('rs179252')].iloc[0]
af=float(fr.MAF if fr.A1==iv.AssessedAllele else 1-fr.MAF)
sx=1/math.sqrt(2*af*(1-af)*(iv.NrSamples+iv.Zscore**2));bx=iv.Zscore*sx
check('exposure beta reconstructed',abs(bx-ex.beta_ref)<1e-12)
check('exposure SE reconstructed',abs(sx-ex.se_ref)<1e-12)
check('effect allele definition',iv.AssessedAllele==ex.effect_allele)
manifests.extend([record(frfile,'existing rounded reference MAF used in manuscript'),record(ROOT/'work/reanalysis_20261002/mr_v1/exposure_comparison.csv','canonical reconstructed exposure')])

canonical=pd.read_csv(ROOT/'outputs/TED_Reanalysis_20261002/Tables/All_2544_genes_five_outcomes.csv')
ocnames={'BBJ':'BBJ_Graves','UKB':'UKB_hyperthyroid','FinnGen_GO':'FinnGen_GO','E4_GRAVES_STRICT':'E4_GRAVES_STRICT','E4_GRAVES_OPHT_STRICT':'E4_GRAVES_OPHT_STRICT'}
gwasrows=[];mrrows=[]
for name,d in data.items():
    if name=='eQTLGen':continue
    for _,row in d.iterrows():
        if name in ['BBJ','UKB']:
            snp=row['rsid' if name=='BBJ' else 'hm_rsid'];ea=row.effect_allele;oa=row.other_allele;pos=int(row.base_pair_location);by=float(row.beta);sy=float(row.standard_error);p=float(row.p_value);eaf=float(row.effect_allele_frequency)
        else:snp=row.rsids;ea=row.alt;oa=row.ref;pos=int(row.pos);by=float(row.beta);sy=float(row.sebeta);p=float(row.pval);eaf=float(row.af_alt)
        original_beta=by;original_se=sy
        if name=='UKB':mu=3731/484598;by/=mu*(1-mu);sy/=mu*(1-mu)
        check(f'{name} {snp} G present',ea=='G' or oa=='G')
        sign=1 if ea=='G' else -1
        by*=sign;eaf_g=eaf if sign==1 else 1-eaf
        gwasrows.append(dict(outcome=name,snp=snp,pos38=pos,source_effect_allele=ea,source_other_allele=oa,source_beta=original_beta,source_se=original_se,source_p=p,source_effect_scale='linear probability' if name=='UKB' else 'log odds',aligned_allele='G',sign_multiplier=sign,aligned_G_beta_logOR=by,aligned_G_se=sy,aligned_G_OR=math.exp(by),aligned_G_CI_low=math.exp(by-norm.ppf(.975)*sy),aligned_G_CI_high=math.exp(by+norm.ppf(.975)*sy),aligned_G_eaf=eaf_g,recomputed_normal_p=2*norm.sf(abs(by/sy))))
        # Existing saved regional extract is an additional independent pipeline check.
        cache=pd.read_csv(ROOT/f'work/followup_validation_20261002/inputs/{name}_TSHR_raw_region.tsv',sep='\t')
        cache=cache[cache.snp.eq(snp)].iloc[0]
        check(f'{name} {snp} fresh GWAS versus regional cache',abs(sign*by-cache.beta)<1e-12 and abs(sy-cache.se)<1e-12 and ea==cache.ea and oa==cache.oa)
        if snp=='rs179252':
            # Earlier BBJ/UKB/GO results use 1.96; the new FinnGen-GD/strict-GO
            # pipeline uses qnorm(.975). Reproduce the declared original pipeline.
            ci_z=1.96 if name in ['BBJ','UKB','FinnGen_GO'] else norm.ppf(.975)
            beta=by/bx;se=sy/abs(bx);pv=2*norm.sf(abs(beta/se));OR=math.exp(beta);lo=math.exp(beta-ci_z*se);hi=math.exp(beta+ci_z*se)
            c=canonical[canonical.gene_symbol.eq('TSHR')&canonical.outcome.eq(ocnames[name])]
            check(f'{name} canonical row unique',len(c)==1,list(c.outcome))
            c=c.iloc[0]
            deltas={k:float(abs(v-c[k])) for k,v in [('beta',beta),('se',se),('OR',OR),('CI_lower',lo),('CI_upper',hi)]}
            check(f'{name} MR full precision reproduced',max(deltas.values())<1e-8,deltas)
            check(f'{name} MR P reproduced',abs(pv-c.pvalue)/c.pvalue<1e-6,{'new':pv,'canonical':c.pvalue})
            mrrows.append(dict(outcome=name,instrument=snp,exposure_allele=iv.AssessedAllele,beta_exposure=bx,se_exposure=sx,exposure_eaf_reference=af,exposure_F=iv.Zscore**2,beta_outcome_G=by,se_outcome=sy,MR_beta=beta,MR_se=se,MR_OR=OR,MR_CI_low=lo,MR_CI_high=hi,MR_p=pv,CI_z_multiplier=ci_z,canonical_OR=float(c.OR),max_numeric_difference=max(deltas.values()),source_gwas_p=p))

pd.DataFrame(gwasrows).to_csv(TABLE/'Table_A3_both_variants_five_disease_GWAS.csv',index=False)
pd.DataFrame(mrrows).to_csv(TABLE/'Table_A4_MR_reproduction.csv',index=False)
pd.DataFrame(checks).to_csv(OUT/'verification_checks.csv',index=False)
(OUT/'SOURCE_MANIFEST.json').write_text(json.dumps(manifests,indent=2),encoding='utf-8')
result={'finished_utc':datetime.now(timezone.utc).isoformat(),'checks':len(checks),'all_pass':all(x['pass_check'] for x in checks),'LD':ldrows,'MR':mrrows,'tissue_analysis_status':'not_started'}
(OUT/'AUDIT_RESULT.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print('COMPLETE',len(checks),'checks passed',flush=True)
