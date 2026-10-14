from pathlib import Path
import gzip,json,time,hashlib,datetime,traceback
import pandas as pd
T=Path(__file__).resolve().parent;W=T.parent;I=T/'inputs';I.mkdir(exist_ok=True)
def log(s):
    text=datetime.datetime.now(datetime.timezone.utc).isoformat()+' '+s
    print(text,flush=True)
    with (T/'run.log').open('a',encoding='utf8') as f:f.write(text+'\n')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    assert not (T/'plan_freeze.json').exists()
    (T/'plan_freeze.json').write_text(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'sha256':sha(T/'PLAN.md')},indent=2),encoding='utf8')
    iv=pd.read_csv(W/'maf_sensitivity_20260905/inputs/instruments_verified.csv')
    genes=['TSHR','APOM','NEU1','IGLV2-18','IGLV3-12','TNF','HIST1H4E','LTA','LTB','C4B','HLA-G']
    ids=iv[iv.gene_symbol.isin(genes)][['gene_symbol','ensembl_gene_id','chr']].drop_duplicates();assert len(ids)==11
    ids.to_csv(I/'gene_identity.csv',index=False);keys={x.encode() for x in ids.ensembl_gene_id}
    src=Path('C:/ProjectTEDGWAS/TrackA_MR/data/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz')
    out=I/'eqtl_11genes.tsv';n=k=0;start=time.time()
    with gzip.open(src,'rb') as f,out.open('wb') as dest:
        head=f.readline();dest.write(head);h=head.decode().strip().split('\t');ix=h.index('Gene')
        for line in f:
            n+=1;c=line.split(b'\t')
            if c[ix] in keys:dest.write(line);k+=1
            if n%10000000==0:log(f'eQTL: scanned {n:,}; selected {k:,}')
    e=pd.read_csv(out,sep='\t');assert set(e.Gene)==set(ids.ensembl_gene_id)
    assert not e.duplicated(['Gene','SNP']).any()
    qc=e.groupby(['GeneSymbol','Gene']).agg(rows=('SNP','size'),gene_chr=('GeneChr','first'),gene_pos_GRCh37=('GenePos','first'),min_pos_GRCh37=('SNPPos','min'),max_pos_GRCh37=('SNPPos','max'),min_p=('Pvalue','min'),max_p=('Pvalue','max'),min_n=('NrSamples','min'),max_n=('NrSamples','max')).reset_index()
    qc['genome_build']='GRCh37';assert qc.max_p.gt(.9).all();qc.to_csv(I/'eqtl_coverage.csv',index=False)
    old=pd.read_csv(W/'maf_sensitivity_20260905/inputs/eqtlgen_9genes.tsv',sep='\t');old=old[old.GeneSymbol.eq('TSHR')]
    fresh=e[e.GeneSymbol.eq('TSHR')];a=fresh.merge(old,on='SNP',suffixes=('_fresh','_old'),validate='one_to_one')
    assert len(a)==len(old)==len(fresh)==7309
    for col in ['SNPChr','SNPPos','AssessedAllele','OtherAllele','Zscore','NrSamples']:
        assert a[col+'_fresh'].eq(a[col+'_old']).all(),col
    # P strings may be rounded independently; signed Z and sample counts define beta reconstruction.
    windows={'source_GenePos_GRCh37':int(fresh.GenePos.iloc[0]),'full_source_cis':True,'half_widths_bp':[250000,500000], 'GWAS_audit_chr14_GRCh38':[80000000,82100000]}
    (T/'frozen_windows.json').write_text(json.dumps(windows,indent=2),encoding='utf8')
    record={'source':str(src),'source_bytes':src.stat().st_size,'source_sha256':sha(src),'rows_scanned':n,'rows_selected':k,'output_sha256':sha(out),'seconds':time.time()-start,'cached_TSHR_identity_rows':len(a)}
    (I/'eqtl_manifest.json').write_text(json.dumps(record,indent=2),encoding='utf8');log(json.dumps(record));print(qc.to_string(index=False))
if __name__=='__main__':
    try:main()
    except Exception:log('FAILED '+traceback.format_exc());raise
