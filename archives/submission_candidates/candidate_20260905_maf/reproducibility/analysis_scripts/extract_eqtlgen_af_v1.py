"""Extract the original eQTLGen cohort AF, verify counts and full gzip CRC."""
from pathlib import Path
import csv,gzip,json,hashlib,datetime
W=Path(__file__).resolve().parent;T=W/'maf_sensitivity_20260905'
source=W/'data/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz'
dest=T/'inputs/eqtlgen_af_selected_v1.csv'
if dest.exists():raise RuntimeError('Versioned AF extract already exists')
coords={}
def register(s,c,p):
    v=(str(c),str(p))
    if s in coords and coords[s]!=v:raise ValueError(f'Inconsistent source hg19 coordinates: {s}')
    coords[s]=v
with (T/'inputs/instruments_verified.csv').open(newline='') as f:
    for r in csv.DictReader(f):register(r['snp'],r['chr'],r['pos_hg19'])
with (T/'inputs/eqtlgen_9genes.tsv').open(newline='') as f:
    for r in csv.DictReader(f,delimiter='\t'):register(r['SNP'],r['SNPChr'],r['SNPPos'])
rows=[];seen=set();total=0;maxerr=0;badcoords=[];invalid=[]
with gzip.open(source,'rt',newline='') as f:
    reader=csv.reader(f,delimiter='\t');header=next(reader)
    expected=['SNP','hg19_chr','hg19_pos','AlleleA','AlleleB','allA_total','allAB_total','allB_total','AlleleB_all']
    assert header==expected,header
    for row in reader:
        total+=1
        s=row[0]
        if s not in coords:continue
        if s in seen:raise ValueError(f'Duplicate requested rsID in official AF: {s}')
        seen.add(s)
        if tuple(row[1:3])!=coords[s]:badcoords.append(dict(snp=s,expected=coords[s],observed=row[1:3]));continue
        aa,ab,bb=map(float,row[5:8]);n=aa+ab+bb;af=float(row[8])
        assert min(aa,ab,bb)>=0 and n>0
        err=abs(af-(ab+2*bb)/(2*n));maxerr=max(maxerr,err)
        assert err<1e-12,(s,err)
        if not 0<af<1:invalid.append(s);continue
        rows.append(dict(snp=s,chr=row[1],pos=row[2],allele_a=row[3],allele_b=row[4],af_b=af,n_genotypes=n,genome_build='GRCh37'))
fields=['snp','chr','pos','allele_a','allele_b','af_b','n_genotypes','genome_build']
with dest.open('w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
manifest=dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_filename=source.name,source_bytes=source.stat().st_size,source_sha256=hashlib.file_digest(source.open('rb'),'sha256').hexdigest(),source_rows=total,requested_unique_snps=len(coords),requested_found=len(seen),usable_selected=len(rows),absent_rsids=sorted(set(coords)-seen),coordinate_mismatches=badcoords,invalid_frequency_rsids=invalid,max_selected_AF_count_error=maxerr,full_gzip_crc='PASS (read to EOF)',selected_sha256=hashlib.file_digest(dest.open('rb'),'sha256').hexdigest(),official_source='https://www.eqtlgen.org/cis-eqtls.html',frequency_cohort='26,609 eQTLGen samples; FHS excluded by provider',effect_frequency_column='AlleleB_all',authorization_manifest='work/data/af_download_manifest_v2.json',local_only=True)
(T/'af_input_verification_v1.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps({k:v for k,v in manifest.items() if k not in ('absent_rsids','coordinate_mismatches','invalid_frequency_rsids')},indent=2))
print('Absent',len(manifest['absent_rsids']),'coordinate mismatches',len(badcoords),'invalid AF',len(invalid))
if badcoords:raise ValueError('Coordinate mismatches require review before analysis')
