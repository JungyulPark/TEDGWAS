from pathlib import Path
import pandas as pd,json,hashlib,zipfile,shutil,re,datetime
T=Path(__file__).resolve().parent;O=T.parent.parent/'outputs/TED_Followup_Validation_20261002'
D=Path('${LOCAL_USER_HOME}/Documents/Codex/TED_Followup_Validation_20261002')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def pp(p):return f'{p:.3f}' if p>=.001 else f'{p:.2e}'
def pval(p):return f'{p:.3g}' if p>=.001 else f'{p:.2e}'
def main():
    assert not D.exists() and not D.with_suffix('.zip').exists()
    df=pd.read_csv(O/'Tables/Table_F1_24_genes_five_outcomes.csv');assert len(df)==120 and df.gene.nunique()==24
    genes=['TSHR','IGF1R','CTLA4','HLA-A','HLA-DQA2','C4A','TUBB','PSMB8','HSD3B7','VKORC1','TNFSF14','PRSS36','MAPKAPK5','IFNGR1','APOM','NEU1','C4B','HIST1H4E','LTA','LTB','TNF','HLA-G','IGLV2-18','IGLV3-12']
    labels={'BBJ_Graves':'BBJ GD','E4_GRAVES_STRICT':'FinnGen GD strict','UKB_hyperthyroid':'UKB hyperthyroidism','FinnGen_GO':'FinnGen GO','E4_GRAVES_OPHT_STRICT':'FinnGen GO strict'}
    lines=['# 24개 유전자 통합표 — 읽기용','', '기본 분석: 1000 Genomes EUR 빈도, p12=10⁻⁵. n IV는 실제 MR 도구 SNP 수다. PP.H3/H4는 공위치화 posterior이며 P값과 다르다.','', '기존 14개와 신규 탐색 후보 10개를 함께 제시했다. HLA-G는 빈도 민감도에서만 발견된 후보이다. MR 추정 불가를 음성 결과로 해석하지 않는다. FinnGen 결과들은 서로 독립적이지 않다.','']
    for oc,label in labels.items():
        lines += [f'## {label}','', '| Gene | IV 수 | MR OR (95% CI) | MR P | PP.H3 | PP.H4 |','|---|---:|---|---:|---:|---:|']
        a=df[df.outcome.eq(oc)].set_index('gene').loc[genes]
        for g,r in a.iterrows():
            fit=pd.notna(r.pvalue);orci=f'{r.OR:.3f} ({r.CI_lower:.3f}–{r.CI_upper:.3f})' if fit else '추정 불가'
            ps=pval(r.pvalue) if fit else '—';iv=str(int(r.n_iv)) if pd.notna(r.n_iv) else '0'
            note=' †' if oc=='BBJ_Graves' and g.startswith('IGLV') else ''
            lines.append(f'| {g}{note} | {iv} | {orci} | {ps} | {pp(r["PP.H3"])} | {pp(r["PP.H4"])} |')
        if oc=='BBJ_Graves':lines += ['', '† IGLV 두 유전자는 핵심 eQTL 부재 및 약한 주변 신호로 공위치화 해석이 제한된다. 계산된 H4는 표시하되 강한 음성 근거로 해석하지 않는다.']
        lines+=['']
    (O/'Tables/Table_F1_review.md').write_text('\n'.join(lines),encoding='utf8')
    report=O/'REPORT_KO.md';text=report.read_text(encoding='utf8');text=text.replace('- [전체 24개 × 5자료 표]', '- [24개 × 5자료 읽기용 통합표](Tables/Table_F1_review.md)\n- [전체 24개 × 5자료 표]');report.write_text(text,encoding='utf8')
    shutil.copy2(Path(__file__),O/'Verification/Scripts'/Path(__file__).name)
    checks=[]
    for f in O.rglob('*.md'):
        for link in re.findall(r'!?\[[^\]]*\]\(([^)]+)\)',f.read_text(encoding='utf8')):
            if re.match(r'https?://',link):continue
            assert (f.parent/link).resolve().is_file(),(f,link)
            checks.append(str((f.parent/link).relative_to(O)))
    for name,n in [('All_24_genes_1440_coloc_settings.csv',1440),('Table_F2_TSHR_windows_priors_frequencies.csv',240),('New_candidates_600_coloc_settings.csv',600),('C4B_HLAG_leave_one_out_72.csv',72)]:assert len(pd.read_csv(O/'Tables'/name))==n
    assert len(list((O/'Figures').glob('*.png')))==4 and len(list((O/'Figures').glob('*.pdf')))==4
    qc=dict(status='PASS',date='2026-10-02',readable_gene_endpoint_rows=120,relative_file_links_checked=len(checks),figures_visually_inspected=4,pdf_figures=4,existing_manuscript_modified=False,TSHR_analysis_version='v2',note='No journal-submission-ready claim. R and independent numeric checks are under Verification.')
    (O/'Verification/package_QC.json').write_text(json.dumps(qc,indent=2),encoding='utf8')
    files=[dict(path=str(p.relative_to(O)).replace('\\','/'),bytes=p.stat().st_size,sha256=sha(p)) for p in sorted(O.rglob('*')) if p.is_file()]
    (O/'PACKAGE_SHA256.json').write_text(json.dumps(dict(created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),files=files),indent=2),encoding='utf8')
    shutil.copytree(O,D)
    for f in files:assert sha(D/f['path'])==f['sha256']
    zpath=D.with_suffix('.zip')
    with zipfile.ZipFile(zpath,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in sorted(D.rglob('*')):
            if p.is_file():z.write(p,arcname=str(Path(D.name)/p.relative_to(D)))
    with zipfile.ZipFile(zpath) as z:assert z.testzip() is None
    result=dict(delivery=str(D),zip=str(zpath),bytes=zpath.stat().st_size,sha256=sha(zpath),files=len(files)+1)
    (T/'delivery.json').write_text(json.dumps(result,indent=2),encoding='utf8');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
