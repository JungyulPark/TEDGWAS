from pathlib import Path
import json,re,shutil,hashlib,zipfile,datetime,subprocess
import pandas as pd
T=Path(__file__).resolve().parent;W=T.parent;O=W.parent/'outputs/TED_Reanalysis_20261002';V=O/'Verification'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
page=(O/'START_HERE.html').read_text(encoding='utf8')
payload=re.search(r'<script id="data" type="application/json">(.*?)</script>',page,re.S).group(1)
d=json.loads(payload);assert len(d['all'])==12720 and len(d['focus'])==70 and len(d['coloc'])==840
lookup={(r['gene'],r['outcome'],r['scenario'],r['p12']):r for r in d['coloc']}
assert len(lookup)==840
for r in d['focus']:
    c=lookup[r['gene_symbol'],r['outcome'],'original_reference',1e-5]
    assert max(abs(r[f'PP.H{i}']-c[f'PP.H{i}']) for i in range(5))<1e-12
    for sc in ['original_reference','available_eqtlgen']:
        for p in [1e-5,5e-6,1e-6]:assert (r['gene_symbol'],r['outcome'],sc,p) in lookup
for href in re.findall(r'(?:href|src)="([^"]+)"',page):
    if '://' not in href:assert (O/href).is_file(),href
script=page.split('</script><script>')[1].split('</script>')[0]
(T/'dashboard_syntax.js').write_text(script,encoding='utf8')
node=Path('${LOCAL_USER_HOME}/.cache/codex-runtimes/codex-primary-runtime/dependencies/node/bin/node.exe')
subprocess.run([str(node),'--check',str(T/'dashboard_syntax.js')],check=True,capture_output=True,text=True)
mr=pd.read_csv(T/'mr_v1/MR_all_genes_all_scenarios.csv')
focus=mr[mr.outcome.eq('E4_GRAVES_STRICT')&mr.BBJ_discovery]
a=set(focus[focus.scenario.eq('original_reference')&focus.focused_GD_corroboration.fillna(False)].gene_symbol)
b=set(focus[focus.scenario.eq('reharmonized_eqtlgen')&focus.focused_GD_corroboration.fillna(False)].gene_symbol)
assert a==b and len(a)==7
shutil.copy2(T/'09_prior_plot_v3.py',V/'scripts/09_prior_plot_v3.py')
shutil.copy2(Path(__file__),V/'scripts'/Path(__file__).name)
(V/'display_verification.json').write_text(json.dumps({'static_payload_rows':{'all':12720,'focus':70,'coloc':840},'MR_coloc_payload_agreement':True,'all_relative_links_exist':True,'JavaScript_syntax_check':'PASS','images':'R1/R2/R4 visually inspected; R3 legend repositioned outside plot','browser_interaction':'NOT VERIFIED: browser security policy rejects local file protocol; no workaround attempted','seven_gene_confirmation_stable_after_frequency_substitution':True},indent=2),encoding='utf8')
(O/'README_KO.txt').write_text('''TED–GD 다섯 질환 자료 재분석 · 2026-10-02

먼저 REPORT_KO.md를 읽으십시오. START_HERE.html은 유전자·prior·빈도 조건을 선택하는 검토 화면입니다. 압축을 풀고 열면 같은 폴더의 Figures를 표시합니다. 브라우저 도구 정책이 file: 열기를 차단하여 실제 클릭 동작은 자동 검증하지 못했습니다. 내장 데이터, 연결 파일, JavaScript 문법은 확인했습니다.

Tables: 전체 2,544개×5개 자료, 주요 14개 통합표, 새 두 자료의 모든 빈도 민감도·leave-one-out·추적 표.
Figures: 검토용 R1–R4의 PNG와 PDF. 기존 원고의 Figure를 교체하지 않았습니다.
Verification: 분석 계획, 원자료·선택파일 해시, 검산, 실행 기록, 공식 표본 수, 코드.
DATA_DICTIONARY.md: 열 이름과 통계 수치의 의미.

핵심: 새 GD에서 TSHR MR 연관성은 강하지만 H4=0.264로 공유변이 지지가 일관되게 재현되지 않았습니다. 기존 결론을 그대로 강화한 결과로 쓰면 안 됩니다.

이 묶음은 분석 검토본입니다. 기존 제출 패키지의 교체본 또는 제출 준비 완료 확인서가 아닙니다.

재현 위치: work/reanalysis_20261002의 버전 스크립트와 기존 로컬 입력을 사용합니다. 원본 GWAS와 이전 입력 파일은 용량 때문에 이 압축에 포함하지 않았습니다. 정확한 경로·해시를 Verification/input_hashes.json 및 extraction_manifest.json에 기록했습니다. 분석 스크립트를 일괄 실행하지 말고 각 완료 단계와 버전별 출력 존재 여부를 확인하십시오.
''',encoding='utf8')
with (T/'analysis.log').open('a',encoding='utf8') as fh:
    fh.write(datetime.datetime.now(datetime.timezone.utc).isoformat()+' Completed: MR primary 5088 rows; sensitivity 20352 rows; new coloc 336 settings; LOO 42 rows; identity 2480 MR / 168 coloc; R/Python validated. Diagnostic plots and delivery assembled. Old manuscript unchanged.\n')
shutil.copy2(T/'analysis.log',V/'analysis.log')
manifest=[{'path':str(p.relative_to(O)).replace('\\','/'),'bytes':p.stat().st_size,'sha256':sha(p)} for p in sorted(O.rglob('*')) if p.is_file() and p.name!='PACKAGE_MANIFEST.json']
(O/'PACKAGE_MANIFEST.json').write_text(json.dumps(manifest,indent=2),encoding='utf8')
dest=Path('${LOCAL_USER_HOME}/Documents/Codex/TED_Reanalysis_20261002');assert not dest.exists(),'Do not overwrite delivery'
shutil.copytree(O,dest)
z=dest.with_suffix('.zip');assert not z.exists()
with zipfile.ZipFile(z,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as f:
    for p in sorted(dest.rglob('*')):
        if p.is_file():f.write(p,arcname=str(p.relative_to(dest.parent)))
with zipfile.ZipFile(z) as f:assert f.testzip() is None
for x in manifest:assert sha(dest/x['path'])==x['sha256']
print(json.dumps({'delivery':str(dest),'zip':str(z),'files':len(manifest)+1,'zip_bytes':z.stat().st_size,'zip_sha256':sha(z),'verification':'PASS numeric/static checks; browser interactions not tested'},indent=2))
