from pathlib import Path
import hashlib
import json
import shutil
import zipfile

W = Path(__file__).resolve().parent
R = W.parents[1]
O = R / 'outputs/TED_Final_Review_20261011'

def sha256(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

source_checks = json.loads((W / 'source_verification.json').read_text(encoding='utf-8'))
assert all(x['all_sha256_match'] for x in source_checks)
pdf_checks = json.loads((O / 'Verification/PDF_assembly_checks.json').read_text(encoding='utf-8'))
assert pdf_checks['status'] == 'PASS' and pdf_checks['pages'] == 13
qa = {
    'delivery_date': '2026-10-11',
    'scientific_content_version': '2026-10-03',
    'new_statistical_analysis_performed_for_this_delivery': False,
    'source_checks': source_checks,
    'combined_pdf': pdf_checks,
    'visual_review': {
        'method': 'Rendered all 13 PDF pages with Poppler and inspected all four contact sheets.',
        'pages_inspected': list(range(1, 14)),
        'clipping_or_overlaps_observed': False,
        'source_page_text_preserved': True,
    },
    'status': 'PASS: delivery assembly checks; author final review remains required.',
}
(O / 'Verification/Delivery_QA_20261011.json').write_text(json.dumps(qa, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')
scripts = O / 'Verification/Delivery_scripts'
scripts.mkdir(exist_ok=True)
for script in W.glob('*.py'):
    shutil.copy2(script, scripts / script.name)

entries = []
for path in sorted(O.rglob('*')):
    if path.is_file() and path != O / 'MANIFEST_SHA256.json':
        entries.append({'path': path.relative_to(O).as_posix(), 'bytes': path.stat().st_size, 'sha256': sha256(path)})
manifest = {
    'delivery_date': '2026-10-11',
    'scientific_content_version': '2026-10-03',
    'scope': 'Every file in this delivery except this manifest itself.',
    'files': entries,
}
(O / 'MANIFEST_SHA256.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')

archive = O.with_suffix('.zip')
assert not archive.exists(), 'Preserve existing delivery archives'
with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for path in sorted(O.rglob('*')):
        if path.is_file():
            z.write(path, (Path(O.name) / path.relative_to(O)).as_posix())
with zipfile.ZipFile(archive) as z:
    assert z.testzip() is None
    for e in entries:
        content = z.read(O.name + '/' + e['path'])
        assert len(content) == e['bytes']
        assert hashlib.sha256(content).hexdigest() == e['sha256']
result = {'zip': str(archive), 'bytes': archive.stat().st_size, 'files': len(entries) + 1, 'crc_and_sha256': 'PASS', 'zip_sha256': sha256(archive)}
(W / 'zip_verification.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
