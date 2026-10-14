#!/usr/bin/env python3
"""Read-only October handoff integrity and tissue-formula verification.

The September entry point is preserved in internal/handoff_before_20261014/.
PASS does not resolve the 13 October review or certify submission readiness.
Python 3.9+ standard library only. Does not change files or line endings.
"""
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parents[1]
CANDIDATE = ROOT / 'submission/candidate_20261011_review'
MASTER = ROOT / 'MANUSCRIPT_TED_TRAP_v5_MASTER.md'


def main():
    manifest = json.loads((CANDIDATE / 'GIT_MANIFEST.json').read_text(encoding='utf-8'))
    for entry in manifest['files']:
        path = (ROOT / entry['path']).resolve()
        assert path.is_relative_to(ROOT.resolve()), entry['path']
        assert path.is_file(), 'Missing file: ' + entry['path']
        assert hashlib.sha256(path.read_bytes()).hexdigest() == entry['sha256'], 'Changed file: ' + entry['path']
    candidates = sorted(p.name for p in (ROOT / 'submission').glob('candidate_*') if p.is_dir())
    assert candidates == [CANDIDATE.name], candidates
    assert not (CANDIDATE / 'Manuscript.md').exists(), 'Keep one editable Markdown master'
    assert MASTER.is_file()
    expected = {entry['path'] for entry in manifest['files']}
    actual = {p.relative_to(ROOT).as_posix() for p in CANDIDATE.rglob('*') if p.is_file() and p.name != 'GIT_MANIFEST.json'}
    assert actual <= expected, 'Unmanifested files: ' + str(actual - expected)
    for p in CANDIDATE.rglob('*.json'):
        json.loads(p.read_text(encoding='utf-8'))
    documents = ['Manuscript.docx', 'Tables.docx', 'Supplementary_Material.docx', 'STROBE_MR_Checklist.docx', 'Cover_Letter_Draft.docx', 'Supplementary_Dataset.xlsx']
    for name in documents:
        with ZipFile(CANDIDATE / name) as archive:
            assert archive.testzip() is None, name
    for name in ['Manuscript.pdf', 'Full_Review.pdf', 'Final_Tables_Figures.pdf']:
        assert (CANDIDATE / name).read_bytes().startswith(b'%PDF-'), name

    def rows(name):
        p = CANDIDATE / 'Verification/Tissue_extension/Tables' / name
        with p.open(encoding='utf-8', newline='') as f:
            return list(csv.DictReader(f))

    planned = rows('Tissue_MR_50_planned_combinations.csv')
    estimated = [r for r in planned if r['status'] == 'estimated']
    assert len(planned) == 50 and len(estimated) == 20
    assert sum(r['status'] == 'not_estimable' for r in planned) == 30
    assert sum(r['family'] == 'primary' for r in planned) == 10
    assert sum(r['family'] == 'secondary' for r in planned) == 40
    key = lambda r: (r['gene'], r['tissue'], r['outcome'])
    harmonized_rows = rows('Harmonized_instruments_all_outcomes.csv')
    harmonized = {key(r): r for r in harmonized_rows}
    assert len(harmonized) == len(harmonized_rows) == 20
    zcrit = NormalDist().inv_cdf(0.975)
    comparisons, max_relative = 0, 0.0
    for r in estimated:
        assert int(r['n_retained_instruments']) == 1
        h = harmonized[key(r)]
        beta = float(h['beta_outcome']) / float(h['beta_exposure'])
        se = float(h['se_outcome']) / abs(float(h['beta_exposure']))
        p = math.erfc(abs(beta / se) / math.sqrt(2))
        n = 10 if r['family'] == 'primary' else 40
        assert int(r['planned_family_tests']) == n
        values = {'beta': beta, 'se': se, 'OR': math.exp(beta),
                  'CI_lower': math.exp(beta - zcrit * se),
                  'CI_upper': math.exp(beta + zcrit * se),
                  'P': p, 'P_adjusted': min(1.0, n * p)}
        for name, calculated in values.items():
            recorded = float(r[name])
            assert math.isclose(calculated, recorded, rel_tol=1e-10, abs_tol=0), (key(r), name)
            max_relative = max(max_relative, abs(calculated - recorded) / max(abs(recorded), 1e-300))
            comparisons += 1
    assert sum(r['gene'] == 'IGF1R' and float(r['P']) > 0.05 for r in estimated) == 10
    categories = [r['statistical_category'] for r in estimated if r['gene'] == 'TSHR']
    assert categories.count('family_corrected') == 8
    assert categories.count('nominal_only') == 1
    assert categories.count('not_nominally_significant') == 1
    for p in CANDIDATE.rglob('*'):
        if p.is_file():
            assert p.suffix.lower() not in {'.raw', '.bed', '.bim', '.fam', '.gz', '.log'}, p
            assert p.name != 'data.txt' and 'Source_extracts' not in p.parts, p
            assert not ('Sources' in p.parts and p.suffix == '.tsv'), p
    print(json.dumps({'status': 'PASS', 'manifest_files_verified': len(manifest['files']),
                      'tissue_estimates_recomputed': len(estimated), 'numeric_comparisons': comparisons,
                      'maximum_relative_difference': max_relative, 'read_only': True,
                      'scope': 'Delivery bytes, document containers, planned tissue grid and Wald formula implementation.',
                      'open_review': 'internal/REVIEW_codex_package_20261011.md; source-dependent regional reproduction and scientific/document review remain open.'}, indent=2))


if __name__ == '__main__':
    main()
