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
    guards = framing_guards()
    print(json.dumps({'status': 'PASS', 'manifest_files_verified': len(manifest['files']),
                      'framing_guards_passed': guards,
                      'tissue_estimates_recomputed': len(estimated), 'numeric_comparisons': comparisons,
                      'maximum_relative_difference': max_relative, 'read_only': True,
                      'scope': 'Delivery bytes, document containers, planned tissue grid and Wald formula implementation.',
                      'review': 'internal/REVIEW_codex_package_20261011.md items addressed on 2026-10-14; author items (journal requirements, all-author approval) remain open.'}, indent=2))


def framing_guards():
    """Wording guards; each encodes a defect that was actually removed (2026-10-14)."""
    import re
    text = MASTER.read_text(encoding='utf-8').replace('\r\n', '\n')
    main_text = text.split('## Declarations')[0]
    abstract = text.split('## Abstract')[1].split('**Keywords')[0]
    conclusions = text.split('### Conclusions')[1].split('## Declarations')[0]
    pkg = {n: (CANDIDATE / n).read_text(encoding='utf-8') for n in
           ['Supplementary_Material.md', 'Figure_Legends.md', 'Tables.md', 'Cover_Letter_Draft.md',
            'SupplementaryTables.md']}
    everything = text + '\n'.join(pkg.values())
    words = lambda s: len([w for w in re.sub(r'[*_`]', '', s).split() if w])
    sections = re.split(r'(?m)^## ', text)
    body = {s.split('\n', 1)[0].strip(): s.split('\n', 1)[1] for s in sections[1:]}
    main_words = sum(words('\n'.join(l.lstrip('#') for l in body[k].split('\n') if l.strip() and l.strip() != '---'))
                     for k in ['Introduction', 'Methods', 'Results', 'Discussion'])
    failures, checked = [], []

    def ok(cond, label):
        checked.append(label)
        if not cond:
            failures.append(label)

    ok(words(abstract) <= 250, 'Abstract <= 250 words')
    ok(main_words <= 5000, 'Main text incl. sub-headings <= 5,000 words')
    ok(not re.search(r'AI-assisted|artificial intelligence|language model|ChatGPT|OpenAI|Codex|Claude', text + pkg['Cover_Letter_Draft.md'] + pkg['Supplementary_Material.md']),
       'No AI-tool declaration (2026-09-05 author decision)')
    ok('0.226' in abstract and '0.264' in abstract, 'Abstract discloses TSHR non-colocalization in UKB and FinnGen GD')
    ok('prior-sensitiv' in abstract and 'strict GO' in abstract.split('prior-sensitiv')[1][:40],
       'Abstract flags prior sensitivity for BBJ and strict GO')
    ok('Mean orbital' in abstract, 'Abstract orbital sentence keeps "Mean" (IGF1R exceeds control in 3 of 4)')
    ok('colocalization' in conclusions.lower() and 'UKB' in conclusions and 'FinnGen GD' in conclusions,
       'Conclusions scope colocalization to the datasets where it failed')
    ok('does not argue against' in conclusions, 'Conclusions keep the IGF1R therapy guard')
    ok(not re.search(r'adipose extension supported|Adipose \*TSHR\*.GD associations passed correction', text),
       'Adipose result not presented as independent support')
    ok('re-express' in text and 'rs12101255' in text, 'Adipose instruments identified as the same regional signal')
    ok('consistently indicated higher odds' not in everything, 'IGF1R not described as consistently indicating higher odds')
    ok(not re.search(r'[Ss]econdary colocalization|secondary interpretive', everything), 'Colocalization not relabelled secondary')
    ok('strengthens the plausibility' not in text, 'Allele audit not overstated')
    ok('0.80 (BBJ), 0.77 (UKB), 0.78 (FinnGen ophthalmopathy), 0.77 (FinnGen GD) and 0.82 (strict ophthalmopathy)' in text,
       'TSHR per-allele ORs explain the per-SD scale')
    ok('5×10⁻⁵ and a PP.H4' not in text and '5×10^−5^ and a PP.H4' not in text, 'Unverifiable Ji et al. prior not asserted')
    ok(not re.search(r'case series|prespecified|pre-specified|first systematic|SuSiE|credible set|DESeq2|padj', everything),
       'Withdrawn or banned phrasing absent')
    ok(not re.search(r'TED-specific (outcome|dataset|cohort|endpoint)', everything), 'FinnGen outcome not labelled TED-specific')
    ok(not re.search(r'(?i)(should be confirmed by|before submission|previously reviewed|unchanged from the earlier|subject to author confirmation|This version is prepared for author review)', everything),
       'No internal to-do or version-history notes in reviewer-facing text')
    ok(not re.search(r'(?<![*\w])P (?:=|<) ', main_text), 'P italic in the main text')
    ok('| PP.H2 |' in pkg['Tables.md'], 'Table 2 shows PP.H2')
    ok('| Q-test *P* |' in pkg['Tables.md'], 'Table 2 shows the Q-test P (2026-10-14 display review)')
    ok(not re.search(r'Table 3[ABC]\b', everything), 'No Table 3A-C after the display review')
    first = {}
    for kind, nums in re.findall(r'\b(Tables?|Figures?) (S?\d(?:(?:, | and |–)S?\d)*)', main_text):
        for n in re.split(r', | and |–', nums):
            seen = first.setdefault((kind.rstrip('s'), n.startswith('S')), [])
            if n not in seen:
                seen.append(n)
    ok(first.get(('Table', False)) == ['1', '2', '3'] and first.get(('Table', True)) == ['S1', 'S2', 'S3', 'S4']
       and first.get(('Figure', False)) == ['1', '2', '3'] and first.get(('Figure', True)) == ['S1', 'S2'],
       'Tables and figures numbered in order of first citation')
    ok(re.findall(r'(?m)^## (Table S\d)', pkg['SupplementaryTables.md']) == ['Table S1', 'Table S2', 'Table S3', 'Table S4'],
       'Supplementary tables S1-S4 in one source file')
    refs = re.findall(r'(?m)^(\d+)\. ', body['References'])
    cited = {int(n) for grp in re.findall(r'\[(\d+(?:[,–-]\s*\d+)*)\]', main_text) for n in re.split(r'[,–-]\s*', grp)}
    ok(len(refs) == 35 and cited == set(range(1, 36)), '35 references, all cited (33 + Zuber 2022, Douglas 2024; 18 October)')
    with ZipFile(CANDIDATE / 'Manuscript.docx') as z:
        doc = re.sub(r'<[^>]+>', '', z.read('word/document.xml').decode('utf-8'))
    ok(all(f'Figure {i}.' in doc for i in (1, 2, 3)), 'Figure 1-3 legends present in Manuscript.docx')
    ok('AI-assisted' not in doc, 'No AI wording in Manuscript.docx')
    derived = CANDIDATE / 'Verification/Review_20261014/derived_values.csv'
    with derived.open(encoding='utf-8', newline='') as f:
        for r in csv.DictReader(f):
            shown = r['displayed'].replace('%', '')
            ok(shown in everything, 'Derived value displayed: ' + r['item'])
    assert not failures, 'Framing guard failures: ' + '; '.join(failures)
    return len(checked)


if __name__ == '__main__':
    main()
