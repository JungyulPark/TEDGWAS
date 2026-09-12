"""Export aggregate alternative-estimator results with explicit uncertainty scales.

The original estimates and P values are preserved. Student-t intervals for
Egger/mode use the same degrees of freedom as their original P-value tests.
No individual-level or raw regional summary statistics are exported.
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import stats

D = Path(__file__).resolve().parents[1]
P = D / 'provenance'
raw = pd.read_csv(P / 'MR_all_estimators_verified.csv')
out = raw.drop(columns=['or_ci_lower', 'or_ci_upper']).copy()
out['test_distribution'] = 'normal'
out['test_df'] = np.nan
out['ci_critical_value'] = stats.norm.ppf(.975)
for method, subtract in [('MR Egger', 2), ('Weighted mode', 1)]:
    ix = out.method == method
    out.loc[ix, 'test_distribution'] = 'Student t'
    out.loc[ix, 'test_df'] = out.loc[ix, 'n_iv'] - subtract
    out.loc[ix, 'ci_critical_value'] = stats.t.ppf(.975, out.loc[ix, 'test_df'])
out['log_or_ci_lower'] = out.beta - out.ci_critical_value * out.se
out['log_or_ci_upper'] = out.beta + out.ci_critical_value * out.se
with np.errstate(over='ignore', under='ignore'):
    out['or_ci_lower'] = np.exp(out.log_or_ci_lower)
    out['or_ci_upper'] = np.exp(out.log_or_ci_upper)
out['frequency_scenario'] = 'original_reference'
out.to_csv(D / 'Supplementary_Data_5_MR_estimators.csv', index=False, na_rep='NA')
checks = []
for method, group in out.groupby('method'):
    test = abs(group.beta / group.se)
    calc = 2 * (stats.t.sf(test, group.test_df) if method in ['MR Egger', 'Weighted mode'] else stats.norm.sf(test))
    assert np.allclose(calc, group.pvalue, rtol=1e-6, atol=1e-300), method
    checks.append({'method': method, 'rows': len(group), 'P_values_match_original_test': True})
assert not out.duplicated(['gene_symbol', 'outcome', 'method']).any()
result = {'status': 'PASS', 'rows': len(out), 'methods': checks,
          'original_estimates_unchanged': True,
          'CI_definition': 'Two-sided 95% intervals use the original test distribution: normal for Wald/IVW/weighted median; Student t with n_iv-2 for MR Egger and n_iv-1 for weighted mode. Log-scale interval endpoints are retained if exponentiation overflows.',
          'sources': ['https://github.com/MRCIEU/TwoSampleMR/blob/master/R/mr.R', 'https://github.com/MRCIEU/TwoSampleMR/blob/master/R/mr_mode.R'],
          'scope': 'Export and arithmetic verification of preserved analyses; this does not rerun estimation, directionality, or multi-signal colocalization.'}
(P / 'complete_statistics_export.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf8')
print(json.dumps(result, indent=2))
