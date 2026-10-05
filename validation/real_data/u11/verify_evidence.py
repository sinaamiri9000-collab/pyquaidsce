"""Read-only integrity/consistency checks of published U11 aggregate evidence."""
from pathlib import Path
import hashlib
import json

import numpy as np
import pandas as pd


def main():
    root = Path(__file__).resolve().parent
    hashes = json.loads((root / 'FILE_HASHES.json').read_text())
    for name, expected in hashes.items():
        path = root / name
        assert path.is_file(), name
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected, name
    result_dir = root / 'results'
    records = {name: json.loads((result_dir / f'{name}.json').read_text())
               for name in ('standard', 'nlsur_criteria')}
    checks = json.loads((result_dir / 'checks.json').read_text())
    assert len(checks) == 30 and all(c['passed'] for c in checks)
    audit = json.loads((result_dir / 'unchanged_engine_audit.json').read_text())
    assert audit['inner_optimization_AST_identical']
    assert audit['outer_nonconvergence_AST_identical']
    history = {}
    for name, record in records.items():
        assert record['converged'] and record['no_curvature_constraint']
        assert record['nobs'] == 37658 and record['n_free'] == 198
        assert sum(r['iterations'] for r in record['inner_history']) == record['n_gn']
        assert len(record['inner_history']) == record['n_outer']
        assert record['n_ifgnls_updates'] == record['n_outer'] - 2
        with np.load(result_dir / f'{name}_parameter_history.npz', allow_pickle=False) as data:
            assert set(data.files) == {'accepted_thetas', 'stage_thetas', 'stage_sigmas'}
            assert data['accepted_thetas'].shape == (record['n_gn'], 198)
            assert data['stage_thetas'].shape == (record['n_outer'], 198)
            assert data['stage_sigmas'].shape == (record['n_outer'], 11, 11)
            assert np.array_equal(data['stage_thetas'][-1], record['theta'])
            assert np.array_equal(data['stage_sigmas'][-1], record['sigma'])
            history[name] = data['accepted_thetas'].copy()
        covariance = pd.read_csv(result_dir / f'{name}_parameter_covariance.csv',
                                 index_col=0, float_precision='round_trip').to_numpy()
        assert covariance.shape == (198, 198) and np.isfinite(covariance).all()
    assert np.array_equal(history['standard'][:11], history['nlsur_criteria'][:11])
    difference = json.loads((result_dir / 'differences.json').read_text())
    delta = np.asarray(records['nlsur_criteria']['theta']) - records['standard']['theta']
    assert np.max(np.abs(delta)) == difference['max_abs_free_coefficient_difference']
    assert (records['nlsur_criteria']['llf'] - records['standard']['llf']) == difference['likelihood_difference']
    print(f'PASS: {len(hashes)} file hashes; 30 recorded checks; source audits; '
          'fit counts, covariance shapes and bitwise-equal initial NLS trajectory.')


if __name__ == '__main__':
    main()
