"""Reproduce the ANSE small4 checks without changing point-estimation options."""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
import pandas as pd
import scipy

from pyquaidsce import quaidsce
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES


CASES = {
    'quadratic-ifgnls': {},
    'quadratic-nls': dict(method='nls'),
    'quadratic-fgnls': dict(method='fgnls'),
    'linear-ifgnls': dict(quadratic=False),
    'uncensored-ifgnls': dict(censor=False),
}


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--bootstrap-reps', type=int, default=400)
    parser.add_argument('--n-jobs', type=int, default=2)
    parser.add_argument('--rep-timeout', type=float, default=60)
    parser.add_argument('--mp-context', choices=['spawn', 'forkserver', 'fork'])
    parser.add_argument('--cases', nargs='+', choices=CASES, default=list(CASES))
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'validation/analytical_se/small4')
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    data_path = Path(BENCH) / 'small4.dta'
    frame = pd.read_stata(data_path)
    common = dict(shares=SHARES, prices=PRICES, expenditure='total',
                  demographics=DEMOS, anot=10.0, start='zero',
                  algorithm='gn', blas_threads=1, verbose=False)
    records = []
    for number, name in enumerate(args.cases):
        options = CASES[name]
        begin = time.perf_counter()
        base = quaidsce(frame, **common, **options)
        point_seconds = time.perf_counter()-begin
        begin = time.perf_counter()
        fitted = quaidsce(frame, **common, **options, analytic=True)
        analytic_seconds = time.perf_counter()-begin
        for field in ['theta', 'b', 'sigma', 'tau']:
            if getattr(base, field) is not None:
                np.testing.assert_array_equal(getattr(base, field), getattr(fitted, field))
        assert base.n_gn == fitted.n_gn and base.n_outer == fitted.n_outer
        inference = fitted.analytical
        vector = fitted.elas.as_stata_vector()
        standard_errors = inference.elasticity_se
        n = fitted.spec.neqn
        names = ([f'income:{i+1}' for i in range(n)]
                 + [f'uncompensated:{i+1},{j+1}' for i in range(n) for j in range(n)]
                 + [f'compensated:{i+1},{j+1}' for i in range(n) for j in range(n)])
        table = pd.DataFrame(dict(name=names, estimate=vector,
                                  analytical_se=standard_errors))
        record = dict(case=name, converged=base.converged, nobs=base.nobs,
                      n_outer=base.n_outer, n_gn=base.n_gn,
                      point_seconds=point_seconds,
                      analytical_fit_seconds=analytic_seconds,
                      inference_seconds=inference.elapsed_seconds,
                      point_estimates_bitwise_equal=True,
                      bread_condition=inference.bread_condition,
                      max_standardized_score=inference.max_standardized_score,
                      notes=fitted.notes)
        table.to_csv(args.output / f'{name}.csv', index=False, float_format='%.17g')
        if args.bootstrap_reps and fitted.spec.censor:
            begin = time.perf_counter()
            boot = quaidsce(frame, **common, **options,
                            reps=args.bootstrap_reps, seed=17000+number,
                            n_jobs=args.n_jobs, bootstrap_start='zero',
                            rep_timeout=args.rep_timeout, mp_context=args.mp_context)
            errors = boot.se[-len(vector):]
            np.savez_compressed(args.output / f'{name}-bootstrap.npz',
                                b_star=boot.boot.b_star, names=boot.names,
                                point=boot.b, se=boot.se)
            table['bootstrap_se'] = errors
            table['analytical_over_bootstrap'] = standard_errors/errors
            record.update(bootstrap_seconds=time.perf_counter()-begin,
                          bootstrap_requested=boot.boot.reps_requested,
                          bootstrap_successful=boot.boot.reps_ok,
                          bootstrap_failures=len(boot.boot.failures),
                          bootstrap_failure_records=list(boot.boot.failures),
                          median_se_ratio=float(np.median(standard_errors/errors)),
                          min_se_ratio=float(np.min(standard_errors/errors)),
                          max_se_ratio=float(np.max(standard_errors/errors)))
        table.to_csv(args.output / f'{name}.csv', index=False, float_format='%.17g')
        np.savez_compressed(args.output / f'{name}-covariance.npz',
                            covariance=inference.covariance,
                            elasticity_covariance=inference.elasticity_covariance,
                            joint_covariance=inference.joint_covariance,
                            bread=inference.bread, meat=inference.meat)
        records.append(record)
        print(json.dumps(record), flush=True)
    report = dict(dataset=str(data_path.relative_to(ROOT)),
                  dataset_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),
                  source_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'],
                                                         cwd=ROOT, text=True).strip(),
                  source_has_local_changes=bool(subprocess.check_output(
                      ['git', 'status', '--porcelain'], cwd=ROOT, text=True)),
                  source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                 for path in ['src/pyquaidsce/inference.py',
                                              'src/pyquaidsce/estimator.py',
                                              'src/pyquaidsce/nlsur.py',
                                              'src/pyquaidsce/probit.py',
                                              'src/pyquaidsce/elasticities.py']},
                  python=platform.python_version(), numpy=np.__version__,
                  scipy=scipy.__version__, pandas=pd.__version__,
                  common_options=common, bootstrap_reps=args.bootstrap_reps,
                  n_jobs=args.n_jobs, rep_timeout=args.rep_timeout,
                  mp_context=args.mp_context, cases=records)
    (args.output / 'summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__':
    main()
