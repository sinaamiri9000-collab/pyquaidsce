"""Check post-estimation analytical derivatives on the small4 benchmark."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
import pandas as pd
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from pyquaidsce import DemandData, quaidsce
from pyquaidsce.inference import _EstimatingSystem, compute_analytical_inference
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'validation/analytical_se/elasticity-derivative-check.json')
    args = parser.parse_args()
    frame = pd.read_stata(Path(BENCH)/'small4.dta')
    options = dict(shares=SHARES, prices=PRICES, expenditure='total',
                   demographics=DEMOS, anot=10., verbose=False, blas_threads=1)
    records = []
    with threadpool_limits(limits=1, user_api='blas'):
        for method in ['fgnls', 'ifgnls']:
            plain = quaidsce(frame, **options, method=method)
            fit = quaidsce(frame, **options, method=method, analytic=True)
            for name in ['theta', 'b', 'sigma', 'tau', 'llf', 'n_outer', 'n_gn', 'converged']:
                np.testing.assert_array_equal(getattr(plain, name), getattr(fit, name))
            prices = np.log(frame[PRICES].to_numpy())
            expenditure = np.log(frame['total'].to_numpy())
            demo = frame[DEMOS].to_numpy()
            design = np.column_stack([prices, expenditure, demo])
            xb = np.column_stack([design, np.ones(len(frame))]) @ fit.tau.reshape(4, -1).T
            data = DemandData(prices, expenditure, frame[SHARES].to_numpy(),
                              demo, norm.cdf(xb), norm.pdf(xb), fit.anot)
            initial = quaidsce(frame, **options, method='nls').theta if method == 'fgnls' else None
            system = _EstimatingSystem(fit, data, design, initial, 2000, None)
            analytical, numerical = [], []
            for _ in range(9):
                begin = time.perf_counter()
                a = np.column_stack(system.elasticity_jacobians())
                analytical.append(time.perf_counter()-begin)
                begin = time.perf_counter()
                b = np.column_stack(system.elasticity_jacobians(relative_step=1e-6))
                numerical.append(time.perf_counter()-begin)
            error = np.max(np.abs(a-b)/(1+np.abs(b)))
            reference = compute_analytical_inference(
                fit, data, design, theta_nls=initial, _relative_step=1e-6)
            relative_se_error = np.max(np.abs(
                fit.analytical.elasticity_se/reference.elasticity_se-1))
            # Relative error of a nearly zero covariance entry is unstable.
            # Scale by the two marginal SEs, retaining every covariance entry.
            covariance_error = np.max(np.abs(
                fit.analytical.elasticity_covariance-reference.elasticity_covariance
            )/np.outer(reference.elasticity_se, reference.elasticity_se))
            assert error < 3e-6
            assert relative_se_error < 5e-5
            assert covariance_error < 5e-5
            records.append(dict(method=method, point_results_bitwise_equal=True,
                                max_scaled_gradient_error=float(error),
                                max_relative_elasticity_se_error=float(relative_se_error),
                                max_covariance_error_in_se_units=float(covariance_error),
                                median_analytical_gradient_seconds=float(np.median(analytical)),
                                median_numerical_gradient_seconds=float(np.median(numerical)),
                                repeats=9, n_outer=fit.n_outer, n_gn=fit.n_gn))
    files = ['src/pyquaidsce/inference.py', 'src/pyquaidsce/_elasticity_derivatives.py',
             'src/pyquaidsce/elasticities.py', 'tools/check_elasticity_derivatives.py']
    report = dict(dataset='small4.dta', observations=len(frame), options=options,
                  cases=records,
                  source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                 for path in files})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
