"""Reproduce two large FGNLS SEs and a limited existing-tolerance diagnostic.

These fits are selected after the Monte Carlo results: they have the largest
income-4 SE in their scenarios. They are not new random validation samples.
The diagnostic does not change the package defaults or the Monte Carlo draws.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from pyquaidsce import DemandData, quaidsce
from pyquaidsce.inference import _EstimatingSystem, compute_analytical_inference
from tools.simulate_separate_participation import generate, OPTIONS, LOG_PRICES, SHARES


def diagnose(observations, correlation, replication, tolerances, numeric_check):
    seed = int(np.random.SeedSequence(27043).generate_state(400, dtype=np.uint32)[replication-1])
    frame = generate(observations, seed, correlation)
    rec = dict(observations=observations, correlation=correlation,
               replication=replication, seed=seed, tolerances=tolerances)
    initial = quaidsce(frame, **OPTIONS, method='nls', **tolerances)
    rec.update(initial_nls_converged=initial.converged, initial_nls_n_gn=initial.n_gn,
               initial_rho=initial.coefs.rho.tolist())
    try:
        fit = quaidsce(frame, **OPTIONS, method='fgnls', analytic=True, **tolerances)
        prices = frame[LOG_PRICES].to_numpy()
        expenditure = frame.lm.to_numpy()
        design = np.column_stack([prices, expenditure, frame[['z', 's']].to_numpy()])
        xb = np.column_stack([design, np.ones(observations)]) @ fit.tau.reshape(4, -1).T
        data = DemandData(prices, expenditure, frame[SHARES].to_numpy(),
                          frame[['z']].to_numpy(), norm.cdf(xb), norm.pdf(xb), fit.anot)
        system = _EstimatingSystem(fit, data, design, initial.theta, 2000, None)
        bread = system.bread()
        rec.update(income4=float(fit.elas.income[3]),
                   income4_se=float(fit.analytical.income_se[3]),
                   min_adjusted_share=float(fit.elas.we.min()),
                   min_eigenvalue_initial_nls_information=float(np.linalg.eigvalsh(
                       -bread[system.nls_slice, system.nls_slice]).min()),
                   min_eigenvalue_final_demand_information=float(np.linalg.eigvalsh(
                       -bread[system.theta_slice, system.theta_slice]).min()),
                   bread_condition=fit.analytical.bread_condition,
                   max_standardized_score=fit.analytical.max_standardized_score,
                   n_outer=fit.n_outer, n_gn=fit.n_gn, final_rho=fit.coefs.rho.tolist())
        if numeric_check:
            reference = compute_analytical_inference(
                fit, data, design, theta_nls=initial.theta, _relative_step=1e-6)
            rec['numerical_gradient_income4_se'] = float(reference.income_se[3])
    except (ValueError, RuntimeError) as exc:
        rec['error'] = f'{type(exc).__name__}: {exc}'
    return rec


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'validation/analytical_se/separate_participation/diagnostics.json')
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api='blas'):
        extreme = [diagnose(8000, .5, 93, {}, True), diagnose(2000, 0., 133, {}, True)]
        accuracy = [extreme[0]]+[
            diagnose(8000, .5, 93, dict(param_tol=p, objective_tol=o, gn_tol=g), False)
            for p, o, g in [(1e-6, 1e-8, 1e-6), (1e-8, 1e-10, 1e-8)]]
    paths = ['src/pyquaidsce/inference.py', 'src/pyquaidsce/_elasticity_derivatives.py',
             'tools/simulate_separate_participation.py', 'tools/diagnose_separate_participation.py']
    report = dict(diagnostic_only=True, package_defaults_changed=False,
                  selection='largest income4 SE in each specified scenario',
                  extreme_cases=extreme, existing_tolerance_probe=accuracy,
                  source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                 for path in paths})
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
