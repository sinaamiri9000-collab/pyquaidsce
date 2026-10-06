"""Known-model Monte Carlo check of the two-step covariance construction.

This is a reference problem, not a simulated QUAIDS dataset: a correctly
specified Probit supplies a generated regressor to a normal linear outcome.
The stages share correlated shocks. Compare Hardin's joint sandwich with the
Murphy--Topel formula in Hardin (2002), equation (1), and empirical dispersion.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from pyquaidsce.inference import _bread_inverse
from pyquaidsce.probit import _lambda_ratio, probit


def run(observations, repetitions, seed):
    rng = np.random.default_rng(seed)
    true_tau = np.array([0.7, -0.8, 0.1])
    true_gamma = 1.4
    true_mean = norm.cdf(true_tau[-1]/np.sqrt(1+true_tau[:-1]@true_tau[:-1]))
    target = true_gamma/true_mean
    rows = []
    for _ in range(repetitions):
        z = np.column_stack([rng.normal(size=(observations, 2)), np.ones(observations)])
        x = rng.normal(size=observations)
        shock = rng.normal(size=observations)
        outcome_shock = 0.6*shock+0.8*rng.normal(size=observations)
        probability = norm.cdf(z @ true_tau)
        outcome = 0.8+0.3*x+true_gamma*probability+outcome_shock
        participation = (z @ true_tau + shock > 0).astype(float)
        first = probit(participation, z[:, :2])
        assert first.converged
        index = z @ first.b
        p, phi = norm.cdf(index), norm.pdf(index)
        design = np.column_stack([np.ones(observations), x, p])
        beta = np.linalg.solve(design.T @ design, design.T @ outcome)
        u = outcome-design @ beta
        q = 2*participation-1
        score_first = z*(q*_lambda_ratio(q*index))[:, None]
        score_second = design*u[:, None]
        a = -np.linalg.inv(first.V)/observations
        d = -(design.T @ design)/observations
        c = -beta[-1]*design.T @ (phi[:, None]*z)/observations
        c[-1] += (phi[:, None]*z).T @ u/observations
        bread = np.block([[a, np.zeros((3, 3))], [c, d]])
        inverse, _ = _bread_inverse(bread, np.ones(6))
        scores = np.column_stack([score_first, score_second])
        scores -= scores.mean(axis=0)
        influence = -scores @ inverse.T
        joint = influence.T @ influence/observations**2
        mean_p = float(p.mean())
        mean_derivative = np.mean(phi[:, None]*z, axis=0)
        mean_influence = p-mean_p+influence[:, :3] @ mean_derivative
        effect_influence = influence[:, -1]/mean_p-beta[-1]*mean_influence/mean_p**2
        effect_se = np.sqrt(effect_influence @ effect_influence)/observations

        # Murphy--Topel with summed scores and inverse observed information.
        sigma2 = float(u @ u/observations)
        s2 = score_second/sigma2
        s2_first = (beta[-1]*u*phi)[:, None]*z/sigma2
        cm, rm = s2.T @ s2_first, s2.T @ score_first
        v1, v2 = first.V, sigma2*np.linalg.inv(design.T @ design)
        mt_second = v2+v2 @ (cm @ v1 @ cm.T-rm @ v1 @ cm.T-cm @ v1 @ rm.T) @ v2
        cross = v2 @ (rm-cm) @ v1
        mt_joint = np.block([[v1, cross.T], [cross, mt_second]])
        # Extend the coefficient covariance to the sample-mean target using
        # the empirical covariance between its raw mean and model influences.
        first_if = observations*score_first @ v1.T
        second_if = observations*(s2-first_if @ (cm/observations).T) @ v2.T
        mt_if = np.column_stack([first_if, second_if])
        raw_mean = p-mean_p
        covariance_mean = mt_if.T @ raw_mean/observations**2
        mean_variance = raw_mean @ raw_mean/observations**2
        gradient = np.zeros(6)
        gradient[:3] = -beta[-1]*mean_derivative/mean_p**2
        gradient[-1] = 1/mean_p
        mean_gradient = -beta[-1]/mean_p**2
        mt_variance = (gradient @ mt_joint @ gradient
                       + 2*mean_gradient*gradient @ covariance_mean
                       + mean_gradient**2*mean_variance)
        assert mt_variance > 0
        rows.append([beta[-1], np.sqrt(joint[-1, -1]),
                     np.sqrt(mt_second[-1, -1]), np.sqrt(v2[-1, -1]),
                     beta[-1]/mean_p, effect_se, np.sqrt(mt_variance)])
    draws = np.asarray(rows)
    output = dict(observations=observations, repetitions=repetitions, seed=seed,
                  true_coefficient=true_gamma, true_population_mean_target=target)
    for label, point, truth, columns in [
            ('coefficient', 0, true_gamma, {'hardin': 1, 'murphy_topel': 2, 'conditional': 3}),
            ('mean_target', 4, target, {'hardin': 5, 'murphy_topel': 6})]:
        output[label] = dict(empirical_sd=float(draws[:, point].std(ddof=1)),
                             bias=float(draws[:, point].mean()-truth), methods={})
        for name, col in columns.items():
            output[label]['methods'][name] = dict(
                mean_se=float(draws[:, col].mean()),
                coverage_95=float(np.mean(np.abs(draws[:, point]-truth) <= 1.96*draws[:, col])))
    return output, draws


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--observations', type=int, default=1000)
    parser.add_argument('--repetitions', type=int, default=500)
    parser.add_argument('--seed', type=int, default=17042)
    parser.add_argument('--output', type=Path, default=ROOT / 'validation/analytical_se/reference')
    args = parser.parse_args()
    with threadpool_limits(limits=1, user_api='blas'):
        report, draws = run(args.observations, args.repetitions, args.seed)
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    np.savetxt(args.output/'draws.csv', draws, delimiter=',', comments='',
               header='coefficient,hardin_se,mt_se,conditional_se,mean_target,hardin_target_se,mt_target_se')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
