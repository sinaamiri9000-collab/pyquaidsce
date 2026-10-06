"""QUAIDS inference checks with distinct, correlated purchase decisions.

The four modeled goods are a partial basket. An outside good closes the budget.
Given X, each purchased share is wstar_i + delta_i*phi_i/Phi_i + bounded noise.
Purchase shocks are correlated Gaussian variables with unit marginal variance.
Thus E[w_i|X] = Phi_i*wstar_i + delta_i*phi_i exactly, for every correlation.
No observed share is clipped, trimmed, or renormalized.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import platform
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
import pandas as pd
import scipy
from scipy.stats import norm, qmc
from threadpoolctl import threadpool_limits

from pyquaidsce import Spec, quaidsce
from pyquaidsce.elasticities import Means, elasticities
from pyquaidsce.params import full_vector, unpack
from pyquaidsce.selection import FirstStageLayout

ALPHA = np.array([.15, .18, .20, .17])
BETA = np.array([.020, -.014, .016, -.022])
GAMMA = -.012*(np.eye(4)-np.ones((4, 4))/4)
LAMBDA = np.array([.003, -.001, .0015, -.0035])
DELTA = np.array([.0025, -.0015, .0010, -.0020])
ETA = np.array([.008, -.015, .015, -.008])
RHO = .35
TAU = np.array([
    [-.5, .2, -.1, .05, .65, .35, .8, -1.6],
    [.1, -.4, .15, -.1, .45, -.3, -.6, -.65],
    [-.2, .05, -.3, .2, .7, .15, .7, -1.1],
    [.15, -.1, .05, -.45, .55, -.25, -.5, -1.35],
])
ANOT = 1.
NOISE = .008
SPEC = Spec(4, 1, quadratic=True, censor=True)
SHARES = ['w1', 'w2', 'w3', 'w4']
LOG_PRICES = ['lp1', 'lp2', 'lp3', 'lp4']
LAYOUT = FirstStageLayout(
    ordered_names=('p1', 'p2', 'p3', 'p4', 'M', 'z', 's'),
    demand_price_names=tuple(LOG_PRICES),
    price_positions={name: j for j, name in enumerate(LOG_PRICES)},
    expenditure_position=4, covariate_positions={'z': 5, 's': 6},
    selection_cf_position=None, constant_position=7,
)
OPTIONS = dict(shares=SHARES, lnprices=LOG_PRICES, lnexpenditure='lm',
               demographics=['z'], selection_covariates=['z', 's'],
               anot=ANOT, start='zero', algorithm='gn', blas_threads=1,
               verbose=False)


def true_theta():
    gamma = np.concatenate([GAMMA[j:3, j] for j in range(3)])
    return np.concatenate([ALPHA, BETA, gamma, LAMBDA, DELTA, ETA[:3], [RHO]])


def regressors(unit):
    return np.column_stack([unit[:, :4]-.5, 2+1.4*(unit[:, 4]-.5),
                            2*unit[:, 5] - 1, 2*unit[:, 6] - 1])


def indices(reg):
    return np.column_stack([reg, np.ones(len(reg))]) @ TAU.T


def latent(reg):
    """Independent implementation of the published QUAIDS mean equation."""
    lnp, lm, z = reg[:, :4], reg[:, 4], reg[:, 5]
    a = ANOT+lnp @ ALPHA+.5*np.einsum('ti,ij,tj->t', lnp, GAMMA, lnp)
    d = lm-a-np.log1p(RHO*z)
    b = BETA+z[:, None]*ETA
    q = np.exp(-np.einsum('ti,ti->t', b, lnp))
    return ALPHA+lnp @ GAMMA.T+b*d[:, None]+LAMBDA*(q*d*d)[:, None]


def mean_model(reg):
    xb = indices(reg)
    return norm.cdf(xb)*latent(reg)+norm.pdf(xb)*DELTA


def support_bounds():
    """Conservative bounds valid over the full covariate and noise support."""
    a_lower = ANOT-.5*np.abs(ALPHA).sum()-.125*np.abs(GAMMA).sum()
    a_upper = ANOT+.5*np.abs(ALPHA).sum()+.125*np.abs(GAMMA).sum()
    d = max(abs(1.3-a_upper-np.log1p(RHO)),
            abs(2.7-a_lower-np.log1p(-RHO)))
    q = np.exp(.5*(np.abs(BETA).sum()+np.abs(ETA).sum()))
    center = 2*TAU[:, 4]+TAU[:, -1]
    radius = (.5*np.abs(TAU[:, :4]).sum(axis=1)
              +.7*np.abs(TAU[:, 4])+np.abs(TAU[:, 5:7]).sum(axis=1))
    xb_min = center-radius
    mills_max = norm.pdf(xb_min)/norm.cdf(xb_min)
    lower_latent = (ALPHA-.5*np.abs(GAMMA).sum(axis=1)
                    -(np.abs(BETA)+np.abs(ETA))*d-np.abs(LAMBDA)*q*d*d)
    lower_positive = lower_latent+np.minimum(DELTA, 0)*mills_max-1.5*NOISE
    # Sum wstar is ALPHA.sum(); all wstar_i are positive. Any purchased
    # subset is bounded by this full sum, positive corrections and max noise.
    outside_lower = 1-ALPHA.sum()-np.maximum(DELTA, 0) @ mills_max-6*NOISE
    return lower_positive, float(outside_lower)


def realize(reg, rng, correlation):
    if not 0 <= correlation < 1:
        raise ValueError('correlation must lie in [0, 1)')
    xb = indices(reg)
    shock = (np.sqrt(correlation)*rng.normal(size=(len(reg), 1))
             +np.sqrt(1-correlation)*rng.normal(size=(len(reg), 4)))
    purchased = xb+shock > 0
    noise = NOISE*(rng.uniform(-1, 1, (len(reg), 4))
                   +.5*rng.uniform(-1, 1, (len(reg), 1)))
    positive = latent(reg)+norm.pdf(xb)/norm.cdf(xb)*DELTA+noise
    shares = purchased*positive
    outside = 1-shares.sum(axis=1)
    if np.any(positive <= 0) or np.any(outside <= 0):
        raise ValueError('DGP support bounds failed')
    return shares, outside


def generate(observations, seed, correlation=.5):
    rng = np.random.default_rng(seed)
    reg = regressors(rng.random((observations, 7)))
    shares, outside = realize(reg, rng, correlation)
    frame = pd.DataFrame(np.column_stack([shares, reg]),
                         columns=SHARES+LOG_PRICES+['lm', 'z', 's'])
    frame['outside_share'] = outside
    return frame


def population_means(power=20):
    engine = qmc.Sobol(7, scramble=False)
    total, cdf, pdf = np.zeros(4), np.zeros(4), np.zeros(4)
    chunk = min(16384, 2**power)
    for _ in range(2**power//chunk):
        reg = regressors(engine.random(chunk))
        xb = indices(reg)
        total += mean_model(reg).sum(axis=0)
        cdf += norm.cdf(xb).sum(axis=0)
        pdf += norm.pdf(xb).sum(axis=0)
    return Means(total/2**power, np.zeros(4), 2., np.zeros(1),
                 cdf/2**power, pdf/2**power, 2*TAU[:, 4]+TAU[:, -1])


def truth(means):
    theta, tau = true_theta(), TAU.ravel()
    e = elasticities(unpack(theta, SPEC), SPEC, means, ANOT,
                     tau=tau, np_prob=LAYOUT.width, layout=LAYOUT).as_stata_vector()
    names = SPEC.full_names(['z'])+LAYOUT.tau_names(4)+SPEC.elas_names()
    functional = (list(names[:-36])+[f'income:{i+1}' for i in range(4)]
                  +[f'uncompensated:{i+1},{j+1}' for i in range(4) for j in range(4)]
                  +[f'compensated:{i+1},{j+1}' for i in range(4) for j in range(4)])
    return np.concatenate([full_vector(theta, SPEC), tau, e]), names, functional


def contrast_matrix():
    """All 630 pairwise differences, fixed independently of simulation results."""
    i, j = np.triu_indices(36, k=1)
    matrix = np.zeros((len(i), 36))
    matrix[np.arange(len(i)), i] = 1
    matrix[np.arange(len(i)), j] = -1
    return matrix, i, j


def coverage_interval(coverage, repetitions):
    """Wilson interval for the Monte Carlo probability, not the model target."""
    z = norm.ppf(.975)
    denom = 1+z*z/repetitions
    middle = (coverage+z*z/(2*repetitions))/denom
    half = z*np.sqrt(coverage*(1-coverage)/repetitions
                     +z*z/(4*repetitions**2))/denom
    return middle-half, middle+half


def calibration(point, se, target, names):
    sd = point.std(axis=0, ddof=1)
    covered = np.mean(np.abs(point-target) <= norm.ppf(.975)*se, axis=0)
    lower, upper = coverage_interval(covered, len(point))
    return pd.DataFrame(dict(
        functional_name=names, truth=target, mean_estimate=point.mean(axis=0),
        bias=point.mean(axis=0)-target, empirical_sd=sd,
        mean_analytical_se=se.mean(axis=0), rms_analytical_se=np.sqrt((se*se).mean(axis=0)),
        mean_se_over_empirical_sd=np.divide(se.mean(axis=0), sd),
        rms_se_over_empirical_sd=np.divide(np.sqrt((se*se).mean(axis=0)), sd),
        coverage_95=covered, coverage_mcse=np.sqrt(covered*(1-covered)/len(point)),
        coverage_wilson_lower=lower, coverage_wilson_upper=upper))


def _one(task):
    number, seed, observations, correlation, methods, timeout = task
    frame = generate(observations, seed, correlation)
    output = []
    for method in methods:
        begin = time.perf_counter()
        record = dict(replication=number, seed=int(seed), method=method)
        try:
            fit = quaidsce(frame, **OPTIONS, method=method, analytic=True,
                           _deadline=begin+timeout)
            record.update(converged=True, seconds=time.perf_counter()-begin,
                          n_outer=fit.n_outer, n_gn=fit.n_gn,
                          inference_seconds=fit.analytical.elapsed_seconds,
                          bread_condition=fit.analytical.bread_condition,
                          max_standardized_score=fit.analytical.max_standardized_score,
                          min_adjusted_share=float(np.min(
                              fit.means.w*fit.means.cdf+fit.coefs.delta*fit.means.pdf)))
            output.append((record, fit.b, fit.analytic_se,
                           fit.analytical.elasticity_covariance))
        except Exception as exc:
            record.update(converged=False, seconds=time.perf_counter()-begin,
                          error=f'{type(exc).__name__}: {exc}')
            output.append((record, None, None, None))
    return output


def summarize(entries, target, names, functional, output, method):
    entries.sort(key=lambda item: item[0])
    if len(entries) < 2:
        return dict(successful=len(entries), error='fewer than two successful fits')
    point = np.stack([v for _, v, _, _ in entries])
    se = np.stack([v for _, _, v, _ in entries])
    cov = np.stack([v for _, _, _, v in entries])
    report = calibration(point, se, target, functional)
    report.insert(0, 'name', names)
    report.to_csv(output/f'{method}.csv', index=False, float_format='%.17g')
    matrix, i, j = contrast_matrix()
    contrast_points = point[:, -36:] @ matrix.T
    contrast_variance = cov[:, i, i]+cov[:, j, j]-2*cov[:, i, j]
    if np.min(contrast_variance) < -1e-12:
        raise ValueError('negative contrast variance')
    contrast_se = np.sqrt(np.maximum(contrast_variance, 0))
    contrast_names = [f'{functional[-36:][a]} minus {functional[-36:][b]}'
                      for a, b in zip(i, j)]
    contrasts = calibration(contrast_points, contrast_se,
                            matrix @ target[-36:], contrast_names)
    contrasts.to_csv(output/f'{method}-contrasts.csv', index=False, float_format='%.17g')
    empirical_cov = np.cov(point[:, -36:], rowvar=False, ddof=1)
    mean_cov = cov.mean(axis=0)
    a, b = np.triu_indices(36)
    sd = np.sqrt(np.diag(empirical_cov))
    normalized_error = (mean_cov-empirical_cov)/np.outer(sd, sd)
    pd.DataFrame(dict(functional_a=np.asarray(functional[-36:])[a],
                      functional_b=np.asarray(functional[-36:])[b],
                      empirical_covariance=empirical_cov[a, b],
                      mean_analytical_covariance=mean_cov[a, b],
                      error_in_empirical_sd_units=normalized_error[a, b])).to_csv(
        output/f'{method}-covariance.csv', index=False, float_format='%.17g')
    np.savez_compressed(output/f'{method}-draws.npz', point=point, se=se,
                        elasticity_covariance=cov, truth=target, names=names,
                        functional_names=functional,
                        replication=np.asarray([i for i, _, _, _ in entries]))
    e = report.iloc[-36:]
    return dict(successful=len(entries),
                median_elasticity_se_ratio=float(e.mean_se_over_empirical_sd.median()),
                min_elasticity_se_ratio=float(e.mean_se_over_empirical_sd.min()),
                max_elasticity_se_ratio=float(e.mean_se_over_empirical_sd.max()),
                median_elasticity_coverage=float(e.coverage_95.median()),
                min_elasticity_coverage=float(e.coverage_95.min()),
                max_elasticity_coverage=float(e.coverage_95.max()),
                median_contrast_se_ratio=float(contrasts.mean_se_over_empirical_sd.median()),
                min_contrast_coverage=float(contrasts.coverage_95.min()),
                max_contrast_coverage=float(contrasts.coverage_95.max()),
                covariance_normalized_rmse=float(np.sqrt(np.mean(normalized_error**2))))


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--observations', type=int, default=2000)
    parser.add_argument('--repetitions', type=int, default=400)
    parser.add_argument('--seed', type=int, default=27043)
    parser.add_argument('--correlation', type=float, default=.5)
    parser.add_argument('--methods', nargs='+', choices=['fgnls', 'ifgnls'],
                        default=['fgnls', 'ifgnls'])
    parser.add_argument('--n-jobs', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'validation/analytical_se/separate_participation')
    args = parser.parse_args()
    positive_lower, outside_lower = support_bounds()
    if np.any(positive_lower <= 0) or outside_lower <= 0:
        raise ValueError('DGP must guarantee positive shares and remaining budget')
    args.output.mkdir(parents=True, exist_ok=True)
    begin = time.perf_counter()
    with threadpool_limits(limits=1, user_api='blas'):
        means = population_means(20)
        target, names, functional = truth(means)
        coarse, _, _ = truth(population_means(18))
    seeds = np.random.SeedSequence(args.seed).generate_state(args.repetitions, dtype=np.uint32)
    tasks = [(i+1, int(seed), args.observations, args.correlation, args.methods, args.timeout)
             for i, seed in enumerate(seeds)]
    records, draws = [], {method: [] for method in args.methods}
    with mp.get_context('spawn').Pool(args.n_jobs) as pool:
        for count, results in enumerate(pool.imap_unordered(_one, tasks), 1):
            for record, point, se, covariance in results:
                records.append(record)
                if point is not None:
                    draws[record['method']].append((record['replication'], point, se, covariance))
            if count % 10 == 0:
                print(f'{count}/{args.repetitions} independent datasets completed; '
                      f'{time.perf_counter()-begin:.1f} seconds', flush=True)
                pd.DataFrame(records).to_csv(args.output/'replications.csv', index=False)
    reports = {method: summarize(entries, target, names, functional, args.output, method)
               for method, entries in draws.items()}
    pd.DataFrame(records).sort_values(['method', 'replication']).to_csv(
        args.output/'replications.csv', index=False)
    source_paths = ['src/pyquaidsce/inference.py', 'tools/simulate_separate_participation.py']
    derivative_path = ROOT/'src/pyquaidsce/_elasticity_derivatives.py'
    if derivative_path.exists():
        source_paths.append(str(derivative_path.relative_to(ROOT)))
    report = dict(observations=args.observations, repetitions=args.repetitions,
                  seed=args.seed, participation_shock_correlation=args.correlation,
                  python=platform.python_version(), numpy=np.__version__,
                  scipy=scipy.__version__, pandas=pd.__version__,
                  methods=reports, options=OPTIONS, timeout=args.timeout, n_jobs=args.n_jobs,
                  tolerances=dict(param_tol=1e-5, objective_tol=1e-7, gn_tol=1e-5,
                                  outer_param_tol=1e-5, max_iter=300, max_outer=200),
                  elapsed_seconds=time.perf_counter()-begin,
                  true_free_parameters=true_theta().tolist(), true_tau=TAU.tolist(),
                  purchased_share_lower_bound=positive_lower.tolist(),
                  outside_share_lower_bound=outside_lower,
                  population_mean_share=means.w.tolist(),
                  population_purchase_probability=means.cdf.tolist(),
                  max_target_integration_difference=float(np.max(np.abs(target-coarse))),
                  git_commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT,
                                                      text=True).strip(),
                  source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                 for path in source_paths})
    (args.output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
