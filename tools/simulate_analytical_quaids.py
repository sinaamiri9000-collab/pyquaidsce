"""Known-parameter Monte Carlo for the censored QUAIDS mean model.

All four modeled goods share a Probit participation event. When they are
purchased, bounded positive shares sum to one; an outside category receives
the budget otherwise. The exact conditional mean is Phi(xb)*wstar+delta*phi(xb).
No generated observation is trimmed, normalized, or clipped after generation.
This controlled case is not a general model of separate purchase decisions.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
import pandas as pd
from scipy.stats import norm, qmc
from threadpoolctl import threadpool_limits

from pyquaidsce import Spec, quaidsce
from pyquaidsce.elasticities import Means, elasticities
from pyquaidsce.params import full_vector, unpack
from pyquaidsce.selection import legacy_layout

ALPHA = np.array([.20, .25, .30, .25])
BETA = np.array([.04, -.02, .01, -.03])
GAMMA = -.025*(np.eye(4)-np.ones((4, 4))/4)
LAMBDA = np.array([.006, -.002, .003, -.007])
DELTA = np.array([.010, -.005, -.003, -.002])
ETA = np.array([.005, -.003, .002, -.004])
RHO = .1
TAU = np.array([-.6, .5, -.3, .4, .8, .25, -1.6])
ANOT = 1.
NOISE = .02
SPEC = Spec(4, 1, quadratic=True, censor=True)
SHARES = ['w1', 'w2', 'w3', 'w4']
LOG_PRICES = ['lp1', 'lp2', 'lp3', 'lp4']
OPTIONS = dict(shares=SHARES, lnprices=LOG_PRICES,
               lnexpenditure='lm', demographics=['z'], anot=ANOT,
               start='zero', algorithm='gn', blas_threads=1, verbose=False)


def true_theta():
    g = np.concatenate([GAMMA[j:3, j] for j in range(3)])
    return np.concatenate([ALPHA, BETA, g, LAMBDA, DELTA, ETA[:3], [RHO]])


def regressors(unit):
    return np.column_stack([unit[:, :4]-.5,
                            2+1.4*(unit[:, 4]-.5), unit[:, 5]-.5])


def latent(reg):
    """Generate from the QUAIDS equations independently of the package model."""
    lnp, lm, z = reg[:, :4], reg[:, 4], reg[:, 5]
    a = ANOT+lnp @ ALPHA+.5*np.einsum('ti,ij,tj->t', lnp, GAMMA, lnp)
    d = lm-a-np.log1p(RHO*z)
    b = BETA+z[:, None]*ETA
    q = np.exp(-np.einsum('ti,ti->t', b, lnp))
    return ALPHA+lnp @ GAMMA.T+b*d[:, None]+LAMBDA*(q*d*d)[:, None]


def mean_model(reg):
    index = np.column_stack([reg, np.ones(len(reg))]) @ TAU
    return norm.cdf(index)[:, None]*latent(reg)+norm.pdf(index)[:, None]*DELTA


def positivity_bound():
    """Conservative bound over the entire covariate and noise support."""
    a_lower = ANOT-.5*np.abs(ALPHA).sum()-.125*np.abs(GAMMA).sum()
    a_upper = ANOT+.5*np.abs(ALPHA).sum()+.125*np.abs(GAMMA).sum()
    d = max(abs(1.3-a_upper-np.log(1.05)), abs(2.7-a_lower-np.log(.95)))
    q = np.exp(.5*np.abs(BETA).sum()+.25*np.abs(ETA).sum())
    max_index = .5*np.abs(TAU[:4]).sum()+.7*abs(TAU[4])+.5*abs(TAU[5])
    mills = norm.pdf(max_index)/norm.cdf(-max_index)
    return (ALPHA-.5*np.abs(GAMMA).sum(axis=1)
            -(np.abs(BETA)+.5*np.abs(ETA))*d-np.abs(LAMBDA)*q*d*d
            +np.minimum(DELTA, 0)*mills-1.5*NOISE)


def generate(observations, seed):
    rng = np.random.default_rng(seed)
    reg = regressors(rng.random((observations, 6)))
    index = np.column_stack([reg, np.ones(observations)]) @ TAU
    participation = index+rng.normal(size=observations) > 0
    noise = rng.uniform(-1, 1, (observations, 4))
    noise = NOISE*(noise-noise.mean(axis=1, keepdims=True))
    positive = latent(reg)+norm.pdf(index)[:, None]/norm.cdf(index)[:, None]*DELTA+noise
    assert np.all(positive > 0) and np.all(positive < 1)
    np.testing.assert_allclose(positive.sum(axis=1), 1, atol=2e-15)
    shares = participation[:, None]*positive
    frame = pd.DataFrame(np.column_stack([shares, reg]),
                         columns=SHARES+LOG_PRICES+['lm', 'z'])
    # This fifth category closes the observed budget, and is not an estimated
    # equation. Modeled shares sum to one conditional on participation.
    frame['outside_share'] = 1-participation.astype(float)
    return frame


def population_means(power=20):
    """Smooth population moments by deterministic Sobol integration.

    Symmetry gives E[lnp]=0, E[lm]=2, E[z]=0, E[xb]=0 and E[Phi(xb)]=1/2
    exactly. Only E[phi(xb)] and the observed share means need integration.
    """
    engine = qmc.Sobol(6, scramble=False)
    total, pdf = np.zeros(4), 0.
    chunk = min(16384, 2**power)
    for _ in range(2**power//chunk):
        reg = regressors(engine.random(chunk))
        index = np.column_stack([reg, np.ones(chunk)]) @ TAU
        total += mean_model(reg).sum(axis=0)
        pdf += norm.pdf(index).sum()
    return Means(total/2**power, np.zeros(4), 2., np.zeros(1),
                 np.full(4, .5), np.full(4, pdf/2**power), np.zeros(4))


def truth(means):
    theta = true_theta()
    tau = np.tile(TAU, 4)
    layout = legacy_layout(LOG_PRICES, ['z'], include_expenditure=True)
    e = elasticities(unpack(theta, SPEC), SPEC, means, ANOT,
                     tau=tau, np_prob=7, layout=layout).as_stata_vector()
    names = SPEC.full_names(['z'])+layout.tau_names(4)+SPEC.elas_names()
    return np.concatenate([full_vector(theta, SPEC), tau, e]), names


def functional_names(names):
    """Identify elasticities by the actual row-major (good, price) order."""
    return (list(names[:-36])+[f'income:{i+1}' for i in range(4)]
            +[f'uncompensated:{i+1},{j+1}' for i in range(4) for j in range(4)]
            +[f'compensated:{i+1},{j+1}' for i in range(4) for j in range(4)])


def _one(task):
    number, seed, observations, methods, timeout, tolerances = task
    frame = generate(observations, seed)
    output = []
    for method in methods:
        begin = time.perf_counter()
        record = dict(replication=number, seed=int(seed), method=method)
        try:
            fit = quaidsce(frame, **OPTIONS, **tolerances, method=method, analytic=True,
                           _deadline=begin+timeout)
            record.update(converged=True, seconds=time.perf_counter()-begin,
                          n_outer=fit.n_outer, n_gn=fit.n_gn,
                          bread_condition=fit.analytical.bread_condition,
                          max_standardized_score=fit.analytical.max_standardized_score)
            output.append((record, fit.b, fit.analytic_se))
        except Exception as exc:
            record.update(converged=False, seconds=time.perf_counter()-begin,
                          error=f'{type(exc).__name__}: {exc}')
            output.append((record, None, None))
    return output


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--observations', type=int, default=2000)
    parser.add_argument('--repetitions', type=int, default=400)
    parser.add_argument('--seed', type=int, default=17043)
    parser.add_argument('--methods', nargs='+', choices=['nls', 'fgnls', 'ifgnls'],
                        default=['nls', 'fgnls', 'ifgnls'])
    parser.add_argument('--n-jobs', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--param-tol', type=float)
    parser.add_argument('--objective-tol', type=float)
    parser.add_argument('--gn-tol', type=float)
    parser.add_argument('--outer-param-tol', type=float)
    parser.add_argument('--output', type=Path,
                        default=ROOT/'validation/analytical_se/quaids_dgp')
    args = parser.parse_args()
    tolerances = {key: getattr(args, key) for key in
                  ['param_tol', 'objective_tol', 'gn_tol', 'outer_param_tol']
                  if getattr(args, key) is not None}
    assert np.all(positivity_bound() > 0), 'DGP must guarantee positive purchased shares'
    args.output.mkdir(parents=True, exist_ok=True)
    with threadpool_limits(limits=1, user_api='blas'):
        means = population_means(20)
        target, names = truth(means)
        coarse, _ = truth(population_means(18))
    seeds = np.random.SeedSequence(args.seed).generate_state(args.repetitions, dtype=np.uint32)
    tasks = [(i+1, int(seed), args.observations, args.methods, args.timeout, tolerances)
             for i, seed in enumerate(seeds)]
    records, draws = [], {method: [] for method in args.methods}
    with mp.get_context('spawn').Pool(args.n_jobs) as pool:
        for count, output in enumerate(pool.imap_unordered(_one, tasks), 1):
            for record, point, se in output:
                records.append(record)
                if point is not None:
                    draws[record['method']].append((record['replication'], point, se))
            if count % 10 == 0:
                print(f'{count}/{args.repetitions} independent datasets completed', flush=True)
                pd.DataFrame(records).to_csv(args.output/'replications.csv', index=False)
    reports = {}
    for method, entries in draws.items():
        entries.sort(key=lambda item: item[0])
        point = np.stack([v for _, v, _ in entries])
        se = np.stack([v for _, _, v in entries])
        sd = point.std(axis=0, ddof=1)
        coverage = np.mean(np.abs(point-target) <= norm.ppf(.975)*se, axis=0)
        ratio = np.divide(se.mean(axis=0), sd, out=np.full_like(sd, np.nan), where=sd > 0)
        pd.DataFrame(dict(name=names, functional_name=functional_names(names),
                          truth=target, mean_estimate=point.mean(axis=0),
                          bias=point.mean(axis=0)-target, empirical_sd=sd,
                          mean_analytical_se=se.mean(axis=0),
                          mean_se_over_empirical_sd=ratio, coverage_95=coverage)).to_csv(
            args.output/f'{method}.csv', index=False, float_format='%.17g')
        np.savez_compressed(args.output/f'{method}-draws.npz', point=point, se=se,
                            truth=target, names=names, functional_names=functional_names(names),
                            replication=np.asarray([i for i, _, _ in entries]))
        reports[method] = dict(successful=len(entries), requested=args.repetitions,
                               median_elasticity_se_ratio=float(np.median(ratio[-36:])),
                               min_elasticity_se_ratio=float(np.min(ratio[-36:])),
                               max_elasticity_se_ratio=float(np.max(ratio[-36:])),
                               median_elasticity_coverage=float(np.median(coverage[-36:])),
                               min_elasticity_coverage=float(np.min(coverage[-36:])),
                               max_elasticity_coverage=float(np.max(coverage[-36:])))
    pd.DataFrame(records).sort_values(['method', 'replication']).to_csv(
        args.output/'replications.csv', index=False)
    report = dict(observations=args.observations, repetitions=args.repetitions,
                  seed=args.seed, methods=reports, options=dict(OPTIONS, **tolerances),
                  timeout=args.timeout, n_jobs=args.n_jobs,
                  true_free_parameters=true_theta().tolist(), true_tau_per_good=TAU.tolist(),
                  purchased_share_lower_bound=positivity_bound().tolist(),
                  population_mean_share=means.w.tolist(),
                  population_mean_phi=means.pdf.tolist(),
                  max_target_integration_difference=float(np.max(np.abs(target-coarse))),
                  source_sha256={path: hashlib.sha256((ROOT/path).read_bytes()).hexdigest()
                                 for path in ['src/pyquaidsce/inference.py',
                                              'tools/simulate_analytical_quaids.py']})
    (args.output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
