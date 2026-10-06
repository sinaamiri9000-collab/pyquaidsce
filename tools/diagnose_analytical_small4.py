"""Diagnose bootstrap tails and starting-value sensitivity on small4.

All converged draws remain in the reported covariance: no trimming or clipping.
The paired cold/warm fits use the same samples and the existing solver settings.
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
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / 'src'))

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES

_CONTEXT = {}


def _init(frame, options, theta, sigma, timeout):
    _CONTEXT.update(frame=frame, options=options, theta=theta, sigma=sigma,
                    timeout=timeout)


def _one(task):
    mode, number, seed = task
    begin = time.perf_counter()
    context = _CONTEXT
    frame = context['frame']
    indices = np.random.default_rng(seed).integers(0, len(frame), size=len(frame))
    options = dict(context['options'])
    if mode == 'warm':
        options.update(initial=context['theta'], sigma_initial=context['sigma'])
    record = dict(mode=mode, replication=number, seed=int(seed))
    try:
        result = quaidsce(frame.iloc[indices].reset_index(drop=True), **options,
                          _deadline=begin+context['timeout'])
        record.update(converged=result.converged, seconds=time.perf_counter()-begin,
                      n_outer=result.n_outer, n_gn=result.n_gn, llf=result.llf,
                      income_1=result.elas.income[0],
                      min_sigma_eigenvalue=float(np.linalg.eigvalsh(result.sigma).min()),
                      parameter_distance=float(np.max(np.abs(result.theta-context['theta'])
                                                       /(1+np.abs(context['theta'])))))
        record.update({f'adjusted_share_{i+1}': float(v) for i, v in enumerate(result.elas.we)})
        vector = result.b if result.converged else None
    except Exception as exc:
        record.update(converged=False, seconds=time.perf_counter()-begin,
                      error=f'{type(exc).__name__}: {exc}')
        vector = None
    return record, vector


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--reps', type=int, default=400)
    parser.add_argument('--seed', type=int, default=17000)
    parser.add_argument('--n-jobs', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=60)
    parser.add_argument('--output', type=Path,
                        default=ROOT / 'validation/analytical_se/diagnostics')
    args = parser.parse_args()
    frame = pd.read_stata(Path(BENCH)/'small4.dta')
    options = dict(shares=SHARES, prices=PRICES, expenditure='total',
                   demographics=DEMOS, anot=10.0, method='ifgnls',
                   start='zero', algorithm='gn', blas_threads=1, verbose=False)
    point = quaidsce(frame, **options, analytic=True)
    seeds = np.random.SeedSequence(args.seed).generate_state(args.reps, dtype=np.uint32)
    tasks = [(mode, i+1, seed) for mode in ['zero', 'warm'] for i, seed in enumerate(seeds)]
    records, vectors = [], {'zero': [], 'warm': []}
    args.output.mkdir(parents=True, exist_ok=True)
    with mp.get_context('spawn').Pool(
            args.n_jobs, initializer=_init,
            initargs=(frame, options, point.theta, point.sigma, args.timeout)) as pool:
        for count, (record, vector) in enumerate(pool.imap_unordered(_one, tasks), 1):
            records.append(record)
            if vector is not None:
                vectors[record['mode']].append((record['replication'], vector))
            if count % 10 == 0:
                print(f'{count}/{len(tasks)} paired bootstrap fits completed', flush=True)
                pd.DataFrame(records).to_csv(args.output/'replications.csv', index=False)
    report = dict(seed=args.seed, reps=args.reps, timeout=args.timeout,
                  n_jobs=args.n_jobs, options=options,
                  inference_source_sha256=hashlib.sha256((ROOT/'src/pyquaidsce/inference.py').read_bytes()).hexdigest(),
                  point_adjusted_shares=point.elas.we.tolist(),
                  point_income_1=float(point.elas.income[0]), modes={})
    for mode, entries in vectors.items():
        entries.sort(key=lambda item: item[0])
        b = np.vstack([v for _, v in entries])
        np.savez_compressed(args.output/f'{mode}-bootstrap.npz', b_star=b,
                            replication=np.asarray([i for i, _ in entries]),
                            names=point.names, point=point.b)
        se = b[:, -36:].std(axis=0, ddof=1)
        income = b[:, -36]
        report['modes'][mode] = dict(successful=len(entries),
                                    income_1_se=float(se[0]),
                                    income_1_quantiles=np.quantile(income, [0,.01,.5,.99,1]).tolist(),
                                    median_se_ratio=float(np.median(point.analytical.elasticity_se/se)),
                                    min_se_ratio=float(np.min(point.analytical.elasticity_se/se)),
                                    max_se_ratio=float(np.max(point.analytical.elasticity_se/se)))
    pd.DataFrame(records).sort_values(['mode', 'replication']).to_csv(
        args.output/'replications.csv', index=False)
    (args.output/'summary.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2), flush=True)


if __name__ == '__main__':
    main()
