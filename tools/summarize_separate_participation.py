"""Combine the four saved separate-participation Monte Carlo scenarios."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import norm

ROOT = Path(__file__).resolve().parents[1]
SCENARIOS = ['n500-corr05', 'n2000-corr05', 'n8000-corr05', 'n2000-corr0']


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--directory', type=Path,
                        default=ROOT/'validation/analytical_se/separate_participation')
    args = parser.parse_args()
    rows, tails, failures = [], [], []
    critical = norm.ppf(.975)
    for scenario in SCENARIOS:
        path = args.directory/scenario
        summary = json.loads((path/'summary.json').read_text())
        records = pd.read_csv(path/'replications.csv')
        requested = summary['repetitions']
        for method, values in summary['methods'].items():
            good = records[(records.method == method) & records.converged]
            row = dict(scenario=scenario, method=method,
                       observations=summary['observations'],
                       shock_correlation=summary['participation_shock_correlation'],
                       requested=requested, **values,
                       median_seconds=float(good.seconds.median()),
                       median_n_outer=float(good.n_outer.median()),
                       median_n_gn=float(good.n_gn.median()),
                       failures=requested-values['successful'])
            row['min_interval_produced_and_covers_rate'] = (
                row['min_elasticity_coverage']*values['successful']/requested)
            rows.append(row)
            with np.load(path/f'{method}-draws.npz') as data:
                point, se, target = data['point'][:, -36:], data['se'][:, -36:], data['truth'][-36:]
                student = (point-target)/se
                for j, name in enumerate(data['functional_names'][-36:]):
                    tails.append(dict(
                        scenario=scenario, method=method, functional_name=str(name),
                        successful=len(point), requested=requested,
                        mean_se=float(se[:, j].mean()), median_se=float(np.median(se[:, j])),
                        max_se=float(se[:, j].max()),
                        max_se_replication=int(data['replication'][se[:, j].argmax()]),
                        studentized_q025=float(np.quantile(student[:, j], .025)),
                        studentized_q975=float(np.quantile(student[:, j], .975)),
                        below_minus_1_96=float(np.mean(student[:, j] < -critical)),
                        above_plus_1_96=float(np.mean(student[:, j] > critical)),
                        interval_produced_and_covers_rate=float(
                            np.sum(np.abs(student[:, j]) <= critical)/requested)))
        bad = records[~records.converged].copy()
        bad.insert(0, 'scenario', scenario)
        failures.append(bad)
    pd.DataFrame(rows).to_csv(args.directory/'comparison.csv', index=False, float_format='%.17g')
    pd.DataFrame(tails).to_csv(args.directory/'studentized-tails.csv', index=False, float_format='%.17g')
    pd.concat(failures, ignore_index=True).to_csv(args.directory/'failures.csv', index=False)
    print(pd.DataFrame(rows)[['scenario', 'method', 'successful', 'median_elasticity_se_ratio',
                             'min_elasticity_coverage', 'median_elasticity_coverage']].to_string(index=False))


if __name__ == '__main__':
    main()
