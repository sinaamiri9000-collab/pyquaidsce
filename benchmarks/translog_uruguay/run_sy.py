"""Reproducible S&Y translog fit using the package's standard SUR settings.

The supplied Stata output is uncensored and is not an external S&Y reference.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parents[1]))
sys.path.insert(0, str(ROOT.parents[1] / "src"))

from pyquaidsce import translog
from pyquaidsce.translog_model import _inner
from benchmarks.translog_uruguay.run_python import DEMOS, SHARES, load_data


def export_result(result, frame, prices, cleaned, elapsed, output, data_path, seed):
    output.mkdir(parents=True, exist_ok=True)
    latent = result.predict(kind="latent_shares")
    predicted = result.predict()
    inner = _inner(result.theta, result.demand_data, result.spec)
    point = pd.DataFrame({name: [result.data[name].mean()]
                          for name in prices+["gasto_total"]+DEMOS})
    h = 1e-6
    plus, minus = point.copy(), point.copy()
    plus["gasto_total"] *= np.exp(h)
    minus["gasto_total"] *= np.exp(-h)
    numerical = (np.log(result.predict(plus, kind="quantities"))
                 - np.log(result.predict(minus, kind="quantities")))/(2*h)
    expenditure_error = float(np.max(abs(numerical[0]-result.elas.income)))
    price_error = 0.0
    for j, price in enumerate(prices):
        plus, minus = point.copy(), point.copy()
        plus[price] *= np.exp(h)
        minus[price] *= np.exp(-h)
        numerical = (np.log(result.predict(plus, kind="quantities"))
                     - np.log(result.predict(minus, kind="quantities")))/(2*h)
        price_error = max(price_error, float(np.max(abs(numerical[0]-result.elas.uncompensated[:, j]))))
    summary = dict(
        model="basic indirect translog + Shonkwiler-Yen + demographic translation",
        nobs=result.nobs, original_rows=len(frame), cleaned_negative_shares=cleaned,
        method=result.method, algorithm=result.algorithm, stop_rule=result.stop_rule,
        start=result.start_method, initialization=result.initialization,
        converged=result.converged, nrtol=result.nrtol, n_outer=result.n_outer,
        n_gn=result.n_gn, runtime_seconds=elapsed, llf=result.llf,
        latent_adding_up=result.spec.latent_adding_up, alpha_sum=float(result.coefs.alpha.sum()),
        n_free=result.spec.n_free, n_probit_coefficients=result.tau.size,
        demographic_names=DEMOS, min_effective_expenditure_ratio=float(inner.b.min()),
        denominator_min=float(inner.denominator.min()), denominator_max=float(inner.denominator.max()),
        fitted_share_sum_min=float(predicted.sum(axis=1).min()),
        fitted_share_sum_max=float(predicted.sum(axis=1).max()),
        latent_share_sum_min=float(latent.sum(axis=1).min()),
        latent_share_sum_max=float(latent.sum(axis=1).max()),
        negative_fitted_share_values=int(np.sum(predicted < 0)),
        sigma_condition_number=float(np.linalg.cond(result.sigma)), r2=result.r2.tolist(),
        predicted_share_evaluation=result.elas.we.tolist(),
        expenditure_elasticity_difference=expenditure_error, price_elasticity_difference=price_error,
        bootstrap_reps=0 if result.boot is None else result.boot.reps_requested,
        bootstrap_successful_reps=0 if result.boot is None else result.boot.reps_ok,
        bootstrap_failures=[] if result.boot is None else result.boot.failures,
        bootstrap_seed=seed if result.boot is not None else None,
        inference="conditional analytical SEs" if result.boot is None else "full two-step bootstrap SEs",
        external_sy_reference_available=False,
        data_sha256=hashlib.sha256(data_path.read_bytes()).hexdigest(),
    )
    summary["passed"] = bool(result.converged and result.nobs == 6573
        and not result.spec.latent_adding_up and result.spec.n_eq_estimated == 14
        and inner.b.min() > 0 and expenditure_error < 1e-7 and price_error < 1e-7
        and np.isfinite(result.b).all() and np.isfinite(result.V).all())
    (output / "sy_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    (output / "python_summary.txt").write_text(result.summary()+"\n"+result.elasticity_tables()+"\n")
    pd.DataFrame(dict(name=result.names, estimate=result.b, se=result.se,
                      analytic_se=result.analytic_se)).to_csv(output / "coefficients_elasticities.csv", index=False)
    np.savez_compressed(output / "python_results.npz", theta=result.theta, tau=result.tau,
        b=result.b, V=result.V, V_est=result.V_est, sigma=result.sigma,
        names=np.array(result.names), r2=result.r2,
        income=result.elas.income, uncompensated=result.elas.uncompensated,
        compensated=result.elas.compensated, predicted_shares=predicted,
        latent_shares=latent, cdf=result.demand_data.cdf,
        selection_index=result.selection_index, estimation_mask=result.estimation_mask)
    if result.boot is not None:
        np.savez_compressed(output / "bootstrap_results.npz", b_star=result.boot.b_star,
                           se=result.boot.se, V=result.boot.V)
    return summary


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path, default=ROOT / "data/bd_uruguay.csv")
    parser.add_argument("--output", type=Path, default=ROOT / "results_sy")
    parser.add_argument("--start", choices=["zero", "mean", "nested"], default="nested")
    parser.add_argument("--reps", type=int, default=0)
    parser.add_argument("--seed", type=int, default=20261001)
    parser.add_argument("--jobs", type=int, default=1)
    args = parser.parse_args()
    frame, prices, cleaned = load_data(args.data)
    started = time.perf_counter()
    result = translog(frame, SHARES, prices=prices, expenditure="gasto_total",
        demographics=DEMOS, censor=True, start=args.start, reps=args.reps, seed=args.seed,
        bootstrap_start="warm", n_jobs=args.jobs, mp_context="spawn")
    summary = export_result(result, frame, prices, cleaned, time.perf_counter()-started,
                            args.output, args.data, args.seed)
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit("The S&Y fit did not pass its convergence/derivative checks")


if __name__ == "__main__":
    main()
