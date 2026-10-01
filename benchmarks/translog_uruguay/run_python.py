"""Reproduce the supplied demandsys benchmark without a Stata installation."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parents[1]
sys.path.insert(0, str(REPO / "src"))

from pyquaidsce import translog
from pyquaidsce.inference import complex_step_jacobian
from pyquaidsce.translog import _point_at_means
from pyquaidsce.translog_model import latent_shares, share_derivatives


SHARES = [f"w{i}_red" for i in range(1, 15)]
DEMOS = ["npersonas", "edad", "Sex", "prim_comp", "sec_comp", "sup_comp", "log_ing"]


def load_data(path=None):
    path = ROOT / "data" / "bd_uruguay.csv" if path is None else Path(path)
    df = pd.read_csv(path, na_values=["."])
    cleaned = int((df[SHARES] < 0).sum().sum())
    df[SHARES] = df[SHARES].clip(lower=0)
    prefix = "P_med" if "P_med1" in df else "p_med"
    return df, [f"{prefix}{i}" for i in range(1, 15)], cleaned


def compare(result, data, reference):
    el, el_se = result.elasticities_at_means(data, standard_errors=True)
    # Diagnostic for the reference's inference convention: its reported price
    # elasticity SEs hold the predicted share denominator fixed. The package's
    # public SEs differentiate the complete elasticity, including that share.
    point = _point_at_means(data, result.price_names, result.expenditure_name,
                            result.demo_names, result.prices_are_logs,
                            result.expenditure_is_log)
    fixed_share = latent_shares(result.theta, point, result.spec)[0]
    conditional_J = complex_step_jacobian(
        lambda t: (share_derivatives(t, point, result.spec)[1][0]
                    / fixed_share[:, None] - np.eye(result.spec.neqn)).ravel(),
        result.theta)
    conditional_se = np.sqrt(np.diag(conditional_J @ result.V_est @ conditional_J.T))
    conditional_se = conditional_se.reshape(result.spec.neqn, result.spec.neqn)
    named_se = dict(zip(result.names, result.se))
    rows = []
    for row in reference["coefficients"]:
        value = result.get(row["name"])
        se = named_se[row["name"]]
        rows.append(dict(group="coefficient", name=row["name"], stata=row["estimate"],
                         python=value, difference=value-row["estimate"],
                         stata_se=row["se"], python_se=se, se_difference=se-row["se"],
                         python_stata_conditional_se=se, stata_conditional_se_difference=se-row["se"]))
    for row in reference["uncompensated"]:
        i, j = row["good"]-1, row["price"]-1
        value, se = el.uncompensated[i, j], el_se.uncompensated[i, j]
        rows.append(dict(group="uncompensated", name=f"good_{i+1}_price_{j+1}",
                         stata=row["estimate"], python=value, difference=value-row["estimate"],
                         stata_se=row["se"], python_se=se, se_difference=se-row["se"],
                         python_stata_conditional_se=conditional_se[i, j],
                         stata_conditional_se_difference=conditional_se[i, j]-row["se"]))
    table = pd.DataFrame(rows)
    c = table[table.group == "coefficient"]
    u = table[table.group == "uncompensated"]
    errors = dict(coefficients=float(c.difference.abs().max()),
                  coefficient_se=float(c.se_difference.abs().max()),
                  uncompensated=float(u.difference.abs().max()),
                  uncompensated_stata_conditional_se=float(u.stata_conditional_se_difference.abs().max()),
                  r2=float(np.max(np.abs(result.r2-reference["r2"]))),
                  llf=abs(result.llf-reference["llf"]))
    # Printed coefficients/elasticities are rounded. Optimization paths and
    # stopping tests also differ; equality beyond these tolerances is not claimed.
    tolerances = dict(coefficients=1e-4, coefficient_se=1e-5,
                      uncompensated=2e-4, uncompensated_stata_conditional_se=1e-5,
                      r2=1e-4, llf=.01)
    summary = dict(nobs=result.nobs, original_rows=len(data), converged=result.converged,
                   llf_python=result.llf, llf_stata=reference["llf"],
                   n_outer=result.n_outer, n_gn=result.n_gn,
                   max_abs_errors=errors, tolerances=tolerances,
                   inference=dict(
                       public_se="full delta method including predicted-share uncertainty",
                       reference_se="conditional delta method with predicted-share denominator fixed",
                       full_delta_vs_reference_max_abs_se_difference=float(u.se_difference.abs().max()),
                       reference_conditional_se_used_only_for_benchmark=True),
                   estimation_start="zero", reference_coefficients_used_as_initial=False,
                   evaluation=dict(prices="arithmetic means on complete price vectors",
                                    expenditure="arithmetic mean on available full-frame rows",
                                    demographics="arithmetic means on complete demographic vectors"))
    summary["passed"] = bool(result.converged and result.nobs == reference["nobs"]
                              and all(errors[k] <= tolerances[k] for k in errors))
    return table, summary, el, el_se


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=Path)
    parser.add_argument("--output", type=Path, default=ROOT / "results")
    args = parser.parse_args()
    df, prices, cleaned = load_data(args.data)
    reference = json.loads((ROOT / "reference" / "stata_reference.json").read_text())
    started = time.perf_counter()
    result = translog(df, SHARES, prices=prices, expenditure="gasto_total",
                       demographics=DEMOS, verbose=True)
    table, summary, el, el_se = compare(result, df, reference)
    summary["runtime_seconds"] = time.perf_counter()-started
    summary["negative_shares_replaced"] = cleaned
    data_path = args.data or ROOT / "data" / "bd_uruguay.csv"
    summary["data_sha256"] = hashlib.sha256(data_path.read_bytes()).hexdigest()
    args.output.mkdir(parents=True, exist_ok=True)
    table.to_csv(args.output / "stata_python_comparison.csv", index=False)
    (args.output / "comparison_summary.json").write_text(json.dumps(summary, indent=2)+"\n")
    (args.output / "python_summary.txt").write_text(result.summary()+"\n"+result.elasticity_tables()+"\n")
    np.savez_compressed(args.output / "python_results.npz", theta=result.theta,
                         b=result.b, V=result.V, V_est=result.V_est, sigma=result.sigma,
                         names=np.array(result.names), r2=result.r2,
                         evaluation_income=el.income, evaluation_uncompensated=el.uncompensated,
                         evaluation_compensated=el.compensated,
                         evaluation_income_se=el_se.income,
                         evaluation_uncompensated_se=el_se.uncompensated,
                         evaluation_compensated_se=el_se.compensated)
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise SystemExit("Uruguay benchmark did not meet its declared tolerances")


if __name__ == "__main__":
    main()
