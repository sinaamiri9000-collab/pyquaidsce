"""Compare the current package with the captured official v1.0.1 baseline."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

from pyquaidsce import quaidsce  # noqa: E402
from pyquaidsce.estimator import first_stage  # noqa: E402
from pyquaidsce.model import DemandData, fitted_shares, jacobian_full  # noqa: E402


SHARES = ["sw1", "sw2", "sw4", "sw9"]
PRICES = ["p1", "p2", "p4", "p9"]
DEMOS = ["x1", "x2"]


def _max_abs(now, old):
    return float(np.max(np.abs(np.asarray(now) - np.asarray(old))))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline",
        type=Path,
        default=ROOT.parents[1] / "baseline_1.0.1" / "small4_ifgnls",
    )
    parser.add_argument("--json-out", type=Path)
    args = parser.parse_args()

    df = pd.read_stata(ROOT / "bench" / "small4.dta")
    res = quaidsce(
        df,
        shares=SHARES,
        prices=PRICES,
        expenditure="total",
        demographics=DEMOS,
        anot=10.0,
        method="ifgnls",
        verbose=False,
    )
    W = df[SHARES].to_numpy(float)
    lnp = np.log(df[PRICES].to_numpy(float))
    lnexp = np.log(df["total"].to_numpy(float))
    demo = df[DEMOS].to_numpy(float)
    fs = first_stage(W, lnp, lnexp, demo, predict="pr")
    d = DemandData(lnp, lnexp, W, demo, fs.cdf, fs.pdf, 10.0)

    arrays = {
        "coefficients": res.b,
        "structural_free_parameters": res.theta,
        "covariance": res.V,
        "structural_covariance": res.V_est,
        "sigma": res.sigma,
        "first_stage_params": res.tau,
        "first_stage_covariance": res.setau,
        "fitted_shares": fitted_shares(res.theta, d, res.spec),
        "jacobian_full": jacobian_full(res.theta, d, res.spec),
        "elasticities_income": res.elas.income,
        "elasticities_marshallian": res.elas.uncompensated,
        "elasticities_hicksian": res.elas.compensated,
    }
    differences = {
        name: _max_abs(value, np.load(args.baseline / f"{name}.npy"))
        for name, value in arrays.items()
    }
    metadata = json.loads((args.baseline / "metadata.json").read_text())
    differences["llf"] = abs(float(res.llf) - float(metadata["llf"]))
    names_equal = res.names == metadata["parameter_names"]
    passed = names_equal and max(differences.values()) <= 1e-10
    report = {
        "baseline_version": "1.0.1",
        "current_version": __import__("pyquaidsce").__version__,
        "dataset": "bench/small4.dta",
        "specification": metadata["specification"],
        "parameter_names_and_order_equal": names_equal,
        "max_absolute_differences": differences,
        "algebraic_tolerance": 1e-10,
        "full_estimation_tolerance": 1e-8,
        "passed": passed,
    }
    rendered = json.dumps(report, indent=2, ensure_ascii=False)
    print(rendered)
    if args.json_out is not None:
        args.json_out.parent.mkdir(parents=True, exist_ok=True)
        args.json_out.write_text(rendered + "\n", encoding="utf-8")
    if not passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
