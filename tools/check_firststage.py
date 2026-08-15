"""Cheap discriminating test: does our first-stage probit reproduce the
`tau` block of log/censor_v2.log?  If it does, we have the right sample.
"""
from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))
from pyquaidsce.estimator import first_stage           # noqa: E402
from tools.statalog import parse_coef_table            # noqa: E402
from tools.paths import find_stata_repo                 # noqa: E402

REPO = str(find_stata_repo())
SHARES = [f"w{i}" for i in range(1, 18)]
PRICES = [f"p{i}" for i in range(1, 18)]
DEMOS = ["x1", "x2", "x3"]


def load(which: str, expvar: str):
    f = (f"{REPO}/data/DS_STATA_3_2__pci2sls_.dta" if which == "A"
         else f"{REPO}/data/DS_STATA_3_2_0_pci2sls_.dta")
    d = pd.read_stata(f, columns=SHARES + PRICES + DEMOS + [expvar])
    d = d[np.isfinite(d.to_numpy()).all(axis=1)]
    d = d[d[expvar] > 0].reset_index(drop=True)
    return d


def run(which: str, expvar: str, predict: str):
    d = load(which, expvar)
    W = d[SHARES].to_numpy(float)
    lnp = np.log(d[PRICES].to_numpy(float))
    lnexp = np.log(d[expvar].to_numpy(float))
    Z = d[DEMOS].to_numpy(float)
    fs = first_stage(W, lnp, lnexp, Z, predict=predict)
    return d, fs


if __name__ == "__main__":
    ref, sc, _ = parse_coef_table(f"{REPO}/log/censor_v2.log")
    tau_ref = np.array(
        [ref[f"tau:{nm}"][0] for i in range(1, 18)
         for nm in ([f"p{j}_{i}" for j in range(1, 18)] + [f"M_{i}"]
                    + [f"x{r}_{i}" for r in (1, 2, 3)] + [f"cons_{i}"])]
    )
    se_ref = np.array(
        [ref[f"tau:{nm}"][1] for i in range(1, 18)
         for nm in ([f"p{j}_{i}" for j in range(1, 18)] + [f"M_{i}"]
                    + [f"x{r}_{i}" for r in (1, 2, 3)] + [f"cons_{i}"])]
    )

    for which in ("A", "B"):
        for expvar in ("total_exp_DS", "total_exp"):
            try:
                d, fs = run(which, expvar, "pr")
            except Exception as exc:  # noqa: BLE001
                print(f"{which:>2} {expvar:<14} FAILED: {exc}")
                continue
            n = len(d)
            dev = np.abs(fs.tau - tau_ref)
            rel = dev / np.maximum(np.abs(tau_ref), 1e-8)
            sedev = np.abs(np.sqrt(np.diag(fs.setau)) - se_ref)
            print(f"{which:>2} {expvar:<14} N={n:<6} "
                  f"max|dtau|={dev.max():.3e}  max rel={rel.max():.3e}  "
                  f"max|dse|={sedev.max():.3e}")
            if dev.max() < 1e-4:
                print("   *** MATCH ***  worst 5 coefficients:")
                order = np.argsort(-dev)[:5]
                for k in order:
                    print(f"     idx {k:4d}  py={fs.tau[k]:+.8f} "
                          f"stata={tau_ref[k]:+.8f}")
