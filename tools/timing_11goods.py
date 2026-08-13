"""
Realistic timing for Sina's own problem size: 11 goods, ~37000 observations,
3 demographics, method(ifgnls).

Rather than guess, build a system of exactly that shape out of the repository's
own household data (so the conditioning is realistic, not synthetic), resample
it to 37000 rows, and run the estimator end to end.

Also exercises `initial=`, both round-tripping our own estimates and the
warm-start path the bootstrap uses.
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyquaidsce import quaidsce  # noqa: E402
from tools.paths import find_stata_repo  # noqa: E402

REPO = str(find_stata_repo())
GOODS = [1, 2, 4, 5, 6, 7, 8, 9, 12, 14, 15]  # 11 goods, all censored
DEMOS = ["x1", "x2", "x3"]


def build(nobs: int, seed: int = 20260726) -> pd.DataFrame:
    f = os.path.join(REPO, "data", "DS_STATA_3_2_0_pci2sls_.dta")
    cols = ([f"w{i}" for i in GOODS] + [f"p{i}" for i in GOODS]
            + ["total_exp"] + DEMOS)
    d = pd.read_stata(f, columns=cols)
    d = d[np.isfinite(d.to_numpy()).all(axis=1)]
    grp = d[[f"w{i}" for i in GOODS]].to_numpy(float).sum(axis=1)
    d = d[grp > 0].copy()
    grp = grp[grp > 0]
    d["total"] = d["total_exp"].to_numpy(float) * grp
    for i in GOODS:
        d[f"sw{i}"] = d[f"w{i}"].to_numpy(float) / grp
    keep = [f"sw{i}" for i in GOODS] + [f"p{i}" for i in GOODS] + ["total"] + DEMOS
    d = d[keep].reset_index(drop=True)
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, len(d), size=nobs)
    return d.iloc[idx].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nobs", type=int, default=37000)
    ap.add_argument("--method", default="ifgnls")
    ap.add_argument("--reps", type=int, default=0)
    ap.add_argument("--n-jobs", type=int, default=2)
    a = ap.parse_args()

    df = build(a.nobs)
    shares = [f"sw{i}" for i in GOODS]
    prices = [f"p{i}" for i in GOODS]
    print(f"{len(df)} observations, {len(shares)} goods, {len(DEMOS)} "
          f"demographics")
    print("zeros per share: "
          + ", ".join(f"{int((df[c] == 0).sum())}" for c in shares))

    t0 = time.time()
    res = quaidsce(df, shares=shares, prices=prices, expenditure="total",
                   demographics=DEMOS, anot=10.0, method=a.method,
                   algorithm="gn", verbose=True)
    el = time.time() - t0
    print(f"\n>>> {a.method}: {el:.1f} s  "
          f"({el / 60:.1f} min)   outer={res.n_outer}  "
          f"Gauss-Newton steps={res.n_gn}  "
          f"({el / max(res.n_gn, 1):.2f} s per step)")
    print(f"    free parameters = {res.spec.n_free}, ll = {res.llf:.4f}")

    # ---- initial= round trip --------------------------------------------- #
    t0 = time.time()
    res2 = quaidsce(df, shares=shares, prices=prices, expenditure="total",
                    demographics=DEMOS, anot=10.0, method=a.method,
                    algorithm="gn", initial=res.theta, verbose=False)
    el2 = time.time() - t0
    print(f"\n>>> restarted from initial=res.theta: {el2:.1f} s "
          f"({res2.n_gn} steps, {el / max(el2, 1e-9):.1f}x faster)")
    print(f"    max |b2 - b| = {np.abs(res2.b - res.b).max():.3e}  "
          f"(same optimum)")

    if a.reps:
        t0 = time.time()
        rb = quaidsce(df, shares=shares, prices=prices, expenditure="total",
                      demographics=DEMOS, anot=10.0, method=a.method,
                      algorithm="gn", reps=a.reps, n_jobs=a.n_jobs,
                      seed=1, verbose=True)
        elb = time.time() - t0
        print(f"\n>>> bootstrap {a.reps} reps on {a.n_jobs} cores: "
              f"{elb:.1f} s ({elb / 60:.1f} min) -> "
              f"{elb / a.reps:.1f} s per replication")
        print(f"    projected for 200 reps on C cores: "
              f"{elb / a.reps * 200 * a.n_jobs / 60:.0f} / C minutes")


if __name__ == "__main__":
    main()
