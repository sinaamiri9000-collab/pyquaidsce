"""
The Python counterpart of the package's own ``example.do``.

Because ``example.do`` builds its data with Stata's random-number generator
(``set seed 123456`` + ``runiform()``/``rnormal()``), which cannot be reproduced
bit for bit outside Stata, this script works on whatever data you point it at.
Two ways to run it:

    # 1. on the repository's own data (a 4-good conditional system)
    python examples/example_quaidsce.py --repo-data ../src/quaidsce-master

    # 2. on your own file
    python examples/example_quaidsce.py --dta mydata.dta \
        --shares w1 w2 w3 w4 --prices p1 p2 p3 p4 \
        --expenditure total --demographics x1 x2 --anot 10 --reps 200

The Stata command it mirrors:

    quaidsce w1 w2 w3 w4, anot(10) reps(200) prices(p1 p2 p3 p4) ///
        expenditure(total) demographics(x1 x2) nolog
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


def repo_data(repo: str, goods=(1, 2, 4, 9), nobs=2000) -> pd.DataFrame:
    """Build a small conditional 4-good example from an external Stata repo.

    This helper is intended for local experimentation when the upstream Stata
    source tree is available; it is not part of the public benchmark.
    """
    f = os.path.join(repo, "data", "DS_STATA_3_2_0_pci2sls_.dta")
    cols = ([f"w{i}" for i in goods] + [f"p{i}" for i in goods]
            + ["total_exp", "x1", "x2"])
    d = pd.read_stata(f, columns=cols)
    grp = d[[f"w{i}" for i in goods]].to_numpy(float).sum(axis=1)
    d = d[grp > 0].copy()
    grp = grp[grp > 0]
    d["total"] = d["total_exp"].to_numpy(float) * grp
    for i in goods:
        d[f"sw{i}"] = d[f"w{i}"].to_numpy(float) / grp
    keep = ([f"sw{i}" for i in goods] + [f"p{i}" for i in goods]
            + ["total", "x1", "x2"])
    return d[keep].iloc[:nobs].reset_index(drop=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--repo-data")
    ap.add_argument("--dta")
    ap.add_argument("--shares", nargs="+")
    ap.add_argument("--prices", nargs="+")
    ap.add_argument("--expenditure")
    ap.add_argument("--demographics", nargs="*", default=[])
    ap.add_argument("--anot", type=float, default=10.0)
    ap.add_argument("--method", default="fgnls",
                    choices=["nls", "fgnls", "ifgnls"])
    ap.add_argument("--reps", type=int, default=0)
    ap.add_argument("--n-jobs", type=int, default=1)
    ap.add_argument("--seed", type=int, default=123456)
    ap.add_argument("--predict", default="pr", choices=["pr", "xb"])
    ap.add_argument("--no-normalize", action="store_true",
                    help="do not rescale the shares to sum to one")
    a = ap.parse_args()

    if a.repo_data:
        goods = (1, 2, 4, 9)
        df = repo_data(a.repo_data)
        shares = [f"sw{i}" for i in goods]
        prices = [f"p{i}" for i in goods]
        expenditure, demos, anot = "total", ["x1", "x2"], a.anot
    else:
        if not (a.dta and a.shares and a.prices and a.expenditure):
            ap.error("need --repo-data, or --dta with --shares/--prices/"
                     "--expenditure")
        df = pd.read_stata(a.dta)
        shares, prices = a.shares, a.prices
        expenditure, demos, anot = a.expenditure, a.demographics, a.anot
        if not a.no_normalize:
            s = df[shares].to_numpy(float).sum(axis=1)
            ok = s > 0
            df = df[ok].reset_index(drop=True)
            s = s[ok]
            for c in shares:
                df[c] = df[c].to_numpy(float) / s

    print(f"data: {len(df)} observations, {len(shares)} goods, "
          f"{len(demos)} demographics")
    zeros = [(int((df[c] == 0).sum()), c) for c in shares]
    print("censoring: " + ", ".join(f"{c}={z}" for z, c in zeros))

    t0 = time.time()
    res = quaidsce(
        df, shares=shares, prices=prices, expenditure=expenditure,
        demographics=demos or None, anot=anot, method=a.method,
        first_stage_predict=a.predict, reps=a.reps, n_jobs=a.n_jobs,
        seed=a.seed, verbose=True,
    )
    print(f"\nestimated in {time.time() - t0:.1f}s "
          f"({res.n_outer} outer, {res.n_gn} Gauss-Newton steps)")

    print(res.summary(elasticities=False))
    print(res.elasticity_tables())

    if res.notes:
        print("\nNotes on Stata-compatibility:")
        for i, n in enumerate(res.notes, 1):
            print(f"  {i}. {n}")

    if res.boot is not None:
        lo, hi = res.boot.percentile_ci(95)
        print(f"\nbootstrap: {res.boot.reps_ok}/{res.boot.reps_requested} "
              f"replications succeeded")
        print("\nExpenditure elasticities with bootstrap percentile CIs")
        k0 = len(res.names) - (2 * len(shares) ** 2 + len(shares))
        for i, c in enumerate(shares):
            j = k0 + i
            print(f"  {c:<14} {res.b[j]:>10.6f}   se {res.boot.se[j]:>9.6f}"
                  f"   [{lo[j]:>9.6f}, {hi[j]:>9.6f}]")


if __name__ == "__main__":
    main()
