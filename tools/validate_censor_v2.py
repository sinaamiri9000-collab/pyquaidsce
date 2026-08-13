"""
Full validation against the author's own Stata run.

Benchmark: quaidsce-master/log/censor_v2.log

    . quaidsce w1-w17, anot(10) prices(p1-p17) expenditure(total) nolog
          demographics(x1-x3)
    (obs = 15,147)
    ... 18 FGNLS iterations ...
    Censored Quadratic AIDS model
    Number of obs = 15147, demographics = 3, Alpha_0 = 10, ll = 482472.11
    649 reported coefficients

The `total` variable of that run is reproduced exactly by
`total_exp_DS` (= total_exp * sum(w1..w17)) in data/DS_STATA_3_2_0_pci2sls_.dta,
restricted to its 15147 strictly positive values; the first-stage probits
confirm this to 4e-6 (Stata's own display precision).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyquaidsce import quaidsce                        # noqa: E402
from tools.statalog import parse_coef_table            # noqa: E402
from tools.paths import find_stata_repo                 # noqa: E402

REPO = str(find_stata_repo())
SHARES = [f"w{i}" for i in range(1, 18)]
PRICES = [f"p{i}" for i in range(1, 18)]
DEMOS = ["x1", "x2", "x3"]
EXPVAR = "total_exp_DS"


def load(normalize: bool = True) -> pd.DataFrame:
    """Rebuild the estimation sample of log/censor_v2.log.

    `total` in that run is `total_exp * sum(w1..w17)`, which the shipped data
    already carries as `total_exp_DS`; its 15147 strictly positive values are
    exactly the reported sample size.  With `normalize=True` the shares are
    rescaled to sum to one within the group, i.e. the run is a *conditional*
    demand system -- which is what the magnitude of the reported alphas implies
    (they range over -0.19 .. 0.32, impossible for raw shares averaging 0.004).
    """
    f = os.path.join(REPO, "data", "DS_STATA_3_2_0_pci2sls_.dta")
    d = pd.read_stata(f, columns=SHARES + PRICES + DEMOS + [EXPVAR])
    d = d[np.isfinite(d.to_numpy()).all(axis=1)]
    d = d[d[EXPVAR] > 0].reset_index(drop=True)
    if normalize:
        grp = d[SHARES].to_numpy(float).sum(axis=1)
        for c in SHARES:
            d[c] = d[c].to_numpy(float) / grp
    return d


def reference():
    ref, sc, _ = parse_coef_table(os.path.join(REPO, "log", "censor_v2.log"))
    return ref, sc


def stata_name_order(neqn=17, demos=DEMOS):
    """The order in which censor_v2.log lists its 649 coefficients."""
    names = [f"alpha:alpha_{i}" for i in range(1, neqn + 1)]
    names += [f"beta:beta_{i}" for i in range(1, neqn + 1)]
    for j in range(1, neqn + 1):
        for i in range(j, neqn + 1):
            names.append(f"gamma:gamma_{i}_{j}")
    names += [f"lambda:lambda_{i}" for i in range(1, neqn + 1)]
    names += [f"delta:delta_{i}" for i in range(1, neqn + 1)]
    for v in demos:
        names += [f"eta:eta_{v}_{i}" for i in range(1, neqn + 1)]
    names += [f"rho:rho_{v}" for v in demos]
    for i in range(1, neqn + 1):
        names += [f"tau:p{j}_{i}" for j in range(1, neqn + 1)]
        names.append(f"tau:M_{i}")
        names += [f"tau:{v}_{i}" for v in demos]
        names.append(f"tau:cons_{i}")
    return names


def compare(res, ref, sc, tag, out_dir):
    order = stata_name_order()
    mine = dict(zip(res.names, res.b))
    mse = dict(zip(res.names, res.se))

    rows = []
    for nm in order:
        if nm not in ref or nm not in mine:
            continue
        b_s, se_s = ref[nm]
        rows.append((nm, mine[nm], b_s, mse[nm], se_s))
    df = pd.DataFrame(rows, columns=["name", "b_py", "b_stata",
                                     "se_py", "se_stata"])
    df["abs_db"] = (df.b_py - df.b_stata).abs()
    df["abs_dse"] = (df.se_py - df.se_stata).abs()
    # Stata prints 7 significant digits, so this is the achievable floor
    df["tol_b"] = np.maximum(np.abs(df.b_stata), 1e-7) * 5e-7
    df["rel_db"] = df.abs_db / np.maximum(np.abs(df.b_stata), 1e-8)

    ll_s = sc["Log-likelihood"]
    print(f"\n===== {tag} =====")
    print(f"N            python {res.nobs:<10}  stata {int(sc['Number of obs'])}")
    print(f"log-lik      python {res.llf:<18.6f} stata {ll_s}  "
          f"|diff| = {abs(res.llf - ll_s):.6g}")
    print(f"coefficients compared: {len(df)} / 649")
    print(f"  max |b_py - b_stata|   = {df.abs_db.max():.3e}"
          f"   (at {df.loc[df.abs_db.idxmax(), 'name']})")
    print(f"  max relative deviation = {df.rel_db.max():.3e}")
    print(f"  within display precision (7 sig digits): "
          f"{int((df.abs_db <= df.tol_b).sum())} / {len(df)}")
    print(f"  max |se_py - se_stata|  = {df.abs_dse.max():.3e}")
    print("\n  worst 12 coefficients by |b diff|:")
    w = df.nlargest(12, "abs_db")
    for _, r in w.iterrows():
        print(f"    {r['name']:<22} py={r.b_py:+.8f}  stata={r.b_stata:+.8f}"
              f"  d={r.abs_db:.2e}   se py={r.se_py:.8f} stata={r.se_stata:.8f}")
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, f"compare_{tag}.csv"), index=False)
    return df, abs(res.llf - ll_s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", default="ifgnls")
    ap.add_argument("--predict", default="pr", choices=["pr", "xb"])
    ap.add_argument("--vce", default="objective", choices=["objective", "final"])
    ap.add_argument("--max-iter", type=int, default=400)
    ap.add_argument("--out", default="out")
    ap.add_argument("--save-summary", action="store_true")
    ap.add_argument("--gn-verbose", action="store_true")
    ap.add_argument("--start", default="zero", choices=["zero", "linear"])
    ap.add_argument("--algorithm", default="lm", choices=["lm", "gn"])
    ap.add_argument("--initial", default="none", choices=["none", "stata"])
    ap.add_argument("--raw-shares", action="store_true")
    a = ap.parse_args()

    d = load(normalize=not a.raw_shares)
    ref, sc = reference()
    init = None
    if a.initial == 'stata':
        from tools.check_at_stata_theta import stata_theta
        init = stata_theta(ref)
    t0 = time.time()
    res = quaidsce(
        d, shares=SHARES, prices=PRICES, expenditure=EXPVAR,
        demographics=DEMOS, anot=10.0, method=a.method,
        first_stage_predict=a.predict, vce_sigma=a.vce, start=a.start, algorithm=a.algorithm,
        max_iter=a.max_iter, verbose=True, gn_verbose=a.gn_verbose,
        initial=init,
    )
    el = time.time() - t0
    tag = f"{a.method}_{a.predict}_{a.vce}_{a.start}_{a.algorithm}_{a.initial}"
    print(f"\nelapsed {el:.1f}s   outer={res.n_outer}  GN steps={res.n_gn}")
    df, dll = compare(res, ref, sc, tag, a.out)

    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, f"summary_{tag}.txt"), "w") as fh:
        fh.write(res.summary(elasticities=False))
        fh.write("\n\n")
        fh.write(res.elasticity_tables())
    with open(os.path.join(a.out, f"meta_{tag}.json"), "w") as fh:
        json.dump(dict(elapsed=el, n_outer=res.n_outer, n_gn=res.n_gn,
                       llf=res.llf, ll_stata=sc["Log-likelihood"],
                       dll=dll, max_db=float(df.abs_db.max()),
                       max_dse=float(df.abs_dse.max()),
                       nobs=res.nobs), fh, indent=2)


if __name__ == "__main__":
    main()
