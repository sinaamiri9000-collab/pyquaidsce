"""
Economic-theory checks on the elasticity formulas, and a smoke test of the two
code paths Stata v2.0 cannot run at all.

The demand-theory identities below must hold *exactly* when the elasticities are
evaluated at shares that are consistent with the model, i.e. at the shares the
estimated parameters imply at the evaluation point:

    Engel aggregation      sum_i w_i e_i        = 1
    Cournot / homogeneity  sum_j eu_ij          = -e_i
    Compensated homogeneity sum_j ec_ij         = 0

They hold only approximately when -- as ``quaidsce_c.ado`` does -- the elasticity
formulas are evaluated at the *observed* mean shares, because those need not
satisfy the model. Both versions are reported.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path

from pyquaidsce import quaidsce
from pyquaidsce.elasticities import Means, elasticities
from pyquaidsce.params import Spec, unpack


def _benchmark_data(path=None):
    if path is None:
        path = (Path(__file__).resolve().parents[1] / "benchmarks" /
                "cquaids_ifgnls_4g_20k" / "data" /
                "benchmark_cquaids_4g_20k.dta")
    return pd.read_stata(path)


SHARES = ["w1", "w2", "w3", "w4"]
PRICES = ["p1", "p2", "p3", "p4"]
DEMOS = ["z1", "z2", "z3"]


def theory_at_model_shares(res, tol=1e-8):
    """Recompute the *latent* elasticities at model-consistent shares."""
    spec = res.spec
    c = res.coefs
    m = res.means
    n = spec.neqn
    L, zbar = m.lnp, m.demo

    lnpindex = res.anot + c.alpha @ L + 0.5 * L @ c.gamma @ L
    lnb = c.beta @ L if spec.quadratic else 0.0
    if spec.ndemo:
        betanz = c.beta + zbar @ c.eta
        lnc = (zbar @ c.eta) @ L
        mbar = 1.0 + c.rho @ zbar
    else:
        betanz, lnc, mbar = c.beta.copy(), 0.0, 1.0
    D = m.lnexp - lnpindex - np.log(mbar)
    q = np.exp(-lnb - lnc)

    # the share the model itself implies at the evaluation point
    w_model = c.alpha + c.gamma @ L + betanz * D
    if spec.quadratic:
        w_model = w_model + c.lam * q * D**2

    mm = Means(w=w_model, lnp=L, lnexp=m.lnexp, demo=zbar,
               cdf=np.ones(n), pdf=np.zeros(n), du=np.ones(n))
    spec_unc = Spec(n, spec.ndemo, spec.quadratic, censor=False)
    el = elasticities(c, spec_unc, mm, a0=res.anot, strict_stata=False)

    out = {
        "sum w_i": float(w_model.sum()),
        "Engel  sum_i w_i e_i - 1": float(w_model @ el.income - 1.0),
        "Cournot max |sum_j eu_ij + e_i|":
            float(np.max(np.abs(el.uncompensated.sum(axis=1) + el.income))),
        "Hicksian max |sum_j ec_ij|":
            float(np.max(np.abs(el.compensated.sum(axis=1)))),
    }
    return out, el


def main():
    df = _benchmark_data()
    print(f"data: {len(df)} obs, shares sum to "
          f"{df[SHARES].to_numpy().sum(axis=1).mean():.10f} on average")

    # ---- 1. the nocensor path: Stata v2.0 errors out here ------------------ #
    print("\n[1] nocensor path (Stata v2.0 raises 'unknown function *()')")
    res = quaidsce(df, shares=SHARES, prices=PRICES, expenditure="total",
                   demographics=DEMOS, anot=10.0, censor=False,
                   method="ifgnls", strict_stata=False, verbose=False)
    c = res.coefs
    print(f"    ran fine: ll = {res.llf:.6f}, {res.n_gn} Gauss-Newton steps")
    print(f"    adding up   sum alpha - 1 = {c.alpha.sum() - 1:+.3e},"
          f"  sum beta = {c.beta.sum():+.3e},"
          f"  sum lambda = {c.lam.sum():+.3e}")
    print(f"    symmetry    max |gamma - gamma'| = "
          f"{np.abs(c.gamma - c.gamma.T).max():.3e}")
    print(f"    homogeneity max |row sum of gamma| = "
          f"{np.abs(c.gamma.sum(axis=1)).max():.3e}")
    print(f"    eta adding up max |row sum| = "
          f"{np.abs(c.eta.sum(axis=1)).max():.3e}")
    ident, _ = theory_at_model_shares(res)
    print("    demand-theory identities at model-consistent shares:")
    for k, v in ident.items():
        print(f"      {k:<34} {v:+.3e}")
    print("    same identities at the OBSERVED mean shares "
          "(what quaidsce_c.ado uses):")
    print(f"      Engel  sum_i w_i e_i - 1           "
          f"{res.means.w @ res.elas.income - 1:+.3e}")
    print(f"      Cournot max |sum_j eu_ij + e_i|    "
          f"{np.max(np.abs(res.elas.uncompensated.sum(axis=1) + res.elas.income)):+.3e}")
    print(f"      Hicksian max |sum_j ec_ij|         "
          f"{np.max(np.abs(res.elas.compensated.sum(axis=1))):+.3e}")

    # ---- 2. noquadratic + demographics + censoring ------------------------- #
    print("\n[2] noquadratic + demographics + censoring: the expenditure "
          "elasticity bug")
    r_bug = quaidsce(df, shares=SHARES, prices=PRICES, expenditure="total",
                     demographics=DEMOS, anot=10.0, quadratic=False,
                     strict_stata=True, verbose=False)
    r_fix = quaidsce(df, shares=SHARES, prices=PRICES, expenditure="total",
                     demographics=DEMOS, anot=10.0, quadratic=False,
                     strict_stata=True, verbose=False)
    r_fix2 = quaidsce(df, shares=SHARES, prices=PRICES, expenditure="total",
                      demographics=DEMOS, anot=10.0, quadratic=False,
                      strict_stata=False, verbose=False)
    print("    good   Stata (strict_stata=True)   corrected "
          "(strict_stata=False)")
    for i, s in enumerate(SHARES):
        print(f"    {s:<7}{r_bug.elas.income[i]:>22.6f}"
              f"{r_fix2.elas.income[i]:>22.6f}")

    # ---- 3. bootstrap ------------------------------------------------------ #
    print("\n[3] bootstrap (Stata: parallel bs, reps())")
    rb = quaidsce(df, shares=SHARES, prices=PRICES, expenditure="total",
                  demographics=DEMOS, anot=10.0, method="fgnls",
                  reps=24, seed=123456, n_jobs=2, verbose=False)
    print(f"    {rb.boot.reps_ok}/{rb.boot.reps_requested} replications "
          f"succeeded")
    k0 = len(rb.names) - (2 * len(SHARES) ** 2 + len(SHARES))
    lo, hi = rb.boot.percentile_ci(95)
    print("    expenditure elasticities: estimate, bootstrap se, "
          "delta-method se, percentile CI")
    for i, s in enumerate(SHARES):
        j = k0 + i
        print(f"      {s:<7}{rb.b[j]:>10.5f}{rb.boot.se[j]:>11.5f}"
              f"{rb.se[j]:>11.5f}   [{lo[j]:>8.4f}, {hi[j]:>8.4f}]")
    print("    (delta-method se is 0 for the elasticity block by construction "
          "in quaidsce;\n     that is exactly why the package bootstraps.)")


if __name__ == "__main__":
    main()
