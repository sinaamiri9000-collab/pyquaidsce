"""
Optimiser-independent validation of the elasticity code and of ``e(V)``.

For each of the four Stata runs in ``bench/small4.log`` we take Stata's *own*
17-digit parameter vector, feed it into our elasticity routine and our
delta-method covariance, and compare against the ``ELAS_INC`` / ``ELAS_UNCOMP``
/ ``ELAS_COMP`` blocks and the standard errors Stata reported for that very
vector.

Because the parameters are held fixed, nothing here depends on the optimiser,
the starting values or the convergence tolerance: any disagreement would be a
genuine difference in the formulas.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))
from pyquaidsce.elasticities import elasticities, sample_means        # noqa: E402
from pyquaidsce.estimator import first_stage                          # noqa: E402
from pyquaidsce.jacfree import make_cache                             # noqa: E402
from pyquaidsce.model import DemandData, residuals                    # noqa: E402
from pyquaidsce.nlsur import (LOG2PI, _normal_equations,              # noqa: E402
                              _whitener)
from pyquaidsce.params import (Spec, delta_matrix, full_vector,        # noqa: E402
                              unpack)
from tools.validate_small4 import (BENCH, DEMOS, PRICES, RUNS, SHARES,  # noqa: E402
                                   parse)

N_G, N_R = 4, 2


def theta_from(ref, quad=True):
    th = [ref[f"alpha:alpha_{i}"][0] for i in range(1, N_G + 1)]
    th += [ref[f"beta:beta_{i}"][0] for i in range(1, N_G + 1)]
    for j in range(1, N_G):
        for i in range(j, N_G):
            th.append(ref[f"gamma:gamma_{i}_{j}"][0])
    if quad:
        th += [ref[f"lambda:lambda_{i}"][0] for i in range(1, N_G + 1)]
    th += [ref[f"delta:delta_{i}"][0] for i in range(1, N_G + 1)]
    for v in DEMOS:
        th += [ref[f"eta:eta_{v}_{i}"][0] for i in range(1, N_G)]
    th += [ref[f"rho:rho_{v}"][0] for v in DEMOS]
    return np.array(th)


def main():
    coefs, lls = parse(os.path.join(BENCH, "small4.log"))
    df = pd.read_stata(os.path.join(BENCH, "small4.dta"))
    W = df[SHARES].to_numpy(float)
    lnp = np.log(df[PRICES].to_numpy(float))
    lnexp = np.log(df["total"].to_numpy(float))
    Z = df[DEMOS].to_numpy(float)
    fs = first_stage(W, lnp, lnexp, Z, predict="pr")
    d = DemandData(lnp=lnp, lnexp=lnexp, shares=W, demo=Z,
                   cdf=fs.cdf, pdf=fs.pdf, a0=10.0)

    print("first-stage probit coefficients (tau), all runs identical")
    ref_tau = np.array(
        [coefs[1][f"tau:{nm}"][0] for i in range(1, N_G + 1)
         for nm in ([f"p{j}_{i}" for j in range(1, N_G + 1)] + [f"M_{i}"]
                    + [f"{v}_{i}" for v in DEMOS] + [f"cons_{i}"])]
    )
    print(f"  max |tau_py - tau_stata| = {np.abs(fs.tau - ref_tau).max():.3e}")
    ref_setau = np.array(
        [coefs[1][f"tau:{nm}"][1] for i in range(1, N_G + 1)
         for nm in ([f"p{j}_{i}" for j in range(1, N_G + 1)] + [f"M_{i}"]
                    + [f"{v}_{i}" for v in DEMOS] + [f"cons_{i}"])]
    )
    print("  max |se(tau)_py - se(tau)_stata| = "
          f"{np.abs(np.sqrt(np.diag(fs.setau)) - ref_setau).max():.3e}")

    rows = []
    for run in (1, 2, 3, 4):
        if run not in coefs:
            continue
        cfg = RUNS[run]
        spec = Spec(N_G, N_R, cfg["quadratic"], cfg["censor"])
        th = theta_from(coefs[run], quad=cfg["quadratic"])
        ref = coefs[run]

        # --- structural coefficients: restriction round trip --------------- #
        bfull_mine = full_vector(th, spec)
        bfull_st = np.array([ref[nm][0] for nm in spec.full_names(DEMOS)])
        d_restr = float(np.abs(bfull_mine - bfull_st).max())

        # --- log-likelihood at Stata's theta ------------------------------- #
        u = residuals(th, d, spec)
        sigma = (u.T @ u) / d.nobs
        _, ld = np.linalg.slogdet(sigma)
        ll = (-(d.nobs * spec.n_eq_estimated / 2.0) * (1.0 + LOG2PI)
              - (d.nobs / 2.0) * ld)
        d_ll = abs(ll - lls[run]) / abs(lls[run])

        # --- delta-method standard errors ---------------------------------- #
        # For NLS and converged IFGNLS, the residual Sigma at the reported
        # theta is the covariance used in e(V). For two-step FGNLS, Stata uses
        # Sigma from the *preceding NLS fit*, which cannot be recovered from the
        # final FGNLS vector alone; deliberately report it as not verified here.
        d_se = np.nan
        if cfg["method"] in {"nls", "ifgnls"}:
            P = _whitener(sigma)
            G, _, _ = _normal_equations(th, d, spec, make_cache(spec), P, 4000)
            Dl = delta_matrix(spec)
            se = np.sqrt(np.diag(Dl @ np.linalg.inv(G) @ Dl.T))
            se_st = np.array([ref[nm][1] for nm in spec.full_names(DEMOS)])
            ok_se = np.isfinite(se_st)
            d_se = (float(np.abs(se[ok_se] - se_st[ok_se]).max())
                    if ok_se.any() else np.nan)

        # --- elasticities at Stata's theta --------------------------------- #
        c = unpack(th, spec)
        means = sample_means(d, fs.du, spec)
        el = elasticities(c, spec, means, a0=10.0, tau=fs.tau,
                          np_prob=fs.np_prob, strict_stata=True)
        ev = el.as_stata_vector()
        ev_st = np.array([ref[nm][0] for nm in spec.elas_names()])
        d_ei = float(np.abs(ev[:N_G] - ev_st[:N_G]).max())
        d_eu = float(np.abs(ev[N_G:N_G + N_G**2]
                            - ev_st[N_G:N_G + N_G**2]).max())
        d_ec = float(np.abs(ev[N_G + N_G**2:] - ev_st[N_G + N_G**2:]).max())

        rows.append(dict(run=run,
                         spec=("fgnls" if run == 1 else
                               "nls" if run == 2 else
                               "ifgnls" if run == 3 else "noquadratic"),
                         restrictions=d_restr, rel_ll=d_ll, se=d_se,
                         ELAS_INC=d_ei, ELAS_UNCOMP=d_eu, ELAS_COMP=d_ec))

    out = pd.DataFrame(rows)
    pd.set_option("display.width", 200)
    print("\nMaximum absolute deviation, evaluated at STATA's own parameter "
          "vector\n(so the optimiser plays no role at all):\n")
    print(out.to_string(index=False, float_format=lambda x: f"{x:.3e}"))
    print("\nA figure around 1e-15..1e-13 means the formula is identical up to "
          "floating-point\nassociativity; 1e-9..1e-8 means Stata's own 17-digit "
          "dump was the limiting factor.")


if __name__ == "__main__":
    main()
