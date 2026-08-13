"""
Why do the ``method(nls)`` estimates differ from Stata's?

Two candidate explanations:

  (a) we converge to a different local minimum of the NLS criterion;
  (b) Stata stops early -- ``nlsur``'s default ``nrtolerance()`` is 1e-5, four
      orders of magnitude looser than what we use.

The test discriminates cleanly.  Evaluate our criterion at *Stata's* reported
parameter vector and look at (i) the objective value and (ii) the relative
gradient.  If Stata's objective is higher than ours *and* its gradient is far
from zero, Stata simply stopped short; if Stata's objective were lower, we would
be the ones in the wrong basin.

The same test is applied to the ``fgnls`` run, where any early stop in the NLS
step propagates through ``Sigma_hat`` into the second-step estimates.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyquaidsce.estimator import first_stage                          # noqa: E402
from pyquaidsce.jacfree import make_cache                             # noqa: E402
from pyquaidsce.model import DemandData, residuals                    # noqa: E402
from pyquaidsce.nlsur import (_normal_equations, _objective,          # noqa: E402
                              _solve_scaled, _whitener, nlsur)
from pyquaidsce.params import Spec                                    # noqa: E402
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES, parse  # noqa: E402

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
    spec = Spec(N_G, N_R, True, True)
    cache = make_cache(spec)
    Im = np.eye(N_G)

    def diag(theta, P, tag):
        obj = _objective(theta, d, spec, P)
        G, g, _ = _normal_equations(theta, d, spec, cache, P, 4000)
        step = _solve_scaled(G, g)
        nrtol = abs(float(step @ g)) / max(abs(obj), 1e-300)
        print(f"    {tag:<34} objective = {obj:.12f}   relative gradient = "
              f"{nrtol:.3e}")
        return obj, nrtol

    print("NLS criterion (Sigma = I)")
    th_s = theta_from(coefs[2])
    res = nlsur(d, spec, method="nls", algorithm="gn", start="zero",
                max_iter=400, chunk=4000)
    o_s, n_s = diag(th_s, Im, "at Stata's method(nls) theta")
    o_p, n_p = diag(res.theta, Im, "at our method(nls) optimum")
    print(f"    our objective is {'LOWER' if o_p < o_s else 'HIGHER'} by "
          f"{abs(o_p - o_s):.6g}  ({abs(o_p - o_s) / o_s:.2e} relative)")
    print("    => " + (
        "Stata stopped early (its nrtolerance default is 1e-5); the two are not"
        " the same point on the criterion."
        if o_p < o_s and n_s > 1e-6 else
        "we may be in a different basin -- investigate."))

    print("\nHow far apart are the two parameter vectors?")
    rel = np.abs(res.theta - th_s) / np.maximum(np.abs(th_s), 1e-8)
    print(f"    max |dtheta| = {np.abs(res.theta - th_s).max():.3e},"
          f" max relative = {rel.max():.3e}")

    print("\nSecond-step consequence (fgnls uses Sigma_hat from the NLS step)")
    for tag, th in (("Stata's NLS residuals", th_s),
                    ("our NLS residuals", res.theta)):
        u = residuals(th, d, spec)
        sig = (u.T @ u) / d.nobs
        _, ld = np.linalg.slogdet(sig)
        print(f"    {tag:<26} ln|Sigma_hat| = {ld:.10f}")

    print("\nRefit fgnls using Stata's own first-step Sigma_hat")
    u = residuals(th_s, d, spec)
    sig_s = (u.T @ u) / d.nobs
    from pyquaidsce.nlsur import gauss_newton
    th2, obj2, it2, ok2 = gauss_newton(th_s, d, spec, sig_s, max_iter=300,
                                       chunk=4000, algorithm="gn")
    th_f = theta_from(coefs[1])
    print(f"    max |dtheta| vs Stata's fgnls estimates = "
          f"{np.abs(th2 - th_f).max():.3e}   (GN steps {it2})")
    res1 = nlsur(d, spec, method="fgnls", algorithm="gn", start="zero",
                 max_iter=400, chunk=4000)
    print(f"    max |dtheta| from our own fgnls run     = "
          f"{np.abs(res1.theta - th_f).max():.3e}")


if __name__ == "__main__":
    main()
