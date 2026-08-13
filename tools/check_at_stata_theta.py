"""
The decisive test of the *model*, independent of the optimiser.

Read Stata's own reported coefficients out of log/censor_v2.log, map them back
into the free parameter vector, and evaluate our objective, Sigma-hat and
log-likelihood there.  If the share equations, the censoring transformation and
the log-likelihood formula are right, e(ll) must come out at 482472.11.

The relative-gradient figure then says whether Stata's point is a stationary
point of *our* criterion -- i.e. whether any remaining disagreement is about the
model or only about which optimum the optimiser walks to.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from pyquaidsce.estimator import first_stage                          # noqa: E402
from pyquaidsce.jacfree import make_cache                             # noqa: E402
from pyquaidsce.model import DemandData, residuals                    # noqa: E402
from pyquaidsce.nlsur import (LOG2PI, _normal_equations, _objective,  # noqa: E402
                              _solve_scaled, _whitener)
from pyquaidsce.params import Spec, delta_matrix, full_vector         # noqa: E402
from tools.statalog import parse_coef_table                           # noqa: E402
from tools.validate_censor_v2 import (DEMOS, EXPVAR, PRICES, REPO,    # noqa: E402
                                      SHARES, load)

N_G, N_R = 17, 3


def stata_theta(ref, neqn=N_G, demos=DEMOS):
    """Reported (full) vector -> free vector.

    Under censoring every alpha, beta, lambda and delta is a free parameter,
    gamma's free part is the vech of its leading (n-1)x(n-1) block, and eta
    drops its last column (adding-up is still imposed there).
    """
    th = []
    th += [ref[f"alpha:alpha_{i}"][0] for i in range(1, neqn + 1)]
    th += [ref[f"beta:beta_{i}"][0] for i in range(1, neqn + 1)]
    for j in range(1, neqn):
        for i in range(j, neqn):
            th.append(ref[f"gamma:gamma_{i}_{j}"][0])
    th += [ref[f"lambda:lambda_{i}"][0] for i in range(1, neqn + 1)]
    th += [ref[f"delta:delta_{i}"][0] for i in range(1, neqn + 1)]
    for v in demos:
        th += [ref[f"eta:eta_{v}_{i}"][0] for i in range(1, neqn)]
    th += [ref[f"rho:rho_{v}"][0] for v in demos]
    return np.array(th)


def main():
    ref, sc, _ = parse_coef_table(os.path.join(REPO, "log", "censor_v2.log"))
    normalize = os.environ.get('NORM', '1') == '1'
    d0 = load(normalize=normalize)
    print('shares normalised to sum to one:', normalize)
    W = d0[SHARES].to_numpy(float)
    lnp = np.log(d0[PRICES].to_numpy(float))
    lnexp = np.log(d0[EXPVAR].to_numpy(float))
    Z = d0[DEMOS].to_numpy(float)

    spec = Spec(N_G, N_R, True, True)
    th = stata_theta(ref)
    print(f"free parameters recovered from the log : {th.size} "
          f"(model expects {spec.n_free})")

    bfull_stata = np.array([ref[nm][0] for nm in spec.full_names(DEMOS)])
    bfull_mine = full_vector(th, spec)
    print("restriction round trip  max |diff| = "
          f"{np.abs(bfull_mine - bfull_stata).max():.3e}"
          "   <- the log's *derived* coefficients (gamma_.,17, gamma_17_17,"
          " eta_.,17) must come out of our restrictions")

    Dl = delta_matrix(spec)
    se_stata = np.array([ref[nm][1] for nm in spec.full_names(DEMOS)])

    for predict in ("pr", "xb"):
        fs = first_stage(W, lnp, lnexp, Z, predict=predict)
        d = DemandData(lnp=lnp, lnexp=lnexp, shares=W, demo=Z,
                       cdf=fs.cdf, pdf=fs.pdf, a0=10.0)
        u = residuals(th, d, spec)
        sigma = (u.T @ u) / d.nobs
        _, logdet = np.linalg.slogdet(sigma)
        ll = -(d.nobs * N_G / 2.0) * (1.0 + LOG2PI) - (d.nobs / 2.0) * logdet
        P = _whitener(sigma)
        obj = _objective(th, d, spec, P)
        cache = make_cache(spec)
        G, g, _ = _normal_equations(th, d, spec, cache, P, 3000)
        step = _solve_scaled(G, g)
        nrtol = abs(float(step @ g)) / max(abs(obj), 1e-300)
        V = np.linalg.inv(G)
        se = np.sqrt(np.diag(Dl @ V @ Dl.T))
        rel = np.abs(se - se_stata) / np.maximum(se_stata, 1e-12)

        print(f"\n--- first_stage_predict = {predict!r}")
        print(f"  log-likelihood at Stata's theta : {ll:.4f}")
        print(f"  Stata's reported e(ll)          : {sc['Log-likelihood']}")
        print(f"  |difference|                    : "
              f"{abs(ll - sc['Log-likelihood']):.4f}  (relative "
              f"{abs(ll - sc['Log-likelihood']) / abs(sc['Log-likelihood']):.2e})")
        print(f"  IFGNLS objective                : {obj:.10f}   "
              f"(N*m = {d.nobs * N_G})")
        print(f"  relative gradient               : {nrtol:.3e}"
              "   <- small => Stata's point is stationary for us too")
        print(f"  delta-method SEs                : max abs dev "
              f"{np.abs(se - se_stata).max():.3e}, median relative dev "
              f"{np.median(rel):.3e}")


if __name__ == "__main__":
    main()
