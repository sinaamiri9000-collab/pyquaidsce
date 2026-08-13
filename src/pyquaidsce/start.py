"""
Starting values for the nonlinear system.

Stata's ``nlsur`` starts every parameter at zero unless ``initial()`` is given.
That works, but on a 17-good censored system the Gauss-Newton path from the
origin is long and runs along a very flat valley.  We therefore build a cheap
consistent starting point from a *linearised* AIDS fit and then let the
optimiser take over; the optimum it converges to is unchanged, only the number
of iterations is.

Linearisation
-------------
Replace the translog price index ``ln a(p)`` by Stone's index
``ln P* = sum_j wbar_j ln p_j`` and set ``b(p) = c(p,z) = mbar(z) = 1``.  The
share equations then become linear,

    w_i / Phi_i  ~  alpha_i + sum_j gamma_ij ln p_j + beta_i * D
                    + sum_r eta_ri z_r * D + lambda_i * D^2
                    + delta_i * phi_i / Phi_i ,      D = ln m - ln P*

so each equation can be fitted by OLS.  The resulting ``gamma`` is then made
symmetric and homogeneous, and the whole thing is mapped into the free
parameter vector with the same restrictions the model imposes.
"""

from __future__ import annotations

import numpy as np

from .model import DemandData
from .params import Spec, free_slices


def linear_start(d: DemandData, spec: Spec) -> np.ndarray:
    n, R = spec.neqn, spec.ndemo
    m = spec.n_eq_estimated
    L, Z = d.lnp, d.demo
    wbar = d.shares[:, :n].mean(axis=0)
    wbar = wbar / max(wbar.sum(), 1e-12)

    lnPs = L @ wbar
    D = d.lnexp - lnPs

    cols = [np.ones(d.nobs), *[L[:, j] for j in range(n)], D]
    if R > 0:
        cols += [Z[:, r] * D for r in range(R)]
    if spec.quadratic:
        cols.append(D**2)
    X0 = np.column_stack(cols)

    alpha = np.zeros(n)
    gamma = np.zeros((n, n))
    beta = np.zeros(n)
    eta = np.zeros((R, n))
    lam = np.zeros(n)
    delta = np.zeros(n)

    for i in range(n):
        if spec.censor:
            phi = np.clip(d.cdf[:, i], 1e-8, None)
            y = d.shares[:, i] / phi
            X = np.column_stack([X0, d.pdf[:, i] / phi])
        else:
            y = d.shares[:, i]
            X = X0
        # ridge-stabilised least squares: we only need a sane starting point
        XtX = X.T @ X
        XtX += 1e-8 * np.trace(XtX) / XtX.shape[0] * np.eye(XtX.shape[0])
        bb = np.linalg.solve(XtX, X.T @ y)
        p = 0
        alpha[i] = bb[p]; p += 1
        gamma[i, :] = bb[p:p + n]; p += n
        beta[i] = bb[p]; p += 1
        if R > 0:
            eta[:, i] = bb[p:p + R]; p += R
        if spec.quadratic:
            lam[i] = bb[p]; p += 1
        if spec.censor:
            delta[i] = bb[p]; p += 1

    # symmetry + homogeneity on gamma
    gamma = 0.5 * (gamma + gamma.T)
    gamma -= gamma.mean(axis=1, keepdims=True)
    gamma = 0.5 * (gamma + gamma.T)

    if not spec.censor:  # adding-up is imposed by the parameterisation
        alpha = alpha + (1.0 - alpha.sum()) / n
        beta = beta - beta.mean()
        lam = lam - lam.mean()
    if R > 0:
        eta = eta - eta.mean(axis=1, keepdims=True)

    # ---- pack into the free vector ---------------------------------------- #
    th = np.zeros(spec.n_free)
    sl = free_slices(spec)
    w = n if spec.censor else n - 1
    th[sl["alpha"]] = alpha[:w]
    th[sl["beta"]] = beta[:w]
    gv = [gamma[i, j] for j in range(n - 1) for i in range(j, n - 1)]
    th[sl["gamma"]] = gv
    if spec.quadratic:
        th[sl["lambda"]] = lam[:w]
    if spec.censor:
        th[sl["delta"]] = delta
    if R > 0:
        ev = []
        for r in range(R):
            ev.extend(eta[r, : n - 1])
        th[sl["eta"]] = ev
        th[sl["rho"]] = 0.0
    return th
