"""
Fast Jacobian with respect to the *free* parameters.

``model.jacobian_full`` differentiates with respect to the unrestricted
parameter vector and then needs a right-multiplication by the delta matrix.  For
the ``gamma`` block that product costs ``O(N m n^2 * n^2)`` -- on a 17-good
system it is the single most expensive operation in the whole estimator, more
than the normal equations themselves.

This module builds the same matrix directly in free coordinates.  The trick is
that the two pieces of the ``gamma`` derivative,

    d w*_i / d gamma_kl  =  ( 1{i=k} L_l + 1{i=l} L_k )  -  S_i * dA_kl ,

both collapse cheaply once the restriction matrix ``Gam = d vech(Gamma)/d theta``
is folded in:

* the second piece is a single row vector per observation,
  ``dA_free = dA_vech @ Gam``, shared by all equations;
* the first piece becomes ``L[t, :] @ C_i`` where
  ``C_i[l, :] = Gam[vech_index(max(i,l), min(i,l)), :]`` -- one small matrix
  product per equation.

Verified against ``model.jacobian_full @ delta_matrix`` in the test suite.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np

from .model import DemandData, _inner
from .params import (Spec, delta_matrix, free_slices, full_slices, unpack,
                     vech_index)


@dataclass
class JacCache:
    spec: Spec
    Gam: np.ndarray  # (n(n+1)/2, n(n-1)/2)
    Ci: List[np.ndarray]  # n matrices of shape (n, n(n-1)/2)
    vidx: List[tuple]
    blk_ab: np.ndarray  # alpha / beta / lambda block, (n, w)
    blk_eta: np.ndarray  # (n, n-1)
    fs: dict
    xs: dict


def make_cache(spec: Spec) -> JacCache:
    n = spec.neqn
    D = delta_matrix(spec)
    fs, xs = full_slices(spec), free_slices(spec)
    Gam = np.ascontiguousarray(D[fs["gamma"], xs["gamma"]])
    vidx = vech_index(n)
    pos = {(i, j): k for k, (i, j) in enumerate(vidx)}
    Ci = []
    for i in range(n):
        Ci.append(np.ascontiguousarray(
            np.stack([Gam[pos[(max(i, l), min(i, l))]] for l in range(n)])
        ))
    blk_ab = np.ascontiguousarray(
        np.eye(n) if spec.censor
        else np.vstack([np.eye(n - 1), -np.ones((1, n - 1))])
    )
    blk_eta = np.ascontiguousarray(
        np.vstack([np.eye(n - 1), -np.ones((1, n - 1))])
    )
    return JacCache(spec, Gam, Ci, vidx, blk_ab, blk_eta, fs, xs)


def jacobian_free(
    theta: np.ndarray,
    d: DemandData,
    spec: Spec,
    cache: JacCache,
    rows: Optional[slice] = None,
) -> np.ndarray:
    """d(fitted)/d(free parameters), shape ``(nc, n_eq, n_free)``."""
    if rows is not None:
        d = DemandData(
            lnp=d.lnp[rows], lnexp=d.lnexp[rows], shares=d.shares[rows],
            demo=d.demo[rows], cdf=d.cdf[rows], pdf=d.pdf[rows], a0=d.a0,
            control_function=d.control_function[rows],
        )
    c = unpack(theta, spec)
    inn = _inner(c, d, spec)
    n, R = spec.neqn, spec.ndemo
    m = spec.n_eq_estimated
    N = d.nobs
    L, Z = d.lnp, d.demo
    xs = cache.xs
    K = spec.n_free

    S = inn.S[:, :m]
    T = inn.T[:, :m]
    D = inn.D

    J = np.zeros((N, m, K))

    # ---- alpha ------------------------------------------------------------- #
    A_full = -S[:, :, None] * L[:, None, :]
    for i in range(m):
        A_full[:, i, i] += 1.0
    J[:, :, xs["alpha"]] = A_full @ cache.blk_ab

    # ---- beta -------------------------------------------------------------- #
    B_full = -T[:, :, None] * L[:, None, :]
    for i in range(m):
        B_full[:, i, i] += D
    J[:, :, xs["beta"]] = B_full @ cache.blk_ab

    # ---- gamma (directly in free coordinates) ------------------------------ #
    dA = np.empty((N, len(cache.vidx)))
    for col, (k, l) in enumerate(cache.vidx):
        dA[:, col] = (0.5 * L[:, k] ** 2) if k == l else (L[:, k] * L[:, l])
    dA_free = dA @ cache.Gam  # (N, n(n-1)/2)
    gsl = xs["gamma"]
    for i in range(m):
        J[:, i, gsl] = L @ cache.Ci[i]
        J[:, i, gsl] -= S[:, i][:, None] * dA_free

    # ---- lambda ------------------------------------------------------------ #
    if spec.quadratic:
        qd2 = inn.q * D**2
        Lm_full = np.zeros((N, m, n))
        for i in range(m):
            Lm_full[:, i, i] = qd2
        J[:, :, xs["lambda"]] = Lm_full @ cache.blk_ab

    # ---- eta and rho ------------------------------------------------------- #
    if R > 0:
        e_free = xs["eta"].start
        ngm1 = n - 1
        for r in range(R):
            zr = Z[:, r]
            blk = -T * zr[:, None]  # (N, m)
            E_full = blk[:, :, None] * L[:, None, :]
            for i in range(m):
                E_full[:, i, i] += zr * D
            J[:, :, e_free + r * ngm1: e_free + (r + 1) * ngm1] = (
                E_full @ cache.blk_eta
            )
        r0 = xs["rho"].start
        for r in range(R):
            J[:, :, r0 + r] = -S * (Z[:, r] / inn.mbar)[:, None]

    # ---- censoring transformation ----------------------------------------- #
    if spec.censor:
        J *= d.cdf[:, :m, None]
        d0 = xs["delta"].start
        J[:, :, d0: d0 + n] = 0.0
        for i in range(m):
            J[:, i, d0 + i] = d.pdf[:, i]
        if spec.control_function:
            c0 = xs["cfcoef"].start
            for i in range(m):
                J[:, i, c0 + i] = d.cdf[:, i] * d.control_function

    return J
