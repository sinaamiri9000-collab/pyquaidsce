"""
The (censored) quadratic almost-ideal demand system.

Share equations, exactly as coded in ``_quaidsce__expshrs`` (Mata):

    ln a(p)   = a0 + sum_i alpha_i ln p_i
                   + 1/2 sum_i sum_j gamma_ij ln p_i ln p_j
    b(p)      = exp( sum_i beta_i ln p_i )
    mbar(z)   = 1 + sum_r rho_r z_r                     (Ray 1983 scaling)
    c(p, z)   = exp( sum_i ( sum_r eta_ri z_r ) ln p_i )
    D         = ln m - ln a(p) - ln mbar(z)

    w*_i      = alpha_i + sum_j gamma_ij ln p_j
                + ( beta_i + sum_r eta_ri z_r ) * D
                + lambda_i / ( b(p) c(p,z) ) * D^2

and the Shonkwiler-Yen (1999) two-step censoring transformation

    w_i       = Phi_i * w*_i + delta_i * phi_i

where ``Phi_i``/``phi_i`` come from the first-stage probit.

IMPORTANT — the ``pdf``/``cdf`` inputs
-------------------------------------
``quaidsce_c.ado`` builds them as::

    quietly predict du`i'                 // after -probit-
    qui replace pdf`i' = normalden(du`i')
    qui replace cdf`i' = normal(du`i')

Stata's ``predict`` after ``probit`` defaults to the **probability**, not the
linear predictor.  So the shipped command actually uses
``Phi_i = Phi(Phi(x'tau))`` and ``phi_i = phi(Phi(x'tau))``.  The
Shonkwiler-Yen estimator calls for ``Phi(x'tau)`` and ``phi(x'tau)``.

Both behaviours are implemented and selected by ``first_stage_predict``:

    "pr"  (default) reproduces the shipped Stata command bit for bit;
    "xb"             is the textbook Shonkwiler-Yen estimator.

The Jacobian below is analytic and exact; it is checked against a
finite-difference Jacobian in the test suite.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .params import Coefs, Spec, full_slices, unpack, vech_index


# --------------------------------------------------------------------------- #
#  Data container
# --------------------------------------------------------------------------- #
@dataclass
class DemandData:
    """Everything the share equations need, already on the estimation sample."""

    lnp: np.ndarray  # (N, n) log prices
    lnexp: np.ndarray  # (N,)  log total expenditure
    shares: np.ndarray  # (N, n) observed budget shares
    demo: np.ndarray  # (N, R) demographics (may be (N, 0))
    cdf: np.ndarray  # (N, n) Phi_i   (ones when nocensor)
    pdf: np.ndarray  # (N, n) phi_i   (zeros when nocensor)
    a0: float = 0.0
    control_function: Optional[np.ndarray] = None  # (N,) external residual

    def __post_init__(self) -> None:
        self.lnp = np.ascontiguousarray(self.lnp, dtype=float)
        self.lnexp = np.ascontiguousarray(self.lnexp, dtype=float).ravel()
        self.shares = np.ascontiguousarray(self.shares, dtype=float)
        self.demo = np.ascontiguousarray(self.demo, dtype=float)
        if self.demo.ndim == 1:
            self.demo = self.demo[:, None]
        self.cdf = np.ascontiguousarray(self.cdf, dtype=float)
        self.pdf = np.ascontiguousarray(self.pdf, dtype=float)
        if self.control_function is None:
            self.control_function = np.zeros(self.lnp.shape[0], dtype=float)
        else:
            self.control_function = np.ascontiguousarray(
                self.control_function, dtype=float
            ).ravel()
        if self.control_function.size != self.lnp.shape[0]:
            raise ValueError("control_function must have one value per observation")
        if not np.isfinite(self.control_function).all():
            raise ValueError("control_function must contain only finite values")

    @property
    def nobs(self) -> int:
        return self.lnp.shape[0]

    def subset(self, idx: np.ndarray) -> "DemandData":
        return DemandData(
            lnp=self.lnp[idx],
            lnexp=self.lnexp[idx],
            shares=self.shares[idx],
            demo=self.demo[idx],
            cdf=self.cdf[idx],
            pdf=self.pdf[idx],
            a0=self.a0,
            control_function=self.control_function[idx],
        )


# --------------------------------------------------------------------------- #
#  Intermediate quantities
# --------------------------------------------------------------------------- #
@dataclass
class _Inner:
    A: np.ndarray  # (N,)   ln a(p)
    lnb: np.ndarray  # (N,)   ln b(p)
    lnc: np.ndarray  # (N,)   ln c(p, z)
    mbar: np.ndarray  # (N,)
    D: np.ndarray  # (N,)   ln m - ln a(p) - ln mbar
    B: np.ndarray  # (N, n) beta_i + sum_r eta_ri z_r
    q: np.ndarray  # (N,)   1 / (b(p) c(p, z))
    wstar: np.ndarray  # (N, n) latent shares
    S: np.ndarray  # (N, n) B + 2 lambda_i D q
    T: np.ndarray  # (N, n) lambda_i D^2 q


def _inner(c: Coefs, d: DemandData, spec: Spec) -> _Inner:
    L, Z = d.lnp, d.demo
    n, R = spec.neqn, spec.ndemo

    A = d.a0 + L @ c.alpha + 0.5 * np.einsum("ti,ij,tj->t", L, c.gamma, L,
                                             optimize=True)
    lnb = L @ c.beta if spec.quadratic else np.zeros(d.nobs)

    if R > 0:
        Ze = Z @ c.eta  # (N, n)
        lnc = np.einsum("tk,tk->t", L, Ze) if spec.quadratic else np.zeros(d.nobs)
        mbar = 1.0 + Z @ c.rho
        B = c.beta[None, :] + Ze
    else:
        lnc = np.zeros(d.nobs)
        mbar = np.ones(d.nobs)
        B = np.repeat(c.beta[None, :], d.nobs, axis=0)

    D = d.lnexp - A - np.log(mbar)
    q = np.exp(-lnb - lnc) if spec.quadratic else np.zeros(d.nobs)

    wstar = (
        c.alpha[None, :]
        + L @ c.gamma.T
        + B * D[:, None]
    )
    if spec.quadratic:
        wstar = wstar + (q * D**2)[:, None] * c.lam[None, :]
        S = B + 2.0 * c.lam[None, :] * (q * D)[:, None]
        T = c.lam[None, :] * (q * D**2)[:, None]
    else:
        S = B
        T = np.zeros_like(B)

    return _Inner(A, lnb, lnc, mbar, D, B, q, wstar, S, T)


# --------------------------------------------------------------------------- #
#  Fitted shares
# --------------------------------------------------------------------------- #
def fitted_shares(theta: np.ndarray, d: DemandData, spec: Spec) -> np.ndarray:
    """Return the fitted values of the *estimated* equations.

    Shape is ``(N, n)`` when censoring, ``(N, n-1)`` otherwise (Stata drops the
    last equation to avoid a singular residual covariance matrix).
    """
    c = unpack(theta, spec)
    inn = _inner(c, d, spec)
    if spec.censor:
        augmented = inn.wstar + d.control_function[:, None] * c.cfcoef[None, :]
        return augmented * d.cdf + c.delta[None, :] * d.pdf
    return inn.wstar[:, : spec.neqn - 1]


def latent_shares(theta: np.ndarray, d: DemandData, spec: Spec) -> np.ndarray:
    """The uncensored (latent) shares w*, all n of them."""
    return _inner(unpack(theta, spec), d, spec).wstar


def augmented_latent_shares(
    theta: np.ndarray, d: DemandData, spec: Spec
) -> np.ndarray:
    """Latent QUAIDS shares augmented by the supplied/generated control residual."""
    c = unpack(theta, spec)
    return _inner(c, d, spec).wstar + d.control_function[:, None] * c.cfcoef


def residuals(theta: np.ndarray, d: DemandData, spec: Spec) -> np.ndarray:
    """``y - f(x, theta)`` for the estimated equations, shape (N, n_eq)."""
    m = spec.n_eq_estimated
    return d.shares[:, :m] - fitted_shares(theta, d, spec)


# --------------------------------------------------------------------------- #
#  Analytic Jacobian
# --------------------------------------------------------------------------- #
def jacobian_full(
    theta: np.ndarray,
    d: DemandData,
    spec: Spec,
    rows: Optional[slice] = None,
) -> np.ndarray:
    """d(fitted)/d(FULL parameter vector), shape (N, n_eq, n_full).

    The full vector is ordered
    ``alpha, beta, vech(gamma), [lambda], [delta], [vec(eta')], [rho]``,
    matching :func:`pyquaidsce.params.full_vector`.  Chain-rule to the free
    parameters is a right-multiplication by the (constant) delta matrix.
    """
    if rows is not None:
        d = DemandData(
            lnp=d.lnp[rows], lnexp=d.lnexp[rows], shares=d.shares[rows],
            demo=d.demo[rows], cdf=d.cdf[rows], pdf=d.pdf[rows], a0=d.a0,
            control_function=d.control_function[rows],
        )
    c = unpack(theta, spec)
    inn = _inner(c, d, spec)
    n, R = spec.neqn, spec.ndemo
    N = d.nobs
    m = spec.n_eq_estimated
    L, Z = d.lnp, d.demo
    sl = full_slices(spec)

    J = np.zeros((N, m, spec.n_full))
    S = inn.S[:, :m]  # (N, m)
    T = inn.T[:, :m]
    D = inn.D

    # ---- alpha:  d w*_i / d alpha_k = 1{i=k} - S_i L_k --------------------- #
    Ja = -S[:, :, None] * L[:, None, :]  # (N, m, n)
    for i in range(m):
        Ja[:, i, i] += 1.0
    J[:, :, sl["alpha"]] = Ja

    # ---- beta:   1{i=k} D - T_i L_k --------------------------------------- #
    Jb = -T[:, :, None] * L[:, None, :]
    for i in range(m):
        Jb[:, i, i] += D
    J[:, :, sl["beta"]] = Jb

    # ---- gamma (vech order) ----------------------------------------------- #
    vidx = vech_index(n)
    g0 = sl["gamma"].start
    for col, (k, l) in enumerate(vidx):
        if k == l:
            dA = 0.5 * L[:, k] ** 2
        else:
            dA = L[:, k] * L[:, l]
        blk = -S * dA[:, None]  # (N, m)
        if k < m:
            blk[:, k] += L[:, l]
        if l < m and l != k:
            blk[:, l] += L[:, k]
        J[:, :, g0 + col] = blk

    # ---- lambda: 1{i=k} q D^2 -------------------------------------------- #
    if spec.quadratic:
        l0 = sl["lambda"].start
        qd2 = inn.q * D**2
        for i in range(m):
            J[:, i, l0 + i] = qd2

    # ---- delta: enters only through the censoring transformation ---------- #
    #   handled below, after the cdf scaling.

    # ---- eta:  1{i=k} z_r D - T_i z_r L_k -------------------------------- #
    if R > 0:
        e0 = sl["eta"].start
        for r in range(R):
            zr = Z[:, r]
            base = -T * (zr[:, None] * 1.0)  # (N, m) placeholder, scaled below
            for k in range(n):
                col = e0 + r * n + k
                blk = base * L[:, k][:, None]
                if k < m:
                    blk[:, k] += zr * D
                J[:, :, col] = blk

        # ---- rho: -S_i z_r / mbar ---------------------------------------- #
        r0 = sl["rho"].start
        for r in range(R):
            J[:, :, r0 + r] = -S * (Z[:, r] / inn.mbar)[:, None]

    # ---- censoring transformation ---------------------------------------- #
    if spec.censor:
        J *= d.cdf[:, :m, None]
        d0 = sl["delta"].start
        for i in range(m):
            J[:, i, d0 + i] = d.pdf[:, i]
        if spec.control_function:
            c0 = sl["cfcoef"].start
            for i in range(m):
                J[:, i, c0 + i] = d.cdf[:, i] * d.control_function

    return J
