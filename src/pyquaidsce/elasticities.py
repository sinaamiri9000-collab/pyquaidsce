"""
At-means expenditure and price elasticities, ported from the ``***Elasticities***``
block of ``quaidsce_c.ado``.

Everything is evaluated at the sample means of the observed shares, log prices,
log expenditure, demographics and -- when censoring is used -- of the
first-stage ``cdf``, ``pdf`` and ``du`` variables, exactly as the ado does.

Latent (uncensored) formulas, Poi (2012) with Ray (1983) scaling
----------------------------------------------------------------
    D        = ln m - ln a(p) - ln mbar
    mu_i     = beta_i + sum_r eta_ri zbar_r                     ("betanz_i")
    e_i      = 1 + (1/w_i) [ mu_i + 2 lambda_i D / (b c) ]
    eu_ij    = -kron_ij + (1/w_i) [ gamma_ij
                 - ( mu_i + 2 lambda_i D/(b c) ) ( alpha_j + sum_k gamma_jk lnp_k )
                 - mu_j lambda_i D^2 / (b c) ]

Shonkwiler-Yen correction actually implemented by the ado
---------------------------------------------------------
    we_i     = wbar_i * cdfbar_i + delta_i * pdfbar_i
    e_i     <- 1 + (1/we_i) [ cdfbar_i (e_i - 1) wbar_i
                              + tau_{i,M} pdfbar_i ( wbar_i - delta_i dubar_i ) ]
    eu_ij   <- -kron_ij + (1/we_i) [ cdfbar_i (eu_ij + kron_ij) wbar_i
                              + tau_{i,j} pdfbar_i ( wbar_i - delta_i dubar_i ) ]
    ec_ij    = eu_ij + e_i * wbar_j                              (Slutsky)

Known deviations in the Stata source, reproduced only when ``strict_stata``
--------------------------------------------------------------------------
1. With **no demographics** and the quadratic term, the last term of the
   uncompensated elasticity uses ``beta_i * lambda_i`` where Poi (2012) has
   ``beta_j * lambda_i``.  ``strict_stata=True`` reproduces the ado.
2. With **demographics but noquadratic**, the ado assigns the expenditure
   elasticity to a *global* macro and then reads an (empty) *local* one.  The
   censoring correction therefore silently uses a zero latent expenditure
   elasticity. ``strict_stata=True`` reproduces that result;
   ``strict_stata=False`` uses the intended ``1 + mu_i / w_i``.
3. The reported ``ELAS_UNCOMP`` / ``ELAS_COMP`` vectors are stored in
   ``i``-major order but *labelled* in ``j``-major order, so ``e(b)``'s
   ``e_a_b`` actually holds ``eu_{b,a}``.  ``as_stata_vector()`` reproduces the
   stored order; the matrices returned here use the natural
   ``[good, price]`` convention.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional
import warnings

import numpy as np
from scipy.stats import norm

from .model import DemandData, _inner
from .params import Coefs, Spec, unpack
from .selection import FirstStageLayout


@dataclass
class Means:
    w: np.ndarray  # (n,) mean observed shares
    lnp: np.ndarray  # (n,)
    lnexp: float
    demo: np.ndarray  # (R,)
    cdf: np.ndarray  # (n,)
    pdf: np.ndarray  # (n,)
    du: np.ndarray  # (n,)
    control_function: float = 0.0


def sample_means(d: DemandData, du: Optional[np.ndarray], spec: Spec) -> Means:
    n = spec.neqn
    return Means(
        w=d.shares[:, :n].mean(axis=0),
        lnp=d.lnp.mean(axis=0),
        lnexp=float(d.lnexp.mean()),
        demo=d.demo.mean(axis=0) if d.demo.shape[1] else np.zeros(0),
        cdf=d.cdf.mean(axis=0),
        pdf=d.pdf.mean(axis=0),
        du=(du.mean(axis=0) if du is not None else np.ones(n)),
        control_function=float(d.control_function.mean()),
    )


@dataclass
class Elasticities:
    income: np.ndarray  # (n,)
    uncompensated: np.ndarray  # (n, n) [good, price]
    compensated: np.ndarray  # (n, n) [good, price]
    income_latent: np.ndarray
    uncompensated_latent: np.ndarray
    we: np.ndarray  # (n,) censoring-adjusted mean share

    # ---- Stata storage order (i-major), for e(b) comparability ------------ #
    def as_stata_vector(self) -> np.ndarray:
        return np.concatenate(
            [
                self.income,
                self.uncompensated.reshape(-1),  # row-major == i-major
                self.compensated.reshape(-1),
            ]
        )


def fitted_share_derivatives(
    theta: np.ndarray,
    d: DemandData,
    spec: Spec,
    *,
    tau: np.ndarray,
    layout: FirstStageLayout,
    selection_index: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Observation-level derivatives of fitted shares.

    Returns ``(dmu_dlnexp, dmu_dlnprices)`` with shapes ``(N, n)`` and
    ``(N, n, n)``.  The control residual is conditional and remains
    fixed.  ``selection_index`` is the Probit linear predictor ``xb``; this
    helper deliberately does not support the legacy double-CDF ``pr`` path.
    """
    if not spec.censor:
        raise ValueError("fitted-share censoring derivatives require censor=True")
    tau = np.asarray(tau, dtype=float).ravel()
    k = np.asarray(selection_index, dtype=float)
    n = spec.neqn
    if k.shape != (d.nobs, n):
        raise ValueError("selection_index must have shape (nobs, neqn)")
    if tau.size != n * layout.width:
        raise ValueError("tau length does not match the first-stage layout")
    if (not np.allclose(d.cdf, norm.cdf(k), atol=1e-11, rtol=1e-11)
            or not np.allclose(d.pdf, norm.pdf(k), atol=1e-11, rtol=1e-11)):
        raise ValueError(
            "DemandData cdf/pdf must be rebuilt from the supplied xb selection_index"
        )

    c = unpack(theta, spec)
    inn = _inner(c, d, spec)
    g = c.alpha[None, :] + d.lnp @ c.gamma.T
    dw_price = (
        c.gamma[None, :, :]
        - inn.S[:, :, None] * g[:, None, :]
        - inn.T[:, :, None] * inn.B[:, None, :]
    )
    augmented = inn.wstar + d.control_function[:, None] * c.cfcoef[None, :]
    selection_term = augmented - c.delta[None, :] * k

    tau_m = np.array([
        layout.coefficient(tau, i, layout.expenditure_position)
        for i in range(n)
    ])
    tau_p = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            tau_p[i, j] = layout.coefficient(
                tau, i, layout.price_position(j)
            )

    d_exp = (
        d.cdf * inn.S
        + d.pdf * tau_m[None, :] * selection_term
    )
    d_price = (
        d.cdf[:, :, None] * dw_price
        + d.pdf[:, :, None] * tau_p[None, :, :]
        * selection_term[:, :, None]
    )
    return d_exp, d_price


def elasticities(
    c: Coefs,
    spec: Spec,
    means: Means,
    a0: float,
    tau: Optional[np.ndarray] = None,
    np_prob: Optional[int] = None,
    strict_stata: bool = False,
    *,
    layout: Optional[FirstStageLayout] = None,
) -> Elasticities:
    """Compute the three elasticity matrices at the sample means."""
    n, R = spec.neqn, spec.ndemo
    L = means.lnp
    wbar = means.w

    # ---- price index and auxiliaries (ado lines 395-454) ------------------ #
    lnpindex = a0 + float(c.alpha @ L) + 0.5 * float(L @ c.gamma @ L)
    gsum = c.gamma @ L  # (n,)  sum_l gamma_jl lnp_l

    lnb = float(c.beta @ L) if spec.quadratic else 0.0  # ado's "bofp"
    if R > 0:
        betanz = c.beta + means.demo @ c.eta  # (n,)
        lnc = float((means.demo @ c.eta) @ L)  # ado's "cofp"
        mbar = 1.0 + float(c.rho @ means.demo)
    else:
        betanz = c.beta.copy()
        lnc = 0.0
        mbar = 1.0

    D = means.lnexp - lnpindex - np.log(mbar)
    q = np.exp(-lnb - lnc)  # 1 / (b(p) c(p, z))

    # ---- predicted shares (ado lines 458-466) ----------------------------- #
    augmented_mean = wbar + c.cfcoef * means.control_function
    if spec.censor:
        we = augmented_mean * means.cdf + c.delta * means.pdf
        cf_denominator = (
            spec.control_function
            or (layout is not None and layout.selection_cf_position is not None)
        )
        if (cf_denominator
                and (not np.isfinite(we).all() or np.any(np.abs(we) <= 1e-12))):
            raise ValueError(
                "censoring-adjusted mean share is nonfinite or too close to zero; "
                "elasticities are undefined"
            )
        if cf_denominator and np.any(we < 0):
            warnings.warn(
                "negative censoring-adjusted mean share(s); elasticities may be "
                "economically difficult to interpret",
                RuntimeWarning,
                stacklevel=2,
            )
    else:
        we = wbar.copy()

    # ---- expenditure elasticities ---------------------------------------- #
    ie_latent = np.empty(n)
    for i in range(n):
        if R == 0:
            val = 1.0 + c.beta[i] / wbar[i]
            if spec.quadratic:
                val = 1.0 + (c.beta[i] + 2.0 * c.lam[i] * q * D) / wbar[i]
        else:
            val = 1.0 + betanz[i] / wbar[i]
            if spec.quadratic:
                val = 1.0 + (betanz[i] + 2.0 * c.lam[i] * q * D) / wbar[i]
        ie_latent[i] = val

    # ------------------------------------------------------------------ #
    # Bug reproduction: in the demographics branch the ado writes the latent
    # expenditure elasticity to a *global* macro,
    #     global ie`i' = 1+`betanz`i''/`w_`i''m
    # and with `noquadratic` the `local` of the same name is never set.  The
    # censoring correction below then expands `(`ie`i''-1)` to the literal
    # "(-1)", i.e. it silently uses a latent elasticity of ZERO.  Stata does not
    # error, it just returns a wrong number -- confirmed against RUN 4 of
    # bench/small4.log to 1e-15.
    # ------------------------------------------------------------------ #
    ie_used = ie_latent
    cf_extension = (
        spec.control_function
        or (layout is not None and layout.selection_cf_position is not None)
    )
    if (strict_stata and R > 0 and not spec.quadratic and spec.censor
            and not cf_extension):
        ie_used = np.zeros(n)

    ie = ie_latent.copy()
    if spec.censor:
        if tau is None or np_prob is None:
            raise ValueError("censored elasticities need tau and np_prob")
        if layout is not None and layout.width != np_prob:
            raise ValueError("first-stage layout width does not match np_prob")
        if tau.size != n * np_prob:
            raise ValueError("tau length does not match the first-stage layout")
        for i in range(n):
            if layout is None:
                exp_pos = n if np_prob == n + R + 2 else None
                tau_m = 0.0 if exp_pos is None else tau[np_prob * i + exp_pos]
            else:
                tau_m = layout.coefficient(
                    tau, i, layout.expenditure_position
                )
            ie[i] = 1.0 + (
                means.cdf[i] * (ie_used[i] - 1.0) * wbar[i]
                + tau_m * means.pdf[i]
                * (augmented_mean[i] - c.delta[i] * means.du[i])
            ) / we[i]
    else:
        ie = ie_used.copy() if ie_used is not ie_latent else ie_latent.copy()

    # ---- uncompensated (Marshallian) ------------------------------------- #
    ue_latent = np.empty((n, n))
    for i in range(n):
        for j in range(n):
            kron = 1.0 if i == j else 0.0
            if R == 0:
                val = -kron + (
                    c.gamma[i, j] - c.beta[i] * (c.alpha[j] + gsum[j])
                ) / wbar[i]
                if spec.quadratic:
                    b_last = c.beta[i] if strict_stata else c.beta[j]
                    val = -kron + (
                        c.gamma[i, j]
                        - (c.beta[i] + 2.0 * c.lam[i] * q * D)
                        * (c.alpha[j] + gsum[j])
                        - b_last * c.lam[i] * q * D**2
                    ) / wbar[i]
            else:
                val = -kron + (
                    c.gamma[i, j] - betanz[i] * (c.alpha[j] + gsum[j])
                ) / wbar[i]
                if spec.quadratic:
                    val = -kron + (
                        c.gamma[i, j]
                        - (betanz[i] + 2.0 * c.lam[i] * q * D)
                        * (c.alpha[j] + gsum[j])
                        - betanz[j] * c.lam[i] * q * D**2
                    ) / wbar[i]
            ue_latent[i, j] = val

    ue = ue_latent.copy()
    if spec.censor:
        for i in range(n):
            for j in range(n):
                kron = 1.0 if i == j else 0.0
                if layout is None:
                    tau_p = tau[np_prob * i + j]
                else:
                    tau_p = layout.coefficient(
                        tau, i, layout.price_position(j)
                    )
                ue[i, j] = -kron + (
                    means.cdf[i] * (ue_latent[i, j] + kron) * wbar[i]
                    + tau_p * means.pdf[i]
                    * (augmented_mean[i] - c.delta[i] * means.du[i])
                ) / we[i]

    # ---- compensated (Hicksian), Slutsky --------------------------------- #
    ce = ue + np.outer(ie, wbar)

    return Elasticities(
        income=ie,
        uncompensated=ue,
        compensated=ce,
        income_latent=ie_latent,
        uncompensated_latent=ue_latent,
        we=we,
    )
