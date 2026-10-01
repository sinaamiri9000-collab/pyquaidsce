"""Basic indirect translog demand with Pollak--Wales demographic translation.

For c = Nu z, b = 1 - p'c/m and r_j = log(p_j/(m b)),

    s_i = (alpha_i + Gamma_i r) / (1 + 1'Gamma r)
    w_i = p_i c_i/m + b s_i.

Uncensored models impose sum(alpha)=1 and Gamma symmetry. S&Y models free all
alphas by default, or impose sum(alpha)=1 with latent_adding=True. Both retain
denominator constant 1 and estimate every equation. There is no independent
committed-quantity intercept (the generalized translog model is not included).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Sequence

import numpy as np

from .data import DemandData
from .elasticities import Elasticities
from .censoring import (corrected_mean, corrected_derivatives,
                        index_coefficients, normal_terms)


@dataclass(frozen=True)
class TranslogSpec:
    neqn: int
    ndemo: int = 0
    censor: bool = False
    latent_adding: bool = False

    def __post_init__(self):
        if not isinstance(self.latent_adding, (bool, np.bool_)):
            raise ValueError("latent_adding must be a boolean")
        if self.latent_adding and not self.censor:
            raise ValueError("latent_adding=True requires censor=True; uncensored adding-up is mandatory")
        if self.neqn < 2 or self.ndemo < 0:
            raise ValueError("translog requires at least two goods and ndemo >= 0")

    @property
    def latent_adding_up(self):
        return not self.censor or self.latent_adding

    @property
    def n_alpha(self):
        return self.neqn - int(self.latent_adding_up)

    @property
    def n_eq_estimated(self):
        return self.neqn if self.censor else self.neqn - 1

    @property
    def n_gamma(self):
        return self.neqn * (self.neqn + 1) // 2

    @property
    def n_free(self):
        return self.n_alpha + self.n_gamma + self.neqn * self.ndemo + (self.neqn if self.censor else 0)

    @property
    def n_full(self):
        return self.n_free + int(self.latent_adding_up)

    def full_names(self, demo_names: Optional[Sequence[str]] = None):
        n = self.neqn
        demos = list(demo_names) if demo_names is not None else [f"z{j+1}" for j in range(self.ndemo)]
        if len(demos) != self.ndemo:
            raise ValueError("demographic names do not match the model dimensions")
        names = [f"alpha:alpha_{i+1}" for i in range(n)]
        names += [f"Gamma:gamma_{i+1}_{j+1}" for i, j in gamma_indices(n)]
        names += [f"Nu:nu_{z}_{i+1}" for z in demos for i in range(n)]
        if self.censor:
            names += [f"delta:delta_{i+1}" for i in range(n)]
        return names

    def elas_names(self):
        n = self.neqn
        names = [f"ELAS_INC:e_{i+1}" for i in range(n)]
        for prefix in ("ELAS_UNCOMP", "ELAS_COMP"):
            names += [f"{prefix}:e_{i+1}_{j+1}" for i in range(n) for j in range(n)]
        return names


@dataclass
class TranslogCoefs:
    alpha: np.ndarray
    gamma: np.ndarray
    nu: np.ndarray  # (good, demographic)
    delta: Optional[np.ndarray] = None


def gamma_indices(n):
    return [(i, j) for i in range(n) for j in range(i, n)]


def unpack(theta, spec):
    theta = np.asarray(theta)
    if theta.ndim != 1 or theta.size != spec.n_free:
        raise ValueError(f"theta must contain {spec.n_free} free translog parameters")
    n = spec.neqn
    a = spec.n_alpha
    alpha = np.empty(n, dtype=theta.dtype)
    alpha[:a] = theta[:a]
    if spec.latent_adding_up:
        alpha[-1] = 1 - alpha[:-1].sum()
    gamma = np.zeros((n, n), dtype=theta.dtype)
    for k, (i, j) in enumerate(gamma_indices(n)):
        gamma[i, j] = gamma[j, i] = theta[a+k]
    offset = a + spec.n_gamma
    end = offset + n * spec.ndemo
    nu = theta[offset:end].reshape(spec.ndemo, n).T.copy()
    delta = theta[end:].copy() if spec.censor else np.zeros(n, dtype=theta.dtype)
    return TranslogCoefs(alpha, gamma, nu, delta)


def pack(coefs, spec):
    if spec.latent_adding_up and not np.allclose(np.sum(coefs.alpha), 1):
        raise ValueError("translog alpha coefficients must sum to one")
    if not np.allclose(coefs.gamma, coefs.gamma.T):
        raise ValueError("translog Gamma must be symmetric")
    parts = [coefs.alpha[:spec.n_alpha],
                           [coefs.gamma[i, j] for i, j in gamma_indices(spec.neqn)],
                           coefs.nu.T.ravel()]
    if spec.censor:
        if coefs.delta is None or np.asarray(coefs.delta).shape != (spec.neqn,):
            raise ValueError("censored translog requires one delta coefficient per good")
        parts.append(coefs.delta)
    return np.concatenate(parts)


def full_vector(theta, spec):
    c = unpack(theta, spec)
    parts = [c.alpha,
                           [c.gamma[i, j] for i, j in gamma_indices(spec.neqn)],
                           c.nu.T.ravel()]
    if spec.censor:
        parts.append(c.delta)
    return np.concatenate(parts)


def delta_matrix(spec):
    if not spec.latent_adding_up:
        return np.eye(spec.n_free)
    out = np.zeros((spec.n_full, spec.n_free))
    n = spec.neqn
    out[:n-1, :n-1] = np.eye(n-1)
    out[n-1, :n-1] = -1
    out[n:, n-1:] = np.eye(spec.n_free - n + 1)
    return out


@dataclass
class _Inner:
    ratio: np.ndarray
    translated: np.ndarray
    b: np.ndarray
    r: np.ndarray
    denominator: np.ndarray
    s: np.ndarray
    w: np.ndarray
    g: np.ndarray
    h: float


def _inner(theta, data, spec):
    c = unpack(theta, spec)
    ratio = np.exp(data.lnp - data.lnexp[:, None])
    translated = ratio * (data.demo @ c.nu.T)
    b = 1 - translated.sum(axis=1)
    if np.any(np.real(b) <= 0) or not np.isfinite(b).all():
        raise ValueError("demographic translation requires positive effective expenditure")
    r = data.lnp - data.lnexp[:, None] - np.log(b)[:, None]
    q = r @ c.gamma.T
    denominator = 1 + q.sum(axis=1)
    if np.any(np.abs(denominator) < 1e-12) or not np.isfinite(denominator).all():
        raise ValueError("translog share denominator is zero or nonfinite")
    s = (c.alpha + q) / denominator[:, None]
    w = translated + b[:, None] * s
    if not np.isfinite(w).all():
        raise ValueError("nonfinite translog shares")
    g = c.gamma.sum(axis=1)
    return _Inner(ratio, translated, b, r, denominator, s, w, g, g.sum())


def latent_shares(theta, data, spec):
    """All goods before censoring, without post-hoc normalization."""
    return _inner(theta, data, spec).w


def predicted_shares(theta, data, spec):
    """Every predicted observed share, using the same mean as estimation."""
    latent = latent_shares(theta, data, spec)
    if spec.censor:
        return corrected_mean(latent, data.cdf, data.pdf, unpack(theta, spec).delta)
    return latent


def fitted_shares(theta, data, spec):
    return predicted_shares(theta, data, spec)[:, :spec.n_eq_estimated]


def jacobian_free(theta, data, spec, rows=None):
    """Analytic derivatives of estimated shares with respect to free parameters."""
    if rows is not None:
        data = data.subset(rows)
    inn = _inner(theta, data, spec)
    n, m = spec.neqn, spec.n_eq_estimated
    out = np.zeros((data.nobs, m, spec.n_free), dtype=np.result_type(np.asarray(theta).dtype, float))
    scale = inn.b / inn.denominator
    a = spec.n_alpha
    alpha_map = np.eye(n)[:, :a]
    if spec.latent_adding_up:
        alpha_map[-1, :] = -1
    out[:, :, :a] = scale[:, None, None] * alpha_map[None, :m, :]
    s = inn.s[:, :m]
    for k, (i, j) in enumerate(gamma_indices(n)):
        term = np.zeros((data.nobs, m), dtype=out.dtype)
        if i < m:
            term[:, i] += inn.r[:, j]
        total = inn.r[:, j].copy()
        if i != j:
            if j < m:
                term[:, j] += inn.r[:, i]
            total += inn.r[:, i]
        out[:, :, a+k] = scale[:, None] * (term - s * total[:, None])
    effect = (inn.g[:m] - s * inn.h) / inn.denominator[:, None] - s
    base = np.eye(n)[:m][None, :, :] + effect[:, :, None]
    offset = a+spec.n_gamma
    for r in range(spec.ndemo):
        out[:, :, offset+r*n:offset+(r+1)*n] = (
            base * inn.ratio[:, None, :] * data.demo[:, r, None, None]
        )
    if spec.censor:
        out *= data.cdf[:, :, None]
        out[:, :, -n:] = data.pdf[:, :, None] * np.eye(n)[None, :, :]
    return out


def share_derivatives(theta, data, spec):
    """Derivatives of every share with respect to log expenditure and log prices."""
    inn = _inner(theta, data, spec)
    c = unpack(theta, spec)
    s, t, b, den = inn.s, inn.translated, inn.b, inn.denominator
    extra = (inn.g - s*inn.h) / den[:, None]
    dm = -t + (1-b)[:, None]*s - extra
    dp = np.eye(spec.neqn)[None] * t[:, :, None]
    dp = dp - s[:, :, None]*t[:, None, :]
    dp += (b/den)[:, None, None] * (c.gamma[None] - s[:, :, None]*inn.g[None, None, :])
    dp += extra[:, :, None]*t[:, None, :]
    return dm, dp


def fitted_share_derivatives(theta, data, spec, *, tau, layout, selection_index):
    """Derivatives of the S&Y observed mean, including participation changes."""
    if not spec.censor:
        return share_derivatives(theta, data, spec)
    k = np.asarray(selection_index)
    if k.shape != (data.nobs, spec.neqn):
        raise ValueError("selection_index must have shape (nobs, neqn)")
    cdf, pdf = normal_terms(k)
    if (not np.allclose(data.cdf, cdf, atol=1e-11, rtol=1e-11)
            or not np.allclose(data.pdf, pdf, atol=1e-11, rtol=1e-11)):
        raise ValueError("DemandData cdf/pdf must match the supplied xb selection_index")
    latent = latent_shares(theta, data, spec)
    dm, dp = share_derivatives(theta, data, spec)
    km, kp = index_coefficients(tau, layout, spec.neqn)
    return corrected_derivatives(latent, dm, dp, k, unpack(theta, spec).delta, km, kp)


def elasticities(theta, data, spec, *, tau=None, layout=None, selection_index=None):
    """Predicted-share expenditure, Marshallian and Hicksian elasticities.

    data may contain one evaluation point or multiple observations. Single-point
    results use the same Elasticities container as AIDS/QUAIDS.
    """
    latent = latent_shares(theta, data, spec)
    dm, dp = share_derivatives(theta, data, spec)
    with np.errstate(divide="ignore", invalid="ignore"):
        latent_income = 1 + dm/latent
        latent_uncomp = dp/latent[:, :, None] - np.eye(spec.neqn)[None]
    if spec.censor:
        if tau is None or layout is None or selection_index is None:
            raise ValueError("S&Y elasticities require tau, layout and the xb selection_index")
        k = np.asarray(selection_index)
        if k.shape != latent.shape:
            raise ValueError("selection_index must have shape (nobs, neqn)")
        cdf, pdf = normal_terms(k)
        w = corrected_mean(latent, cdf, pdf, unpack(theta, spec).delta)
        km, kp = index_coefficients(tau, layout, spec.neqn)
        dm, dp = corrected_derivatives(latent, dm, dp, k, unpack(theta, spec).delta, km, kp)
    else:
        w = latent
    if np.any(np.abs(w) < 1e-12):
        raise ValueError("elasticities are undefined at zero predicted shares")
    income = 1 + dm/w
    uncomp = dp/w[:, :, None] - np.eye(spec.neqn)[None]
    comp = uncomp + income[:, :, None]*w[:, None, :]
    if data.nobs == 1:
        return Elasticities(income[0], uncomp[0], comp[0], latent_income[0],
                            latent_uncomp[0], w[0])
    return Elasticities(income, uncomp, comp, latent_income, latent_uncomp, w)


class TranslogBackend:
    start_methods = ("zero", "mean")

    def make_cache(self, spec):
        return None

    def fitted(self, theta, data, spec):
        return fitted_shares(theta, data, spec)

    def jacobian(self, theta, data, spec, cache, rows):
        return jacobian_free(theta, data, spec, rows)

    def initial(self, data, spec, start):
        theta = np.zeros(spec.n_free)
        if start == "mean":
            theta[:spec.n_alpha] = data.shares[:, :spec.n_alpha].mean(axis=0)
            if spec.censor:
                theta[:spec.n_alpha] /= np.maximum(data.cdf[:, :spec.n_alpha].mean(axis=0), 1e-6)
        return theta


TRANSLOG_BACKEND = TranslogBackend()
