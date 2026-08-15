"""
Stata-compatible probit, used for the Shonkwiler-Yen first stage.

Reproduces ``probit y x1 ... xk`` as called in ``quaidsce_c.ado``:

* maximum likelihood by Newton-Raphson on the exact log-likelihood
  ``sum ln Phi(q_i x_i'b)`` with ``q_i = 2 y_i - 1``;
* coefficient ordering = regressors in the order given, then ``_cons``
  (this matters, because ``quaidsce_c.ado`` indexes into ``tau`` positionally);
* ``e(V)`` = inverse of the *observed* information matrix, which is Stata's
  default ``vce(oim)`` for ``probit``;
* collinear regressors are dropped the way Stata's ``_rmcoll`` does, and the
  corresponding coefficients/variances are reported as 0.

No external dependency beyond numpy/scipy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np
from scipy.stats import norm

from ._timing import check_deadline


@dataclass
class ProbitResult:
    b: np.ndarray  # (k+1,) coefficients, _cons last
    V: np.ndarray  # (k+1, k+1)
    llf: float
    n_obs: int
    n_iter: int
    converged: bool
    dropped: List[int] = field(default_factory=list)  # indices of collinear cols

    def xb(self, X: np.ndarray) -> np.ndarray:
        """Linear predictor for a design matrix *without* the constant."""
        return X @ self.b[:-1] + self.b[-1]


def _log_std_normal_cdf(z: np.ndarray) -> np.ndarray:
    return norm.logcdf(z)


def _lambda_ratio(z: np.ndarray) -> np.ndarray:
    """phi(z)/Phi(z), computed stably for very negative z."""
    out = np.empty_like(z)
    hi = z > -30.0
    out[hi] = np.exp(norm.logpdf(z[hi]) - norm.logcdf(z[hi]))
    # asymptotic expansion phi/Phi ~ -z - 1/z + ... for z -> -inf
    lo = ~hi
    if np.any(lo):
        zl = z[lo]
        out[lo] = -zl - 1.0 / zl + 2.0 / zl**3
    return out


def _drop_collinear(X: np.ndarray, tol: float = 1e-11) -> List[int]:
    """Indices of columns to drop, sweeping left to right like Stata does."""
    keep: List[int] = []
    dropped: List[int] = []
    for j in range(X.shape[1]):
        cand = keep + [j]
        A = X[:, cand]
        # scale-free rank test
        s = np.linalg.svd(A / (np.linalg.norm(A, axis=0, keepdims=True) + 1e-300),
                          compute_uv=False)
        if s[-1] > tol * s[0]:
            keep.append(j)
        else:
            dropped.append(j)
    return dropped


def probit(
    y: np.ndarray,
    X: np.ndarray,
    add_constant: bool = True,
    tol: float = 1e-12,
    max_iter: int = 200,
    deadline: Optional[float] = None,
) -> ProbitResult:
    """Fit a probit model.

    Parameters
    ----------
    y : (N,) 0/1 outcome.
    X : (N, k) regressors, *without* the constant.
    add_constant : append a column of ones (Stata always does unless -nocons-).
    """
    y = np.asarray(y, dtype=float).ravel()
    check_deadline(deadline)
    X = np.asarray(X, dtype=float)
    if X.ndim == 1:
        X = X[:, None]
    n = y.size

    Z = np.hstack([X, np.ones((n, 1))]) if add_constant else X.copy()
    k_all = Z.shape[1]

    dropped = _drop_collinear(Z)
    keep = [j for j in range(k_all) if j not in dropped]
    Zk = Z[:, keep]

    q = 2.0 * y - 1.0
    b = np.zeros(Zk.shape[1])
    # A single OLS-on-scaled-y start is what makes Stata's own start; zero is
    # equally fine since the probit log-likelihood is globally concave.
    llf_old = -np.inf
    converged = False
    it = 0
    for it in range(1, max_iter + 1):
        check_deadline(deadline)
        z = q * (Zk @ b)
        llf = float(np.sum(_log_std_normal_cdf(z)))
        lam = _lambda_ratio(z)
        g = Zk.T @ (q * lam)
        w = lam * (lam + z)
        H = -(Zk * w[:, None]).T @ Zk
        try:
            step = np.linalg.solve(H, -g)
        except np.linalg.LinAlgError:
            step = -np.linalg.lstsq(H, g, rcond=None)[0]
        # Newton with simple backtracking (guards against overshoot in the
        # first couple of steps when a regressor nearly separates the data).
        t = 1.0
        for _ in range(40):
            check_deadline(deadline)
            bt = b + t * step
            zt = q * (Zk @ bt)
            llt = float(np.sum(_log_std_normal_cdf(zt)))
            if np.isfinite(llt) and llt >= llf - 1e-14:
                break
            t *= 0.5
        b_new = b + t * step
        rel = np.max(np.abs(b_new - b) / (np.abs(b) + 1e-8))
        b = b_new
        llf_new = float(np.sum(_log_std_normal_cdf(q * (Zk @ b))))
        if abs(llf_new - llf_old) <= tol * (abs(llf_new) + 1.0) and rel < 1e-10:
            llf_old = llf_new
            converged = True
            break
        llf_old = llf_new

    z = q * (Zk @ b)
    lam = _lambda_ratio(z)
    w = lam * (lam + z)
    H = -(Zk * w[:, None]).T @ Zk
    Vk = np.linalg.inv(-H)

    b_full = np.zeros(k_all)
    V_full = np.zeros((k_all, k_all))
    b_full[keep] = b
    V_full[np.ix_(keep, keep)] = Vk

    return ProbitResult(
        b=b_full,
        V=V_full,
        llf=float(np.sum(_log_std_normal_cdf(z))),
        n_obs=n,
        n_iter=it,
        converged=converged,
        dropped=dropped,
    )
