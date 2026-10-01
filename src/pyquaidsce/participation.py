"""Shared Shonkwiler--Yen participation estimation for demand models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

import numpy as np
from scipy.stats import norm

from ._timing import check_deadline
from .probit import ProbitResult, probit
from .selection import FirstStageLayout, legacy_layout


@dataclass
class FirstStage:
    tau: np.ndarray  # (n * np_prob,) stacked probit coefficients
    setau: np.ndarray  # (n*np_prob, n*np_prob) block diagonal
    cdf: np.ndarray  # (N, n)
    pdf: np.ndarray  # (N, n)
    du: np.ndarray  # (N, n) first-stage linear index X'tau
    np_prob: int
    results: List[ProbitResult]
    layout: FirstStageLayout


def first_stage(
    shares: np.ndarray,
    lnp: np.ndarray,
    lnexp: np.ndarray,
    demo: np.ndarray,
    include_lnexp: bool = True,
    *,
    design: Optional[np.ndarray] = None,
    layout: Optional[FirstStageLayout] = None,
    deadline: Optional[float] = None,
) -> FirstStage:
    """Shonkwiler-Yen step 1: a probit per share.

    The censoring correction always uses the Probit linear index ``X' tau``:
    ``cdf = Phi(X' tau)`` and ``pdf = phi(X' tau)``. This is the textbook
    Shonkwiler-Yen transformation used throughout pyquaidsce >= 1.6.0.

    ``include_lnexp=False`` reproduces the behaviour of the ado when the user
    supplies ``lnexpenditure()`` instead of ``expenditure()``: the local macro
    holding the log-expenditure temp variable is empty in that branch, so log
    expenditure silently drops out of the first-stage probits.
    """
    check_deadline(deadline)
    N, n = shares.shape
    if design is None:
        Z = [lnp]
        if include_lnexp:
            Z.append(lnexp[:, None])
        if demo.shape[1]:
            Z.append(demo)
        X = np.hstack(Z) if Z else np.zeros((N, 0))
        if layout is None:
            layout = legacy_layout(
                [f"p{j + 1}" for j in range(lnp.shape[1])],
                [f"z{r + 1}" for r in range(demo.shape[1])],
                include_expenditure=include_lnexp,
            )
    else:
        X = np.asarray(design, dtype=float)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[0] != N:
            raise ValueError("first-stage design must have one row per observation")
        if layout is None:
            raise ValueError("an explicit first-stage design requires a layout")
    assert layout is not None
    if X.shape[1] != layout.constant_position:
        raise ValueError("first-stage design and layout have inconsistent widths")
    np_prob = layout.width

    tau = np.zeros(n * np_prob)
    setau = np.zeros((n * np_prob, n * np_prob))
    cdf = np.ones((N, n))
    pdf = np.zeros((N, n))
    du = np.ones((N, n))
    res: List[ProbitResult] = []

    for i in range(n):
        check_deadline(deadline)
        w = shares[:, i]
        if w.min() > 0:
            raise ValueError(
                f"no censoring for share {i + 1} found "
                "(the S&Y first stage requires participation variation)"
            )
        z = (w > 0).astype(float)
        if z.max() == z.min():
            raise ValueError(
                f"participation indicator for share {i + 1} has no variation"
            )
        pr = probit(z, X, add_constant=True, deadline=deadline)
        if not pr.converged:
            raise RuntimeError(f"first-stage probit {i + 1} did not converge")
        if (layout.selection_cf_position is not None
                and layout.selection_cf_position in pr.dropped):
            raise ValueError(
                "selection_control_function is collinear with the first-stage "
                f"design in equation {i + 1}; its coefficient was dropped"
            )
        res.append(pr)
        sl = slice(i * np_prob, (i + 1) * np_prob)
        tau[sl] = pr.b
        setau[sl, sl] = pr.V
        xb = pr.xb(X)
        du[:, i] = xb
        pdf[:, i] = norm.pdf(xb)
        cdf[:, i] = norm.cdf(xb)

    return FirstStage(tau, setau, cdf, pdf, du, np_prob, res, layout)

