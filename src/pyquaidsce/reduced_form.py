"""Reduced-form estimation for endogenous total expenditure.

The internal ``ivexp`` path follows the augmented-regression/control-function
construction used for almost-ideal demand systems: log total expenditure is
projected on the exogenous demand-system variables, excluded instruments, and
an intercept.  The OLS residual is then available to both the participation
Probits and the latent demand equations.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
from scipy.stats import f as f_distribution


@dataclass(frozen=True)
class ExpenditureReducedForm:
    """OLS result for the internally generated expenditure control function."""

    outcome_name: str
    regressor_names: tuple[str, ...]
    instrument_names: tuple[str, ...]
    b: np.ndarray
    V: np.ndarray
    fitted: np.ndarray
    residuals: np.ndarray
    nobs: int
    rank: int
    df_resid: int
    rss: float
    r_squared: float
    adjusted_r_squared: float
    excluded_f: float
    excluded_pvalue: float
    excluded_df_num: int
    excluded_df_den: int

    @property
    def se(self) -> np.ndarray:
        diagonal = np.diag(self.V).copy()
        diagonal[diagonal < 0.0] = np.nan
        return np.sqrt(diagonal)

    def named(self) -> dict[str, float]:
        return dict(zip(self.regressor_names, self.b))


def fit_expenditure_reduced_form(
    lnexp: np.ndarray,
    lnp: np.ndarray,
    demographics: np.ndarray,
    instruments: np.ndarray,
    *,
    outcome_name: str,
    price_names: Sequence[str],
    price_inputs_are_logs: bool,
    demographic_names: Sequence[str],
    instrument_names: Sequence[str],
) -> ExpenditureReducedForm:
    """Fit the log-expenditure reduced form and return its OLS residual.

    The unrestricted design is ``[log prices, demographics, instruments, 1]``.
    A classical (homoskedastic) joint F test compares it with the restricted
    design that omits all excluded instruments.
    """

    y = np.asarray(lnexp, dtype=float).reshape(-1)
    lp = np.asarray(lnp, dtype=float)
    demo = np.asarray(demographics, dtype=float)
    iv = np.asarray(instruments, dtype=float)
    if lp.ndim != 2 or demo.ndim != 2 or iv.ndim != 2:
        raise ValueError("reduced-form inputs must be two-dimensional matrices")
    if not (lp.shape[0] == demo.shape[0] == iv.shape[0] == y.size):
        raise ValueError("reduced-form inputs must have the same number of rows")
    if iv.shape[1] == 0:
        raise ValueError("ivexp must contain at least one excluded instrument")
    if not (
        np.isfinite(y).all()
        and np.isfinite(lp).all()
        and np.isfinite(demo).all()
        and np.isfinite(iv).all()
    ):
        raise ValueError("reduced-form variables must contain only finite values")

    controls = np.column_stack([lp, demo])
    unrestricted = np.column_stack([controls, iv, np.ones(y.size)])
    restricted = np.column_stack([controls, np.ones(y.size)])
    nobs, nreg = unrestricted.shape
    if nobs <= nreg:
        raise ValueError(
            "the expenditure reduced form has no residual degrees of freedom; "
            f"need more than {nreg} complete observations"
        )

    b, _, rank, _ = np.linalg.lstsq(unrestricted, y, rcond=None)
    if int(rank) != nreg:
        raise ValueError(
            "the expenditure reduced-form design is rank deficient; check "
            "ivexp, prices, demographics, and duplicate/collinear variables"
        )
    fitted = unrestricted @ b
    residuals = y - fitted
    rss = float(residuals @ residuals)
    df_resid = nobs - nreg
    scale = max(1.0, float(y @ y))
    if rss <= np.finfo(float).eps * scale:
        raise ValueError(
            "the expenditure reduced form fits log expenditure perfectly; "
            "the generated control residual has no usable variation"
        )

    sigma2 = rss / df_resid
    V = sigma2 * np.linalg.inv(unrestricted.T @ unrestricted)
    centered = y - float(np.mean(y))
    tss = float(centered @ centered)
    r_squared = 1.0 - rss / tss if tss > 0.0 else np.nan
    adjusted_r_squared = (
        1.0 - (1.0 - r_squared) * (nobs - 1) / df_resid
        if np.isfinite(r_squared) else np.nan
    )

    b_restricted = np.linalg.lstsq(restricted, y, rcond=None)[0]
    restricted_residuals = y - restricted @ b_restricted
    restricted_rss = float(restricted_residuals @ restricted_residuals)
    q = iv.shape[1]
    numerator = max(0.0, restricted_rss - rss) / q
    excluded_f = numerator / sigma2
    excluded_pvalue = float(
        f_distribution.sf(excluded_f, q, df_resid)
    )

    price_labels = tuple(
        str(name) if price_inputs_are_logs else f"ln({name})"
        for name in price_names
    )
    regressor_names = (
        price_labels
        + tuple(str(name) for name in demographic_names)
        + tuple(str(name) for name in instrument_names)
        + ("_cons",)
    )
    return ExpenditureReducedForm(
        outcome_name=str(outcome_name),
        regressor_names=regressor_names,
        instrument_names=tuple(str(name) for name in instrument_names),
        b=np.asarray(b, dtype=float),
        V=np.asarray(V, dtype=float),
        fitted=np.asarray(fitted, dtype=float),
        residuals=np.asarray(residuals, dtype=float),
        nobs=int(nobs),
        rank=int(rank),
        df_resid=int(df_resid),
        rss=rss,
        r_squared=float(r_squared),
        adjusted_r_squared=float(adjusted_r_squared),
        excluded_f=float(excluded_f),
        excluded_pvalue=excluded_pvalue,
        excluded_df_num=int(q),
        excluded_df_den=int(df_resid),
    )
