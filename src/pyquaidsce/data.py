"""Shared observation container for demand systems."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np


def as_matrix(data, columns):
    """Read numeric columns in their declared order, including an empty design."""
    return (np.column_stack([np.asarray(data[c], dtype=float) for c in columns])
            if columns is not None and len(columns) else np.empty((len(data), 0)))


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
