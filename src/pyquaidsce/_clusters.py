"""Cluster labels and whole-cluster pairs resampling."""

from numbers import Real

import numpy as np
import pandas as pd


def cluster_codes(labels, nobs):
    """Factorize labels without deleting or reordering estimation rows."""
    values = np.asarray(labels)
    if values.ndim != 1 or len(values) != nobs:
        raise ValueError("cluster labels must have one entry per estimation observation")
    if pd.isna(values).any() or any(
        isinstance(x, Real) and not np.isfinite(x) for x in values
    ):
        raise ValueError("cluster labels must not be missing or nonfinite in the estimation sample")
    try:
        codes, groups = pd.factorize(values, sort=False)
    except TypeError as exc:
        raise ValueError("cluster labels must be hashable") from exc
    if len(groups) < 2:
        raise ValueError("cluster inference requires at least two clusters")
    return codes, len(groups)


def cluster_rows(codes, n_clusters):
    """Keep all rows in each cluster, including interleaved observations."""
    order = np.argsort(codes, kind='stable')
    return np.split(order, np.cumsum(np.bincount(codes, minlength=n_clusters))[:-1])


def resample_indices(rng, nobs, groups=None):
    """Draw N rows or G whole clusters with replacement."""
    if groups is None:
        return rng.integers(0, nobs, size=nobs)
    draws = rng.integers(0, len(groups), size=len(groups))
    return np.concatenate([groups[g] for g in draws])
