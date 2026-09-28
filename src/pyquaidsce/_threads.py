"""Temporary BLAS thread limits for pyquaidsce estimators."""

from __future__ import annotations

import inspect
from contextlib import contextmanager
from functools import wraps
from typing import Optional

from threadpoolctl import threadpool_limits


def _validate_blas_threads(n: Optional[int]) -> Optional[int]:
    """Validate and normalize a requested BLAS thread count."""
    if n is None:
        return None
    if isinstance(n, bool):
        raise ValueError("blas_threads must be a positive integer or None")
    try:
        value = int(n)
        exact = float(n) == float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError("blas_threads must be a positive integer or None") from exc
    if not exact or value < 1:
        raise ValueError("blas_threads must be a positive integer or None")
    return value


@contextmanager
def blas_thread_limit(n: Optional[int] = 1):
    """Temporarily limit BLAS threads and restore the prior setting on exit.

    ``None`` leaves the BLAS runtime untouched. The limit is scoped to BLAS so
    unrelated OpenMP thread pools in the caller's process are not modified.
    """
    value = _validate_blas_threads(n)
    if value is None:
        yield
        return
    with threadpool_limits(limits=value, user_api="blas"):
        yield


def with_blas_threads(func):
    """Apply a function's ``blas_threads`` argument for the duration of a call."""
    signature = inspect.signature(func)
    parameter = signature.parameters.get("blas_threads")
    if parameter is None:
        raise TypeError("with_blas_threads requires a blas_threads parameter")
    default = parameter.default

    @wraps(func)
    def wrapped(*args, **kwargs):
        bound = signature.bind_partial(*args, **kwargs)
        value = bound.arguments.get("blas_threads", default)
        with blas_thread_limit(value):
            return func(*args, **kwargs)

    return wrapped

