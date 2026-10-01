"""Model interface consumed by the shared nonlinear SUR optimizer.

Model equations and parameter restrictions stay outside the optimizer. The
default adapter preserves the public QUAIDS solver API and its fast Jacobian.
"""

from __future__ import annotations

from typing import Any, Protocol

import numpy as np

from .data import DemandData


class ModelBackend(Protocol):
    start_methods: tuple

    def make_cache(self, spec: Any) -> Any: ...

    def fitted(self, theta: np.ndarray, data: DemandData, spec: Any) -> np.ndarray: ...

    def jacobian(self, theta: np.ndarray, data: DemandData, spec: Any,
                 cache: Any, rows: slice) -> np.ndarray: ...

    def initial(self, data: DemandData, spec: Any, start: str) -> np.ndarray: ...


class QuaidsBackend:
    start_methods = ("zero", "linear")

    def make_cache(self, spec):
        from .jacfree import make_cache
        return make_cache(spec)

    def fitted(self, theta, data, spec):
        from .model import fitted_shares
        return fitted_shares(theta, data, spec)

    def jacobian(self, theta, data, spec, cache, rows):
        from .jacfree import jacobian_free
        return jacobian_free(theta, data, spec, cache, rows)

    def initial(self, data, spec, start):
        if start == "linear":
            from .start import linear_start
            return linear_start(data, spec)
        return np.zeros(spec.n_free)


QUAIDS_BACKEND = QuaidsBackend()
