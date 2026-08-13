"""Result container mirroring the ``e()`` results left behind by ``quaidsce``."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

import numpy as np

from .elasticities import Elasticities, Means
from .params import Coefs, Spec
from .statafmt import coef_table, g


@dataclass
class QuaidsceResults:
    # ---- specification --------------------------------------------------- #
    spec: Spec
    anot: float
    nobs: int
    method: str
    share_names: List[str]
    price_names: List[str]
    demo_names: List[str]
    expenditure_name: str

    # ---- estimates ------------------------------------------------------- #
    coefs: Coefs
    theta: np.ndarray  # free (estimated) parameter vector
    b: np.ndarray  # e(b): full params + tau + elasticities
    V: np.ndarray  # e(V)
    names: List[str]  # e(b) column names ("eq:name")
    b_est: np.ndarray  # e(best): the free parameters
    V_est: np.ndarray  # e(Vest)

    llf: float
    sigma: np.ndarray
    elas: Elasticities
    means: Means

    tau: Optional[np.ndarray] = None  # stacked probit coefficients
    setau: Optional[np.ndarray] = None
    np_prob: int = 0
    probits: List[object] = field(default_factory=list)

    n_outer: int = 0
    n_gn: int = 0
    converged: bool = True
    boot: Optional["BootResult"] = None  # noqa: F821
    notes: List[str] = field(default_factory=list)

    # ------------------------------------------------------------------ #
    @property
    def se(self) -> np.ndarray:
        d = np.diag(self.V).copy()
        d[d < 0] = np.nan
        return np.sqrt(d)

    def named(self) -> Dict[str, float]:
        return dict(zip(self.names, self.b))

    def get(self, key: str) -> float:
        """Look up a coefficient by ``eq:name`` or bare ``name``."""
        d = self.named()
        if key in d:
            return d[key]
        for k, v in d.items():
            if k.split(":", 1)[-1] == key:
                return v
        raise KeyError(key)

    # ------------------------------------------------------------------ #
    def _title(self) -> str:
        if self.spec.censor and self.spec.quadratic:
            return "Censored Quadratic AIDS model"
        if self.spec.quadratic:
            return "Quadratic AIDS model"
        if self.spec.censor:
            return "Censored AIDS model"
        return "AIDS model"

    def summary(self, level: float = 95.0, elasticities: bool = True) -> str:
        title = self._title()
        head = [
            "",
            title,
            "-" * (20 if len(title) > 10 else 10),
            f"Number of obs          = {g(self.nobs, 10):>10}",
            f"Number of demographics = {g(self.spec.ndemo, 10):>10}",
            f"Alpha_0                = {g(self.anot, 10):>10}",
            f"Log-likelihood         = {g(self.llf, 10):>10}",
            "",
        ]
        if self.boot is not None:
            body = coef_table(self.names, self.b, self.boot.se, level=level,
                              bootstrap=True)
        else:
            has_elasticity_block = self.spec.censor
            n_keep = (len(self.names) if elasticities or not has_elasticity_block
                      else len(self.names) - 2 * self.spec.neqn ** 2
                      - self.spec.neqn)
            body = coef_table(
                self.names[:n_keep], self.b[:n_keep], self.se[:n_keep],
                level=level, hide_zero_se=True,
            )
        return "\n".join(head) + "\n" + body

    def elasticity_tables(self) -> str:
        n = self.spec.neqn
        nm = self.share_names
        out = ["", "Expenditure (income) elasticities, at means", "-" * 44]
        w = max(len(x) for x in nm) + 2
        for i in range(n):
            out.append(f"  {nm[i]:<{w}} {self.elas.income[i]:>12.6f}")
        for lab, M in (
            ("Uncompensated (Marshallian) price elasticities [row = good, "
             "column = price]", self.elas.uncompensated),
            ("Compensated (Hicksian) price elasticities [row = good, "
             "column = price]", self.elas.compensated),
        ):
            out += ["", lab, "-" * min(len(lab), 100)]
            out.append(" " * w + "".join(f"{x[:10]:>12}" for x in nm))
            for i in range(n):
                out.append(
                    f"  {nm[i]:<{w}}" + "".join(f"{M[i, j]:>12.6f}"
                                                for j in range(n))
                )
        return "\n".join(out)

    def __str__(self) -> str:  # pragma: no cover
        return self.summary()
