"""Typed bookkeeping for the first-stage participation equations.

The original implementation inferred coefficient positions from the number of
goods and demographics.  That is only valid for the legacy design
``[all prices, log expenditure, all Ray demographics, constant]``.  The
control-function extension permits an independent selection design, so every
consumer of the first-stage coefficients uses this explicit layout instead.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


@dataclass(frozen=True)
class FirstStageLayout:
    """Column positions within one Probit coefficient vector.

    ``ordered_names`` excludes the constant; the Probit routine appends it.
    ``price_positions`` is keyed by the demand-system price column name, so an
    omitted selection price has no entry and therefore a zero selection
    derivative.  Positions are zero based.
    """

    ordered_names: Tuple[str, ...]
    demand_price_names: Tuple[str, ...]
    price_positions: Dict[str, int]
    expenditure_position: Optional[int]
    covariate_positions: Dict[str, int]
    selection_cf_position: Optional[int]
    constant_position: int

    def __post_init__(self) -> None:
        if self.constant_position != len(self.ordered_names):
            raise ValueError("constant_position must follow all Probit regressors")
        if len(set(self.ordered_names)) != len(self.ordered_names):
            raise ValueError("first-stage coefficient names must be unique")
        positions = list(self.price_positions.values())
        positions += list(self.covariate_positions.values())
        if self.expenditure_position is not None:
            positions.append(self.expenditure_position)
        if self.selection_cf_position is not None:
            positions.append(self.selection_cf_position)
        if len(positions) != len(set(positions)):
            raise ValueError("first-stage layout positions must be unique")
        if positions and (min(positions) < 0 or max(positions) >= self.constant_position):
            raise ValueError("first-stage layout position is out of bounds")

    @property
    def width(self) -> int:
        """Number of coefficients per equation, including the constant."""
        return self.constant_position + 1

    def equation_slice(self, equation: int) -> slice:
        return slice(equation * self.width, (equation + 1) * self.width)

    def tau_names(self, neqn: int) -> List[str]:
        out: List[str] = []
        for i in range(1, neqn + 1):
            out.extend(f"tau:{name}_{i}" for name in self.ordered_names)
            out.append(f"tau:cons_{i}")
        return out

    def coefficient(self, tau, equation: int, position: Optional[int]) -> float:
        if position is None:
            return 0.0
        return float(tau[equation * self.width + position])

    def price_position(self, demand_price_index: int) -> Optional[int]:
        name = self.demand_price_names[demand_price_index]
        return self.price_positions.get(name)


def legacy_layout(
    price_names: Sequence[str],
    demographic_names: Sequence[str],
    *,
    include_expenditure: bool,
) -> FirstStageLayout:
    """Construct the historical full-price/full-demographic layout."""
    ordered: List[str] = []
    price_positions: Dict[str, int] = {}
    for j, price in enumerate(price_names):
        price_positions[str(price)] = len(ordered)
        ordered.append(f"p{j + 1}")
    exp_pos = None
    if include_expenditure:
        exp_pos = len(ordered)
        ordered.append("M")
    cov_positions: Dict[str, int] = {}
    for name in demographic_names:
        cov_positions[str(name)] = len(ordered)
        ordered.append(str(name))
    return FirstStageLayout(
        ordered_names=tuple(ordered),
        demand_price_names=tuple(str(x) for x in price_names),
        price_positions=price_positions,
        expenditure_position=exp_pos,
        covariate_positions=cov_positions,
        selection_cf_position=None,
        constant_position=len(ordered),
    )
