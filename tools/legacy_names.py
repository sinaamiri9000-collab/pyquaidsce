"""Explicit label migration for known historical validation fixtures only."""

import re

_ELASTICITY = re.compile(r"^(ELAS_UNCOMP|ELAS_COMP):e_(\d+)_(\d+)$")


def canonical_legacy_name(name: str) -> str:
    """Swap price-elasticity indices in a verified legacy label.

    This is not automatic version detection. Never apply it to already-correct
    output. Values and covariance positions are deliberately left unchanged.
    """
    match = _ELASTICITY.fullmatch(name)
    if match is None:
        return name
    block, i, j = match.groups()
    return f"{block}:e_{j}_{i}"
