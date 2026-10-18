"""Explicit label migration for known historical validation fixtures only."""

import re

_ELASTICITY = re.compile(r"^(ELAS_UNCOMP|ELAS_COMP):e_(\d+)_(\d+)$")


def canonical_legacy_name(name: str) -> str:
    """Map an archived price-elasticity label to [good, price] indices."""
    match = _ELASTICITY.fullmatch(name)
    if match is None:
        return name
    block, i, j = match.groups()
    return f"{block}:e_{j}_{i}"
