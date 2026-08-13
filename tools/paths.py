"""Paths shared by the validation utilities."""

from __future__ import annotations

import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCH = PROJECT_ROOT / "bench"


def find_stata_repo() -> Path:
    """Locate the separately distributed ``quaidsce-master`` source tree.

    Set ``QUAIDSCE_STATA_REPO`` to avoid auto-discovery.  The common layout is
    two sibling directories, ``pyquaidsce/`` and ``quaidsce-master/``.
    """
    configured = os.environ.get("QUAIDSCE_STATA_REPO")
    candidates = [] if not configured else [Path(configured).expanduser()]
    candidates += [
        PROJECT_ROOT.parent / "quaidsce-master",
        PROJECT_ROOT / "quaidsce-master",
        PROJECT_ROOT.parent / "src" / "quaidsce-master",
    ]
    for candidate in candidates:
        if (candidate / "quaidsce_c.ado").is_file() and (
            candidate / "data" / "DS_STATA_3_2_0_pci2sls_.dta"
        ).is_file():
            return candidate.resolve()
    tried = "\n  ".join(str(x) for x in candidates)
    raise FileNotFoundError(
        "Could not locate quaidsce-master. Set QUAIDSCE_STATA_REPO to its "
        f"directory. Tried:\n  {tried}"
    )
