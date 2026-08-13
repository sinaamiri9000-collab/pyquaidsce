"""Reproduce Stata's ``%w.0g`` display format and its ``_coef_table`` layout.

The point of getting this exactly right is that the Python output can then be
diffed line by line against a Stata ``.log`` file.
"""

from __future__ import annotations

import math
from typing import List, Optional, Sequence

import numpy as np
from scipy.stats import norm


def _fixed(x: float, sig: int) -> str:
    """Fixed-notation string for *x* with *sig* significant digits, trailing
    zeros removed (Stata's ``%g`` never pads)."""
    if x == 0.0:
        return "0"
    mag = math.floor(math.log10(abs(x)))
    dec = max(0, sig - 1 - mag)
    s = f"{x:.{dec}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    if s in ("", "-"):
        s = "0"
    return s


def _strip0(s: str) -> str:
    if s.startswith("0."):
        return s[1:]
    if s.startswith("-0."):
        return "-" + s[2:]
    return s


def g(x, width: int = 9) -> str:
    """Stata's ``%<width>.0g``.

    Stata shows at most ``width - 2`` significant digits, in fixed notation with
    the leading zero dropped, reducing precision until the result fits in
    *width* characters and falling back to exponential notation if it cannot.
    """
    if x is None:
        return "."
    x = float(x)
    if not np.isfinite(x):
        return "."
    if x == 0.0:
        return "0"

    def nsig(s: str) -> int:
        digits = s.split("e")[0].lstrip("-").replace(".", "").lstrip("0")
        return len(digits.rstrip("0")) or 1

    best_fixed = None
    for sig in range(width - 2, 0, -1):
        s = _strip0(_fixed(x, sig))
        if len(s) <= width and float(s) != 0.0:
            best_fixed = s
            break

    best_exp = None
    for dd in range(width - 6, -1, -1):
        s = f"{x:.{dd}e}"
        mant, _, ex = s.partition("e")
        if "." in mant:
            mant = mant.rstrip("0").rstrip(".")
        s = f"{mant}e{ex}"
        if len(s) <= width:
            best_exp = s
            break

    if best_fixed is None:
        return best_exp if best_exp is not None else f"{x:.0e}"
    if best_exp is None:
        return best_fixed
    # Stata keeps whichever representation carries more significant digits,
    # preferring fixed notation on a tie.
    return best_exp if nsig(best_exp) > nsig(best_fixed) else best_fixed


# --------------------------------------------------------------------------- #
#  _coef_table
# --------------------------------------------------------------------------- #
_RULE = "-" * 78
_SEP = "-" * 13 + "+" + "-" * 64


def _header(level: float, bootstrap: bool) -> List[str]:
    lv = f"{level:g}"
    if bootstrap:
        return [
            "             |   Observed   Bootstrap"
            "                         Normal-based",
            "             | coefficient  std. err.      z    P>|z|"
            f"     [{lv}% conf. interval]",
        ]
    return [
        "             | Coefficient  Std. err.      z    P>|z|"
        f"     [{lv}% conf. interval]"
    ]


def coef_table(
    names: Sequence[str],
    b: Sequence[float],
    se: Sequence[float],
    level: float = 95.0,
    bootstrap: bool = False,
    hide_zero_se: bool = False,
    extra: Optional[str] = None,
) -> str:
    """Render a Stata coefficient table, 78 columns wide.

    ``names`` follow the ``equation:name`` convention; an equation header row is
    emitted whenever the equation changes, as Stata does.
    """
    b = np.asarray(b, dtype=float)
    se = np.asarray(se, dtype=float)
    zc = norm.ppf(0.5 + level / 200.0)

    lines: List[str] = []
    if extra:
        lines.append(extra)
    lines.append(_RULE)
    lines += _header(level, bootstrap)
    lines.append(_SEP)

    cur = None
    for nm, bi, si in zip(names, b, se):
        eq, _, sub = nm.partition(":")
        if not sub:
            eq, sub = "", nm
        if eq != cur:
            if cur is not None:
                lines.append(_SEP)
            if eq:
                lines.append(f"{eq:<13}|")
            cur = eq
        if hide_zero_se and (not np.isfinite(si) or si == 0.0):
            lines.append(
                f"{sub:>12} |{g(bi, 9):>11}{'.':>11}{'.':>9}{'.':>8}"
                f"{'.':>13}{'.':>12}"
            )
            continue
        z = bi / si if si > 0 else np.nan
        p = 2.0 * norm.sf(abs(z)) if np.isfinite(z) else np.nan
        lo, hi = bi - zc * si, bi + zc * si
        lines.append(
            f"{sub:>12} |{g(bi, 9):>11}{g(si, 9):>11}"
            f"{z:>9.2f}{p:>8.3f}{g(lo, 9):>13}{g(hi, 9):>12}"
        )
    lines.append(_RULE)
    return "\n".join(lines)
