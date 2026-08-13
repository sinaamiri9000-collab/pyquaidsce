"""Parse Stata ``mat list`` output (including the column-block wrapping)."""

from __future__ import annotations

import re
from typing import Dict

import numpy as np

_HDR = re.compile(r"^(\w+)\[(\d+),(\d+)\]\s*$")
_COLS = re.compile(r"^\s*((?:c\d+\s*)+)$")
_ROW = re.compile(r"^\s*r(\d+)\s+(.*\S)\s*$")
_NUM = re.compile(r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?")


def parse_matrices(path: str) -> Dict[str, np.ndarray]:
    """Return every matrix printed by ``mat list`` in *path*."""
    out: Dict[str, np.ndarray] = {}
    lines = open(path, "r", errors="replace").read().split("\n")
    i = 0
    while i < len(lines):
        m = _HDR.match(lines[i].strip())
        if not m:
            i += 1
            continue
        name, nr, nc = m.group(1), int(m.group(2)), int(m.group(3))
        M = np.full((nr, nc), np.nan)
        i += 1
        cols: list[int] = []
        while i < len(lines):
            ln = lines[i]
            if _HDR.match(ln.strip()):
                break
            cm = _COLS.match(ln)
            if cm:
                cols = [int(c[1:]) - 1 for c in cm.group(1).split()]
                i += 1
                continue
            rm = _ROW.match(ln)
            if rm and cols:
                r = int(rm.group(1)) - 1
                vals = [float(v) for v in _NUM.findall(rm.group(2))]
                for c, v in zip(cols, vals):
                    M[r, c] = v
                i += 1
                continue
            if ln.strip().startswith(".") or ln.strip() == "":
                # a blank line inside a matrix block is just the block gap
                if np.isnan(M).all():
                    i += 1
                    continue
                if ln.strip().startswith("."):
                    break
            i += 1
        out[name] = M
    return out


if __name__ == "__main__":
    import sys

    for k, v in parse_matrices(sys.argv[1]).items():
        print(k, v.shape, "nan:", int(np.isnan(v).sum()))
        print(np.array2string(v[:3, :5], precision=6))
