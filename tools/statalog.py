"""Parse Stata coefficient tables out of a .log file, for benchmarking."""

from __future__ import annotations

import re
from typing import Dict, List, Tuple

_NUM = r"-?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
_ROW = re.compile(
    rf"^\s*(\S+)\s*\|\s*({_NUM})\s+({_NUM})\s+(-?\d+\.\d+|\.)\s+(\d+\.\d+|\.)\s+"
    rf"({_NUM})\s+({_NUM})\s*$"
)
_EQ = re.compile(r"^([A-Za-z_][A-Za-z_0-9]*)\s*\|\s*$")
_SCALAR = re.compile(r"^(Number of obs|Number of demographics|Alpha_0|"
                     r"Log-likelihood)\s*=\s*(\S+)")


def parse_coef_table(path: str) -> Tuple[Dict[str, Tuple[float, float]],
                                         Dict[str, float], List[str]]:
    """Return {"eq:name": (b, se)}, scalars, and the command lines found."""
    coefs: Dict[str, Tuple[float, float]] = {}
    scalars: Dict[str, float] = {}
    cmds: List[str] = []
    eq = ""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(". ") and ("quaidsce" in line or "quaids" in line):
                cmds.append(line[2:].strip())
            m = _SCALAR.match(line.strip())
            if m:
                v = m.group(2).replace(",", "")
                try:
                    scalars[m.group(1)] = float(v)
                except ValueError:
                    pass
                continue
            m = _EQ.match(line.strip())
            if m:
                eq = m.group(1)
                continue
            m = _ROW.match(line)
            if m:
                nm = m.group(1)
                coefs[f"{eq}:{nm}" if eq else nm] = (float(m.group(2)),
                                                     float(m.group(3)))
    return coefs, scalars, cmds


if __name__ == "__main__":
    import sys

    c, s, cmd = parse_coef_table(sys.argv[1])
    print("commands:", cmd)
    print("scalars :", s)
    print("n coefs :", len(c))
    eqs: Dict[str, int] = {}
    for k in c:
        eqs[k.split(":")[0]] = eqs.get(k.split(":")[0], 0) + 1
    print("by eq   :", eqs)
