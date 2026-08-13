"""Controlled same-data runtime benchmark for pyquaidsce 1.0.1."""

from __future__ import annotations

import io
import json
import os
import platform
import sys
import time
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from pyquaidsce import quaidsce


HERE = Path(__file__).resolve().parent
DATA = HERE / "benchmark_11goods_37645.dta"

SHARES = [f"w{i}" for i in range(1, 12)]
PRICES = [
    "tornqvistssb2",
    "tornqvistsweetsnack",
    "tornqvistsweetmeal",
    "tornqvisttea",
    "tornqvistsoursnack",
    "tornqvistfruitveg",
    "tornqvistcereals",
    "tornqvistprotein2",
    "tornqvistdairy",
    "tornqvistoils",
    "tornqvistspices",
]
DEMOGRAPHICS = ["scale", "age", "cfunc"]


def numpy_config_text() -> str:
    buf = io.StringIO()
    with redirect_stdout(buf):
        np.show_config()
    return buf.getvalue()


def main() -> None:
    if not DATA.exists():
        raise FileNotFoundError(
            f"Missing benchmark data: {DATA}. See README.md in this directory."
        )

    # Load before starting the timer so I/O is not part of estimator timing.
    df = pd.read_stata(DATA)
    if len(df) != 37645:
        raise ValueError(f"Expected 37,645 rows; found {len(df):,}.")

    required = SHARES + PRICES + ["tfexp"] + DEMOGRAPHICS
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns: {missing}")

    t0 = time.perf_counter()
    res = quaidsce(
        df,
        shares=SHARES,
        prices=PRICES,
        expenditure="tfexp",
        demographics=DEMOGRAPHICS,
        anot=1.6,
        method="ifgnls",
        algorithm="gn",
        start="zero",
        first_stage_predict="pr",
        strict_stata=True,
        reps=0,
        verbose=False,
    )
    elapsed = time.perf_counter() - t0

    report = {
        "package": "pyquaidsce",
        "package_version": "1.0.1",
        "elapsed_seconds": elapsed,
        "elapsed_minutes": elapsed / 60.0,
        "n": int(len(df)),
        "goods": len(SHARES),
        "demographics": len(DEMOGRAPHICS),
        "method": "ifgnls",
        "algorithm": "gn",
        "start": "zero",
        "anot": 1.6,
        "first_stage_predict": "pr",
        "strict_stata": True,
        "reps": 0,
        "converged": bool(res.converged),
        "n_outer": int(res.n_outer),
        "n_gn": int(res.n_gn),
        "llf": float(res.llf),
        "python": sys.version,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "cpu_count": os.cpu_count(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "pandas": pd.__version__,
        "numpy_config": numpy_config_text(),
    }

    out = HERE / "python_runtime.json"
    out.write_text(json.dumps(report, indent=2), encoding="utf-8")

    print(f"Elapsed: {elapsed:.3f} s ({elapsed / 60:.3f} min)")
    print(f"N: {len(df):,}")
    print(f"Converged: {res.converged}")
    print(f"Outer iterations: {res.n_outer}")
    print(f"Gauss-Newton steps: {res.n_gn}")
    print(f"Log-likelihood: {res.llf:.12f}")
    print(f"Saved: {out}")


if __name__ == "__main__":
    main()
