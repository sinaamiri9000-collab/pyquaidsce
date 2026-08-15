"""Verify deterministic parallel bootstrap ordering for release 1.2.0."""

from __future__ import annotations

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce


DATA = ROOT / "benchmarks/cquaids_ifgnls_4g_20k/data/benchmark_cquaids_4g_20k.dta"
OUTPUT = ROOT / "benchmarks/release_120/results/bootstrap_determinism.json"


def fit(frame):
    return quaidsce(
        frame,
        shares=["w1", "w2", "w3", "w4"],
        prices=["p1", "p2", "p3", "p4"],
        expenditure="total",
        demographics=["z1", "z2", "z3"],
        anot=10.0,
        method="fgnls",
        start="linear",
        first_stage_predict="xb",
        reps=3,
        seed=12077,
        n_jobs=2,
        mp_context="spawn",
        verbose=False,
    )


def main() -> None:
    frame = pd.read_stata(DATA).iloc[:1000].reset_index(drop=True)
    first = fit(frame)
    second = fit(frame)
    assert first.boot is not None and second.boot is not None
    exact = np.array_equal(first.boot.b_star, second.boot.b_star)
    assert exact
    assert first.boot.failures == second.boot.failures
    report = {
        "observations": len(frame),
        "reps": 3,
        "n_jobs": 2,
        "seed": 12077,
        "replications_first": first.boot.reps_ok,
        "replications_second": second.boot.reps_ok,
        "failures": first.boot.failures,
        "failure_records_equal": True,
        "b_star_exactly_equal": exact,
        "max_abs_difference": float(np.max(np.abs(
            first.boot.b_star - second.boot.b_star
        ))),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
