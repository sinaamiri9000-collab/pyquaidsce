"""Reproducible five-replication PR/XB bootstrap smoke test for release 1.3.0."""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce

DATA = ROOT / "benchmarks/cquaids_ifgnls_4g_20k/data/benchmark_cquaids_4g_20k.dta"
OUTPUT = ROOT / "benchmarks/release_120/results/bootstrap_smoke.json"


def fit_mode(frame: pd.DataFrame, predict: str) -> dict:
    started = time.perf_counter()
    res = quaidsce(
        frame,
        shares=["w1", "w2", "w3", "w4"],
        prices=["p1", "p2", "p3", "p4"],
        expenditure="total",
        demographics=["z1", "z2", "z3"],
        anot=10.0,
        method="ifgnls",
        algorithm="gn",
        start="linear",
        first_stage_predict=predict,
        strict_stata=True,
        reps=5,
        seed=12000 + (0 if predict == "pr" else 1),
        bootstrap_start="zero",
        n_jobs=2,
        mp_context="spawn",
        rep_timeout=120.0,
        verbose=True,
    )
    elapsed = time.perf_counter() - started
    n = res.spec.neqn
    n_elasticities = n + 2 * n * n
    assert res.boot is not None and res.boot.reps_ok >= 2
    assert np.isfinite(res.V).all()
    assert np.allclose(res.V, res.boot.V, rtol=0.0, atol=0.0)
    assert np.allclose(res.se, res.boot.se, rtol=1e-12, atol=1e-14)
    assert np.isnan(res.analytic_se[-n_elasticities:]).all()
    assert np.isfinite(res.se[-n_elasticities:]).all()
    return {
        "predict": predict,
        "elapsed_seconds": elapsed,
        "reps_requested": res.boot.reps_requested,
        "reps_ok": res.boot.reps_ok,
        "failures": res.boot.failures,
        "converged": bool(res.converged),
        "bootstrap_se_min": float(np.min(res.se)),
        "bootstrap_se_max": float(np.max(res.se)),
        "elasticity_se_min": float(np.min(res.se[-n_elasticities:])),
        "elasticity_se_max": float(np.max(res.se[-n_elasticities:])),
        "v_se_synchronized": True,
        "analytic_elasticity_se_are_nan": True,
    }


def verify_timeout(frame: pd.DataFrame) -> str:
    try:
        quaidsce(
            frame,
            shares=["w1", "w2", "w3", "w4"],
            prices=["p1", "p2", "p3", "p4"],
            expenditure="total",
            demographics=["z1", "z2", "z3"],
            anot=10.0,
            method="nls",
            start="linear",
            first_stage_predict="xb",
            reps=2,
            seed=12099,
            n_jobs=1,
            rep_timeout=1e-9,
            verbose=False,
        )
    except RuntimeError as exc:
        message = str(exc)
        assert "TimeoutError" in message
        return message
    raise AssertionError("rep_timeout did not stop deliberately expired replications")


def main() -> None:
    frame = pd.read_stata(DATA).iloc[:3000].reset_index(drop=True)
    report = {
        "package_version": "1.3.0",
        "dataset": str(DATA.relative_to(ROOT)),
        "observations": len(frame),
        "runs": [fit_mode(frame, "pr"), fit_mode(frame, "xb")],
        "timeout_check": verify_timeout(frame),
    }
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
