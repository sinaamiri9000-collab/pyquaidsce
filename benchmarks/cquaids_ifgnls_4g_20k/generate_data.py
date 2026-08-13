"""Generate the deterministic synthetic dataset used by the 4-good CQUAIDS benchmark.

The dataset is intentionally independent of any confidential or third-party microdata.
It is designed for implementation-equivalence and runtime benchmarking, not for
recovering a particular set of 'true' CQUAIDS structural parameters.
"""
from __future__ import annotations

from pathlib import Path
from datetime import datetime
import hashlib
import json

import numpy as np
import pandas as pd

NOBS = 20_000
SEED = 20260813
HERE = Path(__file__).resolve().parent
DATA_DIR = HERE / "data"
OUT = DATA_DIR / "benchmark_cquaids_4g_20k.dta"
MANIFEST = HERE / "benchmark_manifest.json"


def make_data(n: int = NOBS, seed: int = SEED) -> pd.DataFrame:
    rng = np.random.default_rng(seed)

    # Three standardized demographic shifters.
    z = rng.normal(size=(n, 3))
    z1, z2, z3 = z.T

    # Moderately correlated positive prices.
    cov = 0.04 * np.eye(4) + 0.01 * np.ones((4, 4))
    lnp = rng.multivariate_normal(
        mean=np.array([0.10, -0.05, 0.00, 0.08]), cov=cov, size=n
    )
    prices = np.exp(lnp)

    # Positive system expenditure with economically plausible covariation.
    lnexp = (
        4.60
        + 0.08 * z1
        - 0.05 * z2
        + 0.04 * z3
        + 0.12 * lnp[:, 0]
        - 0.05 * lnp[:, 2]
        + rng.normal(scale=0.25, size=n)
    )
    total = np.exp(lnexp)

    # Smooth latent budget allocation. Softmax guarantees positive latent shares.
    price_load = np.array(
        [
            [0.70, 0.20, -0.10, 0.05],
            [0.10, 0.60, 0.10, -0.10],
            [-0.20, 0.10, 0.70, 0.05],
            [0.05, -0.10, 0.10, 0.55],
        ]
    )
    demo_load = np.array(
        [
            [0.10, -0.05, 0.04],
            [-0.05, 0.08, -0.02],
            [0.03, -0.04, 0.06],
            [-0.08, 0.02, -0.05],
        ]
    )
    score = (
        np.array([-1.00, -1.25, -1.40, -1.10])[None, :]
        + 0.08 * (lnexp[:, None] - lnexp.mean())
        + 0.15 * (lnp @ price_load.T)
        + z @ demo_load.T
        + rng.normal(scale=0.20, size=(n, 4))
    )
    ex = np.exp(score - score.max(axis=1, keepdims=True))
    latent_share = ex / ex.sum(axis=1, keepdims=True)

    # Participation mechanism creates genuine zero shares for all four goods.
    # At least two goods are forced to participate so no observed share equals one.
    p_intercept = np.array([1.40, 1.00, 0.70, 1.20])
    p_price = np.array(
        [
            [0.40, -0.15, 0.10, 0.00],
            [-0.10, 0.35, -0.05, 0.10],
            [0.00, -0.10, 0.30, -0.15],
            [0.20, 0.00, -0.10, 0.25],
        ]
    )
    p_demo = np.array(
        [
            [0.15, -0.10, 0.05],
            [-0.05, 0.12, 0.08],
            [0.10, -0.06, -0.10],
            [-0.12, 0.04, 0.10],
        ]
    )
    participation_index = (
        p_intercept[None, :]
        + lnp @ p_price.T
        + 0.15 * (lnexp[:, None] - lnexp.mean())
        + z @ p_demo.T
        + rng.normal(size=(n, 4))
    )
    active = participation_index > 0
    too_few = np.flatnonzero(active.sum(axis=1) < 2)
    for row in too_few:
        top2 = np.argpartition(participation_index[row], -2)[-2:]
        active[row, top2] = True

    shares = latent_share * active
    shares /= shares.sum(axis=1, keepdims=True)

    df = pd.DataFrame({f"w{i+1}": shares[:, i] for i in range(4)})
    for i in range(4):
        df[f"p{i+1}"] = prices[:, i]
    df["total"] = total
    df["z1"], df["z2"], df["z3"] = z1, z2, z3
    return df


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    df = make_data()
    shares = ["w1", "w2", "w3", "w4"]
    assert len(df) == NOBS
    assert np.isfinite(df.to_numpy()).all()
    assert (df[["p1", "p2", "p3", "p4", "total"]] > 0).all().all()
    assert (df[shares] >= 0).all().all()
    assert np.max(np.abs(df[shares].sum(axis=1).to_numpy() - 1.0)) < 1e-12
    assert not (df[shares] == 1.0).any().any()
    assert ((df[shares] == 0.0).sum() > 0).all()

    # Stata 14+ format; both Stata and pandas read the same stored doubles.
    df.to_stata(
        OUT, write_index=False, version=118,
        time_stamp=datetime(2026, 8, 13, 0, 0, 0),
        data_label="pyquaidsce controlled CQUAIDS benchmark",
    )
    info = {
        "nobs": NOBS,
        "seed": SEED,
        "goods": 4,
        "demographics": 3,
        "share_variables": shares,
        "price_variables": ["p1", "p2", "p3", "p4"],
        "expenditure_variable": "total",
        "demographic_variables": ["z1", "z2", "z3"],
        "zero_counts": {c: int((df[c] == 0.0).sum()) for c in shares},
        "share_sum_max_abs_error": float(
            np.max(np.abs(df[shares].sum(axis=1).to_numpy() - 1.0))
        ),
    }
    info["sha256"] = sha256(OUT)
    MANIFEST.write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
