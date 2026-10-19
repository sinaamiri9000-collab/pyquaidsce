"""Generate the simulated household data used by the three example scripts."""

from pathlib import Path

import numpy as np
import pandas as pd


def main():
    rng = np.random.default_rng(123456)
    n = 800
    log_prices = rng.normal(0, 0.18, (n, 4))
    hh_size = rng.integers(1, 6, n)
    urban = rng.integers(0, 2, n)
    income = rng.lognormal(7, 0.5, n)
    spending_error = rng.normal(0, 0.3, n)

    # Income helps explain expenditure. The spending error also affects shares.
    log_expenditure = (
        8.3 + 0.00035 * (income - 1000) + 0.035 * hh_size + 0.06 * urban
        + log_prices @ np.array([0.08, -0.05, 0.04, 0.02]) + spending_error
    )
    shares = (
        np.array([0.30, 0.25, 0.20, 0.25])
        + np.outer(log_expenditure - 8.5, [0.025, -0.015, 0.01, -0.02])
        + np.outer(spending_error, [0.06, -0.04, 0.03, -0.05])
        + rng.normal(0, 0.035, (n, 4))
    )

    # Some households do not buy every good. Keep at least one purchase per row.
    participation = (
        0.85 + 0.12 * (log_expenditure[:, None] - 8.5)
        + 0.08 * urban[:, None] + 0.15 * spending_error[:, None]
        + 0.1 * log_prices
    )
    bought = rng.normal(size=(n, 4)) + participation > 0
    bought[~bought.any(axis=1), 0] = True
    shares = np.where(bought, np.maximum(shares, 0.01), 0)
    shares /= shares.sum(axis=1, keepdims=True)

    data = pd.DataFrame({
        **{f"w{i+1}": shares[:, i] for i in range(4)},
        **{f"p{i+1}": np.exp(log_prices[:, i]) for i in range(4)},
        "expenditure": np.exp(log_expenditure),
        "hh_size": hh_size,
        "urban": urban,
        "income": income,
    })
    path = Path(__file__).with_name("household_demand.csv")
    data.to_csv(path, index=False, float_format="%.12g")
    print(f"Saved {len(data)} simulated households to {path}")


if __name__ == "__main__":
    main()
