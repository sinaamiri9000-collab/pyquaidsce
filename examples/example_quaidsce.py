"""Run these examples with: python examples/example_quaidsce.py"""

from pathlib import Path

import pandas as pd
from pyquaidsce import quaidsce


def main():
    # Make sure pyquaidsce is installed in your Python environment.
    # From a terminal: python -m pip install pyquaidsce

    # Read the practice dataset included with these examples.
    data_file = Path(__file__).with_name("data") / "household_demand.csv"
    df = pd.read_csv(data_file)
    shares = ["w1", "w2", "w3", "w4"]
    prices = ["p1", "p2", "p3", "p4"]
    demographics = ["hh_size", "urban"]

    # 1. Basic run: no bootstrap and no expenditure endogeneity correction.
    basic = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, reps=0,
    )
    print(basic.summary())
    print(basic.elasticity_tables())

    # 2. Use income as an instrument to account for endogenous expenditure.
    # The package estimates the expenditure equation and uses its residual.
    iv = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, reps=0, ivexp=["income"],
    )
    print(iv.reduced_form_table())

    # 3. Run bootstrap with the expenditure correction.
    # reps=0 skips bootstrap. reps=10 attempts 10 replications.
    # reps=1 attempts one, but at least two must succeed for standard errors.
    # We use 10 here for practice; increase reps for your research results.
    boot = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, ivexp=["income"],
        reps=10, seed=12345,
    )
    print(boot.elasticity_tables())
    print(f"Successful bootstrap replications: {boot.boot.reps_ok}/10")

    # 4. Run bootstrap on two workers at the same time.
    # warm starts each replication from the full-sample estimates.
    # Each replication may run for at most 900 seconds. A longer replication
    # counts as unsuccessful, and the program moves on to the next one.
    parallel = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, ivexp=["income"],
        reps=10, n_jobs=2, seed=12345, bootstrap_start="warm",
        rep_timeout=900,
    )
    print(parallel.elasticity_tables())

    # 5. Reuse starting values from the same model fitted in example 2.
    # Use theta, not b: theta contains only the free demand parameters.
    restarted = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, ivexp=["income"], reps=0,
        initial=iv.theta, sigma_initial=iv.sigma,
    )
    print(restarted.elasticity_tables())

    # 6. Let the package calculate starting values from a simpler AIDS fit.
    linear_start = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, start="linear", reps=0,
    )
    print(linear_start.elasticity_tables())

    # 7. Fit AIDS by removing the quadratic expenditure term.
    aids = quaidsce(
        df, shares=shares, prices=prices, expenditure="expenditure",
        demographics=demographics, anot=8.0, quadratic=False, reps=0,
    )
    print(aids.elasticity_tables())


# This guard prevents bootstrap workers from running the examples again.
if __name__ == "__main__":
    main()
