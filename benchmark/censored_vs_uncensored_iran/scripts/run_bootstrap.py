"""Run one package-native bootstrap batch for the 10-good U11 benchmark."""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
            "NUMEXPR_NUM_THREADS"):
    os.environ[key] = "1"

HERE = Path(__file__).resolve()
BENCHMARK = HERE.parents[1]
REPO = HERE.parents[3]
sys.path.insert(0, str(REPO / "src"))

import numpy as np
import pandas as pd
import pyquaidsce
from pyquaidsce import quaidsce

GROUPS = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11]
SHARES = [f"w{i}" for i in GROUPS]
PRICES = [f"p{i}" for i in GROUPS]
DEMOGRAPHICS = ["age", "scale", "children"] + [f"region_{i}" for i in range(1, 6)]
INSTRUMENTS = ["log_income_centered", "log_income_centered_sq"]
BATCHES = {
    "first": {"reps": 50, "seed": 20261019, "first_rep": 1},
    "additional": {"reps": 150, "seed": 20261020, "first_rep": 51},
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--model", choices=("censored", "uncensored"), required=True)
    parser.add_argument("--batch", choices=tuple(BATCHES), required=True)
    parser.add_argument("--first-rep", type=int,
                        help="First replication ID for a smaller saved batch")
    parser.add_argument("--reps", type=int,
                        help="Number of replications in this batch")
    parser.add_argument("--seed", type=int,
                        help="Random seed for this batch")
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()
    if args.n_jobs < 1:
        parser.error("--n-jobs must be at least 1")

    defaults = BATCHES[args.batch]
    first_rep = args.first_rep if args.first_rep is not None else defaults["first_rep"]
    reps = args.reps if args.reps is not None else defaults["reps"]
    seed = args.seed if args.seed is not None else defaults["seed"]
    last_rep = first_rep + reps - 1
    if reps < 2:
        parser.error("--reps must be at least 2 for bootstrap inference")
    if args.batch == "first" and (first_rep < 1 or last_rep > 50):
        parser.error("The first batch must use replication IDs 1 through 50")
    if args.batch == "additional" and (first_rep < 51 or last_rep > 200):
        parser.error("Additional batches must use replication IDs 51 through 200")
    censor = args.model == "censored"
    data = pd.read_csv(BENCHMARK / "data" / "iran_food_demand_10goods.csv",
                       float_precision="round_trip")
    if not np.allclose(data[SHARES].sum(axis=1), 1.0, rtol=0, atol=1e-12):
        raise ValueError("The ten retained budget shares must sum to one.")
    log_income = np.log(data["income"].to_numpy(float))
    weight = data["sample_weight"].to_numpy(float)
    centered = log_income - np.average(log_income, weights=weight)
    data[INSTRUMENTS[0]] = centered
    data[INSTRUMENTS[1]] = centered ** 2

    output = BENCHMARK / "results"
    raw = output / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    print(
        f"pyquaidsce={pyquaidsce.__version__}; model={args.model}; "
        f"batch={args.batch}; replication_ids={first_rep}-{last_rep}; "
        f"n={len(data)}; n_jobs={args.n_jobs}; seed={seed}",
        flush=True,
    )
    started = time.perf_counter()
    fit = quaidsce(
        data, shares=SHARES, prices=PRICES, expenditure="expenditure",
        demographics=DEMOGRAPHICS, ivexp=INSTRUMENTS, anot=11,
        quadratic=True, censor=censor, method="ifgnls", start="zero",
        reps=reps, seed=seed, n_jobs=args.n_jobs,
        blas_threads=1, bootstrap_start="warm", max_outer=400,
        max_iter=400, chunk=15000, verbose=True,
    )
    elapsed = time.perf_counter() - started
    if not fit.converged:
        raise RuntimeError("The full-sample model did not converge.")
    if fit.boot is None or fit.boot.reps_ok < 2:
        raise RuntimeError("Fewer than two bootstrap replications succeeded.")

    model = args.model
    if (first_rep, reps, seed) == (defaults["first_rep"], defaults["reps"], defaults["seed"]):
        prefix = f"{model}_batch_{args.batch}"
    else:
        prefix = f"{model}_batch_{args.batch}_{first_rep:03d}_{last_rep:03d}"
    draws = pd.DataFrame(np.asarray(fit.boot.b_star, float), columns=fit.names)
    # BootResult keeps failure messages but returns only successful draws.
    # Its successful rows are sorted by their original replication number.
    failed_local = {
        int(message.split(":", 1)[0].split()[1])
        for message in fit.boot.failures
    }
    successful_local = [i for i in range(1, reps + 1)
                        if i not in failed_local]
    if len(successful_local) != len(draws):
        raise ValueError("Could not map successful draws to their replication IDs.")
    draws.insert(0, "replication", [first_rep + i - 1
                                     for i in successful_local])
    draws.to_csv(raw / f"{prefix}_draws.csv", index=False, float_format="%.17g")

    summary = {
        "package_version": pyquaidsce.__version__,
        "model": model,
        "censor": censor,
        "batch": args.batch,
        "replication_ids": [first_rep, last_rep],
        "first_replication_id": first_rep,
        "requested": int(fit.boot.reps_requested),
        "successful": int(fit.boot.reps_ok),
        "failures": list(fit.boot.failures),
        "seed": int(seed),
        "workers": int(args.n_jobs),
        "resampling": "independent household rows with replacement",
        "reduced_form": "package ivexp ordinary least squares, rebuilt in each draw",
        "excluded_instruments": INSTRUMENTS,
        "observations": int(fit.nobs),
        "original_group_numbers": GROUPS,
        "anot": 11,
        "elapsed_seconds_including_full_sample_fit": elapsed,
        "full_sample_converged": bool(fit.converged),
        "full_sample_outer_iterations": int(fit.n_outer),
        "full_sample_gn_iterations": int(fit.n_gn),
    }
    (raw / f"{prefix}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    # Retain the full-sample point estimates and elasticity matrices. They are
    # deterministic across saved batches; check consistency before saving.
    point = pd.DataFrame({"parameter": fit.names, "estimate": fit.b})
    point_path = output / f"{model}_point_estimates.csv"
    if point_path.exists():
        old = pd.read_csv(point_path)
        if old.parameter.tolist() != point.parameter.tolist() or not np.allclose(
                old.estimate, point.estimate, rtol=1e-10, atol=1e-11):
            raise ValueError("Full-sample point estimates changed between batches.")
    else:
        point.to_csv(point_path, index=False, float_format="%.17g")

    labels = pd.read_csv(BENCHMARK / "data" / "group_labels.csv")
    label_map = dict(zip(labels.original_group, labels.label))
    point_elasticities = pd.DataFrame({
        "original_group": GROUPS,
        "label": [label_map[g] for g in GROUPS],
        "expenditure_elasticity": fit.elas.income,
        "own_price_marshallian": np.diag(fit.elas.uncompensated),
        "own_price_hicksian": np.diag(fit.elas.compensated),
    })
    point_elasticities.to_csv(output / f"{model}_point_elasticities.csv",
                              index=False, float_format="%.17g")
    for kind, values in (("marshallian", fit.elas.uncompensated),
                         ("hicksian", fit.elas.compensated)):
        matrix = pd.DataFrame(values,
            index=[f"{g}: {label_map[g]}" for g in GROUPS],
            columns=[f"{g}: {label_map[g]}" for g in GROUPS])
        matrix.to_csv(output / f"{model}_{kind}_elasticities.csv",
                      index_label="demand_group", float_format="%.17g")
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
