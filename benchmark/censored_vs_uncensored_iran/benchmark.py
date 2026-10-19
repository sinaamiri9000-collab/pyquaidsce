"""Compare QUAIDS elasticities with and without censoring correction."""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import os
from pathlib import Path
import sys
import zipfile

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))

import numpy as np
import pandas as pd
import pyquaidsce
from pyquaidsce import quaidsce

WORKERS = 5
GROUPS = [
    "Sugar-sweetened beverages", "Sweet snacks and added-sugar foods",
    "Non-sugar beverages", "Nuts, dried fruits, dates and honey",
    "Fruits and vegetables", "Bread and other cereal products", "Protein foods",
    "Dairy and eggs", "Oils, fats, spices and condiments", "Rice",
]
SHARES = [f"w{i}" for i in range(1, 11)]
PRICES = [f"p{i}" for i in range(1, 11)]
DEMOGRAPHICS = ["age", "scale", "children"] + [f"region_{i}" for i in range(1, 6)]
INSTRUMENTS = ["log_income_centered", "log_income_centered_sq"]
BATCHES = [(50, 20261019)] + [(25, seed) for seed in range(20261020, 20261026)]
MODELS = ("censored", "uncensored")
DATA = HERE / "iran_food_demand.csv"
ARCHIVE = HERE / "bootstrap_results.zip"
LOG = HERE / "benchmark.log"


def load_data():
    data = pd.read_csv(DATA, float_precision="round_trip")
    columns = SHARES + PRICES + ["expenditure", "income", "sample_weight"] + DEMOGRAPHICS
    if len(data) != 37656 or not np.isfinite(data[columns]).all().all():
        raise ValueError("Expected 37,656 complete household observations.")
    if (data[SHARES] < 0).any().any() or not np.allclose(
        data[SHARES].sum(axis=1), 1, rtol=0, atol=1e-12
    ):
        raise ValueError("Budget shares must be nonnegative and sum to one.")
    if (data[PRICES + ["expenditure", "income", "sample_weight"]] <= 0).any().any():
        raise ValueError("Prices, expenditure, income and survey weights must be positive.")
    log_income = np.log(data.income.to_numpy())
    centered = log_income - np.average(log_income, weights=data.sample_weight)
    data[INSTRUMENTS[0]] = centered
    data[INSTRUMENTS[1]] = centered ** 2
    return data


def estimate(data):
    results = {}
    model_results = {}
    for model in MODELS:
        batches = []
        point = None
        for reps, seed in BATCHES:
            print(f"Estimating {model}: {sum(len(b) for b in batches)}/{sum(r for r, _ in BATCHES)} bootstrap draws completed", flush=True)
            fit = quaidsce(
                data, shares=SHARES, prices=PRICES, expenditure="expenditure",
                demographics=DEMOGRAPHICS, ivexp=INSTRUMENTS, anot=11,
                quadratic=True, censor=model == "censored", method="ifgnls",
                start="zero", reps=reps, seed=seed, n_jobs=WORKERS,
                blas_threads=1, bootstrap_start="warm", max_outer=400,
                max_iter=400, chunk=15000, verbose=True,
            )
            if not fit.converged or fit.boot is None or fit.boot.reps_ok != reps:
                raise RuntimeError(f"{model}: all requested fits must converge for a paired comparison.")
            current = pd.DataFrame({"parameter": fit.names, "estimate": fit.b})
            if point is not None:
                if point.parameter.tolist() != current.parameter.tolist():
                    raise ValueError("Parameter names changed between batches.")
                np.testing.assert_allclose(point.estimate, current.estimate, rtol=1e-10, atol=1e-11)
            point = current
            batches.append(pd.DataFrame(fit.boot.b_star, columns=fit.names))
        draws = pd.concat(batches, ignore_index=True)
        draws.insert(0, "replication", np.arange(1, len(draws) + 1))
        results[model] = (point, draws)
        model_results[model] = {
            "successful_replications": len(draws), "failed_replications": 0,
            "converged": True,
        }
    settings = {
        "observations": len(data), "goods": 10,
        "bootstrap_replications": sum(reps for reps, _ in BATCHES),
        "batches": [{"replications": reps, "seed": seed} for reps, seed in BATCHES],
        "data_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "method": "ifgnls", "anot": 11, "start": "zero", "bootstrap_start": "warm",
        "max_outer": 400, "max_iter": 400, "chunk": 15000, "blas_threads": 1,
        "resampling": "Independent household rows with replacement; matched across models.",
        "instruments": INSTRUMENTS, "model_results": model_results,
        "software": {"pyquaidsce": pyquaidsce.__version__, "numpy": np.__version__, "pandas": pd.__version__},
        "workers": WORKERS,
    }
    return results, settings


def read_results():
    with zipfile.ZipFile(ARCHIVE) as archive:
        settings = json.loads(archive.read("settings.json"))
        if settings["data_sha256"] != hashlib.sha256(DATA.read_bytes()).hexdigest():
            raise ValueError("The saved results refer to a different dataset.")
        results = {}
        for model in MODELS:
            point = pd.read_csv(io.BytesIO(archive.read(f"{model}_point_estimates.csv")), float_precision="round_trip")
            draws = pd.read_csv(io.BytesIO(archive.read(f"{model}_bootstrap_draws.csv")), float_precision="round_trip")
            results[model] = (point, draws)
    return results, settings


def report(results, settings):
    for model, (point, draws) in results.items():
        if draws.replication.tolist() != list(range(1, settings["bootstrap_replications"] + 1)):
            raise ValueError(f"{model}: missing or repeated bootstrap draws.")
        if point.parameter.tolist() != draws.columns[1:].tolist():
            raise ValueError(f"{model}: parameter names do not match the saved draws.")
        if not np.isfinite(draws).all().all() or not np.isfinite(point.estimate).all():
            raise ValueError(f"{model}: non-finite estimates.")
    lines = [
        "Censored vs uncensored QUAIDS", "",
        f"Households: {settings['observations']:,}; food groups: 10.",
        "Estimation: IFGNLS; anot=11; demographic controls; expenditure control function (ivexp).",
        "Both models use the same households and inputs.",
        "Instruments: centered log household income and its square.",
        "Survey weights center log income; both expenditure and demand equations are unweighted.",
        f"Bootstrap: {settings['bootstrap_replications']} successful household draws per model; no failed draws.",
        "Matching draw IDs use the same household resamples for both models.",
        "The expenditure regression is re-estimated in every draw.",
        "Prices are fixed; intervals exclude uncertainty from constructing prices.",
        "Intervals: 95% percentile bootstrap intervals.", "",
    ]
    for title, parameter in [
        ("Expenditure elasticities", lambda i: f"ELAS_INC:e_{i}"),
        ("Uncompensated own-price elasticities", lambda i: f"ELAS_UNCOMP:e_{i}_{i}"),
    ]:
        lines.extend([title, "", f"{'Food group':<39} {'With correction [95% CI]':>27} {'Without correction [95% CI]':>30}", "-" * 100])
        for i, label in enumerate(GROUPS, 1):
            cells = []
            for model in MODELS:
                point, draws = results[model]
                value = point.set_index("parameter").estimate[parameter(i)]
                lo, hi = np.percentile(draws[parameter(i)], [2.5, 97.5])
                cells.append(f"{value:.3f} [{lo:.3f}, {hi:.3f}]")
            lines.append(f"{label:<39} {cells[0]:>27} {cells[1]:>30}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def save_results(results, settings):
    temporary = ARCHIVE.with_suffix(".zip.tmp")
    with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for model, (point, draws) in results.items():
            for suffix, table in [("point_estimates", point), ("bootstrap_draws", draws)]:
                archive.writestr(f"{model}_{suffix}.csv", table.to_csv(index=False, float_format="%.17g"))
        archive.writestr("settings.json", json.dumps(settings, indent=2) + "\n")
    temporary.replace(ARCHIVE)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summarize", action="store_true", help="Rebuild the log from saved results without estimating the models.")
    args = parser.parse_args()
    if args.summarize:
        results, settings = read_results()
    else:
        results, settings = estimate(load_data())
    text = report(results, settings)
    if not args.summarize:
        save_results(results, settings)
    LOG.write_text(text, encoding="utf-8")
    print(text)
    print(f"Results: {ARCHIVE.name}; log: {LOG.name}")


if __name__ == "__main__":
    main()
