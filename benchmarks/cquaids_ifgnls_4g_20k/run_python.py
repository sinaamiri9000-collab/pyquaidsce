"""Run the Python side of the controlled CQUAIDS/IFGNLS benchmark.

The timed interval contains only the estimator call. Data loading and result-file
writing are deliberately outside the timer. By default, common BLAS/OpenMP
thread controls are set to two threads to match the recorded Stata reference
run, which reported two processors available to Stata.
"""
from __future__ import annotations

import os

THREADS = os.environ.get("PYQUAIDSCE_BENCH_THREADS", "2")
for _name in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ[_name] = THREADS

from contextlib import redirect_stdout
from pathlib import Path
import csv
import json
import platform
import sys
import time

import numpy as np
import pandas as pd
import scipy

from pyquaidsce import quaidsce

HERE = Path(__file__).resolve().parent
DATA = HERE / "data" / "benchmark_cquaids_4g_20k.dta"
RESULTS = HERE / "results"
OUTCSV = RESULTS / "python_results.csv"
OUTJSON = RESULTS / "python_runtime.json"
OUTLOG = RESULTS / "python_benchmark.log"

SHARES = ["w1", "w2", "w3", "w4"]
PRICES = ["p1", "p2", "p3", "p4"]
DEMOS = ["z1", "z2", "z3"]


class Tee:
    def __init__(self, *files):
        self.files = files

    def write(self, s):
        for f in self.files:
            f.write(s)
            f.flush()
        return len(s)

    def flush(self):
        for f in self.files:
            f.flush()


def main() -> None:
    RESULTS.mkdir(parents=True, exist_ok=True)
    df = pd.read_stata(DATA)
    zeros = {c: int((df[c] == 0).sum()) for c in SHARES}

    with OUTLOG.open("w", encoding="utf-8", newline="") as lf:
        tee = Tee(sys.stdout, lf)
        with redirect_stdout(tee):
            print("PYQUAIDSCE CONTROLLED CQUAIDS/IFGNLS BENCHMARK")
            print(f"Python: {sys.version.split()[0]}")
            print(f"pyquaidsce: 1.0.1")
            print(f"NumPy: {np.__version__}")
            print(f"SciPy: {scipy.__version__}")
            print(f"pandas: {pd.__version__}")
            print(f"OS/platform: {platform.platform()}")
            print(f"Machine: {platform.machine()}")
            print(f"Processor: {platform.processor()}")
            print(f"Logical CPUs visible to Python: {os.cpu_count()}")
            print(f"Benchmark thread limit: {THREADS}")
            print(f"N: {len(df)}")
            print(f"Zero counts: {zeros}")
            print("Specification: censored QUAIDS, IFGNLS, 4 goods, 3 demographics")
            print("Starting values: zero/default; no user-supplied initial vector")
            print("Bootstrap: disabled")
            print("Stata compatibility: first_stage_predict='pr', strict_stata=True")
            print("\n--- timed estimation starts ---")

            t0 = time.perf_counter()
            res = quaidsce(
                df,
                shares=SHARES,
                prices=PRICES,
                expenditure="total",
                demographics=DEMOS,
                anot=10.0,
                method="ifgnls",
                algorithm="gn",
                start="zero",
                first_stage_predict="pr",
                strict_stata=True,
                reps=0,
                verbose=True,
            )
            runtime = time.perf_counter() - t0

            print("--- timed estimation ends ---")
            print(f"BENCHRUNTIME {runtime:.17g}")
            print(f"BENCHLL {res.llf:.17g}")
            print(f"BENCHN {res.nobs}")
            print(f"CONVERGED {int(bool(res.converged))}")
            print(f"OUTER_ITERATIONS {res.n_outer}")
            print(f"GN_STEPS {res.n_gn}")

    with OUTCSV.open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["name", "python_value", "python_se", "block"])
        for name, value, se in zip(res.names, res.b, res.se):
            w.writerow([name, f"{value:.17g}", f"{se:.17g}", name.split(":", 1)[0]])

    meta = {
        "runtime_seconds": runtime,
        "log_likelihood": float(res.llf),
        "converged": bool(res.converged),
        "outer_iterations": int(res.n_outer),
        "gauss_newton_steps": int(res.n_gn),
        "nobs": int(res.nobs),
        "zero_counts": zeros,
        "python_version": sys.version,
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "pandas_version": pd.__version__,
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor(),
        "logical_cpu_count": os.cpu_count(),
        "benchmark_thread_limit": int(THREADS),
    }
    OUTJSON.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
    print(f"\nSaved benchmark outputs under: {RESULTS}")


if __name__ == "__main__":
    main()
