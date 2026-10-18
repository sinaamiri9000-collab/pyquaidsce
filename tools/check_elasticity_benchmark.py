"""Original fit plus exactly one bootstrap draw for elasticity-label checks.

Example:
    python tools/check_elasticity_benchmark.py --data food_demand_benchmark_4goods.dta \
        --stata-log benchmark_speed_200reps.log --output-dir results

The reference is the observed-coefficient table of the user's existing
quaidsce_c run, not a newly executed Stata run. One draw cannot estimate a
bootstrap covariance; this diagnostic deliberately leaves inference undefined.
"""

from __future__ import annotations

import argparse
import contextlib
from datetime import datetime, timezone
import hashlib
import platform
from pathlib import Path
import re
import subprocess
import sys
import time
import traceback

import numpy as np
import pandas as pd
import scipy

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
import pyquaidsce
from pyquaidsce import quaidsce


def read_stata_observed(path: Path) -> dict:
    """Read the final observed table, whose labels are already canonical."""
    tables, rows, equation = [], {}, None
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^([A-Za-z_]+)\s+\|\s*$", line)
        if match:
            if match[1] == "alpha":
                if rows:
                    tables.append(rows)
                rows = {}
            equation = match[1]
        match = re.match(r"^\s+(\S+)\s+\|\s+([-.\deE+]+)\s+([-.\deE+]+)", line)
        if match and equation:
            rows[f"{equation}:{match[1]}"] = float(match[2])
    if rows:
        tables.append(rows)
    if len(tables) != 2 or len(tables[-1]) != 104:
        raise ValueError("expected Python and quaidsce_c tables with 104 observed values")
    return tables[-1]


def check_labels(result) -> int:
    checks = 0
    n = result.spec.neqn
    for i in range(n):
        assert result.get(f"ELAS_INC:e_{i+1}") == result.elas.income[i]
        checks += 1
        for j in range(n):
            for block, matrix in (("ELAS_UNCOMP", result.elas.uncompensated),
                                  ("ELAS_COMP", result.elas.compensated)):
                assert result.get(f"{block}:e_{i+1}_{j+1}") == matrix[i, j]
                checks += 1
    return checks


class Tee:
    def __init__(self, *streams):
        self.streams = streams

    def write(self, value):
        for stream in self.streams:
            stream.write(value)
        return len(value)

    def flush(self):
        for stream in self.streams:
            stream.flush()


def run(args):
    revision = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
    ).strip()
    dirty = bool(subprocess.check_output(
        ["git", "status", "--porcelain"], cwd=ROOT, text=True
    ).strip())
    print("PYQUAIDSCE PUBLICATION: ORIGINAL FIT + ONE BOOTSTRAP REPLICATION")
    print("UTC:", datetime.now(timezone.utc).isoformat())
    print("Git revision:", revision, "working_tree_dirty:", dirty)
    print("Package:", pyquaidsce.__version__, "module:", pyquaidsce.__file__)
    print("Environment:", platform.platform(), "Python", platform.python_version())
    print("NumPy:", np.__version__, "SciPy:", scipy.__version__, "pandas:", pd.__version__)
    for path in (args.data, args.stata_log):
        print("Input:", path.name, "sha256:", hashlib.sha256(path.read_bytes()).hexdigest())
    print("Bootstrap requested/completed: exactly 1 resampled fit, separately from the original fit.")
    print("Bootstrap SE, covariance, z tests and confidence intervals: NOT ESTIMABLE WITH ONE DRAW.")
    print("Standard summary tables below use conditional analytical SE where implemented; elasticity inference is shown as '.'.")
    print("No Stata process is run. Reference values are rounded observed coefficients in the supplied log.")
    frame = pd.read_stata(args.data)
    kwargs = dict(shares=["w1", "w2", "w3", "w4"], prices=["p1", "p2", "p3", "p4"],
                  expenditure="expenditure", demographics=["age", "scale"], anot=11.0,
                  method="ifgnls", reps=0, start="zero", bootstrap_start="zero",
                  blas_threads=1, verbose=True)
    print("Specification:", kwargs)
    print("Observations:", len(frame), "missing cells:", int(frame.isna().sum().sum()))
    print("\n=== ORIGINAL FIT ===")
    started = time.perf_counter()
    result = quaidsce(frame, **kwargs)
    original_seconds = time.perf_counter() - started
    if not result.converged:
        raise RuntimeError("original fit did not converge")
    checks = check_labels(result)
    print("Converged:", result.converged, "outer:", result.n_outer, "inner:", result.n_gn)
    print("Log-likelihood:", format(result.llf, ".15g"), "seconds:", original_seconds)
    print("Name/matrix assertions passed:", checks)
    if args.baseline_npz is not None:
        with np.load(args.baseline_npz) as before:
            current = dict(b=result.b, theta=result.theta, V=result.V,
                           inc=result.elas.income, u=result.elas.uncompensated,
                           c=result.elas.compensated, sigma=result.sigma, llf=result.llf)
            for name, value in current.items():
                np.testing.assert_array_equal(value, before[name])
        print("Before/after numeric invariance: PASS (b, theta, V, all elasticities, Sigma, llf exactly unchanged).")
    print("\n=== STANDARD PACKAGE COEFFICIENT TABLE: ORIGINAL FIT ===")
    print(result.summary())
    print(result.elasticity_tables())
    print("\n=== EXACTLY ONE BOOTSTRAP REPLICATION ===")
    # Identical seed derivation and row sampling to bootstrap.py's _one_rep.
    seed = int(np.random.SeedSequence(args.seed).generate_state(1, dtype=np.uint32)[0])
    indices = np.random.default_rng(seed).integers(0, len(frame), size=len(frame))
    resampled = frame.iloc[indices].reset_index(drop=True)
    print("Master seed:", args.seed, "replication seed:", seed)
    print("Resample index sha256:", hashlib.sha256(indices.astype('<i8').tobytes()).hexdigest())
    started = time.perf_counter()
    replica = quaidsce(resampled, **kwargs)
    replica_seconds = time.perf_counter() - started
    if not replica.converged:
        raise RuntimeError("single bootstrap replication did not converge")
    replica_checks = check_labels(replica)
    assert replica.names == result.names
    print("Converged:", replica.converged, "outer:", replica.n_outer, "inner:", replica.n_gn)
    print("Log-likelihood:", format(replica.llf, ".15g"), "seconds:", replica_seconds)
    print("Name/matrix assertions passed:", replica_checks)
    print("\n=== STANDARD PACKAGE COEFFICIENT TABLE: SINGLE RESAMPLED FIT ===")
    print(replica.summary())
    print(replica.elasticity_tables())
    reference = read_stata_observed(args.stata_log)
    assert set(result.names) == set(reference)
    rows = []
    print("\n=== ALL REPORTED VALUES (NO BOOTSTRAP STANDARD ERRORS) ===")
    print(f"{'Name':<28} {'Original Python':>19} {'Stata observed':>19} {'Abs difference':>17} {'One bootstrap draw':>20}")
    for name, value, draw in zip(result.names, result.b, replica.b):
        stata = reference[name]
        difference = abs(float(value) - stata)
        rows.append((name, float(value), stata, difference, float(draw)))
        print(f"{name:<28} {value:19.12g} {stata:19.12g} {difference:17.10g} {draw:20.12g}")
    maximum = max(row[3] for row in rows)
    elasticity_maximum = max(row[3] for row in rows if row[0].startswith("ELAS_"))
    print("\nMaximum absolute difference, all observed values:", maximum)
    print("Maximum absolute difference, elasticities:", elasticity_maximum)
    print("All original and resampled name/matrix checks passed:", checks + replica_checks)
    lines = ["# بررسی کشش‌های شاخه publication با یک تکرار", "",
             f"نسخه کد: `{revision}`؛ تغییر ثبت‌نشده در زمان اجرا: `{dirty}`.", "",
             "برآورد اصلی روی ۲٬۶۶۵ مشاهده و دقیقاً یک تکرار بوت‌استرپ اجرا شده است. "
             "ستون Stata از برآوردهای نقطه‌ایِ گزارش ارسالی قبلی خوانده شده؛ Stata دوباره اجرا نشده است.", "",
             "**با یک تکرار، خطای استاندارد و فاصله اطمینان بوت‌استرپ قابل محاسبه نیست.** "
             "ستون تکرار، برآورد یک نمونه بازنمونه‌گیری‌شده است و نباید با برآورد اصلی یکسان باشد.", "",
             f"هر دو برآورد همگرا شدند؛ هر {checks + replica_checks} تطبیق نام و عنصر ماتریس موفق بود.",
             f"بیشترین اختلاف با ارقام گرد‌شده Stata: `{maximum:.10g}` برای تمام خروجی‌ها "
             f"و `{elasticity_maximum:.10g}` برای کشش‌ها.", "",
             f"بذر اصلی: `{args.seed}`؛ بذر تکرار: `{seed}`. "
             "ترتیب اندیس‌ها: کالای مورد تقاضا، سپس قیمت تغییریافته.", "",
             "زمان‌های این محیط فقط برای ثبت اجرا هستند و با زمان‌های دستگاه بنچمارک قبلی قابل مقایسه نیستند.", ""]
    def add_table(title, selected):
        lines.extend([f"## {title}", "",
                      "| نام | Python؛ برآورد اصلی | quaidsce؛ برآورد اصلی قبلی | اختلاف مطلق | یک تکرار بوت‌استرپ Python |",
                      "|---|---:|---:|---:|---:|"])
        for name, value, stata, difference, draw in selected:
            lines.append(f"| `{name}` | {value:.12g} | {stata:.12g} | {difference:.10g} | {draw:.12g} |")
        lines.append("")
    add_table("تمام کشش‌ها", [r for r in rows if r[0].startswith("ELAS_")])
    add_table("ضرایب ساختاری و مرحله اول", [r for r in rows if not r[0].startswith("ELAS_")])
    (args.output_dir / "publication_benchmark_comparison.md").write_text("\n".join(lines), encoding="utf-8")
    print("\nPASS: elasticity labels are consistent in the original fit and the single resampled fit.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, required=True)
    parser.add_argument("--stata-log", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=12345)
    parser.add_argument("--baseline-npz", type=Path,
                        help="optional original-fit arrays captured before the naming fix")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    with (args.output_dir / "publication_benchmark_1rep.log.txt").open("w", encoding="utf-8") as log:
        with contextlib.redirect_stdout(Tee(sys.stdout, log)), contextlib.redirect_stderr(Tee(sys.stderr, log)):
            try:
                run(args)
            except Exception:
                traceback.print_exc()
                raise


if __name__ == "__main__":
    main()
