"""Compare full-precision Stata and Python outputs for the controlled benchmark."""
from __future__ import annotations

from pathlib import Path
import argparse
import json
import re

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
P_LINE = re.compile(r"^BENCHPARM\s+(\S+)\s+(\S+)\s+(\S+)\s*$")
R_LINE = re.compile(r"^BENCHRUNTIME\s+(\S+)\s*$")
LL_LINE = re.compile(r"^BENCHLL\s+(\S+)\s*$")


def parse_stata(path: Path):
    rows, runtime, ll = [], None, None
    with path.open(errors="replace") as f:
        for raw in f:
            line = raw.strip()
            m = P_LINE.match(line)
            if m:
                rows.append((m.group(1), float(m.group(2)), float(m.group(3))))
                continue
            m = R_LINE.match(line)
            if m:
                runtime = float(m.group(1))
                continue
            m = LL_LINE.match(line)
            if m:
                ll = float(m.group(1))
    return pd.DataFrame(rows, columns=["name", "stata_value", "stata_se"]), runtime, ll


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stata-log", default=str(RESULTS / "stata_benchmark.log"))
    ap.add_argument("--python-csv", default=str(RESULTS / "python_results.csv"))
    ap.add_argument("--python-runtime", default=str(RESULTS / "python_runtime.json"))
    ap.add_argument(
        "--same-machine",
        action="store_true",
        help="Assert that the stored Stata and Python timings were produced on the same machine.",
    )
    args = ap.parse_args()

    py = pd.read_csv(args.python_csv)
    st, st_runtime, st_ll = parse_stata(Path(args.stata_log))
    if st.empty:
        raise SystemExit("No BENCHPARM lines found in the Stata log.")
    if st_runtime is None or st_ll is None:
        raise SystemExit("Stata runtime or log-likelihood marker is missing.")

    d = py.merge(st, on="name", how="outer", indicator=True)
    if not (d["_merge"] == "both").all():
        print(d[d["_merge"] != "both"][["name", "_merge"]].to_string(index=False))
        raise SystemExit("Returned-name mismatch between Stata and Python.")

    d["abs_value_diff"] = (d.python_value - d.stata_value).abs()
    d["abs_se_diff"] = (d.python_se - d.stata_se).abs()
    d["rel_value_diff"] = d.abs_value_diff / np.maximum(d.stata_value.abs(), 1e-12)
    d["is_elasticity"] = d.block.isin(["ELAS_INC", "ELAS_UNCOMP", "ELAS_COMP"])
    d["is_first_stage"] = d.block.eq("tau")
    d["is_structural"] = d.block.isin(["alpha", "beta", "gamma", "lambda", "delta", "eta", "rho"])

    meta = json.loads(Path(args.python_runtime).read_text())
    py_runtime = float(meta["runtime_seconds"])
    speed_ratio = st_runtime / py_runtime
    ll_rel = abs(float(meta["log_likelihood"]) - st_ll) / abs(st_ll)

    structural = d[d.is_structural]
    first_stage = d[d.is_first_stage]
    elasticities = d[d.is_elasticity]

    note = (
        "The stored Stata and Python timings were recorded on the same machine. "
        "The ratio is therefore reported as a same-machine wall-clock speedup. "
        "The Python-side log for this recorded run does not contain an explicit BLAS/OpenMP "
        "thread-limit field, so the ratio should not be interpreted as a per-core efficiency comparison."
        if args.same_machine
        else
        "The script was not told that the two timings came from the same machine. "
        "Use --same-machine only when that condition is known to be true."
    )

    summary = {
        "values_compared": int(len(d)),
        "max_abs_structural_parameter_difference": float(structural.abs_value_diff.max()),
        "max_abs_first_stage_parameter_difference": float(first_stage.abs_value_diff.max()),
        "max_abs_elasticity_difference": float(elasticities.abs_value_diff.max()),
        "max_abs_standard_error_difference_nonelasticity": float(d.loc[~d.is_elasticity, "abs_se_diff"].max()),
        "python_log_likelihood": float(meta["log_likelihood"]),
        "stata_log_likelihood": float(st_ll),
        "relative_log_likelihood_difference": float(ll_rel),
        "python_runtime_seconds": py_runtime,
        "stata_runtime_seconds": float(st_runtime),
        "runtime_ratio_stata_over_python": float(speed_ratio),
        "runtime_ratio_is_same_machine": bool(args.same_machine),
        "note": note,
    }

    print("\nCONTROLLED CQUAIDS / IFGNLS NUMERICAL COMPARISON")
    print("--------------------------------------------------")
    print(f"Values compared                         : {len(d)}")
    print(f"Max |structural parameter difference|   : {summary['max_abs_structural_parameter_difference']:.6g}")
    print(f"Max |first-stage parameter difference|  : {summary['max_abs_first_stage_parameter_difference']:.6g}")
    print(f"Max |elasticity difference|             : {summary['max_abs_elasticity_difference']:.6g}")
    print(f"Max |SE difference|, non-elasticities   : {summary['max_abs_standard_error_difference_nonelasticity']:.6g}")
    print(f"Relative log-likelihood difference      : {ll_rel:.6g}")
    print(f"Stata runtime (seconds)                 : {st_runtime:.6f}")
    print(f"Python runtime (seconds)                : {py_runtime:.6f}")
    if args.same_machine:
        print(f"Same-machine wall-clock speedup         : {speed_ratio:.2f}x")
    else:
        print(f"Runtime ratio                           : {speed_ratio:.2f}x [same-machine status not asserted]")

    d = d.drop(columns=["_merge"])
    d.sort_values("abs_value_diff", ascending=False).to_csv(
        RESULTS / "stata_python_comparison.csv", index=False
    )
    (RESULTS / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
