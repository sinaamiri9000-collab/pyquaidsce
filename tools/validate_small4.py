"""
Full-precision validation against ``bench/small4.log``.

The user ran ``deliverables/benchmark_small.do`` in Stata: 4 goods, 2000
observations, 2 demographics, censoring in all four shares, with every
coefficient dumped from Mata at 17 significant digits.  Four specifications
came back:

    RUN 1  quaidsce_c ... (default, i.e. fgnls)
    RUN 2  quaidsce_c ... method(nls)
    RUN 3  quaidsce_c ... method(ifgnls)
    RUN 4  quaidsce_c ... noquadratic
    RUN 5  quaidsce_c ... nocensor          -> errored out in Stata, see below

Unlike ``log/censor_v2.log`` this benchmark also contains the ``ELAS_INC`` /
``ELAS_UNCOMP`` / ``ELAS_COMP`` blocks, so it validates the elasticity code of
v2.0 as well as the estimator.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from typing import Dict, Tuple

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "src"))
from pyquaidsce import quaidsce  # noqa: E402
from tools.paths import BENCH as BENCH_PATH  # noqa: E402

BENCH = str(BENCH_PATH)
SHARES = ["sw1", "sw2", "sw4", "sw9"]
PRICES = ["p1", "p2", "p4", "p9"]
DEMOS = ["x1", "x2"]

RUNS = {
    1: dict(method="fgnls", quadratic=True, censor=True),
    2: dict(method="nls", quadratic=True, censor=True),
    3: dict(method="ifgnls", quadratic=True, censor=True),
    4: dict(method="fgnls", quadratic=False, censor=True),
    5: dict(method="fgnls", quadratic=True, censor=False),
}

_P = re.compile(r"^PARM(\d) (\S+)\s+(\S+)(?:\s+(\S+))?\s*$")
_L = re.compile(r"^LL(\d)\s+(\S+)\s*$")


def parse(path: str):
    coefs: Dict[int, Dict[str, Tuple[float, float]]] = {}
    lls: Dict[int, float] = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            m = _P.match(line)
            if m:
                r = int(m.group(1))
                se = float(m.group(4)) if m.group(4) is not None else np.nan
                coefs.setdefault(r, {})[m.group(2)] = (float(m.group(3)), se)
                continue
            m = _L.match(line)
            if m:
                lls[int(m.group(1))] = float(m.group(2))
    return coefs, lls


def report(run: int, res, ref: Dict[str, Tuple[float, float]], ll_ref: float,
           out_dir: str, label: str):
    mine_b = dict(zip(res.names, res.b))
    mine_se = dict(zip(res.names, res.se))
    rows = []
    for nm, (b_s, se_s) in ref.items():
        if nm not in mine_b:
            rows.append((nm, np.nan, b_s, np.nan, se_s))
            continue
        rows.append((nm, mine_b[nm], b_s, mine_se[nm], se_s))
    df = pd.DataFrame(rows, columns=["name", "b_py", "b_stata",
                                     "se_py", "se_stata"])
    df["block"] = df.name.str.split(":").str[0]
    df["abs_db"] = (df.b_py - df.b_stata).abs()
    df["rel_db"] = df.abs_db / np.maximum(df.b_stata.abs(), 1e-10)
    df["abs_dse"] = (df.se_py - df.se_stata).abs()

    print(f"\n{'=' * 78}\nRUN {run}  ({label})\n{'=' * 78}")
    print(f"  log-likelihood   python {res.llf:.10f}")
    print(f"                   stata  {ll_ref:.10f}")
    print(f"                   relative difference "
          f"{abs(res.llf - ll_ref) / abs(ll_ref):.3e}")
    print(f"  coefficients compared: {int(df.b_py.notna().sum())} "
          f"/ {len(df)}")
    print(f"\n  {'block':<12}{'n':>4}{'max |db|':>12}{'max rel db':>13}"
          f"{'max |dse|':>12}")
    for blk, gg in df.groupby("block", sort=False):
        print(f"  {blk:<12}{len(gg):>4}{gg.abs_db.max():>12.2e}"
              f"{gg.rel_db.max():>13.2e}{gg.abs_dse.max():>12.2e}")
    worst = df.nlargest(6, "rel_db")
    print("\n  worst 6 by relative deviation:")
    for _, r in worst.iterrows():
        print(f"    {r['name']:<24} py={r.b_py:+.14g}")
        print(f"    {'':<24} st={r.b_stata:+.14g}   rel={r.rel_db:.2e}")
    os.makedirs(out_dir, exist_ok=True)
    df.to_csv(os.path.join(out_dir, f"small4_run{run}.csv"), index=False)
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--runs", nargs="*", type=int, default=[1, 2, 3, 4])
    ap.add_argument("--predict", default="pr", choices=["pr", "xb"])
    ap.add_argument("--algorithm", default="gn", choices=["gn", "lm"])
    ap.add_argument("--start", default="zero", choices=["zero", "linear"])
    ap.add_argument("--vce", default="objective", choices=["objective", "final"])
    ap.add_argument("--stop-rule", default="standard", choices=["standard", "tight"])
    ap.add_argument("--out", default="out")
    a = ap.parse_args()

    coefs, lls = parse(os.path.join(BENCH, "small4.log"))
    df0 = pd.read_stata(os.path.join(BENCH, "small4.dta"))
    print(f"data: {len(df0)} obs; zeros per share: "
          + ", ".join(f"{c}={int((df0[c] == 0).sum())}" for c in SHARES))

    summary = []
    for run in a.runs:
        if run not in coefs:
            print(f"\nRUN {run}: not present in the log (Stata errored out)")
            continue
        cfg = RUNS[run]
        res = quaidsce(
            df0, shares=SHARES, prices=PRICES, expenditure="total",
            demographics=DEMOS, anot=10.0,
            method=cfg["method"], quadratic=cfg["quadratic"],
            censor=cfg["censor"], first_stage_predict=a.predict,
            algorithm=a.algorithm, start=a.start, vce_sigma=a.vce,
            stop_rule=a.stop_rule,
            verbose=False,
        )
        lab = (f"method={cfg['method']}"
               + ("" if cfg["quadratic"] else ", noquadratic")
               + ("" if cfg["censor"] else ", nocensor"))
        df = report(run, res, coefs[run], lls[run], a.out, lab)
        summary.append((run, lab, abs(res.llf - lls[run]) / abs(lls[run]),
                        df.abs_db.max(), df.abs_dse.max()))

    print(f"\n{'=' * 78}\nSUMMARY\n{'=' * 78}")
    print(f"{'run':>4}  {'spec':<28}{'rel d ll':>12}{'max |db|':>12}"
          f"{'max |dse|':>12}")
    for run, lab, dll, db, dse in summary:
        print(f"{run:>4}  {lab:<28}{dll:>12.2e}{db:>12.2e}{dse:>12.2e}")


if __name__ == "__main__":
    main()
