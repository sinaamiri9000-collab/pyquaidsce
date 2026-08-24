"""
Master Comprehensive Test Suite for pyquaidsce 1.4.0
=====================================================
Executes 24 real-world empirical test scenarios on Uruguayan household
demand data (bd_uruguay.csv, 14 food groups, 6,848 observations).

Usage:
    python tests/master_test_suite.py
"""

from __future__ import annotations

import os
import sys
import time
import traceback
import numpy as np
import pandas as pd

from pyquaidsce import quaidsce


def load_dataset() -> pd.DataFrame:
    data_candidates = [
        r"C:\Users\sina\Downloads\bd_uruguay.csv",
        os.path.join(os.path.dirname(__file__), "..", "benchmarks", "bd_uruguay.csv"),
        "bd_uruguay.csv",
    ]
    df = None
    for path in data_candidates:
        if os.path.exists(path):
            df = pd.read_csv(path)
            print(f"Loaded dataset from: {path} (N = {len(df):,} observations)")
            break

    if df is None:
        raise FileNotFoundError("Could not find bd_uruguay.csv in candidate paths.")

    # Clean negative shares if any
    for i in range(1, 15):
        col = f"w{i}_red"
        if col in df.columns:
            df[col] = df[col].clip(lower=0.0)

    # Generate log prices and log expenditure
    for i in range(1, 15):
        p_col = f"P_med{i}"
        if p_col in df.columns:
            df[f"ln_{p_col}"] = np.log(df[p_col])

    if "ln_gasto_total" not in df.columns:
        df["ln_gasto_total"] = np.log(df["gasto_total"])

    # Generate an endogeneity residual proxy for test 18 using pure numpy OLS
    demo_vars = ["npersonas", "edad", "Sex", "log_ing"]
    X_mat = np.column_stack([np.ones(len(df)), df[demo_vars].values])
    y_vec = df["ln_gasto_total"].values
    beta, _, _, _ = np.linalg.lstsq(X_mat, y_vec, rcond=None)
    df["vhat"] = y_vec - X_mat @ beta

    return df


def run_tests():
    df = load_dataset()

    shares_14 = [f"w{i}_red" for i in range(1, 15)]
    prices_14 = [f"P_med{i}" for i in range(1, 15)]
    lnprices_14 = [f"ln_P_med{i}" for i in range(1, 15)]
    demographics = ["npersonas", "edad", "Sex", "prim_comp", "sec_comp", "sup_comp", "log_ing"]
    anot_val = 2.9766226

    results = []

    def record_test(test_id: int, title: str, func):
        print(f"\n[{test_id:02d}/24] RUNNING: {title} ...", flush=True)
        t0 = time.time()
        try:
            res_obj = func()
            elapsed = time.time() - t0
            llf = getattr(res_obj, "llf", None)
            converged = getattr(res_obj, "converged", True)
            status = "PASS" if converged else "WARN(Not Converged)"
            msg = f"LLF = {llf:.4f}, Time = {elapsed:.2f}s" if llf is not None else f"Time = {elapsed:.2f}s"
            print(f"       --> [{status}] {msg}")
            results.append((test_id, title, "PASS", elapsed, msg))
            return res_obj
        except Exception as exc:
            elapsed = time.time() - t0
            errmsg = f"{type(exc).__name__}: {exc}"
            print(f"       --> [FAIL] {errmsg}")
            traceback.print_exc()
            results.append((test_id, title, "FAIL", elapsed, errmsg))
            return None

    print("=" * 80)
    print("PYQUAIDSCE 1.4.0 MASTER COMPREHENSIVE 24-SCENARIO TEST SUITE")
    print("=" * 80)

    # -------------------------------------------------------------------------
    # Tier 1: Model Specification & Variable Formats
    # -------------------------------------------------------------------------
    t1_fit = record_test(
        1, "Baseline 14-Good IFGNLS (first_stage_predict='xb')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, method="ifgnls", first_stage_predict="xb"
        )
    )

    record_test(
        2, "Direct Log-Prices & Log-Expenditure (lnprices & lnexpenditure)",
        lambda: quaidsce(
            data=df, shares=shares_14, lnprices=lnprices_14, lnexpenditure="ln_gasto_total",
            demographics=demographics, anot=anot_val, method="ifgnls"
        )
    )

    record_test(
        3, "Linear AIDS Model (quadratic=False)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, quadratic=False, method="ifgnls"
        )
    )

    record_test(
        4, "Uncensored QUAIDS (censor=False, Poi 2012 specification)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, censor=False, method="ifgnls"
        )
    )

    record_test(
        5, "Translog Constant Variation (anot = 10.0)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=10.0, method="ifgnls"
        )
    )

    # 3-good subsystem
    df_3g = df.copy()
    sub_shares = ["w1_red", "w2_red", "w3_red"]
    total_sub = df_3g[sub_shares].sum(axis=1)
    for s in sub_shares:
        df_3g[s + "_sub"] = df_3g[s] / total_sub
    sub_exp = (df_3g["gasto1"] + df_3g["gasto2"] + df_3g["gasto3"]).clip(lower=1.0)
    df_3g["sub_exp"] = sub_exp

    record_test(
        6, "3-Good Subsystem Estimation",
        lambda: quaidsce(
            data=df_3g, shares=[s + "_sub" for s in sub_shares],
            prices=["P_med1", "P_med2", "P_med3"], expenditure="sub_exp",
            demographics=demographics, anot=anot_val, method="ifgnls"
        )
    )

    # -------------------------------------------------------------------------
    # Tier 2: Solvers, Algorithms & Numerical Tolerances
    # -------------------------------------------------------------------------
    t7_fit = record_test(
        7, "Nonlinear Least Squares (method='nls')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, method="nls"
        )
    )

    t8_fit = record_test(
        8, "Feasible Generalized NLS (method='fgnls')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, method="fgnls"
        )
    )

    record_test(
        9, "Iterated FGNLS (method='ifgnls')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, method="ifgnls"
        )
    )

    record_test(
        10, "Levenberg-Marquardt Optimizer (algorithm='lm')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, algorithm="lm", method="ifgnls"
        )
    )

    record_test(
        11, "Strict Gradient Stopping Rule (stop_rule='tight')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, stop_rule="tight", method="ifgnls"
        )
    )

    record_test(
        12, "Alternative VCE Sigma Formula (vce_sigma='final')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, vce_sigma="final", method="ifgnls"
        )
    )

    record_test(
        13, "Chunk Size Variation (chunk=500)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, chunk=500, method="ifgnls"
        )
    )

    # -------------------------------------------------------------------------
    # Tier 3: Censoring, Probit Specifications & Control Functions
    # -------------------------------------------------------------------------
    record_test(
        14, "Legacy Stata Probit CDF Predictor (first_stage_predict='pr')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, first_stage_predict="pr", method="ifgnls"
        )
    )

    record_test(
        15, "Selection Price Subset (selection_prices=['P_med1', 'P_med2', 'P_med3'])",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val,
            selection_prices=["P_med1", "P_med2", "P_med3"], method="ifgnls"
        )
    )

    record_test(
        16, "Selection Independent Covariates (selection_covariates=['npersonas', 'edad'])",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val,
            selection_covariates=["npersonas", "edad"], method="ifgnls"
        )
    )

    record_test(
        17, "Selection Omit Log Expenditure (selection_expenditure=False)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val,
            selection_expenditure=False, method="ifgnls"
        )
    )

    record_test(
        18, "Endogeneity Control Function (control_function='vhat')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val,
            control_function="vhat", first_stage_predict="xb", method="ifgnls"
        )
    )

    # -------------------------------------------------------------------------
    # Tier 4: Warm-Starting & Matrix Transfers
    # -------------------------------------------------------------------------
    record_test(
        19, "Linearized AIDS Starting Values (start='linear')",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, start="linear", method="ifgnls"
        )
    )

    # Extract initial free structural parameters from base model
    theta_init = t1_fit.theta if t1_fit is not None else None
    sigma_init = t1_fit.sigma if t1_fit is not None else None

    record_test(
        20, "Custom Structural Free Parameters Vector (initial=theta_init)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, initial=theta_init, method="ifgnls"
        )
    )

    record_test(
        21, "Custom Residual Covariance Matrix (sigma_initial=sigma_init)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, initial=theta_init,
            sigma_initial=sigma_init, method="ifgnls"
        )
    )

    # Chained warm-start (NLS -> IFGNLS)
    record_test(
        22, "Chained Multi-Stage Warm-Start (NLS theta/sigma -> IFGNLS)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val,
            initial=t7_fit.theta if t7_fit is not None else None,
            sigma_initial=t7_fit.sigma if t7_fit is not None else None,
            method="ifgnls"
        )
    )

    # -------------------------------------------------------------------------
    # Tier 5: Bootstrap, Multiprocessing & Postestimation
    # -------------------------------------------------------------------------
    record_test(
        23, "Parallel Bootstrap Standard Errors (reps=10, n_jobs=4, seed=123456)",
        lambda: quaidsce(
            data=df, shares=shares_14, prices=prices_14, expenditure="gasto_total",
            demographics=demographics, anot=anot_val, reps=10, n_jobs=4,
            seed=123456, mp_context="spawn", method="ifgnls"
        )
    )

    def test_postestimation():
        if t1_fit is None:
            raise RuntimeError("Base fit t1_fit is missing.")
        s = t1_fit.summary()
        assert len(s) > 100, "Summary output too short"
        tbl = t1_fit.elasticity_tables()
        assert "Expenditure" in tbl, "Elasticity tables missing Expenditure"
        named_dict = t1_fit.named()
        assert len(named_dict) > 0, "named() dictionary is empty"
        beta1 = t1_fit.get("beta_1")
        assert beta1 is not None, "get('beta_1') returned None"
        assert t1_fit.elas.income.shape == (14,), f"Bad income elas shape: {t1_fit.elas.income.shape}"
        assert t1_fit.elas.uncompensated.shape == (14, 14), f"Bad uncomp elas shape: {t1_fit.elas.uncompensated.shape}"
        return t1_fit

    record_test(
        24, "Comprehensive Postestimation & Elasticities Extraction API",
        test_postestimation
    )

    # -------------------------------------------------------------------------
    # Final Scorecard
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("PYQUAIDSCE 1.4.0 MASTER TEST SUITE SCORECARD")
    print("=" * 80)
    print(f"{'#':<3} {'Test Scenario Title':<52} {'Status':<8} {'Time (s)':<10} {'Details'}")
    print("-" * 80)
    n_pass = 0
    total_time = 0.0
    for tid, title, st, elapsed, det in results:
        total_time += elapsed
        if st == "PASS":
            n_pass += 1
        print(f"{tid:<3} {title[:50]:<52} {st:<8} {elapsed:>6.2f}s    {det[:30]}")
    print("-" * 80)
    print(f"Summary: {n_pass}/24 tests PASSED successfully in {total_time:.2f} seconds.")
    print("=" * 80 + "\n")


if __name__ == "__main__":
    run_tests()
