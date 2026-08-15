"""
Bridge between Stata and pyquaidsce using Stata Function Interface (sfi).

This module is called directly from within Stata's `pyquaidsce.ado`.
It reads data from Stata's active memory, estimates the censored QUAIDS model,
and populates Stata matrices, scalars, and macros for `ereturn post`.
"""

from __future__ import annotations

import sys
from typing import List, Optional


def run_from_stata(
    shares_str: str,
    prices_str: str,
    expenditure_str: str,
    demographics_str: str,
    anot: float,
    method: str = "ifgnls",
    algorithm: str = "gn",
    reps: int = 0,
    seed: int = -1,
    n_jobs: int = 1,
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
    first_stage_predict: str = "pr",
    strict_stata: bool = True,
    quadratic: bool = True,
    censor: bool = True,
    is_lnprices: bool = False,
    is_lnexp: bool = False,
    control_function: str = "",
    selection_control_function: str = "",
    selection_prices_str: str = "",
    selection_prices_specified: bool = False,
    selection_covariates_str: str = "",
    selection_covariates_specified: bool = False,
    selection_expenditure: bool = True,
    verbose: bool = True,
    b_mat_name: str = "__pyq_b",
    v_mat_name: str = "__pyq_V",
    elas_i_name: str = "__pyq_elas_i",
    elas_u_name: str = "__pyq_elas_u",
    elas_c_name: str = "__pyq_elas_c",
    touse_var: str = "_touse",
) -> None:
    try:
        import sfi
    except ImportError:
        raise RuntimeError("sfi module is only available when running inside Stata (version 16+)")

    import numpy as np
    import pandas as pd

    from .estimator import quaidsce

    # 1. Parse variable lists
    shares = shares_str.split()
    prices = prices_str.split()
    demos = demographics_str.split() if demographics_str.strip() else []
    exp_var = expenditure_str.strip()
    cf_var = control_function.strip() or None
    selection_cf_var = selection_control_function.strip() or None
    selection_prices = (
        selection_prices_str.split()
        if selection_prices_specified else None
    )
    selection_covariates = (
        selection_covariates_str.split()
        if selection_covariates_specified else None
    )

    # 2. Get estimation sample mask from Stata
    touse_raw = sfi.Data.get(var=touse_var)
    touse_mask = [bool(t) for t in touse_raw]
    n_selected = sum(touse_mask)

    if n_selected == 0:
        raise ValueError("No observations selected in estimation sample (all touse == 0).")

    # 3. Read data columns from Stata memory
    all_vars = list(dict.fromkeys(
        shares + prices + demos + [exp_var]
        + ([cf_var] if cf_var is not None else [])
        + ([selection_cf_var] if selection_cf_var is not None else [])
        + ([] if selection_prices is None else selection_prices)
        + ([] if selection_covariates is None else selection_covariates)
    ))
    data_dict = {}
    for v in all_vars:
        col_vals = sfi.Data.get(var=v)
        # Filter to touse observations
        data_dict[v] = [col_vals[i] for i, m in enumerate(touse_mask) if m]

    df = pd.DataFrame(data_dict)

    # 4. Run estimation
    seed_arg = int(seed) if seed >= 0 else None
    res = quaidsce(
        data=df,
        shares=shares,
        prices=None if is_lnprices else prices,
        lnprices=prices if is_lnprices else None,
        expenditure=None if is_lnexp else exp_var,
        lnexpenditure=exp_var if is_lnexp else None,
        demographics=demos if demos else None,
        control_function=cf_var,
        selection_control_function=selection_cf_var,
        selection_prices=selection_prices,
        selection_covariates=selection_covariates,
        selection_expenditure=bool(selection_expenditure),
        anot=float(anot),
        quadratic=quadratic,
        censor=censor,
        method=method,
        algorithm=algorithm,
        reps=int(reps),
        seed=seed_arg,
        n_jobs=int(n_jobs),
        mp_context=mp_context or None,
        rep_timeout=rep_timeout,
        first_stage_predict=first_stage_predict,
        strict_stata=strict_stata,
        verbose=verbose,
    )

    # 5. Populate Stata scalars
    sfi.Scalar.setValue("r_nobs", int(res.nobs))
    sfi.Scalar.setValue("r_llf", float(res.llf))
    sfi.Scalar.setValue("r_anot", float(res.anot))
    sfi.Scalar.setValue("r_ndemo", int(res.spec.ndemo))
    sfi.Scalar.setValue("r_converged", 1 if res.converged else 0)
    sfi.Scalar.setValue("r_n_outer", int(res.n_outer))
    sfi.Scalar.setValue("r_n_gn", int(res.n_gn))

    # 6. Store coefficient vector b and covariance matrix V
    # If bootstrap was run, use bootstrap covariance matrix
    V_matrix = res.boot.V if (res.boot is not None) else res.V
    k_params = len(res.b)

    sfi.Matrix.create(b_mat_name, 1, k_params, 0.0)
    for j, val in enumerate(res.b):
        sfi.Matrix.storeAt(b_mat_name, 0, j, float(val))
    sfi.Matrix.setColNames(b_mat_name, res.names)

    sfi.Matrix.create(v_mat_name, k_params, k_params, 0.0)
    for i in range(k_params):
        for j in range(k_params):
            sfi.Matrix.storeAt(v_mat_name, i, j, float(V_matrix[i, j]))
    sfi.Matrix.setColNames(v_mat_name, res.names)
    sfi.Matrix.setRowNames(v_mat_name, res.names)

    # 7. Store elasticity matrices
    n_goods = len(shares)

    # Expenditure (income) elasticities (1 x n)
    sfi.Matrix.create(elas_i_name, 1, n_goods, 0.0)
    for j, val in enumerate(res.elas.income):
        sfi.Matrix.storeAt(elas_i_name, 0, j, float(val))
    sfi.Matrix.setColNames(elas_i_name, shares)

    # Marshallian (uncompensated) price elasticities (n x n)
    sfi.Matrix.create(elas_u_name, n_goods, n_goods, 0.0)
    for i in range(n_goods):
        for j in range(n_goods):
            sfi.Matrix.storeAt(elas_u_name, i, j, float(res.elas.uncompensated[i, j]))
    sfi.Matrix.setColNames(elas_u_name, shares)
    sfi.Matrix.setRowNames(elas_u_name, shares)

    # Hicksian (compensated) price elasticities (n x n)
    sfi.Matrix.create(elas_c_name, n_goods, n_goods, 0.0)
    for i in range(n_goods):
        for j in range(n_goods):
            sfi.Matrix.storeAt(elas_c_name, i, j, float(res.elas.compensated[i, j]))
    sfi.Matrix.setColNames(elas_c_name, shares)
    sfi.Matrix.setRowNames(elas_c_name, shares)

    # 8. Set title macro
    sfi.Macro.setLocal("model_title", res._title())
