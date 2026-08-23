"""
Bridge between Stata and pyquaidsce using Stata Function Interface (sfi).

This module is called directly from within Stata's `pyquaidsce.ado`.
It reads data from Stata's active memory, estimates the censored QUAIDS model,
and populates Stata matrices, scalars, and macros for `ereturn post`.
"""

from __future__ import annotations

import os
import sys
from typing import Any, Optional

# Module-level state for the out-of-process Stata job.
_BOOT_STATE: dict[str, Any] = {}


def _stata_matrix_to_vector(sfi, name: str) -> Any:
    """Read a Stata matrix by name and return it flattened to a 1-D list."""
    rows = sfi.Matrix.getRowTotal(name)
    cols = sfi.Matrix.getColTotal(name)
    values = []
    for i in range(rows):
        for j in range(cols):
            values.append(float(sfi.Matrix.getAt(name, i, j)))
    if not values:
        raise ValueError(f"Stata matrix {name!r} is empty")
    return values


def _stata_matrix_to_array(sfi, name: str) -> Any:
    """Read a Stata matrix by name and return a nested-list 2-D array."""
    rows = sfi.Matrix.getRowTotal(name)
    cols = sfi.Matrix.getColTotal(name)
    if rows == 0 or cols == 0:
        raise ValueError(f"Stata matrix {name!r} is empty")
    return [
        [float(sfi.Matrix.getAt(name, i, j)) for j in range(cols)]
        for i in range(rows)
    ]


def run_from_stata(
    shares_str: str,
    prices_str: str,
    expenditure_str: str,
    demographics_str: str,
    anot: float,
    method: str = "ifgnls",
    algorithm: str = "gn",
    start: str = "zero",
    reps: int = 0,
    stop_rule: str = "standard",
    bootstrap_start: str = "zero",
    seed: int = -1,
    n_jobs: int = 1,
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
    first_stage_predict: str = "xb",
    strict_stata: bool = False,
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
    vce_sigma: str = "objective",
    initial_mat_name: str = "",
    sigma_initial_mat_name: str = "",
    tol: float = 1e-13,
    max_outer: int = 200,
    max_iter: int = 300,
    chunk: int = 2000,
    nrtol_stop: float = 1e-12,
    inner_nrtol_early: float = 1e-8,
    sigma_tol: float = 1e-11,
    boot_sigma_tol: float = 1e-7,
    gn_verbose: bool = False,
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

    # Optional warm-start matrices passed by name from Stata (sfi.Matrix).
    initial_vec = None
    if initial_mat_name.strip():
        initial_vec = _stata_matrix_to_vector(sfi, initial_mat_name.strip())
    sigma_init_mat = None
    if sigma_initial_mat_name.strip():
        sigma_init_mat = _stata_matrix_to_array(sfi, sigma_initial_mat_name.strip())

    # 4. Save state for potential out-of-process bootstrap
    _BOOT_STATE["df"] = df
    _BOOT_STATE["quaidsce_kwargs"] = dict(
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
        start=start,
        stop_rule=stop_rule,
        bootstrap_start=bootstrap_start,
        first_stage_predict=first_stage_predict,
        strict_stata=strict_stata,
        vce_sigma=vce_sigma,
        tol=float(tol),
        max_outer=int(max_outer),
        max_iter=int(max_iter),
        chunk=int(chunk),
        nrtol_stop=float(nrtol_stop),
        inner_nrtol_early=float(inner_nrtol_early),
        sigma_tol=float(sigma_tol),
        boot_sigma_tol=float(boot_sigma_tol),
        verbose=verbose,
        gn_verbose=bool(gn_verbose),
    )

    # 5. Run estimation
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
        start=start,
        reps=int(reps),
        seed=seed_arg,
        stop_rule=stop_rule,
        bootstrap_start=bootstrap_start,
        n_jobs=int(n_jobs),
        mp_context=mp_context or None,
        rep_timeout=rep_timeout,
        first_stage_predict=first_stage_predict,
        strict_stata=strict_stata,
        vce_sigma=vce_sigma,
        initial=initial_vec,
        sigma_initial=sigma_init_mat,
        tol=float(tol),
        max_outer=int(max_outer),
        max_iter=int(max_iter),
        chunk=int(chunk),
        nrtol_stop=float(nrtol_stop),
        inner_nrtol_early=float(inner_nrtol_early),
        sigma_tol=float(sigma_tol),
        boot_sigma_tol=float(boot_sigma_tol),
        verbose=verbose,
        gn_verbose=bool(gn_verbose),
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


# --------------------------------------------------------------------------- #
#  Out-of-process estimation helpers (called from pyquaidsce.ado)
# --------------------------------------------------------------------------- #

def launch_from_stata(
    shares_str: str,
    prices_str: str,
    expenditure_str: str,
    demographics_str: str,
    anot: float,
    method: str = "ifgnls",
    algorithm: str = "gn",
    start: str = "zero",
    reps: int = 0,
    stop_rule: str = "standard",
    bootstrap_start: str = "zero",
    seed: int = -1,
    n_jobs: int = 1,
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
    first_stage_predict: str = "xb",
    strict_stata: bool = False,
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
    vce_sigma: str = "objective",
    initial_mat_name: str = "",
    sigma_initial_mat_name: str = "",
    tol: float = 1e-13,
    max_outer: int = 200,
    max_iter: int = 300,
    chunk: int = 2000,
    nrtol_stop: float = 1e-12,
    inner_nrtol_early: float = 1e-8,
    sigma_tol: float = 1e-11,
    boot_sigma_tol: float = 1e-7,
    gn_verbose: bool = False,
    touse_var: str = "_touse",
) -> None:
    """Export the Stata sample and launch the complete estimation asynchronously."""
    try:
        import sfi
    except ImportError:
        raise RuntimeError(
            "sfi module is only available when running inside Stata (version 16+)"
        )

    import pandas as pd

    shares = shares_str.split()
    prices = prices_str.split()
    demos = demographics_str.split() if demographics_str.strip() else []
    exp_var = expenditure_str.strip()
    cf_var = control_function.strip() or None
    selection_cf_var = selection_control_function.strip() or None
    selection_prices = (
        selection_prices_str.split() if selection_prices_specified else None
    )
    selection_covariates = (
        selection_covariates_str.split()
        if selection_covariates_specified else None
    )

    touse_mask = [bool(value) for value in sfi.Data.get(var=touse_var)]
    if not any(touse_mask):
        raise ValueError("No observations selected in estimation sample (all touse == 0).")

    all_vars = list(dict.fromkeys(
        shares + prices + demos + [exp_var]
        + ([cf_var] if cf_var is not None else [])
        + ([selection_cf_var] if selection_cf_var is not None else [])
        + ([] if selection_prices is None else selection_prices)
        + ([] if selection_covariates is None else selection_covariates)
    ))
    data_dict = {}
    for variable in all_vars:
        values = sfi.Data.get(var=variable)
        data_dict[variable] = [
            values[index] for index, selected in enumerate(touse_mask) if selected
        ]

    # Optional warm-start matrices passed by name from Stata (sfi.Matrix).
    initial_vec = None
    if initial_mat_name.strip():
        initial_vec = _stata_matrix_to_vector(sfi, initial_mat_name.strip())
    sigma_init_mat = None
    if sigma_initial_mat_name.strip():
        sigma_init_mat = _stata_matrix_to_array(sfi, sigma_initial_mat_name.strip())

    kwargs = dict(
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
        start=start,
        reps=int(reps),
        seed=int(seed) if int(seed) >= 0 else None,
        stop_rule=stop_rule,
        bootstrap_start=bootstrap_start,
        n_jobs=int(n_jobs),
        mp_context=mp_context or None,
        rep_timeout=rep_timeout,
        first_stage_predict=first_stage_predict,
        strict_stata=strict_stata,
        vce_sigma=vce_sigma,
        initial=initial_vec,
        sigma_initial=sigma_init_mat,
        tol=float(tol),
        max_outer=int(max_outer),
        max_iter=int(max_iter),
        chunk=int(chunk),
        nrtol_stop=float(nrtol_stop),
        inner_nrtol_early=float(inner_nrtol_early),
        sigma_tol=float(sigma_tol),
        boot_sigma_tol=float(boot_sigma_tol),
        verbose=verbose,
        gn_verbose=bool(gn_verbose),
    )
    _launch_job(pd.DataFrame(data_dict), kwargs)


def _launch_job(data: Any, kwargs: dict[str, Any]) -> None:
    """Serialize and launch one complete point-estimate/bootstrap job."""
    import pickle
    import subprocess
    import tempfile

    previous = _BOOT_STATE.get("proc")
    if previous is not None and previous.poll() is None:
        previous.terminate()
        try:
            previous.wait(timeout=5)
        except Exception:
            previous.kill()
    _cleanup_boot_files()

    tmpdir = tempfile.gettempdir()
    pid = os.getpid()
    job_path = os.path.join(tmpdir, f"pyq_job_{pid}.pkl")
    result_path = os.path.join(tmpdir, f"pyq_result_{pid}.pkl")
    status_path = os.path.join(tmpdir, f"pyq_status_{pid}.txt")

    for path in (
        job_path,
        result_path,
        status_path,
        result_path + ".tmp",
        result_path + ".error",
    ):
        if os.path.exists(path):
            try:
                os.unlink(path)
            except OSError:
                pass

    with open(job_path, "wb") as file_handle:
        pickle.dump(
            {
                "data": data,
                "kwargs": kwargs,
                "result_path": result_path,
                "status_path": status_path,
            },
            file_handle,
        )

    proc = subprocess.Popen(
        [sys.executable, "-m", "pyquaidsce.bootstrap_runner", "--input", job_path],
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )
    _BOOT_STATE.update(
        proc=proc,
        result_path=result_path,
        status_path=status_path,
        job_path=job_path,
    )


def launch_bootstrap(
    reps: int,
    seed: int = -1,
    n_jobs: int = 1,
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
) -> None:
    """Save a bootstrap job to a temp file and launch a subprocess.

    Uses the DataFrame and kwargs saved by the most recent
    :func:`run_from_stata` call (stored in ``_BOOT_STATE``).
    Returns immediately — Stata's GUI thread is never blocked.
    """
    if "df" not in _BOOT_STATE or "quaidsce_kwargs" not in _BOOT_STATE:
        raise RuntimeError(
            "No point estimate available in memory.  "
            "run_from_stata() must be called before launch_bootstrap()."
        )

    # Compatibility path for callers that first used run_from_stata().
    kwargs = _BOOT_STATE["quaidsce_kwargs"].copy()
    kwargs["reps"] = int(reps)
    kwargs["seed"] = int(seed) if int(seed) >= 0 else None
    kwargs["n_jobs"] = int(n_jobs)
    kwargs["mp_context"] = mp_context or None
    kwargs["rep_timeout"] = rep_timeout

    _launch_job(_BOOT_STATE["df"], kwargs)


def poll_bootstrap() -> None:
    """Check progress of the background bootstrap subprocess.

    Sets Stata scalars and macros that the ado polling loop reads:

    * ``scalar(_pyq_boot_done)`` — 1 when results are ready
    * ``scalar(_pyq_boot_err)``  — 1 when the subprocess has failed
    * ``local _pyq_boot_msg``    — latest progress line
    * ``local _pyq_boot_errmsg`` — error message (when ``_pyq_boot_err == 1``)
    """
    import sfi

    result_path = _BOOT_STATE.get("result_path", "")
    status_path = _BOOT_STATE.get("status_path", "")
    error_path = result_path + ".error" if result_path else ""

    # ---- finished successfully? ------------------------------------------- #
    if result_path and os.path.exists(result_path):
        sfi.Scalar.setValue("_pyq_boot_done", 1)
        sfi.Scalar.setValue("_pyq_boot_err", 0)
        return

    # ---- error file written? ---------------------------------------------- #
    if error_path and os.path.exists(error_path):
        try:
            with open(error_path, "r", encoding="utf-8") as fh:
                errmsg = fh.readline().strip()
        except OSError:
            errmsg = "unknown error"
        sfi.Scalar.setValue("_pyq_boot_done", 0)
        sfi.Scalar.setValue("_pyq_boot_err", 1)
        sfi.Macro.setLocal("_pyq_boot_errmsg", errmsg)
        return

    # ---- subprocess exited without output? -------------------------------- #
    proc = _BOOT_STATE.get("proc")
    if proc is not None and proc.poll() is not None:
        sfi.Scalar.setValue("_pyq_boot_done", 0)
        sfi.Scalar.setValue("_pyq_boot_err", 1)
        sfi.Macro.setLocal(
            "_pyq_boot_errmsg",
            f"Estimation process exited with code {proc.returncode}",
        )
        return

    # ---- still running — relay progress ----------------------------------- #
    sfi.Scalar.setValue("_pyq_boot_done", 0)
    sfi.Scalar.setValue("_pyq_boot_err", 0)
    if status_path and os.path.exists(status_path):
        try:
            with open(status_path, "r", encoding="utf-8") as fh:
                msg = fh.read().strip()
            if msg:
                sfi.Macro.setLocal("_pyq_boot_msg", msg)
        except OSError:
            pass


def load_stata_results(
    b_mat_name: str,
    v_mat_name: str,
    elas_i_name: str,
    elas_u_name: str,
    elas_c_name: str,
) -> None:
    """Load all completed background-estimation results into Stata."""
    import pickle

    import numpy as np
    import sfi

    result_path = _BOOT_STATE.get("result_path", "")
    if not result_path or not os.path.exists(result_path):
        raise FileNotFoundError("Estimation result file not found.")

    with open(result_path, "rb") as file_handle:
        result = pickle.load(file_handle)

    names = list(result["names"])
    b = np.asarray(result["b"], dtype=float)
    V = np.asarray(result["V"], dtype=float)
    shares = list(result["shares"])

    sfi.Matrix.create(b_mat_name, 1, b.size, 0.0)
    for column, value in enumerate(b):
        sfi.Matrix.storeAt(b_mat_name, 0, column, float(value))
    sfi.Matrix.setColNames(b_mat_name, names)

    sfi.Matrix.create(v_mat_name, V.shape[0], V.shape[1], 0.0)
    for row in range(V.shape[0]):
        for column in range(V.shape[1]):
            sfi.Matrix.storeAt(v_mat_name, row, column, float(V[row, column]))
    sfi.Matrix.setColNames(v_mat_name, names)
    sfi.Matrix.setRowNames(v_mat_name, names)

    elasticity_specs = (
        (elas_i_name, np.atleast_2d(result["elas_i"])),
        (elas_u_name, np.asarray(result["elas_u"], dtype=float)),
        (elas_c_name, np.asarray(result["elas_c"], dtype=float)),
    )
    for matrix_name, values in elasticity_specs:
        sfi.Matrix.create(matrix_name, values.shape[0], values.shape[1], 0.0)
        for row in range(values.shape[0]):
            for column in range(values.shape[1]):
                sfi.Matrix.storeAt(matrix_name, row, column, float(values[row, column]))
        sfi.Matrix.setColNames(matrix_name, shares)
        if values.shape[0] == len(shares):
            sfi.Matrix.setRowNames(matrix_name, shares)

    for name in (
        "nobs", "llf", "anot", "ndemo", "converged", "n_outer", "n_gn",
        "boot_reps_ok", "boot_reps_requested",
    ):
        sfi.Scalar.setValue("r_" + name, result[name])
    sfi.Macro.setLocal("model_title", result["title"])
    _cleanup_boot_files()


def load_bootstrap_results(v_mat_name: str) -> None:
    """Read a compatibility bootstrap job and store only its covariance."""
    import pickle

    import numpy as np
    import sfi

    result_path = _BOOT_STATE.get("result_path", "")
    if not result_path or not os.path.exists(result_path):
        raise FileNotFoundError("Bootstrap result file not found.")

    with open(result_path, "rb") as fh:
        result = pickle.load(fh)

    V = np.asarray(result["V"], dtype=float)
    names = list(result["names"])
    k = V.shape[0]

    sfi.Matrix.create(v_mat_name, k, k, 0.0)
    for i in range(k):
        for j in range(k):
            sfi.Matrix.storeAt(v_mat_name, i, j, float(V[i, j]))
    sfi.Matrix.setColNames(v_mat_name, names)
    sfi.Matrix.setRowNames(v_mat_name, names)

    sfi.Scalar.setValue("r_boot_reps_ok", int(result["reps_ok"]))
    sfi.Scalar.setValue(
        "r_boot_reps_requested", int(result["reps_requested"])
    )

    _cleanup_boot_files()


def kill_bootstrap() -> None:
    """Terminate the background bootstrap subprocess and clean up."""
    proc = _BOOT_STATE.get("proc")
    if proc is not None and proc.poll() is None:
        proc.terminate()
        try:
            proc.wait(timeout=5)
        except Exception:
            if hasattr(proc, "kill"):
                proc.kill()
    _cleanup_boot_files()


def _cleanup_boot_files() -> None:
    """Remove temporary IPC files created by a background estimation job."""
    for key in ("job_path", "result_path", "status_path"):
        path = _BOOT_STATE.get(key, "")
        if path:
            for suffix in ("", ".tmp", ".error"):
                p = path + suffix
                if os.path.exists(p):
                    try:
                        os.unlink(p)
                    except OSError:
                        pass
    _BOOT_STATE.pop("proc", None)
    _BOOT_STATE.pop("result_path", None)
    _BOOT_STATE.pop("status_path", None)
    _BOOT_STATE.pop("job_path", None)
