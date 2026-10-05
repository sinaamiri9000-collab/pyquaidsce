"""Public helpers, output behavior, and interface contracts.

The Python estimator is real in bridge tests; only Stata's unavailable sfi
transport is simulated. R checks explicitly inspect the wrapper contract and
are never presented as native R execution.
"""
import ast
import contextlib
import io
import inspect
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time
import traceback
import types
from unittest.mock import patch
import warnings
import numpy as np
from scipy.stats import norm
from pyquaidsce import (quaidsce, first_stage, nlsur, probit, latent_shares,
                       augmented_latent_shares, fitted_share_derivatives,
                       fit_expenditure_reduced_form, full_vector, delta_matrix, sample_means)
from .cases import BASE_ID, PRICES, SHARES
from .oracles import Checks, build_data, independent_fitted, predict


def _blocked_worker(conn, task, payload, blas_threads):
    """Spawn-importable fault injection: no cooperative deadline checks."""
    payload["ready"].set()
    time.sleep(15)
    conn.close()


def contract_tests(out, original, cache, frames, only=None):
    from .run import ROOT, write_json
    records = []
    base = cache[BASE_ID]
    df, kw = frames[BASE_ID]
    d, design, index, _ = build_data(df, kw, base)

    def run(id, group, fn, expectation="public contract holds"):
        if only and id not in only:
            return
        check, log = Checks(), io.StringIO()
        start = time.perf_counter()
        row = dict(id=id, group=group, expectation=expectation, checks=[])
        with contextlib.redirect_stdout(log), contextlib.redirect_stderr(log), warnings.catch_warnings(record=True) as warns:
            warnings.simplefilter("always")
            try:
                fn(check)
            except Exception as exc:
                check.test("unexpected_exception", False, f"{type(exc).__name__}: {exc}", "no exception")
                row["exception"] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
            row["warnings"] = [str(x.message) for x in warns]
        row.update(checks=check.rows, elapsed_seconds=time.perf_counter()-start,
                   status="PASS" if all(x["passed"] for x in check.rows) else "FAIL")
        if log.getvalue():
            row["log"] = log.getvalue()
        write_json(out / "cases" / f"{id}.json", row)
        records.append(row)
        print(f"{id}: {row['status']}", flush=True)

    def no_constant(check):
        X = design[:, :2]
        p = probit((d.shares[:, 0] > 0).astype(float), X, add_constant=False)
        check.test("no_constant_coefficient_dimension", len(p.b) == X.shape[1])
        check.close("xb_without_constant", p.xb(X), X @ p.b)
    run("public-probit_no_constant", "public_helpers", no_constant,
        "ProbitResult.xb respects probit(add_constant=False)")

    def fs(check):
        r = first_stage(d.shares, d.lnp, d.lnexp, d.demo, design=design, layout=base.selection_layout)
        check.close("explicit_design_tau", r.tau, base.tau)
        check.close("canonical_cdf", r.cdf, norm.cdf(index))
        check.close("canonical_pdf", r.pdf, norm.pdf(index))
        check.close("canonical_du", r.du, index)
    run("public-first_stage", "public_helpers", fs)

    def nl(check):
        r = nlsur(d, base.spec, method="ifgnls", max_iter=300)
        check.close("public_nlsur_theta", r.theta, base.theta, atol=1e-6, rtol=1e-5)
        check.close("public_nlsur_sigma", r.sigma, base.sigma, atol=1e-6, rtol=1e-5)
        check.test("public_nlsur_converged", r.converged)
    run("public-nlsur", "public_helpers", nl)

    def means(check):
        r = sample_means(d, index, base.spec)
        for name in ("w", "lnp", "lnexp", "demo", "cdf", "pdf", "du", "control_function"):
            check.close("sample_mean_" + name, getattr(r, name), getattr(base.means, name))
    run("public-sample_means", "public_helpers", means)

    def latent(check):
        check.close("latent_direct_equation", latent_shares(base.theta, d, base.spec),
                    predict(base.coefs, d.lnp, d.lnexp, d.demo, d.a0, base.spec.quadratic), atol=2e-10, rtol=2e-10)
        cfres = cache["cf-both_distinct"] if "cf-both_distinct" in cache else base
        cdf, ckw = frames.get("cf-both_distinct", frames[BASE_ID])
        cd, _, _, _ = build_data(cdf, ckw, cfres)
        check.close("augmented_minus_latent", augmented_latent_shares(cfres.theta, cd, cfres.spec)
                    - latent_shares(cfres.theta, cd, cfres.spec),
                    cd.control_function[:, None] * cfres.coefs.cfcoef, atol=1e-12, rtol=1e-12)
    run("public-latent_shares", "public_helpers", latent)

    def delta(check):
        actual = delta_matrix(base.spec)
        numerical = np.zeros_like(actual)
        for k in range(len(base.theta)):
            plus, minus = base.theta.copy(), base.theta.copy()
            plus[k] += 1e-6
            minus[k] -= 1e-6
            numerical[:, k] = (full_vector(plus, base.spec) - full_vector(minus, base.spec)) / 2e-6
        check.close("delta_finite_difference", actual, numerical, atol=2e-8, rtol=2e-8)
    run("public-delta_mapping", "public_helpers", delta)

    def stale(check):
        stale_index = index + .1
        try:
            fitted_share_derivatives(base.theta, d, base.spec, tau=base.tau,
                                     layout=base.selection_layout, selection_index=stale_index)
        except ValueError as e:
            check.test("stale_design_rejected", "rebuilt" in str(e))
        else:
            check.test("stale_design_rejected", False)
    run("public-stale_derivative_design", "public_helpers", stale)

    def labels(check):
        rdf, _ = frames.get("cf-iv_one", frames[BASE_ID])
        iv = rdf[["iv1"]].to_numpy(float) if "iv1" in rdf else np.random.default_rng(1).normal(size=(len(rdf), 1))
        try:
            fit_expenditure_reduced_form(d.lnexp, d.lnp, d.demo, iv,
                outcome_name="ln(total)", price_names=["p1"], price_inputs_are_logs=False,
                demographic_names=["x1", "x2"], instrument_names=["iv1"])
        except ValueError:
            check.test("mismatched_labels_rejected", True)
        else:
            check.test("mismatched_labels_rejected", False, "returned mislabeled coefficient vector", "ValueError")
    run("public-reduced_form_labels", "public_helpers", labels,
        "reported regressor names must match coefficient dimensions")

    for bootstrap in (False, True):
        id = "output-summary-no_elasticities-" + ("bootstrap" if bootstrap else "analytic")
        target = cache.get("boot-serial") if bootstrap else base
        if target is None:
            continue
        def summary(check, target=target):
            text = target.summary(elasticities=False)
            bare = [x.split(":", 1)[-1] for x in target.spec.elas_names()]
            check.test("elasticity_rows_suppressed", not any(x in text for x in bare),
                       [x for x in bare if x in text], [])
        run(id, "outputs", summary, "summary(elasticities=False) suppresses elasticity coefficient rows")

    def get_missing(check):
        try:
            base.get("does_not_exist")
        except KeyError:
            check.test("missing_coefficient_KeyError", True)
        else:
            check.test("missing_coefficient_KeyError", False)
        try:
            base.reduced_form_table()
        except ValueError:
            check.test("missing_reduced_form_ValueError", True)
        else:
            check.test("missing_reduced_form_ValueError", False)
    run("output-missing_lookup", "outputs", get_missing)

    bootstrap_reference = cache.get("boot-serial")
    if bootstrap_reference is not None:
        for level_id, level in enumerate((-1, 101, float("nan"))):
            for helper in ("ci", "percentile_ci", "summary"):
                def invalid_level(check, helper=helper, level=level):
                    try:
                        if helper == "ci":
                            bootstrap_reference.boot.ci(bootstrap_reference.b, level)
                        elif helper == "percentile_ci":
                            bootstrap_reference.boot.percentile_ci(level)
                        else:
                            bootstrap_reference.summary(level=level)
                    except ValueError:
                        check.test("invalid_confidence_level_rejected", True)
                    else:
                        check.test("invalid_confidence_level_rejected", False, level, "finite confidence percentage within 0 to 100")
                run(f"output-invalid_level-{helper}-{level_id}", "outputs", invalid_level,
                    "invalid confidence percentage is rejected explicitly")
        for level in (0, 100):
            def boundary(check, level=level):
                boot = bootstrap_reference.boot
                lo, hi = boot.percentile_ci(level)
                if level == 0:
                    check.close("zero_percentile_lower", lo, np.median(boot.b_star, axis=0))
                    check.close("zero_percentile_upper", hi, lo)
                    lo, hi = boot.ci(bootstrap_reference.b, level)
                    check.close("zero_normal_lower", lo, bootstrap_reference.b)
                    check.close("zero_normal_upper", hi, bootstrap_reference.b)
                else:
                    check.close("full_percentile_lower", lo, boot.b_star.min(axis=0))
                    check.close("full_percentile_upper", hi, boot.b_star.max(axis=0))
            run(f"output-ci_boundary-{level}", "outputs", boundary,
                "well-defined confidence endpoints are boundary probes, not forced rejection tests")

    def nls_covariance(check):
        result = cache["core-cq-nls-gn-zero"]
        from pyquaidsce import jacobian_full, delta_matrix
        J = np.einsum("tmf,fk->tmk", jacobian_full(result.theta, d, result.spec), delta_matrix(result.spec))
        A = np.einsum("tik,til->kl", J, J)
        B = np.einsum("tik,ij,tjl->kl", J, result.sigma, J)
        inverse = np.linalg.inv(A)
        sandwich = inverse @ B @ inverse
        check.close("conditional_multivariate_NLS_covariance", result.V_est, sandwich, atol=1e-7, rtol=1e-5)
        se_ratio = np.sqrt(np.diag(result.V_est)) / np.sqrt(np.diag(sandwich))
        check.test("NLS_covariance_SE_ratio", np.allclose(se_ratio, 1, atol=1e-5, rtol=1e-5),
                   dict(min_ratio=float(se_ratio.min()), max_ratio=float(se_ratio.max())),
                   "OLS sandwich conditional on first stage; full generated-regressor inference still needs bootstrap")
    run("methodology-NLS_covariance", "methodology", nls_covariance,
        "identity-weighted NLS uses its conditional sandwich covariance, not efficient GLS covariance")

    def tolerance_contract(check):
        source = (ROOT / "src" / "pyquaidsce" / "nlsur.py").read_text()
        tree = ast.parse(source)
        fn = next(x for x in tree.body if isinstance(x, ast.FunctionDef) and x.name == "gauss_newton")
        branches = [x for x in ast.walk(fn) if isinstance(x, ast.If)
                    and isinstance(x.test, ast.Compare) and ast.unparse(x.test) == "stop_rule == 'standard'"]
        names = {x.id for branch in branches for node in branch.body for x in ast.walk(node) if isinstance(x, ast.Name)}
        check.test("standard_uses_public_tolerances", "tol" in names and "nrtol_stop" in names,
                   dict(names=sorted(names), source="src/pyquaidsce/nlsur.py:307"),
                   "documented tol and nrtol_stop affect standard stopping thresholds")
    run("numeric-standard_tolerance_contract", "source_contracts", tolerance_contract,
        "documented convergence tolerances are consulted under the default stop_rule")

    def watchdog(check):
        from pyquaidsce import bootstrap as bootmod
        ctx = mp.get_context("spawn")
        ready = ctx.Event()
        children_before = {x.pid for x in mp.active_children()}
        start = time.perf_counter()
        with patch.object(bootmod, "_one_rep_process", _blocked_worker):
            results = list(bootmod._watchdog_results(ctx, [(1, 123)], {"ready": ready}, 1, 5.0, 1))
        check.test("stuck_worker_started", ready.is_set())
        check.test("hard_watchdog_timeout", len(results) == 1 and "TimeoutError" in results[0][2], results,
                   "parent terminates a worker that never checks its deadline")
        elapsed = time.perf_counter() - start
        check.test("hard_watchdog_wall_bound", elapsed < 9, elapsed, 9)
        check.test("hard_watchdog_no_orphan", {x.pid for x in mp.active_children()} == children_before)
    run("lifecycle-hard_watchdog", "lifecycle", watchdog,
        "terminate an injected noncooperative blocked worker and reap it")

    # Supplemental smoke paths and tighter fits keep the original experiment
    # intact. They are recorded under new IDs, never used to overwrite a failure.
    from .cases import Case
    from .run import capture_fit
    extras = [
        Case("boot-positive_timeout", "bootstrap_supplemental",
             dict(reps=8, seed=16001, n_jobs=2, mp_context="spawn", rep_timeout=10),
             checks=["bootstrap", "bootstrap_compare"], reference="boot-serial",
             expectation="generous watchdog preserves all ordinary draws"),
        Case("boot-default_seed", "bootstrap_supplemental", dict(reps=2, seed=None),
             checks=["bootstrap"], expectation="seed=None supports a finite unseeded bootstrap; realized draws are archived"),
    ]
    if "fork" in mp.get_all_start_methods():
        extras.append(Case("boot-fork", "bootstrap_supplemental", dict(reps=2, seed=16001, n_jobs=2, mp_context="fork"),
                           checks=["bootstrap"], expectation="available fork context works with caller BLAS limited to one"))
    for method in ("nls", "fgnls", "ifgnls"):
        for algorithm in ("gn", "lm"):
            for start in ("zero", "linear"):
                extras.append(Case(f"followup-tight-{method}-{algorithm}-{start}", "followup_tight",
                                   dict(method=method, algorithm=algorithm, start=start, stop_rule="tight"),
                                   checks=["stationarity"],
                                   expectation="tighter stopping removes early-stop ambiguity; original cases remain unchanged"))
    for source_id in ("selection-no_prices", "selection-no_expenditure", "selection-intercept_only", "selection-independent_no_exp"):
        from .cases import scenarios
        original_case = next(x for x in scenarios() if x.id == source_id)
        extras.append(Case("followup-linear-" + source_id, "followup_selection",
                           dict(original_case.kwargs, start="linear"), original_case.fixture,
                           expectation="diagnostic alternative start; not a replacement for the original failure"))
    write_json(out / "supplemental-manifest.json", [dict(id=c.id, group=c.group, kwargs=c.kwargs,
               fixture=c.fixture, expectation=c.expectation) for c in extras])
    for case in extras:
        if only and case.id not in only:
            continue
        if case.reference and case.reference not in cache:
            continue
        row, _ = capture_fit(out, case, original, cache, frames, audit_timeout=20 if case.group.startswith("followup") else 90)
        records.append(row)
        print(f"{case.id}: {row['status']}", flush=True)

    r_source = (ROOT / "R" / "R" / "quaidsce.R").read_text()
    r_methods = (ROOT / "R" / "R" / "methods.R").read_text()
    def r_fitted(check):
        calls = "res_py$fitted_shares()" in r_source
        method_exists = callable(getattr(base, "fitted_shares", None))
        check.test("R_backend_method_exists", not calls or method_exists,
                   dict(wrapper_calls="res_py$fitted_shares()", python_method_exists=method_exists,
                        source="R/R/quaidsce.R"), "R fitted()/residuals() use an available backend API")
    run("R-static-fitted_backend", "interface_static", r_fitted,
        "source contract: R fitted values must not silently become NULL")
    def r_se(check):
        reconstructs = "se <- sqrt(diag(V))" in r_methods
        raw = np.sqrt(np.maximum(np.diag(base.V), 0))
        check.test("R_summary_preserves_undefined_SE", not reconstructs or np.array_equal(np.isnan(raw), np.isnan(base.analytic_se)),
                   dict(raw_zero_elasticity_SE=int(np.sum(raw[-36:] == 0)),
                        python_undefined_elasticity_SE=int(np.isnan(base.analytic_se[-36:]).sum()),
                        source="R/R/methods.R"), "undefined analytical elasticity SEs remain NA/NaN")
    run("R-static-analytic_SE", "interface_static", r_se,
        "source contract: R summary does not convert undefined elasticity uncertainty into zero")
    def r_defaults(check):
        check.test("R_default_method_ifgnls", 'method = "ifgnls"' in r_source)
        check.test("R_forwards_ivexp", "ivexp = if (!is.null(ivexp))" in r_source)
        check.test("R_preserves_empty_selection", "selection_prices = if (!is.null(selection_prices)) as.list(selection_prices) else NULL" in r_source)
    run("R-static-forwarding", "interface_static", r_defaults)

    class FakeSFI:
        def __init__(self, data):
            self.stored = {}
            stored = self.stored
            self.Data = types.SimpleNamespace(get=lambda var: data[var].tolist())
            self.Scalar = types.SimpleNamespace(setValue=lambda name, value: stored.__setitem__(name, value))
            self.Macro = types.SimpleNamespace(setLocal=lambda name, value: stored.__setitem__(name, value))
            self.Matrix = types.SimpleNamespace(
                create=lambda name, rows, cols, value: stored.__setitem__(name, np.full((rows, cols), value, float)),
                storeAt=lambda name, row, col, value: stored[name].__setitem__((row, col), value),
                setColNames=lambda name, names: stored.__setitem__(name+"_cols", list(names)),
                setRowNames=lambda name, names: stored.__setitem__(name+"_rows", list(names)))

    def bridge(check, custom=False):
        from pyquaidsce.stata_bridge import run_from_stata, _BOOT_STATE
        case_id = "cf-iv_one" if custom else BASE_ID
        frame, args = frames.get(case_id, frames[BASE_ID])
        reference = cache.get(case_id, base)
        data = frame.copy()
        data["_touse"] = 1
        fake = FakeSFI(data)
        with patch.dict(sys.modules, {"sfi": fake}):
            run_from_stata(" ".join(SHARES), " ".join(PRICES), "total", "x1 x2", 10,
                           ivexp_str="iv1" if custom else "", verbose=False)
        check.close("bridge_real_estimator_b", fake.stored["__pyq_b"][0], reference.b)
        check.close("bridge_real_estimator_V", fake.stored["__pyq_V"], reference.V)
        check.close("bridge_real_elasticities", fake.stored["__pyq_elas_u"], reference.elas.uncompensated)
        check.test("bridge_parameter_labels", fake.stored["__pyq_b_cols"] == reference.names)
        check.test("bridge_sample_scalar", fake.stored["r_nobs"] == reference.nobs)
        if custom:
            check.close("bridge_reduced_form_b", fake.stored["__pyq_rf_b"][0], reference.reduced_form.b)
        _BOOT_STATE.clear()
    run("Stata-bridge-real_base", "interface_python", lambda c: bridge(c))
    if "cf-iv_one" in cache:
        run("Stata-bridge-real_IV", "interface_python", lambda c: bridge(c, True))

    def async_bridge(check):
        from pyquaidsce.stata_bridge import (launch_from_stata, poll_bootstrap, load_stata_results,
                                             _BOOT_STATE, kill_bootstrap)
        data = original.copy()
        data["_touse"] = 1
        fake = FakeSFI(data)
        # A source checkout is importable by the standalone subprocess even if
        # the package has not been installed into site-packages.
        path = str(ROOT / "src") + os.pathsep + os.environ.get("PYTHONPATH", "")
        with patch.dict(sys.modules, {"sfi": fake}), patch.dict(os.environ, {"PYTHONPATH": path}):
            try:
                launch_from_stata(" ".join(SHARES), " ".join(PRICES), "total", "x1 x2", 10,
                                  reps=8, seed=16001, n_jobs=2, mp_context="spawn", verbose=False)
                deadline = time.perf_counter() + 60
                while time.perf_counter() < deadline:
                    poll_bootstrap()
                    if fake.stored.get("_pyq_boot_done") or fake.stored.get("_pyq_boot_err"):
                        break
                    time.sleep(.05)
                check.test("async_subprocess_complete", fake.stored.get("_pyq_boot_done") == 1,
                           fake.stored.get("_pyq_boot_errmsg"), "completed result")
                if fake.stored.get("_pyq_boot_done"):
                    load_stata_results("b", "V", "ei", "eu", "ec", "theta", "sigma")
                    reference = cache.get("boot-serial")
                    if reference is not None:
                        check.close("async_real_bootstrap_b", fake.stored["b"][0], reference.b)
                        check.close("async_real_bootstrap_V", fake.stored["V"], reference.V)
                    check.close("async_theta", fake.stored["theta"][0], base.theta)
                    check.test("async_bootstrap_counts", fake.stored["r_boot_reps_ok"] == 8)
            finally:
                kill_bootstrap()
                _BOOT_STATE.clear()
    run("Stata-bridge-real_async_bootstrap", "interface_python", async_bridge,
        "real Python subprocess and bootstrap with simulated sfi transport")

    for name in ("R", "Stata"):
        if only and name + "-native" not in only:
            continue
        row = dict(id=name + "-native", group="interface_native", status="NOT_RUN", elapsed_seconds=0,
                   checks=[], reason="native execution is supplemental and no native runtime was used")
        write_json(out / "cases" / f"{name}-native.json", row)
        records.append(row)
    return records
