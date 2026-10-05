"""Run the complete small4 audit without editing the estimator.

    python -m validation.run --output validation/results --large-reps 100

A nonzero exit status means a failed assertion, unexplained nonconvergence,
or an audit error. Missing optional native interfaces are recorded explicitly.
"""
import argparse
import contextlib
import csv
import dataclasses
import hashlib
import inspect
import io
import json
import multiprocessing as mp
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time
import traceback
import unittest
from unittest.mock import patch
import warnings

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))
import numpy as np
import pandas as pd
import scipy
from threadpoolctl import threadpool_limits, threadpool_info
import pyquaidsce
from pyquaidsce import quaidsce
from .cases import BASE, BASE_ID, SHARES, PRICES, scenarios
from .oracles import (Checks, TOLERANCES, check_fit, reduced_form_checks,
                      observation_derivative_checks, model_point_elasticity_checks,
                      independent_probit, bootstrap_checks, stationary_measure)


def plain(x):
    if isinstance(x, np.ndarray):
        return plain(x.tolist())
    if isinstance(x, np.generic):
        return plain(x.item())
    if isinstance(x, dict):
        return {str(k): plain(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [plain(v) for v in x]
    if isinstance(x, float) and not np.isfinite(x):
        return "NaN" if np.isnan(x) else ("Infinity" if x > 0 else "-Infinity")
    if isinstance(x, Path):
        return str(x)
    return x


def write_json(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(plain(obj), ensure_ascii=False, indent=2,
                               allow_nan=False) + "\n", encoding="utf-8")


def write_csv(path, rows, fields=None):
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields, lineterminator="\n")
        w.writeheader()
        for row in rows:
            converted = {}
            for key in fields:
                value = plain(row.get(key))
                converted[key] = json.dumps(value, ensure_ascii=False) if isinstance(value, (dict, list)) else value
            w.writerow(converted)


def fixture(original, name):
    df, change = original.copy(), {}
    if name in ("price_logs", "expenditure_logs", "both_logs", "derived_logs"):
        if name in ("price_logs", "both_logs", "derived_logs"):
            names = []
            for p in PRICES:
                col = "l" + p
                df[col] = np.log(df[p].to_numpy(float))
                names.append(col)
            change.update(prices=None, lnprices=names)
        if name in ("expenditure_logs", "both_logs", "derived_logs"):
            df["lm"] = np.log(df["total"].to_numpy(float))
            change.update(expenditure=None, lnexpenditure="lm")
    if name.startswith("derived") or name in ("cf_nan",):
        lp = np.log(df[PRICES].to_numpy(float))
        lm = np.log(df["total"].to_numpy(float))
        x1, x2 = df["x1"].to_numpy(float), df["x2"].to_numpy(float)
        df["iv1"] = x1**2 + 0.35 * x2**2
        df["iv2"] = np.random.default_rng(2601005).normal(size=len(df))
        X = np.column_stack([lp, df[["x1", "x2"]].to_numpy(float), df[["iv1"]].to_numpy(float), np.ones(len(df))])
        Q, R = np.linalg.qr(X, mode="reduced")
        df["cf"] = lm - X @ np.linalg.solve(R, Q.T @ lm)
        df["cf2"] = 0.7 * df["cf"] + np.random.default_rng(678).normal(0, .03, len(df))
        df["cf_scaled"] = 2 * df["cf"]
        df["constant_cf"] = 0.0
        df["s1"] = (x1 - x1.mean())**2 / max(x1.std()**2, 1e-12)
        df["s2"] = np.random.default_rng(101).normal(size=len(df))
        df["x1_copy"] = x1
        df["lm"] = lm
        if name == "cf_nan":
            df.loc[0, "cf"] = np.nan
    if name == "renamed":
        mapping = {col: "renamed_" + col for col in df.columns}
        df = df.rename(columns=mapping)
        change = dict(shares=[mapping[x] for x in SHARES], prices=[mapping[x] for x in PRICES],
                      expenditure=mapping["total"], demographics=[mapping[x] for x in BASE["demographics"]])
    elif name == "irrelevant":
        df["unused"] = np.random.default_rng(1).normal(size=len(df))
    elif name == "missing_irrelevant":
        df["unused"] = np.nan
    elif name == "reverse_rows":
        df = df.iloc[::-1].reset_index(drop=True)
    elif name == "custom_index":
        df.index = [f"household_{i * 7}" for i in range(len(df))]
    elif name == "permuted_goods":
        change = dict(shares=list(reversed(SHARES)), prices=list(reversed(PRICES)))
    elif name == "three_goods":
        selected = SHARES[:3]
        subtotal = df[selected].sum(axis=1).to_numpy(float)
        keep = subtotal > 0
        df = df.loc[keep].copy()
        subtotal = subtotal[keep]
        df[selected] = df[selected].to_numpy(float) / subtotal[:, None]
        df["total"] = df["total"].to_numpy(float) * subtotal
        change.update(shares=selected, prices=PRICES[:3])
    elif name in ("missing_required", "infinite_required"):
        df.loc[[1, 11, 31], "p1"] = np.nan if name == "missing_required" else np.inf
        df.loc[51, "x1"] = np.nan
    elif name == "negative_share":
        df.loc[0, "sw1"] = -.1
    elif name in ("zero_price", "negative_price"):
        df.loc[0, "p1"] = 0 if name == "zero_price" else -1
    elif name == "zero_expenditure":
        df.loc[0, "total"] = 0
    elif name == "empty_sample":
        df["p1"] = np.nan
    elif name == "no_zero_share":
        df["sw1"] = df["sw1"].clip(lower=.01)
    elif name == "all_zero_share":
        df["sw1"] = 0.0
    elif name == "bad_sum":
        df["sw1"] += .01
    return df, change


def compare_fit(check, res, ref, prefix="equivalent", covariance=True):
    for name in ("theta", "b", "sigma", "llf"):
        check.close(prefix + "_" + name, getattr(res, name), getattr(ref, name), atol=1e-6, rtol=1e-5)
    if covariance:
        check.close(prefix + "_V", res.V, ref.V, atol=1e-6, rtol=1e-5)


def matrix_checks(check, res, ref, perm, demos=False):
    c, b = res.coefs, ref.coefs
    for name in ("alpha", "beta", "lam", "delta"):
        check.close("permutation_" + name, getattr(c, name), getattr(b, name) if demos else getattr(b, name)[perm], atol=1e-6, rtol=1e-5)
    check.close("permutation_gamma", c.gamma, b.gamma if demos else b.gamma[np.ix_(perm, perm)], atol=1e-6, rtol=1e-5)
    check.close("permutation_eta", c.eta, b.eta[perm] if demos else b.eta[:, perm], atol=1e-6, rtol=1e-5)
    check.close("permutation_rho", c.rho, b.rho[perm] if demos else b.rho, atol=1e-6, rtol=1e-5)
    check.close("permutation_income", res.elas.income, ref.elas.income if demos else ref.elas.income[perm], atol=1e-6, rtol=1e-5)
    check.close("permutation_price", res.elas.uncompensated,
                ref.elas.uncompensated if demos else ref.elas.uncompensated[np.ix_(perm, perm)], atol=1e-6, rtol=1e-5)


def capture_fit(out, case, original, cache, frames, audit_timeout=90):
    df, changes = fixture(original, case.fixture)
    kw = dict(BASE, **changes)
    kw.update(case.kwargs)
    # References always belong to a previous frozen scenario.
    ref = cache.get(case.reference)
    if case.reference and ref is None:
        return dict(id=case.id, group=case.group, status="BLOCKED", reason="reference fit unavailable: " + case.reference,
                    checks=[], elapsed_seconds=0.0), None
    if "initial" in case.checks or "initial_sigma" in case.checks:
        kw["initial"] = ref.theta.copy()
    if "initial_sigma" in case.checks or "sigma_only" in case.checks:
        kw["sigma_initial"] = ref.sigma.copy()
    callback = []
    if "callback" in case.checks:
        kw["log"] = lambda *values: callback.append(" ".join(map(str, values)))
    log = io.StringIO()
    check = Checks()
    started = time.perf_counter()
    before = [(x["prefix"], x["num_threads"]) for x in threadpool_info()]
    record = dict(id=case.id, group=case.group, expectation=case.expectation,
                  fixture=case.fixture, kwargs={k: v for k, v in kw.items() if k != "log"},
                  reference=case.reference, checks=[], warnings=[], audit_timeout_seconds=audit_timeout)
    res = None
    reduced_count = [0]
    from pyquaidsce import estimator, bootstrap as bootstrap_module
    original_rf = estimator.fit_expenditure_reduced_form
    original_rep = bootstrap_module._one_rep
    def reduced(*args, **kwargs):
        reduced_count[0] += 1
        return original_rf(*args, **kwargs)
    def rep(task):
        if case.id == "lifecycle-all_failure" or (case.id == "lifecycle-partial_failure" and task[0] == 2):
            return task[0], None, "InjectedValidationFailure: deterministic lifecycle test"
        return original_rep(task)
    with (contextlib.redirect_stdout(log), contextlib.redirect_stderr(log),
          warnings.catch_warnings(record=True) as warns, contextlib.ExitStack() as stack):
        warnings.simplefilter("always")
        if "count_reduced_forms" in case.checks:
            stack.enter_context(patch.object(estimator, "fit_expenditure_reduced_form", side_effect=reduced))
        if case.id in ("lifecycle-partial_failure", "lifecycle-all_failure"):
            stack.enter_context(patch.object(bootstrap_module, "_one_rep", side_effect=rep))
        try:
            # Protect exploratory invalid cases from indefinite iterations;
            # this is a cooperative audit bound, not an alteration of defaults.
            res = quaidsce(df, **kw, _deadline=time.perf_counter() + audit_timeout)
            record["fit_seconds"] = time.perf_counter() - started
            if case.exception:
                check.test("invalid_input_rejected", False, "returned a result", list(case.exception))
            else:
                d, design, index, mask, metrics = check_fit(check, df, kw, res, "nonconvergence" in case.checks)
                record.update(metrics)
                observation_derivative_checks(check, res, d, index)
                model_point_elasticity_checks(check, res)
                if res.reduced_form is not None:
                    reduced_form_checks(check, d, df, kw, res, mask)
                if "equivalent" in case.checks or "equivalent_point" in case.checks:
                    compare_fit(check, res, ref, covariance="equivalent_point" not in case.checks)
                if "initial" in case.checks or "initial_sigma" in case.checks or "sigma_only" in case.checks:
                    compare_fit(check, res, cache[BASE_ID], prefix="warm_start", covariance=True)
                if "complete_case" in case.checks:
                    manual = quaidsce(df.iloc[np.flatnonzero(mask)].copy(), **kw)
                    compare_fit(check, res, manual, prefix="manual_complete_case")
                if "permuted_goods" in case.checks:
                    matrix_checks(check, res, ref, list(reversed(range(res.spec.neqn))))
                if "permuted_demos" in case.checks:
                    matrix_checks(check, res, ref, list(reversed(range(res.spec.ndemo))), demos=True)
                if "selection_reorder" in case.checks:
                    check.close("selection_reorder_structural", res.b[:res.spec.n_full], ref.b[:ref.spec.n_full], atol=1e-6, rtol=1e-5)
                    check.close("selection_reorder_elasticities", res.elas.as_stata_vector(), ref.elas.as_stata_vector(), atol=1e-6, rtol=1e-5)
                if "dropped_column" in case.checks:
                    check.test("collinear_column_reported", all(len(p.dropped) == 1 for p in res.probits))
                if "residual_rescaled" in case.checks:
                    check.close("rescaled_cf_coefficients", 2 * res.coefs.cfcoef, ref.coefs.cfcoef, atol=1e-6, rtol=1e-5)
                    check.close("rescaled_elasticities", res.elas.as_stata_vector(), ref.elas.as_stata_vector(), atol=1e-6, rtol=1e-5)
                if "stationarity" in case.checks:
                    weight = np.eye(res.spec.n_eq_estimated) if res.method == "nls" else None
                    if res.method == "fgnls":
                        preliminary = dict(kw, method="nls")
                        weight = quaidsce(df, **preliminary).sigma
                    measure = stationary_measure(res, d, sigma=weight)
                    check.test("independent_stationarity", measure <= TOLERANCES["stationarity"], measure, TOLERANCES["stationarity"])
                    record["stationarity"] = measure
                if "covariance" in case.checks or case.id in (BASE_ID, "core-cq-nls-gn-zero", "core-cq-fgnls-gn-zero"):
                    from .oracles import independent_fitted
                    from pyquaidsce import jacobian_full, delta_matrix
                    J = np.einsum("tmf,fk->tmk", jacobian_full(res.theta, d, res.spec), delta_matrix(res.spec))
                    sigma = res.sigma
                    if res.method == "fgnls" and kw.get("vce_sigma", "objective") == "objective":
                        sigma = cache["core-cq-nls-gn-zero"].sigma
                    if res.method == "ifgnls" and kw.get("vce_sigma", "objective") == "objective":
                        # At a fixed point, objective/final covariance differ only within sigma_tol.
                        tolerance = 1e-4
                    else:
                        tolerance = 1e-5
                    G = np.einsum("tik,ij,tjl->kl", J, np.linalg.inv(sigma), J)
                    check.close("independent_covariance_formula", res.V_est, np.linalg.inv(G), atol=1e-7, rtol=tolerance)
                if case.id == BASE_ID:
                    independent_probit(check, res, design, d.shares)
                if "bootstrap" in case.checks:
                    bootstrap_checks(check, res)
                    if "partial_failure" in case.checks:
                        check.test("partial_failure_accounting", res.boot.reps_ok == 3 and len(res.boot.failures) == 1)
                        check.test("partial_failure_identity", "rep 2:" in res.boot.failures[0])
                    else:
                        check.test("no_unexplained_bootstrap_failures", res.boot.reps_ok == res.boot.reps_requested,
                                   res.boot.failures, [])
                    if "count_reduced_forms" in case.checks:
                        check.test("reduced_form_rebuilt_each_draw", reduced_count[0] == 1 + res.boot.reps_requested,
                                   reduced_count[0], 1 + res.boot.reps_requested)
                    if "bootstrap_compare" in case.checks:
                        exact = case.id in ("boot-repeat", "boot-parallel", "boot-spawn", "boot-forkserver", "boot-blas_none")
                        if exact:
                            check.test("bootstrap_identical_seed_draws", np.array_equal(res.boot.b_star, ref.boot.b_star))
                        else:
                            check.close("bootstrap_equivalent_seed_draws", res.boot.b_star, ref.boot.b_star, atol=1e-6, rtol=1e-5)
                        check.test("bootstrap_same_failures", res.boot.failures == ref.boot.failures)
                if "verbose" in case.checks:
                    check.test("verbose_progress", "Estimating first-stage" in log.getvalue())
                if "gn_verbose" in case.checks:
                    check.test("gn_iteration_progress", "GN " in log.getvalue())
                if "callback" in case.checks:
                    check.test("callback_received_progress", bool(callback))
                    record["callback_messages"] = callback
                if kw.get("verbose") is False and not kw.get("gn_verbose") and not kw.get("log"):
                    check.test("quiet_stdout", log.getvalue() == "", log.getvalue()[:500], "")
                check.test("public_named_roundtrip", len(res.named()) == len(res.names)
                           and all(res.get(name) == value for name, value in zip(res.names, res.b)))
                check.test("summary_contains_sample", "Number of obs" in res.summary())
                check.test("elasticity_table_labels", "Marshallian" in res.elasticity_tables() and "Hicksian" in res.elasticity_tables())
        except Exception as exc:
            record["exception"] = dict(type=type(exc).__name__, message=str(exc), traceback=traceback.format_exc())
            if case.exception:
                check.test("expected_exception_type", type(exc).__name__ in case.exception, type(exc).__name__, list(case.exception))
                if case.message:
                    check.test("expected_exception_message", case.message.lower() in str(exc).lower(), str(exc), case.message)
            else:
                check.test("unexpected_exception", False, f"{type(exc).__name__}: {exc}", "no exception")
        record["warnings"] = [dict(type=type(w.message).__name__, message=str(w.message)) for w in warns]
    after = [(x["prefix"], x["num_threads"]) for x in threadpool_info()]
    check.test("BLAS_state_restored", before == after, after, before)
    record["elapsed_seconds"] = time.perf_counter() - started
    record["checks"] = check.rows
    good = all(x["passed"] for x in check.rows)
    record["status"] = ("EXPECTED_ERROR" if case.exception else "PASS") if good else "FAIL"
    if res is not None:
        record["result"] = dict(converged=res.converged, nobs=res.nobs, method=res.method,
            n_free=res.spec.n_free, n_reported=len(res.b), n_gn=res.n_gn, n_outer=res.n_outer,
            llf=res.llf, theta=res.theta, b=res.b, names=res.names, sigma=res.sigma,
            V_est=res.V_est, se=res.se, analytic_se=res.analytic_se, notes=res.notes,
            elas=dict(income=res.elas.income, uncompensated=res.elas.uncompensated,
                      compensated=res.elas.compensated),
            probit_dropped=[p.dropped for p in res.probits])
        if not res.converged and "nonconvergence" not in case.checks and not case.exception:
            record["status"] = "NONCONVERGED"
        rows = [dict(name=name, estimate=value, se=res.se[i], analytic_se=res.analytic_se[i])
                for i, (name, value) in enumerate(zip(res.names, res.b))]
        write_csv(out / "coefficients" / f"{case.id}.csv", rows)
        if res.boot is not None:
            record["bootstrap"] = dict(requested=res.boot.reps_requested, successful=res.boot.reps_ok,
                                       failures=res.boot.failures)
            write_csv(out / "bootstrap" / f"{case.id}-draws.csv",
                      [dict(draw_number=i + 1, **dict(zip(res.names, b))) for i, b in enumerate(res.boot.b_star)])
            write_csv(out / "bootstrap" / f"{case.id}-covariance.csv",
                      [dict(name=name, **dict(zip(res.names, row))) for name, row in zip(res.names, res.V)])
        if case.id == BASE_ID:
            write_csv(out / "reference-covariance.csv",
                      [dict(name=name, **dict(zip(res.names, row))) for name, row in zip(res.names, res.V)])
        # Save successful calculation objects even if an assertion failed: later
        # cases must not silently change their reference to a more favorable fit.
        cache[case.id] = res
        frames[case.id] = (df, kw)
    if log.getvalue():
        (out / "logs").mkdir(parents=True, exist_ok=True)
        (out / "logs" / f"{case.id}.log").write_text(log.getvalue(), encoding="utf-8")
    write_json(out / "cases" / f"{case.id}.json", record)
    return record, res


def existing_tests(out):
    buffer = io.StringIO()
    rows = []
    class Result(unittest.TextTestResult):
        def startTest(self, test):
            self.started = time.perf_counter()
            super().startTest(test)
        def addSuccess(self, test):
            rows.append(dict(id=test.id(), status="PASS", elapsed_seconds=time.perf_counter()-self.started))
            super().addSuccess(test)
        def addFailure(self, test, err):
            rows.append(dict(id=test.id(), status="FAIL", elapsed_seconds=time.perf_counter()-self.started))
            super().addFailure(test, err)
        def addError(self, test, err):
            rows.append(dict(id=test.id(), status="ERROR", elapsed_seconds=time.perf_counter()-self.started))
            super().addError(test, err)
        def addSkip(self, test, reason):
            rows.append(dict(id=test.id(), status="SKIP", reason=reason))
            super().addSkip(test, reason)
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"), top_level_dir=str(ROOT))
    with contextlib.redirect_stdout(buffer), contextlib.redirect_stderr(buffer):
        result = unittest.TextTestRunner(stream=buffer, verbosity=2, resultclass=Result).run(suite)
    (out / "existing-tests.log").write_text(buffer.getvalue(), encoding="utf-8")
    write_csv(out / "existing-tests.csv", rows)
    return dict(tests_run=result.testsRun, failures=len(result.failures), errors=len(result.errors),
                skipped=len(result.skipped), successful=result.wasSuccessful())


def metadata(large_reps):
    source = ROOT / "src" / "pyquaidsce"
    return dict(package_version=pyquaidsce.__version__,
        source_base_commit="93ea81e4293ea4d2f0db899588079e1acc22c745",
        checkout_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        branch=subprocess.check_output(["git", "branch", "--show-current"], cwd=ROOT, text=True).strip(),
        created_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        python=platform.python_version(), platform=platform.platform(),
        dependencies=dict(numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__),
        data_sha256=hashlib.sha256((ROOT / "bench" / "small4.dta").read_bytes()).hexdigest(),
        source_sha256={str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted(source.glob("*.py"))},
        blas=threadpool_info(), large_bootstrap_reps=large_reps,
        available_mp_contexts=mp.get_all_start_methods(),
        native_interfaces={k: shutil.which(k) for k in ("Rscript", "stata", "stata-se", "stata-mp")},
        audit_blas_caller_threads=1,
        tolerances=TOLERANCES)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "validation" / "results")
    parser.add_argument("--large-reps", type=int, default=100)
    parser.add_argument("--groups", nargs="*", help="selected groups; reference cases are included automatically")
    parser.add_argument("--only", nargs="*", help="selected case IDs; reference cases are included automatically")
    parser.add_argument("--contract-only", nargs="+", help="selected public/helper/interface contract IDs")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument("--skip-contracts", action="store_true")
    args = parser.parse_args()
    if args.large_reps < 2:
        parser.error("--large-reps must be at least 2")
    if args.contract_only and args.skip_contracts:
        parser.error("--contract-only cannot be combined with --skip-contracts")
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    design = scenarios(args.large_reps)
    wanted = {x.id for x in design if (not args.groups or x.group in args.groups)
              and (not args.only or x.id in args.only)}
    if args.contract_only:
        wanted = set()
        if any(x.startswith("output-") or x == "Stata-bridge-real_async_bootstrap" or x == "boot-positive_timeout" for x in args.contract_only):
            wanted.add("boot-serial")
        if "Stata-bridge-real_IV" in args.contract_only:
            wanted.add("cf-iv_one")
        if "public-latent_shares" in args.contract_only:
            wanted.add("cf-both_distinct")
    if not wanted and not args.contract_only:
        parser.error("no scenarios selected")
    byid = {x.id: x for x in design}
    while True:
        extra = {byid[x].reference for x in wanted if byid[x].reference}
        if extra.issubset(wanted):
            break
        wanted |= extra
    # Bootstrap/contract and warm-start checks also need the explicit baseline.
    wanted |= {BASE_ID, "core-cq-nls-gn-zero", "core-cq-fgnls-gn-zero"}
    selected = [x for x in design if x.id in wanted]
    write_json(out / "manifest.json", [dataclasses.asdict(x) for x in design])
    write_json(out / "environment.json", metadata(args.large_reps))
    original = pd.read_stata(ROOT / "bench" / "small4.dta")
    cache, frames, records = {}, {}, []
    with threadpool_limits(limits=1, user_api="blas"):
        existing = {} if args.skip_existing else existing_tests(out)
        print("Existing suite:", existing, flush=True)
        for number, case in enumerate(selected, 1):
            if case.kwargs.get("mp_context") == "forkserver" and "forkserver" not in mp.get_all_start_methods():
                row = dict(id=case.id, group=case.group, status="NOT_RUN", reason="forkserver unsupported", checks=[], elapsed_seconds=0)
            else:
                row, _ = capture_fit(out, case, original, cache, frames)
            records.append(row)
            bad = [x["name"] for x in row.get("checks", []) if not x["passed"]]
            print(f"[{number}/{len(selected)}] {case.id}: {row['status']} ({row['elapsed_seconds']:.2f}s) {', '.join(bad[:5])}", flush=True)
            if number % 10 == 0 or number == len(selected):
                write_json(out / "results.json", records)
        if not args.skip_contracts:
            from .contracts import contract_tests
            contract_records = contract_tests(out, original, cache, frames, only=args.contract_only)
            records.extend(contract_records)
            if args.contract_only:
                missing = set(args.contract_only) - {r["id"] for r in contract_records}
                if missing:
                    parser.error("unknown or unavailable contract IDs: " + ", ".join(sorted(missing)))
        env_path = out / "environment.json"
        environment = json.loads(env_path.read_text())
        environment["audit_source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                               for p in sorted((ROOT / "validation").glob("*.py"))}
        environment["completed_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        write_json(env_path, environment)
        from .report import generate
        generate(out, records, existing, design, complete=(len(selected) == len(design)))
    failed = any(x["status"] in ("FAIL", "NONCONVERGED", "BLOCKED", "ERROR") for x in records)
    return int(failed or (bool(existing) and not existing["successful"]))


if __name__ == "__main__":
    sys.exit(main())
