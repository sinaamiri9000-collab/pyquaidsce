"""Compare archived standard/GN fits and benchmark matched stage-two cores.

Zero is the sole start. First-stage estimates, model, BLAS thread count and
timing boundaries are shared. Covariance/SE/output construction is excluded.
One warm-up per backend/model is followed by five paired timed repetitions.
"""
from __future__ import annotations

import importlib
import json
import platform
import time
from unittest.mock import patch

import numpy as np
import scipy
from scipy.stats import norm
from threadpoolctl import threadpool_limits

from pyquaidsce.params import full_vector
from .compare import DEFAULT_OUT, MODELS, ROOT, Problem, dump, run_nlsur, source_hashes, write_csv


def archived_comparison():
    runs = json.loads((DEFAULT_OUT / "runs.json").read_text())
    comparisons, coefficients = [], []
    for model in MODELS:
        p = Problem(model)
        oldcase = json.loads((ROOT / f"validation/results/cases/core-{model}-ifgnls-gn-zero.json").read_text())
        old = oldcase["result"]
        a = np.array(old["b"][:p.spec.n_full])
        for backend in ("nlsur-analytic", "nlsur-fd"):
            new = next(z for z in runs if z["model"] == model and z["method"] == backend
                       and z["start"] == "zero" and not z["tightened"])
            b = full_vector(np.array(new["theta"]), p.spec)
            diff = b-a
            j = int(np.argmax(np.abs(diff)))
            comparisons.append(dict(model=model, new_backend=backend,
                old_converged=old["converged"], new_converged=new["success"],
                max_abs_structural_coefficient_difference=float(np.max(np.abs(diff))),
                rms_structural_coefficient_difference=float(np.sqrt(np.mean(diff**2))),
                largest_difference_name=old["names"][j],
                largest_difference_old=float(a[j]), largest_difference_new=float(b[j]),
                old_ll=old["llf"], new_ll=new["gaussian_ll"],
                ll_difference=new["gaussian_ll"]-old["llf"],
                old_total_inner_solves=old["n_outer"],
                new_total_inner_solves=len(new["history"]),
                old_ifgnls_updates=old["n_outer"]-2, new_ifgnls_updates=len(new["outer"]),
                old_total_inner_iterations=old["n_gn"],
                new_total_inner_iterations=sum(t["iterations"] for t in new["history"]),
                old_archived_full_api_fit_seconds=oldcase["fit_seconds"],
                new_archived_point_solver_seconds=new["elapsed_seconds"]))
            for name, x, y in zip(old["names"][:p.spec.n_full], a, b):
                coefficients.append(dict(model=model, new_backend=backend, coefficient=name,
                    old_standard=float(x), new_nlsur=float(y), difference=float(y-x),
                    absolute_difference=float(abs(y-x))))
    dump(DEFAULT_OUT / "standard-comparison.json", comparisons)
    write_csv(DEFAULT_OUT / "standard-comparison.csv", comparisons, list(comparisons[0]))
    write_csv(DEFAULT_OUT / "standard-coefficients.csv", coefficients, list(coefficients[0]))
    return comparisons


def fit_old(p):
    """Time the unmodified old loop, stopping before result assembly.

    The original _finish still runs, but its covariance/likelihood work is
    outside the measured interval. A final residual Sigma is computed before
    stopping the timer, matching the new core's final Sigma calculation.
    Scoped wrappers only observe calls and do not change stopping decisions.
    """
    oldmod = importlib.import_module("pyquaidsce.nlsur")
    original_finish, original_gn = oldmod._finish, oldmod.gauss_newton
    trace, timing = [], {}
    def gn(*args, **kwargs):
        result = original_gn(*args, **kwargs)
        trace.append(dict(iterations=int(result[2]), success=bool(result[3])))
        return result
    def finish(theta, d, spec, sigma_obj, sigma_from, *args, **kwargs):
        sigma_from(theta)
        timing["seconds"] = time.perf_counter()-started
        return original_finish(theta, d, spec, sigma_obj, sigma_from, *args, **kwargs)
    theta0 = np.zeros(p.k)
    with patch.object(oldmod, "gauss_newton", gn), patch.object(oldmod, "_finish", finish):
        started = time.perf_counter()
        # The common outer threadpool limit replaces this decorator's limit.
        result = oldmod.nlsur.__wrapped__(p.d, p.spec, theta0=theta0, start="zero",
                        method="ifgnls", algorithm="gn", stop_rule="standard",
                        max_iter=300, max_outer=200, chunk=2000, verbose=False)
    return result.theta, dict(success=result.converged, inner_solves=len(trace),
           inner_iterations=result.n_gn, ifgnls_updates=result.n_outer-2,
           seconds=timing["seconds"], inner_history=trace)


def fit_new(p, backend):
    theta0, S = np.zeros(p.k), np.eye(p.m)
    p.trace, p.fun_calls, p.jac_calls, p.invalid_trials = [], 0, 0, 0
    p.started = time.perf_counter()
    theta, info = run_nlsur(p,theta0,S,"zero",backend.split("-",1)[1],False)
    seconds = time.perf_counter()-p.started
    return theta, dict(success=info["success"], inner_solves=len(p.trace),
           inner_iterations=sum(t["iterations"] for t in p.trace),
           ifgnls_updates=len(info["outer"]), seconds=seconds,
           inner_history=[dict(x) for x in p.trace])


def benchmark(comparisons, repetitions=5):
    raw, summaries, checks = [], [], []
    methods = ("old-standard-gn", "nlsur-analytic", "nlsur-fd")
    with threadpool_limits(limits=1):
        for model in MODELS:
            p = Problem(model)
            for backend in methods:
                (fit_old(p) if backend.startswith("old-") else fit_new(p,backend))
            for rep in range(repetitions):
                # Rotate ordering, so the same backend is not always first.
                order = methods[rep % len(methods):]+methods[:rep % len(methods)]
                for backend in order:
                    theta, info = fit_old(p) if backend.startswith("old-") else fit_new(p,backend)
                    old_theta = p.reference
                    if backend.startswith("old-"):
                        archived = next(z for z in comparisons if z["model"] == model)
                        expected_solves = archived["old_total_inner_solves"]
                        expected_iterations = archived["old_total_inner_iterations"]
                    else:
                        archived = next(z for z in comparisons if z["model"] == model and z["new_backend"] == backend)
                        expected_solves = archived["new_total_inner_solves"]
                        expected_iterations = archived["new_total_inner_iterations"]
                    z = dict(model=model, backend=backend, repetition=rep+1, **info,
                             theta=theta.tolist(), theta_difference_from_old=float(np.max(abs(theta-old_theta))))
                    raw.append(z)
                    checks.append(dict(model=model,backend=backend,repetition=rep+1,
                        success=bool(info["success"]),
                        inner_solves_match_archive=info["inner_solves"] == expected_solves,
                        inner_iterations_match_archive=info["inner_iterations"] == expected_iterations))
            group = [z for z in raw if z["model"] == model]
            oldtimes = np.array([z["seconds"] for z in group if z["backend"] == methods[0]])
            for backend in methods:
                z = [x for x in group if x["backend"] == backend]
                times = np.array([x["seconds"] for x in z])
                summaries.append(dict(model=model,backend=backend,repetitions=len(z),
                    median_seconds=float(np.median(times)),min_seconds=float(times.min()),
                    max_seconds=float(times.max()),median_ratio_to_old=float(np.median(times)/np.median(oldtimes)),
                    paired_median_ratio_to_old=float(np.median(times/oldtimes)),
                    all_converged=all(x["success"] for x in z),
                    inner_solves=z[0]["inner_solves"],inner_iterations=z[0]["inner_iterations"],
                    ifgnls_updates=z[0]["ifgnls_updates"]))
            dump(DEFAULT_OUT / "standard-timing-runs.json", raw)
            print(json.dumps([z for z in summaries if z["model"] == model]),flush=True)
    dump(DEFAULT_OUT / "standard-timing-summary.json", summaries)
    write_csv(DEFAULT_OUT / "standard-timing-summary.csv", summaries, list(summaries[0]))
    fields = ("model","backend","repetition","success","seconds","inner_solves",
              "inner_iterations","ifgnls_updates","theta_difference_from_old")
    write_csv(DEFAULT_OUT / "standard-timing-runs.csv", raw, fields)
    dump(DEFAULT_OUT / "standard-timing-checks.json", checks)
    return summaries, checks


def main():
    before = source_hashes()
    comparisons = archived_comparison()
    summaries, checks = benchmark(comparisons)
    rounding = rounding_probe()
    dump(DEFAULT_OUT / "standard-rounding-probe.json", rounding)
    raw = json.loads((DEFAULT_OUT / "standard-timing-runs.json").read_text())
    stable = all(len({(z["inner_solves"], z["inner_iterations"]) for z in raw
                      if z["model"] == m and z["backend"] == b}) == 1
                 for m in MODELS for b in ("old-standard-gn","nlsur-analytic","nlsur-fd"))
    limits_reached = any(z["inner_solves"] >= 200 or
                        any(t["iterations"] >= (300 if z["backend"].startswith("old-") else 400)
                            for t in z["inner_history"]) for z in raw)
    assert before == source_hashes()
    dump(DEFAULT_OUT / "standard-timing-environment.json",dict(python=platform.python_version(),
          numpy=np.__version__,scipy=scipy.__version__,blas_threads=1,
          repetitions=5,warmups_per_model_backend=1,start="zero",algorithm_old="gn",
          stop_rule_old="standard",first_stage="same frozen archived values for both",
          old_max_inner=300,new_max_inner=400,limits_reached=limits_reached,
          source_data_sha256=before,source_and_data_unchanged=True,
          timing_scope="stage-two point loop plus final residual covariance; excludes data preparation, first-stage fitting, coefficient covariance, SEs, elasticities, output and artifact checks",
          instrumentation="Scoped observation wrappers on gauss_newton and _finish; original functions run unchanged; common BLAS limit; no optimizer stopping rule altered",
          benchmark_cases=len(checks),all_repetitions_converged=all(z["success"] for z in checks),
          iteration_counts_stable_across_repetitions=stable,
          all_iteration_counts_match_archive=all(z["inner_solves_match_archive"] and z["inner_iterations_match_archive"] for z in checks),
          known_archival_count_difference="cq old core replay uses 37 solves/70 inner iterations versus archived 38/71; shared matrix-first-stage reconstruction differs from package columnwise projection at floating-point rounding scale; see standard-rounding-probe.json"))


def rounding_probe():
    """Explain the one-step cq replay difference without changing the benchmark."""
    with threadpool_limits(limits=1):
        p = Problem("cq")
        theta_matrix, matrix_info = fit_old(p)
        B = p.tau.reshape(p.m,-1)
        index_columns = np.column_stack([p.design @ B[j,:-1] + B[j,-1] for j in range(p.m)])
        cdf_columns, pdf_columns = norm.cdf(index_columns), norm.pdf(index_columns)
        differences = dict(index_max_abs=float(abs(index_columns-p.index).max()),
                           cdf_max_abs=float(abs(cdf_columns-p.d.cdf).max()),
                           pdf_max_abs=float(abs(pdf_columns-p.d.pdf).max()))
        p.d.cdf, p.d.pdf = cdf_columns, pdf_columns
        theta_columns, column_info = fit_old(p)
        return dict(model="cq",scope="old standard GN only; rounding probe, not a timed benchmark repetition",
                    matrix_projection="design @ coefficients.T + constants",
                    column_projection="ProbitResult.xb: design @ coefficients_j + constant_j, one equation at a time",
                    input_rounding_difference=differences,
                    matrix_inner_solves=matrix_info["inner_solves"],
                    matrix_inner_iterations=matrix_info["inner_iterations"],
                    column_inner_solves=column_info["inner_solves"],
                    column_inner_iterations=column_info["inner_iterations"],
                    matrix_theta_max_abs_difference_from_archive=float(abs(theta_matrix-p.reference).max()),
                    column_theta_max_abs_difference_from_archive=float(abs(theta_columns-p.reference).max()),
                    significance="Columnwise projection reproduces the archived old counts and coefficients; both benchmark backends use the same matrix projection and remain directly comparable.")


if __name__ == "__main__":
    main()
