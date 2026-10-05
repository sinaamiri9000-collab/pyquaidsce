"""Frozen audit design: every case has an expectation before execution."""
from dataclasses import dataclass, field
from itertools import product
import multiprocessing as mp

SHARES = ["sw1", "sw2", "sw4", "sw9"]
PRICES = ["p1", "p2", "p4", "p9"]
DEMOS = ["x1", "x2"]
BASE = dict(shares=SHARES, prices=PRICES, expenditure="total",
            demographics=DEMOS, anot=10.0, verbose=False)
BASE_ID = "core-cq-ifgnls-gn-zero"


@dataclass
class Case:
    id: str
    group: str
    kwargs: dict = field(default_factory=dict)
    fixture: str = "original"
    checks: list = field(default_factory=list)
    reference: str = ""
    expectation: str = "converged fit and independent numerical checks"
    exception: tuple = ()
    message: str = ""
    covers: list = field(default_factory=list)


def scenarios(large_reps=100):
    cases = []
    models = {
        "cq": {}, "ca": {"quadratic": False},
        "uq": {"censor": False}, "ua": {"censor": False, "quadratic": False},
        "nq": {"censor": False, "demographics": None},
        "na": {"censor": False, "quadratic": False, "demographics": None},
    }
    for model, settings in models.items():
        for method, algorithm, start in product(
                ("nls", "fgnls", "ifgnls"), ("gn", "lm"), ("zero", "linear")):
            cases.append(Case(f"core-{model}-{method}-{algorithm}-{start}", "core",
                              dict(settings, method=method, algorithm=algorithm, start=start),
                              covers=["quadratic", "censor", "demographics", "method",
                                      "algorithm", "start", "data", "shares", "prices",
                                      "expenditure", "anot", "reps"]))

    def add(id, group, kw=None, fixture="original", checks=None,
            reference="", expectation="converged fit and independent numerical checks",
            exception=(), message="", covers=None):
        cases.append(Case(id, group, kw or {}, fixture, checks or [], reference,
                          expectation, exception, message,
                          list(dict.fromkeys((covers or []) + list((kw or {}).keys())))))

    # 18 input/sample cases. Transform logs only after conversion to float64.
    for censored, reference in ((True, BASE_ID), (False, "core-uq-ifgnls-gn-zero")):
        for encoding in ("price_logs", "expenditure_logs", "both_logs"):
            add(f"input-{encoding}-c{int(censored)}", "inputs",
                dict(censor=censored), encoding, ["equivalent"], reference,
                covers=["lnprices", "lnexpenditure"])
    for id, kw, fixture, check, ref in [
        ("renamed", {}, "renamed", "equivalent", BASE_ID),
        ("irrelevant", {}, "irrelevant", "equivalent", BASE_ID),
        ("reverse_rows", {}, "reverse_rows", "equivalent", BASE_ID),
        ("custom_index", {}, "custom_index", "equivalent", BASE_ID),
        ("permuted_goods", {}, "permuted_goods", "permuted_goods", BASE_ID),
        ("permuted_demos", {"demographics": list(reversed(DEMOS))}, "original", "permuted_demos", BASE_ID),
        ("one_demo", {"demographics": ["x1"]}, "original", "", ""),
        ("three_goods_censored", {}, "three_goods", "", ""),
        ("three_goods_uncensored", {"censor": False}, "three_goods", "", ""),
        ("missing_required", {}, "missing_required", "complete_case", ""),
        ("missing_irrelevant", {}, "missing_irrelevant", "equivalent", BASE_ID),
        ("infinite_required", {}, "infinite_required", "complete_case", ""),
    ]:
        add(f"input-{id}", "inputs", kw, fixture, [check] if check else [], ref,
            covers=["data", "shares", "prices", "demographics"])

    # 18 first-stage designs: None, [], subset, reordering, independent z,
    # True/False/None expenditure, and intercept-only.
    selection = [
        ("explicit_defaults", dict(selection_prices=PRICES, selection_covariates=DEMOS,
                                   selection_expenditure=True), "original", "equivalent", BASE_ID),
        ("price_subset", dict(selection_prices=["p1", "p4"]), "original", "", ""),
        ("price_reorder", dict(selection_prices=list(reversed(PRICES))), "original", "selection_reorder", BASE_ID),
        ("no_prices", dict(selection_prices=[]), "original", "", ""),
        ("covariate_subset", dict(selection_covariates=["x1"]), "original", "", ""),
        ("covariate_reorder", dict(selection_covariates=list(reversed(DEMOS))), "original", "selection_reorder", BASE_ID),
        ("independent_covariate", dict(selection_covariates=["s1"]), "derived", "", ""),
        ("no_covariates", dict(selection_covariates=[]), "original", "", ""),
        ("no_expenditure", dict(selection_expenditure=False), "original", "", ""),
        ("none_raw", dict(selection_expenditure=None), "original", "equivalent", BASE_ID),
        ("none_logged", dict(selection_expenditure=None), "both_logs", "", ""),
        ("false_logged", dict(selection_expenditure=False), "both_logs", "equivalent", "selection-none_logged"),
        ("intercept_only", dict(selection_prices=[], selection_covariates=[], selection_expenditure=False), "original", "", ""),
        ("linear_subset", dict(quadratic=False, selection_prices=["p4", "p1"], selection_expenditure=False), "original", "", ""),
        ("independent_no_exp", dict(selection_covariates=["s1", "s2"], selection_expenditure=False), "derived", "", ""),
        ("price_only", dict(selection_covariates=[], selection_expenditure=False), "original", "", ""),
        ("exp_only", dict(selection_prices=[], selection_covariates=[]), "original", "", ""),
        ("collinear_covariates", dict(selection_covariates=["x1", "x1_copy"]), "derived", "dropped_column", ""),
    ]
    for id, kw, fixture, check, ref in selection:
        add(f"selection-{id}", "selection", kw, fixture, [check] if check else [], ref,
            covers=["selection_prices", "selection_covariates", "selection_expenditure"])

    cf = [
        ("demand", dict(control_function="cf"), "", ""),
        ("participation", dict(selection_control_function="cf"), "", ""),
        ("both_same", dict(control_function="cf", selection_control_function="cf"), "", ""),
        ("both_distinct", dict(control_function="cf", selection_control_function="cf2"), "", ""),
        ("linear_both", dict(quadratic=False, control_function="cf", selection_control_function="cf2"), "", ""),
        ("iv_one", dict(ivexp=["iv1"]), "reduced_form", ""),
        ("iv_two", dict(ivexp=["iv1", "iv2"]), "reduced_form", ""),
        ("manual_one", dict(control_function="cf", selection_control_function="cf"), "equivalent", "cf-iv_one"),
        ("iv_logged", dict(ivexp=["iv1"]), "equivalent", "cf-iv_one"),
        ("iv_custom_selection", dict(ivexp=["iv1"], selection_prices=["p4", "p1"],
                                    selection_expenditure=False, selection_covariates=["s1"]), "reduced_form", ""),
        ("iv_linear", dict(ivexp=["iv1"], quadratic=False), "reduced_form", ""),
        ("residual_rescaled", dict(control_function="cf_scaled", selection_control_function="cf_scaled"),
         "residual_rescaled", "cf-both_same"),
    ]
    for id, kw, check, ref in cf:
        add(f"cf-{id}", "control_function", kw,
            "derived_logs" if id == "iv_logged" else "derived", [check] if check else [], ref,
            covers=["ivexp", "control_function", "selection_control_function"])

    numerical = [
        ("initial", {}, ["initial"], BASE_ID),
        ("initial_sigma", {}, ["initial_sigma"], BASE_ID),
        ("nls_to_ifgnls", {}, ["initial_sigma"], "core-cq-nls-gn-zero"),
        ("initial_precedence", {"start": "linear"}, ["initial"], BASE_ID),
        ("sigma_only", {}, ["sigma_only"], BASE_ID),
        ("chunk_1", {"chunk": 1}, ["equivalent"], BASE_ID),
        ("chunk_127", {"chunk": 127}, ["equivalent"], BASE_ID),
        ("chunk_oversize", {"chunk": 10000}, ["equivalent"], BASE_ID),
        ("blas_2", {"blas_threads": 2}, ["equivalent", "blas_restore"], BASE_ID),
        ("blas_none", {"blas_threads": None}, ["equivalent", "blas_restore"], BASE_ID),
        ("vce_final", {"vce_sigma": "final"}, ["equivalent_point", "covariance"], BASE_ID),
        ("fgnls_vce_final", {"method": "fgnls", "vce_sigma": "final"}, ["equivalent_point", "covariance"], "core-cq-fgnls-gn-zero"),
        ("nls_vce_final", {"method": "nls", "vce_sigma": "final"}, ["equivalent", "covariance"], "core-cq-nls-gn-zero"),
        ("tight", {"stop_rule": "tight"}, ["stationarity"], ""),
        ("tight_tolerances", {"stop_rule": "tight", "tol": 1e-14, "nrtol_stop": 1e-13,
                              "sigma_tol": 1e-6, "inner_nrtol_early": 1e-10}, ["stationarity"], ""),
        ("loose_tolerances", {"tol": 1e-8, "nrtol_stop": 1e-8, "sigma_tol": 1e-4,
                              "inner_nrtol_early": 1e-6}, [], ""),
        ("inner_exhaustion", {"method": "nls", "max_iter": 1}, ["nonconvergence"], ""),
        ("outer_exhaustion", {"max_outer": 2}, ["nonconvergence"], ""),
    ]
    for id, kw, checks, ref in numerical:
        add(f"numeric-{id}", "numerical", kw, checks=checks, reference=ref,
            expectation="nonconverged result with explanatory note" if "nonconvergence" in checks
            else "converged fit and independent numerical checks",
            covers=["initial", "sigma_initial"] if id.startswith("initial") or id == "nls_to_ifgnls" else None)

    interactions = [
        ("iv_lm_linear", dict(ivexp=["iv1"], algorithm="lm", start="linear"), "derived", []),
        ("linear_lm_selection", dict(quadratic=False, algorithm="lm", selection_prices=["p4", "p1"], selection_expenditure=False), "original", []),
        ("cf_lm", dict(control_function="cf", selection_control_function="cf2", algorithm="lm"), "derived", []),
        ("uncensored_logs_tight", dict(censor=False, stop_rule="tight"), "both_logs", []),
        ("no_demo_linear_logs", dict(censor=False, quadratic=False, demographics=None), "both_logs", []),
        ("three_goods_selection", dict(selection_prices=["p4", "p1"], selection_covariates=[], quadratic=False), "three_goods", []),
        ("anot_8", dict(anot=8.0), "original", []),
        ("anot_12", dict(anot=12.0), "original", []),
        ("uppercase", dict(method="IFGNLS", algorithm="GN", start="ZERO", vce_sigma="OBJECTIVE", stop_rule="STANDARD", bootstrap_start="ZERO"), "original", ["equivalent"]),
        ("verbose", dict(verbose=True), "original", ["equivalent", "verbose"]),
        ("gn_verbose", dict(gn_verbose=True), "original", ["equivalent", "gn_verbose"]),
        ("callback", {}, "original", ["equivalent", "callback"]),
    ]
    for id, kw, fixture, checks in interactions:
        add(f"interaction-{id}", "interactions", kw, fixture, checks,
            BASE_ID if "equivalent" in checks else "", covers=["log"] if id == "callback" else None)
    assert len(cases) == 150

    # Exactly 12 smoke batches (8 draws) and 4 publication batches (100 draws).
    smoke = [
        ("serial", {}, "original", ""),
        ("repeat", {}, "original", "boot-serial"),
        ("parallel", {"n_jobs": 2}, "original", "boot-serial"),
        ("spawn", {"n_jobs": 2, "mp_context": "spawn"}, "original", "boot-parallel"),
        ("warm", {"bootstrap_start": "warm"}, "original", "boot-serial"),
        ("blas_2", {"n_jobs": 2, "blas_threads": 2}, "original", "boot-parallel"),
        ("blas_none", {"blas_threads": None}, "original", "boot-serial"),
        ("forkserver", {"n_jobs": 2, "mp_context": "forkserver"}, "original", "boot-parallel"),
        ("iv", {"ivexp": ["iv1"]}, "derived", ""),
        ("uncensored", {"censor": False}, "original", ""),
        ("linear", {"quadratic": False, "method": "fgnls"}, "original", ""),
        ("sigma_tol", {"boot_sigma_tol": 1e-6}, "original", ""),
    ]
    for id, kw, fixture, ref in smoke:
        add(f"boot-{id}", "bootstrap", dict(reps=8, seed=16001, n_jobs=1, **kw)
            if "n_jobs" not in kw else dict(reps=8, seed=16001, **kw), fixture,
            ["bootstrap"] + (["bootstrap_compare"] if ref else []) + (["count_reduced_forms"] if id == "iv" else []), ref,
            covers=["reps", "seed", "n_jobs", "mp_context", "bootstrap_start", "boot_sigma_tol"])
    for model, seed in product(("base", "iv"), (16001, 16002)):
        kw = dict(reps=large_reps, seed=seed, n_jobs=2, mp_context="spawn")
        if model == "iv":
            kw["ivexp"] = ["iv1"]
        add(f"boot-large-{model}-{seed}", "bootstrap_large", kw,
            "derived" if model == "iv" else "original", ["bootstrap"],
            covers=["reps", "seed", "n_jobs", "ivexp"])

    def negative(id, kw=None, fixture="original", message="", exception=("ValueError", "TypeError")):
        add(f"invalid-{id}", "invalid", kw, fixture, expectation="explicit rejection of invalid input",
            exception=exception, message=message,
            covers=["data"] if fixture != "original" else None)

    for arg, values in {
        "method": ["bad"], "algorithm": ["bad"], "start": ["bad"],
        "vce_sigma": ["bad"], "stop_rule": ["bad"], "bootstrap_start": ["bad"],
        "anot": [float("nan"), float("inf")],
        "max_iter": [0, -1, 1.5], "max_outer": [1, -1, 2.5], "chunk": [0, -1, 1.5],
        "tol": [0, -1, float("nan"), float("inf")],
        "nrtol_stop": [0, -1, float("nan"), float("inf")],
        "sigma_tol": [0, -1, float("nan"), float("inf")],
        "inner_nrtol_early": [0, -1, float("nan"), float("inf")],
        "boot_sigma_tol": [0, -1, float("nan"), float("inf")],
        "blas_threads": [0, -1, 1.5, True],
        "rep_timeout": [0, -1, float("nan"), float("inf")],
        "reps": [-1, 2.5],
        "n_jobs": [0, -1, 1.5], "seed": [-1, 1.5],
    }.items():
        for index, value in enumerate(values):
            kw = {arg: value}
            if arg in ("boot_sigma_tol", "n_jobs", "seed"):
                kw["reps"] = 2
            negative(f"{arg}-{index}", kw)
    for id, kw, fixture, msg in [
        ("both_prices", {"lnprices": ["lp1", "lp2", "lp4", "lp9"]}, "original", "exactly one"),
        ("no_prices", {"prices": None}, "original", "exactly one"),
        ("both_expenditure", {"lnexpenditure": "lm"}, "original", "exactly one"),
        ("no_expenditure", {"expenditure": None}, "original", "exactly one"),
        ("two_goods", {"shares": SHARES[:2], "prices": PRICES[:2]}, "original", "at least 3"),
        ("price_count", {"prices": PRICES[:3]}, "original", "number of price"),
        ("duplicate_shares", {"shares": [SHARES[0]] * 4}, "original", "duplicate"),
        ("duplicate_prices", {"prices": [PRICES[0]] * 4}, "original", "duplicate"),
        ("duplicate_demos", {"demographics": ["x1", "x1"]}, "original", "duplicate"),
        ("no_demos", {"demographics": []}, "original", "demographic"),
        ("missing_column", {"expenditure": "absent"}, "original", "not found"),
        ("negative_share", {}, "negative_share", "nonnegative"),
        ("zero_price", {}, "zero_price", "nonpositive"),
        ("negative_price", {}, "negative_price", "nonpositive"),
        ("zero_expenditure", {}, "zero_expenditure", "nonpositive"),
        ("empty_sample", {}, "empty_sample", "no complete"),
        ("no_zero_share", {}, "no_zero_share", "no censoring"),
        ("all_zero_share", {}, "all_zero_share", "no variation"),
        ("uncensored_sum", {"censor": False}, "bad_sum", "sum to one"),
        ("selection_nonprice", {"selection_prices": ["x1"]}, "original", "subset"),
        ("selection_duplicates", {"selection_prices": ["p1", "p1"]}, "original", "duplicate"),
        ("selection_z_duplicates", {"selection_covariates": ["x1", "x1"]}, "original", "duplicate"),
        ("selection_exp_type", {"selection_expenditure": "yes"}, "original", "True, False"),
        ("selection_uncensored", {"censor": False, "selection_prices": []}, "original", "censor=True"),
        ("cf_uncensored", {"censor": False, "control_function": "cf"}, "derived", "censor=True"),
        ("cf_constant", {"control_function": "constant_cf"}, "derived", "variation"),
        ("cf_nan", {"control_function": "cf"}, "cf_nan", "finite"),
        ("cf_demo", {"control_function": "x1"}, "original", "demographic"),
        ("cf_z", {"control_function": "cf", "selection_covariates": ["cf"]}, "derived", "selection_control_function"),
        ("cf_collinear", {"selection_control_function": "lm"}, "derived", "collinear"),
        ("cf_bootstrap", {"control_function": "cf", "reps": 2}, "derived", "disabled"),
        ("selection_cf_bootstrap", {"selection_control_function": "cf", "reps": 2}, "derived", "disabled"),
        ("iv_empty", {"ivexp": []}, "derived", "at least one"),
        ("iv_string", {"ivexp": "iv1"}, "derived", "sequence"),
        ("iv_duplicate", {"ivexp": ["iv1", "iv1"]}, "derived", "duplicate"),
        ("iv_overlap", {"ivexp": ["p1"]}, "derived", "excluded"),
        ("iv_selection_overlap", {"ivexp": ["iv1"], "selection_covariates": ["iv1"]}, "derived", "excluded"),
        ("iv_external", {"ivexp": ["iv1"], "control_function": "cf"}, "derived", "cannot be combined"),
        ("iv_uncensored", {"ivexp": ["iv1"], "censor": False}, "derived", "censor=True"),
        ("iv_rank", {"ivexp": ["x1_copy"]}, "derived", "rank deficient"),
        ("iv_perfect_fit", {"ivexp": ["lm"]}, "derived", "perfectly"),
        ("initial_length", {"initial": [0.0, 0.0]}, "original", "initial"),
        ("initial_nan", {"initial": [float("nan")] * 30}, "original", "finite"),
        ("sigma_shape", {"sigma_initial": [[1.0]]}, "original", "shape"),
        ("sigma_nan", {"sigma_initial": [[float("nan")] * 4] * 4}, "original", "finite"),
        ("mp_context", {"reps": 2, "n_jobs": 2, "mp_context": "absent"}, "original", "not supported"),
    ]:
        negative(id, kw, fixture, msg)
    for removed in ("first_stage_predict", "strict_stata"):
        negative(removed, {removed: True}, exception=("TypeError",))
    # Lifecycle fault scenarios are not input-contract tests.
    add("lifecycle-one_draw", "lifecycle", {"reps": 1, "seed": 16001},
        expectation="one successful draw cannot produce standard errors",
        exception=("RuntimeError",), message="fewer than two", covers=["reps"])
    add("lifecycle-timeout", "lifecycle", {"reps": 2, "seed": 16001, "rep_timeout": 1e-6,
                                          "mp_context": "spawn"}, expectation="watchdog rejects all timed-out draws",
        exception=("RuntimeError",), message="TimeoutError", covers=["rep_timeout", "mp_context"])
    add("lifecycle-partial_failure", "lifecycle", {"reps": 4, "seed": 16001},
        checks=["bootstrap", "partial_failure"], expectation="one injected failed draw is excluded and reported",
        covers=["reps", "seed"])
    add("lifecycle-all_failure", "lifecycle", {"reps": 2, "seed": 16001},
        expectation="all injected failed draws raise an explicit error",
        exception=("RuntimeError",), message="every bootstrap", covers=["reps"])
    return cases
