"""Evidence tables and an honest scorecard, including failed contracts."""
from collections import Counter
from itertools import combinations
import inspect
import json
import re
from pathlib import Path
import numpy as np
from pyquaidsce import quaidsce, __all__ as public_names


def generate(out, records, existing, design, complete):
    from .run import write_csv, write_json, plain
    counts = Counter(x["status"] for x in records)
    flat, scores = [], []
    for r in records:
        failed = [x["name"] for x in r.get("checks", []) if not x["passed"]]
        result = r.get("result", {})
        scores.append(dict(id=r["id"], group=r["group"], status=r["status"],
                           expectation=r.get("expectation", ""),
                           elapsed_seconds=r.get("elapsed_seconds", 0),
                           fit_seconds=r.get("fit_seconds"),
                           converged=result.get("converged"), nobs=result.get("nobs"),
                           llf=result.get("llf"), n_gn=result.get("n_gn"), n_outer=result.get("n_outer"),
                           checks=len(r.get("checks", [])), failed_checks=failed,
                           exception_type=r.get("exception", {}).get("type"),
                           exception_message=r.get("exception", {}).get("message"),
                           reason=r.get("reason", "")))
        for check in r.get("checks", []):
            flat.append(dict(case_id=r["id"], group=r["group"], **check))
    write_csv(out / "scenarios.csv", scores)
    write_csv(out / "checks.csv", flat)
    failures = []
    for r in records:
        for failure in r.get("bootstrap", {}).get("failures", []):
            match = re.match(r"rep (\d+): (.*)", failure)
            failures.append(dict(case_id=r["id"], replication=int(match.group(1)) if match else None,
                                 reason=match.group(2) if match else failure,
                                 seed=r.get("kwargs", {}).get("seed"),
                                 requested=r.get("bootstrap", {}).get("requested"),
                                 injected=r["id"].startswith("lifecycle-")))
    write_csv(out / "bootstrap-failures.csv", failures,
              ["case_id", "replication", "reason", "seed", "requested", "injected"])
    # The case files carry complete raw values. Keep this index compact instead
    # of committing a second copy of every covariance and assertion.
    write_json(out / "results.json", [dict(id=r["id"], group=r["group"], status=r["status"],
                elapsed_seconds=r.get("elapsed_seconds", 0), case_file="cases/" + r["id"] + ".json") for r in records])
    comparisons = []
    for prefix in ("core", "followup-tight"):
        models = ("cq", "ca", "uq", "ua", "nq", "na") if prefix == "core" else ("cq",)
        for model in models:
            for method in ("nls", "fgnls", "ifgnls"):
                ids = [(f"core-{model}-{method}-{a}-{s}" if prefix == "core" else f"followup-tight-{method}-{a}-{s}")
                       for a in ("gn", "lm") for s in ("zero", "linear")]
                for first, second in combinations(ids, 2):
                    r1, r2 = next((r for r in records if r["id"] == first), None), next((r for r in records if r["id"] == second), None)
                    if not r1 or not r2 or "result" not in r1 or "result" not in r2:
                        continue
                    a, b = np.asarray(r1["result"]["b"]), np.asarray(r2["result"]["b"])
                    close = np.allclose(a, b, atol=1e-6, rtol=1e-5)
                    comparisons.append(dict(family=prefix, model=model, method=method, first=first, second=second,
                        both_converged=bool(r1["result"]["converged"] and r2["result"]["converged"]),
                        max_coefficient_difference=float(np.max(np.abs(a-b))),
                        coefficient_agreement="WITHIN_TOLERANCE" if close else "DIFFERS",
                        coefficient_atol=1e-6, coefficient_rtol=1e-5,
                        log_likelihood_difference=abs(r1["result"]["llf"]-r2["result"]["llf"]),
                        relative_SSE_difference=abs(r1.get("residual_sum_squares", 0)-r2.get("residual_sum_squares", 0))
                            / max(abs(r1.get("residual_sum_squares", 0)), 1e-12)))
    write_csv(out / "optimizer-comparisons.csv", comparisons)
    try:
        from .visualize import core_heatmap
        core_heatmap(out, records)
    except ImportError:
        pass  # The audit itself has no matplotlib runtime dependency.
    byid = {r["id"]: r for r in records}
    coverage = []
    guide = (Path(__file__).resolve().parents[1] / "docs" / "user-guide.md").read_text()
    for name, parameter in inspect.signature(quaidsce).parameters.items():
        if name.startswith("_"):
            continue
        mapped = [c.id for c in design if name in c.covers or name in c.kwargs]
        if name in ("data", "shares", "prices", "expenditure", "demographics", "anot", "verbose"):
            mapped = list(dict.fromkeys(["core-cq-ifgnls-gn-zero"] + mapped))
        run = [id for id in mapped if id in byid]
        failed = [id for id in run if byid[id]["status"] in ("FAIL", "NONCONVERGED", "ERROR", "BLOCKED")]
        default = "required" if parameter.default is inspect.Parameter.empty else repr(parameter.default)
        coverage.append(dict(argument=name, default=default, planned_cases=len(mapped),
                             executed_cases=len(run), case_ids=mapped, failed_cases=failed,
                             explicitly_documented_in_user_guide="`" + name + "`" in guide,
                             coverage="EXECUTED" if run else "NOT_RUN",
                             outcome="FAILURES_PRESENT" if failed else "NO_FAILURES_IN_MAPPED_CASES"))
    write_csv(out / "argument-coverage.csv", coverage)
    write_json(out / "argument-coverage.json", coverage)
    api = []
    mapping = {
        "quaidsce": "scenario matrix", "first_stage": "public-first_stage",
        "nlsur": "public-nlsur", "probit": "public-probit_no_constant",
        "Spec": "existing tests.test_suite.InternalMathTests", "Coefs": "matrix/model-point oracles",
        "DemandData": "matrix/oracle reconstruction", "FirstStageLayout": "selection matrix/public-first_stage",
        "ExpenditureReducedForm": "control_function QR oracle",
        "fit_expenditure_reduced_form": "control_function QR oracle/public-reduced_form_labels",
        "QuaidsceResults": "matrix/summary contracts", "Elasticities": "matrix/model-point oracles",
        "elasticities": "model-point finite differences", "fitted_share_derivatives": "observation-level finite differences/public-stale_derivative_design",
        "sample_means": "matrix means/public-sample_means", "unpack": "finite-difference mapping/independent prediction",
        "full_vector": "reported parameter mapping/public-delta_mapping", "delta_matrix": "public-delta_mapping",
        "fitted_shares": "independent direct formula", "latent_shares": "public-latent_shares",
        "augmented_latent_shares": "public-latent_shares", "jacobian_full": "matrix independent numerical Jacobian",
        "__version__": "environment metadata/existing suite",
    }
    for name in public_names:
        api.append(dict(export=name, evidence=mapping.get(name, "UNMAPPED")))
    write_csv(out / "public-api-coverage.csv", api)
    stability = []
    for model in ("base", "iv"):
        paths = [out / "bootstrap" / f"boot-large-{model}-{seed}-draws.csv" for seed in (16001, 16002)]
        if not all(p.exists() for p in paths):
            continue
        import pandas as pd
        frames = [pd.read_csv(p).drop(columns="draw_number") for p in paths]
        se1, se2 = (f.std(ddof=1) for f in frames)
        for name in frames[0].columns:
            denominator = max(float((se1[name] + se2[name]) / 2), 1e-12)
            stability.append(dict(model=model, name=name, se_seed16001=se1[name], se_seed16002=se2[name],
                                  relative_se_difference=abs(se1[name]-se2[name])/denominator))
    write_csv(out / "bootstrap-stability.csv", stability,
              ["model", "name", "se_seed16001", "se_seed16002", "relative_se_difference"])
    summary = dict(complete_design_executed=complete, scenario_counts=dict(counts),
                   scenarios=len(records), checks=len(flat), failed_checks=sum(not x["passed"] for x in flat),
                   existing_suite=existing, public_arguments=len(coverage),
                   arguments_executed=sum(x["coverage"] == "EXECUTED" for x in coverage),
                   total_scenario_seconds=sum(x.get("elapsed_seconds", 0) for x in records),
                   bootstrap_requested=sum(x.get("bootstrap", {}).get("requested", 0) for x in records),
                   bootstrap_successful=sum(x.get("bootstrap", {}).get("successful", 0) for x in records),
                   ordinary_bootstrap_requested=sum(x.get("bootstrap", {}).get("requested", 0) for x in records
                                                    if x["group"] in ("bootstrap", "bootstrap_large", "bootstrap_supplemental")),
                   ordinary_bootstrap_successful=sum(x.get("bootstrap", {}).get("successful", 0) for x in records
                                                     if x["group"] in ("bootstrap", "bootstrap_large", "bootstrap_supplemental")),
                   large_bootstrap_requested=sum(x.get("bootstrap", {}).get("requested", 0) for x in records if x["group"] == "bootstrap_large"),
                   large_bootstrap_successful=sum(x.get("bootstrap", {}).get("successful", 0) for x in records if x["group"] == "bootstrap_large"),
                   native_interfaces="not executed; Python bridge execution and static R contracts are distinct evidence",
                   release_gate_passed=complete and bool(existing.get("successful"))
                       and not any(x["status"] in ("FAIL", "NONCONVERGED", "ERROR", "BLOCKED") for x in records)
                       and all(x["coverage"] == "EXECUTED" for x in coverage))
    write_json(out / "summary.json", summary)
    text = ["# pyquaidsce 1.6.0 small4 validation results", "",
            "This is an audit of the unchanged 1.6.0 estimator, not a comparison with the original Stata ado.", "",
            "The interpreted findings are in [FINDINGS.md](../FINDINGS.md) and [FINDINGS.fa.md](../FINDINGS.fa.md).", "",
            f"- Scenarios: {len(records)}; statuses: `{dict(counts)}`.",
            f"- Assertions: {len(flat)}; failed: {summary['failed_checks']}.",
            f"- Existing suite: `{existing}`.",
            f"- Public arguments with mapped executions: {summary['arguments_executed']}/{len(coverage)}.",
            f"- Full selected design executed: {complete}.",
            f"- Release validation gate passed: **{summary['release_gate_passed']}**.", "",
            "## How to read the evidence", "",
            "`manifest.json` freezes expectations; `scenarios.csv` is the scorecard; `checks.csv` contains observed errors and fixed thresholds. "
            "`cases/` contains configurations, diagnostics, coefficient/SE vectors, theta, Sigma, and structural covariance. "
            "`coefficients/` provides estimates in CSV. `bootstrap/` contains all successful draws and active covariance matrices. "
            "Analytical elasticity SEs are encoded as `NaN` by contract. Small-bootstrap covariance is allowed to be singular.", "",
            "No tolerance was increased to turn a failed case green. A returned nonconverged result is not counted as PASS "
            "unless that exact case deliberately tests exhaustion. Injected lifecycle failures are explicitly labeled.", "",
            "`optimizer-comparisons.csv` contrasts algorithms/starts for the same specification. Its difference labels are "
            "diagnostics, not additional pass/fail assertions: default standard stopping is loose, and FGNLS weights depend "
            "on the preliminary NLS solution. Separate tight-stop follow-ups test stationarity and retain the original failures.", "",
            "## Failed or incomplete cases", "",
            "| Case | Status | Failed checks / reason |", "|---|---|---|"]
    for row in scores:
        if row["status"] in ("FAIL", "NONCONVERGED", "ERROR", "BLOCKED", "NOT_RUN"):
            text.append(f"| {row['id']} | {row['status']} | {', '.join(row['failed_checks']) or row['reason']} |")
    text += ["", "## Bootstrap precision", "",
             "The four large batches use 100 draws each, as requested. `bootstrap-stability.csv` compares SEs across two fixed seeds "
             "for the base and IV models. These are diagnostics, not a proof of frequentist coverage or instrument validity. "
             "No threshold on SE drift is retroactively chosen to declare statistical convergence.", "",
             "A full-success bootstrap stress gate is deliberately stricter than the documented ability to discard failed draws. "
             "A failed stress gate does not imply that covariance/SE accounting was implemented incorrectly; the individual "
             "assertions distinguish the two. Successful-row ordinals are not original failed/successful replication IDs.", "",
             "## Scope", "",
             "small4 provides 2,000 observations, four goods, and two demographics. Subsystems and derived columns are recorded; "
             "the original DTA is unchanged. Existing synthetic/mathematical tests supplement the empirical matrix. "
             "Native R and licensed Stata execution are not represented as tested. This run does not establish behavior on every "
             "dependency version, platform, instrument design, or number of goods.", ""]
    (out / "REPORT.md").write_text("\n".join(text), encoding="utf-8")
    print(json.dumps(plain(summary), indent=2), flush=True)
