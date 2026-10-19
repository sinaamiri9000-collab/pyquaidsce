"""Combine the saved bootstrap batches and write the comparison tables."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
BENCHMARK = HERE.parents[1]
RESULTS = BENCHMARK / "results"
RAW = RESULTS / "raw"
MODELS = ("censored", "uncensored")
GROUPS = [1, 2, 4, 5, 6, 7, 8, 9, 10, 11]


def load_draws(model: str) -> pd.DataFrame:
    batches = []
    jobs = [("first", 1, 50)]
    jobs.extend(("additional", first, min(first + 24, 200))
                for first in range(51, 201, 25))
    for batch_name, first_rep, last_rep in jobs:
        requested = last_rep - first_rep + 1
        if batch_name == "first":
            prefix = f"{model}_batch_first"
        else:
            prefix = f"{model}_batch_additional_{first_rep:03d}_{last_rep:03d}"
        draws = pd.read_csv(RAW / f"{prefix}_draws.csv", float_precision="round_trip")
        summary = json.loads((RAW / f"{prefix}_summary.json").read_text())
        recorded_first = summary.get("first_replication_id", first_rep)
        if (recorded_first != first_rep or summary["requested"] != requested
                or summary["replication_ids"] != [first_rep, last_rep]
                or summary["model"] != model
                or not summary["full_sample_converged"]
                or summary["successful"] != len(draws)):
            raise ValueError(f"{model} {prefix}: batch does not match expected replication IDs")
        failed_local = {
            int(message.split(":", 1)[0].split()[1])
            for message in summary["failures"]
        }
        expected_ids = [first_rep + i - 1 for i in range(1, requested + 1)
                        if i not in failed_local]
        if len(expected_ids) != len(draws):
            raise ValueError(f"{model} {batch_name}: cannot map draws to replication IDs")
        if draws["replication"].tolist() != expected_ids:
            raise ValueError(f"{model} {prefix}: saved replication IDs do not match settings")
        batches.append(draws)
    draws = pd.concat(batches, ignore_index=True)
    if draws.replication.tolist() != list(range(1, 201)):
        raise ValueError(f"{model}: expected replication IDs 1 through 200 once each")
    if draws.isna().any().any() or not np.isfinite(draws.iloc[:, 1:].to_numpy(float)).all():
        raise ValueError(f"{model}: bootstrap draws contain missing or non-finite values")
    draws.to_csv(RESULTS / f"{model}_bootstrap_draws.csv", index=False,
                 float_format="%.17g")
    return draws


def main() -> None:
    all_draws = {model: load_draws(model) for model in MODELS}
    if not all_draws["censored"].replication.equals(all_draws["uncensored"].replication):
        raise ValueError("The two models must use matching replication IDs")
    for censored_path in sorted(RAW.glob("censored_*_summary.json")):
        uncensored_path = censored_path.with_name(censored_path.name.replace("censored_", "uncensored_", 1))
        censored = json.loads(censored_path.read_text())
        uncensored = json.loads(uncensored_path.read_text())
        for key in ("seed", "replication_ids", "requested", "observations",
                    "excluded_instruments", "original_group_numbers", "anot"):
            if censored[key] != uncensored[key]:
                raise ValueError(f"Paired batch settings differ for {key}: {censored_path.name}")
    for model, draws in all_draws.items():
        names = draws.columns.tolist()[1:]
        point = pd.read_csv(RESULTS / f"{model}_point_estimates.csv", float_precision="round_trip")
        if point.parameter.tolist() != names:
            raise ValueError(f"{model}: point and bootstrap parameter orders do not match")
        values = draws[names].to_numpy(float)
        inference = pd.DataFrame({
            "parameter": names,
            "estimate": point.estimate.to_numpy(float),
            "bootstrap_se": values.std(axis=0, ddof=1),
            "percentile_2_5": np.percentile(values, 2.5, axis=0),
            "percentile_97_5": np.percentile(values, 97.5, axis=0),
            "successful_replications": len(values),
        })
        inference.to_csv(RESULTS / f"{model}_bootstrap_inference.csv",
                         index=False, float_format="%.17g")

    labels = pd.read_csv(BENCHMARK / "data" / "group_labels.csv").set_index("original_group")
    point_tables = {
        model: pd.read_csv(RESULTS / f"{model}_point_elasticities.csv", float_precision="round_trip").set_index("original_group")
        for model in MODELS
    }
    rows = []
    dc = all_draws["censored"].set_index("replication")
    du = all_draws["uncensored"].set_index("replication")
    for position, group in enumerate(GROUPS, start=1):
        pc = point_tables["censored"].loc[group]
        pu = point_tables["uncensored"].loc[group]
        inc = f"ELAS_INC:e_{position}"
        own = f"ELAS_UNCOMP:e_{position}_{position}"
        for param in (inc, own):
            if param not in dc.columns or param not in du.columns:
                raise ValueError(f"Missing elasticity parameter {param}")
        c_inc = dc[inc].to_numpy(float)
        u_inc = du[inc].to_numpy(float)
        c_own = dc[own].to_numpy(float)
        u_own = du[own].to_numpy(float)
        c_lo, c_hi = np.percentile(c_inc, [2.5, 97.5])
        u_lo, u_hi = np.percentile(u_inc, [2.5, 97.5])
        cp_lo, cp_hi = np.percentile(c_own, [2.5, 97.5])
        up_lo, up_hi = np.percentile(u_own, [2.5, 97.5])
        diff_lo, diff_hi = np.percentile(c_own - u_own, [2.5, 97.5])
        rows.append({
            "original_group": group,
            "label": labels.loc[group, "label"],
            "expenditure_censored": pc.expenditure_elasticity,
            "expenditure_censored_2_5": c_lo,
            "expenditure_censored_97_5": c_hi,
            "expenditure_uncensored": pu.expenditure_elasticity,
            "expenditure_uncensored_2_5": u_lo,
            "expenditure_uncensored_97_5": u_hi,
            "own_price_censored": pc.own_price_marshallian,
            "own_price_censored_2_5": cp_lo,
            "own_price_censored_97_5": cp_hi,
            "own_price_uncensored": pu.own_price_marshallian,
            "own_price_uncensored_2_5": up_lo,
            "own_price_uncensored_97_5": up_hi,
            "own_price_difference_censored_minus_uncensored": (
                pc.own_price_marshallian - pu.own_price_marshallian),
            "own_price_difference_2_5": diff_lo,
            "own_price_difference_97_5": diff_hi,
        })
    compare = pd.DataFrame(rows)
    compare.to_csv(RESULTS / "elasticity_comparison.csv", index=False,
                   float_format="%.17g")

    def cell(value: float, lower: float, upper: float) -> str:
        return f"{value:.3f} [{lower:.3f}, {upper:.3f}]"

    lines = [
        "# Ten-good censored and uncensored comparison",
        "",
        "Point estimates with 95% percentile intervals from 200 package-native household bootstrap draws.",
        "Own-price elasticities are uncompensated (Marshallian).",
        "",
        "| Original group | Expenditure: censored | Expenditure: uncensored | Own price: censored | Own price: uncensored |",
        "|---|---:|---:|---:|---:|",
    ]
    for r in compare.itertuples(index=False):
        lines.append(
            f"| {r.original_group}. {r.label} | "
            f"{cell(r.expenditure_censored, r.expenditure_censored_2_5, r.expenditure_censored_97_5)} | "
            f"{cell(r.expenditure_uncensored, r.expenditure_uncensored_2_5, r.expenditure_uncensored_97_5)} | "
            f"{cell(r.own_price_censored, r.own_price_censored_2_5, r.own_price_censored_97_5)} | "
            f"{cell(r.own_price_uncensored, r.own_price_uncensored_2_5, r.own_price_uncensored_97_5)} |"
        )
    lines += [
        "",
        "The paired own-price difference is the censored estimate minus the uncensored estimate, using the same household resamples. See `elasticity_comparison.csv` for its point estimate and 95% percentile interval.",
        "",
        "These intervals use 200 replications and describe bootstrap variation for this sample and specification.",
    ]
    (RESULTS / "elasticity_comparison.md").write_text("\n".join(lines) + "\n")

    batch_summaries = {}
    for model in MODELS:
        jobs = [("first", 1, 50)]
        jobs.extend(("additional", first, min(first + 24, 200))
                    for first in range(51, 201, 25))
        records = []
        for batch_name, start, end in jobs:
            if batch_name == "first":
                prefix = f"{model}_batch_first"
            else:
                prefix = f"{model}_batch_additional_{start:03d}_{end:03d}"
            record = json.loads((RAW / f"{prefix}_summary.json").read_text())
            record["first_replication_id"] = start
            record["replication_ids"] = [start, end]
            records.append(record)
        total = sum(r["successful"] for r in records)
        failures = [x for r in records for x in r["failures"]]
        batch_summaries[model] = {
            "censor": model == "censored",
            "requested_total": sum(r["requested"] for r in records),
            "successful_total": total,
            "failures": failures,
            "observations": records[-1]["observations"],
            "package_version": records[-1]["package_version"],
            "resampling": records[-1]["resampling"],
            "reduced_form": records[-1]["reduced_form"],
            "excluded_instruments": records[-1]["excluded_instruments"],
            "batches": records,
        }
        if total != 200:
            raise ValueError(f"{model}: expected 200 successful replications, got {total}")
    (RESULTS / "bootstrap_manifest.json").write_text(
        json.dumps(batch_summaries, indent=2) + "\n")
    comparison_text = (RESULTS / "elasticity_comparison.md").read_text().rstrip()
    final_summary = [
        "Ten-good censored and uncensored QUAIDS demand models",
        "Data: 37,656 Iranian household observations.",
        "Bootstrap: 200 successful package-native draws per model.",
        "Resampling: household rows drawn independently with replacement; no cluster resampling.",
        "The internal ivexp reduced form is rebuilt in every draw.",
        "Intervals: 95% percentile intervals from the 200 draws.",
        "",
        comparison_text,
        "",
    ]
    (RESULTS / "logs" / "final_summary.txt").write_text("\n".join(final_summary))
    print(compare.to_string(index=False))
    print("Both models contain 200 successful replications.")
    checksum_path = RESULTS / "checksums.sha256"
    files = sorted(path for path in BENCHMARK.rglob("*")
                   if path.is_file() and path != checksum_path
                   and "__pycache__" not in path.parts)
    checksum_path.write_text("".join(
        f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(BENCHMARK).as_posix()}\n"
        for path in files
    ))


if __name__ == "__main__":
    main()
