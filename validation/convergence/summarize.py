"""Export settings and summaries without changing any native success flags."""
import json
import pandas as pd
from .compare import DEFAULT_OUT


def main():
    rows = json.loads((DEFAULT_OUT / "runs.json").read_text())
    controls = []
    for r in rows:
        z = {k:r[k] for k in ("model","start","method","tightened","success","status")}
        keys = ("ftol","xtol","gtol","outer_xtol","eps","tau","delta","ifgnlseps",
                "max_inner","max_outer","max_nfev","scaling_anchor","message",
                "earlier_inner_failures")
        z.update({k:r.get(k) for k in keys})
        z.update({"profile_"+k:v for k,v in r.get("options",{}).items()})
        if r["history"]:
            last = r["history"][-1]
            z.update(final_inner_success=last["success"],
                     final_inner_status=last.get("status"),
                     final_inner_message=last.get("message"),
                     final_inner_parameter_pass=last.get("param_ok"),
                     final_inner_rss_pass=last.get("rss_ok"))
        if r.get("outer"):
            last = r["outer"][-1]
            z.update(final_outer_parameter_change=last["parameter_change"],
                     final_outer_covariance_change=last["covariance_change"],
                     final_outer_parameter_pass=last["parameter_pass"],
                     final_outer_covariance_pass=last["covariance_pass"])
        controls.append(z)
    pd.DataFrame(controls).to_csv(DEFAULT_OUT / "controls.csv", index=False)
    data = pd.DataFrame(rows)
    grouped = data.groupby(["method","start","tightened"], sort=True)
    summary = grouped.agg(configurations=("success","size"),native_successes=("success","sum"),
                 median_elapsed_seconds=("elapsed_seconds","median"),
                 max_rms_scaled_score=("rms_scaled_score_inf","max"),
                 max_fitted_difference_from_recorded=("fitted_diff_vs_recorded","max"),
                 max_elasticity_difference_from_recorded=("elasticity_diff_vs_recorded","max"))
    summary.to_csv(DEFAULT_OUT / "summary.csv")
    print(f"Exported settings for {len(controls)} fits and {len(summary)} summary groups")


if __name__ == "__main__":
    main()
