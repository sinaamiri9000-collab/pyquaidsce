"""Post hoc diagnostics: all engines start from the same discovered point.

These probe cases remain separate from the declared primary experiment.
They test whether different terminal points are retained by other engines.
They are not a default initialization strategy or evidence of global maxima.
"""
import json
from threadpoolctl import threadpool_limits
from .compare import DEFAULT_OUT, METHODS, CSV_FIELDS, dump, experiment, write_csv


def main():
    primary = json.loads((DEFAULT_OUT / "runs.json").read_text())
    selections = [("cq", "zero", "profile-BFGS"),
                  ("cq", "zero", "scipy-trf"),
                  ("ca", "zero", "scipy-trf"),
                  ("uq", "zero", "scipy-trf")]
    rows = []
    with threadpool_limits(limits=1):
        for model, start, method in selections:
            source = next(r for r in primary if r["model"] == model
                          and r["start"] == start and r["method"] == method
                          and not r["tightened"])
            for backend in METHODS:
                row = experiment(model, "cross-"+method, backend,
                                 warm_theta=source["theta"])
                row.update(source_start=start, source_method=method,
                           source_success=source["success"], source_ll=source["gaussian_ll"],
                           ll_change_from_shared_start=row["gaussian_ll"]-source["gaussian_ll"])
                rows.append(row)
                print(json.dumps({k:row.get(k) for k in
                     ("model", "start", "method", "success", "gaussian_ll",
                      "ll_change_from_shared_start", "rms_scaled_score_inf")}), flush=True)
    dump(DEFAULT_OUT / "cross-start.json", rows)
    write_csv(DEFAULT_OUT / "cross-start.csv", rows,
              CSV_FIELDS+("source_start", "source_method", "source_success", "source_ll",
                          "ll_change_from_shared_start"))


if __name__ == "__main__":
    main()
