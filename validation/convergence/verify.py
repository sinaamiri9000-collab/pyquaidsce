"""Check experimental artifacts against independent equations and identities.

These checks do not replace any optimizer's convergence decision.
"""
import json
import numpy as np
from threadpoolctl import threadpool_limits
from validation.oracles import independent_fitted
from .compare import DEFAULT_OUT, Problem, dump, source_hashes, whiten, write_csv


def main():
    primary = json.loads((DEFAULT_OUT / "runs.json").read_text())
    cross = json.loads((DEFAULT_OUT / "cross-start.json").read_text())
    env = json.loads((DEFAULT_OUT / "environment.json").read_text())
    probes = json.loads((DEFAULT_OUT / "gradient-probes.json").read_text())
    checks = []
    def check(name, ok, observed=None, tolerance=None):
        checks.append(dict(name=name, passed=bool(ok), observed=observed, tolerance=tolerance))
    check("primary_case_count", len(primary) == 120, len(primary), 120)
    check("cross_start_case_count", len(cross) == 24, len(cross), 24)
    check("linear_start_excluded", all(r["start"] in ("zero","recorded") for r in primary))
    check("package_and_dataset_hashes_unchanged", source_hashes() == env["source_data_sha256"])
    for z in probes:
        if z["h"] == 1e-7:
            check(f"profile_gradient:{z['model']}:{z['point']}:{z['direction']}",
                  z["scaled_error"] < 2e-6, z["scaled_error"], 2e-6)
    with threadpool_limits(limits=1):
        for i, row in enumerate(primary+cross):
            p = Problem(row["model"])
            theta = np.array(row["theta"])
            # Independently written, row-wise model equations on a fixed sample.
            small = p.d.subset(np.linspace(0,p.N-1,12,dtype=int))
            independent = independent_fitted(theta, small, p.spec)
            actual = small.shares[:,:p.m] - p.u(theta)[np.linspace(0,p.N-1,12,dtype=int)]
            err = float(np.max(np.abs(actual-independent)/np.maximum(1,abs(independent))))
            check(f"independent_predictions:{i}", err < 1e-9, err, 1e-9)
            u = p.u(theta)
            S = u.T @ u / p.N
            errS = float(np.max(np.abs(S-np.array(row["sigma"]))))
            check(f"residual_covariance:{i}", errS < 1e-12, errS, 1e-12)
            ll = -p.N/2*(p.m*(1+np.log(2*np.pi))+np.linalg.slogdet(S)[1])
            check(f"concentrated_likelihood:{i}", abs(ll-row["gaussian_ll"]) < 1e-7,
                  abs(ll-row["gaussian_ll"]), 1e-7)
            P = whiten(S)
            Q = float(np.sum((u @ P.T)**2))
            check(f"own_weight_Q_is_Nm:{i}", abs(Q-p.N*p.m) < 1e-7,
                  abs(Q-p.N*p.m), 1e-7)
            Us, sv, Vt = np.linalg.svd(S)
            psi = (Vt.T * (1/np.sqrt(sv))) @ Us.T
            Qsvd = float(np.sum((u @ psi)**2))
            check(f"source_SVD_equals_Cholesky_weighting:{i}", abs(Q-Qsvd) < 1e-7,
                  abs(Q-Qsvd), 1e-7)
    dump(DEFAULT_OUT / "verification.json", checks)
    write_csv(DEFAULT_OUT / "verification.csv", checks, list(checks[0]))
    passed = sum(z["passed"] for z in checks)
    print(f"{passed}/{len(checks)} experimental verification checks passed")
    if passed != len(checks):
        for z in checks:
            if not z["passed"]:
                print(z)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
