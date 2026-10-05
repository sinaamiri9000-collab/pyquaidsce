"""Controlled small4 convergence experiments with three estimator architectures.

Run from the repository root with PYTHONPATH=src:
    python -m validation.convergence.compare

Stata is NOT executed. Its uploaded ado's stopping and update logic is
independently reconstructed, with LAPACK QR and documented Mata tolerance.
SciPy's own convergence decisions are used without replacing their predicates.
Failures, invalid trials and unsuccessful inner solves remain in the output.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import inspect
import json
from pathlib import Path
import platform
import time
from types import SimpleNamespace
import warnings

import numpy as np
import pandas as pd
import scipy
from scipy.linalg import qr, solve_triangular
from scipy.optimize import fixed_point, least_squares, minimize
from threadpoolctl import threadpool_limits

from pyquaidsce.elasticities import elasticities, sample_means
from pyquaidsce.jacfree import jacobian_free, make_cache
from pyquaidsce.model import fitted_shares, residuals
from pyquaidsce.params import Spec, unpack
from validation.oracles import build_data

ROOT = Path(__file__).resolve().parents[2]
DEFAULT_OUT = Path(__file__).resolve().parent / "results"
MODELS = ("cq", "ca", "uq", "ua", "nq", "na")
METHODS = ("scipy-trf", "scipy-lm", "profile-BFGS", "profile-L-BFGS-B",
           "nlsur-fd", "nlsur-analytic")


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def source_hashes():
    paths = sorted((ROOT / "src").rglob("*.py")) + [ROOT / "bench/small4.dta"]
    return {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}


def setup(model):
    record = json.loads((ROOT / f"validation/results/cases/core-{model}-ifgnls-gn-zero.json").read_text())
    r, kw = record["result"], record["kwargs"]
    spec = Spec(4, ndemo=0 if model.startswith("n") else 2,
                quadratic=model.endswith("q"), censor=model.startswith("c"))
    tau = np.array(r["b"][spec.n_full:-36]) if spec.censor else np.empty(0)
    fake = SimpleNamespace(spec=spec, tau=tau, reduced_form=None)
    d, design, index, _ = build_data(pd.read_stata(ROOT / "bench/small4.dta"), kw, fake)
    return d, spec, np.array(r["theta"]), r, tau, index, design


class Problem:
    def __init__(self, model, seconds=45):
        self.d, self.spec, self.reference, self.record, self.tau, self.index, self.design = setup(model)
        self.cache = make_cache(self.spec)
        self.k, self.m, self.N = self.spec.n_free, self.spec.n_eq_estimated, self.d.nobs
        self.started, self.seconds = time.perf_counter(), seconds
        self.fun_calls, self.jac_calls, self.invalid_trials = 0, 0, 0
        self.last = self.reference.copy()
        self.trace = []

    def check(self):
        if time.perf_counter() - self.started > self.seconds:
            raise TimeoutError(f"experiment wall-clock bound {self.seconds:g} seconds")

    def mu(self, theta):
        self.check()
        self.fun_calls += 1
        with np.errstate(all="ignore"):
            mu = fitted_shares(theta, self.d, self.spec)
        if not np.isfinite(mu).all():
            self.invalid_trials += 1
        return mu

    def u(self, theta):
        return self.d.shares[:, :self.m] - self.mu(theta)

    def J(self, theta):
        self.check()
        self.jac_calls += 1
        with np.errstate(all="ignore"):
            return jacobian_free(theta, self.d, self.spec, self.cache)

    def Sigma(self, theta):
        u = self.u(theta)
        return u.T @ u / self.N

    def profile(self, theta):
        u = self.u(theta)
        if not np.isfinite(u).all():
            return np.inf, np.full(self.k, np.nan)
        try:
            L = np.linalg.cholesky(u.T @ u / self.N)
        except np.linalg.LinAlgError:
            return np.inf, np.full(self.k, np.nan)
        v = solve_triangular(L.T, solve_triangular(L, u.T, lower=True), lower=False).T
        g = -np.einsum("tmk,tm->k", self.J(theta), v, optimize=True) / self.N
        return float(np.log(np.diag(L)).sum()), g


def whiten(S):
    return solve_triangular(np.linalg.cholesky(S), np.eye(len(S)), lower=True)


def mata_qrsolve(A, b):
    """QR-pivoted solve using the tolerance documented in Stata [M-5] qrsolve.

    LAPACK's pivot selection and rounding need not equal Mata's implementation.
    Small diagonal rows are set to zero as described in the Mata documentation.
    """
    Q, R, piv = qr(A, mode="economic", pivoting=True)
    eta = 1e-13 * float(np.abs(np.diag(R)).sum()) / len(R)
    rhs, z = Q.T @ b, np.zeros(len(b))
    rank = 0
    for i in range(len(b)-1, -1, -1):
        if abs(R[i, i]) > eta:
            z[i] = (rhs[i] - R[i, i+1:] @ z[i+1:]) / R[i, i]
            rank += 1
    x = np.empty_like(z)
    x[piv] = z
    return x, rank


def nlsur_inner(p, theta, S, derivative, eps, maxiter=400):
    """Unweighted-observation MinSSQ logic, including non-strict Q acceptance."""
    P = whiten(S)  # Algebraically equals the ado's SVD whitening for SPD S.
    mu = p.mu(theta)
    u = p.d.shares[:, :p.m] - mu
    if not np.isfinite(u).all():
        raise ValueError("nlsur initial function evaluation is nonfinite")
    oldQ = float(np.sum((u @ P.T)**2))
    max_halves, deficient, total_halves = 0, 0, 0
    last = {}
    for it in range(1, maxiter+1):
        p.check()
        if derivative == "analytic":
            J = p.J(theta)
        else:
            # The actual FunctionDeriv uses abs(b), unlike the help's text.
            p.jac_calls += 1
            h = 4e-7 * (np.abs(theta) + 4e-7)
            J = np.empty((p.N, p.m, p.k))
            for j in range(p.k):
                trial = theta.copy()
                trial[j] += h[j]
                md = p.mu(trial)
                if not np.isfinite(md).all():
                    raise ValueError("nlsur forward difference is nonfinite")
                J[:, :, j] = (md-mu)/h[j]
        Jw = np.matmul(P, J).reshape(-1, p.k)
        uw = (u @ P.T).ravel()
        G = Jw.T @ Jw
        G = (G+G.T)/2
        step, rank = mata_qrsolve(G, Jw.T @ uw)
        deficient += int(rank < p.k)
        f, halvings = 1.0, 0
        while True:
            p.check()
            trial = theta + f*step
            mut = p.mu(trial)
            ut = p.d.shares[:, :p.m] - mut
            Q = float(np.sum((ut @ P.T)**2)) if np.isfinite(ut).all() else np.inf
            if Q <= oldQ:
                break
            f *= 0.5
            halvings += 1
            # Operational watchdog, explicitly not a convergence decision.
            if halvings > 1100:
                raise RuntimeError("step halving exceeded floating-point watchdog")
        param_ok = bool(np.all(f*np.abs(step) <= eps*(np.abs(theta)+1e-3)))
        rss_ok = bool(abs(oldQ-Q) <= eps*(oldQ+1e-3))
        last = dict(param_ok=param_ok, rss_ok=rss_ok, iterations=it,
                    Q=Q, step_fraction=f, min_rank=rank, derivative=derivative)
        theta, u, oldQ = trial, ut, Q
        mu = mut
        p.last = theta.copy()
        max_halves = max(max_halves, halvings)
        total_halves += halvings
        if param_ok and rss_ok:
            break
    last.update(success=param_ok and rss_ok, max_halvings=max_halves,
                total_halvings=total_halves, rank_deficient_iterations=deficient)
    p.trace.append(last)
    return theta, last["success"]


def run_nlsur(p, theta, S, start, derivative, tightened=False):
    eps, covtol = (1e-8, 1e-12) if tightened else (1e-5, 1e-10)
    outer = []
    if start == "zero":
        theta, _ = nlsur_inner(p, theta, np.eye(p.m), derivative, eps)
        S = p.Sigma(theta)
        theta, _ = nlsur_inner(p, theta, S, derivative, eps)
        S = p.Sigma(theta)
    iiter, limit = 2, 200  # Matches the ado's ordering of its counter.
    success, reason = False, "outer_limit"
    while True:
        oldtheta, oldS = theta.copy(), S.copy()
        theta, inner_ok = nlsur_inner(p, theta, S, derivative, eps)
        iiter += 1
        S = p.Sigma(theta)
        param = float(np.max(np.abs(theta-oldtheta)/(np.abs(theta)+1)))
        cov = float(np.max(np.abs(S-oldS)/(np.abs(S)+1)))
        outer.append(dict(parameter_change=param, covariance_change=cov,
                          parameter_pass=param < eps, covariance_pass=cov < covtol,
                          inner_success=inner_ok))
        if param < eps or cov < covtol or iiter > limit:
            success = inner_ok and iiter <= limit
            reason = ("outer_limit" if iiter > limit else
                      "parameter_OR_covariance" if param < eps and cov < covtol else
                      "parameter" if param < eps else "covariance")
            break
    return theta, dict(success=success, status=reason, outer=outer,
                       eps=eps, tau=1e-3, delta=4e-7, ifgnlseps=covtol,
                       max_inner=400, max_outer=limit,
                       earlier_inner_failures=sum(not x["success"] for x in p.trace))


def run_scipy(p, theta, S, method, tightened=False):
    tol, otol = (1e-10, 1e-10) if tightened else (1e-8, 1e-8)
    tri = np.tril_indices(p.m)
    def update(state):
        theta = state[:p.k]
        S = np.zeros((p.m,p.m))
        S[tri] = state[p.k:]
        S += np.tril(S, -1).T
        P = whiten(S)
        def fun(x):
            return (p.u(x) @ P.T).ravel()
        def jac(x):
            return -np.matmul(P, p.J(x)).reshape(-1,p.k)
        opt = least_squares(fun, theta, jac=jac, method=method,
                            loss="linear", x_scale="jac", ftol=tol, xtol=tol,
                            gtol=tol, max_nfev=400)
        p.last = opt.x.copy()
        p.trace.append(dict(success=bool(opt.success), status=int(opt.status),
                            message=opt.message, nfev=int(opt.nfev), njev=int(opt.njev or 0),
                            cost=float(opt.cost), optimality=float(opt.optimality)))
        if not opt.success:
            raise RuntimeError("inner SciPy failure: " + opt.message)
        return np.concatenate([opt.x, p.Sigma(opt.x)[tri]])
    state = fixed_point(update, np.concatenate([theta, S[tri]]),
                        xtol=otol, maxiter=200, method="iteration")
    return state[:p.k], dict(success=True, status="native_fixed_point",
                            ftol=tol, xtol=tol, gtol=tol, outer_xtol=otol,
                            max_outer=200, max_nfev=400)


def profile_scale(p):
    """Fixed coordinate scaling: inverse RMS whitened-Jacobian column norms.

    Computed at the SAME zero/identity-weight anchor for every optimizer/start.
    This changes coordinates, not the objective or a stopping predicate.
    """
    anchor = np.zeros(p.k)
    Jw = p.J(anchor).reshape(-1,p.k)
    norms = np.linalg.norm(Jw, axis=0) / np.sqrt(p.N)
    norms[norms == 0] = 1.0
    return 1/norms


def run_profile(p, theta, method, tightened=False):
    scale = profile_scale(p)
    origin = theta.copy()
    gtol = 1e-8 if tightened else 1e-5
    def fg(z):
        f, g = p.profile(origin + scale*z)
        return f, g*scale
    def callback(z):
        p.last = origin + scale*z
    opts = dict(gtol=gtol, maxiter=1000)
    if method == "BFGS":
        opts["xrtol"] = 0.0
    else:
        # Explicit native defaults, and a tighter native ftol sensitivity run.
        opts.update(ftol=1e-12 if tightened else 2.220446049250313e-9,
                    maxls=20, maxfun=15000, maxcor=10)
    with warnings.catch_warnings(record=True) as caught:
        opt = minimize(fg, np.zeros(p.k), jac=True, method=method,
                       callback=callback, options=opts)
    theta = origin + scale*opt.x
    p.last = theta.copy()
    return theta, dict(success=bool(opt.success), status=int(opt.status),
                        message=str(opt.message), nfev=int(opt.nfev), njev=int(opt.njev),
                        nit=int(opt.nit), optimizer_gradient_inf=float(np.max(np.abs(opt.jac))),
                        options=opts, scaling_anchor="common zero/identity-weight start",
                        warnings=list(dict.fromkeys(str(w.message) for w in caught)))


def metrics(p, theta):
    # Diagnostics do NOT alter the optimizers' success flags.
    u = residuals(theta, p.d, p.spec)
    S = u.T @ u / p.N
    P = whiten(S)
    J = jacobian_free(theta, p.d, p.spec, p.cache)
    Jw = np.matmul(P, J).reshape(-1,p.k)
    g = Jw.T @ (u @ P.T).ravel()
    colnorm = np.linalg.norm(Jw, axis=0)
    score = np.divide(np.abs(g), colnorm*np.sqrt(p.N),
                      out=np.zeros_like(g), where=colnorm != 0)
    logdet = float(np.linalg.slogdet(S)[1])
    ru = residuals(p.reference, p.d, p.spec)
    rS = ru.T @ ru / p.N
    gain = -p.N/2*(logdet-float(np.linalg.slogdet(rS)[1]))
    means = sample_means(p.d, p.index if p.spec.censor else None, p.spec)
    def elas(t):
        return elasticities(unpack(t,p.spec), p.spec, means, p.d.a0,
                            tau=p.tau if p.spec.censor else None,
                            np_prob=p.design.shape[1]+1 if p.spec.censor else None).as_stata_vector()
    return dict(gaussian_ll=-p.N/2*(p.m*(1+np.log(2*np.pi))+logdet),
                ll_gain_vs_recorded=gain, concentrated_loss=logdet/2,
                profile_gradient_inf=float(np.max(np.abs(g))/p.N),
                rms_scaled_score_inf=float(np.max(score)),
                covariance_condition=float(np.linalg.cond(S)),
                theta_diff_vs_recorded=float(np.max(np.abs(theta-p.reference))),
                fitted_diff_vs_recorded=float(np.max(np.abs(u-ru))),
                elasticity_diff_vs_recorded=float(np.max(np.abs(elas(theta)-elas(p.reference)))),
                ray_min=float(np.min(1+p.d.demo @ unpack(theta,p.spec).rho)),
                theta=theta.tolist(), sigma=S.tolist())


def experiment(model, start, method, tightened=False, seconds=45, warm_theta=None):
    p = Problem(model, seconds)
    theta = np.zeros(p.k) if start == "zero" else p.reference.copy()
    if warm_theta is not None:
        theta = np.asarray(warm_theta).copy()
    S = np.eye(p.m) if start == "zero" else p.Sigma(theta)
    p.last = theta.copy()
    row = dict(model=model, start=start, method=method, tightened=tightened,
               reference_package_converged=bool(p.record["converged"]))
    try:
        if method.startswith("scipy-"):
            theta, info = run_scipy(p,theta,S,method.split("-",1)[1],tightened)
        elif method.startswith("profile-"):
            theta, info = run_profile(p,theta,method.split("-",1)[1],tightened)
        else:
            theta, info = run_nlsur(p,theta,S,start,method.split("-",1)[1],tightened)
        row.update(info)
    except Exception as e:
        theta = p.last
        row.update(success=False, status="exception", error=type(e).__name__+": "+str(e))
    row.update(elapsed_seconds=time.perf_counter()-p.started,
               model_calls=p.fun_calls, jacobian_calls=p.jac_calls,
               invalid_trials=p.invalid_trials, inner_calls=len(p.trace), history=p.trace)
    try:
        with np.errstate(all="ignore"):
            row.update(metrics(p, theta))
    except Exception as e:
        row["metric_error"] = type(e).__name__ + ": " + str(e)
        row["theta"] = theta.tolist()
    return row


def gradient_probes():
    """Independent central directional differences across a range of steps."""
    rows = []
    for model in MODELS:
        p = Problem(model)
        for start in ("zero", "recorded"):
            theta = np.zeros(p.k) if start == "zero" else p.reference.copy()
            # Perturb the recorded fit to test derivatives away from stationarity.
            if start == "recorded":
                theta = theta + 1e-4*np.sin(np.arange(p.k))
            f, g = p.profile(theta)
            rng = np.random.default_rng(24681357)
            for direction in range(3):
                v = rng.normal(size=p.k)
                v /= np.linalg.norm(v)
                exact = float(g @ v)
                for h in (1e-5, 1e-6, 1e-7):
                    plus, _ = p.profile(theta+h*v)
                    minus, _ = p.profile(theta-h*v)
                    fd = (plus-minus)/(2*h)
                    rows.append(dict(model=model, point=start, direction=direction,
                                     h=h, analytic=exact, numerical=fd,
                                     absolute_error=abs(fd-exact),
                                     scaled_error=abs(fd-exact)/max(1,abs(exact))))
    return rows


CSV_FIELDS = ("model", "start", "method", "tightened", "success", "status",
              "gaussian_ll", "ll_gain_vs_recorded", "profile_gradient_inf",
              "rms_scaled_score_inf", "theta_diff_vs_recorded", "fitted_diff_vs_recorded",
              "elasticity_diff_vs_recorded", "ray_min", "covariance_condition",
              "elapsed_seconds", "inner_calls", "model_calls", "jacobian_calls",
              "invalid_trials", "error", "metric_error")


def write_csv(path, rows, fields):
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--models", nargs="+", default=list(MODELS), choices=MODELS)
    parser.add_argument("--methods", nargs="+", default=list(METHODS), choices=METHODS)
    parser.add_argument("--starts", nargs="+", default=["zero","recorded"],
                        choices=["zero","recorded"])
    parser.add_argument("--seconds", type=float, default=45)
    parser.add_argument("--skip-sensitivity", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    before = source_hashes()
    rows = []
    with threadpool_limits(limits=1):
        probes = gradient_probes()
        dump(args.output / "gradient-probes.json", probes)
        write_csv(args.output / "gradient-probes.csv", probes, list(probes[0]))
        configs = [(m,s,a,False) for m in args.models for s in args.starts for a in args.methods]
        if not args.skip_sensitivity:
            configs += [(m,s,a,True) for m in args.models if m in ("cq","ca","uq","ua")
                        for s in args.starts for a in args.methods]
        for model,start,method,tight in configs:
            row = experiment(model,start,method,tight,args.seconds)
            rows.append(row)
            dump(args.output / "runs.json", rows)
            write_csv(args.output / "runs.csv", rows, CSV_FIELDS)
            print(json.dumps({k:row.get(k) for k in
                  ("model","start","method","tightened","success","status",
                   "gaussian_ll","rms_scaled_score_inf","elapsed_seconds","error")}), flush=True)
    after = source_hashes()
    assert before == after, "Package source or dataset changed"
    dump(args.output / "environment.json", dict(python=platform.python_version(),
         numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__,
         source_and_data_unchanged=True, source_data_sha256=before,
         scope="Fixed first stage, small4, stage-two point-estimator comparison; no native Stata; no bootstrap or new inference implementation",
         command_arguments=vars(args)|{"output":str(args.output)},
         blas_threads=1, run_count=len(rows)))
    from scipy.optimize._lsq.common import check_termination
    from scipy.optimize._minpack_py import _fixed_point_helper
    (args.output / "scipy-native-termination.txt").write_text(
        inspect.getsource(check_termination)+"\n"+inspect.getsource(_fixed_point_helper))


if __name__ == "__main__":
    main()
