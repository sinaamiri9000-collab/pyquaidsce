"""
Nonlinear seemingly-unrelated-regression estimator, replicating Stata's
``nlsur`` for the three estimators the package can use.

Given residuals ``u_t(theta)`` (an ``m``-vector per observation):

* **nls**    minimise ``sum_t u_t' u_t``                       (Sigma = I)
* **fgnls**  NLS, then ``Sigma_hat = (1/N) sum_t u_t u_t'``, then minimise
             ``sum_t u_t' Sigma_hat^-1 u_t``                   (two-step)
* **ifgnls** iterate the FGNLS step until the parameters stop moving; this is
             the maximum-likelihood estimator under joint normality

Reported quantities, matching ``e()`` after ``nlsur``:

    e(V)  = ( sum_t X_t' Sigma_hat^-1 X_t )^-1 ,  X_t = d f(x_t, theta)/d theta'
    e(ll) = -(N m / 2) (1 + ln 2*pi) - (N/2) ln |Sigma_hat|

``Sigma_hat`` uses divisor ``N`` (Stata's default; ``dfk`` is not implemented
because ``quaidsce`` never requests it).

The minimiser is a Levenberg-Marquardt-damped Gauss-Newton method using the
exact analytic Jacobian, accumulating the normal equations in observation
blocks so that the ``(N*m) x K`` Jacobian is never materialised in full.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Tuple

import numpy as np

from ._timing import check_deadline
from .jacfree import JacCache, jacobian_free, make_cache
from .model import DemandData, jacobian_full, residuals
from .params import Spec, delta_blocks


LOG2PI = float(np.log(2.0 * np.pi))

try:  # halves the flops of the normal-equations accumulation
    from scipy.linalg.blas import dsyrk as _DSYRK
except Exception:  # pragma: no cover
    _DSYRK = None


# --------------------------------------------------------------------------- #
@dataclass
class NlsurResult:
    theta: np.ndarray
    V: np.ndarray
    sigma: np.ndarray
    llf: float
    obj: float
    nobs: int
    neq: int
    method: str
    n_outer: int
    n_gn: int
    converged: bool
    history: List[float] = field(default_factory=list)


# --------------------------------------------------------------------------- #
def _whitener(sigma: np.ndarray) -> np.ndarray:
    """P with ``P' P = Sigma^-1`` (so ``u' Sigma^-1 u = ||P u||^2``)."""
    L = np.linalg.cholesky(sigma)  # sigma = L L'
    # Sigma^-1 = L^-T L^-1  =>  P = L^-1
    return np.linalg.inv(L)


def _jac_free_reference(
    theta: np.ndarray,
    d: DemandData,
    spec: Spec,
    blocks,
    rows: slice,
) -> np.ndarray:
    """Reference route: differentiate w.r.t. the full vector, then apply Delta.

    Kept only so that the fast path in :mod:`pyquaidsce.jacfree` can be checked
    against it; the estimator itself never calls this.
    """
    Jf = jacobian_full(theta, d, spec, rows=rows)
    nc, m, _ = Jf.shape
    out = np.zeros((nc, m, spec.n_free))
    for full_sl, free_sl, mat in blocks:
        out[:, :, free_sl] = Jf[:, :, full_sl] @ mat
    return out


def _normal_equations(
    theta: np.ndarray,
    d: DemandData,
    spec: Spec,
    cache: JacCache,
    P: np.ndarray,
    chunk: int,
    deadline: Optional[float] = None,
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Accumulate ``G = sum J'Sinv J``, ``g = sum J'Sinv u`` and the objective."""
    N = d.nobs
    m = spec.n_eq_estimated
    K = spec.n_free
    G = np.zeros((K, K))
    g = np.zeros(K)
    obj = 0.0
    identity_P = bool(np.allclose(P, np.eye(m), atol=0.0, rtol=0.0))
    for s in range(0, N, chunk):
        check_deadline(deadline)
        e = min(s + chunk, N)
        sl = slice(s, e)
        u = d.shares[sl, :m] - _fitted_chunk(theta, d, spec, sl)
        J = jacobian_free(theta, d, spec, cache, sl)
        # whiten: batched (m x m) @ (m x K) via BLAS, then one symmetric rank-k
        if identity_P:  # the NLS stage needs no whitening at all
            uw = u
            Jw2 = J.reshape(-1, K)
        else:
            uw = u @ P.T  # (nc, m)
            Jw2 = np.ascontiguousarray(np.matmul(P, J).reshape(-1, K))
        if _DSYRK is not None:
            G += _DSYRK(1.0, Jw2, trans=1, lower=0)
        else:
            G += Jw2.T @ Jw2
        g += Jw2.T @ uw.reshape(-1)
        obj += float(uw.ravel() @ uw.ravel())
    if _DSYRK is not None:  # dsyrk fills only the upper triangle
        G = np.triu(G) + np.triu(G, 1).T
    return G, g, obj


def _fitted_chunk(theta, d: DemandData, spec: Spec, sl: slice) -> np.ndarray:
    from .model import fitted_shares

    sub = DemandData(
        lnp=d.lnp[sl], lnexp=d.lnexp[sl], shares=d.shares[sl],
        demo=d.demo[sl], cdf=d.cdf[sl], pdf=d.pdf[sl], a0=d.a0,
        control_function=d.control_function[sl],
    )
    return fitted_shares(theta, sub, spec)


def _objective(theta, d: DemandData, spec: Spec, P: np.ndarray) -> float:
    u = residuals(theta, d, spec)
    uw = u if P is None else u @ P.T
    return float(uw.ravel() @ uw.ravel())


# --------------------------------------------------------------------------- #
def _safe_obj(
    theta,
    d: DemandData,
    spec: Spec,
    P: np.ndarray,
    deadline: Optional[float] = None,
) -> float:
    check_deadline(deadline)
    try:
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            v = _objective(theta, d, spec, P)
        return v if np.isfinite(v) else np.inf
    except (FloatingPointError, ValueError):
        return np.inf


def _solve_scaled(G: np.ndarray, g: np.ndarray, mu: float = 0.0,
                  rcond: float = 1e-14) -> np.ndarray:
    """Solve ``(G + mu diag(G)) d = g`` in diagonally scaled coordinates.

    The parameter blocks of a QUAIDS system differ in scale by several orders of
    magnitude (think ``rho`` versus ``gamma``), so scaling by ``diag(G)^-1/2``
    before the solve is what keeps the Gauss-Newton direction meaningful.  A
    truncated pseudo-inverse then absorbs the near-singular directions instead
    of letting them blow the step up.
    """
    dg = np.diag(G).copy()
    dg[dg <= 0] = 1.0
    sc = 1.0 / np.sqrt(dg)
    Gs = G * sc[:, None] * sc[None, :]
    if mu:
        Gs = Gs + mu * np.eye(Gs.shape[0])
    gs = g * sc
    Gs = 0.5 * (Gs + Gs.T)
    try:
        ev, Q = np.linalg.eigh(Gs)
        keep = ev > rcond * max(ev.max(), 1e-300)
        ds = Q[:, keep] @ ((Q[:, keep].T @ gs) / ev[keep])
    except np.linalg.LinAlgError:  # pragma: no cover
        ds = np.linalg.lstsq(Gs, gs, rcond=None)[0]
    return ds * sc


def gauss_newton(
    theta0: np.ndarray,
    d: DemandData,
    spec: Spec,
    sigma: np.ndarray,
    tol: float = 1e-13,
    max_iter: int = 200,
    chunk: int = 2000,
    verbose: bool = False,
    say=None,
    algorithm: str = "gn",
    nrtol_stop: float = 1e-12,
    stop_rule: str = "tight",
    deadline: Optional[float] = None,
) -> Tuple[np.ndarray, float, int, bool]:
    """Minimise ``sum_t u_t' sigma^-1 u_t``.

    ``algorithm="lm"`` is a Levenberg-Marquardt trust-region method
    with Nielsen's damping update, working in coordinates scaled by
    ``diag(J'J)^-1/2``.  ``algorithm="gn"`` (the default) is plain
    Gauss-Newton with Hartley
    step halving, i.e. the scheme Stata's ``nl``/``nlsur`` use; it is kept as the
    reference behaviour, but on a large censored system it creeps along the flat
    valley of the criterion while LM turns the corner.  Both target the same
    stationary point.
    """
    cache = make_cache(spec)
    P = _whitener(sigma)
    theta = np.asarray(theta0, dtype=float).copy()
    converged = False
    it = 0
    mu = 1e-3
    nu = 2.0

    for it in range(1, max_iter + 1):
        check_deadline(deadline)
        G, g, obj = _normal_equations(
            theta, d, spec, cache, P, chunk, deadline=deadline
        )

        if algorithm == "gn":
            direction = _solve_scaled(G, g)
            nrtol = abs(float(direction @ g)) / max(abs(obj), 1e-300)
            best, t = None, 1.0
            for _ in range(40):
                check_deadline(deadline)
                cand = theta + t * direction
                oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                if oc < obj:
                    best = (cand, oc, t, 0.0)
                    break
                t *= 0.5
            if best is None:
                m2 = 1e-8
                for _ in range(30):
                    check_deadline(deadline)
                    cand = theta + _solve_scaled(G, g, mu=m2)
                    oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                    if oc < obj:
                        best = (cand, oc, 1.0, m2)
                        break
                    m2 *= 10.0
            if best is None:
                cutoff = 1e-5 if stop_rule == "stata" else nrtol_stop
                converged = bool(nrtol < cutoff)
                break
            cand, oc, t, m2 = best
        else:
            gn = _solve_scaled(G, g)
            nrtol = abs(float(gn @ g)) / max(abs(obj), 1e-300)
            accepted = False
            for _ in range(50):
                check_deadline(deadline)
                step = _solve_scaled(G, g, mu=mu)
                cand = theta + step
                oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                pred = 2.0 * float(step @ g) - float(step @ (G @ step))
                rho = (obj - oc) / pred if pred > 0 else -1.0
                if oc < obj and rho > 0:
                    accepted = True
                    if rho > 0.75:
                        mu = max(mu / 3.0, 1e-14)
                        nu = 2.0
                    elif rho < 0.25:
                        mu = min(mu * 2.0, 1e12)
                    break
                mu = min(mu * nu, 1e14)
                nu *= 2.0
            if not accepted:
                cutoff = 1e-5 if stop_rule == "stata" else nrtol_stop
                converged = bool(nrtol < cutoff)
                break
            t, m2 = 1.0, mu

        rel = (obj - oc) / max(abs(obj), 1e-300)
        smax = float(np.max(np.abs(cand - theta) / (np.abs(theta) + 1e-6)))
        # Stata's mreldif(): the denominator is |b_old| + 1, not |b_old|, which
        # makes the criterion far looser for small coefficients.
        mreldif = float(np.max(np.abs(cand - theta) / (np.abs(theta) + 1.0)))
        theta, obj = cand, oc
        if say is not None:
            say(f"    GN {it:3d}  obj={obj:.12g}  rel={rel:.3e} "
                f"step={smax:.2e} t={t:g} mu={m2:.2e} nrtol={nrtol:.2e}")
        # Stationarity is judged by the scale-free relative gradient (Stata's
        # nrtolerance, whose default is a much looser 1e-5).  The parameter-step
        # test alone never fires on a censored system: the criterion has nearly
        # flat directions, so theta keeps drifting long after the objective has
        # stopped moving.
        if stop_rule == "stata":
            # Stata's nl/nlsur declare convergence when ANY of three criteria is
            # met: tolerance() on the coefficient vector, ltolerance() on the
            # objective, or nrtolerance() on the scaled gradient.  Because the
            # test is a disjunction, the coefficient criterion routinely fires
            # while the gradient is still O(1e-4) -- the criterion of a censored
            # demand system has very flat directions.  That is why Stata's
            # -method(nls)- and -method(fgnls)- results are not fully converged.
            if (mreldif < 1e-5) or (rel < 1e-7) or (nrtol < 1e-5):
                converged = True
                break
        elif nrtol < nrtol_stop or (rel < tol and smax < 1e-9):
            converged = True
            break
    return theta, obj, it, converged


# --------------------------------------------------------------------------- #
def nlsur(
    d: DemandData,
    spec: Spec,
    theta0: Optional[np.ndarray] = None,
    sigma0: Optional[np.ndarray] = None,
    start: str = "zero",
    method: str = "fgnls",
    tol: float = 1e-13,
    sigma_tol: float = 1e-11,
    max_outer: int = 200,
    max_iter: int = 200,
    chunk: int = 2000,
    nrtol_stop: float = 1e-12,
    inner_nrtol_early: float = 1e-8,
    stop_rule: str = "tight",
    vce_sigma: str = "objective",
    algorithm: str = "gn",
    verbose: bool = False,
    log: Optional[Callable[[str], None]] = None,
    gn_log: Optional[Callable[[str], None]] = None,
    deadline: Optional[float] = None,
) -> NlsurResult:
    """Fit the system.

    Parameters
    ----------
    method : {"nls", "fgnls", "ifgnls"}
        ``quaidsce`` v2.0 uses ``fgnls`` unless ``method()`` is given; the
        2021 runs archived in ``log/`` were produced with ``ifgnls``.
    vce_sigma : {"objective", "final"}
        Which ``Sigma_hat`` enters ``e(V)``.  ``"objective"`` is the one used in
        the last minimisation (textbook FGNLS, and what Stata reports);
        ``"final"`` recomputes it from the final residuals.  The two coincide
        for ``ifgnls`` at convergence.
    """
    method = method.lower()
    if method not in ("nls", "fgnls", "ifgnls"):
        raise ValueError(f"unknown method {method!r}")
    if algorithm not in ("gn", "lm"):
        raise ValueError(f"unknown algorithm {algorithm!r}")
    if start not in ("zero", "linear"):
        raise ValueError(f"unknown start {start!r}")
    if stop_rule not in ("tight", "stata"):
        raise ValueError(f"unknown stop_rule {stop_rule!r}")
    if vce_sigma not in ("objective", "final"):
        raise ValueError(f"unknown vce_sigma {vce_sigma!r}")
    say = log or (print if verbose else (lambda *_: None))
    check_deadline(deadline)

    N = d.nobs
    m = spec.n_eq_estimated
    K = spec.n_free
    if theta0 is not None:
        theta = np.asarray(theta0, float).copy()
    elif start == "linear":
        from .start import linear_start
        theta = linear_start(d, spec)
    else:
        theta = np.zeros(K)

    def sigma_from(th: np.ndarray) -> np.ndarray:
        u = residuals(th, d, spec)
        return (u.T @ u) / N

    history: List[float] = []
    total_gn = 0

    # ---- warm start: skip straight to the Sigma fixed-point iteration ----- #
    # A bootstrap replication differs from the full sample only by resampling,
    # so its Sigma_hat is already almost the full-sample one.  Handing both
    # theta and Sigma over means the outer loop starts *at* the fixed point and
    # converges in a handful of iterations instead of rediscovering it from the
    # NLS stage -- which is where nearly all of the cost of IFGNLS sits.
    if sigma0 is not None and theta0 is not None and method == "ifgnls":
        sigma_obj = np.asarray(sigma0, dtype=float)
        history.append(_objective(theta, d, spec, _whitener(sigma_obj)))
        n_outer = 2
        tight = False
        outer_converged = False
        ok = False
        obj = history[-1]
        for outer in range(3, max_outer + 1):
            check_deadline(deadline)
            say(f"FGNLS iteration {outer}...")
            sigma_new = sigma_from(theta)
            theta_prev = theta.copy()
            theta, obj, ngn, ok = gauss_newton(
                theta, d, spec, sigma_new, tol=tol, max_iter=max_iter,
                chunk=chunk, say=gn_log, algorithm=algorithm,
                nrtol_stop=(nrtol_stop if tight
                            else max(nrtol_stop, inner_nrtol_early)),
                stop_rule=stop_rule,
                deadline=deadline,
            )
            total_gn += ngn
            history.append(obj)
            sigma_obj = sigma_new
            n_outer = outer
            rel = float(np.max(np.abs(theta - theta_prev)
                               / (np.abs(theta_prev) + 1e-8)))
            if rel < sigma_tol:
                if tight:
                    outer_converged = True
                    break
                tight = True
            else:
                tight = False
        return _finish(theta, d, spec, sigma_obj, sigma_from, method,
                       vce_sigma, chunk, n_outer, total_gn,
                       bool(ok and outer_converged), obj, history,
                       deadline=deadline)

    # ---- step 1: NLS ------------------------------------------------------ #
    say("Calculating NLS estimates...")
    I_m = np.eye(m)
    theta, obj, ngn, ok = gauss_newton(
        theta, d, spec, I_m, tol=tol, max_iter=max_iter, chunk=chunk,
        say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop,
        stop_rule=stop_rule,
        deadline=deadline,
    )
    total_gn += ngn
    history.append(obj)
    sigma_obj = I_m
    n_outer = 1

    if method != "nls":
        # ---- step 2: FGNLS ----------------------------------------------- #
        say("Calculating FGNLS estimates...")
        sigma = sigma_from(theta)
        theta_prev = theta.copy()
        theta, obj, ngn, ok = gauss_newton(
            theta, d, spec, sigma, tol=tol, max_iter=max_iter, chunk=chunk,
            say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop,
            stop_rule=stop_rule,
            deadline=deadline,
        )
        total_gn += ngn
        history.append(obj)
        sigma_obj = sigma
        n_outer = 2

        if method == "ifgnls":
            # Inexact-outer strategy.  IFGNLS is a fixed-point iteration on
            # Sigma_hat, so there is nothing to gain from solving the early
            # inner problems to machine precision -- their solutions are thrown
            # away by the next Sigma update anyway.  We therefore solve the
            # inner problems loosely (`inner_nrtol_early`) while the outer
            # iteration is still moving, and only tighten to `nrtol_stop` once
            # the outer loop has essentially converged, confirming convergence
            # with a final fully-tight pass.  The fixed point is unchanged; the
            # number of Gauss-Newton steps drops by a large factor.
            tight = False
            outer_converged = False
            for outer in range(3, max_outer + 1):
                check_deadline(deadline)
                say(f"FGNLS iteration {outer}...")
                sigma_new = sigma_from(theta)
                theta_prev = theta.copy()
                theta, obj, ngn, ok = gauss_newton(
                    theta, d, spec, sigma_new, tol=tol, max_iter=max_iter,
                    chunk=chunk, say=gn_log, algorithm=algorithm,
                    nrtol_stop=(nrtol_stop if tight else
                                max(nrtol_stop, inner_nrtol_early)),
                    stop_rule=stop_rule,
                    deadline=deadline,
                )
                total_gn += ngn
                history.append(obj)
                sigma_obj = sigma_new
                n_outer = outer
                rel = float(
                    np.max(np.abs(theta - theta_prev) / (np.abs(theta_prev) + 1e-8))
                )
                if rel < sigma_tol:
                    if tight:
                        outer_converged = True
                        break
                    tight = True  # one more pass, solved to full precision
                else:
                    tight = False

    overall_converged = bool(ok)
    if method == "ifgnls":
        overall_converged = bool(ok and outer_converged)
    return _finish(theta, d, spec, sigma_obj, sigma_from, method, vce_sigma,
                   chunk, n_outer, total_gn, overall_converged, obj, history,
                   deadline=deadline)


def _finish(theta, d, spec, sigma_obj, sigma_from, method, vce_sigma, chunk,
            n_outer, total_gn, ok, obj, history,
            deadline: Optional[float] = None) -> NlsurResult:
    """Assemble e(V), e(ll) and the result object."""
    check_deadline(deadline)
    N, m = d.nobs, spec.n_eq_estimated
    sigma_final = sigma_from(theta)
    sigma_V = sigma_final if vce_sigma == "final" else sigma_obj
    if method == "nls":
        sigma_V = sigma_final

    cache = make_cache(spec)
    P = _whitener(sigma_V)
    G, _, _ = _normal_equations(
        theta, d, spec, cache, P, chunk, deadline=deadline
    )
    V = np.linalg.inv(G)

    _, logdet = np.linalg.slogdet(sigma_final)
    llf = -(N * m / 2.0) * (1.0 + LOG2PI) - (N / 2.0) * logdet

    return NlsurResult(
        theta=theta, V=V, sigma=sigma_final, llf=float(llf), obj=float(obj),
        nobs=N, neq=m, method=method, n_outer=n_outer, n_gn=total_gn,
        converged=bool(ok), history=history,
    )
