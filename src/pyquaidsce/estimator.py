"""
``quaidsce`` — censored QUAIDS estimation, a Python port of Juan C. Caro's
Stata package (Caro, Melo, Molina & Salgado; version 2.0, June 2025).

Pipeline, mirroring ``quaidsce_c.ado`` step by step
---------------------------------------------------
1.  build the estimation sample (Stata's ``marksample`` / ``markout``);
2.  validate: >= 3 shares, prices > 0, expenditure > 0, shares summing to one
    when ``censor=False``, at least one demographic when ``censor=True``;
3.  **first stage** — one probit per share of the 0/1 participation indicator on
    ``[ln p_1..ln p_n, ln m, z_1..z_R]``; stack the coefficients into ``tau``
    and their covariances block-diagonally into ``setau``;
4.  **second stage** — ``nlsur`` on the Shonkwiler-Yen transformed share
    equations;
5.  delta-method map from the estimated to the reported parameter vector,
    ``Vfull = Delta V Delta'``, with ``setau`` appended block-diagonally;
6.  at-means expenditure / uncompensated / compensated elasticities;
7.  optional nonparametric bootstrap (Stata's ``parallel bs, reps()``).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Sequence, Union

import numpy as np
from scipy.stats import norm

from .elasticities import Means, elasticities, sample_means
from .model import DemandData
from .nlsur import nlsur
from .params import Spec, delta_matrix, full_vector, unpack
from .probit import ProbitResult, probit
from .results import QuaidsceResults

ArrayLike = Union[np.ndarray, Sequence[float]]


# --------------------------------------------------------------------------- #
def _as_matrix(df, cols: Optional[Sequence[str]]) -> np.ndarray:
    if cols is None:
        return np.zeros((len(df), 0))
    return np.column_stack([np.asarray(df[c], dtype=float) for c in cols])


@dataclass
class FirstStage:
    tau: np.ndarray  # (n * np_prob,) stacked probit coefficients
    setau: np.ndarray  # (n*np_prob, n*np_prob) block diagonal
    cdf: np.ndarray  # (N, n)
    pdf: np.ndarray  # (N, n)
    du: np.ndarray  # (N, n) whatever `predict` produced
    np_prob: int
    results: List[ProbitResult]


def first_stage(
    shares: np.ndarray,
    lnp: np.ndarray,
    lnexp: np.ndarray,
    demo: np.ndarray,
    predict: str = "pr",
    include_lnexp: bool = True,
) -> FirstStage:
    """Shonkwiler-Yen step 1: a probit per share.

    ``predict`` selects what goes into ``normal()``/``normalden()``:

    * ``"pr"`` — what the shipped ado does, because Stata's ``predict`` after
      ``probit`` defaults to the predicted probability;
    * ``"xb"`` — the linear predictor, i.e. textbook Shonkwiler-Yen.

    ``include_lnexp=False`` reproduces the behaviour of the ado when the user
    supplies ``lnexpenditure()`` instead of ``expenditure()``: the local macro
    holding the log-expenditure temp variable is empty in that branch, so log
    expenditure silently drops out of the first-stage probits.
    """
    predict = str(predict).lower()
    if predict not in {"pr", "xb"}:
        raise ValueError("predict must be 'pr' or 'xb'")
    N, n = shares.shape
    Z = [lnp]
    if include_lnexp:
        Z.append(lnexp[:, None])
    if demo.shape[1]:
        Z.append(demo)
    X = np.hstack(Z)
    np_prob = X.shape[1] + 1  # + intercept

    tau = np.zeros(n * np_prob)
    setau = np.zeros((n * np_prob, n * np_prob))
    cdf = np.ones((N, n))
    pdf = np.zeros((N, n))
    du = np.ones((N, n))
    res: List[ProbitResult] = []

    for i in range(n):
        w = shares[:, i]
        if w.min() > 0:
            raise ValueError(
                f"no censoring for share {i + 1} found "
                "(quaidsce requires zeros in every share)"
            )
        z = (w > 0).astype(float)
        if z.max() == z.min():
            raise ValueError(
                f"participation indicator for share {i + 1} has no variation"
            )
        pr = probit(z, X, add_constant=True)
        if not pr.converged:
            raise RuntimeError(f"first-stage probit {i + 1} did not converge")
        res.append(pr)
        sl = slice(i * np_prob, (i + 1) * np_prob)
        tau[sl] = pr.b
        setau[sl, sl] = pr.V
        xb = pr.xb(X)
        du[:, i] = norm.cdf(xb) if predict == "pr" else xb
        pdf[:, i] = norm.pdf(du[:, i])
        cdf[:, i] = norm.cdf(du[:, i])

    return FirstStage(tau, setau, cdf, pdf, du, np_prob, res)


# --------------------------------------------------------------------------- #
def quaidsce(
    data,
    shares: Sequence[str],
    *,
    prices: Optional[Sequence[str]] = None,
    lnprices: Optional[Sequence[str]] = None,
    expenditure: Optional[str] = None,
    lnexpenditure: Optional[str] = None,
    demographics: Optional[Sequence[str]] = None,
    anot: float,
    quadratic: bool = True,
    censor: bool = True,
    method: str = "fgnls",
    initial: Optional[ArrayLike] = None,
    sigma_initial: Optional[np.ndarray] = None,
    boot_sigma_tol: float = 1e-7,
    start: str = "zero",
    reps: int = 0,
    seed: Optional[int] = None,
    bootstrap_start: str = "zero",
    first_stage_predict: str = "pr",
    strict_stata: bool = True,
    vce_sigma: str = "objective",
    algorithm: str = "gn",
    tol: float = 1e-13,
    max_outer: int = 200,
    max_iter: int = 300,
    chunk: int = 2000,
    nrtol_stop: float = 1e-12,
    inner_nrtol_early: float = 1e-8,
    sigma_tol: float = 1e-11,
    stop_rule: str = "tight",
    n_jobs: int = 1,
    verbose: bool = True,
    gn_verbose: bool = False,
    log=None,
) -> QuaidsceResults:
    """Estimate a (censored) quadratic almost-ideal demand system.

    Parameters follow the Stata syntax as closely as possible.

    Parameters
    ----------
    data : pandas.DataFrame
    shares : the budget-share variables, **all** of them.
    prices / lnprices : exactly one of the two, in the same order as *shares*.
    expenditure / lnexpenditure : exactly one of the two.
    demographics : Ray (1983) scaling variables. Required when ``censor=True``.
    anot : the ``alpha_0`` of the translog price index (Stata's ``anot()``).
    quadratic : ``False`` reproduces ``noquadratic`` (i.e. plain AIDS).
    censor : ``False`` reproduces ``nocensor`` (i.e. Poi's ``quaids``).
    method : ``"nls"``, ``"fgnls"`` (the package default) or ``"ifgnls"``.
    reps : bootstrap replications; 0 disables the bootstrap.
    first_stage_predict : ``"pr"`` reproduces the shipped Stata code, ``"xb"``
        the textbook Shonkwiler-Yen estimator. See :func:`first_stage`.
    strict_stata : keep the (documented) quirks of the original elasticity code.
    """
    say = log or (print if verbose else (lambda *_: None))

    method = str(method).lower()
    algorithm = str(algorithm).lower()
    start = str(start).lower()
    first_stage_predict = str(first_stage_predict).lower()
    vce_sigma = str(vce_sigma).lower()
    stop_rule = str(stop_rule).lower()
    bootstrap_start = str(bootstrap_start).lower()
    if method not in {"nls", "fgnls", "ifgnls"}:
        raise ValueError("method must be 'nls', 'fgnls', or 'ifgnls'")
    if algorithm not in {"gn", "lm"}:
        raise ValueError("algorithm must be 'gn' or 'lm'")
    if start not in {"zero", "linear"}:
        raise ValueError("start must be 'zero' or 'linear'")
    if first_stage_predict not in {"pr", "xb"}:
        raise ValueError("first_stage_predict must be 'pr' or 'xb'")
    if vce_sigma not in {"objective", "final"}:
        raise ValueError("vce_sigma must be 'objective' or 'final'")
    if stop_rule not in {"tight", "stata"}:
        raise ValueError("stop_rule must be 'tight' or 'stata'")
    if bootstrap_start not in {"zero", "warm"}:
        raise ValueError("bootstrap_start must be 'zero' or 'warm'")
    if not np.isfinite(anot):
        raise ValueError("anot must be finite")
    if int(max_outer) < 2 or int(max_iter) < 1 or int(chunk) < 1:
        raise ValueError("max_outer must be >= 2; max_iter and chunk must be >= 1")
    if tol <= 0 or nrtol_stop <= 0 or sigma_tol <= 0:
        raise ValueError("convergence tolerances must be positive")

    if (prices is None) == (lnprices is None):
        raise ValueError("specify exactly one of prices= or lnprices=")
    if (expenditure is None) == (lnexpenditure is None):
        raise ValueError("specify exactly one of expenditure= or lnexpenditure=")

    shares = list(shares)
    neqn = len(shares)
    demo_names = list(demographics or [])
    spec = Spec(neqn=neqn, ndemo=len(demo_names), quadratic=quadratic,
                censor=censor)

    price_names = list(prices or lnprices)
    if len(price_names) != neqn:
        raise ValueError(
            f"number of price variables must equal number of equations ({neqn})"
        )

    # ---- 1. estimation sample (marksample / markout) ---------------------- #
    need = shares + price_names + demo_names + [expenditure or lnexpenditure]
    M = _as_matrix(data, need)
    touse = np.isfinite(M).all(axis=1)

    W = _as_matrix(data, shares)[touse]
    if W.shape[0] == 0:
        raise ValueError("no complete observations remain in the estimation sample")
    if np.any(W < 0):
        raise ValueError("expenditure shares must be nonnegative")
    if prices is not None:
        P = _as_matrix(data, price_names)[touse]
        if P.min() <= 0:
            bad = price_names[int(np.argmin(P.min(axis=0)))]
            raise ValueError(f"nonpositive value(s) for {bad} found")
        lnp = np.log(P)
    else:
        lnp = _as_matrix(data, price_names)[touse]
    if expenditure is not None:
        m = np.asarray(data[expenditure], dtype=float)[touse]
        if m.min() <= 0:
            raise ValueError(f"nonpositive value(s) for {expenditure} found")
        lnexp = np.log(m)
        exp_name = expenditure
    else:
        lnexp = np.asarray(data[lnexpenditure], dtype=float)[touse]
        exp_name = lnexpenditure
    Z = _as_matrix(data, demo_names)[touse] if demo_names else np.zeros(
        (int(touse.sum()), 0)
    )
    N = W.shape[0]

    if not censor:
        s = W.sum(axis=1)
        if np.max(np.abs(s - 1.0) / np.maximum(np.abs(s), 1.0)) >= 1e-4:
            raise ValueError("expenditure shares do not sum to one")

    # ---- 2/3. first stage ------------------------------------------------- #
    notes: List[str] = []
    if censor:
        say("Estimating first-stage probits...")
        fs = first_stage(
            W, lnp, lnexp, Z, predict=first_stage_predict,
            include_lnexp=(expenditure is not None),
        )
        if first_stage_predict == "pr":
            notes.append(
                "First stage uses cdf=Phi(Phi(x'tau)), pdf=phi(Phi(x'tau)) "
                "because Stata's `predict` after `probit` defaults to the "
                "predicted probability. This reproduces quaidsce v2.0 exactly; "
                "pass first_stage_predict='xb' for the textbook "
                "Shonkwiler-Yen transformation."
            )
        if expenditure is None:
            notes.append(
                "lnexpenditure() was used, so - exactly as in quaidsce_c.ado - "
                "log expenditure is omitted from the first-stage probits."
            )
    else:
        fs = FirstStage(
            tau=np.zeros(0), setau=np.zeros((0, 0)),
            cdf=np.ones((N, neqn)), pdf=np.zeros((N, neqn)),
            du=np.ones((N, neqn)), np_prob=0, results=[],
        )

    d = DemandData(lnp=lnp, lnexp=lnexp, shares=W, demo=Z,
                   cdf=fs.cdf, pdf=fs.pdf, a0=float(anot))

    # ---- 4. second stage -------------------------------------------------- #
    theta0 = None if initial is None else np.asarray(initial, float).ravel()
    if theta0 is not None:
        if theta0.size != spec.n_free:
            raise ValueError(
                f"initial must contain {spec.n_free} free parameters; "
                f"got {theta0.size}"
            )
        if not np.isfinite(theta0).all():
            raise ValueError("initial must contain only finite values")
    if sigma_initial is not None:
        sigma_initial = np.asarray(sigma_initial, float)
        want = (spec.n_eq_estimated, spec.n_eq_estimated)
        if sigma_initial.shape != want:
            raise ValueError(f"sigma_initial must have shape {want}")
        if not np.isfinite(sigma_initial).all():
            raise ValueError("sigma_initial must contain only finite values")
    nl = nlsur(
        d, spec, theta0=theta0, sigma0=sigma_initial, start=start,
        method=method, tol=tol,
        max_outer=max_outer,
        max_iter=max_iter, chunk=chunk, nrtol_stop=nrtol_stop, sigma_tol=sigma_tol,
        inner_nrtol_early=inner_nrtol_early, stop_rule=stop_rule,
        vce_sigma=vce_sigma,
        algorithm=algorithm,
        verbose=False, log=say, gn_log=(print if gn_verbose else None),
    )

    # ---- 5. delta method -------------------------------------------------- #
    Delta = delta_matrix(spec)
    bfull = full_vector(nl.theta, spec)
    Vfull = Delta @ nl.V @ Delta.T
    coefs = unpack(nl.theta, spec)

    if censor:
        b = np.concatenate([bfull, fs.tau])
        k1, k2 = Vfull.shape[0], fs.setau.shape[0]
        V = np.zeros((k1 + k2, k1 + k2))
        V[:k1, :k1] = Vfull
        V[k1:, k1:] = fs.setau
    else:
        b, V = bfull, Vfull

    names = spec.full_names(demo_names)
    if censor:
        names = names + spec.tau_names(demo_names)

    # ---- 6. elasticities -------------------------------------------------- #
    means = sample_means(d, fs.du if censor else None, spec)
    el = elasticities(
        coefs, spec, means, a0=float(anot),
        tau=fs.tau if censor else None,
        np_prob=fs.np_prob if censor else None,
        strict_stata=strict_stata,
    )
    if censor:
        ev = el.as_stata_vector()
        b = np.concatenate([b, ev])
        k1, k2 = V.shape[0], ev.size
        Vx = np.zeros((k1 + k2, k1 + k2))
        Vx[:k1, :k1] = V
        V = Vx
        names = names + spec.elas_names()

    if strict_stata and spec.ndemo == 0 and spec.quadratic:
        notes.append(
            "With no demographics, quaidsce_c.ado's uncompensated elasticity "
            "uses beta_i*lambda_i in the last term where Poi (2012) has "
            "beta_j*lambda_i; reproduced here. Pass strict_stata=False to use "
            "the published formula."
        )
    if spec.ndemo > 0 and not spec.quadratic and censor:
        notes.append(
            "quaidsce_c.ado stores the expenditure elasticity in a *global* "
            "macro in the demographics + noquadratic branch and then reads an "
            "empty local. Stata therefore silently uses a zero latent expenditure "
            "elasticity. strict_stata=True reproduces that result; "
            "strict_stata=False uses the intended formula 1 + betanz_i/w_i."
        )
    if not nl.converged:
        notes.append(
            "The nonlinear estimator did not satisfy all requested convergence "
            "criteria. Inspect n_outer/n_gn and refit with larger iteration limits "
            "or different starting values before using the estimates."
        )

    res = QuaidsceResults(
        spec=spec, anot=float(anot), nobs=N, method=method,
        share_names=shares, price_names=price_names, demo_names=demo_names,
        expenditure_name=exp_name,
        coefs=coefs, theta=nl.theta, b=b, V=V, names=names,
        b_est=nl.theta, V_est=nl.V,
        llf=nl.llf, sigma=nl.sigma, elas=el, means=means,
        tau=fs.tau if censor else None,
        setau=fs.setau if censor else None,
        np_prob=fs.np_prob, probits=list(fs.results),
        n_outer=nl.n_outer, n_gn=nl.n_gn, converged=nl.converged,
        notes=notes,
    )

    # ---- 7. bootstrap ----------------------------------------------------- #
    if reps and reps > 0:
        from .bootstrap import bootstrap

        res.boot = bootstrap(
            data=data, shares=shares, prices=prices, lnprices=lnprices,
            expenditure=expenditure, lnexpenditure=lnexpenditure,
            demographics=demographics, anot=anot, quadratic=quadratic,
            censor=censor, method=method, initial=nl.theta,
            sigma_initial=nl.sigma,
            first_stage_predict=first_stage_predict, strict_stata=strict_stata,
            vce_sigma=vce_sigma, sigma_tol=boot_sigma_tol,
            algorithm=algorithm, stop_rule=stop_rule, tol=tol,
            max_outer=max_outer, max_iter=max_iter, chunk=chunk,
            nrtol_stop=nrtol_stop, inner_nrtol_early=inner_nrtol_early,
            bootstrap_start=bootstrap_start,
            reps=int(reps), seed=seed, n_jobs=n_jobs,
            touse=touse, verbose=verbose,
        )
    return res
