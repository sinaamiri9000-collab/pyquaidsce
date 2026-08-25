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

from ._timing import check_deadline
from .elasticities import Means, elasticities, sample_means
from .model import DemandData, fitted_shares
from .nlsur import nlsur
from .params import Spec, delta_matrix, full_vector, unpack
from .probit import ProbitResult, probit
from .reduced_form import fit_expenditure_reduced_form
from .results import QuaidsceResults
from .selection import FirstStageLayout, legacy_layout

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
    layout: FirstStageLayout


def first_stage(
    shares: np.ndarray,
    lnp: np.ndarray,
    lnexp: np.ndarray,
    demo: np.ndarray,
    predict: str = "pr",
    include_lnexp: bool = True,
    *,
    design: Optional[np.ndarray] = None,
    layout: Optional[FirstStageLayout] = None,
    deadline: Optional[float] = None,
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
    check_deadline(deadline)
    if predict not in {"pr", "xb"}:
        raise ValueError("predict must be 'pr' or 'xb'")
    N, n = shares.shape
    if design is None:
        Z = [lnp]
        if include_lnexp:
            Z.append(lnexp[:, None])
        if demo.shape[1]:
            Z.append(demo)
        X = np.hstack(Z) if Z else np.zeros((N, 0))
        if layout is None:
            layout = legacy_layout(
                [f"p{j + 1}" for j in range(lnp.shape[1])],
                [f"z{r + 1}" for r in range(demo.shape[1])],
                include_expenditure=include_lnexp,
            )
    else:
        X = np.asarray(design, dtype=float)
        if X.ndim == 1:
            X = X[:, None]
        if X.shape[0] != N:
            raise ValueError("first-stage design must have one row per observation")
        if layout is None:
            raise ValueError("an explicit first-stage design requires a layout")
    assert layout is not None
    if X.shape[1] != layout.constant_position:
        raise ValueError("first-stage design and layout have inconsistent widths")
    np_prob = layout.width

    tau = np.zeros(n * np_prob)
    setau = np.zeros((n * np_prob, n * np_prob))
    cdf = np.ones((N, n))
    pdf = np.zeros((N, n))
    du = np.ones((N, n))
    res: List[ProbitResult] = []

    for i in range(n):
        check_deadline(deadline)
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
        pr = probit(z, X, add_constant=True, deadline=deadline)
        if not pr.converged:
            raise RuntimeError(f"first-stage probit {i + 1} did not converge")
        if (layout.selection_cf_position is not None
                and layout.selection_cf_position in pr.dropped):
            raise ValueError(
                "selection_control_function is collinear with the first-stage "
                f"design in equation {i + 1}; its coefficient was dropped"
            )
        res.append(pr)
        sl = slice(i * np_prob, (i + 1) * np_prob)
        tau[sl] = pr.b
        setau[sl, sl] = pr.V
        xb = pr.xb(X)
        du[:, i] = norm.cdf(xb) if predict == "pr" else xb
        pdf[:, i] = norm.pdf(du[:, i])
        cdf[:, i] = norm.cdf(du[:, i])

    return FirstStage(tau, setau, cdf, pdf, du, np_prob, res, layout)


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
    ivexp: Optional[Sequence[str]] = None,
    control_function: Optional[str] = None,
    selection_prices: Optional[Sequence[str]] = None,
    selection_expenditure: Optional[bool] = True,
    selection_covariates: Optional[Sequence[str]] = None,
    selection_control_function: Optional[str] = None,
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
    first_stage_predict: str = "xb",
    strict_stata: bool = False,
    vce_sigma: str = "objective",
    algorithm: str = "gn",
    tol: float = 1e-13,
    max_outer: int = 200,
    max_iter: int = 300,
    chunk: int = 2000,
    nrtol_stop: float = 1e-12,
    inner_nrtol_early: float = 1e-8,
    sigma_tol: float = 1e-11,
    stop_rule: str = "standard",
    n_jobs: int = 1,
    verbose: bool = True,
    gn_verbose: bool = False,
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
    log=None,
    _deadline: Optional[float] = None,
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
    control_function : column containing an externally generated reduced-form
        residual. It enters the latent share as ``cfcoef_i * residual``.
    ivexp : excluded instrument column(s) for endogenous log expenditure. The
        package estimates ``ln(m)`` on log prices, Ray demographics, these
        instruments, and a constant. Its residual automatically enters both
        the participation Probits and latent demand equations. It cannot be
        combined with either external control-function argument.
    selection_* : configure the first-stage Probit independently of the Ray
        demographics. ``None`` for prices/covariates preserves their legacy
        sets; ``[]`` selects none. The default includes log expenditure.
    mp_context : multiprocessing start method for bootstrap workers. The safe
        cross-platform default is ``"spawn"``.
    rep_timeout : optional per-replication wall-clock limit in seconds.
        Cooperative Probit/optimizer checks are backed by a parent-side process
        watchdog that can terminate a stuck native call.
    """
    say = log or (print if verbose else (lambda *_: None))
    check_deadline(_deadline)

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
    if stop_rule not in {"tight", "standard"}:
        raise ValueError("stop_rule must be 'tight' or 'standard'")
    if bootstrap_start not in {"zero", "warm"}:
        raise ValueError("bootstrap_start must be 'zero' or 'warm'")
    if not np.isfinite(anot):
        raise ValueError("anot must be finite")
    if int(max_outer) < 2 or int(max_iter) < 1 or int(chunk) < 1:
        raise ValueError("max_outer must be >= 2; max_iter and chunk must be >= 1")
    if tol <= 0 or nrtol_stop <= 0 or sigma_tol <= 0:
        raise ValueError("convergence tolerances must be positive")
    if rep_timeout is not None:
        if not np.isfinite(rep_timeout) or float(rep_timeout) <= 0:
            raise ValueError("rep_timeout must be a finite positive number")

    if (prices is None) == (lnprices is None):
        raise ValueError("specify exactly one of prices= or lnprices=")
    if (expenditure is None) == (lnexpenditure is None):
        raise ValueError("specify exactly one of expenditure= or lnexpenditure=")
    if selection_expenditure not in {None, True, False}:
        raise ValueError("selection_expenditure must be True, False, or None")

    if isinstance(ivexp, str):
        raise ValueError("ivexp must be a sequence of column names, not one string")
    ivexp_names = [] if ivexp is None else list(ivexp)
    if ivexp is not None and not ivexp_names:
        raise ValueError("ivexp must contain at least one excluded instrument")
    if len(ivexp_names) != len(set(ivexp_names)):
        raise ValueError("ivexp must not contain duplicate column names")
    ivexp_active = bool(ivexp_names)
    if ivexp_active and (
        control_function is not None or selection_control_function is not None
    ):
        raise ValueError(
            "ivexp cannot be combined with control_function or "
            "selection_control_function"
        )
    external_cf_active = (
        control_function is not None or selection_control_function is not None
    )
    cf_active = external_cf_active or ivexp_active
    selection_custom = (
        selection_prices is not None
        or selection_covariates is not None
        or selection_expenditure is False
    )
    extension_active = cf_active or selection_custom
    if extension_active and not censor:
        raise ValueError("control-function/selection extensions require censor=True")
    if extension_active and first_stage_predict != "xb":
        raise ValueError(
            "control-function/selection extensions require "
            "first_stage_predict='xb'"
        )
    if external_cf_active and reps and int(reps) > 0:
        raise ValueError(
            "bootstrap with a precomputed control function is disabled: the "
            "reduced form residual must be rebuilt inside every replication; "
            "use ivexp for the internally rebuilt expenditure control function"
        )

    shares = list(shares)
    neqn = len(shares)
    demo_names = list(demographics or [])
    for label, values in (("shares", shares), ("demographics", demo_names)):
        if len(values) != len(set(values)):
            raise ValueError(f"{label} must not contain duplicate column names")
    spec = Spec(
        neqn=neqn,
        ndemo=len(demo_names),
        quadratic=quadratic,
        censor=censor,
        control_function=(control_function is not None or ivexp_active),
    )

    price_names = list(prices or lnprices)
    if len(price_names) != len(set(price_names)):
        raise ValueError("prices/lnprices must not contain duplicate column names")
    if len(price_names) != neqn:
        raise ValueError(
            f"number of price variables must equal number of equations ({neqn})"
        )
    instrument_overlap = sorted(
        set(ivexp_names)
        & set(shares + price_names + demo_names + [expenditure or lnexpenditure])
    )
    if instrument_overlap:
        raise ValueError(
            "ivexp variables must be excluded instruments, not shares, prices, "
            "demographics, or expenditure; overlap: "
            + ", ".join(instrument_overlap)
        )

    selected_prices = price_names if selection_prices is None else list(selection_prices)
    if len(selected_prices) != len(set(selected_prices)):
        raise ValueError("selection_prices must not contain duplicates")
    invalid_selection_prices = [p for p in selected_prices if p not in price_names]
    if invalid_selection_prices:
        raise ValueError(
            "selection_prices must be a subset of demand prices; invalid: "
            + ", ".join(invalid_selection_prices)
        )
    selected_covariates = (
        demo_names if selection_covariates is None else list(selection_covariates)
    )
    if len(selected_covariates) != len(set(selected_covariates)):
        raise ValueError("selection_covariates must not contain duplicates")
    invalid_iv_selection = sorted(set(ivexp_names) & set(selected_covariates))
    if invalid_iv_selection:
        raise ValueError(
            "ivexp variables are excluded instruments and must not enter "
            "selection_covariates; overlap: " + ", ".join(invalid_iv_selection)
        )
    for label, name in (
        ("control_function", control_function),
        ("selection_control_function", selection_control_function),
    ):
        if name is not None and name in demo_names:
            raise ValueError(
                f"{label} must not also be a Ray demographic; use its dedicated API"
            )
    if (control_function is not None
            and control_function in selected_covariates
            and selection_control_function != control_function):
        raise ValueError(
            "a demand control residual in the Probit must be supplied through "
            "selection_control_function, not selection_covariates"
        )
    include_selection_expenditure = (
        expenditure is not None
        if selection_expenditure is None
        else bool(selection_expenditure)
    )

    selection_sources = list(selected_prices)
    if include_selection_expenditure:
        selection_sources.append(expenditure or lnexpenditure)
    selection_sources.extend(selected_covariates)
    if selection_control_function is not None:
        selection_sources.append(selection_control_function)
    if len(selection_sources) != len(set(selection_sources)):
        raise ValueError("the first-stage design must not contain duplicate columns")

    # ---- 1. estimation sample (marksample / markout) ---------------------- #
    need = (
        shares + price_names + demo_names + [expenditure or lnexpenditure]
        + ivexp_names
        + selected_covariates
        + ([control_function] if control_function is not None else [])
        + ([selection_control_function]
           if selection_control_function is not None else [])
    )
    need = list(dict.fromkeys(need))
    missing = [name for name in need if name not in data]
    if missing:
        raise ValueError("column(s) not found in data: " + ", ".join(missing))
    for name in dict.fromkeys(
        x for x in (control_function, selection_control_function) if x is not None
    ):
        raw_residual = np.asarray(data[name], dtype=float)
        if not np.isfinite(raw_residual).all():
            raise ValueError(
                f"control-function column {name!r} must contain only finite values"
            )
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
    reduced_form = None
    if ivexp_active:
        instruments = _as_matrix(data, ivexp_names)[touse]
        reduced_form = fit_expenditure_reduced_form(
            lnexp,
            lnp,
            Z,
            instruments,
            outcome_name=(
                f"ln({expenditure})" if expenditure is not None
                else str(lnexpenditure)
            ),
            price_names=price_names,
            price_inputs_are_logs=lnprices is not None,
            demographic_names=demo_names,
            instrument_names=ivexp_names,
        )
        demand_cf = reduced_form.residuals.copy()
        selection_cf = reduced_form.residuals.copy()
    else:
        demand_cf = (
            np.asarray(data[control_function], dtype=float)[touse]
            if control_function is not None else np.zeros(int(touse.sum()))
        )
        selection_cf = (
            np.asarray(data[selection_control_function], dtype=float)[touse]
            if selection_control_function is not None else None
        )
    for label, values in (
        ("control_function", demand_cf if control_function is not None else None),
        ("selection_control_function", selection_cf),
    ):
        if values is not None:
            scale = max(1.0, float(np.max(np.abs(values))))
            if float(np.ptp(values)) <= 1e-12 * scale:
                raise ValueError(f"{label} residual must have nonzero variation")
    N = W.shape[0]

    if not censor:
        s = W.sum(axis=1)
        if np.max(np.abs(s - 1.0) / np.maximum(np.abs(s), 1.0)) >= 1e-4:
            raise ValueError("expenditure shares do not sum to one")

    # ---- 2/3. first stage ------------------------------------------------- #
    notes: List[str] = []
    if censor:
        check_deadline(_deadline)
        say("Estimating first-stage probits...")
        design_parts = []
        ordered_names: List[str] = []
        price_positions = {}
        for price in selected_prices:
            j = price_names.index(price)
            price_positions[price] = len(ordered_names)
            ordered_names.append(f"p{j + 1}")
            design_parts.append(lnp[:, j])
        expenditure_position = None
        if include_selection_expenditure:
            expenditure_position = len(ordered_names)
            ordered_names.append("M")
            design_parts.append(lnexp)
        selection_Z = (
            _as_matrix(data, selected_covariates)[touse]
            if selected_covariates else np.zeros((N, 0))
        )
        covariate_positions = {}
        for r, name in enumerate(selected_covariates):
            covariate_positions[name] = len(ordered_names)
            # Preserve legacy labels unless the raw name would collide with a
            # synthetic p1..pn/M/constant label.
            report_name = str(name)
            if report_name in ordered_names or report_name == "cons":
                report_name = f"z[{name}]"
            ordered_names.append(report_name)
            design_parts.append(selection_Z[:, r])
        selection_cf_position = None
        if selection_cf is not None:
            selection_cf_position = len(ordered_names)
            if ivexp_active:
                report_name = "cf_ivexp"
                suffix = 2
                while report_name in ordered_names or report_name == "cons":
                    report_name = f"cf_ivexp_{suffix}"
                    suffix += 1
            else:
                report_name = str(selection_control_function)
                if report_name in ordered_names or report_name == "cons":
                    report_name = f"cf[{selection_control_function}]"
            ordered_names.append(report_name)
            design_parts.append(selection_cf)
        layout = FirstStageLayout(
            ordered_names=tuple(ordered_names),
            demand_price_names=tuple(price_names),
            price_positions=price_positions,
            expenditure_position=expenditure_position,
            covariate_positions=covariate_positions,
            selection_cf_position=selection_cf_position,
            constant_position=len(ordered_names),
        )
        selection_design = (
            np.column_stack(design_parts) if design_parts else np.zeros((N, 0))
        )
        fs = first_stage(
            W, lnp, lnexp, Z, predict=first_stage_predict,
            include_lnexp=include_selection_expenditure,
            design=selection_design, layout=layout,
            deadline=_deadline,
        )
        if first_stage_predict == "pr":
            notes.append(
                "First stage uses cdf=Phi(Phi(x'tau)), pdf=phi(Phi(x'tau)) "
                "because Stata's `predict` after `probit` defaults to the "
                "predicted probability. This reproduces quaidsce v2.0 exactly; "
                "pass first_stage_predict='xb' for the textbook "
                "Shonkwiler-Yen transformation."
            )
        if not include_selection_expenditure:
            notes.append(
                "Log expenditure is omitted from the first-stage probits."
            )
    else:
        fs = FirstStage(
            tau=np.zeros(0), setau=np.zeros((0, 0)),
            cdf=np.ones((N, neqn)), pdf=np.zeros((N, neqn)),
            du=np.ones((N, neqn)), np_prob=0, results=[],
            layout=legacy_layout([], [], include_expenditure=False),
        )

    d = DemandData(lnp=lnp, lnexp=lnexp, shares=W, demo=Z,
                   cdf=fs.cdf, pdf=fs.pdf, a0=float(anot),
                   control_function=demand_cf)

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
        deadline=_deadline,
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
        names = names + fs.layout.tau_names(neqn)

    # ---- 6. elasticities -------------------------------------------------- #
    means = sample_means(d, fs.du if censor else None, spec)
    el = elasticities(
        coefs, spec, means, a0=float(anot),
        tau=fs.tau if censor else None,
        np_prob=fs.np_prob if censor else None,
        layout=fs.layout if censor else None,
        strict_stata=strict_stata,
    )
    if censor:
        ev = el.as_stata_vector()
        b = np.concatenate([b, ev])
        k1, k2 = V.shape[0], ev.size
        # The legacy/Stata-compatible result vector contains elasticities, but
        # this release does not claim an analytical delta-method covariance for
        # them. Keep e(V) finite/usable and mask only analytic_se below.
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
    if spec.ndemo > 0 and not spec.quadratic and censor and not cf_active:
        notes.append(
            "quaidsce_c.ado stores the expenditure elasticity in a *global* "
            "macro in the demographics + noquadratic branch and then reads an "
            "empty local. Stata therefore silently uses a zero latent expenditure "
            "elasticity. strict_stata=True reproduces that result; "
            "strict_stata=False uses the intended formula 1 + betanz_i/w_i."
        )
    if spec.ndemo > 0 and not spec.quadratic and censor and cf_active:
        notes.append(
            "The published noquadratic expenditure-elasticity formula is used "
            "for the control-function extension even with strict_stata=True. "
            "Reproducing the legacy Stata local/global-macro bug would be "
            "inconsistent with the derivative of the augmented fitted share."
        )
    if not nl.converged:
        notes.append(
            "The nonlinear estimator did not satisfy all requested convergence "
            "criteria. Inspect n_outer/n_gn and refit with larger iteration limits "
            "or different starting values before using the estimates."
        )
    if ivexp_active:
        notes.append(
            "ivexp internally estimates the log-expenditure reduced form on "
            "log prices, Ray demographics, excluded instruments, and a "
            "constant. Its residual enters both the participation Probits and "
            "latent demand equations. Structural elasticities hold that "
            "residual fixed."
        )
        if not (reps and int(reps) > 0):
            notes.append(
                "Analytical covariance is conditional on the generated ivexp "
                "residual. Use the internal bootstrap for generated-regressor "
                "inference; it re-estimates the reduced form in every replication."
            )
    elif external_cf_active:
        notes.append(
            "Control-function covariance and p-values are conditional on the "
            "externally generated residual; final inference must rebuild the "
            "reduced form in a design-appropriate bootstrap."
        )
    if censor and not (reps and int(reps) > 0):
        notes.append(
            "Analytical elasticity standard errors are not computed. The "
            "reported structural/first-stage analytical covariance is a "
            "block-diagonal conditional approximation; use a full bootstrap "
            "for generated-regressor inference."
        )
    fitted = fitted_shares(nl.theta, d, spec)
    if np.any(fitted < 0):
        notes.append(
            f"Diagnostic: {int(np.sum(fitted < 0))} fitted share values are "
            "negative; values were not clipped."
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
        selection_layout=fs.layout if censor else None,
        control_function_name=("ivexp" if ivexp_active else control_function),
        selection_control_function_name=(
            "ivexp" if ivexp_active else selection_control_function
        ),
        ivexp_names=ivexp_names,
        reduced_form=reduced_form,
        n_outer=nl.n_outer, n_gn=nl.n_gn, converged=nl.converged,
        notes=notes,
    )

    # ---- 7. bootstrap ----------------------------------------------------- #
    if reps and reps > 0:
        from .bootstrap import bootstrap

        res.boot = bootstrap(
            data=data, shares=shares, prices=prices, lnprices=lnprices,
            expenditure=expenditure, lnexpenditure=lnexpenditure,
            demographics=demographics, ivexp=ivexp_names if ivexp_active else None,
            control_function=control_function,
            selection_control_function=selection_control_function,
            selection_prices=selection_prices,
            selection_covariates=selection_covariates,
            selection_expenditure=selection_expenditure,
            anot=anot, quadratic=quadratic,
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
            mp_context=mp_context, rep_timeout=rep_timeout,
        )
        res.V_analytic = res.V.copy()
        res.V = res.boot.V.copy()
        res.notes.append(
            "res.V and res.se contain the successful-replication bootstrap "
            "covariance; res.V_analytic and res.analytic_se retain the "
            "conditional analytical reference."
        )
    return res
