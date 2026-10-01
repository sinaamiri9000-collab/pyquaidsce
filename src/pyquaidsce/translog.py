"""Translog with optional S&Y, using shared participation, SUR and bootstrap."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence

import numpy as np

from ._threads import with_blas_threads
from ._timing import check_deadline
from .data import DemandData, as_matrix as _matrix
from .censoring import normal_terms
from .elasticities import Means
from .inference import complex_step_jacobian
from .nlsur import nlsur
from .result_base import DemandResults
from .participation import first_stage
from .selection import (FirstStageLayout, build_selection_design,
                        selection_layout, selection_index)
from .statafmt import coef_table, g
from .translog_model import (TRANSLOG_BACKEND, TranslogCoefs, TranslogSpec,
                            elasticities, full_vector, latent_shares,
                            predicted_shares, pack, unpack, _inner)


def _input_arrays(data, price_names, expenditure_name, demo_names,
                  prices_are_logs, expenditure_is_log):
    names = price_names + [expenditure_name] + demo_names
    missing = [c for c in names if c not in data]
    if missing:
        raise ValueError("column(s) not found in data: " + ", ".join(missing))
    p = _matrix(data, price_names)
    m = np.asarray(data[expenditure_name], dtype=float)
    z = _matrix(data, demo_names)
    if not prices_are_logs and np.any(p[np.isfinite(p)] <= 0):
        raise ValueError("prices must be strictly positive")
    if not expenditure_is_log and np.any(m[np.isfinite(m)] <= 0):
        raise ValueError("expenditure must be strictly positive")
    with np.errstate(invalid="ignore", divide="ignore"):
        lnp = p if prices_are_logs else np.log(p)
        lnm = m if expenditure_is_log else np.log(m)
    return lnp, lnm, z


def _demand_data(lnp, lnm, z, shares=None):
    shape = lnp.shape
    if shares is None:
        shares = np.zeros(shape)
    return DemandData(lnp, lnm, shares, z, np.ones(shape), np.zeros(shape))


def _point_at_means(data, price_names, expenditure_name, demo_names,
                    prices_are_logs, expenditure_is_log):
    lnp, lnm, z = _input_arrays(data, price_names, expenditure_name, demo_names,
                               prices_are_logs, expenditure_is_log)
    # Stata's atmeans is based on arithmetic means of levels, even when the
    # equations internally use logarithms. Each input group uses its own
    # complete rows: price vectors, demographic vectors, and expenditure.
    # On the estimation sample these masks coincide. An explicit frame can
    # reproduce Stata postestimation over a larger sample with missing prices.
    p = np.exp(lnp)
    m = np.exp(lnm)
    price_rows = np.isfinite(p).all(axis=1)
    demo_rows = np.isfinite(z).all(axis=1)
    if not price_rows.any() or not np.isfinite(m).any():
        raise ValueError("no finite observations for elasticity evaluation")
    if z.shape[1] and not demo_rows.any():
        raise ValueError("no finite demographic observations for elasticity evaluation")
    m = np.where(np.isfinite(m), m, np.nan)
    return _demand_data(np.log(p[price_rows].mean(axis=0))[None],
                        np.array([np.log(np.nanmean(m))]),
                        z[demo_rows].mean(axis=0)[None] if z.shape[1] else np.empty((1, 0)))


def _elasticities(theta, spec, point, tau=None, layout=None, design=None):
    if not spec.censor:
        return elasticities(theta, point, spec)
    k = selection_index(design, tau, layout, spec.neqn)
    return elasticities(theta, point, spec, tau=tau, layout=layout, selection_index=k)


def _reported(theta, spec, point, tau=None, layout=None, design=None):
    parts = [full_vector(theta, spec)]
    if spec.censor:
        parts.append(tau)
    parts.append(_elasticities(theta, spec, point, tau, layout, design).as_stata_vector())
    return np.concatenate(parts)


def _selection_at_means(frame, point, layout, covariate_names):
    cov = _matrix(frame, covariate_names)
    rows = np.isfinite(cov).all(axis=1)
    if cov.shape[1] and not rows.any():
        raise ValueError("no finite selection covariates for elasticity evaluation")
    means = cov[rows].mean(axis=0)[None] if cov.shape[1] else np.empty((1, 0))
    return build_selection_design(point.lnp, point.lnexp, means, layout)


@dataclass
class ElasticityStandardErrors:
    income: np.ndarray
    uncompensated: np.ndarray
    compensated: np.ndarray


@dataclass
class TranslogResults(DemandResults):
    spec: TranslogSpec
    nobs: int
    method: str
    share_names: list
    price_names: list
    demo_names: list
    expenditure_name: str
    coefs: TranslogCoefs
    theta: np.ndarray
    b: np.ndarray
    V: np.ndarray
    names: list
    b_est: np.ndarray
    V_est: np.ndarray
    llf: float
    sigma: np.ndarray
    elas: object
    means: Means
    r2: np.ndarray
    estimation_mask: np.ndarray
    data: object = field(repr=False)
    demand_data: DemandData = field(repr=False)
    evaluation_point: DemandData = field(repr=False)
    prices_are_logs: bool = False
    expenditure_is_log: bool = False
    n_outer: int = 0
    n_gn: int = 0
    converged: bool = True
    boot: Optional[object] = None
    V_analytic: Optional[np.ndarray] = None
    notes: list = field(default_factory=list)
    tau: Optional[np.ndarray] = None
    setau: Optional[np.ndarray] = None
    np_prob: int = 0
    probits: list = field(default_factory=list)
    selection_layout: Optional[FirstStageLayout] = None
    selection_covariate_names: list = field(default_factory=list)
    selection_index: Optional[np.ndarray] = field(default=None, repr=False)
    evaluation_design: Optional[np.ndarray] = field(default=None, repr=False)
    start_method: str = "zero"
    algorithm: str = "gn"
    stop_rule: str = "standard"
    nrtol: float = np.nan
    initialization: dict = field(default_factory=dict)

    @property
    def n_coefficients(self):
        return self.spec.n_full + (0 if self.tau is None else self.tau.size)

    @property
    def compensated_label(self):
        return ("Compensated (Slutsky transformation) price elasticities" if self.spec.censor
                else "Compensated (Hicksian) price elasticities")

    @property
    def analytic_se(self):
        covariance = self.V if self.V_analytic is None else self.V_analytic
        diagonal = np.diag(covariance).copy()
        diagonal[diagonal < 0] = np.nan
        return np.sqrt(diagonal)

    @property
    def se(self):
        return self.boot.se.copy() if self.boot is not None else self.analytic_se

    @property
    def coefficient_covariance(self):
        return self.V[:self.n_coefficients, :self.n_coefficients].copy()

    @property
    def elasticity_covariance(self):
        return self.V[self.n_coefficients:, self.n_coefficients:].copy()

    @property
    def elasticity_se(self):
        n = self.spec.neqn
        se = self.se[self.n_coefficients:]
        return ElasticityStandardErrors(se[:n], se[n:n+n*n].reshape(n, n),
                                        se[n+n*n:].reshape(n, n))

    def _title(self):
        return "Basic translog model" + (" with Shonkwiler-Yen correction" if self.spec.censor else "")

    def summary(self, level=95.0, elasticities=False):
        count = len(self.b) if elasticities else self.n_coefficients
        head = [self._title(), f"Number of obs          = {self.nobs}",
                f"Number of goods        = {self.spec.neqn}",
                f"Demographic method     = Translating",
                f"Number of demographics = {self.spec.ndemo}",
                f"Latent adding-up       = {self.spec.latent_adding_up}",
                f"Starting values        = {self.start_method}",
                f"Scaled gradient        = {g(self.nrtol, 12)}",
                f"Log-likelihood         = {g(self.llf, 12)}",
                f"Converged              = {self.converged}"]
        if self.boot is not None:
            head.append(f"Bootstrap replications = {self.boot.reps_ok}/{self.boot.reps_requested}")
        table = coef_table(self.names[:count], self.b[:count], self.se[:count],
                           level=level, bootstrap=self.boot is not None)
        text = "\n".join(head) + "\n\n" + table
        if self.notes:
            text += "\n\nNotes:\n" + "\n".join(f"- {x}" for x in self.notes)
        return text

    def predict(self, data=None, kind="shares"):
        """Predict every good; preserve supplied row order and missing-input rows."""
        if kind not in {"shares", "quantities", "residuals", "latent_shares", "participation"}:
            raise ValueError("invalid prediction kind")
        if data is None:
            d = self.demand_data
            if kind == "participation":
                return d.cdf.copy()
            fitted = (latent_shares if kind == "latent_shares" else predicted_shares)(self.theta, d, self.spec)
            if kind == "residuals":
                return d.shares - fitted
            return fitted if kind in {"shares", "latent_shares"} else fitted / np.exp(d.lnp-d.lnexp[:, None])
        lnp, lnm, z = _input_arrays(data, self.price_names, self.expenditure_name,
                                    self.demo_names, self.prices_are_logs,
                                    self.expenditure_is_log)
        need_selection = self.spec.censor and kind != "latent_shares"
        cov = _matrix(data, self.selection_covariate_names) if need_selection else np.empty((len(data), 0))
        valid = np.isfinite(np.column_stack([lnp, lnm, z, cov])).all(axis=1)
        out = np.full((len(data), self.spec.neqn), np.nan)
        if valid.any():
            d = _demand_data(lnp[valid], lnm[valid], z[valid])
            if need_selection:
                X = build_selection_design(d.lnp, d.lnexp, cov[valid], self.selection_layout)
                k = selection_index(X, self.tau, self.selection_layout, self.spec.neqn)
                d.cdf, d.pdf = normal_terms(k)
            fitted = (d.cdf if kind == "participation" else
                      (latent_shares if kind == "latent_shares" else predicted_shares)(self.theta, d, self.spec))
            if kind == "quantities":
                fitted = fitted / np.exp(d.lnp-d.lnexp[:, None])
            out[valid] = fitted
        return _matrix(data, self.share_names) - out if kind == "residuals" else out

    def elasticities_at_means(self, data=None, standard_errors=False):
        """Evaluate at level means on the estimation sample or an explicit frame.

        This makes the sample used for Stata postestimation comparisons explicit.
        Delta-method inference treats evaluation means as fixed.
        """
        point = self.evaluation_point if data is None else _point_at_means(
            data, self.price_names, self.expenditure_name, self.demo_names,
            self.prices_are_logs, self.expenditure_is_log)
        design = self.evaluation_design
        if self.spec.censor and data is not None:
            design = _selection_at_means(data, point, self.selection_layout, self.selection_covariate_names)
        result = _elasticities(self.theta, self.spec, point, self.tau, self.selection_layout, design)
        if not standard_errors:
            return result
        if self.spec.censor:
            if self.boot is not None:
                joint = np.r_[self.theta, self.tau]
                J = complex_step_jacobian(
                    lambda t: _elasticities(t[:self.spec.n_free], self.spec, point,
                        t[self.spec.n_free:], self.selection_layout, design).as_stata_vector(), joint)
                free_cov = self.coefficient_covariance
            else:
                J = complex_step_jacobian(lambda t: _elasticities(t, self.spec, point,
                    self.tau, self.selection_layout, design).as_stata_vector(), self.theta)
                free_cov = self.V_est
        else:
            J = complex_step_jacobian(lambda t: elasticities(t, point, self.spec).as_stata_vector(), self.theta)
            free = np.r_[np.arange(self.spec.neqn-1), np.arange(self.spec.neqn, self.spec.n_full)]
            free_cov = self.V[np.ix_(free, free)]
        se = np.sqrt(np.maximum(np.diag(J @ free_cov @ J.T), 0))
        n = self.spec.neqn
        return result, ElasticityStandardErrors(se[:n], se[n:n+n*n].reshape(n, n),
                                                 se[n+n*n:].reshape(n, n))

    def __str__(self):
        return self.summary()


@with_blas_threads
def translog(
    data, shares: Sequence[str], *, prices: Optional[Sequence[str]] = None,
    lnprices: Optional[Sequence[str]] = None, expenditure: Optional[str] = None,
    lnexpenditure: Optional[str] = None, demographics: Optional[Sequence[str]] = None,
    censor: bool = False, selection_prices: Optional[Sequence[str]] = None,
    selection_expenditure: bool = True, selection_covariates: Optional[Sequence[str]] = None,
    method="ifgnls", initial=None, sigma_initial=None,
    start="zero", algorithm="gn", stop_rule="standard", vce_sigma="objective",
    tol=1e-13, max_outer=200, max_iter=300, chunk=1000, nrtol_stop=1e-12,
    inner_nrtol_early=1e-8, sigma_tol=1e-5, share_sum_tol=0.01,
    reps=0, seed=None, bootstrap_start="zero", n_jobs=1, blas_threads=1,
    mp_context=None, rep_timeout=None, verbose=True, gn_verbose=False, log=None,
    _deadline=None,
) -> TranslogResults:
    """Estimate basic translog demand with demographic translation.

    Uncensored SUR omits the last equation and imposes sum(alpha)=1. S&Y uses
    all equations and frees all alphas, retaining denominator constant 1.
    Probits default to all log prices, log expenditure and the demographics;
    their design can be selected independently. No expenditure control function
    or generalized translog intercept is included.

    ``start='nested'`` fits a no-translation S&Y submodel with the same Probits
    and SUR settings to initialize the full model. No final-model restriction
    or optimizer changes. Explicit initial values take precedence.

    Elasticities use predicted shares at arithmetic means of price/expenditure
    levels on the estimation sample. Analytic standard errors use the delta
    method, conditional on those means and, for S&Y, fitted Probits. ``reps>=2``
    re-estimates all stages and supplies full two-step bootstrap inference.
    """
    check_deadline(_deadline)
    method, algorithm, start = (str(x).lower() for x in (method, algorithm, start))
    stop_rule, vce_sigma, bootstrap_start = (
        str(x).lower() for x in (stop_rule, vce_sigma, bootstrap_start))
    if (prices is None) == (lnprices is None):
        raise ValueError("specify exactly one of prices= or lnprices=")
    if (expenditure is None) == (lnexpenditure is None):
        raise ValueError("specify exactly one of expenditure= or lnexpenditure=")
    demo_values = [] if demographics is None else demographics
    for label, values in (("shares", shares), ("prices", prices if prices is not None else lnprices),
                           ("demographics", demo_values)):
        if isinstance(values, str) or len(values) != len(set(values)):
            raise ValueError(f"{label} must be a sequence of distinct column names")
    shares = list(shares)
    price_names = list(prices if prices is not None else lnprices)
    demo_names = list(demo_values)
    exp_name = expenditure if expenditure is not None else lnexpenditure
    spec = TranslogSpec(len(shares), len(demo_names), bool(censor))
    if len(price_names) != spec.neqn:
        raise ValueError("number of prices must match number of shares")
    if not isinstance(censor, (bool, np.bool_)) or not isinstance(selection_expenditure, (bool, np.bool_)):
        raise ValueError("censor and selection_expenditure must be booleans")
    if not censor and (selection_prices is not None or selection_covariates is not None or not selection_expenditure):
        raise ValueError("selection options require censor=True")
    if start == "nested" and not censor:
        raise ValueError("nested starting values require censor=True")
    if start not in {"zero", "mean", "nested"}:
        raise ValueError("unknown start; allowed values are zero, mean and nested")
    selected_prices = price_names if selection_prices is None else selection_prices
    selected_covariates = demo_names if selection_covariates is None else selection_covariates
    for label, cols in (("selection_prices", selected_prices), ("selection_covariates", selected_covariates)):
        if isinstance(cols, str) or len(cols) != len(set(cols)):
            raise ValueError(f"{label} must be a sequence of distinct column names")
    selected_prices, selected_covariates = list(selected_prices), list(selected_covariates)
    if any(p not in price_names for p in selected_prices):
        raise ValueError("selection_prices must be a subset of demand prices")
    if censor and set(selected_covariates) & set(price_names+[exp_name]+shares):
        raise ValueError("selection_covariates must not repeat prices, expenditure or shares")
    missing = [name for name in selected_covariates if name not in data]
    if missing:
        raise ValueError("column(s) not found in data: " + ", ".join(missing))
    for name in shares:
        if name not in data:
            raise ValueError(f"column not found in data: {name}")
    for label, value in (("max_outer", max_outer), ("max_iter", max_iter),
                           ("chunk", chunk), ("n_jobs", n_jobs), ("reps", reps)):
        if isinstance(value, bool) or not np.isfinite(value) or int(value) != value or value < (0 if label == "reps" else 1):
            raise ValueError(f"{label} must be an integer in its valid range")
    if max_outer < 2 or reps == 1:
        raise ValueError("max_outer must be >= 2 and reps must be 0 or >= 2")
    for label, value in (("tol", tol), ("nrtol_stop", nrtol_stop),
                           ("inner_nrtol_early", inner_nrtol_early),
                           ("sigma_tol", sigma_tol), ("share_sum_tol", share_sum_tol)):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be finite and positive")
    if rep_timeout is not None and (not np.isfinite(rep_timeout) or rep_timeout <= 0):
        raise ValueError("rep_timeout must be finite and positive")
    if bootstrap_start not in {"zero", "warm"}:
        raise ValueError("bootstrap_start must be 'zero' or 'warm'")
    lnp, lnm, z = _input_arrays(data, price_names, exp_name, demo_names,
                                prices is None, expenditure is None)
    w = _matrix(data, shares)
    selection_cov = _matrix(data, selected_covariates) if censor else np.empty((len(data), 0))
    mask = np.isfinite(np.column_stack([w, lnp, lnm, z, selection_cov])).all(axis=1)
    if mask.sum() <= spec.n_free:
        raise ValueError("estimation requires more complete observations than free parameters")
    w = w[mask]
    if np.any(w < 0) or np.any(w > 1):
        raise ValueError("shares must be between zero and one; clean negative values explicitly")
    if not censor and np.max(np.abs(w.sum(axis=1)-1)) > share_sum_tol:
        raise ValueError("expenditure shares do not sum to one within share_sum_tol")
    d = _demand_data(lnp[mask], lnm[mask], z[mask], w)
    if initial is not None:
        initial = np.asarray(initial, dtype=float)
        if initial.shape != (spec.n_free,) or not np.isfinite(initial).all():
            raise ValueError(f"initial must contain {spec.n_free} finite free parameters")
    if sigma_initial is not None:
        sigma_initial = np.asarray(sigma_initial, dtype=float)
        m = spec.n_eq_estimated
        if sigma_initial.shape != (m, m) or not np.isfinite(sigma_initial).all() or not np.allclose(sigma_initial, sigma_initial.T):
            raise ValueError("sigma_initial must be a finite symmetric equation covariance")
        np.linalg.cholesky(sigma_initial)
    say = log or (print if verbose else lambda *_: None)
    fs, layout = None, None
    if censor:
        layout = selection_layout(price_names, selected_prices, selected_covariates,
                                   include_expenditure=selection_expenditure)
        design = build_selection_design(d.lnp, d.lnexp, selection_cov[mask], layout)
        say("Estimating first-stage probits...")
        fs = first_stage(w, d.lnp, d.lnexp, d.demo, design=design, layout=layout, deadline=_deadline)
        d.cdf, d.pdf = fs.cdf, fs.pdf
    initialization = {}
    if start == "nested" and initial is None and spec.ndemo:
        say("Obtaining starting values from the S&Y translog without demographic translation...")
        nested_spec = TranslogSpec(spec.neqn, 0, True)
        nested_data = DemandData(d.lnp, d.lnexp, d.shares, np.empty((d.nobs, 0)), d.cdf, d.pdf)
        nested = nlsur(nested_data, nested_spec, start="zero", method=method,
            algorithm=algorithm, stop_rule=stop_rule, vce_sigma=vce_sigma,
            tol=tol, max_outer=int(max_outer), max_iter=int(max_iter), chunk=int(chunk),
            nrtol_stop=nrtol_stop, inner_nrtol_early=inner_nrtol_early, sigma_tol=sigma_tol,
            blas_threads=blas_threads, log=say, gn_log=say if gn_verbose else None,
            deadline=_deadline, model=TRANSLOG_BACKEND)
        c = unpack(nested.theta, nested_spec)
        initial = pack(TranslogCoefs(c.alpha, c.gamma, np.zeros((spec.neqn, spec.ndemo)), c.delta), spec)
        if sigma_initial is None:
            sigma_initial = nested.sigma
        initialization = dict(method=nested.method, n_outer=nested.n_outer,
                              n_gn=nested.n_gn, converged=nested.converged, nrtol=nested.nrtol)
    nl = nlsur(d, spec, theta0=initial, sigma0=sigma_initial, start="zero" if start == "nested" else start,
               method=method, algorithm=algorithm, stop_rule=stop_rule,
               vce_sigma=vce_sigma, tol=tol, max_outer=int(max_outer),
               max_iter=int(max_iter), chunk=int(chunk), nrtol_stop=nrtol_stop,
               inner_nrtol_early=inner_nrtol_early, sigma_tol=sigma_tol,
               blas_threads=blas_threads, log=say,
               gn_log=say if gn_verbose else None, deadline=_deadline,
               model=TRANSLOG_BACKEND)
    sample = data.loc[mask].copy() if hasattr(data, "loc") else data[mask]
    point = _point_at_means(sample, price_names, exp_name, demo_names,
                             prices is None, expenditure is None)
    tau = None if fs is None else fs.tau
    point_design = None if fs is None else _selection_at_means(sample, point, layout, selected_covariates)
    point_index = np.ones((1, spec.neqn))
    if fs is not None:
        point_index = selection_index(point_design, tau, layout, spec.neqn)
        point.cdf, point.pdf = normal_terms(point_index)
    values = _reported(nl.theta, spec, point, tau, layout, point_design)
    # Conditional structural/elasticity inference. First-stage coefficients are
    # reported separately; their uncertainty is included jointly by bootstrap.
    J = complex_step_jacobian(lambda t: _reported(t, spec, point, tau, layout, point_design), nl.theta)
    covariance = J @ nl.V @ J.T
    if fs is not None:
        sl = slice(spec.n_full, spec.n_full+tau.size)
        covariance[sl, sl] = fs.setau
    covariance = (covariance + covariance.T)/2
    el = _elasticities(nl.theta, spec, point, tau, layout, point_design)
    fitted = predicted_shares(nl.theta, d, spec)
    tss = np.sum((w-w.mean(axis=0))**2, axis=0)
    r2 = np.divide(np.sum((w-fitted)**2, axis=0), tss,
                    out=np.full(spec.neqn, np.nan), where=tss > 0)
    means = Means(el.we.copy(), point.lnp[0].copy(), float(point.lnexp[0]),
                   point.demo[0].copy(), point.cdf[0].copy(), point.pdf[0].copy(),
                   point_index[0].copy())
    names = spec.full_names(demo_names)+(layout.tau_names(spec.neqn) if censor else [])+spec.elas_names()
    res = TranslogResults(
        spec, d.nobs, nl.method, shares, price_names, demo_names, exp_name,
        unpack(nl.theta, spec), nl.theta, values, covariance,
        names, nl.theta.copy(), nl.V,
        nl.llf, nl.sigma, el, means, 1-r2, mask.copy(), sample, d, point,
        prices is None, expenditure is None, nl.n_outer, nl.n_gn, nl.converged,
        tau=tau, setau=None if fs is None else fs.setau,
        np_prob=0 if fs is None else fs.np_prob,
        probits=[] if fs is None else list(fs.results), selection_layout=layout,
        selection_covariate_names=selected_covariates if censor else [],
        selection_index=None if fs is None else fs.du.copy(), evaluation_design=point_design,
        start_method=start, algorithm=algorithm, stop_rule=stop_rule, nrtol=nl.nrtol,
        initialization=initialization)
    domain_limited = False
    if censor:
        inner = _inner(nl.theta, d, spec)
        domain_limited = bool(np.min(inner.b) < 1e-6 or np.min(abs(inner.denominator)) < 1e-6)
    if censor and domain_limited and (not np.isfinite(nl.nrtol) or nl.nrtol > 1e-5):
        res.converged = False
        res.notes.append("The fit is near the effective-expenditure/denominator boundary with a large scaled gradient; a domain-limited tiny step does not establish convergence.")
    if not res.converged:
        res.notes.append("The nonlinear SUR optimizer did not converge.")
    if np.any(el.we < 0):
        res.notes.append("Some predicted shares at the evaluation point are negative.")
    if np.any(fitted < 0):
        res.notes.append(f"Diagnostic: {int(np.sum(fitted < 0))} fitted share values are negative; values were not clipped.")
    if censor:
        res.notes.extend([
            "Latent adding-up is relaxed: all alpha coefficients are free; fitted shares are not renormalized.",
            "Compensated elasticities are a Slutsky transformation; utility-consistent Hicksian demand is not guaranteed.",
            "Analytical structural and elasticity SEs condition on fitted Probits and evaluation means; use full two-step bootstrap for generated-regressor inference.",
        ])
    if reps:
        from .bootstrap import bootstrap_estimator
        kwargs = dict(shares=shares, prices=prices, lnprices=lnprices,
                      expenditure=expenditure, lnexpenditure=lnexpenditure,
                      demographics=demo_names, censor=censor, method=method, start=start,
                      algorithm=algorithm, stop_rule=stop_rule, vce_sigma=vce_sigma,
                      tol=tol, max_outer=max_outer, max_iter=max_iter, chunk=chunk,
                      nrtol_stop=nrtol_stop, inner_nrtol_early=inner_nrtol_early,
                      sigma_tol=sigma_tol, share_sum_tol=share_sum_tol,
                      blas_threads=blas_threads,
                      initial=nl.theta if bootstrap_start == "warm" else None,
                      sigma_initial=nl.sigma if bootstrap_start == "warm" else None)
        if censor:
            kwargs.update(selection_prices=selected_prices, selection_expenditure=selection_expenditure,
                          selection_covariates=selected_covariates)
        res.boot = bootstrap_estimator(
            data=sample.reset_index(drop=True), estimator="translog", kwargs=kwargs,
            reps=int(reps), seed=seed, n_jobs=int(n_jobs), blas_threads=blas_threads,
            verbose=verbose, mp_context=mp_context, rep_timeout=rep_timeout)
        res.V_analytic = res.V.copy()
        res.V = res.boot.V.copy()
        if censor:
            res.notes.append("Reported V and SEs use full two-step bootstrap; V_analytic retains conditional inference.")
    return res
