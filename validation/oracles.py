"""Independent equations and finite differences, not archived Stata output.

Tolerances are deliberately declared here before the experiment. A failed
comparison is retained; the runner does not loosen tolerances or reroll seeds.
"""
from dataclasses import replace
import numpy as np
from scipy.stats import norm, f as f_distribution
from scipy.optimize import minimize
from pyquaidsce import (DemandData, fitted_shares, jacobian_full, delta_matrix,
                       full_vector, unpack, elasticities, fitted_share_derivatives)
from pyquaidsce.elasticities import Means
from pyquaidsce._threads import blas_thread_limit

TOLERANCES = {
    "algebra_atol": 1e-10, "algebra_rtol": 1e-10,
    "finite_difference_scaled": 2e-6,
    "equivalent_fit_atol": 1e-6, "equivalent_fit_rtol": 1e-5,
    "stationarity": 1e-9,
    "covariance_atol": 1e-7, "covariance_rtol": 1e-5,
    "probit_coefficient_atol": 2e-5, "probit_coefficient_rtol": 2e-5,
}


class Checks:
    def __init__(self):
        self.rows = []

    def test(self, name, condition, observed=None, expected=None):
        self.rows.append(dict(name=name, passed=bool(condition), observed=observed, expected=expected))

    def close(self, name, actual, expected, atol=1e-10, rtol=1e-10):
        a, b = np.asarray(actual), np.asarray(expected)
        shape_ok = a.shape == b.shape
        error = float(np.max(np.abs(a - b))) if shape_ok and a.size else (0.0 if shape_ok else None)
        self.test(name, shape_ok and np.allclose(a, b, atol=atol, rtol=rtol, equal_nan=True),
                  dict(max_abs_error=error, actual_shape=list(a.shape), expected_shape=list(b.shape)),
                  dict(atol=atol, rtol=rtol))

    def scaled(self, name, actual, expected, tolerance=2e-6):
        a, b = np.asarray(actual), np.asarray(expected)
        err = float(np.max(np.abs(a - b) / np.maximum(1.0, np.abs(b))))
        self.test(name, np.isfinite(err) and err <= tolerance, err, tolerance)


def build_data(frame, kw, res):
    prices = kw.get("prices") or kw.get("lnprices")
    demos = kw.get("demographics") or []
    z_names = kw.get("selection_covariates")
    if z_names is None:
        z_names = demos
    need = list(dict.fromkeys(list(kw["shares"]) + list(prices) + list(demos)
                             + [kw.get("expenditure") or kw.get("lnexpenditure")]
                             + list(kw.get("ivexp") or []) + list(z_names)
                             + [x for x in (kw.get("control_function"), kw.get("selection_control_function")) if x]))
    mask = np.isfinite(frame[need].to_numpy(float)).all(axis=1)
    df = frame.iloc[np.flatnonzero(mask)]
    lp = df[prices].to_numpy(float)
    if kw.get("prices") is not None:
        lp = np.log(lp)
    exp_name = kw.get("expenditure") or kw.get("lnexpenditure")
    lm = df[exp_name].to_numpy(float)
    if kw.get("expenditure") is not None:
        lm = np.log(lm)
    z = df[demos].to_numpy(float) if demos else np.zeros((len(df), 0))
    cf = np.zeros(len(df))
    if res.reduced_form is not None:
        cf = res.reduced_form.residuals
    elif kw.get("control_function"):
        cf = df[kw["control_function"]].to_numpy(float)
    n = len(prices)
    cdf, pdf, index = np.ones((len(df), n)), np.zeros((len(df), n)), np.ones((len(df), n))
    design = np.zeros((len(df), 0))
    if res.spec.censor:
        sp = kw.get("selection_prices")
        if sp is None:
            sp = prices
        parts = [lp[:, prices.index(p)] for p in sp]
        include_exp = kw.get("selection_expenditure", True)
        if include_exp is None:
            include_exp = kw.get("expenditure") is not None
        if include_exp:
            parts.append(lm)
        parts.extend(df[x].to_numpy(float) for x in z_names)
        if res.reduced_form is not None:
            parts.append(res.reduced_form.residuals)
        elif kw.get("selection_control_function"):
            parts.append(df[kw["selection_control_function"]].to_numpy(float))
        design = np.column_stack(parts) if parts else np.zeros((len(df), 0))
        b = res.tau.reshape(n, -1)
        index = design @ b[:, :-1].T + b[:, -1]
        cdf, pdf = norm.cdf(index), norm.pdf(index)
    return DemandData(lp, lm, df[kw["shares"]].to_numpy(float), z, cdf, pdf,
                      float(kw["anot"]), cf), design, index, mask


def predict(c, lp, lm, z, a0, quadratic):
    """Direct published Ray/QUAIDS equations, separate from model._inner."""
    lp, lm, z = np.asarray(lp), np.asarray(lm), np.asarray(z)
    out = np.empty((len(lp), len(c.alpha)))
    # Deliberately use a row-wise formulation rather than model's vectorization.
    for t in range(len(lp)):
        ray = 1.0 + float(z[t] @ c.rho)
        A = a0 + float(c.alpha @ lp[t]) + 0.5 * float(lp[t] @ c.gamma @ lp[t])
        D = lm[t] - A - np.log(ray)
        B = c.beta + z[t] @ c.eta
        out[t] = c.alpha + c.gamma @ lp[t] + B * D
        if quadratic:
            q = np.exp(-float(B @ lp[t]))
            out[t] += c.lam * q * D**2
    return out


def independent_fitted(theta, d, spec):
    c = unpack(theta, spec)
    latent = predict(c, d.lnp, d.lnexp, d.demo, d.a0, spec.quadratic)
    if spec.censor:
        return (latent + d.control_function[:, None] * c.cfcoef) * d.cdf + c.delta * d.pdf
    return latent[:, :spec.neqn - 1]


def numerical_jacobian(theta, d, spec):
    J = np.empty((d.nobs, spec.n_eq_estimated, len(theta)))
    for k in range(len(theta)):
        h = 1e-6 * max(1.0, abs(float(theta[k])))
        plus, minus = theta.copy(), theta.copy()
        plus[k] += h
        minus[k] -= h
        J[:, :, k] = (independent_fitted(plus, d, spec) - independent_fitted(minus, d, spec)) / (2 * h)
    return J


def check_fit(check, frame, kw, res, expected_nonconvergence=False):
    d, design, index, mask = build_data(frame, kw, res)
    n, m, k = res.spec.neqn, res.spec.n_eq_estimated, res.spec.n_free
    check.test("convergence_contract", res.converged != expected_nonconvergence,
               res.converged, not expected_nonconvergence)
    if expected_nonconvergence:
        check.test("nonconvergence_note", any("did not satisfy" in x for x in res.notes))
    check.test("complete_sample_size", res.nobs == int(mask.sum()), res.nobs, int(mask.sum()))
    check.test("parameter_dimensions", len(res.theta) == k and len(res.b) == len(res.names)
               and res.V.shape == (len(res.b), len(res.b)))
    check.test("unique_parameter_names", len(set(res.names)) == len(res.names))
    for label, array in (("b", res.b), ("V", res.V), ("theta", res.theta),
                         ("sigma", res.sigma), ("V_est", res.V_est)):
        check.test(f"finite_{label}", np.isfinite(array).all())
    check.close("V_symmetric", res.V, res.V.T)
    vscale = max(1.0, float(np.max(np.abs(res.V))))
    eigenvalues = np.linalg.eigvalsh(res.V)
    check.test("V_positive_semidefinite", eigenvalues.min() >= -1e-9 * vscale,
               float(eigenvalues.min()), -1e-9 * vscale)
    check.test("sigma_dimensions", res.sigma.shape == (m, m))
    check.test("sigma_positive_definite", np.linalg.eigvalsh(res.sigma).min() > 0)
    check.test("sigma_numerical_rank", np.linalg.matrix_rank(res.sigma) == m,
               int(np.linalg.matrix_rank(res.sigma)), m)
    c = res.coefs
    check.close("gamma_symmetry", c.gamma, c.gamma.T)
    check.close("gamma_homogeneity", c.gamma.sum(axis=1), np.zeros(n))
    if res.spec.ndemo:
        check.close("eta_adding_up", c.eta.sum(axis=1), np.zeros(res.spec.ndemo))
        check.test("ray_positive_domain", np.all(1.0 + d.demo @ c.rho > 0))
    if not res.spec.quadratic:
        check.close("linear_lambda_zero", c.lam, np.zeros(n))
    if not res.spec.censor:
        check.close("uncensored_alpha_adding_up", c.alpha.sum(), 1.0)
        check.close("uncensored_beta_adding_up", c.beta.sum(), 0.0)
        check.close("uncensored_lambda_adding_up", c.lam.sum(), 0.0)
    # Independent published formula for all rows; no clipping is allowed.
    mu = independent_fitted(res.theta, d, res.spec)
    actual = fitted_shares(res.theta, d, res.spec)
    check.close("independent_fitted_shares", actual, mu, atol=2e-10, rtol=2e-10)
    u = d.shares[:, :m] - mu
    sigma = u.T @ u / len(u)
    check.close("residual_covariance", res.sigma, sigma, atol=2e-11, rtol=2e-10)
    _, logdet = np.linalg.slogdet(sigma)
    ll = -len(u) * (m * (1 + np.log(2 * np.pi)) + logdet) / 2
    check.close("gaussian_log_likelihood", res.llf, ll, atol=1e-7, rtol=1e-10)
    check.close("free_to_full_coefficients", res.b[:res.spec.n_full], full_vector(res.theta, res.spec))
    if res.spec.censor:
        check.test("probit_count", len(res.probits) == n)
        check.test("probit_convergence", all(x.converged for x in res.probits))
        check.test("selection_width", res.np_prob == design.shape[1] + 1)
        for i, prob in enumerate(res.probits):
            check.close(f"probit_tau_{i}", res.tau[i * res.np_prob:(i + 1) * res.np_prob], prob.b)
            sl = slice(i * res.np_prob, (i + 1) * res.np_prob)
            check.close(f"probit_covariance_{i}", res.setau[sl, sl], prob.V)
            check.close(f"probit_linear_index_{i}", prob.xb(design), index[:, i])
        check.close("mean_cdf_from_xb", res.means.cdf, d.cdf.mean(axis=0))
        check.close("mean_pdf_from_xb", res.means.pdf, d.pdf.mean(axis=0))
        check.close("mean_du_from_xb", res.means.du, index.mean(axis=0))
    E = res.elas
    check.test("elasticity_shapes", E.income.shape == (n,) and E.uncompensated.shape == (n, n)
               and E.compensated.shape == (n, n))
    check.test("finite_elasticities", all(np.isfinite(x).all() for x in (E.income, E.uncompensated, E.compensated)))
    check.close("Slutsky_identity", E.compensated, E.uncompensated + np.outer(E.income, res.means.w))
    if res.spec.censor:
        count = n + 2 * n * n
        check.close("reported_elasticity_block", res.b[-count:], E.as_stata_vector())
        check.test("analytic_elasticity_se_undefined", np.isnan(res.analytic_se[-count:]).all())
        check.test("analytic_nonelasticity_se_finite", np.isfinite(res.analytic_se[:-count]).all())
    else:
        check.test("uncensored_analytic_se_finite", np.isfinite(res.analytic_se).all())
    D = delta_matrix(res.spec)
    if res.boot is None:
        check.close("delta_covariance_mapping", res.V[:res.spec.n_full, :res.spec.n_full],
                    D @ res.V_est @ D.T, atol=1e-9, rtol=1e-9)
    # Finite differences at deterministic small4 rows, including all free parameters.
    pick = np.linspace(0, d.nobs - 1, 12, dtype=int)
    small = d.subset(pick)
    num = numerical_jacobian(res.theta, small, res.spec)
    analytic = np.einsum("tmf,fk->tmk", jacobian_full(res.theta, small, res.spec), D)
    check.scaled("independent_numerical_jacobian", analytic, num)
    return d, design, index, mask, dict(negative_fitted_values=int(np.sum(mu < 0)),
                                      residual_sum_squares=float(np.sum(u**2)),
                                      llf_independent=float(ll),
                                      sigma_condition_number=float(np.linalg.cond(res.sigma)))


def reduced_form_checks(check, d, frame, kw, res, mask):
    rf = res.reduced_form
    iv = frame.iloc[np.flatnonzero(mask)][kw["ivexp"]].to_numpy(float)
    X = np.column_stack([d.lnp, d.demo, iv, np.ones(d.nobs)])
    Q, R = np.linalg.qr(X, mode="reduced")
    beta = np.linalg.solve(R, Q.T @ d.lnexp)
    u = d.lnexp - X @ beta
    sigma2 = float(u @ u) / (d.nobs - X.shape[1])
    invR = np.linalg.inv(R)
    V = sigma2 * invR @ invR.T
    X0 = np.column_stack([d.lnp, d.demo, np.ones(d.nobs)])
    Q0, R0 = np.linalg.qr(X0, mode="reduced")
    u0 = d.lnexp - X0 @ np.linalg.solve(R0, Q0.T @ d.lnexp)
    q = iv.shape[1]
    F = max(0.0, float(u0 @ u0 - u @ u)) / q / sigma2
    check.close("OLS_QR_coefficients", rf.b, beta, atol=2e-9, rtol=2e-9)
    check.close("OLS_QR_residuals", rf.residuals, u, atol=2e-10, rtol=2e-10)
    check.close("OLS_QR_covariance", rf.V, V, atol=2e-9, rtol=2e-8)
    check.close("OLS_residual_orthogonality", X.T @ rf.residuals, np.zeros(X.shape[1]), atol=2e-7, rtol=0)
    check.close("excluded_instrument_F", rf.excluded_f, F, atol=2e-8, rtol=2e-8)
    check.close("excluded_instrument_p", rf.excluded_pvalue, f_distribution.sf(F, q, rf.df_resid))
    check.test("OLS_degrees_of_freedom", rf.df_resid == d.nobs - X.shape[1]
               and rf.excluded_df_num == q)


def observation_derivative_checks(check, res, d, index):
    if not res.spec.censor:
        return
    pick = np.linspace(0, d.nobs - 1, 8, dtype=int)
    dd = d.subset(pick)
    kk = index[pick].copy()
    exp, price = fitted_share_derivatives(res.theta, dd, res.spec, tau=res.tau,
                                         layout=res.selection_layout, selection_index=kk)
    b = res.tau.reshape(res.spec.neqn, -1)
    h = 1e-5
    def shifted(j, direction):
        lp, lm, k = dd.lnp.copy(), dd.lnexp.copy(), kk.copy()
        pos = res.selection_layout.expenditure_position if j is None else res.selection_layout.price_position(j)
        if j is None:
            lm += direction * h
        else:
            lp[:, j] += direction * h
        if pos is not None:
            k += direction * h * b[:, pos]
        nd = DemandData(lp, lm, dd.shares, dd.demo, norm.cdf(k), norm.pdf(k), dd.a0, dd.control_function)
        return independent_fitted(res.theta, nd, res.spec)
    check.scaled("fitted_expenditure_derivative", exp, (shifted(None, 1) - shifted(None, -1)) / (2 * h))
    for j in range(res.spec.neqn):
        check.scaled(f"fitted_price_derivative_{j}", price[:, :, j], (shifted(j, 1) - shifted(j, -1)) / (2 * h))


def model_point_elasticity_checks(check, res):
    """Finite-difference quantities at one model-consistent point.

    Reported elasticities use observed mean shares and separate mean Probit
    moments. These are not silently equated with mean fitted elasticities.
    This oracle instead constructs a consistent point for the public formulas.
    """
    c, spec, mean = res.coefs, res.spec, res.means
    lp = mean.lnp[None, :].copy()
    lm = np.array([mean.lnexp])
    z = mean.demo[None, :]
    w = predict(c, lp, lm, z, res.anot, spec.quadratic)[0]
    if np.any(np.abs(w) < 1e-8):
        check.test("model_point_nondegenerate", False, w.tolist(), "nonzero latent shares")
        return
    k = mean.du.copy() if spec.censor else np.zeros(spec.neqn)
    cdf, pdf = (norm.cdf(k), norm.pdf(k)) if spec.censor else (np.ones(spec.neqn), np.zeros(spec.neqn))
    mm = Means(w, mean.lnp, mean.lnexp, mean.demo, cdf, pdf, k, mean.control_function)
    e = elasticities(c, spec, mm, res.anot, tau=res.tau, np_prob=res.np_prob,
                     layout=res.selection_layout)
    t = res.tau.reshape(spec.neqn, -1) if spec.censor else None
    def shares(delta_exp=0.0, j=None, delta_price=0.0):
        ll = lp.copy()
        if j is not None:
            ll[0, j] += delta_price
        raw = predict(c, ll, lm + delta_exp, z, res.anot, spec.quadratic)[0]
        if not spec.censor:
            return raw
        ki = k.copy()
        xp = spec.censor and res.selection_layout.expenditure_position
        if xp is not None:
            ki += delta_exp * t[:, xp]
        if j is not None:
            pp = res.selection_layout.price_position(j)
            if pp is not None:
                ki += delta_price * t[:, pp]
        return (raw + c.cfcoef * mean.control_function) * norm.cdf(ki) + c.delta * norm.pdf(ki)
    s0 = shares()
    if np.any(np.abs(s0) < 1e-8):
        check.test("model_point_nondegenerate", False, s0.tolist(), "nonzero fitted shares")
        return
    h = 1e-5
    de = (shares(delta_exp=h) - shares(delta_exp=-h)) / (2 * h)
    eu = np.empty((spec.neqn, spec.neqn))
    for j in range(spec.neqn):
        ds = (shares(j=j, delta_price=h) - shares(j=j, delta_price=-h)) / (2 * h)
        eu[:, j] = ds / s0 - (np.arange(spec.neqn) == j)
    check.scaled("model_point_income_elasticity", e.income, 1.0 + de / s0)
    check.scaled("model_point_price_elasticity", e.uncompensated, eu)
    if not spec.censor:
        check.close("model_point_Engel", w @ e.income, 1.0, atol=1e-8, rtol=1e-8)
        check.close("model_point_Cournot", e.uncompensated.sum(axis=1), -e.income, atol=1e-8, rtol=1e-8)
        check.close("model_point_Hicksian", e.compensated.sum(axis=1), np.zeros(spec.neqn), atol=1e-8, rtol=1e-8)


def stationary_measure(res, d, sigma=None):
    """Scale-free gradient computed from numerical, independent predictions."""
    J = numerical_jacobian(res.theta, d, res.spec)
    u = d.shares[:, :res.spec.n_eq_estimated] - independent_fitted(res.theta, d, res.spec)
    precision = np.linalg.inv(res.sigma if sigma is None else sigma)
    G = np.einsum("tik,ij,tjl->kl", J, precision, J)
    g = np.einsum("tik,ij,tj->k", J, precision, u)
    obj = float(np.einsum("ti,ij,tj->", u, precision, u))
    return abs(float(g @ np.linalg.lstsq(G, g, rcond=None)[0])) / max(obj, 1e-300)


def independent_probit(check, res, design, shares):
    scale = np.std(design, axis=0)
    scale[scale == 0] = 1
    Z = np.column_stack([design / scale, np.ones(len(design))])
    for i, prob in enumerate(res.probits):
        q = 2 * (shares[:, i] > 0).astype(float) - 1
        def loss(b):
            return -float(norm.logcdf(q * (Z @ b)).sum())
        opt = minimize(loss, np.zeros(Z.shape[1]), jac="3-point", method="BFGS",
                       options=dict(gtol=1e-7, maxiter=1500))
        theta = opt.x.copy()
        theta[:-1] /= scale
        check.close(f"independent_probit_ll_{i}", prob.llf, -opt.fun, atol=1e-7, rtol=1e-9)
        check.close(f"independent_probit_b_{i}", prob.b, theta, atol=2e-5, rtol=2e-5)


def bootstrap_checks(check, res):
    boot = res.boot
    check.test("bootstrap_present", boot is not None)
    if boot is None:
        return
    B = boot.b_star
    check.test("bootstrap_dimensions", B.shape == (boot.reps_ok, len(res.b)))
    check.test("bootstrap_accounting", boot.reps_ok + len(boot.failures) == boot.reps_requested)
    check.test("bootstrap_at_least_two", boot.reps_ok >= 2)
    check.test("bootstrap_finite_draws", np.isfinite(B).all())
    check.close("bootstrap_se_ddof1", res.se, B.std(axis=0, ddof=1))
    centered = B - B.mean(axis=0)
    V = centered.T @ centered / (len(B) - 1)
    check.close("bootstrap_covariance_ddof1", res.V, V)
    check.close("bootstrap_covariance_synchronization", res.V, boot.V)
    check.test("analytic_reference_retained", res.V_analytic is not None)
    if res.spec.censor:
        count = res.spec.neqn + 2 * res.spec.neqn ** 2
        check.test("bootstrap_elasticity_se_finite", np.isfinite(res.se[-count:]).all())
    for level in (90, 95, 99):
        lo, hi = boot.ci(res.b, level)
        z = norm.ppf(0.5 + level / 200)
        check.close(f"normal_ci_lower_{level}", lo, res.b - z * res.se)
        check.close(f"normal_ci_upper_{level}", hi, res.b + z * res.se)
        lo, hi = boot.percentile_ci(level)
        check.close(f"percentile_ci_lower_{level}", lo, np.quantile(B, (100-level)/200, axis=0))
        check.close(f"percentile_ci_upper_{level}", hi, np.quantile(B, (100+level)/200, axis=0))
