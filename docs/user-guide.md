# User Guide & Complete API Reference

This document provides a comprehensive reference for all input parameters of `quaidsce(...)` and the structure of the returned `QuaidsceResults` object.

---

## 1. Input Arguments Reference: `quaidsce(...)`

### Data and Variable Specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | `pandas.DataFrame` | *Required* | The dataset containing household expenditure shares, prices, total expenditure, and demographic variables. |
| `shares` | `Sequence[str]` | *Required* | List of column names for budget shares ($w_1, \dots, w_n$). All $n$ goods must be listed. For censored estimation (`censor=True`), zeros must be present. |
| `prices` | `Sequence[str]` | `None` | List of column names for prices ($p_1, \dots, p_n$), in the exact same order as `shares`. Exactly one of `prices` or `lnprices` must be specified. |
| `lnprices` | `Sequence[str]` | `None` | List of column names for log prices ($\ln p_1, \dots, \ln p_n$), in the exact same order as `shares`. |
| `expenditure` | `str` | `None` | Column name for total system expenditure ($m$). Exactly one of `expenditure` or `lnexpenditure` must be specified. |
| `lnexpenditure` | `str` | `None` | Column name for log total expenditure ($\ln m$). |
| `demographics` | `Sequence[str]` | `None` | List of column names for demographic variables ($z_1, \dots, z_R$) used in Ray (1983) demographic scaling. Required when `censor=True`. |
| `anot` | `float` | *Required* | Constant $\alpha_0$ in the translog price index $\ln a(p)$ (equivalent to Stata's `anot(#)`). |

---

### Model Specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `quadratic` | `bool` | `True` | If `True`, estimates the Quadratic AIDS (QUAIDS) model. If `False`, estimates the linear AIDS model (equivalent to Stata's `noquadratic`). |
| `censor` | `bool` | `True` | If `True`, applies the Shonkwiler & Yen (1999) two-step censoring correction. If `False`, estimates uncensored QUAIDS (equivalent to Stata's `nocensor` / Poi 2012). |

---

### Control Function and Selection Design

| Parameter | Type | Default | Description |
|---|---|---|---|
| `control_function` | `str` | `None` | Column containing an externally generated residual. It enters latent demand share `i` as `cfcoef_i * residual` inside the Shonkwiler–Yen CDF multiplier. |
| `selection_prices` | `Sequence[str]` | `None` | Ordered subset of demand-price columns used in each first-stage Probit. `None` keeps all prices; `[]` removes them. |
| `selection_expenditure` | `bool` or `None` | `True` | Include log expenditure in each Probit. `None` reproduces the legacy branch-dependent behavior. |
| `selection_covariates` | `Sequence[str]` | `None` | Independent Probit covariates. `None` reuses Ray demographics; `[]` uses none. |
| `selection_control_function` | `str` | `None` | Residual column used in each Probit with an equation-specific coefficient independent of `cfcoef`. |

These extensions require `censor=True` and `first_stage_predict="xb"`.
Internal `reps>0` is disabled when an external residual is active because a
valid generated-regressor bootstrap must re-estimate the reduced form inside
every replication.

---

### Estimation Methods and Optimizer Settings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `method` | `str` | `"fgnls"` | Second-stage estimation method: `"nls"` (Nonlinear Least Squares), `"fgnls"` (Feasible Generalized NLS), or `"ifgnls"` (Iterated FGNLS). |
| `algorithm` | `str` | `"gn"` | Numerical optimization algorithm: `"gn"` (Gauss-Newton with Hartley step halving, matching Stata) or `"lm"` (Levenberg-Marquardt trust-region damping). |
| `start` | `str` | `"zero"` | Starting value scheme: `"zero"` (starts all parameters at 0, matching Stata default) or `"linear"` (fits a fast linearized AIDS starting point). |
| `initial` | `ArrayLike` | `None` | Optional 1D vector of initial free parameters of length `res.spec.n_free` (useful for restarting or warm-starting). |
| `sigma_initial` | `np.ndarray` | `None` | Optional initial $(m \times m)$ residual covariance matrix $\Sigma$ (used in combination with `initial` for warm-starting IFGNLS). |

---

### Stata Compatibility Switches

| Parameter | Type | Default | Description |
|---|---|---|---|
| `first_stage_predict` | `str` | `"pr"` | First-stage Probit prediction: `"pr"` reproduces Stata's default $\Phi(\Phi(X'\tau))$ and $\phi(\Phi(X'\tau))$; `"xb"` uses the linear index $X'\tau$ (textbook Shonkwiler–Yen). |
| `strict_stata` | `bool` | `False` | If `True`, reproduces Stata's exact documented elasticity calculations (including its quirks). If `False` *(default)*, applies published corrected formulas. |
| `vce_sigma` | `str` | `"objective"` | Residual covariance used in the second-stage standard error formula: `"objective"` (used in the final minimization, matching Stata) or `"final"` (recomputed from final residuals). |

---

### Convergence and Numerical Tolerance

| Parameter | Type | Default | Description |
|---|---|---|---|
| `tol` | `float` | `1e-13` | Objective relative change tolerance. |
| `nrtol_stop` | `float` | `1e-12` | Scaled relative gradient stopping tolerance (Gauss-Newton stationarity). |
| `sigma_tol` | `float` | `1e-11` | Outer fixed-point parameter tolerance for IFGNLS iterations. |
| `inner_nrtol_early` | `float` | `1e-8` | Early-stage inner Gauss-Newton tolerance used during inexact-outer IFGNLS. |
| `stop_rule` | `str` | `"standard"` | Stopping criterion: `"standard"` (disjunctive rule matching Stata's `tolerance` / `ltolerance` / `nrtolerance`) or `"tight"` (stops only on the scaled gradient). |
| `max_iter` | `int` | `300` | Maximum number of inner Gauss-Newton iterations per stage. |
| `max_outer` | `int` | `200` | Maximum number of outer covariance updates for IFGNLS. |
| `chunk` | `int` | `2000` | Observation block size for accumulating normal equations ($J'\Sigma^{-1}J$) without materializing the full Jacobian in RAM. |

---

### Bootstrap and Parallelism

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reps` | `int` | `0` | Number of bootstrap replications. Set `reps=0` to disable bootstrap (point estimates only). |
| `seed` | `int` | `None` | Random seed for reproducible bootstrap resamples. |
| `n_jobs` | `int` | `1` | Number of parallel worker processes for bootstrap replications. |
| `bootstrap_start` | `str` | `"zero"` | Starting values for each bootstrap draw: `"zero"` (starts each draw from zero) or `"warm"` (starts each draw from full-sample estimates). |
| `boot_sigma_tol` | `float` | `1e-7` | Outer covariance tolerance used inside bootstrap replications. |
| `mp_context` | `str` or `None` | `None` | Multiprocessing start method. `None` selects the BLAS-safe cross-platform `"spawn"` default; an available method such as `"forkserver"` can be requested explicitly. |
| `rep_timeout` | `float` or `None` | `None` | Per-replication wall-clock limit. Cooperative Probit/optimizer checks are backed by a parent watchdog that terminates the disposable child process if necessary. Adds process-start overhead when enabled. |

---

### Logging

| Parameter | Type | Default | Description |
|---|---|---|---|
| `verbose` | `bool` | `True` | If `True`, prints high-level estimation and bootstrap progress. |
| `gn_verbose` | `bool` | `False` | If `True`, prints detailed step-by-step Gauss-Newton optimization logs. |

---

## 2. Result Object Reference: `QuaidsceResults`

The object returned by `quaidsce(...)` contains all estimated parameters, standard errors, covariance matrices, elasticities, and diagnostic metadata.

### Overview Table: Python vs. Stata Equivalents

| Attribute / Property | Python Type | Stata Equivalent | Description |
|---|---|---|---|
| `res.b` | `np.ndarray` | `e(b)` | Full parameter vector: structural parameters, Probit parameters ($\tau$), and elasticity estimates. |
| `res.V` | `np.ndarray` | `e(V)` | Active covariance matrix of `res.b`: bootstrap covariance when `reps>0`, otherwise the finite conditional analytical matrix. |
| `res.se` | `np.ndarray` | standard errors | Bootstrap S.E.s when available; otherwise conditional analytical S.E.s with unsupported elasticity entries explicitly set to `NaN`. |
| `res.V_analytic` | `np.ndarray` or `None` | — | Conditional analytical reference retained after bootstrap. |
| `res.analytic_se` | `np.ndarray` | — | Conditional analytical standard errors. Elasticity entries are `NaN` because their analytical delta-method covariance is not implemented. |
| `res.names` | `List[str]` | `colnames e(b)` | Parameter labels in `"equation:name"` format. |
| `res.theta` (or `res.b_est`) | `np.ndarray` | `e(best)` | The vector of estimated free structural parameters. |
| `res.V_est` | `np.ndarray` | `e(Vest)` | Covariance matrix of the free parameters. |
| `res.coefs` | `Coefs` dataclass | `e(alpha)` … `e(rho)` | Unpacked structural parameters ($\alpha, \beta, \gamma, \lambda, \delta, \eta, \rho$). |
| `res.tau` | `np.ndarray` | `e(tau)` | Stacked first-stage Probit coefficients. |
| `res.setau` | `np.ndarray` | `e(V_tau)` | Block-diagonal covariance matrix of first-stage Probits. |
| `res.probits` | `List[ProbitResult]` | `e()` from each probit | Individual Probit estimation result objects. |
| `res.elas` | `Elasticities` | `e(elas_*)` | Demand elasticities evaluated at sample means. |
| `res.llf` | `float` | `e(ll)` | Gaussian log-likelihood value. |
| `res.sigma` | `np.ndarray` | `e(Sigma)` | Residual covariance matrix ($(n \times n)$ or $(n-1 \times n-1)$). |
| `res.nobs` | `int` | `e(N)` | Number of observations in the estimation sample. |
| `res.converged` | `bool` | `e(converged)` | `True` if optimization converged successfully. |
| `res.boot` | `BootResult` or `None` | `r(table)` after `bs` | Bootstrap replicates, standard errors, and percentile confidence intervals (when `reps > 0`). |
| `res.notes` | `List[str]` | — | Informational notes regarding Stata compatibility switches. |

---

### Detailed Access to Structural Coefficients (`res.coefs`)

The structural coefficients are accessible directly as NumPy arrays:

```python
# Individual structural parameter blocks
alpha = res.coefs.alpha   # (n,) array: intercepts
beta  = res.coefs.beta    # (n,) array: price index slopes
gamma = res.coefs.gamma   # (n, n) symmetric array: price interaction terms
lam   = res.coefs.lam     # (n,) array: quadratic expenditure terms
delta = res.coefs.delta   # (n,) array: Shonkwiler-Yen censoring terms
eta   = res.coefs.eta     # (R, n) array: demographic interaction terms
rho   = res.coefs.rho     # (R,) array: Ray demographic scaling factors
```

---

### Detailed Access to Elasticities (`res.elas`)

Elasticities are computed at sample means and stored as natural NumPy arrays:

```python
# 1. Expenditure (Income) Elasticity for each good: shape (n,)
income_elas = res.elas.income

# 2. Uncompensated (Marshallian) Price Elasticities: shape (n, n)
# Row i = good i, Column j = price of good j
uncomp_elas = res.elas.uncompensated

# 3. Compensated (Hicksian) Price Elasticities: shape (n, n)
# Row i = good i, Column j = price of good j
comp_elas = res.elas.compensated
```

---

### Accessing Bootstrap Results (`res.boot`)

When `reps > 0`, `res.boot` contains the distribution of bootstrap replicates:

```python
if res.boot is not None:
    # Standard errors across all parameters in res.names
    se = res.boot.se

    # 95% Normal-based confidence intervals: (lower_bounds, upper_bounds)
    ci_low, ci_high = res.boot.ci(res.b, level=95.0)

    # 95% Empirical Percentile confidence intervals: (lower_bounds, upper_bounds)
    pct_low, pct_high = res.boot.percentile_ci(level=95.0)

    # Full matrix of bootstrap replicates: shape (reps_ok, n_parameters)
    replicate_draws = res.boot.b_star

    # These are synchronized after a successful bootstrap:
    assert np.allclose(res.V, res.boot.V)
    assert np.allclose(res.se, res.boot.se)
```

---

### Helper Methods on `QuaidsceResults`

- **`res.summary(level=95.0)`**: Returns a formatted string with model header statistics and a 78-column coefficient table matching Stata's `_coef_table` format.
- **`res.elasticity_tables()`**: Returns formatted summary tables for expenditure, Marshallian, and Hicksian elasticities.
- **`res.named()`**: Returns a Python dictionary mapping parameter names (`"equation:name"`) to their estimated values.
- **`res.get("beta_1")`**: Look up a specific coefficient by its bare name or equation name.

---

## 3. Practical Code Examples

### Converting Elasticity Matrices to pandas DataFrames

```python
import pandas as pd

# Convert Marshallian elasticities to a labeled DataFrame
good_names = ["Meat", "Dairy", "Cereals", "Other"]

df_uncomp = pd.DataFrame(
    res.elas.uncompensated,
    index=good_names,      # Rows: quantity of good i
    columns=good_names     # Columns: price of good j
)

print("Marshallian Price Elasticities:")
print(df_uncomp)
```

### Saving Estimation Results to CSV

```python
# Build a summary DataFrame
df_results = pd.DataFrame({
    "parameter": res.names,
    "coefficient": res.b,
    "std_error": res.boot.se if res.boot else res.se,
})

# Save to CSV
df_results.to_csv("quaidsce_estimates.csv", index=False)
```
