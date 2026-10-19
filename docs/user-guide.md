# User Guide & Complete API Reference

This guide explains the arguments accepted by `quaidsce(...)` and how to read its results. The defaults below are for Python. For a first run, see [Getting Started](getting-started.md). For complete scripts in Python, Stata, and R, see [Examples](../examples/README.md).

---

## 1. Input Arguments Reference: `quaidsce(...)`

### Data and Variable Specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `data` | `pandas.DataFrame` | *Required* | Your dataset. Use one row per household and numeric columns for the variables used in the model. |
| `shares` | List of column names | *Required* | Spending shares, such as `["w1", "w2", "w3", "w4"]`. A share is spending on a good divided by total spending on the goods in your system. List at least three goods and include them all. With `censor=True`, each share column must have both zero and positive values in the estimation sample. |
| `prices` | List of column names | `None` | Prices before taking logs, such as `["p1", "p2", "p3", "p4"]`. Give one price per share, in the same order. Prices must be positive. Use either `prices` or `lnprices`, not both. |
| `lnprices` | List of column names | `None` | Prices that you have already converted to natural logs. The package uses these values directly. Keep the same good order as `shares`. |
| `expenditure` | Column name | `None` | Total spending on the goods in your system, before taking logs. Values must be positive. Use either `expenditure` or `lnexpenditure`, not both. |
| `lnexpenditure` | Column name | `None` | Total system expenditure that you have already converted to its natural log. The package uses these values directly. |
| `demographics` | List of column names | `None` | Household characteristics, such as `["hh_size", "urban"]`, that enter the demand equations through Ray's demographic scaling. At least one is required when `censor=True`; they are optional when `censor=False`. Code categories as numeric variables, for example `urban=0` or `1`. |
| `ivexp` | List of column names | `None` | Instruments used to account for endogenous expenditure. For example, use `ivexp=["income"]` if household income is your chosen instrument. The package estimates a log-expenditure regression and adds its residual to the demand equations and, when `censor=True`, the participation Probits. Works with either value of `censor`. Leave as `None` to skip this correction. |
| `anot` | Number | *Required* | The fixed constant $\alpha_0$ in the price index, for example `anot=8.0`. You supply this value; the package does not estimate it. It is not the regression intercept or a starting value. |

Rows with missing or nonfinite values in the required columns are left out.
Negative shares and nonpositive prices or expenditure in the remaining sample
cause an error. Externally supplied residual columns must contain only finite
values. See [Methodology](methodology.md) for the model equations and the role of `anot`.

---

### Model Specification

| Parameter | Type | Default | Description |
|---|---|---|---|
| `quadratic` | `bool` | `True` | `True` fits QUAIDS, including the squared log-expenditure term. `False` removes that term and fits AIDS. This corresponds to `noquadratic` in Stata. |
| `censor` | `bool` | `True` | `True` fits a participation Probit for each good and applies the Shonkwiler-Yen correction for zero shares. `False` fits the demand system without this correction and estimates one fewer equation because shares add up to one. This corresponds to `nocensor` in Stata. |

---

### Control Function and Selection Design

| Parameter | Type | Default | Description |
|---|---|---|---|
| `control_function` | Column name | `None` | A residual you have calculated outside the package, for example `"vhat"`. Adds that residual to each demand equation with a separate coefficient (`cfcoef_i`). Leave as `None` to omit it. For the built-in expenditure correction, use `ivexp` instead. |
| `selection_prices` | List of column names | `None` | Prices used in the participation Probits. This is useful when you want your probit equation to have different prices than the ones in the demand stage. `None` uses all demand prices; a list such as `["p1", "p3"]` uses only those prices; `[]` uses no prices. Names must come from your `prices` or `lnprices` list. This does not change the prices in the demand equations. |
| `selection_expenditure` | `bool` or `None` | `True` | `True` includes log expenditure in the participation Probits; `False` leaves it out. `None` includes it when you supply `expenditure` and leaves it out when you supply `lnexpenditure`. Use `True` or `False` to choose explicitly. |
| `selection_covariates` | List of column names | `None` | Other variables used in the participation Probits.This is useful when you want your probit equation to have different(more or less) demographic variables than the ones in the demand stage. `None` uses the variables in `demographics`; a list uses the columns you name; `[]` uses none. These choices do not change `demographics` in the demand equations. |
| `selection_control_function` | Column name | `None` | An externally calculated residual to include in the participation Probits. Each Probit estimates its own coefficient. To use the same residual in both stages, give its column name to both this option and `control_function`. |

`ivexp` and `control_function` work with or without censoring correction.
Without censoring, the control-function coefficients sum to zero so predicted
shares still add to one. The `selection_*` options require `censor=True`:
they configure the participation Probits, which model whether a household buys
each good.

Use either `ivexp` or the external residual options, not both. Variables named
in `ivexp` must not also be shares, prices, expenditure, demographics, or
`selection_covariates`. With `ivexp`, the package can run bootstrap and
re-estimates the expenditure regression in each replication. With an external
residual, use `reps=0`: the package does not know how to rebuild that residual.

---

### Estimation Methods and Optimizer Settings

| Parameter | Type | Default | Description |
|---|---|---|---|
| `method` | `str` | `"ifgnls"` | `"nls"` fits the equations using equally weighted squared errors. `"fgnls"` first runs NLS, then uses the estimated error variances and correlations to fit again. `"ifgnls"` repeats these updates until the estimates converge. These mean nonlinear least squares, feasible generalized nonlinear least squares, and its iterated version. |
| `algorithm` | `str` | `"gn"` | The numerical solver. `"gn"` uses Gauss-Newton and reduces the step size when needed. `"lm"` uses Levenberg-Marquardt, which adjusts the step through damping. Both fit the model selected by `method`. |
| `start` | `str` | `"zero"` | `"zero"` starts all free parameters at zero. `"linear"` calculates starting values from a simpler AIDS fit. If you supply `initial`, that vector is used instead. |
| `initial` | Numeric vector | `None` | Your own starting values for the free demand parameters. To restart the same model from a previous fit, use `initial=res.theta`. Do not use `res.b`, which also includes derived parameters, Probit coefficients, and elasticities. The vector must have `res.spec.n_free` entries in the same parameter order. |
| `sigma_initial` | Numeric matrix | `None` | Starting error covariance matrix for IFGNLS. Use `sigma_initial=res.sigma` together with `initial=res.theta` from the same model. Its shape is $n\times n$ with censoring and $(n-1)\times(n-1)$ without censoring, where $n$ is the number of goods. |

---

### Variance-Covariance Convention

| Parameter | Type | Default | Description |
|---|---|---|---|
| `vce_sigma` | `str` | `"objective"` | Chooses the error covariance used to calculate analytical standard errors after FGNLS/IFGNLS. `"objective"` uses the matrix from the last model-fitting step. `"final"` recalculates it from the final residuals. It does not change the coefficient estimates or replace bootstrap. |

---

### Convergence and Numerical Tolerance

| Parameter | Type | Default | Description |
|---|---|---|---|
| `param_tol` | Positive number | `1e-5` | Stops an inner fit when the relative change in its parameter values is small enough. A smaller value requires a smaller change. |
| `objective_tol` | Positive number | `1e-7` | Stops an inner fit when the relative improvement in its weighted sum of squared errors is small enough. |
| `gn_tol` | Positive number | `1e-5` | Stops an inner fit when the scaled Gauss-Newton criterion is small enough. This checks how much improvement the current search direction suggests. |
| `outer_param_tol` | Positive number | `1e-5` | For IFGNLS, sets how little the parameters must change between successive updates of the error covariance matrix. The criterion must pass twice in a row. |
| `max_iter` | Integer, at least 1 | `300` | Maximum solver iterations within each estimation stage. Raise it if a fit stops because it reaches this limit. |
| `max_outer` | Integer, at least 2 | `200` | Maximum numbered IFGNLS stage. The count includes the initial NLS and FGNLS stages. |
| `chunk` | Integer, at least 1 | `2000` | Number of observations processed together when calculating derivatives. A smaller value uses less working memory; a larger value may reduce processing overhead. It does not change the model. |


See [Methodology](methodology.md#convergence) for the exact formulas.

---

### Bootstrap and Parallelism

| Parameter | Type | Default | Description |
|---|---|---|---|
| `reps` | Nonnegative integer | `0` | Number of bootstrap replications to attempt. `0` skips bootstrap; `200` requests 200 replications. `1` attempts one replication, but the package needs at least two successful replications to calculate standard errors and otherwise raises an error. |
| `seed` | Integer or `None` | `None` | Sets the random seed for bootstrap samples, for example `seed=12345`. Reusing it with the same data and settings gives the same resamples. `None` leaves the seed unset. |
| `n_jobs` | Positive integer | `1` | Number of processes used for parallelized bootstrapping. It is useful when running a high number of `reps` and you want to utilize more of your machine's processing power to finish faster. I recommend using at most one less than the number of logical cores on your system to avoid freezing. This only helps when `reps>0`; it does not split a single fit across workers. |
| `blas_threads` | Positive integer or `None` | `1` | Threads used for matrix calculations within each estimation process. `1` limits each process to one thread; a larger number allows more. `None` keeps your existing setting. This is separate from `n_jobs`, which controls bootstrap workers. |
| `bootstrap_start` | `str` | `"zero"` | `"zero"` starts each bootstrap fit at zero. `"warm"` starts it from the estimates and error covariance of the full-sample fit, which can reduce the work needed in each replication. |
| `mp_context` | `str` or `None` | `None` | How Python starts bootstrap processes. `None` uses `"spawn"`, which works across supported operating systems. You can request `"spawn"` directly; `"forkserver"` and `"fork"` are also accepted where your operating system supports them. Most users can leave the default. |
| `rep_timeout` | Positive number or `None` | `None` | Maximum seconds allowed for each bootstrap replication. If it takes longer, that replication is counted as unsuccessful and the program moves on to the next one. For example, `900` allows 15 minutes per replication. `None` sets no time limit. |

If you use parallel bootstrap (`n_jobs > 1`) in a Python script, put your estimation code inside `main()` and add this at the end of the file:

```python
if __name__ == "__main__":
    main()
```

---

### Logging

| Parameter | Type | Default | Description |
|---|---|---|---|
| `verbose` | `bool` | `True` | `True` prints model-fitting and bootstrap progress. `False` hides these messages. You can still print the final results yourself. |
| `gn_verbose` | `bool` | `False` | `True` prints each step of the numerical solver. `False` hides these detailed messages. This setting is separate from `verbose`. |
| `log` | Callable or `None` | `None` | Optional Python function that receives progress messages instead of the usual output, for example `log=messages.append` to collect them in a list. It overrides `verbose` for those messages. Leave as `None` for normal use. |

---

## 2. Result Object Reference: `QuaidsceResults`

`res = quaidsce(...)` saves the fitted model in `res`. Use the fields and methods below to read its estimates, standard errors, elasticities, and model information.

### Overview Table: Python vs. Stata Equivalents

| Attribute / Property | Python Type | Stata Equivalent | Description |
|---|---|---|---|
| `res.b` | `np.ndarray` | `e(b)` | All reported estimates in one vector: demand parameters, participation Probit coefficients, and elasticities. Use `res.names` to identify each entry. |
| `res.V` | `np.ndarray` | `e(V)` | Covariance matrix for the reported estimates. After bootstrap, this is the bootstrap matrix. Without bootstrap, it is the analytical matrix, conditional on the fitted first stage; its elasticity covariance is not implemented. |
| `res.se` | `np.ndarray` | standard errors | Standard errors in the same order as `res.b`. After bootstrap, these are bootstrap standard errors. Without bootstrap, they are conditional analytical standard errors; elasticity standard errors are `NaN` (not available). |
| `res.V_analytic` | `np.ndarray` or `None` | — | The analytical covariance saved separately when bootstrap replaces `res.V`. `None` means no separate copy was needed. |
| `res.analytic_se` | `np.ndarray` | — | Analytical standard errors, still available after bootstrap. They treat the first-stage estimates as fixed. Elasticity standard errors are `NaN`; use bootstrap for those. |
| `res.names` | `List[str]` | `colnames e(b)` | Labels for `res.b`, such as `"beta:beta_1"`. The same order is used for standard errors and covariance rows and columns. |
| `res.theta` (or `res.b_est`) | `np.ndarray` | `e(b_est)` | Free demand parameters estimated by the solver. Use this vector with `initial` to restart the same model. |
| `res.V_est` | `np.ndarray` | — | Analytical covariance of the free demand parameters. This remains analytical after bootstrap; use `res.V` for the bootstrap covariance of all reported estimates. |
| `res.coefs` | `Coefs` dataclass | Demand entries in `e(b)` | Demand coefficients grouped by their role: $\alpha, \beta, \gamma, \lambda, \delta, \eta, \rho$, and `cfcoef` when a demand control function is used. |
| `res.tau` | `np.ndarray` | Probit entries in `e(b)` | Participation Probit coefficients, stored one good after another. |
| `res.setau` | `np.ndarray` | — | Analytical covariance of the participation Probit coefficients. Each Probit has its own block; covariance between different Probits is not included here. |
| `res.probits` | `List[ProbitResult]` | — | A list containing the separate participation Probit results. |
| `res.elas` | `Elasticities` | `e(elas_*)` | Expenditure and price elasticities at sample means. In each price matrix, row `i` is the good demanded and column `j` is the price that changes. |
| `res.llf` | `float` | `e(ll)` | Gaussian log-likelihood value. |
| `res.sigma` | `np.ndarray` | `e(Sigma)` | Residual covariance matrix ($(n \times n)$ or $(n-1 \times n-1)$). |
| `res.nobs` | `int` | `e(N)` | Number of observations in the estimation sample. |
| `res.converged` | `bool` | `e(converged)` | `True` means the model met its convergence criteria. `False` means it stopped before meeting them. |
| `res.boot` | `BootResult` or `None` | — | Bootstrap results, including successful estimates (`b_star`), standard errors (`se`), successful count (`reps_ok`), requested count (`reps_requested`), and failure messages (`failures`). `None` means bootstrap was not run. |
| `res.notes` | `List[str]` | — | Messages about the fitted model, including convergence and standard errors. |
| `res.ivexp_names` | `List[str]` | `e(ivexp)` | Excluded expenditure instruments used by the internal reduced form. |
| `res.reduced_form` | `ExpenditureReducedForm` or `None` | `e(reduced_form_b)`, `e(reduced_form_V)` | Internal log-expenditure OLS result, including residuals, fitted values, R-squared, and the joint excluded-instrument F test. |

---

### Detailed Access to Structural Coefficients (`res.coefs`)

The structural coefficients are accessible directly as NumPy arrays:

```python
# Individual structural parameter blocks
alpha = res.coefs.alpha   # (n,) array: intercepts
beta  = res.coefs.beta    # (n,) array: linear expenditure coefficients
gamma = res.coefs.gamma   # (n, n) symmetric array: price interaction terms
lam   = res.coefs.lam     # (n,) array: quadratic expenditure terms
delta = res.coefs.delta   # (n,) array: Shonkwiler-Yen censoring terms
eta   = res.coefs.eta     # (R, n) array: demographic interaction terms
rho   = res.coefs.rho     # (R,) array: Ray demographic scaling factors
cfcoef = res.coefs.cfcoef # (n,) array: demand control-function coefficients
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

- **`res.summary(level=95.0)`**: Returns a text table of model information, estimates, standard errors, and 95% confidence intervals. Use `print(res.summary())` to display it, or change `level` for another confidence level.
- **`res.elasticity_tables()`**: Returns formatted summary tables for expenditure, Marshallian, and Hicksian elasticities.
- **`res.reduced_form_table()`**: Returns the internal log-expenditure OLS table and excluded-instrument diagnostic when `ivexp` is active.
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

---

## 4. R Package Reference (`pyquaidsce`)

The R package calls the same Python estimator through `reticulate`. Use `TRUE`, `FALSE`, and `NULL` in place of Python's `True`, `False`, and `None`. Give variable lists as `c("w1", "w2", "w3", "w4")`. The R defaults include `anot=10.0`, `verbose=FALSE`, and `mp_context="spawn"`.

```r
library(pyquaidsce)

# Estimate Censored QUAIDS in R
fit <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  anot = 10.0,
  method = "ifgnls"
)

# Standard S3 methods
summary(fit)
coef(fit)
vcov(fit)

# Custom starting values
fit_warm <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  anot = 10.0,
  initial = fit$theta,
  sigma_initial = fit$sigma
)
```

---

## 5. Stata Interface (`quaidsce` and `pyquaidsce`)

You can use either `quaidsce` or `pyquaidsce` in Stata. Both commands run the same estimator. If you already have the older Stata `quaidsce` installed and want to avoid a name conflict, use `pyquaidsce`.

```stata
* Estimate in Stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) anot(10.0) method(ifgnls)

* Equivalent direct command for the current package
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) anot(10.0) method(ifgnls)

* Warm-start using stored matrices
matrix b_init = e(b_est)
matrix sigma_init = e(Sigma)

quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) anot(10.0) initial(b_init) sigma_initial(sigma_init)
```
