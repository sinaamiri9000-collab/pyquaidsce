# `ivexp` Implementation Report

## Outcome

Version 1.5.0 adds one integrated endogenous-expenditure workflow across all
three public interfaces:

| Interface | Syntax |
|---|---|
| Python | `ivexp=["log_income", "assets"]` |
| R | `ivexp = c("log_income", "assets")` |
| Stata | `ivexp(log_income assets)` |

The existing external `control_function` and
`selection_control_function` paths remain available and backward compatible.

## Econometric contract

If the user supplies expenditure in levels, the package first computes
`ln(expenditure)`; if `lnexpenditure` is supplied, it uses that column directly.
On the final common complete-case estimation sample, it estimates

$$
\ln m_h = \pi_0 + \pi_p'\ln p_h + \pi_z'z_h + \pi_q'q_h + v_h
$$

by OLS, where $q_h$ contains the excluded variables supplied through `ivexp`.
All demand prices and all Ray demographic variables are included as exogenous
controls, together with an intercept.

The generated residual enters the participation equation for every good,

$$
d_{ih}=\mathbb{I}(X_{ih}'\tau_i+\psi_i v_h+u_{ih}>0),
$$

and the latent demand share,

$$
\mu_{ih}=\Phi_{ih}(w^Q_{ih}+\kappa_i v_h)+\delta_i\phi_{ih}.
$$

Thus $\psi_i$ and $\kappa_i$ are distinct, unrestricted, equation-specific
coefficients. This reuses the package's tested control-function parameter map,
Jacobians, Shonkwiler-Yen transformation, and elasticity derivatives.

## Validation rules

- `ivexp` requires at least one unique column.
- Its columns must not overlap shares, demand prices, Ray demographics,
  expenditure, or `selection_covariates`.
- It cannot be combined with either precomputed control-function argument.
- It requires `censor=True` and `first_stage_predict="xb"`.
- The reduced-form matrix must have positive residual degrees of freedom, full
  column rank, finite values, and a nondegenerate residual.
- The package never mutates the caller's DataFrame or creates a visible
  temporary residual column.

## Inference and elasticities

The conditional analytical covariance remains available as a reference, but
it does not account for the generated residual. With `reps>0`, every bootstrap
draw re-estimates the expenditure reduced form before fitting its Probits and
demand system. This is the supported generated-regressor inference path.

Externally supplied residuals remain incompatible with the internal bootstrap,
because their generating equations are not known to the package.

Elasticity formulas are unchanged. The generated residual is held fixed when
computing structural price and expenditure derivatives. Since the OLS reduced
form contains an intercept, its in-sample residual mean is zero up to numerical
precision; the at-means control term is therefore centered naturally.

## Diagnostics and returned objects

The internal reduced form reports:

- coefficient vector and classical OLS covariance;
- fitted values and residuals;
- residual degrees of freedom and rank;
- $R^2$ and adjusted $R^2$;
- the classical joint F test for the excluded instruments.

The F statistic diagnoses relevance under the classical OLS assumptions. It
does not establish the exclusion restriction or instrument validity. The
residual-inclusion interpretation also requires a correctly specified
triangular reduced form and the distributional assumptions supporting the
augmented Probit.

| Interface | Returned diagnostics |
|---|---|
| Python | `res.reduced_form`, `res.reduced_form_table()` |
| R | `fit$reduced_form` |
| Stata | `e(reduced_form_b)`, `e(reduced_form_V)`, `e(reduced_form_r2)`, `e(excluded_iv_F)`, `e(excluded_iv_p)` |

## Verification completed

- All 38 Python unit, theory, regression, bootstrap, and bridge tests pass.
- The internally generated residual matches the same manually generated OLS
  residual, and both paths produce numerically equivalent Probit, structural,
  and elasticity estimates.
- A focused test counts one reduced-form fit for the point estimate plus one in
  every sequential bootstrap replication.
- The Python package compiles successfully.
- A version 1.5.0 wheel builds successfully and passes an isolated import/API
  smoke test.
- The version 1.5.0 source distribution includes the Python, R, and Stata
  interfaces, and all 38 Python tests pass from its extracted tree.
- R and Stata forwarding, option presence, asynchronous bridge transport, and
  stored-result plumbing are covered by repository tests/static checks. Native
  R and Stata executables were not available in the validation environment.
