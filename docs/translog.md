# Basic translog demand with optional Shonkwiler--Yen correction

This is an unreleased addition to the supplied 1.6.0 source. Install this source
with `python -m pip install -e .` before using `translog()`. The existing
`quaidsce()` interface, R wrappers, and Stata commands retain their behavior;
the new entry point currently belongs to the Python API.

```python
from pyquaidsce import translog

fit = translog(
    df,
    shares=["w1", "w2", "w3"],
    prices=["p1", "p2", "p3"],
    expenditure="total",
    demographics=["hh_size", "age"],  # optional
)
print(fit.summary())
print(fit.elasticity_tables())
shares_hat = fit.predict()
quantities_hat = fit.predict(new_data, kind="quantities")
```

All shares and corresponding prices must use the same good ordering. Supply
exactly one of `prices`/`lnprices` and `expenditure`/`lnexpenditure`. Zero
observed shares are allowed without a participation first stage. Negative shares
must be cleaned explicitly; the estimator never changes or renormalizes them.
The default `share_sum_tol=0.01` matches the supplied demandsys share check.

## Equations and parameter restrictions

For demographic vector \(z\), let

\[
c_i=\sum_r \nu_{ir}z_r,\qquad
\widetilde m=m-\sum_k p_kc_k,\qquad
r_j=\log(p_j/\widetilde m).
\]

The expenditure-share equations are

\[
w_i=\frac{p_ic_i}{m}+\frac{\widetilde m}{m}
\frac{\alpha_i+\sum_j\gamma_{ij}r_j}
     {1+\sum_k\sum_j\gamma_{kj}r_j}.
\]

The uncensored implementation imposes \(\sum_i\alpha_i=1\) and
\(\gamma_{ij}=\gamma_{ji}\). It does **not** impose zero Gamma row sums or
zero demographic coefficient sums. Homogeneity follows from the price and
expenditure ratios. This is Pollak--Wales demographic translation, distinct
from the Ray scaling used by the package's AIDS/QUAIDS models.

Without demographics, \(c=0\) and \(\widetilde m=m\). No independent
committed-quantity intercept is estimated: `gtranslog` is outside this phase.
Effective expenditure must be positive and the share denominator must not be
zero. The optimizer rejects trial steps outside that domain.

For \(G\) goods and \(R\) demographics, there are
\((G-1)+G(G+1)/2+GR\) free parameters. The free vector contains the first
\(G-1\) alphas, the symmetric Gamma upper triangle in row order, then Nu
columns in demographic order. The last alpha is recovered by adding-up.
For the Uruguay specification, that is 216 free parameters and 217 reported
structural coefficients.

## Estimation, results, and inference

The default is `method="ifgnls", start="zero", algorithm="gn",
stop_rule="standard"`. `nls` and `fgnls` use the same shared engine.
`start="mean"`, custom `initial`, `sigma_initial`, LM, and convergence controls
are also available. Nonlinear optimization paths can depend on the starts and
stopping rules; the Uruguay comparison records its exact defaults.

The last share equation is omitted from SUR to avoid a singular disturbance
covariance. Predictions still contain all goods. Estimation uses complete rows
of shares, prices, expenditure, and demographics; `fit.estimation_mask` records
that selection. Predictions on a supplied frame retain its row order and
return NaNs where required inputs are missing. Residuals use observed minus
predicted shares.

Useful result fields include:

- `coefs.alpha`, `coefs.gamma`, `coefs.nu`: model-specific coefficient arrays.
- `theta`, `b_est`, `V_est`: free structural estimates and their covariance.
- `b`, `names`, `V`, `se`: full coefficients plus derived elasticities and
  their joint covariance/standard errors.
- `sigma`, `llf`, `r2`, `converged`, `n_outer`, `n_gn`: estimation diagnostics.
- `elas.income`, `elas.uncompensated`, `elas.compensated`: elasticity arrays.
- `elasticity_se`, `elasticity_covariance`, `coefficient_covariance`:
  inference blocks. Price matrices always have rows = goods, columns = prices.

Analytic standard errors differentiate the **complete** elasticity formula,
including the predicted-share denominator, and propagate `V_est` through the
delta method. Parameter derivatives of shares and price/expenditure derivatives
are analytic. The result transformation uses complex-step differentiation,
avoiding finite-difference subtraction error. Evaluation means are treated as
fixed in this analytical inference.

By default elasticities are evaluated at arithmetic means of price and
expenditure levels on the estimation sample, using predicted shares at that
point. They are not averages of observation-level elasticities and do not use
mean observed shares as denominators.

```python
# Explicit postestimation on another frame, for example the full Uruguay data.
elas, se = fit.elasticities_at_means(df, standard_errors=True)
```

For an explicitly supplied frame with missing inputs, prices use rows with a
complete price vector, demographics use rows with a complete demographic
vector, and expenditure uses its available rows. This reproduces the supplied
Uruguay postestimation mean convention. It need not use the same rows as SUR.

## Shared bootstrap

```python
if __name__ == "__main__":  # required for spawn-based scripts
    fit = translog(df, shares=["w1", "w2", "w3"], prices=["p1", "p2", "p3"],
                   expenditure="total", reps=200, seed=123, n_jobs=4,
                   bootstrap_start="warm", mp_context="spawn")
```

Bootstrap jobs reuse the package's resampling, deterministic seeds, worker
supervision, timeouts, and temporary BLAS limits. Each replication re-estimates
the model and its evaluation means. `V` and `se` then describe the bootstrap;
`V_analytic` and `analytic_se` retain the analytical alternative. The serial and
spawn paths are tested with the same seeds.

The model backend separates equations, parameter restrictions, derivatives,
and starting values from SUR. Censoring uses the shared participation and S&Y
composition modules; no expenditure endogeneity correction is exposed.

## Shonkwiler--Yen translog

```python
fit = translog(
    df, shares=["w1", "w2", "w3"], prices=["p1", "p2", "p3"],
    expenditure="total", demographics=["hh_size", "age"],
    censor=True, start="nested",
)
latent = fit.predict(kind="latent_shares")
observed_mean = fit.predict(kind="shares")
participation = fit.predict(kind="participation")
```

S&Y estimates **every** equation and frees **every** alpha by default. Set
`latent_adding=True` to recover the last alpha as `1 - sum(alpha[:-1])`.
It keeps the denominator constant one, Gamma symmetry and the same
demographic translation. Zero Gamma or demographic sums are not required:
sum(alpha)=1 already makes the numerator sum equal to the denominator and
the translated latent shares sum to one. See [latent adding-up](latent-adding.md).
Enabling `latent_adding` without censoring raises an error; the uncensored
restrictions and mandatory adding-up are preserved.
Observed zero shares remain in estimation. Each estimated Probit needs both
positive and zero shares. Demographics are optional for translog censoring.
Observed and predicted shares are not normalized after estimation.

Let \(f_i\) be the latent translated share above. The participation indicator
is \(1(w_i>0)\), with Probit linear predictor \(k_i=x'\tau_i\). Then

\[
W_i=\Phi(k_i)f_i+\delta_i\phi(k_i).
\]

The default Probit uses all log prices, log expenditure and demand demographics.
`selection_prices`, `selection_expenditure` and `selection_covariates` specify
an independent design. `selection_prices=[]` or `selection_covariates=[]`
omits that group; `selection_expenditure=False` omits expenditure. Selection
covariates can differ from translation demographics. Their missing values
participate in estimation-sample selection and observed-mean prediction.
Expenditure is included by default for both level and log input APIs.

For \(G\) goods and \(R\) demographics the structural vector is
`[alpha, Gamma upper triangle, Nu columns, delta]`, with
\(G+G(G+1)/2+GR+G\) free parameters. Probit coefficients are reported after
these structural coefficients and before elasticities. Results expose `tau`,
`setau`, `probits`, `selection_layout`, `selection_index`, `n_coefficients`
and `initialization`; `theta` and `V_est` remain structural quantities.

Define \(t_i=p_ic_i/m\), \(b=1-\sum_i t_i\), \(s_i=N_i/D\),
\(g_i=\sum_j\gamma_{ij}\), \(H=\sum_i g_i\), and
\(T_i=(g_i-s_iH)/D\). The latent log-input derivatives are

\[
F_{i,m}=-t_i+(1-b)s_i-T_i,
\qquad
F_{i,j}=\mathbf1_{i=j}t_i-s_it_j
+\frac bD(\gamma_{ij}-s_ig_j)+T_it_j.
\]

The corrected derivatives and elasticities are

\[
M_i=\Phi(k_i)F_{i,m}+\phi(k_i)\tau_{i,m}(f_i-\delta_i k_i),
\]

\[
B_{ij}=\Phi(k_i)F_{i,j}+\phi(k_i)\tau_{i,p_j}(f_i-\delta_i k_i),
\]

\[
\eta_i=1+M_i/W_i,\qquad
\epsilon^M_{ij}=-\mathbf1_{i=j}+B_{ij}/W_i,\qquad
\epsilon^C_{ij}=\epsilon^M_{ij}+\eta_iW_j.
\]

Selection derivatives for omitted inputs are zero. Every denominator uses the
predicted S&Y share at that evaluation point. The Probit index/CDF/PDF are
recomputed at the level means; means of the individual CDFs are not substituted.
Negative predictions are reported without clipping; zero predicted shares
have undefined elasticities. S&Y compensated results are labeled a **Slutsky
transformation**, because utility-consistent Hicksian demand is not guaranteed
with relaxed latent adding-up and unrestricted selection equations.

### Starting values and the positive-expenditure domain

Defaults remain `method="ifgnls", algorithm="gn", stop_rule="standard",
start="zero"`, matching the standard package engine. `start="mean"` and
explicit `initial`/`sigma_initial` also remain available. For S&Y translation,
`start="nested"` fits the no-translation submodel with the **same** participation
terms and SUR settings, lifts its coefficients with initial Nu equal to zero,
and estimates the complete model. Demographics remain free in the final fit.
Explicit initial values take precedence; no nested fit is then performed.

On Uruguay, zero/mean starts can produce tiny accepted steps near
\(m-\sum_i p_ic_i=0\). A small step alone is not stationarity. `nrtol` records
the final scaled gradient. S&Y translog flags `converged=False` when the fit
is near the effective-expenditure/denominator boundary and that gradient exceeds
1e-5; the standard stopping rules otherwise keep their behavior. The nested
start obtains an interior, stationary fit
with the original GN/IFGNLS settings, without changing data or restrictions.

### Two-step inference

Analytical structural and elasticity SEs are **conditional on fitted Probits**
and fixed evaluation means. The displayed Probit covariance uses observed
information. These analytical blocks are not a complete joint two-step
covariance: generated-regressor and cross-stage uncertainty require bootstrap.
Elasticity derivatives include the predicted-share denominator.

With `reps>=2`, each replication refits all Probits and SUR. `bootstrap_start="warm"`
passes structural/Sigma starting values, while still refitting participation.
`start="nested"` with cold bootstrap refits the initializer too. `V` and `se`
then use the complete bootstrap vector, including cross-stage covariance.
`V_analytic` preserves the conditional analytical alternative. A handful of
replications tests execution only; use enough successful replications for
statistical inference.

`elasticities_at_means(standard_errors=True)` uses the coefficient covariance
and a full delta map at fixed evaluation means. With bootstrap it differentiates
through Probit coefficients as well. Its SEs can differ from `elasticity_se`,
which is the empirical bootstrap SE of elasticities with means re-estimated
inside each replication.

## Stata comparison and references

Run `python benchmarks/translog_uruguay/run_python.py` for the independent
zero-start comparison. The supplied `_ds_translognomata.ado` is an older
alpha-shifting implementation and is not the reference for this translated
specification. The newer wrapper calls `_gtranslog_wrk()`; the equations above
match the official manual and reproduce the reference likelihood and R-squared.

The reported Stata elasticity SEs are numerically reproduced by **holding the
predicted-share denominator fixed** during delta-method propagation. That
conditional calculation is included as a benchmark diagnostic. Public package
SEs use the complete formula and therefore differ. The comparison CSV records
both definitions; it does not substitute the diagnostic into `fit.se`.

- [Stata demandsys manual](https://www.stata.com/manuals/rdemandsys.pdf),
  basic translog and demographic translation, methods and formulas.
- Christensen, Jorgenson, and Lau (1975), *Transcendental Logarithmic Utility Functions*.
- Pollak and Wales (1992), *Demand System Specification and Estimation*.
