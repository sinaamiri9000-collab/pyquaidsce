# Optional latent adding-up with Shonkwiler--Yen

Both Python estimators accept `latent_adding=False`. The default preserves the
existing censored parameterization. To impose adding-up before the S&Y
transformation, use:

```python
from pyquaidsce import quaidsce, translog

q = quaidsce(df, shares, prices=prices, expenditure="total",
             demographics=demos, anot=2.5,
             censor=True, latent_adding=True)
t = translog(df, shares, prices=prices, expenditure="total",
             demographics=demos, censor=True, latent_adding=True)
```

`latent_adding=True` requires `censor=True`. With `censor=False`, adding-up is
always imposed and cannot be disabled; leave `latent_adding` at its default.
The argument must be a boolean. The selected mode is stored in `result.spec`
and retained by all bootstrap refits, including warm starts and spawn workers.

## QUAIDS and AIDS

The last structural coefficients are recovered from

\[
\alpha_n=1-\sum_{i<n}\alpha_i,\qquad
\beta_n=-\sum_{i<n}\beta_i,\qquad
\lambda_n=-\sum_{i<n}\lambda_i.
\]

The lambda block is present only in QUAIDS. Existing Gamma symmetry and zero
row/column sums, and demographic eta sums of zero, are unchanged. If a demand
control function is active, its last coefficient is also recovered so that
`sum(cfcoef)=0`; the augmented latent shares then add up for every residual
value. Ray rho coefficients remain free.

## Basic indirect translog

Write translated shares as

\[
t_i=p_i c_i/m,\quad b=1-\sum_i t_i,\quad
f_i=t_i+b\frac{\alpha_i+(\Gamma r)_i}{1+\sum_i(\Gamma r)_i}.
\]

Imposing `sum(alpha)=1` makes the numerator sum equal to the denominator,
so `sum(f)=sum(t)+b=1`. Only the last alpha is recovered. Gamma remains
symmetric with free row sums, and all demographic translation coefficients
remain free. No generalized-translog intercept is introduced.

## Observed means, derivatives and inference

The observed mean remains

\[
W_i=\Phi(k_i)f_i+\delta_i\phi(k_i).
\]

Its sum is generally not one even when the latent shares add up. Every S&Y
equation is estimated, with unrestricted delta and Probit coefficients.
Neither observed means nor latent shares are normalized after estimation.
Adding-up does not impose positivity or curvature.

The free parameter vector becomes shorter: by three coefficients for QUAIDS,
two for AIDS, and one for translog. An active demand control function removes
one additional QUAIDS/AIDs coefficient. Custom `initial` vectors must match
the selected specification. The analytic Jacobian and covariance mapping
include the recovered last-good coefficients. S&Y elasticity formulas and
their existing evaluation conventions are unchanged.

Translog `start="nested"` retains the selected adding-up mode in its initial
submodel. Convergence still depends on data and starting values; imposing
adding-up does not remove the effective-expenditure domain restriction.

This option is currently exposed in the Python API. The existing R and Stata
entry points retain their defaults.
