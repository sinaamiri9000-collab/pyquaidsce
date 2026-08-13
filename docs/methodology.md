# Methodology and estimator

## Model

The core model is the quadratic almost-ideal demand system (QUAIDS) of Banks,
Blundell and Lewbel (1997), with Ray (1983) demographic scaling in the form used
by Poi (2012).

Let prices be `p`, total system expenditure be `m`, and demographic variables be
`z`. The implementation uses

```text
ln a(p)  = alpha_0 + sum_i alpha_i ln p_i
           + 1/2 sum_i sum_j gamma_ij ln p_i ln p_j
b(p)     = exp(sum_i beta_i ln p_i)
mbar(z)  = 1 + sum_r rho_r z_r
c(p,z)   = exp(sum_i (sum_r eta_ri z_r) ln p_i)
D        = ln m - ln a(p) - ln mbar(z)

w*_i     = alpha_i + sum_j gamma_ij ln p_j
           + (beta_i + sum_r eta_ri z_r) D
           + lambda_i / (b(p)c(p,z)) D^2
```

For censored budget shares, the package applies the Shonkwiler–Yen (1999)
two-step construction. A participation probit is estimated for each good and
the observed-share equation is formed from a normal-CDF multiplier and a
normal-PDF correction term.

## Demand-theory restrictions

The implementation imposes symmetry and homogeneity on the `gamma` matrix. The
Ray demographic coefficients are also completed using their restriction. In the
uncensored model, adding-up restrictions determine the final `alpha`, `beta`,
and `lambda` terms. In the censored branch, the transformed share equations do
not retain the same adding-up identity, so all transformed equations remain in
the estimated system, matching the architecture of `quaidsce`.

## Estimation sequence

### NLS

The nonlinear system is first fit using an identity working covariance. This
provides residuals and an initial estimate of the cross-equation covariance
matrix.

### FGNLS

The residual covariance is used to reweight the nonlinear SUR criterion and the
parameters are re-estimated.

### IFGNLS

With `method="ifgnls"`, the covariance matrix and nonlinear parameters are
updated repeatedly until the outer fixed point satisfies the convergence rule.
Consequently, NLS and FGNLS messages during an IFGNLS run represent required
initial stages, not separate final models.

## Numerical algorithm

The primary optimizer is Gauss–Newton with Hartley step halving and an analytic
Jacobian. A Levenberg–Marquardt alternative exists for difficult numerical
problems. Normal equations are accumulated in observation blocks so the full
stacked Jacobian does not need to be materialized in memory.

`start="zero"` matches the Stata starting convention. `start="linear"` uses a
linearized AIDS starting fit and can be faster, but nonlinear censored demand
systems can have multiple local solutions, so starting values are part of the
reported empirical specification.

## Covariance and likelihood

The residual covariance uses divisor `N` in the compatibility path. The
reported Gaussian-system log likelihood is

```text
-(N*m/2)(1 + ln(2*pi)) - (N/2) ln|Sigma_hat|.
```

The structural covariance is based on the inverse cross-product of the weighted
analytic Jacobian and is then mapped from free to reported parameters by the
restriction Jacobian. First-stage probit covariance blocks are appended in the
Stata-oriented result representation.

## Elasticities

The result object contains expenditure, uncompensated (Marshallian), and
compensated (Hicksian) elasticities. The natural Python matrices use
`[good, price]` indexing. A Stata-order vector is available for direct
comparison with the original command's storage convention.

Because compatibility switches can alter first-stage correction terms or a
published elasticity formula, record `first_stage_predict` and `strict_stata`
with every reported empirical result.

## References

- Banks, J., Blundell, R. and Lewbel, A. (1997). “Quadratic Engel Curves and
  Consumer Demand.” *Review of Economics and Statistics*, 79, 527–539.
- Deaton, A. and Muellbauer, J. (1980). “An Almost Ideal Demand System.”
  *American Economic Review*, 70, 312–326.
- Poi, B. P. (2012). “Easy Demand-System Estimation with quaids.” *Stata
  Journal*, 12, 433–446.
- Ray, R. (1983). “Measuring the Costs of Children: An Alternative Approach.”
  *Journal of Public Economics*, 22, 89–102.
- Shonkwiler, J. S. and Yen, S. T. (1999). “Two-Step Estimation of a Censored
  System of Equations.” *American Journal of Agricultural Economics*, 81,
  972–982.
- Caro, J. C., Melo, G., Molina, J. A. and Salgado, J. C. (2021). “Censored
  QUAIDS estimation with quaidsce.” Boston College Working Papers in Economics
  1045.
