# Stata Compatibility and Validation Guide

`pyquaidsce` retains a Stata-facing interface, coefficient ordering, and output layout so results can be compared conveniently with the original `quaidsce` command. In **v1.7.0**, however, the package no longer offers switches that intentionally reproduce known econometric mistakes in the original ado-file. Python, R, and Stata interfaces now use one canonical implementation.

## Canonical first-stage censoring correction

For the Shonkwiler & Yen (1999) two-step procedure, each participation Probit supplies the linear index $X_i'\tau_i$. The package always computes

$$\Phi_i = \Phi(X_i'\tau_i), \qquad \phi_i = \phi(X_i'\tau_i).$$

The original Stata ado calls `predict` after `probit` without requesting `xb`; because Stata's default prediction is a probability, that code can feed $\Phi(X'\tau)$ back into `normal()` and `normalden()`. `pyquaidsce` intentionally does **not** reproduce that nested transformation.

## Canonical elasticity formulas

`pyquaidsce` also uses the published/theoretically intended elasticity formulas in all specifications. In particular, it does not reproduce two identified ado-file edge-case errors: the no-demographics quadratic Marshallian term that uses the wrong beta index, and the demographics + linear-AIDS censoring branch where a global/local macro mismatch can zero the latent expenditure elasticity.

These choices are fixed behavior in v1.7.0+. There is no compatibility switch in any public interface.

## Comparing pyquaidsce with Stata

### Cross-price names on the publication branch

`ELAS_UNCOMP:e_i_j` and `ELAS_COMP:e_i_j` identify demand for good **i**
with respect to price **j**, matching matrix element `[i-1, j-1]` in Python
and `[i, j]` in R / Stata. Price matrices are flattened row by row.

Before this fix, the same numeric vector used transposed cross-price labels.
The fix changes names only: numerical coefficients, covariance positions,
bootstrap draws, and standalone elasticity matrices are unchanged. Code that
looks up a cross-price elasticity by name now receives the correct direction.

Updating the package does not rewrite stored results. For a verified affected
censored Python result with the standard final elasticity block, rebuild only
that block's names using the corrected package:

```python
n = res.spec.neqn
k = n + 2 * n * n
assert res.spec.censor and len(res.names) == len(res.b)
res.names = list(res.names[:-k]) + res.spec.elas_names()
```

For saved R / Stata output, update coefficient names and both covariance axes
together; also update names attached to bootstrap standard errors. Do not
transpose numeric matrices or reorder bootstrap columns. Archived validation
logs remain unchanged; their readers explicitly normalize known legacy names.

One bootstrap replication can check a fitted result and its labels, but cannot
estimate bootstrap standard errors. The public bootstrap requires at least two
successful replications. `tools/check_elasticity_benchmark.py` performs the
original fit and exactly one resample separately for this diagnostic.

When comparing results, check the following before attributing a difference to the implementation:

1. **Estimation sample** — missing values and positivity filters must select the same observations.
2. **Variable ordering** — shares and prices must use identical good ordering.
3. **Expenditure definition** — use the same level/log definition and same system expenditure.
4. **Starting values** — match `start`, `initial`, and `sigma_initial` where relevant.
5. **Estimation method** — `pyquaidsce` defaults to **IFGNLS**; explicitly match `ifgnls`, `fgnls`, or `nls` in the comparison program.
6. **Known ado deviations** — exact equality is not expected in cases where the original ado uses the nested-Probit transformation or one of the elasticity mistakes described above.

## Current defaults

| Setting | Default |
|---|---|
| Censoring correction | Shonkwiler–Yen using the Probit linear index |
| Elasticity formulas | Corrected theoretical formulas |
| `method` | `"ifgnls"` |
| `start` | `"zero"` |
| `algorithm` | `"gn"` |
| `vce_sigma` | `"objective"` |

## R example

```r
library(pyquaidsce)

fit <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  anot = 10.0
)
```

## Stata example

You can use either `quaidsce` or `pyquaidsce`. Both run the same estimator. If you already have the older Stata `quaidsce` installed and want to avoid a name conflict, use `pyquaidsce`.

```stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) anot(10.0)

* Equivalent direct current-package command
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) anot(10.0)
```

Both examples use IFGNLS by default and the same canonical first-stage and elasticity definitions.
