# Validation and Numerical Evidence

`pyquaidsce` v1.7.0 is validated against econometric identities, numerical derivatives, cross-interface contracts, and archived Stata comparisons. Exact reproduction of known mistakes in the original Stata ado is **not** a v1.7.0 objective.

## 1. Current v1.7.0 validation policy

The current implementation uses two corrected behaviors:

- the Shonkwiler-Yen first stage always uses the Probit linear index, $\Phi(X'\tau)$ and $\phi(X'\tau)$;
- elasticities always use the corrected theoretical formulas.

Therefore, differences from the original ado are expected in specifications affected by its nested-Probit transformation or documented elasticity indexing/macro mistakes. See [Stata Compatibility](stata-compatibility.md).

## 2. Automated Python tests

The repository contains automatic Python tests covering:

- analytic Jacobians against finite differences;
- adding-up, homogeneity, symmetry, and elasticity identities at model-consistent shares;
- Probit estimation and information matrices;
- input validation and parameter mapping;
- control-function and custom-selection extensions;
- integrated `ivexp` reduced forms, excluded-instrument F tests, and bootstrap rebuilding;
- bootstrap covariance synchronization and timeout behavior;
- Stata bridge argument/matrix forwarding;
- verification that the public default estimation method is IFGNLS.

Run them with:

```bash
python -m unittest discover -s tests -v
```

## 3. Comparison with Stata

A four-good censored QUAIDS comparison used 2,665 households, IFGNLS, and
200 bootstrap replications with 8 workers in each program.

| Measure | pyquaidsce | quaidsce_c |
|---|---:|---:|
| Time (seconds) | 122.05 | 10,810.18 |
| Bootstrap replications | 200 | 200 |

pyquaidsce was about **88.6 times faster**. The code, data, and results are in
[benchmark/pyquaidsce vs quaidsce](../benchmark/pyquaidsce%20vs%20quaidsce/).
