# Validation and Numerical Evidence

This page describes the automated tests and the comparison with Stata.

## 1. Automated Python tests

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

## 2. Comparison with Stata

A four-good censored QUAIDS comparison used 2,665 households, IFGNLS, and
200 bootstrap replications with 8 workers in each program.

| Measure | pyquaidsce | quaidsce |
|---|---:|---:|
| Time (seconds) | 122.05 | 10,810.18 |
| Bootstrap replications | 200 | 200 |

pyquaidsce was about **88.6 times faster**. The code, data, and results are in
[benchmark/pyquaidsce vs quaidsce](../benchmark/pyquaidsce%20vs%20quaidsce/).
