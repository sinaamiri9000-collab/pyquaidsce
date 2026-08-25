# Numerical Validation against Stata quaidsce

To verify the numerical precision of `pyquaidsce`, we established a public, reproducible end-to-end benchmark comparing Python against Stata 19.5.

All benchmark data, scripts, and logs can be found in [`benchmarks/cquaids_ifgnls_4g_20k/`](../benchmarks/cquaids_ifgnls_4g_20k/).

---

## 1. Benchmark Specification

The validation benchmark is designed to thoroughly test the full censored QUAIDS estimator under realistic conditions:

- **Sample Size**: 20,000 observations
- **System Dimensions**: 4 goods, 3 demographic scaling variables
- **Censoring**: Zero budget shares present in every good (triggering 4 first-stage participation Probits and Shonkwiler–Yen transformations)
- **Estimator**: Iterated FGNLS (`method="ifgnls"`)
- **Initialization**: Default zero starting values in both programs
- **Dataset**: Identical synthetic `.dta` file read by both Stata and Python

---

## 2. Comparison of Results

A parameter-by-parameter comparison across all 113 reported quantities (structural coefficients, Probit parameters, standard errors, log-likelihood, and elasticities) yields the following differences:

| Parameter Category | Stata vs. Python Max Difference | Note |
|---|---:|---|
| Structural parameters ($\alpha, \beta, \gamma, \lambda, \delta, \eta, \rho$) | $< 1.68 \times 10^{-5}$ | Exact agreement across all coefficients |
| First-stage Probit parameters ($\tau$) | $< 1.21 \times 10^{-7}$ | Matches Stata's Newton-Raphson Probit |
| Income and Price Elasticities | $< 7.34 \times 10^{-7}$ | All Marshallian & Hicksian elasticities |
| Non-elasticity Standard Errors | $< 7.25 \times 10^{-6}$ | Delta-method covariance matrix |
| Log-Likelihood | $< 4.99 \times 10^{-8}$ | Relative difference |

The results show that `pyquaidsce` and Stata produce virtually identical point estimates and elasticities.

---

## 3. Automated Unit & Theory Tests

In addition to the Stata benchmark, `pyquaidsce` includes automated unit tests that verify:
- **Analytic Jacobians**: Checked against numeric finite-difference Jacobians.
- **Economic Theory Restrictions**: Confirming that Slutsky symmetry, price homogeneity, and adding-up restrictions hold at model-consistent shares.
- **Probit Estimator**: Validated against SciPy BFGS optimization and analytic information matrices.
- **Input Validation**: Ensuring appropriate error messages for invalid dimensions or non-positive prices/expenditures.
- **Control Function Extension**: Parameter maps, independent selection layouts,
  analytic/full/fast Jacobians, 220 conditional finite-difference derivatives,
  API gates, and end-to-end recovery checks.
- **Integrated `ivexp` Path**: Reduced-form OLS orthogonality and rank checks,
  excluded-instrument F diagnostics, exact equivalence to the same manually
  generated residual, validation gates, and verification that every bootstrap
  replication re-estimates the reduced form.
- **Release Integration**: Presence of 1.0.2/1.1/1.2 APIs, finite analytical
  covariance handling, real cooperative timeout behavior, and Stata bridge
  argument/matrix forwarding.

Release 1.5.0 collects **38 Python tests**; all pass from the repository, and
all 38 also pass from the built source distribution. The wheel passes an
isolated import/API smoke test. A separate
five-replication `pr`/`xb` bootstrap
smoke test and its JSON results are stored under
`benchmarks/release_120/results/`.

### Running the Test Suite:

```bash
python -m unittest discover -s tests -v
```

---

## 4. Reproducing the Benchmark

To reproduce the benchmark comparison locally:

```bash
cd benchmarks/cquaids_ifgnls_4g_20k/
python generate_data.py       # Generates benchmark_cquaids_4g_20k.dta
python run_python.py          # Runs estimation in Python and saves results
stata -b do run_stata.do      # (Optional) Runs estimation in Stata if installed
python compare_results.py     # Computes differences between Python and Stata
```
