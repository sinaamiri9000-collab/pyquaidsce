# Validation and Numerical Evidence

`pyquaidsce` v1.7.0 is validated against econometric identities, numerical derivatives, cross-interface contracts, and archived Stata comparisons. Exact reproduction of known mistakes in the original Stata ado is **not** a v1.7.0 objective.

## 1. Current v1.7.0 validation policy

The current implementation uses two corrected behaviors:

- the Shonkwiler-Yen first stage always uses the Probit linear index, $\Phi(X'\tau)$ and $\phi(X'\tau)$;
- elasticities always use the corrected theoretical formulas.

Therefore, differences from the original ado are expected in specifications affected by its nested-Probit transformation or documented elasticity indexing/macro mistakes. See [Stata Compatibility](stata-compatibility.md).

## 2. Automated Python tests

The repository contains **38 Python tests** covering:

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

## 3. Archived Stata benchmark

The directory [`benchmarks/cquaids_ifgnls_4g_20k/`](../benchmarks/cquaids_ifgnls_4g_20k/) contains the deterministic 20,000-observation, 4-good, 3-demographic IFGNLS benchmark used during earlier releases. Its stored `results/` files are preserved unchanged for auditability.

For that historical recorded comparison, the maximum reported differences were:

| Quantity | Archived maximum/relative difference |
|---|---:|
| Structural parameters | `1.68e-05` |
| First-stage Probit parameters | `1.21e-07` |
| Expenditure/Marshallian/Hicksian elasticities | `7.34e-07` |
| Non-elasticity standard errors | `7.25e-06` |
| Log-likelihood, relative difference | `4.99e-08` |

These numbers document the historical compatibility implementation; they should not be presented as a fresh v1.7.0 exact-replication result. The current benchmark script uses the v1.7.0 canonical first-stage definition.

## 4. Historical same-machine timing

The archived same-machine wall-clock timings were 1,161.171 seconds for Stata 19.5 and 26.033 seconds for Python, a 44.60x ratio. This is retained as historical performance evidence. A publication that labels the number specifically as a v1.7.0 benchmark should rerun both sides with the final release and report the new measurements.

## 5. Reproducing the benchmark

```bash
cd benchmarks/cquaids_ifgnls_4g_20k/
python generate_data.py
python run_python.py
# optional, when Stata is installed:
# stata -b do run_stata.do
python compare_results.py --same-machine
```

When comparing current v1.7.0 output with the original ado, interpret discrepancies using the documented methodological differences rather than forcing legacy error replication.
