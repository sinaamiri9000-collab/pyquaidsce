# Validation against Stata `quaidsce`

Version 1.0.1 uses one public, reproducible cross-implementation benchmark as
its main end-to-end validation reference:

[`benchmarks/cquaids_ifgnls_4g_20k/`](../benchmarks/cquaids_ifgnls_4g_20k/)

The benchmark deliberately exercises the full censored estimator rather than a
simplified uncensored case.

## Controlled specification

| Item | Setting |
|---|---|
| Observations | 20,000 |
| Goods | 4 |
| Demographic variables | 3 |
| Model | censored QUAIDS (CQUAIDS) |
| Estimator | IFGNLS |
| User-supplied initial vector | none |
| Bootstrap | disabled |
| Data | the same deterministic synthetic `.dta` file in both programs |
| Python compatibility settings | `first_stage_predict="pr"`, `strict_stata=True` |

All four goods contain genuine zero budget shares, so the first-stage
participation probits and censoring correction are active. The benchmark data
are generated from a fixed seed and can be regenerated from the included
script.

## Numerical agreement

The stored Stata and Python outputs contain 113 returned values. The
full-precision comparison gives:

| Quantity | Maximum/relative difference |
|---|---:|
| structural parameters (`alpha` through `rho`) | `1.68e-05` |
| first-stage Probit parameters (`tau`) | `1.21e-07` |
| expenditure/Marshallian/Hicksian elasticities | `7.34e-07` |
| non-elasticity standard errors | `7.25e-06` |
| log-likelihood, relative difference | `4.99e-08` |

The Stata and Python results are **almost identical**. The largest elasticity
difference is below `1e-6`; the two implementations are not claimed to be
bit-for-bit identical.

The repository stores the raw Stata log, Python output, and a
row-by-row comparison file so that this statement can be independently checked.

## Independent internal checks

The automated test suite also checks implementation components that do not
require Stata, including:

- parameter restrictions and the free-to-full parameter map;
- the analytic Jacobian against finite differences;
- the optimized free-space Jacobian against its reference construction;
- the first-stage Probit implementation against an independent optimizer;
- input validation and failure on invalid switches.

The end-to-end regression test uses the same public 20,000-observation benchmark
rather than maintaining a second Stata reference dataset.

Run the tests with:

```bash
python -m unittest discover -s tests -v
```

## Reproducing the cross-program comparison

From `benchmarks/cquaids_ifgnls_4g_20k/`:

```bash
python generate_data.py
python run_python.py
python compare_results.py
```

Run the Stata side with:

```stata
do run_stata.do
```

See the benchmark README for the exact timing boundary and interpretation of
runtime results.
