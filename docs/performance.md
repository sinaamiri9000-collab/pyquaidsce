# Performance and benchmark protocol

Performance is one of the practical motivations for `pyquaidsce`, but a runtime
claim is meaningful only when the data, model specification, initialization,
timing boundary, and computing environment are documented.

## Public benchmark

Version 1.0.1 uses one public performance benchmark:

[`benchmarks/cquaids_ifgnls_4g_20k/`](../benchmarks/cquaids_ifgnls_4g_20k/)

The design exercises the full censored estimator while remaining practical to
reproduce:

- 20,000 observations;
- 4 goods;
- 3 demographic variables;
- genuine zero shares in all four goods;
- censored QUAIDS;
- IFGNLS;
- no user-supplied initial parameter vector;
- no bootstrap;
- the identical stored `.dta` file in Stata and Python.

The data are synthetic and generated from a fixed seed.

## Numerical agreement

The stored Stata and Python outputs contain 113 returned values. The
full-precision comparison gives:

| Quantity | Difference |
|---|---:|
| max absolute structural-parameter difference | `1.67687e-05` |
| max absolute first-stage Probit-parameter difference | `1.21482e-07` |
| max absolute elasticity difference | `7.33670e-07` |
| max absolute non-elasticity standard-error difference | `7.25255e-06` |
| relative log-likelihood difference | `4.98649e-08` |

The Stata and Python results are **almost identical**. In particular, the
largest elasticity difference is below `1e-6`.

## Timing boundary

The Stata script uses Stata's built-in `timer` around the complete
`quaidsce_c ... method(ifgnls)` point-estimation command. Data loading occurs
before `timer on`.

The Python script uses `time.perf_counter()` immediately around `quaidsce(...)`.
Data loading and result-file writing are outside the timer.

Thus both intervals cover the first-stage participation models and the full
IFGNLS point estimator, but not disk I/O or bootstrap inference.

## Recorded same-machine runtimes

The stored final runs were produced on the same Windows machine:

| Implementation | Runtime | Recorded environment |
|---|---:|---|
| Stata | `1161.171 s` | Stata 19.5, Windows, 2 processors reported available |
| Python | `26.032762 s` | Python 3.14.5, Windows Server 2022, Intel64 Family 6 Model 85 |

The resulting same-machine wall-clock ratio is:

```text
1161.171 / 26.03276198497042 = 44.6042
```

or **44.60x** when reported to two decimal places.

The Python log from this particular stored run does not record an explicit
BLAS/OpenMP thread limit. For that reason, the 44.60x figure is a
**same-machine wall-clock speedup**, not a claim about per-core computational
efficiency. The current reproduction script records a thread limit explicitly
to make future reruns even easier to document.

## Bootstrap timing

This benchmark deliberately excludes bootstrap inference. A bootstrap
replication re-estimates the first-stage probits and nonlinear demand system, so
its total runtime depends heavily on the requested number of replications and
parallel-worker configuration. Point-estimation performance is therefore
reported separately.
