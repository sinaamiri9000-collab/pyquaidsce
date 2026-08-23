# Controlled CQUAIDS / IFGNLS benchmark

This is the public numerical-agreement and runtime benchmark for
`pyquaidsce` 1.0.1.

## Design

| Item | Setting |
|---|---|
| Observations | 20,000 |
| Goods | 4 |
| Demographic variables | 3 |
| Model | censored QUAIDS (CQUAIDS) |
| Estimator | IFGNLS |
| User-supplied initial vector | none |
| Bootstrap | disabled |
| Data | identical synthetic `.dta` file in Stata and Python |
| Python Stata-compatibility path | `first_stage_predict="xb"`, `strict_stata=True` |

The dataset is fully synthetic and generated from a fixed seed. Every good
contains genuine zero budget shares, so all four participation probits and the
Shonkwiler–Yen censoring correction are active. Budget shares sum to one in
every observation and no budget share equals one.

The synthetic data-generating process is intended for implementation
comparison and runtime benchmarking. It is not a Monte Carlo design intended
to recover a known vector of structural CQUAIDS parameters.

## Recorded results

The Stata and Python estimates are **almost identical**. The full-precision
comparison covers 113 returned values:

| Quantity | Maximum/relative difference |
|---|---:|
| structural parameters (`alpha` through `rho`) | `1.68e-05` |
| first-stage Probit parameters (`tau`) | `1.21e-07` |
| expenditure/Marshallian/Hicksian elasticities | `7.34e-07` |
| non-elasticity standard errors | `7.25e-06` |
| log-likelihood, relative difference | `4.99e-08` |

The maximum elasticity difference is below `1e-6`. The two implementations are
not claimed to be bit-for-bit identical.

## Same-machine runtime comparison

Both recorded runs were made on the same Windows machine and used the same
stored dataset and model specification.

| Implementation | Runtime |
|---|---:|
| Stata 19.5 `quaidsce_c` | `1161.171 s` |
| Python 3.14.5 `pyquaidsce` | `26.033 s` |
| **Stata / Python runtime ratio** | **`44.60x`** |

Stata's built-in `timer` was placed immediately around the complete
`quaidsce_c ... method(ifgnls)` point-estimation command. Python used
`time.perf_counter()` immediately around `quaidsce(...)`. Data loading and
result-file writing were outside both timing intervals.

The Stata log reports two processors available to Stata. The Python log from
this recorded run does not include an explicit BLAS/OpenMP thread-limit field.
Accordingly, **44.60x is reported as a same-machine wall-clock speedup**, not as
a per-core efficiency comparison.

## Files

```text
cquaids_ifgnls_4g_20k/
├── README.md
├── benchmark_manifest.json
├── generate_data.py
├── run_python.py
├── run_stata.do
├── compare_results.py
├── data/
│   └── benchmark_cquaids_4g_20k.dta
└── results/
    ├── stata_benchmark.log
    ├── python_benchmark.log
    ├── python_results.csv
    ├── python_runtime.json
    ├── stata_python_comparison.csv
    └── comparison_summary.json
```

The files under `results/` are the actual stored outputs used for the headline
comparison. The raw logs are retained so that both timing boundaries and
reported model outputs can be inspected directly.

## Reproduce the benchmark

From this directory, regenerate the deterministic dataset if desired:

```bash
python generate_data.py
```

Install `pyquaidsce==1.0.1` and run the Python side:

```bash
python run_python.py
```

The current `run_python.py` records additional environment metadata and sets a
two-thread limit for common BLAS/OpenMP backends by default. A different limit
can be requested explicitly:

```bash
PYQUAIDSCE_BENCH_THREADS=2 python run_python.py
```

Run the Stata side:

```stata
do run_stata.do
```

Then compare all returned values:

```bash
python compare_results.py --same-machine
```

Use `--same-machine` only when the Stata and Python runtime files being compared
were actually produced on the same computer.

The comparison CSV contains every returned name, the Stata and Python values,
standard errors, absolute differences, and block labels.
