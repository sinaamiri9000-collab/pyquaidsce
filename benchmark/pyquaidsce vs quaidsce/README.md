# pyquaidsce vs quaidsce

This comparison uses the Iranian Household Income and Expenditure Survey (2023–2024): 37,656 households, 6 food groups, and 8 demographic variables. Both commands estimate censored QUAIDS with IFGNLS and the same starting values. Each model was run once, without bootstrap.

| Measure | pyquaidsce | quaidsce_c |
|---|---:|---:|
| Time (seconds) | 119.685 | 21,052.243 |
| Log-likelihood | 379314.2335272935 | 379314.2334614899 |

pyquaidsce took about **176 times less time** in the recorded run. The Stata command completed but did not report a convergence flag.

The files for this comparison are:
- `benchmark.do`: Stata commands for running both models.
- `benchmark_2026-10-10.log`: the original Stata results (to be added).
- `six_commodity_input.dta`: the prepared data (to be added after permission to share it is confirmed).

Run from the root of the repository:

```stata
do "benchmark/pyquaidsce vs quaidsce/benchmark.do"
```

The original log was produced with a temporary fix for the Python executable inside Stata. The new `benchmark.do` does **not** use that fix: the `publication` version of pyquaidsce reads the correct Python setting automatically. The times shown above have not yet been reproduced with this new script.
