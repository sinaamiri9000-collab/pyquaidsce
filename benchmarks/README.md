# Benchmarks

This directory contains **one public performance benchmark** for `pyquaidsce`
1.0.1. Older exploratory timing files are intentionally not kept here.

- [`cquaids_ifgnls_4g_20k/`](cquaids_ifgnls_4g_20k/) — controlled censored
  QUAIDS / IFGNLS benchmark with 20,000 observations, 4 goods, and 3
  demographic variables.

The benchmark has two separate purposes:

1. compare the numerical results returned by Stata `quaidsce` and
   `pyquaidsce`; and
2. measure the point-estimation runtime under an explicitly documented setup.


The stored same-machine run records a **44.60x wall-clock speedup** for Python, while the Stata and Python estimates are **almost identical**.
