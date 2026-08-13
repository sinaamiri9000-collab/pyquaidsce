# Benchmarks

This directory contains the reproducible benchmark comparing `pyquaidsce` against Stata `quaidsce`.

- [`cquaids_ifgnls_4g_20k/`](cquaids_ifgnls_4g_20k/): Controlled benchmark using synthetic household data (20,000 observations, 4 goods, 3 demographics, censored QUAIDS with IFGNLS).

### Key Takeaways:
1. **Numerical Equivalence**: Parameter differences between Stata and Python are $< 1.68 \times 10^{-5}$, and elasticity differences are $< 7.34 \times 10^{-7}$.
2. **Speed Comparison**: Python achieves a **~44.6x wall-clock speedup** (26.0s vs. 1,161.2s in Stata on the same machine).

See [`docs/validation.md`](../docs/validation.md) and [`docs/performance.md`](../docs/performance.md) for full details.
