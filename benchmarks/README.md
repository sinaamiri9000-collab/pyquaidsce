# Benchmarks

This directory contains reproducible benchmark assets for `pyquaidsce`. The
stored numerical-comparison outputs were recorded with an earlier
compatibility-capable release and are retained as **historical validation and
timing evidence**. They are not a fresh v1.7.0 exact-replication claim.

- [`cquaids_ifgnls_4g_20k/`](cquaids_ifgnls_4g_20k/): controlled synthetic
  household benchmark (20,000 observations, 4 goods, 3 demographics, censored
  QUAIDS with IFGNLS). The current Python runner uses the canonical v1.7.0+
  linear-index Shonkwiler-Yen correction.

### Archived recorded results

For the historical recorded comparison, maximum parameter differences were
below `1.68e-05`, maximum elasticity differences were below `7.34e-07`, and the
same-machine wall-clock ratio was approximately **44.6x** (26.0 s in Python
versus 1,161.2 s in Stata). Because v1.7.0 intentionally no longer reproduces
known errors in the original ado, these numerical-agreement figures should be
reported as historical rather than as a current exact-equivalence target.

See [`docs/validation.md`](../docs/validation.md),
[`docs/stata-compatibility.md`](../docs/stata-compatibility.md), and
[`docs/performance.md`](../docs/performance.md) for the current validation policy
and benchmark interpretation.
