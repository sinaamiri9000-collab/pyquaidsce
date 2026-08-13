# Benchmark assets

`small4.dta` and `small4.log` form the public, compact Stata/Python regression
benchmark used by the automated tests and validation utilities.

- `benchmark_small.do` creates/runs the Stata benchmark from the separately
  obtained upstream `quaidsce` source/data tree.
- `small4.log` contains the high-precision Stata reference output.
- `small4.dta` is the corresponding 2,000-observation four-good dataset.
- `timing_11goods_37k.log` records a performance run only; it is not a universal
  speed guarantee.

The full upstream Stata repository is intentionally not vendored here. Set
`QUAIDSCE_STATA_REPO` when running utilities that need it.
