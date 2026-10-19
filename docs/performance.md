# Performance & Benchmarking

One of the primary motivations for creating `pyquaidsce` is estimation speed. In applied empirical demand analysis, running IFGNLS models or bootstrap replications in Stata can become a major time bottleneck.

This document describes the computational design of `pyquaidsce` and the benchmark results.

---

## 1. Why is pyquaidsce Fast?

Several computational strategies make `pyquaidsce` significantly faster than conventional Stata implementations:

1. **Direct Free-Space Analytic Jacobians**:
   Instead of evaluating the full unconstrained Jacobian ($N \times m \times K_{full}$) and then multiplying by large restriction matrices, `pyquaidsce` evaluates the derivative directly in the reduced parameter space ($K_{free}$). This avoids substantial matrix multiplications in every Gauss-Newton iteration.
2. **Symmetric BLAS Operations (`dsyrk`)**:
   Accumulation of normal equations ($J' \Sigma^{-1} J$) is accelerated using low-level symmetric rank-k updates from BLAS.
3. **Memory Chunking**:
   The Jacobian is processed in chunks of observations rather than building the entire matrix in RAM simultaneously, keeping memory usage minimal.
4. **Controlled BLAS Parallelism**:
   Point estimation and each bootstrap worker use one BLAS thread by default (`blas_threads=1`). The limit is temporary in Python and the previous BLAS setting is restored after estimation.

---

## 2. pyquaidsce vs quaidsce

Both programs fit the same four-good censored QUAIDS model using IFGNLS on
2,665 households, with 2 demographic variables and 200 bootstrap replications
using 8 workers on the same computer.

| Program | Time (seconds) |
|---|---:|
| pyquaidsce | 122.05 |
| quaidsce_c | 10,810.18 |

pyquaidsce was about **88.6 times faster** in this run.

The code, data, and results are in [benchmark/pyquaidsce vs quaidsce/](../benchmark/pyquaidsce%20vs%20quaidsce/).

---

## 3. Bootstrap Performance Tips

When running large bootstrap routines (`reps=500+`):

- **Use Multiprocessing**: Specify `n_jobs` equal to the number of available physical CPU cores:
  ```python
  res = quaidsce(..., reps=500, n_jobs=8)
  ```
  If you run parallel bootstrap from a Python script, put the estimation inside `main()` and add this at the end:
  ```python
  if __name__ == "__main__":
      main()
  ```
- **Warm-start Option**: Setting `bootstrap_start="warm"` reuses the full-sample parameter estimates as starting values for each bootstrap draw, further accelerating bootstrap convergence.
- **Time limit per replication**: `rep_timeout=900` allows each bootstrap
  replication to run for up to 900 seconds. If it takes longer, it is counted
  as unsuccessful and the program moves on to the next one.
- **Start method**: leave `mp_context=None` for the safe `spawn` default. An
  explicit context should be used only after testing it on the target operating
  system and BLAS build.
