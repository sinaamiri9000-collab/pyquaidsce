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

We estimated the same six-good censored QUAIDS model using both programs. The data come from the Iranian Household Income and Expenditure Survey (2023–2024).

Both models used IFGNLS, 37,656 observations, eight demographic variables, and the same starting values. Each program was run once without bootstrap.

| Program | Time (seconds) | Log-likelihood |
|---|---:|---:|
| pyquaidsce | 119.685 | 379314.2335272935 |
| quaidsce_c | 21,052.243 | 379314.2334614899 |

In this run, pyquaidsce was about **176 times faster**. The original Stata command did not report a convergence flag, so we do not claim that its convergence was independently confirmed.

The files for this comparison are described in [benchmark/pyquaidsce vs quaidsce/](../benchmark/pyquaidsce%20vs%20quaidsce/). The original timing log was recorded on 10 October 2026. A later run can be made with the updated Stata connection in the publication branch.

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
- **Per-replication deadline**: `rep_timeout=900` drops a replication that passes
  the 15-minute deadline. Cooperative checks normally stop the optimizer first;
  a parent-side watchdog terminates the disposable child process if native code
  remains stuck. Enabling this hard backstop adds process-start overhead.
- **Start method**: leave `mp_context=None` for the safe `spawn` default. An
  explicit context should be used only after testing it on the target operating
  system and BLAS build.
