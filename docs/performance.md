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
4. **Efficient Inexact-Outer IFGNLS**:
   Early outer iterations of the covariance matrix $\Sigma$ do not require solving the inner Gauss-Newton steps to machine precision. Tight tolerances are only applied as the outer fixed-point stabilizes, cutting the total number of optimization steps.
5. **Parallel Bootstrap Worker Management**:
   Each worker process pins its internal BLAS threads to 1, preventing CPU thread oversubscription and ensuring maximum throughput across cores during bootstrap replications.

---

## 2. Controlled Benchmark Results

The benchmark setup consists of **20,000 observations, 4 goods, 3 demographics, and IFGNLS estimation** evaluated on the exact same hardware:

| Platform / Software | Point Estimation Runtime | Speedup |
|---|---:|---:|
| **Stata 19.5** | 1,161.2 seconds (~19 min 21 s) | 1.0x (Baseline) |
| **pyquaidsce (Python 3.14)** | **26.0 seconds** | **~44.6x Faster** |

*Note: Timings measure point-estimation computation time (excluding disk I/O).*

---

## 3. Bootstrap Performance Tips

When running large bootstrap routines (`reps=500+`):

- **Use Multiprocessing**: Specify `n_jobs` equal to the number of available physical CPU cores:
  ```python
  res = quaidsce(..., reps=500, n_jobs=8)
  ```
- **Warm-start Option**: Setting `bootstrap_start="warm"` reuses the full-sample parameter estimates as starting values for each bootstrap draw, further accelerating bootstrap convergence.
