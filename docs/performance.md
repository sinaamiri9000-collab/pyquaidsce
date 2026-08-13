# Performance notes

Performance is one of the practical motivations for `pyquaidsce`, but timing
claims must be tied to a model, convergence path, hardware, BLAS backend, and
sample size.

## Recorded version-1.0.1 timings

The validation workspace recorded the following `ifgnls` runs on two cores with
`scipy-openblas`:

| Problem | Wall time | Iteration information |
|---|---:|---|
| 4 goods, 2,000 observations, 2 demographics | ~1.8 s | small benchmark |
| 17 goods, 15,147 observations, 3 demographics | ~30 min | 42 outer, 417 GN steps |
| 11 goods, 37,000 observations, 3 demographics | ~10.0 min | 35 outer, 200 GN steps |

The raw log for the 11-good benchmark is in
[`../bench/timing_11goods_37k.log`](../bench/timing_11goods_37k.log).

These are **not** universal speed claims against Stata. A fair comparison should
report the same data, same sample, same model, same start, convergence settings,
hardware, Stata version, Python environment, and whether bootstrap inference is
included.

## Why model size matters

The dominant Gauss–Newton work grows with the number of observations, equations,
and especially the number of free parameters. The package accumulates normal
equations in chunks and uses BLAS operations rather than materializing the full
stacked Jacobian, which is the main computational design behind the lower memory
footprint and observed speed.

## Conditioning matters

A model can spend most of its time taking very small or repeatedly halved steps
if variables are poorly scaled. In the 17-good validation exercise, converting
a conditional subsystem to the economically intended normalized shares changed
the numerical conditioning dramatically.

This does **not** mean users should normalize shares mechanically. It means the
estimated system and expenditure definition must be internally coherent. Check
the economic construction first, then numerical conditioning.

## Bootstrap timing

A bootstrap replication re-estimates the first-stage probits and the nonlinear
demand system. A large `reps` value can therefore dominate total runtime even
when one point estimate is fast. Establish the time and convergence rate using a
small diagnostic bootstrap before launching the final inference run.

## Controlled Stata–Python comparison

A same-data, same-machine benchmarking protocol is provided in
[`../bench/runtime_comparison/`](../bench/runtime_comparison/). It intentionally
separates historical timing evidence from the final controlled speed ratio. A
public speedup claim should only be added after both programs have been timed on
the exact same benchmark file and hardware.
