# Convergence architecture comparison for pyquaidsce 1.6.0

This research experiment compares three architectures on `bench/small4.dta`:

1. IFGNLS with native SciPy `least_squares` (TRF / MINPACK LM) inside native
   `fixed_point(method="iteration")` on `[theta, vech(Sigma)]`.
2. Direct concentrated Gaussian likelihood with native SciPy BFGS / L-BFGS-B.
3. An independent Python reconstruction of the uploaded `nlsur.ado` iteration
   and stopping logic, using its numerical derivatives / the package's analytic
   derivatives as separate variants.

**No native Stata executable was run.** The reconstruction uses LAPACK pivoted
QR with the documented Mata rank tolerance, Cholesky whitening, and vectorized
accumulation. These preserve the algebra and stopping logic for positive-definite
Sigma, but do not establish bitwise equivalence to Mata/Stata. The uploaded ado
header is `version 1.6.2 30jan2023`; the package under study is **pyquaidsce 1.6.0**.

The primary experiment contains 120 fits. Linear initialization is excluded at
the user's request. Zero is the sole cold start. Recorded package fits from
`core-*-ifgnls-gn-zero.json` provide common warm starts, not ground truth.
First-stage probit estimates are identical and held fixed within each model.
The six model settings are censored QUAIDS/AIDs with demographics (`cq`,`ca`),
uncensored QUAIDS/AIDS with demographics (`uq`,`ua`), and uncensored
QUAIDS/AIDS without demographics (`nq`,`na`). The sample has 2,000 observations.

There are 72 primary configurations: six model settings × two starts × six
backend variants. A further 48 configurations tighten the native controls for
`cq`,`ca`,`uq`,`ua`, with both starts. The boolean `tightened` describes a
sensitivity experiment, not a new package `stop_rule`.

| Backend | Baseline controls | Tighter sensitivity controls |
|---|---|---|
| SciPy TRF / LM | `ftol=xtol=gtol=1e-8`, `x_scale='jac'`, `max_nfev=400`; outer `xtol=1e-8`, `maxiter=200` | Inner and outer tolerances `1e-10` |
| nlsur reconstruction | `eps=1e-5`, `tau=1e-3`, `ifgnlseps=1e-10`, `delta=4e-7`; inner limit 400, outer limit 200 | `eps=1e-8`, `ifgnlseps=1e-12` |
| BFGS | Per-observation loss; `gtol=1e-5`, `xrtol=0`, `maxiter=1000` | `gtol=1e-8` |
| L-BFGS-B, without parameter bounds | Per-observation loss; `gtol=1e-5`, `ftol=2.220446049250313e-9`, `maxiter=1000`, `maxls=20`, `maxfun=15000`, `maxcor=10` | `gtol=1e-8`, `ftol=1e-12` |

The numerical-derivative increment follows the actual ado:
`4e-7 * (abs(theta_j) + 4e-7)`. Differences are taken between predicted shares,
as in `FunctionDeriv`, rather than between residuals. The analytic variant
changes only the derivative provider. Both nlsur variants use the documented
default pivoted-QR tolerance `1e-13 * trace(abs(R)) / rows(R)`.

The direct-likelihood solvers use a fixed coordinate transform: inverse RMS
Jacobian column norms at the common zero / identity-weight anchor. Exactly zero
columns receive scale 1. This is Jacobian column scaling, with an explicit fixed
anchor; the objective and native termination predicates are unchanged. The
scaling is not dynamically equal to `least_squares(x_scale='jac')`.
Invalid model trials return infinite loss and a nonfinite gradient. Improving
this domain handling is a remaining engineering issue for the profile route.

Each fit has a 45-second operational bound. Observation weights are unity and
BLAS threads are limited to one. Budgets and stopping formulas have different
meanings across engines; equal decimal tolerances would not imply equal accuracy.

## Reproduce

From the repository root, with the package's dependencies plus `threadpoolctl`:

```bash
PYTHONPATH=src python -m validation.convergence.compare
PYTHONPATH=src python -m validation.convergence.cross_start
PYTHONPATH=src python -m validation.convergence.verify
PYTHONPATH=src python -m validation.convergence.summarize
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 PYTHONPATH=src python -m unittest discover -s tests -t . -v
```

The second command runs 24 separately labelled, post hoc shared-point probes.
The third command checks all primary/cross-start outcomes, independent model
equations, the profile gradient, whitening and the concentrated likelihood.

## Findings

At the baseline controls and **zero** start:

| Model | SciPy TRF LL | SciPy LM LL | Profile BFGS LL | nlsur analytic LL |
|---|---:|---:|---:|---:|
| cq | 6913.445850 | 6913.445840 | 10399.155423 **unsuccessful status** | 10234.430476 |
| ca | 6720.399612 | 6720.399612 | 9964.846201 | 9964.846200 |
| uq | 3142.066848 | 3137.557142 | 3137.557144 | 3137.557144 |
| ua | 3136.026726 | 3136.026726 | 3136.026726 | 3136.026726 |
| nq | 3058.095531 | 3058.095531 | 3058.095531 | 3058.095531 |
| na | 3046.122060 | 3046.122060 | 3046.122060 | 3046.122060 |

LL comparisons apply **within the same model and frozen first stage** only.
The larger cq profile likelihood was retained by SciPy IFGNLS and both nlsur
variants when started there. Maximum absolute LL change across the 24 shared-point
probes was `1.27e-9`. BFGS again reported precision loss at the higher cq point;
that status was preserved. The experiment supports multiple retained stationary
points, rather than an explanation based solely on different stopping rules.
It does not prove that any point is a global maximum or classify every point's
curvature. Tightening the controls did not resolve the large cq/ca LL differences.

Native success counts across the 120 primary configurations are TRF 20/20,
LM 20/20, nlsur numerical 20/20, nlsur analytic 20/20, BFGS 13/20, and
L-BFGS-B 14/20. These counts include warm starts and sensitivity runs and are
not accuracy rankings. L-BFGS-B's successful baseline fits can retain measurable
nonstationarity; unsuccessful BFGS fits can be extremely close to stationary.
The native nlsur return flag also does not certify every earlier inner solve:
one tighter warm-start uq analytic run had an earlier inner limit before its
final successful solve.

There are **760/760 passed experimental verification checks** and **42/42 passed
existing unittest tests**. The maximum scaled directional gradient discrepancy
at the declared `h=1e-7` probe was `1.19e-7`, below `2e-6`. Source and dataset
hashes agree before and after the experiment. The existing tests include their
original benchmark checks; they are not the selection criterion for this study.

For the next implementation experiment, the source-defined nlsur IFGNLS loop
with analytic derivatives is a reasonable base: it has explicit reference
controls, preserves the statistical estimator, and performed well from zero
here. SciPy remains a valid standard alternative and a useful independent check.
Direct profile likelihood is especially valuable as an independent check, but
this implementation needs stronger scaling, domain and precision-loss handling
before consideration as the default. The observed larger cq likelihood makes
comparison of independent starting paths important even for the nlsur route.
These are recommendations for further implementation, not production changes.

## Artifacts and scope

- `REPORT.fa.md`: detailed Persian explanation, exact formulas, corrections to
  the attached AI interpretation and architecture proposals.
- `provenance.json`: uploaded-file hashes, roles, header and ado line references.
- `references.json`: primary documentation and implementation references.
- `results/runs.csv` and `runs.json`: all 120 primary fits, including failures,
  native controls, statuses, parameter vectors, covariance and inner histories.
- `results/cross-start.csv` and `.json`: all 24 shared-point probes.
- `results/controls.csv`: native options and final stopping details per primary fit.
- `results/summary.csv`: native status counts and diagnostic summaries by start,
  backend and tolerance setting; these are not accuracy rankings.
- `results/gradient-probes.csv` and `.json`: all 108 directional derivative probes
  at three step sizes, including discrepancies at larger steps.
- `results/verification.csv` and `.json`: all 760 artifact/numerical checks.
- `results/environment.json`: versions, dataset/package hashes and scope.
- `results/scipy-native-termination.txt`: installed SciPy termination source.
- `results/existing-tests.log`: the original 42-test suite's outcomes.

This compares conditional second-stage point estimators. It does not implement
new coefficient covariance estimators, redo the original full API audit, refit
the first stage in bootstrap samples, or validate full censored-model maximum
likelihood. Those are separate tasks. The package estimator, dataset and earlier
validation results are unchanged. Uploaded Stata source remains in the local
attachments; its hashes and independently written reconstruction are recorded.
