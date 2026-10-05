# Findings from the unchanged pyquaidsce 1.6.0

The existing 42-test suite passes, and the independent demand-equation,
Jacobian, first-stage, reduced-form and bootstrap-accounting checks provide
substantial positive evidence. The expanded audit nevertheless does **not**
meet its release validation gate. Exact counts are in `results/summary.json`;
individual observations and fixed thresholds are in `results/checks.csv`.

This report separates mathematical/inference concerns, documented API
contracts, numerical limitations on small4, and supplemental interface source
checks. Nothing in `src/`, `R/`, `stata/`, the old tests, or small4 was patched
to make the experiment pass.

## 1. Conditional analytical covariance of NLS — high priority

Case: `methodology-NLS_covariance`.

The NLS point estimator minimizes identity-weighted residual sum of squares.
Its conditional multivariate NLS covariance, under homogeneous cross-equation
residual covariance Sigma, is

```
A = sum_t J_t' J_t
B = sum_t J_t' Sigma J_t
V_NLS = A^-1 B A^-1
```

The current implementation instead reports the inverse GLS normal matrix
`inv(sum_t J_t' Sigma^-1 J_t)` even for `method='nls'`. That is the efficient
GLS expression, not generally the covariance of the identity-weighted NLS
estimator. On the censored small4 baseline, the maximum absolute covariance
difference is about **0.1721**. Reported free-parameter SEs are approximately
**0.303–0.866 times** the conditional NLS sandwich SEs.

This concerns the conditional second-stage analytical approximation. Neither
expression by itself accounts for estimated first-stage regressors; full
generated-regressor inference still needs the full bootstrap. The audit does
not conflate that separate issue with this NLS weighting mismatch.

Evidence: `results/cases/methodology-NLS_covariance.json`;
implementation: `src/pyquaidsce/nlsur.py`, `_finish`.

## 2. Default stopping rule ignores advertised tolerances — high priority

Case: `numeric-standard_tolerance_contract`.

Under `stop_rule='standard'`, `gauss_newton` uses literal thresholds
`mreldif < 1e-5`, `rel < 1e-7`, and `nrtol < 1e-5`. The advertised `tol` and
`nrtol_stop` values are not consulted in this branch; they do participate in
the separate tight branch. This is a source-backed behavioral contract issue,
not inferred solely from two fits happening to produce similar estimates.

Evidence: `results/cases/numeric-standard_tolerance_contract.json`;
implementation: `src/pyquaidsce/nlsur.py`, `gauss_newton`.

## 3. IFGNLS/LM and restricted selection designs — numerical limitations

Five of the 72 core default-control fits return `converged=False`: censored
QUAIDS/LM with both starts, censored AIDS/LM with zero start, and uncensored
QUAIDS with demographics/LM with both starts. Their GN counterparts converge.
Two tighter baseline IFGNLS/LM follow-ups also fail the convergence gate.
These are not mislabeled as successful estimates. In contrast, the tighter
NLS and FGNLS baseline runs do converge.

Some valid custom selection designs also do not converge. The independent
covariates plus omitted-expenditure design exhausts the **audit's cooperative
90-second point-fit bound**; the alternative-start follow-up uses 20 seconds.
The package's generic deadline error mentions bootstrap even though these
cases are point fits. This does not prove that an unlimited run can never
finish.

The intercept-only Probit design raises `LinAlgError` from a non-positive
definite residual covariance with the default zero start. An alternative
linear start can return an apparently converged fit whose inference checks
fail, including singular/unstable covariance diagnostics. The original four
budget shares sum to one; degeneracy of residual/parameter covariance must
be assessed before treating this specification as usable.

Small absolute covariance asymmetry (around 2.3e-10 in some already
nonconverged omitted-expenditure fits) is recorded by the strict algebra
check and should not be confused with the much larger degeneracy in other
selection fits. No threshold was loosened to erase either observation.

Evidence: core/selection/followup case files, `results/core-matrix.svg`,
`results/optimizer-comparisons.csv`.

## 4. Natural bootstrap failures — optimizer reliability, not accounting

The four requested 100-draw batches have success counts:

| Model | Seed | Successful / requested |
|---|---:|---:|
| Base | 16001 | 99 / 100 |
| Base | 16002 | 98 / 100 |
| Internal IV | 16001 | 99 / 100 |
| Internal IV | 16002 | 98 / 100 |

Thus **394/400** large-batch draws converge. All six failures are second-stage
nonconvergence. The eight-draw uncensored smoke batch has one additional
natural nonconverged draw. The package correctly excludes failed draws and
computes covariance/SE/CI from successful draws; those accounting assertions
pass. The audit's **full-success stress gate is stricter than the documented
ability to discard failed draws**. A flagged batch is not a claim that its
covariance arithmetic is broken.

Serial/parallel ordering, repeatability at fixed seeds, explicit spawn,
available forkserver/fork, reduced-form rebuilding, timeout, noncooperative
watchdog termination, and deterministic injected-failure handling are tested.
The realized `seed=None` draws are archived but not claimed reproducible from
an unspecified seed.

No replacement seeds or retries were used to conceal the failed draws.
`results/bootstrap-failures.csv` identifies them. For example:

```bash
python -m validation.reproduce --failed-draw boot-large-base-16001:30
```

`results/bootstrap-stability.csv` reports SE drift across seeds. With 100
draws, it is a precision diagnostic, not a coverage theorem or an automatic
certificate that every SE is stable.

## 5. Invalid numerical controls — validation gaps

The parameterized negative cases expose acceptance of nonfinite `tol`,
`nrtol_stop`, and `sigma_tol`; nonpositive/nonfinite `inner_nrtol_early`;
negative/fractional `reps`; and zero/negative/fractional `n_jobs`.
`boot_sigma_tol` also lacks consistent upfront validation: some invalid
values fail only later as an all-replications-failed RuntimeError, while
infinity is accepted. Some fractional iteration/chunk values raise TypeError
only when a downstream integer operation is reached.

These findings are separate from the successful rejections for missing or
duplicated variables, nonpositive prices/expenditure, prohibited CF/IV
combinations, rank-deficient reduced forms, and other negative scenarios.

Evidence: `results/scenarios.csv`, group `invalid`, and mapped cases in
`results/argument-coverage.csv`.

The purported complete main-estimator API reference omits the public `log`
callback. It is exercised by the audit and marked as undocumented in the
argument coverage table; this is a documentation gap, not a failed callback.

## 6. Public helper/output contracts — concrete defects

- `public-probit_no_constant`: `probit(..., add_constant=False)` returns a
  coefficient vector without a constant, but `ProbitResult.xb(X)` still
  removes the last coefficient and treats it as the intercept. The test
  raises a matrix-dimension ValueError.
- `public-reduced_form_labels`: wrong lengths of public label arrays are
  accepted, leaving fewer names than coefficients; `named()` would silently
  truncate via `zip`. Label/shape validation is needed.
- `output-summary-no_elasticities-bootstrap`: `summary(elasticities=False)`
  suppresses elasticity coefficient rows without bootstrap but ignores the
  flag when a bootstrap is attached.
- `output-invalid_level-*`: percentages below 0, above 100, or nonfinite are
  not consistently rejected by summary/normal/percentile CI helpers. Some
  produce reversed or undefined intervals. Mathematically well-defined
  0/100 percentile endpoints are tested separately, not declared invalid.

Each case JSON contains the expected and observed behavior and, where
applicable, the exception traceback.

## 7. R wrapper source contracts — no native execution claim

- `R-static-fitted_backend`: the wrapper calls `res_py$fitted_shares()`, but
  `QuaidsceResults` has no such method. The R `tryCatch` fallback makes
  `fitted.values` and residuals NULL, despite the advertised fitted/residual
  methods. This is a source/backend contract finding; no native R result is
  invented.
- `R-static-analytic_SE`: R summary reconstructs SE as `sqrt(diag(V))`.
  For this censored non-bootstrap baseline that yields **36 zero elasticity
  SEs**, while Python intentionally reports 36 undefined/NaN entries. The
  wrapper needs to preserve unsupported uncertainty instead of reconstructing
  it from the placeholder covariance block.

The real Python estimator behind the Stata bridge and its real asynchronous
bootstrap subprocess are exercised through simulated `sfi` transport and
agree with direct Python results. This is not a test of licensed Stata's
ado parser/display/runtime. Native R and Stata are explicitly `NOT_RUN`.

## Reproduction and next work

Run the complete experiment with `python -m validation.run --large-reps 100`.
Use `--only CASE_ID --output /tmp/repro` for a focused case. These flags retain
main-matrix reference fits automatically. For helper/output/interface IDs, use
`--contract-only CONTRACT_ID --output /tmp/helper-repro`. Both commands retain
reference fits automatically; do not use an old mixed output directory for
a subset. The runner exits nonzero for failed assertions or unexplained
ordinary-fit nonconvergence while continuing to collect all evidence.

The package source is intentionally unchanged in this branch. The priority
before publication is to resolve the inference/stopping contracts, assess
singular/weakly identified selection specifications, fix the helper/validation
and wrapper defects, and then rerun the same frozen design. A numerical
limitation may instead need a clearly documented unsupported specification;
it should not be concealed by changing the audit's seed or thresholds.
