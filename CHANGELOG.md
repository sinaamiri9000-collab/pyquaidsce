# Changelog

All notable user-visible changes to `pyquaidsce` are documented here.

## 1.2.0 — 2026-08-15

Consolidated release based on the feature-complete 1.0.2 estimator, with the
Stata integration from 1.1.0 and corrected bootstrap/standard-error work.

- Restored the 1.0.2 external demand control function, independent selection
  design, typed `FirstStageLayout`, control-function Jacobians, and conditional
  price/expenditure derivatives that were absent from the 1.1.0 branch.
- Retained the 1.1.0 `stata_bridge.py`, official Stata command/help/package
  files, documentation tree, benchmark suite, CI, and release workflow.
- Changed parallel bootstrap to a BLAS-safe, cross-platform `spawn` default,
  with configurable `mp_context`, unordered real-time completion reporting,
  and deterministic replicate ordering in the stored `b_star` matrix.
- Added a real cooperative `rep_timeout` checked inside first-stage Probits,
  nonlinear optimizer iterations, and chunked normal-equation accumulation,
  plus a parent-side watchdog that terminates a stuck replication process.
- Kept `res.V` finite when no bootstrap is requested. Unsupported analytical
  elasticity S.E.s are exposed as `NaN` through both `res.se` and
  `res.analytic_se`, avoiding false zero-uncertainty estimates without
  contaminating unrelated covariance contrasts.
- When bootstrap succeeds, synchronized `res.V`/`res.se` to the bootstrap
  covariance while retaining the conditional reference as
  `res.V_analytic`/`res.analytic_se`.
- Added release-integration tests and an explicit source-distribution manifest
  containing the datasets/logs required by the shipped regression tests.
- Made every shipped test and validation entry point import the current
  checkout's `src/` tree before any installed copy of the package.
- Exposed the external control-function/custom-selection design in the Stata
  bridge and ADO, with explicit gates for `first_stage_predict(xb)` and
  reduced-form-aware outer bootstrap inference.

## 1.1.0 — 2026-08-14

Stata integration release and comprehensive documentation expansion.

- **Stata Integration**: Added official `pyquaidsce.ado` and `pyquaidsce.sthlp` commands allowing Stata users to estimate censored QUAIDS directly inside Stata via the Python computational engine (`stata_bridge.py`), supporting full `e()` matrices, scalars, and postestimation commands (`test`, `lincom`, `outreg2`, `esttab`).
- **Stata Package Management**: Added `stata.toc` and `pyquaidsce.pkg` enabling direct 1-line installation in Stata via `net install`.
- **User Guide & API Reference**: Created comprehensive [User Guide](docs/user-guide.md) documenting every input parameter, optimizer setting, and attribute of `QuaidsceResults` with practical extraction recipes.
- **Documentation Refinement**: Simplified all documentation files to clear applied economics terminology.
- **Motivation & Citation**: Documented the empirical motivation for the package and updated citation metadata to 2026.

## 1.0.2 — 2026-08-12

- Added an external demand control function with latent placement
  `Phi * (wQ + cfcoef * residual) + delta * phi`.
- Added independent selection controls for price subsets/order, expenditure,
  covariates, and a selection control-function residual.
- Added typed first-stage layout metadata, the `cfcoef` parameter block,
  analytic full/fast Jacobians, and conditional fitted-share derivatives.
- Added strict gates for incompatible censoring, prediction, bootstrap, and
  collinearity combinations, plus 14 focused control-function tests.
- Preserved the complete 1.0.1 benchmark exactly when the extension is off.

## 1.0.1 — 2026-08-02

Validation and reliability release based on the original Python port.

- aligned the documented Stata-compatible IFGNLS path with the Gauss–Newton
  algorithm used by the benchmark;
- corrected propagation of estimator settings into bootstrap replications;
- stopped non-converged bootstrap replications from entering bootstrap standard
  errors as successful fits;
- corrected IFGNLS convergence reporting at iteration limits;
- made zero-start bootstrap the compatibility default while retaining an
  explicit warm-start option;
- repaired validation-tool paths and Windows-safe parallel bootstrap behavior;
- added collectable regression/theory tests and stronger input validation;
- documented compatibility switches, numerical conditioning, and full-precision
  Stata benchmark evidence.
- consolidated public Stata/Python validation and performance evidence into one
  reproducible controlled CQUAIDS/IFGNLS benchmark.

## 1.0.0

- Initial Python implementation of the censored QUAIDS workflow.
