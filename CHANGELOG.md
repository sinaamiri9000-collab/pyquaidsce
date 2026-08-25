# Changelog

All notable user-visible changes to `pyquaidsce` are documented here.

## 1.5.0 — 2026-08-25

- Added `ivexp` consistently to the Python, R, and Stata interfaces. It fits
  the internal OLS reduced form for log expenditure on log prices, Ray
  demographics, excluded instruments, and a constant.
- The generated residual automatically enters both the participation Probits
  and latent demand equations, using separate equation-specific coefficients.
- Added a typed reduced-form result with coefficients, covariance, fitted
  values, residuals, R-squared, adjusted R-squared, and the classical joint F
  test for excluded instruments. R exposes this as `fit$reduced_form`; Stata
  stores corresponding `e(reduced_form_*)` matrices/scalars.
- Enabled generated-regressor-aware bootstrap inference for `ivexp`: the
  reduced form is re-estimated inside every resample. Precomputed external
  residuals remain bootstrap-disabled because the package cannot rebuild their
  unknown generating equation.
- Preserved the external `control_function` and
  `selection_control_function` APIs and all existing elasticity formulas.
- Added focused numerical, equivalence, validation, bootstrap, R-interface,
  and Stata-bridge regression tests.

## 1.4.0 — 2026-08-24

- **Official R Package (`rquaidsce`)**: Introduced the complete R interface package with native S3 methods (`summary`, `coef`, `vcov`, `print`), full `roxygen2` documentation, CRAN compliance, and `reticulate` backend integration.
- **Enhanced Warm-Starting Matrix Exchange**:
  - In Stata (`pyquaidsce.ado`), exported `e(b_est)` ($1 \times K$ vector of free estimated structural parameters $\theta$) and `e(Sigma)` (residual covariance matrix) to `ereturn matrix`.
  - Stata `initial()` and `sigma_initial()` now seamlessly accept `e(b_est)` and `e(Sigma)` from prior runs.
  - In R (`rquaidsce`), exposed `fit$theta` and `fit$sigma` and added `initial` and `sigma_initial` arguments to `quaidsce()`.
- **Default Predictor Documentation Alignment**: Verified and aligned all documentation and tutorials to reflect `first_stage_predict="xb"` (theoretical textbook Shonkwiler & Yen linear index) as the primary default, while preserving `first_stage_predict="pr"` for legacy Stata compatibility.
- **Author Metadata & Contact**: Updated official author contact to `sinaamiri9000@gmail.com` across all packages, documentation, and metadata files.

- Renamed the `stop_rule` value `"stata"` to `"standard"` and made `"standard"`
  the default in both the Python API and the Stata command (previously
  `"tight"` in Python). The disjunctive Stata-matching behavior itself is
  unchanged; only the label and default moved.
- Changed `strict_stata` to default to `False` (corrected textbook formulas)
  in the Python API, the Stata bridge, and `pyquaidsce.ado`. Pass
  `strict_stata=True` (or `strict_stata(true)` in Stata) for exact replication
  of the original ado's elasticity calculations.
- Ported all remaining Python-only options into `pyquaidsce.ado`:
  `start()`, `initial()`, `sigma_initial()` (Stata matrix names),
  `vce_sigma()`, `tol()`, `nrtol_stop()`, `sigma_tol()`,
  `inner_nrtol_early()`, `max_iter()`, `max_outer()`, `chunk()`,
  `boot_sigma_tol()`, and a new `gnlog` switch mirroring `gn_verbose`.
- Updated `pyquaidsce.sthlp` with the new options and revised defaults.

## 1.3.0 — 2026-08-15

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
