# Changelog

All notable user-visible changes to `pyquaidsce` are documented here.

## Unreleased

- Added Python `latent_adding=False` to `quaidsce()` (QUAIDS/AIDs) and
  `translog()`. With S&Y, `True` imposes adding-up on latent shares through
  the free-parameter map, including QUAIDS control-function coefficients.
  All observed equations and selection-correction coefficients remain free.
  Bootstrap refits and nested/linear/mean starts retain the selected mode.
  Enabling this option without censoring is rejected; uncensored adding-up
  remains mandatory. Existing censored defaults remain relaxed.
- Added Shonkwiler--Yen to Python `translog(censor=True)`: all alphas are
  free, every equation is estimated, and latent/observed shares remain
  unnormalized. Gamma symmetry and demographic quantity translation are kept.
- Added independent participation designs, participation/latent predictions,
  complete S&Y expenditure and price derivatives, conditional analytical SEs,
  and full two-step serial/spawn bootstrap. Compensated S&Y output is labeled
  as a Slutsky transformation.
- Extracted the shared participation engine and censoring transformation from
  QUAIDS orchestration while preserving its existing APIs and behavior.
- Added optional `start="nested"` for S&Y demographic translation. It obtains
  initial values from a no-translation submodel using the same SUR engine;
  final model restrictions and standard GN/IFGNLS defaults are unchanged.
  Translog results now diagnose a large final scaled gradient instead of
  accepting a domain-limited tiny step as sufficient evidence of convergence.

- Added Python `translog()` for basic indirect translog demand,
  with demographic translation, analytic share/parameter derivatives,
  predictions for all goods, and full delta-method elasticity inference.
- Introduced a model backend interface for the existing nonlinear SUR engine,
  shared observation/output helpers, and shared spawn-safe bootstrap execution.
  Existing `quaidsce()` APIs and results retain their behavior.
- Added an independent zero-start Uruguay benchmark against the supplied
  `demandsys` output, including an explicit diagnostic for Stata's
  fixed-predicted-share elasticity standard errors.
- Added translog theory, derivative, missing-data, inference, and serial/spawn
  bootstrap checks. The current R and Stata entry points remain QUAIDS/AIDs
  interfaces; they do not yet expose the new Python translog entry point.

- Unified the main estimator name across interfaces: Python and R use `quaidsce()`, and Stata now uses the `quaidsce` command.
- Renamed the R package from `rquaidsce` to `pyquaidsce`. The Python package name remains `pyquaidsce`.
- Stata accepts both `quaidsce` and `pyquaidsce`. Both run the same estimator; `pyquaidsce` can be used to avoid a name conflict with an older installed `quaidsce` command.


## 1.6.0 — 2026-08-26

- **Breaking API cleanup:** removed `first_stage_predict` and `strict_stata` from Python, R, and Stata. The Shonkwiler-Yen first stage now always uses the Probit linear index, and elasticity calculations always use the corrected theoretical formulas.
- Removed the corresponding legacy branches from the Python estimator, elasticity engine, bootstrap forwarding, Stata bridge/ADO, R wrapper, examples, tests, and user documentation.
- Standardized the package-wide default estimation method to **IFGNLS** across `quaidsce()`, the exported low-level `nlsur()` helper, the R wrapper, the Stata command, examples, and documentation.
- Bumped Python, R, Stata, citation, and release metadata to 1.6.0 and updated the GitHub release workflow default tag.
- Reframed the archived Stata benchmark as historical compatibility/timing evidence so it is not mistaken for a fresh v1.6.0 exact-replication claim.
- Updated the master Python/R/Stata scenario suites so the former legacy-predictor scenario now verifies the implicit IFGNLS default.
- Changed the default IFGNLS outer tolerance to `sigma_tol=1e-5` and the bootstrap default to `boot_sigma_tol=1e-5`.
- Added `blas_threads`; the default is one BLAS thread per estimation process.

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
