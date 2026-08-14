# Changelog

All notable user-visible changes to `pyquaidsce` are documented here.

## 1.1.0 — 2026-08-14

Stata integration release and comprehensive documentation expansion.

- **Stata Integration**: Added official `pyquaidsce.ado` and `pyquaidsce.sthlp` commands allowing Stata users to estimate censored QUAIDS directly inside Stata via the Python computational engine (`stata_bridge.py`), supporting full `e()` matrices, scalars, and postestimation commands (`test`, `lincom`, `outreg2`, `esttab`).
- **Stata Package Management**: Added `stata.toc` and `pyquaidsce.pkg` enabling direct 1-line installation in Stata via `net install`.
- **User Guide & API Reference**: Created comprehensive [User Guide](docs/user-guide.md) documenting every input parameter, optimizer setting, and attribute of `QuaidsceResults` with practical extraction recipes.
- **Documentation Refinement**: Simplified all documentation files to clear applied economics terminology.
- **Motivation & Citation**: Documented the empirical motivation for the package and updated citation metadata to 2026.

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
