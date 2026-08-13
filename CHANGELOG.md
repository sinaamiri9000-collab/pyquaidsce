# Changelog

All notable user-visible changes to `pyquaidsce` are documented here.

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

## 1.0.0

- Initial Python implementation of the censored QUAIDS workflow.
