# pyquaidsce 1.6.0 small4 validation results

This is an audit of the unchanged 1.6.0 estimator, not a comparison with the original Stata ado.

The interpreted findings are in [FINDINGS.md](../FINDINGS.md) and [FINDINGS.fa.md](../FINDINGS.fa.md).

- Scenarios: 322; statuses: `{'PASS': 177, 'NONCONVERGED': 17, 'FAIL': 42, 'EXPECTED_ERROR': 84, 'NOT_RUN': 2}`.
- Assertions: 10998; failed: 71.
- Existing suite: `{'tests_run': 42, 'failures': 0, 'errors': 0, 'skipped': 0, 'successful': True}`.
- Public arguments with mapped executions: 41/41.
- Full selected design executed: True.
- Release validation gate passed: **False**.

## How to read the evidence

`manifest.json` freezes expectations; `scenarios.csv` is the scorecard; `checks.csv` contains observed errors and fixed thresholds. `cases/` contains configurations, diagnostics, coefficient/SE vectors, theta, Sigma, and structural covariance. `coefficients/` provides estimates in CSV. `bootstrap/` contains all successful draws and active covariance matrices. Analytical elasticity SEs are encoded as `NaN` by contract. Small-bootstrap covariance is allowed to be singular.

No tolerance was increased to turn a failed case green. A returned nonconverged result is not counted as PASS unless that exact case deliberately tests exhaustion. Injected lifecycle failures are explicitly labeled.

`optimizer-comparisons.csv` contrasts algorithms/starts for the same specification. Its difference labels are diagnostics, not additional pass/fail assertions: default standard stopping is loose, and FGNLS weights depend on the preliminary NLS solution. Separate tight-stop follow-ups test stationarity and retain the original failures.

## Failed or incomplete cases

| Case | Status | Failed checks / reason |
|---|---|---|
| core-cq-ifgnls-lm-zero | NONCONVERGED | convergence_contract |
| core-cq-ifgnls-lm-linear | NONCONVERGED | convergence_contract |
| core-ca-ifgnls-lm-zero | NONCONVERGED | convergence_contract |
| core-uq-ifgnls-lm-zero | NONCONVERGED | convergence_contract |
| core-uq-ifgnls-lm-linear | NONCONVERGED | convergence_contract |
| selection-no_prices | NONCONVERGED | convergence_contract |
| selection-no_expenditure | NONCONVERGED | convergence_contract, V_symmetric |
| selection-none_logged | NONCONVERGED | convergence_contract, V_symmetric |
| selection-false_logged | NONCONVERGED | convergence_contract, V_symmetric |
| selection-intercept_only | FAIL | unexpected_exception |
| selection-independent_no_exp | FAIL | unexpected_exception |
| selection-price_only | NONCONVERGED | convergence_contract, V_symmetric |
| interaction-iv_lm_linear | NONCONVERGED | convergence_contract |
| interaction-linear_lm_selection | NONCONVERGED | convergence_contract |
| interaction-cf_lm | NONCONVERGED | convergence_contract |
| boot-uncensored | FAIL | no_unexplained_bootstrap_failures |
| boot-large-base-16001 | FAIL | no_unexplained_bootstrap_failures |
| boot-large-base-16002 | FAIL | no_unexplained_bootstrap_failures |
| boot-large-iv-16001 | FAIL | no_unexplained_bootstrap_failures |
| boot-large-iv-16002 | FAIL | no_unexplained_bootstrap_failures |
| invalid-tol-2 | FAIL | invalid_input_rejected |
| invalid-tol-3 | FAIL | invalid_input_rejected |
| invalid-nrtol_stop-2 | FAIL | invalid_input_rejected |
| invalid-nrtol_stop-3 | FAIL | invalid_input_rejected |
| invalid-sigma_tol-2 | FAIL | invalid_input_rejected |
| invalid-sigma_tol-3 | FAIL | invalid_input_rejected |
| invalid-inner_nrtol_early-0 | FAIL | invalid_input_rejected |
| invalid-inner_nrtol_early-1 | FAIL | invalid_input_rejected |
| invalid-inner_nrtol_early-2 | FAIL | invalid_input_rejected |
| invalid-inner_nrtol_early-3 | FAIL | invalid_input_rejected |
| invalid-boot_sigma_tol-0 | FAIL | expected_exception_type |
| invalid-boot_sigma_tol-1 | FAIL | expected_exception_type |
| invalid-boot_sigma_tol-2 | FAIL | expected_exception_type |
| invalid-boot_sigma_tol-3 | FAIL | invalid_input_rejected |
| invalid-reps-0 | FAIL | invalid_input_rejected |
| invalid-reps-1 | FAIL | invalid_input_rejected |
| invalid-n_jobs-0 | FAIL | invalid_input_rejected |
| invalid-n_jobs-1 | FAIL | invalid_input_rejected |
| invalid-n_jobs-2 | FAIL | expected_exception_type |
| public-probit_no_constant | FAIL | unexpected_exception |
| public-reduced_form_labels | FAIL | mismatched_labels_rejected |
| output-summary-no_elasticities-bootstrap | FAIL | elasticity_rows_suppressed |
| output-invalid_level-ci-0 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-percentile_ci-0 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-summary-0 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-ci-1 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-summary-1 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-ci-2 | FAIL | invalid_confidence_level_rejected |
| output-invalid_level-summary-2 | FAIL | invalid_confidence_level_rejected |
| methodology-NLS_covariance | FAIL | conditional_multivariate_NLS_covariance, NLS_covariance_SE_ratio |
| numeric-standard_tolerance_contract | FAIL | standard_uses_public_tolerances |
| followup-tight-ifgnls-lm-zero | NONCONVERGED | convergence_contract, independent_stationarity |
| followup-tight-ifgnls-lm-linear | NONCONVERGED | convergence_contract, independent_stationarity |
| followup-linear-selection-no_prices | NONCONVERGED | convergence_contract |
| followup-linear-selection-no_expenditure | NONCONVERGED | convergence_contract, V_symmetric |
| followup-linear-selection-intercept_only | FAIL | V_symmetric, V_positive_semidefinite, sigma_numerical_rank, gaussian_log_likelihood, analytic_nonelasticity_se_finite |
| followup-linear-selection-independent_no_exp | FAIL | unexpected_exception |
| R-static-fitted_backend | FAIL | R_backend_method_exists |
| R-static-analytic_SE | FAIL | R_summary_preserves_undefined_SE |
| R-native | NOT_RUN | native execution is supplemental and no native runtime was used |
| Stata-native | NOT_RUN | native execution is supplemental and no native runtime was used |

## Bootstrap precision

The four large batches use 100 draws each, as requested. `bootstrap-stability.csv` compares SEs across two fixed seeds for the base and IV models. These are diagnostics, not a proof of frequentist coverage or instrument validity. No threshold on SE drift is retroactively chosen to declare statistical convergence.

A full-success bootstrap stress gate is deliberately stricter than the documented ability to discard failed draws. A failed stress gate does not imply that covariance/SE accounting was implemented incorrectly; the individual assertions distinguish the two. Successful-row ordinals are not original failed/successful replication IDs.

## Scope

small4 provides 2,000 observations, four goods, and two demographics. Subsystems and derived columns are recorded; the original DTA is unchanged. Existing synthetic/mathematical tests supplement the empirical matrix. Native R and licensed Stata execution are not represented as tested. This run does not establish behavior on every dependency version, platform, instrument design, or number of goods.
