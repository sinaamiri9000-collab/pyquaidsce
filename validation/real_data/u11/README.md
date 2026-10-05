# U11: standard versus nlsur stopping predicates

This directory records the **same-engine U11 experiment** on the author's real
data. Both fits use pyquaidsce **1.6.0**, zero initialization, identical prepared
observations and shared first-stage Probits. U11 has **no curvature constraint**.
The second fit changes convergence predicates only; it does not use the QR
backend from the separate small4 architecture study.

Read the [Persian report](REPORT.fa.md), [comparison CSV](results/comparison.csv)
and [complete model specification](specification/model_specification.json).

| Fit | Converged | Outer stages, including NLS/FGNLS | IFGNLS updates | GN steps | Seconds |
|---|---|---:|---:|---:|---:|
| standard | Yes | 48 | 46 | 79 | 117.001 |
| nlsur predicates only | Yes | 17 | 15 | 565 | 822.627 |

The largest full-coefficient difference is `2.10228159174e-5`; the largest fitted
share difference is `4.00423140756e-7`; the largest elasticity difference at the
shared mean point is `2.42934271610e-6`. The likelihood difference is
`-6.30971044302e-8`. The single-run time ratio is approximately `7.03`.
These findings concern this U11 specification and dataset, rather than every
model or package option. Neither stopping contract proves a global optimum.

## Experimental contract

- 37,658 observations, 11 goods, 8 demographics and 198 free demand parameters.
- Quadratic censored demand system, `anot=11`, expenditure control-function
  correction and shared selection equations, as specified in the supplied model.
- Initial NLS with identity covariance, then initial FGNLS and the original
  IFGNLS covariance updates. No warm parameter start, linear start or bootstrap.
- GN, observation chunk size 15,000, one BLAS thread, inner/outer limits 400.
- Direction calculation, analytic derivatives, accumulation, whitening,
  scaling, scaled solve, strict decrease acceptance, halving and damping remain
  unchanged. The AST audits and diffs document the precise allowed differences.
- The experimental inner predicate requires both accepted parameter change
  and objective change to satisfy the uploaded `nlsur.ado` rules, with
  `eps=1e-5`, `tau=1e-3`. The outer predicate uses its parameter OR residual
  covariance change rule, with `ifgnlseps=1e-10`.
- An original-engine failure to accept a step still exits the inner solver.
  The experimental convergence flag then reports False instead of using the
  former gradient predicate. No synthetic zero step or new acceptance rule is
  added. All inner solves succeeded in both completed runs.
- This is not a native Stata run or a replacement package optimizer. The first
  11 accepted NLS parameter vectors are bitwise identical between the runs.

Times include the demand numerical fit, final parameter covariance assembly and
common recording hooks. Shared data preparation and selection Probits are
excluded. There is one fresh numerical fit per backend, not a repeated timing
benchmark. The new standard run bitwise reproduces the previously cached U11
coefficients. All 30 recorded numerical and source checks pass; this does not
replace each fit's convergence outcome.

## Published evidence

- `results/`: estimates, comparison and phase CSVs, full coefficient comparisons,
  residual covariances, free-parameter covariance matrices, elasticity CSVs,
  accepted-step and outer histories, parameter-state NPZs, source diffs and audits.
- `specification/`: exact model and selection layout, shared Probit coefficients
  and diagnostics, elasticity evaluation means, free-parameter order, versions,
  aggregate commodity mapping, fixed weights, price normalization, reduced-form
  coefficients, purchase support and price-imputation counts.
- `logs/`: original solver and comparison-output logs containing numerical
  stage/step diagnostics.
- `run_comparison.py`, `build_report.py`: exact scripts used in the experiment.
  Their original execution hashes are in `results/experiment_manifest.json`.
- `ORIGINAL_BUNDLE_HASHES.json`: inventory of the original local ZIP. Its report
  hash refers to the report before the publication-status sentence was updated;
  it is not the inventory of this directory.
- `FILE_HASHES.json`: SHA256 inventory of the published files, excluding itself.
- `PUBLICATION_MANIFEST.json`: publication scope and original-artifact identity
  checks. The execution commit recorded in the experiment manifest is the
  source snapshot used for estimation, not this evidence-publication commit.

The NPZ histories contain parameter vectors and aggregate residual covariance
matrices, **not observation-level arrays**. Raw household CSVs, identifiers,
prepared household datasets and data-bearing pickles are not included.

## Reproduce with the private inputs

Re-estimation requires the author's original `cldenew.csv`, supplied C11/U11
preparation bundle and the private `common_inputs.pkl` / `auxiliary_u11.pkl`
created during that preparation. Those inputs are deliberately not provided by
this public evidence directory. Their source/data hashes and prepared demand
array fingerprints are recorded in the manifest. Adjust `REPO` and `C11` near
the top of the two original scripts for a different filesystem location.

Use the recorded package/source snapshot and numpy, scipy, pandas and
threadpoolctl versions, then run:

```bash
python run_comparison.py standard
python run_comparison.py nlsur_criteria
python build_report.py
```

The first two commands perform fits and overwrite local result files. The third
reads saved numerical fits and generates the comparisons. It also requires the
private prepared inputs for fitted-share evaluation. No new estimation was run
to publish these files; the original numeric evidence was copied unchanged.

To check the published evidence **without the private data**:

```bash
python validation/real_data/u11/verify_evidence.py
```
