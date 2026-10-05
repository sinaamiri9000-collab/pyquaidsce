# Validation of pyquaidsce 1.6.0 using small4

This directory audits the **unchanged version 1.6.0**, starting from commit
`93ea81e4293ea4d2f0db899588079e1acc22c745`. It does not use archived Stata
agreement as an oracle. It contains the executable design and public evidence,
including failed contracts; it is not a claim that all options have passed.

## Reproduce

From the repository root, with Python >= 3.9 and the package dependencies:

```bash
python -m pip install .
python -m validation.run --large-reps 100 --output validation/results
```

The command exits nonzero when an assertion fails or an ordinary fit does not
converge. All remaining cases still run and evidence is written. Multiprocessing
is guarded by `if __name__ == '__main__'`; the parallel default is tested against
explicit `spawn`. R and licensed Stata runs are supplemental, not assumed.

Use a fresh output directory when selecting a subset so that old files are not
mistaken for current evidence:

```bash
python -m validation.run --groups core inputs --output /tmp/small4-core
python -m validation.run --only invalid-tol-2 --output /tmp/small4-repro
python -m validation.run --contract-only public-probit_no_constant --output /tmp/small4-helper
```

`--only` includes reference fits automatically. The complete matrix is:

`--contract-only` selects helper/output/interface checks and their prerequisite
fits. Pilot harness corrections fixed one overly specific rejection-message
match and applied the already declared fit-equivalence tolerance to separate
optimized fits (rather than the algebra-identity tolerance). The final audit
uses the fixed published thresholds; no estimator source was changed.

- 72 core fits: six valid model specifications × three estimators × two
  algorithms × two starting schemes.
- 78 additional fits: 18 input/sample cases, 18 first-stage designs, 12 control
  function/IV cases, 18 numerical-control cases, 12 interaction cases.
- 12 bootstrap smoke batches with eight draws each.
- Four larger bootstrap batches with **100** draws each: base and IV, each
  at seeds 16001 and 16002. No batch silently defaults back to 500.
- 100 invalid-input cases and four timeout/failure lifecycle cases.
- Public-helper, postestimation, methodological, source-contract, and real
  Python interface-transport checks, plus the pre-existing unittest suite.

## Independent checks

`oracles.py` declares tolerances before running. It reconstructs fitted shares
from the published equations, checks residual covariance/log likelihood,
compares the complete free-parameter Jacobian against finite differences on
deterministic small4 rows, and checks observation-level price/expenditure
derivatives. A separate SciPy optimizer checks the baseline Probits. QR OLS
checks the IV reduced form, covariance, residual orthogonality and F test.

Elasticity finite differences use model-consistent evaluation points. The public
reported elasticities use observed mean shares and separately averaged Probit
moments; they are not silently equated to average fitted elasticities. Adding-up
restrictions are tested only in branches where they are imposed. Negative fitted
shares are recorded as diagnostics and are not clipped or automatically treated
as software failures.

Conditional covariance checks distinguish the implemented inverse GLS normal
matrix from the methodological NLS sandwich check. Generated-regressor
uncertainty still requires the full bootstrap. Analytical elasticity SEs remain
undefined (`NaN`) by the current Python contract. A small bootstrap's covariance
can be singular and is checked for positive semidefiniteness, not invertibility.

## Evidence

- `results/environment.json`: versions, source/data hashes, BLAS and available
  process contexts, and the explicit 100-draw configuration.
- `results/manifest.json`: each scenario's settings and expected behavior.
- `results/scenarios.csv` and `results/checks.csv`: honest status and each assertion.
- `results/argument-coverage.csv`: all 41 public main-estimator inputs and mapped cases.
- `results/public-api-coverage.csv`: public exports and their evidence sources.
- `results/cases/`: per-case settings, diagnostics, estimates, SEs, free-parameter
  covariance, residual covariance, and failed-check evidence.
- `results/coefficients/`: coefficient/SE tables for returned fits.
- `results/bootstrap/`: successful draws in deterministic package order and their
  active covariance matrices. The package does not expose successful draw IDs
  when draws fail; `draw_number` is the **successful-row ordinal**, not an
  invented original replication ID. Failure messages retain the original IDs.
- `results/bootstrap-stability.csv`: SE drift across the two seeds, diagnostic only.
- `results/bootstrap-failures.csv`: original failed replication IDs and reasons.
- `results/optimizer-comparisons.csv`: same-model algorithm/start diagnostics,
  including separate tight-stop follow-ups.
- `results/core-matrix.svg`: optional vector figure generated with matplotlib.
- `results/existing-tests.csv` and `.log`: the existing suite's individual outcomes.
- `results/REPORT.md` and `FINDINGS.md`: scorecard and interpreted findings.

`FINDINGS.fa.md` provides a Persian explanation. `reproduce.py` rebuilds a
natural failed draw with its original sample indices, seed, and solver settings.
Supplemental tests cover confidence-level boundaries, additional process/seed/
timeout paths, noncooperative worker termination, tighter solves, and alternative
selection starts. They use new IDs and never overwrite an original failed case.

`PASS`, `EXPECTED_ERROR`, `FAIL`, `NONCONVERGED`, `BLOCKED`, and `NOT_RUN` have
different meanings. No case is green merely because the estimator returned an
object. Deliberate low-iteration and injected-failure cases state their expected
outcome explicitly. Neither seeds nor numerical thresholds are changed to hide
a failure. The source estimator and `bench/small4.dta` remain unchanged.

The broad suite checks functionality on small4; it does not prove parameter
recovery, confidence-interval coverage, or economic validity of synthetic test
instruments. Existing synthetic/theory tests provide complementary evidence.
