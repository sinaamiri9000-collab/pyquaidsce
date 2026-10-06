# Analytical inference checks on small4

Original small4 implementation: `cf4987e4f7744c38342e81930d5302aedc912f26`. Analytical-gradient and separate-participation implementation: `da34954`. Both are based on main `636609f17e732a57b140cbbad5d2bf4042bc396a`.

The optional Python `analytic=True` calculation runs after the existing point estimate. It includes Probit, demand, estimated GLS-weight, and sample-mean uncertainty. IV and control functions are left for later work. See [Methodology](../../docs/methodology.md#analytical-inference) for the single mathematical definition.

## Small4 comparison

The public `small4.dta` has 2,000 observations. Shares are `sw1 sw2 sw4 sw9`, prices are `p1 p2 p4 p9`, expenditure is `total`, and demographics are `x1 x2`. All cases use `anot=10`, zero starts, GN, and the unchanged default tolerances. Bootstrap fits re-estimate the Probits and demand system and recalculate the sample means.

| Model | Successful bootstrap draws | Point fit (s) | Extra inference (s) | Median analytical / bootstrap elasticity SE | Range |
|---|---:|---:|---:|---:|---:|
| quadratic-ifgnls | 388/400 | 0.410 | 0.057 | 0.954 | 0.148–1.311 |
| quadratic-nls | 399/400 | 0.105 | 0.051 | 0.873 | 0.512–1.039 |
| quadratic-fgnls | 400/400 | 0.150 | 0.172 | 0.997 | 0.355–1.177 |
| linear-ifgnls | 393/400 | 0.231 | 0.039 | 1.100 | 0.618–1.248 |
| uncensored-ifgnls | Not compared | 0.112 | 0.026 | — | — |

The median and range cover all 36 elasticities, including income, uncompensated prices, and compensated prices. FGNLS inference time includes replaying the initial NLS stage. Times are from this machine, not general speed guarantees. Uncensored inference was checked without a bootstrap elasticity comparison.

Each comparison uses 400 requested draws, two workers, a 60-second replication limit, and the package's own bootstrap. Seeds are 17000, 17001, 17002, and 17003 in table order. The saved run used the supported `fork` context on Linux. All successful draws are retained: there is no trimming, clipping, or replacement of unusual estimates. Failure records are in `small4/summary.json`.

## Agreement and limits

The four SE medians are fairly close, but this does **not** establish agreement for every elasticity. In quadratic IFGNLS, the first income elasticity has analytical SE 0.25661 and bootstrap SE 1.73075. The adjusted share in its denominator is 0.06102 at the point estimate. Bootstrap draw 15 has an adjusted share near −0.000751 and an income elasticity near 30.3. The delta method is a local approximation and does not reproduce this near-zero-denominator tail.

NLS and FGNLS also have material differences for some elasticities. For example, FGNLS draw 355 has first adjusted share 0.005681 and first income elasticity −3.3172, compared with point adjusted share 0.06197. Selected draw checks are saved in `extreme-draws.json`. These observations explain a source of finite-sample disagreement; they do not explain every discrepancy or establish QUAIDS confidence-interval coverage.

The same 400 IFGNLS samples were also fitted with a warm start:

| Start | Successful draws | First income elasticity bootstrap SE | Median analytical / bootstrap SE |
|---|---:|---:|---:|
| zero | 388/400 | 1.73075 | 0.954 |
| warm | 390/400 | 1.72282 | 0.949 |

Draw 15 remains near 30.3 with both starts. Warm starts do not remove the tail. The cold-start vectors in this diagnostic match the original package bootstrap bit for bit.

## Implementation checks

All 69 tests passed, including 13 analytical-inference tests and six DGP checks. These check full demand scores against the separate full Jacobian, Probit scores against log-likelihood differences, the complete estimating-equation derivatives for all three methods, generated-mean derivatives, cross-stage covariance, dropped Probit columns, and step/chunk stability. DGP checks cover exact conditional means, positive purchased shares, complete budgets, separate correlated purchases, and the full estimating-equation derivative with an excluded selection covariate. Further cases cover linear and uncensored models, NLS with LM, a custom Probit design, and combined analytical/bootstrap output.

The solver, analytic solver Jacobian, demand model, parameter mapping, Probit optimizer, and elasticity formulas match main byte for byte. Point parameters, residual covariance, Probit coefficients, log-likelihood, and iteration counts are unchanged when analytical inference is enabled. The small4 comparison outputs record exact equality for every fitted case. `CITATION.cff` is unchanged.

## Independent two-step reference

A separate known-model experiment uses 500 independent datasets of 1,000 observations. A correctly specified Probit supplies a generated regressor to a normal linear outcome, with correlated shocks between stages. It compares the joint sandwich with the classical Murphy–Topel expression in Hardin (2002), equation (1). It also checks a target depending on a generated sample mean. This is a reference problem, **not** a QUAIDS simulation.

| Target | Empirical SD | Mean joint-sandwich SE | Mean Murphy–Topel SE | Coverage, joint / MT |
|---|---:|---:|---:|---:|
| coefficient | 0.09682 | 0.10031 | 0.10028 | 0.938 / 0.938 |
| mean_target | 0.19312 | 0.20559 | 0.20567 | 0.958 / 0.960 |

## Known-parameter censored QUAIDS Monte Carlo

A second experiment generates 400 independent datasets of 2,000 observations. Unlike the linear reference above, it uses the censored QUAIDS mean equations, one Ray demographic, nonzero selection corrections, and known demand and Probit parameters. The generator implements the QUAIDS equations separately from the package model.

For covariates X, let D be a shared Probit purchase event. Conditional on D=1, the four shares are w*(X) + delta phi(X tau)/Phi(X tau) + bounded mean-zero noise. Their unconditional mean is exactly Phi(X tau) w*(X) + delta phi(X tau), the fitted model. Purchased shares stay positive over the entire specified support; an outside category receives the budget when D=0. The data therefore contain actual zeros and close the budget without clipping or renormalizing generated observations. This is a controlled common-participation case, not a general simulation of separate purchase decisions.

The target is the package's existing reported elasticity function evaluated at the known parameters and **population** means. No elasticity formula is replaced for this experiment. Smooth population moments use deterministic Sobol integration with 2^20 points; repeating with 2^18 points changed any target by at most 0.000000598.

| Method, default settings | Successful datasets | Median mean SE / empirical SD | Ratio range, all elasticities | Median 95% coverage | Coverage range |
|---|---:|---:|---:|---:|---:|
| nls | 400/400 | 1.209 | 1.051–1.748 | 95.25% | 91.75–97.75% |
| fgnls | 400/400 | 1.054 | 0.978–1.119 | 95.75% | 91.75–97.50% |
| ifgnls | 398/400 | 1.063 | 0.979–1.118 | 95.73% | 91.46–97.49% |

FGNLS and IFGNLS give encouraging calibration in this case. The minimum coverages still fall below 95%; these medians do not establish accurate coverage for every individual elasticity or other DGPs. All coefficient and elasticity biases, SEs, empirical SDs, and coverages are saved in `quaids_dgp/*.csv`; `functional_name` identifies each elasticity by its actual (good, price) order. Failed fits are recorded and are not replaced.

NLS received limited extra accuracy checks. In dataset 385, the default fit stopped after five GN steps, with income elasticity 4 SE of 4.7405. Using the **existing** options `param_tol=1e-8`, `objective_tol=1e-10`, and `gn_tol=1e-8` took eleven steps and gave SE 0.07019. Across the same 400 datasets, these settings reduced its median mean-SE/empirical-SD ratio from 1.209 to 1.064. Separate 10,000-observation runs and the per-fit probe are retained in the NLS folders and `nls-accuracy-probe.json`. These are numerical diagnostics; no solver, tolerance default, or automatic polishing step is changed.

## Separate, correlated purchase decisions

The new generator uses a different Probit coefficient vector for each good and Gaussian purchase shocks with correlation 0.5 within a household. Households remain independent. A covariate `s` enters the Probits but not the demand equations. Each modeled good can be purchased separately; an outside good closes the budget of this partial basket. Purchased shares are `wstar_i + delta_i*phi_i/Phi_i + bounded mean-zero noise`, so each conditional mean is exactly the package's censored QUAIDS equation. Support bounds guarantee positive purchased shares and outside shares without clipping or renormalization. The probability of purchase ranges from 0.399 to 0.605 across goods.

There are 400 independent datasets **per scenario**, with unchanged default tolerances, two workers, zero starts and GN. The 2,000-observation case matches the size of small4. The zero-correlation comparison shares seeds and covariates with its correlated counterpart. Sobol integration with 2^20 versus 2^18 points changes any population target by at most 0.00000242.

| N | Purchase-shock correlation | Method | Successful inference | Median mean SE / empirical SD | Ratio range | Median 95% coverage | Coverage range |
|---:|---:|---|---:|---:|---:|---:|---:|
| 500 | 0.5 | FGNLS | 398/400 | 1.052 | 0.915–1.199 | 95.48% | 91.96–97.74% |
| 500 | 0.5 | IFGNLS | 391/400 | 1.030 | 0.909–1.151 | 95.52% | 91.82–97.70% |
| 2,000 | 0.5 | FGNLS | 400/400 | 1.014 | 0.934–1.151 | 95.25% | 92.50–97.00% |
| 2,000 | 0.5 | IFGNLS | 388/400 | 1.002 | 0.937–1.072 | 95.10% | 91.49–97.16% |
| 8,000 | 0.5 | FGNLS | 400/400 | 1.018 | 0.967–1.175 | 95.00% | 93.25–97.00% |
| 8,000 | 0.5 | IFGNLS | 394/400 | 1.011 | 0.947–1.048 | 94.92% | 93.40–96.70% |
| 2,000 | 0 | FGNLS | 400/400 | 1.048 | 0.950–1.355 | 95.50% | 93.50–98.50% |
| 2,000 | 0 | IFGNLS | 383/400 | 1.013 | 0.930–1.077 | 95.04% | 93.21–98.43% |

Coverage is conditional on obtaining a converged fit and a full-rank inference system. All 46 failures are recorded, including three rank failures; none is replaced. `studentized-tails.csv` also gives the fraction of all 400 requested fits that produce a covering interval. Individual Monte Carlo Wilson intervals and sampling errors are in the per-elasticity CSVs. Around 95% coverage, 400 draws give a Monte Carlo standard error of about 1.1 percentage points. These are pointwise checks, not simultaneous confidence claims across all elasticities.

The fourth income elasticity illustrates a finite-sample limit. With correlated purchases and N=2,000, its IFGNLS coverage is 91.49%; its Wilson interval is below 95%, and most missed intervals fall on the upper side of the studentized distribution. At N=8,000 its coverage rises to 95.18%, its SE ratio is 1.012, and its bias falls from 0.02458 to 0.00482. This supports improving large-sample behavior in this DGP, rather than treating the original coverage loss as only Monte Carlo noise.

Every replication saves the full 36-by-36 elasticity covariance. All 630 pairwise elasticity differences are checked with their full covariance, using contrasts specified before examining results. Their median SE ratios range from 0.999 to 1.062 across methods/scenarios; their lowest coverages range from 90.46% to 92.25%. In IFGNLS, the normalized error between mean analytical and empirical covariance falls from 0.077 at N=500 to 0.050 at N=8,000. The covariance and contrast CSVs retain every entry; no cutoff is used to declare an overall pass.

FGNLS has larger isolated SEs than IFGNLS. In the zero-correlation N=2,000 scenario, its fourth income elasticity has empirical SD 0.11958, mean SE 0.16208 and root-mean-square SE 0.34049. Replication 133 gives SE 4.11745. The existing derivative function checked numerically gives the same SE, and its adjusted share is 0.02023. Its observed final demand information has a negative minimum eigenvalue, so the stopping flag does not certify a well-resolved local minimum for inference. In correlated N=8,000 replication 93, FGNLS SE is 1.60603; ten-times-stricter existing tolerances reduce it to 0.50293, while still stricter settings fail the initial-stage convergence check. This is a limited diagnostic, not a general cure. No polishing step, covariance replacement, solver change, tolerance-default change, or sample deletion is made. These cases remain a limitation of FGNLS analytical inference at the unchanged stopping settings.

Elasticity gradients are now analytical in the post-estimation layer. On small4, their maximum scaled difference from independent central differences is below 0.00000000034, and the largest relative SE difference is below 0.00000091. Median gradient time falls from about 15.6 ms to 0.34 ms over nine evaluations on this machine; this measures the gradient only. Full estimating-equation checks with distinct Probits and the excluded covariate have maximum scaled error below 0.00000000084. The original elasticity function and optimizer remain unchanged.

`separate_participation/` contains the comparison, failure and tail tables, per-target and per-contrast calibration, full covariance comparisons, all successful replication vectors and covariance matrices, source hashes and logs. `diagnostics.json` reproduces the two selected FGNLS extremes and the limited tolerance probe. `elasticity-derivative-check.json` records the small4 numerical comparison and timings.

## Reproduce

Run from the repository root after installing the package and its dependencies:

```bash
python -m unittest discover -s tests -v
python tools/check_analytical_small4.py --bootstrap-reps 400 --n-jobs 2 --rep-timeout 60 --mp-context fork
python tools/diagnose_analytical_small4.py --reps 400 --n-jobs 2 --timeout 60
python tools/check_two_step_reference.py --observations 1000 --repetitions 500 --seed 17042
python tools/simulate_analytical_quaids.py --observations 2000 --repetitions 400 --seed 17043
python tools/simulate_analytical_quaids.py --observations 2000 --repetitions 400 --seed 17043 --methods nls --param-tol 1e-8 --objective-tol 1e-10 --gn-tol 1e-8 --output validation/analytical_se/quaids_dgp_nls_accuracy
python tools/check_elasticity_derivatives.py
python tools/simulate_separate_participation.py --observations 500 --repetitions 400 --correlation .5 --output validation/analytical_se/separate_participation/n500-corr05
python tools/simulate_separate_participation.py --observations 2000 --repetitions 400 --correlation .5 --output validation/analytical_se/separate_participation/n2000-corr05
python tools/simulate_separate_participation.py --observations 8000 --repetitions 400 --correlation .5 --output validation/analytical_se/separate_participation/n8000-corr05
python tools/simulate_separate_participation.py --observations 2000 --repetitions 400 --correlation 0 --output validation/analytical_se/separate_participation/n2000-corr0
python tools/diagnose_separate_participation.py
python tools/summarize_separate_participation.py
```

On Windows, omit `--mp-context fork` to use the package's default process context. The diagnostic script uses the same sampling seeds and point estimator in persistent workers, with cooperative per-fit deadlines; the primary comparison uses the original bootstrap watchdog.

`small4/*.csv` contains estimates, both SEs, and their ratios. The NPZ files retain complete successful bootstrap vectors and covariance matrices. JSON files record versions, dataset/source hashes, settings, and failures; `comparison.csv` is a compact table. The reference folder contains all 500 linear-reference simulation results. For the additional NLS sample-size checks, use the last two Monte Carlo commands with `--observations 10000` and separate output folders.
