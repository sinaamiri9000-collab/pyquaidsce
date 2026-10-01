# Uncensored translog: Uruguay versus Stata demandsys

The supplied data have 6,848 rows, 14 food goods, and seven demographics.
The estimation sample has 6,573 rows after complete-case selection. This
benchmark estimates basic translog with demographic translation and IFGNLS.
The original comparison is uncensored. A separate `run_sy.py` now exercises
Shonkwiler--Yen with relaxed latent adding-up; see below. Generalized translog
and expenditure endogeneity correction remain outside this implementation.

```bash
python -m pip install -e .
python benchmarks/translog_uruguay/run_python.py
```

The primary fit starts at **zero**, with the public GN/IFGNLS defaults. The
printed Stata coefficients are never used as optimizer starting values.
The source CSV and Stata output are stored in `data/` and `reference/`.
`parse_reference.py` builds the machine-readable reference, including small
coefficients printed in scientific notation.

## Data preparation

The runner exactly reproduces the supplied cleaning: replace negative
`w1_red` through `w14_red` values by zero (1,808 replacements in `w8_red`),
then select `P_med1`--`P_med14` if `P_med1` exists, otherwise lower-case
`p_med1`--`p_med14`. It reads Stata's `.` CSV missing values as NaN. Shares
are not renormalized after cleaning.

## Postestimation and inference conventions

The original command `estat elasticities, uncompensated atmeans` reports
6,848 rows, while estimation reports 6,573. Matching its elasticity point
estimates requires price means on complete price vectors, expenditure means on
the full available expenditure sample, and demographic means on complete
demographic vectors in the full frame. The runner uses
`fit.elasticities_at_means(cleaned_full_frame)` explicitly. Ordinary estimator
results continue to use the estimation sample by default.

Stata's printed elasticity standard errors are reproduced by holding the
predicted-share denominator fixed in the parameter delta method. Public Python
SEs instead differentiate the whole elasticity, including that denominator.
Both appear in the CSV:

- `python_se`: public full delta-method inference.
- `python_stata_conditional_se`: a reference diagnostic with the denominator fixed.
- `se_difference` and `stata_conditional_se_difference`: their separate errors.

The public/reference SE difference is reported separately in JSON and is **not**
claimed to vanish. This diagnostic never changes the estimator's public covariance.

## Stored results and tests

`results/comparison_summary.json` records sample sizes, convergence, runtime,
input hash, maximum errors, declared tolerances, and the inference distinction.
`stata_python_comparison.csv` compares 217 structural coefficients and 196
Marshallian elasticities. `python_results.npz` stores the parameter and
covariance arrays, plus all three elasticities. Expenditure and Hicksian
elasticities are checked against demand theory because the supplied Stata
text contains only Marshallian elasticities.

The coefficient and elasticity comparisons are approximate, within declared
tolerances; the reference is printed at finite precision. The log-likelihood
agrees when rounded to the two decimal places displayed by Stata.

Fast external-reference evaluator checks run in the normal test suite. For the
complete independently initialized refit:

```bash
PYQUAIDSCE_RUN_URUGUAY=1 python -m unittest tests.test_translog_reference -v
```

The existing CI matrix exercises Python 3.9--3.13 on Linux and Windows. Local
validation in this task runs Linux, including the shared `spawn` bootstrap;
it does not claim that a Windows CI job was executed here.

## Shonkwiler--Yen with relaxed latent adding-up

```bash
python benchmarks/translog_uruguay/run_sy.py
# Diagnostic three-replication check, not publication-quality bootstrap inference:
python benchmarks/translog_uruguay/run_sy.py --reps 3 --jobs 2
```

This uses the same cleaning and 6,573 observations, all 14 equations, all seven
demographic translations, and independent per-good Probits. All alphas are
free. The denominator constant stays one and fitted shares are not normalized.

The runner selects `start="nested"` while retaining the standard IFGNLS,
Gauss--Newton and stopping settings. Starting at zero with these demographics
can reach the boundary of positive effective expenditure before the gradient
is small. The nested start first fits the model without demographic translation
using identical participation terms, then fits the complete model with every
translation coefficient free. Its final minimum effective-expenditure ratio
is about 0.219, and its scaled gradient is below 1e-5. No observations are
dropped to get that result and no external coefficient estimates are used.

`results_sy/` stores coefficients, elasticities, covariance, sample mask,
participation terms and JSON diagnostics. Input-level finite differences check
the final quantity elasticities through the public prediction function.
The supplied Stata output is **not** a reference for censored translog; only
the original uncensored comparison supports a Stata reproduction claim.

Some predicted shares are negative and several equation R-squared values are
negative in the unrestricted SUR fit. They are reported in the diagnostic
JSON rather than clipped. Convergence alone does not establish economic
regularity or good fit for every equation.
