# Cluster inference: four-good real-data comparison

This check uses the user's supplied `cldenew.csv`, not the later frozen
`cldenewfinal.csv` described in the article brief. The raw-file SHA256 is
`9819184dc68b471247253c74c1f257b98dfc983d3cb6f45218390ea0b49e9fb1`.
Household data and identifiers are not included here.

## Model and sample

- Urban households, seasons 3 and 4: 9,747 rows before the basket filter.
- Groups 6, 7, 8 and 9: fruits/vegetables, non-rice bread, protein, dairy/eggs.
- Expenditure is the sum for these four groups; shares divide group expenditure
  by that sum. Three zero-basket households are excluded, leaving 9,744 rows.
- 1,792 survey clusters identified by the nine-digit household-address prefix;
  cluster sizes range from 1 to 6. This does not establish an official PSU definition.
- Quadratic censored demand, `anot=11`, GN and zero starts, without curvature
  constraints, IVEXP or an external control function.
- Both stages include age, household size and number of children. Region
  dummies are excluded from both stages, as requested by the user.
- The four prices come from the supplied U11 full-sample Majumder / Geometric
  Young preparation. They are not renormalized after subsetting. Both methods
  use these fixed prepared prices; the bootstrap resamples their household rows.
- Probits and demand are unweighted, as in the supplied U11 point fit. This is
  a cluster covariance comparison, not full design-based survey inference.

The package's four default tolerances are unchanged. The iteration limits are
400/400, with chunks of 2,000 and one BLAS thread per fit. Bootstrap uses two
workers and seed 20261006. The first 50 requested draws are retained when
extending to 100; failed draws are recorded and are not replaced.
All Probits, demand parameters, residual covariances and sample means are
re-estimated. Prices are held fixed and there is no expenditure reduced form.

## Results

Ratios below are analytical cluster SE / bootstrap SE over all 36 elasticities.
The analytical covariance uses `cluster_correction=True` (G/(G-1)).

| Method | Requested / successful draws | Median ratio | Ratio range | Bootstrap time (s) |
|---|---:|---:|---:|---:|
| FGNLS | 50 / 43 | 1.112 | 0.799–1.341 | 19.7 |
| FGNLS | 100 / 87 | 1.084 | 0.863–1.386 | 39.7 |
| IFGNLS | 50 / 43 | 1.671 | 0.792–2.344 | 50.9 |
| IFGNLS | 100 / 87 | 1.677 | 0.801–2.375 | 102.5 |

Times are measured on this two-worker machine. Single point fits take roughly
0.67 seconds for FGNLS and 2.31 seconds for IFGNLS, excluding inference.
Extra analytical inference takes 0.50 and 0.21 seconds respectively; FGNLS
includes the initial NLS replay. Point counts are 2 outer / 27 inner for FGNLS,
and 21 outer / 91 inner for IFGNLS. Bootstrap median counts are 2 / 28 and 24 / 81.

For the 100-draw comparison, expenditure elasticities and their SEs are:

| Method | Group | Estimate | Analytical cluster SE | Bootstrap SE |
|---|---:|---:|---:|---:|
| FGNLS | 6 | 0.772957 | 0.012776 | 0.010703 |
| FGNLS | 7 | 0.312451 | 0.017103 | 0.019165 |
| FGNLS | 8 | 1.442893 | 0.014661 | 0.010578 |
| FGNLS | 9 | 0.623567 | 0.013569 | 0.011870 |
| IFGNLS | 6 | 0.748085 | 0.069310 | 0.048879 |
| IFGNLS | 7 | 0.594880 | 0.332107 | 0.154108 |
| IFGNLS | 8 | 1.487674 | 0.036002 | 0.044932 |
| IFGNLS | 9 | 0.495861 | 0.152996 | 0.082323 |

FGNLS shows closer agreement overall, with material differences for individual
elasticities. IFGNLS does not show satisfactory agreement in this sample.
These results do not establish confidence-interval coverage. With 43 and 87
successful draws, bootstrap SEs are still imprecise; unsuccessful draws can
also affect the comparison. No targets or replications are trimmed.

For IFGNLS, define the requested statistic as estimate / SE (null value zero).
The percentage ratio is bootstrap statistic / analytical statistic, times 100;
the same ratio applies for any shared null with a nonzero numerator. Thus
168% means a statistic 68% larger, while 80% means one 20% smaller.
Using the 87 successful fits from 100 requested draws:

| Reported block | Median percentage ratio | Range |
|---|---:|---:|
| Full demand coefficients (41 entries) | 84.4% | 0.192–204.9% |
| Probit coefficients (36 entries) | 92.2% | 61.4–110.3% |
| Elasticities (36 entries) | 167.7% | 80.1–237.5% |

| Expenditure elasticity | Analytical statistic | Bootstrap statistic | Percentage ratio |
|---|---:|---:|---:|
| Group 6 | 10.793 | 15.305 | 141.8% |
| Group 7 | 1.791 | 3.860 | 215.5% |
| Group 8 | 41.322 | 33.110 | 80.1% |
| Group 9 | 3.241 | 6.023 | 185.8% |

The full coefficient range shows that agreement is also poor for some demand
parameters; the elasticity median alone does not summarize the whole model.
`ifgnls-50-t-statistics.csv` and `ifgnls-100-t-statistics.csv` retain every
coefficient and elasticity, both statistics, their percentage ratio and
percentage change. Undefined zero/zero cases, if present, remain missing.

## What the diagnostics establish

Every remaining bootstrap failure is in the group-7 Probit. The final sample
has only 28 zeros for that group, versus 38, 351 and 51 for groups 6, 8 and 9.
A linear-programming support check finds a separating direction in all 13
failed draws, with nonnegative signed margins and a positive mean margin.
The same check finds no separating direction in any of the four original-sample
Probits. This explains a failure mechanism without changing their optimizer.

The initial diagnostic retained region dummies in both stages. Only 11/50 and
26/100 draws succeeded. In 72 of the 74 unsuccessful draws, at least one
good-region cell contained no zero purchase. Its results are retained in
`initial-region-diagnostic/`, separate from the final specification above.

Independent numerical checks on the final real-data fits give:

| Check | FGNLS | IFGNLS |
|---|---:|---:|
| Maximum scaled estimating-equation derivative error | 8.93e-8 | 1.94e-6 |
| Maximum scaled elasticity-gradient error | 4.72e-10 | 1.51e-9 |
| Maximum relative SE difference using numerical elasticity gradients | 8.51e-10 | 8.07e-10 |
| Maximum relative SE difference using the independently differentiated bread | 1.09e-8 | 5.58e-7 |

The independent bread check uses central differences with Richardson
extrapolation, and complex perturbations for the demand-score Sigma block.
The covariance check retains all Probit, weight and mean contributions.
The first three saved bootstrap vectors also match the package's own worker
bit for bit, using the same whole-cluster sampler and seeds.

Thus, the large IFGNLS disagreement is not reproduced as a derivative
implementation error. IFGNLS has a small residual-covariance eigenvalue
(4.32e-7, versus a largest eigenvalue of 0.0274), and a scaled joint-bread
condition number of 1.22e9. Its smallest scaled final-demand information
eigenvalue is 1.20e-7. These indicate sensitivity and weak local curvature;
they do not prove a complete explanation of the bootstrap discrepancy.
The mean-score diagnostic is small (5.63e-6), and adjusted shares are
0.277, 0.105, 0.416 and 0.194, so neither a large residual score nor a
near-zero elasticity denominator is evident at this point estimate.
This comparison remains an unresolved limit for IFGNLS on this specification.

## Software checks and reproduction

All 77 tests passed, including eight cluster checks. These cover independent
manual summation of the full influence functions, cross-stage terms, unequal
interleaved clusters, chunk and label invariance, singleton-cluster CR0/IID
equivalence, estimation-sample alignment, and whole-cluster worker resampling.
The protected optimizer, model, Probit, parameter, Jacobian and elasticity
files, plus `CITATION.cff` and the top-level README, remain byte-identical to
main `636609f17e732a57b140cbbad5d2bf4042bc396a`; hashes are in `core-audit.json`.

See [Methodology](../../../docs/methodology.md#analytical-inference) for the
single covariance definition. `final/` contains settings, every reported SE,
all successful bootstrap vectors, failure records and the diagnostics.
The input is a private prepared CSV, with the columns shown in the command;
reproducing the price preparation requires the user's supplied U11 pipeline.

```bash
python -m unittest discover -s tests -v
python tools/validate_cluster_se.py --data prepared-four-goods.csv --output comparison \
  --shares share4_g6 share4_g7 share4_g8 share4_g9 \
  --prices p_g6 p_g7 p_g8 p_g9 --expenditure expenditure4 \
  --demographics age scale children --cluster source_psu_str --anot 11
python tools/diagnose_cluster_se.py --data prepared-four-goods.csv --comparison comparison
```
