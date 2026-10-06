# Cluster analytical SE: completed real-data checks

Further experiments were stopped at the user's request. This report records
completed checks and keeps unresolved findings explicit. Work is on ANSE;
main and the protected estimation files are unchanged.

## Main findings

1. Joint cluster covariance and whole-cluster bootstrap are implemented. The
   complete influence vector includes Probits, demand, estimated GLS weights
   and sample means. Independent checks support the implementation.
2. Removing children from the participation Probits, while retaining them in
   demand, removes the observed Probit failures: all 100 revised draws succeed.
3. On groups 6–9, revised FGNLS has median analytical / bootstrap elasticity SE
   0.976. Revised IFGNLS remains at 1.690, with material differences by target.
   This discrepancy has not been fixed and cannot be called a validation pass.
4. A separate groups-1–4 diagnostic has IFGNLS median 0.965, but its range is
   0.518–1.142. FGNLS has median 0.802 there. Agreement is model-dependent;
   neither median establishes agreement for every elasticity or interval coverage.
5. Tighter stopping tolerances do not remove the large-parameter root in a
   selected bootstrap sample. Existing warm starts reach another, better root
   in that sample, but do not resolve the overall IFGNLS elasticity-SE gap.

## Data and shared settings

The source is the user's supplied `cldenew.csv`, not the later frozen
`cldenewfinal.csv` described in the article brief. Its SHA256 is
`9819184dc68b471247253c74c1f257b98dfc983d3cb6f45218390ea0b49e9fb1`.
Household rows and identifiers are not published here.

All models use urban households in seasons 3 and 4, quadratic censored demand,
`anot=11`, GN, and zero-start point estimates. There are no region dummies in
the final specifications, curvature constraints, IVEXP or external control
functions. Demand demographics are age, household size and children.
Prices come from the supplied full-sample U11 preparation, are not renormalized
after subsetting, and remain fixed in both inference methods. Each bootstrap
resamples their household rows; it rebuilds Probits, demand, covariance and means.
Neither method includes uncertainty from rebuilding the external price pipeline.

The cluster is identified by the nine-digit household-address prefix. This
matches the supplied article bootstrap grouping, without claiming an official
PSU definition. Probits and demand remain unweighted. Cluster covariance alone
does not implement survey weights, stratification or finite-population correction.

Expenditure and shares are recalculated within each four-good basket. For groups
6–9, 9,747 urban/season rows become 9,744 after excluding three zero-basket
households, in 1,792 clusters. The reference groups 1–4 use 9,282 households in
1,789 clusters. Group prices remain those prepared on the full source sample.

The four package default tolerances remain unchanged. Iteration limits are
400/400, chunks 2,000, one BLAS thread per fit and two bootstrap workers.
Seed is 20261006. The 100-draw run retains the first 50 requested draws and adds
50 more. Failures and unusual estimates are retained in the records; failures
are not replaced and extreme successful draws are not trimmed.

## Completed comparisons

Each ratio is analytical cluster SE / bootstrap SE over all 36 elasticities.
Cluster covariance uses G/(G-1). Draw counts are successful / requested.

| Specification | Method | Draws | Median ratio | Ratio range | Bootstrap time (s) |
|---|---|---:|---:|---:|---:|
| Initial groups 6–9; children in both stages; zero bootstrap | FGNLS | 43/50 | 1.112 | 0.799–1.341 | 19.7 |
| Initial groups 6–9; children in both stages; zero bootstrap | FGNLS | 87/100 | 1.084 | 0.863–1.386 | 39.7 |
| Initial groups 6–9; children in both stages; zero bootstrap | IFGNLS | 43/50 | 1.671 | 0.792–2.344 | 50.9 |
| Initial groups 6–9; children in both stages; zero bootstrap | IFGNLS | 87/100 | 1.677 | 0.801–2.375 | 102.5 |
| Same groups 6–9; warm bootstrap diagnostic | IFGNLS | 43/50 | 1.674 | 0.862–2.329 | 24.8 |
| Same groups 6–9; warm bootstrap diagnostic | IFGNLS | 87/100 | 1.683 | 0.891–2.329 | 49.9 |
| Revised groups 6–9; children only in demand; warm bootstrap | FGNLS | 50/50 | 1.030 | 0.718–1.296 | 34.9 |
| Revised groups 6–9; children only in demand; warm bootstrap | FGNLS | 100/100 | 0.976 | 0.796–1.104 | 64.1 |
| Revised groups 6–9; children only in demand; warm bootstrap | IFGNLS | 50/50 | 1.633 | 0.951–2.542 | 27.4 |
| Revised groups 6–9; children only in demand; warm bootstrap | IFGNLS | 100/100 | 1.690 | 0.922–2.409 | 54.0 |
| Reference groups 1–4; children only in demand; warm bootstrap | FGNLS | 50/50 | 0.778 | 0.476–1.310 | 16.8 |
| Reference groups 1–4; children only in demand; warm bootstrap | FGNLS | 100/100 | 0.802 | 0.437–1.166 | 32.4 |
| Reference groups 1–4; children only in demand; warm bootstrap | IFGNLS | 50/50 | 0.998 | 0.551–1.365 | 17.4 |
| Reference groups 1–4; children only in demand; warm bootstrap | IFGNLS | 100/100 | 0.965 | 0.518–1.142 | 37.2 |

All detailed timings are saved in the corresponding summary JSONs. These are
machine-specific times, not general performance claims. The first diagnostic
also included region dummies in both stages and succeeded in only 11/50 and
26/100 draws; its separate folder preserves those results.

## Requested statistic percentages

The statistic is estimate / SE, using null value zero. The percentage ratio is
bootstrap statistic / analytical statistic, times 100. With the same estimate
and null, it equals 100 times analytical SE / bootstrap SE. Thus 169% means a
statistic 69% larger; 92% means one 8% smaller. This is the existing large-sample
Wald statistic, not a new Student-t degrees-of-freedom correction.

For the latest revised groups-6–9 model, using all 100 successful draws:

| Expenditure elasticity | Analytical statistic | Bootstrap statistic | Percentage ratio |
|---|---:|---:|---:|
| Group 6 | 9.988 | 15.113 | 151.3% |
| Group 7 | 1.733 | 4.088 | 235.9% |
| Group 8 | 38.778 | 35.761 | 92.2% |
| Group 9 | 2.910 | 6.233 | 214.2% |

Across all 36 elasticities, the median is 169.0% and range 92.2–240.9%.
Every coefficient and elasticity is retained in each `*-t-statistics.csv`,
with both statistics, the percentage ratio and percentage change. Zero/zero
cases remain missing. The initial four-good model's coefficient medians differ
from its elasticity median: demand coefficients 84.4%, Probits 92.2%, elasticities
167.7%. A single elasticity summary should not stand in for all coefficients.

## Investigated causes

### Probit failures: a demonstrated support problem

The original groups-6–9 sample has only 38, 28, 351 and 51 zeros respectively.
Only two of the 28 non-buying group-7 households have children, in two distinct
clusters. Every one of the 13 failed initial draws loses both such households.
Then every sampled household with children buys group 7, permitting a separating
Probit direction. Linear programming independently confirms separation in all
13 failures; it finds none in the four original-sample Probits.

Removing children only from the Probits was approved by the user. The revised
100-draw test has no Probit failures. Demand retains the children variable.
This is a specification change for the validation model, not a Probit-optimizer
or package-default change.

### Large demand coefficients: different local roots

On initial bootstrap draw 3, the following runs use the same observations and
first-stage model. Stricter settings are an explicit diagnostic using the existing
API, not changed defaults.

| Start / tolerance | rho for household size | Log-likelihood | Outer / inner |
|---|---:|---:|---:|
| Zero / defaults | 462.063 | 76596.265 | 24 / 132 |
| Zero / stricter | 462.063 | 76596.265 | 31 / 194 |
| Warm / defaults | 0.0217 | 80084.691 | 19 / 27 |
| Warm / stricter | 0.0217 | 80084.691 | 25 / 52 |

Stricter thresholds do not remove the large-parameter root. Warm starts find
a higher-likelihood root in this draw. They do not certify a global optimum,
remove every parameter tail, or fix the overall elasticity-SE discrepancy.
The warm-start 100-draw median ratio remains 1.683.

### IFGNLS sensitivity and the limits of the local approximation

On revised groups 6–9, the smallest residual-covariance eigenvalue is 3.78e-7
versus a largest eigenvalue of 0.0274. The scaled joint-bread condition number
is 1.46e9. Adjusted shares are about 0.277, 0.105, 0.416 and 0.194; a near-zero
elasticity denominator is not evident at the point estimate.

For group-7 expenditure elasticity, the single largest cluster contributes
21.5% of estimated variance, the largest five contribute 51.1%, and the largest
ten contribute 72.6%. Across all elasticities, the median largest-ten share is
54.6%. A nominal total of 1,792 clusters therefore does not imply that every
target's information is spread across many clusters.

A separate leave-one-cluster diagnostic refits 20 samples with warm starts and
stricter existing tolerances. It includes eight clusters with the largest predicted
group-7 elasticity changes and twelve seed-721 random clusters. All fits converge.
The influence prediction and actual change have correlation 0.991 for that
elasticity, but relative RMS error is 20.0% (20.2% across all targets). Some
random-cluster changes are close, while influential-cluster changes are less
linear. This supports finite-sample sensitivity; it is not a full explanation
or a correction of the 100-draw discrepancy.

The reference groups 1–4 were chosen from groups with materially more zero
purchases, without replacing the user's groups 6–9. Their IFGNLS point fit has
condition number 4.92e6 and smallest residual-covariance eigenvalue 2.41e-4.
Its 100-draw median ratio is closer to one, but individual differences remain.
This is an exploratory model comparison, not proof of a general cure. A further
point-only groups-2–5 pilot used 9,230 observations and had a bread condition
number about 1.62e7; no bootstrap comparison was completed for that basket.

## Implementation checks and what remains open

Independent real-data checks on revised groups 6–9 give:

| Check | FGNLS | IFGNLS |
|---|---:|---:|
| Maximum scaled estimating-equation derivative error | 9.71e-8 | 1.08e-6 |
| Maximum scaled elasticity-gradient error | 3.62e-10 | 2.85e-9 |
| Maximum relative SE error using numerical elasticity gradients | 1.38e-9 | 8.09e-10 |
| Maximum relative SE error using independently differentiated bread | 9.69e-9 | 6.94e-7 |

The independent bread check uses central differences with Richardson
extrapolation and complex Sigma perturbations. The first three revised warm
bootstrap vectors match the production worker bit for bit. These checks do not
show a derivative or cluster-summation bug that explains the large IFGNLS gap.

All 77 software tests passed, including eight cluster tests. The last published
implementation check also passed all ten Ubuntu/Windows and Python-3.9–3.13
CI jobs. Point estimates, residual covariance, Probit coefficients, likelihood
and iteration counts are unchanged when analytical clustering is toggled.
The six protected core files, `CITATION.cff` and top-level README are byte-identical
to main `636609f17e732a57b140cbbad5d2bf4042bc396a`; see `core-audit.json`.

No SEs are rescaled to force agreement with bootstrap. No draws are silently
removed, no failed Probits are treated as successful, and no optimizer, stopping
threshold default, covariance-update sequence or model core is changed.
The revised groups-6–9 IFGNLS discrepancy remains open. Bootstrap precision
with 50/100 draws is limited, and no confidence-interval coverage study has
been performed for these real-data models.

## Files and reproduction

| Folder | Contents |
|---|---|
| `initial-region-diagnostic/` | First groups-6–9 check with region dummies |
| `initial-four-good/` | Groups 6–9, no regions, children in both stages, zero bootstrap; root probe |
| `warm-start-diagnostic/` | Same initial specification with warm bootstrap |
| `revised-four-good/` | User-approved children-only-in-demand model; diagnostics and influence probes |
| `reference-four-good/` | Exploratory groups-1–4 comparison |

Each comparison folder saves settings, complete SE tables, replication records,
all successful coefficient vectors and failure lists. Reported numbers above
can be traced to those files. Household input CSVs remain private.
See [Methodology](../../../docs/methodology.md#analytical-inference) for the
single mathematical definition and assumptions.

For the revised four-good model, with the private prepared CSV:

```bash
python tools/validate_cluster_se.py --data prepared-four-goods.csv --output comparison \
  --shares share4_g6 share4_g7 share4_g8 share4_g9 \
  --prices p_g6 p_g7 p_g8 p_g9 --expenditure expenditure4 \
  --demographics age scale children --selection-covariates age scale \
  --cluster source_psu_str --anot 11 --bootstrap-start warm
python tools/diagnose_cluster_se.py --data prepared-four-goods.csv --comparison comparison
```

The archived root and influence diagnostics are in
`tools/probe_cluster_roots.py` and `tools/probe_cluster_influence.py`. Both accept
`--data` and `--comparison`; the root probe defaults to initial draw 3. Use the
initial-four-good settings for the root probe and revised-four-good settings for
the influence probe. No additional fits were run while preparing this report.

Reproducing the original price construction requires the supplied U11 pipeline.
The public tools receive the prepared prices and do not recreate that pipeline.
