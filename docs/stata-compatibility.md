# Stata compatibility and implementation differences

Version 1.0.1 was designed to make direct comparisons with Stata `quaidsce`
possible while also documenting cases where literal reproduction and the
published formula are not the same object.

## First-stage prediction: `pr` versus `xb`

The shipped Stata implementation obtains `predict` after each probit without
requesting the linear predictor. Its subsequent normal CDF/PDF calculations
therefore operate on the predicted probability. In this package:

- `first_stage_predict="pr"` reproduces that shipped Stata path;
- `first_stage_predict="xb"` uses the probit linear index inside the
  Shonkwiler–Yen CDF/PDF terms.

Use `"pr"` when reproducing an existing Stata estimate. Treat `"xb"` as a
methodological choice that should be reported explicitly rather than as a
formatting option.

## `strict_stata`

`strict_stata=True` retains documented Stata-specific elasticity behavior where
needed for direct comparison. `strict_stata=False` uses the corresponding
corrected formula documented during the validation work.

The most important rule is simple: **do not compare two programs while changing
this switch at the same time**. First establish numerical reproduction under the
same convention; only then run an intentional sensitivity analysis.

## Other documented upstream behaviors

The validation work identified several paths in the Stata v2.0 implementation
that require caution, including the interaction of censoring, demographics and
`noquadratic`; the uncensored path; the first-stage use of `lnexpenditure()`;
and the stored orientation/naming of elasticity vectors. The complete evidence,
including benchmark output and code-path notes, is retained in
[`validation.md`](validation.md).

The Python implementation also contains explicit repairs where a literal port
would make a supported Python interface unusable—for example, a valid censored
restart from an already fitted `res.theta`.

## Sample equivalence comes before estimator equivalence

A recurring source of apparently different Python/Stata results is a different
estimation sample. When reproducing a Stata result:

1. apply exactly the same filters before estimation;
2. verify the final observation count;
3. verify the same share and price order;
4. verify the same definition of total system expenditure;
5. then compare optimization and compatibility settings.

Small sample differences can materially alter coefficients and elasticities in a
nonlinear demand system even when both implementations are correct.
