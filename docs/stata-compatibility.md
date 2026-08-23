# Stata Compatibility and Implementation Guide

`pyquaidsce` is designed to provide exact numerical reproduction of the Stata `quaidsce` command (v2.0, Caro et al. 2025), while also offering options to apply textbook econometric formulas where the original Stata ado-file contains documented quirks.

This guide explains the compatibility settings and how to achieve 1-to-1 replication.

---

## 1. First-Stage Probit Predictions: `first_stage_predict`

In the Shonkwiler & Yen (1999) two-step procedure, the participation probability $\Phi_i$ and normal density $\phi_i$ are functions of the Probit index $X_i'\tau_i$:
- **Textbook Shonkwiler–Yen**: $\Phi_i = \Phi(X_i'\tau_i)$ and $\phi_i = \phi(X_i'\tau_i)$

In Stata, the `predict` command immediately following `probit` generates the **predicted probability** by default rather than the linear index ($X_i'\tau_i$). The original `quaidsce_c.ado` file calculates `normal(predict)` and `normalden(predict)`, which computes $\Phi(\Phi(X'\tau))$ and $\phi(\Phi(X'\tau))$.

### Options:
- **`first_stage_predict="xb"`** *(Default)*: Replicates Stata's exact implementation bit-for-bit. Use this if you are comparing results directly against Stata.
- **`first_stage_predict="xb"`**: Uses the linear index $X_i'\tau_i$ inside $\Phi(\cdot)$ and $\phi(\cdot)$, matching standard textbook theory.

---

## 2. Replicating Stata Elasticities: `strict_stata`

During comprehensive code validation against Stata `quaidsce`, two specific behaviors were identified in the Stata ado-file's elasticity calculations:

1. **Quadratic model without demographics (`ndemo=0`)**:
   In the uncompensated price elasticity formula, Stata's ado-file multiplies the quadratic term by $\beta_i$ instead of $\beta_j$.
2. **Linear AIDS with demographics and censoring (`quadratic=False` + `ndemo>0` + `censor=True`)**:
   Stata assigns the latent income elasticity to a global macro and then inadvertently calls an empty local macro, effectively treating the latent elasticity as 0 in the censoring adjustment.

### Options:
- **`strict_stata=False`** *(Default)*: Applies the corrected theoretical formulas (Poi 2012 / Shonkwiler & Yen 1999).
- **`strict_stata=True`**: Replicates Stata's exact returned values so that automated tests and diffs against Stata log files match.

---

## 3. Checklist for Exact Replication against Stata

If you are trying to reproduce an existing Stata estimation in Python and see discrepancies, verify the following steps in order:

1. **Verify the Estimation Sample**:
   Ensure that observations dropped due to missing values (`markout` in Stata) or non-positive values are identical in both programs.
2. **Variable Ordering**:
   Verify that your `prices` and `shares` lists have the exact same ordering of goods.
3. **Total Expenditure Definition**:
   Verify that `expenditure` is the sum of spending across the goods in the system (or log expenditure if `lnexpenditure` is used).
4. **Optimization Starting Values**:
   By default, both Stata `nlsur` and `pyquaidsce` use `start="zero"`. If you provided custom starting values in Stata (`initial(...)`), provide the same vector to `initial=...` in Python.
5. **Estimation Method**:
   Confirm whether you are using `method="ifgnls"`, `method="fgnls"`, or `method="nls"`.

---

## Summary of Defaults

| Parameter | Default | Effect |
|---|---|---|
| `first_stage_predict` | `"pr"` | Matches Stata default probability calculation |
| `strict_stata` | `False` | Corrected textbook elasticity calculations (use `True` to match the ado exactly) |
| `start` | `"zero"` | Starts optimization from zeros, matching Stata `nlsur` |
| `method` | `"fgnls"` | Feasible Generalized NLS (use `"ifgnls"` for iterated) |
