# Validation of `pyquaidsce` against Stata `quaidsce` v2.0

Two independent Stata benchmarks were used, plus a set of checks that need no
Stata at all.

The formulas and the small benchmark were independently re-audited on
2026-08-02. See `docs/development/audit-v1.0.1-fa.md` for the defects found in the original Python
delivery and the repairs in version 1.0.1.

---

## Benchmark A — `bench/small4.log`, full double precision

4 goods, 2000 observations, 2 demographics, censoring in all four shares
(708 / 60 / 187 / 190 zeros). Built by `bench/benchmark_small.do` from
the repository's own `data/DS_STATA_3_2_0_pci2sls_.dta`. Every coefficient was
dumped from Mata with `%21.17g`, so the comparison is at **17 significant
digits**, not Stata's 7-digit display format.

```
quaidsce_c sw1 sw2 sw4 sw9, anot(10) reps(1) prices(p1 p2 p4 p9) ///
    expenditure(total) demographics(x1 x2) nolog [method(...)|noquadratic|nocensor]
```

### A1. Estimated end to end from `theta = 0`, as Stata does

| run | specification | rel. diff in `e(ll)` | max abs diff, 104 coefficients | max abs diff, standard errors |
|---|---|---|---|---|
| 3 | `method(ifgnls)` | **8.5e-09** | **8.4e-07** | **8.6e-08** |
| 1 | default = `method(fgnls)` | **4.1e-07** | **2.1e-05** | **9.0e-06** |
| 4 | `noquadratic` | **1.1e-07** | **6.2e-06** | n/a¹ |
| 2 | `method(nls)` | 5.7e-03 | 4.8e+00 | 3.5e-01 |

¹ the do-file only dumped point estimates for that run.

All four cover the whole reported vector: `alpha`, `beta`, `gamma`, `lambda`,
`delta`, `eta`, `rho`, the 32 first-stage `tau`, and the `ELAS_INC` /
`ELAS_UNCOMP` / `ELAS_COMP` blocks.

**Run 2 is not a failure of the port — Stata stopped before convergence.**
Evaluating *our* NLS criterion at *Stata's* `method(nls)` parameter vector:

```
at Stata's theta        objective = 230.213306755748   relative gradient = 6.9e-04
at our optimum          objective = 230.055832984043   relative gradient = 5.7e-13
```

Our objective is **lower**, and Stata's point has a relative gradient four
orders of magnitude away from zero. Stata's `nl`/`nlsur` declare convergence
when **any** of `tolerance()` (1e-5, on the coefficient vector), `ltolerance()`
(1e-7, on the objective) or `nrtolerance()` (1e-5, on the scaled gradient) is
met. Because the criterion of a censored demand system has very flat
directions, the coefficient test fires while the gradient is still O(1e-4).
`method(nls)` therefore returns a path-dependent, not-fully-converged estimate,
and `method(fgnls)` inherits part of that through `Sigma_hat`.
`method(ifgnls)` is a genuine fixed point and is the only one of the three that
is reproducible across implementations — which is exactly what the table shows.

Passing `stop_rule="stata"` reproduces Stata's disjunctive stopping rule, but it
does **not** improve the agreement (it made run 1 worse, 6.2e-02), because an
early stop depends on the whole iteration path, including Stata's numerical
derivatives. **Recommendation: use `method="ifgnls"`.**

### A2. Evaluated at Stata's own parameter vector — no optimiser involved

Holding `theta` fixed at Stata's 17-digit values isolates the formulas from the
optimiser. Maximum absolute deviations:

| run | parameter restrictions² | rel. `e(ll)` | `ELAS_INC` | `ELAS_UNCOMP` | `ELAS_COMP` |
|---|---|---|---|---|---|
| 1 `fgnls` | **0.0** | 1.2e-09 | **2.1e-10** | **9.7e-11** | **1.0e-10** |
| 2 `nls` | **0.0** | 6.3e-10 | 1.1e-08 | 7.6e-08 | 8.0e-08 |
| 3 `ifgnls` | 1.1e-19 | 8.5e-09 | **1.9e-10** | **1.3e-10** | **1.3e-10** |
| 4 `noquadratic` | **0.0** | 9.8e-09 | **2.2e-10** | **1.0e-10** | **1.1e-10** |

² the *derived* coefficients — `gamma_i,n`, `gamma_n,n` and `eta_r,n`, which
Stata obtains from symmetry, homogeneity and adding-up — recomputed from our
restrictions and compared with Stata's reported values. Exact to the last bit.

First-stage probits, identical across runs:

```
max |tau_py     - tau_stata|      = 1.5e-09
max |se(tau)_py - se(tau)_stata|  = 2.6e-10
```

`e(V)` for `fgnls`: Stata uses the `Sigma_hat` of the **first** step, not of the
final residuals. Substituting each candidate at Stata's `fgnls` point:

| `Sigma_hat` used in `e(V)` | max abs diff in se |
|---|---|
| from the first-step (NLS) residuals — **what Stata does** | **2.1e-03**³ |
| from the final (FGNLS) residuals | 2.8e-01 |
| identity | 3.2e+00 |

³ residual gap only because Stata's *unconverged* first-step `theta` is itself
known only from the log. `pyquaidsce`'s default `vce_sigma="objective"`
implements the correct convention, which is why run 1 above reaches 9.0e-06.

---

## Benchmark B — `log/censor_v2.log`, the author's own 17-good run

```
quaidsce w1-w17, anot(10) prices(p1-p17) expenditure(total) nolog demographics(x1-x3)
(obs = 15,147) ... 18 FGNLS iterations
Censored Quadratic AIDS model, Alpha_0 = 10, Log-likelihood = 482472.11
649 reported coefficients
```

Two things had to be reverse-engineered before this run could be reproduced,
because the log records the command but not the data preparation:

1. **`total`** is not a variable of the shipped data. It is
   `total_exp * (w1 + ... + w17)`, which the file already carries as
   `total_exp_DS`; its 15147 strictly positive values are exactly the reported
   sample size (the full file has 15151 rows).
2. **The shares were rescaled to sum to one.** The reported `alpha`s span
   -0.19 to +0.32, which is impossible for the raw `w1..w17` (means 0.0004 to
   0.037, summing to 0.169). With raw shares our `e(ll)` comes out at 754786
   (56 % off); with `w_i / sum_j w_j` it comes out at **482471.83** against
   Stata's **482472.11**. The run is a *conditional* (within-group) demand
   system.

### B1. Estimated end to end from `theta = 0`

`method(ifgnls)`, Gauss-Newton with step halving, 42 outer iterations / 417
Gauss-Newton steps, 30 minutes on 2 cores:

| quantity | python | stata | deviation |
|---|---|---|---|
| `e(N)` | 15147 | 15147 | — |
| `e(ll)` | 482472.0957 | 482472.11 | **1.4e-02 absolute, 3.0e-08 relative** |
| all 649 coefficients | | | **max 5.5e-05 absolute** |
| all standard errors | | | **max 1.2e-06 absolute** |

The largest coefficient deviations are on `rho_x1` (0.204708 vs 0.204762) and
`alpha_2` (-0.556503 vs -0.556554), i.e. the fourth decimal of coefficients of
order 0.2. They are optimiser-path differences, not formula differences: Stata's
outer loop stopped after 18 iterations, ours ran to 42 to reach a tighter fixed
point, and the criterion is very flat in those directions. Section B2 removes
the optimiser from the comparison entirely and the same quantities then agree to
3.6e-07.

### B2. Evaluated at Stata's reported parameter vector

| quantity | deviation |
|---|---|
| parameter restrictions (`gamma_.,17`, `gamma_17,17`, `eta_.,17`) | 3.6e-07 |
| `e(ll)`: 482471.83 vs 482472.11 | 5.8e-07 relative |
| relative gradient at Stata's point | 2.0e-06 |
| all 275 delta-method standard errors | max 3.6e-07, median rel. 3.5e-05 |
| all 374 first-stage `tau` | max 4.3e-06 |
| all 374 `se(tau)` | max 4.5e-07 |

Every one of these sits at the floor set by the log's 7-significant-digit
display format, i.e. this is a match to display precision.

The same table with `first_stage_predict="xb"` (the textbook Shonkwiler–Yen
transformation) gives `e(ll)` = 387660, a 20 % discrepancy — independent proof
that the shipped command really does feed `predict`'s default *probability*
into `normal()`/`normalden()`.

### Conditioning

The share normalisation is not only a matter of reproducing the log; it
transforms the numerics. On the raw shares the Gauss-Newton step had to be
halved 16 times and after 30 iterations the relative gradient was still 2.4e-02.
On the normalised shares the NLS stage converged in **32 iterations to a
relative gradient of 6e-13 taking full, undamped steps**.

---

## Checks that need no Stata

`python -m unittest discover -s tests -v` (which collects the internal checks
and a Stata regression test):

| check | result |
|---|---|
| delta matrix vs numerical Jacobian of the free → full map, all `(n, R, quadratic, censor)` | 1.4e-10 |
| block decomposition of the delta matrix rebuilds it | exact |
| symmetry, homogeneity, `eta` adding-up; `alpha`/`beta`/`lambda` adding-up under `nocensor` | exact |
| analytic model Jacobian vs central differences | 9.7e-11 relative |
| fast free-space Jacobian vs `Jfull @ Delta` | 3.0e-16 |
| probit vs an independent BFGS maximum-likelihood fit | 1.4e-08 in `b`, 1.1e-13 in the log-likelihood |

`tests/test_theory.py`, on the `nocensor` path (which Stata v2.0 cannot run at
all — see below), with the elasticities evaluated at model-consistent shares:

| identity | residual |
|---|---|
| adding up: `sum alpha - 1`, `sum beta`, `sum lambda` | exactly 0 |
| symmetry `max abs(gamma - gamma')`, homogeneity `max abs(row sums)` | exactly 0 |
| Engel aggregation `sum_i w_i e_i = 1` | 1.1e-16 |
| Cournot/homogeneity `sum_j eu_ij = -e_i` | 2.2e-16 |
| Hicksian homogeneity `sum_j ec_ij = 0` | 2.2e-16 |

---

## Defects found in `quaidsce` v2.0

Ordered by how much they can change a published result.

### 1. `noquadratic` + `demographics()` + censoring silently returns wrong expenditure and compensated elasticities

`quaidsce_c.ado` assigns the latent expenditure elasticity to a **global** macro
and then reads a **local** of the same name:

```stata
global ie`i' = 1+`betanz`i''/`w_`i''m          // <- global
if "`quadratic'" == "" {                       // false under -noquadratic-
    local ie`i' = ...
}
...
local ie`i' = (1+1/we`i'*((cdf`i'm*((`ie`i''-1)*`w_`i''m))+ ...   // <- local, empty
```

With `noquadratic` the `local` is never set, so `(`ie`i''-1)` expands to the
literal `(-1)`: the latent expenditure elasticity is taken to be **zero**.
Stata does not error, it just returns a wrong number, and `ELAS_COMP` inherits
it through the Slutsky identity. On benchmark A:

| good | Stata reports | correct value |
|---|---|---|
| sw1 | 0.398150 | 1.057176 |
| sw2 | 0.289537 | 0.595117 |
| sw4 | 0.285952 | 1.380100 |
| sw9 | 0.242633 | 1.209569 |

Stata's numbers make every good a necessity; correctly, two of the four are
luxuries. `pyquaidsce` reproduces this exactly with `strict_stata=True`
(deviation 2.2e-10, see A2) and gives the correct values with
`strict_stata=False`.

### 2. The first stage uses `Phi(Phi(x'tau))` instead of `Phi(x'tau)`

```stata
quietly predict du`i'                    // after -probit-: default is Pr(y=1)
qui replace pdf`i' = normalden(du`i')
qui replace cdf`i' = normal(du`i')
```

Stata's `predict` after `probit` defaults to the predicted probability, not the
linear predictor. Shonkwiler & Yen (1999) require `Phi(x'tau)` and `phi(x'tau)`.
Because `Phi(x'tau)` lies in `(0,1)`, the shipped code compresses `cdf` into
`(0.500, 0.841)` and `pdf` into `(0.242, 0.399)`: the censoring correction is
attenuated towards a constant, and `delta` absorbs the difference.
`first_stage_predict="pr"` reproduces it; `"xb"` is the textbook estimator.
On benchmark B the two differ by 20 % in `e(ll)`.

### 3. `nocensor` always errors

`np_prob` is only defined inside the censoring branch, but the elasticity
section references it unconditionally:

```stata
local loc = `np_prob'*(`i'-1)+`neqn'+1
```

With `nocensor` the macro is empty, so Stata sees `local loc = *(1-1)+4+1` and
raises `unknown function *()`, `r(133)`. Confirmed on benchmark A, run 5. The
whole `nocensor` path of v2.0 is therefore dead; `pyquaidsce` implements it
(it is Poi's `quaids`) and it satisfies every demand-theory identity exactly.

### 4. `lnexpenditure()` silently drops log expenditure from the first stage

The macro holding the log-expenditure temporary variable is only set in the
`expenditure()` branch, but the first-stage regressor list uses it:

```stata
local zvar `lnprices' `lnexp' `demographics'
```

So `quaidsce ..., lnexpenditure(lnm)` estimates the participation probits
without log expenditure, and `np_prob` is off by one, which shifts every
positional index into `tau` used by the elasticity corrections. Reproduced (and
reported in `res.notes`) for faithfulness.

### 5. `ELAS_UNCOMP` / `ELAS_COMP` are labelled transposed

The row vector is filled with `i` in the outer loop, the name stripe is built
with `j` in the outer loop, so `e(b)`'s `e_a_b` holds `eu_{b,a}`. The diagonal
is unaffected, which is why the package's own `example.do` (which reads only
`e_1_1`, `e_2_2`, …) never trips over it. Anything reading off-diagonal
elasticities out of `e(b)` — e.g. via `parmest` — gets the transpose.

### 6. With no demographics, the uncompensated elasticity uses `beta_i` where Poi (2012) has `beta_j`

```stata
-(`beta'[1,`i']*`lambda'[1,`i']/exp(`bofp')*(`lnexp'm-`lnpindex')^2)
```

The demographics branch correctly uses `betanz_j`. Affects only
`demographics()`-free specifications. `strict_stata=False` uses the published
formula.

### 7. `method(nls)` is not converged

See A1. It is reported as an estimate, but it is an intermediate iterate.

---

## Reproducing all of it

```bash
cd pyquaidsce
python3 -m pip install -e .
python3 -m unittest discover -s tests -v
python3 tests/test_internals.py
python3 tests/test_theory.py
python3 tools/validate_small4.py --runs 1 2 3 4
python3 tools/validate_small4_atstata.py
python3 tools/diagnose_nls.py
python3 tools/check_at_stata_theta.py        # benchmark B
python3 tools/check_firststage.py            # benchmark B, stage 1
python3 tools/validate_censor_v2.py --method ifgnls --algorithm gn
```

The large-data tools auto-discover a sibling `quaidsce-master/` directory. If
the repository lives elsewhere, set `QUAIDSCE_STATA_REPO=/absolute/path/to/quaidsce-master`.
