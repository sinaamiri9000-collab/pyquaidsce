# Stata Package: `pyquaidsce`

Fast Censored Quadratic Almost Ideal Demand System (QUAIDS) estimation in Stata via Python engine.

This package can be called with either `quaidsce` or `pyquaidsce`. Both commands run the same estimator. If you already have the older Stata `quaidsce` installed and want to avoid a name conflict, use `pyquaidsce`.

---

## Prerequisites

1. **Stata 16.0 or newer** (Stata 16, 17, 18, 19 with Python integration).
2. **Python 3.9+** with the `pyquaidsce` package installed.

To verify your Stata Python setup, type in Stata:
```stata
python query
```

If you have not installed the Python core yet, use the Python environment
shown by `python query` to install it from your terminal or command prompt:
```bash
pip install pyquaidsce
```

---

## Installation

### Method 1: Direct Installation from GitHub (Recommended)

Run the following single line inside Stata:

```stata
net install pyquaidsce, from("https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/main/stata") replace
```

Once installed, you can use either command:

```stata
quaidsce ...
pyquaidsce ...
```

Both do the same thing. If you already have the older Stata `quaidsce` installed and want to avoid a name conflict, use `pyquaidsce`.

For help, both names work:

```stata
help quaidsce
help pyquaidsce
```

---

### Method 2: Manual Local Installation

If you cloned the repository locally:

```stata
adopath + "path/to/pyquaidsce/stata"
```
Or copy `quaidsce.ado`, `quaidsce.sthlp`, `pyquaidsce.ado`, and `pyquaidsce.sthlp` directly into your personal Stata PLUS directory (type `sysdir` in Stata to find the exact path, typically `~/ado/plus/p/`).

---

## Files in this Directory

| File | Description |
| :--- | :--- |
| [`pyquaidsce.ado`](pyquaidsce.ado) | The direct implementation of the current Stata interface. Use this command name when you need to distinguish the current package from an older `quaidsce` installation. |
| [`quaidsce.ado`](quaidsce.ado) | Recommended common-name wrapper. It calls the current `pyquaidsce` implementation and keeps the estimator name aligned across Python, R, and Stata. |
| [`quaidsce.sthlp`](quaidsce.sthlp) | Full help for the recommended `quaidsce` command. |
| [`pyquaidsce.sthlp`](pyquaidsce.sthlp) | Full help for the direct `pyquaidsce` command. |
| [`stata.toc`](stata.toc) | Stata package Table of Contents used by Stata's `net` package manager. |
| [`pyquaidsce.pkg`](pyquaidsce.pkg) | Stata package manifest detailing files and metadata for `net install`. |

---

## Quick Start in Stata

```stata
// 1. Load your consumption dataset
use mydata.dta, clear

// 2. Estimate censored QUAIDS with IFGNLS
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)

// 3. or
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)

// 4. Estimate with parallel bootstrap standard errors (e.g. 200 replications across 4 CPU cores)
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) mp_context(spawn) rep_timeout(900) seed(12345)
```

`rep_timeout()` sets the maximum seconds allowed for each bootstrap replication.
If a replication takes longer, it is counted as unsuccessful and the program
moves on to the next one. For example, `rep_timeout(900)` allows 15 minutes.
The default, `rep_timeout(0)`, sets no time limit.


To account for endogenous expenditure, give your instrument variables to
`ivexp()`. For example, use `ivexp(income)` if household income is your chosen
instrument. The package estimates the expenditure equation and adds its
residual to both the participation and demand equations. With bootstrap, it
repeats these steps in every replication:

```stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) ivexp(income) ///
    reps(200) seed(12345)
```

For complete examples with a practice dataset, see
[examples/example_quaidsce.do](../examples/example_quaidsce.do).

---

## Postestimation & Stored Results

Both `quaidsce` and `pyquaidsce` store standard Stata estimation results in `e()`, making them fully compatible with Stata's postestimation toolkit (`test`, `lincom`, `outreg2`, `esttab`):

### Scalars
- `e(N)`: Number of observations in the estimation sample
- `e(ll)`: Log-likelihood
- `e(anot)`: Constant parameter $\alpha_0$ in the translog price index
- `e(ndemo)`: Number of demographic variables
- `e(converged)`: 1 if converged, 0 otherwise
- `e(n_outer)`: Number of outer IFGNLS iterations
- `e(n_gn)`: Number of inner Gauss-Newton steps
- `e(reduced_form_r2)`: R-squared from the internal log-expenditure reduced form
- `e(excluded_iv_F)`, `e(excluded_iv_p)`: Joint excluded-instrument F test and p-value

### Matrices
- `e(b)`: Parameter vector with equation stripes (`alpha`, `beta`, `gamma`, `lambda`, `delta`, `eta`, `rho`, `tau`, `ELAS_INC`, `ELAS_UNCOMP`, `ELAS_COMP`)
- `e(b_est)`: Free estimated structural parameter vector $\theta$ (can be passed to `initial()`)
- `e(V)`: Variance-covariance matrix of estimators (or bootstrap covariance when `reps > 0`)
- `e(Sigma)`: Residual covariance matrix $\Sigma$ (can be passed to `sigma_initial()`)
- `e(elas_i)`: Expenditure (income) elasticities ($1 \times n$)
- `e(elas_u)`: Uncompensated (Marshallian) price elasticities ($n \times n$)
- `e(elas_c)`: Compensated (Hicksian) price elasticities ($n \times n$)
- `e(reduced_form_b)`, `e(reduced_form_V)`: Internal reduced-form coefficients and covariance

### Example: Testing Parameter Restrictions
```stata
// Test equality of income response between goods 1 and 2
test [beta]beta_1 = [beta]beta_2

// View the uncompensated price elasticity matrix
matrix list e(elas_u), format(%10.4f)
```

---

## Troubleshooting

- **Python not found error (`r(7102)`):** Ensure Stata knows where Python is installed. Run `python query` and set your Python path if needed using `set python_exec "C:\path\to\python.exe", permanently`.
- **Package missing:** Run `python query` in Stata, then use that Python executable to run `python -m pip install --upgrade pyquaidsce` from your terminal or command prompt.
