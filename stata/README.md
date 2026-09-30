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

If you have not installed the Python core yet, install it from your terminal or command prompt:
```bash
pip install pyquaidsce
```
*(Note: If `pyquaidsce` is missing when you run the Stata command, `quaidsce.ado` will automatically attempt to install it for you via pip).*

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

// 3. The direct package command is equivalent
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)

// 4. Estimate with parallel bootstrap standard errors (e.g. 200 replications across 4 CPU cores)
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) mp_context(spawn) rep_timeout(900) seed(12345)
```

`rep_timeout()` combines cooperative checks between numerical work units with a
parent-side process watchdog. When a timeout is set, each replication runs in a
disposable child process so a stuck native call can be terminated without
killing unrelated replications.

Advanced control-function and custom-selection options are available from
Stata. The Shonkwiler-Yen correction always uses the Probit linear index. Precomputed residuals require `reps(0)`:

```stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) control_function(vhat) ///
    selection_control_function(vhat_sel) selection_prices(p3 p1) ///
    selection_covariates(urban) selection_noexpenditure reps(0)
```

For endogenous total expenditure, `ivexp()` is the preferred integrated path.
It estimates log expenditure on log prices, Ray demographics, the excluded
instruments, and a constant; the residual enters both participation Probits and
latent demand equations. Its bootstrap rebuilds the reduced form in every
replication:

```stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) ///
    demographics(hh_size urban) ivexp(log_income employment) ///
    reps(200) seed(12345)
```

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
- **Package missing:** Run `pip install pyquaidsce` in your command line or run `python: import subprocess, sys; subprocess.check_call([sys.executable, "-m", "pip", "install", "pyquaidsce"])` directly inside Stata.
