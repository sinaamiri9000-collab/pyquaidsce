# Stata Package: `pyquaidsce`

Fast Censored Quadratic Almost Ideal Demand System (QUAIDS) estimation in Stata via Python engine.

This package provides the official Stata command `pyquaidsce`, allowing Stata users to estimate censored QUAIDS demand systems with the **up to 44.6x speedup** of the `pyquaidsce` Python core while staying completely within the native Stata workflow.

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
*(Note: If `pyquaidsce` is missing when you run the Stata command, `pyquaidsce.ado` will automatically attempt to install it for you via pip).*

---

## Installation

### Method 1: Direct Installation from GitHub (Recommended)

Run the following single line inside Stata:

```stata
net install pyquaidsce, from("https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/main/stata") replace
```

Once installed, Stata will recognize `pyquaidsce` globally and make the interactive help file available via:
```stata
help pyquaidsce
```

---

### Method 2: Manual Local Installation

If you cloned the repository locally:

```stata
adopath + "path/to/pyquaidsce/stata"
```
Or copy `pyquaidsce.ado` and `pyquaidsce.sthlp` directly into your personal Stata PLUS directory (type `sysdir` in Stata to find the exact path, typically `~/ado/plus/p/`).

---

## Files in this Directory

| File | Description |
| :--- | :--- |
| [`pyquaidsce.ado`](pyquaidsce.ado) | The official Stata command program. Parses user syntax, handles sample filtering (`if`/`in`), invokes the Python bridge, posts matrices and scalars to `e()`, and formats the regression and elasticity tables. |
| [`pyquaidsce.sthlp`](pyquaidsce.sthlp) | The official interactive Stata Help file rendered inside Stata's Viewer when running `help pyquaidsce`. |
| [`stata.toc`](stata.toc) | Stata package Table of Contents used by Stata's `net` package manager. |
| [`pyquaidsce.pkg`](pyquaidsce.pkg) | Stata package manifest detailing files and metadata for `net install`. |

---

## Quick Start in Stata

```stata
// 1. Load your consumption dataset
use mydata.dta, clear

// 2. Estimate censored QUAIDS with IFGNLS
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) method(ifgnls)

// 3. Estimate with parallel bootstrap standard errors (e.g. 200 replications across 4 CPU cores)
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) anot(10) reps(200) n_jobs(4) seed(12345)
```

---

## Postestimation & Stored Results

`pyquaidsce` stores standard Stata estimation results in `e()`, making them fully compatible with Stata's postestimation toolkit (`test`, `lincom`, `outreg2`, `esttab`):

### Scalars
- `e(N)`: Number of observations in the estimation sample
- `e(ll)`: Log-likelihood
- `e(anot)`: Constant parameter $\alpha_0$ in the translog price index
- `e(ndemo)`: Number of demographic variables
- `e(converged)`: 1 if converged, 0 otherwise
- `e(n_outer)`: Number of outer IFGNLS iterations
- `e(n_gn)`: Number of inner Gauss-Newton steps

### Matrices
- `e(b)`: Parameter vector with equation stripes (`alpha`, `beta`, `gamma`, `lambda`, `delta`, `eta`, `rho`, `tau`, `ELAS_INC`, `ELAS_UNCOMP`, `ELAS_COMP`)
- `e(V)`: Variance-covariance matrix of estimators (or bootstrap covariance when `reps > 0`)
- `e(elas_i)`: Expenditure (income) elasticities ($1 \times n$)
- `e(elas_u)`: Uncompensated (Marshallian) price elasticities ($n \times n$)
- `e(elas_c)`: Compensated (Hicksian) price elasticities ($n \times n$)

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
