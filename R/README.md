# pyquaidsce: Censored QUAIDS Demand System Estimation in R

An R interface to the high-performance `pyquaidsce` econometric engine, implementing:
- **Quadratic Almost Ideal Demand System (QUAIDS)** (Banks, Blundell, and Lewbel, 1997)
- **Ray (1983) Demographic Scaling** (Poi, 2012)
- **Shonkwiler and Yen (1999) Two-Step Censoring Correction**
- **Iterated Feasible Generalized Nonlinear Least Squares (IFGNLS)**
- **Fast Multiprocessing Nonparametric Bootstrap Standard Errors**
- **Integrated `ivexp` Control Functions for Endogenous Expenditure**
- **External Control Functions & Custom Selection Equations**

---

## Installation

### 1. Requirements
`pyquaidsce` uses `reticulate` to communicate seamlessly with Python. Ensure Python (>= 3.9) and `pyquaidsce` are installed:

```bash
pip install pyquaidsce
```

### 2. Install in R
```r
# Install remotes if not already installed
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install pyquaidsce directly from GitHub
remotes::install_github("sinaamiri9000-collab/pyquaidsce", subdir = "R")
```

---

## Quick Example

```r
library(pyquaidsce)

# Load your household data
# df <- read.csv("household_data.csv")

# 4-good Censored QUAIDS model with demographics
fit <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  anot = 10.0,
  method = "ifgnls",
  censor = TRUE,
  quadratic = TRUE
)

# Standard S3 methods
print(fit)
summary(fit)
coef(fit)
vcov(fit)

# Extract estimated elasticity matrices
fit$elasticities$income          # Expenditure (income) elasticities
fit$elasticities$uncompensated   # Marshallian (uncompensated) price elasticities
fit$elasticities$compensated     # Hicksian (compensated) price elasticities
```

For endogenous total expenditure, pass excluded instruments through `ivexp`:

```r
fit_iv <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  ivexp = c("log_income", "employment_status"),
  anot = 10,
  reps = 200
)

fit_iv$reduced_form$coefficients
fit_iv$reduced_form$excluded.f
```

---

## Convergence controls

| Argument | Default | Criterion |
|---|---|---|
| `param_tol` | `1e-5` | Relative inner parameter change |
| `objective_tol` | `1e-7` | Relative inner weighted-SSR change |
| `gn_tol` | `1e-5` | Scaled Gauss-Newton criterion |
| `outer_param_tol` | `1e-5` | Relative IFGNLS parameter change in two consecutive rounds |

Point estimation and every bootstrap draw use these same four thresholds.
Any one inner criterion may stop the inner solve; IFGNLS requires two
consecutive small outer parameter changes and final inner convergence.
The iteration limits default to `max_iter=300` and `max_outer=200`.
See the [user guide](../docs/user-guide.md#convergence-and-numerical-tolerance)
for the exact formulas. In R, pass the same
argument names, for example `param_tol = 1e-5`.

## Citation & Author
- **Author:** Sina Amiri (Department of Economics, Shiraz University)
- **Email:** sinaamiri9000@gmail.com
- **Repository:** [https://github.com/sinaamiri9000-collab/pyquaidsce](https://github.com/sinaamiri9000-collab/pyquaidsce)
