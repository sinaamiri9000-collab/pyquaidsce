# rquaidsce: Censored QUAIDS Demand System Estimation in R

An R interface to the high-performance `pyquaidsce` econometric engine, implementing:
- **Quadratic Almost Ideal Demand System (QUAIDS)** (Banks, Blundell, and Lewbel, 1997)
- **Ray (1983) Demographic Scaling** (Poi, 2012)
- **Shonkwiler and Yen (1999) Two-Step Censoring Correction**
- **Iterated Feasible Generalized Nonlinear Least Squares (IFGNLS)**
- **Fast Multiprocessing Nonparametric Bootstrap Standard Errors**
- **Control Functions for Endogeneity & Custom Selection Equations**

---

## Installation

### 1. Requirements
`rquaidsce` uses `reticulate` to communicate seamlessly with Python. Ensure Python (>= 3.9) and `pyquaidsce` are installed:

```bash
pip install pyquaidsce
```

### 2. Install in R
```r
# Install devtools if needed
install.packages("devtools")

# Install directly from GitHub
devtools::install_github("sinaamiri9000-collab/pyquaidsce", subdir = "rquaidsce")
```

---

## Quick Example

```r
library(rquaidsce)

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

---

## Citation & Author
- **Author:** Sina Amiri (Department of Economics, Shiraz University)
- **Email:** sinaamiri9000@gmail.com
- **Repository:** [https://github.com/sinaamiri9000-collab/pyquaidsce](https://github.com/sinaamiri9000-collab/pyquaidsce)
