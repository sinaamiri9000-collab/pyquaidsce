# Getting Started with pyquaidsce

This guide walks you through preparing your data, specifying your demand system, running your first estimation, and interpreting the output.

---

## 1. Installation

Install `pyquaidsce` in your Python environment:

```bash
pip install pyquaidsce
```

Or install locally from source:

```bash
git clone https://github.com/sinaamiri9000-collab/pyquaidsce.git
cd pyquaidsce
pip install -e .
```

---

## 2. Preparing Your Data

To estimate a demand system with `quaidsce(...)`, your DataFrame should contain:

1. **Budget shares** ($w_1, \dots, w_n$): The proportion of total expenditure allocated to each good (e.g., $w_i = p_i q_i / m$). For a censored model, these must contain zeros for households that did not consume good $i$.
2. **Prices** ($p_1, \dots, p_n$) or **Log Prices** ($\ln p_1, \dots, \ln p_n$): Must be provided in the exact same order as the shares.
3. **Total Expenditure** ($m$) or **Log Total Expenditure** ($\ln m$): Total spending across all $n$ goods in the system.
4. **Demographic variables** ($z_1, \dots, z_R$): Household characteristics (such as household size, urban dummy, education level) used for Ray (1983) demographic scaling.

> [!TIP]
> **Comparing with Stata?**  
> Make sure your estimation sample in Python matches your Stata sample exactly. Check that the number of observations, variable definitions, and any sample filters (such as dropping missing values or non-positive expenditures) are identical.

---

## 3. Basic Estimation Example

Here is a complete script estimating a 4-good censored QUAIDS model:

```python
import pandas as pd
from pyquaidsce import quaidsce

# 1. Load data
df = pd.read_csv("household_consumption.csv")

# 2. Define variable lists
shares = ["w_meat", "w_dairy", "w_cereal", "w_other"]
prices = ["p_meat", "p_dairy", "p_cereal", "p_other"]
demographics = ["hh_size", "urban"]

# 3. Estimate model
res = quaidsce(
    data=df,
    shares=shares,
    prices=prices,
    expenditure="total_exp",
    demographics=demographics,
    anot=10.0,                  # Price index constant (alpha_0)
    method="ifgnls",            # Iterated FGNLS
    algorithm="gn",             # Gauss-Newton optimizer
    first_stage_predict="xb",   # Textbook Shonkwiler-Yen linear index (default)
    strict_stata=False,         # Corrected textbook formulas (True replicates Stata exactly)
    reps=0,                     # Set reps > 0 to run bootstrap
    verbose=True,
)

# 4. Display results
print(res.summary())
print(res.elasticity_tables())
```

---

## 4. Understanding the Estimation Methods

- **`method="nls"`**: Nonlinear Least Squares (assumes identity error covariance, $\Sigma = I$).
- **`method="fgnls"`**: Feasible Generalized NLS (two-step estimation using the residual covariance from NLS).
- **`method="ifgnls"`** *(Recommended)*: Iterated FGNLS. Continuously updates the residual covariance matrix $\Sigma$ and parameter estimates until convergence. This corresponds to Maximum Likelihood under joint normality.

*Note:* When running `method="ifgnls"`, you will see messages for NLS and initial FGNLS iterations. These are standard starting steps of the IFGNLS algorithm.

---

## 5. Stata-to-Python Option Translation

| Stata Command Option | Python Argument in `quaidsce(...)` |
|---|---|
| `w1 w2 w3 w4` (share varlist) | `shares=["w1", "w2", "w3", "w4"]` |
| `prices(p1 p2 p3 p4)` | `prices=["p1", "p2", "p3", "p4"]` |
| `lnprices(lp1 lp2 lp3 lp4)` | `lnprices=["lp1", "lp2", "lp3", "lp4"]` |
| `expenditure(total)` | `expenditure="total"` |
| `lnexpenditure(ltotal)` | `lnexpenditure="ltotal"` |
| `demographics(z1 z2)` | `demographics=["z1", "z2"]` |
| `anot(10)` | `anot=10.0` |
| `method(ifgnls)` | `method="ifgnls"` |
| `first_stage_predict(xb)` | `first_stage_predict="xb"` |
| `first_stage_predict(pr)` | `first_stage_predict="pr"` |
| `noquadratic` | `quadratic=False` |
| `nocensor` | `censor=False` |
| `initial(b_init)` | `initial=b_init` |
| `sigma_initial(sigma_init)` | `sigma_initial=sigma_init` |
| `reps(200)` | `reps=200` |

---

## 6. Computing Standard Errors via Bootstrap

Because the Shonkwiler–Yen method is a two-step estimator, the conventional standard errors for the second stage do not account for the first-stage Probit estimation error. To obtain valid standard errors and confidence intervals for all parameters and elasticities, use the built-in bootstrap:

```python
if __name__ == "__main__":
    res = quaidsce(
        df,
        shares=shares,
        prices=prices,
        expenditure="total_exp",
        demographics=demographics,
        anot=10.0,
        method="ifgnls",
        reps=200,          # 200 bootstrap replications
        n_jobs=4,          # Parallel execution across 4 CPU cores
        seed=123456,       # Reproducible random seed
        mp_context="spawn", # BLAS-safe cross-platform worker start
        rep_timeout=900,   # Cooperative + hard 15-minute limit per replication
    )

    # Summary table will now include bootstrap standard errors
    print(res.summary())
```

> [!NOTE]
> When using multiprocessing (`n_jobs > 1`) on Windows, always enclose your script within `if __name__ == "__main__":`.

After a successful bootstrap, `res.V` and `res.se` are the bootstrap covariance
and standard errors. The conditional analytical reference is retained in
`res.V_analytic` and `res.analytic_se`. Without bootstrap, analytical elasticity
S.E.s are deliberately `NaN`; the rest of `res.V` remains finite.

`rep_timeout` is checked cooperatively at Probit, nonlinear-iteration, and
chunk boundaries. A parent-side watchdog also terminates the disposable child
process if a native call remains stuck. External control-function models require
an outer bootstrap that rebuilds the reduced form and cannot use the internal
`reps` option. Models using `ivexp` can use internal `reps`, because the package
knows and re-estimates that reduced form in every replication.
