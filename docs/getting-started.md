# Getting started

This guide assumes you are comfortable identifying variables in your dataset but
does not assume advanced Python knowledge.

## 1. Create an environment and install the package

From the repository root:

```bash
python -m venv .venv
```

On Windows:

```powershell
.venv\Scripts\python.exe -m pip install --upgrade pip
.venv\Scripts\python.exe -m pip install -e .
```

On macOS/Linux:

```bash
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -e .
```

## 2. Prepare the data

A call to `quaidsce(...)` needs:

1. all budget-share variables;
2. all price variables, in exactly the same order as the shares;
3. total system expenditure (or its logarithm);
4. Ray demographic-scaling variables when the censored specification is used;
5. the translog price-index constant `anot`.

When comparing Python with an existing Stata run, reproduce the **same estimation
sample** first. Any observations dropped in Stata before estimation must also be
dropped in Python. A different sample can generate economically meaningful
coefficient and elasticity differences even if the two estimators are identical.

For a conditional subsystem, check whether the included budget shares are meant
to sum to one. Rescaling shares without a corresponding economic justification
changes the object being estimated; see `performance.md` for the numerical
conditioning issue that motivated the benchmark normalization check.

## 3. A complete estimation file

Create `estimate.py`:

```python
import pandas as pd
from pyquaidsce import quaidsce


def main():
    df = pd.read_csv("mydata.csv")

    shares = ["w1", "w2", "w3", "w4"]
    prices = ["p1", "p2", "p3", "p4"]
    demographics = ["household_size", "age_head"]

    res = quaidsce(
        df,
        shares=shares,
        prices=prices,
        expenditure="total_expenditure",
        demographics=demographics,
        anot=10.0,
        method="ifgnls",
        algorithm="gn",
        start="zero",
        first_stage_predict="pr",
        strict_stata=True,
        reps=0,
        verbose=True,
    )

    print("Converged:", res.converged)
    print(res.summary())
    print(res.elasticity_tables())


if __name__ == "__main__":
    main()
```

Run it with:

```bash
python estimate.py
```

The `if __name__ == "__main__":` guard is especially important when bootstrap
replications are run in parallel on Windows.

## 4. Stata-to-Python option map

| Stata | Python |
|---|---|
| expenditure-share varlist | `shares=[...]` |
| `prices(...)` | `prices=[...]` |
| `lnprices(...)` | `lnprices=[...]` |
| `expenditure(x)` | `expenditure="x"` |
| `lnexpenditure(x)` | `lnexpenditure="x"` |
| `demographics(...)` | `demographics=[...]` |
| `anot(#)` | `anot=#` |
| `method(ifgnls)` | `method="ifgnls"` |
| `noquadratic` | `quadratic=False` |
| `nocensor` | `censor=False` |
| `initial(...)` | `initial=...` |
| `reps(#)` | `reps=#` |

Give only one of `prices`/`lnprices` and only one of
`expenditure`/`lnexpenditure`.

## 5. NLS, FGNLS, and IFGNLS

`method="ifgnls"` necessarily begins with an NLS estimate, uses its residuals to
estimate the cross-equation covariance matrix, performs FGNLS, and then repeats
that covariance/parameter update until the outer fixed point converges. Messages
such as `Calculating NLS estimates` and `FGNLS iteration ...` are therefore
stages of IFGNLS, not extra unrelated models.

## 6. Starting values and bootstrap

The default `start="zero"` mirrors the Stata starting convention. A successful
result can be reused as an explicit restart:

```python
res2 = quaidsce(
    df,
    shares=shares,
    prices=prices,
    expenditure="total_expenditure",
    demographics=demographics,
    anot=10.0,
    method="ifgnls",
    initial=res.theta,
    sigma_initial=res.sigma,
    reps=0,
)
```

For bootstrap inference, first establish that the main model converges with
`reps=0`. Then try a small diagnostic number of replications before a final
large run. `bootstrap_start="zero"` reproduces the zero-start logic; `"warm"`
reuses the full-sample solution and is a performance/fidelity tradeoff.

Example:

```python
res = quaidsce(
    df,
    shares=shares,
    prices=prices,
    expenditure="total_expenditure",
    demographics=demographics,
    anot=10.0,
    method="ifgnls",
    reps=200,
    bootstrap_start="zero",
    seed=123456,
    n_jobs=4,
)
```

Always inspect the number and nature of failed bootstrap replications rather than
assuming every resample converged.
