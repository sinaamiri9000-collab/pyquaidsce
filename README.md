# pyquaidsce

`pyquaidsce` is a Python package for estimating **Censored Quadratic Almost Ideal Demand Systems (QUAIDS)**, incorporating:

- **Ray (1983) demographic scaling** for household characteristics.
- **Shonkwiler & Yen (1999) two-step correction** for zero expenditure shares (censoring).
- **Nonlinear SUR estimation** supporting NLS, FGNLS, and Iterated FGNLS (IFGNLS).
- **Elasticities**: Expenditure (income), Marshallian (uncompensated), and Hicksian (compensated) price elasticities.
- **Control functions**: Supports integrated expenditure endogeneity correction via `ivexp`, plus externally generated residuals.
- **Parallel bootstrap**: Multi-core bootstrap with progress tracking and runtime safeguards for valid standard errors.
- **​Usable across both Stata and R**: You can call the same Python estimator directly within Stata and R.
---

## Motivation

This project originated from a very practical challenge in my own empirical research. While using the excellent Stata `quaidsce` command developed by Dr. Juan Carlos Caro and colleagues to estimate censored demand systems for my research, the long computation times—especially under Iterated FGNLS (IFGNLS)—became a serious bottleneck for model diagnostics, specification testing, and sensitivity analysis. Furthermore, the long computational time was not due to how quaidsce was implemented, but rather because Stata's `nlsur` command—frequently used for censored QUAIDS estimation—is a general-purpose package that is not optimized for demand systems, especially ones as computationally demanding as censored QUAIDS

At the same time, Python offered a much more flexible and convenient environment for building automated data pipelines and scaling computational workloads on virtual servers. These two factors motivated me to develop `pyquaidsce`: a native Python implementation designed to preserve the rigorous econometric structure of the original Stata command while making censored QUAIDS estimation fast, practical, and easy to integrate into modern empirical workflows.

---

## Why pyquaidsce?

Estimating censored demand systems in empirical research often requires extensive model exploration, sensitivity checks, and bootstrap replications. In Stata, estimating large censored QUAIDS models under Iterated FGNLS (IFGNLS) can take significant computing time (often 20+ minutes for a single run on moderately sized datasets). 

`pyquaidsce` was developed to solve this practical bottleneck:
- **Fast Execution**: Uses analytic Jacobians and NumPy/SciPy linear algebra. In the four-good benchmark below, pyquaidsce was about **88.6 times faster** with 200 bootstrap replications.
- **Pure Python & Lightweight**: Uses `numpy`, `scipy`, `pandas`, and `threadpoolctl`. No complex compilation or heavy external dependencies.
- **Research Workflow Integration**: Easily run demand models in Jupyter notebooks, script automated sensitivity pipelines, and run on cloud servers/clusters.
- **Validated Numerics**: Tested against econometric identities and numerical derivatives; the repository also preserves high-precision historical Stata comparison evidence from earlier compatibility-capable releases.

---

## Installation

Install `pyquaidsce` directly from PyPI via `pip`:

```bash
pip install pyquaidsce
```

For local development or installing from source:

```bash
git clone https://github.com/sinaamiri9000-collab/pyquaidsce.git
cd pyquaidsce
pip install -e .
```

**Requirements**: Python >= 3.9, with `numpy`, `scipy`, `pandas`, and `threadpoolctl`.

---

## Quick Start

Here is a minimal example estimating a 4-good censored QUAIDS model:

```python
import pandas as pd
from pyquaidsce import quaidsce

# Load your household data
df = pd.read_csv("household_data.csv")

# Estimate the model
res = quaidsce(
    df,
    shares=["w1", "w2", "w3", "w4"],          # Budget shares (must contain zeros if censored)
    prices=["p1", "p2", "p3", "p4"],          # Prices corresponding to each good
    expenditure="total_expenditure",           # Total expenditure across the system
    demographics=["hh_size", "urban"],         # Demographic scaling variables (Ray 1983)
    anot=10.0,                                 # Price index constant (alpha_0)
    method="ifgnls",                           # default; alternatives: 'nls' or 'fgnls'
    reps=0,                                    # Set reps=200+ for bootstrap standard errors
    verbose=True,
)

# View summary and elasticity tables
print(res.summary())
print(res.elasticity_tables())
```

To treat total expenditure as endogenous, supply one or more excluded
instruments. The package regresses log expenditure on log prices,
demographics, the excluded instruments, and a constant; its residual enters
both participation and demand equations:

```python
res_iv = quaidsce(
    df,
    shares=["w1", "w2", "w3", "w4"],
    prices=["p1", "p2", "p3", "p4"],
    expenditure="total_expenditure",
    demographics=["hh_size", "urban"],
    ivexp=["log_income", "employment_status"],
    anot=10.0,
    reps=200,  # rebuilds the reduced form in every bootstrap draw
)
print(res_iv.reduced_form_table())
```

### Parallel bootstrap in Python

If you use parallel bootstrap (`n_jobs > 1`) in a Python script, put your estimation code inside `main()` and add this at the end of the file:

```python
if __name__ == "__main__":
    main()
```

This is required by Python multiprocessing when it uses `spawn`.

For a step-by-step tutorial, see [Getting Started](docs/getting-started.md). For a complete reference of all input arguments, control functions, and output attributes, see the [User Guide & API Reference](docs/user-guide.md).

---

## Using pyquaidsce in Stata

you can use `pyquaidsce` in stata. `pyquaidsce` includes a Stata interface (`quaidsce.ado`) for fitting censored QUAIDS models with the Python engine:

```stata
// 1. Install the Stata package directly from GitHub
net install pyquaidsce, from("https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/main/stata") replace

// 2. Estimate your model in Stata with familiar syntax
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) ivexp(log_income employment_status) anot(10) method(ifgnls)

// 3. Postestimation commands work seamlessly
test [beta]beta_1 = [beta]beta_2
matrix list e(elas_u)
```

In Stata, you can use either command name:

```stata
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban)

pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban)
```

Both commands run the same estimator. If you already have the older Stata `quaidsce` installed and want to avoid a name conflict, use `pyquaidsce`.

For help, both names work:

```stata
help quaidsce
help pyquaidsce
```

See the [Stata Package Guide](stata/README.md) for full details, options, and troubleshooting.

---

## Using pyquaidsce in R

Prefer working in R? `pyquaidsce` provides a native R interface with standard S3 methods (`summary`, `coef`, `vcov`, `print`):

```r
# 1. Install remotes (if not already installed) & install the R interface from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("sinaamiri9000-collab/pyquaidsce", subdir = "R")

# 2. Estimate Censored QUAIDS in R
library(pyquaidsce)
fit <- quaidsce(
  data = df,
  shares = c("w1", "w2", "w3", "w4"),
  prices = c("p1", "p2", "p3", "p4"),
  expenditure = "total_exp",
  demographics = c("hh_size", "urban"),
  ivexp = c("log_income", "employment_status"),
  anot = 10.0,
  method = "ifgnls"
)

summary(fit)
fit$elasticities$income
fit$elasticities$uncompensated
```

See the [R Package Documentation](R/README.md) for full details.

---

## Benchmark

We compared **pyquaidsce** with Stata's **quaidsce** on the same computer.
Both programs fit the same censored QUAIDS model using IFGNLS, with 2,665
households, 4 goods, 2 demographic variables, and 200 bootstrap replications
using 8 workers.

| Measure | pyquaidsce | quaidsce_c |
|---|---:|---:|
| Time | 122.05 seconds | 10,810.18 seconds |
| Bootstrap replications | 200 | 200 |

pyquaidsce was about **88.6 times faster** in this run.

The code, data, results, and file descriptions are in [benchmark/pyquaidsce vs quaidsce/](benchmark/pyquaidsce%20vs%20quaidsce/).

---

## Repository Structure

```text
src/pyquaidsce/   Core Python package source code and Stata bridge
stata/            Official Stata package (quaidsce command; pyquaidsce distribution)
R/                R interface to the Python estimator
tests/            Mathematical unit tests, theory checks, and Stata regression tests
examples/         Python, Stata, and R examples with a shared practice dataset
benchmark/        Comparisons of model results and estimation time
docs/             Methodology, getting started, user guide, and validation
tools/            Diagnostic scripts and Stata log comparison utilities
```

---

## Documentation

- [Getting Started Guide](docs/getting-started.md): Practical tutorial, data prep, and Stata-to-Python option map.
- [User Guide & API Reference](docs/user-guide.md): Complete reference for all input parameters, the result object, and code examples.
- [Stata Package Guide](stata/README.md): How to install and run `pyquaidsce` directly within Stata.
- [Methodology & Model Equations](docs/methodology.md): QUAIDS model specification, Shonkwiler-Yen censoring, and elasticity derivations.
- [Examples](examples/README.md): Complete Python, Stata, and R scripts using the same practice dataset.
- [Validation Evidence](docs/validation.md): Automated tests and the comparison with Stata.
- [Performance & Benchmarking](docs/performance.md): Benchmark methodology, timing details, and optimization notes.
- [Contributing](CONTRIBUTING.md): Guidelines for bug reports and contributions.

---

## Citation

If you use `pyquaidsce` in your research, please cite the software:

Amiri, S., Molina Chueca, J. A., Salgado, J. C., & Caro, J. C. (2026). pyquaidsce: Estimation of censored AIDS and QUAIDS demand systems in Python, Stata, and R (Version 1.8.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.23463166

BibTeX:

```bibtex
@software{amiri2026pyquaidsce,
  author = {Amiri, Sina and Molina Chueca, José Alberto and Salgado, Juan Carlos and Caro, Juan Carlos},
  title = {pyquaidsce: Estimation of censored {AIDS} and {QUAIDS} demand systems in {Python}, {Stata}, and {R}},
  year = {2026},
  version = {1.8.0},
  publisher = {Zenodo},
  doi = {10.5281/zenodo.23463166},
  url = {https://github.com/sinaamiri9000-collab/pyquaidsce/releases/tag/v1.8.0}
}
```

The DOI for all versions is [10.5281/zenodo.22990015](https://doi.org/10.5281/zenodo.22990015).

---

## Authors

- **Sina Amiri:** Department of Economics, Shiraz University, Shiraz, Iran.
- **José Alberto Molina Chueca:** Department of Economic Analysis, University of Zaragoza, Zaragoza, Spain; IZA – Institute of Labor Economics, Bonn, Germany.
- **Juan Carlos Salgado:** Instituto Nacional de Salud Pública, Cuernavaca, Mexico; Ministry of Science, Humanities, Technology, and Innovation, Mexico.
- **Juan Carlos Caro:** Departamento de Economía, Universidad de Concepción, Concepción, Chile.

Maintainer: Sina Amiri ([sinaamiri9000@gmail.com](mailto:sinaamiri9000@gmail.com)).

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0-only)**. See [`LICENSE`](LICENSE) for details.
