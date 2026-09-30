# pyquaidsce

`pyquaidsce` is a Python package for estimating **Censored Quadratic Almost Ideal Demand Systems (QUAIDS)**. It provides a fast, native Python alternative to the original Stata `quaidsce` command developed by Dr. Juan Carlos Caro and colleagues, incorporating:

- **Ray (1983) demographic scaling** for household characteristics.
- **Shonkwiler & Yen (1999) two-step correction** for zero expenditure shares (censoring).
- **Nonlinear SUR estimation** supporting NLS, FGNLS, and Iterated FGNLS (IFGNLS).
- **Elasticities**: Expenditure (income), Marshallian (uncompensated), and Hicksian (compensated) price elasticities.
- **Control functions**: Supports integrated expenditure endogeneity correction via `ivexp`, plus externally generated residuals.
- **Parallel bootstrap**: Multi-core bootstrap with progress tracking and runtime safeguards for valid standard errors.
- **Stata integration**: matched coefficient ordering and native Stata-facing output tables, while using the corrected textbook econometric formulas.

---

## Motivation

This project originated from a very practical challenge in my own empirical research. While using the excellent Stata `quaidsce` command developed by Dr. Juan Carlos Caro and colleagues to estimate censored demand systems for my research, the long computation times—especially under Iterated FGNLS (IFGNLS)—became a serious bottleneck for model diagnostics, specification testing, and sensitivity analysis. 

At the same time, Python offered a much more flexible and convenient environment for building automated data pipelines and scaling computational workloads on virtual servers. These two factors motivated me to develop `pyquaidsce`: a native Python implementation designed to preserve the rigorous econometric structure of the original Stata command while making censored QUAIDS estimation fast, practical, and easy to integrate into modern empirical workflows.

---

## Why pyquaidsce?

Estimating censored demand systems in empirical research often requires extensive model exploration, sensitivity checks, and bootstrap replications. In Stata, estimating large censored QUAIDS models under Iterated FGNLS (IFGNLS) can take significant computing time (often 20+ minutes for a single run on moderately sized datasets). 

`pyquaidsce` was developed to solve this practical bottleneck:
- **Fast Execution**: Written with optimized analytic Jacobians and vectorized linear algebra (`numpy`/`scipy`), achieving a **44.6x wall-clock speedup** on the controlled benchmark reported below.
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

Prefer working in Stata? `pyquaidsce` includes an official Stata package (`quaidsce.ado`) that lets you estimate censored QUAIDS models directly inside Stata while harnessing Python's **up to 44.6x speedup**:

```stata
// 1. Install the Stata package directly from GitHub
net install pyquaidsce, from("https://raw.githubusercontent.com/sinaamiri9000-collab/pyquaidsce/main/stata") replace

// 2. Estimate your model in Stata with familiar syntax
quaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban) ivexp(log_income employment_status) anot(10) method(ifgnls)

// 3. Postestimation commands work seamlessly
test [beta]beta_1 = [beta]beta_2
matrix list e(elas_u)
```

The recommended command is `quaidsce`, matching the estimator name used in Python and R. The equivalent command `pyquaidsce` runs the same current implementation. If you still have an older installation of the original Stata `quaidsce` command and want to distinguish the current package explicitly, use:

```stata
pyquaidsce w1 w2 w3 w4, prices(p1 p2 p3 p4) expenditure(total_exp) demographics(hh_size urban)
```

Both help names are available:

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
remotes::install_github("sinaamiri9000-collab/pyquaidsce", subdir = "rquaidsce")

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

See the [R Package Documentation](rquaidsce/README.md) for full details.

---

## Canonical Econometric Behavior

Starting with **v1.6.0**, `pyquaidsce` has one canonical implementation for the two areas where the original Stata `quaidsce` ado contains known deviations from the published formulas:

- The Shonkwiler–Yen first stage always uses the Probit **linear index** $X'\tau$, so the correction uses $\Phi(X'\tau)$ and $\phi(X'\tau)$.
- Elasticities always use the corrected theoretical formulas, including the no-demographics quadratic term and the demographics + linear-AIDS censoring branch.

See [Stata Compatibility](docs/stata-compatibility.md) for details about comparisons with the original Stata command.

---

## Validation & Historical Performance Benchmark

The repository preserves a deterministic 20,000-observation benchmark against Stata 19.5 (4 goods, 3 demographics, censoring, IFGNLS). The stored comparison was produced during an earlier compatibility-capable release and is retained as historical validation and timing evidence. Because v1.6.0 intentionally removes the legacy error-replication paths, exact agreement with the original ado is no longer a release goal in specifications affected by those known deviations.

- **Numerical Agreement**:
  - Structural parameters ($\alpha, \beta, \gamma, \lambda, \delta, \eta, \rho$): Maximum difference $< 1.68 \times 10^{-5}$
  - First-stage Probit parameters ($\tau$): Maximum difference $< 1.21 \times 10^{-7}$
  - Expenditure and price elasticities: Maximum difference $< 7.34 \times 10^{-7}$
  - Log-likelihood: Relative difference $< 4.99 \times 10^{-8}$

- **Speed Comparison (Same-Machine Wall Clock)**:
  - **Stata 19.5**: 1,161.2 seconds (~19 minutes, 21 seconds)
  - **pyquaidsce**: 26.0 seconds
  - **Speedup**: **~44.6x faster**

All raw data, archived logs, scripts, and comparison tables are available in [`benchmarks/cquaids_ifgnls_4g_20k/`](benchmarks/cquaids_ifgnls_4g_20k/). For current v1.6.0 validation policy and the historical timing evidence, see [Validation](docs/validation.md) and [Performance](docs/performance.md).

---

## Repository Structure

```text
src/pyquaidsce/   Core Python package source code and Stata bridge
stata/            Official Stata package (quaidsce command; pyquaidsce distribution)
tests/            Mathematical unit tests, theory checks, and Stata regression tests
examples/         Ready-to-run Python and Stata sample scripts
benchmarks/       Reproducible benchmark data, scripts, and logs
docs/             Methodology, getting started, Stata compatibility, and validation guides
tools/            Diagnostic scripts and Stata log comparison utilities
```

---

## Documentation

- [Getting Started Guide](docs/getting-started.md): Practical tutorial, data prep, and Stata-to-Python option map.
- [User Guide & API Reference](docs/user-guide.md): Complete reference for all input parameters, the result object, and code examples.
- [Stata Package Guide](stata/README.md): How to install and run `pyquaidsce` directly within Stata.
- [Methodology & Model Equations](docs/methodology.md): QUAIDS model specification, Shonkwiler-Yen censoring, and elasticity derivations.
- [Stata Compatibility Guide](docs/stata-compatibility.md): Numerical-comparison guidance and the documented differences from the original Stata ado.
- [Validation Evidence](docs/validation.md): Current v1.6.0 test policy plus archived Stata comparison evidence across 113 reported values.
- [Performance & Benchmarking](docs/performance.md): Benchmark methodology, timing details, and optimization notes.
- [Contributing](CONTRIBUTING.md): Guidelines for bug reports and contributions.

---

## Citation

If you use `pyquaidsce` in your research, please cite both this package and the original Stata `quaidsce` command:

```bibtex
@software{amiri2026pyquaidsce,
  author = {Amiri, Sina},
  title = {pyquaidsce: Fast Censored QUAIDS Demand System Estimation in Python},
  year = {2026},
  url = {https://github.com/sinaamiri9000-collab/pyquaidsce}
}

@techreport{caro2021quaidsce,
  author = {Caro, Juan Carlos and Melo, Grace and Molina, J. A. and Salgado, J. C.},
  title = {Censored QUAIDS estimation with quaidsce},
  institution = {Boston College Department of Economics},
  type = {Boston College Working Papers in Economics},
  number = {1045},
  year = {2021},
  url = {https://ideas.repec.org/p/boc/bocoec/1045.html}
}
```

---

## Author

**Sina Amiri**  
Department of Economics, Shiraz University, Shiraz, Iran.

---

## License

This project is licensed under the **GNU General Public License v3.0 (GPL-3.0-only)**. See [`LICENSE`](LICENSE) for details.
