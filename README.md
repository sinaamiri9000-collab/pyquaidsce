# pyquaidsce

`pyquaidsce` is a Python package for estimating **Censored Quadratic Almost Ideal Demand Systems (QUAIDS)**. It provides a fast, native Python alternative to Stata's official `quaidsce` command, incorporating:

- **Ray (1983) demographic scaling** for household characteristics.
- **Shonkwiler & Yen (1999) two-step correction** for zero expenditure shares (censoring).
- **Nonlinear SUR estimation** supporting NLS, FGNLS, and Iterated FGNLS (IFGNLS).
- **Elasticities**: Expenditure (income), Marshallian (uncompensated), and Hicksian (compensated) price elasticities.
- **Parallel bootstrap** for valid standard errors and confidence intervals.
- **Direct Stata compatibility**: matched coefficient ordering, output tables, and exact numerical reproduction switches.

---

## Motivation

This project originated from a very practical challenge in my own empirical research. While using the excellent Stata `quaidsce` command developed by Dr. Juan Carlos Caro and colleagues to estimate censored demand systems for my research, the long computation times—especially under Iterated FGNLS (IFGNLS)—became a serious bottleneck for model diagnostics, specification testing, and sensitivity analysis. 

At the same time, Python offered a much more flexible and convenient environment for building automated data pipelines and scaling computational workloads on virtual servers. These two factors motivated me to develop `pyquaidsce`: a native Python implementation designed to preserve the rigorous econometric structure of the original Stata command while making censored QUAIDS estimation fast, practical, and easy to integrate into modern empirical workflows.

---

## Why pyquaidsce?

Estimating censored demand systems in empirical research often requires extensive model exploration, sensitivity checks, and bootstrap replications. In Stata, estimating large censored QUAIDS models under Iterated FGNLS (IFGNLS) can take significant computing time (often 20+ minutes for a single run on moderately sized datasets). 

`pyquaidsce` was developed to solve this practical bottleneck:
- **Fast Execution**: Written with optimized analytic Jacobians and vectorized linear algebra (`numpy`/`scipy`), achieving a **~44x speedup** over Stata on standard benchmarks.
- **Pure Python & Lightweight**: Only requires `numpy`, `scipy`, and `pandas`. No complex compilation or heavy external dependencies.
- **Research Workflow Integration**: Easily run demand models in Jupyter notebooks, script automated sensitivity pipelines, and run on cloud servers/clusters.
- **Verified Accuracy**: Delivers parameter estimates and elasticities that match Stata benchmarks up to high numerical precision (within `1e-5` to `1e-7`).

---

## Installation

Clone this repository and install it locally using `pip`:

```bash
git clone https://github.com/sinaamiri9000-collab/pyquaidsce.git
cd pyquaidsce
pip install .
```

For development or editable use:

```bash
pip install -e .
```

**Requirements**: Python >= 3.9, with standard `numpy`, `scipy`, and `pandas`.

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
    method="ifgnls",                           # 'nls', 'fgnls', or 'ifgnls'
    first_stage_predict="pr",                  # 'pr' matches Stata; 'xb' uses linear index
    strict_stata=True,                         # Replicates Stata's exact formulas
    reps=0,                                    # Set reps=200+ for bootstrap standard errors
    verbose=True,
)

# View summary and elasticity tables
print(res.summary())
print(res.elasticity_tables())
```

For a detailed step-by-step tutorial and Stata-to-Python parameter translation, see [Getting Started](docs/getting-started.md).

---

## Key Options & Stata Compatibility

`pyquaidsce` provides two switches to let you choose between literal Stata replication and textbook formulas:

1. **`first_stage_predict`**:
   - `"pr"` *(default)*: Replicates Stata `quaidsce`'s default behavior, which uses the predicted probability $\Phi(X'\tau)$ inside the normal density/distribution terms.
   - `"xb"`: Uses the linear index $X'\tau$, following the standard textbook Shonkwiler–Yen formulation.

2. **`strict_stata`**:
   - `True` *(default)*: Keeps documented Stata conventions and index ordering for 1-to-1 replication against Stata `.log` files.
   - `False`: Uses corrected textbook formulas for documented edge-cases (such as models without demographics or with `noquadratic`).

See [Stata Compatibility](docs/stata-compatibility.md) for full details.

---

## Validation & Performance Benchmark

We evaluated `pyquaidsce` against Stata 19.5 on a deterministic synthetic benchmark dataset with **20,000 observations, 4 goods, 3 demographic variables, censoring across all goods, and IFGNLS estimation**:

- **Numerical Agreement**:
  - Structural parameters ($\alpha, \beta, \gamma, \lambda, \delta, \eta, \rho$): Maximum difference $< 1.68 \times 10^{-5}$
  - First-stage Probit parameters ($\tau$): Maximum difference $< 1.21 \times 10^{-7}$
  - Expenditure and price elasticities: Maximum difference $< 7.34 \times 10^{-7}$
  - Log-likelihood: Relative difference $< 4.99 \times 10^{-8}$

- **Speed Comparison (Same-Machine Wall Clock)**:
  - **Stata 19.5**: 1,161.2 seconds (~19 minutes, 21 seconds)
  - **pyquaidsce**: 26.0 seconds
  - **Speedup**: **~44.6x faster**

All raw data, logs, scripts, and comparison tables are available in [`benchmarks/cquaids_ifgnls_4g_20k/`](benchmarks/cquaids_ifgnls_4g_20k/). For more details, see [Validation](docs/validation.md) and [Performance](docs/performance.md).

---

## Repository Structure

```text
src/pyquaidsce/   Core package source code
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
- [Methodology & Model Equations](docs/methodology.md): QUAIDS model specification, Shonkwiler-Yen censoring, and elasticity derivations.
- [Stata Compatibility Guide](docs/stata-compatibility.md): Explanations of `first_stage_predict`, `strict_stata`, and replication tips.
- [Validation Evidence](docs/validation.md): Numerical comparison across 113 parameters and elasticities against Stata.
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
