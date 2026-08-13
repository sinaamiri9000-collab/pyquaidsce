# pyquaidsce

`pyquaidsce` is a Python implementation of the censored quadratic almost-ideal
demand system (QUAIDS) implemented by the Stata command `quaidsce`. It combines
Ray demographic scaling, the Shonkwiler–Yen two-step correction for zero budget
shares, nonlinear SUR estimation, expenditure and price elasticities, and an
optional bootstrap.

Version **1.0.1** is a compatibility-focused release: the estimator has been
checked against Stata benchmark runs, while documented switches make it
possible to distinguish literal Stata reproduction from selected textbook
formulas.

> **Research-software status.** This package is intended for empirical research.
> Always document your sample construction, model specification, convergence,
> and the compatibility switches used in a published analysis.

## Motivation

This project started from a very practical problem in my own research. While
using the Stata `quaidsce` command developed by Dr. Caro and colleagues to
estimate censored QUAIDS models, I found that the long computation times—
especially under IFGNLS—could become a serious obstacle to model diagnostics
and to exploring how results changed under alternative specifications. At the
same time, Python offered a more convenient environment for building research
pipelines and for running computational work on virtual servers. These two
considerations motivated me to develop `pyquaidsce`: a Python implementation
aimed at preserving the econometric logic of the original command while making
censored QUAIDS estimation more computationally practical and easier to
integrate into modern research workflows.

## Why pyquaidsce?

- Native Python workflow using only `numpy`, `scipy`, and `pandas` at runtime.
- Censored AIDS/QUAIDS with Ray (1983) demographic scaling.
- `nls`, `fgnls`, and iterated `ifgnls` estimation.
- Stata-oriented parameter ordering, result tables, and compatibility options.
- Bootstrap support with parallel replications.
- Analytic Jacobians and regression tests against a full-precision Stata benchmark.
- A public, reproducible Stata/Python benchmark using the same synthetic CQUAIDS
  data and model specification, with raw logs and full-precision result comparison.

## Installation

Clone or download this repository, then install it from the repository root:

```bash
python -m pip install .
```

For development or validation work, use an editable installation:

```bash
python -m pip install -e .
```

The project currently declares support for Python 3.9 or newer.

## Quick start

```python
import pandas as pd
from pyquaidsce import quaidsce

# df = pd.read_csv("mydata.csv")

res = quaidsce(
    df,
    shares=["w1", "w2", "w3", "w4"],
    prices=["p1", "p2", "p3", "p4"],
    expenditure="total_expenditure",
    demographics=["household_size", "age_head"],
    anot=10.0,
    method="ifgnls",
    algorithm="gn",
    reps=0,
    first_stage_predict="pr",
    strict_stata=True,
    verbose=True,
)

print(res.summary())
print(res.elasticity_tables())
```

The price variables must be in the same order as the budget-share variables.
For a censored model, the demographic variables are Ray-scaling variables; they
are not generic additive controls.

For a Stata-to-Python option map and a complete first workflow, see
[`docs/getting-started.md`](docs/getting-started.md).

## Stata compatibility

Two options are especially important:

- `first_stage_predict="pr"` reproduces the first-stage prediction path used by
  the shipped Stata implementation. `first_stage_predict="xb"` uses the linear
  predictor inside the Shonkwiler–Yen normal CDF/PDF terms.
- `strict_stata=True` preserves documented Stata elasticity behavior needed for
  direct reproduction. Use `strict_stata=False` when you deliberately want the
  corresponding corrected formula documented in this project.

These choices are methodological, not cosmetic. See
[`docs/stata-compatibility.md`](docs/stata-compatibility.md) before using the
package in new empirical work.

## Validation and controlled benchmark

The public end-to-end benchmark uses a deterministic synthetic dataset with
**20,000 observations, 4 goods, 3 demographic variables, censoring in every
good, and IFGNLS**. Stata and Python use the same `.dta` file, neither receives
a user-supplied initial parameter vector, and bootstrap is disabled.

The Stata and Python results are **almost identical**:

| Quantity | Maximum/relative difference |
|---|---:|
| structural parameters (`alpha` through `rho`) | `1.68e-05` |
| first-stage Probit parameters (`tau`) | `1.21e-07` |
| expenditure/Marshallian/Hicksian elasticities | `7.34e-07` |
| non-elasticity standard errors | `7.25e-06` |
| log-likelihood, relative difference | `4.99e-08` |

In the recorded **same-machine** runtime comparison, Stata 19.5 took
**1,161.171 seconds (19 min 21.171 s)** and `pyquaidsce` under Python 3.14.5
took **26.033 seconds**. This corresponds to a **44.60x wall-clock speedup**
for this benchmark.

Both timings cover the point-estimation call rather than data loading. The
Stata log reports two processors available to Stata. The recorded Python log
does not contain an explicit BLAS/OpenMP thread-limit field, so the 44.60x
figure is a same-machine wall-clock comparison, not a per-core efficiency
claim.

All data, scripts, raw logs, machine-readable outputs, and the row-by-row
comparison are in
[`benchmarks/cquaids_ifgnls_4g_20k/`](benchmarks/cquaids_ifgnls_4g_20k/). See
[`docs/validation.md`](docs/validation.md) and
[`docs/performance.md`](docs/performance.md) for details.

## Repository guide

```text
src/pyquaidsce/   Python package source
tests/            automated mathematical and regression tests
examples/         user-facing Python and Stata examples
benchmarks/       controlled public performance benchmark and raw results
tools/           validation and diagnostic utilities
docs/            methodology, compatibility, performance, and validation notes
.github/           CI and issue templates
```

Generated wheels and local audit outputs are intentionally not committed. A
wheel should be attached to a versioned GitHub Release instead.

## Documentation

- [Getting started](docs/getting-started.md)
- [Methodology and estimator](docs/methodology.md)
- [Stata compatibility and known implementation differences](docs/stata-compatibility.md)
- [Validation evidence](docs/validation.md)
- [Performance notes](docs/performance.md)
- [Contributing and reproducible bug reports](CONTRIBUTING.md)

## Citation and upstream credit

If you use `pyquaidsce`, cite the software release you used and also cite the
original `quaidsce` work:

> Caro, J. C., Melo, G., Molina, J. A., and Salgado, J. C. (2021),
> *Censored QUAIDS estimation with quaidsce*, Boston College Working Papers in
> Economics 1045.

The Stata module is distributed through SSC/RePEc as Statistical Software
Components S459029. This Python project is a reimplementation and is **not an
official release of the original Stata authors**.

A `CITATION.cff` file is included so GitHub can expose a “Cite this repository”
action.

## Author and maintainer

**Sina Amiri**  
Department of Economics, Shiraz University, Shiraz, Iran.

At the time of the version 1.0.1 public release, Sina Amiri is an undergraduate
student in Economics at Shiraz University.

## Development transparency

The development and review process used AI-assisted coding tools. The public
release is therefore organized around reproducible tests, explicit numerical
benchmarks, and documented econometric choices rather than relying on code
generation provenance as evidence of correctness.

## License

GPL-3.0-only. See [`LICENSE`](LICENSE). The upstream Stata module is also distributed
under GPL v3 through SSC/RePEc.
