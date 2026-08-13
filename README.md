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

## Why pyquaidsce?

- Native Python workflow using only `numpy`, `scipy`, and `pandas` at runtime.
- Censored AIDS/QUAIDS with Ray (1983) demographic scaling.
- `nls`, `fgnls`, and iterated `ifgnls` estimation.
- Stata-oriented parameter ordering, result tables, and compatibility options.
- Bootstrap support with parallel replications.
- Analytic Jacobians and regression tests against a full-precision Stata benchmark.
- A measured 11-good, 37,000-observation, 3-demographic `ifgnls` run completed in
  about 10 minutes in the recorded two-core benchmark. Runtime depends strongly
  on hardware, scaling, convergence path, and the number of free parameters.

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

## Validation

The repository includes a 4-good, 2,000-observation benchmark with a Stata log
printed at high precision. For the `ifgnls` specification, the recorded
comparison gives approximately:

- relative log-likelihood difference: `8.5e-09`;
- maximum absolute coefficient difference: `8.4e-07`;
- maximum absolute standard-error difference: `8.6e-08`.

The test suite also checks the parameter-restriction map, analytic Jacobians,
and the first-stage probit independently of Stata.

Run the tests with:

```bash
python -m unittest discover -s tests -v
```

Full evidence and reproduction instructions are in
[`docs/validation.md`](docs/validation.md).

## Repository guide

```text
src/pyquaidsce/   Python package source
tests/            automated mathematical and regression tests
examples/         user-facing Python and Stata examples
bench/            small reproducible benchmark data and Stata logs
tools/            validation and diagnostic utilities
docs/             methodology, compatibility, performance, and validation notes
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
