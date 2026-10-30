# Software tests

This directory contains developer checks for the correctness, numerical stability,
and consistent behavior of `pyquaidsce`. They check software behavior; they are
not instructions for choosing an empirical specification or statistical tests
for a research application.

If you are a user or researcher who does not intend to develop the package,
start with the [Getting Started guide](../docs/getting-started.md), the
[User Guide and option reference](../docs/user-guide.md), and the runnable
[examples](../examples/README.md). The [Stata guide](../stata/README.md) and
[R guide](../R/README.md) explain those interfaces.

## Automated Python tests

From the repository root, install the package and run:

```bash
python -m pip install -e .
python -m unittest discover -s tests -v
```

The tests cover:

| Files | Checks |
|---|---|
| `test_suite.py`, `test_internals.py` | Parameter restrictions and maps, analytic versus numerical derivatives, Probit calculations, and end-to-end estimation. The suite calls the helper functions in `test_internals.py`. |
| `test_control_function.py`, `test_ivexp.py`, `test_uncensored_control_function.py` | External and integrated expenditure controls, derivatives, input validation, and rebuilding first-stage regressions in bootstrap draws. |
| `test_convergence_and_nls_covariance.py` | Stopping criteria and covariance formulas. |
| `test_elasticity_names.py` | Agreement between elasticity names, matrix cells, and stored outputs. |
| `test_stata_launcher.py`, `test_release_160.py` | Python-side Stata bridge, subprocess handling, and integration behavior. Stata components are simulated where needed; these tests do not launch Stata itself. |
| `test_threads_and_defaults_160.py` | Numerical defaults and scoped thread limits. |

The `160` suffixes are retained filenames: these are active regression tests,
not a requirement to install version 1.6.0.

[GitHub Actions](../.github/workflows/tests.yml) runs this discovery command on
Ubuntu and Windows with Python 3.9 through 3.13. Discovery does not execute every
script in this directory as a standalone program.

## Additional manual checks

- `test_theory.py`: optional diagnostics for demand-theory identities,
  uncensored QUAIDS, censored AIDS, and bootstrap inference. Run
  `python tests/test_theory.py` from the repository root. It prints diagnostics
  and is not a collected `unittest.TestCase` suite.
- `stata/test_stata_small4.do`: end-to-end Stata point-estimate and elasticity
  output checks, with an optional bootstrap check.
- `stata/test_uncensored_ivexp.do`: Stata checks for uncensored expenditure
  controls and adding-up restrictions.

The `.do` files require a licensed Stata installation and the matching Python
backend. Run them from the repository root; the Python CI does not run them.

## Test data and shared helpers

`data/small4.dta` and `data/small4.log` are small reference fixtures used by the
Python tests. `data/benchmark_cquaids_4g_20k.dta` supports the optional bootstrap
diagnostics in [tools](../tools/README.md). These fixtures are for verification,
not a required input format or research dataset for package users.

Some tests import shared helpers from `tools/validate_small4.py`,
`tools/legacy_names.py`, and `tools/paths.py`. Keep those helpers available when
running the suite.
