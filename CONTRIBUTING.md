# Contributing

Contributions are welcome, especially reproducible bug reports, numerical
validation, documentation improvements, and tests that clarify econometric
behavior.

## Development setup

```bash
python -m venv .venv
python -m pip install -e .
python -m unittest discover -s tests -v
```

Use a dedicated virtual environment. Do not commit local data, virtual
environments, `dist/`, or generated audit output.

## Bug reports

For a numerical discrepancy, include as much of the following as possible:

- Python version and operating system;
- `numpy`, `scipy`, and `pandas` versions;
- number of goods, observations, and demographics;
- `method`, `algorithm`, `start`, `first_stage_predict`, `strict_stata`, and
  convergence settings;
- final estimation-sample size and the exact preprocessing filters;
- a minimal synthetic or anonymized dataset if redistribution is permitted;
- the corresponding Stata command and relevant log excerpt when reporting a
  Stata/Python difference.

Do not upload confidential microdata to a public issue.

## Pull requests

A change to model equations, parameter restrictions, Jacobians, elasticities,
first-stage correction, covariance construction, or convergence logic must come
with a test that would fail under the previous incorrect behavior. Compatibility
changes should also state whether they alter direct Stata reproduction.

Run the full test suite before opening a pull request.

## AI-assisted contributions

AI-assisted coding is acceptable, but generated code is not evidence of
correctness. The contributor remains responsible for understanding the change,
checking licensing/provenance, and supplying tests or derivations appropriate to
the scientific claim.
