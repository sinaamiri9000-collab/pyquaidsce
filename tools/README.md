# Developer validation tools

This directory contains software-validation utilities and shared test helpers.
They help maintainers check numerical behavior, bootstrap stability, and stored
reference outputs. Automated regression tests live in [tests](../tests/README.md).

If you are a user or researcher who does not intend to develop `pyquaidsce`,
use the [Getting Started guide](../docs/getting-started.md), the
[User Guide and option reference](../docs/user-guide.md), and the runnable
[examples](../examples/README.md). Interface-specific instructions are in the
[Stata guide](../stata/README.md) and [R guide](../R/README.md).
For the reproducible comparisons used in the paper, see [benchmark](../benchmark/).

## Contents

| File | Purpose |
|---|---|
| `validate_small4.py` | Reads the archived four-good Stata reference log and provides shared fixture settings. Its optional command-line report compares estimates and reports differences; it does not assert that the current estimator reproduces every historical Stata result. |
| `legacy_names.py` | Converts known historical elasticity labels in archived fixtures to the current good-by-price order. Active naming tests use this helper. |
| `paths.py` | Resolves repository-local fixture paths. |
| `bootstrap_smoke.py` | Runs a small parallel bootstrap check, checks covariance/standard-error synchronization, and exercises an intentionally expired replication timeout. |
| `bootstrap_determinism.py` | Repeats a seeded parallel bootstrap and checks that draw ordering and results are identical in the same environment. |
| `__init__.py` | Makes the shared helpers importable by the test suite. |

## Running the optional diagnostics

From the repository root:

```bash
python -m pip install -e .
python tools/bootstrap_smoke.py
python tools/bootstrap_determinism.py
```

Both scripts use `tests/data/benchmark_cquaids_4g_20k.dta`, with a smaller subset
of its observations, and write reports under `out/validation/`. They perform
estimation and take longer than simple algebraic checks. Their results validate
software behavior, not the empirical identification of a research model.

To inspect differences from the archived four-good Stata output:

```bash
python tools/validate_small4.py --runs 3
```

GitHub Actions discovers the Python tests in `tests`; it does not automatically
run every script in `tools`. However, several collected tests import the fixture
and label helpers listed above, so this directory remains part of the test
infrastructure.
