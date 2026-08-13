# Contributing to pyquaidsce

Contributions, feedback, and bug reports are warmly welcome! Whether you are an economist testing a new empirical specification, finding an edge case, or suggesting performance improvements, your input is appreciated.

---

## 1. Development Setup

Clone the repository and set up an editable environment with tests:

```bash
git clone https://github.com/sinaamiri9000-collab/pyquaidsce.git
cd pyquaidsce
pip install -e .
python -m unittest discover -s tests -v
```

---

## 2. Reporting Issues

When submitting a bug report or numerical question, please include:
- Your Python version and OS.
- Installed versions of `numpy`, `scipy`, and `pandas`.
- Model specification: number of goods, demographics, and options used (`method`, `first_stage_predict`, `strict_stata`).
- A minimal reproducible example or anonymized sample data.
- The corresponding Stata command and output if reporting a discrepancy with Stata.

> [!NOTE]
> Please do not upload private or confidential microdata to public issues.

---

## 3. Pull Requests

If proposing modifications to the econometric formulas, Jacobians, or optimization logic:
- Include corresponding unit tests to verify mathematical correctness.
- Ensure all existing tests pass (`python -m unittest discover -s tests`).
- Keep code clean and well-documented.
