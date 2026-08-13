# Release Checklist

Steps before publishing a new release of `pyquaidsce`:

- [x] Author information confirmed (**Sina Amiri**, Department of Economics, Shiraz University).
- [ ] Ensure version number is updated consistently across `pyproject.toml`, `src/pyquaidsce/__init__.py`, and `CITATION.cff`.
- [ ] Set `date-released` in `CITATION.cff`.
- [ ] Verify that all unit and theory tests pass:
  ```bash
  python -m unittest discover -s tests -v
  ```
- [ ] Build distribution packages:
  ```bash
  python -m build
  ```
- [ ] Create a versioned Git tag (e.g., `git tag v1.0.1`).
- [ ] Create a GitHub Release and attach the built wheel/sdist packages.
- [ ] Publish to PyPI (if applicable).
