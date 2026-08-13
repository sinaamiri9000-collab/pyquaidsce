# Public-release checklist

This repository draft is intentionally usable before a GitHub account/repository
URL or personal citation metadata is known. Before making the repository public:

- [ ] Replace `authors: - name: "pyquaidsce project"` in `CITATION.cff` with the
      final software author/maintainer name(s), ORCID(s), and affiliation(s) as
      appropriate.
- [ ] Decide whether version 1.0.1 should be tagged exactly as `v1.0.1` and keep
      the package version synchronized with that tag.
- [ ] Create the GitHub repository with the name `pyquaidsce` if available.
- [ ] Upload this source tree; do **not** commit `dist/` or wheels.
- [ ] Confirm the GitHub Actions test workflow is green on all configured Python
      versions.
- [ ] Create a GitHub Release for `v1.0.1` and attach the built wheel/source
      archive to the release, rather than tracking them in Git.
- [ ] If publishing to PyPI, first verify the project name and release on
      TestPyPI; then use a trusted-publishing workflow rather than a long-lived
      API token.
- [ ] Optionally connect GitHub to Zenodo after the first public release to
      obtain a versioned DOI for research citation.
- [ ] Re-read README claims about performance and compatibility and keep them
      tied to the documented benchmarks rather than presenting a universal speed
      ratio.
