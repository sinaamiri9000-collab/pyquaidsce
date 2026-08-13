# Public-release checklist

Before making `pyquaidsce` 1.0.1 public:

- [x] Software author recorded as **Sina Amiri** with affiliation:
      Department of Economics, Shiraz University, Shiraz, Iran.
- [ ] Add an ORCID to `CITATION.cff` if you want it included in the first public
      release. This is optional and can also be added later.
- [ ] Set `date-released` in `CITATION.cff` to the actual public release date.
- [x] Use version `1.0.1`; create the Git tag as `v1.0.1` when the release is
      ready.
- [x] Keep wheels and other build outputs out of the Git source tree.
- [ ] Confirm the GitHub Actions test workflow is green on the final commit.
- [ ] Create a GitHub Release for `v1.0.1` and attach the distribution files as
      release assets rather than committing them to the repository.
- [ ] Build both a wheel and a source distribution from the final tagged source
      tree. For this pure-Python package, a single `py3-none-any` wheel is
      sufficient across supported platforms.
- [ ] If publishing to PyPI, test the release process first and prefer PyPI
      Trusted Publishing through a dedicated GitHub Actions release workflow.
- [ ] Optionally connect the public GitHub repository to Zenodo after the first
      release to obtain a versioned DOI.
- [ ] Re-read compatibility and performance claims against the stored controlled
      benchmark before publishing.
