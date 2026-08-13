# First publication on GitHub

This project folder is already organized for GitHub. GitHub does not need the
ZIP file itself; the ZIP is only a convenient way to transfer the folder.

## Recommended beginner workflow: GitHub Desktop

1. Extract `pyquaidsce-v1.0.1-github-ready.zip` to a temporary folder.
2. Install GitHub Desktop and sign in to your GitHub account.
3. In GitHub Desktop, create a new local repository named `pyquaidsce` in a
   convenient location such as `Documents/GitHub`.
4. Copy **the contents** of the extracted project folder into the new local
   `pyquaidsce` repository folder. The top-level `README.md`, `pyproject.toml`,
   `src`, `tests`, `docs`, `.github`, etc. should sit directly in the repository
   root; do not create an extra nested `pyquaidsce-v1.0.1-github-ready` folder.
5. Return to GitHub Desktop. It will show the new/changed files.
6. Use the commit summary `Initial public release of pyquaidsce 1.0.1` and
   commit the changes to `main`.
7. Click **Publish repository**.
8. Use repository name `pyquaidsce` and a concise description such as:
   `Censored QUAIDS estimation in Python with Stata-compatible validation.`
9. Keep the repository private while checking the landing page if you want a
   final review before launch. When ready, change visibility to public.
10. Create a version tag/release `v1.0.1` only after the public files and
    citation metadata have been checked.

## Important Git principle

Git does not record every keystroke. It records commits. Changes edited only on
your computer are private. Once a commit is pushed to a public repository, its
committed diff and history can normally be viewed by other people. This is a
feature, not a problem: normal README corrections and documentation improvements
are expected in a healthy software project.

Never commit passwords, API keys, private survey data, or other secrets. Deleting
such material in a later commit does not make the earlier public commit safe.
