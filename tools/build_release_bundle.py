"""Build wheel/sdist and assemble a clean pyquaidsce 1.6.0 release bundle."""

from __future__ import annotations

import hashlib
import shutil
import subprocess
import sys
import zipfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.6.0"
OUT_PARENT = ROOT.parent
RELEASE_NAME = f"PYQUAIDSCE_{VERSION}_FINAL_RELEASE"
RELEASE_DIR = OUT_PARENT / RELEASE_NAME
BUNDLE_ZIP = OUT_PARENT / f"{RELEASE_NAME}_BUNDLE.zip"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def add_tree(archive: zipfile.ZipFile, tree: Path, prefix: str) -> None:
    for path in sorted(tree.rglob("*")):
        if path.is_file():
            archive.write(path, Path(prefix) / path.relative_to(tree))


def main() -> None:
    if RELEASE_DIR.exists() or BUNDLE_ZIP.exists():
        raise FileExistsError("release output already exists; move it aside before rebuilding")

    subprocess.check_call([sys.executable, "-m", "build", "--wheel", "--sdist"], cwd=ROOT)
    dist = ROOT / "dist"
    wheel = dist / f"pyquaidsce-{VERSION}-py3-none-any.whl"
    sdist = dist / f"pyquaidsce-{VERSION}.tar.gz"
    for artifact in (wheel, sdist):
        if not artifact.is_file():
            raise FileNotFoundError(f"missing built artifact: {artifact}")

    clean_source = RELEASE_DIR / "source" / f"pyquaidsce-{VERSION}"
    ignore = shutil.ignore_patterns(
        ".git", "__pycache__", "*.pyc", "*.pyo", "build", "dist",
        "dist_*", "wheel_check", "work_validation", "*.egg-info",
        ".pytest_cache", ".Rcheck",
    )
    shutil.copytree(ROOT, clean_source, ignore=ignore)

    out_dist = RELEASE_DIR / "dist"
    out_dist.mkdir(parents=True)
    shutil.copy2(wheel, out_dist / wheel.name)
    shutil.copy2(sdist, out_dist / sdist.name)

    source_zip = out_dist / f"pyquaidsce-{VERSION}-source.zip"
    with zipfile.ZipFile(source_zip, "w", zipfile.ZIP_DEFLATED) as archive:
        add_tree(archive, clean_source, f"pyquaidsce-{VERSION}")

    artifacts = sorted(out_dist.iterdir())
    (RELEASE_DIR / "SHA256SUMS.txt").write_text(
        "\n".join(f"{sha256(path)}  dist/{path.name}" for path in artifacts) + "\n",
        encoding="utf-8",
    )

    with zipfile.ZipFile(BUNDLE_ZIP, "w", zipfile.ZIP_DEFLATED) as archive:
        add_tree(archive, RELEASE_DIR, RELEASE_NAME)
    print(BUNDLE_ZIP)


if __name__ == "__main__":
    main()
