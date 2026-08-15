"""Assemble a clean, checksumed 1.2.0 release-review bundle."""

from __future__ import annotations

import hashlib
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT_PARENT = ROOT.parent
RELEASE_NAME = "PYQUAIDSCE_1.2.0_FINAL_RELEASE"
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
        raise FileExistsError(
            "release output already exists; move it aside before rebuilding"
        )
    dist = ROOT / "dist"
    wheel = dist / "pyquaidsce-1.2.0-py3-none-any.whl"
    sdist = dist / "pyquaidsce-1.2.0.tar.gz"
    for artifact in (wheel, sdist):
        if not artifact.is_file():
            raise FileNotFoundError(f"missing built artifact: {artifact}")

    clean_source = RELEASE_DIR / "source/pyquaidsce-1.2.0"
    ignore = shutil.ignore_patterns(
        ".git",
        "__pycache__",
        "*.pyc",
        "*.pyo",
        "build",
        "dist",
        "dist_*",
        "wheel_check",
        "work_validation",
        "*.egg-info",
    )
    shutil.copytree(ROOT, clean_source, ignore=ignore)

    out_dist = RELEASE_DIR / "dist"
    out_dist.mkdir(parents=True)
    shutil.copy2(wheel, out_dist / wheel.name)
    shutil.copy2(sdist, out_dist / sdist.name)

    source_zip = out_dist / "pyquaidsce-1.2.0-source.zip"
    with zipfile.ZipFile(source_zip, "w", zipfile.ZIP_DEFLATED) as archive:
        add_tree(archive, clean_source, "pyquaidsce-1.2.0")

    shutil.copy2(ROOT / "MERGE_AUDIT_FA.md", RELEASE_DIR / "MERGE_AUDIT_FA.md")
    qa_dir = RELEASE_DIR / "QA_RESULTS"
    qa_dir.mkdir()
    shutil.copy2(
        ROOT / "benchmarks/release_120/results/bootstrap_smoke.json",
        qa_dir / "bootstrap_smoke.json",
    )
    shutil.copy2(
        ROOT / "benchmarks/release_120/results/bootstrap_determinism.json",
        qa_dir / "bootstrap_determinism.json",
    )

    artifacts = sorted(out_dist.iterdir())
    lines = [f"{sha256(path)}  dist/{path.name}" for path in artifacts]
    (RELEASE_DIR / "SHA256SUMS.txt").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )

    with zipfile.ZipFile(BUNDLE_ZIP, "w", zipfile.ZIP_DEFLATED) as archive:
        add_tree(archive, RELEASE_DIR, RELEASE_NAME)
    print(BUNDLE_ZIP)


if __name__ == "__main__":
    main()
