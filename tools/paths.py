"""Repository-local paths shared by test-fixture and validation helpers."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
BENCH = PROJECT_ROOT / "tests" / "data"
