"""
Standalone estimator runner for out-of-process Stata execution.

Invoked as::

    python -m pyquaidsce.bootstrap_runner --input job.pkl

The input pickle contains:

- ``data``         – pandas DataFrame
- ``kwargs``       – dict of ``quaidsce()`` keyword arguments (including bootstrap controls)
- ``result_path``  – where to write the result pickle
- ``status_path``  – where to write progress updates

This module is launched by :func:`pyquaidsce.stata_bridge.launch_from_stata` as
a separate process so that Stata's GUI thread is never blocked.
"""

from __future__ import annotations

import argparse
import os
import pickle
import sys
import traceback


class _ProgressWriter:
    """Intercept output and mirror the latest progress line to a status file."""

    def __init__(self, status_path: str, original):
        self.status_path = status_path
        self.original = original

    def write(self, s: str) -> int:
        if self.original is not None:
            self.original.write(s)
        stripped = s.strip()
        if stripped:
            try:
                with open(self.status_path, "w", encoding="utf-8") as f:
                    f.write(stripped)
            except OSError:
                pass
        return len(s)

    def flush(self):
        if self.original is not None:
            self.original.flush()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run pyquaidsce estimation in a standalone process."
    )
    parser.add_argument(
        "--input", required=True, help="Path to the job pickle file"
    )
    args = parser.parse_args()

    # ---- load job --------------------------------------------------------- #
    with open(args.input, "rb") as fh:
        job = pickle.load(fh)

    data = job["data"]
    kwargs = job["kwargs"]
    result_path = job["result_path"]
    status_path = job["status_path"]
    error_path = result_path + ".error"

    # ---- redirect stdout for progress ------------------------------------- #
    original_stdout = sys.stdout
    sys.stdout = _ProgressWriter(status_path, original_stdout)

    try:
        from pyquaidsce.estimator import quaidsce

        with open(status_path, "w", encoding="utf-8") as file_handle:
            file_handle.write("Estimating point estimates...")
        res = quaidsce(data=data, **kwargs)

        result = {
            "b": res.b,
            "V": res.V,
            "names": list(res.names),
            "shares": list(res.share_names),
            "elas_i": res.elas.income,
            "elas_u": res.elas.uncompensated,
            "elas_c": res.elas.compensated,
            "nobs": int(res.nobs),
            "llf": float(res.llf),
            "anot": float(res.anot),
            "ndemo": int(res.spec.ndemo),
            "converged": 1 if res.converged else 0,
            "n_outer": int(res.n_outer),
            "n_gn": int(res.n_gn),
            "title": res._title(),
            "boot_reps_ok": res.boot.reps_ok if res.boot else 0,
            "boot_reps_requested": res.boot.reps_requested if res.boot else 0,
            "failures": list(res.boot.failures) if res.boot else [],
            # Compatibility keys used by load_bootstrap_results().
            "reps_ok": res.boot.reps_ok if res.boot else 0,
            "reps_requested": res.boot.reps_requested if res.boot else 0,
        }
        # Atomic-ish write: write to a temp file, then rename
        tmp_path = result_path + ".tmp"
        with open(tmp_path, "wb") as fh:
            pickle.dump(result, fh)
        # On Windows, os.replace is atomic within the same volume.
        os.replace(tmp_path, result_path)

    except Exception as exc:
        with open(error_path, "w", encoding="utf-8") as fh:
            fh.write(f"{type(exc).__name__}: {exc}\n")
            traceback.print_exc(file=fh)
        sys.exit(1)
    finally:
        sys.stdout = original_stdout
        # Clean up the job file (data can be large)
        try:
            os.unlink(args.input)
        except OSError:
            pass


if __name__ == "__main__":
    main()
