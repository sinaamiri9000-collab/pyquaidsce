"""Collectable regression tests (standard-library ``unittest``)."""

from __future__ import annotations

import os
import re
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce
from tests.test_internals import (check_delta_matrix, check_fast_jacobian,
                                  check_jacobian, check_probit,
                                  check_restrictions)

ROOT = Path(__file__).resolve().parents[1]
BENCH = ROOT / "benchmarks" / "cquaids_ifgnls_4g_20k"
DATA = BENCH / "data" / "benchmark_cquaids_4g_20k.dta"
STATA_LOG = BENCH / "results" / "stata_benchmark.log"
SHARES = ["w1", "w2", "w3", "w4"]
PRICES = ["p1", "p2", "p3", "p4"]
DEMOS = ["z1", "z2", "z3"]


def parse_stata_reference(path: Path):
    """Read the full-precision BENCHPARM/BENCHLL records from the stored log."""
    values = {}
    ll = None
    parm = re.compile(r"^BENCHPARM\s+(\S+)\s+([.0-9Ee+\-]+)\s+([.0-9Ee+\-]+)")
    for raw in path.read_text(errors="replace").splitlines():
        line = raw.strip()
        m = parm.match(line)
        if m:
            values[m.group(1)] = (float(m.group(2)), float(m.group(3)))
        elif line.startswith("BENCHLL "):
            ll = float(line.split()[1])
    if not values or ll is None:
        raise RuntimeError(f"could not parse benchmark reference: {path}")
    return values, ll


class InternalMathTests(unittest.TestCase):
    def test_delta_map(self):
        self.assertLess(check_delta_matrix(), 2e-8)

    def test_restrictions(self):
        self.assertTrue(check_restrictions())

    def test_full_jacobian(self):
        self.assertLess(check_jacobian(), 1e-7)

    def test_fast_jacobian(self):
        self.assertLess(check_fast_jacobian(), 1e-12)

    def test_probit(self):
        db, dll = check_probit()
        self.assertLess(db, 1e-6)
        self.assertLess(abs(dll), 1e-9)


class StataBenchmarkTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.ref, cls.ll = parse_stata_reference(STATA_LOG)
        cls.df = pd.read_stata(DATA)

    def test_ifgnls_end_to_end_is_almost_identical_to_stata(self):
        res = quaidsce(
            self.df,
            shares=SHARES,
            prices=PRICES,
            expenditure="total",
            demographics=DEMOS,
            anot=10.0,
            method="ifgnls",
            start="zero",
            first_stage_predict="pr",
            strict_stata=True,
            reps=0,
            verbose=False,
        )
        b_ref = np.array([self.ref[name][0] for name in res.names])
        self.assertTrue(res.converged)
        self.assertLess(abs(res.llf - self.ll) / abs(self.ll), 1e-7)
        self.assertLess(float(np.max(np.abs(res.b - b_ref))), 2e-5)

    def test_invalid_switches_fail_loudly(self):
        common = dict(
            data=self.df.iloc[:200].copy(),
            shares=SHARES,
            prices=PRICES,
            expenditure="total",
            demographics=DEMOS,
            anot=10.0,
            verbose=False,
        )
        with self.assertRaisesRegex(ValueError, "algorithm"):
            quaidsce(**common, algorithm="not-an-algorithm")
        with self.assertRaisesRegex(ValueError, "first_stage_predict"):
            quaidsce(**common, first_stage_predict="probability")
        with self.assertRaisesRegex(ValueError, "initial must contain"):
            quaidsce(**common, initial=np.zeros(2))


if __name__ == "__main__":
    unittest.main()
