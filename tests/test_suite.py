"""Collectable regression tests (standard-library ``unittest``)."""

from __future__ import annotations

import os
import sys
import unittest

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "src")
sys.path.insert(0, ROOT)
sys.path.insert(0, SRC)

from pyquaidsce import quaidsce  # noqa: E402
from tests.test_internals import (check_delta_matrix, check_fast_jacobian,  # noqa: E402
                                  check_jacobian, check_probit,
                                  check_restrictions)
from tools.validate_small4 import (BENCH, DEMOS, PRICES, SHARES, parse)  # noqa: E402


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
        cls.coefs, cls.lls = parse(os.path.join(BENCH, "small4.log"))
        cls.df = pd.read_stata(os.path.join(BENCH, "small4.dta"))

    def test_ifgnls_is_the_public_default_and_runs_end_to_end(self):
        res = quaidsce(
            self.df,
            shares=SHARES,
            prices=PRICES,
            expenditure="total",
            demographics=DEMOS,
            anot=10.0,
            verbose=False,
        )
        self.assertEqual(res.method, "ifgnls")
        self.assertTrue(res.converged)
        self.assertTrue(np.isfinite(res.b).all())
        self.assertTrue(np.isfinite(res.V).all())
        n_elasticities = res.spec.neqn + 2 * res.spec.neqn ** 2
        self.assertTrue(np.isnan(res.se[-n_elasticities:]).all())

    def test_invalid_options_and_removed_legacy_switches_fail_loudly(self):
        common = dict(
            data=self.df,
            shares=SHARES,
            prices=PRICES,
            expenditure="total",
            demographics=DEMOS,
            anot=10.0,
            verbose=False,
        )
        with self.assertRaisesRegex(ValueError, "algorithm"):
            quaidsce(**common, algorithm="not-an-algorithm")
        with self.assertRaises(TypeError):
            quaidsce(**common, first_stage_predict="xb")
        with self.assertRaises(TypeError):
            quaidsce(**common, strict_stata=False)
        with self.assertRaisesRegex(ValueError, "initial must contain"):
            quaidsce(**common, initial=np.zeros(2))


if __name__ == "__main__":
    unittest.main()
