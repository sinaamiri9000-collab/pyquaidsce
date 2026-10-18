"""Check semantic elasticity labels independently of numerical benchmarks."""

from pathlib import Path
import pickle
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from pyquaidsce import quaidsce
from pyquaidsce.bootstrap import BootResult
from pyquaidsce.elasticities import Elasticities
from pyquaidsce.params import Spec
from pyquaidsce.stata_bridge import _BOOT_STATE, load_stata_results
from tools.legacy_names import canonical_legacy_name
from tools.validate_small4 import DEMOS, PRICES, SHARES, parse


class ElasticityNameTests(unittest.TestCase):
    def test_asymmetric_matrices_match_every_named_cell(self):
        for n in (3, 4, 11):
            for quadratic in (False, True):
                for control_function in (False, True):
                    with self.subTest(n=n, quadratic=quadratic, cf=control_function):
                        spec = Spec(n, 1, quadratic, True, control_function)
                        u = np.arange(n * n, dtype=float).reshape(n, n) + 10
                        c = u + 1000
                        income = np.arange(n, dtype=float) + 0.5
                        el = Elasticities(income, u, c, income, u, np.ones(n))
                        named = dict(zip(spec.elas_names(), el.as_stata_vector()))
                        self.assertEqual(len(named), n + 2 * n * n)
                        for i in range(n):
                            self.assertEqual(named[f"ELAS_INC:e_{i+1}"], income[i])
                            for j in range(n):
                                self.assertEqual(named[f"ELAS_UNCOMP:e_{i+1}_{j+1}"], u[i, j])
                                self.assertEqual(named[f"ELAS_COMP:e_{i+1}_{j+1}"], c[i, j])

    def test_legacy_conversion_is_explicit_and_block_specific(self):
        self.assertEqual(canonical_legacy_name("ELAS_UNCOMP:e_2_1"), "ELAS_UNCOMP:e_1_2")
        self.assertEqual(canonical_legacy_name("ELAS_COMP:e_11_2"), "ELAS_COMP:e_2_11")
        for name in ("ELAS_INC:e_1", "gamma:gamma_2_1", "ELAS_COMP:e_2_2"):
            self.assertEqual(canonical_legacy_name(name), name)

    def test_archived_reader_normalizes_values_and_errors_together(self):
        path = str(ROOT / "tests" / "data" / "small4.log")
        raw, _ = parse(path, legacy_elasticity_labels=False)
        canonical, _ = parse(path)
        for run, rows in raw.items():
            for key, value in rows.items():
                self.assertEqual(canonical[run][canonical_legacy_name(key)], value)
        self.assertNotEqual(raw[1]["ELAS_UNCOMP:e_2_1"][0],
                            raw[1]["ELAS_UNCOMP:e_1_2"][0])


class FittedElasticityNameTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame = pd.read_stata(ROOT / "tests" / "data" / "small4.dta")
        cls.kwargs = dict(shares=SHARES, prices=PRICES, expenditure="total",
                          demographics=DEMOS, anot=10.0, verbose=False)
        cls.res = quaidsce(cls.frame, **cls.kwargs)

    def test_result_lookup_matches_matrices(self):
        self.assertTrue(self.res.converged)
        for block, matrix in (("ELAS_UNCOMP", self.res.elas.uncompensated),
                              ("ELAS_COMP", self.res.elas.compensated)):
            for i in range(4):
                for j in range(4):
                    key = f"{block}:e_{i+1}_{j+1}"
                    self.assertEqual(self.res.get(key), matrix[i, j])
                    self.assertEqual(self.res.named()[key], matrix[i, j])

    def test_linear_and_uncensored_layouts(self):
        linear = quaidsce(self.frame, **self.kwargs, quadratic=False)
        for i in range(4):
            for j in range(4):
                self.assertEqual(linear.get(f"ELAS_UNCOMP:e_{i+1}_{j+1}"),
                                 linear.elas.uncompensated[i, j])
        uncensored = quaidsce(self.frame, **self.kwargs, censor=False)
        self.assertEqual(len(uncensored.b), uncensored.spec.n_full)
        self.assertFalse(any(name.startswith("ELAS_") for name in uncensored.names))

    def test_stata_names_and_bootstrap_covariance_share_column_positions(self):
        res = self.res
        n, k = res.spec.neqn, len(res.b)
        # Distinct uncertainty in every column, independently indexed by [i,j].
        draws = res.b + np.arange(5)[:, None] * (np.arange(k)[None, :] + 1) / 100
        boot = BootResult(5, 5, draws, draws.std(axis=0, ddof=1))
        stored = {}

        class Matrix:
            @staticmethod
            def create(name, rows, cols, value):
                stored[name] = np.full((rows, cols), value, dtype=float)

            @staticmethod
            def storeAt(name, row, col, value):
                stored[name][row, col] = value

            @staticmethod
            def setColNames(name, names):
                stored[name + "_cols"] = list(names)

            @staticmethod
            def setRowNames(name, names):
                stored[name + "_rows"] = list(names)

        fake_sfi = types.SimpleNamespace(
            Matrix=Matrix,
            Scalar=types.SimpleNamespace(setValue=lambda name, value: None),
            Macro=types.SimpleNamespace(setLocal=lambda name, value: None),
        )
        payload = dict(b=res.b, V=boot.V, names=res.names, shares=res.share_names,
                       elas_i=res.elas.income, elas_u=res.elas.uncompensated,
                       elas_c=res.elas.compensated, nobs=res.nobs, llf=res.llf,
                       anot=res.anot, ndemo=res.spec.ndemo, converged=1,
                       n_outer=res.n_outer, n_gn=res.n_gn, title=res._title(),
                       boot_reps_ok=5, boot_reps_requested=5)
        with tempfile.TemporaryDirectory() as temp:
            path = str(Path(temp) / "result.pkl")
            with open(path, "wb") as handle:
                pickle.dump(payload, handle)
            with patch.dict(_BOOT_STATE, {"result_path": path}, clear=True), \
                    patch.dict(sys.modules, {"sfi": fake_sfi}):
                load_stata_results("b", "V", "ei", "eu", "ec")
        self.assertEqual(stored["b_cols"], stored["V_cols"])
        self.assertEqual(stored["b_cols"], stored["V_rows"])
        start = k - (n + 2 * n * n)
        for offset, block, matrix_name in ((n, "ELAS_UNCOMP", "eu"),
                                            (n + n * n, "ELAS_COMP", "ec")):
            for i in range(n):
                for j in range(n):
                    key = f"{block}:e_{i+1}_{j+1}"
                    position = start + offset + i * n + j
                    self.assertEqual(stored["b_cols"][position], key)
                    self.assertEqual(stored["b"][0, position], stored[matrix_name][i, j])
                    self.assertAlmostEqual(np.sqrt(stored["V"][position, position]),
                                           boot.se[position], places=12)


if __name__ == "__main__":
    unittest.main()
