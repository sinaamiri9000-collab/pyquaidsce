"""Release-integration checks specific to the merged 1.3.0 tree."""

from __future__ import annotations

import inspect
import os
import pickle
import sys
import tempfile
import time
import types
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

import numpy as np
import pandas as pd

import pyquaidsce
from pyquaidsce import (
    FirstStageLayout,
    fitted_share_derivatives,
    quaidsce,
)
from pyquaidsce.bootstrap import BootResult
from pyquaidsce.stata_bridge import run_from_stata


class MergeIntegrityTests(unittest.TestCase):
    def test_release_exposes_features_from_all_three_branches(self):
        self.assertEqual(pyquaidsce.__version__, "1.4.0")
        self.assertTrue(inspect.isclass(FirstStageLayout))
        self.assertTrue(callable(fitted_share_derivatives))
        self.assertTrue(callable(run_from_stata))
        params = inspect.signature(quaidsce).parameters
        for name in (
            "control_function",
            "selection_prices",
            "selection_expenditure",
            "selection_covariates",
            "selection_control_function",
            "mp_context",
            "rep_timeout",
        ):
            self.assertIn(name, params)
        bridge_params = inspect.signature(run_from_stata).parameters
        for name in (
            "mp_context",
            "rep_timeout",
            "stop_rule",
            "bootstrap_start",
            "control_function",
            "selection_control_function",
            "selection_prices_specified",
            "selection_covariates_specified",
            "selection_expenditure",
        ):
            self.assertIn(name, bridge_params)

    def test_boot_covariance_is_a_regular_matrix(self):
        result = BootResult(
            reps_requested=3,
            reps_ok=3,
            b_star=np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 8.0]]),
            se=np.array([1.0, np.sqrt(28.0 / 3.0)]),
        )
        self.assertEqual(result.V.shape, (2, 2))
        self.assertTrue(np.isfinite(result.V).all())

    def test_timeout_validation_and_expired_internal_deadline(self):
        frame = pd.DataFrame(
            {
                "w1": [0.0, 0.2, 0.3, 0.0],
                "w2": [0.4, 0.0, 0.2, 0.1],
                "w3": [0.6, 0.8, 0.5, 0.9],
                "p1": [1.0, 1.1, 0.9, 1.2],
                "p2": [1.2, 1.0, 1.1, 0.9],
                "p3": [0.8, 0.9, 1.0, 1.1],
                "m": [10.0, 11.0, 9.0, 12.0],
                "z": [0.0, 1.0, 0.5, -0.5],
            }
        )
        common = dict(
            data=frame,
            shares=["w1", "w2", "w3"],
            prices=["p1", "p2", "p3"],
            expenditure="m",
            demographics=["z"],
            anot=1.5,
            verbose=False,
        )
        with self.assertRaisesRegex(ValueError, "rep_timeout"):
            quaidsce(**common, rep_timeout=0)
        with self.assertRaises(TimeoutError):
            quaidsce(**common, _deadline=time.perf_counter() - 1.0)

        watchdog_frame = pd.read_stata(ROOT / "bench" / "small4.dta")
        started = time.perf_counter()
        with self.assertRaisesRegex(RuntimeError, "TimeoutError"):
            quaidsce(
                watchdog_frame,
                shares=["sw1", "sw2", "sw4", "sw9"],
                prices=["p1", "p2", "p4", "p9"],
                expenditure="total",
                demographics=["x1", "x2"],
                anot=10.0,
                method="nls",
                first_stage_predict="xb",
                reps=2,
                seed=12099,
                n_jobs=1,
                mp_context="spawn",
                rep_timeout=1e-6,
                verbose=False,
            )
        # The watchdog must terminate stuck replications promptly, but the
        # absolute wall-clock bound has to scale with the runner: spawning two
        # worker processes on a single-core box contends for the one core and
        # can take ~45 s before either makes progress.
        elapsed_budget = 10.0 if (os.cpu_count() or 1) > 1 else 90.0
        self.assertLess(time.perf_counter() - started, elapsed_budget)

    def test_stata_bridge_forwards_bootstrap_controls(self):
        stored = {}

        class Data:
            values = {
                "_touse": [1, 1],
                "w1": [0.2, 0.3],
                "w2": [0.3, 0.2],
                "w3": [0.5, 0.5],
                "p1": [1.0, 1.1],
                "p2": [1.2, 1.0],
                "p3": [0.9, 1.0],
                "m": [10.0, 11.0],
                "z": [0.0, 1.0],
            }

            @classmethod
            def get(cls, var):
                return cls.values[var]

        class Scalar:
            @staticmethod
            def setValue(name, value):
                stored[name] = value

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

        class Macro:
            @staticmethod
            def setLocal(name, value):
                stored[name] = value

        fake_sfi = types.SimpleNamespace(
            Data=Data, Scalar=Scalar, Matrix=Matrix, Macro=Macro
        )
        captured = {}

        def fake_quaidsce(data, **kwargs):
            captured.update(kwargs)
            return types.SimpleNamespace(
                nobs=2,
                llf=1.0,
                anot=10.0,
                spec=types.SimpleNamespace(ndemo=1),
                converged=True,
                n_outer=2,
                n_gn=3,
                b=np.array([1.0, 2.0]),
                V=np.eye(2),
                boot=None,
                names=["eq1:a", "eq2:b"],
                elas=types.SimpleNamespace(
                    income=np.ones(3),
                    uncompensated=np.eye(3),
                    compensated=np.eye(3),
                ),
                _title=lambda: "Censored QUAIDS",
            )

        with patch.dict(sys.modules, {"sfi": fake_sfi}), patch(
            "pyquaidsce.estimator.quaidsce", side_effect=fake_quaidsce
        ):
            run_from_stata(
                shares_str="w1 w2 w3",
                prices_str="p1 p2 p3",
                expenditure_str="m",
                demographics_str="z",
                anot=10.0,
                stop_rule="standard",
                bootstrap_start="warm",
                n_jobs=2,
                mp_context="spawn",
                rep_timeout=15.0,
                verbose=False,
            )
        self.assertEqual(captured["stop_rule"], "standard")
        self.assertEqual(captured["bootstrap_start"], "warm")
        self.assertEqual(captured["mp_context"], "spawn")
        self.assertEqual(captured["rep_timeout"], 15.0)
        self.assertTrue(np.array_equal(stored["__pyq_V"], np.eye(2)))

    def test_stata_bridge_forwards_control_function_design(self):
        values = {
            "_touse": [1, 1, 0],
            "w1": [0.2, 0.3, 0.1], "w2": [0.3, 0.2, 0.3],
            "w3": [0.5, 0.5, 0.6],
            "p1": [1.0, 1.1, 1.2], "p2": [1.2, 1.0, 0.9],
            "p3": [0.9, 1.0, 1.1], "m": [10.0, 11.0, 12.0],
            "z": [0.0, 1.0, 2.0], "v": [-0.2, 0.3, 0.1],
            "vsel": [0.4, -0.1, 0.2],
        }

        class Data:
            @staticmethod
            def get(var):
                return values[var]

        class Scalar:
            @staticmethod
            def setValue(*_):
                pass

        class Matrix:
            @staticmethod
            def create(*_):
                pass

            @staticmethod
            def storeAt(*_):
                pass

            @staticmethod
            def setColNames(*_):
                pass

            @staticmethod
            def setRowNames(*_):
                pass

        class Macro:
            @staticmethod
            def setLocal(*_):
                pass

        captured = {}

        def fake_quaidsce(data, **kwargs):
            captured.update(kwargs)
            self.assertEqual(list(data.columns), [
                "w1", "w2", "w3", "p1", "p2", "p3", "z", "m",
                "v", "vsel",
            ])
            self.assertEqual(len(data), 2)
            return types.SimpleNamespace(
                nobs=2, llf=1.0, anot=10.0,
                spec=types.SimpleNamespace(ndemo=1), converged=True,
                n_outer=2, n_gn=3, b=np.array([1.0, 2.0]), V=np.eye(2),
                boot=None, names=["eq1:a", "eq2:b"],
                elas=types.SimpleNamespace(
                    income=np.ones(3), uncompensated=np.eye(3),
                    compensated=np.eye(3),
                ),
                _title=lambda: "Censored QUAIDS",
            )

        fake_sfi = types.SimpleNamespace(
            Data=Data, Scalar=Scalar, Matrix=Matrix, Macro=Macro
        )
        with patch.dict(sys.modules, {"sfi": fake_sfi}), patch(
            "pyquaidsce.estimator.quaidsce", side_effect=fake_quaidsce
        ):
            run_from_stata(
                shares_str="w1 w2 w3", prices_str="p1 p2 p3",
                expenditure_str="m", demographics_str="z", anot=10.0,
                control_function="v",
                selection_control_function="vsel",
                selection_prices_str="p3 p1",
                selection_prices_specified=True,
                selection_covariates_str="",
                selection_covariates_specified=True,
                selection_expenditure=False,
                first_stage_predict="xb", verbose=False,
            )
        self.assertEqual(captured["control_function"], "v")
        self.assertEqual(captured["selection_control_function"], "vsel")
        self.assertEqual(captured["selection_prices"], ["p3", "p1"])
        self.assertEqual(captured["selection_covariates"], [])
        self.assertFalse(captured["selection_expenditure"])

    def test_stata_no_options_use_returned_macro_names(self):
        ado = (ROOT / "stata" / "pyquaidsce.ado").read_text(encoding="utf-8")
        self.assertIn('local is_censor = ("`censor\'" == "")', ado)
        self.assertIn('local is_verbose = ("`log\'" == "")', ado)
        self.assertNotIn('local is_censor = ("`nocensor\'"', ado)
        self.assertNotIn('local is_verbose = ("`nolog\'"', ado)
        help_text = (ROOT / "stata" / "pyquaidsce.sthlp").read_text(
            encoding="utf-8"
        )
        for option in (
            "control_function", "selection_control_function",
            "selection_noprices", "selection_nocovariates",
            "selection_noexpenditure",
        ):
            self.assertIn(option, help_text)

    def test_ado_runs_the_complete_estimation_out_of_process(self):
        ado = (ROOT / "stata" / "pyquaidsce.ado").read_text(encoding="utf-8")
        self.assertIn("launch_from_stata", ado)
        self.assertIn("poll_bootstrap", ado)
        self.assertIn("load_stata_results", ado)
        self.assertNotIn("run_from_stata(", ado)
        self.assertNotIn("reps=0,", ado)
        self.assertEqual(ado.count("launch_from_stata("), 1)
        self.assertIn('local method "ifgnls"', ado)

    def test_bootstrap_runner_module_importable(self):
        """bootstrap_runner.py must be importable and expose main()."""
        from pyquaidsce import bootstrap_runner
        self.assertTrue(callable(bootstrap_runner.main))

    def test_bridge_exposes_async_bootstrap_functions(self):
        """stata_bridge must expose launch/poll/load/kill bootstrap functions."""
        from pyquaidsce.stata_bridge import (
            launch_bootstrap,
            launch_from_stata,
            poll_bootstrap,
            load_bootstrap_results,
            load_stata_results,
            kill_bootstrap,
        )
        for fn in (launch_bootstrap, launch_from_stata, poll_bootstrap,
                   load_bootstrap_results, load_stata_results, kill_bootstrap):
            self.assertTrue(callable(fn))

    def test_launch_from_stata_exports_one_complete_job(self):
        from pyquaidsce.stata_bridge import launch_from_stata

        values = {
            "_touse": [1, 0, 1],
            "w1": [0.2, 0.1, 0.3], "w2": [0.3, 0.3, 0.2],
            "w3": [0.5, 0.6, 0.5],
            "p1": [1.0, 1.2, 1.1], "p2": [1.2, 0.9, 1.0],
            "p3": [0.9, 1.1, 1.0], "m": [10.0, 12.0, 11.0],
        }

        class Data:
            @staticmethod
            def get(var):
                return values[var]

        fake_sfi = types.SimpleNamespace(Data=Data)
        with patch.dict(sys.modules, {"sfi": fake_sfi}), patch(
            "pyquaidsce.stata_bridge._launch_job"
        ) as launch:
            launch_from_stata(
                shares_str="w1 w2 w3", prices_str="p1 p2 p3",
                expenditure_str="m", demographics_str="", anot=10.0,
                reps=7, seed=42, n_jobs=2, stop_rule="standard",
                bootstrap_start="warm", verbose=False,
            )

        launch.assert_called_once()
        frame, kwargs = launch.call_args.args
        self.assertEqual(len(frame), 2)
        self.assertEqual(kwargs["reps"], 7)
        self.assertEqual(kwargs["seed"], 42)
        self.assertEqual(kwargs["stop_rule"], "standard")
        self.assertEqual(kwargs["bootstrap_start"], "warm")

    def test_load_stata_results_restores_all_outputs(self):
        from pyquaidsce.stata_bridge import _BOOT_STATE, load_stata_results

        stored = {}

        class Scalar:
            @staticmethod
            def setValue(name, value):
                stored[name] = value

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

        class Macro:
            @staticmethod
            def setLocal(name, value):
                stored[name] = value

        result = {
            "b": np.array([1.0, 2.0]), "V": np.eye(2),
            "names": ["eq1:a", "eq2:b"], "shares": ["w1", "w2"],
            "elas_i": np.array([1.1, 0.9]), "elas_u": np.eye(2),
            "elas_c": np.eye(2) * 2, "nobs": 10, "llf": -3.0,
            "anot": 10.0, "ndemo": 0, "converged": 1,
            "n_outer": 2, "n_gn": 3, "title": "Censored QUAIDS",
            "boot_reps_ok": 6, "boot_reps_requested": 7,
        }
        handle, path = tempfile.mkstemp(suffix=".pkl")
        os.close(handle)
        try:
            with open(path, "wb") as file_handle:
                pickle.dump(result, file_handle)
            _BOOT_STATE["result_path"] = path
            fake_sfi = types.SimpleNamespace(
                Scalar=Scalar, Matrix=Matrix, Macro=Macro
            )
            with patch.dict(sys.modules, {"sfi": fake_sfi}):
                load_stata_results("b", "V", "ei", "eu", "ec")
            self.assertTrue(np.array_equal(stored["b"], [[1.0, 2.0]]))
            self.assertTrue(np.array_equal(stored["V"], np.eye(2)))
            self.assertEqual(stored["r_boot_reps_ok"], 6)
            self.assertEqual(stored["model_title"], "Censored QUAIDS")
        finally:
            if os.path.exists(path):
                os.unlink(path)
            _BOOT_STATE.clear()

    def test_launch_bootstrap_requires_prior_run(self):
        """launch_bootstrap must raise if no point estimate has been run."""
        from pyquaidsce.stata_bridge import _BOOT_STATE, launch_bootstrap
        saved = dict(_BOOT_STATE)
        _BOOT_STATE.clear()
        try:
            with self.assertRaises(RuntimeError):
                launch_bootstrap(reps=10)
        finally:
            _BOOT_STATE.update(saved)


if __name__ == "__main__":
    unittest.main()
