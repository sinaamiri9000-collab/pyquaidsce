"""Explicit stopping controls and covariance for identity-weighted NLS."""

import importlib
import inspect
import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce
from pyquaidsce.model import jacobian_full, residuals
from pyquaidsce.params import Spec, delta_matrix
from tests.test_internals import _fake_data
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES

solver = importlib.import_module("pyquaidsce.nlsur")


class ConvergenceTests(unittest.TestCase):
    def test_all_numerical_entrypoints_expose_the_same_defaults(self):
        from pyquaidsce.stata_bridge import launch_from_stata, run_from_stata

        expected = dict(param_tol=1e-5, objective_tol=1e-7, gn_tol=1e-5)
        for function in (quaidsce, solver.nlsur, solver.gauss_newton,
                         run_from_stata, launch_from_stata):
            params = inspect.signature(function).parameters
            for name, value in expected.items():
                self.assertEqual(params[name].default, value)
            self.assertEqual(params["max_iter"].default, 300)
            for name in ("stop_rule", "tol", "nrtol_stop",
                         "inner_nrtol_early", "sigma_tol", "boot_sigma_tol",
                         "boot_outer_param_tol"):
                self.assertNotIn(name, params)
            if function is not solver.gauss_newton:
                self.assertEqual(params["outer_param_tol"].default, 1e-5)

    def test_each_inner_threshold_can_independently_stop_a_fit(self):
        # A scalar nonlinear least-squares problem: f(theta)=theta^2, y=4.
        # The normal equations and objective are computed independently here.
        def normal(theta, *args, **kwargs):
            x = theta[0]
            u, j = 4 - x*x, 2*x
            return np.array([[j*j]]), np.array([j*u]), float(u*u)

        def objective(theta, *args, **kwargs):
            return float((4 - theta[0]**2)**2)

        strict = dict(param_tol=1e-14, objective_tol=1e-14, gn_tol=1e-14)
        with patch.object(solver, "make_cache", return_value=None), patch.object(
            solver, "_normal_equations", side_effect=normal
        ), patch.object(solver, "_safe_obj", side_effect=objective):
            _, _, count, ok = solver.gauss_newton(
                np.array([1.0]), None, None, np.eye(1), **strict)
            self.assertTrue(ok)
            self.assertGreater(count, 1)
            for name, loose in dict(param_tol=0.8, objective_tol=0.5,
                                    gn_tol=1.1).items():
                with self.subTest(threshold=name):
                    kwargs = dict(strict, **{name: loose})
                    theta, _, count, ok = solver.gauss_newton(
                        np.array([1.0]), None, None, np.eye(1), **kwargs)
                    self.assertTrue(ok)
                    self.assertEqual(count, 1)
                    np.testing.assert_array_equal(theta, [2.5])

    def test_no_improving_step_uses_only_the_gn_threshold(self):
        normal = (np.eye(1), np.array([0.01]), 1.0)
        with patch.object(solver, "make_cache", return_value=None), patch.object(
            solver, "_normal_equations", return_value=normal
        ), patch.object(solver, "_safe_obj", return_value=1.0):
            for algorithm in ("gn", "lm"):
                for gn_tol, expected in ((1e-3, True), (1e-5, False)):
                    with self.subTest(algorithm=algorithm, gn_tol=gn_tol):
                        theta, _, _, ok = solver.gauss_newton(
                            np.array([0.0]), None, None, np.eye(1),
                            param_tol=10, objective_tol=10, gn_tol=gn_tol,
                            algorithm=algorithm)
                        self.assertEqual(ok, expected)
                        np.testing.assert_array_equal(theta, [0.0])

    def test_outer_requires_two_consecutive_small_changes_and_inner_success(self):
        # A large change between two small ones must reset the confirmation.
        for warm in (False, True):
            for final_inner_ok in (False, True):
                with self.subTest(warm=warm, final_inner_ok=final_inner_ok):
                    changes = iter((1e-3, 0.5, 1e-3, 1e-3))
                    calls = []

                    def inner(theta, *args, **kwargs):
                        calls.append(kwargs)
                        if not warm and len(calls) <= 2:
                            return theta.copy(), 1.0, 1, True
                        change = next(changes)
                        new = theta * (1 + change)
                        last = len(calls) == (4 if warm else 6)
                        return new, 1.0, 1, final_inner_ok if last else True

                    data = type("Data", (), {"nobs": 4})()
                    spec = type("Spec", (), {"n_free": 1, "n_eq_estimated": 1})()
                    with patch.object(solver, "gauss_newton", side_effect=inner), \
                            patch.object(solver, "residuals", return_value=np.ones((4, 1))), \
                            patch.object(solver, "_objective", return_value=1.0), \
                            patch.object(solver, "_finish", side_effect=lambda *a, **k: a):
                        result = solver.nlsur(
                            data, spec, theta0=np.ones(1),
                            sigma0=np.eye(1) if warm else None,
                            outer_param_tol=0.01, max_outer=20)
                    self.assertEqual(result[8], 6)  # n_outer
                    self.assertEqual(result[10], final_inner_ok)
                    self.assertTrue(all(c["gn_tol"] == 1e-5 for c in calls))

    def test_invalid_thresholds_fail_before_estimation(self):
        for name in ("param_tol", "objective_tol", "gn_tol",
                     "outer_param_tol"):
            for value in (0, -1, np.nan, np.inf, -np.inf):
                with self.subTest(name=name, value=value):
                    with self.assertRaisesRegex(ValueError, name):
                        quaidsce(None, shares=[], anot=10, **{name: value})
                    with self.assertRaisesRegex(ValueError, name):
                        solver.nlsur(None, None, **{name: value})

    def test_bootstrap_uses_exactly_the_point_estimate_thresholds(self):
        data = pd.read_stata(f"{BENCH}/small4.dta")
        with patch("pyquaidsce.bootstrap.bootstrap") as boot:
            quaidsce(data, shares=SHARES, prices=PRICES, expenditure="total",
                     demographics=DEMOS, anot=10, method="nls", reps=2,
                     param_tol=2e-5, objective_tol=3e-7, gn_tol=4e-5,
                     outer_param_tol=5e-5,
                     verbose=False)
        kwargs = boot.call_args.kwargs
        for name, value in dict(param_tol=2e-5, objective_tol=3e-7,
                                gn_tol=4e-5, outer_param_tol=5e-5).items():
            self.assertEqual(kwargs[name], value)


class NlsCovarianceTests(unittest.TestCase):
    def test_nls_sandwich_matches_independent_full_jacobian_reference(self):
        rng = np.random.default_rng(170)
        for censor in (False, True):
            with self.subTest(censor=censor):
                spec = Spec(3, 1, True, censor)
                data = _fake_data(rng, 180, 3, 1, censored=censor)
                theta = rng.normal(0, 0.05, spec.n_free)
                u = residuals(theta, data, spec)
                sigma = u.T @ u / data.nobs
                J = jacobian_full(theta, data, spec) @ delta_matrix(spec)
                A = np.einsum("nik,nil->kl", J, J)
                B = np.einsum("nik,ij,njl->kl", J, sigma, J)
                expected = np.linalg.solve(A, B) @ np.linalg.inv(A)
                sigma_from = lambda th: sigma
                for vce in ("objective", "final"):
                    result = solver._finish(
                        theta, data, spec, np.eye(spec.n_eq_estimated), sigma_from,
                        "nls", vce, 37, 1, 1, True, 1.0, [1.0])
                    np.testing.assert_allclose(result.V, expected, rtol=1e-8, atol=1e-10)
                    np.testing.assert_array_equal(result.theta, theta)
                gls = np.linalg.inv(np.einsum("nik,ij,njl->kl", J, np.linalg.inv(sigma), J))
                self.assertGreater(np.linalg.norm(expected-gls), 1e-6)

    def test_fgnls_and_ifgnls_covariances_keep_the_gls_formula(self):
        rng = np.random.default_rng(171)
        spec = Spec(3, 1, True, True)
        data = _fake_data(rng, 180, 3, 1)
        theta = rng.normal(0, 0.05, spec.n_free)
        u = residuals(theta, data, spec)
        sigma_final = u.T @ u / data.nobs
        sigma_obj = np.diag([0.1, 0.2, 0.3])
        J = jacobian_full(theta, data, spec) @ delta_matrix(spec)
        for method in ("fgnls", "ifgnls"):
            for vce, sigma in (("objective", sigma_obj), ("final", sigma_final)):
                with self.subTest(method=method, vce=vce):
                    expected = np.linalg.inv(np.einsum(
                        "nik,ij,njl->kl", J, np.linalg.inv(sigma), J))
                    result = solver._finish(
                        theta, data, spec, sigma_obj, lambda th: sigma_final,
                        method, vce, 37, 2, 2, True, 1.0, [1.0])
                    np.testing.assert_allclose(result.V, expected, rtol=1e-8, atol=1e-10)
