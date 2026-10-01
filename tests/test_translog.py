"""Demand-theory, numerical differentiation, and shared-engine translog checks."""

from __future__ import annotations

import os
import sys
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from pyquaidsce import translog
from pyquaidsce.data import DemandData
from pyquaidsce.inference import complex_step_jacobian
from pyquaidsce.translog_model import (TranslogSpec, delta_matrix, elasticities,
                                      fitted_shares, full_vector, jacobian_free,
                                      latent_shares, share_derivatives, unpack)


def synthetic_sample(seed=21, size=420, demographics=True):
    rng = np.random.default_rng(seed)
    spec = TranslogSpec(3, int(demographics))
    theta = np.array([.30, .25, .02, -.004, .002, -.01, .001, .015]
                      + ([.06, -.03, .04] if demographics else []))
    lnp = rng.normal(0, .55, (size, 3))
    lnm = rng.normal(2.3, .5, size)
    z = rng.uniform(0, 2, (size, spec.ndemo))
    d = DemandData(lnp, lnm, np.zeros((size, 3)), z,
                    np.ones((size, 3)), np.zeros((size, 3)))
    w = latent_shares(theta, d, spec)
    noise = rng.normal(0, .002, (size, 2))
    w[:, :2] += noise
    w[:, 2] -= noise.sum(axis=1)
    frame = pd.DataFrame({**{f"p{i+1}": np.exp(lnp[:, i]) for i in range(3)},
                          **{f"w{i+1}": w[:, i] for i in range(3)},
                          "m": np.exp(lnm)})
    if demographics:
        frame["z"] = z[:, 0]
    d.shares = w.copy()
    return frame, d, spec, theta


class TranslogMathTests(unittest.TestCase):
    def setUp(self):
        self.frame, self.data, self.spec, self.theta = synthetic_sample(size=40)

    def test_parameter_map_and_restrictions(self):
        coefs = unpack(self.theta, self.spec)
        self.assertAlmostEqual(coefs.alpha.sum(), 1)
        np.testing.assert_allclose(coefs.gamma, coefs.gamma.T)
        self.assertGreater(np.max(np.abs(coefs.gamma.sum(axis=1))), .001)
        numerical = complex_step_jacobian(lambda t: full_vector(t, self.spec), self.theta)
        np.testing.assert_allclose(numerical, delta_matrix(self.spec), atol=1e-14)

    def test_analytic_parameter_jacobian_and_slicing(self):
        numeric = []
        for k in range(self.spec.n_free):
            step = np.eye(self.spec.n_free)[k]*1e-6
            numeric.append((fitted_shares(self.theta+step, self.data, self.spec)
                            - fitted_shares(self.theta-step, self.data, self.spec))/2e-6)
        analytic = jacobian_free(self.theta, self.data, self.spec)
        np.testing.assert_allclose(analytic, np.stack(numeric, axis=2), rtol=2e-7, atol=2e-9)
        np.testing.assert_allclose(jacobian_free(self.theta, self.data, self.spec, slice(3, 9)),
                                    analytic[3:9])

    def test_price_and_expenditure_derivatives(self):
        dm, dp = share_derivatives(self.theta, self.data, self.spec)
        h = 1e-6
        plus, minus = self.data.subset(slice(None)), self.data.subset(slice(None))
        plus.lnexp = self.data.lnexp+h
        minus.lnexp = self.data.lnexp-h
        numeric = (latent_shares(self.theta, plus, self.spec)-latent_shares(self.theta, minus, self.spec))/(2*h)
        np.testing.assert_allclose(dm, numeric, atol=2e-9)
        for j in range(3):
            plus.lnp = self.data.lnp.copy()
            minus.lnp = self.data.lnp.copy()
            plus.lnexp = minus.lnexp = self.data.lnexp
            plus.lnp[:, j] += h
            minus.lnp[:, j] -= h
            numeric = (latent_shares(self.theta, plus, self.spec)-latent_shares(self.theta, minus, self.spec))/(2*h)
            np.testing.assert_allclose(dp[:, :, j], numeric, atol=2e-9)

    def test_adding_up_homogeneity_engel_and_slutsky(self):
        w = latent_shares(self.theta, self.data, self.spec)
        dm, dp = share_derivatives(self.theta, self.data, self.spec)
        np.testing.assert_allclose(w.sum(axis=1), 1, atol=1e-14)
        np.testing.assert_allclose(dm+dp.sum(axis=2), 0, atol=1e-14)
        el = elasticities(self.theta, self.data, self.spec)
        np.testing.assert_allclose((w*el.income).sum(axis=1), 1, atol=1e-14)
        np.testing.assert_allclose(el.uncompensated.sum(axis=2), -el.income, atol=1e-14)
        scaled = w[:, :, None]*el.compensated
        np.testing.assert_allclose(scaled, scaled.transpose(0, 2, 1), atol=1e-14)
        np.testing.assert_allclose(scaled.sum(axis=1), 0, atol=1e-14)

    def test_invalid_translation_and_denominator(self):
        bad = self.theta.copy()
        bad[-3:] = 1e6
        with self.assertRaisesRegex(ValueError, "effective expenditure"):
            latent_shares(bad, self.data, self.spec)
        spec = TranslogSpec(2)
        d = DemandData(np.zeros((1, 2)), np.ones(1), np.zeros((1, 2)),
                        np.empty((1, 0)), np.ones((1, 2)), np.zeros((1, 2)))
        with self.assertRaisesRegex(ValueError, "denominator"):
            latent_shares(np.array([.4, 1., 0., 0.]), d, spec)


class TranslogEstimationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame, cls.data, cls.spec, cls.theta = synthetic_sample()
        cls.options = dict(shares=["w1", "w2", "w3"], prices=["p1", "p2", "p3"],
                           expenditure="m", demographics=["z"], verbose=False)
        cls.result = translog(cls.frame, **cls.options)

    def test_independent_fit_and_results(self):
        r = self.result
        self.assertTrue(r.converged)
        # Sampling noise permits a tolerance, while a misspecified model does not.
        np.testing.assert_allclose(r.predict(), latent_shares(self.theta, self.data, self.spec), atol=.0012)
        self.assertEqual(r.predict().shape, (len(self.frame), 3))
        np.testing.assert_allclose(r.predict().sum(axis=1), 1, atol=1e-12)
        np.testing.assert_allclose(r.predict(kind="residuals"), self.data.shares-r.predict())
        np.testing.assert_allclose(r.predict(kind="quantities")*self.frame[["p1", "p2", "p3"]],
                                    r.predict()*self.frame.m.to_numpy()[:, None])
        self.assertTrue(np.isfinite(r.se).all())
        self.assertIn("Basic translog", r.summary())
        self.assertIn("Marshallian", r.elasticity_tables())
        self.assertEqual(r.get("alpha_1"), r.coefs.alpha[0])

    def test_levels_logs_missing_rows_and_evaluation(self):
        frame = self.frame.copy()
        frame[["p1", "p2", "p3"]] = np.log(frame[["p1", "p2", "p3"]])
        frame["m"] = np.log(frame.m)
        r = translog(frame, ["w1", "w2", "w3"], lnprices=["p1", "p2", "p3"],
                      lnexpenditure="m", demographics=["z"], initial=self.result.theta,
                      sigma_initial=self.result.sigma, verbose=False)
        np.testing.assert_allclose(r.predict(), self.result.predict(), atol=1e-5)
        supplied = self.frame.copy()
        supplied.loc[3, "p2"] = np.nan
        pred = self.result.predict(supplied)
        self.assertTrue(np.isnan(pred[3]).all())
        self.assertTrue(np.isfinite(pred[np.arange(len(pred)) != 3]).all())
        el, se = self.result.elasticities_at_means(standard_errors=True)
        np.testing.assert_allclose(el.uncompensated, self.result.elas.uncompensated)
        np.testing.assert_allclose(se.uncompensated, self.result.elasticity_se.uncompensated, atol=1e-10)

    def test_no_demographics_and_other_sur_methods(self):
        frame, _, _, theta = synthetic_sample(demographics=False)
        for method in ("nls", "fgnls", "ifgnls"):
            r = translog(frame, ["w1", "w2", "w3"], prices=["p1", "p2", "p3"],
                          expenditure="m", method=method, initial=theta, verbose=False)
            self.assertTrue(r.converged)
            self.assertEqual(r.spec.ndemo, 0)
            self.assertTrue(np.isfinite(r.b).all())

    def test_validation_does_not_clean_or_enable_unimplemented_features(self):
        frame = self.frame.copy()
        frame.loc[0, "w1"] = -1e-9
        with self.assertRaisesRegex(ValueError, "clean negative"):
            translog(frame, **self.options)
        with self.assertRaisesRegex(ValueError, "no censoring"):
            translog(self.frame, **self.options, censor=True)
        with self.assertRaises(TypeError):
            translog(self.frame, **self.options, ivexp=["z"])
        with self.assertRaisesRegex(ValueError, "initial must contain"):
            translog(self.frame, **self.options, initial=[0])

    def test_shared_bootstrap_is_deterministic_across_spawn_and_serial(self):
        options = dict(self.options, reps=3, seed=81, bootstrap_start="warm",
                        initial=self.result.theta, sigma_initial=self.result.sigma)
        serial = translog(self.frame, **options, n_jobs=1)
        spawned = translog(self.frame, **options, n_jobs=2, mp_context="spawn")
        self.assertEqual(serial.boot.reps_ok, 3)
        self.assertEqual(spawned.boot.reps_ok, 3)
        np.testing.assert_allclose(serial.boot.b_star, spawned.boot.b_star, atol=1e-12)
        np.testing.assert_allclose(serial.V, serial.boot.V)
        np.testing.assert_allclose(serial.analytic_se, self.result.analytic_se, atol=2e-6)
        self.assertTrue(np.isfinite(serial.elasticity_se.uncompensated).all())


if __name__ == "__main__":
    unittest.main()
