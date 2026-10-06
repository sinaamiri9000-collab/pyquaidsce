"""Independent mean, support and correlation checks for the second DGP."""

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
from scipy.stats import norm

from pyquaidsce import DemandData, fitted_shares, quaidsce
from pyquaidsce.inference import _EstimatingSystem, _central_jacobian
from pyquaidsce.model import jacobian_full
from pyquaidsce.params import delta_matrix
from tools.simulate_separate_participation import (
    ANOT, LOG_PRICES, OPTIONS, SHARES, SPEC, generate, indices, mean_model, realize, regressors,
    support_bounds, true_theta,
)


class SeparateParticipationDGPTests(unittest.TestCase):
    def test_separate_purchases_positive_shares_and_complete_budget(self):
        positive, outside = support_bounds()
        self.assertTrue(np.all(positive > 0))
        self.assertGreater(outside, 0)
        frame = generate(2000, 27043)
        shares = frame[SHARES].to_numpy()
        purchased = shares > 0
        self.assertTrue(np.all(purchased.any(axis=0)))
        self.assertTrue(np.all((~purchased).any(axis=0)))
        self.assertGreater(np.mean(np.any(purchased != purchased[:, :1], axis=1)), .3)
        self.assertTrue(np.all(shares[purchased] > np.broadcast_to(positive, shares.shape)[purchased]))
        self.assertTrue(np.all(frame.outside_share > outside))
        np.testing.assert_allclose(shares.sum(axis=1)+frame.outside_share, 1, atol=2e-15)

    def test_conditional_mean_matches_existing_package_equations(self):
        reg = regressors(np.random.default_rng(27043).random((2000, 7)))
        xb = indices(reg)
        data = DemandData(reg[:, :4], reg[:, 4], np.zeros((len(reg), 4)),
                          reg[:, 5, None], norm.cdf(xb), norm.pdf(xb), ANOT)
        np.testing.assert_allclose(fitted_shares(true_theta(), data, SPEC),
                                   mean_model(reg), rtol=2e-14, atol=2e-15)

    def test_simulated_conditional_mean_and_correlated_purchase_events(self):
        # Hold X fixed: correlation changes joint purchase events but not
        # marginal Probits or the exact conditional demand mean.
        row = np.array([[.2, -.1, .15, -.3, 2.2, .35, -.25]])
        reg = np.repeat(row, 120000, axis=0)
        expected = mean_model(row)[0]
        event_correlations = []
        for correlation in [0., .5]:
            shares, _ = realize(reg, np.random.default_rng(27044), correlation)
            observed = shares.mean(axis=0)
            mcse = shares.std(axis=0, ddof=1)/np.sqrt(len(reg))
            self.assertTrue(np.all(np.abs(observed-expected) < 5*mcse))
            events = (shares > 0).astype(float)
            event_correlations.append(np.corrcoef(events, rowvar=False)[0, 1])
        self.assertLess(abs(event_correlations[0]), .02)
        self.assertGreater(event_correlations[1], .25)

    def test_complete_bread_with_distinct_probits_and_excluded_covariate(self):
        frame = generate(2000, 2171542997, .5)
        prices = frame[LOG_PRICES].to_numpy()
        expenditure = frame.lm.to_numpy()
        design = np.column_stack([prices, expenditure, frame[['z', 's']].to_numpy()])
        initial = quaidsce(frame, **OPTIONS, method='nls').theta
        for method in ['fgnls', 'ifgnls']:
            with self.subTest(method=method):
                fit = quaidsce(frame, **OPTIONS, method=method)
                xb = np.column_stack([design, np.ones(len(frame))]) @ fit.tau.reshape(4, -1).T
                data = DemandData(prices, expenditure, frame[SHARES].to_numpy(),
                                  frame[['z']].to_numpy(), norm.cdf(xb), norm.pdf(xb), ANOT)
                system = _EstimatingSystem(fit, data, design,
                                           initial if method == 'fgnls' else None, 2000, None)
                step = np.cbrt(np.finfo(float).eps)/4
                small = _central_jacobian(system.mean_score, system.x, system.scales, step)
                large = _central_jacobian(system.mean_score, system.x, system.scales, 2*step)
                reference = (4*small-large)/3
                _, theta, sigma = system.decode(system.x)
                jac = jacobian_full(theta, data, fit.spec) @ delta_matrix(fit.spec)
                u = data.shares-fitted_shares(theta, data, fit.spec)
                # Complex perturbations avoid subtraction of almost equal
                # demand scores in the covariance-weight derivative block.
                for column, (a, b) in enumerate(zip(*system.tril)):
                    h = 1e-20*system.scales[system.sigma_slice.start+column]
                    perturbed = sigma.astype(complex)
                    perturbed[a, b] += 1j*h
                    if a != b:
                        perturbed[b, a] += 1j*h
                    score = np.einsum('tmi,tm->i', jac, u @ np.linalg.inv(perturbed))
                    reference[system.theta_slice, system.sigma_slice.start+column] = (
                        score.imag/h/data.nobs)
                analytical = system.bread()
                error = np.max(np.abs(analytical-reference)/np.maximum(1, np.abs(analytical)))
                self.assertLess(float(error), 5e-6)


if __name__ == '__main__':
    unittest.main()
