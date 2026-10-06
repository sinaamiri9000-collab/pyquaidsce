"""Check the Monte Carlo DGP against the model it is meant to validate."""

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
from scipy.stats import norm

from pyquaidsce import DemandData, fitted_shares
from tools.simulate_analytical_quaids import (ANOT, SHARES, SPEC, TAU, generate,
                                             mean_model, positivity_bound,
                                             regressors, true_theta)


class QuaidsDGPTests(unittest.TestCase):
    def test_positive_purchased_shares_and_complete_budget(self):
        self.assertTrue(np.all(positivity_bound() > 0))
        frame = generate(2000, 17043)
        shares = frame[SHARES].to_numpy()
        participating = frame.outside_share.to_numpy() == 0
        self.assertTrue(participating.any() and (~participating).any())
        self.assertTrue(np.all(shares[participating] > 0))
        self.assertTrue(np.all(shares[participating] < 1))
        np.testing.assert_array_equal(shares[~participating], 0)
        np.testing.assert_allclose(shares.sum(axis=1)+frame.outside_share,
                                   1, atol=2e-15)

    def test_exact_conditional_mean_matches_package_model(self):
        unit = np.random.default_rng(17043).random((2000, 6))
        reg = regressors(unit)
        index = np.column_stack([reg, np.ones(len(reg))]) @ TAU
        data = DemandData(reg[:, :4], reg[:, 4], np.zeros((len(reg), 4)),
                          reg[:, 5, None], np.repeat(norm.cdf(index)[:, None], 4, axis=1),
                          np.repeat(norm.pdf(index)[:, None], 4, axis=1), ANOT)
        np.testing.assert_allclose(fitted_shares(true_theta(), data, SPEC),
                                   mean_model(reg), rtol=2e-14, atol=2e-15)


if __name__ == '__main__':
    unittest.main()
