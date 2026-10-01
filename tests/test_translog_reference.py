"""Fast external-reference checks plus an optional complete Uruguay refit."""

from __future__ import annotations

import json
import os
import sys
import unittest
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from benchmarks.translog_uruguay.run_python import DEMOS, SHARES, compare, load_data
from pyquaidsce import translog
from pyquaidsce.data import DemandData
from pyquaidsce.translog import _point_at_means
from pyquaidsce.translog_model import TranslogSpec, elasticities, latent_shares


class UruguayReferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.df, cls.prices, cls.cleaned = load_data()
        cls.reference = json.loads((ROOT / "benchmarks/translog_uruguay/reference/stata_reference.json").read_text())
        cls.spec = TranslogSpec(14, 7)
        full = np.array([row["estimate"] for row in cls.reference["coefficients"]])
        cls.theta = np.r_[full[:13], full[14:]]
        sample = cls.df.dropna(subset=SHARES+cls.prices+["gasto_total"]+DEMOS)
        cls.d = DemandData(np.log(sample[cls.prices].to_numpy()), np.log(sample.gasto_total.to_numpy()),
                            sample[SHARES].to_numpy(), sample[DEMOS].to_numpy(),
                            np.ones((len(sample), 14)), np.zeros((len(sample), 14)))

    def test_external_coefficients_reproduce_likelihood_and_r2(self):
        self.assertEqual(self.cleaned, 1808)
        self.assertEqual(self.d.nobs, 6573)
        predicted = latent_shares(self.theta, self.d, self.spec)
        u = self.d.shares[:, :13]-predicted[:, :13]
        _, ld = np.linalg.slogdet(u.T @ u/self.d.nobs)
        llf = -self.d.nobs*13/2*(1+np.log(2*np.pi))-self.d.nobs/2*ld
        self.assertLess(abs(llf-self.reference["llf"]), .01)
        r2 = 1-np.sum((self.d.shares-predicted)**2, axis=0)/np.sum((self.d.shares-self.d.shares.mean(axis=0))**2, axis=0)
        np.testing.assert_allclose(r2, self.reference["r2"], atol=1e-4)

    def test_external_elasticities_and_mean_conventions(self):
        point = _point_at_means(self.df, self.prices, "gasto_total", DEMOS, False, False)
        np.testing.assert_allclose(point.lnp[0], np.log(np.exp(self.d.lnp).mean(axis=0)))
        self.assertAlmostEqual(float(np.exp(point.lnexp[0])), self.df.gasto_total.mean())
        expected = np.zeros((14, 14))
        for row in self.reference["uncompensated"]:
            expected[row["good"]-1, row["price"]-1] = row["estimate"]
        # Directly using printed coefficients has additional rounding error.
        np.testing.assert_allclose(elasticities(self.theta, point, self.spec).uncompensated,
                                    expected, atol=2.5e-4)

    @unittest.skipUnless(os.environ.get("PYQUAIDSCE_RUN_URUGUAY") == "1", "set PYQUAIDSCE_RUN_URUGUAY=1 for the complete refit")
    def test_independent_zero_start_refit(self):
        result = translog(self.df, SHARES, prices=self.prices, expenditure="gasto_total",
                           demographics=DEMOS, verbose=False)
        _, summary, _, _ = compare(result, self.df, self.reference)
        self.assertTrue(summary["passed"], summary)
        self.assertFalse(summary["reference_coefficients_used_as_initial"])


if __name__ == "__main__":
    unittest.main()
