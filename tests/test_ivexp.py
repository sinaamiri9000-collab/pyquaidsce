"""Tests for the integrated endogenous-expenditure control-function path."""

from pathlib import Path
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from pyquaidsce import fit_expenditure_reduced_form, quaidsce
from pyquaidsce import estimator as estimator_module


ROOT = Path(__file__).resolve().parents[1]


def _small4_with_instrument() -> pd.DataFrame:
    frame = pd.read_stata(ROOT / "bench" / "small4.dta")
    x1 = frame["x1"].to_numpy(float)
    x2 = frame["x2"].to_numpy(float)
    frame["iv_income"] = x1 * x1 + 0.35 * x2 * x2
    return frame


def _common(frame: pd.DataFrame) -> dict:
    return dict(
        data=frame,
        shares=["sw1", "sw2", "sw4", "sw9"],
        prices=["p1", "p2", "p4", "p9"],
        expenditure="total",
        demographics=["x1", "x2"],
        anot=10.0,
        first_stage_predict="xb",
        method="nls",
        start="linear",
        max_iter=160,
        verbose=False,
    )


class ReducedFormUnitTests(unittest.TestCase):
    def test_ols_residual_is_orthogonal_and_joint_f_is_reported(self):
        rng = np.random.default_rng(20260825)
        n = 700
        lnp = rng.normal(size=(n, 3))
        demo = rng.normal(size=(n, 2))
        iv = rng.normal(size=(n, 2))
        error = rng.normal(scale=0.3, size=n)
        y = (
            1.2
            + lnp @ np.array([0.1, -0.2, 0.05])
            + demo @ np.array([0.15, -0.08])
            + iv @ np.array([0.9, -0.65])
            + error
        )
        result = fit_expenditure_reduced_form(
            y,
            lnp,
            demo,
            iv,
            outcome_name="ln(m)",
            price_names=["p1", "p2", "p3"],
            price_inputs_are_logs=False,
            demographic_names=["z1", "z2"],
            instrument_names=["income", "assets"],
        )
        design = np.column_stack([lnp, demo, iv, np.ones(n)])
        np.testing.assert_allclose(
            design.T @ result.residuals,
            np.zeros(design.shape[1]),
            atol=2e-11,
            rtol=0.0,
        )
        self.assertEqual(result.regressor_names[-1], "_cons")
        self.assertEqual(result.excluded_df_num, 2)
        self.assertGreater(result.excluded_f, 100.0)
        self.assertLess(result.excluded_pvalue, 1e-20)
        self.assertEqual(result.V.shape, (8, 8))

    def test_rank_deficiency_fails_loudly(self):
        rng = np.random.default_rng(9)
        n = 80
        lnp = rng.normal(size=(n, 3))
        demo = rng.normal(size=(n, 1))
        iv = np.column_stack([demo[:, 0], rng.normal(size=n)])
        y = rng.normal(size=n)
        with self.assertRaisesRegex(ValueError, "rank deficient"):
            fit_expenditure_reduced_form(
                y,
                lnp,
                demo,
                iv,
                outcome_name="ln(m)",
                price_names=["p1", "p2", "p3"],
                price_inputs_are_logs=True,
                demographic_names=["z"],
                instrument_names=["bad", "good"],
            )


class IntegratedIvexpTests(unittest.TestCase):
    def test_internal_ivexp_matches_the_same_manually_generated_residual(self):
        frame = _small4_with_instrument()
        lnp = np.log(frame[["p1", "p2", "p4", "p9"]].to_numpy(float))
        demo = frame[["x1", "x2"]].to_numpy(float)
        iv = frame[["iv_income"]].to_numpy(float)
        rf = fit_expenditure_reduced_form(
            np.log(frame["total"].to_numpy(float)),
            lnp,
            demo,
            iv,
            outcome_name="ln(total)",
            price_names=["p1", "p2", "p4", "p9"],
            price_inputs_are_logs=False,
            demographic_names=["x1", "x2"],
            instrument_names=["iv_income"],
        )
        frame["manual_cf"] = rf.residuals

        internal = quaidsce(**_common(frame), ivexp=["iv_income"])
        external = quaidsce(
            **_common(frame),
            control_function="manual_cf",
            selection_control_function="manual_cf",
        )
        self.assertTrue(internal.converged)
        self.assertIsNotNone(internal.reduced_form)
        np.testing.assert_allclose(
            internal.reduced_form.residuals, rf.residuals, atol=2e-14, rtol=2e-14
        )
        # Identical designs can stop one Gauss-Newton iteration apart at the
        # requested numerical tolerance; estimates should still agree tightly.
        np.testing.assert_allclose(internal.theta, external.theta, atol=8e-7, rtol=8e-6)
        np.testing.assert_allclose(internal.tau, external.tau, atol=2e-9, rtol=2e-8)
        np.testing.assert_allclose(
            internal.elas.as_stata_vector(),
            external.elas.as_stata_vector(),
            atol=2e-6,
            rtol=2e-5,
        )
        self.assertIn("cf_ivexp", internal.selection_layout.ordered_names)
        self.assertIn("Excluded-IV F", internal.summary())
        self.assertIn("iv_income", internal.reduced_form_table())

        frame["ln_total"] = np.log(frame["total"].to_numpy(float))
        logged_args = _common(frame)
        logged_args.pop("expenditure")
        logged_args["lnexpenditure"] = "ln_total"
        logged = quaidsce(**logged_args, ivexp=["iv_income"])
        np.testing.assert_allclose(internal.theta, logged.theta, atol=8e-7, rtol=8e-6)
        np.testing.assert_allclose(
            internal.reduced_form.residuals,
            logged.reduced_form.residuals,
            atol=2e-14,
            rtol=2e-14,
        )

    def test_invalid_ivexp_contracts_fail_before_estimation(self):
        frame = _small4_with_instrument()
        common = _common(frame)
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            quaidsce(**common, ivexp=["iv_income"], control_function="x1")
        with self.assertRaisesRegex(ValueError, "excluded instruments"):
            quaidsce(**common, ivexp=["p1"])
        with self.assertRaisesRegex(ValueError, "at least one"):
            quaidsce(**common, ivexp=[])
        with self.assertRaisesRegex(ValueError, "sequence"):
            quaidsce(**common, ivexp="iv_income")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            quaidsce(**common, ivexp=["iv_income", "iv_income"])
        with self.assertRaisesRegex(ValueError, "first_stage_predict='xb'"):
            bad = dict(common)
            bad["first_stage_predict"] = "pr"
            quaidsce(**bad, ivexp=["iv_income"])

    def test_bootstrap_reestimates_reduced_form_in_every_replication(self):
        frame = _small4_with_instrument()
        calls = 0
        original = estimator_module.fit_expenditure_reduced_form

        def counted(*args, **kwargs):
            nonlocal calls
            calls += 1
            return original(*args, **kwargs)

        with patch.object(
            estimator_module, "fit_expenditure_reduced_form", side_effect=counted
        ):
            result = quaidsce(
                **_common(frame),
                ivexp=["iv_income"],
                reps=2,
                seed=711,
                n_jobs=1,
            )
        self.assertEqual(calls, 3)
        self.assertIsNotNone(result.boot)
        self.assertEqual(result.boot.reps_ok, 2)
        self.assertTrue(np.isfinite(result.se).all())


if __name__ == "__main__":
    unittest.main()
