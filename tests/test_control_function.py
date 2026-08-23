"""Acceptance tests for the external control-function extension."""

from __future__ import annotations

import sys
import unittest
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

import numpy as np
from scipy.stats import norm

from pyquaidsce.elasticities import (Means, elasticities,
                                     fitted_share_derivatives, sample_means)
from pyquaidsce.estimator import first_stage, quaidsce
from pyquaidsce.jacfree import jacobian_free, make_cache
from pyquaidsce.model import (DemandData, augmented_latent_shares,
                              fitted_shares, jacobian_full, latent_shares)
from pyquaidsce.params import (Coefs, Spec, delta_matrix, free_slices,
                              full_slices, full_vector, unpack)
from pyquaidsce.selection import FirstStageLayout, legacy_layout


def _data(rng, n=4, R=2, N=18, residual=True):
    lnp = rng.normal(0.0, 0.18, (N, n))
    lnexp = rng.normal(2.0, 0.25, N)
    demo = rng.normal(0.2, 0.15, (N, R))
    shares = rng.uniform(0.05, 0.5, (N, n))
    shares /= shares.sum(axis=1, keepdims=True)
    k = rng.normal(0.0, 0.5, (N, n))
    v = rng.normal(0.0, 0.4, N) if residual else np.zeros(N)
    return DemandData(
        lnp=lnp,
        lnexp=lnexp,
        shares=shares,
        demo=demo,
        cdf=norm.cdf(k),
        pdf=norm.pdf(k),
        a0=1.6,
        control_function=v,
    ), k


def _lift_old_theta(theta, old, new):
    out = np.zeros(new.n_free)
    old_sl, new_sl = free_slices(old), free_slices(new)
    for key in old_sl:
        out[new_sl[key]] = theta[old_sl[key]]
    return out


def _synthetic_cf_dgp(N, seed):
    """Nested sample from a correctly specified CF/selection outcome DGP."""
    import pandas as pd

    rng = np.random.default_rng(seed)
    n = 3
    lnp = rng.normal(0.0, 0.18, (N, n))
    z = rng.normal(0.0, 0.45, N)
    instrument = rng.normal(size=N)
    residual = rng.normal(0.0, 0.30, N)
    # v is the innovation from a correctly specified expenditure reduced form
    # containing prices, z, and an excluded instrument.
    lnexp = (
        2.0 + lnp @ np.array([0.12, -0.08, 0.06])
        + 0.18 * z + 0.45 * instrument + residual
    )
    spec = Spec(n, 1, True, True, True)
    theta = np.zeros(spec.n_free)
    sl = free_slices(spec)
    theta[sl["alpha"]] = np.array([0.31, 0.27, 0.34])
    theta[sl["beta"]] = np.array([0.025, -0.018, 0.012])
    theta[sl["gamma"]] = np.array([0.008, -0.003, 0.006])
    theta[sl["lambda"]] = np.array([0.0015, -0.0010, 0.0005])
    delta = np.array([0.012, 0.010, 0.014])
    kappa = np.array([0.12, -0.10, 0.08])
    theta[sl["delta"]] = delta
    theta[sl["cfcoef"]] = kappa
    theta[sl["eta"]] = np.array([0.008, -0.005])
    theta[sl["rho"]] = np.array([0.025])
    d = DemandData(
        lnp=lnp, lnexp=lnexp, shares=np.zeros((N, n)), demo=z[:, None],
        cdf=np.ones((N, n)), pdf=np.zeros((N, n)), a0=1.6,
        control_function=residual,
    )
    wstar = latent_shares(theta, d, spec)

    tau_p = np.array([
        [0.20, -0.10, 0.05], [-0.10, 0.15, 0.08],
        [0.05, -0.06, 0.12],
    ])
    k = (
        lnp @ tau_p.T
        + lnexp[:, None] * np.array([0.18, 0.12, 0.10])
        + z[:, None] * np.array([0.10, -0.08, 0.05])
        + residual[:, None] * np.array([0.25, -0.20, 0.15])
        + np.array([-0.45, -0.25, -0.05])
    )
    selection_error = rng.normal(size=(N, n))
    selected = k + selection_error > 0
    outcome_noise = rng.normal(0.0, 0.012, (N, n))
    shares = selected * (
        wstar + residual[:, None] * kappa
        + selection_error * delta + outcome_noise
    )
    if np.any(shares < 0):
        raise AssertionError("synthetic DGP calibration produced a negative share")
    frame = pd.DataFrame({
        **{f"w{i + 1}": shares[:, i] for i in range(n)},
        **{f"p{i + 1}": np.exp(lnp[:, i]) for i in range(n)},
        "m": np.exp(lnexp), "z": z, "instrument": instrument, "v": residual,
    })
    return frame, kappa


class ParameterMapAndJacobianTests(unittest.TestCase):
    def test_public_v101_positional_apis_remain_valid(self):
        n = 3
        old_style = Coefs(
            np.full(n, 0.2), np.zeros(n), np.zeros((n, n)),
            np.zeros(n), np.zeros(n), np.zeros((1, n)), np.zeros(1),
        )
        np.testing.assert_array_equal(old_style.cfcoef, np.zeros(n))

        spec = Spec(n, 1, True, True, False)
        c = unpack(np.zeros(spec.n_free), spec)
        means = Means(
            w=np.full(n, 0.25), lnp=np.zeros(n), lnexp=1.0,
            demo=np.zeros(1), cdf=np.full(n, 0.6), pdf=np.full(n, 0.3),
            du=np.zeros(n),
        )
        np_prob = n + 1 + 1 + 1
        tau = np.zeros(n * np_prob)
        positional = elasticities(c, spec, means, 1.6, tau, np_prob, False)
        keyword = elasticities(
            c, spec, means, 1.6, tau=tau, np_prob=np_prob,
            strict_stata=False,
        )
        np.testing.assert_array_equal(
            positional.as_stata_vector(), keyword.as_stata_vector()
        )

    def test_cfcoef_map_dimensions_and_identity(self):
        rng = np.random.default_rng(20260812)
        for n in (3, 4, 10):
            spec = Spec(n, 2, True, True, True)
            theta = rng.normal(0.0, 0.05, spec.n_free)
            c = unpack(theta, spec)
            xs, fs = free_slices(spec), full_slices(spec)
            self.assertEqual(xs["cfcoef"].stop - xs["cfcoef"].start, n)
            self.assertEqual(fs["cfcoef"].stop - fs["cfcoef"].start, n)
            np.testing.assert_array_equal(c.cfcoef, theta[xs["cfcoef"]])
            np.testing.assert_array_equal(
                full_vector(theta, spec)[fs["cfcoef"]], c.cfcoef
            )
            np.testing.assert_array_equal(
                delta_matrix(spec)[fs["cfcoef"], xs["cfcoef"]], np.eye(n)
            )

    def test_active_full_and_fast_jacobians_match_finite_difference(self):
        rng = np.random.default_rng(431)
        for n in (3, 4, 10):
            for residual in (False, True):
                spec = Spec(n, 2, True, True, True)
                d, _ = _data(rng, n=n, R=2, N=7, residual=residual)
                theta = rng.normal(0.0, 0.035, spec.n_free)
                analytic_full = jacobian_full(theta, d, spec)
                analytic = (
                    analytic_full.reshape(-1, spec.n_full) @ delta_matrix(spec)
                )
                fast = jacobian_free(theta, d, spec, make_cache(spec)).reshape(
                    -1, spec.n_free
                )
                np.testing.assert_allclose(fast, analytic, atol=2e-13, rtol=2e-13)
                for h in (1e-4, 1e-5, 1e-6):
                    numeric = np.empty_like(analytic)
                    for q in range(spec.n_free):
                        tp, tm = theta.copy(), theta.copy()
                        tp[q] += h
                        tm[q] -= h
                        numeric[:, q] = (
                            fitted_shares(tp, d, spec).ravel()
                            - fitted_shares(tm, d, spec).ravel()
                        ) / (2.0 * h)
                    scale = max(1.0, float(np.max(np.abs(numeric))))
                    self.assertLess(
                        float(np.max(np.abs(analytic - numeric))) / scale,
                        2e-7,
                    )
                csl = free_slices(spec)["cfcoef"]
                for i in range(n):
                    col = fast.reshape(d.nobs, n, -1)[:, :, csl.start + i]
                    expected = np.zeros_like(col)
                    expected[:, i] = d.cdf[:, i] * d.control_function
                    np.testing.assert_allclose(col, expected, atol=0.0, rtol=0.0)

    def test_cfcoef_zero_preserves_old_blocks_and_elasticities(self):
        rng = np.random.default_rng(901)
        n, R = 4, 2
        old = Spec(n, R, True, True, False)
        new = Spec(n, R, True, True, True)
        d, k = _data(rng, n=n, R=R, N=21, residual=True)
        theta_old = rng.normal(0.0, 0.03, old.n_free)
        theta_new = _lift_old_theta(theta_old, old, new)
        np.testing.assert_allclose(
            fitted_shares(theta_old, d, old),
            fitted_shares(theta_new, d, new),
            atol=0.0,
            rtol=0.0,
        )
        J_old = jacobian_full(theta_old, d, old)
        J_new = jacobian_full(theta_new, d, new)
        for key, old_slice in full_slices(old).items():
            np.testing.assert_allclose(
                J_old[:, :, old_slice],
                J_new[:, :, full_slices(new)[key]],
                atol=0.0,
                rtol=0.0,
            )

        layout = legacy_layout(
            [f"price{j + 1}" for j in range(n)], ["z1", "z2"],
            include_expenditure=True,
        )
        tau = rng.normal(0.0, 0.1, n * layout.width)
        means = sample_means(d, k, old)
        e_old = elasticities(unpack(theta_old, old), old, means, 1.6,
                             tau=tau, np_prob=layout.width, layout=layout,
                             strict_stata=False)
        e_new = elasticities(unpack(theta_new, new), new, means, 1.6,
                             tau=tau, np_prob=layout.width, layout=layout,
                             strict_stata=False)
        np.testing.assert_allclose(e_old.as_stata_vector(), e_new.as_stata_vector(),
                                   atol=0.0, rtol=0.0)


class SelectionLayoutTests(unittest.TestCase):
    def test_layout_supports_subsets_reordering_and_intercept_only(self):
        layout = FirstStageLayout(
            ordered_names=("p3", "p1", "M", "age", "cfunc"),
            demand_price_names=("a", "b", "c"),
            price_positions={"c": 0, "a": 1},
            expenditure_position=2,
            covariate_positions={"age": 3},
            selection_cf_position=4,
            constant_position=5,
        )
        self.assertEqual(layout.price_position(0), 1)
        self.assertIsNone(layout.price_position(1))
        self.assertEqual(layout.price_position(2), 0)
        self.assertEqual(layout.width, 6)
        self.assertEqual(layout.tau_names(1)[-1], "tau:cons_1")

        collision_safe = FirstStageLayout(
            ordered_names=("p1", "z[p1]", "cf[cons]"),
            demand_price_names=("price1", "price2", "price3"),
            price_positions={"price1": 0},
            expenditure_position=None,
            covariate_positions={"p1": 1},
            selection_cf_position=2,
            constant_position=3,
        )
        names = collision_safe.tau_names(1)
        self.assertEqual(len(names), len(set(names)))

        rng = np.random.default_rng(4)
        shares = rng.uniform(0.0, 1.0, (80, 3))
        shares[:20, :] = 0.0
        intercept = FirstStageLayout((), ("a", "b", "c"), {}, None, {}, None, 0)
        fs = first_stage(
            shares, np.zeros((80, 3)), np.ones(80), np.zeros((80, 0)),
            predict="xb", design=np.zeros((80, 0)), layout=intercept,
        )
        self.assertEqual(fs.np_prob, 1)
        self.assertEqual(fs.tau.size, 3)


class PriceExpenditureDerivativeTests(unittest.TestCase):
    def _case(self, residual):
        rng = np.random.default_rng(2026 if residual else 2027)
        n, R, N = 10, 2, 7
        d0, _ = _data(rng, n=n, R=R, N=N, residual=residual)
        spec = Spec(n, R, True, True, True)
        theta = rng.normal(0.0, 0.025, spec.n_free)
        names = tuple(f"p{j + 1}" for j in range(n))
        ordered = tuple([*names, "M", "z1", "cfunc"])
        layout = FirstStageLayout(
            ordered_names=ordered,
            demand_price_names=names,
            price_positions={name: j for j, name in enumerate(names)},
            expenditure_position=n,
            covariate_positions={"z1": n + 1},
            selection_cf_position=n + 2,
            constant_position=n + 3,
        )
        tau = rng.normal(0.0, 0.12, n * layout.width)

        def predict(lnp, lnexp):
            X = np.column_stack(
                [lnp, lnexp, d0.demo[:, 0], d0.control_function]
            )
            k = np.column_stack([
                X @ tau[layout.equation_slice(i)][:-1]
                + tau[layout.equation_slice(i)][-1]
                for i in range(n)
            ])
            d = DemandData(
                lnp=lnp, lnexp=lnexp, shares=d0.shares, demo=d0.demo,
                cdf=norm.cdf(k), pdf=norm.pdf(k), a0=d0.a0,
                control_function=d0.control_function,
            )
            return fitted_shares(theta, d, spec), d, k

        base, d, k = predict(d0.lnp, d0.lnexp)
        self.assertEqual(base.shape, (N, n))
        d_exp, d_price = fitted_share_derivatives(
            theta, d, spec, tau=tau, layout=layout, selection_index=k
        )
        for h in (1e-4, 1e-5, 1e-6):
            p_exp, _, _ = predict(d0.lnp, d0.lnexp + h)
            m_exp, _, _ = predict(d0.lnp, d0.lnexp - h)
            numeric_exp = (p_exp - m_exp) / (2.0 * h)
            self.assertLess(float(np.max(np.abs(d_exp - numeric_exp))), 2e-7)
            for j in range(n):
                lp, lm = d0.lnp.copy(), d0.lnp.copy()
                lp[:, j] += h
                lm[:, j] -= h
                p_price, _, _ = predict(lp, d0.lnexp)
                m_price, _, _ = predict(lm, d0.lnexp)
                numeric = (p_price - m_price) / (2.0 * h)
                self.assertLess(
                    float(np.max(np.abs(d_price[:, :, j] - numeric))),
                    2e-7,
                )

    def test_all_100_price_and_10_expenditure_derivatives_nonzero_residual(self):
        self._case(True)

    def test_all_100_price_and_10_expenditure_derivatives_zero_residual(self):
        self._case(False)

    def test_noquadratic_cf_at_means_matches_direct_expenditure_derivative(self):
        rng = np.random.default_rng(9921)
        n, R = 4, 2
        spec = Spec(n, R, False, True, True)
        d0, k = _data(rng, n=n, R=R, N=1, residual=True)
        theta = rng.normal(0.0, 0.01, spec.n_free)
        sl = free_slices(spec)
        theta[sl["alpha"]] = np.array([0.25, 0.22, 0.28, 0.24])
        theta[sl["delta"]] = np.array([0.02, -0.01, 0.015, 0.005])
        theta[sl["cfcoef"]] = np.array([0.12, -0.08, 0.05, -0.03])
        wstar = latent_shares(theta, d0, spec)
        d = DemandData(
            lnp=d0.lnp, lnexp=d0.lnexp, shares=wstar, demo=d0.demo,
            cdf=d0.cdf, pdf=d0.pdf, a0=d0.a0,
            control_function=d0.control_function,
        )
        layout = legacy_layout(
            [f"p{j + 1}" for j in range(n)], ["z1", "z2"],
            include_expenditure=True,
        )
        tau = rng.normal(0.0, 0.05, n * layout.width)
        derivative, _ = fitted_share_derivatives(
            theta, d, spec, tau=tau, layout=layout, selection_index=k
        )
        e = elasticities(
            unpack(theta, spec), spec, sample_means(d, k, spec), d.a0,
            tau=tau, np_prob=layout.width, strict_stata=True, layout=layout,
        )
        np.testing.assert_allclose(
            e.income, 1.0 + derivative[0] / e.we,
            atol=3e-14, rtol=3e-14,
        )

    def test_subset_reordered_prices_omitted_expenditure_and_distinct_cf(self):
        rng = np.random.default_rng(8642)
        n, R, N = 4, 2, 9
        d0, _ = _data(rng, n=n, R=R, N=N, residual=True)
        selection_cf = rng.normal(0.0, 0.35, N)
        spec = Spec(n, R, True, True, True)
        theta = rng.normal(0.0, 0.025, spec.n_free)
        names = tuple(f"p{j + 1}" for j in range(n))
        layout = FirstStageLayout(
            ordered_names=("p3", "p1", "z1", "selection_cf"),
            demand_price_names=names,
            price_positions={"p3": 0, "p1": 1},
            expenditure_position=None,
            covariate_positions={"z1": 2},
            selection_cf_position=3,
            constant_position=4,
        )
        tau = rng.normal(0.0, 0.08, n * layout.width)

        def predict(lnp, lnexp):
            X = np.column_stack(
                [lnp[:, 2], lnp[:, 0], d0.demo[:, 0], selection_cf]
            )
            k = np.column_stack([
                X @ tau[layout.equation_slice(i)][:-1]
                + tau[layout.equation_slice(i)][-1]
                for i in range(n)
            ])
            d = DemandData(
                lnp=lnp, lnexp=lnexp, shares=d0.shares, demo=d0.demo,
                cdf=norm.cdf(k), pdf=norm.pdf(k), a0=d0.a0,
                control_function=d0.control_function,
            )
            return fitted_shares(theta, d, spec), d, k

        _, d, k = predict(d0.lnp, d0.lnexp)
        d_exp, d_price = fitted_share_derivatives(
            theta, d, spec, tau=tau, layout=layout, selection_index=k
        )
        h = 1e-5
        plus, _, _ = predict(d0.lnp, d0.lnexp + h)
        minus, _, _ = predict(d0.lnp, d0.lnexp - h)
        np.testing.assert_allclose(
            d_exp, (plus - minus) / (2 * h), atol=2e-10, rtol=2e-9
        )
        for j in range(n):
            lp, lm = d0.lnp.copy(), d0.lnp.copy()
            lp[:, j] += h
            lm[:, j] -= h
            plus, _, _ = predict(lp, d0.lnexp)
            minus, _, _ = predict(lm, d0.lnexp)
            np.testing.assert_allclose(
                d_price[:, :, j], (plus - minus) / (2 * h),
                atol=2e-10, rtol=2e-9,
            )


class SyntheticAndGateTests(unittest.TestCase):
    def test_estimator_recovers_known_kappa_as_sample_grows_and_starts_agree(self):
        def fit(frame, start):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", RuntimeWarning)
                return quaidsce(
                    frame,
                    shares=["w1", "w2", "w3"],
                    prices=["p1", "p2", "p3"],
                    expenditure="m",
                    demographics=["z"],
                    control_function="v",
                    selection_control_function="v",
                    anot=1.6,
                    first_stage_predict="xb",
                    reps=0,
                    method="ifgnls",
                    start=start,
                    strict_stata=False,
                    max_outer=80,
                    max_iter=160,
                    verbose=False,
                )

        truth = np.array([0.12, -0.10, 0.08])
        errors = np.empty((4, 2, 3))
        large_seed_17 = None
        for s, seed in enumerate((5, 17, 29, 41)):
            nested, _ = _synthetic_cf_dgp(2400, seed)
            for j, N in enumerate((600, 2400)):
                result = fit(nested.iloc[:N], "zero")
                self.assertTrue(result.converged)
                errors[s, j] = result.coefs.cfcoef - truth
                if seed == 17 and N == 2400:
                    large_seed_17 = result

        rmse_small = np.sqrt(np.mean(errors[:, 0] ** 2))
        rmse_large = np.sqrt(np.mean(errors[:, 1] ** 2))
        large_estimates = errors[:, 1] + truth
        self.assertLess(rmse_large, 0.025)
        self.assertLessEqual(rmse_large, 0.75 * rmse_small)
        np.testing.assert_array_equal(
            np.sign(large_estimates),
            np.broadcast_to(np.sign(truth), large_estimates.shape),
        )

        self.assertIsNotNone(large_seed_17)
        nested_17, _ = _synthetic_cf_dgp(2400, 17)
        large_linear = fit(nested_17, "linear")
        self.assertTrue(large_linear.converged)
        np.testing.assert_allclose(
            large_seed_17.coefs.cfcoef, large_linear.coefs.cfcoef,
            atol=1e-4, rtol=1e-4,
        )
        np.testing.assert_allclose(
            large_seed_17.theta, large_linear.theta,
            atol=1e-4, rtol=1e-4,
        )

    def test_augmented_latent_share_direction_is_exact(self):
        rng = np.random.default_rng(71)
        spec = Spec(4, 2, True, True, True)
        d, _ = _data(rng, residual=True)
        theta = rng.normal(0.0, 0.02, spec.n_free)
        xs = free_slices(spec)["cfcoef"]
        theta[xs] = np.array([0.4, -0.2, 0.1, -0.3])
        difference = augmented_latent_shares(theta, d, spec) - latent_shares(
            theta, d, spec
        )
        np.testing.assert_allclose(
            difference,
            d.control_function[:, None] * theta[xs][None, :],
            atol=2e-16,
            rtol=2e-15,
        )

    def test_invalid_extension_combinations_fail_before_fit(self):
        import pandas as pd

        rng = np.random.default_rng(8)
        N = 40
        frame = pd.DataFrame({
            "w1": np.where(np.arange(N) % 3, 0.3, 0.0),
            "w2": np.where(np.arange(N) % 4, 0.4, 0.0),
            "w3": np.where(np.arange(N) % 5, 0.3, 0.0),
            "p1": np.exp(rng.normal(size=N)),
            "p2": np.exp(rng.normal(size=N)),
            "p3": np.exp(rng.normal(size=N)),
            "m": np.exp(rng.normal(2, 0.2, N)),
            "z": rng.normal(size=N),
            "cf": rng.normal(size=N),
            "constant_cf": 0.0,
        })
        frame["cf_collinear"] = np.log(frame["m"])
        common = dict(data=frame, shares=["w1", "w2", "w3"],
                      prices=["p1", "p2", "p3"], expenditure="m",
                      demographics=["z"], anot=1.6, verbose=False)
        with self.assertRaisesRegex(ValueError, "first_stage_predict='xb'"):
            quaidsce(**common, control_function="cf", first_stage_predict="pr")
        with self.assertRaisesRegex(ValueError, "reduced form residual"):
            quaidsce(**common, selection_control_function="cf",
                      first_stage_predict="xb", reps=2)
        with self.assertRaisesRegex(ValueError, "censor=True"):
            quaidsce(**common, control_function="cf", censor=False,
                      first_stage_predict="xb")
        with self.assertRaisesRegex(ValueError, "subset of demand prices"):
            quaidsce(**common, selection_prices=["not_a_price"],
                      first_stage_predict="xb")
        with self.assertRaisesRegex(ValueError, "nonzero variation"):
            quaidsce(**common, control_function="constant_cf",
                      first_stage_predict="xb")
        with self.assertRaisesRegex(ValueError, "collinear.*dropped"):
            quaidsce(**common, selection_control_function="cf_collinear",
                      first_stage_predict="xb")
        frame.loc[0, "cf"] = np.nan
        with self.assertRaisesRegex(ValueError, "only finite values"):
            quaidsce(**common, control_function="cf",
                      first_stage_predict="xb")

    def test_active_estimator_runs_end_to_end(self):
        import pandas as pd
        from pathlib import Path

        root = Path(__file__).resolve().parents[1]
        frame = pd.read_stata(root / "bench" / "small4.dta")
        lnp = np.log(frame[["p1", "p2", "p4", "p9"]].to_numpy(float))
        lnm = np.log(frame["total"].to_numpy(float))
        x1 = frame["x1"].to_numpy(float)
        reduced_form = np.column_stack([
            np.ones(len(frame)), lnp, frame[["x1", "x2"]].to_numpy(float),
            x1**2,
        ])
        frame["cfunc"] = lnm - reduced_form @ np.linalg.lstsq(
            reduced_form, lnm, rcond=None
        )[0]
        frame["selection_cfunc"] = (
            0.7 * frame["cfunc"].to_numpy(float)
            + np.random.default_rng(678).normal(0.0, 0.03, len(frame))
        )
        res = quaidsce(
            frame,
            shares=["sw1", "sw2", "sw4", "sw9"],
            prices=["p1", "p2", "p4", "p9"],
            expenditure="total",
            demographics=["x1", "x2"],
            control_function="cfunc",
            selection_control_function="selection_cfunc",
            selection_expenditure=True,
            anot=10.0,
            first_stage_predict="xb",
            method="ifgnls",
            start="linear",
            verbose=False,
        )
        self.assertTrue(res.converged)
        self.assertTrue(np.isfinite(res.b).all())
        self.assertEqual(len(res.names), len(res.b))
        self.assertEqual(res.V.shape, (len(res.b), len(res.b)))
        self.assertTrue(np.isfinite(res.V).all())
        n_elasticities = res.spec.neqn + 2 * res.spec.neqn ** 2
        self.assertTrue(np.isnan(res.analytic_se[-n_elasticities:]).all())
        self.assertTrue(np.isfinite(res.analytic_se[:-n_elasticities]).all())
        self.assertEqual(res.selection_layout.selection_cf_position, 7)
        self.assertTrue(any("conditional" in note.lower() for note in res.notes))
        self.assertIn("conditional", res.summary().lower())

        selection_only = quaidsce(
            frame,
            shares=["sw1", "sw2", "sw4", "sw9"],
            prices=["p1", "p2", "p4", "p9"],
            expenditure="total",
            demographics=["x1", "x2"],
            selection_control_function="selection_cfunc",
            selection_expenditure=True,
            anot=10.0,
            first_stage_predict="xb",
            method="ifgnls",
            start="linear",
            verbose=False,
        )
        self.assertTrue(selection_only.converged)
        self.assertFalse(any(
            name.startswith("cfcoef:")
            for name in selection_only.spec.full_names()
        ))
        self.assertIsNotNone(
            selection_only.selection_layout.selection_cf_position
        )

    def test_near_zero_augmented_elasticity_denominator_fails_loudly(self):
        spec = Spec(3, 1, False, True, True)
        coefs = unpack(np.zeros(spec.n_free), spec)
        layout = legacy_layout(["p1", "p2", "p3"], ["z"],
                               include_expenditure=True)
        means = Means(
            w=np.zeros(3), lnp=np.zeros(3), lnexp=1.0, demo=np.zeros(1),
            cdf=np.full(3, 0.5), pdf=np.zeros(3), du=np.zeros(3),
            control_function=0.0,
        )
        with self.assertRaisesRegex(ValueError, "too close to zero"):
            elasticities(
                coefs, spec, means, 1.6,
                tau=np.zeros(3 * layout.width), np_prob=layout.width,
                layout=layout, strict_stata=False,
            )


if __name__ == "__main__":
    unittest.main()
