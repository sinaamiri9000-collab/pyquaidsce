"""Uncensored control functions: restrictions, derivatives and inference."""
from dataclasses import replace
from pathlib import Path
import sys
import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pyquaidsce import quaidsce
from pyquaidsce import estimator as estimator_module
from pyquaidsce.elasticities import elasticities, sample_means
from pyquaidsce.jacfree import jacobian_free, make_cache
from pyquaidsce.model import DemandData, augmented_latent_shares, fitted_shares, jacobian_full
from pyquaidsce.params import Spec, delta_blocks, delta_matrix, free_slices, full_vector, unpack
from pyquaidsce.start import linear_start


def synthetic_data(nobs=350):
    rng = np.random.default_rng(829)
    spec = Spec(3, 1, True, False, True)
    sl = free_slices(spec)
    theta = np.zeros(spec.n_free)
    for name, values in {
        "alpha": [.32, .30], "beta": [.015, -.01],
        "gamma": [.012, -.003, .01], "lambda": [.0008, -.0005],
        "cfcoef": [.08, -.05], "eta": [.006, -.003], "rho": [.02],
    }.items():
        theta[sl[name]] = values
    lnp = rng.normal(0, .18, (nobs, 3))
    z = rng.normal(0, .2, (nobs, 1))
    instrument = rng.normal(0, .4, nobs)
    v = rng.normal(0, .15, nobs)
    lnexp = 2 + lnp @ np.array([.12, -.08, .06]) + .1*z[:, 0] + .5*instrument + v
    d = DemandData(lnp, lnexp, np.zeros((nobs, 3)), z,
                   np.ones((nobs, 3)), np.zeros((nobs, 3)), 1.6, v)
    shares = augmented_latent_shares(theta, d, spec)
    error = rng.normal(0, .001, (nobs, 2))
    shares += np.column_stack([error, -error.sum(axis=1)])
    assert np.min(shares) > 0
    d.shares = shares
    frame = pd.DataFrame({**{f"w{i+1}": shares[:, i] for i in range(3)},
                          **{f"p{i+1}": np.exp(lnp[:, i]) for i in range(3)},
                          "m": np.exp(lnexp), "z": z[:, 0], "income": instrument, "v": v})
    return frame, d, spec, theta


def options(frame):
    return dict(data=frame, shares=["w1", "w2", "w3"], prices=["p1", "p2", "p3"],
                expenditure="m", demographics=["z"], anot=1.6, censor=False,
                method="ifgnls", start="linear", max_outer=150, verbose=False)


class UncensoredControlFunctionTests(unittest.TestCase):
    def test_maps_and_both_jacobians_match_finite_differences(self):
        rng = np.random.default_rng(12)
        for n in [3, 4]:
            for r in [0, 2]:
                for quadratic in [False, True]:
                    with self.subTest(n=n, demographics=r, quadratic=quadratic):
                        spec = Spec(n, r, quadratic, False, True)
                        theta = rng.normal(0, .012, spec.n_free)
                        theta[free_slices(spec)["alpha"]] = 1/n
                        d = DemandData(rng.normal(0, .12, (9, n)), rng.normal(2, .1, 9),
                                       np.full((9,n), 1/n), rng.normal(0, .1, (9,r)),
                                       np.ones((9,n)), np.zeros((9,n)), 1.6, rng.normal(.2, .1, 9))
                        coefs = unpack(theta, spec)
                        self.assertAlmostEqual(coefs.cfcoef.sum(), 0, places=14)
                        np.testing.assert_allclose(augmented_latent_shares(theta, d, spec).sum(axis=1), 1, atol=3e-15)
                        delta = delta_matrix(spec)
                        block_delta = np.zeros_like(delta)
                        for fs, xs, block in delta_blocks(spec): block_delta[fs, xs] = block
                        np.testing.assert_array_equal(delta, block_delta)
                        h = 1e-6
                        numeric = np.empty((9,n-1,spec.n_free))
                        for j in range(spec.n_free):
                            step = np.zeros(spec.n_free); step[j] = h
                            numeric[:,:,j] = (fitted_shares(theta+step,d,spec)-fitted_shares(theta-step,d,spec))/(2*h)
                            np.testing.assert_allclose((full_vector(theta+step,spec)-full_vector(theta-step,spec))/(2*h),delta[:,j],atol=5e-10,rtol=1e-8)
                        fast = jacobian_free(theta,d,spec,make_cache(spec))
                        full = jacobian_full(theta,d,spec) @ delta
                        np.testing.assert_allclose(fast, numeric, atol=3e-9, rtol=3e-7)
                        np.testing.assert_allclose(full, numeric, atol=3e-9, rtol=3e-7)
                        start = linear_start(d,spec)
                        self.assertEqual(start.shape,(spec.n_free,))
                        self.assertTrue(np.isfinite(start).all())

    def test_elasticities_match_quantity_derivatives_with_residual_held_fixed(self):
        _, d, spec, theta = synthetic_data(20)
        point = DemandData(d.lnp.mean(axis=0,keepdims=True), np.array([d.lnexp.mean()]),
                           np.zeros((1,3)), d.demo.mean(axis=0,keepdims=True),
                           np.ones((1,3)), np.zeros((1,3)), 1.6, np.array([.23]))
        point.shares = augmented_latent_shares(theta,point,spec)
        means = sample_means(point,None,spec)
        e = elasticities(unpack(theta,spec),spec,means,1.6)
        def log_quantity(x):
            return np.log(augmented_latent_shares(theta,x,spec)[0])+x.lnexp[0]-x.lnp[0]
        h=1e-6
        numeric_income=(log_quantity(replace(point,lnexp=point.lnexp+h))-log_quantity(replace(point,lnexp=point.lnexp-h)))/(2*h)
        numeric_price=np.zeros((3,3))
        for j in range(3):
            step=np.zeros_like(point.lnp);step[0,j]=h
            numeric_price[:,j]=(log_quantity(replace(point,lnp=point.lnp+step))-log_quantity(replace(point,lnp=point.lnp-step)))/(2*h)
        np.testing.assert_allclose(e.income,numeric_income,atol=2e-9)
        np.testing.assert_allclose(e.uncompensated,numeric_price,atol=2e-9)
        self.assertAlmostEqual(float(means.w @ e.income),1,places=12)
        np.testing.assert_allclose(e.compensated.sum(axis=1),0,atol=2e-12)

    def test_zero_control_coefficients_nest_the_existing_uncensored_model(self):
        _, d, spec, theta = synthetic_data(25)
        old = Spec(3,1,True,False,False)
        sl, old_sl = free_slices(spec),free_slices(old)
        theta[sl['cfcoef']]=0
        old_theta=np.zeros(old.n_free)
        for key in old_sl: old_theta[old_sl[key]]=theta[sl[key]]
        np.testing.assert_array_equal(fitted_shares(theta,d,spec),fitted_shares(old_theta,d,old))
        means=sample_means(d,None,spec)
        np.testing.assert_array_equal(elasticities(unpack(theta,spec),spec,means,1.6).as_stata_vector(),elasticities(unpack(old_theta,old),old,means,1.6).as_stata_vector())

    def test_internal_and_external_controls_match_and_recover_coefficients(self):
        frame,_,_,theta=synthetic_data()
        args=options(frame)
        internal=quaidsce(**args,ivexp=['income'])
        frame['residual']=internal.reduced_form.residuals
        external=quaidsce(**args,control_function='residual')
        known=quaidsce(**args,control_function='v')
        self.assertTrue(internal.converged and external.converged and known.converged)
        np.testing.assert_allclose(internal.theta,external.theta,atol=1e-10)
        np.testing.assert_allclose(known.coefs.cfcoef,[.08,-.05,-.03],atol=.004)
        self.assertAlmostEqual(internal.coefs.cfcoef.sum(),0,places=14)
        self.assertIsNone(internal.selection_layout)
        self.assertIsNone(internal.selection_control_function_name)
        self.assertEqual(internal.probits,[])
        self.assertEqual(len(internal.names),len(internal.b))
        self.assertTrue(np.isnan(internal.analytic_se[-21:]).all())
        self.assertNotIn('ELAS_',internal.summary(elasticities=False))

    def test_bootstrap_rebuilds_reduced_form_and_parallel_results_match(self):
        frame,_,_,_=synthetic_data(180)
        args=options(frame)
        original=estimator_module.fit_expenditure_reduced_form
        with patch.object(estimator_module,'fit_expenditure_reduced_form',wraps=original) as counted:
            serial=quaidsce(**args,ivexp=['income'],reps=3,seed=573,n_jobs=1,bootstrap_start='warm')
        self.assertEqual(counted.call_count,4)
        self.assertEqual(serial.boot.reps_ok,3)
        self.assertTrue(np.isfinite(serial.se).all())
        parallel=quaidsce(**args,ivexp=['income'],reps=3,seed=573,n_jobs=2,bootstrap_start='warm')
        np.testing.assert_allclose(parallel.boot.b_star,serial.boot.b_star,atol=1e-12)
        self.assertEqual(parallel.boot.b_star.shape,(3,len(parallel.names)))

    def test_estimates_do_not_depend_on_which_equation_is_omitted(self):
        frame, _, _, _ = synthetic_data(240)
        args = options(frame)
        original = quaidsce(**args, ivexp=['income'])
        order = [2, 0, 1]
        args.update(shares=[f'w{i+1}' for i in order], prices=[f'p{i+1}' for i in order])
        reordered = quaidsce(**args, ivexp=['income'])
        self.assertTrue(original.converged and reordered.converged)
        reverse = np.argsort(order)
        np.testing.assert_allclose(original.coefs.cfcoef, reordered.coefs.cfcoef[reverse], atol=2e-7)
        np.testing.assert_allclose(original.elas.income, reordered.elas.income[reverse], atol=2e-7)
        np.testing.assert_allclose(original.elas.uncompensated, reordered.elas.uncompensated[np.ix_(reverse, reverse)], atol=2e-7)


if __name__=='__main__': unittest.main()
