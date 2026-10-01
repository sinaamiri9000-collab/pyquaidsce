"""Latent adding-up identities, constrained derivatives and bootstrap refits."""

import unittest
from unittest.mock import patch

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce, translog
from pyquaidsce.data import DemandData
from pyquaidsce.inference import complex_step_jacobian
from pyquaidsce.jacfree import jacobian_free, make_cache
from pyquaidsce.model import (augmented_latent_shares, fitted_shares,
                             jacobian_full)
from pyquaidsce.params import (Spec, delta_blocks, delta_matrix, free_slices,
                              full_slices, full_vector, unpack)
from pyquaidsce.start import linear_start
from pyquaidsce.translog_model import (
    TRANSLOG_BACKEND, TranslogCoefs, TranslogSpec,
    delta_matrix as translog_delta, fitted_shares as translog_fitted,
    full_vector as translog_full, jacobian_free as translog_jacobian,
    latent_shares as translog_latent, pack, share_derivatives,
)


def demand_data():
    rng = np.random.default_rng(891)
    n, size = 4, 19
    return DemandData(
        rng.normal(0, .2, (size, n)), rng.normal(2, .2, size),
        np.full((size, n), 1/n), rng.uniform(.1, .4, (size, 2)),
        rng.uniform(.3, .8, (size, n)), rng.uniform(.1, .3, (size, n)),
        a0=1.6, control_function=rng.normal(0, .2, size),
    )


class LatentAddingMathTests(unittest.TestCase):
    def test_quaids_constraints_and_covariance_maps(self):
        d = demand_data()
        rng = np.random.default_rng(31)
        for quadratic in (False, True):
            for cf in (False, True):
                with self.subTest(quadratic=quadratic, cf=cf):
                    s = Spec(4, 2, quadratic, True, cf, latent_adding=True)
                    old = Spec(4, 2, quadratic, True, cf)
                    self.assertEqual(old.n_free-s.n_free, 2+int(quadratic)+int(cf))
                    self.assertEqual(s.n_eq_estimated, 4)
                    t = rng.normal(0, .015, s.n_free)
                    c = unpack(t, s)
                    self.assertAlmostEqual(c.alpha.sum(), 1)
                    for block in (c.beta, c.lam, c.cfcoef):
                        self.assertAlmostEqual(block.sum(), 0)
                    np.testing.assert_allclose(c.gamma.sum(axis=0), 0, atol=1e-15)
                    np.testing.assert_allclose(c.eta.sum(axis=1), 0, atol=1e-15)
                    np.testing.assert_allclose(augmented_latent_shares(t, d, s).sum(axis=1), 1, atol=1e-14)
                    self.assertGreater(np.max(abs(fitted_shares(t, d, s).sum(axis=1)-1)), .1)
                    D = delta_matrix(s)
                    numerical = np.column_stack([
                        (full_vector(t+row*1e-6, s)-full_vector(t-row*1e-6, s))/2e-6
                        for row in np.eye(s.n_free)
                    ])
                    np.testing.assert_allclose(D, numerical, atol=1e-10)
                    rebuilt = np.zeros_like(D)
                    for fs, xs, block in delta_blocks(s):
                        rebuilt[fs, xs] = block
                    np.testing.assert_array_equal(rebuilt, D)
                    V = D @ D.T
                    fs = full_slices(s)
                    for name in ('alpha', 'beta', 'lambda', 'cfcoef'):
                        if name in fs:
                            selector = np.zeros(s.n_full)
                            selector[fs[name]] = 1
                            np.testing.assert_allclose(selector @ V, 0, atol=1e-14)
                    np.testing.assert_array_equal(D[fs['delta'], free_slices(s)['delta']], np.eye(4))

    def test_quaids_all_equation_jacobian_with_control_function(self):
        d = demand_data()
        rng = np.random.default_rng(14)
        for quadratic in (False, True):
            for cf in (False, True):
                with self.subTest(quadratic=quadratic, cf=cf):
                    s = Spec(4, 2, quadratic, True, cf, latent_adding=True)
                    t = rng.normal(0, .012, s.n_free)
                    analytic = jacobian_free(t, d, s, make_cache(s))
                    full = jacobian_full(t, d, s).reshape(-1, s.n_full) @ delta_matrix(s)
                    np.testing.assert_allclose(analytic.reshape(full.shape), full, atol=1e-14)
                    h = 1e-6
                    numeric = np.stack([
                        (fitted_shares(t+row*h, d, s)-fitted_shares(t-row*h, d, s))/(2*h)
                        for row in np.eye(s.n_free)
                    ], axis=-1)
                    np.testing.assert_allclose(analytic, numeric, atol=3e-10)
                    np.testing.assert_allclose(jacobian_free(t, d, s, make_cache(s), slice(3, 7)), analytic[3:7])

    def test_constrained_linear_and_translog_mean_starts(self):
        d = demand_data()
        for quadratic in (False, True):
            s = Spec(4, 2, quadratic, True, True, latent_adding=True)
            t = linear_start(d, s)
            self.assertTrue(np.isfinite(t).all())
            np.testing.assert_allclose(augmented_latent_shares(t, d, s).sum(axis=1), 1, atol=1e-12)
        s = TranslogSpec(4, 2, True, latent_adding=True)
        t = TRANSLOG_BACKEND.initial(d, s, 'mean')
        self.assertEqual(t.size, s.n_free)
        np.testing.assert_allclose(translog_latent(t, d, s).sum(axis=1), 1)

    def test_translog_constraints_all_equation_derivatives_and_inference(self):
        d = demand_data()
        s = TranslogSpec(4, 2, True, latent_adding=True)
        t = pack(TranslogCoefs(np.array([.2, .3, .25, .25]), np.eye(4)*.003,
                              np.full((4, 2), .006), np.arange(1, 5)*.01), s)
        self.assertEqual(s.n_free, TranslogSpec(4, 2, True).n_free-1)
        self.assertEqual(s.n_eq_estimated, 4)
        np.testing.assert_allclose(translog_latent(t, d, s).sum(axis=1), 1, atol=1e-14)
        dm, dp = share_derivatives(t, d, s)
        np.testing.assert_allclose(dm.sum(axis=1), 0, atol=1e-14)
        np.testing.assert_allclose(dp.sum(axis=1), 0, atol=1e-14)
        analytic = translog_jacobian(t, d, s)
        numeric = complex_step_jacobian(lambda v: translog_fitted(v, d, s), t).reshape(19, 4, -1)
        np.testing.assert_allclose(analytic, numeric, atol=1e-13)
        self.assertTrue(np.all(analytic[:, -1, :3] < 0))
        D = translog_delta(s)
        np.testing.assert_allclose(D, complex_step_jacobian(lambda v: translog_full(v, s), t), atol=1e-14)
        np.testing.assert_allclose(D[:4].sum(axis=0), 0)

    def test_option_validation_and_uncensored_mandatory_adding(self):
        for cls in (Spec, TranslogSpec):
            with self.subTest(cls=cls.__name__):
                with self.assertRaisesRegex(ValueError, 'requires censor=True'):
                    cls(4, 2, censor=False, latent_adding=True)
                with self.assertRaisesRegex(ValueError, 'boolean'):
                    cls(4, 2, censor=True, latent_adding='false')
                self.assertTrue(cls(4, 2, censor=False).latent_adding_up)
                self.assertFalse(cls(4, 2, censor=True).latent_adding_up)
        opts = dict(shares=['w1', 'w2', 'w3'], prices=['p1', 'p2', 'p3'],
                    expenditure='m', demographics=['z'], censor=False, latent_adding=True)
        for fit, extra in ((quaidsce, {'anot': 1.6}), (translog, {})):
            with self.assertRaisesRegex(ValueError, 'requires censor=True'):
                fit(pd.DataFrame(), **opts, **extra)


class LatentAddingEstimationTests(unittest.TestCase):
    def test_quaids_ivexp_and_bootstrap_preserve_latent_constraints(self):
        from tests.test_control_function import _synthetic_cf_dgp
        from pyquaidsce import estimator as module
        frame, _ = _synthetic_cf_dgp(650, 452)
        specs = []
        actual = module.nlsur
        def capture(d, spec, *args, **kwargs):
            specs.append(spec)
            return actual(d, spec, *args, **kwargs)
        with patch.object(module, 'nlsur', side_effect=capture):
            r = quaidsce(frame, ['w1', 'w2', 'w3'], prices=['p1', 'p2', 'p3'],
                expenditure='m', demographics=['z'], ivexp=['instrument'], anot=1.6,
                latent_adding=True, method='nls', max_iter=160, reps=2,
                bootstrap_start='warm', n_jobs=1, seed=812, verbose=False)
        self.assertTrue(r.converged)
        self.assertEqual(r.boot.reps_ok, 2)
        self.assertEqual(len(specs), 3)
        self.assertTrue(all(s.latent_adding and s.control_function for s in specs))
        self.assertAlmostEqual(r.coefs.cfcoef.sum(), 0)
        for names in ('alpha', 'beta', 'lambda', 'cfcoef'):
            idx = [i for i, name in enumerate(r.names) if name.startswith(names+':')]
            expected = 1 if names == 'alpha' else 0
            np.testing.assert_allclose(r.boot.b_star[:, idx].sum(axis=1), expected, atol=1e-12)
        self.assertTrue(np.isfinite(r.boot.V).all())

    def test_translog_nested_and_spawn_bootstrap_preserve_constraint(self):
        from tests.test_translog_sy import censored_sample, OPTIONS
        frame, _, _ = censored_sample(size=450, seed=913)
        r = translog(frame, **OPTIONS, latent_adding=True, start='nested',
                     reps=2, seed=511, bootstrap_start='warm', n_jobs=2,
                     mp_context='spawn')
        self.assertTrue(r.converged)
        self.assertEqual(r.boot.reps_ok, 2)
        self.assertTrue(r.initialization['converged'])
        np.testing.assert_allclose(r.predict(kind='latent_shares').sum(axis=1), 1, atol=1e-13)
        np.testing.assert_allclose(r.boot.b_star[:, :3].sum(axis=1), 1, atol=1e-13)
        np.testing.assert_allclose(r.predict(frame), r.predict(), atol=1e-12)
        self.assertEqual(r.sigma.shape, (3, 3))
        self.assertTrue(np.isfinite(r.se).all())


if __name__ == '__main__':
    unittest.main()
