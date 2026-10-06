"""Independent checks for joint Probit/demand and elasticity inference."""

import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.stats import norm

from pyquaidsce import DemandData, quaidsce
from pyquaidsce.elasticities import elasticities, sample_means
from pyquaidsce.inference import (_EstimatingSystem, _bread_inverse,
                                  _central_jacobian, compute_analytical_inference)
from pyquaidsce.model import fitted_shares, jacobian_full
from pyquaidsce.params import delta_matrix, unpack
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES


def _context(frame, result, covariates=DEMOS):
    lnp = np.log(frame[PRICES].to_numpy())
    lnexp = np.log(frame['total'].to_numpy())
    demo = frame[DEMOS].to_numpy()
    design = np.column_stack([lnp, lnexp, frame[list(covariates)].to_numpy()])
    z = np.column_stack([design, np.ones(len(frame))])
    xb = z @ result.tau.reshape(result.spec.neqn, result.np_prob).T
    data = DemandData(lnp, lnexp, frame[SHARES].to_numpy(), demo,
                      norm.cdf(xb), norm.pdf(xb), result.anot)
    return data, design


class AnalyticalInferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame = pd.read_stata(Path(BENCH) / 'small4.dta')
        cls.kw = dict(shares=SHARES, prices=PRICES, expenditure='total',
                      demographics=DEMOS, anot=10.0, verbose=False)
        cls.base = quaidsce(cls.frame, **cls.kw)
        cls.fit = quaidsce(cls.frame, **cls.kw, analytic=True)
        cls.data, cls.design = _context(cls.frame, cls.fit)

    def test_opt_in_preserves_all_point_estimation_results(self):
        old, new = self.base, self.fit
        for name in ['theta', 'b', 'b_est', 'sigma', 'tau']:
            np.testing.assert_array_equal(getattr(old, name), getattr(new, name))
        for name in ['llf', 'n_outer', 'n_gn', 'converged']:
            self.assertEqual(getattr(old, name), getattr(new, name))
        self.assertIsNone(old.analytical)
        self.assertTrue(np.isnan(old.analytic_se[-36:]).all())
        self.assertTrue(np.isfinite(new.analytic_se).all())
        np.testing.assert_allclose(new.se**2, np.diag(new.V), rtol=1e-14)

    def test_joint_and_elasticity_covariances_have_cross_terms(self):
        result = self.fit
        k = result.spec.n_free
        joint = result.analytical.joint_covariance
        self.assertGreater(np.max(np.abs(joint[:k, k:])), 1e-8)
        width = result.np_prob
        self.assertGreater(np.max(np.abs(joint[k:k+width, k+width:k+2*width])), 1e-8)
        cov = result.analytical.elasticity_covariance
        np.testing.assert_allclose(cov, result.V[-36:, -36:], atol=1e-13)
        self.assertGreater(np.min(np.diag(cov)), 0)
        self.assertGreaterEqual(np.linalg.eigvalsh(cov).min(), -1e-10)
        np.testing.assert_allclose(result.analytical.uncompensated_se.ravel(),
                                   result.analytic_se[-32:-16])

    def test_scores_match_full_jacobian_and_likelihood_differences(self):
        system = _EstimatingSystem(self.fit, self.data, self.design, None, 89, None)
        sl, scores, _ = next(system.chunks(system.x))
        tau, theta, sigma = system.decode(system.x)
        sub = self.data.subset(sl)
        full = jacobian_full(theta, sub, self.fit.spec)
        jac = full @ delta_matrix(self.fit.spec)
        u = sub.shares - fitted_shares(theta, sub, self.fit.spec)
        expected = np.einsum('tmi,tm->ti', jac, u @ np.linalg.inv(sigma))
        np.testing.assert_allclose(scores[:, system.theta_slice], expected,
                                   rtol=1e-10, atol=1e-10)
        z = np.column_stack([self.design[sl], np.ones(sub.nobs)])
        q = 2*(sub.shares[:, 0] > 0)-1
        coef = tau[:self.fit.np_prob]
        for j in [0, self.fit.np_prob-1]:
            step = np.zeros_like(coef)
            step[j] = 1e-5
            numeric = (norm.logcdf(q*(z @ (coef+step)))
                       - norm.logcdf(q*(z @ (coef-step)))) / 2e-5
            np.testing.assert_allclose(scores[:, j], numeric, rtol=1e-6, atol=1e-8)

    def test_probit_bread_matches_observed_information(self):
        inf = self.fit.analytical
        width = self.fit.np_prob
        for i, pr in enumerate(self.fit.probits):
            sl = slice(i*width, (i+1)*width)
            np.testing.assert_allclose(inf.bread[sl, sl],
                                       -np.linalg.inv(pr.V)/self.fit.nobs,
                                       rtol=2e-6, atol=2e-8)

    def test_generated_mean_chain_matches_rebuilding_probit_inputs(self):
        fit = self.fit
        system = _EstimatingSystem(fit, self.data, self.design, None, 2000, None)
        direct, mean_jac = system.elasticity_jacobians()
        total = direct + mean_jac @ system.mean_derivative()
        z = np.column_stack([self.design, np.ones(self.data.nobs)])

        def rebuilt(active_tau):
            tau = fit.tau.copy()
            tau[system.active] = active_tau
            xb = z @ tau.reshape(fit.spec.neqn, fit.np_prob).T
            d = self.data.subset(slice(None))
            d.cdf, d.pdf = norm.cdf(xb), norm.pdf(xb)
            means = sample_means(d, xb, fit.spec)
            return elasticities(unpack(fit.theta, fit.spec), fit.spec, means,
                                fit.anot, tau=tau, np_prob=fit.np_prob,
                                layout=fit.selection_layout).as_stata_vector()

        independent = _central_jacobian(rebuilt, fit.tau[system.active],
                                        relative_step=2e-6)
        np.testing.assert_allclose(total[:, system.tau_slice], independent,
                                   rtol=2e-5, atol=2e-6)
        self.assertGreater(np.max(np.abs(total-direct)), 1e-4)

    def test_step_size_and_chunk_size_stability(self):
        other = compute_analytical_inference(self.fit, self.data, self.design,
                                             chunk=137,
                                             _relative_step=np.cbrt(np.finfo(float).eps)/2)
        np.testing.assert_allclose(other.elasticity_se,
                                   self.fit.analytical.elasticity_se,
                                   rtol=2e-4, atol=2e-7)

    def test_analytical_bread_matches_independent_numerical_derivative(self):
        for method in ['ifgnls', 'nls', 'fgnls']:
            with self.subTest(method=method):
                result = self.fit if method == 'ifgnls' else quaidsce(
                    self.frame, **self.kw, method=method)
                data, design = _context(self.frame, result)
                initial = (quaidsce(self.frame, **self.kw, method='nls').theta
                           if method == 'fgnls' else None)
                system = _EstimatingSystem(result, data, design, initial, 2000, None)
                h = np.cbrt(np.finfo(float).eps)/4
                small = _central_jacobian(system.mean_score, system.x,
                                          system.scales, h)
                large = _central_jacobian(system.mean_score, system.x,
                                          system.scales, 2*h)
                reference = (4*small-large)/3
                analytic = system.bread()
                # Weight derivatives suffer cancellation in real differences
                # when the demand scores are close to zero. Check that block
                # independently with complex perturbations of Sigma instead.
                if system.sigma_slice is not None:
                    _, theta, sigma = system.decode(system.x)
                    jac = jacobian_full(theta, data, result.spec) @ delta_matrix(result.spec)
                    u = data.shares[:, :result.spec.n_eq_estimated] - fitted_shares(theta, data, result.spec)
                    for column, (a, b) in enumerate(zip(*system.tril)):
                        h = 1e-20*system.scales[system.sigma_slice.start+column]
                        perturbed = sigma.astype(complex)
                        perturbed[a, b] += 1j*h
                        if a != b:
                            perturbed[b, a] += 1j*h
                        gradient = np.einsum('tmi,tm->i', jac, u @ np.linalg.inv(perturbed))
                        reference[system.theta_slice, system.sigma_slice.start+column] = gradient.imag/h/data.nobs
                scaled_error = np.max(np.abs(analytic-reference)/np.maximum(1, np.abs(analytic)))
                self.assertLess(scaled_error, 5e-6)

    def test_bootstrap_keeps_the_joint_analytical_reference(self):
        fit = quaidsce(self.frame, **self.kw, analytic=True,
                      reps=4, seed=17031, n_jobs=1)
        self.assertIsNotNone(fit.boot)
        np.testing.assert_array_equal(fit.V_analytic, fit.analytical.covariance)
        np.testing.assert_array_equal(fit.V, fit.boot.V)
        np.testing.assert_allclose(fit.analytic_se, fit.analytical.se)
        np.testing.assert_array_equal(fit.se, fit.boot.se)

    def test_sandwich_matches_independent_two_step_block_formula(self):
        # General two-step block expansion, rather than another matrix solve
        # through the same production covariance code.
        rng = np.random.default_rng(19)
        a = np.array([[-2.0, 0.3], [0.3, -1.0]])
        c = rng.normal(size=(3, 2))
        d = -np.diag([2.0, 1.0, 3.0])
        bread = np.block([[a, np.zeros((2, 3))], [c, d]])
        scores = rng.normal(size=(400, 5)) @ rng.normal(size=(5, 5))
        meat = scores.T @ scores / len(scores)
        inverse, _ = _bread_inverse(bread, np.ones(5))
        full = inverse @ meat @ inverse.T / len(scores)
        ai, di = np.linalg.inv(a), np.linalg.inv(d)
        b11, b12, b22 = meat[:2, :2], meat[:2, 2:], meat[2:, 2:]
        correction = (b22 + c @ ai @ b11 @ ai.T @ c.T
                      - c @ ai @ b12 - b12.T @ ai.T @ c.T)
        expected = di @ correction @ di.T / len(scores)
        np.testing.assert_allclose(full[2:, 2:], expected, rtol=1e-12, atol=1e-14)

    def test_dropped_probit_columns_keep_exact_zero_variance(self):
        frame = self.frame.copy()
        frame['duplicate_x1'] = frame['x1']
        fit = quaidsce(frame, **self.kw, analytic=True,
                      selection_covariates=['x1', 'duplicate_x1'])
        for i, pr in enumerate(fit.probits):
            self.assertTrue(pr.dropped)
            for dropped in pr.dropped:
                index = fit.spec.n_full + i*fit.np_prob + dropped
                np.testing.assert_array_equal(fit.V[index], np.zeros(fit.V.shape[1]))

    def test_nls_fgnls_linear_and_uncensored_support(self):
        for options in [dict(method='nls'), dict(method='fgnls'),
                        dict(quadratic=False), dict(censor=False),
                        dict(method='nls', algorithm='lm'),
                        dict(selection_covariates=['x1'])]:
            with self.subTest(options=options):
                old = quaidsce(self.frame, **self.kw, **options)
                new = quaidsce(self.frame, **self.kw, **options, analytic=True)
                np.testing.assert_array_equal(old.theta, new.theta)
                self.assertEqual(old.n_gn, new.n_gn)
                self.assertTrue(np.isfinite(new.analytical.elasticity_se).all())
                self.assertTrue(np.isfinite(new.analytic_se).all())

    def test_unsupported_stages_and_nonconvergence_are_explicit(self):
        for options in [dict(ivexp=['x2']), dict(control_function='x2'),
                        dict(selection_control_function='x2')]:
            with self.subTest(options=options):
                with self.assertRaisesRegex(NotImplementedError, 'ivexp/control'):
                    quaidsce(self.frame, **self.kw, analytic=True, **options)
        with self.assertRaisesRegex(ValueError, 'converged'):
            quaidsce(self.frame, **self.kw, analytic=True, max_outer=2)
        with self.assertRaisesRegex(ValueError, 'analytic must'):
            quaidsce(self.frame, **self.kw, analytic='yes')


if __name__ == '__main__':
    unittest.main()
