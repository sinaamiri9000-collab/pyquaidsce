"""Independent cluster sandwich and whole-cluster resampling checks."""

import sys
import unittest
from pathlib import Path
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'src'), str(ROOT)]

import numpy as np
import pandas as pd

from pyquaidsce import quaidsce
from pyquaidsce._clusters import cluster_codes, cluster_rows, resample_indices
from pyquaidsce.bootstrap import _one_rep, _WORK
from pyquaidsce.inference import _EstimatingSystem, _pack_means, compute_analytical_inference
from pyquaidsce.params import delta_matrix
from tools.validate_small4 import BENCH, DEMOS, PRICES, SHARES


class ClusterInferenceTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.frame = pd.read_stata(Path(BENCH)/'small4.dta')
        # Unequal, interleaved clusters spanning every computation chunk.
        cls.labels = np.arange(len(cls.frame)) % 37
        cls.labels[::11] = 36
        cls.frame['group'] = cls.labels.astype(str)
        cls.kw = dict(shares=SHARES, prices=PRICES, expenditure='total',
                      demographics=DEMOS, anot=10.0, verbose=False)
        cls.iid = quaidsce(cls.frame, **cls.kw, analytic=True)
        def capture(result, data, design, **kw):
            cls.data, cls.design = data, design
            return compute_analytical_inference(result, data, design, **kw)
        with patch('pyquaidsce.inference.compute_analytical_inference', side_effect=capture):
            cls.fit = quaidsce(cls.frame, **cls.kw, analytic=True, cluster='group')

    def test_independent_full_influence_sandwich(self):
        result = self.fit
        system = _EstimatingSystem(result, self.data, self.design, None, 10000, None)
        _, scores, features = next(system.chunks(system.x, include_means=True))
        centered = scores-scores.mean(axis=0)
        inverse = np.linalg.inv(system.bread())
        parameter_if = -centered @ inverse.T
        direct, means_jac = system.elasticity_jacobians()
        means_if = features-_pack_means(result.means) + parameter_if @ system.mean_derivative().T
        elasticity_if = parameter_if @ direct.T + means_if @ means_jac.T
        k, nt = result.spec.n_free, system.tau.size
        joint_if = np.zeros((len(scores), k+nt))
        joint_if[:, :k] = parameter_if[:, system.theta_slice]
        joint_if[:, k+system.active] = parameter_if[:, system.tau_slice]
        full_if = np.column_stack([joint_if[:, :k] @ delta_matrix(result.spec).T,
                                   joint_if[:, k:], elasticity_if])
        groups = np.unique(self.labels)
        factor = len(groups)/(len(groups)-1)

        def reference(values, denominator):
            # Direct masked sums, independent of the chunked accumulation path.
            sums = np.vstack([values[self.labels == g].sum(axis=0) for g in groups])
            return factor * sums.T @ sums / denominator

        np.testing.assert_allclose(result.V, reference(full_if, len(scores)**2),
                                   rtol=2e-7, atol=1e-9)
        np.testing.assert_allclose(result.analytical.joint_covariance,
                                   reference(joint_if, len(scores)**2), rtol=2e-7, atol=1e-9)
        np.testing.assert_allclose(result.analytical.elasticity_covariance,
                                   reference(elasticity_if, len(scores)**2), rtol=2e-7, atol=1e-9)
        np.testing.assert_allclose(result.analytical.meat,
                                   reference(centered, len(scores)), rtol=2e-9, atol=1e-10)
        np.testing.assert_allclose(result.analytical.bread, self.iid.analytical.bread, rtol=0, atol=0)

    def test_cluster_option_preserves_point_estimates(self):
        for name in ['theta', 'b', 'tau', 'sigma']:
            np.testing.assert_array_equal(getattr(self.iid,name),getattr(self.fit,name))
        for name in ['llf', 'converged', 'n_outer', 'n_gn']:
            self.assertEqual(getattr(self.iid,name),getattr(self.fit,name))
        self.assertEqual(self.fit.n_clusters, 37)
        self.assertEqual(self.fit.analytical.covariance_type, 'cluster')
        self.assertIn('Number of clusters', self.fit.summary())

    def test_singleton_clusters_cr0_equal_iid(self):
        inf = compute_analytical_inference(self.iid, self.data, self.design,
            clusters=np.arange(len(self.frame)), cluster_correction=False)
        np.testing.assert_allclose(inf.covariance, self.iid.V, rtol=2e-12, atol=1e-12)
        corrected = compute_analytical_inference(self.iid, self.data, self.design,
            clusters=np.arange(len(self.frame)))
        n = len(self.frame)
        np.testing.assert_allclose(corrected.covariance, inf.covariance*n/(n-1),
                                   rtol=2e-12, atol=1e-12)

    def test_relabel_chunk_and_row_order_invariance(self):
        inf = compute_analytical_inference(self.fit, self.data, self.design,
            clusters=np.array([f'label-{999-g}' for g in self.labels]), chunk=23)
        np.testing.assert_allclose(inf.covariance, self.fit.V, rtol=2e-7, atol=1e-8)
        order = np.random.default_rng(348).permutation(len(self.frame))
        other = compute_analytical_inference(self.fit, self.data.subset(order),
            self.design[order], clusters=self.labels[order], chunk=67)
        np.testing.assert_allclose(other.covariance, self.fit.V, rtol=2e-7, atol=1e-8)

    def test_estimation_sample_alignment_and_invalid_labels(self):
        frame = self.frame.copy()
        frame.loc[0,'total'] = np.nan
        frame.loc[0,'group'] = None
        masked = quaidsce(frame, **self.kw, analytic=True, cluster='group')
        explicit = quaidsce(frame.iloc[1:], **self.kw, analytic=True, cluster='group')
        np.testing.assert_array_equal(masked.b, explicit.b)
        np.testing.assert_array_equal(masked.V, explicit.V)
        for labels in [np.zeros(len(frame)), [None]*len(frame), [np.inf]*len(frame)]:
            with self.subTest(labels=labels[0]):
                bad = self.frame.assign(group=labels)
                with self.assertRaises(ValueError):
                    quaidsce(bad, **self.kw, analytic=True, cluster='group')
        with self.assertRaisesRegex(ValueError, 'requires'):
            quaidsce(self.frame, **self.kw, cluster='group')
        with self.assertRaisesRegex(ValueError, 'existing'):
            quaidsce(self.frame, **self.kw, analytic=True, cluster='absent')
        with self.assertRaisesRegex(ValueError, 'True or False'):
            quaidsce(self.frame, **self.kw, analytic=True, cluster_correction=1)

    def test_fgnls_and_uncensored_cluster_paths(self):
        for options in [dict(method='fgnls'), dict(censor=False), dict(method='nls')]:
            with self.subTest(options=options):
                iid = quaidsce(self.frame, **self.kw, **options, analytic=True)
                cluster = quaidsce(self.frame, **self.kw, **options, analytic=True, cluster='group')
                np.testing.assert_array_equal(iid.b,cluster.b)
                self.assertEqual(cluster.analytical.n_clusters, 37)
                self.assertTrue(np.isfinite(cluster.analytical.elasticity_se).all())
                self.assertGreater(np.max(np.abs(iid.V-cluster.V)), 1e-8)

    def test_whole_cluster_sampling_and_worker_integration(self):
        labels = np.array(['a','b','a','c','b','b'])
        codes,g = cluster_codes(labels,len(labels))
        groups = cluster_rows(codes,g)
        seed=132
        expected_draws = np.random.default_rng(seed).integers(0, g, size=g)
        expected = np.concatenate([np.flatnonzero(codes==x) for x in expected_draws])
        np.testing.assert_array_equal(resample_indices(np.random.default_rng(seed),len(labels),groups),expected)
        for draw in range(20):
            idx=resample_indices(np.random.default_rng(draw),len(labels),groups)
            # Repeated clusters must contain every household the same number of times.
            counts=np.bincount(idx,minlength=len(labels))
            for rows in groups:
                self.assertTrue(np.all(counts[rows]==counts[rows[0]]))
        frame=pd.DataFrame({'row':np.arange(len(labels))})
        captured=[]
        class Result:
            converged=True
            b=np.array([1.0,2.0])
        def fit(df,**kw):
            captured.append(df['row'].to_numpy())
            return Result()
        with patch.dict(_WORK, dict(df=frame,kw={},rep_timeout=None,cluster_rows=groups),clear=True):
            with patch('pyquaidsce.estimator.quaidsce',side_effect=fit):
                _,b,error=_one_rep((1,seed))
        self.assertIsNone(error)
        np.testing.assert_array_equal(captured[0],expected)
        np.testing.assert_array_equal(b,[1.0,2.0])
        # Legacy IID draws keep their original RNG sequence.
        np.testing.assert_array_equal(resample_indices(np.random.default_rng(seed),6),
                                   np.random.default_rng(seed).integers(0,6,size=6))

    def test_native_bootstrap_keeps_analytical_cluster_reference(self):
        fit=quaidsce(self.frame, **self.kw, analytic=True, cluster='group', reps=4, seed=734)
        self.assertEqual(fit.boot.n_clusters,37)
        self.assertEqual(fit.boot.cluster_name,'group')
        self.assertGreaterEqual(fit.boot.reps_ok,2)
        np.testing.assert_array_equal(fit.V_analytic, self.fit.V)
        np.testing.assert_allclose(fit.V,np.cov(fit.boot.b_star,rowvar=False,ddof=1))
        # The four production worker draws also match explicit whole-cluster fits.
        seeds=np.random.SeedSequence(734).generate_state(4,dtype=np.uint32)
        codes,g=cluster_codes(self.frame['group'],len(self.frame))
        groups=cluster_rows(codes,g)
        expected=[]
        failed=[]
        for i,seed in enumerate(seeds,1):
            idx=resample_indices(np.random.default_rng(int(seed)),len(self.frame),groups)
            result=quaidsce(self.frame.iloc[idx].reset_index(drop=True),**self.kw)
            if result.converged:
                expected.append(result.b)
            else:
                failed.append(f'rep {i}: RuntimeError: second-stage estimator did not converge')
        np.testing.assert_array_equal(fit.boot.b_star,expected)
        self.assertEqual(fit.boot.failures,failed)


if __name__ == '__main__':
    unittest.main()
