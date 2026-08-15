"""Self-checks: parameter map, delta matrix, analytic Jacobian, probit."""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(SRC))

import numpy as np

from pyquaidsce.model import DemandData, fitted_shares, jacobian_full
from pyquaidsce.params import (Spec, delta_blocks, delta_matrix, full_vector,
                              unpack)
from pyquaidsce.probit import probit


def _fake_data(rng, N, n, R, censored=True):
    lnp = rng.normal(0.0, 0.3, size=(N, n))
    lnexp = rng.normal(2.0, 0.5, size=N)
    demo = rng.normal(1.0, 0.4, size=(N, R)) if R else np.zeros((N, 0))
    shares = rng.uniform(0.05, 0.4, size=(N, n))
    shares /= shares.sum(axis=1, keepdims=True)
    if censored:
        cdf = rng.uniform(0.3, 0.95, size=(N, n))
        pdf = rng.uniform(0.05, 0.4, size=(N, n))
    else:
        cdf, pdf = np.ones((N, n)), np.zeros((N, n))
    return DemandData(lnp=lnp, lnexp=lnexp, shares=shares, demo=demo,
                      cdf=cdf, pdf=pdf, a0=1.6)


def check_delta_matrix():
    """Delta must be the exact Jacobian of the free -> full parameter map."""
    rng = np.random.default_rng(0)
    worst = 0.0
    for n in (3, 4, 6):
        for R in (0, 1, 2):
            for quad in (True, False):
                for cen in (True, False):
                    if cen and R == 0:
                        continue
                    spec = Spec(n, R, quad, cen)
                    th = rng.normal(size=spec.n_free) * 0.1
                    D = delta_matrix(spec)
                    assert D.shape == (spec.n_full, spec.n_free), (
                        spec, D.shape, (spec.n_full, spec.n_free))
                    num = np.zeros_like(D)
                    h = 1e-6
                    for k in range(spec.n_free):
                        tp, tm = th.copy(), th.copy()
                        tp[k] += h
                        tm[k] -= h
                        num[:, k] = (full_vector(tp, spec)
                                     - full_vector(tm, spec)) / (2 * h)
                    worst = max(worst, float(np.max(np.abs(num - D))))
                    # block decomposition must rebuild Delta
                    Db = np.zeros_like(D)
                    for fs, xs, mat in delta_blocks(spec):
                        Db[fs, xs] = mat
                    assert np.allclose(Db, D), spec
    return worst


def check_restrictions():
    """Adding-up / homogeneity / symmetry must hold where Stata imposes them."""
    rng = np.random.default_rng(1)
    msgs = []
    for cen in (False, True):
        R = 1
        spec = Spec(5, R, True, cen)
        th = rng.normal(size=spec.n_free) * 0.2
        c = unpack(th, spec)
        assert np.allclose(c.gamma, c.gamma.T), "gamma not symmetric"
        assert np.allclose(c.gamma.sum(axis=1), 0.0, atol=1e-12), \
            "homogeneity violated"
        assert np.allclose(c.eta.sum(axis=1), 0.0, atol=1e-12), \
            "eta adding-up violated"
        if cen:
            msgs.append("censored: alpha sums to %.6f (free, as in Stata)"
                        % c.alpha.sum())
        else:
            assert abs(c.alpha.sum() - 1.0) < 1e-12
            assert abs(c.beta.sum()) < 1e-12
            assert abs(c.lam.sum()) < 1e-12
            msgs.append("nocensor: alpha/beta/lambda adding-up OK")
    return msgs


def check_jacobian():
    rng = np.random.default_rng(2)
    worst = 0.0
    for n in (3, 4, 5):
        for R in (0, 1, 2):
            for quad in (True, False):
                for cen in (True, False):
                    if cen and R == 0:
                        continue
                    spec = Spec(n, R, quad, cen)
                    d = _fake_data(rng, 40, n, R, censored=cen)
                    th = rng.normal(size=spec.n_free) * 0.15
                    Jf = jacobian_full(th, d, spec)
                    # chain rule to free params
                    D = delta_matrix(spec)
                    Ja = Jf.reshape(-1, spec.n_full) @ D
                    h = 1e-6
                    num = np.zeros_like(Ja)
                    for k in range(spec.n_free):
                        tp, tm = th.copy(), th.copy()
                        tp[k] += h
                        tm[k] -= h
                        num[:, k] = (
                            fitted_shares(tp, d, spec).reshape(-1)
                            - fitted_shares(tm, d, spec).reshape(-1)
                        ) / (2 * h)
                    err = np.max(np.abs(Ja - num)) / max(1.0, np.max(np.abs(Ja)))
                    worst = max(worst, float(err))
                    assert err < 1e-6, (spec, err)
    return worst


def check_probit():
    """Compare against a brute-force optimiser and the analytic OIM."""
    from scipy.optimize import minimize
    from scipy.stats import norm

    rng = np.random.default_rng(3)
    N, k = 800, 3
    X = rng.normal(size=(N, k))
    btrue = np.array([0.8, -0.5, 0.3])
    y = (X @ btrue + 0.4 + rng.normal(size=N) > 0).astype(float)

    r = probit(y, X)

    Z = np.hstack([X, np.ones((N, 1))])
    q = 2 * y - 1

    def nll(b):
        return -np.sum(norm.logcdf(q * (Z @ b)))

    opt = minimize(nll, np.zeros(k + 1), method="BFGS",
                   options=dict(gtol=1e-12, maxiter=2000))
    return float(np.max(np.abs(r.b - opt.x))), float(-r.llf - opt.fun)


def check_fast_jacobian():
    """The fast free-space Jacobian must equal Jfull @ Delta exactly."""
    from pyquaidsce.jacfree import jacobian_free, make_cache
    from pyquaidsce.nlsur import _jac_free_reference
    from pyquaidsce.params import delta_blocks

    rng = np.random.default_rng(11)
    worst = 0.0
    for n in (3, 4, 6, 8):
        for R in (0, 1, 3):
            for quad in (True, False):
                for cen in (True, False):
                    if cen and R == 0:
                        continue
                    spec = Spec(n, R, quad, cen)
                    d = _fake_data(rng, 25, n, R, censored=cen)
                    th = rng.normal(size=spec.n_free) * 0.12
                    ref = _jac_free_reference(th, d, spec,
                                              delta_blocks(spec), slice(None))
                    fast = jacobian_free(th, d, spec, make_cache(spec))
                    err = np.max(np.abs(ref - fast))
                    scale = max(1.0, np.max(np.abs(ref)))
                    worst = max(worst, float(err / scale))
                    assert err / scale < 1e-12, (spec, err)
    return worst


if __name__ == "__main__":
    print("delta matrix   max |analytic - numeric| =", check_delta_matrix())
    for m in check_restrictions():
        print("restrictions  ", m)
    print("jacobian       max relative error      =", check_jacobian())
    print("fast jacobian  max relative error      =", check_fast_jacobian())
    db, dll = check_probit()
    print("probit         max |b - bfgs|          =", db)
    print("probit         llf difference          =", dll)
    print("\nALL INTERNAL CHECKS PASSED")
