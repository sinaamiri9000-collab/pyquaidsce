"""Joint estimating-equation covariance and elasticity delta-method inference.

The construction follows Hardin (2002), Stata Journal 2(3), 253--266.
First-stage scores, demand scores and estimated GLS-weight moments are stacked.
Sample-mean influence functions include their dependence on Probit parameters.
Estimating-equation and elasticity derivatives are analytical. Independent
central numerical differences are retained as a verification path. The
covariance is an asymptotic sandwich, without resampling.
"""

from __future__ import annotations

from dataclasses import dataclass
import time
from typing import Optional

import numpy as np
from scipy.stats import norm

from ._clusters import cluster_codes
from ._elasticity_derivatives import elasticity_jacobian
from ._threads import with_blas_threads
from ._timing import check_deadline
from .elasticities import Means, elasticities
from .jacfree import jacobian_free, make_cache
from .model import DemandData, _inner, fitted_shares
from .params import delta_matrix, full_slices, unpack, vech_index
from .probit import _lambda_ratio


@dataclass
class AnalyticalInference:
    """Joint observation-level or cluster sandwich inference at sample means.

    ``covariance`` follows the existing result-vector order. Elasticity arrays
    use the natural [good, price] order, also used by ``result.elas``.
    ``joint_covariance`` orders free demand parameters before the full Probit
    vector, including zero rows/columns for explicitly dropped regressors.
    """

    covariance: np.ndarray
    joint_covariance: np.ndarray
    elasticity_covariance: np.ndarray
    bread: np.ndarray
    meat: np.ndarray
    estimating_parameters: np.ndarray
    score_mean: np.ndarray
    bread_condition: float
    max_standardized_score: float
    elapsed_seconds: float
    neqn: int
    n_clusters: Optional[int] = None
    covariance_type: str = "iid"
    correction_factor: float = 1.0

    @property
    def se(self) -> np.ndarray:
        return np.sqrt(np.diag(self.covariance))

    @property
    def elasticity_se(self) -> np.ndarray:
        return np.sqrt(np.diag(self.elasticity_covariance))

    @property
    def income_se(self) -> np.ndarray:
        return self.elasticity_se[:self.neqn]

    @property
    def uncompensated_se(self) -> np.ndarray:
        n = self.neqn
        return self.elasticity_se[n:n + n*n].reshape(n, n)

    @property
    def compensated_se(self) -> np.ndarray:
        n = self.neqn
        return self.elasticity_se[n + n*n:].reshape(n, n)


def _central_jacobian(function, x, scales=None, relative_step=None):
    """Central derivative with scale-aware steps and domain checks."""
    x = np.asarray(x, dtype=float)
    scales = np.maximum(1.0, np.abs(x)) if scales is None else np.asarray(scales)
    step = np.cbrt(np.finfo(float).eps) if relative_step is None else relative_step
    base = np.asarray(function(x), dtype=float)
    if not np.isfinite(base).all():
        raise ValueError("nonfinite value in analytical inference")
    out = np.empty((base.size, x.size))
    for j in range(x.size):
        h = float(step * scales[j])
        for attempt in range(12):
            plus, minus = x.copy(), x.copy()
            plus[j] += h
            minus[j] -= h
            try:
                with np.errstate(over="raise", invalid="raise", divide="raise"):
                    yp, ym = np.asarray(function(plus)), np.asarray(function(minus))
                valid = np.isfinite(yp).all() and np.isfinite(ym).all()
            except (FloatingPointError, np.linalg.LinAlgError):
                valid = False
            if valid:
                out[:, j] = (yp - ym).ravel() / (2.0 * h)
                break
            h *= 0.5
        else:
            raise ValueError(f"cannot differentiate inference parameter {j}")
    return out


def _bread_inverse(bread, scales):
    """Equilibrated inverse; never regularize unidentified parameters."""
    scaled = bread * scales[None, :]
    rows = np.linalg.norm(scaled, axis=1)
    if np.any(rows == 0) or not np.isfinite(rows).all():
        raise ValueError("analytical inference has an unidentified parameter")
    scaled = scaled / rows[:, None]
    s = np.linalg.svd(scaled, compute_uv=False)
    if s[-1] <= np.finfo(float).eps * len(s) * s[0]:
        raise ValueError("analytical inference requires a full-rank estimating system")
    inverse = np.linalg.solve(scaled, np.eye(len(s)))
    inverse = scales[:, None] * inverse / rows[None, :]
    return inverse, float(s[0] / s[-1])


def _pack_means(means):
    return np.concatenate([means.w, means.lnp, [means.lnexp], means.demo,
                           means.cdf, means.pdf, means.du])


def _unpack_means(values, n, r):
    offset = 2*n + 1 + r
    return Means(values[:n], values[n:2*n], float(values[2*n]),
                 values[2*n + 1:offset], values[offset:offset+n],
                 values[offset+n:offset+2*n], values[offset+2*n:offset+3*n])


class _DemandDerivatives:
    """Exact latent Jacobian and residual-weighted Hessian in free coordinates.

    For D=ln(m)-A-ln(mbar), d2D=dmbar dmbar'/mbar**2.
    For q=exp(-l), d2q=q dl dl'. Apply the product rule to
    w*=alpha+Gamma ln(p)+B D+lambda q D**2. Hessians are contracted
    with the residual weights before allocation, avoiding an N*m*K*K array.
    """

    def __init__(self, spec):
        self.spec = spec
        n, r, k = spec.neqn, spec.ndemo, spec.n_free
        mapping, slices = delta_matrix(spec), full_slices(spec)
        self.alpha = mapping[slices['alpha']]
        self.beta = mapping[slices['beta']]
        self.gamma = np.zeros((n, n, k))
        for row, (i, j) in enumerate(vech_index(n)):
            self.gamma[i, j] = self.gamma[j, i] = mapping[slices['gamma']][row]
        self.lam = mapping[slices['lambda']] if spec.quadratic else np.zeros((n, k))
        self.delta = mapping[slices['delta']] if spec.censor else np.zeros((n, k))
        self.eta = mapping[slices['eta']].reshape(r, n, k) if r else np.zeros((0, n, k))
        self.rho = mapping[slices['rho']] if r else np.zeros((0, k))

    def evaluate(self, theta, data, weights):
        spec, m = self.spec, weights.shape[1]
        coefs = unpack(theta, spec)
        inn = _inner(coefs, data, spec)
        lnp, demo = data.lnp, data.demo
        da = lnp @ self.alpha + 0.5*np.einsum('ti,ijk,tj->tk',
                                             lnp, self.gamma, lnp)
        dm = demo @ self.rho
        dd = -da - dm/inn.mbar[:, None]
        db = self.beta[None, :, :] + np.einsum('tr,rik->tik', demo, self.eta)
        jac = (self.alpha[None, :, :]
               + np.einsum('tj,ijk->tik', lnp, self.gamma)
               + inn.B[:, :, None]*dd[:, None, :]
               + inn.D[:, None, None]*db)

        # Sum_i weights_i H(B_i D), including the residual-Hessian term.
        weighted_db = np.einsum('ti,tik->tk', weights, db[:, :m])
        weighted_b = np.einsum('ti,ti->t', weights, inn.B[:, :m])
        hessian = dd.T @ weighted_db + weighted_db.T @ dd
        hessian += dm.T @ ((weighted_b/inn.mbar**2)[:, None]*dm)
        if spec.quadratic:
            dl = lnp @ self.beta + np.einsum('tr,rik,ti->tk', demo, self.eta, lnp)
            df = inn.q[:, None]*(2*inn.D[:, None]*dd - inn.D[:, None]**2*dl)
            f = inn.q*inn.D**2
            jac += coefs.lam[None, :, None]*df[:, None, :] + f[:, None, None]*self.lam
            weighted_lam = weights @ coefs.lam[:m]
            weighted_dlam = weights @ self.lam[:m]
            hessian += df.T @ weighted_dlam + weighted_dlam.T @ df
            qlam = inn.q*weighted_lam
            hessian += dd.T @ ((2*qlam)[:, None]*dd)
            hessian += dm.T @ ((2*qlam*inn.D/inn.mbar**2)[:, None]*dm)
            cross = dd.T @ ((-2*qlam*inn.D)[:, None]*dl)
            hessian += cross + cross.T
            hessian += dl.T @ ((qlam*inn.D**2)[:, None]*dl)
        return jac[:, :m], hessian


class _EstimatingSystem:
    """Scores for the existing estimator, with no optimization calls."""

    def __init__(self, result, data, design, theta_nls, chunk, deadline):
        self.result, self.data, self.spec = result, data, result.spec
        self.chunk, self.deadline = chunk, deadline
        self.cache = make_cache(self.spec)
        self.derivatives = _DemandDerivatives(self.spec)
        n, k = self.spec.neqn, self.spec.n_free
        self.tau = np.asarray(result.tau if self.spec.censor else [], dtype=float)
        self.design = None
        self.active = []
        self.good_columns = []
        if self.spec.censor:
            self.design = np.column_stack([design, np.ones(data.nobs)])
            for i, pr in enumerate(result.probits):
                keep = [j for j in range(result.np_prob) if j not in pr.dropped]
                self.good_columns.append(keep)
                self.active.extend(i*result.np_prob + j for j in keep)
        self.active = np.asarray(self.active, dtype=int)
        self.tau_slice = slice(0, len(self.active))
        offset = len(self.active)
        parts = [self.tau[self.active]]
        self.nls_slice = None
        if result.method == "fgnls":
            if theta_nls is None:
                raise ValueError("FGNLS inference needs the initial NLS estimates")
            self.nls_slice = slice(offset, offset+k)
            offset += k
            parts.append(np.asarray(theta_nls))
        self.theta_slice = slice(offset, offset+k)
        offset += k
        parts.append(result.theta)
        self.sigma_slice = None
        self.tril = np.tril_indices(self.spec.n_eq_estimated)
        if result.method != "nls":
            source = result.theta if theta_nls is None else theta_nls
            u = data.shares[:, :self.spec.n_eq_estimated] - fitted_shares(
                source, data, self.spec)
            sigma = u.T @ u / data.nobs
            parts.append(sigma[self.tril])
            self.sigma_slice = slice(offset, offset + len(self.tril[0]))
        self.x = np.concatenate(parts)
        self.scales = np.maximum(1.0, np.abs(self.x))
        if self.sigma_slice is not None:
            sigma_scale = np.sqrt(np.outer(np.diag(sigma), np.diag(sigma)))
            self.scales[self.sigma_slice] = sigma_scale[self.tril]

    def decode(self, vector):
        tau = self.tau.copy()
        tau[self.active] = vector[self.tau_slice]
        theta = vector[self.theta_slice]
        m = self.spec.n_eq_estimated
        sigma = np.eye(m)
        if self.sigma_slice is not None:
            sigma[self.tril] = vector[self.sigma_slice]
            sigma[(self.tril[1], self.tril[0])] = vector[self.sigma_slice]
            np.linalg.cholesky(sigma)
        return tau, theta, sigma

    def chunks(self, vector, include_means=False):
        tau, theta, sigma = self.decode(vector)
        weight = np.linalg.inv(sigma)
        n, m = self.spec.neqn, self.spec.n_eq_estimated
        for start in range(0, self.data.nobs, self.chunk):
            check_deadline(self.deadline)
            sl = slice(start, min(start+self.chunk, self.data.nobs))
            sub = self.data.subset(sl)
            scores = []
            if self.spec.censor:
                z = self.design[sl]
                xb = z @ tau.reshape(n, self.result.np_prob).T
                sub.cdf, sub.pdf = norm.cdf(xb), norm.pdf(xb)
                q = 2.0 * (sub.shares > 0) - 1.0
                lam = _lambda_ratio(q * xb)
                for i, keep in enumerate(self.good_columns):
                    scores.append(z[:, keep] * (q[:, i]*lam[:, i])[:, None])
            else:
                xb = np.ones((sub.nobs, n))
            if self.nls_slice is not None:
                initial = vector[self.nls_slice]
                u_initial = sub.shares[:, :m] - fitted_shares(initial, sub, self.spec)
                j_initial = jacobian_free(initial, sub, self.spec, self.cache)
                scores.append(np.einsum("tmi,tm->ti", j_initial, u_initial))
            u = sub.shares[:, :m] - fitted_shares(theta, sub, self.spec)
            j = jacobian_free(theta, sub, self.spec, self.cache)
            scores.append(np.einsum("tmi,tm->ti", j, u @ weight))
            if self.sigma_slice is not None:
                source_u = u if self.nls_slice is None else u_initial
                a, b = self.tril
                scores.append(source_u[:, a] * source_u[:, b] - sigma[self.tril])
            features = None
            if include_means:
                features = np.column_stack([sub.shares, sub.lnp, sub.lnexp,
                                            sub.demo, sub.cdf, sub.pdf, xb])
            yield sl, np.column_stack(scores), features

    def mean_score(self, vector):
        out = np.zeros_like(vector)
        for _, scores, _ in self.chunks(vector):
            out += scores.sum(axis=0)
        return out / self.data.nobs

    def bread(self):
        """Mean derivative of every score, including residual Hessians."""
        out = np.zeros((self.x.size, self.x.size))
        tau, theta, sigma = self.decode(self.x)
        weight = np.linalg.inv(sigma)
        n, m = self.spec.neqn, self.spec.n_eq_estimated
        for start in range(0, self.data.nobs, self.chunk):
            check_deadline(self.deadline)
            sl = slice(start, min(start+self.chunk, self.data.nobs))
            sub = self.data.subset(sl)
            z = self.design[sl] if self.spec.censor else None
            xb = (z @ tau.reshape(n, self.result.np_prob).T
                  if self.spec.censor else np.ones((sub.nobs, n)))
            if self.spec.censor:
                sub.cdf, sub.pdf = norm.cdf(xb), norm.pdf(xb)
                position = 0
                for i, keep in enumerate(self.good_columns):
                    q = 2.0*(sub.shares[:, i] > 0)-1.0
                    index = q*xb[:, i]
                    lam = _lambda_ratio(index)
                    information = lam*(lam+index)
                    # Match the derivative of the backend's tail expansion.
                    tail = index <= -30
                    information[tail] = 1-1/index[tail]**2+6/index[tail]**4
                    block = slice(position, position+len(keep))
                    zk = z[:, keep]
                    out[block, block] -= zk.T @ (information[:, None]*zk)
                    position += len(keep)

            def demand_blocks(parameters, w):
                u = sub.shares[:, :m] - fitted_shares(parameters, sub, self.spec)
                j = jacobian_free(parameters, sub, self.spec, self.cache)
                wu = u @ w
                jw = np.matmul(w, j)
                residual_weights = wu*sub.cdf[:, :m] if self.spec.censor else wu
                jl, hessian = self.derivatives.evaluate(parameters, sub, residual_weights)
                hh = hessian - np.einsum('tmi,tmj->ij', j, jw)
                cross = np.zeros((self.spec.n_free, len(self.active)))
                fitted_tau = []
                if self.spec.censor:
                    coefs = unpack(parameters, self.spec)
                    latent = _inner(coefs, sub, self.spec).wstar
                    position = 0
                    for i, keep in enumerate(self.good_columns):
                        correction = latent[:, i]-coefs.delta[i]*xb[:, i]
                        derivative = sub.pdf[:, i, None]*(
                            wu[:, i, None]*(jl[:, i]-xb[:, i, None]*self.derivatives.delta[i])
                            - jw[:, i]*correction[:, None])
                        block = slice(position, position+len(keep))
                        cross[:, block] = derivative.T @ z[:, keep]
                        fitted_tau.append(sub.pdf[:, i, None]*correction[:, None]*z[:, keep])
                        position += len(keep)
                return hh, cross, j, u, jw, wu, fitted_tau

            blocks = demand_blocks(theta, weight)
            hh, cross, j, u, jw, wu, fitted_tau = blocks
            out[self.theta_slice, self.theta_slice] += hh
            out[self.theta_slice, self.tau_slice] += cross
            if self.nls_slice is not None:
                initial_blocks = demand_blocks(self.x[self.nls_slice], np.eye(m))
                out[self.nls_slice, self.nls_slice] += initial_blocks[0]
                out[self.nls_slice, self.tau_slice] += initial_blocks[1]
            if self.sigma_slice is not None:
                for column, (a, b) in enumerate(zip(*self.tril)):
                    value = -jw[:, a].T @ wu[:, b]
                    if a != b:
                        value -= jw[:, b].T @ wu[:, a]
                    out[self.theta_slice, self.sigma_slice.start+column] += value
                source = blocks if self.nls_slice is None else initial_blocks
                js, us, fitted_tau = source[2], source[3], source[6]
                source_slice = self.theta_slice if self.nls_slice is None else self.nls_slice
                for row, (a, b) in enumerate(zip(*self.tril)):
                    dest = self.sigma_slice.start+row
                    out[dest, source_slice] -= (js[:, a]*us[:, b, None]
                                                + js[:, b]*us[:, a, None]).sum(axis=0)
                    position = 0
                    for i, keep in enumerate(self.good_columns):
                        block = slice(position, position+len(keep))
                        if i == a:
                            out[dest, block] -= us[:, b] @ fitted_tau[i]
                        if i == b:
                            out[dest, block] -= us[:, a] @ fitted_tau[i]
                        position += len(keep)
                out[self.sigma_slice, self.sigma_slice] -= sub.nobs*np.eye(len(self.tril[0]))
        return out / self.data.nobs

    def mean_derivative(self):
        """Dependence of generated sample means on the active Probit vector."""
        n, r = self.spec.neqn, self.spec.ndemo
        out = np.zeros((5*n + r + 1, self.x.size))
        if not self.spec.censor:
            return out
        offset = 2*n + 1 + r
        xb = self.design @ self.tau.reshape(n, self.result.np_prob).T
        phi = norm.pdf(xb)
        for column, flat in enumerate(self.active):
            i, j = divmod(int(flat), self.result.np_prob)
            z = self.design[:, j]
            out[offset+i, column] = np.mean(phi[:, i] * z)
            out[offset+n+i, column] = np.mean(-xb[:, i] * phi[:, i] * z)
            out[offset+2*n+i, column] = np.mean(z)
        return out

    def elasticity_jacobians(self, relative_step=None):
        n, r = self.spec.neqn, self.spec.ndemo
        nt = len(self.active)
        theta = self.result.theta
        means = _pack_means(self.result.means)
        x = np.concatenate([self.tau[self.active], theta, means])

        def function(values):
            check_deadline(self.deadline)
            tau = self.tau.copy()
            tau[self.active] = values[:nt]
            coefs = unpack(values[nt:nt+theta.size], self.spec)
            evaluation = _unpack_means(values[nt+theta.size:], n, r)
            return elasticities(coefs, self.spec, evaluation, self.result.anot,
                                tau=tau if self.spec.censor else None,
                                np_prob=self.result.np_prob,
                                layout=self.result.selection_layout).as_stata_vector()

        if relative_step is None:
            check_deadline(self.deadline)
            jac = elasticity_jacobian(
                theta, self.spec, self.result.means, self.result.anot,
                self.tau, self.result.np_prob, self.result.selection_layout,
                self.active,
            )
        else:
            # Private independent numerical reference; never used by fitting.
            jac = _central_jacobian(function, x, relative_step=relative_step)
        core = np.zeros((jac.shape[0], self.x.size))
        core[:, self.tau_slice] = jac[:, :nt]
        core[:, self.theta_slice] = jac[:, nt:nt+theta.size]
        return core, jac[:, nt+theta.size:]


@with_blas_threads
def compute_analytical_inference(
    result,
    data: DemandData,
    selection_design: Optional[np.ndarray] = None,
    *,
    theta_nls: Optional[np.ndarray] = None,
    chunk: int = 2000,
    blas_threads: Optional[int] = 1,
    deadline: Optional[float] = None,
    clusters: Optional[np.ndarray] = None,
    cluster_correction: bool = True,
    _relative_step: Optional[float] = None,
) -> AnalyticalInference:
    """Compute joint sandwich and delta-method covariance after fitting.

    Observations, or supplied clusters, must be independent. The estimators must be identified and
    satisfy their estimating equations to adequate numerical accuracy. This
    implementation excludes internal/external control functions and weighted
    survey designs. FGNLS additionally requires the original initial NLS estimates.
    """
    started = time.perf_counter()
    if data.nobs != result.nobs or int(chunk) < 1:
        raise ValueError("use the estimation sample and a positive chunk size")
    if not isinstance(cluster_correction, (bool, np.bool_)):
        raise ValueError("cluster_correction must be True or False")
    codes, n_clusters = (None, None) if clusters is None else cluster_codes(clusters, data.nobs)
    factor = n_clusters/(n_clusters-1) if n_clusters is not None and cluster_correction else 1.0
    if result.spec.control_function or result.control_function_name is not None \
            or result.selection_control_function_name is not None \
            or result.reduced_form is not None:
        raise NotImplementedError("analytical inference for ivexp/control functions is not implemented")
    if not result.converged:
        raise ValueError("analytical inference requires a converged estimate")
    if result.spec.censor:
        if selection_design is None or np.asarray(selection_design).shape != (
                data.nobs, result.np_prob-1):
            raise ValueError("supply the first-stage design on the estimation sample")
        if len(result.probits) != result.spec.neqn:
            raise ValueError("analytical inference needs every fitted Probit")
    system = _EstimatingSystem(result, data, selection_design, theta_nls,
                               int(chunk), deadline)
    bread = system.bread()
    inverse, condition = _bread_inverse(bread, system.scales)
    direct, mean_jac = system.elasticity_jacobians(_relative_step)
    total = direct + mean_jac @ system.mean_derivative()
    spec, nobs = result.spec, data.nobs
    k, nt, ne = spec.n_free, system.tau.size, total.shape[0]
    joint_map = np.zeros((k+nt, system.x.size))
    joint_map[:k, system.theta_slice] = np.eye(k)
    if nt:
        joint_map[k+system.active, np.arange(len(system.active))] = 1.0
    report_map = np.zeros((spec.n_full+nt, system.x.size))
    report_map[:spec.n_full, system.theta_slice] = delta_matrix(spec)
    if nt:
        report_map[spec.n_full+system.active, np.arange(len(system.active))] = 1.0
    reported = np.vstack([report_map, total]) if spec.censor else report_map
    raw_mean_map = np.vstack([np.zeros((len(report_map), mean_jac.shape[1])),
                              mean_jac]) if spec.censor else np.zeros(
                                  (len(report_map), mean_jac.shape[1]))
    transform = -reported @ inverse
    joint_transform = -joint_map @ inverse
    elasticity_transform = -total @ inverse
    covariance = np.zeros((len(reported), len(reported)))
    joint_covariance = np.zeros((k+nt, k+nt))
    elasticity_covariance = np.zeros((ne, ne))
    meat = np.zeros_like(bread)
    mean_values = _pack_means(result.means)
    score_mean = system.mean_score(system.x)
    if codes is not None:
        sums = [np.zeros((n_clusters, width)) for width in
                (len(reported), k+nt, ne, len(system.x))]
    for sl, scores, features in system.chunks(system.x, include_means=True):
        centered = scores - score_mean
        means_residual = features - mean_values
        influence = centered @ transform.T + means_residual @ raw_mean_map.T
        joint_influence = centered @ joint_transform.T
        elasticity_influence = (centered @ elasticity_transform.T
                                + means_residual @ mean_jac.T)
        if codes is None:
            covariance += influence.T @ influence
            joint_covariance += joint_influence.T @ joint_influence
            elasticity_covariance += elasticity_influence.T @ elasticity_influence
            meat += centered.T @ centered
        else:
            check_deadline(deadline)
            for target, values in zip(sums, (influence, joint_influence,
                                             elasticity_influence, centered)):
                np.add.at(target, codes[sl], values)
    if codes is not None:
        covariance, joint_covariance, elasticity_covariance, meat = [s.T @ s for s in sums]
    meat /= nobs
    if factor != 1.0:
        meat *= factor
    standard = np.sqrt(np.diag(meat)/nobs)
    score_size = np.divide(np.abs(score_mean), standard,
                           out=np.zeros_like(score_mean), where=standard > 0)
    covariance /= nobs*nobs
    joint_covariance /= nobs*nobs
    elasticity_covariance /= nobs*nobs
    if factor != 1.0:
        covariance *= factor
        joint_covariance *= factor
        elasticity_covariance *= factor
    if covariance.shape != result.V.shape:
        raise ValueError("analytical covariance does not match the reported vector")
    return AnalyticalInference(
        covariance, joint_covariance, elasticity_covariance, bread, meat,
        system.x.copy(), score_mean, condition, float(np.max(score_size)),
        time.perf_counter()-started, spec.neqn,
        n_clusters, "cluster" if codes is not None else "iid", factor,
    )
