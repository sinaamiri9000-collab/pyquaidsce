"""Analytical first derivatives of the existing at-means elasticities.

Only post-estimation inference uses this module. The point elasticity function
in elasticities.py remains authoritative and unchanged. Product, chain and
quotient rules give derivatives in [active tau, free theta, sample means]
coordinates; independent finite differences check the complete result.
"""

from __future__ import annotations

import numpy as np

from .params import delta_matrix, full_slices, unpack, vech_index


def elasticity_jacobian(theta, spec, means, a0, tau, np_prob, layout, active):
    """Return the complete elasticity Jacobian in natural row-major order.

    Means follow [w, lnprices, lnexpenditure, demographics, cdf, pdf, index].
    They are differentiated independently here. Their dependence on tau is
    added separately through the sample-mean influence functions.
    """
    if spec.control_function:
        raise NotImplementedError('analytical elasticity derivatives exclude control functions')
    n, r, k = spec.neqn, spec.ndemo, spec.n_free
    active = np.asarray(active, dtype=int)
    nt = len(active)
    width = nt+k+5*n+r+1
    coefs = unpack(theta, spec)
    mapping, slices = delta_matrix(spec), full_slices(spec)

    def coefficient_derivative(name, shape):
        derivative = np.zeros((*shape, width))
        derivative[..., nt:nt+k] = mapping[slices[name]].reshape(*shape, k)
        return derivative

    da = coefficient_derivative('alpha', (n,))
    db = coefficient_derivative('beta', (n,))
    dg = np.zeros((n, n, width))
    for row, (i, j) in enumerate(vech_index(n)):
        dg[i, j, nt:nt+k] = dg[j, i, nt:nt+k] = mapping[slices['gamma']][row]
    dl = coefficient_derivative('lambda', (n,)) if spec.quadratic else np.zeros((n, width))
    dd = coefficient_derivative('delta', (n,)) if spec.censor else np.zeros((n, width))
    de = coefficient_derivative('eta', (r, n)) if r else np.zeros((0, n, width))
    dr = coefficient_derivative('rho', (r,)) if r else np.zeros((0, width))

    identity = np.eye(width)
    offset = nt+k
    dw = identity[offset:offset+n]
    dprice = identity[offset+n:offset+2*n]
    dexp = identity[offset+2*n]
    dz = identity[offset+2*n+1:offset+2*n+1+r]
    offset += 2*n+1+r
    dcdf, dpdf, dindex = (identity[offset:offset+n], identity[offset+n:offset+2*n],
                         identity[offset+2*n:offset+3*n])
    w, price, z = means.w, means.lnp, means.demo
    alpha, beta, gamma = coefs.alpha, coefs.beta, coefs.gamma
    lam, delta, eta, rho = coefs.lam, coefs.delta, coefs.eta, coefs.rho

    price_index = a0+alpha @ price+.5*price @ gamma @ price
    dprice_index = (price @ da+alpha @ dprice
                    +.5*np.einsum('i,ijk,j->k', price, dg, price)
                    +(gamma @ price) @ dprice)
    G = alpha+gamma @ price
    dG = da+np.einsum('ijk,j->ik', dg, price)+gamma @ dprice
    B = beta+z @ eta
    dB = db+eta.T @ dz+np.einsum('r,rik->ik', z, de)
    mbar = 1+rho @ z
    dmbar = z @ dr+rho @ dz
    if mbar <= 0:
        raise ValueError('elasticity derivatives require positive demographic scale')
    D = means.lnexp-price_index-np.log(mbar)
    dD = dexp-dprice_index-dmbar/mbar
    if spec.quadratic:
        q = np.exp(-B @ price)
        dq = -q*(price @ dB+B @ dprice)
        S = B+2*lam*q*D
        dS = dB+2*(q*D*dl+lam[:, None]*(D*dq+q*dD))
        T = lam*q*D**2
        dT = q*D**2*dl+lam[:, None]*(D**2*dq+2*q*D*dD)
    else:
        S, dS = B, dB
        T, dT = np.zeros(n), np.zeros((n, width))
    H = gamma-S[:, None]*G[None, :]-T[:, None]*B[None, :]
    dH = (dg-dS[:, None, :]*G[None, :, None]-S[:, None, None]*dG[None, :, :]
          -dT[:, None, :]*B[None, :, None]-T[:, None, None]*dB[None, :, :])

    if spec.censor:
        tau = np.asarray(tau)
        dtau = np.zeros((tau.size, width))
        dtau[active, np.arange(nt)] = 1
        tm, tp = np.zeros(n), np.zeros((n, n))
        dtm, dtp = np.zeros((n, width)), np.zeros((n, n, width))
        for i in range(n):
            exp_position = (layout.expenditure_position if layout is not None
                            else n if np_prob == n+r+2 else None)
            if exp_position is not None:
                tm[i] = tau[i*np_prob+exp_position]
                dtm[i] = dtau[i*np_prob+exp_position]
            for j in range(n):
                price_position = layout.price_position(j) if layout is not None else j
                if price_position is not None:
                    tp[i, j] = tau[i*np_prob+price_position]
                    dtp[i, j] = dtau[i*np_prob+price_position]
        cdf, pdf, index = means.cdf, means.pdf, means.du
        we = w*cdf+delta*pdf
        dwe = cdf[:, None]*dw+w[:, None]*dcdf+pdf[:, None]*dd+delta[:, None]*dpdf
        F = w-delta*index
        dF = dw-index[:, None]*dd-delta[:, None]*dindex
        P = pdf*F
        dP = F[:, None]*dpdf+pdf[:, None]*dF
        income_numerator = cdf*S+tm*P
        dincome_numerator = (cdf[:, None]*dS+S[:, None]*dcdf
                             +tm[:, None]*dP+P[:, None]*dtm)
        price_numerator = cdf[:, None]*H+tp*P[:, None]
        dprice_numerator = (cdf[:, None, None]*dH+H[:, :, None]*dcdf[:, None, :]
                            +tp[:, :, None]*dP[:, None, :]+P[:, None, None]*dtp)
    else:
        we, dwe = w, dw
        income_numerator, dincome_numerator = S, dS
        price_numerator, dprice_numerator = H, dH
    if not np.isfinite(we).all() or np.any(we == 0):
        raise ValueError('elasticity derivatives require finite nonzero mean shares')
    income = 1+income_numerator/we
    dincome = (dincome_numerator-(income_numerator/we)[:, None]*dwe)/we[:, None]
    dprice_elasticity = (dprice_numerator
                         -(price_numerator/we[:, None])[:, :, None]*dwe[:, None, :]
                         )/we[:, None, None]
    dcompensated = (dprice_elasticity+dincome[:, None, :]*w[None, :, None]
                    +income[:, None, None]*dw[None, :, :])
    jac = np.concatenate([dincome, dprice_elasticity.reshape(n*n, width),
                          dcompensated.reshape(n*n, width)])
    if not np.isfinite(jac).all():
        raise ValueError('nonfinite analytical elasticity derivative')
    return jac
