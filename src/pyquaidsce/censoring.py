"""Model-independent Shonkwiler--Yen means and input derivatives.

Latent demand formulas remain in their model backends. This composition also
accepts control-function-augmented latent means with their residual held fixed.
"""

import numpy as np
from scipy.special import ndtr


def normal_terms(index):
    """CDF/PDF of the Probit xb index; preserve complex-step inputs."""
    index = np.asarray(index)
    return ndtr(index), np.exp(-0.5 * index**2) / np.sqrt(2 * np.pi)


def corrected_mean(latent, cdf, pdf, delta):
    return cdf * latent + pdf * np.asarray(delta)[None, :]


def corrected_derivatives(latent, dm, dp, index, delta, index_dm, index_dp):
    """Log-expenditure/log-price derivatives, shapes (N,G) and (N,G,G)."""
    cdf, pdf = normal_terms(index)
    selection = pdf * (latent - np.asarray(delta)[None, :] * index)
    return (cdf * dm + selection * np.asarray(index_dm)[None, :],
            cdf[:, :, None] * dp
            + selection[:, :, None] * np.asarray(index_dp)[None, :, :])


def index_coefficients(tau, layout, neqn):
    """Use explicit positions, with zero derivatives for omitted inputs."""
    coefficients = np.asarray(tau).reshape(neqn, layout.width)
    dm = np.zeros(neqn, dtype=coefficients.dtype)
    dp = np.zeros((neqn, neqn), dtype=coefficients.dtype)
    if layout.expenditure_position is not None:
        dm = coefficients[:, layout.expenditure_position]
    for j in range(neqn):
        pos = layout.price_position(j)
        if pos is not None:
            dp[:, j] = coefficients[:, pos]
    return dm, dp
