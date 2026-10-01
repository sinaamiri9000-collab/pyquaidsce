"""Delta-method utilities shared by model-specific result transformations."""

import numpy as np


def complex_step_jacobian(function, theta, step=1e-25):
    """Differentiate an analytic transformation without subtraction error.

    The supplied transformation must preserve complex inputs and avoid clipping
    or taking their absolute values in its returned quantities.
    """
    theta = np.asarray(theta, dtype=float)
    value = np.asarray(function(theta)).ravel()
    out = np.empty((value.size, theta.size))
    for k in range(theta.size):
        perturbed = theta.astype(complex)
        perturbed[k] += step * 1j
        out[:, k] = np.imag(np.asarray(function(perturbed)).ravel()) / step
    if not np.isfinite(out).all():
        raise ValueError("nonfinite derivatives in the delta-method transformation")
    return out
