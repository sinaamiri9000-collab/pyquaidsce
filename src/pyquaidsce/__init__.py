"""
pyquaidsce — a faithful Python port of the Stata package ``quaidsce`` v2.0.

    Caro, J. C., Melo, G., Molina, J. A. and Salgado, J. C. (2025)
    "quaidsce: Censored QUAIDS estimation", https://github.com/juancaros/quaidsce

The model is Banks, Blundell & Lewbel's (1997) quadratic AIDS with Ray's (1983)
demographic scaling (Poi 2012), estimated by nonlinear SUR after the
Shonkwiler & Yen (1999) two-step correction for zero expenditure shares.

Quick start
-----------
>>> from pyquaidsce import quaidsce
>>> res = quaidsce(df,
...                shares=["w1", "w2", "w3", "w4"],
...                prices=["p1", "p2", "p3", "p4"],
...                expenditure="total",
...                demographics=["income"],
...                anot=1.6, reps=200)
>>> print(res.summary())
>>> print(res.elasticity_tables())
"""

from .elasticities import (Elasticities, elasticities,
                           fitted_share_derivatives, sample_means)
from .estimator import first_stage, quaidsce
from .model import (DemandData, augmented_latent_shares, fitted_shares,
                    jacobian_full, latent_shares)
from .nlsur import nlsur
from .params import Coefs, Spec, delta_matrix, full_vector, unpack
from .probit import probit
from .results import QuaidsceResults
from .selection import FirstStageLayout

__version__ = "1.2.0"

__all__ = [
    "quaidsce",
    "first_stage",
    "nlsur",
    "probit",
    "Spec",
    "Coefs",
    "DemandData",
    "FirstStageLayout",
    "QuaidsceResults",
    "Elasticities",
    "elasticities",
    "fitted_share_derivatives",
    "sample_means",
    "unpack",
    "full_vector",
    "delta_matrix",
    "fitted_shares",
    "latent_shares",
    "augmented_latent_shares",
    "jacobian_full",
    "__version__",
]
