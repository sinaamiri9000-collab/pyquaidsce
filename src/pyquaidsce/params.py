"""
Parameter mapping for the (censored) QUAIDS model.

This is a line-by-line port of the Mata routines in ``utils__quaidsce.mata``:

    _quaidsce__getcoefs_wrk   ->  unpack()
    _quaidsce__getcoefs       ->  unpack() + Coefs
    _quaidsce__fullvector     ->  full_vector()
    _quaidsce__delta          ->  delta_matrix()

Everything here is deliberately *literal*: the ordering of the free parameter
vector, the way the last good's parameters are recovered from the adding-up /
homogeneity / symmetry restrictions, and the ordering of the reported "full"
parameter vector all follow the Stata/Mata source exactly, because that is what
makes ``e(b)`` comparable element by element.

Notation
--------
n       number of goods / equations                (Mata: neqn)
R       number of demographic variables            (Mata: ndemo)
quad    True  -> quadratic (QUAIDS) term included  (Mata: quadratics == "")
        False -> AIDS, i.e. Stata's -noquadratic-
censor  True  -> Shonkwiler-Yen two-step censoring (Mata: censor == "")
        False -> Stata's -nocensor-
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Tuple

import numpy as np


# --------------------------------------------------------------------------- #
#  Specification
# --------------------------------------------------------------------------- #
@dataclass(frozen=True)
class Spec:
    """Model dimensions and switches."""

    neqn: int
    ndemo: int = 0
    quadratic: bool = True
    censor: bool = True
    control_function: bool = False

    def __post_init__(self) -> None:
        if self.neqn < 3:
            raise ValueError("must specify at least 3 expenditure shares")
        if self.censor and self.ndemo == 0:
            # quaidsce_c.ado, line ~122
            raise ValueError(
                "at least one demographic variable is needed for censoring "
                "correction"
            )
        if self.control_function and not self.censor:
            raise ValueError("control function is supported only with censor=True")

    # ---- counts ---------------------------------------------------------- #
    @property
    def n_free(self) -> int:
        """Number of *estimated* parameters, i.e. Stata's ``nparam()``."""
        n, R = self.neqn, self.ndemo
        if self.censor:
            #  alpha n, beta n, gamma n(n-1)/2, [lambda n], delta n
            k = 2 * n + n * (n - 1) // 2
            if self.quadratic:
                k += n
            k += n  # delta
        else:
            #  alpha n-1, beta n-1, gamma n(n-1)/2, [lambda n-1]
            k = 2 * (n - 1) + n * (n - 1) // 2
            if self.quadratic:
                k += n - 1
        if R > 0:
            k += R * (n - 1) + R
        if self.control_function:
            k += n
        return k

    @property
    def n_eq_estimated(self) -> int:
        """Equations handed to ``nlsur``: all n when censoring, else n-1."""
        return self.neqn if self.censor else self.neqn - 1

    @property
    def n_full(self) -> int:
        """Length of the reported (unrestricted) parameter vector."""
        n, R = self.neqn, self.ndemo
        k = 2 * n + n * (n + 1) // 2  # alpha, beta, vech(gamma)
        if self.quadratic:
            k += n
        if self.censor:
            k += n
        if self.control_function:
            k += n
        if R > 0:
            k += R * n + R
        return k

    # ---- names ----------------------------------------------------------- #
    def full_names(self, demo_names: List[str] | None = None) -> List[str]:
        """Reported coefficient names, exactly as quaidsce_c.ado's namestripe."""
        n, R = self.neqn, self.ndemo
        demo_names = list(demo_names or [f"z{r + 1}" for r in range(R)])
        out = [f"alpha:alpha_{i + 1}" for i in range(n)]
        out += [f"beta:beta_{i + 1}" for i in range(n)]
        for j in range(1, n + 1):  # column-major vech
            for i in range(j, n + 1):
                out.append(f"gamma:gamma_{i}_{j}")
        if self.quadratic:
            out += [f"lambda:lambda_{i + 1}" for i in range(n)]
        if self.censor:
            out += [f"delta:delta_{i + 1}" for i in range(n)]
        if self.control_function:
            out += [f"cfcoef:cfcoef_{i + 1}" for i in range(n)]
        if R > 0:
            for v in demo_names:
                out += [f"eta:eta_{v}_{i + 1}" for i in range(n)]
            out += [f"rho:rho_{v}" for v in demo_names]
        return out

    def tau_names(self, demo_names: List[str] | None = None) -> List[str]:
        """Probit coefficient names: p1_i .. pn_i, M_i, <demos>_i, cons_i."""
        n, R = self.neqn, self.ndemo
        demo_names = list(demo_names or [f"z{r + 1}" for r in range(R)])
        out: List[str] = []
        for i in range(1, n + 1):
            out += [f"tau:p{j}_{i}" for j in range(1, n + 1)]
            out.append(f"tau:M_{i}")
            out += [f"tau:{v}_{i}" for v in demo_names]
            out.append(f"tau:cons_{i}")
        return out

    def elas_names(self) -> List[str]:
        n = self.neqn
        out = [f"ELAS_INC:e_{i + 1}" for i in range(n)]
        for j in range(1, n + 1):
            for i in range(1, n + 1):
                out.append(f"ELAS_UNCOMP:e_{i}_{j}")
        for j in range(1, n + 1):
            for i in range(1, n + 1):
                out.append(f"ELAS_COMP:e_{i}_{j}")
        return out


# --------------------------------------------------------------------------- #
#  Coefficient container
# --------------------------------------------------------------------------- #
@dataclass
class Coefs:
    alpha: np.ndarray  # (n,)
    beta: np.ndarray  # (n,)
    gamma: np.ndarray  # (n, n) symmetric
    lam: np.ndarray  # (n,)   zeros when noquadratic
    delta: np.ndarray  # (n,)   ones  when nocensor
    eta: np.ndarray  # (R, n)  empty when R == 0
    rho: np.ndarray  # (R,)    empty when R == 0
    # Keep a default so the public seven-argument constructor from version
    # 1.0.1 remains valid. Internal code always passes this block explicitly.
    cfcoef: Optional[np.ndarray] = None

    def __post_init__(self) -> None:
        if self.cfcoef is None:
            self.cfcoef = np.zeros_like(self.alpha, dtype=float)


# --------------------------------------------------------------------------- #
#  getcoefs_wrk
# --------------------------------------------------------------------------- #
def unpack(theta: np.ndarray, spec: Spec) -> Coefs:
    """Free parameter vector -> structural coefficients.

    Direct port of ``_quaidsce__getcoefs_wrk``.  Note the two Mata quirks that
    are reproduced faithfully:

    * ``alpha`` is initialised to **1** and ``beta``/``lambda`` to **0**, so that
      under ``nocensor`` the n-th good's parameters satisfy the adding-up
      restrictions sum(alpha)=1, sum(beta)=0, sum(lambda)=0.
    * When censoring **is** used the n-th element is *overwritten* by a free
      parameter, i.e. adding-up is **not** imposed on alpha/beta/lambda (it
      cannot be, because the Shonkwiler-Yen transformation destroys it).  It
      *is* still imposed on ``eta``.
    """
    theta = np.asarray(theta, dtype=float).ravel()
    n, R = spec.neqn, spec.ndemo
    if theta.size != spec.n_free:
        raise ValueError(
            f"expected {spec.n_free} free parameters, got {theta.size}"
        )
    col = 0

    alpha = np.ones(n)
    for i in range(n - 1):
        alpha[i] = theta[col]
        alpha[n - 1] -= alpha[i]
        col += 1
    if spec.censor:
        alpha[n - 1] = theta[col]
        col += 1

    beta = np.zeros(n)
    for i in range(n - 1):
        beta[i] = theta[col]
        beta[n - 1] -= beta[i]
        col += 1
    if spec.censor:
        beta[n - 1] = theta[col]
        col += 1

    gamma = np.zeros((n, n))
    # j outer so that the ordering matches vech()/invvech()
    for j in range(n - 1):
        for i in range(j, n - 1):
            gamma[i, j] = theta[col]
            if i != j:
                gamma[j, i] = theta[col]
            col += 1
    for i in range(n - 1):
        for j in range(n - 1):
            gamma[i, n - 1] -= gamma[i, j]
        gamma[n - 1, i] = gamma[i, n - 1]
    for i in range(n - 1):
        gamma[n - 1, n - 1] -= gamma[i, n - 1]

    lam = np.zeros(n)
    if spec.quadratic:
        for i in range(n - 1):
            lam[i] = theta[col]
            lam[n - 1] -= lam[i]
            col += 1
        if spec.censor:
            lam[n - 1] = theta[col]
            col += 1

    delta = np.ones(n)
    if spec.censor:
        for i in range(n):
            delta[i] = theta[col]
            col += 1

    cfcoef = np.zeros(n)
    if spec.control_function:
        cfcoef[:] = theta[col:col + n]
        col += n

    eta = np.zeros((R, n))
    rho = np.zeros(R)
    if R > 0:
        for i in range(R):
            for j in range(n - 1):
                eta[i, j] = theta[col]
                eta[i, n - 1] -= eta[i, j]
                col += 1
        for i in range(R):
            rho[i] = theta[col]
            col += 1

    assert col == theta.size, (col, theta.size)
    return Coefs(alpha, beta, gamma, lam, delta, eta, rho, cfcoef)


# --------------------------------------------------------------------------- #
#  fullvector
# --------------------------------------------------------------------------- #
def vech(a: np.ndarray) -> np.ndarray:
    """Mata's vech(): column-major lower triangle of a square matrix."""
    n = a.shape[0]
    return np.concatenate([a[j:, j] for j in range(n)])


def full_vector(theta: np.ndarray, spec: Spec) -> np.ndarray:
    """Free vector -> reported vector, i.e. ``_quaidsce__fullvector``.

    Order: alpha, beta, vech(gamma), [lambda], [delta], [vec(eta')], [rho].
    """
    c = unpack(theta, spec)
    parts = [c.alpha, c.beta, vech(c.gamma)]
    if spec.quadratic:
        parts.append(c.lam)
    if spec.censor:
        parts.append(c.delta)
    if spec.control_function:
        parts.append(c.cfcoef)
    if spec.ndemo > 0:
        parts.append(c.eta.reshape(-1))  # vec(eta')' == row-major eta
        parts.append(c.rho)
    return np.concatenate(parts)


# --------------------------------------------------------------------------- #
#  delta matrix
# --------------------------------------------------------------------------- #
def _blockdiag(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    out = np.zeros((a.shape[0] + b.shape[0], a.shape[1] + b.shape[1]))
    out[: a.shape[0], : a.shape[1]] = a
    out[a.shape[0]:, a.shape[1]:] = b
    return out


def delta_matrix(spec: Spec) -> np.ndarray:
    """Jacobian d(full vector)/d(free vector); port of ``_quaidsce__delta``.

    The map is affine, so this constant matrix is both the delta-method
    transform for the covariance (``Vfull = Delta V Delta'``) and the chain-rule
    factor for the model Jacobian (``J_free = J_full @ Delta``).
    """
    ng, R = spec.neqn, spec.ndemo
    ngm1 = ng - 1

    if spec.censor:
        block = np.eye(ng)
    else:
        block = np.vstack([np.eye(ngm1), -np.ones((1, ngm1))])
    blockd = np.vstack([np.eye(ngm1), -np.ones((1, ngm1))])

    Delta = block
    Delta = _blockdiag(Delta, block)  # beta

    # ---- Gamma ---------------------------------------------------------- #
    Gamma = np.zeros(((ng + 1) * ng // 2, (ngm1 + 1) * ngm1 // 2))
    m = 0
    for j in range(1, ng + 1):
        for i in range(j, ng + 1):
            nn = 0
            for jc in range(1, ng):
                for ic in range(jc, ng):
                    if j < ng and i < ng:  # Case 1
                        if jc == j and ic == i:
                            Gamma[m, nn] = 1
                    elif j < ng and i == ng:  # Case 2
                        if jc == j or ic == j:
                            Gamma[m, nn] = -1
                    elif j == ng and i == ng:  # Case 3
                        Gamma[m, nn] += 1
                        if ic != jc:
                            Gamma[m, nn] += 1
                    nn += 1
            m += 1
    Delta = _blockdiag(Delta, Gamma)

    if spec.quadratic:
        Delta = _blockdiag(Delta, block)
    if spec.censor:
        Delta = _blockdiag(Delta, block)
    if spec.control_function:
        Delta = _blockdiag(Delta, np.eye(ng))
    if R > 0:
        for _ in range(R):
            Delta = _blockdiag(Delta, blockd)
        Delta = _blockdiag(Delta, np.eye(R))
    return Delta


# --------------------------------------------------------------------------- #
#  index bookkeeping for the *full* vector (used by the Jacobian)
# --------------------------------------------------------------------------- #
def full_slices(spec: Spec) -> dict:
    """Slices of the full parameter vector, keyed by block name."""
    n, R = spec.neqn, spec.ndemo
    out, p = {}, 0
    out["alpha"] = slice(p, p + n); p += n
    out["beta"] = slice(p, p + n); p += n
    ng2 = n * (n + 1) // 2
    out["gamma"] = slice(p, p + ng2); p += ng2
    if spec.quadratic:
        out["lambda"] = slice(p, p + n); p += n
    if spec.censor:
        out["delta"] = slice(p, p + n); p += n
    if spec.control_function:
        out["cfcoef"] = slice(p, p + n); p += n
    if R > 0:
        out["eta"] = slice(p, p + R * n); p += R * n
        out["rho"] = slice(p, p + R); p += R
    assert p == spec.n_full
    return out


def vech_index(n: int) -> List[Tuple[int, int]]:
    """(i, j) 0-based row/col pairs in vech order (column-major lower tri)."""
    return [(i, j) for j in range(n) for i in range(j, n)]


def free_slices(spec: Spec) -> dict:
    """Slices of the *free* parameter vector, keyed by block name."""
    n, R = spec.neqn, spec.ndemo
    w = n if spec.censor else n - 1  # width of the alpha/beta/lambda blocks
    out, p = {}, 0
    out["alpha"] = slice(p, p + w); p += w
    out["beta"] = slice(p, p + w); p += w
    ng2 = n * (n - 1) // 2
    out["gamma"] = slice(p, p + ng2); p += ng2
    if spec.quadratic:
        out["lambda"] = slice(p, p + w); p += w
    if spec.censor:
        out["delta"] = slice(p, p + n); p += n
    if spec.control_function:
        out["cfcoef"] = slice(p, p + n); p += n
    if R > 0:
        out["eta"] = slice(p, p + R * (n - 1)); p += R * (n - 1)
        out["rho"] = slice(p, p + R); p += R
    assert p == spec.n_free, (p, spec.n_free)
    return out


def delta_blocks(spec: Spec):
    """Block decomposition of :func:`delta_matrix`.

    Returns a list of ``(full_slice, free_slice, matrix)`` triples.  Because the
    delta matrix is block diagonal, applying it block by block turns an
    ``O(n_full * n_free)`` dense product into a handful of tiny ones -- which
    matters a lot inside the Gauss-Newton loop.
    """
    n, R = spec.neqn, spec.ndemo
    ngm1 = n - 1
    fs, xs = full_slices(spec), free_slices(spec)
    block = np.eye(n) if spec.censor else np.vstack(
        [np.eye(ngm1), -np.ones((1, ngm1))]
    )
    blockd = np.vstack([np.eye(ngm1), -np.ones((1, ngm1))])

    D = delta_matrix(spec)
    gfull, gfree = fs["gamma"], xs["gamma"]
    Gamma = D[gfull, gfree]

    out = [
        (fs["alpha"], xs["alpha"], block),
        (fs["beta"], xs["beta"], block),
        (gfull, gfree, Gamma),
    ]
    if spec.quadratic:
        out.append((fs["lambda"], xs["lambda"], block))
    if spec.censor:
        out.append((fs["delta"], xs["delta"], np.eye(n)))
    if spec.control_function:
        out.append((fs["cfcoef"], xs["cfcoef"], np.eye(n)))
    if R > 0:
        e_full, e_free = fs["eta"].start, xs["eta"].start
        for r in range(R):
            out.append(
                (slice(e_full + r * n, e_full + (r + 1) * n),
                 slice(e_free + r * ngm1, e_free + (r + 1) * ngm1),
                 blockd)
            )
        out.append((fs["rho"], xs["rho"], np.eye(R)))
    return out
