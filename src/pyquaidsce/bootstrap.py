"""
Nonparametric bootstrap corresponding to
``parallel bs, reps(#): quaidsce_c ...``.

``quaidsce.ado`` is nothing but a bootstrap wrapper around ``quaidsce_c``::

    parallel bs, reps(`reps'): quaidsce_c `0'

Stata's ``bootstrap`` draws ``N`` observations with replacement from the
estimation sample, re-runs the *whole* two-step procedure (probits included --
which is the entire point, since that is what makes the second-stage standard
errors valid) and reports the standard deviation of the replicates together with
normal-based confidence intervals.

Replicates that fail (a share with no zeros in the resample, a probit that will
not converge, a singular residual covariance, or second-stage nonconvergence)
are dropped, as Stata drops replications that error out. Python and Stata use
different random-number generators, so equal seeds do not imply identical
resamples.
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field
from typing import List, Optional

import numpy as np


@dataclass
class BootResult:
    reps_requested: int
    reps_ok: int
    b_star: np.ndarray  # (reps_ok, k)
    se: np.ndarray  # (k,)
    failures: List[str] = field(default_factory=list)

    @property
    def V(self) -> np.ndarray:
        return np.cov(self.b_star, rowvar=False, ddof=1)

    def ci(self, b: np.ndarray, level: float = 95.0):
        from scipy.stats import norm

        z = norm.ppf(0.5 + level / 200.0)
        return b - z * self.se, b + z * self.se

    def percentile_ci(self, level: float = 95.0):
        lo = (100.0 - level) / 2.0
        return (
            np.percentile(self.b_star, lo, axis=0),
            np.percentile(self.b_star, 100.0 - lo, axis=0),
        )


_WORK = {}


def set_blas_threads(n: int) -> bool:
    """Set the BLAS thread count at run time.

    Bootstrap workers are CPU-bound in BLAS.  If each of ``n_jobs`` processes
    also starts its own pool of BLAS threads the machine is oversubscribed by a
    factor of ``n_jobs`` and everything slows to a crawl -- in testing, 2
    workers x 2 BLAS threads on 2 cores was **4.5x slower** than 2 workers x 1
    thread.  So each worker pins BLAS to a single thread and parallelism is
    taken across replications instead, which is the efficient way round.

    Returns True if the runtime knob was found.  If it was not, export
    ``OMP_NUM_THREADS=1`` (and ``OPENBLAS_NUM_THREADS=1``) before starting
    Python.
    """
    import ctypes
    import glob

    try:
        import numpy as _np
    except Exception:  # pragma: no cover
        return False
    root = os.path.dirname(os.path.dirname(_np.__file__))
    pats = ("**/libscipy_openblas*.so*", "**/libopenblas*.so*", "**/libmkl_rt*.so*")
    for pat in pats:
        for path in glob.glob(os.path.join(root, pat), recursive=True)[:6]:
            try:
                lib = ctypes.CDLL(path)
            except OSError:
                continue
            for nm in ("scipy_openblas_set_num_threads",
                       "openblas_set_num_threads",
                       "goto_set_num_threads",
                       "MKL_Set_Num_Threads"):
                fn = getattr(lib, nm, None)
                if fn is not None:
                    fn(ctypes.c_int(int(n)))
                    return True
    return False


def _init_worker(payload, pin_blas: bool = True):
    _WORK.update(payload)
    if pin_blas:
        for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                    "MKL_NUM_THREADS"):
            os.environ[var] = "1"
        set_blas_threads(1)


def _one_rep(rep_seed: int):
    from .estimator import quaidsce

    df = _WORK["df"]
    kw = _WORK["kw"]
    n = len(df)
    rng = np.random.default_rng(rep_seed)
    idx = rng.integers(0, n, size=n)
    boot_df = df.iloc[idx].reset_index(drop=True)
    try:
        r = quaidsce(boot_df, reps=0, verbose=False, **kw)
        if not r.converged:
            return None, "RuntimeError: second-stage estimator did not converge"
        return np.asarray(r.b, dtype=float), None
    except Exception as exc:  # noqa: BLE001
        return None, f"{type(exc).__name__}: {exc}"


def bootstrap(
    *,
    data,
    shares,
    prices,
    lnprices,
    expenditure,
    lnexpenditure,
    demographics,
    anot,
    quadratic,
    censor,
    method,
    initial,
    sigma_initial,
    first_stage_predict,
    strict_stata,
    vce_sigma,
    algorithm,
    stop_rule,
    tol: float,
    max_outer: int,
    max_iter: int,
    chunk: int,
    nrtol_stop: float,
    inner_nrtol_early: float,
    sigma_tol: float,
    bootstrap_start: str,
    reps: int,
    seed: Optional[int],
    n_jobs: int,
    touse: np.ndarray,
    verbose: bool = True,
) -> BootResult:
    df = data.loc[touse].reset_index(drop=True) if hasattr(data, "loc") else data
    kw = dict(
        shares=list(shares),
        prices=None if prices is None else list(prices),
        lnprices=None if lnprices is None else list(lnprices),
        expenditure=expenditure,
        lnexpenditure=lnexpenditure,
        demographics=None if demographics is None else list(demographics),
        anot=anot,
        quadratic=quadratic,
        censor=censor,
        method=method,
        initial=(np.asarray(initial, float)
                 if bootstrap_start == "warm" else None),
        sigma_initial=(np.asarray(sigma_initial, float)
                       if bootstrap_start == "warm" else None),
        first_stage_predict=first_stage_predict,
        strict_stata=strict_stata,
        vce_sigma=vce_sigma,
        algorithm=algorithm,
        stop_rule=stop_rule,
        tol=tol,
        max_outer=max_outer,
        max_iter=max_iter,
        chunk=chunk,
        nrtol_stop=nrtol_stop,
        inner_nrtol_early=inner_nrtol_early,
        sigma_tol=sigma_tol,
    )
    payload = {"df": df, "kw": kw}

    ss = np.random.SeedSequence(seed)
    seeds = [int(x) for x in ss.generate_state(reps, dtype=np.uint32)]

    out: List[np.ndarray] = []
    fails: List[str] = []

    n_jobs = max(1, min(int(n_jobs), os.cpu_count() or 1))
    if n_jobs == 1:
        # one worker: leave BLAS multi-threaded, it has the machine to itself
        _init_worker(payload, pin_blas=False)
        for i, s in enumerate(seeds, 1):
            b, err = _one_rep(s)
            if b is None:
                fails.append(err)
            else:
                out.append(b)
            if verbose and (i % max(1, reps // 20) == 0 or i == reps):
                print(f"  bootstrap {i}/{reps} ({len(fails)} failed)")
    else:
        import multiprocessing as mp

        methods = mp.get_all_start_methods()
        ctx = mp.get_context("fork" if "fork" in methods else "spawn")
        with ctx.Pool(n_jobs, initializer=_init_worker,
                      initargs=(payload, True)) as pool:
            for i, (b, err) in enumerate(
                pool.imap(_one_rep, seeds), 1
            ):
                if b is None:
                    fails.append(err)
                else:
                    out.append(b)
                if verbose and (i % max(1, reps // 20) == 0 or i == reps):
                    print(f"  bootstrap {i}/{reps} ({len(fails)} failed)")

    if not out:
        raise RuntimeError(
            "every bootstrap replication failed; first error: "
            + (fails[0] if fails else "unknown")
        )
    if len(out) < 2:
        raise RuntimeError(
            "fewer than two bootstrap replications converged; standard errors "
            "cannot be computed"
        )
    B = np.vstack(out)
    se = B.std(axis=0, ddof=1)
    return BootResult(reps_requested=reps, reps_ok=B.shape[0], b_star=B,
                      se=se, failures=fails)
