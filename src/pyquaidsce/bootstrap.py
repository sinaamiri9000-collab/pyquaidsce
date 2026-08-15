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
import time
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
    """Set the BLAS thread count at run time across platforms.

    Bootstrap workers are CPU-bound in BLAS.  If each of ``n_jobs`` processes
    also starts its own pool of BLAS threads the machine is oversubscribed by a
    factor of ``n_jobs`` and everything slows to a crawl -- in testing, 2
    workers x 2 BLAS threads on 2 cores was **4.5x slower** than 2 workers x 1
    thread.  So each worker pins BLAS to a single thread and parallelism is
    taken across replications instead, which is the efficient way round.

    Returns ``True`` if a runtime knob was found and set.  Environment
    variables are also pinned in every worker before estimation starts.
    """
    try:
        from threadpoolctl import threadpool_limits

        threadpool_limits(limits=int(n))
        return True
    except Exception:
        # threadpoolctl is optional; fall back to vendor runtime symbols.
        pass

    import ctypes
    import glob

    try:
        import numpy as _np
    except Exception:  # pragma: no cover
        return False
    root = os.path.dirname(os.path.dirname(_np.__file__))
    pats = (
        "**/libscipy_openblas*.so*",
        "**/libopenblas*.so*",
        "**/libmkl_rt*.so*",
        "**/libscipy_openblas*.dll*",
        "**/libopenblas*.dll*",
        "**/mkl_rt*.dll*",
        "**/*openblas*.dylib*",
        "**/*mkl*.dylib*",
    )
    found = False
    for pat in pats:
        for path in glob.glob(os.path.join(root, pat), recursive=True)[:10]:
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
                    found = True
    return found


def _init_worker(payload, pin_blas: bool = True):
    _WORK.update(payload)
    if pin_blas:
        for var in (
            "OMP_NUM_THREADS",
            "OPENBLAS_NUM_THREADS",
            "MKL_NUM_THREADS",
            "VECLIB_MAXIMUM_THREADS",
            "NUMEXPR_NUM_THREADS",
        ):
            os.environ[var] = "1"
        set_blas_threads(1)


def _one_rep(task):
    from .estimator import quaidsce

    rep_index, rep_seed = task
    df = _WORK["df"]
    kw = _WORK["kw"]
    rep_timeout = _WORK.get("rep_timeout")
    n = len(df)
    rng = np.random.default_rng(rep_seed)
    idx = rng.integers(0, n, size=n)
    boot_df = df.iloc[idx].reset_index(drop=True)
    try:
        deadline = (
            None if rep_timeout is None
            else time.perf_counter() + float(rep_timeout)
        )
        r = quaidsce(
            boot_df, reps=0, verbose=False, _deadline=deadline, **kw
        )
        if deadline is not None and time.perf_counter() > deadline:
            raise TimeoutError(
                f"bootstrap replication exceeded {rep_timeout:g} seconds"
            )
        if not r.converged:
            return (
                rep_index,
                None,
                "RuntimeError: second-stage estimator did not converge",
            )
        return rep_index, np.asarray(r.b, dtype=float), None
    except Exception as exc:  # noqa: BLE001
        return rep_index, None, f"{type(exc).__name__}: {exc}"


def _one_rep_process(conn, task, payload, pin_blas: bool):
    """Run one replication in a disposable process for hard time limits."""
    try:
        _init_worker(payload, pin_blas=pin_blas)
        conn.send(_one_rep(task))
    except BaseException as exc:  # noqa: BLE001
        try:
            conn.send((task[0], None, f"{type(exc).__name__}: {exc}"))
        except Exception:
            pass
    finally:
        conn.close()


def _watchdog_results(ctx, tasks, payload, n_jobs: int, timeout: float):
    """Yield replications while enforcing a parent-side wall-clock timeout.

    Each active replication owns a disposable child process.  This costs more
    process-start overhead than a reusable pool, but lets the parent terminate
    a genuinely stuck BLAS/native call without killing unrelated replications.
    Cooperative deadlines inside Probit/NLSUR remain the fast first line of
    defence; this watchdog is the hard backstop.
    """
    from multiprocessing.connection import wait

    active = {}
    next_task = 0

    def launch(task):
        parent_conn, child_conn = ctx.Pipe(duplex=False)
        proc = ctx.Process(
            target=_one_rep_process,
            args=(child_conn, task, payload, n_jobs > 1),
        )
        proc.start()
        child_conn.close()
        active[task[0]] = (proc, parent_conn, time.perf_counter())

    def refill():
        nonlocal next_task
        while next_task < len(tasks) and len(active) < n_jobs:
            launch(tasks[next_task])
            next_task += 1

    refill()
    try:
        while active:
            by_conn = {record[1]: rep_index
                       for rep_index, record in active.items()}
            for conn in wait(list(by_conn), timeout=0.05):
                rep_index = by_conn[conn]
                proc, parent_conn, _ = active.pop(rep_index)
                try:
                    result = parent_conn.recv()
                except EOFError:
                    result = (
                        rep_index,
                        None,
                        f"RuntimeError: bootstrap worker exited with code "
                        f"{proc.exitcode} without returning a result",
                    )
                finally:
                    parent_conn.close()
                    proc.join(timeout=1.0)
                yield result

            now = time.perf_counter()
            expired = [
                rep_index for rep_index, (_, _, started) in active.items()
                if now - started >= timeout
            ]
            for rep_index in expired:
                proc, parent_conn, _ = active.pop(rep_index)
                proc.terminate()
                proc.join(timeout=1.0)
                if proc.is_alive() and hasattr(proc, "kill"):
                    proc.kill()
                    proc.join(timeout=1.0)
                parent_conn.close()
                yield (
                    rep_index,
                    None,
                    f"TimeoutError: bootstrap replication exceeded "
                    f"{timeout:g} seconds",
                )
            refill()
    finally:
        for proc, parent_conn, _ in active.values():
            if proc.is_alive():
                proc.terminate()
            proc.join(timeout=1.0)
            parent_conn.close()


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
    mp_context: Optional[str] = None,
    rep_timeout: Optional[float] = None,
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
    payload = {"df": df, "kw": kw, "rep_timeout": rep_timeout}

    ss = np.random.SeedSequence(seed)
    seeds = [int(x) for x in ss.generate_state(reps, dtype=np.uint32)]
    tasks = list(enumerate(seeds, 1))

    out = []
    fail_records = []

    n_jobs = max(1, min(int(n_jobs), os.cpu_count() or 1))
    started = time.perf_counter()

    def record(rep_index, b, err, completed):
        if b is None:
            fail_records.append((rep_index, err))
        else:
            out.append((rep_index, b))
        if verbose:
            elapsed = time.perf_counter() - started
            pct = 100.0 * completed / reps
            rate = completed / elapsed if elapsed > 0 else 0.0
            failed = len(fail_records)
            suffix = f", {failed} failed" if failed else ""
            print(
                f"  bootstrap {completed}/{reps} [{pct:5.1f}%] "
                f"({elapsed:.1f}s, {rate:.2f} reps/s{suffix})",
                flush=True,
            )

    if rep_timeout is not None:
        import multiprocessing as mp

        methods = mp.get_all_start_methods()
        if mp_context is None:
            context_name = "spawn" if "spawn" in methods else methods[0]
        else:
            context_name = str(mp_context).lower()
            if context_name not in methods:
                raise ValueError(
                    f"multiprocessing context {context_name!r} is not supported; "
                    f"available methods are {methods}"
                )
        ctx = mp.get_context(context_name)
        for completed, (rep_index, b, err) in enumerate(
            _watchdog_results(
                ctx, tasks, payload, n_jobs, float(rep_timeout)
            ),
            1,
        ):
            record(rep_index, b, err, completed)
    elif n_jobs == 1:
        # one worker: leave BLAS multi-threaded, it has the machine to itself
        _init_worker(payload, pin_blas=False)
        for completed, task in enumerate(tasks, 1):
            rep_index, b, err = _one_rep(task)
            record(rep_index, b, err, completed)
    else:
        import multiprocessing as mp

        methods = mp.get_all_start_methods()
        if mp_context is None:
            # ``spawn`` avoids inheriting BLAS/OpenMP locks and is available on
            # every supported platform.  Users may opt into forkserver/fork.
            context_name = "spawn" if "spawn" in methods else methods[0]
        else:
            context_name = str(mp_context).lower()
            if context_name not in methods:
                raise ValueError(
                    f"multiprocessing context {context_name!r} is not supported; "
                    f"available methods are {methods}"
                )
        ctx = mp.get_context(context_name)
        with ctx.Pool(n_jobs, initializer=_init_worker,
                      initargs=(payload, True)) as pool:
            for completed, (rep_index, b, err) in enumerate(
                pool.imap_unordered(_one_rep, tasks), 1
            ):
                record(rep_index, b, err, completed)

    if not out:
        raise RuntimeError(
            "every bootstrap replication failed; first error: "
            + (fail_records[0][1] if fail_records else "unknown")
        )
    if len(out) < 2:
        raise RuntimeError(
            "fewer than two bootstrap replications converged; standard errors "
            "cannot be computed"
        )
    out.sort(key=lambda item: item[0])
    fail_records.sort(key=lambda item: item[0])
    B = np.vstack([b for _, b in out])
    fails = [f"rep {i}: {err}" for i, err in fail_records]
    se = B.std(axis=0, ddof=1)
    return BootResult(reps_requested=reps, reps_ok=B.shape[0], b_star=B,
                      se=se, failures=fails)
