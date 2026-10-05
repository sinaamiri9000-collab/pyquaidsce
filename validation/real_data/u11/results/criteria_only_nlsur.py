@with_blas_threads
def nlsur(d: DemandData, spec: Spec, theta0: Optional[np.ndarray]=None, sigma0: Optional[np.ndarray]=None, start: str='zero', method: str='ifgnls', tol: float=1e-13, sigma_tol: float=1e-05, max_outer: int=200, max_iter: int=200, chunk: int=2000, nrtol_stop: float=1e-12, inner_nrtol_early: float=1e-08, stop_rule: str='standard', vce_sigma: str='objective', algorithm: str='gn', blas_threads: Optional[int]=1, verbose: bool=False, log: Optional[Callable[[str], None]]=None, gn_log: Optional[Callable[[str], None]]=None, deadline: Optional[float]=None) -> NlsurResult:
    """Fit the system.

    Parameters
    ----------
    method : {"nls", "fgnls", "ifgnls"}
        Estimation method. ``"ifgnls"`` is the pyquaidsce default.
    blas_threads : positive integer or None
        Temporarily limits BLAS threads during this fit. ``None`` leaves the
        caller's BLAS runtime unchanged.
    vce_sigma : {"objective", "final"}
        Which ``Sigma_hat`` enters ``e(V)``.  ``"objective"`` is the one used in
        the last minimisation (textbook FGNLS, and what Stata reports);
        ``"final"`` recomputes it from the final residuals.  The two coincide
        for ``ifgnls`` at convergence.
    """
    method = method.lower()
    if method not in ('nls', 'fgnls', 'ifgnls'):
        raise ValueError(f'unknown method {method!r}')
    if algorithm not in ('gn', 'lm'):
        raise ValueError(f'unknown algorithm {algorithm!r}')
    if start not in ('zero', 'linear'):
        raise ValueError(f'unknown start {start!r}')
    if stop_rule not in ('tight', 'standard'):
        raise ValueError(f'unknown stop_rule {stop_rule!r}')
    if vce_sigma not in ('objective', 'final'):
        raise ValueError(f'unknown vce_sigma {vce_sigma!r}')
    say = log or (print if verbose else lambda *_: None)
    check_deadline(deadline)
    N = d.nobs
    m = spec.n_eq_estimated
    K = spec.n_free
    if theta0 is not None:
        theta = np.asarray(theta0, float).copy()
    elif start == 'linear':
        from .start import linear_start
        theta = linear_start(d, spec)
    else:
        theta = np.zeros(K)

    def sigma_from(th: np.ndarray) -> np.ndarray:
        u = residuals(th, d, spec)
        return u.T @ u / N
    history: List[float] = []
    total_gn = 0
    if sigma0 is not None and theta0 is not None and (method == 'ifgnls'):
        sigma_obj = np.asarray(sigma0, dtype=float)
        history.append(_objective(theta, d, spec, _whitener(sigma_obj)))
        n_outer = 2
        tight = False
        outer_converged = False
        ok = False
        obj = history[-1]
        for outer in range(3, max_outer + 1):
            check_deadline(deadline)
            say(f'FGNLS iteration {outer}...')
            sigma_new = sigma_from(theta)
            theta_prev = theta.copy()
            theta, obj, ngn, ok = gauss_newton(theta, d, spec, sigma_new, tol=tol, max_iter=max_iter, chunk=chunk, say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop if tight else max(nrtol_stop, inner_nrtol_early), stop_rule=stop_rule, deadline=deadline)
            total_gn += ngn
            history.append(obj)
            sigma_obj = sigma_new
            n_outer = outer
            rel = float(np.max(np.abs(theta - theta_prev) / (np.abs(theta_prev) + 1e-08)))
            if observe_outer(theta_prev, theta, sigma_new, outer, ok):
                outer_converged = True
                break
        return _finish(theta, d, spec, sigma_obj, sigma_from, method, vce_sigma, chunk, n_outer, total_gn, bool(ok and outer_converged), obj, history, deadline=deadline)
    say('Calculating NLS estimates...')
    I_m = np.eye(m)
    theta, obj, ngn, ok = gauss_newton(theta, d, spec, I_m, tol=tol, max_iter=max_iter, chunk=chunk, say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop, stop_rule=stop_rule, deadline=deadline)
    total_gn += ngn
    history.append(obj)
    sigma_obj = I_m
    n_outer = 1
    if method != 'nls':
        say('Calculating FGNLS estimates...')
        sigma = sigma_from(theta)
        theta_prev = theta.copy()
        theta, obj, ngn, ok = gauss_newton(theta, d, spec, sigma, tol=tol, max_iter=max_iter, chunk=chunk, say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop, stop_rule=stop_rule, deadline=deadline)
        total_gn += ngn
        history.append(obj)
        sigma_obj = sigma
        n_outer = 2
        if method == 'ifgnls':
            tight = False
            outer_converged = False
            for outer in range(3, max_outer + 1):
                check_deadline(deadline)
                say(f'FGNLS iteration {outer}...')
                sigma_new = sigma_from(theta)
                theta_prev = theta.copy()
                theta, obj, ngn, ok = gauss_newton(theta, d, spec, sigma_new, tol=tol, max_iter=max_iter, chunk=chunk, say=gn_log, algorithm=algorithm, nrtol_stop=nrtol_stop if tight else max(nrtol_stop, inner_nrtol_early), stop_rule=stop_rule, deadline=deadline)
                total_gn += ngn
                history.append(obj)
                sigma_obj = sigma_new
                n_outer = outer
                rel = float(np.max(np.abs(theta - theta_prev) / (np.abs(theta_prev) + 1e-08)))
                if observe_outer(theta_prev, theta, sigma_new, outer, ok):
                    outer_converged = True
                    break
    overall_converged = bool(ok)
    if method == 'ifgnls':
        overall_converged = bool(ok and outer_converged)
    return _finish(theta, d, spec, sigma_obj, sigma_from, method, vce_sigma, chunk, n_outer, total_gn, overall_converged, obj, history, deadline=deadline)
