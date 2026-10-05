def gauss_newton(theta0: np.ndarray, d: DemandData, spec: Spec, sigma: np.ndarray, tol: float=1e-13, max_iter: int=200, chunk: int=2000, verbose: bool=False, say=None, algorithm: str='gn', nrtol_stop: float=1e-12, stop_rule: str='standard', deadline: Optional[float]=None) -> Tuple[np.ndarray, float, int, bool]:
    """Minimise ``sum_t u_t' sigma^-1 u_t``.

    ``algorithm="lm"`` is a Levenberg-Marquardt trust-region method
    with Nielsen's damping update, working in coordinates scaled by
    ``diag(J'J)^-1/2``.  ``algorithm="gn"`` (the default) is plain
    Gauss-Newton with Hartley
    step halving, i.e. the scheme Stata's ``nl``/``nlsur`` use; it is kept as the
    reference behaviour, but on a large censored system it creeps along the flat
    valley of the criterion while LM turns the corner.  Both target the same
    stationary point.
    """
    cache = make_cache(spec)
    P = _whitener(sigma)
    theta = np.asarray(theta0, dtype=float).copy()
    converged = False
    it = 0
    mu = 0.001
    nu = 2.0
    for it in range(1, max_iter + 1):
        check_deadline(deadline)
        G, g, obj = _normal_equations(theta, d, spec, cache, P, chunk, deadline=deadline)
        if algorithm == 'gn':
            direction = _solve_scaled(G, g)
            nrtol = abs(float(direction @ g)) / max(abs(obj), 1e-300)
            best, t = (None, 1.0)
            for _ in range(40):
                check_deadline(deadline)
                cand = theta + t * direction
                oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                if oc < obj:
                    best = (cand, oc, t, 0.0)
                    break
                t *= 0.5
            if best is None:
                m2 = 1e-08
                for _ in range(30):
                    check_deadline(deadline)
                    cand = theta + _solve_scaled(G, g, mu=m2)
                    oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                    if oc < obj:
                        best = (cand, oc, 1.0, m2)
                        break
                    m2 *= 10.0
            if best is None:
                cutoff = 1e-05 if stop_rule == 'standard' else nrtol_stop
                converged = bool(nrtol < cutoff)
                break
            cand, oc, t, m2 = best
        else:
            gn = _solve_scaled(G, g)
            nrtol = abs(float(gn @ g)) / max(abs(obj), 1e-300)
            accepted = False
            for _ in range(50):
                check_deadline(deadline)
                step = _solve_scaled(G, g, mu=mu)
                cand = theta + step
                oc = _safe_obj(cand, d, spec, P, deadline=deadline)
                pred = 2.0 * float(step @ g) - float(step @ (G @ step))
                rho = (obj - oc) / pred if pred > 0 else -1.0
                if oc < obj and rho > 0:
                    accepted = True
                    if rho > 0.75:
                        mu = max(mu / 3.0, 1e-14)
                        nu = 2.0
                    elif rho < 0.25:
                        mu = min(mu * 2.0, 1000000000000.0)
                    break
                mu = min(mu * nu, 100000000000000.0)
                nu *= 2.0
            if not accepted:
                cutoff = 1e-05 if stop_rule == 'standard' else nrtol_stop
                converged = bool(nrtol < cutoff)
                break
            t, m2 = (1.0, mu)
        rel = (obj - oc) / max(abs(obj), 1e-300)
        smax = float(np.max(np.abs(cand - theta) / (np.abs(theta) + 1e-06)))
        mreldif = float(np.max(np.abs(cand - theta) / (np.abs(theta) + 1.0)))
        criterion_param, criterion_rss = observe_accepted(theta, cand, obj, oc, it, mreldif, nrtol, t, m2)
        theta, obj = (cand, oc)
        if say is not None:
            say(f'    GN {it:3d}  obj={obj:.12g}  rel={rel:.3e} step={smax:.2e} t={t:g} mu={m2:.2e} nrtol={nrtol:.2e}')
        if stop_rule == 'standard':
            if mreldif < 1e-05 or rel < 1e-07 or nrtol < 1e-05:
                converged = True
                break
        elif nrtol < nrtol_stop or (rel < tol and smax < 1e-09):
            converged = True
            break
    return (theta, obj, it, converged)
