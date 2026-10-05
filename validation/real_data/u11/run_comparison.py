"""Compare U11 stopping predicates with identical pyquaidsce 1.6.0 numerics.

Both runs start at zero, include original NLS and FGNLS initialization, and use
the same prepared demand data / already estimated Probits from the C11 inputs.
AST copies add observations; the second changes convergence predicates only.
No production source file or scientific preparation code is edited.
"""
from pathlib import Path
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS',
            'NUMEXPR_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    os.environ[key] = '1'
import sys, ast, inspect, copy, importlib, time, pickle, json, hashlib, difflib, argparse
from unittest.mock import patch
ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/pyquaidsce')
C11 = Path('/workspace/real-data-c11')
sys.path[:0] = [str(REPO / 'src'), str(C11 / 'D1_C11_FINAL_PSU_BOOTSTRAP_V2')]
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
solver = importlib.import_module('pyquaidsce.nlsur')
OUT = ROOT / 'results'
OUT.mkdir(exist_ok=True)
EPS, TAU, COVTOL = 1e-5, 1e-3, 1e-10
STATE = {}

def save(name, obj):
    def fallback(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x).__name__)
    (OUT / name).write_text(json.dumps(obj, default=fallback, indent=2, allow_nan=False) + '\n')

def observe_accepted(old, new, oldQ, newQ, it, mreldif, nrtol, fraction, damping):
    ratios = np.abs(new - old) / (EPS * (np.abs(old) + TAU))
    param = bool(np.all(np.abs(new - old) <= EPS * (np.abs(old) + TAU)))
    rss = bool(abs(newQ - oldQ) <= EPS * (oldQ + TAU))
    idx = int(np.argmax(ratios))
    row = dict(stage=STATE['stage'], iteration=it, Q_old=float(oldQ), Q=float(newQ),
               param_ratio=float(np.max(ratios)), rss_ratio=float(abs(newQ-oldQ)/(EPS*(oldQ+TAU))),
               param_pass=param, rss_pass=rss, worst_theta_index=idx,
               worst_theta_component=STATE['labels'][idx],
               worst_old=float(old[idx]), worst_new=float(new[idx]),
               theta_step_max_abs=float(np.max(np.abs(new-old))),
               standard_mreldif=float(mreldif), nrtol=float(nrtol),
               accepted_fraction=float(fraction), damping=float(damping))
    STATE['steps'].append(row)
    STATE['accepted_thetas'].append(new.copy())
    if it == 1 or it % 10 == 0 or (param and rss):
        print(f"{STATE['backend']} stage {STATE['stage']} GN {it}: "
              f"Q={newQ:.12g}, parameter ratio={row['param_ratio']:.4g}, "
              f"RSS ratio={row['rss_ratio']:.4g}, elapsed={time.perf_counter()-STATE['started']:.1f}s", flush=True)
    return param, rss

def observe_outer(old_theta, new_theta, old_sigma, outer, inner_ok):
    new_sigma = STATE['last_sigma']
    param = float(np.max(np.abs(new_theta-old_theta)/(1+np.abs(new_theta))))
    cov = float(np.max(np.abs(new_sigma-old_sigma)/(1+np.abs(new_sigma))))
    passed = param < EPS or cov < COVTOL
    STATE['outer'].append(dict(stage=outer, parameter_change=param, covariance_change=cov,
                              parameter_pass=param < EPS, covariance_pass=cov < COVTOL,
                              inner_success=bool(inner_ok), stop=passed))
    print(f'  nlsur criteria outer {outer}: theta={param:.4g}, Sigma={cov:.4g}, '
          f'inner success={inner_ok}, stop={passed}', flush=True)
    return passed

original_gn = inspect.getsource(solver.gauss_newton)
original_outer = inspect.getsource(solver.nlsur)

class ObserveGN(ast.NodeTransformer):
    def visit_Assign(self, node):
        if (len(node.targets) == 1 and ast.unparse(node.targets[0]) == '(theta, obj)'
                and ast.unparse(node.value) == '(cand, oc)'):
            observer = ast.parse('criterion_param, criterion_rss = observe_accepted(theta, cand, obj, oc, it, mreldif, nrtol, t, m2)').body[0]
            return [observer, node]
        return node

class CriteriaGN(ast.NodeTransformer):
    def visit_Assign(self, node):
        if (len(node.targets) == 1 and ast.unparse(node.targets[0]) == 'converged'
                and ast.unparse(node.value) == 'bool(nrtol < cutoff)'):
            node.value = ast.Constant(False)
        return node
    def visit_If(self, node):
        if ast.unparse(node.test) == "stop_rule == 'standard'":
            return ast.parse('if criterion_param and criterion_rss:\n    converged = True\n    break').body[0]
        return self.generic_visit(node)

class CriteriaOuter(ast.NodeTransformer):
    def visit_If(self, node):
        if ast.unparse(node.test) == 'rel < sigma_tol':
            return ast.parse('if observe_outer(theta_prev, theta, sigma_new, outer, ok):\n    outer_converged = True\n    break').body[0]
        return self.generic_visit(node)

class NormalizeGN(ast.NodeTransformer):
    def visit_Assign(self, node):
        if isinstance(node.value, ast.Call) and isinstance(node.value.func, ast.Name) and node.value.func.id == 'observe_accepted':
            return None
        if len(node.targets) == 1 and ast.unparse(node.targets[0]) == 'converged' and ast.unparse(node.value) in ('bool(nrtol < cutoff)', 'False'):
            node.value = ast.Constant('CONVERGENCE_BOOLEAN')
        return node
    def visit_If(self, node):
        if ast.unparse(node.test) in ("stop_rule == 'standard'", 'criterion_param and criterion_rss'):
            return None
        return self.generic_visit(node)

class NormalizeOuter(ast.NodeTransformer):
    def visit_If(self, node):
        if (ast.unparse(node.test) == 'rel < sigma_tol' or
                isinstance(node.test, ast.Call) and isinstance(node.test.func, ast.Name) and node.test.func.id == 'observe_outer'):
            return None
        return self.generic_visit(node)

standard_gn = ast.fix_missing_locations(ObserveGN().visit(ast.parse(original_gn)))
criteria_gn = ast.fix_missing_locations(CriteriaGN().visit(copy.deepcopy(standard_gn)))
criteria_outer = ast.fix_missing_locations(CriteriaOuter().visit(ast.parse(original_outer)))

def dump_normalized(tree, transformer):
    return ast.dump(transformer().visit(copy.deepcopy(tree)), include_attributes=False)

gn0 = dump_normalized(ast.parse(original_gn), NormalizeGN)
gn1 = dump_normalized(standard_gn, NormalizeGN)
gn2 = dump_normalized(criteria_gn, NormalizeGN)
outer0 = dump_normalized(ast.parse(original_outer), NormalizeOuter)
outer1 = dump_normalized(criteria_outer, NormalizeOuter)
assert gn0 == gn1 == gn2, 'Non-convergence inner optimization arithmetic changed'
assert outer0 == outer1, 'Non-convergence IFGNLS structure changed'
audit = dict(inner_optimization_AST_identical=True, outer_nonconvergence_AST_identical=True,
             standard_has_only_observations=True,
             normalized_inner_sha256=hashlib.sha256(gn0.encode()).hexdigest(),
             normalized_outer_sha256=hashlib.sha256(outer0.encode()).hexdigest(),
             no_accepted_step_policy='Keep original optimizer exit; report False instead of applying old nrtol flag. No new step or acceptance rule.',
             initial_NLS_and_FGNLS_unchanged=True, eps=EPS, tau=TAU, ifgnlseps=COVTOL)
save('unchanged_engine_audit.json', audit)
for name, tree in [('standard_observed_gauss_newton.py', standard_gn),
                   ('criteria_only_gauss_newton.py', criteria_gn),
                   ('criteria_only_nlsur.py', criteria_outer)]:
    (OUT / name).write_text(ast.unparse(tree) + '\n')
for label, old, new in [('inner', ast.parse(original_gn), criteria_gn), ('outer', ast.parse(original_outer), criteria_outer)]:
    (OUT / f'{label}_criteria_diff.patch').write_text(''.join(difflib.unified_diff(
        ast.unparse(old).splitlines(keepends=True), ast.unparse(new).splitlines(keepends=True),
        fromfile='original AST formatted', tofile='criteria-only AST formatted')))

def run(backend):
    with (C11 / 'results/common_inputs.pkl').open('rb') as f:
        inputs = pickle.load(f)
    d, spec = inputs['d'], inputs['spec']
    from pyquaidsce.params import free_slices
    sl = free_slices(spec)
    labels = [None] * spec.n_free
    for name, block in sl.items():
        labels[block] = [f'{name}[{j+1}]' for j in range(block.stop-block.start)]
    STATE.clear()
    STATE.update(backend=backend, labels=labels, stage=0, steps=[], inner=[], outer=[],
                 accepted_thetas=[], stage_thetas=[], stage_sigmas=[])
    namespace = dict(vars(solver))
    namespace.update(observe_accepted=observe_accepted, observe_outer=observe_outer)
    tree = standard_gn if backend == 'standard' else criteria_gn
    exec(compile(tree, '<observed GN; numerics unchanged>', 'exec'), namespace)
    inner_function = namespace['gauss_newton']

    def monitored_inner(*args, **kwargs):
        STATE['stage'] += 1
        old_theta = np.asarray(args[0]).copy()
        old_sigma = np.asarray(args[3]).copy()
        begun = time.perf_counter()
        result = inner_function(*args, **kwargs)
        theta, Q, count, ok = result
        numeric_seconds = time.perf_counter() - begun
        u = solver.residuals(theta, d, spec)
        sigma = u.T @ u / d.nobs
        STATE['last_sigma'] = sigma
        llf = -(d.nobs*spec.n_eq_estimated/2)*(1+solver.LOG2PI)-d.nobs/2*np.linalg.slogdet(sigma)[1]
        legacy_outer_change = float(np.max(np.abs(theta-old_theta)/(np.abs(old_theta)+1e-8)))
        param = float(np.max(np.abs(theta-old_theta)/(1+np.abs(theta))))
        cov = float(np.max(np.abs(sigma-old_sigma)/(1+np.abs(sigma))))
        row = dict(stage=STATE['stage'], phase='NLS' if STATE['stage']==1 else 'FGNLS' if STATE['stage']==2 else 'IFGNLS',
                   iterations=int(count), success=bool(ok), Q=float(Q), seconds=numeric_seconds,
                   llf=float(llf), legacy_outer_change=legacy_outer_change,
                   nlsur_outer_parameter_change=param, nlsur_outer_covariance_change=cov)
        STATE['inner'].append(row)
        STATE['stage_thetas'].append(theta.copy())
        STATE['stage_sigmas'].append(sigma.copy())
        np.savez(OUT / f'{backend}_latest.npz', theta=theta, sigma=sigma,
                 stage=STATE['stage'], n_gn=sum(x['iterations'] for x in STATE['inner']))
        pd.DataFrame(STATE['inner']).to_csv(OUT / f'{backend}_inner_history.csv', index=False)
        pd.DataFrame(STATE['steps']).to_csv(OUT / f'{backend}_accepted_steps.csv', index=False)
        print(f"  {backend} {row['phase']} stage {STATE['stage']}: {count} GN, "
              f"success={ok}, LL={llf:.9f}", flush=True)
        return result

    namespace['gauss_newton'] = monitored_inner
    if backend == 'nlsur_criteria':
        exec(compile(criteria_outer, '<IFGNLS: convergence changes only>', 'exec'), namespace)
        fit = namespace['nlsur']
    else:
        fit = solver.nlsur
    with threadpool_limits(limits=1):
        STATE['started'] = time.perf_counter()
        with patch.object(solver, 'gauss_newton', monitored_inner):
            result = fit(d, spec, theta0=None, sigma0=None, start='zero', method='ifgnls',
                         max_outer=400, max_iter=400, chunk=15000, sigma_tol=1e-5,
                         stop_rule='standard', algorithm='gn', blas_threads=1,
                         log=lambda line: print(f'{backend}: {line}', flush=True))
        seconds = time.perf_counter() - STATE['started']
    with (OUT / f'{backend}_result.pkl').open('wb') as f:
        pickle.dump(result, f, protocol=pickle.HIGHEST_PROTOCOL)
    record = dict(backend=backend, nobs=d.nobs, n_free=spec.n_free,
                  converged=result.converged, n_outer=result.n_outer,
                  n_ifgnls_updates=max(result.n_outer-2, 0), n_gn=result.n_gn,
                  seconds=seconds, llf=result.llf, obj=result.obj,
                  theta=result.theta, sigma=result.sigma, inner_history=STATE['inner'],
                  outer_history=STATE['outer'],
                  earlier_inner_failures=sum(not r['success'] for r in STATE['inner'][:-1]),
                  timing_scope='Demand numerical fit, final covariance assembly and identical observation hooks. Shared preparation and Probits excluded.',
                  n_outer_convention='Includes initial NLS and initial FGNLS, as in package 1.6.0.',
                  no_curvature_constraint=True, start='zero', bootstrap_reps=0)
    save(f'{backend}.json', record)
    pd.DataFrame(STATE['outer']).to_csv(OUT / f'{backend}_outer_history.csv', index=False)
    np.savez(OUT / f'{backend}_parameter_history.npz',
             accepted_thetas=np.asarray(STATE['accepted_thetas']),
             stage_thetas=np.asarray(STATE['stage_thetas']), stage_sigmas=np.asarray(STATE['stage_sigmas']))
    print(json.dumps({k:record[k] for k in ('backend','converged','n_outer','n_gn','seconds','llf')}), flush=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('backend', choices=['standard', 'nlsur_criteria'])
    run(parser.parse_args().backend)
