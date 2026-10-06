"""Compare joint cluster S.E.s with nested whole-cluster bootstrap samples.

Input prices are prepared once and held fixed. No IV or control function is
used. All Probits, demand parameters, residual covariances and sample means
are re-estimated in every draw. Output contains aggregate estimates only.
"""

from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
import multiprocessing as mp
import platform
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'src'))

import numpy as np
import pandas as pd
from pyquaidsce import quaidsce
from pyquaidsce._clusters import cluster_codes, cluster_rows, resample_indices
from pyquaidsce.bootstrap import set_blas_threads

_STATE = {}


def _init(frame, kwargs, groups):
    set_blas_threads(1)
    _STATE.update(frame=frame, kwargs=kwargs, groups=groups)


def _draw(task):
    index, seed = task
    frame, kw, groups = (_STATE[k] for k in ('frame','kwargs','groups'))
    rows = resample_indices(np.random.default_rng(seed), len(frame), groups)
    started = time.perf_counter()
    try:
        result = quaidsce(frame.iloc[rows].reset_index(drop=True), **kw,
                         _deadline=time.perf_counter()+120)
        record = dict(rep=index, seed=seed, nobs=len(rows), converged=result.converged,
                      n_outer=result.n_outer, n_gn=result.n_gn,
                      seconds=time.perf_counter()-started)
        if not result.converged:
            record['error'] = 'second-stage estimator did not converge'
            return record, None
        return record, result.b
    except Exception as exc:
        return dict(rep=index, seed=seed, nobs=len(rows), converged=False,
                    seconds=time.perf_counter()-started,
                    error=f'{type(exc).__name__}: {exc}'), None


def _snapshot(out, method, count, fit, iid, records, estimates, elapsed):
    ordered = sorted(records, key=lambda row: row['rep'])
    good = sorted(estimates)
    vectors = np.vstack([estimates[i] for i in good])
    boot_se = vectors.std(axis=0, ddof=1)
    table = pd.DataFrame(dict(name=fit.names, estimate=fit.b,
                              analytical_iid_se=iid.analytic_se,
                              analytical_cluster_se=fit.analytic_se,
                              bootstrap_se=boot_se))
    table['cluster_to_bootstrap'] = np.divide(fit.analytic_se, boot_se,
        out=np.full_like(boot_se,np.nan), where=boot_se>0)
    table.to_csv(out/f'{method}-{count}-comparison.csv',index=False)
    pd.DataFrame(ordered).to_csv(out/f'{method}-{count}-replications.csv',index=False)
    np.savez_compressed(out/f'{method}-{count}-draws.npz',replications=good,b=vectors)
    n=fit.spec.neqn
    e=table.iloc[-(n+2*n*n):]
    summaries={}
    for kind,block in [('all',e),('income',e.iloc[:n]),('uncompensated',e.iloc[n:n+n*n]),
                        ('compensated',e.iloc[n+n*n:])]:
        ratio=block.cluster_to_bootstrap.to_numpy()
        summaries[kind]=dict(median=float(np.nanmedian(ratio)),min=float(np.nanmin(ratio)),
                             max=float(np.nanmax(ratio)))
    result=dict(method=method,reps_requested=count,reps_ok=len(good),
                failures=[r for r in ordered if not r['converged']],
                bootstrap_seconds=elapsed,elasticity_se_ratios=summaries,
                point_outer=fit.n_outer,point_gn=fit.n_gn,
                inference_seconds=fit.analytical.elapsed_seconds,
                score_diagnostic=fit.analytical.max_standardized_score,
                bread_condition=fit.analytical.bread_condition,
                bootstrap_outer_median=float(np.median([r['n_outer'] for r in ordered if r['converged']])),
                bootstrap_gn_median=float(np.median([r['n_gn'] for r in ordered if r['converged']])))
    (out/f'{method}-{count}-summary.json').write_text(json.dumps(result,indent=2)+'\n')
    # Full failure records are saved above; keep progress output compact.
    print(json.dumps(dict(result,failures=len(result['failures']))),flush=True)
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--shares',nargs='+',required=True)
    parser.add_argument('--prices',nargs='+',required=True)
    parser.add_argument('--expenditure',required=True)
    parser.add_argument('--demographics',nargs='+',required=True)
    parser.add_argument('--selection-covariates',nargs='*',default=None,
                        help='Probit covariates; default uses all demand demographics')
    parser.add_argument('--cluster',required=True)
    parser.add_argument('--anot',type=float,required=True)
    parser.add_argument('--methods',nargs='+',default=['ifgnls','fgnls'],choices=['ifgnls','fgnls','nls'])
    parser.add_argument('--checkpoints',nargs='+',type=int,default=[50,100])
    parser.add_argument('--seed',type=int,default=20261006)
    parser.add_argument('--workers',type=int,default=2)
    args=parser.parse_args()
    checkpoints=sorted(set(args.checkpoints))
    if min(checkpoints)<2 or args.workers<1:
        parser.error('checkpoints must be >=2 and workers >=1')
    args.output.mkdir(parents=True,exist_ok=True)
    frame=pd.read_csv(args.data,dtype={args.cluster:str})
    codes,g=cluster_codes(frame[args.cluster],len(frame))
    groups=cluster_rows(codes,g)
    seeds=np.random.SeedSequence(args.seed).generate_state(max(checkpoints),dtype=np.uint32)
    settings=dict(shares=args.shares,prices=args.prices,expenditure=args.expenditure,
                  demographics=args.demographics,anot=args.anot,censor=True,quadratic=True,
                  start='zero',algorithm='gn',param_tol=1e-5,objective_tol=1e-7,
                  gn_tol=1e-5,outer_param_tol=1e-5,max_outer=400,max_iter=400,
                  chunk=2000,blas_threads=1,verbose=False)
    if args.selection_covariates is not None:
        settings['selection_covariates']=args.selection_covariates
    metadata=dict(nobs=len(frame),n_clusters=g,cluster=args.cluster,seed=args.seed,checkpoints=checkpoints,
                  input_sha256=hashlib.sha256(args.data.read_bytes()).hexdigest(),
                  prices_fixed=True,ivexp=None,control_function=None,cluster_correction='G/(G-1)',
                  settings=settings,workers=args.workers,python=platform.python_version(),
                  numpy=np.__version__,pandas=pd.__version__,
                  source_sha256={name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest()
                                 for name in ['src/pyquaidsce/inference.py','src/pyquaidsce/_clusters.py',
                                              'tools/validate_cluster_se.py']})
    (args.output/'settings.json').write_text(json.dumps(metadata,indent=2)+'\n')
    for method in args.methods:
        kw=dict(settings,method=method)
        started=time.perf_counter()
        iid=quaidsce(frame,**kw,analytic=True)
        fit=quaidsce(frame,**kw,analytic=True,cluster=args.cluster)
        for name in ('theta','b','tau','sigma'):
            np.testing.assert_array_equal(getattr(iid,name),getattr(fit,name))
        for name in ('llf','n_outer','n_gn','converged'):
            assert getattr(iid,name)==getattr(fit,name)
        print(json.dumps(dict(method=method,point_seconds=time.perf_counter()-started,
                              outer=fit.n_outer,inner=fit.n_gn)),flush=True)
        records,estimates=[],{}
        bootstrap_started=time.perf_counter()
        with ProcessPoolExecutor(max_workers=args.workers,mp_context=mp.get_context('spawn'),
                initializer=_init,initargs=(frame,kw,groups)) as pool:
            previous=0
            for count in checkpoints:
                tasks=[(i+1,int(seeds[i])) for i in range(previous,count)]
                futures=[pool.submit(_draw,task) for task in tasks]
                for future in as_completed(futures):
                    record,b=future.result()
                    records.append(record)
                    if b is not None:
                        estimates[record['rep']]=b
                    if len(records)%10==0:
                        print(f'{method}: {len(records)}/{count}, successful={len(estimates)}',flush=True)
                if len(estimates)<2:
                    raise RuntimeError('fewer than two bootstrap fits converged')
                _snapshot(args.output,method,count,fit,iid,records,estimates,
                          time.perf_counter()-bootstrap_started)
                previous=count


if __name__=='__main__':
    main()
