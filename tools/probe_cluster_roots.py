"""Archived local-root diagnostic for a completed cluster-bootstrap comparison."""

import argparse
from pathlib import Path
import sys,time,json
from unittest.mock import patch
REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'src'))
import numpy as np,pandas as pd
from pyquaidsce import quaidsce
from pyquaidsce._clusters import cluster_codes,cluster_rows,resample_indices
from pyquaidsce.inference import compute_analytical_inference,_EstimatingSystem
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--comparison',type=Path,required=True)
parser.add_argument('--draw',type=int,default=3)
args=parser.parse_args()
if args.draw<1:
    parser.error('draw must be positive')
ROOT=args.comparison
settings=json.loads((ROOT/'settings.json').read_text())
cluster=settings['cluster']
frame=pd.read_csv(args.data,dtype={cluster:str})
kw=settings['settings']
base=quaidsce(frame,**kw,method='ifgnls')
codes,g=cluster_codes(frame[cluster],len(frame));groups=cluster_rows(codes,g)
seeds=np.random.SeedSequence(settings['seed']).generate_state(max(100,args.draw),dtype=np.uint32)
rows=resample_indices(np.random.default_rng(int(seeds[args.draw-1])),len(frame),groups)
rep=frame.iloc[rows].reset_index(drop=True)
cases=[('zero-default',{},None,None),('warm-default',{},base.theta,base.sigma),
       ('zero-stricter',dict(param_tol=1e-8,objective_tol=1e-10,gn_tol=1e-8,outer_param_tol=1e-8),None,None),
       ('warm-stricter',dict(param_tol=1e-8,objective_tol=1e-10,gn_tol=1e-8,outer_param_tol=1e-8),base.theta,base.sigma)]
records=[]
for name,tolerances,initial,sigma_initial in cases:
    options=dict(kw,**tolerances);ctx={}
    def capture(result,data,design,**options):
        ctx.update(data=data,design=design)
        return compute_analytical_inference(result,data,design,**options)
    start=time.perf_counter()
    try:
        with patch('pyquaidsce.inference.compute_analytical_inference',side_effect=capture):
            fit=quaidsce(rep,**options,method='ifgnls',initial=initial,sigma_initial=sigma_initial,
                         analytic=True,_deadline=time.perf_counter()+60)
        system=_EstimatingSystem(fit,ctx['data'],ctx['design'],None,2000,None)
        bread=system.bread();info=-bread[system.theta_slice,system.theta_slice]
        profile=-(bread[system.theta_slice,system.theta_slice]-bread[system.theta_slice,system.sigma_slice]@np.linalg.solve(bread[system.sigma_slice,system.sigma_slice],bread[system.sigma_slice,system.theta_slice]))
        profile=(profile+profile.T)/2
        rec=dict(case=name,seconds=time.perf_counter()-start,converged=fit.converged,outer=fit.n_outer,
                 inner=fit.n_gn,llf=fit.llf,rho=fit.coefs.rho.tolist(),income=fit.elas.income.tolist(),
                 income_se=fit.analytical.income_se.tolist(),condition=fit.analytical.bread_condition,
                 standardized_score=fit.analytical.max_standardized_score,
                 min_fixed_information_eigenvalue=float(np.linalg.eigvalsh(info).min()),
                 min_profile_information_eigenvalue=float(np.linalg.eigvalsh(profile).min()))
    except Exception as exc:
        rec=dict(case=name,seconds=time.perf_counter()-start,error=f'{type(exc).__name__}: {exc}')
    records.append(rec);print(json.dumps(rec),flush=True)
(ROOT/'root-probe.json').write_text(json.dumps(records,indent=2)+'\n')
