"""Archived influence and leave-one-cluster diagnostic for a completed comparison."""

import argparse
from pathlib import Path
import sys,json,time
from unittest.mock import patch
REPO=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(REPO/'src'))
import numpy as np,pandas as pd
from pyquaidsce import quaidsce
from pyquaidsce._clusters import cluster_codes
from pyquaidsce.inference import compute_analytical_inference,_EstimatingSystem,_bread_inverse,_pack_means
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--data',type=Path,required=True)
parser.add_argument('--comparison',type=Path,required=True)
args=parser.parse_args()
ROOT=args.comparison
settings=json.loads((ROOT/'settings.json').read_text());kw=settings['settings']
cluster=settings['cluster']
frame=pd.read_csv(args.data,dtype={cluster:str});codes,g=cluster_codes(frame[cluster],len(frame))
ctx={}
def capture(result,data,design,**kw):
    ctx.update(data=data,design=design)
    return compute_analytical_inference(result,data,design,**kw)
with patch('pyquaidsce.inference.compute_analytical_inference',side_effect=capture):
    base=quaidsce(frame,**kw,method='ifgnls',analytic=True,cluster=cluster)
s=_EstimatingSystem(base,ctx['data'],ctx['design'],None,2000,None)
inv,_=_bread_inverse(s.bread(),s.scales);direct,mj=s.elasticity_jacobians();total=direct+mj@s.mean_derivative()
sums=np.zeros((g,36));scoremean=s.mean_score(s.x);means=_pack_means(base.means)
for sl,scores,features in s.chunks(s.x,include_means=True):
    values=-(scores-scoremean)@inv.T@total.T+(features-means)@mj.T
    np.add.at(sums,codes[sl],values)
squared=sums*sums
ranked=np.sort(squared,axis=0)[::-1]
fractions=np.cumsum(ranked,axis=0)/squared.sum(axis=0)
concentration=dict(income2_top1=float(fractions[0,1]),
                   income2_top5=float(fractions[4,1]),
                   income2_top10=float(fractions[9,1]),
                   median_elasticity_top10=float(np.median(fractions[9])),
                   maximum_elasticity_top10=float(np.max(fractions[9])))
(ROOT/'influence-concentration.json').write_text(json.dumps(concentration,indent=2)+'\n')
sizes=np.bincount(codes)
prediction=-sums/(len(frame)-sizes)[:,None]
high=np.argsort(abs(prediction[:,1]))[-8:]
random=np.random.default_rng(721).choice(g,size=12,replace=False)
selected=np.unique(np.r_[high,random])
point=base.elas.as_stata_vector();rows=[]
for code in selected:
    f=frame.loc[codes!=code].reset_index(drop=True);started=time.perf_counter()
    options=dict(kw,param_tol=1e-8,objective_tol=1e-10,gn_tol=1e-8,outer_param_tol=1e-8)
    fit=quaidsce(f,**options,method='ifgnls',initial=base.theta,sigma_initial=base.sigma,
                _deadline=time.perf_counter()+60)
    actual=fit.elas.as_stata_vector()-point
    row=dict(cluster_code=int(code),size=int(sizes[code]),converged=fit.converged,
             seconds=time.perf_counter()-started,outer=fit.n_outer,inner=fit.n_gn,
             predicted=prediction[code].tolist(),actual=actual.tolist())
    rows.append(row)
    print(json.dumps(dict(code=int(code),predicted_income2=float(prediction[code,1]),
                         actual_income2=float(actual[1]),converged=fit.converged)),flush=True)
actual=np.array([r['actual'] for r in rows]);pred=np.array([r['predicted'] for r in rows])
report=dict(records=rows,selection='eight largest predicted income-2 changes plus twelve seed-721 random clusters',
            income2_relative_rms_error=float(np.linalg.norm(actual[:,1]-pred[:,1])/np.linalg.norm(pred[:,1])),
            total_relative_rms_error=float(np.linalg.norm(actual-pred)/np.linalg.norm(pred)),
            income2_correlation=float(np.corrcoef(actual[:,1],pred[:,1])[0,1]))
(ROOT/'delete-cluster-probe.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:v for k,v in report.items() if k!='records'}),flush=True)
