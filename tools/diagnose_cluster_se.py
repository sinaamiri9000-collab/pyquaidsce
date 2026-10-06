"""Independent real-data derivative, separation and bootstrap-worker checks.

Use settings saved by tools/validate_cluster_se.py and its prepared input CSV.
This diagnostic does not change the optimizer or any estimation thresholds.
"""

import argparse,sys,time,json,re
from pathlib import Path
from unittest.mock import patch
REPO=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(REPO/'src'),str(REPO)]
import numpy as np,pandas as pd
from scipy.optimize import linprog
from pyquaidsce import quaidsce
from pyquaidsce.inference import _EstimatingSystem,_central_jacobian,_bread_inverse,_pack_means,compute_analytical_inference
from pyquaidsce.model import jacobian_full,fitted_shares
from pyquaidsce.params import delta_matrix
from pyquaidsce.probit import probit,_lambda_ratio
from pyquaidsce._clusters import cluster_codes,cluster_rows,resample_indices


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data',type=Path,required=True)
    parser.add_argument('--comparison',type=Path,required=True)
    args=parser.parse_args()
    ROOT=args.comparison
    settings=json.loads((ROOT/'settings.json').read_text())
    cluster=settings['cluster']
    frame=pd.read_csv(args.data,dtype={cluster:str})
    kw=settings['settings']; codes,g=cluster_codes(frame[cluster],len(frame))
    groups=cluster_rows(codes,g)
    report={}
    for method in ['ifgnls','fgnls']:
        captured={}
        def capture(result,data,design,**kwargs):
            captured.update(data=data,design=design,kwargs=kwargs)
            return compute_analytical_inference(result,data,design,**kwargs)
        started=time.perf_counter()
        with patch('pyquaidsce.inference.compute_analytical_inference',side_effect=capture):
            fit=quaidsce(frame,**kw,method=method,analytic=True,cluster=cluster)
        elapsed=time.perf_counter()-started
        system=_EstimatingSystem(fit,captured['data'],captured['design'],captured['kwargs']['theta_nls'],2000,None)
        bread=system.bread()
        h=np.cbrt(np.finfo(float).eps)/4
        a=_central_jacobian(system.mean_score,system.x,system.scales,h)
        b=_central_jacobian(system.mean_score,system.x,system.scales,h*2)
        ref=(4*a-b)/3
        tau,theta,sigma=system.decode(system.x)
        data=captured['data']
        jac=jacobian_full(theta,data,fit.spec) @ delta_matrix(fit.spec)
        u=data.shares[:,:fit.spec.n_eq_estimated]-fitted_shares(theta,data,fit.spec)
        for column,(aa,bb) in enumerate(zip(*system.tril)):
            step=1e-20*system.scales[system.sigma_slice.start+column]
            perturb=sigma.astype(complex)
            perturb[aa,bb]+=1j*step
            if aa!=bb: perturb[bb,aa]+=1j*step
            grad=np.einsum('tmi,tm->i',jac,u@np.linalg.inv(perturb))
            ref[system.theta_slice,system.sigma_slice.start+column]=grad.imag/step/data.nobs
        direct,means=system.elasticity_jacobians()
        nd,nm=system.elasticity_jacobians(relative_step=1e-6)
        numerical=compute_analytical_inference(fit,data,captured['design'],**captured['kwargs'],_relative_step=1e-6)
        information=-bread[system.theta_slice,system.theta_slice]
        sd=np.sqrt(np.abs(np.diag(information)))
        scaled_information=information/np.outer(sd,sd)
        coefvar=fit.analytical.joint_covariance[:fit.spec.n_free,:fit.spec.n_free]
        stats=dict(point_and_inference_seconds=elapsed,inference_seconds=fit.analytical.elapsed_seconds,
                   point_seconds=elapsed-fit.analytical.elapsed_seconds,outer=fit.n_outer,inner=fit.n_gn,
                   scaled_bread_error=float(np.max(np.abs(bread-ref)/np.maximum(1,np.abs(bread)))),
                   scaled_elasticity_gradient_error=float(np.max(np.abs(np.column_stack([direct,means])-np.column_stack([nd,nm]))/(1+np.abs(np.column_stack([nd,nm]))))),
                   max_relative_elasticity_se_error=float(np.max(np.abs(numerical.elasticity_se/fit.analytical.elasticity_se-1))),
                   min_final_demand_information_eigenvalue=float(np.linalg.eigvalsh(information).min()),
                   min_scaled_information_eigenvalue=float(np.linalg.eigvalsh(scaled_information).min()),
                   sigma_eigenvalues=np.linalg.eigvalsh(sigma).tolist(),adjusted_shares=fit.elas.we.tolist(),
                   rho=fit.coefs.rho.tolist(),bread_condition=fit.analytical.bread_condition,
                   score_diagnostic=fit.analytical.max_standardized_score)
        # Compare S.E.s using the independently differentiated full bread as well,
        # since a small entrywise error can matter in an ill-conditioned system.
        inverse,_=_bread_inverse(ref,system.scales)
        total=direct+means@system.mean_derivative()
        sums=np.zeros((g,len(fit.analytical.elasticity_se)))
        score_mean=system.mean_score(system.x)
        mean_values=_pack_means(fit.means)
        for sl,scores,features in system.chunks(system.x,include_means=True):
            values=-(scores-score_mean)@inverse.T@total.T+(features-mean_values)@means.T
            np.add.at(sums,codes[sl],values)
        se_reference=np.sqrt(np.sum(sums*sums,axis=0)*g/(g-1))/len(frame)
        stats['max_relative_se_error_numerical_bread']=float(np.max(np.abs(
            se_reference/fit.analytical.elasticity_se-1)))
        report[method]=stats
        if method=='ifgnls':
            point_ifgnls=fit
        print(method,json.dumps(stats),flush=True)

    # Examine all remaining Probit failures with the unchanged backend.
    replications=pd.read_csv(ROOT/'ifgnls-100-replications.csv')
    failed=replications.loc[~replications.converged]
    checks=[]
    for _,r in failed.iterrows():
        match=re.search(r'first-stage probit (\d+)', str(r.error))
        if match is None:
            continue
        equation=int(match.group(1))-1
        rows=resample_indices(np.random.default_rng(int(r.seed)),len(frame),groups)
        f=frame.iloc[rows]
        covariates=kw.get('selection_covariates',kw['demographics'])
        X=np.column_stack([np.log(f[kw['prices']]),np.log(f[kw['expenditure']]),f[covariates]])
        y=(f[kw['shares'][equation]].to_numpy()>0).astype(float)
        pr=probit(y,X)
        Z=np.column_stack([X,np.ones(len(X))])
        q=2*y-1; ratio=_lambda_ratio(q*pr.xb(X))
        scores=Z*(q*ratio)[:,None]
        gradient=scores.sum(axis=0)
        standardized_gradient=np.max(np.abs(gradient)/np.sqrt((scores**2).sum(axis=0)))
        # A separating direction has nonnegative signed margins and at least one
        # strictly positive margin. Bounded LP coordinates avoid scale ambiguity.
        scaled=Z/np.maximum(np.linalg.norm(Z,axis=0),1)
        signed=q[:,None]*scaled
        lp=linprog(-signed.mean(axis=0),A_ub=-signed,b_ub=np.zeros(len(Z)),
                   bounds=[(-1,1)]*Z.shape[1],method='highs')
        checks.append(dict(rep=int(r.rep),equation=equation+1,zeros=int((y==0).sum()),iterations=pr.n_iter,
                           converged=pr.converged,max_standardized_score=float(standardized_gradient),
                           lp_success=bool(lp.success),separating_margin=float(-lp.fun) if lp.success else None,
                           min_margin=float((signed@lp.x).min()) if lp.success else None))
    report['failed_probit_checks']=checks

    # The diagnostic runner uses the production sampler/seed convention.
    seeds=np.random.SeedSequence(settings['seed']).generate_state(100,dtype=np.uint32)
    from pyquaidsce.bootstrap import _one_rep,_WORK
    worker_kw=dict(kw,method='ifgnls')
    worker_kw.pop('verbose',None)
    if settings.get('bootstrap_start','zero')=='warm':
        worker_kw.update(initial=point_ifgnls.theta,sigma_initial=point_ifgnls.sigma)
    _WORK.update(df=frame,kw=worker_kw,rep_timeout=None,cluster_rows=groups)
    saved=np.load(ROOT/'ifgnls-100-draws.npz')
    confirmed=[]
    for rep in [int(x) for x in saved['replications'][:3]]:
        index,vector,error=_one_rep((rep,int(seeds[rep-1])))
        assert error is None
        expected=saved['b'][np.flatnonzero(saved['replications']==rep)[0]]
        np.testing.assert_array_equal(vector,expected)
        confirmed.append(rep)
    report['native_worker_draws_bitwise_equal']=confirmed
    (ROOT/'diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(checks),flush=True)


if __name__ == "__main__":
    main()
