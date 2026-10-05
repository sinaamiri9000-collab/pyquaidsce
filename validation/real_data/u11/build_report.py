"""Compare saved U11 fits and create a private aggregate report, without refits."""
from pathlib import Path
import sys, json, pickle, hashlib, zipfile, subprocess
ROOT = Path(__file__).resolve().parent
REPO = Path('/workspace/pyquaidsce')
C11 = Path('/workspace/real-data-c11')
sys.path[:0] = [str(REPO/'src'), str(C11/'D1_C11_FINAL_PSU_BOOTSTRAP_V2')]
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
from pyquaidsce.params import unpack, full_vector, free_slices
from pyquaidsce.model import fitted_shares, residuals
from pyquaidsce.elasticities import elasticities
OUT = ROOT/'results'

def save(name, obj):
    def fallback(x):
        if isinstance(x, np.ndarray): return x.tolist()
        if isinstance(x, np.generic): return x.item()
        raise TypeError(type(x).__name__)
    (OUT/name).write_text(json.dumps(obj, default=fallback, indent=2, allow_nan=False)+'\n')

with (C11/'results/common_inputs.pkl').open('rb') as f: inputs=pickle.load(f)
with (C11/'results/auxiliary_u11.pkl').open('rb') as f: aux=pickle.load(f)
fits={}
for name in ('standard','nlsur_criteria'):
    with (OUT/f'{name}_result.pkl').open('rb') as f: fits[name]=pickle.load(f)
records={name:json.loads((OUT/f'{name}.json').read_text()) for name in fits}
d,spec=inputs['d'],inputs['spec']
names=spec.full_names(aux.demo_names)
vectors={name:full_vector(r.theta,spec) for name,r in fits.items()}
diffcoef=vectors['nlsur_criteria']-vectors['standard']
coefficients=pd.DataFrame(dict(coefficient=names,standard=vectors['standard'],
    nlsur_criteria=vectors['nlsur_criteria'],difference=diffcoef,absolute_difference=np.abs(diffcoef)))
coefficients.to_csv(OUT/'coefficient_comparison.csv',index=False)
blocks=[]
for block,group in coefficients.groupby(coefficients.coefficient.str.split(':').str[0],sort=False):
    blocks.append(dict(block=block,max_abs_difference=float(group.absolute_difference.max()),
                       rms_difference=float(np.sqrt(np.mean(group.difference**2)))))
pd.DataFrame(blocks).to_csv(OUT/'coefficient_blocks.csv',index=False)
keys=('backend','nobs','n_free','converged','n_outer','n_ifgnls_updates','n_gn','seconds','llf','earlier_inner_failures')
summary=pd.DataFrame([{k:r[k] for k in keys} for r in records.values()])
summary.to_csv(OUT/'comparison.csv',index=False)
phase_rows=[]
for name,r in records.items():
    history=pd.DataFrame(r['inner_history'])
    for phase,group in history.groupby('phase',sort=False):
        steps=int(group.iterations.sum())
        seconds=float(group.seconds.sum())
        phase_rows.append(dict(backend=name,phase=phase,inner_calls=len(group),
            GN_steps=steps,inner_seconds=seconds,seconds_per_GN=seconds/steps))
pd.DataFrame(phase_rows).to_csv(OUT/'phase_summary.csv',index=False)
end_rows=[]
for name in records:
    trace=pd.read_csv(OUT/f'{name}_accepted_steps.csv',float_precision='round_trip')
    for stage,group in trace.groupby('stage',sort=False):
        last=group.iloc[-1]
        end_rows.append(dict(backend=name,stage=int(stage),last_accepted_iteration=int(last.iteration),
            param_ratio=float(last.param_ratio),rss_ratio=float(last.rss_ratio),
            worst_theta_component=last.worst_theta_component,nrtol=float(last.nrtol),
            standard_mreldif=float(last.standard_mreldif),
            standard_gradient_pass=last.nrtol<1e-5,
            standard_parameter_pass=last.standard_mreldif<1e-5,
            standard_objective_pass=(last.Q_old-last.Q)/abs(last.Q_old)<1e-7))
pd.DataFrame(end_rows).to_csv(OUT/'stage_end_diagnostics.csv',index=False)

checks=[]
def check(name,passed,value=None):
    checks.append(dict(check=name,passed=bool(passed),value=value))
with threadpool_limits(limits=1):
    predictions={name:fitted_shares(r.theta,d,spec) for name,r in fits.items()}
    els={name:elasticities(unpack(r.theta,spec),spec,aux.means,a0=d.a0,
            tau=aux.tau,np_prob=aux.np_prob,layout=aux.selection_layout) for name,r in fits.items()}
    e0,e1=els['standard'],els['nlsur_criteria']
    delta_share=predictions['nlsur_criteria']-predictions['standard']
    v0,v1=fits['standard'].V,fits['nlsur_criteria'].V
    difference=dict(max_abs_free_coefficient_difference=float(np.max(np.abs(fits['nlsur_criteria'].theta-fits['standard'].theta))),
        max_abs_full_coefficient_difference=float(np.max(np.abs(diffcoef))),
        rms_full_coefficient_difference=float(np.sqrt(np.mean(diffcoef**2))),
        fitted_share_max_abs_difference=float(np.max(np.abs(delta_share))),
        fitted_share_rms_difference=float(np.sqrt(np.mean(delta_share**2))),
        likelihood_difference=fits['nlsur_criteria'].llf-fits['standard'].llf,
        residual_covariance_relative_frobenius_difference=float(np.linalg.norm(fits['nlsur_criteria'].sigma-fits['standard'].sigma)/np.linalg.norm(fits['standard'].sigma)),
        parameter_covariance_relative_frobenius_difference=float(np.linalg.norm(v1-v0)/np.linalg.norm(v0)),
        all_elasticities_max_abs_difference=float(np.max(np.abs(e1.as_stata_vector()-e0.as_stata_vector()))),
        marshallian_max_abs_difference=float(np.max(np.abs(e1.uncompensated-e0.uncompensated))),
        hicksian_max_abs_difference=float(np.max(np.abs(e1.compensated-e0.compensated))),
        income_max_abs_difference=float(np.max(np.abs(e1.income-e0.income))),
        nlsur_criteria_to_standard_time_ratio=records['nlsur_criteria']['seconds']/records['standard']['seconds'],
        standard_vs_cached_auxiliary_max_abs_theta_difference=float(np.max(np.abs(fits['standard'].theta-aux.theta))),
        standard_vs_cached_auxiliary_likelihood_difference=fits['standard'].llf-aux.llf,
        largest_coefficient_changes=coefficients.nlargest(12,'absolute_difference').to_dict('records'))
    for name,r in fits.items():
        co=unpack(r.theta,spec)
        u=residuals(r.theta,d,spec)
        sigma=u.T@u/d.nobs
        steps=pd.read_csv(OUT/f'{name}_accepted_steps.csv',float_precision='round_trip')
        h=pd.read_csv(OUT/f'{name}_inner_history.csv')
        check(f'{name}: finite theta',np.isfinite(r.theta).all())
        check(f'{name}: finite parameter covariance',np.isfinite(r.V).all())
        check(f'{name}: positive definite residual covariance',np.linalg.eigvalsh(r.sigma).min()>0,float(np.linalg.eigvalsh(r.sigma).min()))
        check(f'{name}: covariance equals residual covariance',np.max(np.abs(sigma-r.sigma))<1e-12,float(np.max(np.abs(sigma-r.sigma))))
        check(f'{name}: Gamma symmetry',np.max(np.abs(co.gamma-co.gamma.T))<1e-12)
        check(f'{name}: Gamma homogeneity',np.max(np.abs(co.gamma.sum(axis=1)))<1e-10)
        check(f'{name}: eta adding-up',np.max(np.abs(co.eta.sum(axis=1)))<1e-10)
        check(f'{name}: positive demographic scaling',np.min(1+d.demo@co.rho)>0,float(np.min(1+d.demo@co.rho)))
        check(f'{name}: GN trace counts',int(h.iterations.sum())==r.n_gn)
        check(f'{name}: outer trace counts',len(h)==r.n_outer)
        check(f'{name}: accepted objective decreases',(steps.Q<steps.Q_old).all())
        check(f'{name}: finite elasticities',np.isfinite(els[name].as_stata_vector()).all())
        for kind,arr in [('marshallian',els[name].uncompensated),('hicksian',els[name].compensated),('income',els[name].income)]:
            pd.DataFrame(arr).to_csv(OUT/f'{name}_{kind}_elasticities.csv',index=False)
    history0=np.load(OUT/'standard_parameter_history.npz')['accepted_thetas']
    history1=np.load(OUT/'nlsur_criteria_parameter_history.npz')['accepted_thetas']
    a,b=history0[0],history1[0]
    first_difference=float(np.max(np.abs(a-b)))
    check('exactly equal first accepted theta',np.array_equal(a,b),first_difference)
    prefix=min(records['standard']['inner_history'][0]['iterations'],records['nlsur_criteria']['inner_history'][0]['iterations'])
    check('exactly equal common initial NLS trajectory',np.array_equal(history0[:prefix],history1[:prefix]),int(prefix))
    check('standard reproduces cached unrestricted U11',difference['standard_vs_cached_auxiliary_max_abs_theta_difference']<1e-8,difference['standard_vs_cached_auxiliary_max_abs_theta_difference'])
audit=json.loads((OUT/'unchanged_engine_audit.json').read_text())
check('unchanged inner optimization AST',audit['inner_optimization_AST_identical'])
check('unchanged outer nonconvergence AST',audit['outer_nonconvergence_AST_identical'])
dirty=subprocess.check_output(['git','diff','--name-only','HEAD','--','src/pyquaidsce'],cwd=REPO,text=True)
check('production source files unmodified',not dirty.strip())
save('checks.json',checks)
save('differences.json',difference)
assert all(c['passed'] for c in checks),[c for c in checks if not c['passed']]

arrays={key:dict(shape=getattr(d,key).shape,dtype=str(getattr(d,key).dtype),
         sha256=hashlib.sha256(np.ascontiguousarray(getattr(d,key)).tobytes()).hexdigest())
         for key in ('lnp','lnexp','shares','demo','cdf','pdf','control_function')}
manifest=dict(package='pyquaidsce 1.6.0',repository_branch=subprocess.check_output(['git','branch','--show-current'],cwd=REPO,text=True).strip(),
    repository_HEAD=subprocess.check_output(['git','rev-parse','HEAD'],cwd=REPO,text=True).strip(),
    data_sha256=json.loads((C11/'results/preparation.json').read_text())['data_sha256'],
    demand_inputs=arrays,nobs=d.nobs,n_free=spec.n_free,n_eq_estimated=spec.n_eq_estimated,
    shared_selection=True,first_stage_estimates='Previously fitted full-sample U11 Probits, unchanged for both runs.',
    initial_theta='All zeros, no theta0/sigma0; original initial NLS and FGNLS stages preserved.',
    no_curvature_constraint=True,no_linear_start=True,bootstrap_reps=0,fits_per_backend=1,
    max_inner=400,max_outer=400,chunk=15000,BLAS_threads=1,eps=1e-5,tau=1e-3,ifgnlseps=1e-10,
    native_stata_executed=False,optimization_engine_changed=False,checks_passed=len(checks),
    timing_scope='New single numerical fit for each backend, including final covariance assembly and common recording hooks. Shared preparation and first-stage Probits excluded.',
    scripts={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.glob('*.py')},
    production_nlsur_sha256=hashlib.sha256((REPO/'src/pyquaidsce/nlsur.py').read_bytes()).hexdigest(),
    publication='Private local results. No raw household arrays or pickles in the shareable bundle.')
save('experiment_manifest.json',manifest)

def table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(map(str,row))+' |' for row in rows])
main=table(['روش','همگرا','کل مراحل بیرونی','به‌روزرسانی IFGNLS','گام داخلی','ثانیه','LL'],
    [[r.backend,str(r.converged),r.n_outer,r.n_ifgnls_updates,r.n_gn,f'{r.seconds:.3f}',f'{r.llf:.9f}'] for r in summary.itertuples()])
phase_table=table(['روش','مرحله','تعداد حل داخلی','کل گام GN','ثانیهٔ حل‌های داخلی','ثانیه به‌ازای GN'],
    [[r['backend'],r['phase'],r['inner_calls'],r['GN_steps'],f"{r['inner_seconds']:.3f}",f"{r['seconds_per_GN']:.3f}"] for r in phase_rows])
top=difference['largest_coefficient_changes'][0]
outcome=('هر دو اجرای U11 همگرا شده‌اند.' if all(r.converged for r in fits.values())
         else 'وضعیت همگرایی دو اجرا در جدول گزارش شده است.')
report=f'''# مقایسهٔ U11 بدون قید انحنا: standard و فقط معیارهای nlsur

{outcome} در این داده، نتیجهٔ standard از نظر ضرایب، سهم‌ها، درست‌نمایی و کشش‌ها بسیار نزدیک به نتیجهٔ معیارهای nlsur است. بیشترین اختلاف ضریب حدود {difference['max_abs_full_coefficient_difference']:.3g} و قدرمطلق اختلاف LL حدود {abs(difference['likelihood_difference']):.3g} است. افزایش زمان حدود {difference['nlsur_criteria_to_standard_time_ratio']:.2f} برابر، همراه با افزایش گام‌های داخلی رخ داده است. این نتیجه شاهدی دربارهٔ همین U11 است و به‌تنهایی جای بررسی همهٔ مدل‌ها و گزینه‌های پکیج را نمی‌گیرد.

{main}

هر دو اجرای جدید روی همان دادهٔ اصلی آماده‌شده، ۳۷٬۶۵۸ مشاهده، ۱۱ کالا، ۸ متغیر جمعیتی و ۱۹۸ پارامتر آزاد انجام شده‌اند. `anot=11`، مدل درجه‌دو، تصحیح سانسور و تابع کنترل همان مشخصات U11 در بستهٔ ارسالی‌اند. U11 قید انحنای C11 را ندارد و در مختصات معمول پارامترهای پکیج تخمین زده می‌شود.

شروع هر دو صفر است؛ نقطهٔ تخمین قبلی به‌عنوان شروع استفاده نشده است. NLS اولیه با ماتریس واحد و سپس FGNLS اولیه و IFGNLS خود نسخهٔ ۱٫۶٫۰ حفظ شده‌اند. عدد n_outer پکیج شامل دو مرحلهٔ اولیه است؛ ستون به‌روزرسانی IFGNLS آن دو را کنار می‌گذارد. بوت‌استرپ و شروع خطی انجام نشده‌اند.

## زمان و اختلاف نتایج

زمان جدول فقط خود تخمین تقاضا، ساخت کوواریانس نهایی و ثبت تاریخچهٔ یکسان است. آماده‌سازی داده و Probitها مشترک و خارج از این زمان‌اند. هر روش فقط یک بار تخمین زده شده؛ نسبت زمانی این اجرا benchmark چندباره نیست.

- نسبت زمان معیارهای nlsur به standard: {difference['nlsur_criteria_to_standard_time_ratio']:.6g}.
- بیشترین اختلاف ضریب آزاد: {difference['max_abs_free_coefficient_difference']:.12g}.
- بیشترین اختلاف ضریب کامل: {difference['max_abs_full_coefficient_difference']:.12g}؛ برای `{top['coefficient']}`، از {top['standard']:.12g} به {top['nlsur_criteria']:.12g}.
- RMS اختلاف ضرایب کامل: {difference['rms_full_coefficient_difference']:.12g}.
- بیشترین اختلاف سهم پیش‌بینی‌شده: {difference['fitted_share_max_abs_difference']:.12g}؛ RMS: {difference['fitted_share_rms_difference']:.12g}.
- تغییر LL: {difference['likelihood_difference']:+.12g}.
- بیشترین اختلاف کشش‌ها در میانگین مشترک: {difference['all_elasticities_max_abs_difference']:.12g}؛ درآمدی: {difference['income_max_abs_difference']:.12g}.
- اختلاف نسبی فروبنیوس کوواریانس باقیمانده: {difference['residual_covariance_relative_frobenius_difference']:.12g}؛ کوواریانس پارامترها: {difference['parameter_covariance_relative_frobenius_difference']:.12g}.

تفکیک هزینه و گام‌ها به مرحلهٔ اولیه و IFGNLS:

{phase_table}

ثانیهٔ حل‌های داخلی جمع زمان ثبت‌شدهٔ فراخوانی‌های GN است؛ ثبت snapshot پس از حل و ساخت کوواریانس نهایی در این ستون نیستند، ولی در زمان کل جدول اصلی هستند. هزینهٔ یک گام می‌تواند با تعداد نصف‌کردن گام تغییر کند؛ مقایسهٔ زمان کل باید همراه شمار گام‌ها و نتیجهٔ همگرایی خوانده شود.

کشش‌های مارشالی، هیکسی و درآمدی با تابع عمومی نسخهٔ ۱٫۶٫۰ و میانگین‌ها، tau و مشخصات انتخاب یکسان محاسبه شده‌اند؛ این مقایسهٔ کشش در نقطهٔ میانگین است، نه در تک‌تک خانوارها. CSV همهٔ ضرایب و کشش‌ها همراه گزارش است. پرچم همگرایی هر روش فقط به قرارداد توقف خودش مربوط است و بهینهٔ سراسری را اثبات نمی‌کند.

## شروط مورد مقایسه

standard همان شروط نسخهٔ ۱٫۶٫۰ را نگه می‌دارد: توقف داخلی با OR تغییر پارامتر، تغییر هدف و معیار GN؛ توقف بیرونی با تغییر نسبی پارامتر کمتر از sigma_tol در دو دور متوالی.

در روش دوم، تنها شروط اعلام همگرایی از nlsur.ado ارسالی گرفته شده‌اند:

```
برای تمام پارامترهای آزاد:
abs(theta_new - theta_old) <= 1e-5 * (abs(theta_old) + 1e-3)
و هم‌زمان:
abs(Q_new - Q_old) <= 1e-5 * (Q_old + 1e-3)

بیرونی:
max(abs(theta_new-theta_old)/(1+abs(theta_new))) < 1e-5
OR
max(abs(Sigma_new-Sigma_old)/(1+abs(Sigma_new))) < 1e-10
```

اگر همان موتور اصلی هیچ گامی نپذیرد، به همان شیوه از حل خارج می‌شود؛ فقط پرچم همگرایی داخلی روش دوم False است، به‌جای استفاده از شرط قدیمی nrtol. هیچ گام صفر، پذیرش جدید، تغییر حل خطی یا افزایش سقف تکرار افزوده نشده است. اگر حلقهٔ بیرونی متوقف شود ولی آخرین حل داخلی موفق نباشد، خروجی دوم converged=False خواهد داشت.

## شواهد حفظ موتور و ورودی

AST هر دو نسخهٔ حل داخلی، پس از کنارگذاشتن فقط شروط همگرایی و ثبت مشاهدات، با اصل دقیقاً یکسان است. AST حلقهٔ بیرونی نیز پس از کنارگذاشتن شرط توقف یکسان است. جهت GN، ژاکوبین تحلیلی، chunking، معادلات نرمال، whitening، scaling، حل مقیاس‌بندی‌شده، پذیرش با کاهش سخت‌گیرانهٔ هدف، نصف‌کردن گام و fallback میرایی حفظ شده‌اند. خود src ریپو و کد علمی ZIP ارسالی تغییر نکرده‌اند. آزمایش دارای QR استفاده نشده و native Stata اجرا نشده است.

اولین {prefix} بردار پارامتر پذیرفته‌شدهٔ NLS در دو اجرا بیت‌به‌بیت برابرند؛ مسیر تا توقف NLS اولیهٔ standard یکسان است. standard جدید، U11 کمکی ذخیره‌شدهٔ قبلی را نیز با بیشترین اختلاف پارامتر {difference['standard_vs_cached_auxiliary_max_abs_theta_difference']:.12g} و اختلاف LL {difference['standard_vs_cached_auxiliary_likelihood_difference']:+.12g} بازتولید کرده است.

{len(checks)} کنترل عددی، هویتی و ممیزی منبع پاس شده‌اند. این کنترل‌ها جای پرچم همگرایی را نمی‌گیرند. تعداد حل‌های داخلی ناموفق پیش از حل آخر: standard برابر {records['standard']['earlier_inner_failures']} و روش معیارهای nlsur برابر {records['nlsur_criteria']['earlier_inner_failures']}. جزئیات هر مرحله و هر گام پذیرفته‌شده در CSVها ثبت شده‌اند.

## فایل‌ها و بازتولید

comparison.csv، coefficient_comparison.csv، coefficient_blocks.csv، differences.json، تاریخچه‌های هر دو روش، پارامترهای تاریخی، CSV کشش‌ها، unchanged_engine_audit.json، diff تغییر شروط و experiment_manifest.json همراه گزارش‌اند.

اسکریپت‌ها از ورودی خصوصی common_inputs.pkl و auxiliary_u11.pkl ساخته‌شده در آزمایش C11 استفاده می‌کنند. برای بازتولید، ابتدا آماده‌سازی قبلی روی دادهٔ اصلی و Probitهای آن لازم است؛ این فایل‌های حاوی داده و CSV خام در بستهٔ نتایج قرار ندارند. مسیرهای REPO و C11 در ابتدای اسکریپت‌ها قابل تنظیم‌اند. سپس:

```
python run_comparison.py standard
python run_comparison.py nlsur_criteria
python build_report.py
```

دو دستور اول نتایج محلی را بازنویسی و مدل را تخمین می‌زنند؛ دستور سوم فقط نتیجه‌های ذخیره‌شده را مقایسه می‌کند. گزارش و دادهٔ اصلی منتشر نشده‌اند.
'''
(ROOT/'REPORT.fa.md').write_text(report)
target=Path('/workspace/outputs/u11-convergence-comparison.zip')
target.parent.mkdir(exist_ok=True)
files=[(ROOT/'REPORT.fa.md','REPORT.fa.md')]
files += [(p,p.name) for p in ROOT.glob('*.py')]
for pattern in ('*.json','*.csv','*.patch','*.py','*_parameter_history.npz'):
    files += [(p,'results/'+p.name) for p in OUT.glob(pattern)]
with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED) as z:
    for p,name in files: z.write(p,name)
    z.writestr('FILE_HASHES.json',json.dumps({name:hashlib.sha256(p.read_bytes()).hexdigest() for p,name in files},indent=2))
with zipfile.ZipFile(target) as z:
    assert z.testzip() is None
    assert not any(name.endswith('.pkl') or 'cldenew.csv' in name for name in z.namelist())
print(summary.to_string(index=False))
print(json.dumps(difference,indent=2))
print('Checks:',len(checks),'passed. Bundle:',target,'bytes:',target.stat().st_size)
